/*
* Copyright(c) 2019 Intel Corporation
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/
#include <stdlib.h>

#include "EbEncHandle.h"
#include "EbPacketizationProcess.h"
#include "EbEntropyCodingResults.h"
#include "EbSequenceControlSet.h"
#include "EbPictureControlSet.h"
#include "EbEntropyCoding.h"
#include "EbRateControlTasks.h"
#include "EbTime.h"
#include "EbModeDecisionProcess.h"
#include "EbPictureDemuxResults.h"
#include "EbLog.h"
#include "EbSvtAv1ErrorCodes.h"
#include "EbPictureDecisionResults.h"
#include "EbRestoration.h" // RDCOST_DBL

#define RDCOST_DBL_WITH_NATIVE_BD_DIST(RM, R, D, BD) \
    RDCOST_DBL((RM), (R), (double)((D) >> (2 * (BD - 8))))

/**************************************
 * Type Declarations
 **************************************/
typedef struct EbPPSConfig {
    uint8_t pps_id;
    uint8_t constrained_flag;
} EbPPSConfig;

/**************************************
 * Context
 **************************************/
typedef struct PacketizationContext {
    EbDctor      dctor;
    EbFifo      *entropy_coding_input_fifo_ptr;
    EbFifo      *rate_control_tasks_output_fifo_ptr;
    EbFifo      *picture_decision_results_output_fifo_ptr; // to motion estimation process
    EbPPSConfig *pps_config;
    EbFifo      *picture_demux_fifo_ptr; // to picture manager process
    uint64_t     dpb_disp_order[8], dpb_dec_order[8];
    uint64_t     tot_shown_frames;
    uint64_t     disp_order_continuity_count;
} PacketizationContext;

static EbBool is_passthrough_data(EbLinkedListNode *data_node) { return data_node->passthrough; }

int  compute_rdmult_sse(PictureControlSet *pcs_ptr, uint8_t q_index, uint8_t bit_depth);
void free_temporal_filtering_buffer(PictureControlSet *pcs_ptr, SequenceControlSet *scs_ptr);
void recon_output(PictureControlSet *pcs_ptr, SequenceControlSet *scs_ptr);
void init_resize_picture(SequenceControlSet *scs_ptr, PictureParentControlSet *pcs_ptr);
void pad_ref_and_set_flags(PictureControlSet *pcs_ptr, SequenceControlSet *scs_ptr);
void update_rc_counts(PictureParentControlSet *ppcs_ptr);
void ssim_calculations(PictureControlSet *pcs_ptr, SequenceControlSet *scs_ptr, EbBool free_memory);

// Extracts passthrough data from a linked list. The extracted data nodes are removed from the original linked list and
// returned as a linked list. Does not gaurantee the original order of the nodes.
static EbLinkedListNode *extract_passthrough_data(EbLinkedListNode **ll_ptr_ptr) {
    EbLinkedListNode *ll_rest_ptr = NULL;
    EbLinkedListNode *ll_pass_ptr = split_eb_linked_list(
        *ll_ptr_ptr, &ll_rest_ptr, is_passthrough_data);
    *ll_ptr_ptr = ll_rest_ptr;
    return ll_pass_ptr;
}

static void packetization_context_dctor(EbPtr p) {
    EbThreadContext      *thread_context_ptr = (EbThreadContext *)p;
    PacketizationContext *obj                = (PacketizationContext *)thread_context_ptr->priv;
    EB_FREE_ARRAY(obj->pps_config);
    EB_FREE_ARRAY(obj);
}

EbErrorType packetization_context_ctor(EbThreadContext   *thread_context_ptr,
                                       const EbEncHandle *enc_handle_ptr, int rate_control_index,
                                       int demux_index, int me_port_index) {
    PacketizationContext *context_ptr;
    EB_CALLOC_ARRAY(context_ptr, 1);
    thread_context_ptr->priv  = context_ptr;
    thread_context_ptr->dctor = packetization_context_dctor;

    context_ptr->dctor                         = packetization_context_dctor;
    context_ptr->entropy_coding_input_fifo_ptr = svt_system_resource_get_consumer_fifo(
        enc_handle_ptr->entropy_coding_results_resource_ptr, 0);
    context_ptr->rate_control_tasks_output_fifo_ptr = svt_system_resource_get_producer_fifo(
        enc_handle_ptr->rate_control_tasks_resource_ptr, rate_control_index);
    context_ptr->picture_demux_fifo_ptr = svt_system_resource_get_producer_fifo(
        enc_handle_ptr->picture_demux_results_resource_ptr, demux_index);
    context_ptr->picture_decision_results_output_fifo_ptr = svt_system_resource_get_producer_fifo(
        enc_handle_ptr->picture_decision_results_resource_ptr, me_port_index);
    EB_MALLOC_ARRAY(context_ptr->pps_config, 1);

    return EB_ErrorNone;
}
static inline int get_reorder_queue_pos(const EncodeContext *encode_context_ptr, int delta) {
    return (encode_context_ptr->packetization_reorder_queue_head_index + delta) %
        PACKETIZATION_REORDER_QUEUE_MAX_DEPTH;
}

static inline PacketizationReorderEntry *get_reorder_queue_entry(
    const EncodeContext *encode_context_ptr, int delta) {
    int pos = get_reorder_queue_pos(encode_context_ptr, delta);
    return encode_context_ptr->packetization_reorder_queue[pos];
}

static uint32_t count_frames_in_next_tu(const EncodeContext *encode_context_ptr,
                                        uint32_t            *data_size) {
    int i      = 0;
    *data_size = 0;
    do {
        const PacketizationReorderEntry *queue_entry_ptr = get_reorder_queue_entry(
            encode_context_ptr, i);
        const EbObjectWrapper *wrapper = queue_entry_ptr->output_stream_wrapper_ptr;
        //not a completed td
        if (!wrapper)
            return 0;

        const EbBufferHeaderType *output_stream_ptr = (EbBufferHeaderType *)wrapper->object_ptr;
        *data_size += output_stream_ptr->n_filled_len;

        i++;
        //we have a td when we got a displable frame
        if (queue_entry_ptr->show_frame)
            break;
    } while (i < PACKETIZATION_REORDER_QUEUE_MAX_DEPTH);
    return i;
}

static int pts_descend(const void *pa, const void *pb) {
    EbObjectWrapper    *a  = *(EbObjectWrapper **)pa;
    EbObjectWrapper    *b  = *(EbObjectWrapper **)pb;
    EbBufferHeaderType *ba = (EbBufferHeaderType *)(a->object_ptr);
    EbBufferHeaderType *bb = (EbBufferHeaderType *)(b->object_ptr);
    return (int)(bb->pts - ba->pts);
}

static void push_undisplayed_frame(EncodeContext *encode_context_ptr, EbObjectWrapper *wrapper) {
    if (encode_context_ptr->picture_decision_undisplayed_queue_count >= UNDISP_QUEUE_SIZE) {
        SVT_ERROR("bug, too many frames in undisplayed queue");
        return;
    }
    uint32_t count = encode_context_ptr->picture_decision_undisplayed_queue_count;
    encode_context_ptr->picture_decision_undisplayed_queue[count++] = wrapper;
    encode_context_ptr->picture_decision_undisplayed_queue_count    = count;
}

static EbObjectWrapper *pop_undisplayed_frame(EncodeContext *encode_context_ptr) {
    uint32_t count = encode_context_ptr->picture_decision_undisplayed_queue_count;
    if (!count) {
        SVT_ERROR("bug, no frame in undisplayed queue");
        return NULL;
    }
    count--;
    EbObjectWrapper *ret = encode_context_ptr->picture_decision_undisplayed_queue[count];
    encode_context_ptr->picture_decision_undisplayed_queue_count = count;
    return ret;
}

static void sort_undisplayed_frame(EncodeContext *encode_context_ptr) {
    qsort(&encode_context_ptr->picture_decision_undisplayed_queue[0],
          encode_context_ptr->picture_decision_undisplayed_queue_count,
          sizeof(EbObjectWrapper *),
          pts_descend);
}

#if DETAILED_FRAME_OUTPUT
static void print_detailed_frame_info(PacketizationContext            *context_ptr,
                                      const PacketizationReorderEntry *queue_entry_ptr) {
    int32_t i;
    uint8_t showTab[] = {'H', 'V'};

    //Countinuity count check of visible frames
    if (queue_entry_ptr->show_frame) {
        if (context_ptr->disp_order_continuity_count == queue_entry_ptr->poc)
            context_ptr->disp_order_continuity_count++;
        else {
            SVT_ERROR("disp_order_continuity_count Error1 POC:%i\n", (int32_t)queue_entry_ptr->poc);
            exit(0);
        }
    }

    if (queue_entry_ptr->has_show_existing) {
        if (context_ptr->disp_order_continuity_count ==
            context_ptr->dpb_disp_order[queue_entry_ptr->show_existing_frame])
            context_ptr->disp_order_continuity_count++;
        else {
            SVT_ERROR("disp_order_continuity_count Error2 POC:%i\n", (int32_t)queue_entry_ptr->poc);
            exit(0);
        }
    }

    //update total number of shown frames
    if (queue_entry_ptr->show_frame)
        context_ptr->tot_shown_frames++;
    if (queue_entry_ptr->has_show_existing)
        context_ptr->tot_shown_frames++;

    //implement the GOP here - Serial dec order
    if (queue_entry_ptr->frame_type == KEY_FRAME) {
        //reset the DPB on a Key frame
        for (i = 0; i < 8; i++) {
            context_ptr->dpb_dec_order[i]  = queue_entry_ptr->picture_number;
            context_ptr->dpb_disp_order[i] = queue_entry_ptr->poc;
        }
        SVT_LOG("%i  %i  %c ****KEY***** %i frames\n",
                (int32_t)queue_entry_ptr->picture_number,
                (int32_t)queue_entry_ptr->poc,
                showTab[queue_entry_ptr->show_frame],
                (int32_t)context_ptr->tot_shown_frames);
    } else {
        int32_t LASTrefIdx = queue_entry_ptr->av1_ref_signal.ref_dpb_index[0];
        int32_t BWDrefIdx  = queue_entry_ptr->av1_ref_signal.ref_dpb_index[4];

        if (queue_entry_ptr->frame_type == INTER_FRAME) {
            if (queue_entry_ptr->has_show_existing)
                SVT_LOG("%i (%i  %i)    %i  (%i  %i)   %c  showEx: %i   %i frames\n",
                        (int32_t)queue_entry_ptr->picture_number,
                        (int32_t)context_ptr->dpb_dec_order[LASTrefIdx],
                        (int32_t)context_ptr->dpb_dec_order[BWDrefIdx],
                        (int32_t)queue_entry_ptr->poc,
                        (int32_t)context_ptr->dpb_disp_order[LASTrefIdx],
                        (int32_t)context_ptr->dpb_disp_order[BWDrefIdx],
                        showTab[queue_entry_ptr->show_frame],
                        (int32_t)context_ptr->dpb_disp_order[queue_entry_ptr->show_existing_frame],
                        (int32_t)context_ptr->tot_shown_frames);
            else
                SVT_LOG("%i (%i  %i)    %i  (%i  %i)   %c  %i frames\n",
                        (int32_t)queue_entry_ptr->picture_number,
                        (int32_t)context_ptr->dpb_dec_order[LASTrefIdx],
                        (int32_t)context_ptr->dpb_dec_order[BWDrefIdx],
                        (int32_t)queue_entry_ptr->poc,
                        (int32_t)context_ptr->dpb_disp_order[LASTrefIdx],
                        (int32_t)context_ptr->dpb_disp_order[BWDrefIdx],
                        showTab[queue_entry_ptr->show_frame],
                        (int32_t)context_ptr->tot_shown_frames);

            if (queue_entry_ptr->ref_poc_list0 != context_ptr->dpb_disp_order[LASTrefIdx]) {
                SVT_LOG("L0 MISMATCH POC:%i\n", (int32_t)queue_entry_ptr->poc);
                exit(0);
            }

            for (int rr = 0; rr < 7; rr++) {
                uint8_t dpb_spot = queue_entry_ptr->av1_ref_signal.ref_dpb_index[rr];

                if (queue_entry_ptr->ref_poc_array[rr] != context_ptr->dpb_disp_order[dpb_spot])
                    SVT_LOG("REF_POC MISMATCH POC:%i  ref:%i\n", (int32_t)queue_entry_ptr->poc, rr);
            }
        } else {
            if (queue_entry_ptr->has_show_existing)
                SVT_LOG("%i  %i  %c   showEx: %i ----INTRA---- %i frames \n",
                        (int32_t)queue_entry_ptr->picture_number,
                        (int32_t)queue_entry_ptr->poc,
                        showTab[queue_entry_ptr->show_frame],
                        (int32_t)context_ptr->dpb_disp_order[queue_entry_ptr->show_existing_frame],
                        (int32_t)context_ptr->tot_shown_frames);
            else
                SVT_LOG("%i  %i  %c   ----INTRA---- %i frames\n",
                        (int32_t)queue_entry_ptr->picture_number,
                        (int32_t)queue_entry_ptr->poc,
                        (int32_t)showTab[queue_entry_ptr->show_frame],
                        (int32_t)context_ptr->tot_shown_frames);
        }

        //Update the DPB
        for (i = 0; i < 8; i++) {
            if ((queue_entry_ptr->av1_ref_signal.refresh_frame_mask >> i) & 1) {
                context_ptr->dpb_dec_order[i]  = queue_entry_ptr->picture_number;
                context_ptr->dpb_disp_order[i] = queue_entry_ptr->poc;
            }
        }
    }
}
#endif

static void collect_frames_info(PacketizationContext *context_ptr,
                                const EncodeContext *encode_context_ptr, int frames) {
    for (int i = 0; i < frames; i++) {
        PacketizationReorderEntry *queue_entry_ptr = get_reorder_queue_entry(encode_context_ptr, i);
        EbBufferHeaderType        *output_stream_ptr =
            (EbBufferHeaderType *)queue_entry_ptr->output_stream_wrapper_ptr->object_ptr;
#if DETAILED_FRAME_OUTPUT
        print_detailed_frame_info(context_ptr, queue_entry_ptr);
#else
        (void)context_ptr;

#endif
        // Calculate frame latency in milliseconds
        uint64_t finish_time_seconds   = 0;
        uint64_t finish_time_u_seconds = 0;
        svt_av1_get_time(&finish_time_seconds, &finish_time_u_seconds);

        output_stream_ptr->n_tick_count = (uint32_t)svt_av1_compute_overall_elapsed_time_ms(
            queue_entry_ptr->start_time_seconds,
            queue_entry_ptr->start_time_u_seconds,
            finish_time_seconds,
            finish_time_u_seconds);
        output_stream_ptr->p_app_private = queue_entry_ptr->out_meta_data;
        if (queue_entry_ptr->is_alt_ref)
            output_stream_ptr->flags |= (uint32_t)EB_BUFFERFLAG_IS_ALT_REF;
    }
}

#define TD_SIZE 2

//a tu start with a td, + 0 more not displable frame, + 1 display frame
static EbErrorType encode_tu(EncodeContext *encode_context_ptr, int frames, uint32_t total_bytes,
                             EbBufferHeaderType *output_stream_ptr) {
    total_bytes += TD_SIZE;
    if (total_bytes > output_stream_ptr->n_alloc_len) {
        uint8_t *pbuff;
        EB_MALLOC(pbuff, total_bytes);
        if (!pbuff) {
            SVT_ERROR("failed to allocate more memory in encode_tu");
            return EB_ErrorInsufficientResources;
        }
        EB_MEMCPY(pbuff,
                  output_stream_ptr->p_buffer,
                  output_stream_ptr->n_alloc_len > total_bytes ? total_bytes
                                                               : output_stream_ptr->n_alloc_len);
        EB_FREE(output_stream_ptr->p_buffer);
        output_stream_ptr->p_buffer    = pbuff;
        output_stream_ptr->n_alloc_len = total_bytes;
    }
    uint8_t *dst = output_stream_ptr->p_buffer + total_bytes;
    //we use last frame's output_stream_ptr to hold entire tu, so we need copy backward.
    for (int i = frames - 1; i >= 0; i--) {
        PacketizationReorderEntry *queue_entry_ptr = get_reorder_queue_entry(encode_context_ptr, i);
        EbObjectWrapper           *wrapper         = queue_entry_ptr->output_stream_wrapper_ptr;
        EbBufferHeaderType        *src_stream_ptr  = (EbBufferHeaderType *)wrapper->object_ptr;
        uint32_t                   size            = src_stream_ptr->n_filled_len;
        dst -= size;
        memmove(dst, src_stream_ptr->p_buffer, size);
        //1. The last frame is a displayable frame, others are undisplayed.
        //2. We do not push alt ref frame since the overlay frame will carry the pts.
        //3. Release alt ref stream buffer here for it will not be sent out
        if (i != frames - 1 && !queue_entry_ptr->is_alt_ref)
            push_undisplayed_frame(encode_context_ptr, wrapper);
        else if (queue_entry_ptr->is_alt_ref)
            EB_FREE(src_stream_ptr->p_buffer);
    }
    if (frames > 1)
        sort_undisplayed_frame(encode_context_ptr);
    dst -= TD_SIZE;
    encode_td_av1(dst);
    output_stream_ptr->n_filled_len = total_bytes;
    output_stream_ptr->flags |= EB_BUFFERFLAG_HAS_TD;
    return EB_ErrorNone;
}

static EbErrorType copy_data_from_bitstream(EncodeContext      *encode_context_ptr,
                                            Bitstream          *bitstream_ptr,
                                            EbBufferHeaderType *output_stream_ptr) {
    EbErrorType return_error = EB_ErrorNone;
    int         size         = bitstream_get_bytes_count(bitstream_ptr);

    CHECK_REPORT_ERROR((size + output_stream_ptr->n_filled_len < output_stream_ptr->n_alloc_len),
                       encode_context_ptr->app_callback_ptr,
                       EB_ENC_EC_ERROR2);

    void *dest = output_stream_ptr->p_buffer + output_stream_ptr->n_filled_len;
    bitstream_copy(bitstream_ptr, dest, size);
    output_stream_ptr->n_filled_len += size;

    return return_error;
}

static void encode_show_existing(EncodeContext             *encode_context_ptr,
                                 PacketizationReorderEntry *queue_entry_ptr,
                                 EbBufferHeaderType        *output_stream_ptr) {
    uint8_t *dst = output_stream_ptr->p_buffer;

    encode_td_av1(dst);
    output_stream_ptr->n_filled_len = TD_SIZE;

    copy_data_from_bitstream(encode_context_ptr, queue_entry_ptr->bitstream_ptr, output_stream_ptr);

    output_stream_ptr->flags |= (EB_BUFFERFLAG_SHOW_EXT | EB_BUFFERFLAG_HAS_TD);
}

static void release_frames(EncodeContext *encode_context_ptr, int frames) {
    for (int i = 0; i < frames; i++) {
        PacketizationReorderEntry *queue_entry_ptr = get_reorder_queue_entry(encode_context_ptr, i);
        queue_entry_ptr->out_meta_data             = (EbLinkedListNode *)NULL;
        // Reset the Reorder Queue Entry
        queue_entry_ptr->picture_number += PACKETIZATION_REORDER_QUEUE_MAX_DEPTH;
        queue_entry_ptr->output_stream_wrapper_ptr = (EbObjectWrapper *)NULL;
    }
    encode_context_ptr->packetization_reorder_queue_head_index = get_reorder_queue_pos(
        encode_context_ptr, frames);
}

inline static void clear_eos_flag(EbBufferHeaderType *output_stream_ptr) {
    output_stream_ptr->flags &= ~EB_BUFFERFLAG_EOS;
}

inline static void set_eos_flag(EbBufferHeaderType *output_stream_ptr) {
    output_stream_ptr->flags |= EB_BUFFERFLAG_EOS;
}

/* Wrapper function to capture the return of EB_MALLOC */
static inline EbErrorType malloc_p_buffer(EbBufferHeaderType *output_stream_ptr) {
    EB_MALLOC(output_stream_ptr->p_buffer, output_stream_ptr->n_alloc_len);
    return EB_ErrorNone;
}

/* Realloc when bitstream pointer size is not enough to write data of size sz */
static EbErrorType realloc_output_bitstream(Bitstream *bitstream_ptr, uint32_t sz) {
    if (bitstream_ptr && sz > 0) {
        OutputBitstreamUnit *output_bitstream_ptr = bitstream_ptr->output_bitstream_ptr;
        if (output_bitstream_ptr) {
            output_bitstream_ptr->size = sz;
            EB_REALLOC_ARRAY(output_bitstream_ptr->buffer_begin_av1, output_bitstream_ptr->size);
            output_bitstream_ptr->buffer_av1 = output_bitstream_ptr->buffer_begin_av1;
        }
    }
    return EB_ErrorNone;
}
void *packetization_kernel(void *input_ptr) {
    // Context
    EbThreadContext      *thread_context_ptr = (EbThreadContext *)input_ptr;
    PacketizationContext *context_ptr        = (PacketizationContext *)thread_context_ptr->priv;

    // Input
    EbObjectWrapper *entropy_coding_results_wrapper_ptr;

    // Output
    EbObjectWrapper *rate_control_tasks_wrapper_ptr;
    EbObjectWrapper *picture_manager_results_wrapper_ptr;

    context_ptr->tot_shown_frames            = 0;
    context_ptr->disp_order_continuity_count = 0;

    for (;;) {
        // Get EntropyCoding Results
        EB_GET_FULL_OBJECT(context_ptr->entropy_coding_input_fifo_ptr,
                           &entropy_coding_results_wrapper_ptr);

        EntropyCodingResults *entropy_coding_results_ptr =
            (EntropyCodingResults *)entropy_coding_results_wrapper_ptr->object_ptr;
        PictureControlSet *pcs_ptr = (PictureControlSet *)
                                         entropy_coding_results_ptr->pcs_wrapper_ptr->object_ptr;
        SequenceControlSet *scs_ptr = (SequenceControlSet *)pcs_ptr->scs_wrapper_ptr->object_ptr;
        EncodeContext      *encode_context_ptr = scs_ptr->encode_context_ptr;
        FrameHeader        *frm_hdr            = &pcs_ptr->parent_pcs_ptr->frm_hdr;
        Av1Common *const    cm                 = pcs_ptr->parent_pcs_ptr->av1_cm;
        uint16_t            tile_cnt = cm->tiles_info.tile_rows * cm->tiles_info.tile_cols;
        PictureParentControlSet *parent_pcs_ptr = (PictureParentControlSet *)
                                                      pcs_ptr->parent_pcs_ptr;

        if (parent_pcs_ptr->superres_total_recode_loop > 0 &&
            parent_pcs_ptr->superres_recode_loop < parent_pcs_ptr->superres_total_recode_loop) {
            // Reset the Bitstream before writing to it
            bitstream_reset(pcs_ptr->bitstream_ptr);
            write_frame_header_av1(pcs_ptr->bitstream_ptr, scs_ptr, pcs_ptr, 0);
            int64_t bits = (int64_t)bitstream_get_bytes_count(pcs_ptr->bitstream_ptr) << 3;
            int64_t rate = bits << 5; // To match scale.
            bitstream_reset(pcs_ptr->bitstream_ptr);
            int64_t sse       = parent_pcs_ptr->luma_sse;
            uint8_t bit_depth = pcs_ptr->hbd_mode_decision ? 10 : 8;
            uint8_t qindex    = parent_pcs_ptr->frm_hdr.quantization_params.base_q_idx;
            int32_t rdmult    = compute_rdmult_sse(pcs_ptr, qindex, bit_depth);

            double rdcost = RDCOST_DBL_WITH_NATIVE_BD_DIST(
                rdmult, rate, sse, scs_ptr->static_config.encoder_bit_depth);

#if DEBUG_SUPERRES_RECODE
            printf(
                "\n####### %s - frame %d, loop %d/%d, denom %d, rate %I64d, sse %I64d, rdcost "
                "%.2f, qindex %d, rdmult %d\n",
                __FUNCTION__,
                (int)parent_pcs_ptr->picture_number,
                parent_pcs_ptr->superres_recode_loop,
                parent_pcs_ptr->superres_total_recode_loop,
                parent_pcs_ptr->superres_denom,
                rate,
                sse,
                rdcost,
                qindex,
                rdmult);
#endif

            assert(parent_pcs_ptr->superres_total_recode_loop <= SCALE_NUMERATOR + 1);
            parent_pcs_ptr->superres_rdcost[parent_pcs_ptr->superres_recode_loop] = rdcost;
            ++parent_pcs_ptr->superres_recode_loop;

            if (parent_pcs_ptr->superres_recode_loop <=
                parent_pcs_ptr->superres_total_recode_loop) {
                EbBool do_recode = EB_FALSE;
                if (parent_pcs_ptr->superres_recode_loop ==
                    parent_pcs_ptr->superres_total_recode_loop) {
                    // compare rdcosts to determine whether need to recode again
                    // rdcost is the smaller the better
                    int best_index = 0;
                    for (int i = 1; i < parent_pcs_ptr->superres_total_recode_loop; ++i) {
                        double rdcost1 = parent_pcs_ptr->superres_rdcost[best_index];
                        double rdcost2 = parent_pcs_ptr->superres_rdcost[i];
                        if (rdcost2 < rdcost1) {
                            best_index = i;
                        }
                    }

                    if (best_index != parent_pcs_ptr->superres_total_recode_loop - 1) {
                        do_recode = EB_TRUE;
                        parent_pcs_ptr->superres_denom =
                            parent_pcs_ptr->superres_denom_array[best_index];
#if DEBUG_SUPERRES_RECODE
                        printf("\n####### %s - frame %d, extra loop, pick denom %d\n",
                               __FUNCTION__,
                               (int)parent_pcs_ptr->picture_number,
                               parent_pcs_ptr->superres_denom);
#endif
                    }
                } else {
                    do_recode = EB_TRUE;
                }

                if (do_recode) {
                    init_resize_picture(scs_ptr, parent_pcs_ptr);

                    // reset gm based on super-res on/off
                    set_gm_controls(parent_pcs_ptr, derive_gm_level(parent_pcs_ptr));

                    // Initialize Segments as picture decision process
                    parent_pcs_ptr->me_segments_completion_count = 0;
                    parent_pcs_ptr->me_processed_b64_count       = 0;

                    if (parent_pcs_ptr->reference_picture_wrapper_ptr != NULL) {
                        // update mi_rows and mi_cols for the reference pic wrapper (used in mfmv for other pictures)
                        EbReferenceObject *reference_object =
                            parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr;
                        svt_reference_object_reset(reference_object, scs_ptr);
                    }
#if DEBUG_SUPERRES_RECODE
                    printf("\n%s - send superres recode task to open loop ME. Frame %d, denom %d\n",
                           __FUNCTION__,
                           (int)parent_pcs_ptr->picture_number,
                           parent_pcs_ptr->superres_denom);
#endif

                    for (uint32_t segment_index = 0;
                         segment_index < parent_pcs_ptr->me_segments_total_count;
                         ++segment_index) {
                        // Get Empty Results Object
                        EbObjectWrapper *out_results_wrapper;
                        svt_get_empty_object(context_ptr->picture_decision_results_output_fifo_ptr,
                                             &out_results_wrapper);

                        PictureDecisionResults *out_results = (PictureDecisionResults *)
                                                                  out_results_wrapper->object_ptr;
                        out_results->pcs_wrapper_ptr = parent_pcs_ptr->p_pcs_wrapper_ptr;
                        out_results->segment_index   = segment_index;
                        out_results->task_type       = TASK_SUPERRES_RE_ME;
                        //Post the Full Results Object
                        svt_post_full_object(out_results_wrapper);
                    }

                    // Release the Entropy Coding Result
                    svt_release_object(entropy_coding_results_wrapper_ptr);
                    continue;
                }
            }
        }

        if (parent_pcs_ptr->superres_total_recode_loop > 0) {
            // Release pa_ref_objs
            // Delayed call from Initial Rate Control process / Source Based Operations process
            if (parent_pcs_ptr->tpl_ctrls.enable) {
                if (parent_pcs_ptr->temporal_layer_index == 0) {
                    for (uint32_t i = 0; i < parent_pcs_ptr->tpl_group_size; i++) {
                        if (parent_pcs_ptr->tpl_group[i]->slice_type == P_SLICE) {
                            if (parent_pcs_ptr->tpl_group[i]->ext_mg_id ==
                                parent_pcs_ptr->ext_mg_id + 1) {
                                release_pa_reference_objects(scs_ptr, parent_pcs_ptr->tpl_group[i]);
                            }
                        } else {
                            if (parent_pcs_ptr->tpl_group[i]->ext_mg_id ==
                                parent_pcs_ptr->ext_mg_id) {
                                release_pa_reference_objects(scs_ptr, parent_pcs_ptr->tpl_group[i]);
                            }
                        }
                    }
                }
            } else {
                release_pa_reference_objects(scs_ptr, parent_pcs_ptr);
            }

            // Delayed call from Rate Control process for multiple coding loop frames
            if (scs_ptr->static_config.pass == ENC_MIDDLE_PASS ||
                scs_ptr->static_config.pass == ENC_LAST_PASS || scs_ptr->lap_enabled ||
                (!(scs_ptr->static_config.pass == ENC_MIDDLE_PASS ||
                   scs_ptr->static_config.pass == ENC_LAST_PASS) &&
                 scs_ptr->static_config.pass != ENC_FIRST_PASS &&
                 scs_ptr->static_config.rate_control_mode == 2))
                update_rc_counts(parent_pcs_ptr);

            // Release pa me ptr. For non-superres-recode, it's released in mode_decision_kernel
            assert(pcs_ptr->parent_pcs_ptr->me_data_wrapper_ptr != NULL);
            assert(pcs_ptr->parent_pcs_ptr->pa_me_data != NULL);
            svt_release_object(pcs_ptr->parent_pcs_ptr->me_data_wrapper_ptr);
            pcs_ptr->parent_pcs_ptr->me_data_wrapper_ptr = NULL;
            pcs_ptr->parent_pcs_ptr->pa_me_data          = NULL;

            // Delayed call from Rest process
            {
                if (scs_ptr->static_config.stat_report) {
                    // memory is freed in the ssim_calculations call
                    ssim_calculations(pcs_ptr, scs_ptr, EB_TRUE);
                } else {
                    // free memory used by psnr_calculations
                    free_temporal_filtering_buffer(pcs_ptr, scs_ptr);
                }

                if (scs_ptr->static_config.recon_enabled) {
                    recon_output(pcs_ptr, scs_ptr);
                }

                if (parent_pcs_ptr->is_used_as_reference_flag) {
                    EbObjectWrapper     *picture_demux_results_wrapper_ptr;
                    PictureDemuxResults *picture_demux_results_rtr;

                    // Get Empty PicMgr Results
                    svt_get_empty_object(context_ptr->picture_demux_fifo_ptr,
                                         &picture_demux_results_wrapper_ptr);

                    picture_demux_results_rtr = (PictureDemuxResults *)
                                                    picture_demux_results_wrapper_ptr->object_ptr;
                    picture_demux_results_rtr->reference_picture_wrapper_ptr =
                        parent_pcs_ptr->reference_picture_wrapper_ptr;
                    picture_demux_results_rtr->scs_wrapper_ptr = parent_pcs_ptr->scs_wrapper_ptr;
                    picture_demux_results_rtr->picture_number  = parent_pcs_ptr->picture_number;
                    picture_demux_results_rtr->picture_type    = EB_PIC_REFERENCE;

                    // Post Reference Picture
                    svt_post_full_object(picture_demux_results_wrapper_ptr);
                }
            }

            // Delayed call from Entropy Coding process
            {
                // Release the List 0 Reference Pictures
                for (uint32_t ref_idx = 0; ref_idx < pcs_ptr->parent_pcs_ptr->ref_list0_count;
                     ++ref_idx) {
                    if (pcs_ptr->ref_pic_ptr_array[0][ref_idx] != NULL) {
                        svt_release_object(pcs_ptr->ref_pic_ptr_array[0][ref_idx]);
                    }
                }

                // Release the List 1 Reference Pictures
                for (uint32_t ref_idx = 0; ref_idx < pcs_ptr->parent_pcs_ptr->ref_list1_count;
                     ++ref_idx) {
                    if (pcs_ptr->ref_pic_ptr_array[1][ref_idx] != NULL) {
                        svt_release_object(pcs_ptr->ref_pic_ptr_array[1][ref_idx]);
                    }
                }

                //free palette data
                if (pcs_ptr->tile_tok[0][0])
                    EB_FREE_ARRAY(pcs_ptr->tile_tok[0][0]);
            }
        }

        //****************************************************
        // Input Entropy Results into Reordering Queue
        //****************************************************
        //get a new entry spot
        int32_t queue_entry_index = pcs_ptr->parent_pcs_ptr->decode_order %
            PACKETIZATION_REORDER_QUEUE_MAX_DEPTH;
        PacketizationReorderEntry *queue_entry_ptr =
            encode_context_ptr->packetization_reorder_queue[queue_entry_index];
        queue_entry_ptr->start_time_seconds   = pcs_ptr->parent_pcs_ptr->start_time_seconds;
        queue_entry_ptr->start_time_u_seconds = pcs_ptr->parent_pcs_ptr->start_time_u_seconds;
        queue_entry_ptr->is_alt_ref           = pcs_ptr->parent_pcs_ptr->is_alt_ref;
        svt_get_empty_object(scs_ptr->encode_context_ptr->stream_output_fifo_ptr,
                             &pcs_ptr->parent_pcs_ptr->output_stream_wrapper_ptr);
        EbObjectWrapper *output_stream_wrapper_ptr =
            pcs_ptr->parent_pcs_ptr->output_stream_wrapper_ptr;
        EbBufferHeaderType *output_stream_ptr = (EbBufferHeaderType *)
                                                    output_stream_wrapper_ptr->object_ptr;

        if (frm_hdr->frame_type == KEY_FRAME) {
            if (scs_ptr->static_config.mastering_display.max_luma)
                svt_add_metadata(pcs_ptr->parent_pcs_ptr->input_ptr,
                                 EB_AV1_METADATA_TYPE_HDR_MDCV,
                                 (const uint8_t *)&scs_ptr->static_config.mastering_display,
                                 sizeof(scs_ptr->static_config.mastering_display));
            if (scs_ptr->static_config.content_light_level.max_cll)
                svt_add_metadata(pcs_ptr->parent_pcs_ptr->input_ptr,
                                 EB_AV1_METADATA_TYPE_HDR_CLL,
                                 (const uint8_t *)&scs_ptr->static_config.content_light_level,
                                 sizeof(scs_ptr->static_config.content_light_level));
        }

        output_stream_ptr->flags = 0;
        if (pcs_ptr->parent_pcs_ptr->end_of_sequence_flag) {
            output_stream_ptr->flags |= EB_BUFFERFLAG_EOS;
        }
        output_stream_ptr->n_filled_len = 0;
        output_stream_ptr->pts          = pcs_ptr->parent_pcs_ptr->input_ptr->pts;
        //we output one temporal unit a time, so dts alwasy equals to pts.
        output_stream_ptr->dts           = output_stream_ptr->pts;
        output_stream_ptr->pic_type      = pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag
                 ? pcs_ptr->parent_pcs_ptr->idr_flag ? EB_AV1_KEY_PICTURE : pcs_ptr->slice_type
                 : EB_AV1_NON_REF_PICTURE;
        output_stream_ptr->p_app_private = pcs_ptr->parent_pcs_ptr->input_ptr->p_app_private;
        output_stream_ptr->qp            = pcs_ptr->parent_pcs_ptr->picture_qp;

        if (scs_ptr->static_config.stat_report) {
            output_stream_ptr->luma_sse  = pcs_ptr->parent_pcs_ptr->luma_sse;
            output_stream_ptr->cr_sse    = pcs_ptr->parent_pcs_ptr->cr_sse;
            output_stream_ptr->cb_sse    = pcs_ptr->parent_pcs_ptr->cb_sse;
            output_stream_ptr->luma_ssim = pcs_ptr->parent_pcs_ptr->luma_ssim;
            output_stream_ptr->cr_ssim   = pcs_ptr->parent_pcs_ptr->cr_ssim;
            output_stream_ptr->cb_ssim   = pcs_ptr->parent_pcs_ptr->cb_ssim;
        } else {
            output_stream_ptr->luma_sse  = 0;
            output_stream_ptr->cr_sse    = 0;
            output_stream_ptr->cb_sse    = 0;
            output_stream_ptr->luma_ssim = 0;
            output_stream_ptr->cr_ssim   = 0;
            output_stream_ptr->cb_ssim   = 0;
        }

        // Get Empty Rate Control Input Tasks
        svt_get_empty_object(context_ptr->rate_control_tasks_output_fifo_ptr,
                             &rate_control_tasks_wrapper_ptr);
        RateControlTasks *rate_control_tasks_ptr = (RateControlTasks *)
                                                       rate_control_tasks_wrapper_ptr->object_ptr;
        rate_control_tasks_ptr->pcs_wrapper_ptr = pcs_ptr->picture_parent_control_set_wrapper_ptr;
        rate_control_tasks_ptr->task_type       = RC_PACKETIZATION_FEEDBACK_RESULT;
        if (scs_ptr->enable_dec_order ||
            (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE &&
             pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr)) {
            if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE &&
                // Force each frame to update their data so future frames can use it,
                // even if the current frame did not use it.  This enables REF frames to
                // have the feature off, while NREF frames can have it on.  Used for multi-threading.
                pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr) {
                for (uint16_t tile_idx = 0; tile_idx < tile_cnt; tile_idx++) {
                    svt_av1_reset_cdf_symbol_counters(
                        pcs_ptr->entropy_coding_info[tile_idx]->entropy_coder_ptr->fc);
                    ((EbReferenceObject *)
                         pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
                        ->frame_context =
                        (*pcs_ptr->entropy_coding_info[tile_idx]->entropy_coder_ptr->fc);
                }
            }
            // Get Empty Results Object
            svt_get_empty_object(context_ptr->picture_demux_fifo_ptr,
                                 &picture_manager_results_wrapper_ptr);

            PictureDemuxResults *picture_manager_results_ptr =
                (PictureDemuxResults *)picture_manager_results_wrapper_ptr->object_ptr;
            picture_manager_results_ptr->picture_number  = pcs_ptr->picture_number;
            picture_manager_results_ptr->picture_type    = EB_PIC_FEEDBACK;
            picture_manager_results_ptr->decode_order    = pcs_ptr->parent_pcs_ptr->decode_order;
            picture_manager_results_ptr->scs_wrapper_ptr = pcs_ptr->scs_wrapper_ptr;
        }
        // Reset the Bitstream before writing to it
        bitstream_reset(pcs_ptr->bitstream_ptr);

        size_t metadata_sz = 0;

        // Code the SPS
        if (frm_hdr->frame_type == KEY_FRAME) {
            encode_sps_av1(pcs_ptr->bitstream_ptr, scs_ptr);
            // Add CLL and MDCV meta when frame is keyframe and SPS is written
            write_metadata_av1(pcs_ptr->bitstream_ptr,
                               pcs_ptr->parent_pcs_ptr->input_ptr->metadata,
                               EB_AV1_METADATA_TYPE_HDR_CLL);
            write_metadata_av1(pcs_ptr->bitstream_ptr,
                               pcs_ptr->parent_pcs_ptr->input_ptr->metadata,
                               EB_AV1_METADATA_TYPE_HDR_MDCV);
        }

        if (frm_hdr->show_frame) {
            // Add HDR10+ dynamic metadata when show frame flag is enabled
            write_metadata_av1(pcs_ptr->bitstream_ptr,
                               pcs_ptr->parent_pcs_ptr->input_ptr->metadata,
                               EB_AV1_METADATA_TYPE_ITUT_T35);
            svt_metadata_array_free(&pcs_ptr->parent_pcs_ptr->input_ptr->metadata);
        } else {
            // Copy metadata pointer to the queue entry related to current frame number
            uint64_t                   current_picture_number = pcs_ptr->picture_number;
            PacketizationReorderEntry *temp_entry =
                encode_context_ptr
                    ->packetization_reorder_queue[current_picture_number %
                                                  PACKETIZATION_REORDER_QUEUE_MAX_DEPTH];
            temp_entry->metadata = pcs_ptr->parent_pcs_ptr->input_ptr->metadata;
            pcs_ptr->parent_pcs_ptr->input_ptr->metadata = NULL;
            metadata_sz = svt_metadata_size(temp_entry->metadata, EB_AV1_METADATA_TYPE_ITUT_T35);
        }

        write_frame_header_av1(pcs_ptr->bitstream_ptr, scs_ptr, pcs_ptr, 0);

        output_stream_ptr->n_alloc_len =
            (uint32_t)(bitstream_get_bytes_count(pcs_ptr->bitstream_ptr) + TD_SIZE + metadata_sz);
        malloc_p_buffer(output_stream_ptr);

        assert(output_stream_ptr->p_buffer != NULL && "bit-stream memory allocation failure");

        copy_data_from_bitstream(encode_context_ptr, pcs_ptr->bitstream_ptr, output_stream_ptr);

        if (pcs_ptr->parent_pcs_ptr->has_show_existing) {
            uint64_t                   next_picture_number = pcs_ptr->picture_number + 1;
            PacketizationReorderEntry *temp_entry =
                encode_context_ptr
                    ->packetization_reorder_queue[next_picture_number %
                                                  PACKETIZATION_REORDER_QUEUE_MAX_DEPTH];
            // Check if the temporal entry has metadata
            if (temp_entry->metadata) {
                // Get bitstream from queue entry
                Bitstream *bitstream_ptr       = queue_entry_ptr->bitstream_ptr;
                size_t     current_metadata_sz = svt_metadata_size(temp_entry->metadata,
                                                               EB_AV1_METADATA_TYPE_ITUT_T35);
                // 16 bytes for frame header
                realloc_output_bitstream(bitstream_ptr, (uint32_t)(current_metadata_sz + 16));
            }
            // Reset the Bitstream before writing to it
            bitstream_reset(queue_entry_ptr->bitstream_ptr);
            write_metadata_av1(queue_entry_ptr->bitstream_ptr,
                               temp_entry->metadata,
                               EB_AV1_METADATA_TYPE_ITUT_T35);
            write_frame_header_av1(queue_entry_ptr->bitstream_ptr, scs_ptr, pcs_ptr, 1);
            svt_metadata_array_free(&temp_entry->metadata);
        }

        // Send the number of bytes per frame to RC
        pcs_ptr->parent_pcs_ptr->total_num_bits = output_stream_ptr->n_filled_len << 3;
        if (scs_ptr->passes == 3 && scs_ptr->static_config.pass == ENC_MIDDLE_PASS) {
            StatStruct stat_struct;
            stat_struct.poc = pcs_ptr->picture_number;
            if (scs_ptr->mid_pass_ctrls.ds)
                stat_struct.total_num_bits = pcs_ptr->parent_pcs_ptr->total_num_bits * DS_SC_FACT /
                    10;
            else
                stat_struct.total_num_bits = pcs_ptr->parent_pcs_ptr->total_num_bits;
            stat_struct.qindex       = frm_hdr->quantization_params.base_q_idx;
            stat_struct.worst_qindex = quantizer_to_qindex[(uint8_t)scs_ptr->static_config.qp];
            if (scs_ptr->rc_stat_gen_pass_mode &&
                !pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag &&
                !pcs_ptr->parent_pcs_ptr->first_frame_in_minigop)
                stat_struct.total_num_bits = 0;
            stat_struct.temporal_layer_index = pcs_ptr->temporal_layer_index;
            (scs_ptr->twopass.stats_buf_ctx->stats_in_start +
             pcs_ptr->parent_pcs_ptr->picture_number)
                ->stat_struct = stat_struct;
        }
        queue_entry_ptr->total_num_bits = pcs_ptr->parent_pcs_ptr->total_num_bits;
        queue_entry_ptr->frame_type     = frm_hdr->frame_type;
        queue_entry_ptr->poc            = pcs_ptr->picture_number;
        svt_memcpy(&queue_entry_ptr->av1_ref_signal,
                   &pcs_ptr->parent_pcs_ptr->av1_ref_signal,
                   sizeof(Av1RpsNode));

        queue_entry_ptr->slice_type = pcs_ptr->slice_type;
#if DETAILED_FRAME_OUTPUT
        queue_entry_ptr->ref_poc_list0 = pcs_ptr->parent_pcs_ptr->ref_pic_poc_array[REF_LIST_0][0];
        queue_entry_ptr->ref_poc_list1 = pcs_ptr->parent_pcs_ptr->ref_pic_poc_array[REF_LIST_1][0];
        svt_memcpy(queue_entry_ptr->ref_poc_array,
                   pcs_ptr->parent_pcs_ptr->av1_ref_signal.ref_poc_array,
                   7 * sizeof(uint64_t));
#endif
        queue_entry_ptr->show_frame          = frm_hdr->show_frame;
        queue_entry_ptr->has_show_existing   = pcs_ptr->parent_pcs_ptr->has_show_existing;
        queue_entry_ptr->show_existing_frame = frm_hdr->show_existing_frame;

        //Store the output buffer in the Queue
        queue_entry_ptr->output_stream_wrapper_ptr = output_stream_wrapper_ptr;

        // Note: last chance here to add more output meta data for an encoded picture -->

        // collect output meta data
        queue_entry_ptr->out_meta_data = concat_eb_linked_list(
            extract_passthrough_data(&(pcs_ptr->parent_pcs_ptr->data_ll_head_ptr)),
            pcs_ptr->parent_pcs_ptr->app_out_data_ll_head_ptr);
        pcs_ptr->parent_pcs_ptr->app_out_data_ll_head_ptr = (EbLinkedListNode *)NULL;

        // Calling callback functions to release the memory allocated for data linked list in the application
        while (pcs_ptr->parent_pcs_ptr->data_ll_head_ptr != NULL) {
            EbLinkedListNode *app_data_ll_head_temp_ptr =
                pcs_ptr->parent_pcs_ptr->data_ll_head_ptr->next;
            if (pcs_ptr->parent_pcs_ptr->data_ll_head_ptr->release_cb_fnc_ptr != NULL)
                pcs_ptr->parent_pcs_ptr->data_ll_head_ptr->release_cb_fnc_ptr(
                    pcs_ptr->parent_pcs_ptr->data_ll_head_ptr);
            pcs_ptr->parent_pcs_ptr->data_ll_head_ptr = app_data_ll_head_temp_ptr;
        }

        if (scs_ptr->speed_control_flag) {
            // update speed control variables
            svt_block_on_mutex(encode_context_ptr->sc_buffer_mutex);
            encode_context_ptr->sc_frame_out++;
            svt_release_mutex(encode_context_ptr->sc_buffer_mutex);
        }

        // Post Rate Control Taks
        svt_post_full_object(rate_control_tasks_wrapper_ptr);
        if (scs_ptr->enable_dec_order ||
            (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE &&
             pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr))
            // Post the Full Results Object
            svt_post_full_object(picture_manager_results_wrapper_ptr);
        else
            // Since feedback is not set to PM, life count of is reduced here instead of PM
            svt_release_object(pcs_ptr->scs_wrapper_ptr);
        svt_release_object(pcs_ptr->parent_pcs_ptr->enc_dec_ptr->enc_dec_wrapper_ptr); //Child
        //Release the Parent PCS then the Child PCS
        assert(entropy_coding_results_ptr->pcs_wrapper_ptr->live_count == 1);
        svt_release_object(entropy_coding_results_ptr->pcs_wrapper_ptr); //Child
        // Release the Entropy Coding Result
        svt_release_object(entropy_coding_results_wrapper_ptr);

        //****************************************************
        // Process the head of the queue
        //****************************************************
        // Look at head of queue and see if we got a td
        uint32_t frames, total_bytes;
        while ((frames = count_frames_in_next_tu(encode_context_ptr, &total_bytes))) {
            collect_frames_info(context_ptr, encode_context_ptr, frames);
            //last frame in a termporal unit is a displable frame. only the last frame has pts.
            //so we collect all bitstream to last frames' output buffer.
            queue_entry_ptr           = get_reorder_queue_entry(encode_context_ptr, frames - 1);
            output_stream_wrapper_ptr = queue_entry_ptr->output_stream_wrapper_ptr;
            output_stream_ptr         = (EbBufferHeaderType *)output_stream_wrapper_ptr->object_ptr;
            EbBool eos                = output_stream_ptr->flags & EB_BUFFERFLAG_EOS;

            encode_tu(encode_context_ptr, frames, total_bytes, output_stream_ptr);

            if (eos && queue_entry_ptr->has_show_existing)
                clear_eos_flag(output_stream_ptr);

            svt_post_full_object(output_stream_wrapper_ptr);
            if (queue_entry_ptr->has_show_existing) {
                EbObjectWrapper *existed = pop_undisplayed_frame(encode_context_ptr);
                if (existed) {
                    EbBufferHeaderType *existed_output_stream_ptr = (EbBufferHeaderType *)
                                                                        existed->object_ptr;
                    encode_show_existing(
                        encode_context_ptr, queue_entry_ptr, existed_output_stream_ptr);
                    if (eos)
                        set_eos_flag(existed_output_stream_ptr);
                    svt_post_full_object(existed);
                }
            }
            release_frames(encode_context_ptr, frames);
        }
    }
    return NULL;
}
