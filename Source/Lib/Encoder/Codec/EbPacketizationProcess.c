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
#define DETAILED_FRAME_OUTPUT 0

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
    EbFifo *     entropy_coding_input_fifo_ptr;
    EbFifo *     rate_control_tasks_output_fifo_ptr;
    EbPPSConfig *pps_config;
    EbFifo *     picture_manager_input_fifo_ptr; // to picture-manager
    uint64_t     dpb_disp_order[8], dpb_dec_order[8];
    uint64_t     tot_shown_frames;
    uint64_t     disp_order_continuity_count;
} PacketizationContext;

static EbBool is_passthrough_data(EbLinkedListNode *data_node) { return data_node->passthrough; }

// Extracts passthrough data from a linked list. The extracted data nodes are removed from the original linked list and
// returned as a linked list. Does not gaurantee the original order of the nodes.
static EbLinkedListNode *extract_passthrough_data(EbLinkedListNode **ll_ptr_ptr) {
    EbLinkedListNode *ll_rest_ptr = NULL;
    EbLinkedListNode *ll_pass_ptr =
        split_eb_linked_list(*ll_ptr_ptr, &ll_rest_ptr, is_passthrough_data);
    *ll_ptr_ptr = ll_rest_ptr;
    return ll_pass_ptr;
}

static void packetization_context_dctor(EbPtr p) {
    EbThreadContext *     thread_context_ptr = (EbThreadContext *)p;
    PacketizationContext *obj                = (PacketizationContext *)thread_context_ptr->priv;
    EB_FREE_ARRAY(obj->pps_config);
    EB_FREE_ARRAY(obj);
}

EbErrorType packetization_context_ctor(EbThreadContext *  thread_context_ptr,
                                       const EbEncHandle *enc_handle_ptr, int rate_control_index,
                                       int demux_index) {
    PacketizationContext *context_ptr;
    EB_CALLOC_ARRAY(context_ptr, 1);
    thread_context_ptr->priv  = context_ptr;
    thread_context_ptr->dctor = packetization_context_dctor;

    context_ptr->dctor                         = packetization_context_dctor;
    context_ptr->entropy_coding_input_fifo_ptr = svt_system_resource_get_consumer_fifo(
        enc_handle_ptr->entropy_coding_results_resource_ptr, 0);
    context_ptr->rate_control_tasks_output_fifo_ptr = svt_system_resource_get_producer_fifo(
        enc_handle_ptr->rate_control_tasks_resource_ptr, rate_control_index);
    context_ptr->picture_manager_input_fifo_ptr = svt_system_resource_get_producer_fifo(
        enc_handle_ptr->picture_demux_results_resource_ptr, demux_index);
    EB_MALLOC_ARRAY(context_ptr->pps_config, 1);

    return EB_ErrorNone;
}

void update_rc_rate_tables(PictureControlSet *pcs_ptr, SequenceControlSet *scs_ptr) {
    if (use_input_stat(scs_ptr) && scs_ptr->static_config.rate_control_mode == 1)
        return; //skip update for 2pass VBR
    // SB Loop
    if (scs_ptr->static_config.rate_control_mode > 0) {
        EncodeContext *encode_context_ptr = (EncodeContext *)scs_ptr->encode_context_ptr;

        int32_t      qp_index;
        FrameHeader *frm_hdr = &pcs_ptr->parent_pcs_ptr->frm_hdr;
        uint32_t sad_interval_index;
        uint64_t sad_bits[NUMBER_OF_SAD_INTERVALS] = {0};
        uint32_t count[NUMBER_OF_SAD_INTERVALS]    = {0};

        encode_context_ptr->rate_control_tables_array_updated = EB_TRUE;

        for (uint32_t sb_index = 0; sb_index < pcs_ptr->sb_total_count; ++sb_index) {
            SuperBlock *sb_ptr = pcs_ptr->sb_ptr_array[sb_index];
            SbParams *  sb_params_ptr = &pcs_ptr->parent_pcs_ptr->sb_params_array[sb_index];

            if (sb_params_ptr->is_complete_sb) {
                if (pcs_ptr->slice_type == I_SLICE) {
                    uint32_t intra_sad_interval_index =
                        pcs_ptr->parent_pcs_ptr->intra_sad_interval_index[sb_index];

                    sad_bits[intra_sad_interval_index] += sb_ptr->total_bits;
                    count[intra_sad_interval_index]++;
                } else {
                    sad_interval_index =
                        pcs_ptr->parent_pcs_ptr->inter_sad_interval_index[sb_index];

                    sad_bits[sad_interval_index] += sb_ptr->total_bits;
                    count[sad_interval_index]++;
                }
            }
        }
        svt_block_on_mutex(encode_context_ptr->rate_table_update_mutex);

        uint64_t ref_qindex_dequant =
            (uint64_t)pcs_ptr->parent_pcs_ptr->deq_bd
                .y_dequant_qtx[frm_hdr->quantization_params.base_q_idx][1];
        uint64_t sad_bits_ref_dequant;
        uint64_t weight               = 0;
        if (pcs_ptr->slice_type == I_SLICE) {
            if (scs_ptr->input_resolution < INPUT_SIZE_4K_RANGE) {
                for (sad_interval_index = 0;
                        sad_interval_index < NUMBER_OF_INTRA_SAD_INTERVALS;
                        sad_interval_index++) {
                    if (count[sad_interval_index] > 5)
                        weight = 8;
                    else if (count[sad_interval_index] > 1)
                        weight = 5;
                    else if (count[sad_interval_index] == 1)
                        weight = 2;
                    if (count[sad_interval_index] > 0) {
                        sad_bits[sad_interval_index] /= count[sad_interval_index];
                        sad_bits_ref_dequant =
                            sad_bits[sad_interval_index] * ref_qindex_dequant;
                        for (qp_index = scs_ptr->static_config.min_qp_allowed;
                                qp_index <= (int32_t)scs_ptr->static_config.max_qp_allowed;
                                qp_index++) {
                            encode_context_ptr->rate_control_tables_array[qp_index]
                                .intra_sad_bits_array[pcs_ptr->temporal_layer_index]
                                                        [sad_interval_index] = (EbBitNumber)(
                                ((weight * sad_bits_ref_dequant /
                                    pcs_ptr->parent_pcs_ptr->deq_bd
                                        .y_dequant_qtx[quantizer_to_qindex[qp_index]][1]) +
                                    (10 - weight) *
                                        (uint32_t)encode_context_ptr
                                            ->rate_control_tables_array[qp_index]
                                            .intra_sad_bits_array[pcs_ptr
                                                                    ->temporal_layer_index]
                                                                [sad_interval_index] +
                                    5) /
                                10);

                            encode_context_ptr->rate_control_tables_array[qp_index]
                                .intra_sad_bits_array[pcs_ptr->temporal_layer_index]
                                                        [sad_interval_index] =
                                MIN((uint16_t)encode_context_ptr
                                        ->rate_control_tables_array[qp_index]
                                        .intra_sad_bits_array[pcs_ptr->temporal_layer_index]
                                                                [sad_interval_index],
                                    (uint16_t)((1 << 15) - 1));
                        }
                    }
                }
            } else {
                for (sad_interval_index = 0;
                        sad_interval_index < NUMBER_OF_INTRA_SAD_INTERVALS;
                        sad_interval_index++) {
                    if (count[sad_interval_index] > 10)
                        weight = 8;
                    else if (count[sad_interval_index] > 5)
                        weight = 5;
                    else if (count[sad_interval_index] >= 1)
                        weight = 1;
                    if (count[sad_interval_index] > 0) {
                        sad_bits[sad_interval_index] /= count[sad_interval_index];
                        sad_bits_ref_dequant =
                            sad_bits[sad_interval_index] * ref_qindex_dequant;
                        for (qp_index = scs_ptr->static_config.min_qp_allowed;
                                qp_index <= (int32_t)scs_ptr->static_config.max_qp_allowed;
                                qp_index++) {
                            encode_context_ptr->rate_control_tables_array[qp_index]
                                .intra_sad_bits_array[pcs_ptr->temporal_layer_index]
                                                        [sad_interval_index] = (EbBitNumber)(
                                ((weight * sad_bits_ref_dequant /
                                    pcs_ptr->parent_pcs_ptr->deq_bd
                                        .y_dequant_qtx[quantizer_to_qindex[qp_index]][1]) +
                                    (10 - weight) *
                                        (uint32_t)encode_context_ptr
                                            ->rate_control_tables_array[qp_index]
                                            .intra_sad_bits_array[pcs_ptr
                                                                    ->temporal_layer_index]
                                                                [sad_interval_index] +
                                    5) /
                                10);

                            encode_context_ptr->rate_control_tables_array[qp_index]
                                .intra_sad_bits_array[pcs_ptr->temporal_layer_index]
                                                        [sad_interval_index] =
                                MIN((uint16_t)encode_context_ptr
                                        ->rate_control_tables_array[qp_index]
                                        .intra_sad_bits_array[pcs_ptr->temporal_layer_index]
                                                                [sad_interval_index],
                                    (uint16_t)((1 << 15) - 1));
                        }
                    }
                }
            }
        } else {
            if (scs_ptr->input_resolution < INPUT_SIZE_4K_RANGE) {
                for (sad_interval_index = 0; sad_interval_index < NUMBER_OF_SAD_INTERVALS;
                        sad_interval_index++) {
                    if (count[sad_interval_index] > 5)
                        weight = 8;
                    else if (count[sad_interval_index] > 1)
                        weight = 5;
                    else if (count[sad_interval_index] == 1)
                        weight = 1;
                    if (count[sad_interval_index] > 0) {
                        sad_bits[sad_interval_index] /= count[sad_interval_index];
                        sad_bits_ref_dequant =
                            sad_bits[sad_interval_index] * ref_qindex_dequant;
                        for (qp_index = scs_ptr->static_config.min_qp_allowed;
                                qp_index <= (int32_t)scs_ptr->static_config.max_qp_allowed;
                                qp_index++) {
                            encode_context_ptr->rate_control_tables_array[qp_index]
                                .sad_bits_array[pcs_ptr->temporal_layer_index]
                                                [sad_interval_index] = (EbBitNumber)(
                                ((weight * sad_bits_ref_dequant /
                                    pcs_ptr->parent_pcs_ptr->deq_bd
                                        .y_dequant_qtx[quantizer_to_qindex[qp_index]][1]) +
                                    (10 - weight) *
                                        (uint32_t)encode_context_ptr
                                            ->rate_control_tables_array[qp_index]
                                            .sad_bits_array[pcs_ptr->temporal_layer_index]
                                                        [sad_interval_index] +
                                    5) /
                                10);
                            encode_context_ptr->rate_control_tables_array[qp_index]
                                .sad_bits_array[pcs_ptr->temporal_layer_index]
                                                [sad_interval_index] =
                                MIN((uint16_t)encode_context_ptr
                                        ->rate_control_tables_array[qp_index]
                                        .sad_bits_array[pcs_ptr->temporal_layer_index]
                                                        [sad_interval_index],
                                    (uint16_t)((1 << 15) - 1));
                        }
                    }
                }
            } else {
                for (sad_interval_index = 0; sad_interval_index < NUMBER_OF_SAD_INTERVALS;
                        sad_interval_index++) {
                    if (count[sad_interval_index] > 10)
                        weight = 7;
                    else if (count[sad_interval_index] > 5)
                        weight = 5;
                    else if (sad_interval_index > ((NUMBER_OF_SAD_INTERVALS >> 1) - 1) &&
                                count[sad_interval_index] > 1)
                        weight = 1;
                    if (count[sad_interval_index] > 0) {
                        sad_bits[sad_interval_index] /= count[sad_interval_index];
                        sad_bits_ref_dequant =
                            sad_bits[sad_interval_index] * ref_qindex_dequant;
                        for (qp_index = scs_ptr->static_config.min_qp_allowed;
                                qp_index <= (int32_t)scs_ptr->static_config.max_qp_allowed;
                                qp_index++) {
                            encode_context_ptr->rate_control_tables_array[qp_index]
                                .sad_bits_array[pcs_ptr->temporal_layer_index]
                                                [sad_interval_index] = (EbBitNumber)(
                                ((weight * sad_bits_ref_dequant /
                                    pcs_ptr->parent_pcs_ptr->deq_bd
                                        .y_dequant_qtx[quantizer_to_qindex[qp_index]][1]) +
                                    (10 - weight) *
                                        (uint32_t)encode_context_ptr
                                            ->rate_control_tables_array[qp_index]
                                            .sad_bits_array[pcs_ptr->temporal_layer_index]
                                                        [sad_interval_index] +
                                    5) /
                                10);
                            encode_context_ptr->rate_control_tables_array[qp_index]
                                .sad_bits_array[pcs_ptr->temporal_layer_index]
                                                [sad_interval_index] =
                                MIN((uint16_t)encode_context_ptr
                                        ->rate_control_tables_array[qp_index]
                                        .sad_bits_array[pcs_ptr->temporal_layer_index]
                                                        [sad_interval_index],
                                    (uint16_t)((1 << 15) - 1));
                        }
                    }
                }
            }
        }
        svt_release_mutex(encode_context_ptr->rate_table_update_mutex);
    }
}

static inline int get_reorder_queue_pos(const EncodeContext *encode_context_ptr, int delta) {
    return (encode_context_ptr->packetization_reorder_queue_head_index + delta) % PACKETIZATION_REORDER_QUEUE_MAX_DEPTH;
}

static inline PacketizationReorderEntry* get_reorder_queue_entry(const EncodeContext *encode_context_ptr, int delta) {
    int pos = get_reorder_queue_pos(encode_context_ptr, delta);
    return encode_context_ptr->packetization_reorder_queue[pos];
}

static uint32_t count_frames_in_next_tu(const EncodeContext *encode_context_ptr, uint32_t *data_size) {
    int                              i = 0;
    *data_size                         = 0;
    do {
        const PacketizationReorderEntry *queue_entry_ptr = get_reorder_queue_entry(
            encode_context_ptr, i);
        const EbObjectWrapper *wrapper = queue_entry_ptr->output_stream_wrapper_ptr;
        //not a completed td
        if (!wrapper) return 0;

        const EbBufferHeaderType *output_stream_ptr = (EbBufferHeaderType *)wrapper->object_ptr;
        *data_size += output_stream_ptr->n_filled_len;

        i++;
        //we have a td when we got a displable frame
        if (queue_entry_ptr->show_frame) break;
    } while (i < PACKETIZATION_REORDER_QUEUE_MAX_DEPTH);
    return i;
}

static int pts_descend(const void *pa, const void*pb) {
    EbObjectWrapper* a = *(EbObjectWrapper**)pa;
    EbObjectWrapper* b = *(EbObjectWrapper**)pb;
    EbBufferHeaderType *ba = (EbBufferHeaderType *)(a->object_ptr);
    EbBufferHeaderType *bb = (EbBufferHeaderType *)(b->object_ptr);
    return (int)(bb->pts - ba->pts);
}

static void push_undisplayed_frame(EncodeContext *encode_context_ptr, EbObjectWrapper *wrapper) {
    if (encode_context_ptr->picture_decision_undisplayed_queue_count >= REF_FRAMES) {
        SVT_ERROR("bug, too many frames in undisplayed queue");
        return;
    }
    uint32_t count = encode_context_ptr->picture_decision_undisplayed_queue_count;
    encode_context_ptr->picture_decision_undisplayed_queue[count++] = wrapper;
    encode_context_ptr->picture_decision_undisplayed_queue_count = count;
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

static void sort_undisplayed_frame(EncodeContext* encode_context_ptr) {
    qsort(&encode_context_ptr->picture_decision_undisplayed_queue[0],
          encode_context_ptr->picture_decision_undisplayed_queue_count,
          sizeof(EbObjectWrapper *),
          pts_descend);
}

#if DETAILED_FRAME_OUTPUT
static void print_detailed_frame_info(PacketizationContext* context_ptr, const PacketizationReorderEntry *queue_entry_ptr) {
    int32_t i;
    uint8_t showTab[] = {'H', 'V'};

    //Countinuity count check of visible frames
    if (queue_entry_ptr->show_frame) {
        if (context_ptr->disp_order_continuity_count == queue_entry_ptr->poc)
            context_ptr->disp_order_continuity_count++;
        else {
            SVT_LOG("SVT [ERROR]: disp_order_continuity_count Error1 POC:%i\n",
                    (int32_t)queue_entry_ptr->poc);
            exit(0);
        }
    }

    if (queue_entry_ptr->has_show_existing) {
        if (context_ptr->disp_order_continuity_count ==
            context_ptr->dpb_disp_order[queue_entry_ptr->show_existing_frame])
            context_ptr->disp_order_continuity_count++;
        else {
            SVT_LOG("SVT [ERROR]: disp_order_continuity_count Error2 POC:%i\n",
                    (int32_t)queue_entry_ptr->poc);
            exit(0);
        }
    }

    //update total number of shown frames
    if (queue_entry_ptr->show_frame) context_ptr->tot_shown_frames++;
    if (queue_entry_ptr->has_show_existing) context_ptr->tot_shown_frames++;

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
                        (int32_t)context_ptr
                            ->dpb_disp_order[queue_entry_ptr->show_existing_frame],
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

            if (queue_entry_ptr->ref_poc_list0 !=
                context_ptr->dpb_disp_order[LASTrefIdx]) {
                SVT_LOG("L0 MISMATCH POC:%i\n", (int32_t)queue_entry_ptr->poc);
                exit(0);
            }

            for (int rr = 0; rr < 7; rr++) {
                uint8_t dpb_spot = queue_entry_ptr->av1_ref_signal.ref_dpb_index[rr];

                if (queue_entry_ptr->ref_poc_array[rr] !=
                    context_ptr->dpb_disp_order[dpb_spot])
                    SVT_LOG("REF_POC MISMATCH POC:%i  ref:%i\n",
                            (int32_t)queue_entry_ptr->poc,
                            rr);
            }
        } else {
            if (queue_entry_ptr->has_show_existing)
                SVT_LOG("%i  %i  %c   showEx: %i ----INTRA---- %i frames \n",
                        (int32_t)queue_entry_ptr->picture_number,
                        (int32_t)queue_entry_ptr->poc,
                        showTab[queue_entry_ptr->show_frame],
                        (int32_t)context_ptr
                            ->dpb_disp_order[queue_entry_ptr->show_existing_frame],
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

static void collect_frames_info(PacketizationContext* context_ptr, const EncodeContext *encode_context_ptr, int frames) {
    for (int i = 0; i < frames; i++) {
        PacketizationReorderEntry *queue_entry_ptr = get_reorder_queue_entry(encode_context_ptr, i);
        EbBufferHeaderType        *output_stream_ptr = (EbBufferHeaderType *) queue_entry_ptr->output_stream_wrapper_ptr->object_ptr;
#if DETAILED_FRAME_OUTPUT
        print_detailed_frame_info(context_ptr, queue_entry_ptr) ;
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
    uint8_t *dst                    = output_stream_ptr->p_buffer + total_bytes;
    //we use last frame's output_stream_ptr to hold entire tu, so we need copy backward.
    for (int i = frames - 1; i >= 0; i--) {
        PacketizationReorderEntry *queue_entry_ptr = get_reorder_queue_entry(encode_context_ptr, i);
        EbObjectWrapper* wrapper = queue_entry_ptr->output_stream_wrapper_ptr;
        EbBufferHeaderType *       src_stream_ptr = (EbBufferHeaderType *)wrapper->object_ptr;
        uint32_t size = src_stream_ptr->n_filled_len;
        dst -= size;
        memmove(dst, src_stream_ptr->p_buffer, size);
        //1. The last frame is a displayable frame, others are undisplayed.
        //2. We do not push alt ref frame since the overlay frame will carry the pts.
        if (i != frames - 1 && !queue_entry_ptr->is_alt_ref)
            push_undisplayed_frame(encode_context_ptr, wrapper);
    }
    if (frames > 1)
        sort_undisplayed_frame(encode_context_ptr);
    dst -= TD_SIZE;
    encode_td_av1(dst);
    output_stream_ptr->n_filled_len = total_bytes;
    output_stream_ptr->flags |= EB_BUFFERFLAG_HAS_TD;
    return EB_ErrorNone;
}

static EbErrorType copy_data_from_bitstream(EncodeContext *encode_context_ptr, Bitstream *bitstream_ptr, EbBufferHeaderType *output_stream_ptr) {
    EbErrorType          return_error = EB_ErrorNone;
    int size = bitstream_get_bytes_count(bitstream_ptr);

    CHECK_REPORT_ERROR((size + output_stream_ptr->n_filled_len
                       < output_stream_ptr->n_alloc_len),
                       encode_context_ptr->app_callback_ptr,
                       EB_ENC_EC_ERROR2);

    void* dest = output_stream_ptr->p_buffer + output_stream_ptr->n_filled_len;
    bitstream_copy(bitstream_ptr, dest, size);
    output_stream_ptr->n_filled_len += size;

    return return_error;
}

static void encode_show_existing(EncodeContext *encode_context_ptr,
                                 PacketizationReorderEntry *queue_entry_ptr,
                                 EbBufferHeaderType        *output_stream_ptr) {
    uint8_t* dst = output_stream_ptr->p_buffer;

    encode_td_av1(dst);
    output_stream_ptr->n_filled_len = TD_SIZE;

    copy_data_from_bitstream(encode_context_ptr,
                 queue_entry_ptr->bitstream_ptr,
                 output_stream_ptr);

    output_stream_ptr->flags |= (EB_BUFFERFLAG_SHOW_EXT | EB_BUFFERFLAG_HAS_TD);
}

static void release_frames(EncodeContext *encode_context_ptr, int frames) {
    for (int i = 0; i < frames; i++) {
        PacketizationReorderEntry *queue_entry_ptr = get_reorder_queue_entry(encode_context_ptr, i);
        queue_entry_ptr->out_meta_data = (EbLinkedListNode *)NULL;
        // Reset the Reorder Queue Entry
        queue_entry_ptr->picture_number += PACKETIZATION_REORDER_QUEUE_MAX_DEPTH;
        queue_entry_ptr->output_stream_wrapper_ptr = (EbObjectWrapper *)NULL;
    }
    encode_context_ptr->packetization_reorder_queue_head_index = get_reorder_queue_pos(encode_context_ptr, frames);
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

void *packetization_kernel(void *input_ptr) {
    // Context
    EbThreadContext *     thread_context_ptr = (EbThreadContext *)input_ptr;
    PacketizationContext *context_ptr        = (PacketizationContext *)thread_context_ptr->priv;

    // Input
    EbObjectWrapper *     entropy_coding_results_wrapper_ptr;

    // Output
    EbObjectWrapper *    rate_control_tasks_wrapper_ptr;
    EbObjectWrapper *    picture_manager_results_wrapper_ptr;

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
        EncodeContext *     encode_context_ptr = scs_ptr->encode_context_ptr;
        FrameHeader *    frm_hdr    = &pcs_ptr->parent_pcs_ptr->frm_hdr;
        Av1Common *const cm = pcs_ptr->parent_pcs_ptr->av1_cm;
        uint16_t            tile_cnt = cm->tiles_info.tile_rows * cm->tiles_info.tile_cols;
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

        output_stream_ptr->flags = 0;
        output_stream_ptr->flags |=
            (encode_context_ptr->terminating_sequence_flag_received == EB_TRUE &&
             pcs_ptr->parent_pcs_ptr->decode_order ==
                 encode_context_ptr->terminating_picture_number)
                ? EB_BUFFERFLAG_EOS
                : 0;
        output_stream_ptr->n_filled_len = 0;
        output_stream_ptr->pts          = pcs_ptr->parent_pcs_ptr->input_ptr->pts;
        //we output one temporal unit a time, so dts alwasy equals to pts.
        output_stream_ptr->dts          = output_stream_ptr->pts;
        output_stream_ptr->pic_type =
            pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag
                ? pcs_ptr->parent_pcs_ptr->idr_flag ? EB_AV1_KEY_PICTURE : pcs_ptr->slice_type
                : EB_AV1_NON_REF_PICTURE;
        output_stream_ptr->p_app_private = pcs_ptr->parent_pcs_ptr->input_ptr->p_app_private;
        output_stream_ptr->qp            = pcs_ptr->parent_pcs_ptr->picture_qp;

        if (scs_ptr->static_config.stat_report) {
            output_stream_ptr->luma_sse = pcs_ptr->parent_pcs_ptr->luma_sse;
            output_stream_ptr->cr_sse   = pcs_ptr->parent_pcs_ptr->cr_sse;
            output_stream_ptr->cb_sse   = pcs_ptr->parent_pcs_ptr->cb_sse;
            output_stream_ptr->luma_ssim = pcs_ptr->parent_pcs_ptr->luma_ssim;
            output_stream_ptr->cr_ssim   = pcs_ptr->parent_pcs_ptr->cr_ssim;
            output_stream_ptr->cb_ssim   = pcs_ptr->parent_pcs_ptr->cb_ssim;
        } else {
            output_stream_ptr->luma_sse = 0;
            output_stream_ptr->cr_sse   = 0;
            output_stream_ptr->cb_sse   = 0;
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

        if(use_input_stat(scs_ptr) ||
            (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE &&
            pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr)) {
            if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE &&
                pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr &&
                pcs_ptr->parent_pcs_ptr->frame_end_cdf_update_mode) {
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
            svt_get_empty_object(context_ptr->picture_manager_input_fifo_ptr,
                                 &picture_manager_results_wrapper_ptr);

            PictureDemuxResults *picture_manager_results_ptr =
                (PictureDemuxResults *)picture_manager_results_wrapper_ptr->object_ptr;
            picture_manager_results_ptr->picture_number  = pcs_ptr->picture_number;
            picture_manager_results_ptr->picture_type    = EB_PIC_FEEDBACK;
            picture_manager_results_ptr->decode_order = pcs_ptr->parent_pcs_ptr->decode_order;
            picture_manager_results_ptr->scs_wrapper_ptr = pcs_ptr->scs_wrapper_ptr;
        }
        // Reset the Bitstream before writing to it
        bitstream_reset(pcs_ptr->bitstream_ptr);

        // Code the SPS
        if (frm_hdr->frame_type == KEY_FRAME) { encode_sps_av1(pcs_ptr->bitstream_ptr, scs_ptr); }

        write_frame_header_av1(pcs_ptr->bitstream_ptr, scs_ptr, pcs_ptr, 0);

        output_stream_ptr->n_alloc_len = bitstream_get_bytes_count(pcs_ptr->bitstream_ptr) + TD_SIZE;
        malloc_p_buffer(output_stream_ptr);

        assert(output_stream_ptr->p_buffer != NULL && "bit-stream memory allocation failure");

        copy_data_from_bitstream(encode_context_ptr,
                    pcs_ptr->bitstream_ptr,
                    output_stream_ptr);

        if (pcs_ptr->parent_pcs_ptr->has_show_existing) {
            // Reset the Bitstream before writing to it
            bitstream_reset(queue_entry_ptr->bitstream_ptr);
            write_frame_header_av1(queue_entry_ptr->bitstream_ptr, scs_ptr, pcs_ptr, 1);
        }

        // Send the number of bytes per frame to RC
        pcs_ptr->parent_pcs_ptr->total_num_bits = output_stream_ptr->n_filled_len << 3;
        queue_entry_ptr->total_num_bits         = pcs_ptr->parent_pcs_ptr->total_num_bits;
        // update the rate tables used in RC based on the encoded bits of each sb
        update_rc_rate_tables(pcs_ptr, scs_ptr);
        queue_entry_ptr->frame_type = frm_hdr->frame_type;
        queue_entry_ptr->poc        = pcs_ptr->picture_number;
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

        if (scs_ptr->static_config.speed_control_flag) {
            // update speed control variables
            svt_block_on_mutex(encode_context_ptr->sc_buffer_mutex);
            encode_context_ptr->sc_frame_out++;
            svt_release_mutex(encode_context_ptr->sc_buffer_mutex);
        }

        // Post Rate Control Taks
        svt_post_full_object(rate_control_tasks_wrapper_ptr);
        if (use_input_stat(scs_ptr) || (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE &&
            pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr))
            // Post the Full Results Object
            svt_post_full_object(picture_manager_results_wrapper_ptr);
        else
            // Since feedback is not set to PM, life count of is reduced here instead of PM
            svt_release_object(pcs_ptr->scs_wrapper_ptr);
        //Release the Parent PCS then the Child PCS
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
            EbBool eos                = output_stream_ptr->flags &  EB_BUFFERFLAG_EOS;

            encode_tu(encode_context_ptr, frames, total_bytes, output_stream_ptr);

            if (eos && queue_entry_ptr->has_show_existing)
                clear_eos_flag(output_stream_ptr);

            svt_post_full_object(output_stream_wrapper_ptr);
            if (queue_entry_ptr->has_show_existing) {
                EbObjectWrapper *existed = pop_undisplayed_frame(encode_context_ptr);
                if (existed) {
                    EbBufferHeaderType *existed_output_stream_ptr = (EbBufferHeaderType *)existed->object_ptr;
                    encode_show_existing(encode_context_ptr, queue_entry_ptr, existed_output_stream_ptr);
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
