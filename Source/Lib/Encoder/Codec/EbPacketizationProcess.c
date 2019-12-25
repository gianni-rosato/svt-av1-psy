// clang-format off
/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
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
#define DETAILED_FRAME_OUTPUT 0

/**************************************
 * Type Declarations
 **************************************/
typedef struct EbPPSConfig
{
    uint8_t pps_id;
    uint8_t constrained_flag;
} EbPPSConfig;

/**************************************
 * Context
 **************************************/
typedef struct PacketizationContext
{
    EbDctor      dctor;
    EbFifo      *entropy_coding_input_fifo_ptr;
    EbFifo      *rate_control_tasks_output_fifo_ptr;
    EbPPSConfig *pps_config;
    EbFifo      *picture_manager_input_fifo_ptr;   // to picture-manager
    uint64_t   dpb_disp_order[8], dpb_dec_order[8];
    uint64_t   tot_shown_frames;
    uint64_t   disp_order_continuity_count;
} PacketizationContext;

static EbBool is_passthrough_data(EbLinkedListNode* data_node)
{
    return data_node->passthrough;
}

// Extracts passthrough data from a linked list. The extracted data nodes are removed from the original linked list and
// returned as a linked list. Does not gaurantee the original order of the nodes.
static EbLinkedListNode* extract_passthrough_data(EbLinkedListNode** ll_ptr_ptr)
{
    EbLinkedListNode* ll_rest_ptr = NULL;
    EbLinkedListNode* ll_pass_ptr = split_eb_linked_list(*ll_ptr_ptr, &ll_rest_ptr, is_passthrough_data);
    *ll_ptr_ptr = ll_rest_ptr;
    return ll_pass_ptr;
}

static void packetization_context_dctor(EbPtr p)
{
    EbThreadContext     *thread_context_ptr = (EbThreadContext *)p;
    PacketizationContext *obj = (PacketizationContext*)thread_context_ptr->priv;
    EB_FREE_ARRAY(obj->pps_config);
    EB_FREE_ARRAY(obj);
}

EbErrorType packetization_context_ctor(
    EbThreadContext     *thread_context_ptr,
    const EbEncHandle   *enc_handle_ptr,
    int                 rate_control_index,
    int                 demux_index)
{
    PacketizationContext  *context_ptr;
    EB_CALLOC_ARRAY(context_ptr, 1);
    thread_context_ptr->priv = context_ptr;
    thread_context_ptr->dctor = packetization_context_dctor;

    context_ptr->dctor = packetization_context_dctor;
    context_ptr->entropy_coding_input_fifo_ptr =
        eb_system_resource_get_consumer_fifo(enc_handle_ptr->entropy_coding_results_resource_ptr, 0);
    context_ptr->rate_control_tasks_output_fifo_ptr =
        eb_system_resource_get_producer_fifo(enc_handle_ptr->rate_control_tasks_resource_ptr, rate_control_index);
    context_ptr->picture_manager_input_fifo_ptr =
        eb_system_resource_get_producer_fifo(enc_handle_ptr->picture_demux_results_resource_ptr, demux_index);
    EB_MALLOC_ARRAY(context_ptr->pps_config, 1);

    return EB_ErrorNone;
}
#define TD_SIZE                     2
#define OBU_FRAME_HEADER_SIZE       3
#define TILES_GROUP_SIZE            1

// Write TD after offsetting the stream buffer
static void write_td (
    EbBufferHeaderType  *out_str_ptr,
    EbBool               show_ex,
    EbBool               has_tiles){
    uint8_t  td_buff[TD_SIZE] = { 0,0 };
    uint8_t  obu_frame_header_size = has_tiles ? OBU_FRAME_HEADER_SIZE + TILES_GROUP_SIZE : OBU_FRAME_HEADER_SIZE;
    if (out_str_ptr &&
        (out_str_ptr->n_alloc_len > (out_str_ptr->n_filled_len + 2))) {
        uint8_t *src_address = (show_ex == EB_FALSE) ?  out_str_ptr->p_buffer :
                out_str_ptr->p_buffer + out_str_ptr->n_filled_len - (obu_frame_header_size);

        uint8_t *dst_address = src_address + TD_SIZE;

        uint32_t move_size   = (show_ex == EB_FALSE) ? out_str_ptr->n_filled_len :
                               (obu_frame_header_size);

        memmove(dst_address,
                src_address,
                move_size);

        encode_td_av1((uint8_t*)(&td_buff));

        EB_MEMCPY(src_address,
                  &td_buff,
                  TD_SIZE);
    }
}

void update_rc_rate_tables(
    PictureControlSet            *picture_control_set_ptr,
    SequenceControlSet           *sequence_control_set_ptr) {
    EncodeContext* encode_context_ptr = (EncodeContext*)sequence_control_set_ptr->encode_context_ptr;

    uint32_t  intra_sad_interval_index;
    uint32_t  sad_interval_index;

    SuperBlock   *sb_ptr;
    SbParams *sb_params_ptr;

    uint32_t  sb_index;
    int32_t   qp_index;
    FrameHeader *frm_hdr = &picture_control_set_ptr->parent_pcs_ptr->frm_hdr;

    // LCU Loop
    if (sequence_control_set_ptr->static_config.rate_control_mode > 0) {
        uint64_t  sadBits[NUMBER_OF_SAD_INTERVALS] = { 0 };
        uint32_t  count[NUMBER_OF_SAD_INTERVALS] = { 0 };

        encode_context_ptr->rate_control_tables_array_updated = EB_TRUE;

        for (sb_index = 0; sb_index < picture_control_set_ptr->sb_total_count; ++sb_index) {
            sb_ptr = picture_control_set_ptr->sb_ptr_array[sb_index];
            sb_params_ptr = &sequence_control_set_ptr->sb_params_array[sb_index];

            if (sb_params_ptr->is_complete_sb) {
                if (picture_control_set_ptr->slice_type == I_SLICE) {
                    intra_sad_interval_index = picture_control_set_ptr->parent_pcs_ptr->intra_sad_interval_index[sb_index];

                    sadBits[intra_sad_interval_index] += sb_ptr->total_bits;
                    count[intra_sad_interval_index] ++;
                }
                else {
                    sad_interval_index = picture_control_set_ptr->parent_pcs_ptr->inter_sad_interval_index[sb_index];

                    sadBits[sad_interval_index] += sb_ptr->total_bits;
                    count[sad_interval_index] ++;
                }
            }
        }
        {
            eb_block_on_mutex(encode_context_ptr->rate_table_update_mutex);

            uint64_t ref_qindex_dequant = (uint64_t)picture_control_set_ptr->parent_pcs_ptr->deq.y_dequant_QTX[frm_hdr->quantization_params.base_q_idx][1];
            uint64_t sad_bits_ref_dequant = 0;
            uint64_t weight = 0;
            {
                if (picture_control_set_ptr->slice_type == I_SLICE) {
                    if (sequence_control_set_ptr->input_resolution < INPUT_SIZE_4K_RANGE) {
                        for (sad_interval_index = 0; sad_interval_index < NUMBER_OF_INTRA_SAD_INTERVALS; sad_interval_index++) {
                            if (count[sad_interval_index] > 5)
                                weight = 8;
                            else if (count[sad_interval_index] > 1)
                                weight = 5;
                            else if (count[sad_interval_index] == 1)
                                weight = 2;
                            if (count[sad_interval_index] > 0) {
                                sadBits[sad_interval_index] /= count[sad_interval_index];
                                sad_bits_ref_dequant = sadBits[sad_interval_index] * ref_qindex_dequant;
                                for (qp_index = sequence_control_set_ptr->static_config.min_qp_allowed; qp_index <= (int32_t)sequence_control_set_ptr->static_config.max_qp_allowed; qp_index++) {
                                    encode_context_ptr->rate_control_tables_array[qp_index].intra_sad_bits_array[picture_control_set_ptr->temporal_layer_index][sad_interval_index] =
                                        (EbBitNumber)(((weight * sad_bits_ref_dequant / picture_control_set_ptr->parent_pcs_ptr->deq.y_dequant_QTX[quantizer_to_qindex[qp_index]][1])
                                            + (10 - weight) * (uint32_t)encode_context_ptr->rate_control_tables_array[qp_index].intra_sad_bits_array[picture_control_set_ptr->temporal_layer_index][sad_interval_index] + 5) / 10);

                                    encode_context_ptr->rate_control_tables_array[qp_index].intra_sad_bits_array[picture_control_set_ptr->temporal_layer_index][sad_interval_index] =
                                        MIN((uint16_t)encode_context_ptr->rate_control_tables_array[qp_index].intra_sad_bits_array[picture_control_set_ptr->temporal_layer_index][sad_interval_index], (uint16_t)((1 << 15) - 1));
                                }
                            }
                        }
                    }
                    else {
                        for (sad_interval_index = 0; sad_interval_index < NUMBER_OF_INTRA_SAD_INTERVALS; sad_interval_index++) {
                            if (count[sad_interval_index] > 10)
                                weight = 8;
                            else if (count[sad_interval_index] > 5)
                                weight = 5;
                            else if (count[sad_interval_index] >= 1)
                                weight = 1;
                            if (count[sad_interval_index] > 0) {
                                sadBits[sad_interval_index] /= count[sad_interval_index];
                                sad_bits_ref_dequant = sadBits[sad_interval_index] * ref_qindex_dequant;
                                for (qp_index = sequence_control_set_ptr->static_config.min_qp_allowed; qp_index <= (int32_t)sequence_control_set_ptr->static_config.max_qp_allowed; qp_index++) {
                                    encode_context_ptr->rate_control_tables_array[qp_index].intra_sad_bits_array[picture_control_set_ptr->temporal_layer_index][sad_interval_index] =
                                        (EbBitNumber)(((weight * sad_bits_ref_dequant / picture_control_set_ptr->parent_pcs_ptr->deq.y_dequant_QTX[quantizer_to_qindex[qp_index]][1])
                                            + (10 - weight) * (uint32_t)encode_context_ptr->rate_control_tables_array[qp_index].intra_sad_bits_array[picture_control_set_ptr->temporal_layer_index][sad_interval_index] + 5) / 10);

                                    encode_context_ptr->rate_control_tables_array[qp_index].intra_sad_bits_array[picture_control_set_ptr->temporal_layer_index][sad_interval_index] =
                                        MIN((uint16_t)encode_context_ptr->rate_control_tables_array[qp_index].intra_sad_bits_array[picture_control_set_ptr->temporal_layer_index][sad_interval_index], (uint16_t)((1 << 15) - 1));
                                }
                            }
                        }
                    }
                }
                else {
                    if (sequence_control_set_ptr->input_resolution < INPUT_SIZE_4K_RANGE) {
                        for (sad_interval_index = 0; sad_interval_index < NUMBER_OF_SAD_INTERVALS; sad_interval_index++) {
                            if (count[sad_interval_index] > 5)
                                weight = 8;
                            else if (count[sad_interval_index] > 1)
                                weight = 5;
                            else if (count[sad_interval_index] == 1)
                                weight = 1;
                            if (count[sad_interval_index] > 0) {
                                sadBits[sad_interval_index] /= count[sad_interval_index];
                                sad_bits_ref_dequant = sadBits[sad_interval_index] * ref_qindex_dequant;
                                for (qp_index = sequence_control_set_ptr->static_config.min_qp_allowed; qp_index <= (int32_t)sequence_control_set_ptr->static_config.max_qp_allowed; qp_index++) {
                                    encode_context_ptr->rate_control_tables_array[qp_index].sad_bits_array[picture_control_set_ptr->temporal_layer_index][sad_interval_index] =
                                        (EbBitNumber)(((weight * sad_bits_ref_dequant / picture_control_set_ptr->parent_pcs_ptr->deq.y_dequant_QTX[quantizer_to_qindex[qp_index]][1])
                                            + (10 - weight) * (uint32_t)encode_context_ptr->rate_control_tables_array[qp_index].sad_bits_array[picture_control_set_ptr->temporal_layer_index][sad_interval_index] + 5) / 10);
                                    encode_context_ptr->rate_control_tables_array[qp_index].sad_bits_array[picture_control_set_ptr->temporal_layer_index][sad_interval_index] =
                                        MIN((uint16_t)encode_context_ptr->rate_control_tables_array[qp_index].sad_bits_array[picture_control_set_ptr->temporal_layer_index][sad_interval_index], (uint16_t)((1 << 15) - 1));
                                }
                            }
                        }
                    }
                    else {
                        for (sad_interval_index = 0; sad_interval_index < NUMBER_OF_SAD_INTERVALS; sad_interval_index++) {
                            if (count[sad_interval_index] > 10)
                                weight = 7;
                            else if (count[sad_interval_index] > 5)
                                weight = 5;
                            else if (sad_interval_index > ((NUMBER_OF_SAD_INTERVALS >> 1) - 1) && count[sad_interval_index] > 1)
                                weight = 1;
                            if (count[sad_interval_index] > 0) {
                                sadBits[sad_interval_index] /= count[sad_interval_index];
                                sad_bits_ref_dequant = sadBits[sad_interval_index] * ref_qindex_dequant;
                                for (qp_index = sequence_control_set_ptr->static_config.min_qp_allowed; qp_index <= (int32_t)sequence_control_set_ptr->static_config.max_qp_allowed; qp_index++) {
                                    encode_context_ptr->rate_control_tables_array[qp_index].sad_bits_array[picture_control_set_ptr->temporal_layer_index][sad_interval_index] =
                                        (EbBitNumber)(((weight * sad_bits_ref_dequant / picture_control_set_ptr->parent_pcs_ptr->deq.y_dequant_QTX[quantizer_to_qindex[qp_index]][1])
                                            + (10 - weight) * (uint32_t)encode_context_ptr->rate_control_tables_array[qp_index].sad_bits_array[picture_control_set_ptr->temporal_layer_index][sad_interval_index] + 5) / 10);
                                    encode_context_ptr->rate_control_tables_array[qp_index].sad_bits_array[picture_control_set_ptr->temporal_layer_index][sad_interval_index] =
                                        MIN((uint16_t)encode_context_ptr->rate_control_tables_array[qp_index].sad_bits_array[picture_control_set_ptr->temporal_layer_index][sad_interval_index], (uint16_t)((1 << 15) - 1));
                                }
                            }
                        }
                    }
                }
            }
            eb_release_mutex(encode_context_ptr->rate_table_update_mutex);
        }
    }
}
void* packetization_kernel(void *input_ptr)
{
    // Context
    EbThreadContext              *thread_context_ptr = (EbThreadContext*)input_ptr;
    PacketizationContext         *context_ptr = (PacketizationContext*)thread_context_ptr->priv;

    PictureControlSet            *picture_control_set_ptr;

    // Config
    SequenceControlSet           *sequence_control_set_ptr;

    FrameHeader                  *frm_hdr;

    // Encoding Context
    EncodeContext                *encode_context_ptr;

    // Input
    EbObjectWrapper              *entropyCodingResultsWrapperPtr;
    EntropyCodingResults         *entropyCodingResultsPtr;

    // Output
    EbObjectWrapper              *output_stream_wrapper_ptr;
    EbBufferHeaderType             *output_stream_ptr;
    EbObjectWrapper              *rateControlTasksWrapperPtr;
    RateControlTasks             *rateControlTasksPtr;
    EbObjectWrapper               *picture_manager_results_wrapper_ptr;
    PictureDemuxResults           *picture_manager_results_ptr;
    // Queue variables
    int32_t                         queueEntryIndex;
    PacketizationReorderEntry    *queueEntryPtr;
    EbLinkedListNode               *appDataLLHeadTempPtr;

    context_ptr->tot_shown_frames = 0;
    context_ptr->disp_order_continuity_count = 0;

    for (;;) {
        // Get EntropyCoding Results
        eb_get_full_object(
            context_ptr->entropy_coding_input_fifo_ptr,
            &entropyCodingResultsWrapperPtr);
        entropyCodingResultsPtr = (EntropyCodingResults*)entropyCodingResultsWrapperPtr->object_ptr;
        picture_control_set_ptr = (PictureControlSet*)entropyCodingResultsPtr->picture_control_set_wrapper_ptr->object_ptr;
        sequence_control_set_ptr = (SequenceControlSet*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr;
        encode_context_ptr = (EncodeContext*)sequence_control_set_ptr->encode_context_ptr;
        frm_hdr = &picture_control_set_ptr->parent_pcs_ptr->frm_hdr;
        //****************************************************
        // Input Entropy Results into Reordering Queue
        //****************************************************
        //get a new entry spot
        queueEntryIndex = picture_control_set_ptr->parent_pcs_ptr->decode_order % PACKETIZATION_REORDER_QUEUE_MAX_DEPTH;
        queueEntryPtr = encode_context_ptr->packetization_reorder_queue[queueEntryIndex];
        queueEntryPtr->start_time_seconds = picture_control_set_ptr->parent_pcs_ptr->start_time_seconds;
        queueEntryPtr->start_time_u_seconds = picture_control_set_ptr->parent_pcs_ptr->start_time_u_seconds;
        queueEntryPtr->is_alt_ref = picture_control_set_ptr->parent_pcs_ptr->is_alt_ref;
        eb_get_empty_object(
            sequence_control_set_ptr->encode_context_ptr->stream_output_fifo_ptr,
            &picture_control_set_ptr->parent_pcs_ptr->output_stream_wrapper_ptr);
        output_stream_wrapper_ptr = picture_control_set_ptr->parent_pcs_ptr->output_stream_wrapper_ptr;
        output_stream_ptr = (EbBufferHeaderType*)output_stream_wrapper_ptr->object_ptr;
        output_stream_ptr->p_buffer =(uint8_t*)malloc(output_stream_ptr->n_alloc_len);
        assert(output_stream_ptr->p_buffer != NULL && "bit-stream memory allocation failure");

        output_stream_ptr->flags = 0;
        output_stream_ptr->flags |= (encode_context_ptr->terminating_sequence_flag_received == EB_TRUE && picture_control_set_ptr->parent_pcs_ptr->decode_order == encode_context_ptr->terminating_picture_number) ? EB_BUFFERFLAG_EOS : 0;
        output_stream_ptr->n_filled_len = 0;
        output_stream_ptr->pts = picture_control_set_ptr->parent_pcs_ptr->input_ptr->pts;
        output_stream_ptr->dts = picture_control_set_ptr->parent_pcs_ptr->decode_order - (uint64_t)(1 << picture_control_set_ptr->parent_pcs_ptr->hierarchical_levels) + 1;
        output_stream_ptr->pic_type = picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag ?
            picture_control_set_ptr->parent_pcs_ptr->idr_flag ? EB_AV1_KEY_PICTURE :
            picture_control_set_ptr->slice_type : EB_AV1_NON_REF_PICTURE;
        output_stream_ptr->p_app_private = picture_control_set_ptr->parent_pcs_ptr->input_ptr->p_app_private;
        output_stream_ptr->qp            = picture_control_set_ptr->parent_pcs_ptr->picture_qp;

        if (sequence_control_set_ptr->static_config.stat_report){
            output_stream_ptr->luma_sse = picture_control_set_ptr->parent_pcs_ptr->luma_sse;
            output_stream_ptr->cr_sse   = picture_control_set_ptr->parent_pcs_ptr->cr_sse;
            output_stream_ptr->cb_sse   = picture_control_set_ptr->parent_pcs_ptr->cb_sse;
        } else {
            output_stream_ptr->luma_sse = 0;
            output_stream_ptr->cr_sse   = 0;
            output_stream_ptr->cb_sse   = 0;
        }

        // Get Empty Rate Control Input Tasks
        eb_get_empty_object(
            context_ptr->rate_control_tasks_output_fifo_ptr,
            &rateControlTasksWrapperPtr);
        rateControlTasksPtr = (RateControlTasks*)rateControlTasksWrapperPtr->object_ptr;
        rateControlTasksPtr->picture_control_set_wrapper_ptr = picture_control_set_ptr->picture_parent_control_set_wrapper_ptr;
        rateControlTasksPtr->task_type = RC_PACKETIZATION_FEEDBACK_RESULT;

        if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE &&
            picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr) {

            if (picture_control_set_ptr->parent_pcs_ptr->frame_end_cdf_update_mode) {
                eb_av1_reset_cdf_symbol_counters(picture_control_set_ptr->entropy_coder_ptr->fc);
                ((EbReferenceObject*)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->frame_context
                    = (*picture_control_set_ptr->entropy_coder_ptr->fc);
            }
            // Get Empty Results Object
            eb_get_empty_object(
                context_ptr->picture_manager_input_fifo_ptr,
                &picture_manager_results_wrapper_ptr);

            picture_manager_results_ptr = (PictureDemuxResults*)picture_manager_results_wrapper_ptr->object_ptr;
            picture_manager_results_ptr->picture_number = picture_control_set_ptr->picture_number;
            picture_manager_results_ptr->picture_type = EB_PIC_FEEDBACK;
            picture_manager_results_ptr->sequence_control_set_wrapper_ptr = picture_control_set_ptr->sequence_control_set_wrapper_ptr;
        }
        else {
            picture_manager_results_wrapper_ptr = EB_NULL;
            (void)picture_manager_results_ptr;
            (void)picture_manager_results_wrapper_ptr;
        }
        // Reset the bitstream before writing to it
        reset_bitstream(
            picture_control_set_ptr->bitstream_ptr->output_bitstream_ptr);

        // Code the SPS
        if (frm_hdr->frame_type == KEY_FRAME) {
            encode_sps_av1(
                picture_control_set_ptr->bitstream_ptr,
                sequence_control_set_ptr);
        }

        write_frame_header_av1(
            picture_control_set_ptr->bitstream_ptr,
            sequence_control_set_ptr,
            picture_control_set_ptr,
            0);

        // Copy Slice Header to the Output Bitstream
        copy_rbsp_bitstream_to_payload(
            picture_control_set_ptr->bitstream_ptr,
            output_stream_ptr->p_buffer,
            (uint32_t*) &(output_stream_ptr->n_filled_len),
            (uint32_t*) &(output_stream_ptr->n_alloc_len),
            encode_context_ptr);
        if (picture_control_set_ptr->parent_pcs_ptr->has_show_existing) {
            // Reset the bitstream before writing to it
            reset_bitstream(
                picture_control_set_ptr->bitstream_ptr->output_bitstream_ptr);
            write_frame_header_av1(
                picture_control_set_ptr->bitstream_ptr,
                sequence_control_set_ptr,
                picture_control_set_ptr,
                1);

            // Copy Slice Header to the Output Bitstream
            copy_rbsp_bitstream_to_payload(
                picture_control_set_ptr->bitstream_ptr,
                output_stream_ptr->p_buffer,
                (uint32_t*)&(output_stream_ptr->n_filled_len),
                (uint32_t*)&(output_stream_ptr->n_alloc_len),
                encode_context_ptr);

            output_stream_ptr->flags |= EB_BUFFERFLAG_SHOW_EXT;
        }

        // Send the number of bytes per frame to RC
        picture_control_set_ptr->parent_pcs_ptr->total_num_bits = output_stream_ptr->n_filled_len << 3;
        queueEntryPtr->total_num_bits = picture_control_set_ptr->parent_pcs_ptr->total_num_bits;
        // update the rate tables used in RC based on the encoded bits of each sb
        update_rc_rate_tables(
            picture_control_set_ptr,
            sequence_control_set_ptr);
        queueEntryPtr->frame_type = frm_hdr->frame_type;
        queueEntryPtr->poc = picture_control_set_ptr->picture_number;
        memcpy(&queueEntryPtr->av1_ref_signal, &picture_control_set_ptr->parent_pcs_ptr->av1_ref_signal, sizeof(Av1RpsNode));

        queueEntryPtr->slice_type = picture_control_set_ptr->slice_type;
#if DETAILED_FRAME_OUTPUT
        queueEntryPtr->ref_poc_list0 = picture_control_set_ptr->parent_pcs_ptr->ref_pic_poc_array[REF_LIST_0][0];
        queueEntryPtr->ref_poc_list1 = picture_control_set_ptr->parent_pcs_ptr->ref_pic_poc_array[REF_LIST_1][0];
        memcpy(queueEntryPtr->ref_poc_array, picture_control_set_ptr->parent_pcs_ptr->av1RefSignal.ref_poc_array, 7 * sizeof(uint64_t));
#endif
        queueEntryPtr->show_frame = frm_hdr->show_frame;
        queueEntryPtr->has_show_existing = picture_control_set_ptr->parent_pcs_ptr->has_show_existing;
        queueEntryPtr->show_existing_frame = frm_hdr->show_existing_frame;

        //Store the output buffer in the Queue
        queueEntryPtr->output_stream_wrapper_ptr = output_stream_wrapper_ptr;

        // Note: last chance here to add more output meta data for an encoded picture -->

        // collect output meta data
        queueEntryPtr->out_meta_data = concat_eb_linked_list(extract_passthrough_data(&(picture_control_set_ptr->parent_pcs_ptr->data_ll_head_ptr)),
            picture_control_set_ptr->parent_pcs_ptr->app_out_data_ll_head_ptr);
        picture_control_set_ptr->parent_pcs_ptr->app_out_data_ll_head_ptr = (EbLinkedListNode *)EB_NULL;

        // Calling callback functions to release the memory allocated for data linked list in the application
        while (picture_control_set_ptr->parent_pcs_ptr->data_ll_head_ptr != EB_NULL) {
            appDataLLHeadTempPtr = picture_control_set_ptr->parent_pcs_ptr->data_ll_head_ptr->next;
            if (picture_control_set_ptr->parent_pcs_ptr->data_ll_head_ptr->release_cb_fnc_ptr != EB_NULL)
                picture_control_set_ptr->parent_pcs_ptr->data_ll_head_ptr->release_cb_fnc_ptr(picture_control_set_ptr->parent_pcs_ptr->data_ll_head_ptr);
            picture_control_set_ptr->parent_pcs_ptr->data_ll_head_ptr = appDataLLHeadTempPtr;
        }

        if (sequence_control_set_ptr->static_config.speed_control_flag) {
            // update speed control variables
            eb_block_on_mutex(encode_context_ptr->sc_buffer_mutex);
            encode_context_ptr->sc_frame_out++;
            eb_release_mutex(encode_context_ptr->sc_buffer_mutex);
        }

        // Post Rate Control Taks
        eb_post_full_object(rateControlTasksWrapperPtr);
        if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE &&
            picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr)
            // Post the Full Results Object
            eb_post_full_object(picture_manager_results_wrapper_ptr);
        else
            // Since feedback is not set to PM, life count of is reduced here instead of PM
            eb_release_object(picture_control_set_ptr->sequence_control_set_wrapper_ptr);
        //Release the Parent PCS then the Child PCS
        eb_release_object(entropyCodingResultsPtr->picture_control_set_wrapper_ptr);//Child

        // Release the Entropy Coding Result
        eb_release_object(entropyCodingResultsWrapperPtr);

        //****************************************************
        // Process the head of the queue
        //****************************************************
        // Look at head of queue and see if any picture is ready to go
        queueEntryPtr = encode_context_ptr->packetization_reorder_queue[encode_context_ptr->packetization_reorder_queue_head_index];

        while (queueEntryPtr->output_stream_wrapper_ptr != EB_NULL) {
            EbBool has_tiles = (EbBool)(sequence_control_set_ptr->static_config.tile_columns || sequence_control_set_ptr->static_config.tile_rows);
            output_stream_wrapper_ptr = queueEntryPtr->output_stream_wrapper_ptr;
            output_stream_ptr = (EbBufferHeaderType*)output_stream_wrapper_ptr->object_ptr;

            if (queueEntryPtr->has_show_existing) {
                write_td(output_stream_ptr, EB_TRUE, has_tiles);
                output_stream_ptr->n_filled_len += TD_SIZE;
            }

            if (encode_context_ptr->td_needed == EB_TRUE){
                output_stream_ptr->flags |= (uint32_t)EB_BUFFERFLAG_HAS_TD;
                write_td(output_stream_ptr, EB_FALSE, has_tiles);
                encode_context_ptr->td_needed = EB_FALSE;
                output_stream_ptr->n_filled_len += TD_SIZE;
            }

            if (queueEntryPtr->has_show_existing || queueEntryPtr->show_frame)
                encode_context_ptr->td_needed = EB_TRUE;

#if DETAILED_FRAME_OUTPUT
            {
                int32_t i;
                uint8_t  showTab[] = { 'H', 'V' };

                //Countinuity count check of visible frames
                if (queueEntryPtr->show_frame) {
                    if (context_ptr->disp_order_continuity_count == queueEntryPtr->poc)
                        context_ptr->disp_order_continuity_count++;
                    else
                    {
                        SVT_LOG("SVT [ERROR]: disp_order_continuity_count Error1 POC:%i\n", (int32_t)queueEntryPtr->poc);
                        exit(0);
                    }
                }

                if (queueEntryPtr->has_show_existing) {
                    if (context_ptr->disp_order_continuity_count == context_ptr->dpb_disp_order[queueEntryPtr->show_existing_frame])
                        context_ptr->disp_order_continuity_count++;
                    else
                    {
                        SVT_LOG("SVT [ERROR]: disp_order_continuity_count Error2 POC:%i\n", (int32_t)queueEntryPtr->poc);
                        exit(0);
                    }
                }

                //update total number of shown frames
                if (queueEntryPtr->show_frame)
                    context_ptr->tot_shown_frames++;
                if (queueEntryPtr->has_show_existing)
                    context_ptr->tot_shown_frames++;

                //implement the GOP here - Serial dec order
                if (queueEntryPtr->frame_type == KEY_FRAME)
                {
                    //reset the DPB on a Key frame
                    for (i = 0; i < 8; i++)
                    {
                        context_ptr->dpb_dec_order[i] = queueEntryPtr->picture_number;
                        context_ptr->dpb_disp_order[i] = queueEntryPtr->poc;
                    }
                    SVT_LOG("%i  %i  %c ****KEY***** %i frames\n", (int32_t)queueEntryPtr->picture_number, (int32_t)queueEntryPtr->poc, showTab[queueEntryPtr->show_frame], (int32_t)context_ptr->tot_shown_frames);
                }
                else
                {
                    int32_t LASTrefIdx = queueEntryPtr->av1_ref_signal.ref_dpb_index[0];
                    int32_t BWDrefIdx = queueEntryPtr->av1_ref_signal.ref_dpb_index[4];

                    if (queueEntryPtr->frame_type == INTER_FRAME)
                    {
                        if (queueEntryPtr->has_show_existing)
                            SVT_LOG("%i (%i  %i)    %i  (%i  %i)   %c  showEx: %i   %i frames\n",
                            (int32_t)queueEntryPtr->picture_number, (int32_t)context_ptr->dpb_dec_order[LASTrefIdx], (int32_t)context_ptr->dpb_dec_order[BWDrefIdx],
                                (int32_t)queueEntryPtr->poc, (int32_t)context_ptr->dpb_disp_order[LASTrefIdx], (int32_t)context_ptr->dpb_disp_order[BWDrefIdx],
                                showTab[queueEntryPtr->show_frame], (int32_t)context_ptr->dpb_disp_order[queueEntryPtr->show_existing_frame], (int32_t)context_ptr->tot_shown_frames);
                        else
                            SVT_LOG("%i (%i  %i)    %i  (%i  %i)   %c  %i frames\n",
                            (int32_t)queueEntryPtr->picture_number, (int32_t)context_ptr->dpb_dec_order[LASTrefIdx], (int32_t)context_ptr->dpb_dec_order[BWDrefIdx],
                                (int32_t)queueEntryPtr->poc, (int32_t)context_ptr->dpb_disp_order[LASTrefIdx], (int32_t)context_ptr->dpb_disp_order[BWDrefIdx],
                                showTab[queueEntryPtr->show_frame], (int32_t)context_ptr->tot_shown_frames);

                        if (queueEntryPtr->ref_poc_list0 != context_ptr->dpb_disp_order[LASTrefIdx])
                        {
                            SVT_LOG("L0 MISMATCH POC:%i\n", (int32_t)queueEntryPtr->poc);
                            exit(0);
                        }
                        if (sequence_control_set_ptr->static_config.hierarchical_levels == 3 && queueEntryPtr->slice_type == B_SLICE && queueEntryPtr->ref_poc_list1 != context_ptr->dpb_disp_order[BWDrefIdx])
                        {
                            SVT_LOG("L1 MISMATCH POC:%i\n", (int32_t)queueEntryPtr->poc);
                            exit(0);
                        }

                        for (int rr = 0; rr < 7; rr++)
                        {
                            uint8_t dpb_spot = queueEntryPtr->av1RefSignal.refDpbIndex[rr];

                            if (queueEntryPtr->ref_poc_array[rr] != context_ptr->dpbDispOrder[dpb_spot])
                                SVT_LOG("REF_POC MISMATCH POC:%i  ref:%i\n", (int32_t)queueEntryPtr->poc, rr);
                        }
                    }
                    else
                    {
                        if (queueEntryPtr->has_show_existing)
                            SVT_LOG("%i  %i  %c   showEx: %i ----INTRA---- %i frames \n", (int32_t)queueEntryPtr->picture_number, (int32_t)queueEntryPtr->poc,
                                showTab[queueEntryPtr->show_frame], (int32_t)context_ptr->dpb_disp_order[queueEntryPtr->show_existing_frame], (int32_t)context_ptr->tot_shown_frames);
                        else
                            SVT_LOG("%i  %i  %c   ----INTRA---- %i frames\n", (int32_t)queueEntryPtr->picture_number, (int32_t)queueEntryPtr->poc,
                            (int32_t)showTab[queueEntryPtr->show_frame], (int32_t)context_ptr->tot_shown_frames);
                    }

                    //Update the DPB
                    for (i = 0; i < 8; i++)
                    {
                        if ((queueEntryPtr->av1_ref_signal.refresh_frame_mask >> i) & 1)
                        {
                            context_ptr->dpb_dec_order[i] = queueEntryPtr->picture_number;
                            context_ptr->dpb_disp_order[i] = queueEntryPtr->poc;
                        }
                    }
                }
            }
#endif
            // Calculate frame latency in milliseconds
            double latency = 0.0;
            uint64_t finishTimeSeconds = 0;
            uint64_t finishTimeuSeconds = 0;
            EbFinishTime(&finishTimeSeconds, &finishTimeuSeconds);

            EbComputeOverallElapsedTimeMs(
                queueEntryPtr->start_time_seconds,
                queueEntryPtr->start_time_u_seconds,
                finishTimeSeconds,
                finishTimeuSeconds,
                &latency);

            output_stream_ptr->n_tick_count = (uint32_t)latency;
            output_stream_ptr->p_app_private = queueEntryPtr->out_meta_data;
            if (queueEntryPtr->is_alt_ref)
                output_stream_ptr->flags |= (uint32_t)EB_BUFFERFLAG_IS_ALT_REF;

            eb_post_full_object(output_stream_wrapper_ptr);
            queueEntryPtr->out_meta_data = (EbLinkedListNode *)EB_NULL;

            // Reset the Reorder Queue Entry
            queueEntryPtr->picture_number += PACKETIZATION_REORDER_QUEUE_MAX_DEPTH;
            queueEntryPtr->output_stream_wrapper_ptr = (EbObjectWrapper *)EB_NULL;

            if (encode_context_ptr->statistics_port_active)
                queueEntryPtr->outputStatisticsWrapperPtr = (EbObjectWrapper *)EB_NULL;
            // Increment the Reorder Queue head Ptr
            encode_context_ptr->packetization_reorder_queue_head_index =
                (encode_context_ptr->packetization_reorder_queue_head_index == PACKETIZATION_REORDER_QUEUE_MAX_DEPTH - 1) ? 0 : encode_context_ptr->packetization_reorder_queue_head_index + 1;

            queueEntryPtr = encode_context_ptr->packetization_reorder_queue[encode_context_ptr->packetization_reorder_queue_head_index];
        }
    }
    return EB_NULL;
}
// clang-format on
