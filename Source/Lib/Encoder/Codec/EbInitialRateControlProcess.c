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

#include "EbEncHandle.h"
#include "EbPictureControlSet.h"
#include "EbSequenceControlSet.h"
#include "EbMotionEstimationResults.h"
#include "EbInitialRateControlProcess.h"
#include "EbInitialRateControlResults.h"
#include "EbMotionEstimationContext.h"
#include "EbUtility.h"
#include "EbReferenceObject.h"
#include "EbResize.h"
#include "common_dsp_rtcd.h"
#if FTR_LAD_MG
#include "EbLog.h"
#include "EbPictureDecisionProcess.h"
#endif
/**************************************
 * Context
 **************************************/
#if FTR_LAD_MG
typedef struct LadQueueEntry {
    EbDctor dctor;
    PictureParentControlSet *pcs;
} LadQueueEntry;

typedef struct LadQueue {
    LadQueueEntry **cir_buf;  //circular buffer holding the entries
    uint32_t          head;
    uint32_t          tail;
} LadQueue;

/* look ahead queue constructor*/
EbErrorType lad_queue_entry_ctor(LadQueueEntry *entry_ptr) {

    entry_ptr->pcs = NULL;
    return EB_ErrorNone;
}
#endif


typedef struct InitialRateControlContext {
    EbFifo *motion_estimation_results_input_fifo_ptr;
    EbFifo *initialrate_control_results_output_fifo_ptr;
#if FTR_LAD_MG
    LadQueue         *lad_queue;
#endif

} InitialRateControlContext;

/**************************************
* Macros
**************************************/
static void initial_rate_control_context_dctor(EbPtr p) {
    EbThreadContext *          thread_context_ptr = (EbThreadContext *)p;
    InitialRateControlContext *obj = (InitialRateControlContext *)thread_context_ptr->priv;

#if FTR_LAD_MG
    EB_DELETE_PTR_ARRAY(obj->lad_queue->cir_buf, REFERENCE_QUEUE_MAX_DEPTH);
    EB_FREE(obj->lad_queue);
#endif
    EB_FREE_ARRAY(obj);
}

/************************************************
* Initial Rate Control Context Constructor
************************************************/
EbErrorType initial_rate_control_context_ctor(EbThreadContext *  thread_context_ptr,
                                              const EbEncHandle *enc_handle_ptr) {
    InitialRateControlContext *context_ptr;
    EB_CALLOC_ARRAY(context_ptr, 1);
    thread_context_ptr->priv  = context_ptr;
    thread_context_ptr->dctor = initial_rate_control_context_dctor;

    context_ptr->motion_estimation_results_input_fifo_ptr = svt_system_resource_get_consumer_fifo(
        enc_handle_ptr->motion_estimation_results_resource_ptr, 0);
    context_ptr->initialrate_control_results_output_fifo_ptr =
        svt_system_resource_get_producer_fifo(
            enc_handle_ptr->initial_rate_control_results_resource_ptr, 0);

#if FTR_LAD_MG
    EB_MALLOC(context_ptr->lad_queue, sizeof(LadQueue));

    EB_ALLOC_PTR_ARRAY(context_ptr->lad_queue->cir_buf, REFERENCE_QUEUE_MAX_DEPTH);
    for (uint32_t picture_index = 0; picture_index < REFERENCE_QUEUE_MAX_DEPTH; ++picture_index) {
        EB_NEW(context_ptr->lad_queue->cir_buf[picture_index], lad_queue_entry_ctor);
    }
    context_ptr->lad_queue->head = 0;
    context_ptr->lad_queue->tail = 0;
#endif

    return EB_ErrorNone;
}

/************************************************
* Update BEA Information Based on Lookahead
** Average zzCost of Collocated SB throughout lookahead frames
** Set isMostOfPictureNonMoving based on number of non moving SBs
** LAD Window: min (2xmgpos+1 or sliding window size)
************************************************/

void update_bea_info_over_time(EncodeContext *          encode_context_ptr,
                               PictureParentControlSet *pcs_ptr) {
    uint64_t non_moving_index_sum = 0;

    SequenceControlSet *scs_ptr = (SequenceControlSet *)pcs_ptr->scs_wrapper_ptr->object_ptr;
    // Update motionIndexArray of the current picture by averaging the motionIndexArray of the N future pictures
    // Determine number of frames to check N
    uint32_t update_non_moving_index_array_frames_to_check = MIN(
        MIN(((pcs_ptr->pred_struct_ptr->pred_struct_period << 1) + 1), pcs_ptr->frames_in_sw),
        scs_ptr->static_config.look_ahead_distance);
    uint64_t me_dist           = 0;
    uint8_t  me_dist_pic_count = 0;
    uint32_t complete_sb_count = 0;
    // SB Loop
    for (uint16_t sb_idx = 0; sb_idx < pcs_ptr->sb_total_count; ++sb_idx) {
        uint16_t  non_moving_index_over_sliding_window = pcs_ptr->non_moving_index_array[sb_idx];
        uint16_t  frames_to_check_index;
        SbParams *sb_params = &pcs_ptr->sb_params_array[sb_idx];
        complete_sb_count++;

        // Walk the first N entries in the sliding window starting picture + 1
        uint32_t input_queue_index =
            encode_context_ptr->initial_rate_control_reorder_queue_head_index ==
                INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH - 1
            ? 0
            : encode_context_ptr->initial_rate_control_reorder_queue_head_index + 1;
        for (frames_to_check_index = 0;
             frames_to_check_index < update_non_moving_index_array_frames_to_check - 1;
             frames_to_check_index++) {
            InitialRateControlReorderEntry *temp_queue_entry_ptr =
                encode_context_ptr->initial_rate_control_reorder_queue[input_queue_index];
            PictureParentControlSet *temp_pcs_ptr =
                (PictureParentControlSet *)(temp_queue_entry_ptr->parent_pcs_wrapper_ptr)
                    ->object_ptr;

            if (temp_pcs_ptr->slice_type == I_SLICE || temp_pcs_ptr->end_of_sequence_flag)
                break;
            // Limit the distortion to lower layers 0, 1 and 2 only. Higher layers have close temporal distance and lower distortion that might contaminate the data
            if (sb_params->is_complete_sb &&
                temp_pcs_ptr->temporal_layer_index <
                    MAX((int8_t)pcs_ptr->hierarchical_levels - 1, 2)) {
                if (sb_idx == 0)
                    me_dist_pic_count++;
                me_dist += (temp_pcs_ptr->slice_type == I_SLICE)
                    ? 0
                    : (uint64_t)temp_pcs_ptr->rc_me_distortion[sb_idx];
            }
#if !TUNE_REDESIGN_TF_CTRLS
            // Store the filtered_sse of next ALT_REF picture in the I slice to be used in QP Scaling
            if (pcs_ptr->slice_type == I_SLICE && pcs_ptr->filtered_sse == 0 && sb_idx == 0 &&
                temp_pcs_ptr->temporal_layer_index == 0) {
                pcs_ptr->filtered_sse    = temp_pcs_ptr->filtered_sse;
                pcs_ptr->filtered_sse_uv = temp_pcs_ptr->filtered_sse_uv;
            }
#endif
            non_moving_index_over_sliding_window += temp_pcs_ptr->non_moving_index_array[sb_idx];

            // Increment the input_queue_index Iterator
            input_queue_index = (input_queue_index ==
                                 INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH - 1)
                ? 0
                : input_queue_index + 1;
        }
        pcs_ptr->non_moving_index_array[sb_idx] = (uint8_t)(non_moving_index_over_sliding_window /
                                                            (frames_to_check_index + 1));

        non_moving_index_sum += pcs_ptr->non_moving_index_array[sb_idx];
    }
    assert(complete_sb_count > 0);
    pcs_ptr->non_moving_index_average = (uint16_t)non_moving_index_sum / pcs_ptr->sb_total_count;
    me_dist_pic_count                 = MAX(me_dist_pic_count, 1);
    pcs_ptr->qp_scaling_average_complexity = (uint16_t)((uint64_t)me_dist / complete_sb_count /
                                                        256 / me_dist_pic_count);
    return;
}

/****************************************
* Init ZZ Cost array to default values
** Used when no Lookahead is available
****************************************/
void init_zz_cost_info(PictureParentControlSet *pcs_ptr) {
    uint16_t sb_idx;
    pcs_ptr->non_moving_index_average = INVALID_ZZ_COST;

    // SB Loop
    for (sb_idx = 0; sb_idx < pcs_ptr->sb_total_count; ++sb_idx)
        pcs_ptr->non_moving_index_array[sb_idx] = INVALID_ZZ_COST;
    return;
}

/************************************************
* Update uniform motion field
** Update Uniformly moving SBs using
** collocated SBs infor in lookahead pictures
** LAD Window: min (2xmgpos+1 or sliding window size)
************************************************/
void update_motion_field_uniformity_over_time(EncodeContext *          encode_context_ptr,
                                              SequenceControlSet *     scs_ptr,
                                              PictureParentControlSet *pcs_ptr) {
    //SVT_LOG("To update POC %d\tframesInSw = %d\n", pcs_ptr->picture_number, pcs_ptr->frames_in_sw);
    // Determine number of frames to check N
    uint32_t no_frames_to_check = MIN(
        MIN(((pcs_ptr->pred_struct_ptr->pred_struct_period << 1) + 1), pcs_ptr->frames_in_sw),
        scs_ptr->static_config.look_ahead_distance);

    // Walk the first N entries in the sliding window starting picture + 1
    uint32_t input_queue_index =
        (encode_context_ptr->initial_rate_control_reorder_queue_head_index ==
         INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH - 1)
        ? 0
        : encode_context_ptr->initial_rate_control_reorder_queue_head_index;
    for (uint32_t frames_to_check_index = 0; frames_to_check_index < no_frames_to_check - 1;
         frames_to_check_index++) {
        InitialRateControlReorderEntry *temp_queue_entry_ptr =
            encode_context_ptr->initial_rate_control_reorder_queue[input_queue_index];
        PictureParentControlSet *temp_pcs_ptr =
            ((PictureParentControlSet *)(temp_queue_entry_ptr->parent_pcs_wrapper_ptr)->object_ptr);

        if (temp_pcs_ptr->end_of_sequence_flag)
            break;
        // Increment the input_queue_index Iterator
        input_queue_index = (input_queue_index == INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH - 1)
            ? 0
            : input_queue_index + 1;
    }
    return;
}
InitialRateControlReorderEntry *determine_picture_offset_in_queue(
    EncodeContext *encode_context_ptr, PictureParentControlSet *pcs_ptr,
    MotionEstimationResults *in_results_ptr) {
    InitialRateControlReorderEntry *queue_entry_ptr;
    int32_t                         queue_entry_index;

    queue_entry_index = (int32_t)(
        pcs_ptr->picture_number -
        encode_context_ptr
            ->initial_rate_control_reorder_queue
                [encode_context_ptr->initial_rate_control_reorder_queue_head_index]
            ->picture_number);
    queue_entry_index += encode_context_ptr->initial_rate_control_reorder_queue_head_index;
    queue_entry_index = (queue_entry_index > INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH - 1)
        ? queue_entry_index - INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH
        : queue_entry_index;
    queue_entry_ptr   = encode_context_ptr->initial_rate_control_reorder_queue[queue_entry_index];
    queue_entry_ptr->parent_pcs_wrapper_ptr = in_results_ptr->pcs_wrapper_ptr;
    queue_entry_ptr->picture_number         = pcs_ptr->picture_number;

    return queue_entry_ptr;
}

void get_histogram_queue_data(SequenceControlSet *scs_ptr, EncodeContext *encode_context_ptr,
                              PictureParentControlSet *pcs_ptr) {
    HlRateControlHistogramEntry *histogram_queue_entry_ptr;
    int32_t                      histogram_queue_entry_index;

    // Determine offset from the Head Ptr for HLRC histogram queue
    svt_block_on_mutex(scs_ptr->encode_context_ptr->hl_rate_control_historgram_queue_mutex);
    histogram_queue_entry_index = (int32_t)(
        pcs_ptr->picture_number -
        encode_context_ptr
            ->hl_rate_control_historgram_queue[encode_context_ptr
                                                   ->hl_rate_control_historgram_queue_head_index]
            ->picture_number);
    histogram_queue_entry_index += encode_context_ptr->hl_rate_control_historgram_queue_head_index;
    histogram_queue_entry_index = (histogram_queue_entry_index >
                                   HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1)
        ? histogram_queue_entry_index - HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH
        : (histogram_queue_entry_index < 0)
        ? histogram_queue_entry_index + HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH
        : histogram_queue_entry_index;
    histogram_queue_entry_ptr =
        encode_context_ptr->hl_rate_control_historgram_queue[histogram_queue_entry_index];

    //histogram_queue_entry_ptr->parent_pcs_wrapper_ptr  = in_results_ptr->pcs_wrapper_ptr;
    histogram_queue_entry_ptr->picture_number       = pcs_ptr->picture_number;
    histogram_queue_entry_ptr->end_of_sequence_flag = pcs_ptr->end_of_sequence_flag;
    histogram_queue_entry_ptr->slice_type           = pcs_ptr->slice_type;
    histogram_queue_entry_ptr->temporal_layer_index = pcs_ptr->temporal_layer_index;
    histogram_queue_entry_ptr->full_sb_count        = pcs_ptr->full_sb_count;
    histogram_queue_entry_ptr->life_count           = 0;
    histogram_queue_entry_ptr->passed_to_hlrc       = EB_FALSE;
    histogram_queue_entry_ptr->is_coded             = EB_FALSE;
    histogram_queue_entry_ptr->total_num_bits_coded = 0;
    histogram_queue_entry_ptr->frames_in_sw         = 0;
    svt_memcpy(histogram_queue_entry_ptr->me_distortion_histogram,
               pcs_ptr->me_distortion_histogram,
               sizeof(uint16_t) * NUMBER_OF_SAD_INTERVALS);

    svt_memcpy(histogram_queue_entry_ptr->ois_distortion_histogram,
               pcs_ptr->ois_distortion_histogram,
               sizeof(uint16_t) * NUMBER_OF_INTRA_SAD_INTERVALS);

    svt_release_mutex(scs_ptr->encode_context_ptr->hl_rate_control_historgram_queue_mutex);
    //SVT_LOG("Test1 POC: %d\t POC: %d\t LifeCount: %d\n", histogram_queue_entry_ptr->picture_number, pcs_ptr->picture_number,  histogram_queue_entry_ptr->life_count);

    return;
}

void update_histogram_queue_entry(SequenceControlSet *scs_ptr, EncodeContext *encode_context_ptr,
                                  PictureParentControlSet *pcs_ptr, uint32_t frames_in_sw) {
    HlRateControlHistogramEntry *histogram_queue_entry_ptr;
    int32_t                      histogram_queue_entry_index;

    svt_block_on_mutex(scs_ptr->encode_context_ptr->hl_rate_control_historgram_queue_mutex);

    histogram_queue_entry_index = (int32_t)(
        pcs_ptr->picture_number -
        encode_context_ptr
            ->hl_rate_control_historgram_queue[encode_context_ptr
                                                   ->hl_rate_control_historgram_queue_head_index]
            ->picture_number);
    histogram_queue_entry_index += encode_context_ptr->hl_rate_control_historgram_queue_head_index;
    histogram_queue_entry_index = (histogram_queue_entry_index >
                                   HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1)
        ? histogram_queue_entry_index - HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH
        : (histogram_queue_entry_index < 0)
        ? histogram_queue_entry_index + HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH
        : histogram_queue_entry_index;
    histogram_queue_entry_ptr =
        encode_context_ptr->hl_rate_control_historgram_queue[histogram_queue_entry_index];
    histogram_queue_entry_ptr->passed_to_hlrc = EB_TRUE;

    if (scs_ptr->static_config.rate_control_mode == 2)
        histogram_queue_entry_ptr->life_count +=
            (int16_t)(scs_ptr->static_config.intra_period_length + 1) -
            3; // FramelevelRC does not decrease the life count for first picture in each temporal layer
    else
        histogram_queue_entry_ptr->life_count += pcs_ptr->historgram_life_count;

    histogram_queue_entry_ptr->frames_in_sw = frames_in_sw;
    svt_release_mutex(scs_ptr->encode_context_ptr->hl_rate_control_historgram_queue_mutex);

    return;
}

void svt_av1_build_quantizer(AomBitDepth bit_depth, int32_t y_dc_delta_q, int32_t u_dc_delta_q,
                             int32_t u_ac_delta_q, int32_t v_dc_delta_q, int32_t v_ac_delta_q,
                             Quants *const quants, Dequants *const deq);



#if FTR_LAD_MG
/*
 get size  of the  lad queue
*/
uint32_t  get_lad_q_size(InitialRateControlContext *ctx)
{
    uint32_t size = 0;
    uint32_t idx = ctx->lad_queue->head;
    LadQueueEntry   *queue_entry = ctx->lad_queue->cir_buf[idx];
    while (queue_entry->pcs != NULL) {
        queue_entry = ctx->lad_queue->cir_buf[++idx];
        size++;
    }
    return size;
}
#if LAD_MG_PRINT
/*
 dump the content of the  queue for debug purpose
*/
void print_lad_queue(InitialRateControlContext *ctx, uint8_t log)
{
    if (log) {
        LadQueue* queue = ctx->lad_queue;
        uint32_t idx = queue->head;
        LadQueueEntry   *queue_entry = queue->cir_buf[idx];

        SVT_LOG("\n lad_queue size:%i  ", get_lad_q_size(ctx));

        while (queue_entry->pcs != NULL) {
            SVT_LOG("%i-%lld ", queue_entry->pcs->ext_mg_id, queue_entry->pcs->picture_number);
            idx = OUT_Q_ADVANCE(idx);
            queue_entry = queue->cir_buf[idx];
        }
        SVT_LOG("\n");
    }
}
#endif
/*
 store pictures in the lad queue
*/
void push_to_lad_queue(
    PictureParentControlSet *pcs,
    InitialRateControlContext *ctx)
{
    LadQueue *queue = ctx->lad_queue;
    uint32_t entry_idx = pcs->decode_order % REFERENCE_QUEUE_MAX_DEPTH;
    LadQueueEntry   *queue_entry  = queue->cir_buf[entry_idx];
    assert_err(queue_entry->pcs == NULL, "lad queue overflow");
    if (queue_entry->pcs == NULL)
       queue_entry->pcs = pcs;
    else
        SVT_ERROR("\n lad queue overflow \n");

}

/* send picture out from irc process */
void irc_send_picture_out(InitialRateControlContext *ctx, PictureParentControlSet *pcs)
{
    EbObjectWrapper *out_results_wrapper_ptr;
    // Get Empty Results Object
    svt_get_empty_object(ctx->initialrate_control_results_output_fifo_ptr,
           &out_results_wrapper_ptr);
    InitialRateControlResults *out_results_ptr =
          (InitialRateControlResults *)out_results_wrapper_ptr->object_ptr;
    //SVT_LOG("iRC Out:%lld\n",pcs->picture_number);
    out_results_ptr->pcs_wrapper_ptr = pcs->p_pcs_wrapper_ptr;
    svt_post_full_object(out_results_wrapper_ptr);
}
#if FTR_USE_LAD_TPL
uint8_t is_frame_already_exists(PictureParentControlSet *pcs, uint32_t end_index, uint64_t pic_num) {
    for (uint32_t i = 0; i < end_index; i++)
        if (pcs->tpl_group[i]->picture_number == pic_num)
            return 1;
    return 0;
}
#endif

#if SIM_OLD_REF
/* decide to inject a frame into the tpl group*/
uint8_t  inj_to_tpl_group( PictureParentControlSet* pcs)
{
    uint8_t inj = 0;
    if (pcs->hierarchical_levels != 4) {
        if (pcs->is_used_as_reference_flag)
            inj = 1;
        else
            inj = 0;
    }
    else {

        if (pcs->temporal_layer_index < 4)
            inj = 1;
        else if (pcs->is_used_as_reference_flag && (pcs->pic_index == 6 || pcs->pic_index == 10))//TODO: could be removed once TPL r0 adapts dyncamically  to TPL group size
            inj = 1;
        else
            inj = 0;
    }

    return inj;
}
#endif
/*
 copy the number of pcs entries from the the output queue to extended  buffer
*/
void store_extended_group(
    PictureParentControlSet    *pcs,
    InitialRateControlContext  *ctx,
    uint32_t                    start_idx,
    int64_t                     end_mg)
{
    LadQueue* queue = ctx->lad_queue;
    uint32_t pic_i = 0;
    uint32_t q_idx = start_idx;
    LadQueueEntry  *entry = queue->cir_buf[q_idx];

    while (entry->pcs != NULL) {

        if (entry->pcs->ext_mg_id <= end_mg) {

            assert_err(pic_i < MAX_TPL_EXT_GROUP_SIZE,"exceeding size of ext group");
            pcs->ext_group[pic_i++] = queue->cir_buf[q_idx]->pcs;
        }

        //Increment the queue_index Iterator
        q_idx = OUT_Q_ADVANCE(q_idx);
        //get the next entry
        entry = queue->cir_buf[q_idx];
    }

    pcs->ext_group_size = pic_i;
#if LAD_MG_PRINT
    const uint8_t log = 0;
    if (log) {
        SVT_LOG("\n EXT group Pic:%lld  size:%i  \n", pcs->picture_number, pcs->ext_group_size);
        for (uint32_t i = 0; i < pcs->ext_group_size; i++) {
            if(pcs->ext_group[i]->temporal_layer_index==0)
                SVT_LOG(" | ");
            SVT_LOG("%lld ", pcs->ext_group[i]->picture_number);

        }
        SVT_LOG("\n");
    }
#endif

#if 0 //force test
   if( pcs->tpl_group_size != pcs->ext_group_size)
       printf("asfafs");
    for (uint32_t i = 0; i < pcs->ext_group_size; i++)
       if(pcs->tpl_group[i]!=pcs->ext_group[i])
           printf("asfafs");


    pcs->tpl_group_size = pcs->ext_group_size;

    for (uint32_t i = 0; i < pcs->ext_group_size; i++)
        pcs->tpl_group[i] = pcs->ext_group[i];
#endif


    //new tpl group needs to stop at the second I
    pcs->ntpl_group_size = 0;
    uint8_t is_gop_end = 0;
    int64_t last_intra_mg_id;

    for (uint32_t i = 0; i < pcs->ext_group_size; i++) {

        PictureParentControlSet *cur_pcs = pcs->ext_group[i];
        if (cur_pcs->slice_type == I_SLICE) {
            if (is_delayed_intra(cur_pcs)) {
                if (i == 0) {
                    pcs->ntpl_group[pcs->ntpl_group_size++] = cur_pcs;
                }
                else
                    break;
            }
            else {
                if (i == 0) {
                    pcs->ntpl_group[pcs->ntpl_group_size++] = cur_pcs;
                }
                else{
                    pcs->ntpl_group[pcs->ntpl_group_size++] = cur_pcs;
                    last_intra_mg_id = cur_pcs->ext_mg_id;
                    is_gop_end = 1;
                }
            }
        }
        else {

            if (is_gop_end == 0)
                pcs->ntpl_group[pcs->ntpl_group_size++] = cur_pcs;
            else if (cur_pcs->ext_mg_id == last_intra_mg_id)
                pcs->ntpl_group[pcs->ntpl_group_size++] = cur_pcs;
            else
                break;
        }
    }
#if FTR_USE_LAD_TPL

    SequenceControlSet *scs_ptr = (SequenceControlSet *)pcs->scs_wrapper_ptr->object_ptr;
    if (scs_ptr->lad_mg) {
        pcs->tpl_group_size = pcs->ntpl_group_size;
        memset(pcs->tpl_valid_pic, 0, MAX_TPL_EXT_GROUP_SIZE * sizeof(uint8_t));
        pcs->tpl_valid_pic[0] = 1;
        pcs->used_tpl_frame_num = 0;
        for (uint32_t i = 0; i < pcs->ntpl_group_size; i++) {
            pcs->tpl_group[i] = pcs->ntpl_group[i];
            // Check wether the i-th pic already exists in the tpl group
            if (!is_frame_already_exists(pcs, i, pcs->tpl_group[i]->picture_number)) {
                // Discard non-ref pic from the tpl group
#if SIM_OLD_REF
                uint8_t inject_frame = inj_to_tpl_group(pcs->tpl_group[i]);
                if (inject_frame) {
#else
                if (pcs->tpl_group[i]->is_used_as_reference_flag) {
#endif
                    if (pcs->slice_type != I_SLICE) {
                        // Discard low important pictures from tpl group
#if CLN_TPL_CONNECT_FLAG
                        if (pcs->tpl_ctrls.tpl_opt_flag && pcs->tpl_ctrls.reduced_tpl_group) {
#else
                        if (pcs->tpl_ctrls.reduced_tpl_group) {
#endif
                            if (pcs->tpl_group[i]->temporal_layer_index <= pcs->tpl_ctrls.reduced_tpl_group) {
                                pcs->tpl_valid_pic[i] = 1;
                                pcs->used_tpl_frame_num++;
                            }
                        }
                        else {
                            pcs->tpl_valid_pic[i] = 1;
                            pcs->used_tpl_frame_num++;
                        }
                    }
                    else {
                        pcs->tpl_valid_pic[i] = 1;
                        pcs->used_tpl_frame_num++;
                    }
                }
            }
        }
#if !FTR_USE_LAD_TPL
        uint16_t original_tpl_group = pcs->tpl_group_size;
        for (pic_i = original_tpl_group; pic_i < pcs->ntpl_group_size; ++pic_i) {
            PictureParentControlSet* pcs_tpl_ptr = (PictureParentControlSet *)pcs->tpl_group[pic_i];
            pcs_tpl_ptr->num_tpl_grps++;
        }
#endif
    }
#endif

#if LAD_MG_PRINT
    if (log) {
        SVT_LOG("\n NEW TPL group Pic:%lld  size:%i  \n", pcs->picture_number, pcs->ntpl_group_size);
        for (uint32_t i = 0; i < pcs->ntpl_group_size; i++) {
            if (pcs->ext_group[i]->temporal_layer_index == 0)
                SVT_LOG(" | ");
            SVT_LOG("%lld ", pcs->ntpl_group[i]->picture_number);
        }
        SVT_LOG("\n");
    }
#endif

}

/*
 scan the queue and determine if pictures can go outside
 pictures are stored in dec order.
 only base pictures are hold. the rest including LDP ones are pass-thru
*/
void process_lad_queue(
    InitialRateControlContext *ctx, uint8_t pass_thru)
{
    LadQueue *queue = ctx->lad_queue;
    LadQueueEntry  *head_entry = queue->cir_buf[queue->head];

    while (head_entry->pcs != NULL) {

        PictureParentControlSet *head_pcs = head_entry->pcs;


        uint8_t send_out;
        if (!pass_thru) {
            if (head_pcs->temporal_layer_index == 0) {

                //delayed Intra can use the whole window relative to the next Base
                uint8_t target_mgs = is_delayed_intra(head_pcs) ? head_pcs->scs_ptr->lad_mg + 1 : head_pcs->scs_ptr->lad_mg;
                target_mgs++;//add one for the MG including the head
                {

                    uint8_t num_mgs = 0;//number of MGs accumulated so far in the queue including the MG where the head belongs
                    int64_t cur_mg = head_pcs->ext_mg_id;

                    uint32_t  tmp_idx = queue->head;
                    LadQueueEntry  *tmp_entry = queue->cir_buf[tmp_idx];
                    uint32_t tot_acc_frames_in_cur_mg = 0;
                    send_out = 0;
                    while (tmp_entry->pcs != NULL) {
                        PictureParentControlSet *tmp_pcs = tmp_entry->pcs;

                        assert_err(tmp_pcs->ext_mg_id >= head_pcs->ext_mg_id, "err in mg id");
                        //adjust the lad if we hit an EOS
                        if (tmp_pcs->end_of_sequence_flag)
                            target_mgs = MIN(target_mgs, (uint8_t)(tmp_pcs->ext_mg_id - head_pcs->ext_mg_id + 1));//+1: to include the MG where the head belongs

                        if (tmp_pcs->ext_mg_id >= cur_mg) {

                            if (tmp_pcs->ext_mg_id > cur_mg)
                               assert_err(tmp_pcs->ext_mg_id == cur_mg+1, "err continuity in mg id");

                            tot_acc_frames_in_cur_mg++;

                            if (tot_acc_frames_in_cur_mg == tmp_pcs->ext_mg_size) {
                                num_mgs++;
                                cur_mg = tmp_pcs->ext_mg_id;
                                tot_acc_frames_in_cur_mg = 0;
                            }

                            if (num_mgs == target_mgs) {
                                store_extended_group(head_pcs, ctx, queue->head, tmp_pcs->ext_mg_id);
                                send_out = 1;
                                break;
                            }
                        }

                        tmp_idx = OUT_Q_ADVANCE(tmp_idx);
                        tmp_entry = queue->cir_buf[tmp_idx];
                    }
                }

            }
            else {
                send_out = 1;
            }
        }
        else {
            send_out = 1;
        }

        if (send_out) {
            //take the picture out from iRc process
            irc_send_picture_out(ctx, head_pcs);
            //advance the head
            head_entry->pcs = NULL;
            queue->head = OUT_Q_ADVANCE(queue->head);
            head_entry = queue->cir_buf[queue->head];
        }
        else {
            break;
        }
    }

}
#endif


/* Initial Rate Control Kernel */

/*********************************************************************************
*
* @brief
*  The Initial Rate Control process determines the initial bit budget for each picture
*  depending on the data gathered in the Picture Analysis and Motion Estimation processes
*  as well as the settings determined in the Picture Decision process.
*
* @par Description:
*  The Initial Rate Control process employs a sliding window buffer to analyze
*  multiple pictures if a delay is allowed. Note that no reference picture data is
*  used in this process.
*
* @param[in] Picture
*  The Initial Rate Control Kernel takes a picture and determines the initial bit budget
*  for each picture depending on the data that was gathered in Picture Analysis and
*  Motion Estimation processes
*
* @param[out] Bit Budget
*  Bit Budget is the amount of budgetted bits for a picture
*
* @remarks
*  Temporal noise reduction is currently performed in Initial Rate Control Process.
*  In the future we might decide to move it to Motion Analysis Process.
*
********************************************************************************/
void *initial_rate_control_kernel(void *input_ptr) {
    EbThreadContext *          thread_context_ptr = (EbThreadContext *)input_ptr;
    InitialRateControlContext *context_ptr = (InitialRateControlContext *)thread_context_ptr->priv;

    EbObjectWrapper *in_results_wrapper_ptr;

    EbObjectWrapper *out_results_wrapper_ptr;

    EbObjectWrapper *reference_picture_wrapper_ptr;

    // Segments
    for (;;) {
        // Get Input Full Object
        EB_GET_FULL_OBJECT(context_ptr->motion_estimation_results_input_fifo_ptr,
                           &in_results_wrapper_ptr);

        MotionEstimationResults *in_results_ptr = (MotionEstimationResults *)
                                                      in_results_wrapper_ptr->object_ptr;
        PictureParentControlSet *pcs_ptr = (PictureParentControlSet *)
                                               in_results_ptr->pcs_wrapper_ptr->object_ptr;

        // Set the segment counter
#if FTR_TPL_TR
        if (in_results_ptr->task_type == TASK_TPL_TR_ME)
        {
            pcs_ptr->me_trailing_segments_completion_count++;
            if (pcs_ptr->me_trailing_segments_completion_count == pcs_ptr->me_segments_total_count) {
#if DEBUG_TR
                printf("[%ld]: iRC in \n", pcs_ptr->picture_number);
#endif
                svt_post_semaphore(pcs_ptr->pame_trail_done_semaphore);

            }
        }
        else {

#endif
        pcs_ptr->me_segments_completion_count++;

        // If the picture is complete, proceed
        if (pcs_ptr->me_segments_completion_count == pcs_ptr->me_segments_total_count) {
            SequenceControlSet *scs_ptr = (SequenceControlSet *)
                                              pcs_ptr->scs_wrapper_ptr->object_ptr;
            EncodeContext *encode_context_ptr = (EncodeContext *)scs_ptr->encode_context_ptr;
            if (pcs_ptr->picture_number == 0) {
                Quants *const   quants_8bit = &scs_ptr->quants_8bit;
                Dequants *const deq_8bit    = &scs_ptr->deq_8bit;
                svt_av1_build_quantizer(
                    AOM_BITS_8,
                    pcs_ptr->frm_hdr.quantization_params.delta_q_dc[AOM_PLANE_Y],
                    pcs_ptr->frm_hdr.quantization_params.delta_q_dc[AOM_PLANE_U],
                    pcs_ptr->frm_hdr.quantization_params.delta_q_ac[AOM_PLANE_U],
                    pcs_ptr->frm_hdr.quantization_params.delta_q_dc[AOM_PLANE_V],
                    pcs_ptr->frm_hdr.quantization_params.delta_q_ac[AOM_PLANE_V],
                    quants_8bit,
                    deq_8bit);

                if (scs_ptr->static_config.encoder_bit_depth == AOM_BITS_10) {
                    Quants *const   quants_bd = &scs_ptr->quants_bd;
                    Dequants *const deq_bd    = &scs_ptr->deq_bd;
                    svt_av1_build_quantizer(
                        AOM_BITS_10,
                        pcs_ptr->frm_hdr.quantization_params.delta_q_dc[AOM_PLANE_Y],
                        pcs_ptr->frm_hdr.quantization_params.delta_q_dc[AOM_PLANE_U],
                        pcs_ptr->frm_hdr.quantization_params.delta_q_ac[AOM_PLANE_U],
                        pcs_ptr->frm_hdr.quantization_params.delta_q_dc[AOM_PLANE_V],
                        pcs_ptr->frm_hdr.quantization_params.delta_q_ac[AOM_PLANE_V],
                        quants_bd,
                        deq_bd);
                }
            }
            if (scs_ptr->static_config.enable_tpl_la && scs_ptr->in_loop_me == 0) {
#if FIX_DDL
                svt_set_cond_var(&pcs_ptr->me_ready, 1);
#else
                svt_post_semaphore(pcs_ptr->pame_done_semaphore);
                atomic_set_u32(&pcs_ptr->pame_done, 1);
#endif
            }
            if (scs_ptr->static_config.look_ahead_distance == 0 ||
                scs_ptr->static_config.enable_tpl_la == 0) {
                // Release Pa Ref pictures when not needed
                // Release Pa ref after when TPL is OFF
                if (!scs_ptr->in_loop_me && scs_ptr->static_config.enable_tpl_la == 0)
                    release_pa_reference_objects(scs_ptr, pcs_ptr);
            }
            /*In case Look-Ahead is zero there is no need to place pictures in the
              re-order queue. this will cause an artificial delay since pictures come in dec-order*/
            if (scs_ptr->static_config.look_ahead_distance == 0) {
                for (uint8_t temporal_layer_index = 0;
                     temporal_layer_index < EB_MAX_TEMPORAL_LAYERS;
                     temporal_layer_index++)
                    pcs_ptr->frames_in_interval[temporal_layer_index] = 0;

                pcs_ptr->frames_in_sw           = 0;
                pcs_ptr->historgram_life_count  = 0;
                pcs_ptr->scene_change_in_gop    = EB_FALSE;
                pcs_ptr->end_of_sequence_region = EB_FALSE;

                init_zz_cost_info(pcs_ptr);


#if FTR_LAD_MG
                push_to_lad_queue(pcs_ptr, context_ptr);
#if LAD_MG_PRINT
                print_lad_queue(context_ptr,0);
#endif
                uint8_t lad_queue_pass_thru = !scs_ptr->static_config.enable_tpl_la;
                process_lad_queue(context_ptr, lad_queue_pass_thru);
#else

                // Get Empty Results Object
                svt_get_empty_object(context_ptr->initialrate_control_results_output_fifo_ptr,
                                     &out_results_wrapper_ptr);

                InitialRateControlResults *out_results_ptr =
                    (InitialRateControlResults *)out_results_wrapper_ptr->object_ptr;

                out_results_ptr->pcs_wrapper_ptr = pcs_ptr->p_pcs_wrapper_ptr;
                svt_post_full_object(out_results_wrapper_ptr);
#endif

            } else {
                //****************************************************
                // Input Motion Analysis Results into Reordering Queue
                //****************************************************

                if (!pcs_ptr->is_overlay)
                    // Determine offset from the Head Ptr
                    determine_picture_offset_in_queue(encode_context_ptr, pcs_ptr, in_results_ptr);
                if (scs_ptr->static_config.rate_control_mode && !use_input_stat(scs_ptr) &&
                    !scs_ptr->lap_enabled)
                    // Getting the Histogram Queue Data
                    get_histogram_queue_data(scs_ptr, encode_context_ptr, pcs_ptr);
                for (uint8_t temporal_layer_index = 0;
                     temporal_layer_index < EB_MAX_TEMPORAL_LAYERS;
                     temporal_layer_index++)
                    pcs_ptr->frames_in_interval[temporal_layer_index] = 0;
                pcs_ptr->frames_in_sw          = 0;
                pcs_ptr->historgram_life_count = 0;
                pcs_ptr->scene_change_in_gop   = EB_FALSE;
                EbBool move_slide_window_flag  = EB_TRUE;
                while (move_slide_window_flag) {
                    // Check if the sliding window condition is valid
                    uint32_t queue_entry_index_temp =
                        encode_context_ptr->initial_rate_control_reorder_queue_head_index;
                    EbBool end_of_sequence_flag =
                        encode_context_ptr
                            ->initial_rate_control_reorder_queue[queue_entry_index_temp]
                            ->parent_pcs_wrapper_ptr
                        ? (((PictureParentControlSet
                                 *)(encode_context_ptr
                                        ->initial_rate_control_reorder_queue[queue_entry_index_temp]
                                        ->parent_pcs_wrapper_ptr)
                                ->object_ptr))
                              ->end_of_sequence_flag
                        : EB_FALSE;
                    uint8_t frames_in_sw = 0;
                    while (move_slide_window_flag && !end_of_sequence_flag &&
                           queue_entry_index_temp <=
                               encode_context_ptr->initial_rate_control_reorder_queue_head_index +
                                   scs_ptr->static_config.look_ahead_distance) {
                        // frames_in_sw <= scs_ptr->static_config.look_ahead_distance){
                        frames_in_sw++;

                        uint32_t queue_entry_index_temp2 =
                            (queue_entry_index_temp >
                             INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH - 1)
                            ? queue_entry_index_temp - INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH
                            : queue_entry_index_temp;

                        if (encode_context_ptr
                                ->initial_rate_control_reorder_queue[queue_entry_index_temp2]
                                ->parent_pcs_wrapper_ptr != NULL) {
                            PictureParentControlSet *pcs =
                                (PictureParentControlSet *)(encode_context_ptr
                                                                ->initial_rate_control_reorder_queue
                                                                    [queue_entry_index_temp2]
                                                                ->parent_pcs_wrapper_ptr)
                                    ->object_ptr;
                            if (pcs->is_next_frame_intra)
                                break;
                        }

                        move_slide_window_flag = (EbBool)(
                            move_slide_window_flag &&
                            (encode_context_ptr
                                 ->initial_rate_control_reorder_queue[queue_entry_index_temp2]
                                 ->parent_pcs_wrapper_ptr != NULL));
                        if (encode_context_ptr
                                ->initial_rate_control_reorder_queue[queue_entry_index_temp2]
                                ->parent_pcs_wrapper_ptr != NULL) {
                            // check if it is the last frame. If we have reached the last frame, we would output the buffered frames in the Queue.
                            end_of_sequence_flag = ((PictureParentControlSet
                                                         *)(encode_context_ptr
                                                                ->initial_rate_control_reorder_queue
                                                                    [queue_entry_index_temp2]
                                                                ->parent_pcs_wrapper_ptr)
                                                        ->object_ptr)
                                                       ->end_of_sequence_flag;
                        } else
                            end_of_sequence_flag = EB_FALSE;
                        queue_entry_index_temp++;
                    }

                    if (move_slide_window_flag) {
                        //get a new entry spot
                        InitialRateControlReorderEntry *queue_entry_ptr =
                            encode_context_ptr->initial_rate_control_reorder_queue
                                [encode_context_ptr->initial_rate_control_reorder_queue_head_index];
                        pcs_ptr =
                            ((PictureParentControlSet *)(queue_entry_ptr->parent_pcs_wrapper_ptr)
                                 ->object_ptr);
                        scs_ptr = (SequenceControlSet *)pcs_ptr->scs_wrapper_ptr->object_ptr;
                        // overlay picture was not added to the queue. For the alt_ref picture with an overlay picture, it loops on both alt ref and overlay pictures
                        uint8_t has_overlay = pcs_ptr->is_alt_ref ? 1 : 0;
                        for (uint8_t loop_index = 0; loop_index <= has_overlay; loop_index++) {
                            if (loop_index)
                                pcs_ptr = pcs_ptr->overlay_ppcs_ptr;
                            pcs_ptr->frames_in_sw = frames_in_sw;
                            queue_entry_index_temp =
                                encode_context_ptr->initial_rate_control_reorder_queue_head_index;
                            end_of_sequence_flag = EB_FALSE;
                            // find the frames_in_interval for the peroid I frames
                            while (!end_of_sequence_flag &&
                                   queue_entry_index_temp <=
                                       encode_context_ptr
                                               ->initial_rate_control_reorder_queue_head_index +
                                           scs_ptr->static_config.look_ahead_distance) {
                                uint32_t queue_entry_index_temp2 =
                                    (queue_entry_index_temp >
                                     INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH - 1)
                                    ? queue_entry_index_temp -
                                        INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH
                                    : queue_entry_index_temp;
                                //exit if we hit a non valid entry
                                if (encode_context_ptr
                                        ->initial_rate_control_reorder_queue
                                            [queue_entry_index_temp2]
                                        ->parent_pcs_wrapper_ptr == NULL)
                                    break;
                                PictureParentControlSet *pcs_ptr_temp =
                                    ((PictureParentControlSet
                                          *)(encode_context_ptr
                                                 ->initial_rate_control_reorder_queue
                                                     [queue_entry_index_temp2]
                                                 ->parent_pcs_wrapper_ptr)
                                         ->object_ptr);
                                if (scs_ptr->intra_period_length != -1) {
                                    if (pcs_ptr->picture_number %
                                            ((scs_ptr->intra_period_length + 1)) ==
                                        0) {
                                        pcs_ptr->frames_in_interval[pcs_ptr_temp
                                                                        ->temporal_layer_index]++;
                                        if (pcs_ptr_temp->scene_change_flag)
                                            pcs_ptr->scene_change_in_gop = EB_TRUE;
                                    }
                                }

                                pcs_ptr_temp->historgram_life_count++;
                                end_of_sequence_flag = pcs_ptr_temp->end_of_sequence_flag;
                                queue_entry_index_temp++;
                            }

                            if ((scs_ptr->static_config.look_ahead_distance != 0) &&
                                (frames_in_sw < (scs_ptr->static_config.look_ahead_distance + 1)))
                                pcs_ptr->end_of_sequence_region = EB_TRUE;
                            else
                                pcs_ptr->end_of_sequence_region = EB_FALSE;
                            if (scs_ptr->static_config.rate_control_mode &&
                                !use_input_stat(scs_ptr) && !scs_ptr->lap_enabled)
                                // Determine offset from the Head Ptr for HLRC histogram queue and set the life count
                                if (scs_ptr->static_config.look_ahead_distance != 0)
                                    // Update Histogram Queue Entry Life count
                                    update_histogram_queue_entry(
                                        scs_ptr, encode_context_ptr, pcs_ptr, frames_in_sw);
                            // BACKGROUND ENHANCEMENT Part II
                            if (!pcs_ptr->end_of_sequence_flag &&
                                scs_ptr->static_config.look_ahead_distance != 0) {
                                // Update BEA information based on Lookahead information
                                update_bea_info_over_time(encode_context_ptr, pcs_ptr);
                            } else {
                                // Reset zzCost information to default When there's no lookahead available
                                init_zz_cost_info(pcs_ptr);
                            }

                            // Use the temporal layer 0 is_sb_motion_field_non_uniform array for all the other layer pictures in the mini GOP
                            if (!pcs_ptr->end_of_sequence_flag &&
                                scs_ptr->static_config.look_ahead_distance != 0) {
                                // Updat uniformly moving SBs based on Collocated SBs in LookAhead window
                                update_motion_field_uniformity_over_time(
                                    encode_context_ptr, scs_ptr, pcs_ptr);
                            }
                            if (pcs_ptr->is_used_as_reference_flag) {
                                // Get Empty Reference Picture Object
                                svt_get_empty_object(
                                    scs_ptr->encode_context_ptr->reference_picture_pool_fifo_ptr,
                                    &reference_picture_wrapper_ptr);
                                if (loop_index) {
                                    pcs_ptr->reference_picture_wrapper_ptr =
                                        reference_picture_wrapper_ptr;
                                    // Give the new Reference a nominal live_count of 1
                                    svt_object_inc_live_count(
                                        pcs_ptr->reference_picture_wrapper_ptr, 1);
                                } else {
                                    ((PictureParentControlSet
                                          *)(queue_entry_ptr->parent_pcs_wrapper_ptr->object_ptr))
                                        ->reference_picture_wrapper_ptr =
                                        reference_picture_wrapper_ptr;
                                    // Give the new Reference a nominal live_count of 1
                                    svt_object_inc_live_count(
                                        ((PictureParentControlSet *)(queue_entry_ptr
                                                                         ->parent_pcs_wrapper_ptr
                                                                         ->object_ptr))
                                            ->reference_picture_wrapper_ptr,
                                        1);
                                }
                            } else {
                                pcs_ptr->reference_picture_wrapper_ptr = NULL;
                            }
                            // Get Empty Results Object
                            svt_get_empty_object(
                                context_ptr->initialrate_control_results_output_fifo_ptr,
                                &out_results_wrapper_ptr);

                            InitialRateControlResults *out_results_ptr =
                                (InitialRateControlResults *)out_results_wrapper_ptr->object_ptr;

                            if (loop_index)
                                out_results_ptr->pcs_wrapper_ptr = pcs_ptr->p_pcs_wrapper_ptr;
                            else
                                out_results_ptr->pcs_wrapper_ptr =
                                    queue_entry_ptr->parent_pcs_wrapper_ptr;
                            if (scs_ptr->static_config.look_ahead_distance != 0 &&
                                scs_ptr->static_config.enable_tpl_la &&
                                ((has_overlay == 0 && loop_index == 0) ||
                                 (has_overlay == 1 && loop_index == 1))) {
                                // Release Pa Ref pictures when not needed
                                release_pa_reference_objects(scs_ptr, pcs_ptr);
                                //loop_index ? pcs_ptr : queueEntryPtr);
                            }
                            // Post the Full Results Object
                            svt_post_full_object(out_results_wrapper_ptr);
                        }
                        // Reset the Reorder Queue Entry
                        queue_entry_ptr->picture_number +=
                            INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH;
                        queue_entry_ptr->parent_pcs_wrapper_ptr = (EbObjectWrapper *)NULL;

                        // Increment the Reorder Queue head Ptr
                        encode_context_ptr->initial_rate_control_reorder_queue_head_index =
                            (encode_context_ptr->initial_rate_control_reorder_queue_head_index ==
                             INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH - 1)
                            ? 0
                            : encode_context_ptr->initial_rate_control_reorder_queue_head_index + 1;

                        queue_entry_ptr =
                            encode_context_ptr->initial_rate_control_reorder_queue
                                [encode_context_ptr->initial_rate_control_reorder_queue_head_index];
                    }
                }
            }
        }
#if FTR_TPL_TR
        }
#endif
        // Release the Input Results
        svt_release_object(in_results_wrapper_ptr);
    }
    return NULL;
}
