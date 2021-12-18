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
#include "EbLog.h"
#include "EbPictureDecisionProcess.h"
#if FTR_1PAS_VBR
#include "firstpass.h"
#endif
/**************************************
 * Context
 **************************************/
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


typedef struct InitialRateControlContext {
    EbFifo *motion_estimation_results_input_fifo_ptr;
    EbFifo *initialrate_control_results_output_fifo_ptr;
    LadQueue         *lad_queue;

} InitialRateControlContext;

/**************************************
* Macros
**************************************/
static void initial_rate_control_context_dctor(EbPtr p) {
    EbThreadContext *          thread_context_ptr = (EbThreadContext *)p;
    InitialRateControlContext *obj = (InitialRateControlContext *)thread_context_ptr->priv;

    EB_DELETE_PTR_ARRAY(obj->lad_queue->cir_buf, REFERENCE_QUEUE_MAX_DEPTH);
    EB_FREE(obj->lad_queue);
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

    EB_MALLOC(context_ptr->lad_queue, sizeof(LadQueue));

    EB_ALLOC_PTR_ARRAY(context_ptr->lad_queue->cir_buf, REFERENCE_QUEUE_MAX_DEPTH);
    for (uint32_t picture_index = 0; picture_index < REFERENCE_QUEUE_MAX_DEPTH; ++picture_index) {
        EB_NEW(context_ptr->lad_queue->cir_buf[picture_index], lad_queue_entry_ctor);
    }
    context_ptr->lad_queue->head = 0;
    context_ptr->lad_queue->tail = 0;

    return EB_ErrorNone;
}
#if !CLN_NON_MOVING_CNT
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
#endif

void svt_av1_build_quantizer(AomBitDepth bit_depth, int32_t y_dc_delta_q, int32_t u_dc_delta_q,
                             int32_t u_ac_delta_q, int32_t v_dc_delta_q, int32_t v_ac_delta_q,
                             Quants *const quants, Dequants *const deq);



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
void irc_send_picture_out(InitialRateControlContext *ctx, PictureParentControlSet *pcs, EbBool superres_recode)
{
    EbObjectWrapper *out_results_wrapper_ptr;
    // Get Empty Results Object
    svt_get_empty_object(ctx->initialrate_control_results_output_fifo_ptr,
           &out_results_wrapper_ptr);
    InitialRateControlResults *out_results_ptr =
          (InitialRateControlResults *)out_results_wrapper_ptr->object_ptr;
    //SVT_LOG("iRC Out:%lld\n",pcs->picture_number);
    out_results_ptr->pcs_wrapper_ptr = pcs->p_pcs_wrapper_ptr;
    out_results_ptr->superres_recode = superres_recode;
    svt_post_full_object(out_results_wrapper_ptr);
}
uint8_t is_frame_already_exists(PictureParentControlSet *pcs, uint32_t end_index, uint64_t pic_num) {
    for (uint32_t i = 0; i < end_index; i++)
        if (pcs->tpl_group[i]->picture_number == pic_num)
            return 1;
    return 0;
}

/* decide to inject a frame into the tpl group*/
uint8_t  inj_to_tpl_group( PictureParentControlSet* pcs)
{
    uint8_t inj = 0;
    if (pcs->hierarchical_levels != 4) {
        if (pcs->temporal_layer_index < pcs->hierarchical_levels)
            inj = 1;
        else
            inj = 0;
    }
    else {
#if CLN_MERGE_MRP_SIG
#if CLN_MRP_ENC_CONFIG
        if (pcs->scs_ptr->mrp_ctrls.referencing_scheme == 1) {
#else
        if (pcs->scs_ptr->static_config.mrp_ctrls.referencing_scheme == 1) {
#endif
#else
        if (pcs->scs_ptr->mrp_init_level == 1) {
#endif
            if (pcs->temporal_layer_index < 4)
                inj = 1;
            else if (pcs->is_used_as_reference_flag && (pcs->pic_index == 6 || pcs->pic_index == 10))//TODO: could be removed once TPL r0 adapts dyncamically  to TPL group size
                inj = 1;
            else
                inj = 0;

        }
        else {
            if (pcs->temporal_layer_index < 4)
                inj = 1;
            else if (pcs->is_used_as_reference_flag && (pcs->pic_index == 0 || pcs->pic_index == 8))
                inj = 1;
            else
                inj = 0;

        }
    }

    return inj;
}
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
    if (pcs->tpl_group_size != pcs->ext_group_size)
        printf("asfafs");
    for (uint32_t i = 0; i < pcs->ext_group_size; i++)
        if (pcs->tpl_group[i] != pcs->ext_group[i])
            printf("asfafs");


    pcs->tpl_group_size = pcs->ext_group_size;

    for (uint32_t i = 0; i < pcs->ext_group_size; i++)
        pcs->tpl_group[i] = pcs->ext_group[i];
#endif


    //new tpl group needs to stop at the second I
    pcs->ntpl_group_size = 0;
    uint8_t is_gop_end = 0;
    int64_t last_intra_mg_id;
#if FTR_1PAS_VBR
#if FTR_LAD_INPUT
#if FIX_DATA_RACE_2PASS
    uint32_t mg_size;
    if (pcs->scs_ptr->static_config.enable_adaptive_mini_gop == 0) {
        mg_size = 1 << pcs->scs_ptr->static_config.hierarchical_levels;
    }
    else {
#if FIX_DG
        mg_size = 1 << pcs->scs_ptr->static_config.max_heirachical_level;
#else
        mg_size = 1 << pcs->hierarchical_levels;
#endif
    }
#else
    uint32_t mg_size                = 1 << pcs->scs_ptr->static_config.hierarchical_levels;
#endif
    uint32_t limited_tpl_group_size = pcs->slice_type == I_SLICE
        ? MIN(1 + (pcs->scs_ptr->tpl_lad_mg + 1) * mg_size, pcs->ext_group_size)
        : MIN((pcs->scs_ptr->tpl_lad_mg + 1) * mg_size, pcs->ext_group_size);
#else
    uint32_t limited_tpl_group_size = pcs->slice_type == I_SLICE ? MIN(33, pcs->ext_group_size) : MIN(32, pcs->ext_group_size);
#endif
    for (uint32_t i = 0; i < limited_tpl_group_size; i++) {
#else
    for (uint32_t i = 0; i < pcs->ext_group_size; i++) {
#endif
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
#if !OPT_COMBINE_TPL_FOR_LAD
    SequenceControlSet *scs_ptr = (SequenceControlSet *)pcs->scs_wrapper_ptr->object_ptr;
    if (scs_ptr->lad_mg)
#endif
    {
        pcs->tpl_group_size = pcs->ntpl_group_size;
        memset(pcs->tpl_valid_pic, 0, MAX_TPL_EXT_GROUP_SIZE * sizeof(uint8_t));
        pcs->tpl_valid_pic[0] = 1;
        pcs->used_tpl_frame_num = 0;
        for (uint32_t i = 0; i < pcs->ntpl_group_size; i++) {
            pcs->tpl_group[i] = pcs->ntpl_group[i];
            // Check wether the i-th pic already exists in the tpl group
            if (!is_frame_already_exists(pcs, i, pcs->tpl_group[i]->picture_number)) {
                // Discard non-ref pic from the tpl group
                uint8_t inject_frame = inj_to_tpl_group(pcs->tpl_group[i]);
                if (inject_frame) {
                    if (pcs->slice_type != I_SLICE) {
                        // Discard low important pictures from tpl group
#if OPT5_TPL_REDUCE_PIC
                        if (pcs->tpl_ctrls.tpl_opt_flag && (pcs->tpl_ctrls.reduced_tpl_group >= 0)) {
#else
                        if (pcs->tpl_ctrls.tpl_opt_flag && pcs->tpl_ctrls.reduced_tpl_group) {
#endif
#if !TPL_CLEANUP
                            if (pcs->tpl_group[i]->temporal_layer_index <= pcs->tpl_ctrls.reduced_tpl_group + (pcs->hierarchical_levels == 5 ? 1 : 0)) {
#else
                            if (pcs->tpl_group[i]->temporal_layer_index <= pcs->tpl_ctrls.reduced_tpl_group) {
#endif
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
    }

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
#if FTR_1PAS_VBR
                        if (tmp_pcs->end_of_sequence_flag)
                            head_pcs->end_of_sequence_region = EB_TRUE;
#endif
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
#if FTR_1PAS_VBR
#if CLN_ENC_CONFIG_SIG
            if(head_pcs->scs_ptr->static_config.pass == ENC_MIDDLE_PASS || head_pcs->scs_ptr->static_config.pass == ENC_LAST_PASS || head_pcs->scs_ptr->lap_enabled) {
#else
            if (head_pcs->scs_ptr->lap_enabled || use_input_stat(head_pcs->scs_ptr)) {
#endif
                head_pcs->stats_in_offset = head_pcs->decode_order;
#if CLN_ENC_CONFIG_SIG
                head_pcs->stats_in_end_offset = head_pcs->ext_group_size && !(head_pcs->scs_ptr->static_config.pass == ENC_MIDDLE_PASS || head_pcs->scs_ptr->static_config.pass == ENC_LAST_PASS) ?
#else
                head_pcs->stats_in_end_offset = head_pcs->ext_group_size && !use_input_stat(head_pcs->scs_ptr) ?
#endif
                    MIN((uint64_t)(head_pcs->scs_ptr->twopass.stats_buf_ctx->stats_in_end_write - head_pcs->scs_ptr->twopass.stats_buf_ctx->stats_in_start),
                      head_pcs->stats_in_offset + (uint64_t)head_pcs->ext_group_size + 1) :
                    (uint64_t)(head_pcs->scs_ptr->twopass.stats_buf_ctx->stats_in_end_write - head_pcs->scs_ptr->twopass.stats_buf_ctx->stats_in_start);
                head_pcs->frames_in_sw = (int)(head_pcs->stats_in_end_offset - head_pcs->stats_in_offset);
#if FIX_VBR_R2R
                if (head_pcs->scs_ptr->enable_dec_order == 0 && head_pcs->scs_ptr->lap_enabled && head_pcs->temporal_layer_index == 0) {
#else
                if (head_pcs->scs_ptr->lap_enabled && head_pcs->temporal_layer_index == 0) {
#endif
                    for (uint64_t num_frames = head_pcs->stats_in_offset; num_frames < head_pcs->stats_in_end_offset; ++num_frames) {
                        FIRSTPASS_STATS *cur_frame =
                            head_pcs->scs_ptr->twopass.stats_buf_ctx->stats_in_start + num_frames;
                        if ((int64_t)cur_frame->frame > head_pcs->scs_ptr->twopass.stats_buf_ctx->last_frame_accumulated) {
                            svt_av1_accumulate_stats(head_pcs->scs_ptr->twopass.stats_buf_ctx->total_stats, cur_frame);
                            head_pcs->scs_ptr->twopass.stats_buf_ctx->last_frame_accumulated = (int64_t)cur_frame->frame;
                        }
                    }
                }
            }
#endif
            //take the picture out from iRc process
            irc_send_picture_out(ctx, head_pcs, EB_FALSE);
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
        pcs_ptr->me_segments_completion_count++;

        // If the picture is complete, proceed
        if (pcs_ptr->me_segments_completion_count == pcs_ptr->me_segments_total_count) {
            SequenceControlSet *scs_ptr = (SequenceControlSet *)
                                              pcs_ptr->scs_wrapper_ptr->object_ptr;

            if (in_results_ptr->task_type == TASK_SUPERRES_RE_ME) {
                // do necessary steps as normal routine
                {
                    // Release Pa Ref pictures when not needed
                    // Don't release if superres recode loop is actived (auto-dual or auto-all mode)
                    if (pcs_ptr->superres_total_recode_loop == 0) {  // QThreshold or auto-solo mode
                        if (scs_ptr->static_config.enable_tpl_la) {
                            for (uint32_t i = 0; i < pcs_ptr->tpl_group_size; i++) {
                                if (pcs_ptr->tpl_group[i]->slice_type == P_SLICE) {
                                    if (pcs_ptr->tpl_group[i]->ext_mg_id == pcs_ptr->ext_mg_id + 1)
                                        release_pa_reference_objects(scs_ptr, pcs_ptr->tpl_group[i]);
                                }
                                else {
                                    if (pcs_ptr->tpl_group[i]->ext_mg_id == pcs_ptr->ext_mg_id)
                                        release_pa_reference_objects(scs_ptr, pcs_ptr->tpl_group[i]);
                                }
                                if (pcs_ptr->tpl_group[i]->non_tf_input)
                                    EB_DELETE(pcs_ptr->tpl_group[i]->non_tf_input);
                            }
                        } else {
                            release_pa_reference_objects(scs_ptr, pcs_ptr);
                        }
                    }

                    /*In case Look-Ahead is zero there is no need to place pictures in the
                      re-order queue. this will cause an artificial delay since pictures come in dec-order*/
                    pcs_ptr->frames_in_sw = 0;
                    pcs_ptr->end_of_sequence_region = EB_FALSE;
#if !CLN_NON_MOVING_CNT
                    init_zz_cost_info(pcs_ptr);
#endif
                }

                // post to downstream process
                irc_send_picture_out(context_ptr, pcs_ptr, EB_TRUE);

                // Release the Input Results
                svt_release_object(in_results_wrapper_ptr);
                continue;
            }

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
            // perform tpl_la on unscaled frames only
            if (scs_ptr->static_config.enable_tpl_la && !pcs_ptr->frame_superres_enabled) {
                svt_set_cond_var(&pcs_ptr->me_ready, 1);
            }

            // Release Pa Ref pictures when not needed
            // Release Pa ref when
            //   1. TPL is OFF and
            //   2. super-res mode is NONE or FIXED or RANDOM.
            //     For other super-res modes, pa_ref_objs are needed in TASK_SUPERRES_RE_ME task
#if FIX_LAD_MG_0_HANG && !FIX_I51
            // Release Pa Ref if mg_lad is 0 and P slice (not belonging to any TPL group)
            uint8_t release_pa_ref = 0;
            if (scs_ptr->static_config.enable_tpl_la == 0 && scs_ptr->static_config.superres_mode <= SUPERRES_RANDOM)
                release_pa_ref = 1;
#if FTR_LAD_INPUT && !FIX_VBR_R2R
            else if (scs_ptr->tpl_lad_mg == 0 && pcs_ptr->slice_type == P_SLICE)
#else
            else if (scs_ptr->lad_mg == 0 && pcs_ptr->slice_type == P_SLICE)
#endif
                release_pa_ref = 1;

            if (release_pa_ref)
#else
            if (scs_ptr->static_config.enable_tpl_la == 0 && scs_ptr->static_config.superres_mode <= SUPERRES_RANDOM)
#endif
                release_pa_reference_objects(scs_ptr, pcs_ptr);

            /*In case Look-Ahead is zero there is no need to place pictures in the
              re-order queue. this will cause an artificial delay since pictures come in dec-order*/
            pcs_ptr->frames_in_sw           = 0;
            pcs_ptr->end_of_sequence_region = EB_FALSE;
#if !CLN_NON_MOVING_CNT
            init_zz_cost_info(pcs_ptr);
#endif

            push_to_lad_queue(pcs_ptr, context_ptr);
#if LAD_MG_PRINT
            print_lad_queue(context_ptr,0);
#endif
            uint8_t lad_queue_pass_thru = !(scs_ptr->static_config.enable_tpl_la && !pcs_ptr->frame_superres_enabled);
            process_lad_queue(context_ptr, lad_queue_pass_thru);

        }
        // Release the Input Results
        svt_release_object(in_results_wrapper_ptr);
    }
    return NULL;
}
