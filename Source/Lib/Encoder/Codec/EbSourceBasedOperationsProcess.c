/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbPictureControlSet.h"
#include "EbSequenceControlSet.h"

#include "EbSourceBasedOperationsProcess.h"
#include "EbInitialRateControlResults.h"
#include "EbPictureDemuxResults.h"
#ifdef ARCH_X86
#include "emmintrin.h"
#endif
#include "EbEncHandle.h"
#include "EbUtility.h"

/**************************************
 * Context
 **************************************/

typedef struct SourceBasedOperationsContext {
    EbDctor dctor;
    EbFifo *initial_rate_control_results_input_fifo_ptr;
    EbFifo *picture_demux_results_output_fifo_ptr;
    // local zz cost array
    uint32_t complete_sb_count;
    uint8_t *y_mean_ptr;
    uint8_t *cr_mean_ptr;
    uint8_t *cb_mean_ptr;
} SourceBasedOperationsContext;

static void source_based_operations_context_dctor(EbPtr p) {
    EbThreadContext *             thread_context_ptr = (EbThreadContext *)p;
    SourceBasedOperationsContext *obj = (SourceBasedOperationsContext *)thread_context_ptr->priv;
    EB_FREE_ARRAY(obj);
}

/************************************************
* Source Based Operation Context Constructor
************************************************/
EbErrorType source_based_operations_context_ctor(EbThreadContext *  thread_context_ptr,
                                                 const EbEncHandle *enc_handle_ptr, int index) {
    SourceBasedOperationsContext *context_ptr;
    EB_CALLOC_ARRAY(context_ptr, 1);
    thread_context_ptr->priv  = context_ptr;
    thread_context_ptr->dctor = source_based_operations_context_dctor;

    context_ptr->initial_rate_control_results_input_fifo_ptr = eb_system_resource_get_consumer_fifo(
        enc_handle_ptr->initial_rate_control_results_resource_ptr, index);
    context_ptr->picture_demux_results_output_fifo_ptr = eb_system_resource_get_producer_fifo(
        enc_handle_ptr->picture_demux_results_resource_ptr, index);
    return EB_ErrorNone;
}

/***************************************************
* Derives BEA statistics and set activity flags
***************************************************/
void derive_picture_activity_statistics(PictureParentControlSet *pcs_ptr)

{
    uint64_t non_moving_index_min = ~0u;
    uint64_t non_moving_index_max = 0;
    uint64_t non_moving_index_sum = 0;
    uint32_t complete_sb_count    = 0;
    uint32_t non_moving_sb_count  = 0;
    uint32_t sb_total_count       = pcs_ptr->sb_total_count;
    uint32_t tot_nmv_idx          = 0;

    uint32_t sb_index;
    for (sb_index = 0; sb_index < sb_total_count; ++sb_index) {
        SbParams *sb_params = &pcs_ptr->sb_params_array[sb_index];
        if (sb_params->is_complete_sb) {
            non_moving_index_min = pcs_ptr->non_moving_index_array[sb_index] < non_moving_index_min
                                       ? pcs_ptr->non_moving_index_array[sb_index]
                                       : non_moving_index_min;

            non_moving_index_max = pcs_ptr->non_moving_index_array[sb_index] > non_moving_index_max
                                       ? pcs_ptr->non_moving_index_array[sb_index]
                                       : non_moving_index_max;
            if (pcs_ptr->non_moving_index_array[sb_index] < NON_MOVING_SCORE_1)
                non_moving_sb_count++;
            complete_sb_count++;

            non_moving_index_sum += pcs_ptr->non_moving_index_array[sb_index];

            if (pcs_ptr->non_moving_index_array[sb_index] < NON_MOVING_SCORE_1) tot_nmv_idx++;
        }
    }

    if (complete_sb_count > 0) {
        pcs_ptr->non_moving_index_average = (uint16_t)(non_moving_index_sum / complete_sb_count);
        pcs_ptr->kf_zeromotion_pct        = (non_moving_sb_count * 100) / complete_sb_count;
    }
    pcs_ptr->non_moving_index_min_distance = (uint16_t)(
        ABS((int32_t)(pcs_ptr->non_moving_index_average) - (int32_t)non_moving_index_min));
    pcs_ptr->non_moving_index_max_distance = (uint16_t)(
        ABS((int32_t)(pcs_ptr->non_moving_index_average) - (int32_t)non_moving_index_max));
    return;
}

/************************************************
 * Source Based Operations Kernel
 * Source-based operations process involves a number of analysis algorithms
 * to identify spatiotemporal characteristics of the input pictures.
 ************************************************/
void *source_based_operations_kernel(void *input_ptr) {
    EbThreadContext *             thread_context_ptr = (EbThreadContext *)input_ptr;
    SourceBasedOperationsContext *context_ptr =
        (SourceBasedOperationsContext *)thread_context_ptr->priv;
    PictureParentControlSet *  pcs_ptr;
    EbObjectWrapper *          in_results_wrapper_ptr;
    InitialRateControlResults *in_results_ptr;
    EbObjectWrapper *          out_results_wrapper_ptr;
    PictureDemuxResults *      out_results_ptr;

    for (;;) {
        // Get Input Full Object
        EB_GET_FULL_OBJECT(context_ptr->initial_rate_control_results_input_fifo_ptr,
                           &in_results_wrapper_ptr);

        in_results_ptr = (InitialRateControlResults *)in_results_wrapper_ptr->object_ptr;
        pcs_ptr        = (PictureParentControlSet *)in_results_ptr->pcs_wrapper_ptr->object_ptr;
        context_ptr->complete_sb_count             = 0;
        uint32_t sb_total_count                    = pcs_ptr->sb_total_count;
        uint32_t sb_index;

        /***********************************************SB-based operations************************************************************/
        for (sb_index = 0; sb_index < sb_total_count; ++sb_index) {
            SbParams *sb_params      = &pcs_ptr->sb_params_array[sb_index];
            EbBool    is_complete_sb = sb_params->is_complete_sb;
            uint8_t * y_mean_ptr     = pcs_ptr->y_mean[sb_index];
#ifdef ARCH_X86
            _mm_prefetch((const char *)y_mean_ptr, _MM_HINT_T0);
#endif
            uint8_t *cr_mean_ptr = pcs_ptr->cr_mean[sb_index];
            uint8_t *cb_mean_ptr = pcs_ptr->cb_mean[sb_index];
#ifdef ARCH_X86
            _mm_prefetch((const char *)cr_mean_ptr, _MM_HINT_T0);
            _mm_prefetch((const char *)cb_mean_ptr, _MM_HINT_T0);
#endif
            context_ptr->y_mean_ptr  = y_mean_ptr;
            context_ptr->cr_mean_ptr = cr_mean_ptr;
            context_ptr->cb_mean_ptr = cb_mean_ptr;

            if (is_complete_sb) { context_ptr->complete_sb_count++; }
        }
        /*********************************************Picture-based operations**********************************************************/

        // Activity statistics derivation
        derive_picture_activity_statistics(pcs_ptr);

        // Get Empty Results Object
        eb_get_empty_object(context_ptr->picture_demux_results_output_fifo_ptr,
                            &out_results_wrapper_ptr);

        out_results_ptr = (PictureDemuxResults *)out_results_wrapper_ptr->object_ptr;
        out_results_ptr->pcs_wrapper_ptr = in_results_ptr->pcs_wrapper_ptr;
        out_results_ptr->picture_type    = EB_PIC_INPUT;

        // Release the Input Results
        eb_release_object(in_results_wrapper_ptr);

        // Post the Full Results Object
        eb_post_full_object(out_results_wrapper_ptr);
    }
    return NULL;
}
