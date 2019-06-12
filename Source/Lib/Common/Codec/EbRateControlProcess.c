/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

/*
* Copyright (c) 2016, Alliance for Open Media. All rights reserved
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at www.aomedia.org/license/software. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at www.aomedia.org/license/patent.
*/
#include <stdlib.h>

#include "EbDefinitions.h"
#include "EbRateControlProcess.h"
#include "EbSystemResourceManager.h"
#include "EbSequenceControlSet.h"
#include "EbPictureControlSet.h"
#include "EbUtility.h"
#include "EbSvtAv1ErrorCodes.h"
#include "EbEntropyCoding.h"

#include "EbRateControlResults.h"
#include "EbRateControlTasks.h"
#include "RateControlModel.h"

#include "EbSegmentation.h"

// calculate the QP based on the QP scaling
uint32_t qp_scaling_calc(
    SequenceControlSet *sequence_control_set_ptr,
    EB_SLICE            slice_type,
    uint32_t            temporal_layer_index,
    uint32_t            base_qp);

/*****************************
* Internal Typedefs
*****************************/
void rate_control_layer_reset(
    RateControlLayerContext *rate_control_layer_ptr,
    PictureControlSet       *picture_control_set_ptr,
    RateControlContext      *rate_control_context_ptr,
    uint32_t                 picture_area_in_pixel,
    EbBool                  was_used)
{
    SequenceControlSet *sequence_control_set_ptr = (SequenceControlSet *)picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr;
    uint32_t            slice_num;
    uint32_t            temporal_layer_index;
    uint64_t            total_frame_in_interval;
    uint64_t            sum_bits_per_sw = 0;

    rate_control_layer_ptr->target_bit_rate = picture_control_set_ptr->parent_pcs_ptr->target_bit_rate*(uint64_t)rate_percentage_layer_array[sequence_control_set_ptr->static_config.hierarchical_levels][rate_control_layer_ptr->temporal_index] / 100;
    // update this based on temporal layers
    rate_control_layer_ptr->frame_rate = sequence_control_set_ptr->frame_rate;

    total_frame_in_interval = sequence_control_set_ptr->static_config.intra_period_length + 1;

    if (sequence_control_set_ptr->static_config.look_ahead_distance != 0 && sequence_control_set_ptr->intra_period_length != -1) {
        if (picture_control_set_ptr->picture_number % ((sequence_control_set_ptr->intra_period_length + 1)) == 0) {
            total_frame_in_interval = 0;
            for (temporal_layer_index = 0; temporal_layer_index < EB_MAX_TEMPORAL_LAYERS; temporal_layer_index++) {
                rate_control_context_ptr->frames_in_interval[temporal_layer_index] = picture_control_set_ptr->parent_pcs_ptr->frames_in_interval[temporal_layer_index];
                total_frame_in_interval += picture_control_set_ptr->parent_pcs_ptr->frames_in_interval[temporal_layer_index];
                sum_bits_per_sw += picture_control_set_ptr->parent_pcs_ptr->bits_per_sw_per_layer[temporal_layer_index];
            }
#if ADAPTIVE_PERCENTAGE
            rate_control_layer_ptr->target_bit_rate = picture_control_set_ptr->parent_pcs_ptr->target_bit_rate* picture_control_set_ptr->parent_pcs_ptr->bits_per_sw_per_layer[rate_control_layer_ptr->temporal_index] / sum_bits_per_sw;
#endif
        }
    }

    if (sequence_control_set_ptr->static_config.intra_period_length != -1)
        rate_control_layer_ptr->frame_rate = sequence_control_set_ptr->frame_rate * rate_control_context_ptr->frames_in_interval[rate_control_layer_ptr->temporal_index] / total_frame_in_interval;
    else {
        switch (picture_control_set_ptr->parent_pcs_ptr->hierarchical_levels) {
        case 0:
            break;
        case 1:
            if (sequence_control_set_ptr->static_config.intra_period_length == -1)
                rate_control_layer_ptr->frame_rate = rate_control_layer_ptr->frame_rate >> 1;
            break;
        case 2:
            if (rate_control_layer_ptr->temporal_index == 0)
                rate_control_layer_ptr->frame_rate = rate_control_layer_ptr->frame_rate >> 2;
            else
                rate_control_layer_ptr->frame_rate = rate_control_layer_ptr->frame_rate >> (3 - rate_control_layer_ptr->temporal_index);
            break;
        case 3:
            if (rate_control_layer_ptr->temporal_index == 0)
                rate_control_layer_ptr->frame_rate = rate_control_layer_ptr->frame_rate >> 3;
            else
                rate_control_layer_ptr->frame_rate = rate_control_layer_ptr->frame_rate >> (4 - rate_control_layer_ptr->temporal_index);
            break;
        case 4:
            if (rate_control_layer_ptr->temporal_index == 0)
                rate_control_layer_ptr->frame_rate = rate_control_layer_ptr->frame_rate >> 4;
            else
                rate_control_layer_ptr->frame_rate = rate_control_layer_ptr->frame_rate >> (5 - rate_control_layer_ptr->temporal_index);
            break;
        case 5:
            if (rate_control_layer_ptr->temporal_index == 0)
                rate_control_layer_ptr->frame_rate = rate_control_layer_ptr->frame_rate >> 5;
            else
                rate_control_layer_ptr->frame_rate = rate_control_layer_ptr->frame_rate >> (6 - rate_control_layer_ptr->temporal_index);
            break;

        default:
            break;
        }
    }

    rate_control_layer_ptr->coeff_averaging_weight1 = 5;

    rate_control_layer_ptr->coeff_averaging_weight2 = 16 - rate_control_layer_ptr->coeff_averaging_weight1;
    if (rate_control_layer_ptr->frame_rate == 0) { // no frame in that layer
        rate_control_layer_ptr->frame_rate = 1 << RC_PRECISION;
    }
    rate_control_layer_ptr->channel_bit_rate = (((rate_control_layer_ptr->target_bit_rate << (2 * RC_PRECISION)) / rate_control_layer_ptr->frame_rate) + RC_PRECISION_OFFSET) >> RC_PRECISION;
    rate_control_layer_ptr->channel_bit_rate = (uint64_t)MAX((int64_t)1, (int64_t)rate_control_layer_ptr->channel_bit_rate);
    rate_control_layer_ptr->ec_bit_constraint = rate_control_layer_ptr->channel_bit_rate;

    // This is only for the initial frame, because the feedback is from packetization now and all of these are considered
    // considering the bits for slice header
    // *Note - only one-slice-per picture is supported for UHD
    slice_num = 1;

    rate_control_layer_ptr->ec_bit_constraint -= SLICE_HEADER_BITS_NUM * slice_num;

    rate_control_layer_ptr->ec_bit_constraint = MAX(1, rate_control_layer_ptr->ec_bit_constraint);

    rate_control_layer_ptr->previous_bit_constraint = rate_control_layer_ptr->channel_bit_rate;
    rate_control_layer_ptr->bit_constraint = rate_control_layer_ptr->channel_bit_rate;
    rate_control_layer_ptr->dif_total_and_ec_bits = 0;

    rate_control_layer_ptr->frame_same_distortion_min_qp_count = 0;
    rate_control_layer_ptr->max_qp = picture_control_set_ptr->picture_qp;

    rate_control_layer_ptr->alpha = 1 << (RC_PRECISION - 1);
    {
        if (!was_used) {
            rate_control_layer_ptr->same_distortion_count = 0;
            rate_control_layer_ptr->k_coeff = 3 << RC_PRECISION;
            rate_control_layer_ptr->previous_k_coeff = 3 << RC_PRECISION;
            rate_control_layer_ptr->c_coeff = (rate_control_layer_ptr->channel_bit_rate << (2 * RC_PRECISION)) / picture_area_in_pixel / CCOEFF_INIT_FACT;
            rate_control_layer_ptr->previous_c_coeff = (rate_control_layer_ptr->channel_bit_rate << (2 * RC_PRECISION)) / picture_area_in_pixel / CCOEFF_INIT_FACT;
            // These are for handling Pred structure 2, when for higher temporal layer, frames can arrive in different orders
            // They should be modifed in a way that gets these from previous layers
            rate_control_layer_ptr->previous_frame_qp = 32;
            rate_control_layer_ptr->previous_frame_bit_actual = 1200;
            rate_control_layer_ptr->previous_framequantized_coeff_bit_actual = 1000;
            rate_control_layer_ptr->previous_frame_distortion_me = 10000000;
            rate_control_layer_ptr->previous_frame_qp = picture_control_set_ptr->picture_qp;
            rate_control_layer_ptr->delta_qp_fraction = 0;
            rate_control_layer_ptr->previous_frame_average_qp = picture_control_set_ptr->picture_qp;
            rate_control_layer_ptr->previous_calculated_frame_qp = picture_control_set_ptr->picture_qp;
            rate_control_layer_ptr->calculated_frame_qp = picture_control_set_ptr->picture_qp;
            rate_control_layer_ptr->critical_states = 0;
        }
        else {
            rate_control_layer_ptr->same_distortion_count = 0;
            rate_control_layer_ptr->critical_states = 0;
        }
    }
}

void rate_control_layer_reset_part2(
    RateControlContext      *context_ptr,
    RateControlLayerContext *rate_control_layer_ptr,
    PictureControlSet       *picture_control_set_ptr)
{
    // update this based on temporal layers
    rate_control_layer_ptr->max_qp = (uint32_t)CLIP3(0, 63, (int32_t)context_ptr->qp_scaling_map[rate_control_layer_ptr->temporal_index][picture_control_set_ptr->picture_qp]);
    // These are for handling Pred structure 2, when for higher temporal layer, frames can arrive in different orders
    // They should be modifed in a way that gets these from previous layers
    rate_control_layer_ptr->previous_frame_qp = rate_control_layer_ptr->max_qp;
    rate_control_layer_ptr->previous_frame_average_qp = rate_control_layer_ptr->max_qp;
    rate_control_layer_ptr->previous_calculated_frame_qp = rate_control_layer_ptr->max_qp;
    rate_control_layer_ptr->calculated_frame_qp = rate_control_layer_ptr->max_qp;
}

EbErrorType high_level_rate_control_context_ctor(
    HighLevelRateControlContext *entry_ptr) {
    (void)entry_ptr;

    return EB_ErrorNone;
}

EbErrorType rate_control_layer_context_ctor(
    RateControlLayerContext *entry_ptr) {

    entry_ptr->first_frame = 1;
    entry_ptr->first_non_intra_frame = 1;

    return EB_ErrorNone;
}

void rate_control_interval_param_context_dctor(EbPtr p)
{
    RateControlIntervalParamContext* obj = (RateControlIntervalParamContext*)p;
    EB_DELETE_PTR_ARRAY(obj->rate_control_layer_array, EB_MAX_TEMPORAL_LAYERS);
}

EbErrorType rate_control_interval_param_context_ctor(
    RateControlIntervalParamContext *entry_ptr) {
    uint32_t temporal_index;

    entry_ptr->dctor = rate_control_interval_param_context_dctor;

    EB_ALLOC_PTR_ARRAY(entry_ptr->rate_control_layer_array, EB_MAX_TEMPORAL_LAYERS);

    for (temporal_index = 0; temporal_index < EB_MAX_TEMPORAL_LAYERS; temporal_index++) {
        EB_NEW(
            entry_ptr->rate_control_layer_array[temporal_index],
            rate_control_layer_context_ctor);
        entry_ptr->rate_control_layer_array[temporal_index]->temporal_index = temporal_index;
        entry_ptr->rate_control_layer_array[temporal_index]->frame_rate = 1 << RC_PRECISION;
    }

    return EB_ErrorNone;
}

EbErrorType rate_control_coded_frames_stats_context_ctor(
    CodedFramesStatsEntry  *entry_ptr,
    uint64_t                picture_number) {

    entry_ptr->picture_number = picture_number;
    entry_ptr->frame_total_bit_actual = -1;

    return EB_ErrorNone;
}

void rate_control_context_dctor(EbPtr p)
{
    RateControlContext* obj = (RateControlContext*)p;
#if OVERSHOOT_STAT_PRINT
    EB_DELETE_PTR_ARRAY(obj->coded_frames_stat_queue, CODED_FRAMES_STAT_QUEUE_MAX_DEPTH);
#endif
    EB_DELETE_PTR_ARRAY(obj->rate_control_param_queue, PARALLEL_GOP_MAX_NUMBER);
    EB_DELETE(obj->high_level_rate_control_ptr);
    EB_DELETE(obj->rc_model_ptr);

}

EbErrorType rate_control_context_ctor(
    RateControlContext *context_ptr,
    EbFifo             *rate_control_input_tasks_fifo_ptr,
    EbFifo             *rate_control_output_results_fifo_ptr,
    int32_t             intra_period)
{
    uint32_t interval_index;

#if OVERSHOOT_STAT_PRINT
    uint32_t picture_index;
#endif

    context_ptr->dctor = rate_control_context_dctor;
    context_ptr->rate_control_input_tasks_fifo_ptr = rate_control_input_tasks_fifo_ptr;
    context_ptr->rate_control_output_results_fifo_ptr = rate_control_output_results_fifo_ptr;

    // High level RC
    EB_NEW(
        context_ptr->high_level_rate_control_ptr,
        high_level_rate_control_context_ctor);

    EB_NEW(context_ptr->rc_model_ptr, rate_control_model_ctor);

    EB_ALLOC_PTR_ARRAY(context_ptr->rate_control_param_queue, PARALLEL_GOP_MAX_NUMBER);

    for (interval_index = 0; interval_index < PARALLEL_GOP_MAX_NUMBER; interval_index++) {
        EB_NEW(
            context_ptr->rate_control_param_queue[interval_index],
            rate_control_interval_param_context_ctor);
        context_ptr->rate_control_param_queue[interval_index]->first_poc = (interval_index*(uint32_t)(intra_period + 1));
        context_ptr->rate_control_param_queue[interval_index]->last_poc = ((interval_index + 1)*(uint32_t)(intra_period + 1)) - 1;
    }

#if OVERSHOOT_STAT_PRINT
    EB_ALLOC_PTR_ARRAY(context_ptr->coded_frames_stat_queue, CODED_FRAMES_STAT_QUEUE_MAX_DEPTH);

    for (picture_index = 0; picture_index < CODED_FRAMES_STAT_QUEUE_MAX_DEPTH; ++picture_index) {
        EB_NEW(context_ptr->coded_frames_stat_queue[picture_index],
            rate_control_coded_frames_stats_context_ctor,
            picture_index);
    }
    context_ptr->min_bit_actual_per_gop = 0xfffffffffffff;
#endif
    context_ptr->intra_coef_rate = 4;

    return EB_ErrorNone;
}
uint64_t predict_bits(
    EncodeContext                 *encode_context_ptr,
    HlRateControlHistogramEntry   *hl_rate_control_histogram_ptr_temp,
    uint32_t                         qp,
    uint32_t                         area_in_pixel)
{
    uint64_t total_bits = 0;

    if (hl_rate_control_histogram_ptr_temp->is_coded) {
        // If the frame is already coded, use the actual number of bits
        total_bits = hl_rate_control_histogram_ptr_temp->total_num_bits_coded;
    }
    else {
        RateControlTables     *rate_control_tables_ptr = &encode_context_ptr->rate_control_tables_array[qp];
        EbBitNumber             *sad_bits_array_ptr = rate_control_tables_ptr->sad_bits_array[hl_rate_control_histogram_ptr_temp->temporal_layer_index];
        EbBitNumber             *intra_sad_bits_array_ptr = rate_control_tables_ptr->intra_sad_bits_array[0];
        uint32_t                pred_bits_ref_qp = 0;

        if (hl_rate_control_histogram_ptr_temp->slice_type == I_SLICE) {
            // Loop over block in the frame and calculated the predicted bits at reg QP
            unsigned i;
            uint32_t accum = 0;
            for (i = 0; i < NUMBER_OF_INTRA_SAD_INTERVALS; ++i)
                accum += (uint32_t)((uint32_t)hl_rate_control_histogram_ptr_temp->ois_distortion_histogram[i] * (uint32_t)intra_sad_bits_array_ptr[i]);
            pred_bits_ref_qp = accum;
            total_bits += pred_bits_ref_qp;
        }
        else {
            unsigned i;
            uint32_t accum = 0;
            uint32_t accum_intra = 0;
            for (i = 0; i < NUMBER_OF_SAD_INTERVALS; ++i)
            {
                accum += (uint32_t)((uint32_t)hl_rate_control_histogram_ptr_temp->me_distortion_histogram[i] * (uint32_t)sad_bits_array_ptr[i]);
                accum_intra += (uint32_t)((uint32_t)hl_rate_control_histogram_ptr_temp->ois_distortion_histogram[i] * (uint32_t)intra_sad_bits_array_ptr[i]);
            }
            if (accum > accum_intra * 3)
                pred_bits_ref_qp = accum_intra;
            else
                pred_bits_ref_qp = accum;
            total_bits += pred_bits_ref_qp;
        }

        // Scale for in complete LCSs
        //  total_bits is normalized based on the area because of the sbs at the picture boundries
        total_bits = total_bits * (uint64_t)area_in_pixel / (hl_rate_control_histogram_ptr_temp->full_sb_count << 12);
    }
    return total_bits;
}

void high_level_rc_input_picture_vbr(
    PictureParentControlSet     *picture_control_set_ptr,
    SequenceControlSet          *sequence_control_set_ptr,
    EncodeContext               *encode_context_ptr,
    RateControlContext          *context_ptr,
    HighLevelRateControlContext *high_level_rate_control_ptr)
{
    EbBool                      end_of_sequence_flag = EB_TRUE;

    HlRateControlHistogramEntry *hl_rate_control_histogram_ptr_temp;
    // Queue variables
    uint32_t                     queue_entry_index_temp;
    uint32_t                     queue_entry_index_temp2;
    uint32_t                     queue_entry_index_head_temp;

    uint64_t                     min_la_bit_distance;
    uint32_t                     selected_ref_qp_table_index;
    uint32_t                     selected_ref_qp;
#if RC_UPDATE_TARGET_RATE
    uint32_t                     selected_org_ref_qp;
#endif
    uint32_t                     previous_selected_ref_qp = encode_context_ptr->previous_selected_ref_qp;
    uint64_t                     max_coded_poc = encode_context_ptr->max_coded_poc;
    uint32_t                     max_coded_poc_selected_ref_qp = encode_context_ptr->max_coded_poc_selected_ref_qp;

    uint32_t                     ref_qp_index;
    uint32_t                     ref_qp_index_temp;
    uint32_t                     ref_qp_table_index;

    uint32_t                     area_in_pixel;
    uint32_t                     num_of_full_sbs;
    uint32_t                     qp_search_min;
    uint32_t                     qp_search_max;
    int32_t                      qp_step = 1;
    EbBool                      best_qp_found;
    uint32_t                     temporal_layer_index;
    EbBool                      tables_updated;

    uint64_t                     bit_constraint_per_sw = 0;

    RateControlTables           *rate_control_tables_ptr;
    EbBitNumber                 *sad_bits_array_ptr;
    EbBitNumber                 *intra_sad_bits_array_ptr;
    uint32_t                     pred_bits_ref_qp;

    for (temporal_layer_index = 0; temporal_layer_index < EB_MAX_TEMPORAL_LAYERS; temporal_layer_index++)
        picture_control_set_ptr->bits_per_sw_per_layer[temporal_layer_index] = 0;
    picture_control_set_ptr->total_bits_per_gop = 0;

    area_in_pixel = sequence_control_set_ptr->seq_header.max_frame_width * sequence_control_set_ptr->seq_header.max_frame_height;;

    eb_block_on_mutex(sequence_control_set_ptr->encode_context_ptr->rate_table_update_mutex);

    tables_updated = sequence_control_set_ptr->encode_context_ptr->rate_control_tables_array_updated;
    picture_control_set_ptr->percentage_updated = EB_FALSE;
    if (sequence_control_set_ptr->static_config.look_ahead_distance != 0) {
        // Increamenting the head of the hl_rate_control_historgram_queue and clean up the entores
        hl_rate_control_histogram_ptr_temp = (encode_context_ptr->hl_rate_control_historgram_queue[encode_context_ptr->hl_rate_control_historgram_queue_head_index]);
        while ((hl_rate_control_histogram_ptr_temp->life_count == 0) && hl_rate_control_histogram_ptr_temp->passed_to_hlrc) {
            eb_block_on_mutex(sequence_control_set_ptr->encode_context_ptr->hl_rate_control_historgram_queue_mutex);
            // Reset the Reorder Queue Entry
            hl_rate_control_histogram_ptr_temp->picture_number += INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH;
            hl_rate_control_histogram_ptr_temp->life_count = -1;
            hl_rate_control_histogram_ptr_temp->passed_to_hlrc = EB_FALSE;
            hl_rate_control_histogram_ptr_temp->is_coded = EB_FALSE;
            hl_rate_control_histogram_ptr_temp->total_num_bits_coded = 0;

            // Increment the Reorder Queue head Ptr
            encode_context_ptr->hl_rate_control_historgram_queue_head_index =
                (encode_context_ptr->hl_rate_control_historgram_queue_head_index == HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1) ? 0 : encode_context_ptr->hl_rate_control_historgram_queue_head_index + 1;
            eb_release_mutex(sequence_control_set_ptr->encode_context_ptr->hl_rate_control_historgram_queue_mutex);
            hl_rate_control_histogram_ptr_temp = encode_context_ptr->hl_rate_control_historgram_queue[encode_context_ptr->hl_rate_control_historgram_queue_head_index];
        }
        // For the case that number of frames in the sliding window is less than size of the look ahead or intra Refresh. i.e. end of sequence
        if ((picture_control_set_ptr->frames_in_sw < MIN(sequence_control_set_ptr->static_config.look_ahead_distance + 1, (uint32_t)sequence_control_set_ptr->intra_period_length + 1))) {
            selected_ref_qp = max_coded_poc_selected_ref_qp;

            // Update the QP for the sliding window based on the status of RC
            if ((context_ptr->extra_bits_gen > (int64_t)(context_ptr->virtual_buffer_size << 3)))
                selected_ref_qp = (uint32_t)MAX((int32_t)selected_ref_qp - 2, 0);
            else if ((context_ptr->extra_bits_gen > (int64_t)(context_ptr->virtual_buffer_size << 2)))
                selected_ref_qp = (uint32_t)MAX((int32_t)selected_ref_qp - 1, 0);
            if ((context_ptr->extra_bits_gen < -(int64_t)(context_ptr->virtual_buffer_size << 2)))
                selected_ref_qp += 2;
            else if ((context_ptr->extra_bits_gen < -(int64_t)(context_ptr->virtual_buffer_size << 1)))
                selected_ref_qp += 1;
            if ((picture_control_set_ptr->frames_in_sw < (uint32_t)(sequence_control_set_ptr->intra_period_length + 1)) &&
                (picture_control_set_ptr->picture_number % ((sequence_control_set_ptr->intra_period_length + 1)) == 0)) {
                selected_ref_qp++;
            }

            selected_ref_qp = (uint32_t)CLIP3(
                sequence_control_set_ptr->static_config.min_qp_allowed,
                sequence_control_set_ptr->static_config.max_qp_allowed,
                selected_ref_qp);

            queue_entry_index_head_temp = (int32_t)(picture_control_set_ptr->picture_number - encode_context_ptr->hl_rate_control_historgram_queue[encode_context_ptr->hl_rate_control_historgram_queue_head_index]->picture_number);
            queue_entry_index_head_temp += encode_context_ptr->hl_rate_control_historgram_queue_head_index;
            queue_entry_index_head_temp = (queue_entry_index_head_temp > HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1) ?
                queue_entry_index_head_temp - HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH :
                queue_entry_index_head_temp;

            queue_entry_index_temp = queue_entry_index_head_temp;
            {
                hl_rate_control_histogram_ptr_temp = (encode_context_ptr->hl_rate_control_historgram_queue[queue_entry_index_temp]);

                if (hl_rate_control_histogram_ptr_temp->slice_type == I_SLICE)
                    ref_qp_index_temp = context_ptr->qp_scaling_map_I_SLICE[selected_ref_qp];
                else
                    ref_qp_index_temp = context_ptr->qp_scaling_map[hl_rate_control_histogram_ptr_temp->temporal_layer_index][selected_ref_qp];

                ref_qp_index_temp = (uint32_t)CLIP3(
                    sequence_control_set_ptr->static_config.min_qp_allowed,
                    sequence_control_set_ptr->static_config.max_qp_allowed,
                    ref_qp_index_temp);

                hl_rate_control_histogram_ptr_temp->pred_bits_ref_qp[ref_qp_index_temp] = 0;
                rate_control_tables_ptr = &encode_context_ptr->rate_control_tables_array[ref_qp_index_temp];
                sad_bits_array_ptr = rate_control_tables_ptr->sad_bits_array[hl_rate_control_histogram_ptr_temp->temporal_layer_index];
                intra_sad_bits_array_ptr = rate_control_tables_ptr->intra_sad_bits_array[hl_rate_control_histogram_ptr_temp->temporal_layer_index];
                pred_bits_ref_qp = 0;
                num_of_full_sbs = 0;

                if (hl_rate_control_histogram_ptr_temp->slice_type == I_SLICE) {
                    // Loop over block in the frame and calculated the predicted bits at reg QP
                    {
                        unsigned i;
                        uint32_t accum = 0;
                        for (i = 0; i < NUMBER_OF_INTRA_SAD_INTERVALS; ++i)
                            accum += (uint32_t)(hl_rate_control_histogram_ptr_temp->ois_distortion_histogram[i] * intra_sad_bits_array_ptr[i]);
                        pred_bits_ref_qp = accum;
                        num_of_full_sbs = hl_rate_control_histogram_ptr_temp->full_sb_count;
                    }
                    hl_rate_control_histogram_ptr_temp->pred_bits_ref_qp[ref_qp_index_temp] += pred_bits_ref_qp;
                }

                else {
                    {
                        unsigned i;
                        uint32_t accum = 0;
                        for (i = 0; i < NUMBER_OF_SAD_INTERVALS; ++i)
                            accum += (uint32_t)(hl_rate_control_histogram_ptr_temp->me_distortion_histogram[i] * sad_bits_array_ptr[i]);
                        pred_bits_ref_qp = accum;
                        num_of_full_sbs = hl_rate_control_histogram_ptr_temp->full_sb_count;
                    }
                    hl_rate_control_histogram_ptr_temp->pred_bits_ref_qp[ref_qp_index_temp] += pred_bits_ref_qp;
                }

                // Scale for in complete
                //  pred_bits_ref_qp is normalized based on the area because of the LCUs at the picture boundries
                hl_rate_control_histogram_ptr_temp->pred_bits_ref_qp[ref_qp_index_temp] = hl_rate_control_histogram_ptr_temp->pred_bits_ref_qp[ref_qp_index_temp] * (uint64_t)area_in_pixel / (num_of_full_sbs << 12);

                // Store the pred_bits_ref_qp for the first frame in the window to PCS
                picture_control_set_ptr->pred_bits_ref_qp[ref_qp_index_temp] = hl_rate_control_histogram_ptr_temp->pred_bits_ref_qp[ref_qp_index_temp];
            }
        }
        else {
            // Loop over the QPs and find the best QP
            min_la_bit_distance = MAX_UNSIGNED_VALUE;
            qp_search_min = (uint8_t)CLIP3(
                sequence_control_set_ptr->static_config.min_qp_allowed,
                MAX_REF_QP_NUM,//sequence_control_set_ptr->static_config.max_qp_allowed,
                (uint32_t)MAX((int32_t)sequence_control_set_ptr->qp - 40, 0));

            qp_search_max = (uint8_t)CLIP3(
                sequence_control_set_ptr->static_config.min_qp_allowed,
                MAX_REF_QP_NUM,//sequence_control_set_ptr->static_config.max_qp_allowed,
                sequence_control_set_ptr->qp + 40);

            for (ref_qp_table_index = qp_search_min; ref_qp_table_index < qp_search_max; ref_qp_table_index++)
                high_level_rate_control_ptr->pred_bits_ref_qpPerSw[ref_qp_table_index] = 0;
            bit_constraint_per_sw = high_level_rate_control_ptr->bit_constraint_per_sw * picture_control_set_ptr->frames_in_sw / (sequence_control_set_ptr->static_config.look_ahead_distance + 1);

            // Update the target rate for the sliding window based on the status of RC
            if ((context_ptr->extra_bits_gen > (int64_t)(context_ptr->virtual_buffer_size * 10)))
                bit_constraint_per_sw = bit_constraint_per_sw * 130 / 100;
            else if ((context_ptr->extra_bits_gen > (int64_t)(context_ptr->virtual_buffer_size << 3)))
                bit_constraint_per_sw = bit_constraint_per_sw * 120 / 100;
            else if ((context_ptr->extra_bits_gen > (int64_t)(context_ptr->virtual_buffer_size << 2)))
                bit_constraint_per_sw = bit_constraint_per_sw * 110 / 100;
            if ((context_ptr->extra_bits_gen < -(int64_t)(context_ptr->virtual_buffer_size << 3)))
                bit_constraint_per_sw = bit_constraint_per_sw * 80 / 100;
            else if ((context_ptr->extra_bits_gen < -(int64_t)(context_ptr->virtual_buffer_size << 2)))
                bit_constraint_per_sw = bit_constraint_per_sw * 90 / 100;
            // Loop over proper QPs and find the Predicted bits for that QP. Find the QP with the closest total predicted rate to target bits for the sliding window.
            previous_selected_ref_qp = CLIP3(
                qp_search_min,
                qp_search_max,
                previous_selected_ref_qp);
            ref_qp_table_index = previous_selected_ref_qp;
            selected_ref_qp_table_index = ref_qp_table_index;
            selected_ref_qp = selected_ref_qp_table_index;
            best_qp_found = EB_FALSE;
            while (ref_qp_table_index >= qp_search_min && ref_qp_table_index <= qp_search_max && !best_qp_found) {
                ref_qp_index = CLIP3(
                    sequence_control_set_ptr->static_config.min_qp_allowed,
                    MAX_REF_QP_NUM,//sequence_control_set_ptr->static_config.max_qp_allowed,
                    ref_qp_table_index);
                high_level_rate_control_ptr->pred_bits_ref_qpPerSw[ref_qp_index] = 0;

                // Finding the predicted bits for each frame in the sliding window at the reference Qp(s)
                queue_entry_index_head_temp = (int32_t)(picture_control_set_ptr->picture_number - encode_context_ptr->hl_rate_control_historgram_queue[encode_context_ptr->hl_rate_control_historgram_queue_head_index]->picture_number);
                queue_entry_index_head_temp += encode_context_ptr->hl_rate_control_historgram_queue_head_index;
                queue_entry_index_head_temp = (queue_entry_index_head_temp > HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1) ?
                    queue_entry_index_head_temp - HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH :
                    queue_entry_index_head_temp;

                queue_entry_index_temp = queue_entry_index_head_temp;
                // This is set to false, so the last frame would go inside the loop
                end_of_sequence_flag = EB_FALSE;

                while (!end_of_sequence_flag &&
                    queue_entry_index_temp <= queue_entry_index_head_temp + sequence_control_set_ptr->static_config.look_ahead_distance) {
                    queue_entry_index_temp2 = (queue_entry_index_temp > HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1) ? queue_entry_index_temp - HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH : queue_entry_index_temp;
                    hl_rate_control_histogram_ptr_temp = (encode_context_ptr->hl_rate_control_historgram_queue[queue_entry_index_temp2]);

                    if (hl_rate_control_histogram_ptr_temp->slice_type == I_SLICE)
                        ref_qp_index_temp = context_ptr->qp_scaling_map_I_SLICE[ref_qp_index];
                    else
                        ref_qp_index_temp = context_ptr->qp_scaling_map[hl_rate_control_histogram_ptr_temp->temporal_layer_index][ref_qp_index];

                    ref_qp_index_temp = (uint32_t)CLIP3(
                        sequence_control_set_ptr->static_config.min_qp_allowed,
                        sequence_control_set_ptr->static_config.max_qp_allowed,
                        ref_qp_index_temp);

                    hl_rate_control_histogram_ptr_temp->pred_bits_ref_qp[ref_qp_index_temp] = 0;

                    if (ref_qp_table_index == previous_selected_ref_qp) {
                        eb_block_on_mutex(sequence_control_set_ptr->encode_context_ptr->hl_rate_control_historgram_queue_mutex);
                        hl_rate_control_histogram_ptr_temp->life_count--;
                        eb_release_mutex(sequence_control_set_ptr->encode_context_ptr->hl_rate_control_historgram_queue_mutex);
                    }
                    hl_rate_control_histogram_ptr_temp->pred_bits_ref_qp[ref_qp_index_temp] = predict_bits(
                        encode_context_ptr,
                        hl_rate_control_histogram_ptr_temp,
                        ref_qp_index_temp,
                        area_in_pixel);

                    high_level_rate_control_ptr->pred_bits_ref_qpPerSw[ref_qp_index] += hl_rate_control_histogram_ptr_temp->pred_bits_ref_qp[ref_qp_index_temp];

                    // Store the pred_bits_ref_qp for the first frame in the window to PCS
                    if (queue_entry_index_head_temp == queue_entry_index_temp2)
                        picture_control_set_ptr->pred_bits_ref_qp[ref_qp_index_temp] = hl_rate_control_histogram_ptr_temp->pred_bits_ref_qp[ref_qp_index_temp];

                    end_of_sequence_flag = hl_rate_control_histogram_ptr_temp->end_of_sequence_flag;
                    queue_entry_index_temp++;
                }

                if (min_la_bit_distance >= (uint64_t)ABS((int64_t)high_level_rate_control_ptr->pred_bits_ref_qpPerSw[ref_qp_index] - (int64_t)bit_constraint_per_sw)) {
                    min_la_bit_distance = (uint64_t)ABS((int64_t)high_level_rate_control_ptr->pred_bits_ref_qpPerSw[ref_qp_index] - (int64_t)bit_constraint_per_sw);
                    selected_ref_qp_table_index = ref_qp_table_index;
                    selected_ref_qp = ref_qp_index;
                }
                else
                    best_qp_found = EB_TRUE;
                if (ref_qp_table_index == previous_selected_ref_qp) {
                    if (high_level_rate_control_ptr->pred_bits_ref_qpPerSw[ref_qp_index] > bit_constraint_per_sw)
                        qp_step = +1;
                    else
                        qp_step = -1;
                }
                ref_qp_table_index = (uint32_t)(ref_qp_table_index + qp_step);
            }
        }

#if RC_UPDATE_TARGET_RATE
        selected_org_ref_qp = selected_ref_qp;
        if (sequence_control_set_ptr->intra_period_length != -1 && picture_control_set_ptr->picture_number % ((sequence_control_set_ptr->intra_period_length + 1)) == 0 &&
            (int32_t)picture_control_set_ptr->frames_in_sw > sequence_control_set_ptr->intra_period_length) {
            if (picture_control_set_ptr->picture_number > 0)
                picture_control_set_ptr->intra_selected_org_qp = (uint8_t)selected_ref_qp;
            ref_qp_index = selected_ref_qp;
            high_level_rate_control_ptr->pred_bits_ref_qpPerSw[ref_qp_index] = 0;

            if (high_level_rate_control_ptr->pred_bits_ref_qpPerSw[ref_qp_index] == 0) {
                // Finding the predicted bits for each frame in the sliding window at the reference Qp(s)
                //queue_entry_index_temp = encode_context_ptr->hl_rate_control_historgram_queue_head_index;
                queue_entry_index_head_temp = (int32_t)(picture_control_set_ptr->picture_number - encode_context_ptr->hl_rate_control_historgram_queue[encode_context_ptr->hl_rate_control_historgram_queue_head_index]->picture_number);
                queue_entry_index_head_temp += encode_context_ptr->hl_rate_control_historgram_queue_head_index;
                queue_entry_index_head_temp = (queue_entry_index_head_temp > HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1) ?
                    queue_entry_index_head_temp - HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH :
                    queue_entry_index_head_temp;

                queue_entry_index_temp = queue_entry_index_head_temp;

                // This is set to false, so the last frame would go inside the loop
                end_of_sequence_flag = EB_FALSE;

                while (!end_of_sequence_flag &&
                    //queue_entry_index_temp <= encode_context_ptr->hl_rate_control_historgram_queue_head_index+sequence_control_set_ptr->static_config.look_ahead_distance){
                    queue_entry_index_temp <= queue_entry_index_head_temp + sequence_control_set_ptr->static_config.look_ahead_distance) {
                    queue_entry_index_temp2 = (queue_entry_index_temp > HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1) ? queue_entry_index_temp - HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH : queue_entry_index_temp;
                    hl_rate_control_histogram_ptr_temp = (encode_context_ptr->hl_rate_control_historgram_queue[queue_entry_index_temp2]);

                    if (hl_rate_control_histogram_ptr_temp->slice_type == I_SLICE)
                        ref_qp_index_temp = context_ptr->qp_scaling_map_I_SLICE[ref_qp_index];
                    else
                        ref_qp_index_temp = context_ptr->qp_scaling_map[hl_rate_control_histogram_ptr_temp->temporal_layer_index][ref_qp_index];

                    ref_qp_index_temp = (uint32_t)CLIP3(
                        sequence_control_set_ptr->static_config.min_qp_allowed,
                        sequence_control_set_ptr->static_config.max_qp_allowed,
                        ref_qp_index_temp);

                    hl_rate_control_histogram_ptr_temp->pred_bits_ref_qp[ref_qp_index_temp] = predict_bits(
                        encode_context_ptr,
                        hl_rate_control_histogram_ptr_temp,
                        ref_qp_index_temp,
                        area_in_pixel);

                    high_level_rate_control_ptr->pred_bits_ref_qpPerSw[ref_qp_index] += hl_rate_control_histogram_ptr_temp->pred_bits_ref_qp[ref_qp_index_temp];
                    // Store the pred_bits_ref_qp for the first frame in the window to PCS
                    //  if(encode_context_ptr->hl_rate_control_historgram_queue_head_index == queue_entry_index_temp2)
                    if (queue_entry_index_head_temp == queue_entry_index_temp2)
                        picture_control_set_ptr->pred_bits_ref_qp[ref_qp_index_temp] = hl_rate_control_histogram_ptr_temp->pred_bits_ref_qp[ref_qp_index_temp];

                    end_of_sequence_flag = hl_rate_control_histogram_ptr_temp->end_of_sequence_flag;
                    queue_entry_index_temp++;
                }
            }
        }
#endif
        picture_control_set_ptr->tables_updated = tables_updated;
        EbBool expensive_i_slice = EB_FALSE;
        // Looping over the window to find the percentage of bit allocation in each layer
        if ((sequence_control_set_ptr->intra_period_length != -1) &&
            ((int32_t)picture_control_set_ptr->frames_in_sw > sequence_control_set_ptr->intra_period_length) &&
            ((int32_t)picture_control_set_ptr->frames_in_sw > sequence_control_set_ptr->intra_period_length)) {
            uint64_t i_slice_bits = 0;

            if (picture_control_set_ptr->picture_number % ((sequence_control_set_ptr->intra_period_length + 1)) == 0) {
                queue_entry_index_head_temp = (int32_t)(picture_control_set_ptr->picture_number - encode_context_ptr->hl_rate_control_historgram_queue[encode_context_ptr->hl_rate_control_historgram_queue_head_index]->picture_number);
                queue_entry_index_head_temp += encode_context_ptr->hl_rate_control_historgram_queue_head_index;
                queue_entry_index_head_temp = (queue_entry_index_head_temp > HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1) ?
                    queue_entry_index_head_temp - HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH :
                    queue_entry_index_head_temp;

                queue_entry_index_temp = queue_entry_index_head_temp;

                // This is set to false, so the last frame would go inside the loop
                end_of_sequence_flag = EB_FALSE;

                while (!end_of_sequence_flag &&
                    queue_entry_index_temp <= queue_entry_index_head_temp + sequence_control_set_ptr->intra_period_length) {
                    queue_entry_index_temp2 = (queue_entry_index_temp > HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1) ? queue_entry_index_temp - HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH : queue_entry_index_temp;
                    hl_rate_control_histogram_ptr_temp = (encode_context_ptr->hl_rate_control_historgram_queue[queue_entry_index_temp2]);

                    if (hl_rate_control_histogram_ptr_temp->slice_type == I_SLICE)
                        ref_qp_index_temp = context_ptr->qp_scaling_map_I_SLICE[selected_ref_qp];
                    else
                        ref_qp_index_temp = context_ptr->qp_scaling_map[hl_rate_control_histogram_ptr_temp->temporal_layer_index][selected_ref_qp];

                    ref_qp_index_temp = (uint32_t)CLIP3(
                        sequence_control_set_ptr->static_config.min_qp_allowed,
                        sequence_control_set_ptr->static_config.max_qp_allowed,
                        ref_qp_index_temp);

                    if (queue_entry_index_temp == queue_entry_index_head_temp)
                        i_slice_bits = hl_rate_control_histogram_ptr_temp->pred_bits_ref_qp[ref_qp_index_temp];
                    picture_control_set_ptr->total_bits_per_gop += hl_rate_control_histogram_ptr_temp->pred_bits_ref_qp[ref_qp_index_temp];
                    picture_control_set_ptr->bits_per_sw_per_layer[hl_rate_control_histogram_ptr_temp->temporal_layer_index] += hl_rate_control_histogram_ptr_temp->pred_bits_ref_qp[ref_qp_index_temp];
                    picture_control_set_ptr->percentage_updated = EB_TRUE;

                    end_of_sequence_flag = hl_rate_control_histogram_ptr_temp->end_of_sequence_flag;
                    queue_entry_index_temp++;
                }
                if (i_slice_bits * 100 > 85 * picture_control_set_ptr->total_bits_per_gop)
                    expensive_i_slice = EB_TRUE;
                if (picture_control_set_ptr->total_bits_per_gop == 0) {
                    for (temporal_layer_index = 0; temporal_layer_index < EB_MAX_TEMPORAL_LAYERS; temporal_layer_index++)
                        picture_control_set_ptr->bits_per_sw_per_layer[temporal_layer_index] = rate_percentage_layer_array[sequence_control_set_ptr->static_config.hierarchical_levels][temporal_layer_index];
                }
            }
        }
        else {
            for (temporal_layer_index = 0; temporal_layer_index < EB_MAX_TEMPORAL_LAYERS; temporal_layer_index++)
                picture_control_set_ptr->bits_per_sw_per_layer[temporal_layer_index] = rate_percentage_layer_array[sequence_control_set_ptr->static_config.hierarchical_levels][temporal_layer_index];
        }
        if (expensive_i_slice) {
            if (tables_updated)
                selected_ref_qp = (uint32_t)MAX((int32_t)selected_ref_qp - 1, 0);
            else
                selected_ref_qp = (uint32_t)MAX((int32_t)selected_ref_qp - 3, 0);
            selected_ref_qp = (uint32_t)CLIP3(
                sequence_control_set_ptr->static_config.min_qp_allowed,
                sequence_control_set_ptr->static_config.max_qp_allowed,
                selected_ref_qp);
        }
        // Set the QP
        previous_selected_ref_qp = selected_ref_qp;
        if (picture_control_set_ptr->picture_number > max_coded_poc && picture_control_set_ptr->temporal_layer_index < 2 && !picture_control_set_ptr->end_of_sequence_region) {
            max_coded_poc = picture_control_set_ptr->picture_number;
            max_coded_poc_selected_ref_qp = previous_selected_ref_qp;
            encode_context_ptr->previous_selected_ref_qp = previous_selected_ref_qp;
            encode_context_ptr->max_coded_poc = max_coded_poc;
            encode_context_ptr->max_coded_poc_selected_ref_qp = max_coded_poc_selected_ref_qp;
        }

        if (picture_control_set_ptr->slice_type == I_SLICE)
            picture_control_set_ptr->best_pred_qp = (uint8_t)context_ptr->qp_scaling_map_I_SLICE[selected_ref_qp];
        else
            picture_control_set_ptr->best_pred_qp = (uint8_t)context_ptr->qp_scaling_map[picture_control_set_ptr->temporal_layer_index][selected_ref_qp];

        picture_control_set_ptr->best_pred_qp = (uint8_t)CLIP3(
            sequence_control_set_ptr->static_config.min_qp_allowed,
            sequence_control_set_ptr->static_config.max_qp_allowed,
            picture_control_set_ptr->best_pred_qp);

#if RC_UPDATE_TARGET_RATE
        if (picture_control_set_ptr->picture_number == 0) {
            high_level_rate_control_ptr->prev_intra_selected_ref_qp = selected_ref_qp;
            high_level_rate_control_ptr->prev_intra_org_selected_ref_qp = selected_ref_qp;
        }
        if (sequence_control_set_ptr->intra_period_length != -1) {
            if (picture_control_set_ptr->picture_number % ((sequence_control_set_ptr->intra_period_length + 1)) == 0) {
                high_level_rate_control_ptr->prev_intra_selected_ref_qp = selected_ref_qp;
                high_level_rate_control_ptr->prev_intra_org_selected_ref_qp = selected_org_ref_qp;
            }
        }
#endif
        picture_control_set_ptr->target_bits_best_pred_qp = picture_control_set_ptr->pred_bits_ref_qp[picture_control_set_ptr->best_pred_qp];
#if RC_PRINTS
        if (picture_control_set_ptr->slice_type == 2)
        {
            SVT_LOG("\nTID: %d\t", picture_control_set_ptr->temporal_layer_index);
            SVT_LOG("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t\n",
                picture_control_set_ptr->picture_number,
                picture_control_set_ptr->best_pred_qp,
                selected_ref_qp,
                (int)picture_control_set_ptr->target_bits_best_pred_qp,
                (int)high_level_rate_control_ptr->pred_bits_ref_qpPerSw[selected_ref_qp - 1],
                (int)high_level_rate_control_ptr->pred_bits_ref_qpPerSw[selected_ref_qp],
                (int)high_level_rate_control_ptr->pred_bits_ref_qpPerSw[selected_ref_qp + 1],
                (int)high_level_rate_control_ptr->bit_constraint_per_sw,
                (int)bit_constraint_per_sw/*,
                (int)high_level_rate_control_ptr->virtual_buffer_level*/);
        }
#endif
    }
    eb_release_mutex(sequence_control_set_ptr->encode_context_ptr->rate_table_update_mutex);
}
void frame_level_rc_input_picture_vbr(
    PictureControlSet               *picture_control_set_ptr,
    SequenceControlSet              *sequence_control_set_ptr,
    RateControlContext              *context_ptr,
    RateControlLayerContext         *rate_control_layer_ptr,
    RateControlIntervalParamContext *rate_control_param_ptr)
{
    RateControlLayerContext *rate_control_layer_temp_ptr;

    // Tiles
    uint32_t                 picture_area_in_pixel;
    uint32_t                 area_in_pixel;

    // SB Loop variables
    SbParams               *sb_params_ptr;
    uint32_t                 sb_index;
    uint64_t                 temp_qp;
    uint32_t                 area_in_sbs;

    picture_area_in_pixel = sequence_control_set_ptr->seq_header.max_frame_height*sequence_control_set_ptr->seq_header.max_frame_width;

    if (rate_control_layer_ptr->first_frame == 1) {
        rate_control_layer_ptr->first_frame = 0;
        picture_control_set_ptr->parent_pcs_ptr->first_frame_in_temporal_layer = 1;
    }
    else
        picture_control_set_ptr->parent_pcs_ptr->first_frame_in_temporal_layer = 0;
    if (picture_control_set_ptr->slice_type != I_SLICE) {
        if (rate_control_layer_ptr->first_non_intra_frame == 1) {
            rate_control_layer_ptr->first_non_intra_frame = 0;
            picture_control_set_ptr->parent_pcs_ptr->first_non_intra_frame_in_temporal_layer = 1;
        }
        else
            picture_control_set_ptr->parent_pcs_ptr->first_non_intra_frame_in_temporal_layer = 0;
    }
    else
        picture_control_set_ptr->parent_pcs_ptr->first_non_intra_frame_in_temporal_layer = 0;

    picture_control_set_ptr->parent_pcs_ptr->target_bits_rc = 0;

    // ***Rate Control***
    area_in_sbs = 0;
    area_in_pixel = 0;

    for (sb_index = 0; sb_index < picture_control_set_ptr->sb_total_count; ++sb_index) {
        sb_params_ptr = &sequence_control_set_ptr->sb_params_array[sb_index];

        if (sb_params_ptr->is_complete_sb) {
            // add the area of one LCU (64x64=4096) to the area of the tile
            area_in_pixel += 4096;
            area_in_sbs++;
        }
        else {
            // add the area of the LCU to the area of the tile
            area_in_pixel += sb_params_ptr->width * sb_params_ptr->height;
        }
    }
    rate_control_layer_ptr->area_in_pixel = area_in_pixel;

    if (picture_control_set_ptr->parent_pcs_ptr->first_frame_in_temporal_layer || (picture_control_set_ptr->picture_number == rate_control_param_ptr->first_poc)) {
        if (sequence_control_set_ptr->static_config.enable_qp_scaling_flag && (picture_control_set_ptr->picture_number != rate_control_param_ptr->first_poc)) {
            picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                (int32_t)sequence_control_set_ptr->static_config.min_qp_allowed,
                (int32_t)sequence_control_set_ptr->static_config.max_qp_allowed,
                (int32_t)(rate_control_param_ptr->intra_frames_qp + context_ptr->qp_scaling_map[picture_control_set_ptr->temporal_layer_index][rate_control_param_ptr->intra_frames_qp_bef_scal] - context_ptr->qp_scaling_map_I_SLICE[rate_control_param_ptr->intra_frames_qp_bef_scal]));
        }

        if (picture_control_set_ptr->picture_number == 0) {
            rate_control_param_ptr->intra_frames_qp = sequence_control_set_ptr->qp;
            rate_control_param_ptr->intra_frames_qp_bef_scal = (uint8_t)sequence_control_set_ptr->qp;
        }

        if (picture_control_set_ptr->picture_number == rate_control_param_ptr->first_poc) {
            uint32_t temporal_layer_idex;
            rate_control_param_ptr->previous_virtual_buffer_level = context_ptr->virtual_buffer_level_initial_value;
            rate_control_param_ptr->virtual_buffer_level = context_ptr->virtual_buffer_level_initial_value;
            rate_control_param_ptr->extra_ap_bit_ratio_i = 0;
            if (picture_control_set_ptr->parent_pcs_ptr->end_of_sequence_region) {
                rate_control_param_ptr->last_poc = MAX(rate_control_param_ptr->first_poc + picture_control_set_ptr->parent_pcs_ptr->frames_in_sw - 1, rate_control_param_ptr->first_poc);
                rate_control_param_ptr->last_gop = EB_TRUE;
            }

            if ((context_ptr->extra_bits > (int64_t)(context_ptr->virtual_buffer_size >> 8)) ||
                (context_ptr->extra_bits < -(int64_t)(context_ptr->virtual_buffer_size >> 8))) {
                int64_t extra_bits_per_gop = 0;

                if (picture_control_set_ptr->parent_pcs_ptr->end_of_sequence_region) {
                    if ((context_ptr->extra_bits > (int64_t)(context_ptr->virtual_buffer_size << 4)) ||
                        (context_ptr->extra_bits < -(int64_t)(context_ptr->virtual_buffer_size << 4))) {
                        extra_bits_per_gop = context_ptr->extra_bits;
                        extra_bits_per_gop = CLIP3(
                            -(int64_t)(context_ptr->vb_fill_threshold2 << 3),
                            (int64_t)(context_ptr->vb_fill_threshold2 << 3),
                            extra_bits_per_gop);
                    }
                    else
                        if ((context_ptr->extra_bits > (int64_t)(context_ptr->virtual_buffer_size << 3)) ||
                            (context_ptr->extra_bits < -(int64_t)(context_ptr->virtual_buffer_size << 3))) {
                            extra_bits_per_gop = context_ptr->extra_bits;
                            extra_bits_per_gop = CLIP3(
                                -(int64_t)(context_ptr->vb_fill_threshold2 << 2),
                                (int64_t)(context_ptr->vb_fill_threshold2 << 2),
                                extra_bits_per_gop);
                        }
                        else if ((context_ptr->extra_bits > (int64_t)(context_ptr->virtual_buffer_size << 2)) ||
                            (context_ptr->extra_bits < -(int64_t)(context_ptr->virtual_buffer_size << 2))) {
                            extra_bits_per_gop = CLIP3(
                                -(int64_t)context_ptr->vb_fill_threshold2 << 1,
                                (int64_t)context_ptr->vb_fill_threshold2 << 1,
                                extra_bits_per_gop);
                        }
                        else {
                            extra_bits_per_gop = CLIP3(
                                -(int64_t)context_ptr->vb_fill_threshold1,
                                (int64_t)context_ptr->vb_fill_threshold1,
                                extra_bits_per_gop);
                        }
                }
                else {
                    if ((context_ptr->extra_bits > (int64_t)(context_ptr->virtual_buffer_size << 3)) ||
                        (context_ptr->extra_bits < -(int64_t)(context_ptr->virtual_buffer_size << 3))) {
                        extra_bits_per_gop = context_ptr->extra_bits;
                        extra_bits_per_gop = CLIP3(
                            -(int64_t)(context_ptr->vb_fill_threshold2 << 2),
                            (int64_t)(context_ptr->vb_fill_threshold2 << 2),
                            extra_bits_per_gop);
                    }
                    else if ((context_ptr->extra_bits > (int64_t)(context_ptr->virtual_buffer_size << 2)) ||
                        (context_ptr->extra_bits < -(int64_t)(context_ptr->virtual_buffer_size << 2))) {
                        extra_bits_per_gop = CLIP3(
                            -(int64_t)context_ptr->vb_fill_threshold2 << 1,
                            (int64_t)context_ptr->vb_fill_threshold2 << 1,
                            extra_bits_per_gop);
                    }
                }

                rate_control_param_ptr->virtual_buffer_level -= extra_bits_per_gop;
                rate_control_param_ptr->previous_virtual_buffer_level -= extra_bits_per_gop;
                context_ptr->extra_bits -= extra_bits_per_gop;
            }

            for (temporal_layer_idex = 0; temporal_layer_idex < EB_MAX_TEMPORAL_LAYERS; temporal_layer_idex++) {
                rate_control_layer_temp_ptr = rate_control_param_ptr->rate_control_layer_array[temporal_layer_idex];
                rate_control_layer_reset(
                    rate_control_layer_temp_ptr,
                    picture_control_set_ptr,
                    context_ptr,
                    picture_area_in_pixel,
                    rate_control_param_ptr->was_used);
            }
        }

        picture_control_set_ptr->parent_pcs_ptr->sad_me = 0;
        // Finding the QP of the Intra frame by using variance tables
        if (picture_control_set_ptr->slice_type == I_SLICE) {
            uint32_t         selected_ref_qp;

            if (sequence_control_set_ptr->static_config.look_ahead_distance == 0)
                printf("ERROR: LAD=0 is not supported\n");
            else {
                selected_ref_qp = picture_control_set_ptr->parent_pcs_ptr->best_pred_qp;
                picture_control_set_ptr->picture_qp = (uint8_t)selected_ref_qp;
                picture_control_set_ptr->parent_pcs_ptr->calculated_qp = picture_control_set_ptr->picture_qp;
            }

            // Update the QP based on the VB
            if (picture_control_set_ptr->parent_pcs_ptr->end_of_sequence_region) {
                if (rate_control_param_ptr->virtual_buffer_level >= context_ptr->vb_fill_threshold2 << 1)
                    picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp + (uint8_t)THRESHOLD2QPINCREASE + 2;
                else if (rate_control_param_ptr->virtual_buffer_level >= context_ptr->vb_fill_threshold2)
                    picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp + (uint8_t)THRESHOLD2QPINCREASE;
                else if (rate_control_param_ptr->virtual_buffer_level >= context_ptr->vb_fill_threshold1 &&
                    rate_control_param_ptr->virtual_buffer_level < context_ptr->vb_fill_threshold2) {
                    picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp + (uint8_t)THRESHOLD1QPINCREASE;
                }
                if (rate_control_param_ptr->virtual_buffer_level <= -(context_ptr->vb_fill_threshold2 << 2))
                    picture_control_set_ptr->picture_qp = (uint8_t)MAX((int32_t)picture_control_set_ptr->picture_qp - (int32_t)THRESHOLD2QPINCREASE - (int32_t)2, 0);
                else
                    if (rate_control_param_ptr->virtual_buffer_level <= -(context_ptr->vb_fill_threshold2 << 1))
                        picture_control_set_ptr->picture_qp = (uint8_t)MAX((int32_t)picture_control_set_ptr->picture_qp - (int32_t)THRESHOLD2QPINCREASE - (int32_t)1, 0);
                    else if (rate_control_param_ptr->virtual_buffer_level <= 0)
                        picture_control_set_ptr->picture_qp = (uint8_t)MAX((int32_t)picture_control_set_ptr->picture_qp - (int32_t)THRESHOLD2QPINCREASE, 0);
            }
            else {
                if (rate_control_param_ptr->virtual_buffer_level >= context_ptr->vb_fill_threshold2)
                    picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp + (uint8_t)THRESHOLD2QPINCREASE;
                if (rate_control_param_ptr->virtual_buffer_level <= -(context_ptr->vb_fill_threshold2 << 2))
                    picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp - (uint8_t)THRESHOLD2QPINCREASE - (int32_t)2;
                else if (rate_control_param_ptr->virtual_buffer_level <= -(context_ptr->vb_fill_threshold2 << 1))
                    picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp - (uint8_t)THRESHOLD2QPINCREASE - (int32_t)1;
                else
                    if (rate_control_param_ptr->virtual_buffer_level <= 0)
                        picture_control_set_ptr->picture_qp = (uint8_t)MAX((int32_t)picture_control_set_ptr->picture_qp - (int32_t)THRESHOLD2QPINCREASE, 0);
            }
            picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                sequence_control_set_ptr->static_config.min_qp_allowed,
                sequence_control_set_ptr->static_config.max_qp_allowed,
                picture_control_set_ptr->picture_qp);
        }
        else {
            // LCU Loop
            for (sb_index = 0; sb_index < picture_control_set_ptr->sb_total_count; ++sb_index) {
                sb_params_ptr = &sequence_control_set_ptr->sb_params_array[sb_index];

                if (sb_params_ptr->is_complete_sb)
                    picture_control_set_ptr->parent_pcs_ptr->sad_me += picture_control_set_ptr->parent_pcs_ptr->rc_me_distortion[sb_index];
            }

            //  tileSadMe is normalized based on the area because of the LCUs at the tile boundries
            picture_control_set_ptr->parent_pcs_ptr->sad_me = MAX((picture_control_set_ptr->parent_pcs_ptr->sad_me*rate_control_layer_ptr->area_in_pixel / (area_in_sbs << 12)), 1);

            // totalSquareMad has RC_PRECISION precision
            picture_control_set_ptr->parent_pcs_ptr->sad_me <<= RC_PRECISION;
        }

        temp_qp = picture_control_set_ptr->picture_qp;

        if (picture_control_set_ptr->picture_number == rate_control_param_ptr->first_poc) {
            uint32_t temporal_layer_idex;
            for (temporal_layer_idex = 0; temporal_layer_idex < EB_MAX_TEMPORAL_LAYERS; temporal_layer_idex++) {
                rate_control_layer_temp_ptr = rate_control_param_ptr->rate_control_layer_array[temporal_layer_idex];
                rate_control_layer_reset_part2(
                    context_ptr,
                    rate_control_layer_temp_ptr,
                    picture_control_set_ptr);
            }
        }

        if (picture_control_set_ptr->picture_number == 0) {
            context_ptr->base_layer_frames_avg_qp = picture_control_set_ptr->picture_qp + 1;
            context_ptr->base_layer_intra_frames_avg_qp = picture_control_set_ptr->picture_qp;
        }
    }
    else {
        picture_control_set_ptr->parent_pcs_ptr->sad_me = 0;

        // if the pixture is an I slice, for now we set the QP as the QP of the previous frame
        if (picture_control_set_ptr->slice_type == I_SLICE) {
            uint32_t         selected_ref_qp;

            if (sequence_control_set_ptr->static_config.look_ahead_distance == 0)
                printf("ERROR: LAD=0 is not supported\n");
            else {
                selected_ref_qp = picture_control_set_ptr->parent_pcs_ptr->best_pred_qp;
                picture_control_set_ptr->picture_qp = (uint8_t)selected_ref_qp;
                picture_control_set_ptr->parent_pcs_ptr->calculated_qp = picture_control_set_ptr->picture_qp;
            }

            picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                sequence_control_set_ptr->static_config.min_qp_allowed,
                sequence_control_set_ptr->static_config.max_qp_allowed,
                picture_control_set_ptr->picture_qp);

            temp_qp = picture_control_set_ptr->picture_qp;
        }

        else { // Not an I slice
            // combining the target rate from initial RC and frame level RC
            if (sequence_control_set_ptr->static_config.look_ahead_distance != 0) {
                picture_control_set_ptr->parent_pcs_ptr->target_bits_rc = rate_control_layer_ptr->bit_constraint;
                rate_control_layer_ptr->ec_bit_constraint = (rate_control_layer_ptr->alpha * picture_control_set_ptr->parent_pcs_ptr->target_bits_best_pred_qp +
                    ((1 << RC_PRECISION) - rate_control_layer_ptr->alpha) * picture_control_set_ptr->parent_pcs_ptr->target_bits_rc + RC_PRECISION_OFFSET) >> RC_PRECISION;

                rate_control_layer_ptr->ec_bit_constraint = (uint64_t)MAX((int64_t)rate_control_layer_ptr->ec_bit_constraint - (int64_t)rate_control_layer_ptr->dif_total_and_ec_bits, 1);

                picture_control_set_ptr->parent_pcs_ptr->target_bits_rc = rate_control_layer_ptr->ec_bit_constraint;
            }

            // LCU Loop
            for (sb_index = 0; sb_index < picture_control_set_ptr->sb_total_count; ++sb_index) {
                sb_params_ptr = &sequence_control_set_ptr->sb_params_array[sb_index];

                if (sb_params_ptr->is_complete_sb)
                    picture_control_set_ptr->parent_pcs_ptr->sad_me += picture_control_set_ptr->parent_pcs_ptr->rc_me_distortion[sb_index];
            }

            //  tileSadMe is normalized based on the area because of the LCUs at the tile boundries
            picture_control_set_ptr->parent_pcs_ptr->sad_me = MAX((picture_control_set_ptr->parent_pcs_ptr->sad_me*rate_control_layer_ptr->area_in_pixel / (area_in_sbs << 12)), 1);
            picture_control_set_ptr->parent_pcs_ptr->sad_me <<= RC_PRECISION;
            if (rate_control_layer_ptr->area_in_pixel > 0)
                rate_control_layer_ptr->total_mad = MAX((picture_control_set_ptr->parent_pcs_ptr->sad_me / rate_control_layer_ptr->area_in_pixel), 1);
            if (!rate_control_layer_ptr->feedback_arrived)
                rate_control_layer_ptr->previous_frame_distortion_me = picture_control_set_ptr->parent_pcs_ptr->sad_me;
            {
                uint64_t qp_calc_temp1, qp_calc_temp2, qp_calc_temp3;

                qp_calc_temp1 = picture_control_set_ptr->parent_pcs_ptr->sad_me *rate_control_layer_ptr->total_mad;
                qp_calc_temp2 =
                    MAX((int64_t)(rate_control_layer_ptr->ec_bit_constraint << (2 * RC_PRECISION)) - (int64_t)rate_control_layer_ptr->c_coeff*(int64_t)rate_control_layer_ptr->area_in_pixel,
                    (int64_t)(rate_control_layer_ptr->ec_bit_constraint << (2 * RC_PRECISION - 2)));

                // This is a more complex but with higher precision implementation
                if (qp_calc_temp1 > qp_calc_temp2)
                    qp_calc_temp3 = (uint64_t)((qp_calc_temp1 / qp_calc_temp2)*rate_control_layer_ptr->k_coeff);
                else
                    qp_calc_temp3 = (uint64_t)(qp_calc_temp1*rate_control_layer_ptr->k_coeff / qp_calc_temp2);
                temp_qp = (uint64_t)(log2f_high_precision(MAX(((qp_calc_temp3 + RC_PRECISION_OFFSET) >> RC_PRECISION)*((qp_calc_temp3 + RC_PRECISION_OFFSET) >> RC_PRECISION)*((qp_calc_temp3 + RC_PRECISION_OFFSET) >> RC_PRECISION), 1), RC_PRECISION));

                rate_control_layer_ptr->calculated_frame_qp = (uint8_t)(CLIP3(1, 63, (uint32_t)(temp_qp + RC_PRECISION_OFFSET) >> RC_PRECISION));
                picture_control_set_ptr->parent_pcs_ptr->calculated_qp = (uint8_t)(CLIP3(1, 63, (uint32_t)(temp_qp + RC_PRECISION_OFFSET) >> RC_PRECISION));
            }

            temp_qp += rate_control_layer_ptr->delta_qp_fraction;
            picture_control_set_ptr->picture_qp = (uint8_t)((temp_qp + RC_PRECISION_OFFSET) >> RC_PRECISION);
            // Use the QP of HLRC instead of calculated one in FLRC
            if (picture_control_set_ptr->parent_pcs_ptr->hierarchical_levels > 1) {
                picture_control_set_ptr->picture_qp = picture_control_set_ptr->parent_pcs_ptr->best_pred_qp;
                picture_control_set_ptr->parent_pcs_ptr->calculated_qp = picture_control_set_ptr->parent_pcs_ptr->best_pred_qp;
            }
        }
        if (picture_control_set_ptr->parent_pcs_ptr->first_non_intra_frame_in_temporal_layer && picture_control_set_ptr->temporal_layer_index == 0 && picture_control_set_ptr->slice_type != I_SLICE)
            picture_control_set_ptr->picture_qp = (uint8_t)(rate_control_param_ptr->intra_frames_qp + context_ptr->qp_scaling_map[picture_control_set_ptr->temporal_layer_index][rate_control_param_ptr->intra_frames_qp_bef_scal] - context_ptr->qp_scaling_map_I_SLICE[rate_control_param_ptr->intra_frames_qp_bef_scal]);
        if (!rate_control_layer_ptr->feedback_arrived && picture_control_set_ptr->slice_type != I_SLICE) {
            picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                (int32_t)sequence_control_set_ptr->static_config.min_qp_allowed,
                (int32_t)sequence_control_set_ptr->static_config.max_qp_allowed,
                (int32_t)(rate_control_param_ptr->intra_frames_qp + context_ptr->qp_scaling_map[picture_control_set_ptr->temporal_layer_index][rate_control_param_ptr->intra_frames_qp_bef_scal] - context_ptr->qp_scaling_map_I_SLICE[rate_control_param_ptr->intra_frames_qp_bef_scal]));
        }

        if (picture_control_set_ptr->parent_pcs_ptr->end_of_sequence_region) {
            if (rate_control_param_ptr->virtual_buffer_level > context_ptr->vb_fill_threshold2 << 2)
                picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp + (uint8_t)THRESHOLD2QPINCREASE + 4;
            else if (rate_control_param_ptr->virtual_buffer_level > context_ptr->vb_fill_threshold2 << 1)
                picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp + (uint8_t)THRESHOLD2QPINCREASE + 3;
            else if (rate_control_param_ptr->virtual_buffer_level > context_ptr->vb_fill_threshold2)
                picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp + (uint8_t)THRESHOLD2QPINCREASE + 2;
            else if (rate_control_param_ptr->virtual_buffer_level > context_ptr->vb_fill_threshold1 &&
                rate_control_param_ptr->virtual_buffer_level < context_ptr->vb_fill_threshold2) {
                picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp + (uint8_t)THRESHOLD1QPINCREASE + 2;
            }
        }
        else {
            if (rate_control_param_ptr->virtual_buffer_level > context_ptr->vb_fill_threshold2 << 2)
                picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp + (uint8_t)THRESHOLD2QPINCREASE + 2;
            else if (rate_control_param_ptr->virtual_buffer_level > context_ptr->vb_fill_threshold2 << 1)
                picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp + (uint8_t)THRESHOLD2QPINCREASE + 1;
            else if (rate_control_param_ptr->virtual_buffer_level > context_ptr->vb_fill_threshold2)
                picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp + (uint8_t)THRESHOLD2QPINCREASE + 1;
            else if (rate_control_param_ptr->virtual_buffer_level > context_ptr->vb_fill_threshold1 &&
                rate_control_param_ptr->virtual_buffer_level < context_ptr->vb_fill_threshold2) {
                picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp + (uint8_t)THRESHOLD1QPINCREASE;
            }
        }
        if (picture_control_set_ptr->parent_pcs_ptr->end_of_sequence_region) {
            if (rate_control_param_ptr->virtual_buffer_level < -(context_ptr->vb_fill_threshold2 << 2))
                picture_control_set_ptr->picture_qp = (uint8_t)MAX((int32_t)picture_control_set_ptr->picture_qp - (int32_t)THRESHOLD2QPINCREASE - 2, 0);
            else if (rate_control_param_ptr->virtual_buffer_level < -(context_ptr->vb_fill_threshold2 << 1))
                picture_control_set_ptr->picture_qp = (uint8_t)MAX((int32_t)picture_control_set_ptr->picture_qp - (int32_t)THRESHOLD2QPINCREASE - 1, 0);
            else if (rate_control_param_ptr->virtual_buffer_level < 0)
                picture_control_set_ptr->picture_qp = (uint8_t)MAX((int32_t)picture_control_set_ptr->picture_qp - (int32_t)THRESHOLD2QPINCREASE, 0);
        }
        else {
            if (rate_control_param_ptr->virtual_buffer_level < -(context_ptr->vb_fill_threshold2 << 2))
                picture_control_set_ptr->picture_qp = (uint8_t)MAX((int32_t)picture_control_set_ptr->picture_qp - (int32_t)THRESHOLD2QPINCREASE - 1, 0);
            else if (rate_control_param_ptr->virtual_buffer_level < -context_ptr->vb_fill_threshold2)
                picture_control_set_ptr->picture_qp = (uint8_t)MAX((int32_t)picture_control_set_ptr->picture_qp - (int32_t)THRESHOLD2QPINCREASE, 0);
        }

        // limiting the QP based on the predicted QP
        if (sequence_control_set_ptr->static_config.look_ahead_distance != 0) {
            if (picture_control_set_ptr->parent_pcs_ptr->end_of_sequence_region) {
                picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                    (uint32_t)MAX((int32_t)picture_control_set_ptr->parent_pcs_ptr->best_pred_qp - 8, 0),
                    (uint32_t)picture_control_set_ptr->parent_pcs_ptr->best_pred_qp + 8,
                    (uint32_t)picture_control_set_ptr->picture_qp);
            }
            else {
                picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                    (uint32_t)MAX((int32_t)picture_control_set_ptr->parent_pcs_ptr->best_pred_qp - 8, 0),
                    (uint32_t)picture_control_set_ptr->parent_pcs_ptr->best_pred_qp + 8,
                    (uint32_t)picture_control_set_ptr->picture_qp);
            }
        }
        if (picture_control_set_ptr->picture_number != rate_control_param_ptr->first_poc &&
            picture_control_set_ptr->picture_qp == picture_control_set_ptr->parent_pcs_ptr->best_pred_qp && rate_control_param_ptr->virtual_buffer_level > context_ptr->vb_fill_threshold1) {
            if (rate_control_param_ptr->extra_ap_bit_ratio_i > 200)
                picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp + 3;
            else if (rate_control_param_ptr->extra_ap_bit_ratio_i > 100)
                picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp + 2;
            else if (rate_control_param_ptr->extra_ap_bit_ratio_i > 50)
                picture_control_set_ptr->picture_qp++;
        }
        //Limiting the QP based on the QP of the Reference frame

        uint32_t ref_qp;
        if ((int32_t)picture_control_set_ptr->temporal_layer_index == 0 && picture_control_set_ptr->slice_type != I_SLICE) {
            if (picture_control_set_ptr->ref_slice_type_array[0][0] == I_SLICE) {
                picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                    (uint32_t)picture_control_set_ptr->ref_pic_qp_array[0][0],
                    (uint32_t)picture_control_set_ptr->picture_qp,
                    picture_control_set_ptr->picture_qp);
            }
            else {
                picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                    (uint32_t)MAX((int32_t)picture_control_set_ptr->ref_pic_qp_array[0][0] - 1, 0),
                    (uint32_t)picture_control_set_ptr->picture_qp,
                    picture_control_set_ptr->picture_qp);
            }
        }
        else {
            ref_qp = 0;
            if (picture_control_set_ptr->ref_slice_type_array[0][0] != I_SLICE)
                ref_qp = MAX(ref_qp, picture_control_set_ptr->ref_pic_qp_array[0][0]);
            if ((picture_control_set_ptr->slice_type == B_SLICE) && (picture_control_set_ptr->ref_slice_type_array[1][0] != I_SLICE))
                ref_qp = MAX(ref_qp, picture_control_set_ptr->ref_pic_qp_array[1][0]);
            if (ref_qp > 0) {
                picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                    (uint32_t)ref_qp - 1,
                    picture_control_set_ptr->picture_qp,
                    picture_control_set_ptr->picture_qp);
            }
        }
        // limiting the QP between min Qp allowed and max Qp allowed
        picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
            sequence_control_set_ptr->static_config.min_qp_allowed,
            sequence_control_set_ptr->static_config.max_qp_allowed,
            picture_control_set_ptr->picture_qp);

        rate_control_layer_ptr->delta_qp_fraction = CLIP3(-RC_PRECISION_OFFSET, RC_PRECISION_OFFSET, -((int64_t)temp_qp - (int64_t)(picture_control_set_ptr->picture_qp << RC_PRECISION)));

        if (picture_control_set_ptr->parent_pcs_ptr->sad_me == rate_control_layer_ptr->previous_frame_distortion_me &&
            (rate_control_layer_ptr->previous_frame_distortion_me != 0))
            rate_control_layer_ptr->same_distortion_count++;
        else
            rate_control_layer_ptr->same_distortion_count = 0;
    }

    rate_control_layer_ptr->previous_c_coeff = rate_control_layer_ptr->c_coeff;
    rate_control_layer_ptr->previous_k_coeff = rate_control_layer_ptr->k_coeff;
    rate_control_layer_ptr->previous_calculated_frame_qp = rate_control_layer_ptr->calculated_frame_qp;
}

void frame_level_rc_feedback_picture_vbr(
    PictureParentControlSet *parentpicture_control_set_ptr,
    SequenceControlSet      *sequence_control_set_ptr,
    RateControlContext      *context_ptr)
{
    RateControlLayerContext         *rate_control_layer_temp_ptr;
    RateControlIntervalParamContext *rate_control_param_ptr;
    RateControlLayerContext         *rate_control_layer_ptr;
    // SB Loop variables
    uint32_t                         slice_num;
    uint64_t                         previous_frame_bit_actual;

    if (sequence_control_set_ptr->intra_period_length == -1)
        rate_control_param_ptr = context_ptr->rate_control_param_queue[0];
    else {
        uint32_t interval_index_temp = 0;
        while ((!(parentpicture_control_set_ptr->picture_number >= context_ptr->rate_control_param_queue[interval_index_temp]->first_poc &&
            parentpicture_control_set_ptr->picture_number <= context_ptr->rate_control_param_queue[interval_index_temp]->last_poc)) &&
            (interval_index_temp < PARALLEL_GOP_MAX_NUMBER)) {
            interval_index_temp++;
        }
        CHECK_REPORT_ERROR(
            interval_index_temp != PARALLEL_GOP_MAX_NUMBER,
            sequence_control_set_ptr->encode_context_ptr->app_callback_ptr,
            EB_ENC_RC_ERROR2);
        rate_control_param_ptr = context_ptr->rate_control_param_queue[interval_index_temp];
    }

    rate_control_layer_ptr = rate_control_param_ptr->rate_control_layer_array[parentpicture_control_set_ptr->temporal_layer_index];

    rate_control_layer_ptr->max_qp = 0;

    rate_control_layer_ptr->feedback_arrived = EB_TRUE;
    rate_control_layer_ptr->max_qp = MAX(rate_control_layer_ptr->max_qp, parentpicture_control_set_ptr->picture_qp);

    rate_control_layer_ptr->previous_frame_qp = parentpicture_control_set_ptr->picture_qp;
    rate_control_layer_ptr->previous_frame_bit_actual = parentpicture_control_set_ptr->total_num_bits;
    if (parentpicture_control_set_ptr->quantized_coeff_num_bits == 0)
        parentpicture_control_set_ptr->quantized_coeff_num_bits = 1;
    rate_control_layer_ptr->previous_framequantized_coeff_bit_actual = parentpicture_control_set_ptr->quantized_coeff_num_bits;

    // Setting Critical states for adjusting the averaging weights on C and K
    if ((parentpicture_control_set_ptr->sad_me > (3 * rate_control_layer_ptr->previous_frame_distortion_me) >> 1) &&
        (rate_control_layer_ptr->previous_frame_distortion_me != 0)) {
        rate_control_layer_ptr->critical_states = 3;
    }
    else if (rate_control_layer_ptr->critical_states)
        rate_control_layer_ptr->critical_states--;
    else
        rate_control_layer_ptr->critical_states = 0;
    if (parentpicture_control_set_ptr->slice_type != I_SLICE) {
        // Updating c_coeff
        rate_control_layer_ptr->c_coeff = (((int64_t)rate_control_layer_ptr->previous_frame_bit_actual - (int64_t)rate_control_layer_ptr->previous_framequantized_coeff_bit_actual) << (2 * RC_PRECISION))
            / rate_control_layer_ptr->area_in_pixel;
        rate_control_layer_ptr->c_coeff = MAX(rate_control_layer_ptr->c_coeff, 1);

        // Updating k_coeff
        if ((parentpicture_control_set_ptr->sad_me + RC_PRECISION_OFFSET) >> RC_PRECISION > 5) {
            {
                uint64_t test1, test2, test3;
                test1 = rate_control_layer_ptr->previous_framequantized_coeff_bit_actual*(two_to_power_qp_over_three[parentpicture_control_set_ptr->picture_qp]);
                test2 = MAX(parentpicture_control_set_ptr->sad_me / rate_control_layer_ptr->area_in_pixel, 1);
                test3 = test1 * 65536 / test2 * 65536 / parentpicture_control_set_ptr->sad_me;

                rate_control_layer_ptr->k_coeff = test3;
            }
        }

        if (rate_control_layer_ptr->critical_states) {
            rate_control_layer_ptr->k_coeff = (8 * rate_control_layer_ptr->k_coeff + 8 * rate_control_layer_ptr->previous_k_coeff + 8) >> 4;
            rate_control_layer_ptr->c_coeff = (8 * rate_control_layer_ptr->c_coeff + 8 * rate_control_layer_ptr->previous_c_coeff + 8) >> 4;
        }
        else {
            rate_control_layer_ptr->k_coeff = (rate_control_layer_ptr->coeff_averaging_weight1*rate_control_layer_ptr->k_coeff + rate_control_layer_ptr->coeff_averaging_weight2*rate_control_layer_ptr->previous_k_coeff + 8) >> 4;
            rate_control_layer_ptr->c_coeff = (rate_control_layer_ptr->coeff_averaging_weight1*rate_control_layer_ptr->c_coeff + rate_control_layer_ptr->coeff_averaging_weight2*rate_control_layer_ptr->previous_c_coeff + 8) >> 4;
        }
        rate_control_layer_ptr->k_coeff = MIN(rate_control_layer_ptr->k_coeff, rate_control_layer_ptr->previous_k_coeff * 4);
        rate_control_layer_ptr->c_coeff = MIN(rate_control_layer_ptr->c_coeff, rate_control_layer_ptr->previous_c_coeff * 4);
        if (parentpicture_control_set_ptr->slice_type != I_SLICE)
            rate_control_layer_ptr->previous_frame_distortion_me = parentpicture_control_set_ptr->sad_me;
        else
            rate_control_layer_ptr->previous_frame_distortion_me = 0;
    }

    if (sequence_control_set_ptr->static_config.look_ahead_distance != 0) {
        if (parentpicture_control_set_ptr->slice_type == I_SLICE) {
            if (parentpicture_control_set_ptr->total_num_bits < parentpicture_control_set_ptr->target_bits_best_pred_qp << 1)
                context_ptr->base_layer_intra_frames_avg_qp = (3 * context_ptr->base_layer_intra_frames_avg_qp + parentpicture_control_set_ptr->picture_qp + 2) >> 2;
            else if (parentpicture_control_set_ptr->total_num_bits > parentpicture_control_set_ptr->target_bits_best_pred_qp << 2)
                context_ptr->base_layer_intra_frames_avg_qp = (3 * context_ptr->base_layer_intra_frames_avg_qp + parentpicture_control_set_ptr->picture_qp + 4 + 2) >> 2;
            else if (parentpicture_control_set_ptr->total_num_bits > parentpicture_control_set_ptr->target_bits_best_pred_qp << 1)
                context_ptr->base_layer_intra_frames_avg_qp = (3 * context_ptr->base_layer_intra_frames_avg_qp + parentpicture_control_set_ptr->picture_qp + 2 + 2) >> 2;
        }
    }

    {
        uint64_t previous_frame_ec_bits = 0;
        EbBool picture_min_qp_allowed = EB_TRUE;
        rate_control_layer_ptr->previous_frame_average_qp = 0;
        rate_control_layer_ptr->previous_frame_average_qp += rate_control_layer_ptr->previous_frame_qp;
        previous_frame_ec_bits += rate_control_layer_ptr->previous_frame_bit_actual;
        if (rate_control_layer_ptr->same_distortion_count == 0 ||
            parentpicture_control_set_ptr->picture_qp != sequence_control_set_ptr->static_config.min_qp_allowed) {
            picture_min_qp_allowed = EB_FALSE;
        }
        if (picture_min_qp_allowed)
            rate_control_layer_ptr->frame_same_distortion_min_qp_count++;
        else
            rate_control_layer_ptr->frame_same_distortion_min_qp_count = 0;

        rate_control_layer_ptr->previous_ec_bits = previous_frame_ec_bits;
        previous_frame_bit_actual = parentpicture_control_set_ptr->total_num_bits;
        if (parentpicture_control_set_ptr->first_frame_in_temporal_layer)
            rate_control_layer_ptr->dif_total_and_ec_bits = (previous_frame_bit_actual - previous_frame_ec_bits);
        else
            rate_control_layer_ptr->dif_total_and_ec_bits = ((previous_frame_bit_actual - previous_frame_ec_bits) + rate_control_layer_ptr->dif_total_and_ec_bits) >> 1;
        // update bitrate of different layers in the interval based on the rate of the I frame
        if (parentpicture_control_set_ptr->picture_number == rate_control_param_ptr->first_poc &&
            (parentpicture_control_set_ptr->slice_type == I_SLICE) &&
            sequence_control_set_ptr->static_config.intra_period_length != -1) {
            uint32_t temporal_layer_idex;
            uint64_t target_bit_rate;
            uint64_t channel_bit_rate;
            uint64_t sum_bits_per_sw = 0;
#if ADAPTIVE_PERCENTAGE
            if (sequence_control_set_ptr->static_config.look_ahead_distance != 0) {
                if (parentpicture_control_set_ptr->tables_updated && parentpicture_control_set_ptr->percentage_updated) {
                    parentpicture_control_set_ptr->bits_per_sw_per_layer[0] =
                        (uint64_t)MAX((int64_t)parentpicture_control_set_ptr->bits_per_sw_per_layer[0] + (int64_t)parentpicture_control_set_ptr->total_num_bits - (int64_t)parentpicture_control_set_ptr->target_bits_best_pred_qp, 1);
                }
            }
#endif

            if (sequence_control_set_ptr->static_config.look_ahead_distance != 0 && sequence_control_set_ptr->intra_period_length != -1) {
                for (temporal_layer_idex = 0; temporal_layer_idex < EB_MAX_TEMPORAL_LAYERS; temporal_layer_idex++)
                    sum_bits_per_sw += parentpicture_control_set_ptr->bits_per_sw_per_layer[temporal_layer_idex];
            }

            for (temporal_layer_idex = 0; temporal_layer_idex < EB_MAX_TEMPORAL_LAYERS; temporal_layer_idex++) {
                rate_control_layer_temp_ptr = rate_control_param_ptr->rate_control_layer_array[temporal_layer_idex];

                target_bit_rate = (uint64_t)((int64_t)parentpicture_control_set_ptr->target_bit_rate -
                    MIN((int64_t)parentpicture_control_set_ptr->target_bit_rate * 3 / 4, (int64_t)(parentpicture_control_set_ptr->total_num_bits*context_ptr->frame_rate / (sequence_control_set_ptr->static_config.intra_period_length + 1)) >> RC_PRECISION))
                    *rate_percentage_layer_array[sequence_control_set_ptr->static_config.hierarchical_levels][temporal_layer_idex] / 100;

#if ADAPTIVE_PERCENTAGE
                if (sequence_control_set_ptr->static_config.look_ahead_distance != 0 && sequence_control_set_ptr->intra_period_length != -1) {
                    target_bit_rate = (uint64_t)((int64_t)parentpicture_control_set_ptr->target_bit_rate -
                        MIN((int64_t)parentpicture_control_set_ptr->target_bit_rate * 3 / 4, (int64_t)(parentpicture_control_set_ptr->total_num_bits*context_ptr->frame_rate / (sequence_control_set_ptr->static_config.intra_period_length + 1)) >> RC_PRECISION))
                        *parentpicture_control_set_ptr->bits_per_sw_per_layer[temporal_layer_idex] / sum_bits_per_sw;
                }
#endif
                // update this based on temporal layers
                if (temporal_layer_idex == 0)
                    channel_bit_rate = (((target_bit_rate << (2 * RC_PRECISION)) / MAX(1, rate_control_layer_temp_ptr->frame_rate - (1 * context_ptr->frame_rate / (sequence_control_set_ptr->static_config.intra_period_length + 1)))) + RC_PRECISION_OFFSET) >> RC_PRECISION;
                else
                    channel_bit_rate = (((target_bit_rate << (2 * RC_PRECISION)) / rate_control_layer_temp_ptr->frame_rate) + RC_PRECISION_OFFSET) >> RC_PRECISION;
                channel_bit_rate = (uint64_t)MAX((int64_t)1, (int64_t)channel_bit_rate);
                rate_control_layer_temp_ptr->ec_bit_constraint = channel_bit_rate;

                slice_num = 1;
                rate_control_layer_temp_ptr->ec_bit_constraint -= SLICE_HEADER_BITS_NUM * slice_num;

                rate_control_layer_temp_ptr->previous_bit_constraint = channel_bit_rate;
                rate_control_layer_temp_ptr->bit_constraint = channel_bit_rate;
                rate_control_layer_temp_ptr->channel_bit_rate = channel_bit_rate;
            }
            if ((int64_t)parentpicture_control_set_ptr->target_bit_rate * 3 / 4 < (int64_t)(parentpicture_control_set_ptr->total_num_bits*context_ptr->frame_rate / (sequence_control_set_ptr->static_config.intra_period_length + 1)) >> RC_PRECISION) {
                rate_control_param_ptr->previous_virtual_buffer_level += (int64_t)((parentpicture_control_set_ptr->total_num_bits*context_ptr->frame_rate / (sequence_control_set_ptr->static_config.intra_period_length + 1)) >> RC_PRECISION) - (int64_t)parentpicture_control_set_ptr->target_bit_rate * 3 / 4;
                context_ptr->extra_bits_gen -= (int64_t)((parentpicture_control_set_ptr->total_num_bits*context_ptr->frame_rate / (sequence_control_set_ptr->static_config.intra_period_length + 1)) >> RC_PRECISION) - (int64_t)parentpicture_control_set_ptr->target_bit_rate * 3 / 4;
            }
        }

        if (previous_frame_bit_actual) {
            uint64_t bit_changes_rate;
            // Updating virtual buffer level and it can be negative
            if ((parentpicture_control_set_ptr->picture_number == rate_control_param_ptr->first_poc) &&
                (parentpicture_control_set_ptr->slice_type == I_SLICE) &&
                (rate_control_param_ptr->last_gop == EB_FALSE) &&
                sequence_control_set_ptr->static_config.intra_period_length != -1) {
                rate_control_param_ptr->virtual_buffer_level =
                    (int64_t)rate_control_param_ptr->previous_virtual_buffer_level;
            }
            else {
                rate_control_param_ptr->virtual_buffer_level =
                    (int64_t)rate_control_param_ptr->previous_virtual_buffer_level +
                    (int64_t)previous_frame_bit_actual - (int64_t)rate_control_layer_ptr->channel_bit_rate;
                context_ptr->extra_bits_gen -= (int64_t)previous_frame_bit_actual - (int64_t)rate_control_layer_ptr->channel_bit_rate;
            }
            if (parentpicture_control_set_ptr->hierarchical_levels > 1 && rate_control_layer_ptr->frame_same_distortion_min_qp_count > 10) {
                rate_control_layer_ptr->previous_bit_constraint = (int64_t)rate_control_layer_ptr->channel_bit_rate;
                rate_control_param_ptr->virtual_buffer_level = ((int64_t)context_ptr->virtual_buffer_size >> 1);
            }
            // Updating bit difference
            rate_control_layer_ptr->bit_diff = (int64_t)rate_control_param_ptr->virtual_buffer_level
                //- ((int64_t)context_ptr->virtual_buffer_size>>1);
                - ((int64_t)rate_control_layer_ptr->channel_bit_rate >> 1);

            // Limit the bit difference
            rate_control_layer_ptr->bit_diff = CLIP3(-(int64_t)(rate_control_layer_ptr->channel_bit_rate), (int64_t)(rate_control_layer_ptr->channel_bit_rate >> 1), rate_control_layer_ptr->bit_diff);
            bit_changes_rate = rate_control_layer_ptr->frame_rate;

            // Updating bit Constraint
            rate_control_layer_ptr->bit_constraint = MAX((int64_t)rate_control_layer_ptr->previous_bit_constraint - ((rate_control_layer_ptr->bit_diff << RC_PRECISION) / ((int64_t)bit_changes_rate)), 1);

            // Limiting the bit_constraint
            if (parentpicture_control_set_ptr->temporal_layer_index == 0) {
                rate_control_layer_ptr->bit_constraint = CLIP3(rate_control_layer_ptr->channel_bit_rate >> 2,
                    rate_control_layer_ptr->channel_bit_rate * 200 / 100,
                    rate_control_layer_ptr->bit_constraint);
            }
            else {
                rate_control_layer_ptr->bit_constraint = CLIP3(rate_control_layer_ptr->channel_bit_rate >> 1,
                    rate_control_layer_ptr->channel_bit_rate * 200 / 100,
                    rate_control_layer_ptr->bit_constraint);
            }
            rate_control_layer_ptr->ec_bit_constraint = (uint64_t)MAX((int64_t)rate_control_layer_ptr->bit_constraint - (int64_t)rate_control_layer_ptr->dif_total_and_ec_bits, 1);
            rate_control_param_ptr->previous_virtual_buffer_level = rate_control_param_ptr->virtual_buffer_level;
            rate_control_layer_ptr->previous_bit_constraint = rate_control_layer_ptr->bit_constraint;
        }

        rate_control_param_ptr->processed_frames_number++;
        rate_control_param_ptr->in_use = EB_TRUE;
        // check if all the frames in the interval have arrived
        if (rate_control_param_ptr->processed_frames_number == (rate_control_param_ptr->last_poc - rate_control_param_ptr->first_poc + 1) &&
            sequence_control_set_ptr->intra_period_length != -1) {
            uint32_t temporal_index;
            int64_t extra_bits;
            rate_control_param_ptr->first_poc += PARALLEL_GOP_MAX_NUMBER * (uint32_t)(sequence_control_set_ptr->intra_period_length + 1);
            rate_control_param_ptr->last_poc += PARALLEL_GOP_MAX_NUMBER * (uint32_t)(sequence_control_set_ptr->intra_period_length + 1);
            rate_control_param_ptr->processed_frames_number = 0;
            rate_control_param_ptr->extra_ap_bit_ratio_i = 0;
            rate_control_param_ptr->in_use = EB_FALSE;
            rate_control_param_ptr->was_used = EB_TRUE;
            rate_control_param_ptr->last_gop = EB_FALSE;
            rate_control_param_ptr->first_pic_actual_qp_assigned = EB_FALSE;
            for (temporal_index = 0; temporal_index < EB_MAX_TEMPORAL_LAYERS; temporal_index++) {
                rate_control_param_ptr->rate_control_layer_array[temporal_index]->first_frame = 1;
                rate_control_param_ptr->rate_control_layer_array[temporal_index]->first_non_intra_frame = 1;
                rate_control_param_ptr->rate_control_layer_array[temporal_index]->feedback_arrived = EB_FALSE;
            }
            extra_bits = ((int64_t)context_ptr->virtual_buffer_size >> 1) - (int64_t)rate_control_param_ptr->virtual_buffer_level;

            rate_control_param_ptr->virtual_buffer_level = context_ptr->virtual_buffer_size >> 1;
            context_ptr->extra_bits += extra_bits;
        }
        // Allocate the extra_bits among other GOPs
        if ((parentpicture_control_set_ptr->temporal_layer_index <= 2) &&
            ((context_ptr->extra_bits > (int64_t)(context_ptr->virtual_buffer_size >> 8)) ||
            (context_ptr->extra_bits < -(int64_t)(context_ptr->virtual_buffer_size >> 8)))) {
            uint32_t interval_index_temp, interval_in_use_count;
            int64_t extra_bits_per_gop;
            int64_t extra_bits = context_ptr->extra_bits;
            int32_t clip_coef1, clip_coef2;
            if (parentpicture_control_set_ptr->end_of_sequence_region) {
                clip_coef1 = -1;
                clip_coef2 = -1;
            }
            else {
                if (context_ptr->extra_bits > (int64_t)(context_ptr->virtual_buffer_size << 3) ||
                    context_ptr->extra_bits < -(int64_t)(context_ptr->virtual_buffer_size << 3)) {
                    clip_coef1 = 0;
                    clip_coef2 = 0;
                }
                else {
                    clip_coef1 = 2;
                    clip_coef2 = 4;
                }
            }

            interval_in_use_count = 0;

            if (extra_bits > 0) {
                // Extra bits to be distributed
                // Distribute it among those that are consuming more
                for (interval_index_temp = 0; interval_index_temp < PARALLEL_GOP_MAX_NUMBER; interval_index_temp++) {
                    if (context_ptr->rate_control_param_queue[interval_index_temp]->in_use &&
                        context_ptr->rate_control_param_queue[interval_index_temp]->virtual_buffer_level > ((int64_t)context_ptr->virtual_buffer_size >> 1)) {
                        interval_in_use_count++;
                    }
                }
                // Distribute the rate among them
                if (interval_in_use_count) {
                    extra_bits_per_gop = extra_bits / interval_in_use_count;
                    if (clip_coef1 > 0)
                        extra_bits_per_gop = CLIP3(
                            -(int64_t)context_ptr->virtual_buffer_size >> clip_coef1,
                            (int64_t)context_ptr->virtual_buffer_size >> clip_coef1,
                            extra_bits_per_gop);
                    else
                        extra_bits_per_gop = CLIP3(
                            -(int64_t)context_ptr->virtual_buffer_size << (-clip_coef1),
                            (int64_t)context_ptr->virtual_buffer_size << (-clip_coef1),
                            extra_bits_per_gop);

                    for (interval_index_temp = 0; interval_index_temp < PARALLEL_GOP_MAX_NUMBER; interval_index_temp++) {
                        if (context_ptr->rate_control_param_queue[interval_index_temp]->in_use &&
                            context_ptr->rate_control_param_queue[interval_index_temp]->virtual_buffer_level > ((int64_t)context_ptr->virtual_buffer_size >> 1)) {
                            context_ptr->rate_control_param_queue[interval_index_temp]->virtual_buffer_level -= extra_bits_per_gop;
                            context_ptr->rate_control_param_queue[interval_index_temp]->previous_virtual_buffer_level -= extra_bits_per_gop;
                            context_ptr->extra_bits -= extra_bits_per_gop;
                        }
                    }
                }
                // if no interval with more consuming was found, allocate it to ones with consuming less
                else {
                    interval_in_use_count = 0;
                    // Distribute it among those that are consuming less
                    for (interval_index_temp = 0; interval_index_temp < PARALLEL_GOP_MAX_NUMBER; interval_index_temp++) {
                        if (context_ptr->rate_control_param_queue[interval_index_temp]->in_use &&
                            context_ptr->rate_control_param_queue[interval_index_temp]->virtual_buffer_level <= ((int64_t)context_ptr->virtual_buffer_size >> 1)) {
                            interval_in_use_count++;
                        }
                    }
                    if (interval_in_use_count) {
                        extra_bits_per_gop = extra_bits / interval_in_use_count;
                        if (clip_coef2 > 0)
                            extra_bits_per_gop = CLIP3(
                                -(int64_t)context_ptr->virtual_buffer_size >> clip_coef2,
                                (int64_t)context_ptr->virtual_buffer_size >> clip_coef2,
                                extra_bits_per_gop);
                        else
                            extra_bits_per_gop = CLIP3(
                                -(int64_t)context_ptr->virtual_buffer_size << (-clip_coef2),
                                (int64_t)context_ptr->virtual_buffer_size << (-clip_coef2),
                                extra_bits_per_gop);
                        // Distribute the rate among them
                        for (interval_index_temp = 0; interval_index_temp < PARALLEL_GOP_MAX_NUMBER; interval_index_temp++) {
                            if (context_ptr->rate_control_param_queue[interval_index_temp]->in_use &&
                                context_ptr->rate_control_param_queue[interval_index_temp]->virtual_buffer_level <= ((int64_t)context_ptr->virtual_buffer_size >> 1)) {
                                context_ptr->rate_control_param_queue[interval_index_temp]->virtual_buffer_level -= extra_bits_per_gop;
                                context_ptr->rate_control_param_queue[interval_index_temp]->previous_virtual_buffer_level -= extra_bits_per_gop;
                                context_ptr->extra_bits -= extra_bits_per_gop;
                            }
                        }
                    }
                }
            }
            else {
                // Distribute it among those that are consuming less
                for (interval_index_temp = 0; interval_index_temp < PARALLEL_GOP_MAX_NUMBER; interval_index_temp++) {
                    if (context_ptr->rate_control_param_queue[interval_index_temp]->in_use &&
                        context_ptr->rate_control_param_queue[interval_index_temp]->virtual_buffer_level < ((int64_t)context_ptr->virtual_buffer_size >> 1)) {
                        interval_in_use_count++;
                    }
                }
                if (interval_in_use_count) {
                    extra_bits_per_gop = extra_bits / interval_in_use_count;
                    if (clip_coef1 > 0)
                        extra_bits_per_gop = CLIP3(
                            -(int64_t)context_ptr->virtual_buffer_size >> clip_coef1,
                            (int64_t)context_ptr->virtual_buffer_size >> clip_coef1,
                            extra_bits_per_gop);
                    else
                        extra_bits_per_gop = CLIP3(
                            -(int64_t)context_ptr->virtual_buffer_size << (-clip_coef1),
                            (int64_t)context_ptr->virtual_buffer_size << (-clip_coef1),
                            extra_bits_per_gop);
                    // Distribute the rate among them
                    for (interval_index_temp = 0; interval_index_temp < PARALLEL_GOP_MAX_NUMBER; interval_index_temp++) {
                        if (context_ptr->rate_control_param_queue[interval_index_temp]->in_use &&
                            context_ptr->rate_control_param_queue[interval_index_temp]->virtual_buffer_level < ((int64_t)context_ptr->virtual_buffer_size >> 1)) {
                            context_ptr->rate_control_param_queue[interval_index_temp]->virtual_buffer_level -= extra_bits_per_gop;
                            context_ptr->rate_control_param_queue[interval_index_temp]->previous_virtual_buffer_level -= extra_bits_per_gop;
                            context_ptr->extra_bits -= extra_bits_per_gop;
                        }
                    }
                }
                // if no interval with less consuming was found, allocate it to ones with consuming more
                else {
                    interval_in_use_count = 0;
                    for (interval_index_temp = 0; interval_index_temp < PARALLEL_GOP_MAX_NUMBER; interval_index_temp++) {
                        if (context_ptr->rate_control_param_queue[interval_index_temp]->in_use &&
                            context_ptr->rate_control_param_queue[interval_index_temp]->virtual_buffer_level < (int64_t)(context_ptr->virtual_buffer_size)) {
                            interval_in_use_count++;
                        }
                    }
                    if (interval_in_use_count) {
                        extra_bits_per_gop = extra_bits / interval_in_use_count;
                        if (clip_coef2 > 0)
                            extra_bits_per_gop = CLIP3(
                                -(int64_t)context_ptr->virtual_buffer_size >> clip_coef2,
                                (int64_t)context_ptr->virtual_buffer_size >> clip_coef2,
                                extra_bits_per_gop);
                        else
                            extra_bits_per_gop = CLIP3(
                                -(int64_t)context_ptr->virtual_buffer_size << (-clip_coef2),
                                (int64_t)context_ptr->virtual_buffer_size << (-clip_coef2),
                                extra_bits_per_gop);
                        // Distribute the rate among them
                        for (interval_index_temp = 0; interval_index_temp < PARALLEL_GOP_MAX_NUMBER; interval_index_temp++) {
                            if (context_ptr->rate_control_param_queue[interval_index_temp]->in_use &&
                                context_ptr->rate_control_param_queue[interval_index_temp]->virtual_buffer_level < (int64_t)(context_ptr->virtual_buffer_size)) {
                                context_ptr->rate_control_param_queue[interval_index_temp]->virtual_buffer_level -= extra_bits_per_gop;
                                context_ptr->rate_control_param_queue[interval_index_temp]->previous_virtual_buffer_level -= extra_bits_per_gop;
                                context_ptr->extra_bits -= extra_bits_per_gop;
                            }
                        }
                    }
                }
            }
        }
    }

#if RC_PRINTS
    if (parentpicture_control_set_ptr->temporal_layer_index == 0)
    {
        SVT_LOG("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.0f\t%.0f\t%.0f\t%.0f\t%d\t%d\n",
            (int)parentpicture_control_set_ptr->slice_type,
            (int)parentpicture_control_set_ptr->picture_number,
            (int)parentpicture_control_set_ptr->temporal_layer_index,
            (int)parentpicture_control_set_ptr->picture_qp, (int)parentpicture_control_set_ptr->calculated_qp, (int)parentpicture_control_set_ptr->best_pred_qp,
            (int)previous_frame_bit_actual,
            (int)parentpicture_control_set_ptr->target_bits_best_pred_qp,
            (int)parentpicture_control_set_ptr->target_bits_rc,
            (int)rate_control_layer_ptr->channel_bit_rate,
            (int)rate_control_layer_ptr->bit_constraint,
            (double)rate_control_layer_ptr->c_coeff,
            (double)rate_control_layer_ptr->k_coeff,
            (double)parentpicture_control_set_ptr->sad_me,
            (double)context_ptr->extra_bits_gen,
            (int)rate_control_param_ptr->virtual_buffer_level,
            (int)context_ptr->extra_bits);
    }
#endif
}
void high_level_rc_input_picture_cvbr(
    PictureParentControlSet     *picture_control_set_ptr,
    SequenceControlSet          *sequence_control_set_ptr,
    EncodeContext               *encode_context_ptr,
    RateControlContext          *context_ptr,
    HighLevelRateControlContext *high_level_rate_control_ptr)
{
    EbBool                      end_of_sequence_flag = EB_TRUE;

    HlRateControlHistogramEntry *hl_rate_control_histogram_ptr_temp;
    // Queue variables
    uint32_t                     queue_entry_index_temp;
    uint32_t                     queue_entry_index_temp2;
    uint32_t                     queue_entry_index_head_temp;

    uint64_t                     min_la_bit_distance;
    uint32_t                     selected_ref_qp_table_index;
    uint32_t                     selected_ref_qp;
#if RC_UPDATE_TARGET_RATE
    uint32_t                     selected_org_ref_qp;
#endif
    uint32_t                     previous_selected_ref_qp = encode_context_ptr->previous_selected_ref_qp;
    uint64_t                     max_coded_poc = encode_context_ptr->max_coded_poc;
    uint32_t                     max_coded_poc_selected_ref_qp = encode_context_ptr->max_coded_poc_selected_ref_qp;

    uint32_t                     ref_qp_index;
    uint32_t                     ref_qp_index_temp;
    uint32_t                     ref_qp_table_index;

    uint32_t                     area_in_pixel;
    uint32_t                     num_of_full_sbs;
    uint32_t                     qp_search_min;
    uint32_t                     qp_search_max;
    int32_t                      qp_step = 1;
    EbBool                      best_qp_found;
    uint32_t                     temporal_layer_index;
    EbBool                      tables_updated;

    uint64_t                     bit_constraint_per_sw = 0;

    RateControlTables           *rate_control_tables_ptr;
    EbBitNumber                 *sad_bits_array_ptr;
    EbBitNumber                 *intra_sad_bits_array_ptr;
    uint32_t                     pred_bits_ref_qp;
    int delta_qp = 0;

    for (temporal_layer_index = 0; temporal_layer_index < EB_MAX_TEMPORAL_LAYERS; temporal_layer_index++)
        picture_control_set_ptr->bits_per_sw_per_layer[temporal_layer_index] = 0;
    picture_control_set_ptr->total_bits_per_gop = 0;

    area_in_pixel = sequence_control_set_ptr->seq_header.max_frame_width * sequence_control_set_ptr->seq_header.max_frame_height;;

    eb_block_on_mutex(sequence_control_set_ptr->encode_context_ptr->rate_table_update_mutex);

    tables_updated = sequence_control_set_ptr->encode_context_ptr->rate_control_tables_array_updated;
    picture_control_set_ptr->percentage_updated = EB_FALSE;
    if (sequence_control_set_ptr->static_config.look_ahead_distance != 0) {
        // Increamenting the head of the hl_rate_control_historgram_queue and clean up the entores
        hl_rate_control_histogram_ptr_temp = (encode_context_ptr->hl_rate_control_historgram_queue[encode_context_ptr->hl_rate_control_historgram_queue_head_index]);

        while ((hl_rate_control_histogram_ptr_temp->life_count == 0) && hl_rate_control_histogram_ptr_temp->passed_to_hlrc) {
            eb_block_on_mutex(sequence_control_set_ptr->encode_context_ptr->hl_rate_control_historgram_queue_mutex);
            // Reset the Reorder Queue Entry
            hl_rate_control_histogram_ptr_temp->picture_number += INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH;
            hl_rate_control_histogram_ptr_temp->life_count = -1;
            hl_rate_control_histogram_ptr_temp->passed_to_hlrc = EB_FALSE;
            hl_rate_control_histogram_ptr_temp->is_coded = EB_FALSE;
            hl_rate_control_histogram_ptr_temp->total_num_bits_coded = 0;

            // Increment the Reorder Queue head Ptr
            encode_context_ptr->hl_rate_control_historgram_queue_head_index =
                (encode_context_ptr->hl_rate_control_historgram_queue_head_index == HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1) ? 0 : encode_context_ptr->hl_rate_control_historgram_queue_head_index + 1;
            eb_release_mutex(sequence_control_set_ptr->encode_context_ptr->hl_rate_control_historgram_queue_mutex);
            hl_rate_control_histogram_ptr_temp = encode_context_ptr->hl_rate_control_historgram_queue[encode_context_ptr->hl_rate_control_historgram_queue_head_index];
        }
        // For the case that number of frames in the sliding window is less than size of the look ahead or intra Refresh. i.e. end of sequence
        if ((picture_control_set_ptr->frames_in_sw < MIN(sequence_control_set_ptr->static_config.look_ahead_distance + 1, (uint32_t)sequence_control_set_ptr->intra_period_length + 1))) {
            selected_ref_qp = max_coded_poc_selected_ref_qp;

            // Update the QP for the sliding window based on the status of RC
            if ((context_ptr->extra_bits_gen > (int64_t)(context_ptr->virtual_buffer_size << 3)))
                selected_ref_qp = (uint32_t)MAX((int32_t)selected_ref_qp - 2, 0);
            else if ((context_ptr->extra_bits_gen > (int64_t)(context_ptr->virtual_buffer_size << 2)))
                selected_ref_qp = (uint32_t)MAX((int32_t)selected_ref_qp - 1, 0);
            if ((context_ptr->extra_bits_gen < -(int64_t)(context_ptr->virtual_buffer_size << 2)))
                selected_ref_qp += 2;
            else if ((context_ptr->extra_bits_gen < -(int64_t)(context_ptr->virtual_buffer_size << 1)))
                selected_ref_qp += 1;
            if ((picture_control_set_ptr->frames_in_sw < (uint32_t)(sequence_control_set_ptr->intra_period_length + 1)) &&
                (picture_control_set_ptr->picture_number % ((sequence_control_set_ptr->intra_period_length + 1)) == 0)) {
                selected_ref_qp++;
            }

            selected_ref_qp = (uint32_t)CLIP3(
                sequence_control_set_ptr->static_config.min_qp_allowed,
                sequence_control_set_ptr->static_config.max_qp_allowed,
                selected_ref_qp);

            queue_entry_index_head_temp = (int32_t)(picture_control_set_ptr->picture_number - encode_context_ptr->hl_rate_control_historgram_queue[encode_context_ptr->hl_rate_control_historgram_queue_head_index]->picture_number);
            queue_entry_index_head_temp += encode_context_ptr->hl_rate_control_historgram_queue_head_index;
            queue_entry_index_head_temp = (queue_entry_index_head_temp > HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1) ?
                queue_entry_index_head_temp - HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH :
                queue_entry_index_head_temp;

            queue_entry_index_temp = queue_entry_index_head_temp;
            {
                hl_rate_control_histogram_ptr_temp = (encode_context_ptr->hl_rate_control_historgram_queue[queue_entry_index_temp]);

                if (hl_rate_control_histogram_ptr_temp->slice_type == I_SLICE)
                    ref_qp_index_temp = context_ptr->qp_scaling_map_I_SLICE[selected_ref_qp];
                else
                    ref_qp_index_temp = context_ptr->qp_scaling_map[hl_rate_control_histogram_ptr_temp->temporal_layer_index][selected_ref_qp];

                ref_qp_index_temp = (uint32_t)CLIP3(
                    sequence_control_set_ptr->static_config.min_qp_allowed,
                    sequence_control_set_ptr->static_config.max_qp_allowed,
                    ref_qp_index_temp);

                hl_rate_control_histogram_ptr_temp->pred_bits_ref_qp[ref_qp_index_temp] = 0;
                rate_control_tables_ptr = &encode_context_ptr->rate_control_tables_array[ref_qp_index_temp];
                sad_bits_array_ptr = rate_control_tables_ptr->sad_bits_array[hl_rate_control_histogram_ptr_temp->temporal_layer_index];
                intra_sad_bits_array_ptr = rate_control_tables_ptr->intra_sad_bits_array[hl_rate_control_histogram_ptr_temp->temporal_layer_index];
                pred_bits_ref_qp = 0;
                num_of_full_sbs = 0;

                if (hl_rate_control_histogram_ptr_temp->slice_type == I_SLICE) {
                    // Loop over block in the frame and calculated the predicted bits at reg QP
                    {
                        unsigned i;
                        uint32_t accum = 0;
                        for (i = 0; i < NUMBER_OF_INTRA_SAD_INTERVALS; ++i)
                            accum += (uint32_t)(hl_rate_control_histogram_ptr_temp->ois_distortion_histogram[i] * intra_sad_bits_array_ptr[i]);
                        pred_bits_ref_qp = accum;
                        num_of_full_sbs = hl_rate_control_histogram_ptr_temp->full_sb_count;
                    }
                    hl_rate_control_histogram_ptr_temp->pred_bits_ref_qp[ref_qp_index_temp] += pred_bits_ref_qp;
                }

                else {
                    {
                        unsigned i;
                        uint32_t accum = 0;
                        for (i = 0; i < NUMBER_OF_SAD_INTERVALS; ++i)
                            accum += (uint32_t)(hl_rate_control_histogram_ptr_temp->me_distortion_histogram[i] * sad_bits_array_ptr[i]);
                        pred_bits_ref_qp = accum;
                        num_of_full_sbs = hl_rate_control_histogram_ptr_temp->full_sb_count;
                    }
                    hl_rate_control_histogram_ptr_temp->pred_bits_ref_qp[ref_qp_index_temp] += pred_bits_ref_qp;
                }

                // Scale for in complete
                //  pred_bits_ref_qp is normalized based on the area because of the LCUs at the picture boundries
                hl_rate_control_histogram_ptr_temp->pred_bits_ref_qp[ref_qp_index_temp] = hl_rate_control_histogram_ptr_temp->pred_bits_ref_qp[ref_qp_index_temp] * (uint64_t)area_in_pixel / (num_of_full_sbs << 12);

                // Store the pred_bits_ref_qp for the first frame in the window to PCS
                picture_control_set_ptr->pred_bits_ref_qp[ref_qp_index_temp] = hl_rate_control_histogram_ptr_temp->pred_bits_ref_qp[ref_qp_index_temp];
            }
        }
        else {
            // Loop over the QPs and find the best QP
            min_la_bit_distance = MAX_UNSIGNED_VALUE;
            qp_search_min = (uint8_t)CLIP3(
                sequence_control_set_ptr->static_config.min_qp_allowed,
                MAX_REF_QP_NUM,//sequence_control_set_ptr->static_config.max_qp_allowed,
                (uint32_t)MAX((int32_t)sequence_control_set_ptr->qp - 40, 0));

            qp_search_max = (uint8_t)CLIP3(
                sequence_control_set_ptr->static_config.min_qp_allowed,
                MAX_REF_QP_NUM,//sequence_control_set_ptr->static_config.max_qp_allowed,
                sequence_control_set_ptr->qp + 40);

            for (ref_qp_table_index = qp_search_min; ref_qp_table_index < qp_search_max; ref_qp_table_index++)
                high_level_rate_control_ptr->pred_bits_ref_qpPerSw[ref_qp_table_index] = 0;
            bit_constraint_per_sw = high_level_rate_control_ptr->bit_constraint_per_sw * picture_control_set_ptr->frames_in_sw / (sequence_control_set_ptr->static_config.look_ahead_distance + 1);

            // Update the target rate for the sliding window based on the status of RC
            if ((context_ptr->extra_bits_gen > (int64_t)(context_ptr->virtual_buffer_size * 10)))
                bit_constraint_per_sw = bit_constraint_per_sw * 130 / 100;
            else if ((context_ptr->extra_bits_gen > (int64_t)(context_ptr->virtual_buffer_size << 3)))
                bit_constraint_per_sw = bit_constraint_per_sw * 120 / 100;
            else if ((context_ptr->extra_bits_gen > (int64_t)(context_ptr->virtual_buffer_size << 2)))
                bit_constraint_per_sw = bit_constraint_per_sw * 110 / 100;
            if ((context_ptr->extra_bits_gen < -(int64_t)(context_ptr->virtual_buffer_size << 3)))
                bit_constraint_per_sw = bit_constraint_per_sw * 80 / 100;
            else if ((context_ptr->extra_bits_gen < -(int64_t)(context_ptr->virtual_buffer_size << 2)))
                bit_constraint_per_sw = bit_constraint_per_sw * 90 / 100;
            // Loop over proper QPs and find the Predicted bits for that QP. Find the QP with the closest total predicted rate to target bits for the sliding window.
            previous_selected_ref_qp = CLIP3(
                qp_search_min + 1,
                qp_search_max - 1,
                previous_selected_ref_qp);
            ref_qp_table_index = previous_selected_ref_qp;
            ref_qp_index = ref_qp_table_index;
            selected_ref_qp_table_index = ref_qp_table_index;
            selected_ref_qp = selected_ref_qp_table_index;
            if (sequence_control_set_ptr->intra_period_length != -1 && picture_control_set_ptr->picture_number % ((sequence_control_set_ptr->intra_period_length + 1)) == 0 &&
                (int32_t)picture_control_set_ptr->frames_in_sw > sequence_control_set_ptr->intra_period_length) {
                best_qp_found = EB_FALSE;
                while (ref_qp_table_index >= qp_search_min && ref_qp_table_index <= qp_search_max && !best_qp_found) {
                    ref_qp_index = CLIP3(
                        sequence_control_set_ptr->static_config.min_qp_allowed,
                        MAX_REF_QP_NUM,//sequence_control_set_ptr->static_config.max_qp_allowed,
                        ref_qp_table_index);
                    high_level_rate_control_ptr->pred_bits_ref_qpPerSw[ref_qp_index] = 0;

                    // Finding the predicted bits for each frame in the sliding window at the reference Qp(s)
                    queue_entry_index_head_temp = (int32_t)(picture_control_set_ptr->picture_number - encode_context_ptr->hl_rate_control_historgram_queue[encode_context_ptr->hl_rate_control_historgram_queue_head_index]->picture_number);
                    queue_entry_index_head_temp += encode_context_ptr->hl_rate_control_historgram_queue_head_index;
                    queue_entry_index_head_temp = (queue_entry_index_head_temp > HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1) ?
                        queue_entry_index_head_temp - HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH :
                        queue_entry_index_head_temp;

                    queue_entry_index_temp = queue_entry_index_head_temp;
                    // This is set to false, so the last frame would go inside the loop
                    end_of_sequence_flag = EB_FALSE;

                    while (!end_of_sequence_flag &&
                        queue_entry_index_temp <= queue_entry_index_head_temp + sequence_control_set_ptr->static_config.look_ahead_distance) {
                        queue_entry_index_temp2 = (queue_entry_index_temp > HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1) ? queue_entry_index_temp - HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH : queue_entry_index_temp;
                        hl_rate_control_histogram_ptr_temp = (encode_context_ptr->hl_rate_control_historgram_queue[queue_entry_index_temp2]);

                        if (hl_rate_control_histogram_ptr_temp->slice_type == I_SLICE)
                            ref_qp_index_temp = context_ptr->qp_scaling_map_I_SLICE[ref_qp_index];
                        else
                            ref_qp_index_temp = context_ptr->qp_scaling_map[hl_rate_control_histogram_ptr_temp->temporal_layer_index][ref_qp_index];

                        ref_qp_index_temp = (uint32_t)CLIP3(
                            sequence_control_set_ptr->static_config.min_qp_allowed,
                            sequence_control_set_ptr->static_config.max_qp_allowed,
                            ref_qp_index_temp);

                        hl_rate_control_histogram_ptr_temp->pred_bits_ref_qp[ref_qp_index_temp] = 0;

                        if (ref_qp_table_index == previous_selected_ref_qp) {
                            eb_block_on_mutex(sequence_control_set_ptr->encode_context_ptr->hl_rate_control_historgram_queue_mutex);
                            hl_rate_control_histogram_ptr_temp->life_count--;
                            eb_release_mutex(sequence_control_set_ptr->encode_context_ptr->hl_rate_control_historgram_queue_mutex);
                        }
                        hl_rate_control_histogram_ptr_temp->pred_bits_ref_qp[ref_qp_index_temp] = predict_bits(
                            encode_context_ptr,
                            hl_rate_control_histogram_ptr_temp,
                            ref_qp_index_temp,
                            area_in_pixel);

                        high_level_rate_control_ptr->pred_bits_ref_qpPerSw[ref_qp_index] += hl_rate_control_histogram_ptr_temp->pred_bits_ref_qp[ref_qp_index_temp];
                        // Store the pred_bits_ref_qp for the first frame in the window to PCS
                        if (queue_entry_index_head_temp == queue_entry_index_temp2)
                            picture_control_set_ptr->pred_bits_ref_qp[ref_qp_index_temp] = hl_rate_control_histogram_ptr_temp->pred_bits_ref_qp[ref_qp_index_temp];

                        end_of_sequence_flag = hl_rate_control_histogram_ptr_temp->end_of_sequence_flag;
                        queue_entry_index_temp++;
                    }

                    if (min_la_bit_distance >= (uint64_t)ABS((int64_t)high_level_rate_control_ptr->pred_bits_ref_qpPerSw[ref_qp_index] - (int64_t)bit_constraint_per_sw)) {
                        min_la_bit_distance = (uint64_t)ABS((int64_t)high_level_rate_control_ptr->pred_bits_ref_qpPerSw[ref_qp_index] - (int64_t)bit_constraint_per_sw);
                        selected_ref_qp_table_index = ref_qp_table_index;
                        selected_ref_qp = ref_qp_index;
                    }
                    else
                        best_qp_found = EB_TRUE;
                    if (ref_qp_table_index == previous_selected_ref_qp) {
                        if (high_level_rate_control_ptr->pred_bits_ref_qpPerSw[ref_qp_index] > bit_constraint_per_sw)
                            qp_step = +1;
                        else
                            qp_step = -1;
                    }
                    ref_qp_table_index = (uint32_t)(ref_qp_table_index + qp_step);
                }

                if (ref_qp_index == sequence_control_set_ptr->static_config.max_qp_allowed && high_level_rate_control_ptr->pred_bits_ref_qpPerSw[ref_qp_index] > bit_constraint_per_sw) {
                    delta_qp = (int)((high_level_rate_control_ptr->pred_bits_ref_qpPerSw[ref_qp_index] - bit_constraint_per_sw) * 100 / (high_level_rate_control_ptr->pred_bits_ref_qpPerSw[ref_qp_index - 1] - high_level_rate_control_ptr->pred_bits_ref_qpPerSw[ref_qp_index]));
                    delta_qp = (delta_qp + 50) / 100;
                }
            }
        }
#if RC_UPDATE_TARGET_RATE
        selected_org_ref_qp = selected_ref_qp;
        if (sequence_control_set_ptr->intra_period_length != -1 && picture_control_set_ptr->picture_number % ((sequence_control_set_ptr->intra_period_length + 1)) == 0 &&
            (int32_t)picture_control_set_ptr->frames_in_sw > sequence_control_set_ptr->intra_period_length) {
            if (picture_control_set_ptr->picture_number > 0)
                picture_control_set_ptr->intra_selected_org_qp = (uint8_t)selected_ref_qp;
            ref_qp_index = selected_ref_qp;
            high_level_rate_control_ptr->pred_bits_ref_qpPerSw[ref_qp_index] = 0;

            if (high_level_rate_control_ptr->pred_bits_ref_qpPerSw[ref_qp_index] == 0) {
                // Finding the predicted bits for each frame in the sliding window at the reference Qp(s)
                //queue_entry_index_temp = encode_context_ptr->hl_rate_control_historgram_queue_head_index;
                queue_entry_index_head_temp = (int32_t)(picture_control_set_ptr->picture_number - encode_context_ptr->hl_rate_control_historgram_queue[encode_context_ptr->hl_rate_control_historgram_queue_head_index]->picture_number);
                queue_entry_index_head_temp += encode_context_ptr->hl_rate_control_historgram_queue_head_index;
                queue_entry_index_head_temp = (queue_entry_index_head_temp > HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1) ?
                    queue_entry_index_head_temp - HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH :
                    queue_entry_index_head_temp;

                queue_entry_index_temp = queue_entry_index_head_temp;

                // This is set to false, so the last frame would go inside the loop
                end_of_sequence_flag = EB_FALSE;

                while (!end_of_sequence_flag &&
                    //queue_entry_index_temp <= encode_context_ptr->hl_rate_control_historgram_queue_head_index+sequence_control_set_ptr->static_config.look_ahead_distance){
                    queue_entry_index_temp <= queue_entry_index_head_temp + sequence_control_set_ptr->static_config.look_ahead_distance) {
                    queue_entry_index_temp2 = (queue_entry_index_temp > HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1) ? queue_entry_index_temp - HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH : queue_entry_index_temp;
                    hl_rate_control_histogram_ptr_temp = (encode_context_ptr->hl_rate_control_historgram_queue[queue_entry_index_temp2]);

                    if (hl_rate_control_histogram_ptr_temp->slice_type == I_SLICE)
                        ref_qp_index_temp = context_ptr->qp_scaling_map_I_SLICE[ref_qp_index];
                    else
                        ref_qp_index_temp = context_ptr->qp_scaling_map[hl_rate_control_histogram_ptr_temp->temporal_layer_index][ref_qp_index];

                    ref_qp_index_temp = (uint32_t)CLIP3(
                        sequence_control_set_ptr->static_config.min_qp_allowed,
                        sequence_control_set_ptr->static_config.max_qp_allowed,
                        ref_qp_index_temp);

                    hl_rate_control_histogram_ptr_temp->pred_bits_ref_qp[ref_qp_index_temp] = predict_bits(
                        encode_context_ptr,
                        hl_rate_control_histogram_ptr_temp,
                        ref_qp_index_temp,
                        area_in_pixel);

                    high_level_rate_control_ptr->pred_bits_ref_qpPerSw[ref_qp_index] += hl_rate_control_histogram_ptr_temp->pred_bits_ref_qp[ref_qp_index_temp];
                    // Store the pred_bits_ref_qp for the first frame in the window to PCS
                    //  if(encode_context_ptr->hl_rate_control_historgram_queue_head_index == queue_entry_index_temp2)
                    if (queue_entry_index_head_temp == queue_entry_index_temp2)
                        picture_control_set_ptr->pred_bits_ref_qp[ref_qp_index_temp] = hl_rate_control_histogram_ptr_temp->pred_bits_ref_qp[ref_qp_index_temp];

                    end_of_sequence_flag = hl_rate_control_histogram_ptr_temp->end_of_sequence_flag;
                    queue_entry_index_temp++;
                }
            }
        }
#endif
        picture_control_set_ptr->tables_updated = tables_updated;

        // Looping over the window to find the percentage of bit allocation in each layer
        if ((sequence_control_set_ptr->intra_period_length != -1) &&
            ((int32_t)picture_control_set_ptr->frames_in_sw > sequence_control_set_ptr->intra_period_length) &&
            ((int32_t)picture_control_set_ptr->frames_in_sw > sequence_control_set_ptr->intra_period_length)) {
            if (picture_control_set_ptr->picture_number % ((sequence_control_set_ptr->intra_period_length + 1)) == 0) {
                queue_entry_index_head_temp = (int32_t)(picture_control_set_ptr->picture_number - encode_context_ptr->hl_rate_control_historgram_queue[encode_context_ptr->hl_rate_control_historgram_queue_head_index]->picture_number);
                queue_entry_index_head_temp += encode_context_ptr->hl_rate_control_historgram_queue_head_index;
                queue_entry_index_head_temp = (queue_entry_index_head_temp > HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1) ?
                    queue_entry_index_head_temp - HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH :
                    queue_entry_index_head_temp;

                queue_entry_index_temp = queue_entry_index_head_temp;

                // This is set to false, so the last frame would go inside the loop
                end_of_sequence_flag = EB_FALSE;

                while (!end_of_sequence_flag &&
                    queue_entry_index_temp <= queue_entry_index_head_temp + sequence_control_set_ptr->intra_period_length) {
                    queue_entry_index_temp2 = (queue_entry_index_temp > HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1) ? queue_entry_index_temp - HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH : queue_entry_index_temp;
                    hl_rate_control_histogram_ptr_temp = (encode_context_ptr->hl_rate_control_historgram_queue[queue_entry_index_temp2]);

                    if (hl_rate_control_histogram_ptr_temp->slice_type == I_SLICE)
                        ref_qp_index_temp = context_ptr->qp_scaling_map_I_SLICE[selected_ref_qp];
                    else
                        ref_qp_index_temp = context_ptr->qp_scaling_map[hl_rate_control_histogram_ptr_temp->temporal_layer_index][selected_ref_qp];

                    ref_qp_index_temp = (uint32_t)CLIP3(
                        sequence_control_set_ptr->static_config.min_qp_allowed,
                        sequence_control_set_ptr->static_config.max_qp_allowed,
                        ref_qp_index_temp);

                    picture_control_set_ptr->total_bits_per_gop += hl_rate_control_histogram_ptr_temp->pred_bits_ref_qp[ref_qp_index_temp];
                    picture_control_set_ptr->bits_per_sw_per_layer[hl_rate_control_histogram_ptr_temp->temporal_layer_index] += hl_rate_control_histogram_ptr_temp->pred_bits_ref_qp[ref_qp_index_temp];
                    picture_control_set_ptr->percentage_updated = EB_TRUE;

                    end_of_sequence_flag = hl_rate_control_histogram_ptr_temp->end_of_sequence_flag;
                    queue_entry_index_temp++;
                }
                if (picture_control_set_ptr->total_bits_per_gop == 0) {
                    for (temporal_layer_index = 0; temporal_layer_index < EB_MAX_TEMPORAL_LAYERS; temporal_layer_index++)
                        picture_control_set_ptr->bits_per_sw_per_layer[temporal_layer_index] = rate_percentage_layer_array[sequence_control_set_ptr->static_config.hierarchical_levels][temporal_layer_index];
                }
            }
        }
        else {
            for (temporal_layer_index = 0; temporal_layer_index < EB_MAX_TEMPORAL_LAYERS; temporal_layer_index++)
                picture_control_set_ptr->bits_per_sw_per_layer[temporal_layer_index] = rate_percentage_layer_array[sequence_control_set_ptr->static_config.hierarchical_levels][temporal_layer_index];
        }

        // Set the QP
        previous_selected_ref_qp = selected_ref_qp;
        if (picture_control_set_ptr->picture_number > max_coded_poc && picture_control_set_ptr->temporal_layer_index < 2 && !picture_control_set_ptr->end_of_sequence_region) {
            max_coded_poc = picture_control_set_ptr->picture_number;
            max_coded_poc_selected_ref_qp = previous_selected_ref_qp;
            encode_context_ptr->previous_selected_ref_qp = previous_selected_ref_qp;
            encode_context_ptr->max_coded_poc = max_coded_poc;
            encode_context_ptr->max_coded_poc_selected_ref_qp = max_coded_poc_selected_ref_qp;
        }

        if (picture_control_set_ptr->slice_type == I_SLICE)
            picture_control_set_ptr->best_pred_qp = (uint8_t)context_ptr->qp_scaling_map_I_SLICE[selected_ref_qp];
        else
            picture_control_set_ptr->best_pred_qp = (uint8_t)context_ptr->qp_scaling_map[picture_control_set_ptr->temporal_layer_index][selected_ref_qp];

        picture_control_set_ptr->target_bits_best_pred_qp = picture_control_set_ptr->pred_bits_ref_qp[picture_control_set_ptr->best_pred_qp];
        picture_control_set_ptr->best_pred_qp = (uint8_t)CLIP3(
            sequence_control_set_ptr->static_config.min_qp_allowed,
            sequence_control_set_ptr->static_config.max_qp_allowed,
            (uint8_t)((int)picture_control_set_ptr->best_pred_qp + delta_qp));

#if RC_UPDATE_TARGET_RATE
        if (picture_control_set_ptr->picture_number == 0) {
            high_level_rate_control_ptr->prev_intra_selected_ref_qp = selected_ref_qp;
            high_level_rate_control_ptr->prev_intra_org_selected_ref_qp = selected_ref_qp;
        }
        if (sequence_control_set_ptr->intra_period_length != -1) {
            if (picture_control_set_ptr->picture_number % ((sequence_control_set_ptr->intra_period_length + 1)) == 0) {
                high_level_rate_control_ptr->prev_intra_selected_ref_qp = selected_ref_qp;
                high_level_rate_control_ptr->prev_intra_org_selected_ref_qp = selected_org_ref_qp;
            }
        }
#endif
#if RC_PRINTS
        ////if (picture_control_set_ptr->slice_type == 2)
        {
            SVT_LOG("\nTID: %d\t", picture_control_set_ptr->temporal_layer_index);
            SVT_LOG("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t\n",
                picture_control_set_ptr->picture_number,
                picture_control_set_ptr->best_pred_qp,
                selected_ref_qp,
                (int)picture_control_set_ptr->target_bits_best_pred_qp,
                (int)high_level_rate_control_ptr->pred_bits_ref_qpPerSw[selected_ref_qp - 1],
                (int)high_level_rate_control_ptr->pred_bits_ref_qpPerSw[selected_ref_qp],
                (int)high_level_rate_control_ptr->pred_bits_ref_qpPerSw[selected_ref_qp + 1],
                (int)high_level_rate_control_ptr->bit_constraint_per_sw,
                (int)bit_constraint_per_sw/*,
                (int)high_level_rate_control_ptr->virtual_buffer_level*/);
        }
#endif
    }
    eb_release_mutex(sequence_control_set_ptr->encode_context_ptr->rate_table_update_mutex);
}
void frame_level_rc_input_picture_cvbr(
    PictureControlSet               *picture_control_set_ptr,
    SequenceControlSet              *sequence_control_set_ptr,
    RateControlContext              *context_ptr,
    RateControlLayerContext         *rate_control_layer_ptr,
    RateControlIntervalParamContext *rate_control_param_ptr)
{
    RateControlLayerContext *rate_control_layer_temp_ptr;

    // Tiles
    uint32_t                 picture_area_in_pixel;
    uint32_t                 area_in_pixel;

    // LCU Loop variables
    SbParams               *sb_params_ptr;
    uint32_t                 sb_index;
    uint64_t                 temp_qp;
    uint32_t                 area_in_sbs;

    picture_area_in_pixel = sequence_control_set_ptr->seq_header.max_frame_height*sequence_control_set_ptr->seq_header.max_frame_width;

    if (rate_control_layer_ptr->first_frame == 1) {
        rate_control_layer_ptr->first_frame = 0;
        picture_control_set_ptr->parent_pcs_ptr->first_frame_in_temporal_layer = 1;
    }
    else
        picture_control_set_ptr->parent_pcs_ptr->first_frame_in_temporal_layer = 0;
    if (picture_control_set_ptr->slice_type != I_SLICE) {
        if (rate_control_layer_ptr->first_non_intra_frame == 1) {
            rate_control_layer_ptr->first_non_intra_frame = 0;
            picture_control_set_ptr->parent_pcs_ptr->first_non_intra_frame_in_temporal_layer = 1;
        }
        else
            picture_control_set_ptr->parent_pcs_ptr->first_non_intra_frame_in_temporal_layer = 0;
    }
    else
        picture_control_set_ptr->parent_pcs_ptr->first_non_intra_frame_in_temporal_layer = 0;

    picture_control_set_ptr->parent_pcs_ptr->target_bits_rc = 0;

    // ***Rate Control***
    area_in_sbs = 0;
    area_in_pixel = 0;

    for (sb_index = 0; sb_index < picture_control_set_ptr->sb_total_count; ++sb_index) {
        sb_params_ptr = &sequence_control_set_ptr->sb_params_array[sb_index];

        if (sb_params_ptr->is_complete_sb) {
            // add the area of one LCU (64x64=4096) to the area of the tile
            area_in_pixel += 4096;
            area_in_sbs++;
        }
        else {
            // add the area of the LCU to the area of the tile
            area_in_pixel += sb_params_ptr->width * sb_params_ptr->height;
        }
    }
    rate_control_layer_ptr->area_in_pixel = area_in_pixel;

    if (picture_control_set_ptr->parent_pcs_ptr->first_frame_in_temporal_layer || (picture_control_set_ptr->picture_number == rate_control_param_ptr->first_poc)) {
        if (sequence_control_set_ptr->static_config.enable_qp_scaling_flag && (picture_control_set_ptr->picture_number != rate_control_param_ptr->first_poc)) {
            picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                (int32_t)sequence_control_set_ptr->static_config.min_qp_allowed,
                (int32_t)sequence_control_set_ptr->static_config.max_qp_allowed,
                (int32_t)(rate_control_param_ptr->intra_frames_qp + context_ptr->qp_scaling_map[picture_control_set_ptr->temporal_layer_index][rate_control_param_ptr->intra_frames_qp_bef_scal] - context_ptr->qp_scaling_map_I_SLICE[rate_control_param_ptr->intra_frames_qp_bef_scal]));
        }

        if (picture_control_set_ptr->picture_number == 0) {
            rate_control_param_ptr->intra_frames_qp = sequence_control_set_ptr->qp;
            rate_control_param_ptr->intra_frames_qp_bef_scal = (uint8_t)sequence_control_set_ptr->qp;
        }

        if (picture_control_set_ptr->picture_number == rate_control_param_ptr->first_poc) {
            uint32_t temporal_layer_idex;
            rate_control_param_ptr->previous_virtual_buffer_level = context_ptr->virtual_buffer_level_initial_value;
            rate_control_param_ptr->virtual_buffer_level = context_ptr->virtual_buffer_level_initial_value;
            rate_control_param_ptr->extra_ap_bit_ratio_i = 0;
            if (picture_control_set_ptr->parent_pcs_ptr->end_of_sequence_region) {
                rate_control_param_ptr->last_poc = MAX(rate_control_param_ptr->first_poc + picture_control_set_ptr->parent_pcs_ptr->frames_in_sw - 1, rate_control_param_ptr->first_poc);
                rate_control_param_ptr->last_gop = EB_TRUE;
            }

            for (temporal_layer_idex = 0; temporal_layer_idex < EB_MAX_TEMPORAL_LAYERS; temporal_layer_idex++) {
                rate_control_layer_temp_ptr = rate_control_param_ptr->rate_control_layer_array[temporal_layer_idex];
                rate_control_layer_reset(
                    rate_control_layer_temp_ptr,
                    picture_control_set_ptr,
                    context_ptr,
                    picture_area_in_pixel,
                    rate_control_param_ptr->was_used);
            }
        }

        picture_control_set_ptr->parent_pcs_ptr->sad_me = 0;
        // Finding the QP of the Intra frame by using variance tables
        if (picture_control_set_ptr->slice_type == I_SLICE) {
            uint32_t         selected_ref_qp;

            if (sequence_control_set_ptr->static_config.look_ahead_distance == 0)
                printf("ERROR: LAD=0 is not supported\n");
            else {
                selected_ref_qp = picture_control_set_ptr->parent_pcs_ptr->best_pred_qp;
                picture_control_set_ptr->picture_qp = (uint8_t)selected_ref_qp;
                picture_control_set_ptr->parent_pcs_ptr->calculated_qp = picture_control_set_ptr->picture_qp;
            }

            picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                sequence_control_set_ptr->static_config.min_qp_allowed,
                sequence_control_set_ptr->static_config.max_qp_allowed,
                picture_control_set_ptr->picture_qp);
        }
        else {
            // LCU Loop
            for (sb_index = 0; sb_index < picture_control_set_ptr->sb_total_count; ++sb_index) {
                sb_params_ptr = &sequence_control_set_ptr->sb_params_array[sb_index];

                if (sb_params_ptr->is_complete_sb)
                    picture_control_set_ptr->parent_pcs_ptr->sad_me += picture_control_set_ptr->parent_pcs_ptr->rc_me_distortion[sb_index];
            }

            //  tileSadMe is normalized based on the area because of the LCUs at the tile boundries
            picture_control_set_ptr->parent_pcs_ptr->sad_me = MAX((picture_control_set_ptr->parent_pcs_ptr->sad_me*rate_control_layer_ptr->area_in_pixel / (area_in_sbs << 12)), 1);

            // totalSquareMad has RC_PRECISION precision
            picture_control_set_ptr->parent_pcs_ptr->sad_me <<= RC_PRECISION;
        }

        temp_qp = picture_control_set_ptr->picture_qp;

        if (picture_control_set_ptr->picture_number == rate_control_param_ptr->first_poc) {
            uint32_t temporal_layer_idex;
            for (temporal_layer_idex = 0; temporal_layer_idex < EB_MAX_TEMPORAL_LAYERS; temporal_layer_idex++) {
                rate_control_layer_temp_ptr = rate_control_param_ptr->rate_control_layer_array[temporal_layer_idex];
                rate_control_layer_reset_part2(
                    context_ptr,
                    rate_control_layer_temp_ptr,
                    picture_control_set_ptr);
            }
        }

        if (picture_control_set_ptr->picture_number == 0) {
            context_ptr->base_layer_frames_avg_qp = picture_control_set_ptr->picture_qp + 1;
            context_ptr->base_layer_intra_frames_avg_qp = picture_control_set_ptr->picture_qp;
        }
    }
    else {
        picture_control_set_ptr->parent_pcs_ptr->sad_me = 0;

        HighLevelRateControlContext *high_level_rate_control_ptr = context_ptr->high_level_rate_control_ptr;
        EncodeContext               *encode_context_ptr = sequence_control_set_ptr->encode_context_ptr;
        HlRateControlHistogramEntry *hl_rate_control_histogram_ptr_temp;
        // Queue variables
        uint32_t                     queue_entry_index_temp;
        uint32_t                     queue_entry_index_temp2;
        uint32_t                     queue_entry_index_head_temp;

        uint64_t                     min_la_bit_distance;
        uint32_t                     selected_ref_qp_table_index;
        uint32_t                     selected_ref_qp;
        uint32_t                     previous_selected_ref_qp = encode_context_ptr->previous_selected_ref_qp;

        uint32_t                     ref_qp_index;
        uint32_t                     ref_qp_index_temp;
        uint32_t                     ref_qp_table_index;

        uint32_t                     qp_search_min;
        uint32_t                     qp_search_max;
        int32_t                      qp_step = 1;
        EbBool                      best_qp_found;

        uint64_t                     bit_constraint_per_sw = 0;
        EbBool                      end_of_sequence_flag = EB_TRUE;

        // Loop over the QPs and find the best QP
        min_la_bit_distance = MAX_UNSIGNED_VALUE;
        qp_search_min = (uint8_t)CLIP3(
            sequence_control_set_ptr->static_config.min_qp_allowed,
            MAX_REF_QP_NUM,//sequence_control_set_ptr->static_config.max_qp_allowed,
            (uint32_t)MAX((int32_t)sequence_control_set_ptr->qp - 40, 0));

        qp_search_max = (uint8_t)CLIP3(
            sequence_control_set_ptr->static_config.min_qp_allowed,
            MAX_REF_QP_NUM,
            sequence_control_set_ptr->qp + 40);

        for (ref_qp_table_index = qp_search_min; ref_qp_table_index < qp_search_max; ref_qp_table_index++)
            high_level_rate_control_ptr->pred_bits_ref_qpPerSw[ref_qp_table_index] = 0;
        // Finding the predicted bits for each frame in the sliding window at the reference Qp(s)
        ///queue_entry_index_head_temp = (int32_t)(picture_control_set_ptr->picture_number - encode_context_ptr->hl_rate_control_historgram_queue[encode_context_ptr->hl_rate_control_historgram_queue_head_index]->picture_number);
        queue_entry_index_head_temp = (int32_t)(rate_control_param_ptr->first_poc - encode_context_ptr->hl_rate_control_historgram_queue[encode_context_ptr->hl_rate_control_historgram_queue_head_index]->picture_number);
        queue_entry_index_head_temp += encode_context_ptr->hl_rate_control_historgram_queue_head_index;
        queue_entry_index_head_temp = (queue_entry_index_head_temp > HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1) ?
            queue_entry_index_head_temp - HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH :
            queue_entry_index_head_temp;

        if (picture_control_set_ptr->parent_pcs_ptr->picture_number + picture_control_set_ptr->parent_pcs_ptr->frames_in_sw > rate_control_param_ptr->first_poc + sequence_control_set_ptr->static_config.intra_period_length + 1)
            bit_constraint_per_sw = high_level_rate_control_ptr->bit_constraint_per_sw;
        else
            bit_constraint_per_sw = high_level_rate_control_ptr->bit_constraint_per_sw * encode_context_ptr->hl_rate_control_historgram_queue[queue_entry_index_head_temp]->frames_in_sw / (sequence_control_set_ptr->static_config.look_ahead_distance + 1);

        // Loop over proper QPs and find the Predicted bits for that QP. Find the QP with the closest total predicted rate to target bits for the sliding window.
        previous_selected_ref_qp = CLIP3(
            qp_search_min + 1,
            qp_search_max - 1,
            previous_selected_ref_qp);
        ref_qp_table_index = previous_selected_ref_qp;
        ref_qp_index = ref_qp_table_index;
        selected_ref_qp_table_index = ref_qp_table_index;
        selected_ref_qp = selected_ref_qp_table_index;
        best_qp_found = EB_FALSE;
        while (ref_qp_table_index >= qp_search_min && ref_qp_table_index <= qp_search_max && !best_qp_found) {
            ref_qp_index = CLIP3(
                sequence_control_set_ptr->static_config.min_qp_allowed,
                MAX_REF_QP_NUM,
                ref_qp_table_index);
            high_level_rate_control_ptr->pred_bits_ref_qpPerSw[ref_qp_index] = 0;

            queue_entry_index_temp = queue_entry_index_head_temp;
            // This is set to false, so the last frame would go inside the loop
            end_of_sequence_flag = EB_FALSE;

            while (!end_of_sequence_flag &&
                //queue_entry_index_temp <= queue_entry_index_head_temp + sequence_control_set_ptr->static_config.look_ahead_distance) {
                queue_entry_index_temp <= queue_entry_index_head_temp + encode_context_ptr->hl_rate_control_historgram_queue[queue_entry_index_head_temp]->frames_in_sw - 1) {
                queue_entry_index_temp2 = (queue_entry_index_temp > HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1) ? queue_entry_index_temp - HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH : queue_entry_index_temp;
                hl_rate_control_histogram_ptr_temp = (encode_context_ptr->hl_rate_control_historgram_queue[queue_entry_index_temp2]);

                if (hl_rate_control_histogram_ptr_temp->slice_type == I_SLICE)
                    ref_qp_index_temp = context_ptr->qp_scaling_map_I_SLICE[ref_qp_index];
                else
                    ref_qp_index_temp = context_ptr->qp_scaling_map[hl_rate_control_histogram_ptr_temp->temporal_layer_index][ref_qp_index];

                ref_qp_index_temp = (uint32_t)CLIP3(
                    sequence_control_set_ptr->static_config.min_qp_allowed,
                    sequence_control_set_ptr->static_config.max_qp_allowed,
                    ref_qp_index_temp);

                hl_rate_control_histogram_ptr_temp->pred_bits_ref_qp[ref_qp_index_temp] = 0;

                if (ref_qp_table_index == previous_selected_ref_qp) {
                    eb_block_on_mutex(sequence_control_set_ptr->encode_context_ptr->hl_rate_control_historgram_queue_mutex);
                    hl_rate_control_histogram_ptr_temp->life_count--;
                    eb_release_mutex(sequence_control_set_ptr->encode_context_ptr->hl_rate_control_historgram_queue_mutex);
                }

                hl_rate_control_histogram_ptr_temp->pred_bits_ref_qp[ref_qp_index_temp] = predict_bits(
                    encode_context_ptr,
                    hl_rate_control_histogram_ptr_temp,
                    ref_qp_index_temp,
                    area_in_pixel);

                high_level_rate_control_ptr->pred_bits_ref_qpPerSw[ref_qp_index] += hl_rate_control_histogram_ptr_temp->pred_bits_ref_qp[ref_qp_index_temp];
                // Store the pred_bits_ref_qp for the first frame in the window to PCS
                if (queue_entry_index_head_temp == queue_entry_index_temp2)
                    picture_control_set_ptr->parent_pcs_ptr->pred_bits_ref_qp[ref_qp_index_temp] = hl_rate_control_histogram_ptr_temp->pred_bits_ref_qp[ref_qp_index_temp];

                end_of_sequence_flag = hl_rate_control_histogram_ptr_temp->end_of_sequence_flag;
                queue_entry_index_temp++;
            }

            if (min_la_bit_distance >= (uint64_t)ABS((int64_t)high_level_rate_control_ptr->pred_bits_ref_qpPerSw[ref_qp_index] - (int64_t)bit_constraint_per_sw)) {
                min_la_bit_distance = (uint64_t)ABS((int64_t)high_level_rate_control_ptr->pred_bits_ref_qpPerSw[ref_qp_index] - (int64_t)bit_constraint_per_sw);
                selected_ref_qp_table_index = ref_qp_table_index;
                selected_ref_qp = ref_qp_index;
            }
            else
                best_qp_found = EB_TRUE;
            if (ref_qp_table_index == previous_selected_ref_qp) {
                if (high_level_rate_control_ptr->pred_bits_ref_qpPerSw[ref_qp_index] > bit_constraint_per_sw)
                    qp_step = +1;
                else
                    qp_step = -1;
            }
            ref_qp_table_index = (uint32_t)(ref_qp_table_index + qp_step);
        }

        int delta_qp = 0;
        if (ref_qp_index == sequence_control_set_ptr->static_config.max_qp_allowed && high_level_rate_control_ptr->pred_bits_ref_qpPerSw[ref_qp_index] > bit_constraint_per_sw) {
            delta_qp = (int)((high_level_rate_control_ptr->pred_bits_ref_qpPerSw[ref_qp_index] - bit_constraint_per_sw) * 100 / (high_level_rate_control_ptr->pred_bits_ref_qpPerSw[ref_qp_index - 1] - high_level_rate_control_ptr->pred_bits_ref_qpPerSw[ref_qp_index]));
            delta_qp = (delta_qp + 50) / 100;
        }

        if (picture_control_set_ptr->slice_type == I_SLICE)
            picture_control_set_ptr->parent_pcs_ptr->best_pred_qp = (uint8_t)context_ptr->qp_scaling_map_I_SLICE[selected_ref_qp];
        else
            picture_control_set_ptr->parent_pcs_ptr->best_pred_qp = (uint8_t)context_ptr->qp_scaling_map[picture_control_set_ptr->temporal_layer_index][selected_ref_qp];

        picture_control_set_ptr->parent_pcs_ptr->best_pred_qp = (uint8_t)CLIP3(
            sequence_control_set_ptr->static_config.min_qp_allowed,
            sequence_control_set_ptr->static_config.max_qp_allowed,
            (uint8_t)((int)picture_control_set_ptr->parent_pcs_ptr->best_pred_qp + delta_qp));

        // if the pixture is an I slice, for now we set the QP as the QP of the previous frame
        if (picture_control_set_ptr->slice_type == I_SLICE) {
            uint32_t         selected_ref_qp;

            if (sequence_control_set_ptr->static_config.look_ahead_distance == 0)
                printf("ERROR: LAD=0 is not supported\n");
            else {
                selected_ref_qp = picture_control_set_ptr->parent_pcs_ptr->best_pred_qp;
                picture_control_set_ptr->picture_qp = (uint8_t)selected_ref_qp;
                picture_control_set_ptr->parent_pcs_ptr->calculated_qp = picture_control_set_ptr->picture_qp;
            }

            picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                sequence_control_set_ptr->static_config.min_qp_allowed,
                sequence_control_set_ptr->static_config.max_qp_allowed,
                picture_control_set_ptr->picture_qp);

            temp_qp = picture_control_set_ptr->picture_qp;
        }

        else { // Not an I slice
            // combining the target rate from initial RC and frame level RC
            if (sequence_control_set_ptr->static_config.look_ahead_distance != 0) {
                picture_control_set_ptr->parent_pcs_ptr->target_bits_rc = rate_control_layer_ptr->bit_constraint;
                rate_control_layer_ptr->ec_bit_constraint = (rate_control_layer_ptr->alpha * picture_control_set_ptr->parent_pcs_ptr->target_bits_best_pred_qp +
                    ((1 << RC_PRECISION) - rate_control_layer_ptr->alpha) * picture_control_set_ptr->parent_pcs_ptr->target_bits_rc + RC_PRECISION_OFFSET) >> RC_PRECISION;

                rate_control_layer_ptr->ec_bit_constraint = (uint64_t)MAX((int64_t)rate_control_layer_ptr->ec_bit_constraint - (int64_t)rate_control_layer_ptr->dif_total_and_ec_bits, 1);

                picture_control_set_ptr->parent_pcs_ptr->target_bits_rc = rate_control_layer_ptr->ec_bit_constraint;
            }

            // LCU Loop
            for (sb_index = 0; sb_index < picture_control_set_ptr->sb_total_count; ++sb_index) {
                sb_params_ptr = &sequence_control_set_ptr->sb_params_array[sb_index];

                if (sb_params_ptr->is_complete_sb)
                    picture_control_set_ptr->parent_pcs_ptr->sad_me += picture_control_set_ptr->parent_pcs_ptr->rc_me_distortion[sb_index];
            }

            //  tileSadMe is normalized based on the area because of the LCUs at the tile boundries
            picture_control_set_ptr->parent_pcs_ptr->sad_me = MAX((picture_control_set_ptr->parent_pcs_ptr->sad_me*rate_control_layer_ptr->area_in_pixel / (area_in_sbs << 12)), 1);
            picture_control_set_ptr->parent_pcs_ptr->sad_me <<= RC_PRECISION;
            if (rate_control_layer_ptr->area_in_pixel > 0)
                rate_control_layer_ptr->total_mad = MAX((picture_control_set_ptr->parent_pcs_ptr->sad_me / rate_control_layer_ptr->area_in_pixel), 1);
            if (!rate_control_layer_ptr->feedback_arrived)
                rate_control_layer_ptr->previous_frame_distortion_me = picture_control_set_ptr->parent_pcs_ptr->sad_me;
            {
                uint64_t qp_calc_temp1, qp_calc_temp2, qp_calc_temp3;

                qp_calc_temp1 = picture_control_set_ptr->parent_pcs_ptr->sad_me *rate_control_layer_ptr->total_mad;
                qp_calc_temp2 =
                    MAX((int64_t)(rate_control_layer_ptr->ec_bit_constraint << (2 * RC_PRECISION)) - (int64_t)rate_control_layer_ptr->c_coeff*(int64_t)rate_control_layer_ptr->area_in_pixel,
                    (int64_t)(rate_control_layer_ptr->ec_bit_constraint << (2 * RC_PRECISION - 2)));

                // This is a more complex but with higher precision implementation
                if (qp_calc_temp1 > qp_calc_temp2)
                    qp_calc_temp3 = (uint64_t)((qp_calc_temp1 / qp_calc_temp2)*rate_control_layer_ptr->k_coeff);
                else
                    qp_calc_temp3 = (uint64_t)(qp_calc_temp1*rate_control_layer_ptr->k_coeff / qp_calc_temp2);
                temp_qp = (uint64_t)(log2f_high_precision(MAX(((qp_calc_temp3 + RC_PRECISION_OFFSET) >> RC_PRECISION)*((qp_calc_temp3 + RC_PRECISION_OFFSET) >> RC_PRECISION)*((qp_calc_temp3 + RC_PRECISION_OFFSET) >> RC_PRECISION), 1), RC_PRECISION));

                rate_control_layer_ptr->calculated_frame_qp = (uint8_t)(CLIP3(1, 63, (uint32_t)(temp_qp + RC_PRECISION_OFFSET) >> RC_PRECISION));
                picture_control_set_ptr->parent_pcs_ptr->calculated_qp = (uint8_t)(CLIP3(1, 63, (uint32_t)(temp_qp + RC_PRECISION_OFFSET) >> RC_PRECISION));
            }

            temp_qp += rate_control_layer_ptr->delta_qp_fraction;
            picture_control_set_ptr->picture_qp = (uint8_t)((temp_qp + RC_PRECISION_OFFSET) >> RC_PRECISION);
            // Use the QP of HLRC instead of calculated one in FLRC
            if (picture_control_set_ptr->parent_pcs_ptr->hierarchical_levels > 1) {
                picture_control_set_ptr->picture_qp = picture_control_set_ptr->parent_pcs_ptr->best_pred_qp;
                picture_control_set_ptr->parent_pcs_ptr->calculated_qp = picture_control_set_ptr->parent_pcs_ptr->best_pred_qp;
            }
        }
        if (picture_control_set_ptr->parent_pcs_ptr->first_non_intra_frame_in_temporal_layer && picture_control_set_ptr->temporal_layer_index == 0 && picture_control_set_ptr->slice_type != I_SLICE)
            picture_control_set_ptr->picture_qp = (uint8_t)(rate_control_param_ptr->intra_frames_qp + context_ptr->qp_scaling_map[picture_control_set_ptr->temporal_layer_index][rate_control_param_ptr->intra_frames_qp_bef_scal] - context_ptr->qp_scaling_map_I_SLICE[rate_control_param_ptr->intra_frames_qp_bef_scal]);
        if (!rate_control_layer_ptr->feedback_arrived && picture_control_set_ptr->slice_type != I_SLICE) {
            picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                (int32_t)sequence_control_set_ptr->static_config.min_qp_allowed,
                (int32_t)sequence_control_set_ptr->static_config.max_qp_allowed,
                (int32_t)(rate_control_param_ptr->intra_frames_qp + context_ptr->qp_scaling_map[picture_control_set_ptr->temporal_layer_index][rate_control_param_ptr->intra_frames_qp_bef_scal] - context_ptr->qp_scaling_map_I_SLICE[rate_control_param_ptr->intra_frames_qp_bef_scal]));
        }

        if (picture_control_set_ptr->parent_pcs_ptr->end_of_sequence_region) {
            if (rate_control_param_ptr->virtual_buffer_level > context_ptr->vb_fill_threshold2 << 2)
                picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp + (uint8_t)THRESHOLD2QPINCREASE + 4;
            else if (rate_control_param_ptr->virtual_buffer_level > context_ptr->vb_fill_threshold2 << 1)
                picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp + (uint8_t)THRESHOLD2QPINCREASE + 3;
            else if (rate_control_param_ptr->virtual_buffer_level > context_ptr->vb_fill_threshold2)
                picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp + (uint8_t)THRESHOLD2QPINCREASE + 2;
            else if (rate_control_param_ptr->virtual_buffer_level > context_ptr->vb_fill_threshold1 &&
                rate_control_param_ptr->virtual_buffer_level < context_ptr->vb_fill_threshold2) {
                picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp + (uint8_t)THRESHOLD1QPINCREASE + 2;
            }
        }
        else {
            //if (rate_control_param_ptr->virtual_buffer_level > context_ptr->vb_fill_threshold2 << 2){
            if (rate_control_param_ptr->virtual_buffer_level > context_ptr->vb_fill_threshold2 + (int64_t)(context_ptr->virtual_buffer_size * 2 / 3))
                picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp + (uint8_t)THRESHOLD2QPINCREASE + 2;
            //else if (rate_control_param_ptr->virtual_buffer_level > context_ptr->vb_fill_threshold2 << 1){
            else if (rate_control_param_ptr->virtual_buffer_level > context_ptr->vb_fill_threshold2 + (int64_t)(context_ptr->virtual_buffer_size / 3))
                picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp + (uint8_t)THRESHOLD2QPINCREASE + 1;
            else if (rate_control_param_ptr->virtual_buffer_level > context_ptr->vb_fill_threshold2)
                picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp + (uint8_t)THRESHOLD2QPINCREASE + 1;
            else if (rate_control_param_ptr->virtual_buffer_level > context_ptr->vb_fill_threshold1 &&
                rate_control_param_ptr->virtual_buffer_level < context_ptr->vb_fill_threshold2) {
                picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp + (uint8_t)THRESHOLD1QPINCREASE;
            }
        }
        if (picture_control_set_ptr->parent_pcs_ptr->end_of_sequence_region) {
            if (rate_control_param_ptr->virtual_buffer_level < -(context_ptr->vb_fill_threshold2 << 2))
                picture_control_set_ptr->picture_qp = (uint8_t)MAX((int32_t)picture_control_set_ptr->picture_qp - (int32_t)THRESHOLD2QPINCREASE - 2, 0);
            else if (rate_control_param_ptr->virtual_buffer_level < -(context_ptr->vb_fill_threshold2 << 1))
                picture_control_set_ptr->picture_qp = (uint8_t)MAX((int32_t)picture_control_set_ptr->picture_qp - (int32_t)THRESHOLD2QPINCREASE - 1, 0);
            else if (rate_control_param_ptr->virtual_buffer_level < 0)
                picture_control_set_ptr->picture_qp = (uint8_t)MAX((int32_t)picture_control_set_ptr->picture_qp - (int32_t)THRESHOLD2QPINCREASE, 0);
        }
        else {
            if (rate_control_param_ptr->virtual_buffer_level < -(context_ptr->vb_fill_threshold2 << 2))
                picture_control_set_ptr->picture_qp = (uint8_t)MAX((int32_t)picture_control_set_ptr->picture_qp - (int32_t)THRESHOLD2QPINCREASE - 1, 0);
            else if (rate_control_param_ptr->virtual_buffer_level < -context_ptr->vb_fill_threshold2)
                picture_control_set_ptr->picture_qp = (uint8_t)MAX((int32_t)picture_control_set_ptr->picture_qp - (int32_t)THRESHOLD2QPINCREASE, 0);
        }

        uint32_t ref_qp;
        if ((int32_t)picture_control_set_ptr->temporal_layer_index == 0 && picture_control_set_ptr->slice_type != I_SLICE) {
            if (picture_control_set_ptr->ref_slice_type_array[0][0] == I_SLICE) {
                /*    picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                        (uint32_t)picture_control_set_ptr->ref_pic_qp_array[0],
                        (uint32_t)picture_control_set_ptr->picture_qp,
                        picture_control_set_ptr->picture_qp);*/
            }
            else {
                picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                    (uint32_t)MAX((int32_t)picture_control_set_ptr->ref_pic_qp_array[0][0] - 1, 0),
                    (uint32_t)picture_control_set_ptr->ref_pic_qp_array[0][0] + 3,
                    picture_control_set_ptr->picture_qp);
            }
        }
        else {
            ref_qp = 0;
            if (picture_control_set_ptr->ref_slice_type_array[0][0] != I_SLICE)
                ref_qp = MAX(ref_qp, picture_control_set_ptr->ref_pic_qp_array[0][0]);
            if ((picture_control_set_ptr->slice_type == B_SLICE) && (picture_control_set_ptr->ref_slice_type_array[1][0] != I_SLICE))
                ref_qp = MAX(ref_qp, picture_control_set_ptr->ref_pic_qp_array[1][0]);
            if (ref_qp > 0) {
                picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                    (uint32_t)ref_qp - 1,
                    picture_control_set_ptr->picture_qp,
                    picture_control_set_ptr->picture_qp);
            }
        }
        // limiting the QP between min Qp allowed and max Qp allowed
        picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
            sequence_control_set_ptr->static_config.min_qp_allowed,
            sequence_control_set_ptr->static_config.max_qp_allowed,
            picture_control_set_ptr->picture_qp);

        rate_control_layer_ptr->delta_qp_fraction = CLIP3(-RC_PRECISION_OFFSET, RC_PRECISION_OFFSET, -((int64_t)temp_qp - (int64_t)(picture_control_set_ptr->picture_qp << RC_PRECISION)));

        if (picture_control_set_ptr->parent_pcs_ptr->sad_me == rate_control_layer_ptr->previous_frame_distortion_me &&
            (rate_control_layer_ptr->previous_frame_distortion_me != 0))
            rate_control_layer_ptr->same_distortion_count++;
        else
            rate_control_layer_ptr->same_distortion_count = 0;
    }

    rate_control_layer_ptr->previous_c_coeff = rate_control_layer_ptr->c_coeff;
    rate_control_layer_ptr->previous_k_coeff = rate_control_layer_ptr->k_coeff;
    rate_control_layer_ptr->previous_calculated_frame_qp = rate_control_layer_ptr->calculated_frame_qp;
}

void frame_level_rc_feedback_picture_cvbr(
    PictureParentControlSet *parentpicture_control_set_ptr,
    SequenceControlSet      *sequence_control_set_ptr,
    RateControlContext      *context_ptr)
{
    RateControlLayerContext             *rate_control_layer_temp_ptr;
    RateControlIntervalParamContext     *rate_control_param_ptr;
    RateControlLayerContext             *rate_control_layer_ptr;
    // LCU Loop variables
    uint32_t                       slice_num;
    uint64_t                       previous_frame_bit_actual;

    if (sequence_control_set_ptr->intra_period_length == -1)
        rate_control_param_ptr = context_ptr->rate_control_param_queue[0];
    else {
        uint32_t interval_index_temp = 0;
        while ((!(parentpicture_control_set_ptr->picture_number >= context_ptr->rate_control_param_queue[interval_index_temp]->first_poc &&
            parentpicture_control_set_ptr->picture_number <= context_ptr->rate_control_param_queue[interval_index_temp]->last_poc)) &&
            (interval_index_temp < PARALLEL_GOP_MAX_NUMBER)) {
            interval_index_temp++;
        }
        CHECK_REPORT_ERROR(
            interval_index_temp != PARALLEL_GOP_MAX_NUMBER,
            sequence_control_set_ptr->encode_context_ptr->app_callback_ptr,
            EB_ENC_RC_ERROR2);
        rate_control_param_ptr = context_ptr->rate_control_param_queue[interval_index_temp];
    }

    rate_control_layer_ptr = rate_control_param_ptr->rate_control_layer_array[parentpicture_control_set_ptr->temporal_layer_index];

    rate_control_layer_ptr->max_qp = 0;

    rate_control_layer_ptr->feedback_arrived = EB_TRUE;
    rate_control_layer_ptr->max_qp = MAX(rate_control_layer_ptr->max_qp, parentpicture_control_set_ptr->picture_qp);

    rate_control_layer_ptr->previous_frame_qp = parentpicture_control_set_ptr->picture_qp;
    rate_control_layer_ptr->previous_frame_bit_actual = parentpicture_control_set_ptr->total_num_bits;
    if (parentpicture_control_set_ptr->quantized_coeff_num_bits == 0)
        parentpicture_control_set_ptr->quantized_coeff_num_bits = 1;
    rate_control_layer_ptr->previous_framequantized_coeff_bit_actual = parentpicture_control_set_ptr->quantized_coeff_num_bits;

    // Setting Critical states for adjusting the averaging weights on C and K
    if ((parentpicture_control_set_ptr->sad_me > (3 * rate_control_layer_ptr->previous_frame_distortion_me) >> 1) &&
        (rate_control_layer_ptr->previous_frame_distortion_me != 0)) {
        rate_control_layer_ptr->critical_states = 3;
    }
    else if (rate_control_layer_ptr->critical_states)
        rate_control_layer_ptr->critical_states--;
    else
        rate_control_layer_ptr->critical_states = 0;
    if (parentpicture_control_set_ptr->slice_type != I_SLICE) {
        // Updating c_coeff
        rate_control_layer_ptr->c_coeff = (((int64_t)rate_control_layer_ptr->previous_frame_bit_actual - (int64_t)rate_control_layer_ptr->previous_framequantized_coeff_bit_actual) << (2 * RC_PRECISION))
            / rate_control_layer_ptr->area_in_pixel;
        rate_control_layer_ptr->c_coeff = MAX(rate_control_layer_ptr->c_coeff, 1);

        // Updating k_coeff
        if ((parentpicture_control_set_ptr->sad_me + RC_PRECISION_OFFSET) >> RC_PRECISION > 5) {
            {
                uint64_t test1, test2, test3;
                test1 = rate_control_layer_ptr->previous_framequantized_coeff_bit_actual*(two_to_power_qp_over_three[parentpicture_control_set_ptr->picture_qp]);
                test2 = MAX(parentpicture_control_set_ptr->sad_me / rate_control_layer_ptr->area_in_pixel, 1);
                test3 = test1 * 65536 / test2 * 65536 / parentpicture_control_set_ptr->sad_me;

                rate_control_layer_ptr->k_coeff = test3;
            }
        }

        if (rate_control_layer_ptr->critical_states) {
            rate_control_layer_ptr->k_coeff = (8 * rate_control_layer_ptr->k_coeff + 8 * rate_control_layer_ptr->previous_k_coeff + 8) >> 4;
            rate_control_layer_ptr->c_coeff = (8 * rate_control_layer_ptr->c_coeff + 8 * rate_control_layer_ptr->previous_c_coeff + 8) >> 4;
        }
        else {
            rate_control_layer_ptr->k_coeff = (rate_control_layer_ptr->coeff_averaging_weight1*rate_control_layer_ptr->k_coeff + rate_control_layer_ptr->coeff_averaging_weight2*rate_control_layer_ptr->previous_k_coeff + 8) >> 4;
            rate_control_layer_ptr->c_coeff = (rate_control_layer_ptr->coeff_averaging_weight1*rate_control_layer_ptr->c_coeff + rate_control_layer_ptr->coeff_averaging_weight2*rate_control_layer_ptr->previous_c_coeff + 8) >> 4;
        }
        rate_control_layer_ptr->k_coeff = MIN(rate_control_layer_ptr->k_coeff, rate_control_layer_ptr->previous_k_coeff * 4);
        rate_control_layer_ptr->c_coeff = MIN(rate_control_layer_ptr->c_coeff, rate_control_layer_ptr->previous_c_coeff * 4);
        if (parentpicture_control_set_ptr->slice_type != I_SLICE)
            rate_control_layer_ptr->previous_frame_distortion_me = parentpicture_control_set_ptr->sad_me;
        else
            rate_control_layer_ptr->previous_frame_distortion_me = 0;
    }

    if (sequence_control_set_ptr->static_config.look_ahead_distance != 0) {
        if (parentpicture_control_set_ptr->slice_type == I_SLICE) {
            if (parentpicture_control_set_ptr->total_num_bits < parentpicture_control_set_ptr->target_bits_best_pred_qp << 1)
                context_ptr->base_layer_intra_frames_avg_qp = (3 * context_ptr->base_layer_intra_frames_avg_qp + parentpicture_control_set_ptr->picture_qp + 2) >> 2;
            else if (parentpicture_control_set_ptr->total_num_bits > parentpicture_control_set_ptr->target_bits_best_pred_qp << 2)
                context_ptr->base_layer_intra_frames_avg_qp = (3 * context_ptr->base_layer_intra_frames_avg_qp + parentpicture_control_set_ptr->picture_qp + 4 + 2) >> 2;
            else if (parentpicture_control_set_ptr->total_num_bits > parentpicture_control_set_ptr->target_bits_best_pred_qp << 1)
                context_ptr->base_layer_intra_frames_avg_qp = (3 * context_ptr->base_layer_intra_frames_avg_qp + parentpicture_control_set_ptr->picture_qp + 2 + 2) >> 2;
        }
    }

    {
        uint64_t previous_frame_ec_bits = 0;
        EbBool picture_min_qp_allowed = EB_TRUE;
        rate_control_layer_ptr->previous_frame_average_qp = 0;
        rate_control_layer_ptr->previous_frame_average_qp += rate_control_layer_ptr->previous_frame_qp;
        previous_frame_ec_bits += rate_control_layer_ptr->previous_frame_bit_actual;
        if (rate_control_layer_ptr->same_distortion_count == 0 ||
            parentpicture_control_set_ptr->picture_qp != sequence_control_set_ptr->static_config.min_qp_allowed) {
            picture_min_qp_allowed = EB_FALSE;
        }
        if (picture_min_qp_allowed)
            rate_control_layer_ptr->frame_same_distortion_min_qp_count++;
        else
            rate_control_layer_ptr->frame_same_distortion_min_qp_count = 0;

        rate_control_layer_ptr->previous_ec_bits = previous_frame_ec_bits;
        previous_frame_bit_actual = parentpicture_control_set_ptr->total_num_bits;
        if (parentpicture_control_set_ptr->first_frame_in_temporal_layer)
            rate_control_layer_ptr->dif_total_and_ec_bits = (previous_frame_bit_actual - previous_frame_ec_bits);
        else
            rate_control_layer_ptr->dif_total_and_ec_bits = ((previous_frame_bit_actual - previous_frame_ec_bits) + rate_control_layer_ptr->dif_total_and_ec_bits) >> 1;
        // update bitrate of different layers in the interval based on the rate of the I frame
        if (parentpicture_control_set_ptr->picture_number == rate_control_param_ptr->first_poc &&
            (parentpicture_control_set_ptr->slice_type == I_SLICE) &&
            sequence_control_set_ptr->static_config.intra_period_length != -1) {
            uint32_t temporal_layer_idex;
            uint64_t target_bit_rate;
            uint64_t channel_bit_rate;
            uint64_t sum_bits_per_sw = 0;
#if ADAPTIVE_PERCENTAGE
            if (sequence_control_set_ptr->static_config.look_ahead_distance != 0) {
                if (parentpicture_control_set_ptr->tables_updated && parentpicture_control_set_ptr->percentage_updated) {
                    parentpicture_control_set_ptr->bits_per_sw_per_layer[0] =
                        (uint64_t)MAX((int64_t)parentpicture_control_set_ptr->bits_per_sw_per_layer[0] + (int64_t)parentpicture_control_set_ptr->total_num_bits - (int64_t)parentpicture_control_set_ptr->target_bits_best_pred_qp, 1);
                }
            }
#endif

            if (sequence_control_set_ptr->static_config.look_ahead_distance != 0 && sequence_control_set_ptr->intra_period_length != -1) {
                for (temporal_layer_idex = 0; temporal_layer_idex < EB_MAX_TEMPORAL_LAYERS; temporal_layer_idex++)
                    sum_bits_per_sw += parentpicture_control_set_ptr->bits_per_sw_per_layer[temporal_layer_idex];
            }

            for (temporal_layer_idex = 0; temporal_layer_idex < EB_MAX_TEMPORAL_LAYERS; temporal_layer_idex++) {
                rate_control_layer_temp_ptr = rate_control_param_ptr->rate_control_layer_array[temporal_layer_idex];

                target_bit_rate = (uint64_t)((int64_t)parentpicture_control_set_ptr->target_bit_rate -
                    MIN((int64_t)parentpicture_control_set_ptr->target_bit_rate * 3 / 4, (int64_t)(parentpicture_control_set_ptr->total_num_bits*context_ptr->frame_rate / (sequence_control_set_ptr->static_config.intra_period_length + 1)) >> RC_PRECISION))
                    *rate_percentage_layer_array[sequence_control_set_ptr->static_config.hierarchical_levels][temporal_layer_idex] / 100;

#if ADAPTIVE_PERCENTAGE
                if (sequence_control_set_ptr->static_config.look_ahead_distance != 0 && sequence_control_set_ptr->intra_period_length != -1) {
                    target_bit_rate = (uint64_t)((int64_t)parentpicture_control_set_ptr->target_bit_rate -
                        MIN((int64_t)parentpicture_control_set_ptr->target_bit_rate * 3 / 4, (int64_t)(parentpicture_control_set_ptr->total_num_bits*context_ptr->frame_rate / (sequence_control_set_ptr->static_config.intra_period_length + 1)) >> RC_PRECISION))
                        *parentpicture_control_set_ptr->bits_per_sw_per_layer[temporal_layer_idex] / sum_bits_per_sw;
                }
#endif
                // update this based on temporal layers
                if (temporal_layer_idex == 0)
                    channel_bit_rate = (((target_bit_rate << (2 * RC_PRECISION)) / MAX(1, rate_control_layer_temp_ptr->frame_rate - (1 * context_ptr->frame_rate / (sequence_control_set_ptr->static_config.intra_period_length + 1)))) + RC_PRECISION_OFFSET) >> RC_PRECISION;
                else
                    channel_bit_rate = (((target_bit_rate << (2 * RC_PRECISION)) / rate_control_layer_temp_ptr->frame_rate) + RC_PRECISION_OFFSET) >> RC_PRECISION;
                channel_bit_rate = (uint64_t)MAX((int64_t)1, (int64_t)channel_bit_rate);
                rate_control_layer_temp_ptr->ec_bit_constraint = channel_bit_rate;

                slice_num = 1;
                rate_control_layer_temp_ptr->ec_bit_constraint -= SLICE_HEADER_BITS_NUM * slice_num;

                rate_control_layer_temp_ptr->previous_bit_constraint = channel_bit_rate;
                rate_control_layer_temp_ptr->bit_constraint = channel_bit_rate;
                rate_control_layer_temp_ptr->channel_bit_rate = channel_bit_rate;
            }
            if ((int64_t)parentpicture_control_set_ptr->target_bit_rate * 3 / 4 < (int64_t)(parentpicture_control_set_ptr->total_num_bits*context_ptr->frame_rate / (sequence_control_set_ptr->static_config.intra_period_length + 1)) >> RC_PRECISION) {
                rate_control_param_ptr->previous_virtual_buffer_level += (int64_t)((parentpicture_control_set_ptr->total_num_bits*context_ptr->frame_rate / (sequence_control_set_ptr->static_config.intra_period_length + 1)) >> RC_PRECISION) - (int64_t)parentpicture_control_set_ptr->target_bit_rate * 3 / 4;
                context_ptr->extra_bits_gen = 0;
            }
        }

        if (previous_frame_bit_actual) {
            uint64_t bit_changes_rate;
            // Updating virtual buffer level and it can be negative
            context_ptr->extra_bits_gen = 0;
            rate_control_param_ptr->virtual_buffer_level =
                (int64_t)rate_control_param_ptr->previous_virtual_buffer_level +
                (int64_t)previous_frame_bit_actual - (int64_t)context_ptr->high_level_rate_control_ptr->channel_bit_rate_per_frame;
            if (parentpicture_control_set_ptr->hierarchical_levels > 1 && rate_control_layer_ptr->frame_same_distortion_min_qp_count > 10) {
                rate_control_layer_ptr->previous_bit_constraint = (int64_t)rate_control_layer_ptr->channel_bit_rate;
                rate_control_param_ptr->virtual_buffer_level = ((int64_t)context_ptr->virtual_buffer_size >> 1);
            }
            // Updating bit difference
            rate_control_layer_ptr->bit_diff = (int64_t)rate_control_param_ptr->virtual_buffer_level
                //- ((int64_t)context_ptr->virtual_buffer_size>>1);
                - ((int64_t)rate_control_layer_ptr->channel_bit_rate >> 1);

            // Limit the bit difference
            rate_control_layer_ptr->bit_diff = CLIP3(-(int64_t)(rate_control_layer_ptr->channel_bit_rate), (int64_t)(rate_control_layer_ptr->channel_bit_rate >> 1), rate_control_layer_ptr->bit_diff);
            bit_changes_rate = rate_control_layer_ptr->frame_rate;

            // Updating bit Constraint
            rate_control_layer_ptr->bit_constraint = MAX((int64_t)rate_control_layer_ptr->previous_bit_constraint - ((rate_control_layer_ptr->bit_diff << RC_PRECISION) / ((int64_t)bit_changes_rate)), 1);

            // Limiting the bit_constraint
            if (parentpicture_control_set_ptr->temporal_layer_index == 0) {
                rate_control_layer_ptr->bit_constraint = CLIP3(rate_control_layer_ptr->channel_bit_rate >> 2,
                    rate_control_layer_ptr->channel_bit_rate * 200 / 100,
                    rate_control_layer_ptr->bit_constraint);
            }
            else {
                rate_control_layer_ptr->bit_constraint = CLIP3(rate_control_layer_ptr->channel_bit_rate >> 1,
                    rate_control_layer_ptr->channel_bit_rate * 200 / 100,
                    rate_control_layer_ptr->bit_constraint);
            }
            rate_control_layer_ptr->ec_bit_constraint = (uint64_t)MAX((int64_t)rate_control_layer_ptr->bit_constraint - (int64_t)rate_control_layer_ptr->dif_total_and_ec_bits, 1);
            rate_control_param_ptr->previous_virtual_buffer_level = rate_control_param_ptr->virtual_buffer_level;
            rate_control_layer_ptr->previous_bit_constraint = rate_control_layer_ptr->bit_constraint;
        }

        rate_control_param_ptr->processed_frames_number++;
        rate_control_param_ptr->in_use = EB_TRUE;
        // check if all the frames in the interval have arrived
        if (rate_control_param_ptr->processed_frames_number == (rate_control_param_ptr->last_poc - rate_control_param_ptr->first_poc + 1) &&
            sequence_control_set_ptr->intra_period_length != -1) {
            uint32_t temporal_index;
            int64_t extra_bits;
            rate_control_param_ptr->first_poc += PARALLEL_GOP_MAX_NUMBER * (uint32_t)(sequence_control_set_ptr->intra_period_length + 1);
            rate_control_param_ptr->last_poc += PARALLEL_GOP_MAX_NUMBER * (uint32_t)(sequence_control_set_ptr->intra_period_length + 1);
            rate_control_param_ptr->processed_frames_number = 0;
            rate_control_param_ptr->extra_ap_bit_ratio_i = 0;
            rate_control_param_ptr->in_use = EB_FALSE;
            rate_control_param_ptr->was_used = EB_TRUE;
            rate_control_param_ptr->last_gop = EB_FALSE;
            rate_control_param_ptr->first_pic_actual_qp_assigned = EB_FALSE;
            for (temporal_index = 0; temporal_index < EB_MAX_TEMPORAL_LAYERS; temporal_index++) {
                rate_control_param_ptr->rate_control_layer_array[temporal_index]->first_frame = 1;
                rate_control_param_ptr->rate_control_layer_array[temporal_index]->first_non_intra_frame = 1;
                rate_control_param_ptr->rate_control_layer_array[temporal_index]->feedback_arrived = EB_FALSE;
            }
            extra_bits = ((int64_t)context_ptr->virtual_buffer_size >> 1) - (int64_t)rate_control_param_ptr->virtual_buffer_level;

            rate_control_param_ptr->virtual_buffer_level = context_ptr->virtual_buffer_size >> 1;
            context_ptr->extra_bits += extra_bits;
        }
        context_ptr->extra_bits = 0;
    }

#if RC_PRINTS
    ///if (parentpicture_control_set_ptr->temporal_layer_index == 0)
    {
        SVT_LOG("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.0f\t%.0f\t%.0f\t%.0f\t%d\t%d\n",
            (int)parentpicture_control_set_ptr->slice_type,
            (int)parentpicture_control_set_ptr->picture_number,
            (int)parentpicture_control_set_ptr->temporal_layer_index,
            (int)parentpicture_control_set_ptr->picture_qp, (int)parentpicture_control_set_ptr->calculated_qp, (int)parentpicture_control_set_ptr->best_pred_qp,
            (int)previous_frame_bit_actual,
            (int)parentpicture_control_set_ptr->target_bits_best_pred_qp,
            (int)parentpicture_control_set_ptr->target_bits_rc,
            (int)rate_control_layer_ptr->channel_bit_rate,
            (int)rate_control_layer_ptr->bit_constraint,
            (double)rate_control_layer_ptr->c_coeff,
            (double)rate_control_layer_ptr->k_coeff,
            (double)parentpicture_control_set_ptr->sad_me,
            (double)context_ptr->extra_bits_gen,
            (int)rate_control_param_ptr->virtual_buffer_level,
            (int)context_ptr->extra_bits);
    }
#endif
}

void high_level_rc_feed_back_picture(
    PictureParentControlSet *picture_control_set_ptr,
    SequenceControlSet      *sequence_control_set_ptr)
{
    // Queue variables
    HlRateControlHistogramEntry *hl_rate_control_histogram_ptr_temp;
    uint32_t                     queue_entry_index_head_temp;

    //SVT_LOG("\nOut:%d Slidings: ",picture_control_set_ptr->picture_number);
    if (sequence_control_set_ptr->static_config.look_ahead_distance != 0) {
        // Update the coded rate in the histogram queue
        if (picture_control_set_ptr->picture_number >= sequence_control_set_ptr->encode_context_ptr->hl_rate_control_historgram_queue[sequence_control_set_ptr->encode_context_ptr->hl_rate_control_historgram_queue_head_index]->picture_number) {
            queue_entry_index_head_temp = (int32_t)(picture_control_set_ptr->picture_number - sequence_control_set_ptr->encode_context_ptr->hl_rate_control_historgram_queue[sequence_control_set_ptr->encode_context_ptr->hl_rate_control_historgram_queue_head_index]->picture_number);
            queue_entry_index_head_temp += sequence_control_set_ptr->encode_context_ptr->hl_rate_control_historgram_queue_head_index;
            queue_entry_index_head_temp = (queue_entry_index_head_temp > HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1) ?
                queue_entry_index_head_temp - HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH :
                queue_entry_index_head_temp;

            hl_rate_control_histogram_ptr_temp = (sequence_control_set_ptr->encode_context_ptr->hl_rate_control_historgram_queue[queue_entry_index_head_temp]);
            if (hl_rate_control_histogram_ptr_temp->picture_number == picture_control_set_ptr->picture_number &&
                hl_rate_control_histogram_ptr_temp->passed_to_hlrc) {
                eb_block_on_mutex(sequence_control_set_ptr->encode_context_ptr->hl_rate_control_historgram_queue_mutex);
                hl_rate_control_histogram_ptr_temp->total_num_bits_coded = picture_control_set_ptr->total_num_bits;
                hl_rate_control_histogram_ptr_temp->is_coded = EB_TRUE;
                eb_release_mutex(sequence_control_set_ptr->encode_context_ptr->hl_rate_control_historgram_queue_mutex);
            }
        }
    }
}
// rate control QP refinement
void rate_control_refinement(
    PictureControlSet               *picture_control_set_ptr,
    SequenceControlSet              *sequence_control_set_ptr,
    RateControlIntervalParamContext *rate_control_param_ptr,
    RateControlIntervalParamContext *prev_gop_rate_control_param_ptr,
    RateControlIntervalParamContext *next_gop_rate_control_param_ptr) {
    if (picture_control_set_ptr->picture_number == rate_control_param_ptr->first_poc && picture_control_set_ptr->picture_number != 0 && !prev_gop_rate_control_param_ptr->scene_change_in_gop) {
        int16_t deltaApQp = (int16_t)prev_gop_rate_control_param_ptr->first_pic_actual_qp - (int16_t)prev_gop_rate_control_param_ptr->first_pic_pred_qp;
        int64_t extraApBitRatio = (prev_gop_rate_control_param_ptr->first_pic_pred_bits != 0) ?
            (((int64_t)prev_gop_rate_control_param_ptr->first_pic_actual_bits - (int64_t)prev_gop_rate_control_param_ptr->first_pic_pred_bits) * 100) / ((int64_t)prev_gop_rate_control_param_ptr->first_pic_pred_bits) :
            0;
        extraApBitRatio += (int64_t)deltaApQp * 15;
        if (extraApBitRatio > 200)
            picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp + 3;
        else if (extraApBitRatio > 100)
            picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp + 2;
        else if (extraApBitRatio > 50)
            picture_control_set_ptr->picture_qp++;
    }

    if (picture_control_set_ptr->picture_number == rate_control_param_ptr->first_poc && picture_control_set_ptr->picture_number != 0) {
        uint8_t qpIncAllowed = 3;
        uint8_t qpDecAllowed = 4;
        if (picture_control_set_ptr->parent_pcs_ptr->intra_selected_org_qp + 10 <= prev_gop_rate_control_param_ptr->first_pic_actual_qp)
            qpDecAllowed = (uint8_t)(prev_gop_rate_control_param_ptr->first_pic_actual_qp - picture_control_set_ptr->parent_pcs_ptr->intra_selected_org_qp) >> 1;
        if (picture_control_set_ptr->parent_pcs_ptr->intra_selected_org_qp >= prev_gop_rate_control_param_ptr->first_pic_actual_qp + 10)
        {
            qpIncAllowed = (uint8_t)(picture_control_set_ptr->parent_pcs_ptr->intra_selected_org_qp - prev_gop_rate_control_param_ptr->first_pic_actual_qp) * 2 / 3;
            if (prev_gop_rate_control_param_ptr->first_pic_actual_qp <= 15)
                qpIncAllowed += 5;
            else if (prev_gop_rate_control_param_ptr->first_pic_actual_qp <= 20)
                qpIncAllowed += 4;
            else if (prev_gop_rate_control_param_ptr->first_pic_actual_qp <= 25)
                qpIncAllowed += 3;
        }
        else if (prev_gop_rate_control_param_ptr->scene_change_in_gop)
            qpIncAllowed = 5;
        if (picture_control_set_ptr->parent_pcs_ptr->end_of_sequence_region) {
            qpIncAllowed += 2;
            qpDecAllowed += 4;
        }
        picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
            (uint32_t)MAX((int32_t)prev_gop_rate_control_param_ptr->first_pic_actual_qp - (int32_t)qpDecAllowed, 0),
            (uint32_t)prev_gop_rate_control_param_ptr->first_pic_actual_qp + qpIncAllowed,
            picture_control_set_ptr->picture_qp);
    }

    // Scene change
    if (picture_control_set_ptr->slice_type == I_SLICE && picture_control_set_ptr->picture_number != rate_control_param_ptr->first_poc) {
        if (next_gop_rate_control_param_ptr->first_pic_actual_qp_assigned) {
            picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                (uint32_t)MAX((int32_t)next_gop_rate_control_param_ptr->first_pic_actual_qp - (int32_t)1, 0),
                (uint32_t)next_gop_rate_control_param_ptr->first_pic_actual_qp + 8,
                picture_control_set_ptr->picture_qp);
        }
        else {
            if (rate_control_param_ptr->first_pic_actual_qp < 20) {
                picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                    (uint32_t)MAX((int32_t)rate_control_param_ptr->first_pic_actual_qp - (int32_t)4, 0),
                    (uint32_t)rate_control_param_ptr->first_pic_actual_qp + 10,
                    picture_control_set_ptr->picture_qp);
            }
            else {
                picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                    (uint32_t)MAX((int32_t)rate_control_param_ptr->first_pic_actual_qp - (int32_t)4, 0),
                    (uint32_t)rate_control_param_ptr->first_pic_actual_qp + 8,
                    picture_control_set_ptr->picture_qp);
            }
        }
    }

    if (sequence_control_set_ptr->intra_period_length != -1 && picture_control_set_ptr->parent_pcs_ptr->hierarchical_levels < 2 && (int32_t)picture_control_set_ptr->temporal_layer_index == 0 && picture_control_set_ptr->slice_type != I_SLICE) {
        if (next_gop_rate_control_param_ptr->first_pic_actual_qp_assigned || next_gop_rate_control_param_ptr->was_used) {
            picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                (uint32_t)MAX((int32_t)next_gop_rate_control_param_ptr->first_pic_actual_qp - (int32_t)4, 0),
                (uint32_t)picture_control_set_ptr->picture_qp,
                picture_control_set_ptr->picture_qp);
        }
        else {
            picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                (uint32_t)MAX((int32_t)rate_control_param_ptr->first_pic_actual_qp - (int32_t)4, 0),
                (uint32_t)picture_control_set_ptr->picture_qp,
                picture_control_set_ptr->picture_qp);
        }
    }
}
// initialize the rate control parameter at the beginning
void init_rc(
    RateControlContext *context_ptr,
    PictureControlSet  *picture_control_set_ptr,
    SequenceControlSet *sequence_control_set_ptr) {
    context_ptr->high_level_rate_control_ptr->target_bit_rate = sequence_control_set_ptr->static_config.target_bit_rate;
    context_ptr->high_level_rate_control_ptr->frame_rate = sequence_control_set_ptr->frame_rate;
    context_ptr->high_level_rate_control_ptr->channel_bit_rate_per_frame = (uint64_t)MAX((int64_t)1, (int64_t)((context_ptr->high_level_rate_control_ptr->target_bit_rate << RC_PRECISION) / context_ptr->high_level_rate_control_ptr->frame_rate));

    context_ptr->high_level_rate_control_ptr->channel_bit_rate_per_sw = context_ptr->high_level_rate_control_ptr->channel_bit_rate_per_frame * (sequence_control_set_ptr->static_config.look_ahead_distance + 1);
    context_ptr->high_level_rate_control_ptr->bit_constraint_per_sw = context_ptr->high_level_rate_control_ptr->channel_bit_rate_per_sw;

#if RC_UPDATE_TARGET_RATE
    context_ptr->high_level_rate_control_ptr->previous_updated_bit_constraint_per_sw = context_ptr->high_level_rate_control_ptr->channel_bit_rate_per_sw;
#endif

    int32_t total_frame_in_interval = sequence_control_set_ptr->intra_period_length;
    uint32_t gopPeriod = (1 << picture_control_set_ptr->parent_pcs_ptr->hierarchical_levels);
    context_ptr->frame_rate = sequence_control_set_ptr->frame_rate;
    while (total_frame_in_interval >= 0) {
        if (total_frame_in_interval % (gopPeriod) == 0)
            context_ptr->frames_in_interval[0] ++;
        else if (total_frame_in_interval % (gopPeriod >> 1) == 0)
            context_ptr->frames_in_interval[1] ++;
        else if (total_frame_in_interval % (gopPeriod >> 2) == 0)
            context_ptr->frames_in_interval[2] ++;
        else if (total_frame_in_interval % (gopPeriod >> 3) == 0)
            context_ptr->frames_in_interval[3] ++;
        else if (total_frame_in_interval % (gopPeriod >> 4) == 0)
            context_ptr->frames_in_interval[4] ++;
        else if (total_frame_in_interval % (gopPeriod >> 5) == 0)
            context_ptr->frames_in_interval[5] ++;
        total_frame_in_interval--;
    }
    if (sequence_control_set_ptr->static_config.rate_control_mode == 2) { // VBR
        context_ptr->virtual_buffer_size = (((uint64_t)sequence_control_set_ptr->static_config.target_bit_rate * 3) << RC_PRECISION) / (context_ptr->frame_rate);
        context_ptr->rate_average_periodin_frames = (uint64_t)sequence_control_set_ptr->static_config.intra_period_length + 1;
        context_ptr->virtual_buffer_level_initial_value = context_ptr->virtual_buffer_size >> 1;
        context_ptr->virtual_buffer_level = context_ptr->virtual_buffer_size >> 1;
        context_ptr->previous_virtual_buffer_level = context_ptr->virtual_buffer_size >> 1;
        context_ptr->vb_fill_threshold1 = (context_ptr->virtual_buffer_size * 6) >> 3;
        context_ptr->vb_fill_threshold2 = (context_ptr->virtual_buffer_size << 3) >> 3;
        context_ptr->base_layer_frames_avg_qp = sequence_control_set_ptr->qp;
        context_ptr->base_layer_intra_frames_avg_qp = sequence_control_set_ptr->qp;
    }
    else if (sequence_control_set_ptr->static_config.rate_control_mode == 3) {
        context_ptr->virtual_buffer_size = ((uint64_t)sequence_control_set_ptr->static_config.target_bit_rate);// vbv_buf_size);
        context_ptr->rate_average_periodin_frames = (uint64_t)sequence_control_set_ptr->static_config.intra_period_length + 1;
        context_ptr->virtual_buffer_level_initial_value = context_ptr->virtual_buffer_size >> 1;
        context_ptr->virtual_buffer_level = context_ptr->virtual_buffer_size >> 1;
        context_ptr->previous_virtual_buffer_level = context_ptr->virtual_buffer_size >> 1;
        context_ptr->vb_fill_threshold1 = context_ptr->virtual_buffer_level_initial_value + (context_ptr->virtual_buffer_size / 4);
        context_ptr->vb_fill_threshold2 = context_ptr->virtual_buffer_level_initial_value + (context_ptr->virtual_buffer_size / 3);
        context_ptr->base_layer_frames_avg_qp = sequence_control_set_ptr->qp;
        context_ptr->base_layer_intra_frames_avg_qp = sequence_control_set_ptr->qp;
    }

    for (uint32_t base_qp = 0; base_qp < MAX_REF_QP_NUM; base_qp++) {
        if (base_qp < 64) {
            context_ptr->qp_scaling_map_I_SLICE[base_qp] = qp_scaling_calc(
                sequence_control_set_ptr,
                I_SLICE,
                0,
                base_qp);
        }
        else
            context_ptr->qp_scaling_map_I_SLICE[base_qp] = (uint32_t)CLIP3(0, 63, (int)base_qp - (63 - (int)context_ptr->qp_scaling_map_I_SLICE[63]));
        for (uint32_t temporal_layer_index = 0; temporal_layer_index < sequence_control_set_ptr->static_config.hierarchical_levels+1; temporal_layer_index++) {
            if (base_qp < 64) {
                context_ptr->qp_scaling_map[temporal_layer_index][base_qp] = qp_scaling_calc(
                    sequence_control_set_ptr,
                    0,
                    temporal_layer_index,
                    base_qp);
            }
            else
                context_ptr->qp_scaling_map[temporal_layer_index][base_qp] = (uint32_t)CLIP3(0, 63, (int)base_qp - (63 - (int)context_ptr->qp_scaling_map[temporal_layer_index][63]));
        }
    }
}

#define MAX_Q_INDEX 255
#define MIN_Q_INDEX 0

extern int16_t av1_ac_quant_Q3(int32_t qindex, int32_t delta, AomBitDepth bit_depth);
// These functions use formulaic calculations to make playing with the
// quantizer tables easier. If necessary they can be replaced by lookup
// tables if and when things settle down in the experimental bitstream

double av1_convert_qindex_to_q(int32_t qindex, AomBitDepth bit_depth) {
    // Convert the index to a real Q value (scaled down to match old Q values)
    switch (bit_depth) {
    case AOM_BITS_8: return av1_ac_quant_Q3(qindex, 0, bit_depth) / 4.0;
    case AOM_BITS_10: return av1_ac_quant_Q3(qindex, 0, bit_depth) / 16.0;
    case AOM_BITS_12: return av1_ac_quant_Q3(qindex, 0, bit_depth) / 64.0;
    default:
        assert(0 && "bit_depth should be AOM_BITS_8, AOM_BITS_10 or AOM_BITS_12");
        return -1.0;
    }
}
int32_t av1_compute_qdelta(double qstart, double qtarget,
    AomBitDepth bit_depth) {
    int32_t start_index = MAX_Q_INDEX;
    int32_t target_index = MAX_Q_INDEX;
    int32_t i;

    // Convert the average q value to an index.
    for (i = MIN_Q_INDEX; i < MAX_Q_INDEX; ++i) {
        start_index = i;
        if (av1_convert_qindex_to_q(i, bit_depth) >= qstart) break;
    }

    // Convert the q target to an index
    for (i = MIN_Q_INDEX; i < MAX_Q_INDEX; ++i) {
        target_index = i;
        if (av1_convert_qindex_to_q(i, bit_depth) >= qtarget) break;
    }

    return target_index - start_index;
}
// calculate the QP based on the QP scaling
uint32_t qp_scaling_calc(
    SequenceControlSet *sequence_control_set_ptr,
    EB_SLICE            slice_type,
    uint32_t            temporal_layer_index,
    uint32_t            base_qp)
{
    // AMIR to fix
    uint32_t    scaled_qp = 0;
    int         base_qindex;

    const  double delta_rate_new[2][6] =
    { { 0.40, 0.7, 0.85, 1.0, 1.0, 1.0 },
    { 0.35, 0.6, 0.8,  0.9, 1.0, 1.0 } };

    int qindex = quantizer_to_qindex[base_qp];
    const double q = av1_convert_qindex_to_q(qindex, (AomBitDepth)sequence_control_set_ptr->static_config.encoder_bit_depth);
    int delta_qindex;

    if (slice_type == I_SLICE) {
        delta_qindex = av1_compute_qdelta(
            q,
            q* 0.25,
            (AomBitDepth)sequence_control_set_ptr->static_config.encoder_bit_depth);
    }
    else {
        delta_qindex = av1_compute_qdelta(
            q,
            q* delta_rate_new[sequence_control_set_ptr->static_config.hierarchical_levels == 4][temporal_layer_index], // RC does not support 5L
            //q* delta_rate_new[0][temporal_layer_index], // RC does not support 5L
            (AomBitDepth)sequence_control_set_ptr->static_config.encoder_bit_depth);
    }

    base_qindex = MAX(qindex + delta_qindex, MIN_Q_INDEX);
    scaled_qp = (uint32_t)(base_qindex) >> 2;

    return scaled_qp;
}
typedef struct {
    // Rate targetting variables
    int base_frame_target;  // A baseline frame target before adjustment
                            // for previous under or over shoot.
    int this_frame_target;  // Actual frame target after rc adjustment.
    int projected_frame_size;
    int sb64_target_rate;
    int last_q[FRAME_TYPES];  // Separate values for Intra/Inter
    int last_boosted_qindex;  // Last boosted GF/KF/ARF q
    int last_kf_qindex;       // Q index of the last key frame coded.

    int gfu_boost;
    int kf_boost;

    // double rate_correction_factors[RATE_FACTOR_LEVELS];

    int frames_since_golden;
    int frames_till_gf_update_due;
    int min_gf_interval;
    int max_gf_interval;
    int static_scene_max_gf_interval;
    int baseline_gf_interval;
    int constrained_gf_group;
    int frames_to_key;
    int frames_since_key;
    int this_key_frame_forced;
    int next_key_frame_forced;
    int source_alt_ref_pending;
    int source_alt_ref_active;
    int is_src_frame_alt_ref;
    int sframe_due;

    // Length of the bi-predictive frame group interval
    int bipred_group_interval;

    // NOTE: Different types of frames may have different bits allocated
    //       accordingly, aiming to achieve the overall optimal RD performance.
    int is_bwd_ref_frame;
    int is_last_bipred_frame;
    int is_bipred_frame;
    int is_src_frame_ext_arf;

    int avg_frame_bandwidth;  // Average frame size target for clip
    int min_frame_bandwidth;  // Minimum allocation used for any frame
    int max_frame_bandwidth;  // Maximum burst rate allowed for a frame.

    int ni_av_qi;
    int ni_tot_qi;
    int ni_frames;
    int avg_frame_qindex[FRAME_TYPES];
    double tot_q;
    double avg_q;

    int64_t buffer_level;
    int64_t bits_off_target;
    int64_t vbr_bits_off_target;
    int64_t vbr_bits_off_target_fast;

    int decimation_factor;
    int decimation_count;

    int rolling_target_bits;
    int rolling_actual_bits;

    int long_rolling_target_bits;
    int long_rolling_actual_bits;

    int rate_error_estimate;

    int64_t total_actual_bits;
    int64_t total_target_bits;
    int64_t total_target_vs_actual;

    int worst_quality;
    int best_quality;

    int64_t starting_buffer_level;
    int64_t optimal_buffer_level;
    int64_t maximum_buffer_size;

    // rate control history for last frame(1) and the frame before(2).
    // -1: undershot
    //  1: overshoot
    //  0: not initialized.
    int rc_1_frame;
    int rc_2_frame;
    int q_1_frame;
    int q_2_frame;

    // Auto frame-scaling variables.
 //   int rf_level_maxq[RATE_FACTOR_LEVELS];
    float_t arf_boost_factor;
    // Q index used for ALT frame
    int arf_q;
} RATE_CONTROL;
#define STATIC_MOTION_THRESH 95

enum {
    INTER_NORMAL = 0,
    INTER_LOW = 1,
    INTER_HIGH = 2,
    GF_ARF_LOW = 3,
    GF_ARF_STD = 4,
    KF_STD = 5,
    RATE_FACTOR_LEVELS = 6
} RATE_FACTOR_LEVEL;

enum {
    KF_UPDATE = 0,
    LF_UPDATE = 1,
    GF_UPDATE = 2,
    ARF_UPDATE = 3,
    OVERLAY_UPDATE = 4,
    BRF_UPDATE = 5,            // Backward Reference Frame
    LAST_BIPRED_UPDATE = 6,    // Last Bi-predictive Frame
    BIPRED_UPDATE = 7,         // Bi-predictive Frame, but not the last one
    INTNL_OVERLAY_UPDATE = 8,  // Internal Overlay Frame
    INTNL_ARF_UPDATE = 9,      // Internal Altref Frame (candidate for ALTREF2)
    FRAME_UPDATE_TYPES = 10
} FRAME_UPDATE_TYPE;

// that are not marked as coded with 0,0 motion in the first pass.
#define FAST_MOVING_KF_GROUP_THRESH 5
#if QPS_TUNING
#define MEDIUM_MOVING_KF_GROUP_THRESH  30
#define STATIC_KF_GROUP_THRESH         80
#define MAX_QPS_COMP_I                100
#define MAX_QPS_COMP_NONI             300
#define HIGH_QPS_COMP_THRESHOLD        80
#define LOW_QPS_COMP_THRESHOLD         40
#define HIGH_FILTERED_THRESHOLD     (4<<8) // 8 bit precision
#define LOW_FILTERED_THRESHOLD      (1<<8) // 8 bit precision
#else
#define STATIC_KF_GROUP_THRESH 99
#define MAX_QPS_COMP_I        60
#define MAX_QPS_COMP_NONI    200
#endif
#define QPS_SW_THRESH          8

#define ASSIGN_MINQ_TABLE(bit_depth, name)                   \
  do {                                                       \
    switch (bit_depth) {                                     \
      case AOM_BITS_8: name = name##_8; break;               \
      case AOM_BITS_10: name = name##_10; break;             \
      case AOM_BITS_12: name = name##_12; break;             \
      default:                                               \
        assert(0 &&                                          \
               "bit_depth should be AOM_BITS_8, AOM_BITS_10" \
               " or AOM_BITS_12");                           \
        name = NULL;                                         \
    }                                                        \
  } while (0)

// Tables relating active max Q to active min Q
static int kf_low_motion_minq_8[QINDEX_RANGE];
static int kf_high_motion_minq_8[QINDEX_RANGE];
static int arfgf_low_motion_minq_8[QINDEX_RANGE];
static int arfgf_high_motion_minq_8[QINDEX_RANGE];
static int inter_minq_8[QINDEX_RANGE];
static int rtc_minq_8[QINDEX_RANGE];

static int kf_low_motion_minq_10[QINDEX_RANGE];
static int kf_high_motion_minq_10[QINDEX_RANGE];
static int arfgf_low_motion_minq_10[QINDEX_RANGE];
static int arfgf_high_motion_minq_10[QINDEX_RANGE];
static int inter_minq_10[QINDEX_RANGE];
static int rtc_minq_10[QINDEX_RANGE];
static int kf_low_motion_minq_12[QINDEX_RANGE];
static int kf_high_motion_minq_12[QINDEX_RANGE];
static int arfgf_low_motion_minq_12[QINDEX_RANGE];
static int arfgf_high_motion_minq_12[QINDEX_RANGE];
static int inter_minq_12[QINDEX_RANGE];
static int rtc_minq_12[QINDEX_RANGE];

static int gf_high = 2000;
static int gf_low = 400;
static int kf_high = 5000;
static int kf_low = 400;

// Functions to compute the active minq lookup table entries based on a
// formulaic approach to facilitate easier adjustment of the Q tables.
// The formulae were derived from computing a 3rd order polynomial best
// fit to the original data (after plotting real maxq vs minq (not q index))
static int get_minq_index(double maxq, double x3, double x2, double x1,
    AomBitDepth bit_depth) {
    int i;
    const double minqtarget = AOMMIN(((x3 * maxq + x2) * maxq + x1) * maxq, maxq);

    // Special case handling to deal with the step from q2.0
    // down to lossless mode represented by q 1.0.
    if (minqtarget <= 2.0) return 0;

    for (i = 0; i < QINDEX_RANGE; i++)
        if (minqtarget <= av1_convert_qindex_to_q(i, bit_depth)) return i;
    return QINDEX_RANGE - 1;
}

static void init_minq_luts(int *kf_low_m, int *kf_high_m, int *arfgf_low,
    int *arfgf_high, int *inter, int *rtc,
    AomBitDepth bit_depth) {
    int i;
    for (i = 0; i < QINDEX_RANGE; i++) {
        const double maxq = av1_convert_qindex_to_q(i, bit_depth);
        kf_low_m[i] = get_minq_index(maxq, 0.000001, -0.0004, 0.150, bit_depth);
        kf_high_m[i] = get_minq_index(maxq, 0.0000021, -0.00125, 0.45, bit_depth);
        arfgf_low[i] = get_minq_index(maxq, 0.0000015, -0.0009, 0.30, bit_depth);
        arfgf_high[i] = get_minq_index(maxq, 0.0000021, -0.00125, 0.55, bit_depth);
        inter[i] = get_minq_index(maxq, 0.00000271, -0.00113, 0.90, bit_depth);
        rtc[i] = get_minq_index(maxq, 0.00000271, -0.00113, 0.70, bit_depth);
    }
}

void av1_rc_init_minq_luts(void) {
    init_minq_luts(kf_low_motion_minq_8, kf_high_motion_minq_8,
        arfgf_low_motion_minq_8, arfgf_high_motion_minq_8,
        inter_minq_8, rtc_minq_8, AOM_BITS_8);
    init_minq_luts(kf_low_motion_minq_10, kf_high_motion_minq_10,
        arfgf_low_motion_minq_10, arfgf_high_motion_minq_10,
        inter_minq_10, rtc_minq_10, AOM_BITS_10);
    init_minq_luts(kf_low_motion_minq_12, kf_high_motion_minq_12,
        arfgf_low_motion_minq_12, arfgf_high_motion_minq_12,
        inter_minq_12, rtc_minq_12, AOM_BITS_12);
}

static int get_active_quality(int q, int gfu_boost, int low, int high,
    int *low_motion_minq, int *high_motion_minq) {
    if (gfu_boost > high)
        return low_motion_minq[q];
    else if (gfu_boost < low)
        return high_motion_minq[q];
    else {
        const int gap = high - low;
        const int offset = high - gfu_boost;
        const int qdiff = high_motion_minq[q] - low_motion_minq[q];
        const int adjustment = ((offset * qdiff) + (gap >> 1)) / gap;
        return low_motion_minq[q] + adjustment;
    }
}

static int get_kf_active_quality(const RATE_CONTROL *const rc, int q,
    AomBitDepth bit_depth) {
    int *kf_low_motion_minq;
    int *kf_high_motion_minq;
    ASSIGN_MINQ_TABLE(bit_depth, kf_low_motion_minq);
    ASSIGN_MINQ_TABLE(bit_depth, kf_high_motion_minq);
    return get_active_quality(q, rc->kf_boost, kf_low, kf_high,
        kf_low_motion_minq, kf_high_motion_minq);
}

static int get_gf_active_quality(const RATE_CONTROL *const rc, int q,
    AomBitDepth bit_depth) {
    int *arfgf_low_motion_minq;
    int *arfgf_high_motion_minq;
    ASSIGN_MINQ_TABLE(bit_depth, arfgf_low_motion_minq);
    ASSIGN_MINQ_TABLE(bit_depth, arfgf_high_motion_minq);
    return get_active_quality(q, rc->gfu_boost, gf_low, gf_high,
        arfgf_low_motion_minq, arfgf_high_motion_minq);
}

static int get_gf_high_motion_quality(int q, AomBitDepth bit_depth) {
    int *arfgf_high_motion_minq;
    ASSIGN_MINQ_TABLE(bit_depth, arfgf_high_motion_minq);
    return arfgf_high_motion_minq[q];
}

static int adaptive_qindex_calc(
    PictureControlSet         *picture_control_set_ptr,
    RATE_CONTROL                *rc,
    int                          qindex) {
    SequenceControlSet        *sequence_control_set_ptr = picture_control_set_ptr->parent_pcs_ptr->sequence_control_set_ptr;
    const Av1Common  *const cm = picture_control_set_ptr->parent_pcs_ptr->av1_cm;

    const int cq_level = qindex;
    int active_best_quality = 0;
    int active_worst_quality = qindex;
    rc->arf_q = 0;
    int q;
    int is_src_frame_alt_ref, refresh_golden_frame, refresh_alt_ref_frame, is_intrl_arf_boost, rf_level, update_type;
    is_src_frame_alt_ref = 0;
    refresh_golden_frame = frame_is_intra_only(picture_control_set_ptr->parent_pcs_ptr) ? 1 : 0;
    refresh_alt_ref_frame = (picture_control_set_ptr->parent_pcs_ptr->temporal_layer_index == 0) ? 1 : 0;
    is_intrl_arf_boost = (picture_control_set_ptr->parent_pcs_ptr->temporal_layer_index > 0 && picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag) ? 1 : 0;
    rf_level = (frame_is_intra_only(picture_control_set_ptr->parent_pcs_ptr)) ? KF_STD :
        (picture_control_set_ptr->parent_pcs_ptr->temporal_layer_index == 0) ? GF_ARF_STD :
        picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag ? GF_ARF_LOW : INTER_NORMAL;

    update_type = (frame_is_intra_only(picture_control_set_ptr->parent_pcs_ptr)) ? KF_UPDATE :
        (picture_control_set_ptr->parent_pcs_ptr->temporal_layer_index == 0) ? ARF_UPDATE :
        picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag ? INTNL_ARF_UPDATE : LF_UPDATE;
    const int bit_depth = sequence_control_set_ptr->static_config.encoder_bit_depth;
    // Since many frames can be processed at the same time, storing/using arf_q in rc param is not sufficient and will create a run to run.
    // So, for each frame, arf_q is updated based on the qp of its references.
    if (picture_control_set_ptr->ref_slice_type_array[0][0] != I_SLICE)
        rc->arf_q = MAX(rc->arf_q, ((picture_control_set_ptr->ref_pic_qp_array[0][0] << 2) + 2));
    if ((picture_control_set_ptr->slice_type == B_SLICE) && (picture_control_set_ptr->ref_slice_type_array[1][0] != I_SLICE))
        rc->arf_q = MAX(rc->arf_q, ((picture_control_set_ptr->ref_pic_qp_array[1][0] << 2) + 2));
    if (frame_is_intra_only(picture_control_set_ptr->parent_pcs_ptr)) {
        // Not forced keyframe.
        double q_adj_factor = 1.0;
        double q_val;

        rc->worst_quality = MAXQ;
        rc->best_quality = MINQ;
        int max_qp_scaling_avg_comp_I = sequence_control_set_ptr->input_resolution < 2 ? (MAX_QPS_COMP_I >> 1) : MAX_QPS_COMP_I;
#if QPS_TUNING
        // Update the complexity for very fast moving content.
        if (picture_control_set_ptr->parent_pcs_ptr->kf_zeromotion_pct <= FAST_MOVING_KF_GROUP_THRESH)
            picture_control_set_ptr->parent_pcs_ptr->qp_scaling_average_complexity = max_qp_scaling_avg_comp_I;

        // For the low filtered ALT_REF pictures (next ALT_REF) where complexity is low and picture is static, decrease the complexity/QP of the I_SLICE.
        // The improved area will be propagated to future frames
        if (picture_control_set_ptr->parent_pcs_ptr->qp_scaling_average_complexity <= LOW_QPS_COMP_THRESHOLD &&
            picture_control_set_ptr->parent_pcs_ptr->filtered_sse < LOW_FILTERED_THRESHOLD && picture_control_set_ptr->parent_pcs_ptr->filtered_sse_uv < LOW_FILTERED_THRESHOLD &&
            picture_control_set_ptr->parent_pcs_ptr->kf_zeromotion_pct > STATIC_KF_GROUP_THRESH)
            picture_control_set_ptr->parent_pcs_ptr->qp_scaling_average_complexity >>= 1;

        // For the highly filtered ALT_REF pictures (next ALT_REF), increase the complexity/QP of the I_SLICE to save on rate
        if (picture_control_set_ptr->parent_pcs_ptr->filtered_sse + picture_control_set_ptr->parent_pcs_ptr->filtered_sse_uv >= HIGH_FILTERED_THRESHOLD)
            picture_control_set_ptr->parent_pcs_ptr->qp_scaling_average_complexity = max_qp_scaling_avg_comp_I;
#else
        // Update the complexity for very fast moving content
        if (picture_control_set_ptr->parent_pcs_ptr->kf_zeromotion_pct <= FAST_MOVING_KF_GROUP_THRESH)
            picture_control_set_ptr->parent_pcs_ptr->qp_scaling_average_complexity <<= 1;
#endif
        picture_control_set_ptr->parent_pcs_ptr->qp_scaling_average_complexity = MIN(max_qp_scaling_avg_comp_I, picture_control_set_ptr->parent_pcs_ptr->qp_scaling_average_complexity);

        // cross multiplication to derive kf_boost from non_moving_average_score; kf_boost range is [kf_low,kf_high], and non_moving_average_score range [0,max_qp_scaling_avg_comp_I]
        rc->kf_boost = (((max_qp_scaling_avg_comp_I - (picture_control_set_ptr->parent_pcs_ptr->qp_scaling_average_complexity))  * (kf_high - kf_low)) / max_qp_scaling_avg_comp_I) + kf_low;
        // Baseline value derived from cpi->active_worst_quality and kf boost.
        active_best_quality =
            get_kf_active_quality(rc, active_worst_quality, bit_depth);
        // Allow somewhat lower kf minq with small image formats.
        if ((cm->frm_size.frame_width * cm->frm_size.frame_height) <= (352 * 288))
            q_adj_factor -= 0.25;
        // Make a further adjustment based on the kf zero motion measure.
        q_adj_factor += 0.05 - (0.001 * (double)picture_control_set_ptr->parent_pcs_ptr->kf_zeromotion_pct/*(double)cpi->twopass.kf_zeromotion_pct*/);

        // Convert the adjustment factor to a qindex delta
        // on active_best_quality.
        q_val = av1_convert_qindex_to_q(active_best_quality, bit_depth);
        active_best_quality +=
            av1_compute_qdelta(q_val, q_val * q_adj_factor, bit_depth);
    }
    else if (!is_src_frame_alt_ref &&
        (refresh_golden_frame || is_intrl_arf_boost ||
            refresh_alt_ref_frame)) {
#if QPS_TUNING
        // Clip the complexity of highly complex pictures to maximum.
        if (picture_control_set_ptr->parent_pcs_ptr->qp_scaling_average_complexity > HIGH_QPS_COMP_THRESHOLD)
            picture_control_set_ptr->parent_pcs_ptr->qp_scaling_average_complexity = MAX_QPS_COMP_NONI;
#endif
        rc->gfu_boost = (((MAX_QPS_COMP_NONI - (picture_control_set_ptr->parent_pcs_ptr->qp_scaling_average_complexity))  * (gf_high - gf_low)) / MAX_QPS_COMP_NONI) + gf_low;
#if QPS_TUNING
        // For the highly filtered ALT_REF pictures or where complexity is medium or picture is medium moving, add a boost to decrease the QP of the ALT_REF.
        // The improved area will be propagated to future frames
        rc->arf_boost_factor = (picture_control_set_ptr->parent_pcs_ptr->qp_scaling_average_complexity > LOW_QPS_COMP_THRESHOLD || picture_control_set_ptr->parent_pcs_ptr->kf_zeromotion_pct < MEDIUM_MOVING_KF_GROUP_THRESH || picture_control_set_ptr->parent_pcs_ptr->filtered_sse >= HIGH_FILTERED_THRESHOLD) ?
            (float_t)1.3 : 1;
#else
        rc->arf_boost_factor = 1;
#endif
        q = active_worst_quality;

        // non ref frame or repeated frames with re-encode
        if (!refresh_alt_ref_frame && !is_intrl_arf_boost)
            active_best_quality = cq_level;
        else {
            // base layer
            if (update_type == ARF_UPDATE) {
                active_best_quality = get_gf_active_quality(rc, q, bit_depth);
                //*arf_q = active_best_quality;
                rc->arf_q = active_best_quality;
                const int min_boost = get_gf_high_motion_quality(q, bit_depth);
                const int boost = min_boost - active_best_quality;

                active_best_quality = min_boost - (int)(boost * rc->arf_boost_factor);
            }
            else
                active_best_quality = rc->arf_q;
            // active_best_quality is updated with the q index of the reference
            if (rf_level == GF_ARF_LOW)
                active_best_quality = (active_best_quality + cq_level + 1) / 2;
        }
    }
    else
        active_best_quality = cq_level;
    q = active_best_quality;
    clamp(q, active_best_quality, active_worst_quality);

    return q;
}

void* rate_control_kernel(void *input_ptr)
{
    // Context
    RateControlContext                *context_ptr = (RateControlContext  *)input_ptr;

    RateControlIntervalParamContext   *rate_control_param_ptr;

    RateControlIntervalParamContext   *prev_gop_rate_control_param_ptr;
    RateControlIntervalParamContext   *next_gop_rate_control_param_ptr;

    PictureControlSet                 *picture_control_set_ptr;
    PictureParentControlSet           *parentpicture_control_set_ptr;

    // Config
    SequenceControlSet                *sequence_control_set_ptr;

    // Input
    EbObjectWrapper                   *rate_control_tasks_wrapper_ptr;
    RateControlTasks                  *rate_control_tasks_ptr;

    // Output
    EbObjectWrapper                   *rate_control_results_wrapper_ptr;
    RateControlResults                *rate_control_results_ptr;

    RateControlLayerContext           *rate_control_layer_ptr;

    uint64_t                           total_number_of_fb_frames = 0;

    RateControlTaskTypes               task_type;
    EbRateControlModel          *rc_model_ptr;
    RATE_CONTROL                 rc;

    rc_model_ptr = context_ptr->rc_model_ptr;

    for (;;) {
        // Get RateControl Task
        eb_get_full_object(
            context_ptr->rate_control_input_tasks_fifo_ptr,
            &rate_control_tasks_wrapper_ptr);

        rate_control_tasks_ptr = (RateControlTasks*)rate_control_tasks_wrapper_ptr->object_ptr;
        task_type = rate_control_tasks_ptr->task_type;

        // Modify these for different temporal layers later
        switch (task_type) {
        case RC_PICTURE_MANAGER_RESULT:

            picture_control_set_ptr = (PictureControlSet  *)rate_control_tasks_ptr->picture_control_set_wrapper_ptr->object_ptr;
            sequence_control_set_ptr = (SequenceControlSet *)picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr;
            FrameHeader *frm_hdr = &picture_control_set_ptr->parent_pcs_ptr->frm_hdr;

            if (picture_control_set_ptr->picture_number == 0) {
                rate_control_model_init(rc_model_ptr, sequence_control_set_ptr);

                av1_rc_init_minq_luts();
                //init rate control parameters
                init_rc(
                    context_ptr,
                    picture_control_set_ptr,
                    sequence_control_set_ptr);
            }
            if (sequence_control_set_ptr->static_config.rate_control_mode)
            {
                picture_control_set_ptr->parent_pcs_ptr->intra_selected_org_qp = 0;
                // High level RC
                if (sequence_control_set_ptr->static_config.rate_control_mode == 2)
                    high_level_rc_input_picture_vbr(
                        picture_control_set_ptr->parent_pcs_ptr,
                        sequence_control_set_ptr,
                        sequence_control_set_ptr->encode_context_ptr,
                        context_ptr,
                        context_ptr->high_level_rate_control_ptr);
                else if (sequence_control_set_ptr->static_config.rate_control_mode == 3)
                    high_level_rc_input_picture_cvbr(
                        picture_control_set_ptr->parent_pcs_ptr,
                        sequence_control_set_ptr,
                        sequence_control_set_ptr->encode_context_ptr,
                        context_ptr,
                        context_ptr->high_level_rate_control_ptr);
            }

            // Frame level RC. Find the ParamPtr for the current GOP
            if (sequence_control_set_ptr->intra_period_length == -1 || sequence_control_set_ptr->static_config.rate_control_mode == 0) {
                rate_control_param_ptr = context_ptr->rate_control_param_queue[0];
                prev_gop_rate_control_param_ptr = context_ptr->rate_control_param_queue[0];
                next_gop_rate_control_param_ptr = context_ptr->rate_control_param_queue[0];
            }
            else {
                uint32_t interval_index_temp = 0;
                EbBool intervalFound = EB_FALSE;
                while ((interval_index_temp < PARALLEL_GOP_MAX_NUMBER) && !intervalFound) {
                    if (picture_control_set_ptr->picture_number >= context_ptr->rate_control_param_queue[interval_index_temp]->first_poc &&
                        picture_control_set_ptr->picture_number <= context_ptr->rate_control_param_queue[interval_index_temp]->last_poc) {
                        intervalFound = EB_TRUE;
                    }
                    else
                        interval_index_temp++;
                }
                CHECK_REPORT_ERROR(
                    interval_index_temp != PARALLEL_GOP_MAX_NUMBER,
                    sequence_control_set_ptr->encode_context_ptr->app_callback_ptr,
                    EB_ENC_RC_ERROR2);

                rate_control_param_ptr = context_ptr->rate_control_param_queue[interval_index_temp];

                prev_gop_rate_control_param_ptr = (interval_index_temp == 0) ?
                    context_ptr->rate_control_param_queue[PARALLEL_GOP_MAX_NUMBER - 1] :
                    context_ptr->rate_control_param_queue[interval_index_temp - 1];
                next_gop_rate_control_param_ptr = (interval_index_temp == PARALLEL_GOP_MAX_NUMBER - 1) ?
                    context_ptr->rate_control_param_queue[0] :
                    context_ptr->rate_control_param_queue[interval_index_temp + 1];
            }

            rate_control_layer_ptr = rate_control_param_ptr->rate_control_layer_array[picture_control_set_ptr->temporal_layer_index];

            if (sequence_control_set_ptr->static_config.rate_control_mode == 0) {
                // if RC mode is 0,  fixed QP is used
                // QP scaling based on POC number for Flat IPPP structure
                frm_hdr->quantization_params.base_q_idx = quantizer_to_qindex[picture_control_set_ptr->picture_qp];

                if (sequence_control_set_ptr->static_config.enable_qp_scaling_flag && picture_control_set_ptr->parent_pcs_ptr->qp_on_the_fly == EB_FALSE) {
                    const int32_t qindex = quantizer_to_qindex[(uint8_t)sequence_control_set_ptr->qp];
                    const double q_val = av1_convert_qindex_to_q(qindex, (AomBitDepth)sequence_control_set_ptr->static_config.encoder_bit_depth);
                    // if there are need enough pictures in the LAD/SlidingWindow, the adaptive QP scaling is not used
                    if (picture_control_set_ptr->parent_pcs_ptr->frames_in_sw >= QPS_SW_THRESH) {
                        int32_t new_qindex = adaptive_qindex_calc(
                            picture_control_set_ptr,
                            &rc,
                            qindex);

                        frm_hdr->quantization_params.base_q_idx =
                            (uint8_t)CLIP3(
                            (int32_t)quantizer_to_qindex[sequence_control_set_ptr->static_config.min_qp_allowed],
                                (int32_t)quantizer_to_qindex[sequence_control_set_ptr->static_config.max_qp_allowed],
                                (int32_t)(new_qindex));
                    }
                    else if (picture_control_set_ptr->slice_type == I_SLICE) {
                        const int32_t delta_qindex = av1_compute_qdelta(
                            q_val,
                            q_val * 0.25,
                            (AomBitDepth)sequence_control_set_ptr->static_config.encoder_bit_depth);
                        frm_hdr->quantization_params.base_q_idx =
                            (uint8_t)CLIP3(
                            (int32_t)quantizer_to_qindex[sequence_control_set_ptr->static_config.min_qp_allowed],
                                (int32_t)quantizer_to_qindex[sequence_control_set_ptr->static_config.max_qp_allowed],
                                (int32_t)(qindex + delta_qindex));
                    }
                    else {
                        const  double delta_rate_new[2][6] =
                        { { 0.40, 0.7, 0.85, 1.0, 1.0, 1.0 },
                        { 0.35, 0.6, 0.8,  0.9, 1.0, 1.0 } };

                        const int32_t delta_qindex = av1_compute_qdelta(
                            q_val,
                            q_val * delta_rate_new[picture_control_set_ptr->parent_pcs_ptr->hierarchical_levels == 4][picture_control_set_ptr->parent_pcs_ptr->temporal_layer_index],
                            (AomBitDepth)sequence_control_set_ptr->static_config.encoder_bit_depth);

                        frm_hdr->quantization_params.base_q_idx =
                            (uint8_t)CLIP3(
                            (int32_t)quantizer_to_qindex[sequence_control_set_ptr->static_config.min_qp_allowed],
                                (int32_t)quantizer_to_qindex[sequence_control_set_ptr->static_config.max_qp_allowed],
                                (int32_t)(qindex + delta_qindex));
                    }
                    picture_control_set_ptr->picture_qp = (uint8_t)CLIP3((int32_t)sequence_control_set_ptr->static_config.min_qp_allowed, (int32_t)sequence_control_set_ptr->static_config.max_qp_allowed, frm_hdr->quantization_params.base_q_idx >> 2);
                }

                else if (picture_control_set_ptr->parent_pcs_ptr->qp_on_the_fly == EB_TRUE) {
                    picture_control_set_ptr->picture_qp = (uint8_t)CLIP3((int32_t)sequence_control_set_ptr->static_config.min_qp_allowed, (int32_t)sequence_control_set_ptr->static_config.max_qp_allowed, picture_control_set_ptr->parent_pcs_ptr->picture_qp);
                    frm_hdr->quantization_params.base_q_idx = quantizer_to_qindex[picture_control_set_ptr->picture_qp];
                }

                picture_control_set_ptr->parent_pcs_ptr->picture_qp = picture_control_set_ptr->picture_qp;
                setup_segmentation(
                        picture_control_set_ptr,
                        sequence_control_set_ptr,
                        rate_control_layer_ptr
                );
            }
            else {
                // ***Rate Control***
                if (sequence_control_set_ptr->static_config.rate_control_mode == 1)
                    picture_control_set_ptr->picture_qp = rate_control_get_quantizer(rc_model_ptr, picture_control_set_ptr->parent_pcs_ptr);
                else if (sequence_control_set_ptr->static_config.rate_control_mode == 2) {
                    frame_level_rc_input_picture_vbr(
                        picture_control_set_ptr,
                        sequence_control_set_ptr,
                        context_ptr,
                        rate_control_layer_ptr,
                        rate_control_param_ptr);

                    // rate control QP refinement
                    rate_control_refinement(
                        picture_control_set_ptr,
                        sequence_control_set_ptr,
                        rate_control_param_ptr,
                        prev_gop_rate_control_param_ptr,
                        next_gop_rate_control_param_ptr);
                }
                else if (sequence_control_set_ptr->static_config.rate_control_mode == 3) {
                    frame_level_rc_input_picture_cvbr(
                        picture_control_set_ptr,
                        sequence_control_set_ptr,
                        context_ptr,
                        rate_control_layer_ptr,
                        rate_control_param_ptr);
                }
                picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                    sequence_control_set_ptr->static_config.min_qp_allowed,
                    sequence_control_set_ptr->static_config.max_qp_allowed,
                    picture_control_set_ptr->picture_qp);
                frm_hdr->quantization_params.base_q_idx = quantizer_to_qindex[picture_control_set_ptr->picture_qp];
            }

            picture_control_set_ptr->parent_pcs_ptr->picture_qp = picture_control_set_ptr->picture_qp;

            if (picture_control_set_ptr->parent_pcs_ptr->temporal_layer_index == 0 && sequence_control_set_ptr->static_config.look_ahead_distance != 0)
                context_ptr->base_layer_frames_avg_qp = (3 * context_ptr->base_layer_frames_avg_qp + picture_control_set_ptr->picture_qp + 2) >> 2;
            if (picture_control_set_ptr->slice_type == I_SLICE) {
                if (picture_control_set_ptr->picture_number == rate_control_param_ptr->first_poc) {
                    rate_control_param_ptr->first_pic_pred_qp = (uint16_t)picture_control_set_ptr->parent_pcs_ptr->best_pred_qp;
                    rate_control_param_ptr->first_pic_actual_qp = (uint16_t)picture_control_set_ptr->picture_qp;
                    rate_control_param_ptr->scene_change_in_gop = picture_control_set_ptr->parent_pcs_ptr->scene_change_in_gop;
                    rate_control_param_ptr->first_pic_actual_qp_assigned = EB_TRUE;
                }
                {
                    if (picture_control_set_ptr->picture_number == rate_control_param_ptr->first_poc) {
                        if (sequence_control_set_ptr->static_config.look_ahead_distance != 0)
                            context_ptr->base_layer_intra_frames_avg_qp = (3 * context_ptr->base_layer_intra_frames_avg_qp + picture_control_set_ptr->picture_qp + 2) >> 2;
                    }

                    if (picture_control_set_ptr->picture_number == rate_control_param_ptr->first_poc) {
                        rate_control_param_ptr->intra_frames_qp = picture_control_set_ptr->picture_qp;
                        rate_control_param_ptr->next_gop_intra_frame_qp = picture_control_set_ptr->picture_qp;
                        rate_control_param_ptr->intra_frames_qp_bef_scal = (uint8_t)sequence_control_set_ptr->static_config.max_qp_allowed;
                        for (uint32_t qindex = sequence_control_set_ptr->static_config.min_qp_allowed; qindex <= sequence_control_set_ptr->static_config.max_qp_allowed; qindex++) {
                            if (rate_control_param_ptr->intra_frames_qp <= context_ptr->qp_scaling_map_I_SLICE[qindex]) {
                                rate_control_param_ptr->intra_frames_qp_bef_scal = (uint8_t)qindex;
                                break;
                            }
                        }
                    }
                }
            }

            picture_control_set_ptr->parent_pcs_ptr->average_qp = 0;
            LargestCodingUnit         *sb_ptr;
            uint32_t                       lcuCodingOrder;
            for (lcuCodingOrder = 0; lcuCodingOrder < sequence_control_set_ptr->sb_tot_cnt; ++lcuCodingOrder) {
                sb_ptr = picture_control_set_ptr->sb_ptr_array[lcuCodingOrder];
#if ADD_DELTA_QP_SUPPORT

                sb_ptr->qp = quantizer_to_qindex[picture_control_set_ptr->picture_qp];
#else
                sb_ptr->qp = (uint8_t)picture_control_set_ptr->picture_qp;
#endif
                picture_control_set_ptr->parent_pcs_ptr->average_qp += sb_ptr->qp;
            }

            // Get Empty Rate Control Results Buffer
            eb_get_empty_object(
                context_ptr->rate_control_output_results_fifo_ptr,
                &rate_control_results_wrapper_ptr);
            rate_control_results_ptr = (RateControlResults*)rate_control_results_wrapper_ptr->object_ptr;
            rate_control_results_ptr->picture_control_set_wrapper_ptr = rate_control_tasks_ptr->picture_control_set_wrapper_ptr;

            // Post Full Rate Control Results
            eb_post_full_object(rate_control_results_wrapper_ptr);

            // Release Rate Control Tasks
            eb_release_object(rate_control_tasks_wrapper_ptr);

            break;

        case RC_PACKETIZATION_FEEDBACK_RESULT:

            parentpicture_control_set_ptr = (PictureParentControlSet  *)rate_control_tasks_ptr->picture_control_set_wrapper_ptr->object_ptr;
            sequence_control_set_ptr = (SequenceControlSet *)parentpicture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr;
            if (sequence_control_set_ptr->static_config.rate_control_mode) {
                if (sequence_control_set_ptr->static_config.rate_control_mode == 1)
                    rate_control_update_model(rc_model_ptr, parentpicture_control_set_ptr);
                ReferenceQueueEntry           *reference_entry_ptr;
                uint32_t                          reference_queue_index;
                EncodeContext             *encode_context_ptr = sequence_control_set_ptr->encode_context_ptr;
                reference_queue_index = encode_context_ptr->reference_picture_queue_head_index;
                // Find the Reference in the Reference Queue
                do {
                    reference_entry_ptr = encode_context_ptr->reference_picture_queue[reference_queue_index];
                    if (reference_entry_ptr->picture_number == parentpicture_control_set_ptr->picture_number) {
                        // Set the feedback arrived
                        reference_entry_ptr->feedback_arrived = EB_TRUE;
                    }

                    // Increment the reference_queue_index Iterator
                    reference_queue_index = (reference_queue_index == REFERENCE_QUEUE_MAX_DEPTH - 1) ? 0 : reference_queue_index + 1;
                } while ((reference_queue_index != encode_context_ptr->reference_picture_queue_tail_index) && (reference_entry_ptr->picture_number != parentpicture_control_set_ptr->picture_number));
            }
            // Frame level RC
            if (sequence_control_set_ptr->intra_period_length == -1 || sequence_control_set_ptr->static_config.rate_control_mode == 0) {
                rate_control_param_ptr = context_ptr->rate_control_param_queue[0];
                prev_gop_rate_control_param_ptr = context_ptr->rate_control_param_queue[0];
                if (parentpicture_control_set_ptr->slice_type == I_SLICE) {
                    if (parentpicture_control_set_ptr->total_num_bits > MAX_BITS_PER_FRAME)
                        context_ptr->max_rate_adjust_delta_qp++;
                    else if (context_ptr->max_rate_adjust_delta_qp > 0 && parentpicture_control_set_ptr->total_num_bits < MAX_BITS_PER_FRAME * 85 / 100)
                        context_ptr->max_rate_adjust_delta_qp--;
                    context_ptr->max_rate_adjust_delta_qp = CLIP3(0, 63, context_ptr->max_rate_adjust_delta_qp);
                    context_ptr->max_rate_adjust_delta_qp = 0;
                }
            }
            else {
                uint32_t interval_index_temp = 0;
                EbBool intervalFound = EB_FALSE;
                while ((interval_index_temp < PARALLEL_GOP_MAX_NUMBER) && !intervalFound) {
                    if (parentpicture_control_set_ptr->picture_number >= context_ptr->rate_control_param_queue[interval_index_temp]->first_poc &&
                        parentpicture_control_set_ptr->picture_number <= context_ptr->rate_control_param_queue[interval_index_temp]->last_poc) {
                        intervalFound = EB_TRUE;
                    }
                    else
                        interval_index_temp++;
                }
                CHECK_REPORT_ERROR(
                    interval_index_temp != PARALLEL_GOP_MAX_NUMBER,
                    sequence_control_set_ptr->encode_context_ptr->app_callback_ptr,
                    EB_ENC_RC_ERROR2);

                rate_control_param_ptr = context_ptr->rate_control_param_queue[interval_index_temp];

                prev_gop_rate_control_param_ptr = (interval_index_temp == 0) ?
                    context_ptr->rate_control_param_queue[PARALLEL_GOP_MAX_NUMBER - 1] :
                    context_ptr->rate_control_param_queue[interval_index_temp - 1];
            }
            if (sequence_control_set_ptr->static_config.rate_control_mode != 0) {
                context_ptr->previous_virtual_buffer_level = context_ptr->virtual_buffer_level;

                context_ptr->virtual_buffer_level =
                    (int64_t)context_ptr->previous_virtual_buffer_level +
                    (int64_t)parentpicture_control_set_ptr->total_num_bits - (int64_t)context_ptr->high_level_rate_control_ptr->channel_bit_rate_per_frame;

                high_level_rc_feed_back_picture(
                    parentpicture_control_set_ptr,
                    sequence_control_set_ptr);
                if (sequence_control_set_ptr->static_config.rate_control_mode == 2)
                    frame_level_rc_feedback_picture_vbr(
                        parentpicture_control_set_ptr,
                        sequence_control_set_ptr,
                        context_ptr);
                else if (sequence_control_set_ptr->static_config.rate_control_mode == 3)
                    frame_level_rc_feedback_picture_cvbr(
                        parentpicture_control_set_ptr,
                        sequence_control_set_ptr,
                        context_ptr);
                if (parentpicture_control_set_ptr->picture_number == rate_control_param_ptr->first_poc) {
                    rate_control_param_ptr->first_pic_pred_bits = parentpicture_control_set_ptr->target_bits_best_pred_qp;
                    rate_control_param_ptr->first_pic_actual_bits = parentpicture_control_set_ptr->total_num_bits;
                    {
                        int16_t deltaApQp = (int16_t)rate_control_param_ptr->first_pic_actual_qp - (int16_t)rate_control_param_ptr->first_pic_pred_qp;
                        rate_control_param_ptr->extra_ap_bit_ratio_i = (rate_control_param_ptr->first_pic_pred_bits != 0) ?
                            (((int64_t)rate_control_param_ptr->first_pic_actual_bits - (int64_t)rate_control_param_ptr->first_pic_pred_bits) * 100) / ((int64_t)rate_control_param_ptr->first_pic_pred_bits) :
                            0;
                        rate_control_param_ptr->extra_ap_bit_ratio_i += (int64_t)deltaApQp * 15;
                    }
                }
            }

            // Queue variables
#if OVERSHOOT_STAT_PRINT
            if (sequence_control_set_ptr->intra_period_length != -1) {
                int32_t                       queueEntryIndex;
                uint32_t                       queueEntryIndexTemp;
                uint32_t                       queueEntryIndexTemp2;
                CodedFramesStatsEntry     *queueEntryPtr;
                EbBool                      moveSlideWondowFlag = EB_TRUE;
                EbBool                      end_of_sequence_flag = EB_TRUE;
                uint32_t                       frames_in_sw;

                // Determine offset from the Head Ptr
                queue_entry_index = (int32_t)(parentpicture_control_set_ptr->picture_number - context_ptr->coded_frames_stat_queue[context_ptr->coded_frames_stat_queue_head_index]->picture_number);
                queue_entry_index += context_ptr->coded_frames_stat_queue_head_index;
                queue_entry_index = (queue_entry_index > CODED_FRAMES_STAT_QUEUE_MAX_DEPTH - 1) ? queue_entry_index - CODED_FRAMES_STAT_QUEUE_MAX_DEPTH : queue_entry_index;
                queue_entry_ptr = context_ptr->coded_frames_stat_queue[queue_entry_index];

                queue_entry_ptr->frame_total_bit_actual = (uint64_t)parentpicture_control_set_ptr->total_num_bits;
                queue_entry_ptr->picture_number = parentpicture_control_set_ptr->picture_number;
                queue_entry_ptr->end_of_sequence_flag = parentpicture_control_set_ptr->end_of_sequence_flag;
                context_ptr->rate_average_periodin_frames = (uint64_t)sequence_control_set_ptr->static_config.intra_period + 1;

                //SVT_LOG("\n0_POC: %d\n",
                //    queue_entry_ptr->picture_number);
                move_slide_wondow_flag = EB_TRUE;
                while (move_slide_wondow_flag) {
                    //  SVT_LOG("\n1_POC: %d\n",
                    //      queue_entry_ptr->picture_number);
                      // Check if the sliding window condition is valid
                    queue_entry_index_temp = context_ptr->coded_frames_stat_queue_head_index;
                    if (context_ptr->coded_frames_stat_queue[queue_entry_index_temp]->frame_total_bit_actual != -1)
                        end_of_sequence_flag = context_ptr->coded_frames_stat_queue[queue_entry_index_temp]->end_of_sequence_flag;
                    else
                        end_of_sequence_flag = EB_FALSE;
                    while (move_slide_wondow_flag && !end_of_sequence_flag &&
                        queue_entry_index_temp < context_ptr->coded_frames_stat_queue_head_index + context_ptr->rate_average_periodin_frames) {
                        // SVT_LOG("\n2_POC: %d\n",
                        //     queue_entry_ptr->picture_number);

                        queue_entry_index_temp2 = (queue_entry_index_temp > CODED_FRAMES_STAT_QUEUE_MAX_DEPTH - 1) ? queue_entry_index_temp - CODED_FRAMES_STAT_QUEUE_MAX_DEPTH : queue_entry_index_temp;

                        move_slide_wondow_flag = (EbBool)(move_slide_wondow_flag && (context_ptr->coded_frames_stat_queue[queue_entry_index_temp2]->frame_total_bit_actual != -1));

                        if (context_ptr->coded_frames_stat_queue[queue_entry_index_temp2]->frame_total_bit_actual != -1) {
                            // check if it is the last frame. If we have reached the last frame, we would output the buffered frames in the Queue.
                            end_of_sequence_flag = context_ptr->coded_frames_stat_queue[queue_entry_index_temp]->end_of_sequence_flag;
                        }
                        else
                            end_of_sequence_flag = EB_FALSE;
                        queue_entry_index_temp =
                            (queue_entry_index_temp == CODED_FRAMES_STAT_QUEUE_MAX_DEPTH - 1) ? 0 : queue_entry_index_temp + 1;
                    }

                    if (move_slide_wondow_flag) {
                        //get a new entry spot
                        queue_entry_ptr = (context_ptr->coded_frames_stat_queue[context_ptr->coded_frames_stat_queue_head_index]);
                        queue_entry_index_temp = context_ptr->coded_frames_stat_queue_head_index;
                        // This is set to false, so the last frame would go inside the loop
                        end_of_sequence_flag = EB_FALSE;
                        frames_in_sw = 0;
                        context_ptr->total_bit_actual_per_sw = 0;

                        while (!end_of_sequence_flag &&
                            queue_entry_index_temp < context_ptr->coded_frames_stat_queue_head_index + context_ptr->rate_average_periodin_frames) {
                            frames_in_sw++;

                            queue_entry_index_temp2 = (queue_entry_index_temp > CODED_FRAMES_STAT_QUEUE_MAX_DEPTH - 1) ? queue_entry_index_temp - CODED_FRAMES_STAT_QUEUE_MAX_DEPTH : queue_entry_index_temp;

                            context_ptr->total_bit_actual_per_sw += context_ptr->coded_frames_stat_queue[queue_entry_index_temp2]->frame_total_bit_actual;
                            end_of_sequence_flag = context_ptr->coded_frames_stat_queue[queue_entry_index_temp2]->end_of_sequence_flag;

                            queue_entry_index_temp =
                                (queue_entry_index_temp == CODED_FRAMES_STAT_QUEUE_MAX_DEPTH - 1) ? 0 : queue_entry_index_temp + 1;
                        }
                        //

                        //if(frames_in_sw == context_ptr->rate_average_periodin_frames)
                        //    SVT_LOG("POC:%d\t %.3f\n", queue_entry_ptr->picture_number, (double)context_ptr->total_bit_actual_per_sw*(sequence_control_set_ptr->frame_rate>> RC_PRECISION)/(double)frames_in_sw/1000);
                        if (frames_in_sw == (uint32_t)sequence_control_set_ptr->intra_period_length + 1) {
                            context_ptr->max_bit_actual_per_sw = MAX(context_ptr->max_bit_actual_per_sw, context_ptr->total_bit_actual_per_sw*(sequence_control_set_ptr->frame_rate >> RC_PRECISION) / frames_in_sw / 1000);
                            if (queue_entry_ptr->picture_number % ((sequence_control_set_ptr->intra_period_length + 1)) == 0) {
                                context_ptr->max_bit_actual_per_gop = MAX(context_ptr->max_bit_actual_per_gop, context_ptr->total_bit_actual_per_sw*(sequence_control_set_ptr->frame_rate >> RC_PRECISION) / frames_in_sw / 1000);
                                context_ptr->min_bit_actual_per_gop = MIN(context_ptr->min_bit_actual_per_gop, context_ptr->total_bit_actual_per_sw*(sequence_control_set_ptr->frame_rate >> RC_PRECISION) / frames_in_sw / 1000);
                                if (1) {
                                    //if (context_ptr->total_bit_actual_per_sw > sequence_control_set_ptr->static_config.max_buffersize){
                                    SVT_LOG("POC:%d\t%.0f\t%.2f%% \n",
                                        (int)queue_entry_ptr->picture_number,
                                        (double)((int64_t)context_ptr->total_bit_actual_per_sw*(sequence_control_set_ptr->frame_rate >> RC_PRECISION) / frames_in_sw / 1000),
                                        (double)(100 * (double)context_ptr->total_bit_actual_per_sw*(sequence_control_set_ptr->frame_rate >> RC_PRECISION) / frames_in_sw / (double)sequence_control_set_ptr->static_config.target_bit_rate) - 100);
                                }
                            }
                        }
                        if (frames_in_sw == context_ptr->rate_average_periodin_frames - 1) {
                            //SVT_LOG("\n%d MAX\n", (int32_t)context_ptr->max_bit_actual_per_sw);
                            SVT_LOG("\n%d GopMax\t", (int32_t)context_ptr->max_bit_actual_per_gop);
                            SVT_LOG("%d GopMin\n", (int32_t)context_ptr->min_bit_actual_per_gop);
                        }
                        // Reset the Queue Entry
                        queue_entry_ptr->picture_number += CODED_FRAMES_STAT_QUEUE_MAX_DEPTH;
                        queue_entry_ptr->frame_total_bit_actual = -1;

                        // Increment the Reorder Queue head Ptr
                        context_ptr->coded_frames_stat_queue_head_index =
                            (context_ptr->coded_frames_stat_queue_head_index == CODED_FRAMES_STAT_QUEUE_MAX_DEPTH - 1) ? 0 : context_ptr->coded_frames_stat_queue_head_index + 1;

                        queue_entry_ptr = (context_ptr->coded_frames_stat_queue[context_ptr->coded_frames_stat_queue_head_index]);
                    }
                }
            }
#endif
            total_number_of_fb_frames++;

            // Release the SequenceControlSet
            eb_release_object(parentpicture_control_set_ptr->sequence_control_set_wrapper_ptr);
            // Release the ParentPictureControlSet
            eb_release_object(parentpicture_control_set_ptr->input_picture_wrapper_ptr);
            eb_release_object(rate_control_tasks_ptr->picture_control_set_wrapper_ptr);

            // Release Rate Control Tasks
            eb_release_object(rate_control_tasks_wrapper_ptr);
            break;

        case RC_ENTROPY_CODING_ROW_FEEDBACK_RESULT:

            // Extract bits-per-lcu-row

            // Release Rate Control Tasks
            eb_release_object(rate_control_tasks_wrapper_ptr);

            break;

        default:
            picture_control_set_ptr = (PictureControlSet*)rate_control_tasks_ptr->picture_control_set_wrapper_ptr->object_ptr;
            sequence_control_set_ptr = (SequenceControlSet*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr;

            break;
        }
    }
    EB_DELETE(rc_model_ptr);
    return EB_NULL;
}
