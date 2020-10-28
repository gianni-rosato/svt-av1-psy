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
#include "EbTransforms.h"
#include "aom_dsp_rtcd.h"
#include "EbRateDistortionCost.h"
#include "EbLog.h"
#include "EbIntraPrediction.h"
#include "EbMotionEstimation.h"
/**************************************
 * Context
 **************************************/
typedef struct InitialRateControlContext {
    EbFifo *motion_estimation_results_input_fifo_ptr;
    EbFifo *initialrate_control_results_output_fifo_ptr;
} InitialRateControlContext;

/**************************************
* Macros
**************************************/
static void initial_rate_control_context_dctor(EbPtr p) {
    EbThreadContext *          thread_context_ptr = (EbThreadContext *)p;
    InitialRateControlContext *obj = (InitialRateControlContext *)thread_context_ptr->priv;
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
    context_ptr->initialrate_control_results_output_fifo_ptr = svt_system_resource_get_producer_fifo(
        enc_handle_ptr->initial_rate_control_results_resource_ptr, 0);

    return EB_ErrorNone;
}

/************************************************
* Release Pa Reference Objects
** Check if reference pictures are needed
** release them when appropriate
************************************************/
void release_pa_reference_objects(SequenceControlSet *scs_ptr, PictureParentControlSet *pcs_ptr) {
    // PA Reference Pictures
    if (pcs_ptr->slice_type != I_SLICE) {
        uint32_t num_of_list_to_search = (pcs_ptr->slice_type == P_SLICE) ? REF_LIST_0 : REF_LIST_1;

        // List Loop
        for (uint32_t list_index = REF_LIST_0; list_index <= num_of_list_to_search; ++list_index) {
            // Release PA Reference Pictures
            uint8_t num_of_ref_pic_to_search =
                (pcs_ptr->slice_type == P_SLICE)
                    ? MIN(pcs_ptr->ref_list0_count, scs_ptr->reference_count)
                    : (list_index == REF_LIST_0)
                          ? MIN(pcs_ptr->ref_list0_count, scs_ptr->reference_count)
                          : MIN(pcs_ptr->ref_list1_count, scs_ptr->reference_count);

            for (uint32_t ref_pic_index = 0; ref_pic_index < num_of_ref_pic_to_search;
                 ++ref_pic_index) {
                if (pcs_ptr->ref_pa_pic_ptr_array[list_index][ref_pic_index] != NULL) {
                    svt_release_object(pcs_ptr->ref_pa_pic_ptr_array[list_index][ref_pic_index]);
                }
            }
        }
    }

    if (pcs_ptr->pa_reference_picture_wrapper_ptr != NULL) {
        svt_release_object(pcs_ptr->pa_reference_picture_wrapper_ptr);
    }

    return;
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
        uint16_t non_moving_index_over_sliding_window = pcs_ptr->non_moving_index_array[sb_idx];
        uint16_t frames_to_check_index;
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

            if (temp_pcs_ptr->slice_type == I_SLICE || temp_pcs_ptr->end_of_sequence_flag) break;
            // Limit the distortion to lower layers 0, 1 and 2 only. Higher layers have close temporal distance and lower distortion that might contaminate the data
            if (sb_params->is_complete_sb && temp_pcs_ptr->temporal_layer_index <
                MAX((int8_t)pcs_ptr->hierarchical_levels - 1, 2)) {
                if (sb_idx == 0) me_dist_pic_count++;
                me_dist += (temp_pcs_ptr->slice_type == I_SLICE)
                               ? 0
                               : (uint64_t)temp_pcs_ptr->rc_me_distortion[sb_idx];
            }
            // Store the filtered_sse of next ALT_REF picture in the I slice to be used in QP Scaling
            if (pcs_ptr->slice_type == I_SLICE && pcs_ptr->filtered_sse == 0 && sb_idx == 0 &&
                temp_pcs_ptr->temporal_layer_index == 0) {
                pcs_ptr->filtered_sse    = temp_pcs_ptr->filtered_sse;
                pcs_ptr->filtered_sse_uv = temp_pcs_ptr->filtered_sse_uv;
            }
            non_moving_index_over_sliding_window += temp_pcs_ptr->non_moving_index_array[sb_idx];

            // Increment the input_queue_index Iterator
            input_queue_index =
                (input_queue_index == INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH - 1)
                    ? 0
                    : input_queue_index + 1;
        }
        pcs_ptr->non_moving_index_array[sb_idx] =
            (uint8_t)(non_moving_index_over_sliding_window / (frames_to_check_index + 1));

        non_moving_index_sum += pcs_ptr->non_moving_index_array[sb_idx];
    }
    assert(complete_sb_count > 0);
    pcs_ptr->non_moving_index_average = (uint16_t)non_moving_index_sum / pcs_ptr->sb_total_count;
    me_dist_pic_count                 = MAX(me_dist_pic_count, 1);
    pcs_ptr->qp_scaling_average_complexity =
        (uint16_t)((uint64_t)me_dist / complete_sb_count / 256 / me_dist_pic_count);
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

        if (temp_pcs_ptr->end_of_sequence_flag) break;
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

    queue_entry_index =
        (int32_t)(pcs_ptr->picture_number -
                  encode_context_ptr
                      ->initial_rate_control_reorder_queue
                          [encode_context_ptr->initial_rate_control_reorder_queue_head_index]
                      ->picture_number);
    queue_entry_index += encode_context_ptr->initial_rate_control_reorder_queue_head_index;
    queue_entry_index = (queue_entry_index > INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH - 1)
                            ? queue_entry_index - INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH
                            : queue_entry_index;
    queue_entry_ptr = encode_context_ptr->initial_rate_control_reorder_queue[queue_entry_index];
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
    histogram_queue_entry_index =
        (histogram_queue_entry_index > HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1)
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
    histogram_queue_entry_index =
        (histogram_queue_entry_index > HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1)
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

static void generate_lambda_scaling_factor(PictureParentControlSet         *pcs_ptr, int64_t mc_dep_cost_base)
{
    Av1Common *cm = pcs_ptr->av1_cm;
    const int step = 1 << (pcs_ptr->is_720p_or_larger ? 2 : 1);
    const int mi_cols_sr = ((pcs_ptr->aligned_width + 15) / 16) << 2;

    const int block_size = BLOCK_16X16;
    const int num_mi_w = mi_size_wide[block_size];
    const int num_mi_h = mi_size_high[block_size];
    const int num_cols = (mi_cols_sr + num_mi_w - 1) / num_mi_w;
    const int num_rows = (cm->mi_rows + num_mi_h - 1) / num_mi_h;
    const int stride   = mi_cols_sr >> (1 + pcs_ptr->is_720p_or_larger);
    const double c = 1.2;

    for (int row = 0; row < num_rows; row++) {
        for (int col = 0; col < num_cols; col++) {
            double intra_cost = 0.0;
            double mc_dep_cost = 0.0;
            const int index = row * num_cols + col;
            for (int mi_row = row * num_mi_h; mi_row < (row + 1) * num_mi_h; mi_row += step) {
                for (int mi_col = col * num_mi_w; mi_col < (col + 1) * num_mi_w; mi_col += step) {
                    if (mi_row >= cm->mi_rows || mi_col >= mi_cols_sr) continue;

                    const int index1 = (mi_row >> (1 + pcs_ptr->is_720p_or_larger)) * stride + (mi_col >> (1 + pcs_ptr->is_720p_or_larger));
                    TplStats *tpl_stats_ptr = pcs_ptr->tpl_stats[index1];
                    int64_t mc_dep_delta =
                        RDCOST(pcs_ptr->base_rdmult, tpl_stats_ptr->mc_dep_rate, tpl_stats_ptr->mc_dep_dist);
                    intra_cost  += (double)(tpl_stats_ptr->recrf_dist << RDDIV_BITS);
                    mc_dep_cost += (double)(tpl_stats_ptr->recrf_dist << RDDIV_BITS) + mc_dep_delta;
                }
            }
            double rk = 0;
            if (mc_dep_cost > 0 && intra_cost > 0) {
                rk = intra_cost / mc_dep_cost;
            }

            pcs_ptr->tpl_rdmult_scaling_factors[index] = (mc_dep_cost_base) ? rk / pcs_ptr->r0 + c : c;
        }
    }

    return;
}

static AOM_INLINE void get_quantize_error(MacroblockPlane *p,
                                          const TranLow *coeff, TranLow *qcoeff,
                                          TranLow *dqcoeff, TxSize tx_size,
                                          uint16_t *eob, int64_t *recon_error,
                                          int64_t *sse) {
  const ScanOrder *const scan_order = &av1_scan_orders[tx_size][DCT_DCT]; //&av1_default_scan_orders[tx_size]
  int pix_num = 1 << num_pels_log2_lookup[txsize_to_bsize[tx_size]];
  const int shift = tx_size == TX_32X32 ? 0 : 2;

  svt_av1_quantize_fp(coeff, pix_num, p->zbin_qtx, p->round_fp_qtx, p->quant_fp_qtx,
                      p->quant_shift_qtx, qcoeff, dqcoeff, p->dequant_qtx, eob,
                      scan_order->scan, scan_order->iscan);

  *recon_error = svt_av1_block_error(coeff, dqcoeff, pix_num, sse) >> shift;
  *recon_error = AOMMAX(*recon_error, 1);

  *sse = (*sse) >> shift;
  *sse = AOMMAX(*sse, 1);
}

static int rate_estimator(TranLow *qcoeff, int eob, TxSize tx_size) {
  const ScanOrder *const scan_order = &av1_scan_orders[tx_size][DCT_DCT]; //&av1_default_scan_orders[tx_size]

  assert((1 << num_pels_log2_lookup[txsize_to_bsize[tx_size]]) >= eob);

  int rate_cost = 1;

  for (int idx = 0; idx < eob; ++idx) {
    int abs_level = abs(qcoeff[scan_order->scan[idx]]);
    rate_cost += (int)(log1p(abs_level) / log(2.0)) + 1;
  }

  return (rate_cost << AV1_PROB_COST_SHIFT);
}

static void result_model_store(PictureParentControlSet *pcs_ptr, TplStats  *tpl_stats_ptr,
        uint32_t mb_origin_x, uint32_t mb_origin_y) {
  const int mi_height = mi_size_high[BLOCK_16X16];
  const int mi_width = mi_size_wide[BLOCK_16X16];
  const int step = 1 << (pcs_ptr->is_720p_or_larger ? 2 : 1);
  const int shift = 3 + pcs_ptr->is_720p_or_larger;
  const int aligned16_width = ((pcs_ptr->aligned_width + 15) / 16) << 4;

  int64_t srcrf_dist = tpl_stats_ptr->srcrf_dist / (mi_height * mi_width);
  int64_t recrf_dist = tpl_stats_ptr->recrf_dist / (mi_height * mi_width);
  int64_t srcrf_rate = tpl_stats_ptr->srcrf_rate / (mi_height * mi_width);
  int64_t recrf_rate = tpl_stats_ptr->recrf_rate / (mi_height * mi_width);

  srcrf_dist = AOMMAX(1, srcrf_dist);
  recrf_dist = AOMMAX(1, recrf_dist);
  srcrf_rate = AOMMAX(1, srcrf_rate);
  recrf_rate = AOMMAX(1, recrf_rate);

  for (int idy = 0; idy < mi_height; idy += step) {
    TplStats *dst_ptr = pcs_ptr->tpl_stats[((mb_origin_y >> shift) + (idy >> 1)) * (aligned16_width >> shift) + (mb_origin_x >> shift)];
    for (int idx = 0; idx < mi_width; idx += step) {
      dst_ptr->srcrf_dist = srcrf_dist;
      dst_ptr->recrf_dist = recrf_dist;
      dst_ptr->srcrf_rate = srcrf_rate;
      dst_ptr->recrf_rate = recrf_rate;
      dst_ptr->mv = tpl_stats_ptr->mv;
      dst_ptr->ref_frame_poc = tpl_stats_ptr->ref_frame_poc;
      ++dst_ptr;
    }
  }
}

static const int16_t dc_qlookup_QTX[QINDEX_RANGE] = {
  4,    8,    8,    9,    10,  11,  12,  12,  13,  14,  15,   16,   17,   18,
  19,   19,   20,   21,   22,  23,  24,  25,  26,  26,  27,   28,   29,   30,
  31,   32,   32,   33,   34,  35,  36,  37,  38,  38,  39,   40,   41,   42,
  43,   43,   44,   45,   46,  47,  48,  48,  49,  50,  51,   52,   53,   53,
  54,   55,   56,   57,   57,  58,  59,  60,  61,  62,  62,   63,   64,   65,
  66,   66,   67,   68,   69,  70,  70,  71,  72,  73,  74,   74,   75,   76,
  77,   78,   78,   79,   80,  81,  81,  82,  83,  84,  85,   85,   87,   88,
  90,   92,   93,   95,   96,  98,  99,  101, 102, 104, 105,  107,  108,  110,
  111,  113,  114,  116,  117, 118, 120, 121, 123, 125, 127,  129,  131,  134,
  136,  138,  140,  142,  144, 146, 148, 150, 152, 154, 156,  158,  161,  164,
  166,  169,  172,  174,  177, 180, 182, 185, 187, 190, 192,  195,  199,  202,
  205,  208,  211,  214,  217, 220, 223, 226, 230, 233, 237,  240,  243,  247,
  250,  253,  257,  261,  265, 269, 272, 276, 280, 284, 288,  292,  296,  300,
  304,  309,  313,  317,  322, 326, 330, 335, 340, 344, 349,  354,  359,  364,
  369,  374,  379,  384,  389, 395, 400, 406, 411, 417, 423,  429,  435,  441,
  447,  454,  461,  467,  475, 482, 489, 497, 505, 513, 522,  530,  539,  549,
  559,  569,  579,  590,  602, 614, 626, 640, 654, 668, 684,  700,  717,  736,
  755,  775,  796,  819,  843, 869, 896, 925, 955, 988, 1022, 1058, 1098, 1139,
  1184, 1232, 1282, 1336,
};

static const int16_t dc_qlookup_10_QTX[QINDEX_RANGE] = {
  4,    9,    10,   13,   15,   17,   20,   22,   25,   28,   31,   34,   37,
  40,   43,   47,   50,   53,   57,   60,   64,   68,   71,   75,   78,   82,
  86,   90,   93,   97,   101,  105,  109,  113,  116,  120,  124,  128,  132,
  136,  140,  143,  147,  151,  155,  159,  163,  166,  170,  174,  178,  182,
  185,  189,  193,  197,  200,  204,  208,  212,  215,  219,  223,  226,  230,
  233,  237,  241,  244,  248,  251,  255,  259,  262,  266,  269,  273,  276,
  280,  283,  287,  290,  293,  297,  300,  304,  307,  310,  314,  317,  321,
  324,  327,  331,  334,  337,  343,  350,  356,  362,  369,  375,  381,  387,
  394,  400,  406,  412,  418,  424,  430,  436,  442,  448,  454,  460,  466,
  472,  478,  484,  490,  499,  507,  516,  525,  533,  542,  550,  559,  567,
  576,  584,  592,  601,  609,  617,  625,  634,  644,  655,  666,  676,  687,
  698,  708,  718,  729,  739,  749,  759,  770,  782,  795,  807,  819,  831,
  844,  856,  868,  880,  891,  906,  920,  933,  947,  961,  975,  988,  1001,
  1015, 1030, 1045, 1061, 1076, 1090, 1105, 1120, 1137, 1153, 1170, 1186, 1202,
  1218, 1236, 1253, 1271, 1288, 1306, 1323, 1342, 1361, 1379, 1398, 1416, 1436,
  1456, 1476, 1496, 1516, 1537, 1559, 1580, 1601, 1624, 1647, 1670, 1692, 1717,
  1741, 1766, 1791, 1817, 1844, 1871, 1900, 1929, 1958, 1990, 2021, 2054, 2088,
  2123, 2159, 2197, 2236, 2276, 2319, 2363, 2410, 2458, 2508, 2561, 2616, 2675,
  2737, 2802, 2871, 2944, 3020, 3102, 3188, 3280, 3375, 3478, 3586, 3702, 3823,
  3953, 4089, 4236, 4394, 4559, 4737, 4929, 5130, 5347,
};

static const int16_t dc_qlookup_12_QTX[QINDEX_RANGE] = {
  4,     12,    18,    25,    33,    41,    50,    60,    70,    80,    91,
  103,   115,   127,   140,   153,   166,   180,   194,   208,   222,   237,
  251,   266,   281,   296,   312,   327,   343,   358,   374,   390,   405,
  421,   437,   453,   469,   484,   500,   516,   532,   548,   564,   580,
  596,   611,   627,   643,   659,   674,   690,   706,   721,   737,   752,
  768,   783,   798,   814,   829,   844,   859,   874,   889,   904,   919,
  934,   949,   964,   978,   993,   1008,  1022,  1037,  1051,  1065,  1080,
  1094,  1108,  1122,  1136,  1151,  1165,  1179,  1192,  1206,  1220,  1234,
  1248,  1261,  1275,  1288,  1302,  1315,  1329,  1342,  1368,  1393,  1419,
  1444,  1469,  1494,  1519,  1544,  1569,  1594,  1618,  1643,  1668,  1692,
  1717,  1741,  1765,  1789,  1814,  1838,  1862,  1885,  1909,  1933,  1957,
  1992,  2027,  2061,  2096,  2130,  2165,  2199,  2233,  2267,  2300,  2334,
  2367,  2400,  2434,  2467,  2499,  2532,  2575,  2618,  2661,  2704,  2746,
  2788,  2830,  2872,  2913,  2954,  2995,  3036,  3076,  3127,  3177,  3226,
  3275,  3324,  3373,  3421,  3469,  3517,  3565,  3621,  3677,  3733,  3788,
  3843,  3897,  3951,  4005,  4058,  4119,  4181,  4241,  4301,  4361,  4420,
  4479,  4546,  4612,  4677,  4742,  4807,  4871,  4942,  5013,  5083,  5153,
  5222,  5291,  5367,  5442,  5517,  5591,  5665,  5745,  5825,  5905,  5984,
  6063,  6149,  6234,  6319,  6404,  6495,  6587,  6678,  6769,  6867,  6966,
  7064,  7163,  7269,  7376,  7483,  7599,  7715,  7832,  7958,  8085,  8214,
  8352,  8492,  8635,  8788,  8945,  9104,  9275,  9450,  9639,  9832,  10031,
  10245, 10465, 10702, 10946, 11210, 11482, 11776, 12081, 12409, 12750, 13118,
  13501, 13913, 14343, 14807, 15290, 15812, 16356, 16943, 17575, 18237, 18949,
  19718, 20521, 21387,
};

int16_t av1_dc_quant_qtx(int qindex, int delta, AomBitDepth bit_depth) {
  const int q_clamped = clamp(qindex + delta, 0, MAXQ);
  switch (bit_depth) {
    case AOM_BITS_8: return dc_qlookup_QTX[q_clamped];
    case AOM_BITS_10: return dc_qlookup_10_QTX[q_clamped];
    case AOM_BITS_12: return dc_qlookup_12_QTX[q_clamped];
    default:
      assert(0 && "bit_depth should be AOM_BITS_8, AOM_BITS_10 or AOM_BITS_12");
      return -1;
  }
}

int svt_av1_compute_rd_mult_based_on_qindex(AomBitDepth bit_depth, int qindex) {
  const int q = av1_dc_quant_qtx(qindex, 0, bit_depth);
  //const int q = svt_av1_dc_quant_Q3(qindex, 0, bit_depth);
  int rdmult = q * q;
  rdmult = rdmult * 3 + (rdmult * 2 / 3);
  switch (bit_depth) {
    case AOM_BITS_8: break;
    case AOM_BITS_10: rdmult = ROUND_POWER_OF_TWO(rdmult, 4); break;
    case AOM_BITS_12: rdmult = ROUND_POWER_OF_TWO(rdmult, 8); break;
    default:
      assert(0 && "bit_depth should be AOM_BITS_8, AOM_BITS_10 or AOM_BITS_12");
      return -1;
  }
  return rdmult > 0 ? rdmult : 1;
}

void svt_av1_set_quantizer(PictureParentControlSet *pcs_ptr, int32_t q);
void svt_av1_build_quantizer(AomBitDepth bit_depth, int32_t y_dc_delta_q, int32_t u_dc_delta_q,
                             int32_t u_ac_delta_q, int32_t v_dc_delta_q, int32_t v_ac_delta_q,
                             Quants *const quants, Dequants *const deq);

#define TPL_DEP_COST_SCALE_LOG2 4
double svt_av1_convert_qindex_to_q(int32_t qindex, AomBitDepth bit_depth);
int32_t svt_av1_compute_qdelta(double qstart, double qtarget, AomBitDepth bit_depth);
extern void filter_intra_edge(OisMbResults *ois_mb_results_ptr, uint8_t mode, uint16_t max_frame_width, uint16_t max_frame_height,
                            int32_t p_angle, int32_t cu_origin_x, int32_t cu_origin_y, uint8_t *above_row, uint8_t *left_col);

//Given one reference frame identified by the pair (list_index,ref_index)
//indicate if ME data is valid
static uint8_t is_me_data_valid(
    const MeSbResults           *me_results,
    uint32_t                     me_mb_offset,
    uint8_t                      list_idx,
    uint8_t                      ref_idx) {

    uint8_t total_me_cnt = me_results->total_me_candidate_index[me_mb_offset];
    const MeCandidate *me_block_results = &me_results->me_candidate_array[me_mb_offset*MAX_PA_ME_CAND];

    for (uint32_t me_cand_i = 0; me_cand_i < total_me_cnt; ++me_cand_i) {
        const MeCandidate *me_cand = &me_block_results[me_cand_i];
        assert(/*me_cand->direction >= 0 && */me_cand->direction <= 2);
        if (me_cand->direction == 0 || me_cand->direction == 2) {
            if (list_idx == me_cand->ref0_list && ref_idx == me_cand->ref_idx_l0)
                return 1;
        }
        if (me_cand->direction == 1 || me_cand->direction == 2) {
            if (list_idx == me_cand->ref1_list && ref_idx == me_cand->ref_idx_l1)
                return 1;
        }
    }
    return 0;
}
/************************************************
* Genrate TPL MC Flow Dispenser  Based on Lookahead
** LAD Window: sliding window size
************************************************/
void tpl_mc_flow_dispenser(
    EncodeContext                   *encode_context_ptr,
    SequenceControlSet              *scs_ptr,
    PictureParentControlSet         *pcs_ptr,
    int32_t                          frame_idx)
{
    uint32_t    picture_width_in_sb = (pcs_ptr->enhanced_picture_ptr->width + BLOCK_SIZE_64 - 1) / BLOCK_SIZE_64;
    uint32_t    picture_width_in_mb = (pcs_ptr->enhanced_picture_ptr->width + 16 - 1) / 16;
    uint32_t    picture_height_in_sb = (pcs_ptr->enhanced_picture_ptr->height + BLOCK_SIZE_64 - 1) / BLOCK_SIZE_64;
    int16_t     x_curr_mv = 0;
    int16_t     y_curr_mv = 0;
    uint32_t    me_mb_offset = 0;
    TxSize      tx_size = TX_16X16;
    EbPictureBufferDesc  *ref_pic_ptr;
    struct      ScaleFactors sf;
    BlockGeom   blk_geom;
    uint32_t    kernel = (EIGHTTAP_REGULAR << 16) | EIGHTTAP_REGULAR;
    EbPictureBufferDesc *input_picture_ptr = pcs_ptr->enhanced_picture_ptr;
    TplStats  tpl_stats;

    (void)scs_ptr;

    DECLARE_ALIGNED(32, uint8_t, predictor8[256 * 2]);
    DECLARE_ALIGNED(32, int16_t, src_diff[256]);
    DECLARE_ALIGNED(32, TranLow, coeff[256]);
    DECLARE_ALIGNED(32, TranLow, qcoeff[256]);
    DECLARE_ALIGNED(32, TranLow, dqcoeff[256]);
    DECLARE_ALIGNED(32, TranLow, best_coeff[256]);
    uint8_t *predictor = predictor8;

    blk_geom.bwidth  = 16;
    blk_geom.bheight = 16;

    svt_av1_setup_scale_factors_for_frame(
                &sf, picture_width_in_sb * BLOCK_SIZE_64,
                picture_height_in_sb * BLOCK_SIZE_64,
                picture_width_in_sb * BLOCK_SIZE_64,
                picture_height_in_sb * BLOCK_SIZE_64);

    MacroblockPlane mb_plane;
    int32_t qIndex = quantizer_to_qindex[(uint8_t)scs_ptr->static_config.qp];

    const  double delta_rate_new[7][6] =
    {
        { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 }, // 1L
        { 0.6, 1.0, 1.0, 1.0, 1.0, 1.0 }, // 2L
        { 0.6, 0.8, 1.0, 1.0, 1.0, 1.0 }, // 3L
        { 0.6 , 0.8, 0.9, 1.0, 1.0, 1.0 }, // 4L
        { 0.35, 0.6, 0.8,  0.9, 1.0, 1.0},  //5L
        { 0.35, 0.6, 0.8,  0.9, 0.95, 1.0}  //6L
    };
    double q_val;  q_val = svt_av1_convert_qindex_to_q(qIndex, 8);
    int32_t delta_qindex;
    if(pcs_ptr->slice_type == I_SLICE)
        delta_qindex = svt_av1_compute_qdelta(
            q_val,
            q_val * 0.25,
            8);
    else
         delta_qindex = svt_av1_compute_qdelta(
            q_val,
            q_val * delta_rate_new[pcs_ptr->hierarchical_levels]
            [pcs_ptr->temporal_layer_index],
            8);
    qIndex =
       (qIndex + delta_qindex);

    Quants *const quants_bd = &pcs_ptr->quants_bd;
    Dequants *const deq_bd = &pcs_ptr->deq_bd;
    svt_av1_set_quantizer(
        pcs_ptr,
        pcs_ptr->frm_hdr.quantization_params.base_q_idx);
    svt_av1_build_quantizer(
        /*pcs_ptr->hbd_mode_decision ? AOM_BITS_10 :*/ AOM_BITS_8,
        pcs_ptr->frm_hdr.quantization_params.delta_q_dc[AOM_PLANE_Y],
        pcs_ptr->frm_hdr.quantization_params.delta_q_dc[AOM_PLANE_U],
        pcs_ptr->frm_hdr.quantization_params.delta_q_ac[AOM_PLANE_U],
        pcs_ptr->frm_hdr.quantization_params.delta_q_dc[AOM_PLANE_V],
        pcs_ptr->frm_hdr.quantization_params.delta_q_ac[AOM_PLANE_V],
        quants_bd,
        deq_bd);
    mb_plane.quant_qtx       = pcs_ptr->quants_bd.y_quant[qIndex];
    mb_plane.quant_fp_qtx    = pcs_ptr->quants_bd.y_quant_fp[qIndex];
    mb_plane.round_fp_qtx    = pcs_ptr->quants_bd.y_round_fp[qIndex];
    mb_plane.quant_shift_qtx = pcs_ptr->quants_bd.y_quant_shift[qIndex];
    mb_plane.zbin_qtx        = pcs_ptr->quants_bd.y_zbin[qIndex];
    mb_plane.round_qtx       = pcs_ptr->quants_bd.y_round[qIndex];
    mb_plane.dequant_qtx     = pcs_ptr->deq_bd.y_dequant_qtx[qIndex];
    pcs_ptr->base_rdmult = svt_av1_compute_rd_mult_based_on_qindex((AomBitDepth)8/*scs_ptr->static_config.encoder_bit_depth*/, qIndex) / 6;

    // Walk the first N entries in the sliding window
    for (uint32_t sb_index = 0; sb_index < pcs_ptr->sb_total_count; ++sb_index) {
        {
            SbParams *sb_params = &scs_ptr->sb_params_array[sb_index];
            uint32_t pa_blk_index = 0;
            while (pa_blk_index < CU_MAX_COUNT) {
                const CodedBlockStats *blk_stats_ptr;
                blk_stats_ptr = get_coded_blk_stats(pa_blk_index);
                uint8_t bsize = blk_stats_ptr->size;
                EbBool small_boundary_blk = EB_FALSE;

                //if(sb_params->raster_scan_blk_validity[md_scan_to_raster_scan[pa_blk_index]])
                {
                    uint32_t cu_origin_x = sb_params->origin_x + blk_stats_ptr->origin_x;
                    uint32_t cu_origin_y = sb_params->origin_y + blk_stats_ptr->origin_y;
                    if ((blk_stats_ptr->origin_x % 16) == 0 && (blk_stats_ptr->origin_y % 16) == 0 &&
                            ((pcs_ptr->enhanced_picture_ptr->width - cu_origin_x) < 16 || (pcs_ptr->enhanced_picture_ptr->height - cu_origin_y) < 16))
                        small_boundary_blk = EB_TRUE;
                }
                if(bsize != 16 && !small_boundary_blk) {
                    pa_blk_index++;
                    continue;
                }
                if (sb_params->raster_scan_blk_validity[md_scan_to_raster_scan[pa_blk_index]]) {
                    uint32_t mb_origin_x = sb_params->origin_x + blk_stats_ptr->origin_x;
                    uint32_t mb_origin_y = sb_params->origin_y + blk_stats_ptr->origin_y;
                    const int dst_buffer_stride = input_picture_ptr->stride_y;
                    const int dst_mb_offset = mb_origin_y * dst_buffer_stride + mb_origin_x;
                    const int dst_basic_offset = input_picture_ptr->origin_y * input_picture_ptr->stride_y + input_picture_ptr->origin_x;
                    uint8_t *dst_buffer = encode_context_ptr->mc_flow_rec_picture_buffer[frame_idx] + dst_basic_offset + dst_mb_offset;
                    int64_t inter_cost;
                    int64_t recon_error = 1, sse = 1;
                    uint64_t best_ref_poc = 0;
                    int32_t best_rf_idx = -1;
                    int64_t best_inter_cost = INT64_MAX;
                    MV final_best_mv = {0, 0};
                    uint32_t max_inter_ref = MAX_PA_ME_MV;
                    OisMbResults *ois_mb_results_ptr = pcs_ptr->ois_mb_results[(mb_origin_y >> 4) * picture_width_in_mb + (mb_origin_x >> 4)];
                    int64_t best_intra_cost = ois_mb_results_ptr->intra_cost;
                    uint8_t best_mode = DC_PRED;
                    uint8_t *src_mb = input_picture_ptr->buffer_y + input_picture_ptr->origin_x + mb_origin_x +
                                     (input_picture_ptr->origin_y + mb_origin_y) * input_picture_ptr->stride_y;
                    memset(&tpl_stats, 0, sizeof(tpl_stats));
                    blk_geom.origin_x = blk_stats_ptr->origin_x;
                    blk_geom.origin_y = blk_stats_ptr->origin_y;
                    me_mb_offset = get_me_info_index(pcs_ptr->max_number_of_pus_per_sb, &blk_geom, 0, 0);
                    for(uint32_t rf_idx = 0; rf_idx < max_inter_ref; rf_idx++) {
                        uint32_t list_index = rf_idx < 4 ? 0 : 1;
                        uint32_t ref_pic_index = rf_idx >= 4 ? (rf_idx - 4) : rf_idx;

                        if( (list_index == 0 && (ref_pic_index+1) > pcs_ptr->ref_list0_count_try) ||
                            (list_index == 1 && (ref_pic_index+1) > pcs_ptr->ref_list1_count_try) )
                            continue;
                        if( !is_me_data_valid( pcs_ptr->pa_me_data->me_results[sb_index], me_mb_offset, list_index, ref_pic_index))
                            continue;
                        if(!pcs_ptr->ref_pa_pic_ptr_array[list_index][ref_pic_index])
                            continue;
                        uint64_t ref_poc = pcs_ptr->ref_pic_poc_array[list_index][ref_pic_index];
                        uint32_t ref_frame_idx = 0;
                        while(ref_frame_idx < MAX_TPL_LA_SW && encode_context_ptr->poc_map_idx[ref_frame_idx] != ref_poc)
                            ref_frame_idx++;
                        if(ref_frame_idx == MAX_TPL_LA_SW || (int32_t)ref_frame_idx >= frame_idx) {
                            continue;
                        }

                        EbPaReferenceObject * referenceObject = (EbPaReferenceObject*)pcs_ptr->ref_pa_pic_ptr_array[list_index][ref_pic_index]->object_ptr;
                        ref_pic_ptr = (EbPictureBufferDesc*)referenceObject->input_padded_picture_ptr;

                        const int ref_basic_offset = ref_pic_ptr->origin_y * ref_pic_ptr->stride_y + ref_pic_ptr->origin_x;
                        const int ref_mb_offset = mb_origin_y * ref_pic_ptr->stride_y + mb_origin_x;
                        uint8_t *ref_mb = ref_pic_ptr->buffer_y + ref_basic_offset + ref_mb_offset;

                        struct Buf2D ref_buf = { NULL, ref_pic_ptr->buffer_y + ref_basic_offset,
                                                  ref_pic_ptr->width, ref_pic_ptr->height,
                                                  ref_pic_ptr->stride_y };
                        const MeSbResults *me_results = pcs_ptr->pa_me_data->me_results[sb_index];
                        x_curr_mv = me_results->me_mv_array[me_mb_offset * MAX_PA_ME_MV + (list_index ? 4 : 0) + ref_pic_index].x_mv << 1;
                        y_curr_mv = me_results->me_mv_array[me_mb_offset * MAX_PA_ME_MV + (list_index ? 4 : 0) + ref_pic_index].y_mv << 1;
                        InterPredParams inter_pred_params;
                        svt_av1_init_inter_params(&inter_pred_params, 16, 16, mb_origin_y,
                                mb_origin_x, 0, 0, 8, 0, 0,
                                &sf, &ref_buf, kernel);

                        inter_pred_params.conv_params = get_conv_params(0, 0, 0, 8);

                        MV best_mv = {y_curr_mv, x_curr_mv};
                        av1_build_inter_predictor(pcs_ptr->av1_cm,
                                                  ref_mb,
                                                  input_picture_ptr->stride_y,
                                                  predictor,
                                                  16,
                                                  &best_mv,
                                                  mb_origin_x,
                                                  mb_origin_y,
                                                  &inter_pred_params);
                        svt_aom_subtract_block(16, 16, src_diff, 16, src_mb, input_picture_ptr->stride_y, predictor, 16);

                        svt_av1_wht_fwd_txfm(src_diff, 16, coeff, tx_size, 8, 0);

                        inter_cost = svt_aom_satd(coeff, 256);
                        if (inter_cost < best_inter_cost) {
                            memcpy(best_coeff, coeff, sizeof(best_coeff));
                            best_ref_poc = ref_poc;
                            best_rf_idx = rf_idx;
                            best_inter_cost = inter_cost;
                            final_best_mv = best_mv;

                            if (best_inter_cost < best_intra_cost) best_mode = NEWMV;
                        }
                    } // rf_idx
                    if(best_inter_cost < INT64_MAX) {
                        uint16_t eob = 0;
                        get_quantize_error(&mb_plane, best_coeff, qcoeff, dqcoeff, tx_size, &eob, &recon_error, &sse);
                        int rate_cost = pcs_ptr->tpl_opt_flag? 0 : rate_estimator(qcoeff, eob, tx_size);
                        tpl_stats.srcrf_rate = rate_cost << TPL_DEP_COST_SCALE_LOG2;
                    }
                    best_intra_cost = AOMMAX(best_intra_cost, 1);
                    if (frame_is_intra_only(pcs_ptr))
                        best_inter_cost = 0;
                    else
                        best_inter_cost = AOMMIN(best_intra_cost, best_inter_cost);

                    tpl_stats.srcrf_dist = recon_error << (TPL_DEP_COST_SCALE_LOG2);

                    if (best_mode == NEWMV) {
                        // inter recon with rec_picture as reference pic
                        uint64_t ref_poc = best_ref_poc;
                        uint32_t ref_frame_idx = 0;
                        while(ref_frame_idx < MAX_TPL_LA_SW && encode_context_ptr->poc_map_idx[ref_frame_idx] != ref_poc)
                            ref_frame_idx++;
                        assert(ref_frame_idx != MAX_TPL_LA_SW);

                        const int ref_basic_offset = input_picture_ptr->origin_y * input_picture_ptr->stride_y + input_picture_ptr->origin_x;
                        const int ref_mb_offset = mb_origin_y * input_picture_ptr->stride_y + mb_origin_x;
                        uint8_t *ref_mb = encode_context_ptr->mc_flow_rec_picture_buffer[ref_frame_idx] + ref_basic_offset + ref_mb_offset;

                        struct Buf2D ref_buf = { NULL, encode_context_ptr->mc_flow_rec_picture_buffer[ref_frame_idx] + ref_basic_offset,
                                                  input_picture_ptr->width, input_picture_ptr->height,
                                                  input_picture_ptr->stride_y};
                        InterPredParams inter_pred_params;
                        svt_av1_init_inter_params(&inter_pred_params, 16, 16, mb_origin_y,
                            mb_origin_x, 0, 0, 8, 0, 0,
                            &sf, &ref_buf, kernel);

                        inter_pred_params.conv_params = get_conv_params(0, 0, 0, 8);
                        av1_build_inter_predictor(pcs_ptr->av1_cm,
                                                  ref_mb,
                                                  input_picture_ptr->stride_y,
                                                  dst_buffer,
                                                  dst_buffer_stride,
                                                  &final_best_mv,
                                                  mb_origin_x,
                                                  mb_origin_y,
                                                  &inter_pred_params);
                    } else {
                        // intra recon
                        uint8_t *above_row;
                        uint8_t *left_col;
                        DECLARE_ALIGNED(16, uint8_t, left_data[MAX_TX_SIZE * 2 + 32]);
                        DECLARE_ALIGNED(16, uint8_t, above_data[MAX_TX_SIZE * 2 + 32]);

                        above_row = above_data + 16;
                        left_col = left_data + 16;
                        uint8_t *recon_buffer =
                            encode_context_ptr->mc_flow_rec_picture_buffer[frame_idx] +
                            dst_basic_offset;
                        update_neighbor_samples_array_open_loop_mb_recon(above_row - 1,
                                                                         left_col - 1,
                                                                         recon_buffer,
                                                                         dst_buffer_stride,
                                                                         mb_origin_x,
                                                                         mb_origin_y,
                                                                         16,
                                                                         16,
                                                                         input_picture_ptr->width,
                                                                         input_picture_ptr->height);
                        uint8_t ois_intra_mode = ois_mb_results_ptr->intra_mode;
                        int32_t p_angle = av1_is_directional_mode((PredictionMode)ois_intra_mode) ? mode_to_angle_map[(PredictionMode)ois_intra_mode] : 0;
                        // Edge filter
                        if(av1_is_directional_mode((PredictionMode)ois_intra_mode) && 1/*scs_ptr->seq_header.enable_intra_edge_filter*/) {
                            filter_intra_edge(ois_mb_results_ptr, ois_intra_mode, scs_ptr->seq_header.max_frame_width, scs_ptr->seq_header.max_frame_height, p_angle, mb_origin_x, mb_origin_y, above_row, left_col);
                        }
                        // PRED
                        intra_prediction_open_loop_mb(p_angle,
                                                      ois_intra_mode,
                                                      mb_origin_x,
                                                      mb_origin_y,
                                                      TX_16X16,
                                                      above_row,
                                                      left_col,
                                                      dst_buffer,
                                                      dst_buffer_stride);
                    }

                    svt_aom_subtract_block(16, 16, src_diff, 16, src_mb, input_picture_ptr->stride_y, dst_buffer, dst_buffer_stride);
                    svt_av1_wht_fwd_txfm(src_diff, 16, coeff, tx_size, 8, 0);

                    uint16_t eob = 0;

                    get_quantize_error(&mb_plane, coeff, qcoeff, dqcoeff, tx_size, &eob, &recon_error, &sse);
                    int rate_cost = pcs_ptr->tpl_opt_flag ? 0 : rate_estimator(qcoeff, eob, tx_size);

                    if(eob) {
                        av1_inv_transform_recon8bit((int32_t*)dqcoeff, dst_buffer, dst_buffer_stride, dst_buffer, dst_buffer_stride, TX_16X16, DCT_DCT, PLANE_TYPE_Y, eob, 0);
                    }

                    tpl_stats.recrf_dist = recon_error << (TPL_DEP_COST_SCALE_LOG2);
                    tpl_stats.recrf_rate = rate_cost << TPL_DEP_COST_SCALE_LOG2;
                    if (best_mode != NEWMV) {
                        tpl_stats.srcrf_dist = recon_error << (TPL_DEP_COST_SCALE_LOG2);
                        tpl_stats.srcrf_rate = rate_cost << TPL_DEP_COST_SCALE_LOG2;
                    }
                    tpl_stats.recrf_dist = AOMMAX(tpl_stats.srcrf_dist, tpl_stats.recrf_dist);
                    tpl_stats.recrf_rate = AOMMAX(tpl_stats.srcrf_rate, tpl_stats.recrf_rate);
                    if (!frame_is_intra_only(pcs_ptr) && best_rf_idx != -1) {
                        tpl_stats.mv = final_best_mv;
                        tpl_stats.ref_frame_poc = best_ref_poc;
                    }
                    // Motion flow dependency dispenser.
                    result_model_store(pcs_ptr, &tpl_stats, mb_origin_x, mb_origin_y);
                }
                pa_blk_index++;
            }
        }
    }
    // padding current recon picture
    generate_padding(
        encode_context_ptr->mc_flow_rec_picture_buffer[frame_idx],
        input_picture_ptr->stride_y,
        input_picture_ptr->width,
        input_picture_ptr->height,
        input_picture_ptr->origin_x,
        input_picture_ptr->origin_y);

    return;
}

static int get_overlap_area(int grid_pos_row, int grid_pos_col, int ref_pos_row,
                            int ref_pos_col, int block, int/*BLOCK_SIZE*/ bsize) {
  int width = 0, height = 0;
  int bw = 4 << mi_size_wide_log2[bsize];
  int bh = 4 << mi_size_high_log2[bsize];

  switch (block) {
    case 0:
      width = grid_pos_col + bw - ref_pos_col;
      height = grid_pos_row + bh - ref_pos_row;
      break;
    case 1:
      width = ref_pos_col + bw - grid_pos_col;
      height = grid_pos_row + bh - ref_pos_row;
      break;
    case 2:
      width = grid_pos_col + bw - ref_pos_col;
      height = ref_pos_row + bh - grid_pos_row;
      break;
    case 3:
      width = ref_pos_col + bw - grid_pos_col;
      height = ref_pos_row + bh - grid_pos_row;
      break;
    default: assert(0);
  }

  return width * height;
}

static int round_floor(int ref_pos, int bsize_pix) {
  int round;
  if (ref_pos < 0)
    round = -(1 + (-ref_pos - 1) / bsize_pix);
  else
    round = ref_pos / bsize_pix;

  return round;
}

static int64_t delta_rate_cost(int64_t delta_rate, int64_t recrf_dist,
                               int64_t srcrf_dist, int pix_num) {
  double beta = (double)srcrf_dist / recrf_dist;
  int64_t rate_cost = delta_rate;

  if (srcrf_dist <= 128) return rate_cost;

  double dr =
      (double)(delta_rate >> (TPL_DEP_COST_SCALE_LOG2 + AV1_PROB_COST_SHIFT)) /
      pix_num;

  double log_den = log(beta) / log(2.0) + 2.0 * dr;

  if (log_den > log(10.0) / log(2.0)) {
    rate_cost = (int64_t)((log(1.0 / beta) * pix_num) / log(2.0) / 2.0);
    rate_cost <<= (TPL_DEP_COST_SCALE_LOG2 + AV1_PROB_COST_SHIFT);
    return rate_cost;
  }

  double num = pow(2.0, log_den);
  double den = num * beta + (1 - beta) * beta;

  rate_cost = (int64_t)((pix_num * log(num / den)) / log(2.0) / 2.0);

  rate_cost <<= (TPL_DEP_COST_SCALE_LOG2 + AV1_PROB_COST_SHIFT);

  return rate_cost;
}

static AOM_INLINE void tpl_model_update_b(PictureParentControlSet *ref_pcs_ptr, PictureParentControlSet *pcs_ptr,
                                          TplStats *tpl_stats_ptr,
                                          int mi_row, int mi_col,
                                          const int/*BLOCK_SIZE*/ bsize) {
  Av1Common *ref_cm = ref_pcs_ptr->av1_cm;
  TplStats *ref_tpl_stats_ptr;

  const FULLPEL_MV full_mv = get_fullmv_from_mv(&tpl_stats_ptr->mv);
  const int ref_pos_row = mi_row * MI_SIZE + full_mv.row;
  const int ref_pos_col = mi_col * MI_SIZE + full_mv.col;

  const int bw = 4 << mi_size_wide_log2[bsize];
  const int bh = 4 << mi_size_high_log2[bsize];
  const int mi_height = mi_size_high[bsize];
  const int mi_width = mi_size_wide[bsize];
  const int pix_num = bw * bh;
  const int shift = pcs_ptr->is_720p_or_larger ? 2 : 1;
  const int mi_cols_sr = ((ref_pcs_ptr->aligned_width + 15) / 16) << 2;

  // top-left on grid block location in pixel
  int grid_pos_row_base = round_floor(ref_pos_row, bh) * bh;
  int grid_pos_col_base = round_floor(ref_pos_col, bw) * bw;
  int block;

  int64_t cur_dep_dist = tpl_stats_ptr->recrf_dist - tpl_stats_ptr->srcrf_dist;
  int64_t mc_dep_dist = (int64_t)(
      tpl_stats_ptr->mc_dep_dist *
      ((double)(tpl_stats_ptr->recrf_dist - tpl_stats_ptr->srcrf_dist) /
       tpl_stats_ptr->recrf_dist));
  int64_t delta_rate = tpl_stats_ptr->recrf_rate - tpl_stats_ptr->srcrf_rate;
  int64_t mc_dep_rate =
      delta_rate_cost(tpl_stats_ptr->mc_dep_rate, tpl_stats_ptr->recrf_dist,
                      tpl_stats_ptr->srcrf_dist, pix_num);

  for (block = 0; block < 4; ++block) {
    int grid_pos_row = grid_pos_row_base + bh * (block >> 1);
    int grid_pos_col = grid_pos_col_base + bw * (block & 0x01);

    if (grid_pos_row >= 0 && grid_pos_row < ref_cm->mi_rows * MI_SIZE &&
        grid_pos_col >= 0 && grid_pos_col < ref_cm->mi_cols * MI_SIZE) {
      int overlap_area = get_overlap_area(
          grid_pos_row, grid_pos_col, ref_pos_row, ref_pos_col, block, bsize);
      int ref_mi_row = round_floor(grid_pos_row, bh) * mi_height;
      int ref_mi_col = round_floor(grid_pos_col, bw) * mi_width;
      const int step = 1 << (pcs_ptr->is_720p_or_larger ? 2 : 1);

      for (int idy = 0; idy < mi_height; idy += step) {
        for (int idx = 0; idx < mi_width; idx += step) {
          ref_tpl_stats_ptr = ref_pcs_ptr->tpl_stats[((ref_mi_row + idy) >> shift) * (mi_cols_sr >> shift)  + ((ref_mi_col + idx) >> shift)];
          ref_tpl_stats_ptr->mc_dep_dist +=
              ((cur_dep_dist + mc_dep_dist) * overlap_area) / pix_num;
          ref_tpl_stats_ptr->mc_dep_rate +=
              ((delta_rate + mc_dep_rate) * overlap_area) / pix_num;
          assert(overlap_area >= 0);
        }
      }
    }
  }
}

static AOM_INLINE void tpl_model_update(PictureParentControlSet *pcs_array[MAX_TPL_LA_SW], int32_t frame_idx, int mi_row, int mi_col, const int/*BLOCK_SIZE*/ bsize, uint8_t frames_in_sw) {
  const int mi_height = mi_size_high[bsize];
  const int mi_width = mi_size_wide[bsize];
  PictureParentControlSet *pcs_ptr = pcs_array[frame_idx];
  const int/*BLOCK_SIZE*/ block_size = pcs_ptr->is_720p_or_larger ? BLOCK_16X16 : BLOCK_8X8;
  const int step = 1 << (pcs_ptr->is_720p_or_larger ? 2 : 1);
  const int shift = pcs_ptr->is_720p_or_larger ? 2 : 1;
  const int mi_cols_sr = ((pcs_ptr->aligned_width + 15) / 16) << 2;
  int i = 0;

  for (int idy = 0; idy < mi_height; idy += step) {
    for (int idx = 0; idx < mi_width; idx += step) {
      TplStats *tpl_stats_ptr = pcs_ptr->tpl_stats[(((mi_row + idy) >> shift) * (mi_cols_sr >> shift)) + ((mi_col + idx) >> shift)];

      while(i<frames_in_sw && pcs_array[i]->picture_number != tpl_stats_ptr->ref_frame_poc)
        i++;
      if(i<frames_in_sw)
        tpl_model_update_b(pcs_array[i], pcs_ptr, tpl_stats_ptr, mi_row + idy, mi_col + idx, block_size);
    }
  }
}

/************************************************
* Genrate TPL MC Flow Synthesizer Based on Lookahead
** LAD Window: sliding window size
************************************************/
void tpl_mc_flow_synthesizer(
    PictureParentControlSet         *pcs_array[MAX_TPL_LA_SW],
    int32_t                          frame_idx,
    uint8_t                          frames_in_sw)
{
    Av1Common *cm = pcs_array[frame_idx]->av1_cm;
    const int/*BLOCK_SIZE*/ bsize = BLOCK_16X16;
    const int mi_height = mi_size_high[bsize];
    const int mi_width = mi_size_wide[bsize];

    for (int mi_row = 0; mi_row < cm->mi_rows; mi_row += mi_height) {
        for (int mi_col = 0; mi_col < cm->mi_cols; mi_col += mi_width) {
            tpl_model_update(pcs_array, frame_idx, mi_row, mi_col, bsize, frames_in_sw);
        }
    }
    return;
}

static void generate_r0beta(PictureParentControlSet *pcs_ptr)
{
    Av1Common *cm = pcs_ptr->av1_cm;
    SequenceControlSet *scs_ptr = pcs_ptr->scs_ptr;
    int64_t intra_cost_base = 0;
    int64_t mc_dep_cost_base = 0;
    const int step = 1 << (pcs_ptr->is_720p_or_larger ? 2 : 1);
    const int mi_cols_sr = ((pcs_ptr->aligned_width + 15) / 16) << 2;
    const int shift = pcs_ptr->is_720p_or_larger ? 2 : 1;

    for (int row = 0; row < cm->mi_rows; row += step) {
        for (int col = 0; col < mi_cols_sr; col += step) {
            TplStats *tpl_stats_ptr = pcs_ptr->tpl_stats[(row >> shift) * (mi_cols_sr >> shift) + (col >> shift)];
            int64_t mc_dep_delta =
                RDCOST(pcs_ptr->base_rdmult, tpl_stats_ptr->mc_dep_rate, tpl_stats_ptr->mc_dep_dist);
            intra_cost_base  += (tpl_stats_ptr->recrf_dist << RDDIV_BITS);
            mc_dep_cost_base += (tpl_stats_ptr->recrf_dist << RDDIV_BITS) + mc_dep_delta;
        }
    }

    if (mc_dep_cost_base != 0) {
      pcs_ptr->r0 = (double)intra_cost_base / mc_dep_cost_base;
    }

    //SVT_LOG("generate_r0beta ------> poc %ld\t%.0f\t%.0f \t%.5f base_rdmult=%d\n", pcs_ptr->picture_number, (double)intra_cost_base, (double)mc_dep_cost_base, pcs_ptr->r0, pcs_ptr->base_rdmult);
    generate_lambda_scaling_factor(pcs_ptr, mc_dep_cost_base);

    const uint32_t sb_sz = scs_ptr->seq_header.sb_size == BLOCK_128X128 ? 128 : 64;
    const uint32_t picture_sb_width  = (uint32_t)((scs_ptr->seq_header.max_frame_width  + sb_sz - 1) / sb_sz);
    const uint32_t picture_sb_height = (uint32_t)((scs_ptr->seq_header.max_frame_height + sb_sz - 1) / sb_sz);
    const uint32_t picture_width_in_mb  = (scs_ptr->seq_header.max_frame_width  + 16 - 1) / 16;
    const uint32_t picture_height_in_mb = (scs_ptr->seq_header.max_frame_height + 16 - 1) / 16;
    const uint32_t blks = scs_ptr->seq_header.sb_size == BLOCK_128X128 ? (128 >> (3 + pcs_ptr->is_720p_or_larger))
                                                                       : (64  >> (3 + pcs_ptr->is_720p_or_larger));
    for (uint32_t sb_y = 0; sb_y < picture_sb_height; ++sb_y) {
        for (uint32_t sb_x = 0; sb_x < picture_sb_width; ++sb_x) {
            int64_t intra_cost = 0;
            int64_t mc_dep_cost = 0;
            for (uint32_t blky_offset = 0; blky_offset < blks; blky_offset++) {
                for (uint32_t blkx_offset = 0; blkx_offset < blks; blkx_offset++) {
                    uint32_t blkx = ((sb_x * sb_sz) >> (3 + pcs_ptr->is_720p_or_larger)) + blkx_offset;
                    uint32_t blky = ((sb_y * sb_sz) >> (3 + pcs_ptr->is_720p_or_larger)) + blky_offset;
                    if((blkx >> (1 - pcs_ptr->is_720p_or_larger)) >= picture_width_in_mb ||
                       (blky >> (1 - pcs_ptr->is_720p_or_larger)) >= picture_height_in_mb)
                        continue;
                    TplStats *tpl_stats_ptr = pcs_ptr->tpl_stats[blky * (mi_cols_sr >> shift) + blkx];
                    int64_t mc_dep_delta =
                        RDCOST(pcs_ptr->base_rdmult, tpl_stats_ptr->mc_dep_rate, tpl_stats_ptr->mc_dep_dist);
                    intra_cost  += (tpl_stats_ptr->recrf_dist << RDDIV_BITS);
                    mc_dep_cost += (tpl_stats_ptr->recrf_dist << RDDIV_BITS) + mc_dep_delta;
                }
            }
            double beta = 1.0;
            if (mc_dep_cost > 0 && intra_cost > 0) {
            //if (mc_dep_cost > ((scs_ptr->seq_header.sb_size == BLOCK_128X128 ? 8 : 4)*
            //    (scs_ptr->seq_header.sb_size == BLOCK_128X128 ? 8 : 4) << RDDIV_BITS)
            //    && intra_cost > ((scs_ptr->seq_header.sb_size == BLOCK_128X128 ? 8 : 4)*
            //    (scs_ptr->seq_header.sb_size == BLOCK_128X128 ? 8 : 4) << RDDIV_BITS) ){
                double rk = (double)intra_cost / mc_dep_cost;
                beta = (pcs_ptr->r0 / rk);
                //if (pcs_ptr->picture_number == 16)
                //SVT_LOG(
                //        "generate_r0beta ---> poc%ld sb_x=%d sb_y=%d r0=%f rk=%f\t intraC=%lld "
                //        "mc_dep=%lld\tbeta=%f\n",
                //        pcs_ptr->picture_number,
                //        sb_x,
                //        sb_y,
                //        pcs_ptr->r0,
                //        rk,
                //        intra_cost,
                //        mc_dep_cost,
                //        beta);
                assert(beta > 0.0);
            }
            //SVT_LOG("generate_r0beta sbxy=%d %d, rk=%f beta=%f\n", sb_x, sb_y, rk, beta);
            pcs_ptr->tpl_beta[sb_y * picture_sb_width + sb_x] = beta;
        }
    }
    return;
}

/************************************************
* Genrate TPL MC Flow Based on Lookahead
** LAD Window: sliding window size
************************************************/
EbErrorType tpl_mc_flow(
    EncodeContext                   *encode_context_ptr,
    SequenceControlSet              *scs_ptr,
    PictureParentControlSet         *pcs_ptr)
{
    PictureParentControlSet          *pcs_array[MAX_TPL_LA_SW] = {NULL, };

    uint32_t                         inputQueueIndex;
    int32_t                          frames_in_sw = MIN(MAX_TPL_LA_SW, pcs_ptr->frames_in_sw);
    uint32_t                         shift = pcs_ptr->is_720p_or_larger ? 0 : 1;
    uint32_t picture_width_in_mb  = (pcs_ptr->enhanced_picture_ptr->width  + 16 - 1) / 16;
    uint32_t picture_height_in_mb = (pcs_ptr->enhanced_picture_ptr->height + 16 - 1) / 16;
    EbBool                           start_is_intra = EB_FALSE;
    EbBool                           got_intra_in_sw = EB_FALSE;
    EbByte                           mc_flow_rec_picture_buffer_noref = NULL;

    (void)scs_ptr;
    pcs_ptr->r0 = 0;
    pcs_array[0] = pcs_ptr;
    // Walk the first N entries in the sliding window
    inputQueueIndex = encode_context_ptr->initial_rate_control_reorder_queue_head_index;
    for (int32_t frame_idx = 0; frame_idx < pcs_ptr->frames_in_sw; frame_idx++) {
        InitialRateControlReorderEntry   *temporaryQueueEntryPtr = encode_context_ptr->initial_rate_control_reorder_queue[inputQueueIndex];
        PictureParentControlSet          *temp_pcs_ptr = ((PictureParentControlSet*)(temporaryQueueEntryPtr->parent_pcs_wrapper_ptr)->object_ptr);

        // sort to be decode order
        if(frame_idx == 0) {
            pcs_array[0] = temp_pcs_ptr;
        } else {
            int32_t i = 0;
            for (; i < frame_idx; i++) {
                if (temp_pcs_ptr->decode_order < pcs_array[i]->decode_order) {
                    for (int32_t j = frame_idx; j > i; j--)
                        pcs_array[j] = pcs_array[j-1];
                    pcs_array[i] = temp_pcs_ptr;
                    break;
                }
            }
            if (i == frame_idx)
                pcs_array[i] = temp_pcs_ptr;
        }

        // Increment the inputQueueIndex Iterator
        inputQueueIndex = (inputQueueIndex == INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH - 1) ? 0 : inputQueueIndex + 1;
    }

    for (int32_t frame_idx = 1; frame_idx < MIN(16, frames_in_sw); frame_idx++) {
        if(frame_is_intra_only(pcs_array[frame_idx])) {
            got_intra_in_sw = EB_TRUE;
            break;
        }
    }
    if (got_intra_in_sw) {
        return EB_ErrorNone;
    }
    start_is_intra = frame_is_intra_only(pcs_array[0]);

    for(int32_t frame_idx = 0; frame_idx < MAX_TPL_LA_SW; frame_idx++) {
        encode_context_ptr->poc_map_idx[frame_idx] = -1;
        encode_context_ptr->mc_flow_rec_picture_buffer[frame_idx] = NULL;
    }
    EB_MALLOC_ARRAY(mc_flow_rec_picture_buffer_noref, pcs_ptr->enhanced_picture_ptr->luma_size);
    for(int32_t frame_idx = 0; frame_idx < frames_in_sw; frame_idx++) {
        if (pcs_array[frame_idx]->is_used_as_reference_flag) {
            EB_MALLOC_ARRAY(encode_context_ptr->mc_flow_rec_picture_buffer[frame_idx], pcs_ptr->enhanced_picture_ptr->luma_size);
        } else {
            encode_context_ptr->mc_flow_rec_picture_buffer[frame_idx] = mc_flow_rec_picture_buffer_noref;
        }
    }
    if (!encode_context_ptr->mc_flow_rec_picture_buffer_saved)
        EB_MALLOC_ARRAY(encode_context_ptr->mc_flow_rec_picture_buffer_saved, pcs_ptr->enhanced_picture_ptr->luma_size);

    if (pcs_array[0]->temporal_layer_index == 0) {
        // dispenser I0 or frame_idx0 pic in LA1
        int32_t sw_length = MIN(17, (frames_in_sw));
        EbPictureBufferDesc *input_picture_ptr = pcs_array[0]->enhanced_picture_ptr;
        encode_context_ptr->poc_map_idx[0] = pcs_array[0]->picture_number;
        for(int32_t frame_idx = 0; frame_idx < sw_length; frame_idx++) {
            if (!start_is_intra && frame_idx == 0) {
                uint8_t *dst_buffer = encode_context_ptr->mc_flow_rec_picture_buffer[0];
                memcpy(dst_buffer, encode_context_ptr->mc_flow_rec_picture_buffer_saved, input_picture_ptr->stride_y * (input_picture_ptr->origin_y * 2 + input_picture_ptr->height));
                continue;
            }
            encode_context_ptr->poc_map_idx[frame_idx] = pcs_array[frame_idx]->picture_number;
            for (uint32_t blky = 0; blky < (picture_height_in_mb << shift); blky++) {
                memset(pcs_array[frame_idx]->tpl_stats[blky * (picture_width_in_mb << shift)], 0, (picture_width_in_mb << shift) * sizeof(TplStats));
            }

            tpl_mc_flow_dispenser(encode_context_ptr, scs_ptr, pcs_array[frame_idx], frame_idx);

        }

        // synthesizer
        for(int32_t frame_idx = sw_length - 1; frame_idx >= 0; frame_idx--)
            tpl_mc_flow_synthesizer(pcs_array, frame_idx, sw_length);
        generate_r0beta(pcs_array[0]);

        // The second part is for the next base layer frame to use the available pictures.
        // i.e. POC 16 have access to picture 1,2,...15. So dispenser and synthesizer are called.
        // In the next call, the stats for POC 16 is updated using pictures 17,... 32
        encode_context_ptr->poc_map_idx[0] = pcs_array[0]->picture_number;
        for (int32_t frame_idx = 1; frame_idx < sw_length; frame_idx++) {
            encode_context_ptr->poc_map_idx[frame_idx] = pcs_array[frame_idx]->picture_number;
            if (frame_idx == 1 && pcs_array[frame_idx]->temporal_layer_index == 0) {
                // copy frame_idx0 input to rec before run frame_idx1 picture
                memcpy(encode_context_ptr->mc_flow_rec_picture_buffer[0],
                       pcs_array[0]->enhanced_picture_ptr->buffer_y,
                       pcs_array[0]->enhanced_picture_ptr->stride_y *
                           (pcs_array[0]->enhanced_picture_ptr->origin_y * 2 +
                            pcs_array[0]->enhanced_picture_ptr->height));
            }
            for (uint32_t blky = 0; blky < (picture_height_in_mb << shift); blky++) {
                memset(pcs_array[frame_idx]->tpl_stats[blky * (picture_width_in_mb << shift)], 0, (picture_width_in_mb << shift) * sizeof(TplStats));
            }

            tpl_mc_flow_dispenser(encode_context_ptr, scs_ptr, pcs_array[frame_idx], frame_idx);
            if (frame_idx == 1 && pcs_array[frame_idx]->temporal_layer_index == 0) {
                // save frame_idx1 picture buffer for next LA
                memcpy(encode_context_ptr->mc_flow_rec_picture_buffer_saved, encode_context_ptr->mc_flow_rec_picture_buffer[frame_idx], input_picture_ptr->stride_y * (input_picture_ptr->origin_y * 2 + input_picture_ptr->height));
            }
        }
        // synthesizer
        for (int32_t frame_idx = sw_length - 1; frame_idx >= 1; frame_idx--)
            // to make sure synthesizer is not called more than one time
            if(pcs_array[frame_idx]->picture_number <= pcs_array[1]->picture_number)
            tpl_mc_flow_synthesizer(pcs_array, frame_idx, sw_length);
    }

    for(int32_t frame_idx = 0; frame_idx < frames_in_sw; frame_idx++) {
        if ( encode_context_ptr->mc_flow_rec_picture_buffer[frame_idx] &&
             encode_context_ptr->mc_flow_rec_picture_buffer[frame_idx] != mc_flow_rec_picture_buffer_noref ) {
            EB_FREE_ARRAY(encode_context_ptr->mc_flow_rec_picture_buffer[frame_idx]);
        }
        encode_context_ptr->mc_flow_rec_picture_buffer[frame_idx] = NULL;
    }
    EB_FREE_ARRAY(mc_flow_rec_picture_buffer_noref);

    return EB_ErrorNone;
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

    EbObjectWrapper *        in_results_wrapper_ptr;

    EbObjectWrapper *          out_results_wrapper_ptr;

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

        uint32_t segment_index = in_results_ptr->segment_index;

        // Set the segment mask
        SEGMENT_COMPLETION_MASK_SET(pcs_ptr->me_segments_completion_mask, segment_index);

        // If the picture is complete, proceed
        if (SEGMENT_COMPLETION_MASK_TEST(pcs_ptr->me_segments_completion_mask,
                                         pcs_ptr->me_segments_total_count)) {
            SequenceControlSet *scs_ptr = (SequenceControlSet *)
                                              pcs_ptr->scs_wrapper_ptr->object_ptr;
            EncodeContext *encode_context_ptr = (EncodeContext *)scs_ptr->encode_context_ptr;
            if (scs_ptr->static_config.look_ahead_distance == 0 || scs_ptr->static_config.enable_tpl_la == 0) {
                // Release Pa Ref pictures when not needed
                release_pa_reference_objects(scs_ptr, pcs_ptr);
            }
            /*In case Look-Ahead is zero there is no need to place pictures in the
              re-order queue. this will cause an artificial delay since pictures come in dec-order*/
            if (scs_ptr->static_config.look_ahead_distance == 0) {


                for (uint8_t temporal_layer_index = 0; temporal_layer_index < EB_MAX_TEMPORAL_LAYERS;
                    temporal_layer_index++)
                    pcs_ptr->frames_in_interval[temporal_layer_index] = 0;

                pcs_ptr->frames_in_sw = 0;
                pcs_ptr->historgram_life_count = 0;
                pcs_ptr->scene_change_in_gop = EB_FALSE;
                pcs_ptr->end_of_sequence_region = EB_FALSE;

                init_zz_cost_info(pcs_ptr);

                // Get Empty Results Object
                svt_get_empty_object(
                    context_ptr->initialrate_control_results_output_fifo_ptr,
                    &out_results_wrapper_ptr);

                InitialRateControlResults * out_results_ptr =
                    (InitialRateControlResults *)out_results_wrapper_ptr->object_ptr;

                out_results_ptr->pcs_wrapper_ptr = pcs_ptr->p_pcs_wrapper_ptr;
                svt_post_full_object(out_results_wrapper_ptr);

            }
            else {
            //****************************************************
            // Input Motion Analysis Results into Reordering Queue
            //****************************************************

            if (!pcs_ptr->is_overlay)
                // Determine offset from the Head Ptr
                determine_picture_offset_in_queue(
                    encode_context_ptr, pcs_ptr, in_results_ptr);

            if (use_input_stat(scs_ptr) && scs_ptr->static_config.rate_control_mode == 1)
                ; //skip 2pass VBR
            else
            if (scs_ptr->static_config.rate_control_mode) {
                // Getting the Histogram Queue Data
                get_histogram_queue_data(scs_ptr, encode_context_ptr, pcs_ptr);
            }

            for (uint8_t temporal_layer_index = 0; temporal_layer_index < EB_MAX_TEMPORAL_LAYERS;
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
                    encode_context_ptr->initial_rate_control_reorder_queue[queue_entry_index_temp]
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
                        (queue_entry_index_temp > INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH - 1)
                        ? queue_entry_index_temp - INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH
                        : queue_entry_index_temp;

                    move_slide_window_flag =
                        (EbBool)(move_slide_window_flag &&
                                 (encode_context_ptr
                                      ->initial_rate_control_reorder_queue[queue_entry_index_temp2]
                                      ->parent_pcs_wrapper_ptr != NULL));
                    if (encode_context_ptr
                            ->initial_rate_control_reorder_queue[queue_entry_index_temp2]
                            ->parent_pcs_wrapper_ptr != NULL) {
                        // check if it is the last frame. If we have reached the last frame, we would output the buffered frames in the Queue.
                        end_of_sequence_flag =
                            ((PictureParentControlSet *)(encode_context_ptr
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
                    pcs_ptr = ((PictureParentControlSet *)(queue_entry_ptr->parent_pcs_wrapper_ptr)
                                   ->object_ptr);
                    scs_ptr = (SequenceControlSet *)pcs_ptr->scs_wrapper_ptr->object_ptr;
                    // overlay picture was not added to the queue. For the alt_ref picture with an overlay picture, it loops on both alt ref and overlay pictures
                    uint8_t has_overlay = pcs_ptr->is_alt_ref ? 1 : 0;
                    for (uint8_t loop_index = 0; loop_index <= has_overlay; loop_index++) {
                        if (loop_index) pcs_ptr = pcs_ptr->overlay_ppcs_ptr;
                        pcs_ptr->frames_in_sw = frames_in_sw;
                        queue_entry_index_temp =
                            encode_context_ptr->initial_rate_control_reorder_queue_head_index;
                        end_of_sequence_flag = EB_FALSE;
                        // find the frames_in_interval for the peroid I frames
                        while (
                            !end_of_sequence_flag &&
                            queue_entry_index_temp <=
                                encode_context_ptr->initial_rate_control_reorder_queue_head_index +
                                    scs_ptr->static_config.look_ahead_distance) {
                            uint32_t queue_entry_index_temp2 =
                                (queue_entry_index_temp >
                                 INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH - 1)
                                ? queue_entry_index_temp -
                                    INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH
                                : queue_entry_index_temp;
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
                                    pcs_ptr
                                        ->frames_in_interval[pcs_ptr_temp->temporal_layer_index]++;
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

                        if (use_input_stat(scs_ptr) && scs_ptr->static_config.rate_control_mode == 1)
                            ; //skip 2pass VBR
                        else
                        if (scs_ptr->static_config.rate_control_mode) {
                            // Determine offset from the Head Ptr for HLRC histogram queue and set the life count
                            if (scs_ptr->static_config.look_ahead_distance != 0) {
                                // Update Histogram Queue Entry Life count
                                update_histogram_queue_entry(
                                    scs_ptr, encode_context_ptr, pcs_ptr, frames_in_sw);
                            }
                        }
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
                        // Get Empty Reference Picture Object
                        svt_get_empty_object(
                            scs_ptr->encode_context_ptr->reference_picture_pool_fifo_ptr,
                            &reference_picture_wrapper_ptr);
                        if (loop_index) {
                            pcs_ptr->reference_picture_wrapper_ptr = reference_picture_wrapper_ptr;
                            // Give the new Reference a nominal live_count of 1
                            svt_object_inc_live_count(pcs_ptr->reference_picture_wrapper_ptr, 1);
                        } else {
                            ((PictureParentControlSet *)(queue_entry_ptr->parent_pcs_wrapper_ptr
                                                             ->object_ptr))
                                ->reference_picture_wrapper_ptr = reference_picture_wrapper_ptr;
                            // Give the new Reference a nominal live_count of 1
                            svt_object_inc_live_count(
                                ((PictureParentControlSet *)(queue_entry_ptr->parent_pcs_wrapper_ptr
                                                                 ->object_ptr))
                                    ->reference_picture_wrapper_ptr,
                                1);
                        }
                        if (scs_ptr->static_config.look_ahead_distance != 0 &&
                            scs_ptr->static_config.enable_tpl_la &&
                            pcs_ptr->temporal_layer_index == 0) {
                            tpl_mc_flow(encode_context_ptr, scs_ptr, pcs_ptr);
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
                        if (scs_ptr->static_config.look_ahead_distance != 0 && scs_ptr->static_config.enable_tpl_la
                            && ((has_overlay == 0 && loop_index == 0) || (has_overlay == 1 && loop_index == 1))) {
                            // Release Pa Ref pictures when not needed
                            release_pa_reference_objects(scs_ptr, pcs_ptr);
                                //loop_index ? pcs_ptr : queueEntryPtr);
                        }
                        // Post the Full Results Object
                        svt_post_full_object(out_results_wrapper_ptr);
                    }
                    // Reset the Reorder Queue Entry
                    queue_entry_ptr->picture_number += INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH;
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

        // Release the Input Results
        svt_release_object(in_results_wrapper_ptr);
    }
    return NULL;
}
