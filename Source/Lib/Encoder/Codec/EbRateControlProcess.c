/*
* Copyright(c) 2019 Intel Corporation
* Copyright (c) 2016, Alliance for Open Media. All rights reserved
*
* This source code is subject to the terms of the BSD 3-Clause Clear License and
* the Alliance for Open Media Patent License 1.0. If the BSD 3-Clause Clear License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/
#include <stdlib.h>

#include "EbDefinitions.h"
#include "EbEncHandle.h"
#include "EbRateControlProcess.h"
#include "EbSequenceControlSet.h"
#include "EbPictureControlSet.h"
#include "EbUtility.h"
#include "EbSvtAv1ErrorCodes.h"
#include "EbEntropyCoding.h"

#include "EbRateControlResults.h"
#include "EbRateControlTasks.h"

#include "EbSegmentation.h"
#include "EbLog.h"
#include "EbRateDistortionCost.h"
#include "EbLambdaRateTables.h"
#include "pass2_strategy.h"

#include "EbTransforms.h"
#include "aom_dsp_rtcd.h"
#include "EbLog.h"
#include "EbIntraPrediction.h"
#include "EbMotionEstimation.h"

#include "EbPictureDecisionResults.h"
#include "EbResize.h"

static const double tpl_hl_islice_div_factor[EB_MAX_TEMPORAL_LAYERS]     = {1, 1, 1, 2, 1, 0.7};
static const double tpl_hl_base_frame_div_factor[EB_MAX_TEMPORAL_LAYERS] = {1, 1, 1, 3, 1, 0.9};
#define KB 400

typedef struct RateControlContext {
    EbFifo *rate_control_input_tasks_fifo_ptr;
    EbFifo *rate_control_output_results_fifo_ptr;
    EbFifo *picture_decision_results_output_fifo_ptr;

    RateControlIntervalParamContext **rate_control_param_queue;
} RateControlContext;
EbErrorType rate_control_coded_frames_stats_context_ctor(coded_frames_stats_entry *entry_ptr,
                                                         uint64_t                  picture_number) {
    entry_ptr->picture_number         = picture_number;
    entry_ptr->frame_total_bit_actual = -1;

    return EB_ErrorNone;
}
static void rate_control_context_dctor(EbPtr p) {
    EbThreadContext    *thread_context_ptr = (EbThreadContext *)p;
    RateControlContext *obj                = (RateControlContext *)thread_context_ptr->priv;
    if (obj->rate_control_param_queue)
        EB_FREE_2D(obj->rate_control_param_queue);
    EB_FREE_ARRAY(obj);
}

EbErrorType rate_control_context_ctor(EbThreadContext   *thread_context_ptr,
                                      const EbEncHandle *enc_handle_ptr, int me_port_index) {
    uint32_t interval_index;

    int32_t intra_period =
        enc_handle_ptr->scs_instance_array[0]->scs_ptr->static_config.intra_period_length;

    RateControlContext *context_ptr;
    EB_CALLOC_ARRAY(context_ptr, 1);
    thread_context_ptr->priv  = context_ptr;
    thread_context_ptr->dctor = rate_control_context_dctor;

    context_ptr->rate_control_input_tasks_fifo_ptr = svt_system_resource_get_consumer_fifo(
        enc_handle_ptr->rate_control_tasks_resource_ptr, 0);
    context_ptr->rate_control_output_results_fifo_ptr = svt_system_resource_get_producer_fifo(
        enc_handle_ptr->rate_control_results_resource_ptr, 0);
    context_ptr->picture_decision_results_output_fifo_ptr = svt_system_resource_get_producer_fifo(
        enc_handle_ptr->picture_decision_results_resource_ptr, me_port_index);

    EB_MALLOC_2D(context_ptr->rate_control_param_queue, (int32_t)PARALLEL_GOP_MAX_NUMBER, 1);

    for (interval_index = 0; interval_index < PARALLEL_GOP_MAX_NUMBER; interval_index++) {
        context_ptr->rate_control_param_queue[interval_index]->first_poc =
            (interval_index * (uint32_t)(intra_period + 1));
        context_ptr->rate_control_param_queue[interval_index]->last_poc =
            ((interval_index + 1) * (uint32_t)(intra_period + 1)) - 1;
        context_ptr->rate_control_param_queue[interval_index]->processed_frame_number = 0;
        context_ptr->rate_control_param_queue[interval_index]->end_of_seq_seen        = 0;
        context_ptr->rate_control_param_queue[interval_index]->last_i_qp              = 0;
    }

    return EB_ErrorNone;
}

#define MAX_Q_INDEX 255
#define MIN_Q_INDEX 0

extern int16_t svt_av1_ac_quant_q3(int32_t qindex, int32_t delta, AomBitDepth bit_depth);
// These functions use formulaic calculations to make playing with the
// quantizer tables easier. If necessary they can be replaced by lookup
// tables if and when things settle down in the experimental Bitstream

int32_t svt_av1_convert_qindex_to_q_fp8(int32_t qindex, AomBitDepth bit_depth) {
    // Convert the index to a real Q value (scaled down to match old Q values)
    switch (bit_depth) {
    case AOM_BITS_8: return svt_av1_ac_quant_q3(qindex, 0, bit_depth) << 6; // / 4.0;
    case AOM_BITS_10: return svt_av1_ac_quant_q3(qindex, 0, bit_depth) << 4; // / 16.0;
    case AOM_BITS_12: return svt_av1_ac_quant_q3(qindex, 0, bit_depth) << 3; // / 64.0;
    default: assert(0 && "bit_depth should be AOM_BITS_8, AOM_BITS_10 or AOM_BITS_12"); return -1;
    }
}

int32_t svt_av1_compute_qdelta_fp(int32_t qstart_fp8, int32_t qtarget_fp8, AomBitDepth bit_depth) {
    int32_t start_index  = MAX_Q_INDEX;
    int32_t target_index = MAX_Q_INDEX;
    int32_t i;

    // Convert the average q value to an index.
    for (i = MIN_Q_INDEX; i < MAX_Q_INDEX; ++i) {
        start_index = i;
        if (svt_av1_convert_qindex_to_q_fp8(i, bit_depth) >= qstart_fp8)
            break;
    }

    // Convert the q target to an index
    for (i = MIN_Q_INDEX; i < MAX_Q_INDEX; ++i) {
        target_index = i;
        if (svt_av1_convert_qindex_to_q_fp8(i, bit_depth) >= qtarget_fp8)
            break;
    }

    return target_index - start_index;
}
double svt_av1_convert_qindex_to_q(int32_t qindex, AomBitDepth bit_depth) {
    // Convert the index to a real Q value (scaled down to match old Q values)
    switch (bit_depth) {
    case AOM_BITS_8: return svt_av1_ac_quant_q3(qindex, 0, bit_depth) / 4.0;
    case AOM_BITS_10: return svt_av1_ac_quant_q3(qindex, 0, bit_depth) / 16.0;
    case AOM_BITS_12: return svt_av1_ac_quant_q3(qindex, 0, bit_depth) / 64.0;
    default: assert(0 && "bit_depth should be AOM_BITS_8, AOM_BITS_10 or AOM_BITS_12"); return -1.0;
    }
}
int32_t svt_av1_compute_qdelta(double qstart, double qtarget, AomBitDepth bit_depth) {
    int32_t start_index  = MAX_Q_INDEX;
    int32_t target_index = MAX_Q_INDEX;
    int32_t i;

    // Convert the average q value to an index.
    for (i = MIN_Q_INDEX; i < MAX_Q_INDEX; ++i) {
        start_index = i;
        if (svt_av1_convert_qindex_to_q(i, bit_depth) >= qstart)
            break;
    }

    // Convert the q target to an index
    for (i = MIN_Q_INDEX; i < MAX_Q_INDEX; ++i) {
        target_index = i;
        if (svt_av1_convert_qindex_to_q(i, bit_depth) >= qtarget)
            break;
    }

    return target_index - start_index;
}
#define STATIC_MOTION_THRESH 95

// that are not marked as coded with 0,0 motion in the first pass.
#define FAST_MOVING_KF_GROUP_THRESH 5
#define MEDIUM_MOVING_KF_GROUP_THRESH 30
#define MAX_QPS_COMP_I 150
#define MAX_QPS_COMP_I_LR 42
#define MAX_QPS_COMP_NONI 300
#define HIGH_QPS_COMP_THRESHOLD 80
#define LOW_QPS_COMP_THRESHOLD 40
#define HIGH_FILTERED_THRESHOLD (4 << 8) // 8 bit precision
#define LOW_FILTERED_THRESHOLD (2 << 8) // 8 bit precision
#define MAX_REF_AREA_I 50 // Max ref area for I slice
#define MAX_REF_AREA_NONI 50 // Max ref area for Non I slice
#define MAX_REF_AREA_NONI_LOW_RES 40 // Max ref area for Non I slice in low resolution
#define REF_AREA_DIF_THRESHOLD 10 // Difference threshold for ref area between two frames
#define REF_AREA_LOW_THRESHOLD 8 // Low threshold for ref area
#define REF_AREA_MED_THRESHOLD 16 // Medium threshold for ref area
#define ME_SAD_LOW_THRESHOLD1 15 // Low sad_ threshold1 for me distortion (very low)
#define ME_SAD_LOW_THRESHOLD2 25 // Low sad_ threshold2 for me distortion (low)
#define ME_SAD_HIGH_THRESHOLD 80 // High sad_ threshold2 for me distortion (high)

#define SUPERRES_QADJ_PER_DENOM_KEYFRAME_SOLO 0
#define SUPERRES_QADJ_PER_DENOM_KEYFRAME 2
#define SUPERRES_QADJ_PER_DENOM_ARFFRAME 0

#define ASSIGN_MINQ_TABLE(bit_depth, name)                   \
    do {                                                     \
        name = NULL;                                         \
        switch (bit_depth) {                                 \
        case AOM_BITS_10: name = name##_10; break;           \
        case AOM_BITS_12: name = name##_12; break;           \
        case AOM_BITS_8: name = name##_8; break;             \
        }                                                    \
        assert(name &&                                       \
               "bit_depth should be AOM_BITS_8, AOM_BITS_10" \
               " or AOM_BITS_12");                           \
        if (!name)                                           \
            name = name##_8;                                 \
    } while (0)
static int kf_low_motion_minq_cqp_8[QINDEX_RANGE] = {
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   2,   2,   2,   2,   2,   2,   2,   3,
    3,   3,   3,   3,   3,   3,   4,   4,   4,   4,   4,   4,   4,   4,   5,   5,   5,   5,   5,
    5,   5,   6,   6,   6,   6,   6,   6,   6,   7,   7,   7,   7,   7,   7,   7,   7,   8,   8,
    8,   8,   8,   9,   9,   9,   9,   10,  10,  10,  10,  11,  11,  11,  11,  12,  12,  12,  12,
    13,  13,  13,  13,  14,  14,  14,  15,  15,  15,  16,  16,  16,  17,  17,  18,  18,  18,  19,
    19,  19,  20,  20,  20,  21,  21,  22,  22,  23,  23,  24,  24,  24,  25,  25,  26,  26,  27,
    27,  28,  28,  29,  30,  30,  31,  31,  32,  32,  33,  34,  34,  35,  36,  36,  37,  37,  38,
    39,  39,  40,  41,  42,  42,  43,  44,  45,  45,  46,  47,  48,  49,  50,  51,  51,  52,  53,
    54,  55,  56,  57,  58,  59,  60,  61,  62,  64,  65,  66,  67,  69,  70,  71,  72,  74,  75,
    77,  78,  80,  82,  83,  85,  87,  89,  91,  93,  95,  96,  97,  99,  100, 101, 103, 104, 105,
    107, 109, 110, 112, 114, 116, 118, 120, 122, 124, 125, 127, 129, 131, 134, 136, 138, 140, 142,
    144, 147, 149, 151, 154, 156, 158, 161, 163};
static int kf_low_motion_minq_cqp_10[QINDEX_RANGE] = {
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   10,  10,  11,
    11,  11,  11,  11,  11,  12,  12,  12,  12,  12,  13,  13,  13,  13,  13,  13,  13,  14,  14,
    14,  14,  14,  14,  14,  15,  15,  15,  15,  15,  16,  16,  16,  16,  16,  16,  16,  17,  17,
    17,  17,  17,  18,  18,  18,  18,  19,  19,  19,  19,  20,  20,  20,  21,  21,  21,  21,  22,
    22,  22,  22,  23,  23,  23,  23,  24,  24,  24,  25,  25,  25,  26,  26,  26,  27,  27,  27,
    28,  28,  28,  29,  29,  29,  30,  30,  31,  31,  32,  32,  32,  33,  33,  34,  34,  34,  35,
    35,  36,  36,  37,  37,  38,  38,  39,  39,  40,  40,  41,  41,  42,  42,  43,  44,  44,  45,
    46,  46,  47,  47,  48,  49,  49,  50,  51,  51,  52,  53,  54,  54,  55,  56,  57,  58,  58,
    59,  60,  61,  62,  63,  64,  65,  66,  67,  68,  69,  70,  71,  72,  73,  74,  76,  77,  78,
    80,  81,  83,  84,  86,  87,  89,  91,  93,  95,  96,  97,  98,  100, 101, 102, 103, 105, 106,
    108, 109, 111, 113, 115, 117, 119, 121, 122, 124, 126, 128, 130, 132, 134, 136, 138, 140, 142,
    144, 147, 149, 151, 154, 156, 159, 161, 163};
static int kf_low_motion_minq_cqp_12[QINDEX_RANGE] = {
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   13,
    13,  13,  13,  14,  14,  14,  14,  14,  14,  15,  15,  15,  15,  15,  16,  16,  16,  16,  16,
    16,  16,  17,  17,  17,  17,  17,  17,  18,  18,  18,  18,  18,  18,  18,  19,  19,  19,  19,
    19,  19,  20,  20,  20,  20,  21,  21,  21,  21,  22,  22,  22,  22,  23,  23,  23,  23,  24,
    24,  24,  24,  25,  25,  25,  25,  26,  26,  26,  27,  27,  27,  28,  28,  28,  29,  29,  29,
    30,  30,  30,  31,  31,  31,  32,  32,  33,  33,  33,  34,  34,  35,  35,  35,  36,  36,  37,
    37,  38,  38,  39,  39,  39,  40,  40,  41,  41,  42,  42,  43,  44,  44,  45,  45,  46,  46,
    47,  48,  48,  49,  49,  50,  51,  51,  52,  53,  53,  54,  55,  56,  56,  57,  58,  59,  59,
    60,  61,  62,  63,  64,  65,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,  76,  78,  79,
    80,  82,  83,  85,  86,  88,  90,  91,  93,  95,  96,  97,  99,  100, 101, 102, 104, 105, 106,
    108, 110, 111, 113, 115, 117, 119, 121, 122, 124, 126, 128, 130, 132, 134, 136, 138, 140, 142,
    144, 147, 149, 152, 154, 156, 159, 161, 163};
static int kf_high_motion_minq_8[QINDEX_RANGE] = {
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   2,   2,   3,   3,   4,   4,   5,
    5,   5,   6,   6,   7,   7,   8,   8,   8,   9,   9,   10,  10,  11,  11,  11,  12,  12,  13,
    13,  14,  14,  14,  15,  15,  16,  16,  16,  17,  17,  18,  18,  19,  19,  19,  20,  20,  21,
    21,  21,  22,  22,  23,  23,  24,  24,  24,  25,  25,  26,  26,  26,  27,  27,  28,  28,  28,
    29,  29,  30,  30,  30,  31,  31,  32,  32,  32,  33,  33,  34,  34,  34,  35,  35,  36,  36,
    36,  37,  38,  39,  39,  40,  41,  42,  42,  43,  44,  45,  46,  46,  47,  48,  49,  49,  50,
    51,  51,  52,  53,  54,  54,  55,  56,  57,  58,  59,  61,  62,  63,  64,  65,  66,  67,  68,
    69,  70,  71,  72,  73,  74,  76,  77,  78,  80,  81,  82,  84,  85,  86,  88,  89,  90,  92,
    93,  95,  96,  97,  97,  98,  99,  100, 100, 101, 102, 103, 104, 105, 106, 107, 107, 108, 109,
    110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 121, 122, 123, 124, 124, 125, 126,
    128, 128, 129, 130, 131, 131, 132, 134, 135, 136, 137, 138, 139, 140, 141, 143, 143, 144, 146,
    146, 147, 149, 150, 151, 152, 153, 155, 156, 158, 158, 160, 161, 163, 164, 166, 166, 168, 170,
    171, 173, 174, 176, 178, 179, 181, 183, 185, 187, 189, 191, 193, 195, 197, 200, 201, 204, 206,
    209, 212, 214, 216, 219, 222, 224, 227, 230};

static int arfgf_low_motion_minq_8[QINDEX_RANGE] = {
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   2,   2,   2,   3,   3,   3,   3,   4,   4,   4,   5,   5,   5,   5,   6,   6,   6,
    7,   7,   7,   7,   8,   8,   8,   9,   9,   9,   9,   10,  10,  10,  10,  11,  11,  11,  12,
    12,  12,  12,  13,  13,  13,  13,  14,  14,  14,  15,  15,  15,  15,  16,  16,  16,  16,  17,
    17,  17,  17,  18,  18,  18,  18,  19,  19,  19,  20,  20,  20,  20,  21,  21,  21,  21,  22,
    22,  22,  23,  23,  24,  24,  25,  25,  26,  26,  27,  27,  28,  28,  29,  29,  30,  30,  31,
    31,  32,  32,  33,  33,  34,  34,  35,  36,  36,  37,  38,  38,  39,  40,  41,  41,  42,  43,
    43,  44,  45,  45,  46,  47,  48,  49,  49,  50,  51,  52,  53,  54,  54,  55,  56,  57,  58,
    59,  60,  61,  62,  63,  64,  65,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,  77,  78,
    79,  80,  81,  83,  84,  85,  86,  87,  89,  90,  91,  92,  94,  95,  96,  97,  97,  98,  100,
    100, 101, 102, 102, 103, 105, 106, 106, 107, 109, 110, 110, 112, 113, 114, 116, 116, 118, 119,
    120, 122, 123, 125, 125, 127, 128, 130, 132, 133, 134, 135, 137, 139, 140, 141, 143, 145, 146,
    148, 150, 152, 154, 155, 158, 160, 162, 164, 166, 168, 171, 173, 176, 178, 181, 183, 186, 188,
    191, 194, 197, 200, 203, 206, 210, 213, 216};

static int arfgf_high_motion_minq_8[QINDEX_RANGE] = {
    0,   0,   0,   0,   0,   0,   0,   0,   2,   2,   3,   3,   4,   4,   5,   5,   6,   7,   7,
    8,   8,   9,   9,   10,  10,  11,  11,  12,  12,  13,  13,  14,  14,  15,  16,  16,  17,  17,
    18,  18,  19,  19,  20,  20,  21,  21,  22,  22,  23,  23,  24,  24,  25,  25,  26,  26,  27,
    27,  28,  28,  29,  29,  30,  31,  31,  32,  32,  33,  33,  34,  34,  35,  35,  36,  36,  37,
    37,  38,  38,  39,  39,  40,  40,  41,  41,  42,  42,  43,  43,  44,  44,  45,  45,  46,  46,
    46,  47,  48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,  64,
    65,  66,  67,  68,  68,  69,  70,  72,  73,  74,  76,  77,  79,  80,  81,  83,  84,  85,  87,
    88,  89,  91,  92,  93,  95,  96,  97,  98,  99,  100, 100, 101, 102, 103, 104, 105, 106, 107,
    108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 123, 124, 125,
    126, 127, 128, 129, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 139, 140, 141, 142,
    144, 144, 145, 146, 147, 148, 149, 151, 151, 152, 153, 155, 156, 156, 157, 159, 160, 161, 162,
    163, 164, 166, 167, 169, 169, 170, 172, 173, 175, 176, 178, 179, 180, 181, 183, 184, 186, 188,
    189, 191, 192, 194, 196, 197, 199, 201, 202, 204, 206, 209, 210, 212, 214, 217, 218, 220, 223,
    225, 228, 230, 232, 234, 237, 239, 242, 245};
static int inter_minq_8[QINDEX_RANGE] = {
    0,   0,   2,   2,   3,   4,   5,   6,   7,   8,   9,   10,  10,  11,  12,  13,  14,  15,  16,
    17,  18,  18,  19,  20,  21,  22,  23,  24,  25,  26,  26,  27,  28,  29,  30,  31,  32,  33,
    33,  34,  35,  36,  37,  38,  39,  40,  40,  41,  42,  43,  44,  45,  46,  47,  47,  48,  49,
    50,  51,  52,  53,  53,  54,  55,  56,  57,  58,  59,  59,  60,  61,  62,  63,  64,  65,  65,
    66,  67,  68,  69,  70,  71,  71,  72,  73,  74,  75,  76,  77,  77,  78,  79,  80,  81,  82,
    83,  84,  86,  88,  89,  91,  93,  94,  96,  97,  97,  98,  99,  100, 101, 102, 102, 103, 104,
    105, 106, 107, 107, 108, 109, 110, 111, 112, 114, 115, 116, 117, 119, 120, 121, 122, 122, 123,
    124, 125, 126, 127, 127, 128, 129, 131, 132, 133, 134, 135, 136, 137, 138, 139, 139, 140, 141,
    142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 157, 158, 159,
    161, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 176, 177,
    178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196,
    196, 197, 199, 199, 200, 201, 203, 203, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215,
    216, 217, 219, 220, 221, 222, 223, 225, 226, 227, 228, 230, 231, 232, 234, 235, 236, 238, 239,
    240, 242, 243, 245, 246, 248, 250, 251, 253};
static int rtc_minq_8[QINDEX_RANGE] = {
    0,   0,   0,   0,   0,   2,   3,   3,   4,   5,   5,   6,   7,   7,   8,   9,   9,   10,  11,
    12,  12,  13,  14,  14,  15,  16,  16,  17,  18,  18,  19,  20,  20,  21,  22,  22,  23,  24,
    24,  25,  26,  26,  27,  28,  28,  29,  30,  31,  31,  32,  33,  33,  34,  35,  35,  36,  37,
    37,  38,  39,  39,  40,  41,  41,  42,  42,  43,  44,  44,  45,  46,  46,  47,  48,  48,  49,
    50,  50,  51,  52,  52,  53,  54,  54,  55,  56,  56,  57,  58,  58,  59,  60,  60,  61,  61,
    62,  63,  65,  66,  67,  69,  70,  71,  72,  74,  75,  76,  78,  79,  80,  81,  83,  84,  85,
    86,  88,  89,  90,  91,  93,  94,  96,  97,  98,  98,  99,  100, 101, 102, 103, 104, 105, 106,
    107, 108, 109, 110, 110, 112, 113, 114, 115, 116, 118, 119, 120, 121, 122, 123, 123, 124, 125,
    126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 142, 143,
    144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162,
    162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180,
    181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 199, 200,
    201, 202, 203, 205, 206, 207, 208, 210, 211, 212, 214, 215, 216, 218, 219, 221, 222, 224, 225,
    227, 229, 230, 232, 234, 235, 237, 239, 241};

static int kf_high_motion_minq_10[QINDEX_RANGE] = {
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   11,  11,  11,  12,  13,  13,  14,  14,  15,  15,  16,  16,  17,  17,  18,  18,  19,
    19,  20,  20,  21,  21,  22,  22,  22,  23,  23,  24,  24,  25,  25,  26,  26,  27,  27,  27,
    28,  28,  29,  29,  29,  30,  30,  31,  31,  32,  32,  32,  33,  33,  33,  34,  34,  35,  35,
    35,  36,  36,  37,  37,  37,  38,  38,  39,  39,  39,  40,  40,  41,  41,  41,  42,  42,  42,
    43,  43,  44,  45,  45,  46,  47,  48,  48,  49,  50,  50,  51,  52,  52,  53,  54,  54,  55,
    56,  56,  57,  58,  58,  59,  60,  61,  62,  63,  64,  64,  66,  67,  67,  69,  69,  70,  71,
    72,  73,  74,  75,  76,  77,  79,  80,  81,  82,  84,  85,  86,  87,  88,  90,  91,  92,  94,
    95,  96,  97,  97,  98,  99,  100, 101, 101, 102, 103, 104, 105, 105, 106, 107, 108, 109, 110,
    111, 112, 113, 114, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 123, 124, 125, 125, 126,
    128, 129, 129, 130, 131, 132, 133, 134, 135, 136, 137, 139, 139, 140, 141, 143, 143, 144, 146,
    147, 147, 149, 150, 151, 152, 153, 155, 156, 158, 159, 160, 161, 163, 164, 166, 166, 168, 170,
    171, 173, 174, 176, 178, 179, 181, 184, 185, 187, 189, 191, 193, 195, 197, 200, 201, 204, 206,
    209, 212, 214, 216, 219, 222, 224, 227, 230};

static int arfgf_low_motion_minq_10[QINDEX_RANGE] = {
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   10,  11,  11,  11,  12,  12,  12,  13,  13,
    13,  14,  14,  14,  15,  15,  16,  16,  16,  17,  17,  17,  17,  18,  18,  18,  19,  19,  19,
    20,  20,  20,  21,  21,  21,  21,  22,  22,  22,  23,  23,  23,  24,  24,  24,  24,  25,  25,
    25,  25,  26,  26,  26,  26,  27,  27,  27,  28,  28,  28,  28,  28,  29,  29,  29,  30,  30,
    30,  30,  31,  31,  32,  32,  33,  33,  34,  34,  35,  35,  36,  36,  37,  37,  37,  38,  38,
    39,  39,  40,  40,  41,  41,  41,  42,  42,  43,  44,  44,  45,  46,  46,  47,  48,  48,  49,
    49,  50,  50,  51,  52,  52,  53,  54,  55,  56,  56,  57,  58,  59,  59,  60,  61,  62,  62,
    63,  64,  65,  66,  67,  68,  69,  69,  70,  72,  72,  73,  74,  75,  77,  77,  78,  79,  80,
    82,  83,  84,  85,  86,  87,  88,  90,  91,  92,  93,  94,  95,  96,  97,  98,  98,  99,  101,
    101, 102, 103, 103, 104, 106, 106, 107, 108, 110, 110, 111, 113, 114, 114, 116, 117, 119, 120,
    121, 122, 123, 125, 126, 128, 129, 131, 132, 133, 135, 136, 137, 139, 140, 142, 144, 145, 146,
    148, 150, 152, 154, 156, 158, 160, 162, 164, 166, 169, 171, 173, 176, 178, 181, 184, 186, 189,
    191, 194, 197, 200, 203, 206, 210, 213, 216};

static int arfgf_high_motion_minq_10[QINDEX_RANGE] = {
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   10,  11,
    11,  12,  13,  13,  14,  14,  15,  16,  16,  17,  18,  18,  19,  19,  20,  20,  21,  22,  22,
    23,  23,  24,  24,  25,  26,  26,  27,  27,  28,  28,  29,  30,  30,  30,  31,  32,  32,  33,
    33,  34,  34,  35,  35,  36,  36,  37,  37,  38,  38,  39,  39,  40,  40,  41,  41,  42,  42,
    42,  43,  44,  44,  44,  45,  45,  46,  46,  47,  47,  48,  48,  49,  49,  50,  50,  51,  51,
    52,  52,  53,  54,  55,  56,  57,  58,  59,  60,  60,  61,  62,  63,  64,  65,  66,  67,  67,
    68,  69,  70,  71,  72,  72,  73,  75,  76,  77,  78,  80,  81,  82,  84,  85,  86,  87,  89,
    90,  91,  92,  94,  95,  96,  97,  98,  99,  99,  100, 101, 102, 103, 104, 105, 105, 106, 107,
    108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 120, 121, 121, 122, 123, 124, 125, 125,
    126, 127, 128, 129, 130, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 140, 141, 142,
    144, 145, 145, 146, 147, 148, 149, 151, 152, 152, 153, 155, 156, 156, 157, 159, 160, 161, 163,
    163, 164, 166, 167, 169, 170, 170, 172, 173, 175, 176, 178, 179, 181, 181, 183, 184, 186, 188,
    189, 191, 192, 194, 196, 197, 199, 201, 202, 204, 206, 209, 210, 212, 214, 217, 218, 220, 223,
    225, 228, 230, 232, 234, 237, 240, 242, 245};

static int inter_minq_10[QINDEX_RANGE] = {
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   11,  11,  12,  13,  14,  15,  16,  17,
    18,  19,  20,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  29,  30,  31,  32,  33,  34,
    35,  36,  37,  37,  39,  39,  40,  41,  42,  43,  44,  44,  45,  46,  47,  48,  49,  50,  51,
    51,  52,  53,  54,  55,  56,  57,  58,  58,  59,  60,  61,  62,  62,  63,  64,  65,  66,  67,
    68,  69,  69,  70,  71,  72,  73,  73,  74,  75,  76,  77,  78,  79,  79,  80,  81,  82,  83,
    84,  85,  87,  88,  90,  92,  93,  95,  96,  97,  98,  99,  99,  100, 101, 102, 103, 104, 104,
    105, 106, 107, 108, 109, 109, 110, 111, 113, 114, 115, 116, 118, 119, 120, 121, 122, 123, 123,
    124, 125, 126, 127, 127, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 140, 141,
    142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 158, 160,
    161, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 177,
    178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196,
    196, 197, 199, 199, 200, 201, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215,
    216, 218, 219, 220, 221, 222, 223, 225, 226, 227, 228, 230, 231, 232, 234, 235, 236, 238, 239,
    240, 242, 243, 245, 246, 248, 250, 251, 253};
static int rtc_minq_10[QINDEX_RANGE] = {
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   11,  11,  12,  13,  13,
    14,  15,  16,  16,  17,  18,  19,  19,  20,  21,  22,  22,  23,  24,  24,  25,  26,  27,  28,
    28,  29,  29,  30,  31,  32,  32,  33,  34,  34,  35,  36,  36,  37,  38,  38,  39,  40,  41,
    41,  42,  42,  43,  44,  44,  45,  46,  46,  47,  48,  48,  49,  50,  50,  51,  51,  52,  53,
    53,  54,  55,  55,  56,  56,  57,  58,  58,  59,  60,  60,  61,  62,  62,  63,  63,  64,  64,
    65,  67,  68,  69,  70,  71,  72,  74,  75,  76,  77,  78,  80,  81,  82,  83,  84,  86,  87,
    88,  89,  90,  91,  93,  94,  95,  96,  97,  98,  99,  100, 101, 102, 103, 104, 105, 105, 106,
    107, 108, 109, 110, 111, 112, 113, 114, 116, 117, 118, 119, 120, 121, 122, 123, 124, 124, 125,
    126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144,
    144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162,
    163, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180,
    181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 198, 199, 200,
    201, 202, 203, 205, 206, 207, 208, 210, 211, 212, 214, 215, 216, 218, 219, 221, 222, 224, 225,
    227, 229, 230, 232, 234, 235, 237, 239, 241};

static int kf_high_motion_minq_12[QINDEX_RANGE] = {
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   13,  14,  14,  15,  15,  16,  16,  17,  17,  18,  18,  19,  19,  20,  20,
    21,  21,  22,  22,  23,  23,  23,  24,  24,  25,  25,  26,  26,  27,  27,  28,  28,  28,  29,
    29,  30,  30,  31,  31,  31,  32,  32,  33,  33,  33,  34,  34,  35,  35,  35,  36,  36,  37,
    37,  37,  38,  38,  39,  39,  39,  40,  40,  40,  41,  41,  41,  42,  42,  43,  43,  43,  44,
    44,  45,  45,  46,  47,  47,  48,  49,  49,  50,  51,  51,  52,  53,  53,  54,  55,  55,  56,
    57,  57,  58,  59,  59,  60,  61,  62,  63,  64,  64,  65,  66,  67,  68,  69,  70,  71,  72,
    73,  74,  75,  76,  77,  78,  79,  80,  82,  83,  84,  85,  86,  88,  89,  90,  91,  92,  94,
    95,  96,  97,  98,  98,  99,  100, 101, 101, 102, 103, 104, 105, 106, 107, 107, 108, 109, 110,
    111, 112, 113, 114, 115, 115, 116, 117, 118, 119, 120, 121, 122, 123, 123, 124, 125, 125, 126,
    127, 128, 128, 129, 130, 131, 132, 132, 133, 134, 135, 136, 137, 137, 138, 139, 139, 140, 141,
    142, 142, 143, 144, 145, 145, 146, 147, 148, 149, 150, 151, 151, 152, 153, 154, 155, 155, 156,
    157, 158, 159, 160, 161, 162, 163, 165, 166, 167, 168, 170, 171, 172, 173, 175, 176, 178, 179,
    181, 183, 184, 186, 188, 190, 191, 193, 195};

static int arfgf_low_motion_minq_12[QINDEX_RANGE] = {
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   13,  13,  14,  14,  14,  15,  15,
    15,  16,  16,  16,  17,  17,  17,  18,  18,  18,  19,  19,  19,  20,  20,  20,  21,  21,  21,
    22,  22,  22,  22,  23,  23,  23,  24,  24,  24,  25,  25,  25,  25,  26,  26,  26,  26,  27,
    27,  27,  28,  28,  28,  28,  29,  29,  29,  29,  30,  30,  30,  30,  31,  31,  31,  31,  32,
    32,  32,  33,  33,  34,  34,  35,  35,  35,  36,  36,  37,  37,  38,  38,  39,  39,  39,  40,
    40,  41,  41,  42,  42,  42,  43,  43,  44,  45,  45,  46,  46,  47,  48,  48,  49,  49,  50,
    51,  51,  52,  52,  53,  54,  54,  55,  56,  57,  57,  58,  59,  60,  60,  61,  62,  63,  63,
    64,  65,  66,  67,  68,  69,  70,  70,  71,  72,  73,  74,  75,  76,  77,  78,  79,  80,  81,
    82,  83,  84,  86,  87,  88,  89,  90,  91,  92,  94,  95,  96,  96,  97,  98,  98,  99,  100,
    100, 101, 102, 102, 103, 104, 105, 105, 106, 107, 108, 108, 109, 110, 111, 111, 112, 113, 114,
    115, 115, 116, 117, 118, 119, 120, 121, 122, 122, 123, 124, 124, 125, 126, 127, 128, 129, 129,
    130, 131, 132, 134, 135, 136, 137, 138, 139, 141, 142, 143, 144, 146, 147, 149, 151, 152, 154,
    155, 157, 159, 161, 163, 165, 167, 169, 171};

static int arfgf_high_motion_minq_12[QINDEX_RANGE] = {
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   13,  14,  14,  15,  16,  16,  17,  17,  18,  19,  19,  20,  20,  21,  22,  22,  23,  23,
    24,  25,  25,  26,  26,  27,  27,  28,  28,  29,  30,  30,  31,  31,  32,  32,  33,  33,  34,
    34,  35,  35,  36,  36,  37,  37,  38,  38,  39,  39,  40,  40,  41,  41,  42,  42,  43,  43,
    44,  44,  45,  45,  46,  46,  47,  47,  48,  48,  49,  49,  49,  50,  50,  51,  51,  52,  52,
    53,  53,  54,  55,  56,  57,  58,  59,  59,  60,  61,  62,  63,  64,  65,  65,  66,  67,  68,
    69,  70,  71,  71,  72,  73,  74,  75,  77,  78,  79,  80,  82,  83,  84,  85,  87,  88,  89,
    90,  92,  93,  94,  95,  96,  97,  98,  99,  100, 101, 101, 102, 103, 104, 105, 106, 106, 107,
    108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 119, 120, 121, 122, 122, 123, 124, 125, 125,
    126, 127, 128, 129, 130, 131, 132, 132, 133, 134, 135, 136, 137, 138, 139, 140, 140, 141, 142,
    143, 144, 144, 145, 146, 147, 148, 149, 150, 150, 151, 152, 153, 154, 154, 155, 156, 157, 158,
    158, 159, 160, 161, 162, 163, 163, 164, 165, 166, 167, 168, 169, 170, 170, 171, 172, 173, 174,
    175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 187, 188, 189, 190, 192, 193, 194, 196,
    197, 199, 200, 202, 203, 205, 207, 208, 210};

static int inter_minq_12[QINDEX_RANGE] = {
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   13,  14,  15,  16,  17,
    18,  19,  20,  21,  22,  23,  23,  24,  25,  26,  27,  28,  29,  30,  31,  32,  32,  33,  34,
    35,  36,  37,  38,  39,  40,  40,  41,  42,  43,  44,  45,  46,  47,  47,  48,  49,  50,  51,
    52,  53,  53,  54,  55,  56,  57,  58,  59,  59,  60,  61,  62,  63,  64,  65,  65,  66,  67,
    68,  69,  70,  70,  71,  72,  73,  74,  75,  76,  76,  77,  78,  79,  80,  80,  81,  82,  83,
    84,  85,  87,  89,  90,  92,  93,  95,  96,  97,  98,  99,  99,  100, 101, 102, 103, 104, 104,
    105, 106, 107, 108, 109, 109, 110, 111, 113, 114, 115, 116, 118, 119, 120, 121, 122, 123, 123,
    124, 125, 126, 127, 127, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 140, 141,
    142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 158, 160,
    161, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 177,
    178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196,
    196, 197, 199, 199, 200, 201, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215,
    216, 217, 219, 220, 221, 222, 223, 225, 226, 227, 228, 230, 231, 232, 234, 235, 236, 238, 239,
    240, 242, 243, 245, 246, 248, 250, 251, 253};
static int rtc_minq_12[QINDEX_RANGE] = {
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   13,  14,
    15,  16,  16,  17,  18,  19,  19,  20,  21,  22,  22,  23,  24,  25,  25,  26,  27,  28,  28,
    29,  30,  30,  31,  32,  32,  33,  34,  34,  35,  36,  37,  37,  38,  39,  39,  40,  41,  41,
    42,  43,  43,  44,  45,  45,  46,  46,  47,  48,  48,  49,  50,  50,  51,  52,  52,  53,  54,
    54,  55,  55,  56,  57,  57,  58,  58,  59,  60,  60,  61,  62,  62,  63,  63,  64,  65,  65,
    66,  67,  68,  69,  71,  72,  73,  74,  75,  76,  78,  79,  80,  81,  82,  84,  85,  86,  87,
    88,  90,  91,  92,  93,  94,  95,  96,  97,  98,  99,  100, 101, 102, 103, 104, 105, 106, 107,
    107, 108, 109, 110, 111, 112, 113, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 124, 125,
    126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144,
    145, 146, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162,
    163, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180,
    181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 199, 200,
    201, 202, 203, 205, 206, 207, 208, 210, 211, 212, 214, 215, 216, 218, 219, 221, 222, 224, 225,
    227, 229, 230, 232, 234, 235, 237, 239, 241};
static int gf_high_tpl_la = 2400;
static int gf_low_tpl_la  = 300;
static int kf_high        = 5000;
static int kf_low         = 400;
static int get_active_quality(int q, int gfu_boost, int low, int high, int *low_motion_minq,
                              int *high_motion_minq) {
    if (gfu_boost > high)
        return low_motion_minq[q];
    else if (gfu_boost < low)
        return high_motion_minq[q];
    else {
        const int gap        = high - low;
        const int offset     = high - gfu_boost;
        const int qdiff      = high_motion_minq[q] - low_motion_minq[q];
        const int adjustment = ((offset * qdiff) + (gap >> 1)) / gap;
        return low_motion_minq[q] + adjustment;
    }
}
static int get_kf_active_quality_tpl(const RATE_CONTROL *const rc, int q, AomBitDepth bit_depth) {
    int *kf_low_motion_minq_cqp;
    int *kf_high_motion_minq;
    ASSIGN_MINQ_TABLE(bit_depth, kf_low_motion_minq_cqp);
    ASSIGN_MINQ_TABLE(bit_depth, kf_high_motion_minq);
    return get_active_quality(
        q, rc->kf_boost, kf_low, kf_high, kf_low_motion_minq_cqp, kf_high_motion_minq);
}
static int get_gf_active_quality_tpl_la(const RATE_CONTROL *const rc, int q,
                                        AomBitDepth bit_depth) {
    int *arfgf_low_motion_minq;
    int *arfgf_high_motion_minq;
    ASSIGN_MINQ_TABLE(bit_depth, arfgf_low_motion_minq);
    ASSIGN_MINQ_TABLE(bit_depth, arfgf_high_motion_minq);
    return get_active_quality(q,
                              rc->gfu_boost,
                              gf_low_tpl_la,
                              gf_high_tpl_la,
                              arfgf_low_motion_minq,
                              arfgf_high_motion_minq);
}
static int get_gf_high_motion_quality(int q, AomBitDepth bit_depth) {
    int *arfgf_high_motion_minq;
    ASSIGN_MINQ_TABLE(bit_depth, arfgf_high_motion_minq);
    return arfgf_high_motion_minq[q];
}
int16_t svt_av1_dc_quant_qtx(int32_t qindex, int32_t delta, AomBitDepth bit_depth);

static int get_cqp_kf_boost_from_r0(double r0, int frames_to_key,
                                    EbInputResolution input_resolution) {
    double factor;
    // when frames_to_key not available, it is set to -1. In this case the factor is set to average of min and max
    if (frames_to_key == -1)
        factor = (10.0 + 4.0) / 2;
    else {
        factor = sqrt((double)frames_to_key);
        factor = AOMMIN(factor, 10.0);
        factor = AOMMAX(factor, 4.0);
    }
    const int is_720p_or_smaller = input_resolution <= INPUT_SIZE_720p_RANGE;
    const int boost = is_720p_or_smaller ? (int)rint(3 * (75.0 + 17.0 * factor) / 2 / r0)
                                         : (int)rint(2 * (75.0 + 17.0 * factor) / r0);
    return boost;
}

double svt_av1_get_gfu_boost_projection_factor(double min_factor, double max_factor,
                                               int frame_count) {
    double factor = sqrt((double)frame_count);
    factor        = AOMMIN(factor, max_factor);
    factor        = AOMMAX(factor, min_factor);
    factor        = (200.0 + 10.0 * factor);
    return factor;
}

#define MAX_GFUBOOST_FACTOR 10.0
//#define MIN_GFUBOOST_FACTOR 4.0
static int get_gfu_boost_from_r0_lap(double min_factor, double max_factor, double r0,
                                     int frames_to_key) {
    double factor = svt_av1_get_gfu_boost_projection_factor(min_factor, max_factor, frames_to_key);
    const int boost = (int)rint(factor / r0);
    return boost;
}

int svt_av1_get_deltaq_offset(AomBitDepth bit_depth, int qindex, double beta, EB_SLICE slice_type) {
    assert(beta > 0.0);
    int q = svt_av1_dc_quant_qtx(qindex, 0, bit_depth);
    int newq;
    // use a less aggressive action when lowering the q for non I_slice
    if (slice_type != I_SLICE && beta > 1)
        newq = (int)rint(q / sqrt(sqrt(beta)));
    else
        newq = (int)rint(q / sqrt(beta));
    int orig_qindex = qindex;
    if (newq < q) {
        do {
            qindex--;
            q = svt_av1_dc_quant_qtx(qindex, 0, bit_depth);
        } while (newq < q && qindex > 0);
    } else {
        do {
            qindex++;
            q = svt_av1_dc_quant_qtx(qindex, 0, bit_depth);
        } while (newq > q && qindex < MAXQ);
    }
    return qindex - orig_qindex;
}

#define MIN_BPB_FACTOR 0.005
#define MAX_BPB_FACTOR 50
int svt_av1_rc_bits_per_mb(FrameType frame_type, int qindex, double correction_factor,
                           const int bit_depth, const int is_screen_content_type,
                           int onepass_cbr_mode) {
    const double q          = svt_av1_convert_qindex_to_q(qindex, bit_depth);
    int          enumerator = frame_type == KEY_FRAME ? 1400000 : 1000000;
    if (onepass_cbr_mode == 2) {
        enumerator = frame_type == KEY_FRAME ? 7000000 : 5250000;
    } else if (onepass_cbr_mode == 1) {
        enumerator = frame_type == KEY_FRAME ? 2000000 : 1500000;
    }
    if (is_screen_content_type) {
        enumerator = frame_type == KEY_FRAME ? 1000000 : 750000;
    }
    assert(correction_factor <= MAX_BPB_FACTOR && correction_factor >= MIN_BPB_FACTOR);

    // q based adjustment to baseline enumerator
    return (int)(enumerator * correction_factor / q);
}

static int find_qindex_by_rate(int desired_bits_per_mb, const int bit_depth, FrameType frame_type,
                               const int is_screen_content_type, int onepass_cbr_mode,
                               int best_qindex, int worst_qindex) {
    assert(best_qindex <= worst_qindex);
    int low  = best_qindex;
    int high = worst_qindex;
    while (low < high) {
        const int mid             = (low + high) >> 1;
        const int mid_bits_per_mb = svt_av1_rc_bits_per_mb(
            frame_type, mid, 1.0, bit_depth, is_screen_content_type, onepass_cbr_mode);
        if (mid_bits_per_mb > desired_bits_per_mb) {
            low = mid + 1;
        } else {
            high = mid;
        }
    }
    assert(low == high);
    assert(svt_av1_rc_bits_per_mb(
               frame_type, low, 1.0, bit_depth, is_screen_content_type, onepass_cbr_mode) <=
               desired_bits_per_mb ||
           low == worst_qindex);
    return low;
}

int svt_av1_compute_qdelta_by_rate(const RATE_CONTROL *rc, FrameType frame_type, int qindex,
                                   double rate_target_ratio, const int bit_depth,
                                   const int is_screen_content_type) {
    // Look up the current projected bits per block for the base index
    const int base_bits_per_mb = svt_av1_rc_bits_per_mb(
        frame_type, qindex, 1.0, bit_depth, is_screen_content_type, rc->onepass_cbr_mode);

    // Find the target bits per mb based on the base value and given ratio.
    const int target_bits_per_mb = (int)(rate_target_ratio * base_bits_per_mb);

    const int target_index = find_qindex_by_rate(target_bits_per_mb,
                                                 bit_depth,
                                                 frame_type,
                                                 is_screen_content_type,
                                                 rc->onepass_cbr_mode,
                                                 rc->best_quality,
                                                 rc->worst_quality);
    return target_index - qindex;
}

static const double rate_factor_deltas[RATE_FACTOR_LEVELS] = {
    1.00, // INTER_NORMAL
    1.00, // INTER_LOW
    1.00, // INTER_HIGH
    1.50, // GF_ARF_LOW
    2.00, // GF_ARF_STD
    2.00, // KF_STD
};

int svt_av1_frame_type_qdelta(RATE_CONTROL *rc, int rf_level, int q, const int bit_depth,
                              const int sc_content_detected) {
    const int /*rate_factor_level*/ rf_lvl     = rf_level; //get_rate_factor_level(&cpi->gf_group);
    const FrameType                 frame_type = (rf_lvl == KF_STD) ? KEY_FRAME : INTER_FRAME;
    double                          rate_factor;

    rate_factor = rate_factor_deltas[rf_lvl];
    if (rf_lvl == GF_ARF_LOW) {
        rate_factor -= (0 /*cpi->gf_group.layer_depth[cpi->gf_group.index]*/ - 2) * 0.1;
        rate_factor = AOMMAX(rate_factor, 1.0);
    }
    return svt_av1_compute_qdelta_by_rate(
        rc, frame_type, q, rate_factor, bit_depth, sc_content_detected);
}

static const rate_factor_level rate_factor_levels[FRAME_UPDATE_TYPES] = {
    KF_STD, // KF_UPDATE
    INTER_NORMAL, // LF_UPDATE
    GF_ARF_STD, // GF_UPDATE
    GF_ARF_STD, // ARF_UPDATE
    INTER_NORMAL, // OVERLAY_UPDATE
    INTER_NORMAL, // INTNL_OVERLAY_UPDATE
    GF_ARF_LOW, // INTNL_ARF_UPDATE
};

static rate_factor_level get_rate_factor_level(const GF_GROUP *const gf_group,
                                               unsigned char         gf_group_index) {
    const FRAME_UPDATE_TYPE update_type = gf_group->update_type[gf_group_index];
    assert(update_type < FRAME_UPDATE_TYPES);
    return rate_factor_levels[update_type];
}

int av1_frame_type_qdelta_org(RATE_CONTROL *rc, GF_GROUP *gf_group, unsigned char gf_group_index,
                              int q, const int bit_depth, uint8_t sc_content_detected) {
    const rate_factor_level rf_lvl     = get_rate_factor_level(gf_group, gf_group_index);
    const FrameType         frame_type = (rf_lvl == KF_STD) ? KEY_FRAME : INTER_FRAME;
    double                  rate_factor;

    rate_factor = rate_factor_deltas[rf_lvl];
    //anaghdin: to remove?
    if (rf_lvl == GF_ARF_LOW) {
        rate_factor -= (gf_group->layer_depth[gf_group_index] - 2) * 0.1;
        rate_factor = AOMMAX(rate_factor, 1.0);
    }
    return svt_av1_compute_qdelta_by_rate(
        rc, frame_type, q, rate_factor, bit_depth, sc_content_detected);
}

static void adjust_active_best_and_worst_quality_org(PictureControlSet *pcs_ptr, RATE_CONTROL *rc,
                                                     int *active_worst, int *active_best) {
    int                 active_best_quality  = *active_best;
    int                 active_worst_quality = *active_worst;
    SequenceControlSet *scs_ptr              = pcs_ptr->parent_pcs_ptr->scs_ptr;
    const int           bit_depth            = scs_ptr->static_config.encoder_bit_depth;

    EncodeContext         *encode_context_ptr = scs_ptr->encode_context_ptr;
    TWO_PASS *const        twopass            = &scs_ptr->twopass;
    const enum aom_rc_mode rc_mode            = encode_context_ptr->rc_cfg.mode;
    GF_GROUP              *gf_group           = &encode_context_ptr->gf_group;
    // Extension to max or min Q if undershoot or overshoot is outside
    // the permitted range.
    if (rc_mode != AOM_Q) {
        if (frame_is_intra_only(pcs_ptr->parent_pcs_ptr) ||
            (pcs_ptr->parent_pcs_ptr->temporal_layer_index < 2 && scs_ptr->is_short_clip) ||
            (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag && !scs_ptr->is_short_clip)) {
            active_best_quality -= (twopass->extend_minq + twopass->extend_minq_fast);
            active_worst_quality += (twopass->extend_maxq / 2);
        } else {
            active_best_quality -= (twopass->extend_minq + twopass->extend_minq_fast) / 2;
            active_worst_quality += twopass->extend_maxq;
        }
    }

    // Static forced key frames Q restrictions dealt with elsewhere.
    const int qdelta     = av1_frame_type_qdelta_org(rc,
                                                 gf_group,
                                                 pcs_ptr->parent_pcs_ptr->gf_group_index,
                                                 active_worst_quality,
                                                 bit_depth,
                                                 pcs_ptr->parent_pcs_ptr->sc_class1);
    active_worst_quality = AOMMAX(active_worst_quality + qdelta, active_best_quality);
    active_best_quality  = clamp(active_best_quality, rc->best_quality, rc->worst_quality);
    active_worst_quality = clamp(active_worst_quality, active_best_quality, rc->worst_quality);

    *active_best  = active_best_quality;
    *active_worst = active_worst_quality;
}

static void adjust_active_best_and_worst_quality(PictureControlSet *pcs_ptr, RATE_CONTROL *rc,
                                                 int rf_level, int *active_worst,
                                                 int *active_best) {
    int active_best_quality  = *active_best;
    int active_worst_quality = *active_worst;
    ;
    SequenceControlSet *scs_ptr   = pcs_ptr->parent_pcs_ptr->scs_ptr;
    const int           bit_depth = scs_ptr->static_config.encoder_bit_depth;

    // Static forced key frames Q restrictions dealt with elsewhere.
    if (!frame_is_intra_only(pcs_ptr->parent_pcs_ptr)
        /*|| (cpi->twopass.last_kfgroup_zeromotion_pct < STATIC_MOTION_THRESH)*/) {
        const int qdelta = svt_av1_frame_type_qdelta(
            rc, rf_level, active_worst_quality, bit_depth, pcs_ptr->parent_pcs_ptr->sc_class1);
        active_worst_quality = AOMMAX(active_worst_quality + qdelta, active_best_quality);
    }

    active_best_quality  = clamp(active_best_quality, rc->best_quality, rc->worst_quality);
    active_worst_quality = clamp(active_worst_quality, active_best_quality, rc->worst_quality);

    *active_best  = active_best_quality;
    *active_worst = active_worst_quality;
}

/******************************************************
 * cqp_qindex_calc_tpl_la
 * Assign the q_index per frame.
 * Used in the one pass encoding with tpl stats
 ******************************************************/
static int cqp_qindex_calc_tpl_la(PictureControlSet *pcs_ptr, RATE_CONTROL *rc, int qindex) {
    SequenceControlSet *scs_ptr              = pcs_ptr->parent_pcs_ptr->scs_ptr;
    const int           cq_level             = qindex;
    int                 active_best_quality  = 0;
    int                 active_worst_quality = qindex;
    rc->arf_q                                = 0;
    int q;
    int refresh_golden_frame, refresh_alt_ref_frame, is_intrl_arf_boost, rf_level;
    refresh_golden_frame  = frame_is_intra_only(pcs_ptr->parent_pcs_ptr) ? 1 : 0;
    refresh_alt_ref_frame = (pcs_ptr->parent_pcs_ptr->temporal_layer_index == 0) ? 1 : 0;
    is_intrl_arf_boost    = (pcs_ptr->parent_pcs_ptr->temporal_layer_index > 0 &&
                          pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag)
           ? 1
           : 0;
    rf_level              = (frame_is_intra_only(pcs_ptr->parent_pcs_ptr))  ? KF_STD
                     : (pcs_ptr->parent_pcs_ptr->temporal_layer_index == 0) ? GF_ARF_STD
                     : pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag   ? GF_ARF_LOW
                                                                            : INTER_NORMAL;

    const int bit_depth = scs_ptr->static_config.encoder_bit_depth;
    // Since many frames can be processed at the same time, storing/using arf_q in rc param is not sufficient and will create a run to run.
    // So, for each frame, arf_q is updated based on the qp of its references.
    rc->arf_q = MAX(rc->arf_q, ((pcs_ptr->ref_pic_qp_array[0][0] << 2) + 2));
    if (pcs_ptr->slice_type == B_SLICE)
        rc->arf_q = MAX(rc->arf_q, ((pcs_ptr->ref_pic_qp_array[1][0] << 2) + 2));

    if (frame_is_intra_only(pcs_ptr->parent_pcs_ptr)) {
        // Not forced keyframe.
        double q_adj_factor = 1.0;
        double q_val;
        rc->worst_quality = MAXQ;
        rc->best_quality  = MINQ;
        // The new tpl only looks at pictures in tpl group, which is fewer than before,
        // As a results, we defined a factor to adjust r0
        if (pcs_ptr->parent_pcs_ptr->frm_hdr.frame_type != KEY_FRAME) {
            double factor;
            factor                      = 1;
            pcs_ptr->parent_pcs_ptr->r0 = pcs_ptr->parent_pcs_ptr->r0 / factor;
        }
        pcs_ptr->parent_pcs_ptr->r0 = pcs_ptr->parent_pcs_ptr->r0 /
            tpl_hl_islice_div_factor[scs_ptr->max_heirachical_level];
        if (pcs_ptr->parent_pcs_ptr->frm_hdr.frame_type == KEY_FRAME) {
            {
                double factor = 1.0;
                if (pcs_ptr->parent_pcs_ptr->r0 < 0.2) {
                    double mult = 1.0;
                    factor      = (double)(mult * 255.0) / qindex;
                }
                pcs_ptr->parent_pcs_ptr->r0 = pcs_ptr->parent_pcs_ptr->r0 / factor;
            }
        }
        if (!scs_ptr->lad_mg) {
            double div_factor           = pcs_ptr->parent_pcs_ptr->tpl_ctrls.r0_adjust_factor
                          ? pcs_ptr->parent_pcs_ptr->used_tpl_frame_num * 0.2
                          : 1.0;
            pcs_ptr->parent_pcs_ptr->r0 = pcs_ptr->parent_pcs_ptr->r0 / div_factor;
        }
        // when frames_to_key not available, i.e. in 1 pass encoding
        rc->kf_boost = get_cqp_kf_boost_from_r0(
            pcs_ptr->parent_pcs_ptr->r0, -1, scs_ptr->input_resolution);
        int max_boost = pcs_ptr->parent_pcs_ptr->used_tpl_frame_num * KB;
        rc->kf_boost  = AOMMIN(rc->kf_boost, max_boost);
        // Baseline value derived from cpi->active_worst_quality and kf boost.
        active_best_quality = get_kf_active_quality_tpl(rc, active_worst_quality, bit_depth);
        // Allow somewhat lower kf minq with small image formats.
        if (pcs_ptr->parent_pcs_ptr->input_resolution == INPUT_SIZE_240p_RANGE)
            q_adj_factor -= (pcs_ptr->parent_pcs_ptr->tune_tpl_for_chroma) ? 0.2 : 0.15;
        // Make a further adjustment based on the kf zero motion measure.

        // Convert the adjustment factor to a qindex delta
        // on active_best_quality.
        q_val = svt_av1_convert_qindex_to_q(active_best_quality, bit_depth);
        active_best_quality += svt_av1_compute_qdelta(q_val, q_val * q_adj_factor, bit_depth);
    } else if (refresh_golden_frame || is_intrl_arf_boost || refresh_alt_ref_frame) {
        //double min_boost_factor = sqrt(1 << pcs_ptr->parent_pcs_ptr->hierarchical_levels);
        double min_boost_factor = (int32_t)1 << (pcs_ptr->parent_pcs_ptr->hierarchical_levels >> 1);
        if (pcs_ptr->parent_pcs_ptr->hierarchical_levels & 1) {
            min_boost_factor *= CONST_SQRT2;
        }
        // The new tpl only looks at pictures in tpl group, which is fewer than before,
        // As a results, we defined a factor to adjust r0
        if (pcs_ptr->parent_pcs_ptr->temporal_layer_index == 0) {
            double div_factor = 1;
            if (pcs_ptr->parent_pcs_ptr->tpl_ctrls.r0_adjust_factor) {
                div_factor = pcs_ptr->parent_pcs_ptr->used_tpl_frame_num *
                    pcs_ptr->parent_pcs_ptr->tpl_ctrls.r0_adjust_factor;
            }
            pcs_ptr->parent_pcs_ptr->r0 = pcs_ptr->parent_pcs_ptr->r0 / div_factor;
            pcs_ptr->parent_pcs_ptr->r0 = pcs_ptr->parent_pcs_ptr->r0 /
                tpl_hl_base_frame_div_factor[scs_ptr->max_heirachical_level];
        }

        int num_stats_required_for_gfu_boost = pcs_ptr->parent_pcs_ptr->tpl_group_size +
            (1 << pcs_ptr->parent_pcs_ptr->hierarchical_levels);

        rc->gfu_boost        = get_gfu_boost_from_r0_lap(min_boost_factor,
                                                  MAX_GFUBOOST_FACTOR,
                                                  pcs_ptr->parent_pcs_ptr->r0,
                                                  num_stats_required_for_gfu_boost);
        rc->arf_boost_factor = (pcs_ptr->ref_slice_type_array[0][0] == I_SLICE &&
                                pcs_ptr->ref_pic_r0[0][0] - pcs_ptr->parent_pcs_ptr->r0 >= 0.08)
            ? (float_t)1.3
            : (float_t)1;
        q                    = active_worst_quality;

        // non ref frame or repeated frames with re-encode
        if (!refresh_alt_ref_frame && !is_intrl_arf_boost)
            active_best_quality = cq_level;
        else {
            if (!is_intrl_arf_boost) {
                active_best_quality = get_gf_active_quality_tpl_la(rc, q, bit_depth);
                rc->arf_q           = active_best_quality;
                const int min_boost = get_gf_high_motion_quality(q, bit_depth);
                const int boost     = min_boost - active_best_quality;
                active_best_quality = min_boost - (int)(boost * rc->arf_boost_factor);
            } else {
                EbReferenceObject *ref_obj_l0 =
                    (EbReferenceObject *)pcs_ptr->ref_pic_ptr_array[REF_LIST_0][0]->object_ptr;
                EbReferenceObject *ref_obj_l1 = NULL;
                if (pcs_ptr->slice_type == B_SLICE)
                    ref_obj_l1 =
                        (EbReferenceObject *)pcs_ptr->ref_pic_ptr_array[REF_LIST_1][0]->object_ptr;

                uint8_t ref_tmp_layer = ref_obj_l0->tmp_layer_idx;
                if (pcs_ptr->slice_type == B_SLICE)
                    ref_tmp_layer = MAX(ref_tmp_layer, ref_obj_l1->tmp_layer_idx);
                active_best_quality    = rc->arf_q;
                int8_t tmp_layer_delta = (int8_t)pcs_ptr->parent_pcs_ptr->temporal_layer_index -
                    (int8_t)ref_tmp_layer;
                // active_best_quality is updated with the q index of the reference
                if (rf_level == GF_ARF_LOW) {
                    while (tmp_layer_delta--)
                        active_best_quality = (active_best_quality + cq_level + 1) / 2;
                }
            }
            // For alt_ref and GF frames (including internal arf frames) adjust the
            // worst allowed quality as well. This insures that even on hard
            // sections we dont clamp the Q at the same value for arf frames and
            // leaf (non arf) frames. This is important to the TPL model which assumes
            // Q drops with each arf level.
            active_worst_quality = (active_best_quality + (3 * active_worst_quality) + 2) / 4;
        }
    } else
        active_best_quality = cq_level;

    if (pcs_ptr->parent_pcs_ptr->temporal_layer_index)
        active_best_quality = MAX(active_best_quality, rc->arf_q);
    adjust_active_best_and_worst_quality(
        pcs_ptr, rc, rf_level, &active_worst_quality, &active_best_quality);
    q = active_best_quality;
    clamp(q, active_best_quality, active_worst_quality);
    pcs_ptr->parent_pcs_ptr->top_index    = active_worst_quality;
    pcs_ptr->parent_pcs_ptr->bottom_index = active_best_quality;
    assert(pcs_ptr->parent_pcs_ptr->top_index <= rc->worst_quality &&
           pcs_ptr->parent_pcs_ptr->top_index >= rc->best_quality);
    assert(pcs_ptr->parent_pcs_ptr->bottom_index <= rc->worst_quality &&
           pcs_ptr->parent_pcs_ptr->bottom_index >= rc->best_quality);
    return q;
}

#define DEFAULT_KF_BOOST 2700
#define DEFAULT_GF_BOOST 1350

int svt_av1_get_q_index_from_qstep_ratio(int leaf_qindex, double qstep_ratio, const int bit_depth) {
    const double leaf_qstep   = svt_av1_dc_quant_qtx(leaf_qindex, 0, bit_depth);
    const double target_qstep = leaf_qstep * qstep_ratio;
    int          qindex       = leaf_qindex;
    for (qindex = leaf_qindex; qindex > 0; --qindex) {
        const double qstep = svt_av1_dc_quant_qtx(qindex, 0, bit_depth);
        if (qstep + 0.1 <= target_qstep)
            break;
    }
    return qindex;
}
/******************************************************
 * non_base_boost
 * Compute a non-base frame boost.
 ******************************************************/
int8_t non_base_boost(PictureControlSet *pcs_ptr) {
    int8_t             q_boost = 0;
    EbReferenceObject *ref_obj_l0 =
        (EbReferenceObject *)pcs_ptr->ref_pic_ptr_array[REF_LIST_0][0]->object_ptr;
    uint32_t l0_was_intra = 0;
    if (ref_obj_l0->slice_type != I_SLICE) {
        for (uint32_t sb_index = 0; sb_index < pcs_ptr->sb_total_count_pix; sb_index++) {
            l0_was_intra += ref_obj_l0->sb_intra[sb_index];
        }
    }
    if (l0_was_intra) {
        int8_t intra_percentage = (l0_was_intra * 100) / pcs_ptr->sb_total_count_pix;
        q_boost                 = intra_percentage >> 2;
    }
    return q_boost;
}
/******************************************************
 * cqp_qindex_calc
 * Assign the q_index per frame.
 * Used in the one pass encoding with no look ahead
 ******************************************************/
static int cqp_qindex_calc(PictureControlSet *pcs_ptr, int qindex) {
    SequenceControlSet *scs_ptr = pcs_ptr->parent_pcs_ptr->scs_ptr;
    int                 q;
    const int           bit_depth = scs_ptr->static_config.encoder_bit_depth;

    int active_best_quality  = 0;
    int active_worst_quality = qindex;

    double q_val = svt_av1_convert_qindex_to_q(qindex, bit_depth);

    int offset_idx = -1;
    if (!pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag)
        offset_idx = -1;
    else if (pcs_ptr->parent_pcs_ptr->idr_flag)
        offset_idx = 0;
    else
        offset_idx = MIN(pcs_ptr->temporal_layer_index + 1, FIXED_QP_OFFSET_COUNT - 1);

    double q_val_target = (offset_idx == -1)
        ? q_val
        : MAX(q_val -
                  (q_val * percents[pcs_ptr->parent_pcs_ptr->hierarchical_levels <= 4][offset_idx] /
                   100),
              0.0);

    if (scs_ptr->static_config.pred_structure == EB_PRED_LOW_DELAY_P ||
        scs_ptr->static_config.pred_structure == EB_PRED_LOW_DELAY_B) {
        if (pcs_ptr->parent_pcs_ptr->temporal_layer_index) {
            int8_t boost = non_base_boost(pcs_ptr);
            if (boost)
                q_val_target = MAX(0, q_val_target - (boost * q_val_target) / 100);
        }
    }
    const int32_t delta_qindex = svt_av1_compute_qdelta(q_val, q_val_target, bit_depth);

    active_best_quality = (int32_t)(qindex + delta_qindex);
    q                   = active_best_quality;
    clamp(q, active_best_quality, active_worst_quality);

    return q;
}
const int64_t q_factor[2][6] = {{100, 110, 120, 138, 140, 150}, {100, 110, 112, 125, 135, 140}};
// The table we use is modified from libaom; here is the original, from libaom:
// static const int rd_frame_type_factor[FRAME_UPDATE_TYPES] = { 128, 144, 128,
//                                                               128, 144, 144,
//                                                               128 };
static const int rd_frame_type_factor[FRAME_UPDATE_TYPES] = {128, 164, 128, 128, 164, 164, 128};
/*
 * Set the sse lambda based on the bit_depth, then update based on frame position.
 */
int compute_rdmult_sse(PictureControlSet *pcs_ptr, uint8_t q_index, uint8_t bit_depth) {
    FrameType frame_type = pcs_ptr->parent_pcs_ptr->frm_hdr.frame_type;
    // To set gf_update_type based on current TL vs. the max TL (e.g. for 5L, max TL is 4)
    uint8_t temporal_layer_index = pcs_ptr->temporal_layer_index;
    uint8_t max_temporal_layer   = pcs_ptr->parent_pcs_ptr->hierarchical_levels;

    int64_t rdmult = bit_depth == 8 ? av1_lambda_mode_decision8_bit_sse[q_index]
        : bit_depth == 10           ? av1lambda_mode_decision10_bit_sse[q_index]
                                    : av1lambda_mode_decision12_bit_sse[q_index];

    // Update rdmult based on the frame's position in the miniGOP
    if (frame_type != KEY_FRAME) {
        uint8_t gf_update_type = temporal_layer_index == 0 ? ARF_UPDATE
            : temporal_layer_index < max_temporal_layer    ? INTNL_ARF_UPDATE
                                                           : LF_UPDATE;
        rdmult                 = (rdmult * rd_frame_type_factor[gf_update_type]) >> 7;
    }
    assert(pcs_ptr->parent_pcs_ptr->tpl_ctrls.vq_adjust_lambda_sb > 0);
    int    qdiff = q_index - quantizer_to_qindex[pcs_ptr->parent_pcs_ptr->picture_qp];
    int8_t qidx  = (qdiff < -8) ? 0 : (qdiff < -4) ? 1 : (qdiff < -2) ? 2 : (qdiff <= 4) ? 3 : 4;

    if (pcs_ptr->parent_pcs_ptr->tpl_ctrls.enable) {
        rdmult = (rdmult *
                  q_factor[pcs_ptr->parent_pcs_ptr->tpl_ctrls.vq_adjust_lambda_sb - 1][qidx]) >>
            7;
    }
    return (int)rdmult;
}
static void sb_setup_lambda(PictureControlSet *pcs_ptr, SuperBlock *sb_ptr) {
    const Av1Common *const   cm       = pcs_ptr->parent_pcs_ptr->av1_cm;
    PictureParentControlSet *ppcs_ptr = pcs_ptr->parent_pcs_ptr;
    SequenceControlSet      *scs_ptr  = ppcs_ptr->scs_ptr;

    int mi_col = sb_ptr->origin_x / 4;
    int mi_row = sb_ptr->origin_y / 4;

    const int mi_col_sr = coded_to_superres_mi(mi_col, ppcs_ptr->superres_denom);
    assert(ppcs_ptr->enhanced_unscaled_picture_ptr);
    // ALIGN_POWER_OF_TWO(pixels, 3) >> 2 ??
    const int mi_cols_sr     = ((ppcs_ptr->enhanced_unscaled_picture_ptr->width + 15) / 16) << 2;
    const int sb_mi_width_sr = coded_to_superres_mi(mi_size_wide[scs_ptr->seq_header.sb_size],
                                                    ppcs_ptr->superres_denom);
    const int bsize_base     = ppcs_ptr->tpl_ctrls.synth_blk_size == 32 ? BLOCK_32X32 : BLOCK_16X16;
    const int num_mi_w       = mi_size_wide[bsize_base];
    const int num_mi_h       = mi_size_high[bsize_base];
    const int num_cols       = (mi_cols_sr + num_mi_w - 1) / num_mi_w;
    const int num_rows       = (cm->mi_rows + num_mi_h - 1) / num_mi_h;
    const int num_bcols      = (sb_mi_width_sr + num_mi_w - 1) / num_mi_w;
    const int num_brows = (mi_size_high[scs_ptr->seq_header.sb_size] + num_mi_h - 1) / num_mi_h;

    int row, col;

    int32_t base_block_count = 0;
    double  log_sum          = 0.0;

    for (row = mi_row / num_mi_w; row < num_rows && row < mi_row / num_mi_w + num_brows; ++row) {
        for (col = mi_col_sr / num_mi_h; col < num_cols && col < mi_col_sr / num_mi_h + num_bcols;
             ++col) {
            const int index = row * num_cols + col;
            log_sum += log(ppcs_ptr->pa_me_data->tpl_rdmult_scaling_factors[index]);
            ++base_block_count;
        }
    }
    assert(base_block_count > 0);

    uint8_t bit_depth = pcs_ptr->hbd_mode_decision ? 10 : 8;

    const int orig_rdmult = compute_rdmult_sse(
        pcs_ptr, ppcs_ptr->frm_hdr.quantization_params.base_q_idx, bit_depth);
    const int    new_rdmult     = compute_rdmult_sse(pcs_ptr, sb_ptr->qindex, bit_depth);
    const double scaling_factor = (double)new_rdmult / (double)orig_rdmult;
    //double scale_adj = exp(log(scaling_factor) - log_sum / base_block_count);
    double scale_adj = scaling_factor / exp(log_sum / base_block_count);

    for (row = mi_row / num_mi_w; row < num_rows && row < mi_row / num_mi_w + num_brows; ++row) {
        for (col = mi_col_sr / num_mi_h; col < num_cols && col < mi_col_sr / num_mi_h + num_bcols;
             ++col) {
            const int index                                            = row * num_cols + col;
            ppcs_ptr->pa_me_data->tpl_sb_rdmult_scaling_factors[index] = scale_adj *
                ppcs_ptr->pa_me_data->tpl_rdmult_scaling_factors[index];
        }
    }
    ppcs_ptr->blk_lambda_tuning = EB_TRUE;
}

/******************************************************
 * sb_qp_derivation_tpl_la
 * Calculates the QP per SB based on the tpl statistics
 * used in one pass and second pass of two pass encoding
 ******************************************************/
void sb_qp_derivation_tpl_la(PictureControlSet *pcs_ptr) {
    PictureParentControlSet *ppcs_ptr = pcs_ptr->parent_pcs_ptr;
    SequenceControlSet      *scs_ptr  = pcs_ptr->parent_pcs_ptr->scs_ptr;
    SuperBlock              *sb_ptr;
    uint32_t                 sb_addr;
    if (pcs_ptr->temporal_layer_index == 0)
        pcs_ptr->parent_pcs_ptr->frm_hdr.delta_q_params.delta_q_present = 1;
    else
        pcs_ptr->parent_pcs_ptr->frm_hdr.delta_q_params.delta_q_present = 0;
    // super res pictures scaled with different sb count, should use sb_total_count for each picture
    uint16_t sb_cnt = scs_ptr->sb_tot_cnt;
    if (scs_ptr->static_config.superres_mode > SUPERRES_NONE)
        sb_cnt = ppcs_ptr->sb_total_count;
    if ((pcs_ptr->parent_pcs_ptr->frm_hdr.delta_q_params.delta_q_present) &&
        (pcs_ptr->parent_pcs_ptr->tpl_is_valid == 1)) {
        for (sb_addr = 0; sb_addr < sb_cnt; ++sb_addr) {
            sb_ptr        = pcs_ptr->sb_ptr_array[sb_addr];
            double beta   = ppcs_ptr->pa_me_data->tpl_beta[sb_addr];
            int    offset = svt_av1_get_deltaq_offset(scs_ptr->static_config.encoder_bit_depth,
                                                   ppcs_ptr->frm_hdr.quantization_params.base_q_idx,
                                                   beta,
                                                   pcs_ptr->parent_pcs_ptr->slice_type);
            offset        = AOMMIN(
                offset, pcs_ptr->parent_pcs_ptr->frm_hdr.delta_q_params.delta_q_res * 9 * 4 - 1);
            offset = AOMMAX(
                offset, -pcs_ptr->parent_pcs_ptr->frm_hdr.delta_q_params.delta_q_res * 9 * 4 + 1);
            sb_ptr->qindex = CLIP3(
                pcs_ptr->parent_pcs_ptr->frm_hdr.delta_q_params.delta_q_res,
                255 - pcs_ptr->parent_pcs_ptr->frm_hdr.delta_q_params.delta_q_res,
                ((int16_t)ppcs_ptr->frm_hdr.quantization_params.base_q_idx + (int16_t)offset));

            sb_setup_lambda(pcs_ptr, sb_ptr);
        }

    } else {
        for (sb_addr = 0; sb_addr < sb_cnt; ++sb_addr) {
            sb_ptr         = pcs_ptr->sb_ptr_array[sb_addr];
            sb_ptr->qindex = quantizer_to_qindex[pcs_ptr->picture_qp];
        }
    }
}

static int av1_find_qindex(double desired_q, aom_bit_depth_t bit_depth, int best_qindex,
                           int worst_qindex) {
    assert(best_qindex <= worst_qindex);
    int low  = best_qindex;
    int high = worst_qindex;
    while (low < high) {
        const int    mid   = (low + high) >> 1;
        const double mid_q = svt_av1_convert_qindex_to_q(mid, bit_depth);
        if (mid_q < desired_q) {
            low = mid + 1;
        } else {
            high = mid;
        }
    }
    assert(low == high);
    assert(svt_av1_convert_qindex_to_q(low, bit_depth) >= desired_q || low == worst_qindex);
    return low;
}
static int find_fp_qindex(aom_bit_depth_t bit_depth) {
#ifdef ARCH_X86_64
    aom_clear_system_state();
#endif
    return av1_find_qindex(FIRST_PASS_Q, bit_depth, 0, QINDEX_RANGE - 1);
}
int svt_av1_rc_get_default_min_gf_interval(int width, int height, double framerate) {
    // Assume we do not need any constraint lower than 4K 20 fps
    static const double factor_safe = 3840 * 2160 * 20.0;
    const double        factor      = width * height * framerate;
    const int default_interval = clamp((int)(framerate * 0.125), MIN_GF_INTERVAL, MAX_GF_INTERVAL);

    if (factor <= factor_safe)
        return default_interval;
    else
        return AOMMAX(default_interval, (int)(MIN_GF_INTERVAL * factor / factor_safe + 0.5));
    // Note this logic makes:
    // 4K24: 5
    // 4K30: 6
    // 4K60: 12
}
void set_rc_buffer_sizes(SequenceControlSet *scs_ptr) {
    EncodeContext        *encode_context_ptr = scs_ptr->encode_context_ptr;
    RATE_CONTROL         *rc                 = &encode_context_ptr->rc;
    RateControlCfg *const rc_cfg             = &encode_context_ptr->rc_cfg;
    const int64_t         bandwidth          = scs_ptr->static_config.target_bit_rate;
    const int64_t         starting           = rc_cfg->starting_buffer_level_ms;
    const int64_t         optimal            = rc_cfg->optimal_buffer_level_ms;
    const int64_t         maximum            = rc_cfg->maximum_buffer_size_ms;

    rc->starting_buffer_level = starting * bandwidth / 1000;
    rc->optimal_buffer_level  = (optimal == 0) ? bandwidth / 8 : optimal * bandwidth / 1000;
    rc->maximum_buffer_size   = (maximum == 0) ? bandwidth / 8 : maximum * bandwidth / 1000;
}

int svt_av1_rc_get_default_max_gf_interval(double framerate, int min_gf_interval) {
    int interval = AOMMIN(MAX_GF_INTERVAL, (int)(framerate * 0.75));
    interval += (interval & 0x01); // Round to even value
    interval = AOMMAX(MAX_GF_INTERVAL, interval);
    return AOMMAX(interval, min_gf_interval);
}

//#define INT_MAX 0x7fffffff
#define BPER_MB_NORMBITS 9
#define FRAME_OVERHEAD_BITS 200
static void av1_rc_init(SequenceControlSet *scs_ptr) {
    EncodeContext              *encode_context_ptr = scs_ptr->encode_context_ptr;
    RATE_CONTROL               *rc                 = &encode_context_ptr->rc;
    const RateControlCfg *const rc_cfg             = &encode_context_ptr->rc_cfg;
    const uint32_t              width              = scs_ptr->seq_header.max_frame_width;
    const uint32_t              height             = scs_ptr->seq_header.max_frame_height;
    int                         i;

    if ((!(scs_ptr->static_config.pass == ENC_MIDDLE_PASS ||
           scs_ptr->static_config.pass == ENC_LAST_PASS) &&
         scs_ptr->static_config.pass != ENC_FIRST_PASS && rc_cfg->mode == AOM_CBR)) {
        rc->avg_frame_qindex[KEY_FRAME]   = rc_cfg->worst_allowed_q;
        rc->avg_frame_qindex[INTER_FRAME] = rc_cfg->worst_allowed_q;
        rc->frames_till_gf_update_due     = 0;
    } else {
        rc->avg_frame_qindex[KEY_FRAME]   = (rc_cfg->worst_allowed_q + rc_cfg->best_allowed_q) / 2;
        rc->avg_frame_qindex[INTER_FRAME] = (rc_cfg->worst_allowed_q + rc_cfg->best_allowed_q) / 2;
    }
    rc->buffer_level    = rc->starting_buffer_level;
    rc->bits_off_target = rc->starting_buffer_level;

    rc->rolling_target_bits    = rc->avg_frame_bandwidth;
    rc->rolling_actual_bits    = rc->avg_frame_bandwidth;
    rc->total_actual_bits      = 0;
    rc->total_target_bits      = 0;
    rc->total_target_vs_actual = 0;

    rc->frames_since_key      = 8; // Sensible default for first frame.
    rc->this_key_frame_forced = 0;
    for (i = 0; i < MAX_TEMPORAL_LAYERS + 1; ++i) { rc->rate_correction_factors[i] = 0.7; }
    rc->rate_correction_factors[KF_STD] = 1.0;
    rc->min_gf_interval                 = encode_context_ptr->gf_cfg.min_gf_interval;
    rc->max_gf_interval                 = encode_context_ptr->gf_cfg.max_gf_interval;
    if (rc->min_gf_interval == 0)
        rc->min_gf_interval = svt_av1_rc_get_default_min_gf_interval(
            width,
            height, /*oxcf->frm_dim_cfg.width, oxcf->frm_dim_cfg.height,*/
            scs_ptr->double_frame_rate /*oxcf->input_cfg.init_framerate*/);
    if (rc->max_gf_interval == 0)
        rc->max_gf_interval = svt_av1_rc_get_default_max_gf_interval(
            scs_ptr->double_frame_rate /*oxcf->input_cfg.init_framerate*/, rc->min_gf_interval);
    rc->baseline_gf_interval = (rc->min_gf_interval + rc->max_gf_interval) / 2;
    //rc->avg_frame_low_motion = 0;

    // Set absolute upper and lower quality limits
    rc->worst_quality = rc_cfg->worst_allowed_q;
    rc->best_quality  = rc_cfg->best_allowed_q;
    if (!(scs_ptr->static_config.pass == ENC_MIDDLE_PASS ||
          scs_ptr->static_config.pass == ENC_LAST_PASS) &&
        scs_ptr->static_config.pass != ENC_FIRST_PASS && rc_cfg->mode == AOM_CBR) {
        if (scs_ptr->static_config.pred_structure == EB_PRED_LOW_DELAY_P) {
            rc->onepass_cbr_mode = 2;
        } else {
            rc->onepass_cbr_mode = 1;
        }
    } else {
        rc->onepass_cbr_mode = 0;
    }
    if (scs_ptr->lap_enabled ||
        (!(scs_ptr->static_config.pass == ENC_MIDDLE_PASS ||
           scs_ptr->static_config.pass == ENC_LAST_PASS) &&
         scs_ptr->static_config.pass != ENC_FIRST_PASS && rc_cfg->mode == AOM_CBR)) {
        double frame_rate = (double)scs_ptr->static_config.frame_rate_numerator /
            (double)scs_ptr->static_config.frame_rate_denominator;
        // Each frame can have a different duration, as the frame rate in the source
        // isn't guaranteed to be constant. The frame rate prior to the first frame
        // encoded in the second pass is a guess. However, the sum duration is not.
        // It is calculated based on the actual durations of all frames from the
        // first pass.
        svt_av1_new_framerate(scs_ptr, frame_rate);
    }
}

#define MIN_BOOST_COMBINE_FACTOR 4.0
#define MAX_BOOST_COMBINE_FACTOR 12.0
void process_tpl_stats_frame_kf_gfu_boost(PictureControlSet *pcs_ptr) {
    SequenceControlSet *scs_ptr            = pcs_ptr->parent_pcs_ptr->scs_ptr;
    EncodeContext      *encode_context_ptr = scs_ptr->encode_context_ptr;
    RATE_CONTROL *const rc                 = &encode_context_ptr->rc;
    // The new tpl only looks at pictures in tpl group, which is fewer than before,
    // As a results, we defined a factor to adjust r0
    if (pcs_ptr->parent_pcs_ptr->slice_type != 2) {
        double div_factor = 1;
        if (pcs_ptr->parent_pcs_ptr->tpl_ctrls.r0_adjust_factor) {
            div_factor = pcs_ptr->parent_pcs_ptr->used_tpl_frame_num *
                pcs_ptr->parent_pcs_ptr->tpl_ctrls.r0_adjust_factor;
        }

        pcs_ptr->parent_pcs_ptr->r0 = pcs_ptr->parent_pcs_ptr->r0 / div_factor;
    } else if (pcs_ptr->parent_pcs_ptr->frm_hdr.frame_type != KEY_FRAME) {
        double factor               = 2;
        pcs_ptr->parent_pcs_ptr->r0 = pcs_ptr->parent_pcs_ptr->r0 / factor;
    }
    rc->gfu_boost = get_gfu_boost_from_r0_lap(MIN_BOOST_COMBINE_FACTOR,
                                              MAX_GFUBOOST_FACTOR,
                                              pcs_ptr->parent_pcs_ptr->r0,
                                              rc->frames_to_key);
    if (pcs_ptr->parent_pcs_ptr->frm_hdr.frame_type == KEY_FRAME) {
        pcs_ptr->parent_pcs_ptr->r0 = pcs_ptr->parent_pcs_ptr->r0 /
            tpl_hl_islice_div_factor[scs_ptr->static_config.hierarchical_levels];
        if (pcs_ptr->parent_pcs_ptr->frm_hdr.frame_type == KEY_FRAME) {
            {
                double factor = 1.0;
                if (pcs_ptr->parent_pcs_ptr->r0 < 0.2) {
                    double mult = 1.0;
                    factor      = (double)(mult * 255.0) / rc->active_worst_quality;
                }
                pcs_ptr->parent_pcs_ptr->r0 = pcs_ptr->parent_pcs_ptr->r0 / factor;
            }
        }
        // when frames_to_key not available, i.e. in 1 pass encoding
        rc->kf_boost = get_cqp_kf_boost_from_r0(
            pcs_ptr->parent_pcs_ptr->r0, rc->frames_to_key, scs_ptr->input_resolution);
        int max_boost = pcs_ptr->parent_pcs_ptr->used_tpl_frame_num * KB;
        rc->kf_boost  = AOMMIN(rc->kf_boost, max_boost);
    }
}

static void get_intra_q_and_bounds(PictureControlSet *pcs_ptr, int *active_best, int *active_worst,
                                   int cq_level, int is_fwd_kf) {
    SequenceControlSet *scs_ptr            = pcs_ptr->parent_pcs_ptr->scs_ptr;
    EncodeContext      *encode_context_ptr = scs_ptr->encode_context_ptr;
    RATE_CONTROL       *rc                 = &encode_context_ptr->rc;
    TWO_PASS *const     twopass            = &scs_ptr->twopass;
    int                 active_best_quality;
    int                 active_worst_quality = *active_worst;
    const int           bit_depth            = scs_ptr->static_config.encoder_bit_depth;

    if (rc->frames_to_key <= 1 && encode_context_ptr->rc_cfg.mode == AOM_Q) {
        // If the next frame is also a key frame or the current frame is the
        // only frame in the sequence in AOM_Q mode, just use the cq_level
        // as q.
        active_best_quality  = cq_level;
        active_worst_quality = cq_level;
    } else if (is_fwd_kf) {
        // Handle the special case for forward reference key frames.
        // Increase the boost because this keyframe is used as a forward and
        // backward reference.
        const int    qindex         = rc->last_boosted_qindex;
        const double last_boosted_q = svt_av1_convert_qindex_to_q(qindex, bit_depth);
        const int    delta_qindex   = svt_av1_compute_qdelta(
            last_boosted_q, last_boosted_q * 0.25, bit_depth);
        active_best_quality = AOMMAX(qindex + delta_qindex, rc->best_quality);
    } else {
        // Not forced keyframe.
        double q_adj_factor = 1.0;
        double q_val;
        rc->worst_quality = MAXQ;
        rc->best_quality  = MINQ;

        // Baseline value derived from cpi->active_worst_quality and kf boost.
        active_best_quality = get_kf_active_quality_tpl(rc, active_worst_quality, bit_depth);
        if (twopass->kf_zeromotion_pct >= STATIC_KF_GROUP_THRESH) {
            active_best_quality /= 3;
        }
        if (pcs_ptr->parent_pcs_ptr->sc_class1 && scs_ptr->passes == 1 &&
            (encode_context_ptr->rc_cfg.mode == AOM_VBR ||
             encode_context_ptr->rc_cfg.mode == AOM_CBR))
            active_best_quality /= 2;
        // Allow somewhat lower kf minq with small image formats.
        if (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_240p_RANGE)
            q_adj_factor -= (pcs_ptr->parent_pcs_ptr->tune_tpl_for_chroma) ? 0.2 : 0.15;
        // Make a further adjustment based on the kf zero motion measure.

        // Convert the adjustment factor to a qindex delta
        // on active_best_quality.
        q_val = svt_av1_convert_qindex_to_q(active_best_quality, bit_depth);
        active_best_quality += svt_av1_compute_qdelta(q_val, q_val * q_adj_factor, bit_depth);

        // Tweak active_best_quality for AOM_Q mode when superres is on, as this
        // will be used directly as 'q' later.
        if (encode_context_ptr->rc_cfg.mode == AOM_Q &&
            (scs_ptr->static_config.superres_mode == SUPERRES_QTHRESH ||
             scs_ptr->static_config.superres_mode == SUPERRES_AUTO) &&
            pcs_ptr->parent_pcs_ptr->superres_denom != SCALE_NUMERATOR) {
            active_best_quality = AOMMAX(
                active_best_quality -
                    ((pcs_ptr->parent_pcs_ptr->superres_denom - SCALE_NUMERATOR) *
                     SUPERRES_QADJ_PER_DENOM_KEYFRAME),
                0);
        }
    }

    *active_best  = active_best_quality;
    *active_worst = active_worst_quality;
    return;
}

// Returns |active_best_quality| for an inter frame.
// The |active_best_quality| depends on different rate control modes:
// VBR, Q, CQ, CBR.
// The returning active_best_quality could further be adjusted in
// adjust_active_best_and_worst_quality().
static int get_active_best_quality(PictureControlSet *pcs_ptr, const int active_worst_quality,
                                   const int cq_level) {
    SequenceControlSet    *scs_ptr            = pcs_ptr->parent_pcs_ptr->scs_ptr;
    EncodeContext         *encode_context_ptr = scs_ptr->encode_context_ptr;
    RATE_CONTROL          *rc                 = &encode_context_ptr->rc;
    const enum aom_rc_mode rc_mode            = encode_context_ptr->rc_cfg.mode;
    GF_GROUP *const        gf_group           = &encode_context_ptr->gf_group;
    const int              bit_depth          = scs_ptr->static_config.encoder_bit_depth;
    //int rf_level, update_type;
    const RefreshFrameFlagsInfo *const refresh_frame_flags =
        &pcs_ptr->parent_pcs_ptr->refresh_frame;
    const int is_intrl_arf_boost = gf_group->update_type[pcs_ptr->parent_pcs_ptr->gf_group_index] ==
        INTNL_ARF_UPDATE;
    int *inter_minq;
    ASSIGN_MINQ_TABLE(bit_depth, inter_minq);
    int       active_best_quality = 0;
    const int is_leaf_frame       = !(refresh_frame_flags->golden_frame ||
                                refresh_frame_flags->alt_ref_frame || is_intrl_arf_boost);
    const int is_overlay_frame    = pcs_ptr->parent_pcs_ptr->is_overlay; //rc->is_src_frame_alt_ref;

    if (is_leaf_frame || is_overlay_frame) {
        if (rc_mode == AOM_Q)
            return cq_level;

        active_best_quality = inter_minq[active_worst_quality];
        // For the constrained quality mode we don't want
        // q to fall below the cq level.
        if ((rc_mode == AOM_CQ) && (active_best_quality < cq_level)) {
            active_best_quality = cq_level;
        }
        return active_best_quality;
    }
    // Determine active_best_quality for frames that are not leaf or overlay.
    int q = active_worst_quality;
    // Use the lower of active_worst_quality and recent
    // average Q as basis for GF/ARF best Q limit unless last frame was
    // a key frame.
    if ((rc_mode == AOM_VBR || rc_mode == AOM_CBR) && rc->frames_since_key > 1 &&
        rc->avg_frame_qindex[INTER_FRAME] < active_worst_quality) {
        q = rc->avg_frame_qindex[INTER_FRAME];
    }
    if (rc_mode == AOM_CQ && q < cq_level)
        q = cq_level;
    active_best_quality = get_gf_active_quality_tpl_la(rc, q, bit_depth);
    // Constrained quality use slightly lower active best.
    if (rc_mode == AOM_CQ)
        active_best_quality = active_best_quality * 15 / 16;
    const int min_boost = get_gf_high_motion_quality(q, bit_depth);
    const int boost     = min_boost - active_best_quality;

    rc->arf_boost_factor = (pcs_ptr->ref_slice_type_array[0][0] == I_SLICE &&
                            pcs_ptr->ref_pic_r0[0][0] - pcs_ptr->parent_pcs_ptr->r0 >= 0.08)
        ? (float_t)1.3
        : (float_t)1;
    active_best_quality  = min_boost - (int)(boost * rc->arf_boost_factor);
    if (!is_intrl_arf_boost)
        return active_best_quality;

    if (rc_mode == AOM_Q || rc_mode == AOM_CQ)
        active_best_quality = rc->arf_q;
    int this_height =
        gf_group->layer_depth[pcs_ptr->parent_pcs_ptr
                                  ->gf_group_index]; //gf_group_pyramid_level(gf_group, gf_index);
    while (this_height > 1) {
        active_best_quality = (active_best_quality + active_worst_quality + 1) / 2;
        --this_height;
    }
    return active_best_quality;
}

static double get_rate_correction_factor(PictureParentControlSet *ppcs_ptr/*,
                                         int width, int height*/) {
    SequenceControlSet                *scs_ptr             = ppcs_ptr->scs_ptr;
    EncodeContext                     *encode_context_ptr  = scs_ptr->encode_context_ptr;
    RATE_CONTROL                      *rc                  = &encode_context_ptr->rc;
    const RefreshFrameFlagsInfo *const refresh_frame_flags = &ppcs_ptr->refresh_frame;
    double                             rcf;

    if (ppcs_ptr->frm_hdr.frame_type == KEY_FRAME) {
        rcf = rc->rate_correction_factors[KF_STD];
    } else if (scs_ptr->static_config.pass == ENC_MIDDLE_PASS ||
               scs_ptr->static_config.pass == ENC_LAST_PASS || scs_ptr->lap_enabled) {
        const rate_factor_level rf_lvl = ppcs_ptr->temporal_layer_index + 1;
        rcf                            = rc->rate_correction_factors[rf_lvl];
    } else {
        if ((refresh_frame_flags->alt_ref_frame || refresh_frame_flags->golden_frame) &&
            !rc->is_src_frame_alt_ref && //!cpi->use_svc &&
            (encode_context_ptr->rc_cfg.mode != AOM_CBR ||
             encode_context_ptr->rc_cfg.gf_cbr_boost_pct > 20))
            rcf = rc->rate_correction_factors[GF_ARF_STD];
        else
            rcf = rc->rate_correction_factors[INTER_NORMAL];
    }
    //rcf *= resize_rate_factor(&cpi->oxcf.frm_dim_cfg, width, height);
    return fclamp(rcf, MIN_BPB_FACTOR, MAX_BPB_FACTOR);
}

static void set_rate_correction_factor(PictureParentControlSet *ppcs_ptr, double factor/*,
                                       int width, int height*/) {
    SequenceControlSet                *scs_ptr             = ppcs_ptr->scs_ptr;
    EncodeContext                     *encode_context_ptr  = scs_ptr->encode_context_ptr;
    RATE_CONTROL                      *rc                  = &encode_context_ptr->rc;
    const RefreshFrameFlagsInfo *const refresh_frame_flags = &ppcs_ptr->refresh_frame;

    // Normalize RCF to account for the size-dependent scaling factor.
    //factor /= resize_rate_factor(&cpi->oxcf.frm_dim_cfg, width, height);

    factor = fclamp(factor, MIN_BPB_FACTOR, MAX_BPB_FACTOR);

    if (ppcs_ptr->frm_hdr.frame_type == KEY_FRAME) {
        rc->rate_correction_factors[KF_STD] = factor;
    } else if (scs_ptr->static_config.pass == ENC_MIDDLE_PASS ||
               scs_ptr->static_config.pass == ENC_LAST_PASS || scs_ptr->lap_enabled) {
        const rate_factor_level rf_lvl      = ppcs_ptr->temporal_layer_index + 1;
        rc->rate_correction_factors[rf_lvl] = factor;
    } else {
        if ((refresh_frame_flags->alt_ref_frame || refresh_frame_flags->golden_frame) &&
            !rc->is_src_frame_alt_ref && //!cpi->use_svc &&
            (encode_context_ptr->rc_cfg.mode != AOM_CBR ||
             encode_context_ptr->rc_cfg.gf_cbr_boost_pct > 20))
            rc->rate_correction_factors[GF_ARF_STD] = factor;
        else
            rc->rate_correction_factors[INTER_NORMAL] = factor;
    }
}

// Calculate rate for the given 'q'.
static int get_bits_per_mb(PictureParentControlSet *ppcs_ptr, int use_cyclic_refresh,
                           double correction_factor, int q) {
    SequenceControlSet *scs_ptr = ppcs_ptr->scs_ptr;
    return use_cyclic_refresh
        ? 0 /*av1_cyclic_refresh_rc_bits_per_mb(cpi, q, correction_factor)*/
        : svt_av1_rc_bits_per_mb(ppcs_ptr->frm_hdr.frame_type,
                                 q,
                                 correction_factor,
                                 scs_ptr->static_config.encoder_bit_depth,
                                 ppcs_ptr->sc_class1,
                                 scs_ptr->encode_context_ptr->rc.onepass_cbr_mode);
}
//anaghdin --> update this if desired is less than cur_bit_diff by far
// Similar to find_qindex_by_rate() function in ratectrl.c, but returns the q
// index with rate just above or below the desired rate, depending on which of
// the two rates is closer to the desired rate.
// Also, respects the selected aq_mode when computing the rate.
static int find_closest_qindex_by_rate(int desired_bits_per_mb, PictureParentControlSet *ppcs_ptr,
                                       double correction_factor, int best_qindex,
                                       int worst_qindex) {
    const int use_cyclic_refresh = 0 /*cpi->oxcf.q_cfg.aq_mode == CYCLIC_REFRESH_AQ &&
                                 cpi->cyclic_refresh->apply_cyclic_refresh*/
        ;

    // Find 'qindex' based on 'desired_bits_per_mb'.
    assert(best_qindex <= worst_qindex);
    int low  = best_qindex;
    int high = worst_qindex;
    while (low < high) {
        const int mid             = (low + high) >> 1;
        const int mid_bits_per_mb = get_bits_per_mb(
            ppcs_ptr, use_cyclic_refresh, correction_factor, mid);
        if (mid_bits_per_mb > desired_bits_per_mb) {
            low = mid + 1;
        } else {
            high = mid;
        }
    }
    assert(low == high);

    // Calculate rate difference of this q index from the desired rate.
    const int curr_q           = low;
    const int curr_bits_per_mb = get_bits_per_mb(
        ppcs_ptr, use_cyclic_refresh, correction_factor, curr_q);
    const int curr_bit_diff = (curr_bits_per_mb <= desired_bits_per_mb)
        ? desired_bits_per_mb - curr_bits_per_mb
        : INT_MAX;
    assert((curr_bit_diff != INT_MAX && curr_bit_diff >= 0) || curr_q == worst_qindex);

    // Calculate rate difference for previous q index too.
    const int prev_q = curr_q - 1;
    int       prev_bit_diff;
    if (curr_bit_diff == INT_MAX || curr_q == best_qindex) {
        prev_bit_diff = INT_MAX;
    } else {
        const int prev_bits_per_mb = get_bits_per_mb(
            ppcs_ptr, use_cyclic_refresh, correction_factor, prev_q);
        assert(prev_bits_per_mb > desired_bits_per_mb);
        prev_bit_diff = prev_bits_per_mb - desired_bits_per_mb;
    }

    // Pick one of the two q indices, depending on which one has rate closer to
    // the desired rate.
    return (curr_bit_diff <= prev_bit_diff) ? curr_q : prev_q;
}

static int adjust_q_cbr(PictureParentControlSet *ppcs_ptr, int q /*, int active_worst_quality*/) {
    SequenceControlSet *scs_ptr            = ppcs_ptr->scs_ptr;
    EncodeContext      *encode_context_ptr = scs_ptr->encode_context_ptr;
    RATE_CONTROL       *rc                 = &encode_context_ptr->rc;
    //const AV1_COMMON *const cm = &cpi->common;
    const RefreshFrameFlagsInfo *const refresh_frame_flags        = &ppcs_ptr->refresh_frame;
    const int                          max_delta                  = 16;
    const int                          change_avg_frame_bandwidth = abs(rc->avg_frame_bandwidth -
                                               rc->prev_avg_frame_bandwidth) >
        0.1 * (rc->avg_frame_bandwidth);
    // If resolution changes or avg_frame_bandwidth significantly changed,
    // then set this flag to indicate change in target bits per macroblock.
    const int change_target_bits_mb =
        /*cm->prev_frame &&
      (cm->width != cm->prev_frame->width ||
       cm->height != cm->prev_frame->height || */
        change_avg_frame_bandwidth;
    // Apply some control/clamp to QP under certain conditions.
    if (ppcs_ptr->frm_hdr.frame_type != KEY_FRAME && /*!cpi->use_svc &&*/
        rc->frames_since_key > 1 && !change_target_bits_mb &&
        (!encode_context_ptr->rc_cfg.gf_cbr_boost_pct ||
         !(refresh_frame_flags->alt_ref_frame || refresh_frame_flags->golden_frame))) {
        // Make sure q is between oscillating Qs to prevent resonance.
        if (rc->rc_1_frame * rc->rc_2_frame == -1 && rc->q_1_frame != rc->q_2_frame) {
            q = clamp(
                q, AOMMIN(rc->q_1_frame, rc->q_2_frame), AOMMAX(rc->q_1_frame, rc->q_2_frame));
        }
        // Limit the decrease in Q from previous frame.
        if (rc->q_1_frame - q > max_delta)
            q = rc->q_1_frame - max_delta;
    }

    return AOMMAX(AOMMIN(q, rc->worst_quality), rc->best_quality);
}

static int av1_rc_regulate_q(PictureParentControlSet *ppcs_ptr, int target_bits_per_frame,
                             int active_best_quality, int active_worst_quality, int width,
                             int height) {
    const int    MBs = ((width + 15) / 16) * ((height + 15) / 16); //av1_get_MBs(width, height);
    const double correction_factor  = get_rate_correction_factor(ppcs_ptr /*, width, height*/);
    const int    target_bits_per_mb = (int)(((uint64_t)target_bits_per_frame << BPER_MB_NORMBITS) /
                                         MBs);

    int q = find_closest_qindex_by_rate(
        target_bits_per_mb, ppcs_ptr, correction_factor, active_best_quality, active_worst_quality);
    SequenceControlSet *scs_ptr            = ppcs_ptr->scs_ptr;
    EncodeContext      *encode_context_ptr = scs_ptr->encode_context_ptr;
    if (!(scs_ptr->static_config.pass == ENC_MIDDLE_PASS ||
          scs_ptr->static_config.pass == ENC_LAST_PASS) &&
        scs_ptr->static_config.pass != ENC_FIRST_PASS &&
        encode_context_ptr->rc_cfg.mode == AOM_CBR) {
        return adjust_q_cbr(ppcs_ptr, q /*, active_worst_quality*/);
    }

    return q;
}

static int get_q(PictureControlSet *pcs_ptr, const int active_worst_quality,
                 const int active_best_quality) {
    SequenceControlSet    *scs_ptr            = pcs_ptr->parent_pcs_ptr->scs_ptr;
    EncodeContext         *encode_context_ptr = scs_ptr->encode_context_ptr;
    RATE_CONTROL          *rc                 = &encode_context_ptr->rc;
    TWO_PASS *const        twopass            = &scs_ptr->twopass;
    const enum aom_rc_mode rc_mode            = encode_context_ptr->rc_cfg.mode;
    const int              width  = pcs_ptr->parent_pcs_ptr->av1_cm->frm_size.frame_width;
    const int              height = pcs_ptr->parent_pcs_ptr->av1_cm->frm_size.frame_height;
    int                    q;
    if (rc_mode == AOM_Q ||
        (frame_is_intra_only(pcs_ptr->parent_pcs_ptr) &&
         twopass->kf_zeromotion_pct >= STATIC_KF_GROUP_THRESH && rc->frames_to_key > 1)) {
        q = active_best_quality;
    } else {
        q = av1_rc_regulate_q(pcs_ptr->parent_pcs_ptr,
                              pcs_ptr->parent_pcs_ptr->this_frame_target,
                              active_best_quality,
                              active_worst_quality,
                              width,
                              height);
        if (q > active_worst_quality) {
            // Special case when we are targeting the max allowed rate.
            if (pcs_ptr->parent_pcs_ptr->this_frame_target < rc->max_frame_bandwidth) {
                q = active_worst_quality;
            }
        }
        q = AOMMAX(q, active_best_quality);
    }
    return q;
}

// Adjust active_worst_quality level based on buffer level.
static int calc_active_worst_quality_no_stats_cbr(PictureParentControlSet *ppcs_ptr) {
    // Adjust active_worst_quality: If buffer is above the optimal/target level,
    // bring active_worst_quality down depending on fullness of buffer.
    // If buffer is below the optimal level, let the active_worst_quality go from
    // ambient Q (at buffer = optimal level) to worst_quality level
    // (at buffer = critical level).
    SequenceControlSet *scs_ptr            = ppcs_ptr->scs_ptr;
    EncodeContext      *encode_context_ptr = scs_ptr->encode_context_ptr;
    RATE_CONTROL       *rc                 = &encode_context_ptr->rc;
    // Buffer level below which we push active_worst to worst_quality.
    int64_t critical_level = rc->optimal_buffer_level >> 3;
    int64_t buff_lvl_step  = 0;
    int     adjustment     = 0;
    int     active_worst_quality;
    int     ambient_qp;
    if (ppcs_ptr->frm_hdr.frame_type == KEY_FRAME)
        return rc->worst_quality;
    // For ambient_qp we use minimum of avg_frame_qindex[KEY_FRAME/INTER_FRAME]
    // for the first few frames following key frame. These are both initialized
    // to worst_quality and updated with (3/4, 1/4) average in postencode_update.
    // So for first few frames following key, the qp of that key frame is weighted
    // into the active_worst_quality setting.
    int32_t frame_updated = 0;
    svt_block_on_mutex(encode_context_ptr->frame_updated_mutex);
    frame_updated = encode_context_ptr->frame_updated;
    svt_release_mutex(encode_context_ptr->frame_updated_mutex);
    ambient_qp           = (frame_updated < 4)
                  ? AOMMIN(rc->avg_frame_qindex[INTER_FRAME], rc->avg_frame_qindex[KEY_FRAME])
                  : rc->avg_frame_qindex[INTER_FRAME];
    active_worst_quality = AOMMIN(rc->worst_quality, ambient_qp * 5 / 4);
    if (rc->buffer_level > rc->optimal_buffer_level) {
        // Adjust down.
        // Maximum limit for down adjustment, ~30%.
        int max_adjustment_down = active_worst_quality / 3;
        if (max_adjustment_down) {
            buff_lvl_step = ((rc->maximum_buffer_size - rc->optimal_buffer_level) /
                             max_adjustment_down);
            if (buff_lvl_step)
                adjustment = (int)((rc->buffer_level - rc->optimal_buffer_level) / buff_lvl_step);
            active_worst_quality -= adjustment;
        }
    } else if (rc->buffer_level > critical_level) {
        // Adjust up from ambient Q.
        if (critical_level) {
            buff_lvl_step = (rc->optimal_buffer_level - critical_level);
            if (buff_lvl_step) {
                adjustment = (int)((rc->worst_quality - ambient_qp) *
                                   (rc->optimal_buffer_level - rc->buffer_level) / buff_lvl_step);
            }
            active_worst_quality = ambient_qp + adjustment;
        }
    } else {
        // Set to worst_quality if buffer is below critical level.
        active_worst_quality = rc->worst_quality;
    }
    return active_worst_quality;
}

// Calculate the active_best_quality level.
static int calc_active_best_quality_no_stats_cbr(PictureParentControlSet *ppcs_ptr,
                                                 int active_worst_quality, int width, int height) {
    SequenceControlSet                *scs_ptr             = ppcs_ptr->scs_ptr;
    EncodeContext                     *encode_context_ptr  = scs_ptr->encode_context_ptr;
    RATE_CONTROL                      *rc                  = &encode_context_ptr->rc;
    const RefreshFrameFlagsInfo *const refresh_frame_flags = &ppcs_ptr->refresh_frame;
    int                               *rtc_minq;
    const int                          bit_depth = scs_ptr->static_config.encoder_bit_depth;
    int                                active_best_quality = rc->best_quality;
    ASSIGN_MINQ_TABLE(bit_depth, rtc_minq);

    if (frame_is_intra_only(ppcs_ptr)) {
        if (ppcs_ptr->frame_offset > 0) {
            // not first frame of one pass and kf_boost is set
            double q_adj_factor = 1.0;
            double q_val;
            active_best_quality =
                //get_kf_active_quality_cqp(rc, rc->avg_frame_qindex[KEY_FRAME], bit_depth);
                get_kf_active_quality_tpl(rc, rc->avg_frame_qindex[KEY_FRAME], bit_depth);
            // Allow somewhat lower kf minq with small image formats.
            if ((width * height) <= (352 * 288)) {
                q_adj_factor -= 0.25;
            }
            // Convert the adjustment factor to a qindex delta
            // on active_best_quality.
            q_val = svt_av1_convert_qindex_to_q(active_best_quality, bit_depth);
            //active_best_quality +=
            //    av1_compute_qdelta(rc, q_val, q_val * q_adj_factor, bit_depth);
            active_best_quality += svt_av1_compute_qdelta(q_val, q_val * q_adj_factor, bit_depth);
        }
    } else if (!rc->is_src_frame_alt_ref /*&& !cpi->use_svc */ &&
               encode_context_ptr->rc_cfg.gf_cbr_boost_pct &&
               (refresh_frame_flags->golden_frame || refresh_frame_flags->alt_ref_frame)) {
        // Use the lower of active_worst_quality and recent
        // average Q as basis for GF/ARF best Q limit unless last frame was
        // a key frame.
        int q = active_worst_quality;
        if (rc->frames_since_key > 1 && rc->avg_frame_qindex[INTER_FRAME] < active_worst_quality) {
            q = rc->avg_frame_qindex[INTER_FRAME];
        }
        //active_best_quality = get_gf_active_quality(rc, q, bit_depth);
        active_best_quality = get_gf_active_quality_tpl_la(rc, q, bit_depth);
    } else {
        // Use the lower of active_worst_quality and recent/average Q.
        svt_block_on_mutex(encode_context_ptr->frame_updated_mutex);
        FRAME_TYPE frame_type = (encode_context_ptr->frame_updated > 0) ? INTER_FRAME : KEY_FRAME;
        svt_release_mutex(encode_context_ptr->frame_updated_mutex);
        if (rc->avg_frame_qindex[frame_type] < active_worst_quality)
            active_best_quality = rtc_minq[rc->avg_frame_qindex[frame_type]];
        else
            active_best_quality = rtc_minq[active_worst_quality];
    }
    return active_best_quality;
}

static int rc_pick_q_and_bounds_no_stats_cbr(PictureControlSet *pcs_ptr) {
    SequenceControlSet *scs_ptr            = pcs_ptr->parent_pcs_ptr->scs_ptr;
    EncodeContext      *encode_context_ptr = scs_ptr->encode_context_ptr;
    RATE_CONTROL       *rc                 = &encode_context_ptr->rc;
    GF_GROUP *const     gf_group           = &encode_context_ptr->gf_group;
    int                 q;
    const int           bit_depth = scs_ptr->static_config.encoder_bit_depth;
    const int           width     = pcs_ptr->parent_pcs_ptr->av1_cm->frm_size.frame_width;
    const int           height    = pcs_ptr->parent_pcs_ptr->av1_cm->frm_size.frame_height;
    int active_worst_quality      = calc_active_worst_quality_no_stats_cbr(pcs_ptr->parent_pcs_ptr);
    int active_best_quality       = calc_active_best_quality_no_stats_cbr(
        pcs_ptr->parent_pcs_ptr, active_worst_quality, width, height);
    assert(encode_context_ptr->rc_cfg.mode == AOM_CBR);

    // Clip the active best and worst quality values to limits
    active_best_quality  = clamp(active_best_quality, rc->best_quality, rc->worst_quality);
    active_worst_quality = clamp(active_worst_quality, active_best_quality, rc->worst_quality);

    pcs_ptr->parent_pcs_ptr->top_index    = active_worst_quality;
    pcs_ptr->parent_pcs_ptr->bottom_index = active_best_quality;

    // Limit Q range for the adaptive loop.
    if (pcs_ptr->parent_pcs_ptr->frm_hdr.frame_type == KEY_FRAME && !rc->this_key_frame_forced &&
        pcs_ptr->parent_pcs_ptr->frame_offset != 0) {
        int qdelta = 0;
#ifdef ARCH_X86_64
        aom_clear_system_state();
#endif
        qdelta                             = svt_av1_compute_qdelta_by_rate(rc,
                                                pcs_ptr->parent_pcs_ptr->frm_hdr.frame_type,
                                                active_worst_quality,
                                                2.0,
                                                bit_depth,
                                                pcs_ptr->parent_pcs_ptr->sc_class1);
        pcs_ptr->parent_pcs_ptr->top_index = active_worst_quality + qdelta;
        pcs_ptr->parent_pcs_ptr->top_index = AOMMAX(pcs_ptr->parent_pcs_ptr->top_index,
                                                    pcs_ptr->parent_pcs_ptr->bottom_index);
    }

    // Special case code to try and match quality with forced key frames
    if (pcs_ptr->parent_pcs_ptr->frm_hdr.frame_type == KEY_FRAME && rc->this_key_frame_forced) {
        q = rc->last_boosted_qindex;
    } else {
        q = av1_rc_regulate_q(pcs_ptr->parent_pcs_ptr,
                              pcs_ptr->parent_pcs_ptr->this_frame_target,
                              active_best_quality,
                              active_worst_quality,
                              width,
                              height);
        if (q > pcs_ptr->parent_pcs_ptr->top_index) {
            // Special case when we are targeting the max allowed rate
            if (pcs_ptr->parent_pcs_ptr->this_frame_target >= rc->max_frame_bandwidth)
                pcs_ptr->parent_pcs_ptr->top_index = q;
            else
                q = pcs_ptr->parent_pcs_ptr->top_index;
        }
    }
    assert(pcs_ptr->parent_pcs_ptr->top_index <= rc->worst_quality &&
           pcs_ptr->parent_pcs_ptr->top_index >= rc->best_quality);
    assert(pcs_ptr->parent_pcs_ptr->bottom_index <= rc->worst_quality &&
           pcs_ptr->parent_pcs_ptr->bottom_index >= rc->best_quality);
    assert(q <= rc->worst_quality && q >= rc->best_quality);

    if (gf_group->update_type[pcs_ptr->parent_pcs_ptr->gf_group_index] == ARF_UPDATE)
        rc->arf_q = q;
    return q;
}

/******************************************************
 * rc_pick_q_and_bounds
 * assigns the q_index per frame using first pass statistics per frame.
 * used in the second pass of two pass encoding
 ******************************************************/
static int rc_pick_q_and_bounds(PictureControlSet *pcs_ptr) {
    SequenceControlSet                *scs_ptr            = pcs_ptr->parent_pcs_ptr->scs_ptr;
    EncodeContext                     *encode_context_ptr = scs_ptr->encode_context_ptr;
    RATE_CONTROL                      *rc                 = &encode_context_ptr->rc;
    GF_GROUP *const                    gf_group           = &encode_context_ptr->gf_group;
    const RefreshFrameFlagsInfo *const refresh_frame_flags =
        &pcs_ptr->parent_pcs_ptr->refresh_frame;
    const enum aom_rc_mode rc_mode              = encode_context_ptr->rc_cfg.mode;
    const int              cq_level             = encode_context_ptr->rc_cfg.cq_level;
    int                    active_best_quality  = 0;
    int                    active_worst_quality = rc->active_worst_quality;
    int                    q;
    int is_intrl_arf_boost = gf_group->update_type[pcs_ptr->parent_pcs_ptr->gf_group_index] ==
        INTNL_ARF_UPDATE;

    if (frame_is_intra_only(pcs_ptr->parent_pcs_ptr)) {
        const int is_fwd_kf = pcs_ptr->parent_pcs_ptr->frm_hdr.frame_type == KEY_FRAME &&
            pcs_ptr->parent_pcs_ptr->frm_hdr.show_frame == 0;
        get_intra_q_and_bounds(
            pcs_ptr, &active_best_quality, &active_worst_quality, cq_level, is_fwd_kf);
    } else {
        const int pyramid_level = gf_group->layer_depth[pcs_ptr->parent_pcs_ptr->gf_group_index];
        if ((pyramid_level <= 1) || (pyramid_level > MAX_ARF_LAYERS) || rc_mode == AOM_Q) {
            active_best_quality = get_active_best_quality(pcs_ptr, active_worst_quality, cq_level);
        } else {
            active_best_quality = rc->active_best_quality[pyramid_level - 1] + 1;
            active_best_quality = AOMMIN(active_best_quality, active_worst_quality);
            active_best_quality += (active_worst_quality - active_best_quality) / 2;
        }
        // For alt_ref and GF frames (including internal arf frames) adjust the
        // worst allowed quality as well. This insures that even on hard
        // sections we dont clamp the Q at the same value for arf frames and
        // leaf (non arf) frames. This is important to the TPL model which assumes
        // Q drops with each arf level.
        if (!(rc->is_src_frame_alt_ref) &&
            (refresh_frame_flags->golden_frame || refresh_frame_flags->alt_ref_frame ||
             is_intrl_arf_boost)) {
            active_worst_quality = (active_best_quality + (3 * active_worst_quality) + 2) / 4;
        }
    }

    adjust_active_best_and_worst_quality_org(
        pcs_ptr, rc, &active_worst_quality, &active_best_quality);

    q = get_q(pcs_ptr, active_worst_quality, active_best_quality);

    // Special case when we are targeting the max allowed rate.
    if (pcs_ptr->parent_pcs_ptr->this_frame_target >= rc->max_frame_bandwidth &&
        q > active_worst_quality) {
        active_worst_quality = q;
    }
    pcs_ptr->parent_pcs_ptr->top_index    = active_worst_quality;
    pcs_ptr->parent_pcs_ptr->bottom_index = active_best_quality;
    assert(pcs_ptr->parent_pcs_ptr->top_index <= rc->worst_quality &&
           pcs_ptr->parent_pcs_ptr->top_index >= rc->best_quality);
    assert(pcs_ptr->parent_pcs_ptr->bottom_index <= rc->worst_quality &&
           pcs_ptr->parent_pcs_ptr->bottom_index >= rc->best_quality);

    assert(q <= rc->worst_quality && q >= rc->best_quality);

    if (gf_group->update_type[pcs_ptr->parent_pcs_ptr->gf_group_index] == ARF_UPDATE)
        rc->arf_q = q;

    return q;
}

static int av1_estimate_bits_at_q(FrameType frame_type, int q, int mbs, double correction_factor,
                                  AomBitDepth bit_depth, uint8_t sc_content_detected,
                                  int onepass_cbr_mode) {
    const int bpm = (int)(svt_av1_rc_bits_per_mb(
        frame_type, q, correction_factor, bit_depth, sc_content_detected, onepass_cbr_mode));
    return AOMMAX(FRAME_OVERHEAD_BITS, (int)((uint64_t)bpm * mbs) >> BPER_MB_NORMBITS);
}

static void av1_rc_update_rate_correction_factors(PictureParentControlSet *ppcs_ptr, int width,
                                                  int height) {
    SequenceControlSet *scs_ptr            = ppcs_ptr->scs_ptr;
    EncodeContext      *encode_context_ptr = scs_ptr->encode_context_ptr;
    RATE_CONTROL       *rc                 = &encode_context_ptr->rc;
    int                 correction_factor  = 100;
    double rate_correction_factor = get_rate_correction_factor(ppcs_ptr /*, width, height*/);
    double adjustment_limit;
    //const int MBs = av1_get_MBs(width, height);
    const int MBs = ((width + 15) / 16) * ((height + 15) / 16); //av1_get_MBs(width, height);

    int projected_size_based_on_q = 0;

    // Do not update the rate factors for arf overlay frames.
    if (rc->is_src_frame_alt_ref)
        return;

    // Clear down mmx registers to allow floating point in what follows
    //aom_clear_system_state();

    // Work out how big we would have expected the frame to be at this Q given
    // the current correction factor.
    // Stay in double to avoid int overflow when values are large

    {
        projected_size_based_on_q = av1_estimate_bits_at_q(
            ppcs_ptr->frm_hdr.frame_type,
            ppcs_ptr->frm_hdr.quantization_params.base_q_idx /*cm->quant_params.base_qindex*/,
            MBs,
            rate_correction_factor,
            scs_ptr->static_config.encoder_bit_depth,
            ppcs_ptr->sc_class1,
            rc->onepass_cbr_mode);
    }
    // Work out a size correction factor.
    if (projected_size_based_on_q > FRAME_OVERHEAD_BITS)
        correction_factor = (int)((100 * (int64_t)ppcs_ptr->projected_frame_size) /
                                  projected_size_based_on_q);

    // More heavily damped adjustment used if we have been oscillating either side
    // of target.
    if (correction_factor > 0) {
        adjustment_limit = 0.25 + 0.5 * AOMMIN(1, fabs(log10(0.01 * correction_factor)));
    } else {
        adjustment_limit = 0.75;
    }

    rc->q_2_frame = rc->q_1_frame;
    rc->q_1_frame =
        ppcs_ptr->frm_hdr.quantization_params.base_q_idx; //cm->quant_params.base_qindex;
    rc->rc_2_frame = rc->rc_1_frame;
    if (correction_factor > 110)
        rc->rc_1_frame = -1;
    else if (correction_factor < 90)
        rc->rc_1_frame = 1;
    else
        rc->rc_1_frame = 0;

    if (correction_factor > 102) {
        // We are not already at the worst allowable quality
        correction_factor      = (int)(100 + ((correction_factor - 100) * adjustment_limit));
        rate_correction_factor = (rate_correction_factor * correction_factor) / 100;
        // Keep rate_correction_factor within limits
        if (rate_correction_factor > MAX_BPB_FACTOR)
            rate_correction_factor = MAX_BPB_FACTOR;
    } else if (correction_factor < 99) {
        // We are not already at the best allowable quality
        correction_factor      = (int)(100 - ((100 - correction_factor) * adjustment_limit));
        rate_correction_factor = (rate_correction_factor * correction_factor) / 100;

        // Keep rate_correction_factor within limits
        if (rate_correction_factor < MIN_BPB_FACTOR)
            rate_correction_factor = MIN_BPB_FACTOR;
    }

    set_rate_correction_factor(ppcs_ptr, rate_correction_factor /*, width, height*/);
}

// Update the buffer level: leaky bucket model.
static void update_buffer_level(PictureParentControlSet *ppcs_ptr, int encoded_frame_size) {
    SequenceControlSet *scs_ptr            = ppcs_ptr->scs_ptr;
    EncodeContext      *encode_context_ptr = scs_ptr->encode_context_ptr;
    RATE_CONTROL       *rc                 = &encode_context_ptr->rc;

    // Non-viewable frames are a special case and are treated as pure overhead.
    if (!ppcs_ptr->frm_hdr.showable_frame)
        rc->bits_off_target -= encoded_frame_size;
    else
        rc->bits_off_target += rc->avg_frame_bandwidth - encoded_frame_size;

    // Clip the buffer level to the maximum specified buffer size.
    rc->bits_off_target = AOMMIN(rc->bits_off_target, rc->maximum_buffer_size);
    rc->buffer_level    = rc->bits_off_target;

    //if (cpi->use_svc) update_layer_buffer_level(&cpi->svc, encoded_frame_size);
}

static void av1_rc_postencode_update(PictureParentControlSet *ppcs_ptr) {
    const RefreshFrameFlagsInfo *const refresh_frame_flags = &ppcs_ptr->refresh_frame;
    SequenceControlSet                *scs_ptr             = ppcs_ptr->scs_ptr;
    EncodeContext                     *encode_context_ptr  = scs_ptr->encode_context_ptr;
    RATE_CONTROL                      *rc                  = &encode_context_ptr->rc;
    GF_GROUP *const                    gf_group            = &encode_context_ptr->gf_group;
    CurrentFrame *const                current_frame       = &ppcs_ptr->av1_cm->current_frame;
    current_frame->frame_type                              = ppcs_ptr->frm_hdr.frame_type;
    FrameHeader *frm_hdr                                   = &ppcs_ptr->frm_hdr;
    const int    width                                     = ppcs_ptr->av1_cm->frm_size.frame_width;
    const int    height = ppcs_ptr->av1_cm->frm_size.frame_height;

    const int is_intrnl_arf = gf_group->update_type[ppcs_ptr->gf_group_index] == INTNL_ARF_UPDATE;

    const int qindex = frm_hdr->quantization_params.base_q_idx; //cm->quant_params.base_qindex;

    // Update rate control heuristics
    ppcs_ptr->projected_frame_size = (int)ppcs_ptr->total_num_bits;
    // Post encode loop adjustment of Q prediction.
    av1_rc_update_rate_correction_factors(ppcs_ptr, width, height);

    // Keep a record of last Q and ambient average Q.
    if (current_frame->frame_type == KEY_FRAME) {
        rc->avg_frame_qindex[KEY_FRAME] = ROUND_POWER_OF_TWO(
            3 * rc->avg_frame_qindex[KEY_FRAME] + qindex, 2);
        svt_block_on_mutex(encode_context_ptr->frame_updated_mutex);
        encode_context_ptr->frame_updated = 0;
        svt_release_mutex(encode_context_ptr->frame_updated_mutex);
    } else {
        svt_block_on_mutex(encode_context_ptr->frame_updated_mutex);
        encode_context_ptr->frame_updated++;
        svt_release_mutex(encode_context_ptr->frame_updated_mutex);
        if (/*(cpi->use_svc && cpi->oxcf.rc_cfg.mode == AOM_CBR) ||*/
            (!rc->is_src_frame_alt_ref &&
             !(refresh_frame_flags->golden_frame || is_intrnl_arf ||
               refresh_frame_flags->alt_ref_frame))) {
            rc->avg_frame_qindex[INTER_FRAME] = ROUND_POWER_OF_TWO(
                3 * rc->avg_frame_qindex[INTER_FRAME] + qindex, 2);
        }
    }

    // Keep record of last boosted (KF/GF/ARF) Q value.
    // If the current frame is coded at a lower Q then we also update it.
    // If all mbs in this group are skipped only update if the Q value is
    // better than that already stored.
    // This is used to help set quality in forced key frames to reduce popping
    if ((qindex < rc->last_boosted_qindex) || (current_frame->frame_type == KEY_FRAME) ||
        (!rc->constrained_gf_group &&
         (refresh_frame_flags->alt_ref_frame || is_intrnl_arf ||
          (refresh_frame_flags->golden_frame && !rc->is_src_frame_alt_ref)))) {
        rc->last_boosted_qindex = qindex;
    }
    update_buffer_level(ppcs_ptr, ppcs_ptr->projected_frame_size);
    rc->prev_avg_frame_bandwidth = rc->avg_frame_bandwidth;

    // Rolling monitors of whether we are over or underspending used to help
    // regulate min and Max Q in two pass.
    if (current_frame->frame_type != KEY_FRAME) {
        rc->rolling_target_bits = (int)ROUND_POWER_OF_TWO_64(
            rc->rolling_target_bits * 3 + ppcs_ptr->this_frame_target, 2);
        rc->rolling_actual_bits = (int)ROUND_POWER_OF_TWO_64(
            rc->rolling_actual_bits * 3 + ppcs_ptr->projected_frame_size, 2);
    }

    // Actual bits spent
    rc->total_actual_bits += ppcs_ptr->projected_frame_size;
    rc->total_target_bits += ppcs_ptr->frm_hdr.showable_frame ? rc->avg_frame_bandwidth : 0;

    rc->total_target_vs_actual = rc->total_actual_bits - rc->total_target_bits;

    if (current_frame->frame_type == KEY_FRAME)
        rc->frames_since_key = 0;
}
void update_rc_counts(PictureParentControlSet *ppcs_ptr) {
    SequenceControlSet *scs_ptr            = ppcs_ptr->scs_ptr;
    EncodeContext      *encode_context_ptr = scs_ptr->encode_context_ptr;
    RATE_CONTROL       *rc                 = &encode_context_ptr->rc;
    GF_GROUP *const     gf_group           = &encode_context_ptr->gf_group;

    if (ppcs_ptr->frm_hdr.showable_frame) {
        // If this is a show_existing_frame with a source other than altref,
        // or if it is not a displayed forward keyframe, the keyframe update
        // counters were incremented when it was originally encoded.
        rc->frames_since_key++;
        rc->frames_to_key--;
        rc->frames_till_gf_update_due--;
    }
    //update_gf_group_index(cpi);
    // Increment the gf group index ready for the next frame. If this is
    // a show_existing_frame with a source other than altref, or if it is not
    // a displayed forward keyframe, the index was incremented when it was
    // originally encoded.
    ++gf_group->index;
}
#define VBR_PCT_ADJUSTMENT_LIMIT 50
// For VBR...adjustment to the frame target based on error from previous frames
static void vbr_rate_correction(PictureControlSet *pcs_ptr, int *this_frame_target) {
    SequenceControlSet *scs_ptr             = pcs_ptr->parent_pcs_ptr->scs_ptr;
    EncodeContext      *encode_context_ptr  = scs_ptr->encode_context_ptr;
    RATE_CONTROL       *rc                  = &encode_context_ptr->rc;
    TWO_PASS *const     twopass             = &scs_ptr->twopass;
    int64_t             vbr_bits_off_target = rc->vbr_bits_off_target;
    const int           stats_count         = twopass->stats_buf_ctx->total_stats != NULL
                          ? (int)twopass->stats_buf_ctx->total_stats->count
                          : 0;
    const int frame_window = AOMMIN(16, (int)(stats_count - (int)pcs_ptr->picture_number));
    assert(VBR_PCT_ADJUSTMENT_LIMIT <= 100);
    if (frame_window > 0) {
        const int max_delta = (int)AOMMIN(
            abs((int)(vbr_bits_off_target / frame_window)),
            ((int64_t)(*this_frame_target) * VBR_PCT_ADJUSTMENT_LIMIT) / 100);

        // vbr_bits_off_target > 0 means we have extra bits to spend
        // vbr_bits_off_target < 0 we are currently overshooting
        *this_frame_target += (vbr_bits_off_target >= 0) ? max_delta : -max_delta;
    }

    // Fast redistribution of bits arising from massive local undershoot.
    // Dont do it for kf,arf,gf or overlay frames.
    if (!frame_is_kf_gf_arf(pcs_ptr->parent_pcs_ptr) && !rc->is_src_frame_alt_ref &&
        rc->vbr_bits_off_target_fast) {
        int one_frame_bits = AOMMAX(rc->avg_frame_bandwidth, *this_frame_target);
        int fast_extra_bits;
        fast_extra_bits = (int)AOMMIN(rc->vbr_bits_off_target_fast, one_frame_bits);
        fast_extra_bits = (int)AOMMIN(fast_extra_bits,
                                      AOMMAX(one_frame_bits / 8, rc->vbr_bits_off_target_fast / 8));
        *this_frame_target += (int)fast_extra_bits;
        rc->vbr_bits_off_target_fast -= fast_extra_bits;
    }
}

static INLINE void set_refresh_frame_flags(RefreshFrameFlagsInfo *const refresh_frame_flags,
                                           bool refresh_gf, bool refresh_bwdref, bool refresh_arf) {
    refresh_frame_flags->golden_frame  = refresh_gf;
    refresh_frame_flags->bwd_ref_frame = refresh_bwdref;
    refresh_frame_flags->alt_ref_frame = refresh_arf;
}

static void av1_configure_buffer_updates(PictureControlSet           *pcs_ptr,
                                         RefreshFrameFlagsInfo *const refresh_frame_flags,
                                         int                          force_refresh_all) {
    // NOTE(weitinglin): Should we define another function to take care of
    // cpi->rc.is_$Source_Type to make this function as it is in the comment?
    SequenceControlSet     *scs_ptr            = pcs_ptr->parent_pcs_ptr->scs_ptr;
    EncodeContext          *encode_context_ptr = scs_ptr->encode_context_ptr;
    RATE_CONTROL           *rc                 = &encode_context_ptr->rc;
    GF_GROUP               *gf_group           = &encode_context_ptr->gf_group;
    const FRAME_UPDATE_TYPE type = gf_group->update_type[pcs_ptr->parent_pcs_ptr->gf_group_index];
    rc->is_src_frame_alt_ref     = 0;

    switch (type) {
    case KF_UPDATE: set_refresh_frame_flags(refresh_frame_flags, true, true, true); break;

    case LF_UPDATE: set_refresh_frame_flags(refresh_frame_flags, false, false, false); break;

    case GF_UPDATE: set_refresh_frame_flags(refresh_frame_flags, true, false, false); break;

    case OVERLAY_UPDATE:
        set_refresh_frame_flags(refresh_frame_flags, true, false, false);
        rc->is_src_frame_alt_ref = 1;
        break;

    case ARF_UPDATE:
        // NOTE: BWDREF does not get updated along with ALTREF_FRAME.
        set_refresh_frame_flags(refresh_frame_flags, false, false, true);
        break;

    case INTNL_OVERLAY_UPDATE:
        set_refresh_frame_flags(refresh_frame_flags, false, false, false);
        rc->is_src_frame_alt_ref = 1;
        break;

    case INTNL_ARF_UPDATE: set_refresh_frame_flags(refresh_frame_flags, false, true, false); break;

    default: assert(0); break;
    }

    if (force_refresh_all)
        set_refresh_frame_flags(refresh_frame_flags, true, true, true);
}
static void av1_set_target_rate(PictureControlSet *pcs_ptr) {
    SequenceControlSet         *scs_ptr            = pcs_ptr->parent_pcs_ptr->scs_ptr;
    EncodeContext              *encode_context_ptr = scs_ptr->encode_context_ptr;
    int                         target_rate        = pcs_ptr->parent_pcs_ptr->base_frame_target;
    const RateControlCfg *const rc_cfg             = &encode_context_ptr->rc_cfg;
    // Correction to rate target based on prior over or under shoot.
    if (rc_cfg->mode == AOM_VBR || rc_cfg->mode == AOM_CQ)
        vbr_rate_correction(pcs_ptr, &target_rate);
    pcs_ptr->parent_pcs_ptr->this_frame_target = target_rate;
}
static double av1_get_compression_ratio(PictureParentControlSet *ppcs_ptr,
                                        size_t                   encoded_frame_size) {
    const int upscaled_width = ppcs_ptr->av1_cm->frm_size.superres_upscaled_width;
    const int height         = ppcs_ptr->av1_cm->frm_size.frame_height; //cm->height;
    const int luma_pic_size  = upscaled_width * height;
    const /*BITSTREAM_PROFILE*/ EbAv1SeqProfile profile = ppcs_ptr->scs_ptr->seq_header.seq_profile;
    const int pic_size_profile_factor                   = profile == /*PROFILE_0*/ MAIN_PROFILE
                          ? 15
                          : (profile == /*PROFILE_1*/ HIGH_PROFILE ? 30 : 36);
    encoded_frame_size = (encoded_frame_size > 129 ? encoded_frame_size - 128 : 1);
    const size_t uncompressed_frame_size = (luma_pic_size * pic_size_profile_factor) >> 3;
    return uncompressed_frame_size / (double)encoded_frame_size;
}
/**************************************************************************************************************
* get_kf_q_tpl()
* This function finds the q for a selected active quality for key frame. The functionality is the
* reverse of get_kf_active_quality_tpl()
**************************************************************************************************************/
static int get_kf_q_tpl(const RATE_CONTROL *const rc, int target_active_quality,
                        AomBitDepth bit_depth) {
    int *kf_low_motion_minq_cqp;
    int *kf_high_motion_minq;
    ASSIGN_MINQ_TABLE(bit_depth, kf_low_motion_minq_cqp);
    ASSIGN_MINQ_TABLE(bit_depth, kf_high_motion_minq);
    int q              = rc->active_worst_quality;
    int active_quality = get_active_quality(
        q, rc->kf_boost, kf_low, kf_high, kf_low_motion_minq_cqp, kf_high_motion_minq);
    int prev_dif = abs(target_active_quality - active_quality);
    while (abs(target_active_quality - active_quality) > 4 &&
           abs(target_active_quality - active_quality) <= prev_dif) {
        if (active_quality > target_active_quality)
            q--;
        else
            q++;
        active_quality = get_active_quality(
            q, rc->kf_boost, kf_low, kf_high, kf_low_motion_minq_cqp, kf_high_motion_minq);
    }
    return q;
}
/**************************************************************************************************************
*This function finds the q for a selected active quality for base layer frames. The functionality is the reverse of get_kf_active_quality_tpl()
**************************************************************************************************************/
static int get_gfu_q_tpl(const RATE_CONTROL *const rc, int target_active_quality,
                         AomBitDepth bit_depth) {
    int *arfgf_low_motion_minq;
    int *arfgf_high_motion_minq;
    ASSIGN_MINQ_TABLE(bit_depth, arfgf_low_motion_minq);
    ASSIGN_MINQ_TABLE(bit_depth, arfgf_high_motion_minq);

    int q              = rc->active_worst_quality;
    int active_quality = get_active_quality(q,
                                            rc->gfu_boost,
                                            gf_low_tpl_la,
                                            gf_high_tpl_la,
                                            arfgf_low_motion_minq,
                                            arfgf_high_motion_minq);

    int prev_dif = abs(target_active_quality - active_quality);
    while (abs(target_active_quality - active_quality) > 4 &&
           abs(target_active_quality - active_quality) <= prev_dif) {
        if (active_quality > target_active_quality)
            q--;
        else
            q++;
        active_quality = get_active_quality(q,
                                            rc->gfu_boost,
                                            gf_low_tpl_la,
                                            gf_high_tpl_la,
                                            arfgf_low_motion_minq,
                                            arfgf_high_motion_minq);
    }
    return q;
}
/**************************************************************************************************************
* capped_crf_reencode()
* This function performs re-encoding for capped CRF. It adjusts the QP, and active_worst_quality
**************************************************************************************************************/
static void capped_crf_reencode(PictureParentControlSet *ppcs_ptr, int *const q) {
    SequenceControlSet *scs_ptr            = ppcs_ptr->scs_ptr;
    EncodeContext      *encode_context_ptr = scs_ptr->encode_context_ptr;
    RATE_CONTROL *const rc                 = &encode_context_ptr->rc;

    if (ppcs_ptr->temporal_layer_index == 0 && ppcs_ptr->slice_type != 2) {
        ppcs_ptr->gf_group_index                                                    = 1;
        scs_ptr->encode_context_ptr->gf_group.update_type[ppcs_ptr->gf_group_index] = 6;
    }
    uint32_t frame_rate   = ((scs_ptr->frame_rate + (1 << (RC_PRECISION - 1))) >> RC_PRECISION);
    int      frames_in_sw = (int)rc->rate_average_periodin_frames;

    int64_t spent_bits_sw       = 0, available_bit_sw;
    int     coded_frames_num_sw = 0;
    // Find the start and the end of the sliding window
    int32_t start_index = ((ppcs_ptr->picture_number / frames_in_sw) * frames_in_sw) %
        CODED_FRAMES_STAT_QUEUE_MAX_DEPTH;
    int32_t end_index   = start_index + frames_in_sw;
    frames_in_sw        = (scs_ptr->passes > 1)
               ? MIN(end_index, (int32_t)scs_ptr->twopass.stats_buf_ctx->total_stats->count) - start_index
               : frames_in_sw;
    int64_t max_bits_sw = (int64_t)scs_ptr->static_config.max_bit_rate * (int32_t)frames_in_sw /
        frame_rate;
    max_bits_sw += (max_bits_sw * scs_ptr->static_config.mbr_over_shoot_pct / 100);
    // Loop over the sliding window and calculated the spent bits
    for (int index = start_index; index < end_index; index++) {
        int32_t queue_entry_index                 = (index > CODED_FRAMES_STAT_QUEUE_MAX_DEPTH - 1)
                            ? index - CODED_FRAMES_STAT_QUEUE_MAX_DEPTH
                            : index;
        coded_frames_stats_entry *queue_entry_ptr = rc->coded_frames_stat_queue[queue_entry_index];
        spent_bits_sw += (queue_entry_ptr->frame_total_bit_actual > 0)
            ? queue_entry_ptr->frame_total_bit_actual
            : 0;
        coded_frames_num_sw += (queue_entry_ptr->frame_total_bit_actual > 0) ? 1 : 0;
    }
    available_bit_sw = MAX(max_bits_sw - spent_bits_sw, 0);

    int remaining_frames       = frames_in_sw - coded_frames_num_sw;
    int available_bit_ratio    = (int)(100 * available_bit_sw / max_bits_sw);
    int available_frames_ratio = 100 * remaining_frames / frames_in_sw;

    int worst_quality = (int32_t)quantizer_to_qindex[scs_ptr->static_config.max_qp_allowed];
    if (*q < worst_quality && ppcs_ptr->projected_frame_size > ppcs_ptr->max_frame_size &&
        ppcs_ptr->temporal_layer_index == 0) {
        int          tmp_q;
        int          ref_qindex = rc->active_worst_quality;
        const double ref_q    = svt_av1_convert_qindex_to_q(ref_qindex, scs_ptr->encoder_bit_depth);
        int64_t      ref_bits = (int64_t)(ppcs_ptr->projected_frame_size);
        int64_t      target_bits = ppcs_ptr->max_frame_size;
        int          low         = rc->best_quality;
        int          high        = rc->worst_quality;

        while (low < high) {
            const int    mid      = (low + high) >> 1;
            const double q_tmp1   = svt_av1_convert_qindex_to_q(mid, scs_ptr->encoder_bit_depth);
            const int    mid_bits = (int)(ref_bits * ref_q / q_tmp1);

            if (mid_bits > target_bits)
                low = mid + 1;
            else
                high = mid;
        }
        tmp_q = low;

        rc->active_worst_quality = CLIP3(
            (int32_t)quantizer_to_qindex[scs_ptr->static_config.min_qp_allowed],
            (int32_t)quantizer_to_qindex[scs_ptr->static_config.max_qp_allowed],
            tmp_q);
#if DEBUG_RC_CAP_LOG
        if (ppcs_ptr->temporal_layer_index <= 0)
            printf("Reencode POC:%lld\tQindex:%d\t%d\t%d\tWorseActive%d\t%d\t%d\n",
                   ppcs_ptr->picture_number,
                   ppcs_ptr->frm_hdr.quantization_params.base_q_idx,
                   ppcs_ptr->projected_frame_size,
                   ppcs_ptr->max_frame_size,
                   rc->active_worst_quality,
                   ppcs_ptr->bottom_index,
                   ppcs_ptr->top_index);
#endif
        ppcs_ptr->top_index = rc->active_worst_quality;
        ppcs_ptr->q_high    = rc->active_worst_quality;
    }
    // Decrease the active worse quality based on the projected frame size and max frame size
    else if (ppcs_ptr->projected_frame_size < ppcs_ptr->max_frame_size &&
             ppcs_ptr->temporal_layer_index == 0 && ppcs_ptr->loop_count == 0 &&
             rc->active_worst_quality > quantizer_to_qindex[(uint8_t)scs_ptr->static_config.qp] &&
             (available_bit_ratio > available_frames_ratio)) {
        if (ppcs_ptr->projected_frame_size < ppcs_ptr->max_frame_size / 3)
            rc->active_worst_quality -= (rc->active_worst_quality / 5);
        else if (ppcs_ptr->projected_frame_size < ppcs_ptr->max_frame_size / 2)
            rc->active_worst_quality -= (rc->active_worst_quality / 8);
        else if (ppcs_ptr->projected_frame_size < 2 * ppcs_ptr->max_frame_size / 3)
            rc->active_worst_quality -= (rc->active_worst_quality / 12);

        rc->active_worst_quality = CLIP3(
            (int32_t)quantizer_to_qindex[(uint8_t)scs_ptr->static_config.qp],
            (int32_t)quantizer_to_qindex[scs_ptr->static_config.max_qp_allowed],
            rc->active_worst_quality);
    }
}
static void av1_rc_compute_frame_size_bounds(PictureParentControlSet *ppcs_ptr, int frame_target,
                                             int *frame_under_shoot_limit,
                                             int *frame_over_shoot_limit) {
    EncodeContext *const        encode_context_ptr = ppcs_ptr->scs_ptr->encode_context_ptr;
    RATE_CONTROL *const         rc                 = &(encode_context_ptr->rc);
    const RateControlCfg *const rc_cfg             = &encode_context_ptr->rc_cfg;
    if (rc_cfg->mode == AOM_Q) {
        const int tolerance = (int)AOMMAX(
            100, ((int64_t)encode_context_ptr->recode_tolerance * ppcs_ptr->max_frame_size) / 100);
        *frame_under_shoot_limit = ppcs_ptr->loop_count
            ? AOMMAX(ppcs_ptr->max_frame_size - tolerance, 0)
            : 0;
        *frame_over_shoot_limit  = AOMMIN(ppcs_ptr->max_frame_size + tolerance, INT_MAX);
    } else {
        // For very small rate targets where the fractional adjustment
        // may be tiny make sure there is at least a minimum range.
        assert(encode_context_ptr->recode_tolerance <= 100);
        const int tolerance = (int)AOMMAX(
            100, ((int64_t)encode_context_ptr->recode_tolerance * frame_target) / 100);
        *frame_under_shoot_limit = AOMMAX(frame_target - tolerance, 0);
        *frame_over_shoot_limit  = AOMMIN(frame_target + tolerance, rc->max_frame_bandwidth);
    }
}

// Function to test for conditions that indicate we should loop
// back and recode a frame.
static AOM_INLINE int recode_loop_test(PictureParentControlSet *ppcs_ptr, int high_limit,
                                       int low_limit, int q, int maxq, int minq) {
    EncodeContext *const        encode_context_ptr = ppcs_ptr->scs_ptr->encode_context_ptr;
    RATE_CONTROL *const         rc                 = &(encode_context_ptr->rc);
    const RateControlCfg *const rc_cfg             = &encode_context_ptr->rc_cfg;
    const int                   frame_is_kfgfarf   = frame_is_kf_gf_arf(ppcs_ptr);
    int                         force_recode       = 0;
    if ((ppcs_ptr->projected_frame_size >= rc->max_frame_bandwidth) ||
        (encode_context_ptr->recode_loop == ALLOW_RECODE) ||
        (frame_is_kfgfarf && (encode_context_ptr->recode_loop >= ALLOW_RECODE_KFMAXBW))) {
        // TODO(agrange) high_limit could be greater than the scale-down threshold.
        if ((ppcs_ptr->projected_frame_size > high_limit && q < maxq) ||
            (ppcs_ptr->projected_frame_size < low_limit && q > minq)) {
            force_recode = 1;
        } else if (rc_cfg->mode == AOM_CQ) {
            // Deal with frame undershoot and whether or not we are
            // below the automatically set cq level.
            if (q > rc_cfg->cq_level &&
                ppcs_ptr->projected_frame_size < ((ppcs_ptr->this_frame_target * 7) >> 3)) {
                force_recode = 1;
            }
        }
    }
    return force_recode;
}

// get overshoot regulated q based on q_low
static int get_regulated_q_overshoot(PictureParentControlSet *ppcs_ptr, int q_low, int q_high,
                                     int top_index, int bottom_index) {
    const int width  = ppcs_ptr->av1_cm->frm_size.frame_width;
    const int height = ppcs_ptr->av1_cm->frm_size.frame_height;

    av1_rc_update_rate_correction_factors(ppcs_ptr, width, height);

    int q_regulated = av1_rc_regulate_q(ppcs_ptr,
                                        ppcs_ptr->this_frame_target,
                                        bottom_index,
                                        AOMMAX(q_high, top_index),
                                        width,
                                        height);
    int retries     = 0;
    while (q_regulated < q_low && retries < 10) {
        av1_rc_update_rate_correction_factors(ppcs_ptr, width, height);
        q_regulated = av1_rc_regulate_q(ppcs_ptr,
                                        ppcs_ptr->this_frame_target,
                                        bottom_index,
                                        AOMMAX(q_high, top_index),
                                        width,
                                        height);
        retries++;
    }
    return q_regulated;
}

// get undershoot regulated q based on q_high
static AOM_INLINE int get_regulated_q_undershoot(PictureParentControlSet *ppcs_ptr, int q_high,
                                                 int top_index, int bottom_index) {
    const int width  = ppcs_ptr->av1_cm->frm_size.frame_width;
    const int height = ppcs_ptr->av1_cm->frm_size.frame_height;

    av1_rc_update_rate_correction_factors(ppcs_ptr, width, height);
    int q_regulated = av1_rc_regulate_q(
        ppcs_ptr, ppcs_ptr->this_frame_target, bottom_index, top_index, width, height);

    int retries = 0;
    while (q_regulated > q_high && retries < 10) {
        av1_rc_update_rate_correction_factors(ppcs_ptr, width, height);
        q_regulated = av1_rc_regulate_q(
            ppcs_ptr, ppcs_ptr->this_frame_target, bottom_index, top_index, width, height);
        retries++;
    }
    return q_regulated;
}

// This function works out whether we under- or over-shot
// our bitrate target and adjusts q as appropriate.  Also decides whether
// or not we should do another recode loop, indicated by *loop
void recode_loop_update_q(PictureParentControlSet *ppcs_ptr, int *const loop, int *const q,
                          int *const q_low, int *const q_high, const int top_index,
                          const int bottom_index, int *const undershoot_seen,
                          int *const overshoot_seen, int *const low_cr_seen, const int loop_count) {
    SequenceControlSet *const   scs_ptr            = ppcs_ptr->scs_ptr;
    EncodeContext *const        encode_context_ptr = scs_ptr->encode_context_ptr;
    RATE_CONTROL *const         rc                 = &(encode_context_ptr->rc);
    const RateControlCfg *const rc_cfg             = &encode_context_ptr->rc_cfg;
    const int do_dummy_pack = (scs_ptr->encode_context_ptr->recode_loop >= ALLOW_RECODE_KFMAXBW &&
                               !(rc_cfg->mode == AOM_Q &&
                                 scs_ptr->static_config.max_bit_rate == 0)) ||
        rc_cfg->min_cr > 0;
    if (do_dummy_pack) {
        svt_block_on_mutex(ppcs_ptr->pcs_total_rate_mutex);
        ppcs_ptr->projected_frame_size =
            (int)(((ppcs_ptr->pcs_total_rate + (1 << (AV1_PROB_COST_SHIFT - 1))) >>
                   AV1_PROB_COST_SHIFT) +
                  ((ppcs_ptr->frm_hdr.frame_type == KEY_FRAME) ? 13 : 0));
        svt_release_mutex(ppcs_ptr->pcs_total_rate_mutex);
    } else {
        ppcs_ptr->projected_frame_size = 0;
    }
    if (ppcs_ptr->loop_count && rc_cfg->mode != AOM_Q) {
        // scale rc->projected_frame_size with *0.8 for loop_count>=1
    }
    *loop = 0;
    if (scs_ptr->encode_context_ptr->recode_loop == ALLOW_RECODE_KFMAXBW &&
        ppcs_ptr->frm_hdr.frame_type != KEY_FRAME) {
        // skip re-encode for inter frame when setting -recode-loop 1
        return;
    }

    const int min_cr = rc_cfg->min_cr;
    if (min_cr > 0) {
        //aom_clear_system_state();
        const double compression_ratio = av1_get_compression_ratio(
            ppcs_ptr, ppcs_ptr->projected_frame_size >> 3);
        const double target_cr = min_cr / 100.0;
        if (compression_ratio < target_cr) {
            *low_cr_seen = 1;
            if (*q < rc->worst_quality) {
                const double cr_ratio    = target_cr / compression_ratio;
                const int    projected_q = AOMMAX(*q + 1, (int)(*q * cr_ratio * cr_ratio));
                *q                       = AOMMIN(AOMMIN(projected_q, *q + 32), rc->worst_quality);
                *q_low                   = AOMMAX(*q, *q_low);
                *q_high                  = AOMMAX(*q, *q_high);
                *loop                    = 1;
            }
        }
        if (*low_cr_seen)
            return;
    }
    // Used for capped CRF. Update the active worse quality
    if (rc_cfg->mode == AOM_Q && scs_ptr->static_config.max_bit_rate) {
        if (ppcs_ptr->temporal_layer_index > 0)
            return;
        else
            capped_crf_reencode(ppcs_ptr, q);
    }
    const int last_q                 = *q;
    int       frame_over_shoot_limit = 0, frame_under_shoot_limit = 0;
    av1_rc_compute_frame_size_bounds(
        ppcs_ptr, ppcs_ptr->this_frame_target, &frame_under_shoot_limit, &frame_over_shoot_limit);
    if (frame_over_shoot_limit == 0)
        frame_over_shoot_limit = 1;

    if (recode_loop_test(ppcs_ptr,
                         frame_over_shoot_limit,
                         frame_under_shoot_limit,
                         *q,
                         AOMMAX(*q_high, top_index),
                         bottom_index)) {
        const int width  = ppcs_ptr->av1_cm->frm_size.frame_width;
        const int height = ppcs_ptr->av1_cm->frm_size.frame_height;
        // Is the projected frame size out of range and are we allowed
        // to attempt to recode.

        // Frame size out of permitted range:
        // Update correction factor & compute new Q to try...
        // Frame is too large
        if (ppcs_ptr->projected_frame_size > ppcs_ptr->this_frame_target) {
            // Special case if the projected size is > the max allowed.
            if (*q == *q_high && ppcs_ptr->projected_frame_size >= rc->max_frame_bandwidth) {
                const double q_val_high_current = svt_av1_convert_qindex_to_q(
                    *q_high, scs_ptr->static_config.encoder_bit_depth);
                const double q_val_high_new = q_val_high_current *
                    ((double)ppcs_ptr->projected_frame_size / rc->max_frame_bandwidth);
                *q_high = av1_find_qindex(q_val_high_new,
                                          scs_ptr->static_config.encoder_bit_depth,
                                          rc->best_quality,
                                          rc->worst_quality);
            }
            // Raise Qlow as to at least the current value
            *q_low = AOMMIN(*q + 1, *q_high);

            if (*undershoot_seen || loop_count > 2 ||
                (loop_count == 2 && !frame_is_intra_only(ppcs_ptr))) {
                av1_rc_update_rate_correction_factors(ppcs_ptr, width, height);

                *q = (*q_high + *q_low + 1) / 2;
            } else if (loop_count == 2 && frame_is_intra_only(ppcs_ptr)) {
                const int q_mid       = (*q_high + *q_low + 1) / 2;
                const int q_regulated = get_regulated_q_overshoot(
                    ppcs_ptr, *q_low, *q_high, top_index, bottom_index);
                // Get 'q' in-between 'q_mid' and 'q_regulated' for a smooth
                // transition between loop_count < 2 and loop_count > 2.
                *q = (q_mid + q_regulated + 1) / 2;
            } else {
                *q = get_regulated_q_overshoot(ppcs_ptr, *q_low, *q_high, top_index, bottom_index);
            }

            *overshoot_seen = 1;
        } else {
            // Frame is too small
            *q_high = AOMMAX(*q - 1, *q_low);

            if (*overshoot_seen || loop_count > 2 ||
                (loop_count == 2 && !frame_is_intra_only(ppcs_ptr))) {
                av1_rc_update_rate_correction_factors(ppcs_ptr, width, height);
                *q = (*q_high + *q_low) / 2;
            } else if (loop_count == 2 && frame_is_intra_only(ppcs_ptr)) {
                const int q_mid       = (*q_high + *q_low) / 2;
                const int q_regulated = get_regulated_q_undershoot(
                    ppcs_ptr, *q_high, top_index, bottom_index);
                // Get 'q' in-between 'q_mid' and 'q_regulated' for a smooth
                // transition between loop_count < 2 and loop_count > 2.
                *q = (q_mid + q_regulated) / 2;

                // Special case reset for qlow for constrained quality.
                // This should only trigger where there is very substantial
                // undershoot on a frame and the auto cq level is above
                // the user passsed in value.
                if (rc_cfg->mode == AOM_CQ && q_regulated < *q_low) {
                    *q_low = *q;
                }
            } else {
                *q = get_regulated_q_undershoot(ppcs_ptr, *q_high, top_index, bottom_index);

                // Special case reset for qlow for constrained quality.
                // This should only trigger where there is very substantial
                // undershoot on a frame and the auto cq level is above
                // the user passsed in value.
                if (rc_cfg->mode == AOM_CQ && *q < *q_low) {
                    *q_low = *q;
                }
            }

            *undershoot_seen = 1;
        }

        // Clamp Q to upper and lower limits:
        *q = clamp(*q, *q_low, *q_high);
    }

    *q    = (uint8_t)CLIP3((int32_t)quantizer_to_qindex[scs_ptr->static_config.min_qp_allowed],
                        (int32_t)quantizer_to_qindex[scs_ptr->static_config.max_qp_allowed],
                        *q);
    *loop = (*q != last_q);
    // Used for capped CRF. Update the active worse quality based on the final assigned qindex
    if (rc_cfg->mode == AOM_Q && scs_ptr->static_config.max_bit_rate && *loop == 0 &&
        ppcs_ptr->temporal_layer_index == 0 && ppcs_ptr->loop_count > 0) {
        if (ppcs_ptr->slice_type == I_SLICE)
            rc->active_worst_quality = get_kf_q_tpl(
                rc, *q, scs_ptr->static_config.encoder_bit_depth);
        else
            rc->active_worst_quality = get_gfu_q_tpl(
                rc, *q, scs_ptr->static_config.encoder_bit_depth);

        rc->active_worst_quality = CLIP3(
            (int32_t)quantizer_to_qindex[scs_ptr->static_config.min_qp_allowed],
            (int32_t)quantizer_to_qindex[scs_ptr->static_config.max_qp_allowed],
            rc->active_worst_quality);
    }
}
/************************************************************************************************
* Populate the required parameters in two_pass structure from other structures
*************************************************************************************************/
static void restore_two_pass_param(PictureParentControlSet         *ppcs_ptr,
                                   RateControlIntervalParamContext *rate_control_param_ptr) {
    SequenceControlSet *scs_ptr = ppcs_ptr->scs_ptr;
    TWO_PASS *const     twopass = &scs_ptr->twopass;
    if (ppcs_ptr->scs_ptr->enable_dec_order == 1 && ppcs_ptr->scs_ptr->lap_enabled &&
        ppcs_ptr->temporal_layer_index == 0) {
        for (uint64_t num_frames = ppcs_ptr->stats_in_offset;
             num_frames < ppcs_ptr->stats_in_end_offset;
             ++num_frames) {
            FIRSTPASS_STATS *cur_frame = ppcs_ptr->scs_ptr->twopass.stats_buf_ctx->stats_in_start +
                num_frames;
            if ((int64_t)cur_frame->frame >
                ppcs_ptr->scs_ptr->twopass.stats_buf_ctx->last_frame_accumulated) {
                svt_av1_accumulate_stats(ppcs_ptr->scs_ptr->twopass.stats_buf_ctx->total_stats,
                                         cur_frame);
                ppcs_ptr->scs_ptr->twopass.stats_buf_ctx->last_frame_accumulated =
                    (int64_t)cur_frame->frame;
            }
        }
    }
    twopass->stats_in = scs_ptr->twopass.stats_buf_ctx->stats_in_start + ppcs_ptr->stats_in_offset;
    twopass->stats_buf_ctx->stats_in_end = scs_ptr->twopass.stats_buf_ctx->stats_in_start +
        ppcs_ptr->stats_in_end_offset;
    twopass->kf_group_bits       = rate_control_param_ptr->kf_group_bits;
    twopass->kf_group_error_left = rate_control_param_ptr->kf_group_error_left;
}
/************************************************************************************************
* Populate the required parameters in gf_group structure from other structures
*************************************************************************************************/
static void restore_gf_group_param(PictureParentControlSet *ppcs_ptr) {
    SequenceControlSet *scs_ptr            = ppcs_ptr->scs_ptr;
    EncodeContext      *encode_context_ptr = scs_ptr->encode_context_ptr;
    GF_GROUP *const     gf_group           = &encode_context_ptr->gf_group;
    gf_group->index                        = ppcs_ptr->gf_group_index;
    gf_group->size                         = ppcs_ptr->gf_group_size;
    gf_group->update_type[gf_group->index] = ppcs_ptr->update_type;
    gf_group->layer_depth[gf_group->index] = ppcs_ptr->layer_depth;
    gf_group->arf_boost[gf_group->index]   = ppcs_ptr->arf_boost;
}
/************************************************************************************************
* Populate the required parameters in rc, twopass and gf_group structures from other structures
*************************************************************************************************/
static void restore_param(PictureParentControlSet         *ppcs_ptr,
                          RateControlIntervalParamContext *rate_control_param_ptr) {
    SequenceControlSet *scs_ptr            = ppcs_ptr->scs_ptr;
    EncodeContext      *encode_context_ptr = scs_ptr->encode_context_ptr;
    if (!(!(scs_ptr->static_config.pass == ENC_MIDDLE_PASS ||
            scs_ptr->static_config.pass == ENC_LAST_PASS) &&
          scs_ptr->static_config.pass != ENC_FIRST_PASS &&
          scs_ptr->static_config.rate_control_mode == 2))
        restore_two_pass_param(ppcs_ptr, rate_control_param_ptr);

    const KeyFrameCfg *const kf_cfg = &encode_context_ptr->kf_cfg;
    ppcs_ptr->frames_since_key      = (int)(ppcs_ptr->decode_order - ppcs_ptr->last_idr_picture);
    int key_max                     = kf_cfg->key_freq_max;
    if (scs_ptr->lap_enabled) {
        if (scs_ptr->static_config.hierarchical_levels != ppcs_ptr->hierarchical_levels ||
            ppcs_ptr->end_of_sequence_region)
            key_max = (int)MIN(
                kf_cfg->key_freq_max,
                (int)((int64_t)((scs_ptr->twopass.stats_buf_ctx->stats_in_end - 1)->frame) -
                      ppcs_ptr->last_idr_picture + 1));
        else
            key_max = kf_cfg->key_freq_max;
    } else {
        if (!(!(scs_ptr->static_config.pass == ENC_MIDDLE_PASS ||
                scs_ptr->static_config.pass == ENC_LAST_PASS) &&
              scs_ptr->static_config.pass != ENC_FIRST_PASS &&
              scs_ptr->static_config.rate_control_mode == 2))
            key_max = (int)MIN(
                kf_cfg->key_freq_max,
                (int)((int64_t)((scs_ptr->twopass.stats_buf_ctx->stats_in_end - 1)->frame) -
                      ppcs_ptr->last_idr_picture + 1));
    }
    if (!(!(scs_ptr->static_config.pass == ENC_MIDDLE_PASS ||
            scs_ptr->static_config.pass == ENC_LAST_PASS) &&
          scs_ptr->static_config.pass != ENC_FIRST_PASS &&
          scs_ptr->static_config.rate_control_mode == 2)) {
        ppcs_ptr->frames_to_key = key_max - ppcs_ptr->frames_since_key;
        TWO_PASS *const twopass = &scs_ptr->twopass;
        RATE_CONTROL   *rc      = &encode_context_ptr->rc;
        // For the last minigop of the sequence, when look ahead is not long enough to find the GOP size, the GOP size is set
        // to kf_cfg->key_freq_max and the kf_group_bits is calculated based on that. However, when we get closer to the end, the
        // end of sequence will be in the look ahead and frames_to_key is updated. In this case, kf_group_bits is calculated based
        // on the new GOP size
        if (scs_ptr->lap_enabled && (kf_cfg->key_freq_max != ppcs_ptr->frames_since_key) &&
            (scs_ptr->lad_mg + 1) * (1 << scs_ptr->static_config.hierarchical_levels) <
                scs_ptr->static_config.intra_period_length &&
            (scs_ptr->static_config.hierarchical_levels != ppcs_ptr->hierarchical_levels ||
             ppcs_ptr->end_of_sequence_region) &&
            !rate_control_param_ptr->end_of_seq_seen) {
            twopass->kf_group_bits = (ppcs_ptr->frames_to_key) * twopass->kf_group_bits /
                (kf_cfg->key_freq_max - ppcs_ptr->frames_since_key);
            rate_control_param_ptr->end_of_seq_seen = 1;
        }
        rc->frames_to_key    = ppcs_ptr->frames_to_key;
        rc->frames_since_key = ppcs_ptr->frames_since_key;
    }
    restore_gf_group_param(ppcs_ptr);
}
/************************************************************************************************
* Store the required parameters from rc structure to other structures
*************************************************************************************************/
static void store_rc_param(PictureParentControlSet *ppcs_ptr) {
    SequenceControlSet *scs_ptr            = ppcs_ptr->scs_ptr;
    EncodeContext      *encode_context_ptr = scs_ptr->encode_context_ptr;
    RATE_CONTROL       *rc                 = &encode_context_ptr->rc;

    ppcs_ptr->is_src_frame_alt_ref = ppcs_ptr->is_overlay;
    if (ppcs_ptr->is_new_gf_group) {
        for (uint8_t frame_idx = 0; frame_idx < (int32_t)ppcs_ptr->gf_interval; frame_idx++) {
            ppcs_ptr->gf_group[frame_idx]->num_stats_used_for_gfu_boost =
                rc->num_stats_used_for_gfu_boost;
            ppcs_ptr->gf_group[frame_idx]->num_stats_required_for_gfu_boost =
                rc->num_stats_required_for_gfu_boost;
        }
    }
}
/************************************************************************************************
* Store the required parameters in two_pass structure from other structures
*************************************************************************************************/
static void store_two_pass_param(PictureParentControlSet         *ppcs_ptr,
                                 RateControlIntervalParamContext *rate_control_param_ptr) {
    SequenceControlSet *scs_ptr = ppcs_ptr->scs_ptr;
    TWO_PASS *const     twopass = &scs_ptr->twopass;

    rate_control_param_ptr->kf_group_bits       = twopass->kf_group_bits;
    rate_control_param_ptr->kf_group_error_left = twopass->kf_group_error_left;
}
/************************************************************************************************
* Store the required parameters from gf_group structure to other structures
*************************************************************************************************/
static void store_gf_group_param(PictureParentControlSet *ppcs_ptr) {
    SequenceControlSet *scs_ptr            = ppcs_ptr->scs_ptr;
    EncodeContext      *encode_context_ptr = scs_ptr->encode_context_ptr;
    GF_GROUP *const     gf_group           = &encode_context_ptr->gf_group;
    if (ppcs_ptr->is_new_gf_group) {
        for (uint8_t frame_idx = 0; frame_idx < (int32_t)ppcs_ptr->gf_interval; frame_idx++) {
            uint8_t gf_group_index = ppcs_ptr->slice_type == I_SLICE ? frame_idx : frame_idx + 1;
            ppcs_ptr->gf_group[frame_idx]->gf_group_index = gf_group_index;
            ppcs_ptr->gf_group[frame_idx]->gf_group_size  = MAX(gf_group->size,
                                                               ppcs_ptr->gf_interval);
            ppcs_ptr->gf_group[frame_idx]->update_type    = gf_group->update_type[gf_group_index];
            ppcs_ptr->gf_group[frame_idx]->layer_depth    = gf_group->layer_depth[gf_group_index];
            ppcs_ptr->gf_group[frame_idx]->arf_boost      = gf_group->arf_boost[gf_group_index];
            ppcs_ptr->gf_group[frame_idx]->base_frame_target =
                gf_group->bit_allocation[gf_group_index];
        }
    }
}
/************************************************************************************************
* Store the required parameters from rc, twopass and gf_group structures to other structures
*************************************************************************************************/
static void store_param(PictureParentControlSet         *ppcs_ptr,
                        RateControlIntervalParamContext *rate_control_param_ptr) {
    store_rc_param(ppcs_ptr);
    store_two_pass_param(ppcs_ptr, rate_control_param_ptr);
    store_gf_group_param(ppcs_ptr);
}
/************************************************************************************************
* Calculates the stat of coded frames over the averaging period
*************************************************************************************************/
static void coded_frames_stat_calc(PictureParentControlSet *ppcs_ptr) {
    int32_t                   queue_entry_index;
    coded_frames_stats_entry *queue_entry_ptr;
    EbBool                    move_slide_window_flag = EB_TRUE;
    EbBool                    end_of_sequence_flag   = EB_TRUE;

    SequenceControlSet *scs_ptr = (SequenceControlSet *)ppcs_ptr->scs_wrapper_ptr->object_ptr;
    EncodeContext      *encode_context_ptr = scs_ptr->encode_context_ptr;
    RATE_CONTROL       *rc                 = &encode_context_ptr->rc;
    // Determine offset from the Head Ptr
    queue_entry_index =
        (int32_t)(ppcs_ptr->picture_number -
                  rc->coded_frames_stat_queue[rc->coded_frames_stat_queue_head_index]
                      ->picture_number);
    queue_entry_index += rc->coded_frames_stat_queue_head_index;
    queue_entry_index = (queue_entry_index > CODED_FRAMES_STAT_QUEUE_MAX_DEPTH - 1)
        ? queue_entry_index - CODED_FRAMES_STAT_QUEUE_MAX_DEPTH
        : queue_entry_index;
    queue_entry_ptr   = rc->coded_frames_stat_queue[queue_entry_index];

    queue_entry_ptr->frame_total_bit_actual = (uint64_t)ppcs_ptr->total_num_bits;
    queue_entry_ptr->picture_number         = ppcs_ptr->picture_number;
    queue_entry_ptr->end_of_sequence_flag   = ppcs_ptr->end_of_sequence_flag;

    move_slide_window_flag = EB_TRUE;
    while (move_slide_window_flag) {
        // Check if the sliding window condition is valid
        uint32_t queue_entry_index_temp = rc->coded_frames_stat_queue_head_index;
        if (rc->coded_frames_stat_queue[queue_entry_index_temp]->frame_total_bit_actual != -1)
            end_of_sequence_flag =
                rc->coded_frames_stat_queue[queue_entry_index_temp]->end_of_sequence_flag;
        else
            end_of_sequence_flag = EB_FALSE;
        while (move_slide_window_flag && !end_of_sequence_flag &&
               queue_entry_index_temp <
                   rc->coded_frames_stat_queue_head_index + rc->rate_average_periodin_frames) {
            uint32_t queue_entry_index_temp2 = (queue_entry_index_temp >
                                                CODED_FRAMES_STAT_QUEUE_MAX_DEPTH - 1)
                ? queue_entry_index_temp - CODED_FRAMES_STAT_QUEUE_MAX_DEPTH
                : queue_entry_index_temp;

            move_slide_window_flag = (EbBool)(move_slide_window_flag &&
                                              (rc->coded_frames_stat_queue[queue_entry_index_temp2]
                                                   ->frame_total_bit_actual != -1));

            if (rc->coded_frames_stat_queue[queue_entry_index_temp2]->frame_total_bit_actual !=
                -1) {
                // check if it is the last frame. If we have reached the last frame, we would output the buffered frames in the Queue.
                end_of_sequence_flag =
                    rc->coded_frames_stat_queue[queue_entry_index_temp2]->end_of_sequence_flag;
            } else
                end_of_sequence_flag = EB_FALSE;
            queue_entry_index_temp++;
        }

        if (move_slide_window_flag) {
            //get a new entry spot
            queue_entry_ptr = (rc->coded_frames_stat_queue[rc->coded_frames_stat_queue_head_index]);
            queue_entry_index_temp = rc->coded_frames_stat_queue_head_index;
            // This is set to false, so the last frame would go inside the loop
            end_of_sequence_flag        = EB_FALSE;
            uint32_t frames_in_sw       = 0;
            rc->total_bit_actual_per_sw = 0;

            while (!end_of_sequence_flag &&
                   queue_entry_index_temp <
                       rc->coded_frames_stat_queue_head_index + rc->rate_average_periodin_frames) {
                frames_in_sw++;

                uint32_t queue_entry_index_temp2 = (queue_entry_index_temp >
                                                    CODED_FRAMES_STAT_QUEUE_MAX_DEPTH - 1)
                    ? queue_entry_index_temp - CODED_FRAMES_STAT_QUEUE_MAX_DEPTH
                    : queue_entry_index_temp;

                rc->total_bit_actual_per_sw +=
                    rc->coded_frames_stat_queue[queue_entry_index_temp2]->frame_total_bit_actual;
                end_of_sequence_flag =
                    rc->coded_frames_stat_queue[queue_entry_index_temp2]->end_of_sequence_flag;

                queue_entry_index_temp++;
            }
            uint32_t frame_rate = ((scs_ptr->frame_rate + (1 << (RC_PRECISION - 1))) >>
                                   RC_PRECISION);
            assert(frames_in_sw > 0);
            if (frames_in_sw == (uint32_t)rc->rate_average_periodin_frames) {
                rc->max_bit_actual_per_sw = MAX(
                    rc->max_bit_actual_per_sw,
                    rc->total_bit_actual_per_sw * frame_rate / frames_in_sw / 1000);
                if (queue_entry_ptr->picture_number % ((rc->rate_average_periodin_frames)) == 0) {
                    rc->max_bit_actual_per_gop = MAX(
                        rc->max_bit_actual_per_gop,
                        rc->total_bit_actual_per_sw * frame_rate / frames_in_sw / 1000);
                    rc->min_bit_actual_per_gop = MIN(
                        rc->min_bit_actual_per_gop,
                        rc->total_bit_actual_per_sw * frame_rate / frames_in_sw / 1000);
#if DEBUG_RC_CAP_LOG
                    SVT_LOG("POC:%d\t%.0f\t%.2f%% \n",
                            (int)queue_entry_ptr->picture_number,
                            (double)((int64_t)rc->total_bit_actual_per_sw * frame_rate /
                                     frames_in_sw / 1000),
                            (double)(100 * (double)rc->total_bit_actual_per_sw * frame_rate /
                                     frames_in_sw /
                                     MAX((double)scs_ptr->static_config.max_bit_rate, 1.0)) -
                                100);
#endif
                }
            }
#if DEBUG_RC_CAP_LOG
            if (frames_in_sw == rc->rate_average_periodin_frames - 1) {
                SVT_LOG("\n%d GopMax\t", (int32_t)rc->max_bit_actual_per_gop);
                SVT_LOG("%d GopMin\n", (int32_t)rc->min_bit_actual_per_gop);
            }
#endif
            // Reset the Queue Entry
            queue_entry_ptr->picture_number += CODED_FRAMES_STAT_QUEUE_MAX_DEPTH;
            queue_entry_ptr->frame_total_bit_actual = -1;

            // Increment the Queue head Ptr
            rc->coded_frames_stat_queue_head_index = (rc->coded_frames_stat_queue_head_index ==
                                                      CODED_FRAMES_STAT_QUEUE_MAX_DEPTH - 1)
                ? 0
                : rc->coded_frames_stat_queue_head_index + 1;

            queue_entry_ptr = (rc->coded_frames_stat_queue[rc->coded_frames_stat_queue_head_index]);
        }
    }
}
/****************************************************************************************
* reset_rc_param
* reset RC related variable in PPCS
*****************************************************************************************/
void reset_rc_param(PictureParentControlSet *ppcs) {
    ppcs->loop_count      = 0;
    ppcs->overshoot_seen  = 0;
    ppcs->undershoot_seen = 0;
}

void *rate_control_kernel(void *input_ptr) {
    // Context
    EbThreadContext    *thread_context_ptr = (EbThreadContext *)input_ptr;
    RateControlContext *context_ptr        = (RateControlContext *)thread_context_ptr->priv;

    RateControlIntervalParamContext *rate_control_param_ptr;

    PictureControlSet       *pcs_ptr;
    PictureParentControlSet *parentpicture_control_set_ptr;

    // Config
    SequenceControlSet *scs_ptr;

    // Input
    EbObjectWrapper  *rate_control_tasks_wrapper_ptr;
    RateControlTasks *rate_control_tasks_ptr;

    // Output
    EbObjectWrapper    *rate_control_results_wrapper_ptr;
    RateControlResults *rate_control_results_ptr;

    uint64_t total_number_of_fb_frames = 0;

    RateControlTaskTypes task_type;
    RATE_CONTROL        *rc;
    for (;;) {
        // Get RateControl Task
        EB_GET_FULL_OBJECT(context_ptr->rate_control_input_tasks_fifo_ptr,
                           &rate_control_tasks_wrapper_ptr);

        rate_control_tasks_ptr = (RateControlTasks *)rate_control_tasks_wrapper_ptr->object_ptr;
        task_type              = rate_control_tasks_ptr->task_type;
        EbBool is_superres_recode_task = (task_type == RC_INPUT_SUPERRES_RECODE) ? EB_TRUE
                                                                                 : EB_FALSE;

        // Modify these for different temporal layers later
        switch (task_type) {
        case RC_INPUT_SUPERRES_RECODE:
            assert(scs_ptr->static_config.superres_mode == SUPERRES_QTHRESH ||
                   scs_ptr->static_config.superres_mode == SUPERRES_AUTO);
            // intentionally reuse code in RC_INPUT
        case RC_INPUT:
            pcs_ptr = (PictureControlSet *)rate_control_tasks_ptr->pcs_wrapper_ptr->object_ptr;

            scs_ptr              = (SequenceControlSet *)pcs_ptr->scs_wrapper_ptr->object_ptr;
            FrameHeader *frm_hdr = &pcs_ptr->parent_pcs_ptr->frm_hdr;
            rc                   = &scs_ptr->encode_context_ptr->rc;
            if (scs_ptr->passes > 1 && scs_ptr->static_config.max_bit_rate)
                rc->rate_average_periodin_frames =
                    (uint64_t)scs_ptr->twopass.stats_buf_ctx->total_stats->count;
            else
                rc->rate_average_periodin_frames = 60;
            if (!is_superres_recode_task) {
                pcs_ptr->parent_pcs_ptr->blk_lambda_tuning = EB_FALSE;
            }
            reset_rc_param(pcs_ptr->parent_pcs_ptr);
            // Frame level RC. Find the ParamPtr for the current GOP
            {
                uint32_t interval_index_temp = 0;
                EbBool   interval_found      = EB_FALSE;
                while ((interval_index_temp < PARALLEL_GOP_MAX_NUMBER) && !interval_found) {
                    if (pcs_ptr->picture_number >=
                            context_ptr->rate_control_param_queue[interval_index_temp]->first_poc &&
                        pcs_ptr->picture_number <=
                            context_ptr->rate_control_param_queue[interval_index_temp]->last_poc) {
                        interval_found = EB_TRUE;
                    } else
                        interval_index_temp++;
                }
                CHECK_REPORT_ERROR(interval_index_temp != PARALLEL_GOP_MAX_NUMBER,
                                   scs_ptr->encode_context_ptr->app_callback_ptr,
                                   EB_ENC_RC_ERROR2);

                rate_control_param_ptr = context_ptr->rate_control_param_queue[interval_index_temp];
            }
            pcs_ptr->parent_pcs_ptr->rate_control_param_ptr = rate_control_param_ptr;
            if (!is_superres_recode_task) {
                if (scs_ptr->static_config.rate_control_mode) {
                    if (pcs_ptr->picture_number == 0) {
                        set_rc_buffer_sizes(scs_ptr);
                        av1_rc_init(scs_ptr);
                    }
                    restore_param(pcs_ptr->parent_pcs_ptr, rate_control_param_ptr);
                    if (scs_ptr->static_config.rate_control_mode == 2 &&
                        scs_ptr->static_config.pred_structure == EB_PRED_LOW_DELAY_P)
                        svt_av1_get_one_pass_rt_params(pcs_ptr->parent_pcs_ptr);
                    else
                        svt_av1_get_second_pass_params(pcs_ptr->parent_pcs_ptr);

                    av1_configure_buffer_updates(
                        pcs_ptr, &(pcs_ptr->parent_pcs_ptr->refresh_frame), 0);
                    av1_set_target_rate(pcs_ptr);
                    store_param(pcs_ptr->parent_pcs_ptr, rate_control_param_ptr);
                }
            }

            if (scs_ptr->static_config.rate_control_mode == 0) {
                // if RC mode is 0,  fixed QP is used
                // QP scaling based on POC number for Flat IPPP structure
                frm_hdr->quantization_params.base_q_idx = quantizer_to_qindex[pcs_ptr->picture_qp];

                if (scs_ptr->static_config.use_fixed_qindex_offsets == 1) {
                    pcs_ptr->picture_qp = scs_ptr->static_config.qp;
                    int32_t qindex      = quantizer_to_qindex[(uint8_t)scs_ptr->static_config.qp];
                    if (!frame_is_intra_only(pcs_ptr->parent_pcs_ptr)) {
                        qindex +=
                            scs_ptr->static_config.qindex_offsets[pcs_ptr->temporal_layer_index];
                    } else {
                        qindex += scs_ptr->static_config.key_frame_qindex_offset;
                    }
                    qindex = CLIP3(quantizer_to_qindex[scs_ptr->static_config.min_qp_allowed],
                                   quantizer_to_qindex[scs_ptr->static_config.max_qp_allowed],
                                   qindex);
                    int32_t chroma_qindex = qindex;
                    if (frame_is_intra_only(pcs_ptr->parent_pcs_ptr)) {
                        chroma_qindex += scs_ptr->static_config.key_frame_chroma_qindex_offset;
                    } else {
                        chroma_qindex += scs_ptr->static_config
                                             .chroma_qindex_offsets[pcs_ptr->temporal_layer_index];
                    }

                    chroma_qindex = CLIP3(
                        quantizer_to_qindex[scs_ptr->static_config.min_qp_allowed],
                        quantizer_to_qindex[scs_ptr->static_config.max_qp_allowed],
                        chroma_qindex);
                    frm_hdr->quantization_params.base_q_idx = qindex;
                    frm_hdr->quantization_params.delta_q_dc[1] =
                        frm_hdr->quantization_params.delta_q_dc[2] = frm_hdr->quantization_params
                                                                         .delta_q_ac[1] =
                            frm_hdr->quantization_params.delta_q_ac[2] = (chroma_qindex - qindex);
                    pcs_ptr->picture_qp                                = (uint8_t)CLIP3(
                        (int32_t)scs_ptr->static_config.min_qp_allowed,
                        (int32_t)scs_ptr->static_config.max_qp_allowed,
                        (frm_hdr->quantization_params.base_q_idx + 2) >> 2);
                    /*
                    printf("\nSVT: Frame Type = %s, PicNumber = %lld, DecoderOrder = %lld, Temp Layer "
                           "Index = %d, Picture QP = %d, Picture Qindex = %d, Picture Chroma Qindex = %d\n",
                           frame_is_intra_only(pcs_ptr->parent_pcs_ptr) ? "INTRA" : "INTER",
                           pcs_ptr->picture_number,
                           pcs_ptr->parent_pcs_ptr->decode_order,
                           pcs_ptr->temporal_layer_index,
                           pcs_ptr->picture_qp,
                           frm_hdr->quantization_params.base_q_idx,
                           chroma_qindex);
*/
                } else if (scs_ptr->enable_qp_scaling_flag &&
                           pcs_ptr->parent_pcs_ptr->qp_on_the_fly == EB_FALSE) {
                    const int32_t qindex = quantizer_to_qindex[(uint8_t)scs_ptr->static_config.qp];
                    // if there are need enough pictures in the LAD/SlidingWindow, the adaptive QP scaling is not used
                    int32_t new_qindex;
                    if (scs_ptr->static_config.pass != ENC_FIRST_PASS) {
                        // Content adaptive qp assignment
                        if (pcs_ptr->parent_pcs_ptr->tpl_ctrls.enable) {
                            if (pcs_ptr->picture_number == 0) {
                                rc->active_worst_quality =
                                    quantizer_to_qindex[(uint8_t)scs_ptr->static_config.qp];
                                av1_rc_init(scs_ptr);
                            }
                            new_qindex = cqp_qindex_calc_tpl_la(
                                pcs_ptr, rc, rc->active_worst_quality);
                        } else
                            new_qindex = cqp_qindex_calc(pcs_ptr, qindex);
                    } else {
                        new_qindex = find_fp_qindex(
                            (AomBitDepth)scs_ptr->static_config.encoder_bit_depth);
                    }
                    frm_hdr->quantization_params.base_q_idx = (uint8_t)CLIP3(
                        (int32_t)quantizer_to_qindex[scs_ptr->static_config.min_qp_allowed],
                        (int32_t)quantizer_to_qindex[scs_ptr->static_config.max_qp_allowed],
                        (int32_t)(new_qindex));

                    pcs_ptr->picture_qp = (uint8_t)CLIP3(
                        (int32_t)scs_ptr->static_config.min_qp_allowed,
                        (int32_t)scs_ptr->static_config.max_qp_allowed,
                        (frm_hdr->quantization_params.base_q_idx + 2) >> 2);
                    // max bit rate is only active for 1 pass CRF
                    if (scs_ptr->static_config.rate_control_mode == 0 &&
                        scs_ptr->static_config.max_bit_rate)
                        crf_assign_max_rate(pcs_ptr->parent_pcs_ptr);
                } else if (pcs_ptr->parent_pcs_ptr->qp_on_the_fly == EB_TRUE) {
                    pcs_ptr->picture_qp = (uint8_t)CLIP3(
                        (int32_t)scs_ptr->static_config.min_qp_allowed,
                        (int32_t)scs_ptr->static_config.max_qp_allowed,
                        pcs_ptr->parent_pcs_ptr->picture_qp);
                    frm_hdr->quantization_params.base_q_idx =
                        quantizer_to_qindex[pcs_ptr->picture_qp];
                }

                pcs_ptr->parent_pcs_ptr->picture_qp = pcs_ptr->picture_qp;
                setup_segmentation(pcs_ptr, scs_ptr);
            } else {
                // ***Rate Control***
                int32_t new_qindex;
                int32_t update_type = scs_ptr->encode_context_ptr->gf_group
                                          .update_type[pcs_ptr->parent_pcs_ptr->gf_group_index];

                if (pcs_ptr->parent_pcs_ptr->tpl_ctrls.enable && pcs_ptr->parent_pcs_ptr->r0 != 0 &&
                    (update_type == KF_UPDATE || update_type == GF_UPDATE ||
                     update_type == ARF_UPDATE)) {
                    process_tpl_stats_frame_kf_gfu_boost(pcs_ptr);
                }
                // Qindex calculating
                if (scs_ptr->static_config.rate_control_mode == 2 &&
                    scs_ptr->static_config.pass == ENC_SINGLE_PASS)
                    new_qindex = rc_pick_q_and_bounds_no_stats_cbr(pcs_ptr);
                else
                    new_qindex = rc_pick_q_and_bounds(pcs_ptr);
                frm_hdr->quantization_params.base_q_idx = (uint8_t)CLIP3(
                    (int32_t)quantizer_to_qindex[scs_ptr->static_config.min_qp_allowed],
                    (int32_t)quantizer_to_qindex[scs_ptr->static_config.max_qp_allowed],
                    (int32_t)(new_qindex));

                pcs_ptr->picture_qp = (uint8_t)CLIP3(
                    (int32_t)scs_ptr->static_config.min_qp_allowed,
                    (int32_t)scs_ptr->static_config.max_qp_allowed,
                    (frm_hdr->quantization_params.base_q_idx + 2) >> 2);

                //Limiting the QP based on the QP of the Reference frame
                if ((int32_t)pcs_ptr->temporal_layer_index != 0) {
                    uint32_t ref_qp = 0;
                    if (pcs_ptr->ref_slice_type_array[0][0] != I_SLICE)
                        ref_qp = pcs_ptr->ref_pic_qp_array[0][0];
                    if ((pcs_ptr->slice_type == B_SLICE) &&
                        (pcs_ptr->ref_slice_type_array[1][0] != I_SLICE))
                        ref_qp = MAX(ref_qp, pcs_ptr->ref_pic_qp_array[1][0]);
                    if (ref_qp > 0 && pcs_ptr->picture_qp < ref_qp) {
                        pcs_ptr->picture_qp = (uint8_t)CLIP3(scs_ptr->static_config.min_qp_allowed,
                                                             scs_ptr->static_config.max_qp_allowed,
                                                             (uint8_t)(ref_qp));
                    }
                }
                frm_hdr->quantization_params.base_q_idx = quantizer_to_qindex[pcs_ptr->picture_qp];
            }
            if (pcs_ptr->parent_pcs_ptr->slice_type == I_SLICE) {
                rate_control_param_ptr->last_i_qp = pcs_ptr->picture_qp;
            } else if (pcs_ptr->parent_pcs_ptr->transition_present) {
                uint32_t min_ref_qp = rate_control_param_ptr->last_i_qp;
                if (pcs_ptr->ref_slice_type_array[0][0] != I_SLICE)
                    min_ref_qp = pcs_ptr->ref_pic_qp_array[0][0];
                if ((pcs_ptr->slice_type == B_SLICE) &&
                    (pcs_ptr->ref_slice_type_array[1][0] != I_SLICE))
                    min_ref_qp = MIN(min_ref_qp, pcs_ptr->ref_pic_qp_array[1][0]);

                pcs_ptr->picture_qp = (uint8_t)CLIP3(
                    scs_ptr->static_config.min_qp_allowed,
                    scs_ptr->static_config.max_qp_allowed,
                    (min_ref_qp + rate_control_param_ptr->last_i_qp) / 2);

                frm_hdr->quantization_params.base_q_idx = quantizer_to_qindex[pcs_ptr->picture_qp];
            }
            pcs_ptr->parent_pcs_ptr->picture_qp = pcs_ptr->picture_qp;

            if (!is_superres_recode_task) {
                // Determine superres parameters for 1-pass encoding or 2nd pass of 2-pass encoding
                // if superres_mode is SUPERRES_QTHRESH or SUPERRES_AUTO.
                // SUPERRES_FIXED and SUPERRES_RANDOM modes are handled in picture decision process.
                if (scs_ptr->static_config.pass == ENC_SINGLE_PASS) {
                    if (scs_ptr->static_config.superres_mode > SUPERRES_RANDOM) {
                        // determine denom and scale down picture by selected denom
                        init_resize_picture(scs_ptr, pcs_ptr->parent_pcs_ptr);
                        if (pcs_ptr->parent_pcs_ptr->frame_superres_enabled) {
                            // reset gm based on super-res on/off
                            set_gm_controls(pcs_ptr->parent_pcs_ptr,
                                            derive_gm_level(pcs_ptr->parent_pcs_ptr));

                            // Initialize Segments as picture decision process
                            pcs_ptr->parent_pcs_ptr->me_segments_completion_count = 0;
                            pcs_ptr->parent_pcs_ptr->me_processed_b64_count       = 0;

                            for (uint32_t segment_index = 0;
                                 segment_index < pcs_ptr->parent_pcs_ptr->me_segments_total_count;
                                 ++segment_index) {
                                // Get Empty Results Object
                                EbObjectWrapper *out_results_wrapper;
                                svt_get_empty_object(
                                    context_ptr->picture_decision_results_output_fifo_ptr,
                                    &out_results_wrapper);

                                PictureDecisionResults *out_results =
                                    (PictureDecisionResults *)out_results_wrapper->object_ptr;
                                out_results->pcs_wrapper_ptr =
                                    pcs_ptr->parent_pcs_ptr->p_pcs_wrapper_ptr;
                                out_results->segment_index = segment_index;
                                out_results->task_type     = TASK_SUPERRES_RE_ME;
                                //Post the Full Results Object
                                svt_post_full_object(out_results_wrapper);
                            }

                            // Release Rate Control Tasks
                            svt_release_object(rate_control_tasks_wrapper_ptr);

                            break;
                        } else {
                            // pa_ref_objs are no longer needed if super-res isn't performed on current frame
                            if (pcs_ptr->parent_pcs_ptr->tpl_ctrls.enable) {
                                if (pcs_ptr->parent_pcs_ptr->temporal_layer_index == 0) {
                                    for (uint32_t i = 0;
                                         i < pcs_ptr->parent_pcs_ptr->tpl_group_size;
                                         i++) {
                                        if (pcs_ptr->parent_pcs_ptr->tpl_group[i]->slice_type ==
                                            P_SLICE) {
                                            if (pcs_ptr->parent_pcs_ptr->tpl_group[i]->ext_mg_id ==
                                                pcs_ptr->parent_pcs_ptr->ext_mg_id + 1) {
                                                release_pa_reference_objects(
                                                    scs_ptr, pcs_ptr->parent_pcs_ptr->tpl_group[i]);
                                            }
                                        } else {
                                            if (pcs_ptr->parent_pcs_ptr->tpl_group[i]->ext_mg_id ==
                                                pcs_ptr->parent_pcs_ptr->ext_mg_id) {
                                                release_pa_reference_objects(
                                                    scs_ptr, pcs_ptr->parent_pcs_ptr->tpl_group[i]);
                                            }
                                        }
                                    }
                                }
                            } else {
                                release_pa_reference_objects(scs_ptr, pcs_ptr->parent_pcs_ptr);
                            }
                        }
                    }
                }
            }

            // QPM with tpl_la
            if (scs_ptr->static_config.enable_adaptive_quantization == 2 &&
                pcs_ptr->parent_pcs_ptr->tpl_ctrls.enable && pcs_ptr->parent_pcs_ptr->r0 != 0) {
                sb_qp_derivation_tpl_la(pcs_ptr);
            } else {
                pcs_ptr->parent_pcs_ptr->frm_hdr.delta_q_params.delta_q_present = 0;
                for (int sb_addr = 0; sb_addr < pcs_ptr->sb_total_count_pix; ++sb_addr) {
                    SuperBlock *sb_ptr = pcs_ptr->sb_ptr_array[sb_addr];
                    sb_ptr->qindex     = quantizer_to_qindex[pcs_ptr->picture_qp];
                }
            }
            if (scs_ptr->static_config.rate_control_mode && !is_superres_recode_task) {
                update_rc_counts(pcs_ptr->parent_pcs_ptr);
            }
            // Get Empty Rate Control Results Buffer
            svt_get_empty_object(context_ptr->rate_control_output_results_fifo_ptr,
                                 &rate_control_results_wrapper_ptr);
            rate_control_results_ptr = (RateControlResults *)
                                           rate_control_results_wrapper_ptr->object_ptr;
            rate_control_results_ptr->pcs_wrapper_ptr = rate_control_tasks_ptr->pcs_wrapper_ptr;
            rate_control_results_ptr->superres_recode = is_superres_recode_task;

            // Post Full Rate Control Results
            svt_post_full_object(rate_control_results_wrapper_ptr);

            // Release Rate Control Tasks
            svt_release_object(rate_control_tasks_wrapper_ptr);

            break;

        case RC_PACKETIZATION_FEEDBACK_RESULT:

            parentpicture_control_set_ptr = (PictureParentControlSet *)
                                                rate_control_tasks_ptr->pcs_wrapper_ptr->object_ptr;
            scs_ptr = (SequenceControlSet *)
                          parentpicture_control_set_ptr->scs_wrapper_ptr->object_ptr;

            if (scs_ptr->static_config.intra_period_length != -1 &&
                scs_ptr->static_config.rate_control_mode != 0) {
                uint32_t interval_index_temp = 0;
                for (RateControlIntervalParamContext **rc_param_queue =
                         context_ptr->rate_control_param_queue;
                     interval_index_temp < PARALLEL_GOP_MAX_NUMBER &&
                     (rc_param_queue[interval_index_temp]->first_poc >
                          parentpicture_control_set_ptr->picture_number ||
                      rc_param_queue[interval_index_temp]->last_poc <
                          parentpicture_control_set_ptr->picture_number);
                     ++interval_index_temp) {}
                CHECK_REPORT_ERROR(interval_index_temp != PARALLEL_GOP_MAX_NUMBER,
                                   scs_ptr->encode_context_ptr->app_callback_ptr,
                                   EB_ENC_RC_ERROR2);
                rate_control_param_ptr = context_ptr->rate_control_param_queue[interval_index_temp];
                rate_control_param_ptr->processed_frame_number++;
                // check if all the frames in the interval have arrived
                if (scs_ptr->static_config.intra_period_length + 1 ==
                    rate_control_param_ptr->processed_frame_number) {
                    uint64_t shift_count = (uint64_t)(scs_ptr->static_config.intra_period_length +
                                                      1) *
                        PARALLEL_GOP_MAX_NUMBER;
                    rate_control_param_ptr->first_poc += shift_count;
                    rate_control_param_ptr->last_poc += shift_count;
                    rate_control_param_ptr->processed_frame_number = 0;
                }
            }
            if (scs_ptr->static_config.rate_control_mode) {
                restore_gf_group_param(parentpicture_control_set_ptr);
                av1_rc_postencode_update(parentpicture_control_set_ptr);
                // Qindex calculating
                if (!(scs_ptr->static_config.rate_control_mode == 2 &&
                      scs_ptr->static_config.pass == ENC_SINGLE_PASS))
                    svt_av1_twopass_postencode_update(parentpicture_control_set_ptr);
            }
            // Queue variables
            if (scs_ptr->static_config.max_bit_rate)
                coded_frames_stat_calc(parentpicture_control_set_ptr);
            total_number_of_fb_frames++;

            // Release the SequenceControlSet
            svt_release_object(parentpicture_control_set_ptr->scs_wrapper_ptr);
            // Release the ParentPictureControlSet
            //y8b needs to get decremented at the same time of regular input
            //  svt_release_object_with_call_stack(parentpicture_control_set_ptr->eb_y8b_wrapper_ptr, 3000, parentpicture_control_set_ptr->picture_number);
            svt_release_object(parentpicture_control_set_ptr->eb_y8b_wrapper_ptr);

            svt_release_object(parentpicture_control_set_ptr->input_picture_wrapper_ptr);
            svt_release_object(rate_control_tasks_ptr->pcs_wrapper_ptr);

            // Release Rate Control Tasks
            svt_release_object(rate_control_tasks_wrapper_ptr);
            break;

        default:
            pcs_ptr = (PictureControlSet *)rate_control_tasks_ptr->pcs_wrapper_ptr->object_ptr;
            scs_ptr = (SequenceControlSet *)pcs_ptr->scs_wrapper_ptr->object_ptr;

            break;
        }
    }

    return NULL;
}
