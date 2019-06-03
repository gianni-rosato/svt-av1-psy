/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbRateControl_h
#define EbRateControl_h

#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"
#include "EbSvtAv1Enc.h"
#include "EbPictureControlSet.h"

#define CCOEFF_INIT_FACT              2
#define SAD_CLIP_COEFF                5
// 88 + 3*16*8
#define SLICE_HEADER_BITS_NUM       104
#define RC_PRECISION                16
#define RC_PRECISION_OFFSET         (1 << (RC_PRECISION - 1))

#define OVERSHOOT_STAT_PRINT             0 // Do not remove.
                                           // For printing overshooting percentages for both RC and fixed QP.
                                           // Target rate and and max buffer size should be set properly even for fixed QP.
                                           // Disabled by default.
#if OVERSHOOT_STAT_PRINT
#define CODED_FRAMES_STAT_QUEUE_MAX_DEPTH   10000
#endif
#define RC_PRINTS                   0
#define ADAPTIVE_PERCENTAGE         1
#define RC_UPDATE_TARGET_RATE       1

#define RC_QPMOD_MAXQP             54

static const uint32_t  rate_percentage_layer_array[EB_MAX_TEMPORAL_LAYERS][EB_MAX_TEMPORAL_LAYERS] =
{
    {100,  0,  0,  0,  0,  0 },
    { 70, 30,  0,  0,  0,  0 },
    { 70, 15, 15,  0,  0,  0 },
    { 55, 15, 15, 15,  0,  0 },
    { 40, 15, 15, 15, 15,  0 },
    { 30, 10, 15, 15, 15, 15 }
};

// range from 0 to 51
// precision is 16 bits
static const uint64_t two_to_power_qp_over_three[] =
{
         0x10000,      0x1428A,     0x19660,     0x20000,
         0x28514,      0x32CC0,     0x40000,     0x50A29,
         0x65980,      0x80000,     0xA1451,     0xCB2FF,
        0x100000,     0x1428A3,    0x1965FF,    0x200000,
        0x285146,     0x32CBFD,    0x400000,    0x50A28C,
        0x6597FB,     0x800000,    0xA14518,    0xCB2FF5,
       0x1000000,    0x1428A30,   0x1965FEA,   0x2000000,
       0x285145F,    0x32CBFD5,   0x4000000,   0x50A28BE,
       0x6597FA9,    0x8000000,   0xA14517D,   0xCB2FF53,
      0x10000000,   0x1428A2FA,  0x1965FEA5,  0x20000000,
      0x285145F3,   0x32CBFD4A,  0x40000000,  0x50A28BE6,
      0x6597FA95,   0x80000000,  0xA14517CC,  0xCB2FF52A,
     0x100000000,  0x1428A2F99, 0x1965FEA54, 0x200000000
};
/**************************************
 * Input Port Types
 **************************************/
typedef enum RateControlInputPortTypes
{
    RATE_CONTROL_INPUT_PORT_PICTURE_MANAGER = 0,
    RATE_CONTROL_INPUT_PORT_PACKETIZATION = 1,
    RATE_CONTROL_INPUT_PORT_ENTROPY_CODING = 2,
    RATE_CONTROL_INPUT_PORT_TOTAL_COUNT = 3,
    RATE_CONTROL_INPUT_PORT_INVALID = ~0,
} RateControlInputPortTypes;

/**************************************
 * Input Port Config
 **************************************/
typedef struct RateControlPorts
{
    RateControlInputPortTypes    type;
    uint32_t                           count;
} RateControlPorts;

/**************************************
 * Coded Frames Stats
 **************************************/
typedef struct CodedFramesStatsEntry
{
    uint64_t               picture_number;
    int64_t               frame_total_bit_actual;
    EbBool              end_of_sequence_flag;
} CodedFramesStatsEntry;
/**************************************
 * Context
 **************************************/

#if RC
typedef struct RateControlLayerContext
{
    uint64_t   previous_frame_distortion_me;
    uint64_t   previous_frame_bit_actual;
    uint64_t   previous_framequantized_coeff_bit_actual;
    EbBool     feedback_arrived;

    uint64_t   target_bit_rate;
    uint64_t   frame_rate;
    uint64_t   channel_bit_rate;

    uint64_t   previous_bit_constraint;
    uint64_t   bit_constraint;
    uint64_t   ec_bit_constraint;
    uint64_t   previous_ec_bits;
    int64_t    dif_total_and_ec_bits;

    int64_t    bit_diff;
    uint32_t   coeff_averaging_weight1;
    uint32_t   coeff_averaging_weight2; // coeff_averaging_weight2 = 16- coeff_averaging_weight1
    //Ccoeffs have 2*RC_PRECISION precision
    int64_t    c_coeff;
    int64_t    previous_c_coeff;
    //Kcoeffs have RC_PRECISION precision
    uint64_t   k_coeff;
    uint64_t   previous_k_coeff;

    //delta_qp_fraction has RC_PRECISION precision
    int64_t    delta_qp_fraction;
    uint32_t   previous_frame_qp;
    uint32_t   calculated_frame_qp;
    uint32_t   previous_calculated_frame_qp;
    uint32_t   area_in_pixel;
    uint32_t   previous_frame_average_qp;

    //total_mad has RC_PRECISION precision
    uint64_t   total_mad;

    uint32_t   first_frame;
    uint32_t   first_non_intra_frame;
    uint32_t   same_distortion_count;
    uint32_t   frame_same_distortion_min_qp_count;
    uint32_t   critical_states;

    uint32_t   max_qp;
    uint32_t   temporal_index;

    uint64_t   alpha;
} RateControlLayerContext;

typedef struct RateControlIntervalParamContext
{
    uint64_t                     first_poc;
    uint64_t                     last_poc;
    EbBool                       in_use;
    EbBool                       was_used;
    uint64_t                     processed_frames_number;
    EbBool                       last_gop;
    RateControlLayerContext   **rate_control_layer_array;

    int64_t                      virtual_buffer_level;
    int64_t                      previous_virtual_buffer_level;
    uint32_t                     intra_frames_qp;
    uint8_t                      intra_frames_qp_bef_scal;

    uint32_t                     next_gop_intra_frame_qp;
    uint64_t                     first_pic_pred_bits;
    uint64_t                     first_pic_actual_bits;
    uint16_t                     first_pic_pred_qp;
    uint16_t                     first_pic_actual_qp;
    EbBool                       first_pic_actual_qp_assigned;
    EbBool                       scene_change_in_gop;
    EbBool                       min_target_rate_assigned;
    int64_t                      extra_ap_bit_ratio_i;
} RateControlIntervalParamContext;

typedef struct HighLevelRateControlContext
{
    uint64_t target_bit_rate;
    uint64_t frame_rate;
    uint64_t channel_bit_rate_per_frame;
    uint64_t channel_bit_rate_per_sw;
    uint64_t bit_constraint_per_sw;
    uint64_t pred_bits_ref_qpPerSw[MAX_REF_QP_NUM];
#if RC_UPDATE_TARGET_RATE
    uint32_t prev_intra_selected_ref_qp;
    uint32_t prev_intra_org_selected_ref_qp;
    uint64_t previous_updated_bit_constraint_per_sw;
#endif
} HighLevelRateControlContext;

typedef struct RateControlContext
{
    EbFifo                            *rate_control_input_tasks_fifo_ptr;
    EbFifo                            *rate_control_output_results_fifo_ptr;

    HighLevelRateControlContext       *high_level_rate_control_ptr;

    RateControlIntervalParamContext  **rate_control_param_queue;
    uint64_t                           rate_control_param_queue_head_index;

    uint64_t                           frame_rate;

    uint64_t                           virtual_buffer_size;

    int64_t                            virtual_buffer_level_initial_value;
    int64_t                            previous_virtual_buffer_level;

    int64_t                            virtual_buffer_level;

    //Virtual Buffer Thresholds
    int64_t                            vb_fill_threshold1;
    int64_t                            vb_fill_threshold2;

    // Rate Control Previous Bits Queue
#if OVERSHOOT_STAT_PRINT
    CodedFramesStatsEntry            **coded_frames_stat_queue;
    uint32_t                           coded_frames_stat_queue_head_index;
    uint32_t                           coded_frames_stat_queue_tail_index;

    uint64_t                           total_bit_actual_per_sw;
    uint64_t                           max_bit_actual_per_sw;
    uint64_t                           max_bit_actual_per_gop;
    uint64_t                           min_bit_actual_per_gop;
    uint64_t                           avg_bit_actual_per_gop;

#endif

    uint64_t                           rate_average_periodin_frames;
    uint32_t                           base_layer_frames_avg_qp;
    uint32_t                           base_layer_intra_frames_avg_qp;

    EbBool                            end_of_sequence_region;

    uint32_t                           intra_coef_rate;

    uint64_t                           frames_in_interval[EB_MAX_TEMPORAL_LAYERS];
    int64_t                            extra_bits;
    int64_t                            extra_bits_gen;
    int16_t                            max_rate_adjust_delta_qp;

    uint32_t                           qp_scaling_map[EB_MAX_TEMPORAL_LAYERS][MAX_REF_QP_NUM];
    uint32_t                           qp_scaling_map_I_SLICE[MAX_REF_QP_NUM];
} RateControlContext;
#else
typedef struct RateControlLayerContext
{
    uint64_t                  previousFrameSadMe;
    uint64_t                  previousFrameBitActual;
    uint64_t                  previousFrameQuantizedCoeffBitActual;
    EbBool                 feedbackArrived;

    uint64_t                  target_bit_rate;
    uint64_t                  frame_rate;
    uint64_t                  channelBitRate;

    uint64_t                  previousBitConstraint;
    uint64_t                  bitConstraint;
    uint64_t                  ecBitConstraint;
    uint64_t                  previousEcBits;
    int64_t                  difTotalAndEcBits;
    int64_t                  prevDifTotalAndEcBits;

    int64_t                  bitDiff;
    uint32_t                  coeffAveragingWeight1;
    uint32_t                  coeffAveragingWeight2; // coeffAveragingWeight2 = 16- coeffAveragingWeight1
    //Ccoeffs have 2*RC_PRECISION precision
    int64_t                  cCoeff;
    int64_t                  previousCCoeff;
    //Kcoeffs have RC_PRECISION precision
    uint64_t                  kCoeff;
    uint64_t                  previousKCoeff;
    uint64_t                  coeffWeight;

    //deltaQpFraction has RC_PRECISION precision
    int64_t                  deltaQpFraction;
    uint32_t                  previousFrameQp;
    uint32_t                  calculatedFrameQp;
    uint32_t                  previousCalculatedFrameQp;
    uint32_t                  areaInPixel;
    uint32_t                  previousFrameAverageQp;

    //totalMad has RC_PRECISION precision
    uint64_t                  totalMad;

    uint32_t                  firstFrame;
    uint32_t                  firstNonIntraFrame;
    uint32_t                  sameSADCount;
    uint32_t                  frameSameSADMinQpCount;
    uint32_t                  criticalStates;

    uint32_t                  maxQp;
    uint32_t                  temporalIndex;

    uint64_t                  alpha;
} RateControlLayerContext;

typedef struct RateControlIntervalParamContext
{
    uint64_t                       firstPoc;
    uint64_t                       lastPoc;
    EbBool                      inUse;
    EbBool                      wasUsed;
    uint64_t                       processedFramesNumber;
    EbBool                      lastGop;
    RateControlLayerContext  **rateControlLayerArray;

    int64_t                       virtualBufferLevel;
    int64_t                       previousVirtualBufferLevel;
    uint32_t                       intraFramesQp;

    uint32_t                       nextGopIntraFrameQp;
    int64_t                       totalExtraBits;
    uint64_t                       firstPicPredBits;
    uint64_t                       firstPicActualBits;
    uint16_t                       firstPicPredQp;
    uint16_t                       firstPicActualQp;
    EbBool                       firstPicActualQpAssigned;
    EbBool                      scene_change_in_gop;
    EbBool                      min_target_rate_assigned;
    int64_t                       extraApBitRatioI;
} RateControlIntervalParamContext;

typedef struct HighLevelRateControlContext
{
    uint64_t                       targetBitsPerSlidingWindow;
    uint64_t                       target_bit_rate;
    uint64_t                       frame_rate;
    uint64_t                       channelBitRatePerFrame;
    uint64_t                       channelBitRatePerSw;
    uint64_t                       bitConstraintPerSw;
    uint64_t                       predBitsRefQpPerSw[MAX_REF_QP_NUM];
#if RC_UPDATE_TARGET_RATE
    uint32_t                       prevIntraSelectedRefQp;
    uint32_t                       prevIntraOrgSelectedRefQp;
    uint64_t                       previousUpdatedBitConstraintPerSw;
#endif
} HighLevelRateControlContext;

typedef struct RateControlContext
{
    EbFifo                    *rate_control_input_tasks_fifo_ptr;
    EbFifo                    *rate_control_output_results_fifo_ptr;

    HighLevelRateControlContext *highLevelRateControlPtr;

    RateControlIntervalParamContext **rateControlParamQueue;
    uint64_t                       rateControlParamQueueHeadIndex;

    uint64_t                       frame_rate;

    uint64_t                       virtualBufferSize;

    int64_t                       virtualBufferLevelInitialValue;
    int64_t                       previousVirtualBufferLevel;

    int64_t                       virtualBufferLevel;

    //Virtual Buffer Thresholds
    int64_t                       vbFillThreshold1;
    int64_t                       vbFillThreshold2;

    // Rate Control Previous Bits Queue
#if OVERSHOOT_STAT_PRINT
    CodedFramesStatsEntry    **codedFramesStatQueue;
    uint32_t                       codedFramesStatQueueHeadIndex;
    uint32_t                       codedFramesStatQueueTailIndex;

    uint64_t                       totalBitActualPerSw;
    uint64_t                       maxBitActualPerSw;
    uint64_t                       maxBitActualPerGop;

#endif

    uint64_t                       rateAveragePeriodinFrames;
    uint32_t                       baseLayerFramesAvgQp;
    uint32_t                       baseLayerIntraFramesAvgQp;

    uint32_t                       baseLayerIntraFramesAvgQpFloat;
    EbBool                      end_of_sequence_region;

    uint32_t                       intraCoefRate;
    uint32_t                       nonPeriodicIntraCoefRate;

    uint64_t                       frames_in_interval[EB_MAX_TEMPORAL_LAYERS];
    int64_t                       extraBits;
    int64_t                       extraBitsGen;
    int16_t                      maxRateAdjustDeltaQP;
} RateControlContext;
#endif
/**************************************
 * Extern Function Declarations
 **************************************/
extern EbErrorType rate_control_layer_context_ctor(
    RateControlLayerContext **entry_dbl_ptr);

extern EbErrorType rate_control_interval_param_context_ctor(
    RateControlIntervalParamContext **entry_dbl_ptr);

extern EbErrorType rate_control_coded_frames_stats_context_ctor(
    CodedFramesStatsEntry **entry_dbl_ptr,
    uint64_t                  picture_number);

extern EbErrorType rate_control_context_ctor(
    RateControlContext **context_dbl_ptr,
    EbFifo              *rate_control_input_tasks_fifo_ptr,
    EbFifo              *rate_control_output_results_fifo_ptr,
    int32_t                intra_period_length);

extern void* rate_control_kernel(void *input_ptr);

#endif // EbRateControl_h
