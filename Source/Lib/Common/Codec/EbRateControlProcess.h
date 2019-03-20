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
#define PPS_BITS_NUM                 80
#define SPS_BITS_NUM                296
#define RC_PRECISION                16
#define RC_PRECISION_OFFSET         (1 << (RC_PRECISION - 1))

#define OVERSHOOT_STAT_PRINT             0 // Do not remove.
                                           // For printing overshooting percentages for both RC and fixed QP.
                                           // Target rate and and max buffer size should be set properly even for fixed QP.
                                           // Disabled by default.
#if OVERSHOOT_STAT_PRINT
#define CODED_FRAMES_STAT_QUEUE_MAX_DEPTH   10000
#endif

#define ADAPTIVE_PERCENTAGE 1
#define     RC_INTRA_QP_OFFSET              (-1)

#define     RC_DISABLE_FLRC                 0
#define     RC_DISABLE_FLRC_RATE_UPDATE     0
#define     RC_DISABLE_HLRC_RATE_INPUT      1
#define     RC_DISABLE_PRED_QP_LIMIT        1
#define     RC_DISABLE_EXTRA_BITS           0
#define     RC_NEW_EXTRA_BITS               1
#define     RC_UPDATE_TARGET_RATE           1

#define        RC_QPMOD_MAXQP                    42

static const uint32_t  rate_percentage_layer_array[EB_MAX_TEMPORAL_LAYERS][EB_MAX_TEMPORAL_LAYERS] = {
    {100,  0,  0,  0,  0,  0 },
    { 70, 30,  0,  0,  0,  0 },
    { 70, 15, 15,  0,  0,  0 },
    { 55, 15, 15, 15,  0,  0 },
    { 40, 15, 15, 15, 15,  0 },
    { 30, 10, 15, 15, 15, 15 }
};

// range from 0 to 51
// precision is 16 bits
static const uint64_t two_to_power_qp_over_three[] = {
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
// range is from -51 to 51 (0 to 102)
static const uint64_t two_to_power_x_over_six[] = {
         0xB5,       0xCB,        0xE4,        0xFF,        0x11F,        0x142,
        0x16A,      0x196,       0x1C8,       0x1FF,        0x23E,        0x285,
        0x2D4,      0x32C,       0x390,       0x3FF,        0x47D,        0x50A,
        0x5A8,      0x659,       0x720,       0x7FF,        0x8FA,        0xA14,
        0xB50,      0xCB2,       0xE41,       0xFFF,       0x11F5,       0x1428,
       0x16A0,     0x1965,      0x1C82,      0x1FFF,       0x23EB,       0x2851,
       0x2D41,     0x32CB,      0x3904,      0x3FFF,       0x47D6,       0x50A2,
       0x5A82,     0x6597,      0x7208,      0x7FFF,       0x8FAC,       0xA144,
       0xB504,     0xCB2F,      0xE411,      0xFFFF,      0x11F58,      0x14288,
      0x16A08,    0x1965E,     0x1C822,     0x1FFFE,      0x23EB1,    0x28511,
      0x2D410,    0x32CBC,     0x39044,     0x3FFFC,      0x47D62,    0x50A23,
      0x5A821,    0x65979,     0x72088,     0x7FFF8,      0x8FAC4,    0xA1447,
      0xB5043,    0xCB2F2,     0xE4110,     0xFFFF0,     0x11F588,   0x14288E,
     0x16A087,   0x1965E5,    0x1C8221,    0x1FFFE0,     0x23EB11,   0x28511D,
     0x2D410F,   0x32CBCA,    0x390443,    0x3FFFC0,     0x47D623,   0x50A23B,
     0x5A821F,   0x659794,    0x720886,    0x7FFF80,     0x8FAC46,   0xA14476,
     0xB5043E,   0xCB2F29,    0xE4110C,    0xFFFF00,    0x11F588C,  0x14288ED,
     0x16A087C
};
/**************************************
 * Input Port Types
 **************************************/
typedef enum RATE_CONTROL_INPUT_PORT_TYPES {
    RATE_CONTROL_INPUT_PORT_PICTURE_MANAGER = 0,
    RATE_CONTROL_INPUT_PORT_PACKETIZATION = 1,
    RATE_CONTROL_INPUT_PORT_ENTROPY_CODING = 2,
    RATE_CONTROL_INPUT_PORT_TOTAL_COUNT = 3,
    RATE_CONTROL_INPUT_PORT_INVALID = ~0,
} RATE_CONTROL_INPUT_PORT_TYPES;

/**************************************
 * Input Port Config
 **************************************/
typedef struct RateControlPorts_s {
    RATE_CONTROL_INPUT_PORT_TYPES    type;
    uint32_t                           count;
} RateControlPorts_t;

/**************************************
 * Coded Frames Stats
 **************************************/
typedef struct CodedFramesStatsEntry_s {
    uint64_t               picture_number;
    int64_t               frameTotalBitActual;
    EbBool              end_of_sequence_flag;
} CodedFramesStatsEntry_t;
/**************************************
 * Context
 **************************************/
typedef struct RateControlLayerContext_s
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

} RateControlLayerContext_t;

typedef struct RateControlIntervalParamContext_s
{
    uint64_t                       firstPoc;
    uint64_t                       lastPoc;
    EbBool                      inUse;
    EbBool                      wasUsed;
    uint64_t                       processedFramesNumber;
    EbBool                      lastGop;
    RateControlLayerContext_t  **rateControlLayerArray;

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

} RateControlIntervalParamContext_t;

typedef struct HighLevelRateControlContext_s
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


} HighLevelRateControlContext_t;

typedef struct RateControlContext_s
{
    EbFifo_t                    *rate_control_input_tasks_fifo_ptr;
    EbFifo_t                    *rate_control_output_results_fifo_ptr;

    HighLevelRateControlContext_t *highLevelRateControlPtr;

    RateControlIntervalParamContext_t **rateControlParamQueue;
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
    CodedFramesStatsEntry_t    **codedFramesStatQueue;
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


} RateControlContext_t;

/**************************************
 * Extern Function Declarations
 **************************************/
extern EbErrorType rate_control_layer_context_ctor(
    RateControlLayerContext_t **entry_dbl_ptr);

extern EbErrorType rate_control_interval_param_context_ctor(
    RateControlIntervalParamContext_t **entry_dbl_ptr);

extern EbErrorType rate_control_coded_frames_stats_context_ctor(
    CodedFramesStatsEntry_t **entry_dbl_ptr,
    uint64_t                  picture_number);

extern EbErrorType rate_control_context_ctor(
    RateControlContext_t **context_dbl_ptr,
    EbFifo_t              *rate_control_input_tasks_fifo_ptr,
    EbFifo_t              *rate_control_output_results_fifo_ptr,
    int32_t                intra_period_length);

extern void* rate_control_kernel(void *input_ptr);

#endif // EbRateControl_h