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
#include "EbObject.h"

#define CCOEFF_INIT_FACT 2
#define SAD_CLIP_COEFF 5
// 88 + 3*16*8
#define SLICE_HEADER_BITS_NUM 104
#define RC_PRECISION 16
#define RC_PRECISION_OFFSET (1 << (RC_PRECISION - 1))

#define OVERSHOOT_STAT_PRINT 0
/* Do not remove
 * For printing overshooting percentages for both RC and fixed QP.
 * Target rate and and max buffer size should be set properly even for fixed QP.
 * Disabled by default.
*/
#if OVERSHOOT_STAT_PRINT
#define CODED_FRAMES_STAT_QUEUE_MAX_DEPTH 10000
#endif
#define RC_PRINTS 0
#define ADAPTIVE_PERCENTAGE 1

#define RC_QPMOD_MAXQP 54

/**************************************
 * Input Port Types
 **************************************/
typedef enum RateControlInputPortTypes {
    RATE_CONTROL_INPUT_PORT_PICTURE_MANAGER = 0,
    RATE_CONTROL_INPUT_PORT_PACKETIZATION   = 1,
    RATE_CONTROL_INPUT_PORT_ENTROPY_CODING  = 2,
    RATE_CONTROL_INPUT_PORT_TOTAL_COUNT     = 3,
    RATE_CONTROL_INPUT_PORT_INVALID         = ~0,
} RateControlInputPortTypes;

/**************************************
 * Input Port Config
 **************************************/
typedef struct RateControlPorts {
    RateControlInputPortTypes type;
    uint32_t                  count;
} RateControlPorts;

/**************************************
 * Context
 **************************************/

typedef struct RateControlLayerContext {
    EbDctor  dctor;
    uint64_t previous_frame_distortion_me;
    uint64_t previous_frame_bit_actual;
    uint64_t previous_framequantized_coeff_bit_actual;
    EbBool   feedback_arrived;

    uint64_t target_bit_rate;
    uint64_t frame_rate;
    uint64_t channel_bit_rate;

    uint64_t previous_bit_constraint;
    uint64_t bit_constraint;
    uint64_t ec_bit_constraint;
    uint64_t previous_ec_bits;
    int64_t  dif_total_and_ec_bits;

    int64_t  bit_diff;
    uint32_t coeff_averaging_weight1;
    uint32_t coeff_averaging_weight2; // coeff_averaging_weight2 = 16- coeff_averaging_weight1
    //Ccoeffs have 2*RC_PRECISION precision
    int64_t c_coeff;
    int64_t previous_c_coeff;
    //Kcoeffs have RC_PRECISION precision
    uint64_t k_coeff;
    uint64_t previous_k_coeff;

    //delta_qp_fraction has RC_PRECISION precision
    int64_t  delta_qp_fraction;
    uint32_t previous_frame_qp;
    uint32_t calculated_frame_qp;
    uint32_t previous_calculated_frame_qp;
    uint32_t area_in_pixel;
    uint32_t previous_frame_average_qp;

    //total_mad has RC_PRECISION precision
    uint64_t total_mad;

    uint32_t first_frame;
    uint32_t first_non_intra_frame;
    uint32_t same_distortion_count;
    uint32_t frame_same_distortion_min_qp_count;
    uint32_t critical_states;

    uint32_t max_qp;
    uint32_t temporal_index;

    uint64_t alpha;
    //segmentation
    int32_t prev_segment_qps[MAX_SEGMENTS];

} RateControlLayerContext;

/**************************************
 * Extern Function Declarations
 **************************************/
EbErrorType rate_control_context_ctor(EbThreadContext *  thread_context_ptr,
                                      const EbEncHandle *enc_handle_ptr);

extern void *rate_control_kernel(void *input_ptr);

#endif // EbRateControl_h
