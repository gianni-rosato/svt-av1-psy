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

// Threshold used to define if a KF group is static (e.g. a slide show).
// Essentially, this means that no frame in the group has more than 1% of MBs
// that are not marked as coded with 0,0 motion in the first pass.
#define STATIC_KF_GROUP_THRESH 99
#define STATIC_KF_GROUP_FLOAT_THRESH 0.99

// Minimum and maximum height for the new pyramid structure.
// (Old structure supports height = 1, but does NOT support height = 4).
#define MIN_PYRAMID_LVL 0
#define MAX_PYRAMID_LVL 4

#define MIN_GF_INTERVAL 4
#define MAX_GF_INTERVAL 32
#define FIXED_GF_INTERVAL 8  // Used in some testing modes only
#define MAX_GF_LENGTH_LAP 16

#define MAX_NUM_GF_INTERVALS 15

#define MAX_ARF_LAYERS 6

enum {
    KF_UPDATE,
    LF_UPDATE,
    GF_UPDATE,
    ARF_UPDATE,
    OVERLAY_UPDATE,
    INTNL_OVERLAY_UPDATE,  // Internal Overlay Frame
    INTNL_ARF_UPDATE,      // Internal Altref Frame
    FRAME_UPDATE_TYPES
} UENUM1BYTE(FRAME_UPDATE_TYPE);

typedef enum rate_factor_level {
    INTER_NORMAL       = 0,
    INTER_LOW          = 1,
    INTER_HIGH         = 2,
    GF_ARF_LOW         = 3,
    GF_ARF_STD         = 4,
    KF_STD             = 5,
    RATE_FACTOR_LEVELS = 6
} rate_factor_level;

typedef struct {
    // Rate targetting variables
    int base_frame_target; // A baseline frame target before adjustment
        // for previous under or over shoot.
    int this_frame_target; // Actual frame target after rc adjustment.
    int projected_frame_size;
    int sb64_target_rate;
    int last_q[FRAME_TYPES]; // Separate values for Intra/Inter
    int last_boosted_qindex; // Last boosted GF/KF/ARF q
    int last_kf_qindex; // Q index of the last key frame coded.

    int gfu_boost;
    int kf_boost;

    double rate_correction_factors[RATE_FACTOR_LEVELS];

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

    int avg_frame_bandwidth; // Average frame size target for clip
    int min_frame_bandwidth; // Minimum allocation used for any frame
    int max_frame_bandwidth; // Maximum burst rate allowed for a frame.

    int    ni_av_qi;
    int    ni_tot_qi;
    int    ni_frames;
    int    avg_frame_qindex[FRAME_TYPES];
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

    // real for TWOPASS_RC
    int prev_avg_frame_bandwidth; //only for CBR?
    int active_worst_quality;
    int active_best_quality[MAX_ARF_LAYERS + 1];
    int base_layer_qp;

    // number of determined gf group length left
    int intervals_till_gf_calculate_due;
    // stores gf group length intervals
    int gf_intervals[MAX_NUM_GF_INTERVALS];
    // the current index in gf_intervals
    int cur_gf_index;

    // gop bit budget
    int64_t gf_group_bits;

    // Total number of stats used only for kf_boost calculation.
    int num_stats_used_for_kf_boost;
    // Total number of stats used only for gfu_boost calculation.
    int num_stats_used_for_gfu_boost;
    // Total number of stats required by gfu_boost calculation.
    int num_stats_required_for_gfu_boost;
    int enable_scenecut_detection;
    int use_arf_in_this_kf_group;
    int next_is_fwd_key;
} RATE_CONTROL;

/**************************************
 * Input Port Types
 **************************************/
typedef enum RateControlInputPortTypes {
#if FEATURE_INL_ME
    RATE_CONTROL_INPUT_PORT_INLME = 0,
#else
    RATE_CONTROL_INPUT_PORT_PICTURE_MANAGER = 0,
#endif
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
double svt_av1_convert_qindex_to_q(int32_t qindex, AomBitDepth bit_depth);
int svt_av1_rc_get_default_min_gf_interval(int width, int height, double framerate);
int svt_av1_rc_get_default_max_gf_interval(double framerate, int min_gf_interval);
double svt_av1_get_gfu_boost_projection_factor(double min_factor, double max_factor, int frame_count);

EbErrorType rate_control_context_ctor(EbThreadContext *  thread_context_ptr,
                                      const EbEncHandle *enc_handle_ptr);

extern void *rate_control_kernel(void *input_ptr);

#endif // EbRateControl_h
