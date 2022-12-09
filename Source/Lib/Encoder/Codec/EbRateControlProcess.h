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
#include "EbInvTransforms.h"

#define MINQ_ADJ_LIMIT 48
#define HIGH_UNDERSHOOT_RATIO 2
#define CCOEFF_INIT_FACT 2
#define SAD_CLIP_COEFF 5
// 88 + 3*16*8
#define SLICE_HEADER_BITS_NUM 104
#define RC_PRECISION 16
#define RC_PRECISION_OFFSET (1 << (RC_PRECISION - 1))

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
#define FIXED_GF_INTERVAL 8 // Used in some testing modes only
#define MAX_GF_LENGTH_LAP 16
#define MAX_ARF_LAYERS 6

enum {
    KF_UPDATE,
    LF_UPDATE,
    GF_UPDATE,
    ARF_UPDATE,
    OVERLAY_UPDATE,
    INTNL_OVERLAY_UPDATE, // Internal Overlay Frame
    INTNL_ARF_UPDATE, // Internal Altref Frame
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
#define CODED_FRAMES_STAT_QUEUE_MAX_DEPTH 2000
// max bit rate average period
#define MAX_RATE_AVG_PERIOD (CODED_FRAMES_STAT_QUEUE_MAX_DEPTH >> 1)
#define CRITICAL_BUFFER_LEVEL 15
#define OPTIMAL_BUFFER_LEVEL 70
/**************************************
 * Coded Frames Stats
 **************************************/
typedef struct coded_frames_stats_entry {
    EbDctor  dctor;
    uint64_t picture_number;
    int64_t  frame_total_bit_actual;
    Bool     end_of_sequence_flag;
} coded_frames_stats_entry;

typedef enum {
    NO_RESIZE      = 0,
    DOWN_THREEFOUR = 1, // From orig to 3/4.
    DOWN_ONEHALF   = 2, // From orig or 3/4 to 1/2.
    UP_THREEFOUR   = -1, // From 1/2 to 3/4.
    UP_ORIG        = -2, // From 1/2 or 3/4 to orig.
} RESIZE_ACTION;

typedef enum { ORIG = 0, THREE_QUARTER = 1, ONE_HALF = 2 } RESIZE_STATE;

/*!
 * \brief Desired dimensions for an externally triggered resize.
 *
 * When resize is triggered externally, the desired dimensions are stored in
 * this struct until used in the next frame to be coded. These values are
 * effective only for one frame and are reset after they are used.
 */
typedef struct {
    RESIZE_STATE resize_state;
    uint8_t      resize_denom;
} ResizePendingParams;

extern EbErrorType rate_control_coded_frames_stats_context_ctor(coded_frames_stats_entry *entry_ptr,
                                                                uint64_t picture_number);
typedef struct {
    int     last_boosted_qindex; // Last boosted GF/KF/ARF q
    int     gfu_boost;
    int     kf_boost;
    double  rate_correction_factors[MAX_TEMPORAL_LAYERS + 1];
    int     onepass_cbr_mode; // 0: not 1pass cbr, 1: 1pass cbr for low delay
    int     baseline_gf_interval;
    int     constrained_gf_group;
    int     frames_to_key;
    int     frames_since_key;
    int     this_key_frame_forced;
    int     avg_frame_bandwidth; // Average frame size target for clip
    int     max_frame_bandwidth; // Maximum burst rate allowed for a frame.
    int     avg_frame_qindex[FRAME_TYPES];
    int64_t buffer_level;
    int64_t bits_off_target;
    int64_t vbr_bits_off_target;
    int64_t vbr_bits_off_target_fast;
    int     rolling_target_bits;
    int     rolling_actual_bits;
    int     rate_error_estimate;

    int64_t total_actual_bits;
    int64_t total_target_bits;

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

    // gop bit budget
    int64_t gf_group_bits;
    // Total number of stats used only for gfu_boost calculation.
    int num_stats_used_for_gfu_boost;
    // Total number of stats required by gfu_boost calculation.
    int num_stats_required_for_gfu_boost;
    // Rate Control stat Queue
    coded_frames_stats_entry **coded_frames_stat_queue;
    uint32_t                   coded_frames_stat_queue_head_index;
    uint32_t                   coded_frames_stat_queue_tail_index;

    uint64_t total_bit_actual_per_sw;
    uint64_t max_bit_actual_per_sw;
    uint64_t max_bit_actual_per_gop;
    uint64_t min_bit_actual_per_gop;
    uint64_t avg_bit_actual_per_gop;
    uint64_t rate_average_periodin_frames;

    EbHandle rc_mutex;
    // For dynamic resize, 1 pass cbr.
    RESIZE_STATE resize_state;
    int32_t      resize_avg_qp;
    int32_t      resize_buffer_underflow;
    int32_t      resize_count;
    int32_t      last_q[FRAME_TYPES]; // Q used on last encoded frame of the given type.

} RATE_CONTROL;

/**************************************
 * Input Port Types
 **************************************/
typedef enum RateControlInputPortTypes {
    RATE_CONTROL_INPUT_PORT_INLME          = 0,
    RATE_CONTROL_INPUT_PORT_PACKETIZATION  = 1,
    RATE_CONTROL_INPUT_PORT_ENTROPY_CODING = 2,
    RATE_CONTROL_INPUT_PORT_TOTAL_COUNT    = 3,
    RATE_CONTROL_INPUT_PORT_INVALID        = ~0,
} RateControlInputPortTypes;

/**************************************
 * Input Port Config
 **************************************/
typedef struct RateControlPorts {
    RateControlInputPortTypes type;
    uint32_t                  count;
} RateControlPorts;

typedef enum PicMgrInputPortTypes {
    PIC_MGR_INPUT_PORT_SOP           = 0,
    PIC_MGR_INPUT_PORT_PACKETIZATION = 1,
    PIC_MGR_INPUT_PORT_REST          = 2,
    PIC_MGR_INPUT_PORT_TOTAL_COUNT   = 3,
    PIC_MGR_INPUT_PORT_INVALID       = ~0,
} PicMgrInputPortTypes;
typedef struct PicMgrPorts {
    PicMgrInputPortTypes type;
    uint32_t             count;
} PicMgrPorts;
/**************************************
 * Context
 **************************************/

/**************************************
 * Extern Function Declarations
 **************************************/
int32_t svt_av1_convert_qindex_to_q_fp8(int32_t qindex, EbBitDepth bit_depth);
double  svt_av1_convert_qindex_to_q(int32_t qindex, EbBitDepth bit_depth);
int     svt_av1_rc_get_default_min_gf_interval(int width, int height, double framerate);
int     svt_av1_rc_get_default_max_gf_interval(double framerate, int min_gf_interval);
double  svt_av1_get_gfu_boost_projection_factor(double min_factor, double max_factor,
                                                int frame_count);

EbErrorType rate_control_context_ctor(EbThreadContext   *thread_context_ptr,
                                      const EbEncHandle *enc_handle_ptr, int me_port_index);

extern void *rate_control_kernel(void *input_ptr);
int svt_aom_compute_rd_mult_based_on_qindex(EbBitDepth bit_depth, FRAME_UPDATE_TYPE update_type,
                                            int qindex);
struct PictureControlSet;
int svt_aom_compute_rd_mult(struct PictureControlSet *pcs, uint8_t q_index, uint8_t me_q_index,
                            uint8_t bit_depth);
#endif // EbRateControl_h
