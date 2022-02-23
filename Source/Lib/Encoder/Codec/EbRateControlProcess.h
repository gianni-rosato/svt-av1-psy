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

#define MINQ_ADJ_LIMIT 48
#define MINQ_ADJ_LIMIT_CQ 20
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
// max bit rate average period in second. default is set to 2 second
#define MAX_RATE_AVG_PERIOD_IN_SEC 2
#define CODED_FRAMES_STAT_QUEUE_MAX_DEPTH 2000
#define CRITICAL_BUFFER_LEVEL 15
#define OPTIMAL_BUFFER_LEVEL 70
/**************************************
 * Coded Frames Stats
 **************************************/
typedef struct coded_frames_stats_entry {
    EbDctor  dctor;
    uint64_t picture_number;
    int64_t  frame_total_bit_actual;
    EbBool   end_of_sequence_flag;
} coded_frames_stats_entry;

extern EbErrorType rate_control_coded_frames_stats_context_ctor(coded_frames_stats_entry *entry_ptr,
                                                                uint64_t picture_number);
typedef struct {
    int    last_boosted_qindex; // Last boosted GF/KF/ARF q
    int    gfu_boost;
    int    kf_boost;
    double rate_correction_factors[MAX_TEMPORAL_LAYERS + 1];
    int    min_gf_interval;
    int    max_gf_interval;
    int    frames_till_gf_update_due;
    int    onepass_cbr_mode; // 0: not 1pass cbr, 1: 1pass cbr normal, 2: 1pass cbr real time
    int    baseline_gf_interval;
    int    constrained_gf_group;
    int    frames_to_key;
    int    frames_since_key;
    int    this_key_frame_forced;
    int    is_src_frame_alt_ref;

    int     avg_frame_bandwidth; // Average frame size target for clip
    int     min_frame_bandwidth; // Minimum allocation used for any frame
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
    int gf_interval;

    // gop bit budget
    int64_t gf_group_bits;

    // Total number of stats used only for kf_boost calculation.
    int num_stats_used_for_kf_boost;
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
} RATE_CONTROL;
typedef struct RateControlIntervalParamContext {
    EbDctor  dctor;
    uint64_t first_poc;
    uint64_t last_poc;

    // Projected total bits available for a key frame group of frames
    int64_t kf_group_bits;
    // Error score of frames still to be coded in kf group
    int64_t kf_group_error_left;
    uint8_t end_of_seq_seen;
    int32_t processed_frame_number;
    uint8_t last_i_qp;
} RateControlIntervalParamContext;
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
int32_t svt_av1_convert_qindex_to_q_fp8(int32_t qindex, AomBitDepth bit_depth);
double  svt_av1_convert_qindex_to_q(int32_t qindex, AomBitDepth bit_depth);
int     svt_av1_rc_get_default_min_gf_interval(int width, int height, double framerate);
int     svt_av1_rc_get_default_max_gf_interval(double framerate, int min_gf_interval);
double  svt_av1_get_gfu_boost_projection_factor(double min_factor, double max_factor,
                                                int frame_count);

EbErrorType rate_control_context_ctor(EbThreadContext   *thread_context_ptr,
                                      const EbEncHandle *enc_handle_ptr, int me_port_index);

extern void *rate_control_kernel(void *input_ptr);
#endif // EbRateControl_h
