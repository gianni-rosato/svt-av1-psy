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

#ifndef EbEncodeContext_h
#define EbEncodeContext_h

#include <stdio.h>

#include "EbDefinitions.h"
#include "EbSvtAv1Enc.h"
#include "EbPictureDecisionReorderQueue.h"
#include "EbPictureDecisionQueue.h"
#include "EbPictureManagerQueue.h"
#include "EbPacketizationReorderQueue.h"
#include "EbInitialRateControlReorderQueue.h"
#include "EbPictureManagerReorderQueue.h"
#include "EbCabacContextModel.h"
#include "EbMdRateEstimation.h"
#include "EbPredictionStructure.h"
#include "EbRateControlTables.h"
#include "EbObject.h"
#include "encoder.h"
#include "firstpass.h"

// *Note - the queues are small for testing purposes.  They should be increased when they are done.
#define PRE_ASSIGNMENT_MAX_DEPTH 128 // should be large enough to hold an entire prediction period
#define INPUT_QUEUE_MAX_DEPTH 5000
#define REFERENCE_QUEUE_MAX_DEPTH 5000
#define PICTURE_DECISION_PA_REFERENCE_QUEUE_MAX_DEPTH 5000

#define PICTURE_DECISION_REORDER_QUEUE_MAX_DEPTH 2048
#define INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH 2048
#define PICTURE_MANAGER_REORDER_QUEUE_MAX_DEPTH 2048
#define HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH 2048
#define PACKETIZATION_REORDER_QUEUE_MAX_DEPTH 2048

// RC Groups: They should be a power of 2, so we can replace % by &.
// Instead of using x % y, we use x && (y-1)
#define PARALLEL_GOP_MAX_NUMBER 256
#define RC_GROUP_IN_GOP_MAX_NUMBER 512
#define PICTURE_IN_RC_GROUP_MAX_NUMBER 64

typedef struct DpbDependentList
{
    int32_t                              list[1 << MAX_TEMPORAL_LAYERS];
    uint32_t                             list_count;
} DpbDependentList;

typedef struct DPBInfo {
    uint64_t picture_number;
    int32_t dep_count;
    int32_t dep_list0_count;
    int32_t dep_list1_count;
    uint8_t temporal_layer_index;
    EbBool  is_displayed;
    EbBool  is_used;
    EbBool  is_alt_ref;
    DpbDependentList dep_list0;
    DpbDependentList dep_list1;
} DPBInfo;

typedef struct FirstPassStatsOut {
    FIRSTPASS_STATS* stat;
    size_t size;
    size_t capability;
} FirstPassStatsOut;

typedef struct EncodeContext {
    EbDctor dctor;
    // Callback Functions
    EbCallback *app_callback_ptr;

    EbHandle total_number_of_recon_frame_mutex;
    uint64_t total_number_of_recon_frames;

    // Overlay input picture fifo
    EbFifo *overlay_input_picture_pool_fifo_ptr;
    // Output Buffer Fifos
    EbFifo *stream_output_fifo_ptr;
    EbFifo *recon_output_fifo_ptr;

    // Picture Buffer Fifos
    EbFifo *reference_picture_pool_fifo_ptr;
    EbFifo *pa_reference_picture_pool_fifo_ptr;
#if FEATURE_INL_ME
    EbFifo *down_scaled_picture_pool_fifo_ptr;
#endif

    // Picture Decision Reorder Queue
    PictureDecisionReorderEntry **picture_decision_reorder_queue;
    uint32_t                      picture_decision_reorder_queue_head_index;
    //hold undisplayed frame for show existing frame. It's ordered with pts Descend.
    EbObjectWrapper              *picture_decision_undisplayed_queue[REF_FRAMES];
    uint32_t                      picture_decision_undisplayed_queue_count;
    // Picture Manager Pre-Assignment Buffer
    uint32_t          pre_assignment_buffer_intra_count;
    uint32_t          pre_assignment_buffer_idr_count;
    uint32_t          pre_assignment_buffer_scene_change_count;
    uint32_t          pre_assignment_buffer_scene_change_index;
    uint32_t          pre_assignment_buffer_eos_flag;
    uint64_t          decode_base_number;
    EbObjectWrapper **pre_assignment_buffer;
    uint32_t          pre_assignment_buffer_count;

    // Picture Decision Circular Queues
    PaReferenceQueueEntry **picture_decision_pa_reference_queue;
    uint32_t                picture_decision_pa_reference_queue_head_index;
    uint32_t                picture_decision_pa_reference_queue_tail_index;

    // Picture Manager Circular Queues
    InputQueueEntry **    input_picture_queue;
    uint32_t              input_picture_queue_head_index;
    uint32_t              input_picture_queue_tail_index;
    ReferenceQueueEntry **reference_picture_queue;
    uint32_t              reference_picture_queue_head_index;
    uint32_t              reference_picture_queue_tail_index;

    // Initial Rate Control Reorder Queue
    InitialRateControlReorderEntry **initial_rate_control_reorder_queue;
    uint32_t                         initial_rate_control_reorder_queue_head_index;
    uint32_t        dep_q_head;
    uint32_t        dep_q_tail;
    PicQueueEntry **dep_cnt_picture_queue; //buffer to sotre all pictures needing dependent-count clean-up in PicMgr

    // High Level Rate Control Histogram Queue
    HlRateControlHistogramEntry **hl_rate_control_historgram_queue;
    uint32_t                      hl_rate_control_historgram_queue_head_index;
    EbHandle                      hl_rate_control_historgram_queue_mutex;

    // Packetization Reorder Queue
    PacketizationReorderEntry **packetization_reorder_queue;
    uint32_t                    packetization_reorder_queue_head_index;

    // GOP Counters
    uint32_t intra_period_position; // Current position in intra period
    uint32_t pred_struct_position; // Current position within a prediction structure
    uint32_t elapsed_non_idr_count;
    uint32_t elapsed_non_cra_count;
    int64_t  current_input_poc;
    EbBool   initial_picture;
    uint64_t last_idr_picture; // the most recently occured IDR picture (in decode order)

    // Sequence Termination Flags
    uint64_t terminating_picture_number;
    EbBool   terminating_sequence_flag_received;

    // Signalling the need for a td structure to be written in the Bitstream - only used in the PK process so no need for a mutex
    EbBool td_needed;

    // Prediction Structure
    PredictionStructureGroup *prediction_structure_group_ptr;
    // Rate Control Bit Tables
    RateControlTables *rate_control_tables_array;
    EbBool             rate_control_tables_array_updated;
    EbHandle           rate_table_update_mutex;

    // Speed Control
    int64_t   sc_buffer;
    int64_t   sc_frame_in;
    int64_t   sc_frame_out;
    EbHandle  sc_buffer_mutex;
    EbEncMode enc_mode;

    // Rate Control
    uint32_t previous_selected_ref_qp;
    uint64_t max_coded_poc;
    uint32_t max_coded_poc_selected_ref_qp;

    // Dynamic GOP
    uint32_t         previous_mini_gop_hierarchical_levels;
    EbObjectWrapper *previous_picture_control_set_wrapper_ptr;
    EbHandle         shared_reference_mutex;
    uint64_t picture_number_alt; // The picture number overlay includes all the overlay frames

    EbHandle stat_file_mutex;

    //DPB list management
    DPBInfo dpb_list[REF_FRAMES];
    uint64_t display_picture_number;
    EbBool  is_mini_gop_changed;
    EbBool  is_i_slice_in_last_mini_gop;
    uint64_t i_slice_picture_number_in_last_mini_gop;
    uint64_t poc_map_idx[MAX_TPL_LA_SW];
#if FEATURE_IN_LOOP_TPL
    EbPictureBufferDesc *mc_flow_rec_picture_buffer[MAX_TPL_LA_SW];
    EbPictureBufferDesc *mc_flow_rec_picture_buffer_noref;
#else
    EbByte  mc_flow_rec_picture_buffer[MAX_TPL_LA_SW];
    EbByte  mc_flow_rec_picture_buffer_saved;
#endif
    FrameInfo      frame_info;
    TwoPassCfg     two_pass_cfg; // two pass datarate control
    RATE_CONTROL   rc;
    RateControlCfg rc_cfg;
    GF_GROUP       gf_group;
    KeyFrameCfg    kf_cfg;
    GFConfig       gf_cfg;
    FIRSTPASS_STATS *frame_stats_buffer;
    // Number of stats buffers required for look ahead
    int num_lap_buffers;
    STATS_BUFFER_CTX stats_buf_context;
    SvtAv1FixedBuf rc_twopass_stats_in; // replaced oxcf->two_pass_cfg.stats_in in aom
    FirstPassStatsOut stats_out;
} EncodeContext;

typedef struct EncodeContextInitData {
    int32_t junk;
} EncodeContextInitData;

/**************************************
 * Extern Function Declarations
 **************************************/
extern EbErrorType encode_context_ctor(EncodeContext *encode_context_ptr,
                                       EbPtr          object_init_data_ptr);
#endif // EbEncodeContext_h
