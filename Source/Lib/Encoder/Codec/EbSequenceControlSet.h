/*
* Copyright(c) 2019 Intel Corporation
* Copyright (c) 2016, Alliance for Open Media. All rights reserved
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

#ifndef EbSequenceControlSet_h
#define EbSequenceControlSet_h

#include "EbDefinitions.h"
#include "EbAv1Structs.h"
#include "EbEncodeContext.h"
#include "EbObject.h"
#include "firstpass.h"

#ifdef __cplusplus
extern "C" {
#endif
#if CLIP_BASED_DYNAMIC_MINIGOP
    typedef struct MiniGopSizeCtrls {
        EbBool adptive_enable; // 0: Off, 1: Switch at clip level, 2: Switch at GOP level
        double short_shot_th; // Threshold to determine short scene.
        double animation_type_th; // Threshold to determine animation scene
        double lm_th; // Threshold to determine low motion scene
        double hm_th; // Threshold to determine high motion scene
        double lfr_th; // Threshold to determine low frame rate scene
        double hsa_th; // Threshold to determine high static area scene
        double hmv_di_th; // Threshold to determine high mv direction scene
        double lmv_di_th; // Threshold to determine low mv direction scene
    } MiniGopSizeCtrls;
#endif
/************************************
     * Sequence Control Set
     ************************************/
typedef struct SequenceControlSet {
    /*!< Pointer to the dtor of the struct*/
    EbDctor dctor;
    /*!< Encoding context pointer containing the handle pointer */
    EncodeContext *encode_context_ptr;
     /*!< 2ndpass enc mode, available at firstpass encoder */
     int8_t enc_mode_2ndpass;
    /*!< API structure */
    EbSvtAv1EncConfiguration static_config;
    /*!< Pointer to prediction structure containing the mini-gop information */
    PredictionStructure *pred_struct_ptr;
    /*!< Super block geomerty pointer */
    SbGeom *sb_geom;
    /*!< Array of superblock parameters computed at the resource coordination stage */
    SbParams *sb_params_array;
    /*!< Bitstream level */
    BitstreamLevel level[MAX_NUM_OPERATING_POINTS];
    /*!< Sequence header structure, common between the encoder and decoder */
    SeqHeader seq_header;

    /*!< Sequence coding parameters
            parameters/features are set to be set for the full stream
            but encoding decisions may still be taken at a picture / sub picture level
    */

    /*!< Maximum number of references that a picture can have within the stream needs to be cleaned up*/
    uint32_t max_ref_count;
    /*!< Maximum number of references that a picture can have within the stream */
    uint32_t reference_count;
    /*!< Maximum number of allowed temporal layers */
    uint32_t max_temporal_layers;
    /*!< Overflow bits used for the picture order count increments */
    uint32_t bits_for_picture_order_count;
    /*!< Screen change detection mode 0=OFF, 1=use decimated picture, 2=use full picture */
    EbScdMode scd_mode;
    /*!< Number of delay frames needed to implement future window
         for algorithms such as SceneChange or TemporalFiltering */
    uint32_t scd_delay;
    /*!<  */
    EbBlockMeanPrec block_mean_calc_prec;
    /*!< CDF (The signal changes per preset; 0: CDF update, 1: no CDF update) Default is 0.*/
    uint8_t cdf_mode;
    /*!< Down-sampling method @ ME and alt-ref temporal filtering
        (The signal changes per preset; 0: filtering, 1: decimation) Default is 0. */
    uint8_t down_sampling_method_me_search;
#if CLN_GEOM
    uint32_t geom_idx;   //geometry type
#endif

    /*  1..15    | 17..31  | 33..47  |
              16 |       32|       48|
      lad mg=2: delay the first MG (1-16) until the next 2 MGs(17-48) are gop , TF, and ME ready
    */
    uint8_t lad_mg;   //delay all pictures within a given MG, until N future MGs are  gop , TF, and ME ready
#if FTR_LAD_INPUT
    uint8_t tpl_lad_mg;   //delay all pictures within a given MG, until N future MGs are  gop , TF, and ME ready used for tpl
#endif
    /*!< 1: Specifies that loop restoration filter should use boundary pixels in the search.  Must be
            set at the sequence level because it requires a buffer allocation to copy the pixels
            to be used in the search.
         0: Specifies that loop restoration filter should not use boundary pixels in the search.*/
    uint8_t use_boundaries_in_rest_search;
    uint8_t enable_pic_mgr_dec_order; // if enabled: pic mgr starts pictures in dec order
    uint8_t enable_dec_order; // if enabled: encoding are in dec order
    /*!< Use in loop motion OIS
         Default is 1. */
    uint8_t in_loop_ois;
    /*!< Allow the usage of motion field motion vector in the stream
        (The signal changes per preset; 0: Enabled, 1: Disabled) Default is 1. */
    uint8_t mfmv_enabled;
    /*!< Film grain seed */
    uint16_t film_grain_random_seed;
    /*!< over_boundary_block: pad resolution to a multiple of SB for smaller overhead
        (The signal changes per preset; 0: No over boundary blk allowed, 1: over boundary blk allowed) Default is 1.
        to enable when md_skip_blk is on */
    uint8_t over_boundary_block_mode;
    /*!< Enable compound prediction to be used in the stream, decisions will be taken at a picture level subsequently
    (The signal changes per preset; 0: compound disabled, 1: compound enabled) Default is 1. */
    uint8_t compound_mode;

    /*!< Sequence resolution parameters */
    uint32_t          chroma_format_idc;
    uint16_t          subsampling_x; // add chroma subsampling parameters
    uint16_t          subsampling_y;
    uint16_t          max_input_luma_width;
    uint16_t          max_input_luma_height;
    uint16_t          max_input_chroma_width;
    uint16_t          max_input_chroma_height;
    uint16_t          max_input_pad_bottom;
    uint16_t          max_input_pad_right;
    uint32_t          chroma_width;
    uint32_t          chroma_height;
    uint32_t          pad_right;
    uint32_t          pad_bottom;
    uint16_t          left_padding;
    uint16_t          top_padding;
    uint16_t          right_padding;
    uint16_t          bot_padding;
    uint32_t          frame_rate;
    uint32_t          encoder_bit_depth;
    EbInputResolution input_resolution;

    /*!< Super block parameters set for the stream */
    uint8_t  sb_sz;
    uint8_t  max_sb_depth;
#if FTR_16K
    uint16_t  pic_width_in_sb;
    uint16_t  picture_height_in_sb;
#else
    uint8_t  pic_width_in_sb;
    uint8_t  picture_height_in_sb;
#endif
    uint16_t sb_total_count;
    uint16_t sb_size_pix;
    uint16_t sb_tot_cnt;
    uint16_t max_block_cnt;
#if FTR_NEW_WN_LVLS
    /*!< Restoration Unit parameters set for the stream */
    int32_t rest_units_per_tile;
#endif
    /*!< Block limits */
    uint8_t max_blk_size;
    uint8_t min_blk_size;
    uint8_t max_intra_size;
    uint8_t min_intra_size;

    /*!< Sub picture reagions for picture analysis */
    uint32_t picture_analysis_number_of_regions_per_width;
    uint32_t picture_analysis_number_of_regions_per_height;

    /*!< Tile groups per hierarchical layers */
    uint8_t tile_group_col_count_array[MAX_TEMPORAL_LAYERS];
    uint8_t tile_group_row_count_array[MAX_TEMPORAL_LAYERS];

    /*!< Segements (sub picture) count for different processes */
    uint32_t me_segment_column_count_array[MAX_TEMPORAL_LAYERS];
    uint32_t me_segment_row_count_array[MAX_TEMPORAL_LAYERS];
#if OPT_1P
    uint32_t fpass_segment_column_count;
    uint32_t fpass_segment_row_count;
#endif
    uint32_t enc_dec_segment_col_count_array[MAX_TEMPORAL_LAYERS];
    uint32_t enc_dec_segment_row_count_array[MAX_TEMPORAL_LAYERS];
    uint32_t tpl_segment_col_count_array;
    uint32_t tpl_segment_row_count_array;
    uint32_t cdef_segment_column_count;
    uint32_t cdef_segment_row_count;
    uint32_t rest_segment_column_count;
    uint32_t rest_segment_row_count;
    uint32_t tf_segment_column_count;
    uint32_t tf_segment_row_count;

    /*!< Picture, reference, recon and input output buffer count */
    uint32_t picture_control_set_pool_init_count;
    uint32_t me_pool_init_count;
    uint32_t picture_control_set_pool_init_count_child;
    uint32_t enc_dec_pool_init_count;
    uint32_t pa_reference_picture_buffer_init_count;
    uint32_t reference_picture_buffer_init_count;
    uint32_t input_buffer_fifo_init_count;
    uint32_t overlay_input_picture_buffer_init_count;
    uint32_t output_stream_buffer_fifo_init_count;
    uint32_t output_recon_buffer_fifo_init_count;

    /*!< Inter processes fifos count */
    uint32_t resource_coordination_fifo_init_count;
    uint32_t picture_analysis_fifo_init_count;
    uint32_t picture_decision_fifo_init_count;
    uint32_t motion_estimation_fifo_init_count;
    uint32_t initial_rate_control_fifo_init_count;
    uint32_t picture_demux_fifo_init_count;
    uint32_t tpl_disp_fifo_init_count;
    uint32_t rate_control_tasks_fifo_init_count;
    uint32_t rate_control_fifo_init_count;
    uint32_t mode_decision_configuration_fifo_init_count;
    uint32_t enc_dec_fifo_init_count;
    uint32_t entropy_coding_fifo_init_count;
    uint32_t dlf_fifo_init_count;
    uint32_t cdef_fifo_init_count;
    uint32_t rest_fifo_init_count;

    /*!< Thread count for each process */
    uint32_t picture_analysis_process_init_count;
    uint32_t motion_estimation_process_init_count;
    uint32_t source_based_operations_process_init_count;
    uint32_t mode_decision_configuration_process_init_count;
    uint32_t enc_dec_process_init_count;
    uint32_t entropy_coding_process_init_count;
    uint32_t dlf_process_init_count;
    uint32_t cdef_process_init_count;
    uint32_t rest_process_init_count;
    uint32_t tpl_disp_process_init_count;
    uint32_t total_process_init_count;
    int32_t  lap_enabled;
    TWO_PASS twopass;
    double   double_frame_rate;
    Quants   quants_bd; // follows input bit depth
    Dequants deq_bd; // follows input bit depth
    Quants   quants_8bit; // 8bit
    Dequants deq_8bit; // 8bit
    ScaleFactors sf_identity;
    uint8_t  mrp_init_level; //sequence based MRP level
    int32_t nmv_vec_cost[MV_JOINTS];
    int32_t nmv_costs[2][MV_VALS];
    uint8_t mvrate_set;
#if CLIP_BASED_DYNAMIC_MINIGOP
    MiniGopSizeCtrls mgs_ctls;
#endif
#if FTR_OPT_MPASS
    /*!< The RC stat generation pass mode (0: The default, 1: optimized)*/
    uint8_t rc_stat_gen_pass_mode;
#endif
#if FTR_NEW_QPS
    int cqp_base_q_tf;
    int cqp_base_q;
#endif
} SequenceControlSet;

typedef struct EbSequenceControlSetInitData {
    EncodeContext *encode_context_ptr;
    int32_t        sb_size;
} EbSequenceControlSetInitData;

typedef struct EbSequenceControlSetInstance {
    EbDctor             dctor;
    EncodeContext *     encode_context_ptr;
    SequenceControlSet *scs_ptr;
    EbHandle            config_mutex;
} EbSequenceControlSetInstance;

/**************************************
     * Extern Function Declarations
     **************************************/
extern EbErrorType svt_sequence_control_set_creator(EbPtr *object_dbl_ptr,
                                                    EbPtr  object_init_data_ptr);

extern EbErrorType svt_sequence_control_set_ctor(SequenceControlSet *object,
                                                 EbPtr               object_init_data_ptr);

extern EbErrorType copy_sequence_control_set(SequenceControlSet *dst, SequenceControlSet *src);

extern EbErrorType svt_sequence_control_set_instance_ctor(EbSequenceControlSetInstance *object_ptr);

extern EbErrorType sb_params_init(SequenceControlSet *scs_ptr);

extern EbErrorType derive_input_resolution(EbInputResolution *input_resolution,
                                           uint32_t           input_size);

EbErrorType sb_geom_init(SequenceControlSet *scs_ptr);

inline static EbBool use_input_stat(const SequenceControlSet *scs_ptr) {
    return !!scs_ptr->static_config.rc_twopass_stats_in.sz;
}

inline static EbBool use_output_stat(const SequenceControlSet *scs_ptr) {
    return scs_ptr->static_config.rc_firstpass_stats_out;
}
#if FTR_MULTI_PASS_API
inline static EbBool is_middle_pass(const SequenceControlSet *scs_ptr) {
    return scs_ptr->static_config.rc_middlepass_stats_out;
}
#if FTR_OPT_MPASS_DOWN_SAMPLE
inline static EbBool is_middle_pass_ds(const SequenceControlSet *scs_ptr) {
    return scs_ptr->static_config.rc_middlepass_ds_stats_out;
}
#endif
#endif
#ifdef __cplusplus
}
#endif
#endif // EbSequenceControlSet_h
