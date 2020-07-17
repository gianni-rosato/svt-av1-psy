/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbAppConfig_h
#define EbAppConfig_h

#include <stdio.h>

#include "EbSvtAv1Enc.h"

#ifdef _WIN32
#define fseeko _fseeki64
#define ftello _ftelli64
#endif
// Define Cross-Platform 64-bit fseek() and ftell()

/** The AppExitConditionType type is used to define the App main loop exit
conditions.
*/
typedef enum AppExitConditionType {
    APP_ExitConditionNone = 0,
    APP_ExitConditionFinished,
    APP_ExitConditionError
} AppExitConditionType;

/** The AppPortActiveType type is used to define the state of output ports in
the App.
*/
typedef enum AppPortActiveType { APP_PortActive = 0, APP_PortInactive } AppPortActiveType;

/** The EbPtr type is intended to be used to pass pointers to and from the svt
API.  This is a 32 bit pointer and is aligned on a 32 bit word boundary.
*/
typedef void *EbPtr;

#define WARNING_LENGTH 100

// memory map to be removed and replaced by malloc / free
typedef enum EbPtrType {
    EB_N_PTR     = 0, // malloc'd pointer
    EB_A_PTR     = 1, // malloc'd pointer aligned
    EB_MUTEX     = 2, // mutex
    EB_SEMAPHORE = 3, // semaphore
    EB_THREAD    = 4 // thread handle
} EbPtrType;
typedef struct EbMemoryMapEntry {
    EbPtr     ptr; // points to a memory pointer
    EbPtrType ptr_type; // pointer type
} EbMemoryMapEntry;

extern EbMemoryMapEntry *app_memory_map; // App Memory table
extern uint32_t *        app_memory_map_index; // App Memory index
extern uint64_t *        total_app_memory; // App Memory malloc'd
extern uint32_t          app_malloc_count;

#define MAX_APP_NUM_PTR (0x186A0 << 2) // Maximum number of pointers to be allocated for the app

#define EB_APP_MALLOC(type, pointer, n_elements, pointer_class, return_type) \
    pointer = (type)malloc(n_elements);                                      \
    if (pointer == (type)NULL) {                                             \
        return return_type;                                                  \
    } else {                                                                 \
        app_memory_map[*(app_memory_map_index)].ptr_type = pointer_class;    \
        app_memory_map[(*(app_memory_map_index))++].ptr  = pointer;          \
        if (n_elements % 8 == 0) {                                           \
            *total_app_memory += (n_elements);                               \
        } else {                                                             \
            *total_app_memory += ((n_elements) + (8 - ((n_elements) % 8)));  \
        }                                                                    \
    }                                                                        \
    if (*(app_memory_map_index) >= MAX_APP_NUM_PTR) { return return_type; }  \
    app_malloc_count++;

#define EB_APP_MALLOC_NR(type, pointer, n_elements, pointer_class, return_type) \
    (void)return_type;                                                          \
    pointer = (type)malloc(n_elements);                                         \
    if (pointer == (type)NULL) {                                                \
        return_type = EB_ErrorInsufficientResources;                            \
        fprintf(stderr, "Malloc has failed due to insuffucient resources");     \
        return;                                                                 \
    } else {                                                                    \
        app_memory_map[*(app_memory_map_index)].ptr_type = pointer_class;       \
        app_memory_map[(*(app_memory_map_index))++].ptr  = pointer;             \
        if (n_elements % 8 == 0) {                                              \
            *total_app_memory += (n_elements);                                  \
        } else {                                                                \
            *total_app_memory += ((n_elements) + (8 - ((n_elements) % 8)));     \
        }                                                                       \
    }                                                                           \
    if (*(app_memory_map_index) >= MAX_APP_NUM_PTR) {                           \
        return_type = EB_ErrorInsufficientResources;                            \
        fprintf(stderr, "Malloc has failed due to insuffucient resources");     \
        return;                                                                 \
    }                                                                           \
    app_malloc_count++;

#define EB_APP_MEMORY()                                                        \
    fprintf(stderr, "Total Number of Mallocs in App: %d\n", app_malloc_count); \
    fprintf(stderr, "Total App Memory: %.2lf KB\n\n", *total_app_memory / (double)1024);

#define MAX_CHANNEL_NUMBER 6
#define MAX_NUM_TOKENS 210

#ifdef _WIN32
#define FOPEN(f, s, m) fopen_s(&f, s, m)
#else
#define FOPEN(f, s, m) f = fopen(s, m)
#endif

typedef struct EbPerformanceContext {
    /****************************************
     * Computational Performance Data
     ****************************************/
    uint64_t lib_start_time[2]; // [sec, micro_sec] including init time
    uint64_t encode_start_time[2]; // [sec, micro_sec] first frame sent

    double total_execution_time; // includes init
    double total_encode_time; // not including init

    uint64_t total_latency;
    uint32_t max_latency;

    uint64_t starts_time;
    uint64_t startu_time;
    uint64_t frame_count;

    double average_speed;
    double average_latency;

    uint64_t byte_count;
    double   sum_luma_psnr;
    double   sum_cr_psnr;
    double   sum_cb_psnr;

    double sum_luma_sse;
    double sum_cr_sse;
    double sum_cb_sse;

    double sum_luma_ssim;
    double sum_cr_ssim;
    double sum_cb_ssim;

    uint64_t sum_qp;

} EbPerformanceContext;

typedef struct EbConfig {
    /****************************************
     * File I/O
     ****************************************/
    FILE *        config_file;
    FILE *        input_file;
    EbBool        input_file_is_fifo;
    FILE *        bitstream_file;
    FILE *        recon_file;
    FILE *        error_log_file;
    FILE *        stat_file;
    FILE *        buffer_file;
    FILE *        qp_file;
    FILE *        input_stat_file;
    FILE *        output_stat_file;
    FILE *        input_pred_struct_file;
    char *        input_pred_struct_filename;
    EbBool        use_input_stat_file;
    EbBool        use_output_stat_file;
    EbBool        y4m_input;
    unsigned char y4m_buf[9];
    EbBool        use_qp_file;
    EbBool        no_progress;
    uint8_t       stat_report;
    uint32_t      frame_rate;
    uint32_t      frame_rate_numerator;
    uint32_t      frame_rate_denominator;
    uint32_t      injector_frame_rate;
    uint32_t      injector;
    uint32_t      speed_control_flag;
    uint32_t      encoder_bit_depth;
    EbBool        is_16bit_pipeline;
    uint32_t      encoder_color_format;
    uint32_t      compressed_ten_bit_format;
    uint32_t      source_width;
    uint32_t      source_height;

    uint32_t input_padded_width;
    uint32_t input_padded_height;

    // -1 indicates unknown (auto-detect at earliest opportunity)
    // auto-detect is performed on load for files and at end of stream for pipes
    int64_t   frames_to_be_encoded;
    int32_t   frames_encoded;
    int32_t   buffered_input;
    uint8_t **sequence_buffer;

    uint8_t latency_mode;

    /*****************************************
     * Coding Structure
     *****************************************/
    uint8_t  enc_mode;
    uint8_t  snd_pass_enc_mode;
    int32_t  intra_period;
    uint32_t intra_refresh_type;
    uint32_t hierarchical_levels;
    uint32_t pred_structure;

    /****************************************
     * Quantization
     ****************************************/
    uint32_t qp;

    /****************************************
     * Film Grain
     ****************************************/
    uint32_t film_grain_denoise_strength;
    /****************************************
     * DLF
     ****************************************/
    EbBool disable_dlf_flag;

    /****************************************
     * Local Warped Motion
     ****************************************/
    int enable_warped_motion;

    /****************************************
     * Global Motion
     ****************************************/
    EbBool enable_global_motion;

    /****************************************
     * CDEF Mode
     * 0         OFF
     * 1         1 step refinement
     * 2         4 step refinement
     * 3         8 step refinement
     * 4         16 step refinement
     * 5         64 step refinement
    ****************************************/
    int cdef_mode;

    /****************************************
     * Restoration filtering
    ****************************************/
    int enable_restoration_filtering;
    int sg_filter_mode;
    int wn_filter_mode;

#if 0//!REMOVE_COMBINE_CLASS12

    /****************************************
     * class12
    ****************************************/
    int combine_class_12;
#endif
    /****************************************
     * edge based skip angle intra
    ****************************************/
    int edge_skp_angle_intra;
    /****************************************
     * intra angle delta
    ****************************************/
    int intra_angle_delta;
    /****************************************
     * intra inter compoound
    ****************************************/
    int inter_intra_compound;
    /****************************************
     * paeth
    ****************************************/
    int enable_paeth;
    /****************************************
     * smooth
    ****************************************/
    int enable_smooth;
    /****************************************
     * motion field motion vector
    ****************************************/
    int enable_mfmv;
    /****************************************
     * redundant block
    ****************************************/
    int enable_redundant_blk;
    /****************************************
      * spatial sse in full loop
     ****************************************/
    int spatial_sse_fl;
#if 1//!REMOVE_ME_SUBPEL_CODE
    /****************************************
      * subpel
     ****************************************/
    int enable_subpel;
#endif
    /****************************************
      * over boundry block
     ****************************************/
    int over_bndry_blk;
    /****************************************
      * new nearest comb injection
     ****************************************/
    int new_nearest_comb_inject;
    /****************************************
      * prune unipred at me
     ****************************************/
    int prune_unipred_me;
    /****************************************
      * prune ref frame for rec partitions
     ****************************************/
    int prune_ref_rec_part;
    /****************************************
      * nsq table
     ****************************************/
    int nsq_table;
    /****************************************
      * frame end cdf update
     ****************************************/
    int frame_end_cdf_update;
    /****************************************
      * predictive me
     ****************************************/
    int pred_me;
    /****************************************
      * bipred 3x3 injection
     ****************************************/
    int bipred_3x3_inject;
    /****************************************
      * compound level
     ****************************************/
    int compound_level;

    /****************************************
     * Chroma
     *
     * Level                Settings
     * CHROMA_MODE_0  0     Full chroma search @ MD
     * CHROMA_MODE_1  1     Fast chroma search @ MD
     * CHROMA_MODE_2  2     Chroma blind @ MD + CFL @ EP
     * CHROMA_MODE_3  3     Chroma blind @ MD + no CFL @ EP
     *
     * Default is -1 (AUTO)  */
    int set_chroma_mode;

    /* Disable chroma from luma (CFL)
     *
     * Default is -1 (auto) */
    int disable_cfl_flag;

    /****************************************
     * OBMC
     ****************************************/
    EbBool enable_obmc;

    /****************************************
     * RDOQ
     * ****************************************/
    int enable_rdoq;

    /****************************************
     * Filter intra prediction
     ****************************************/
    EbBool enable_filter_intra;

    /****************************************
     * Intra Edge Filter
     ****************************************/
    int enable_intra_edge_filter;

    /****************************************
     * Picture based rate estimation
     ****************************************/
    int pic_based_rate_est;

    /****************************************
     * ME Tools
     ****************************************/
    EbBool use_default_me_hme;
    EbBool enable_hme_flag;
    EbBool enable_hme_level0_flag;
    EbBool enable_hme_level1_flag;
    EbBool enable_hme_level2_flag;
    EbBool ext_block_flag;

    /****************************************
     * ME Parameters
     ****************************************/
    uint32_t search_area_width;
    uint32_t search_area_height;

    /****************************************
     * HME Parameters
     ****************************************/
    uint32_t number_hme_search_region_in_width;
    uint32_t number_hme_search_region_in_height;
    uint32_t hme_level0_total_search_area_width;
    uint32_t hme_level0_total_search_area_height;
    uint32_t hme_level0_column_index;
    uint32_t hme_level0_row_index;
    uint32_t hme_level1_column_index;
    uint32_t hme_level1_row_index;
    uint32_t hme_level2_column_index;
    uint32_t hme_level2_row_index;
    uint32_t hme_level0_search_area_in_width_array[EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT];
    uint32_t hme_level0_search_area_in_height_array[EB_HME_SEARCH_AREA_ROW_MAX_COUNT];
    uint32_t hme_level1_search_area_in_width_array[EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT];
    uint32_t hme_level1_search_area_in_height_array[EB_HME_SEARCH_AREA_ROW_MAX_COUNT];
    uint32_t hme_level2_search_area_in_width_array[EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT];
    uint32_t hme_level2_search_area_in_height_array[EB_HME_SEARCH_AREA_ROW_MAX_COUNT];

    /****************************************
     * MD Parameters
     ****************************************/
    int8_t  enable_hbd_mode_decision;
    int32_t enable_palette;
    int32_t tile_columns;
    int32_t tile_rows;

    /****************************************
     * Rate Control
     ****************************************/
    uint32_t scene_change_detection;
    uint32_t rate_control_mode;
    uint32_t look_ahead_distance;
#if 1//TPL_LA
    uint32_t enable_tpl_la;
#endif
    uint32_t target_bit_rate;
    uint32_t max_qp_allowed;
    uint32_t min_qp_allowed;
    uint32_t vbv_bufsize;

    EbBool enable_adaptive_quantization;

    /****************************************
     * Optional Features
     ****************************************/

    uint32_t screen_content_mode;
    int      intrabc_mode;
    uint32_t high_dynamic_range_input;
    EbBool   unrestricted_motion_vector;

    /****************************************
     * Annex A Parameters
     ****************************************/
    uint32_t profile;
    uint32_t tier;
    uint32_t level;

    /****************************************
     * On-the-fly Testing
     ****************************************/
    EbBool eos_flag;

    /****************************************
    * CPU FLAGS available
    ****************************************/
    CPU_FLAGS cpu_flags_limit;

    /****************************************
     * Computational Performance Data
     ****************************************/
    EbPerformanceContext performance_context;

    /****************************************
    * Instance Info
    ****************************************/
    uint32_t channel_id;
    uint32_t active_channel_count;
    uint32_t logical_processors;
    uint32_t unpin;
    int32_t  target_socket;
    EbBool   stop_encoder; // to signal CTRL+C Event, need to stop encoding.

    uint64_t processed_frame_count;
    uint64_t processed_byte_count;

    uint64_t byte_count_since_ivf;
    uint64_t ivf_count;

    // --- start: ALTREF_FILTERING_SUPPORT
    /****************************************
     * ALT-REF related Parameters
     ****************************************/
    EbBool  enable_altrefs;
    uint8_t altref_strength;
    uint8_t altref_nframes;
    EbBool  enable_overlays;
    // --- end: ALTREF_FILTERING_SUPPORT

    /****************************************
     * Super-resolution related Parameters
     ****************************************/
    SUPERRES_MODE superres_mode;
    uint8_t       superres_denom;
    uint8_t       superres_kf_denom;
    uint8_t       superres_qthres;

    // square cost weighting for deciding if a/b shapes could be skipped
    uint32_t sq_weight;

    // inter/intra class pruning costs before MD stage 1/2
    uint64_t md_stage_1_class_prune_th;
    uint64_t md_stage_1_cand_prune_th;
    uint64_t md_stage_2_3_class_prune_th;
    uint64_t md_stage_2_3_cand_prune_th;

    // prediction structure
    PredictionStructureConfigEntry pred_struct[1 << (MAX_HIERARCHICAL_LEVEL - 1)];
    EbBool                         enable_manual_pred_struct;
    int32_t                        manual_pred_struct_entry_num;
#if 1//ON_OFF_FEATURE_MRP
    int                 mrp_level;
#endif
} EbConfig;

extern void eb_config_ctor(EbConfig *config_ptr);
extern void eb_config_dtor(EbConfig *config_ptr);

extern EbErrorType read_command_line(int32_t argc, char *const argv[], EbConfig **config,
                                     uint32_t num_channels, EbErrorType *return_errors,
                                     char *warning_str[WARNING_LENGTH]);
extern uint32_t    get_help(int32_t argc, char *const argv[]);
extern uint32_t    get_number_of_channels(int32_t argc, char *const argv[]);

#endif //EbAppConfig_h
