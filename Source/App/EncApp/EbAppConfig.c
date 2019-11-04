/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "EbAppString.h"
#include "EbAppConfig.h"
#include "EbAppInputy4m.h"

#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

/**********************************
 * Defines
 **********************************/
#define HELP_TOKEN                      "-help"
#define CHANNEL_NUMBER_TOKEN            "-nch"
#define COMMAND_LINE_MAX_SIZE           2048
#define CONFIG_FILE_TOKEN               "-c"
#define INPUT_FILE_TOKEN                "-i"
#define OUTPUT_BITSTREAM_TOKEN          "-b"
#define OUTPUT_RECON_TOKEN              "-o"
#define ERROR_FILE_TOKEN                "-errlog"
#define QP_FILE_TOKEN                   "-qp-file"
#if TWO_PASS
#define INPUT_STAT_FILE_TOKEN           "-input-stat-file"
#define OUTPUT_STAT_FILE_TOKEN          "-output-stat-file"
#endif
#define STAT_FILE_TOKEN                 "-stat-file"
#define WIDTH_TOKEN                     "-w"
#define HEIGHT_TOKEN                    "-h"
#define NUMBER_OF_PICTURES_TOKEN        "-n"
#define BUFFERED_INPUT_TOKEN            "-nb"
#define BASE_LAYER_SWITCH_MODE_TOKEN    "-base-layer-switch-mode" // no Eval
#define QP_TOKEN                        "-q"
#define USE_QP_FILE_TOKEN               "-use-q-file"
#define STAT_REPORT_TOKEN               "-stat-report"
#define FRAME_RATE_TOKEN                "-fps"
#define FRAME_RATE_NUMERATOR_TOKEN      "-fps-num"
#define FRAME_RATE_DENOMINATOR_TOKEN    "-fps-denom"
#define ENCODER_BIT_DEPTH               "-bit-depth"
#define ENCODER_COLOR_FORMAT            "-color-format"
#define INPUT_COMPRESSED_TEN_BIT_FORMAT "-compressed-ten-bit-format"
#define ENCMODE_TOKEN                   "-enc-mode"
#if TWO_PASS_USE_2NDP_ME_IN_1STP
#define ENCMODE2P_TOKEN                 "-enc-mode-2p"
#endif
#define HIERARCHICAL_LEVELS_TOKEN       "-hierarchical-levels" // no Eval
#define PRED_STRUCT_TOKEN               "-pred-struct"
#define INTRA_PERIOD_TOKEN              "-intra-period"
#define PROFILE_TOKEN                   "-profile"
#define TIER_TOKEN                      "-tier"
#define LEVEL_TOKEN                     "-level"
#define LATENCY_MODE                    "-latency-mode" // no Eval
#define FILM_GRAIN_TOKEN                "-film-grain"
#define INTERLACED_VIDEO_TOKEN          "-interlaced-video"
#define SEPERATE_FILDS_TOKEN            "-separate-fields"
#define INTRA_REFRESH_TYPE_TOKEN        "-irefresh-type" // no Eval
#define LOOP_FILTER_DISABLE_TOKEN       "-dlf"
#define LOCAL_WARPED_ENABLE_TOKEN       "-local-warp"
#define GLOBAL_MOTION_ENABLE_TOKEN      "-global-motion"
#define OBMC_TOKEN                      "-obmc"
#define RDOQ_TOKEN                      "-rdoq"
#define FILTER_INTRA_TOKEN              "-filter-intra"
#define USE_DEFAULT_ME_HME_TOKEN        "-use-default-me-hme"
#define HME_ENABLE_TOKEN                "-hme"
#define HME_L0_ENABLE_TOKEN             "-hme-l0"
#define HME_L1_ENABLE_TOKEN             "-hme-l1"
#define HME_L2_ENABLE_TOKEN             "-hme-l2"
#define EXT_BLOCK                       "-ext-block"
#define IN_LOOP_ME                      "-in-loop-me"
#define SEARCH_AREA_WIDTH_TOKEN         "-search-w"
#define SEARCH_AREA_HEIGHT_TOKEN        "-search-h"
#define NUM_HME_SEARCH_WIDTH_TOKEN      "-num-hme-w"
#define NUM_HME_SEARCH_HEIGHT_TOKEN     "-num-hme-h"
#define HME_SRCH_T_L0_WIDTH_TOKEN       "-hme-tot-l0-w"
#define HME_SRCH_T_L0_HEIGHT_TOKEN      "-hme-tot-l0-h"
#define HME_LEVEL0_WIDTH                "-hme-l0-w"
#define HME_LEVEL0_HEIGHT               "-hme-l0-h"
#define HME_LEVEL1_WIDTH                "-hme-l1-w"
#define HME_LEVEL1_HEIGHT               "-hme-l1-h"
#define HME_LEVEL2_WIDTH                "-hme-l2-w"
#define HME_LEVEL2_HEIGHT               "-hme-l2-h"
#define SCREEN_CONTENT_TOKEN            "-scm"
// --- start: ALTREF_FILTERING_SUPPORT
#define ENABLE_ALTREFS                  "-enable-altrefs"
#define ALTREF_STRENGTH                 "-altref-strength"
#define ALTREF_NFRAMES                  "-altref-nframes"
#define ENABLE_OVERLAYS                 "-enable-overlays"
// --- end: ALTREF_FILTERING_SUPPORT
#define HBD_MD_ENABLE_TOKEN             "-hbd-md"
#define PALETTE_TOKEN                   "-palette"
#define OLPD_REFINEMENT_TOKEN           "-olpd-refinement"
#define CONSTRAINED_INTRA_ENABLE_TOKEN  "-constrd-intra"
#define HDR_INPUT_TOKEN                 "-hdr"
#define RATE_CONTROL_ENABLE_TOKEN       "-rc"
#define TARGET_BIT_RATE_TOKEN           "-tbr"
#define MAX_QP_TOKEN                    "-max-qp"
#define MIN_QP_TOKEN                    "-min-qp"
#define ADAPTIVE_QP_ENABLE_TOKEN        "-adaptive-quantization"
#define LOOK_AHEAD_DIST_TOKEN           "-lad"
#define SUPER_BLOCK_SIZE_TOKEN          "-sb-size"
#define TILE_ROW_TOKEN                   "-tile-rows"
#define TILE_COL_TOKEN                   "-tile-columns"

#define SQ_WEIGHT_TOKEN                 "-sqw"
#define ENABLE_AMP_TOKEN                "-enable-amp"

#define SCENE_CHANGE_DETECTION_TOKEN    "-scd"
#define INJECTOR_TOKEN                  "-inj"  // no Eval
#define INJECTOR_FRAMERATE_TOKEN        "-inj-frm-rt" // no Eval
#define SPEED_CONTROL_TOKEN             "-speed-ctrl"
#define ASM_TYPE_TOKEN                  "-asm"
#define THREAD_MGMNT                    "-lp"
#define TARGET_SOCKET                   "-ss"
#define UNRESTRICTED_MOTION_VECTOR      "-umv"
#define CONFIG_FILE_COMMENT_CHAR    '#'
#define CONFIG_FILE_NEWLINE_CHAR    '\n'
#define CONFIG_FILE_RETURN_CHAR     '\r'
#define CONFIG_FILE_VALUE_SPLIT     ':'
#define CONFIG_FILE_SPACE_CHAR      ' '
#define CONFIG_FILE_ARRAY_SEP_CHAR  CONFIG_FILE_SPACE_CHAR
#define CONFIG_FILE_TAB_CHAR        '\t'
#define CONFIG_FILE_NULL_CHAR       '\0'
#define CONFIG_FILE_MAX_ARG_COUNT   256
#define CONFIG_FILE_MAX_VAR_LEN     128
#define EVENT_FILE_MAX_ARG_COUNT    20
#define EVENT_FILE_MAX_VAR_LEN      256
#define BUFFER_FILE_MAX_ARG_COUNT   320
#define BUFFER_FILE_MAX_VAR_LEN     128

#define MDS1_PRUNE_C_TH             "-mds1p-class-th"
#define MDS1_PRUNE_S_TH             "-mds1p-cand-th"
#define MDS2_PRUNE_C_TH             "-mds2p-class-th"
#define MDS2_PRUNE_S_TH             "-mds2p-cand-th"

/**********************************
 * Set Cfg Functions
 **********************************/
static void SetCfgInputFile(const char *filename, EbConfig *cfg)
{
    if (cfg->input_file && !cfg->input_file_is_fifo)
        fclose(cfg->input_file);

    if (!filename) {
        cfg->input_file = NULL;
        return;
    }

    if (!strcmp(filename, "stdin"))
        cfg->input_file = stdin;
    else
        FOPEN(cfg->input_file, filename, "rb");

    if (cfg->input_file == NULL) {
        return;
    }

#ifdef _WIN32
    cfg->input_file_is_fifo =
    GetFileType(cfg->input_file) == FILE_TYPE_PIPE;
#else
    int fd = fileno(cfg->input_file);
    struct stat statbuf;
    fstat(fd, &statbuf);
    cfg->input_file_is_fifo = S_ISFIFO(statbuf.st_mode);
#endif

    cfg->y4m_input = check_if_y4m(cfg);
};

static void SetCfgStreamFile                    (const char *value, EbConfig *cfg)
{
    if (cfg->bitstream_file) { fclose(cfg->bitstream_file); }
    FOPEN(cfg->bitstream_file,value, "wb");
};
static void SetCfgErrorFile                     (const char *value, EbConfig *cfg)
{
    if (cfg->error_log_file) { fclose(cfg->error_log_file); }
    FOPEN(cfg->error_log_file,value, "w+");
};
static void SetCfgReconFile                     (const char *value, EbConfig *cfg)
{
    if (cfg->recon_file) { fclose(cfg->recon_file); }
    FOPEN(cfg->recon_file,value, "wb");
};
static void SetCfgQpFile                        (const char *value, EbConfig *cfg)
{
    if (cfg->qp_file) { fclose(cfg->qp_file); }
    FOPEN(cfg->qp_file,value, "r");
};
#if TWO_PASS
static void set_input_stat_file(const char *value, EbConfig *cfg)
{
    if (cfg->input_stat_file) { fclose(cfg->input_stat_file); }
    FOPEN(cfg->input_stat_file, value, "rb");
};
static void set_output_stat_file(const char *value, EbConfig *cfg)
{
    if (cfg->output_stat_file) { fclose(cfg->output_stat_file); }
    FOPEN(cfg->output_stat_file, value, "wb");
};
#if TWO_PASS_USE_2NDP_ME_IN_1STP
static void set_snd_pass_enc_mode(const char *value, EbConfig *cfg) { cfg->snd_pass_enc_mode = (uint8_t)strtoul(value, NULL, 0); };
#endif
#endif
static void SetCfgStatFile(const char *value, EbConfig *cfg)
{
    if (cfg->stat_file) { fclose(cfg->stat_file); }
    FOPEN(cfg->stat_file, value, "wb");
};
static void SetStatReport                       (const char *value, EbConfig *cfg) {cfg->stat_report = (uint8_t) strtoul(value, NULL, 0);};
static void SetCfgSourceWidth                   (const char *value, EbConfig *cfg) {cfg->source_width = strtoul(value, NULL, 0);};
static void SetInterlacedVideo                  (const char *value, EbConfig *cfg) {cfg->interlaced_video  = (EbBool) strtoul(value, NULL, 0);};
static void SetSeperateFields                   (const char *value, EbConfig *cfg) {cfg->separate_fields = (EbBool) strtoul(value, NULL, 0);};
static void SetCfgSourceHeight                  (const char *value, EbConfig *cfg) {cfg->source_height = strtoul(value, NULL, 0) >> cfg->separate_fields;};
static void SetCfgFramesToBeEncoded             (const char *value, EbConfig *cfg) {cfg->frames_to_be_encoded = strtol(value,  NULL, 0) << cfg->separate_fields;};
static void SetBufferedInput                    (const char *value, EbConfig *cfg) {cfg->buffered_input = (strtol(value, NULL, 0) != -1 && cfg->separate_fields) ? strtol(value, NULL, 0) << cfg->separate_fields : strtol(value, NULL, 0);};
static void SetFrameRate                        (const char *value, EbConfig *cfg) {
    cfg->frame_rate = strtoul(value, NULL, 0);
    if (cfg->frame_rate > 1000 )
        cfg->frame_rate = cfg->frame_rate;
    else
        cfg->frame_rate = cfg->frame_rate << 16;
}

static void SetFrameRateNumerator               (const char *value, EbConfig *cfg) { cfg->frame_rate_numerator = strtoul(value, NULL, 0);};
static void SetFrameRateDenominator             (const char *value, EbConfig *cfg) { cfg->frame_rate_denominator = strtoul(value, NULL, 0);};
static void SetEncoderBitDepth                  (const char *value, EbConfig *cfg) {cfg->encoder_bit_depth = strtoul(value, NULL, 0);}
static void SetEncoderColorFormat               (const char *value, EbConfig *cfg) {cfg->encoder_color_format = strtoul(value, NULL, 0);}
static void SetcompressedTenBitFormat           (const char *value, EbConfig *cfg) {cfg->compressed_ten_bit_format = strtoul(value, NULL, 0);}
static void SetBaseLayerSwitchMode              (const char *value, EbConfig *cfg) {cfg->base_layer_switch_mode = (EbBool) strtoul(value, NULL, 0);};
static void SetencMode                          (const char *value, EbConfig *cfg) {cfg->enc_mode = (uint8_t)strtoul(value, NULL, 0);};
static void SetCfgIntraPeriod                   (const char *value, EbConfig *cfg) {cfg->intra_period = strtol(value,  NULL, 0);};
static void SetCfgIntraRefreshType              (const char *value, EbConfig *cfg) {cfg->intra_refresh_type = strtol(value,  NULL, 0);};
static void SetHierarchicalLevels               (const char *value, EbConfig *cfg) { cfg->hierarchical_levels = strtol(value, NULL, 0); };
static void SetCfgPredStructure                 (const char *value, EbConfig *cfg) { cfg->pred_structure = strtol(value, NULL, 0); };
static void SetCfgQp                            (const char *value, EbConfig *cfg) {cfg->qp = strtoul(value, NULL, 0);};
static void SetCfgUseQpFile                     (const char *value, EbConfig *cfg) {cfg->use_qp_file = (EbBool)strtol(value, NULL, 0); };
static void SetCfgFilmGrain                     (const char *value, EbConfig *cfg) { cfg->film_grain_denoise_strength = strtol(value, NULL, 0); };  //not bool to enable possible algorithm extension in the future
static void SetDisableDlfFlag                   (const char *value, EbConfig *cfg) {cfg->disable_dlf_flag = (EbBool)strtoul(value, NULL, 0);};
static void SetEnableLocalWarpedMotionFlag      (const char *value, EbConfig *cfg) {cfg->enable_warped_motion = (EbBool)strtoul(value, NULL, 0);};
static void SetEnableGlobalMotionFlag           (const char *value, EbConfig *cfg) {cfg->enable_global_motion = (EbBool)strtoul(value, NULL, 0);};
static void SetEnableObmcFlag                   (const char *value, EbConfig *cfg) {cfg->enable_obmc = (EbBool)strtoul(value, NULL, 0);};
static void SetEnableRdoqFlag                   (const char *value, EbConfig *cfg) {cfg->enable_rdoq = (int8_t)strtol(value, NULL, 0);};

static void SetEnableFilterIntraFlag            (const char *value, EbConfig *cfg) {cfg->enable_filter_intra = (EbBool)strtoul(value, NULL, 0);};
static void SetEnableHmeFlag                    (const char *value, EbConfig *cfg) {cfg->enable_hme_flag = (EbBool)strtoul(value, NULL, 0);};
static void SetEnableHmeLevel0Flag              (const char *value, EbConfig *cfg) {cfg->enable_hme_level0_flag = (EbBool)strtoul(value, NULL, 0);};
static void SetTileRow                          (const char *value, EbConfig *cfg) { cfg->tile_rows = strtoul(value, NULL, 0); };
static void SetTileCol                          (const char *value, EbConfig *cfg) { cfg->tile_columns = strtoul(value, NULL, 0); };

static void SetSceneChangeDetection             (const char *value, EbConfig *cfg) {cfg->scene_change_detection = strtoul(value, NULL, 0);};
static void SetLookAheadDistance                (const char *value, EbConfig *cfg) {cfg->look_ahead_distance = strtoul(value, NULL, 0);};
static void SetRateControlMode                  (const char *value, EbConfig *cfg) {cfg->rate_control_mode = strtoul(value, NULL, 0);};
static void SetTargetBitRate                    (const char *value, EbConfig *cfg) {cfg->target_bit_rate = 1000*strtoul(value, NULL, 0);};
static void SetMaxQpAllowed                     (const char *value, EbConfig *cfg) {cfg->max_qp_allowed = strtoul(value, NULL, 0);};
static void SetMinQpAllowed                     (const char *value, EbConfig *cfg) {cfg->min_qp_allowed = strtoul(value, NULL, 0);};
static void SetAdaptiveQuantization             (const char *value, EbConfig *cfg) {cfg->enable_adaptive_quantization = (EbBool)strtol(value,  NULL, 0);};
static void SetEnableHmeLevel1Flag              (const char *value, EbConfig *cfg) {cfg->enable_hme_level1_flag  = (EbBool)strtoul(value, NULL, 0);};
static void SetEnableHmeLevel2Flag              (const char *value, EbConfig *cfg) {cfg->enable_hme_level2_flag  = (EbBool)strtoul(value, NULL, 0);};
static void SetCfgSearchAreaWidth               (const char *value, EbConfig *cfg) {cfg->search_area_width = strtoul(value, NULL, 0);};
static void SetCfgSearchAreaHeight              (const char *value, EbConfig *cfg) {cfg->search_area_height = strtoul(value, NULL, 0);};
static void SetCfgNumberHmeSearchRegionInWidth  (const char *value, EbConfig *cfg) {cfg->number_hme_search_region_in_width = strtoul(value, NULL, 0);};
static void SetCfgNumberHmeSearchRegionInHeight (const char *value, EbConfig *cfg) {cfg->number_hme_search_region_in_height = strtoul(value, NULL, 0);};
static void SetCfgHmeLevel0TotalSearchAreaWidth (const char *value, EbConfig *cfg) {cfg->hme_level0_total_search_area_width = strtoul(value, NULL, 0);};
static void SetCfgHmeLevel0TotalSearchAreaHeight(const char *value, EbConfig *cfg) {cfg->hme_level0_total_search_area_height = strtoul(value, NULL, 0);};
static void SetCfgUseDefaultMeHme               (const char *value, EbConfig *cfg) {cfg->use_default_me_hme = (EbBool)strtol(value, NULL, 0); };
static void SetEnableExtBlockFlag(const char *value, EbConfig *cfg) { cfg->ext_block_flag = (EbBool)strtoul(value, NULL, 0); };
static void SetEnableInLoopMeFlag(const char *value, EbConfig *cfg) { cfg->in_loop_me_flag = (EbBool)strtoul(value, NULL, 0); };
static void SetHmeLevel0SearchAreaInWidthArray  (const char *value, EbConfig *cfg) {cfg->hme_level0_search_area_in_width_array[cfg->hme_level0_column_index++] = strtoul(value, NULL, 0);};
static void SetHmeLevel0SearchAreaInHeightArray (const char *value, EbConfig *cfg) {cfg->hme_level0_search_area_in_height_array[cfg->hme_level0_row_index++] = strtoul(value, NULL, 0);};
static void SetHmeLevel1SearchAreaInWidthArray  (const char *value, EbConfig *cfg) {cfg->hme_level1_search_area_in_width_array[cfg->hme_level1_column_index++] = strtoul(value, NULL, 0);};
static void SetHmeLevel1SearchAreaInHeightArray (const char *value, EbConfig *cfg) {cfg->hme_level1_search_area_in_height_array[cfg->hme_level1_row_index++] = strtoul(value, NULL, 0);};
static void SetHmeLevel2SearchAreaInWidthArray  (const char *value, EbConfig *cfg) {cfg->hme_level2_search_area_in_width_array[cfg->hme_level2_column_index++] = strtoul(value, NULL, 0);};
static void SetHmeLevel2SearchAreaInHeightArray (const char *value, EbConfig *cfg) {cfg->hme_level2_search_area_in_height_array[cfg->hme_level2_row_index++] = strtoul(value, NULL, 0);};
static void SetScreenContentMode                (const char *value, EbConfig *cfg) {cfg->screen_content_mode                                                 = strtoul(value, NULL, 0);};
// --- start: ALTREF_FILTERING_SUPPORT
static void SetEnableAltRefs                    (const char *value, EbConfig *cfg) {cfg->enable_altrefs = (EbBool)strtoul(value, NULL, 0);};
static void SetAltRefStrength                   (const char *value, EbConfig *cfg) {cfg->altref_strength = (uint8_t)strtoul(value, NULL, 0);};
static void SetAltRefNFrames                    (const char *value, EbConfig *cfg) {cfg->altref_nframes = (uint8_t)strtoul(value, NULL, 0);};
static void SetEnableOverlays                   (const char *value, EbConfig *cfg) { cfg->enable_overlays = (EbBool)strtoul(value, NULL, 0); };
// --- end: ALTREF_FILTERING_SUPPORT
static void SetEnableHBDModeDecision            (const char *value, EbConfig *cfg) {cfg->enable_hbd_mode_decision = (uint8_t)strtoul(value, NULL, 0);};
static void SetEnablePalette                    (const char *value, EbConfig *cfg) { cfg->enable_palette = (int32_t)strtoul(value, NULL, 0); };
static void SetEnableOlpdRefinement              (const char *value, EbConfig *cfg) { cfg->olpd_refinement = (int32_t)strtoul(value, NULL, 0); };
static void SetEnableConstrainedIntra           (const char *value, EbConfig *cfg) {cfg->constrained_intra                                             = (EbBool)strtoul(value, NULL, 0);};
static void SetHighDynamicRangeInput            (const char *value, EbConfig *cfg) {cfg->high_dynamic_range_input            = strtol(value,  NULL, 0);};
static void SetProfile                          (const char *value, EbConfig *cfg) {cfg->profile                          = strtol(value,  NULL, 0);};
static void SetTier                             (const char *value, EbConfig *cfg) {cfg->tier                             = strtol(value,  NULL, 0);};
static void SetLevel                            (const char *value, EbConfig *cfg) {

    if (strtoul( value, NULL,0) != 0 || EB_STRCMP(value, "0") == 0 )
        cfg->level = (uint32_t)(10*strtod(value,  NULL));
    else
        cfg->level = 9999999;
};
static void SetInjector                         (const char *value, EbConfig *cfg) {cfg->injector                         = strtol(value,  NULL, 0);};
static void SpeedControlFlag                    (const char *value, EbConfig *cfg) { cfg->speed_control_flag = strtol(value, NULL, 0); };
static void SetInjectorFrameRate                (const char *value, EbConfig *cfg) {
    cfg->injector_frame_rate = strtoul(value, NULL, 0);
    if (cfg->injector_frame_rate > 1000 )
        cfg->injector_frame_rate = cfg->injector_frame_rate;
    else
        cfg->injector_frame_rate = cfg->injector_frame_rate << 16;
}
static void SetLatencyMode                      (const char *value, EbConfig *cfg)  {cfg->latency_mode               = (uint8_t)strtol(value, NULL, 0);};
static void SetAsmType                          (const char *value, EbConfig *cfg)  {
    const struct {
        char *name;
        CPU_FLAGS flags;
    } param_maps[] ={
        {"c",       0},                             {"0",  0},
        {"mmx",     (CPU_FLAGS_MMX << 1) - 1},      {"1",  (CPU_FLAGS_MMX << 1) - 1},
        {"sse",     (CPU_FLAGS_SSE << 1) - 1},      {"2",  (CPU_FLAGS_SSE << 1) - 1},
        {"sse2",    (CPU_FLAGS_SSE2 << 1) - 1},     {"3",  (CPU_FLAGS_SSE2 << 1) - 1},
        {"sse3",    (CPU_FLAGS_SSE3 << 1) - 1},     {"4",  (CPU_FLAGS_SSE3 << 1) - 1},
        {"ssse3",   (CPU_FLAGS_SSSE3 << 1) - 1},    {"5",  (CPU_FLAGS_SSSE3 << 1) - 1},
        {"sse4_1",  (CPU_FLAGS_SSE4_1 << 1) - 1},   {"6",  (CPU_FLAGS_SSE4_1 << 1) - 1},
        {"sse4_2",  (CPU_FLAGS_SSE4_2 << 1) - 1},   {"7",  (CPU_FLAGS_SSE4_2 << 1) - 1},
        {"avx",     (CPU_FLAGS_AVX << 1) - 1},      {"8",  (CPU_FLAGS_AVX << 1) - 1},
        {"avx2",    (CPU_FLAGS_AVX2 << 1) - 1},     {"9",  (CPU_FLAGS_AVX2 << 1) - 1},
        {"avx512",  (CPU_FLAGS_AVX512VL << 1) - 1}, {"10", (CPU_FLAGS_AVX512VL << 1) - 1},
        {"max",     CPU_FLAGS_ALL},                 {"11", CPU_FLAGS_ALL},
    };
    const uint32_t para_map_size = sizeof(param_maps) / sizeof(param_maps[0]);
    uint32_t i;

    for (i = 0; i < para_map_size; ++i) {
        if (EB_STRCMP(value, param_maps[i].name) == 0) {
            cfg->cpu_flags_limit = param_maps[i].flags;
            return;
        }
    }

    cfg->cpu_flags_limit = CPU_FLAGS_INVALID;
};
static void SetLogicalProcessors                (const char *value, EbConfig *cfg)  {cfg->logical_processors         = (uint32_t)strtoul(value, NULL, 0);};
static void SetTargetSocket                     (const char *value, EbConfig *cfg)  {cfg->target_socket              = (int32_t)strtol(value, NULL, 0);};
static void SetUnrestrictedMotionVector         (const char *value, EbConfig *cfg)  {cfg->unrestricted_motion_vector = (EbBool)strtol(value, NULL, 0);};

static void SetSquareWeight                     (const char *value, EbConfig *cfg)  {cfg->sq_weight                  = (uint64_t)strtoul(value, NULL, 0);
        if (cfg->sq_weight == 0)
            cfg->sq_weight = (uint32_t)~0;
}

static void SetMDS1_PRUNE_C_TH(const char *value, EbConfig *cfg) {
    cfg->md_stage_1_class_prune_th = (uint64_t)strtoul(value, NULL, 0);
    if (cfg->md_stage_1_class_prune_th == 0)
        cfg->md_stage_1_class_prune_th = (uint64_t)~0;
}
static void SetMDS1_PRUNE_S_TH(const char *value, EbConfig *cfg) {
    cfg->md_stage_1_cand_prune_th = (uint64_t)strtoul(value, NULL, 0);
    if (cfg->md_stage_1_cand_prune_th == 0)
        cfg->md_stage_1_cand_prune_th = (uint64_t)~0;
}
static void SetMDS2_PRUNE_C_TH(const char *value, EbConfig *cfg) {
    cfg->md_stage_2_class_prune_th = (uint64_t)strtoul(value, NULL, 0);
    if (cfg->md_stage_2_class_prune_th == 0)
        cfg->md_stage_2_class_prune_th = (uint64_t)~0;
}
static void SetMDS2_PRUNE_S_TH(const char *value, EbConfig *cfg) {
    cfg->md_stage_2_cand_prune_th = (uint64_t)strtoul(value, NULL, 0);
    if (cfg->md_stage_2_cand_prune_th == 0)
        cfg->md_stage_2_cand_prune_th = (uint64_t)~0;
}

static void set_enable_auto_max_partition           (const char *value, EbConfig *cfg) { cfg->enable_auto_max_partition = (uint8_t)strtol(value, NULL, 0); };

enum cfg_type{
    SINGLE_INPUT,   // Configuration parameters that have only 1 value input
    ARRAY_INPUT     // Configuration parameters that have multiple values as input
};

/**********************************
 * Config Entry Struct
 **********************************/
typedef struct config_entry_s {
    enum  cfg_type type;
    const char *token;
    const char *name;
    void (*scf)(const char *, EbConfig *);
} config_entry_t;

/**********************************
 * Config Entry Array
 **********************************/
config_entry_t config_entry[] = {
    // File I/O
    { SINGLE_INPUT, INPUT_FILE_TOKEN, "InputFile", SetCfgInputFile },
    { SINGLE_INPUT, OUTPUT_BITSTREAM_TOKEN,   "StreamFile",       SetCfgStreamFile },
    { SINGLE_INPUT, ERROR_FILE_TOKEN, "ErrorFile", SetCfgErrorFile },
    { SINGLE_INPUT, OUTPUT_RECON_TOKEN, "ReconFile", SetCfgReconFile },
    { SINGLE_INPUT, QP_FILE_TOKEN, "QpFile", SetCfgQpFile },
    { SINGLE_INPUT, STAT_FILE_TOKEN, "StatFile", SetCfgStatFile },
#if TWO_PASS
    { SINGLE_INPUT, INPUT_STAT_FILE_TOKEN, "input_stat_file", set_input_stat_file },
    { SINGLE_INPUT, OUTPUT_STAT_FILE_TOKEN, "output_stat_file", set_output_stat_file },
#endif

    // Interlaced Video
    { SINGLE_INPUT, INTERLACED_VIDEO_TOKEN , "InterlacedVideo" , SetInterlacedVideo },
    { SINGLE_INPUT, SEPERATE_FILDS_TOKEN, "SeperateFields", SetSeperateFields },
    // Picture Dimensions
    { SINGLE_INPUT, WIDTH_TOKEN, "SourceWidth", SetCfgSourceWidth },
    { SINGLE_INPUT, HEIGHT_TOKEN, "SourceHeight", SetCfgSourceHeight },
    // Prediction Structure
    { SINGLE_INPUT, NUMBER_OF_PICTURES_TOKEN, "FrameToBeEncoded", SetCfgFramesToBeEncoded },
    { SINGLE_INPUT, BUFFERED_INPUT_TOKEN, "BufferedInput", SetBufferedInput },
    { SINGLE_INPUT, BASE_LAYER_SWITCH_MODE_TOKEN, "BaseLayerSwitchMode", SetBaseLayerSwitchMode },
    { SINGLE_INPUT, ENCMODE_TOKEN, "EncoderMode", SetencMode},
#if TWO_PASS_USE_2NDP_ME_IN_1STP
    { SINGLE_INPUT, ENCMODE2P_TOKEN, "EncoderMode2p", set_snd_pass_enc_mode},
#endif
    { SINGLE_INPUT, INTRA_PERIOD_TOKEN, "IntraPeriod", SetCfgIntraPeriod },
    { SINGLE_INPUT, INTRA_REFRESH_TYPE_TOKEN, "IntraRefreshType", SetCfgIntraRefreshType },
    { SINGLE_INPUT, FRAME_RATE_TOKEN, "FrameRate", SetFrameRate },
    { SINGLE_INPUT, FRAME_RATE_NUMERATOR_TOKEN, "FrameRateNumerator", SetFrameRateNumerator },
    { SINGLE_INPUT, FRAME_RATE_DENOMINATOR_TOKEN, "FrameRateDenominator", SetFrameRateDenominator },
    { SINGLE_INPUT, ENCODER_BIT_DEPTH, "EncoderBitDepth", SetEncoderBitDepth },
    { SINGLE_INPUT, ENCODER_COLOR_FORMAT, "EncoderColorFormat", SetEncoderColorFormat},
    { SINGLE_INPUT, INPUT_COMPRESSED_TEN_BIT_FORMAT, "CompressedTenBitFormat", SetcompressedTenBitFormat },
    { SINGLE_INPUT, HIERARCHICAL_LEVELS_TOKEN, "HierarchicalLevels", SetHierarchicalLevels },
    { SINGLE_INPUT, PRED_STRUCT_TOKEN, "PredStructure", SetCfgPredStructure },
     { SINGLE_INPUT, TILE_ROW_TOKEN, "TileRow", SetTileRow},
     { SINGLE_INPUT, TILE_COL_TOKEN, "TileCol", SetTileCol},
    // Rate Control
    { SINGLE_INPUT, SCENE_CHANGE_DETECTION_TOKEN, "SceneChangeDetection", SetSceneChangeDetection},
    { SINGLE_INPUT, QP_TOKEN, "QP", SetCfgQp },
    { SINGLE_INPUT, USE_QP_FILE_TOKEN, "UseQpFile", SetCfgUseQpFile },
    { SINGLE_INPUT, STAT_REPORT_TOKEN, "StatReport", SetStatReport },
    { SINGLE_INPUT, RATE_CONTROL_ENABLE_TOKEN, "RateControlMode", SetRateControlMode },
    { SINGLE_INPUT, LOOK_AHEAD_DIST_TOKEN, "LookAheadDistance",                             SetLookAheadDistance},
    { SINGLE_INPUT, TARGET_BIT_RATE_TOKEN, "TargetBitRate", SetTargetBitRate },
    { SINGLE_INPUT, MAX_QP_TOKEN, "MaxQpAllowed", SetMaxQpAllowed },
    { SINGLE_INPUT, MIN_QP_TOKEN, "MinQpAllowed", SetMinQpAllowed },
    { SINGLE_INPUT, ADAPTIVE_QP_ENABLE_TOKEN, "AdaptiveQuantization", SetAdaptiveQuantization },

    // DLF
    { SINGLE_INPUT, LOOP_FILTER_DISABLE_TOKEN, "LoopFilterDisable", SetDisableDlfFlag },
    // LOCAL WARPED MOTION
    { SINGLE_INPUT, LOCAL_WARPED_ENABLE_TOKEN, "LocalWarpedMotion", SetEnableLocalWarpedMotionFlag },
    // GLOBAL MOTION
    { SINGLE_INPUT, GLOBAL_MOTION_ENABLE_TOKEN, "GlobalMotion", SetEnableGlobalMotionFlag },
    // OBMC
    { SINGLE_INPUT, OBMC_TOKEN, "Obmc", SetEnableObmcFlag },
    // RDOQ
    { SINGLE_INPUT, RDOQ_TOKEN, "RDOQ", SetEnableRdoqFlag },
    // Filter Intra
    { SINGLE_INPUT, FILTER_INTRA_TOKEN, "FilterIntra", SetEnableFilterIntraFlag },
    // ME Tools
    { SINGLE_INPUT, USE_DEFAULT_ME_HME_TOKEN, "UseDefaultMeHme", SetCfgUseDefaultMeHme },
    { SINGLE_INPUT, HME_ENABLE_TOKEN, "HME", SetEnableHmeFlag },
    { SINGLE_INPUT, HME_L0_ENABLE_TOKEN, "HMELevel0", SetEnableHmeLevel0Flag },
    { SINGLE_INPUT, HME_L1_ENABLE_TOKEN, "HMELevel1", SetEnableHmeLevel1Flag },
    { SINGLE_INPUT, HME_L2_ENABLE_TOKEN, "HMELevel2", SetEnableHmeLevel2Flag },
    { SINGLE_INPUT, EXT_BLOCK, "ExtBlockFlag", SetEnableExtBlockFlag },
    { SINGLE_INPUT, IN_LOOP_ME, "InLoopMeFlag", SetEnableInLoopMeFlag },
    // ME Parameters
    { SINGLE_INPUT, SEARCH_AREA_WIDTH_TOKEN, "SearchAreaWidth", SetCfgSearchAreaWidth },
    { SINGLE_INPUT, SEARCH_AREA_HEIGHT_TOKEN, "SearchAreaHeight", SetCfgSearchAreaHeight },
    // HME Parameters
    { SINGLE_INPUT, NUM_HME_SEARCH_WIDTH_TOKEN, "NumberHmeSearchRegionInWidth", SetCfgNumberHmeSearchRegionInWidth },
    { SINGLE_INPUT, NUM_HME_SEARCH_HEIGHT_TOKEN, "NumberHmeSearchRegionInHeight", SetCfgNumberHmeSearchRegionInHeight },
    { SINGLE_INPUT, HME_SRCH_T_L0_WIDTH_TOKEN, "HmeLevel0TotalSearchAreaWidth", SetCfgHmeLevel0TotalSearchAreaWidth },
    { SINGLE_INPUT, HME_SRCH_T_L0_HEIGHT_TOKEN, "HmeLevel0TotalSearchAreaHeight", SetCfgHmeLevel0TotalSearchAreaHeight },
    // MD Parameters
    { SINGLE_INPUT, SCREEN_CONTENT_TOKEN, "ScreenContentMode", SetScreenContentMode},
    { SINGLE_INPUT, HBD_MD_ENABLE_TOKEN, "HighBitDepthModeDecision", SetEnableHBDModeDecision },
    { SINGLE_INPUT, PALETTE_TOKEN, "PaletteMode", SetEnablePalette },
    { SINGLE_INPUT, OLPD_REFINEMENT_TOKEN, "OlpdRefinement", SetEnableOlpdRefinement },
    { SINGLE_INPUT, CONSTRAINED_INTRA_ENABLE_TOKEN, "ConstrainedIntra", SetEnableConstrainedIntra},
    // Thread Management
    { SINGLE_INPUT, THREAD_MGMNT, "logicalProcessors", SetLogicalProcessors },
    { SINGLE_INPUT, TARGET_SOCKET, "TargetSocket", SetTargetSocket },
    // Optional Features
    { SINGLE_INPUT, UNRESTRICTED_MOTION_VECTOR, "UnrestrictedMotionVector", SetUnrestrictedMotionVector },

//    { SINGLE_INPUT, BITRATE_REDUCTION_TOKEN, "bit_rate_reduction", SetBitRateReduction },
    { SINGLE_INPUT, HDR_INPUT_TOKEN, "HighDynamicRangeInput", SetHighDynamicRangeInput },
    // Latency
    { SINGLE_INPUT, INJECTOR_TOKEN, "Injector", SetInjector },
    { SINGLE_INPUT, INJECTOR_FRAMERATE_TOKEN, "InjectorFrameRate", SetInjectorFrameRate },
    { SINGLE_INPUT, SPEED_CONTROL_TOKEN, "SpeedControlFlag", SpeedControlFlag },
    // Annex A parameters
    { SINGLE_INPUT, PROFILE_TOKEN, "Profile", SetProfile },
    { SINGLE_INPUT, TIER_TOKEN, "Tier", SetTier },
    { SINGLE_INPUT, LEVEL_TOKEN, "Level", SetLevel },
    { SINGLE_INPUT, LATENCY_MODE, "LatencyMode", SetLatencyMode },
    { SINGLE_INPUT, FILM_GRAIN_TOKEN, "FilmGrain", SetCfgFilmGrain },
    // Asm Type
    { SINGLE_INPUT, ASM_TYPE_TOKEN, "Asm", SetAsmType },
    // HME
    { ARRAY_INPUT,HME_LEVEL0_WIDTH, "HmeLevel0SearchAreaInWidth", SetHmeLevel0SearchAreaInWidthArray },
    { ARRAY_INPUT,HME_LEVEL0_HEIGHT, "HmeLevel0SearchAreaInHeight", SetHmeLevel0SearchAreaInHeightArray },
    { ARRAY_INPUT,HME_LEVEL1_WIDTH, "HmeLevel1SearchAreaInWidth", SetHmeLevel1SearchAreaInWidthArray },
    { ARRAY_INPUT,HME_LEVEL1_HEIGHT, "HmeLevel1SearchAreaInHeight", SetHmeLevel1SearchAreaInHeightArray },
    { ARRAY_INPUT,HME_LEVEL2_WIDTH, "HmeLevel2SearchAreaInWidth", SetHmeLevel2SearchAreaInWidthArray },
    { ARRAY_INPUT,HME_LEVEL2_HEIGHT, "HmeLevel2SearchAreaInHeight", SetHmeLevel2SearchAreaInHeightArray },
    // --- start: ALTREF_FILTERING_SUPPORT
    { SINGLE_INPUT, ENABLE_ALTREFS, "EnableAltRefs", SetEnableAltRefs },
    { SINGLE_INPUT, ALTREF_STRENGTH, "AltRefStrength", SetAltRefStrength },
    { SINGLE_INPUT, ALTREF_NFRAMES, "AltRefNframes", SetAltRefNFrames },
    { SINGLE_INPUT, ENABLE_OVERLAYS, "EnableOverlays", SetEnableOverlays },
    // --- end: ALTREF_FILTERING_SUPPORT

    { SINGLE_INPUT, SQ_WEIGHT_TOKEN, "SquareWeight", SetSquareWeight },
    { SINGLE_INPUT, ENABLE_AMP_TOKEN, "AutomaxPartition", set_enable_auto_max_partition },

    { SINGLE_INPUT, MDS1_PRUNE_C_TH, "MDStage1PruneClassThreshold", SetMDS1_PRUNE_C_TH },
    { SINGLE_INPUT, MDS1_PRUNE_S_TH, "MDStage1PruneCandThreshold", SetMDS1_PRUNE_S_TH },
    { SINGLE_INPUT, MDS2_PRUNE_C_TH, "MDStage2PruneClassThreshold", SetMDS2_PRUNE_C_TH },
    { SINGLE_INPUT, MDS2_PRUNE_S_TH, "MDStage2PruneCandThreshold", SetMDS2_PRUNE_S_TH },

    // Termination
    {SINGLE_INPUT,NULL,  NULL,                                NULL}
};

/**********************************
 * Constructor
 **********************************/
void eb_config_ctor(EbConfig *config_ptr)
{
    memset(config_ptr, 0, sizeof(*config_ptr));
    config_ptr->error_log_file                         = stderr;
    config_ptr->frame_rate                            = 30 << 16;
    config_ptr->encoder_bit_depth                      = 8;
    config_ptr->encoder_color_format                   = 1; //EB_YUV420
    config_ptr->buffered_input                        = -1;

    config_ptr->qp                                   = 50;
    config_ptr->use_qp_file                          = EB_FALSE;
    config_ptr->look_ahead_distance                  = (uint32_t)~0;
    config_ptr->target_bit_rate                        = 7000000;
    config_ptr->max_qp_allowed                       = 63;
    config_ptr->min_qp_allowed                       = 10;

    config_ptr->enable_adaptive_quantization         = 2;
    config_ptr->enc_mode                              = MAX_ENC_PRESET;
#if TWO_PASS_USE_2NDP_ME_IN_1STP
    config_ptr->snd_pass_enc_mode                     = MAX_ENC_PRESET + 1;
#endif
    config_ptr->intra_period                          = -2;
    config_ptr->intra_refresh_type                     = 1;
    config_ptr->hierarchical_levels                   = 4;
    config_ptr->pred_structure                        = 2;
    config_ptr->enable_global_motion                 = EB_TRUE;
    config_ptr->enable_obmc                          = EB_TRUE;
    config_ptr->enable_rdoq                          = -1;
    config_ptr->enable_filter_intra                  = EB_TRUE;
    config_ptr->in_loop_me_flag                      = EB_TRUE;
    config_ptr->use_default_me_hme                   = EB_TRUE;
    config_ptr->enable_hme_flag                        = EB_TRUE;
    config_ptr->enable_hme_level0_flag                  = EB_TRUE;
    config_ptr->search_area_width                      = 16;
    config_ptr->search_area_height                     = 7;
    config_ptr->number_hme_search_region_in_width         = 2;
    config_ptr->number_hme_search_region_in_height        = 2;
    config_ptr->hme_level0_total_search_area_width        = 64;
    config_ptr->hme_level0_total_search_area_height       = 25;
    config_ptr->hme_level0_search_area_in_width_array[0]   = 32;
    config_ptr->hme_level0_search_area_in_width_array[1]   = 32;
    config_ptr->hme_level0_search_area_in_height_array[0]  = 12;
    config_ptr->hme_level0_search_area_in_height_array[1]  = 13;
    config_ptr->hme_level1_search_area_in_width_array[0]   = 1;
    config_ptr->hme_level1_search_area_in_width_array[1]   = 1;
    config_ptr->hme_level1_search_area_in_height_array[0]  = 1;
    config_ptr->hme_level1_search_area_in_height_array[1]  = 1;
    config_ptr->hme_level2_search_area_in_width_array[0]   = 1;
    config_ptr->hme_level2_search_area_in_width_array[1]   = 1;
    config_ptr->hme_level2_search_area_in_height_array[0]  = 1;
    config_ptr->hme_level2_search_area_in_height_array[1]  = 1;
    config_ptr->screen_content_mode                  = 2;
    config_ptr->enable_hbd_mode_decision             = 1;
    config_ptr->enable_palette                       = -1;
    config_ptr->olpd_refinement                      = -1;
    config_ptr->injector_frame_rate                    = 60 << 16;

    // ASM Type
    config_ptr->cpu_flags_limit                         = CPU_FLAGS_ALL;

    config_ptr->target_socket                         = -1;

    config_ptr->unrestricted_motion_vector           = EB_TRUE;

    // --- start: ALTREF_FILTERING_SUPPORT
    config_ptr->enable_altrefs                       = EB_TRUE;
    config_ptr->altref_strength                      = 5;
    config_ptr->altref_nframes                       = 7;
    // --- end: ALTREF_FILTERING_SUPPORT

    config_ptr->sq_weight                            = 100;
    config_ptr->enable_auto_max_partition             = 1;

    config_ptr->md_stage_1_cand_prune_th                = 75;
    config_ptr->md_stage_1_class_prune_th                = 100;
    config_ptr->md_stage_2_cand_prune_th                = 15;
    config_ptr->md_stage_2_class_prune_th                = 25;

    return;
}

/**********************************
 * Destructor
 **********************************/
void eb_config_dtor(EbConfig *config_ptr)
{
    // Close any files that are open
    if (config_ptr->config_file) {
        fclose(config_ptr->config_file);
        config_ptr->config_file = (FILE *) NULL;
    }

    if (config_ptr->input_file) {
        if (!config_ptr->input_file_is_fifo)
            fclose(config_ptr->input_file);
        config_ptr->input_file = (FILE *) NULL;
    }

    if (config_ptr->bitstream_file) {
        fclose(config_ptr->bitstream_file);
        config_ptr->bitstream_file = (FILE *) NULL;
    }

    if (config_ptr->recon_file) {
        fclose(config_ptr->recon_file);
        config_ptr->recon_file = (FILE *)NULL;
    }

    if (config_ptr->error_log_file) {
        fclose(config_ptr->error_log_file);
        config_ptr->error_log_file = (FILE *) NULL;
    }

    if (config_ptr->qp_file) {
        fclose(config_ptr->qp_file);
        config_ptr->qp_file = (FILE *)NULL;
    }

    if (config_ptr->stat_file) {
        fclose(config_ptr->stat_file);
        config_ptr->stat_file = (FILE *) NULL;
    }
#if TWO_PASS
    if (config_ptr->input_stat_file) {
        fclose(config_ptr->input_stat_file);
        config_ptr->input_stat_file = (FILE *)NULL;
    }
    if (config_ptr->output_stat_file) {
        fclose(config_ptr->output_stat_file);
        config_ptr->output_stat_file = (FILE *)NULL;
    }
#endif
    return;
}

/**********************************
 * File Size
 **********************************/
static int32_t findFileSize(
    FILE * const pFile)
{
    int32_t fileSize;

    fseek(pFile, 0, SEEK_END);
    fileSize = ftell(pFile);
    rewind(pFile);

    return fileSize;
}

/**********************************
 * Line Split
 **********************************/
static void lineSplit(
    uint32_t       *argc,
    char           *argv  [CONFIG_FILE_MAX_ARG_COUNT],
    uint32_t        argLen[CONFIG_FILE_MAX_ARG_COUNT],
    char           *linePtr)
{
    uint32_t i=0;
    *argc = 0;

    while((*linePtr != CONFIG_FILE_NEWLINE_CHAR) &&
            (*linePtr != CONFIG_FILE_RETURN_CHAR) &&
            (*linePtr != CONFIG_FILE_COMMENT_CHAR) &&
            (*argc < CONFIG_FILE_MAX_ARG_COUNT)) {
        // Increment past whitespace
        while((*linePtr == CONFIG_FILE_SPACE_CHAR || *linePtr == CONFIG_FILE_TAB_CHAR) && (*linePtr != CONFIG_FILE_NEWLINE_CHAR))
            ++linePtr;
        // Set arg
        if ((*linePtr != CONFIG_FILE_NEWLINE_CHAR) &&
                (*linePtr != CONFIG_FILE_RETURN_CHAR) &&
                (*linePtr != CONFIG_FILE_COMMENT_CHAR) &&
                (*argc < CONFIG_FILE_MAX_ARG_COUNT)) {
            argv[*argc] = linePtr;

            // Increment to next whitespace
            while(*linePtr != CONFIG_FILE_SPACE_CHAR &&
                    *linePtr != CONFIG_FILE_TAB_CHAR &&
                    *linePtr != CONFIG_FILE_NEWLINE_CHAR &&
                    *linePtr != CONFIG_FILE_RETURN_CHAR) {
                ++linePtr;
                ++i;
            }

            // Set arg length
            argLen[(*argc)++] = i;

            i=0;
        }
    }

    return;
}

/**********************************
* Set Config value
**********************************/
static void SetConfigValue(EbConfig *config, const char *name,
                           const char *value) {
    int32_t i=0;

    while(config_entry[i].name != NULL) {
        if(EB_STRCMP(config_entry[i].name, name) == 0)
            (*config_entry[i].scf)((const char *) value, config);
        ++i;
    }

    return;
}

/**********************************
* Parse Config File
**********************************/
static void ParseConfigFile(
    EbConfig *config,
    char       *buffer,
    int32_t         size)
{
    uint32_t argc;
    char *argv[CONFIG_FILE_MAX_ARG_COUNT];
    uint32_t argLen[CONFIG_FILE_MAX_ARG_COUNT];

    char varName[CONFIG_FILE_MAX_VAR_LEN];
    char varValue[CONFIG_FILE_MAX_ARG_COUNT][CONFIG_FILE_MAX_VAR_LEN];

    uint32_t valueIndex;

    uint32_t commentSectionFlag = 0;
    uint32_t newLineFlag = 0;

    // Keep looping until we process the entire file
    while(size--) {
        commentSectionFlag = ((*buffer == CONFIG_FILE_COMMENT_CHAR) || (commentSectionFlag != 0)) ? 1 : commentSectionFlag;

        // At the beginning of each line
        if ((newLineFlag == 1) && (commentSectionFlag == 0)) {
            // Do an argc/argv split for the line
            lineSplit(&argc, argv, argLen, buffer);

            if ((argc > 2) && (*argv[1] == CONFIG_FILE_VALUE_SPLIT)) {
                // ***NOTE - We're assuming that the variable name is the first arg and
                // the variable value is the third arg.

                // Cap the length of the variable name
                argLen[0] = (argLen[0] > CONFIG_FILE_MAX_VAR_LEN - 1) ? CONFIG_FILE_MAX_VAR_LEN - 1 : argLen[0];
                // Copy the variable name
                EB_STRNCPY(varName, CONFIG_FILE_MAX_VAR_LEN, argv[0], argLen[0]);
                // Null terminate the variable name
                varName[argLen[0]] = CONFIG_FILE_NULL_CHAR;

                for(valueIndex=0; (valueIndex < CONFIG_FILE_MAX_ARG_COUNT - 2) && (valueIndex < (argc - 2)); ++valueIndex) {
                    // Cap the length of the variable
                    argLen[valueIndex+2] = (argLen[valueIndex+2] > CONFIG_FILE_MAX_VAR_LEN - 1) ? CONFIG_FILE_MAX_VAR_LEN - 1 : argLen[valueIndex+2];
                    // Copy the variable name
                    EB_STRNCPY(varValue[valueIndex], CONFIG_FILE_MAX_VAR_LEN, argv[valueIndex+2], argLen[valueIndex+2]);
                    // Null terminate the variable name
                    varValue[valueIndex][argLen[valueIndex+2]] = CONFIG_FILE_NULL_CHAR;

                    SetConfigValue(config, varName, varValue[valueIndex]);
                }
            }
        }

        commentSectionFlag = (*buffer == CONFIG_FILE_NEWLINE_CHAR) ? 0 : commentSectionFlag;
        newLineFlag = (*buffer == CONFIG_FILE_NEWLINE_CHAR) ? 1 : 0;
        ++buffer;
    }

    return;
}

/******************************************
* Find Token
******************************************/
static int32_t FindToken(
    int32_t         argc,
    char * const    argv[],
    char const *    token,
    char*           configStr)
{
    int32_t return_error = -1;

    while((argc > 0) && (return_error != 0)) {
        return_error = EB_STRCMP(argv[--argc], token);
        if (return_error == 0) {
            EB_STRCPY(configStr, COMMAND_LINE_MAX_SIZE, argv[argc + 1]);
        }
    }

    return return_error;
}

/**********************************
* Read Config File
**********************************/
static int32_t ReadConfigFile(
    EbConfig  *config,
    char        *configPath,
    uint32_t     instance_idx)
{
    int32_t return_error = 0;

    // Open the config file
    FOPEN(config->config_file, configPath, "rb");

    if (config->config_file != (FILE*) NULL) {
        int32_t configFileSize = findFileSize(config->config_file);
        char *configFileBuffer = (char*) malloc(configFileSize);

        if (configFileBuffer != (char *) NULL) {
            int32_t resultSize = (int32_t) fread(configFileBuffer, 1, configFileSize, config->config_file);

            if (resultSize == configFileSize) {
                ParseConfigFile(config, configFileBuffer, configFileSize);
            } else {
                printf("Error channel %u: File Read Failed\n",instance_idx+1);
                return_error = -1;
            }
        } else {
            printf("Error channel %u: Memory Allocation Failed\n",instance_idx+1);
            return_error = -1;
        }

        free(configFileBuffer);
        fclose(config->config_file);
        config->config_file = (FILE*) NULL;
    } else {
        printf("Error channel %u: Couldn't open Config File: %s\n", instance_idx+1,configPath);
        return_error = -1;
    }

    return return_error;
}

/******************************************
* Verify Settings
******************************************/
static EbErrorType VerifySettings(EbConfig *config, uint32_t channelNumber)
{
    EbErrorType return_error = EB_ErrorNone;

    // Check Input File
    if(config->input_file == (FILE*) NULL) {
        fprintf(config->error_log_file, "Error instance %u: Invalid Input File\n",channelNumber+1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->frames_to_be_encoded <= -1) {
        fprintf(config->error_log_file, "Error instance %u: FrameToBeEncoded must be greater than 0\n",channelNumber+1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->buffered_input < -1) {
        fprintf(config->error_log_file, "Error instance %u: Invalid buffered_input. buffered_input must greater or equal to -1\n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->buffered_input > config->frames_to_be_encoded) {
        fprintf(config->error_log_file, "Error instance %u: Invalid buffered_input. buffered_input must be less or equal to the number of frames to be encoded\n",channelNumber+1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->use_qp_file == EB_TRUE && config->qp_file == NULL) {
        fprintf(config->error_log_file, "Error instance %u: Could not find QP file, UseQpFile is set to 1\n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }
    if (config->separate_fields > 1) {
        fprintf(config->error_log_file, "Error Instance %u: Invalid SeperateFields Input\n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->encoder_bit_depth == 10 && config->separate_fields == 1)
    {
        fprintf(config->error_log_file, "Error instance %u: Separate fields is not supported for 10 bit input \n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->encoder_color_format != 1) {
        fprintf(config->error_log_file, "Error instance %u: Only support 420 now \n", channelNumber + 1);
        return_error = EB_ErrorBadParameter;
    }
    if (config->injector > 1 ){
        fprintf(config->error_log_file, "Error Instance %u: Invalid injector [0 - 1]\n",channelNumber+1);
        return_error = EB_ErrorBadParameter;
    }

    if(config->injector_frame_rate > (240<<16) && config->injector){
        fprintf(config->error_log_file, "Error Instance %u: The maximum allowed injector_frame_rate is 240 fps\n",channelNumber+1);
        return_error = EB_ErrorBadParameter;
    }
    // Check that the injector frame_rate is non-zero
    if(config->injector_frame_rate <= 0 && config->injector) {
        fprintf(config->error_log_file, "Error Instance %u: The injector frame rate should be greater than 0 fps \n",channelNumber+1);
        return_error = EB_ErrorBadParameter;
    }

    // target_socket
    if (config->target_socket != -1 && config->target_socket != 0 && config->target_socket != 1) {
        fprintf(config->error_log_file, "Error instance %u: Invalid target_socket [-1 - 1], your input: %d\n", channelNumber + 1, config->target_socket);
        return_error = EB_ErrorBadParameter;
    }

    // Local Warped Motion
    if (config->enable_warped_motion != 0 && config->enable_warped_motion != 1) {
        fprintf(config->error_log_file, "Error instance %u: Invalid warped motion flag [0 - 1], your input: %d\n", channelNumber + 1, config->target_socket);
        return_error = EB_ErrorBadParameter;
    }

    // Global Motion
    if (config->enable_global_motion != 0 && config->enable_global_motion != 1) {
        fprintf(config->error_log_file, "Error instance %u: Invalid global motion flag [0 - 1], your input: %d\n", channelNumber + 1, config->target_socket);
        return_error = EB_ErrorBadParameter;
    }

    // OBMC
    if (config->enable_obmc != 0 && config->enable_obmc != 1) {
        fprintf(config->error_log_file, "Error instance %u: Invalid OBMC flag [0 - 1], your input: %d\n", channelNumber + 1, config->target_socket);
        return_error = EB_ErrorBadParameter;
    }

    // Filter Intra prediction
    if (config->enable_filter_intra != 0 && config->enable_filter_intra != 1) {
        fprintf(config->error_log_file, "Error instance %u: Invalid Filter Intra flag [0 - 1], your input: %d\n", channelNumber + 1, config->target_socket);
        return_error = EB_ErrorBadParameter;
    }

    // HBD mode decision
    if (config->enable_hbd_mode_decision != 0 && config->enable_hbd_mode_decision != 1 && config->enable_hbd_mode_decision != 2) {
        fprintf(config->error_log_file, "Error instance %u: Invalid HBD mode decision flag [0 - 2], your input: %d\n", channelNumber + 1, config->target_socket);
        return_error = EB_ErrorBadParameter;
    }



    return return_error;
}

/******************************************
 * Find Token for multiple inputs
 ******************************************/
int32_t FindTokenMultipleInputs(
    int32_t         argc,
    char* const     argv[],
    const char*     token,
    char**          configStr)
{
    int32_t return_error = -1;
    int32_t done = 0;
    while((argc > 0) && (return_error != 0)) {
        return_error = EB_STRCMP(argv[--argc], token);
        if (return_error == 0) {
            int32_t count;
            for (count=0; count < MAX_CHANNEL_NUMBER  ; ++count){
                if (done ==0){
                    if (argv[argc + count + 1] ){
                        if (strtoul(argv[argc + count + 1], NULL,0) != 0 || EB_STRCMP(argv[argc + count + 1], "0") == 0 ){
                            EB_STRCPY(configStr[count], COMMAND_LINE_MAX_SIZE, argv[argc + count + 1]);
                        }else if (argv[argc + count + 1][0] != '-'){
                            EB_STRCPY(configStr[count], COMMAND_LINE_MAX_SIZE, argv[argc + count + 1]);
                        }else {
                            EB_STRCPY(configStr[count], COMMAND_LINE_MAX_SIZE," ");
                            done = 1;
                        }
                    }else{
                        EB_STRCPY(configStr[count], COMMAND_LINE_MAX_SIZE, " ");
                        done =1;
                        //return return_error;
                    }
                }else
                    EB_STRCPY(configStr[count], COMMAND_LINE_MAX_SIZE, " ");
            }
        }
    }

    return return_error;
}

uint32_t get_help(
    int32_t     argc,
    char *const argv[])
{
    char config_string[COMMAND_LINE_MAX_SIZE];
    if (FindToken(argc, argv, HELP_TOKEN, config_string) == 0) {
        int32_t token_index = -1;

        printf("\n%-25s\t%-25s\t%s\n\n" ,"TOKEN", "DESCRIPTION", "INPUT TYPE");
        printf("%-25s\t%-25s\t%s\n" ,"-nch", "NumberOfChannels", "Single input");
        while (config_entry[++token_index].token != NULL)
            printf("%-25s\t%-25s\t%s\n", config_entry[token_index].token, config_entry[token_index].name, config_entry[token_index].type ? "Array input": "Single input");
        return 1;
    }
    else
        return 0;
}

/******************************************************
* Get the number of channels and validate it with input
******************************************************/
uint32_t get_number_of_channels(
    int32_t     argc,
    char *const argv[])
{
    char config_string[COMMAND_LINE_MAX_SIZE];
    uint32_t channelNumber;
    if (FindToken(argc, argv, CHANNEL_NUMBER_TOKEN, config_string) == 0) {
        // Set the input file
        channelNumber = strtol(config_string,  NULL, 0);
        if ((channelNumber > (uint32_t) MAX_CHANNEL_NUMBER) || channelNumber == 0){
            printf("Error: The number of channels has to be within the range [1,%u]\n",(uint32_t) MAX_CHANNEL_NUMBER);
            return 0;
        }else{
            return channelNumber;
        }
    }
    return 1;
}

void mark_token_as_read(
    const char  * token,
    char        * cmd_copy[],
    int32_t     * cmd_token_cnt
    )
{
    int32_t cmd_copy_index;
    for (cmd_copy_index = 0; cmd_copy_index < *(cmd_token_cnt); ++cmd_copy_index) {
        if (!EB_STRCMP(cmd_copy[cmd_copy_index], token))
            cmd_copy[cmd_copy_index] = cmd_copy[--(*cmd_token_cnt)];
    }
}

EbBool is_negative_number(
    const char* string) {
    int32_t length = (int32_t)strlen(string);
    int32_t index = 0;
    if (string[0] != '-') return EB_FALSE;
    for (index = 1; index < length; index++)
    {
        if (string[index] < '0' || string[index] > '9')
            return EB_FALSE;
    }
    return EB_TRUE;
}

// Computes the number of frames in the input file
int32_t ComputeFramesToBeEncoded(
    EbConfig   *config)
{
    uint64_t fileSize   = 0;
    int32_t  frame_count = 0;
    uint32_t frameSize;
    uint64_t currLoc;

    if (config->input_file) {
        currLoc = ftello(config->input_file); // get current fp location
        fseeko(config->input_file, 0L, SEEK_END);
        fileSize = ftello(config->input_file);
        fseeko(config->input_file, currLoc, SEEK_SET); // seek back to that location
    }

    frameSize = config->input_padded_width * config->input_padded_height; // Luma
    frameSize += 2 * (frameSize >> (3 - config->encoder_color_format)); // Add Chroma
    frameSize = frameSize << ((config->encoder_bit_depth == 10) ? 1 : 0);

    if (frameSize == 0)
        return -1;

    if (config->encoder_bit_depth == 10 && config->compressed_ten_bit_format == 1)
        frame_count = (int32_t)(2 * ((double)fileSize / frameSize) / 1.25);
    else
        frame_count = (int32_t)(fileSize / frameSize);

    if (frame_count == 0)
        return -1;

    return frame_count;
}

/******************************************
* Read Command Line
******************************************/
EbErrorType read_command_line(
    int32_t        argc,
    char *const    argv[],
    EbConfig   **configs,
    uint32_t       num_channels,
    EbErrorType  *return_errors)
{
    EbErrorType return_error = EB_ErrorBadParameter;
    char            config_string[COMMAND_LINE_MAX_SIZE];        // for one input options
    char           *config_strings[MAX_CHANNEL_NUMBER]; // for multiple input options
    char           *cmd_copy[MAX_NUM_TOKENS];                 // keep track of extra tokens
    uint32_t    index           = 0;
    int32_t             cmd_token_cnt   = 0;                        // total number of tokens
    int32_t             token_index     = -1;
    int32_t ret_y4m;

    for (index = 0; index < MAX_CHANNEL_NUMBER; ++index)
        config_strings[index] = (char*)malloc(sizeof(char)*COMMAND_LINE_MAX_SIZE);
    // Copy tokens (except for CHANNEL_NUMBER_TOKEN ) into a temp token buffer hosting all tokens that are passed through the command line
    size_t len = EB_STRLEN(CHANNEL_NUMBER_TOKEN, COMMAND_LINE_MAX_SIZE);
    for (token_index = 0; token_index < argc; ++token_index) {
        if ((argv[token_index][0] == '-') && strncmp(argv[token_index], CHANNEL_NUMBER_TOKEN, len) && !is_negative_number(argv[token_index]))
                cmd_copy[cmd_token_cnt++] = argv[token_index];
    }

    /***************************************************************************************************/
    /****************  Find configuration files tokens and call respective functions  ******************/
    /***************************************************************************************************/

    // Find the Config File Path in the command line
    if (FindTokenMultipleInputs(argc, argv, CONFIG_FILE_TOKEN, config_strings) == 0) {
        mark_token_as_read(CONFIG_FILE_TOKEN, cmd_copy, &cmd_token_cnt);
        // Parse the config file
        for (index = 0; index < num_channels; ++index){
            return_errors[index] = (EbErrorType)ReadConfigFile(configs[index], config_strings[index], index);
            return_error = (EbErrorType)(return_error &  return_errors[index]);
        }
    }
    else {
        if (FindToken(argc, argv, CONFIG_FILE_TOKEN, config_string) == 0) {
            printf("Error: Config File Token Not Found\n");
            return EB_ErrorBadParameter;
        }
        else
            return_error = EB_ErrorNone;
    }

    /***************************************************************************************************/
    /***********   Find SINGLE_INPUT configuration parameter tokens and call respective functions  **********/
    /***************************************************************************************************/
    token_index = -1;
    // Parse command line for tokens
    while (config_entry[++token_index].name != NULL){
        if (config_entry[token_index].type == SINGLE_INPUT){
            if (FindTokenMultipleInputs(argc, argv, config_entry[token_index].token, config_strings) == 0) {
                // When a token is found mark it as found in the temp token buffer
                mark_token_as_read(config_entry[token_index].token, cmd_copy, &cmd_token_cnt);

                // Fill up the values corresponding to each channel
                for (index = 0; index < num_channels; ++index) {
                    if (EB_STRCMP(config_strings[index], " "))
                        (*config_entry[token_index].scf)(config_strings[index], configs[index]);
                    else
                        break;
                }
            }
        }
    }

    /***************************************************************************************************/
    /********************** Parse parameters from input file if in y4m format **************************/
    /********************** overriding config file and command line inputs    **************************/
    /***************************************************************************************************/

    for (index = 0; index < num_channels; ++index) {
        if ((configs[index])->y4m_input == EB_TRUE){
            ret_y4m = read_y4m_header(configs[index]);
            if(ret_y4m == EB_ErrorBadParameter){
                printf("Error found when reading the y4m file parameters.\n");
                return EB_ErrorBadParameter;
            }
        }
    }

    /***************************************************************************************************/
    /***********   Find SPECIAL configuration parameter tokens and call respective functions  **********/
    /***************************************************************************************************/
    // Parse command line for search region at level 0 width token
    if (FindTokenMultipleInputs(argc, argv, HME_LEVEL0_WIDTH, config_strings) == 0) {
        uint32_t inputIndex = 0, lastIndex = 0;
        uint32_t done = 1;

        mark_token_as_read(HME_LEVEL0_WIDTH, cmd_copy, &cmd_token_cnt);

        for (index = 0; done && (index < num_channels); ++index){
            configs[index]->hme_level0_column_index = 0;
            for (inputIndex = lastIndex; inputIndex < configs[index]->number_hme_search_region_in_width + lastIndex; ++inputIndex){
                if (EB_STRCMP(config_strings[inputIndex], " "))
                    SetHmeLevel0SearchAreaInWidthArray(config_strings[inputIndex], configs[index]);
                else{
                    done = 0;
                    break;
                }
            }
            lastIndex += configs[index]->number_hme_search_region_in_width;
        }
    }

    //// Parse command line for search region at level 0 height token
    if (FindTokenMultipleInputs(argc, argv, HME_LEVEL0_HEIGHT, config_strings) == 0) {
        uint32_t inputIndex = 0, lastIndex = 0;
        uint32_t done = 1;

        mark_token_as_read(HME_LEVEL0_HEIGHT, cmd_copy, &cmd_token_cnt);

        for (index = 0; done && (index < num_channels); ++index){
            configs[index]->hme_level0_row_index = 0;
            for (inputIndex = lastIndex; inputIndex < configs[index]->number_hme_search_region_in_height + lastIndex; ++inputIndex){
                if (EB_STRCMP(config_strings[inputIndex], " "))
                    SetHmeLevel0SearchAreaInHeightArray(config_strings[inputIndex], configs[index]);
                else{
                    done = 0;
                    break;
                }
            }
            lastIndex += configs[index]->number_hme_search_region_in_height;
        }
    }

    // Parse command line for search region at level 1 Height token
    if (FindTokenMultipleInputs(argc, argv, HME_LEVEL1_HEIGHT, config_strings) == 0) {
        uint32_t inputIndex = 0, lastIndex = 0;
        uint32_t done = 1;

        mark_token_as_read(HME_LEVEL1_HEIGHT, cmd_copy, &cmd_token_cnt);

        for (index = 0; done && (index < num_channels); ++index){
            configs[index]->hme_level1_row_index = 0;
            for (inputIndex = lastIndex; inputIndex < configs[index]->number_hme_search_region_in_height + lastIndex; ++inputIndex){
                if (EB_STRCMP(config_strings[inputIndex], " "))
                    SetHmeLevel1SearchAreaInHeightArray(config_strings[inputIndex], configs[index]);
                else{
                    done = 0;
                    break;
                }
            }
            lastIndex += configs[index]->number_hme_search_region_in_height;
        }
    }

    // Parse command line for search region at level 1 width token
    if (FindTokenMultipleInputs(argc, argv, HME_LEVEL1_WIDTH, config_strings) == 0) {
        uint32_t inputIndex = 0, lastIndex = 0;
        uint32_t done = 1;

        mark_token_as_read(HME_LEVEL1_WIDTH, cmd_copy, &cmd_token_cnt);

        for (index = 0; done && (index < num_channels); ++index){
            configs[index]->hme_level1_column_index = 0;
            for (inputIndex = lastIndex; inputIndex < configs[index]->number_hme_search_region_in_width + lastIndex; ++inputIndex){
                if (EB_STRCMP(config_strings[inputIndex], " "))
                    SetHmeLevel1SearchAreaInWidthArray(config_strings[inputIndex], configs[index]);
                else{
                    done = 0;
                    break;
                }
            }
            lastIndex += configs[index]->number_hme_search_region_in_width;
        }
    }

    // Parse command line for search region at level 2 width token
    if (FindTokenMultipleInputs(argc, argv, HME_LEVEL2_WIDTH, config_strings) == 0) {
        uint32_t inputIndex = 0, lastIndex = 0;
        uint32_t done = 1;

        mark_token_as_read(HME_LEVEL2_WIDTH, cmd_copy, &cmd_token_cnt);

        for (index = 0; done && (index < num_channels); ++index){
            configs[index]->hme_level2_column_index = 0;
            for (inputIndex = lastIndex; inputIndex < configs[index]->number_hme_search_region_in_width + lastIndex; ++inputIndex){
                if (EB_STRCMP(config_strings[inputIndex], " "))
                    SetHmeLevel2SearchAreaInWidthArray(config_strings[inputIndex], configs[index]);
                else{
                    done = 0;
                    break;
                }
            }
            lastIndex += configs[index]->number_hme_search_region_in_width;
        }
    }

    // Parse command line for search region at level 2 height token
    if (FindTokenMultipleInputs(argc, argv, HME_LEVEL2_HEIGHT, config_strings) == 0) {
        uint32_t inputIndex = 0, lastIndex = 0;
        uint32_t done = 1;

        mark_token_as_read(HME_LEVEL2_HEIGHT, cmd_copy, &cmd_token_cnt);

        for (index = 0; done && (index < num_channels); ++index){
            configs[index]->hme_level2_row_index = 0;
            for (inputIndex = lastIndex; inputIndex < configs[index]->number_hme_search_region_in_height + lastIndex; ++inputIndex){
                if (EB_STRCMP(config_strings[inputIndex], " "))
                    SetHmeLevel2SearchAreaInHeightArray(config_strings[inputIndex], configs[index]);
                else{
                    done = 0;
                    break;
                }
            }
            lastIndex += configs[index]->number_hme_search_region_in_height;
        }
    }

    /***************************************************************************************************/
    /**************************************   Verify configuration parameters   ************************/
    /***************************************************************************************************/
    // Verify the config values
    if (return_error == 0) {
        return_error = EB_ErrorBadParameter;
        for (index = 0; index < num_channels; ++index){
            if (return_errors[index] == EB_ErrorNone){
                return_errors[index] = VerifySettings(configs[index], index);

                // Assuming no errors, add padding to width and height
                if (return_errors[index] == EB_ErrorNone) {
                    configs[index]->input_padded_width = configs[index]->source_width + configs[index]->source_width % 8;
                    configs[index]->input_padded_height = configs[index]->source_height + configs[index]->source_width % 8;
                }

                // Assuming no errors, set the frames to be encoded to the number of frames in the input yuv
                if (return_errors[index] == EB_ErrorNone && configs[index]->frames_to_be_encoded == 0)
                    configs[index]->frames_to_be_encoded = ComputeFramesToBeEncoded(configs[index]);

                if (configs[index]->frames_to_be_encoded == -1) {
                    fprintf(configs[index]->error_log_file, "Error instance %u: Input yuv does not contain enough frames \n", index + 1);
                    return_errors[index] = EB_ErrorBadParameter;
                }

                // Force the injector latency mode, and injector frame rate when speed control is on
                if (return_errors[index] == EB_ErrorNone && configs[index]->speed_control_flag == 1)
                    configs[index]->injector    = 1;
            }
            return_error = (EbErrorType)(return_error & return_errors[index]);
        }
    }

    // Print message for unprocessed tokens
    if (cmd_token_cnt > 0) {
        int32_t cmd_copy_index;
        printf("Unprocessed tokens: ");
        for (cmd_copy_index = 0; cmd_copy_index < cmd_token_cnt; ++cmd_copy_index)
            printf(" %s ", cmd_copy[cmd_copy_index]);
        printf("\n\n");
        return_error = EB_ErrorBadParameter;
    }

    for (index = 0; index < MAX_CHANNEL_NUMBER; ++index)
        free(config_strings[index]);
    return return_error;
}
