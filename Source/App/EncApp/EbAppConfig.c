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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "EbAppConfig.h"
#include "EbAppContext.h"
#include "EbAppInputy4m.h"
#ifdef _WIN32
#include <windows.h>
#include <io.h>
#else
#include <unistd.h>
#include <sys/file.h>
#endif

#if !defined(_WIN32) || !defined(HAVE_STRNLEN_S)
#include "safe_str_lib.h"
#endif

/**********************************
 * Defines
 **********************************/
#define HELP_TOKEN "-help"
#define HELP_LONG_TOKEN "--help"
#define CHANNEL_NUMBER_TOKEN "-nch"
#define COMMAND_LINE_MAX_SIZE 2048
#define CONFIG_FILE_TOKEN "-c"
#define CONFIG_FILE_LONG_TOKEN "--config"
#define INPUT_FILE_TOKEN "-i"
#define OUTPUT_BITSTREAM_TOKEN "-b"
#define OUTPUT_RECON_TOKEN "-o"
#define ERROR_FILE_TOKEN "-errlog"
#define QP_FILE_TOKEN "-qp-file"

/* two pass token */
#define PASS_TOKEN "--pass"
#define TWO_PASS_STATS_TOKEN "--stats"
#define PASSES_TOKEN "--passes"
#define INPUT_STAT_FILE_TOKEN "-input-stat-file"
#define OUTPUT_STAT_FILE_TOKEN "-output-stat-file"
#define STAT_FILE_TOKEN "-stat-file"
#define INPUT_PREDSTRUCT_FILE_TOKEN "-pred-struct-file"
#define WIDTH_TOKEN "-w"
#define HEIGHT_TOKEN "-h"
#define NUMBER_OF_PICTURES_TOKEN "-n"
#define BUFFERED_INPUT_TOKEN "-nb"
#define NO_PROGRESS_TOKEN "--no-progress" // tbd if it should be removed
#define PROGRESS_TOKEN "--progress"
#define BASE_LAYER_SWITCH_MODE_TOKEN "-base-layer-switch-mode" // no Eval
#define QP_TOKEN "-q"
#define USE_QP_FILE_TOKEN "-use-q-file"
#define STAT_REPORT_TOKEN "-stat-report"
#define FRAME_RATE_TOKEN "-fps"
#define FRAME_RATE_NUMERATOR_TOKEN "-fps-num"
#define FRAME_RATE_DENOMINATOR_TOKEN "-fps-denom"
#define ENCODER_BIT_DEPTH "-bit-depth"
#define ENCODER_16BIT_PIPELINE "-16bit-pipeline"
#define ENCODER_COLOR_FORMAT "-color-format"
#define INPUT_COMPRESSED_TEN_BIT_FORMAT "-compressed-ten-bit-format"
#define ENCMODE_TOKEN "-enc-mode"
#define HIERARCHICAL_LEVELS_TOKEN "-hierarchical-levels" // no Eval
#define PRED_STRUCT_TOKEN "-pred-struct"
#define INTRA_PERIOD_TOKEN "-intra-period"
#define PROFILE_TOKEN "-profile"
#define TIER_TOKEN "-tier"
#define LEVEL_TOKEN "-level"
#define FILM_GRAIN_TOKEN "-film-grain"
#define INTRA_REFRESH_TYPE_TOKEN "-irefresh-type" // no Eval
#define LOOP_FILTER_DISABLE_TOKEN "-dlf"
#define CDEF_LEVEL_TOKEN "-cdef-level"
#define RESTORATION_ENABLE_TOKEN "-restoration-filtering"
#define SG_FILTER_MODE_TOKEN "-sg-filter-mode"
#define WN_FILTER_MODE_TOKEN "-wn-filter-mode"
#define INTRA_ANGLE_DELTA_TOKEN "-intra-angle-delta"
#define INTER_INTRA_COMPOUND_TOKEN "-interintra-comp"
#define PAETH_TOKEN "-paeth"
#define SMOOTH_TOKEN "-smooth"
#define MFMV_ENABLE_TOKEN "-mfmv"
#define REDUNDANT_BLK_TOKEN "-redundant-blk"
#define SPATIAL_SSE_FL_TOKEN "-spatial-sse-full-loop-level"
#define OVR_BNDRY_BLK_TOKEN "-over-bndry-blk"
#define NEW_NEAREST_COMB_INJECT_TOKEN "-new-nrst-near-comb"
#define NSQ_TABLE_TOKEN "-nsq-table-use"
#define FRAME_END_CDF_UPDATE_TOKEN "-framend-cdf-upd-mode"
#define LOCAL_WARPED_ENABLE_TOKEN "-local-warp"
#define GLOBAL_MOTION_ENABLE_TOKEN "-global-motion"
#define OBMC_TOKEN "-obmc-level"
#define RDOQ_TOKEN "-rdoq-level"
#define PRED_ME_TOKEN "-pred-me"
#define BIPRED_3x3_TOKEN "-bipred-3x3"
#define COMPOUND_LEVEL_TOKEN "-compound"
#define FILTER_INTRA_TOKEN "-filter-intra-level"
#define INTRA_EDGE_FILTER_TOKEN "-intra-edge-filter"
#define PIC_BASED_RATE_EST_TOKEN "-pic-based-rate-est"
#define PIC_BASED_RATE_EST_NEW_TOKEN "--enable-pic-based-rate-est"
#define USE_DEFAULT_ME_HME_TOKEN "-use-default-me-hme"
#define HME_ENABLE_TOKEN "-hme"
#define HME_L0_ENABLE_TOKEN "-hme-l0"
#define HME_L1_ENABLE_TOKEN "-hme-l1"
#define HME_L2_ENABLE_TOKEN "-hme-l2"
#define EXT_BLOCK "-ext-block"
#define SEARCH_AREA_WIDTH_TOKEN "-search-w"
#define SEARCH_AREA_HEIGHT_TOKEN "-search-h"
#define SCREEN_CONTENT_TOKEN "-scm"
#define INTRABC_MODE_TOKEN "-intrabc-mode"
// --- start: ALTREF_FILTERING_SUPPORT
#define TF_LEVEL "-tf-level"
#define ALTREF_STRENGTH "-altref-strength"
#define ALTREF_NFRAMES "-altref-nframes"
#define ENABLE_OVERLAYS "-enable-overlays"
// --- end: ALTREF_FILTERING_SUPPORT
// --- start: SUPER-RESOLUTION SUPPORT
#define SUPERRES_MODE_INPUT "-superres-mode"
#define SUPERRES_DENOM "-superres-denom"
#define SUPERRES_KF_DENOM "-superres-kf-denom"
#define SUPERRES_QTHRES "-superres-qthres"
// --- end: SUPER-RESOLUTION SUPPORT
#define HBD_MD_ENABLE_TOKEN "-hbd-md"
#define PALETTE_TOKEN "-palette-level"
#define HDR_INPUT_TOKEN "-hdr"
#define RATE_CONTROL_ENABLE_TOKEN "-rc"
#define TARGET_BIT_RATE_TOKEN "-tbr"
#define MAX_QP_TOKEN "-max-qp"
#define VBV_BUFSIZE_TOKEN "-vbv-bufsize"
#define MIN_QP_TOKEN "-min-qp"
#define VBR_BIAS_PCT_TOKEN "-bias-pct"
#define VBR_MIN_SECTION_PCT_TOKEN "-minsection-pct"
#define VBR_MAX_SECTION_PCT_TOKEN "-maxsection-pct"
#define UNDER_SHOOT_PCT_TOKEN "-undershoot-pct"
#define OVER_SHOOT_PCT_TOKEN "-overshoot-pct"
#define ADAPTIVE_QP_ENABLE_TOKEN "-adaptive-quantization"
#define LOOK_AHEAD_DIST_TOKEN "-lad"
#define ENABLE_TPL_LA_TOKEN "-enable-tpl-la"
#define SUPER_BLOCK_SIZE_TOKEN "-sb-size"
#define TILE_ROW_TOKEN "-tile-rows"
#define TILE_COL_TOKEN "-tile-columns"

#define SQ_WEIGHT_TOKEN "-sqw"
#define CHROMA_MODE_TOKEN "-chroma-mode"
#define DISABLE_CFL_TOKEN "-dcfl"

#define SCENE_CHANGE_DETECTION_TOKEN "-scd"
#define INJECTOR_TOKEN "-inj" // no Eval
#define INJECTOR_FRAMERATE_TOKEN "-inj-frm-rt" // no Eval
#define SPEED_CONTROL_TOKEN "-speed-ctrl"
#define ASM_TYPE_TOKEN "-asm"
#define THREAD_MGMNT "-lp"
#define UNPIN_TOKEN "-unpin"
#define TARGET_SOCKET "-ss"
#define UNRESTRICTED_MOTION_VECTOR "-umv"
#define CONFIG_FILE_COMMENT_CHAR '#'
#define CONFIG_FILE_NEWLINE_CHAR '\n'
#define CONFIG_FILE_RETURN_CHAR '\r'
#define CONFIG_FILE_VALUE_SPLIT ':'
#define CONFIG_FILE_SPACE_CHAR ' '
#define CONFIG_FILE_ARRAY_SEP_CHAR CONFIG_FILE_SPACE_CHAR
#define CONFIG_FILE_TAB_CHAR '\t'
#define CONFIG_FILE_NULL_CHAR '\0'
#define CONFIG_FILE_MAX_ARG_COUNT 256
#define CONFIG_FILE_MAX_VAR_LEN 128
#define EVENT_FILE_MAX_ARG_COUNT 20
#define EVENT_FILE_MAX_VAR_LEN 256
#define BUFFER_FILE_MAX_ARG_COUNT 320
#define BUFFER_FILE_MAX_VAR_LEN 128

#define MDS_1_PRUNE_C_TH "-mds-1-class-th"
#define MDS_1_PRUNE_S_TH "-mds-1-cand-th"
#define MDS_2_3_PRUNE_C_TH "-mds-2-3-class-th"
#define MDS_2_3_PRUNE_S_TH "-mds-2-3-cand-th"
//double dash
#define PRESET_TOKEN "--preset"
#define QP_FILE_NEW_TOKEN "--qpfile"
#define INPUT_DEPTH_TOKEN "--input-depth"
#define KEYINT_TOKEN "--keyint"
#define LOOKAHEAD_NEW_TOKEN "--lookahead"

#define STAT_REPORT_NEW_TOKEN "--enable-stat-report"
#define RESTORATION_ENABLE_NEW_TOKEN "--enable-restoration-filtering"
#define INTER_INTRA_COMPOUND_NEW_TOKEN "--enable-interintra-comp"
#define FRAC_SEARCH_64_NEW_TOKEN "--enable-frac-search-64"
#define MFMV_ENABLE_NEW_TOKEN "--enable-mfmv"
#define REDUNDANT_BLK_NEW_TOKEN "--enable-redundant-blk"
#define SPATIAL_SSE_FL_NEW_TOKEN "--enable-spatial-sse-full-loop-level"
#define OVR_BNDRY_BLK_NEW_TOKEN "--enable-over-bndry-blk"
#define NEW_NEAREST_COMB_INJECT_NEW_TOKEN "--enable-new-nrst-near-comb"
#define NX4_4XN_MV_INJECT_NEW_TOKEN "--enable-nx4-4xn-mv-inject"
#define NSQ_TABLE_NEW_TOKEN "--enable-nsq-table-use"
#define FRAME_END_CDF_UPDATE_NEW_TOKEN "--enable-framend-cdf-upd-mode"
#define LOCAL_WARPED_ENABLE_NEW_TOKEN "--enable-local-warp"
#define GLOBAL_MOTION_ENABLE_NEW_TOKEN "--enable-global-motion"
#define RDOQ_NEW_TOKEN "--rdoq-level"
#define FILTER_INTRA_NEW_TOKEN "--filter-intra-level"
#define HDR_INPUT_NEW_TOKEN "--enable-hdr"
#define ADAPTIVE_QP_ENABLE_NEW_TOKEN "--aq-mode"
#define INPUT_FILE_LONG_TOKEN "--input"
#define OUTPUT_BITSTREAM_LONG_TOKEN "--output"
#define OUTPUT_RECON_LONG_TOKEN "--recon"
#define WIDTH_LONG_TOKEN "--width"
#define HEIGHT_LONG_TOKEN "--height"
#define NUMBER_OF_PICTURES_LONG_TOKEN "--frames"
#define QP_LONG_TOKEN "--qp"
#define CLASS_12_NEW_TOKEN "--enable-class-12"
#define LOOP_FILTER_DISABLE_NEW_TOKEN "--disable-dlf"
#define DISABLE_CFL_NEW_TOKEN "--disable-cfl"
#define INTRA_EDGE_FILTER_NEW_TOKEN "--enable-intra-edge-filter"
#define INTRA_ANGLE_DELTA_NEW_TOKEN "--enable-intra-angle-delta"
#define PAETH_NEW_TOKEN "--enable-paeth"
#define SMOOTH_NEW_TOKEN "--enable-smooth"
#define MRP_LEVEL_TOKEN "--mrp-level"

#ifdef _WIN32
static HANDLE get_file_handle(FILE* fp)
{
    return (HANDLE)_get_osfhandle(_fileno(fp));
}
#endif

static EbBool fopen_and_lock(FILE** file, const char* name, EbBool write)
{
    if (!file || !name)
        return EB_FALSE;

    const char* mode = write ? "wb" : "rb";
    FOPEN(*file, name, mode);
    if (!*file)
        return EB_FALSE;

#ifdef _WIN32
    HANDLE handle = get_file_handle(*file);
    if (handle == INVALID_HANDLE_VALUE) return EB_FALSE;
    if (LockFile(handle, 0, 0, MAXDWORD, MAXDWORD))
        return EB_TRUE;
#else
    int fd = fileno(*file);
    if (flock(fd, LOCK_EX | LOCK_NB) == 0)
        return EB_TRUE;
#endif
    fprintf(stderr, "ERROR: locking %s failed, is it used by other encoder?\n", name);
    return EB_FALSE;
}

/**********************************
 * Set Cfg Functions
 **********************************/
static void set_cfg_input_file(const char *filename, EbConfig *cfg) {
    if (cfg->input_file && !cfg->input_file_is_fifo) fclose(cfg->input_file);

    if (!filename) {
        cfg->input_file = NULL;
        return;
    }

    if (!strcmp(filename, "stdin"))
        cfg->input_file = stdin;
    else
        FOPEN(cfg->input_file, filename, "rb");

    if (cfg->input_file == NULL) { return; }

#ifdef _WIN32
    HANDLE handle = (HANDLE)_get_osfhandle(_fileno(cfg->input_file));
    if (handle == INVALID_HANDLE_VALUE) return;
    cfg->input_file_is_fifo = GetFileType(handle) == FILE_TYPE_PIPE;
#else
    int         fd = fileno(cfg->input_file);
    struct stat statbuf;
    fstat(fd, &statbuf);
    cfg->input_file_is_fifo = S_ISFIFO(statbuf.st_mode);
#endif

    cfg->y4m_input = check_if_y4m(cfg);
};

static void set_pred_struct_file(const char *value, EbConfig *cfg) {

    if (cfg->input_pred_struct_filename) { free(cfg->input_pred_struct_filename); }
    cfg->input_pred_struct_filename = (char *)malloc(strlen(value) + 1);
    strcpy_s(cfg->input_pred_struct_filename, strlen(value) + 1, value);

    cfg->enable_manual_pred_struct = EB_TRUE;
};

static void set_cfg_stream_file(const char *value, EbConfig *cfg) {
    if (cfg->bitstream_file && cfg->bitstream_file != stdout) { fclose(cfg->bitstream_file); }

    if (!strcmp(value, "stdout")) {
        cfg->bitstream_file = stdout;
    } else {
        FOPEN(cfg->bitstream_file, value, "wb");
    }
};
static void set_cfg_error_file(const char *value, EbConfig *cfg) {
    if (cfg->error_log_file && cfg->error_log_file != stderr) { fclose(cfg->error_log_file); }
    FOPEN(cfg->error_log_file, value, "w+");
};
static void set_cfg_recon_file(const char *value, EbConfig *cfg) {
    if (cfg->recon_file) { fclose(cfg->recon_file); }
    FOPEN(cfg->recon_file, value, "wb");
};
static void set_cfg_qp_file(const char *value, EbConfig *cfg) {
    if (cfg->qp_file) { fclose(cfg->qp_file); }
    FOPEN(cfg->qp_file, value, "r");
};

static void set_pass(const char* value, EbConfig *cfg) {
    cfg->pass = strtol(value, NULL, 0);
}

static void set_two_pass_stats(const char* value, EbConfig *cfg) {
#ifndef _WIN32
    cfg->stats = strdup(value);
#else
    cfg->stats = _strdup(value);
#endif
}

static void set_passes(const char* value, EbConfig *cfg) {
    (void)value;
    (void)cfg;
    /* empty function, we will handle passes at higher level*/
    return;
}

static void set_cfg_stat_file(const char *value, EbConfig *cfg) {
    if (cfg->stat_file) { fclose(cfg->stat_file); }
    FOPEN(cfg->stat_file, value, "wb");
};
static void set_stat_report(const char *value, EbConfig *cfg) {
    cfg->stat_report = (uint8_t)strtoul(value, NULL, 0);
};
static void set_cfg_source_width(const char *value, EbConfig *cfg) {
    cfg->source_width = strtoul(value, NULL, 0);
};
static void set_cfg_source_height(const char *value, EbConfig *cfg) {
    cfg->source_height = strtoul(value, NULL, 0);
};
static void set_cfg_frames_to_be_encoded(const char *value, EbConfig *cfg) {
    cfg->frames_to_be_encoded = strtol(value, NULL, 0);
};
static void set_buffered_input(const char *value, EbConfig *cfg) {
    cfg->buffered_input = strtol(value, NULL, 0);
};
static void set_no_progress(const char *value, EbConfig *cfg) {
    switch (value ? *value : '1') {
    case '0': cfg->progress = 1; break; // equal to --progress 1
    default: cfg->progress = 0; break; // equal to --progress 0
    }
}
static void set_progress(const char *value, EbConfig *cfg) {
    switch (value ? *value : '1') {
    case '0': cfg->progress = 0; break; // no progress printed
    case '2': cfg->progress = 2; break; // aomenc style progress
    default: cfg->progress = 1; break; // default progress
    }
}
static void set_frame_rate(const char *value, EbConfig *cfg) {
    cfg->frame_rate = strtoul(value, NULL, 0);
    if (cfg->frame_rate <= 1000)
        cfg->frame_rate <<= 16;
}

static void set_frame_rate_numerator(const char *value, EbConfig *cfg) {
    cfg->frame_rate_numerator = strtoul(value, NULL, 0);
};
static void set_frame_rate_denominator(const char *value, EbConfig *cfg) {
    cfg->frame_rate_denominator = strtoul(value, NULL, 0);
};
static void set_encoder_bit_depth(const char *value, EbConfig *cfg) {
    cfg->encoder_bit_depth = strtoul(value, NULL, 0);
}
static void set_encoder_16bit_pipeline(const char *value, EbConfig *cfg) {
    cfg->is_16bit_pipeline = (EbBool)strtoul(value, NULL, 0);
}
static void set_encoder_color_format(const char *value, EbConfig *cfg) {
    cfg->encoder_color_format = strtoul(value, NULL, 0);
}
static void set_compressed_ten_bit_format(const char *value, EbConfig *cfg) {
    cfg->compressed_ten_bit_format = strtoul(value, NULL, 0);
}
static void set_enc_mode(const char *value, EbConfig *cfg) {
    cfg->enc_mode = (uint8_t)strtoul(value, NULL, 0);
};
static void set_cfg_intra_period(const char *value, EbConfig *cfg) {
    cfg->intra_period = strtol(value, NULL, 0);
};
static void set_cfg_intra_refresh_type(const char *value, EbConfig *cfg) {
    cfg->intra_refresh_type = strtol(value, NULL, 0);
};
static void set_hierarchical_levels(const char *value, EbConfig *cfg) {
    cfg->hierarchical_levels = strtol(value, NULL, 0);
};
static void set_cfg_pred_structure(const char *value, EbConfig *cfg) {
    cfg->pred_structure = strtol(value, NULL, 0);
};
static void set_cfg_qp(const char *value, EbConfig *cfg) { cfg->qp = strtoul(value, NULL, 0); };
static void set_cfg_use_qp_file(const char *value, EbConfig *cfg) {
    cfg->use_qp_file = (EbBool)strtol(value, NULL, 0);
};
static void set_cfg_film_grain(const char *value, EbConfig *cfg) {
    cfg->film_grain_denoise_strength = strtol(value, NULL, 0);
}; //not bool to enable possible algorithm extension in the future
static void set_disable_dlf_flag(const char *value, EbConfig *cfg) {
    cfg->disable_dlf_flag = (EbBool)strtoul(value, NULL, 0);
};
static void set_enable_local_warped_motion_flag(const char *value, EbConfig *cfg) {
    cfg->enable_warped_motion = strtol(value, NULL, 0);
};
static void set_enable_global_motion_flag(const char *value, EbConfig *cfg) {
    cfg->enable_global_motion = (EbBool)strtoul(value, NULL, 0);
};
static void set_cdef_level(const char *value, EbConfig *cfg) {
    cfg->cdef_level = strtol(value, NULL, 0);
};
static void set_enable_restoration_filter_flag(const char *value, EbConfig *cfg) {
    cfg->enable_restoration_filtering = strtol(value, NULL, 0);
};
static void set_sg_filter_mode(const char *value, EbConfig *cfg) {
    cfg->sg_filter_mode = strtol(value, NULL, 0);
};
static void set_wn_filter_mode(const char *value, EbConfig *cfg) {
    cfg->wn_filter_mode = strtol(value, NULL, 0);
};
static void set_intra_angle_delta_flag(const char *value, EbConfig *cfg) {
    cfg->intra_angle_delta = strtol(value, NULL, 0);
};
static void set_interintra_compound_flag(const char *value, EbConfig *cfg) {
    cfg->inter_intra_compound = strtol(value, NULL, 0);
};
static void set_enable_paeth_flag(const char *value, EbConfig *cfg) {
    cfg->enable_paeth = strtol(value, NULL, 0);
};
static void set_enable_smooth_flag(const char *value, EbConfig *cfg) {
    cfg->enable_smooth = strtol(value, NULL, 0);
};
static void set_mrp_level(const char *value, EbConfig *cfg) {
    cfg->mrp_level = strtol(value, NULL, 0);
};
static void set_enable_mfmv_flag(const char *value, EbConfig *cfg) {
    cfg->enable_mfmv = strtol(value, NULL, 0);
};
static void set_enable_redundant_blk_flag(const char *value, EbConfig *cfg) {
    cfg->enable_redundant_blk = strtol(value, NULL, 0);
};
static void set_spatial_sse_full_loop_level_flag(const char *value, EbConfig *cfg) {
    cfg->spatial_sse_full_loop_level = strtol(value, NULL, 0);
};
static void set_over_bndry_blk_flag(const char *value, EbConfig *cfg) {
    cfg->over_bndry_blk = strtol(value, NULL, 0);
};
static void set_new_nearest_comb_inject_flag(const char *value, EbConfig *cfg) {
    cfg->new_nearest_comb_inject = strtol(value, NULL, 0);
};
static void set_nsq_table_flag(const char *value, EbConfig *cfg) {
    cfg->nsq_table = strtol(value, NULL, 0);
};
static void set_frame_end_cdf_update_flag(const char *value, EbConfig *cfg) {
    cfg->frame_end_cdf_update = strtol(value, NULL, 0);
};
static void set_chroma_mode(const char *value, EbConfig *cfg) {
    cfg->set_chroma_mode = strtol(value, NULL, 0);
};
static void set_disable_cfl_flag(const char *value, EbConfig *cfg) {
    cfg->disable_cfl_flag = strtol(value, NULL, 0);
};
static void set_obmc_level_flag(const char *value, EbConfig *cfg) {
    cfg->obmc_level = (EbBool)strtoul(value, NULL, 0);
};
static void set_rdoq_level_flag(const char *value, EbConfig *cfg) {
    cfg->rdoq_level = strtol(value, NULL, 0);
};
static void set_predictive_me_flag(const char *value, EbConfig *cfg) {
    cfg->pred_me = strtol(value, NULL, 0);
};
static void set_bipred3x3inject_flag(const char *value, EbConfig *cfg) {
    cfg->bipred_3x3_inject = strtol(value, NULL, 0);
};
static void set_compound_level_flag(const char *value, EbConfig *cfg) {
    cfg->compound_level = strtol(value, NULL, 0);
};
static void set_filter_intra_level_flag(const char *value, EbConfig *cfg) {
    cfg->filter_intra_level = (int8_t)strtoul(value, NULL, 0);
};
static void set_enable_intra_edge_filter_flag(const char *value, EbConfig *cfg) {
    cfg->enable_intra_edge_filter = strtol(value, NULL, 0);
};
static void set_pic_based_rate_est(const char *value, EbConfig *cfg) {
    cfg->pic_based_rate_est = strtol(value, NULL, 0);
};
static void set_enable_hme_flag(const char *value, EbConfig *cfg) {
    cfg->enable_hme_flag = (EbBool)strtoul(value, NULL, 0);
};
static void set_enable_hme_level_0_flag(const char *value, EbConfig *cfg) {
    cfg->enable_hme_level0_flag = (EbBool)strtoul(value, NULL, 0);
};
static void set_tile_row(const char *value, EbConfig *cfg) {
    cfg->tile_rows = strtoul(value, NULL, 0);
};
static void set_tile_col(const char *value, EbConfig *cfg) {
    cfg->tile_columns = strtoul(value, NULL, 0);
};
static void set_scene_change_detection(const char *value, EbConfig *cfg) {
    cfg->scene_change_detection = strtoul(value, NULL, 0);
}
static void set_look_ahead_distance(const char *value, EbConfig *cfg) {
    cfg->look_ahead_distance = strtoul(value, NULL, 0);
};
static void set_enable_tpl_la(const char *value, EbConfig *cfg) {
    cfg->enable_tpl_la = strtoul(value, NULL, 0);
};
static void set_rate_control_mode(const char *value, EbConfig *cfg) {
    cfg->rate_control_mode = strtoul(value, NULL, 0);
};
static void set_target_bit_rate(const char *value, EbConfig *cfg) {
    cfg->target_bit_rate = 1000 * strtoul(value, NULL, 0);
};
static void set_vbv_buf_size(const char *value, EbConfig *cfg) {
    cfg->vbv_bufsize = 1000 * strtoul(value, NULL, 0);
};
static void set_max_qp_allowed(const char *value, EbConfig *cfg) {
    cfg->max_qp_allowed = strtoul(value, NULL, 0);
};
static void set_min_qp_allowed(const char *value, EbConfig *cfg) {
    cfg->min_qp_allowed = strtoul(value, NULL, 0);
};
static void set_vbr_bias_pct(const char *value, EbConfig *cfg) {
    cfg->vbr_bias_pct = strtoul(value, NULL, 0);
};
static void set_vbr_min_section_pct(const char *value, EbConfig *cfg) {
    cfg->vbr_min_section_pct = strtoul(value, NULL, 0);
};
static void set_vbr_max_section_pct(const char *value, EbConfig *cfg) {
    cfg->vbr_max_section_pct = strtoul(value, NULL, 0);
};
static void set_under_shoot_pct(const char *value, EbConfig *cfg) {
    cfg->under_shoot_pct = strtoul(value, NULL, 0);
};
static void set_over_shoot_pct(const char *value, EbConfig *cfg) {
    cfg->over_shoot_pct = strtoul(value, NULL, 0);
};
static void set_adaptive_quantization(const char *value, EbConfig *cfg) {
    cfg->enable_adaptive_quantization = (EbBool)strtol(value, NULL, 0);
};
static void set_enable_hme_level_1_flag(const char *value, EbConfig *cfg) {
    cfg->enable_hme_level1_flag = (EbBool)strtoul(value, NULL, 0);
};
static void set_enable_hme_level_2_flag(const char *value, EbConfig *cfg) {
    cfg->enable_hme_level2_flag = (EbBool)strtoul(value, NULL, 0);
};
static void set_cfg_search_area_width(const char *value, EbConfig *cfg) {
    cfg->search_area_width = strtoul(value, NULL, 0);
};
static void set_cfg_search_area_height(const char *value, EbConfig *cfg) {
    cfg->search_area_height = strtoul(value, NULL, 0);
};
static void set_cfg_use_default_me_hme(const char *value, EbConfig *cfg) {
    cfg->use_default_me_hme = (EbBool)strtol(value, NULL, 0);
};
static void set_enable_ext_block_flag(const char *value, EbConfig *cfg) {
    cfg->ext_block_flag = (EbBool)strtoul(value, NULL, 0);
};
static void set_screen_content_mode(const char *value, EbConfig *cfg) {
    cfg->screen_content_mode = strtoul(value, NULL, 0);
};
static void set_intrabc_mode(const char *value, EbConfig *cfg) {
    cfg->intrabc_mode = strtol(value, NULL, 0);
};
// --- start: ALTREF_FILTERING_SUPPORT
static void set_tf_level(const char *value, EbConfig *cfg) {
    cfg->tf_level = (int8_t)strtoul(value, NULL, 0);
};
static void set_altref_strength(const char *value, EbConfig *cfg) {
    cfg->altref_strength = (uint8_t)strtoul(value, NULL, 0);
};
static void set_altref_n_frames(const char *value, EbConfig *cfg) {
    cfg->altref_nframes = (uint8_t)strtoul(value, NULL, 0);
};
static void set_enable_overlays(const char *value, EbConfig *cfg) {
    cfg->enable_overlays = (EbBool)strtoul(value, NULL, 0);
};
// --- end: ALTREF_FILTERING_SUPPORT
// --- start: SUPER-RESOLUTION SUPPORT
static void set_superres_mode(const char *value, EbConfig *cfg) {
    cfg->superres_mode = (SUPERRES_MODE)strtoul(value, NULL, 0);
};
static void set_superres_denom(const char *value, EbConfig *cfg) {
    cfg->superres_denom = (uint8_t)strtoul(value, NULL, 0);
};
static void set_superres_kf_denom(const char *value, EbConfig *cfg) {
    cfg->superres_kf_denom = (uint8_t)strtoul(value, NULL, 0);
};
static void set_superres_qthres(const char *value, EbConfig *cfg) {
    cfg->superres_qthres = (uint8_t)strtoul(value, NULL, 0);
};
// --- end: SUPER-RESOLUTION SUPPORT
static void set_enable_hbd_mode_decision(const char *value, EbConfig *cfg) {
    cfg->enable_hbd_mode_decision = (uint8_t)strtoul(value, NULL, 0);
};
static void set_palette_level(const char *value, EbConfig *cfg) {
    cfg->palette_level = (int32_t)strtoul(value, NULL, 0);
};
static void set_high_dynamic_range_input(const char *value, EbConfig *cfg) {
    cfg->high_dynamic_range_input = strtol(value, NULL, 0);
};
static void set_profile(const char *value, EbConfig *cfg) {
    cfg->profile = strtol(value, NULL, 0);
};
static void set_tier(const char *value, EbConfig *cfg) { cfg->tier = strtol(value, NULL, 0); };
static void set_level(const char *value, EbConfig *cfg) {
    if (strtoul(value, NULL, 0) != 0 || strcmp(value, "0") == 0)
        cfg->level = (uint32_t)(10 * strtod(value, NULL));
    else
        cfg->level = 9999999;
};
static void set_injector(const char *value, EbConfig *cfg) {
    cfg->injector = strtol(value, NULL, 0);
};
static void speed_control_flag(const char *value, EbConfig *cfg) {
    cfg->speed_control_flag = strtol(value, NULL, 0);
};
static void set_injector_frame_rate(const char *value, EbConfig *cfg) {
    cfg->injector_frame_rate = strtoul(value, NULL, 0);
    if (cfg->injector_frame_rate <= 1000)
        cfg->injector_frame_rate <<= 16;
}
static void set_asm_type(const char *value, EbConfig *cfg) {
    const struct {
        const char *name;
        CPU_FLAGS   flags;
    } param_maps[] = {
        {"c", 0},
        {"0", 0},
        {"mmx", (CPU_FLAGS_MMX << 1) - 1},
        {"1", (CPU_FLAGS_MMX << 1) - 1},
        {"sse", (CPU_FLAGS_SSE << 1) - 1},
        {"2", (CPU_FLAGS_SSE << 1) - 1},
        {"sse2", (CPU_FLAGS_SSE2 << 1) - 1},
        {"3", (CPU_FLAGS_SSE2 << 1) - 1},
        {"sse3", (CPU_FLAGS_SSE3 << 1) - 1},
        {"4", (CPU_FLAGS_SSE3 << 1) - 1},
        {"ssse3", (CPU_FLAGS_SSSE3 << 1) - 1},
        {"5", (CPU_FLAGS_SSSE3 << 1) - 1},
        {"sse4_1", (CPU_FLAGS_SSE4_1 << 1) - 1},
        {"6", (CPU_FLAGS_SSE4_1 << 1) - 1},
        {"sse4_2", (CPU_FLAGS_SSE4_2 << 1) - 1},
        {"7", (CPU_FLAGS_SSE4_2 << 1) - 1},
        {"avx", (CPU_FLAGS_AVX << 1) - 1},
        {"8", (CPU_FLAGS_AVX << 1) - 1},
        {"avx2", (CPU_FLAGS_AVX2 << 1) - 1},
        {"9", (CPU_FLAGS_AVX2 << 1) - 1},
        {"avx512", (CPU_FLAGS_AVX512VL << 1) - 1},
        {"10", (CPU_FLAGS_AVX512VL << 1) - 1},
        {"max", CPU_FLAGS_ALL},
        {"11", CPU_FLAGS_ALL},
    };
    const uint32_t para_map_size = sizeof(param_maps) / sizeof(param_maps[0]);
    uint32_t       i;

    for (i = 0; i < para_map_size; ++i) {
        if (strcmp(value, param_maps[i].name) == 0) {
            cfg->cpu_flags_limit = param_maps[i].flags;
            return;
        }
    }

    cfg->cpu_flags_limit = CPU_FLAGS_INVALID;
};
static void set_logical_processors(const char *value, EbConfig *cfg) {
    cfg->logical_processors = (uint32_t)strtoul(value, NULL, 0);
};
static void set_unpin_execution(const char *value, EbConfig *cfg) {
    cfg->unpin = (uint32_t)strtoul(value, NULL, 0);
};
static void set_target_socket(const char *value, EbConfig *cfg) {
    cfg->target_socket = (int32_t)strtol(value, NULL, 0);
};
static void set_unrestricted_motion_vector(const char *value, EbConfig *cfg) {
    cfg->unrestricted_motion_vector = (EbBool)strtol(value, NULL, 0);
};

enum CfgType {
    SINGLE_INPUT, // Configuration parameters that have only 1 value input
    ARRAY_INPUT // Configuration parameters that have multiple values as input
};

/**********************************
 * Config Entry Struct
 **********************************/
typedef struct config_entry_s {
    enum CfgType type;
    const char * token;
    const char * name;
    void (*scf)(const char *, EbConfig *);
} ConfigEntry;

/**********************************
 * Config Entry Array
 **********************************/
ConfigEntry config_entry_options[] = {
    // File I/O
    {SINGLE_INPUT, HELP_TOKEN, "Show usage options and exit", set_cfg_input_file},
    {SINGLE_INPUT, INPUT_FILE_TOKEN, "Input filename", set_cfg_input_file},
    {SINGLE_INPUT, INPUT_FILE_LONG_TOKEN, "Input filename", set_cfg_input_file},

    {SINGLE_INPUT, OUTPUT_BITSTREAM_TOKEN, "Output filename", set_cfg_stream_file},
    {SINGLE_INPUT, OUTPUT_BITSTREAM_LONG_TOKEN, "Output filename", set_cfg_stream_file},

    {SINGLE_INPUT, ERROR_FILE_TOKEN, "Error filename", set_cfg_error_file},
    {SINGLE_INPUT, OUTPUT_RECON_TOKEN, "Recon filename", set_cfg_recon_file},
    {SINGLE_INPUT, OUTPUT_RECON_LONG_TOKEN, "Recon filename", set_cfg_recon_file},

    {SINGLE_INPUT, STAT_FILE_TOKEN, "Stat filename", set_cfg_stat_file},
    {SINGLE_INPUT, NULL, NULL, NULL}};

ConfigEntry config_entry_global_options[] = {
    // Picture Dimensions
    {SINGLE_INPUT, WIDTH_TOKEN, "Frame width", set_cfg_source_width},
    {SINGLE_INPUT, WIDTH_LONG_TOKEN, "Frame width", set_cfg_source_width},

    {SINGLE_INPUT, HEIGHT_TOKEN, "Frame height", set_cfg_source_height},
    {SINGLE_INPUT, HEIGHT_LONG_TOKEN, "Frame height", set_cfg_source_height},

    {SINGLE_INPUT,
     NUMBER_OF_PICTURES_TOKEN,
     "Stop encoding after n input frames",
     set_cfg_frames_to_be_encoded},
    {SINGLE_INPUT,
     NUMBER_OF_PICTURES_LONG_TOKEN,
     "Stop encoding after n input frames",
     set_cfg_frames_to_be_encoded},

    {SINGLE_INPUT, BUFFERED_INPUT_TOKEN, "Buffer n input frames", set_buffered_input},
    {SINGLE_INPUT,
     PROGRESS_TOKEN,
     "Change verbosity of the output (0: no progress is printed, 1: default, 2: aomenc style "
     "machine parsable output)",
     set_progress},
    {SINGLE_INPUT,
     NO_PROGRESS_TOKEN,
     "Do not print out progress, if set to 1 it is equivalent to `" PROGRESS_TOKEN
     " 0`, else  `" PROGRESS_TOKEN " 1`",
     set_no_progress},
    {SINGLE_INPUT,
     ENCODER_COLOR_FORMAT,
     "Set encoder color format(EB_YUV400, EB_YUV420, EB_YUV422, EB_YUV444)",
     set_encoder_color_format},
    {SINGLE_INPUT,
     PROFILE_TOKEN,
     "Bitstream profile number to use(0: main profile[default], 1: high profile, 2: professional "
     "profile) ",
     set_profile},
    {SINGLE_INPUT, FRAME_RATE_TOKEN, "Stream frame rate (rate/scale)", set_frame_rate},
    {SINGLE_INPUT,
     FRAME_RATE_NUMERATOR_TOKEN,
     "Stream frame rate numerator",
     set_frame_rate_numerator},
    {SINGLE_INPUT,
     FRAME_RATE_DENOMINATOR_TOKEN,
     "Stream frame rate denominator",
     set_frame_rate_denominator},
    //{SINGLE_INPUT, ENCODER_BIT_DEPTH, "Bit depth for codec(8 or 10)", set_encoder_bit_depth},
    {SINGLE_INPUT, INPUT_DEPTH_TOKEN, "Bit depth for codec(8 or 10)", set_encoder_bit_depth},
    {SINGLE_INPUT,
     ENCODER_16BIT_PIPELINE,
     "Bit depth for enc-dec(0: lbd[default], 1: hbd)",
     set_encoder_16bit_pipeline},
    //{SINGLE_INPUT, LEVEL_TOKEN, "Level", set_level},
    {SINGLE_INPUT,
     HIERARCHICAL_LEVELS_TOKEN,
     "Set hierarchical levels(3 or 4[default])",
     set_hierarchical_levels},
    {SINGLE_INPUT,
     PRED_STRUCT_TOKEN,
     "Set prediction structure( 0: low delay P, 1: low delay B, 2: random access [default])",
     set_cfg_pred_structure},
    //{SINGLE_INPUT,
    // HDR_INPUT_TOKEN,
    // "Enable high dynamic range(0: OFF[default], ON: 1)",
    // set_high_dynamic_range_input},
    {SINGLE_INPUT,
     HDR_INPUT_NEW_TOKEN,
     "Enable high dynamic range(0: OFF[default], ON: 1)",
     set_high_dynamic_range_input},
    // Asm Type
    {SINGLE_INPUT,
     ASM_TYPE_TOKEN,
     "Limit assembly instruction set [0 - 11] or [c, mmx, sse, sse2, sse3, ssse3, sse4_1, sse4_2,"
     " avx, avx2, avx512, max], by default highest level supported by CPU",
     set_asm_type},
    {SINGLE_INPUT, THREAD_MGMNT, "number of logical processors to be used", set_logical_processors},
    {SINGLE_INPUT,
     UNPIN_TOKEN,
     "Allows the execution to be pined/unpined to/from a specific number of cores \n"
     "The combinational use of --unpin with --lp results in memory reduction while allowing the "
     "execution to work on any of the cores and not restrict it to specific cores \n"
     "--unpin is overwritten to 0 when --ss is set to 0 or 1. ( 0: OFF ,1: ON [default]) \n"
     "Example: 72 core machine: \n"
     "72 jobs x -- lp 1 -- unpin 1 \n"
     "36 jobs x -- lp 2 -- unpin 1 \n"
     "18 jobs x -- lp 4 -- unpin 1 ",
     set_unpin_execution},
    {SINGLE_INPUT,
     TARGET_SOCKET,
     "Specify  which socket the encoder runs on"
     "--unpin is overwritten to 0 when --ss is set to 0 or 1",
     set_target_socket},
    // Termination
    {SINGLE_INPUT, NULL, NULL, NULL}};

ConfigEntry config_entry_rc[] = {
    // Rate Control
    {SINGLE_INPUT,
     RATE_CONTROL_ENABLE_TOKEN,
     "Rate control mode(0 = CQP , 1 = VBR , 2 = CVBR)",
     set_rate_control_mode},
    {SINGLE_INPUT, TARGET_BIT_RATE_TOKEN, "Target Bitrate (kbps)", set_target_bit_rate},
    {SINGLE_INPUT,
     USE_QP_FILE_TOKEN,
     "Overwrite QP assignment using qp values in QP file",
     set_cfg_use_qp_file},
    //{SINGLE_INPUT, QP_FILE_TOKEN, "Path to Qp file", set_cfg_qp_file},
    {SINGLE_INPUT, QP_FILE_NEW_TOKEN, "Path to Qp file", set_cfg_qp_file},
    {SINGLE_INPUT, MAX_QP_TOKEN, "Maximum (worst) quantizer[0-63]", set_max_qp_allowed},
    {SINGLE_INPUT, MIN_QP_TOKEN, "Minimum (best) quantizer[0-63]", set_min_qp_allowed},
    //{SINGLE_INPUT,
    // ADAPTIVE_QP_ENABLE_TOKEN,
    // "Set adaptive QP level(0: OFF ,1: variance base using segments ,2: Deltaq pred efficiency)",
    // set_adaptive_quantization},
    {SINGLE_INPUT,
     ADAPTIVE_QP_ENABLE_TOKEN,
     "Set adaptive QP level(0: OFF ,1: variance base using segments ,2: Deltaq pred efficiency)",
     set_adaptive_quantization},
    {SINGLE_INPUT, VBV_BUFSIZE_TOKEN, "VBV buffer size", set_vbv_buf_size},
    {SINGLE_INPUT, UNDER_SHOOT_PCT_TOKEN, "Datarate undershoot (min) target (%)", set_under_shoot_pct},
    {SINGLE_INPUT, OVER_SHOOT_PCT_TOKEN, "Datarate overshoot (max) target (%)", set_over_shoot_pct},
    // Termination
    {SINGLE_INPUT, NULL, NULL, NULL}};
ConfigEntry config_entry_2p[] = {
    // 2 pass
    {SINGLE_INPUT, PASS_TOKEN, "Multipass bitrate control (1: first pass, generates stats file , 2: second pass, uses stats file)", set_pass},
    {SINGLE_INPUT, TWO_PASS_STATS_TOKEN, "Filename for 2 pass stats(\"svtav1_2pass.log\" : [Default])", set_two_pass_stats},
    {SINGLE_INPUT, PASSES_TOKEN, "Number of passes (1: one pass encode, 2: two passes encode)", set_passes},
    {SINGLE_INPUT, VBR_BIAS_PCT_TOKEN, "CBR/VBR bias (0=CBR, 100=VBR)", set_vbr_bias_pct},
    {SINGLE_INPUT, VBR_MIN_SECTION_PCT_TOKEN, "GOP min bitrate (% of target)", set_vbr_min_section_pct},
    {SINGLE_INPUT, VBR_MAX_SECTION_PCT_TOKEN, "GOP max bitrate (% of target)", set_vbr_max_section_pct},
    // Termination
    {SINGLE_INPUT, NULL, NULL, NULL}};
ConfigEntry config_entry_intra_refresh[] = {
    {SINGLE_INPUT,
     KEYINT_TOKEN,
     "Intra period interval(frames) (-2: default intra period , -1: No intra update or [0 - 2^31-2]; [-2-255] if RateControlMode >= 1)",
     set_cfg_intra_period},
    {SINGLE_INPUT,
     INTRA_REFRESH_TYPE_TOKEN,
     "Intra refresh type (1: CRA (Open GOP)2: IDR (Closed GOP))",
     set_tile_row},
    // Termination
    {SINGLE_INPUT, NULL, NULL, NULL}};
ConfigEntry config_entry_specific[] = {
    // Prediction Structure
    //{SINGLE_INPUT, ENCMODE_TOKEN, "Encoder mode/Preset used[0-8]", set_enc_mode},
    {SINGLE_INPUT, PRESET_TOKEN, "Encoder mode/Preset used[-2,-1,0,..,8]", set_enc_mode},
    {SINGLE_INPUT,
     INPUT_COMPRESSED_TEN_BIT_FORMAT,
     "Offline packing of the 2bits: requires two bits packed input (0: OFF[default], 1: ON)",
     set_compressed_ten_bit_format},
    {SINGLE_INPUT, TILE_ROW_TOKEN, "Number of tile rows to use, log2[0-6]", set_tile_row},
    {SINGLE_INPUT, TILE_COL_TOKEN, "Number of tile columns to use, log2[0-4]", set_tile_col},
    {SINGLE_INPUT, QP_TOKEN, "Constant/Constrained Quality level", set_cfg_qp},
    {SINGLE_INPUT, QP_LONG_TOKEN, "Constant/Constrained Quality level", set_cfg_qp},

    {SINGLE_INPUT,
     LOOKAHEAD_NEW_TOKEN,
     "When RC is ON , it is best to set this parameter to be equal to the intra period value",
     set_look_ahead_distance},
    // DLF
    {SINGLE_INPUT,
     LOOP_FILTER_DISABLE_NEW_TOKEN,
     "Disable loop filter(0: loop filter enabled[default] ,1: loop filter disabled)",
     set_disable_dlf_flag},
    // CDEF
     {SINGLE_INPUT,
     CDEF_LEVEL_TOKEN,
     "CDEF Level, 0: OFF, 1-5: ON with 64,16,8,4,1 step refinement, -1: DEFAULT",
     set_cdef_level},
    // RESTORATION
    {SINGLE_INPUT,
     RESTORATION_ENABLE_NEW_TOKEN,
     "Enable the loop restoration filter(0: OFF ,1: ON ,-1:DEFAULT)",
     set_enable_restoration_filter_flag},
    {SINGLE_INPUT,
     SG_FILTER_MODE_TOKEN,
     "Self-guided filter mode (0:OFF, 1: step 0, 2: step 1, 3: step 4, 4: step 16, -1: DEFAULT)",
     set_sg_filter_mode},
    {SINGLE_INPUT,
     WN_FILTER_MODE_TOKEN,
     "Wiener filter mode (0:OFF, 1: 3-Tap luma/ 3-Tap chroma, 2: 5-Tap luma/ 5-Tap chroma, 3: "
     "7-Tap luma/ 7-Tap chroma, -1: DEFAULT)",
     set_wn_filter_mode},
    {SINGLE_INPUT,
     MRP_LEVEL_TOKEN,
     "Multi reference frame levels( 0: OFF, 1: FULL, 2: Level1 .. 9: Level8,  -1: DEFAULT)",
     set_mrp_level},
    {SINGLE_INPUT, LOOK_AHEAD_DIST_TOKEN,
    "Set look ahead distance",
    set_look_ahead_distance},
    {SINGLE_INPUT,
    ENABLE_TPL_LA_TOKEN,
    "RDO based on frame temporal dependency (0: off, 1: backward source based)",
    set_enable_tpl_la},
    {SINGLE_INPUT,
     MFMV_ENABLE_NEW_TOKEN,
     "Enable motion field motion vector( 0: OFF, 1: ON, -1: DEFAULT)",
     set_enable_mfmv_flag},
    {SINGLE_INPUT,
     REDUNDANT_BLK_NEW_TOKEN,
     "Use the same md results(mode, residual , cost,etc..)as the previously processed identical "
     "block(0: OFF, 1: ON, -1: DEFAULT)",
     set_enable_redundant_blk_flag},
     {SINGLE_INPUT,
      SPATIAL_SSE_FL_NEW_TOKEN,
      "Enable spatial sse full loop(0: OFF, 1: ON, -1: DEFAULT)",
      set_spatial_sse_full_loop_level_flag},
    {SINGLE_INPUT,
     OVR_BNDRY_BLK_NEW_TOKEN,
     "Enable over boundary block mode (0: OFF, 1: ON, -1: DEFAULT)",
     set_over_bndry_blk_flag},
    {SINGLE_INPUT,
     NEW_NEAREST_COMB_INJECT_NEW_TOKEN,
     "Enable new nearest near comb injection (0: OFF, 1: ON, -1: DEFAULT)",
     set_new_nearest_comb_inject_flag},
    {SINGLE_INPUT,
     NSQ_TABLE_NEW_TOKEN,
     "Enable nsq table (0: OFF, 1: ON, -1: DEFAULT)",
     set_nsq_table_flag},
    {SINGLE_INPUT,
     FRAME_END_CDF_UPDATE_NEW_TOKEN,
     "Enable frame end cdf update mode (0: OFF, 1: ON, -1: DEFAULT)",
     set_frame_end_cdf_update_flag},

    // CHROMA
    {SINGLE_INPUT, CHROMA_MODE_TOKEN, "Select chroma mode([0-3], -1: DEFAULT)", set_chroma_mode},
    {SINGLE_INPUT,
     DISABLE_CFL_NEW_TOKEN,
     "Disable chroma from luma (CFL) flag (0: OFF (do not disable), 1: ON (disable), -1: DEFAULT)",
     set_disable_cfl_flag},

    // LOCAL WARPED MOTION
    {SINGLE_INPUT,
     LOCAL_WARPED_ENABLE_NEW_TOKEN,
     "Enable warped motion use , 0 = OFF, 1 = ON, -1 = DEFAULT",
     set_enable_local_warped_motion_flag},
    // GLOBAL MOTION
    {SINGLE_INPUT,
     GLOBAL_MOTION_ENABLE_NEW_TOKEN,
     "Enable global motion (0: OFF, 1: ON [default])",
     set_enable_global_motion_flag},
    // INTRA ANGLE DELTA
    {SINGLE_INPUT,
     INTRA_ANGLE_DELTA_TOKEN,
     "Enable intra angle delta filtering filtering (0: OFF, 1: ON, -1: DEFAULT)",
     set_intra_angle_delta_flag},
    // INTER INTRA COMPOUND
    {SINGLE_INPUT,
     INTER_INTRA_COMPOUND_NEW_TOKEN,
     "Enable interintra compound (0: OFF, 1: ON (default))",
     set_interintra_compound_flag},
    // PAETH
    {SINGLE_INPUT,
     PAETH_NEW_TOKEN,
     "Enable paeth (0: OFF, 1: ON, -1: DEFAULT)",
     set_enable_paeth_flag},
    // SMOOTH
    {SINGLE_INPUT,
     SMOOTH_NEW_TOKEN,
     "Enable smooth (0: OFF, 1: ON, -1: DEFAULT)",
     set_enable_smooth_flag},
    // OBMC
     {SINGLE_INPUT, OBMC_TOKEN, "OBMC Level(0: OFF, 1: Fully ON, 2 and 3 are faster levels, -1: DEFAULT)", set_obmc_level_flag},
    // RDOQ
    {SINGLE_INPUT,
     RDOQ_NEW_TOKEN,
     "Enable RDOQ (0: OFF, 1: ON, -1: DEFAULT)",
     set_rdoq_level_flag},

    // Filter Intra
    {SINGLE_INPUT,
     FILTER_INTRA_NEW_TOKEN,
     "Enable filter intra prediction mode (0: OFF, 1: ON [default])",
     set_filter_intra_level_flag},

    // Edge Intra Filter
    {SINGLE_INPUT,
     INTRA_EDGE_FILTER_NEW_TOKEN,
     "Enable intra edge filter (0: OFF, 1: ON, -1: DEFAULT)",
     set_enable_intra_edge_filter_flag},
    // Picture based rate estimation
    {SINGLE_INPUT,
     PIC_BASED_RATE_EST_NEW_TOKEN,
     "Enable picture based rate estimation (0: OFF, 1: ON, -1: DEFAULT)",
     set_pic_based_rate_est},

    // PREDICTIVE ME
    {SINGLE_INPUT,
     PRED_ME_TOKEN,
     "Set predictive motion estimation level(-1: default, [0-5])",
     set_predictive_me_flag},
    // BIPRED 3x3 INJECTION
    {SINGLE_INPUT,
     BIPRED_3x3_TOKEN,
     "Set bipred3x3 injection (0: OFF, 1: ON FULL, 2: Reduced set, -1: DEFAULT)",
     set_bipred3x3inject_flag},
    // COMPOUND MODE
    {SINGLE_INPUT,
     COMPOUND_LEVEL_TOKEN,
     "Enable compound mode(0: OFF, 1:ON[AVG/DIST/DIFF], 2: ON[AVG/DIST/DIFF/WEDGE], -1: default)",
     set_compound_level_flag},
    // ME Tools
    {SINGLE_INPUT,
     USE_DEFAULT_ME_HME_TOKEN,
     "Use default motion estimation/hierarchical motion estimation settings(0: OFF, 1: "
     "ON[default])",
     set_cfg_use_default_me_hme},
    {SINGLE_INPUT,
     HME_ENABLE_TOKEN,
     "Enable hierarchical motion estimation(0: OFF, 1: ON)",
     set_enable_hme_flag},
    {SINGLE_INPUT,
     HME_L0_ENABLE_TOKEN,
     "Enable hierarchical motion estimation Level 0 (0: OFF, 1: ON)",
     set_enable_hme_level_0_flag},
    {SINGLE_INPUT,
     HME_L1_ENABLE_TOKEN,
     "Enable hierarchical motion estimation Level 1 (0: OFF, 1: ON)",
     set_enable_hme_level_1_flag},
    {SINGLE_INPUT,
     HME_L2_ENABLE_TOKEN,
     "Enable hierarchical motion estimation Level 2 (0: OFF, 1: ON)",
     set_enable_hme_level_2_flag},
    {SINGLE_INPUT,
     EXT_BLOCK,
     "Enable the rectangular and asymetric block (0: OFF, 1: ON)",
     set_enable_ext_block_flag},
    // ME Parameters
    {SINGLE_INPUT,
     SEARCH_AREA_WIDTH_TOKEN,
     "Set search area in width[1-256]",
     set_cfg_search_area_width},
    {SINGLE_INPUT,
     SEARCH_AREA_HEIGHT_TOKEN,
     "Set search area in height[1-256]",
     set_cfg_search_area_height},
    // MD Parameters
    {SINGLE_INPUT,
     SCREEN_CONTENT_TOKEN,
     "Set screen content detection level([0-2], 0: DEFAULT)",
     set_screen_content_mode},
    {SINGLE_INPUT,
     INTRABC_MODE_TOKEN,
     "Set intraBC mode (0: OFF, 1: ON slow, 2: ON faster, 3: ON fastest, -1: DEFAULT)",
     set_intrabc_mode},
    {SINGLE_INPUT,
     HBD_MD_ENABLE_TOKEN,
     "Enable high bit depth mode decision(0: OFF, 1: ON partially[default],2: fully ON)",
     set_enable_hbd_mode_decision},
    {SINGLE_INPUT,
     PALETTE_TOKEN,
     "Set palette prediction mode(-1: default or [0-6])",
     set_palette_level},
    // Optional Features
    {SINGLE_INPUT,
     UNRESTRICTED_MOTION_VECTOR,
     "Allow motion vectors to reach outside of the picture boundary(O: OFF, 1: ON[default])",
     set_unrestricted_motion_vector},

    //{ SINGLE_INPUT, BITRATE_REDUCTION_TOKEN, "bit_rate_reduction", SetBitRateReduction },
    // Latency
    {SINGLE_INPUT,
     INJECTOR_TOKEN,
     "Inject pictures at defined frame rate(0: OFF[default],1: ON)",
     set_injector},
    {SINGLE_INPUT, INJECTOR_FRAMERATE_TOKEN, "Set injector frame rate", set_injector_frame_rate},
    {SINGLE_INPUT,
     SPEED_CONTROL_TOKEN,
     "Enable speed control(0: OFF[default], 1: ON)",
     speed_control_flag},
    // Annex A parameters
    {SINGLE_INPUT,
     FILM_GRAIN_TOKEN,
     "Enable film grain(0: OFF[default], 1-50: ON, film-grain denoising strength)",
     set_cfg_film_grain},
    // --- start: ALTREF_FILTERING_SUPPORT
     {SINGLE_INPUT,
      TF_LEVEL,
      "Set altref level(-1: Default; 0: OFF; 1: ON; 2 and 3: Faster levels)",
      set_tf_level},
    {SINGLE_INPUT,
     ALTREF_STRENGTH,
     "AltRef filter strength([0-6], default: 5)",
     set_altref_strength},
    {SINGLE_INPUT, ALTREF_NFRAMES, "AltRef max frames([0-10], default: 7)", set_altref_n_frames},
    {SINGLE_INPUT,
     ENABLE_OVERLAYS,
     "Enable the insertion of an extra picture called overlayer picture which will be used as an "
     "extra reference frame for the base-layer picture(0: OFF[default], 1: ON)",
     set_enable_overlays},
    // --- end: ALTREF_FILTERING_SUPPORT
    {SINGLE_INPUT, STAT_REPORT_NEW_TOKEN, "Stat Report", set_stat_report},
    {SINGLE_INPUT,
     INTRA_ANGLE_DELTA_NEW_TOKEN,
     "Enable intra angle delta filtering filtering (0: OFF, 1: ON, -1: DEFAULT)",
     set_intra_angle_delta_flag},

    // double dash
    //{SINGLE_INPUT,
    //NX4_4XN_MV_INJECT_NEW_TOKEN,
    // "nx4ParentMvInjection",
    // set_nx4_4xn_parent_mv_inject_flag},

    // Termination
    {SINGLE_INPUT, NULL, NULL, NULL}};

ConfigEntry config_entry[] = {
    // File I/O
    {SINGLE_INPUT, INPUT_FILE_TOKEN, "InputFile", set_cfg_input_file},
    {SINGLE_INPUT, OUTPUT_BITSTREAM_TOKEN, "StreamFile", set_cfg_stream_file},
    {SINGLE_INPUT, ERROR_FILE_TOKEN, "ErrorFile", set_cfg_error_file},
    {SINGLE_INPUT, OUTPUT_RECON_TOKEN, "ReconFile", set_cfg_recon_file},
    {SINGLE_INPUT, QP_FILE_TOKEN, "QpFile", set_cfg_qp_file},
    {SINGLE_INPUT, STAT_FILE_TOKEN, "StatFile", set_cfg_stat_file},

    // two pass
    {SINGLE_INPUT, PASS_TOKEN, "Pass", set_pass},
    {SINGLE_INPUT, TWO_PASS_STATS_TOKEN, "Two pass stat", set_two_pass_stats},

    {SINGLE_INPUT, INPUT_PREDSTRUCT_FILE_TOKEN, "PredStructFile", set_pred_struct_file},
    // Picture Dimensions
    {SINGLE_INPUT, WIDTH_TOKEN, "SourceWidth", set_cfg_source_width},
    {SINGLE_INPUT, HEIGHT_TOKEN, "SourceHeight", set_cfg_source_height},
    // Prediction Structure
    {SINGLE_INPUT, NUMBER_OF_PICTURES_TOKEN, "FrameToBeEncoded", set_cfg_frames_to_be_encoded},
    {SINGLE_INPUT, BUFFERED_INPUT_TOKEN, "BufferedInput", set_buffered_input},
    {SINGLE_INPUT, PROGRESS_TOKEN, "Progress", set_progress},
    {SINGLE_INPUT, NO_PROGRESS_TOKEN, "NoProgress", set_no_progress},
    {SINGLE_INPUT, ENCMODE_TOKEN, "EncoderMode", set_enc_mode},
    {SINGLE_INPUT, INTRA_PERIOD_TOKEN, "IntraPeriod", set_cfg_intra_period},
    {SINGLE_INPUT, INTRA_REFRESH_TYPE_TOKEN, "IntraRefreshType", set_cfg_intra_refresh_type},
    {SINGLE_INPUT, FRAME_RATE_TOKEN, "FrameRate", set_frame_rate},
    {SINGLE_INPUT, FRAME_RATE_NUMERATOR_TOKEN, "FrameRateNumerator", set_frame_rate_numerator},
    {SINGLE_INPUT,
     FRAME_RATE_DENOMINATOR_TOKEN,
     "FrameRateDenominator",
     set_frame_rate_denominator},
    {SINGLE_INPUT, ENCODER_BIT_DEPTH, "EncoderBitDepth", set_encoder_bit_depth},
    {SINGLE_INPUT, ENCODER_16BIT_PIPELINE, "Encoder16BitPipeline", set_encoder_16bit_pipeline},
    {SINGLE_INPUT, ENCODER_COLOR_FORMAT, "EncoderColorFormat", set_encoder_color_format},
    {SINGLE_INPUT,
     INPUT_COMPRESSED_TEN_BIT_FORMAT,
     "CompressedTenBitFormat",
     set_compressed_ten_bit_format},
    {SINGLE_INPUT, HIERARCHICAL_LEVELS_TOKEN, "HierarchicalLevels", set_hierarchical_levels},
    {SINGLE_INPUT, PRED_STRUCT_TOKEN, "PredStructure", set_cfg_pred_structure},
    {SINGLE_INPUT, TILE_ROW_TOKEN, "TileRow", set_tile_row},
    {SINGLE_INPUT, TILE_COL_TOKEN, "TileCol", set_tile_col},
    // Rate Control
    {SINGLE_INPUT,
     SCENE_CHANGE_DETECTION_TOKEN,
     "SceneChangeDetection",
     set_scene_change_detection},
    {SINGLE_INPUT, QP_TOKEN, "QP", set_cfg_qp},
    {SINGLE_INPUT, USE_QP_FILE_TOKEN, "UseQpFile", set_cfg_use_qp_file},
    {SINGLE_INPUT, STAT_REPORT_TOKEN, "StatReport", set_stat_report},
    {SINGLE_INPUT, RATE_CONTROL_ENABLE_TOKEN, "RateControlMode", set_rate_control_mode},
    {SINGLE_INPUT, LOOK_AHEAD_DIST_TOKEN, "LookAheadDistance", set_look_ahead_distance},
    {SINGLE_INPUT, ENABLE_TPL_LA_TOKEN, "EnableTplLA", set_enable_tpl_la},
    {SINGLE_INPUT, TARGET_BIT_RATE_TOKEN, "TargetBitRate", set_target_bit_rate},
    {SINGLE_INPUT, MAX_QP_TOKEN, "MaxQpAllowed", set_max_qp_allowed},
    {SINGLE_INPUT, MIN_QP_TOKEN, "MinQpAllowed", set_min_qp_allowed},
    {SINGLE_INPUT, VBV_BUFSIZE_TOKEN, "VBVBufSize", set_vbv_buf_size},
    {SINGLE_INPUT, ADAPTIVE_QP_ENABLE_TOKEN, "AdaptiveQuantization", set_adaptive_quantization},
    {SINGLE_INPUT, VBR_BIAS_PCT_TOKEN, "CBR/VBR bias (0=CBR, 100=VBR)", set_vbr_bias_pct},
    {SINGLE_INPUT, VBR_MIN_SECTION_PCT_TOKEN, "GOP min bitrate (% of target)", set_vbr_min_section_pct},
    {SINGLE_INPUT, VBR_MAX_SECTION_PCT_TOKEN, "GOP max bitrate (% of target)", set_vbr_max_section_pct},
    {SINGLE_INPUT, UNDER_SHOOT_PCT_TOKEN, "Datarate undershoot (min) target (%)", set_under_shoot_pct},
    {SINGLE_INPUT, OVER_SHOOT_PCT_TOKEN, "Datarate overshoot (max) target (%)", set_over_shoot_pct},

    // DLF
    {SINGLE_INPUT, LOOP_FILTER_DISABLE_TOKEN, "LoopFilterDisable", set_disable_dlf_flag},

    // CDEF
    {SINGLE_INPUT, CDEF_LEVEL_TOKEN, "CDEFLevel", set_cdef_level},

    // RESTORATION
    {SINGLE_INPUT,
     RESTORATION_ENABLE_TOKEN,
     "RestorationFilter",
     set_enable_restoration_filter_flag},
    {SINGLE_INPUT, SG_FILTER_MODE_TOKEN, "SelfGuidedFilterMode", set_sg_filter_mode},
    {SINGLE_INPUT, WN_FILTER_MODE_TOKEN, "WienerFilterMode", set_wn_filter_mode},
    {SINGLE_INPUT, MRP_LEVEL_TOKEN, "MrpLevel", set_mrp_level},
    {SINGLE_INPUT, MFMV_ENABLE_TOKEN, "Mfmv", set_enable_mfmv_flag},
    {SINGLE_INPUT, REDUNDANT_BLK_TOKEN, "RedundantBlock", set_enable_redundant_blk_flag},
    {SINGLE_INPUT, SPATIAL_SSE_FL_TOKEN, "SpatialSSEfl", set_spatial_sse_full_loop_level_flag},
    {SINGLE_INPUT, OVR_BNDRY_BLK_TOKEN, "OverBoundryBlock", set_over_bndry_blk_flag},
    {SINGLE_INPUT,
     NEW_NEAREST_COMB_INJECT_TOKEN,
     "NewNearestCombInjection",
     set_new_nearest_comb_inject_flag},
    {SINGLE_INPUT, NSQ_TABLE_TOKEN, "NsqTable", set_nsq_table_flag},
    {SINGLE_INPUT, FRAME_END_CDF_UPDATE_TOKEN, "FrameEndCdfUpdate", set_frame_end_cdf_update_flag},

    // CHROMA
    {SINGLE_INPUT, CHROMA_MODE_TOKEN, "ChromaMode", set_chroma_mode},
    {SINGLE_INPUT, DISABLE_CFL_TOKEN, "DisableCFL", set_disable_cfl_flag},

    // LOCAL WARPED MOTION
    {SINGLE_INPUT,
     LOCAL_WARPED_ENABLE_TOKEN,
     "LocalWarpedMotion",
     set_enable_local_warped_motion_flag},
    // GLOBAL MOTION
    {SINGLE_INPUT, GLOBAL_MOTION_ENABLE_TOKEN, "GlobalMotion", set_enable_global_motion_flag},
    // INTRA ANGLE DELTA
    {SINGLE_INPUT, INTRA_ANGLE_DELTA_TOKEN, "IntraAngleDelta", set_intra_angle_delta_flag},

    // INTER INTRA COMPOUND
    {SINGLE_INPUT, INTER_INTRA_COMPOUND_TOKEN, "InterIntraCompound", set_interintra_compound_flag},
    // PAETH
    {SINGLE_INPUT, PAETH_TOKEN, "Paeth", set_enable_paeth_flag},
    // SMOOTH
    {SINGLE_INPUT, SMOOTH_TOKEN, "Smooth", set_enable_smooth_flag},
    // OBMC
    {SINGLE_INPUT, OBMC_TOKEN, "Obmc", set_obmc_level_flag},
    // RDOQ
    {SINGLE_INPUT, RDOQ_TOKEN, "RDOQ", set_rdoq_level_flag},
    // Filter Intra
    {SINGLE_INPUT, FILTER_INTRA_TOKEN, "FilterIntra", set_filter_intra_level_flag},
    // Edge Intra Filter
    {SINGLE_INPUT, INTRA_EDGE_FILTER_TOKEN, "IntraEdgeFilter", set_enable_intra_edge_filter_flag},

    // Picture based rate estimation
    {SINGLE_INPUT, PIC_BASED_RATE_EST_TOKEN, "PicBasedRateEst", set_pic_based_rate_est},
    {SINGLE_INPUT, PIC_BASED_RATE_EST_NEW_TOKEN, "PicBasedRateEst", set_pic_based_rate_est},

    // PREDICTIVE ME
    {SINGLE_INPUT, PRED_ME_TOKEN, "PredMe", set_predictive_me_flag},
    // BIPRED 3x3 INJECTION
    {SINGLE_INPUT, BIPRED_3x3_TOKEN, "Bipred3x3", set_bipred3x3inject_flag},
    // COMPOUND MODE
    {SINGLE_INPUT, COMPOUND_LEVEL_TOKEN, "CompoundLevel", set_compound_level_flag},

    // ME Tools
    {SINGLE_INPUT, USE_DEFAULT_ME_HME_TOKEN, "UseDefaultMeHme", set_cfg_use_default_me_hme},
    {SINGLE_INPUT, HME_ENABLE_TOKEN, "HME", set_enable_hme_flag},
    {SINGLE_INPUT, HME_L0_ENABLE_TOKEN, "HMELevel0", set_enable_hme_level_0_flag},
    {SINGLE_INPUT, HME_L1_ENABLE_TOKEN, "HMELevel1", set_enable_hme_level_1_flag},
    {SINGLE_INPUT, HME_L2_ENABLE_TOKEN, "HMELevel2", set_enable_hme_level_2_flag},
    {SINGLE_INPUT, EXT_BLOCK, "ExtBlockFlag", set_enable_ext_block_flag},
    // ME Parameters
    {SINGLE_INPUT, SEARCH_AREA_WIDTH_TOKEN, "SearchAreaWidth", set_cfg_search_area_width},
    {SINGLE_INPUT, SEARCH_AREA_HEIGHT_TOKEN, "SearchAreaHeight", set_cfg_search_area_height},

    // MD Parameters
    {SINGLE_INPUT, SCREEN_CONTENT_TOKEN, "ScreenContentMode", set_screen_content_mode},
    {SINGLE_INPUT, INTRABC_MODE_TOKEN, "IntraBCMode", set_intrabc_mode},
    {SINGLE_INPUT, HBD_MD_ENABLE_TOKEN, "HighBitDepthModeDecision", set_enable_hbd_mode_decision},
    {SINGLE_INPUT, PALETTE_TOKEN, "PaletteLevel", set_palette_level},
    // Thread Management
    {SINGLE_INPUT, THREAD_MGMNT, "LogicalProcessors", set_logical_processors},
    {SINGLE_INPUT, UNPIN_TOKEN, "UnpinExecution", set_unpin_execution},
    {SINGLE_INPUT, TARGET_SOCKET, "TargetSocket", set_target_socket},
    // Optional Features
    {SINGLE_INPUT,
     UNRESTRICTED_MOTION_VECTOR,
     "UnrestrictedMotionVector",
     set_unrestricted_motion_vector},

    //    { SINGLE_INPUT, BITRATE_REDUCTION_TOKEN, "bit_rate_reduction", SetBitRateReduction },
    {SINGLE_INPUT, HDR_INPUT_TOKEN, "HighDynamicRangeInput", set_high_dynamic_range_input},
    // Latency
    {SINGLE_INPUT, INJECTOR_TOKEN, "Injector", set_injector},
    {SINGLE_INPUT, INJECTOR_FRAMERATE_TOKEN, "InjectorFrameRate", set_injector_frame_rate},
    {SINGLE_INPUT, SPEED_CONTROL_TOKEN, "SpeedControlFlag", speed_control_flag},
    // Annex A parameters
    {SINGLE_INPUT, PROFILE_TOKEN, "Profile", set_profile},
    {SINGLE_INPUT, TIER_TOKEN, "Tier", set_tier},
    {SINGLE_INPUT, LEVEL_TOKEN, "Level", set_level},
    {SINGLE_INPUT, FILM_GRAIN_TOKEN, "FilmGrain", set_cfg_film_grain},
    // Asm Type
    {SINGLE_INPUT, ASM_TYPE_TOKEN, "Asm", set_asm_type},
    // --- start: ALTREF_FILTERING_SUPPORT
    {SINGLE_INPUT, TF_LEVEL, "TfLevel", set_tf_level},
    {SINGLE_INPUT, ALTREF_STRENGTH, "AltRefStrength", set_altref_strength},
    {SINGLE_INPUT, ALTREF_NFRAMES, "AltRefNframes", set_altref_n_frames},
    {SINGLE_INPUT, ENABLE_OVERLAYS, "EnableOverlays", set_enable_overlays},
    // --- end: ALTREF_FILTERING_SUPPORT
    // Super-resolution support
    {SINGLE_INPUT, SUPERRES_MODE_INPUT, "SuperresMode", set_superres_mode},
    {SINGLE_INPUT, SUPERRES_DENOM, "SuperresDenom", set_superres_denom},
    {SINGLE_INPUT, SUPERRES_KF_DENOM, "SuperresKfDenom", set_superres_kf_denom},
    {SINGLE_INPUT, SUPERRES_QTHRES, "SuperresQthres", set_superres_qthres},

    // double dash
    {SINGLE_INPUT, PRESET_TOKEN, "Encoder mode/Preset used[-2,-1,0,..,8]", set_enc_mode},
    {SINGLE_INPUT, QP_FILE_NEW_TOKEN, "Path to Qp file", set_cfg_qp_file},
    {SINGLE_INPUT, INPUT_DEPTH_TOKEN, "Bit depth for codec(8 or 10)", set_encoder_bit_depth},
    {SINGLE_INPUT,
     KEYINT_TOKEN,
     "Intra period interval(frames) (-2: default intra period, -1: No intra update or [0 - 2^31-2]; [-2 - 255] if RateControlMode>=1)",
     set_cfg_intra_period},
    {SINGLE_INPUT,
     LOOKAHEAD_NEW_TOKEN,
     "When RC is ON , it is best to set this parameter to be equal to the intra period value",
     set_look_ahead_distance},

    {SINGLE_INPUT, STAT_REPORT_NEW_TOKEN, "Stat Report", set_stat_report},
    {SINGLE_INPUT,
     RESTORATION_ENABLE_NEW_TOKEN,
     "Restoration Filter",
     set_enable_restoration_filter_flag},
    {SINGLE_INPUT,
     INTER_INTRA_COMPOUND_NEW_TOKEN,
     "Inter Intra Compound",
     set_interintra_compound_flag},
    //{SINGLE_INPUT, FRAC_SEARCH_64_NEW_TOKEN, "FracSearch64", ??}, // todo
    {SINGLE_INPUT, MFMV_ENABLE_NEW_TOKEN, "Mfmv token with double dash", set_enable_mfmv_flag},
    {SINGLE_INPUT, REDUNDANT_BLK_NEW_TOKEN, "Redundant Block", set_enable_redundant_blk_flag},
    {SINGLE_INPUT, SPATIAL_SSE_FL_NEW_TOKEN, "Spatial SSE fl", set_spatial_sse_full_loop_level_flag},
    {SINGLE_INPUT, OVR_BNDRY_BLK_NEW_TOKEN, "Over Boundry Block", set_over_bndry_blk_flag},
    {SINGLE_INPUT,
     NEW_NEAREST_COMB_INJECT_NEW_TOKEN,
     "New Nearest Comb Injection",
     set_new_nearest_comb_inject_flag},
    {SINGLE_INPUT, NSQ_TABLE_NEW_TOKEN, "Nsq Table", set_nsq_table_flag},
    {SINGLE_INPUT,
     FRAME_END_CDF_UPDATE_NEW_TOKEN,
     "Frame End Cdf Update",
     set_frame_end_cdf_update_flag},
    {SINGLE_INPUT,
     LOCAL_WARPED_ENABLE_NEW_TOKEN,
     "Local Warped Motion [0 = OFF, 1 = ON, -1 = DEFAULT]",
     set_enable_local_warped_motion_flag},
    {SINGLE_INPUT, GLOBAL_MOTION_ENABLE_NEW_TOKEN, "Global Motion", set_enable_global_motion_flag},
    {SINGLE_INPUT, RDOQ_NEW_TOKEN, "RDOQ double dash token", set_rdoq_level_flag},
    {SINGLE_INPUT, FILTER_INTRA_NEW_TOKEN, "Filter Intra", set_filter_intra_level_flag},
    {SINGLE_INPUT, HDR_INPUT_NEW_TOKEN, "High Dynamic Range Input", set_high_dynamic_range_input},
    {SINGLE_INPUT,
     ADAPTIVE_QP_ENABLE_NEW_TOKEN,
     "Adaptive Quantization",
     set_adaptive_quantization},

    {SINGLE_INPUT, INPUT_FILE_LONG_TOKEN, "Input File", set_cfg_input_file},
    {SINGLE_INPUT, OUTPUT_BITSTREAM_LONG_TOKEN, "Stream File", set_cfg_stream_file},
    {SINGLE_INPUT, OUTPUT_RECON_LONG_TOKEN, "Recon File", set_cfg_recon_file},

    {SINGLE_INPUT, WIDTH_LONG_TOKEN, "Source Width", set_cfg_source_width},
    {SINGLE_INPUT, HEIGHT_LONG_TOKEN, "Source Height", set_cfg_source_height},
    {SINGLE_INPUT,
     NUMBER_OF_PICTURES_LONG_TOKEN,
     "Frame To Be Encoded",
     set_cfg_frames_to_be_encoded},
    {SINGLE_INPUT, QP_LONG_TOKEN, "QP double dash token", set_cfg_qp},
    {SINGLE_INPUT, LOOP_FILTER_DISABLE_NEW_TOKEN, "Loop Filter Disable", set_disable_dlf_flag},

    {SINGLE_INPUT, DISABLE_CFL_NEW_TOKEN, "Disable CFL", set_disable_cfl_flag},
    {SINGLE_INPUT,
     INTRA_EDGE_FILTER_NEW_TOKEN,
     "Intra Edge Filter",
     set_enable_intra_edge_filter_flag},
    {SINGLE_INPUT, INTRA_ANGLE_DELTA_NEW_TOKEN, "Intra Angle Delta", set_intra_angle_delta_flag},
    {SINGLE_INPUT, PAETH_NEW_TOKEN, "Paeth New Token", set_enable_paeth_flag},
    {SINGLE_INPUT, SMOOTH_NEW_TOKEN, "Smooth New Token", set_enable_smooth_flag},

    // Termination
    {SINGLE_INPUT, NULL, NULL, NULL}};

/**********************************
 * Constructor
 **********************************/
EbConfig * eb_config_ctor(EncodePass pass) {
    EbConfig *config_ptr = (EbConfig *)calloc(1, sizeof(EbConfig));
    if (!config_ptr)
        return NULL;

    if (pass == ENCODE_FIRST_PASS)
        config_ptr->pass = 1;
    else if (pass == ENCODE_LAST_PASS)
        config_ptr->pass = 2;

    config_ptr->error_log_file         = stderr;
    config_ptr->frame_rate             = 30 << 16;
    config_ptr->encoder_bit_depth      = 8;
    config_ptr->is_16bit_pipeline = 0;
    config_ptr->encoder_color_format   = 1; //EB_YUV420
    config_ptr->buffered_input         = -1;
    config_ptr->qp                  = 50;
    config_ptr->use_qp_file         = EB_FALSE;
    config_ptr->look_ahead_distance = (uint32_t)~0;
    config_ptr->enable_tpl_la       = 1;
    config_ptr->target_bit_rate     = 7000000;
    config_ptr->max_qp_allowed      = 63;
    config_ptr->min_qp_allowed      = 10;

    config_ptr->enable_adaptive_quantization              = 2;
    config_ptr->enc_mode                                  = MAX_ENC_PRESET;
    config_ptr->vbr_bias_pct        = 50;
    config_ptr->vbr_min_section_pct = 0;
    config_ptr->vbr_max_section_pct = 2000;
    config_ptr->under_shoot_pct     = 25;
    config_ptr->over_shoot_pct      = 25;
    config_ptr->intra_period                              = -2;
    config_ptr->intra_refresh_type                        = 1;
    config_ptr->hierarchical_levels                       = 4;
    config_ptr->pred_structure                            = 2;
    config_ptr->disable_dlf_flag                          = EB_FALSE;
    config_ptr->enable_global_motion                      = EB_TRUE;
    config_ptr->progress                                  = 1;
    config_ptr->enable_warped_motion                      = DEFAULT;
    config_ptr->cdef_level                                = DEFAULT;
    config_ptr->enable_restoration_filtering              = DEFAULT;
    config_ptr->sg_filter_mode                            = DEFAULT;
    config_ptr->wn_filter_mode                            = DEFAULT;
    config_ptr->intra_angle_delta                         = DEFAULT;
    config_ptr->inter_intra_compound                      = DEFAULT;
    config_ptr->enable_paeth                              = DEFAULT;
    config_ptr->enable_smooth                             = DEFAULT;
    config_ptr->mrp_level                                 = DEFAULT;
    config_ptr->enable_mfmv                               = DEFAULT;
    config_ptr->enable_redundant_blk                      = DEFAULT;
    config_ptr->spatial_sse_full_loop_level               = DEFAULT;
    config_ptr->over_bndry_blk                            = DEFAULT;
    config_ptr->new_nearest_comb_inject                   = DEFAULT;
    config_ptr->nsq_table                                 = DEFAULT;
    config_ptr->frame_end_cdf_update                      = DEFAULT;
    config_ptr->set_chroma_mode                           = DEFAULT;
    config_ptr->disable_cfl_flag                          = DEFAULT;
    config_ptr->obmc_level                                = DEFAULT;
    config_ptr->rdoq_level                                = DEFAULT;
    config_ptr->pred_me                                   = DEFAULT;
    config_ptr->bipred_3x3_inject                         = DEFAULT;
    config_ptr->compound_level                            = DEFAULT;
    config_ptr->filter_intra_level                        = DEFAULT;
    config_ptr->enable_intra_edge_filter                  = DEFAULT;
    config_ptr->pic_based_rate_est                        = DEFAULT;
    config_ptr->ext_block_flag                            = EB_FALSE;
    config_ptr->use_default_me_hme                        = EB_TRUE;
    config_ptr->enable_hme_flag                           = EB_TRUE;
    config_ptr->enable_hme_level0_flag                    = EB_TRUE;
    config_ptr->enable_hme_level1_flag                    = EB_FALSE;
    config_ptr->enable_hme_level2_flag                    = EB_FALSE;
    config_ptr->search_area_width                         = 16;
    config_ptr->search_area_height                        = 7;
    config_ptr->number_hme_search_region_in_width         = 2;
    config_ptr->number_hme_search_region_in_height        = 2;
    config_ptr->hme_level0_total_search_area_width        = 64;
    config_ptr->hme_level0_total_search_area_height       = 25;
    config_ptr->hme_level0_search_area_in_width_array[0]  = 32;
    config_ptr->hme_level0_search_area_in_width_array[1]  = 32;
    config_ptr->hme_level0_search_area_in_height_array[0] = 12;
    config_ptr->hme_level0_search_area_in_height_array[1] = 13;
    config_ptr->hme_level1_search_area_in_width_array[0]  = 1;
    config_ptr->hme_level1_search_area_in_width_array[1]  = 1;
    config_ptr->hme_level1_search_area_in_height_array[0] = 1;
    config_ptr->hme_level1_search_area_in_height_array[1] = 1;
    config_ptr->hme_level2_search_area_in_width_array[0]  = 1;
    config_ptr->hme_level2_search_area_in_width_array[1]  = 1;
    config_ptr->hme_level2_search_area_in_height_array[0] = 1;
    config_ptr->hme_level2_search_area_in_height_array[1] = 1;
    config_ptr->screen_content_mode                       = 2;
    config_ptr->enable_hbd_mode_decision                  = DEFAULT;
    config_ptr->intrabc_mode                              = DEFAULT;
    config_ptr->palette_level                             = DEFAULT;
    config_ptr->injector_frame_rate                       = 60 << 16;
    config_ptr->speed_control_flag                        = 0;

    // ASM Type
    config_ptr->cpu_flags_limit = CPU_FLAGS_ALL;

    config_ptr->unpin     = 1;
    config_ptr->target_socket = -1;

    config_ptr->unrestricted_motion_vector = EB_TRUE;

    // --- start: ALTREF_FILTERING_SUPPORT
    config_ptr->tf_level = DEFAULT;
    config_ptr->altref_strength = 5;
    config_ptr->altref_nframes = 13;
    // --- end: ALTREF_FILTERING_SUPPORT

    // start - super-resolution support
    config_ptr->superres_mode     = SUPERRES_NONE; // disabled
    config_ptr->superres_denom    = 8; // no scaling
    config_ptr->superres_kf_denom = 8; // no scaling
    config_ptr->superres_qthres   = 43; // random threshold for now
    // end - super-resolution support

    config_ptr->pass = DEFAULT;

    return config_ptr;
}

/**********************************
 * Destructor
 **********************************/
void eb_config_dtor(EbConfig *config_ptr) {
    if (!config_ptr)
        return;
    // Close any files that are open
    if (config_ptr->config_file) {
        fclose(config_ptr->config_file);
        config_ptr->config_file = (FILE *)NULL;
    }

    if (config_ptr->input_file) {
        if (!config_ptr->input_file_is_fifo) fclose(config_ptr->input_file);
        config_ptr->input_file = (FILE *)NULL;
    }

    if (config_ptr->bitstream_file) {
        fclose(config_ptr->bitstream_file);
        config_ptr->bitstream_file = (FILE *)NULL;
    }

    if (config_ptr->recon_file) {
        fclose(config_ptr->recon_file);
        config_ptr->recon_file = (FILE *)NULL;
    }

    if (config_ptr->input_pred_struct_file) {
        fclose(config_ptr->input_pred_struct_file);
        config_ptr->input_pred_struct_file = (FILE *)NULL;
    }

    if (config_ptr->input_pred_struct_filename) {
        free(config_ptr->input_pred_struct_filename);
        config_ptr->input_pred_struct_filename = NULL;
    }

    if (config_ptr->error_log_file && config_ptr->error_log_file != stderr) {
        fclose(config_ptr->error_log_file);
        config_ptr->error_log_file = (FILE *)NULL;
    }

    if (config_ptr->qp_file) {
        fclose(config_ptr->qp_file);
        config_ptr->qp_file = (FILE *)NULL;
    }

    if (config_ptr->stat_file) {
        fclose(config_ptr->stat_file);
        config_ptr->stat_file = (FILE *)NULL;
    }
    free((void*)config_ptr->stats);
    free(config_ptr);
    return;
}

EbErrorType enc_channel_ctor(EncChannel* c, EncodePass pass) {
    c->config = eb_config_ctor(pass);
    if (!c->config)
        return EB_ErrorInsufficientResources;
    c->app_callback = (EbAppContext *)malloc(sizeof(EbAppContext));
    if (!c->app_callback)
        return EB_ErrorInsufficientResources;
    memset(c->app_callback, 0, sizeof(EbAppContext));
    c->exit_cond        = APP_ExitConditionError;
    c->exit_cond_output = APP_ExitConditionError;
    c->exit_cond_recon  = APP_ExitConditionError;
    c->exit_cond_input  = APP_ExitConditionError;
    c->active = EB_FALSE;
    return EB_ErrorNone;
}

void enc_channel_dctor(EncChannel* c)
{
    eb_config_dtor(c->config);
    free(c->app_callback);
}

/**********************************
 * File Size
 **********************************/
static int32_t find_file_size(FILE *const pFile) {
    int32_t file_size;

    fseek(pFile, 0, SEEK_END);
    file_size = ftell(pFile);
    rewind(pFile);

    return file_size;
}

/**********************************
 * Line Split
 **********************************/
static void line_split(uint32_t *argc, char *argv[CONFIG_FILE_MAX_ARG_COUNT],
                       uint32_t arg_len[CONFIG_FILE_MAX_ARG_COUNT], char *linePtr) {
    uint32_t i = 0;
    *argc      = 0;

    while ((*linePtr != CONFIG_FILE_NEWLINE_CHAR) && (*linePtr != CONFIG_FILE_RETURN_CHAR) &&
           (*linePtr != CONFIG_FILE_COMMENT_CHAR) && (*argc < CONFIG_FILE_MAX_ARG_COUNT)) {
        // Increment past whitespace
        while ((*linePtr == CONFIG_FILE_SPACE_CHAR || *linePtr == CONFIG_FILE_TAB_CHAR) &&
               (*linePtr != CONFIG_FILE_NEWLINE_CHAR))
            ++linePtr;
        // Set arg
        if ((*linePtr != CONFIG_FILE_NEWLINE_CHAR) && (*linePtr != CONFIG_FILE_RETURN_CHAR) &&
            (*linePtr != CONFIG_FILE_COMMENT_CHAR) && (*argc < CONFIG_FILE_MAX_ARG_COUNT)) {
            argv[*argc] = linePtr;

            // Increment to next whitespace
            while (*linePtr != CONFIG_FILE_SPACE_CHAR && *linePtr != CONFIG_FILE_TAB_CHAR &&
                   *linePtr != CONFIG_FILE_NEWLINE_CHAR && *linePtr != CONFIG_FILE_RETURN_CHAR) {
                ++linePtr;
                ++i;
            }

            // Set arg length
            arg_len[(*argc)++] = i;

            i = 0;
        }
    }

    return;
}

/**********************************
* Set Config value
**********************************/
static void set_config_value(EbConfig *config, const char *name, const char *value) {
    int32_t i = 0;

    while (config_entry[i].name != NULL) {
        if (strcmp(config_entry[i].name, name) == 0)
            (*config_entry[i].scf)((const char *)value, config);
        ++i;
    }

    return;
}

/**********************************
* Parse Config File
**********************************/
static void parse_config_file(EbConfig *config, char *buffer, int32_t size) {
    uint32_t argc;
    char *   argv[CONFIG_FILE_MAX_ARG_COUNT];
    uint32_t arg_len[CONFIG_FILE_MAX_ARG_COUNT];

    char var_name[CONFIG_FILE_MAX_VAR_LEN];
    char var_value[CONFIG_FILE_MAX_ARG_COUNT][CONFIG_FILE_MAX_VAR_LEN];

    uint32_t value_index;

    uint32_t comment_section_flag = 0;
    uint32_t new_line_flag        = 0;

    // Keep looping until we process the entire file
    while (size--) {
        comment_section_flag =
            ((*buffer == CONFIG_FILE_COMMENT_CHAR) || (comment_section_flag != 0))
                ? 1
                : comment_section_flag;

        // At the beginning of each line
        if ((new_line_flag == 1) && (comment_section_flag == 0)) {
            // Do an argc/argv split for the line
            line_split(&argc, argv, arg_len, buffer);

            if ((argc > 2) && (*argv[1] == CONFIG_FILE_VALUE_SPLIT)) {
                // ***NOTE - We're assuming that the variable name is the first arg and
                // the variable value is the third arg.

                // Cap the length of the variable name
                arg_len[0] = (arg_len[0] > CONFIG_FILE_MAX_VAR_LEN - 1)
                                 ? CONFIG_FILE_MAX_VAR_LEN - 1
                                 : arg_len[0];
                // Copy the variable name
                strncpy_s(var_name, CONFIG_FILE_MAX_VAR_LEN, argv[0], arg_len[0]);
                // Null terminate the variable name
                var_name[arg_len[0]] = CONFIG_FILE_NULL_CHAR;

                for (value_index = 0;
                     (value_index < CONFIG_FILE_MAX_ARG_COUNT - 2) && (value_index < (argc - 2));
                     ++value_index) {
                    // Cap the length of the variable
                    arg_len[value_index + 2] =
                        (arg_len[value_index + 2] > CONFIG_FILE_MAX_VAR_LEN - 1)
                            ? CONFIG_FILE_MAX_VAR_LEN - 1
                            : arg_len[value_index + 2];
                    // Copy the variable name
                    strncpy_s(var_value[value_index],
                               CONFIG_FILE_MAX_VAR_LEN,
                               argv[value_index + 2],
                               arg_len[value_index + 2]);
                    // Null terminate the variable name
                    var_value[value_index][arg_len[value_index + 2]] = CONFIG_FILE_NULL_CHAR;

                    set_config_value(config, var_name, var_value[value_index]);
                }
            }
        }

        comment_section_flag = (*buffer == CONFIG_FILE_NEWLINE_CHAR) ? 0 : comment_section_flag;
        new_line_flag        = (*buffer == CONFIG_FILE_NEWLINE_CHAR) ? 1 : 0;
        ++buffer;
    }

    return;
}

/******************************************
* Find Token
******************************************/
static int32_t find_token(int32_t argc, char *const argv[], char const *token, char *configStr) {
    int32_t return_error = -1;

    while ((argc > 0) && (return_error != 0)) {
        return_error = strcmp(argv[--argc], token);
        if (return_error == 0 && configStr)
            strcpy_s(configStr, COMMAND_LINE_MAX_SIZE, argv[argc + 1]);
    }

    return return_error;
}

/**********************************
* Read Config File
**********************************/
static int32_t read_config_file(EbConfig *config, char *config_path, uint32_t instance_idx) {
    int32_t return_error = 0;

    // Open the config file
    FOPEN(config->config_file, config_path, "rb");

    if (config->config_file != (FILE *)NULL) {
        int32_t config_file_size   = find_file_size(config->config_file);
        char *  config_file_buffer = (char *)malloc(config_file_size);

        if (config_file_buffer != (char *)NULL) {
            int32_t result_size =
                (int32_t)fread(config_file_buffer, 1, config_file_size, config->config_file);

            if (result_size == config_file_size) {
                parse_config_file(config, config_file_buffer, config_file_size);
            } else {
                fprintf(stderr, "Error channel %u: File Read Failed\n", instance_idx + 1);
                return_error = -1;
            }
        } else {
            fprintf(stderr, "Error channel %u: Memory Allocation Failed\n", instance_idx + 1);
            return_error = -1;
        }

        free(config_file_buffer);
        fclose(config->config_file);
        config->config_file = (FILE *)NULL;
    } else {
        fprintf(stderr,
                "Error channel %u: Couldn't open Config File: %s\n",
                instance_idx + 1,
                config_path);
        return_error = -1;
    }

    return return_error;
}

/* get config->rc_twopass_stats_in from config->input_stat_file */
EbBool load_twopass_stats_in(EbConfig *config)
{
#ifdef _WIN32
    int fd = _fileno(config->input_stat_file);
    struct _stat file_stat;
    int ret = _fstat(fd, &file_stat);
#else
    int fd = fileno(config->input_stat_file);
    struct stat file_stat;
    int ret = fstat(fd, &file_stat);
#endif
    if (ret) {
        return EB_FALSE;
    }
    config->rc_twopass_stats_in.buf = malloc(file_stat.st_size);
    if (config->rc_twopass_stats_in.buf) {
        config->rc_twopass_stats_in.sz = (uint64_t)file_stat.st_size;
        if (fread(config->rc_twopass_stats_in.buf, 1, file_stat.st_size,
            config->input_stat_file) != (size_t)file_stat.st_size) {
            return EB_FALSE;
        }
    }
    return config->rc_twopass_stats_in.buf != NULL;
}

/* set two passes stats information to EbConfig
 */
EbErrorType set_two_passes_stats(EbConfig *config, EncodePass pass,
    const SvtAv1FixedBuf* rc_twopass_stats_in, uint32_t channel_number)
{
    switch (pass) {
        case ENCODE_SINGLE_PASS: {
            const char* stats = config->stats ? config->stats : "svtav1_2pass.log";
            if (config->pass == 1) {
                if (!fopen_and_lock(&config->output_stat_file, stats, EB_TRUE)) {
                    fprintf(config->error_log_file,
                        "Error instance %u: can't open stats file %s for write \n",
                        channel_number + 1, stats);
                    return EB_ErrorBadParameter;
                }
                config->rc_firstpass_stats_out = EB_TRUE;
            } else if (config->pass == 2) {
                if (!fopen_and_lock(&config->input_stat_file, stats, EB_FALSE)) {
                    fprintf(config->error_log_file,
                        "Error instance %u: can't read stats file %s for read\n",
                        channel_number + 1, stats);
                    return EB_ErrorBadParameter;
                }
                if (!load_twopass_stats_in(config)) {
                    fprintf(config->error_log_file,
                        "Error instance %u: can't load file %s\n",
                        channel_number + 1, stats);
                    return EB_ErrorBadParameter;
                }
            }
            break;
        }
        case ENCODE_FIRST_PASS: {
            // for combined two passes,
            // we only ouptut first pass stats when user explicitly set the --stats
            if (config->stats) {
                if (!fopen_and_lock(&config->output_stat_file, config->stats, EB_TRUE)) {
                    fprintf(config->error_log_file,
                        "Error instance %u: can't open stats file %s for write \n",
                        channel_number + 1, config->stats);
                    return EB_ErrorBadParameter;
                }
            }
            config->rc_firstpass_stats_out = EB_TRUE;
            break;
        }
        case ENCODE_LAST_PASS: {
            if (!rc_twopass_stats_in->sz) {
                fprintf(config->error_log_file,
                        "Error instance %u: bug, combined 2passes need stats in for second pass \n",
                        channel_number + 1);
                return EB_ErrorBadParameter;
            }
            config->rc_twopass_stats_in = *rc_twopass_stats_in;
            break;
        }
        default: {
            assert(0);
            break;
        }
    }
    return EB_ErrorNone;
}

/******************************************
* Verify Settings
******************************************/
static EbErrorType verify_settings(EbConfig *config, uint32_t channel_number) {
    EbErrorType return_error = EB_ErrorNone;

    // Check Input File
    if (config->input_file == (FILE *)NULL) {
        fprintf(
            config->error_log_file, "Error instance %u: Invalid Input File\n", channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->frames_to_be_encoded <= -1) {
        fprintf(config->error_log_file,
                "Error instance %u: FrameToBeEncoded must be greater than 0\n",
                channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->buffered_input < -1) {
        fprintf(config->error_log_file,
                "Error instance %u: Invalid buffered_input. buffered_input must greater or equal "
                "to -1\n",
                channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->buffered_input > config->frames_to_be_encoded) {
        fprintf(config->error_log_file,
                "Error instance %u: Invalid buffered_input. buffered_input must be less or equal "
                "to the number of frames to be encoded\n",
                channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->use_qp_file == EB_TRUE && config->qp_file == NULL) {
        fprintf(config->error_log_file,
                "Error instance %u: Could not find QP file, UseQpFile is set to 1\n",
                channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->encoder_color_format != 1) {
        fprintf(config->error_log_file,
                "Error instance %u: Only support 420 now \n",
                channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }
    if (config->injector > 1) {
        fprintf(config->error_log_file,
                "Error Instance %u: Invalid injector [0 - 1]\n",
                channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }

    if (config->injector_frame_rate > (240 << 16) && config->injector) {
        fprintf(config->error_log_file,
                "Error Instance %u: The maximum allowed injector_frame_rate is 240 fps\n",
                channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }
    // Check that the injector frame_rate is non-zero
    if (!config->injector_frame_rate && config->injector) {
        fprintf(config->error_log_file,
                "Error Instance %u: The injector frame rate should be greater than 0 fps \n",
                channel_number + 1);
        return_error = EB_ErrorBadParameter;
    }

    // target_socket
    if (config->target_socket != -1 && config->target_socket != 0 && config->target_socket != 1) {
        fprintf(config->error_log_file,
                "Error instance %u: Invalid target_socket [-1 - 1], your input: %d\n",
                channel_number + 1,
                config->target_socket);
        return_error = EB_ErrorBadParameter;
    }

    if (config->input_stat_file && config->output_stat_file) {
        fprintf(config->error_log_file,
                "Error instance %u: do not set input_stat_file and output_stat_file at same time\n",
                channel_number + 1);
        return EB_ErrorBadParameter;
    }
    int pass = config->pass;

    if (pass != 2 && pass != 1 && pass != DEFAULT) {
        fprintf(config->error_log_file,
                "Error instance %u: %d pass encode is not supported\n",
                channel_number + 1, config->pass);
        return EB_ErrorBadParameter;
    }

    if (pass != DEFAULT  && (config->input_stat_file || config->output_stat_file)) {
        fprintf(config->error_log_file,
                "Error instance %u: --pass can't work with -input-stat-file or -output-stat-file \n",
                channel_number + 1);
        return EB_ErrorBadParameter;
    }

    if ((pass != DEFAULT || config->input_stat_file || config->output_stat_file) && channel_number > 0) {
        fprintf(config->error_log_file,
                "Error instance %u: 2 pass encode for multi instance is not supported\n",
                channel_number + 1);
        return EB_ErrorBadParameter;
    }
    if (pass != DEFAULT || config->input_stat_file || config->output_stat_file) {
        if (config->hierarchical_levels != 4) {
            fprintf(config->error_log_file,
                "Error instance %u: 2 pass encode for hierarchical_levels %u is not supported\n",
                channel_number + 1, config->hierarchical_levels);
            return EB_ErrorBadParameter;
        }
        if (config->enable_overlays) {
            fprintf(config->error_log_file,
                "Error instance %u: 2 pass encode for overlays is not supported\n",
                channel_number + 1);
            return EB_ErrorBadParameter;
        }
        if (config->intra_refresh_type != 2) {
            fprintf(config->error_log_file,
                "Error instance %u: 2 pass encode for intra_refresh_type %u is not supported\n",
                channel_number + 1,config->intra_refresh_type);
            return EB_ErrorBadParameter;
        }
    }
    return return_error;
}

/******************************************
 * Find Token for multiple inputs
 ******************************************/
int32_t find_token_multiple_inputs(int32_t argc, char *const argv[], const char *token,
                                   char **configStr) {
    int32_t return_error = -1;
    int32_t done         = 0;
    while ((argc > 0) && (return_error != 0)) {
        return_error = strcmp(argv[--argc], token);
        if (return_error == 0) {
            for (unsigned count = 0; count < MAX_CHANNEL_NUMBER; ++count) {
                if (done == 0) {
                    if (argv[argc + count + 1]) {
                        if (strtoul(argv[argc + count + 1], NULL, 0) != 0 ||
                            strcmp(argv[argc + count + 1], "0") == 0) {
                            strcpy_s(
                                configStr[count], COMMAND_LINE_MAX_SIZE, argv[argc + count + 1]);
                        } else if (argv[argc + count + 1][0] != '-') {
                            strcpy_s(
                                configStr[count], COMMAND_LINE_MAX_SIZE, argv[argc + count + 1]);
                        } else {
                            strcpy_s(configStr[count], COMMAND_LINE_MAX_SIZE, " ");
                            done = 1;
                        }
                    } else {
                        strcpy_s(configStr[count], COMMAND_LINE_MAX_SIZE, " ");
                        done = 1;
                        //return return_error;
                    }
                } else
                    strcpy_s(configStr[count], COMMAND_LINE_MAX_SIZE, " ");
            }
        }
    }

    return return_error;
}

static int check_long(ConfigEntry cfg_entry, ConfigEntry cfg_entry_next) {
    return cfg_entry_next.name ? !strcmp(cfg_entry.name, cfg_entry_next.name) : 0;
}

uint32_t get_help(int32_t argc, char *const argv[]) {
    char config_string[COMMAND_LINE_MAX_SIZE];
    if (find_token(argc, argv, HELP_TOKEN, config_string) &&
        find_token(argc, argv, HELP_LONG_TOKEN, config_string))
        return 0;

    printf(
        "Usage: SvtAv1EncApp <options> -b dst_filename -i src_filename\n\n"
        "Examples:\n"
        "Two passes encode:\n"
        "    SvtAv1EncApp <--stats svtav1_2pass.log> --pass 1 -b dst_filename -i src_filename\n"
        "    SvtAv1EncApp <--stats svtav1_2pass.log> --pass 2 -b dst_filename -i src_filename\n"
        "Or a combined cli:\n"
        "    SvtAv1EncApp <--stats svtav1_2pass.log> --passes 2 -b dst_filename -i src_filename\n"
        "\nOptions:\n");
    for (ConfigEntry *options_token_index = config_entry_options; options_token_index->token;
         ++options_token_index) {
        // this only works if short and long token are one after another
        switch (check_long(*options_token_index, options_token_index[1])) {
        case 1:
            printf("  %s, %-25s    %-25s\n",
                   options_token_index->token,
                   options_token_index[1].token,
                   options_token_index->name);
            ++options_token_index;
            break;
        default:
            printf(options_token_index->token[1] == '-' ? "      %-25s    %-25s\n"
                                                        : "      -%-25s   %-25s\n",
                   options_token_index->token,
                   options_token_index->name);
        }
    }
    printf("\nEncoder Global Options:\n");
    for (ConfigEntry *global_options_token_index = config_entry_global_options;
         global_options_token_index->token;
         ++global_options_token_index) {
        switch (check_long(*global_options_token_index, global_options_token_index[1])) {
        case 1:
            printf("  %s, %-25s    %-25s\n",
                   global_options_token_index->token,
                   global_options_token_index[1].token,
                   global_options_token_index->name);
            ++global_options_token_index;
            break;
        default:
            printf(global_options_token_index->token[1] == '-' ? "      %-25s    %-25s\n"
                                                               : "      -%-25s   %-25s\n",
                   global_options_token_index->token,
                   global_options_token_index->name);
        }
    }
    printf("\nRate Control Options:\n");
    for (ConfigEntry *rc_token_index = config_entry_rc; rc_token_index->token; ++rc_token_index) {
        switch (check_long(*rc_token_index, rc_token_index[1])) {
        case 1:
            printf("  %s, %-25s    %-25s\n",
                   rc_token_index->token,
                   rc_token_index[1].token,
                   rc_token_index->name);
            ++rc_token_index;
            break;
        default:
            printf(rc_token_index->token[1] == '-' ? "      %-25s    %-25s\n"
                                                   : "      -%-25s   %-25s\n",
                   rc_token_index->token,
                   rc_token_index->name);
        }
    }
    printf("\nTwopass Options:\n");
    for (ConfigEntry *two_p_token_index = config_entry_2p; two_p_token_index->token;
         ++two_p_token_index) {
        switch (check_long(*two_p_token_index, two_p_token_index[1])) {
        case 1:
            printf("  %s, %-25s    %-25s\n",
                   two_p_token_index->token,
                   two_p_token_index[1].token,
                   two_p_token_index->name);
            ++two_p_token_index;
            break;
        default:
            printf(two_p_token_index->token[1] == '-' ? "      %-25s    %-25s\n"
                                                      : "      -%-25s   %-25s\n",
                   two_p_token_index->token,
                   two_p_token_index->name);
        }
    }
    printf("\nKeyframe Placement Options:\n");
    for (ConfigEntry *kf_token_index = config_entry_intra_refresh; kf_token_index->token;
         ++kf_token_index) {
        switch (check_long(*kf_token_index, kf_token_index[1])) {
        case 1:
            printf("  %s, %-25s    %-25s\n",
                   kf_token_index->token,
                   kf_token_index[1].token,
                   kf_token_index->name);
            ++kf_token_index;
            break;
        default:
            printf(kf_token_index->token[1] == '-' ? "      %-25s    %-25s\n"
                                                   : "      -%-25s   %-25s\n",
                   kf_token_index->token,
                   kf_token_index->name);
        }
    }
    printf("\nAV1 Specific Options:\n");
    for (ConfigEntry *sp_token_index = config_entry_specific; sp_token_index->token;
         ++sp_token_index) {
        switch (check_long(*sp_token_index, sp_token_index[1])) {
        case 1:
            printf("  %s, %-25s    %-25s\n",
                   sp_token_index->token,
                   sp_token_index[1].token,
                   sp_token_index->name);
            ++sp_token_index;
            break;
        default:
            printf(sp_token_index->token[1] == '-' ? "      %-25s    %-25s\n"
                                                   : "      -%-25s   %-25s\n",
                   sp_token_index->token,
                   sp_token_index->name);
        }
    }
    return 1;
}

/******************************************************
* Get the number of channels and validate it with input
******************************************************/
uint32_t get_number_of_channels(int32_t argc, char *const argv[]) {
    char config_string[COMMAND_LINE_MAX_SIZE];
    if (find_token(argc, argv, CHANNEL_NUMBER_TOKEN, config_string) == 0) {
        // Set the input file
        uint32_t channel_number = strtol(config_string, NULL, 0);
        if ((channel_number > MAX_CHANNEL_NUMBER) || channel_number == 0) {
            fprintf(stderr,
                    "Error: The number of channels has to be within the range [1,%u]\n",
                    MAX_CHANNEL_NUMBER);
            return 0;
        }
        return channel_number;
    }
    return 1;
}

static EbBool check_two_pass_conflicts(int32_t argc, char *const argv[])
{
    char     config_string[COMMAND_LINE_MAX_SIZE];
    const char* conflicts[] = {
        PASS_TOKEN,
        INPUT_STAT_FILE_TOKEN,
        OUTPUT_STAT_FILE_TOKEN,
        NULL,
    };
    int i = 0;
    const char* token;
    while ((token = conflicts[i])) {
        if (find_token(argc, argv, token, config_string) == 0) {
            fprintf(stderr,
                "--passes 2 conflicts with %s\n", token);
            return EB_TRUE;
        }
        i++;
    }
    return EB_FALSE;
}

uint32_t get_passes(int32_t argc, char *const argv[], EncodePass pass[MAX_ENCODE_PASS]) {
    char     config_string[COMMAND_LINE_MAX_SIZE];
    uint32_t passes;
    if (find_token(argc, argv, PASSES_TOKEN, config_string) != 0) {
        pass[0] = ENCODE_SINGLE_PASS;
        return 1;
    }
    passes = strtol(config_string, NULL, 0);
    if (passes == 0 || passes > 2) {
        fprintf(stderr,
            "Error: The number of passes has to be within the range [1,%u]\n",
            (uint32_t)MAX_ENCODE_PASS);
        return 0;
    }
    if (passes == 1) {
        pass[0] = ENCODE_SINGLE_PASS;
        return 1;
    }

    int preset = MAX_ENC_PRESET;
    if (find_token(argc, argv, PRESET_TOKEN, config_string) == 0
        || find_token(argc, argv, ENCMODE_TOKEN, config_string) == 0) {
        preset = strtol(config_string, NULL, 0);
    }
    int rc_mode = 0;
#if FIX_RC_TOKEN
    if (find_token(argc, argv, RATE_CONTROL_ENABLE_TOKEN, config_string) == 0 ||
        find_token(argc, argv, "--rc", config_string) == 0 )
#else
    if (find_token(argc, argv, RATE_CONTROL_ENABLE_TOKEN, config_string) == 0 )
#endif
        rc_mode = strtol(config_string, NULL, 0);

    if (preset > 3 && rc_mode == 0) {
        fprintf(stderr,
            "\nWarn: --passes 2 CRF for preset > 3 is not supported yet, force single pass\n\n");
        pass[0] = ENCODE_SINGLE_PASS;
        return 1;
    }
    if (check_two_pass_conflicts(argc, argv))
        return 0;

    pass[0] = ENCODE_FIRST_PASS;
    pass[1] = ENCODE_LAST_PASS;
    return 2;
}

void mark_token_as_read(const char *token, char *cmd_copy[], int32_t *cmd_token_cnt) {
    int32_t cmd_copy_index;
    for (cmd_copy_index = 0; cmd_copy_index < *(cmd_token_cnt); ++cmd_copy_index) {
        if (!strcmp(cmd_copy[cmd_copy_index], token))
            cmd_copy[cmd_copy_index] = cmd_copy[--(*cmd_token_cnt)];
    }
}

EbBool is_negative_number(const char *string) {
    int32_t length = (int32_t)strlen(string);
    int32_t index  = 0;
    if (string[0] != '-') return EB_FALSE;
    for (index = 1; index < length; index++) {
        if (string[index] < '0' || string[index] > '9') return EB_FALSE;
    }
    return EB_TRUE;
}

// Computes the number of frames in the input file
int32_t compute_frames_to_be_encoded(EbConfig *config) {
    uint64_t file_size   = 0;
    int32_t  frame_count = 0;
    uint32_t frame_size;

    // Pipes contain data streams whose end we cannot know before we reach it.
    // For pipes, we leave it up to the eof logic to detect how many frames to eventually encode.
    if (config->input_file == stdin || config->input_file_is_fifo) return -1;

    if (config->input_file) {
        uint64_t curr_loc = ftello(config->input_file); // get current fp location
        fseeko(config->input_file, 0L, SEEK_END);
        file_size = ftello(config->input_file);
        fseeko(config->input_file, curr_loc, SEEK_SET); // seek back to that location
    }

    frame_size = config->input_padded_width * config->input_padded_height; // Luma
    frame_size += 2 * (frame_size >> (3 - config->encoder_color_format)); // Add Chroma
    frame_size = frame_size << ((config->encoder_bit_depth == 10) ? 1 : 0);

    if (frame_size == 0) return -1;

    if (config->encoder_bit_depth == 10 && config->compressed_ten_bit_format == 1)
        frame_count = (int32_t)(2 * ((double)file_size / frame_size) / 1.25);
    else
        frame_count = (int32_t)(file_size / frame_size);

    if (frame_count == 0) return -1;

    return frame_count;
}

/**********************************
* Parse Pred Struct File
**********************************/
static int32_t parse_pred_struct_file(EbConfig *config, char *buffer, int32_t size) {
    uint32_t argc;
    char *   argv[CONFIG_FILE_MAX_ARG_COUNT];
    uint32_t arg_len[CONFIG_FILE_MAX_ARG_COUNT];

    char var_name[CONFIG_FILE_MAX_VAR_LEN];
    char var_value[CONFIG_FILE_MAX_ARG_COUNT][CONFIG_FILE_MAX_VAR_LEN];

    uint32_t value_index;

    uint32_t comment_section_flag = 0;
    uint32_t new_line_flag        = 0;
    int32_t  entry_num            = 0;
    int32_t  display_order = 0, num_ref_list0 = 0, num_ref_list1 = 0;
    int32_t  idx_ref_list0 = 0, idx_ref_list1 = 0;

    // Keep looping until we process the entire file
    while (size--) {
        comment_section_flag =
            ((*buffer == CONFIG_FILE_COMMENT_CHAR) || (comment_section_flag != 0))
                ? 1
                : comment_section_flag;

        // At the beginning of each line
        if ((new_line_flag == 1) && (comment_section_flag == 0)) {
            // Do an argc/argv split for the line
            line_split(&argc, argv, arg_len, buffer);

            if ((argc > 2) && (*argv[1] == CONFIG_FILE_VALUE_SPLIT)) {
                // ***NOTE - We're assuming that the variable name is the first arg and
                // the variable value is the third arg.

                // Cap the length of the variable name
                arg_len[0] = (arg_len[0] > CONFIG_FILE_MAX_VAR_LEN - 1)
                                 ? CONFIG_FILE_MAX_VAR_LEN - 1
                                 : arg_len[0];
                // Copy the variable name
                strncpy_s(var_name, CONFIG_FILE_MAX_VAR_LEN, argv[0], arg_len[0]);
                // Null terminate the variable name
                var_name[arg_len[0]] = CONFIG_FILE_NULL_CHAR;
                if (strcmp(var_name, "PredStructEntry")) { continue; }

                ++entry_num;
                idx_ref_list0 = idx_ref_list1 = num_ref_list0 = num_ref_list1 = 0;

                for (value_index = 0;
                     (value_index < CONFIG_FILE_MAX_ARG_COUNT - 2) && (value_index < (argc - 2));
                     ++value_index) {
                    // Cap the length of the variable
                    arg_len[value_index + 2] =
                        (arg_len[value_index + 2] > CONFIG_FILE_MAX_VAR_LEN - 1)
                            ? CONFIG_FILE_MAX_VAR_LEN - 1
                            : arg_len[value_index + 2];
                    // Copy the variable name
                    strncpy_s(var_value[value_index],
                               CONFIG_FILE_MAX_VAR_LEN,
                               argv[value_index + 2],
                               arg_len[value_index + 2]);
                    // Null terminate the variable name
                    var_value[value_index][arg_len[value_index + 2]] = CONFIG_FILE_NULL_CHAR;

                    switch (value_index) {
                    case 0:
                        display_order = strtoul(var_value[value_index], NULL, 0);
                        if (display_order >= (1 << (MAX_HIERARCHICAL_LEVEL - 1))) { return -1; }
                        break;
                    case 1:
                        config->pred_struct[display_order].decode_order =
                            strtoul(var_value[value_index], NULL, 0);
                        if (config->pred_struct[display_order].decode_order >=
                            (1 << (MAX_HIERARCHICAL_LEVEL - 1))) {
                            return -1;
                        }
                        break;
                    case 2:
                        config->pred_struct[display_order].temporal_layer_index =
                            strtoul(var_value[value_index], NULL, 0);
                        break;
                    case 3: num_ref_list0 = strtoul(var_value[value_index], NULL, 0); break;
                    case 4: num_ref_list1 = strtoul(var_value[value_index], NULL, 0); break;
                    default:
                        if (idx_ref_list0 < num_ref_list0) {
                            if (idx_ref_list0 < REF_LIST_MAX_DEPTH) {
                                config->pred_struct[display_order].ref_list0[idx_ref_list0] =
                                    strtoul(var_value[value_index], NULL, 0);
                            }
                            ++idx_ref_list0;
                        } else if (idx_ref_list1 < num_ref_list1) {
                            if (idx_ref_list1 < REF_LIST_MAX_DEPTH - 1) {
                                config->pred_struct[display_order].ref_list1[idx_ref_list1] =
                                    strtoul(var_value[value_index], NULL, 0);
                            }
                            ++idx_ref_list1;
                        }
                        break;
                    }
                }
                for (; num_ref_list0 < REF_LIST_MAX_DEPTH; ++num_ref_list0) {
                    config->pred_struct[display_order].ref_list0[num_ref_list0] = 0;
                }
                for (; num_ref_list1 < REF_LIST_MAX_DEPTH - 1; ++num_ref_list1) {
                    config->pred_struct[display_order].ref_list1[num_ref_list1] = 0;
                }
            }
        }

        comment_section_flag = (*buffer == CONFIG_FILE_NEWLINE_CHAR) ? 0 : comment_section_flag;
        new_line_flag        = (*buffer == CONFIG_FILE_NEWLINE_CHAR) ? 1 : 0;
        ++buffer;
    }
    config->manual_pred_struct_entry_num = entry_num;

    return 0;
}

/**********************************
* Read Prediction Structure File
**********************************/
static int32_t read_pred_struct_file(EbConfig *config, char *PredStructPath,
                                     uint32_t instance_idx) {
    int32_t return_error = 0;

    FOPEN(config->input_pred_struct_file, PredStructPath, "rb");

    if (config->input_pred_struct_file != (FILE *)NULL) {
        int32_t config_file_size   = find_file_size(config->input_pred_struct_file);
        char *  config_file_buffer = (char *)malloc(config_file_size);

        if (config_file_buffer != (char *)NULL) {
            int32_t result_size = (int32_t)fread(
                config_file_buffer, 1, config_file_size, config->input_pred_struct_file);

            if (result_size == config_file_size) {
                parse_pred_struct_file(config, config_file_buffer, config_file_size);
            } else {
                fprintf(stderr, "Error channel %u: File Read Failed\n", instance_idx + 1);
                return_error = -1;
            }
        } else {
            fprintf(stderr, "Error channel %u: Memory Allocation Failed\n", instance_idx + 1);
            return_error = -1;
        }

        free(config_file_buffer);
        fclose(config->input_pred_struct_file);
        config->input_pred_struct_file = (FILE *)NULL;
    } else {
        fprintf(stderr,
                "Error channel %u: Couldn't open Manual Prediction Structure File: %s\n",
                instance_idx + 1,
                PredStructPath);
        return_error = -1;
    }

    return return_error;
}

EbErrorType handle_short_tokens(char *string) {
    char *short_token_only = string + 2;
    if (*(string + 3) == '\0') {
        fprintf(stderr,
                "warning: please use -%s instead of --%s\n",
                short_token_only,
                short_token_only);
        return EB_ErrorBadParameter;
    }
    return EB_ErrorNone;
}

const char *handle_warnings(const char *token, char *print_message, uint8_t double_dash_token) {
    const char *linked_token = "";

    if (strcmp(token, ENCMODE_TOKEN) == 0) linked_token = PRESET_TOKEN;
    if (strcmp(token, ENCODER_BIT_DEPTH) == 0) linked_token = INPUT_DEPTH_TOKEN;
    if (strcmp(token, INTRA_PERIOD_TOKEN) == 0) linked_token = KEYINT_TOKEN;
    if (strcmp(token, QP_FILE_TOKEN) == 0) linked_token = QP_FILE_NEW_TOKEN;
    if (strcmp(token, LOOK_AHEAD_DIST_TOKEN) == 0) linked_token = LOOKAHEAD_NEW_TOKEN;

    if (strcmp(token, STAT_REPORT_TOKEN) == 0) linked_token = STAT_REPORT_NEW_TOKEN;
    if (strcmp(token, RESTORATION_ENABLE_TOKEN) == 0)
        linked_token = RESTORATION_ENABLE_NEW_TOKEN;
    if (strcmp(token, INTER_INTRA_COMPOUND_TOKEN) == 0)
        linked_token = INTER_INTRA_COMPOUND_NEW_TOKEN;
    if (strcmp(token, MFMV_ENABLE_TOKEN) == 0) linked_token = MFMV_ENABLE_NEW_TOKEN;
    if (strcmp(token, REDUNDANT_BLK_TOKEN) == 0) linked_token = REDUNDANT_BLK_NEW_TOKEN;
    if (strcmp(token, SPATIAL_SSE_FL_TOKEN) == 0) linked_token = SPATIAL_SSE_FL_NEW_TOKEN;
    if (strcmp(token, OVR_BNDRY_BLK_TOKEN) == 0) linked_token = OVR_BNDRY_BLK_NEW_TOKEN;
    if (strcmp(token, NEW_NEAREST_COMB_INJECT_TOKEN) == 0)
        linked_token = NEW_NEAREST_COMB_INJECT_NEW_TOKEN;
    if (strcmp(token, NSQ_TABLE_TOKEN) == 0) linked_token = NSQ_TABLE_NEW_TOKEN;
    if (strcmp(token, FRAME_END_CDF_UPDATE_TOKEN) == 0)
        linked_token = FRAME_END_CDF_UPDATE_NEW_TOKEN;
    if (strcmp(token, LOCAL_WARPED_ENABLE_TOKEN) == 0)
        linked_token = LOCAL_WARPED_ENABLE_NEW_TOKEN;
    if (strcmp(token, GLOBAL_MOTION_ENABLE_TOKEN) == 0)
        linked_token = GLOBAL_MOTION_ENABLE_NEW_TOKEN;
    if (strcmp(token, RDOQ_TOKEN) == 0) linked_token = RDOQ_NEW_TOKEN;
    if (strcmp(token, FILTER_INTRA_TOKEN) == 0) linked_token = FILTER_INTRA_NEW_TOKEN;
    if (strcmp(token, HDR_INPUT_TOKEN) == 0) linked_token = HDR_INPUT_NEW_TOKEN;
    if (strcmp(token, ADAPTIVE_QP_ENABLE_TOKEN) == 0)
        linked_token = ADAPTIVE_QP_ENABLE_NEW_TOKEN;

    if (strcmp(token, DISABLE_CFL_TOKEN) == 0) linked_token = DISABLE_CFL_NEW_TOKEN;
    if (strcmp(token, INTRA_EDGE_FILTER_TOKEN) == 0) linked_token = INTRA_EDGE_FILTER_NEW_TOKEN;
    if (strcmp(token, INTRA_ANGLE_DELTA_TOKEN) == 0) linked_token = INTRA_ANGLE_DELTA_NEW_TOKEN;
    if (strcmp(token, PAETH_TOKEN) == 0) linked_token = PAETH_NEW_TOKEN;
    if (strcmp(token, SMOOTH_TOKEN) == 0) linked_token = SMOOTH_NEW_TOKEN;
    if (strcmp(token, PIC_BASED_RATE_EST_TOKEN) == 0) linked_token = PIC_BASED_RATE_EST_NEW_TOKEN;

    if (strnlen_s(linked_token, WARNING_LENGTH) > 1) {
        const char *message_str = " will be deprecated soon, please use ";
        size_t offset;
        strcpy_s(print_message, WARNING_LENGTH, token);
        offset = strnlen_s(print_message, WARNING_LENGTH);
        strcpy_s(print_message + offset, WARNING_LENGTH - offset, message_str);
        offset = strnlen_s(print_message, WARNING_LENGTH);
        strcpy_s(print_message + offset, WARNING_LENGTH - offset, linked_token);
        return print_message;
    } else if (double_dash_token == 0) {
       const char *message_str = " will be deprecated soon, please use -";
        size_t offset;
        strcpy_s(print_message, WARNING_LENGTH, token);
        offset = strnlen_s(print_message, WARNING_LENGTH);
        strcpy_s(print_message + offset, WARNING_LENGTH - offset, message_str);
        offset = strnlen_s(print_message, WARNING_LENGTH);
        strcpy_s(print_message + offset, WARNING_LENGTH - offset, token);
        return print_message;
    }
    return "";
}

/******************************************
* Read Command Line
******************************************/
EbErrorType read_command_line(int32_t argc, char *const argv[], EncChannel *channels,
                              uint32_t num_channels,
                              char *warning_str[WARNING_LENGTH]) {
    EbErrorType return_error = EB_ErrorNone;
    char        config_string[COMMAND_LINE_MAX_SIZE]; // for one input options
    char *      config_strings[MAX_CHANNEL_NUMBER]; // for multiple input options
    char *      cmd_copy[MAX_NUM_TOKENS]; // keep track of extra tokens
    uint32_t    index         = 0;
    int32_t     cmd_token_cnt = 0; // total number of tokens
    int32_t     token_index   = -1;
    int32_t     ret_y4m;

    for (index = 0; index < MAX_CHANNEL_NUMBER; ++index)
        config_strings[index] = (char *)malloc(sizeof(char) * COMMAND_LINE_MAX_SIZE);
    // Copy tokens (except for CHANNEL_NUMBER_TOKEN and PASSES_TOKEN ) into a temp token buffer hosting all tokens that are passed through the command line
    size_t len = COMMAND_LINE_MAX_SIZE;
    for (token_index = 0; token_index < argc; ++token_index) {
        if ((argv[token_index][0] == '-') &&
            strncmp(argv[token_index], CHANNEL_NUMBER_TOKEN, len) &&
            strncmp(argv[token_index], PASSES_TOKEN, len) &&
            !is_negative_number(argv[token_index]))
            cmd_copy[cmd_token_cnt++] = argv[token_index];
    }

    /***************************************************************************************************/
    /****************  Find configuration files tokens and call respective functions  ******************/
    /***************************************************************************************************/

    // Find the Config File Path in the command line
    if (find_token_multiple_inputs(argc, argv, CONFIG_FILE_TOKEN, config_strings) == 0) {
        mark_token_as_read(CONFIG_FILE_TOKEN, cmd_copy, &cmd_token_cnt);
        // Parse the config file
        for (index = 0; index < num_channels; ++index) {
            EncChannel* c = channels + index;
            c->return_error =
                (EbErrorType)read_config_file(c->config, config_strings[index], index);
            return_error = (EbErrorType)(return_error & c->return_error);
        }
    } else if (find_token_multiple_inputs(argc, argv, CONFIG_FILE_LONG_TOKEN, config_strings) ==
               0) {
        mark_token_as_read(CONFIG_FILE_LONG_TOKEN, cmd_copy, &cmd_token_cnt);
        // Parse the config file
        for (index = 0; index < num_channels; ++index) {
            EncChannel* c = channels + index;
            c->return_error =
                (EbErrorType)read_config_file(c->config, config_strings[index], index);
            return_error = (EbErrorType)(return_error & c->return_error);
        }
    } else {
        if (find_token(argc, argv, CONFIG_FILE_TOKEN, config_string) == 0) {
            fprintf(stderr, "Error: Config File Token Not Found\n");
            return EB_ErrorBadParameter;
        } else
            return_error = EB_ErrorNone;
    }

    /***************************************************************************************************/
    /***********   Find SINGLE_INPUT configuration parameter tokens and call respective functions  **********/
    /***************************************************************************************************/
    token_index                     = -1;
    EbErrorType return_result_error = EB_ErrorNone;
    uint32_t    warning_index       = -1;
    // Parse command line for tokens
    while (config_entry[++token_index].name != NULL) {
        if (config_entry[token_index].type == SINGLE_INPUT) {
            char message[WARNING_LENGTH] = "";
            // concat strings with '-'
            char concat_str[WARNING_LENGTH] = "-";
            strcpy_s(concat_str + 1, sizeof(concat_str) - 1, config_entry[token_index].token);
            if (find_token_multiple_inputs(
                    argc, argv, config_entry[token_index].token, config_strings) == 0) {
                //Warning for one dash
                if (*(config_entry[token_index].token + 2) != '\0' && // check for small token
                    *(config_entry[token_index].token + 1) != '-' && // check for --
                    return_result_error != EB_ErrorBadParameter) {
                    handle_warnings(config_entry[token_index].token, message, 0);
                    strcpy_s(warning_str[++warning_index], WARNING_LENGTH, message);
                }
                // When a token is found mark it as found in the temp token buffer
                mark_token_as_read(config_entry[token_index].token, cmd_copy, &cmd_token_cnt);

                // Fill up the values corresponding to each channel
                for (index = 0; index < num_channels; ++index) {
                    EncChannel* c = channels + index;
                    if (strcmp(config_strings[index], " "))
                        (*config_entry[token_index].scf)(config_strings[index], c->config);
                    else
                        break;
                }
            } else if (find_token_multiple_inputs(argc, argv, concat_str, config_strings) ==
                       0) { // handle double dash
                handle_warnings(config_entry[token_index].token, message, 1);
                // handle warnings for new tokens
                if (strnlen_s(message, sizeof(message) > 1)) {
                    char double_dash_warning[WARNING_LENGTH] = "-";
                    strcpy_s(double_dash_warning + 1,
                              sizeof(double_dash_warning), message);
                    strcpy_s(warning_str[++warning_index], WARNING_LENGTH, double_dash_warning);
                    warning_index++;
                }
                return_result_error = handle_short_tokens(concat_str);
                // When a token is found mark it as found in the temp token buffer
                mark_token_as_read(concat_str, cmd_copy, &cmd_token_cnt);

                // Fill up the values corresponding to each channel
                for (index = 0; index < num_channels; ++index) {
                    EncChannel* c = channels + index;
                    if (strcmp(config_strings[index], " "))
                        (*config_entry[token_index].scf)(config_strings[index], c->config);
                    else
                        break;
                }
            }
        }
    }
    if (return_result_error == EB_ErrorBadParameter) return EB_ErrorBadParameter;

    /***************************************************************************************************/
    /********************** Parse parameters from input file if in y4m format **************************/
    /********************** overriding config file and command line inputs    **************************/
    /***************************************************************************************************/

    for (index = 0; index < num_channels; ++index) {
        EncChannel* c = channels + index;
        if (c->config->y4m_input == EB_TRUE) {
            ret_y4m = read_y4m_header(c->config);
            if (ret_y4m == EB_ErrorBadParameter) {
                fprintf(stderr, "Error found when reading the y4m file parameters.\n");
                return EB_ErrorBadParameter;
            }
        }
    }
    /***************************************************************************************************/
    /*******************************   Parse manual prediction structure  ******************************/
    /***************************************************************************************************/
    for (index = 0; index < num_channels; ++index) {
        EncChannel* c = channels + index;
        EbConfig* config = c->config;
        if (config->enable_manual_pred_struct == EB_TRUE) {
            c->return_error = (EbErrorType)read_pred_struct_file(
                config, config->input_pred_struct_filename, index);
            return_error = (EbErrorType)(return_error & c->return_error);
        }
    }
    /***************************************************************************************************/
    /**************************************   Verify configuration parameters   ************************/
    /***************************************************************************************************/
    // Verify the config values
    if (return_error == 0) {
        return_error = EB_ErrorBadParameter;
        for (index = 0; index < num_channels; ++index) {
            EncChannel* c = channels + index;
            if (c->return_error == EB_ErrorNone) {
                EbConfig* config = c->config;
                c->return_error = verify_settings(config, index);

                // Assuming no errors, add padding to width and height
                if (c->return_error == EB_ErrorNone) {
                    config->input_padded_width = config->source_width;
                    config->input_padded_height = config->source_height;
                }

                // Assuming no errors, set the frames to be encoded to the number of frames in the input yuv
                if (c->return_error == EB_ErrorNone &&
                    config->frames_to_be_encoded == 0)
                    config->frames_to_be_encoded =
                        compute_frames_to_be_encoded(config);

                // For pipe input it is fine if we have -1 here (we will update on end of stream)
                if (config->frames_to_be_encoded == -1
                    && config->input_file != stdin
                    && !config->input_file_is_fifo) {
                    fprintf(config->error_log_file,
                            "Error instance %u: Input yuv does not contain enough frames \n",
                            index + 1);
                    c->return_error = EB_ErrorBadParameter;
                }

                // Force the injector latency mode, and injector frame rate when speed control is on
                if (c->return_error == EB_ErrorNone && config->speed_control_flag == 1)
                    config->injector = 1;
            }
            return_error = (EbErrorType)(return_error & c->return_error);
        }
    }

    // Print message for unprocessed tokens
    if (cmd_token_cnt > 0) {
        int32_t cmd_copy_index;
        fprintf(stderr, "Unprocessed tokens: ");
        for (cmd_copy_index = 0; cmd_copy_index < cmd_token_cnt; ++cmd_copy_index)
            fprintf(stderr, " %s ", cmd_copy[cmd_copy_index]);
        fprintf(stderr, "\n\n");
        return_error = EB_ErrorBadParameter;
    }

    for (index = 0; index < MAX_CHANNEL_NUMBER; ++index) free(config_strings[index]);
    return return_error;
}
