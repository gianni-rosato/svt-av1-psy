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
#include <io.h>
#else
#include <unistd.h>
#endif

/**********************************
 * Defines
 **********************************/
#define HELP_TOKEN "-help"
#define CHANNEL_NUMBER_TOKEN "-nch"
#define COMMAND_LINE_MAX_SIZE 2048
#define CONFIG_FILE_TOKEN "-c"
#define INPUT_FILE_TOKEN "-i"
#define OUTPUT_BITSTREAM_TOKEN "-b"
#define OUTPUT_RECON_TOKEN "-o"
#define ERROR_FILE_TOKEN "-errlog"
#define QP_FILE_TOKEN "-qp-file"
#define INPUT_STAT_FILE_TOKEN "-input-stat-file"
#define OUTPUT_STAT_FILE_TOKEN "-output-stat-file"
#define STAT_FILE_TOKEN "-stat-file"
#define WIDTH_TOKEN "-w"
#define HEIGHT_TOKEN "-h"
#define NUMBER_OF_PICTURES_TOKEN "-n"
#define BUFFERED_INPUT_TOKEN "-nb"
#define BASE_LAYER_SWITCH_MODE_TOKEN "-base-layer-switch-mode" // no Eval
#define QP_TOKEN "-q"
#define USE_QP_FILE_TOKEN "-use-q-file"
#define STAT_REPORT_TOKEN "-stat-report"
#define FRAME_RATE_TOKEN "-fps"
#define FRAME_RATE_NUMERATOR_TOKEN "-fps-num"
#define FRAME_RATE_DENOMINATOR_TOKEN "-fps-denom"
#define ENCODER_BIT_DEPTH "-bit-depth"
#define ENCODER_COLOR_FORMAT "-color-format"
#define INPUT_COMPRESSED_TEN_BIT_FORMAT "-compressed-ten-bit-format"
#define ENCMODE_TOKEN "-enc-mode"
#define ENCMODE2P_TOKEN "-enc-mode-2p"
#define HIERARCHICAL_LEVELS_TOKEN "-hierarchical-levels" // no Eval
#define PRED_STRUCT_TOKEN "-pred-struct"
#define INTRA_PERIOD_TOKEN "-intra-period"
#define PROFILE_TOKEN "-profile"
#define TIER_TOKEN "-tier"
#define LEVEL_TOKEN "-level"
#define LATENCY_MODE "-latency-mode" // no Eval
#define FILM_GRAIN_TOKEN "-film-grain"
#define INTRA_REFRESH_TYPE_TOKEN "-irefresh-type" // no Eval
#define LOOP_FILTER_DISABLE_TOKEN "-dlf"
#define RESTORATION_ENABLE_TOKEN "-restoration-filtering"
#define CLASS_12_TOKEN "-class-12"
#define EDGE_SKIP_ANGLE_INTRA_TOKEN "-intra-edge-skp"
#define INTER_INTRA_COMPOUND_TOKEN "-interintra-comp"
#define FRAC_SEARCH_64_TOKEN "-frac-search-64"
#define MFMV_ENABLE_TOKEN "-mfmv"
#define REDUNDANT_BLK_TOKEN "-redundant-blk"
#define TRELLIS_ENABLE_TOKEN "-trellis"
#define SPATIAL_SSE_FL_TOKEN "-spatial-sse-fl"
#define SUBPEL_TOKEN "-subpel"
#define OVR_BNDRY_BLK_TOKEN "-over-bndry-blk"
#define NEW_NEAREST_COMB_INJECT_TOKEN "-new-nrst-near-comb"
#define NX4_4XN_MV_INJECT_TOKEN "-nx4-4xn-mv-inject"
#define PRUNE_UNIPRED_ME_TOKEN "-prune-unipred-me"
#define PRUNE_REF_REC_PART_TOKEN "-prune-ref-rec-part"
#define NSQ_TABLE_TOKEN "-nsq-table-use"
#define FRAME_END_CDF_UPDATE_TOKEN "-framend-cdf-upd-mode"
#define LOCAL_WARPED_ENABLE_TOKEN "-local-warp"
#define GLOBAL_MOTION_ENABLE_TOKEN "-global-motion"
#define OBMC_TOKEN "-obmc"
#define RDOQ_TOKEN "-rdoq"
#define PRED_ME_TOKEN "-pred-me"
#define BIPRED_3x3_TOKEN "-bipred-3x3"
#define COMPOUND_LEVEL_TOKEN "-compound"
#define FILTER_INTRA_TOKEN "-filter-intra"
#define USE_DEFAULT_ME_HME_TOKEN "-use-default-me-hme"
#define HME_ENABLE_TOKEN "-hme"
#define HME_L0_ENABLE_TOKEN "-hme-l0"
#define HME_L1_ENABLE_TOKEN "-hme-l1"
#define HME_L2_ENABLE_TOKEN "-hme-l2"
#define EXT_BLOCK "-ext-block"
#define SEARCH_AREA_WIDTH_TOKEN "-search-w"
#define SEARCH_AREA_HEIGHT_TOKEN "-search-h"
#define NUM_HME_SEARCH_WIDTH_TOKEN "-num-hme-w"
#define NUM_HME_SEARCH_HEIGHT_TOKEN "-num-hme-h"
#define HME_SRCH_T_L0_WIDTH_TOKEN "-hme-tot-l0-w"
#define HME_SRCH_T_L0_HEIGHT_TOKEN "-hme-tot-l0-h"
#define HME_LEVEL0_WIDTH "-hme-l0-w"
#define HME_LEVEL0_HEIGHT "-hme-l0-h"
#define HME_LEVEL1_WIDTH "-hme-l1-w"
#define HME_LEVEL1_HEIGHT "-hme-l1-h"
#define HME_LEVEL2_WIDTH "-hme-l2-w"
#define HME_LEVEL2_HEIGHT "-hme-l2-h"
#define SCREEN_CONTENT_TOKEN "-scm"
// --- start: ALTREF_FILTERING_SUPPORT
#define ENABLE_ALTREFS "-enable-altrefs"
#define ALTREF_STRENGTH "-altref-strength"
#define ALTREF_NFRAMES "-altref-nframes"
#define ENABLE_OVERLAYS "-enable-overlays"
// --- end: ALTREF_FILTERING_SUPPORT
#define HBD_MD_ENABLE_TOKEN "-hbd-md"
#define PALETTE_TOKEN "-palette"
#define OLPD_REFINEMENT_TOKEN "-olpd-refinement"
#define HDR_INPUT_TOKEN "-hdr"
#define RATE_CONTROL_ENABLE_TOKEN "-rc"
#define TARGET_BIT_RATE_TOKEN "-tbr"
#define MAX_QP_TOKEN "-max-qp"
#define VBV_BUFSIZE_TOKEN "-vbv-bufsize"
#define MIN_QP_TOKEN "-min-qp"
#define ADAPTIVE_QP_ENABLE_TOKEN "-adaptive-quantization"
#define LOOK_AHEAD_DIST_TOKEN "-lad"
#define SUPER_BLOCK_SIZE_TOKEN "-sb-size"
#define TILE_ROW_TOKEN "-tile-rows"
#define TILE_COL_TOKEN "-tile-columns"

#define SQ_WEIGHT_TOKEN "-sqw"
#define ENABLE_AMP_TOKEN "-enable-amp"
#define CHROMA_MODE_TOKEN "-chroma-mode"

#define SCENE_CHANGE_DETECTION_TOKEN "-scd"
#define INJECTOR_TOKEN "-inj" // no Eval
#define INJECTOR_FRAMERATE_TOKEN "-inj-frm-rt" // no Eval
#define SPEED_CONTROL_TOKEN "-speed-ctrl"
#define ASM_TYPE_TOKEN "-asm"
#define THREAD_MGMNT "-lp"
#define UNPIN_LP1_TOKEN "-unpin-lp1"
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

#define MDS1_PRUNE_C_TH "-mds1p-class-th"
#define MDS1_PRUNE_S_TH "-mds1p-cand-th"
#define MDS2_PRUNE_C_TH "-mds2p-class-th"
#define MDS2_PRUNE_S_TH "-mds2p-cand-th"

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
    if (handle == INVALID_HANDLE_VALUE)
        return;
    cfg->input_file_is_fifo = GetFileType(handle) == FILE_TYPE_PIPE;
#else
    int         fd = fileno(cfg->input_file);
    struct stat statbuf;
    fstat(fd, &statbuf);
    cfg->input_file_is_fifo = S_ISFIFO(statbuf.st_mode);
#endif

    cfg->y4m_input = check_if_y4m(cfg);
};

static void set_cfg_stream_file(const char *value, EbConfig *cfg) {
    if (cfg->bitstream_file) { fclose(cfg->bitstream_file); }
    FOPEN(cfg->bitstream_file, value, "wb");
};
static void set_cfg_error_file(const char *value, EbConfig *cfg) {
    if (cfg->error_log_file) { fclose(cfg->error_log_file); }
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
static void set_input_stat_file(const char *value, EbConfig *cfg) {
    if (cfg->input_stat_file) { fclose(cfg->input_stat_file); }
    FOPEN(cfg->input_stat_file, value, "rb");
};
static void set_output_stat_file(const char *value, EbConfig *cfg) {
    if (cfg->output_stat_file) { fclose(cfg->output_stat_file); }
    FOPEN(cfg->output_stat_file, value, "wb");
};
static void set_snd_pass_enc_mode(const char *value, EbConfig *cfg) {
    cfg->snd_pass_enc_mode = (uint8_t)strtoul(value, NULL, 0);
};
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
static void set_frame_rate(const char *value, EbConfig *cfg) {
    cfg->frame_rate = strtoul(value, NULL, 0);
    if (cfg->frame_rate > 1000)
        cfg->frame_rate = cfg->frame_rate;
    else
        cfg->frame_rate = cfg->frame_rate << 16;
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
static void set_encoder_color_format(const char *value, EbConfig *cfg) {
    cfg->encoder_color_format = strtoul(value, NULL, 0);
}
static void set_compressed_ten_bit_format(const char *value, EbConfig *cfg) {
    cfg->compressed_ten_bit_format = strtoul(value, NULL, 0);
}
static void set_base_layer_switch_mode(const char *value, EbConfig *cfg) {
    cfg->base_layer_switch_mode = (EbBool)strtoul(value, NULL, 0);
};
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
    cfg->enable_warped_motion = (EbBool)strtoul(value, NULL, 0);
};
static void set_enable_global_motion_flag(const char *value, EbConfig *cfg) {
    cfg->enable_global_motion = (EbBool)strtoul(value, NULL, 0);
};
static void set_enable_restoration_filter_flag(const char *value, EbConfig *cfg) {
    cfg->enable_restoration_filtering = strtol(value, NULL, 0);
};
static void set_class_12_flag(const char *value, EbConfig *cfg) {
    cfg->combine_class_12 = strtol(value, NULL, 0);
};
static void set_edge_skip_angle_intra_flag(const char *value, EbConfig *cfg) {
    cfg->edge_skp_angle_intra = strtol(value, NULL, 0);
};
static void set_interintra_compound_flag(const char *value, EbConfig *cfg) {
    cfg->inter_intra_compound = strtol(value, NULL, 0);
};
static void set_fractional_search_64_flag(const char *value, EbConfig *cfg) {
    cfg->fract_search_64 = strtol(value, NULL, 0);
};
static void set_enable_mfmv_flag(const char *value, EbConfig *cfg) {
    cfg->enable_mfmv = strtol(value, NULL, 0);
};
static void set_enable_redundant_blk_flag(const char *value, EbConfig *cfg) {
    cfg->enable_redundant_blk = strtol(value, NULL, 0);
};
static void set_enable_trellis_flag(const char *value, EbConfig *cfg) {
    cfg->enable_trellis = strtol(value, NULL, 0);
};
static void set_spatial_sse_fl_flag(const char *value, EbConfig *cfg) {
    cfg->spatial_sse_fl = strtol(value, NULL, 0);
};
static void set_enable_sub_pel_flag(const char *value, EbConfig *cfg) {
    cfg->enable_subpel = strtol(value, NULL, 0);
};
static void set_over_bndry_blk_flag(const char *value, EbConfig *cfg) {
    cfg->over_bndry_blk = strtol(value, NULL, 0);
};
static void set_new_nearest_comb_inject_flag(const char *value, EbConfig *cfg) {
    cfg->new_nearest_comb_inject = strtol(value, NULL, 0);
};
static void set_nx4_4xn_parent_mv_inject_flag(const char *value, EbConfig *cfg) {
    cfg->nx4_4xn_parent_mv_inject = strtol(value, NULL, 0);
};
static void set_prune_unipred_me_flag(const char *value, EbConfig *cfg) {
    cfg->prune_unipred_me = strtol(value, NULL, 0);
};
static void set_prune_ref_rec_part_flag(const char *value, EbConfig *cfg) {
    cfg->prune_ref_rec_part = strtol(value, NULL, 0);
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
static void set_enable_obmc_flag(const char *value, EbConfig *cfg) {
    cfg->enable_obmc = (EbBool)strtoul(value, NULL, 0);
};
static void set_enable_rdoq_flag(const char *value, EbConfig *cfg) {
    cfg->enable_rdoq = strtol(value, NULL, 0);
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
static void set_enable_filter_intra_flag(const char *value, EbConfig *cfg) {
    cfg->enable_filter_intra = (EbBool)strtoul(value, NULL, 0);
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
};
static void set_look_ahead_distance(const char *value, EbConfig *cfg) {
    cfg->look_ahead_distance = strtoul(value, NULL, 0);
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
static void set_cfg_number_hme_search_region_in_width(const char *value, EbConfig *cfg) {
    cfg->number_hme_search_region_in_width = strtoul(value, NULL, 0);
};
static void set_cfg_number_hme_search_region_in_height(const char *value, EbConfig *cfg) {
    cfg->number_hme_search_region_in_height = strtoul(value, NULL, 0);
};
static void set_cfg_hme_level_0_total_search_area_width(const char *value, EbConfig *cfg) {
    cfg->hme_level0_total_search_area_width = strtoul(value, NULL, 0);
};
static void set_cfg_hme_level_0_total_search_area_height(const char *value, EbConfig *cfg) {
    cfg->hme_level0_total_search_area_height = strtoul(value, NULL, 0);
};
static void set_cfg_use_default_me_hme(const char *value, EbConfig *cfg) {
    cfg->use_default_me_hme = (EbBool)strtol(value, NULL, 0);
};
static void set_enable_ext_block_flag(const char *value, EbConfig *cfg) {
    cfg->ext_block_flag = (EbBool)strtoul(value, NULL, 0);
};
static void set_hme_level_0_search_area_in_width_array(const char *value, EbConfig *cfg) {
    cfg->hme_level0_search_area_in_width_array[cfg->hme_level0_column_index++] =
        strtoul(value, NULL, 0);
};
static void set_hme_level_0_search_area_in_height_array(const char *value, EbConfig *cfg) {
    cfg->hme_level0_search_area_in_height_array[cfg->hme_level0_row_index++] =
        strtoul(value, NULL, 0);
};
static void set_hme_level_1_search_area_in_width_array(const char *value, EbConfig *cfg) {
    cfg->hme_level1_search_area_in_width_array[cfg->hme_level1_column_index++] =
        strtoul(value, NULL, 0);
};
static void set_hme_level_1_search_area_in_height_array(const char *value, EbConfig *cfg) {
    cfg->hme_level1_search_area_in_height_array[cfg->hme_level1_row_index++] =
        strtoul(value, NULL, 0);
};
static void set_hme_level_2_search_area_in_width_array(const char *value, EbConfig *cfg) {
    cfg->hme_level2_search_area_in_width_array[cfg->hme_level2_column_index++] =
        strtoul(value, NULL, 0);
};
static void set_hme_level_2_search_area_in_height_array(const char *value, EbConfig *cfg) {
    cfg->hme_level2_search_area_in_height_array[cfg->hme_level2_row_index++] =
        strtoul(value, NULL, 0);
};
static void set_screen_content_mode(const char *value, EbConfig *cfg) {
    cfg->screen_content_mode = strtoul(value, NULL, 0);
};
// --- start: ALTREF_FILTERING_SUPPORT
static void set_enable_altrefs(const char *value, EbConfig *cfg) {
    cfg->enable_altrefs = (EbBool)strtoul(value, NULL, 0);
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
static void set_enable_hbd_mode_decision(const char *value, EbConfig *cfg) {
    cfg->enable_hbd_mode_decision = (uint8_t)strtoul(value, NULL, 0);
};
static void set_enable_palette(const char *value, EbConfig *cfg) {
    cfg->enable_palette = (int32_t)strtoul(value, NULL, 0);
};
static void set_enable_olpd_refinement(const char *value, EbConfig *cfg) {
    cfg->olpd_refinement = (int32_t)strtoul(value, NULL, 0);
};
static void set_high_dynamic_range_input(const char *value, EbConfig *cfg) {
    cfg->high_dynamic_range_input = strtol(value, NULL, 0);
};
static void set_profile(const char *value, EbConfig *cfg) {
    cfg->profile = strtol(value, NULL, 0);
};
static void set_tier(const char *value, EbConfig *cfg) { cfg->tier = strtol(value, NULL, 0); };
static void set_level(const char *value, EbConfig *cfg) {
    if (strtoul(value, NULL, 0) != 0 || EB_STRCMP(value, "0") == 0)
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
    if (cfg->injector_frame_rate > 1000)
        cfg->injector_frame_rate = cfg->injector_frame_rate;
    else
        cfg->injector_frame_rate = cfg->injector_frame_rate << 16;
}
static void set_latency_mode(const char *value, EbConfig *cfg) {
    cfg->latency_mode = (uint8_t)strtol(value, NULL, 0);
};
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
        if (EB_STRCMP(value, param_maps[i].name) == 0) {
            cfg->cpu_flags_limit = param_maps[i].flags;
            return;
        }
    }

    cfg->cpu_flags_limit = CPU_FLAGS_INVALID;
};
static void set_logical_processors(const char *value, EbConfig *cfg) {
    cfg->logical_processors = (uint32_t)strtoul(value, NULL, 0);
};
static void set_unpin_single_core_execution(const char *value, EbConfig *cfg) {
    cfg->unpin_lp1 = (uint32_t)strtoul(value, NULL, 0);
};
static void set_target_socket(const char *value, EbConfig *cfg) {
    cfg->target_socket = (int32_t)strtol(value, NULL, 0);
};
static void set_unrestricted_motion_vector(const char *value, EbConfig *cfg) {
    cfg->unrestricted_motion_vector = (EbBool)strtol(value, NULL, 0);
};

static void set_square_weight(const char *value, EbConfig *cfg) {
    cfg->sq_weight = (uint64_t)strtoul(value, NULL, 0);
    if (cfg->sq_weight == 0) cfg->sq_weight = (uint32_t)~0;
}

static void set_mds1_prune_c_th(const char *value, EbConfig *cfg) {
    cfg->md_stage_1_class_prune_th = (uint64_t)strtoul(value, NULL, 0);
    if (cfg->md_stage_1_class_prune_th == 0) cfg->md_stage_1_class_prune_th = (uint64_t)~0;
}
static void set_mds1_prune_s_th(const char *value, EbConfig *cfg) {
    cfg->md_stage_1_cand_prune_th = (uint64_t)strtoul(value, NULL, 0);
    if (cfg->md_stage_1_cand_prune_th == 0) cfg->md_stage_1_cand_prune_th = (uint64_t)~0;
}
static void set_mds2_prune_c_th(const char *value, EbConfig *cfg) {
    cfg->md_stage_2_class_prune_th = (uint64_t)strtoul(value, NULL, 0);
    if (cfg->md_stage_2_class_prune_th == 0) cfg->md_stage_2_class_prune_th = (uint64_t)~0;
}
static void set_mds2_prune_s_th(const char *value, EbConfig *cfg) {
    cfg->md_stage_2_cand_prune_th = (uint64_t)strtoul(value, NULL, 0);
    if (cfg->md_stage_2_cand_prune_th == 0) cfg->md_stage_2_cand_prune_th = (uint64_t)~0;
}

static void set_enable_auto_max_partition(const char *value, EbConfig *cfg) {
    cfg->enable_auto_max_partition = (uint8_t)strtol(value, NULL, 0);
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
ConfigEntry config_entry[] = {
    // File I/O
    {SINGLE_INPUT, INPUT_FILE_TOKEN, "InputFile", set_cfg_input_file},
    {SINGLE_INPUT, OUTPUT_BITSTREAM_TOKEN, "StreamFile", set_cfg_stream_file},
    {SINGLE_INPUT, ERROR_FILE_TOKEN, "ErrorFile", set_cfg_error_file},
    {SINGLE_INPUT, OUTPUT_RECON_TOKEN, "ReconFile", set_cfg_recon_file},
    {SINGLE_INPUT, QP_FILE_TOKEN, "QpFile", set_cfg_qp_file},
    {SINGLE_INPUT, STAT_FILE_TOKEN, "StatFile", set_cfg_stat_file},
    {SINGLE_INPUT, INPUT_STAT_FILE_TOKEN, "input_stat_file", set_input_stat_file},
    {SINGLE_INPUT, OUTPUT_STAT_FILE_TOKEN, "output_stat_file", set_output_stat_file},
    // Picture Dimensions
    {SINGLE_INPUT, WIDTH_TOKEN, "SourceWidth", set_cfg_source_width},
    {SINGLE_INPUT, HEIGHT_TOKEN, "SourceHeight", set_cfg_source_height},
    // Prediction Structure
    {SINGLE_INPUT, NUMBER_OF_PICTURES_TOKEN, "FrameToBeEncoded", set_cfg_frames_to_be_encoded},
    {SINGLE_INPUT, BUFFERED_INPUT_TOKEN, "BufferedInput", set_buffered_input},
    {SINGLE_INPUT, BASE_LAYER_SWITCH_MODE_TOKEN, "BaseLayerSwitchMode", set_base_layer_switch_mode},
    {SINGLE_INPUT, ENCMODE_TOKEN, "EncoderMode", set_enc_mode},
    {SINGLE_INPUT, ENCMODE2P_TOKEN, "EncoderMode2p", set_snd_pass_enc_mode},
    {SINGLE_INPUT, INTRA_PERIOD_TOKEN, "IntraPeriod", set_cfg_intra_period},
    {SINGLE_INPUT, INTRA_REFRESH_TYPE_TOKEN, "IntraRefreshType", set_cfg_intra_refresh_type},
    {SINGLE_INPUT, FRAME_RATE_TOKEN, "FrameRate", set_frame_rate},
    {SINGLE_INPUT, FRAME_RATE_NUMERATOR_TOKEN, "FrameRateNumerator", set_frame_rate_numerator},
    {SINGLE_INPUT, FRAME_RATE_DENOMINATOR_TOKEN, "FrameRateDenominator", set_frame_rate_denominator},
    {SINGLE_INPUT, ENCODER_BIT_DEPTH, "EncoderBitDepth", set_encoder_bit_depth},
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
    {SINGLE_INPUT, TARGET_BIT_RATE_TOKEN, "TargetBitRate", set_target_bit_rate},
    {SINGLE_INPUT, MAX_QP_TOKEN, "MaxQpAllowed", set_max_qp_allowed},
    {SINGLE_INPUT, MIN_QP_TOKEN, "MinQpAllowed", set_min_qp_allowed},
    {SINGLE_INPUT, VBV_BUFSIZE_TOKEN, "VBVBufSize", set_vbv_buf_size},
    {SINGLE_INPUT, ADAPTIVE_QP_ENABLE_TOKEN, "AdaptiveQuantization", set_adaptive_quantization},

    // DLF
    {SINGLE_INPUT, LOOP_FILTER_DISABLE_TOKEN, "LoopFilterDisable", set_disable_dlf_flag},

    // RESTORATION
    {SINGLE_INPUT,
     RESTORATION_ENABLE_TOKEN,
     "RestorationFilter",
     set_enable_restoration_filter_flag},

    {SINGLE_INPUT, MFMV_ENABLE_TOKEN, "Mfmv", set_enable_mfmv_flag},
    {SINGLE_INPUT, REDUNDANT_BLK_TOKEN, "RedundantBlock", set_enable_redundant_blk_flag},
    {SINGLE_INPUT, TRELLIS_ENABLE_TOKEN, "Trellis", set_enable_trellis_flag},
    {SINGLE_INPUT, SPATIAL_SSE_FL_TOKEN, "SpatialSSEfl", set_spatial_sse_fl_flag},
    {SINGLE_INPUT, SUBPEL_TOKEN, "Subpel", set_enable_sub_pel_flag},
    {SINGLE_INPUT, OVR_BNDRY_BLK_TOKEN, "OverBoundryBlock", set_over_bndry_blk_flag},
    {SINGLE_INPUT,
     NEW_NEAREST_COMB_INJECT_TOKEN,
     "NewNearestCombInjection",
     set_new_nearest_comb_inject_flag},
    {SINGLE_INPUT,
     NX4_4XN_MV_INJECT_TOKEN,
     "nx4ParentMvInjection",
     set_nx4_4xn_parent_mv_inject_flag},
    {SINGLE_INPUT, PRUNE_UNIPRED_ME_TOKEN, "PruneUnipredMe", set_prune_unipred_me_flag},
    {SINGLE_INPUT, PRUNE_REF_REC_PART_TOKEN, "PruneRefRecPart", set_prune_ref_rec_part_flag},
    {SINGLE_INPUT, NSQ_TABLE_TOKEN, "NsqTable", set_nsq_table_flag},
    {SINGLE_INPUT, FRAME_END_CDF_UPDATE_TOKEN, "FrameEndCdfUpdate", set_frame_end_cdf_update_flag},

    // CHROMA
    {SINGLE_INPUT, CHROMA_MODE_TOKEN, "ChromaMode", set_chroma_mode},

    // LOCAL WARPED MOTION
    {SINGLE_INPUT,
     LOCAL_WARPED_ENABLE_TOKEN,
     "LocalWarpedMotion",
     set_enable_local_warped_motion_flag},
    // GLOBAL MOTION
    {SINGLE_INPUT, GLOBAL_MOTION_ENABLE_TOKEN, "GlobalMotion", set_enable_global_motion_flag},

    // CLASS 12
    {SINGLE_INPUT, CLASS_12_TOKEN, "CombineClass12", set_class_12_flag},
    // EDGE SKIP ANGLE INTRA
    {SINGLE_INPUT,
     EDGE_SKIP_ANGLE_INTRA_TOKEN,
     "EdgeSkipAngleIntra",
     set_edge_skip_angle_intra_flag},
    // INTER INTRA COMPOUND
    {SINGLE_INPUT, INTER_INTRA_COMPOUND_TOKEN, "InterIntraCompound", set_interintra_compound_flag},
    // FRACTIONAL SEARCH 64x64
    {SINGLE_INPUT, FRAC_SEARCH_64_TOKEN, "FractionalSearch64", set_fractional_search_64_flag},

    // OBMC
    {SINGLE_INPUT, OBMC_TOKEN, "Obmc", set_enable_obmc_flag},
    // RDOQ
    {SINGLE_INPUT, RDOQ_TOKEN, "RDOQ", set_enable_rdoq_flag},

    // Filter Intra
    {SINGLE_INPUT, FILTER_INTRA_TOKEN, "FilterIntra", set_enable_filter_intra_flag},

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
    // HME Parameters
    {SINGLE_INPUT,
     NUM_HME_SEARCH_WIDTH_TOKEN,
     "NumberHmeSearchRegionInWidth",
     set_cfg_number_hme_search_region_in_width},
    {SINGLE_INPUT,
     NUM_HME_SEARCH_HEIGHT_TOKEN,
     "NumberHmeSearchRegionInHeight",
     set_cfg_number_hme_search_region_in_height},
    {SINGLE_INPUT,
     HME_SRCH_T_L0_WIDTH_TOKEN,
     "HmeLevel0TotalSearchAreaWidth",
     set_cfg_hme_level_0_total_search_area_width},
    {SINGLE_INPUT,
     HME_SRCH_T_L0_HEIGHT_TOKEN,
     "HmeLevel0TotalSearchAreaHeight",
     set_cfg_hme_level_0_total_search_area_height},
    // MD Parameters
    {SINGLE_INPUT, SCREEN_CONTENT_TOKEN, "ScreenContentMode", set_screen_content_mode},
    {SINGLE_INPUT, HBD_MD_ENABLE_TOKEN, "HighBitDepthModeDecision", set_enable_hbd_mode_decision},
    {SINGLE_INPUT, PALETTE_TOKEN, "PaletteMode", set_enable_palette},
    {SINGLE_INPUT, OLPD_REFINEMENT_TOKEN, "OlpdRefinement", set_enable_olpd_refinement},
    // Thread Management
    {SINGLE_INPUT, THREAD_MGMNT, "LogicalProcessors", set_logical_processors},
    {SINGLE_INPUT, UNPIN_LP1_TOKEN, "UnpinSingleCoreExecution", set_unpin_single_core_execution},
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
    {SINGLE_INPUT, SPEED_CONTROL_TOKEN, "speed_control_flag", speed_control_flag},
    // Annex A parameters
    {SINGLE_INPUT, PROFILE_TOKEN, "Profile", set_profile},
    {SINGLE_INPUT, TIER_TOKEN, "Tier", set_tier},
    {SINGLE_INPUT, LEVEL_TOKEN, "Level", set_level},
    {SINGLE_INPUT, LATENCY_MODE, "LatencyMode", set_latency_mode},
    {SINGLE_INPUT, FILM_GRAIN_TOKEN, "FilmGrain", set_cfg_film_grain},
    // Asm Type
    {SINGLE_INPUT, ASM_TYPE_TOKEN, "Asm", set_asm_type},
    // HME
    {ARRAY_INPUT,
     HME_LEVEL0_WIDTH,
     "HmeLevel0SearchAreaInWidth",
     set_hme_level_0_search_area_in_width_array},
    {ARRAY_INPUT,
     HME_LEVEL0_HEIGHT,
     "HmeLevel0SearchAreaInHeight",
     set_hme_level_0_search_area_in_height_array},
    {ARRAY_INPUT,
     HME_LEVEL1_WIDTH,
     "HmeLevel1SearchAreaInWidth",
     set_hme_level_1_search_area_in_width_array},
    {ARRAY_INPUT,
     HME_LEVEL1_HEIGHT,
     "HmeLevel1SearchAreaInHeight",
     set_hme_level_1_search_area_in_height_array},
    {ARRAY_INPUT,
     HME_LEVEL2_WIDTH,
     "HmeLevel2SearchAreaInWidth",
     set_hme_level_2_search_area_in_width_array},
    {ARRAY_INPUT,
     HME_LEVEL2_HEIGHT,
     "HmeLevel2SearchAreaInHeight",
     set_hme_level_2_search_area_in_height_array},
    // --- start: ALTREF_FILTERING_SUPPORT
    {SINGLE_INPUT, ENABLE_ALTREFS, "EnableAltRefs", set_enable_altrefs},
    {SINGLE_INPUT, ALTREF_STRENGTH, "AltRefStrength", set_altref_strength},
    {SINGLE_INPUT, ALTREF_NFRAMES, "AltRefNframes", set_altref_n_frames},
    {SINGLE_INPUT, ENABLE_OVERLAYS, "EnableOverlays", set_enable_overlays},
    // --- end: ALTREF_FILTERING_SUPPORT

    {SINGLE_INPUT, SQ_WEIGHT_TOKEN, "SquareWeight", set_square_weight},
    {SINGLE_INPUT, ENABLE_AMP_TOKEN, "AutomaxPartition", set_enable_auto_max_partition},

    {SINGLE_INPUT, MDS1_PRUNE_C_TH, "MDStage1PruneClassThreshold", set_mds1_prune_c_th},
    {SINGLE_INPUT, MDS1_PRUNE_S_TH, "MDStage1PruneCandThreshold", set_mds1_prune_s_th},
    {SINGLE_INPUT, MDS2_PRUNE_C_TH, "MDStage2PruneClassThreshold", set_mds2_prune_c_th},
    {SINGLE_INPUT, MDS2_PRUNE_S_TH, "MDStage2PruneCandThreshold", set_mds2_prune_s_th},

    // Termination
    {SINGLE_INPUT, NULL, NULL, NULL}};

/**********************************
 * Constructor
 **********************************/
void eb_config_ctor(EbConfig *config_ptr) {
    memset(config_ptr, 0, sizeof(*config_ptr));
    config_ptr->error_log_file       = stderr;
    config_ptr->frame_rate           = 30 << 16;
    config_ptr->encoder_bit_depth    = 8;
    config_ptr->encoder_color_format = 1; //EB_YUV420
    config_ptr->buffered_input       = -1;

    config_ptr->qp                  = 50;
    config_ptr->use_qp_file         = EB_FALSE;
    config_ptr->look_ahead_distance = (uint32_t)~0;
    config_ptr->target_bit_rate     = 7000000;
    config_ptr->max_qp_allowed      = 63;
    config_ptr->min_qp_allowed      = 10;

    config_ptr->enable_adaptive_quantization              = 2;
    config_ptr->enc_mode                                  = MAX_ENC_PRESET;
    config_ptr->snd_pass_enc_mode                         = MAX_ENC_PRESET + 1;
    config_ptr->intra_period                              = -2;
    config_ptr->intra_refresh_type                        = 1;
    config_ptr->hierarchical_levels                       = 4;
    config_ptr->pred_structure                            = 2;
    config_ptr->enable_global_motion                      = EB_TRUE;
    config_ptr->enable_restoration_filtering              = DEFAULT;
    config_ptr->combine_class_12                          = DEFAULT;
    config_ptr->edge_skp_angle_intra                      = DEFAULT;
    config_ptr->inter_intra_compound                      = DEFAULT;
    config_ptr->fract_search_64                           = DEFAULT;
    config_ptr->enable_mfmv                               = DEFAULT;
    config_ptr->enable_redundant_blk                      = DEFAULT;
    config_ptr->enable_trellis                            = DEFAULT;
    config_ptr->spatial_sse_fl                            = DEFAULT;
    config_ptr->enable_subpel                             = DEFAULT;
    config_ptr->over_bndry_blk                            = DEFAULT;
    config_ptr->new_nearest_comb_inject                   = DEFAULT;
    config_ptr->nx4_4xn_parent_mv_inject                  = DEFAULT;
    config_ptr->prune_unipred_me                          = DEFAULT;
    config_ptr->prune_ref_rec_part                        = DEFAULT;
    config_ptr->nsq_table                                 = DEFAULT;
    config_ptr->frame_end_cdf_update                      = DEFAULT;
    config_ptr->set_chroma_mode                           = DEFAULT;
    config_ptr->enable_obmc                               = EB_TRUE;
    config_ptr->enable_rdoq                               = DEFAULT;
    config_ptr->pred_me                                   = DEFAULT;
    config_ptr->bipred_3x3_inject                         = DEFAULT;
    config_ptr->compound_level                            = DEFAULT;
    config_ptr->enable_filter_intra                       = EB_TRUE;
    config_ptr->use_default_me_hme                        = EB_TRUE;
    config_ptr->enable_hme_flag                           = EB_TRUE;
    config_ptr->enable_hme_level0_flag                    = EB_TRUE;
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
    config_ptr->enable_hbd_mode_decision                  = 2;
    config_ptr->enable_palette                            = -1;
    config_ptr->olpd_refinement                           = -1;
    config_ptr->injector_frame_rate                       = 60 << 16;

    // ASM Type
    config_ptr->cpu_flags_limit = CPU_FLAGS_ALL;

    config_ptr->unpin_lp1     = 1;
    config_ptr->target_socket = -1;

    config_ptr->unrestricted_motion_vector = EB_TRUE;

    // --- start: ALTREF_FILTERING_SUPPORT
    config_ptr->enable_altrefs  = EB_TRUE;
    config_ptr->altref_strength = 5;
    config_ptr->altref_nframes  = 7;
    // --- end: ALTREF_FILTERING_SUPPORT

    config_ptr->sq_weight                 = 100;
    config_ptr->enable_auto_max_partition = 1;

    config_ptr->md_stage_1_cand_prune_th  = 75;
    config_ptr->md_stage_1_class_prune_th = 100;
    config_ptr->md_stage_2_cand_prune_th  = 15;
    config_ptr->md_stage_2_class_prune_th = 25;

    return;
}

/**********************************
 * Destructor
 **********************************/
void eb_config_dtor(EbConfig *config_ptr) {
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
    if (config_ptr->input_stat_file) {
        fclose(config_ptr->input_stat_file);
        config_ptr->input_stat_file = (FILE *)NULL;
    }
    if (config_ptr->output_stat_file) {
        fclose(config_ptr->output_stat_file);
        config_ptr->output_stat_file = (FILE *)NULL;
    }
    return;
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
        if (EB_STRCMP(config_entry[i].name, name) == 0)
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
                EB_STRNCPY(var_name, CONFIG_FILE_MAX_VAR_LEN, argv[0], arg_len[0]);
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
                    EB_STRNCPY(var_value[value_index],
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
        return_error = EB_STRCMP(argv[--argc], token);
        if (return_error == 0) { EB_STRCPY(configStr, COMMAND_LINE_MAX_SIZE, argv[argc + 1]); }
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
    if (config->injector_frame_rate <= 0 && config->injector) {
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
        return_error = EB_STRCMP(argv[--argc], token);
        if (return_error == 0) {
            int32_t count;
            for (count = 0; count < MAX_CHANNEL_NUMBER; ++count) {
                if (done == 0) {
                    if (argv[argc + count + 1]) {
                        if (strtoul(argv[argc + count + 1], NULL, 0) != 0 ||
                            EB_STRCMP(argv[argc + count + 1], "0") == 0) {
                            EB_STRCPY(
                                configStr[count], COMMAND_LINE_MAX_SIZE, argv[argc + count + 1]);
                        } else if (argv[argc + count + 1][0] != '-') {
                            EB_STRCPY(
                                configStr[count], COMMAND_LINE_MAX_SIZE, argv[argc + count + 1]);
                        } else {
                            EB_STRCPY(configStr[count], COMMAND_LINE_MAX_SIZE, " ");
                            done = 1;
                        }
                    } else {
                        EB_STRCPY(configStr[count], COMMAND_LINE_MAX_SIZE, " ");
                        done = 1;
                        //return return_error;
                    }
                } else
                    EB_STRCPY(configStr[count], COMMAND_LINE_MAX_SIZE, " ");
            }
        }
    }

    return return_error;
}

uint32_t get_help(int32_t argc, char *const argv[]) {
    char config_string[COMMAND_LINE_MAX_SIZE];
    if (find_token(argc, argv, HELP_TOKEN, config_string) == 0) {
        int32_t token_index = -1;

        fprintf(stderr, "\n%-25s\t%-25s\t%s\n\n", "TOKEN", "DESCRIPTION", "INPUT TYPE");
        fprintf(stderr, "%-25s\t%-25s\t%s\n", "-nch", "NumberOfChannels", "Single input");
        while (config_entry[++token_index].token != NULL)
            fprintf(stderr,
                    "%-25s\t%-25s\t%s\n",
                    config_entry[token_index].token,
                    config_entry[token_index].name,
                    config_entry[token_index].type ? "Array input" : "Single input");
        return 1;
    } else
        return 0;
}

/******************************************************
* Get the number of channels and validate it with input
******************************************************/
uint32_t get_number_of_channels(int32_t argc, char *const argv[]) {
    char     config_string[COMMAND_LINE_MAX_SIZE];
    uint32_t channel_number;
    if (find_token(argc, argv, CHANNEL_NUMBER_TOKEN, config_string) == 0) {
        // Set the input file
        channel_number = strtol(config_string, NULL, 0);
        if ((channel_number > (uint32_t)MAX_CHANNEL_NUMBER) || channel_number == 0) {
            fprintf(stderr,
                    "Error: The number of channels has to be within the range [1,%u]\n",
                    (uint32_t)MAX_CHANNEL_NUMBER);
            return 0;
        } else {
            return channel_number;
        }
    }
    return 1;
}

void mark_token_as_read(const char *token, char *cmd_copy[], int32_t *cmd_token_cnt) {
    int32_t cmd_copy_index;
    for (cmd_copy_index = 0; cmd_copy_index < *(cmd_token_cnt); ++cmd_copy_index) {
        if (!EB_STRCMP(cmd_copy[cmd_copy_index], token))
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
    uint64_t curr_loc;

    if (config->input_file) {
        curr_loc = ftello(config->input_file); // get current fp location
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

/******************************************
* Read Command Line
******************************************/
EbErrorType read_command_line(int32_t argc, char *const argv[], EbConfig **configs,
                              uint32_t num_channels, EbErrorType *return_errors) {
    EbErrorType return_error = EB_ErrorBadParameter;
    char        config_string[COMMAND_LINE_MAX_SIZE]; // for one input options
    char *      config_strings[MAX_CHANNEL_NUMBER]; // for multiple input options
    char *      cmd_copy[MAX_NUM_TOKENS]; // keep track of extra tokens
    uint32_t    index         = 0;
    int32_t     cmd_token_cnt = 0; // total number of tokens
    int32_t     token_index   = -1;
    int32_t     ret_y4m;

    for (index = 0; index < MAX_CHANNEL_NUMBER; ++index)
        config_strings[index] = (char *)malloc(sizeof(char) * COMMAND_LINE_MAX_SIZE);
    // Copy tokens (except for CHANNEL_NUMBER_TOKEN ) into a temp token buffer hosting all tokens that are passed through the command line
    size_t len = EB_STRLEN(CHANNEL_NUMBER_TOKEN, COMMAND_LINE_MAX_SIZE);
    for (token_index = 0; token_index < argc; ++token_index) {
        if ((argv[token_index][0] == '-') &&
            strncmp(argv[token_index], CHANNEL_NUMBER_TOKEN, len) &&
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
            return_errors[index] =
                (EbErrorType)read_config_file(configs[index], config_strings[index], index);
            return_error = (EbErrorType)(return_error & return_errors[index]);
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
    token_index = -1;
    // Parse command line for tokens
    while (config_entry[++token_index].name != NULL) {
        if (config_entry[token_index].type == SINGLE_INPUT) {
            if (find_token_multiple_inputs(
                    argc, argv, config_entry[token_index].token, config_strings) == 0) {
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
        if ((configs[index])->y4m_input == EB_TRUE) {
            ret_y4m = read_y4m_header(configs[index]);
            if (ret_y4m == EB_ErrorBadParameter) {
                fprintf(stderr, "Error found when reading the y4m file parameters.\n");
                return EB_ErrorBadParameter;
            }
        }
    }

    /***************************************************************************************************/
    /***********   Find SPECIAL configuration parameter tokens and call respective functions  **********/
    /***************************************************************************************************/
    // Parse command line for search region at level 0 width token
    if (find_token_multiple_inputs(argc, argv, HME_LEVEL0_WIDTH, config_strings) == 0) {
        uint32_t input_index = 0, last_index = 0;
        uint32_t done = 1;

        mark_token_as_read(HME_LEVEL0_WIDTH, cmd_copy, &cmd_token_cnt);

        for (index = 0; done && (index < num_channels); ++index) {
            configs[index]->hme_level0_column_index = 0;
            for (input_index = last_index;
                 input_index < configs[index]->number_hme_search_region_in_width + last_index;
                 ++input_index) {
                if (EB_STRCMP(config_strings[input_index], " "))
                    set_hme_level_0_search_area_in_width_array(config_strings[input_index],
                                                               configs[index]);
                else {
                    done = 0;
                    break;
                }
            }
            last_index += configs[index]->number_hme_search_region_in_width;
        }
    }

    //// Parse command line for search region at level 0 height token
    if (find_token_multiple_inputs(argc, argv, HME_LEVEL0_HEIGHT, config_strings) == 0) {
        uint32_t input_index = 0, last_index = 0;
        uint32_t done = 1;

        mark_token_as_read(HME_LEVEL0_HEIGHT, cmd_copy, &cmd_token_cnt);

        for (index = 0; done && (index < num_channels); ++index) {
            configs[index]->hme_level0_row_index = 0;
            for (input_index = last_index;
                 input_index < configs[index]->number_hme_search_region_in_height + last_index;
                 ++input_index) {
                if (EB_STRCMP(config_strings[input_index], " "))
                    set_hme_level_0_search_area_in_height_array(config_strings[input_index],
                                                                configs[index]);
                else {
                    done = 0;
                    break;
                }
            }
            last_index += configs[index]->number_hme_search_region_in_height;
        }
    }

    // Parse command line for search region at level 1 Height token
    if (find_token_multiple_inputs(argc, argv, HME_LEVEL1_HEIGHT, config_strings) == 0) {
        uint32_t input_index = 0, last_index = 0;
        uint32_t done = 1;

        mark_token_as_read(HME_LEVEL1_HEIGHT, cmd_copy, &cmd_token_cnt);

        for (index = 0; done && (index < num_channels); ++index) {
            configs[index]->hme_level1_row_index = 0;
            for (input_index = last_index;
                 input_index < configs[index]->number_hme_search_region_in_height + last_index;
                 ++input_index) {
                if (EB_STRCMP(config_strings[input_index], " "))
                    set_hme_level_1_search_area_in_height_array(config_strings[input_index],
                                                                configs[index]);
                else {
                    done = 0;
                    break;
                }
            }
            last_index += configs[index]->number_hme_search_region_in_height;
        }
    }

    // Parse command line for search region at level 1 width token
    if (find_token_multiple_inputs(argc, argv, HME_LEVEL1_WIDTH, config_strings) == 0) {
        uint32_t input_index = 0, last_index = 0;
        uint32_t done = 1;

        mark_token_as_read(HME_LEVEL1_WIDTH, cmd_copy, &cmd_token_cnt);

        for (index = 0; done && (index < num_channels); ++index) {
            configs[index]->hme_level1_column_index = 0;
            for (input_index = last_index;
                 input_index < configs[index]->number_hme_search_region_in_width + last_index;
                 ++input_index) {
                if (EB_STRCMP(config_strings[input_index], " "))
                    set_hme_level_1_search_area_in_width_array(config_strings[input_index],
                                                               configs[index]);
                else {
                    done = 0;
                    break;
                }
            }
            last_index += configs[index]->number_hme_search_region_in_width;
        }
    }

    // Parse command line for search region at level 2 width token
    if (find_token_multiple_inputs(argc, argv, HME_LEVEL2_WIDTH, config_strings) == 0) {
        uint32_t input_index = 0, last_index = 0;
        uint32_t done = 1;

        mark_token_as_read(HME_LEVEL2_WIDTH, cmd_copy, &cmd_token_cnt);

        for (index = 0; done && (index < num_channels); ++index) {
            configs[index]->hme_level2_column_index = 0;
            for (input_index = last_index;
                 input_index < configs[index]->number_hme_search_region_in_width + last_index;
                 ++input_index) {
                if (EB_STRCMP(config_strings[input_index], " "))
                    set_hme_level_2_search_area_in_width_array(config_strings[input_index],
                                                               configs[index]);
                else {
                    done = 0;
                    break;
                }
            }
            last_index += configs[index]->number_hme_search_region_in_width;
        }
    }

    // Parse command line for search region at level 2 height token
    if (find_token_multiple_inputs(argc, argv, HME_LEVEL2_HEIGHT, config_strings) == 0) {
        uint32_t input_index = 0, last_index = 0;
        uint32_t done = 1;

        mark_token_as_read(HME_LEVEL2_HEIGHT, cmd_copy, &cmd_token_cnt);

        for (index = 0; done && (index < num_channels); ++index) {
            configs[index]->hme_level2_row_index = 0;
            for (input_index = last_index;
                 input_index < configs[index]->number_hme_search_region_in_height + last_index;
                 ++input_index) {
                if (EB_STRCMP(config_strings[input_index], " "))
                    set_hme_level_2_search_area_in_height_array(config_strings[input_index],
                                                                configs[index]);
                else {
                    done = 0;
                    break;
                }
            }
            last_index += configs[index]->number_hme_search_region_in_height;
        }
    }

    /***************************************************************************************************/
    /**************************************   Verify configuration parameters   ************************/
    /***************************************************************************************************/
    // Verify the config values
    if (return_error == 0) {
        return_error = EB_ErrorBadParameter;
        for (index = 0; index < num_channels; ++index) {
            if (return_errors[index] == EB_ErrorNone) {
                return_errors[index] = verify_settings(configs[index], index);

                // Assuming no errors, add padding to width and height
                if (return_errors[index] == EB_ErrorNone) {
                    configs[index]->input_padded_width =
                        configs[index]->source_width + configs[index]->source_width % 8;
                    configs[index]->input_padded_height =
                        configs[index]->source_height + configs[index]->source_width % 8;
                }

                // Assuming no errors, set the frames to be encoded to the number of frames in the input yuv
                if (return_errors[index] == EB_ErrorNone &&
                    configs[index]->frames_to_be_encoded == 0)
                    configs[index]->frames_to_be_encoded =
                        compute_frames_to_be_encoded(configs[index]);

                if (configs[index]->frames_to_be_encoded == -1) {
                    fprintf(configs[index]->error_log_file,
                            "Error instance %u: Input yuv does not contain enough frames \n",
                            index + 1);
                    return_errors[index] = EB_ErrorBadParameter;
                }

                // Force the injector latency mode, and injector frame rate when speed control is on
                if (return_errors[index] == EB_ErrorNone && configs[index]->speed_control_flag == 1)
                    configs[index]->injector = 1;
            }
            return_error = (EbErrorType)(return_error & return_errors[index]);
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
