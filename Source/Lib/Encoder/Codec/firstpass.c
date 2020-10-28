/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <limits.h>
#include <math.h>
#include <stdio.h>

#include "aom_dsp_rtcd.h"
#include "EbDefinitions.h"
#include "EbRateControlProcess.h"
#include "EbSequenceControlSet.h"
#include "EbPictureControlSet.h"
#include "firstpass.h"
#include "EbLog.h"
#include "EbModeDecisionProcess.h"
#include "EbCodingLoop.h"
#include "dwt.h" // to move to firstpass.c
#include "EbPictureDecisionProcess.h"
#include "EbModeDecisionConfigurationProcess.h"
#include "mv.h"

#define INCRMENT_CAND_TOTAL_COUNT(cnt)                                                     \
    MULTI_LINE_MACRO_BEGIN cnt++;                                                          \
    if (cnt >= MODE_DECISION_CANDIDATE_MAX_COUNT_Y)                                          \
        SVT_LOG(" ERROR: reaching limit for MODE_DECISION_CANDIDATE_MAX_COUNT %i\n", cnt); \
    MULTI_LINE_MACRO_END

#define OUTPUT_FPF 0

#define INTRA_MODE_PENALTY 1024
#define NEW_MV_MODE_PENALTY 32
#define DARK_THRESH 64

#define NCOUNT_INTRA_THRESH 8192
#define NCOUNT_INTRA_FACTOR 3

#define STATS_CAPABILITY_INIT 100
//1.5 times larger than request.
#define STATS_CAPABILITY_GROW(s) (s * 3 /2)

static EbErrorType realloc_stats_out(FirstPassStatsOut* out, uint64_t frame_number) {
    if (frame_number < out->size)
        return EB_ErrorNone;
    if (frame_number >= out->capability) {
        size_t capability = frame_number >= STATS_CAPABILITY_INIT ?
            STATS_CAPABILITY_GROW(frame_number) : STATS_CAPABILITY_INIT;
        EB_REALLOC_ARRAY(out->stat, capability);
        out->capability = capability;
    }
    out->size = frame_number + 1;
    return EB_ErrorNone;
}

static AOM_INLINE void output_stats(SequenceControlSet *scs_ptr, FIRSTPASS_STATS *stats,
                                    uint64_t frame_number) {
    FirstPassStatsOut* stats_out = &scs_ptr->encode_context_ptr->stats_out;
    svt_block_on_mutex(scs_ptr->encode_context_ptr->stat_file_mutex);
    if (realloc_stats_out(stats_out, frame_number) != EB_ErrorNone) {
        SVT_ERROR("realloc_stats_out request %d entries failed failed\n", frame_number);
    } else {
        stats_out->stat[frame_number] = *stats;
    }
// TEMP debug code
#if OUTPUT_FPF
    {
        FILE *fpfile;
        if (frame_number == 0)
            fpfile = fopen("firstpass.stt", "w");
        else
            fpfile = fopen("firstpass.stt", "a");
        fprintf(fpfile,
                "%12.0lf %12.4lf %12.0lf %12.0lf %12.0lf %12.4lf %12.4lf"
                "%12.4lf %12.4lf %12.4lf %12.4lf %12.4lf %12.4lf %12.4lf %12.4lf"
                "%12.4lf %12.4lf %12.0lf %12.4lf %12.4lf %12.4lf %12.4lf\n",
                stats->frame,
                stats->weight,
                stats->intra_error,
                stats->coded_error,
                stats->sr_coded_error,
                stats->pcnt_inter,
                stats->pcnt_motion,
                stats->pcnt_second_ref,
                stats->pcnt_neutral,
                stats->intra_skip_pct,
                stats->inactive_zone_rows,
                stats->inactive_zone_cols,
                stats->MVr,
                stats->mvr_abs,
                stats->MVc,
                stats->mvc_abs,
                stats->MVrv,
                stats->MVcv,
                stats->mv_in_out_count,
                stats->new_mv_count,
                stats->count,
                stats->duration);
        fclose(fpfile);
    }
#endif
    svt_release_mutex(scs_ptr->encode_context_ptr->stat_file_mutex);
}
void svt_av1_twopass_zero_stats(FIRSTPASS_STATS *section) {
    section->frame                    = 0.0;
    section->weight                   = 0.0;
    section->intra_error              = 0.0;
    section->frame_avg_wavelet_energy = 0.0;
    section->coded_error              = 0.0;
    section->sr_coded_error           = 0.0;
    section->pcnt_inter               = 0.0;
    section->pcnt_motion              = 0.0;
    section->pcnt_second_ref          = 0.0;
    section->pcnt_neutral             = 0.0;
    section->intra_skip_pct           = 0.0;
    section->inactive_zone_rows       = 0.0;
    section->inactive_zone_cols       = 0.0;
    section->MVr                      = 0.0;
    section->mvr_abs                  = 0.0;
    section->MVc                      = 0.0;
    section->mvc_abs                  = 0.0;
    section->MVrv                     = 0.0;
    section->MVcv                     = 0.0;
    section->mv_in_out_count          = 0.0;
    section->new_mv_count             = 0.0;
    section->count                    = 0.0;
    section->duration                 = 1.0;
}
void svt_av1_accumulate_stats(FIRSTPASS_STATS *section, const FIRSTPASS_STATS *frame) {
    section->frame += frame->frame;
    section->weight += frame->weight;
    section->intra_error += frame->intra_error;
    section->frame_avg_wavelet_energy += frame->frame_avg_wavelet_energy;
    section->coded_error += frame->coded_error;
    section->sr_coded_error += frame->sr_coded_error;
    section->pcnt_inter += frame->pcnt_inter;
    section->pcnt_motion += frame->pcnt_motion;
    section->pcnt_second_ref += frame->pcnt_second_ref;
    section->pcnt_neutral += frame->pcnt_neutral;
    section->intra_skip_pct += frame->intra_skip_pct;
    section->inactive_zone_rows += frame->inactive_zone_rows;
    section->inactive_zone_cols += frame->inactive_zone_cols;
    section->MVr += frame->MVr;
    section->mvr_abs += frame->mvr_abs;
    section->MVc += frame->MVc;
    section->mvc_abs += frame->mvc_abs;
    section->MVrv += frame->MVrv;
    section->MVcv += frame->MVcv;
    section->mv_in_out_count += frame->mv_in_out_count;
    section->new_mv_count += frame->new_mv_count;
    section->count += frame->count;
    section->duration += frame->duration;
}
void svt_av1_end_first_pass(PictureParentControlSet *pcs_ptr) {
    SequenceControlSet *scs_ptr = pcs_ptr->scs_ptr;
    TWO_PASS *          twopass = &scs_ptr->twopass;

    if (twopass->stats_buf_ctx->total_stats) {
        // add the total to the end of the file
        output_stats(scs_ptr, twopass->stats_buf_ctx->total_stats, pcs_ptr->picture_number + 1);
    }
}
static double raw_motion_error_stdev(int *raw_motion_err_list, int raw_motion_err_counts) {
    int64_t sum_raw_err   = 0;
    double  raw_err_avg   = 0;
    double  raw_err_stdev = 0;
    if (raw_motion_err_counts == 0) return 0;

    int i;
    for (i = 0; i < raw_motion_err_counts; i++) { sum_raw_err += raw_motion_err_list[i]; }
    raw_err_avg = (double)sum_raw_err / raw_motion_err_counts;
    for (i = 0; i < raw_motion_err_counts; i++) {
        raw_err_stdev +=
            (raw_motion_err_list[i] - raw_err_avg) * (raw_motion_err_list[i] - raw_err_avg);
    }
    // Calculate the standard deviation for the motion error of all the inter
    // blocks of the 0,0 motion using the last source
    // frame as the reference.
    raw_err_stdev = sqrt(raw_err_stdev / raw_motion_err_counts);
    return raw_err_stdev;
}
#define UL_INTRA_THRESH 50
#define INVALID_ROW -1
// Accumulates motion vector stats.
// Modifies member variables of "stats".
void accumulate_mv_stats(const MV best_mv, const FULLPEL_MV mv, const int mb_row, const int mb_col,
                         const int mb_rows, const int mb_cols, MV *last_mv, FRAME_STATS *stats) {
    if (is_zero_mv(&best_mv)) return;

    ++stats->mv_count;
    // Non-zero vector, was it different from the last non zero vector?
    if (!is_equal_mv(&best_mv, last_mv)) ++stats->new_mv_count;
    *last_mv = best_mv;

    // Does the row vector point inwards or outwards?
    if (mb_row < mb_rows / 2) {
        if (mv.row > 0) {
            --stats->sum_in_vectors;
        } else if (mv.row < 0) {
            ++stats->sum_in_vectors;
        }
    } else if (mb_row > mb_rows / 2) {
        if (mv.row > 0) {
            ++stats->sum_in_vectors;
        } else if (mv.row < 0) {
            --stats->sum_in_vectors;
        }
    }

    // Does the col vector point inwards or outwards?
    if (mb_col < mb_cols / 2) {
        if (mv.col > 0) {
            --stats->sum_in_vectors;
        } else if (mv.col < 0) {
            ++stats->sum_in_vectors;
        }
    } else if (mb_col > mb_cols / 2) {
        if (mv.col > 0) {
            ++stats->sum_in_vectors;
        } else if (mv.col < 0) {
            --stats->sum_in_vectors;
        }
    }
}
// Updates the first pass stats of this frame.
// Input:
//   cpi: the encoder setting. Only a few params in it will be used.
//   stats: stats accumulated for this frame.
//   raw_err_stdev: the statndard deviation for the motion error of all the
//                  inter blocks of the (0,0) motion using the last source
//                  frame as the reference.
//   frame_number: current frame number.
//   ts_duration: Duration of the frame / collection of frames.
// Updates:
//   twopass->total_stats: the accumulated stats.
//   twopass->stats_buf_ctx->stats_in_end: the pointer to the current stats,
//                                         update its value and its position
//                                         in the buffer.
static void update_firstpass_stats(PictureParentControlSet *pcs_ptr,
                                   const FRAME_STATS *const stats, const double raw_err_stdev,
                                   const int frame_number, const int64_t ts_duration) {
    SequenceControlSet *scs_ptr = pcs_ptr->scs_ptr;
    TWO_PASS *          twopass = &scs_ptr->twopass;

    const uint32_t   mb_cols          = (scs_ptr->seq_header.max_frame_width + 16 - 1) / 16;
    const uint32_t   mb_rows          = (scs_ptr->seq_header.max_frame_height + 16 - 1) / 16;
    FIRSTPASS_STATS *this_frame_stats = twopass->stats_buf_ctx->stats_in_end;
    FIRSTPASS_STATS  fps;
    // The minimum error here insures some bit allocation to frames even
    // in static regions. The allocation per MB declines for larger formats
    // where the typical "real" energy per MB also falls.
    // Initial estimate here uses sqrt(mbs) to define the min_err, where the
    // number of mbs is proportional to the image area.
    const int num_mbs = mb_rows * mb_cols;
    //(cpi->oxcf.resize_cfg.resize_mode != RESIZE_NONE)
    //    ? cpi->initial_mbs
    //    : mi_params->MBs;
    const double min_err = 200 * sqrt(num_mbs);

    fps.weight                   = stats->intra_factor * stats->brightness_factor;
    fps.frame                    = frame_number;
    fps.coded_error              = (double)(stats->coded_error >> 8) + min_err;
    fps.sr_coded_error           = (double)(stats->sr_coded_error >> 8) + min_err;
    fps.tr_coded_error           = (double)(stats->tr_coded_error >> 8) + min_err;
    fps.intra_error              = (double)(stats->intra_error >> 8) + min_err;
    fps.frame_avg_wavelet_energy = (double)stats->frame_avg_wavelet_energy;
    fps.count                    = 1.0;
    fps.pcnt_inter               = (double)stats->inter_count / num_mbs;
    fps.pcnt_second_ref          = (double)stats->second_ref_count / num_mbs;
    fps.pcnt_third_ref           = (double)stats->third_ref_count / num_mbs;
    fps.pcnt_neutral             = (double)stats->neutral_count / num_mbs;
    fps.intra_skip_pct           = (double)stats->intra_skip_count / num_mbs;
    fps.inactive_zone_rows       = (double)stats->image_data_start_row;
    fps.inactive_zone_cols       = (double)0; // TODO(paulwilkins): fix
    fps.raw_error_stdev          = raw_err_stdev;

    if (stats->mv_count > 0) {
        fps.MVr     = (double)stats->sum_mvr / stats->mv_count;
        fps.mvr_abs = (double)stats->sum_mvr_abs / stats->mv_count;
        fps.MVc     = (double)stats->sum_mvc / stats->mv_count;
        fps.mvc_abs = (double)stats->sum_mvc_abs / stats->mv_count;
        fps.MVrv    = ((double)stats->sum_mvrs -
                    ((double)stats->sum_mvr * stats->sum_mvr / stats->mv_count)) /
                   stats->mv_count;
        fps.MVcv = ((double)stats->sum_mvcs -
                    ((double)stats->sum_mvc * stats->sum_mvc / stats->mv_count)) /
                   stats->mv_count;
        fps.mv_in_out_count = (double)stats->sum_in_vectors / (stats->mv_count * 2);
        fps.new_mv_count    = stats->new_mv_count;
        fps.pcnt_motion     = (double)stats->mv_count / num_mbs;
    } else {
        fps.MVr             = 0.0;
        fps.mvr_abs         = 0.0;
        fps.MVc             = 0.0;
        fps.mvc_abs         = 0.0;
        fps.MVrv            = 0.0;
        fps.MVcv            = 0.0;
        fps.mv_in_out_count = 0.0;
        fps.new_mv_count    = 0.0;
        fps.pcnt_motion     = 0.0;
    }

    // TODO(paulwilkins):  Handle the case when duration is set to 0, or
    // something less than the full time between subsequent values of
    // cpi->source_time_stamp.
    fps.duration = (double)ts_duration;

    // We will store the stats inside the persistent twopass struct (and NOT the
    // local variable 'fps'), and then cpi->output_pkt_list will point to it.
    *this_frame_stats = fps;
    output_stats(scs_ptr, this_frame_stats, pcs_ptr->picture_number);
    if (twopass->stats_buf_ctx->total_stats != NULL) {
        svt_av1_accumulate_stats(twopass->stats_buf_ctx->total_stats, &fps);
    }
    /*In the case of two pass, first pass uses it as a circular buffer,
   * when LAP is enabled it is used as a linear buffer*/
    twopass->stats_buf_ctx->stats_in_end++;
    if ((use_output_stat(scs_ptr)) &&
        (twopass->stats_buf_ctx->stats_in_end >= twopass->stats_buf_ctx->stats_in_buf_end)) {
        twopass->stats_buf_ctx->stats_in_end = twopass->stats_buf_ctx->stats_in_start;
    }
}

static FRAME_STATS accumulate_frame_stats(FRAME_STATS *mb_stats, int mb_rows, int mb_cols) {
    FRAME_STATS stats = {0};
    int         i, j;

    stats.image_data_start_row = INVALID_ROW;
    for (j = 0; j < mb_rows; j++) {
        for (i = 0; i < mb_cols; i++) {
            FRAME_STATS mb_stat = mb_stats[j * mb_cols + i];
            stats.brightness_factor += mb_stat.brightness_factor;
            stats.coded_error += mb_stat.coded_error;
            stats.frame_avg_wavelet_energy += mb_stat.frame_avg_wavelet_energy;
            if (stats.image_data_start_row == INVALID_ROW &&
                mb_stat.image_data_start_row != INVALID_ROW) {
                stats.image_data_start_row = mb_stat.image_data_start_row;
            }
            stats.inter_count += mb_stat.inter_count;
            stats.intra_error += mb_stat.intra_error;
            stats.intra_factor += mb_stat.intra_factor;
            stats.intra_skip_count += mb_stat.intra_skip_count;
            stats.mv_count += mb_stat.mv_count;
            stats.neutral_count += mb_stat.neutral_count;
            stats.new_mv_count += mb_stat.new_mv_count;
            stats.second_ref_count += mb_stat.second_ref_count;
            stats.sr_coded_error += mb_stat.sr_coded_error;
            stats.sum_in_vectors += mb_stat.sum_in_vectors;
            stats.sum_mvc += mb_stat.sum_mvc;
            stats.sum_mvc_abs += mb_stat.sum_mvc_abs;
            stats.sum_mvcs += mb_stat.sum_mvcs;
            stats.sum_mvr += mb_stat.sum_mvr;
            stats.sum_mvr_abs += mb_stat.sum_mvr_abs;
            stats.sum_mvrs += mb_stat.sum_mvrs;
            stats.third_ref_count += mb_stat.third_ref_count;
            stats.tr_coded_error += mb_stat.tr_coded_error;
        }
    }
    return stats;
}
/**************************************************
* average_non_16x16_stats
* Handle stat for non 16x16 blocks. For non 16x16 blocks, some of the stats are increased multiple times
* First find the last block in the 16x16 area and then devide the stats by the number of small blocks
 **************************************************/
// Handle stat for non 16x16 blocks. For non 16x16 blocks, some of the stats are increased multiple times
// First find the last block in the 16x16 area and then devide the stats by the number of small blocks
void average_non_16x16_stats(FRAME_STATS *mb_stats, int blk_num) {
    mb_stats->brightness_factor /= blk_num;
    mb_stats->frame_avg_wavelet_energy /= blk_num;
    mb_stats->inter_count /= blk_num;
    mb_stats->intra_skip_count /= blk_num;
    mb_stats->mv_count /= blk_num;
    mb_stats->neutral_count /= blk_num;
    mb_stats->new_mv_count /= blk_num;
    mb_stats->second_ref_count /= blk_num;
    mb_stats->sum_in_vectors /= blk_num;
    mb_stats->sum_mvc /= blk_num;
    mb_stats->sum_mvc_abs /= blk_num;
    mb_stats->sum_mvcs /= blk_num;
    mb_stats->sum_mvr /= blk_num;
    mb_stats->sum_mvr_abs /= blk_num;
    mb_stats->sum_mvrs /= blk_num;
    mb_stats->third_ref_count /= blk_num;
    mb_stats->intra_factor /= blk_num;
}
/**************************************************
 * Reset first pass stat
 **************************************************/
void setup_firstpass_data(PictureParentControlSet *pcs_ptr) {
    SequenceControlSet *scs_ptr        = pcs_ptr->scs_ptr;
    FirstPassData *     firstpass_data = &pcs_ptr->firstpass_data;
    const uint32_t      mb_cols        = (scs_ptr->seq_header.max_frame_width + 16 - 1) / 16;
    const uint32_t      mb_rows        = (scs_ptr->seq_header.max_frame_height + 16 - 1) / 16;
    const uint32_t      num_mbs        = mb_cols * mb_rows;
    memset(firstpass_data->mb_stats, 0, sizeof(*firstpass_data->mb_stats) * num_mbs);
    for (uint32_t i = 0; i < num_mbs; i++)
        firstpass_data->mb_stats[i].image_data_start_row = INVALID_ROW;
}

void first_pass_frame_end(PictureParentControlSet *pcs_ptr, const int64_t ts_duration) {
    SequenceControlSet *scs_ptr = pcs_ptr->scs_ptr;
    const uint32_t      mb_cols = (scs_ptr->seq_header.max_frame_width + 16 - 1) / 16;
    const uint32_t      mb_rows = (scs_ptr->seq_header.max_frame_height + 16 - 1) / 16;

    int *        raw_motion_err_list = pcs_ptr->firstpass_data.raw_motion_err_list;
    FRAME_STATS *mb_stats            = pcs_ptr->firstpass_data.mb_stats;

    FRAME_STATS  stats                      = accumulate_frame_stats(mb_stats, mb_rows, mb_cols);
    int          total_raw_motion_err_count = frame_is_intra_only(pcs_ptr) ? 0 : mb_rows * mb_cols;
    const double raw_err_stdev =
        raw_motion_error_stdev(raw_motion_err_list, total_raw_motion_err_count);
    // Clamp the image start to rows/2. This number of rows is discarded top
    // and bottom as dead data so rows / 2 means the frame is blank.
    if ((stats.image_data_start_row > (int)mb_rows / 2) ||
        (stats.image_data_start_row == INVALID_ROW)) {
        stats.image_data_start_row = mb_rows / 2;
    }
    // Exclude any image dead zone
    if (stats.image_data_start_row > 0) {
        stats.intra_skip_count =
            AOMMAX(0, stats.intra_skip_count - (stats.image_data_start_row * (int)mb_cols * 2));
    }
    const int num_mbs = mb_rows * mb_cols;
    /*(cpi->oxcf.resize_cfg.resize_mode != RESIZE_NONE)
        ? cpi->initial_mbs
        : mi_params->MBs;*/
    stats.intra_factor      = stats.intra_factor / (double)num_mbs;
    stats.brightness_factor = stats.brightness_factor / (double)num_mbs;
    update_firstpass_stats(
        pcs_ptr, &stats, raw_err_stdev, (const int)pcs_ptr->picture_number, ts_duration);
}
/******************************************************
* Derive Pre-Analysis settings for first pass
Input   : encoder mode and tune
Output  : Pre-Analysis signal(s)
******************************************************/
extern EbErrorType first_pass_signal_derivation_pre_analysis(SequenceControlSet *     scs_ptr,
                                                             PictureParentControlSet *pcs_ptr) {
    EbErrorType return_error = EB_ErrorNone;
    // Derive HME Flag
    pcs_ptr->enable_hme_flag        = 1;
    pcs_ptr->enable_hme_level0_flag = 1;
    pcs_ptr->enable_hme_level1_flag = 1;
    pcs_ptr->enable_hme_level2_flag = 1;

    //// Set here to allocate resources for the downsampled pictures used in HME (generated in PictureAnalysis)
    //// Will be later updated for SC/NSC in PictureDecisionProcess
    pcs_ptr->tf_enable_hme_flag                  = 0;
    pcs_ptr->tf_enable_hme_level0_flag           = 0;
    pcs_ptr->tf_enable_hme_level1_flag           = 0;
    pcs_ptr->tf_enable_hme_level2_flag           = 0;
    scs_ptr->seq_header.enable_intra_edge_filter = 0;
    scs_ptr->seq_header.pic_based_rate_est       = 0;
    scs_ptr->seq_header.enable_restoration       = 0;
    scs_ptr->seq_header.cdef_level/*enable_cdef*/= 0;
    scs_ptr->seq_header.enable_warped_motion     = 0;

    return return_error;
}
extern EbErrorType av1_intra_full_cost(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                                       struct ModeDecisionCandidateBuffer *candidate_buffer_ptr,
                                       BlkStruct *blk_ptr, uint64_t *y_distortion,
                                       uint64_t *cb_distortion, uint64_t *cr_distortion,
                                       uint64_t lambda, uint64_t *y_coeff_bits,
                                       uint64_t *cb_coeff_bits, uint64_t *cr_coeff_bits,
                                       BlockSize bsize);

extern EbErrorType av1_inter_full_cost(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                                       struct ModeDecisionCandidateBuffer *candidate_buffer_ptr,
                                       BlkStruct *blk_ptr, uint64_t *y_distortion,
                                       uint64_t *cb_distortion, uint64_t *cr_distortion,
                                       uint64_t lambda, uint64_t *y_coeff_bits,
                                       uint64_t *cb_coeff_bits, uint64_t *cr_coeff_bits,
                                       BlockSize bsize);
void first_pass_perform_tx_partitioning(ModeDecisionCandidateBuffer *candidate_buffer,
    ModeDecisionContext *context_ptr, PictureControlSet *pcs_ptr,
    uint8_t start_tx_depth, uint8_t end_tx_depth,
    uint64_t *y_coeff_bits,
    uint64_t *y_full_distortion);
extern const EbPredictionFunc svt_product_prediction_fun_table[3];

extern void first_pass_loop_core(PictureControlSet *pcs_ptr,
                                 ModeDecisionContext *        context_ptr,
                                 ModeDecisionCandidateBuffer *candidate_buffer,
                                 ModeDecisionCandidate *      candidate_ptr,
                                 EbPictureBufferDesc *        input_picture_ptr,
                                 uint32_t input_origin_index, uint32_t blk_origin_index) {
    uint64_t y_full_distortion[DIST_CALC_TOTAL];

    uint64_t y_coeff_bits;

    int32_t is_inter = (candidate_buffer->candidate_ptr->type == INTER_MODE ||
                        candidate_buffer->candidate_ptr->use_intrabc)
                           ? EB_TRUE
                           : EB_FALSE;

    // initialize TU Split
    y_full_distortion[DIST_CALC_RESIDUAL]   = 0;
    y_full_distortion[DIST_CALC_PREDICTION] = 0;
    y_coeff_bits                            = 0;

    candidate_ptr->full_distortion = 0;
    memset(candidate_ptr->eob[0], 0, sizeof(uint16_t));
    memset(candidate_ptr->eob[1], 0, sizeof(uint16_t));
    memset(candidate_ptr->eob[2], 0, sizeof(uint16_t));

    candidate_ptr->chroma_distortion             = 0;
    candidate_ptr->chroma_distortion_inter_depth = 0;
    // Set Skip Flag
    candidate_ptr->skip_flag = EB_FALSE;
    if (is_inter)
    svt_product_prediction_fun_table[candidate_ptr->type](
        context_ptr->hbd_mode_decision, context_ptr, pcs_ptr, candidate_buffer);

    // Initialize luma CBF
    candidate_ptr->y_has_coeff = 0;
    candidate_ptr->u_has_coeff = 0;
    candidate_ptr->v_has_coeff = 0;

    // Initialize tx type
    for (int tu_index = 0; tu_index < MAX_TXB_COUNT; tu_index++)
        candidate_ptr->transform_type[tu_index] = DCT_DCT;
    uint8_t start_tx_depth = 0;
    uint8_t end_tx_depth   = 0;
    if (context_ptr->md_tx_size_search_mode == 0) {
        start_tx_depth = end_tx_depth;
    } else if (context_ptr->md_staging_tx_size_mode == 0) {
        start_tx_depth = end_tx_depth = candidate_buffer->candidate_ptr->tx_depth;
    }
    //Y Residual: residual for INTRA is computed inside the TU loop
    if (is_inter)
        //Y Residual
        residual_kernel(input_picture_ptr->buffer_y,
                        input_origin_index,
                        input_picture_ptr->stride_y,
                        candidate_buffer->prediction_ptr->buffer_y,
                        blk_origin_index,
                        candidate_buffer->prediction_ptr->stride_y,
                        (int16_t *)candidate_buffer->residual_ptr->buffer_y,
                        blk_origin_index,
                        candidate_buffer->residual_ptr->stride_y,
                        context_ptr->hbd_mode_decision,
                        context_ptr->blk_geom->bwidth,
                        context_ptr->blk_geom->bheight);
    first_pass_perform_tx_partitioning(
        candidate_buffer,
        context_ptr,
        pcs_ptr,
        start_tx_depth,
        end_tx_depth,
        &y_coeff_bits,
        &y_full_distortion[0]);
    candidate_ptr->chroma_distortion_inter_depth = 0;
    candidate_ptr->chroma_distortion             = 0;

    candidate_ptr->block_has_coeff =
        (candidate_ptr->y_has_coeff | candidate_ptr->u_has_coeff | candidate_ptr->v_has_coeff)
            ? EB_TRUE
            : EB_FALSE;

    //ALL PLANE
    *(candidate_buffer->full_cost_ptr) = 0;
}
#define LOW_MOTION_ERROR_THRESH 25
/***************************************************************************
* Computes and returns the intra pred error of a block.
* intra pred error: sum of squared error of the intra predicted residual.
* Modifies:
*   stats->intra_skip_count
*   stats->image_data_start_row
*   stats->intra_factor
*   stats->brightness_factor
*   stats->intra_error
*   stats->frame_avg_wavelet_energy
* Returns:
*   this_intra_error.
***************************************************************************/
static int firstpass_intra_prediction(PictureControlSet *pcs_ptr,
                                      ModeDecisionContext *        context_ptr,
                                      ModeDecisionCandidateBuffer *candidate_buffer,
                                      ModeDecisionCandidate *      candidate_ptr,
                                      EbPictureBufferDesc *        input_picture_ptr,
                                      uint32_t input_origin_index, uint32_t blk_origin_index,
                                      FRAME_STATS *const stats) {
    int32_t         mb_row      = context_ptr->blk_origin_y >> 4;
    int32_t         mb_col      = context_ptr->blk_origin_x >> 4;
    const int       use_dc_pred = (mb_col || mb_row) && (!mb_col || !mb_row);
    const BlockSize bsize       = context_ptr->blk_geom->bsize;

    // Initialize tx_depth
    candidate_buffer->candidate_ptr->tx_depth =
        use_dc_pred ? 0 : (bsize == BLOCK_16X16 ? 2 : bsize == BLOCK_8X8 ? 1 : 0);
    candidate_buffer->candidate_ptr->fast_luma_rate   = 0;
    candidate_buffer->candidate_ptr->fast_chroma_rate = 0;
    context_ptr->md_staging_skip_interpolation_search = EB_TRUE;
    context_ptr->md_staging_skip_chroma_pred          = EB_FALSE;
    context_ptr->md_staging_tx_size_mode              = 0;
    context_ptr->md_staging_skip_full_chroma          = EB_FALSE;
    context_ptr->md_staging_skip_rdoq                 = EB_TRUE;
    context_ptr->md_staging_spatial_sse_full_loop_level     = context_ptr->spatial_sse_full_loop_level;

    first_pass_loop_core(pcs_ptr,
                         context_ptr,
                         candidate_buffer,
                         candidate_ptr,
                         input_picture_ptr,
                         input_origin_index,
                         blk_origin_index);

    EbSpatialFullDistType spatial_full_dist_type_fun = context_ptr->hbd_mode_decision
                                                           ? svt_full_distortion_kernel16_bits
                                                           : svt_spatial_full_distortion_kernel;

    int this_intra_error =
        (uint32_t)(spatial_full_dist_type_fun(input_picture_ptr->buffer_y,
                                              input_origin_index,
                                              input_picture_ptr->stride_y,
                                              candidate_buffer->prediction_ptr->buffer_y,
                                              blk_origin_index,
                                              candidate_buffer->prediction_ptr->stride_y,
                                              context_ptr->blk_geom->bwidth,
                                              context_ptr->blk_geom->bheight));

    if (this_intra_error < UL_INTRA_THRESH) {
        ++stats->intra_skip_count;
    } else if ((mb_col > 0) && (stats->image_data_start_row == INVALID_ROW)) {
        stats->image_data_start_row = mb_row;
    }

    if (context_ptr->hbd_mode_decision) {
        switch (pcs_ptr->parent_pcs_ptr->av1_cm->bit_depth) {
        case AOM_BITS_8: break;
        case AOM_BITS_10: this_intra_error >>= 4; break;
        case AOM_BITS_12: this_intra_error >>= 8; break;
        default:
            assert(0 &&
                   "seq_params->bit_depth should be AOM_BITS_8, "
                   "AOM_BITS_10 or AOM_BITS_12");
            return -1;
        }
    }

    // aom_clear_system_state();
    double log_intra = log1p((double)this_intra_error);
    if (log_intra < 10.0)
        stats->intra_factor += 1.0 + ((10.0 - log_intra) * 0.05);
    else
        stats->intra_factor += 1.0;

    int level_sample;
    if (context_ptr->hbd_mode_decision)
        level_sample = ((uint16_t *)input_picture_ptr->buffer_y)[input_origin_index];
    else
        level_sample = input_picture_ptr->buffer_y[input_origin_index];

    if ((level_sample < DARK_THRESH) && (log_intra < 9.0))
        stats->brightness_factor += 1.0 + (0.01 * (DARK_THRESH - level_sample));
    else
        stats->brightness_factor += 1.0;
    // Intrapenalty below deals with situations where the intra and inter
    // error scores are very low (e.g. a plain black frame).
    // We do not have special cases in first pass for 0,0 and nearest etc so
    // all inter modes carry an overhead cost estimate for the mv.
    // When the error score is very low this causes us to pick all or lots of
    // INTRA modes and throw lots of key frames.
    // This penalty adds a cost matching that of a 0,0 mv to the intra case.
    this_intra_error += INTRA_MODE_PENALTY;

    const int hbd    = context_ptr->hbd_mode_decision;
    const int stride = input_picture_ptr->stride_y;
    if (hbd) {
        uint16_t *buf = &((uint16_t *)input_picture_ptr->buffer_y)[input_origin_index];
        for (int r8 = 0; r8 < 2; ++r8) {
            for (int c8 = 0; c8 < 2; ++c8) {
                stats->frame_avg_wavelet_energy += svt_av1_haar_ac_sad_8x8_uint8_input(
                    CONVERT_TO_BYTEPTR(buf) + c8 * 8 + r8 * 8 * stride, stride, hbd);
            }
        }
    } else {
        uint8_t *buf = &(input_picture_ptr->buffer_y)[input_origin_index];
        for (int r8 = 0; r8 < 2; ++r8) {
            for (int c8 = 0; c8 < 2; ++c8) {
                stats->frame_avg_wavelet_energy +=
                    svt_av1_haar_ac_sad_8x8_uint8_input(buf + c8 * 8 + r8 * 8 * stride, stride, hbd);
            }
        }
    }
    // Accumulate the intra error.
    stats->intra_error += (int64_t)this_intra_error;
    return this_intra_error;
}
/***************************************************************************
* Computes and returns the inter prediction error from the last frame.
* Computes inter prediction errors from the golden and alt ref frames and
* Updates stats accordingly.
* Modifies:
*    stats: many member params in it.
*  Returns:
*    this_inter_error
***************************************************************************/
static int firstpass_inter_prediction(
    PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
    ModeDecisionCandidateBuffer *candidate_buffer, ModeDecisionCandidate *candidate_ptr,
    EbPictureBufferDesc *input_picture_ptr, uint32_t input_origin_index, uint32_t blk_origin_index,
    uint32_t fast_candidate_total_count, const int this_intra_error,
    /*int *raw_motion_err_list, */ FRAME_STATS *stats) {

    int32_t        mb_row = context_ptr->blk_origin_y >> 4;
    int32_t        mb_col = context_ptr->blk_origin_x >> 4;
    const uint32_t mb_cols =
        (pcs_ptr->parent_pcs_ptr->scs_ptr->seq_header.max_frame_width + 16 - 1) / 16;
    const uint32_t mb_rows =
        (pcs_ptr->parent_pcs_ptr->scs_ptr->seq_header.max_frame_height + 16 - 1) / 16;
    int this_inter_error = this_intra_error;
    FULLPEL_MV mv = kZeroFullMv;
    MV last_mv;

    uint32_t full_lambda = context_ptr->hbd_mode_decision
                               ? context_ptr->full_lambda_md[EB_10_BIT_MD]
                               : context_ptr->full_lambda_md[EB_8_BIT_MD];
    int errorperbit = full_lambda >> RD_EPB_SHIFT;
    errorperbit += (errorperbit == 0);
    EbSpatialFullDistType spatial_full_dist_type_fun = context_ptr->hbd_mode_decision
                                                           ? svt_full_distortion_kernel16_bits
                                                           : svt_spatial_full_distortion_kernel;

    int motion_error = 0;
    // TODO(pengchong): Replace the hard-coded threshold
    // anaghdin: to add support
    if (1) //(raw_motion_error > LOW_MOTION_ERROR_THRESH)
    {
        uint32_t                      cand_index = 1;
        ModeDecisionCandidateBuffer **candidate_buffer_ptr_array_base =
            context_ptr->candidate_buffer_ptr_array;
        ModeDecisionCandidateBuffer **candidate_buffer_ptr_array =
            &(candidate_buffer_ptr_array_base[0]);

        candidate_buffer = candidate_buffer_ptr_array[cand_index];
        candidate_ptr    = candidate_buffer->candidate_ptr =
            &context_ptr->fast_candidate_array[cand_index];
        context_ptr->best_candidate_index_array[cand_index] = cand_index;
        // Initialize tx_depth
        candidate_buffer->candidate_ptr->tx_depth = 0;
        candidate_buffer->candidate_ptr->fast_luma_rate   = 0;
        candidate_buffer->candidate_ptr->fast_chroma_rate = 0;
        candidate_buffer->candidate_ptr->interp_filters   = 0;
        svt_product_prediction_fun_table[candidate_ptr->type](
            context_ptr->hbd_mode_decision, context_ptr, pcs_ptr, candidate_buffer);
        *(candidate_buffer->full_cost_ptr) = 0;
        // To convert full-pel MV
        mv.col = candidate_buffer->candidate_ptr->motion_vector_xl0 >> 3;
        mv.row = candidate_buffer->candidate_ptr->motion_vector_yl0 >> 3;

        last_mv.col = context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds].ref_mvs[1][0].as_mv.col;
        last_mv.row = context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds].ref_mvs[1][0].as_mv.row;

        motion_error =
            (uint32_t)(spatial_full_dist_type_fun(input_picture_ptr->buffer_y,
                                                  input_origin_index,
                                                  input_picture_ptr->stride_y,
                                                  candidate_buffer->prediction_ptr->buffer_y,
                                                  blk_origin_index,
                                                  candidate_buffer->prediction_ptr->stride_y,
                                                  context_ptr->blk_geom->bwidth,
                                                  context_ptr->blk_geom->bheight));

        // Assume 0,0 motion with no mv overhead.
        if (mv.col != 0 && mv.row != 0) {
            const MV temp_full_mv = get_mv_from_fullmv(&mv);
            motion_error += mv_err_cost(&temp_full_mv,
                                        &last_mv,
                                        context_ptr->md_rate_estimation_ptr->nmv_vec_cost,
                                        context_ptr->md_rate_estimation_ptr->nmvcoststack,
                                        errorperbit) +
                            NEW_MV_MODE_PENALTY;
        }

        // Motion search in 2nd reference frame.
        int gf_motion_error = motion_error;
        if (fast_candidate_total_count > 2) {
            cand_index++;
            candidate_buffer = candidate_buffer_ptr_array[cand_index];
            candidate_ptr    = candidate_buffer->candidate_ptr =
                &context_ptr->fast_candidate_array[cand_index];
            context_ptr->best_candidate_index_array[cand_index] = cand_index;
            // Initialize tx_depth
            candidate_buffer->candidate_ptr->tx_depth = 0;
            candidate_buffer->candidate_ptr->interp_filters   = 0;
            svt_product_prediction_fun_table[candidate_ptr->type](
                context_ptr->hbd_mode_decision, context_ptr, pcs_ptr, candidate_buffer);

            gf_motion_error =
                (uint32_t)(spatial_full_dist_type_fun(input_picture_ptr->buffer_y,
                                                      input_origin_index,
                                                      input_picture_ptr->stride_y,
                                                      candidate_buffer->prediction_ptr->buffer_y,
                                                      blk_origin_index,
                                                      candidate_buffer->prediction_ptr->stride_y,
                                                      context_ptr->blk_geom->bwidth,
                                                      context_ptr->blk_geom->bheight));
            // To convert full-pel MV
            FULLPEL_MV gf_mv;
            gf_mv.col = candidate_buffer->candidate_ptr->motion_vector_xl1 >> 3;
            gf_mv.row = candidate_buffer->candidate_ptr->motion_vector_yl1 >> 3;

            // Assume 0,0 motion with no mv overhead.
            if (gf_mv.col != 0 && gf_mv.row != 0) {
                const MV temp_full_mv = get_mv_from_fullmv(&gf_mv);
                gf_motion_error += mv_err_cost(&temp_full_mv,
                                               &kZeroMv,
                                               context_ptr->md_rate_estimation_ptr->nmv_vec_cost,
                                               context_ptr->md_rate_estimation_ptr->nmvcoststack,
                                               errorperbit) +
                                   NEW_MV_MODE_PENALTY;
            }
        }

        if (gf_motion_error < motion_error && gf_motion_error < this_intra_error) {
            ++stats->second_ref_count;
        }
        // In accumulating a score for the 2nd reference frame take the
        // best of the motion predicted score and the intra coded error
        // (just as will be done for) accumulation of "coded_error" for
        // the last frame.
        if (fast_candidate_total_count > 2) {
            stats->sr_coded_error += AOMMIN(gf_motion_error, this_intra_error);
        } else {
            stats->sr_coded_error += motion_error;
        }

        // Motion search in 3rd reference frame.
        int alt_motion_error = motion_error;
        if (alt_motion_error < motion_error && alt_motion_error < gf_motion_error &&
            alt_motion_error < this_intra_error) {
            ++stats->third_ref_count;
        }
        // In accumulating a score for the 3rd reference frame take the
        // best of the motion predicted score and the intra coded error
        // (just as will be done for) accumulation of "coded_error" for
        // the last frame.
        // alt_ref_frame is not supported yet
        if (0 /*alt_ref_frame != NULL*/) {
            stats->tr_coded_error += AOMMIN(alt_motion_error, this_intra_error);
        } else {
            stats->tr_coded_error += motion_error;
        }
    } else {
        stats->sr_coded_error += motion_error;
        stats->tr_coded_error += motion_error;
    }

    // Start by assuming that intra mode is best.
    if (motion_error <= this_intra_error) {
#ifdef ARCH_X86_64
        aom_clear_system_state();
#endif
        // Keep a count of cases where the inter and intra were very close
        // and very low. This helps with scene cut detection for example in
        // cropped clips with black bars at the sides or top and bottom.
        if (((this_intra_error - INTRA_MODE_PENALTY) * 9 <= motion_error * 10) &&
            (this_intra_error < (2 * INTRA_MODE_PENALTY))) {
            stats->neutral_count += 1.0;
            // Also track cases where the intra is not much worse than the inter
            // and use this in limiting the GF/arf group length.
        } else if ((this_intra_error > NCOUNT_INTRA_THRESH) &&
                   (this_intra_error < (NCOUNT_INTRA_FACTOR * motion_error))) {
            stats->neutral_count +=
                (double)motion_error / DOUBLE_DIVIDE_CHECK((double)this_intra_error);
        }
        const MV best_mv = get_mv_from_fullmv(&mv);
        this_inter_error = motion_error;
        stats->sum_mvr += best_mv.row;
        stats->sum_mvr_abs += abs(best_mv.row);
        stats->sum_mvc += best_mv.col;
        stats->sum_mvc_abs += abs(best_mv.col);
        stats->sum_mvrs += best_mv.row * best_mv.row;
        stats->sum_mvcs += best_mv.col * best_mv.col;
        ++stats->inter_count;
        accumulate_mv_stats(best_mv, mv, mb_row, mb_col, mb_rows, mb_cols, &last_mv, stats);
    }

    return this_inter_error;
}

void set_inter_comp_controls(ModeDecisionContext *mdctxt, uint8_t inter_comp_mode);
void set_dist_based_ref_pruning_controls(ModeDecisionContext *mdctxt, uint8_t dist_based_ref_pruning_level);

/******************************************************
* Derive md Settings(feature signals) that could be
  changed  at the block level
******************************************************/
extern EbErrorType first_pass_signal_derivation_block(
    ModeDecisionContext   *context_ptr) {

    EbErrorType return_error = EB_ErrorNone;


    // Set dist_based_ref_pruning
    context_ptr->dist_based_ref_pruning = 0;

    set_dist_based_ref_pruning_controls(context_ptr, context_ptr->dist_based_ref_pruning);

    // set compound_types_to_try
    set_inter_comp_controls(context_ptr, 0);

    context_ptr->compound_types_to_try = MD_COMP_AVG;

    // Do not add MD_COMP_WEDGE  beyond this point
    if (get_wedge_params_bits(context_ptr->blk_geom->bsize) == 0)
        context_ptr->compound_types_to_try = MIN(context_ptr->compound_types_to_try, MD_COMP_DIFF0);
    context_ptr->inject_inter_candidates = 1;

    return return_error;
}

void product_coding_loop_init_fast_loop(ModeDecisionContext *context_ptr,
                                        NeighborArrayUnit *  skip_coeff_neighbor_array,
                                        NeighborArrayUnit *  inter_pred_dir_neighbor_array,
                                        NeighborArrayUnit *  ref_frame_type_neighbor_array,
                                        NeighborArrayUnit *  intra_luma_mode_neighbor_array,
                                        NeighborArrayUnit *  skip_flag_neighbor_array,
                                        NeighborArrayUnit *  mode_type_neighbor_array,
                                        NeighborArrayUnit *  leaf_depth_neighbor_array,
                                        NeighborArrayUnit *  leaf_partition_neighbor_array);
// inject intra candidates for first pass
void  first_pass_inject_intra_candidates(
    ModeDecisionContext          *context_ptr,
    uint32_t                     *candidate_total_cnt){


    uint32_t                    cand_total_cnt = 0;
    ModeDecisionCandidate    *cand_array = context_ptr->fast_candidate_array;

    cand_array[cand_total_cnt].type = INTRA_MODE;
    cand_array[cand_total_cnt].merge_flag = EB_FALSE;

    cand_array[cand_total_cnt].palette_info = NULL;

    cand_array[cand_total_cnt].intra_luma_mode = DC_PRED;
    cand_array[cand_total_cnt].distortion_ready = 0;
    cand_array[cand_total_cnt].use_intrabc = 0;
    cand_array[cand_total_cnt].filter_intra_mode = FILTER_INTRA_MODES;
    cand_array[cand_total_cnt].is_directional_mode_flag = 0;
    cand_array[cand_total_cnt].angle_delta[PLANE_TYPE_Y] = 0;

    cand_array[cand_total_cnt].intra_chroma_mode = UV_DC_PRED;

    cand_array[cand_total_cnt].is_directional_chroma_mode_flag = 0;
    cand_array[cand_total_cnt].angle_delta[PLANE_TYPE_UV] = 0;

    cand_array[cand_total_cnt].cfl_alpha_signs = 0;
    cand_array[cand_total_cnt].cfl_alpha_idx = 0;
    cand_array[cand_total_cnt].transform_type[0] = DCT_DCT;
    cand_array[cand_total_cnt].transform_type_uv = DCT_DCT;

    cand_array[cand_total_cnt].ref_frame_type = INTRA_FRAME;
    cand_array[cand_total_cnt].pred_mode = DC_PRED;
    cand_array[cand_total_cnt].motion_mode = SIMPLE_TRANSLATION;
    cand_array[cand_total_cnt].is_interintra_used = 0;
    INCRMENT_CAND_TOTAL_COUNT(cand_total_cnt);

    // update the total number of candidates injected
    (*candidate_total_cnt) = cand_total_cnt;

    return;
}
void choose_best_av1_mv_pred(ModeDecisionContext *           context_ptr,
                             struct MdRateEstimationContext *md_rate_estimation_ptr,
                             BlkStruct *blk_ptr, MvReferenceFrame ref_frame, uint8_t is_compound,
                             PredictionMode mode, //NEW or NEW_NEW
                             int16_t mv0x, int16_t mv0y, int16_t mv1x, int16_t mv1y,
                             uint8_t *bestDrlIndex, // output
                             IntMv    best_pred_mv[2] ) ;// output
EbBool mrp_is_already_injected_mv_l0(ModeDecisionContext *context_ptr, int16_t mv_x, int16_t mv_y,
                                     uint8_t ref_type);

EbBool mrp_is_already_injected_mv_l1(ModeDecisionContext *context_ptr, int16_t mv_x, int16_t mv_y,
                                     uint8_t ref_type);

EbBool is_valid_unipred_ref(
    struct ModeDecisionContext *context_ptr,
    uint8_t inter_cand_group,
    uint8_t list_idx, uint8_t ref_idx);


// inject new candidates for first pass
void first_pass_inject_new_candidates(const SequenceControlSet *  scs_ptr,
                           struct ModeDecisionContext *context_ptr, PictureControlSet *pcs_ptr,
                           EbBool is_compound_enabled, uint32_t me_sb_addr,
                           uint32_t me_block_offset, uint32_t *candidate_total_cnt) {
    ModeDecisionCandidate *cand_array      = context_ptr->fast_candidate_array;
    uint32_t               cand_total_cnt  = (*candidate_total_cnt);

    const MeSbResults *me_results =
        pcs_ptr->parent_pcs_ptr->pa_me_data->me_results[me_sb_addr];

    uint8_t            total_me_cnt     = me_results->total_me_candidate_index[me_block_offset];

    const MeCandidate *me_block_results = &me_results->me_candidate_array[context_ptr->me_cand_offset];

    MacroBlockD *      xd               = context_ptr->blk_ptr->av1xd;
    int                inside_tile      = 1;
    int                umv0tile         = (scs_ptr->static_config.unrestricted_motion_vector == 0);
    uint32_t           mi_row           = context_ptr->blk_origin_y >> MI_SIZE_LOG2;
    uint32_t           mi_col           = context_ptr->blk_origin_x >> MI_SIZE_LOG2;

    for (uint8_t me_candidate_index = 0; me_candidate_index < total_me_cnt; ++me_candidate_index) {
        const MeCandidate *me_block_results_ptr = &me_block_results[me_candidate_index];
        const uint8_t      inter_direction      = me_block_results_ptr->direction;
        const uint8_t      list0_ref_index      = me_block_results_ptr->ref_idx_l0;
        const uint8_t      list1_ref_index      = me_block_results_ptr->ref_idx_l1;

        /**************
            NEWMV L0
        ************* */
        if (inter_direction == 0) {

            if (!is_valid_unipred_ref(context_ptr, MIN(TOT_INTER_GROUP-1,PA_ME_GROUP), REF_LIST_0, list0_ref_index))
                continue;

            if (list0_ref_index > context_ptr->md_max_ref_count - 1) continue;
            int16_t to_inject_mv_x =
                context_ptr
                ->sb_me_mv[context_ptr->blk_geom->blkidx_mds][REF_LIST_0][list0_ref_index][0];
            int16_t to_inject_mv_y =
                context_ptr
                ->sb_me_mv[context_ptr->blk_geom->blkidx_mds][REF_LIST_0][list0_ref_index][1];
            uint8_t to_inject_ref_type = svt_get_ref_frame_type(REF_LIST_0, list0_ref_index);

            inside_tile = 1;
            if (umv0tile)
                inside_tile = is_inside_tile_boundary(&(xd->tile),
                                                      to_inject_mv_x,
                                                      to_inject_mv_y,
                                                      mi_col,
                                                      mi_row,
                                                      context_ptr->blk_geom->bsize);
            uint8_t skip_cand =  (!inside_tile);

            if (!skip_cand &&
                (context_ptr->injected_mv_count_l0 == 0 ||
                 mrp_is_already_injected_mv_l0(
                     context_ptr, to_inject_mv_x, to_inject_mv_y, to_inject_ref_type) ==
                     EB_FALSE)) {

                cand_array[cand_total_cnt].type                    = INTER_MODE;
                cand_array[cand_total_cnt].distortion_ready        = 0;
                cand_array[cand_total_cnt].use_intrabc             = 0;
                cand_array[cand_total_cnt].merge_flag              = EB_FALSE;
                cand_array[cand_total_cnt].prediction_direction[0] = (EbPredDirection)0;
                cand_array[cand_total_cnt].inter_mode              = NEWMV;
                cand_array[cand_total_cnt].pred_mode               = NEWMV;
                cand_array[cand_total_cnt].motion_mode             = SIMPLE_TRANSLATION;
                cand_array[cand_total_cnt].is_compound             = 0;
                cand_array[cand_total_cnt].is_new_mv               = 1;
                cand_array[cand_total_cnt].drl_index               = 0;

                // Set the MV to ME result
                cand_array[cand_total_cnt].motion_vector_xl0 = to_inject_mv_x;
                cand_array[cand_total_cnt].motion_vector_yl0 = to_inject_mv_y;

                // will be needed later by the rate estimation
                cand_array[cand_total_cnt].ref_mv_index   = 0;
                cand_array[cand_total_cnt].ref_frame_type =
                    svt_get_ref_frame_type(REF_LIST_0, list0_ref_index);
                cand_array[cand_total_cnt].ref_frame_index_l0 = list0_ref_index;
                cand_array[cand_total_cnt].ref_frame_index_l1 = -1;

                cand_array[cand_total_cnt].transform_type[0] = DCT_DCT;
                cand_array[cand_total_cnt].transform_type_uv = DCT_DCT;

                cand_array[cand_total_cnt].motion_vector_pred_x[REF_LIST_0] = 0;
                cand_array[cand_total_cnt].motion_vector_pred_y[REF_LIST_0] = 0;

                    cand_array[cand_total_cnt].is_interintra_used = 0;
                    cand_array[cand_total_cnt].motion_mode        = SIMPLE_TRANSLATION;

                    INCRMENT_CAND_TOTAL_COUNT(cand_total_cnt);

                context_ptr->injected_mv_x_l0_array[context_ptr->injected_mv_count_l0] =
                    to_inject_mv_x;
                context_ptr->injected_mv_y_l0_array[context_ptr->injected_mv_count_l0] =
                    to_inject_mv_y;
                context_ptr->injected_ref_type_l0_array[context_ptr->injected_mv_count_l0] =
                    to_inject_ref_type;
                ++context_ptr->injected_mv_count_l0;

            }
        }

        if (is_compound_enabled) {
            /**************
               NEWMV L1
           ************* */
            if (inter_direction == 1) {
                if (!is_valid_unipred_ref(context_ptr, MIN(TOT_INTER_GROUP-1,PA_ME_GROUP), REF_LIST_1, list1_ref_index))
                    continue;

                if (list1_ref_index > context_ptr->md_max_ref_count - 1) continue;
                int16_t to_inject_mv_x = context_ptr->sb_me_mv[context_ptr->blk_geom->blkidx_mds]
                    [REF_LIST_1][list1_ref_index][0];
                int16_t to_inject_mv_y = context_ptr->sb_me_mv[context_ptr->blk_geom->blkidx_mds]
                    [REF_LIST_1][list1_ref_index][1];
                uint8_t to_inject_ref_type = svt_get_ref_frame_type(REF_LIST_1, list1_ref_index);

                inside_tile = 1;
                if (umv0tile)
                    inside_tile = is_inside_tile_boundary(&(xd->tile),
                                                          to_inject_mv_x,
                                                          to_inject_mv_y,
                                                          mi_col,
                                                          mi_row,
                                                          context_ptr->blk_geom->bsize);
                uint8_t skip_cand =  !inside_tile;

                if (!skip_cand &&
                    (context_ptr->injected_mv_count_l1 == 0 ||
                     mrp_is_already_injected_mv_l1(
                         context_ptr, to_inject_mv_x, to_inject_mv_y, to_inject_ref_type) ==
                         EB_FALSE)) {

                    cand_array[cand_total_cnt].type                    = INTER_MODE;
                    cand_array[cand_total_cnt].distortion_ready        = 0;
                    cand_array[cand_total_cnt].use_intrabc             = 0;
                    cand_array[cand_total_cnt].merge_flag              = EB_FALSE;
                    cand_array[cand_total_cnt].prediction_direction[0] = (EbPredDirection)1;
                    cand_array[cand_total_cnt].inter_mode              = NEWMV;
                    cand_array[cand_total_cnt].pred_mode               = NEWMV;
                    cand_array[cand_total_cnt].motion_mode             = SIMPLE_TRANSLATION;
                    cand_array[cand_total_cnt].is_compound             = 0;
                    cand_array[cand_total_cnt].is_new_mv               = 1;
                    cand_array[cand_total_cnt].drl_index               = 0;

                    // Set the MV to ME result
                    cand_array[cand_total_cnt].motion_vector_xl1 = to_inject_mv_x;
                    cand_array[cand_total_cnt].motion_vector_yl1 = to_inject_mv_y;

                    // will be needed later by the rate estimation
                    cand_array[cand_total_cnt].ref_mv_index   = 0;
                    cand_array[cand_total_cnt].ref_frame_type =
                        svt_get_ref_frame_type(REF_LIST_1, list1_ref_index);
                    cand_array[cand_total_cnt].ref_frame_index_l0 = -1;
                    cand_array[cand_total_cnt].ref_frame_index_l1 = list1_ref_index;

                    cand_array[cand_total_cnt].transform_type[0] = DCT_DCT;
                    cand_array[cand_total_cnt].transform_type_uv = DCT_DCT;
                    cand_array[cand_total_cnt].motion_vector_pred_x[REF_LIST_1] = 0;
                    cand_array[cand_total_cnt].motion_vector_pred_y[REF_LIST_1] = 0;
                        cand_array[cand_total_cnt].is_interintra_used = 0;
                        cand_array[cand_total_cnt].motion_mode        = SIMPLE_TRANSLATION;

                    INCRMENT_CAND_TOTAL_COUNT(cand_total_cnt);

                    context_ptr->injected_mv_x_l1_array[context_ptr->injected_mv_count_l1] =
                        to_inject_mv_x;
                    context_ptr->injected_mv_y_l1_array[context_ptr->injected_mv_count_l1] =
                        to_inject_mv_y;
                    context_ptr->injected_ref_type_l1_array[context_ptr->injected_mv_count_l1] =
                        to_inject_ref_type;
                    ++context_ptr->injected_mv_count_l1;

                }
            }
        }
    }
    // update the total number of candidates injected
    (*candidate_total_cnt) = cand_total_cnt;
}

// inject inter candidates for first pass
void first_pass_inject_inter_candidates(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                             const SequenceControlSet *scs_ptr, uint32_t *candidate_total_cnt) {


    FrameHeader *          frm_hdr        = &pcs_ptr->parent_pcs_ptr->frm_hdr;
    uint32_t               cand_total_cnt = *candidate_total_cnt;

    EbBool       is_compound_enabled      = (frm_hdr->reference_mode == SINGLE_REFERENCE) ? 0 : 1;

    first_pass_inject_new_candidates(scs_ptr,
                          context_ptr,
                          pcs_ptr,
                          is_compound_enabled,
                          context_ptr->me_sb_addr,
                          context_ptr->me_block_offset,
                          &cand_total_cnt);

    // update the total number of candidates injected
    *candidate_total_cnt = cand_total_cnt;
}
void check_mv_validity(int16_t x_mv, int16_t y_mv, uint8_t need_shift);

// generate candidates for first pass
EbErrorType first_pass_generate_md_stage_0_cand(
    ModeDecisionContext *context_ptr,
    uint32_t            *candidate_total_count_ptr,
    PictureControlSet   *pcs_ptr)
{
    const SequenceControlSet *scs_ptr = (SequenceControlSet*)pcs_ptr->scs_wrapper_ptr->object_ptr;
    const EB_SLICE slice_type = pcs_ptr->slice_type;
    uint32_t cand_total_cnt = 0;
    // Reset duplicates variables
    context_ptr->injected_mv_count_l0 = 0;
    context_ptr->injected_mv_count_l1 = 0;
    context_ptr->injected_mv_count_bipred = 0;

    //----------------------
    // Intra
    if (context_ptr->blk_geom->sq_size < 128) {
                first_pass_inject_intra_candidates(
                    context_ptr,
                &cand_total_cnt);
    }
    if (slice_type != I_SLICE && context_ptr->inject_inter_candidates) {
            first_pass_inject_inter_candidates(
                pcs_ptr,
                context_ptr,
                scs_ptr,
                &cand_total_cnt);
    }
//#if INJECT_BACKUP_CANDIDATE
//    // For I_SLICE, DC is always injected, and therefore there is no a risk of no candidates @ md_syage_0()
//    // For non I_SLICE, there is a risk of no candidates @ md_stage_0() because of the INTER candidates pruning techniques
//    if (slice_type != I_SLICE && cand_total_cnt == 0) {
//        inject_zz_backup_candidate(
//            context_ptr,
//            &cand_total_cnt);
//    }
//#endif
    *candidate_total_count_ptr = cand_total_cnt;
    CandClass  cand_class_it;
    memset(context_ptr->md_stage_0_count, 0, CAND_CLASS_TOTAL * sizeof(uint32_t));

    uint32_t cand_i;
    for (cand_i = 0; cand_i < cand_total_cnt; cand_i++)
    {
        ModeDecisionCandidate * cand_ptr = &context_ptr->fast_candidate_array[cand_i];

        if (cand_ptr->type == INTRA_MODE) {
            // Intra prediction
            if (cand_ptr->palette_info == NULL ||
                    cand_ptr->palette_info->pmi.palette_size[0] == 0) {
            cand_ptr->cand_class = CAND_CLASS_0;
            context_ptr->md_stage_0_count[CAND_CLASS_0]++;
            }
            else {
                // Palette Prediction
                cand_ptr->cand_class = CAND_CLASS_3;
                context_ptr->md_stage_0_count[CAND_CLASS_3]++;
            }
        }
        else  if (cand_ptr->is_new_mv) {
            // MV Prediction
            cand_ptr->cand_class = CAND_CLASS_1;
            context_ptr->md_stage_0_count[CAND_CLASS_1]++;
        }
        else {
            //MVP Prediction
            cand_ptr->cand_class = CAND_CLASS_2;
            context_ptr->md_stage_0_count[CAND_CLASS_2]++;
        }
    }
    uint32_t fast_accum = 0;
    for (cand_class_it = CAND_CLASS_0; cand_class_it < CAND_CLASS_TOTAL; cand_class_it++) {
        fast_accum += context_ptr->md_stage_0_count[cand_class_it];
    }
    assert(fast_accum == cand_total_cnt);

    //check if final MV is within AV1 limits
    for (cand_i = 0; cand_i < cand_total_cnt; cand_i++)
    {
        ModeDecisionCandidate * cand_ptr = &context_ptr->fast_candidate_array[cand_i];

        if (cand_ptr->type == INTER_MODE) {
            if (cand_ptr->prediction_direction[0] == UNI_PRED_LIST_0 ||
                cand_ptr->prediction_direction[0] == BI_PRED)
                check_mv_validity(cand_ptr->motion_vector_xl0,
                    cand_ptr->motion_vector_yl0, 0);

            if (cand_ptr->prediction_direction[0] == UNI_PRED_LIST_1 ||
                cand_ptr->prediction_direction[0] == BI_PRED)
                check_mv_validity(cand_ptr->motion_vector_xl1,
                    cand_ptr->motion_vector_yl1, 0);

        }
    }

    return EB_ErrorNone;
}
void distortion_based_modulator(ModeDecisionContext *context_ptr,
    EbPictureBufferDesc *input_picture_ptr, uint32_t input_origin_index,
    EbPictureBufferDesc *recon_ptr, uint32_t blk_origin_index);

extern void first_pass_md_encode_block(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                                       EbPictureBufferDesc *        input_picture_ptr,
                                       ModeDecisionCandidateBuffer *bestcandidate_buffers[5]) {
    ModeDecisionCandidateBuffer **candidate_buffer_ptr_array_base =
        context_ptr->candidate_buffer_ptr_array;
    ModeDecisionCandidateBuffer **candidate_buffer_ptr_array;
    const BlockGeom *             blk_geom = context_ptr->blk_geom;
    ModeDecisionCandidateBuffer * candidate_buffer;
    ModeDecisionCandidate *       fast_candidate_array = context_ptr->fast_candidate_array;
    uint32_t                      candidate_index;
    uint32_t                      fast_candidate_total_count;
    uint32_t                      best_intra_mode = EB_INTRA_MODE_INVALID;
    const uint32_t                input_origin_index =
        (context_ptr->blk_origin_y + input_picture_ptr->origin_y) * input_picture_ptr->stride_y +
        (context_ptr->blk_origin_x + input_picture_ptr->origin_x);

    const uint32_t blk_origin_index =
        blk_geom->origin_x + blk_geom->origin_y * context_ptr->sb_size;
    BlkStruct *blk_ptr         = context_ptr->blk_ptr;
    candidate_buffer_ptr_array = &(candidate_buffer_ptr_array_base[0]);
    first_pass_signal_derivation_block(context_ptr);

    blk_ptr->av1xd->tile.mi_col_start = context_ptr->sb_ptr->tile_info.mi_col_start;
    blk_ptr->av1xd->tile.mi_col_end   = context_ptr->sb_ptr->tile_info.mi_col_end;
    blk_ptr->av1xd->tile.mi_row_start = context_ptr->sb_ptr->tile_info.mi_row_start;
    blk_ptr->av1xd->tile.mi_row_end   = context_ptr->sb_ptr->tile_info.mi_row_end;

    product_coding_loop_init_fast_loop(context_ptr,
                                       context_ptr->skip_coeff_neighbor_array,
                                       context_ptr->inter_pred_dir_neighbor_array,
                                       context_ptr->ref_frame_type_neighbor_array,
                                       context_ptr->intra_luma_mode_neighbor_array,
                                       context_ptr->skip_flag_neighbor_array,
                                       context_ptr->mode_type_neighbor_array,
                                       context_ptr->leaf_depth_neighbor_array,
                                       context_ptr->leaf_partition_neighbor_array);

    FrameHeader *frm_hdr = &pcs_ptr->parent_pcs_ptr->frm_hdr;
    // Generate MVP(s)
    if (!context_ptr->shut_fast_rate) {
        if (frm_hdr->allow_intrabc)
            generate_av1_mvp_table(&context_ptr->sb_ptr->tile_info,
                                   context_ptr,
                                   context_ptr->blk_ptr,
                                   context_ptr->blk_geom,
                                   context_ptr->blk_origin_x,
                                   context_ptr->blk_origin_y,
                                   pcs_ptr->parent_pcs_ptr->ref_frame_type_arr,
                                   1,
                                   pcs_ptr);
        else if (pcs_ptr->slice_type != I_SLICE)
            generate_av1_mvp_table(&context_ptr->sb_ptr->tile_info,
                                   context_ptr,
                                   context_ptr->blk_ptr,
                                   context_ptr->blk_geom,
                                   context_ptr->blk_origin_x,
                                   context_ptr->blk_origin_y,
                                   pcs_ptr->parent_pcs_ptr->ref_frame_type_arr,
                                   1,
                                   pcs_ptr);
    } else {
        init_xd(pcs_ptr, context_ptr);
    }
    // Read and (if needed) perform 1/8 Pel ME MVs refinement
    if (pcs_ptr->slice_type != I_SLICE)
        read_refine_me_mvs(pcs_ptr, context_ptr, input_picture_ptr);
    if (context_ptr->ref_pruning_ctrls.enabled)
        // Perform md reference pruning
        perform_md_reference_pruning(pcs_ptr, context_ptr, input_picture_ptr);
    context_ptr->inject_inter_candidates = 1;
    first_pass_generate_md_stage_0_cand(
     context_ptr, &fast_candidate_total_count, pcs_ptr);

    int32_t        mb_row = context_ptr->blk_origin_y >> 4;
    int32_t        mb_col = context_ptr->blk_origin_x >> 4;
    const uint32_t mb_cols =
        (pcs_ptr->parent_pcs_ptr->scs_ptr->seq_header.max_frame_width + 16 - 1) / 16;
    FRAME_STATS *mb_stats =
        pcs_ptr->parent_pcs_ptr->firstpass_data.mb_stats + mb_row * mb_cols + mb_col;
    //int *raw_motion_err_list = pcs_ptr->parent_pcs_ptr->firstpass_data.raw_motion_err_list +
    //    mb_row * mb_cols + mb_col;

    ModeDecisionCandidate *candidate_ptr;
    uint32_t               cand_index = 0;
    candidate_buffer                  = candidate_buffer_ptr_array[cand_index];
    candidate_ptr = candidate_buffer->candidate_ptr = &fast_candidate_array[cand_index];

    int this_intra_error = firstpass_intra_prediction(pcs_ptr,
                                                      context_ptr,
                                                      candidate_buffer,
                                                      candidate_ptr,
                                                      input_picture_ptr,
                                                      input_origin_index,
                                                      blk_origin_index,
                                                      mb_stats);

    int this_inter_error = this_intra_error;
    if (pcs_ptr->slice_type != I_SLICE && fast_candidate_total_count > 1) {
        this_inter_error = firstpass_inter_prediction(pcs_ptr,
                                                      context_ptr,
                                                      candidate_buffer,
                                                      candidate_ptr,
                                                      input_picture_ptr,
                                                      input_origin_index,
                                                      blk_origin_index,
                                                      fast_candidate_total_count,
                                                      this_intra_error,
                                                      //raw_motion_err_list,
                                                      mb_stats);

        mb_stats->coded_error += this_inter_error;
    } else {
        mb_stats->sr_coded_error += this_intra_error;
        mb_stats->tr_coded_error += this_intra_error;
        mb_stats->coded_error += this_intra_error;
    }
    // choose between Intra and inter LAST based on inter/intra error
    if (this_inter_error < this_intra_error)
        context_ptr->best_candidate_index_array[0] = 1;
    else
        context_ptr->best_candidate_index_array[0] = 0;
    // Handle stat for non 16x16 blocks. For non 16x16 blocks, some of the stats are increased multiple times
    // First find the last block in the 16x16 area and then devide the stats by the number of small blocks
    if (context_ptr->blk_geom->bsize != BLOCK_16X16 &&
        (context_ptr->blk_origin_x + context_ptr->blk_geom->bwidth ==
             pcs_ptr->parent_pcs_ptr->aligned_width ||
         (context_ptr->blk_geom->origin_x + context_ptr->blk_geom->bwidth) % FORCED_BLK_SIZE ==
             0) &&
        (context_ptr->blk_origin_y + context_ptr->blk_geom->bheight ==
             pcs_ptr->parent_pcs_ptr->aligned_height ||
         (context_ptr->blk_geom->origin_y + context_ptr->blk_geom->bheight) % FORCED_BLK_SIZE ==
             0)) {
        int blk_num =
            (((context_ptr->blk_geom->origin_x % FORCED_BLK_SIZE) + context_ptr->blk_geom->bwidth) /
             context_ptr->blk_geom->bwidth) *
            (((context_ptr->blk_geom->origin_y % FORCED_BLK_SIZE) +
              context_ptr->blk_geom->bheight) /
             context_ptr->blk_geom->bheight);
        average_non_16x16_stats(mb_stats, blk_num);
    }

    // Full Mode Decision (choose the best mode)
    candidate_index = product_full_mode_decision(
        context_ptr,
        blk_ptr,
        candidate_buffer_ptr_array,
        1,
        context_ptr->best_candidate_index_array,
        &best_intra_mode);
    candidate_buffer = candidate_buffer_ptr_array[candidate_index];

    bestcandidate_buffers[0] = candidate_buffer;
    uint8_t sq_index         = svt_log2f(context_ptr->blk_geom->sq_size) - 2;
    if (context_ptr->blk_geom->shape == PART_N) {
        context_ptr->parent_sq_type[sq_index] = candidate_buffer->candidate_ptr->type;

        context_ptr->parent_sq_has_coeff[sq_index] =
            (candidate_buffer->candidate_ptr->y_has_coeff ||
             candidate_buffer->candidate_ptr->u_has_coeff ||
             candidate_buffer->candidate_ptr->v_has_coeff)
                ? 1
                : 0;

        context_ptr->parent_sq_pred_mode[sq_index] = candidate_buffer->candidate_ptr->pred_mode;
    }
        if (!context_ptr->blk_geom->has_uv) {
        // Store the luma data for 4x* and *x4 blocks to be used for CFL
        EbPictureBufferDesc *recon_ptr = candidate_buffer->recon_ptr;
        uint32_t             rec_luma_offset =
            context_ptr->blk_geom->origin_x + context_ptr->blk_geom->origin_y * recon_ptr->stride_y;
        if (context_ptr->hbd_mode_decision) {
            for (uint32_t j = 0; j < context_ptr->blk_geom->bheight; ++j)
                memcpy(
                    context_ptr->cfl_temp_luma_recon16bit + rec_luma_offset +
                        j * recon_ptr->stride_y,
                    ((uint16_t *)recon_ptr->buffer_y) + (rec_luma_offset + j * recon_ptr->stride_y),
                    sizeof(uint16_t) * context_ptr->blk_geom->bwidth);
        } else {
            for (uint32_t j = 0; j < context_ptr->blk_geom->bheight; ++j)
                memcpy(&context_ptr->cfl_temp_luma_recon[rec_luma_offset + j * recon_ptr->stride_y],
                       recon_ptr->buffer_y + rec_luma_offset + j * recon_ptr->stride_y,
                       context_ptr->blk_geom->bwidth);
        }
    }
#if FIX_10BIT
        if (pcs_ptr->parent_pcs_ptr->scs_ptr->encoder_bit_depth == EB_8BIT) {
#endif
            if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag) {
                EbPictureBufferDesc *ref_pic = ((EbReferenceObject *)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->reference_picture;
                uint8_t *src_ptr = input_picture_ptr->buffer_y + input_origin_index;//recon_ptr->buffer_y + rec_luma_offset;
                uint8_t *dst_ptr =
                    ref_pic->buffer_y + ref_pic->origin_x + context_ptr->blk_origin_x +
                    (ref_pic->origin_y + context_ptr->blk_origin_y) * ref_pic->stride_y;
                for (uint32_t j = 0; j < context_ptr->blk_geom->bheight; j++)
                    svt_memcpy(dst_ptr + j * ref_pic->stride_y,
                        src_ptr + j * input_picture_ptr->stride_y,
                        context_ptr->blk_geom->bwidth * sizeof(uint8_t));
            }
#if FIX_10BIT
        }
        else {
            EbPictureBufferDesc *input_picture = pcs_ptr->input_frame16bit;

            EbPictureBufferDesc *ref_pic = ((EbReferenceObject *)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->reference_picture16bit;

             uint32_t src_block_index = input_picture->origin_x + context_ptr->blk_origin_x +
                    (input_picture->origin_y + context_ptr->blk_origin_y) * input_picture->stride_y;

             uint32_t dst_block_index = ref_pic->origin_x + context_ptr->blk_origin_x +
                    (ref_pic->origin_y + context_ptr->blk_origin_y) * ref_pic->stride_y;

            uint16_t *src = &(((uint16_t *)input_picture->buffer_y)[src_block_index]);
            uint16_t *dst = &(((uint16_t *)ref_pic->buffer_y)[dst_block_index]);

            for (int i = 0; i < context_ptr->blk_geom->bheight; i++) {
                svt_memcpy(dst, src, context_ptr->blk_geom->bwidth * sizeof(uint16_t));
                src += input_picture->stride_y;
                dst += ref_pic->stride_y;
            }
        }
#endif
#if !FIX_VBR_BUG
    //copy neigh recon data in blk_ptr
    {
        uint32_t             j;
        EbPictureBufferDesc *recon_ptr = candidate_buffer->recon_ptr;
        uint32_t             rec_luma_offset =
            context_ptr->blk_geom->origin_x + context_ptr->blk_geom->origin_y * recon_ptr->stride_y;

        uint32_t rec_cb_offset = ((((context_ptr->blk_geom->origin_x >> 3) << 3) +
                                   ((context_ptr->blk_geom->origin_y >> 3) << 3) *
                                       candidate_buffer->recon_ptr->stride_cb) >>
                                  1);
        uint32_t rec_cr_offset = ((((context_ptr->blk_geom->origin_x >> 3) << 3) +
                                   ((context_ptr->blk_geom->origin_y >> 3) << 3) *
                                       candidate_buffer->recon_ptr->stride_cr) >>
                                  1);

        if (!context_ptr->hbd_mode_decision) {
            memcpy(context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                       .neigh_top_recon[0],
                   recon_ptr->buffer_y + rec_luma_offset +
                       (context_ptr->blk_geom->bheight - 1) * recon_ptr->stride_y,
                   context_ptr->blk_geom->bwidth);
            if (context_ptr->blk_geom->has_uv && context_ptr->chroma_level <= CHROMA_MODE_1) {
                memcpy(context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                           .neigh_top_recon[1],
                       recon_ptr->buffer_cb + rec_cb_offset +
                           (context_ptr->blk_geom->bheight_uv - 1) * recon_ptr->stride_cb,
                       context_ptr->blk_geom->bwidth_uv);
                memcpy(context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                           .neigh_top_recon[2],
                       recon_ptr->buffer_cr + rec_cr_offset +
                           (context_ptr->blk_geom->bheight_uv - 1) * recon_ptr->stride_cr,
                       context_ptr->blk_geom->bwidth_uv);
            }

            for (j = 0; j < context_ptr->blk_geom->bheight; ++j)
                context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                    .neigh_left_recon[0][j] =
                    recon_ptr->buffer_y[rec_luma_offset + context_ptr->blk_geom->bwidth - 1 +
                                        j * recon_ptr->stride_y];

            if (context_ptr->blk_geom->has_uv && context_ptr->chroma_level <= CHROMA_MODE_1) {
                for (j = 0; j < context_ptr->blk_geom->bheight_uv; ++j) {
                    context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                        .neigh_left_recon[1][j] =
                        recon_ptr->buffer_cb[rec_cb_offset + context_ptr->blk_geom->bwidth_uv - 1 +
                                             j * recon_ptr->stride_cb];
                    context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                        .neigh_left_recon[2][j] =
                        recon_ptr->buffer_cr[rec_cr_offset + context_ptr->blk_geom->bwidth_uv - 1 +
                                             j * recon_ptr->stride_cr];
                }
            }
            if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag) {
                EbPictureBufferDesc *ref_pic = ((EbReferenceObject *)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->reference_picture;
                uint8_t *src_ptr = recon_ptr->buffer_y + rec_luma_offset;
                uint8_t *dst_ptr =
                    ref_pic->buffer_y + ref_pic->origin_x + context_ptr->blk_origin_x +
                    (ref_pic->origin_y + context_ptr->blk_origin_y) * ref_pic->stride_y;
                for (j = 0; j < context_ptr->blk_geom->bheight; j++)
                    svt_memcpy(dst_ptr + j * ref_pic->stride_y,
                        src_ptr + j * recon_ptr->stride_y,
                        context_ptr->blk_geom->bwidth * sizeof(uint8_t));
            }

        } else {
            uint16_t sz = sizeof(uint16_t);
            memcpy(
                context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                    .neigh_top_recon_16bit[0],
                recon_ptr->buffer_y + sz * (rec_luma_offset + (context_ptr->blk_geom->bheight - 1) *
                                                                  recon_ptr->stride_y),
                sz * context_ptr->blk_geom->bwidth);
            if (context_ptr->blk_geom->has_uv && context_ptr->chroma_level <= CHROMA_MODE_1) {
                memcpy(context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                           .neigh_top_recon_16bit[1],
                       recon_ptr->buffer_cb +
                           sz * (rec_cb_offset +
                                 (context_ptr->blk_geom->bheight_uv - 1) * recon_ptr->stride_cb),
                       sz * context_ptr->blk_geom->bwidth_uv);
                memcpy(context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                           .neigh_top_recon_16bit[2],
                       recon_ptr->buffer_cr +
                           sz * (rec_cr_offset +
                                 (context_ptr->blk_geom->bheight_uv - 1) * recon_ptr->stride_cr),
                       sz * context_ptr->blk_geom->bwidth_uv);
            }

            for (j = 0; j < context_ptr->blk_geom->bheight; ++j)
                context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                    .neigh_left_recon_16bit[0][j] =
                    ((uint16_t *)
                         recon_ptr->buffer_y)[rec_luma_offset + context_ptr->blk_geom->bwidth - 1 +
                                              j * recon_ptr->stride_y];

            if (context_ptr->blk_geom->has_uv && context_ptr->chroma_level <= CHROMA_MODE_1) {
                for (j = 0; j < context_ptr->blk_geom->bheight_uv; ++j) {
                    context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                        .neigh_left_recon_16bit[1][j] =
                        ((uint16_t *)recon_ptr
                             ->buffer_cb)[rec_cb_offset + context_ptr->blk_geom->bwidth_uv - 1 +
                                          j * recon_ptr->stride_cb];
                    context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                        .neigh_left_recon_16bit[2][j] =
                        ((uint16_t *)recon_ptr
                             ->buffer_cr)[rec_cr_offset + context_ptr->blk_geom->bwidth_uv - 1 +
                                          j * recon_ptr->stride_cr];
                }
            }
        }
    }
#endif
    context_ptr->md_local_blk_unit[blk_ptr->mds_idx].avail_blk_flag = EB_TRUE;
}

void set_tf_controls(PictureDecisionContext *context_ptr, uint8_t tf_level);
/******************************************************
* Derive Multi-Processes Settings for first pass
Input   : encoder mode and tune
Output  : Multi-Processes signal(s)
******************************************************/
EbErrorType first_pass_signal_derivation_multi_processes(SequenceControlSet *     scs_ptr,
                                                         PictureParentControlSet *pcs_ptr,
                                                         PictureDecisionContext * context_ptr) {
    EbErrorType return_error = EB_ErrorNone;
    FrameHeader *frm_hdr = &pcs_ptr->frm_hdr;
    // If enabled here, the hme enable flags should also be enabled in ResourceCoordinationProcess
    // to ensure that resources are allocated for the downsampled pictures used in HME
    pcs_ptr->enable_hme_flag = 1;
    pcs_ptr->enable_hme_level0_flag = 1;
    pcs_ptr->enable_hme_level1_flag = 1;
    pcs_ptr->enable_hme_level2_flag = 1;

    pcs_ptr->tf_enable_hme_flag = 0;
    pcs_ptr->tf_enable_hme_level0_flag = 0;
    pcs_ptr->tf_enable_hme_level1_flag = 0;
    pcs_ptr->tf_enable_hme_level2_flag = 0;

    // Set the Multi-Pass PD level
    pcs_ptr->multi_pass_pd_level = MULTI_PASS_PD_OFF;

    // Set disallow_nsq
    pcs_ptr->disallow_nsq = EB_TRUE;

    pcs_ptr->max_number_of_pus_per_sb = SQUARE_PU_COUNT;
    pcs_ptr->disallow_all_nsq_blocks_below_8x8 = EB_TRUE;

    // Set disallow_all_nsq_blocks_below_16x16: 16x8, 8x16, 16x4, 4x16
    pcs_ptr->disallow_all_nsq_blocks_below_16x16 = EB_TRUE;

    pcs_ptr->disallow_all_nsq_blocks_below_64x64 = EB_TRUE;
    pcs_ptr->disallow_all_nsq_blocks_below_32x32 = EB_TRUE;
    pcs_ptr->disallow_all_nsq_blocks_above_64x64 = EB_TRUE;
    pcs_ptr->disallow_all_nsq_blocks_above_32x32 = EB_TRUE;
    // disallow_all_nsq_blocks_above_16x16
    pcs_ptr->disallow_all_nsq_blocks_above_16x16 = EB_TRUE;

    pcs_ptr->disallow_HVA_HVB_HV4 = EB_TRUE;
    pcs_ptr->disallow_HV4 = EB_TRUE;

    // Set disallow_all_non_hv_nsq_blocks_below_16x16
    pcs_ptr->disallow_all_non_hv_nsq_blocks_below_16x16 = EB_TRUE;

    // Set disallow_all_h4_v4_blocks_below_16x16
    pcs_ptr->disallow_all_h4_v4_blocks_below_16x16 = EB_TRUE;

    frm_hdr->allow_screen_content_tools = 0;
    frm_hdr->allow_intrabc = 0;

    // Palette Modes:
    //    0:OFF
    //    1:Slow    NIC=7/4/4
    //    2:        NIC=7/2/2
    //    3:        NIC=7/2/2 + No K means for non ref
    //    4:        NIC=4/2/1
    //    5:        NIC=4/2/1 + No K means for Inter frame
    //    6:        Fastest NIC=4/2/1 + No K means for non base + step for non base for most dominent
    pcs_ptr->palette_level = 0;
    // Loop filter Level                            Settings
    // 0                                            OFF
    // 1                                            CU-BASED
    // 2                                            LIGHT FRAME-BASED
    // 3                                            FULL FRAME-BASED
    pcs_ptr->loop_filter_mode = 0;

    // CDEF Level                                   Settings
    // 0                                            OFF
    // 1                                            1 step refinement
    // 2                                            4 step refinement
    // 3                                            8 step refinement
    // 4                                            16 step refinement
    // 5                                            64 step refinement
    pcs_ptr->cdef_level = 0;

    // SG Level                                    Settings
    // 0                                            OFF
    // 1                                            0 step refinement
    // 2                                            1 step refinement
    // 3                                            4 step refinement
    // 4                                            16 step refinement
    Av1Common *cm = pcs_ptr->av1_cm;
    cm->sg_filter_mode = 0;

    // WN Level                                     Settings
    // 0                                            OFF
    // 1                                            3-Tap luma/ 3-Tap chroma
    // 2                                            5-Tap luma/ 5-Tap chroma
    // 3                                            7-Tap luma/ 5-Tap chroma
    cm->wn_filter_mode = 0;

    // Intra prediction modes                       Settings
    // 0                                            FULL
    // 1                                            LIGHT per block : disable_z2_prediction && disable_angle_refinement  for 64/32/4
    // 2                                            OFF per block : disable_angle_prediction for 64/32/4
    // 3                                            OFF : disable_angle_prediction
    // 4                                            OIS based Intra
    // 5                                            Light OIS based Intra
    pcs_ptr->intra_pred_mode = 3;

    // Set Tx Search     Settings
    // 0                 OFF
    // 1                 ON
    pcs_ptr->tx_size_search_mode = 1;

    // Set frame end cdf update mode      Settings
    // 0                                  OFF
    // 1                                  ON
    if (scs_ptr->static_config.frame_end_cdf_update == DEFAULT)
        pcs_ptr->frame_end_cdf_update_mode = 0;
    else
        pcs_ptr->frame_end_cdf_update_mode =
        scs_ptr->static_config.frame_end_cdf_update;

     pcs_ptr->frm_hdr.use_ref_frame_mvs = 0;

    // Global motion level                        Settings
    // GM_FULL                                    Exhaustive search mode.
    // GM_DOWN                                    Downsampled resolution with a
    // downsampling factor of 2 in each dimension GM_TRAN_ONLY Translation only
    // using ME MV.
    pcs_ptr->gm_level = GM_DOWN;

    // Exit TX size search when all coefficients are zero
    // 0: OFF
    // 1: ON
    pcs_ptr->tx_size_early_exit = 0;


    context_ptr->tf_level = 0;
    set_tf_controls(context_ptr, context_ptr->tf_level);
    // MRP control
    // 0: OFF (1,1)  ; override features
    // 1: FULL (4,3) ; override features
    // 2: (4,3) ; No-override features
    // 3: (3,3) ; No-override features
    // 4: (3,2) ; No-override features
    // 5: (2,3) ; No-override features
    // 6: (2,2) ; No-override features
    // 7: (2,1) ; No-override features
    // 8: (1,2) ; No-override features
    // 9: (1,1) ; No-override features
    // Level 0 , 1  : set ref_list0_count_try and ref_list1_count_try and Override MRP-related features
    // Level 2 .. 9 : Only set ref_list0_count_try and ref_list1_count_try

    pcs_ptr->tpl_opt_flag = 1;
    return return_error;
}
void set_txt_cycle_reduction_controls(ModeDecisionContext *mdctxt, uint8_t txt_cycles_red_mode);
void set_nsq_cycle_redcution_controls(ModeDecisionContext *mdctxt, uint16_t nsq_cycles_red_mode);
void set_depth_cycle_redcution_controls(ModeDecisionContext *mdctxt, uint8_t depth_cycles_red_mode) ;
void adaptive_md_cycles_redcution_controls(ModeDecisionContext *mdctxt, uint8_t adaptive_md_cycles_red_mode);
void set_obmc_controls(ModeDecisionContext *mdctxt, uint8_t obmc_mode) ;
void set_txs_cycle_reduction_controls(ModeDecisionContext *mdctxt, uint8_t txs_cycles_red_mode);

void coeff_based_switch_md_controls(ModeDecisionContext *mdctxt, uint8_t switch_md_mode_based_on_sq_coeff_level);
void md_subpel_me_controls(ModeDecisionContext *mdctxt, uint8_t md_subpel_me_level);
void md_subpel_pme_controls(ModeDecisionContext *mdctxt, uint8_t md_subpel_pme_level);
void md_nsq_motion_search_controls(ModeDecisionContext *mdctxt, uint8_t md_nsq_mv_search_level);
void md_pme_search_controls(ModeDecisionContext *mdctxt, uint8_t md_pme_level);

/******************************************************
* Derive EncDec Settings for first pass
Input   : encoder mode and pd pass
Output  : EncDec Kernel signal(s)
******************************************************/
EbErrorType first_pass_signal_derivation_enc_dec_kernel(
    PictureControlSet *pcs_ptr,
    ModeDecisionContext *context_ptr) {
    EbErrorType return_error = EB_ErrorNone;

    uint8_t pd_pass = context_ptr->pd_pass;
    // sb_classifier levels
    // Level                Settings
    // 0                    Off
    // 1                    TH 80%
    // 2                    TH 70%
    // 3                    TH 60%
    // 4                    TH 50%
    // 5                    TH 40%
    context_ptr->enable_area_based_cycles_allocation = 0;
    // Tx_search Level for Luma                       Settings
    // TX_SEARCH_DCT_DCT_ONLY                         DCT_DCT only
    // TX_SEARCH_DCT_TX_TYPES                         Tx search DCT type(s): DCT_DCT, V_DCT, H_DCT
    // TX_SEARCH_ALL_TX_TYPES                         Tx search all type(s)
    context_ptr->tx_search_level = TX_SEARCH_DCT_DCT_ONLY;
    uint8_t txt_cycles_reduction_level = 0;
    set_txt_cycle_reduction_controls(context_ptr, txt_cycles_reduction_level);
    context_ptr->interpolation_search_level = IFS_OFF;
    // Set Chroma Mode
    // Level                Settings
    // CHROMA_MODE_0  0     Full chroma search @ MD
    // CHROMA_MODE_1  1     Fast chroma search @ MD
    // CHROMA_MODE_2  2     Chroma blind @ MD + CFL @ EP
    // CHROMA_MODE_3  3     Chroma blind @ MD + no CFL @ EP
    context_ptr->chroma_level = CHROMA_MODE_2;

    // Chroma independent modes search
    // Level                Settings
    // 0                    post first md_stage
    // 1                    post last md_stage
    context_ptr->chroma_at_last_md_stage = 0;
    context_ptr->chroma_at_last_md_stage_intra_th = (uint64_t)~0;
    context_ptr->chroma_at_last_md_stage_cfl_th = (uint64_t)~0;

    // Cfl level
    // Level                Settings
    // 0                    Allow cfl
    // 1                    Disable cfl

    context_ptr->md_disable_cfl = EB_TRUE;

    // Set disallow_4x4
    context_ptr->disallow_4x4 = EB_FALSE;

    context_ptr->md_disallow_nsq = pcs_ptr->parent_pcs_ptr->disallow_nsq;

    // Set global MV injection
    // Level                Settings
    // 0                    Injection off
    // 1                    On
    context_ptr->global_mv_injection = 0;

    context_ptr->new_nearest_injection = 0;
    context_ptr->new_nearest_near_comb_injection = 0;

    // Set warped motion injection
    // Level                Settings
    // 0                    OFF
    // 1                    On
    context_ptr->warped_motion_injection = 0;

    // Set unipred3x3 injection
    // Level                Settings
    // 0                    OFF
    // 1                    ON FULL
    // 2                    Reduced set
    context_ptr->unipred3x3_injection = 0;

    // Set bipred3x3 injection
    // Level                Settings
    // 0                    OFF
    // 1                    ON FULL
    // 2                    Reduced set
    context_ptr->bipred3x3_injection = 0;

    // Level   Settings
    // 0       OFF: No compound mode search : AVG only
    // 1       ON: Full - AVG/DIST/DIFF/WEDGE
    // 2       ON: Fast - Use AVG only for non-closest ref frames or ref frames with high distortion
    context_ptr->inter_compound_mode = 0;
    // Set PME level
    context_ptr->md_pme_level = 0;
    md_pme_search_controls(context_ptr, context_ptr->md_pme_level);
    if (pd_pass == PD_PASS_0) {
        context_ptr->md_staging_mode = MD_STAGING_MODE_0;
    }
    else if (pd_pass == PD_PASS_1) {
        context_ptr->md_staging_mode = MD_STAGING_MODE_1;
    }
    else
        context_ptr->md_staging_mode = MD_STAGING_MODE_1;


    // Set md staging count level
    // Level 0              minimum count = 1
    // Level 1              set towards the best possible partitioning (to further optimize)
    // Level 2              HG: breack down or look up-table(s) are required !
    if (pd_pass == PD_PASS_0) {
        context_ptr->md_staging_count_level = 0;
    }
    else if (pd_pass == PD_PASS_1) {
        context_ptr->md_staging_count_level = 1;
    }
    else {
        context_ptr->md_staging_count_level = 2;
    }

    // Derive Spatial SSE Flag
    context_ptr->spatial_sse_full_loop_level = EB_TRUE;

    context_ptr->blk_skip_decision = EB_FALSE;

    // Derive Trellis Quant Coeff Optimization Flag
    context_ptr->rdoq_level = EB_FALSE;

    // Derive redundant block
    context_ptr->redundant_blk = EB_FALSE;

    // md_stage_1_cand_prune_th (for single candidate removal per class)
    // Remove candidate if deviation to the best is higher than md_stage_1_cand_prune_th
    context_ptr->md_stage_1_cand_prune_th = (uint64_t)~0;
    // md_stage_1_class_prune_th (for class removal)
    // Remove class if deviation to the best higher than TH_C

    context_ptr->md_stage_1_class_prune_th = (uint64_t)~0;

    // md_stage_2_3_cand_prune_th (for single candidate removal per class)
    // Remove candidate if deviation to the best is higher than
    // md_stage_2_3_cand_prune_th
    context_ptr->md_stage_2_3_cand_prune_th = (uint64_t)~0;

    // md_stage_2_3_class_prune_th (for class removal)
    // Remove class if deviation to the best is higher than
    // md_stage_2_3_class_prune_th

    context_ptr->md_stage_2_3_class_prune_th = (uint64_t)~0;

    context_ptr->coeff_area_based_bypass_nsq_th = 0;

    uint8_t adaptive_md_cycles_level = 0;
    adaptive_md_cycles_redcution_controls(context_ptr, adaptive_md_cycles_level);
    // Weighting (expressed as a percentage) applied to
    // square shape costs for determining if a and b
    // shapes should be skipped. Namely:
    // skip HA, HB, and H4 if h_cost > (weighted sq_cost)
    // skip VA, VB, and V4 if v_cost > (weighted sq_cost)
    context_ptr->sq_weight = (uint32_t)~0;
    // Set coeff_based_nsq_cand_reduction
    context_ptr->switch_md_mode_based_on_sq_coeff = 0;
    coeff_based_switch_md_controls(context_ptr, context_ptr->switch_md_mode_based_on_sq_coeff);

    // Set pic_obmc_level @ MD
    context_ptr->md_pic_obmc_level = 0;
    set_obmc_controls(context_ptr, context_ptr->md_pic_obmc_level);

    // Set enable_inter_intra @ MD
    context_ptr->md_inter_intra_level = 0;

    // Set enable_paeth @ MD
    context_ptr->md_enable_paeth = 0;

    // Set enable_smooth @ MD
    context_ptr->md_enable_smooth = 0;

    // Set md_tx_size_search_mode @ MD
    context_ptr->md_tx_size_search_mode = pcs_ptr->parent_pcs_ptr->tx_size_search_mode;

    // Assign whether to use TXS in inter classes (if TXS is ON)
    // 0 OFF - Use TXS for intra candidates only
    // 1 ON  - Use TXS for all candidates
    // 2 ON  - INTER TXS restricted to max 1 depth
    context_ptr->txs_in_inter_classes = 0;


    // Each NIC scaling level corresponds to a scaling factor, given by the below {x,y}
    // combinations, where x is the numerator, and y is the denominator.  e.g. {1,8} corresponds
    // to 1/8x scaling of the base NICs, which are set in set_md_stage_counts().
    //{10,8 },    // level0
    //{ 8,8 },    // level1
    //{ 7,8 },    // level2
    //{ 6,8 },    // level3
    //{ 5,8 },    // level4
    //{ 4,8 },    // level5
    //{ 3,8 },    // level6
    //{ 2,8 },    // level7
    //{ 3,16},    // level8
    //{ 1,8 },    // level9
    //{ 1,16}     // level10
    context_ptr->nic_scaling_level = 9;

    uint8_t txs_cycles_reduction_level = 0;
    set_txs_cycle_reduction_controls(context_ptr, txs_cycles_reduction_level);

    // Set md_filter_intra_mode @ MD
    // md_filter_intra_level specifies whether filter intra would be active
    // for a given prediction candidate in mode decision.
    // md_filter_intra_level | Settings
    // 0                      | OFF
    // 1                      | ON
    context_ptr->md_filter_intra_level = 0;

    // Set md_allow_intrabc @ MD
    context_ptr->md_allow_intrabc = 0;

    // Set md_palette_level @ MD
    context_ptr->md_palette_level = 0;

    context_ptr->md_nsq_mv_search_level = 0;
    md_nsq_motion_search_controls(context_ptr, context_ptr->md_nsq_mv_search_level);

    context_ptr->md_nsq_mv_search_level = 4;
    md_nsq_motion_search_controls(context_ptr, context_ptr->md_nsq_mv_search_level);

    context_ptr->md_subpel_me_level = 0;
    md_subpel_me_controls(context_ptr, context_ptr->md_subpel_me_level);


    context_ptr->md_subpel_pme_level = 2;
    md_subpel_pme_controls(context_ptr, context_ptr->md_subpel_pme_level);

    // Set max_ref_count @ MD
    context_ptr->md_max_ref_count = 4;
    // Set dc_cand_only_flag
    context_ptr->dc_cand_only_flag = EB_TRUE;

    // Set intra_angle_delta @ MD
    context_ptr->md_intra_angle_delta = 0;

    // Set disable_angle_z2_prediction_flag
    context_ptr->disable_angle_z2_intra_flag = EB_TRUE;
    // Use coeff rate and slit flag rate only (i.e. no fast rate)
    context_ptr->shut_fast_rate = EB_FALSE;
    context_ptr->skip_intra = 0;

    context_ptr->mds3_intra_prune_th = (uint16_t)~0;
    context_ptr->skip_cfl_cost_dev_th = (uint16_t)~0;

    return return_error;
}
/******************************************************
* Derive Mode Decision Config Settings for first pass
Input   : encoder mode and tune
Output  : EncDec Kernel signal(s)
******************************************************/
EbErrorType first_pass_signal_derivation_mode_decision_config_kernel(
    PictureControlSet *pcs_ptr,
    ModeDecisionConfigurationContext *context_ptr) {

    EbErrorType return_error = EB_ErrorNone;

    // ADP
    context_ptr->adp_level = pcs_ptr->parent_pcs_ptr->enc_mode;

    // CDF
    pcs_ptr->update_cdf = 0;

    // Filter INTRA
    // pic_filter_intra_level specifies whether filter intra would be active
    // for a given picture.
    // pic_filter_intra_level | Settings
    // 0                      | OFF
    // 1                      | ON
    pcs_ptr->pic_filter_intra_level = 0;

    // High Precision
    FrameHeader *frm_hdr = &pcs_ptr->parent_pcs_ptr->frm_hdr;
    frm_hdr->allow_high_precision_mv = 0;

    // Warped
    frm_hdr->allow_warped_motion = 0;
    frm_hdr->is_motion_mode_switchable = frm_hdr->allow_warped_motion;

    // pic_obmc_level - pic_obmc_level is used to define md_pic_obmc_level.
    // The latter determines the OBMC settings in the function set_obmc_controls.
    // Please check the definitions of the flags/variables in the function
    // set_obmc_controls corresponding to the pic_obmc_level settings.
    //  pic_obmc_level  |              Default Encoder Settings             |     Command Line Settings
    //         0        | OFF subject to possible constraints               | OFF everywhere in encoder
    //         1        | ON subject to possible constraints                | Fully ON in PD_PASS_2
    //         2        | Faster level subject to possible constraints      | Level 2 everywhere in PD_PASS_2
    //         3        | Even faster level subject to possible constraints | Level 3 everywhere in PD_PASS_3
    pcs_ptr->parent_pcs_ptr->pic_obmc_level = 0;

    // Switchable Motion Mode
    frm_hdr->is_motion_mode_switchable = frm_hdr->is_motion_mode_switchable ||
        pcs_ptr->parent_pcs_ptr->pic_obmc_level;

    // HBD Mode
    pcs_ptr->hbd_mode_decision = EB_8_BIT_MD; //first pass hard coded to 8bit
    return return_error;
}
void* set_me_hme_params_oq(
    MeContext                     *me_context_ptr,
    PictureParentControlSet       *pcs_ptr,
    SequenceControlSet            *scs_ptr,
    EbInputResolution             input_resolution);
void *set_me_hme_params_from_config(SequenceControlSet *scs_ptr, MeContext *me_context_ptr) ;
void set_me_hme_ref_prune_ctrls(MeContext* context_ptr, uint8_t prune_level) ;
void set_me_sr_adjustment_ctrls(MeContext* context_ptr, uint8_t sr_adjustment_level);
/******************************************************
* Derive ME Settings for first pass
  Input   : encoder mode and tune
  Output  : ME Kernel signal(s)
******************************************************/
EbErrorType first_pass_signal_derivation_me_kernel(
    SequenceControlSet        *scs_ptr,
    PictureParentControlSet   *pcs_ptr,
    MotionEstimationContext_t   *context_ptr) {
    EbErrorType return_error = EB_ErrorNone;

    // Set ME/HME search regions

    if (scs_ptr->static_config.use_default_me_hme)
        set_me_hme_params_oq(
            context_ptr->me_context_ptr,
            pcs_ptr,
            scs_ptr,
            scs_ptr->input_resolution);
    else
        set_me_hme_params_from_config(
            scs_ptr,
            context_ptr->me_context_ptr);


    // Set HME flags
    context_ptr->me_context_ptr->enable_hme_flag = pcs_ptr->enable_hme_flag;
    context_ptr->me_context_ptr->enable_hme_level0_flag = pcs_ptr->enable_hme_level0_flag;
    context_ptr->me_context_ptr->enable_hme_level1_flag = pcs_ptr->enable_hme_level1_flag;
    context_ptr->me_context_ptr->enable_hme_level2_flag = pcs_ptr->enable_hme_level2_flag;

    // HME Search Method
    context_ptr->me_context_ptr->hme_search_method = SUB_SAD_SEARCH;

    // ME Search Method
    context_ptr->me_context_ptr->me_search_method = SUB_SAD_SEARCH;

    context_ptr->me_context_ptr->compute_global_motion = EB_FALSE;

    // Set hme/me based reference pruning level (0-4)
    set_me_hme_ref_prune_ctrls(context_ptr->me_context_ptr, 0);

    // Set hme-based me sr adjustment level
    set_me_sr_adjustment_ctrls(context_ptr->me_context_ptr, 0);

    return return_error;
};
