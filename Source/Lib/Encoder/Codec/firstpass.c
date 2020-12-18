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
#ifdef ARCH_X86_64
#include <xmmintrin.h>
#endif
#include "EbMotionEstimation.h"
//#include "EbMotionEstimationProcess.h"
#undef _MM_HINT_T2
#define _MM_HINT_T2 1

#define OUTPUT_FPF 0

#define INTRA_MODE_PENALTY 1024
#define NEW_MV_MODE_PENALTY 32
#define DARK_THRESH 64

#define NCOUNT_INTRA_THRESH 8192
#define NCOUNT_INTRA_FACTOR 3

#define STATS_CAPABILITY_INIT 100
//1.5 times larger than request.
#define STATS_CAPABILITY_GROW(s) (s * 3 / 2)
static EbErrorType realloc_stats_out(SequenceControlSet *scs_ptr, FirstPassStatsOut *out,
                                     uint64_t frame_number) {
    if (frame_number < out->size)
        return EB_ErrorNone;

    if ((int64_t)frame_number >= (int64_t)out->capability - 1) {
        size_t capability = (int64_t)frame_number >= (int64_t)STATS_CAPABILITY_INIT - 1
            ? STATS_CAPABILITY_GROW(frame_number)
            : STATS_CAPABILITY_INIT;
        if (scs_ptr->lap_enabled) {
            //store the data points before re-allocation
            uint64_t stats_in_start_offset = 0;
            uint64_t stats_in_offset       = 0;
            uint64_t stats_in_end_offset   = 0;
            if (frame_number) {
                stats_in_start_offset = scs_ptr->twopass.stats_buf_ctx->stats_in_start - out->stat;
                stats_in_offset       = scs_ptr->twopass.stats_in - out->stat;
                stats_in_end_offset   = scs_ptr->twopass.stats_buf_ctx->stats_in_end - out->stat;
            }
            EB_REALLOC_ARRAY(out->stat, capability);
            // restore the pointers after re-allocation is done
            scs_ptr->twopass.stats_buf_ctx->stats_in_start = out->stat + stats_in_start_offset;
            scs_ptr->twopass.stats_in                      = out->stat + stats_in_offset;
            scs_ptr->twopass.stats_buf_ctx->stats_in_end   = out->stat + stats_in_end_offset;
        } else {
            EB_REALLOC_ARRAY(out->stat, capability);
        }
        out->capability = capability;
    }
    out->size = frame_number + 1;
    return EB_ErrorNone;
}

static AOM_INLINE void output_stats(SequenceControlSet *scs_ptr, FIRSTPASS_STATS *stats,
                                    uint64_t frame_number) {
    FirstPassStatsOut *stats_out = &scs_ptr->encode_context_ptr->stats_out;
    svt_block_on_mutex(scs_ptr->encode_context_ptr->stat_file_mutex);
    if (realloc_stats_out(scs_ptr, stats_out, frame_number) != EB_ErrorNone) {
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
    section->raw_error_stdev          = 0.0;
    section->pcnt_third_ref           = 0.0;
    section->tr_coded_error           = 0.0;
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
    if (raw_motion_err_counts == 0)
        return 0;

    int i;
    for (i = 0; i < raw_motion_err_counts; i++) { sum_raw_err += raw_motion_err_list[i]; }
    raw_err_avg = (double)sum_raw_err / raw_motion_err_counts;
    for (i = 0; i < raw_motion_err_counts; i++) {
        raw_err_stdev += (raw_motion_err_list[i] - raw_err_avg) *
            (raw_motion_err_list[i] - raw_err_avg);
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
    if (is_zero_mv(&best_mv))
        return;

    ++stats->mv_count;
    // Non-zero vector, was it different from the last non zero vector?
    if (!is_equal_mv(&best_mv, last_mv))
        ++stats->new_mv_count;
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
static void update_firstpass_stats(PictureParentControlSet *pcs_ptr, const FRAME_STATS *const stats,
                                   const double raw_err_stdev, const int frame_number,
                                   const int64_t ts_duration) {
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
    const double raw_err_stdev              = raw_motion_error_stdev(raw_motion_err_list,
                                                        total_raw_motion_err_count);
    // Clamp the image start to rows/2. This number of rows is discarded top
    // and bottom as dead data so rows / 2 means the frame is blank.
    if ((stats.image_data_start_row > (int)mb_rows / 2) ||
        (stats.image_data_start_row == INVALID_ROW)) {
        stats.image_data_start_row = mb_rows / 2;
    }
    // Exclude any image dead zone
    if (stats.image_data_start_row > 0) {
        stats.intra_skip_count = AOMMAX(
            0, stats.intra_skip_count - (stats.image_data_start_row * (int)mb_cols * 2));
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
* Derive Pre-Analysis settings for first pass for pcs
Input   : encoder mode and tune
Output  : Pre-Analysis signal(s)
******************************************************/
extern EbErrorType first_pass_signal_derivation_pre_analysis_pcs(PictureParentControlSet *pcs_ptr) {
    EbErrorType return_error = EB_ErrorNone;
    // Derive HME Flag
    pcs_ptr->enable_hme_flag        = 1;
    pcs_ptr->enable_hme_level0_flag = 1;
    pcs_ptr->enable_hme_level1_flag = 1;
    pcs_ptr->enable_hme_level2_flag = 1;

    //// Set here to allocate resources for the downsampled pictures used in HME (generated in PictureAnalysis)
    //// Will be later updated for SC/NSC in PictureDecisionProcess
    pcs_ptr->tf_enable_hme_flag                    = 0;
    pcs_ptr->tf_enable_hme_level0_flag             = 0;
    pcs_ptr->tf_enable_hme_level1_flag             = 0;
    pcs_ptr->tf_enable_hme_level2_flag             = 0;

    return return_error;
}

/******************************************************
* Derive Pre-Analysis settings for first pass for scs
Input   : encoder mode and tune
Output  : Pre-Analysis signal(s)
******************************************************/
extern EbErrorType first_pass_signal_derivation_pre_analysis_scs(SequenceControlSet * scs_ptr) {
    EbErrorType return_error = EB_ErrorNone;
    scs_ptr->seq_header.enable_intra_edge_filter   = 0;
    scs_ptr->seq_header.pic_based_rate_est         = 0;
    scs_ptr->seq_header.enable_restoration         = 0;
    scs_ptr->seq_header.cdef_level /*enable_cdef*/ = 0;
    scs_ptr->seq_header.enable_warped_motion       = 0;

    return return_error;
}

#define LOW_MOTION_ERROR_THRESH 25
void set_tf_controls(PictureParentControlSet *pcs_ptr, uint8_t tf_level);
/******************************************************
* Derive Multi-Processes Settings for first pass
Input   : encoder mode and tune
Output  : Multi-Processes signal(s)
******************************************************/
EbErrorType first_pass_signal_derivation_multi_processes(SequenceControlSet *     scs_ptr,
                                                         PictureParentControlSet *pcs_ptr,
                                                         PictureDecisionContext * context_ptr) {
    EbErrorType  return_error = EB_ErrorNone;
    FrameHeader *frm_hdr      = &pcs_ptr->frm_hdr;
    // If enabled here, the hme enable flags should also be enabled in ResourceCoordinationProcess
    // to ensure that resources are allocated for the downsampled pictures used in HME
    pcs_ptr->enable_hme_flag        = 1;
    pcs_ptr->enable_hme_level0_flag = 1;
    pcs_ptr->enable_hme_level1_flag = 1;
    pcs_ptr->enable_hme_level2_flag = 1;

    pcs_ptr->tf_enable_hme_flag        = 0;
    pcs_ptr->tf_enable_hme_level0_flag = 0;
    pcs_ptr->tf_enable_hme_level1_flag = 0;
    pcs_ptr->tf_enable_hme_level2_flag = 0;

    // Set the Multi-Pass PD level
    pcs_ptr->multi_pass_pd_level = MULTI_PASS_PD_OFF;

    // Set disallow_nsq
    pcs_ptr->disallow_nsq = EB_TRUE;

    pcs_ptr->max_number_of_pus_per_sb          = SQUARE_PU_COUNT;
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
    pcs_ptr->disallow_HV4         = EB_TRUE;

    // Set disallow_all_non_hv_nsq_blocks_below_16x16
    pcs_ptr->disallow_all_non_hv_nsq_blocks_below_16x16 = EB_TRUE;

    // Set disallow_all_h4_v4_blocks_below_16x16
    pcs_ptr->disallow_all_h4_v4_blocks_below_16x16 = EB_TRUE;

    frm_hdr->allow_screen_content_tools = 0;
    frm_hdr->allow_intrabc              = 0;

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
    Av1Common *cm      = pcs_ptr->av1_cm;
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
        pcs_ptr->frame_end_cdf_update_mode = scs_ptr->static_config.frame_end_cdf_update;

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
    set_tf_controls(pcs_ptr, context_ptr->tf_level);
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
    pcs_ptr->tpl_trailing_frame_count = 0;
    return return_error;
}
/******************************************************
* Derive Mode Decision Config Settings for first pass
Input   : encoder mode and tune
Output  : EncDec Kernel signal(s)
******************************************************/
EbErrorType first_pass_signal_derivation_mode_decision_config_kernel(PictureControlSet *pcs_ptr) {
    EbErrorType return_error = EB_ErrorNone;
    // CDF
    pcs_ptr->cdf_ctrl.enabled = pcs_ptr->cdf_ctrl.update_coef = 0;
    pcs_ptr->cdf_ctrl.update_mv = pcs_ptr->cdf_ctrl.update_se = 0;

    // Filter INTRA
    // pic_filter_intra_level specifies whether filter intra would be active
    // for a given picture.
    // pic_filter_intra_level | Settings
    // 0                      | OFF
    // 1                      | ON
    pcs_ptr->pic_filter_intra_level = 0;

    // High Precision
    FrameHeader *frm_hdr             = &pcs_ptr->parent_pcs_ptr->frm_hdr;
    frm_hdr->allow_high_precision_mv = 0;

    // Warped
    frm_hdr->allow_warped_motion       = 0;
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
void *set_me_hme_params_oq(MeContext *me_context_ptr, PictureParentControlSet *pcs_ptr,
                           SequenceControlSet *scs_ptr, EbInputResolution input_resolution);
void *set_me_hme_params_from_config(SequenceControlSet *scs_ptr, MeContext *me_context_ptr);
void  set_me_hme_ref_prune_ctrls(MeContext *context_ptr, uint8_t prune_level);
void  set_me_sr_adjustment_ctrls(MeContext *context_ptr, uint8_t sr_adjustment_level);
void  set_gm_controls(PictureParentControlSet *pcs_ptr, uint8_t gm_level);
/******************************************************
* Derive ME Settings for first pass
  Input   : encoder mode and tune
  Output  : ME Kernel signal(s)
******************************************************/
EbErrorType first_pass_signal_derivation_me_kernel(SequenceControlSet *       scs_ptr,
                                                   PictureParentControlSet *  pcs_ptr,
                                                   MotionEstimationContext_t *context_ptr) {
    EbErrorType return_error = EB_ErrorNone;

    // Set ME/HME search regions

    if (scs_ptr->static_config.use_default_me_hme)
        set_me_hme_params_oq(
            context_ptr->me_context_ptr, pcs_ptr, scs_ptr, scs_ptr->input_resolution);
    else
        set_me_hme_params_from_config(scs_ptr, context_ptr->me_context_ptr);

    // Set HME flags
    context_ptr->me_context_ptr->enable_hme_flag        = pcs_ptr->enable_hme_flag;
    context_ptr->me_context_ptr->enable_hme_level0_flag = pcs_ptr->enable_hme_level0_flag;
    context_ptr->me_context_ptr->enable_hme_level1_flag = pcs_ptr->enable_hme_level1_flag;
    context_ptr->me_context_ptr->enable_hme_level2_flag = pcs_ptr->enable_hme_level2_flag;

    // HME Search Method
    context_ptr->me_context_ptr->hme_search_method = SUB_SAD_SEARCH;

    // ME Search Method
    context_ptr->me_context_ptr->me_search_method = SUB_SAD_SEARCH;
    uint8_t gm_level                              = 0;
    set_gm_controls(pcs_ptr, gm_level);

    // Set hme/me based reference pruning level (0-4)
    set_me_hme_ref_prune_ctrls(context_ptr->me_context_ptr, 0);

    // Set hme-based me sr adjustment level
    set_me_sr_adjustment_ctrls(context_ptr->me_context_ptr, 0);

    return return_error;
};

/***************************************************************************
* Computes and returns the intra pred error of a block using src.
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
static int open_loop_firstpass_intra_prediction(uint32_t blk_origin_x, uint32_t blk_origin_y,
                                                uint8_t bwidth, uint8_t bheight,
                                                EbPictureBufferDesc *input_picture_ptr,
                                                uint32_t             input_origin_index,
                                                FRAME_STATS *const   stats) {
    int32_t   mb_row      = blk_origin_y >> 4;
    int32_t   mb_col      = blk_origin_x >> 4;
    const int use_dc_pred = (mb_col || mb_row) && (!mb_col || !mb_row);

    uint32_t sub_blk_origin_x, sub_blk_origin_y;
    uint8_t *above_row;
    uint8_t *left_col;

    DECLARE_ALIGNED(16, uint8_t, left_data[MAX_TX_SIZE * 2 + 32]);
    DECLARE_ALIGNED(16, uint8_t, above_data[MAX_TX_SIZE * 2 + 32]);
    DECLARE_ALIGNED(32, uint8_t, predictor8[256 * 2]);
    uint8_t *predictor = predictor8;

    uint8_t sub_blk_rows = use_dc_pred
        ? (bheight == FORCED_BLK_SIZE && bwidth == FORCED_BLK_SIZE) ? 1 : bheight / 8
        : bheight / 4;
    uint8_t sub_blk_cols = use_dc_pred
        ? (bheight == FORCED_BLK_SIZE && bwidth == FORCED_BLK_SIZE) ? 1 : bwidth / 8
        : bwidth / 4;

    for (uint32_t sub_blk_index_y = 0; sub_blk_index_y < sub_blk_rows; ++sub_blk_index_y) {
        for (uint32_t sub_blk_index_x = 0; sub_blk_index_x < sub_blk_cols; ++sub_blk_index_x) {
            TxSize tx_size   = use_dc_pred
                  ? ((bheight == FORCED_BLK_SIZE && bwidth == FORCED_BLK_SIZE) ? TX_16X16 : TX_8X8)
                  : TX_4X4;
            sub_blk_origin_x = blk_origin_x + sub_blk_index_x * bwidth / sub_blk_cols;
            sub_blk_origin_y = blk_origin_y + sub_blk_index_y * bheight / sub_blk_rows;
            above_row        = above_data + 16;
            left_col         = left_data + 16;

            // Fill Neighbor Arrays
            update_neighbor_samples_array_open_loop_mb(above_row - 1,
                                                       left_col - 1,
                                                       input_picture_ptr,
                                                       input_picture_ptr->stride_y,
                                                       sub_blk_origin_x,
                                                       sub_blk_origin_y,
                                                       bwidth / sub_blk_cols,
                                                       bheight / sub_blk_rows);
            // PRED
            predictor = &predictor8[(sub_blk_origin_x - blk_origin_x) +
                                    (sub_blk_origin_y - blk_origin_y) * FORCED_BLK_SIZE];
            intra_prediction_open_loop_mb(0,
                                          DC_PRED,
                                          sub_blk_origin_x,
                                          sub_blk_origin_y,
                                          tx_size,
                                          above_row,
                                          left_col,
                                          predictor,
                                          FORCED_BLK_SIZE);
        }
    }

    EbSpatialFullDistType spatial_full_dist_type_fun = svt_spatial_full_distortion_kernel;
    int this_intra_error = (uint32_t)(spatial_full_dist_type_fun(input_picture_ptr->buffer_y,
                                                                 input_origin_index,
                                                                 input_picture_ptr->stride_y,
                                                                 predictor8,
                                                                 0,
                                                                 FORCED_BLK_SIZE,
                                                                 bwidth,
                                                                 bheight));

    if (this_intra_error < UL_INTRA_THRESH) {
        ++stats->intra_skip_count;
    } else if ((mb_col > 0) && (stats->image_data_start_row == INVALID_ROW)) {
        stats->image_data_start_row = mb_row;
    }
    // aom_clear_system_state();
    double log_intra = log1p((double)this_intra_error);
    if (log_intra < 10.0)
        stats->intra_factor += 1.0 + ((10.0 - log_intra) * 0.05);
    else
        stats->intra_factor += 1.0;

    int level_sample = input_picture_ptr->buffer_y[input_origin_index];

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

    const int stride = input_picture_ptr->stride_y;
    uint8_t * buf    = &(input_picture_ptr->buffer_y)[input_origin_index];
    for (int r8 = 0; r8 < 2; ++r8) {
        for (int c8 = 0; c8 < 2; ++c8) {
            stats->frame_avg_wavelet_energy += svt_av1_haar_ac_sad_8x8_uint8_input(
                buf + c8 * 8 + r8 * 8 * stride, stride, 0);
        }
    }

    // Accumulate the intra error.
    stats->intra_error += (int64_t)this_intra_error;
    return this_intra_error;
}
/***************************************************************************
* Computes and returns the inter prediction error from the src last frame.
* Computes inter prediction errors from the golden and updates stats accordingly.
* Modifies:
*    stats: many member params in it.
*  Returns:
*    this_inter_error
***************************************************************************/
static int open_loop_firstpass_inter_prediction(
    PictureParentControlSet *ppcs_ptr, uint32_t me_sb_addr, uint32_t blk_origin_x,
    uint32_t blk_origin_y, uint8_t bwidth, uint8_t bheight, EbPictureBufferDesc *input_picture_ptr,
    uint32_t input_origin_index, const int this_intra_error, MV *last_mv, int raw_motion_err,
    FRAME_STATS *stats) {
    int32_t        mb_row  = blk_origin_y >> 4;
    int32_t        mb_col  = blk_origin_x >> 4;
    const uint32_t mb_cols = (ppcs_ptr->scs_ptr->seq_header.max_frame_width + FORCED_BLK_SIZE - 1) /
        FORCED_BLK_SIZE;
    const uint32_t mb_rows = (ppcs_ptr->scs_ptr->seq_header.max_frame_height + FORCED_BLK_SIZE -
                              1) /
        FORCED_BLK_SIZE;
    int                   this_inter_error           = this_intra_error;
    FULLPEL_MV            mv                         = kZeroFullMv;
    EbSpatialFullDistType spatial_full_dist_type_fun = svt_spatial_full_distortion_kernel;

    int motion_error = 0;
    // TODO(pengchong): Replace the hard-coded threshold
    if (raw_motion_err > LOW_MOTION_ERROR_THRESH) {
        uint32_t           me_mb_offset = 0;
        BlockGeom          blk_geom;
        const MeSbResults *me_results = ppcs_ptr->pa_me_data->me_results[me_sb_addr];
        uint32_t           me_sb_size = ppcs_ptr->scs_ptr->sb_sz;
        blk_geom.origin_x             = blk_origin_x - (blk_origin_x / me_sb_size) * me_sb_size;
        blk_geom.origin_y             = blk_origin_y - (blk_origin_y / me_sb_size) * me_sb_size;
        blk_geom.bwidth               = FORCED_BLK_SIZE;
        blk_geom.bheight              = FORCED_BLK_SIZE;
        me_mb_offset       = get_me_info_index(ppcs_ptr->max_number_of_pus_per_sb, &blk_geom, 0, 0);
        uint8_t list_index = 0;
        uint8_t ref_pic_index = 0;
        mv.col =
            me_results
                ->me_mv_array[me_mb_offset * MAX_PA_ME_MV + (list_index ? 4 : 0) + ref_pic_index]
                .x_mv >>
            2;
        mv.row =
            me_results
                ->me_mv_array[me_mb_offset * MAX_PA_ME_MV + (list_index ? 4 : 0) + ref_pic_index]
                .y_mv >>
            2;

        EbPictureBufferDesc *last_input_picture_ptr = ppcs_ptr->first_pass_ref_count
            ? ppcs_ptr->first_pass_ref_ppcs_ptr[0]->enhanced_picture_ptr
            : NULL;
        int32_t ref_origin_index = last_input_picture_ptr->origin_x + (blk_origin_x + mv.col) +
            (blk_origin_y + mv.row + last_input_picture_ptr->origin_y) *
                last_input_picture_ptr->stride_y;

        motion_error = (uint32_t)(spatial_full_dist_type_fun(input_picture_ptr->buffer_y,
                                                             input_origin_index,
                                                             input_picture_ptr->stride_y,
                                                             last_input_picture_ptr->buffer_y,
                                                             ref_origin_index,
                                                             last_input_picture_ptr->stride_y,
                                                             bwidth,
                                                             bheight));

        // Assume 0,0 motion with no mv overhead.
        if (mv.col != 0 && mv.row != 0) {
            motion_error += NEW_MV_MODE_PENALTY;
        }
        // Motion search in 2nd reference frame.
        int gf_motion_error = motion_error;
        if (ppcs_ptr->first_pass_ref_count > 1 &&
            me_results->total_me_candidate_index[me_mb_offset] > 1) {
            // To convert full-pel MV
            list_index    = 0;
            ref_pic_index = 1;
            FULLPEL_MV gf_mv;
            gf_mv.col = me_results
                            ->me_mv_array[me_mb_offset * MAX_PA_ME_MV + (list_index ? 4 : 0) +
                                          ref_pic_index]
                            .x_mv >>
                2;
            gf_mv.row = me_results
                            ->me_mv_array[me_mb_offset * MAX_PA_ME_MV + (list_index ? 4 : 0) +
                                          ref_pic_index]
                            .y_mv >>
                2;
            EbPictureBufferDesc *golden_input_picture_ptr =
                ppcs_ptr->first_pass_ref_ppcs_ptr[1]->enhanced_picture_ptr;
            ref_origin_index = golden_input_picture_ptr->origin_x + (blk_origin_x + gf_mv.col) +
                (blk_origin_y + gf_mv.row + golden_input_picture_ptr->origin_y) *
                    golden_input_picture_ptr->stride_y;

            gf_motion_error = (uint32_t)(
                spatial_full_dist_type_fun(input_picture_ptr->buffer_y,
                                           input_origin_index,
                                           input_picture_ptr->stride_y,
                                           golden_input_picture_ptr->buffer_y,
                                           ref_origin_index,
                                           golden_input_picture_ptr->stride_y,
                                           bwidth,
                                           bheight));

            // Assume 0,0 motion with no mv overhead.
            if (gf_mv.col != 0 && gf_mv.row != 0) {
                gf_motion_error += NEW_MV_MODE_PENALTY;
            }
        }

        if (gf_motion_error < motion_error && gf_motion_error < this_intra_error) {
            ++stats->second_ref_count;
            motion_error = gf_motion_error;
        }
        // In accumulating a score for the 2nd reference frame take the
        // best of the motion predicted score and the intra coded error
        // (just as will be done for) accumulation of "coded_error" for
        // the last frame.
        if (ppcs_ptr->first_pass_ref_count > 1 && (gf_motion_error < motion_error * 3)) {
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
        stats->tr_coded_error += motion_error;
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
            stats->neutral_count += (double)motion_error /
                DOUBLE_DIVIDE_CHECK((double)this_intra_error);
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
        accumulate_mv_stats(best_mv, mv, mb_row, mb_col, mb_rows, mb_cols, last_mv, stats);
    }

    return this_inter_error;
}
/***************************************************************************
* Perform the processing for first pass.
* For each 16x16 blocks performs DC and ME results from LAST frame and store
* the required data.
***************************************************************************/
static EbErrorType first_pass_frame(PictureParentControlSet *ppcs_ptr) {
    EbPictureBufferDesc *input_picture_ptr      = ppcs_ptr->enhanced_picture_ptr;
    EbPictureBufferDesc *last_input_picture_ptr = ppcs_ptr->first_pass_ref_count
        ? ppcs_ptr->first_pass_ref_ppcs_ptr[0]->enhanced_picture_ptr
        : NULL;

    const uint32_t blk_cols = (uint32_t)(input_picture_ptr->width + FORCED_BLK_SIZE - 1) /
        FORCED_BLK_SIZE;
    const uint32_t blk_rows = (uint32_t)(input_picture_ptr->height + FORCED_BLK_SIZE - 1) /
        FORCED_BLK_SIZE;

    uint32_t me_sb_size         = ppcs_ptr->scs_ptr->sb_sz;
    uint32_t me_pic_width_in_sb = (ppcs_ptr->aligned_width + me_sb_size - 1) / me_sb_size;
    uint32_t me_sb_x, me_sb_y, me_sb_addr;

    uint32_t blk_width, blk_height, blk_origin_x, blk_origin_y;
    MV       first_top_mv = kZeroMv;
    MV       last_mv;
    uint32_t input_origin_index;

    EbSpatialFullDistType spatial_full_dist_type_fun = svt_spatial_full_distortion_kernel;

    for (uint32_t blk_index_y = 0; blk_index_y < blk_rows; ++blk_index_y) {
        for (uint32_t blk_index_x = 0; blk_index_x < blk_cols; ++blk_index_x) {
            blk_origin_x = blk_index_x * FORCED_BLK_SIZE;
            blk_origin_y = blk_index_y * FORCED_BLK_SIZE;
            me_sb_x      = blk_origin_x / me_sb_size;
            me_sb_y      = blk_origin_y / me_sb_size;
            me_sb_addr   = me_sb_x + me_sb_y * me_pic_width_in_sb;

            blk_width  = (ppcs_ptr->aligned_width - blk_origin_x) < FORCED_BLK_SIZE
                 ? ppcs_ptr->aligned_width - blk_origin_x
                 : FORCED_BLK_SIZE;
            blk_height = (ppcs_ptr->aligned_height - blk_origin_y) < FORCED_BLK_SIZE
                ? ppcs_ptr->aligned_height - blk_origin_y
                : FORCED_BLK_SIZE;

            input_origin_index = (input_picture_ptr->origin_y + blk_origin_y) *
                    input_picture_ptr->stride_y +
                (input_picture_ptr->origin_x + blk_origin_x);

            FRAME_STATS *mb_stats = ppcs_ptr->firstpass_data.mb_stats + blk_index_y * blk_cols +
                blk_index_x;

            int this_intra_error = open_loop_firstpass_intra_prediction(blk_origin_x,
                                                                        blk_origin_y,
                                                                        blk_width,
                                                                        blk_height,
                                                                        input_picture_ptr,
                                                                        input_origin_index,
                                                                        mb_stats);
            int this_inter_error = this_intra_error;

            if (blk_origin_x == 0)
                last_mv = first_top_mv;

            if (ppcs_ptr->first_pass_ref_count) {
                ppcs_ptr->firstpass_data.raw_motion_err_list[blk_index_y * blk_cols + blk_index_x] =
                    (uint32_t)(spatial_full_dist_type_fun(input_picture_ptr->buffer_y,
                                                          input_origin_index,
                                                          input_picture_ptr->stride_y,
                                                          last_input_picture_ptr->buffer_y,
                                                          input_origin_index,
                                                          input_picture_ptr->stride_y,
                                                          blk_width,
                                                          blk_height));

                this_inter_error = open_loop_firstpass_inter_prediction(
                    ppcs_ptr,
                    me_sb_addr,
                    blk_origin_x,
                    blk_origin_y,
                    blk_width,
                    blk_height,
                    input_picture_ptr,
                    input_origin_index,
                    this_intra_error,
                    &last_mv,
                    ppcs_ptr->firstpass_data
                        .raw_motion_err_list[blk_index_y * blk_cols + blk_index_x],
                    mb_stats);

                if (blk_origin_x == 0)
                    first_top_mv = last_mv;

                mb_stats->coded_error += this_inter_error;
            } else {
                mb_stats->sr_coded_error += this_intra_error;
                mb_stats->tr_coded_error += this_intra_error;
                mb_stats->coded_error += this_intra_error;
            }
        }
    }

    return EB_ErrorNone;
}
/***************************************************************************
* Prepare the me context for performing first pass me.
***************************************************************************/
static void first_pass_setup_me_context(MotionEstimationContext_t *context_ptr,
                                        PictureParentControlSet *  ppcs_ptr,
                                        EbPictureBufferDesc *input_picture_ptr, int blk_row,
                                        int blk_col, uint32_t ss_x, uint32_t ss_y) {
    // setup the references
    context_ptr->me_context_ptr->num_of_list_to_search       = 0;
    context_ptr->me_context_ptr->num_of_ref_pic_to_search[0] = 0;
    context_ptr->me_context_ptr->num_of_ref_pic_to_search[1] = 0;
    context_ptr->me_context_ptr->temporal_layer_index        = 0;
    context_ptr->me_context_ptr->is_used_as_reference_flag   = 1;

    if (ppcs_ptr->first_pass_ref_count) {
        context_ptr->me_context_ptr->me_ds_ref_array[0][0] =
            ppcs_ptr->first_pass_ref_ppcs_ptr[0]->ds_pics;
        context_ptr->me_context_ptr->num_of_ref_pic_to_search[0]++;
    }
    if (ppcs_ptr->first_pass_ref_count > 1) {
        context_ptr->me_context_ptr->me_ds_ref_array[0][1] =
            ppcs_ptr->first_pass_ref_ppcs_ptr[1]->ds_pics;
        context_ptr->me_context_ptr->num_of_ref_pic_to_search[0]++;
    }

    context_ptr->me_context_ptr->me_type = ME_FIRST_PASS;
    // Set 1/4 and 1/16 ME reference buffer(s); filtered or decimated
    EbPictureBufferDesc *quarter_pic_ptr = ppcs_ptr->ds_pics.quarter_picture_ptr;

    EbPictureBufferDesc *sixteenth_pic_ptr = ppcs_ptr->ds_pics.sixteenth_picture_ptr;
    // Parts from MotionEstimationKernel()
    uint32_t sb_origin_x = (uint32_t)(blk_col * BLOCK_SIZE_64);
    uint32_t sb_origin_y = (uint32_t)(blk_row * BLOCK_SIZE_64);

    uint32_t sb_width  = (input_picture_ptr->width - sb_origin_x) < BLOCK_SIZE_64
         ? input_picture_ptr->width - sb_origin_x
         : BLOCK_SIZE_64;
    uint32_t sb_height = (input_picture_ptr->height - sb_origin_y) < BLOCK_SIZE_64
        ? input_picture_ptr->height - sb_origin_y
        : BLOCK_SIZE_64;
    // Load the SB from the input to the intermediate SB buffer
    int buffer_index = (input_picture_ptr->origin_y + sb_origin_y) * input_picture_ptr->stride_y +
        input_picture_ptr->origin_x + sb_origin_x;

    // set search method
    context_ptr->me_context_ptr->hme_search_method = SUB_SAD_SEARCH;

#ifdef ARCH_X86_64
    uint8_t *src_ptr = &(input_picture_ptr->buffer_y[buffer_index]);
    //_MM_HINT_T0     //_MM_HINT_T1    //_MM_HINT_T2    //_MM_HINT_NTA
    uint32_t i;
    for (i = 0; i < sb_height; i++) {
        char const *p = (char const *)(src_ptr + i * input_picture_ptr->stride_y);
        _mm_prefetch(p, _MM_HINT_T2);
    }
#endif
    context_ptr->me_context_ptr->sb_src_ptr    = &(input_picture_ptr->buffer_y[buffer_index]);
    context_ptr->me_context_ptr->sb_src_stride = input_picture_ptr->stride_y;

    // Load the 1/4 decimated SB from the 1/4 decimated input to the 1/4 intermediate SB buffer
    buffer_index = (quarter_pic_ptr->origin_y + (sb_origin_y >> ss_y)) * quarter_pic_ptr->stride_y +
        quarter_pic_ptr->origin_x + (sb_origin_x >> ss_x);

    for (uint32_t sb_row = 0; sb_row < (sb_height >> ss_y); sb_row++) {
        EB_MEMCPY((&(context_ptr->me_context_ptr->quarter_sb_buffer
                         [sb_row * context_ptr->me_context_ptr->quarter_sb_buffer_stride])),
                  (&(quarter_pic_ptr->buffer_y[buffer_index + sb_row * quarter_pic_ptr->stride_y])),
                  (sb_width >> ss_x) * sizeof(uint8_t));
    }

    // Load the 1/16 decimated SB from the 1/16 decimated input to the 1/16 intermediate SB buffer
    buffer_index = (sixteenth_pic_ptr->origin_y + (sb_origin_y >> 2)) *
            sixteenth_pic_ptr->stride_y +
        sixteenth_pic_ptr->origin_x + (sb_origin_x >> 2);

    uint8_t *frame_ptr = &(sixteenth_pic_ptr->buffer_y[buffer_index]);
    uint8_t *local_ptr = context_ptr->me_context_ptr->sixteenth_sb_buffer;

    if (context_ptr->me_context_ptr->hme_search_method == FULL_SAD_SEARCH) {
        for (uint32_t sb_row = 0; sb_row < (sb_height >> 2); sb_row += 1) {
            EB_MEMCPY(local_ptr, frame_ptr, (sb_width >> 2) * sizeof(uint8_t));
            local_ptr += 16;
            frame_ptr += sixteenth_pic_ptr->stride_y;
        }
    } else {
        for (uint32_t sb_row = 0; sb_row < (sb_height >> 2); sb_row += 2) {
            EB_MEMCPY(local_ptr, frame_ptr, (sb_width >> 2) * sizeof(uint8_t));
            local_ptr += 16;
            frame_ptr += sixteenth_pic_ptr->stride_y << 1;
        }
    }
}
/***************************************************************************
* Perform the motion estimation for first pass.
***************************************************************************/
static EbErrorType first_pass_me(PictureParentControlSet *  ppcs_ptr,
                                 MotionEstimationContext_t *me_context_ptr, int32_t segment_index) {
    EbPictureBufferDesc *input_picture_ptr = ppcs_ptr->enhanced_picture_ptr;

    uint32_t blk_cols = (uint32_t)(input_picture_ptr->width + BLOCK_SIZE_64 - 1) / BLOCK_SIZE_64;
    uint32_t blk_rows = (uint32_t)(input_picture_ptr->height + BLOCK_SIZE_64 - 1) / BLOCK_SIZE_64;
    uint32_t ss_x     = ppcs_ptr->scs_ptr->subsampling_x;
    uint32_t ss_y     = ppcs_ptr->scs_ptr->subsampling_y;

    MeContext *context_ptr = me_context_ptr->me_context_ptr;

    uint32_t x_seg_idx;
    uint32_t y_seg_idx;
    uint32_t picture_width_in_b64  = blk_cols;
    uint32_t picture_height_in_b64 = blk_rows;
    SEGMENT_CONVERT_IDX_TO_XY(
        segment_index, x_seg_idx, y_seg_idx, ppcs_ptr->first_pass_seg_column_count);
    uint32_t x_b64_start_idx = SEGMENT_START_IDX(
        x_seg_idx, picture_width_in_b64, ppcs_ptr->first_pass_seg_column_count);
    uint32_t x_b64_end_idx = SEGMENT_END_IDX(
        x_seg_idx, picture_width_in_b64, ppcs_ptr->first_pass_seg_column_count);
    uint32_t y_b64_start_idx = SEGMENT_START_IDX(
        y_seg_idx, picture_height_in_b64, ppcs_ptr->first_pass_seg_row_count);
    uint32_t y_b64_end_idx = SEGMENT_END_IDX(
        y_seg_idx, picture_height_in_b64, ppcs_ptr->first_pass_seg_row_count);

    for (uint32_t blk_row = y_b64_start_idx; blk_row < y_b64_end_idx; blk_row++) {
        for (uint32_t blk_col = x_b64_start_idx; blk_col < x_b64_end_idx; blk_col++) {
            // Initialize ME context
            first_pass_setup_me_context(
                me_context_ptr, ppcs_ptr, input_picture_ptr, blk_row, blk_col, ss_x, ss_y);
            // Perform ME - context_ptr will store the outputs (MVs, buffers, etc)
            // Block-based MC using open-loop HME + refinement
            motion_estimate_sb(ppcs_ptr, // source picture control set -> references come from here
                               (uint32_t)blk_row * blk_cols + blk_col,
                               (uint32_t)blk_col * BLOCK_SIZE_64, // x block
                               (uint32_t)blk_row * BLOCK_SIZE_64, // y block
                               context_ptr,
                               input_picture_ptr); // source picture
        }
    }
    return EB_ErrorNone;
}

/************************************************************************************
* Performs the first pass based on open loop data.
* Source frames are used for Intra and Inter prediction.
* ME is done per segment but the remaining parts performed per frame.
************************************************************************************/
void open_loop_first_pass(PictureParentControlSet *  ppcs_ptr,
                          MotionEstimationContext_t *me_context_ptr, int32_t segment_index) {
    me_context_ptr->me_context_ptr->min_frame_size = MIN(ppcs_ptr->aligned_height,
                                                         ppcs_ptr->aligned_width);
    // Perform the me for the first pass for each segment
    if (ppcs_ptr->first_pass_ref_count)
        first_pass_me(ppcs_ptr, me_context_ptr, segment_index);

    svt_block_on_mutex(ppcs_ptr->first_pass_mutex);
    ppcs_ptr->first_pass_seg_acc++;
    if (ppcs_ptr->first_pass_seg_acc == ppcs_ptr->first_pass_seg_total_count) {
        setup_firstpass_data(ppcs_ptr);
        // Perform the processing of the frame for each frame after me is done for all blocks
        first_pass_frame(ppcs_ptr);

        first_pass_frame_end(ppcs_ptr, ppcs_ptr->ts_duration);
        if (ppcs_ptr->end_of_sequence_flag && !ppcs_ptr->scs_ptr->lap_enabled)
            svt_av1_end_first_pass(ppcs_ptr);
        // Signal that the first pass is done
        svt_post_semaphore(ppcs_ptr->first_pass_done_semaphore);
    }

    svt_release_mutex(ppcs_ptr->first_pass_mutex);
}
