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
#include "EbPictureDecisionResults.h"

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
#if FIX_INVALID_PTR_AFTER_REALLOC
static EbErrorType realloc_stats_out(SequenceControlSet *scs_ptr, uint64_t frame_number) {
    FirstPassStatsOut *out = &scs_ptr->encode_context_ptr->stats_out;
#else
static EbErrorType realloc_stats_out(SequenceControlSet *scs_ptr, FirstPassStatsOut *out,
                                     uint64_t frame_number) {
#endif
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
                stats_in_end_offset   = scs_ptr->twopass.stats_buf_ctx->stats_in_end_write - out->stat;
            }
            EB_REALLOC_ARRAY(out->stat, capability);
            // restore the pointers after re-allocation is done
            scs_ptr->twopass.stats_buf_ctx->stats_in_start = out->stat + stats_in_start_offset;
            scs_ptr->twopass.stats_in                      = out->stat + stats_in_offset;
            scs_ptr->twopass.stats_buf_ctx->stats_in_end_write = out->stat + stats_in_end_offset;
        } else {
            EB_REALLOC_ARRAY(out->stat, capability);
        }
#if FIX_INVALID_PTR_AFTER_REALLOC
        memset(out->stat + out->capability, 0, sizeof(*out->stat) * (capability - out->capability));
#endif
        out->capability = capability;
    }
    out->size = frame_number + 1;
    return EB_ErrorNone;
}

#if FIX_INVALID_PTR_AFTER_REALLOC
static AOM_INLINE void output_stats(SequenceControlSet *scs_ptr, STATS_BUFFER_CTX *stats_buf_ctx,
                                    size_t offset, uint64_t frame_number) {
#else
static AOM_INLINE void output_stats(SequenceControlSet *scs_ptr, FIRSTPASS_STATS *stats,
                                    uint64_t frame_number) {
    FirstPassStatsOut *stats_out = &scs_ptr->encode_context_ptr->stats_out;
#endif
    svt_block_on_mutex(scs_ptr->encode_context_ptr->stat_file_mutex);
#if FIX_INVALID_PTR_AFTER_REALLOC
    if (realloc_stats_out(scs_ptr, frame_number) != EB_ErrorNone) {
#else
    if (realloc_stats_out(scs_ptr, stats_out, frame_number) != EB_ErrorNone) {
#endif
        SVT_ERROR("realloc_stats_out request %d entries failed failed\n", frame_number);
#if !FIX_INVALID_PTR_AFTER_REALLOC
    } else {
        stats_out->stat[frame_number] = *stats;
#endif
    }
#if FIX_INVALID_PTR_AFTER_REALLOC
    FIRSTPASS_STATS *stats        = *(FIRSTPASS_STATS **)((uint8_t *)stats_buf_ctx + offset);
    scs_ptr->encode_context_ptr->stats_out.stat[frame_number] = *stats;
#endif

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
    section->coded_error              = 0.0;
    section->sr_coded_error           = 0.0;
    section->pcnt_inter               = 0.0;
    section->pcnt_motion              = 0.0;
    section->pcnt_second_ref          = 0.0;
    section->pcnt_neutral             = 0.0;
    section->intra_skip_pct           = 0.0;
    section->inactive_zone_rows       = 0.0;
    section->inactive_zone_cols       = 0.0;
#if !CLN_2PASS
    section->MVr                      = 0.0;
#endif
    section->mvr_abs                  = 0.0;
#if !CLN_2PASS
    section->MVc                      = 0.0;
#endif
    section->mvc_abs                  = 0.0;
#if !CLN_2PASS
    section->MVrv                     = 0.0;
    section->MVcv                     = 0.0;
#endif
    section->mv_in_out_count          = 0.0;
#if !CLN_2PASS
    section->new_mv_count             = 0.0;
#endif
    section->count                    = 0.0;
    section->duration                 = 1.0;
#if !CLN_2PASS
    section->raw_error_stdev          = 0.0;
    section->pcnt_third_ref           = 0.0;
    section->tr_coded_error           = 0.0;
#endif
#if FTR_NEW_MULTI_PASS
    memset(&section->stat_struct, 0, sizeof(StatStruct));
#endif
}
void svt_av1_accumulate_stats(FIRSTPASS_STATS *section, const FIRSTPASS_STATS *frame) {
    section->frame += frame->frame;
    section->weight += frame->weight;
    section->intra_error += frame->intra_error;
    section->coded_error += frame->coded_error;
    section->sr_coded_error += frame->sr_coded_error;
    section->pcnt_inter += frame->pcnt_inter;
    section->pcnt_motion += frame->pcnt_motion;
    section->pcnt_second_ref += frame->pcnt_second_ref;
    section->pcnt_neutral += frame->pcnt_neutral;
    section->intra_skip_pct += frame->intra_skip_pct;
    section->inactive_zone_rows += frame->inactive_zone_rows;
    section->inactive_zone_cols += frame->inactive_zone_cols;
#if !CLN_2PASS
    section->MVr += frame->MVr;
#endif
    section->mvr_abs += frame->mvr_abs;
#if !CLN_2PASS
    section->MVc += frame->MVc;
#endif
    section->mvc_abs += frame->mvc_abs;
#if !CLN_2PASS
    section->MVrv += frame->MVrv;
    section->MVcv += frame->MVcv;
#endif
    section->mv_in_out_count += frame->mv_in_out_count;
#if !CLN_2PASS
    section->new_mv_count += frame->new_mv_count;
#endif
    section->count += frame->count;
    section->duration += frame->duration;
}
void svt_av1_end_first_pass(PictureParentControlSet *pcs_ptr) {
#if FIX_INVALID_PTR_AFTER_REALLOC
    if (pcs_ptr->scs_ptr->twopass.stats_buf_ctx->total_stats) {
#else
    SequenceControlSet *scs_ptr = pcs_ptr->scs_ptr;
    TWO_PASS *          twopass = &scs_ptr->twopass;

    if (twopass->stats_buf_ctx->total_stats) {
#endif
        // add the total to the end of the file
#if FIX_INVALID_PTR_AFTER_REALLOC
        output_stats(pcs_ptr->scs_ptr,
                     pcs_ptr->scs_ptr->twopass.stats_buf_ctx,
                     offsetof(STATS_BUFFER_CTX, total_stats),
                     pcs_ptr->picture_number + 1);
#else
        output_stats(scs_ptr, twopass->stats_buf_ctx->total_stats, pcs_ptr->picture_number + 1);
#endif
    }
}
#if !CLN_2PASS
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
#endif
#define UL_INTRA_THRESH 50
#define INVALID_ROW -1
// Accumulates motion vector stats.
// Modifies member variables of "stats".
void accumulate_mv_stats(const MV best_mv, const FULLPEL_MV mv, const int mb_row, const int mb_col,
                         const int mb_rows, const int mb_cols, MV *last_mv, FRAME_STATS *stats) {
    if (is_zero_mv(&best_mv))
        return;

    ++stats->mv_count;
#if !CLN_2PASS
    // Non-zero vector, was it different from the last non zero vector?
    if (!is_equal_mv(&best_mv, last_mv))
        ++stats->new_mv_count;
#endif

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
#if CLN_2PASS
#if !FIX_INVALID_PTR_AFTER_REALLOC
                                   const int frame_number,
#endif
#else
#if !FIX_INVALID_PTR_AFTER_REALLOC
                                   const double raw_err_stdev, const int frame_number,
#else
                                   const double raw_err_stdev,
#endif
#endif
#if FTR_2PASS_1PASS_UNIFICATION
                                   const double ts_duration) {
#else
                                   const int64_t ts_duration) {
#endif
    SequenceControlSet *scs_ptr = pcs_ptr->scs_ptr;
    TWO_PASS *          twopass = &scs_ptr->twopass;

    const uint32_t   mb_cols          = (scs_ptr->seq_header.max_frame_width + 16 - 1) / 16;
    const uint32_t   mb_rows          = (scs_ptr->seq_header.max_frame_height + 16 - 1) / 16;
    FIRSTPASS_STATS *this_frame_stats = twopass->stats_buf_ctx->stats_in_end_write;
    FIRSTPASS_STATS  fps;
    // The minimum error here insures some bit allocation to frames even
    // in static regions. The allocation per MB declines for larger formats
    // where the typical "real" energy per MB also falls.
    // Initial estimate here uses sqrt(mbs) to define the min_err, where the
    // number of mbs is proportional to the image area.
#if OPT_FIRST_PASS3
    uint32_t step = (uint32_t)pcs_ptr->bypass_blk_step;
    const int num_mbs = (mb_rows / step) * (mb_cols/ step);
#else
    const int num_mbs = mb_rows * mb_cols;
#endif
    //(cpi->oxcf.resize_cfg.resize_mode != RESIZE_NONE)
    //    ? cpi->initial_mbs
    //    : mi_params->MBs;
    const double min_err = 200 * sqrt(num_mbs);

    if (pcs_ptr->skip_frame) {
        FirstPassStatsOut *stats_out = &scs_ptr->encode_context_ptr->stats_out;
#if FIX_INVALID_PTR_AFTER_REALLOC
        fps       = stats_out->stat[pcs_ptr->picture_number - 1];
        fps.frame = (double)pcs_ptr->picture_number;
#else
        fps       = stats_out->stat[frame_number - 1];
        fps.frame = frame_number;
#endif
    }else{
    fps.weight                   = stats->intra_factor * stats->brightness_factor;
#if FIX_INVALID_PTR_AFTER_REALLOC
    fps.frame                    = (double)pcs_ptr->picture_number;
#else
    fps.frame                    = frame_number;
#endif
    fps.coded_error              = (double)(stats->coded_error >> 8) + min_err;
    fps.sr_coded_error           = (double)(stats->sr_coded_error >> 8) + min_err;
#if !CLN_2PASS
    fps.tr_coded_error           = (double)(stats->tr_coded_error >> 8) + min_err;
#endif
    fps.intra_error              = (double)(stats->intra_error >> 8) + min_err;
    fps.count                    = 1.0;
    fps.pcnt_inter               = (double)stats->inter_count / num_mbs;
    fps.pcnt_second_ref          = (double)stats->second_ref_count / num_mbs;
#if !CLN_2PASS
    fps.pcnt_third_ref           = (double)stats->third_ref_count / num_mbs;
#endif
    fps.pcnt_neutral             = (double)stats->neutral_count / num_mbs;
    fps.intra_skip_pct           = (double)stats->intra_skip_count / num_mbs;
    fps.inactive_zone_rows       = (double)stats->image_data_start_row;
    fps.inactive_zone_cols       = (double)0; // TODO(paulwilkins): fix
#if !CLN_2PASS
    fps.raw_error_stdev          = raw_err_stdev;
#endif

    if (stats->mv_count > 0) {
#if !CLN_2PASS
        fps.MVr     = (double)stats->sum_mvr / stats->mv_count;
#endif
        fps.mvr_abs = (double)stats->sum_mvr_abs / stats->mv_count;
#if !CLN_2PASS
        fps.MVc     = (double)stats->sum_mvc / stats->mv_count;
#endif
        fps.mvc_abs = (double)stats->sum_mvc_abs / stats->mv_count;
#if !CLN_2PASS
        fps.MVrv    = ((double)stats->sum_mvrs -
                    ((double)stats->sum_mvr * stats->sum_mvr / stats->mv_count)) /
            stats->mv_count;
        fps.MVcv = ((double)stats->sum_mvcs -
                    ((double)stats->sum_mvc * stats->sum_mvc / stats->mv_count)) /
            stats->mv_count;
#endif
        fps.mv_in_out_count = (double)stats->sum_in_vectors / (stats->mv_count * 2);
#if !CLN_2PASS
        fps.new_mv_count    = stats->new_mv_count;
#endif
        fps.pcnt_motion     = (double)stats->mv_count / num_mbs;
    } else {
#if !CLN_2PASS
        fps.MVr             = 0.0;
#endif
        fps.mvr_abs         = 0.0;
#if !CLN_2PASS
        fps.MVc             = 0.0;
#endif
        fps.mvc_abs         = 0.0;
#if !CLN_2PASS
        fps.MVrv            = 0.0;
        fps.MVcv            = 0.0;
#endif
        fps.mv_in_out_count = 0.0;
#if !CLN_2PASS
        fps.new_mv_count    = 0.0;
#endif
        fps.pcnt_motion     = 0.0;
    }
#if 0//FTR_OPT_IPP_DOWN_SAMPLE
   // fps.weight = stats->intra_factor * stats->brightness_factor;
   // fps.frame = frame_number;
    fps.coded_error *= 3;// (double)(stats->coded_error >> 8) + min_err;
    fps.sr_coded_error *= 3;// (double)(stats->sr_coded_error >> 8) + min_err;
    fps.tr_coded_error *= 3;//(double)(stats->tr_coded_error >> 8) + min_err;
    fps.intra_error *= 3;// (double)(stats->intra_error >> 8) + min_err;
   // fps.count = 1.0;
   // fps.pcnt_inter = (double)stats->inter_count / num_mbs;
  //  fps.pcnt_second_ref = (double)stats->second_ref_count / num_mbs;
   // fps.pcnt_third_ref = (double)stats->third_ref_count / num_mbs;
  //  fps.pcnt_neutral = (double)stats->neutral_count / num_mbs;
  //  fps.intra_skip_pct = (double)stats->intra_skip_count / num_mbs;
    fps.inactive_zone_rows *= 2;// (double)stats->image_data_start_row;
    fps.inactive_zone_cols *= 2;// (double)0; // TODO(paulwilkins): fix
   // fps.raw_error_stdev = raw_err_stdev;

    if (stats->mv_count > 0) {
       // fps.MVr = (double)stats->sum_mvr / stats->mv_count;
      //  fps.mvr_abs = (double)stats->sum_mvr_abs / stats->mv_count;
       // fps.MVc = (double)stats->sum_mvc / stats->mv_count;
     //   fps.mvc_abs = (double)stats->sum_mvc_abs / stats->mv_count;
      //  fps.MVrv = ((double)stats->sum_mvrs -
      //      ((double)stats->sum_mvr * stats->sum_mvr / stats->mv_count)) /
       //     stats->mv_count;
      //  fps.MVcv = ((double)stats->sum_mvcs -
     //       ((double)stats->sum_mvc * stats->sum_mvc / stats->mv_count)) /
     //       stats->mv_count;
     //   fps.mv_in_out_count = (double)stats->sum_in_vectors / (stats->mv_count * 2);
        fps.new_mv_count *= 3;// stats->new_mv_count;
    //    fps.pcnt_motion = (double)stats->mv_count / num_mbs;
    }
    //else {
    //    fps.MVr = 0.0;
    //    fps.mvr_abs = 0.0;
    //    fps.MVc = 0.0;
    //    fps.mvc_abs = 0.0;
    //    fps.MVrv = 0.0;
    //    fps.MVcv = 0.0;
    //    fps.mv_in_out_count = 0.0;
    //    fps.new_mv_count = 0.0;
    //    fps.pcnt_motion = 0.0;
    //}
#endif
    // TODO(paulwilkins):  Handle the case when duration is set to 0, or
    // something less than the full time between subsequent values of
    // cpi->source_time_stamp.
    fps.duration = (double)ts_duration;
    }
#if FTR_NEW_MULTI_PASS
    memset(&fps.stat_struct, 0, sizeof(StatStruct));
#endif
    // We will store the stats inside the persistent twopass struct (and NOT the
    // local variable 'fps'), and then cpi->output_pkt_list will point to it.
    *this_frame_stats = fps;
#if FIX_INVALID_PTR_AFTER_REALLOC
    output_stats(scs_ptr,
                 twopass->stats_buf_ctx,
                 offsetof(STATS_BUFFER_CTX, stats_in_end_write),
                 pcs_ptr->picture_number);
#else
    output_stats(scs_ptr, this_frame_stats, pcs_ptr->picture_number);
#endif
#if FTR_1PAS_VBR
    if (twopass->stats_buf_ctx->total_stats != NULL && use_output_stat(scs_ptr)) {
#else
    if (twopass->stats_buf_ctx->total_stats != NULL) {
#endif
        svt_av1_accumulate_stats(twopass->stats_buf_ctx->total_stats, &fps);
    }
    /*In the case of two pass, first pass uses it as a circular buffer,
   * when LAP is enabled it is used as a linear buffer*/
    twopass->stats_buf_ctx->stats_in_end_write++;

    if ((use_output_stat(scs_ptr)) &&
        (twopass->stats_buf_ctx->stats_in_end_write >= twopass->stats_buf_ctx->stats_in_buf_end)) {
        twopass->stats_buf_ctx->stats_in_end_write = twopass->stats_buf_ctx->stats_in_start;
    }
}

#if OPT_FIRST_PASS3
static FRAME_STATS accumulate_frame_stats(FRAME_STATS *mb_stats, int mb_rows, int mb_cols, PictureParentControlSet *pcs_ptr) {
#else
static FRAME_STATS accumulate_frame_stats(FRAME_STATS *mb_stats, int mb_rows, int mb_cols) {
#endif
    FRAME_STATS stats = {0};
    int         i, j;
#if OPT_FIRST_PASS3
    uint32_t blks_in_b64 = BLOCK_SIZE_64 / FORCED_BLK_SIZE;
#endif
    stats.image_data_start_row = INVALID_ROW;
    for (j = 0; j < mb_rows; j++) {
        for (i = 0; i < mb_cols; i++) {
            FRAME_STATS mb_stat = mb_stats[j * mb_cols + i];
#if OPT_FIRST_PASS3
            int sb_index_x = i / blks_in_b64;
            int sb_index_y = j / blks_in_b64;
            if(pcs_ptr->bypass_blk_step > 1)
                if((sb_index_x% (uint32_t)pcs_ptr->bypass_blk_step != 0) || (sb_index_y % (uint32_t)pcs_ptr->bypass_blk_step != 0))
                    continue;
#endif
            stats.brightness_factor += mb_stat.brightness_factor;
            stats.coded_error += mb_stat.coded_error;
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
#if !CLN_2PASS
            stats.new_mv_count += mb_stat.new_mv_count;
#endif
            stats.second_ref_count += mb_stat.second_ref_count;
            stats.sr_coded_error += mb_stat.sr_coded_error;
            stats.sum_in_vectors += mb_stat.sum_in_vectors;
            stats.sum_mvc += mb_stat.sum_mvc;
            stats.sum_mvc_abs += mb_stat.sum_mvc_abs;
#if !CLN_2PASS
            stats.sum_mvcs += mb_stat.sum_mvcs;
            stats.sum_mvr += mb_stat.sum_mvr;
            stats.sum_mvr_abs += mb_stat.sum_mvr_abs;
            stats.sum_mvrs += mb_stat.sum_mvrs;
            stats.third_ref_count += mb_stat.third_ref_count;
            stats.tr_coded_error += mb_stat.tr_coded_error;
#endif
        }
    }
    return stats;
}
/**************************************************
 * Reset first pass stat
 **************************************************/
void setup_firstpass_data_seg(PictureParentControlSet *ppcs_ptr, int32_t segment_index) {
    SequenceControlSet *scs_ptr        = ppcs_ptr->scs_ptr;
    FirstPassData *     firstpass_data = &ppcs_ptr->firstpass_data;
    const uint32_t      mb_cols        = (scs_ptr->seq_header.max_frame_width + 16 - 1) / 16;
    const uint32_t      mb_rows        = (scs_ptr->seq_header.max_frame_height + 16 - 1) / 16;
    EbPictureBufferDesc *input_picture_ptr = ppcs_ptr->enhanced_picture_ptr;

    uint32_t blk_cols = (uint32_t)(input_picture_ptr->width + BLOCK_SIZE_64 - 1) / BLOCK_SIZE_64;
    uint32_t blk_rows = (uint32_t)(input_picture_ptr->height + BLOCK_SIZE_64 - 1) / BLOCK_SIZE_64;

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

    const uint32_t mb_y_end = (y_b64_end_idx << 2) > mb_rows ? mb_rows : (y_b64_end_idx << 2);
    const uint32_t mb_x_end = (x_b64_end_idx << 2) > mb_cols ? mb_cols : (x_b64_end_idx << 2);

    for (uint32_t mb_y = (y_b64_start_idx << 2); mb_y < mb_y_end; mb_y++) {
        for (uint32_t mb_x = (x_b64_start_idx << 2); mb_x < mb_x_end; mb_x++) {
            memset(firstpass_data->mb_stats + mb_x + mb_y * mb_cols, 0, sizeof(*firstpass_data->mb_stats));
            firstpass_data->mb_stats[mb_x + mb_y * mb_cols].image_data_start_row = INVALID_ROW;
        }
    }
}
#if FTR_2PASS_1PASS_UNIFICATION
void first_pass_frame_end(PictureParentControlSet *pcs_ptr, const double ts_duration) {
#else
void first_pass_frame_end(PictureParentControlSet *pcs_ptr, const int64_t ts_duration) {
#endif
    SequenceControlSet *scs_ptr = pcs_ptr->scs_ptr;
    const uint32_t      mb_cols = (scs_ptr->seq_header.max_frame_width + 16 - 1) / 16;
    const uint32_t      mb_rows = (scs_ptr->seq_header.max_frame_height + 16 - 1) / 16;

#if !CLN_2PASS
    int *        raw_motion_err_list = pcs_ptr->firstpass_data.raw_motion_err_list;
#endif
    FRAME_STATS *mb_stats            = pcs_ptr->firstpass_data.mb_stats;

    FRAME_STATS  stats;
#if !CLN_2PASS
    double raw_err_stdev = 0;
#endif
    if (!pcs_ptr->skip_frame) {
#if OPT_FIRST_PASS3
        stats                      = accumulate_frame_stats(mb_stats, mb_rows, mb_cols, pcs_ptr);
#else
        stats                      = accumulate_frame_stats(mb_stats, mb_rows, mb_cols);
#endif
#if !CLN_2PASS
        int total_raw_motion_err_count = frame_is_intra_only(pcs_ptr) ? 0 : mb_rows * mb_cols;
        raw_err_stdev              = raw_motion_error_stdev(raw_motion_err_list,
                                                        total_raw_motion_err_count);
#endif
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
#if OPT_FIRST_PASS3
    uint32_t step = (uint32_t)pcs_ptr->bypass_blk_step;
    const int num_mbs = (mb_rows/ step) * (mb_cols/ step);
#else
    const int num_mbs = mb_rows * mb_cols;
#endif
    /*(cpi->oxcf.resize_cfg.resize_mode != RESIZE_NONE)
        ? cpi->initial_mbs
        : mi_params->MBs;*/
    stats.intra_factor      = stats.intra_factor / (double)num_mbs;
    stats.brightness_factor = stats.brightness_factor / (double)num_mbs;
    }
    update_firstpass_stats(
#if CLN_2PASS
#if FIX_INVALID_PTR_AFTER_REALLOC
        pcs_ptr, &stats, ts_duration);
#else
        pcs_ptr, &stats, (const int)pcs_ptr->picture_number, ts_duration);
#endif
#else
#if FIX_INVALID_PTR_AFTER_REALLOC
        pcs_ptr, &stats, raw_err_stdev, ts_duration);
#else
        pcs_ptr, &stats, raw_err_stdev, (const int)pcs_ptr->picture_number, ts_duration);
#endif
#endif
}
#if TUNE_MULTI_PASS
void samepred_pass_frame_end(PictureParentControlSet *pcs_ptr, const double ts_duration) {
    SequenceControlSet *scs_ptr = pcs_ptr->scs_ptr;
    TWO_PASS *          twopass = &scs_ptr->twopass;

    FIRSTPASS_STATS *this_frame_stats = twopass->stats_buf_ctx->stats_in_end_write;
    FIRSTPASS_STATS  fps ;

    fps.intra_error= 0;
    fps.coded_error= 0;
    fps.sr_coded_error= 0;
    fps.pcnt_inter= 0;
    fps.pcnt_motion= 0;
    fps.pcnt_second_ref= 0;
    fps.pcnt_neutral= 0;
    fps.intra_skip_pct= 0;
    fps.inactive_zone_rows= 0;
    fps.inactive_zone_cols= 0;
#if !CLN_2PASS
    fps.MVr= 0;
#endif
    fps.mvr_abs= 0;
#if !CLN_2PASS
    fps.MVc= 0;
#endif
    fps.mvc_abs= 0;
#if !CLN_2PASS
    fps.MVrv= 0;
    fps.MVcv= 0;
#endif
    fps.mv_in_out_count= 0;
#if !CLN_2PASS
    fps.new_mv_count= 0;
#endif

    fps.weight = 0;
    fps.duration = (double)ts_duration;
    fps.count = 1.0;
    fps.frame = (const int)pcs_ptr->picture_number;
    memset(&fps.stat_struct, 0, sizeof(StatStruct));
    fps.stat_struct.poc = pcs_ptr->picture_number;
#if FTR_OPT_MPASS_DOWN_SAMPLE
    if (is_middle_pass_ds(scs_ptr))
        fps.stat_struct.total_num_bits = pcs_ptr->total_num_bits *DS_SC_FACT / 10;
    else
        fps.stat_struct.total_num_bits = pcs_ptr->total_num_bits;
#else
    fps.stat_struct.total_num_bits = pcs_ptr->total_num_bits;
#endif
    fps.stat_struct.qindex = pcs_ptr->frm_hdr.quantization_params.base_q_idx;
    fps.stat_struct.worst_qindex = quantizer_to_qindex[(uint8_t)scs_ptr->static_config.qp];
#if  FTR_OPT_MPASS_BYPASS_FRAMES
    fps.stat_struct.temporal_layer_index = pcs_ptr->temporal_layer_index;
#endif

    *this_frame_stats = fps;
#if FIX_INVALID_PTR_AFTER_REALLOC
    output_stats(scs_ptr,
                 twopass->stats_buf_ctx,
                 offsetof(STATS_BUFFER_CTX, stats_in_end_write),
                 pcs_ptr->picture_number);
#else
    output_stats(scs_ptr, this_frame_stats, pcs_ptr->picture_number);
#endif
    if (twopass->stats_buf_ctx->total_stats != NULL) {
        svt_av1_accumulate_stats(twopass->stats_buf_ctx->total_stats, &fps);
    }
    /*In the case of two pass, first pass uses it as a circular buffer,
   * when LAP is enabled it is used as a linear buffer*/
    twopass->stats_buf_ctx->stats_in_end_write++;

    if ((twopass->stats_buf_ctx->stats_in_end_write >= twopass->stats_buf_ctx->stats_in_buf_end))
        twopass->stats_buf_ctx->stats_in_end_write = twopass->stats_buf_ctx->stats_in_start;
    if (pcs_ptr->end_of_sequence_flag){
        if (twopass->stats_buf_ctx->total_stats) {
            // add the total to the end of the file
#if FIX_INVALID_PTR_AFTER_REALLOC
            output_stats(scs_ptr,
                         twopass->stats_buf_ctx,
                         offsetof(STATS_BUFFER_CTX, total_stats),
                         scs_ptr->encode_context_ptr->terminating_picture_number + 1);
#else
            output_stats(scs_ptr,
                         twopass->stats_buf_ctx->total_stats,
                         scs_ptr->encode_context_ptr->terminating_picture_number + 1);
#endif
        }
    }
}
#endif
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
#if OPT_FIRST_PASS3
    pcs_ptr->bypass_blk_step = pcs_ptr->scs_ptr->static_config.ipp_ctrls.bypass_blk_step ? ((pcs_ptr->scs_ptr->input_resolution < INPUT_SIZE_360p_RANGE) ? 1 : 2) : 1;
#endif
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
#define MOTION_ERROR_THRESH 500
void set_tf_controls(PictureParentControlSet *pcs_ptr, uint8_t tf_level);
#if FTR_NEW_WN_LVLS
void set_wn_filter_ctrls(Av1Common* cm, uint8_t wn_filter_lvl);
#endif
#if CLN_DLF_SIGNALS
void set_dlf_controls(PictureParentControlSet *pcs_ptr, uint8_t dlf_level);
#endif
/******************************************************
* Derive Multi-Processes Settings for first pass
Input   : encoder mode and tune
Output  : Multi-Processes signal(s)
******************************************************/

EbErrorType first_pass_signal_derivation_multi_processes(SequenceControlSet *     scs_ptr,
                                                         PictureParentControlSet *pcs_ptr) {
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
#if CLN_DLF_SIGNALS
    set_dlf_controls(pcs_ptr, 0);
#else
    // Loop filter Level                            Settings
    // 0                                            OFF
    // 1                                            CU-BASED
    // 2                                            LIGHT FRAME-BASED
    // 3                                            FULL FRAME-BASED
    pcs_ptr->loop_filter_mode = 0;
#endif
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

#if FTR_NEW_WN_LVLS
    set_wn_filter_ctrls(cm, 0);
#else
    // WN Level                                     Settings
    // 0                                            OFF
    // 1                                            3-Tap luma/ 3-Tap chroma
    // 2                                            5-Tap luma/ 5-Tap chroma
    // 3                                            7-Tap luma/ 5-Tap chroma
    cm->wn_filter_mode = 0;
#endif

#if !TUNE_INTRA_LEVELS
    // Intra prediction modes                       Settings
    // 0                                            FULL
    // 1                                            LIGHT per block : disable_z2_prediction && disable_angle_refinement  for 64/32/4
    // 2                                            OFF per block : disable_angle_prediction for 64/32/4
    // 3                                            OFF : disable_angle_prediction
    // 4                                            OIS based Intra
    // 5                                            Light OIS based Intra
#endif
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
#if !OPT_TXS_SEARCH
    // Exit TX size search when all coefficients are zero
    // 0: OFF
    // 1: ON
    pcs_ptr->tx_size_early_exit = 0;
#endif
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
#if FIX_REMOVE_PD1
    //  pic_obmc_level  | Default Encoder Settings
    //         0        | OFF subject to possible constraints
    //       > 1        | Faster level subject to possible constraints
#else
    //  pic_obmc_level  |              Default Encoder Settings             |     Command Line Settings
    //         0        | OFF subject to possible constraints               | OFF everywhere in encoder
    //         1        | ON subject to possible constraints                | Fully ON in PD_PASS_2
    //         2        | Faster level subject to possible constraints      | Level 2 everywhere in PD_PASS_2
    //         3        | Even faster level subject to possible constraints | Level 3 everywhere in PD_PASS_3
#endif
    pcs_ptr->parent_pcs_ptr->pic_obmc_level = 0;

    // Switchable Motion Mode
    frm_hdr->is_motion_mode_switchable = frm_hdr->is_motion_mode_switchable ||
        pcs_ptr->parent_pcs_ptr->pic_obmc_level;

    // HBD Mode
    pcs_ptr->hbd_mode_decision = EB_8_BIT_MD; //first pass hard coded to 8bit
#if OPT9_RATE_ESTIMATION
    pcs_ptr->parent_pcs_ptr->partition_contexts = PARTITION_CONTEXTS;
#endif
    pcs_ptr->parent_pcs_ptr->bypass_cost_table_gen = 0;
#if  FTR_SIMPLIFIED_MV_COST
#if CLN_RATE_EST_CTRLS
    pcs_ptr->approx_inter_rate = 0;
#else
    pcs_ptr->use_low_precision_cost_estimation = 0;
#endif
#endif
    return return_error;
}
#if FTR_2PASS_1PASS_UNIFICATION
/************************************************
 * Set ME/HME Params for the first pass encoding
 ************************************************/
void *set_first_pass_me_hme_params_oq(MeContext *me_context_ptr, PictureParentControlSet *pcs_ptr,
    SequenceControlSet *scs_ptr, EbInputResolution input_resolution) {
    UNUSED(scs_ptr);

#if FTR_BIAS_STAT
    me_context_ptr->stat_factor = 100;
#endif
    // HME/ME default settings
    me_context_ptr->number_hme_search_region_in_width = 2;
    me_context_ptr->number_hme_search_region_in_height = 2;
    me_context_ptr->reduce_hme_l0_sr_th_min = 0;
    me_context_ptr->reduce_hme_l0_sr_th_max = 0;

    // Set the minimum ME search area
#if IPP_CTRL
    if (!scs_ptr->static_config.ipp_ctrls.reduce_me_search) {
#else
    if (pcs_ptr->scs_ptr->enc_mode_2ndpass <= ENC_M4) {
#endif
        me_context_ptr->search_area_width = me_context_ptr->search_area_height = 8;
        me_context_ptr->max_me_search_width = me_context_ptr->max_me_search_height = 8;
    }
    else {
        me_context_ptr->search_area_width = 8;
        me_context_ptr->search_area_height = 3;
        me_context_ptr->max_me_search_width = 8;
        me_context_ptr->max_me_search_height = 5;
    }

    if (pcs_ptr->input_resolution < INPUT_SIZE_1080p_RANGE) {
        me_context_ptr->hme_level0_total_search_area_width = me_context_ptr->hme_level0_total_search_area_height = 8;
        me_context_ptr->hme_level0_max_total_search_area_width = me_context_ptr->hme_level0_max_total_search_area_height = 192;
    }
    else {
        me_context_ptr->hme_level0_total_search_area_width = me_context_ptr->hme_level0_total_search_area_height = 16;
        me_context_ptr->hme_level0_max_total_search_area_width = me_context_ptr->hme_level0_max_total_search_area_height = 192;
    }
    me_context_ptr->hme_level0_total_search_area_width = me_context_ptr->hme_level0_total_search_area_height = me_context_ptr->hme_level0_total_search_area_width / 2;
    me_context_ptr->hme_level0_max_total_search_area_width = me_context_ptr->hme_level0_max_total_search_area_height = me_context_ptr->hme_level0_max_total_search_area_width / 2;

    me_context_ptr->hme_level0_max_search_area_in_width_array[0] =
        me_context_ptr->hme_level0_max_search_area_in_width_array[1] =
        me_context_ptr->hme_level0_max_total_search_area_width / me_context_ptr->number_hme_search_region_in_width;
    me_context_ptr->hme_level0_max_search_area_in_height_array[0] =
        me_context_ptr->hme_level0_max_search_area_in_height_array[1] =
        me_context_ptr->hme_level0_max_total_search_area_height / me_context_ptr->number_hme_search_region_in_height;
    me_context_ptr->hme_level0_search_area_in_width_array[0] =
        me_context_ptr->hme_level0_search_area_in_width_array[1] =
        me_context_ptr->hme_level0_total_search_area_width / me_context_ptr->number_hme_search_region_in_width;
    me_context_ptr->hme_level0_search_area_in_height_array[0] =
        me_context_ptr->hme_level0_search_area_in_height_array[1] =
        me_context_ptr->hme_level0_total_search_area_height / me_context_ptr->number_hme_search_region_in_height;

    me_context_ptr->hme_level1_search_area_in_width_array[0] =
        me_context_ptr->hme_level1_search_area_in_width_array[1] = 8;
    me_context_ptr->hme_level1_search_area_in_height_array[0] =
        me_context_ptr->hme_level1_search_area_in_height_array[1] = 3;

    me_context_ptr->hme_level2_search_area_in_width_array[0] =
        me_context_ptr->hme_level2_search_area_in_width_array[1] = 8;

    me_context_ptr->hme_level2_search_area_in_height_array[0] =
        me_context_ptr->hme_level2_search_area_in_height_array[1] = 3;

    me_context_ptr->hme_level1_search_area_in_width_array[0] =
        me_context_ptr->hme_level1_search_area_in_width_array[1] =
        me_context_ptr->hme_level1_search_area_in_height_array[0] =
        me_context_ptr->hme_level1_search_area_in_height_array[1] = 16 / 2;

    me_context_ptr->hme_level2_search_area_in_width_array[0] =
        me_context_ptr->hme_level2_search_area_in_width_array[1] =
        me_context_ptr->hme_level2_search_area_in_height_array[0] =
        me_context_ptr->hme_level2_search_area_in_height_array[1] = 16 / 2;

    if (input_resolution <= INPUT_SIZE_720p_RANGE)
        me_context_ptr->hme_decimation = pcs_ptr->enc_mode <= ENC_MRS ? ONE_DECIMATION_HME : TWO_DECIMATION_HME;
    else
        me_context_ptr->hme_decimation = TWO_DECIMATION_HME;

    // Scale up the MIN ME area if low frame rate
    uint8_t  low_frame_rate_flag = (scs_ptr->static_config.frame_rate >> 16) < 50 ? 1 : 0;
    if (low_frame_rate_flag) {
        me_context_ptr->search_area_width = (me_context_ptr->search_area_width * 3) / 2;
        me_context_ptr->search_area_height = (me_context_ptr->search_area_height * 3) / 2;
    }
#if FTR_HME_ME_EARLY_EXIT
    me_context_ptr->me_early_exit_th = 0;
#endif
    return NULL;
};
#else
void *set_me_hme_params_oq(MeContext *me_context_ptr, PictureParentControlSet *pcs_ptr,
                           SequenceControlSet *scs_ptr, EbInputResolution input_resolution);
#endif
void *set_me_hme_params_from_config(SequenceControlSet *scs_ptr, MeContext *me_context_ptr);
void  set_me_hme_ref_prune_ctrls(MeContext *context_ptr, uint8_t prune_level);
void  set_me_sr_adjustment_ctrls(MeContext *context_ptr, uint8_t sr_adjustment_level);
void  set_gm_controls(PictureParentControlSet *pcs_ptr, uint8_t gm_level);
void  set_prehme_ctrls(MeContext* context, uint8_t level);
/******************************************************
* Derive ME Settings for first pass
  Input   : encoder mode and tune
  Output  : ME Kernel signal(s)
******************************************************/
EbErrorType first_pass_signal_derivation_me_kernel(SequenceControlSet *       scs_ptr,
                                                   PictureParentControlSet *  pcs_ptr,
                                                   MotionEstimationContext_t *context_ptr) {
    EbErrorType return_error = EB_ErrorNone;
#if FTR_TUNE_PRUNING
    EbInputResolution input_resolution = scs_ptr->input_resolution;
    context_ptr->me_context_ptr->input_resolution = input_resolution;
    context_ptr->me_context_ptr->clip_class = pcs_ptr->sc_class1;
#endif
    // Set ME/HME search regions

    if (scs_ptr->static_config.use_default_me_hme)
#if FTR_2PASS_1PASS_UNIFICATION
        set_first_pass_me_hme_params_oq(
            context_ptr->me_context_ptr, pcs_ptr, scs_ptr, scs_ptr->input_resolution);
#else
        set_me_hme_params_oq(
            context_ptr->me_context_ptr, pcs_ptr, scs_ptr, scs_ptr->input_resolution);
#endif
    else
        set_me_hme_params_from_config(scs_ptr, context_ptr->me_context_ptr);

    // Set HME flags
    context_ptr->me_context_ptr->enable_hme_flag        = pcs_ptr->enable_hme_flag;
    context_ptr->me_context_ptr->enable_hme_level0_flag = pcs_ptr->enable_hme_level0_flag;
    context_ptr->me_context_ptr->enable_hme_level1_flag = pcs_ptr->enable_hme_level1_flag;
#if TUNE_FIRSTPASS_HME2
#if IPP_CTRL
    context_ptr->me_context_ptr->enable_hme_level2_flag = !scs_ptr->static_config.ipp_ctrls.reduce_me_search ? pcs_ptr->enable_hme_level2_flag : 0;
#else
    context_ptr->me_context_ptr->enable_hme_level2_flag = scs_ptr->enc_mode_2ndpass <= ENC_M4 ? pcs_ptr->enable_hme_level2_flag : 0;
#endif
#else
    context_ptr->me_context_ptr->enable_hme_level2_flag = scs_ptr->enc_mode_2ndpass <= ENC_M7 ? pcs_ptr->enable_hme_level2_flag : 0;
#endif

    // HME Search Method
    context_ptr->me_context_ptr->hme_search_method = SUB_SAD_SEARCH;

    // ME Search Method
    context_ptr->me_context_ptr->me_search_method = SUB_SAD_SEARCH;

#if TUNE_HME_SUB
    // skip search line hme level0
    context_ptr->me_context_ptr->skip_search_line_hme0 = 0;
#endif

    uint8_t gm_level                              = 0;
    set_gm_controls(pcs_ptr, gm_level);

    // Set pre-hme level (0-2)
    uint8_t prehme_level = 0;
    set_prehme_ctrls(context_ptr->me_context_ptr, prehme_level);

    // Set hme/me based reference pruning level (0-4)
    set_me_hme_ref_prune_ctrls(context_ptr->me_context_ptr, 0);

    // Set hme-based me sr adjustment level
    set_me_sr_adjustment_ctrls(context_ptr->me_context_ptr, 0);
    context_ptr->me_context_ptr->prune_me_candidates_th = 0; // No impact on tf
#if FTR_LIMIT_ME_CANDS
    context_ptr->me_context_ptr->use_best_unipred_cand_only = 0; // No impact on tf
#endif
    set_prehme_ctrls(context_ptr->me_context_ptr, 0);
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
static int open_loop_firstpass_intra_prediction(PictureParentControlSet *ppcs_ptr,int raw_motion_err , uint32_t blk_origin_x, uint32_t blk_origin_y,
                                                uint8_t bwidth, uint8_t bheight,
                                                EbPictureBufferDesc *input_picture_ptr,
                                                uint32_t             input_origin_index,
                                                FRAME_STATS *const   stats) {
    int32_t   mb_row      = blk_origin_y >> 4;
    int32_t   mb_col      = blk_origin_x >> 4;
    const int use_dc_pred = (mb_col || mb_row) && (!mb_col || !mb_row);
    uint8_t  use8blk = 0;
#if FIX_PRESET_TUNING
#if IPP_CTRL
    if (!ppcs_ptr->scs_ptr->static_config.ipp_ctrls.use8blk) {
#else
    if (ppcs_ptr->scs_ptr->enc_mode_2ndpass <= ENC_M4) {
#endif
#else
    if (ppcs_ptr->scs_ptr->enc_mode_2ndpass <= ENC_M7) {
#endif
        use8blk = 0;
    }
    else {
        if (ppcs_ptr->first_pass_ref_count)
            if (raw_motion_err > MOTION_ERROR_THRESH)
                use8blk = 0;
            else
                use8blk = 1;
        else
            use8blk = 0;
    }

    uint32_t sub_blk_origin_x, sub_blk_origin_y;
    uint8_t *above_row;
    uint8_t *left_col;

    DECLARE_ALIGNED(16, uint8_t, left_data[MAX_TX_SIZE * 2 + 32]);
    DECLARE_ALIGNED(16, uint8_t, above_data[MAX_TX_SIZE * 2 + 32]);
    DECLARE_ALIGNED(32, uint8_t, predictor8[256 * 2]);
    uint8_t *predictor = predictor8;
    uint8_t sub_blk_rows = use_dc_pred
        ? (bheight == FORCED_BLK_SIZE && bwidth == FORCED_BLK_SIZE) ? 1 : bheight / 8
        : use8blk ? bheight / 8 : bheight / 4 ;
    uint8_t sub_blk_cols = use_dc_pred
        ? (bheight == FORCED_BLK_SIZE && bwidth == FORCED_BLK_SIZE) ? 1 : bwidth / 8
        : use8blk ? bwidth / 8 : bwidth / 4 ;
    for (uint32_t sub_blk_index_y = 0; sub_blk_index_y < sub_blk_rows; ++sub_blk_index_y) {
        for (uint32_t sub_blk_index_x = 0; sub_blk_index_x < sub_blk_cols; ++sub_blk_index_x) {
            TxSize tx_size   = use_dc_pred
                  ? ((bheight == FORCED_BLK_SIZE && bwidth == FORCED_BLK_SIZE) ? TX_16X16 : TX_8X8)
                  : use8blk ? TX_8X8 : TX_4X4 ;
            sub_blk_origin_x = blk_origin_x + sub_blk_index_x * bwidth / sub_blk_cols;
            sub_blk_origin_y = blk_origin_y + sub_blk_index_y * bheight / sub_blk_rows;
            above_row        = above_data + 16;
            left_col         = left_data + 16;

            // Fill Neighbor Arrays
            update_neighbor_samples_array_open_loop_mb(
                                                       0, // use_top_righ_bottom_left
                                                       0, // update_top_neighbor
                                                       above_row - 1,
                                                       left_col - 1,
                                                       input_picture_ptr,
                                                       input_picture_ptr->stride_y,
                                                       sub_blk_origin_x,
                                                       sub_blk_origin_y,
                                                       bwidth / sub_blk_cols,
                                                       bheight / sub_blk_rows);
            // point to  top_neighbor at input buffer
            if (sub_blk_origin_y != 0 ) {
                (above_row)  = ((input_picture_ptr->buffer_y + (((sub_blk_origin_y + input_picture_ptr->origin_y) * input_picture_ptr->stride_y) + (sub_blk_origin_x + input_picture_ptr->origin_x)) ) - input_picture_ptr->stride_y);
            }
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
#if OPT_FP_LOG
    //22025 is (exp(10) -1) equivalent to (log1p((double)this_intra_error) < 10.0)
    if (this_intra_error < 22025)
        stats->intra_factor += 1.5 - log1p((double)this_intra_error) * 0.05;
    else
        stats->intra_factor += 1.0;

    int level_sample = input_picture_ptr->buffer_y[input_origin_index];

    //8102 is (exp(9) -1) equivalent to (log1p((double)this_intra_error) < 9.0)
    if ((level_sample < DARK_THRESH) && (this_intra_error < 8102))
        stats->brightness_factor += 1.0 + (0.01 * (DARK_THRESH - level_sample));
    else
        stats->brightness_factor += 1.0;
#else /*OPT_FP_LOG*/
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
#endif /*OPT_FP_LOG*/
    // Intrapenalty below deals with situations where the intra and inter
    // error scores are very low (e.g. a plain black frame).
    // We do not have special cases in first pass for 0,0 and nearest etc so
    // all inter modes carry an overhead cost estimate for the mv.
    // When the error score is very low this causes us to pick all or lots of
    // INTRA modes and throw lots of key frames.
    // This penalty adds a cost matching that of a 0,0 mv to the intra case.
    this_intra_error += INTRA_MODE_PENALTY;

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
#if OPT_FIRST_PASS4
    FRAME_STATS *stats, int down_step) {
#else
    FRAME_STATS *stats) {
#endif
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
#if !OPT_ME
        uint8_t list_index = 0;
#endif
        uint8_t ref_pic_index = 0;

#if OPT_ME
#if ME_8X8
#if FTR_M13
        if (!ppcs_ptr->enable_me_8x8) {
            if (me_mb_offset >= MAX_SB64_PU_COUNT_NO_8X8)
                me_mb_offset = me_idx_85_8x8_to_16x16_conversion[me_mb_offset - MAX_SB64_PU_COUNT_NO_8X8];
            if (!ppcs_ptr->enable_me_16x16)
                if (me_mb_offset >= MAX_SB64_PU_COUNT_WO_16X16)
                    me_mb_offset = me_idx_16x16_to_parent_32x32_conversion[me_mb_offset - MAX_SB64_PU_COUNT_WO_16X16];
        }
#else
        if (!ppcs_ptr->enable_me_8x8)
            if (me_mb_offset >= MAX_SB64_PU_COUNT_NO_8X8)
                me_mb_offset = me_idx_85_8x8_to_16x16_conversion[me_mb_offset - MAX_SB64_PU_COUNT_NO_8X8];
#endif
#endif
            mv.col = (me_results->me_mv_array[me_mb_offset*ppcs_ptr->pa_me_data->max_refs + ref_pic_index].x_mv) >> 2;
            mv.row = (me_results->me_mv_array[me_mb_offset*ppcs_ptr->pa_me_data->max_refs + ref_pic_index].y_mv) >> 2;

#else
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
#endif
        EbPictureBufferDesc *last_input_picture_ptr = ppcs_ptr->first_pass_ref_count
            ? ppcs_ptr->first_pass_ref_ppcs_ptr[0]->enhanced_picture_ptr
            : NULL;
        int32_t ref_origin_index;
        if (last_input_picture_ptr != NULL)
        {
            ref_origin_index = last_input_picture_ptr->origin_x + (blk_origin_x + mv.col) +
                (blk_origin_y + mv.row + last_input_picture_ptr->origin_y) *
                last_input_picture_ptr->stride_y;
#if OPT_FIRST_PASS4
            motion_error = (uint32_t)(spatial_full_dist_type_fun(input_picture_ptr->buffer_y,
                input_origin_index,
                input_picture_ptr->stride_y << down_step,
                last_input_picture_ptr->buffer_y,
                ref_origin_index,
                last_input_picture_ptr->stride_y << down_step,
                bwidth,
                bheight >> down_step) << down_step);
#else
            motion_error = (uint32_t)(spatial_full_dist_type_fun(input_picture_ptr->buffer_y,
                input_origin_index,
                input_picture_ptr->stride_y,
                last_input_picture_ptr->buffer_y,
                ref_origin_index,
                last_input_picture_ptr->stride_y,
                bwidth,
                bheight));
#endif
        }

        // Assume 0,0 motion with no mv overhead.
        if (mv.col != 0 && mv.row != 0) {
            motion_error += NEW_MV_MODE_PENALTY;
        }
        // Motion search in 2nd reference frame.
        int gf_motion_error = motion_error;
        if (ppcs_ptr->first_pass_ref_count > 1 &&
            me_results->total_me_candidate_index[me_mb_offset] > 1) {
            // To convert full-pel MV
#if !OPT_ME
            list_index    = 0;
#endif
            ref_pic_index = 1;
            FULLPEL_MV gf_mv;
#if OPT_ME

            gf_mv.col = (me_results->me_mv_array[me_mb_offset*ppcs_ptr->pa_me_data->max_refs + ref_pic_index].x_mv) >> 2;
            gf_mv.row = (me_results->me_mv_array[me_mb_offset*ppcs_ptr->pa_me_data->max_refs + ref_pic_index].y_mv) >> 2;

#else
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
#endif
            EbPictureBufferDesc *golden_input_picture_ptr =
                ppcs_ptr->first_pass_ref_ppcs_ptr[1]->enhanced_picture_ptr;
            ref_origin_index = golden_input_picture_ptr->origin_x + (blk_origin_x + gf_mv.col) +
                (blk_origin_y + gf_mv.row + golden_input_picture_ptr->origin_y) *
                    golden_input_picture_ptr->stride_y;
#if OPT_FIRST_PASS4
            gf_motion_error = (uint32_t)(
                spatial_full_dist_type_fun(input_picture_ptr->buffer_y,
                    input_origin_index,
                    input_picture_ptr->stride_y << down_step,
                    golden_input_picture_ptr->buffer_y,
                    ref_origin_index,
                    golden_input_picture_ptr->stride_y << down_step,
                    bwidth,
                    bheight >> down_step) << down_step);
#else
            gf_motion_error = (uint32_t)(
                spatial_full_dist_type_fun(input_picture_ptr->buffer_y,
                                           input_origin_index,
                                           input_picture_ptr->stride_y,
                                           golden_input_picture_ptr->buffer_y,
                                           ref_origin_index,
                                           golden_input_picture_ptr->stride_y,
                                           bwidth,
                                           bheight));
#endif

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

#if !CLN_2PASS
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
#endif
    } else {
        stats->sr_coded_error += motion_error;
#if !CLN_2PASS
        stats->tr_coded_error += motion_error;
#endif
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
#if !CLN_2PASS
        stats->sum_mvr += best_mv.row;
#endif
        stats->sum_mvr_abs += abs(best_mv.row);
        stats->sum_mvc += best_mv.col;
        stats->sum_mvc_abs += abs(best_mv.col);
#if !CLN_2PASS
        stats->sum_mvrs += best_mv.row * best_mv.row;
        stats->sum_mvcs += best_mv.col * best_mv.col;
#endif
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
static EbErrorType first_pass_frame_seg(PictureParentControlSet *ppcs_ptr, int32_t segment_index) {
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

    uint32_t blks_in_b64 = BLOCK_SIZE_64 / FORCED_BLK_SIZE;
    uint32_t picture_width_in_b64  = (uint32_t)(input_picture_ptr->width + BLOCK_SIZE_64 - 1) / BLOCK_SIZE_64;
    uint32_t picture_height_in_b64 = (uint32_t)(input_picture_ptr->height + BLOCK_SIZE_64 - 1) / BLOCK_SIZE_64;

    uint32_t x_seg_idx;
    uint32_t y_seg_idx;
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

    const uint32_t blk_index_y_end = (y_b64_end_idx * blks_in_b64) > blk_rows ? blk_rows : (y_b64_end_idx * blks_in_b64);
    const uint32_t blk_index_x_end = (x_b64_end_idx * blks_in_b64) > blk_cols ? blk_cols : (x_b64_end_idx * blks_in_b64);
    EbSpatialFullDistType spatial_full_dist_type_fun = svt_spatial_full_distortion_kernel;
#if OPT_FIRST_PASS4
    int down_step = ppcs_ptr->scs_ptr->static_config.ipp_ctrls.dist_ds;
#endif
    for (uint32_t blk_index_y = (y_b64_start_idx * blks_in_b64); blk_index_y < blk_index_y_end; blk_index_y++) {
        for (uint32_t blk_index_x = (x_b64_start_idx * blks_in_b64); blk_index_x < blk_index_x_end; blk_index_x++) {
#if OPT_FIRST_PASS3
            int sb_index_x = blk_index_x / blks_in_b64;
            int sb_index_y = blk_index_y / blks_in_b64;
            if (ppcs_ptr->bypass_blk_step > 1)
                if ((sb_index_x % (uint32_t)ppcs_ptr->bypass_blk_step != 0) || (sb_index_y % (uint32_t)ppcs_ptr->bypass_blk_step != 0))
                    continue;
#endif
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

            if (ppcs_ptr->first_pass_ref_count)
                ppcs_ptr->firstpass_data.raw_motion_err_list[blk_index_y * blk_cols + blk_index_x] =
#if OPT_FIRST_PASS4
            (uint32_t)(spatial_full_dist_type_fun(input_picture_ptr->buffer_y,
                input_origin_index,
                input_picture_ptr->stride_y << down_step,
                last_input_picture_ptr->buffer_y,
                input_origin_index,
                input_picture_ptr->stride_y << down_step,
                blk_width,
                blk_height >> down_step) << down_step);
#else
                    (uint32_t)(spatial_full_dist_type_fun(input_picture_ptr->buffer_y,
                                                          input_origin_index,
                                                          input_picture_ptr->stride_y,
                                                          last_input_picture_ptr->buffer_y,
                                                          input_origin_index,
                                                          input_picture_ptr->stride_y,
                                                          blk_width,
                                                          blk_height));
#endif

            int this_intra_error = open_loop_firstpass_intra_prediction(ppcs_ptr, ppcs_ptr->firstpass_data.raw_motion_err_list[blk_index_y * blk_cols + blk_index_x] ,blk_origin_x,
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
#if OPT_FIRST_PASS4
                    mb_stats,
                    down_step);
#else
                    mb_stats);
#endif

                if (blk_origin_x == 0)
                    first_top_mv = last_mv;

                mb_stats->coded_error += this_inter_error;
            } else {
                mb_stats->sr_coded_error += this_intra_error;
#if !CLN_2PASS
                mb_stats->tr_coded_error += this_intra_error;
#endif
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

    // Load the SB from the input to the intermediate SB buffer
    int buffer_index = (input_picture_ptr->origin_y + sb_origin_y) * input_picture_ptr->stride_y +
        input_picture_ptr->origin_x + sb_origin_x;

    // set search method
    context_ptr->me_context_ptr->hme_search_method = SUB_SAD_SEARCH;

#ifdef ARCH_X86_64
    uint8_t *src_ptr = &(input_picture_ptr->buffer_y[buffer_index]);

    uint32_t sb_height = (input_picture_ptr->height - sb_origin_y) < BLOCK_SIZE_64
        ? input_picture_ptr->height - sb_origin_y
        : BLOCK_SIZE_64;
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

    context_ptr->me_context_ptr->quarter_sb_buffer = &quarter_pic_ptr->buffer_y[buffer_index];
    context_ptr->me_context_ptr->quarter_sb_buffer_stride = quarter_pic_ptr->stride_y;

    // Load the 1/16 decimated SB from the 1/16 decimated input to the 1/16 intermediate SB buffer
    buffer_index = (sixteenth_pic_ptr->origin_y + (sb_origin_y >> 2)) *
            sixteenth_pic_ptr->stride_y +
        sixteenth_pic_ptr->origin_x + (sb_origin_x >> 2);

    context_ptr->me_context_ptr->sixteenth_sb_buffer = &sixteenth_pic_ptr->buffer_y[buffer_index];
    context_ptr->me_context_ptr->sixteenth_sb_buffer_stride = sixteenth_pic_ptr->stride_y;
}
/***************************************************************************
* Perform the motion estimation for first pass.
***************************************************************************/
#if OPT_FIRST_PASS3
static EbErrorType first_pass_me(PictureParentControlSet *  ppcs_ptr,
    MotionEstimationContext_t *me_context_ptr, int32_t segment_index) {
    EbPictureBufferDesc *input_picture_ptr = ppcs_ptr->enhanced_picture_ptr;
    uint32_t ss_x = ppcs_ptr->scs_ptr->subsampling_x;
    uint32_t ss_y = ppcs_ptr->scs_ptr->subsampling_y;
    const uint32_t blk_cols = (uint32_t)(input_picture_ptr->width + FORCED_BLK_SIZE - 1) /
        FORCED_BLK_SIZE;
    const uint32_t blk_rows = (uint32_t)(input_picture_ptr->height + FORCED_BLK_SIZE - 1) /
        FORCED_BLK_SIZE;
    uint32_t sb_cols = (uint32_t)(input_picture_ptr->width + BLOCK_SIZE_64 - 1) / BLOCK_SIZE_64;
    uint32_t blks_in_b64 = BLOCK_SIZE_64 / FORCED_BLK_SIZE;
    uint32_t picture_width_in_b64 = (uint32_t)(input_picture_ptr->width + BLOCK_SIZE_64 - 1) / BLOCK_SIZE_64;
    uint32_t picture_height_in_b64 = (uint32_t)(input_picture_ptr->height + BLOCK_SIZE_64 - 1) / BLOCK_SIZE_64;

    uint32_t x_seg_idx;
    uint32_t y_seg_idx;
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

    const uint32_t blk_index_y_end = (y_b64_end_idx * blks_in_b64) > blk_rows ? blk_rows : (y_b64_end_idx * blks_in_b64);
    const uint32_t blk_index_x_end = (x_b64_end_idx * blks_in_b64) > blk_cols ? blk_cols : (x_b64_end_idx * blks_in_b64);
    MeContext *context_ptr = me_context_ptr->me_context_ptr;
#if OPT_ME
    // assume max 2 references for first pass
    ppcs_ptr->pa_me_data->max_cand = 3;
    ppcs_ptr->pa_me_data->max_refs = 2;
    ppcs_ptr->pa_me_data->max_l0 = 1;
#endif
    for (uint32_t blk_index_y = (y_b64_start_idx * blks_in_b64); blk_index_y < blk_index_y_end; blk_index_y+= blks_in_b64) {
        for (uint32_t blk_index_x = (x_b64_start_idx * blks_in_b64); blk_index_x < blk_index_x_end; blk_index_x+= blks_in_b64) {
            int sb_index_x = blk_index_x / blks_in_b64;
            int sb_index_y = blk_index_y / blks_in_b64;
            if (ppcs_ptr->bypass_blk_step > 1)
                if ((sb_index_x % (uint32_t)ppcs_ptr->bypass_blk_step != 0) || (sb_index_y % (uint32_t)ppcs_ptr->bypass_blk_step != 0))
                    continue;

            // Initialize ME context
            first_pass_setup_me_context(
                me_context_ptr, ppcs_ptr, input_picture_ptr, sb_index_y, sb_index_x, ss_x, ss_y);
            // Perform ME - context_ptr will store the outputs (MVs, buffers, etc)
            // Block-based MC using open-loop HME + refinement
            motion_estimate_sb(
                ppcs_ptr, // source picture control set -> references come from here
                (uint32_t)sb_index_y * sb_cols + sb_index_x,
                (uint32_t)sb_index_x * BLOCK_SIZE_64, // x block
                (uint32_t)sb_index_y * BLOCK_SIZE_64, // y block
                context_ptr,
                input_picture_ptr); // source picture
        }
    }
    return EB_ErrorNone;
}
#else
static EbErrorType first_pass_me(PictureParentControlSet *  ppcs_ptr,
                                 MotionEstimationContext_t *me_context_ptr, int32_t segment_index) {
    EbPictureBufferDesc *input_picture_ptr = ppcs_ptr->enhanced_picture_ptr;

    uint32_t blk_cols = (uint32_t)(input_picture_ptr->width + BLOCK_SIZE_64 - 1) / BLOCK_SIZE_64;
    uint32_t blk_rows = (uint32_t)(input_picture_ptr->height + BLOCK_SIZE_64 - 1) / BLOCK_SIZE_64;
    uint32_t ss_x     = ppcs_ptr->scs_ptr->subsampling_x;
    uint32_t ss_y     = ppcs_ptr->scs_ptr->subsampling_y;

    MeContext *context_ptr = me_context_ptr->me_context_ptr;
#if OPT_ME
     // assume max 2 references for first pass
     ppcs_ptr->pa_me_data->max_cand = 3;
     ppcs_ptr->pa_me_data->max_refs = 2;
     ppcs_ptr->pa_me_data->max_l0 = 1;
#endif
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
#if OPT_FIRST_PASS3
    uint32_t step = (uint32_t)ppcs_ptr->bypass_blk_step;
    for (uint32_t blk_row = y_b64_start_idx; blk_row < y_b64_end_idx; blk_row+= step) {
        for (uint32_t blk_col = x_b64_start_idx; blk_col < x_b64_end_idx; blk_col+= step) {
#else
    for (uint32_t blk_row = y_b64_start_idx; blk_row < y_b64_end_idx; blk_row++) {
        for (uint32_t blk_col = x_b64_start_idx; blk_col < x_b64_end_idx; blk_col++) {
#endif
            // Initialize ME context
            first_pass_setup_me_context(
                me_context_ptr, ppcs_ptr, input_picture_ptr, blk_row, blk_col, ss_x, ss_y);
            // Perform ME - context_ptr will store the outputs (MVs, buffers, etc)
            // Block-based MC using open-loop HME + refinement
            motion_estimate_sb(
                                ppcs_ptr, // source picture control set -> references come from here
                               (uint32_t)blk_row * blk_cols + blk_col,
                               (uint32_t)blk_col * BLOCK_SIZE_64, // x block
                               (uint32_t)blk_row * BLOCK_SIZE_64, // y block
                               context_ptr,
                               input_picture_ptr); // source picture
        }
    }
    return EB_ErrorNone;
}
#endif
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
#if IPP_CTRL
    if (!ppcs_ptr->scs_ptr->static_config.ipp_ctrls.skip_frame_first_pass)
#else
    if (ppcs_ptr->scs_ptr->enc_mode_2ndpass <= ENC_M4)
#endif
        ppcs_ptr->skip_frame =0;
    else {
#if OPT_FIRST_PASS
#if ENBLE_SKIP_FRAME_IN_VBR_MODE
#if FIX_ISSUE_50
        if ((ppcs_ptr->scs_ptr->static_config.ipp_ctrls.skip_frame_first_pass == 1) && (ppcs_ptr->picture_number % 8 > 0))
#else
        if ((ppcs_ptr->scs_ptr->static_config.final_pass_rc_mode == 0) && (ppcs_ptr->picture_number % 8 > 0))
#endif
            ppcs_ptr->skip_frame = 1;
        else if ((ppcs_ptr->scs_ptr->static_config.ipp_ctrls.skip_frame_first_pass == 2) && (ppcs_ptr->picture_number > 7 && ppcs_ptr->picture_number % 8 > 0))
#else
        if ((ppcs_ptr->scs_ptr->static_config.final_pass_rc_mode == 0) && (ppcs_ptr->picture_number % 8 > 0))
#endif
            ppcs_ptr->skip_frame = 1;
        else
#endif
        if (ppcs_ptr->picture_number > 3 && ppcs_ptr->picture_number % 4 > 0)
            ppcs_ptr->skip_frame = 1;
        else
            ppcs_ptr->skip_frame = 0;
    }
#if OPT_FIRST_PASS3
    ppcs_ptr->bypass_blk_step = ppcs_ptr->scs_ptr->static_config.ipp_ctrls.bypass_blk_step ? ((ppcs_ptr->scs_ptr->input_resolution < INPUT_SIZE_360p_RANGE) ? 1 : 2) : 1;
#endif
    if (!ppcs_ptr->skip_frame)
    if (ppcs_ptr->first_pass_ref_count)
        first_pass_me(ppcs_ptr, me_context_ptr, segment_index);

    if (!ppcs_ptr->skip_frame){
    setup_firstpass_data_seg(ppcs_ptr, segment_index);
    // Perform the processing of the segment for each frame after me is done for all blocks
    first_pass_frame_seg(ppcs_ptr, segment_index);
    }
    svt_block_on_mutex(ppcs_ptr->first_pass_mutex);
    ppcs_ptr->first_pass_seg_acc++;
    if (ppcs_ptr->first_pass_seg_acc == ppcs_ptr->first_pass_seg_total_count) {

        first_pass_frame_end(ppcs_ptr, ppcs_ptr->ts_duration);
        if (ppcs_ptr->end_of_sequence_flag && !ppcs_ptr->scs_ptr->lap_enabled)
            svt_av1_end_first_pass(ppcs_ptr);
        // Signal that the first pass is done
        svt_post_semaphore(ppcs_ptr->first_pass_done_semaphore);
    }

    svt_release_mutex(ppcs_ptr->first_pass_mutex);
}
