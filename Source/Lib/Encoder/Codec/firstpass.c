// clang-format off
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
static EbErrorType realloc_stats_out(SequenceControlSet *scs, FirstPassStatsOut *out,
                                     uint64_t frame_number) {
    if (frame_number < out->size)
        return EB_ErrorNone;

    if ((int64_t)frame_number >= (int64_t)out->capability - 1) {
        size_t capability = (int64_t)frame_number >= (int64_t)STATS_CAPABILITY_INIT - 1
            ? STATS_CAPABILITY_GROW(frame_number)
            : STATS_CAPABILITY_INIT;
        if (scs->lap_rc) {
            //store the data points before re-allocation
            uint64_t stats_in_start_offset = 0;
            uint64_t stats_in_offset       = 0;
            uint64_t stats_in_end_offset   = 0;
            if (frame_number) {
                stats_in_start_offset = scs->twopass.stats_buf_ctx->stats_in_start - out->stat;
                stats_in_offset       = scs->twopass.stats_in - out->stat;
                stats_in_end_offset   = scs->twopass.stats_buf_ctx->stats_in_end_write -
                    out->stat;
            }
            EB_REALLOC_ARRAY(out->stat, capability);
            // restore the pointers after re-allocation is done
            scs->twopass.stats_buf_ctx->stats_in_start     = out->stat + stats_in_start_offset;
            scs->twopass.stats_in                          = out->stat + stats_in_offset;
            scs->twopass.stats_buf_ctx->stats_in_end_write = out->stat + stats_in_end_offset;
        } else {
            EB_REALLOC_ARRAY(out->stat, capability);
        }
        out->capability = capability;
    }
    out->size = frame_number + 1;
    return EB_ErrorNone;
}

static AOM_INLINE void output_stats(SequenceControlSet *scs, const FIRSTPASS_STATS *stats,
                                    uint64_t frame_number) {
    FirstPassStatsOut *stats_out = &scs->enc_ctx->stats_out;
    svt_block_on_mutex(scs->enc_ctx->stat_file_mutex);
    if (realloc_stats_out(scs, stats_out, frame_number) != EB_ErrorNone) {
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
    svt_release_mutex(scs->enc_ctx->stat_file_mutex);
}
void svt_av1_twopass_zero_stats(FIRSTPASS_STATS *section) {
    section->frame              = 0.0;
    section->weight             = 0.0;
    section->intra_error        = 0.0;
    section->coded_error        = 0.0;
    section->sr_coded_error     = 0.0;
    section->pcnt_inter         = 0.0;
    section->pcnt_motion        = 0.0;
    section->pcnt_second_ref    = 0.0;
    section->pcnt_neutral       = 0.0;
    section->intra_skip_pct     = 0.0;
    section->inactive_zone_rows = 0.0;
    section->inactive_zone_cols = 0.0;
    section->mvr_abs = 0.0;
    section->mvc_abs = 0.0;
    section->mv_in_out_count = 0.0;
    section->count    = 0.0;
    section->duration = 1.0;
    memset(&section->stat_struct, 0, sizeof(StatStruct));
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
    section->mvr_abs += frame->mvr_abs;
    section->mvc_abs += frame->mvc_abs;
    section->mv_in_out_count += frame->mv_in_out_count;
    section->count += frame->count;
    section->duration += frame->duration;
}
void svt_av1_end_first_pass(PictureParentControlSet *pcs) {
    SequenceControlSet *scs = pcs->scs;
    TWO_PASS *          twopass = &scs->twopass;

    if (twopass->stats_buf_ctx->total_stats) {
        // add the total to the end of the file
        svt_block_on_mutex(twopass->stats_buf_ctx->stats_in_write_mutex);

        FIRSTPASS_STATS total_stats = *twopass->stats_buf_ctx->total_stats;
        output_stats(scs, &total_stats, pcs->picture_number + 1);

        svt_release_mutex(twopass->stats_buf_ctx->stats_in_write_mutex);
    }
}
#define UL_INTRA_THRESH 50
#define INVALID_ROW -1
// Accumulates motion vector stats.
// Modifies member variables of "stats".
static void accumulate_mv_stats(const MV best_mv, const FULLPEL_MV mv, const int mb_row, const int mb_col,
                         const int mb_rows, const int mb_cols, MV *last_mv, FRAME_STATS *stats) {
    if (is_zero_mv(&best_mv))
        return;

    ++stats->mv_count;

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
static void update_firstpass_stats(PictureParentControlSet *pcs, const FRAME_STATS *const stats,
                                   const int frame_number,
                                   uint8_t skip_frame, uint8_t bypass_blk_step,
                                   const double ts_duration) {
    SequenceControlSet *scs = pcs->scs;
    TWO_PASS *          twopass = &scs->twopass;

    const uint32_t   mb_cols          = (scs->max_input_luma_width + 16 - 1) / 16;
    const uint32_t   mb_rows          = (scs->max_input_luma_height + 16 - 1) / 16;

    svt_block_on_mutex(twopass->stats_buf_ctx->stats_in_write_mutex);

    FIRSTPASS_STATS *this_frame_stats = twopass->stats_buf_ctx->stats_in_end_write;
    FIRSTPASS_STATS  fps;
    // The minimum error here insures some bit allocation to frames even
    // in static regions. The allocation per MB declines for larger formats
    // where the typical "real" energy per MB also falls.
    // Initial estimate here uses sqrt(mbs) to define the min_err, where the
    // number of mbs is proportional to the image area.
    uint32_t step = (uint32_t)bypass_blk_step;
    const int num_mbs = (mb_rows / step) * (mb_cols / step);
    //(cpi->oxcf.resize_cfg.resize_mode != RESIZE_NONE)
    //    ? cpi->initial_mbs
    //    : mi_params->MBs;
    const double min_err = 200 * sqrt(num_mbs);
    if (skip_frame) {
        FirstPassStatsOut *stats_out = &scs->enc_ctx->stats_out;
        fps       = stats_out->stat[frame_number - 1];
        fps.frame = frame_number;
    } else {
        fps.weight = stats->intra_factor * stats->brightness_factor;
        fps.frame = frame_number;
        fps.coded_error    = (double)(stats->coded_error >> 8) + min_err;
        fps.sr_coded_error = (double)(stats->sr_coded_error >> 8) + min_err;
        fps.intra_error     = (double)(stats->intra_error >> 8) + min_err;
        // if blocks are skipped, the errors need to be updated
        if (bypass_blk_step == 2) {
            fps.coded_error *= 3;
            fps.sr_coded_error *= 3;
            fps.intra_error *= 3;
        }
        fps.count           = 1.0;
        fps.pcnt_inter      = (double)stats->inter_count / num_mbs;
        fps.pcnt_second_ref = (double)stats->second_ref_count / num_mbs;
        fps.pcnt_neutral       = (double)stats->neutral_count / num_mbs;
        fps.intra_skip_pct     = (double)stats->intra_skip_count / num_mbs;
        fps.inactive_zone_rows = (double)stats->image_data_start_row;
        fps.inactive_zone_cols = (double)0; // TODO(paulwilkins): fix

        if (stats->mv_count > 0) {
            fps.mvr_abs = (double)stats->sum_mvr_abs / stats->mv_count;
            fps.mvc_abs = (double)stats->sum_mvc_abs / stats->mv_count;
            fps.mv_in_out_count = (double)stats->sum_in_vectors / (stats->mv_count * 2);
            fps.pcnt_motion = (double)stats->mv_count / num_mbs;
        } else {
            fps.mvr_abs = 0.0;
            fps.mvc_abs = 0.0;
            fps.mv_in_out_count = 0.0;
            fps.pcnt_motion = 0.0;
        }
        // TODO(paulwilkins):  Handle the case when duration is set to 0, or
        // something less than the full time between subsequent values of
        // cpi->source_time_stamp.
        fps.duration = (double)ts_duration;
    }
    memset(&fps.stat_struct, 0, sizeof(StatStruct));
    // We will store the stats inside the persistent twopass struct (and NOT the
    // local variable 'fps'), and then cpi->output_pkt_list will point to it.
    *this_frame_stats = fps;
    output_stats(scs, &fps, pcs->picture_number);
    if (twopass->stats_buf_ctx->total_stats != NULL &&
        scs->static_config.pass == ENC_FIRST_PASS) {
        svt_av1_accumulate_stats(twopass->stats_buf_ctx->total_stats, &fps);
    }
    /*In the case of two pass, first pass uses it as a circular buffer,
   * when LAP is enabled it is used as a linear buffer*/
    twopass->stats_buf_ctx->stats_in_end_write++;
    if (scs->static_config.pass == ENC_FIRST_PASS &&
        (twopass->stats_buf_ctx->stats_in_end_write >= twopass->stats_buf_ctx->stats_in_buf_end)) {
        twopass->stats_buf_ctx->stats_in_end_write = twopass->stats_buf_ctx->stats_in_start;
    }
    svt_release_mutex(twopass->stats_buf_ctx->stats_in_write_mutex);
}

static FRAME_STATS accumulate_frame_stats(FRAME_STATS *mb_stats, int mb_rows, int mb_cols,
                                          uint8_t bypass_blk_step) {
    FRAME_STATS stats = {0};
    int         i, j;
    uint32_t blks_in_b64 = BLOCK_SIZE_64 / FORCED_BLK_SIZE;
    stats.image_data_start_row = INVALID_ROW;
    for (j = 0; j < mb_rows; j++) {
        for (i = 0; i < mb_cols; i++) {
            FRAME_STATS mb_stat = mb_stats[j * mb_cols + i];
            int sb_index_x = i / blks_in_b64;
            int sb_index_y = j / blks_in_b64;
            if (bypass_blk_step > 1)
                if ((sb_index_x % (uint32_t)bypass_blk_step != 0) ||
                    (sb_index_y % (uint32_t)bypass_blk_step != 0))
                    continue;
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
            stats.second_ref_count += mb_stat.second_ref_count;
            stats.sr_coded_error += mb_stat.sr_coded_error;
            stats.sum_in_vectors += mb_stat.sum_in_vectors;
            stats.sum_mvc += mb_stat.sum_mvc;
            stats.sum_mvc_abs += mb_stat.sum_mvc_abs;
        }
    }
    return stats;
}
/**************************************************
 * Reset first pass stat
 **************************************************/
void setup_firstpass_data_seg(PictureParentControlSet *ppcs, int32_t segment_index) {
    SequenceControlSet * scs           = ppcs->scs;
    FirstPassData *      firstpass_data    = &ppcs->firstpass_data;
    const uint32_t       mb_cols           = (scs->max_input_luma_width + 16 - 1) / 16;
    const uint32_t       mb_rows           = (scs->max_input_luma_height + 16 - 1) / 16;
    EbPictureBufferDesc *input_pic = ppcs->enhanced_pic;

    uint32_t blk_cols = (uint32_t)(input_pic->width + BLOCK_SIZE_64 - 1) / BLOCK_SIZE_64;
    uint32_t blk_rows = (uint32_t)(input_pic->height + BLOCK_SIZE_64 - 1) / BLOCK_SIZE_64;

    uint32_t x_seg_idx;
    uint32_t y_seg_idx;
    uint32_t picture_width_in_b64  = blk_cols;
    uint32_t picture_height_in_b64 = blk_rows;
    SEGMENT_CONVERT_IDX_TO_XY(
        segment_index, x_seg_idx, y_seg_idx, ppcs->first_pass_seg_column_count);
    uint32_t x_b64_start_idx = SEGMENT_START_IDX(
        x_seg_idx, picture_width_in_b64, ppcs->first_pass_seg_column_count);
    uint32_t x_b64_end_idx = SEGMENT_END_IDX(
        x_seg_idx, picture_width_in_b64, ppcs->first_pass_seg_column_count);
    uint32_t y_b64_start_idx = SEGMENT_START_IDX(
        y_seg_idx, picture_height_in_b64, ppcs->first_pass_seg_row_count);
    uint32_t y_b64_end_idx = SEGMENT_END_IDX(
        y_seg_idx, picture_height_in_b64, ppcs->first_pass_seg_row_count);

    const uint32_t mb_y_end = (y_b64_end_idx << 2) > mb_rows ? mb_rows : (y_b64_end_idx << 2);
    const uint32_t mb_x_end = (x_b64_end_idx << 2) > mb_cols ? mb_cols : (x_b64_end_idx << 2);

    for (uint32_t mb_y = (y_b64_start_idx << 2); mb_y < mb_y_end; mb_y++) {
        for (uint32_t mb_x = (x_b64_start_idx << 2); mb_x < mb_x_end; mb_x++) {
            memset(firstpass_data->mb_stats + mb_x + mb_y * mb_cols,
                   0,
                   sizeof(*firstpass_data->mb_stats));
            firstpass_data->mb_stats[mb_x + mb_y * mb_cols].image_data_start_row = INVALID_ROW;
        }
    }
}
static void first_pass_frame_end(PictureParentControlSet *pcs, uint8_t skip_frame,
                          uint8_t bypass_blk_step, const double ts_duration) {
    SequenceControlSet *scs = pcs->scs;
    const uint32_t      mb_cols = (scs->max_input_luma_width + 16 - 1) / 16;
    const uint32_t      mb_rows = (scs->max_input_luma_height + 16 - 1) / 16;

    FRAME_STATS *mb_stats = pcs->firstpass_data.mb_stats;

    FRAME_STATS stats;

    if (!skip_frame) {
        stats = accumulate_frame_stats(mb_stats, mb_rows, mb_cols, bypass_blk_step);
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
        uint32_t step = (uint32_t)bypass_blk_step;
        const int num_mbs = (mb_rows / step) * (mb_cols / step);
        /*(cpi->oxcf.resize_cfg.resize_mode != RESIZE_NONE)
        ? cpi->initial_mbs
        : mi_params->MBs;*/
        stats.intra_factor      = stats.intra_factor / (double)num_mbs;
        stats.brightness_factor = stats.brightness_factor / (double)num_mbs;
    } else {
        memset(&stats, 0, sizeof(stats));
    }
    update_firstpass_stats(
        pcs,
        &stats,
        (const int)pcs->picture_number,
        skip_frame,
        bypass_blk_step,
        ts_duration);

}

#define LOW_MOTION_ERROR_THRESH 25
#define MOTION_ERROR_THRESH 500
void set_tf_controls(PictureParentControlSet *pcs, uint8_t tf_level);

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
static int open_loop_firstpass_intra_prediction(
    PictureParentControlSet *ppcs, int raw_motion_err, uint32_t blk_org_x,
    uint32_t blk_org_y, uint8_t bwidth, uint8_t bheight, EbPictureBufferDesc *input_pic,
    uint32_t input_origin_index, FRAME_STATS *const stats) {
    int32_t   mb_row      = blk_org_y >> 4;
    int32_t   mb_col      = blk_org_x >> 4;
    const int use_dc_pred = (mb_col || mb_row) && (!mb_col || !mb_row);
    uint8_t   use8blk     = 0;
    if (!ppcs->scs->ipp_pass_ctrls.use8blk) {
        use8blk = 0;
    } else {
        if (ppcs->first_pass_ref_count)
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
    uint8_t *predictor    = predictor8;
    uint8_t  sub_blk_rows = use_dc_pred
         ? (bheight == FORCED_BLK_SIZE && bwidth == FORCED_BLK_SIZE) ? 1 : bheight / 8
         : use8blk ? bheight / 8
                   : bheight / 4;
    uint8_t  sub_blk_cols = use_dc_pred
         ? (bheight == FORCED_BLK_SIZE && bwidth == FORCED_BLK_SIZE) ? 1 : bwidth / 8
         : use8blk ? bwidth / 8
                   : bwidth / 4;
    for (uint32_t sub_blk_index_y = 0; sub_blk_index_y < sub_blk_rows; ++sub_blk_index_y) {
        for (uint32_t sub_blk_index_x = 0; sub_blk_index_x < sub_blk_cols; ++sub_blk_index_x) {
            TxSize tx_size   = use_dc_pred
                  ? ((bheight == FORCED_BLK_SIZE && bwidth == FORCED_BLK_SIZE) ? TX_16X16 : TX_8X8)
                  : use8blk ? TX_8X8
                            : TX_4X4;
            sub_blk_origin_x = blk_org_x + sub_blk_index_x * bwidth / sub_blk_cols;
            sub_blk_origin_y = blk_org_y + sub_blk_index_y * bheight / sub_blk_rows;
            above_row        = above_data + 16;
            left_col         = left_data + 16;

            // Fill Neighbor Arrays
            svt_aom_update_neighbor_samples_array_open_loop_mb(0, // use_top_righ_bottom_left
                                                       0, // update_top_neighbor
                                                       above_row - 1,
                                                       left_col - 1,
                                                       input_pic,
                                                       input_pic->stride_y,
                                                       sub_blk_origin_x,
                                                       sub_blk_origin_y,
                                                       bwidth / sub_blk_cols,
                                                       bheight / sub_blk_rows);
            // point to  top_neighbor at input buffer
            if (sub_blk_origin_y != 0) {
                (above_row) = ((input_pic->buffer_y +
                                (((sub_blk_origin_y + input_pic->org_y) *
                                  input_pic->stride_y) +
                                 (sub_blk_origin_x + input_pic->org_x))) -
                               input_pic->stride_y);
            }
            // PRED
            predictor = &predictor8[(sub_blk_origin_x - blk_org_x) +
                                    (sub_blk_origin_y - blk_org_y) * FORCED_BLK_SIZE];
            svt_aom_intra_prediction_open_loop_mb(0,
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
    int this_intra_error = (uint32_t)(spatial_full_dist_type_fun(input_pic->buffer_y,
                                                                 input_origin_index,
                                                                 input_pic->stride_y,
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
    //22025 is (exp(10) -1) equivalent to (log1p((double)this_intra_error) < 10.0)
    if (this_intra_error < 22025)
        stats->intra_factor += 1.5 - log1p((double)this_intra_error) * 0.05;
    else
        stats->intra_factor += 1.0;

    int level_sample = input_pic->buffer_y[input_origin_index];

    //8102 is (exp(9) -1) equivalent to (log1p((double)this_intra_error) < 9.0)
    if ((level_sample < DARK_THRESH) && (this_intra_error < 8102))
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
    PictureParentControlSet *ppcs, uint32_t me_sb_addr, uint32_t blk_org_x,
    uint32_t blk_org_y, uint8_t bwidth, uint8_t bheight, EbPictureBufferDesc *input_pic,
    uint32_t input_origin_index, const int this_intra_error, MV *last_mv, int raw_motion_err,
    FRAME_STATS *stats, int down_step) {
    int32_t        mb_row  = blk_org_y >> 4;
    int32_t        mb_col  = blk_org_x >> 4;
    const uint32_t mb_cols = (ppcs->scs->max_input_luma_width + FORCED_BLK_SIZE - 1) /
        FORCED_BLK_SIZE;
    const uint32_t mb_rows = (ppcs->scs->max_input_luma_height + FORCED_BLK_SIZE -
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
        const MeSbResults *me_results = ppcs->pa_me_data->me_results[me_sb_addr];
        uint32_t           me_sb_size = ppcs->scs->b64_size;
        blk_geom.org_x             = blk_org_x - (blk_org_x / me_sb_size) * me_sb_size;
        blk_geom.org_y             = blk_org_y - (blk_org_y / me_sb_size) * me_sb_size;
        blk_geom.bwidth               = FORCED_BLK_SIZE;
        blk_geom.bheight              = FORCED_BLK_SIZE;
        me_mb_offset = svt_aom_get_me_info_index(ppcs->max_number_of_pus_per_sb, &blk_geom, 0, 0);
        uint8_t ref_pic_index = 0;

        if (!ppcs->enable_me_8x8) {
            if (me_mb_offset >= MAX_SB64_PU_COUNT_NO_8X8)
                me_mb_offset =
                    me_idx_85_8x8_to_16x16_conversion[me_mb_offset - MAX_SB64_PU_COUNT_NO_8X8];
            if (!ppcs->enable_me_16x16)
                if (me_mb_offset >= MAX_SB64_PU_COUNT_WO_16X16)
                    me_mb_offset =
                        me_idx_16x16_to_parent_32x32_conversion[me_mb_offset -
                                                                MAX_SB64_PU_COUNT_WO_16X16];
        }
        mv.col = (me_results
                      ->me_mv_array[me_mb_offset * ppcs->pa_me_data->max_refs + ref_pic_index]
                      .x_mv) >>
            2;
        mv.row = (me_results
                      ->me_mv_array[me_mb_offset * ppcs->pa_me_data->max_refs + ref_pic_index]
                      .y_mv) >>
            2;

        EbPictureBufferDesc *last_input_picture_ptr = ppcs->first_pass_ref_count
            ? ppcs->first_pass_ref_ppcs_ptr[0]->enhanced_pic
            : NULL;
        int32_t              ref_origin_index;
        if (last_input_picture_ptr != NULL) {
            ref_origin_index = (int32_t)last_input_picture_ptr->org_x +
                ((int32_t)blk_org_x + mv.col) +
                ((int32_t)blk_org_y + mv.row + (int32_t)last_input_picture_ptr->org_y) *
                    (int32_t)last_input_picture_ptr->stride_y;
            motion_error = (uint32_t)(
                spatial_full_dist_type_fun(input_pic->buffer_y,
                                           input_origin_index,
                                           input_pic->stride_y << down_step,
                                           last_input_picture_ptr->buffer_y,
                                           ref_origin_index,
                                           last_input_picture_ptr->stride_y << down_step,
                                           bwidth,
                                           bheight >> down_step)
                << down_step);
        }

        // Assume 0,0 motion with no mv overhead.
        if (mv.col != 0 && mv.row != 0) {
            motion_error += NEW_MV_MODE_PENALTY;
        }
        // Motion search in 2nd reference frame.
        int gf_motion_error = motion_error;
        if (ppcs->first_pass_ref_count > 1 &&
            me_results->total_me_candidate_index[me_mb_offset] > 1) {
            // To convert full-pel MV
            ref_pic_index = 1;
            FULLPEL_MV gf_mv;

            gf_mv.col =
                (me_results
                     ->me_mv_array[me_mb_offset * ppcs->pa_me_data->max_refs + ref_pic_index]
                     .x_mv) >>
                2;
            gf_mv.row =
                (me_results
                     ->me_mv_array[me_mb_offset * ppcs->pa_me_data->max_refs + ref_pic_index]
                     .y_mv) >>
                2;

            EbPictureBufferDesc *golden_input_picture_ptr =
                ppcs->first_pass_ref_ppcs_ptr[1]->enhanced_pic;
            ref_origin_index = (int32_t)golden_input_picture_ptr->org_x +
                ((int32_t)blk_org_x + gf_mv.col) +
                ((int32_t)blk_org_y + gf_mv.row + (int32_t)golden_input_picture_ptr->org_y) *
                    (int32_t)golden_input_picture_ptr->stride_y;
            gf_motion_error = (uint32_t)(
                spatial_full_dist_type_fun(input_pic->buffer_y,
                                           input_origin_index,
                                           input_pic->stride_y << down_step,
                                           golden_input_picture_ptr->buffer_y,
                                           ref_origin_index,
                                           golden_input_picture_ptr->stride_y << down_step,
                                           bwidth,
                                           bheight >> down_step)
                << down_step);

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
        if (ppcs->first_pass_ref_count > 1 && (gf_motion_error < motion_error * 3)) {
            stats->sr_coded_error += AOMMIN(gf_motion_error, this_intra_error);
        } else {
            stats->sr_coded_error += motion_error;
        }

    } else {
        stats->sr_coded_error += motion_error;
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
        stats->sum_mvr_abs += abs(best_mv.row);
        stats->sum_mvc += best_mv.col;
        stats->sum_mvc_abs += abs(best_mv.col);
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
static EbErrorType first_pass_frame_seg(PictureParentControlSet *ppcs, int32_t segment_index,
                                        uint8_t bypass_blk_step) {
    EbPictureBufferDesc *input_pic      = ppcs->enhanced_pic;
    EbPictureBufferDesc *last_input_picture_ptr = ppcs->first_pass_ref_count
        ? ppcs->first_pass_ref_ppcs_ptr[0]->enhanced_pic
        : NULL;

    const uint32_t blk_cols = (uint32_t)(input_pic->width + FORCED_BLK_SIZE - 1) /
        FORCED_BLK_SIZE;
    const uint32_t blk_rows = (uint32_t)(input_pic->height + FORCED_BLK_SIZE - 1) /
        FORCED_BLK_SIZE;

    uint32_t me_sb_size         = ppcs->scs->b64_size;
    uint32_t me_pic_width_in_sb = (ppcs->aligned_width + me_sb_size - 1) / me_sb_size;
    uint32_t me_sb_x, me_sb_y, me_sb_addr;

    uint32_t blk_width, blk_height, blk_org_x, blk_org_y;
    MV       first_top_mv = kZeroMv;
    MV       last_mv;
    uint32_t input_origin_index;

    uint32_t blks_in_b64          = BLOCK_SIZE_64 / FORCED_BLK_SIZE;
    uint32_t picture_width_in_b64 = (uint32_t)(input_pic->width + BLOCK_SIZE_64 - 1) /
        BLOCK_SIZE_64;
    uint32_t picture_height_in_b64 = (uint32_t)(input_pic->height + BLOCK_SIZE_64 - 1) /
        BLOCK_SIZE_64;

    uint32_t x_seg_idx;
    uint32_t y_seg_idx;
    SEGMENT_CONVERT_IDX_TO_XY(
        segment_index, x_seg_idx, y_seg_idx, ppcs->first_pass_seg_column_count);
    uint32_t x_b64_start_idx = SEGMENT_START_IDX(
        x_seg_idx, picture_width_in_b64, ppcs->first_pass_seg_column_count);
    uint32_t x_b64_end_idx = SEGMENT_END_IDX(
        x_seg_idx, picture_width_in_b64, ppcs->first_pass_seg_column_count);
    uint32_t y_b64_start_idx = SEGMENT_START_IDX(
        y_seg_idx, picture_height_in_b64, ppcs->first_pass_seg_row_count);
    uint32_t y_b64_end_idx = SEGMENT_END_IDX(
        y_seg_idx, picture_height_in_b64, ppcs->first_pass_seg_row_count);

    const uint32_t        blk_index_y_end            = (y_b64_end_idx * blks_in_b64) > blk_rows
                          ? blk_rows
                          : (y_b64_end_idx * blks_in_b64);
    const uint32_t        blk_index_x_end            = (x_b64_end_idx * blks_in_b64) > blk_cols
                          ? blk_cols
                          : (x_b64_end_idx * blks_in_b64);
    EbSpatialFullDistType spatial_full_dist_type_fun = svt_spatial_full_distortion_kernel;
    int down_step = ppcs->scs->ipp_pass_ctrls.dist_ds;
    for (uint32_t blk_index_y = (y_b64_start_idx * blks_in_b64); blk_index_y < blk_index_y_end;
         blk_index_y++) {
        for (uint32_t blk_index_x = (x_b64_start_idx * blks_in_b64); blk_index_x < blk_index_x_end;
             blk_index_x++) {
            int sb_index_x = blk_index_x / blks_in_b64;
            int sb_index_y = blk_index_y / blks_in_b64;
            if (bypass_blk_step > 1)
                if ((sb_index_x % (uint32_t)bypass_blk_step != 0) ||
                    (sb_index_y % (uint32_t)bypass_blk_step != 0))
                    continue;
            blk_org_x = blk_index_x * FORCED_BLK_SIZE;
            blk_org_y = blk_index_y * FORCED_BLK_SIZE;
            me_sb_x      = blk_org_x / me_sb_size;
            me_sb_y      = blk_org_y / me_sb_size;
            me_sb_addr   = me_sb_x + me_sb_y * me_pic_width_in_sb;

            blk_width  = (ppcs->aligned_width - blk_org_x) < FORCED_BLK_SIZE
                 ? ppcs->aligned_width - blk_org_x
                 : FORCED_BLK_SIZE;
            blk_height = (ppcs->aligned_height - blk_org_y) < FORCED_BLK_SIZE
                ? ppcs->aligned_height - blk_org_y
                : FORCED_BLK_SIZE;

            input_origin_index = (input_pic->org_y + blk_org_y) *
                    input_pic->stride_y +
                (input_pic->org_x + blk_org_x);

            FRAME_STATS *mb_stats = ppcs->firstpass_data.mb_stats + blk_index_y * blk_cols +
                blk_index_x;

            if (ppcs->first_pass_ref_count)
                ppcs->firstpass_data.raw_motion_err_list[blk_index_y * blk_cols + blk_index_x] =
                    (uint32_t)(spatial_full_dist_type_fun(input_pic->buffer_y,
                                                          input_origin_index,
                                                          input_pic->stride_y << down_step,
                                                          last_input_picture_ptr->buffer_y,
                                                          input_origin_index,
                                                          input_pic->stride_y << down_step,
                                                          blk_width,
                                                          blk_height >> down_step)
                               << down_step);

            int this_intra_error = open_loop_firstpass_intra_prediction(
                ppcs,
                ppcs->firstpass_data.raw_motion_err_list[blk_index_y * blk_cols + blk_index_x],
                blk_org_x,
                blk_org_y,
                blk_width,
                blk_height,
                input_pic,
                input_origin_index,
                mb_stats);

            int this_inter_error = this_intra_error;

            if (blk_org_x == 0)
                last_mv = first_top_mv;

            if (ppcs->first_pass_ref_count) {
                this_inter_error = open_loop_firstpass_inter_prediction(
                    ppcs,
                    me_sb_addr,
                    blk_org_x,
                    blk_org_y,
                    blk_width,
                    blk_height,
                    input_pic,
                    input_origin_index,
                    this_intra_error,
                    &last_mv,
                    ppcs->firstpass_data
                        .raw_motion_err_list[blk_index_y * blk_cols + blk_index_x],
                    mb_stats,
                    down_step);

                if (blk_org_x == 0)
                    first_top_mv = last_mv;

                mb_stats->coded_error += this_inter_error;
            } else {
                mb_stats->sr_coded_error += this_intra_error;
                mb_stats->coded_error += this_intra_error;
            }
        }
    }

    return EB_ErrorNone;
}
/***************************************************************************
* Prepare the me context for performing first pass me.
***************************************************************************/
static void first_pass_setup_me_context(MotionEstimationContext_t *me_context_ptr,
                                        PictureParentControlSet *  ppcs,
                                        EbPictureBufferDesc *input_pic, int blk_row,
                                        int blk_col, uint32_t ss_x, uint32_t ss_y) {
    // setup the references
    me_context_ptr->me_ctx->num_of_list_to_search = 1;
    me_context_ptr->me_ctx->num_of_ref_pic_to_search[0] = 0;
    me_context_ptr->me_ctx->num_of_ref_pic_to_search[1] = 0;
    me_context_ptr->me_ctx->temporal_layer_index        = 0;
    me_context_ptr->me_ctx->is_ref   = 1;

    if (ppcs->first_pass_ref_count) {
        me_context_ptr->me_ctx->me_ds_ref_array[0][0] =
            ppcs->first_pass_ref_ppcs_ptr[0]->ds_pics;
        me_context_ptr->me_ctx->num_of_ref_pic_to_search[0]++;
    }
    if (ppcs->first_pass_ref_count > 1) {
        me_context_ptr->me_ctx->me_ds_ref_array[0][1] =
            ppcs->first_pass_ref_ppcs_ptr[1]->ds_pics;
        me_context_ptr->me_ctx->num_of_ref_pic_to_search[0]++;
    }

    me_context_ptr->me_ctx->me_type = ME_FIRST_PASS;
    // Set 1/4 and 1/16 ME reference buffer(s); filtered or decimated
    EbPictureBufferDesc *quarter_pic_ptr = ppcs->ds_pics.quarter_picture_ptr;

    EbPictureBufferDesc *sixteenth_pic_ptr = ppcs->ds_pics.sixteenth_picture_ptr;
    // Parts from MotionEstimationKernel()
    uint32_t b64_origin_x = (uint32_t)(blk_col * BLOCK_SIZE_64);
    uint32_t b64_origin_y = (uint32_t)(blk_row * BLOCK_SIZE_64);

    // Load the SB from the input to the intermediate SB buffer
    int buffer_index = (input_pic->org_y + b64_origin_y) * input_pic->stride_y +
        input_pic->org_x + b64_origin_x;

    // set search method
    me_context_ptr->me_ctx->hme_search_method = SUB_SAD_SEARCH;

#ifdef ARCH_X86_64
    uint8_t *src_ptr = &(input_pic->buffer_y[buffer_index]);

    uint32_t b64_height = (input_pic->height - b64_origin_y) < BLOCK_SIZE_64
        ? input_pic->height - b64_origin_y
        : BLOCK_SIZE_64;
    //_MM_HINT_T0     //_MM_HINT_T1    //_MM_HINT_T2    //_MM_HINT_NTA
    uint32_t i;
    for (i = 0; i < b64_height; i++) {
        char const *p = (char const *)(src_ptr + i * input_pic->stride_y);
        _mm_prefetch(p, _MM_HINT_T2);
    }
#endif
    me_context_ptr->me_ctx->b64_src_ptr    = &(input_pic->buffer_y[buffer_index]);
    me_context_ptr->me_ctx->b64_src_stride = input_pic->stride_y;

    // Load the 1/4 decimated SB from the 1/4 decimated input to the 1/4 intermediate SB buffer
    buffer_index = (quarter_pic_ptr->org_y + (b64_origin_y >> ss_y)) * quarter_pic_ptr->stride_y +
        quarter_pic_ptr->org_x + (b64_origin_x >> ss_x);

    me_context_ptr->me_ctx->quarter_b64_buffer = &quarter_pic_ptr->buffer_y[buffer_index];
    me_context_ptr->me_ctx->quarter_b64_buffer_stride = quarter_pic_ptr->stride_y;

    // Load the 1/16 decimated SB from the 1/16 decimated input to the 1/16 intermediate SB buffer
    buffer_index = (sixteenth_pic_ptr->org_y + (b64_origin_y >> 2)) * sixteenth_pic_ptr->stride_y +
        sixteenth_pic_ptr->org_x + (b64_origin_x >> 2);

    me_context_ptr->me_ctx->sixteenth_b64_buffer = &sixteenth_pic_ptr->buffer_y[buffer_index];
    me_context_ptr->me_ctx->sixteenth_b64_buffer_stride = sixteenth_pic_ptr->stride_y;
}
/***************************************************************************
* Perform the motion estimation for first pass.
***************************************************************************/
static EbErrorType first_pass_me(PictureParentControlSet *  ppcs,
                                 MotionEstimationContext_t *me_context_ptr, int32_t segment_index) {
    EbPictureBufferDesc *input_pic = ppcs->enhanced_pic;
    uint32_t             ss_x              = ppcs->scs->subsampling_x;
    uint32_t             ss_y              = ppcs->scs->subsampling_y;
    const uint32_t       blk_cols = (uint32_t)(input_pic->width + FORCED_BLK_SIZE - 1) /
        FORCED_BLK_SIZE;
    const uint32_t blk_rows = (uint32_t)(input_pic->height + FORCED_BLK_SIZE - 1) /
        FORCED_BLK_SIZE;
    uint32_t sb_cols     = (uint32_t)(input_pic->width + BLOCK_SIZE_64 - 1) / BLOCK_SIZE_64;
    uint32_t blks_in_b64 = BLOCK_SIZE_64 / FORCED_BLK_SIZE;
    uint32_t picture_width_in_b64 = (uint32_t)(input_pic->width + BLOCK_SIZE_64 - 1) /
        BLOCK_SIZE_64;
    uint32_t picture_height_in_b64 = (uint32_t)(input_pic->height + BLOCK_SIZE_64 - 1) /
        BLOCK_SIZE_64;

    uint32_t x_seg_idx;
    uint32_t y_seg_idx;
    SEGMENT_CONVERT_IDX_TO_XY(
        segment_index, x_seg_idx, y_seg_idx, ppcs->first_pass_seg_column_count);
    uint32_t x_b64_start_idx = SEGMENT_START_IDX(
        x_seg_idx, picture_width_in_b64, ppcs->first_pass_seg_column_count);
    uint32_t x_b64_end_idx = SEGMENT_END_IDX(
        x_seg_idx, picture_width_in_b64, ppcs->first_pass_seg_column_count);
    uint32_t y_b64_start_idx = SEGMENT_START_IDX(
        y_seg_idx, picture_height_in_b64, ppcs->first_pass_seg_row_count);
    uint32_t y_b64_end_idx = SEGMENT_END_IDX(
        y_seg_idx, picture_height_in_b64, ppcs->first_pass_seg_row_count);

    const uint32_t blk_index_y_end = (y_b64_end_idx * blks_in_b64) > blk_rows
        ? blk_rows
        : (y_b64_end_idx * blks_in_b64);
    const uint32_t blk_index_x_end = (x_b64_end_idx * blks_in_b64) > blk_cols
        ? blk_cols
        : (x_b64_end_idx * blks_in_b64);
    MeContext *    me_ctx     = me_context_ptr->me_ctx;
    // assume max 2 references for first pass
    ppcs->pa_me_data->max_cand = 3;
    ppcs->pa_me_data->max_refs = 2;
    ppcs->pa_me_data->max_l0   = 1;
    for (uint32_t blk_index_y = (y_b64_start_idx * blks_in_b64); blk_index_y < blk_index_y_end; blk_index_y += blks_in_b64) {
        for (uint32_t blk_index_x = (x_b64_start_idx * blks_in_b64); blk_index_x < blk_index_x_end; blk_index_x += blks_in_b64) {
            int b64_index_x = blk_index_x / blks_in_b64;
            int b64_index_y = blk_index_y / blks_in_b64;

            uint8_t bypass_blk_step = me_context_ptr->me_ctx->bypass_blk_step;
            if (bypass_blk_step > 1)
                if ((b64_index_x % (uint32_t)bypass_blk_step != 0) || (b64_index_y % (uint32_t)bypass_blk_step != 0))
                    continue;

            // Initialize ME context
            first_pass_setup_me_context(me_context_ptr, ppcs, input_pic, b64_index_y, b64_index_x, ss_x, ss_y);
            // Perform ME - me_ctx will store the outputs (MVs, buffers, etc)
            // Block-based MC using open-loop HME + refinement

            svt_aom_motion_estimation_b64(ppcs, // source picture control set -> references come from here
                (uint32_t)b64_index_y * sb_cols + b64_index_x,
                (uint32_t)b64_index_x * BLOCK_SIZE_64, // x block
                (uint32_t)b64_index_y * BLOCK_SIZE_64, // y block
                me_ctx,
                input_pic); // source picture
        }
    }
    return EB_ErrorNone;
}
/************************************************************************************
* Performs the first pass based on open loop data.
* Source frames are used for Intra and Inter prediction.
* ME is done per segment but the remaining parts performed per frame.
************************************************************************************/
void open_loop_first_pass(PictureParentControlSet *  ppcs,
                          MotionEstimationContext_t *me_context_ptr, int32_t segment_index) {
    me_context_ptr->me_ctx->tf_mv_dist_th = CLIP3(64, 450, (int)((int) MIN(ppcs->aligned_height,ppcs->aligned_width) - 150));
    // Perform the me for the first pass for each segment

    me_context_ptr->me_ctx->bypass_blk_step =
        ppcs->scs->ipp_pass_ctrls.bypass_blk_step
        ?
        ((ppcs->scs->input_resolution < INPUT_SIZE_360p_RANGE) ? 1 : 2)
        : 1;
    if (!me_context_ptr->me_ctx->skip_frame)
        if (ppcs->first_pass_ref_count)
            first_pass_me(ppcs, me_context_ptr, segment_index);
    if (!me_context_ptr->me_ctx->skip_frame) {
        setup_firstpass_data_seg(ppcs, segment_index);
        // Perform the processing of the segment for each frame after me is done for all blocks
        first_pass_frame_seg(
            ppcs, segment_index, me_context_ptr->me_ctx->bypass_blk_step);
    }
    svt_block_on_mutex(ppcs->first_pass_mutex);
    ppcs->first_pass_seg_acc++;
    if (ppcs->first_pass_seg_acc == ppcs->first_pass_seg_total_count) {
        first_pass_frame_end(ppcs,
                             me_context_ptr->me_ctx->skip_frame,
                             me_context_ptr->me_ctx->bypass_blk_step,
                             ppcs->ts_duration);
        if (ppcs->end_of_sequence_flag && !ppcs->scs->lap_rc)
            svt_av1_end_first_pass(ppcs);
        // Signal that the first pass is done
        svt_post_semaphore(ppcs->first_pass_done_semaphore);
    }

    svt_release_mutex(ppcs->first_pass_mutex);
}
// clang-format on
