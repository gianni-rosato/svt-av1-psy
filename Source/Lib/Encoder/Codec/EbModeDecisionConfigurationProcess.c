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

#include <stdlib.h>

#include "EbEncHandle.h"
#include "EbUtility.h"
#include "EbPictureControlSet.h"
#include "EbModeDecisionConfigurationProcess.h"
#include "EbRateControlResults.h"
#include "EbEncDecTasks.h"
#include "EbReferenceObject.h"
#include "EbModeDecisionProcess.h"
#include "av1me.h"
#include "EbQMatrices.h"
#include "EbLog.h"
#include "EbCoefficients.h"
#include "EbCommonUtils.h"
#include "EbResize.h"
#include "EbInvTransforms.h"
#include "EncModeConfig.h"
#include "EbGlobalMotionEstimation.h"
#include "aom_dsp_rtcd.h"
#define MAX_MESH_SPEED 5 // Max speed setting for mesh motion method
static MeshPattern good_quality_mesh_patterns[MAX_MESH_SPEED + 1][MAX_MESH_STEP] = {
    {{64, 8}, {28, 4}, {15, 1}, {7, 1}},
    {{64, 8}, {28, 4}, {15, 1}, {7, 1}},
    {{64, 8}, {14, 2}, {7, 1}, {7, 1}},
    {{64, 16}, {24, 8}, {12, 4}, {7, 1}},
    {{64, 16}, {24, 8}, {12, 4}, {7, 1}},
    {{64, 16}, {24, 8}, {12, 4}, {7, 1}},
};
// TODO: These settings are pretty relaxed, tune them for
// each speed setting
static MeshPattern intrabc_mesh_patterns[MAX_MESH_SPEED + 1][MAX_MESH_STEP] = {
    {{256, 1}, {256, 1}, {0, 0}, {0, 0}},
    {{256, 1}, {256, 1}, {0, 0}, {0, 0}},
    {{64, 1}, {64, 1}, {0, 0}, {0, 0}},
    {{64, 1}, {64, 1}, {0, 0}, {0, 0}},
    {{64, 4}, {16, 1}, {0, 0}, {0, 0}},
    {{64, 4}, {16, 1}, {0, 0}, {0, 0}},
};
void set_global_motion_field(PictureControlSet *pcs) {
    // Init Global Motion Vector
    uint8_t frame_index;
    for (frame_index = INTRA_FRAME; frame_index <= ALTREF_FRAME; ++frame_index) {
        pcs->ppcs->global_motion[frame_index].wmtype   = IDENTITY;
        pcs->ppcs->global_motion[frame_index].alpha    = 0;
        pcs->ppcs->global_motion[frame_index].beta     = 0;
        pcs->ppcs->global_motion[frame_index].delta    = 0;
        pcs->ppcs->global_motion[frame_index].gamma    = 0;
        pcs->ppcs->global_motion[frame_index].invalid  = 0;
        pcs->ppcs->global_motion[frame_index].wmmat[0] = 0;
        pcs->ppcs->global_motion[frame_index].wmmat[1] = 0;
        pcs->ppcs->global_motion[frame_index].wmmat[2] = (1 << WARPEDMODEL_PREC_BITS);
        pcs->ppcs->global_motion[frame_index].wmmat[3] = 0;
        pcs->ppcs->global_motion[frame_index].wmmat[4] = 0;
        pcs->ppcs->global_motion[frame_index].wmmat[5] = (1 << WARPEDMODEL_PREC_BITS);
        pcs->ppcs->global_motion[frame_index].wmmat[6] = 0;
        pcs->ppcs->global_motion[frame_index].wmmat[7] = 0;
    }

    //Update MV
    PictureParentControlSet *ppcs = pcs->ppcs;
    for (frame_index = INTRA_FRAME; frame_index <= ALTREF_FRAME; ++frame_index) {
        if (ppcs->is_global_motion[get_list_idx(frame_index)][get_ref_frame_idx(frame_index)])
            ppcs->global_motion[frame_index] =
                ppcs->svt_aom_global_motion_estimation[get_list_idx(frame_index)][get_ref_frame_idx(frame_index)];
        uint8_t sf = ppcs->gm_downsample_level == GM_DOWN ? 2 : ppcs->gm_downsample_level == GM_DOWN16 ? 4 : 1;
        svt_aom_upscale_wm_params(&ppcs->global_motion[frame_index], sf);
    }
}

void svt_av1_build_quantizer(EbBitDepth bit_depth, int32_t y_dc_delta_q, int32_t u_dc_delta_q, int32_t u_ac_delta_q,
                             int32_t v_dc_delta_q, int32_t v_ac_delta_q, Quants *const quants, Dequants *const deq, PictureParentControlSet *pcs) {
    int32_t i, q, quant_qtx;

    for (q = 0; q < QINDEX_RANGE; q++) {
        int32_t qzbin_factor     = svt_aom_get_qzbin_factor(q, bit_depth);
        int32_t qrounding_factor = q == 0 ? 64 : 48;
        int32_t qrounding_factor_fp = pcs->scs->static_config.tune != 3 ? 64 : 48;
        int diff = q - pcs->frm_hdr.quantization_params.base_q_idx; // q-range diff based on current quantizer

        for (i = 0; i < 2; ++i) {
            quant_qtx                   = i == 0 ? svt_aom_dc_quant_qtx(q, y_dc_delta_q, bit_depth)
                                                 : svt_aom_ac_quant_qtx(q, 0, bit_depth);
            if (pcs->scs->static_config.sharpness != 0) {
                if (pcs->scs->static_config.sharpness > 0) { // If sharpness is positive
                    if (diff < 0) { // If we're going to lower quant, decrease zbin factor and increase rounding
                        qzbin_factor -= MAX(pcs->scs->static_config.sharpness << 1, abs(diff)); // Example range with --sharpness 4: Subtract qzbin_factor by diff, or by (4 << 1), preferring which is closer to 0
                        qrounding_factor += MAX(pcs->scs->static_config.sharpness << 1, abs(diff));
                        qrounding_factor_fp += MAX(pcs->scs->static_config.sharpness << 1, abs(diff));
                    } else if (diff > 0) { // If we're going higher quant, reset both zbin and rounding factor
                        qzbin_factor     = svt_aom_get_qzbin_factor(q, bit_depth);
                        qrounding_factor = q == 0 ? 64 : 48;
                        qrounding_factor_fp = pcs->scs->static_config.tune != 3 ? 64 : 48;
                    }
                } else if (pcs->scs->static_config.sharpness < 0) { // If sharpness is negative
                    if (diff < 0) { // If we're going to lower quant, reset both zbin and rounding factor
                        qzbin_factor     = svt_aom_get_qzbin_factor(q, bit_depth);
                        qrounding_factor = q == 0 ? 64 : 48;
                        qrounding_factor_fp = pcs->scs->static_config.tune != 3 ? 64 : 48;
                    } else if (diff > 0) { // If we're going higher quant, increase zbin and decrease rounding factors
                        qzbin_factor += MIN(abs(pcs->scs->static_config.sharpness) << 1, diff); // Choose minimum here since diff > 0
                        qrounding_factor -= MIN(abs(pcs->scs->static_config.sharpness) << 1, diff);
                        qrounding_factor_fp -= MIN(abs(pcs->scs->static_config.sharpness) << 1, diff);
                    }
                }
                qzbin_factor = MIN(MAX(qzbin_factor, 1), 256); // Ensure we don't go too low or high
                qrounding_factor = MIN(MAX(qrounding_factor, 1), 256);
                qrounding_factor_fp = MIN(MAX(qrounding_factor_fp, 1), 256);
            }
            svt_aom_invert_quant(&quants->y_quant[q][i], &quants->y_quant_shift[q][i], quant_qtx);
            quants->y_quant_fp[q][i] = (int16_t)((1 << 16) / quant_qtx);
            quants->y_round_fp[q][i] = (int16_t)((qrounding_factor_fp * quant_qtx) >> 7);
            quants->y_zbin[q][i]     = (int16_t)ROUND_POWER_OF_TWO(qzbin_factor * quant_qtx, 7);
            quants->y_round[q][i]    = (int16_t)((qrounding_factor * quant_qtx) >> 7);
            deq->y_dequant_qtx[q][i] = (int16_t)quant_qtx;
            quant_qtx                = i == 0 ? svt_aom_dc_quant_qtx(q, u_dc_delta_q, bit_depth)
                                              : svt_aom_ac_quant_qtx(q, u_ac_delta_q, bit_depth);
            svt_aom_invert_quant(&quants->u_quant[q][i], &quants->u_quant_shift[q][i], quant_qtx);
            quants->u_quant_fp[q][i] = (int16_t)((1 << 16) / quant_qtx);
            quants->u_round_fp[q][i] = (int16_t)((qrounding_factor_fp * quant_qtx) >> 7);
            quants->u_zbin[q][i]     = (int16_t)ROUND_POWER_OF_TWO(qzbin_factor * quant_qtx, 7);
            quants->u_round[q][i]    = (int16_t)((qrounding_factor * quant_qtx) >> 7);
            deq->u_dequant_qtx[q][i] = (int16_t)quant_qtx;
            quant_qtx                = i == 0 ? svt_aom_dc_quant_qtx(q, v_dc_delta_q, bit_depth)
                                              : svt_aom_ac_quant_qtx(q, v_ac_delta_q, bit_depth);
            svt_aom_invert_quant(&quants->v_quant[q][i], &quants->v_quant_shift[q][i], quant_qtx);
            quants->v_quant_fp[q][i] = (int16_t)((1 << 16) / quant_qtx);
            quants->v_round_fp[q][i] = (int16_t)((qrounding_factor_fp * quant_qtx) >> 7);
            quants->v_zbin[q][i]     = (int16_t)ROUND_POWER_OF_TWO(qzbin_factor * quant_qtx, 7);
            quants->v_round[q][i]    = (int16_t)((qrounding_factor * quant_qtx) >> 7);
            deq->v_dequant_qtx[q][i] = (int16_t)quant_qtx;
        }

        for (i = 2; i < 8; i++) { // 8: SIMD width
            quants->y_quant[q][i]       = quants->y_quant[q][1];
            quants->y_quant_fp[q][i]    = quants->y_quant_fp[q][1];
            quants->y_round_fp[q][i]    = quants->y_round_fp[q][1];
            quants->y_quant_shift[q][i] = quants->y_quant_shift[q][1];
            quants->y_zbin[q][i]        = quants->y_zbin[q][1];
            quants->y_round[q][i]       = quants->y_round[q][1];
            deq->y_dequant_qtx[q][i]    = deq->y_dequant_qtx[q][1];

            quants->u_quant[q][i]       = quants->u_quant[q][1];
            quants->u_quant_fp[q][i]    = quants->u_quant_fp[q][1];
            quants->u_round_fp[q][i]    = quants->u_round_fp[q][1];
            quants->u_quant_shift[q][i] = quants->u_quant_shift[q][1];
            quants->u_zbin[q][i]        = quants->u_zbin[q][1];
            quants->u_round[q][i]       = quants->u_round[q][1];
            deq->u_dequant_qtx[q][i]    = deq->u_dequant_qtx[q][1];
            quants->v_quant[q][i]       = quants->u_quant[q][1];
            quants->v_quant_fp[q][i]    = quants->v_quant_fp[q][1];
            quants->v_round_fp[q][i]    = quants->v_round_fp[q][1];
            quants->v_quant_shift[q][i] = quants->v_quant_shift[q][1];
            quants->v_zbin[q][i]        = quants->v_zbin[q][1];
            quants->v_round[q][i]       = quants->v_round[q][1];
            deq->v_dequant_qtx[q][i]    = deq->v_dequant_qtx[q][1];
        }
    }
}

// Reduce the large number of quantizers to a smaller number of levels for which
// different matrices may be defined
static INLINE int aom_get_qmlevel(int qindex, int first, int last) {
    // mapping qindex(0, 255) to QM level(first, last)
    return first + (qindex * (last + 1 - first)) / QINDEX_RANGE;
}

void svt_av1_qm_init(PictureParentControlSet *pcs) {
    const uint8_t num_planes = 3; // MAX_MB_PLANE;// NM- No monochroma
    uint8_t       q, c, t;
    int32_t       current;
    for (q = 0; q < NUM_QM_LEVELS; ++q) {
        for (c = 0; c < num_planes; ++c) {
            current = 0;
            for (t = 0; t < TX_SIZES_ALL; ++t) {
                const int32_t size       = tx_size_2d[t];
                const TxSize  qm_tx_size = av1_get_adjusted_tx_size(t);
                if (q == NUM_QM_LEVELS - 1) {
                    pcs->gqmatrix[q][c][t]  = NULL;
                    pcs->giqmatrix[q][c][t] = NULL;
                } else if (t != qm_tx_size) { // Reuse matrices for 'qm_tx_size'
                    pcs->gqmatrix[q][c][t]  = pcs->gqmatrix[q][c][qm_tx_size];
                    pcs->giqmatrix[q][c][t] = pcs->giqmatrix[q][c][qm_tx_size];
                } else {
                    assert(current + size <= QM_TOTAL_SIZE);
                    pcs->gqmatrix[q][c][t]  = &wt_matrix_ref[q][c >= 1][current];
                    pcs->giqmatrix[q][c][t] = &iwt_matrix_ref[q][c >= 1][current];
                    current += size;
                }
            }
        }
    }

    if (pcs->frm_hdr.quantization_params.using_qmatrix) {
        const int32_t min_qmlevel = pcs->scs->static_config.min_qm_level;
        const int32_t max_qmlevel = pcs->scs->static_config.max_qm_level;
        const int32_t base_qindex = pcs->frm_hdr.quantization_params.base_q_idx;

        pcs->frm_hdr.quantization_params.qm[AOM_PLANE_Y] = aom_get_qmlevel(base_qindex, min_qmlevel, max_qmlevel);
        pcs->frm_hdr.quantization_params.qm[AOM_PLANE_U] = aom_get_qmlevel(
            base_qindex + pcs->frm_hdr.quantization_params.delta_q_ac[AOM_PLANE_U], min_qmlevel, max_qmlevel);
        pcs->frm_hdr.quantization_params.qm[AOM_PLANE_V] = aom_get_qmlevel(
            base_qindex + pcs->frm_hdr.quantization_params.delta_q_ac[AOM_PLANE_V], min_qmlevel, max_qmlevel);
#if DEBUG_QM_LEVEL
        SVT_LOG("\n[svt_av1_qm_init] Frame %d - qindex %d, qmlevel %d %d %d\n",
                (int)pcs->picture_number,
                base_qindex,
                pcs->frm_hdr.quantization_params.qm[AOM_PLANE_Y],
                pcs->frm_hdr.quantization_params.qm[AOM_PLANE_U],
                pcs->frm_hdr.quantization_params.qm[AOM_PLANE_V]);
#endif
    }
}

/******************************************************
* Set the reference sg ep for a given picture
******************************************************/
void set_reference_sg_ep(PictureControlSet *pcs) {
    Av1Common         *cm = pcs->ppcs->av1_cm;
    EbReferenceObject *ref_obj_l0, *ref_obj_l1;
    memset(cm->sg_frame_ep_cnt, 0, SGRPROJ_PARAMS * sizeof(int32_t));
    cm->sg_frame_ep = 0;

    // NADER: set cm->sg_ref_frame_ep[0] = cm->sg_ref_frame_ep[1] = -1 to perform all iterations
    switch (pcs->slice_type) {
    case I_SLICE:
        cm->sg_ref_frame_ep[0] = -1;
        cm->sg_ref_frame_ep[1] = -1;
        break;
    case B_SLICE:
        ref_obj_l0             = (EbReferenceObject *)pcs->ref_pic_ptr_array[REF_LIST_0][0]->object_ptr;
        ref_obj_l1             = (EbReferenceObject *)pcs->ref_pic_ptr_array[REF_LIST_1][0]->object_ptr;
        cm->sg_ref_frame_ep[0] = ref_obj_l0->sg_frame_ep;
        cm->sg_ref_frame_ep[1] = ref_obj_l1->sg_frame_ep;
        break;
    case P_SLICE:
        ref_obj_l0             = (EbReferenceObject *)pcs->ref_pic_ptr_array[REF_LIST_0][0]->object_ptr;
        cm->sg_ref_frame_ep[0] = ref_obj_l0->sg_frame_ep;
        cm->sg_ref_frame_ep[1] = 0;
        break;
    default: SVT_LOG("SG: Not supported picture type"); break;
    }
}

void mode_decision_configuration_init_qp_update(PictureControlSet *pcs) {
    FrameHeader *frm_hdr  = &pcs->ppcs->frm_hdr;
    pcs->intra_coded_area = 0;
    pcs->skip_coded_area  = 0;
    pcs->hp_coded_area    = 0;
    // Init block selection
    // Set reference sg ep
    set_reference_sg_ep(pcs);
    set_global_motion_field(pcs);

    svt_av1_qm_init(pcs->ppcs);
    MdRateEstimationContext *md_rate_est_ctx;

    md_rate_est_ctx = pcs->md_rate_est_ctx;

    if (pcs->ppcs->frm_hdr.primary_ref_frame != PRIMARY_REF_NONE)
        memcpy(&pcs->md_frame_context,
               &pcs->ref_frame_context[pcs->ppcs->frm_hdr.primary_ref_frame],
               sizeof(FRAME_CONTEXT));
    else {
        svt_av1_default_coef_probs(&pcs->md_frame_context, frm_hdr->quantization_params.base_q_idx);
        svt_aom_init_mode_probs(&pcs->md_frame_context);
    }
    // Initial Rate Estimation of the syntax elements
    svt_aom_estimate_syntax_rate(md_rate_est_ctx,
                                 pcs->slice_type == I_SLICE ? TRUE : FALSE,
                                 pcs->ppcs->scs->seq_header.filter_intra_level,
                                 pcs->ppcs->frm_hdr.allow_screen_content_tools,
                                 pcs->ppcs->enable_restoration,
                                 pcs->ppcs->frm_hdr.allow_intrabc,
                                 &pcs->md_frame_context);
    // Initial Rate Estimation of the Motion vectors
    svt_aom_estimate_mv_rate(pcs, md_rate_est_ctx, &pcs->md_frame_context);
    // Initial Rate Estimation of the quantized coefficients
    svt_aom_estimate_coefficients_rate(md_rate_est_ctx, &pcs->md_frame_context);
}

/******************************************************
 * Compute Tc, and Beta offsets for a given picture
 ******************************************************/

static void mode_decision_configuration_context_dctor(EbPtr p) {
    EbThreadContext                  *thread_ctx = (EbThreadContext *)p;
    ModeDecisionConfigurationContext *obj        = (ModeDecisionConfigurationContext *)thread_ctx->priv;

    EB_FREE_ARRAY(obj);
}
/******************************************************
 * Mode Decision Configuration Context Constructor
 ******************************************************/
EbErrorType svt_aom_mode_decision_configuration_context_ctor(EbThreadContext   *thread_ctx,
                                                             const EbEncHandle *enc_handle_ptr, int input_index,
                                                             int output_index) {
    ModeDecisionConfigurationContext *context_ptr;
    EB_CALLOC_ARRAY(context_ptr, 1);
    thread_ctx->priv  = context_ptr;
    thread_ctx->dctor = mode_decision_configuration_context_dctor;

    // Input/Output System Resource Manager FIFOs
    context_ptr->rate_control_input_fifo_ptr = svt_system_resource_get_consumer_fifo(
        enc_handle_ptr->rate_control_results_resource_ptr, input_index);
    context_ptr->mode_decision_configuration_output_fifo_ptr = svt_system_resource_get_producer_fifo(
        enc_handle_ptr->enc_dec_tasks_resource_ptr, output_index);
    return EB_ErrorNone;
}

/******************************************************
* Derive Mode Decision Config Settings for first pass
Input   : encoder mode and tune
Output  : EncDec Kernel signal(s)
******************************************************/
static INLINE int get_relative_dist(const OrderHintInfo *oh, int a, int b) {
    if (!oh->enable_order_hint)
        return 0;

    const int bits = oh->order_hint_bits;

    assert(bits >= 1);
    assert(a >= 0 && a < (1 << bits));
    assert(b >= 0 && b < (1 << bits));

    int       diff = a - b;
    const int m    = 1 << (bits - 1);
    diff           = (diff & (m - 1)) - (diff & m);
    return diff;
}

static int get_block_position(Av1Common *cm, int *mi_r, int *mi_c, int blk_row, int blk_col, MV mv, int sign_bias) {
    const int base_blk_row = (blk_row >> 3) << 3;
    const int base_blk_col = (blk_col >> 3) << 3;

    const int row_offset = (mv.row >= 0) ? (mv.row >> (4 + MI_SIZE_LOG2)) : -((-mv.row) >> (4 + MI_SIZE_LOG2));

    const int col_offset = (mv.col >= 0) ? (mv.col >> (4 + MI_SIZE_LOG2)) : -((-mv.col) >> (4 + MI_SIZE_LOG2));

    const int row = (sign_bias == 1) ? blk_row - row_offset : blk_row + row_offset;
    const int col = (sign_bias == 1) ? blk_col - col_offset : blk_col + col_offset;

    if (row < 0 || row >= (cm->mi_rows >> 1) || col < 0 || col >= (cm->mi_cols >> 1))
        return 0;

    if (row < base_blk_row - (MAX_OFFSET_HEIGHT >> 3) || row >= base_blk_row + 8 + (MAX_OFFSET_HEIGHT >> 3) ||
        col < base_blk_col - (MAX_OFFSET_WIDTH >> 3) || col >= base_blk_col + 8 + (MAX_OFFSET_WIDTH >> 3))
        return 0;

    *mi_r = row;
    *mi_c = col;

    return 1;
}

#define MFMV_STACK_SIZE 3

// Note: motion_filed_projection finds motion vectors of current frame's
// reference frame, and projects them to current frame. To make it clear,
// let's call current frame's reference frame as start frame.
// Call Start frame's reference frames as reference frames.
// Call ref_offset as frame distances between start frame and its reference
// frames.
static int motion_field_projection(Av1Common *cm, PictureControlSet *pcs, MvReferenceFrame start_frame, int dir) {
    TPL_MV_REF *tpl_mvs_base           = pcs->tpl_mvs;
    int         ref_offset[REF_FRAMES] = {0};

    uint8_t list_idx0, ref_idx_l0;
    list_idx0                          = get_list_idx(start_frame);
    ref_idx_l0                         = get_ref_frame_idx(start_frame);
    EbReferenceObject *start_frame_buf = (EbReferenceObject *)pcs->ref_pic_ptr_array[list_idx0][ref_idx_l0]->object_ptr;

    if (start_frame_buf == NULL)
        return 0;

    if (start_frame_buf->frame_type == KEY_FRAME || start_frame_buf->frame_type == INTRA_ONLY_FRAME)
        return 0;

    // MFMV is not applied when the reference picture is of a different spatial resolution
    // (described in the AV1 spec section 7.9.2.)
    if (start_frame_buf->mi_rows != cm->mi_rows || start_frame_buf->mi_cols != cm->mi_cols) {
        return 0;
    }

    const int                 start_frame_order_hint        = start_frame_buf->order_hint;
    const unsigned int *const ref_order_hints               = &start_frame_buf->ref_order_hint[0];
    int                       start_to_current_frame_offset = get_relative_dist(
        &pcs->ppcs->scs->seq_header.order_hint_info, start_frame_order_hint, pcs->ppcs->cur_order_hint);

    for (int i = LAST_FRAME; i <= INTER_REFS_PER_FRAME; ++i)
        ref_offset[i] = get_relative_dist(
            &pcs->ppcs->scs->seq_header.order_hint_info, start_frame_order_hint, ref_order_hints[i - LAST_FRAME]);

    if (dir == 2)
        start_to_current_frame_offset = -start_to_current_frame_offset;

    const MV_REF *const mv_ref_base = start_frame_buf->mvs;
    const int           mvs_rows    = (cm->mi_rows + 1) >> 1;
    const int           mvs_cols    = (cm->mi_cols + 1) >> 1;

    for (int blk_row = 0; blk_row < mvs_rows; ++blk_row) {
        for (int blk_col = 0; blk_col < mvs_cols; ++blk_col) {
            const MV_REF *const mv_ref = &mv_ref_base[blk_row * mvs_cols + blk_col];
            MV                  fwd_mv = mv_ref->mv.as_mv;

            if (mv_ref->ref_frame > INTRA_FRAME) {
                MV        this_mv;
                int       mi_r, mi_c;
                const int ref_frame_offset = ref_offset[mv_ref->ref_frame];

                int pos_valid = abs(ref_frame_offset) <= MAX_FRAME_DISTANCE && ref_frame_offset > 0 &&
                    abs(start_to_current_frame_offset) <= MAX_FRAME_DISTANCE;

                if (pos_valid) {
                    get_mv_projection(&this_mv, fwd_mv, start_to_current_frame_offset, ref_frame_offset);
                    pos_valid = get_block_position(cm, &mi_r, &mi_c, blk_row, blk_col, this_mv, dir >> 1);
                }

                if (pos_valid) {
                    const int mi_offset = mi_r * (cm->mi_stride >> 1) + mi_c;

                    tpl_mvs_base[mi_offset].mfmv0.as_mv.row  = fwd_mv.row;
                    tpl_mvs_base[mi_offset].mfmv0.as_mv.col  = fwd_mv.col;
                    tpl_mvs_base[mi_offset].ref_frame_offset = ref_frame_offset;
                }
            }
        }
    }

    return 1;
}
static void av1_setup_motion_field(Av1Common *cm, PictureControlSet *pcs) {
    const OrderHintInfo *const order_hint_info = &pcs->ppcs->scs->seq_header.order_hint_info;
    memset(pcs->ref_frame_side, 0, sizeof(pcs->ref_frame_side));
    if (!order_hint_info->enable_order_hint)
        return;

    TPL_MV_REF *tpl_mvs_base = pcs->tpl_mvs;
    int         size         = ((cm->mi_rows + MAX_MIB_SIZE) >> 1) * (cm->mi_stride >> 1);

    const int                cur_order_hint = pcs->ppcs->cur_order_hint;
    const EbReferenceObject *ref_buf[INTER_REFS_PER_FRAME];
    int                      ref_order_hint[INTER_REFS_PER_FRAME];

    for (int ref_frame = LAST_FRAME; ref_frame <= ALTREF_FRAME; ref_frame++) {
        const int ref_idx    = ref_frame - LAST_FRAME;
        int       order_hint = 0;
        uint8_t   list_idx0, ref_idx_l0;
        list_idx0              = get_list_idx(ref_frame);
        ref_idx_l0             = get_ref_frame_idx(ref_frame);
        EbReferenceObject *buf = (EbReferenceObject *)pcs->ref_pic_ptr_array[list_idx0][ref_idx_l0]->object_ptr;

        if (buf != NULL)
            order_hint = buf->order_hint;

        ref_buf[ref_idx]        = buf;
        ref_order_hint[ref_idx] = order_hint;

        if (get_relative_dist(order_hint_info, order_hint, cur_order_hint) > 0)
            pcs->ref_frame_side[ref_frame] = 1;
        else if (order_hint == cur_order_hint)
            pcs->ref_frame_side[ref_frame] = -1;
    }

    //for a frame based mfmv, we need to keep computing the ref_frame_side regardless mfmv is used or no
    if (!pcs->ppcs->frm_hdr.use_ref_frame_mvs)
        return;

    for (int idx = 0; idx < size; ++idx) {
        tpl_mvs_base[idx].mfmv0.as_int     = INVALID_MV;
        tpl_mvs_base[idx].ref_frame_offset = 0;
    }

    int ref_stamp = MFMV_STACK_SIZE - 1;
    if (ref_buf[0 /*LAST_FRAME - LAST_FRAME*/] != NULL) {
        const int alt_of_lst_order_hint =
            ref_buf[0 /*LAST_FRAME - LAST_FRAME*/]->ref_order_hint[ALTREF_FRAME - LAST_FRAME];
        const int is_lst_overlay = (alt_of_lst_order_hint == ref_order_hint[GOLDEN_FRAME - LAST_FRAME]);
        if (!is_lst_overlay)
            motion_field_projection(cm, pcs, LAST_FRAME, 2);

        --ref_stamp;
    }

    if (get_relative_dist(order_hint_info, ref_order_hint[BWDREF_FRAME - LAST_FRAME], cur_order_hint) > 0) {
        if (motion_field_projection(cm, pcs, BWDREF_FRAME, 0))
            --ref_stamp;
    }

    if (get_relative_dist(order_hint_info, ref_order_hint[ALTREF2_FRAME - LAST_FRAME], cur_order_hint) > 0) {
        if (motion_field_projection(cm, pcs, ALTREF2_FRAME, 0))
            --ref_stamp;
    }

    if (get_relative_dist(order_hint_info, ref_order_hint[ALTREF_FRAME - LAST_FRAME], cur_order_hint) > 0 &&
        ref_stamp >= 0)
        if (motion_field_projection(cm, pcs, ALTREF_FRAME, 0))
            --ref_stamp;

    if (ref_stamp >= 0)
        motion_field_projection(cm, pcs, LAST2_FRAME, 2);
}
EbErrorType svt_av1_hash_table_create(HashTable *p_hash_table);
void       *rtime_alloc_block_hash_block_is_same(size_t size) { return malloc(size); }
int32_t     svt_aom_noise_log1p_fp16(int32_t noise_level_fp16);
/* Determine the frame complexity level (stored under pcs->coeff_lvl) based
on the ME distortion and QP. */
static void set_frame_coeff_lvl(PictureControlSet *pcs) {
    // Derive the input nois level
    EbPictureBufferDesc *input_pic = pcs->ppcs->enhanced_pic;

    EbByte buffer_y = input_pic->buffer_y + input_pic->org_y * input_pic->stride_y + input_pic->org_x;

    int32_t noise_level_fp16 = svt_estimate_noise_fp16(buffer_y, // Y
                                                       input_pic->width,
                                                       input_pic->height,
                                                       input_pic->stride_y);

    noise_level_fp16         = svt_aom_noise_log1p_fp16(noise_level_fp16);
    uint64_t tot_me_8x8_dist = 0;
    for (uint32_t b64_idx = 0; b64_idx < pcs->b64_total_count; b64_idx++) {
        tot_me_8x8_dist += pcs->ppcs->me_8x8_distortion[b64_idx];
    }
    uint64_t me_8x8_dist_per_sb  = tot_me_8x8_dist / pcs->b64_total_count;
    uint64_t cmplx               = me_8x8_dist_per_sb / MAX(1, pcs->scs->static_config.qp);
    uint64_t coeff_low_level_th  = COEFF_LVL_TH_0;
    uint64_t coeff_high_level_th = COEFF_LVL_TH_1;
    if (pcs->ppcs->input_resolution <= INPUT_SIZE_240p_RANGE) {
        coeff_low_level_th  = (uint64_t)((double)coeff_low_level_th * 1.7);
        coeff_high_level_th = (uint64_t)((double)coeff_high_level_th * 1.7);
    } else if (pcs->ppcs->input_resolution <= INPUT_SIZE_480p_RANGE) {
        coeff_low_level_th  = (uint64_t)((double)coeff_low_level_th * 1.3);
        coeff_high_level_th = (uint64_t)((double)coeff_high_level_th * 1.3);
    } else if (pcs->ppcs->input_resolution <= INPUT_SIZE_720p_RANGE) {
        coeff_low_level_th  = (uint64_t)((double)coeff_low_level_th * 1.2);
        coeff_high_level_th = (uint64_t)((double)coeff_high_level_th * 1.2);
    }

    if (noise_level_fp16 < 26572 /*FLOAT2FP(log1p(0.5), 16, int32_t)*/) {
        coeff_low_level_th  = (uint64_t)((double)coeff_low_level_th * 0.7);
        coeff_high_level_th = (uint64_t)((double)coeff_high_level_th * 0.7);
    } else if (noise_level_fp16 > 45426 /*FLOAT2FP(log1p(1.0), 16, int32_t)*/) {
        coeff_low_level_th  = (uint64_t)((double)coeff_low_level_th * 1.05);
        coeff_high_level_th = (uint64_t)((double)coeff_high_level_th * 1.05);
    }
    pcs->coeff_lvl = NORMAL_LVL;
    if (cmplx < coeff_low_level_th) {
        pcs->coeff_lvl = LOW_LVL;
    } else if (cmplx > coeff_high_level_th) {
        pcs->coeff_lvl = HIGH_LVL;
    }
}

/* Mode Decision Configuration Kernel */

/*********************************************************************************
 *
 * @brief
 *  The Mode Decision Configuration Process involves a number of initialization steps,
 *  setting flags for a number of features, and determining the blocks to be considered
 *  in subsequent MD stages.
 *
 * @par Description:
 *  The Mode Decision Configuration Process involves a number of initialization steps,
 *  setting flags for a number of features, and determining the blocks to be considered
 *  in subsequent MD stages. Examples of flags that are set are the flags for filter intra,
 *  eighth-pel, OBMC and warped motion and flags for updating the cumulative density functions
 *  Examples of initializations include initializations for picture chroma QP offsets,
 *  CDEF strength, self-guided restoration filter parameters, quantization parameters,
 *  lambda arrays, mv and coefficient rate estimation arrays.
 *
 *  The set of blocks to be processed in subsequent MD stages is decided in this process as a
 *  function of the picture depth mode (pic_depth_mode).
 *
 * @param[in] Configurations
 *  Configuration flags that are to be set
 *
 * @param[out] Initializations
 *  Initializations for various flags and variables
 *
 ********************************************************************************/
void *svt_aom_mode_decision_configuration_kernel(void *input_ptr) {
    // Context & SCS & PCS
    EbThreadContext                  *thread_ctx  = (EbThreadContext *)input_ptr;
    ModeDecisionConfigurationContext *context_ptr = (ModeDecisionConfigurationContext *)thread_ctx->priv;
    // Input
    EbObjectWrapper *rc_results_wrapper;

    // Output
    EbObjectWrapper *enc_dec_tasks_wrapper;

    for (;;) {
        // Get RateControl Results
        EB_GET_FULL_OBJECT(context_ptr->rate_control_input_fifo_ptr, &rc_results_wrapper);

        RateControlResults *rc_results = (RateControlResults *)rc_results_wrapper->object_ptr;
        PictureControlSet  *pcs        = (PictureControlSet *)rc_results->pcs_wrapper->object_ptr;
        SequenceControlSet *scs        = pcs->scs;
        pcs->min_me_clpx               = 0;
        pcs->max_me_clpx               = 0;
        pcs->avg_me_clpx               = 0;
        if (pcs->slice_type != I_SLICE) {
            uint32_t b64_idx;
            uint64_t avg_me_clpx = 0;
            uint64_t min_dist    = (uint64_t)~0;
            uint64_t max_dist    = 0;

            for (b64_idx = 0; b64_idx < pcs->ppcs->b64_total_count; ++b64_idx) {
                avg_me_clpx += pcs->ppcs->me_8x8_cost_variance[b64_idx];
                min_dist = MIN(pcs->ppcs->me_8x8_cost_variance[b64_idx], min_dist);
                max_dist = MAX(pcs->ppcs->me_8x8_cost_variance[b64_idx], max_dist);
            }
            pcs->min_me_clpx = min_dist;
            pcs->max_me_clpx = max_dist;
            pcs->avg_me_clpx = avg_me_clpx / pcs->ppcs->b64_total_count;
        }
        pcs->coeff_lvl = INVALID_LVL;

        if (pcs->slice_type != I_SLICE && !pcs->ppcs->sc_class1) {
            set_frame_coeff_lvl(pcs);
        }
        // -------
        // Scale references if resolution of the reference is different than the input
        // super-res reference frame size is same as original input size, only check current frame scaled flag;
        // reference scaling resizes reference frame to different size, need check each reference frame for scaling
        // -------
        if ((pcs->ppcs->frame_superres_enabled == 1 || scs->static_config.resize_mode != RESIZE_NONE) &&
            pcs->slice_type != I_SLICE) {
            if (pcs->ppcs->is_ref == TRUE && pcs->ppcs->ref_pic_wrapper != NULL) {
                // update mi_rows and mi_cols for the reference pic wrapper (used in mfmv for other
                // pictures)
                EbReferenceObject *ref_object = pcs->ppcs->ref_pic_wrapper->object_ptr;
                ref_object->mi_rows           = pcs->ppcs->aligned_height >> MI_SIZE_LOG2;
                ref_object->mi_cols           = pcs->ppcs->aligned_width >> MI_SIZE_LOG2;
            }

            svt_aom_scale_rec_references(pcs, pcs->ppcs->enhanced_pic, pcs->hbd_md);
        }

        FrameHeader *frm_hdr = &pcs->ppcs->frm_hdr;
        pcs->rtc_tune        = (scs->static_config.pred_structure == SVT_AV1_PRED_LOW_DELAY_B) ? true : false;
        // Mode Decision Configuration Kernel Signal(s) derivation
        svt_aom_sig_deriv_mode_decision_config(scs, pcs);

        if (pcs->slice_type != I_SLICE && scs->mfmv_enabled)
            av1_setup_motion_field(pcs->ppcs->av1_cm, pcs);

        pcs->intra_coded_area = 0;
        pcs->skip_coded_area  = 0;
        pcs->hp_coded_area    = 0;
        // Init block selection
        // Set reference sg ep
        set_reference_sg_ep(pcs);
        if (pcs->ppcs->gm_ctrls.use_ref_info) {
            assert(pcs->slice_type != I_SLICE);
            svt_aom_global_motion_estimation(pcs->ppcs, pcs->ppcs->enhanced_pic);
        }
        set_global_motion_field(pcs);

        svt_av1_qm_init(pcs->ppcs);
        MdRateEstimationContext *md_rate_est_ctx;

        // QP
        context_ptr->qp = pcs->picture_qp;

        // QP Index
        context_ptr->qp_index = (uint8_t)frm_hdr->quantization_params.base_q_idx;

        md_rate_est_ctx = pcs->md_rate_est_ctx;
        if (pcs->ppcs->frm_hdr.primary_ref_frame != PRIMARY_REF_NONE)
            memcpy(&pcs->md_frame_context,
                   &pcs->ref_frame_context[pcs->ppcs->frm_hdr.primary_ref_frame],
                   sizeof(FRAME_CONTEXT));
        else {
            svt_av1_default_coef_probs(&pcs->md_frame_context, frm_hdr->quantization_params.base_q_idx);
            svt_aom_init_mode_probs(&pcs->md_frame_context);
        }
        // Initial Rate Estimation of the syntax elements
        svt_aom_estimate_syntax_rate(md_rate_est_ctx,
                                     pcs->slice_type == I_SLICE ? TRUE : FALSE,
                                     scs->seq_header.filter_intra_level,
                                     pcs->ppcs->frm_hdr.allow_screen_content_tools,
                                     pcs->ppcs->enable_restoration,
                                     pcs->ppcs->frm_hdr.allow_intrabc,
                                     &pcs->md_frame_context);
        // Initial Rate Estimation of the Motion vectors
        svt_aom_estimate_mv_rate(pcs, md_rate_est_ctx, &pcs->md_frame_context);
        // Initial Rate Estimation of the quantized coefficients
        svt_aom_estimate_coefficients_rate(md_rate_est_ctx, &pcs->md_frame_context);
        if (frm_hdr->allow_intrabc) {
            int            i;
            int            speed = 1;
            SpeedFeatures *sf    = &pcs->sf;

            const int mesh_speed           = AOMMIN(speed, MAX_MESH_SPEED);
            sf->exhaustive_searches_thresh = (1 << 25);
            if (mesh_speed > 0)
                sf->exhaustive_searches_thresh = sf->exhaustive_searches_thresh << 1;

            for (i = 0; i < MAX_MESH_STEP; ++i) {
                sf->mesh_patterns[i].range    = good_quality_mesh_patterns[mesh_speed][i].range;
                sf->mesh_patterns[i].interval = good_quality_mesh_patterns[mesh_speed][i].interval;
            }

            if (pcs->slice_type == I_SLICE) {
                for (i = 0; i < MAX_MESH_STEP; ++i) {
                    sf->mesh_patterns[i].range    = intrabc_mesh_patterns[mesh_speed][i].range;
                    sf->mesh_patterns[i].interval = intrabc_mesh_patterns[mesh_speed][i].interval;
                }
            }

            {
                // add to hash table
                const int pic_width  = pcs->ppcs->aligned_width;
                const int pic_height = pcs->ppcs->aligned_height;

                uint32_t *block_hash_values[2][2];
                int8_t   *is_block_same[2][3];
                int       k, j;

                for (k = 0; k < 2; k++) {
                    for (j = 0; j < 2; j++)
                        block_hash_values[k][j] = rtime_alloc_block_hash_block_is_same(sizeof(uint32_t) * pic_width *
                                                                                       pic_height);
                    for (j = 0; j < 3; j++)
                        is_block_same[k][j] = rtime_alloc_block_hash_block_is_same(sizeof(int8_t) * pic_width *
                                                                                   pic_height);
                }
                svt_aom_rtime_alloc_svt_av1_hash_table_create(&pcs->hash_table);
                Yv12BufferConfig cpi_source;
                svt_aom_link_eb_to_aom_buffer_desc_8bit(pcs->ppcs->enhanced_pic, &cpi_source);

                svt_av1_crc_calculator_init(&pcs->crc_calculator1, 24, 0x5D6DCB);
                svt_av1_crc_calculator_init(&pcs->crc_calculator2, 24, 0x864CFB);

                svt_av1_generate_block_2x2_hash_value(&cpi_source, block_hash_values[0], is_block_same[0], pcs);
                uint8_t       src_idx     = 0;
                const uint8_t max_sb_size = pcs->ppcs->intraBC_ctrls.max_block_size_hash;
                for (int size = 4; size <= max_sb_size; size <<= 1, src_idx = !src_idx) {
                    const uint8_t dst_idx = !src_idx;
                    svt_av1_generate_block_hash_value(&cpi_source,
                                                      size,
                                                      block_hash_values[src_idx],
                                                      block_hash_values[dst_idx],
                                                      is_block_same[src_idx],
                                                      is_block_same[dst_idx],
                                                      pcs);
                    if (size != 4 || pcs->ppcs->intraBC_ctrls.hash_4x4_blocks)
                        svt_aom_rtime_alloc_svt_av1_add_to_hash_map_by_row_with_precal_data(&pcs->hash_table,
                                                                                            block_hash_values[dst_idx],
                                                                                            is_block_same[dst_idx][2],
                                                                                            pic_width,
                                                                                            pic_height,
                                                                                            size);
                }
                for (k = 0; k < 2; k++) {
                    for (j = 0; j < 2; j++) free(block_hash_values[k][j]);
                    for (j = 0; j < 3; j++) free(is_block_same[k][j]);
                }
            }

            svt_av1_init3smotion_compensation(&pcs->ss_cfg, pcs->ppcs->enhanced_pic->stride_y);
        }
        CdefControls *cdef_ctrls = &pcs->ppcs->cdef_ctrls;
        uint8_t       skip_perc  = pcs->ref_skip_percentage;
        if ((skip_perc > 75 && cdef_ctrls->use_skip_detector) ||
            (scs->vq_ctrls.sharpness_ctrls.cdef && pcs->ppcs->is_noise_level))
            pcs->ppcs->cdef_level = 0;
        else {
            if (cdef_ctrls->use_reference_cdef_fs) {
                if (pcs->slice_type != I_SLICE) {
                    uint8_t lowest_sg  = TOTAL_STRENGTHS - 1;
                    uint8_t highest_sg = 0;
                    // Determine luma pred filter
                    // Add filter from list0
                    EbReferenceObject *ref_obj_l0 =
                        (EbReferenceObject *)pcs->ref_pic_ptr_array[REF_LIST_0][0]->object_ptr;
                    for (uint8_t fs = 0; fs < ref_obj_l0->ref_cdef_strengths_num; fs++) {
                        if (ref_obj_l0->ref_cdef_strengths[0][fs] < lowest_sg)
                            lowest_sg = ref_obj_l0->ref_cdef_strengths[0][fs];
                        if (ref_obj_l0->ref_cdef_strengths[0][fs] > highest_sg)
                            highest_sg = ref_obj_l0->ref_cdef_strengths[0][fs];
                    }
                    if (pcs->slice_type == B_SLICE) {
                        // Add filter from list1
                        EbReferenceObject *ref_obj_l1 =
                            (EbReferenceObject *)pcs->ref_pic_ptr_array[REF_LIST_1][0]->object_ptr;
                        for (uint8_t fs = 0; fs < ref_obj_l1->ref_cdef_strengths_num; fs++) {
                            if (ref_obj_l1->ref_cdef_strengths[0][fs] < lowest_sg)
                                lowest_sg = ref_obj_l1->ref_cdef_strengths[0][fs];
                            if (ref_obj_l1->ref_cdef_strengths[0][fs] > highest_sg)
                                highest_sg = ref_obj_l1->ref_cdef_strengths[0][fs];
                        }
                    }
                    if (pcs->rtc_tune) {
                        int8_t mid_filter     = MIN(63, MAX(0, MAX(lowest_sg, highest_sg)));
                        cdef_ctrls->pred_y_f  = mid_filter;
                        cdef_ctrls->pred_uv_f = 0;
                    } else {
                        int8_t mid_filter     = MIN(63, MAX(0, (lowest_sg + highest_sg) / 2));
                        cdef_ctrls->pred_y_f  = mid_filter;
                        cdef_ctrls->pred_uv_f = 0;
                    }
                    cdef_ctrls->first_pass_fs_num          = 0;
                    cdef_ctrls->default_second_pass_fs_num = 0;
                    // Set cdef to off if pred is.
                    if ((cdef_ctrls->pred_y_f == 0) && (cdef_ctrls->pred_uv_f == 0))
                        pcs->ppcs->cdef_level = 0;
                }
            } else if (cdef_ctrls->search_best_ref_fs) {
                if (pcs->slice_type != I_SLICE) {
                    cdef_ctrls->first_pass_fs_num          = 1;
                    cdef_ctrls->default_second_pass_fs_num = 0;

                    // Add filter from list0, if not the same as the default
                    EbReferenceObject *ref_obj_l0 =
                        (EbReferenceObject *)pcs->ref_pic_ptr_array[REF_LIST_0][0]->object_ptr;
                    if (ref_obj_l0->ref_cdef_strengths[0][0] != cdef_ctrls->default_first_pass_fs[0]) {
                        cdef_ctrls->default_first_pass_fs[1] = ref_obj_l0->ref_cdef_strengths[0][0];
                        (cdef_ctrls->first_pass_fs_num)++;
                    }

                    if (pcs->slice_type == B_SLICE) {
                        EbReferenceObject *ref_obj_l1 =
                            (EbReferenceObject *)pcs->ref_pic_ptr_array[REF_LIST_1][0]->object_ptr;
                        // Add filter from list1, if different from default filter and list0 filter
                        if (ref_obj_l1->ref_cdef_strengths[0][0] != cdef_ctrls->default_first_pass_fs[0] &&
                            ref_obj_l1->ref_cdef_strengths[0][0] !=
                                cdef_ctrls->default_first_pass_fs[cdef_ctrls->first_pass_fs_num - 1]) {
                            cdef_ctrls->default_first_pass_fs[cdef_ctrls->first_pass_fs_num] =
                                ref_obj_l1->ref_cdef_strengths[0][0];
                            (cdef_ctrls->first_pass_fs_num)++;

                            // Chroma
                            if (ref_obj_l0->ref_cdef_strengths[1][0] == cdef_ctrls->default_first_pass_fs_uv[0] &&
                                ref_obj_l1->ref_cdef_strengths[1][0] == cdef_ctrls->default_first_pass_fs_uv[0]) {
                                cdef_ctrls->default_first_pass_fs_uv[0] = -1;
                                cdef_ctrls->default_first_pass_fs_uv[1] = -1;
                            }
                        }
                        // if list0/list1 filters are the same, skip CDEF search, and use the filter selected by the ref frames
                        else if (cdef_ctrls->first_pass_fs_num == 2 &&
                                 ref_obj_l0->ref_cdef_strengths[0][0] == ref_obj_l1->ref_cdef_strengths[0][0]) {
                            cdef_ctrls->use_reference_cdef_fs = 1;

                            cdef_ctrls->pred_y_f  = ref_obj_l0->ref_cdef_strengths[0][0];
                            cdef_ctrls->pred_uv_f = MIN(
                                63,
                                MAX(0,
                                    (ref_obj_l0->ref_cdef_strengths[1][0] + ref_obj_l1->ref_cdef_strengths[1][0]) / 2));
                            cdef_ctrls->first_pass_fs_num          = 0;
                            cdef_ctrls->default_second_pass_fs_num = 0;
                        }
                    }
                    // Chroma
                    else if (ref_obj_l0->ref_cdef_strengths[1][0] == cdef_ctrls->default_first_pass_fs_uv[0]) {
                        cdef_ctrls->default_first_pass_fs_uv[0] = -1;
                        cdef_ctrls->default_first_pass_fs_uv[1] = -1;
                    }

                    // Set cdef to off if pred luma is.
                    if (cdef_ctrls->first_pass_fs_num == 1)
                        pcs->ppcs->cdef_level = 0;
                }
            }
        }

        if (scs->vq_ctrls.sharpness_ctrls.restoration && pcs->ppcs->is_noise_level) {
            pcs->ppcs->enable_restoration = 0;
        }

        // Post the results to the MD processes
        uint16_t tg_count = pcs->ppcs->tile_group_cols * pcs->ppcs->tile_group_rows;
        for (uint16_t tile_group_idx = 0; tile_group_idx < tg_count; tile_group_idx++) {
            svt_get_empty_object(context_ptr->mode_decision_configuration_output_fifo_ptr, &enc_dec_tasks_wrapper);

            EncDecTasks *enc_dec_tasks      = (EncDecTasks *)enc_dec_tasks_wrapper->object_ptr;
            enc_dec_tasks->pcs_wrapper      = rc_results->pcs_wrapper;
            enc_dec_tasks->input_type       = rc_results->superres_recode ? ENCDEC_TASKS_SUPERRES_INPUT
                                                                          : ENCDEC_TASKS_MDC_INPUT;
            enc_dec_tasks->tile_group_index = tile_group_idx;

            // Post the Full Results Object
            svt_post_full_object(enc_dec_tasks_wrapper);

            if (rc_results->superres_recode) {
                // for superres input, only send one task
                break;
            }
        }

        // Release Rate Control Results
        svt_release_object(rc_results_wrapper);
    }

    return NULL;
}
