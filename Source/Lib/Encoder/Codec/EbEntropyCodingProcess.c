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
#include <stdio.h>
#include "EbEncHandle.h"
#include "EbEntropyCodingProcess.h"
#include "EbEncDecResults.h"
#include "EbEntropyCodingResults.h"
#include "EbRateControlTasks.h"
#include "EbCabacContextModel.h"
#include "EbLog.h"
#include "common_dsp_rtcd.h"
#define AV1_MIN_TILE_SIZE_BYTES 1
void svt_av1_reset_loop_restoration(PictureControlSet *piCSetPtr, uint16_t tile_idx);

static void rest_context_dctor(EbPtr p) {
    EbThreadContext      *thread_context_ptr = (EbThreadContext *)p;
    EntropyCodingContext *obj                = (EntropyCodingContext *)thread_context_ptr->priv;
    EB_FREE_ARRAY(obj);
}

/******************************************************
 * Enc Dec Context Constructor
 ******************************************************/
EbErrorType entropy_coding_context_ctor(EbThreadContext   *thread_context_ptr,
                                        const EbEncHandle *enc_handle_ptr, int index,
                                        int rate_control_index) {
    EntropyCodingContext *context_ptr;
    EB_CALLOC_ARRAY(context_ptr, 1);
    thread_context_ptr->priv  = context_ptr;
    thread_context_ptr->dctor = rest_context_dctor;

    context_ptr->is_16bit =
        (Bool)(enc_handle_ptr->scs_instance_array[0]->scs_ptr->static_config.encoder_bit_depth >
               EB_EIGHT_BIT);
    ;

    // Input/Output System Resource Manager FIFOs
    context_ptr->enc_dec_input_fifo_ptr = svt_system_resource_get_consumer_fifo(
        enc_handle_ptr->rest_results_resource_ptr, index);
    context_ptr->entropy_coding_output_fifo_ptr = svt_system_resource_get_producer_fifo(
        enc_handle_ptr->entropy_coding_results_resource_ptr, index);
    context_ptr->rate_control_output_fifo_ptr = svt_system_resource_get_producer_fifo(
        enc_handle_ptr->rate_control_tasks_resource_ptr, rate_control_index);

    return EB_ErrorNone;
}

/***********************************************
 * Entropy Coding Reset Neighbor Arrays
 ***********************************************/
static void entropy_coding_reset_neighbor_arrays(PictureControlSet *pcs_ptr, uint16_t tile_idx) {
    neighbor_array_unit_reset(pcs_ptr->mode_type_neighbor_array[tile_idx]);

    neighbor_array_unit_reset(pcs_ptr->partition_context_neighbor_array[tile_idx]);

    neighbor_array_unit_reset(pcs_ptr->skip_flag_neighbor_array[tile_idx]);

    neighbor_array_unit_reset(pcs_ptr->skip_coeff_neighbor_array[tile_idx]);
    neighbor_array_unit_reset(pcs_ptr->luma_dc_sign_level_coeff_neighbor_array[tile_idx]);
    neighbor_array_unit_reset(pcs_ptr->cb_dc_sign_level_coeff_neighbor_array[tile_idx]);
    neighbor_array_unit_reset(pcs_ptr->cr_dc_sign_level_coeff_neighbor_array[tile_idx]);
    neighbor_array_unit_reset(pcs_ptr->ref_frame_type_neighbor_array[tile_idx]);

    neighbor_array_unit_reset(pcs_ptr->intra_luma_mode_neighbor_array[tile_idx]);
    neighbor_array_unit_reset32(pcs_ptr->interpolation_type_neighbor_array[tile_idx]);
    neighbor_array_unit_reset(pcs_ptr->txfm_context_array[tile_idx]);
    neighbor_array_unit_reset(pcs_ptr->segmentation_id_pred_array[tile_idx]);
    return;
}

void av1_get_syntax_rate_from_cdf(int32_t *costs, const AomCdfProb *cdf, const int32_t *inv_map);

void svt_av1_cost_tokens_from_cdf(int32_t *costs, const AomCdfProb *cdf, const int32_t *inv_map) {
    // int32_t i;
    // AomCdfProb prev_cdf = 0;
    // for (i = 0;; ++i) {
    //     AomCdfProb p15 = AOM_ICDF(cdf[i]) - prev_cdf;
    //     p15 = (p15 < EC_MIN_PROB) ? EC_MIN_PROB : p15;
    //     prev_cdf = AOM_ICDF(cdf[i]);
    //
    //     if (inv_map)
    //         costs[inv_map[i]] = av1_cost_symbol(p15);
    //     else
    //         costs[i] = av1_cost_symbol(p15);
    //
    //     // Stop once we reach the end of the CDF
    //     if (cdf[i] == AOM_ICDF(CDF_PROB_TOP)) break;
    // }

    av1_get_syntax_rate_from_cdf(costs, cdf, inv_map);
}

static void build_nmv_component_cost_table(int32_t *mvcost, const NmvComponent *const mvcomp,
                                           MvSubpelPrecision precision) {
    int32_t i, v;
    int32_t sign_cost[2], class_cost[MV_CLASSES], class0_cost[CLASS0_SIZE];
    int32_t bits_cost[MV_OFFSET_BITS][2];
    int32_t class0_fp_cost[CLASS0_SIZE][MV_FP_SIZE], fp_cost[MV_FP_SIZE];
    int32_t class0_hp_cost[2], hp_cost[2];

    svt_av1_cost_tokens_from_cdf(sign_cost, mvcomp->sign_cdf, NULL);
    svt_av1_cost_tokens_from_cdf(class_cost, mvcomp->classes_cdf, NULL);
    svt_av1_cost_tokens_from_cdf(class0_cost, mvcomp->class0_cdf, NULL);
    for (i = 0; i < MV_OFFSET_BITS; ++i)
        svt_av1_cost_tokens_from_cdf(bits_cost[i], mvcomp->bits_cdf[i], NULL);
    for (i = 0; i < CLASS0_SIZE; ++i)
        svt_av1_cost_tokens_from_cdf(class0_fp_cost[i], mvcomp->class0_fp_cdf[i], NULL);
    svt_av1_cost_tokens_from_cdf(fp_cost, mvcomp->fp_cdf, NULL);

    if (precision > MV_SUBPEL_LOW_PRECISION) {
        svt_av1_cost_tokens_from_cdf(class0_hp_cost, mvcomp->class0_hp_cdf, NULL);
        svt_av1_cost_tokens_from_cdf(hp_cost, mvcomp->hp_cdf, NULL);
    }
    mvcost[0] = 0;
    for (v = 1; v <= MV_MAX; ++v) {
        int32_t z, c, o, d, e, f, cost = 0;
        z = v - 1;
        c = svt_av1_get_mv_class(z, &o);
        cost += class_cost[c];
        d = (o >> 3); /* int32_t mv data */
        f = (o >> 1) & 3; /* fractional pel mv data */
        e = (o & 1); /* high precision mv data */
        if (c == MV_CLASS_0)
            cost += class0_cost[d];
        else {
            const int32_t b = c + CLASS0_BITS - 1; /* number of bits */
            for (i = 0; i < b; ++i) cost += bits_cost[i][((d >> i) & 1)];
        }
        if (precision > MV_SUBPEL_NONE) {
            if (c == MV_CLASS_0)
                cost += class0_fp_cost[d][f];
            else
                cost += fp_cost[f];
            if (precision > MV_SUBPEL_LOW_PRECISION) {
                if (c == MV_CLASS_0)
                    cost += class0_hp_cost[e];
                else
                    cost += hp_cost[e];
            }
        }
        mvcost[v]  = cost + sign_cost[0];
        mvcost[-v] = cost + sign_cost[1];
    }
}
void svt_av1_build_nmv_cost_table(int32_t *mvjoint, int32_t *mvcost[2], const NmvContext *ctx,
                                  MvSubpelPrecision precision) {
    svt_av1_cost_tokens_from_cdf(mvjoint, ctx->joints_cdf, NULL);
    build_nmv_component_cost_table(mvcost[0], &ctx->comps[0], precision);
    build_nmv_component_cost_table(mvcost[1], &ctx->comps[1], precision);
}

/**************************************************
 * Reset Entropy Coding Picture
 **************************************************/
static void reset_entropy_coding_picture(EntropyCodingContext *context_ptr,
                                         PictureControlSet *pcs_ptr, SequenceControlSet *scs_ptr) {
    uint16_t tile_cnt = pcs_ptr->parent_pcs_ptr->av1_cm->tiles_info.tile_rows *
        pcs_ptr->parent_pcs_ptr->av1_cm->tiles_info.tile_cols;
    uint16_t tile_idx = 0;
    uint32_t entropy_coding_qp;

    context_ptr->is_16bit = (Bool)(scs_ptr->static_config.encoder_bit_depth > EB_EIGHT_BIT);
    FrameHeader *frm_hdr  = &pcs_ptr->parent_pcs_ptr->frm_hdr;
    // Asuming cb and cr offset to be the same for chroma QP in both slice and pps for lambda computation
    entropy_coding_qp = pcs_ptr->parent_pcs_ptr->frm_hdr.quantization_params.base_q_idx;

    for (tile_idx = 0; tile_idx < tile_cnt; tile_idx++) {
        pcs_ptr->parent_pcs_ptr->prev_qindex[tile_idx] =
            pcs_ptr->parent_pcs_ptr->frm_hdr.quantization_params.base_q_idx;
    }
    if (pcs_ptr->parent_pcs_ptr->frm_hdr.allow_intrabc)
        assert(pcs_ptr->parent_pcs_ptr->frm_hdr.delta_lf_params.delta_lf_present == 0);
    if (pcs_ptr->parent_pcs_ptr->frm_hdr.delta_lf_params.delta_lf_present) {
        pcs_ptr->parent_pcs_ptr->prev_delta_lf_from_base = 0;
        const int32_t frame_lf_count                     = pcs_ptr->parent_pcs_ptr->monochrome == 0
                                ? FRAME_LF_COUNT
                                : FRAME_LF_COUNT - 2;
        for (int32_t lf_id = 0; lf_id < frame_lf_count; ++lf_id)
            pcs_ptr->parent_pcs_ptr->prev_delta_lf[lf_id] = 0;
    }

    // pass the ent
    for (tile_idx = 0; tile_idx < tile_cnt; tile_idx++) {
        OutputBitstreamUnit *output_bitstream_ptr =
            (OutputBitstreamUnit *)(pcs_ptr->entropy_coding_info[tile_idx]
                                        ->entropy_coder_ptr->ec_output_bitstream_ptr);
        //****************************************************************//
        pcs_ptr->entropy_coding_info[tile_idx]->entropy_coder_ptr->ec_writer.allow_update_cdf =
            !pcs_ptr->parent_pcs_ptr->large_scale_tile;
        pcs_ptr->entropy_coding_info[tile_idx]->entropy_coder_ptr->ec_writer.allow_update_cdf =
            pcs_ptr->entropy_coding_info[tile_idx]->entropy_coder_ptr->ec_writer.allow_update_cdf &&
            !frm_hdr->disable_cdf_update;
        aom_start_encode(&pcs_ptr->entropy_coding_info[tile_idx]->entropy_coder_ptr->ec_writer,
                         output_bitstream_ptr);
        // ADD Reset here
        if (pcs_ptr->parent_pcs_ptr->frm_hdr.primary_ref_frame != PRIMARY_REF_NONE)
            svt_memcpy(
                pcs_ptr->entropy_coding_info[tile_idx]->entropy_coder_ptr->fc,
                &pcs_ptr->ref_frame_context[pcs_ptr->parent_pcs_ptr->frm_hdr.primary_ref_frame],
                sizeof(FRAME_CONTEXT));
        else
            reset_entropy_coder(scs_ptr->encode_context_ptr,
                                pcs_ptr->entropy_coding_info[tile_idx]->entropy_coder_ptr,
                                entropy_coding_qp,
                                pcs_ptr->slice_type);

        entropy_coding_reset_neighbor_arrays(pcs_ptr, tile_idx);
    }
    return;
}

/* Entropy Coding */

/*********************************************************************************
*
* @brief
*  The Entropy Coding process is responsible for producing an AV1 conformant bitstream for each frame.
*
* @par Description:
*  The entropy coder is a frame-based process and is based on multi-symbol arithmetic range coding.
*  It takes as input the coding decisions and information for each block and produces as output the bitstream
*  for each frame.
*
* @param[in] Coding Decisions
*  Coding decisions and information for each block.
*
* @param[out] bitstream
*  Bitstream for each block
*
********************************************************************************/
void *entropy_coding_kernel(void *input_ptr) {
    // Context & SCS & PCS
    EbThreadContext      *thread_context_ptr = (EbThreadContext *)input_ptr;
    EntropyCodingContext *context_ptr        = (EntropyCodingContext *)thread_context_ptr->priv;

    // Input
    EbObjectWrapper *rest_results_wrapper_ptr;

    // Output
    EbObjectWrapper      *entropy_coding_results_wrapper_ptr;
    EntropyCodingResults *entropy_coding_results_ptr;

    for (;;) {
        // Get Mode Decision Results
        EB_GET_FULL_OBJECT(context_ptr->enc_dec_input_fifo_ptr, &rest_results_wrapper_ptr);

        RestResults       *rest_results_ptr = (RestResults *)rest_results_wrapper_ptr->object_ptr;
        PictureControlSet *pcs_ptr          = (PictureControlSet *)
                                         rest_results_ptr->pcs_wrapper_ptr->object_ptr;
        SequenceControlSet *scs_ptr = pcs_ptr->scs_ptr;
        // SB Constants

        uint8_t sb_size = (uint8_t)scs_ptr->sb_size;

        uint8_t  sb_size_log2    = (uint8_t)svt_log2f(sb_size);
        uint32_t pic_width_in_sb = (pcs_ptr->parent_pcs_ptr->aligned_width + sb_size - 1) >>
            sb_size_log2;
        uint16_t         tile_idx        = rest_results_ptr->tile_index;
        Av1Common *const cm              = pcs_ptr->parent_pcs_ptr->av1_cm;
        const uint16_t   tile_cnt        = cm->tiles_info.tile_rows * cm->tiles_info.tile_cols;
        const uint16_t   tile_col        = tile_idx % cm->tiles_info.tile_cols;
        const uint16_t   tile_row        = tile_idx / cm->tiles_info.tile_cols;
        const uint16_t   tile_sb_start_x = cm->tiles_info.tile_col_start_mi[tile_col] >>
            scs_ptr->seq_header.sb_size_log2;
        const uint16_t tile_sb_start_y = cm->tiles_info.tile_row_start_mi[tile_row] >>
            scs_ptr->seq_header.sb_size_log2;

        uint16_t tile_width_in_sb = (cm->tiles_info.tile_col_start_mi[tile_col + 1] -
                                     cm->tiles_info.tile_col_start_mi[tile_col]) >>
            scs_ptr->seq_header.sb_size_log2;
        uint16_t tile_height_in_sb = (cm->tiles_info.tile_row_start_mi[tile_row + 1] -
                                      cm->tiles_info.tile_row_start_mi[tile_row]) >>
            scs_ptr->seq_header.sb_size_log2;

        Bool frame_entropy_done = FALSE;

        svt_block_on_mutex(pcs_ptr->entropy_coding_pic_mutex);
        if (pcs_ptr->entropy_coding_pic_reset_flag) {
            pcs_ptr->entropy_coding_pic_reset_flag = FALSE;

            reset_entropy_coding_picture(context_ptr, pcs_ptr, scs_ptr);
        }
        svt_release_mutex(pcs_ptr->entropy_coding_pic_mutex);

#if TURN_OFF_EC_FIRST_PASS
        if (scs_ptr->static_config.pass != ENC_FIRST_PASS &&
            !is_pic_skipped(pcs_ptr->parent_pcs_ptr)) {
#endif
            for (uint32_t y_sb_index = 0; y_sb_index < tile_height_in_sb; ++y_sb_index) {
                for (uint32_t x_sb_index = 0; x_sb_index < tile_width_in_sb; ++x_sb_index) {
                    uint16_t    sb_index = (uint16_t)((x_sb_index + tile_sb_start_x) +
                                                   (y_sb_index + tile_sb_start_y) *
                                                       pic_width_in_sb);
                    SuperBlock *sb_ptr   = pcs_ptr->sb_ptr_array[sb_index];

                    context_ptr->sb_origin_x = (x_sb_index + tile_sb_start_x) << sb_size_log2;
                    context_ptr->sb_origin_y = (y_sb_index + tile_sb_start_y) << sb_size_log2;
                    if (x_sb_index == 0 && y_sb_index == 0) {
                        svt_av1_reset_loop_restoration(pcs_ptr, tile_idx);
                        context_ptr->tok = pcs_ptr->tile_tok[tile_row][tile_col];
                    }

                    EbPictureBufferDesc *coeff_picture_ptr =
                        pcs_ptr->parent_pcs_ptr->enc_dec_ptr->quantized_coeff[sb_index];
                    write_sb(context_ptr,
                             sb_ptr,
                             pcs_ptr,
                             tile_idx,
                             pcs_ptr->entropy_coding_info[tile_idx]->entropy_coder_ptr,
                             coeff_picture_ptr);
                }
            }
#if TURN_OFF_EC_FIRST_PASS
        }
#endif
        Bool pic_ready = TRUE;

        // Current tile ready
        encode_slice_finish(pcs_ptr->entropy_coding_info[tile_idx]->entropy_coder_ptr);

        svt_block_on_mutex(pcs_ptr->entropy_coding_pic_mutex);
        pcs_ptr->entropy_coding_info[tile_idx]->entropy_coding_tile_done = TRUE;
        for (uint16_t i = 0; i < tile_cnt; i++) {
            if (pcs_ptr->entropy_coding_info[i]->entropy_coding_tile_done == FALSE) {
                pic_ready = FALSE;
                break;
            }
        }
        svt_release_mutex(pcs_ptr->entropy_coding_pic_mutex);
        if (pic_ready) {
            if (pcs_ptr->parent_pcs_ptr->superres_total_recode_loop == 0) {
                // Release the List 0 Reference Pictures
                for (uint32_t ref_idx = 0; ref_idx < pcs_ptr->parent_pcs_ptr->ref_list0_count;
                     ++ref_idx) {
                    if (pcs_ptr->ref_pic_ptr_array[0][ref_idx] != NULL) {
                        svt_release_object(pcs_ptr->ref_pic_ptr_array[0][ref_idx]);
                    }
                }
                // Release the List 1 Reference Pictures
                for (uint32_t ref_idx = 0; ref_idx < pcs_ptr->parent_pcs_ptr->ref_list1_count;
                     ++ref_idx) {
                    if (pcs_ptr->ref_pic_ptr_array[1][ref_idx] != NULL) {
                        svt_release_object(pcs_ptr->ref_pic_ptr_array[1][ref_idx]);
                    }
                }

                //free palette data
                if (pcs_ptr->tile_tok[0][0])
                    EB_FREE_ARRAY(pcs_ptr->tile_tok[0][0]);
            }
            frame_entropy_done = TRUE;
        }

        if (frame_entropy_done) {
            // Get Empty Entropy Coding Results
            svt_get_empty_object(context_ptr->entropy_coding_output_fifo_ptr,
                                 &entropy_coding_results_wrapper_ptr);
            entropy_coding_results_ptr = (EntropyCodingResults *)
                                             entropy_coding_results_wrapper_ptr->object_ptr;
            entropy_coding_results_ptr->pcs_wrapper_ptr = rest_results_ptr->pcs_wrapper_ptr;

            // Post EntropyCoding Results
            svt_post_full_object(entropy_coding_results_wrapper_ptr);
        }

        // Release Mode Decision Results
        svt_release_object(rest_results_wrapper_ptr);
    }

    return NULL;
}
