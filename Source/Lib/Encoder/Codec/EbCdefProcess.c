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
#include "aom_dsp_rtcd.h"
#include "EbDefinitions.h"
#include "EbEncHandle.h"
#include "EbCdefProcess.h"
#include "EbEncDecResults.h"
#include "EbThreads.h"
#include "EbReferenceObject.h"
#include "EbEncCdef.h"
#include "EbEncDecProcess.h"
#include "EbPictureBufferDesc.h"
#include "EbSequenceControlSet.h"
#include "EbUtility.h"
#include "EbPictureControlSet.h"

void copy_sb8_16(uint16_t *dst, int32_t dstride, const uint8_t *src, int32_t src_voffset,
                 int32_t src_hoffset, int32_t sstride, int32_t vsize, int32_t hsize);

void   *svt_aom_memalign(size_t align, size_t size);
void    svt_aom_free(void *memblk);
void   *svt_aom_malloc(size_t size);
int32_t svt_sb_all_skip(PictureControlSet *pcs_ptr, const Av1Common *const cm, int32_t mi_row,
                        int32_t mi_col);
int32_t svt_sb_compute_cdef_list(PictureControlSet *pcs_ptr, const Av1Common *const cm,
                                 int32_t mi_row, int32_t mi_col, CdefList *dlist, BlockSize bs);
void    finish_cdef_search(PictureControlSet *pcs_ptr);
void    av1_cdef_frame16bit(EncDecContext *context_ptr, SequenceControlSet *scs_ptr,
                            PictureControlSet *pCs);
void    svt_av1_cdef_frame(EncDecContext *context_ptr, SequenceControlSet *scs_ptr,
                           PictureControlSet *pCs);
void    svt_av1_loop_restoration_save_boundary_lines(const Yv12BufferConfig *frame, Av1Common *cm,
                                                     int32_t after_cdef);
void    svt_av1_superres_upscale_frame(struct Av1Common *cm, PictureControlSet *pcs_ptr,
                                       SequenceControlSet *scs_ptr);
void    set_unscaled_input_16bit(PictureControlSet *pcs_ptr);

void get_recon_pic(PictureControlSet *pcs_ptr, EbPictureBufferDesc **recon_ptr, EbBool is_highbd);

/**************************************
 * Cdef Context
 **************************************/
typedef struct CdefContext {
    EbFifo *cdef_input_fifo_ptr;
    EbFifo *cdef_output_fifo_ptr;
} CdefContext;

static void cdef_context_dctor(EbPtr p) {
    EbThreadContext *thread_context_ptr = (EbThreadContext *)p;
    CdefContext     *obj                = (CdefContext *)thread_context_ptr->priv;
    EB_FREE_ARRAY(obj);
}

/******************************************************
 * Cdef Context Constructor
 ******************************************************/
EbErrorType cdef_context_ctor(EbThreadContext   *thread_context_ptr,
                              const EbEncHandle *enc_handle_ptr, int index) {
    CdefContext *context_ptr;
    EB_CALLOC_ARRAY(context_ptr, 1);
    thread_context_ptr->priv  = context_ptr;
    thread_context_ptr->dctor = cdef_context_dctor;

    // Input/Output System Resource Manager FIFOs
    context_ptr->cdef_input_fifo_ptr = svt_system_resource_get_consumer_fifo(
        enc_handle_ptr->dlf_results_resource_ptr, index);
    context_ptr->cdef_output_fifo_ptr = svt_system_resource_get_producer_fifo(
        enc_handle_ptr->cdef_results_resource_ptr, index);

    return EB_ErrorNone;
}

#define default_mse_uv 1040400
/* Search for the best filter strength pair for each 64x64 filter block.
 *
 * For each 64x64 filter block and each plane, search the allowable filter strength pairs.
 * Call cdef_filter_fb() to perform filtering, then compute the MSE for each pair.
*/
void cdef_seg_search(PictureControlSet *pcs_ptr, SequenceControlSet *scs_ptr,
                     uint32_t segment_index) {
    struct PictureParentControlSet *ppcs    = pcs_ptr->parent_pcs_ptr;
    FrameHeader                    *frm_hdr = &ppcs->frm_hdr;
    Av1Common                      *cm      = pcs_ptr->parent_pcs_ptr->av1_cm;
    uint32_t                        x_seg_idx;
    uint32_t                        y_seg_idx;
    uint32_t picture_width_in_b64  = (pcs_ptr->parent_pcs_ptr->aligned_width + 64 - 1) / 64;
    uint32_t picture_height_in_b64 = (pcs_ptr->parent_pcs_ptr->aligned_height + 64 - 1) / 64;
    SEGMENT_CONVERT_IDX_TO_XY(
        segment_index, x_seg_idx, y_seg_idx, pcs_ptr->cdef_segments_column_count);
    uint32_t x_b64_start_idx = SEGMENT_START_IDX(
        x_seg_idx, picture_width_in_b64, pcs_ptr->cdef_segments_column_count);
    uint32_t x_b64_end_idx = SEGMENT_END_IDX(
        x_seg_idx, picture_width_in_b64, pcs_ptr->cdef_segments_column_count);
    uint32_t y_b64_start_idx = SEGMENT_START_IDX(
        y_seg_idx, picture_height_in_b64, pcs_ptr->cdef_segments_row_count);
    uint32_t y_b64_end_idx = SEGMENT_END_IDX(
        y_seg_idx, picture_height_in_b64, pcs_ptr->cdef_segments_row_count);

    const int32_t mi_rows = cm->mi_rows;
    const int32_t mi_cols = cm->mi_cols;
    uint32_t      fbr, fbc;
    int32_t       gi;
    int32_t       nvb, nhb;
    CdefControls *cdef_ctrls                 = &pcs_ptr->parent_pcs_ptr->cdef_ctrls;
    const int     first_pass_fs_num          = cdef_ctrls->first_pass_fs_num;
    const int     default_second_pass_fs_num = cdef_ctrls->default_second_pass_fs_num;
    int32_t       pri_strength;
    uint64_t      curr_mse;
    int32_t       sec_strength;
    uint8_t      *src[3];
    uint8_t      *ref_coeff[3];
    CdefList      dlist[MI_SIZE_128X128 * MI_SIZE_128X128];
    int32_t       stride_src[3];
    int32_t       stride_ref[3];
    int32_t       bsize[3];
    int32_t       mi_wide_l2[3];
    int32_t       mi_high_l2[3];
    int32_t       xdec[3];
    int32_t       ydec[3];
    int32_t       pli;
    int32_t       cdef_count;
    int32_t       coeff_shift = AOMMAX(scs_ptr->static_config.encoder_bit_depth - 8, 0);
    int32_t       nvfb        = (mi_rows + MI_SIZE_64X64 - 1) / MI_SIZE_64X64;
    int32_t       nhfb        = (mi_cols + MI_SIZE_64X64 - 1) / MI_SIZE_64X64;
    int32_t       pri_damping = 3 + (frm_hdr->quantization_params.base_q_idx >> 6);
    int32_t       sec_damping = pri_damping;
    const int32_t num_planes  = 3;

    DECLARE_ALIGNED(32, uint16_t, inbuf[CDEF_INBUF_SIZE]);
    uint16_t *in;
    DECLARE_ALIGNED(32, uint8_t, tmp_dst[1 << (MAX_SB_SIZE_LOG2 * 2)]);
    EbPictureBufferDesc *input_picture_ptr = (EbPictureBufferDesc *)
                                                 pcs_ptr->parent_pcs_ptr->enhanced_picture_ptr;
    EbPictureBufferDesc *recon_picture_ptr;
    if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
        recon_picture_ptr = ((EbReferenceObject *)
                                 pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
                                ->reference_picture;
    else
        recon_picture_ptr = pcs_ptr->parent_pcs_ptr->enc_dec_ptr->recon_picture_ptr;

    for (pli = 0; pli < num_planes; pli++) {
        int32_t subsampling_x = (pli == 0) ? 0 : 1;
        int32_t subsampling_y = (pli == 0) ? 0 : 1;
        xdec[pli]             = subsampling_x;
        ydec[pli]             = subsampling_y;
        bsize[pli]            = ydec[pli] ? (xdec[pli] ? BLOCK_4X4 : BLOCK_8X4)
                                          : (xdec[pli] ? BLOCK_4X8 : BLOCK_8X8);
        mi_wide_l2[pli]       = MI_SIZE_LOG2 - subsampling_x;
        mi_high_l2[pli]       = MI_SIZE_LOG2 - subsampling_y;

        src[pli]        = (uint8_t *)pcs_ptr->src[pli];
        ref_coeff[pli]  = (uint8_t *)pcs_ptr->ref_coeff[pli];
        stride_src[pli] = pli == 0
            ? recon_picture_ptr->stride_y
            : (pli == 1 ? recon_picture_ptr->stride_cb : recon_picture_ptr->stride_cr);
        stride_ref[pli] = pli == 0
            ? input_picture_ptr->stride_y
            : (pli == 1 ? input_picture_ptr->stride_cb : input_picture_ptr->stride_cr);
    }

    in = inbuf + CDEF_VBORDER * CDEF_BSTRIDE + CDEF_HBORDER;
    // This is the SB loop
    for (fbr = y_b64_start_idx; fbr < y_b64_end_idx; ++fbr) {
        for (fbc = x_b64_start_idx; fbc < x_b64_end_idx; ++fbc) {
            int32_t        dirinit    = 0;
            const uint32_t lc         = MI_SIZE_64X64 * fbc;
            const uint32_t lr         = MI_SIZE_64X64 * fbr;
            nhb                       = AOMMIN(MI_SIZE_64X64, mi_cols - lc);
            nvb                       = AOMMIN(MI_SIZE_64X64, mi_rows - lr);
            int32_t           hb_step = 1; //these should be all time with 64x64 SBs
            int32_t           vb_step = 1;
            BlockSize         bs      = BLOCK_64X64;
            ModeInfo        **mi      = pcs_ptr->mi_grid_base + lr * cm->mi_stride + lc;
            const MbModeInfo *mbmi    = &mi[0]->mbmi;
            const BlockSize   sb_type = mbmi->block_mi.sb_type;
            if (((fbc & 1) && (sb_type == BLOCK_128X128 || sb_type == BLOCK_128X64)) ||
                ((fbr & 1) && (sb_type == BLOCK_128X128 || sb_type == BLOCK_64X128)))
                continue;
            if (sb_type == BLOCK_128X128 || sb_type == BLOCK_128X64 || sb_type == BLOCK_64X128)
                bs = sb_type;

            if (bs == BLOCK_128X128 || bs == BLOCK_128X64) {
                nhb     = AOMMIN(MI_SIZE_128X128, cm->mi_cols - lc);
                hb_step = 2;
            }
            if (bs == BLOCK_128X128 || bs == BLOCK_64X128) {
                nvb     = AOMMIN(MI_SIZE_128X128, cm->mi_rows - lr);
                vb_step = 2;
            }
            const uint32_t fb_idx = fbr * nhfb + fbc;
            // No filtering if the entire filter block is skipped
            cdef_count = svt_sb_compute_cdef_list(pcs_ptr, cm, lr, lc, dlist, bs);
            if (cdef_count == 0) {
                pcs_ptr->skip_cdef_seg[fb_idx] = 1;
                continue;
            } else {
                pcs_ptr->skip_cdef_seg[fb_idx] = 0;
            }
            uint8_t(*dir)[CDEF_NBLOCKS][CDEF_NBLOCKS] = &pcs_ptr->cdef_dir_data[fb_idx].dir;
            int32_t(*var)[CDEF_NBLOCKS][CDEF_NBLOCKS] = &pcs_ptr->cdef_dir_data[fb_idx].var;
            for (pli = 0; pli < num_planes; pli++) {
                /* We avoid filtering the pixels for which some of the pixels to
                   average are outside the frame. We could change the filter instead,
                   but it would add special cases for any future vectorization.
                   No need to set pli == 2 because the copy size will be the same as for pli == 1. */
                if (pli < 2)
                    memset(inbuf, (uint8_t)CDEF_VERY_LARGE, sizeof(inbuf[0]) * CDEF_INBUF_SIZE);
                int32_t yoff  = CDEF_VBORDER * (fbr != 0);
                int32_t xoff  = CDEF_HBORDER * (fbc != 0);
                int32_t ysize = (nvb << mi_high_l2[pli]) +
                    CDEF_VBORDER * ((int32_t)fbr + vb_step < nvfb) + yoff;
                int32_t xsize = (nhb << mi_wide_l2[pli]) +
                    CDEF_HBORDER * ((int32_t)fbc + hb_step < nhfb) + xoff;

                copy_sb8_16(&in[(-yoff * CDEF_BSTRIDE - xoff)],
                            CDEF_BSTRIDE,
                            src[pli],
                            (lr << mi_high_l2[pli]) - yoff,
                            (lc << mi_wide_l2[pli]) - xoff,
                            stride_src[pli],
                            ysize,
                            xsize);
                uint8_t subsampling_factor = cdef_ctrls->subsampling_factor;
                /*
                Cap the subsampling for certain block sizes.

                The intrinsics process several lines simultaneously, so blocks can only be subsampled
                a finite amount before there is no more speed gain.  If the space between processed lines
                is too large, the intrinsics will begin accessing memory outside the block.
                */
                if (bsize[pli] == BLOCK_8X8)
                    subsampling_factor = MIN(subsampling_factor, 4);
                else if (bsize[pli] == BLOCK_8X4 || bsize[pli] == BLOCK_4X8)
                    subsampling_factor = MIN(subsampling_factor, 2);
                else if (bsize[pli] == BLOCK_4X4)
                    subsampling_factor = MIN(subsampling_factor, 1);
                /* first cdef stage
                 * Perform the pri_filter strength search for the current sub_block
                 */
                for (gi = 0; gi < first_pass_fs_num; gi++) {
                    if (!pli || (cdef_ctrls->default_first_pass_fs_uv[gi] != -1)) {
                        pri_strength = cdef_ctrls->default_first_pass_fs[gi] / CDEF_SEC_STRENGTHS;
                        sec_strength = cdef_ctrls->default_first_pass_fs[gi] % CDEF_SEC_STRENGTHS;

                        svt_cdef_filter_fb(tmp_dst,
                                           NULL,
                                           0,
                                           in,
                                           xdec[pli],
                                           ydec[pli],
                                           *dir,
                                           &dirinit,
                                           *var,
                                           pli,
                                           dlist,
                                           cdef_count,
                                           pri_strength,
                                           sec_strength + (sec_strength == 3),
                                           pri_damping,
                                           sec_damping,
                                           coeff_shift,
                                           subsampling_factor);
                        curr_mse = svt_compute_cdef_dist_8bit(
                            ref_coeff[pli] + (lr << mi_high_l2[pli]) * stride_ref[pli] +
                                (lc << mi_wide_l2[pli]),
                            stride_ref[pli],
                            tmp_dst,
                            dlist,
                            cdef_count,
                            (BlockSize)bsize[pli],
                            coeff_shift,
                            pli,
                            subsampling_factor);

                        if (pli < 2)
                            pcs_ptr->mse_seg[pli][fb_idx][gi] = curr_mse *
                                subsampling_factor; // Phoenix - *2 b/c using every other line
                        else
                            pcs_ptr->mse_seg[1][fb_idx][gi] += (curr_mse * subsampling_factor);
                    } else
                        pcs_ptr->mse_seg[1][fb_idx][gi] = default_mse_uv * 64;
                }

                /* second cdef stage
                 * Perform the sec_filter strength search for the current sub_block
                 */
                for (gi = first_pass_fs_num; gi < first_pass_fs_num + default_second_pass_fs_num;
                     gi++) {
                    if (!pli ||
                        (cdef_ctrls->default_second_pass_fs_uv[gi - first_pass_fs_num] != -1)) {
                        pri_strength = cdef_ctrls->default_second_pass_fs[gi - first_pass_fs_num] /
                            CDEF_SEC_STRENGTHS;
                        sec_strength = cdef_ctrls->default_second_pass_fs[gi - first_pass_fs_num] %
                            CDEF_SEC_STRENGTHS;

                        svt_cdef_filter_fb(tmp_dst,
                                           NULL,
                                           0,
                                           in,
                                           xdec[pli],
                                           ydec[pli],
                                           *dir,
                                           &dirinit,
                                           *var,
                                           pli,
                                           dlist,
                                           cdef_count,
                                           pri_strength,
                                           sec_strength + (sec_strength == 3),
                                           pri_damping,
                                           sec_damping,
                                           coeff_shift,
                                           subsampling_factor);
                        curr_mse = svt_compute_cdef_dist_8bit(
                            ref_coeff[pli] + (lr << mi_high_l2[pli]) * stride_ref[pli] +
                                (lc << mi_wide_l2[pli]),
                            stride_ref[pli],
                            tmp_dst,
                            dlist,
                            cdef_count,
                            (BlockSize)bsize[pli],
                            coeff_shift,
                            pli,
                            subsampling_factor);

                        if (pli < 2)
                            pcs_ptr->mse_seg[pli][fb_idx][gi] = curr_mse *
                                subsampling_factor; // Phoenix - *2 b/c using every other line
                        else
                            pcs_ptr->mse_seg[1][fb_idx][gi] += (curr_mse * subsampling_factor);
                    } else
                        pcs_ptr->mse_seg[1][fb_idx][gi] = default_mse_uv * 64;
                }
            }
        }
    }
}
/* Search for the best filter strength pair for each 64x64 filter block.
 *
 * For each 64x64 filter block and each plane, search the allowable filter strength pairs.
 * Call cdef_filter_fb() to perform filtering, then compute the MSE for each pair.
*/
void cdef_seg_search16bit(PictureControlSet *pcs_ptr, SequenceControlSet *scs_ptr,
                          uint32_t segment_index) {
    EbPictureBufferDesc *input_pic_ptr = pcs_ptr->input_frame16bit;

    EbPictureBufferDesc *recon_pic_ptr;
    get_recon_pic(pcs_ptr, &recon_pic_ptr, 1);
    struct PictureParentControlSet *ppcs    = pcs_ptr->parent_pcs_ptr;
    FrameHeader                    *frm_hdr = &ppcs->frm_hdr;
    Av1Common                      *cm      = pcs_ptr->parent_pcs_ptr->av1_cm;
    uint32_t                        x_seg_idx;
    uint32_t                        y_seg_idx;
    uint32_t picture_width_in_b64  = (pcs_ptr->parent_pcs_ptr->aligned_width + 64 - 1) / 64;
    uint32_t picture_height_in_b64 = (pcs_ptr->parent_pcs_ptr->aligned_height + 64 - 1) / 64;
    SEGMENT_CONVERT_IDX_TO_XY(
        segment_index, x_seg_idx, y_seg_idx, pcs_ptr->cdef_segments_column_count);
    uint32_t x_b64_start_idx = SEGMENT_START_IDX(
        x_seg_idx, picture_width_in_b64, pcs_ptr->cdef_segments_column_count);
    uint32_t x_b64_end_idx = SEGMENT_END_IDX(
        x_seg_idx, picture_width_in_b64, pcs_ptr->cdef_segments_column_count);
    uint32_t y_b64_start_idx = SEGMENT_START_IDX(
        y_seg_idx, picture_height_in_b64, pcs_ptr->cdef_segments_row_count);
    uint32_t y_b64_end_idx = SEGMENT_END_IDX(
        y_seg_idx, picture_height_in_b64, pcs_ptr->cdef_segments_row_count);

    int32_t       mi_rows = cm->mi_rows;
    int32_t       mi_cols = cm->mi_cols;
    uint32_t      fbr, fbc;
    int32_t       gi;
    int32_t       pri_strength;
    uint64_t      curr_mse;
    int32_t       sec_strength;
    CdefControls *cdef_ctrls                 = &pcs_ptr->parent_pcs_ptr->cdef_ctrls;
    const int     first_pass_fs_num          = cdef_ctrls->first_pass_fs_num;
    const int     default_second_pass_fs_num = cdef_ctrls->default_second_pass_fs_num;
    uint16_t     *src[3];
    uint16_t     *ref_coeff[3];
    CdefList      dlist[MI_SIZE_128X128 * MI_SIZE_128X128];
    int32_t       stride_src[3];
    int32_t       stride_ref[3];
    int32_t       bsize[3];
    int32_t       mi_wide_l2[3];
    int32_t       mi_high_l2[3];
    int32_t       xdec[3];
    int32_t       ydec[3];
    int32_t       pli;
    int32_t       cdef_count;
    int32_t       coeff_shift = AOMMAX(scs_ptr->static_config.encoder_bit_depth - 8, 0);
    int32_t       nvfb        = (mi_rows + MI_SIZE_64X64 - 1) / MI_SIZE_64X64;
    int32_t       nhfb        = (mi_cols + MI_SIZE_64X64 - 1) / MI_SIZE_64X64;
    int32_t       pri_damping = 3 + (frm_hdr->quantization_params.base_q_idx >> 6);
    int32_t       sec_damping = pri_damping;
    const int32_t num_planes  = 3;
    DECLARE_ALIGNED(32, uint16_t, inbuf[CDEF_INBUF_SIZE]);
    uint16_t *in;
    DECLARE_ALIGNED(32, uint16_t, tmp_dst[1 << (MAX_SB_SIZE_LOG2 * 2)]);

    for (pli = 0; pli < num_planes; pli++) {
        int32_t subsampling_x = (pli == 0) ? 0 : 1;
        int32_t subsampling_y = (pli == 0) ? 0 : 1;
        xdec[pli]             = subsampling_x;
        ydec[pli]             = subsampling_y;
        bsize[pli]            = ydec[pli] ? (xdec[pli] ? BLOCK_4X4 : BLOCK_8X4)
                                          : (xdec[pli] ? BLOCK_4X8 : BLOCK_8X8);

        mi_wide_l2[pli] = MI_SIZE_LOG2 - subsampling_x;
        mi_high_l2[pli] = MI_SIZE_LOG2 - subsampling_y;

        src[pli]        = pcs_ptr->src[pli];
        ref_coeff[pli]  = pcs_ptr->ref_coeff[pli];
        stride_src[pli] = pli == 0
            ? recon_pic_ptr->stride_y
            : (pli == 1 ? recon_pic_ptr->stride_cb : recon_pic_ptr->stride_cr);
        stride_ref[pli] = pli == 0
            ? input_pic_ptr->stride_y
            : (pli == 1 ? input_pic_ptr->stride_cb : input_pic_ptr->stride_cr);
    }

    in = inbuf + CDEF_VBORDER * CDEF_BSTRIDE + CDEF_HBORDER;

    for (fbr = y_b64_start_idx; fbr < y_b64_end_idx; ++fbr) {
        for (fbc = x_b64_start_idx; fbc < x_b64_end_idx; ++fbc) {
            const uint32_t lc = MI_SIZE_64X64 * fbc;
            const uint32_t lr = MI_SIZE_64X64 * fbr;
            int32_t        nvb, nhb;
            int32_t        dirinit    = 0;
            nhb                       = AOMMIN(MI_SIZE_64X64, mi_cols - lc);
            nvb                       = AOMMIN(MI_SIZE_64X64, mi_rows - lr);
            int32_t           hb_step = 1; //these should be all time with 64x64 SBs
            int32_t           vb_step = 1;
            BlockSize         bs      = BLOCK_64X64;
            ModeInfo        **mi      = pcs_ptr->mi_grid_base + lr * cm->mi_stride + lc;
            const MbModeInfo *mbmi    = &mi[0]->mbmi;
            const BlockSize   sb_type = mbmi->block_mi.sb_type;
            if (((fbc & 1) && (sb_type == BLOCK_128X128 || sb_type == BLOCK_128X64)) ||
                ((fbr & 1) && (sb_type == BLOCK_128X128 || sb_type == BLOCK_64X128)))
                continue;
            if (sb_type == BLOCK_128X128 || sb_type == BLOCK_128X64 || sb_type == BLOCK_64X128)
                bs = sb_type;
            if (bs == BLOCK_128X128 || bs == BLOCK_128X64) {
                nhb     = AOMMIN(MI_SIZE_128X128, mi_cols - lc);
                hb_step = 2;
            }
            if (bs == BLOCK_128X128 || bs == BLOCK_64X128) {
                nvb     = AOMMIN(MI_SIZE_128X128, mi_rows - lr);
                vb_step = 2;
            }
            const uint32_t fb_idx = fbr * nhfb + fbc;
            cdef_count            = svt_sb_compute_cdef_list(
                pcs_ptr, cm, fbr * MI_SIZE_64X64, fbc * MI_SIZE_64X64, dlist, bs);

            if (cdef_count == 0) {
                pcs_ptr->skip_cdef_seg[fb_idx] = 1;
                continue;
            } else {
                pcs_ptr->skip_cdef_seg[fb_idx] = 0;
            }
            uint8_t(*dir)[CDEF_NBLOCKS][CDEF_NBLOCKS] = &pcs_ptr->cdef_dir_data[fb_idx].dir;
            int32_t(*var)[CDEF_NBLOCKS][CDEF_NBLOCKS] = &pcs_ptr->cdef_dir_data[fb_idx].var;
            for (pli = 0; pli < num_planes; pli++) {
                /* We avoid filtering the pixels for which some of the pixels to
                   average are outside the frame. We could change the filter instead,
                   but it would add special cases for any future vectorization.
                   No need to set pli == 2 because the copy size will be the same as for pli == 1.*/
                if (pli < 2)
                    memset(inbuf, (uint8_t)CDEF_VERY_LARGE, sizeof(inbuf[0]) * CDEF_INBUF_SIZE);
                int32_t yoff  = CDEF_VBORDER * (fbr != 0);
                int32_t xoff  = CDEF_HBORDER * (fbc != 0);
                int32_t ysize = (nvb << mi_high_l2[pli]) +
                    CDEF_VBORDER * ((int32_t)fbr + vb_step < nvfb) + yoff;
                int32_t xsize = (nhb << mi_wide_l2[pli]) +
                    CDEF_HBORDER * ((int32_t)fbc + hb_step < nhfb) + xoff;

                copy_sb16_16(&in[(-yoff * CDEF_BSTRIDE - xoff)],
                             CDEF_BSTRIDE,
                             src[pli],
                             (lr << mi_high_l2[pli]) - yoff,
                             (lc << mi_wide_l2[pli]) - xoff,
                             stride_src[pli],
                             ysize,
                             xsize);

                uint8_t subsampling_factor = cdef_ctrls->subsampling_factor;
                /*
                Cap the subsampling for certain block sizes.

                The intrinsics process several lines simultaneously, so blocks can only be subsampled
                a finite amount before there is no more speed gain.  If the space between processed lines
                is too large, the intrinsics will begin accessing memory outside the block.
                */
                if (bsize[pli] == BLOCK_8X8)
#ifndef NON_AVX512_SUPPORT
                    subsampling_factor = MIN(subsampling_factor, 2);
#else
                    subsampling_factor = MIN(subsampling_factor, 4);
#endif
                else if (bsize[pli] == BLOCK_8X4 || bsize[pli] == BLOCK_4X8)
                    subsampling_factor = MIN(subsampling_factor, 2);
                else if (bsize[pli] == BLOCK_4X4)
                    subsampling_factor = MIN(subsampling_factor, 1);
                /* first cdef stage
                 * Perform the pri_filter strength search for the current sub_block
                 */
                for (gi = 0; gi < first_pass_fs_num; gi++) {
                    if (!pli || (cdef_ctrls->default_first_pass_fs_uv[gi] != -1)) {
                        pri_strength = cdef_ctrls->default_first_pass_fs[gi] / CDEF_SEC_STRENGTHS;
                        sec_strength = cdef_ctrls->default_first_pass_fs[gi] % CDEF_SEC_STRENGTHS;

                        svt_cdef_filter_fb(NULL,
                                           tmp_dst,
                                           0,
                                           in,
                                           xdec[pli],
                                           ydec[pli],
                                           *dir,
                                           &dirinit,
                                           *var,
                                           pli,
                                           dlist,
                                           cdef_count,
                                           pri_strength,
                                           sec_strength + (sec_strength == 3),
                                           pri_damping,
                                           sec_damping,
                                           coeff_shift,
                                           subsampling_factor);
                        curr_mse = svt_compute_cdef_dist_16bit(
                            ref_coeff[pli] + (lr << mi_high_l2[pli]) * stride_ref[pli] +
                                (lc << mi_wide_l2[pli]),
                            stride_ref[pli],
                            tmp_dst,
                            dlist,
                            cdef_count,
                            (BlockSize)bsize[pli],
                            coeff_shift,
                            pli,
                            subsampling_factor);

                        if (pli < 2)
                            pcs_ptr->mse_seg[pli][fb_idx][gi] = curr_mse * subsampling_factor;
                        else
                            pcs_ptr->mse_seg[1][fb_idx][gi] += (curr_mse * subsampling_factor);
                    } else
                        pcs_ptr->mse_seg[1][fb_idx][gi] = default_mse_uv * 64;
                }

                /* second cdef stage
                 * Perform the sec_filter strength search for the current sub_block
                 */

                for (gi = first_pass_fs_num; gi < first_pass_fs_num + default_second_pass_fs_num;
                     gi++) {
                    if (!pli ||
                        (cdef_ctrls->default_second_pass_fs_uv[gi - first_pass_fs_num] != 1)) {
                        pri_strength = cdef_ctrls->default_second_pass_fs[gi - first_pass_fs_num] /
                            CDEF_SEC_STRENGTHS;
                        sec_strength = cdef_ctrls->default_second_pass_fs[gi - first_pass_fs_num] %
                            CDEF_SEC_STRENGTHS;

                        svt_cdef_filter_fb(NULL,
                                           tmp_dst,
                                           0,
                                           in,
                                           xdec[pli],
                                           ydec[pli],
                                           *dir,
                                           &dirinit,
                                           *var,
                                           pli,
                                           dlist,
                                           cdef_count,
                                           pri_strength,
                                           sec_strength + (sec_strength == 3),
                                           pri_damping,
                                           sec_damping,
                                           coeff_shift,
                                           subsampling_factor);
                        curr_mse = svt_compute_cdef_dist_16bit(
                            ref_coeff[pli] + (lr << mi_high_l2[pli]) * stride_ref[pli] +
                                (lc << mi_wide_l2[pli]),
                            stride_ref[pli],
                            tmp_dst,
                            dlist,
                            cdef_count,
                            (BlockSize)bsize[pli],
                            coeff_shift,
                            pli,
                            subsampling_factor);

                        if (pli < 2)
                            pcs_ptr->mse_seg[pli][fb_idx][gi] = curr_mse * subsampling_factor;
                        else
                            pcs_ptr->mse_seg[1][fb_idx][gi] += (curr_mse * subsampling_factor);
                    } else
                        pcs_ptr->mse_seg[1][fb_idx][gi] = default_mse_uv * 64;
                }
            }
        }
    }
}

/******************************************************
 * CDEF Kernel
 ******************************************************/
void *cdef_kernel(void *input_ptr) {
    // Context & SCS & PCS
    EbThreadContext    *thread_context_ptr = (EbThreadContext *)input_ptr;
    CdefContext        *context_ptr        = (CdefContext *)thread_context_ptr->priv;
    PictureControlSet  *pcs_ptr;
    SequenceControlSet *scs_ptr;

    //// Input
    EbObjectWrapper *dlf_results_wrapper_ptr;
    DlfResults      *dlf_results_ptr;

    //// Output
    EbObjectWrapper *cdef_results_wrapper_ptr;
    CdefResults     *cdef_results_ptr;

    // SB Loop variables

    for (;;) {
        FrameHeader *frm_hdr;

        // Get DLF Results
        EB_GET_FULL_OBJECT(context_ptr->cdef_input_fifo_ptr, &dlf_results_wrapper_ptr);

        dlf_results_ptr = (DlfResults *)dlf_results_wrapper_ptr->object_ptr;
        pcs_ptr         = (PictureControlSet *)dlf_results_ptr->pcs_wrapper_ptr->object_ptr;
        scs_ptr         = (SequenceControlSet *)pcs_ptr->scs_wrapper_ptr->object_ptr;

        EbBool     is_16bit      = scs_ptr->is_16bit_pipeline;
        Av1Common *cm            = pcs_ptr->parent_pcs_ptr->av1_cm;
        frm_hdr                  = &pcs_ptr->parent_pcs_ptr->frm_hdr;
        CdefControls *cdef_ctrls = &pcs_ptr->parent_pcs_ptr->cdef_ctrls;
        if (!cdef_ctrls->use_reference_cdef_fs) {
            if (scs_ptr->seq_header.cdef_level && pcs_ptr->parent_pcs_ptr->cdef_level) {
                if (is_16bit)
                    cdef_seg_search16bit(pcs_ptr, scs_ptr, dlf_results_ptr->segment_index);
                else
                    cdef_seg_search(pcs_ptr, scs_ptr, dlf_results_ptr->segment_index);
            }
        }
        //all seg based search is done. update total processed segments. if all done, finish the search and perfrom application.
        svt_block_on_mutex(pcs_ptr->cdef_search_mutex);

        pcs_ptr->tot_seg_searched_cdef++;
        if (pcs_ptr->tot_seg_searched_cdef == pcs_ptr->cdef_segments_total_count) {
            // SVT_LOG("    CDEF all seg here  %i\n", pcs_ptr->picture_number);
            if (scs_ptr->seq_header.cdef_level && pcs_ptr->parent_pcs_ptr->cdef_level) {
                finish_cdef_search(pcs_ptr);

                if (scs_ptr->seq_header.enable_restoration != 0 ||
                    pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ||
                    scs_ptr->static_config.recon_enabled) {
                    // Do application iff there are non-zero filters
                    if (frm_hdr->cdef_params.cdef_y_strength[0] != 0 ||
                        frm_hdr->cdef_params.cdef_uv_strength[0] != 0 ||
                        pcs_ptr->parent_pcs_ptr->nb_cdef_strengths != 1) {
                        if (is_16bit)
                            av1_cdef_frame16bit(0, scs_ptr, pcs_ptr);
                        else
                            svt_av1_cdef_frame(0, scs_ptr, pcs_ptr);
                    }
                }
            } else {
                frm_hdr->cdef_params.cdef_bits             = 0;
                frm_hdr->cdef_params.cdef_y_strength[0]    = 0;
                pcs_ptr->parent_pcs_ptr->nb_cdef_strengths = 1;
                frm_hdr->cdef_params.cdef_uv_strength[0]   = 0;
            }

            //restoration prep
            if (scs_ptr->seq_header.enable_restoration) {
                svt_av1_loop_restoration_save_boundary_lines(cm->frame_to_show, cm, 1);
            }

            // ------- start: Normative upscaling - super-resolution tool
            if (frm_hdr->allow_intrabc == 0 && !av1_superres_unscaled(&cm->frm_size)) {
                svt_av1_superres_upscale_frame(cm, pcs_ptr, scs_ptr);

                if (is_16bit) {
                    set_unscaled_input_16bit(pcs_ptr);
                }
            }
            // ------- end: Normative upscaling - super-resolution tool

            pcs_ptr->rest_segments_column_count = scs_ptr->rest_segment_column_count;
            pcs_ptr->rest_segments_row_count    = scs_ptr->rest_segment_row_count;
            pcs_ptr->rest_segments_total_count  = (uint16_t)(pcs_ptr->rest_segments_column_count *
                                                            pcs_ptr->rest_segments_row_count);
            pcs_ptr->tot_seg_searched_rest      = 0;
            pcs_ptr->parent_pcs_ptr->av1_cm->use_boundaries_in_rest_search =
                scs_ptr->use_boundaries_in_rest_search;
            pcs_ptr->rest_extend_flag[0] = EB_FALSE;
            pcs_ptr->rest_extend_flag[1] = EB_FALSE;
            pcs_ptr->rest_extend_flag[2] = EB_FALSE;

            uint32_t segment_index;
            for (segment_index = 0; segment_index < pcs_ptr->rest_segments_total_count;
                 ++segment_index) {
                // Get Empty Cdef Results to Rest
                svt_get_empty_object(context_ptr->cdef_output_fifo_ptr, &cdef_results_wrapper_ptr);
                cdef_results_ptr = (struct CdefResults *)cdef_results_wrapper_ptr->object_ptr;
                cdef_results_ptr->pcs_wrapper_ptr = dlf_results_ptr->pcs_wrapper_ptr;
                cdef_results_ptr->segment_index   = segment_index;
                // Post Cdef Results
                svt_post_full_object(cdef_results_wrapper_ptr);
            }
        }
        svt_release_mutex(pcs_ptr->cdef_search_mutex);

        // Release Dlf Results
        svt_release_object(dlf_results_wrapper_ptr);
    }

    return NULL;
}
