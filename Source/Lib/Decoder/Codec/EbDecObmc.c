/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

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
#include <stdlib.h>
#include <string.h>

#include "EbCodingUnit.h"
#include "EbSvtAv1Dec.h"
#include "EbDecHandle.h"
#include "EbDecProcessFrame.h"
#include "EbDecObmc.h"
#include "EbDecNbr.h"
#include "EbDecUtils.h"
#include "EbDecInverseQuantize.h"
#include "EbDecInterPrediction.h"

//This function is present in encoder also, but encoder structures & decoder structures are different.
static INLINE int dec_is_neighbor_overlappable(const ModeInfo_t *mbmi){
    // TODO: currently intrabc  is not supporting
    return mbmi->use_intrabc || mbmi->ref_frame[0] > INTRA_FRAME;
}

//static INLINE int dec_is_motion_variation_allowed_bsize(BlockSize bsize) {
//    return AOMMIN(block_size_wide[bsize], block_size_high[bsize]) >= 8;
//}

void av1_modify_neighbor_predictor_for_obmc(ModeInfo_t *mbmi) {
    mbmi->ref_frame[1] = NONE_FRAME;
    mbmi->inter_compound.type = COMPOUND_AVERAGE;
    return;
}

// HW does not support < 4x4 prediction. To limit the bandwidth requirement, if
// block-size of current plane is smaller than 8x8, always only blend with the
// left neighbor(s) (skip blending with the above side).
#define DISABLE_CHROMA_U8X8_OBMC 0  // 0: one-sided obmc; 1: disable

int av1_skip_u4x4_pred_in_obmc(BlockSize bsize, int32_t sub_x, int32_t sub_y, int dir) {

    assert(is_motion_variation_allowed_bsize(bsize));

    const BlockSize bsize_plane =
        get_plane_block_size(bsize, sub_x, sub_y);
    switch (bsize_plane) {
#if DISABLE_CHROMA_U8X8_OBMC
    case BLOCK_4X4:
    case BLOCK_8X4:
    case BLOCK_4X8: return 1; break;
#else
    case BLOCK_4X4:
    case BLOCK_8X4:
    case BLOCK_4X8: return dir == 0; break;
#endif
    default: return 0;
    }
}

// obmc_mask_N[overlap_position]
static const uint8_t obmc_mask_1[1] = { 64 };
static const uint8_t obmc_mask_2[2] = { 45, 64 };

static const uint8_t obmc_mask_4[4] = { 39, 50, 59, 64 };

static const uint8_t obmc_mask_8[8] = { 36, 42, 48, 53, 57, 61, 64, 64 };

static const uint8_t obmc_mask_16[16] = { 34, 37, 40, 43, 46, 49, 52, 54,
                                          56, 58, 60, 61, 64, 64, 64, 64 };

static const uint8_t obmc_mask_32[32] = { 33, 35, 36, 38, 40, 41, 43, 44,
                                          45, 47, 48, 50, 51, 52, 53, 55,
                                          56, 57, 58, 59, 60, 60, 61, 62,
                                          64, 64, 64, 64, 64, 64, 64, 64 };

static const uint8_t obmc_mask_64[64] = {
  33, 34, 35, 35, 36, 37, 38, 39, 40, 40, 41, 42, 43, 44, 44, 44,
  45, 46, 47, 47, 48, 49, 50, 51, 51, 51, 52, 52, 53, 54, 55, 56,
  56, 56, 57, 57, 58, 58, 59, 60, 60, 60, 60, 60, 61, 62, 62, 62,
  62, 62, 63, 63, 63, 63, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
};

const uint8_t *av1_get_obmc_mask(int length) {
    switch (length) {
    case 1: return obmc_mask_1;
    case 2: return obmc_mask_2;
    case 4: return obmc_mask_4;
    case 8: return obmc_mask_8;
    case 16: return obmc_mask_16;
    case 32: return obmc_mask_32;
    case 64: return obmc_mask_64;
    default: assert(0); return NULL;
    }
}

void aom_highbd_blend_a64_vmask_c(uint8_t *dst_8, uint32_t dst_stride,
    const uint8_t *src0_8, uint32_t src0_stride,
    const uint8_t *src1_8, uint32_t src1_stride,
    const uint8_t *mask, int w, int h, int bd)
{
    int i, j;
    uint16_t *dst = CONVERT_TO_SHORTPTR(dst_8);
    const uint16_t *src0 = CONVERT_TO_SHORTPTR(src0_8);
    const uint16_t *src1 = CONVERT_TO_SHORTPTR(src1_8);
    (void)bd;

    assert(IMPLIES(src0 == dst, src0_stride == dst_stride));
    assert(IMPLIES(src1 == dst, src1_stride == dst_stride));

    assert(h >= 1);
    assert(w >= 1);
    assert(IS_POWER_OF_TWO(h));
    assert(IS_POWER_OF_TWO(w));

    assert(bd == 8 || bd == 10 || bd == 12);

    for (i = 0; i < h; ++i) {
        const int m = mask[i];
        for (j = 0; j < w; ++j) {
            dst[i * dst_stride + j] = AOM_BLEND_A64(m, src0[i * src0_stride + j],
                src1[i * src1_stride + j]);
        }
    }
}

void aom_blend_a64_vmask_c(uint8_t *dst, uint32_t dst_stride,
    const uint8_t *src0, uint32_t src0_stride,
    const uint8_t *src1, uint32_t src1_stride,
    const uint8_t *mask, int w, int h)
{
    int i, j;

    assert(IMPLIES(src0 == dst, src0_stride == dst_stride));
    assert(IMPLIES(src1 == dst, src1_stride == dst_stride));

    assert(h >= 1);
    assert(w >= 1);
    assert(IS_POWER_OF_TWO(h));
    assert(IS_POWER_OF_TWO(w));

    for (i = 0; i < h; ++i) {
        const int m = mask[i];
        for (j = 0; j < w; ++j) {
            dst[i * dst_stride + j] = AOM_BLEND_A64(m, src0[i * src0_stride + j],
                src1[i * src1_stride + j]);
        }
    }
}

void aom_highbd_blend_a64_hmask_c(uint8_t *dst_8, uint32_t dst_stride,
    const uint8_t *src0_8, uint32_t src0_stride,
    const uint8_t *src1_8, uint32_t src1_stride,
    const uint8_t *mask, int w, int h, int bd)
{
    int i, j;
    uint16_t *dst = CONVERT_TO_SHORTPTR(dst_8);
    const uint16_t *src0 = CONVERT_TO_SHORTPTR(src0_8);
    const uint16_t *src1 = CONVERT_TO_SHORTPTR(src1_8);
    (void)bd;

    assert(IMPLIES(src0 == dst, src0_stride == dst_stride));
    assert(IMPLIES(src1 == dst, src1_stride == dst_stride));

    assert(h >= 1);
    assert(w >= 1);
    assert(IS_POWER_OF_TWO(h));
    assert(IS_POWER_OF_TWO(w));

    assert(bd == 8 || bd == 10 || bd == 12);

    for (i = 0; i < h; ++i) {
        for (j = 0; j < w; ++j) {
            dst[i * dst_stride + j] = AOM_BLEND_A64(
                mask[j], src0[i * src0_stride + j], src1[i * src1_stride + j]);
        }
    }
}

void aom_blend_a64_hmask_c(uint8_t *dst, uint32_t dst_stride,
    const uint8_t *src0, uint32_t src0_stride,
    const uint8_t *src1, uint32_t src1_stride,
    const uint8_t *mask, int w, int h) {
    int i, j;

    assert(IMPLIES(src0 == dst, src0_stride == dst_stride));
    assert(IMPLIES(src1 == dst, src1_stride == dst_stride));

    assert(h >= 1);
    assert(w >= 1);
    assert(IS_POWER_OF_TWO(h));
    assert(IS_POWER_OF_TWO(w));

    for (i = 0; i < h; ++i) {
        for (j = 0; j < w; ++j) {
            dst[i * dst_stride + j] = AOM_BLEND_A64(
                mask[j], src0[i * src0_stride + j], src1[i * src1_stride + j]);
        }
    }
}

static INLINE void build_obmc_inter_pred_above(EbDecHandle *dec_handle,
    PartitionInfo_t *pi, BlockSize bsize, int rel_mi_col, uint8_t above_mi_width,
    uint8_t *above_tmp_buf[MAX_MB_PLANE], int above_tmp_stride[MAX_MB_PLANE],
    uint8_t *curr_blk_recon_buf[MAX_MB_PLANE],
    int32_t curr_recon_stride[MAX_MB_PLANE], const int num_planes)
{

    EbPictureBufferDesc *recon_picture_buf = dec_handle->cur_pic_buf[0]->
        ps_pic_buf;
    const int is_hbd = (recon_picture_buf->bit_depth != EB_8BIT) ? 1 : 0;
    const int overlap =
        AOMMIN(block_size_high[bsize], block_size_high[BLOCK_64X64]) >> 1;
    uint8_t *above_buf;
    int32_t above_stride;
    uint8_t *tmp_recon_buf;
    int32_t tmp_recon_stride;

    for (int plane = 0; plane < num_planes; ++plane) {
        int32_t sub_x = (plane > 0) ? pi->subsampling_x : 0;
        int32_t sub_y = (plane > 0) ? pi->subsampling_y : 0;
        const int bw = (above_mi_width * MI_SIZE) >> sub_x;
        const int bh = overlap >> sub_y;

        if (av1_skip_u4x4_pred_in_obmc(bsize, sub_x, sub_y, 0)) continue;

        if (recon_picture_buf->bit_depth != EB_8BIT) {
            above_buf = (uint8_t *)((uint16_t *)above_tmp_buf[plane] +
                ((rel_mi_col*MI_SIZE) >> sub_x) +
                0/*No y-offset for obmc above pred*/);
            above_stride = above_tmp_stride[plane];
            tmp_recon_buf = (uint8_t *)((uint16_t *)curr_blk_recon_buf[plane] +
                ((rel_mi_col*MI_SIZE) >> sub_x) +
                0/*No y-offset for obmc above pred*/);
            tmp_recon_stride = curr_recon_stride[plane];

        }
        else {
            above_buf = above_tmp_buf[plane] + ((rel_mi_col*MI_SIZE) >> sub_x) +
                0/*No y-offset for obmc above pred*/;
            above_stride = above_tmp_stride[plane];
            tmp_recon_buf = curr_blk_recon_buf[plane] +
                ((rel_mi_col*MI_SIZE) >> sub_x) +
                0/*No y-offset for obmc above pred*/;
            tmp_recon_stride = curr_recon_stride[plane];
        }

        const uint8_t *const mask = av1_get_obmc_mask(bh);

        if (is_hbd)
            aom_highbd_blend_a64_vmask_c(CONVERT_TO_BYTEPTR(tmp_recon_buf),
                tmp_recon_stride, CONVERT_TO_BYTEPTR(tmp_recon_buf),
                tmp_recon_stride, CONVERT_TO_BYTEPTR(above_buf), above_stride,
                mask, bw, bh, recon_picture_buf->bit_depth);

        else
            aom_blend_a64_vmask_c(tmp_recon_buf, tmp_recon_stride, tmp_recon_buf,
                tmp_recon_stride, above_buf, above_stride, mask, bw, bh);
    }

}


static INLINE void build_obmc_inter_pred_left(EbDecHandle *dec_handle,
    PartitionInfo_t *pi, BlockSize bsize,
    int rel_mi_row, uint8_t left_mi_height,
    uint8_t *left_tmp_buf[MAX_MB_PLANE], int left_tmp_stride[MAX_MB_PLANE],
    uint8_t *curr_blk_recon_buf[MAX_MB_PLANE],
    int32_t curr_recon_stride[MAX_MB_PLANE], const int num_planes)
{
    EbPictureBufferDesc *recon_picture_buf = dec_handle->cur_pic_buf[0]->
        ps_pic_buf;
    const int is_hbd = (recon_picture_buf->bit_depth != EB_8BIT) ? 1 : 0;
    const int overlap =
        AOMMIN(block_size_wide[bsize], block_size_wide[BLOCK_64X64]) >> 1;

    uint8_t *left_buf;
    int32_t left_stride;
    uint8_t *tmp_recon_buf;
    int32_t tmp_recon_stride;
    for (int plane = 0; plane < num_planes; ++plane) {

        int32_t sub_x = (plane > 0) ? pi->subsampling_x : 0;
        int32_t sub_y = (plane > 0) ? pi->subsampling_y : 0;
        const int bw = overlap >> sub_x;
        const int bh = (left_mi_height * MI_SIZE) >> sub_y;


        if (av1_skip_u4x4_pred_in_obmc(bsize, sub_x, sub_y, 1)) continue;

        if (recon_picture_buf->bit_depth != EB_8BIT) {
            left_buf = (uint8_t *)((uint16_t *)left_tmp_buf[plane] +
                ((MI_SIZE*rel_mi_row*left_tmp_stride[plane]) >> sub_y) +
                0/*No x offst for left obmc pred*/);
            left_stride = left_tmp_stride[plane];

            tmp_recon_buf = (uint8_t *)((uint16_t *)curr_blk_recon_buf[plane] +
                ((MI_SIZE*rel_mi_row*curr_recon_stride[plane]) >> sub_y) +
                0/*No y-offset for obmc above pred*/);
            tmp_recon_stride = curr_recon_stride[plane];
        }
        else {
            left_buf = left_tmp_buf[plane] +
                ((MI_SIZE*rel_mi_row*left_tmp_stride[plane]) >> sub_y) +
                0/*No x offst for left obmc pred*/;
            left_stride = left_tmp_stride[plane];

            tmp_recon_buf = curr_blk_recon_buf[plane] +
                ((MI_SIZE*rel_mi_row*curr_recon_stride[plane]) >> sub_y) +
                0/*No y-offset for obmc above pred*/;
            tmp_recon_stride = curr_recon_stride[plane];
        }

        const uint8_t *const mask = av1_get_obmc_mask(bw);

        if (is_hbd)
            aom_highbd_blend_a64_hmask_c(CONVERT_TO_BYTEPTR(tmp_recon_buf),
                tmp_recon_stride, CONVERT_TO_BYTEPTR(tmp_recon_buf),
                tmp_recon_stride,CONVERT_TO_BYTEPTR(left_buf), left_stride,
                mask, bw, bh, recon_picture_buf->bit_depth);
        else
            aom_blend_a64_hmask_c(tmp_recon_buf, tmp_recon_stride, tmp_recon_buf,
                tmp_recon_stride, left_buf, left_stride, mask, bw, bh);
    }

}

static INLINE void dec_build_prediction_by_above_pred(EbDecHandle *dec_handle,
    PartitionInfo_t *backup_pi, BlockSize bsize, int bw4, int mi_row, int mi_col,
    int rel_mi_col, uint8_t above_mi_width, ModeInfo_t *above_mbmi,
    uint8_t *tmp_buf[MAX_MB_PLANE], int tmp_stride[MAX_MB_PLANE],
    const int num_planes)
{
    EbPictureBufferDesc *recon_picture_buf = dec_handle->cur_pic_buf[0]->
        ps_pic_buf;
    const int above_mi_col = mi_col + rel_mi_col;
    int mi_x, mi_y;
    uint8_t *tmp_recon_buf;
    int32_t tmp_recon_stride;
    ModeInfo_t bakup_abv_mbmi = *above_mbmi;
    backup_pi->mi = &bakup_abv_mbmi;
    av1_modify_neighbor_predictor_for_obmc(backup_pi->mi);

    backup_pi->mb_to_left_edge = 8 * MI_SIZE * (-above_mi_col);
    backup_pi->mb_to_right_edge += (bw4 - rel_mi_col - above_mi_width)
                                    * MI_SIZE * 8;

    mi_x = above_mi_col << MI_SIZE_LOG2;
    mi_y = mi_row << MI_SIZE_LOG2;
    backup_pi->wpx[0] = (above_mi_width * MI_SIZE);

    for (int plane = 0; plane < num_planes; ++plane) {
        uint8_t sub_x = (plane > 0) ? backup_pi->subsampling_x : 0;
        uint8_t sub_y = (plane > 0) ? backup_pi->subsampling_y : 0;

        if (recon_picture_buf->bit_depth != EB_8BIT) {
            tmp_recon_buf = (uint8_t *)((uint16_t *)tmp_buf[plane]+
                ((rel_mi_col*MI_SIZE) >> sub_x) +
                0/*No y-offset for obmc above pred*/);

            tmp_recon_stride = tmp_stride[plane];
        }
        else {
            tmp_recon_buf = tmp_buf[plane] + ((rel_mi_col*MI_SIZE) >> sub_x) +
                0/*No y-offset for obmc above pred*/;
            tmp_recon_stride = tmp_stride[plane];
        }

        if (av1_skip_u4x4_pred_in_obmc(bsize, sub_x, sub_y, 0)) continue;
        svtav1_predict_inter_block_plane(dec_handle, backup_pi, plane,
            1/*obmc*/, mi_x, mi_y, (void *)tmp_recon_buf, tmp_recon_stride,
            0/*some_use_intra*/, recon_picture_buf->bit_depth);


    }
}

static void dec_build_prediction_by_above_preds(EbDecHandle *dec_handle,
    PartitionInfo_t *pi, int mi_row, int mi_col,
    uint8_t *above_dst_buf[MAX_MB_PLANE],
    int above_dst_stride[MAX_MB_PLANE])
{
    if (!pi->up_available) return;
    PartitionInfo_t backup_pi = *pi;

    // Adjust mb_to_bottom_edge to have the correct value for the OBMC
    // prediction block. This is half the height of the original block,
    // except for 128-wide blocks, where we only use a height of 32.
    BlockSize bsize = pi->mi->sb_type;
    int bh4 =  mi_size_high[bsize];
    int bw4 =  mi_size_wide[bsize];
    int currblock_height = bh4 * MI_SIZE;
    int pred_height = AOMMIN(currblock_height / 2, MAX_OBMC_LEN);
    backup_pi.mb_to_bottom_edge += (currblock_height - pred_height) * 8;

    int bh = clamp(block_size_high[bsize] >>  1, 4,
        block_size_high[BLOCK_64X64] >> 1);
    backup_pi.hpx[0] = bh;

    EbColorConfig *cc = &dec_handle->seq_header.color_config;
    FrameHeader *frame_info = &dec_handle->frame_header;
    const int num_planes = av1_num_planes(cc);
    int nb_count = 0;
    uint8_t mi_step;
    int nb_max = max_neighbor_obmc[mi_size_wide_log2[bsize]];

    const int end_col = AOMMIN(mi_col + bw4, (int)frame_info->mi_cols);

    //Calculating buffers for current block i.e getting recon_buffer for blending
    uint8_t *curr_blk_recon_buf[MAX_MB_PLANE];
    int32_t curr_recon_stride[MAX_MB_PLANE];
    int32_t sub_x, sub_y;
    EbPictureBufferDesc *recon_picture_buf = dec_handle->cur_pic_buf[0]->
        ps_pic_buf;
    for (int plane = 0; plane < num_planes; ++plane) {
        sub_x = (plane > 0) ? pi->subsampling_x : 0;
        sub_y = (plane > 0) ? pi->subsampling_y : 0;

        derive_blk_pointers(recon_picture_buf, plane,
            mi_col*MI_SIZE >> sub_x, mi_row*MI_SIZE >> sub_y,
            (void *)&curr_blk_recon_buf[plane], &curr_recon_stride[plane],
            sub_x, sub_y);
    }

    for (int above_mi_col = mi_col; above_mi_col < end_col && nb_count < nb_max;
        above_mi_col += mi_step) {
        ModeInfo_t *above_mi = get_cur_mode_info(dec_handle, mi_row-1,
            above_mi_col, NULL);

        mi_step =
            AOMMIN(mi_size_wide[above_mi->sb_type], mi_size_wide[BLOCK_64X64]);

        // If we're considering a block with width 4, it should be treated as
        // half of a pair of blocks with chroma information in the second. Move
        // above_mi_col back to the start of the pair if needed, set above_mbmi
        // to point at the block with chroma information, and set mi_step to 2 to
        // step over the entire pair at the end of the iteration.
        if (mi_step == 1) {
            above_mi = get_cur_mode_info(dec_handle, mi_row - 1,
                above_mi_col | 1, NULL);
            mi_step = 2;
        }
        if (dec_is_neighbor_overlappable(above_mi)) {
            ++nb_count;
            /*OBMC above prediction*/
            dec_build_prediction_by_above_pred(dec_handle, &backup_pi, bsize,
                bw4, mi_row, mi_col, above_mi_col - mi_col,
                AOMMIN((uint8_t)bw4, mi_step), above_mi, above_dst_buf,
                above_dst_stride, num_planes);
            /*OBMC blending for above prediction*/
            build_obmc_inter_pred_above(dec_handle, &backup_pi, bsize,
                above_mi_col - mi_col, AOMMIN((uint8_t)bw4, mi_step),
                above_dst_buf, above_dst_stride, curr_blk_recon_buf,
                curr_recon_stride, num_planes);
        }
    }

}


static INLINE void dec_build_prediction_by_left_pred(EbDecHandle *dec_handle,
    PartitionInfo_t *backup_pi, BlockSize bsize, int bh4, int mi_row, int mi_col,
    int rel_mi_row, uint8_t left_mi_height, ModeInfo_t *left_mbmi,
    uint8_t *tmp_buf[MAX_MB_PLANE], int tmp_stride[MAX_MB_PLANE],
    const int num_planes)
{
    EbPictureBufferDesc *recon_picture_buf = dec_handle->cur_pic_buf[0]->
        ps_pic_buf;
    const int left_mi_row = mi_row + rel_mi_row;
    int mi_x, mi_y;
    uint8_t *tmp_recon_buf;
    int32_t tmp_recon_stride;
    ModeInfo_t bakup_left_mbmi = *left_mbmi;
    backup_pi->mi = &bakup_left_mbmi;
    av1_modify_neighbor_predictor_for_obmc(backup_pi->mi);

    backup_pi->mb_to_top_edge = 8 * MI_SIZE * (-left_mi_row);
    backup_pi->mb_to_bottom_edge +=
        (bh4 - rel_mi_row - left_mi_height) * MI_SIZE * 8;

    mi_x = mi_col << MI_SIZE_LOG2;
    mi_y = left_mi_row << MI_SIZE_LOG2;
    backup_pi->hpx[0] = (left_mi_height << MI_SIZE_LOG2);

    for (int plane = 0; plane < num_planes; ++plane) {
        int32_t sub_x = (plane > 0) ? backup_pi->subsampling_x : 0;
        int32_t sub_y = (plane > 0) ? backup_pi->subsampling_y : 0;

        if (recon_picture_buf->bit_depth != EB_8BIT) {
            tmp_recon_buf = (uint8_t *)((uint16_t *)tmp_buf[plane]+
            ((MI_SIZE*rel_mi_row*tmp_stride[plane]) >> sub_y) +
                0/*No x offst for left obmc pred*/);

            tmp_recon_stride = tmp_stride[plane];
        }
        else {
            tmp_recon_buf = tmp_buf[plane] +
                ((MI_SIZE*rel_mi_row*tmp_stride[plane]) >> sub_y) +
                0/*No x offst for left obmc pred*/;
            tmp_recon_stride = tmp_stride[plane];
        }

        if (av1_skip_u4x4_pred_in_obmc(bsize, sub_x, sub_y, 1)) continue;
       // dec_build_inter_predictors(ctxt->cm, pi, j, &backup_mbmi, 1, bw, bh, mi_x,
       //                            mi_y);
        svtav1_predict_inter_block_plane(dec_handle, backup_pi, plane,
            1/*obmc*/, mi_x, mi_y, (void *)tmp_recon_buf, tmp_recon_stride,
            0/*some_use_intra*/, recon_picture_buf->bit_depth);

    }
}

static void dec_build_prediction_by_left_preds(EbDecHandle *dec_handle,
    PartitionInfo_t *pi, int mi_row, int mi_col,
    uint8_t *left_dst_buf[MAX_MB_PLANE],
    int left_dst_stride[MAX_MB_PLANE])
{
    if (!pi->left_available) return;
    PartitionInfo_t backup_pi = *pi;

    // Adjust mb_to_right_edge to have the correct value for the OBMC
    // prediction block. This is half the width of the original block,
    // except for 128-wide blocks, where we only use a width of 32.
    BlockSize bsize = pi->mi->sb_type;
    int bh4 = mi_size_high[bsize];
    int bw4 = mi_size_wide[bsize];
    int currblock_width = bw4 * MI_SIZE;
    int pred_width = AOMMIN(currblock_width / 2, MAX_OBMC_LEN);
    backup_pi.mb_to_right_edge += (currblock_width - pred_width) * 8;

    int bw = clamp(block_size_wide[bsize] >>  1, 4,
        block_size_wide[BLOCK_64X64] >> 1);
    backup_pi.wpx[0] = bw;

    EbColorConfig *cc = &dec_handle->seq_header.color_config;
    FrameHeader *frame_info = &dec_handle->frame_header;
    const int num_planes = av1_num_planes(cc);
    int nb_count = 0;
    uint8_t mi_step;
    int nb_max = max_neighbor_obmc[mi_size_high_log2[bsize]];

    const int end_row = AOMMIN(mi_row + bh4, (int)frame_info->mi_rows);

    //Calculating buffers for current block i.e getting recon_buffer
    uint8_t *curr_blk_recon_buf[MAX_MB_PLANE];
    int32_t curr_recon_stride[MAX_MB_PLANE];
    int32_t sub_x, sub_y;
    EbPictureBufferDesc *recon_picture_buf = dec_handle->cur_pic_buf[0]->
        ps_pic_buf;
    for (int plane = 0; plane < num_planes; ++plane) {
        sub_x = (plane > 0) ? pi->subsampling_x : 0;
        sub_y = (plane > 0) ? pi->subsampling_y : 0;

        derive_blk_pointers(recon_picture_buf, plane,
            mi_col*MI_SIZE >> sub_x, mi_row*MI_SIZE >> sub_y,
            (void *)&curr_blk_recon_buf[plane], &curr_recon_stride[plane], sub_x, sub_y);
    }


    for (int left_mi_row = mi_row; left_mi_row < end_row && nb_count < nb_max;
        left_mi_row += mi_step) {
        ModeInfo_t *left_mi =
            get_cur_mode_info(dec_handle, left_mi_row, mi_col-1, NULL);
        mi_step =
            AOMMIN(mi_size_high[left_mi->sb_type], mi_size_high[BLOCK_64X64]);
        if (mi_step == 1) {
            left_mi = get_cur_mode_info(dec_handle, left_mi_row | 1,
                 mi_col - 1, NULL);
            mi_step = 2;
        }
        if (dec_is_neighbor_overlappable(left_mi)) {
            ++nb_count;
            /*OBMC left prediction*/
            dec_build_prediction_by_left_pred(dec_handle, &backup_pi, bsize,
                bh4, mi_row, mi_col, left_mi_row - mi_row,
                AOMMIN((uint8_t)bh4, mi_step), left_mi, left_dst_buf,
                left_dst_stride, num_planes);
            /*OBMC blending for left prediction*/
            build_obmc_inter_pred_left(dec_handle, &backup_pi, bsize,
                left_mi_row - mi_row, AOMMIN((uint8_t)bh4, mi_step),
                left_dst_buf, left_dst_stride, curr_blk_recon_buf,
                curr_recon_stride, num_planes);
        }
    }
}


void dec_build_obmc_inter_predictors_sb( EbDecHandle *dec_handle,
    PartitionInfo_t *pi, int mi_row, int mi_col)
{
    DecModCtxt *dec_mod_ctxt = (DecModCtxt *)dec_handle->pv_dec_mod_ctxt;
    uint8_t *dst_buf[MAX_MB_PLANE];
    dec_mod_ctxt->obmc_ctx.dst_stride[AOM_PLANE_Y] = MAX_SB_SIZE;
    dec_mod_ctxt->obmc_ctx.dst_stride[AOM_PLANE_U] = MAX_SB_SIZE;
    dec_mod_ctxt->obmc_ctx.dst_stride[AOM_PLANE_V] = MAX_SB_SIZE;
    dst_buf[0] = dec_mod_ctxt->obmc_ctx.tmp_obmc_bufs[AOM_PLANE_Y];
    dst_buf[1] = dec_mod_ctxt->obmc_ctx.tmp_obmc_bufs[AOM_PLANE_U];
    dst_buf[2] = dec_mod_ctxt->obmc_ctx.tmp_obmc_bufs[AOM_PLANE_V];
    /*OBMC above prediction followed by Blending happen in below fun call*/
    dec_build_prediction_by_above_preds(dec_handle, pi, mi_row, mi_col,
        dst_buf, dec_mod_ctxt->obmc_ctx.dst_stride);

    /*OBMC left prediction followed by Blending happen in below fun call*/
    dec_build_prediction_by_left_preds(dec_handle, pi, mi_row, mi_col,
        dst_buf, dec_mod_ctxt->obmc_ctx.dst_stride);
}

