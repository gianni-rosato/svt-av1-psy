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

#ifndef EbTransforms_h
#define EbTransforms_h

#include "EbDefinitions.h"
#include "EbCoefficients.h"
#include "EbInvTransforms.h"

#ifdef __cplusplus
extern "C" {
#endif

#include "EbEncDecProcess.h"

static const int8_t fwd_shift_4x4[3]   = {2, 0, 0};
static const int8_t fwd_shift_8x8[3]   = {2, -1, 0};
static const int8_t fwd_shift_16x16[3] = {2, -2, 0};
static const int8_t fwd_shift_32x32[3] = {2, -4, 0};
static const int8_t fwd_shift_64x64[3] = {0, -2, -2};
static const int8_t fwd_shift_4x8[3]   = {2, -1, 0};
static const int8_t fwd_shift_8x4[3]   = {2, -1, 0};
static const int8_t fwd_shift_8x16[3]  = {2, -2, 0};
static const int8_t fwd_shift_16x8[3]  = {2, -2, 0};
static const int8_t fwd_shift_16x32[3] = {2, -4, 0};
static const int8_t fwd_shift_32x16[3] = {2, -4, 0};
static const int8_t fwd_shift_32x64[3] = {0, -2, -2};
static const int8_t fwd_shift_64x32[3] = {2, -4, -2};
static const int8_t fwd_shift_4x16[3]  = {2, -1, 0};
static const int8_t fwd_shift_16x4[3]  = {2, -1, 0};
static const int8_t fwd_shift_8x32[3]  = {2, -2, 0};
static const int8_t fwd_shift_32x8[3]  = {2, -2, 0};
static const int8_t fwd_shift_16x64[3] = {0, -2, 0};
static const int8_t fwd_shift_64x16[3] = {2, -4, 0};

static const int8_t fwd_cos_bit_col[MAX_TXWH_IDX /*txw_idx*/][MAX_TXWH_IDX /*txh_idx*/] = {
    {13, 13, 13, 0, 0},
    {13, 13, 13, 12, 0},
    {13, 13, 13, 12, 13},
    {0, 13, 13, 12, 13},
    {0, 0, 13, 12, 13}};
static const int8_t fwd_cos_bit_row[MAX_TXWH_IDX /*txw_idx*/][MAX_TXWH_IDX /*txh_idx*/] = {
    {13, 13, 12, 0, 0},
    {13, 13, 13, 12, 0},
    {13, 13, 12, 13, 12},
    {0, 12, 13, 12, 11},
    {0, 0, 12, 11, 10}};



static const int8_t fdct4_range_mult2[4]       = {0, 2, 3, 3};
static const int8_t fdct8_range_mult2[6]       = {0, 2, 4, 5, 5, 5};
static const int8_t fdct16_range_mult2[8]      = {0, 2, 4, 6, 7, 7, 7, 7};
static const int8_t fdct32_range_mult2[10]     = {0, 2, 4, 6, 8, 9, 9, 9, 9, 9};
static const int8_t fdct64_range_mult2[12]     = {0, 2, 4, 6, 8, 10, 11, 11, 11, 11, 11, 11};
static const int8_t fadst4_range_mult2[7]      = {0, 2, 4, 3, 3, 3, 3};
static const int8_t fadst8_range_mult2[8]      = {0, 0, 1, 3, 3, 5, 5, 5};
static const int8_t fadst16_range_mult2[10]    = {0, 0, 1, 3, 3, 5, 5, 7, 7, 7};
static const int8_t fadst32_range_mult2[12]    = {0, 0, 1, 3, 3, 5, 5, 7, 7, 9, 9, 9};
static const int8_t max_fwd_range_mult2_col[5] = {3, 5, 7, 9, 11};
static const int8_t fidtx4_range_mult2[1]      = {1};
static const int8_t fidtx8_range_mult2[1]      = {2};
static const int8_t fidtx16_range_mult2[1]     = {3};
static const int8_t fidtx32_range_mult2[1]     = {4};
static const int8_t fidtx64_range_mult2[1]     = {5};

#define BLOCK_SIZES_ALL 22
static INLINE int is_rect_tx(TxSize tx_size) { return tx_size >= TX_SIZES; }
static INLINE int is_rect_tx_allowed_bsize(BlockSize bsize) {
    static const char lut[BLOCK_SIZES_ALL] = {
        0, // BLOCK_4X4
        1, // BLOCK_4X8
        1, // BLOCK_8X4
        0, // BLOCK_8X8
        1, // BLOCK_8X16
        1, // BLOCK_16X8
        0, // BLOCK_16X16
        1, // BLOCK_16X32
        1, // BLOCK_32X16
        0, // BLOCK_32X32
        1, // BLOCK_32X64
        1, // BLOCK_64X32
        0, // BLOCK_64X64
        0, // BLOCK_64X128
        0, // BLOCK_128X64
        0, // BLOCK_128X128
        1, // BLOCK_4X16
        1, // BLOCK_16X4
        1, // BLOCK_8X32
        1, // BLOCK_32X8
        1, // BLOCK_16X64
        1, // BLOCK_64X16
    };

    return lut[bsize];
}
static INLINE int is_rect_tx_allowed(/*const MacroBlockD *xd,*/
                                     const MbModeInfo *mbmi) {
    return is_rect_tx_allowed_bsize(mbmi->block_mi.sb_type) /*&&
            !xd->lossless[mbmi->segment_id]*/
        ;
}


////////////////////// QUANTIZATION//////////////
typedef struct QuantParam {
    int32_t      log_scale;
    TxSize       tx_size;
    const QmVal *qmatrix;
    const QmVal *iqmatrix;
} QuantParam;


static const uint32_t q_func[] = {26214, 23302, 20560, 18396, 16384, 14564};

extern EbErrorType av1_estimate_transform(int16_t *residual_buffer, uint32_t residual_stride,
                                          int32_t *coeff_buffer, uint32_t coeff_stride,
                                          TxSize transform_size, uint64_t *three_quad_energy,
                                          uint32_t bit_depth, TxType transform_type,
                                          PlaneType            component_type,
                                          EB_TRANS_COEFF_SHAPE trans_coeff_shape);

extern int32_t av1_quantize_inv_quantize(
        PictureControlSet *pcs_ptr, ModeDecisionContext *md_context, int32_t *coeff,
        const uint32_t coeff_stride, int32_t *quant_coeff, int32_t *recon_coeff, uint32_t qp,
        int32_t segmentation_qp_offset, uint32_t width, uint32_t height, TxSize txsize, uint16_t *eob,
        uint32_t *y_count_non_zero_coeffs, uint32_t component_type, uint32_t bit_increment,
        TxType tx_type, ModeDecisionCandidateBuffer *candidate_buffer, int16_t txb_skip_context,
        int16_t dc_sign_context, PredictionMode pred_mode, EbBool is_intra_bc, uint32_t lambda,EbBool is_encode_pass);


#ifdef __cplusplus
}
#endif

#endif // EbTransforms_h
