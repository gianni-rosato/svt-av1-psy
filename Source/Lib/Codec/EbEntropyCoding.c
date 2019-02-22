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

#include "EbEntropyCoding.h"
#include "EbEntropyCodingUtil.h"
#include "EbUtility.h"
#include "EbCodingUnit.h"
#include "EbErrorCodes.h"
#include "EbTransforms.h"
#include "EbEntropyCodingProcess.h"

#include "aom_dsp_rtcd.h"

#define S32 32*32
#define S16 16*16
#define S8  8*8
#define S4  4*4

int32_t av1_loop_restoration_corners_in_sb(Av1Common *cm, int32_t plane,
    int32_t mi_row, int32_t mi_col, BlockSize bsize,
    int32_t *rcol0, int32_t *rcol1, int32_t *rrow0,
    int32_t *rrow1, int32_t *tile_tl_idx);

#define CHAR_BIT      8         /* number of bits in a char */
#if ADD_DELTA_QP_SUPPORT
#define OD_CLZ0 (1)
#define OD_CLZ(x) (-get_msb(x))
#define OD_ILOG_NZ(x) (OD_CLZ0 - OD_CLZ(x))
#endif

static INLINE int32_t is_comp_ref_allowed(BlockSize bsize) {
    return AOMMIN(block_size_wide[bsize], block_size_high[bsize]) >= 8;
}
int32_t av1_loop_restoration_corners_in_sb(Av1Common *cm, int32_t plane,
    int32_t mi_row, int32_t mi_col, BlockSize bsize,
    int32_t *rcol0, int32_t *rcol1, int32_t *rrow0,
    int32_t *rrow1, int32_t *tile_tl_idx);

/*****************************
* Enums
*****************************/

enum COEFF_SCAN_TYPE
{
    SCAN_ZIGZAG = 0,      // zigzag scan
    SCAN_HOR,             // first scan is horizontal
    SCAN_VER,             // first scan is vertical
    SCAN_DIAG             // up-right diagonal scan
};

extern void av1_set_ref_frame(MvReferenceFrame *rf,
    int8_t ref_frame_type);

/************************************************
* CABAC Encoder Constructor
************************************************/
void CabacCtor(
    CabacEncodeContext_t *cabacEncContextPtr)
{
    EB_MEMSET(cabacEncContextPtr, 0, sizeof(CabacEncodeContext_t));

    return;
}

static INLINE int32_t does_level_match(int32_t width, int32_t height, double fps,
    int32_t lvl_width, int32_t lvl_height,
    double lvl_fps, int32_t lvl_dim_mult) {
    const int64_t lvl_luma_pels = lvl_width * lvl_height;
    const double lvl_display_sample_rate = lvl_luma_pels * lvl_fps;
    const int64_t luma_pels = width * height;
    const double display_sample_rate = luma_pels * fps;
    return luma_pels <= lvl_luma_pels &&
        display_sample_rate <= lvl_display_sample_rate &&
        width <= lvl_width * lvl_dim_mult &&
        height <= lvl_height * lvl_dim_mult;
}

static void SetBitstreamLevelTier(SequenceControlSet_t *scsPtr) {
    // TODO(any): This is a placeholder function that only addresses dimensions
    // and max display sample rates.
    // Need to add checks for max bit rate, max decoded luma sample rate, header
    // rate, etc. that are not covered by this function.

    BitstreamLevel bl = { 9, 3 };
    if (does_level_match(scsPtr->luma_width, scsPtr->luma_height, (scsPtr->frame_rate >> 16), 512,
        288, 30.0, 4)) {
        bl.major = 2;
        bl.minor = 0;
    }
    else if (does_level_match(scsPtr->luma_width, scsPtr->luma_height, (scsPtr->frame_rate >> 16),
        704, 396, 30.0, 4)) {
        bl.major = 2;
        bl.minor = 1;
    }
    else if (does_level_match(scsPtr->luma_width, scsPtr->luma_height, (scsPtr->frame_rate >> 16),
        1088, 612, 30.0, 4)) {
        bl.major = 3;
        bl.minor = 0;
    }
    else if (does_level_match(scsPtr->luma_width, scsPtr->luma_height, (scsPtr->frame_rate >> 16),
        1376, 774, 30.0, 4)) {
        bl.major = 3;
        bl.minor = 1;
    }
    else if (does_level_match(scsPtr->luma_width, scsPtr->luma_height, (scsPtr->frame_rate >> 16),
        2048, 1152, 30.0, 3)) {
        bl.major = 4;
        bl.minor = 0;
    }
    else if (does_level_match(scsPtr->luma_width, scsPtr->luma_height, (scsPtr->frame_rate >> 16),
        2048, 1152, 60.0, 3)) {
        bl.major = 4;
        bl.minor = 1;
    }
    else if (does_level_match(scsPtr->luma_width, scsPtr->luma_height, (scsPtr->frame_rate >> 16),
        4096, 2176, 30.0, 2)) {
        bl.major = 5;
        bl.minor = 0;
    }
    else if (does_level_match(scsPtr->luma_width, scsPtr->luma_height, (scsPtr->frame_rate >> 16),
        4096, 2176, 60.0, 2)) {
        bl.major = 5;
        bl.minor = 1;
    }
    else if (does_level_match(scsPtr->luma_width, scsPtr->luma_height, (scsPtr->frame_rate >> 16),
        4096, 2176, 120.0, 2)) {
        bl.major = 5;
        bl.minor = 2;
    }
    else if (does_level_match(scsPtr->luma_width, scsPtr->luma_height, (scsPtr->frame_rate >> 16),
        8192, 4352, 30.0, 2)) {
        bl.major = 6;
        bl.minor = 0;
    }
    else if (does_level_match(scsPtr->luma_width, scsPtr->luma_height, (scsPtr->frame_rate >> 16),
        8192, 4352, 60.0, 2)) {
        bl.major = 6;
        bl.minor = 1;
    }
    else if (does_level_match(scsPtr->luma_width, scsPtr->luma_height, (scsPtr->frame_rate >> 16),
        8192, 4352, 120.0, 2)) {
        bl.major = 6;
        bl.minor = 2;
    }
    else if (does_level_match(scsPtr->luma_width, scsPtr->luma_height, (scsPtr->frame_rate >> 16),
        16384, 8704, 30.0, 2)) {
        bl.major = 7;
        bl.minor = 0;
    }
    else if (does_level_match(scsPtr->luma_width, scsPtr->luma_height, (scsPtr->frame_rate >> 16),
        16384, 8704, 60.0, 2)) {
        bl.major = 7;
        bl.minor = 1;
    }
    else if (does_level_match(scsPtr->luma_width, scsPtr->luma_height, (scsPtr->frame_rate >> 16),
        16384, 8704, 120.0, 2)) {
        bl.major = 7;
        bl.minor = 2;
    }
    for (int32_t i = 0; i < MAX_NUM_OPERATING_POINTS; ++i) {
        scsPtr->level[i] = bl;
        scsPtr->tier[i] = 0;  // setting main tier by default

    }
}





const uint8_t KEobOffsetBits[12] = { 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

static void WriteGolomb(aom_writer *w, int32_t level) {
    int32_t x = level + 1;
    int32_t i = x;
    int32_t length = 0;

    while (i) {
        i >>= 1;
        ++length;
    }
    assert(length > 0);

    for (i = 0; i < length - 1; ++i) aom_write_bit(w, 0);

    for (i = length - 1; i >= 0; --i) aom_write_bit(w, (x >> i) & 0x01);
}

static const uint8_t EobToPosSmall[33] = {
    0, 1, 2,                                        // 0-2
    3, 3,                                           // 3-4
    4, 4, 4, 4,                                     // 5-8
    5, 5, 5, 5, 5, 5, 5, 5,                         // 9-16
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6  // 17-32
};

static const uint8_t EobToPosLarge[17] = {
    6,                               // place holder
    7,                               // 33-64
    8, 8,                           // 65-128
    9, 9, 9, 9,                   // 129-256
    10, 10, 10, 10, 10, 10, 10, 10,  // 257-512
    11                               // 513-
};
const int16_t KEobGroupStart[12] = { 0,  1,  2,  3,   5,   9,
                                        17, 33, 65, 129, 257, 513 };

static INLINE int16_t GetEobPosToken(const int16_t eob, int16_t *const extra) {
    int16_t t;

    if (eob < 33) {
        t = EobToPosSmall[eob];
    }
    else {
        const int16_t e = MIN((eob - 1) >> 5, 16);
        t = EobToPosLarge[e];
    }

    *extra = eob - KEobGroupStart[t];

    return t;
}

static INLINE uint8_t *SetLevels(uint8_t *const levelsBuf, const int32_t width) {
    return levelsBuf + TX_PAD_TOP * (width + TX_PAD_HOR);
}

static INLINE void av1TxbInitLevels(
    int32_t         *coeffBufferPtr,
    const uint32_t    coeffStride,
    const int16_t    width,
    const int16_t    height,
    uint8_t *const    levels) {

    const int16_t stride = width + TX_PAD_HOR;
    uint8_t *ls = levels;

    memset(levels - TX_PAD_TOP * stride, 0,
        sizeof(*levels) * TX_PAD_TOP * stride);
    memset(levels + stride * height, 0,
        sizeof(*levels) * (TX_PAD_BOTTOM * stride + TX_PAD_END));

    for (int16_t i = 0; i < height; i++) {
        for (int16_t j = 0; j < width; j++) {
            *ls++ = (uint8_t)CLIP3(0, INT8_MAX, abs(coeffBufferPtr[i * coeffStride + j]));
        }
        for (int16_t j = 0; j < TX_PAD_HOR; j++) {
            *ls++ = 0;
        }
    }
}

/************************************************************************************************/
// blockd.h

static INLINE int32_t get_txb_wide(TxSize tx_size) {
    tx_size = av1_get_adjusted_tx_size(tx_size);
    return tx_size_wide[tx_size];
}
static INLINE int32_t get_txb_high(TxSize tx_size) {
    tx_size = av1_get_adjusted_tx_size(tx_size);
    return tx_size_high[tx_size];
}
static INLINE int32_t get_txb_bwl(TxSize tx_size) {
    tx_size = av1_get_adjusted_tx_size(tx_size);
    return tx_size_wide_log2[tx_size];
}

static INLINE int16_t GetBrCtx(
    const uint8_t *const levels,
    const int16_t c,  // raster order
    const int16_t bwl,
    const TxType tx_type

) {
    const int16_t row = c >> bwl;
    const int16_t col = c - (row << bwl);
    const int16_t stride = (1 << bwl) + TX_PAD_HOR;
    const TX_CLASS tx_class = tx_type_to_class[tx_type];
    const int16_t pos = row * stride + col;
    int16_t mag = levels[pos + 1];
    mag += levels[pos + stride];
    switch (tx_class) {
    case TX_CLASS_2D:
        mag += levels[pos + stride + 1];
        mag = MIN((mag + 1) >> 1, 6);
        if (c == 0) return mag;
        if ((row < 2) && (col < 2)) return mag + 7;
        break;

    case TX_CLASS_HORIZ:
        mag += levels[pos + 2];
        mag = AOMMIN((mag + 1) >> 1, 6);
        if (c == 0) return mag;
        if (col == 0) return mag + 7;
        break;
    case TX_CLASS_VERT:
        mag += levels[pos + (stride << 1)];
        mag = AOMMIN((mag + 1) >> 1, 6);
        if (c == 0) return mag;
        if (row == 0) return mag + 7;
        break;


    default: break;
    }

    return mag + 14;
}

void GetTxbCtx(
    const int32_t               plane,
    NeighborArrayUnit_t     *dcSignLevelCoeffNeighborArray,
    uint32_t                  cu_origin_x,
    uint32_t                  cu_origin_y,
    const BlockSize        plane_bsize,
    const TxSize           tx_size,
    int16_t *const           txb_skip_ctx,
    int16_t *const           dc_sign_ctx) {

    uint32_t dcSignLevelCoeffLeftNeighborIndex = GetNeighborArrayUnitLeftIndex(
        dcSignLevelCoeffNeighborArray,
        cu_origin_y);
    uint32_t dcSignLevelCoeffTopNeighborIndex = GetNeighborArrayUnitTopIndex(
        dcSignLevelCoeffNeighborArray,
        cu_origin_x);

#define MAX_TX_SIZE_UNIT 16
    static const int8_t signs[3] = { 0, -1, 1 };
    const int32_t txb_w_unit = tx_size_wide_unit[tx_size];
    const int32_t txb_h_unit = tx_size_high_unit[tx_size];
    int16_t dc_sign = 0;
    uint16_t k = 0;

    uint8_t sign;

#if TILES
    if (dcSignLevelCoeffNeighborArray->topArray[dcSignLevelCoeffTopNeighborIndex] != INVALID_NEIGHBOR_DATA){
#else
    if (cu_origin_y != 0) {//dcSignLevelCoeffNeighborArray->topArray[dcSignLevelCoeffTopNeighborIndex] != INVALID_NEIGHBOR_DATA){ AMIR
#endif
        do {
            sign = ((uint8_t)dcSignLevelCoeffNeighborArray->topArray[k + dcSignLevelCoeffTopNeighborIndex] >> COEFF_CONTEXT_BITS);
            assert(sign <= 2);
            dc_sign += signs[sign];
        } while (++k < txb_w_unit);
    }

#if TILES
    if (dcSignLevelCoeffNeighborArray->leftArray[dcSignLevelCoeffLeftNeighborIndex] != INVALID_NEIGHBOR_DATA){ 
#else
    if (cu_origin_x != 0) {// dcSignLevelCoeffNeighborArray->leftArray[dcSignLevelCoeffLeftNeighborIndex] != INVALID_NEIGHBOR_DATA){ AMIR
#endif
        k = 0;
        do {
            sign = ((uint8_t)dcSignLevelCoeffNeighborArray->leftArray[k + dcSignLevelCoeffLeftNeighborIndex] >> COEFF_CONTEXT_BITS);
            assert(sign <= 2);
            dc_sign += signs[sign];
        } while (++k < txb_h_unit);
    }

    if (dc_sign > 0)
        *dc_sign_ctx = 2;
    else if (dc_sign < 0)
        *dc_sign_ctx = 1;
    else
        *dc_sign_ctx = 0;

    if (plane == 0) {
        if (plane_bsize == txsize_to_bsize[tx_size]) {
            *txb_skip_ctx = 0;
        }
        else {
            static const uint8_t skip_contexts[5][5] = { { 1, 2, 2, 2, 3 },
            { 1, 4, 4, 4, 5 },
            { 1, 4, 4, 4, 5 },
            { 1, 4, 4, 4, 5 },
            { 1, 4, 4, 4, 6 } };
            int32_t top = 0;
            int32_t left = 0;

            k = 0;
#if TILES
            if (dcSignLevelCoeffNeighborArray->topArray[dcSignLevelCoeffTopNeighborIndex] != INVALID_NEIGHBOR_DATA) {
#else
            if (cu_origin_y != 0) {
#endif
                do {
                    top |= (int32_t)(dcSignLevelCoeffNeighborArray->topArray[k + dcSignLevelCoeffTopNeighborIndex]);
                } while (++k < txb_w_unit);
            }
            top &= COEFF_CONTEXT_MASK;

#if TILES
            if (dcSignLevelCoeffNeighborArray->leftArray[dcSignLevelCoeffLeftNeighborIndex] != INVALID_NEIGHBOR_DATA) {
#else
            if (cu_origin_x != 0) {
#endif
                k = 0;
                do {
                    left |= (int32_t)(dcSignLevelCoeffNeighborArray->leftArray[k + dcSignLevelCoeffLeftNeighborIndex]);
                } while (++k < txb_h_unit);
            }
            left &= COEFF_CONTEXT_MASK;
            //do {
            //    top |= a[k];
            //} while (++k < txb_w_unit);
            //top &= COEFF_CONTEXT_MASK;

            //k = 0;
            //do {
            //    left |= l[k];
            //} while (++k < txb_h_unit);
            //left &= COEFF_CONTEXT_MASK;
            const int32_t max = AOMMIN(top | left, 4);
            const int32_t min = AOMMIN(AOMMIN(top, left), 4);

            *txb_skip_ctx = skip_contexts[min][max];
        }
    }
    else {
        //const int32_t ctx_base = get_entropy_context(tx_size, a, l);
        int16_t ctx_base_left = 0;
        int16_t ctx_base_top = 0;

#if TILES
        if (dcSignLevelCoeffNeighborArray->topArray[dcSignLevelCoeffTopNeighborIndex] != INVALID_NEIGHBOR_DATA) {
#else
        if (cu_origin_y != 0) {
#endif
            k = 0;
            do {
                ctx_base_top += (dcSignLevelCoeffNeighborArray->topArray[k + dcSignLevelCoeffTopNeighborIndex] != 0);
            } while (++k < txb_w_unit);
        }
#if TILES
        if (dcSignLevelCoeffNeighborArray->leftArray[dcSignLevelCoeffLeftNeighborIndex] != INVALID_NEIGHBOR_DATA) {
#else
        if (cu_origin_x != 0) {
#endif
            k = 0;
            do {
                ctx_base_left += (dcSignLevelCoeffNeighborArray->leftArray[k + dcSignLevelCoeffLeftNeighborIndex] != 0);
            } while (++k < txb_h_unit);
        }
        const int32_t ctx_base = ((ctx_base_left != 0) + (ctx_base_top != 0));
        const int32_t ctx_offset = (num_pels_log2_lookup[plane_bsize] >
            num_pels_log2_lookup[txsize_to_bsize[tx_size]])
            ? 10
            : 7;
        *txb_skip_ctx = (int16_t)(ctx_base + ctx_offset);
    }
#undef MAX_TX_SIZE_UNIT
}

void Av1WriteTxType(
    PictureParentControlSet_t   *pcsPtr,
    FRAME_CONTEXT               *frameContext,
    aom_writer                  *ecWriter,
    CodingUnit_t                *cu_ptr,
    uint32_t                      intraDir,
    TxType                     txType,
    TxSize                      txSize) {

    const int32_t isInter = (cu_ptr->prediction_mode_flag == INTER_MODE);

    const TxSize squareTxSize = txsize_sqr_map[txSize];

    if (get_ext_tx_types(txSize, isInter, pcsPtr->reduced_tx_set_used) > 1 &&
        (pcsPtr->base_qindex > 0)) {
        const TxSetType txSetType =
            get_ext_tx_set_type(txSize, isInter, pcsPtr->reduced_tx_set_used);
        const int32_t eset = get_ext_tx_set(txSize, isInter, pcsPtr->reduced_tx_set_used);
        // eset == 0 should correspond to a set with only DCT_DCT and there
        // is no need to send the tx_type
        assert(eset > 0);
        assert(av1_ext_tx_used[txSetType][txType]);
        if (isInter) {
            aom_write_symbol(ecWriter, av1_ext_tx_ind[txSetType][txType],
                frameContext->inter_ext_tx_cdf[eset][squareTxSize],
                av1_num_ext_tx_set[txSetType]);
        }
        else {


            aom_write_symbol(
                ecWriter, av1_ext_tx_ind[txSetType][txType],
                frameContext->intra_ext_tx_cdf[eset][squareTxSize][intraDir],
                av1_num_ext_tx_set[txSetType]);
        }
    }
}


static INLINE void set_dc_sign(int32_t *cul_level, int32_t dc_val) {
    if (dc_val < 0)
        *cul_level |= 1 << COEFF_CONTEXT_BITS;
    else if (dc_val > 0)
        *cul_level += 2 << COEFF_CONTEXT_BITS;
}



int32_t  Av1WriteCoeffsTxb1D(
    PictureParentControlSet_t   *parentPcsPtr,
    FRAME_CONTEXT               *frameContext,
    aom_writer                  *ecWriter,
    CodingUnit_t                *cu_ptr,
    TxSize                     txSize,
    uint32_t                       pu_index,
    uint32_t                       tu_index,
    uint32_t                       intraLumaDir,
    int32_t                     *coeffBufferPtr,
    const uint16_t                coeffStride,
    COMPONENT_TYPE              componentType,
    int16_t                      txbSkipCtx,
    int16_t                      dcSignCtx,
    int16_t                      eob)
{
    (void)pu_index;
    (void)coeffStride;
    const TxSize txs_ctx = (TxSize)((txsize_sqr_map[txSize] + txsize_sqr_up_map[txSize] + 1) >> 1);
    TxType txType = cu_ptr->transform_unit_array[tu_index].transform_type[componentType];
    const SCAN_ORDER *const scan_order = &av1_scan_orders[txSize][txType];
    const int16_t *const scan = scan_order->scan;
    int32_t c;
    const int16_t bwl = (const uint16_t)get_txb_bwl(txSize);
    const uint16_t width = (const uint16_t)get_txb_wide(txSize);
    const uint16_t height = (const uint16_t)get_txb_high(txSize);

    uint8_t levelsBuf[TX_PAD_2D];
    uint8_t *const levels = SetLevels(levelsBuf, width);
    DECLARE_ALIGNED(16, int8_t, coeffContexts[MAX_TX_SQUARE]);

    aom_write_symbol(ecWriter, eob == 0,

        frameContext->txb_skip_cdf[txs_ctx][txbSkipCtx], 2);

    if (componentType == 0 && eob == 0) {
        // INTER. Chroma follows Luma in transform type
        if (cu_ptr->prediction_mode_flag == INTER_MODE) {
            txType = cu_ptr->transform_unit_array[tu_index].transform_type[PLANE_TYPE_Y] = DCT_DCT;
            cu_ptr->transform_unit_array[tu_index].transform_type[PLANE_TYPE_UV] = DCT_DCT;
        }
        else { // INTRA
            txType = cu_ptr->transform_unit_array[tu_index].transform_type[PLANE_TYPE_Y] = DCT_DCT;
        }
        assert(txType == DCT_DCT);
    }
    if (eob == 0) return 0;




    av1TxbInitLevels(
        coeffBufferPtr,
        width,
        width,
        height,
        levels);

    if (componentType == COMPONENT_LUMA) {
        Av1WriteTxType(
            parentPcsPtr,
            frameContext,
            ecWriter,
            cu_ptr,
            intraLumaDir,
            txType,
            txSize);
    }

    int16_t eobExtra;
    const int16_t eobPt = GetEobPosToken(eob, &eobExtra);
    const int16_t eobMultiSize = txsize_log2_minus4[txSize];
    const int16_t eobMultiCtx = (tx_type_to_class[txType] == TX_CLASS_2D) ? 0 : 1;
    switch (eobMultiSize) {
    case 0:
        aom_write_symbol(ecWriter, eobPt - 1,
            frameContext->eob_flag_cdf16[componentType][eobMultiCtx], 5);
        break;
    case 1:
        aom_write_symbol(ecWriter, eobPt - 1,
            frameContext->eob_flag_cdf32[componentType][eobMultiCtx], 6);
        break;
    case 2:
        aom_write_symbol(ecWriter, eobPt - 1,
            frameContext->eob_flag_cdf64[componentType][eobMultiCtx], 7);
        break;
    case 3:
        aom_write_symbol(ecWriter, eobPt - 1,
            frameContext->eob_flag_cdf128[componentType][eobMultiCtx], 8);
        break;
    case 4:
        aom_write_symbol(ecWriter, eobPt - 1,
            frameContext->eob_flag_cdf256[componentType][eobMultiCtx], 9);
        break;
    case 5:
        aom_write_symbol(ecWriter, eobPt - 1,
            frameContext->eob_flag_cdf512[componentType][eobMultiCtx], 10);
        break;
    default:
        aom_write_symbol(ecWriter, eobPt - 1,
            frameContext->eob_flag_cdf1024[componentType][eobMultiCtx], 11);
        break;
    }

    uint8_t eobOffsetBits = KEobOffsetBits[eobPt];
    if (eobOffsetBits > 0) {
        int32_t eobShift = eobOffsetBits - 1;
        int32_t bit = (eobExtra & (1 << eobShift)) ? 1 : 0;
        aom_write_symbol(ecWriter, bit, frameContext->eob_extra_cdf[txs_ctx][componentType][eobPt], 2);
        for (int32_t i = 1; i < eobOffsetBits; i++) {
            eobShift = eobOffsetBits - 1 - i;
            bit = (eobExtra & (1 << eobShift)) ? 1 : 0;
            aom_write_bit(ecWriter, bit);
        }
    }

    av1_get_nz_map_contexts(levels, scan, eob, txSize, tx_type_to_class[txType], coeffContexts);

    for (c = eob - 1; c >= 0; --c) {

        const int16_t pos = scan[c];
        const int32_t v = coeffBufferPtr[pos];
        const int16_t coeffCtx = coeffContexts[pos];
        int32_t level = ABS(v);

        if (c == eob - 1) {
            aom_write_symbol(
                ecWriter, AOMMIN(level, 3) - 1,
                frameContext->coeff_base_eob_cdf[txs_ctx][componentType][coeffCtx], 3);
        }
        else {
            aom_write_symbol(ecWriter, AOMMIN(level, 3),
                frameContext->coeff_base_cdf[txs_ctx][componentType][coeffCtx],
                4);
        }
        if (level > NUM_BASE_LEVELS) {
            // level is above 1.
            int32_t base_range = level - 1 - NUM_BASE_LEVELS;
            int16_t brCtx = GetBrCtx(levels, pos, bwl, txType);

            for (int32_t idx = 0; idx < COEFF_BASE_RANGE; idx += BR_CDF_SIZE - 1) {
                const int32_t k = AOMMIN(base_range - idx, BR_CDF_SIZE - 1);
                aom_write_symbol(
                    ecWriter, k,
                    frameContext->coeff_br_cdf[AOMMIN(txs_ctx, TX_32X32)][componentType][brCtx],
                    BR_CDF_SIZE);
                if (k < BR_CDF_SIZE - 1) break;
            }
        }
    }
    // Loop to code all signs in the transform block,
    // starting with the sign of DC (if applicable)


    int32_t cul_level = 0;
    for (c = 0; c < eob; ++c) {

        const int16_t pos = scan[c];
        const int32_t v = coeffBufferPtr[pos];
        int32_t level = ABS(v);
        cul_level += level;

        const int32_t sign = (v < 0) ? 1 : 0;
        if (level) {
            if (c == 0) {
                aom_write_symbol(
                    ecWriter, sign, frameContext->dc_sign_cdf[componentType][dcSignCtx], 2);
            }
            else {
                aom_write_bit(ecWriter, sign);
            }
            if (level > COEFF_BASE_RANGE + NUM_BASE_LEVELS) {
                WriteGolomb(ecWriter,
                    level - COEFF_BASE_RANGE - 1 - NUM_BASE_LEVELS);
            }
        }
    }

    cul_level = AOMMIN(COEFF_CONTEXT_MASK, cul_level);
    // DC value
    set_dc_sign(&cul_level, coeffBufferPtr[0]);
    return cul_level;



}

/************************************
******* Av1EncodeTuCoeff
**************************************/
static EbErrorType Av1EncodeCoeff1D(
    PictureControlSet_t     *pcsPtr,
    EntropyCodingContext_t  *context_ptr,
    FRAME_CONTEXT           *frameContext,
    aom_writer              *ecWriter,
    CodingUnit_t           *cu_ptr,
    uint32_t                  cu_origin_x,
    uint32_t                  cu_origin_y,
    uint32_t                  intraLumaDir,
    BlockSize              plane_bsize,
    EbPictureBufferDesc_t  *coeffPtr,
    NeighborArrayUnit_t     *luma_dc_sign_level_coeff_neighbor_array,
    NeighborArrayUnit_t     *cr_dc_sign_level_coeff_neighbor_array,
    NeighborArrayUnit_t     *cb_dc_sign_level_coeff_neighbor_array)
{

    EbErrorType return_error = EB_ErrorNone;

    const BlockGeom *blk_geom = Get_blk_geom_mds(cu_ptr->mds_idx);
    int32_t cul_level_y, cul_level_cb = 0 , cul_level_cr = 0;




    uint16_t txb_count = blk_geom->txb_count;
    uint8_t txb_itr = 0;

    for (txb_itr = 0; txb_itr < txb_count; txb_itr++) {

        const TxSize tx_size = blk_geom->txsize[txb_itr];
        const TxSize chroma_tx_size = blk_geom->txsize_uv[txb_itr];
        int32_t *coeffBuffer;

        const uint32_t coeff1dOffset = context_ptr->coded_area_sb;

        coeffBuffer = (int32_t*)coeffPtr->bufferY + coeff1dOffset;

        {
            int16_t txbSkipCtx = 0;
            int16_t dcSignCtx = 0;

            GetTxbCtx(
                COMPONENT_LUMA,
                luma_dc_sign_level_coeff_neighbor_array,
                cu_origin_x + blk_geom->tx_org_x[txb_itr] - blk_geom->origin_x,
                cu_origin_y + blk_geom->tx_org_y[txb_itr] - blk_geom->origin_y,
                plane_bsize,
                tx_size,
                &txbSkipCtx,
                &dcSignCtx);


            cul_level_y =
                Av1WriteCoeffsTxb1D(
                    pcsPtr->parent_pcs_ptr,
                    frameContext,
                    ecWriter,
                    cu_ptr,
                    tx_size,
                    0,
                    txb_itr,
                    intraLumaDir,
                    coeffBuffer,
                    coeffPtr->strideY,
                    COMPONENT_LUMA,
                    txbSkipCtx,
                    dcSignCtx,
                    cu_ptr->transform_unit_array[txb_itr].nz_coef_count[0]);
        }

        if (blk_geom->has_uv) {

            // cb
            coeffBuffer = (int32_t*)coeffPtr->bufferCb + context_ptr->coded_area_sb_uv;
            {
                int16_t txbSkipCtx = 0;
                int16_t dcSignCtx = 0;

                GetTxbCtx(
                    COMPONENT_CHROMA,
                    cb_dc_sign_level_coeff_neighbor_array,
                    ROUND_UV(cu_origin_x + blk_geom->tx_org_x[txb_itr] - blk_geom->origin_x) >> 1,
                    ROUND_UV(cu_origin_y + blk_geom->tx_org_y[txb_itr] - blk_geom->origin_y) >> 1,
                    blk_geom->bsize_uv,
                    chroma_tx_size,
                    &txbSkipCtx,
                    &dcSignCtx);


                cul_level_cb =
                    Av1WriteCoeffsTxb1D(
                        pcsPtr->parent_pcs_ptr,
                        frameContext,
                        ecWriter,
                        cu_ptr,
                        chroma_tx_size,
                        0,
                        txb_itr,
                        intraLumaDir,
                        coeffBuffer,
                        coeffPtr->strideCb,
                        COMPONENT_CHROMA,
                        txbSkipCtx,
                        dcSignCtx,
                        cu_ptr->transform_unit_array[txb_itr].nz_coef_count[1]);

            }

            // cr
            coeffBuffer = (int32_t*)coeffPtr->bufferCr + context_ptr->coded_area_sb_uv;
            {

                int16_t txbSkipCtx = 0;
                int16_t dcSignCtx = 0;

                GetTxbCtx(
                    COMPONENT_CHROMA,
                    cr_dc_sign_level_coeff_neighbor_array,
                    ROUND_UV(cu_origin_x + blk_geom->tx_org_x[txb_itr] - blk_geom->origin_x) >> 1,
                    ROUND_UV(cu_origin_y + blk_geom->tx_org_y[txb_itr] - blk_geom->origin_y) >> 1,
                    blk_geom->bsize_uv,
                    chroma_tx_size,
                    &txbSkipCtx,
                    &dcSignCtx);


                cul_level_cr =
                    Av1WriteCoeffsTxb1D(
                        pcsPtr->parent_pcs_ptr,
                        frameContext,
                        ecWriter,
                        cu_ptr,
                        chroma_tx_size,
                        0,
                        txb_itr,
                        intraLumaDir,
                        coeffBuffer,
                        coeffPtr->strideCr,
                        COMPONENT_CHROMA,
                        txbSkipCtx,
                        dcSignCtx,
                        cu_ptr->transform_unit_array[txb_itr].nz_coef_count[2]);
            }
        }

        // Update the luma Dc Sign Level Coeff Neighbor Array
        {
            uint8_t dcSignLevelCoeff = (uint8_t)cul_level_y;
            //if (!txb_ptr->lumaCbf)
            //    dcSignLevelCoeff = 0;
            NeighborArrayUnitModeWrite(
                luma_dc_sign_level_coeff_neighbor_array,
                (uint8_t*)&dcSignLevelCoeff,
                cu_origin_x + blk_geom->tx_org_x[txb_itr] - blk_geom->origin_x,
                cu_origin_y + blk_geom->tx_org_y[txb_itr] - blk_geom->origin_y,
                blk_geom->tx_width[txb_itr],
                blk_geom->tx_height[txb_itr],
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
        }

        // Update the cb Dc Sign Level Coeff Neighbor Array

        if (blk_geom->has_uv)
        {
            uint8_t dcSignLevelCoeff = (uint8_t)cul_level_cb;
            NeighborArrayUnitModeWrite(
                cb_dc_sign_level_coeff_neighbor_array,
                (uint8_t*)&dcSignLevelCoeff,
                ROUND_UV(cu_origin_x + blk_geom->tx_org_x[txb_itr] - blk_geom->origin_x) >> 1,
                ROUND_UV(cu_origin_y + blk_geom->tx_org_y[txb_itr] - blk_geom->origin_y) >> 1,
                blk_geom->tx_width_uv[txb_itr],
                blk_geom->tx_height_uv[txb_itr],
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

        }

        if (blk_geom->has_uv)
            // Update the cr DC Sign Level Coeff Neighbor Array
        {
            uint8_t dcSignLevelCoeff = (uint8_t)cul_level_cr;
            NeighborArrayUnitModeWrite(
                cr_dc_sign_level_coeff_neighbor_array,
                (uint8_t*)&dcSignLevelCoeff,
                ROUND_UV(cu_origin_x + blk_geom->tx_org_x[txb_itr] - blk_geom->origin_x) >> 1,
                ROUND_UV(cu_origin_y + blk_geom->tx_org_y[txb_itr] - blk_geom->origin_y) >> 1,
                blk_geom->tx_width_uv[txb_itr],
                blk_geom->tx_height_uv[txb_itr],
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
        }


        context_ptr->coded_area_sb += blk_geom->tx_width[txb_itr] * blk_geom->tx_height[txb_itr];

        if (blk_geom->has_uv)
            context_ptr->coded_area_sb_uv += blk_geom->tx_width_uv[txb_itr] * blk_geom->tx_height_uv[txb_itr];

    }

    return return_error;
}

/*********************************************************************
* EncodePartitionAv1
*   Encodes the partition
*********************************************************************/
// Return the number of elements in the partition CDF when
// partitioning the (square) block with luma block size of bsize.
static INLINE int32_t partition_cdf_length(BlockSize bsize) {
    if (bsize <= BLOCK_8X8)
        return PARTITION_TYPES;
    else if (bsize == BLOCK_128X128)
        return EXT_PARTITION_TYPES - 2;
    else
        return EXT_PARTITION_TYPES;
}
static int32_t cdf_element_prob(const aom_cdf_prob *const cdf,
    size_t element) {
    assert(cdf != NULL);
    return (element > 0 ? cdf[element - 1] : CDF_PROB_TOP) - cdf[element];
}
static void partition_gather_horz_alike(aom_cdf_prob *out,
    const aom_cdf_prob *const in,
    BlockSize bsize) {


    out[0] = CDF_PROB_TOP;
    out[0] -= cdf_element_prob(in, PARTITION_HORZ);
    out[0] -= cdf_element_prob(in, PARTITION_SPLIT);
    out[0] -= cdf_element_prob(in, PARTITION_HORZ_A);
    out[0] -= cdf_element_prob(in, PARTITION_HORZ_B);
    out[0] -= cdf_element_prob(in, PARTITION_VERT_A);
    if (bsize != BLOCK_128X128) out[0] -= cdf_element_prob(in, PARTITION_HORZ_4);
    out[0] = AOM_ICDF(out[0]);
    out[1] = AOM_ICDF(CDF_PROB_TOP);
}
static void partition_gather_vert_alike(aom_cdf_prob *out,
    const aom_cdf_prob *const in,
    BlockSize bsize) {
    out[0] = CDF_PROB_TOP;
    out[0] -= cdf_element_prob(in, PARTITION_VERT);
    out[0] -= cdf_element_prob(in, PARTITION_SPLIT);
    out[0] -= cdf_element_prob(in, PARTITION_HORZ_A);
    out[0] -= cdf_element_prob(in, PARTITION_VERT_A);
    out[0] -= cdf_element_prob(in, PARTITION_VERT_B);
    if (bsize != BLOCK_128X128) out[0] -= cdf_element_prob(in, PARTITION_VERT_4);
    out[0] = AOM_ICDF(out[0]);
    out[1] = AOM_ICDF(CDF_PROB_TOP);
}
static void EncodePartitionAv1(
    SequenceControlSet_t    *sequence_control_set_ptr,
    FRAME_CONTEXT           *frameContext,
    aom_writer              *ecWriter,
    BlockSize              bsize,
    PartitionType          p,
    uint32_t                  cu_origin_x,
    uint32_t                  cu_origin_y,
    NeighborArrayUnit_t    *partition_context_neighbor_array)
{

    const int32_t is_partition_point = bsize >= BLOCK_8X8;

    if (!is_partition_point) return;

    const int32_t hbs = (mi_size_wide[bsize] << 2) >> 1;
    const int32_t hasRows = (cu_origin_y + hbs) < sequence_control_set_ptr->luma_height;
    const int32_t hasCols = (cu_origin_x + hbs) < sequence_control_set_ptr->luma_width;

    uint32_t partitionContextLeftNeighborIndex = GetNeighborArrayUnitLeftIndex(
        partition_context_neighbor_array,
        cu_origin_y);
    uint32_t partitionContextTopNeighborIndex = GetNeighborArrayUnitTopIndex(
        partition_context_neighbor_array,
        cu_origin_x);

    uint32_t contextIndex = 0;

    const PARTITION_CONTEXT above_ctx = (((PartitionContext*)partition_context_neighbor_array->topArray)[partitionContextTopNeighborIndex].above == (int8_t)INVALID_NEIGHBOR_DATA) ?
        0 : ((PartitionContext*)partition_context_neighbor_array->topArray)[partitionContextTopNeighborIndex].above;
    const PARTITION_CONTEXT left_ctx = (((PartitionContext*)partition_context_neighbor_array->leftArray)[partitionContextLeftNeighborIndex].left == (int8_t)INVALID_NEIGHBOR_DATA) ?
        0 : ((PartitionContext*)partition_context_neighbor_array->leftArray)[partitionContextLeftNeighborIndex].left;

    const int32_t bsl = mi_size_wide_log2[bsize] - mi_size_wide_log2[BLOCK_8X8];
    int32_t above = (above_ctx >> bsl) & 1, left = (left_ctx >> bsl) & 1;

    assert(mi_size_wide_log2[bsize] == mi_size_high_log2[bsize]);
    assert(bsl >= 0);

    contextIndex = (left * 2 + above) + bsl * PARTITION_PLOFFSET;

    if (!hasRows && !hasCols) {
        assert(p == PARTITION_SPLIT);
        return;
    }

    if (hasRows && hasCols) {

        aom_write_symbol(
            ecWriter,
            p,
            frameContext->partition_cdf[contextIndex],
            partition_cdf_length(bsize));

    }
    else if (!hasRows && hasCols) {
        aom_cdf_prob cdf[2];
        partition_gather_vert_alike(cdf, frameContext->partition_cdf[contextIndex], bsize);
        aom_write_symbol(
            ecWriter,
            p == PARTITION_SPLIT,
            cdf,
            2);
    }
    else {
        aom_cdf_prob cdf[2];
        partition_gather_horz_alike(cdf, frameContext->partition_cdf[contextIndex], bsize);
        aom_write_symbol(
            ecWriter,
            p == PARTITION_SPLIT,
            cdf,
            2);
    }

    return;
}

/*********************************************************************
* EncodeSkipCoeffAv1
*   Encodes the skip coefficient flag
*********************************************************************/
static void EncodeSkipCoeffAv1(
    FRAME_CONTEXT           *frameContext,
    aom_writer              *ecWriter,
    EbBool                 skipCoeffFlag,
    uint32_t                  cu_origin_x,
    uint32_t                  cu_origin_y,
    NeighborArrayUnit_t    *skip_coeff_neighbor_array)
{

    uint32_t skipCoeffLeftNeighborIndex = GetNeighborArrayUnitLeftIndex(
        skip_coeff_neighbor_array,
        cu_origin_y);
    uint32_t skipCoeffTopNeighborIndex = GetNeighborArrayUnitTopIndex(
        skip_coeff_neighbor_array,
        cu_origin_x);

    uint32_t contextIndex = 0;

    contextIndex =
        (skip_coeff_neighbor_array->leftArray[skipCoeffLeftNeighborIndex] == (uint8_t)INVALID_NEIGHBOR_DATA) ? 0 :
        (skip_coeff_neighbor_array->leftArray[skipCoeffLeftNeighborIndex]) ? 1 : 0;


    contextIndex +=
        (skip_coeff_neighbor_array->topArray[skipCoeffTopNeighborIndex] == (uint8_t)INVALID_NEIGHBOR_DATA) ? 0 :
        (skip_coeff_neighbor_array->topArray[skipCoeffTopNeighborIndex]) ? 1 : 0;

    aom_write_symbol(
        ecWriter,
        skipCoeffFlag ? 1 : 0,
        frameContext->skip_cdfs[contextIndex],
        2);

    return;
}
/*********************************************************************
* EncodeIntraLumaModeAv1
*   Encodes the Intra Luma Mode
*********************************************************************/
static void EncodeIntraLumaModeAv1(
    FRAME_CONTEXT           *frameContext,
    aom_writer              *ecWriter,
    CodingUnit_t            *cu_ptr,
    uint32_t                  cu_origin_x,
    uint32_t                  cu_origin_y,
    uint32_t                  lumaMode,
    NeighborArrayUnit_t    *mode_type_neighbor_array,
    NeighborArrayUnit_t    *intra_luma_mode_neighbor_array)
{
    uint32_t modeTypeLeftNeighborIndex = GetNeighborArrayUnitLeftIndex(
        mode_type_neighbor_array,
        cu_origin_y);
    uint32_t modeTypeTopNeighborIndex = GetNeighborArrayUnitTopIndex(
        mode_type_neighbor_array,
        cu_origin_x);
    uint32_t intraLumaModeLeftNeighborIndex = GetNeighborArrayUnitLeftIndex(
        intra_luma_mode_neighbor_array,
        cu_origin_y);
    uint32_t intraLumaModeTopNeighborIndex = GetNeighborArrayUnitTopIndex(
        intra_luma_mode_neighbor_array,
        cu_origin_x);


    uint32_t topContext = 0, leftContext = 0;

    uint32_t left_neighbor_mode = (uint32_t)(
        (mode_type_neighbor_array->leftArray[modeTypeLeftNeighborIndex] != INTRA_MODE) ? (uint32_t)DC_PRED :
        intra_luma_mode_neighbor_array->leftArray[intraLumaModeLeftNeighborIndex]);

    uint32_t top_neighbor_mode = (uint32_t)(
        (mode_type_neighbor_array->topArray[modeTypeTopNeighborIndex] != INTRA_MODE) ? (uint32_t)DC_PRED :
        intra_luma_mode_neighbor_array->topArray[intraLumaModeTopNeighborIndex]);

    topContext = intra_mode_context[top_neighbor_mode];
    leftContext = intra_mode_context[left_neighbor_mode];

    aom_write_symbol(
        ecWriter,
        lumaMode,
        frameContext->kf_y_cdf[topContext][leftContext],
        INTRA_MODES);

    if (cu_ptr->pred_mode != INTRA_MODE_4x4)

        if (cu_ptr->prediction_unit_array[0].use_angle_delta && cu_ptr->prediction_unit_array[0].is_directional_mode_flag) {
            aom_write_symbol(ecWriter,
                cu_ptr->prediction_unit_array[0].angle_delta[PLANE_TYPE_Y] + MAX_ANGLE_DELTA,
                frameContext->angle_delta_cdf[lumaMode - V_PRED],
                2 * MAX_ANGLE_DELTA + 1);

        }

    return;
}
/*********************************************************************
* EncodeIntraLumaModeNonKeyAv1
*   Encodes the Intra Luma Mode for non Key frames
*********************************************************************/
static void EncodeIntraLumaModeNonKeyAv1(
    FRAME_CONTEXT           *frameContext,
    aom_writer              *ecWriter,
    CodingUnit_t            *cu_ptr,

    BlockSize                bsize,
    uint32_t                  lumaMode)
{

    aom_write_symbol(
        ecWriter,
        lumaMode,
        frameContext->y_mode_cdf[size_group_lookup[bsize]],
        INTRA_MODES);

    if (cu_ptr->pred_mode != INTRA_MODE_4x4)

        if (cu_ptr->prediction_unit_array[0].use_angle_delta && cu_ptr->prediction_unit_array[0].is_directional_mode_flag) {
            aom_write_symbol(ecWriter,
                cu_ptr->prediction_unit_array[0].angle_delta[PLANE_TYPE_Y] + MAX_ANGLE_DELTA,
                frameContext->angle_delta_cdf[lumaMode - V_PRED],
                2 * MAX_ANGLE_DELTA + 1);
        }

    return;
}

static void write_cfl_alphas(FRAME_CONTEXT *const ec_ctx, int32_t idx,
    int32_t joint_sign, aom_writer *w) {
    aom_write_symbol(w, joint_sign, ec_ctx->cfl_sign_cdf, CFL_JOINT_SIGNS);
    // Magnitudes are only signaled for nonzero codes.
    if (CFL_SIGN_U(joint_sign) != CFL_SIGN_ZERO) {
        aom_cdf_prob *cdf_u = ec_ctx->cfl_alpha_cdf[CFL_CONTEXT_U(joint_sign)];
        aom_write_symbol(w, CFL_IDX_U(idx), cdf_u, CFL_ALPHABET_SIZE);
    }
    if (CFL_SIGN_V(joint_sign) != CFL_SIGN_ZERO) {
        aom_cdf_prob *cdf_v = ec_ctx->cfl_alpha_cdf[CFL_CONTEXT_V(joint_sign)];
        aom_write_symbol(w, CFL_IDX_V(idx), cdf_v, CFL_ALPHABET_SIZE);
    }
}


/*********************************************************************
* EncodeIntraChromaModeAv1
*   Encodes the Intra Chroma Mode
*********************************************************************/
static void EncodeIntraChromaModeAv1(
    FRAME_CONTEXT           *frameContext,
    aom_writer              *ecWriter,
    CodingUnit_t            *cu_ptr,
    uint32_t                  lumaMode,
    uint32_t                  chroma_mode,
    uint8_t                   cflAllowed)
{
    aom_write_symbol(
        ecWriter,
        chroma_mode,
        frameContext->uv_mode_cdf[cflAllowed][lumaMode],
        UV_INTRA_MODES - !cflAllowed);

    if (chroma_mode == UV_CFL_PRED)
        write_cfl_alphas(frameContext, cu_ptr->prediction_unit_array->cfl_alpha_idx, cu_ptr->prediction_unit_array->cfl_alpha_signs, ecWriter);

    if (cu_ptr->pred_mode != INTRA_MODE_4x4)
        if (cu_ptr->prediction_unit_array[0].use_angle_delta && cu_ptr->prediction_unit_array[0].is_directional_chroma_mode_flag) {
            aom_write_symbol(ecWriter,
                cu_ptr->prediction_unit_array[0].angle_delta[PLANE_TYPE_UV] + MAX_ANGLE_DELTA,
                frameContext->angle_delta_cdf[chroma_mode - V_PRED],
                2 * MAX_ANGLE_DELTA + 1);

        }

    return;
}

/*********************************************************************
* EncodeSkipModeAv1
*   Encodes the skip Mode flag
*********************************************************************/
static void EncodeSkipModeAv1(
    FRAME_CONTEXT           *frameContext,
    aom_writer              *ecWriter,
    EbBool                 skipModeFlag,
    uint32_t                  cu_origin_x,
    uint32_t                  cu_origin_y,
    NeighborArrayUnit_t    *skip_flag_neighbor_array)
{

    uint32_t skipFlagLeftNeighborIndex = GetNeighborArrayUnitLeftIndex(
        skip_flag_neighbor_array,
        cu_origin_y);
    uint32_t skipFlagTopNeighborIndex = GetNeighborArrayUnitTopIndex(
        skip_flag_neighbor_array,
        cu_origin_x);

    uint32_t contextIndex = 0;

    contextIndex =
        (skip_flag_neighbor_array->leftArray[skipFlagLeftNeighborIndex] == (uint8_t)INVALID_NEIGHBOR_DATA) ? 0 :
        (skip_flag_neighbor_array->leftArray[skipFlagLeftNeighborIndex]) ? 1 : 0;


    contextIndex +=
        (skip_flag_neighbor_array->topArray[skipFlagTopNeighborIndex] == (uint8_t)INVALID_NEIGHBOR_DATA) ? 0 :
        (skip_flag_neighbor_array->topArray[skipFlagTopNeighborIndex]) ? 1 : 0;

    aom_write_symbol(
        ecWriter,
        skipModeFlag ? 1 : 0,
        frameContext->skip_mode_cdfs[contextIndex],
        2);

    return;
}
/*********************************************************************
* EncodePredModeAv1
*   Encodes the Prediction Mode
*********************************************************************/
static void EncodePredModeAv1(
    FRAME_CONTEXT           *frameContext,
    aom_writer              *ecWriter,
    EbBool                 predModeFlag,
    uint32_t                  cu_origin_x,
    uint32_t                  cu_origin_y,
    NeighborArrayUnit_t    *mode_type_neighbor_array)
{
    uint32_t modeTypeLeftNeighborIndex = GetNeighborArrayUnitLeftIndex(
        mode_type_neighbor_array,
        cu_origin_y);
    uint32_t modeTypeTopNeighborIndex = GetNeighborArrayUnitTopIndex(
        mode_type_neighbor_array,
        cu_origin_x);

    uint32_t contextIndex = 0;

    if (mode_type_neighbor_array->leftArray[modeTypeLeftNeighborIndex] != (uint8_t)INVALID_MODE && mode_type_neighbor_array->topArray[modeTypeTopNeighborIndex] != (uint8_t)INVALID_MODE) {
        contextIndex = (mode_type_neighbor_array->leftArray[modeTypeLeftNeighborIndex] == (uint8_t)INTRA_MODE && mode_type_neighbor_array->topArray[modeTypeTopNeighborIndex] == (uint8_t)INTRA_MODE) ? 3 :
            (mode_type_neighbor_array->leftArray[modeTypeLeftNeighborIndex] == (uint8_t)INTRA_MODE || mode_type_neighbor_array->topArray[modeTypeTopNeighborIndex] == (uint8_t)INTRA_MODE) ? 1 : 0;

    }
    else  if (mode_type_neighbor_array->leftArray[modeTypeLeftNeighborIndex] != (uint8_t)INVALID_MODE) {
        contextIndex = (mode_type_neighbor_array->leftArray[modeTypeLeftNeighborIndex] == (uint8_t)INTRA_MODE) ? 2 : 0;
    }
    else if (mode_type_neighbor_array->topArray[modeTypeTopNeighborIndex] != (uint8_t)INVALID_MODE) {
        contextIndex = (mode_type_neighbor_array->topArray[modeTypeTopNeighborIndex] == (uint8_t)INTRA_MODE) ? 2 : 0;
    }
    else {
        contextIndex = 0;
    }

    aom_write_symbol(
        ecWriter,
        predModeFlag == INTER_MODE ? 1 : 0,
        frameContext->intra_inter_cdf[contextIndex],
        2);

    return;
}

//****************************************************************************************************//

/*********************************************************************
* motion_mode_allowed
*   checks the motion modes that are allowed for the current block
*********************************************************************/

static INLINE int is_inter_mode(PredictionMode mode)
{
    return mode >= SINGLE_INTER_MODE_START && mode < SINGLE_INTER_MODE_END;
}

static INLINE int is_global_mv_block(
    const PredictionMode          mode,
    const BlockSize               bsize,
    TransformationType            type)
{
    return (mode == GLOBALMV || mode == GLOBAL_GLOBALMV)
            && type > TRANSLATION
            && is_motion_variation_allowed_bsize(bsize);
}

MOTION_MODE motion_mode_allowed(
    const PictureControlSet_t       *picture_control_set_ptr,
    const CodingUnit_t              *cu_ptr,
    const BlockSize                 bsize,
    MvReferenceFrame                rf0,
    MvReferenceFrame                rf1,
    PredictionMode                  mode)
{
    if(!picture_control_set_ptr->parent_pcs_ptr->switchable_motion_mode)
        return SIMPLE_TRANSLATION;

    if (picture_control_set_ptr->parent_pcs_ptr->cur_frame_force_integer_mv == 0) {
        const TransformationType gm_type =
            picture_control_set_ptr->parent_pcs_ptr->global_motion[rf0].wmtype;
        if (is_global_mv_block(mode, bsize, gm_type))
            return SIMPLE_TRANSLATION;
    }

    if (is_motion_variation_allowed_bsize(bsize) &&
        is_inter_mode(mode) &&
        rf1 != INTRA_FRAME &&
        !(rf1 > INTRA_FRAME) ) // is_motion_variation_allowed_compound
    {
        if (!has_overlappable_candidates(cu_ptr)) // check_num_overlappable_neighbors
            return SIMPLE_TRANSLATION;

        if (cu_ptr->prediction_unit_array[0].num_proj_ref >= 1 &&
           (picture_control_set_ptr->parent_pcs_ptr->allow_warped_motion)) // TODO(JS): when scale is added, put: && !av1_is_scaled(&(xd->block_refs[0]->sf))
        {
            if (picture_control_set_ptr->parent_pcs_ptr->cur_frame_force_integer_mv)
                return OBMC_CAUSAL;
            return WARPED_CAUSAL;
        }
        return OBMC_CAUSAL;
    } else
        return SIMPLE_TRANSLATION;
}

/*********************************************************************
* write_motion_mode
*   Encodes the Motion Mode (obmc or warped)
*********************************************************************/
static void write_motion_mode(
    FRAME_CONTEXT            *frame_context,
    aom_writer               *ec_writer,
    BlockSize                 bsize,
    MOTION_MODE               motion_mode,
    MvReferenceFrame          rf0,
    MvReferenceFrame          rf1,
    CodingUnit_t             *cu_ptr,
    PictureControlSet_t      *picture_control_set_ptr)
{
    const PredictionMode mode = cu_ptr->prediction_unit_array[0].inter_mode;
    MOTION_MODE last_motion_mode_allowed =
        motion_mode_allowed(picture_control_set_ptr, cu_ptr, bsize, rf0, rf1, mode);

    switch (last_motion_mode_allowed) {
    case SIMPLE_TRANSLATION: break;
    case OBMC_CAUSAL:
        aom_write_symbol(
            ec_writer,
            SIMPLE_TRANSLATION, // motion_mode == OBMC_CAUSAL, TODO: support OBMC
            frame_context->obmc_cdf[bsize],
            2);
        break;
    default:
        aom_write_symbol(
            ec_writer,
            motion_mode,
            frame_context->motion_mode_cdf[bsize],
            MOTION_MODES);
    }

    return;
}

//****************************************************************************************************//
extern  int8_t av1_ref_frame_type(const MvReferenceFrame *const rf);
uint16_t compound_mode_ctx_map[3][COMP_NEWMV_CTXS] = {
   { 0, 1, 1, 1, 1 },
   { 1, 2, 3, 4, 4 },
   { 4, 4, 5, 6, 7 },
};

static int16_t Av1ModeContextAnalyzer(
    const int16_t *const mode_context, const MvReferenceFrame *const rf) {
    const int8_t ref_frame = av1_ref_frame_type(rf);

    if (rf[1] <= INTRA_FRAME) return mode_context[ref_frame];

    const int16_t newmv_ctx = mode_context[ref_frame] & NEWMV_CTX_MASK;
    const int16_t refmv_ctx =
        (mode_context[ref_frame] >> REFMV_OFFSET) & REFMV_CTX_MASK;
    ASSERT((refmv_ctx >> 1) < 3);
    const int16_t comp_ctx = compound_mode_ctx_map[refmv_ctx >> 1][AOMMIN(
        newmv_ctx, COMP_NEWMV_CTXS - 1)];
    return comp_ctx;
}




EbErrorType EncodeSliceFinish(
    EntropyCoder_t        *entropy_coder_ptr)
{
    EbErrorType return_error = EB_ErrorNone;

    aom_stop_encode(&entropy_coder_ptr->ecWriter);

    return return_error;
}

EbErrorType ResetBitstream(
    EbPtr bitstreamPtr)
{
    EbErrorType return_error = EB_ErrorNone;
    OutputBitstreamUnit_t *outputBitstreamPtr = (OutputBitstreamUnit_t*)bitstreamPtr;

    output_bitstream_reset(
        outputBitstreamPtr);

    return return_error;
}

EbErrorType ResetEntropyCoder(
    EncodeContext_t            *encode_context_ptr,
    EntropyCoder_t             *entropy_coder_ptr,
    uint32_t                      qp,
    EB_SLICE                    slice_type)
{
    EbErrorType return_error = EB_ErrorNone;

    (void)encode_context_ptr;
    (void)slice_type;
    av1_default_coef_probs(entropy_coder_ptr->fc, qp);
    init_mode_probs(entropy_coder_ptr->fc);


    return return_error;
}

EbErrorType CopyRbspBitstreamToPayload(
    Bitstream_t *bitstreamPtr,
    EbByte      outputBuffer,
    uint32_t      *outputBufferIndex,
    uint32_t      *outputBufferSize,
    EncodeContext_t         *encode_context_ptr)
{
    EbErrorType return_error = EB_ErrorNone;
    OutputBitstreamUnit_t *outputBitstreamPtr = (OutputBitstreamUnit_t*)bitstreamPtr->outputBitstreamPtr;


    CHECK_REPORT_ERROR(
        ((outputBitstreamPtr->writtenBitsCount >> 3) + (*outputBufferIndex) < (*outputBufferSize)),
        encode_context_ptr->app_callback_ptr,
        EB_ENC_EC_ERROR2);



    output_bitstream_rbsp_to_payload(
        outputBitstreamPtr,
        outputBuffer,
        outputBufferIndex,
        outputBufferSize,
        0);

    return return_error;
}


EbErrorType BitstreamCtor(
    Bitstream_t **bitstreamDblPtr,
    uint32_t bufferSize)
{
    EbErrorType return_error = EB_ErrorNone;
    EB_MALLOC(Bitstream_t*, *bitstreamDblPtr, sizeof(Bitstream_t), EB_N_PTR);

    EB_MALLOC(EbPtr, (*bitstreamDblPtr)->outputBitstreamPtr, sizeof(OutputBitstreamUnit_t), EB_N_PTR);

    return_error = output_bitstream_unit_ctor(
        (OutputBitstreamUnit_t *)(*bitstreamDblPtr)->outputBitstreamPtr,
        bufferSize);

    return return_error;
}




EbErrorType EntropyCoderCtor(
    EntropyCoder_t **entropyCoderDblPtr,
    uint32_t bufferSize)
{
    EbErrorType return_error = EB_ErrorNone;
    EB_MALLOC(EntropyCoder_t*, *entropyCoderDblPtr, sizeof(EntropyCoder_t), EB_N_PTR);

    EB_MALLOC(EbPtr, (*entropyCoderDblPtr)->cabacEncodeContextPtr, sizeof(CabacEncodeContext_t), EB_N_PTR);

    EB_MALLOC(FRAME_CONTEXT*, (*entropyCoderDblPtr)->fc, sizeof(FRAME_CONTEXT), EB_N_PTR);

    EB_MALLOC(EbPtr, (*entropyCoderDblPtr)->ecOutputBitstreamPtr, sizeof(OutputBitstreamUnit_t), EB_N_PTR);

    return_error = output_bitstream_unit_ctor(
        (OutputBitstreamUnit_t *)(*entropyCoderDblPtr)->ecOutputBitstreamPtr,
        bufferSize);

    CabacCtor(
        (CabacEncodeContext_t *)(*entropyCoderDblPtr)->cabacEncodeContextPtr);


    return_error = output_bitstream_unit_ctor(
        &((((CabacEncodeContext_t*)(*entropyCoderDblPtr)->cabacEncodeContextPtr)->bacEncContext).m_pcTComBitIf),
        bufferSize);

    return return_error;
}




EbPtr EntropyCoderGetBitstreamPtr(
    EntropyCoder_t *entropy_coder_ptr)
{
    CabacEncodeContext_t *cabacEncCtxPtr = (CabacEncodeContext_t*)entropy_coder_ptr->cabacEncodeContextPtr;
    EbPtr bitstreamPtr = (EbPtr)&(cabacEncCtxPtr->bacEncContext.m_pcTComBitIf);

    return bitstreamPtr;
}

//*******************************************************************************************//
//*******************************************************************************************//
//*******************************************************************************************//
//*******************************************************************************************//
// aom_integer.c
static const size_t kMaximumLeb128Size = 8;
static const uint64_t kMaximumLeb128Value = 0xFFFFFFFFFFFFFF;  // 2 ^ 56 - 1

size_t aom_uleb_size_in_bytes(uint64_t value) {
    size_t size = 0;
    do {
        ++size;
    } while ((value >>= 7) != 0);
    return size;
}

int32_t aom_uleb_encode(uint64_t value, size_t available, uint8_t *coded_value,
    size_t *coded_size) {
    const size_t leb_size = aom_uleb_size_in_bytes(value);
    if (value > kMaximumLeb128Value || leb_size > kMaximumLeb128Size ||
        leb_size > available || !coded_value || !coded_size) {
        return -1;
    }

    for (size_t i = 0; i < leb_size; ++i) {
        uint8_t byte = value & 0x7f;
        value >>= 7;

        if (value != 0) byte |= 0x80;  // Signal that more bytes follow.

        *(coded_value + i) = byte;
    }

    *coded_size = leb_size;
    return 0;
}


int32_t aom_wb_is_byte_aligned(const struct aom_write_bit_buffer *wb) {
    return (wb->bit_offset % CHAR_BIT == 0);
}


uint32_t aom_wb_bytes_written(const struct aom_write_bit_buffer *wb) {
    return wb->bit_offset / CHAR_BIT + (wb->bit_offset % CHAR_BIT > 0);
}

void aom_wb_write_bit(struct aom_write_bit_buffer *wb, int32_t bit) {
    const int32_t off = (int32_t)wb->bit_offset;
    const int32_t p = off / CHAR_BIT;
    const int32_t q = CHAR_BIT - 1 - off % CHAR_BIT;
    if (q == CHAR_BIT - 1) {
        // Zero next char and write bit
        wb->bit_buffer[p] = (uint8_t)(bit << q);
    }
    else {
        wb->bit_buffer[p] &= ~(1 << q);
        wb->bit_buffer[p] |= bit << q;
    }
    wb->bit_offset = off + 1;
}

void aom_wb_overwrite_bit(struct aom_write_bit_buffer *wb, int32_t bit) {
    // Do not zero bytes but overwrite exisiting values
    const int32_t off = (int32_t)wb->bit_offset;
    const int32_t p = off / CHAR_BIT;
    const int32_t q = CHAR_BIT - 1 - off % CHAR_BIT;
    wb->bit_buffer[p] &= ~(1 << q);
    wb->bit_buffer[p] |= bit << q;
    wb->bit_offset = off + 1;
}

void aom_wb_write_literal(struct aom_write_bit_buffer *wb, int32_t data, int32_t bits) {
    int32_t bit;
    for (bit = bits - 1; bit >= 0; bit--) aom_wb_write_bit(wb, (data >> bit) & 1);
}

void aom_wb_write_inv_signed_literal(struct aom_write_bit_buffer *wb, int32_t data,
    int32_t bits) {
    aom_wb_write_literal(wb, data, bits + 1);
}


//*******************************************************************************************//

static void WriteInterMode(
    FRAME_CONTEXT       *frameContext,
    aom_writer          *ecWriter,
    PredictionMode     mode,
    const int16_t       mode_ctx,
    uint32_t                  cu_origin_x,
    uint32_t                  cu_origin_y)
{

    (void)cu_origin_x;
    (void)cu_origin_y;
    int16_t newmv_ctx = mode_ctx & NEWMV_CTX_MASK;
    ASSERT(newmv_ctx<NEWMV_MODE_CONTEXTS);
    aom_write_symbol(ecWriter, mode != NEWMV, frameContext->newmv_cdf[newmv_ctx], 2);

    if (mode != NEWMV) {
        const int16_t zeromv_ctx =
            (mode_ctx >> GLOBALMV_OFFSET) & GLOBALMV_CTX_MASK;
        aom_write_symbol(ecWriter, mode != GLOBALMV, frameContext->zeromv_cdf[zeromv_ctx], 2);

        if (mode != GLOBALMV) {
            int16_t refmv_ctx = (mode_ctx >> REFMV_OFFSET) & REFMV_CTX_MASK;
            ASSERT(refmv_ctx<REFMV_MODE_CONTEXTS);            
            aom_write_symbol(ecWriter, mode != NEARESTMV, frameContext->refmv_cdf[refmv_ctx], 2);
        }
    }
}


//extern INLINE int8_t av1_ref_frame_type(const MvReferenceFrame *const rf);
extern uint8_t av1_drl_ctx(const CandidateMv *ref_mv_stack, int32_t ref_idx);

void WriteDrlIdx(
    FRAME_CONTEXT       *frameContext,
    aom_writer          *ecWriter,
    CodingUnit_t        *cu_ptr) {

    //uint8_t ref_frame_type = av1_ref_frame_type(mbmi->ref_frame);
    uint8_t ref_frame_type = cu_ptr->prediction_unit_array[0].ref_frame_type;
    MacroBlockD*  xd = cu_ptr->av1xd;
    //assert(mbmi->ref_mv_idx < 3);

    //xd->ref_mv_stack[ref_frame][ref_mv_idx].comp_mv;

    const int32_t new_mv = cu_ptr->pred_mode == NEWMV || cu_ptr->pred_mode == NEW_NEWMV;
    if (new_mv) {
        int32_t idx;
        for (idx = 0; idx < 2; ++idx) {
            if (xd->ref_mv_count[ref_frame_type] > idx + 1) {
                uint8_t drl_ctx =
                    av1_drl_ctx(xd->final_ref_mv_stack, idx);

                aom_write_symbol(ecWriter, cu_ptr->drl_index != idx, frameContext->drl_cdf[drl_ctx],
                    2);

                if (cu_ptr->drl_index == idx) return;
            }
        }
        return;
    }

    if (have_nearmv_in_inter_mode(cu_ptr->pred_mode)) {
        int32_t idx;
        // TODO(jingning): Temporary solution to compensate the NEARESTMV offset.
        for (idx = 1; idx < 3; ++idx) {
            if (xd->ref_mv_count[ref_frame_type] > idx + 1) {
                uint8_t drl_ctx =
                    av1_drl_ctx(xd->final_ref_mv_stack, idx);
    
                aom_write_symbol(ecWriter, cu_ptr->drl_index != (idx - 1),
                    frameContext->drl_cdf[drl_ctx], 2);

                if (cu_ptr->drl_index == (idx - 1)) return;
            }
        }
        return;
    }
}

extern MV_JOINT_TYPE av1_get_mv_joint(int32_t diff[2]);
static void encode_mv_component(aom_writer *w, int32_t comp, nmv_component *mvcomp,
    MvSubpelPrecision precision) {
    int32_t offset;
    const int32_t sign = comp < 0;
    const int32_t mag = sign ? -comp : comp;
    const int32_t mv_class = av1_get_mv_class(mag - 1, &offset);
    const int32_t d = offset >> 3;         // int32_t mv data
    const int32_t fr = (offset >> 1) & 3;  // fractional mv data
    const int32_t hp = offset & 1;         // high precision mv data

    assert(comp != 0);

    // Sign
    aom_write_symbol(w, sign, mvcomp->sign_cdf, 2);

    // Class
    aom_write_symbol(w, mv_class, mvcomp->classes_cdf, MV_CLASSES);

    // Integer bits
    if (mv_class == MV_CLASS_0) {
        aom_write_symbol(w, d, mvcomp->class0_cdf, CLASS0_SIZE);
    }
    else {
        int32_t i;
        const int32_t n = mv_class + CLASS0_BITS - 1;  // number of bits
        for (i = 0; i < n; ++i)
            aom_write_symbol(w, (d >> i) & 1, mvcomp->bits_cdf[i], 2);
    }
    // Fractional bits
    if (precision > MV_SUBPEL_NONE) {
        aom_write_symbol(
            w, fr,
            mv_class == MV_CLASS_0 ? mvcomp->class0_fp_cdf[d] : mvcomp->fp_cdf,
            MV_FP_SIZE);
    }

    // High precision bit
    if (precision > MV_SUBPEL_LOW_PRECISION)
        aom_write_symbol(
            w, hp, mv_class == MV_CLASS_0 ? mvcomp->class0_hp_cdf : mvcomp->hp_cdf,
            2);
}

MV_JOINT_TYPE av1_get_mv_joint_diff(int32_t diff[2]) {
    if (diff[0] == 0) {
        return diff[1] == 0 ? MV_JOINT_ZERO : MV_JOINT_HNZVZ;
    }
    else {
        return diff[1] == 0 ? MV_JOINT_HZVNZ : MV_JOINT_HNZVNZ;
    }
}


static INLINE int32_t is_mv_valid(const MV *mv) {
    return mv->row > MV_LOW && mv->row < MV_UPP && mv->col > MV_LOW &&
        mv->col < MV_UPP;
}

void av1_encode_mv(
    PictureParentControlSet_t   *pcsPtr,
    aom_writer                  *ecWriter,
    const MV *mv,
    const MV *ref,
    nmv_context *mvctx,
    int32_t usehp) {



    if (is_mv_valid(mv) == 0)
        printf("Corrupted MV\n");


    int32_t diff[2] = { mv->row - ref->row, mv->col - ref->col };
    const MV_JOINT_TYPE j = av1_get_mv_joint_diff(diff);

    if (pcsPtr->cur_frame_force_integer_mv) {
        usehp = MV_SUBPEL_NONE;
    }
    aom_write_symbol(ecWriter, j, mvctx->joints_cdf, MV_JOINTS);
    if (mv_joint_vertical(j))
        encode_mv_component(ecWriter, diff[0], &mvctx->comps[0], (MvSubpelPrecision)usehp);

    if (mv_joint_horizontal(j))
        encode_mv_component(ecWriter, diff[1], &mvctx->comps[1], (MvSubpelPrecision)usehp);

    // If auto_mv_step_size is enabled then keep track of the largest
    // motion vector component used.
    //if (cpi->sf.mv.auto_mv_step_size) {
    //    uint32_t maxv = AOMMAX(abs(mv->row), abs(mv->col)) >> 3;
    //    cpi->max_mv_magnitude = AOMMAX(maxv, cpi->max_mv_magnitude);
    //}
}

///InterpFilter av1_extract_interp_filter(uint32_t filters,
//    int32_t x_filter) {
//    return (InterpFilter)((filters >> (x_filter ? 16 : 0)) & 0xffff);
//}
#define INTER_FILTER_COMP_OFFSET (SWITCHABLE_FILTERS + 1)
#define INTER_FILTER_DIR_OFFSET ((SWITCHABLE_FILTERS + 1) * 2)

int32_t av1_get_pred_context_switchable_interp(
    NeighborArrayUnit_t     *ref_frame_type_neighbor_array,
    MvReferenceFrame rf0,
    MvReferenceFrame rf1,
    NeighborArrayUnit32_t     *interpolation_type_neighbor_array,

    uint32_t cu_origin_x,
    uint32_t cu_origin_y,
    //const MacroBlockD *xd,
    int32_t dir
) {

    uint32_t interpolationTypeLeftNeighborIndex = GetNeighborArrayUnitLeftIndex32(
        interpolation_type_neighbor_array,
        cu_origin_y);
    uint32_t interpolationTypeTopNeighborIndex = GetNeighborArrayUnitTopIndex32(
        interpolation_type_neighbor_array,
        cu_origin_x);

    uint32_t rfLeftNeighborIndex = GetNeighborArrayUnitLeftIndex(
        ref_frame_type_neighbor_array,
        cu_origin_y);
    uint32_t rfTopNeighborIndex = GetNeighborArrayUnitTopIndex(
        ref_frame_type_neighbor_array,
        cu_origin_x);

    //const MbModeInfo *const mbmi = xd->mi[0];
    const int32_t ctx_offset =
        (rf1 > INTRA_FRAME) * INTER_FILTER_COMP_OFFSET;
    MvReferenceFrame ref_frame =
        (dir < 2) ? rf0 : rf1;
    // Note:
    // The mode info data structure has a one element border above and to the
    // left of the entries corresponding to real macroblocks.
    // The prediction flags in these dummy entries are initialized to 0.
    int32_t filter_type_ctx = ctx_offset + (dir & 0x01) * INTER_FILTER_DIR_OFFSET;
    int32_t left_type = SWITCHABLE_FILTERS;
    int32_t above_type = SWITCHABLE_FILTERS;

    if (cu_origin_x != 0 /*&& interpolation_type_neighbor_array->leftArray[interpolationTypeLeftNeighborIndex] != SWITCHABLE_FILTERS*/) {

        MvReferenceFrame rf_left[2];
        av1_set_ref_frame(rf_left, ref_frame_type_neighbor_array->leftArray[rfLeftNeighborIndex]);
        uint32_t leftNeigh = (uint32_t)interpolation_type_neighbor_array->leftArray[interpolationTypeLeftNeighborIndex];
        left_type = (rf_left[0] == ref_frame || rf_left[1] == ref_frame) ? av1_extract_interp_filter(leftNeigh, dir & 0x01) : SWITCHABLE_FILTERS;
    }
    //get_ref_filter_type(xd->mi[-1], xd, dir, ref_frame);



    if (cu_origin_y != 0 /*&& interpolation_type_neighbor_array->topArray[interpolationTypeTopNeighborIndex] != SWITCHABLE_FILTERS*/) {
        MvReferenceFrame rf_above[2];
        av1_set_ref_frame(rf_above, ref_frame_type_neighbor_array->topArray[rfTopNeighborIndex]);
        uint32_t aboveNeigh = (uint32_t)interpolation_type_neighbor_array->topArray[interpolationTypeTopNeighborIndex];

        above_type = (rf_above[0] == ref_frame || rf_above[1] == ref_frame) ? av1_extract_interp_filter(aboveNeigh, dir & 0x01) : SWITCHABLE_FILTERS;
        //get_ref_filter_type(xd->mi[-xd->mi_stride], xd, dir, ref_frame);
    }

    if (left_type == above_type) {
        filter_type_ctx += left_type;
    }
    else if (left_type == SWITCHABLE_FILTERS) {
        assert(above_type != SWITCHABLE_FILTERS);
        filter_type_ctx += above_type;
    }
    else if (above_type == SWITCHABLE_FILTERS) {
        assert(left_type != SWITCHABLE_FILTERS);
        filter_type_ctx += left_type;
    }
    else {
        filter_type_ctx += SWITCHABLE_FILTERS;
    }
    //  printf("\t %d\t %d\t %d",filter_type_ctx,left_type ,above_type);

    return filter_type_ctx;
}


/*INLINE*/ int32_t is_nontrans_global_motion_EC(
    MvReferenceFrame             rf0,
    MvReferenceFrame             rf1,
    CodingUnit_t                  *cu_ptr,
    BlockSize sb_type,
    PictureParentControlSet_t   *pcsPtr) {
    int32_t ref;

    // First check if all modes are GLOBALMV
    if (cu_ptr->pred_mode != GLOBALMV && cu_ptr->pred_mode != GLOBAL_GLOBALMV)
        return 0;

    if (MIN(mi_size_wide[sb_type], mi_size_high[sb_type]) < 2)
        return 0;

    // Now check if all global motion is non translational
    for (ref = 0; ref < 1 + cu_ptr->prediction_unit_array->is_compound; ++ref) {
        if (pcsPtr->global_motion[ref ? rf1 : rf0].wmtype == TRANSLATION)
            return 0;
    }
    return 1;
}
static int32_t av1_is_interp_needed(
    MvReferenceFrame            rf0,
    MvReferenceFrame            rf1,
    CodingUnit_t               *cu_ptr,
    BlockSize                   bsize,
    PictureParentControlSet_t  *pcsPtr)
{
    if (cu_ptr->skip_flag)
        return 0;

    if (cu_ptr->prediction_unit_array[0].motion_mode == WARPED_CAUSAL)
        return 0;

    if (is_nontrans_global_motion_EC(rf0, rf1, cu_ptr, bsize, pcsPtr))
        return 0;

    return 1;
}

void write_mb_interp_filter(
    NeighborArrayUnit_t     *ref_frame_type_neighbor_array,
    BlockSize bsize,
    MvReferenceFrame rf0,
    MvReferenceFrame rf1,
    PictureParentControlSet_t   *pcsPtr,
    aom_writer                  *ecWriter,
    CodingUnit_t             *cu_ptr,
    EntropyCoder_t *entropy_coder_ptr,
    NeighborArrayUnit32_t     *interpolation_type_neighbor_array,
    uint32_t blkOriginX,
    uint32_t blkOriginY)
{
    Av1Common *const cm = pcsPtr->av1_cm; //&cpi->common;
    //const MbModeInfo *const mbmi = xd->mi[0];
    //FRAME_CONTEXT *ec_ctx = xd->tile_ctx;

    if (!av1_is_interp_needed(rf0, rf1, cu_ptr, bsize, pcsPtr)) {

        /* assert(mbmi->interp_filters ==
                av1_broadcast_interp_filter(
                    av1_unswitchable_filter(cm->interp_filter)));*/
        return;
    }
    if (cm->interp_filter == SWITCHABLE) {
        int32_t dir;
        for (dir = 0; dir < 2; ++dir) {
            // printf("\nP:%d\tX: %d\tY:%d\t %d",(pcsPtr)->picture_number,blkOriginX,blkOriginY ,((ecWriter)->ec).rng);
            const int32_t ctx = av1_get_pred_context_switchable_interp(
                ref_frame_type_neighbor_array,

                rf0,
                rf1,
                interpolation_type_neighbor_array,

                blkOriginX,
                blkOriginY,
                //xd,
                dir
            );
            InterpFilter filter = av1_extract_interp_filter(cu_ptr->interp_filters, dir);
            ASSERT(ctx < SWITCHABLE_FILTER_CONTEXTS);
            aom_write_symbol(ecWriter, filter, /*ec_ctx*/entropy_coder_ptr->fc->switchable_interp_cdf[ctx],
                SWITCHABLE_FILTERS);

            // ++pcsPtr->interp_filter_selected[0][filter];
            if (pcsPtr->sequence_control_set_ptr->enable_dual_filter == 0) return;
        }
    }
}

static void WriteInterCompoundMode(
    FRAME_CONTEXT       *frameContext,
    aom_writer          *ecWriter,
    PredictionMode     mode,
    const int16_t       mode_ctx) {
    assert(is_inter_compound_mode(mode));
    aom_write_symbol(ecWriter, INTER_COMPOUND_OFFSET(mode),
        frameContext->inter_compound_mode_cdf[mode_ctx],
        INTER_COMPOUND_MODES);

}

int32_t Av1GetReferenceModeContext(
    uint32_t                  cu_origin_x,
    uint32_t                  cu_origin_y,
    NeighborArrayUnit_t    *mode_type_neighbor_array,
    NeighborArrayUnit_t    *inter_pred_dir_neighbor_array)
{


    uint32_t modeTypeLeftNeighborIndex = GetNeighborArrayUnitLeftIndex(
        mode_type_neighbor_array,
        cu_origin_y);
    uint32_t modeTypeTopNeighborIndex = GetNeighborArrayUnitTopIndex(
        mode_type_neighbor_array,
        cu_origin_x);

    int32_t ctx = 0;

    // Note:
    // The mode info data structure has a one element border above and to the
    // left of the entries corresponding to real macroblocks.
    // The prediction flags in these dummy entries are initialized to 0.
    if (mode_type_neighbor_array->leftArray[modeTypeLeftNeighborIndex] != (uint8_t)INVALID_MODE && mode_type_neighbor_array->topArray[modeTypeTopNeighborIndex] != (uint8_t)INVALID_MODE) {  // both edges available
        const int32_t topIntra = (mode_type_neighbor_array->topArray[modeTypeTopNeighborIndex] == (uint8_t)INTRA_MODE);
        const int32_t leftIntra = (mode_type_neighbor_array->leftArray[modeTypeLeftNeighborIndex] == (uint8_t)INTRA_MODE);
        //if (has_above && has_left) {  // both edges available
        if (!(inter_pred_dir_neighbor_array->topArray[modeTypeTopNeighborIndex] == BI_PRED && !topIntra) && !(inter_pred_dir_neighbor_array->leftArray[modeTypeLeftNeighborIndex] == BI_PRED && !leftIntra)) {
            // neither edge uses comp pred (0/1)
            ctx = (inter_pred_dir_neighbor_array->topArray[modeTypeTopNeighborIndex] == UNI_PRED_LIST_1) ^
                (inter_pred_dir_neighbor_array->leftArray[modeTypeLeftNeighborIndex] == UNI_PRED_LIST_1);
        }
        else if (!(inter_pred_dir_neighbor_array->topArray[modeTypeTopNeighborIndex] == BI_PRED && !topIntra)/*has_second_ref(above_mbmi)*/) {
            // one of two edges uses comp pred (2/3)
            ctx = 2 + ((inter_pred_dir_neighbor_array->topArray[modeTypeTopNeighborIndex] == UNI_PRED_LIST_1) ||
                topIntra);
        }
        else if (!(inter_pred_dir_neighbor_array->leftArray[modeTypeLeftNeighborIndex] == BI_PRED && !leftIntra)/*has_second_ref(left_mbmi)*/) {
            // one of two edges uses comp pred (2/3)
            ctx = 2 + ((inter_pred_dir_neighbor_array->leftArray[modeTypeLeftNeighborIndex] == UNI_PRED_LIST_1) ||
                leftIntra);
        }
        else {  // both edges use comp pred (4){

            ctx = 4;
        }
    }
    else if (mode_type_neighbor_array->leftArray[modeTypeLeftNeighborIndex] != (uint8_t)INVALID_MODE) {  // one edge available

        if (!(inter_pred_dir_neighbor_array->leftArray[modeTypeLeftNeighborIndex] == BI_PRED && mode_type_neighbor_array->leftArray[modeTypeLeftNeighborIndex] != (uint8_t)INTRA_MODE)/*has_second_ref(edge_mbmi)*/) {
            // edge does not use comp pred (0/1)
            ctx = (inter_pred_dir_neighbor_array->leftArray[modeTypeLeftNeighborIndex] == UNI_PRED_LIST_1);
        }
        else {
            // edge uses comp pred (3)
            ctx = 3;
        }
    }
    else if (mode_type_neighbor_array->topArray[modeTypeTopNeighborIndex] != (uint8_t)INVALID_MODE) {  // one edge available

        if (!(inter_pred_dir_neighbor_array->topArray[modeTypeTopNeighborIndex] == BI_PRED && mode_type_neighbor_array->topArray[modeTypeTopNeighborIndex] != (uint8_t)INTRA_MODE)/*has_second_ref(edge_mbmi)*/) {
            // edge does not use comp pred (0/1)
            ctx = (inter_pred_dir_neighbor_array->topArray[modeTypeTopNeighborIndex] == UNI_PRED_LIST_1);
        }
        else {
            // edge uses comp pred (3)
            ctx = 3;
        }
    }
    else {  // no edges available (1)
        ctx = 1;
    }

    assert(ctx >= 0 && ctx < COMP_INTER_CONTEXTS);
    return ctx;
}


int32_t Av1GetCompReferenceTypeContext(
    uint32_t                  cu_origin_x,
    uint32_t                  cu_origin_y,
    NeighborArrayUnit_t    *mode_type_neighbor_array,
    NeighborArrayUnit_t     *inter_pred_dir_neighbor_array) {

    int32_t pred_context = 0;

    uint32_t modeTypeLeftNeighborIndex = GetNeighborArrayUnitLeftIndex(
        mode_type_neighbor_array,
        cu_origin_y);
    uint32_t modeTypeTopNeighborIndex = GetNeighborArrayUnitTopIndex(
        mode_type_neighbor_array,
        cu_origin_x);


    if (mode_type_neighbor_array->leftArray[modeTypeLeftNeighborIndex] != (uint8_t)INVALID_MODE && mode_type_neighbor_array->topArray[modeTypeTopNeighborIndex] != (uint8_t)INVALID_MODE) {  // both edges available
        const int32_t above_intra = (mode_type_neighbor_array->topArray[modeTypeTopNeighborIndex] == (uint8_t)INTRA_MODE);
        const int32_t left_intra = (mode_type_neighbor_array->leftArray[modeTypeLeftNeighborIndex] == (uint8_t)INTRA_MODE);

        if (above_intra && left_intra) {  // intra/intra
            pred_context = 2;
        }
        else if (left_intra) {  // Intra & Inter. Left is Intra, check Top

            if (inter_pred_dir_neighbor_array->topArray[modeTypeTopNeighborIndex] != BI_PRED)  // single pred
                pred_context = 2;
            else  // comp pred
                pred_context = 1 + 2 * 0/* has_uni_comp_refs(inter_mbmi)*/;
        }
        else if (above_intra) {  // Intra & Inter. Top is Intra, check Left

            if (inter_pred_dir_neighbor_array->leftArray[modeTypeLeftNeighborIndex] != BI_PRED)  // single pred
                pred_context = 2;
            else  // comp pred
                pred_context = 1 + 2 * 0/* has_uni_comp_refs(inter_mbmi)*/;
        }
        else {  // inter/inter
            const int32_t a_sg = (inter_pred_dir_neighbor_array->topArray[modeTypeTopNeighborIndex] != BI_PRED);
            const int32_t l_sg = (inter_pred_dir_neighbor_array->leftArray[modeTypeLeftNeighborIndex] != BI_PRED);
            //const MvReferenceFrame frfa = above_mbmi->ref_frame[0];
            //const MvReferenceFrame frfl = left_mbmi->ref_frame[0];

            if (a_sg && l_sg) {  // single/single
                pred_context = 1 + 2 * (!((inter_pred_dir_neighbor_array->topArray[modeTypeTopNeighborIndex] == UNI_PRED_LIST_1) ^
                    (inter_pred_dir_neighbor_array->leftArray[modeTypeLeftNeighborIndex] == UNI_PRED_LIST_1)));
            }
            else if (l_sg || a_sg) {  // single/comp
                const int32_t uni_rfc =
                    a_sg ? 0 /*has_uni_comp_refs(left_mbmi) */ : 0 /*has_uni_comp_refs(above_mbmi)*/;

                if (!uni_rfc)  // comp bidir
                    pred_context = 1;
                else  // comp unidir
                    pred_context = 3 + (!((inter_pred_dir_neighbor_array->topArray[modeTypeTopNeighborIndex] == UNI_PRED_LIST_1) ^
                    (inter_pred_dir_neighbor_array->leftArray[modeTypeLeftNeighborIndex] == UNI_PRED_LIST_1)));
            }
            else {  // comp/comp
                const int32_t a_uni_rfc = 0;// has_uni_comp_refs(above_mbmi);
                const int32_t l_uni_rfc = 0;// has_uni_comp_refs(left_mbmi);
                pred_context = 0;

                if (!a_uni_rfc && !l_uni_rfc)  // bidir/bidir
                    pred_context = 0;
                else if (!a_uni_rfc || !l_uni_rfc)  // unidir/bidir
                    pred_context = 2;
                else  // unidir/unidir
                    pred_context =
                    3 + (!(0 ^ 0));
            }
        }
    }
    else if (mode_type_neighbor_array->leftArray[modeTypeLeftNeighborIndex] != (uint8_t)INVALID_MODE) {  // one edge available


        if (mode_type_neighbor_array->leftArray[modeTypeLeftNeighborIndex] == (uint8_t)INTRA_MODE) {  // intra
            pred_context = 2;
        }
        else {                           // inter
            if (inter_pred_dir_neighbor_array->leftArray[modeTypeLeftNeighborIndex] != BI_PRED)  // single pred
                pred_context = 2;
            else  // comp pred
                pred_context = 4 * 0/*has_uni_comp_refs(edge_mbmi)*/;
        }
    }
    else if (mode_type_neighbor_array->topArray[modeTypeTopNeighborIndex] != (uint8_t)INVALID_MODE) {  // one edge available


        if (mode_type_neighbor_array->topArray[modeTypeTopNeighborIndex] == (uint8_t)INTRA_MODE) {  // intra
            pred_context = 2;
        }
        else {                           // inter
            if (inter_pred_dir_neighbor_array->topArray[modeTypeTopNeighborIndex] != BI_PRED)  // single pred
                pred_context = 2;
            else  // comp pred
                pred_context = 4 * 0/*has_uni_comp_refs(edge_mbmi)*/;
        }
    }
    else {  // no edges available
        pred_context = 2;
    }

    //assert(pred_context >= 0 && pred_context < COMP_REF_TYPE_CONTEXTS);
    return pred_context;
}

void Av1CollectNeighborsRefCounts(
    CodingUnit_t            *cu_ptr,
    uint32_t                   cu_origin_x,
    uint32_t                   cu_origin_y,
    NeighborArrayUnit_t     *mode_type_neighbor_array,
    NeighborArrayUnit_t     *inter_pred_dir_neighbor_array,
    NeighborArrayUnit_t     *ref_frame_type_neighbor_array) {


    av1_zero(cu_ptr->av1xd->neighbors_ref_counts);

    uint8_t *const ref_counts = cu_ptr->av1xd->neighbors_ref_counts;


    uint32_t modeTypeLeftNeighborIndex = GetNeighborArrayUnitLeftIndex(
        mode_type_neighbor_array,
        cu_origin_y);
    uint32_t modeTypeTopNeighborIndex = GetNeighborArrayUnitTopIndex(
        mode_type_neighbor_array,
        cu_origin_x);


    const int32_t topInter = (mode_type_neighbor_array->topArray[modeTypeTopNeighborIndex] == (uint8_t)INTER_MODE);
    const int32_t leftInter = (mode_type_neighbor_array->leftArray[modeTypeLeftNeighborIndex] == (uint8_t)INTER_MODE);

    const int32_t topInImage = (mode_type_neighbor_array->topArray[modeTypeTopNeighborIndex] != (uint8_t)INVALID_MODE);
    const int32_t leftInImage = (mode_type_neighbor_array->leftArray[modeTypeLeftNeighborIndex] != (uint8_t)INVALID_MODE);

    // Above neighbor
    if (topInImage && topInter) {
        MvReferenceFrame topRfType[2];
        av1_set_ref_frame(topRfType, ref_frame_type_neighbor_array->topArray[modeTypeTopNeighborIndex]);
        switch (inter_pred_dir_neighbor_array->topArray[modeTypeTopNeighborIndex]) {
        case UNI_PRED_LIST_0:
            ref_counts[topRfType[0]]++;

            break;
        case UNI_PRED_LIST_1:
            ref_counts[topRfType[0]]++;
            break;

        case BI_PRED:
            ref_counts[topRfType[0]]++;
            ref_counts[topRfType[1]]++;
            break;
        default:
            printf("ERROR[AN]: Invalid Pred Direction\n");
            break;
        }
    }

    // Left neighbor
    if (leftInImage && leftInter) {
        MvReferenceFrame leftRfType[2];
        av1_set_ref_frame(leftRfType, ref_frame_type_neighbor_array->leftArray[modeTypeLeftNeighborIndex]);
        switch (inter_pred_dir_neighbor_array->leftArray[modeTypeLeftNeighborIndex]) {
        case UNI_PRED_LIST_0:
            ref_counts[leftRfType[0]]++;

            break;
        case UNI_PRED_LIST_1:
            ref_counts[leftRfType[0]]++;
            break;

        case BI_PRED:
            ref_counts[leftRfType[0]]++;
            ref_counts[leftRfType[1]]++;
            break;
        default:
            printf("ERROR[AN]: Invalid Pred Direction\n");
            break;
        }
    }
}

#define WRITE_REF_BIT(bname, pname) \
  aom_write_symbol(ecWriter, bname, av1_get_pred_cdf_##pname(), 2)
/***************************************************************************************/

// == Common context functions for both comp and single ref ==
//
// Obtain contexts to signal a reference frame to be either LAST/LAST2 or
// LAST3/GOLDEN.
static int32_t get_pred_context_ll2_or_l3gld(const MacroBlockD *xd) {
    const uint8_t *const ref_counts = &xd->neighbors_ref_counts[0];

    // Count of LAST + LAST2
    const int32_t last_last2_count = ref_counts[LAST_FRAME] + ref_counts[LAST2_FRAME];
    // Count of LAST3 + GOLDEN
    const int32_t last3_gld_count =
        ref_counts[LAST3_FRAME] + ref_counts[GOLDEN_FRAME];

    const int32_t pred_context = (last_last2_count == last3_gld_count)
        ? 1
        : ((last_last2_count < last3_gld_count) ? 0 : 2);

    assert(pred_context >= 0 && pred_context < REF_CONTEXTS);
    return pred_context;
}

// Obtain contexts to signal a reference frame to be either LAST or LAST2.
static int32_t get_pred_context_last_or_last2(const MacroBlockD *xd) {
    const uint8_t *const ref_counts = &xd->neighbors_ref_counts[0];

    // Count of LAST
    const int32_t last_count = ref_counts[LAST_FRAME];
    // Count of LAST2
    const int32_t last2_count = ref_counts[LAST2_FRAME];

    const int32_t pred_context =
        (last_count == last2_count) ? 1 : ((last_count < last2_count) ? 0 : 2);

    assert(pred_context >= 0 && pred_context < REF_CONTEXTS);
    return pred_context;
}

// Obtain contexts to signal a reference frame to be either LAST3 or GOLDEN.
static int32_t get_pred_context_last3_or_gld(const MacroBlockD *xd) {
    const uint8_t *const ref_counts = &xd->neighbors_ref_counts[0];

    // Count of LAST3
    const int32_t last3_count = ref_counts[LAST3_FRAME];
    // Count of GOLDEN
    const int32_t gld_count = ref_counts[GOLDEN_FRAME];

    const int32_t pred_context =
        (last3_count == gld_count) ? 1 : ((last3_count < gld_count) ? 0 : 2);

    assert(pred_context >= 0 && pred_context < REF_CONTEXTS);
    return pred_context;
}

// Obtain contexts to signal a reference frame be either BWDREF/ALTREF2, or
// ALTREF.
static int32_t get_pred_context_brfarf2_or_arf(const MacroBlockD *xd) {
    const uint8_t *const ref_counts = &xd->neighbors_ref_counts[0];

    // Counts of BWDREF, ALTREF2, or ALTREF frames (B, A2, or A)
    const int32_t brfarf2_count =
        ref_counts[BWDREF_FRAME] + ref_counts[ALTREF2_FRAME];
    const int32_t arf_count = ref_counts[ALTREF_FRAME];

    const int32_t pred_context =
        (brfarf2_count == arf_count) ? 1 : ((brfarf2_count < arf_count) ? 0 : 2);

    assert(pred_context >= 0 && pred_context < REF_CONTEXTS);
    return pred_context;
}

// Obtain contexts to signal a reference frame be either BWDREF or ALTREF2.
static int32_t get_pred_context_brf_or_arf2(const MacroBlockD *xd) {
    const uint8_t *const ref_counts = &xd->neighbors_ref_counts[0];

    // Count of BWDREF frames (B)
    const int32_t brf_count = ref_counts[BWDREF_FRAME];
    // Count of ALTREF2 frames (A2)
    const int32_t arf2_count = ref_counts[ALTREF2_FRAME];

    const int32_t pred_context =
        (brf_count == arf2_count) ? 1 : ((brf_count < arf2_count) ? 0 : 2);

    assert(pred_context >= 0 && pred_context < REF_CONTEXTS);
    return pred_context;
}

// == Context functions for comp ref ==
//
// Returns a context number for the given MB prediction signal
// Signal the first reference frame for a compound mode be either
// GOLDEN/LAST3, or LAST/LAST2.
int32_t av1_get_pred_context_comp_ref_p(const MacroBlockD *xd) {
    return get_pred_context_ll2_or_l3gld(xd);
}

// Returns a context number for the given MB prediction signal
// Signal the first reference frame for a compound mode be LAST,
// conditioning on that it is known either LAST/LAST2.
int32_t av1_get_pred_context_comp_ref_p1(const MacroBlockD *xd) {
    return get_pred_context_last_or_last2(xd);
}

// Returns a context number for the given MB prediction signal
// Signal the first reference frame for a compound mode be GOLDEN,
// conditioning on that it is known either GOLDEN or LAST3.
int32_t av1_get_pred_context_comp_ref_p2(const MacroBlockD *xd) {
    return get_pred_context_last3_or_gld(xd);
}

// Signal the 2nd reference frame for a compound mode be either
// ALTREF, or ALTREF2/BWDREF.
int32_t av1_get_pred_context_comp_bwdref_p(const MacroBlockD *xd) {
    return get_pred_context_brfarf2_or_arf(xd);
}

// Signal the 2nd reference frame for a compound mode be either
// ALTREF2 or BWDREF.
int32_t av1_get_pred_context_comp_bwdref_p1(const MacroBlockD *xd) {
    return get_pred_context_brf_or_arf2(xd);
}

// == Context functions for single ref ==
//
// For the bit to signal whether the single reference is a forward reference
// frame or a backward reference frame.
int32_t av1_get_pred_context_single_ref_p1(const MacroBlockD *xd) {
    const uint8_t *const ref_counts = &xd->neighbors_ref_counts[0];

    // Count of forward reference frames
    const int32_t fwd_count = ref_counts[LAST_FRAME] + ref_counts[LAST2_FRAME] +
        ref_counts[LAST3_FRAME] + ref_counts[GOLDEN_FRAME];
    // Count of backward reference frames
    const int32_t bwd_count = ref_counts[BWDREF_FRAME] + ref_counts[ALTREF2_FRAME] +
        ref_counts[ALTREF_FRAME];

    const int32_t pred_context =
        (fwd_count == bwd_count) ? 1 : ((fwd_count < bwd_count) ? 0 : 2);

    assert(pred_context >= 0 && pred_context < REF_CONTEXTS);
    return pred_context;
}

// For the bit to signal whether the single reference is ALTREF_FRAME or
// non-ALTREF backward reference frame, knowing that it shall be either of
// these 2 choices.
int32_t av1_get_pred_context_single_ref_p2(const MacroBlockD *xd) {
    return get_pred_context_brfarf2_or_arf(xd);
}

// For the bit to signal whether the single reference is LAST3/GOLDEN or
// LAST2/LAST, knowing that it shall be either of these 2 choices.
int32_t av1_get_pred_context_single_ref_p3(const MacroBlockD *xd) {
    return get_pred_context_ll2_or_l3gld(xd);
}

// For the bit to signal whether the single reference is LAST2_FRAME or
// LAST_FRAME, knowing that it shall be either of these 2 choices.
int32_t av1_get_pred_context_single_ref_p4(const MacroBlockD *xd) {
    return get_pred_context_last_or_last2(xd);
}

// For the bit to signal whether the single reference is GOLDEN_FRAME or
// LAST3_FRAME, knowing that it shall be either of these 2 choices.
int32_t av1_get_pred_context_single_ref_p5(const MacroBlockD *xd) {
    return get_pred_context_last3_or_gld(xd);
}

// For the bit to signal whether the single reference is ALTREF2_FRAME or
// BWDREF_FRAME, knowing that it shall be either of these 2 choices.
int32_t av1_get_pred_context_single_ref_p6(const MacroBlockD *xd) {
    return get_pred_context_brf_or_arf2(xd);
}
/***************************************************************************************/



// This function encodes the reference frame
static void WriteRefFrames(
    FRAME_CONTEXT               *frameContext,
    aom_writer                  *ecWriter,
    PictureParentControlSet_t   *pcsPtr,
    CodingUnit_t                *cu_ptr,
    BlockSize                   bsize,
    uint32_t                       cu_origin_x,
    uint32_t                       cu_origin_y,
    NeighborArrayUnit_t         *mode_type_neighbor_array,
    NeighborArrayUnit_t         *inter_pred_dir_neighbor_array)
{
    //const MbModeInfo *const mbmi = &xd->mi[0]->mbmi;
    const int32_t is_compound = (cu_ptr->prediction_unit_array[0].inter_pred_direction_index == BI_PRED); //is_compound;// has_second_ref(mbmi);
    //const int32_t segment_id = mbmi->segment_id;

    // If segment level coding of this signal is disabled...
    // or the segment allows multiple reference frame options
    /*if (segfeature_active(&cm->seg, segment_id, SEG_LVL_REF_FRAME)) {
    assert(!is_compound);
    assert(mbmi->ref_frame[0] ==
    get_segdata(&cm->seg, segment_id, SEG_LVL_REF_FRAME));
    }
    else if (segfeature_active(&cm->seg, segment_id, SEG_LVL_SKIP) ||
    segfeature_active(&cm->seg, segment_id, SEG_LVL_GLOBALMV)) {
    assert(!is_compound);
    assert(mbmi->ref_frame[0] == LAST_FRAME);
    }
    else*/ {
    // does the feature use compound prediction or not
    // (if not specified at the frame/segment level)
        if (pcsPtr->reference_mode == REFERENCE_MODE_SELECT) {

            if (is_comp_ref_allowed(bsize)) {
                int32_t context = 0;
                context = Av1GetReferenceModeContext(
                    cu_origin_x,
                    cu_origin_y,
                    mode_type_neighbor_array,
                    inter_pred_dir_neighbor_array);
                aom_write_symbol(ecWriter, is_compound, frameContext->comp_inter_cdf[context], 2);


            }
        }
        else {
            assert((!is_compound) == (pcsPtr->reference_mode == SINGLE_REFERENCE));
        }


        int32_t context = 0;
        if (is_compound) {
            const COMP_REFERENCE_TYPE comp_ref_type = /*has_uni_comp_refs(mbmi)
                                                      ? UNIDIR_COMP_REFERENCE
                                                      : */BIDIR_COMP_REFERENCE;
            MvReferenceFrame refType[2];
            av1_set_ref_frame(refType, cu_ptr->prediction_unit_array[0].ref_frame_type);


            context = Av1GetCompReferenceTypeContext(
                cu_origin_x,
                cu_origin_y,
                mode_type_neighbor_array,
                inter_pred_dir_neighbor_array);

            aom_write_symbol(ecWriter, comp_ref_type, frameContext->comp_ref_type_cdf[context],
                2);

            if (comp_ref_type == UNIDIR_COMP_REFERENCE) {
                printf("ERROR[AN]: UNIDIR_COMP_REFERENCE not supported\n");
                //const int32_t bit = mbmi->ref_frame[0] == BWDREF_FRAME;
                //WRITE_REF_BIT(bit, uni_comp_ref_p);

                //if (!bit) {
                //    assert(mbmi->ref_frame[0] == LAST_FRAME);
                //    const int32_t bit1 = mbmi->ref_frame[1] == LAST3_FRAME ||
                //        mbmi->ref_frame[1] == GOLDEN_FRAME;
                //    WRITE_REF_BIT(bit1, uni_comp_ref_p1);
                //    if (bit1) {
                //        const int32_t bit2 = mbmi->ref_frame[1] == GOLDEN_FRAME;
                //        WRITE_REF_BIT(bit2, uni_comp_ref_p2);
                //    }
                //}
                //else {
                //    assert(mbmi->ref_frame[1] == ALTREF_FRAME);
                //}

                //return;
            }

            assert(comp_ref_type == BIDIR_COMP_REFERENCE);

            const int32_t bit = (refType[0] == GOLDEN_FRAME ||
                refType[0] == LAST3_FRAME);

            context = av1_get_pred_context_comp_ref_p(cu_ptr->av1xd);
            aom_write_symbol(ecWriter, bit, frameContext->comp_ref_cdf[context][0],
                2);
            //            WRITE_REF_BIT(bit, comp_ref_p);


            if (!bit) {
                const int32_t bit1 = (refType[0] == LAST2_FRAME);
                context = av1_get_pred_context_comp_ref_p1(cu_ptr->av1xd);
                aom_write_symbol(ecWriter, bit1, frameContext->comp_ref_cdf[context][1],
                    2);
                //WRITE_REF_BIT(bit1, comp_ref_p1);
            }
            else {
                const int32_t bit2 = (refType[0] == GOLDEN_FRAME);
                context = av1_get_pred_context_comp_ref_p2(cu_ptr->av1xd);
                aom_write_symbol(ecWriter, bit2, frameContext->comp_ref_cdf[context][2],
                    2);

                //WRITE_REF_BIT(bit2, comp_ref_p2);
            }

            const int32_t bit_bwd = (refType[1] == ALTREF_FRAME);
            context = av1_get_pred_context_comp_bwdref_p(cu_ptr->av1xd);
            aom_write_symbol(ecWriter, bit_bwd, frameContext->comp_bwdref_cdf[context][0],
                2);
            //WRITE_REF_BIT(bit_bwd, comp_bwdref_p);

            if (!bit_bwd) {
                context = av1_get_pred_context_comp_bwdref_p1(cu_ptr->av1xd);
                aom_write_symbol(ecWriter, refType[1] == ALTREF2_FRAME, frameContext->comp_bwdref_cdf[context][1],
                    2);
                //WRITE_REF_BIT(mbmi->ref_frame[1] == ALTREF2_FRAME, comp_bwdref_p1);
            }

        }
        else {
            const int32_t bit0 = (cu_ptr->prediction_unit_array[0].ref_frame_type <= ALTREF_FRAME &&
                cu_ptr->prediction_unit_array[0].ref_frame_type >= BWDREF_FRAME);//0

            context = av1_get_pred_context_single_ref_p1(cu_ptr->av1xd);
            aom_write_symbol(ecWriter, bit0, frameContext->single_ref_cdf[context][0],
                2);
            //WRITE_REF_BIT(bit0, single_ref_p1);

            if (bit0) {
                const int32_t bit1 = (cu_ptr->prediction_unit_array[0].ref_frame_type == ALTREF_FRAME);
                context = av1_get_pred_context_single_ref_p2(cu_ptr->av1xd);
                aom_write_symbol(ecWriter, bit1, frameContext->single_ref_cdf[context][1],
                    2);
                //WRITE_REF_BIT(bit1, single_ref_p2);

                if (!bit1) {
                    context = av1_get_pred_context_single_ref_p6(cu_ptr->av1xd);
                    aom_write_symbol(ecWriter, cu_ptr->prediction_unit_array[0].ref_frame_type == ALTREF2_FRAME, frameContext->single_ref_cdf[context][5],
                        2);
                    //WRITE_REF_BIT(mbmi->ref_frame[0] == ALTREF2_FRAME, single_ref_p6);
                }
            }
            else {
                const int32_t bit2 = (cu_ptr->prediction_unit_array[0].ref_frame_type == LAST3_FRAME ||
                    cu_ptr->prediction_unit_array[0].ref_frame_type == GOLDEN_FRAME); //0
                context = av1_get_pred_context_single_ref_p3(cu_ptr->av1xd);
                aom_write_symbol(ecWriter, bit2, frameContext->single_ref_cdf[context][2],
                    2);
                //WRITE_REF_BIT(bit2, single_ref_p3);
                if (!bit2) {
                    const int32_t bit3 = (cu_ptr->prediction_unit_array[0].ref_frame_type != LAST_FRAME); //0;
                    context = av1_get_pred_context_single_ref_p4(cu_ptr->av1xd);

                    aom_write_symbol(ecWriter, bit3, frameContext->single_ref_cdf[context][3],
                        2);
                    //WRITE_REF_BIT(bit3, single_ref_p4);
                }
                else {
                    const int32_t bit4 = (cu_ptr->prediction_unit_array[0].ref_frame_type != LAST3_FRAME);
                    context = av1_get_pred_context_single_ref_p5(cu_ptr->av1xd);
                    aom_write_symbol(ecWriter, bit4, frameContext->single_ref_cdf[context][4],
                        2);
                    //WRITE_REF_BIT(bit4, single_ref_p5);
                }
            }
        }
    }
}


static void encode_restoration_mode(PictureParentControlSet_t *pcsPtr,
    struct aom_write_bit_buffer *wb) {

    Av1Common* cm = pcsPtr->av1_cm;

    //printf("ERROR[AN]: encode_restoration_mode might not work. Double check the reference code\n");
    assert(!pcsPtr->all_lossless);
    // move out side of the function
    //if (!cm->seq_params.enable_restoration) return;

    if (pcsPtr->allow_intrabc) return;

    const int32_t num_planes = 3;// av1_num_planes(cm);
    int32_t all_none = 1, chroma_none = 1;
    for (int32_t p = 0; p < num_planes; ++p) {
        //   aom_wb_write_bit(wb, 0);
        //   aom_wb_write_bit(wb, 0);

        RestorationInfo *rsi = &pcsPtr->av1_cm->rst_info[p];

        //if (p==0)
        //   printf("POC:%i Luma restType:%i\n", pcsPtr->picture_number, rsi->frame_restoration_type);
        //if (p == 1)
        //    printf("POC:%i Cb restType:%i\n", pcsPtr->picture_number, rsi->frame_restoration_type);
        //if (p == 2)
        //    printf("POC:%i Cr restType:%i\n", pcsPtr->picture_number, rsi->frame_restoration_type);


        if (rsi->frame_restoration_type != RESTORE_NONE) {
            all_none = 0;
            chroma_none &= (int32_t)(p == 0);
        }
        switch (rsi->frame_restoration_type) {
        case RESTORE_NONE:
            aom_wb_write_bit(wb, 0);
            aom_wb_write_bit(wb, 0);
            break;
        case RESTORE_WIENER:
            aom_wb_write_bit(wb, 1);
            aom_wb_write_bit(wb, 0);
            break;
        case RESTORE_SGRPROJ:
            aom_wb_write_bit(wb, 1);
            aom_wb_write_bit(wb, 1);
            break;
        case RESTORE_SWITCHABLE:
            aom_wb_write_bit(wb, 0);
            aom_wb_write_bit(wb, 1);
            break;
        default: assert(0);
        }
    }
    if (!all_none) {

        //  assert(cm->seq_params.sb_size == BLOCK_64X64 ||
        //      cm->seq_params.sb_size == BLOCK_128X128);
        const int32_t sb_size = pcsPtr->sequence_control_set_ptr->sb_size == BLOCK_128X128 ? 128 : 64;;
        RestorationInfo *rsi = &pcsPtr->av1_cm->rst_info[0];

        assert(rsi->restoration_unit_size >= sb_size);
        assert(RESTORATION_UNITSIZE_MAX == 256);

        if (sb_size == 64) {
            aom_wb_write_bit(wb, rsi->restoration_unit_size > 64);
        }
        if (rsi->restoration_unit_size > 64) {
            aom_wb_write_bit(wb, rsi->restoration_unit_size > 128);
        }
    }

    if (num_planes > 1) {
        int32_t s = 1;// AOMMIN(cm->subsampling_x, cm->subsampling_y);
        if (s && !chroma_none) {
            aom_wb_write_bit(wb, cm->rst_info[1].restoration_unit_size !=
                cm->rst_info[0].restoration_unit_size);
            assert(cm->rst_info[1].restoration_unit_size ==
                cm->rst_info[0].restoration_unit_size ||
                cm->rst_info[1].restoration_unit_size ==
                (cm->rst_info[0].restoration_unit_size >> s));
            assert(cm->rst_info[2].restoration_unit_size ==
                cm->rst_info[1].restoration_unit_size);
        }
        else if (!s) {
            assert(cm->rst_info[1].restoration_unit_size ==
                cm->rst_info[0].restoration_unit_size);
            assert(cm->rst_info[2].restoration_unit_size ==
                cm->rst_info[1].restoration_unit_size);
        }
    }


}


static void encode_loopfilter(PictureParentControlSet_t *pcsPtr, struct aom_write_bit_buffer *wb) {
    assert(!pcsPtr->coded_lossless);
    if (pcsPtr->allow_intrabc) return;
    const int32_t num_planes = 3;// av1_num_planes(pcsPtr);

    struct loopfilter *lf = &pcsPtr->lf;

    // Encode the loop filter level and type
    aom_wb_write_literal(wb, lf->filter_level[0], 6);
    aom_wb_write_literal(wb, lf->filter_level[1], 6);
    if (num_planes > 1) {
        if (lf->filter_level[0] || lf->filter_level[1]) {
            aom_wb_write_literal(wb, lf->filter_level_u, 6);
            aom_wb_write_literal(wb, lf->filter_level_v, 6);
        }
    }
    aom_wb_write_literal(wb, lf->sharpness_level, 3);

    // Write out loop filter deltas applied at the MB level based on mode or
    // ref frame (if they are enabled).
    aom_wb_write_bit(wb, lf->mode_ref_delta_enabled);
    if (lf->mode_ref_delta_enabled) {
        printf("ERROR[AN]: Loop Filter is not supported yet \n");
        /* aom_wb_write_bit(wb, lf->mode_ref_delta_update);
        if (lf->mode_ref_delta_update) {
        const int32_t prime_idx = pcsPtr->primary_ref_frame;
        const int32_t buf_idx =
        prime_idx == PRIMARY_REF_NONE ? -1 : cm->frame_refs[prime_idx].idx;
        int8_t last_ref_deltas[TOTAL_REFS_PER_FRAME];
        if (prime_idx == PRIMARY_REF_NONE || buf_idx < 0) {
        av1_set_default_ref_deltas(last_ref_deltas);
        } else {
        memcpy(last_ref_deltas, cm->buffer_pool->frame_bufs[buf_idx].ref_deltas,
        TOTAL_REFS_PER_FRAME);
        }
        for (i = 0; i < TOTAL_REFS_PER_FRAME; i++) {
        const int32_t delta = lf->ref_deltas[i];
        const int32_t changed = delta != last_ref_deltas[i];
        aom_wb_write_bit(wb, changed);
        if (changed) aom_wb_write_inv_signed_literal(wb, delta, 6);
        }
        int8_t last_mode_deltas[MAX_MODE_LF_DELTAS];
        if (prime_idx == PRIMARY_REF_NONE || buf_idx < 0) {
        av1_set_default_mode_deltas(last_mode_deltas);
        } else {
        memcpy(last_mode_deltas,
        cm->buffer_pool->frame_bufs[buf_idx].mode_deltas,
        MAX_MODE_LF_DELTAS);
        }

        for (i = 0; i < MAX_MODE_LF_DELTAS; i++) {
        const int32_t delta = lf->mode_deltas[i];
        const int32_t changed = delta != last_mode_deltas[i];
        aom_wb_write_bit(wb, changed);
        if (changed) aom_wb_write_inv_signed_literal(wb, delta, 6);
        }
        }*/
    }


}

static void encode_cdef(const PictureParentControlSet_t *pcsPtr, struct aom_write_bit_buffer *wb) {
    //assert(!cm->coded_lossless);
    // moved out side
    //if (!cm->seq_params.enable_cdef) return;

    if (pcsPtr->allow_intrabc) return;

    const int32_t num_planes = 3;// av1_num_planes(pcsPtr);
    int32_t i;
    aom_wb_write_literal(wb, pcsPtr->cdef_pri_damping - 3, 2);
    assert(pcsPtr->cdef_pri_damping == pcsPtr->cdef_sec_damping);
    aom_wb_write_literal(wb, pcsPtr->cdef_bits, 2);
    for (i = 0; i < pcsPtr->nb_cdef_strengths; i++) {
        aom_wb_write_literal(wb, pcsPtr->cdef_strengths[i], CDEF_STRENGTH_BITS);
        if (num_planes > 1)
            aom_wb_write_literal(wb, pcsPtr->cdef_uv_strengths[i], CDEF_STRENGTH_BITS);
    }
}

static void write_delta_q(struct aom_write_bit_buffer *wb, int32_t delta_q) {
    if (delta_q != 0) {
        aom_wb_write_bit(wb, 1);
        aom_wb_write_inv_signed_literal(wb, delta_q, 6);
    }
    else {
        aom_wb_write_bit(wb, 0);
    }
}

static void encode_quantization(const PictureParentControlSet_t *const pcsPtr,
    struct aom_write_bit_buffer *wb) {
    const int32_t num_planes = 3;//av1_num_planes(cm);

    aom_wb_write_literal(wb, pcsPtr->base_qindex, QINDEX_BITS);
    write_delta_q(wb, pcsPtr->y_dc_delta_q);
    if (num_planes > 1) {
        int32_t diff_uv_delta = (pcsPtr->u_dc_delta_q != pcsPtr->v_dc_delta_q) ||
            (pcsPtr->u_ac_delta_q != pcsPtr->v_ac_delta_q);
        if (pcsPtr->separate_uv_delta_q) aom_wb_write_bit(wb, diff_uv_delta);
        write_delta_q(wb, pcsPtr->u_dc_delta_q);
        write_delta_q(wb, pcsPtr->u_ac_delta_q);
        if (diff_uv_delta) {
            write_delta_q(wb, pcsPtr->v_dc_delta_q);
            write_delta_q(wb, pcsPtr->v_ac_delta_q);
        }
    }
    aom_wb_write_bit(wb, pcsPtr->using_qmatrix);
    if (pcsPtr->using_qmatrix) {

        aom_wb_write_literal(wb, pcsPtr->qm_y, QM_LEVEL_BITS);
        aom_wb_write_literal(wb, pcsPtr->qm_u, QM_LEVEL_BITS);
        if (!pcsPtr->separate_uv_delta_q)
            assert(pcsPtr->qm_u == pcsPtr->qm_v);
        else
            aom_wb_write_literal(wb, pcsPtr->qm_v, QM_LEVEL_BITS);

    }
}

static void write_tile_info_max_tile(const PictureParentControlSet_t *const pcsPtr,
    struct aom_write_bit_buffer *wb) {


    Av1Common * cm = pcsPtr->av1_cm;
    aom_wb_write_bit(wb, pcsPtr->uniform_tile_spacing_flag);

    if (pcsPtr->uniform_tile_spacing_flag) {

#if !TILES
        //CHKN: no tiles
        cm->log2_tile_cols = cm->min_log2_tile_cols;
#endif
        // Uniform spaced tiles with power-of-two number of rows and columns
        // tile columns
        int32_t ones = cm->log2_tile_cols - cm->min_log2_tile_cols;
        while (ones--) {
            aom_wb_write_bit(wb, 1);
        }
        if (cm->log2_tile_cols < cm->max_log2_tile_cols) {
            aom_wb_write_bit(wb, 0);
        }

        // rows
#if ! TILES
        //CHKN: no tiles
        cm->log2_tile_rows = cm->min_log2_tile_rows;
#endif
        ones = cm->log2_tile_rows - cm->min_log2_tile_rows;
        while (ones--) {
            aom_wb_write_bit(wb, 1);
        }
        if (cm->log2_tile_rows < cm->max_log2_tile_rows) {
            aom_wb_write_bit(wb, 0);
        }

    }
    else {
        // Explicit tiles with configurable tile widths and heights
        printf("ERROR[AN]:  NON uniform_tile_spacing_flag not supported yet\n");
        //// columns
        //for (i = 0; i < cm->tile_cols; i++) {
        //    size_sb = cm->tile_col_start_sb[i + 1] - cm->tile_col_start_sb[i];
        //    wb_write_uniform(wb, AOMMIN(width_sb, cm->max_tile_width_sb),
        //        size_sb - 1);
        //    width_sb -= size_sb;
        //}
        //assert(width_sb == 0);

        //// rows
        //for (i = 0; i < cm->tile_rows; i++) {
        //    size_sb = cm->tile_row_start_sb[i + 1] - cm->tile_row_start_sb[i];
        //    wb_write_uniform(wb, AOMMIN(height_sb, cm->max_tile_height_sb),
        //        size_sb - 1);
        //    height_sb -= size_sb;
        //}
        //assert(height_sb == 0);
    }
}

#define MAX_TILE_WIDTH (4096)        // Max Tile width in pixels
#define MAX_TILE_AREA (4096 * 2304)  // Maximum tile area in pixels
// Find smallest k>=0 such that (blk_size << k) >= target
static int32_t tile_log2(int32_t blk_size, int32_t target) {
    int32_t k;
    for (k = 0; (blk_size << k) < target; k++) {
    }
    return k;
}
void av1_get_tile_limits(PictureParentControlSet_t * pcsPtr) {

    Av1Common * cm = pcsPtr->av1_cm;

    int32_t mi_cols = ALIGN_POWER_OF_TWO(cm->mi_cols, pcsPtr->sequence_control_set_ptr->mib_size_log2);
    int32_t mi_rows = ALIGN_POWER_OF_TWO(cm->mi_rows, pcsPtr->sequence_control_set_ptr->mib_size_log2);
    int32_t sb_cols = mi_cols >> pcsPtr->sequence_control_set_ptr->mib_size_log2;
    int32_t sb_rows = mi_rows >> pcsPtr->sequence_control_set_ptr->mib_size_log2;

    int32_t sb_size_log2 = pcsPtr->sequence_control_set_ptr->mib_size_log2 + MI_SIZE_LOG2;
    cm->max_tile_width_sb = MAX_TILE_WIDTH >> sb_size_log2;
    int32_t max_tile_area_sb = MAX_TILE_AREA >> (2 * sb_size_log2);

    cm->min_log2_tile_cols = tile_log2(cm->max_tile_width_sb, sb_cols);
    cm->max_log2_tile_cols = tile_log2(1, AOMMIN(sb_cols, MAX_TILE_COLS));
    cm->max_log2_tile_rows = tile_log2(1, AOMMIN(sb_rows, MAX_TILE_ROWS));
    cm->min_log2_tile_rows = 0; // CHKN Tiles 
    cm->min_log2_tiles = tile_log2(max_tile_area_sb, sb_cols * sb_rows);
    cm->min_log2_tiles = AOMMAX(cm->min_log2_tiles, cm->min_log2_tile_cols);
}

#if TILES
void av1_calculate_tile_cols(PictureParentControlSet_t * pcs_ptr) {

    Av1Common *const cm = pcs_ptr->av1_cm;

    int mi_cols = ALIGN_POWER_OF_TWO(cm->mi_cols, pcs_ptr->sequence_control_set_ptr->mib_size_log2);
    int mi_rows = ALIGN_POWER_OF_TWO(cm->mi_rows, pcs_ptr->sequence_control_set_ptr->mib_size_log2);
    int sb_cols = mi_cols >> pcs_ptr->sequence_control_set_ptr->mib_size_log2;
    int sb_rows = mi_rows >> pcs_ptr->sequence_control_set_ptr->mib_size_log2;
    int i;

    if (cm->uniform_tile_spacing_flag) {
        int start_sb;
        int size_sb = ALIGN_POWER_OF_TWO(sb_cols, cm->log2_tile_cols);
        size_sb >>= cm->log2_tile_cols;
        assert(size_sb > 0);
        for (i = 0, start_sb = 0; start_sb < sb_cols; i++) {
            cm->tile_col_start_sb[i] = start_sb;
            start_sb += size_sb;
        }
        cm->tile_cols = i;
        cm->tile_col_start_sb[i] = sb_cols;
        cm->min_log2_tile_rows = AOMMAX(cm->min_log2_tiles - cm->log2_tile_cols, 0);
        cm->max_tile_height_sb = sb_rows >> cm->min_log2_tile_rows;

        cm->tile_width = size_sb << pcs_ptr->sequence_control_set_ptr->mib_size_log2;
        cm->tile_width = AOMMIN(cm->tile_width, cm->mi_cols);
    }
    else {
        int max_tile_area_sb = (sb_rows * sb_cols);
        int widest_tile_sb = 1;
        cm->log2_tile_cols = tile_log2(1, cm->tile_cols);
        for (i = 0; i < cm->tile_cols; i++) {
            int size_sb = cm->tile_col_start_sb[i + 1] - cm->tile_col_start_sb[i];
            widest_tile_sb = AOMMAX(widest_tile_sb, size_sb);
        }
        if (cm->min_log2_tiles)
            max_tile_area_sb >>= (cm->min_log2_tiles + 1);
        
        cm->max_tile_height_sb = AOMMAX(max_tile_area_sb / widest_tile_sb, 1);
    }
}

void av1_calculate_tile_rows(PictureParentControlSet_t * pcsPtr)
{

    Av1Common *const cm = pcsPtr->av1_cm;

    int mi_rows = ALIGN_POWER_OF_TWO(cm->mi_rows, pcsPtr->sequence_control_set_ptr->mib_size_log2);
    int sb_rows = mi_rows >> pcsPtr->sequence_control_set_ptr->mib_size_log2;
    int start_sb, size_sb, i;

    if (cm->uniform_tile_spacing_flag) {
        size_sb = ALIGN_POWER_OF_TWO(sb_rows, cm->log2_tile_rows);
        size_sb >>= cm->log2_tile_rows;
        assert(size_sb > 0);
        for (i = 0, start_sb = 0; start_sb < sb_rows; i++) {
            cm->tile_row_start_sb[i] = start_sb;
            start_sb += size_sb;
        }
        cm->tile_rows = i;
        cm->tile_row_start_sb[i] = sb_rows;

        cm->tile_height = size_sb << pcsPtr->sequence_control_set_ptr->mib_size_log2;
        cm->tile_height = AOMMIN(cm->tile_height, cm->mi_rows);
    }
    else {
        cm->log2_tile_rows = tile_log2(1, cm->tile_rows);
    }
}

 void set_tile_info(PictureParentControlSet_t * pcs_ptr)
{

     /*  Tiling algorithm:
        input : log2_tile_count ==> tile_count = 1<<log2_tile_count
        
        step1) compute pic_size_in_sb
        step2) then round up to the closed n.tile_count. 
        step3) tile_size = rounded_pic_size_in_sb / tile_count.
        step4) we fill tiles of size tile_size until we reach the end of the pic

        Note that: the last tile could have smaller size, and the final number 
        of tiles could be less than tile_count
     */
   
    Av1Common * cm = pcs_ptr->av1_cm;
    int i, start_sb;
    //to connect later if non uniform tile spacing is needed.
    int tile_width_count = 0;
    int tile_height_count= 0;
    int tile_widths[MAX_TILE_COLS] = {0};
    int tile_heights[MAX_TILE_ROWS] = { 0 };

    av1_get_tile_limits(pcs_ptr);


    // configure tile columns
    if (tile_width_count == 0 || tile_height_count == 0) 
    {
        cm->uniform_tile_spacing_flag = 1;
        cm->log2_tile_cols = AOMMAX(pcs_ptr->sequence_control_set_ptr->static_config.tile_columns, cm->min_log2_tile_cols);
        cm->log2_tile_cols = AOMMIN(cm->log2_tile_cols, cm->max_log2_tile_cols);
    }
    else {
        int mi_cols = ALIGN_POWER_OF_TWO(cm->mi_cols, pcs_ptr->sequence_control_set_ptr->mib_size_log2);
        int sb_cols = mi_cols >> pcs_ptr->sequence_control_set_ptr->mib_size_log2;
        int size_sb, j = 0;
        cm->uniform_tile_spacing_flag = 0;
        for (i = 0, start_sb = 0; start_sb < sb_cols && i < MAX_TILE_COLS; i++) {
            cm->tile_col_start_sb[i] = start_sb;
            size_sb = tile_widths[j++];
            if (j >= tile_width_count) j = 0;
            start_sb += AOMMIN(size_sb, cm->max_tile_width_sb);
        }
        cm->tile_cols = i;
        cm->tile_col_start_sb[i] = sb_cols;
    }
    av1_calculate_tile_cols(pcs_ptr);

    // configure tile rows
    if (cm->uniform_tile_spacing_flag) {
        cm->log2_tile_rows = AOMMAX(pcs_ptr->sequence_control_set_ptr->static_config.tile_rows, cm->min_log2_tile_rows);
        cm->log2_tile_rows = AOMMIN(cm->log2_tile_rows, cm->max_log2_tile_rows);
    }
    else {
        int mi_rows = ALIGN_POWER_OF_TWO(cm->mi_rows, pcs_ptr->sequence_control_set_ptr->mib_size_log2);
        int sb_rows = mi_rows >> pcs_ptr->sequence_control_set_ptr->mib_size_log2;
        int size_sb, j = 0;
        for (i = 0, start_sb = 0; start_sb < sb_rows && i < MAX_TILE_ROWS; i++) {
            cm->tile_row_start_sb[i] = start_sb;
            size_sb = tile_heights[j++];
            if (j >= tile_height_count) j = 0;
            start_sb += AOMMIN(size_sb, cm->max_tile_height_sb);
        }
        cm->tile_rows = i;
        cm->tile_row_start_sb[i] = sb_rows;
    }
    av1_calculate_tile_rows(pcs_ptr);
}

 void av1_tile_set_row(TileInfo *tile, PictureParentControlSet_t * pcs_ptr, int row)
 {

     Av1Common *const cm = pcs_ptr->av1_cm;

     assert(row < cm->tile_rows);
     int mi_row_start = cm->tile_row_start_sb[row]    << pcs_ptr->sequence_control_set_ptr->mib_size_log2;
     int mi_row_end  = cm->tile_row_start_sb[row + 1] << pcs_ptr->sequence_control_set_ptr->mib_size_log2;
     tile->tile_row = row;
     tile->mi_row_start = mi_row_start;
     tile->mi_row_end = AOMMIN(mi_row_end, cm->mi_rows);
     assert(tile->mi_row_end > tile->mi_row_start);
 }

 void av1_tile_set_col(TileInfo *tile, PictureParentControlSet_t * pcs_ptr, int col) {

     Av1Common *const cm = pcs_ptr->av1_cm;

     assert(col < cm->tile_cols);
     int mi_col_start = cm->tile_col_start_sb[col] << pcs_ptr->sequence_control_set_ptr->mib_size_log2;
     int mi_col_end = cm->tile_col_start_sb[col + 1]
         << pcs_ptr->sequence_control_set_ptr->mib_size_log2;
     tile->tile_col = col;
     tile->mi_col_start = mi_col_start;
     tile->mi_col_end = AOMMIN(mi_col_end, cm->mi_cols);
     assert(tile->mi_col_end > tile->mi_col_start);
 }
#endif


static void write_tile_info(const PictureParentControlSet_t *const pcs_ptr,
    //struct aom_write_bit_buffer *saved_wb,
    struct aom_write_bit_buffer *wb) {

    av1_get_tile_limits((PictureParentControlSet_t *)pcs_ptr);
#if AV1_UPGRADE
    write_tile_info_max_tile(pcs_ptr, wb);

#if TILES
    if (pcs_ptr->av1_cm->tile_rows * pcs_ptr->av1_cm->tile_cols > 1) {

        // tile id used for cdf update
        aom_wb_write_literal(wb, 0, pcs_ptr->av1_cm->log2_tile_cols + pcs_ptr->av1_cm->log2_tile_rows);
        // Number of bytes in tile size - 1
        aom_wb_write_literal(wb, 3, 2);
    }
#endif

#else
    if (pcs_ptr->large_scale_tile) {
        printf("ERROR[AN]: large_scale_tile not supported yet\n");
        //const int32_t tile_width =
        //    ALIGN_POWER_OF_TWO(cm->tile_width, cm->seq_params.mib_size_log2) >>
        //    cm->seq_params.mib_size_log2;
        //const int32_t tile_height =
        //    ALIGN_POWER_OF_TWO(cm->tile_height, cm->seq_params.mib_size_log2) >>
        //    cm->seq_params.mib_size_log2;

        //assert(tile_width > 0);
        //assert(tile_height > 0);

        //// Write the tile sizes
        //if (cm->seq_params.sb_size == BLOCK_128X128) {
        //    assert(tile_width <= 32);
        //    assert(tile_height <= 32);
        //    aom_wb_write_literal(wb, tile_width - 1, 5);
        //    aom_wb_write_literal(wb, tile_height - 1, 5);
        //} else {
        //    assert(tile_width <= 64);
        //    assert(tile_height <= 64);
        //    aom_wb_write_literal(wb, tile_width - 1, 6);
        //    aom_wb_write_literal(wb, tile_height - 1, 6);
        //}
    }
    else {
        write_tile_info_max_tile(pcs_ptr, wb);
    }
#endif
    /**saved_wb = *wb;*/
#if !AV1_UPGRADE
    if (pcs_ptr->large_scale_tile) {

        printf("ERROR[AN]: large_scale_tile not supported yet\n");
#endif
        //if (cm->tile_rows * cm->tile_cols > 1) {
        //    // Note that the last item in the uncompressed header is the data
        //    // describing tile configuration.
        //    // Number of bytes in tile column size - 1
        //    aom_wb_write_literal(wb, 0, 2);
        //    // Number of bytes in tile size - 1
        //    aom_wb_write_literal(wb, 0, 2);
        //}
        //return;
#if !AV1_UPGRADE
    }
#endif
    //if (cm->tile_rows * cm->tile_cols > 1) {
    //    // Number of bytes in tile size - 1
    //    aom_wb_write_literal(wb, 3, 2);
    //}
}

static void write_frame_size(PictureParentControlSet_t *pcsPtr,
    int32_t frame_size_override,
    struct aom_write_bit_buffer *wb) {
    SequenceControlSet_t *scsPtr = (SequenceControlSet_t*)pcsPtr->sequence_control_set_wrapper_ptr->objectPtr;
    (void)(*pcsPtr);
    (void)frame_size_override;
    //const int32_t coded_width = cm->superres_upscaled_width - 1;
    //const int32_t coded_height = cm->superres_upscaled_height - 1;

    //if (frame_size_override) {
    //    const SequenceHeader *seq_params = &cm->seq_params;
    //    int32_t num_bits_width = seq_params->num_bits_width;
    //    int32_t num_bits_height = seq_params->num_bits_height;
    //    aom_wb_write_literal(wb, coded_width, num_bits_width);
    //    aom_wb_write_literal(wb, coded_height, num_bits_height);
    //}
    if (scsPtr->enable_superres) {
        printf("ERROR[AN]: enable_superres not supported yet\n");
        //write_superres_scale(cm, wb);
    }

    aom_wb_write_bit(wb, 0);
    //write_render_size(cm, wb);

}

static void WriteProfile(BITSTREAM_PROFILE profile,
    struct aom_write_bit_buffer *wb) {
    assert(profile >= PROFILE_0 && profile < MAX_PROFILES);
    aom_wb_write_literal(wb, profile, PROFILE_BITS);

}

static void write_bitdepth(SequenceControlSet_t *scsPtr/*Av1Common *const cm*/,
    struct aom_write_bit_buffer *wb) {
    // Profile 0/1: [0] for 8 bit, [1]  10-bit
    // Profile   2: [0] for 8 bit, [10] 10-bit, [11] - 12-bit
    aom_wb_write_bit(wb, scsPtr->static_config.encoder_bit_depth == EB_8BIT ? 0 : 1);
    if (scsPtr->static_config.profile == PROFILE_2 && scsPtr->static_config.encoder_bit_depth != EB_8BIT) {
        printf("ERROR[AN]: Profile 2 Not supported\n");
        aom_wb_write_bit(wb, scsPtr->static_config.encoder_bit_depth == EB_10BIT ? 0 : 1);
    }
}
#if AV1_UPGRADE
static void write_color_config(
    SequenceControlSet_t *scsPtr/*Av1Common *const cm*/, struct aom_write_bit_buffer *wb) {

    //write_bitdepth(cm, wb);
    write_bitdepth(scsPtr, wb);
    const int32_t is_monochrome = 0;// cm->seq_params.monochrome;
    // monochrome bit
    aom_wb_write_bit(wb, is_monochrome);
    //if (cm->profile != PROFILE_1)
    //    aom_wb_write_bit(wb, is_monochrome);
    //else
    //    assert(!is_monochrome);
    if (1/*cm->color_primaries == AOM_CICP_CP_UNSPECIFIED &&
         cm->transfer_characteristics == AOM_CICP_TC_UNSPECIFIED &&
         cm->matrix_coefficients == AOM_CICP_MC_UNSPECIFIED*/) {
        aom_wb_write_bit(wb, 0);  // No color description present
    }
    else {
        //aom_wb_write_bit(wb, 1);  // Color description present
        //aom_wb_write_literal(wb, cm->color_primaries, 8);
        //aom_wb_write_literal(wb, cm->transfer_characteristics, 8);
        //aom_wb_write_literal(wb, cm->matrix_coefficients, 8);
    }
    if (is_monochrome) {
        printf("ERROR[AN]: is_monochrome not supported yet\n");
        return;
    }
    if (0/*cm->color_primaries == AOM_CICP_CP_BT_709 &&
         cm->transfer_characteristics == AOM_CICP_TC_SRGB &&
         cm->matrix_coefficients ==
         AOM_CICP_MC_IDENTITY*/) {  // it would be better to remove this
         // dependency too
         //assert(cm->subsampling_x == 0 && cm->subsampling_y == 0);
         //assert(cm->profile == PROFILE_1 ||
         //    (cm->profile == PROFILE_2 && cm->bit_depth == AOM_BITS_12));
    }
    else {
        // 0: [16, 235] (i.e. xvYCC), 1: [0, 255]
        aom_wb_write_bit(wb, 0);
        //aom_wb_write_bit(wb, cm->color_range);

        //if (cm->profile == PROFILE_0) {
        //    // 420 only
        //    assert(cm->subsampling_x == 1 && cm->subsampling_y == 1);
        //}
        //else if (cm->profile == PROFILE_1) {
        //    // 444 only
        //    assert(cm->subsampling_x == 0 && cm->subsampling_y == 0);
        //}
        //else if (cm->profile == PROFILE_2) {
        //    if (cm->bit_depth == AOM_BITS_12) {
        //        // 420, 444 or 422
        //        aom_wb_write_bit(wb, cm->subsampling_x);
        //        if (cm->subsampling_x == 0) {
        //            assert(cm->subsampling_y == 0 &&
        //                "4:4:0 subsampling not allowed in AV1");
        //        }
        //        else {
        //            aom_wb_write_bit(wb, cm->subsampling_y);
        //        }
        //    }
        //    else {
        //        // 422 only
        //        assert(cm->subsampling_x == 1 && cm->subsampling_y == 0);
        //    }
        //}
        //if (cm->matrix_coefficients == AOM_CICP_MC_IDENTITY) {
        //    assert(cm->subsampling_x == 0 && cm->subsampling_y == 0);
        //}
        aom_wb_write_literal(wb, 0, 2);
        //if (cm->subsampling_x == 1 && cm->subsampling_y == 1) {
        //    aom_wb_write_literal(wb, cm->chroma_sample_position, 2);
        //}
    }
    aom_wb_write_bit(wb, 0);
    //aom_wb_write_bit(wb, cm->separate_uv_delta_q);
}

#else
static void WriteBitdepthColorspaceSampling(
    SequenceControlSet_t *scsPtr/*Av1Common *const cm*/, struct aom_write_bit_buffer *wb) {
    write_bitdepth(scsPtr, wb);

    const int32_t is_monochrome = 0;// cm->seq_params.monochrome;
    // monochrome bit
    aom_wb_write_bit(wb, is_monochrome);
    //if (cm->profile != PROFILE_1)
    //    aom_wb_write_bit(wb, is_monochrome);
    //else
    //    assert(!is_monochrome);


    if (1/*cm->color_primaries == AOM_CICP_CP_UNSPECIFIED &&
         cm->transfer_characteristics == AOM_CICP_TC_UNSPECIFIED &&
         cm->matrix_coefficients == AOM_CICP_MC_UNSPECIFIED*/) {
        aom_wb_write_bit(wb, 0);  // No color description present
    }
    else {
        //aom_wb_write_bit(wb, 1);  // Color description present
        //aom_wb_write_literal(wb, cm->color_primaries, 8);
        //aom_wb_write_literal(wb, cm->transfer_characteristics, 8);
        //aom_wb_write_literal(wb, cm->matrix_coefficients, 8);
    }

    if (is_monochrome) {
        printf("ERROR[AN]: is_monochrome not supported yet\n");
        return;
    }



    if (0/*cm->color_primaries == AOM_CICP_CP_BT_709 &&
         cm->transfer_characteristics == AOM_CICP_TC_SRGB &&
         cm->matrix_coefficients ==
         AOM_CICP_MC_IDENTITY*/) {  // it would be better to remove this
         // dependency too
         //assert(cm->subsampling_x == 0 && cm->subsampling_y == 0);
         //assert(cm->profile == PROFILE_1 ||
         //    (cm->profile == PROFILE_2 && cm->bit_depth == AOM_BITS_12));
    }
    else {
        // 0: [16, 235] (i.e. xvYCC), 1: [0, 255]
        aom_wb_write_bit(wb, 0);
        //aom_wb_write_bit(wb, cm->color_range);

        //if (cm->profile == PROFILE_0) {
        //    // 420 only
        //    assert(cm->subsampling_x == 1 && cm->subsampling_y == 1);
        //}
        //else if (cm->profile == PROFILE_1) {
        //    // 444 only
        //    assert(cm->subsampling_x == 0 && cm->subsampling_y == 0);
        //}
        //else if (cm->profile == PROFILE_2) {
        //    if (cm->bit_depth == AOM_BITS_12) {
        //        // 420, 444 or 422
        //        aom_wb_write_bit(wb, cm->subsampling_x);
        //        if (cm->subsampling_x == 0) {
        //            assert(cm->subsampling_y == 0 &&
        //                "4:4:0 subsampling not allowed in AV1");
        //        }
        //        else {
        //            aom_wb_write_bit(wb, cm->subsampling_y);
        //        }
        //    }
        //    else {
        //        // 422 only
        //        assert(cm->subsampling_x == 1 && cm->subsampling_y == 0);
        //    }
        //}
        aom_wb_write_literal(wb, 0, 2);
        //if (cm->subsampling_x == 1 && cm->subsampling_y == 1) {
        //    aom_wb_write_literal(wb, cm->chroma_sample_position, 2);
        //}
    }
    aom_wb_write_bit(wb, 0);
    //aom_wb_write_bit(wb, cm->separate_uv_delta_q);
}
#endif
void WriteSequenceHeader(SequenceControlSet_t *scsPtr/*AV1_COMP *cpi*/, struct aom_write_bit_buffer *wb) {
    //    Av1Common *const cm = &cpi->common;
    //    SequenceHeader *seq_params = &cm->seq_params;
    //



    int32_t max_frame_width = scsPtr->luma_width;
    /*        cpi->oxcf.forced_max_frame_width
    ? cpi->oxcf.forced_max_frame_width
    : cpi->oxcf.width;*/
    int32_t max_frame_height = scsPtr->luma_height;
    /*cpi->oxcf.forced_max_frame_height
    ? cpi->oxcf.forced_max_frame_height
    : cpi->oxcf.height;*/

    aom_wb_write_literal(wb, scsPtr->num_bits_width - 1, 4);
    aom_wb_write_literal(wb, scsPtr->num_bits_height - 1, 4);
    aom_wb_write_literal(wb, max_frame_width - 1, scsPtr->num_bits_width);
    aom_wb_write_literal(wb, max_frame_height - 1, scsPtr->num_bits_height);
    //
    /* Placeholder for actually writing to the bitstream */

    if (!scsPtr->reduced_still_picture_hdr) {

        //scsPtr->frame_id_numbers_present_flag = 0;
        //    cm->large_scale_tile ? 0 : cm->error_resilient_mode;

        aom_wb_write_bit(wb, scsPtr->frame_id_numbers_present_flag);
        if (scsPtr->frame_id_numbers_present_flag) {
            // We must always have delta_frame_id_length < frame_id_length,
            // in order for a frame to be referenced with a unique delta.
            // Avoid wasting bits by using a coding that enforces this restriction.
            aom_wb_write_literal(wb, scsPtr->delta_frame_id_length - 2, 4);
            aom_wb_write_literal(
                wb, scsPtr->frame_id_length - scsPtr->delta_frame_id_length - 1,
                3);
        }

    }

    aom_wb_write_bit(wb, scsPtr->sb_size == BLOCK_128X128 ? 1 : 0);
    //    write_sb_size(seq_params, wb);
    aom_wb_write_bit(wb, scsPtr->enable_filter_intra);
#if !INTRA_10BIT_SUPPORT
    if (scsPtr->static_config.encoder_bit_depth == 8)
#endif

#if  DIS_EDGE_FIL
        scsPtr->enable_intra_edge_filter = 0;
#else
        scsPtr->enable_intra_edge_filter = 1;
#endif
    aom_wb_write_bit(wb, scsPtr->enable_intra_edge_filter);

    if (!scsPtr->reduced_still_picture_hdr) {

        aom_wb_write_bit(wb, scsPtr->enable_interintra_compound);
        aom_wb_write_bit(wb, scsPtr->enable_masked_compound);
        aom_wb_write_bit(wb, scsPtr->static_config.enable_warped_motion);
        aom_wb_write_bit(wb, scsPtr->enable_dual_filter);

        aom_wb_write_bit(wb, scsPtr->enable_order_hint);

        if (scsPtr->enable_order_hint) {
            aom_wb_write_bit(wb, scsPtr->enable_jnt_comp);
            aom_wb_write_bit(wb, scsPtr->enable_ref_frame_mvs);
        }

        if (scsPtr->force_screen_content_tools == 2) {
            aom_wb_write_bit(wb, 1);
        }
        else {
            aom_wb_write_bit(wb, 0);
            aom_wb_write_bit(wb, scsPtr->force_screen_content_tools);
        }
        //
        if (scsPtr->force_screen_content_tools > 0) {
            if (scsPtr->force_integer_mv == 2) {
                aom_wb_write_bit(wb, 1);
            }
            else {
                aom_wb_write_bit(wb, 0);
                aom_wb_write_bit(wb, scsPtr->force_integer_mv);
            }
        }
        else {
            assert(scsPtr->force_integer_mv == 2);
        }

        if (scsPtr->enable_order_hint)
            aom_wb_write_literal(wb, scsPtr->order_hint_bits_minus1, 3);
    }

    aom_wb_write_bit(wb, scsPtr->enable_superres);
    aom_wb_write_bit(wb, scsPtr->enable_cdef);
    aom_wb_write_bit(wb, scsPtr->enable_restoration);
}

// Recenters a non-negative literal v around a reference r
static uint16_t recenter_nonneg(uint16_t r, uint16_t v) {
    if (v > (r << 1))
        return v;
    else if (v >= r)
        return ((v - r) << 1);
    else
        return ((r - v) << 1) - 1;
}

// Recenters a non-negative literal v in [0, n-1] around a
// reference r also in [0, n-1]
static uint16_t recenter_finite_nonneg(uint16_t n, uint16_t r, uint16_t v) {
    if ((r << 1) <= n) {
        return recenter_nonneg(r, v);
    }
    else {
        return recenter_nonneg(n - 1 - r, n - 1 - v);
    }
}

int32_t aom_count_primitive_symmetric(int16_t v, uint32_t abs_bits) {
    return (v == 0 ? 1 : abs_bits + 2);
}

// Encodes a value v in [0, n-1] quasi-uniformly
void aom_write_primitive_quniform(aom_writer *w, uint16_t n, uint16_t v) {
    if (n <= 1) return;
    const int32_t l = get_msb(n - 1) + 1;
    const int32_t m = (1 << l) - n;
    if (v < m) {
        aom_write_literal(w, v, l - 1);
    }
    else {
        aom_write_literal(w, m + ((v - m) >> 1), l - 1);
        aom_write_bit(w, (v - m) & 1);
    }
}

static void aom_wb_write_primitive_quniform(struct aom_write_bit_buffer *wb,
    uint16_t n, uint16_t v) {
    if (n <= 1) return;
    const int32_t l = get_msb(n - 1) + 1;
    const int32_t m = (1 << l) - n;
    if (v < m) {
        aom_wb_write_literal(wb, v, l - 1);
    }
    else {
        aom_wb_write_literal(wb, m + ((v - m) >> 1), l - 1);
        aom_wb_write_bit(wb, (v - m) & 1);
    }
}

int32_t aom_count_primitive_quniform(uint16_t n, uint16_t v) {
    if (n <= 1) return 0;
    const int32_t l = get_msb(n - 1) + 1;
    const int32_t m = (1 << l) - n;
    return v < m ? l - 1 : l;
}

// Finite subexponential code that codes a symbol v in [0, n-1] with parameter k
void aom_write_primitive_subexpfin(aom_writer *w, uint16_t n, uint16_t k,
    uint16_t v) {
    int32_t i = 0;
    int32_t mk = 0;
    while (1) {
        int32_t b = (i ? k + i - 1 : k);
        int32_t a = (1 << b);
        if (n <= mk + 3 * a) {
            aom_write_primitive_quniform(w, (uint16_t)(n - mk), (uint16_t)(v - mk));
            break;
        }
        else {
            int32_t t = (v >= mk + a);
            aom_write_bit(w, t);
            if (t) {
                i = i + 1;
                mk += a;
            }
            else {
                aom_write_literal(w, v - mk, b);
                break;
            }
        }
    }
}

static void aom_wb_write_primitive_subexpfin(struct aom_write_bit_buffer *wb,
    uint16_t n, uint16_t k,
    uint16_t v) {
    int32_t i = 0;
    int32_t mk = 0;
    while (1) {
        int32_t b = (i ? k + i - 1 : k);
        int32_t a = (1 << b);
        if (n <= mk + 3 * a) {
            aom_wb_write_primitive_quniform(wb, (uint16_t)(n - mk), (uint16_t)(v - mk));
            break;
        }
        else {
            int32_t t = (v >= mk + a);
            aom_wb_write_bit(wb, t);
            if (t) {
                i = i + 1;
                mk += a;
            }
            else {
                aom_wb_write_literal(wb, v - mk, b);
                break;
            }
        }
    }
}

int32_t aom_count_primitive_subexpfin(uint16_t n, uint16_t k, uint16_t v) {
    int32_t count = 0;
    int32_t i = 0;
    int32_t mk = 0;
    while (1) {
        int32_t b = (i ? k + i - 1 : k);
        int32_t a = (1 << b);
        if (n <= mk + 3 * a) {
            count += aom_count_primitive_quniform((uint16_t)(n - mk), (uint16_t)(v - mk));
            break;
        }
        else {
            int32_t t = (v >= mk + a);
            count++;
            if (t) {
                i = i + 1;
                mk += a;
            }
            else {
                count += b;
                break;
            }
        }
    }
    return count;
}
// Finite subexponential code that codes a symbol v in[0, n - 1] with parameter k
// based on a reference ref also in [0, n-1].
// Recenters symbol around r first and then uses a finite subexponential code.
void aom_write_primitive_refsubexpfin(aom_writer *w, uint16_t n, uint16_t k,
    uint16_t ref, uint16_t v) {
    aom_write_primitive_subexpfin(w, n, k, recenter_finite_nonneg(n, ref, v));
}

static void aom_wb_write_primitive_refsubexpfin(struct aom_write_bit_buffer *wb,
    uint16_t n, uint16_t k,
    uint16_t ref, uint16_t v) {
    aom_wb_write_primitive_subexpfin(wb, n, k, recenter_finite_nonneg(n, ref, v));
}

void aom_wb_write_signed_primitive_refsubexpfin(struct aom_write_bit_buffer *wb,
    uint16_t n, uint16_t k,
    int16_t ref, int16_t v) {
    ref += n - 1;
    v += n - 1;
    const uint16_t scaled_n = (n << 1) - 1;
    aom_wb_write_primitive_refsubexpfin(wb, scaled_n, k, ref, v);
}

int32_t aom_count_primitive_refsubexpfin(uint16_t n, uint16_t k, uint16_t ref,
    uint16_t v) {
    return aom_count_primitive_subexpfin(n, k, recenter_finite_nonneg(n, ref, v));
}

static void write_global_motion_params(const EbWarpedMotionParams *params,
    const EbWarpedMotionParams *ref_params,
    struct aom_write_bit_buffer *wb,
    int32_t allow_hp) {
    const TransformationType type = params->wmtype;
    assert(type == TRANSLATION || type == IDENTITY);
    aom_wb_write_bit(wb, type != IDENTITY);
    if (type != IDENTITY) {
#if GLOBAL_TRANS_TYPES > 4
        aom_wb_write_literal(wb, type - 1, GLOBAL_TYPE_BITS);
#else
        aom_wb_write_bit(wb, type == ROTZOOM);
        if (type != ROTZOOM) aom_wb_write_bit(wb, type == TRANSLATION);
#endif  // GLOBAL_TRANS_TYPES > 4
    }

    if (type >= ROTZOOM) {

        int16_t ref2 = (int16_t)((ref_params->wmmat[2] >> GM_ALPHA_PREC_DIFF) - (1 << GM_ALPHA_PREC_BITS));
        int16_t v2 = (int16_t)((params->wmmat[2] >> GM_ALPHA_PREC_DIFF) - (1 << GM_ALPHA_PREC_BITS));

        int16_t ref3 = (int16_t)(ref_params->wmmat[3] >> GM_ALPHA_PREC_DIFF);
        int16_t v3 = (int16_t)(params->wmmat[3] >> GM_ALPHA_PREC_DIFF);

        aom_wb_write_signed_primitive_refsubexpfin(
            wb, GM_ALPHA_MAX + 1, SUBEXPFIN_K,
            ref2/*(ref_params->wmmat[2] >> GM_ALPHA_PREC_DIFF) - (1 << GM_ALPHA_PREC_BITS)*/,
            v2/*(int16_t)((params->wmmat[2] >> GM_ALPHA_PREC_DIFF) - (1 << GM_ALPHA_PREC_BITS))*/);
        aom_wb_write_signed_primitive_refsubexpfin(
            wb, GM_ALPHA_MAX + 1, SUBEXPFIN_K,
            ref3/*(ref_params->wmmat[3] >> GM_ALPHA_PREC_DIFF)*/,
            v3/*(int16_t)(params->wmmat[3] >> GM_ALPHA_PREC_DIFF)*/);
    }

    if (type >= AFFINE) {
        int16_t ref4 = (int16_t)(ref_params->wmmat[4] >> GM_ALPHA_PREC_DIFF);
        int16_t v4 = (int16_t)(params->wmmat[4] >> GM_ALPHA_PREC_DIFF);

        int16_t ref5 = (int16_t)((ref_params->wmmat[5] >> GM_ALPHA_PREC_DIFF) - (1 << GM_ALPHA_PREC_BITS));
        int16_t v5 = (int16_t)((params->wmmat[5] >> GM_ALPHA_PREC_DIFF) - (1 << GM_ALPHA_PREC_BITS));

        aom_wb_write_signed_primitive_refsubexpfin(
            wb, GM_ALPHA_MAX + 1, SUBEXPFIN_K,
            ref4/*(ref_params->wmmat[4] >> GM_ALPHA_PREC_DIFF)*/,
            v4/*(int16_t)(params->wmmat[4] >> GM_ALPHA_PREC_DIFF)*/);
        aom_wb_write_signed_primitive_refsubexpfin(
            wb, GM_ALPHA_MAX + 1, SUBEXPFIN_K,
            ref5/*(ref_params->wmmat[5] >> GM_ALPHA_PREC_DIFF) -    (1 << GM_ALPHA_PREC_BITS)*/,
            v5/*(int16_t)(params->wmmat[5] >> GM_ALPHA_PREC_DIFF) - (1 << GM_ALPHA_PREC_BITS)*/);
    }

    if (type >= TRANSLATION) {

        const int32_t trans_bits = (type == TRANSLATION)
            ? GM_ABS_TRANS_ONLY_BITS - !allow_hp
            : GM_ABS_TRANS_BITS;
        const int32_t trans_prec_diff = (type == TRANSLATION)
            ? GM_TRANS_ONLY_PREC_DIFF + !allow_hp
            : GM_TRANS_PREC_DIFF;
        aom_wb_write_signed_primitive_refsubexpfin(
            wb, (1 << trans_bits) + 1, SUBEXPFIN_K,
            (int16_t)(ref_params->wmmat[0] >> trans_prec_diff),
            (int16_t)(params->wmmat[0] >> trans_prec_diff));
        aom_wb_write_signed_primitive_refsubexpfin(
            wb, (1 << trans_bits) + 1, SUBEXPFIN_K,
            (int16_t)(ref_params->wmmat[1] >> trans_prec_diff),
            (int16_t)(params->wmmat[1] >> trans_prec_diff));
    }
}
static void WriteGlobalMotion(
    PictureParentControlSet_t *pcsPtr,
    struct aom_write_bit_buffer *wb)

{
    int32_t frame;
    for (frame = LAST_FRAME; frame <= ALTREF_FRAME; ++frame) {
        const EbWarpedMotionParams *ref_params = &default_warp_params;
        //pcsPtr->prev_frame ? &pcsPtr->prev_frame->global_motion[frame] : &default_warp_params;

        write_global_motion_params(&pcsPtr->global_motion[frame], ref_params, wb,
            pcsPtr->allow_high_precision_mv);
        // TODO(sarahparker, debargha): The logic in the commented out code below
        // does not work currently and causes mismatches when resize is on.
        // Fix it before turning the optimization back on.
        /*
        Yv12BufferConfig *ref_buf = get_ref_frame_buffer(cpi, frame);
        if (cpi->source->y_crop_width == ref_buf->y_crop_width &&
        cpi->source->y_crop_height == ref_buf->y_crop_height) {
        write_global_motion_params(&cm->global_motion[frame],
        &cm->prev_frame->global_motion[frame], wb,
        cm->allow_high_precision_mv);
        } else {
        assert(cm->global_motion[frame].wmtype == IDENTITY &&
        "Invalid warp type for frames of different resolutions");
        }
        */
        /*
        printf("Frame %d/%d: Enc Ref %d: %d %d %d %d\n",
        cm->current_video_frame, cm->show_frame, frame,
        cm->global_motion[frame].wmmat[0],
        cm->global_motion[frame].wmmat[1], cm->global_motion[frame].wmmat[2],
        cm->global_motion[frame].wmmat[3]);
        */
    }
}

static void write_film_grain_params(PictureParentControlSet_t *pcsPtr,
    struct aom_write_bit_buffer *wb) {

    aom_film_grain_t *pars = &pcsPtr->film_grain_params;

    aom_wb_write_bit(wb, pars->apply_grain);
    if (!pars->apply_grain) return;

    aom_wb_write_literal(wb, pars->random_seed, 16);

    if (pcsPtr->av1FrameType == INTER_FRAME) {
        EbReferenceObject_t* refObj0 = (EbReferenceObject_t*)pcsPtr->childPcs->ref_pic_ptr_array[REF_LIST_0]->objectPtr;
        int32_t ref_idx = 0;
        pars->update_parameters = 1;
        if (film_grain_params_equal(&refObj0->film_grain_params, pars)) {
            pars->update_parameters = 0;
            ref_idx = get_ref_frame_map_idx(pcsPtr, LAST_FRAME);
        }
        else if (pcsPtr->childPcs->slice_type == B_SLICE)
        {
            EbReferenceObject_t* refObj1 = (EbReferenceObject_t*)pcsPtr->childPcs->ref_pic_ptr_array[REF_LIST_1]->objectPtr;
            if (film_grain_params_equal(&refObj1->film_grain_params, pars)) {
                pars->update_parameters = 0;
                ref_idx = get_ref_frame_map_idx(pcsPtr, ALTREF_FRAME);  //todo: will it always be ALF_REF in L1?
            }
        }
        aom_wb_write_bit(wb, pars->update_parameters);
        if (!pars->update_parameters) {
            aom_wb_write_literal(wb, ref_idx, 3);
            return;
        }
    }
    else
        pars->update_parameters = 1;


    // Scaling functions parameters
    aom_wb_write_literal(wb, pars->num_y_points, 4);  // max 14
    for (int32_t i = 0; i < pars->num_y_points; i++) {
        aom_wb_write_literal(wb, pars->scaling_points_y[i][0], 8);
        aom_wb_write_literal(wb, pars->scaling_points_y[i][1], 8);
    }

    if (!pcsPtr->sequence_control_set_ptr->monochrome)
        aom_wb_write_bit(wb, pars->chroma_scaling_from_luma);
    else
        pars->chroma_scaling_from_luma = 0;  // for monochrome override to 0

    if (pcsPtr->sequence_control_set_ptr->monochrome || pars->chroma_scaling_from_luma ||
        // todo: add corresponding check when subsampling variables are present
        ((pcsPtr->sequence_control_set_ptr->subsampling_x == 1) &&
        (pcsPtr->sequence_control_set_ptr->subsampling_y == 1) &&
            (pars->num_y_points == 0))) {
        pars->num_cb_points = 0;
        pars->num_cr_points = 0;
    }
    else {
        aom_wb_write_literal(wb, pars->num_cb_points, 4);  // max 10
        for (int32_t i = 0; i < pars->num_cb_points; i++) {
            aom_wb_write_literal(wb, pars->scaling_points_cb[i][0], 8);
            aom_wb_write_literal(wb, pars->scaling_points_cb[i][1], 8);
        }

        aom_wb_write_literal(wb, pars->num_cr_points, 4);  // max 10
        for (int32_t i = 0; i < pars->num_cr_points; i++) {
            aom_wb_write_literal(wb, pars->scaling_points_cr[i][0], 8);
            aom_wb_write_literal(wb, pars->scaling_points_cr[i][1], 8);
        }
    }

    aom_wb_write_literal(wb, pars->scaling_shift - 8, 2);  // 8 + value

    // AR coefficients
    // Only sent if the corresponsing scaling function has
    // more than 0 points

    aom_wb_write_literal(wb, pars->ar_coeff_lag, 2);

    int32_t num_pos_luma = 2 * pars->ar_coeff_lag * (pars->ar_coeff_lag + 1);
    int32_t num_pos_chroma = num_pos_luma;
    if (pars->num_y_points > 0) ++num_pos_chroma;

    if (pars->num_y_points)
        for (int32_t i = 0; i < num_pos_luma; i++)
            aom_wb_write_literal(wb, pars->ar_coeffs_y[i] + 128, 8);

    if (pars->num_cb_points || pars->chroma_scaling_from_luma)
        for (int32_t i = 0; i < num_pos_chroma; i++)
            aom_wb_write_literal(wb, pars->ar_coeffs_cb[i] + 128, 8);

    if (pars->num_cr_points || pars->chroma_scaling_from_luma)
        for (int32_t i = 0; i < num_pos_chroma; i++)
            aom_wb_write_literal(wb, pars->ar_coeffs_cr[i] + 128, 8);

    aom_wb_write_literal(wb, pars->ar_coeff_shift - 6, 2);  // 8 + value

    aom_wb_write_literal(wb, pars->grain_scale_shift, 2);

    if (pars->num_cb_points) {
        aom_wb_write_literal(wb, pars->cb_mult, 8);
        aom_wb_write_literal(wb, pars->cb_luma_mult, 8);
        aom_wb_write_literal(wb, pars->cb_offset, 9);
    }

    if (pars->num_cr_points) {
        aom_wb_write_literal(wb, pars->cr_mult, 8);
        aom_wb_write_literal(wb, pars->cr_luma_mult, 8);
        aom_wb_write_literal(wb, pars->cr_offset, 9);
    }

    aom_wb_write_bit(wb, pars->overlap_flag);

    aom_wb_write_bit(wb, pars->clip_to_restricted_range);
}

// New function based on HLS R18
static void WriteUncompressedHeaderObu(SequenceControlSet_t *scsPtr/*AV1_COMP *cpi*/,
    PictureParentControlSet_t *pcsPtr,
    //struct aom_write_bit_buffer *saved_wb,
    struct aom_write_bit_buffer *wb,
    uint8_t showExisting) {
    // Av1Common *const cm = &cpi->common;
    // MacroBlockD *const xd = &cpi->td.mb.e_mbd;

    // NOTE: By default all coded frames to be used as a reference
    pcsPtr->is_reference_frame = 1;

    if (!scsPtr->reduced_still_picture_hdr) {

        if (showExisting) {
            //printf("ERROR[AN]: show_existing_frame not supported yet\n");
            //RefCntBuffer *const frame_bufs = cm->buffer_pool->frame_bufs;
            //const int32_t frame_to_show = cm->ref_frame_map[cpi->showExistingLoc];

            //if (frame_to_show < 0 || frame_bufs[frame_to_show].ref_count < 1) {
            //    aom_internal_error(&cm->error, AOM_CODEC_UNSUP_BITSTREAM,
            //        "Buffer %d does not contain a reconstructed frame",
            //        frame_to_show);
            //}
            //ref_cnt_fb(frame_bufs, &cm->new_fb_idx, frame_to_show);

            aom_wb_write_bit(wb, 1);  // show_existing_frame
            aom_wb_write_literal(wb, pcsPtr->showExistingLoc, 3);

            if (scsPtr->frame_id_numbers_present_flag) {
                printf("ERROR[AN]: frame_id_numbers_present_flag not supported yet\n");
                /*int32_t frame_id_len = cm->seq_params.frame_id_length;
                int32_t display_frame_id = cm->ref_frame_id[cpi->showExistingLoc];
                aom_wb_write_literal(wb, display_frame_id, frame_id_len);*/
            }

            //        if (cm->reset_decoder_state &&
            //            frame_bufs[frame_to_show].av1FrameType != KEY_FRAME) {
            //            aom_internal_error(
            //                &cm->error, AOM_CODEC_UNSUP_BITSTREAM,
            //                "show_existing_frame to reset state on KEY_FRAME only");
            //        }


            return;
        }
        else {
            aom_wb_write_bit(wb, 0);  // show_existing_frame
        }


        //pcsPtr->av1FrameType = pcsPtr->intra_only ? INTRA_ONLY_FRAME : pcsPtr->av1FrameType;

        aom_wb_write_literal(wb, pcsPtr->av1FrameType, 2);

        // if (pcsPtr->intra_only) pcsPtr->av1FrameType = INTRA_ONLY_FRAME;

        aom_wb_write_bit(wb, pcsPtr->showFrame);

        if (!pcsPtr->showFrame) {
            aom_wb_write_bit(wb, pcsPtr->showable_frame);
        }
        if (pcsPtr->av1FrameType == S_FRAME) {
            assert(pcsPtr->error_resilient_mode);
        }

        else if (!(pcsPtr->av1FrameType == KEY_FRAME && pcsPtr->showFrame)) {

            aom_wb_write_bit(wb, pcsPtr->error_resilient_mode);
        }

    }

    aom_wb_write_bit(wb, pcsPtr->disable_cdf_update);

    if (scsPtr->force_screen_content_tools == 2) {
        aom_wb_write_bit(wb, pcsPtr->allow_screen_content_tools);
    }
    else {
        assert(pcsPtr->allow_screen_content_tools ==
            scsPtr->force_screen_content_tools);
    }

    if (pcsPtr->allow_screen_content_tools) {
        if (scsPtr->force_integer_mv == 2) {
            aom_wb_write_bit(wb, pcsPtr->cur_frame_force_integer_mv);
        }
        else {
            assert(pcsPtr->cur_frame_force_integer_mv == scsPtr->force_integer_mv);
        }
    }
    else {
        assert(pcsPtr->cur_frame_force_integer_mv == 0);
    }

    if (!scsPtr->reduced_still_picture_hdr) {

        if (scsPtr->frame_id_numbers_present_flag) {
            int32_t frame_id_len = scsPtr->frame_id_length;
            aom_wb_write_literal(wb, pcsPtr->current_frame_id, frame_id_len);
        }

        //if (cm->width > cm->seq_params.max_frame_width ||
        //    cm->height > cm->seq_params.max_frame_height) {
        //    aom_internal_error(&cm->error, AOM_CODEC_UNSUP_BITSTREAM,
        //        "Frame dimensions are larger than the maximum values");
        //}

        int32_t frame_size_override_flag = 0;
        /*        (pcsPtr->av1FrameType == S_FRAME) ? 1
        : (cm->width != cm->seq_params.max_frame_width ||
        cm->height != cm->seq_params.max_frame_height);*/
        if (pcsPtr->av1FrameType != S_FRAME) aom_wb_write_bit(wb, frame_size_override_flag);

        if (scsPtr->enable_order_hint)
            aom_wb_write_literal(wb, (int32_t)pcsPtr->frame_offset,
                scsPtr->order_hint_bits_minus1 + 1);


        if (!pcsPtr->error_resilient_mode && !frame_is_intra_only(pcsPtr)) {
            aom_wb_write_literal(wb, pcsPtr->primary_ref_frame, PRIMARY_REF_BITS);
        }

    }
    int32_t frame_size_override_flag = 0;
    if (pcsPtr->av1FrameType == KEY_FRAME) {
        if (!pcsPtr->showFrame) {
            aom_wb_write_literal(wb, pcsPtr->av1RefSignal.refreshFrameMask, REF_FRAMES);
        }
    }
    else {

        if (pcsPtr->av1FrameType == INTRA_ONLY_FRAME) {
            // pcsPtr->refresh_frame_mask = get_refresh_mask(cpi);
            int32_t updated_fb = -1;
            for (int32_t i = 0; i < REF_FRAMES; i++) {
                // If more than one frame is refreshed, it doesn't matter which one
                // we pick, so pick the first.
                if (pcsPtr->av1RefSignal.refreshFrameMask & (1 << i)) {
                    updated_fb = i;
                    break;
                }
            }
            assert(updated_fb >= 0);
            pcsPtr->fb_of_context_type[pcsPtr->frame_context_idx] = updated_fb;

            aom_wb_write_literal(wb, pcsPtr->av1RefSignal.refreshFrameMask, REF_FRAMES);

        }
        else if (pcsPtr->av1FrameType == INTER_FRAME || frame_is_sframe(pcsPtr)) {

            //pcsPtr->refresh_frame_mask = get_refresh_mask(cpi);
            if (pcsPtr->av1FrameType == INTER_FRAME) {
                aom_wb_write_literal(wb, pcsPtr->av1RefSignal.refreshFrameMask, REF_FRAMES);
            }
            else {
                assert(frame_is_sframe(pcsPtr) && pcsPtr->av1RefSignal.refreshFrameMask == 0xFF);
            }

            int32_t updated_fb = -1;
            for (int32_t i = 0; i < REF_FRAMES; i++) {
                // If more than one frame is refreshed, it doesn't matter which one
                // we pick, so pick the first.
                if (pcsPtr->av1RefSignal.refreshFrameMask & (1 << i)) {
                    updated_fb = i;
                    break;
                }
            }
            // large scale tile sometimes won't refresh any fbs
            if (updated_fb >= 0) {
                pcsPtr->fb_of_context_type[pcsPtr->frame_context_idx] = updated_fb;
            }

            if (!pcsPtr->av1RefSignal.refreshFrameMask) {
                // NOTE: "cpi->refresh_frame_mask == 0" indicates that the coded frame
                //       will not be used as a reference
                pcsPtr->is_reference_frame = 0;
            }

        }
    }

    if (pcsPtr->av1FrameType == KEY_FRAME) {
        write_frame_size(pcsPtr, frame_size_override_flag, wb);
        //assert(av1_superres_unscaled(cm) ||
        //    !(cm->allow_intrabc && NO_FILTER_FOR_IBC));
        if (pcsPtr->allow_screen_content_tools &&
            0 /*(av1_superres_unscaled(cm) || !NO_FILTER_FOR_IBC)*/)
            aom_wb_write_bit(wb, pcsPtr->allow_intrabc);
        // all eight fbs are refreshed, pick one that will live long enough
        pcsPtr->fb_of_context_type[REGULAR_FRAME] = 0;
    }
    else {
        if (pcsPtr->av1FrameType == INTRA_ONLY_FRAME) {
            write_frame_size(pcsPtr, frame_size_override_flag, wb);
            if (pcsPtr->allow_screen_content_tools &&
                0 /*(av1_superres_unscaled(cm) || !NO_FILTER_FOR_IBC)*/)
                aom_wb_write_bit(wb, pcsPtr->allow_intrabc);
        }
        else if (pcsPtr->av1FrameType == INTER_FRAME || frame_is_sframe(pcsPtr)) {
            MvReferenceFrame ref_frame;

            assert(pcsPtr->frame_refs_short_signaling == 0);
            // NOTE: Error resilient mode turns off frame_refs_short_signaling
            //       automatically.
            if (scsPtr->enable_order_hint)
                aom_wb_write_bit(wb, pcsPtr->frame_refs_short_signaling);

            if (pcsPtr->frame_refs_short_signaling) {
                aom_wb_write_literal(wb, get_ref_frame_map_idx(pcsPtr, LAST_FRAME),
                    REF_FRAMES_LOG2);
                aom_wb_write_literal(wb, get_ref_frame_map_idx(pcsPtr, GOLDEN_FRAME),
                    REF_FRAMES_LOG2);
            }
            for (ref_frame = LAST_FRAME; ref_frame <= ALTREF_FRAME; ++ref_frame) {
                assert(get_ref_frame_map_idx(pcsPtr, ref_frame) != INVALID_IDX);
                if (!pcsPtr->frame_refs_short_signaling)
                    aom_wb_write_literal(wb, get_ref_frame_map_idx(pcsPtr, ref_frame),
                        REF_FRAMES_LOG2);

                if (scsPtr->frame_id_numbers_present_flag) {
                    printf("ERROR[AN]: frame_id_numbers_present_flag not supported yet\n");
                    //int32_t i = get_ref_frame_map_idx(cpi, ref_frame);
                    //int32_t frame_id_len = cm->seq_params.frame_id_length;
                    //int32_t diff_len = cm->seq_params.delta_frame_id_length;
                    //int32_t delta_frame_id_minus1 =
                    //    ((cm->current_frame_id - cm->ref_frame_id[i] +
                    //    (1 << frame_id_len)) %
                    //    (1 << frame_id_len)) -
                    //    1;
                    //if (delta_frame_id_minus1 < 0 ||
                    //    delta_frame_id_minus1 >= (1 << diff_len))
                    //    cm->invalid_delta_frame_id_minus1 = 1;
                    //aom_wb_write_literal(wb, delta_frame_id_minus1, diff_len);
                }
            }

            if (!pcsPtr->error_resilient_mode && frame_size_override_flag) {
                printf("ERROR[AN]: frame_size_override_flag not supported yet\n");
                //write_frame_size_with_refs(pcsPtr, wb);
            }
            else {
                write_frame_size(pcsPtr, frame_size_override_flag, wb);
            }

            if (pcsPtr->cur_frame_force_integer_mv) {
                pcsPtr->allow_high_precision_mv = 0;
            }
            else {
                aom_wb_write_bit(wb, pcsPtr->allow_high_precision_mv);
            }

#define LOG_SWITCHABLE_FILTERS  2

            // OPT fix_interp_filter(cm, cpi->td.counts);
            // write_frame_interp_filter(cm->interp_filter, wb);
            aom_wb_write_bit(wb, pcsPtr->av1_cm->interp_filter == SWITCHABLE);
            if (pcsPtr->av1_cm->interp_filter != SWITCHABLE)
                aom_wb_write_literal(wb, pcsPtr->av1_cm->interp_filter, LOG_SWITCHABLE_FILTERS);

            aom_wb_write_bit(wb, pcsPtr->switchable_motion_mode);
            if (frame_might_allow_ref_frame_mvs(pcsPtr, scsPtr)) {

                aom_wb_write_bit(wb, pcsPtr->allow_ref_frame_mvs);
            }
        }
    }

    //if (scsPtr->frame_id_numbers_present_flag)
    //    pcsPtr->refresh_mask = get_refresh_mask(pcsPtr);
#if AV1_UPGRADE
    const int32_t might_bwd_adapt =
        !(scsPtr->reduced_still_picture_hdr) && !(pcsPtr->disable_cdf_update);
    if (pcsPtr->large_scale_tile)
        pcsPtr->refresh_frame_context = REFRESH_FRAME_CONTEXT_DISABLED;
#else
    const int32_t might_bwd_adapt =
        !(scsPtr->reduced_still_picture_hdr) && !(pcsPtr->large_scale_tile) && !(pcsPtr->disable_cdf_update);
#endif
    if (might_bwd_adapt) {
        aom_wb_write_bit(
            wb, pcsPtr->refresh_frame_context == REFRESH_FRAME_CONTEXT_DISABLED);
    }

    write_tile_info(pcsPtr, /*saved_wb,*/ wb);

    encode_quantization(pcsPtr, wb);

    aom_wb_write_bit(wb, 0);
    //encode_segmentation(cm, xd, wb);
        //if (pcsPtr->delta_q_present_flag)
           // assert(delta_q_allowed == 1 && pcsPtr->base_qindex > 0);

    if (pcsPtr->base_qindex > 0) {

        aom_wb_write_bit(wb, pcsPtr->delta_q_present_flag);
        if (pcsPtr->delta_q_present_flag) {

#if ADD_DELTA_QP_SUPPORT //PART 0
            aom_wb_write_literal(wb, OD_ILOG_NZ(pcsPtr->delta_q_res) - 1, 2);
            pcsPtr->prev_qindex = pcsPtr->base_qindex;
            if (pcsPtr->allow_intrabc)
                assert(pcsPtr->delta_lf_present_flag == 0);
            else
                aom_wb_write_bit(wb, pcsPtr->delta_lf_present_flag);
            if (pcsPtr->delta_lf_present_flag) {
                aom_wb_write_literal(wb, OD_ILOG_NZ(pcsPtr->delta_lf_res) - 1, 2);
                pcsPtr->prev_delta_lf_from_base = 0;
                aom_wb_write_bit(wb, pcsPtr->delta_lf_multi);
                const int32_t frame_lf_count =
                    pcsPtr->monochrome == 0 ? FRAME_LF_COUNT : FRAME_LF_COUNT - 2;
                for (int32_t lf_id = 0; lf_id < frame_lf_count; ++lf_id)
                    pcsPtr->prev_delta_lf[lf_id] = 0;
            }
#else
            printf("ERROR[AN]: delta_q_present_flag not supported yet\n");
#endif
            //                aom_wb_write_literal(wb, OD_ILOG_NZ(cm->delta_q_res) - 1, 2);
            //                xd->prev_qindex = cm->base_qindex;
            //                if (cm->allow_intrabc)
            //                    assert(cm->delta_lf_present_flag == 0);
            //                else
            //                    aom_wb_write_bit(wb, cm->delta_lf_present_flag);
            //                if (cm->delta_lf_present_flag) {
            //                    aom_wb_write_literal(wb, OD_ILOG_NZ(cm->delta_lf_res) - 1, 2);
            //                    xd->prev_delta_lf_from_base = 0;
            //                    aom_wb_write_bit(wb, cm->delta_lf_multi);
            //                    const int32_t frame_lf_count =
            //                        av1_num_planes(cm) > 1 ? FRAME_LF_COUNT : FRAME_LF_COUNT - 2;
            //                    for (int32_t lf_id = 0; lf_id < frame_lf_count; ++lf_id)
            //                        xd->prev_delta_lf[lf_id] = 0;
            //                }
        }
    }
    if (pcsPtr->all_lossless) {

        printf("ERROR[AN]: all_lossless\n");
        //assert(av1_superres_unscaled(pcsPtr));
    }
    else {
        if (!pcsPtr->coded_lossless) {
            encode_loopfilter(pcsPtr, wb);
            if (scsPtr->enable_cdef) {
                encode_cdef(pcsPtr, wb);
            }
        }
        if (scsPtr->enable_restoration) {
            encode_restoration_mode(pcsPtr, wb);
        }
    }

    aom_wb_write_bit(wb, pcsPtr->tx_mode == TX_MODE_SELECT);
    //write_tx_mode(cm, &pcsPtr->tx_mode, wb);


    if (pcsPtr->allow_comp_inter_inter) {
        const int32_t use_hybrid_pred = pcsPtr->reference_mode == REFERENCE_MODE_SELECT;

        aom_wb_write_bit(wb, use_hybrid_pred);
    }


    if (pcsPtr->is_skip_mode_allowed) aom_wb_write_bit(wb, pcsPtr->skip_mode_flag);

    if (frame_might_allow_warped_motion(pcsPtr, scsPtr))
        aom_wb_write_bit(wb, pcsPtr->allow_warped_motion);
    else
        assert(!pcsPtr->allow_warped_motion);

    aom_wb_write_bit(wb, pcsPtr->reduced_tx_set_used);

    if (!frame_is_intra_only(pcsPtr)) {
        //  printf("ERROR[AN]: Global motion not supported yet\n");
        WriteGlobalMotion(pcsPtr, wb);
    }



    if (scsPtr->film_grain_params_present && (pcsPtr->showFrame || pcsPtr->showable_frame)) {
        write_film_grain_params(pcsPtr, wb);
    }
}

uint32_t WriteObuHeader(obuType obuType, int32_t obuExtension,
    uint8_t *const dst) {
    struct aom_write_bit_buffer wb = { dst, 0 };
    uint32_t size = 0;

    aom_wb_write_literal(&wb, 0, 1);  // forbidden bit.
    aom_wb_write_literal(&wb, (int32_t)obuType, 4);
    aom_wb_write_literal(&wb, obuExtension ? 1 : 0, 1);
    aom_wb_write_literal(&wb, 1, 1);  // obu_has_payload_length_field
    aom_wb_write_literal(&wb, 0, 1);  // reserved

    if (obuExtension) {
        aom_wb_write_literal(&wb, obuExtension & 0xFF, 8);
    }

    size = aom_wb_bytes_written(&wb);
    return size;
}
int32_t WriteUlebObuSize(uint32_t obuHeaderSize, uint32_t obuPayloadSize,
    uint8_t *dest) {
    const uint32_t obuSize = obuPayloadSize;
    const uint32_t offset = obuHeaderSize;
    size_t codedObuSize = 0;

    if (aom_uleb_encode(obuSize, sizeof(obuSize), dest + offset,
        &codedObuSize) != 0) {
        return AOM_CODEC_ERROR;
    }

    return AOM_CODEC_OK;
}
static size_t ObuMemMove(
    uint32_t obuHeaderSize,
    uint32_t obuPayloadSize,
    uint8_t *data)
{
    const size_t lengthFieldSize = aom_uleb_size_in_bytes(obuPayloadSize);
    const uint32_t moveDstOffset =
        (uint32_t)lengthFieldSize + obuHeaderSize;
    const uint32_t moveSrcOffset = obuHeaderSize;
    const uint32_t moveSize = obuPayloadSize;
    memmove(data + moveDstOffset, data + moveSrcOffset, moveSize);
    return lengthFieldSize;
}

static void add_trailing_bits(struct aom_write_bit_buffer *wb) {
    if (aom_wb_is_byte_aligned(wb)) {
        aom_wb_write_literal(wb, 0x80, 8);
    }
    else {
        // assumes that the other bits are already 0s
        aom_wb_write_bit(wb, 1);
    }
}
static void write_bitstream_level(BitstreamLevel bl,
    struct aom_write_bit_buffer *wb) {
    uint8_t seq_level_idx = major_minor_to_seq_level_idx(bl);
    assert(is_valid_seq_level_idx(seq_level_idx));
    aom_wb_write_literal(wb, seq_level_idx, LEVEL_BITS);
}
static uint32_t WriteSequenceHeaderObu(
    SequenceControlSet_t *scsPtr,
    uint8_t *const dst,
    uint8_t numberSpatialLayers)
{
    struct aom_write_bit_buffer wb = { dst, 0 };
    uint32_t size = 0;

    SetBitstreamLevelTier(scsPtr);

    WriteProfile((BITSTREAM_PROFILE)scsPtr->static_config.profile, &wb);

    // Still picture or not
    aom_wb_write_bit(&wb, scsPtr->still_picture);
    assert(IMPLIES(!scsPtr->still_picture,
        !scsPtr->reduced_still_picture_hdr));

    // whether to use reduced still picture header
    aom_wb_write_bit(&wb, scsPtr->reduced_still_picture_hdr);

    if (scsPtr->reduced_still_picture_hdr) {
        printf("ERROR[AN]: reduced_still_picture_hdr not supported\n");
        //write_bitstream_level(cm->seq_params.level[0], &wb);
    }
    else {
#if AV1_UPGRADE
        aom_wb_write_bit(&wb, scsPtr->timing_info_present);  // timing info present flag

        if (scsPtr->timing_info_present) {
            // timing_info
            printf("ERROR[AN]: timing_info_present not supported\n");
            /*write_timing_info_header(cm, &wb);
            aom_wb_write_bit(&wb, cm->decoder_model_info_present_flag);
            if (cm->decoder_model_info_present_flag) write_decoder_model_info(cm, &wb);*/
        }
        aom_wb_write_bit(&wb, scsPtr->display_model_info_present_flag);


        uint8_t operating_points_cnt_minus_1 =
            numberSpatialLayers > 1 ? numberSpatialLayers - 1 : 0;
        aom_wb_write_literal(&wb, operating_points_cnt_minus_1,
            OP_POINTS_CNT_MINUS_1_BITS);
        int32_t i;
        for (i = 0; i < operating_points_cnt_minus_1 + 1; i++) {
            aom_wb_write_literal(&wb, scsPtr->operating_point_idc[i],
                OP_POINTS_IDC_BITS);
            write_bitstream_level(scsPtr->level[i], &wb);
            if (scsPtr->level[i].major > 3)
                aom_wb_write_bit(&wb, scsPtr->tier[i]);
            if (scsPtr->decoder_model_info_present_flag) {
                printf("ERROR[AN]: decoder_model_info_present_flag not supported\n");
                //aom_wb_write_bit(&wb,
                //    cm->op_params[i].decoder_model_param_present_flag);
                //if (cm->op_params[i].decoder_model_param_present_flag)
                //    write_dec_model_op_parameters(cm, &wb, i);
            }
            if (scsPtr->display_model_info_present_flag) {
                printf("ERROR[AN]: display_model_info_present_flag not supported\n");
                //aom_wb_write_bit(&wb,
                //    cm->op_params[i].display_model_param_present_flag);
                //if (cm->op_params[i].display_model_param_present_flag) {
                //    assert(cm->op_params[i].initial_display_delay <= 10);
                //    aom_wb_write_literal(&wb, cm->op_params[i].initial_display_delay - 1,
                //        4);
                //}
            }
        }
    }
    WriteSequenceHeader(scsPtr, &wb);

    write_color_config(scsPtr, &wb);
#else

        uint8_t operating_points_cnt_minus_1 =
            numberSpatialLayers > 1 ? numberSpatialLayers - 1 : 0;
        aom_wb_write_literal(&wb, operating_points_cnt_minus_1,
            OP_POINTS_CNT_MINUS_1_BITS);
        int32_t i;
        if (operating_points_cnt_minus_1 == 0) {
            scsPtr->operating_point_idc[0] = 0;
        }
        else {
            printf("ERROR[AN]: more than 1 operating_point_idc not supported\n");
            // Set operating_point_idc[] such that for the i-th operating point the
            // first (operating_points_cnt-i) spatial layers and the first temporal
            // layer are decoded Note that highest quality operating point should come
            // first
            /*for (i = 0; i < operating_points_cnt_minus_1 + 1; i++)
                cm->seq_params.operating_point_idc[i] =
                (~(~0u << (operating_points_cnt_minus_1 + 1 - i)) << 8) | 1;*/
        }

        for (i = 0; i < operating_points_cnt_minus_1 + 1; i++) {
            aom_wb_write_literal(&wb, scsPtr->operating_point_idc[i],
                OP_POINTS_IDC_BITS);
            write_bitstream_level(scsPtr->level[i], &wb);
            if (scsPtr->level[i].major > 3)
                aom_wb_write_bit(&wb, scsPtr->tier[i]);
        }
}


    WriteSequenceHeader(scsPtr, &wb);

    // color_config
    WriteBitdepthColorspaceSampling(scsPtr, &wb);

    // timing_info
    if (!scsPtr->reduced_still_picture_hdr)
        aom_wb_write_bit(&wb, scsPtr->timing_info_present);  // timing info present flag
    else
        assert(scsPtr->timing_info_present == 0);

    if (scsPtr->timing_info_present) {
        // timing_info
        printf("ERROR[AN]: timing_info_present not supported\n");
        /*write_timing_info_header(cm, &wb);
        aom_wb_write_bit(&wb, cm->decoder_model_info_present_flag);
        if (cm->decoder_model_info_present_flag) write_decoder_model_info(cm, &wb);*/
    }

    if (scsPtr->operating_points_decoder_model_cnt > 0) {
        printf("ERROR[AN]: operating_points_decoder_model_cnt not supported\n");
        //aom_wb_write_bit(&wb, 1);
        //aom_wb_write_literal(&wb, cm->operating_points_decoder_model_cnt - 1, 5);
    }
    else {
        aom_wb_write_bit(&wb, 0);
    }
    /*for (int32_t op_num = 0; op_num < scsPtr->operating_points_decoder_model_cnt;
       ++op_num) {
        aom_wb_write_literal(
            &wb, cm->op_params[op_num].decoder_model_operating_point_idc, 12);
        aom_wb_write_bit(&wb,
                     cm->op_params[op_num].display_model_param_present_flag);
        if (cm->op_params[op_num].display_model_param_present_flag)
        aom_wb_write_literal(&wb, cm->op_params[op_num].initial_display_delay - 1,
                           4);
        if (cm->decoder_model_info_present_flag) {
            aom_wb_write_bit(&wb,
                       cm->op_params[op_num].decoder_model_param_present_flag);
            if (cm->op_params[op_num].decoder_model_param_present_flag)
                write_dec_model_op_parameters(cm, &wb, op_num);
        }
  }*/


#endif
    aom_wb_write_bit(&wb, scsPtr->film_grain_params_present);

    add_trailing_bits(&wb);

    size = aom_wb_bytes_written(&wb);
    return size;
}
#if TILES
static uint32_t write_tile_group_header(uint8_t *const dst, int startTile,
    int endTile, int tiles_log2,
    int tile_start_and_end_present_flag)
{
    struct aom_write_bit_buffer wb = { dst, 0 };
    uint32_t size = 0;

    if (!tiles_log2) return size;

    aom_wb_write_bit(&wb, tile_start_and_end_present_flag);

    if (tile_start_and_end_present_flag) {
        aom_wb_write_literal(&wb, startTile, tiles_log2);
        aom_wb_write_literal(&wb, endTile, tiles_log2);
    }

    size = aom_wb_bytes_written(&wb);
    return size;
}
#endif
static uint32_t WriteFrameHeaderObu(
    SequenceControlSet_t      *scsPtr,
    PictureParentControlSet_t *pcsPtr,
    uint8_t                   *const dst,
    uint8_t showExisting,
    int32_t appendTrailingBits

)
{
    struct aom_write_bit_buffer wb = { dst, 0 };
    uint32_t totalSize = 0;

    WriteUncompressedHeaderObu(scsPtr, pcsPtr,/* saved_wb,*/ &wb, showExisting);

    if (appendTrailingBits)
        add_trailing_bits(&wb);

    if (showExisting) {
        totalSize = aom_wb_bytes_written(&wb);
        return totalSize;
    }


    totalSize = aom_wb_bytes_written(&wb);
    return totalSize;
}

/**************************************************
* EncodeFrameHeaderHeader
**************************************************/
EbErrorType WriteFrameHeaderAv1(
    Bitstream_t *bitstreamPtr,
    SequenceControlSet_t *scsPtr,
    PictureControlSet_t *pcsPtr,
    uint8_t showExisting)
{
    EbErrorType                 return_error = EB_ErrorNone;
    OutputBitstreamUnit_t       *outputBitstreamPtr = (OutputBitstreamUnit_t*)bitstreamPtr->outputBitstreamPtr;
    PictureParentControlSet_t   *parentPcsPtr = pcsPtr->parent_pcs_ptr;
    uint8_t                     *data = outputBitstreamPtr->bufferAv1;
    uint32_t obuHeaderSize = 0;

    int32_t currDataSize = 0;

    const uint8_t obuExtensionHeader = 0;

    // A new tile group begins at this tile.  Write the obu header and
    // tile group header
    const obuType obuType = showExisting ? OBU_FRAME_HEADER : OBU_FRAME;
    currDataSize =
        WriteObuHeader(obuType, obuExtensionHeader, data);
    obuHeaderSize = currDataSize;

    currDataSize +=
        WriteFrameHeaderObu(scsPtr, parentPcsPtr, /*saved_wb,*/ data + currDataSize, showExisting, showExisting);

#if TILES   
    const int n_log2_tiles = parentPcsPtr->av1_cm->log2_tile_rows + parentPcsPtr->av1_cm->log2_tile_cols;
    int tile_start_and_end_present_flag = 0;

   currDataSize += write_tile_group_header(data + currDataSize,0,
        0, n_log2_tiles, tile_start_and_end_present_flag
        );
#else
    //currDataSize += write_tile_group_header(
    //    data + currDataSize, tile_idx,
    //    AOMMIN(tile_idx + tg_size - 1, tile_cols * tile_rows - 1),
    //    n_log2_tiles, cm->num_tg > 1);
#endif

    if (!showExisting) {
        // Add data from EC stream to Picture Stream.
#if TILES
        int32_t frameSize = parentPcsPtr->av1_cm->tile_cols*parentPcsPtr->av1_cm->tile_rows==1 ? pcsPtr->entropy_coder_ptr->ecWriter.pos : pcsPtr->entropy_coder_ptr->ec_frame_size;
#else
        int32_t frameSize = pcsPtr->entropy_coder_ptr->ecWriter.pos;
#endif
        OutputBitstreamUnit_t *ecOutputBitstreamPtr = (OutputBitstreamUnit_t*)pcsPtr->entropy_coder_ptr->ecOutputBitstreamPtr;
        //****************************************************************//
        // Copy from EC stream to frame stream
        memcpy(data + currDataSize, ecOutputBitstreamPtr->bufferBeginAv1, frameSize);
        currDataSize += (frameSize);
    }
    const uint32_t obuPayloadSize = currDataSize - obuHeaderSize;
    const size_t lengthFieldSize =
        ObuMemMove(obuHeaderSize, obuPayloadSize, data);
    if (WriteUlebObuSize(obuHeaderSize, obuPayloadSize, data) !=
        AOM_CODEC_OK) {
        assert(0);
    }
    currDataSize += (int32_t)lengthFieldSize;
    data += currDataSize;

    outputBitstreamPtr->bufferAv1 = data;
    return return_error;
}

/**************************************************
* EncodeSPSAv1
**************************************************/
EbErrorType EncodeSPSAv1(
    Bitstream_t *bitstreamPtr,
    SequenceControlSet_t *scsPtr)
{
    EbErrorType            return_error = EB_ErrorNone;
    OutputBitstreamUnit_t  *outputBitstreamPtr = (OutputBitstreamUnit_t*)bitstreamPtr->outputBitstreamPtr;
    uint8_t                *data = outputBitstreamPtr->bufferAv1;
    uint32_t                obuHeaderSize = 0;
    uint32_t                obuPayloadSize = 0;
    const uint8_t enhancementLayersCnt = 0;// cm->enhancementLayersCnt;

    // write sequence header obu if KEY_FRAME, preceded by 4-byte size
    obuHeaderSize = WriteObuHeader(OBU_SEQUENCE_HEADER, 0, data);

    obuPayloadSize = WriteSequenceHeaderObu(scsPtr,/*cpi,*/ data + obuHeaderSize,
        enhancementLayersCnt);

    const size_t lengthFieldSize = ObuMemMove(obuHeaderSize, obuPayloadSize, data);
    if (WriteUlebObuSize(obuHeaderSize, obuPayloadSize, data) !=
        AOM_CODEC_OK) {
        // return AOM_CODEC_ERROR;
    }

    data += obuHeaderSize + obuPayloadSize + lengthFieldSize;
    outputBitstreamPtr->bufferAv1 = data;
    return return_error;
}
/**************************************************
* EncodeTDAv1
**************************************************/
EbErrorType EncodeTDAv1(
    Bitstream_t *bitstreamPtr)
{
    EbErrorType            return_error = EB_ErrorNone;
    OutputBitstreamUnit_t   *outputBitstreamPtr = (OutputBitstreamUnit_t*)bitstreamPtr->outputBitstreamPtr;
    uint8_t                 *data = outputBitstreamPtr->bufferAv1;

    // move data and insert OBU_TD preceded by optional 4 byte size
    uint32_t obuHeaderSize = 1;
    const uint32_t obuPayloadSize = 0;
    const size_t lengthFieldSize =
        aom_uleb_size_in_bytes(obuPayloadSize);

    /*if (ctx->pending_cx_data) {
    const size_t move_offset = lengthFieldSize + 1;
    memmove(ctx->pending_cx_data + move_offset, ctx->pending_cx_data,
    ctx->pending_cx_data_sz);
    }*/
    obuHeaderSize = WriteObuHeader(
        OBU_TEMPORAL_DELIMITER, 0, data);

    // OBUs are preceded/succeeded by an unsigned leb128 coded integer.
    if (WriteUlebObuSize(obuHeaderSize, obuPayloadSize, data) != AOM_CODEC_OK) {
        //return AOM_CODEC_ERROR;
    }
    data += obuHeaderSize + obuPayloadSize + lengthFieldSize;
    outputBitstreamPtr->bufferAv1 = data;
    return return_error;
}

#if ADD_DELTA_QP_SUPPORT
static void Av1writeDeltaQindex(
    FRAME_CONTEXT *frameContext,
    int32_t           delta_qindex,
    aom_writer    *w)
{
    int32_t sign = delta_qindex < 0;
    int32_t abs = sign ? -delta_qindex : delta_qindex;
    int32_t rem_bits, thr;
    int32_t smallval = abs < DELTA_Q_SMALL ? 1 : 0;
    //FRAME_CONTEXT *ec_ctx = xd->tile_ctx;

    aom_write_symbol(w, AOMMIN(abs, DELTA_Q_SMALL), frameContext->delta_q_cdf,
        DELTA_Q_PROBS + 1);

    if (!smallval) {
        rem_bits = OD_ILOG_NZ(abs - 1) - 1;
        thr = (1 << rem_bits) + 1;
        aom_write_literal(w, rem_bits - 1, 3);
        aom_write_literal(w, abs - thr, rem_bits);
    }
    if (abs > 0) {
        aom_write_bit(w, sign);
    }
}
#endif


static void write_cdef(
    SequenceControlSet_t     *seqCSetPtr,
    PictureControlSet_t     *p_pcs_ptr,
    //Av1Common *cm,
    MacroBlockD *const xd,
    aom_writer *w,
    int32_t skip, int32_t mi_col, int32_t mi_row)
{

    (void)xd;
    Av1Common *cm = p_pcs_ptr->parent_pcs_ptr->av1_cm;


    if (p_pcs_ptr->parent_pcs_ptr->coded_lossless || p_pcs_ptr->parent_pcs_ptr->allow_intrabc) {
        // Initialize to indicate no CDEF for safety.
        p_pcs_ptr->parent_pcs_ptr->cdef_bits = 0;
        p_pcs_ptr->parent_pcs_ptr->cdef_strengths[0] = 0;
        p_pcs_ptr->parent_pcs_ptr->nb_cdef_strengths = 1;
        p_pcs_ptr->parent_pcs_ptr->cdef_uv_strengths[0] = 0;
        return;
    }

    const int32_t m = ~((1 << (6 - MI_SIZE_LOG2)) - 1);
    const ModeInfo *mi =
        p_pcs_ptr->mi_grid_base[(mi_row & m) * cm->mi_stride + (mi_col & m)];
    //cm->mi_grid_visible[(mi_row & m) * cm->mi_stride + (mi_col & m)];

// Initialise when at top left part of the superblock
    if (!(mi_row & (seqCSetPtr->mib_size - 1)) &&
        !(mi_col & (seqCSetPtr->mib_size - 1))) {  // Top left?
        p_pcs_ptr->cdef_preset[0] = p_pcs_ptr->cdef_preset[1] = p_pcs_ptr->cdef_preset[2] =
            p_pcs_ptr->cdef_preset[3] = -1;
    }

    // Emit CDEF param at first non-skip coding block
    const int32_t mask = 1 << (6 - MI_SIZE_LOG2);
    const int32_t index = seqCSetPtr->sb_size == BLOCK_128X128
        ? !!(mi_col & mask) + 2 * !!(mi_row & mask)
        : 0;

    if (p_pcs_ptr->cdef_preset[index] == -1 && !skip) {
        aom_write_literal(w, mi->mbmi.cdef_strength, p_pcs_ptr->parent_pcs_ptr->cdef_bits);
        p_pcs_ptr->cdef_preset[index] = mi->mbmi.cdef_strength;


    }

}


void av1_reset_loop_restoration(PictureControlSet_t     *piCSetPtr) {
    for (int32_t p = 0; p < 3; ++p) {
        set_default_wiener(piCSetPtr->wiener_info + p);
        set_default_sgrproj(piCSetPtr->sgrproj_info + p);
    }
}
static void write_wiener_filter(int32_t wiener_win, const WienerInfo *wiener_info,
    WienerInfo *ref_wiener_info, aom_writer *wb) {
    if (wiener_win == WIENER_WIN)
        aom_write_primitive_refsubexpfin(
            wb, WIENER_FILT_TAP0_MAXV - WIENER_FILT_TAP0_MINV + 1,
            WIENER_FILT_TAP0_SUBEXP_K,
            ref_wiener_info->vfilter[0] - WIENER_FILT_TAP0_MINV,
            wiener_info->vfilter[0] - WIENER_FILT_TAP0_MINV);
    else
        assert(wiener_info->vfilter[0] == 0 &&
            wiener_info->vfilter[WIENER_WIN - 1] == 0);
    aom_write_primitive_refsubexpfin(
        wb, WIENER_FILT_TAP1_MAXV - WIENER_FILT_TAP1_MINV + 1,
        WIENER_FILT_TAP1_SUBEXP_K,
        ref_wiener_info->vfilter[1] - WIENER_FILT_TAP1_MINV,
        wiener_info->vfilter[1] - WIENER_FILT_TAP1_MINV);
    aom_write_primitive_refsubexpfin(
        wb, WIENER_FILT_TAP2_MAXV - WIENER_FILT_TAP2_MINV + 1,
        WIENER_FILT_TAP2_SUBEXP_K,
        ref_wiener_info->vfilter[2] - WIENER_FILT_TAP2_MINV,
        wiener_info->vfilter[2] - WIENER_FILT_TAP2_MINV);
    if (wiener_win == WIENER_WIN)
        aom_write_primitive_refsubexpfin(
            wb, WIENER_FILT_TAP0_MAXV - WIENER_FILT_TAP0_MINV + 1,
            WIENER_FILT_TAP0_SUBEXP_K,
            ref_wiener_info->hfilter[0] - WIENER_FILT_TAP0_MINV,
            wiener_info->hfilter[0] - WIENER_FILT_TAP0_MINV);
    else
        assert(wiener_info->hfilter[0] == 0 &&
            wiener_info->hfilter[WIENER_WIN - 1] == 0);
    aom_write_primitive_refsubexpfin(
        wb, WIENER_FILT_TAP1_MAXV - WIENER_FILT_TAP1_MINV + 1,
        WIENER_FILT_TAP1_SUBEXP_K,
        ref_wiener_info->hfilter[1] - WIENER_FILT_TAP1_MINV,
        wiener_info->hfilter[1] - WIENER_FILT_TAP1_MINV);
    aom_write_primitive_refsubexpfin(
        wb, WIENER_FILT_TAP2_MAXV - WIENER_FILT_TAP2_MINV + 1,
        WIENER_FILT_TAP2_SUBEXP_K,
        ref_wiener_info->hfilter[2] - WIENER_FILT_TAP2_MINV,
        wiener_info->hfilter[2] - WIENER_FILT_TAP2_MINV);
    memcpy(ref_wiener_info, wiener_info, sizeof(*wiener_info));
}

static void write_sgrproj_filter(const SgrprojInfo *sgrproj_info,
    SgrprojInfo *ref_sgrproj_info,
    aom_writer *wb) {
    aom_write_literal(wb, sgrproj_info->ep, SGRPROJ_PARAMS_BITS);
    const sgr_params_type *params = &sgr_params[sgrproj_info->ep];

    if (params->r[0] == 0) {
        assert(sgrproj_info->xqd[0] == 0);
        aom_write_primitive_refsubexpfin(
            wb, SGRPROJ_PRJ_MAX1 - SGRPROJ_PRJ_MIN1 + 1, SGRPROJ_PRJ_SUBEXP_K,
            (uint16_t)(ref_sgrproj_info->xqd[1] - SGRPROJ_PRJ_MIN1),
            (uint16_t)(sgrproj_info->xqd[1] - SGRPROJ_PRJ_MIN1));
    }
    else if (params->r[1] == 0) {
        aom_write_primitive_refsubexpfin(
            wb, SGRPROJ_PRJ_MAX0 - SGRPROJ_PRJ_MIN0 + 1, SGRPROJ_PRJ_SUBEXP_K,
            (uint16_t)(ref_sgrproj_info->xqd[0] - SGRPROJ_PRJ_MIN0),
            (uint16_t)(sgrproj_info->xqd[0] - SGRPROJ_PRJ_MIN0));
    }
    else {
        aom_write_primitive_refsubexpfin(
            wb, SGRPROJ_PRJ_MAX0 - SGRPROJ_PRJ_MIN0 + 1, SGRPROJ_PRJ_SUBEXP_K,
            (uint16_t)(ref_sgrproj_info->xqd[0] - SGRPROJ_PRJ_MIN0),
            (uint16_t)(sgrproj_info->xqd[0] - SGRPROJ_PRJ_MIN0));
        aom_write_primitive_refsubexpfin(
            wb, SGRPROJ_PRJ_MAX1 - SGRPROJ_PRJ_MIN1 + 1, SGRPROJ_PRJ_SUBEXP_K,
            (uint16_t)(ref_sgrproj_info->xqd[1] - SGRPROJ_PRJ_MIN1),
            (uint16_t)(sgrproj_info->xqd[1] - SGRPROJ_PRJ_MIN1));
    }

    memcpy(ref_sgrproj_info, sgrproj_info, sizeof(*sgrproj_info));
}
static void loop_restoration_write_sb_coeffs(PictureControlSet_t     *piCSetPtr, FRAME_CONTEXT           *frameContext, const Av1Common *const cm,
    //MacroBlockD *xd,
    const RestorationUnitInfo *rui,
    aom_writer *const w, int32_t plane/*,
    FRAME_COUNTS *counts*/)
{
    const RestorationInfo *rsi = cm->rst_info + plane;
    RestorationType frame_rtype = rsi->frame_restoration_type;
    if (frame_rtype == RESTORE_NONE) return;

    //(void)counts;
//    assert(!cm->all_lossless);

    const int32_t wiener_win = (plane > 0) ? WIENER_WIN_CHROMA : WIENER_WIN;
    WienerInfo *wiener_info = piCSetPtr->wiener_info + plane;
    SgrprojInfo *sgrproj_info = piCSetPtr->sgrproj_info + plane;
    RestorationType unit_rtype = rui->restoration_type;



    if (frame_rtype == RESTORE_SWITCHABLE) {
        aom_write_symbol(w, unit_rtype, /*xd->tile_ctx->*/frameContext->switchable_restore_cdf,
            RESTORE_SWITCHABLE_TYPES);
#if CONFIG_ENTROPY_STATS
        ++counts->switchable_restore[unit_rtype];
#endif
        switch (unit_rtype) {
        case RESTORE_WIENER:
            write_wiener_filter(wiener_win, &rui->wiener_info, wiener_info, w);
            //printf("POC:%i plane:%i v:%i %i %i  h:%i %i %i\n", piCSetPtr->picture_number, plane, rui->wiener_info.vfilter[0], rui->wiener_info.vfilter[1], rui->wiener_info.vfilter[2], rui->wiener_info.hfilter[0], rui->wiener_info.hfilter[1], rui->wiener_info.hfilter[2]);
            break;
        case RESTORE_SGRPROJ:
            write_sgrproj_filter(&rui->sgrproj_info, sgrproj_info, w);
            //printf("POC:%i plane:%i ep:%i xqd_0:%i  xqd_1:%i\n", piCSetPtr->picture_number, plane, rui->sgrproj_info.ep, rui->sgrproj_info.xqd[0], rui->sgrproj_info.xqd[1]);
            break;
        default: assert(unit_rtype == RESTORE_NONE);// printf("POC:%i plane:%i OFF\n", piCSetPtr->picture_number, plane);
            break;
        }
    }
    else if (frame_rtype == RESTORE_WIENER) {
        aom_write_symbol(w, unit_rtype != RESTORE_NONE,
            /*xd->tile_ctx->*/frameContext->wiener_restore_cdf, 2);
#if CONFIG_ENTROPY_STATS
        ++counts->wiener_restore[unit_rtype != RESTORE_NONE];
#endif
        if (unit_rtype != RESTORE_NONE) {
            write_wiener_filter(wiener_win, &rui->wiener_info, wiener_info, w);
            //printf("POC:%i plane:%i v:%i %i %i  h:%i %i %i\n", piCSetPtr->picture_number, plane, rui->wiener_info.vfilter[0], rui->wiener_info.vfilter[1], rui->wiener_info.vfilter[2], rui->wiener_info.hfilter[0], rui->wiener_info.hfilter[1], rui->wiener_info.hfilter[2]);
        }
        //else
            //printf("POC:%i plane:%i OFF\n", piCSetPtr->picture_number, plane);
    }
    else if (frame_rtype == RESTORE_SGRPROJ) {
        aom_write_symbol(w, unit_rtype != RESTORE_NONE,
            /*xd->tile_ctx->*/frameContext->sgrproj_restore_cdf, 2);
#if CONFIG_ENTROPY_STATS
        ++counts->sgrproj_restore[unit_rtype != RESTORE_NONE];
#endif
        if (unit_rtype != RESTORE_NONE) {
            write_sgrproj_filter(&rui->sgrproj_info, sgrproj_info, w);
            //printf("POC:%i plane:%i ep:%i xqd_0:%i  xqd_1:%i\n", piCSetPtr->picture_number, plane, rui->sgrproj_info.ep, rui->sgrproj_info.xqd[0], rui->sgrproj_info.xqd[1]);
        }
        //else
        //    printf("POC:%i plane:%i OFF\n", piCSetPtr->picture_number, plane);
    }
}

EbErrorType ec_update_neighbors(
    PictureControlSet_t     *picture_control_set_ptr,
    EntropyCodingContext_t  *context_ptr,
    uint32_t                 blkOriginX,
    uint32_t                 blkOriginY,
    CodingUnit_t            *cu_ptr,
    BlockSize                bsize,
    EbPictureBufferDesc_t   *coeffPtr)
{
    UNUSED(coeffPtr);
    EbErrorType return_error = EB_ErrorNone;
    NeighborArrayUnit_t     *mode_type_neighbor_array = picture_control_set_ptr->mode_type_neighbor_array;
    NeighborArrayUnit_t     *partition_context_neighbor_array = picture_control_set_ptr->partition_context_neighbor_array;
    NeighborArrayUnit_t     *skip_flag_neighbor_array = picture_control_set_ptr->skip_flag_neighbor_array;
    NeighborArrayUnit_t     *skip_coeff_neighbor_array = picture_control_set_ptr->skip_coeff_neighbor_array;
    NeighborArrayUnit_t     *luma_dc_sign_level_coeff_neighbor_array = picture_control_set_ptr->luma_dc_sign_level_coeff_neighbor_array;
    NeighborArrayUnit_t     *cr_dc_sign_level_coeff_neighbor_array = picture_control_set_ptr->cr_dc_sign_level_coeff_neighbor_array;
    NeighborArrayUnit_t     *cb_dc_sign_level_coeff_neighbor_array = picture_control_set_ptr->cb_dc_sign_level_coeff_neighbor_array;
    NeighborArrayUnit_t     *inter_pred_dir_neighbor_array = picture_control_set_ptr->inter_pred_dir_neighbor_array;
    NeighborArrayUnit_t     *ref_frame_type_neighbor_array = picture_control_set_ptr->ref_frame_type_neighbor_array;
    NeighborArrayUnit32_t   *interpolation_type_neighbor_array = picture_control_set_ptr->interpolation_type_neighbor_array;
    const BlockGeom         *blk_geom = Get_blk_geom_mds(cu_ptr->mds_idx);
    EbBool                   skipCoeff = EB_FALSE;
    PartitionContext         partition;

    skipCoeff = cu_ptr->block_has_coeff ? 0 : 1;
    // Update the Leaf Depth Neighbor Array
    partition.above = partition_context_lookup[bsize].above;
    partition.left = partition_context_lookup[bsize].left;

    NeighborArrayUnitModeWrite(
        partition_context_neighbor_array,
        (uint8_t*)&partition,
        blkOriginX,
        blkOriginY,
        blk_geom->bwidth,
        blk_geom->bheight,
        NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

    // Update the Mode Type Neighbor Array
    {
        uint8_t prediction_mode_flag = (uint8_t)cu_ptr->prediction_mode_flag;
        NeighborArrayUnitModeWrite(
            mode_type_neighbor_array,
            &prediction_mode_flag,
            blkOriginX,
            blkOriginY,
            blk_geom->bwidth,
            blk_geom->bheight,
            NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
    }

    // Update the Skip Flag Neighbor Array
    {
        uint8_t skip_flag = (uint8_t)cu_ptr->skip_flag;
        NeighborArrayUnitModeWrite(
            skip_flag_neighbor_array,
            (uint8_t*)&skip_flag,
            blkOriginX,
            blkOriginY,
            blk_geom->bwidth,
            blk_geom->bheight,
            NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
    }

    // Update the Skip Coeff Neighbor Array
    {
        //
        NeighborArrayUnitModeWrite(
            skip_coeff_neighbor_array,
            (uint8_t*)&skipCoeff,
            blkOriginX,
            blkOriginY,
            blk_geom->bwidth,
            blk_geom->bheight,
            NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
    }


    if (skipCoeff)
    {
        uint8_t dcSignLevelCoeff = 0;

        NeighborArrayUnitModeWrite(
            luma_dc_sign_level_coeff_neighbor_array,
            (uint8_t*)&dcSignLevelCoeff,
            blkOriginX,
            blkOriginY,
            blk_geom->bwidth,
            blk_geom->bheight,
            NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

        if (blk_geom->has_uv)
            NeighborArrayUnitModeWrite(
                cb_dc_sign_level_coeff_neighbor_array,
                (uint8_t*)&dcSignLevelCoeff,
                ((blkOriginX >> 3) << 3) >> 1,
                ((blkOriginY >> 3) << 3) >> 1,
                blk_geom->bwidth_uv,
                blk_geom->bheight_uv,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

        if (blk_geom->has_uv)
            NeighborArrayUnitModeWrite(
                cr_dc_sign_level_coeff_neighbor_array,
                (uint8_t*)&dcSignLevelCoeff,
                ((blkOriginX >> 3) << 3) >> 1,
                ((blkOriginY >> 3) << 3) >> 1,
                blk_geom->bwidth_uv,
                blk_geom->bheight_uv,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

        context_ptr->coded_area_sb += blk_geom->bwidth * blk_geom->bheight;

        if (blk_geom->has_uv)
            context_ptr->coded_area_sb_uv += blk_geom->bwidth_uv * blk_geom->bheight_uv;

    }

    // Update the Inter Pred Type Neighbor Array
    {
        uint8_t inter_pred_direction_index = (uint8_t)cu_ptr->prediction_unit_array->inter_pred_direction_index;
        NeighborArrayUnitModeWrite(
            inter_pred_dir_neighbor_array,
            (uint8_t*)&(inter_pred_direction_index),
            blkOriginX,
            blkOriginY,
            blk_geom->bwidth,
            blk_geom->bheight,
            NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
    }

    // Update the refFrame Type Neighbor Array
    {
        uint8_t ref_frame_type = (uint8_t)cu_ptr->prediction_unit_array[0].ref_frame_type;
        NeighborArrayUnitModeWrite(
            ref_frame_type_neighbor_array,
            (uint8_t*)&(ref_frame_type),
            blkOriginX,
            blkOriginY,
            blk_geom->bwidth,
            blk_geom->bheight,
            NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
    }

    // Update the interpolation Type Neighbor Array
    {
        uint32_t interpolationType = cu_ptr->interp_filters;
        NeighborArrayUnitModeWrite32(
            interpolation_type_neighbor_array,
            interpolationType,
            blkOriginX,
            blkOriginY,
            blk_geom->bwidth,
            blk_geom->bheight,
            NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
    }
    return return_error;

}

EbErrorType write_modes_b(
    PictureControlSet_t     *picture_control_set_ptr,
    EntropyCodingContext_t  *context_ptr,
    EntropyCoder_t          *entropy_coder_ptr,
    LargestCodingUnit_t     *tbPtr,
    CodingUnit_t            *cu_ptr,
    EbPictureBufferDesc_t   *coeffPtr)
{
    UNUSED(tbPtr);
    EbErrorType return_error = EB_ErrorNone;
    FRAME_CONTEXT           *frameContext = entropy_coder_ptr->fc;
    aom_writer              *ecWriter = &entropy_coder_ptr->ecWriter;
    SequenceControlSet_t     *sequence_control_set_ptr = (SequenceControlSet_t*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->objectPtr;

    NeighborArrayUnit_t     *mode_type_neighbor_array = picture_control_set_ptr->mode_type_neighbor_array;
    NeighborArrayUnit_t     *intra_luma_mode_neighbor_array = picture_control_set_ptr->intra_luma_mode_neighbor_array;
    NeighborArrayUnit_t     *skip_flag_neighbor_array = picture_control_set_ptr->skip_flag_neighbor_array;
    NeighborArrayUnit_t     *skip_coeff_neighbor_array = picture_control_set_ptr->skip_coeff_neighbor_array;
    NeighborArrayUnit_t     *luma_dc_sign_level_coeff_neighbor_array = picture_control_set_ptr->luma_dc_sign_level_coeff_neighbor_array;
    NeighborArrayUnit_t     *cr_dc_sign_level_coeff_neighbor_array = picture_control_set_ptr->cr_dc_sign_level_coeff_neighbor_array;
    NeighborArrayUnit_t     *cb_dc_sign_level_coeff_neighbor_array = picture_control_set_ptr->cb_dc_sign_level_coeff_neighbor_array;
    NeighborArrayUnit_t     *inter_pred_dir_neighbor_array = picture_control_set_ptr->inter_pred_dir_neighbor_array;
    NeighborArrayUnit_t     *ref_frame_type_neighbor_array = picture_control_set_ptr->ref_frame_type_neighbor_array;
    NeighborArrayUnit32_t   *interpolation_type_neighbor_array = picture_control_set_ptr->interpolation_type_neighbor_array;

    const BlockGeom          *blk_geom = Get_blk_geom_mds(cu_ptr->mds_idx);
    uint32_t blkOriginX = context_ptr->sb_origin_x + blk_geom->origin_x;
    uint32_t blkOriginY = context_ptr->sb_origin_y + blk_geom->origin_y;
    BlockSize bsize = blk_geom->bsize;
    EbBool                   skipCoeff = EB_FALSE;
    skipCoeff = cu_ptr->block_has_coeff ? 0 : 1;

    if (picture_control_set_ptr->slice_type == I_SLICE) {
        //const int32_t skip = write_skip(cm, xd, mbmi->segment_id, mi, w);
        EncodeSkipCoeffAv1(
            frameContext,
            ecWriter,
            skipCoeff,
            blkOriginX,
            blkOriginY,
            skip_coeff_neighbor_array);

        write_cdef(
            sequence_control_set_ptr,
            picture_control_set_ptr,
            cu_ptr->av1xd,
            ecWriter,
            skipCoeff,
            blkOriginX >> MI_SIZE_LOG2,
            blkOriginY >> MI_SIZE_LOG2);

#if ADD_DELTA_QP_SUPPORT //PART 1
        if (picture_control_set_ptr->parent_pcs_ptr->delta_q_present_flag) {

            int32_t current_q_index = cu_ptr->qp;

            int32_t super_block_upper_left = (((blkOriginY >> 2) & (sequence_control_set_ptr->mib_size - 1)) == 0) && (((blkOriginX >> 2) & (sequence_control_set_ptr->mib_size - 1)) == 0);
            /*((mi_row & (cm->seq_params.mib_size - 1)) == 0) && ((mi_col & (cm->seq_params.mib_size - 1)) == 0);*/

            BlockSize bsize = cu_size == 64 ? BLOCK_64X64 : cu_size == 32 ? BLOCK_32X32 : cu_size == 16 ? BLOCK_16X16 : cu_size == 8 ? BLOCK_8X8 : BLOCK_4X4;
            if (cu_size == 8 && cu_ptr->prediction_mode_flag == INTRA_MODE && cu_ptr->pred_mode == INTRA_MODE_4x4) {
                bsize = BLOCK_4X4;
            }
            if ((bsize != sequence_control_set_ptr->sb_size || skipCoeff == 0) && super_block_upper_left) {
                assert(current_q_index > 0);
                int32_t reduced_delta_qindex = (current_q_index - picture_control_set_ptr->parent_pcs_ptr->prev_qindex) / picture_control_set_ptr->parent_pcs_ptr->delta_q_res;

                //write_delta_qindex(xd, reduced_delta_qindex, w);
                Av1writeDeltaQindex(
                    frameContext,
                    reduced_delta_qindex,
                    ecWriter);
                /*if (picture_control_set_ptr->picture_number == 0){
                printf("%d\t%d\t%d\t%d\n",
                blkOriginX,
                blkOriginY,
                current_q_index,
                picture_control_set_ptr->parent_pcs_ptr->prev_qindex);
                }*/
                picture_control_set_ptr->parent_pcs_ptr->prev_qindex = current_q_index;

            }
        }
#endif


        {
            const uint32_t intra_luma_mode = cu_ptr->pred_mode;
            uint32_t intra_chroma_mode = cu_ptr->prediction_unit_array->intra_chroma_mode;

            EncodeIntraLumaModeAv1(
                frameContext,
                ecWriter,
                cu_ptr,
                blkOriginX,
                blkOriginY,
                intra_luma_mode,
                mode_type_neighbor_array,
                intra_luma_mode_neighbor_array);

            NeighborArrayUnitModeWrite(
                intra_luma_mode_neighbor_array,
                (uint8_t*)&intra_luma_mode,
                blkOriginX,
                blkOriginY,
                blk_geom->bwidth,
                blk_geom->bheight,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);


            if (blk_geom->has_uv)
                EncodeIntraChromaModeAv1(
                    frameContext,
                    ecWriter,
                    cu_ptr,
                    intra_luma_mode,
                    intra_chroma_mode,
                    blk_geom->bwidth <= 32 && blk_geom->bheight <= 32);

            if (!skipCoeff) {
                Av1EncodeCoeff1D(
                    picture_control_set_ptr,
                    context_ptr,
                    frameContext,
                    ecWriter,
                    cu_ptr,
                    blkOriginX,
                    blkOriginY,
                    intra_luma_mode,
                    bsize,
                    coeffPtr,
                    luma_dc_sign_level_coeff_neighbor_array,
                    cr_dc_sign_level_coeff_neighbor_array,
                    cb_dc_sign_level_coeff_neighbor_array);
            }
        }
    }
    else {
        if (picture_control_set_ptr->parent_pcs_ptr->skip_mode_flag && is_comp_ref_allowed(bsize)) {
            EncodeSkipModeAv1(
                frameContext,
                ecWriter,
                cu_ptr->skip_flag,
                blkOriginX,
                blkOriginY,
                skip_flag_neighbor_array);
        }
        if (!picture_control_set_ptr->parent_pcs_ptr->skip_mode_flag && cu_ptr->skip_flag) {
            printf("ERROR[AN]: SKIP not supported\n");
        }

        if (!cu_ptr->skip_flag) {
            //const int32_t skip = write_skip(cm, xd, mbmi->segment_id, mi, w);
            EncodeSkipCoeffAv1(
                frameContext,
                ecWriter,
                skipCoeff,
                blkOriginX,
                blkOriginY,
                skip_coeff_neighbor_array);
        }

        write_cdef(
            sequence_control_set_ptr,
            picture_control_set_ptr, /*cm,*/
            cu_ptr->av1xd,
            ecWriter,
            cu_ptr->skip_flag ? 1 : skipCoeff,
            blkOriginX >> MI_SIZE_LOG2,
            blkOriginY >> MI_SIZE_LOG2);
#if ADD_DELTA_QP_SUPPORT//PART 2
        if (picture_control_set_ptr->parent_pcs_ptr->delta_q_present_flag) {
            int32_t current_q_index = cu_ptr->qp;

            int32_t super_block_upper_left = (((blkOriginY >> 2) & (sequence_control_set_ptr->mib_size - 1)) == 0) && (((blkOriginX >> 2) & (sequence_control_set_ptr->mib_size - 1)) == 0);
            /*((mi_row & (cm->seq_params.mib_size - 1)) == 0) && ((mi_col & (cm->seq_params.mib_size - 1)) == 0);*/

            BlockSize bsize = cu_size == 64 ? BLOCK_64X64 : cu_size == 32 ? BLOCK_32X32 : cu_size == 16 ? BLOCK_16X16 : cu_size == 8 ? BLOCK_8X8 : BLOCK_4X4;
            if (cu_size == 8 && cu_ptr->prediction_mode_flag == INTRA_MODE && cu_ptr->pred_mode == INTRA_MODE_4x4) {
                bsize = BLOCK_4X4;
            }

            if ((bsize != sequence_control_set_ptr->sb_size || skipCoeff == 0) && super_block_upper_left) {
                assert(current_q_index > 0);

                int32_t reduced_delta_qindex = (current_q_index - picture_control_set_ptr->parent_pcs_ptr->prev_qindex) / picture_control_set_ptr->parent_pcs_ptr->delta_q_res;

                //write_delta_qindex(xd, reduced_delta_qindex, w);

                Av1writeDeltaQindex(
                    frameContext,
                    reduced_delta_qindex,
                    ecWriter);

                picture_control_set_ptr->parent_pcs_ptr->prev_qindex = current_q_index;
            }
        }

#endif
        if (!cu_ptr->skip_flag) {


            //write_is_inter(cm, xd, mbmi->segment_id, w, is_inter)
            EncodePredModeAv1(
                frameContext,
                ecWriter,
                cu_ptr->prediction_mode_flag,
                blkOriginX,
                blkOriginY,
                mode_type_neighbor_array);

            if (cu_ptr->prediction_mode_flag == INTRA_MODE) {

                uint32_t intra_luma_mode = cu_ptr->pred_mode;
                uint32_t intra_chroma_mode = cu_ptr->prediction_unit_array->intra_chroma_mode;

                EncodeIntraLumaModeNonKeyAv1(
                    frameContext,
                    ecWriter,
                    cu_ptr,
                    bsize,
                    intra_luma_mode);

                NeighborArrayUnitModeWrite(
                    intra_luma_mode_neighbor_array,
                    (uint8_t*)&intra_luma_mode,
                    blkOriginX,
                    blkOriginY,
                    blk_geom->bwidth,
                    blk_geom->bheight,
                    NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);


                if (blk_geom->has_uv)
                    EncodeIntraChromaModeAv1(
                        frameContext,
                        ecWriter,
                        cu_ptr,
                        intra_luma_mode,
                        intra_chroma_mode,
                        blk_geom->bwidth <= 32 && blk_geom->bheight <= 32);
            }
            else {

                Av1CollectNeighborsRefCounts(
                    cu_ptr,
                    blkOriginX,
                    blkOriginY,
                    mode_type_neighbor_array,
                    inter_pred_dir_neighbor_array,
                    ref_frame_type_neighbor_array);


                WriteRefFrames(
                    frameContext,
                    ecWriter,
                    picture_control_set_ptr->parent_pcs_ptr,
                    cu_ptr,
                    bsize,
                    blkOriginX,
                    blkOriginY,
                    mode_type_neighbor_array,
                    inter_pred_dir_neighbor_array);

                MvReferenceFrame rf[2];
                av1_set_ref_frame(rf, cu_ptr->prediction_unit_array[0].ref_frame_type);
                int16_t mode_ctx = Av1ModeContextAnalyzer(cu_ptr->inter_mode_ctx, rf);
                PredictionMode inter_mode = (PredictionMode)cu_ptr->prediction_unit_array[0].inter_mode;
                const int32_t is_compound = (cu_ptr->prediction_unit_array[0].inter_pred_direction_index == BI_PRED);

                // If segment skip is not enabled code the mode.
                if (1) {
                    if (cu_ptr->prediction_unit_array[0].is_compound) {
                        WriteInterCompoundMode(
                            frameContext,
                            ecWriter,
                            inter_mode,
                            mode_ctx);
                    }
                    else /*if (is_inter_singleref_mode(mode))*/
                        WriteInterMode(
                            frameContext,
                            ecWriter,
                            inter_mode,
                            mode_ctx,
                            blkOriginX,
                            blkOriginY);


                    if (inter_mode == NEWMV || inter_mode == NEW_NEWMV || have_nearmv_in_inter_mode(inter_mode)) {
                        WriteDrlIdx(
                            frameContext,
                            ecWriter,
                            cu_ptr);

                    }
                }

                if (inter_mode == NEWMV || inter_mode == NEW_NEWMV) {
                    IntMv ref_mv;

                    for (uint8_t ref = 0; ref < 1 + is_compound; ++ref) {
                        nmv_context *nmvc = &frameContext->nmvc;
                        ref_mv = cu_ptr->predmv[ref];

                        MV mv;
                        mv.row = cu_ptr->prediction_unit_array[0].mv[ref].y;
                        mv.col = cu_ptr->prediction_unit_array[0].mv[ref].x;
                        if (cu_ptr->prediction_unit_array[0].inter_pred_direction_index == UNI_PRED_LIST_1)
                        {
                            mv.row = cu_ptr->prediction_unit_array[0].mv[1].y;
                            mv.col = cu_ptr->prediction_unit_array[0].mv[1].x;
                        }

                        av1_encode_mv(
                            picture_control_set_ptr->parent_pcs_ptr,
                            ecWriter,
                            &mv,
                            &ref_mv.as_mv,
                            nmvc,
                            picture_control_set_ptr->parent_pcs_ptr->allow_high_precision_mv);
                    }
                }
                else if (inter_mode == NEAREST_NEWMV || inter_mode == NEAR_NEWMV) {
                    nmv_context *nmvc = &frameContext->nmvc;
                    IntMv ref_mv = cu_ptr->predmv[1];

                    MV mv;
                    mv.row = cu_ptr->prediction_unit_array[0].mv[1].y;
                    mv.col = cu_ptr->prediction_unit_array[0].mv[1].x;

                    av1_encode_mv(
                        picture_control_set_ptr->parent_pcs_ptr,
                        ecWriter,
                        &mv,
                        &ref_mv.as_mv,
                        nmvc,
                        picture_control_set_ptr->parent_pcs_ptr->allow_high_precision_mv);

                }
                else if (inter_mode == NEW_NEARESTMV || inter_mode == NEW_NEARMV) {
                    nmv_context *nmvc = &frameContext->nmvc;
                    IntMv ref_mv = cu_ptr->predmv[0];

                    MV mv;
                    mv.row = cu_ptr->prediction_unit_array[0].mv[0].y;
                    mv.col = cu_ptr->prediction_unit_array[0].mv[0].x;

                    av1_encode_mv(
                        picture_control_set_ptr->parent_pcs_ptr,
                        ecWriter,
                        &mv,
                        &ref_mv.as_mv,
                        nmvc,
                        picture_control_set_ptr->parent_pcs_ptr->allow_high_precision_mv);

                }

                if (picture_control_set_ptr->parent_pcs_ptr->switchable_motion_mode
                    && rf[1] != INTRA_FRAME) {
                    write_motion_mode(  
                        frameContext,
                        ecWriter,
                        bsize,
                        cu_ptr->prediction_unit_array[0].motion_mode,
                        rf[0],
                        rf[1],
                        cu_ptr,
                        picture_control_set_ptr);
                }

                if (sequence_control_set_ptr->enable_masked_compound || sequence_control_set_ptr->enable_jnt_comp)
                    printf("ERROR[AN]: masked_compound_used and enable_jnt_comp not supported\n");

                // No filter for Global MV
                write_mb_interp_filter(
                    ref_frame_type_neighbor_array,
                    bsize,
                    rf[0],
                    rf[1],
                    picture_control_set_ptr->parent_pcs_ptr,
                    ecWriter,
                    cu_ptr,
                    entropy_coder_ptr,
                    interpolation_type_neighbor_array,
                    blkOriginX,
                    blkOriginY
                );

            }

            if (!skipCoeff) {
                uint32_t intra_luma_mode = DC_PRED;
                if (cu_ptr->prediction_mode_flag == INTRA_MODE)
                    intra_luma_mode = (uint32_t)cu_ptr->pred_mode;

                {
                    Av1EncodeCoeff1D(
                        picture_control_set_ptr,
                        context_ptr,
                        frameContext,
                        ecWriter,
                        cu_ptr,
                        blkOriginX,
                        blkOriginY,
                        intra_luma_mode,
                        bsize,
                        coeffPtr,
                        luma_dc_sign_level_coeff_neighbor_array,
                        cr_dc_sign_level_coeff_neighbor_array,
                        cb_dc_sign_level_coeff_neighbor_array);
                }
            }
        }
    }
    // Update the neighbors
    ec_update_neighbors(
        picture_control_set_ptr,
        context_ptr,
        blkOriginX,
        blkOriginY,
        cu_ptr,
        bsize,
        coeffPtr);



    return return_error;

}
/**********************************************
* Write sb
**********************************************/
EB_EXTERN EbErrorType write_sb(
    EntropyCodingContext_t  *context_ptr,
    LargestCodingUnit_t     *tbPtr,
    PictureControlSet_t     *picture_control_set_ptr,
    EntropyCoder_t          *entropy_coder_ptr,
    EbPictureBufferDesc_t   *coeffPtr)
{
    EbErrorType return_error = EB_ErrorNone;
    FRAME_CONTEXT           *frameContext = entropy_coder_ptr->fc;
    aom_writer              *ecWriter = &entropy_coder_ptr->ecWriter;
    SequenceControlSet_t     *sequence_control_set_ptr = (SequenceControlSet_t*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->objectPtr;
    NeighborArrayUnit_t     *partition_context_neighbor_array = picture_control_set_ptr->partition_context_neighbor_array;

    // CU Varaiables
    const BlockGeom          *blk_geom;
    CodingUnit_t             *cu_ptr;
    uint32_t                    cu_index = 0;
    uint32_t                    final_cu_index = 0;
    uint32_t                    cu_origin_x;
    uint32_t                    cu_origin_y;
    BlockSize                bsize;

    context_ptr->coded_area_sb = 0;
    context_ptr->coded_area_sb_uv = 0;

    tbPtr->quantized_coeffs_bits = 0;
    EbBool checkCuOutOfBound = EB_FALSE;

    SbGeom_t * sb_geom = &sequence_control_set_ptr->sb_geom[tbPtr->index];// .block_is_inside_md_scan[blk_index])

#if !TILES 
    if (context_ptr->sb_origin_x == 0 && context_ptr->sb_origin_y == 0)

        av1_reset_loop_restoration(picture_control_set_ptr);
#endif
    if (!(sb_geom->is_complete_sb)) {

        checkCuOutOfBound = EB_TRUE;
    }
    do {
        EbBool codeCuCond = EB_TRUE; // Code cu only if it is inside the picture

        cu_ptr = &tbPtr->final_cu_arr[final_cu_index];

        blk_geom = Get_blk_geom_mds(cu_index); // AMIR to be replaced with /*cu_ptr->mds_idx*/

        bsize = blk_geom->bsize;
        ASSERT(bsize < BlockSizeS_ALL);
        cu_origin_x = context_ptr->sb_origin_x + blk_geom->origin_x;
        cu_origin_y = context_ptr->sb_origin_y + blk_geom->origin_y;
        if (checkCuOutOfBound) {
            codeCuCond = (EbBool)sb_geom->block_is_inside_md_scan[cu_index]; // check if cu is inside the picture

            if ((cu_origin_x < sequence_control_set_ptr->luma_width) && (cu_origin_y < sequence_control_set_ptr->luma_height))
                codeCuCond = EB_TRUE;
        }


        if (codeCuCond) {
            uint32_t blkOriginX = cu_origin_x;
            uint32_t blkOriginY = cu_origin_y;

            const int32_t hbs = mi_size_wide[bsize] >> 1;
            const int32_t quarter_step = mi_size_wide[bsize] >> 2;
            Av1Common* cm = picture_control_set_ptr->parent_pcs_ptr->av1_cm;
            int32_t mi_row = blkOriginY >> MI_SIZE_LOG2;
            int32_t mi_col = blkOriginX >> MI_SIZE_LOG2;

            if (bsize >= BLOCK_8X8) {

                for (int32_t plane = 0; plane < 3; ++plane) {
                    int32_t rcol0, rcol1, rrow0, rrow1, tile_tl_idx;
                    if (av1_loop_restoration_corners_in_sb(cm, plane, mi_row, mi_col, bsize,
                        &rcol0, &rcol1, &rrow0, &rrow1,
                        &tile_tl_idx)) {
                        const int32_t rstride = cm->rst_info[plane].horz_units_per_tile;
                        for (int32_t rrow = rrow0; rrow < rrow1; ++rrow) {
                            for (int32_t rcol = rcol0; rcol < rcol1; ++rcol) {
                                const int32_t runit_idx = tile_tl_idx + rcol + rrow * rstride;
                                const RestorationUnitInfo *rui =
                                    &cm->rst_info[plane].unit_info[runit_idx];
                                loop_restoration_write_sb_coeffs(picture_control_set_ptr, frameContext, cm, /*xd,*/ rui, ecWriter, plane);
                            }
                        }
                    }
                }


                // Code Split Flag
                EncodePartitionAv1(
                    sequence_control_set_ptr,
                    frameContext,
                    ecWriter,
                    bsize,
                    tbPtr->cu_partition_array[cu_index],
                    blkOriginX,
                    blkOriginY,
                    partition_context_neighbor_array);

            }

            switch (tbPtr->cu_partition_array[cu_index]) {
            case PARTITION_NONE:
                write_modes_b(

                    picture_control_set_ptr,
                    context_ptr,
                    entropy_coder_ptr,
                    tbPtr,
                    cu_ptr,
                    coeffPtr);
                break;

            case PARTITION_HORZ:
                write_modes_b(
                    picture_control_set_ptr,
                    context_ptr,
                    entropy_coder_ptr,
                    tbPtr,
                    cu_ptr,
                    coeffPtr);

                if (mi_row + hbs < cm->mi_rows) {
                    final_cu_index++;
                    cu_ptr = &tbPtr->final_cu_arr[final_cu_index];
                    write_modes_b(
                        picture_control_set_ptr,
                        context_ptr,
                        entropy_coder_ptr,
                        tbPtr,
                        cu_ptr,
                        coeffPtr);
                }
                break;

            case PARTITION_VERT:
                write_modes_b(
                    picture_control_set_ptr,
                    context_ptr,
                    entropy_coder_ptr,
                    tbPtr,
                    cu_ptr,
                    coeffPtr);
                if (mi_col + hbs < cm->mi_cols) {
                    final_cu_index++;
                    cu_ptr = &tbPtr->final_cu_arr[final_cu_index];
                    write_modes_b(
                        picture_control_set_ptr,
                        context_ptr,
                        entropy_coder_ptr,
                        tbPtr,
                        cu_ptr,
                        coeffPtr);
                }
                break;
            case PARTITION_SPLIT:
                break;
            case PARTITION_HORZ_A:
                write_modes_b(
                    picture_control_set_ptr,
                    context_ptr,
                    entropy_coder_ptr,
                    tbPtr,
                    cu_ptr,
                    coeffPtr);

                final_cu_index++;
                cu_ptr = &tbPtr->final_cu_arr[final_cu_index];
                write_modes_b(
                    picture_control_set_ptr,
                    context_ptr,
                    entropy_coder_ptr,
                    tbPtr,
                    cu_ptr,
                    coeffPtr);

                final_cu_index++;
                cu_ptr = &tbPtr->final_cu_arr[final_cu_index];
                write_modes_b(
                    picture_control_set_ptr,
                    context_ptr,
                    entropy_coder_ptr,
                    tbPtr,
                    cu_ptr,
                    coeffPtr);

                break;
            case PARTITION_HORZ_B:
                write_modes_b(
                    picture_control_set_ptr,
                    context_ptr,
                    entropy_coder_ptr,
                    tbPtr,
                    cu_ptr,
                    coeffPtr);

                final_cu_index++;
                cu_ptr = &tbPtr->final_cu_arr[final_cu_index];
                write_modes_b(
                    picture_control_set_ptr,
                    context_ptr,
                    entropy_coder_ptr,
                    tbPtr,
                    cu_ptr,
                    coeffPtr);

                final_cu_index++;
                cu_ptr = &tbPtr->final_cu_arr[final_cu_index];
                write_modes_b(
                    picture_control_set_ptr,
                    context_ptr,
                    entropy_coder_ptr,
                    tbPtr,
                    cu_ptr,
                    coeffPtr);

                break;
            case PARTITION_VERT_A:
                write_modes_b(
                    picture_control_set_ptr,
                    context_ptr,
                    entropy_coder_ptr,
                    tbPtr,
                    cu_ptr,
                    coeffPtr);

                final_cu_index++;
                cu_ptr = &tbPtr->final_cu_arr[final_cu_index];
                write_modes_b(
                    picture_control_set_ptr,
                    context_ptr,
                    entropy_coder_ptr,
                    tbPtr,
                    cu_ptr,
                    coeffPtr);

                final_cu_index++;
                cu_ptr = &tbPtr->final_cu_arr[final_cu_index];
                write_modes_b(
                    picture_control_set_ptr,
                    context_ptr,
                    entropy_coder_ptr,
                    tbPtr,
                    cu_ptr,
                    coeffPtr);

                break;
            case PARTITION_VERT_B:
                write_modes_b(
                    picture_control_set_ptr,
                    context_ptr,
                    entropy_coder_ptr,
                    tbPtr,
                    cu_ptr,
                    coeffPtr);

                final_cu_index++;
                cu_ptr = &tbPtr->final_cu_arr[final_cu_index];
                write_modes_b(
                    picture_control_set_ptr,
                    context_ptr,
                    entropy_coder_ptr,
                    tbPtr,
                    cu_ptr,
                    coeffPtr);

                final_cu_index++;
                cu_ptr = &tbPtr->final_cu_arr[final_cu_index];
                write_modes_b(
                    picture_control_set_ptr,
                    context_ptr,
                    entropy_coder_ptr,
                    tbPtr,
                    cu_ptr,
                    coeffPtr);

                break;
            case PARTITION_HORZ_4:
                for (int32_t i = 0; i < 4; ++i) {
                    int32_t this_mi_row = mi_row + i * quarter_step;
                    if (i > 0 && this_mi_row >= cm->mi_rows) break;

                    if (i > 0) {
                        final_cu_index++;
                        cu_ptr = &tbPtr->final_cu_arr[final_cu_index];
                    }
                    write_modes_b(
                        picture_control_set_ptr,
                        context_ptr,
                        entropy_coder_ptr,
                        tbPtr,
                        cu_ptr,
                        coeffPtr);
                }
                break;
            case PARTITION_VERT_4:
                for (int32_t i = 0; i < 4; ++i) {
                    int32_t this_mi_col = mi_col + i * quarter_step;
                    if (i > 0 && this_mi_col >= cm->mi_cols) break;
                    if (i > 0) {
                        final_cu_index++;
                        cu_ptr = &tbPtr->final_cu_arr[final_cu_index];
                    }
                    write_modes_b(
                        picture_control_set_ptr,
                        context_ptr,
                        entropy_coder_ptr,
                        tbPtr,
                        cu_ptr,
                        coeffPtr);
                }
                break;
            default: assert(0);
            }

            if (tbPtr->cu_partition_array[cu_index] != PARTITION_SPLIT) {
                final_cu_index++;
                cu_index += ns_depth_offset[sequence_control_set_ptr->sb_size == BLOCK_128X128][blk_geom->depth];

            }
            else {
                cu_index += d1_depth_offset[sequence_control_set_ptr->sb_size == BLOCK_128X128][blk_geom->depth];

            }

        }
        else
            ++cu_index;
    } while (cu_index < sequence_control_set_ptr->max_block_cnt);
    return return_error;
}
