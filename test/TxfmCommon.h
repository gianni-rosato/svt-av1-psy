/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/******************************************************************************
 * @file TxfmCommon.h
 *
 * @brief declaration of fwd/inv txfm functions, including:
 * - 1d fwd/inv txfm function
 * - map type and size to function;
 *
 * @author Cidana-Edmond, Cidana-Wenyao
 *
 ******************************************************************************/
#ifndef _TXFM_TEST_H_
#define _TXFM_TEST_H_

#include "aom_dsp_rtcd.h"

#ifdef __cplusplus
extern "C" {
#endif

static const int8_t test_txfm_range[12] = {
    20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20};

// export forward transform functions
void eb_av1_fdct4_new(const int32_t *input, int32_t *output, int8_t cos_bit,
                   const int8_t *stage_range);
void eb_av1_fdct8_new(const int32_t *input, int32_t *output, int8_t cos_bit,
                   const int8_t *stage_range);
void eb_av1_fdct16_new(const int32_t *input, int32_t *output, int8_t cos_bit,
                    const int8_t *stage_range);
void eb_av1_fdct32_new(const int32_t *input, int32_t *output, int8_t cos_bit,
                    const int8_t *stage_range);
void eb_av1_fdct64_new(const int32_t *input, int32_t *output, int8_t cos_bit,
                    const int8_t *stage_range);
void eb_av1_fadst4_new(const int32_t *input, int32_t *output, int8_t cos_bit,
                    const int8_t *stage_range);
void eb_av1_fadst8_new(const int32_t *input, int32_t *output, int8_t cos_bit,
                    const int8_t *stage_range);
void eb_av1_fadst16_new(const int32_t *input, int32_t *output, int8_t cos_bit,
                     const int8_t *stage_range);
void av1_fadst32_new(const int32_t *input, int32_t *output, int8_t cos_bit,
                     const int8_t *stage_range);
void eb_av1_fidentity4_c(const int32_t *input, int32_t *output, int8_t cos_bit,
                      const int8_t *stage_range);
void eb_av1_fidentity8_c(const int32_t *input, int32_t *output, int8_t cos_bit,
                      const int8_t *stage_range);
void eb_av1_fidentity16_c(const int32_t *input, int32_t *output, int8_t cos_bit,
                       const int8_t *stage_range);
void eb_av1_fidentity32_c(const int32_t *input, int32_t *output, int8_t cos_bit,
                       const int8_t *stage_range);
void av1_fidentity64_c(const int32_t *input, int32_t *output, int8_t cos_bit,
                       const int8_t *stage_range);

// export inverse transform functions
void eb_av1_idct4_new(const int32_t *input, int32_t *output, int8_t cos_bit,
                   const int8_t *stage_range);
void eb_av1_idct8_new(const int32_t *input, int32_t *output, int8_t cos_bit,
                   const int8_t *stage_range);
void eb_av1_idct16_new(const int32_t *input, int32_t *output, int8_t cos_bit,
                    const int8_t *stage_range);
void eb_av1_idct32_new(const int32_t *input, int32_t *output, int8_t cos_bit,
                    const int8_t *stage_range);
void eb_av1_idct64_new(const int32_t *input, int32_t *output, int8_t cos_bit,
                    const int8_t *stage_range);
void eb_av1_iadst4_new(const int32_t *input, int32_t *output, int8_t cos_bit,
                    const int8_t *stage_range);
void eb_av1_iadst8_new(const int32_t *input, int32_t *output, int8_t cos_bit,
                    const int8_t *stage_range);
void eb_av1_iadst16_new(const int32_t *input, int32_t *output, int8_t cos_bit,
                     const int8_t *stage_range);
void av1_iadst32_new(const int32_t *input, int32_t *output, int8_t cos_bit,
                     const int8_t *stage_range);
void eb_av1_iidentity4_c(const int32_t *input, int32_t *output, int8_t cos_bit,
                      const int8_t *stage_range);
void eb_av1_iidentity8_c(const int32_t *input, int32_t *output, int8_t cos_bit,
                      const int8_t *stage_range);
void eb_av1_iidentity16_c(const int32_t *input, int32_t *output, int8_t cos_bit,
                       const int8_t *stage_range);
void eb_av1_iidentity32_c(const int32_t *input, int32_t *output, int8_t cos_bit,
                       const int8_t *stage_range);
void av1_iidentity64_c(const int32_t *input, int32_t *output, int8_t cos_bit,
                       const int8_t *stage_range);

void Av1TransformConfig(TxType tx_type, TxSize tx_size, Txfm2DFlipCfg *cfg);

typedef void (*Txfm1dFunc)(const int32_t *input, int32_t *output,
                           int8_t cos_bit, const int8_t *stage_range);

typedef void (*TxfmFwd2dFunc)(int16_t *input, int32_t *output,
                              uint32_t input_stride, TxType transform_type,
                              uint8_t bit_depth);

static INLINE Txfm1dFunc fwd_txfm_type_to_func(TxfmType txfm_type) {
    switch (txfm_type) {
    case TXFM_TYPE_DCT4: return eb_av1_fdct4_new;
    case TXFM_TYPE_DCT8: return eb_av1_fdct8_new;
    case TXFM_TYPE_DCT16: return eb_av1_fdct16_new;
    case TXFM_TYPE_DCT32: return eb_av1_fdct32_new;
    case TXFM_TYPE_DCT64: return eb_av1_fdct64_new;
    case TXFM_TYPE_ADST4: return eb_av1_fadst4_new;
    case TXFM_TYPE_ADST8: return eb_av1_fadst8_new;
    case TXFM_TYPE_ADST16: return eb_av1_fadst16_new;
    case TXFM_TYPE_ADST32: return av1_fadst32_new;
    case TXFM_TYPE_IDENTITY4: return eb_av1_fidentity4_c;
    case TXFM_TYPE_IDENTITY8: return eb_av1_fidentity8_c;
    case TXFM_TYPE_IDENTITY16: return eb_av1_fidentity16_c;
    case TXFM_TYPE_IDENTITY32: return eb_av1_fidentity32_c;
    case TXFM_TYPE_IDENTITY64: return av1_fidentity64_c;
    default: assert(0); return NULL;
    }
}

static INLINE Txfm1dFunc inv_txfm_type_to_func(TxfmType txfm_type) {
    switch (txfm_type) {
    case TXFM_TYPE_DCT4: return eb_av1_idct4_new;
    case TXFM_TYPE_DCT8: return eb_av1_idct8_new;
    case TXFM_TYPE_DCT16: return eb_av1_idct16_new;
    case TXFM_TYPE_DCT32: return eb_av1_idct32_new;
    case TXFM_TYPE_DCT64: return eb_av1_idct64_new;
    case TXFM_TYPE_ADST4: return eb_av1_iadst4_new;
    case TXFM_TYPE_ADST8: return eb_av1_iadst8_new;
    case TXFM_TYPE_ADST16: return eb_av1_iadst16_new;
    case TXFM_TYPE_ADST32: return av1_iadst32_new;
    case TXFM_TYPE_IDENTITY4: return eb_av1_iidentity4_c;
    case TXFM_TYPE_IDENTITY8: return eb_av1_iidentity8_c;
    case TXFM_TYPE_IDENTITY16: return eb_av1_iidentity16_c;
    case TXFM_TYPE_IDENTITY32: return eb_av1_iidentity32_c;
    case TXFM_TYPE_IDENTITY64: return av1_iidentity64_c;
    default: assert(0); return NULL;
    }
}

using FwdTxfm2dFunc = void (*)(int16_t *input, int32_t *output, uint32_t stride,
                               TxType tx_type, uint8_t bd);
static const FwdTxfm2dFunc fwd_txfm_2d_c_func[TX_SIZES_ALL] = {
    Av1TransformTwoD_4x4_c,   Av1TransformTwoD_8x8_c,
    Av1TransformTwoD_16x16_c, Av1TransformTwoD_32x32_c,
    Av1TransformTwoD_64x64_c, eb_av1_fwd_txfm2d_4x8_c,
    eb_av1_fwd_txfm2d_8x4_c,     eb_av1_fwd_txfm2d_8x16_c,
    eb_av1_fwd_txfm2d_16x8_c,    eb_av1_fwd_txfm2d_16x32_c,
    eb_av1_fwd_txfm2d_32x16_c,   eb_av1_fwd_txfm2d_32x64_c,
    eb_av1_fwd_txfm2d_64x32_c,   eb_av1_fwd_txfm2d_4x16_c,
    eb_av1_fwd_txfm2d_16x4_c,    eb_av1_fwd_txfm2d_8x32_c,
    eb_av1_fwd_txfm2d_32x8_c,    eb_av1_fwd_txfm2d_16x64_c,
    eb_av1_fwd_txfm2d_64x16_c,
};

static INLINE int get_txfm1d_size(TxfmType txfm_type) {
    switch (txfm_type) {
    case TXFM_TYPE_DCT4:
    case TXFM_TYPE_ADST4:
    case TXFM_TYPE_IDENTITY4: return 4;
    case TXFM_TYPE_DCT8:
    case TXFM_TYPE_ADST8:
    case TXFM_TYPE_IDENTITY8: return 8;
    case TXFM_TYPE_DCT16:
    case TXFM_TYPE_ADST16:
    case TXFM_TYPE_IDENTITY16: return 16;
    case TXFM_TYPE_DCT32:
    case TXFM_TYPE_ADST32:
    case TXFM_TYPE_IDENTITY32: return 32;
    case TXFM_TYPE_DCT64:
    case TXFM_TYPE_IDENTITY64: return 64;
    default: assert(0); return 0;
    }
}

static INLINE bool is_txfm_allowed(TxType tx_type, TxSize tx_size) {
    /* According to 5.11 */
    const TxSize sqr_up = txsize_sqr_up_map[tx_size];
    if (sqr_up > TX_32X32 && tx_type != DCT_DCT)
        return false;
    else if (sqr_up == TX_32X32 && (tx_type != DCT_DCT && tx_type != IDTX))
        return false;
    else /* <=TX_16X16, all type are supported */
        return true;
}

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

using IsTxTypeImpFunc = bool (*)(const TxType tx_type);
static INLINE bool all_txtype_imp(const TxType tx_type) {
    (void)tx_type;
    return true;
}

static INLINE bool dct_adst_combine_imp(const TxType tx_type) {
    switch (tx_type) {
    case DCT_DCT:
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
    case FLIPADST_DCT:
    case DCT_FLIPADST:
    case FLIPADST_FLIPADST:
    case ADST_FLIPADST:
    case FLIPADST_ADST: return true;
    default: return false;
    }
}

#ifdef __cplusplus
}
#endif

#endif  // _TXFM_TEST_H_
