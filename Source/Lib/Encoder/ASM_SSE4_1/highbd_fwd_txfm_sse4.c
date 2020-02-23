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

#include <assert.h>
#include <smmintrin.h> /* SSE4.1 */

#include "EbDefinitions.h"
#include "aom_dsp_rtcd.h"
#include "emmintrin.h"
#include "EbTransforms.h"

static const int8_t *fwd_txfm_shift_ls[TX_SIZES_ALL] = {
    fwd_shift_4x4,   fwd_shift_8x8,   fwd_shift_16x16, fwd_shift_32x32, fwd_shift_64x64,
    fwd_shift_4x8,   fwd_shift_8x4,   fwd_shift_8x16,  fwd_shift_16x8,  fwd_shift_16x32,
    fwd_shift_32x16, fwd_shift_32x64, fwd_shift_64x32, fwd_shift_4x16,  fwd_shift_16x4,
    fwd_shift_8x32,  fwd_shift_32x8,  fwd_shift_16x64, fwd_shift_64x16,
};

static INLINE void load_buffer_4x4(const int16_t *input, __m128i *in, int32_t stride,
                                   int32_t flipud, int32_t fliplr, int32_t shift) {
    if (!flipud) {
        in[0] = _mm_loadl_epi64((const __m128i *)(input + 0 * stride));
        in[1] = _mm_loadl_epi64((const __m128i *)(input + 1 * stride));
        in[2] = _mm_loadl_epi64((const __m128i *)(input + 2 * stride));
        in[3] = _mm_loadl_epi64((const __m128i *)(input + 3 * stride));
    } else {
        in[0] = _mm_loadl_epi64((const __m128i *)(input + 3 * stride));
        in[1] = _mm_loadl_epi64((const __m128i *)(input + 2 * stride));
        in[2] = _mm_loadl_epi64((const __m128i *)(input + 1 * stride));
        in[3] = _mm_loadl_epi64((const __m128i *)(input + 0 * stride));
    }

    if (fliplr) {
        in[0] = _mm_shufflelo_epi16(in[0], 0x1b);
        in[1] = _mm_shufflelo_epi16(in[1], 0x1b);
        in[2] = _mm_shufflelo_epi16(in[2], 0x1b);
        in[3] = _mm_shufflelo_epi16(in[3], 0x1b);
    }

    in[0] = _mm_cvtepi16_epi32(in[0]);
    in[1] = _mm_cvtepi16_epi32(in[1]);
    in[2] = _mm_cvtepi16_epi32(in[2]);
    in[3] = _mm_cvtepi16_epi32(in[3]);

    in[0] = _mm_slli_epi32(in[0], shift);
    in[1] = _mm_slli_epi32(in[1], shift);
    in[2] = _mm_slli_epi32(in[2], shift);
    in[3] = _mm_slli_epi32(in[3], shift);
}

static void fidtx4x4_sse4_1(__m128i *in, __m128i *out, int32_t bit, int32_t col_num) {
    (void)bit;
    __m128i fact   = _mm_set1_epi32(new_sqrt2);
    __m128i offset = _mm_set1_epi32(1 << (new_sqrt2_bits - 1));
    __m128i a_low;
    __m128i v[4];

    for (int32_t i = 0; i < 4; i++) {
        a_low  = _mm_mullo_epi32(in[i * col_num], fact);
        a_low  = _mm_add_epi32(a_low, offset);
        out[i] = _mm_srai_epi32(a_low, new_sqrt2_bits);
    }

    // Transpose for 4x4
    v[0] = _mm_unpacklo_epi32(out[0], out[1]);
    v[1] = _mm_unpackhi_epi32(out[0], out[1]);
    v[2] = _mm_unpacklo_epi32(out[2], out[3]);
    v[3] = _mm_unpackhi_epi32(out[2], out[3]);

    out[0] = _mm_unpacklo_epi64(v[0], v[2]);
    out[1] = _mm_unpackhi_epi64(v[0], v[2]);
    out[2] = _mm_unpacklo_epi64(v[1], v[3]);
    out[3] = _mm_unpackhi_epi64(v[1], v[3]);
}

// We only use stage-2 bit;
// shift[0] is used in load_buffer_4x4()
// shift[1] is used in txfm_func_col()
// shift[2] is used in txfm_func_row()
static void fdct4x4_sse4_1(__m128i *in, __m128i *out, int32_t bit, const int32_t num_col) {
    const int32_t *cospi   = cospi_arr(bit);
    const __m128i  cospi32 = _mm_set1_epi32(cospi[32]);
    const __m128i  cospi48 = _mm_set1_epi32(cospi[48]);
    const __m128i  cospi16 = _mm_set1_epi32(cospi[16]);
    const __m128i  rnding  = _mm_set1_epi32(1 << (bit - 1));
    __m128i        s0, s1, s2, s3;
    __m128i        u0, u1, u2, u3;
    __m128i        v0, v1, v2, v3;

    int32_t endidx = 3 * num_col;
    s0             = _mm_add_epi32(in[0], in[endidx]);
    s3             = _mm_sub_epi32(in[0], in[endidx]);
    endidx -= num_col;
    s1 = _mm_add_epi32(in[num_col], in[endidx]);
    s2 = _mm_sub_epi32(in[num_col], in[endidx]);

    // btf_32_sse4_1_type0(cospi32, cospi32, s[01], u[02], bit);
    u0 = _mm_mullo_epi32(s0, cospi32);
    u1 = _mm_mullo_epi32(s1, cospi32);
    u2 = _mm_add_epi32(u0, u1);
    v0 = _mm_sub_epi32(u0, u1);

    u3 = _mm_add_epi32(u2, rnding);
    v1 = _mm_add_epi32(v0, rnding);

    u0 = _mm_srai_epi32(u3, bit);
    u2 = _mm_srai_epi32(v1, bit);

    // btf_32_sse4_1_type1(cospi48, cospi16, s[23], u[13], bit);
    v0 = _mm_mullo_epi32(s2, cospi48);
    v1 = _mm_mullo_epi32(s3, cospi16);
    v2 = _mm_add_epi32(v0, v1);

    v3 = _mm_add_epi32(v2, rnding);
    u1 = _mm_srai_epi32(v3, bit);

    v0 = _mm_mullo_epi32(s2, cospi16);
    v1 = _mm_mullo_epi32(s3, cospi48);
    v2 = _mm_sub_epi32(v1, v0);

    v3 = _mm_add_epi32(v2, rnding);
    u3 = _mm_srai_epi32(v3, bit);

    // Note: shift[1] and shift[2] are zeros

    // Transpose 4x4 32-bit
    v0 = _mm_unpacklo_epi32(u0, u1);
    v1 = _mm_unpackhi_epi32(u0, u1);
    v2 = _mm_unpacklo_epi32(u2, u3);
    v3 = _mm_unpackhi_epi32(u2, u3);

    out[0] = _mm_unpacklo_epi64(v0, v2);
    out[1] = _mm_unpackhi_epi64(v0, v2);
    out[2] = _mm_unpacklo_epi64(v1, v3);
    out[3] = _mm_unpackhi_epi64(v1, v3);
}

static INLINE void write_buffer_4x4(__m128i *res, int32_t *output) {
    _mm_storeu_si128((__m128i *)(output + 0 * 4), res[0]);
    _mm_storeu_si128((__m128i *)(output + 1 * 4), res[1]);
    _mm_storeu_si128((__m128i *)(output + 2 * 4), res[2]);
    _mm_storeu_si128((__m128i *)(output + 3 * 4), res[3]);
}

static void fadst4x4_sse4_1(__m128i *in, __m128i *out, int32_t bit, const int32_t num_col) {
    const int32_t *sinpi  = sinpi_arr(bit);
    const __m128i  rnding = _mm_set1_epi32(1 << (bit - 1));
    const __m128i  sinpi1 = _mm_set1_epi32((int32_t)sinpi[1]);
    const __m128i  sinpi2 = _mm_set1_epi32((int32_t)sinpi[2]);
    const __m128i  sinpi3 = _mm_set1_epi32((int32_t)sinpi[3]);
    const __m128i  sinpi4 = _mm_set1_epi32((int32_t)sinpi[4]);
    __m128i        t;
    __m128i        s0, s1, s2, s3, s4, s5, s6, s7;
    __m128i        x0, x1, x2, x3;
    __m128i        u0, u1, u2, u3;
    __m128i        v0, v1, v2, v3;

    int32_t idx = 0 * num_col;
    s0          = _mm_mullo_epi32(in[idx], sinpi1);
    s1          = _mm_mullo_epi32(in[idx], sinpi4);
    t           = _mm_add_epi32(in[idx], in[idx + num_col]);
    idx += num_col;
    s2 = _mm_mullo_epi32(in[idx], sinpi2);
    s3 = _mm_mullo_epi32(in[idx], sinpi1);
    idx += num_col;
    s4 = _mm_mullo_epi32(in[idx], sinpi3);
    idx += num_col;
    s5 = _mm_mullo_epi32(in[idx], sinpi4);
    s6 = _mm_mullo_epi32(in[idx], sinpi2);
    s7 = _mm_sub_epi32(t, in[idx]);

    t  = _mm_add_epi32(s0, s2);
    x0 = _mm_add_epi32(t, s5);
    x1 = _mm_mullo_epi32(s7, sinpi3);
    t  = _mm_sub_epi32(s1, s3);
    x2 = _mm_add_epi32(t, s6);
    x3 = s4;

    s0 = _mm_add_epi32(x0, x3);
    s1 = x1;
    s2 = _mm_sub_epi32(x2, x3);
    t  = _mm_sub_epi32(x2, x0);
    s3 = _mm_add_epi32(t, x3);

    u0 = _mm_add_epi32(s0, rnding);
    u0 = _mm_srai_epi32(u0, bit);

    u1 = _mm_add_epi32(s1, rnding);
    u1 = _mm_srai_epi32(u1, bit);

    u2 = _mm_add_epi32(s2, rnding);
    u2 = _mm_srai_epi32(u2, bit);

    u3 = _mm_add_epi32(s3, rnding);
    u3 = _mm_srai_epi32(u3, bit);

    v0 = _mm_unpacklo_epi32(u0, u1);
    v1 = _mm_unpackhi_epi32(u0, u1);
    v2 = _mm_unpacklo_epi32(u2, u3);
    v3 = _mm_unpackhi_epi32(u2, u3);

    out[0] = _mm_unpacklo_epi64(v0, v2);
    out[1] = _mm_unpackhi_epi64(v0, v2);
    out[2] = _mm_unpacklo_epi64(v1, v3);
    out[3] = _mm_unpackhi_epi64(v1, v3);
}

void eb_av1_fwd_txfm2d_4x4_sse4_1(int16_t *input, int32_t *coeff, uint32_t stride, TxType tx_type,
                                  uint8_t bd) {
    __m128i       in[4];
    const int8_t *shift   = fwd_txfm_shift_ls[TX_4X4];
    const int32_t txw_idx = get_txw_idx(TX_4X4);
    const int32_t txh_idx = get_txh_idx(TX_4X4);

    switch (tx_type) {
    case DCT_DCT:
        load_buffer_4x4(input, in, stride, 0, 0, shift[0]);
        fdct4x4_sse4_1(in, in, fwd_cos_bit_col[txw_idx][txh_idx], 1);
        fdct4x4_sse4_1(in, in, fwd_cos_bit_row[txw_idx][txh_idx], 1);
        write_buffer_4x4(in, coeff);
        break;
    case ADST_DCT:
        load_buffer_4x4(input, in, stride, 0, 0, shift[0]);
        fadst4x4_sse4_1(in, in, fwd_cos_bit_col[txw_idx][txh_idx], 1);
        fdct4x4_sse4_1(in, in, fwd_cos_bit_row[txw_idx][txh_idx], 1);
        write_buffer_4x4(in, coeff);
        break;
    case DCT_ADST:
        load_buffer_4x4(input, in, stride, 0, 0, shift[0]);
        fdct4x4_sse4_1(in, in, fwd_cos_bit_col[txw_idx][txh_idx], 1);
        fadst4x4_sse4_1(in, in, fwd_cos_bit_row[txw_idx][txh_idx], 1);
        write_buffer_4x4(in, coeff);
        break;
    case ADST_ADST:
        load_buffer_4x4(input, in, stride, 0, 0, shift[0]);
        fadst4x4_sse4_1(in, in, fwd_cos_bit_col[txw_idx][txh_idx], 1);
        fadst4x4_sse4_1(in, in, fwd_cos_bit_row[txw_idx][txh_idx], 1);
        write_buffer_4x4(in, coeff);
        break;
    case FLIPADST_DCT:
        load_buffer_4x4(input, in, stride, 1, 0, shift[0]);
        fadst4x4_sse4_1(in, in, fwd_cos_bit_col[txw_idx][txh_idx], 1);
        fdct4x4_sse4_1(in, in, fwd_cos_bit_row[txw_idx][txh_idx], 1);
        write_buffer_4x4(in, coeff);
        break;
    case DCT_FLIPADST:
        load_buffer_4x4(input, in, stride, 0, 1, shift[0]);
        fdct4x4_sse4_1(in, in, fwd_cos_bit_col[txw_idx][txh_idx], 1);
        fadst4x4_sse4_1(in, in, fwd_cos_bit_row[txw_idx][txh_idx], 1);
        write_buffer_4x4(in, coeff);
        break;
    case FLIPADST_FLIPADST:
        load_buffer_4x4(input, in, stride, 1, 1, shift[0]);
        fadst4x4_sse4_1(in, in, fwd_cos_bit_col[txw_idx][txh_idx], 1);
        fadst4x4_sse4_1(in, in, fwd_cos_bit_row[txw_idx][txh_idx], 1);
        write_buffer_4x4(in, coeff);
        break;
    case ADST_FLIPADST:
        load_buffer_4x4(input, in, stride, 0, 1, shift[0]);
        fadst4x4_sse4_1(in, in, fwd_cos_bit_col[txw_idx][txh_idx], 1);
        fadst4x4_sse4_1(in, in, fwd_cos_bit_row[txw_idx][txh_idx], 1);
        write_buffer_4x4(in, coeff);
        break;
    case FLIPADST_ADST:
        load_buffer_4x4(input, in, stride, 1, 0, shift[0]);
        fadst4x4_sse4_1(in, in, fwd_cos_bit_col[txw_idx][txh_idx], 1);
        fadst4x4_sse4_1(in, in, fwd_cos_bit_row[txw_idx][txh_idx], 1);
        write_buffer_4x4(in, coeff);
        break;
    case IDTX:
        load_buffer_4x4(input, in, stride, 0, 0, shift[0]);
        fidtx4x4_sse4_1(in, in, fwd_cos_bit_col[txw_idx][txh_idx], 1);
        fidtx4x4_sse4_1(in, in, fwd_cos_bit_row[txw_idx][txh_idx], 1);
        write_buffer_4x4(in, coeff);
        break;
    case V_DCT:
        load_buffer_4x4(input, in, stride, 0, 0, shift[0]);
        fdct4x4_sse4_1(in, in, fwd_cos_bit_col[txw_idx][txh_idx], 1);
        fidtx4x4_sse4_1(in, in, fwd_cos_bit_row[txw_idx][txh_idx], 1);
        write_buffer_4x4(in, coeff);
        break;
    case H_DCT:
        load_buffer_4x4(input, in, stride, 0, 0, shift[0]);
        fidtx4x4_sse4_1(in, in, fwd_cos_bit_row[txw_idx][txh_idx], 1);
        fdct4x4_sse4_1(in, in, fwd_cos_bit_col[txw_idx][txh_idx], 1);
        write_buffer_4x4(in, coeff);
        break;
    case V_ADST:
        load_buffer_4x4(input, in, stride, 0, 0, shift[0]);
        fadst4x4_sse4_1(in, in, fwd_cos_bit_col[txw_idx][txh_idx], 1);
        fidtx4x4_sse4_1(in, in, fwd_cos_bit_row[txw_idx][txh_idx], 1);
        write_buffer_4x4(in, coeff);
        break;
    case H_ADST:
        load_buffer_4x4(input, in, stride, 0, 0, shift[0]);
        fidtx4x4_sse4_1(in, in, fwd_cos_bit_row[txw_idx][txh_idx], 1);
        fadst4x4_sse4_1(in, in, fwd_cos_bit_col[txw_idx][txh_idx], 1);
        write_buffer_4x4(in, coeff);
        break;
    case V_FLIPADST:
        load_buffer_4x4(input, in, stride, 1, 0, shift[0]);
        fadst4x4_sse4_1(in, in, fwd_cos_bit_row[txw_idx][txh_idx], 1);
        fidtx4x4_sse4_1(in, in, fwd_cos_bit_row[txw_idx][txh_idx], 1);
        write_buffer_4x4(in, coeff);
        break;
    case H_FLIPADST:
        load_buffer_4x4(input, in, stride, 0, 1, shift[0]);
        fidtx4x4_sse4_1(in, in, fwd_cos_bit_row[txw_idx][txh_idx], 1);
        fadst4x4_sse4_1(in, in, fwd_cos_bit_row[txw_idx][txh_idx], 1);
        write_buffer_4x4(in, coeff);
        break;
    default: assert(0);
    }
    (void)bd;
}
