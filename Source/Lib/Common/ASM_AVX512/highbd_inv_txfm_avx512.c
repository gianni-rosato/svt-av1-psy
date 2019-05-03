#include <assert.h>
#include <immintrin.h>
#include "EbDefinitions.h"
#include "aom_dsp_rtcd.h"
#include "EbTransforms.h"

extern const int8_t *eb_inv_txfm_shift_ls[];
const int32_t *cospi_arr(int32_t n);

#define TRANSPOSE_4X4_AVX512(x0, x1, x2, x3, y0, y1, y2, y3) \
    do {                                                    \
    __m512i u0, u1, u2, u3;                             \
    u0 = _mm512_unpacklo_epi32(x0, x1);                    \
    u1 = _mm512_unpackhi_epi32(x0, x1);                    \
    u2 = _mm512_unpacklo_epi32(x2, x3);                    \
    u3 = _mm512_unpackhi_epi32(x2, x3);                    \
    y0 = _mm512_unpacklo_epi64(u0, u2);                    \
    y1 = _mm512_unpackhi_epi64(u0, u2);                    \
    y2 = _mm512_unpacklo_epi64(u1, u3);                    \
    y3 = _mm512_unpackhi_epi64(u1, u3);                    \
    } while (0)


static INLINE void transpose_16x16_avx512(int32_t stride, const __m512i *in, __m512i *out) {
    __m512i out1[16];
    TRANSPOSE_4X4_AVX512(in[0 * stride], in[1 * stride], in[2 * stride], in[3 * stride], out1[0], out1[1], out1[2], out1[3]);
    TRANSPOSE_4X4_AVX512(in[4 * stride], in[5 * stride], in[6 * stride], in[7 * stride], out1[4], out1[5], out1[6], out1[7]);
    TRANSPOSE_4X4_AVX512(in[8 * stride], in[9 * stride], in[10 * stride], in[11 * stride], out1[8], out1[9], out1[10], out1[11]);
    TRANSPOSE_4X4_AVX512(in[12 * stride], in[13 * stride], in[14 * stride], in[15 * stride], out1[12], out1[13], out1[14], out1[15]);

    __m128i * outptr = (__m128i *)(out + 0 * stride);

    //will get first row of transpose matrix from corresponding 4 vectors in out1
    outptr[0] = _mm512_extracti32x4_epi32(out1[0], 0);
    outptr[1] = _mm512_extracti32x4_epi32(out1[4], 0);
    outptr[2] = _mm512_extracti32x4_epi32(out1[8], 0);
    outptr[3] = _mm512_extracti32x4_epi32(out1[12], 0);

    //will get second row of transpose matrix from corresponding 4 vectors in out1
    outptr = (__m128i *)(out + 1 * stride);
    outptr[0] = _mm512_extracti32x4_epi32(out1[1], 0);
    outptr[1] = _mm512_extracti32x4_epi32(out1[5], 0);
    outptr[2] = _mm512_extracti32x4_epi32(out1[9], 0);
    outptr[3] = _mm512_extracti32x4_epi32(out1[13], 0);

    //will get 3rd row of transpose matrix from corresponding 4 vectors in out1
    outptr = (__m128i *)(out + 2 * stride);
    outptr[0] = _mm512_extracti32x4_epi32(out1[2], 0);
    outptr[1] = _mm512_extracti32x4_epi32(out1[6], 0);
    outptr[2] = _mm512_extracti32x4_epi32(out1[10], 0);
    outptr[3] = _mm512_extracti32x4_epi32(out1[14], 0);

    //will get 4th row of transpose matrix from corresponding 4 vectors in out1
    outptr = (__m128i *)(out + 3 * stride);
    outptr[0] = _mm512_extracti32x4_epi32(out1[3], 0);
    outptr[1] = _mm512_extracti32x4_epi32(out1[7], 0);
    outptr[2] = _mm512_extracti32x4_epi32(out1[11], 0);
    outptr[3] = _mm512_extracti32x4_epi32(out1[15], 0);

    //will get 5th row of transpose matrix from corresponding 4 vectors in out1
    outptr = (__m128i *)(out + 4 * stride);
    outptr[0] = _mm512_extracti32x4_epi32(out1[0], 1);
    outptr[1] = _mm512_extracti32x4_epi32(out1[4], 1);
    outptr[2] = _mm512_extracti32x4_epi32(out1[8], 1);
    outptr[3] = _mm512_extracti32x4_epi32(out1[12], 1);

    //will get 6th row of transpose matrix from corresponding 4 vectors in out1
    outptr = (__m128i *)(out + 5 * stride);
    outptr[0] = _mm512_extracti32x4_epi32(out1[1], 1);
    outptr[1] = _mm512_extracti32x4_epi32(out1[5], 1);
    outptr[2] = _mm512_extracti32x4_epi32(out1[9], 1);
    outptr[3] = _mm512_extracti32x4_epi32(out1[13], 1);

    //will get 7th row of transpose matrix from corresponding 4 vectors in out1
    outptr = (__m128i *)(out + 6 * stride);
    outptr[0] = _mm512_extracti32x4_epi32(out1[2], 1);
    outptr[1] = _mm512_extracti32x4_epi32(out1[6], 1);
    outptr[2] = _mm512_extracti32x4_epi32(out1[10], 1);
    outptr[3] = _mm512_extracti32x4_epi32(out1[14], 1);

    //will get 8th row of transpose matrix from corresponding 4 vectors in out1
    outptr = (__m128i *)(out + 7 * stride);
    outptr[0] = _mm512_extracti32x4_epi32(out1[3], 1);
    outptr[1] = _mm512_extracti32x4_epi32(out1[7], 1);
    outptr[2] = _mm512_extracti32x4_epi32(out1[11], 1);
    outptr[3] = _mm512_extracti32x4_epi32(out1[15], 1);

    //will get 9th row of transpose matrix from corresponding 4 vectors in out1
    outptr = (__m128i *)(out + 8 * stride);
    outptr[0] = _mm512_extracti32x4_epi32(out1[0], 2);
    outptr[1] = _mm512_extracti32x4_epi32(out1[4], 2);
    outptr[2] = _mm512_extracti32x4_epi32(out1[8], 2);
    outptr[3] = _mm512_extracti32x4_epi32(out1[12], 2);

    //will get 10th row of transpose matrix from corresponding 4 vectors in out1
    outptr = (__m128i *)(out + 9 * stride);
    outptr[0] = _mm512_extracti32x4_epi32(out1[1], 2);
    outptr[1] = _mm512_extracti32x4_epi32(out1[5], 2);
    outptr[2] = _mm512_extracti32x4_epi32(out1[9], 2);
    outptr[3] = _mm512_extracti32x4_epi32(out1[13], 2);

    //will get 11th row of transpose matrix from corresponding 4 vectors in out1
    outptr = (__m128i *)(out + 10 * stride);
    outptr[0] = _mm512_extracti32x4_epi32(out1[2], 2);
    outptr[1] = _mm512_extracti32x4_epi32(out1[6], 2);
    outptr[2] = _mm512_extracti32x4_epi32(out1[10], 2);
    outptr[3] = _mm512_extracti32x4_epi32(out1[14], 2);

    //will get 12th row of transpose matrix from corresponding 4 vectors in out1
    outptr = (__m128i *)(out + 11 * stride);
    outptr[0] = _mm512_extracti32x4_epi32(out1[3], 2);
    outptr[1] = _mm512_extracti32x4_epi32(out1[7], 2);
    outptr[2] = _mm512_extracti32x4_epi32(out1[11], 2);
    outptr[3] = _mm512_extracti32x4_epi32(out1[15], 2);

    //will get 13th row of transpose matrix from corresponding 4 vectors in out1
    outptr = (__m128i *)(out + 12 * stride);
    outptr[0] = _mm512_extracti32x4_epi32(out1[0], 3);
    outptr[1] = _mm512_extracti32x4_epi32(out1[4], 3);
    outptr[2] = _mm512_extracti32x4_epi32(out1[8], 3);
    outptr[3] = _mm512_extracti32x4_epi32(out1[12], 3);

    //will get 14th row of transpose matrix from corresponding 4 vectors in out1
    outptr = (__m128i *)(out + 13 * stride);
    outptr[0] = _mm512_extracti32x4_epi32(out1[1], 3);
    outptr[1] = _mm512_extracti32x4_epi32(out1[5], 3);
    outptr[2] = _mm512_extracti32x4_epi32(out1[9], 3);
    outptr[3] = _mm512_extracti32x4_epi32(out1[13], 3);

    //will get 15th row of transpose matrix from corresponding 4 vectors in out1
    outptr = (__m128i *)(out + 14 * stride);
    outptr[0] = _mm512_extracti32x4_epi32(out1[2], 3);
    outptr[1] = _mm512_extracti32x4_epi32(out1[6], 3);
    outptr[2] = _mm512_extracti32x4_epi32(out1[10], 3);
    outptr[3] = _mm512_extracti32x4_epi32(out1[14], 3);

    //will get 16th row of transpose matrix from corresponding 4 vectors in out1
    outptr = (__m128i *)(out + 15 * stride);
    outptr[0] = _mm512_extracti32x4_epi32(out1[3], 3);
    outptr[1] = _mm512_extracti32x4_epi32(out1[7], 3);
    outptr[2] = _mm512_extracti32x4_epi32(out1[11], 3);
    outptr[3] = _mm512_extracti32x4_epi32(out1[15], 3);
}

static INLINE __m512i half_btf_avx512(const __m512i *w0, const __m512i *n0,
    const __m512i *w1, const __m512i *n1,
    const __m512i *rounding, int32_t bit) {
    __m512i x, y;

    x = _mm512_mullo_epi32(*w0, *n0);
    y = _mm512_mullo_epi32(*w1, *n1);
    x = _mm512_add_epi32(x, y);
    x = _mm512_add_epi32(x, *rounding);
    x = _mm512_srai_epi32(x, (uint8_t)bit);
    return x;
}

static INLINE void highbd_clamp_epi32_avx512(__m512i *x, int32_t bd) {
    const __m512i zero = _mm512_setzero_si512();
    const __m512i max = _mm512_set1_epi32((1 << bd) - 1);

    *x = _mm512_min_epi32(*x, max);
    *x = _mm512_max_epi32(*x, zero);
}


static void load_buffer_16x16_avx512(const int32_t *coeff, __m512i *in) {
    int32_t i;
    for (i = 0; i < 16; ++i) {
        in[i] = _mm512_load_si512((const __m512i *)coeff);
        coeff += 16;
    }
}

static INLINE void write_buffer_16x16_avx512(__m512i *in, uint16_t *output,
    int32_t stride, int32_t fliplr,
    int32_t flipud, int32_t bd) {
    __m512i u1, v0, index;
    __m256i u0;
    int32_t i = 0;
    uint32_t idx[] = { 15 , 14 , 13 , 12 , 11 , 10 , 9 , 8 , 7 , 6 , 5 , 4 , 3 , 2 , 1 , 0 };
    index = _mm512_loadu_si512((const __m256i *)idx);

    if (flipud) {
        output += stride * 15;
        stride = -stride;
    }

    while (i < 16) {
        u0 = _mm256_loadu_si256((const __m256i *)output);

        u1 = _mm512_cvtepi16_epi32(u0);
        v0 = in[i];
        if (fliplr) {
        v0 = _mm512_permutexvar_epi32(index,v0);
        }

        v0 = _mm512_add_epi32(v0, u1);
        highbd_clamp_epi32_avx512(&v0, bd);

        u0 = _mm512_cvtepi32_epi16(v0);

        _mm256_storeu_si256((__m256i *)output, u0);

        output += stride;
        i += 1;
    }
}

static INLINE void round_shift_16x16_avx512(__m512i *in, int32_t shift) {
    uint8_t ushift = (uint8_t)shift;
    __m512i rnding = _mm512_set1_epi32(1 << (ushift - 1));
    int32_t i = 0;

    while (i < 16) {
        in[i] = _mm512_add_epi32(in[i], rnding);
        in[i] = _mm512_srai_epi32(in[i], ushift);
        i++;
    }
}

static INLINE void iidentity16_and_round_shift_avx512(__m512i *input, int32_t shift)
{
    const __m512i scalar = _mm512_set1_epi32(NewSqrt2);
    const __m512i rnding = _mm512_set1_epi32((1 << (NewSqrt2Bits - 2)) +
        (!!(shift) << (shift + NewSqrt2Bits - 2)));

    for (int32_t i = 0; i < 16; i++) {
        input[i] = _mm512_mullo_epi32(input[i], scalar);
        input[i] = _mm512_add_epi32(input[i], rnding);
        input[i] = _mm512_srai_epi32(input[i], NewSqrt2Bits - 1 + shift);
    }
}

static INLINE void idct16_col_avx512(__m512i *in, __m512i *out, int32_t bit)
{

    const int32_t *cospi = cospi_arr(bit);
    const __m512i cospi60 = _mm512_set1_epi32(cospi[60]);
    const __m512i cospi28 = _mm512_set1_epi32(cospi[28]);
    const __m512i cospi44 = _mm512_set1_epi32(cospi[44]);
    const __m512i cospi12 = _mm512_set1_epi32(cospi[12]);
    const __m512i cospi52 = _mm512_set1_epi32(cospi[52]);
    const __m512i cospi20 = _mm512_set1_epi32(cospi[20]);
    const __m512i cospi36 = _mm512_set1_epi32(cospi[36]);
    const __m512i cospi4 = _mm512_set1_epi32(cospi[4]);
    const __m512i cospi56 = _mm512_set1_epi32(cospi[56]);
    const __m512i cospi24 = _mm512_set1_epi32(cospi[24]);
    const __m512i cospi40 = _mm512_set1_epi32(cospi[40]);
    const __m512i cospi8 = _mm512_set1_epi32(cospi[8]);
    const __m512i cospi32 = _mm512_set1_epi32(cospi[32]);
    const __m512i cospi48 = _mm512_set1_epi32(cospi[48]);
    const __m512i cospi16 = _mm512_set1_epi32(cospi[16]);
    const __m512i cospim4 = _mm512_set1_epi32(-cospi[4]);
    const __m512i cospim36 = _mm512_set1_epi32(-cospi[36]);
    const __m512i cospim20 = _mm512_set1_epi32(-cospi[20]);
    const __m512i cospim52 = _mm512_set1_epi32(-cospi[52]);
    const __m512i cospim8 = _mm512_set1_epi32(-cospi[8]);
    const __m512i cospim40 = _mm512_set1_epi32(-cospi[40]);
    const __m512i cospim32 = _mm512_set1_epi32(-cospi[32]);
    const __m512i cospim48 = _mm512_set1_epi32(-cospi[48]);
    const __m512i cospim16 = _mm512_set1_epi32(-cospi[16]);
    const __m512i rounding = _mm512_set1_epi32(1 << (bit - 1));
    __m512i tmp[16], tmp2[16];
    //stage 1

    //stage 2
    tmp[8] = half_btf_avx512(&cospi60, &in[1],
    &cospim4, &in[15], &rounding, bit);
    tmp[9] = half_btf_avx512(&cospi28, &in[9],
    &cospim36, &in[7], &rounding, bit);
    tmp[10] = half_btf_avx512(&cospi44, &in[5],
    &cospim20, &in[11], &rounding, bit);
    tmp[11] = half_btf_avx512(&cospi12, &in[13],
    &cospim52, &in[3], &rounding, bit);
    tmp[12] = half_btf_avx512(&cospi52, &in[13],
    &cospi12, &in[3], &rounding, bit);
    tmp[13] = half_btf_avx512(&cospi20, &in[5],
    &cospi44, &in[11], &rounding, bit);
    tmp[14] = half_btf_avx512(&cospi36, &in[9],
    &cospi28, &in[7], &rounding, bit);
    tmp[15] = half_btf_avx512(&cospi4, &in[1],
    &cospi60, &in[15], &rounding, bit);

    //stage 3
    tmp2[0] = half_btf_avx512(&cospi56, &in[2],
    &cospim8, &in[14], &rounding, bit);
    tmp2[1] = half_btf_avx512(&cospi24, &in[10],
    &cospim40, &in[6], &rounding, bit);
    tmp2[2] = half_btf_avx512(&cospi40, &in[10],
    &cospi24, &in[6], &rounding, bit);
    tmp2[3] = half_btf_avx512(&cospi8, &in[2],
    &cospi56, &in[14], &rounding, bit);
    tmp2[4] = _mm512_add_epi32(tmp[8], tmp[9]);
    tmp2[5] = _mm512_sub_epi32(tmp[8], tmp[9]);
    tmp2[6] = _mm512_sub_epi32(tmp[11], tmp[10]);
    tmp2[7] = _mm512_add_epi32(tmp[10], tmp[11]);
    tmp2[8] = _mm512_add_epi32(tmp[12], tmp[13]);
    tmp2[9] = _mm512_sub_epi32(tmp[12], tmp[13]);
    tmp2[10] = _mm512_sub_epi32(tmp[15], tmp[14]);
    tmp2[11] = _mm512_add_epi32(tmp[14], tmp[15]);

    //stage 4
    tmp[0] = half_btf_avx512(&cospi32, &in[0],
    &cospi32, &in[8], &rounding, bit);
    tmp[1] = half_btf_avx512(&cospi32, &in[0],
    &cospim32, &in[8], &rounding, bit);
    tmp[2] = half_btf_avx512(&cospi48, &in[4],
    &cospim16, &in[12], &rounding, bit);
    tmp[3] = half_btf_avx512(&cospi16, &in[4],
    &cospi48, &in[12], &rounding, bit);
    tmp[4] = _mm512_add_epi32(tmp2[0], tmp2[1]);
    tmp[5] = _mm512_sub_epi32(tmp2[0], tmp2[1]);
    tmp[6] = _mm512_sub_epi32(tmp2[3], tmp2[2]);
    tmp[7] = _mm512_add_epi32(tmp2[2], tmp2[3]);
    tmp[9] = half_btf_avx512(&cospim16, &tmp2[5],
    &cospi48, &tmp2[10], &rounding, bit);
    tmp[10] = half_btf_avx512(&cospim48, &tmp2[6],
    &cospim16, &tmp2[9], &rounding, bit);
    tmp[13] = half_btf_avx512(&cospim16, &tmp2[6],
    &cospi48, &tmp2[9], &rounding, bit);
    tmp[14] = half_btf_avx512(&cospi48, &tmp2[5],
    &cospi16, &tmp2[10], &rounding, bit);

    //stage 5
    tmp2[12] = _mm512_sub_epi32(tmp2[11], tmp2[8]);
    tmp2[15] = _mm512_add_epi32(tmp2[8], tmp2[11]);
    tmp2[8] = _mm512_add_epi32(tmp2[4], tmp2[7]);
    tmp2[11] = _mm512_sub_epi32(tmp2[4], tmp2[7]);
    tmp2[0] = _mm512_add_epi32(tmp[0], tmp[3]);
    tmp2[1] = _mm512_add_epi32(tmp[1], tmp[2]);
    tmp2[2] = _mm512_sub_epi32(tmp[1], tmp[2]);
    tmp2[3] = _mm512_sub_epi32(tmp[0], tmp[3]);
    tmp2[5] = half_btf_avx512(&cospim32, &tmp[5],
    &cospi32, &tmp[6], &rounding, bit);
    tmp2[6] = half_btf_avx512(&cospi32, &tmp[5],
    &cospi32, &tmp[6], &rounding, bit);
    tmp2[9] = _mm512_add_epi32(tmp[9], tmp[10]);
    tmp2[10] = _mm512_sub_epi32(tmp[9], tmp[10]);
    tmp2[13] = _mm512_sub_epi32(tmp[14], tmp[13]);
    tmp2[14] = _mm512_add_epi32(tmp[13], tmp[14]);

    //stage 6
    tmp[0] = _mm512_add_epi32(tmp2[0], tmp[7]);
    tmp[1] = _mm512_add_epi32(tmp2[1], tmp2[6]);
    tmp[2] = _mm512_add_epi32(tmp2[2], tmp2[5]);
    tmp[3] = _mm512_add_epi32(tmp2[3], tmp[4]);
    tmp[4] = _mm512_sub_epi32(tmp2[3], tmp[4]);
    tmp[5] = _mm512_sub_epi32(tmp2[2], tmp2[5]);
    tmp[6] = _mm512_sub_epi32(tmp2[1], tmp2[6]);
    tmp[7] = _mm512_sub_epi32(tmp2[0], tmp[7]);
    tmp[10] = half_btf_avx512(&cospim32, &tmp2[10],
    &cospi32, &tmp2[13], &rounding, bit);
    tmp[11] = half_btf_avx512(&cospim32, &tmp2[11],
    &cospi32, &tmp2[12], &rounding, bit);
    tmp[12] = half_btf_avx512(&cospi32, &tmp2[11],
    &cospi32, &tmp2[12], &rounding, bit);
    tmp[13] = half_btf_avx512(&cospi32, &tmp2[10],
    &cospi32, &tmp2[13], &rounding, bit);

    //stage 7
    out[0] = _mm512_add_epi32(tmp[0], tmp2[15]);
    out[1] = _mm512_add_epi32(tmp[1], tmp2[14]);
    out[2] = _mm512_add_epi32(tmp[2], tmp[13]);
    out[3] = _mm512_add_epi32(tmp[3], tmp[12]);
    out[4] = _mm512_add_epi32(tmp[4], tmp[11]);
    out[5] = _mm512_add_epi32(tmp[5], tmp[10]);
    out[6] = _mm512_add_epi32(tmp[6], tmp2[9]);
    out[7] = _mm512_add_epi32(tmp[7], tmp2[8]);
    out[8] = _mm512_sub_epi32(tmp[7], tmp2[8]);
    out[9] = _mm512_sub_epi32(tmp[6], tmp2[9]);
    out[10] = _mm512_sub_epi32(tmp[5], tmp[10]);
    out[11] = _mm512_sub_epi32(tmp[4], tmp[11]);
    out[12] = _mm512_sub_epi32(tmp[3], tmp[12]);
    out[13] = _mm512_sub_epi32(tmp[2], tmp[13]);
    out[14] = _mm512_sub_epi32(tmp[1], tmp2[14]);
    out[15] = _mm512_sub_epi32(tmp[0], tmp2[15]);
}


static INLINE void iadst16_col_avx512(__m512i *in, __m512i *out,
    int8_t cos_bit) {
    const int32_t *cospi = cospi_arr(cos_bit);
    const __m512i cospi2 = _mm512_set1_epi32(cospi[2]);
    const __m512i cospi62 = _mm512_set1_epi32(cospi[62]);
    const __m512i cospi10 = _mm512_set1_epi32(cospi[10]);
    const __m512i cospi54 = _mm512_set1_epi32(cospi[54]);
    const __m512i cospi18 = _mm512_set1_epi32(cospi[18]);
    const __m512i cospi46 = _mm512_set1_epi32(cospi[46]);
    const __m512i cospi26 = _mm512_set1_epi32(cospi[26]);
    const __m512i cospi38 = _mm512_set1_epi32(cospi[38]);
    const __m512i cospi34 = _mm512_set1_epi32(cospi[34]);
    const __m512i cospi30 = _mm512_set1_epi32(cospi[30]);
    const __m512i cospi42 = _mm512_set1_epi32(cospi[42]);
    const __m512i cospi22 = _mm512_set1_epi32(cospi[22]);
    const __m512i cospi50 = _mm512_set1_epi32(cospi[50]);
    const __m512i cospi14 = _mm512_set1_epi32(cospi[14]);
    const __m512i cospi58 = _mm512_set1_epi32(cospi[58]);
    const __m512i cospi6 = _mm512_set1_epi32(cospi[6]);

    const __m512i cospi8 = _mm512_set1_epi32(cospi[8]);
    const __m512i cospi56 = _mm512_set1_epi32(cospi[56]);
    const __m512i cospi40 = _mm512_set1_epi32(cospi[40]);
    const __m512i cospi24 = _mm512_set1_epi32(cospi[24]);

    const __m512i cospi16 = _mm512_set1_epi32(cospi[16]);
    const __m512i cospi48 = _mm512_set1_epi32(cospi[48]);
    const __m512i cospi32 = _mm512_set1_epi32(cospi[32]);

    const __m512i cospim2 = _mm512_set1_epi32(-cospi[2]);
    const __m512i cospim10 = _mm512_set1_epi32(-cospi[10]);
    const __m512i cospim18 = _mm512_set1_epi32(-cospi[18]);
    const __m512i cospim26 = _mm512_set1_epi32(-cospi[26]);
    const __m512i cospim34 = _mm512_set1_epi32(-cospi[34]);
    const __m512i cospim42 = _mm512_set1_epi32(-cospi[42]);
    const __m512i cospim50 = _mm512_set1_epi32(-cospi[50]);
    const __m512i cospim58 = _mm512_set1_epi32(-cospi[58]);
    const __m512i cospim56 = _mm512_set1_epi32(-cospi[56]);
    const __m512i cospim24 = _mm512_set1_epi32(-cospi[24]);
    const __m512i cospim8 = _mm512_set1_epi32(-cospi[8]);
    const __m512i cospim40 = _mm512_set1_epi32(-cospi[40]);
    const __m512i cospim48 = _mm512_set1_epi32(-cospi[48]);
    const __m512i cospim16 = _mm512_set1_epi32(-cospi[16]);
    const __m512i cospim32 = _mm512_set1_epi32(-cospi[32]);

    const __m256i negative = _mm256_set1_epi32(-1);
    const __m512i rounding = _mm512_set1_epi32(1 << (cos_bit - 1));

    __m512i tmp[16], tmp2[16], tmp3[16];
    __m256i temp1, temp2;
    //stage 1

    //stage 2
    tmp[0] = half_btf_avx512(&cospi2, &in[15],
    &cospi62, &in[0], &rounding, cos_bit);
    tmp[1] = half_btf_avx512(&cospi62, &in[15],
    &cospim2, &in[0], &rounding, cos_bit);
    tmp[2] = half_btf_avx512(&cospi10, &in[13],
    &cospi54, &in[2], &rounding, cos_bit);
    tmp[3] = half_btf_avx512(&cospi54, &in[13],
    &cospim10, &in[2], &rounding, cos_bit);
    tmp[4] = half_btf_avx512(&cospi18, &in[11],
    &cospi46, &in[4], &rounding, cos_bit);
    tmp[5] = half_btf_avx512(&cospi46, &in[11],
    &cospim18, &in[4], &rounding, cos_bit);
    tmp[6] = half_btf_avx512(&cospi26, &in[9],
    &cospi38, &in[6], &rounding, cos_bit);
    tmp[7] = half_btf_avx512(&cospi38, &in[9],
    &cospim26, &in[6], &rounding, cos_bit);
    tmp[8] = half_btf_avx512(&cospi34, &in[7],
    &cospi30, &in[8], &rounding, cos_bit);
    tmp[9] = half_btf_avx512(&cospi30, &in[7],
    &cospim34, &in[8], &rounding, cos_bit);
    tmp[10] = half_btf_avx512(&cospi42, &in[5],
    &cospi22, &in[10], &rounding, cos_bit);
    tmp[11] = half_btf_avx512(&cospi22, &in[5],
    &cospim42, &in[10], &rounding, cos_bit);
    tmp[12] = half_btf_avx512(&cospi50, &in[3],
    &cospi14, &in[12], &rounding, cos_bit);
    tmp[13] = half_btf_avx512(&cospi14, &in[3],
    &cospim50, &in[12], &rounding, cos_bit);
    tmp[14] = half_btf_avx512(&cospi58, &in[1],
    &cospi6, &in[14], &rounding, cos_bit);
    tmp[15] = half_btf_avx512(&cospi6, &in[1],
    &cospim58, &in[14], &rounding, cos_bit);

    //stage 3
    tmp3[0] = _mm512_add_epi32(tmp[0], tmp[8]);
    tmp3[1] = _mm512_add_epi32(tmp[1], tmp[9]);
    tmp3[2] = _mm512_add_epi32(tmp[2], tmp[10]);
    tmp3[3] = _mm512_add_epi32(tmp[3], tmp[11]);
    tmp3[4] = _mm512_add_epi32(tmp[4], tmp[12]);
    tmp3[5] = _mm512_add_epi32(tmp[5], tmp[13]);
    tmp3[6] = _mm512_add_epi32(tmp[6], tmp[14]);
    tmp3[7] = _mm512_add_epi32(tmp[7], tmp[15]);
    tmp2[8] = _mm512_sub_epi32(tmp[0], tmp[8]);
    tmp2[9] = _mm512_sub_epi32(tmp[1], tmp[9]);
    tmp2[10] = _mm512_sub_epi32(tmp[2], tmp[10]);
    tmp2[11] = _mm512_sub_epi32(tmp[3], tmp[11]);
    tmp2[12] = _mm512_sub_epi32(tmp[4], tmp[12]);
    tmp2[13] = _mm512_sub_epi32(tmp[5], tmp[13]);
    tmp2[14] = _mm512_sub_epi32(tmp[6], tmp[14]);
    tmp2[15] = _mm512_sub_epi32(tmp[7], tmp[15]);


    //stage 4
    tmp[8] = half_btf_avx512(
    &cospi8, &tmp2[8], &cospi56, &tmp2[9], &rounding, cos_bit);
    tmp[9] = half_btf_avx512(
    &cospi56, &tmp2[8], &cospim8, &tmp2[9], &rounding, cos_bit);
    tmp[10] = half_btf_avx512(
    &cospi40, &tmp2[10], &cospi24, &tmp2[11], &rounding, cos_bit);
    tmp[11] = half_btf_avx512(
    &cospi24, &tmp2[10], &cospim40, &tmp2[11], &rounding, cos_bit);
    tmp[12] = half_btf_avx512(
    &cospim56, &tmp2[12], &cospi8, &tmp2[13], &rounding, cos_bit);
    tmp[13] = half_btf_avx512(
    &cospi8, &tmp2[12], &cospi56, &tmp2[13], &rounding, cos_bit);
    tmp[14] = half_btf_avx512(
    &cospim24, &tmp2[14], &cospi40, &tmp2[15], &rounding, cos_bit);
    tmp[15] = half_btf_avx512(
    &cospi40, &tmp2[14], &cospi24, &tmp2[15], &rounding, cos_bit);

    //stage 5
    tmp3[8] = _mm512_add_epi32(tmp3[0], tmp3[4]);
    tmp3[9] = _mm512_add_epi32(tmp3[1], tmp3[5]);
    tmp3[10] = _mm512_add_epi32(tmp3[2], tmp3[6]);
    tmp3[11] = _mm512_add_epi32(tmp3[3], tmp3[7]);
    tmp2[4] = _mm512_sub_epi32(tmp3[0], tmp3[4]);
    tmp2[5] = _mm512_sub_epi32(tmp3[1], tmp3[5]);
    tmp2[6] = _mm512_sub_epi32(tmp3[2], tmp3[6]);
    tmp2[7] = _mm512_sub_epi32(tmp3[3], tmp3[7]);
    tmp3[12] = _mm512_add_epi32(tmp[8], tmp[12]);
    tmp3[13] = _mm512_add_epi32(tmp[9], tmp[13]);
    tmp3[14] = _mm512_add_epi32(tmp[10], tmp[14]);
    tmp3[15] = _mm512_add_epi32(tmp[11], tmp[15]);
    tmp2[12] = _mm512_sub_epi32(tmp[8], tmp[12]);
    tmp2[13] = _mm512_sub_epi32(tmp[9], tmp[13]);
    tmp2[14] = _mm512_sub_epi32(tmp[10], tmp[14]);
    tmp2[15] = _mm512_sub_epi32(tmp[11], tmp[15]);

    //stage 6
    tmp[4] = half_btf_avx512(
    &cospi16, &tmp2[4], &cospi48, &tmp2[5], &rounding, cos_bit);
    tmp[5] = half_btf_avx512(
    &cospi48, &tmp2[4], &cospim16, &tmp2[5], &rounding, cos_bit);
    tmp[6] = half_btf_avx512(
    &cospim48, &tmp2[6], &cospi16, &tmp2[7], &rounding, cos_bit);
    tmp[7] = half_btf_avx512(
    &cospi16, &tmp2[6], &cospi48, &tmp2[7], &rounding, cos_bit);
    tmp[12] = half_btf_avx512(
    &cospi16, &tmp2[12], &cospi48, &tmp2[13], &rounding, cos_bit);
    tmp[13] = half_btf_avx512(
    &cospi48, &tmp2[12], &cospim16, &tmp2[13], &rounding, cos_bit);
    tmp[14] = half_btf_avx512(
    &cospim48, &tmp2[14], &cospi16, &tmp2[15], &rounding, cos_bit);
    tmp[15] = half_btf_avx512(
    &cospi16, &tmp2[14], &cospi48, &tmp2[15], &rounding, cos_bit);

    //stage 7
    out[0] = _mm512_add_epi32(tmp3[8], tmp3[10]);
    out[2] = _mm512_add_epi32(tmp[12], tmp[14]);
    out[12] = _mm512_add_epi32(tmp[5], tmp[7]);
    out[14] = _mm512_add_epi32(tmp3[13], tmp3[15]);
    tmp2[1] = _mm512_add_epi32(tmp3[9], tmp3[11]);
    tmp2[2] = _mm512_sub_epi32(tmp3[8], tmp3[10]);
    tmp2[3] = _mm512_sub_epi32(tmp3[9], tmp3[11]);
    tmp2[4] = _mm512_add_epi32(tmp[4], tmp[6]);
    tmp2[6] = _mm512_sub_epi32(tmp[4], tmp[6]);
    tmp2[7] = _mm512_sub_epi32(tmp[5], tmp[7]);
    tmp2[8] = _mm512_add_epi32(tmp3[12], tmp3[14]);
    tmp2[10] = _mm512_sub_epi32(tmp3[12], tmp3[14]);
    tmp2[11] = _mm512_sub_epi32(tmp3[13], tmp3[15]);
    tmp2[13] = _mm512_add_epi32(tmp[13], tmp[15]);
    tmp2[14] = _mm512_sub_epi32(tmp[12], tmp[14]);
    tmp2[15] = _mm512_sub_epi32(tmp[13], tmp[15]);

    //stage 8
    out[4] = half_btf_avx512(
    &cospi32, &tmp2[6], &cospi32, &tmp2[7], &rounding, cos_bit);
    out[6] = half_btf_avx512(
    &cospi32, &tmp2[10], &cospi32, &tmp2[11], &rounding, cos_bit);
    out[8] = half_btf_avx512(
    &cospi32, &tmp2[2], &cospim32, &tmp2[3], &rounding, cos_bit);
    out[10] = half_btf_avx512(
    &cospi32, &tmp2[14], &cospim32, &tmp2[15], &rounding, cos_bit);
    tmp[2] = half_btf_avx512(
    &cospi32, &tmp2[2], &cospi32, &tmp2[3], &rounding, cos_bit);
    tmp[7] = half_btf_avx512(
    &cospi32, &tmp2[6], &cospim32, &tmp2[7], &rounding, cos_bit);
    tmp[11] = half_btf_avx512(
    &cospi32, &tmp2[10], &cospim32, &tmp2[11], &rounding, cos_bit);
    tmp[14] = half_btf_avx512(
    &cospi32, &tmp2[14], &cospi32, &tmp2[15], &rounding, cos_bit);

    //stage 9
    temp1 = _mm512_extracti64x4_epi64(tmp2[8], 0);
    temp2 = _mm512_extracti64x4_epi64(tmp2[8], 1);
    temp1 = _mm256_sign_epi32(temp1, negative);
    temp2 = _mm256_sign_epi32(temp2, negative);
    out[1] = _mm512_inserti64x4(out[1], temp1, 0);
    out[1] = _mm512_inserti64x4(out[1], temp2, 1);

    temp1 = _mm512_extracti64x4_epi64(tmp2[4], 0);
    temp2 = _mm512_extracti64x4_epi64(tmp2[4], 1);
    temp1 = _mm256_sign_epi32(temp1, negative);
    temp2 = _mm256_sign_epi32(temp2, negative);
    out[3] = _mm512_inserti64x4(out[3], temp1, 0);
    out[3] = _mm512_inserti64x4(out[3], temp2, 1);

    temp1 = _mm512_extracti64x4_epi64(tmp[14], 0);
    temp2 = _mm512_extracti64x4_epi64(tmp[14], 1);
    temp1 = _mm256_sign_epi32(temp1, negative);
    temp2 = _mm256_sign_epi32(temp2, negative);
    out[5] = _mm512_inserti64x4(out[5], temp1, 0);
    out[5] = _mm512_inserti64x4(out[5], temp2, 1);

    temp1 = _mm512_extracti64x4_epi64(tmp[2], 0);
    temp2 = _mm512_extracti64x4_epi64(tmp[2], 1);
    temp1 = _mm256_sign_epi32(temp1, negative);
    temp2 = _mm256_sign_epi32(temp2, negative);
    out[7] = _mm512_inserti64x4(out[7], temp1, 0);
    out[7] = _mm512_inserti64x4(out[7], temp2, 1);

    temp1 = _mm512_extracti64x4_epi64(tmp[11], 0);
    temp2 = _mm512_extracti64x4_epi64(tmp[11], 1);
    temp1 = _mm256_sign_epi32(temp1, negative);
    temp2 = _mm256_sign_epi32(temp2, negative);
    out[9] = _mm512_inserti64x4(out[9], temp1, 0);
    out[9] = _mm512_inserti64x4(out[9], temp2, 1);

    temp1 = _mm512_extracti64x4_epi64(tmp[7], 0);
    temp2 = _mm512_extracti64x4_epi64(tmp[7], 1);
    temp1 = _mm256_sign_epi32(temp1, negative);
    temp2 = _mm256_sign_epi32(temp2, negative);
    out[11] = _mm512_inserti64x4(out[11], temp1, 0);
    out[11] = _mm512_inserti64x4(out[11], temp2, 1);

    temp1 = _mm512_extracti64x4_epi64(tmp2[13], 0);
    temp2 = _mm512_extracti64x4_epi64(tmp2[13], 1);
    temp1 = _mm256_sign_epi32(temp1, negative);
    temp2 = _mm256_sign_epi32(temp2, negative);
    out[13] = _mm512_inserti64x4(out[13], temp1, 0);
    out[13] = _mm512_inserti64x4(out[13], temp2, 1);

    temp1 = _mm512_extracti64x4_epi64(tmp2[1], 0);
    temp2 = _mm512_extracti64x4_epi64(tmp2[1], 1);
    temp1 = _mm256_sign_epi32(temp1, negative);
    temp2 = _mm256_sign_epi32(temp2, negative);
    out[15] = _mm512_inserti64x4(out[15], temp1, 0);
    out[15] = _mm512_inserti64x4(out[15], temp2, 1);
}

void av1_inv_txfm2d_add_16x16_avx512(const int32_t *coeff, uint16_t *output,
    int32_t stride, TxType tx_type, int32_t bd) {
    __m512i in[16], out[16];
    const int8_t *shift = eb_inv_txfm_shift_ls[TX_16X16];
    const int32_t txw_idx = get_txw_idx(TX_16X16);
    const int32_t txh_idx = get_txh_idx(TX_16X16);

    switch (tx_type) {
    case IDTX:
        load_buffer_16x16_avx512(coeff, in);
        iidentity16_and_round_shift_avx512(in, -shift[0]);
        iidentity16_and_round_shift_avx512(in, -shift[1]);
        write_buffer_16x16_avx512(in, output, stride, 0, 0, bd);
        break;

    case V_DCT:
        load_buffer_16x16_avx512(coeff, in);
        iidentity16_and_round_shift_avx512(in, -shift[0]);
        idct16_col_avx512(in, out, inv_cos_bit_col[txw_idx][txh_idx]);
        round_shift_16x16_avx512(out, -shift[1]);
        write_buffer_16x16_avx512(out, output, stride, 0, 0, bd);
        break;

    case H_DCT:
        load_buffer_16x16_avx512(coeff, in);
        transpose_16x16_avx512(1 , in, out);
        idct16_col_avx512(out, in, inv_cos_bit_row[txw_idx][txh_idx]);
        transpose_16x16_avx512(1 , in, out);
        round_shift_16x16_avx512(out, -shift[0]);
        iidentity16_and_round_shift_avx512(out, -shift[1]);
        write_buffer_16x16_avx512(out, output, stride, 0, 0, bd);
        break;

    case V_ADST:
        load_buffer_16x16_avx512(coeff, in);
        iidentity16_and_round_shift_avx512(in, -shift[0]);
        iadst16_col_avx512(in, out, inv_cos_bit_col[txw_idx][txh_idx]);
        round_shift_16x16_avx512(out, -shift[1]);
        write_buffer_16x16_avx512(out, output, stride, 0, 0, bd);
        break;

    case H_ADST:
        load_buffer_16x16_avx512(coeff, in);
        transpose_16x16_avx512(1 , in, out);
        iadst16_col_avx512(out, in, inv_cos_bit_row[txw_idx][txh_idx]);
        transpose_16x16_avx512(1 , in, out);
        round_shift_16x16_avx512(out, -shift[0]);
        iidentity16_and_round_shift_avx512(out, -shift[1]);
        write_buffer_16x16_avx512(out, output, stride, 0, 0, bd);
        break;

    case V_FLIPADST:
        load_buffer_16x16_avx512(coeff, in);
        iidentity16_and_round_shift_avx512(in, -shift[0]);
        iadst16_col_avx512(in, out, inv_cos_bit_col[txw_idx][txh_idx]);
        round_shift_16x16_avx512(out, -shift[1]);
        write_buffer_16x16_avx512(out, output, stride, 0, 1, bd);
        break;

    case H_FLIPADST:
        load_buffer_16x16_avx512(coeff, in);
        transpose_16x16_avx512(1 , in, out);
        iadst16_col_avx512(out, in, inv_cos_bit_row[txw_idx][txh_idx]);
        transpose_16x16_avx512(1 , in, out);
        round_shift_16x16_avx512(out, -shift[0]);
        iidentity16_and_round_shift_avx512(out, -shift[1]);
        write_buffer_16x16_avx512(out, output, stride, 1, 0, bd);
        break;

    case DCT_DCT:
        load_buffer_16x16_avx512(coeff, in);
        transpose_16x16_avx512(1 , in, out);
        idct16_col_avx512(out, in, inv_cos_bit_row[txw_idx][txh_idx]);
        round_shift_16x16_avx512(in, -shift[0]);
        transpose_16x16_avx512(1 , in, out);
        idct16_col_avx512(out, in, inv_cos_bit_col[txw_idx][txh_idx]);
        round_shift_16x16_avx512(in, -shift[1]);
        write_buffer_16x16_avx512(in, output, stride, 0, 0, bd);
        break;

    case DCT_ADST:
        load_buffer_16x16_avx512(coeff, in);
        transpose_16x16_avx512(1 , in, out);
        iadst16_col_avx512(out, in, inv_cos_bit_row[txw_idx][txh_idx]);
        round_shift_16x16_avx512(in, -shift[0]);
        transpose_16x16_avx512(1 , in, out);
        idct16_col_avx512(out, in, inv_cos_bit_col[txw_idx][txh_idx]);
        round_shift_16x16_avx512(in, -shift[1]);
        write_buffer_16x16_avx512(in, output, stride, 0, 0, bd);
        break;

    case ADST_DCT:
        load_buffer_16x16_avx512(coeff, in);
        transpose_16x16_avx512(1 , in, out);
        idct16_col_avx512(out, in, inv_cos_bit_row[txw_idx][txh_idx]);
        round_shift_16x16_avx512(in, -shift[0]);
        transpose_16x16_avx512(1 , in, out);
        iadst16_col_avx512(out, in, inv_cos_bit_col[txw_idx][txh_idx]);
        round_shift_16x16_avx512(in, -shift[1]);
        write_buffer_16x16_avx512(in, output, stride, 0, 0, bd);
        break;

    case ADST_ADST:
        load_buffer_16x16_avx512(coeff, in);
        transpose_16x16_avx512(1 , in, out);
        iadst16_col_avx512(out, in, inv_cos_bit_row[txw_idx][txh_idx]);
        round_shift_16x16_avx512(in, -shift[0]);
        transpose_16x16_avx512(1 , in, out);
        iadst16_col_avx512(out, in, inv_cos_bit_col[txw_idx][txh_idx]);
        round_shift_16x16_avx512(in, -shift[1]);
        write_buffer_16x16_avx512(in, output, stride, 0, 0, bd);
        break;

    case FLIPADST_DCT:
        load_buffer_16x16_avx512(coeff, in);
        transpose_16x16_avx512(1 , in, out);
        idct16_col_avx512(out, in, inv_cos_bit_row[txw_idx][txh_idx]);
        round_shift_16x16_avx512(in, -shift[0]);
        transpose_16x16_avx512(1 , in, out);
        iadst16_col_avx512(out, in, inv_cos_bit_col[txw_idx][txh_idx]);
        round_shift_16x16_avx512(in, -shift[1]);
        write_buffer_16x16_avx512(in, output, stride, 0, 1, bd);
        break;

    case DCT_FLIPADST:
        load_buffer_16x16_avx512(coeff, in);
        transpose_16x16_avx512(1 , in, out);
        iadst16_col_avx512(out, in, inv_cos_bit_row[txw_idx][txh_idx]);
        round_shift_16x16_avx512(in, -shift[0]);
        transpose_16x16_avx512(1 , in, out);
        idct16_col_avx512(out, in, inv_cos_bit_col[txw_idx][txh_idx]);
        round_shift_16x16_avx512(in, -shift[1]);
        write_buffer_16x16_avx512(in, output, stride, 1, 0, bd);
        break;

    case ADST_FLIPADST:
        load_buffer_16x16_avx512(coeff, in);
        transpose_16x16_avx512(1 , in, out);
        iadst16_col_avx512(out, in, inv_cos_bit_row[txw_idx][txh_idx]);
        round_shift_16x16_avx512(in, -shift[0]);
        transpose_16x16_avx512(1 , in, out);
        iadst16_col_avx512(out, in, inv_cos_bit_col[txw_idx][txh_idx]);
        round_shift_16x16_avx512(in, -shift[1]);
        write_buffer_16x16_avx512(in, output, stride, 1, 0, bd);
        break;

    case FLIPADST_FLIPADST:
        load_buffer_16x16_avx512(coeff, in);
        transpose_16x16_avx512(1 , in, out);
        iadst16_col_avx512(out, in, inv_cos_bit_row[txw_idx][txh_idx]);
        round_shift_16x16_avx512(in, -shift[0]);
        transpose_16x16_avx512(1 , in, out);
        iadst16_col_avx512(out, in, inv_cos_bit_col[txw_idx][txh_idx]);
        round_shift_16x16_avx512(in, -shift[1]);
        write_buffer_16x16_avx512(in, output, stride, 1, 1, bd);
        break;

    case FLIPADST_ADST:
        load_buffer_16x16_avx512(coeff, in);
        transpose_16x16_avx512(1 , in, out);
        iadst16_col_avx512(out, in, inv_cos_bit_row[txw_idx][txh_idx]);
        round_shift_16x16_avx512(in, -shift[0]);
        transpose_16x16_avx512(1 , in, out);
        iadst16_col_avx512(out, in, inv_cos_bit_col[txw_idx][txh_idx]);
        round_shift_16x16_avx512(in, -shift[1]);
        write_buffer_16x16_avx512(in, output, stride, 0, 1, bd);
        break;

    default:assert(0);
    }
}
