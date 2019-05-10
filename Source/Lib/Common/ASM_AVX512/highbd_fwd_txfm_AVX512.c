#include <assert.h>
#include "EbDefinitions.h"
#include "aom_dsp_rtcd.h"
#include "EbTransforms.h"
#include <immintrin.h>

const int32_t *cospi_arr(int32_t n);
const int32_t *sinpi_arr(int32_t n);

#define TRANSPOSE_4X4_AVX512(x0, x1, x2, x3, y0, y1, y2, y3) \
  do {                                                \
    __m512i u0, u1, u2, u3;                           \
    u0 = _mm512_unpacklo_epi32(x0, x1);                  \
    u1 = _mm512_unpackhi_epi32(x0, x1);                  \
    u2 = _mm512_unpacklo_epi32(x2, x3);                  \
    u3 = _mm512_unpackhi_epi32(x2, x3);                  \
    y0 = _mm512_unpacklo_epi64(u0, u2);                  \
    y1 = _mm512_unpackhi_epi64(u0, u2);                  \
    y2 = _mm512_unpacklo_epi64(u1, u3);                  \
    y3 = _mm512_unpackhi_epi64(u1, u3);                  \
    } while (0)

static const int8_t *fwd_txfm_shift_ls[TX_SIZES_ALL] = {
    fwd_shift_4x4, fwd_shift_8x8, fwd_shift_16x16, fwd_shift_32x32,
    fwd_shift_64x64, fwd_shift_4x8, fwd_shift_8x4, fwd_shift_8x16,
    fwd_shift_16x8, fwd_shift_16x32, fwd_shift_32x16, fwd_shift_32x64,
    fwd_shift_64x32, fwd_shift_4x16, fwd_shift_16x4, fwd_shift_8x32,
    fwd_shift_32x8, fwd_shift_16x64, fwd_shift_64x16,
};

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

static INLINE void load_buffer_16x16_avx512(const int16_t *input, __m512i *out,
    int32_t stride, int32_t flipud, int32_t fliplr, int32_t shift) {
    __m256i temp[16];
    uint8_t ushift = (uint8_t)shift;
    if (flipud) {
        /* load rows upside down (bottom to top) */
        for (int32_t i = 0; i < 16; i++) {
            int idx = 15 - i;
            temp[idx] = _mm256_loadu_si256((const __m256i *)(input + i * stride));
            out[idx] = _mm512_cvtepi16_epi32(temp[idx]);
            out[idx] = _mm512_slli_epi32(out[idx], ushift);
        }
    }
    else {
        /* load rows normally */
        for (int32_t i = 0; i < 16; i++) {
            temp[i] = _mm256_loadu_si256((const __m256i *)(input + i * stride));
            out[i] = _mm512_cvtepi16_epi32(temp[i]);
            out[i] = _mm512_slli_epi32(out[i], ushift);
        }
    }

    if (fliplr) {
      /*flip columns left to right*/
        uint32_t idx[] = { 15 , 14 , 13 , 12 , 11 , 10 , 9 , 8 , 7 , 6 , 5 , 4 , 3 , 2 , 1 , 0 };
        __m512i index = _mm512_loadu_si512(idx);

        for (int32_t i = 0; i < 16; i++) {
            out[i] = _mm512_permutexvar_epi32(index, out[i]);
        }
    }

}

static void fidtx16x16_avx512(const __m512i *in, __m512i *out, int8_t bit, int32_t col_num) {
    (void)bit;
    const uint8_t bits = 12;       // NewSqrt2Bits = 12
    const int32_t sqrt = 2 * 5793; // 2 * NewSqrt2
    const __m512i newsqrt = _mm512_set1_epi32(sqrt);
    const __m512i rounding = _mm512_set1_epi32(1 << (bits - 1));
    __m512i temp;
    int32_t num_iters = 16 * col_num;
    for (int32_t i = 0; i < num_iters; i++) {
        temp = _mm512_mullo_epi32(in[i], newsqrt);
        temp = _mm512_add_epi32(temp, rounding);
        out[i] = _mm512_srai_epi32(temp, bits);
    }
}

static INLINE void col_txfm_16x16_rounding_avx512(__m512i *in, int32_t shift) {
    uint8_t ushift = (uint8_t)shift;
    const __m512i rounding = _mm512_set1_epi32(1 << (ushift - 1));
    for (int32_t i = 0; i < 16; i++) {
        in[i] = _mm512_add_epi32(in[i], rounding);
        in[i] = _mm512_srai_epi32(in[i], ushift);
    }
}

static INLINE void write_buffer_16x16(const __m512i *res, int32_t *output) {
    int32_t fact = -1, index = -1;
    for (int32_t i = 0; i < 8; i++)
    {
        _mm512_storeu_si512((__m512i *)(output + (++fact) * 32), res[++index]);
        _mm512_storeu_si512((__m512i *)(output + (fact) * 32 + 16), res[++index]);
      
    }
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

static void fadst16x16_avx512(const __m512i *in, __m512i *out, int8_t bit, const int32_t col_num) {
    const int32_t *cospi = cospi_arr(bit);
    const __m512i cospi32 = _mm512_set1_epi32(cospi[32]);
    const __m512i cospi48 = _mm512_set1_epi32(cospi[48]);
    const __m512i cospi16 = _mm512_set1_epi32(cospi[16]);
    const __m512i cospim16 = _mm512_set1_epi32(-cospi[16]);
    const __m512i cospim48 = _mm512_set1_epi32(-cospi[48]);
    const __m512i cospi8 = _mm512_set1_epi32(cospi[8]);
    const __m512i cospi56 = _mm512_set1_epi32(cospi[56]);
    const __m512i cospim56 = _mm512_set1_epi32(-cospi[56]);
    const __m512i cospim8 = _mm512_set1_epi32(-cospi[8]);
    const __m512i cospi24 = _mm512_set1_epi32(cospi[24]);
    const __m512i cospim24 = _mm512_set1_epi32(-cospi[24]);
    const __m512i cospim40 = _mm512_set1_epi32(-cospi[40]);
    const __m512i cospi40 = _mm512_set1_epi32(cospi[40]);
    const __m512i cospi2 = _mm512_set1_epi32(cospi[2]);
    const __m512i cospi62 = _mm512_set1_epi32(cospi[62]);
    const __m512i cospim2 = _mm512_set1_epi32(-cospi[2]);
    const __m512i cospi10 = _mm512_set1_epi32(cospi[10]);
    const __m512i cospi54 = _mm512_set1_epi32(cospi[54]);
    const __m512i cospim10 = _mm512_set1_epi32(-cospi[10]);
    const __m512i cospi18 = _mm512_set1_epi32(cospi[18]);
    const __m512i cospi46 = _mm512_set1_epi32(cospi[46]);
    const __m512i cospim18 = _mm512_set1_epi32(-cospi[18]);
    const __m512i cospi26 = _mm512_set1_epi32(cospi[26]);
    const __m512i cospi38 = _mm512_set1_epi32(cospi[38]);
    const __m512i cospim26 = _mm512_set1_epi32(-cospi[26]);
    const __m512i cospi34 = _mm512_set1_epi32(cospi[34]);
    const __m512i cospi30 = _mm512_set1_epi32(cospi[30]);
    const __m512i cospim34 = _mm512_set1_epi32(-cospi[34]);
    const __m512i cospi42 = _mm512_set1_epi32(cospi[42]);
    const __m512i cospi22 = _mm512_set1_epi32(cospi[22]);
    const __m512i cospim42 = _mm512_set1_epi32(-cospi[42]);
    const __m512i cospi50 = _mm512_set1_epi32(cospi[50]);
    const __m512i cospi14 = _mm512_set1_epi32(cospi[14]);
    const __m512i cospim50 = _mm512_set1_epi32(-cospi[50]);
    const __m512i cospi58 = _mm512_set1_epi32(cospi[58]);
    const __m512i cospi6 = _mm512_set1_epi32(cospi[6]);
    const __m512i cospim58 = _mm512_set1_epi32(-cospi[58]);
    const __m512i rnding = _mm512_set1_epi32(1 << (bit - 1));
    const __m512i zero = _mm512_setzero_si512();

    __m512i u[16], v[16], x, y;
    int32_t col;

    for (col = 0; col < col_num; ++col) {
        // stage 0
        // stage 1
        u[0] = in[0 * col_num + col];
        u[1] = _mm512_sub_epi32(zero, in[15 * col_num + col]);
        u[2] = _mm512_sub_epi32(zero, in[7 * col_num + col]);
        u[3] = in[8 * col_num + col];
        u[4] = _mm512_sub_epi32(zero, in[3 * col_num + col]);
        u[5] = in[12 * col_num + col];
        u[6] = in[4 * col_num + col];
        u[7] = _mm512_sub_epi32(zero, in[11 * col_num + col]);
        u[8] = _mm512_sub_epi32(zero, in[1 * col_num + col]);
        u[9] = in[14 * col_num + col];
        u[10] = in[6 * col_num + col];
        u[11] = _mm512_sub_epi32(zero, in[9 * col_num + col]);
        u[12] = in[2 * col_num + col];
        u[13] = _mm512_sub_epi32(zero, in[13 * col_num + col]);
        u[14] = _mm512_sub_epi32(zero, in[5 * col_num + col]);
        u[15] = in[10 * col_num + col];

        // stage 2
        v[0] = u[0];
        v[1] = u[1];

        x = _mm512_mullo_epi32(u[2], cospi32);
        y = _mm512_mullo_epi32(u[3], cospi32);
        v[2] = _mm512_add_epi32(x, y);
        v[2] = _mm512_add_epi32(v[2], rnding);
        v[2] = _mm512_srai_epi32(v[2], bit);

        v[3] = _mm512_sub_epi32(x, y);
        v[3] = _mm512_add_epi32(v[3], rnding);
        v[3] = _mm512_srai_epi32(v[3], bit);

        v[4] = u[4];
        v[5] = u[5];

        x = _mm512_mullo_epi32(u[6], cospi32);
        y = _mm512_mullo_epi32(u[7], cospi32);
        v[6] = _mm512_add_epi32(x, y);
        v[6] = _mm512_add_epi32(v[6], rnding);
        v[6] = _mm512_srai_epi32(v[6], bit);

        v[7] = _mm512_sub_epi32(x, y);
        v[7] = _mm512_add_epi32(v[7], rnding);
        v[7] = _mm512_srai_epi32(v[7], bit);

        v[8] = u[8];
        v[9] = u[9];

        x = _mm512_mullo_epi32(u[10], cospi32);
        y = _mm512_mullo_epi32(u[11], cospi32);
        v[10] = _mm512_add_epi32(x, y);
        v[10] = _mm512_add_epi32(v[10], rnding);
        v[10] = _mm512_srai_epi32(v[10], bit);

        v[11] = _mm512_sub_epi32(x, y);
        v[11] = _mm512_add_epi32(v[11], rnding);
        v[11] = _mm512_srai_epi32(v[11], bit);

        v[12] = u[12];
        v[13] = u[13];

        x = _mm512_mullo_epi32(u[14], cospi32);
        y = _mm512_mullo_epi32(u[15], cospi32);
        v[14] = _mm512_add_epi32(x, y);
        v[14] = _mm512_add_epi32(v[14], rnding);
        v[14] = _mm512_srai_epi32(v[14], bit);

        v[15] = _mm512_sub_epi32(x, y);
        v[15] = _mm512_add_epi32(v[15], rnding);
        v[15] = _mm512_srai_epi32(v[15], bit);

        // stage 3
        u[0] = _mm512_add_epi32(v[0], v[2]);
        u[1] = _mm512_add_epi32(v[1], v[3]);
        u[2] = _mm512_sub_epi32(v[0], v[2]);
        u[3] = _mm512_sub_epi32(v[1], v[3]);
        u[4] = _mm512_add_epi32(v[4], v[6]);
        u[5] = _mm512_add_epi32(v[5], v[7]);
        u[6] = _mm512_sub_epi32(v[4], v[6]);
        u[7] = _mm512_sub_epi32(v[5], v[7]);
        u[8] = _mm512_add_epi32(v[8], v[10]);
        u[9] = _mm512_add_epi32(v[9], v[11]);
        u[10] = _mm512_sub_epi32(v[8], v[10]);
        u[11] = _mm512_sub_epi32(v[9], v[11]);
        u[12] = _mm512_add_epi32(v[12], v[14]);
        u[13] = _mm512_add_epi32(v[13], v[15]);
        u[14] = _mm512_sub_epi32(v[12], v[14]);
        u[15] = _mm512_sub_epi32(v[13], v[15]);

        // stage 4
        v[0] = u[0];
        v[1] = u[1];
        v[2] = u[2];
        v[3] = u[3];
        v[4] = half_btf_avx512(&cospi16, &u[4], &cospi48, &u[5], &rnding, bit);
        v[5] = half_btf_avx512(&cospi48, &u[4], &cospim16, &u[5], &rnding, bit);
        v[6] = half_btf_avx512(&cospim48, &u[6], &cospi16, &u[7], &rnding, bit);
        v[7] = half_btf_avx512(&cospi16, &u[6], &cospi48, &u[7], &rnding, bit);
        v[8] = u[8];
        v[9] = u[9];
        v[10] = u[10];
        v[11] = u[11];
        v[12] = half_btf_avx512(&cospi16, &u[12], &cospi48, &u[13], &rnding, bit);
        v[13] = half_btf_avx512(&cospi48, &u[12], &cospim16, &u[13], &rnding, bit);
        v[14] = half_btf_avx512(&cospim48, &u[14], &cospi16, &u[15], &rnding, bit);
        v[15] = half_btf_avx512(&cospi16, &u[14], &cospi48, &u[15], &rnding, bit);

        // stage 5
        u[0] = _mm512_add_epi32(v[0], v[4]);
        u[1] = _mm512_add_epi32(v[1], v[5]);
        u[2] = _mm512_add_epi32(v[2], v[6]);
        u[3] = _mm512_add_epi32(v[3], v[7]);
        u[4] = _mm512_sub_epi32(v[0], v[4]);
        u[5] = _mm512_sub_epi32(v[1], v[5]);
        u[6] = _mm512_sub_epi32(v[2], v[6]);
        u[7] = _mm512_sub_epi32(v[3], v[7]);
        u[8] = _mm512_add_epi32(v[8], v[12]);
        u[9] = _mm512_add_epi32(v[9], v[13]);
        u[10] = _mm512_add_epi32(v[10], v[14]);
        u[11] = _mm512_add_epi32(v[11], v[15]);
        u[12] = _mm512_sub_epi32(v[8], v[12]);
        u[13] = _mm512_sub_epi32(v[9], v[13]);
        u[14] = _mm512_sub_epi32(v[10], v[14]);
        u[15] = _mm512_sub_epi32(v[11], v[15]);

        // stage 6
        v[0] = u[0];
        v[1] = u[1];
        v[2] = u[2];
        v[3] = u[3];
        v[4] = u[4];
        v[5] = u[5];
        v[6] = u[6];
        v[7] = u[7];
        v[8] = half_btf_avx512(&cospi8, &u[8], &cospi56, &u[9], &rnding, bit);
        v[9] = half_btf_avx512(&cospi56, &u[8], &cospim8, &u[9], &rnding, bit);
        v[10] = half_btf_avx512(&cospi40, &u[10], &cospi24, &u[11], &rnding, bit);
        v[11] = half_btf_avx512(&cospi24, &u[10], &cospim40, &u[11], &rnding, bit);
        v[12] = half_btf_avx512(&cospim56, &u[12], &cospi8, &u[13], &rnding, bit);
        v[13] = half_btf_avx512(&cospi8, &u[12], &cospi56, &u[13], &rnding, bit);
        v[14] = half_btf_avx512(&cospim24, &u[14], &cospi40, &u[15], &rnding, bit);
        v[15] = half_btf_avx512(&cospi40, &u[14], &cospi24, &u[15], &rnding, bit);

        // stage 7
        u[0] = _mm512_add_epi32(v[0], v[8]);
        u[1] = _mm512_add_epi32(v[1], v[9]);
        u[2] = _mm512_add_epi32(v[2], v[10]);
        u[3] = _mm512_add_epi32(v[3], v[11]);
        u[4] = _mm512_add_epi32(v[4], v[12]);
        u[5] = _mm512_add_epi32(v[5], v[13]);
        u[6] = _mm512_add_epi32(v[6], v[14]);
        u[7] = _mm512_add_epi32(v[7], v[15]);
        u[8] = _mm512_sub_epi32(v[0], v[8]);
        u[9] = _mm512_sub_epi32(v[1], v[9]);
        u[10] = _mm512_sub_epi32(v[2], v[10]);
        u[11] = _mm512_sub_epi32(v[3], v[11]);
        u[12] = _mm512_sub_epi32(v[4], v[12]);
        u[13] = _mm512_sub_epi32(v[5], v[13]);
        u[14] = _mm512_sub_epi32(v[6], v[14]);
        u[15] = _mm512_sub_epi32(v[7], v[15]);

        // stage 8
        v[0] = half_btf_avx512(&cospi2, &u[0], &cospi62, &u[1], &rnding, bit);
        v[1] = half_btf_avx512(&cospi62, &u[0], &cospim2, &u[1], &rnding, bit);
        v[2] = half_btf_avx512(&cospi10, &u[2], &cospi54, &u[3], &rnding, bit);
        v[3] = half_btf_avx512(&cospi54, &u[2], &cospim10, &u[3], &rnding, bit);
        v[4] = half_btf_avx512(&cospi18, &u[4], &cospi46, &u[5], &rnding, bit);
        v[5] = half_btf_avx512(&cospi46, &u[4], &cospim18, &u[5], &rnding, bit);
        v[6] = half_btf_avx512(&cospi26, &u[6], &cospi38, &u[7], &rnding, bit);
        v[7] = half_btf_avx512(&cospi38, &u[6], &cospim26, &u[7], &rnding, bit);
        v[8] = half_btf_avx512(&cospi34, &u[8], &cospi30, &u[9], &rnding, bit);
        v[9] = half_btf_avx512(&cospi30, &u[8], &cospim34, &u[9], &rnding, bit);
        v[10] = half_btf_avx512(&cospi42, &u[10], &cospi22, &u[11], &rnding, bit);
        v[11] = half_btf_avx512(&cospi22, &u[10], &cospim42, &u[11], &rnding, bit);
        v[12] = half_btf_avx512(&cospi50, &u[12], &cospi14, &u[13], &rnding, bit);
        v[13] = half_btf_avx512(&cospi14, &u[12], &cospim50, &u[13], &rnding, bit);
        v[14] = half_btf_avx512(&cospi58, &u[14], &cospi6, &u[15], &rnding, bit);
        v[15] = half_btf_avx512(&cospi6, &u[14], &cospim58, &u[15], &rnding, bit);

        // stage 9
        out[0 * col_num + col] = v[1];
        out[1 * col_num + col] = v[14];
        out[2 * col_num + col] = v[3];
        out[3 * col_num + col] = v[12];
        out[4 * col_num + col] = v[5];
        out[5 * col_num + col] = v[10];
        out[6 * col_num + col] = v[7];
        out[7 * col_num + col] = v[8];
        out[8 * col_num + col] = v[9];
        out[9 * col_num + col] = v[6];
        out[10 * col_num + col] = v[11];
        out[11 * col_num + col] = v[4];
        out[12 * col_num + col] = v[13];
        out[13 * col_num + col] = v[2];
        out[14 * col_num + col] = v[15];
        out[15 * col_num + col] = v[0];
    }
}
static void fdct16x16_avx512(const __m512i *in, __m512i *out, int8_t bit, const int32_t col_num) {
    const int32_t *cospi = cospi_arr(bit);
    const __m512i cospi32 = _mm512_set1_epi32(cospi[32]);
    const __m512i cospim32 = _mm512_set1_epi32(-cospi[32]);
    const __m512i cospi48 = _mm512_set1_epi32(cospi[48]);
    const __m512i cospi16 = _mm512_set1_epi32(cospi[16]);
    const __m512i cospim48 = _mm512_set1_epi32(-cospi[48]);
    const __m512i cospim16 = _mm512_set1_epi32(-cospi[16]);
    const __m512i cospi56 = _mm512_set1_epi32(cospi[56]);
    const __m512i cospi8 = _mm512_set1_epi32(cospi[8]);
    const __m512i cospi24 = _mm512_set1_epi32(cospi[24]);
    const __m512i cospi40 = _mm512_set1_epi32(cospi[40]);
    const __m512i cospi60 = _mm512_set1_epi32(cospi[60]);
    const __m512i cospi4 = _mm512_set1_epi32(cospi[4]);
    const __m512i cospi28 = _mm512_set1_epi32(cospi[28]);
    const __m512i cospi36 = _mm512_set1_epi32(cospi[36]);
    const __m512i cospi44 = _mm512_set1_epi32(cospi[44]);
    const __m512i cospi20 = _mm512_set1_epi32(cospi[20]);
    const __m512i cospi12 = _mm512_set1_epi32(cospi[12]);
    const __m512i cospi52 = _mm512_set1_epi32(cospi[52]);
    const __m512i rnding = _mm512_set1_epi32(1 << (bit - 1));
    __m512i u[16], v[16], x;
    int32_t col;

    for (col = 0; col < col_num; ++col) {
        // stage 0
        // stage 1
        u[0] = _mm512_add_epi32(in[0 * col_num + col], in[15 * col_num + col]);
        u[15] = _mm512_sub_epi32(in[0 * col_num + col], in[15 * col_num + col]);
        u[1] = _mm512_add_epi32(in[1 * col_num + col], in[14 * col_num + col]);
        u[14] = _mm512_sub_epi32(in[1 * col_num + col], in[14 * col_num + col]);
        u[2] = _mm512_add_epi32(in[2 * col_num + col], in[13 * col_num + col]);
        u[13] = _mm512_sub_epi32(in[2 * col_num + col], in[13 * col_num + col]);
        u[3] = _mm512_add_epi32(in[3 * col_num + col], in[12 * col_num + col]);
        u[12] = _mm512_sub_epi32(in[3 * col_num + col], in[12 * col_num + col]);
        u[4] = _mm512_add_epi32(in[4 * col_num + col], in[11 * col_num + col]);
        u[11] = _mm512_sub_epi32(in[4 * col_num + col], in[11 * col_num + col]);
        u[5] = _mm512_add_epi32(in[5 * col_num + col], in[10 * col_num + col]);
        u[10] = _mm512_sub_epi32(in[5 * col_num + col], in[10 * col_num + col]);
        u[6] = _mm512_add_epi32(in[6 * col_num + col], in[9 * col_num + col]);
        u[9] = _mm512_sub_epi32(in[6 * col_num + col], in[9 * col_num + col]);
        u[7] = _mm512_add_epi32(in[7 * col_num + col], in[8 * col_num + col]);
        u[8] = _mm512_sub_epi32(in[7 * col_num + col], in[8 * col_num + col]);

        // stage 2
        v[0] = _mm512_add_epi32(u[0], u[7]);
        v[7] = _mm512_sub_epi32(u[0], u[7]);
        v[1] = _mm512_add_epi32(u[1], u[6]);
        v[6] = _mm512_sub_epi32(u[1], u[6]);
        v[2] = _mm512_add_epi32(u[2], u[5]);
        v[5] = _mm512_sub_epi32(u[2], u[5]);
        v[3] = _mm512_add_epi32(u[3], u[4]);
        v[4] = _mm512_sub_epi32(u[3], u[4]);
        v[8] = u[8];
        v[9] = u[9];

        v[10] = _mm512_mullo_epi32(u[10], cospim32);
        x = _mm512_mullo_epi32(u[13], cospi32);
        v[10] = _mm512_add_epi32(v[10], x);
        v[10] = _mm512_add_epi32(v[10], rnding);
        v[10] = _mm512_srai_epi32(v[10], bit);

        v[13] = _mm512_mullo_epi32(u[10], cospi32);
        x = _mm512_mullo_epi32(u[13], cospim32);
        v[13] = _mm512_sub_epi32(v[13], x);
        v[13] = _mm512_add_epi32(v[13], rnding);
        v[13] = _mm512_srai_epi32(v[13], bit);

        v[11] = _mm512_mullo_epi32(u[11], cospim32);
        x = _mm512_mullo_epi32(u[12], cospi32);
        v[11] = _mm512_add_epi32(v[11], x);
        v[11] = _mm512_add_epi32(v[11], rnding);
        v[11] = _mm512_srai_epi32(v[11], bit);

        v[12] = _mm512_mullo_epi32(u[11], cospi32);
        x = _mm512_mullo_epi32(u[12], cospim32);
        v[12] = _mm512_sub_epi32(v[12], x);
        v[12] = _mm512_add_epi32(v[12], rnding);
        v[12] = _mm512_srai_epi32(v[12], bit);
        v[14] = u[14];
        v[15] = u[15];

        // stage 3
        u[0] = _mm512_add_epi32(v[0], v[3]);
        u[3] = _mm512_sub_epi32(v[0], v[3]);
        u[1] = _mm512_add_epi32(v[1], v[2]);
        u[2] = _mm512_sub_epi32(v[1], v[2]);
        u[4] = v[4];

        u[5] = _mm512_mullo_epi32(v[5], cospim32);
        x = _mm512_mullo_epi32(v[6], cospi32);
        u[5] = _mm512_add_epi32(u[5], x);
        u[5] = _mm512_add_epi32(u[5], rnding);
        u[5] = _mm512_srai_epi32(u[5], bit);

        u[6] = _mm512_mullo_epi32(v[5], cospi32);
        x = _mm512_mullo_epi32(v[6], cospim32);
        u[6] = _mm512_sub_epi32(u[6], x);
        u[6] = _mm512_add_epi32(u[6], rnding);
        u[6] = _mm512_srai_epi32(u[6], bit);

        u[7] = v[7];
        u[8] = _mm512_add_epi32(v[8], v[11]);
        u[11] = _mm512_sub_epi32(v[8], v[11]);
        u[9] = _mm512_add_epi32(v[9], v[10]);
        u[10] = _mm512_sub_epi32(v[9], v[10]);
        u[12] = _mm512_sub_epi32(v[15], v[12]);
        u[15] = _mm512_add_epi32(v[15], v[12]);
        u[13] = _mm512_sub_epi32(v[14], v[13]);
        u[14] = _mm512_add_epi32(v[14], v[13]);

        // stage 4
        u[0] = _mm512_mullo_epi32(u[0], cospi32);
        u[1] = _mm512_mullo_epi32(u[1], cospi32);
        v[0] = _mm512_add_epi32(u[0], u[1]);
        v[0] = _mm512_add_epi32(v[0], rnding);
        v[0] = _mm512_srai_epi32(v[0], bit);

        v[1] = _mm512_sub_epi32(u[0], u[1]);
        v[1] = _mm512_add_epi32(v[1], rnding);
        v[1] = _mm512_srai_epi32(v[1], bit);

        v[2] = _mm512_mullo_epi32(u[2], cospi48);
        x = _mm512_mullo_epi32(u[3], cospi16);
        v[2] = _mm512_add_epi32(v[2], x);
        v[2] = _mm512_add_epi32(v[2], rnding);
        v[2] = _mm512_srai_epi32(v[2], bit);

        v[3] = _mm512_mullo_epi32(u[2], cospi16);
        x = _mm512_mullo_epi32(u[3], cospi48);
        v[3] = _mm512_sub_epi32(x, v[3]);
        v[3] = _mm512_add_epi32(v[3], rnding);
        v[3] = _mm512_srai_epi32(v[3], bit);

        v[4] = _mm512_add_epi32(u[4], u[5]);
        v[5] = _mm512_sub_epi32(u[4], u[5]);
        v[6] = _mm512_sub_epi32(u[7], u[6]);
        v[7] = _mm512_add_epi32(u[7], u[6]);
        v[8] = u[8];

        v[9] = _mm512_mullo_epi32(u[9], cospim16);
        x = _mm512_mullo_epi32(u[14], cospi48);
        v[9] = _mm512_add_epi32(v[9], x);
        v[9] = _mm512_add_epi32(v[9], rnding);
        v[9] = _mm512_srai_epi32(v[9], bit);

        v[14] = _mm512_mullo_epi32(u[9], cospi48);
        x = _mm512_mullo_epi32(u[14], cospim16);
        v[14] = _mm512_sub_epi32(v[14], x);
        v[14] = _mm512_add_epi32(v[14], rnding);
        v[14] = _mm512_srai_epi32(v[14], bit);

        v[10] = _mm512_mullo_epi32(u[10], cospim48);
        x = _mm512_mullo_epi32(u[13], cospim16);
        v[10] = _mm512_add_epi32(v[10], x);
        v[10] = _mm512_add_epi32(v[10], rnding);
        v[10] = _mm512_srai_epi32(v[10], bit);

        v[13] = _mm512_mullo_epi32(u[10], cospim16);
        x = _mm512_mullo_epi32(u[13], cospim48);
        v[13] = _mm512_sub_epi32(v[13], x);
        v[13] = _mm512_add_epi32(v[13], rnding);
        v[13] = _mm512_srai_epi32(v[13], bit);

        v[11] = u[11];
        v[12] = u[12];
        v[15] = u[15];

        // stage 5
        u[0] = v[0];
        u[1] = v[1];
        u[2] = v[2];
        u[3] = v[3];

        u[4] = _mm512_mullo_epi32(v[4], cospi56);
        x = _mm512_mullo_epi32(v[7], cospi8);
        u[4] = _mm512_add_epi32(u[4], x);
        u[4] = _mm512_add_epi32(u[4], rnding);
        u[4] = _mm512_srai_epi32(u[4], bit);

        u[7] = _mm512_mullo_epi32(v[4], cospi8);
        x = _mm512_mullo_epi32(v[7], cospi56);
        u[7] = _mm512_sub_epi32(x, u[7]);
        u[7] = _mm512_add_epi32(u[7], rnding);
        u[7] = _mm512_srai_epi32(u[7], bit);

        u[5] = _mm512_mullo_epi32(v[5], cospi24);
        x = _mm512_mullo_epi32(v[6], cospi40);
        u[5] = _mm512_add_epi32(u[5], x);
        u[5] = _mm512_add_epi32(u[5], rnding);
        u[5] = _mm512_srai_epi32(u[5], bit);

        u[6] = _mm512_mullo_epi32(v[5], cospi40);
        x = _mm512_mullo_epi32(v[6], cospi24);
        u[6] = _mm512_sub_epi32(x, u[6]);
        u[6] = _mm512_add_epi32(u[6], rnding);
        u[6] = _mm512_srai_epi32(u[6], bit);

        u[8] = _mm512_add_epi32(v[8], v[9]);
        u[9] = _mm512_sub_epi32(v[8], v[9]);
        u[10] = _mm512_sub_epi32(v[11], v[10]);
        u[11] = _mm512_add_epi32(v[11], v[10]);
        u[12] = _mm512_add_epi32(v[12], v[13]);
        u[13] = _mm512_sub_epi32(v[12], v[13]);
        u[14] = _mm512_sub_epi32(v[15], v[14]);
        u[15] = _mm512_add_epi32(v[15], v[14]);

        // stage 6
        v[0] = u[0];
        v[1] = u[1];
        v[2] = u[2];
        v[3] = u[3];
        v[4] = u[4];
        v[5] = u[5];
        v[6] = u[6];
        v[7] = u[7];

        v[8] = _mm512_mullo_epi32(u[8], cospi60);
        x = _mm512_mullo_epi32(u[15], cospi4);
        v[8] = _mm512_add_epi32(v[8], x);
        v[8] = _mm512_add_epi32(v[8], rnding);
        v[8] = _mm512_srai_epi32(v[8], bit);

        v[15] = _mm512_mullo_epi32(u[8], cospi4);
        x = _mm512_mullo_epi32(u[15], cospi60);
        v[15] = _mm512_sub_epi32(x, v[15]);
        v[15] = _mm512_add_epi32(v[15], rnding);
        v[15] = _mm512_srai_epi32(v[15], bit);

        v[9] = _mm512_mullo_epi32(u[9], cospi28);
        x = _mm512_mullo_epi32(u[14], cospi36);
        v[9] = _mm512_add_epi32(v[9], x);
        v[9] = _mm512_add_epi32(v[9], rnding);
        v[9] = _mm512_srai_epi32(v[9], bit);

        v[14] = _mm512_mullo_epi32(u[9], cospi36);
        x = _mm512_mullo_epi32(u[14], cospi28);
        v[14] = _mm512_sub_epi32(x, v[14]);
        v[14] = _mm512_add_epi32(v[14], rnding);
        v[14] = _mm512_srai_epi32(v[14], bit);

        v[10] = _mm512_mullo_epi32(u[10], cospi44);
        x = _mm512_mullo_epi32(u[13], cospi20);
        v[10] = _mm512_add_epi32(v[10], x);
        v[10] = _mm512_add_epi32(v[10], rnding);
        v[10] = _mm512_srai_epi32(v[10], bit);

        v[13] = _mm512_mullo_epi32(u[10], cospi20);
        x = _mm512_mullo_epi32(u[13], cospi44);
        v[13] = _mm512_sub_epi32(x, v[13]);
        v[13] = _mm512_add_epi32(v[13], rnding);
        v[13] = _mm512_srai_epi32(v[13], bit);

        v[11] = _mm512_mullo_epi32(u[11], cospi12);
        x = _mm512_mullo_epi32(u[12], cospi52);
        v[11] = _mm512_add_epi32(v[11], x);
        v[11] = _mm512_add_epi32(v[11], rnding);
        v[11] = _mm512_srai_epi32(v[11], bit);

        v[12] = _mm512_mullo_epi32(u[11], cospi52);
        x = _mm512_mullo_epi32(u[12], cospi12);
        v[12] = _mm512_sub_epi32(x, v[12]);
        v[12] = _mm512_add_epi32(v[12], rnding);
        v[12] = _mm512_srai_epi32(v[12], bit);

        out[0 * col_num + col] = v[0];
        out[1 * col_num + col] = v[8];
        out[2 * col_num + col] = v[4];
        out[3 * col_num + col] = v[12];
        out[4 * col_num + col] = v[2];
        out[5 * col_num + col] = v[10];
        out[6 * col_num + col] = v[6];
        out[7 * col_num + col] = v[14];
        out[8 * col_num + col] = v[1];
        out[9 * col_num + col] = v[9];
        out[10 * col_num + col] = v[5];
        out[11 * col_num + col] = v[13];
        out[12 * col_num + col] = v[3];
        out[13 * col_num + col] = v[11];
        out[14 * col_num + col] = v[7];
        out[15 * col_num + col] = v[15];
    }
}

void av1_fwd_txfm2d_16x16_avx512(int16_t *input, int32_t *coeff, uint32_t stride, TxType tx_type, uint8_t  bd)
{
    __m512i in[16], out[16];
    const int8_t *shift = fwd_txfm_shift_ls[TX_16X16];
    const int32_t txw_idx = get_txw_idx(TX_16X16);
    const int32_t txh_idx = get_txh_idx(TX_16X16);
    const int32_t col_num = 1;
    switch (tx_type) {
    case IDTX:
        load_buffer_16x16_avx512(input, in, stride, 0, 0, shift[0]);
        fidtx16x16_avx512(in, out, fwd_cos_bit_col[txw_idx][txh_idx], col_num);
        col_txfm_16x16_rounding_avx512(out, -shift[1]);
        fidtx16x16_avx512(out, out, fwd_cos_bit_row[txw_idx][txh_idx], col_num);
        write_buffer_16x16(out, coeff);
        break;
    case DCT_DCT:
        load_buffer_16x16_avx512(input, in, stride, 0, 0, shift[0]);
        fdct16x16_avx512(in, out, fwd_cos_bit_col[txw_idx][txh_idx], col_num);
        col_txfm_16x16_rounding_avx512(out, -shift[1]);
        transpose_16x16_avx512(1, out, in);
        fdct16x16_avx512(in, out, fwd_cos_bit_row[txw_idx][txh_idx], col_num);
        transpose_16x16_avx512(1, out, in);
        write_buffer_16x16(in, coeff);
        break;
    case ADST_DCT:
        load_buffer_16x16_avx512(input, in, stride, 0, 0, shift[0]);
        fadst16x16_avx512(in, out, fwd_cos_bit_col[txw_idx][txh_idx], col_num);
        col_txfm_16x16_rounding_avx512(out, -shift[1]);
        transpose_16x16_avx512(1, out, in);
        fdct16x16_avx512(in, out, fwd_cos_bit_row[txw_idx][txh_idx], col_num);
        transpose_16x16_avx512(1, out, in);
        write_buffer_16x16(in, coeff);
        break;
    case DCT_ADST:
        load_buffer_16x16_avx512(input, in, stride, 0, 0, shift[0]);
        fdct16x16_avx512(in, out, fwd_cos_bit_col[txw_idx][txh_idx], col_num);
        col_txfm_16x16_rounding_avx512(out, -shift[1]);
        transpose_16x16_avx512(1, out, in);
        fadst16x16_avx512(in, out, fwd_cos_bit_row[txw_idx][txh_idx], col_num);
        transpose_16x16_avx512(1, out, in);
        write_buffer_16x16(in, coeff);
        break;
    case ADST_ADST:
        load_buffer_16x16_avx512(input, in, stride, 0, 0, shift[0]);
        fadst16x16_avx512(in, out, fwd_cos_bit_col[txw_idx][txh_idx], col_num);
        col_txfm_16x16_rounding_avx512(out, -shift[1]);
        transpose_16x16_avx512(1, out, in);
        fadst16x16_avx512(in, out, fwd_cos_bit_row[txw_idx][txh_idx], col_num);
        transpose_16x16_avx512(1, out, in);
        write_buffer_16x16(in, coeff);
        break;
    case DCT_FLIPADST:
        load_buffer_16x16_avx512(input, in, stride, 0, 1, shift[0]);
        fdct16x16_avx512(in, out, fwd_cos_bit_col[txw_idx][txh_idx], col_num);
        col_txfm_16x16_rounding_avx512(out, -shift[1]);
        transpose_16x16_avx512(1, out, in);
        fadst16x16_avx512(in, out, fwd_cos_bit_row[txw_idx][txh_idx], col_num);
        transpose_16x16_avx512(1, out, in);
        write_buffer_16x16(in, coeff);
        break;
    case FLIPADST_DCT:
        load_buffer_16x16_avx512(input, in, stride, 1, 0, shift[0]);
        fadst16x16_avx512(in, out, fwd_cos_bit_col[txw_idx][txh_idx], col_num);
        col_txfm_16x16_rounding_avx512(out, -shift[1]);
        transpose_16x16_avx512(1, out, in);
        fdct16x16_avx512(in, out, fwd_cos_bit_row[txw_idx][txh_idx], col_num);
        transpose_16x16_avx512(1, out, in);
        write_buffer_16x16(in, coeff);
        break;
    case FLIPADST_FLIPADST:
        load_buffer_16x16_avx512(input, in, stride, 1, 1, shift[0]);
        fadst16x16_avx512(in, out, fwd_cos_bit_col[txw_idx][txh_idx], col_num);
        col_txfm_16x16_rounding_avx512(out, -shift[1]);
        transpose_16x16_avx512(1, out, in);
        fadst16x16_avx512(in, out, fwd_cos_bit_row[txw_idx][txh_idx], col_num);
        transpose_16x16_avx512(1, out, in);
        write_buffer_16x16(in, coeff);
        break;
    case ADST_FLIPADST:
        load_buffer_16x16_avx512(input, in, stride, 0, 1, shift[0]);
        fadst16x16_avx512(in, out, fwd_cos_bit_col[txw_idx][txh_idx], col_num);
        col_txfm_16x16_rounding_avx512(out, -shift[1]);
        transpose_16x16_avx512(1, out, in);
        fadst16x16_avx512(in, out, fwd_cos_bit_row[txw_idx][txh_idx], col_num);
        transpose_16x16_avx512(1, out, in);
        write_buffer_16x16(in, coeff);
        break;
    case FLIPADST_ADST:
        load_buffer_16x16_avx512(input, in, stride, 1, 0, shift[0]);
        fadst16x16_avx512(in, out, fwd_cos_bit_col[txw_idx][txh_idx], col_num);
        col_txfm_16x16_rounding_avx512(out, -shift[1]);
        transpose_16x16_avx512(1, out, in);
        fadst16x16_avx512(in, out, fwd_cos_bit_row[txw_idx][txh_idx], col_num);
        transpose_16x16_avx512(1, out, in);
        write_buffer_16x16(in, coeff);
        break;
    case V_DCT:
        load_buffer_16x16_avx512(input, in, stride, 0, 0, shift[0]);
        fdct16x16_avx512(in, out, fwd_cos_bit_col[txw_idx][txh_idx], col_num);
        col_txfm_16x16_rounding_avx512(out, -shift[1]);
        fidtx16x16_avx512(out, out, fwd_cos_bit_row[txw_idx][txh_idx], col_num);
        write_buffer_16x16(out, coeff);
        break;
    case H_DCT:
        load_buffer_16x16_avx512(input, in, stride, 0, 0, shift[0]);
        fidtx16x16_avx512(in, in, fwd_cos_bit_col[txw_idx][txh_idx], col_num);
        col_txfm_16x16_rounding_avx512(in, -shift[1]);
        transpose_16x16_avx512(1, in, out);
        fdct16x16_avx512(out, in, fwd_cos_bit_row[txw_idx][txh_idx], col_num);
        transpose_16x16_avx512(1, in, out);
        write_buffer_16x16(out, coeff);
        break;
    case V_ADST:
        load_buffer_16x16_avx512(input, in, stride, 0, 0, shift[0]);
        fadst16x16_avx512(in, out, fwd_cos_bit_col[txw_idx][txh_idx], col_num);
        col_txfm_16x16_rounding_avx512(out, -shift[1]);
        fidtx16x16_avx512(out, out, fwd_cos_bit_row[txw_idx][txh_idx], col_num);
        write_buffer_16x16(out, coeff);
        break;
    case H_ADST:
        load_buffer_16x16_avx512(input, in, stride, 0, 0, shift[0]);
        fidtx16x16_avx512(in, in, fwd_cos_bit_col[txw_idx][txh_idx], col_num);
        col_txfm_16x16_rounding_avx512(in, -shift[1]);
        transpose_16x16_avx512(1, in, out);
        fadst16x16_avx512(out, in, fwd_cos_bit_row[txw_idx][txh_idx], col_num);
        transpose_16x16_avx512(1, in, out);
        write_buffer_16x16(out, coeff);
        break;
    case V_FLIPADST:
        load_buffer_16x16_avx512(input, in, stride, 1, 0, shift[0]);
        fadst16x16_avx512(in, out, fwd_cos_bit_col[txw_idx][txh_idx], col_num);
        col_txfm_16x16_rounding_avx512(out, -shift[1]);
        fidtx16x16_avx512(out, out, fwd_cos_bit_row[txw_idx][txh_idx], col_num);
        write_buffer_16x16(out, coeff);
        break;
    case H_FLIPADST:
        load_buffer_16x16_avx512(input, in, stride, 0, 1, shift[0]);
        fidtx16x16_avx512(in, in, fwd_cos_bit_col[txw_idx][txh_idx], col_num);
        col_txfm_16x16_rounding_avx512(in, -shift[1]);
        transpose_16x16_avx512(1, in, out);
        fadst16x16_avx512(out, in, fwd_cos_bit_row[txw_idx][txh_idx], col_num);
        transpose_16x16_avx512(1, in, out);
        write_buffer_16x16(out, coeff);
        break;
    default: assert(0);
    }
    (void)bd;

}

