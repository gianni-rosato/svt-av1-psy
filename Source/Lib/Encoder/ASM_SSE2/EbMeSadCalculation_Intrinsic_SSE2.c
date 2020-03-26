/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbDefinitions.h"
#include <emmintrin.h>
#include "stdint.h"

static INLINE void sad8x4x2_sse2_intrin(const uint8_t *src, const uint32_t src_stride,
                                        const uint8_t *ref, const uint32_t ref_stride,
                                        __m128i *sad8x8) {
    *sad8x8 = _mm_add_epi32(*sad8x8,
                            _mm_sad_epu8(_mm_loadu_si128((__m128i *)(src + 0 * src_stride)),
                                         _mm_loadu_si128((__m128i *)(ref + 0 * ref_stride))));
    *sad8x8 = _mm_add_epi32(*sad8x8,
                            _mm_sad_epu8(_mm_loadu_si128((__m128i *)(src + 2 * src_stride)),
                                         _mm_loadu_si128((__m128i *)(ref + 2 * ref_stride))));
    *sad8x8 = _mm_add_epi32(*sad8x8,
                            _mm_sad_epu8(_mm_loadu_si128((__m128i *)(src + 4 * src_stride)),
                                         _mm_loadu_si128((__m128i *)(ref + 4 * ref_stride))));
    *sad8x8 = _mm_add_epi32(*sad8x8,
                            _mm_sad_epu8(_mm_loadu_si128((__m128i *)(src + 6 * src_stride)),
                                         _mm_loadu_si128((__m128i *)(ref + 6 * ref_stride))));
}

void sad_calculation_8x8_16x16_sse2_intrin(uint8_t *src, uint32_t src_stride, uint8_t *ref,
                                           uint32_t ref_stride, uint32_t *p_best_sad_8x8,
                                           uint32_t *p_best_sad_16x16, uint32_t *p_best_mv8x8,
                                           uint32_t *p_best_mv16x16, uint32_t mv,
                                           uint32_t *p_sad16x16, EbBool sub_sad) {
    __m128i xmm_sad16x16, xmm_sad8x8[2], xmm_sad16x16_total, sad8x8_0_3, sad8x8_less_than_bitmask,
        xmm_n1;
    __m128i sad8x8_greater_or_eq_bitmask, best_mv_8x8, best_sad8x8, xmm_p_best_sad_8x8,
        xmm_p_best_mv_8x8, xmm_mv;

    xmm_sad8x8[0] = xmm_sad8x8[1] = _mm_setzero_si128();

    //sad8x8_0, sad8x8_1
    sad8x4x2_sse2_intrin(
        src + 0 * src_stride, src_stride, ref + 0 * ref_stride, ref_stride, &xmm_sad8x8[0]);

    //sad8x8_2, sad8x8_3
    sad8x4x2_sse2_intrin(
        src + 8 * src_stride, src_stride, ref + 8 * ref_stride, ref_stride, &xmm_sad8x8[1]);

    if (sub_sad) {
        xmm_sad8x8[0] = _mm_slli_epi32(xmm_sad8x8[0], 1);
        xmm_sad8x8[1] = _mm_slli_epi32(xmm_sad8x8[1], 1);
    } else {
        //sad8x8_0, sad8x8_1
        sad8x4x2_sse2_intrin(
            src + 1 * src_stride, src_stride, ref + 1 * ref_stride, ref_stride, &xmm_sad8x8[0]);

        //sad8x8_2, sad8x8_3
        sad8x4x2_sse2_intrin(
            src + 9 * src_stride, src_stride, ref + 9 * ref_stride, ref_stride, &xmm_sad8x8[1]);
    }

    xmm_sad16x16       = _mm_add_epi32(xmm_sad8x8[0], xmm_sad8x8[1]);
    xmm_sad16x16_total = _mm_add_epi32(_mm_srli_si128(xmm_sad16x16, 8), xmm_sad16x16);

    *p_sad16x16 = _mm_cvtsi128_si32(xmm_sad16x16_total);

    sad8x8_0_3 = _mm_packs_epi32(xmm_sad8x8[0], xmm_sad8x8[1]);

    xmm_mv = _mm_cvtsi64_si128(mv);
    xmm_mv = _mm_unpacklo_epi32(xmm_mv, xmm_mv);
    xmm_mv = _mm_unpacklo_epi64(xmm_mv, xmm_mv);

    xmm_p_best_sad_8x8 = _mm_loadu_si128((__m128i *)p_best_sad_8x8);
    xmm_p_best_mv_8x8  = _mm_loadu_si128((__m128i *)p_best_mv8x8);

    // sad8x8_0 < p_best_sad_8x8[0] for 0 to 3
    sad8x8_less_than_bitmask = _mm_cmplt_epi32(sad8x8_0_3, xmm_p_best_sad_8x8);

    xmm_n1 = _mm_cmpeq_epi8(xmm_sad8x8[0], xmm_sad8x8[0]);

    sad8x8_greater_or_eq_bitmask = _mm_sub_epi32(xmm_n1, sad8x8_less_than_bitmask);

    best_sad8x8 = _mm_or_si128(_mm_and_si128(xmm_p_best_sad_8x8, sad8x8_greater_or_eq_bitmask),
                               _mm_and_si128(sad8x8_less_than_bitmask, sad8x8_0_3));
    best_mv_8x8 = _mm_or_si128(_mm_and_si128(xmm_p_best_mv_8x8, sad8x8_greater_or_eq_bitmask),
                               _mm_and_si128(sad8x8_less_than_bitmask, xmm_mv));

    _mm_storeu_si128((__m128i *)p_best_sad_8x8, best_sad8x8);
    _mm_storeu_si128((__m128i *)p_best_mv8x8, best_mv_8x8);

    uint64_t sad16x16 = _mm_cvtsi128_si64(xmm_sad16x16_total);
    if (sad16x16 < p_best_sad_16x16[0]) {
        p_best_sad_16x16[0] = (uint32_t)sad16x16;
        p_best_mv16x16[0]   = _mm_cvtsi128_si32(xmm_mv);
    }
}

void sad_calculation_32x32_64x64_sse2_intrin(uint32_t *p_sad16x16, uint32_t *p_best_sad_32x32,
                                             uint32_t *p_best_sad_64x64, uint32_t *p_best_mv32x32,
                                             uint32_t *p_best_mv64x64, uint32_t mv) {
    __m128i xmm_n1, sad32x32_greater_than_bitmask, sad32x32_less_than_or_eq_bitmask, best_sad32x32,
        best_mv_32x32, xmm_mv;
    __m128i sad_16x16_0_7_lo, sad_16x16_0_7_hi, sad_16x16_8_15_lo, sad_16x16_8_15_hi, xmm_sad64x64,
        xmm_sad64x64_total, xmm_p_best_sad_32x32, xmm_p_best_mv_32x32;

    sad_16x16_0_7_lo  = _mm_unpacklo_epi32(_mm_loadu_si128((__m128i *)p_sad16x16),
                                          _mm_loadu_si128((__m128i *)(p_sad16x16 + 4)));
    sad_16x16_0_7_hi  = _mm_unpackhi_epi32(_mm_loadu_si128((__m128i *)p_sad16x16),
                                          _mm_loadu_si128((__m128i *)(p_sad16x16 + 4)));
    sad_16x16_8_15_lo = _mm_unpacklo_epi32(_mm_loadu_si128((__m128i *)(p_sad16x16 + 8)),
                                           _mm_loadu_si128((__m128i *)(p_sad16x16 + 12)));
    sad_16x16_8_15_hi = _mm_unpackhi_epi32(_mm_loadu_si128((__m128i *)(p_sad16x16 + 8)),
                                           _mm_loadu_si128((__m128i *)(p_sad16x16 + 12)));

    xmm_sad64x64 =
        _mm_add_epi32(_mm_add_epi32(_mm_unpacklo_epi64(sad_16x16_0_7_lo, sad_16x16_8_15_lo),
                                    _mm_unpackhi_epi64(sad_16x16_0_7_lo, sad_16x16_8_15_lo)),
                      _mm_add_epi32(_mm_unpacklo_epi64(sad_16x16_0_7_hi, sad_16x16_8_15_hi),
                                    _mm_unpackhi_epi64(sad_16x16_0_7_hi, sad_16x16_8_15_hi)));

    xmm_sad64x64_total = _mm_add_epi32(_mm_srli_si128(xmm_sad64x64, 8), xmm_sad64x64);

    xmm_sad64x64_total = _mm_add_epi32(_mm_srli_si128(xmm_sad64x64_total, 4), xmm_sad64x64_total);

    xmm_mv = _mm_cvtsi32_si128(mv);
    xmm_mv = _mm_unpacklo_epi32(xmm_mv, xmm_mv);
    xmm_mv = _mm_unpacklo_epi64(xmm_mv, xmm_mv);

    xmm_p_best_sad_32x32 = _mm_loadu_si128((__m128i *)p_best_sad_32x32);
    xmm_p_best_mv_32x32  = _mm_loadu_si128((__m128i *)p_best_mv32x32);

    sad32x32_greater_than_bitmask = _mm_cmpgt_epi32(
        xmm_p_best_sad_32x32, xmm_sad64x64); // _mm_cmplt_epi32(xmm_p_best_sad_32x32, xmm_sad64x64);

    xmm_n1 =
        _mm_cmpeq_epi8(xmm_mv, xmm_mv); // anything compared to itself is equal (get 0xFFFFFFFF)
    sad32x32_less_than_or_eq_bitmask = _mm_sub_epi32(xmm_n1, sad32x32_greater_than_bitmask);

    best_sad32x32 =
        _mm_or_si128(_mm_and_si128(xmm_p_best_sad_32x32, sad32x32_less_than_or_eq_bitmask),
                     _mm_and_si128(xmm_sad64x64, sad32x32_greater_than_bitmask));
    best_mv_32x32 =
        _mm_or_si128(_mm_and_si128(xmm_p_best_mv_32x32, sad32x32_less_than_or_eq_bitmask),
                     _mm_and_si128(xmm_mv, sad32x32_greater_than_bitmask));

    _mm_storeu_si128((__m128i *)p_best_sad_32x32, best_sad32x32);
    _mm_storeu_si128((__m128i *)p_best_mv32x32, best_mv_32x32);

    uint32_t sad64x64 = _mm_cvtsi128_si32(xmm_sad64x64_total);
    if (sad64x64 < p_best_sad_64x64[0]) {
        p_best_sad_64x64[0] = sad64x64;
        p_best_mv64x64[0]   = _mm_cvtsi128_si32(xmm_mv);
    }
}

void initialize_buffer_32bits_sse2_intrin(uint32_t *pointer, uint32_t count128, uint32_t count32,
                                          uint32_t value) {
    __m128i  xmm1, xmm2;
    uint32_t index128;
    xmm2 = _mm_cvtsi32_si128(value);
    xmm1 = _mm_or_si128(_mm_slli_si128(xmm2, 4), xmm2);
    xmm2 = _mm_or_si128(_mm_slli_si128(xmm1, 8), xmm1);

    for (index128 = 0; index128 < count128; ++index128) {
        _mm_storeu_si128((__m128i *)pointer, xmm2);
        pointer += 4;
    }
    if (count32 == 3) { //Initialize 96 bits
        _mm_storel_epi64((__m128i *)(pointer), xmm2);
        *(pointer + 2) = _mm_cvtsi128_si32(xmm2);
    } else if (count32 == 2) { // Initialize 64 bits
        _mm_storel_epi64((__m128i *)pointer, xmm2);
    } else if (count32 == 1) { // Initialize 32 bits
        *(pointer) = _mm_cvtsi128_si32(xmm2);
    }
}
