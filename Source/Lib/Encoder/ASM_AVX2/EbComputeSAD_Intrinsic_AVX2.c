/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <assert.h>

#include "EbComputeSAD_AVX2.h"
#include "EbDefinitions.h"
#include <immintrin.h>
#include "EbMemory_AVX2.h"
#include "EbComputeSAD.h"

#define UPDATE_BEST(s, k, offset)        \
    tem_sum_1 = _mm_extract_epi32(s, k); \
    if (tem_sum_1 < low_sum) {           \
        low_sum = tem_sum_1;             \
        x_best  = j + offset + k;        \
        y_best  = i;                     \
    }

void ext_sad_calculation_8x8_16x16_avx2_intrin(uint8_t *src, uint32_t src_stride, uint8_t *ref,
                                               uint32_t ref_stride, uint32_t *p_best_sad_8x8,
                                               uint32_t *p_best_sad_16x16, uint32_t *p_best_mv8x8,
                                               uint32_t *p_best_mv16x16, uint32_t mv,
                                               uint32_t *p_sad16x16, uint32_t *p_sad8x8,
                                               EbBool sub_sad) {
    __m128i xmm_sad16x16, xmm_sad16x16_total, sad8x8_0_3;
    __m128i sad8x8_less_than_bitmask, best_mv8x8;
    __m128i best_sad8x8, xmm_best_sad8x8, xmm_best_mv8x8;
    __m256i sad8x8_0_3_256, src_256, ref_256;
    __m128i xmm_mv = _mm_set1_epi32(mv);

    src_256        = _mm256_setr_m128i(_mm_loadu_si128((__m128i const *)(src + 0 * src_stride)),
                                _mm_loadu_si128((__m128i const *)(src + 8 * src_stride)));
    ref_256        = _mm256_setr_m128i(_mm_loadu_si128((__m128i const *)(ref + 0 * ref_stride)),
                                _mm_loadu_si128((__m128i const *)(ref + 8 * ref_stride)));
    sad8x8_0_3_256 = _mm256_sad_epu8(src_256, ref_256);

    src_256        = _mm256_setr_m128i(_mm_loadu_si128((__m128i const *)(src + 2 * src_stride)),
                                _mm_loadu_si128((__m128i const *)(src + 10 * src_stride)));
    ref_256        = _mm256_setr_m128i(_mm_loadu_si128((__m128i const *)(ref + 2 * ref_stride)),
                                _mm_loadu_si128((__m128i const *)(ref + 10 * ref_stride)));
    sad8x8_0_3_256 = _mm256_add_epi32(_mm256_sad_epu8(src_256, ref_256), sad8x8_0_3_256);

    src_256        = _mm256_setr_m128i(_mm_loadu_si128((__m128i const *)(src + 4 * src_stride)),
                                _mm_loadu_si128((__m128i const *)(src + 12 * src_stride)));
    ref_256        = _mm256_setr_m128i(_mm_loadu_si128((__m128i const *)(ref + 4 * ref_stride)),
                                _mm_loadu_si128((__m128i const *)(ref + 12 * ref_stride)));
    sad8x8_0_3_256 = _mm256_add_epi32(_mm256_sad_epu8(src_256, ref_256), sad8x8_0_3_256);

    src_256        = _mm256_setr_m128i(_mm_loadu_si128((__m128i const *)(src + 6 * src_stride)),
                                _mm_loadu_si128((__m128i const *)(src + 14 * src_stride)));
    ref_256        = _mm256_setr_m128i(_mm_loadu_si128((__m128i const *)(ref + 6 * ref_stride)),
                                _mm_loadu_si128((__m128i const *)(ref + 14 * ref_stride)));
    sad8x8_0_3_256 = _mm256_add_epi32(_mm256_sad_epu8(src_256, ref_256), sad8x8_0_3_256);

    if (sub_sad)
        sad8x8_0_3_256 = _mm256_slli_epi32(sad8x8_0_3_256, 1);
    else {
        src_256        = _mm256_setr_m128i(_mm_loadu_si128((__m128i const *)(src + 1 * src_stride)),
                                    _mm_loadu_si128((__m128i const *)(src + 9 * src_stride)));
        ref_256        = _mm256_setr_m128i(_mm_loadu_si128((__m128i const *)(ref + 1 * ref_stride)),
                                    _mm_loadu_si128((__m128i const *)(ref + 9 * ref_stride)));
        sad8x8_0_3_256 = _mm256_add_epi32(_mm256_sad_epu8(src_256, ref_256), sad8x8_0_3_256);

        src_256        = _mm256_setr_m128i(_mm_loadu_si128((__m128i const *)(src + 3 * src_stride)),
                                    _mm_loadu_si128((__m128i const *)(src + 11 * src_stride)));
        ref_256        = _mm256_setr_m128i(_mm_loadu_si128((__m128i const *)(ref + 3 * ref_stride)),
                                    _mm_loadu_si128((__m128i const *)(ref + 11 * ref_stride)));
        sad8x8_0_3_256 = _mm256_add_epi32(_mm256_sad_epu8(src_256, ref_256), sad8x8_0_3_256);

        src_256        = _mm256_setr_m128i(_mm_loadu_si128((__m128i const *)(src + 5 * src_stride)),
                                    _mm_loadu_si128((__m128i const *)(src + 13 * src_stride)));
        ref_256        = _mm256_setr_m128i(_mm_loadu_si128((__m128i const *)(ref + 5 * ref_stride)),
                                    _mm_loadu_si128((__m128i const *)(ref + 13 * ref_stride)));
        sad8x8_0_3_256 = _mm256_add_epi32(_mm256_sad_epu8(src_256, ref_256), sad8x8_0_3_256);

        src_256        = _mm256_setr_m128i(_mm_loadu_si128((__m128i const *)(src + 7 * src_stride)),
                                    _mm_loadu_si128((__m128i const *)(src + 15 * src_stride)));
        ref_256        = _mm256_setr_m128i(_mm_loadu_si128((__m128i const *)(ref + 7 * ref_stride)),
                                    _mm_loadu_si128((__m128i const *)(ref + 15 * ref_stride)));
        sad8x8_0_3_256 = _mm256_add_epi32(_mm256_sad_epu8(src_256, ref_256), sad8x8_0_3_256);
    }

    sad8x8_0_3 = _mm_packs_epi32(_mm256_castsi256_si128(sad8x8_0_3_256),
                                 _mm256_extracti128_si256(sad8x8_0_3_256, 1));
    _mm_storeu_si128((__m128i *)p_sad8x8, sad8x8_0_3);

    xmm_best_sad8x8 = _mm_loadu_si128((__m128i const *)p_best_sad_8x8);
    xmm_best_mv8x8  = _mm_loadu_si128((__m128i const *)p_best_mv8x8);

    // sad8x8_0 < p_best_sad_8x8[0] for 0 to 3
    sad8x8_less_than_bitmask = _mm_cmplt_epi32(sad8x8_0_3, xmm_best_sad8x8);
    best_sad8x8 = _mm_blendv_epi8(xmm_best_sad8x8, sad8x8_0_3, sad8x8_less_than_bitmask);
    best_mv8x8  = _mm_blendv_epi8(xmm_best_mv8x8, xmm_mv, sad8x8_less_than_bitmask);
    _mm_storeu_si128((__m128i *)p_best_sad_8x8, best_sad8x8);
    _mm_storeu_si128((__m128i *)p_best_mv8x8, best_mv8x8);

    xmm_sad16x16       = _mm_add_epi32(_mm256_castsi256_si128(sad8x8_0_3_256),
                                 _mm256_extracti128_si256(sad8x8_0_3_256, 1));
    xmm_sad16x16_total = _mm_add_epi32(_mm_srli_si128(xmm_sad16x16, 8), xmm_sad16x16);
    p_sad16x16[0]      = _mm_cvtsi128_si32(xmm_sad16x16_total);

    if (p_sad16x16[0] < p_best_sad_16x16[0]) {
        p_best_sad_16x16[0] = p_sad16x16[0];
        p_best_mv16x16[0]   = mv;
    }
}
#if !REMOVE_UNUSED_CODE
void sad_loop_kernel_sparse_avx2_intrin(
    uint8_t * src, // input parameter, source samples Ptr
    uint32_t  src_stride, // input parameter, source stride
    uint8_t * ref, // input parameter, reference samples Ptr
    uint32_t  ref_stride, // input parameter, reference stride
    uint32_t  block_height, // input parameter, block height (M)
    uint32_t  block_width, // input parameter, block width (N)
    uint64_t *best_sad, int16_t *x_search_center, int16_t *y_search_center,
    uint32_t src_stride_raw, // input parameter, source stride (no line skipping)
    int16_t search_area_width, int16_t search_area_height) {
    int16_t        x_best = *x_search_center, y_best = *y_search_center;
    uint32_t       low_sum   = 0xffffff;
    uint32_t       tem_sum_1 = 0;
    int16_t        i, j;
    uint32_t       k, l;
    uint32_t       leftover = search_area_width & 7;
    const uint8_t *p_ref, *p_src;
    __m128i        s0, s1, s2, s3, s4, s5, s6, s8 = _mm_set1_epi32(-1);
    __m256i        ss0, ss1, ss2, ss3, ss4, ss5, ss6, ss7, ss8;

    if (leftover) {
        for (k = 0; k < leftover; k++) s8 = _mm_slli_si128(s8, 2);
    }

    switch (block_width) {
    case 4:

        if (!(block_height % 4)) {
            uint32_t src_stride_t = 3 * src_stride;
            uint32_t ref_stride_t = 3 * ref_stride;
            for (i = 0; i < search_area_height; i++) {
                uint32_t start_w = (i & 1) << 3;
                for (j = start_w; j <= search_area_width - 8; j += 16) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss5 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k += 4) {
                        ss0 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_ref)),
                            _mm_loadu_si128((__m128i *)(p_ref + 2 * ref_stride)),
                            0x1);
                        ss1 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(
                                _mm_loadu_si128((__m128i *)(p_ref + ref_stride))),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride_t)),
                            0x1);
                        ss2 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_unpacklo_epi64(
                                _mm_cvtsi32_si128(*(uint32_t *)p_src),
                                _mm_cvtsi32_si128(*(uint32_t *)(p_src + src_stride)))),
                            _mm_unpacklo_epi64(
                                _mm_cvtsi32_si128(*(uint32_t *)(p_src + 2 * src_stride)),
                                _mm_cvtsi32_si128(*(uint32_t *)(p_src + src_stride_t))),
                            0x1);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        p_src += src_stride << 2;
                        p_ref += ref_stride << 2;
                    }
                    ss3       = _mm256_adds_epu16(ss3, ss5);
                    s3        = _mm_adds_epu16(_mm256_castsi256_si128(ss3),
                                        _mm256_extracti128_si256(ss3, 1));
                    s3        = _mm_minpos_epu16(s3);
                    tem_sum_1 = _mm_extract_epi16(s3, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s3, 1));
                        y_best  = i;
                    }
                }

                if (leftover && j < search_area_width) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss5 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k += 4) {
                        ss0 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_ref)),
                            _mm_loadu_si128((__m128i *)(p_ref + 2 * ref_stride)),
                            0x1);
                        ss1 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(
                                _mm_loadu_si128((__m128i *)(p_ref + ref_stride))),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride_t)),
                            0x1);
                        ss2 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_unpacklo_epi64(
                                _mm_cvtsi32_si128(*(uint32_t *)p_src),
                                _mm_cvtsi32_si128(*(uint32_t *)(p_src + src_stride)))),
                            _mm_unpacklo_epi64(
                                _mm_cvtsi32_si128(*(uint32_t *)(p_src + 2 * src_stride)),
                                _mm_cvtsi32_si128(*(uint32_t *)(p_src + src_stride_t))),
                            0x1);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        p_src += src_stride << 2;
                        p_ref += ref_stride << 2;
                    }
                    ss3       = _mm256_adds_epu16(ss3, ss5);
                    s3        = _mm_adds_epu16(_mm256_castsi256_si128(ss3),
                                        _mm256_extracti128_si256(ss3, 1));
                    s3        = _mm_or_si128(s3, s8);
                    s3        = _mm_minpos_epu16(s3);
                    tem_sum_1 = _mm_extract_epi16(s3, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s3, 1));
                        y_best  = i;
                    }
                }
                ref += src_stride_raw;
            }
        } else {
            for (i = 0; i < search_area_height; i++) {
                uint32_t start_w = (i & 1) << 3;
                for (j = start_w; j <= search_area_width - 8; j += 16) {
                    p_src = src;
                    p_ref = ref + j;
                    s3    = _mm_setzero_si128();
                    for (k = 0; k < block_height; k += 2) {
                        s0 = _mm_loadu_si128((__m128i *)p_ref);
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + ref_stride));
                        s2 = _mm_cvtsi32_si128(*(uint32_t *)p_src);
                        s5 = _mm_cvtsi32_si128(*(uint32_t *)(p_src + src_stride));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s1, s5, 0));
                        p_src += src_stride << 1;
                        p_ref += ref_stride << 1;
                    }
                    s3        = _mm_minpos_epu16(s3);
                    tem_sum_1 = _mm_extract_epi16(s3, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s3, 1));
                        y_best  = i;
                    }
                }

                if (leftover && j < search_area_width) {
                    p_src = src;
                    p_ref = ref + j;
                    s3    = _mm_setzero_si128();
                    for (k = 0; k < block_height; k += 2) {
                        s0 = _mm_loadu_si128((__m128i *)p_ref);
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + ref_stride));
                        s2 = _mm_cvtsi32_si128(*(uint32_t *)p_src);
                        s5 = _mm_cvtsi32_si128(*(uint32_t *)(p_src + src_stride));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s1, s5, 0));
                        p_src += src_stride << 1;
                        p_ref += ref_stride << 1;
                    }
                    s3        = _mm_or_si128(s3, s8);
                    s3        = _mm_minpos_epu16(s3);
                    tem_sum_1 = _mm_extract_epi16(s3, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s3, 1));
                        y_best  = i;
                    }
                }
                ref += src_stride_raw;
            }
        }

        break;

    case 8:
        if (!(block_height % 4)) {
            uint32_t src_stride_t = 3 * src_stride;
            uint32_t ref_stride_t = 3 * ref_stride;
            for (i = 0; i < search_area_height; i++) {
                uint32_t start_w = (i & 1) << 3;
                for (j = start_w; j <= search_area_width - 8; j += 16) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k += 4) {
                        ss0 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_ref)),
                            _mm_loadu_si128((__m128i *)(p_ref + 2 * ref_stride)),
                            0x1);
                        ss1 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(
                                _mm_loadu_si128((__m128i *)(p_ref + ref_stride))),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride_t)),
                            0x1);
                        ss2 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_unpacklo_epi64(
                                _mm_loadl_epi64((__m128i *)p_src),
                                _mm_loadl_epi64((__m128i *)(p_src + src_stride)))),
                            _mm_unpacklo_epi64(_mm_loadl_epi64((__m128i *)(p_src + 2 * src_stride)),
                                               _mm_loadl_epi64((__m128i *)(p_src + src_stride_t))),
                            0x1);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride << 2;
                        p_ref += ref_stride << 2;
                    }
                    ss3 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    s3        = _mm_adds_epu16(_mm256_castsi256_si128(ss3),
                                        _mm256_extracti128_si256(ss3, 1));
                    s3        = _mm_minpos_epu16(s3);
                    tem_sum_1 = _mm_extract_epi16(s3, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s3, 1));
                        y_best  = i;
                    }
                }

                if (leftover && j < search_area_width) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k += 4) {
                        ss0 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_ref)),
                            _mm_loadu_si128((__m128i *)(p_ref + 2 * ref_stride)),
                            0x1);
                        ss1 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(
                                _mm_loadu_si128((__m128i *)(p_ref + ref_stride))),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride_t)),
                            0x1);
                        ss2 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_unpacklo_epi64(
                                _mm_loadl_epi64((__m128i *)p_src),
                                _mm_loadl_epi64((__m128i *)(p_src + src_stride)))),
                            _mm_unpacklo_epi64(_mm_loadl_epi64((__m128i *)(p_src + 2 * src_stride)),
                                               _mm_loadl_epi64((__m128i *)(p_src + src_stride_t))),
                            0x1);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride << 2;
                        p_ref += ref_stride << 2;
                    }
                    ss3 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    s3        = _mm_adds_epu16(_mm256_castsi256_si128(ss3),
                                        _mm256_extracti128_si256(ss3, 1));
                    s3        = _mm_or_si128(s3, s8);
                    s3        = _mm_minpos_epu16(s3);
                    tem_sum_1 = _mm_extract_epi16(s3, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s3, 1));
                        y_best  = i;
                    }
                }
                ref += src_stride_raw;
            }
        } else {
            for (i = 0; i < search_area_height; i++) {
                uint32_t start_w = (i & 1) << 3;
                for (j = start_w; j <= search_area_width - 8; j += 16) {
                    p_src = src;
                    p_ref = ref + j;
                    s3 = s4 = _mm_setzero_si128();
                    for (k = 0; k < block_height; k += 2) {
                        s0 = _mm_loadu_si128((__m128i *)p_ref);
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + ref_stride));
                        s2 = _mm_loadl_epi64((__m128i *)p_src);
                        s5 = _mm_loadl_epi64((__m128i *)(p_src + src_stride));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s1, s5, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s1, s5, 5));
                        p_src += src_stride << 1;
                        p_ref += ref_stride << 1;
                    }
                    s3        = _mm_adds_epu16(s3, s4);
                    s3        = _mm_minpos_epu16(s3);
                    tem_sum_1 = _mm_extract_epi16(s3, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s3, 1));
                        y_best  = i;
                    }
                }

                if (leftover && j < search_area_width) {
                    p_src = src;
                    p_ref = ref + j;
                    s3 = s4 = _mm_setzero_si128();
                    for (k = 0; k < block_height; k += 2) {
                        s0 = _mm_loadu_si128((__m128i *)p_ref);
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + ref_stride));
                        s2 = _mm_loadl_epi64((__m128i *)p_src);
                        s5 = _mm_loadl_epi64((__m128i *)(p_src + src_stride));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s1, s5, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s1, s5, 5));
                        p_src += src_stride << 1;
                        p_ref += ref_stride << 1;
                    }
                    s3        = _mm_adds_epu16(s3, s4);
                    s3        = _mm_or_si128(s3, s8);
                    s3        = _mm_minpos_epu16(s3);
                    tem_sum_1 = _mm_extract_epi16(s3, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s3, 1));
                        y_best  = i;
                    }
                }
                ref += src_stride_raw;
            }
        }
        break;

    case 16:
        if (block_height <= 16) {
            for (i = 0; i < search_area_height; i++) {
                uint32_t start_w = (i & 1) << 3;
                for (j = start_w; j <= search_area_width - 8; j += 16) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k += 2) {
                        ss0 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_ref)),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride)),
                            0x1);
                        ss1 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)(p_ref + 8))),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride + 8)),
                            0x1);
                        ss2 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_src)),
                            _mm_loadu_si128((__m128i *)(p_src + src_stride)),
                            0x1);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += 2 * src_stride;
                        p_ref += 2 * ref_stride;
                    }
                    ss3 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    s3        = _mm_adds_epu16(_mm256_castsi256_si128(ss3),
                                        _mm256_extracti128_si256(ss3, 1));
                    s3        = _mm_minpos_epu16(s3);
                    tem_sum_1 = _mm_extract_epi16(s3, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s3, 1));
                        y_best  = i;
                    }
                }

                if (leftover && j < search_area_width) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k += 2) {
                        ss0 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_ref)),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride)),
                            0x1);
                        ss1 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)(p_ref + 8))),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride + 8)),
                            0x1);
                        ss2 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_src)),
                            _mm_loadu_si128((__m128i *)(p_src + src_stride)),
                            0x1);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += 2 * src_stride;
                        p_ref += 2 * ref_stride;
                    }
                    ss3 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    s3        = _mm_adds_epu16(_mm256_castsi256_si128(ss3),
                                        _mm256_extracti128_si256(ss3, 1));
                    s3        = _mm_or_si128(s3, s8);
                    s3        = _mm_minpos_epu16(s3);
                    tem_sum_1 = _mm_extract_epi16(s3, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s3, 1));
                        y_best  = i;
                    }
                }
                ref += src_stride_raw;
            }
        } else if (block_height <= 32) {
            for (i = 0; i < search_area_height; i++) {
                uint32_t start_w = (i & 1) << 3;
                for (j = start_w; j <= search_area_width - 8; j += 16) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k += 2) {
                        ss0 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_ref)),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride)),
                            0x1);
                        ss1 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)(p_ref + 8))),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride + 8)),
                            0x1);
                        ss2 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_src)),
                            _mm_loadu_si128((__m128i *)(p_src + src_stride)),
                            0x1);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += 2 * src_stride;
                        p_ref += 2 * ref_stride;
                    }
                    ss3 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    s3        = _mm256_castsi256_si128(ss3);
                    s5        = _mm256_extracti128_si256(ss3, 1);
                    s4        = _mm_minpos_epu16(s3);
                    s6        = _mm_minpos_epu16(s5);
                    s4        = _mm_unpacklo_epi16(s4, s4);
                    s4        = _mm_unpacklo_epi32(s4, s4);
                    s4        = _mm_unpacklo_epi64(s4, s4);
                    s6        = _mm_unpacklo_epi16(s6, s6);
                    s6        = _mm_unpacklo_epi32(s6, s6);
                    s6        = _mm_unpacklo_epi64(s6, s6);
                    s3        = _mm_sub_epi16(s3, s4);
                    s5        = _mm_adds_epu16(s5, s3);
                    s5        = _mm_sub_epi16(s5, s6);
                    s5        = _mm_minpos_epu16(s5);
                    tem_sum_1 = _mm_extract_epi16(s5, 0);
                    tem_sum_1 += _mm_extract_epi16(s4, 0);
                    tem_sum_1 += _mm_extract_epi16(s6, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s5, 1));
                        y_best  = i;
                    }
                }

                if (leftover && j < search_area_width) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k += 2) {
                        ss0 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_ref)),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride)),
                            0x1);
                        ss1 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)(p_ref + 8))),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride + 8)),
                            0x1);
                        ss2 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_src)),
                            _mm_loadu_si128((__m128i *)(p_src + src_stride)),
                            0x1);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += 2 * src_stride;
                        p_ref += 2 * ref_stride;
                    }
                    ss3 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    s3        = _mm256_castsi256_si128(ss3);
                    s5        = _mm256_extracti128_si256(ss3, 1);
                    s3        = _mm_or_si128(s3, s8);
                    s5        = _mm_or_si128(s5, s8);
                    s4        = _mm_minpos_epu16(s3);
                    s6        = _mm_minpos_epu16(s5);
                    s4        = _mm_unpacklo_epi16(s4, s4);
                    s4        = _mm_unpacklo_epi32(s4, s4);
                    s4        = _mm_unpacklo_epi64(s4, s4);
                    s6        = _mm_unpacklo_epi16(s6, s6);
                    s6        = _mm_unpacklo_epi32(s6, s6);
                    s6        = _mm_unpacklo_epi64(s6, s6);
                    s3        = _mm_sub_epi16(s3, s4);
                    s5        = _mm_adds_epu16(s5, s3);
                    s5        = _mm_sub_epi16(s5, s6);
                    s5        = _mm_minpos_epu16(s5);
                    tem_sum_1 = _mm_extract_epi16(s5, 0);
                    tem_sum_1 += _mm_extract_epi16(s4, 0);
                    tem_sum_1 += _mm_extract_epi16(s6, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s5, 1));
                        y_best  = i;
                    }
                }
                ref += src_stride_raw;
            }
        } else {
            for (i = 0; i < search_area_height; i++) {
                uint32_t start_w = (i & 1) << 3;
                for (j = start_w; j <= search_area_width - 8; j += 16) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k += 2) {
                        ss0 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_ref)),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride)),
                            0x1);
                        ss1 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)(p_ref + 8))),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride + 8)),
                            0x1);
                        ss2 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_src)),
                            _mm_loadu_si128((__m128i *)(p_src + src_stride)),
                            0x1);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += 2 * src_stride;
                        p_ref += 2 * ref_stride;
                    }
                    ss3       = _mm256_adds_epu16(ss3, ss4);
                    ss5       = _mm256_adds_epu16(ss5, ss6);
                    ss0       = _mm256_adds_epu16(ss3, ss5);
                    s0        = _mm_adds_epu16(_mm256_castsi256_si128(ss0),
                                        _mm256_extracti128_si256(ss0, 1));
                    s0        = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            ss4 = _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256());
                            ss3 = _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256());
                            ss6 = _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256());
                            ss5 = _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256());
                            ss4 = _mm256_add_epi32(ss4, ss6);
                            ss3 = _mm256_add_epi32(ss3, ss5);
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss4),
                                               _mm256_extracti128_si256(ss4, 1));
                            s3  = _mm_add_epi32(_mm256_castsi256_si128(ss3),
                                               _mm256_extracti128_si256(ss3, 1));
                            UPDATE_BEST(s0, 0, 0);
                            UPDATE_BEST(s0, 1, 0);
                            UPDATE_BEST(s0, 2, 0);
                            UPDATE_BEST(s0, 3, 0);
                            UPDATE_BEST(s3, 0, 4);
                            UPDATE_BEST(s3, 1, 4);
                            UPDATE_BEST(s3, 2, 4);
                            UPDATE_BEST(s3, 3, 4);
                        }
                    }
                }

                if (leftover && j < search_area_width) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k += 2) {
                        ss0 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_ref)),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride)),
                            0x1);
                        ss1 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)(p_ref + 8))),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride + 8)),
                            0x1);
                        ss2 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_src)),
                            _mm_loadu_si128((__m128i *)(p_src + src_stride)),
                            0x1);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += 2 * src_stride;
                        p_ref += 2 * ref_stride;
                    }
                    ss3       = _mm256_adds_epu16(ss3, ss4);
                    ss5       = _mm256_adds_epu16(ss5, ss6);
                    ss0       = _mm256_adds_epu16(ss3, ss5);
                    s0        = _mm_adds_epu16(_mm256_castsi256_si128(ss0),
                                        _mm256_extracti128_si256(ss0, 1));
                    s0        = _mm_or_si128(s0, s8);
                    s0        = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            ss4 = _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256());
                            ss3 = _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256());
                            ss6 = _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256());
                            ss5 = _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256());
                            ss4 = _mm256_add_epi32(ss4, ss6);
                            ss3 = _mm256_add_epi32(ss3, ss5);
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss4),
                                               _mm256_extracti128_si256(ss4, 1));
                            s3  = _mm_add_epi32(_mm256_castsi256_si128(ss3),
                                               _mm256_extracti128_si256(ss3, 1));
                            k   = leftover;
                            while (k > 0) {
                                for (l = 0; l < 4 && k; l++, k--) {
                                    tem_sum_1 = _mm_extract_epi32(s0, 0);
                                    s0        = _mm_srli_si128(s0, 4);
                                    if (tem_sum_1 < low_sum) {
                                        low_sum = tem_sum_1;
                                        x_best  = (int16_t)(j + leftover - k);
                                        y_best  = i;
                                    }
                                }
                                s0 = s3;
                            }
                        }
                    }
                }
                ref += src_stride_raw;
            }
        }
        break;

    case 24:
        if (block_height <= 16) {
            for (i = 0; i < search_area_height; i++) {
                uint32_t start_w = (i & 1) << 3;
                for (j = start_w; j <= search_area_width - 8; j += 16) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    ss3       = _mm256_adds_epu16(ss3, ss4);
                    ss5       = _mm256_adds_epu16(ss5, ss6);
                    s3        = _mm_adds_epu16(_mm256_castsi256_si128(ss3),
                                        _mm256_extracti128_si256(ss3, 1));
                    s5        = _mm256_castsi256_si128(ss5);
                    s4        = _mm_minpos_epu16(s3);
                    s6        = _mm_minpos_epu16(s5);
                    s4        = _mm_unpacklo_epi16(s4, s4);
                    s4        = _mm_unpacklo_epi32(s4, s4);
                    s4        = _mm_unpacklo_epi64(s4, s4);
                    s6        = _mm_unpacklo_epi16(s6, s6);
                    s6        = _mm_unpacklo_epi32(s6, s6);
                    s6        = _mm_unpacklo_epi64(s6, s6);
                    s3        = _mm_sub_epi16(s3, s4);
                    s5        = _mm_adds_epu16(s5, s3);
                    s5        = _mm_sub_epi16(s5, s6);
                    s5        = _mm_minpos_epu16(s5);
                    tem_sum_1 = _mm_extract_epi16(s5, 0);
                    tem_sum_1 += _mm_extract_epi16(s4, 0);
                    tem_sum_1 += _mm_extract_epi16(s6, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s5, 1));
                        y_best  = i;
                    }
                }

                if (leftover && j < search_area_width) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    ss3       = _mm256_adds_epu16(ss3, ss4);
                    ss5       = _mm256_adds_epu16(ss5, ss6);
                    s3        = _mm_adds_epu16(_mm256_castsi256_si128(ss3),
                                        _mm256_extracti128_si256(ss3, 1));
                    s5        = _mm256_castsi256_si128(ss5);
                    s3        = _mm_or_si128(s3, s8);
                    s5        = _mm_or_si128(s5, s8);
                    s4        = _mm_minpos_epu16(s3);
                    s6        = _mm_minpos_epu16(s5);
                    s4        = _mm_unpacklo_epi16(s4, s4);
                    s4        = _mm_unpacklo_epi32(s4, s4);
                    s4        = _mm_unpacklo_epi64(s4, s4);
                    s6        = _mm_unpacklo_epi16(s6, s6);
                    s6        = _mm_unpacklo_epi32(s6, s6);
                    s6        = _mm_unpacklo_epi64(s6, s6);
                    s3        = _mm_sub_epi16(s3, s4);
                    s5        = _mm_adds_epu16(s5, s3);
                    s5        = _mm_sub_epi16(s5, s6);
                    s5        = _mm_minpos_epu16(s5);
                    tem_sum_1 = _mm_extract_epi16(s5, 0);
                    tem_sum_1 += _mm_extract_epi16(s4, 0);
                    tem_sum_1 += _mm_extract_epi16(s6, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s5, 1));
                        y_best  = i;
                    }
                }
                ref += src_stride_raw;
            }
        } else {
            for (i = 0; i < search_area_height; i++) {
                uint32_t start_w = (i & 1) << 3;
                for (j = start_w; j <= search_area_width - 8; j += 16) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    ss3       = _mm256_adds_epu16(ss3, ss4);
                    ss5       = _mm256_adds_epu16(ss5, ss6);
                    s3        = _mm256_castsi256_si128(ss3);
                    s4        = _mm256_extracti128_si256(ss3, 1);
                    s5        = _mm256_castsi256_si128(ss5);
                    s0        = _mm_adds_epu16(_mm_adds_epu16(s3, s4), s5);
                    s0        = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            s0 = _mm_unpacklo_epi16(s3, _mm_setzero_si128());
                            s3 = _mm_unpackhi_epi16(s3, _mm_setzero_si128());
                            s1 = _mm_unpacklo_epi16(s4, _mm_setzero_si128());
                            s4 = _mm_unpackhi_epi16(s4, _mm_setzero_si128());
                            s2 = _mm_unpacklo_epi16(s5, _mm_setzero_si128());
                            s5 = _mm_unpackhi_epi16(s5, _mm_setzero_si128());
                            s0 = _mm_add_epi32(_mm_add_epi32(s0, s1), s2);
                            s3 = _mm_add_epi32(_mm_add_epi32(s3, s4), s5);
                            UPDATE_BEST(s0, 0, 0);
                            UPDATE_BEST(s0, 1, 0);
                            UPDATE_BEST(s0, 2, 0);
                            UPDATE_BEST(s0, 3, 0);
                            UPDATE_BEST(s3, 0, 4);
                            UPDATE_BEST(s3, 1, 4);
                            UPDATE_BEST(s3, 2, 4);
                            UPDATE_BEST(s3, 3, 4);
                        }
                    }
                }

                if (leftover && j < search_area_width) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    ss3       = _mm256_adds_epu16(ss3, ss4);
                    ss5       = _mm256_adds_epu16(ss5, ss6);
                    s3        = _mm256_castsi256_si128(ss3);
                    s4        = _mm256_extracti128_si256(ss3, 1);
                    s5        = _mm256_castsi256_si128(ss5);
                    s0        = _mm_adds_epu16(_mm_adds_epu16(s3, s4), s5);
                    s0        = _mm_or_si128(s0, s8);
                    s0        = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            s0 = _mm_unpacklo_epi16(s3, _mm_setzero_si128());
                            s3 = _mm_unpackhi_epi16(s3, _mm_setzero_si128());
                            s1 = _mm_unpacklo_epi16(s4, _mm_setzero_si128());
                            s4 = _mm_unpackhi_epi16(s4, _mm_setzero_si128());
                            s2 = _mm_unpacklo_epi16(s5, _mm_setzero_si128());
                            s5 = _mm_unpackhi_epi16(s5, _mm_setzero_si128());
                            s0 = _mm_add_epi32(_mm_add_epi32(s0, s1), s2);
                            s3 = _mm_add_epi32(_mm_add_epi32(s3, s4), s5);
                            k  = leftover;
                            while (k > 0) {
                                for (l = 0; l < 4 && k; l++, k--) {
                                    tem_sum_1 = _mm_extract_epi32(s0, 0);
                                    s0        = _mm_srli_si128(s0, 4);
                                    if (tem_sum_1 < low_sum) {
                                        low_sum = tem_sum_1;
                                        x_best  = (int16_t)(j + leftover - k);
                                        y_best  = i;
                                    }
                                }
                                s0 = s3;
                            }
                        }
                    }
                }
                ref += src_stride_raw;
            }
        }
        break;

    case 32:
        if (block_height <= 16) {
            for (i = 0; i < search_area_height; i++) {
                uint32_t start_w = (i & 1) << 3;
                for (j = start_w; j <= search_area_width - 8; j += 16) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    ss3 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    s3        = _mm256_castsi256_si128(ss3);
                    s5        = _mm256_extracti128_si256(ss3, 1);
                    s4        = _mm_minpos_epu16(s3);
                    s6        = _mm_minpos_epu16(s5);
                    s4        = _mm_unpacklo_epi16(s4, s4);
                    s4        = _mm_unpacklo_epi32(s4, s4);
                    s4        = _mm_unpacklo_epi64(s4, s4);
                    s6        = _mm_unpacklo_epi16(s6, s6);
                    s6        = _mm_unpacklo_epi32(s6, s6);
                    s6        = _mm_unpacklo_epi64(s6, s6);
                    s3        = _mm_sub_epi16(s3, s4);
                    s5        = _mm_adds_epu16(s5, s3);
                    s5        = _mm_sub_epi16(s5, s6);
                    s5        = _mm_minpos_epu16(s5);
                    tem_sum_1 = _mm_extract_epi16(s5, 0);
                    tem_sum_1 += _mm_extract_epi16(s4, 0);
                    tem_sum_1 += _mm_extract_epi16(s6, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s5, 1));
                        y_best  = i;
                    }
                }

                if (leftover && j < search_area_width) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    ss3 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    s3        = _mm256_castsi256_si128(ss3);
                    s5        = _mm256_extracti128_si256(ss3, 1);
                    s3        = _mm_or_si128(s3, s8);
                    s5        = _mm_or_si128(s5, s8);
                    s4        = _mm_minpos_epu16(s3);
                    s6        = _mm_minpos_epu16(s5);
                    s4        = _mm_unpacklo_epi16(s4, s4);
                    s4        = _mm_unpacklo_epi32(s4, s4);
                    s4        = _mm_unpacklo_epi64(s4, s4);
                    s6        = _mm_unpacklo_epi16(s6, s6);
                    s6        = _mm_unpacklo_epi32(s6, s6);
                    s6        = _mm_unpacklo_epi64(s6, s6);
                    s3        = _mm_sub_epi16(s3, s4);
                    s5        = _mm_adds_epu16(s5, s3);
                    s5        = _mm_sub_epi16(s5, s6);
                    s5        = _mm_minpos_epu16(s5);
                    tem_sum_1 = _mm_extract_epi16(s5, 0);
                    tem_sum_1 += _mm_extract_epi16(s4, 0);
                    tem_sum_1 += _mm_extract_epi16(s6, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s5, 1));
                        y_best  = i;
                    }
                }
                ref += src_stride_raw;
            }
        } else if (block_height <= 32) {
            for (i = 0; i < search_area_height; i++) {
                uint32_t start_w = (i & 1) << 3;
                for (j = start_w; j <= search_area_width - 8; j += 16) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    ss3       = _mm256_adds_epu16(ss3, ss4);
                    ss5       = _mm256_adds_epu16(ss5, ss6);
                    ss6       = _mm256_adds_epu16(ss3, ss5);
                    s3        = _mm256_castsi256_si128(ss6);
                    s4        = _mm256_extracti128_si256(ss6, 1);
                    s0        = _mm_adds_epu16(s3, s4);
                    s0        = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            ss4 = _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256());
                            ss3 = _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256());
                            ss6 = _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256());
                            ss5 = _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256());
                            ss4 = _mm256_add_epi32(ss4, ss6);
                            ss3 = _mm256_add_epi32(ss3, ss5);
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss4),
                                               _mm256_extracti128_si256(ss4, 1));
                            s3  = _mm_add_epi32(_mm256_castsi256_si128(ss3),
                                               _mm256_extracti128_si256(ss3, 1));
                            UPDATE_BEST(s0, 0, 0);
                            UPDATE_BEST(s0, 1, 0);
                            UPDATE_BEST(s0, 2, 0);
                            UPDATE_BEST(s0, 3, 0);
                            UPDATE_BEST(s3, 0, 4);
                            UPDATE_BEST(s3, 1, 4);
                            UPDATE_BEST(s3, 2, 4);
                            UPDATE_BEST(s3, 3, 4);
                        }
                    }
                }

                if (leftover && j < search_area_width) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    ss3 = _mm256_adds_epu16(ss3, ss4);
                    ss5 = _mm256_adds_epu16(ss5, ss6);
                    ss6 = _mm256_adds_epu16(ss3, ss5);
                    s3  = _mm256_castsi256_si128(ss6);
                    s4  = _mm256_extracti128_si256(ss6, 1);
                    s0  = _mm_adds_epu16(s3, s4);
                    //s0 = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0        = _mm_or_si128(s0, s8);
                    s0        = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            ss4 = _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256());
                            ss3 = _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256());
                            ss6 = _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256());
                            ss5 = _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256());
                            ss4 = _mm256_add_epi32(ss4, ss6);
                            ss3 = _mm256_add_epi32(ss3, ss5);
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss4),
                                               _mm256_extracti128_si256(ss4, 1));
                            s3  = _mm_add_epi32(_mm256_castsi256_si128(ss3),
                                               _mm256_extracti128_si256(ss3, 1));
                            k   = leftover;
                            while (k > 0) {
                                for (l = 0; l < 4 && k; l++, k--) {
                                    tem_sum_1 = _mm_extract_epi32(s0, 0);
                                    s0        = _mm_srli_si128(s0, 4);
                                    if (tem_sum_1 < low_sum) {
                                        low_sum = tem_sum_1;
                                        x_best  = (int16_t)(j + leftover - k);
                                        y_best  = i;
                                    }
                                }
                                s0 = s3;
                            }
                        }
                    }
                }
                ref += src_stride_raw;
            }
        } else {
            for (i = 0; i < search_area_height; i++) {
                uint32_t start_w = (i & 1) << 3;
                for (j = start_w; j <= search_area_width - 8; j += 16) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    ss7 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    s3        = _mm256_castsi256_si128(ss7);
                    s4        = _mm256_extracti128_si256(ss7, 1);
                    s0        = _mm_adds_epu16(s3, s4);
                    s0        = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            ss0 = _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256());
                            ss3 = _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256());
                            ss1 = _mm256_unpacklo_epi16(ss4, _mm256_setzero_si256());
                            ss4 = _mm256_unpackhi_epi16(ss4, _mm256_setzero_si256());
                            ss2 = _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256());
                            ss5 = _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256());
                            ss7 = _mm256_unpacklo_epi16(ss6, _mm256_setzero_si256());
                            ss6 = _mm256_unpackhi_epi16(ss6, _mm256_setzero_si256());
                            ss0 = _mm256_add_epi32(_mm256_add_epi32(ss0, ss1),
                                                   _mm256_add_epi32(ss2, ss7));
                            ss3 = _mm256_add_epi32(_mm256_add_epi32(ss3, ss4),
                                                   _mm256_add_epi32(ss5, ss6));
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss0),
                                               _mm256_extracti128_si256(ss0, 1));
                            s3  = _mm_add_epi32(_mm256_castsi256_si128(ss3),
                                               _mm256_extracti128_si256(ss3, 1));
                            UPDATE_BEST(s0, 0, 0);
                            UPDATE_BEST(s0, 1, 0);
                            UPDATE_BEST(s0, 2, 0);
                            UPDATE_BEST(s0, 3, 0);
                            UPDATE_BEST(s3, 0, 4);
                            UPDATE_BEST(s3, 1, 4);
                            UPDATE_BEST(s3, 2, 4);
                            UPDATE_BEST(s3, 3, 4);
                        }
                    }
                }

                if (leftover && j < search_area_width) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    ss7 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    s3        = _mm256_castsi256_si128(ss7);
                    s4        = _mm256_extracti128_si256(ss7, 1);
                    s0        = _mm_adds_epu16(s3, s4);
                    s0        = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            ss0 = _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256());
                            ss3 = _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256());
                            ss1 = _mm256_unpacklo_epi16(ss4, _mm256_setzero_si256());
                            ss4 = _mm256_unpackhi_epi16(ss4, _mm256_setzero_si256());
                            ss2 = _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256());
                            ss5 = _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256());
                            ss7 = _mm256_unpacklo_epi16(ss6, _mm256_setzero_si256());
                            ss6 = _mm256_unpackhi_epi16(ss6, _mm256_setzero_si256());
                            ss0 = _mm256_add_epi32(_mm256_add_epi32(ss0, ss1),
                                                   _mm256_add_epi32(ss2, ss7));
                            ss3 = _mm256_add_epi32(_mm256_add_epi32(ss3, ss4),
                                                   _mm256_add_epi32(ss5, ss6));
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss0),
                                               _mm256_extracti128_si256(ss0, 1));
                            s3  = _mm_add_epi32(_mm256_castsi256_si128(ss3),
                                               _mm256_extracti128_si256(ss3, 1));
                            k   = leftover;
                            while (k > 0) {
                                for (l = 0; l < 4 && k; l++, k--) {
                                    tem_sum_1 = _mm_extract_epi32(s0, 0);
                                    s0        = _mm_srli_si128(s0, 4);
                                    if (tem_sum_1 < low_sum) {
                                        low_sum = tem_sum_1;
                                        x_best  = (int16_t)(j + leftover - k);
                                        y_best  = i;
                                    }
                                }
                                s0 = s3;
                            }
                        }
                    }
                }
                ref += src_stride_raw;
            }
        }
        break;

    case 48:
        if (block_height <= 32) {
            for (i = 0; i < search_area_height; i++) {
                uint32_t start_w = (i & 1) << 3;
                for (j = start_w; j <= search_area_width - 8; j += 16) {
                    p_src = src;
                    p_ref = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        s0  = _mm_loadu_si128((__m128i *)(p_ref + 32));
                        s1  = _mm_loadu_si128((__m128i *)(p_ref + 40));
                        s2  = _mm_loadu_si128((__m128i *)(p_src + 32));
                        s3  = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4  = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5  = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6  = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    s3        = _mm_adds_epu16(s3, s4);
                    s5        = _mm_adds_epu16(s5, s6);
                    s0        = _mm_adds_epu16(s3, s5);
                    ss3       = _mm256_adds_epu16(ss3, ss4);
                    ss5       = _mm256_adds_epu16(ss5, ss6);
                    ss6       = _mm256_adds_epu16(ss3, ss5);
                    s0        = _mm_adds_epu16(s0,
                                        _mm_adds_epu16(_mm256_castsi256_si128(ss6),
                                                       _mm256_extracti128_si256(ss6, 1)));
                    s0        = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            ss4 = _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256());
                            ss3 = _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256());
                            ss6 = _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256());
                            ss5 = _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256());
                            ss4 = _mm256_add_epi32(ss4, ss6);
                            ss3 = _mm256_add_epi32(ss3, ss5);
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss4),
                                               _mm256_extracti128_si256(ss4, 1));
                            s1  = _mm_add_epi32(_mm256_castsi256_si128(ss3),
                                               _mm256_extracti128_si256(ss3, 1));
                            s0  = _mm_add_epi32(s0, _mm_unpacklo_epi16(s3, _mm_setzero_si128()));
                            s0  = _mm_add_epi32(s0, _mm_unpacklo_epi16(s5, _mm_setzero_si128()));
                            s1  = _mm_add_epi32(s1, _mm_unpackhi_epi16(s3, _mm_setzero_si128()));
                            s1  = _mm_add_epi32(s1, _mm_unpackhi_epi16(s5, _mm_setzero_si128()));
                            UPDATE_BEST(s0, 0, 0);
                            UPDATE_BEST(s0, 1, 0);
                            UPDATE_BEST(s0, 2, 0);
                            UPDATE_BEST(s0, 3, 0);
                            UPDATE_BEST(s1, 0, 4);
                            UPDATE_BEST(s1, 1, 4);
                            UPDATE_BEST(s1, 2, 4);
                            UPDATE_BEST(s1, 3, 4);
                        }
                    }
                }

                if (leftover && j < search_area_width) {
                    p_src = src;
                    p_ref = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        s0  = _mm_loadu_si128((__m128i *)(p_ref + 32));
                        s1  = _mm_loadu_si128((__m128i *)(p_ref + 40));
                        s2  = _mm_loadu_si128((__m128i *)(p_src + 32));
                        s3  = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4  = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5  = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6  = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    s3        = _mm_adds_epu16(s3, s4);
                    s5        = _mm_adds_epu16(s5, s6);
                    s0        = _mm_adds_epu16(s3, s5);
                    ss3       = _mm256_adds_epu16(ss3, ss4);
                    ss5       = _mm256_adds_epu16(ss5, ss6);
                    ss6       = _mm256_adds_epu16(ss3, ss5);
                    s0        = _mm_adds_epu16(s0,
                                        _mm_adds_epu16(_mm256_castsi256_si128(ss6),
                                                       _mm256_extracti128_si256(ss6, 1)));
                    s0        = _mm_or_si128(s0, s8);
                    s0        = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            ss4 = _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256());
                            ss3 = _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256());
                            ss6 = _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256());
                            ss5 = _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256());
                            ss4 = _mm256_add_epi32(ss4, ss6);
                            ss3 = _mm256_add_epi32(ss3, ss5);
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss4),
                                               _mm256_extracti128_si256(ss4, 1));
                            s1  = _mm_add_epi32(_mm256_castsi256_si128(ss3),
                                               _mm256_extracti128_si256(ss3, 1));
                            s0  = _mm_add_epi32(s0, _mm_unpacklo_epi16(s3, _mm_setzero_si128()));
                            s0  = _mm_add_epi32(s0, _mm_unpacklo_epi16(s5, _mm_setzero_si128()));
                            s1  = _mm_add_epi32(s1, _mm_unpackhi_epi16(s3, _mm_setzero_si128()));
                            s1  = _mm_add_epi32(s1, _mm_unpackhi_epi16(s5, _mm_setzero_si128()));
                            k   = leftover;
                            while (k > 0) {
                                for (l = 0; l < 4 && k; l++, k--) {
                                    tem_sum_1 = _mm_extract_epi32(s0, 0);
                                    s0        = _mm_srli_si128(s0, 4);
                                    if (tem_sum_1 < low_sum) {
                                        low_sum = tem_sum_1;
                                        x_best  = (int16_t)(j + leftover - k);
                                        y_best  = i;
                                    }
                                }
                                s0 = s1;
                            }
                        }
                    }
                }
                ref += src_stride_raw;
            }
        } else {
            for (i = 0; i < search_area_height; i++) {
                uint32_t start_w = (i & 1) << 3;
                for (j = start_w; j <= search_area_width - 8; j += 16) {
                    p_src = src;
                    p_ref = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        s0  = _mm_loadu_si128((__m128i *)(p_ref + 32));
                        s1  = _mm_loadu_si128((__m128i *)(p_ref + 40));
                        s2  = _mm_loadu_si128((__m128i *)(p_src + 32));
                        s3  = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4  = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5  = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6  = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    s0        = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    ss7       = _mm256_adds_epu16(ss3, ss4);
                    ss8       = _mm256_adds_epu16(ss5, ss6);
                    ss7       = _mm256_adds_epu16(ss7, ss8);
                    s0        = _mm_adds_epu16(s0,
                                        _mm_adds_epu16(_mm256_castsi256_si128(ss7),
                                                       _mm256_extracti128_si256(ss7, 1)));
                    s0        = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            ss0 = _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256());
                            ss3 = _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256());
                            ss1 = _mm256_unpacklo_epi16(ss4, _mm256_setzero_si256());
                            ss4 = _mm256_unpackhi_epi16(ss4, _mm256_setzero_si256());
                            ss2 = _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256());
                            ss5 = _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256());
                            ss7 = _mm256_unpacklo_epi16(ss6, _mm256_setzero_si256());
                            ss6 = _mm256_unpackhi_epi16(ss6, _mm256_setzero_si256());
                            ss0 = _mm256_add_epi32(_mm256_add_epi32(ss0, ss1),
                                                   _mm256_add_epi32(ss2, ss7));
                            ss3 = _mm256_add_epi32(_mm256_add_epi32(ss3, ss4),
                                                   _mm256_add_epi32(ss5, ss6));
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss0),
                                               _mm256_extracti128_si256(ss0, 1));
                            s1  = _mm_add_epi32(_mm256_castsi256_si128(ss3),
                                               _mm256_extracti128_si256(ss3, 1));
                            s0  = _mm_add_epi32(s0, _mm_unpacklo_epi16(s3, _mm_setzero_si128()));
                            s0  = _mm_add_epi32(s0, _mm_unpacklo_epi16(s4, _mm_setzero_si128()));
                            s0  = _mm_add_epi32(s0, _mm_unpacklo_epi16(s5, _mm_setzero_si128()));
                            s0  = _mm_add_epi32(s0, _mm_unpacklo_epi16(s6, _mm_setzero_si128()));
                            s1  = _mm_add_epi32(s1, _mm_unpackhi_epi16(s3, _mm_setzero_si128()));
                            s1  = _mm_add_epi32(s1, _mm_unpackhi_epi16(s4, _mm_setzero_si128()));
                            s1  = _mm_add_epi32(s1, _mm_unpackhi_epi16(s5, _mm_setzero_si128()));
                            s1  = _mm_add_epi32(s1, _mm_unpackhi_epi16(s6, _mm_setzero_si128()));
                            UPDATE_BEST(s0, 0, 0);
                            UPDATE_BEST(s0, 1, 0);
                            UPDATE_BEST(s0, 2, 0);
                            UPDATE_BEST(s0, 3, 0);
                            UPDATE_BEST(s1, 0, 4);
                            UPDATE_BEST(s1, 1, 4);
                            UPDATE_BEST(s1, 2, 4);
                            UPDATE_BEST(s1, 3, 4);
                        }
                    }
                }

                if (leftover && j < search_area_width) {
                    p_src = src;
                    p_ref = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        s0  = _mm_loadu_si128((__m128i *)(p_ref + 32));
                        s1  = _mm_loadu_si128((__m128i *)(p_ref + 40));
                        s2  = _mm_loadu_si128((__m128i *)(p_src + 32));
                        s3  = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4  = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5  = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6  = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    s0        = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    ss7       = _mm256_adds_epu16(ss3, ss4);
                    ss8       = _mm256_adds_epu16(ss5, ss6);
                    ss7       = _mm256_adds_epu16(ss7, ss8);
                    s0        = _mm_adds_epu16(s0,
                                        _mm_adds_epu16(_mm256_castsi256_si128(ss7),
                                                       _mm256_extracti128_si256(ss7, 1)));
                    s0        = _mm_or_si128(s0, s8);
                    s0        = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            ss0 = _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256());
                            ss3 = _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256());
                            ss1 = _mm256_unpacklo_epi16(ss4, _mm256_setzero_si256());
                            ss4 = _mm256_unpackhi_epi16(ss4, _mm256_setzero_si256());
                            ss2 = _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256());
                            ss5 = _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256());
                            ss7 = _mm256_unpacklo_epi16(ss6, _mm256_setzero_si256());
                            ss6 = _mm256_unpackhi_epi16(ss6, _mm256_setzero_si256());
                            ss0 = _mm256_add_epi32(_mm256_add_epi32(ss0, ss1),
                                                   _mm256_add_epi32(ss2, ss7));
                            ss3 = _mm256_add_epi32(_mm256_add_epi32(ss3, ss4),
                                                   _mm256_add_epi32(ss5, ss6));
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss0),
                                               _mm256_extracti128_si256(ss0, 1));
                            s1  = _mm_add_epi32(_mm256_castsi256_si128(ss3),
                                               _mm256_extracti128_si256(ss3, 1));
                            s0  = _mm_add_epi32(s0, _mm_unpacklo_epi16(s3, _mm_setzero_si128()));
                            s0  = _mm_add_epi32(s0, _mm_unpacklo_epi16(s4, _mm_setzero_si128()));
                            s0  = _mm_add_epi32(s0, _mm_unpacklo_epi16(s5, _mm_setzero_si128()));
                            s0  = _mm_add_epi32(s0, _mm_unpacklo_epi16(s6, _mm_setzero_si128()));
                            s1  = _mm_add_epi32(s1, _mm_unpackhi_epi16(s3, _mm_setzero_si128()));
                            s1  = _mm_add_epi32(s1, _mm_unpackhi_epi16(s4, _mm_setzero_si128()));
                            s1  = _mm_add_epi32(s1, _mm_unpackhi_epi16(s5, _mm_setzero_si128()));
                            s1  = _mm_add_epi32(s1, _mm_unpackhi_epi16(s6, _mm_setzero_si128()));
                            k   = leftover;
                            while (k > 0) {
                                for (l = 0; l < 4 && k; l++, k--) {
                                    tem_sum_1 = _mm_extract_epi32(s0, 0);
                                    s0        = _mm_srli_si128(s0, 4);
                                    if (tem_sum_1 < low_sum) {
                                        low_sum = tem_sum_1;
                                        x_best  = (int16_t)(j + leftover - k);
                                        y_best  = i;
                                    }
                                }
                                s0 = s1;
                            }
                        }
                    }
                }
                ref += src_stride_raw;
            }
        }
        break;

    case 64:
        if (block_height <= 32) {
            for (i = 0; i < search_area_height; i++) {
                uint32_t start_w = (i & 1) << 3;
                for (j = start_w; j <= search_area_width - 8; j += 16) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        ss0 = _mm256_loadu_si256((__m256i *)(p_ref + 32));
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 40));
                        ss2 = _mm256_loadu_si256((__m256i *)(p_src + 32));
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    ss7 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    s3        = _mm256_castsi256_si128(ss7);
                    s4        = _mm256_extracti128_si256(ss7, 1);
                    s0        = _mm_adds_epu16(s3, s4);
                    s0        = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            ss0 = _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256());
                            ss3 = _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256());
                            ss1 = _mm256_unpacklo_epi16(ss4, _mm256_setzero_si256());
                            ss4 = _mm256_unpackhi_epi16(ss4, _mm256_setzero_si256());
                            ss2 = _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256());
                            ss5 = _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256());
                            ss7 = _mm256_unpacklo_epi16(ss6, _mm256_setzero_si256());
                            ss6 = _mm256_unpackhi_epi16(ss6, _mm256_setzero_si256());
                            ss0 = _mm256_add_epi32(_mm256_add_epi32(ss0, ss1),
                                                   _mm256_add_epi32(ss2, ss7));
                            ss3 = _mm256_add_epi32(_mm256_add_epi32(ss3, ss4),
                                                   _mm256_add_epi32(ss5, ss6));
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss0),
                                               _mm256_extracti128_si256(ss0, 1));
                            s3  = _mm_add_epi32(_mm256_castsi256_si128(ss3),
                                               _mm256_extracti128_si256(ss3, 1));
                            UPDATE_BEST(s0, 0, 0);
                            UPDATE_BEST(s0, 1, 0);
                            UPDATE_BEST(s0, 2, 0);
                            UPDATE_BEST(s0, 3, 0);
                            UPDATE_BEST(s3, 0, 4);
                            UPDATE_BEST(s3, 1, 4);
                            UPDATE_BEST(s3, 2, 4);
                            UPDATE_BEST(s3, 3, 4);
                        }
                    }
                }

                if (leftover && j < search_area_width) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        ss0 = _mm256_loadu_si256((__m256i *)(p_ref + 32));
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 40));
                        ss2 = _mm256_loadu_si256((__m256i *)(p_src + 32));
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    ss7 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    s3        = _mm256_castsi256_si128(ss7);
                    s4        = _mm256_extracti128_si256(ss7, 1);
                    s0        = _mm_adds_epu16(s3, s4);
                    s0        = _mm_or_si128(s0, s8);
                    s0        = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            ss0 = _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256());
                            ss3 = _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256());
                            ss1 = _mm256_unpacklo_epi16(ss4, _mm256_setzero_si256());
                            ss4 = _mm256_unpackhi_epi16(ss4, _mm256_setzero_si256());
                            ss2 = _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256());
                            ss5 = _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256());
                            ss7 = _mm256_unpacklo_epi16(ss6, _mm256_setzero_si256());
                            ss6 = _mm256_unpackhi_epi16(ss6, _mm256_setzero_si256());
                            ss0 = _mm256_add_epi32(_mm256_add_epi32(ss0, ss1),
                                                   _mm256_add_epi32(ss2, ss7));
                            ss3 = _mm256_add_epi32(_mm256_add_epi32(ss3, ss4),
                                                   _mm256_add_epi32(ss5, ss6));
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss0),
                                               _mm256_extracti128_si256(ss0, 1));
                            s3  = _mm_add_epi32(_mm256_castsi256_si128(ss3),
                                               _mm256_extracti128_si256(ss3, 1));
                            k   = leftover;
                            while (k > 0) {
                                for (l = 0; l < 4 && k; l++, k--) {
                                    tem_sum_1 = _mm_extract_epi32(s0, 0);
                                    s0        = _mm_srli_si128(s0, 4);
                                    if (tem_sum_1 < low_sum) {
                                        low_sum = tem_sum_1;
                                        x_best  = (int16_t)(j + leftover - k);
                                        y_best  = i;
                                    }
                                }
                                s0 = s3;
                            }
                        }
                    }
                }
                ref += src_stride_raw;
            }
        } else {
            __m256i ss9, ss10;
            for (i = 0; i < search_area_height; i++) {
                uint32_t start_w = (i & 1) << 3;
                for (j = start_w; j <= search_area_width - 8; j += 16) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = ss7 = ss8 = ss9 = ss10 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        ss0 = _mm256_loadu_si256((__m256i *)(p_ref + 32));
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 40));
                        ss2 = _mm256_loadu_si256((__m256i *)(p_src + 32));
                        ss7 = _mm256_adds_epu16(ss7, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss8 = _mm256_adds_epu16(ss8, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss9 = _mm256_adds_epu16(ss9, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss10 =
                            _mm256_adds_epu16(ss10, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    ss0 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    ss0       = _mm256_adds_epu16(ss0,
                                            _mm256_adds_epu16(_mm256_adds_epu16(ss7, ss8),
                                                              _mm256_adds_epu16(ss9, ss10)));
                    s0        = _mm_adds_epu16(_mm256_castsi256_si128(ss0),
                                        _mm256_extracti128_si256(ss0, 1));
                    s0        = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            ss0 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss4, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss6, _mm256_setzero_si256())));
                            ss1 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss4, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss6, _mm256_setzero_si256())));
                            ss2 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss7, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss8, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss9, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss10, _mm256_setzero_si256())));
                            ss3 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss7, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss8, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss9, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss10, _mm256_setzero_si256())));
                            ss0 = _mm256_add_epi32(ss0, ss2);
                            ss1 = _mm256_add_epi32(ss1, ss3);
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss0),
                                               _mm256_extracti128_si256(ss0, 1));
                            s3  = _mm_add_epi32(_mm256_castsi256_si128(ss1),
                                               _mm256_extracti128_si256(ss1, 1));
                            UPDATE_BEST(s0, 0, 0);
                            UPDATE_BEST(s0, 1, 0);
                            UPDATE_BEST(s0, 2, 0);
                            UPDATE_BEST(s0, 3, 0);
                            UPDATE_BEST(s3, 0, 4);
                            UPDATE_BEST(s3, 1, 4);
                            UPDATE_BEST(s3, 2, 4);
                            UPDATE_BEST(s3, 3, 4);
                        }
                    }
                }

                if (leftover && j < search_area_width) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = ss7 = ss8 = ss9 = ss10 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        ss0 = _mm256_loadu_si256((__m256i *)(p_ref + 32));
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 40));
                        ss2 = _mm256_loadu_si256((__m256i *)(p_src + 32));
                        ss7 = _mm256_adds_epu16(ss7, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss8 = _mm256_adds_epu16(ss8, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss9 = _mm256_adds_epu16(ss9, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss10 =
                            _mm256_adds_epu16(ss10, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    ss0 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    ss0       = _mm256_adds_epu16(ss0,
                                            _mm256_adds_epu16(_mm256_adds_epu16(ss7, ss8),
                                                              _mm256_adds_epu16(ss9, ss10)));
                    s0        = _mm_adds_epu16(_mm256_castsi256_si128(ss0),
                                        _mm256_extracti128_si256(ss0, 1));
                    s0        = _mm_or_si128(s0, s8);
                    s0        = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            ss0 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss4, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss6, _mm256_setzero_si256())));
                            ss1 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss4, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss6, _mm256_setzero_si256())));
                            ss2 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss7, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss8, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss9, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss10, _mm256_setzero_si256())));
                            ss3 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss7, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss8, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss9, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss10, _mm256_setzero_si256())));
                            ss0 = _mm256_add_epi32(ss0, ss2);
                            ss1 = _mm256_add_epi32(ss1, ss3);
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss0),
                                               _mm256_extracti128_si256(ss0, 1));
                            s3  = _mm_add_epi32(_mm256_castsi256_si128(ss1),
                                               _mm256_extracti128_si256(ss1, 1));
                            k   = leftover;
                            while (k > 0) {
                                for (l = 0; l < 4 && k; l++, k--) {
                                    tem_sum_1 = _mm_extract_epi32(s0, 0);
                                    s0        = _mm_srli_si128(s0, 4);
                                    if (tem_sum_1 < low_sum) {
                                        low_sum = tem_sum_1;
                                        x_best  = (int16_t)(j + leftover - k);
                                        y_best  = i;
                                    }
                                }
                                s0 = s3;
                            }
                        }
                    }
                }
                ref += src_stride_raw;
            }
        }
        break;

    default: assert(0); break;
    }

    *best_sad        = low_sum;
    *x_search_center = x_best;
    *y_search_center = y_best;
}
#endif
/*******************************************************************************
 * Requirement: width   = 4, 8, 16, 24, 32, 48 or 64
 * Requirement: block_height <= 64
 * Requirement: block_height % 2 = 0 when width = 4 or 8
*******************************************************************************/
void sad_loop_kernel_avx2_intrin(
    uint8_t * src, // input parameter, source samples Ptr
    uint32_t  src_stride, // input parameter, source stride
    uint8_t * ref, // input parameter, reference samples Ptr
    uint32_t  ref_stride, // input parameter, reference stride
    uint32_t  block_height, // input parameter, block height (M)
    uint32_t  block_width, // input parameter, block width (N)
    uint64_t *best_sad, int16_t *x_search_center, int16_t *y_search_center,
    uint32_t src_stride_raw, // input parameter, source stride (no line skipping)
    int16_t search_area_width, int16_t search_area_height) {
    int16_t        x_best = *x_search_center, y_best = *y_search_center;
    uint32_t       low_sum   = 0xffffff;
    uint32_t       tem_sum_1 = 0;
    int16_t        i, j;
    uint32_t       k, l;
    uint32_t       leftover = search_area_width & 7;
    const uint8_t *p_ref, *p_src;
    __m128i        s0, s1, s2, s3, s4, s5, s6, s8 = _mm_set1_epi32(-1);
    __m256i        ss0, ss1, ss2, ss3, ss4, ss5, ss6, ss7, ss8;

    if (leftover) {
        for (k = 0; k < leftover; k++) s8 = _mm_slli_si128(s8, 2);
    }

    switch (block_width) {
    case 4:

        if (!(block_height % 4)) {
            uint32_t src_stride_t = 3 * src_stride;
            uint32_t ref_stride_t = 3 * ref_stride;
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss5 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k += 4) {
                        ss0 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_ref)),
                            _mm_loadu_si128((__m128i *)(p_ref + 2 * ref_stride)),
                            0x1);
                        ss1 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(
                                _mm_loadu_si128((__m128i *)(p_ref + ref_stride))),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride_t)),
                            0x1);
                        ss2 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_unpacklo_epi64(
                                _mm_cvtsi32_si128(*(uint32_t *)p_src),
                                _mm_cvtsi32_si128(*(uint32_t *)(p_src + src_stride)))),
                            _mm_unpacklo_epi64(
                                _mm_cvtsi32_si128(*(uint32_t *)(p_src + 2 * src_stride)),
                                _mm_cvtsi32_si128(*(uint32_t *)(p_src + src_stride_t))),
                            0x1);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        p_src += src_stride << 2;
                        p_ref += ref_stride << 2;
                    }
                    ss3       = _mm256_adds_epu16(ss3, ss5);
                    s3        = _mm_adds_epu16(_mm256_castsi256_si128(ss3),
                                        _mm256_extracti128_si256(ss3, 1));
                    s3        = _mm_minpos_epu16(s3);
                    tem_sum_1 = _mm_extract_epi16(s3, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s3, 1));
                        y_best  = i;
                    }
                }

                if (leftover) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss5 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k += 4) {
                        ss0 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_ref)),
                            _mm_loadu_si128((__m128i *)(p_ref + 2 * ref_stride)),
                            0x1);
                        ss1 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(
                                _mm_loadu_si128((__m128i *)(p_ref + ref_stride))),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride_t)),
                            0x1);
                        ss2 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_unpacklo_epi64(
                                _mm_cvtsi32_si128(*(uint32_t *)p_src),
                                _mm_cvtsi32_si128(*(uint32_t *)(p_src + src_stride)))),
                            _mm_unpacklo_epi64(
                                _mm_cvtsi32_si128(*(uint32_t *)(p_src + 2 * src_stride)),
                                _mm_cvtsi32_si128(*(uint32_t *)(p_src + src_stride_t))),
                            0x1);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        p_src += src_stride << 2;
                        p_ref += ref_stride << 2;
                    }
                    ss3       = _mm256_adds_epu16(ss3, ss5);
                    s3        = _mm_adds_epu16(_mm256_castsi256_si128(ss3),
                                        _mm256_extracti128_si256(ss3, 1));
                    s3        = _mm_or_si128(s3, s8);
                    s3        = _mm_minpos_epu16(s3);
                    tem_sum_1 = _mm_extract_epi16(s3, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s3, 1));
                        y_best  = i;
                    }
                }
                ref += src_stride_raw;
            }
        } else {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    s3    = _mm_setzero_si128();
                    for (k = 0; k < block_height; k += 2) {
                        s0 = _mm_loadu_si128((__m128i *)p_ref);
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + ref_stride));
                        s2 = _mm_cvtsi32_si128(*(uint32_t *)p_src);
                        s5 = _mm_cvtsi32_si128(*(uint32_t *)(p_src + src_stride));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s1, s5, 0));
                        p_src += src_stride << 1;
                        p_ref += ref_stride << 1;
                    }
                    s3        = _mm_minpos_epu16(s3);
                    tem_sum_1 = _mm_extract_epi16(s3, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s3, 1));
                        y_best  = i;
                    }
                }

                if (leftover) {
                    p_src = src;
                    p_ref = ref + j;
                    s3    = _mm_setzero_si128();
                    for (k = 0; k < block_height; k += 2) {
                        s0 = _mm_loadu_si128((__m128i *)p_ref);
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + ref_stride));
                        s2 = _mm_cvtsi32_si128(*(uint32_t *)p_src);
                        s5 = _mm_cvtsi32_si128(*(uint32_t *)(p_src + src_stride));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s1, s5, 0));
                        p_src += src_stride << 1;
                        p_ref += ref_stride << 1;
                    }
                    s3        = _mm_or_si128(s3, s8);
                    s3        = _mm_minpos_epu16(s3);
                    tem_sum_1 = _mm_extract_epi16(s3, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s3, 1));
                        y_best  = i;
                    }
                }
                ref += src_stride_raw;
            }
        }

        break;

    case 8:
        if (!(block_height % 4)) {
            uint32_t src_stride_t = 3 * src_stride;
            uint32_t ref_stride_t = 3 * ref_stride;
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k += 4) {
                        ss0 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_ref)),
                            _mm_loadu_si128((__m128i *)(p_ref + 2 * ref_stride)),
                            0x1);
                        ss1 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(
                                _mm_loadu_si128((__m128i *)(p_ref + ref_stride))),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride_t)),
                            0x1);
                        ss2 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_unpacklo_epi64(
                                _mm_loadl_epi64((__m128i *)p_src),
                                _mm_loadl_epi64((__m128i *)(p_src + src_stride)))),
                            _mm_unpacklo_epi64(_mm_loadl_epi64((__m128i *)(p_src + 2 * src_stride)),
                                               _mm_loadl_epi64((__m128i *)(p_src + src_stride_t))),
                            0x1);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride << 2;
                        p_ref += ref_stride << 2;
                    }
                    ss3 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    s3        = _mm_adds_epu16(_mm256_castsi256_si128(ss3),
                                        _mm256_extracti128_si256(ss3, 1));
                    s3        = _mm_minpos_epu16(s3);
                    tem_sum_1 = _mm_extract_epi16(s3, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s3, 1));
                        y_best  = i;
                    }
                }

                if (leftover) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k += 4) {
                        ss0 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_ref)),
                            _mm_loadu_si128((__m128i *)(p_ref + 2 * ref_stride)),
                            0x1);
                        ss1 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(
                                _mm_loadu_si128((__m128i *)(p_ref + ref_stride))),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride_t)),
                            0x1);
                        ss2 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_unpacklo_epi64(
                                _mm_loadl_epi64((__m128i *)p_src),
                                _mm_loadl_epi64((__m128i *)(p_src + src_stride)))),
                            _mm_unpacklo_epi64(_mm_loadl_epi64((__m128i *)(p_src + 2 * src_stride)),
                                               _mm_loadl_epi64((__m128i *)(p_src + src_stride_t))),
                            0x1);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride << 2;
                        p_ref += ref_stride << 2;
                    }
                    ss3 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    s3        = _mm_adds_epu16(_mm256_castsi256_si128(ss3),
                                        _mm256_extracti128_si256(ss3, 1));
                    s3        = _mm_or_si128(s3, s8);
                    s3        = _mm_minpos_epu16(s3);
                    tem_sum_1 = _mm_extract_epi16(s3, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s3, 1));
                        y_best  = i;
                    }
                }
                ref += src_stride_raw;
            }
        } else {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    s3 = s4 = _mm_setzero_si128();
                    for (k = 0; k < block_height; k += 2) {
                        s0 = _mm_loadu_si128((__m128i *)p_ref);
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + ref_stride));
                        s2 = _mm_loadl_epi64((__m128i *)p_src);
                        s5 = _mm_loadl_epi64((__m128i *)(p_src + src_stride));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s1, s5, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s1, s5, 5));
                        p_src += src_stride << 1;
                        p_ref += ref_stride << 1;
                    }
                    s3        = _mm_adds_epu16(s3, s4);
                    s3        = _mm_minpos_epu16(s3);
                    tem_sum_1 = _mm_extract_epi16(s3, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s3, 1));
                        y_best  = i;
                    }
                }

                if (leftover) {
                    p_src = src;
                    p_ref = ref + j;
                    s3 = s4 = _mm_setzero_si128();
                    for (k = 0; k < block_height; k += 2) {
                        s0 = _mm_loadu_si128((__m128i *)p_ref);
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + ref_stride));
                        s2 = _mm_loadl_epi64((__m128i *)p_src);
                        s5 = _mm_loadl_epi64((__m128i *)(p_src + src_stride));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s1, s5, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s1, s5, 5));
                        p_src += src_stride << 1;
                        p_ref += ref_stride << 1;
                    }
                    s3        = _mm_adds_epu16(s3, s4);
                    s3        = _mm_or_si128(s3, s8);
                    s3        = _mm_minpos_epu16(s3);
                    tem_sum_1 = _mm_extract_epi16(s3, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s3, 1));
                        y_best  = i;
                    }
                }
                ref += src_stride_raw;
            }
        }
        break;

    case 16:
        if (block_height <= 16) {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k += 2) {
                        ss0 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_ref)),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride)),
                            0x1);
                        ss1 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)(p_ref + 8))),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride + 8)),
                            0x1);
                        ss2 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_src)),
                            _mm_loadu_si128((__m128i *)(p_src + src_stride)),
                            0x1);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += 2 * src_stride;
                        p_ref += 2 * ref_stride;
                    }
                    ss3 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    s3        = _mm_adds_epu16(_mm256_castsi256_si128(ss3),
                                        _mm256_extracti128_si256(ss3, 1));
                    s3        = _mm_minpos_epu16(s3);
                    tem_sum_1 = _mm_extract_epi16(s3, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s3, 1));
                        y_best  = i;
                    }
                }

                if (leftover) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k += 2) {
                        ss0 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_ref)),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride)),
                            0x1);
                        ss1 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)(p_ref + 8))),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride + 8)),
                            0x1);
                        ss2 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_src)),
                            _mm_loadu_si128((__m128i *)(p_src + src_stride)),
                            0x1);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += 2 * src_stride;
                        p_ref += 2 * ref_stride;
                    }
                    ss3 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    s3        = _mm_adds_epu16(_mm256_castsi256_si128(ss3),
                                        _mm256_extracti128_si256(ss3, 1));
                    s3        = _mm_or_si128(s3, s8);
                    s3        = _mm_minpos_epu16(s3);
                    tem_sum_1 = _mm_extract_epi16(s3, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s3, 1));
                        y_best  = i;
                    }
                }
                ref += src_stride_raw;
            }
        } else if (block_height <= 32) {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k += 2) {
                        ss0 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_ref)),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride)),
                            0x1);
                        ss1 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)(p_ref + 8))),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride + 8)),
                            0x1);
                        ss2 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_src)),
                            _mm_loadu_si128((__m128i *)(p_src + src_stride)),
                            0x1);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += 2 * src_stride;
                        p_ref += 2 * ref_stride;
                    }
                    ss3 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    s3        = _mm256_castsi256_si128(ss3);
                    s5        = _mm256_extracti128_si256(ss3, 1);
                    s4        = _mm_minpos_epu16(s3);
                    s6        = _mm_minpos_epu16(s5);
                    s4        = _mm_unpacklo_epi16(s4, s4);
                    s4        = _mm_unpacklo_epi32(s4, s4);
                    s4        = _mm_unpacklo_epi64(s4, s4);
                    s6        = _mm_unpacklo_epi16(s6, s6);
                    s6        = _mm_unpacklo_epi32(s6, s6);
                    s6        = _mm_unpacklo_epi64(s6, s6);
                    s3        = _mm_sub_epi16(s3, s4);
                    s5        = _mm_adds_epu16(s5, s3);
                    s5        = _mm_sub_epi16(s5, s6);
                    s5        = _mm_minpos_epu16(s5);
                    tem_sum_1 = _mm_extract_epi16(s5, 0);
                    tem_sum_1 += _mm_extract_epi16(s4, 0);
                    tem_sum_1 += _mm_extract_epi16(s6, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s5, 1));
                        y_best  = i;
                    }
                }

                if (leftover) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k += 2) {
                        ss0 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_ref)),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride)),
                            0x1);
                        ss1 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)(p_ref + 8))),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride + 8)),
                            0x1);
                        ss2 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_src)),
                            _mm_loadu_si128((__m128i *)(p_src + src_stride)),
                            0x1);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += 2 * src_stride;
                        p_ref += 2 * ref_stride;
                    }
                    ss3 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    s3        = _mm256_castsi256_si128(ss3);
                    s5        = _mm256_extracti128_si256(ss3, 1);
                    s3        = _mm_or_si128(s3, s8);
                    s5        = _mm_or_si128(s5, s8);
                    s4        = _mm_minpos_epu16(s3);
                    s6        = _mm_minpos_epu16(s5);
                    s4        = _mm_unpacklo_epi16(s4, s4);
                    s4        = _mm_unpacklo_epi32(s4, s4);
                    s4        = _mm_unpacklo_epi64(s4, s4);
                    s6        = _mm_unpacklo_epi16(s6, s6);
                    s6        = _mm_unpacklo_epi32(s6, s6);
                    s6        = _mm_unpacklo_epi64(s6, s6);
                    s3        = _mm_sub_epi16(s3, s4);
                    s5        = _mm_adds_epu16(s5, s3);
                    s5        = _mm_sub_epi16(s5, s6);
                    s5        = _mm_minpos_epu16(s5);
                    tem_sum_1 = _mm_extract_epi16(s5, 0);
                    tem_sum_1 += _mm_extract_epi16(s4, 0);
                    tem_sum_1 += _mm_extract_epi16(s6, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s5, 1));
                        y_best  = i;
                    }
                }
                ref += src_stride_raw;
            }
        } else {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k += 2) {
                        ss0 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_ref)),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride)),
                            0x1);
                        ss1 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)(p_ref + 8))),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride + 8)),
                            0x1);
                        ss2 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_src)),
                            _mm_loadu_si128((__m128i *)(p_src + src_stride)),
                            0x1);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += 2 * src_stride;
                        p_ref += 2 * ref_stride;
                    }
                    ss3       = _mm256_adds_epu16(ss3, ss4);
                    ss5       = _mm256_adds_epu16(ss5, ss6);
                    ss0       = _mm256_adds_epu16(ss3, ss5);
                    s0        = _mm_adds_epu16(_mm256_castsi256_si128(ss0),
                                        _mm256_extracti128_si256(ss0, 1));
                    s0        = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            ss4 = _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256());
                            ss3 = _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256());
                            ss6 = _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256());
                            ss5 = _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256());
                            ss4 = _mm256_add_epi32(ss4, ss6);
                            ss3 = _mm256_add_epi32(ss3, ss5);
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss4),
                                               _mm256_extracti128_si256(ss4, 1));
                            s3  = _mm_add_epi32(_mm256_castsi256_si128(ss3),
                                               _mm256_extracti128_si256(ss3, 1));
                            UPDATE_BEST(s0, 0, 0);
                            UPDATE_BEST(s0, 1, 0);
                            UPDATE_BEST(s0, 2, 0);
                            UPDATE_BEST(s0, 3, 0);
                            UPDATE_BEST(s3, 0, 4);
                            UPDATE_BEST(s3, 1, 4);
                            UPDATE_BEST(s3, 2, 4);
                            UPDATE_BEST(s3, 3, 4);
                        }
                    }
                }

                if (leftover) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k += 2) {
                        ss0 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_ref)),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride)),
                            0x1);
                        ss1 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)(p_ref + 8))),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride + 8)),
                            0x1);
                        ss2 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_src)),
                            _mm_loadu_si128((__m128i *)(p_src + src_stride)),
                            0x1);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += 2 * src_stride;
                        p_ref += 2 * ref_stride;
                    }
                    ss3       = _mm256_adds_epu16(ss3, ss4);
                    ss5       = _mm256_adds_epu16(ss5, ss6);
                    ss0       = _mm256_adds_epu16(ss3, ss5);
                    s0        = _mm_adds_epu16(_mm256_castsi256_si128(ss0),
                                        _mm256_extracti128_si256(ss0, 1));
                    s0        = _mm_or_si128(s0, s8);
                    s0        = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            ss4 = _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256());
                            ss3 = _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256());
                            ss6 = _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256());
                            ss5 = _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256());
                            ss4 = _mm256_add_epi32(ss4, ss6);
                            ss3 = _mm256_add_epi32(ss3, ss5);
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss4),
                                               _mm256_extracti128_si256(ss4, 1));
                            s3  = _mm_add_epi32(_mm256_castsi256_si128(ss3),
                                               _mm256_extracti128_si256(ss3, 1));
                            k   = leftover;
                            while (k > 0) {
                                for (l = 0; l < 4 && k; l++, k--) {
                                    tem_sum_1 = _mm_extract_epi32(s0, 0);
                                    s0        = _mm_srli_si128(s0, 4);
                                    if (tem_sum_1 < low_sum) {
                                        low_sum = tem_sum_1;
                                        x_best  = (int16_t)(j + leftover - k);
                                        y_best  = i;
                                    }
                                }
                                s0 = s3;
                            }
                        }
                    }
                }
                ref += src_stride_raw;
            }
        }
        break;

    case 24:
        if (block_height <= 16) {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    ss3       = _mm256_adds_epu16(ss3, ss4);
                    ss5       = _mm256_adds_epu16(ss5, ss6);
                    s3        = _mm_adds_epu16(_mm256_castsi256_si128(ss3),
                                        _mm256_extracti128_si256(ss3, 1));
                    s5        = _mm256_castsi256_si128(ss5);
                    s4        = _mm_minpos_epu16(s3);
                    s6        = _mm_minpos_epu16(s5);
                    s4        = _mm_unpacklo_epi16(s4, s4);
                    s4        = _mm_unpacklo_epi32(s4, s4);
                    s4        = _mm_unpacklo_epi64(s4, s4);
                    s6        = _mm_unpacklo_epi16(s6, s6);
                    s6        = _mm_unpacklo_epi32(s6, s6);
                    s6        = _mm_unpacklo_epi64(s6, s6);
                    s3        = _mm_sub_epi16(s3, s4);
                    s5        = _mm_adds_epu16(s5, s3);
                    s5        = _mm_sub_epi16(s5, s6);
                    s5        = _mm_minpos_epu16(s5);
                    tem_sum_1 = _mm_extract_epi16(s5, 0);
                    tem_sum_1 += _mm_extract_epi16(s4, 0);
                    tem_sum_1 += _mm_extract_epi16(s6, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s5, 1));
                        y_best  = i;
                    }
                }

                if (leftover) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    ss3       = _mm256_adds_epu16(ss3, ss4);
                    ss5       = _mm256_adds_epu16(ss5, ss6);
                    s3        = _mm_adds_epu16(_mm256_castsi256_si128(ss3),
                                        _mm256_extracti128_si256(ss3, 1));
                    s5        = _mm256_castsi256_si128(ss5);
                    s3        = _mm_or_si128(s3, s8);
                    s5        = _mm_or_si128(s5, s8);
                    s4        = _mm_minpos_epu16(s3);
                    s6        = _mm_minpos_epu16(s5);
                    s4        = _mm_unpacklo_epi16(s4, s4);
                    s4        = _mm_unpacklo_epi32(s4, s4);
                    s4        = _mm_unpacklo_epi64(s4, s4);
                    s6        = _mm_unpacklo_epi16(s6, s6);
                    s6        = _mm_unpacklo_epi32(s6, s6);
                    s6        = _mm_unpacklo_epi64(s6, s6);
                    s3        = _mm_sub_epi16(s3, s4);
                    s5        = _mm_adds_epu16(s5, s3);
                    s5        = _mm_sub_epi16(s5, s6);
                    s5        = _mm_minpos_epu16(s5);
                    tem_sum_1 = _mm_extract_epi16(s5, 0);
                    tem_sum_1 += _mm_extract_epi16(s4, 0);
                    tem_sum_1 += _mm_extract_epi16(s6, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s5, 1));
                        y_best  = i;
                    }
                }
                ref += src_stride_raw;
            }
        } else {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    ss3       = _mm256_adds_epu16(ss3, ss4);
                    ss5       = _mm256_adds_epu16(ss5, ss6);
                    s3        = _mm256_castsi256_si128(ss3);
                    s4        = _mm256_extracti128_si256(ss3, 1);
                    s5        = _mm256_castsi256_si128(ss5);
                    s0        = _mm_adds_epu16(_mm_adds_epu16(s3, s4), s5);
                    s0        = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            s0 = _mm_unpacklo_epi16(s3, _mm_setzero_si128());
                            s3 = _mm_unpackhi_epi16(s3, _mm_setzero_si128());
                            s1 = _mm_unpacklo_epi16(s4, _mm_setzero_si128());
                            s4 = _mm_unpackhi_epi16(s4, _mm_setzero_si128());
                            s2 = _mm_unpacklo_epi16(s5, _mm_setzero_si128());
                            s5 = _mm_unpackhi_epi16(s5, _mm_setzero_si128());
                            s0 = _mm_add_epi32(_mm_add_epi32(s0, s1), s2);
                            s3 = _mm_add_epi32(_mm_add_epi32(s3, s4), s5);
                            UPDATE_BEST(s0, 0, 0);
                            UPDATE_BEST(s0, 1, 0);
                            UPDATE_BEST(s0, 2, 0);
                            UPDATE_BEST(s0, 3, 0);
                            UPDATE_BEST(s3, 0, 4);
                            UPDATE_BEST(s3, 1, 4);
                            UPDATE_BEST(s3, 2, 4);
                            UPDATE_BEST(s3, 3, 4);
                        }
                    }
                }

                if (leftover) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    ss3       = _mm256_adds_epu16(ss3, ss4);
                    ss5       = _mm256_adds_epu16(ss5, ss6);
                    s3        = _mm256_castsi256_si128(ss3);
                    s4        = _mm256_extracti128_si256(ss3, 1);
                    s5        = _mm256_castsi256_si128(ss5);
                    s0        = _mm_adds_epu16(_mm_adds_epu16(s3, s4), s5);
                    s0        = _mm_or_si128(s0, s8);
                    s0        = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            s0 = _mm_unpacklo_epi16(s3, _mm_setzero_si128());
                            s3 = _mm_unpackhi_epi16(s3, _mm_setzero_si128());
                            s1 = _mm_unpacklo_epi16(s4, _mm_setzero_si128());
                            s4 = _mm_unpackhi_epi16(s4, _mm_setzero_si128());
                            s2 = _mm_unpacklo_epi16(s5, _mm_setzero_si128());
                            s5 = _mm_unpackhi_epi16(s5, _mm_setzero_si128());
                            s0 = _mm_add_epi32(_mm_add_epi32(s0, s1), s2);
                            s3 = _mm_add_epi32(_mm_add_epi32(s3, s4), s5);
                            k  = leftover;
                            while (k > 0) {
                                for (l = 0; l < 4 && k; l++, k--) {
                                    tem_sum_1 = _mm_extract_epi32(s0, 0);
                                    s0        = _mm_srli_si128(s0, 4);
                                    if (tem_sum_1 < low_sum) {
                                        low_sum = tem_sum_1;
                                        x_best  = (int16_t)(j + leftover - k);
                                        y_best  = i;
                                    }
                                }
                                s0 = s3;
                            }
                        }
                    }
                }
                ref += src_stride_raw;
            }
        }
        break;

    case 32:
        if (block_height <= 16) {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    ss3 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    s3        = _mm256_castsi256_si128(ss3);
                    s5        = _mm256_extracti128_si256(ss3, 1);
                    s4        = _mm_minpos_epu16(s3);
                    s6        = _mm_minpos_epu16(s5);
                    s4        = _mm_unpacklo_epi16(s4, s4);
                    s4        = _mm_unpacklo_epi32(s4, s4);
                    s4        = _mm_unpacklo_epi64(s4, s4);
                    s6        = _mm_unpacklo_epi16(s6, s6);
                    s6        = _mm_unpacklo_epi32(s6, s6);
                    s6        = _mm_unpacklo_epi64(s6, s6);
                    s3        = _mm_sub_epi16(s3, s4);
                    s5        = _mm_adds_epu16(s5, s3);
                    s5        = _mm_sub_epi16(s5, s6);
                    s5        = _mm_minpos_epu16(s5);
                    tem_sum_1 = _mm_extract_epi16(s5, 0);
                    tem_sum_1 += _mm_extract_epi16(s4, 0);
                    tem_sum_1 += _mm_extract_epi16(s6, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s5, 1));
                        y_best  = i;
                    }
                }

                if (leftover) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    ss3 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    s3        = _mm256_castsi256_si128(ss3);
                    s5        = _mm256_extracti128_si256(ss3, 1);
                    s3        = _mm_or_si128(s3, s8);
                    s5        = _mm_or_si128(s5, s8);
                    s4        = _mm_minpos_epu16(s3);
                    s6        = _mm_minpos_epu16(s5);
                    s4        = _mm_unpacklo_epi16(s4, s4);
                    s4        = _mm_unpacklo_epi32(s4, s4);
                    s4        = _mm_unpacklo_epi64(s4, s4);
                    s6        = _mm_unpacklo_epi16(s6, s6);
                    s6        = _mm_unpacklo_epi32(s6, s6);
                    s6        = _mm_unpacklo_epi64(s6, s6);
                    s3        = _mm_sub_epi16(s3, s4);
                    s5        = _mm_adds_epu16(s5, s3);
                    s5        = _mm_sub_epi16(s5, s6);
                    s5        = _mm_minpos_epu16(s5);
                    tem_sum_1 = _mm_extract_epi16(s5, 0);
                    tem_sum_1 += _mm_extract_epi16(s4, 0);
                    tem_sum_1 += _mm_extract_epi16(s6, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s5, 1));
                        y_best  = i;
                    }
                }
                ref += src_stride_raw;
            }
        } else if (block_height <= 32) {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    ss3       = _mm256_adds_epu16(ss3, ss4);
                    ss5       = _mm256_adds_epu16(ss5, ss6);
                    ss6       = _mm256_adds_epu16(ss3, ss5);
                    s3        = _mm256_castsi256_si128(ss6);
                    s4        = _mm256_extracti128_si256(ss6, 1);
                    s0        = _mm_adds_epu16(s3, s4);
                    s0        = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            ss4 = _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256());
                            ss3 = _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256());
                            ss6 = _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256());
                            ss5 = _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256());
                            ss4 = _mm256_add_epi32(ss4, ss6);
                            ss3 = _mm256_add_epi32(ss3, ss5);
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss4),
                                               _mm256_extracti128_si256(ss4, 1));
                            s3  = _mm_add_epi32(_mm256_castsi256_si128(ss3),
                                               _mm256_extracti128_si256(ss3, 1));
                            UPDATE_BEST(s0, 0, 0);
                            UPDATE_BEST(s0, 1, 0);
                            UPDATE_BEST(s0, 2, 0);
                            UPDATE_BEST(s0, 3, 0);
                            UPDATE_BEST(s3, 0, 4);
                            UPDATE_BEST(s3, 1, 4);
                            UPDATE_BEST(s3, 2, 4);
                            UPDATE_BEST(s3, 3, 4);
                        }
                    }
                }

                if (leftover) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    ss3 = _mm256_adds_epu16(ss3, ss4);
                    ss5 = _mm256_adds_epu16(ss5, ss6);
                    ss6 = _mm256_adds_epu16(ss3, ss5);
                    s3  = _mm256_castsi256_si128(ss6);
                    s4  = _mm256_extracti128_si256(ss6, 1);
                    s0  = _mm_adds_epu16(s3, s4);
                    //s0 = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0        = _mm_or_si128(s0, s8);
                    s0        = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            ss4 = _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256());
                            ss3 = _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256());
                            ss6 = _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256());
                            ss5 = _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256());
                            ss4 = _mm256_add_epi32(ss4, ss6);
                            ss3 = _mm256_add_epi32(ss3, ss5);
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss4),
                                               _mm256_extracti128_si256(ss4, 1));
                            s3  = _mm_add_epi32(_mm256_castsi256_si128(ss3),
                                               _mm256_extracti128_si256(ss3, 1));
                            k   = leftover;
                            while (k > 0) {
                                for (l = 0; l < 4 && k; l++, k--) {
                                    tem_sum_1 = _mm_extract_epi32(s0, 0);
                                    s0        = _mm_srli_si128(s0, 4);
                                    if (tem_sum_1 < low_sum) {
                                        low_sum = tem_sum_1;
                                        x_best  = (int16_t)(j + leftover - k);
                                        y_best  = i;
                                    }
                                }
                                s0 = s3;
                            }
                        }
                    }
                }
                ref += src_stride_raw;
            }
        } else {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    ss7 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    s3        = _mm256_castsi256_si128(ss7);
                    s4        = _mm256_extracti128_si256(ss7, 1);
                    s0        = _mm_adds_epu16(s3, s4);
                    s0        = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            ss0 = _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256());
                            ss3 = _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256());
                            ss1 = _mm256_unpacklo_epi16(ss4, _mm256_setzero_si256());
                            ss4 = _mm256_unpackhi_epi16(ss4, _mm256_setzero_si256());
                            ss2 = _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256());
                            ss5 = _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256());
                            ss7 = _mm256_unpacklo_epi16(ss6, _mm256_setzero_si256());
                            ss6 = _mm256_unpackhi_epi16(ss6, _mm256_setzero_si256());
                            ss0 = _mm256_add_epi32(_mm256_add_epi32(ss0, ss1),
                                                   _mm256_add_epi32(ss2, ss7));
                            ss3 = _mm256_add_epi32(_mm256_add_epi32(ss3, ss4),
                                                   _mm256_add_epi32(ss5, ss6));
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss0),
                                               _mm256_extracti128_si256(ss0, 1));
                            s3  = _mm_add_epi32(_mm256_castsi256_si128(ss3),
                                               _mm256_extracti128_si256(ss3, 1));
                            UPDATE_BEST(s0, 0, 0);
                            UPDATE_BEST(s0, 1, 0);
                            UPDATE_BEST(s0, 2, 0);
                            UPDATE_BEST(s0, 3, 0);
                            UPDATE_BEST(s3, 0, 4);
                            UPDATE_BEST(s3, 1, 4);
                            UPDATE_BEST(s3, 2, 4);
                            UPDATE_BEST(s3, 3, 4);
                        }
                    }
                }

                if (leftover) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    ss7 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    s3        = _mm256_castsi256_si128(ss7);
                    s4        = _mm256_extracti128_si256(ss7, 1);
                    s0        = _mm_adds_epu16(s3, s4);
                    s0        = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            ss0 = _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256());
                            ss3 = _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256());
                            ss1 = _mm256_unpacklo_epi16(ss4, _mm256_setzero_si256());
                            ss4 = _mm256_unpackhi_epi16(ss4, _mm256_setzero_si256());
                            ss2 = _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256());
                            ss5 = _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256());
                            ss7 = _mm256_unpacklo_epi16(ss6, _mm256_setzero_si256());
                            ss6 = _mm256_unpackhi_epi16(ss6, _mm256_setzero_si256());
                            ss0 = _mm256_add_epi32(_mm256_add_epi32(ss0, ss1),
                                                   _mm256_add_epi32(ss2, ss7));
                            ss3 = _mm256_add_epi32(_mm256_add_epi32(ss3, ss4),
                                                   _mm256_add_epi32(ss5, ss6));
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss0),
                                               _mm256_extracti128_si256(ss0, 1));
                            s3  = _mm_add_epi32(_mm256_castsi256_si128(ss3),
                                               _mm256_extracti128_si256(ss3, 1));
                            k   = leftover;
                            while (k > 0) {
                                for (l = 0; l < 4 && k; l++, k--) {
                                    tem_sum_1 = _mm_extract_epi32(s0, 0);
                                    s0        = _mm_srli_si128(s0, 4);
                                    if (tem_sum_1 < low_sum) {
                                        low_sum = tem_sum_1;
                                        x_best  = (int16_t)(j + leftover - k);
                                        y_best  = i;
                                    }
                                }
                                s0 = s3;
                            }
                        }
                    }
                }
                ref += src_stride_raw;
            }
        }
        break;

    case 48:
        if (block_height <= 32) {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        s0  = _mm_loadu_si128((__m128i *)(p_ref + 32));
                        s1  = _mm_loadu_si128((__m128i *)(p_ref + 40));
                        s2  = _mm_loadu_si128((__m128i *)(p_src + 32));
                        s3  = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4  = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5  = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6  = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    s3        = _mm_adds_epu16(s3, s4);
                    s5        = _mm_adds_epu16(s5, s6);
                    s0        = _mm_adds_epu16(s3, s5);
                    ss3       = _mm256_adds_epu16(ss3, ss4);
                    ss5       = _mm256_adds_epu16(ss5, ss6);
                    ss6       = _mm256_adds_epu16(ss3, ss5);
                    s0        = _mm_adds_epu16(s0,
                                        _mm_adds_epu16(_mm256_castsi256_si128(ss6),
                                                       _mm256_extracti128_si256(ss6, 1)));
                    s0        = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            ss4 = _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256());
                            ss3 = _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256());
                            ss6 = _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256());
                            ss5 = _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256());
                            ss4 = _mm256_add_epi32(ss4, ss6);
                            ss3 = _mm256_add_epi32(ss3, ss5);
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss4),
                                               _mm256_extracti128_si256(ss4, 1));
                            s1  = _mm_add_epi32(_mm256_castsi256_si128(ss3),
                                               _mm256_extracti128_si256(ss3, 1));
                            s0  = _mm_add_epi32(s0, _mm_unpacklo_epi16(s3, _mm_setzero_si128()));
                            s0  = _mm_add_epi32(s0, _mm_unpacklo_epi16(s5, _mm_setzero_si128()));
                            s1  = _mm_add_epi32(s1, _mm_unpackhi_epi16(s3, _mm_setzero_si128()));
                            s1  = _mm_add_epi32(s1, _mm_unpackhi_epi16(s5, _mm_setzero_si128()));
                            UPDATE_BEST(s0, 0, 0);
                            UPDATE_BEST(s0, 1, 0);
                            UPDATE_BEST(s0, 2, 0);
                            UPDATE_BEST(s0, 3, 0);
                            UPDATE_BEST(s1, 0, 4);
                            UPDATE_BEST(s1, 1, 4);
                            UPDATE_BEST(s1, 2, 4);
                            UPDATE_BEST(s1, 3, 4);
                        }
                    }
                }

                if (leftover) {
                    p_src = src;
                    p_ref = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        s0  = _mm_loadu_si128((__m128i *)(p_ref + 32));
                        s1  = _mm_loadu_si128((__m128i *)(p_ref + 40));
                        s2  = _mm_loadu_si128((__m128i *)(p_src + 32));
                        s3  = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4  = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5  = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6  = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    s3        = _mm_adds_epu16(s3, s4);
                    s5        = _mm_adds_epu16(s5, s6);
                    s0        = _mm_adds_epu16(s3, s5);
                    ss3       = _mm256_adds_epu16(ss3, ss4);
                    ss5       = _mm256_adds_epu16(ss5, ss6);
                    ss6       = _mm256_adds_epu16(ss3, ss5);
                    s0        = _mm_adds_epu16(s0,
                                        _mm_adds_epu16(_mm256_castsi256_si128(ss6),
                                                       _mm256_extracti128_si256(ss6, 1)));
                    s0        = _mm_or_si128(s0, s8);
                    s0        = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            ss4 = _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256());
                            ss3 = _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256());
                            ss6 = _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256());
                            ss5 = _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256());
                            ss4 = _mm256_add_epi32(ss4, ss6);
                            ss3 = _mm256_add_epi32(ss3, ss5);
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss4),
                                               _mm256_extracti128_si256(ss4, 1));
                            s1  = _mm_add_epi32(_mm256_castsi256_si128(ss3),
                                               _mm256_extracti128_si256(ss3, 1));
                            s0  = _mm_add_epi32(s0, _mm_unpacklo_epi16(s3, _mm_setzero_si128()));
                            s0  = _mm_add_epi32(s0, _mm_unpacklo_epi16(s5, _mm_setzero_si128()));
                            s1  = _mm_add_epi32(s1, _mm_unpackhi_epi16(s3, _mm_setzero_si128()));
                            s1  = _mm_add_epi32(s1, _mm_unpackhi_epi16(s5, _mm_setzero_si128()));
                            k   = leftover;
                            while (k > 0) {
                                for (l = 0; l < 4 && k; l++, k--) {
                                    tem_sum_1 = _mm_extract_epi32(s0, 0);
                                    s0        = _mm_srli_si128(s0, 4);
                                    if (tem_sum_1 < low_sum) {
                                        low_sum = tem_sum_1;
                                        x_best  = (int16_t)(j + leftover - k);
                                        y_best  = i;
                                    }
                                }
                                s0 = s1;
                            }
                        }
                    }
                }
                ref += src_stride_raw;
            }
        } else {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        s0  = _mm_loadu_si128((__m128i *)(p_ref + 32));
                        s1  = _mm_loadu_si128((__m128i *)(p_ref + 40));
                        s2  = _mm_loadu_si128((__m128i *)(p_src + 32));
                        s3  = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4  = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5  = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6  = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    s0        = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    ss7       = _mm256_adds_epu16(ss3, ss4);
                    ss8       = _mm256_adds_epu16(ss5, ss6);
                    ss7       = _mm256_adds_epu16(ss7, ss8);
                    s0        = _mm_adds_epu16(s0,
                                        _mm_adds_epu16(_mm256_castsi256_si128(ss7),
                                                       _mm256_extracti128_si256(ss7, 1)));
                    s0        = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            ss0 = _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256());
                            ss3 = _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256());
                            ss1 = _mm256_unpacklo_epi16(ss4, _mm256_setzero_si256());
                            ss4 = _mm256_unpackhi_epi16(ss4, _mm256_setzero_si256());
                            ss2 = _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256());
                            ss5 = _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256());
                            ss7 = _mm256_unpacklo_epi16(ss6, _mm256_setzero_si256());
                            ss6 = _mm256_unpackhi_epi16(ss6, _mm256_setzero_si256());
                            ss0 = _mm256_add_epi32(_mm256_add_epi32(ss0, ss1),
                                                   _mm256_add_epi32(ss2, ss7));
                            ss3 = _mm256_add_epi32(_mm256_add_epi32(ss3, ss4),
                                                   _mm256_add_epi32(ss5, ss6));
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss0),
                                               _mm256_extracti128_si256(ss0, 1));
                            s1  = _mm_add_epi32(_mm256_castsi256_si128(ss3),
                                               _mm256_extracti128_si256(ss3, 1));
                            s0  = _mm_add_epi32(s0, _mm_unpacklo_epi16(s3, _mm_setzero_si128()));
                            s0  = _mm_add_epi32(s0, _mm_unpacklo_epi16(s4, _mm_setzero_si128()));
                            s0  = _mm_add_epi32(s0, _mm_unpacklo_epi16(s5, _mm_setzero_si128()));
                            s0  = _mm_add_epi32(s0, _mm_unpacklo_epi16(s6, _mm_setzero_si128()));
                            s1  = _mm_add_epi32(s1, _mm_unpackhi_epi16(s3, _mm_setzero_si128()));
                            s1  = _mm_add_epi32(s1, _mm_unpackhi_epi16(s4, _mm_setzero_si128()));
                            s1  = _mm_add_epi32(s1, _mm_unpackhi_epi16(s5, _mm_setzero_si128()));
                            s1  = _mm_add_epi32(s1, _mm_unpackhi_epi16(s6, _mm_setzero_si128()));
                            UPDATE_BEST(s0, 0, 0);
                            UPDATE_BEST(s0, 1, 0);
                            UPDATE_BEST(s0, 2, 0);
                            UPDATE_BEST(s0, 3, 0);
                            UPDATE_BEST(s1, 0, 4);
                            UPDATE_BEST(s1, 1, 4);
                            UPDATE_BEST(s1, 2, 4);
                            UPDATE_BEST(s1, 3, 4);
                        }
                    }
                }

                if (leftover) {
                    p_src = src;
                    p_ref = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        s0  = _mm_loadu_si128((__m128i *)(p_ref + 32));
                        s1  = _mm_loadu_si128((__m128i *)(p_ref + 40));
                        s2  = _mm_loadu_si128((__m128i *)(p_src + 32));
                        s3  = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4  = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5  = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6  = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    s0        = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    ss7       = _mm256_adds_epu16(ss3, ss4);
                    ss8       = _mm256_adds_epu16(ss5, ss6);
                    ss7       = _mm256_adds_epu16(ss7, ss8);
                    s0        = _mm_adds_epu16(s0,
                                        _mm_adds_epu16(_mm256_castsi256_si128(ss7),
                                                       _mm256_extracti128_si256(ss7, 1)));
                    s0        = _mm_or_si128(s0, s8);
                    s0        = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            ss0 = _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256());
                            ss3 = _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256());
                            ss1 = _mm256_unpacklo_epi16(ss4, _mm256_setzero_si256());
                            ss4 = _mm256_unpackhi_epi16(ss4, _mm256_setzero_si256());
                            ss2 = _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256());
                            ss5 = _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256());
                            ss7 = _mm256_unpacklo_epi16(ss6, _mm256_setzero_si256());
                            ss6 = _mm256_unpackhi_epi16(ss6, _mm256_setzero_si256());
                            ss0 = _mm256_add_epi32(_mm256_add_epi32(ss0, ss1),
                                                   _mm256_add_epi32(ss2, ss7));
                            ss3 = _mm256_add_epi32(_mm256_add_epi32(ss3, ss4),
                                                   _mm256_add_epi32(ss5, ss6));
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss0),
                                               _mm256_extracti128_si256(ss0, 1));
                            s1  = _mm_add_epi32(_mm256_castsi256_si128(ss3),
                                               _mm256_extracti128_si256(ss3, 1));
                            s0  = _mm_add_epi32(s0, _mm_unpacklo_epi16(s3, _mm_setzero_si128()));
                            s0  = _mm_add_epi32(s0, _mm_unpacklo_epi16(s4, _mm_setzero_si128()));
                            s0  = _mm_add_epi32(s0, _mm_unpacklo_epi16(s5, _mm_setzero_si128()));
                            s0  = _mm_add_epi32(s0, _mm_unpacklo_epi16(s6, _mm_setzero_si128()));
                            s1  = _mm_add_epi32(s1, _mm_unpackhi_epi16(s3, _mm_setzero_si128()));
                            s1  = _mm_add_epi32(s1, _mm_unpackhi_epi16(s4, _mm_setzero_si128()));
                            s1  = _mm_add_epi32(s1, _mm_unpackhi_epi16(s5, _mm_setzero_si128()));
                            s1  = _mm_add_epi32(s1, _mm_unpackhi_epi16(s6, _mm_setzero_si128()));
                            k   = leftover;
                            while (k > 0) {
                                for (l = 0; l < 4 && k; l++, k--) {
                                    tem_sum_1 = _mm_extract_epi32(s0, 0);
                                    s0        = _mm_srli_si128(s0, 4);
                                    if (tem_sum_1 < low_sum) {
                                        low_sum = tem_sum_1;
                                        x_best  = (int16_t)(j + leftover - k);
                                        y_best  = i;
                                    }
                                }
                                s0 = s1;
                            }
                        }
                    }
                }
                ref += src_stride_raw;
            }
        }
        break;

    case 64:
        if (block_height <= 32) {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        ss0 = _mm256_loadu_si256((__m256i *)(p_ref + 32));
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 40));
                        ss2 = _mm256_loadu_si256((__m256i *)(p_src + 32));
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    ss7 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    s3        = _mm256_castsi256_si128(ss7);
                    s4        = _mm256_extracti128_si256(ss7, 1);
                    s0        = _mm_adds_epu16(s3, s4);
                    s0        = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            ss0 = _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256());
                            ss3 = _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256());
                            ss1 = _mm256_unpacklo_epi16(ss4, _mm256_setzero_si256());
                            ss4 = _mm256_unpackhi_epi16(ss4, _mm256_setzero_si256());
                            ss2 = _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256());
                            ss5 = _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256());
                            ss7 = _mm256_unpacklo_epi16(ss6, _mm256_setzero_si256());
                            ss6 = _mm256_unpackhi_epi16(ss6, _mm256_setzero_si256());
                            ss0 = _mm256_add_epi32(_mm256_add_epi32(ss0, ss1),
                                                   _mm256_add_epi32(ss2, ss7));
                            ss3 = _mm256_add_epi32(_mm256_add_epi32(ss3, ss4),
                                                   _mm256_add_epi32(ss5, ss6));
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss0),
                                               _mm256_extracti128_si256(ss0, 1));
                            s3  = _mm_add_epi32(_mm256_castsi256_si128(ss3),
                                               _mm256_extracti128_si256(ss3, 1));
                            UPDATE_BEST(s0, 0, 0);
                            UPDATE_BEST(s0, 1, 0);
                            UPDATE_BEST(s0, 2, 0);
                            UPDATE_BEST(s0, 3, 0);
                            UPDATE_BEST(s3, 0, 4);
                            UPDATE_BEST(s3, 1, 4);
                            UPDATE_BEST(s3, 2, 4);
                            UPDATE_BEST(s3, 3, 4);
                        }
                    }
                }

                if (leftover) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        ss0 = _mm256_loadu_si256((__m256i *)(p_ref + 32));
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 40));
                        ss2 = _mm256_loadu_si256((__m256i *)(p_src + 32));
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    ss7 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    s3        = _mm256_castsi256_si128(ss7);
                    s4        = _mm256_extracti128_si256(ss7, 1);
                    s0        = _mm_adds_epu16(s3, s4);
                    s0        = _mm_or_si128(s0, s8);
                    s0        = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            ss0 = _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256());
                            ss3 = _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256());
                            ss1 = _mm256_unpacklo_epi16(ss4, _mm256_setzero_si256());
                            ss4 = _mm256_unpackhi_epi16(ss4, _mm256_setzero_si256());
                            ss2 = _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256());
                            ss5 = _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256());
                            ss7 = _mm256_unpacklo_epi16(ss6, _mm256_setzero_si256());
                            ss6 = _mm256_unpackhi_epi16(ss6, _mm256_setzero_si256());
                            ss0 = _mm256_add_epi32(_mm256_add_epi32(ss0, ss1),
                                                   _mm256_add_epi32(ss2, ss7));
                            ss3 = _mm256_add_epi32(_mm256_add_epi32(ss3, ss4),
                                                   _mm256_add_epi32(ss5, ss6));
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss0),
                                               _mm256_extracti128_si256(ss0, 1));
                            s3  = _mm_add_epi32(_mm256_castsi256_si128(ss3),
                                               _mm256_extracti128_si256(ss3, 1));
                            k   = leftover;
                            while (k > 0) {
                                for (l = 0; l < 4 && k; l++, k--) {
                                    tem_sum_1 = _mm_extract_epi32(s0, 0);
                                    s0        = _mm_srli_si128(s0, 4);
                                    if (tem_sum_1 < low_sum) {
                                        low_sum = tem_sum_1;
                                        x_best  = (int16_t)(j + leftover - k);
                                        y_best  = i;
                                    }
                                }
                                s0 = s3;
                            }
                        }
                    }
                }
                ref += src_stride_raw;
            }
        } else {
            __m256i ss9, ss10;
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = ss7 = ss8 = ss9 = ss10 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        ss0 = _mm256_loadu_si256((__m256i *)(p_ref + 32));
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 40));
                        ss2 = _mm256_loadu_si256((__m256i *)(p_src + 32));
                        ss7 = _mm256_adds_epu16(ss7, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss8 = _mm256_adds_epu16(ss8, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss9 = _mm256_adds_epu16(ss9, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss10 =
                            _mm256_adds_epu16(ss10, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    ss0 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    ss0       = _mm256_adds_epu16(ss0,
                                            _mm256_adds_epu16(_mm256_adds_epu16(ss7, ss8),
                                                              _mm256_adds_epu16(ss9, ss10)));
                    s0        = _mm_adds_epu16(_mm256_castsi256_si128(ss0),
                                        _mm256_extracti128_si256(ss0, 1));
                    s0        = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            ss0 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss4, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss6, _mm256_setzero_si256())));
                            ss1 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss4, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss6, _mm256_setzero_si256())));
                            ss2 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss7, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss8, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss9, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss10, _mm256_setzero_si256())));
                            ss3 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss7, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss8, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss9, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss10, _mm256_setzero_si256())));
                            ss0 = _mm256_add_epi32(ss0, ss2);
                            ss1 = _mm256_add_epi32(ss1, ss3);
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss0),
                                               _mm256_extracti128_si256(ss0, 1));
                            s3  = _mm_add_epi32(_mm256_castsi256_si128(ss1),
                                               _mm256_extracti128_si256(ss1, 1));
                            UPDATE_BEST(s0, 0, 0);
                            UPDATE_BEST(s0, 1, 0);
                            UPDATE_BEST(s0, 2, 0);
                            UPDATE_BEST(s0, 3, 0);
                            UPDATE_BEST(s3, 0, 4);
                            UPDATE_BEST(s3, 1, 4);
                            UPDATE_BEST(s3, 2, 4);
                            UPDATE_BEST(s3, 3, 4);
                        }
                    }
                }

                if (leftover) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = ss7 = ss8 = ss9 = ss10 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        ss0 = _mm256_loadu_si256((__m256i *)(p_ref + 32));
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 40));
                        ss2 = _mm256_loadu_si256((__m256i *)(p_src + 32));
                        ss7 = _mm256_adds_epu16(ss7, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss8 = _mm256_adds_epu16(ss8, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss9 = _mm256_adds_epu16(ss9, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss10 =
                            _mm256_adds_epu16(ss10, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    ss0 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    ss0       = _mm256_adds_epu16(ss0,
                                            _mm256_adds_epu16(_mm256_adds_epu16(ss7, ss8),
                                                              _mm256_adds_epu16(ss9, ss10)));
                    s0        = _mm_adds_epu16(_mm256_castsi256_si128(ss0),
                                        _mm256_extracti128_si256(ss0, 1));
                    s0        = _mm_or_si128(s0, s8);
                    s0        = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            ss0 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss4, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss6, _mm256_setzero_si256())));
                            ss1 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss4, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss6, _mm256_setzero_si256())));
                            ss2 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss7, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss8, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss9, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss10, _mm256_setzero_si256())));
                            ss3 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss7, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss8, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss9, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss10, _mm256_setzero_si256())));
                            ss0 = _mm256_add_epi32(ss0, ss2);
                            ss1 = _mm256_add_epi32(ss1, ss3);
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss0),
                                               _mm256_extracti128_si256(ss0, 1));
                            s3  = _mm_add_epi32(_mm256_castsi256_si128(ss1),
                                               _mm256_extracti128_si256(ss1, 1));
                            k   = leftover;
                            while (k > 0) {
                                for (l = 0; l < 4 && k; l++, k--) {
                                    tem_sum_1 = _mm_extract_epi32(s0, 0);
                                    s0        = _mm_srli_si128(s0, 4);
                                    if (tem_sum_1 < low_sum) {
                                        low_sum = tem_sum_1;
                                        x_best  = (int16_t)(j + leftover - k);
                                        y_best  = i;
                                    }
                                }
                                s0 = s3;
                            }
                        }
                    }
                }
                ref += src_stride_raw;
            }
        }
        break;

    default: assert(0); break;
    }

    *best_sad        = low_sum;
    *x_search_center = x_best;
    *y_search_center = y_best;
}

/*******************************************************************************
* Requirement: height % 4 = 0
*******************************************************************************/
SIMD_INLINE uint32_t
compute4x_m_sad_avx2(const uint8_t *src, // input parameter, source samples Ptr
                     uint32_t       src_stride, // input parameter, source stride
                     const uint8_t *ref, // input parameter, reference samples Ptr
                     uint32_t       ref_stride, // input parameter, reference stride
                     uint32_t       height) // input parameter, block height (M)
{
    uint32_t y = height;
    __m128i  xmm0;
    __m256i  ymm = _mm256_setzero_si256();

    do {
        const __m256i src0123 = load_u8_4x4_avx2(src, src_stride);
        const __m256i ref0123 = load_u8_4x4_avx2(ref, ref_stride);
        ymm                   = _mm256_add_epi32(ymm, _mm256_sad_epu8(src0123, ref0123));
        src += src_stride << 2;
        ref += ref_stride << 2;
        y -= 4;
    } while (y);

    xmm0 = _mm_add_epi32(_mm256_castsi256_si128(ymm), _mm256_extracti128_si256(ymm, 1));

    return (uint32_t)_mm_cvtsi128_si32(xmm0);
}

uint32_t eb_compute4x_m_sad_avx2_intrin(
    const uint8_t *src, // input parameter, source samples Ptr
    uint32_t       src_stride, // input parameter, source stride
    const uint8_t *ref, // input parameter, reference samples Ptr
    uint32_t       ref_stride, // input parameter, reference stride
    uint32_t       height, // input parameter, block height (M)
    uint32_t       width) // input parameter, block width (N)
{
    (void)width;
    return compute4x_m_sad_avx2(src, src_stride, ref, ref_stride, height);
}

static INLINE uint32_t sad_final_4_val_avx2(const __m256i sad) {
    const __m128i sad_lo = _mm256_castsi256_si128(sad);
    const __m128i sad_hi = _mm256_extracti128_si256(sad, 1);
    const __m128i sum    = _mm_add_epi32(sad_lo, sad_hi);
    const __m128i sum_hi = _mm_srli_si128(sum, 8);
    const __m128i d      = _mm_add_epi32(sum, sum_hi);
    return (uint32_t)_mm_cvtsi128_si32(d);
}

static INLINE uint32_t sad_final_8_val_avx2(const __m256i sad) {
    const __m128i sad_lo = _mm256_castsi256_si128(sad);
    const __m128i sad_hi = _mm256_extracti128_si256(sad, 1);
    const __m128i sum    = _mm_add_epi32(sad_lo, sad_hi);
    const __m128i sum_hi = _mm_srli_si128(sum, 8);
    const __m128i d      = _mm_add_epi32(sum, sum_hi);
    return (uint32_t)_mm_cvtsi128_si32(d) + _mm_extract_epi32(d, 1);
}

/*******************************************************************************
* Requirement: height % 4 = 0
*******************************************************************************/
SIMD_INLINE uint32_t
compute8x_m_sad_avx2(const uint8_t *src, // input parameter, source samples Ptr
                     uint32_t       src_stride, // input parameter, source stride
                     const uint8_t *ref, // input parameter, reference samples Ptr
                     uint32_t       ref_stride, // input parameter, reference stride
                     uint32_t       height) // input parameter, block height (M)
{
    uint32_t y   = height;
    __m256i  sad = _mm256_setzero_si256();

    do {
        const __m256i src0123 = load_u8_8x4_avx2(src, src_stride);
        const __m256i ref0123 = load_u8_8x4_avx2(ref, ref_stride);
        sad                   = _mm256_add_epi32(sad, _mm256_sad_epu8(src0123, ref0123));
        src += src_stride << 2;
        ref += ref_stride << 2;
        y -= 4;
    } while (y);

    return sad_final_4_val_avx2(sad);
}

uint32_t eb_compute8x_m_sad_avx2_intrin(
    const uint8_t *src, // input parameter, source samples Ptr
    uint32_t       src_stride, // input parameter, source stride
    const uint8_t *ref, // input parameter, reference samples Ptr
    uint32_t       ref_stride, // input parameter, reference stride
    uint32_t       height, // input parameter, block height (M)
    uint32_t       width) // input parameter, block width (N)
{
    (void)width;
    return compute8x_m_sad_avx2(src, src_stride, ref, ref_stride, height);
}

static INLINE void compute16_sad_avx2(const uint8_t *const src, const uint32_t src_stride,
                                      const uint8_t *const ref, const uint32_t ref_stride,
                                      __m256i *const sum) {
    const __m256i s   = loadu_u8_16x2_avx2(src, src_stride);
    const __m256i r   = loadu_u8_16x2_avx2(ref, ref_stride);
    const __m256i sad = _mm256_sad_epu8(s, r);
    *sum              = _mm256_add_epi32(*sum, sad);
}

static INLINE void compute32_sad_avx2(const uint8_t *const src, const uint8_t *const ref,
                                      __m256i *const sum) {
    const __m256i s   = _mm256_loadu_si256((__m256i *)src);
    const __m256i r   = _mm256_loadu_si256((__m256i *)ref);
    const __m256i sad = _mm256_sad_epu8(s, r);
    *sum              = _mm256_add_epi32(*sum, sad);
}

/*******************************************************************************
* Requirement: height % 4 = 0
*******************************************************************************/
SIMD_INLINE uint32_t
compute16x_m_sad_avx2(const uint8_t *src, // input parameter, source samples Ptr
                      uint32_t       src_stride, // input parameter, source stride
                      const uint8_t *ref, // input parameter, reference samples Ptr
                      uint32_t       ref_stride, // input parameter, reference stride
                      uint32_t       height) // input parameter, block height (M)
{
    uint32_t y   = height;
    __m256i  sad = _mm256_setzero_si256();

    do {
        compute16_sad_avx2(
            src + 0 * src_stride, src_stride, ref + 0 * ref_stride, ref_stride, &sad);
        compute16_sad_avx2(
            src + 2 * src_stride, src_stride, ref + 2 * ref_stride, ref_stride, &sad);
        src += src_stride << 2;
        ref += ref_stride << 2;
        y -= 4;
    } while (y);

    return sad_final_4_val_avx2(sad);
}

uint32_t eb_compute16x_m_sad_avx2_intrin(
    const uint8_t *src, // input parameter, source samples Ptr
    uint32_t       src_stride, // input parameter, source stride
    const uint8_t *ref, // input parameter, reference samples Ptr
    uint32_t       ref_stride, // input parameter, reference stride
    uint32_t       height, // input parameter, block height (M)
    uint32_t       width) // input parameter, block width (N)
{
    (void)width;
    return compute16x_m_sad_avx2(src, src_stride, ref, ref_stride, height);
}

/*******************************************************************************
* Requirement: height % 2 = 0
*******************************************************************************/
uint32_t eb_compute24x_m_sad_avx2_intrin(
    const uint8_t *src, // input parameter, source samples Ptr
    uint32_t       src_stride, // input parameter, source stride
    const uint8_t *ref, // input parameter, reference samples Ptr
    uint32_t       ref_stride, // input parameter, reference stride
    uint32_t       height, // input parameter, block height (M)
    uint32_t       width) // input parameter, block width (N)
{
    uint32_t y   = height;
    __m256i  sad = _mm256_setzero_si256();

    (void)width;

    do {
        compute32_sad_avx2(src + 0 * src_stride, ref + 0 * ref_stride, &sad);
        compute32_sad_avx2(src + 1 * src_stride, ref + 1 * ref_stride, &sad);
        src += src_stride << 1;
        ref += ref_stride << 1;
        y -= 2;
    } while (y);

    const __m128i sad_lo = _mm256_castsi256_si128(sad);
    const __m128i sad_hi = _mm256_extracti128_si256(sad, 1);
    const __m128i sum_lo = _mm_add_epi32(sad_lo, _mm_srli_si128(sad_lo, 8));
    const __m128i sum    = _mm_add_epi32(sum_lo, sad_hi);

    return (uint32_t)_mm_cvtsi128_si32(sum);
}

/*******************************************************************************
* Requirement: height % 2 = 0
*******************************************************************************/
static INLINE uint32_t
compute32x_m_sad_avx2(const uint8_t *src, // input parameter, source samples Ptr
                      uint32_t       src_stride, // input parameter, source stride
                      const uint8_t *ref, // input parameter, reference samples Ptr
                      uint32_t       ref_stride, // input parameter, reference stride
                      uint32_t       height) // input parameter, block height (M)
{
    uint32_t y   = height;
    __m256i  sad = _mm256_setzero_si256();

    do {
        compute32_sad_avx2(src + 0 * src_stride, ref + 0 * ref_stride, &sad);
        compute32_sad_avx2(src + 1 * src_stride, ref + 1 * ref_stride, &sad);
        src += src_stride << 1;
        ref += ref_stride << 1;
        y -= 2;
    } while (y);

    return sad_final_4_val_avx2(sad);
}

uint32_t eb_compute32x_m_sad_avx2_intrin(
    const uint8_t *src, // input parameter, source samples Ptr
    uint32_t       src_stride, // input parameter, source stride
    const uint8_t *ref, // input parameter, reference samples Ptr
    uint32_t       ref_stride, // input parameter, reference stride
    uint32_t       height, // input parameter, block height (M)
    uint32_t       width) // input parameter, block width (N)
{
    (void)width;
    return compute32x_m_sad_avx2(src, src_stride, ref, ref_stride, height);
}

/*******************************************************************************
* Requirement: height % 2 = 0
*******************************************************************************/
uint32_t eb_compute48x_m_sad_avx2_intrin(
    const uint8_t *src, // input parameter, source samples Ptr
    uint32_t       src_stride, // input parameter, source stride
    const uint8_t *ref, // input parameter, reference samples Ptr
    uint32_t       ref_stride, // input parameter, reference stride
    uint32_t       height, // input parameter, block height (M)
    uint32_t       width) // input parameter, block width (N)
{
    uint32_t y   = height;
    __m256i  sad = _mm256_setzero_si256();

    (void)width;

    do {
        compute32_sad_avx2(src + 0 * src_stride, ref + 0 * ref_stride, &sad);
        compute32_sad_avx2(src + 1 * src_stride, ref + 1 * ref_stride, &sad);
        compute16_sad_avx2(src + 32, src_stride, ref + 32, ref_stride, &sad);
        src += src_stride << 1;
        ref += ref_stride << 1;
        y -= 2;
    } while (y);

    return sad_final_4_val_avx2(sad);
}

/*******************************************************************************
* Requirement: height % 2 = 0
*******************************************************************************/
SIMD_INLINE uint32_t
compute64x_m_sad_avx2(const uint8_t *src, // input parameter, source samples Ptr
                      uint32_t       src_stride, // input parameter, source stride
                      const uint8_t *ref, // input parameter, reference samples Ptr
                      uint32_t       ref_stride, // input parameter, reference stride
                      uint32_t       height) // input parameter, block height (M)
{
    uint32_t y   = height;
    __m256i  sad = _mm256_setzero_si256();

    do {
        compute32_sad_avx2(src, ref, &sad);
        compute32_sad_avx2(src + 32, ref + 32, &sad);
        compute32_sad_avx2(src + src_stride, ref + ref_stride, &sad);
        compute32_sad_avx2(src + src_stride + 32, ref + ref_stride + 32, &sad);
        src += src_stride << 1;
        ref += ref_stride << 1;
        y -= 2;
    } while (y);

    return sad_final_4_val_avx2(sad);
}

SIMD_INLINE uint32_t
compute128x_m_sad_avx2(const uint8_t *src, // input parameter, source samples Ptr
                       uint32_t       src_stride, // input parameter, source stride
                       const uint8_t *ref, // input parameter, reference samples Ptr
                       uint32_t       ref_stride, // input parameter, reference stride
                       uint32_t       height) // input parameter, block height (M)
{
    uint32_t y   = height;
    __m256i  sad = _mm256_setzero_si256();

    do {
        compute32_sad_avx2(src + 0 * 32, ref + 0 * 32, &sad);
        compute32_sad_avx2(src + 1 * 32, ref + 1 * 32, &sad);
        compute32_sad_avx2(src + 2 * 32, ref + 2 * 32, &sad);
        compute32_sad_avx2(src + 3 * 32, ref + 3 * 32, &sad);
        src += src_stride;
        ref += ref_stride;
    } while (--y);

    return sad_final_4_val_avx2(sad);
}

uint32_t eb_compute64x_m_sad_avx2_intrin(
    const uint8_t *src, // input parameter, source samples Ptr
    uint32_t       src_stride, // input parameter, source stride
    const uint8_t *ref, // input parameter, reference samples Ptr
    uint32_t       ref_stride, // input parameter, reference stride
    uint32_t       height, // input parameter, block height (M)
    uint32_t       width) // input parameter, block width (N)
{
    (void)width;
    return compute64x_m_sad_avx2(src, src_stride, ref, ref_stride, height);
}

uint32_t eb_compute128x_m_sad_avx2_intrin(
    const uint8_t *src, // input parameter, source samples Ptr
    uint32_t       src_stride, // input parameter, source stride
    const uint8_t *ref, // input parameter, reference samples Ptr
    uint32_t       ref_stride, // input parameter, reference stride
    uint32_t       height, // input parameter, block height (M)
    uint32_t       width) // input parameter, block width (N)
{
    (void)width;
    return compute128x_m_sad_avx2(src, src_stride, ref, ref_stride, height);
}

#if !REMOVE_UNUSED_CODE
static INLINE void sad_eight_8x2x2_avx2(const uint8_t *src, const uint32_t src_stride,
                                        const uint8_t *ref, const uint32_t ref_stride,
                                        __m256i d[2]) {
    const __m256i s  = loadu_u8_16x2_avx2(src, src_stride);
    const __m256i r0 = loadu_u8_16x2_avx2(ref + 0, ref_stride);
    const __m256i r1 = loadu_u8_16x2_avx2(ref + 8, ref_stride);
    d[0]             = _mm256_adds_epu16(d[0], _mm256_mpsadbw_epu8(r0, s, 0x00));
    d[0]             = _mm256_adds_epu16(d[0], _mm256_mpsadbw_epu8(r0, s, 0x2D)); // 101 101
    d[1]             = _mm256_adds_epu16(d[1], _mm256_mpsadbw_epu8(r1, s, 0x12)); // 010 010
    d[1]             = _mm256_adds_epu16(d[1], _mm256_mpsadbw_epu8(r1, s, 0x3F)); // 111 111
}
#endif
/*
NOTE: Slower than avx2. Abandoned.
static INLINE void sad_eight_8x2x2_avx512(const uint8_t *const src,
    const uint32_t src_stride, const uint8_t *const ref,
    const uint32_t ref_stride, __m512i *const d)
{
    const __m512i c0 =
        _mm512_setr_epi32(0, 0, 0, 0, 2, 2, 2, 2, 4, 4, 4, 4, 6, 6, 6, 6);
    const __m512i c1 =
        _mm512_setr_epi32(1, 1, 1, 1, 3, 3, 3, 3, 5, 5, 5, 5, 7, 7, 7, 7);
    const __m512i c2 = _mm512_setr_epi64(0, 1, 1, 2, 4, 5, 5, 6);

    const __m256i s01 = loadu_u8_16x2_avx2(src, src_stride);
    const __m512i s = _mm512_castsi256_si512(s01);
    const __m512i ss0 = _mm512_permutexvar_epi32(c0, s);
    const __m512i ss1 = _mm512_permutexvar_epi32(c1, s);

    const __m256i r0 = _mm256_loadu_si256((__m256i*)ref);
    const __m256i r1 = _mm256_loadu_si256((__m256i*)(ref + ref_stride));
    const __m512i r = _mm512_inserti64x4(_mm512_castsi256_si512(r0), r1, 1);
    const __m512i rr0 = _mm512_permutexvar_epi64(c2, r);

    *d = _mm512_adds_epu16(*d, _mm512_dbsad_epu8(ss0, rr0, 0x94)); //10 01 01 00
    *d = _mm512_adds_epu16(*d, _mm512_dbsad_epu8(ss1, rr0, 0xE9)); //11 10 10 01
}

Usage:
s_512 = _mm512_setzero_si512();
sad_eight_8x2x2_avx512(src + 0 * src_stride, 8 * src_stride,
    ref + 0 * ref_stride, 8 * ref_stride, &s_512);
sad_eight_8x2x2_avx512(src + 2 * src_stride, 8 * src_stride,
     ref + 2 * ref_stride, 8 * ref_stride, &s_512);
sad_eight_8x2x2_avx512(src + 4 * src_stride, 8 * src_stride,
     ref + 4 * ref_stride, 8 * ref_stride, &s_512);
sad_eight_8x2x2_avx512(src + 6 * src_stride, 8 * src_stride,
     ref + 6 * ref_stride, 8 * ref_stride, &s_512);
sad_eight_8x2x2_avx512(src + 1 * src_stride, 8 * src_stride,
     ref + 1 * ref_stride, 8 * ref_stride, &s_512);
sad_eight_8x2x2_avx512(src + 3 * src_stride, 8 * src_stride,
     ref + 3 * ref_stride, 8 * ref_stride, &s_512);
sad_eight_8x2x2_avx512(src + 5 * src_stride, 8 * src_stride,
     ref + 5 * ref_stride, 8 * ref_stride, &s_512);
sad_eight_8x2x2_avx512(src + 7 * src_stride, 8 * src_stride,
     ref + 7 * ref_stride, 8 * ref_stride, &s_512);
s_256[0] = _mm512_castsi512_si256(s_512);
s_256[1] = _mm512_extracti64x4_epi64(s_512, 1);
*/

#if !REMOVE_UNUSED_CODE
/*******************************************************************************
* Requirement: p_best_sad_8x8[i] must be less than 0x7FFFFFFF because signed comparison is used.
*******************************************************************************/
void get_eight_horizontal_search_point_results_8x8_16x16_pu_avx2_intrin(
    uint8_t *src, uint32_t src_stride, uint8_t *ref, uint32_t ref_stride, uint32_t *p_best_sad_8x8,
    uint32_t *p_best_mv8x8, uint32_t *p_best_sad_16x16, uint32_t *p_best_mv16x16, uint32_t mv,
    uint16_t *p_sad16x16, EbBool sub_sad) {
    __m128i s_128[4], s16_128, best_sad_128, best_mv_128, mv_128;
    __m256i s_256[2];

    /*
    -------------------------------------   -----------------------------------
    | 8x8_00 | 8x8_01 | 8x8_04 | 8x8_05 |   8x8_16 | 8x8_17 | 8x8_20 | 8x8_21 |
    -------------------------------------   -----------------------------------
    | 8x8_02 | 8x8_03 | 8x8_06 | 8x8_07 |   8x8_18 | 8x8_19 | 8x8_22 | 8x8_23 |
    -----------------------   -----------   ----------------------   ----------
    | 8x8_08 | 8x8_09 | 8x8_12 | 8x8_13 |   8x8_24 | 8x8_25 | 8x8_29 | 8x8_29 |
    ----------------------    -----------   ---------------------    ----------
    | 8x8_10 | 8x8_11 | 8x8_14 | 8x8_15 |   8x8_26 | 8x8_27 | 8x8_30 | 8x8_31 |
    -------------------------------------   -----------------------------------

    -------------------------------------   -----------------------------------
    | 8x8_32 | 8x8_33 | 8x8_36 | 8x8_37 |   8x8_48 | 8x8_49 | 8x8_52 | 8x8_53 |
    -------------------------------------   -----------------------------------
    | 8x8_34 | 8x8_35 | 8x8_38 | 8x8_39 |   8x8_50 | 8x8_51 | 8x8_54 | 8x8_55 |
    -----------------------   -----------   ----------------------   ----------
    | 8x8_40 | 8x8_41 | 8x8_44 | 8x8_45 |   8x8_56 | 8x8_57 | 8x8_60 | 8x8_61 |
    ----------------------    -----------   ---------------------    ----------
    | 8x8_42 | 8x8_43 | 8x8_46 | 8x8_48 |   8x8_58 | 8x8_59 | 8x8_62 | 8x8_63 |
    -------------------------------------   -----------------------------------
    */

    /*
    ----------------------    ----------------------
    |  16x16_0  |  16x16_1  |  16x16_4  |  16x16_5  |
    ----------------------    ----------------------
    |  16x16_2  |  16x16_3  |  16x16_6  |  16x16_7  |
    -----------------------   -----------------------
    |  16x16_8  |  16x16_9  |  16x16_12 |  16x16_13 |
    ----------------------    ----------------------
    |  16x16_10 |  16x16_11 |  16x16_14 |  16x16_15 |
    -----------------------   -----------------------
    */

    // 8x8_0, 8x8_1, 8x8_2 & 8x8_3
    s_256[0] = s_256[1] = _mm256_setzero_si256();
    sad_eight_8x2x2_avx2(
        src + 0 * src_stride, 8 * src_stride, ref + 0 * ref_stride, 8 * ref_stride, s_256);
    sad_eight_8x2x2_avx2(
        src + 2 * src_stride, 8 * src_stride, ref + 2 * ref_stride, 8 * ref_stride, s_256);
    sad_eight_8x2x2_avx2(
        src + 4 * src_stride, 8 * src_stride, ref + 4 * ref_stride, 8 * ref_stride, s_256);
    sad_eight_8x2x2_avx2(
        src + 6 * src_stride, 8 * src_stride, ref + 6 * ref_stride, 8 * ref_stride, s_256);

    //16x16
    if (sub_sad) {
        s_256[0] = _mm256_slli_epi16(s_256[0], 1);
        s_256[1] = _mm256_slli_epi16(s_256[1], 1);
    } else {
        // 8x8_0, 8x8_1, 8x8_2 & 8x8_3
        sad_eight_8x2x2_avx2(
            src + 1 * src_stride, 8 * src_stride, ref + 1 * ref_stride, 8 * ref_stride, s_256);
        sad_eight_8x2x2_avx2(
            src + 3 * src_stride, 8 * src_stride, ref + 3 * ref_stride, 8 * ref_stride, s_256);
        sad_eight_8x2x2_avx2(
            src + 5 * src_stride, 8 * src_stride, ref + 5 * ref_stride, 8 * ref_stride, s_256);
        sad_eight_8x2x2_avx2(
            src + 7 * src_stride, 8 * src_stride, ref + 7 * ref_stride, 8 * ref_stride, s_256);
    }

    s_128[0] = _mm256_castsi256_si128(s_256[0]);
    s_128[1] = _mm256_castsi256_si128(s_256[1]);
    s_128[2] = _mm256_extracti128_si256(s_256[0], 1);
    s_128[3] = _mm256_extracti128_si256(s_256[1], 1);
    s16_128  = _mm_adds_epu16(s_128[0], s_128[1]);
    s16_128  = _mm_adds_epu16(s16_128, s_128[2]);
    s16_128  = _mm_adds_epu16(s16_128, s_128[3]);

    //sotore the 8 SADs(16x16 SADs)
    _mm_storeu_si128((__m128i *)p_sad16x16, s16_128);

    //find the best for 16x16
    const __m128i  minpos16 = _mm_minpos_epu16(s16_128);
    const uint32_t sad      = _mm_extract_epi16(minpos16, 0);

    if (sad < p_best_sad_16x16[0]) {
        const int32_t best_mv = 4 * _mm_extract_epi16(minpos16, 1);
        const int16_t x_mv    = _MVXT(mv) + (int16_t)best_mv;
        const int16_t y_mv    = _MVYT(mv);
        p_best_sad_16x16[0]   = sad;
        *p_best_mv16x16       = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
    }

    //find the best for 8x8_0, 8x8_1, 8x8_2 & 8x8_3
    s_128[0]                 = _mm_minpos_epu16(s_128[0]);
    s_128[1]                 = _mm_minpos_epu16(s_128[1]);
    s_128[2]                 = _mm_minpos_epu16(s_128[2]);
    s_128[3]                 = _mm_minpos_epu16(s_128[3]);
    s_128[0]                 = _mm_unpacklo_epi16(s_128[0], s_128[1]);
    s_128[2]                 = _mm_unpacklo_epi16(s_128[2], s_128[3]);
    s_128[0]                 = _mm_unpacklo_epi32(s_128[0], s_128[2]);
    s_128[1]                 = _mm_unpackhi_epi16(s_128[0], _mm_setzero_si128());
    const __m128i sad_128    = _mm_unpacklo_epi16(s_128[0], _mm_setzero_si128());
    const __m128i offset_128 = _mm_slli_epi16(s_128[1], 2);

    best_sad_128       = _mm_loadu_si128((__m128i *)p_best_sad_8x8);
    const __m128i mask = _mm_cmpgt_epi32(best_sad_128, sad_128);
    best_sad_128       = _mm_min_epu32(sad_128, best_sad_128);
    _mm_storeu_si128((__m128i *)p_best_sad_8x8, best_sad_128);

    mv_128      = _mm_set1_epi32(mv);
    mv_128      = _mm_add_epi16(mv_128, offset_128);
    mv_128      = _mm_and_si128(mv_128, mask);
    best_mv_128 = _mm_loadu_si128((__m128i *)p_best_mv8x8);
    best_mv_128 = _mm_andnot_si128(mask, best_mv_128);
    best_mv_128 = _mm_or_si128(best_mv_128, mv_128);
    _mm_storeu_si128((__m128i *)p_best_mv8x8, best_mv_128);
}
static INLINE __m256i add_4_sad_avx2(const uint16_t *const p_sad16x16) {
    const __m128i s0 = _mm_loadu_si128((__m128i *)(p_sad16x16 + 0 * 8));
    const __m128i s1 = _mm_loadu_si128((__m128i *)(p_sad16x16 + 1 * 8));
    const __m128i s2 = _mm_loadu_si128((__m128i *)(p_sad16x16 + 2 * 8));
    const __m128i s3 = _mm_loadu_si128((__m128i *)(p_sad16x16 + 3 * 8));
    __m256i       s_256[4];

    s_256[0] = _mm256_cvtepu16_epi32(s0);
    s_256[1] = _mm256_cvtepu16_epi32(s1);
    s_256[2] = _mm256_cvtepu16_epi32(s2);
    s_256[3] = _mm256_cvtepu16_epi32(s3);

    s_256[0] = _mm256_add_epi32(s_256[0], s_256[1]);
    s_256[2] = _mm256_add_epi32(s_256[2], s_256[3]);
    return _mm256_add_epi32(s_256[0], s_256[2]);
}

#endif
#if !REMOVE_UNUSED_CODE
static INLINE void update_best_sse2(const __m128i sad, const uint32_t mv,
                                    uint32_t *const p_best_sad_32x32,
                                    uint32_t *const p_best_mv32x32) {
    const __m128i  minpos  = _mm_minpos_epu16(sad);
    const uint32_t min_sad = (uint16_t)_mm_extract_epi16(minpos, 0);

    if (min_sad < *p_best_sad_32x32) {
        const int32_t best_mv = 4 * _mm_extract_epi16(minpos, 1);
        const int16_t x_mv    = _MVXT(mv) + (int16_t)best_mv;
        const int16_t y_mv    = _MVYT(mv);
        *p_best_sad_32x32     = min_sad;
        *p_best_mv32x32       = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
    }
}

/*******************************************
Calculate SAD for 32x32,64x64 from 16x16
and check if there is improvement, if yes keep
the best SAD+MV
*******************************************/
void get_eight_horizontal_search_point_results_32x32_64x64_pu_avx2_intrin(
    uint16_t *p_sad16x16, uint32_t *p_best_sad_32x32, uint32_t *p_best_sad_64x64,
    uint32_t *p_best_mv32x32, uint32_t *p_best_mv64x64, uint32_t mv) {
    __m256i s_256[4], sum[2];

    /*--------------------
    |  32x32_0  |  32x32_1
    ----------------------
    |  32x32_2  |  32x32_3
    ----------------------*/

    /*  data ordering in p_sad16x16 buffer

    Search    Search            Search
    Point 0   Point 1           Point 7
    ---------------------------------------
    16x16_0    |    x    |    x    | ...... |    x    |
    ---------------------------------------
    16x16_1    |    x    |    x    | ...... |    x    |

    16x16_n    |    x    |    x    | ...... |    x    |

    ---------------------------------------
    16x16_15   |    x    |    x    | ...... |    x    |
    ---------------------------------------
    */

    s_256[0] = add_4_sad_avx2(p_sad16x16 + 0 * 8);
    s_256[1] = add_4_sad_avx2(p_sad16x16 + 4 * 8);
    s_256[2] = add_4_sad_avx2(p_sad16x16 + 8 * 8);
    s_256[3] = add_4_sad_avx2(p_sad16x16 + 12 * 8);

    sum[0]              = _mm256_add_epi32(s_256[0], s_256[1]);
    sum[1]              = _mm256_add_epi32(s_256[2], s_256[3]);
    sum[0]              = _mm256_add_epi32(sum[0], sum[1]);
    const __m128i sad_0 = _mm256_castsi256_si128(sum[0]);
    const __m128i sad_1 = _mm256_extracti128_si256(sum[0], 1);

    const __m128i  ss      = _mm_packus_epi32(sad_0, sad_1);
    const __m128i  minpos  = _mm_minpos_epu16(ss);
    const uint32_t min_sad = (uint16_t)_mm_extract_epi16(minpos, 0);

    if (min_sad != 0xFFFF) {
        if (min_sad < *p_best_sad_64x64) {
            const uint32_t best_mv_64x64 = 4 * _mm_extract_epi16(minpos, 1);
            const int16_t  x_mv        = _MVXT(mv) + (int16_t)best_mv_64x64;
            const int16_t  y_mv        = _MVYT(mv);
            *p_best_sad_64x64          = min_sad;
            *p_best_mv64x64            = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        s_256[0] = _mm256_packus_epi32(s_256[0], s_256[1]);
        s_256[2] = _mm256_packus_epi32(s_256[2], s_256[3]);
        s_256[0] = _mm256_permute4x64_epi64(s_256[0], 0xD8);
        s_256[2] = _mm256_permute4x64_epi64(s_256[2], 0xD8);

        const __m128i s0_128 = _mm256_castsi256_si128(s_256[0]);
        const __m128i s1_128 = _mm256_extracti128_si256(s_256[0], 1);
        const __m128i s2_128 = _mm256_castsi256_si128(s_256[2]);
        const __m128i s3_128 = _mm256_extracti128_si256(s_256[2], 1);

        update_best_sse2(s0_128, mv, &p_best_sad_32x32[0], &p_best_mv32x32[0]);
        update_best_sse2(s1_128, mv, &p_best_sad_32x32[1], &p_best_mv32x32[1]);
        update_best_sse2(s2_128, mv, &p_best_sad_32x32[2], &p_best_mv32x32[2]);
        update_best_sse2(s3_128, mv, &p_best_sad_32x32[3], &p_best_mv32x32[3]);
    } else {
        __m128i        s0, s1, s2, s3, s4, s5;
        __m256i        ss0, ss1, ss2, ss3, ss4, ss5, ss6, ss7;
        const __m128i  min_val0    = _mm_min_epu32(sad_0, sad_1);
        const __m128i  min_val0_hi = _mm_srli_si128(min_val0, 8);
        const __m128i  min_val1    = _mm_min_epu32(min_val0, min_val0_hi);
        const __m128i  min_val1_hi = _mm_srli_si128(min_val1, 4);
        const __m128i  min_val2    = _mm_min_epu32(min_val1, min_val1_hi);
        const uint32_t min_val     = _mm_cvtsi128_si32(min_val2);

        if (min_val < *p_best_sad_64x64) {
            uint32_t sad_val[8];
            _mm256_storeu_si256((__m256i *)sad_val, sum[0]);

            for (int i = 0; i < 8; ++i) {
                if (min_val == sad_val[i]) {
                    const int16_t x_mv = _MVXT(mv) + (int16_t)i * 4;
                    const int16_t y_mv = _MVYT(mv);
                    *p_best_sad_64x64  = min_val;
                    *p_best_mv64x64    = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
                    break;
                }
            }
        }

        // XY
        // X: 32x32 block [0..3]
        // Y: Search position [0..7]
        ss0 = s_256[0]; // 07 06 05 04  03 02 01 00
        ss1 = s_256[1]; // 17 16 15 14  13 12 11 10
        ss2 = s_256[2]; // 27 26 25 24  23 22 21 20
        ss3 = s_256[3]; // 37 36 35 34  33 32 31 30
        ss4 = _mm256_unpacklo_epi32(ss0, ss1); // 15 05 14 04  11 01 10 00
        ss5 = _mm256_unpacklo_epi32(ss2, ss3); // 35 25 34 24  31 21 30 20
        ss6 = _mm256_unpackhi_epi32(ss0, ss1); // 17 07 16 06  13 03 12 02
        ss7 = _mm256_unpackhi_epi32(ss2, ss3); // 37 27 36 26  33 23 32 22
        ss0 = _mm256_unpacklo_epi64(ss4, ss5); // 34 24 14 04  30 20 10 00
        ss1 = _mm256_unpackhi_epi64(ss4, ss5); // 35 25 15 05  31 21 11 01
        ss2 = _mm256_unpacklo_epi64(ss6, ss7); // 36 26 16 06  32 22 12 02
        ss3 = _mm256_unpackhi_epi64(ss6, ss7); // 37 27 17 07  33 23 13 03

        //   ss4   |  ss5-2  |                ss6        |
        // a0 : a1 | a2 : a3 | min(a0, a1) : min(a2, a3) |    | (ss4 & !ss6) | ((ss5-2) & ss6) | ((ss4 & !ss6) | ((ss5-2) & ss6)) |
        // > (-1)  | >  (-3) |         >     (-1)        | a3 |       0      |       -3        |              -3                  |
        // > (-1)  | >  (-3) |         <=     (0)        | a1 |      -1      |        0        |              -1                  |
        // > (-1)  | <= (-2) |         >     (-1)        | a2 |       0      |       -2        |              -2                  |
        // > (-1)  | <= (-2) |         <=     (0)        | a1 |      -1      |        0        |              -1                  |
        // <= (0)  | >  (-3) |         >     (-1)        | a3 |       0      |       -3        |              -3                  |
        // <= (0)  | >  (-3) |         <=     (0)        | a0 |       0      |        0        |               0                  |
        // <= (0)  | <= (-2) |         >     (-1)        | a2 |       0      |       -2        |              -2                  |
        // <= (0)  | <= (-2) |         <=     (0)        | a0 |       0      |        0        |               0                  |

        // *** 8 search points per position ***

        // ss0: Search Pos 0,4 for blocks 0,1,2,3
        // ss1: Search Pos 1,5 for blocks 0,1,2,3
        // ss2: Search Pos 2,6 for blocks 0,1,2,3
        // ss3: Search Pos 3,7 for blocks 0,1,2,3

        ss4 = _mm256_cmpgt_epi32(ss0, ss1);
        ss0 = _mm256_min_epi32(ss0, ss1);
        ss5 = _mm256_cmpgt_epi32(ss2, ss3);
        ss2 = _mm256_min_epi32(ss2, ss3);
        ss5 = _mm256_sub_epi32(ss5, _mm256_set1_epi32(2)); // ss5-2

        // *** 4 search points per position ***
        ss6 = _mm256_cmpgt_epi32(ss0, ss2);
        ss0 = _mm256_min_epi32(ss0, ss2);
        ss5 = _mm256_and_si256(ss5, ss6); // (ss5-2) & ss6
        ss4 = _mm256_andnot_si256(ss6, ss4); // ss4 & !ss6
        ss4 = _mm256_or_si256(ss4, ss5); // (ss4 & !ss6) | ((ss5-2) & ss6)

        // *** 2 search points per position ***

        // ss0 = 8 SADs, two search points for each 32x32
        // ss4 = 8 MVs, two search points for each 32x32
        //
        // XY
        // X: 32x32 block [0..3]
        // Y: search position [0..1]
        // Format: 00 10 20 30  01 11 21 31

        // Each 128 bits contains 4 32x32 32bit block results
        // SAD
        s0 = _mm256_castsi256_si128(ss0);
        s1 = _mm256_extracti128_si256(ss0, 1);
        // MV
        s2 = _mm256_castsi256_si128(ss4);
        s3 = _mm256_extracti128_si256(ss4, 1);

        // Choose the best MV out of the two, use s4 to hold results of min
        s4 = _mm_cmpgt_epi32(s0, s1);
        s0 = _mm_min_epi32(s0, s1);

        // Extract MV's based on the blocks to s2
        s3 = _mm_sub_epi32(s3, _mm_set1_epi32(4)); // s3-4
        // Remove the MV's are not used from s2
        s2 = _mm_andnot_si128(s4, s2);
        // Remove the MV's that are not used from s3 (inverse from s2 above operation)
        s3 = _mm_and_si128(s4, s3);
        // Combine the remaining candidates into s2
        s2 = _mm_or_si128(s2, s3);
        // Convert MV's into encoders format
        s2 = _mm_slli_epi32(s2, 2); // mv info (negative)

        // ***SAD***
        // s0: current SAD candidates for each 32x32
        // s1: best SAD's for 32x32

        // Load best SAD's
        s1 = _mm_loadu_si128((__m128i *)p_best_sad_32x32);

        // Determine which candidates are better than the current best SAD's.
        // s4 is used to determine the MV's of the new best SAD's
        s4 = _mm_cmpgt_epi32(s1, s0);
        // Combine old and new min SAD's
        s0 = _mm_min_epu32(s0, s1);
        // Store new best SAD's back to memory
        _mm_storeu_si128((__m128i *)p_best_sad_32x32, s0);

        // ***Motion Vectors***
        // Load best MV's
        // s3: candidate MV's
        // s4: results of comparing SAD's
        // s5: previous best MV's

        // Load previous best MV's
        s5 = _mm_loadu_si128((__m128i *)p_best_mv32x32);
        // Remove the MV's that are being replaced
        s5 = _mm_andnot_si128(s4, s5);

        //Set mvx to s0 and mvy to s1
        s0 = _mm_set1_epi32(mv & 0xFFFF);
        s1 = _mm_set1_epi32(mv & 0xFFFF0000);

        // Add candidate MV's to base MV
        s3 = _mm_sub_epi32(s0, s2);
        // Limit to int16_t
        s3 = _mm_and_si128(s3, _mm_set1_epi32(0xFFFF));
        //mvy | mvx
        s3 = _mm_or_si128(s3, s1);

        // Remove non-candidate's
        s3 = _mm_and_si128(s3, s4);
        // Combine remaining candidates with remaining best MVs
        s3 = _mm_or_si128(s3, s5);
        // Store back to memory
        _mm_storeu_si128((__m128i *)p_best_mv32x32, s3);
    }
}
#endif
#if !REMOVE_UNUSED_CODE_PH2
/*******************************************************************************
* Requirement: width   = 4, 8, 16, 24, 32, 48 or 64
* Requirement: block_height <= 64
* Requirement: block_height % 2 = 0 when width = 4 or 8
*******************************************************************************/
void sad_loop_kernel_avx2_hme_l0_intrin(
    uint8_t * src, // input parameter, source samples Ptr
    uint32_t  src_stride, // input parameter, source stride
    uint8_t * ref, // input parameter, reference samples Ptr
    uint32_t  ref_stride, // input parameter, reference stride
    uint32_t  block_height, // input parameter, block height (M)
    uint32_t  block_width, // input parameter, block width (N)
    uint64_t *best_sad, int16_t *x_search_center, int16_t *y_search_center,
    uint32_t src_stride_raw, // input parameter, source stride (no line skipping)
    int16_t search_area_width, int16_t search_area_height) {
    int16_t        x_best = *x_search_center, y_best = *y_search_center;
    uint32_t       low_sum   = 0xffffff;
    uint32_t       tem_sum_1 = 0;
    int16_t        i, j;
    uint32_t       k;
    const uint8_t *p_ref, *p_src;
    __m128i        s0, s1, s2, s3, s4, s5, s6, s7 = _mm_set1_epi32(-1);
    __m256i        ss0, ss1, ss2, ss3, ss4, ss5, ss6, ss7, ss8, ss9, ss10, ss11;

    switch (block_width) {
    case 4:

        if (!(block_height % 4)) {
            uint32_t src_stride_t = 3 * src_stride;
            uint32_t ref_stride_t = 3 * ref_stride;
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss5 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k += 4) {
                        ss0 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_ref)),
                            _mm_loadu_si128((__m128i *)(p_ref + 2 * ref_stride)),
                            0x1);
                        ss1 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(
                                _mm_loadu_si128((__m128i *)(p_ref + ref_stride))),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride_t)),
                            0x1);
                        ss2 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_unpacklo_epi64(
                                _mm_cvtsi32_si128(*(uint32_t *)p_src),
                                _mm_cvtsi32_si128(*(uint32_t *)(p_src + src_stride)))),
                            _mm_unpacklo_epi64(
                                _mm_cvtsi32_si128(*(uint32_t *)(p_src + 2 * src_stride)),
                                _mm_cvtsi32_si128(*(uint32_t *)(p_src + src_stride_t))),
                            0x1);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        p_src += src_stride << 2;
                        p_ref += ref_stride << 2;
                    }
                    ss3       = _mm256_adds_epu16(ss3, ss5);
                    s3        = _mm_adds_epu16(_mm256_castsi256_si128(ss3),
                                        _mm256_extracti128_si256(ss3, 1));
                    s3        = _mm_minpos_epu16(s3);
                    tem_sum_1 = _mm_extract_epi16(s3, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s3, 1));
                        y_best  = i;
                    }
                }
                ref += src_stride_raw;
            }
        } else {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    s3    = _mm_setzero_si128();
                    for (k = 0; k < block_height; k += 2) {
                        s0 = _mm_loadu_si128((__m128i *)p_ref);
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + ref_stride));
                        s2 = _mm_cvtsi32_si128(*(uint32_t *)p_src);
                        s5 = _mm_cvtsi32_si128(*(uint32_t *)(p_src + src_stride));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s1, s5, 0));
                        p_src += src_stride << 1;
                        p_ref += ref_stride << 1;
                    }
                    s3        = _mm_minpos_epu16(s3);
                    tem_sum_1 = _mm_extract_epi16(s3, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s3, 1));
                        y_best  = i;
                    }
                }
                ref += src_stride_raw;
            }
        }

        break;

    case 8:
        if (!(block_height % 4)) {
            uint32_t src_stride_t = 3 * src_stride;
            uint32_t ref_stride_t = 3 * ref_stride;
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k += 4) {
                        ss0 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_ref)),
                            _mm_loadu_si128((__m128i *)(p_ref + 2 * ref_stride)),
                            0x1);
                        ss1 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(
                                _mm_loadu_si128((__m128i *)(p_ref + ref_stride))),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride_t)),
                            0x1);
                        ss2 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_unpacklo_epi64(
                                _mm_loadl_epi64((__m128i *)p_src),
                                _mm_loadl_epi64((__m128i *)(p_src + src_stride)))),
                            _mm_unpacklo_epi64(_mm_loadl_epi64((__m128i *)(p_src + 2 * src_stride)),
                                               _mm_loadl_epi64((__m128i *)(p_src + src_stride_t))),
                            0x1);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride << 2;
                        p_ref += ref_stride << 2;
                    }
                    ss3 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    s3        = _mm_adds_epu16(_mm256_castsi256_si128(ss3),
                                        _mm256_extracti128_si256(ss3, 1));
                    s3        = _mm_minpos_epu16(s3);
                    tem_sum_1 = _mm_extract_epi16(s3, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s3, 1));
                        y_best  = i;
                    }
                }
                ref += src_stride_raw;
            }
        } else {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    s3 = s4 = _mm_setzero_si128();
                    for (k = 0; k < block_height; k += 2) {
                        s0 = _mm_loadu_si128((__m128i *)p_ref);
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + ref_stride));
                        s2 = _mm_loadl_epi64((__m128i *)p_src);
                        s5 = _mm_loadl_epi64((__m128i *)(p_src + src_stride));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s1, s5, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s1, s5, 5));
                        p_src += src_stride << 1;
                        p_ref += ref_stride << 1;
                    }
                    s3        = _mm_adds_epu16(s3, s4);
                    s3        = _mm_minpos_epu16(s3);
                    tem_sum_1 = _mm_extract_epi16(s3, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s3, 1));
                        y_best  = i;
                    }
                }
                ref += src_stride_raw;
            }
        }
        break;

    case 16:
        if (block_height <= 16) {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 16; j += 16) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    ss7 = ss9 = ss10 = ss11 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k += 2) {
                        ss0 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_ref)),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride)),
                            0x1);
                        ss1 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)(p_ref + 8))),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride + 8)),
                            0x1);
                        ss2 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_src)),
                            _mm_loadu_si128((__m128i *)(p_src + src_stride)),
                            0x1);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111

                        ss0 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)(p_ref + 16))),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride + 16)),
                            0x1);
                        ss7  = _mm256_adds_epu16(ss7, _mm256_mpsadbw_epu8(ss1, ss2, 0));
                        ss11 = _mm256_adds_epu16(ss11, _mm256_mpsadbw_epu8(ss1, ss2, 45));
                        ss9  = _mm256_adds_epu16(ss9, _mm256_mpsadbw_epu8(ss0, ss2, 18));
                        ss10 = _mm256_adds_epu16(ss10, _mm256_mpsadbw_epu8(ss0, ss2, 63));

                        p_src += 2 * src_stride;
                        p_ref += 2 * ref_stride;
                    }
                    ss3 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    s3        = _mm_adds_epu16(_mm256_castsi256_si128(ss3),
                                        _mm256_extracti128_si256(ss3, 1));
                    s3        = _mm_minpos_epu16(s3);
                    tem_sum_1 = _mm_extract_epi16(s3, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s3, 1));
                        y_best  = i;
                    }

                    ss7       = _mm256_adds_epu16(_mm256_adds_epu16(ss7, ss11),
                                            _mm256_adds_epu16(ss9, ss10));
                    s7        = _mm_adds_epu16(_mm256_castsi256_si128(ss7),
                                        _mm256_extracti128_si256(ss7, 1));
                    s7        = _mm_minpos_epu16(s7);
                    tem_sum_1 = _mm_extract_epi16(s7, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = (int16_t)(j + 8 + _mm_extract_epi16(s7, 1));
                        y_best  = i;
                    }
                }
                ref += src_stride_raw;
            }
        } else if (block_height <= 32) {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k += 2) {
                        ss0 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_ref)),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride)),
                            0x1);
                        ss1 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)(p_ref + 8))),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride + 8)),
                            0x1);
                        ss2 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_src)),
                            _mm_loadu_si128((__m128i *)(p_src + src_stride)),
                            0x1);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += 2 * src_stride;
                        p_ref += 2 * ref_stride;
                    }
                    ss3 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    s3        = _mm256_castsi256_si128(ss3);
                    s5        = _mm256_extracti128_si256(ss3, 1);
                    s4        = _mm_minpos_epu16(s3);
                    s6        = _mm_minpos_epu16(s5);
                    s4        = _mm_unpacklo_epi16(s4, s4);
                    s4        = _mm_unpacklo_epi32(s4, s4);
                    s4        = _mm_unpacklo_epi64(s4, s4);
                    s6        = _mm_unpacklo_epi16(s6, s6);
                    s6        = _mm_unpacklo_epi32(s6, s6);
                    s6        = _mm_unpacklo_epi64(s6, s6);
                    s3        = _mm_sub_epi16(s3, s4);
                    s5        = _mm_adds_epu16(s5, s3);
                    s5        = _mm_sub_epi16(s5, s6);
                    s5        = _mm_minpos_epu16(s5);
                    tem_sum_1 = _mm_extract_epi16(s5, 0);
                    tem_sum_1 += _mm_extract_epi16(s4, 0);
                    tem_sum_1 += _mm_extract_epi16(s6, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s5, 1));
                        y_best  = i;
                    }
                }
                ref += src_stride_raw;
            }
        } else {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k += 2) {
                        ss0 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_ref)),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride)),
                            0x1);
                        ss1 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)(p_ref + 8))),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride + 8)),
                            0x1);
                        ss2 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_src)),
                            _mm_loadu_si128((__m128i *)(p_src + src_stride)),
                            0x1);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += 2 * src_stride;
                        p_ref += 2 * ref_stride;
                    }
                    ss3       = _mm256_adds_epu16(ss3, ss4);
                    ss5       = _mm256_adds_epu16(ss5, ss6);
                    ss0       = _mm256_adds_epu16(ss3, ss5);
                    s0        = _mm_adds_epu16(_mm256_castsi256_si128(ss0),
                                        _mm256_extracti128_si256(ss0, 1));
                    s0        = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            ss4 = _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256());
                            ss3 = _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256());
                            ss6 = _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256());
                            ss5 = _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256());
                            ss4 = _mm256_add_epi32(ss4, ss6);
                            ss3 = _mm256_add_epi32(ss3, ss5);
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss4),
                                               _mm256_extracti128_si256(ss4, 1));
                            s3  = _mm_add_epi32(_mm256_castsi256_si128(ss3),
                                               _mm256_extracti128_si256(ss3, 1));
                            UPDATE_BEST(s0, 0, 0);
                            UPDATE_BEST(s0, 1, 0);
                            UPDATE_BEST(s0, 2, 0);
                            UPDATE_BEST(s0, 3, 0);
                            UPDATE_BEST(s3, 0, 4);
                            UPDATE_BEST(s3, 1, 4);
                            UPDATE_BEST(s3, 2, 4);
                            UPDATE_BEST(s3, 3, 4);
                        }
                    }
                }
                ref += src_stride_raw;
            }
        }
        break;

    case 24:
        if (block_height <= 16) {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    ss3       = _mm256_adds_epu16(ss3, ss4);
                    ss5       = _mm256_adds_epu16(ss5, ss6);
                    s3        = _mm_adds_epu16(_mm256_castsi256_si128(ss3),
                                        _mm256_extracti128_si256(ss3, 1));
                    s5        = _mm256_castsi256_si128(ss5);
                    s4        = _mm_minpos_epu16(s3);
                    s6        = _mm_minpos_epu16(s5);
                    s4        = _mm_unpacklo_epi16(s4, s4);
                    s4        = _mm_unpacklo_epi32(s4, s4);
                    s4        = _mm_unpacklo_epi64(s4, s4);
                    s6        = _mm_unpacklo_epi16(s6, s6);
                    s6        = _mm_unpacklo_epi32(s6, s6);
                    s6        = _mm_unpacklo_epi64(s6, s6);
                    s3        = _mm_sub_epi16(s3, s4);
                    s5        = _mm_adds_epu16(s5, s3);
                    s5        = _mm_sub_epi16(s5, s6);
                    s5        = _mm_minpos_epu16(s5);
                    tem_sum_1 = _mm_extract_epi16(s5, 0);
                    tem_sum_1 += _mm_extract_epi16(s4, 0);
                    tem_sum_1 += _mm_extract_epi16(s6, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s5, 1));
                        y_best  = i;
                    }
                }
                ref += src_stride_raw;
            }
        } else {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    ss3       = _mm256_adds_epu16(ss3, ss4);
                    ss5       = _mm256_adds_epu16(ss5, ss6);
                    s3        = _mm256_castsi256_si128(ss3);
                    s4        = _mm256_extracti128_si256(ss3, 1);
                    s5        = _mm256_castsi256_si128(ss5);
                    s0        = _mm_adds_epu16(_mm_adds_epu16(s3, s4), s5);
                    s0        = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            s0 = _mm_unpacklo_epi16(s3, _mm_setzero_si128());
                            s3 = _mm_unpackhi_epi16(s3, _mm_setzero_si128());
                            s1 = _mm_unpacklo_epi16(s4, _mm_setzero_si128());
                            s4 = _mm_unpackhi_epi16(s4, _mm_setzero_si128());
                            s2 = _mm_unpacklo_epi16(s5, _mm_setzero_si128());
                            s5 = _mm_unpackhi_epi16(s5, _mm_setzero_si128());
                            s0 = _mm_add_epi32(_mm_add_epi32(s0, s1), s2);
                            s3 = _mm_add_epi32(_mm_add_epi32(s3, s4), s5);
                            UPDATE_BEST(s0, 0, 0);
                            UPDATE_BEST(s0, 1, 0);
                            UPDATE_BEST(s0, 2, 0);
                            UPDATE_BEST(s0, 3, 0);
                            UPDATE_BEST(s3, 0, 4);
                            UPDATE_BEST(s3, 1, 4);
                            UPDATE_BEST(s3, 2, 4);
                            UPDATE_BEST(s3, 3, 4);
                        }
                    }
                }
                ref += src_stride_raw;
            }
        }
        break;

    case 32:
        if (block_height <= 16) {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    ss3 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    s3        = _mm256_castsi256_si128(ss3);
                    s5        = _mm256_extracti128_si256(ss3, 1);
                    s4        = _mm_minpos_epu16(s3);
                    s6        = _mm_minpos_epu16(s5);
                    s4        = _mm_unpacklo_epi16(s4, s4);
                    s4        = _mm_unpacklo_epi32(s4, s4);
                    s4        = _mm_unpacklo_epi64(s4, s4);
                    s6        = _mm_unpacklo_epi16(s6, s6);
                    s6        = _mm_unpacklo_epi32(s6, s6);
                    s6        = _mm_unpacklo_epi64(s6, s6);
                    s3        = _mm_sub_epi16(s3, s4);
                    s5        = _mm_adds_epu16(s5, s3);
                    s5        = _mm_sub_epi16(s5, s6);
                    s5        = _mm_minpos_epu16(s5);
                    tem_sum_1 = _mm_extract_epi16(s5, 0);
                    tem_sum_1 += _mm_extract_epi16(s4, 0);
                    tem_sum_1 += _mm_extract_epi16(s6, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s5, 1));
                        y_best  = i;
                    }
                }
                ref += src_stride_raw;
            }
        } else if (block_height <= 32) {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    ss3       = _mm256_adds_epu16(ss3, ss4);
                    ss5       = _mm256_adds_epu16(ss5, ss6);
                    ss6       = _mm256_adds_epu16(ss3, ss5);
                    s3        = _mm256_castsi256_si128(ss6);
                    s4        = _mm256_extracti128_si256(ss6, 1);
                    s0        = _mm_adds_epu16(s3, s4);
                    s0        = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            ss4 = _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256());
                            ss3 = _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256());
                            ss6 = _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256());
                            ss5 = _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256());
                            ss4 = _mm256_add_epi32(ss4, ss6);
                            ss3 = _mm256_add_epi32(ss3, ss5);
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss4),
                                               _mm256_extracti128_si256(ss4, 1));
                            s3  = _mm_add_epi32(_mm256_castsi256_si128(ss3),
                                               _mm256_extracti128_si256(ss3, 1));
                            UPDATE_BEST(s0, 0, 0);
                            UPDATE_BEST(s0, 1, 0);
                            UPDATE_BEST(s0, 2, 0);
                            UPDATE_BEST(s0, 3, 0);
                            UPDATE_BEST(s3, 0, 4);
                            UPDATE_BEST(s3, 1, 4);
                            UPDATE_BEST(s3, 2, 4);
                            UPDATE_BEST(s3, 3, 4);
                        }
                    }
                }
                ref += src_stride_raw;
            }
        } else {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    ss7 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    s3        = _mm256_castsi256_si128(ss7);
                    s4        = _mm256_extracti128_si256(ss7, 1);
                    s0        = _mm_adds_epu16(s3, s4);
                    s0        = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            ss0 = _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256());
                            ss3 = _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256());
                            ss1 = _mm256_unpacklo_epi16(ss4, _mm256_setzero_si256());
                            ss4 = _mm256_unpackhi_epi16(ss4, _mm256_setzero_si256());
                            ss2 = _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256());
                            ss5 = _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256());
                            ss7 = _mm256_unpacklo_epi16(ss6, _mm256_setzero_si256());
                            ss6 = _mm256_unpackhi_epi16(ss6, _mm256_setzero_si256());
                            ss0 = _mm256_add_epi32(_mm256_add_epi32(ss0, ss1),
                                                   _mm256_add_epi32(ss2, ss7));
                            ss3 = _mm256_add_epi32(_mm256_add_epi32(ss3, ss4),
                                                   _mm256_add_epi32(ss5, ss6));
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss0),
                                               _mm256_extracti128_si256(ss0, 1));
                            s3  = _mm_add_epi32(_mm256_castsi256_si128(ss3),
                                               _mm256_extracti128_si256(ss3, 1));
                            UPDATE_BEST(s0, 0, 0);
                            UPDATE_BEST(s0, 1, 0);
                            UPDATE_BEST(s0, 2, 0);
                            UPDATE_BEST(s0, 3, 0);
                            UPDATE_BEST(s3, 0, 4);
                            UPDATE_BEST(s3, 1, 4);
                            UPDATE_BEST(s3, 2, 4);
                            UPDATE_BEST(s3, 3, 4);
                        }
                    }
                }
                ref += src_stride_raw;
            }
        }
        break;

    case 48:
        if (block_height <= 32) {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        s0  = _mm_loadu_si128((__m128i *)(p_ref + 32));
                        s1  = _mm_loadu_si128((__m128i *)(p_ref + 40));
                        s2  = _mm_loadu_si128((__m128i *)(p_src + 32));
                        s3  = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4  = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5  = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6  = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    s3        = _mm_adds_epu16(s3, s4);
                    s5        = _mm_adds_epu16(s5, s6);
                    s0        = _mm_adds_epu16(s3, s5);
                    ss3       = _mm256_adds_epu16(ss3, ss4);
                    ss5       = _mm256_adds_epu16(ss5, ss6);
                    ss6       = _mm256_adds_epu16(ss3, ss5);
                    s0        = _mm_adds_epu16(s0,
                                        _mm_adds_epu16(_mm256_castsi256_si128(ss6),
                                                       _mm256_extracti128_si256(ss6, 1)));
                    s0        = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            ss4 = _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256());
                            ss3 = _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256());
                            ss6 = _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256());
                            ss5 = _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256());
                            ss4 = _mm256_add_epi32(ss4, ss6);
                            ss3 = _mm256_add_epi32(ss3, ss5);
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss4),
                                               _mm256_extracti128_si256(ss4, 1));
                            s1  = _mm_add_epi32(_mm256_castsi256_si128(ss3),
                                               _mm256_extracti128_si256(ss3, 1));
                            s0  = _mm_add_epi32(s0, _mm_unpacklo_epi16(s3, _mm_setzero_si128()));
                            s0  = _mm_add_epi32(s0, _mm_unpacklo_epi16(s5, _mm_setzero_si128()));
                            s1  = _mm_add_epi32(s1, _mm_unpackhi_epi16(s3, _mm_setzero_si128()));
                            s1  = _mm_add_epi32(s1, _mm_unpackhi_epi16(s5, _mm_setzero_si128()));
                            UPDATE_BEST(s0, 0, 0);
                            UPDATE_BEST(s0, 1, 0);
                            UPDATE_BEST(s0, 2, 0);
                            UPDATE_BEST(s0, 3, 0);
                            UPDATE_BEST(s1, 0, 4);
                            UPDATE_BEST(s1, 1, 4);
                            UPDATE_BEST(s1, 2, 4);
                            UPDATE_BEST(s1, 3, 4);
                        }
                    }
                }
                ref += src_stride_raw;
            }
        } else {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        s0  = _mm_loadu_si128((__m128i *)(p_ref + 32));
                        s1  = _mm_loadu_si128((__m128i *)(p_ref + 40));
                        s2  = _mm_loadu_si128((__m128i *)(p_src + 32));
                        s3  = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4  = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5  = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6  = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    s0        = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    ss7       = _mm256_adds_epu16(ss3, ss4);
                    ss8       = _mm256_adds_epu16(ss5, ss6);
                    ss7       = _mm256_adds_epu16(ss7, ss8);
                    s0        = _mm_adds_epu16(s0,
                                        _mm_adds_epu16(_mm256_castsi256_si128(ss7),
                                                       _mm256_extracti128_si256(ss7, 1)));
                    s0        = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            ss0 = _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256());
                            ss3 = _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256());
                            ss1 = _mm256_unpacklo_epi16(ss4, _mm256_setzero_si256());
                            ss4 = _mm256_unpackhi_epi16(ss4, _mm256_setzero_si256());
                            ss2 = _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256());
                            ss5 = _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256());
                            ss7 = _mm256_unpacklo_epi16(ss6, _mm256_setzero_si256());
                            ss6 = _mm256_unpackhi_epi16(ss6, _mm256_setzero_si256());
                            ss0 = _mm256_add_epi32(_mm256_add_epi32(ss0, ss1),
                                                   _mm256_add_epi32(ss2, ss7));
                            ss3 = _mm256_add_epi32(_mm256_add_epi32(ss3, ss4),
                                                   _mm256_add_epi32(ss5, ss6));
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss0),
                                               _mm256_extracti128_si256(ss0, 1));
                            s1  = _mm_add_epi32(_mm256_castsi256_si128(ss3),
                                               _mm256_extracti128_si256(ss3, 1));
                            s0  = _mm_add_epi32(s0, _mm_unpacklo_epi16(s3, _mm_setzero_si128()));
                            s0  = _mm_add_epi32(s0, _mm_unpacklo_epi16(s4, _mm_setzero_si128()));
                            s0  = _mm_add_epi32(s0, _mm_unpacklo_epi16(s5, _mm_setzero_si128()));
                            s0  = _mm_add_epi32(s0, _mm_unpacklo_epi16(s6, _mm_setzero_si128()));
                            s1  = _mm_add_epi32(s1, _mm_unpackhi_epi16(s3, _mm_setzero_si128()));
                            s1  = _mm_add_epi32(s1, _mm_unpackhi_epi16(s4, _mm_setzero_si128()));
                            s1  = _mm_add_epi32(s1, _mm_unpackhi_epi16(s5, _mm_setzero_si128()));
                            s1  = _mm_add_epi32(s1, _mm_unpackhi_epi16(s6, _mm_setzero_si128()));
                            UPDATE_BEST(s0, 0, 0);
                            UPDATE_BEST(s0, 1, 0);
                            UPDATE_BEST(s0, 2, 0);
                            UPDATE_BEST(s0, 3, 0);
                            UPDATE_BEST(s1, 0, 4);
                            UPDATE_BEST(s1, 1, 4);
                            UPDATE_BEST(s1, 2, 4);
                            UPDATE_BEST(s1, 3, 4);
                        }
                    }
                }
                ref += src_stride_raw;
            }
        }
        break;

    case 64:
        if (block_height <= 32) {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        ss0 = _mm256_loadu_si256((__m256i *)(p_ref + 32));
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 40));
                        ss2 = _mm256_loadu_si256((__m256i *)(p_src + 32));
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    ss7 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    s3        = _mm256_castsi256_si128(ss7);
                    s4        = _mm256_extracti128_si256(ss7, 1);
                    s0        = _mm_adds_epu16(s3, s4);
                    s0        = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            ss0 = _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256());
                            ss3 = _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256());
                            ss1 = _mm256_unpacklo_epi16(ss4, _mm256_setzero_si256());
                            ss4 = _mm256_unpackhi_epi16(ss4, _mm256_setzero_si256());
                            ss2 = _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256());
                            ss5 = _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256());
                            ss7 = _mm256_unpacklo_epi16(ss6, _mm256_setzero_si256());
                            ss6 = _mm256_unpackhi_epi16(ss6, _mm256_setzero_si256());
                            ss0 = _mm256_add_epi32(_mm256_add_epi32(ss0, ss1),
                                                   _mm256_add_epi32(ss2, ss7));
                            ss3 = _mm256_add_epi32(_mm256_add_epi32(ss3, ss4),
                                                   _mm256_add_epi32(ss5, ss6));
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss0),
                                               _mm256_extracti128_si256(ss0, 1));
                            s3  = _mm_add_epi32(_mm256_castsi256_si128(ss3),
                                               _mm256_extracti128_si256(ss3, 1));
                            UPDATE_BEST(s0, 0, 0);
                            UPDATE_BEST(s0, 1, 0);
                            UPDATE_BEST(s0, 2, 0);
                            UPDATE_BEST(s0, 3, 0);
                            UPDATE_BEST(s3, 0, 4);
                            UPDATE_BEST(s3, 1, 4);
                            UPDATE_BEST(s3, 2, 4);
                            UPDATE_BEST(s3, 3, 4);
                        }
                    }
                }
                ref += src_stride_raw;
            }
        } else {
            __m256i ss9, ss10;
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = ss7 = ss8 = ss9 = ss10 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        ss0 = _mm256_loadu_si256((__m256i *)(p_ref + 32));
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 40));
                        ss2 = _mm256_loadu_si256((__m256i *)(p_src + 32));
                        ss7 = _mm256_adds_epu16(ss7, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss8 = _mm256_adds_epu16(ss8, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss9 = _mm256_adds_epu16(ss9, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss10 =
                            _mm256_adds_epu16(ss10, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    ss0 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    ss0       = _mm256_adds_epu16(ss0,
                                            _mm256_adds_epu16(_mm256_adds_epu16(ss7, ss8),
                                                              _mm256_adds_epu16(ss9, ss10)));
                    s0        = _mm_adds_epu16(_mm256_castsi256_si128(ss0),
                                        _mm256_extracti128_si256(ss0, 1));
                    s0        = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            ss0 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss4, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss6, _mm256_setzero_si256())));
                            ss1 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss4, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss6, _mm256_setzero_si256())));
                            ss2 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss7, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss8, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss9, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss10, _mm256_setzero_si256())));
                            ss3 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss7, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss8, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss9, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss10, _mm256_setzero_si256())));
                            ss0 = _mm256_add_epi32(ss0, ss2);
                            ss1 = _mm256_add_epi32(ss1, ss3);
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss0),
                                               _mm256_extracti128_si256(ss0, 1));
                            s3  = _mm_add_epi32(_mm256_castsi256_si128(ss1),
                                               _mm256_extracti128_si256(ss1, 1));
                            UPDATE_BEST(s0, 0, 0);
                            UPDATE_BEST(s0, 1, 0);
                            UPDATE_BEST(s0, 2, 0);
                            UPDATE_BEST(s0, 3, 0);
                            UPDATE_BEST(s3, 0, 4);
                            UPDATE_BEST(s3, 1, 4);
                            UPDATE_BEST(s3, 2, 4);
                            UPDATE_BEST(s3, 3, 4);
                        }
                    }
                }
                ref += src_stride_raw;
            }
        }
        break;

    default: assert(0); break;
    }

    *best_sad        = low_sum;
    *x_search_center = x_best;
    *y_search_center = y_best;
}
#endif
#if RESTRUCTURE_SAD
#define UPDATE_BEST_PME(s, k, offset)                                               \
    tem_sum_1 = _mm_extract_epi32(s, k);                                            \
    if (tem_sum_1 < low_sum) {                                                      \
        low_sum = tem_sum_1;                                                        \
        x_best  = mvx + (search_position_start_x + (j + offset + k)) * search_step; \
        y_best  = mvy + (search_position_start_y + i) * search_step;                \
    }

/*******************************************************************************
* performs sad search for given block and search area
* return best sad with its best mvx/mvy
* Requirement: width   = 4, 8, 16, 24, 32, 48, 64, 128
* Requirement: block_height <= 128 when width = 64 or 128, <= 64 otherwise
* Requirement: block_height % 2 = 0 when width = 4 or 8
*******************************************************************************/
void pme_sad_loop_kernel_avx2(uint8_t * src, // input parameter, source samples Ptr
                              uint32_t  src_stride, // input parameter, source stride
                              uint8_t * ref, // input parameter, reference samples Ptr
                              uint32_t  ref_stride, // input parameter, reference stride
                              uint32_t  block_height, // input parameter, block height (M)
                              uint32_t  block_width, // input parameter, block width (N)
                              uint32_t *best_sad, int16_t *best_mvx, int16_t *best_mvy,
                              int16_t search_position_start_x, int16_t search_position_start_y,
                              int16_t search_area_width, int16_t search_area_height,
                              int16_t search_step, int16_t mvx, int16_t mvy) {
    int16_t        x_best = *best_mvx, y_best = *best_mvy;
    uint32_t       low_sum  = *best_sad;
    uint32_t       tem_sum_1 = 0;
    int16_t        i, j;
    uint32_t       k;
    const uint8_t *p_ref, *p_src;
    __m128i        s0, s1, s2, s3, s4, s5, s6, s7 = _mm_set1_epi32(-1);
    __m256i        ss0, ss1, ss2, ss3, ss4, ss5, ss6, ss7, ss8, ss9, ss10, ss11;
    __m256i        ss12, ss13, ss14, ss15, ss16, ss17, ss18;

    switch (block_width) {
    case 4:
        if (block_height == 4) {
            uint32_t src_stride_t = 3 * src_stride;
            uint32_t ref_stride_t = 3 * ref_stride;
            ss2                 = _mm256_insertf128_si256(
                _mm256_castsi128_si256(
                    _mm_unpacklo_epi64(_mm_cvtsi32_si128(*(uint32_t *)src),
                                       _mm_cvtsi32_si128(*(uint32_t *)(src + src_stride)))),
                _mm_unpacklo_epi64(_mm_cvtsi32_si128(*(uint32_t *)(src + 2 * src_stride)),
                                   _mm_cvtsi32_si128(*(uint32_t *)(src + src_stride_t))),
                0x1);
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_ref = ref + j;
                    ss3 = ss5 = _mm256_setzero_si256();
                    ss0       = _mm256_insertf128_si256(
                        _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_ref)),
                        _mm_loadu_si128((__m128i *)(p_ref + 2 * ref_stride)),
                        0x1);
                    ss1 = _mm256_insertf128_si256(
                        _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)(p_ref + ref_stride))),
                        _mm_loadu_si128((__m128i *)(p_ref + ref_stride_t)),
                        0x1);
                    ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                    ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                    p_ref += ref_stride << 2;
                    ss3     = _mm256_adds_epu16(ss3, ss5);
                    s3      = _mm_adds_epu16(_mm256_castsi256_si128(ss3),
                                        _mm256_extracti128_si256(ss3, 1));
                    s3      = _mm_minpos_epu16(s3);
                    tem_sum_1 = _mm_extract_epi16(s3, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = mvx + (search_position_start_x +
                                       (int16_t)(j + _mm_extract_epi16(s3, 1))) *
                                          search_step;
                        y_best = mvy + (search_position_start_y + i) * search_step;
                    }
                }
                ref += ref_stride;
            }
        } else if (block_height == 8) {
            uint32_t src_stride_t = 3 * src_stride;
            uint32_t ref_stride_t = 3 * ref_stride;
            p_src                = src;
            ss2                 = _mm256_insertf128_si256(
                _mm256_castsi128_si256(
                    _mm_unpacklo_epi64(_mm_cvtsi32_si128(*(uint32_t *)p_src),
                                       _mm_cvtsi32_si128(*(uint32_t *)(p_src + src_stride)))),
                _mm_unpacklo_epi64(_mm_cvtsi32_si128(*(uint32_t *)(p_src + 2 * src_stride)),
                                   _mm_cvtsi32_si128(*(uint32_t *)(p_src + src_stride_t))),
                0x1);
            p_src += src_stride << 2;
            ss6 = _mm256_insertf128_si256(
                _mm256_castsi128_si256(
                    _mm_unpacklo_epi64(_mm_cvtsi32_si128(*(uint32_t *)p_src),
                                       _mm_cvtsi32_si128(*(uint32_t *)(p_src + src_stride)))),
                _mm_unpacklo_epi64(_mm_cvtsi32_si128(*(uint32_t *)(p_src + 2 * src_stride)),
                                   _mm_cvtsi32_si128(*(uint32_t *)(p_src + src_stride_t))),
                0x1);
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_ref = ref + j;
                    ss3 = ss5 = _mm256_setzero_si256();
                    ss0       = _mm256_insertf128_si256(
                        _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_ref)),
                        _mm_loadu_si128((__m128i *)(p_ref + 2 * ref_stride)),
                        0x1);
                    ss1 = _mm256_insertf128_si256(
                        _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)(p_ref + ref_stride))),
                        _mm_loadu_si128((__m128i *)(p_ref + ref_stride_t)),
                        0x1);
                    ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                    ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                    p_ref += ref_stride << 2;
                    ss0 = _mm256_insertf128_si256(
                        _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_ref)),
                        _mm_loadu_si128((__m128i *)(p_ref + 2 * ref_stride)),
                        0x1);
                    ss1 = _mm256_insertf128_si256(
                        _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)(p_ref + ref_stride))),
                        _mm_loadu_si128((__m128i *)(p_ref + ref_stride_t)),
                        0x1);
                    ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss6, 0));
                    ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss6, 18)); // 010 010

                    ss3     = _mm256_adds_epu16(ss3, ss5);
                    s3      = _mm_adds_epu16(_mm256_castsi256_si128(ss3),
                                        _mm256_extracti128_si256(ss3, 1));
                    s3      = _mm_minpos_epu16(s3);
                    tem_sum_1 = _mm_extract_epi16(s3, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = mvx + (search_position_start_x +
                                       (int16_t)(j + _mm_extract_epi16(s3, 1))) *
                                          search_step;
                        y_best = mvy + (search_position_start_y + i) * search_step;
                    }
                }
                ref += ref_stride;
            }
        } else if (block_height == 16) {
            uint32_t src_stride_t = 3 * src_stride;
            uint32_t ref_stride_t = 3 * ref_stride;
            p_src                = src;
            ss2                 = _mm256_insertf128_si256(
                _mm256_castsi128_si256(
                    _mm_unpacklo_epi64(_mm_cvtsi32_si128(*(uint32_t *)p_src),
                                       _mm_cvtsi32_si128(*(uint32_t *)(p_src + src_stride)))),
                _mm_unpacklo_epi64(_mm_cvtsi32_si128(*(uint32_t *)(p_src + 2 * src_stride)),
                                   _mm_cvtsi32_si128(*(uint32_t *)(p_src + src_stride_t))),
                0x1);
            p_src += src_stride << 2;
            ss6 = _mm256_insertf128_si256(
                _mm256_castsi128_si256(
                    _mm_unpacklo_epi64(_mm_cvtsi32_si128(*(uint32_t *)p_src),
                                       _mm_cvtsi32_si128(*(uint32_t *)(p_src + src_stride)))),
                _mm_unpacklo_epi64(_mm_cvtsi32_si128(*(uint32_t *)(p_src + 2 * src_stride)),
                                   _mm_cvtsi32_si128(*(uint32_t *)(p_src + src_stride_t))),
                0x1);
            p_src += src_stride << 2;
            ss7 = _mm256_insertf128_si256(
                _mm256_castsi128_si256(
                    _mm_unpacklo_epi64(_mm_cvtsi32_si128(*(uint32_t *)p_src),
                                       _mm_cvtsi32_si128(*(uint32_t *)(p_src + src_stride)))),
                _mm_unpacklo_epi64(_mm_cvtsi32_si128(*(uint32_t *)(p_src + 2 * src_stride)),
                                   _mm_cvtsi32_si128(*(uint32_t *)(p_src + src_stride_t))),
                0x1);
            p_src += src_stride << 2;
            ss8 = _mm256_insertf128_si256(
                _mm256_castsi128_si256(
                    _mm_unpacklo_epi64(_mm_cvtsi32_si128(*(uint32_t *)p_src),
                                       _mm_cvtsi32_si128(*(uint32_t *)(p_src + src_stride)))),
                _mm_unpacklo_epi64(_mm_cvtsi32_si128(*(uint32_t *)(p_src + 2 * src_stride)),
                                   _mm_cvtsi32_si128(*(uint32_t *)(p_src + src_stride_t))),
                0x1);
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_ref = ref + j;
                    ss3 = ss5 = _mm256_setzero_si256();
                    ss0       = _mm256_insertf128_si256(
                        _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_ref)),
                        _mm_loadu_si128((__m128i *)(p_ref + 2 * ref_stride)),
                        0x1);
                    ss1 = _mm256_insertf128_si256(
                        _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)(p_ref + ref_stride))),
                        _mm_loadu_si128((__m128i *)(p_ref + ref_stride_t)),
                        0x1);
                    ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                    ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                    p_ref += ref_stride << 2;
                    ss0 = _mm256_insertf128_si256(
                        _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_ref)),
                        _mm_loadu_si128((__m128i *)(p_ref + 2 * ref_stride)),
                        0x1);
                    ss1 = _mm256_insertf128_si256(
                        _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)(p_ref + ref_stride))),
                        _mm_loadu_si128((__m128i *)(p_ref + ref_stride_t)),
                        0x1);
                    ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss6, 0));
                    ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss6, 18)); // 010 010
                    p_ref += ref_stride << 2;
                    ss0 = _mm256_insertf128_si256(
                        _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_ref)),
                        _mm_loadu_si128((__m128i *)(p_ref + 2 * ref_stride)),
                        0x1);
                    ss1 = _mm256_insertf128_si256(
                        _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)(p_ref + ref_stride))),
                        _mm_loadu_si128((__m128i *)(p_ref + ref_stride_t)),
                        0x1);
                    ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss7, 0));
                    ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss7, 18)); // 010 010
                    p_ref += ref_stride << 2;
                    ss0 = _mm256_insertf128_si256(
                        _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_ref)),
                        _mm_loadu_si128((__m128i *)(p_ref + 2 * ref_stride)),
                        0x1);
                    ss1 = _mm256_insertf128_si256(
                        _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)(p_ref + ref_stride))),
                        _mm_loadu_si128((__m128i *)(p_ref + ref_stride_t)),
                        0x1);
                    ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss8, 0));
                    ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss8, 18)); // 010 010

                    ss3     = _mm256_adds_epu16(ss3, ss5);
                    s3      = _mm_adds_epu16(_mm256_castsi256_si128(ss3),
                                        _mm256_extracti128_si256(ss3, 1));
                    s3      = _mm_minpos_epu16(s3);
                    tem_sum_1 = _mm_extract_epi16(s3, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = mvx + (search_position_start_x +
                                       (int16_t)(j + _mm_extract_epi16(s3, 1))) *
                                          search_step;
                        y_best = mvy + (search_position_start_y + i) * search_step;
                    }
                }
                ref += ref_stride;
            }
        }
        else if (!(block_height % 4)) {
            uint32_t src_stride_t = 3 * src_stride;
            uint32_t ref_stride_t = 3 * ref_stride;
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss5 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k += 4) {
                        ss0 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_ref)),
                            _mm_loadu_si128((__m128i *)(p_ref + 2 * ref_stride)),
                            0x1);
                        ss1 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)(p_ref + ref_stride))),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride_t)),
                            0x1);
                        ss2 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_unpacklo_epi64(
                                _mm_cvtsi32_si128(*(uint32_t *)p_src),
                                _mm_cvtsi32_si128(*(uint32_t *)(p_src + src_stride)))),
                            _mm_unpacklo_epi64(
                                _mm_cvtsi32_si128(*(uint32_t *)(p_src + 2 * src_stride)),
                                _mm_cvtsi32_si128(*(uint32_t *)(p_src + src_stride_t))),
                            0x1);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        p_src += src_stride << 2;
                        p_ref += ref_stride << 2;
                    }
                    ss3     = _mm256_adds_epu16(ss3, ss5);
                    s3      = _mm_adds_epu16(_mm256_castsi256_si128(ss3),
                                        _mm256_extracti128_si256(ss3, 1));
                    s3      = _mm_minpos_epu16(s3);
                    tem_sum_1 = _mm_extract_epi16(s3, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = mvx + (search_position_start_x +
                                       (int16_t)(j + _mm_extract_epi16(s3, 1))) *
                                          search_step;
                        y_best = mvy + (search_position_start_y + i) * search_step;
                    }
                }
                ref += ref_stride;
            }
        } else {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    s3   = _mm_setzero_si128();
                    for (k = 0; k < block_height; k += 2) {
                        s0 = _mm_loadu_si128((__m128i *)p_ref);
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + ref_stride));
                        s2 = _mm_cvtsi32_si128(*(uint32_t *)p_src);
                        s5 = _mm_cvtsi32_si128(*(uint32_t *)(p_src + src_stride));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s1, s5, 0));
                        p_src += src_stride << 1;
                        p_ref += ref_stride << 1;
                    }
                    s3      = _mm_minpos_epu16(s3);
                    tem_sum_1 = _mm_extract_epi16(s3, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = mvx + (search_position_start_x +
                                       (int16_t)(j + _mm_extract_epi16(s3, 1))) *
                                          search_step;
                        y_best = mvy + (search_position_start_y + i) * search_step;
                    }
                }
                ref += ref_stride;
            }
        }
        break;

    case 8:
        if (block_height == 4) {
            uint32_t src_stride_t = 3 * src_stride;
            uint32_t ref_stride_t = 3 * ref_stride;
            p_src                = src;
            ss2                 = _mm256_insertf128_si256(
                _mm256_castsi128_si256(
                    _mm_unpacklo_epi64(_mm_loadl_epi64((__m128i *)p_src),
                                       _mm_loadl_epi64((__m128i *)(p_src + src_stride)))),
                _mm_unpacklo_epi64(_mm_loadl_epi64((__m128i *)(p_src + 2 * src_stride)),
                                   _mm_loadl_epi64((__m128i *)(p_src + src_stride_t))),
                0x1);
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    ss0                   = _mm256_insertf128_si256(
                        _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_ref)),
                        _mm_loadu_si128((__m128i *)(p_ref + 2 * ref_stride)),
                        0x1);
                    ss1 = _mm256_insertf128_si256(
                        _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)(p_ref + ref_stride))),
                        _mm_loadu_si128((__m128i *)(p_ref + ref_stride_t)),
                        0x1);
                    ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                    ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                    ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                    ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                    ss3 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    s3      = _mm_adds_epu16(_mm256_castsi256_si128(ss3),
                                        _mm256_extracti128_si256(ss3, 1));
                    s3      = _mm_minpos_epu16(s3);
                    tem_sum_1 = _mm_extract_epi16(s3, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = mvx + (search_position_start_x +
                                       (int16_t)(j + _mm_extract_epi16(s3, 1))) *
                                          search_step;
                        y_best = mvy + (search_position_start_y + i) * search_step;
                    }
                }
                ref += ref_stride;
            }
        } else if (block_height == 8) {
            uint32_t src_stride_t = 3 * src_stride;
            uint32_t ref_stride_t = 3 * ref_stride;
            p_src               = src;
            ss2                 = _mm256_insertf128_si256(
                _mm256_castsi128_si256(
                    _mm_unpacklo_epi64(_mm_loadl_epi64((__m128i *)p_src),
                                       _mm_loadl_epi64((__m128i *)(p_src + src_stride)))),
                _mm_unpacklo_epi64(_mm_loadl_epi64((__m128i *)(p_src + 2 * src_stride)),
                                   _mm_loadl_epi64((__m128i *)(p_src + src_stride_t))),
                0x1);
            p_src += src_stride << 2;
            ss7 = _mm256_insertf128_si256(
                _mm256_castsi128_si256(
                    _mm_unpacklo_epi64(_mm_loadl_epi64((__m128i *)p_src),
                                       _mm_loadl_epi64((__m128i *)(p_src + src_stride)))),
                _mm_unpacklo_epi64(_mm_loadl_epi64((__m128i *)(p_src + 2 * src_stride)),
                                   _mm_loadl_epi64((__m128i *)(p_src + src_stride_t))),
                0x1);
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    ss0                   = _mm256_insertf128_si256(
                        _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_ref)),
                        _mm_loadu_si128((__m128i *)(p_ref + 2 * ref_stride)),
                        0x1);
                    ss1 = _mm256_insertf128_si256(
                        _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)(p_ref + ref_stride))),
                        _mm_loadu_si128((__m128i *)(p_ref + ref_stride_t)),
                        0x1);
                    ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                    ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                    ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                    ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                    p_ref += ref_stride << 2;
                    ss0 = _mm256_insertf128_si256(
                        _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_ref)),
                        _mm_loadu_si128((__m128i *)(p_ref + 2 * ref_stride)),
                        0x1);
                    ss1 = _mm256_insertf128_si256(
                        _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)(p_ref + ref_stride))),
                        _mm_loadu_si128((__m128i *)(p_ref + ref_stride_t)),
                        0x1);
                    ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss7, 0));
                    ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss7, 45)); // 101 101
                    ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss7, 18)); // 010 010
                    ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss7, 63)); // 111 111
                    ss3 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    s3      = _mm_adds_epu16(_mm256_castsi256_si128(ss3),
                                        _mm256_extracti128_si256(ss3, 1));
                    s3      = _mm_minpos_epu16(s3);
                    tem_sum_1 = _mm_extract_epi16(s3, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = mvx + (search_position_start_x +
                                       (int16_t)(j + _mm_extract_epi16(s3, 1))) *
                                          search_step;
                        y_best = mvy + (search_position_start_y + i) * search_step;
                    }
                }
                ref += ref_stride;
            }
        }
        else if (!(block_height % 4)) {
            uint32_t src_stride_t = 3 * src_stride;
            uint32_t ref_stride_t = 3 * ref_stride;
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k += 4) {
                        ss0 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_ref)),
                            _mm_loadu_si128((__m128i *)(p_ref + 2 * ref_stride)),
                            0x1);
                        ss1 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)(p_ref + ref_stride))),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride_t)),
                            0x1);
                        ss2 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_unpacklo_epi64(
                                _mm_loadl_epi64((__m128i *)p_src),
                                _mm_loadl_epi64((__m128i *)(p_src + src_stride)))),
                            _mm_unpacklo_epi64(_mm_loadl_epi64((__m128i *)(p_src + 2 * src_stride)),
                                               _mm_loadl_epi64((__m128i *)(p_src + src_stride_t))),
                            0x1);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride << 2;
                        p_ref += ref_stride << 2;
                    }
                    ss3 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    s3      = _mm_adds_epu16(_mm256_castsi256_si128(ss3),
                                        _mm256_extracti128_si256(ss3, 1));
                    s3      = _mm_minpos_epu16(s3);
                    tem_sum_1 = _mm_extract_epi16(s3, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = mvx + (search_position_start_x +
                                       (int16_t)(j + _mm_extract_epi16(s3, 1))) *
                                          search_step;
                        y_best = mvy + (search_position_start_y + i) * search_step;
                    }
                }
                ref += ref_stride;
            }
        } else {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    s3 = s4 = _mm_setzero_si128();
                    for (k = 0; k < block_height; k += 2) {
                        s0 = _mm_loadu_si128((__m128i *)p_ref);
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + ref_stride));
                        s2 = _mm_loadl_epi64((__m128i *)p_src);
                        s5 = _mm_loadl_epi64((__m128i *)(p_src + src_stride));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s1, s5, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s1, s5, 5));
                        p_src += src_stride << 1;
                        p_ref += ref_stride << 1;
                    }
                    s3      = _mm_adds_epu16(s3, s4);
                    s3      = _mm_minpos_epu16(s3);
                    tem_sum_1 = _mm_extract_epi16(s3, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = mvx + (search_position_start_x +
                                       (int16_t)(j + _mm_extract_epi16(s3, 1))) *
                                          search_step;
                        y_best = mvy + (search_position_start_y + i) * search_step;
                    }
                }
                ref += ref_stride;
            }
        }
        break;

    case 16:
        if (block_height <= 16 && !(search_area_width & 15)) {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 16; j += 16) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    ss7 = ss9 = ss10 = ss11 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k += 2) {
                        ss0 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_ref)),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride)),
                            0x1);
                        ss1 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)(p_ref + 8))),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride + 8)),
                            0x1);
                        ss2 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_src)),
                            _mm_loadu_si128((__m128i *)(p_src + src_stride)),
                            0x1);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111

                        ss0 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)(p_ref + 16))),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride + 16)),
                            0x1);
                        ss7  = _mm256_adds_epu16(ss7, _mm256_mpsadbw_epu8(ss1, ss2, 0));
                        ss11 = _mm256_adds_epu16(ss11, _mm256_mpsadbw_epu8(ss1, ss2, 45));
                        ss9  = _mm256_adds_epu16(ss9, _mm256_mpsadbw_epu8(ss0, ss2, 18));
                        ss10 = _mm256_adds_epu16(ss10, _mm256_mpsadbw_epu8(ss0, ss2, 63));

                        p_src += 2 * src_stride;
                        p_ref += 2 * ref_stride;
                    }
                    ss3 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    s3      = _mm_adds_epu16(_mm256_castsi256_si128(ss3),
                                        _mm256_extracti128_si256(ss3, 1));
                    s3      = _mm_minpos_epu16(s3);
                    tem_sum_1 = _mm_extract_epi16(s3, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = mvx + (search_position_start_x +
                                       (int16_t)(j + _mm_extract_epi16(s3, 1))) *
                                          search_step;
                        y_best = mvy + (search_position_start_y + i) * search_step;
                        ;
                    }

                    ss7     = _mm256_adds_epu16(_mm256_adds_epu16(ss7, ss11),
                                            _mm256_adds_epu16(ss9, ss10));
                    s7      = _mm_adds_epu16(_mm256_castsi256_si128(ss7),
                                        _mm256_extracti128_si256(ss7, 1));
                    s7      = _mm_minpos_epu16(s7);
                    tem_sum_1 = _mm_extract_epi16(s7, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = mvx + (search_position_start_x +
                                       (int16_t)(j + 8 + _mm_extract_epi16(s7, 1))) *
                                          search_step;
                        y_best = mvy + (search_position_start_y + i) * search_step;
                    }
                }
                ref += ref_stride;
            }
        } else if (block_height <= 32) {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k += 2) {
                        ss0 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_ref)),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride)),
                            0x1);
                        ss1 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)(p_ref + 8))),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride + 8)),
                            0x1);
                        ss2 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_src)),
                            _mm_loadu_si128((__m128i *)(p_src + src_stride)),
                            0x1);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += 2 * src_stride;
                        p_ref += 2 * ref_stride;
                    }
                    ss3 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    s3      = _mm256_castsi256_si128(ss3);
                    s5      = _mm256_extracti128_si256(ss3, 1);
                    s4      = _mm_minpos_epu16(s3);
                    s6      = _mm_minpos_epu16(s5);
                    s4      = _mm_unpacklo_epi16(s4, s4);
                    s4      = _mm_unpacklo_epi32(s4, s4);
                    s4      = _mm_unpacklo_epi64(s4, s4);
                    s6      = _mm_unpacklo_epi16(s6, s6);
                    s6      = _mm_unpacklo_epi32(s6, s6);
                    s6      = _mm_unpacklo_epi64(s6, s6);
                    s3      = _mm_sub_epi16(s3, s4);
                    s5      = _mm_adds_epu16(s5, s3);
                    s5      = _mm_sub_epi16(s5, s6);
                    s5      = _mm_minpos_epu16(s5);
                    tem_sum_1 = _mm_extract_epi16(s5, 0);
                    tem_sum_1 += _mm_extract_epi16(s4, 0);
                    tem_sum_1 += _mm_extract_epi16(s6, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = mvx + (search_position_start_x +
                                       (int16_t)(j + _mm_extract_epi16(s5, 1))) *
                                          search_step;
                        y_best = mvy + (search_position_start_y + i) * search_step;
                    }
                }
                ref += ref_stride;
            }
        } else {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k += 2) {
                        ss0 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_ref)),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride)),
                            0x1);
                        ss1 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)(p_ref + 8))),
                            _mm_loadu_si128((__m128i *)(p_ref + ref_stride + 8)),
                            0x1);
                        ss2 = _mm256_insertf128_si256(
                            _mm256_castsi128_si256(_mm_loadu_si128((__m128i *)p_src)),
                            _mm_loadu_si128((__m128i *)(p_src + src_stride)),
                            0x1);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += 2 * src_stride;
                        p_ref += 2 * ref_stride;
                    }
                    ss3     = _mm256_adds_epu16(ss3, ss4);
                    ss5     = _mm256_adds_epu16(ss5, ss6);
                    ss0     = _mm256_adds_epu16(ss3, ss5);
                    s0      = _mm_adds_epu16(_mm256_castsi256_si128(ss0),
                                        _mm256_extracti128_si256(ss0, 1));
                    s0      = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = mvx + (search_position_start_x +
                                           (int16_t)(j + _mm_extract_epi16(s0, 1))) *
                                              search_step;
                            y_best = mvy + (search_position_start_y + i) * search_step;
                        } else {
                            ss4 = _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256());
                            ss3 = _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256());
                            ss6 = _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256());
                            ss5 = _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256());
                            ss4 = _mm256_add_epi32(ss4, ss6);
                            ss3 = _mm256_add_epi32(ss3, ss5);
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss4),
                                               _mm256_extracti128_si256(ss4, 1));
                            s3  = _mm_add_epi32(_mm256_castsi256_si128(ss3),
                                               _mm256_extracti128_si256(ss3, 1));
                            UPDATE_BEST_PME(s0, 0, 0);
                            UPDATE_BEST_PME(s0, 1, 0);
                            UPDATE_BEST_PME(s0, 2, 0);
                            UPDATE_BEST_PME(s0, 3, 0);
                            UPDATE_BEST_PME(s3, 0, 4);
                            UPDATE_BEST_PME(s3, 1, 4);
                            UPDATE_BEST_PME(s3, 2, 4);
                            UPDATE_BEST_PME(s3, 3, 4);
                        }
                    }
                }
                ref += ref_stride;
            }
        }
        break;

    case 24:
        if (block_height <= 16) {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    ss3     = _mm256_adds_epu16(ss3, ss4);
                    ss5     = _mm256_adds_epu16(ss5, ss6);
                    s3      = _mm_adds_epu16(_mm256_castsi256_si128(ss3),
                                        _mm256_extracti128_si256(ss3, 1));
                    s5      = _mm256_castsi256_si128(ss5);
                    s4      = _mm_minpos_epu16(s3);
                    s6      = _mm_minpos_epu16(s5);
                    s4      = _mm_unpacklo_epi16(s4, s4);
                    s4      = _mm_unpacklo_epi32(s4, s4);
                    s4      = _mm_unpacklo_epi64(s4, s4);
                    s6      = _mm_unpacklo_epi16(s6, s6);
                    s6      = _mm_unpacklo_epi32(s6, s6);
                    s6      = _mm_unpacklo_epi64(s6, s6);
                    s3      = _mm_sub_epi16(s3, s4);
                    s5      = _mm_adds_epu16(s5, s3);
                    s5      = _mm_sub_epi16(s5, s6);
                    s5      = _mm_minpos_epu16(s5);
                    tem_sum_1 = _mm_extract_epi16(s5, 0);
                    tem_sum_1 += _mm_extract_epi16(s4, 0);
                    tem_sum_1 += _mm_extract_epi16(s6, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = mvx + (search_position_start_x +
                                       (int16_t)(j + _mm_extract_epi16(s5, 1))) *
                                          search_step;
                        y_best = mvy + (search_position_start_y + i) * search_step;
                    }
                }
                ref += ref_stride;
            }
        } else {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    ss3     = _mm256_adds_epu16(ss3, ss4);
                    ss5     = _mm256_adds_epu16(ss5, ss6);
                    s3      = _mm256_castsi256_si128(ss3);
                    s4      = _mm256_extracti128_si256(ss3, 1);
                    s5      = _mm256_castsi256_si128(ss5);
                    s0      = _mm_adds_epu16(_mm_adds_epu16(s3, s4), s5);
                    s0      = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = mvx + (search_position_start_x +
                                           (int16_t)(j + _mm_extract_epi16(s0, 1))) *
                                              search_step;
                            y_best = mvy + (search_position_start_y + i) * search_step;
                        } else {
                            s0 = _mm_unpacklo_epi16(s3, _mm_setzero_si128());
                            s3 = _mm_unpackhi_epi16(s3, _mm_setzero_si128());
                            s1 = _mm_unpacklo_epi16(s4, _mm_setzero_si128());
                            s4 = _mm_unpackhi_epi16(s4, _mm_setzero_si128());
                            s2 = _mm_unpacklo_epi16(s5, _mm_setzero_si128());
                            s5 = _mm_unpackhi_epi16(s5, _mm_setzero_si128());
                            s0 = _mm_add_epi32(_mm_add_epi32(s0, s1), s2);
                            s3 = _mm_add_epi32(_mm_add_epi32(s3, s4), s5);
                            UPDATE_BEST_PME(s0, 0, 0);
                            UPDATE_BEST_PME(s0, 1, 0);
                            UPDATE_BEST_PME(s0, 2, 0);
                            UPDATE_BEST_PME(s0, 3, 0);
                            UPDATE_BEST_PME(s3, 0, 4);
                            UPDATE_BEST_PME(s3, 1, 4);
                            UPDATE_BEST_PME(s3, 2, 4);
                            UPDATE_BEST_PME(s3, 3, 4);
                        }
                    }
                }
                ref += ref_stride;
            }
        }
        break;

    case 32:
        if (block_height <= 16) {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    ss3 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    s3      = _mm256_castsi256_si128(ss3);
                    s5      = _mm256_extracti128_si256(ss3, 1);
                    s4      = _mm_minpos_epu16(s3);
                    s6      = _mm_minpos_epu16(s5);
                    s4      = _mm_unpacklo_epi16(s4, s4);
                    s4      = _mm_unpacklo_epi32(s4, s4);
                    s4      = _mm_unpacklo_epi64(s4, s4);
                    s6      = _mm_unpacklo_epi16(s6, s6);
                    s6      = _mm_unpacklo_epi32(s6, s6);
                    s6      = _mm_unpacklo_epi64(s6, s6);
                    s3      = _mm_sub_epi16(s3, s4);
                    s5      = _mm_adds_epu16(s5, s3);
                    s5      = _mm_sub_epi16(s5, s6);
                    s5      = _mm_minpos_epu16(s5);
                    tem_sum_1 = _mm_extract_epi16(s5, 0);
                    tem_sum_1 += _mm_extract_epi16(s4, 0);
                    tem_sum_1 += _mm_extract_epi16(s6, 0);
                    if (tem_sum_1 < low_sum) {
                        low_sum = tem_sum_1;
                        x_best  = mvx + (search_position_start_x +
                                       (int16_t)(j + _mm_extract_epi16(s5, 1))) *
                                          search_step;
                        y_best = mvy + (search_position_start_y + i) * search_step;
                    }
                }
                ref += ref_stride;
            }
        } else if (block_height <= 32) {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    ss3     = _mm256_adds_epu16(ss3, ss4);
                    ss5     = _mm256_adds_epu16(ss5, ss6);
                    ss6     = _mm256_adds_epu16(ss3, ss5);
                    s3      = _mm256_castsi256_si128(ss6);
                    s4      = _mm256_extracti128_si256(ss6, 1);
                    s0      = _mm_adds_epu16(s3, s4);
                    s0      = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = mvx + (search_position_start_x +
                                           (int16_t)(j + _mm_extract_epi16(s0, 1))) *
                                              search_step;
                            y_best = mvy + (search_position_start_y + i) * search_step;
                        } else {
                            ss4 = _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256());
                            ss3 = _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256());
                            ss6 = _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256());
                            ss5 = _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256());
                            ss4 = _mm256_add_epi32(ss4, ss6);
                            ss3 = _mm256_add_epi32(ss3, ss5);
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss4),
                                               _mm256_extracti128_si256(ss4, 1));
                            s3  = _mm_add_epi32(_mm256_castsi256_si128(ss3),
                                               _mm256_extracti128_si256(ss3, 1));
                            UPDATE_BEST_PME(s0, 0, 0);
                            UPDATE_BEST_PME(s0, 1, 0);
                            UPDATE_BEST_PME(s0, 2, 0);
                            UPDATE_BEST_PME(s0, 3, 0);
                            UPDATE_BEST_PME(s3, 0, 4);
                            UPDATE_BEST_PME(s3, 1, 4);
                            UPDATE_BEST_PME(s3, 2, 4);
                            UPDATE_BEST_PME(s3, 3, 4);
                        }
                    }
                }
                ref += ref_stride;
            }
        } else {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    ss7 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    s3      = _mm256_castsi256_si128(ss7);
                    s4      = _mm256_extracti128_si256(ss7, 1);
                    s0      = _mm_adds_epu16(s3, s4);
                    s0      = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = mvx + (search_position_start_x +
                                           (int16_t)(j + _mm_extract_epi16(s0, 1))) *
                                              search_step;
                            y_best = mvy + (search_position_start_y + i) * search_step;
                        } else {
                            ss0 = _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256());
                            ss3 = _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256());
                            ss1 = _mm256_unpacklo_epi16(ss4, _mm256_setzero_si256());
                            ss4 = _mm256_unpackhi_epi16(ss4, _mm256_setzero_si256());
                            ss2 = _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256());
                            ss5 = _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256());
                            ss7 = _mm256_unpacklo_epi16(ss6, _mm256_setzero_si256());
                            ss6 = _mm256_unpackhi_epi16(ss6, _mm256_setzero_si256());
                            ss0 = _mm256_add_epi32(_mm256_add_epi32(ss0, ss1),
                                                   _mm256_add_epi32(ss2, ss7));
                            ss3 = _mm256_add_epi32(_mm256_add_epi32(ss3, ss4),
                                                   _mm256_add_epi32(ss5, ss6));
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss0),
                                               _mm256_extracti128_si256(ss0, 1));
                            s3  = _mm_add_epi32(_mm256_castsi256_si128(ss3),
                                               _mm256_extracti128_si256(ss3, 1));
                            UPDATE_BEST_PME(s0, 0, 0);
                            UPDATE_BEST_PME(s0, 1, 0);
                            UPDATE_BEST_PME(s0, 2, 0);
                            UPDATE_BEST_PME(s0, 3, 0);
                            UPDATE_BEST_PME(s3, 0, 4);
                            UPDATE_BEST_PME(s3, 1, 4);
                            UPDATE_BEST_PME(s3, 2, 4);
                            UPDATE_BEST_PME(s3, 3, 4);
                        }
                    }
                }
                ref += ref_stride;
            }
        }
        break;

    case 48:
        if (block_height <= 32) {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        s0  = _mm_loadu_si128((__m128i *)(p_ref + 32));
                        s1  = _mm_loadu_si128((__m128i *)(p_ref + 40));
                        s2  = _mm_loadu_si128((__m128i *)(p_src + 32));
                        s3  = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4  = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5  = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6  = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    s3      = _mm_adds_epu16(s3, s4);
                    s5      = _mm_adds_epu16(s5, s6);
                    s0      = _mm_adds_epu16(s3, s5);
                    ss3     = _mm256_adds_epu16(ss3, ss4);
                    ss5     = _mm256_adds_epu16(ss5, ss6);
                    ss6     = _mm256_adds_epu16(ss3, ss5);
                    s0      = _mm_adds_epu16(s0,
                                        _mm_adds_epu16(_mm256_castsi256_si128(ss6),
                                                       _mm256_extracti128_si256(ss6, 1)));
                    s0      = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = mvx + (search_position_start_x +
                                           (int16_t)(j + _mm_extract_epi16(s0, 1))) *
                                              search_step;
                            y_best = mvy + (search_position_start_y + i) * search_step;
                        } else {
                            ss4 = _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256());
                            ss3 = _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256());
                            ss6 = _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256());
                            ss5 = _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256());
                            ss4 = _mm256_add_epi32(ss4, ss6);
                            ss3 = _mm256_add_epi32(ss3, ss5);
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss4),
                                               _mm256_extracti128_si256(ss4, 1));
                            s1  = _mm_add_epi32(_mm256_castsi256_si128(ss3),
                                               _mm256_extracti128_si256(ss3, 1));
                            s0  = _mm_add_epi32(s0, _mm_unpacklo_epi16(s3, _mm_setzero_si128()));
                            s0  = _mm_add_epi32(s0, _mm_unpacklo_epi16(s5, _mm_setzero_si128()));
                            s1  = _mm_add_epi32(s1, _mm_unpackhi_epi16(s3, _mm_setzero_si128()));
                            s1  = _mm_add_epi32(s1, _mm_unpackhi_epi16(s5, _mm_setzero_si128()));
                            UPDATE_BEST_PME(s0, 0, 0);
                            UPDATE_BEST_PME(s0, 1, 0);
                            UPDATE_BEST_PME(s0, 2, 0);
                            UPDATE_BEST_PME(s0, 3, 0);
                            UPDATE_BEST_PME(s1, 0, 4);
                            UPDATE_BEST_PME(s1, 1, 4);
                            UPDATE_BEST_PME(s1, 2, 4);
                            UPDATE_BEST_PME(s1, 3, 4);
                        }
                    }
                }
                ref += ref_stride;
            }
        } else {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        s0  = _mm_loadu_si128((__m128i *)(p_ref + 32));
                        s1  = _mm_loadu_si128((__m128i *)(p_ref + 40));
                        s2  = _mm_loadu_si128((__m128i *)(p_src + 32));
                        s3  = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4  = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5  = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6  = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    s0      = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    ss7     = _mm256_adds_epu16(ss3, ss4);
                    ss8     = _mm256_adds_epu16(ss5, ss6);
                    ss7     = _mm256_adds_epu16(ss7, ss8);
                    s0      = _mm_adds_epu16(s0,
                                        _mm_adds_epu16(_mm256_castsi256_si128(ss7),
                                                       _mm256_extracti128_si256(ss7, 1)));
                    s0      = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = mvx + (search_position_start_x +
                                           (int16_t)(j + _mm_extract_epi16(s0, 1))) *
                                              search_step;
                            y_best = mvy + (search_position_start_y + i) * search_step;
                        } else {
                            ss0 = _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256());
                            ss3 = _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256());
                            ss1 = _mm256_unpacklo_epi16(ss4, _mm256_setzero_si256());
                            ss4 = _mm256_unpackhi_epi16(ss4, _mm256_setzero_si256());
                            ss2 = _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256());
                            ss5 = _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256());
                            ss7 = _mm256_unpacklo_epi16(ss6, _mm256_setzero_si256());
                            ss6 = _mm256_unpackhi_epi16(ss6, _mm256_setzero_si256());
                            ss0 = _mm256_add_epi32(_mm256_add_epi32(ss0, ss1),
                                                   _mm256_add_epi32(ss2, ss7));
                            ss3 = _mm256_add_epi32(_mm256_add_epi32(ss3, ss4),
                                                   _mm256_add_epi32(ss5, ss6));
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss0),
                                               _mm256_extracti128_si256(ss0, 1));
                            s1  = _mm_add_epi32(_mm256_castsi256_si128(ss3),
                                               _mm256_extracti128_si256(ss3, 1));
                            s0  = _mm_add_epi32(s0, _mm_unpacklo_epi16(s3, _mm_setzero_si128()));
                            s0  = _mm_add_epi32(s0, _mm_unpacklo_epi16(s4, _mm_setzero_si128()));
                            s0  = _mm_add_epi32(s0, _mm_unpacklo_epi16(s5, _mm_setzero_si128()));
                            s0  = _mm_add_epi32(s0, _mm_unpacklo_epi16(s6, _mm_setzero_si128()));
                            s1  = _mm_add_epi32(s1, _mm_unpackhi_epi16(s3, _mm_setzero_si128()));
                            s1  = _mm_add_epi32(s1, _mm_unpackhi_epi16(s4, _mm_setzero_si128()));
                            s1  = _mm_add_epi32(s1, _mm_unpackhi_epi16(s5, _mm_setzero_si128()));
                            s1  = _mm_add_epi32(s1, _mm_unpackhi_epi16(s6, _mm_setzero_si128()));
                            UPDATE_BEST_PME(s0, 0, 0);
                            UPDATE_BEST_PME(s0, 1, 0);
                            UPDATE_BEST_PME(s0, 2, 0);
                            UPDATE_BEST_PME(s0, 3, 0);
                            UPDATE_BEST_PME(s1, 0, 4);
                            UPDATE_BEST_PME(s1, 1, 4);
                            UPDATE_BEST_PME(s1, 2, 4);
                            UPDATE_BEST_PME(s1, 3, 4);
                        }
                    }
                }
                ref += ref_stride;
            }
        }
        break;

    case 64:
        if (block_height <= 32) {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        ss0 = _mm256_loadu_si256((__m256i *)(p_ref + 32));
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 40));
                        ss2 = _mm256_loadu_si256((__m256i *)(p_src + 32));
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    ss7 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    s3      = _mm256_castsi256_si128(ss7);
                    s4      = _mm256_extracti128_si256(ss7, 1);
                    s0      = _mm_adds_epu16(s3, s4);
                    s0      = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = mvx + (search_position_start_x +
                                           (int16_t)(j + _mm_extract_epi16(s0, 1))) *
                                              search_step;
                            y_best = mvy + (search_position_start_y + i) * search_step;
                        } else {
                            ss0 = _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256());
                            ss3 = _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256());
                            ss1 = _mm256_unpacklo_epi16(ss4, _mm256_setzero_si256());
                            ss4 = _mm256_unpackhi_epi16(ss4, _mm256_setzero_si256());
                            ss2 = _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256());
                            ss5 = _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256());
                            ss7 = _mm256_unpacklo_epi16(ss6, _mm256_setzero_si256());
                            ss6 = _mm256_unpackhi_epi16(ss6, _mm256_setzero_si256());
                            ss0 = _mm256_add_epi32(_mm256_add_epi32(ss0, ss1),
                                                   _mm256_add_epi32(ss2, ss7));
                            ss3 = _mm256_add_epi32(_mm256_add_epi32(ss3, ss4),
                                                   _mm256_add_epi32(ss5, ss6));
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss0),
                                               _mm256_extracti128_si256(ss0, 1));
                            s3  = _mm_add_epi32(_mm256_castsi256_si128(ss3),
                                               _mm256_extracti128_si256(ss3, 1));
                            UPDATE_BEST_PME(s0, 0, 0);
                            UPDATE_BEST_PME(s0, 1, 0);
                            UPDATE_BEST_PME(s0, 2, 0);
                            UPDATE_BEST_PME(s0, 3, 0);
                            UPDATE_BEST_PME(s3, 0, 4);
                            UPDATE_BEST_PME(s3, 1, 4);
                            UPDATE_BEST_PME(s3, 2, 4);
                            UPDATE_BEST_PME(s3, 3, 4);
                        }
                    }
                }
                ref += ref_stride;
            }
        } else if (block_height <= 64) {
            __m256i ss9, ss10;
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = ss7 = ss8 = ss9 = ss10 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        ss0 = _mm256_loadu_si256((__m256i *)(p_ref + 32));
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 40));
                        ss2 = _mm256_loadu_si256((__m256i *)(p_src + 32));
                        ss7 = _mm256_adds_epu16(ss7, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss8 = _mm256_adds_epu16(ss8, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss9 = _mm256_adds_epu16(ss9, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss10 =
                            _mm256_adds_epu16(ss10, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    ss0 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    ss0     = _mm256_adds_epu16(ss0,
                                            _mm256_adds_epu16(_mm256_adds_epu16(ss7, ss8),
                                                              _mm256_adds_epu16(ss9, ss10)));
                    s0      = _mm_adds_epu16(_mm256_castsi256_si128(ss0),
                                        _mm256_extracti128_si256(ss0, 1));
                    s0      = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = mvx + (search_position_start_x +
                                           (int16_t)(j + _mm_extract_epi16(s0, 1))) *
                                              search_step;
                            y_best = mvy + (search_position_start_y + i) * search_step;
                        } else {
                            ss0 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss4, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss6, _mm256_setzero_si256())));
                            ss1 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss4, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss6, _mm256_setzero_si256())));
                            ss2 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss7, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss8, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss9, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss10, _mm256_setzero_si256())));
                            ss3 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss7, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss8, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss9, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss10, _mm256_setzero_si256())));
                            ss0 = _mm256_add_epi32(ss0, ss2);
                            ss1 = _mm256_add_epi32(ss1, ss3);
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss0),
                                               _mm256_extracti128_si256(ss0, 1));
                            s3  = _mm_add_epi32(_mm256_castsi256_si128(ss1),
                                               _mm256_extracti128_si256(ss1, 1));
                            UPDATE_BEST_PME(s0, 0, 0);
                            UPDATE_BEST_PME(s0, 1, 0);
                            UPDATE_BEST_PME(s0, 2, 0);
                            UPDATE_BEST_PME(s0, 3, 0);
                            UPDATE_BEST_PME(s3, 0, 4);
                            UPDATE_BEST_PME(s3, 1, 4);
                            UPDATE_BEST_PME(s3, 2, 4);
                            UPDATE_BEST_PME(s3, 3, 4);
                        }
                    }
                }
                ref += ref_stride;
            }
        } else {
            __m256i ss9, ss10;
            __m256i ssa1, ssa2, ssa3, ssa4, ssa5, ssa6, ssa7, ssa8;
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = ss7 = ss8 = ss9 = ss10 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        ss0 = _mm256_loadu_si256((__m256i *)(p_ref + 32));
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 40));
                        ss2 = _mm256_loadu_si256((__m256i *)(p_src + 32));
                        ss7 = _mm256_adds_epu16(ss7, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss8 = _mm256_adds_epu16(ss8, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss9 = _mm256_adds_epu16(ss9, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss10 =
                            _mm256_adds_epu16(ss10, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                        if (k == 63) {
                            ssa1 = ss3;
                            ssa2 = ss4;
                            ssa3 = ss5;
                            ssa4 = ss6;
                            ssa5 = ss7;
                            ssa6 = ss8;
                            ssa7 = ss9;
                            ssa8 = ss10;
                            ss3 = ss4 = ss5 = ss6 = ss7 = ss8 = ss9 = ss10 = _mm256_setzero_si256();
                        }
                    }
                    ss0 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    ss0     = _mm256_adds_epu16(ss0,
                                            _mm256_adds_epu16(_mm256_adds_epu16(ss7, ss8),
                                                              _mm256_adds_epu16(ss9, ss10)));
                    ss1     = _mm256_adds_epu16(_mm256_adds_epu16(ssa1, ssa2),
                                            _mm256_adds_epu16(ssa3, ssa4));
                    ss1     = _mm256_adds_epu16(ss1,
                                            _mm256_adds_epu16(_mm256_adds_epu16(ssa5, ssa6),
                                                              _mm256_adds_epu16(ssa7, ssa8)));
                    ss0     = _mm256_adds_epu16(ss0, ss1);
                    s0      = _mm_adds_epu16(_mm256_castsi256_si128(ss0),
                                        _mm256_extracti128_si256(ss0, 1));
                    s0      = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = mvx + (search_position_start_x +
                                           (int16_t)(j + _mm_extract_epi16(s0, 1))) *
                                              search_step;
                            y_best = mvy + (search_position_start_y + i) * search_step;
                        } else {
                            ss0 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss4, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss6, _mm256_setzero_si256())));
                            ss1 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss4, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss6, _mm256_setzero_si256())));
                            ss2 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss7, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss8, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss9, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss10, _mm256_setzero_si256())));
                            ss3 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss7, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss8, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss9, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss10, _mm256_setzero_si256())));
                            ss4 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ssa1, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ssa2, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ssa3, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ssa4, _mm256_setzero_si256())));
                            ss5 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ssa1, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ssa2, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ssa3, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ssa4, _mm256_setzero_si256())));
                            ss6 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ssa5, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ssa6, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ssa7, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ssa8, _mm256_setzero_si256())));
                            ss7 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ssa5, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ssa6, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ssa7, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ssa8, _mm256_setzero_si256())));
                            ss0 = _mm256_add_epi32(ss0, ss2);
                            ss1 = _mm256_add_epi32(ss1, ss3);
                            ss2 = _mm256_add_epi32(ss4, ss6);
                            ss3 = _mm256_add_epi32(ss5, ss7);
                            ss0 = _mm256_add_epi32(ss0, ss2);
                            ss1 = _mm256_add_epi32(ss1, ss3);
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss0),
                                               _mm256_extracti128_si256(ss0, 1));
                            s3  = _mm_add_epi32(_mm256_castsi256_si128(ss1),
                                               _mm256_extracti128_si256(ss1, 1));
                            UPDATE_BEST_PME(s0, 0, 0);
                            UPDATE_BEST_PME(s0, 1, 0);
                            UPDATE_BEST_PME(s0, 2, 0);
                            UPDATE_BEST_PME(s0, 3, 0);
                            UPDATE_BEST_PME(s3, 0, 4);
                            UPDATE_BEST_PME(s3, 1, 4);
                            UPDATE_BEST_PME(s3, 2, 4);
                            UPDATE_BEST_PME(s3, 3, 4);
                        }
                    }
                }
                ref += ref_stride;
            }
        }
        break;
    case 128:
        if (block_height <= 64) {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = ss7 = ss8 = ss9 = ss10 = ss11 = ss12 = ss13 = ss14 =
                        ss15 = ss16 = ss17 = ss18 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        ss0 = _mm256_loadu_si256((__m256i *)(p_ref + 32));
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 40));
                        ss2 = _mm256_loadu_si256((__m256i *)(p_src + 32));
                        ss7 = _mm256_adds_epu16(ss7, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss8 = _mm256_adds_epu16(ss8, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss9 = _mm256_adds_epu16(ss9, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss10 =
                            _mm256_adds_epu16(ss10, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        ss0  = _mm256_loadu_si256((__m256i *)(p_ref + 64));
                        ss1  = _mm256_loadu_si256((__m256i *)(p_ref + 72));
                        ss2  = _mm256_loadu_si256((__m256i *)(p_src + 64));
                        ss11 = _mm256_adds_epu16(ss11, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss12 =
                            _mm256_adds_epu16(ss12, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss13 =
                            _mm256_adds_epu16(ss13, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss14 =
                            _mm256_adds_epu16(ss14, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        ss0  = _mm256_loadu_si256((__m256i *)(p_ref + 96));
                        ss1  = _mm256_loadu_si256((__m256i *)(p_ref + 104));
                        ss2  = _mm256_loadu_si256((__m256i *)(p_src + 96));
                        ss15 = _mm256_adds_epu16(ss15, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss16 =
                            _mm256_adds_epu16(ss16, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss17 =
                            _mm256_adds_epu16(ss17, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss18 =
                            _mm256_adds_epu16(ss18, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    ss0 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    ss0     = _mm256_adds_epu16(ss0,
                                            _mm256_adds_epu16(_mm256_adds_epu16(ss7, ss8),
                                                              _mm256_adds_epu16(ss9, ss10)));
                    ss0     = _mm256_adds_epu16(ss0,
                                            _mm256_adds_epu16(_mm256_adds_epu16(ss11, ss12),
                                                              _mm256_adds_epu16(ss13, ss14)));
                    ss0     = _mm256_adds_epu16(ss0,
                                            _mm256_adds_epu16(_mm256_adds_epu16(ss15, ss16),
                                                              _mm256_adds_epu16(ss17, ss18)));
                    s0      = _mm_adds_epu16(_mm256_castsi256_si128(ss0),
                                        _mm256_extracti128_si256(ss0, 1));
                    s0      = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = mvx + (search_position_start_x +
                                           (int16_t)(j + _mm_extract_epi16(s0, 1))) *
                                              search_step;
                            y_best = mvy + (search_position_start_y + i) * search_step;
                        } else {
                            ss0 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss4, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss6, _mm256_setzero_si256())));
                            ss1 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss4, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss6, _mm256_setzero_si256())));
                            ss2 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss7, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss8, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss9, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss10, _mm256_setzero_si256())));
                            ss3 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss7, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss8, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss9, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss10, _mm256_setzero_si256())));
                            ss4 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss11, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss12, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss13, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss14, _mm256_setzero_si256())));
                            ss5 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss11, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss12, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss13, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss14, _mm256_setzero_si256())));
                            ss6 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss15, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss16, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss17, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss18, _mm256_setzero_si256())));
                            ss7 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss15, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss16, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss17, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss18, _mm256_setzero_si256())));
                            ss0 = _mm256_add_epi32(ss0, ss2);
                            ss1 = _mm256_add_epi32(ss1, ss3);
                            ss2 = _mm256_add_epi32(ss4, ss6);
                            ss3 = _mm256_add_epi32(ss5, ss7);
                            ss0 = _mm256_add_epi32(ss0, ss2);
                            ss1 = _mm256_add_epi32(ss1, ss3);
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss0),
                                               _mm256_extracti128_si256(ss0, 1));
                            s3  = _mm_add_epi32(_mm256_castsi256_si128(ss1),
                                               _mm256_extracti128_si256(ss1, 1));
                            UPDATE_BEST_PME(s0, 0, 0);
                            UPDATE_BEST_PME(s0, 1, 0);
                            UPDATE_BEST_PME(s0, 2, 0);
                            UPDATE_BEST_PME(s0, 3, 0);
                            UPDATE_BEST_PME(s3, 0, 4);
                            UPDATE_BEST_PME(s3, 1, 4);
                            UPDATE_BEST_PME(s3, 2, 4);
                            UPDATE_BEST_PME(s3, 3, 4);
                        }
                    }
                }
                ref += ref_stride;
            }
        } else {
            __m256i ssa1, ssa2, ssa3, ssa4, ssa5, ssa6, ssa7, ssa8, ssa9, ssa10, ssa11, ssa12,
                ssa13, ssa14, ssa15, ssa16;
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    ss3 = ss4 = ss5 = ss6 = ss7 = ss8 = ss9 = ss10 = ss11 = ss12 = ss13 = ss14 =
                        ss15 = ss16 = ss17 = ss18 = _mm256_setzero_si256();
                    for (k = 0; k < block_height; k++) {
                        ss0 = _mm256_loadu_si256((__m256i *)p_ref);
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 8));
                        ss2 = _mm256_loadu_si256((__m256i *)p_src);
                        ss3 = _mm256_adds_epu16(ss3, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss4 = _mm256_adds_epu16(ss4, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss5 = _mm256_adds_epu16(ss5, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss6 = _mm256_adds_epu16(ss6, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        ss0 = _mm256_loadu_si256((__m256i *)(p_ref + 32));
                        ss1 = _mm256_loadu_si256((__m256i *)(p_ref + 40));
                        ss2 = _mm256_loadu_si256((__m256i *)(p_src + 32));
                        ss7 = _mm256_adds_epu16(ss7, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss8 = _mm256_adds_epu16(ss8, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss9 = _mm256_adds_epu16(ss9, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss10 =
                            _mm256_adds_epu16(ss10, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        ss0  = _mm256_loadu_si256((__m256i *)(p_ref + 64));
                        ss1  = _mm256_loadu_si256((__m256i *)(p_ref + 72));
                        ss2  = _mm256_loadu_si256((__m256i *)(p_src + 64));
                        ss11 = _mm256_adds_epu16(ss11, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss12 =
                            _mm256_adds_epu16(ss12, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss13 =
                            _mm256_adds_epu16(ss13, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss14 =
                            _mm256_adds_epu16(ss14, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        ss0  = _mm256_loadu_si256((__m256i *)(p_ref + 96));
                        ss1  = _mm256_loadu_si256((__m256i *)(p_ref + 104));
                        ss2  = _mm256_loadu_si256((__m256i *)(p_src + 96));
                        ss15 = _mm256_adds_epu16(ss15, _mm256_mpsadbw_epu8(ss0, ss2, 0));
                        ss16 =
                            _mm256_adds_epu16(ss16, _mm256_mpsadbw_epu8(ss0, ss2, 45)); // 101 101
                        ss17 =
                            _mm256_adds_epu16(ss17, _mm256_mpsadbw_epu8(ss1, ss2, 18)); // 010 010
                        ss18 =
                            _mm256_adds_epu16(ss18, _mm256_mpsadbw_epu8(ss1, ss2, 63)); // 111 111
                        p_src += src_stride;
                        p_ref += ref_stride;
                        if (k == 63) {
                            ssa1  = ss3;
                            ssa2  = ss4;
                            ssa3  = ss5;
                            ssa4  = ss6;
                            ssa5  = ss7;
                            ssa6  = ss8;
                            ssa7  = ss9;
                            ssa8  = ss10;
                            ssa9  = ss11;
                            ssa10 = ss12;
                            ssa11 = ss13;
                            ssa12 = ss14;
                            ssa13 = ss15;
                            ssa14 = ss16;
                            ssa15 = ss17;
                            ssa16 = ss18;
                            ss3 = ss4 = ss5 = ss6 = ss7 = ss8 = ss9 = ss10 = ss11 = ss12 = ss13 =
                                ss14 = ss15 = ss16 = ss17 = ss18 = _mm256_setzero_si256();
                        }
                    }
                    ss0 =
                        _mm256_adds_epu16(_mm256_adds_epu16(ss3, ss4), _mm256_adds_epu16(ss5, ss6));
                    ss0 = _mm256_adds_epu16(ss0,
                                            _mm256_adds_epu16(_mm256_adds_epu16(ss7, ss8),
                                                              _mm256_adds_epu16(ss9, ss10)));
                    ss0 = _mm256_adds_epu16(ss0,
                                            _mm256_adds_epu16(_mm256_adds_epu16(ss11, ss12),
                                                              _mm256_adds_epu16(ss13, ss14)));
                    ss0 = _mm256_adds_epu16(ss0,
                                            _mm256_adds_epu16(_mm256_adds_epu16(ss15, ss16),
                                                              _mm256_adds_epu16(ss17, ss18)));
                    ss1 = _mm256_adds_epu16(_mm256_adds_epu16(ssa1, ssa2),
                                            _mm256_adds_epu16(ssa3, ssa4));
                    ss1 = _mm256_adds_epu16(ss1,
                                            _mm256_adds_epu16(_mm256_adds_epu16(ssa5, ssa6),
                                                              _mm256_adds_epu16(ssa7, ssa8)));
                    ss1 = _mm256_adds_epu16(ss1,
                                            _mm256_adds_epu16(_mm256_adds_epu16(ssa9, ssa10),
                                                              _mm256_adds_epu16(ssa11, ssa12)));
                    ss1 = _mm256_adds_epu16(ss1,
                                            _mm256_adds_epu16(_mm256_adds_epu16(ssa13, ssa14),
                                                              _mm256_adds_epu16(ssa15, ssa16)));
                    ss0 = _mm256_adds_epu16(ss0, ss1);

                    s0      = _mm_adds_epu16(_mm256_castsi256_si128(ss0),
                                        _mm256_extracti128_si256(ss0, 1));
                    s0      = _mm_minpos_epu16(s0);
                    tem_sum_1 = _mm_extract_epi16(s0, 0);
                    if (tem_sum_1 < low_sum) {
                        if (tem_sum_1 != 0xFFFF) { // no overflow
                            low_sum = tem_sum_1;
                            x_best  = mvx + (search_position_start_x +
                                           (int16_t)(j + _mm_extract_epi16(s0, 1))) *
                                              search_step;
                            y_best = mvy + (search_position_start_y + i) * search_step;
                        } else {
                            ss0 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss3, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss4, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss5, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss6, _mm256_setzero_si256())));
                            ss1 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss3, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss4, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss5, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss6, _mm256_setzero_si256())));
                            ss2 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss7, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss8, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss9, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss10, _mm256_setzero_si256())));
                            ss3 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss7, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss8, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss9, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss10, _mm256_setzero_si256())));
                            ss4 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss11, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss12, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss13, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss14, _mm256_setzero_si256())));
                            ss5 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss11, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss12, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss13, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss14, _mm256_setzero_si256())));
                            ss6 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss15, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss16, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ss17, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ss18, _mm256_setzero_si256())));
                            ss7 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss15, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss16, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ss17, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ss18, _mm256_setzero_si256())));
                            ss8 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ssa1, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ssa2, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ssa3, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ssa4, _mm256_setzero_si256())));
                            ss9 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ssa1, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ssa2, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ssa3, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ssa4, _mm256_setzero_si256())));
                            ss10 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ssa5, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ssa6, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ssa7, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ssa8, _mm256_setzero_si256())));
                            ss11 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ssa5, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ssa6, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ssa7, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ssa8, _mm256_setzero_si256())));
                            ss12 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ssa9, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ssa10, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ssa11, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ssa12, _mm256_setzero_si256())));
                            ss13 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ssa9, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ssa10, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ssa11, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ssa12, _mm256_setzero_si256())));
                            ss14 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ssa13, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ssa14, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpacklo_epi16(ssa15, _mm256_setzero_si256()),
                                    _mm256_unpacklo_epi16(ssa16, _mm256_setzero_si256())));
                            ss15 = _mm256_add_epi32(
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ssa13, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ssa14, _mm256_setzero_si256())),
                                _mm256_add_epi32(
                                    _mm256_unpackhi_epi16(ssa15, _mm256_setzero_si256()),
                                    _mm256_unpackhi_epi16(ssa16, _mm256_setzero_si256())));
                            ss0 = _mm256_add_epi32(ss0, ss2);
                            ss1 = _mm256_add_epi32(ss1, ss3);
                            ss2 = _mm256_add_epi32(ss4, ss6);
                            ss3 = _mm256_add_epi32(ss5, ss7);
                            ss4 = _mm256_add_epi32(ss8, ss10);
                            ss5 = _mm256_add_epi32(ss9, ss11);
                            ss6 = _mm256_add_epi32(ss12, ss14);
                            ss7 = _mm256_add_epi32(ss13, ss15);
                            ss0 = _mm256_add_epi32(ss0, ss2);
                            ss1 = _mm256_add_epi32(ss1, ss3);
                            ss2 = _mm256_add_epi32(ss4, ss6);
                            ss3 = _mm256_add_epi32(ss5, ss7);
                            ss0 = _mm256_add_epi32(ss0, ss2);
                            ss1 = _mm256_add_epi32(ss1, ss3);
                            s0  = _mm_add_epi32(_mm256_castsi256_si128(ss0),
                                               _mm256_extracti128_si256(ss0, 1));
                            s3  = _mm_add_epi32(_mm256_castsi256_si128(ss1),
                                               _mm256_extracti128_si256(ss1, 1));
                            UPDATE_BEST_PME(s0, 0, 0);
                            UPDATE_BEST_PME(s0, 1, 0);
                            UPDATE_BEST_PME(s0, 2, 0);
                            UPDATE_BEST_PME(s0, 3, 0);
                            UPDATE_BEST_PME(s3, 0, 4);
                            UPDATE_BEST_PME(s3, 1, 4);
                            UPDATE_BEST_PME(s3, 2, 4);
                            UPDATE_BEST_PME(s3, 3, 4);
                        }
                    }
                }
                ref += ref_stride;
            }
        }
        break;

    default: assert(0); break;
    }

    *best_sad = low_sum;
    *best_mvx = x_best;
    *best_mvy = y_best;
}
#endif

uint32_t eb_aom_sad4x4_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr,
                            int ref_stride) {
    return compute4x_m_sad_avx2(src_ptr, src_stride, ref_ptr, ref_stride, 4);
}

uint32_t eb_aom_sad4x8_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr,
                            int ref_stride) {
    return compute4x_m_sad_avx2(src_ptr, src_stride, ref_ptr, ref_stride, 8);
}

uint32_t eb_aom_sad4x16_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr,
                             int ref_stride) {
    return compute4x_m_sad_avx2(src_ptr, src_stride, ref_ptr, ref_stride, 16);
}

uint32_t eb_aom_sad8x4_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr,
                            int ref_stride) {
    return compute8x_m_sad_avx2(src_ptr, src_stride, ref_ptr, ref_stride, 4);
}

uint32_t eb_aom_sad8x8_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr,
                            int ref_stride) {
    return compute8x_m_sad_avx2(src_ptr, src_stride, ref_ptr, ref_stride, 8);
}

uint32_t eb_aom_sad8x16_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr,
                             int ref_stride) {
    return compute8x_m_sad_avx2(src_ptr, src_stride, ref_ptr, ref_stride, 16);
}

uint32_t eb_aom_sad8x32_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr,
                             int ref_stride) {
    return compute8x_m_sad_avx2(src_ptr, src_stride, ref_ptr, ref_stride, 32);
}

uint32_t eb_aom_sad16x4_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr,
                             int ref_stride) {
    return compute16x_m_sad_avx2(src_ptr, src_stride, ref_ptr, ref_stride, 4);
}

uint32_t eb_aom_sad16x8_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr,
                             int ref_stride) {
    return compute16x_m_sad_avx2(src_ptr, src_stride, ref_ptr, ref_stride, 8);
}

uint32_t eb_aom_sad16x16_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr,
                              int ref_stride) {
    return compute16x_m_sad_avx2(src_ptr, src_stride, ref_ptr, ref_stride, 16);
}

uint32_t eb_aom_sad16x32_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr,
                              int ref_stride) {
    return compute16x_m_sad_avx2(src_ptr, src_stride, ref_ptr, ref_stride, 32);
}

uint32_t eb_aom_sad16x64_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr,
                              int ref_stride) {
    return compute16x_m_sad_avx2(src_ptr, src_stride, ref_ptr, ref_stride, 64);
}

uint32_t eb_aom_sad32x8_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr,
                             int ref_stride) {
    return compute32x_m_sad_avx2(src_ptr, src_stride, ref_ptr, ref_stride, 8);
}

uint32_t eb_aom_sad32x16_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr,
                              int ref_stride) {
    return compute32x_m_sad_avx2(src_ptr, src_stride, ref_ptr, ref_stride, 16);
}

uint32_t eb_aom_sad32x32_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr,
                              int ref_stride) {
    return compute32x_m_sad_avx2(src_ptr, src_stride, ref_ptr, ref_stride, 32);
}

uint32_t eb_aom_sad32x64_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr,
                              int ref_stride) {
    return compute32x_m_sad_avx2(src_ptr, src_stride, ref_ptr, ref_stride, 64);
}

uint32_t eb_aom_sad64x16_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr,
                              int ref_stride) {
    return compute64x_m_sad_avx2(src_ptr, src_stride, ref_ptr, ref_stride, 16);
}

uint32_t eb_aom_sad64x32_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr,
                              int ref_stride) {
    return compute64x_m_sad_avx2(src_ptr, src_stride, ref_ptr, ref_stride, 32);
}

uint32_t eb_aom_sad64x64_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr,
                              int ref_stride) {
    return compute64x_m_sad_avx2(src_ptr, src_stride, ref_ptr, ref_stride, 64);
}

uint32_t eb_aom_sad128x64_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr,
                               int ref_stride) {
    return compute128x_m_sad_avx2(src_ptr, src_stride, ref_ptr, ref_stride, 64);
}

uint32_t eb_aom_sad128x128_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr,
                                int ref_stride) {
    return compute128x_m_sad_avx2(src_ptr, src_stride, ref_ptr, ref_stride, 128);
}

SIMD_INLINE void eb_aom_sad4xhx4d_calc_avx2(const uint8_t *src, int src_stride,
                                            const uint8_t *const ref[4], int ref_stride,
                                            uint32_t res[4], uint32_t height) {
    uint32_t       y;
    uint32_t       hsteps = height >> 2;
    const uint8_t *ref0, *ref1, *ref2, *ref3;
    assert(!(height % 4));

    ref0 = ref[0];
    ref1 = ref[1];
    ref2 = ref[2];
    ref3 = ref[3];

    __m256i mm256_sad0 = _mm256_setzero_si256();
    __m256i mm256_sad1 = _mm256_setzero_si256();
    __m256i mm256_sad2 = _mm256_setzero_si256();
    __m256i mm256_sad3 = _mm256_setzero_si256();
    __m256i mm256_src;
    __m256i mm256_ref0;
    __m256i mm256_ref1;
    __m256i mm256_ref2;
    __m256i mm256_ref3;

    for (y = 0; y < hsteps; y++) {
        mm256_src = _mm256_setr_epi32(*(uint32_t *)src,
                                      *(uint32_t *)(src + 1 * src_stride),
                                      *(uint32_t *)(src + 2 * src_stride),
                                      *(uint32_t *)(src + 3 * src_stride),
                                      0,
                                      0,
                                      0,
                                      0);

        mm256_ref0 = _mm256_setr_epi32(*(uint32_t *)ref0,
                                       *(uint32_t *)(ref0 + 1 * ref_stride),
                                       *(uint32_t *)(ref0 + 2 * ref_stride),
                                       *(uint32_t *)(ref0 + 3 * ref_stride),
                                       0,
                                       0,
                                       0,
                                       0);

        mm256_ref1 = _mm256_setr_epi32(*(uint32_t *)ref1,
                                       *(uint32_t *)(ref1 + 1 * ref_stride),
                                       *(uint32_t *)(ref1 + 2 * ref_stride),
                                       *(uint32_t *)(ref1 + 3 * ref_stride),
                                       0,
                                       0,
                                       0,
                                       0);

        mm256_ref2 = _mm256_setr_epi32(*(uint32_t *)ref2,
                                       *(uint32_t *)(ref2 + 1 * ref_stride),
                                       *(uint32_t *)(ref2 + 2 * ref_stride),
                                       *(uint32_t *)(ref2 + 3 * ref_stride),
                                       0,
                                       0,
                                       0,
                                       0);

        mm256_ref3 = _mm256_setr_epi32(*(uint32_t *)ref3,
                                       *(uint32_t *)(ref3 + 1 * ref_stride),
                                       *(uint32_t *)(ref3 + 2 * ref_stride),
                                       *(uint32_t *)(ref3 + 3 * ref_stride),
                                       0,
                                       0,
                                       0,
                                       0);

        mm256_sad0 = _mm256_add_epi16(mm256_sad0, _mm256_sad_epu8(mm256_src, mm256_ref0));

        mm256_sad1 = _mm256_add_epi16(mm256_sad1, _mm256_sad_epu8(mm256_src, mm256_ref1));

        mm256_sad2 = _mm256_add_epi16(mm256_sad2, _mm256_sad_epu8(mm256_src, mm256_ref2));

        mm256_sad3 = _mm256_add_epi16(mm256_sad3, _mm256_sad_epu8(mm256_src, mm256_ref3));

        src += src_stride << 2;
        ref0 += ref_stride << 2;
        ref1 += ref_stride << 2;
        ref2 += ref_stride << 2;
        ref3 += ref_stride << 2;
    }

    res[0] = _mm256_extract_epi32(mm256_sad0, 0) + _mm256_extract_epi32(mm256_sad0, 2);
    res[1] = _mm256_extract_epi32(mm256_sad1, 0) + _mm256_extract_epi32(mm256_sad1, 2);
    res[2] = _mm256_extract_epi32(mm256_sad2, 0) + _mm256_extract_epi32(mm256_sad2, 2);
    res[3] = _mm256_extract_epi32(mm256_sad3, 0) + _mm256_extract_epi32(mm256_sad3, 2);
}

void eb_aom_sad4x4x4d_avx2(const uint8_t *src, int src_stride, const uint8_t *const ref[4],
                           int ref_stride, uint32_t res[4]) {
    eb_aom_sad4xhx4d_calc_avx2(src, src_stride, ref, ref_stride, res, 4);
}

void eb_aom_sad4x8x4d_avx2(const uint8_t *src, int src_stride, const uint8_t *const ref[4],
                           int ref_stride, uint32_t res[4]) {
    eb_aom_sad4xhx4d_calc_avx2(src, src_stride, ref, ref_stride, res, 8);
}

void eb_aom_sad4x16x4d_avx2(const uint8_t *src, int src_stride, const uint8_t *const ref[4],
                            int ref_stride, uint32_t res[4]) {
    eb_aom_sad4xhx4d_calc_avx2(src, src_stride, ref, ref_stride, res, 16);
}

SIMD_INLINE void eb_aom_sad8xhx4d_calc_avx2(const uint8_t *src, int src_stride,
                                            const uint8_t *const ref[4], int ref_stride,
                                            uint32_t res[4], uint32_t height) {
    __m128i  xmm0;
    __m256i  ymm0 = _mm256_setzero_si256();
    __m256i  ymm1 = _mm256_setzero_si256();
    __m256i  ymm2 = _mm256_setzero_si256();
    __m256i  ymm3 = _mm256_setzero_si256();

    const uint8_t *ref0, *ref1, *ref2, *ref3;

    ref0 = ref[0];
    ref1 = ref[1];
    ref2 = ref[2];
    ref3 = ref[3];

    for (uint32_t y = 0; y < height; y += 4) {
        __m128i src01, src23;
        src01 = _mm_loadl_epi64((__m128i *)(src + 0 * src_stride));
        src01 = _mm_castpd_si128(_mm_loadh_pd(
            _mm_castsi128_pd(src01), (const double *)(const void *)(src + 1 * src_stride)));
        src23 = _mm_loadl_epi64((__m128i *)(src + 2 * src_stride));
        src23                 = _mm_castpd_si128(_mm_loadh_pd(
            _mm_castsi128_pd(src23), (const double *)(const void *)(src + 3 * src_stride)));
        const __m256i src0123 = _mm256_setr_m128i(src01, src23);

        src01 = _mm_loadl_epi64((__m128i *)(ref0 + 0 * ref_stride));
        src01 = _mm_castpd_si128(
            _mm_loadh_pd(_mm_castsi128_pd(src01), (const double *)(const void*)(ref0 + 1 * ref_stride)));
        src23 = _mm_loadl_epi64((__m128i *)(ref0 + 2 * ref_stride));
        src23 = _mm_castpd_si128(
            _mm_loadh_pd(_mm_castsi128_pd(src23), (const double *)(const void*)(ref0 + 3 * ref_stride)));
        const __m256i ref0123_0 = _mm256_setr_m128i(src01, src23);

        src01 = _mm_loadl_epi64((__m128i *)(ref1 + 0 * ref_stride));
        src01 = _mm_castpd_si128(
            _mm_loadh_pd(_mm_castsi128_pd(src01), (const double *)(const void*)(ref1 + 1 * ref_stride)));
        src23 = _mm_loadl_epi64((__m128i *)(ref1 + 2 * ref_stride));
        src23 = _mm_castpd_si128(
            _mm_loadh_pd(_mm_castsi128_pd(src23), (const double *)(const void*)(ref1 + 3 * ref_stride)));
        const __m256i ref0123_1 = _mm256_setr_m128i(src01, src23);

        src01 = _mm_loadl_epi64((__m128i *)(ref2 + 0 * ref_stride));
        src01 = _mm_castpd_si128(
            _mm_loadh_pd(_mm_castsi128_pd(src01), (const double *)(const void*)(ref2 + 1 * ref_stride)));
        src23 = _mm_loadl_epi64((__m128i *)(ref2 + 2 * ref_stride));
        src23 = _mm_castpd_si128(
            _mm_loadh_pd(_mm_castsi128_pd(src23), (const double *)(const void*)(ref2 + 3 * ref_stride)));
        const __m256i ref0123_2 = _mm256_setr_m128i(src01, src23);

        src01 = _mm_loadl_epi64((__m128i *)(ref3 + 0 * ref_stride));
        src01 = _mm_castpd_si128(
            _mm_loadh_pd(_mm_castsi128_pd(src01), (const double *)(const void*)(ref3 + 1 * ref_stride)));
        src23 = _mm_loadl_epi64((__m128i *)(ref3 + 2 * ref_stride));
        src23 = _mm_castpd_si128(
            _mm_loadh_pd(_mm_castsi128_pd(src23), (const double *)(const void*)(ref3 + 3 * ref_stride)));
        const __m256i ref0123_3 = _mm256_setr_m128i(src01, src23);

        ymm0 = _mm256_add_epi32(ymm0, _mm256_sad_epu8(src0123, ref0123_0));
        ymm1 = _mm256_add_epi32(ymm1, _mm256_sad_epu8(src0123, ref0123_1));
        ymm2 = _mm256_add_epi32(ymm2, _mm256_sad_epu8(src0123, ref0123_2));
        ymm3 = _mm256_add_epi32(ymm3, _mm256_sad_epu8(src0123, ref0123_3));
        src += src_stride << 2;
        ref0 += ref_stride << 2;
        ref1 += ref_stride << 2;
        ref2 += ref_stride << 2;
        ref3 += ref_stride << 2;
    }

    xmm0   = _mm_add_epi32(_mm256_castsi256_si128(ymm0), _mm256_extracti128_si256(ymm0, 1));
    xmm0   = _mm_add_epi32(xmm0, _mm_srli_si128(xmm0, 8));
    res[0] = (uint32_t)_mm_cvtsi128_si32(xmm0);

    xmm0   = _mm_add_epi32(_mm256_castsi256_si128(ymm1), _mm256_extracti128_si256(ymm1, 1));
    xmm0   = _mm_add_epi32(xmm0, _mm_srli_si128(xmm0, 8));
    res[1] = (uint32_t)_mm_cvtsi128_si32(xmm0);

    xmm0   = _mm_add_epi32(_mm256_castsi256_si128(ymm2), _mm256_extracti128_si256(ymm2, 1));
    xmm0   = _mm_add_epi32(xmm0, _mm_srli_si128(xmm0, 8));
    res[2] = (uint32_t)_mm_cvtsi128_si32(xmm0);

    xmm0   = _mm_add_epi32(_mm256_castsi256_si128(ymm3), _mm256_extracti128_si256(ymm3, 1));
    xmm0   = _mm_add_epi32(xmm0, _mm_srli_si128(xmm0, 8));
    res[3] = (uint32_t)_mm_cvtsi128_si32(xmm0);
}

void eb_aom_sad8x4x4d_avx2(const uint8_t *src, int src_stride, const uint8_t *const ref[4],
                           int ref_stride, uint32_t res[4]) {
    eb_aom_sad8xhx4d_calc_avx2(src, src_stride, ref, ref_stride, res, 4);
}

void eb_aom_sad8x8x4d_avx2(const uint8_t *src, int src_stride, const uint8_t *const ref[4],
                           int ref_stride, uint32_t res[4]) {
    eb_aom_sad8xhx4d_calc_avx2(src, src_stride, ref, ref_stride, res, 8);
}

void eb_aom_sad8x16x4d_avx2(const uint8_t *src, int src_stride, const uint8_t *const ref[4],
                            int ref_stride, uint32_t res[4]) {
    eb_aom_sad8xhx4d_calc_avx2(src, src_stride, ref, ref_stride, res, 16);
}

void eb_aom_sad8x32x4d_avx2(const uint8_t *src, int src_stride, const uint8_t *const ref[4],
                            int ref_stride, uint32_t res[4]) {
    eb_aom_sad8xhx4d_calc_avx2(src, src_stride, ref, ref_stride, res, 32);
}

SIMD_INLINE void eb_aom_sad16xhx4d_calc_avx2(const uint8_t *src, int src_stride,
                                             const uint8_t *const ref[4], int ref_stride,
                                             uint32_t res[4], uint32_t height) {
    __m128i  xmm0;
    __m256i  ymm0 = _mm256_setzero_si256();
    __m256i  ymm1 = _mm256_setzero_si256();
    __m256i  ymm2 = _mm256_setzero_si256();
    __m256i  ymm3 = _mm256_setzero_si256();
    uint32_t y;

    const uint8_t *ref0, *ref1, *ref2, *ref3;

    ref0 = ref[0];
    ref1 = ref[1];
    ref2 = ref[2];
    ref3 = ref[3];

    for (y = 0; y < height; y += 2) {
        __m256i src0 = _mm256_setr_m128i(_mm_loadu_si128((__m128i *)(src + 0 * src_stride)),
                                         _mm_loadu_si128((__m128i *)(src + 1 * src_stride)));

        __m256i r0 = _mm256_setr_m128i(_mm_loadu_si128((__m128i *)(ref0 + 0 * ref_stride)),
                                       _mm_loadu_si128((__m128i *)(ref0 + 1 * ref_stride)));

        __m256i r1 = _mm256_setr_m128i(_mm_loadu_si128((__m128i *)(ref1 + 0 * ref_stride)),
                                       _mm_loadu_si128((__m128i *)(ref1 + 1 * ref_stride)));

        __m256i r2 = _mm256_setr_m128i(_mm_loadu_si128((__m128i *)(ref2 + 0 * ref_stride)),
                                       _mm_loadu_si128((__m128i *)(ref2 + 1 * ref_stride)));

        __m256i r3 = _mm256_setr_m128i(_mm_loadu_si128((__m128i *)(ref3 + 0 * ref_stride)),
                                       _mm_loadu_si128((__m128i *)(ref3 + 1 * ref_stride)));

        ymm0 = _mm256_add_epi32(ymm0, _mm256_sad_epu8(src0, r0));
        ymm1 = _mm256_add_epi32(ymm1, _mm256_sad_epu8(src0, r1));
        ymm2 = _mm256_add_epi32(ymm2, _mm256_sad_epu8(src0, r2));
        ymm3 = _mm256_add_epi32(ymm3, _mm256_sad_epu8(src0, r3));

        src += src_stride << 1;
        ref0 += ref_stride << 1;
        ref1 += ref_stride << 1;
        ref2 += ref_stride << 1;
        ref3 += ref_stride << 1;
    }

    xmm0   = _mm_add_epi32(_mm256_castsi256_si128(ymm0), _mm256_extracti128_si256(ymm0, 1));
    xmm0   = _mm_add_epi32(xmm0, _mm_srli_si128(xmm0, 8));
    res[0] = (uint32_t)_mm_cvtsi128_si32(xmm0);

    xmm0   = _mm_add_epi32(_mm256_castsi256_si128(ymm1), _mm256_extracti128_si256(ymm1, 1));
    xmm0   = _mm_add_epi32(xmm0, _mm_srli_si128(xmm0, 8));
    res[1] = (uint32_t)_mm_cvtsi128_si32(xmm0);

    xmm0   = _mm_add_epi32(_mm256_castsi256_si128(ymm2), _mm256_extracti128_si256(ymm2, 1));
    xmm0   = _mm_add_epi32(xmm0, _mm_srli_si128(xmm0, 8));
    res[2] = (uint32_t)_mm_cvtsi128_si32(xmm0);

    xmm0   = _mm_add_epi32(_mm256_castsi256_si128(ymm3), _mm256_extracti128_si256(ymm3, 1));
    xmm0   = _mm_add_epi32(xmm0, _mm_srli_si128(xmm0, 8));
    res[3] = (uint32_t)_mm_cvtsi128_si32(xmm0);
}

void eb_aom_sad16x4x4d_avx2(const uint8_t *src, int src_stride, const uint8_t *const ref[4],
                            int ref_stride, uint32_t res[4]) {
    eb_aom_sad16xhx4d_calc_avx2(src, src_stride, ref, ref_stride, res, 4);
}

void eb_aom_sad16x8x4d_avx2(const uint8_t *src, int src_stride, const uint8_t *const ref[4],
                            int ref_stride, uint32_t res[4]) {
    eb_aom_sad16xhx4d_calc_avx2(src, src_stride, ref, ref_stride, res, 8);
}

void eb_aom_sad16x16x4d_avx2(const uint8_t *src, int src_stride, const uint8_t *const ref[4],
                             int ref_stride, uint32_t res[4]) {
    eb_aom_sad16xhx4d_calc_avx2(src, src_stride, ref, ref_stride, res, 16);
}

void eb_aom_sad16x32x4d_avx2(const uint8_t *src, int src_stride, const uint8_t *const ref[4],
                             int ref_stride, uint32_t res[4]) {
    eb_aom_sad16xhx4d_calc_avx2(src, src_stride, ref, ref_stride, res, 32);
}

void eb_aom_sad16x64x4d_avx2(const uint8_t *src, int src_stride, const uint8_t *const ref[4],
                             int ref_stride, uint32_t res[4]) {
    eb_aom_sad16xhx4d_calc_avx2(src, src_stride, ref, ref_stride, res, 64);
}

void eb_aom_sad32x8x4d_avx2(const uint8_t *src, int src_stride, const uint8_t *const ref[4],
                            int ref_stride, uint32_t res[4]) {
    __m256i        src_reg, ref0_reg, ref1_reg, ref2_reg, ref3_reg;
    __m256i        sum_ref0, sum_ref1, sum_ref2, sum_ref3;
    __m256i        sum_mlow, sum_mhigh;
    int            i;
    const uint8_t *ref0, *ref1, *ref2, *ref3;

    ref0     = ref[0];
    ref1     = ref[1];
    ref2     = ref[2];
    ref3     = ref[3];
    sum_ref0 = _mm256_set1_epi16(0);
    sum_ref1 = _mm256_set1_epi16(0);
    sum_ref2 = _mm256_set1_epi16(0);
    sum_ref3 = _mm256_set1_epi16(0);
    for (i = 0; i < 8; i++) {
        // load src and all refs
        src_reg  = _mm256_loadu_si256((const __m256i *)src);
        ref0_reg = _mm256_loadu_si256((const __m256i *)ref0);
        ref1_reg = _mm256_loadu_si256((const __m256i *)ref1);
        ref2_reg = _mm256_loadu_si256((const __m256i *)ref2);
        ref3_reg = _mm256_loadu_si256((const __m256i *)ref3);
        // sum of the absolute differences between every ref-i to src
        ref0_reg = _mm256_sad_epu8(ref0_reg, src_reg);
        ref1_reg = _mm256_sad_epu8(ref1_reg, src_reg);
        ref2_reg = _mm256_sad_epu8(ref2_reg, src_reg);
        ref3_reg = _mm256_sad_epu8(ref3_reg, src_reg);
        // sum every ref-i
        sum_ref0 = _mm256_add_epi32(sum_ref0, ref0_reg);
        sum_ref1 = _mm256_add_epi32(sum_ref1, ref1_reg);
        sum_ref2 = _mm256_add_epi32(sum_ref2, ref2_reg);
        sum_ref3 = _mm256_add_epi32(sum_ref3, ref3_reg);

        src += src_stride;
        ref0 += ref_stride;
        ref1 += ref_stride;
        ref2 += ref_stride;
        ref3 += ref_stride;
    }
    {
        __m128i sum;
        // in sum_ref-i the result is saved in the first 4 bytes
        // the other 4 bytes are zeroed.
        // sum_ref1 and sum_ref3 are shifted left by 4 bytes
        sum_ref1 = _mm256_slli_si256(sum_ref1, 4);
        sum_ref3 = _mm256_slli_si256(sum_ref3, 4);

        // merge sum_ref0 and sum_ref1 also sum_ref2 and sum_ref3
        sum_ref0 = _mm256_or_si256(sum_ref0, sum_ref1);
        sum_ref2 = _mm256_or_si256(sum_ref2, sum_ref3);

        // merge every 64 bit from each sum_ref-i
        sum_mlow  = _mm256_unpacklo_epi64(sum_ref0, sum_ref2);
        sum_mhigh = _mm256_unpackhi_epi64(sum_ref0, sum_ref2);

        // add the low 64 bit to the high 64 bit
        sum_mlow = _mm256_add_epi32(sum_mlow, sum_mhigh);

        // add the low 128 bit to the high 128 bit
        sum =
            _mm_add_epi32(_mm256_castsi256_si128(sum_mlow), _mm256_extractf128_si256(sum_mlow, 1));

        _mm_storeu_si128((__m128i *)(res), sum);
    }
    _mm256_zeroupper();
}

void eb_aom_sad32x16x4d_avx2(const uint8_t *src, int src_stride, const uint8_t *const ref[4],
                             int ref_stride, uint32_t res[4]) {
    __m256i        src_reg, ref0_reg, ref1_reg, ref2_reg, ref3_reg;
    __m256i        sum_ref0, sum_ref1, sum_ref2, sum_ref3;
    __m256i        sum_mlow, sum_mhigh;
    int            i;
    const uint8_t *ref0, *ref1, *ref2, *ref3;

    ref0     = ref[0];
    ref1     = ref[1];
    ref2     = ref[2];
    ref3     = ref[3];
    sum_ref0 = _mm256_set1_epi16(0);
    sum_ref1 = _mm256_set1_epi16(0);
    sum_ref2 = _mm256_set1_epi16(0);
    sum_ref3 = _mm256_set1_epi16(0);
    for (i = 0; i < 16; i++) {
        // load src and all refs
        src_reg  = _mm256_loadu_si256((const __m256i *)src);
        ref0_reg = _mm256_loadu_si256((const __m256i *)ref0);
        ref1_reg = _mm256_loadu_si256((const __m256i *)ref1);
        ref2_reg = _mm256_loadu_si256((const __m256i *)ref2);
        ref3_reg = _mm256_loadu_si256((const __m256i *)ref3);
        // sum of the absolute differences between every ref-i to src
        ref0_reg = _mm256_sad_epu8(ref0_reg, src_reg);
        ref1_reg = _mm256_sad_epu8(ref1_reg, src_reg);
        ref2_reg = _mm256_sad_epu8(ref2_reg, src_reg);
        ref3_reg = _mm256_sad_epu8(ref3_reg, src_reg);
        // sum every ref-i
        sum_ref0 = _mm256_add_epi32(sum_ref0, ref0_reg);
        sum_ref1 = _mm256_add_epi32(sum_ref1, ref1_reg);
        sum_ref2 = _mm256_add_epi32(sum_ref2, ref2_reg);
        sum_ref3 = _mm256_add_epi32(sum_ref3, ref3_reg);

        src += src_stride;
        ref0 += ref_stride;
        ref1 += ref_stride;
        ref2 += ref_stride;
        ref3 += ref_stride;
    }
    {
        __m128i sum;
        // in sum_ref-i the result is saved in the first 4 bytes
        // the other 4 bytes are zeroed.
        // sum_ref1 and sum_ref3 are shifted left by 4 bytes
        sum_ref1 = _mm256_slli_si256(sum_ref1, 4);
        sum_ref3 = _mm256_slli_si256(sum_ref3, 4);

        // merge sum_ref0 and sum_ref1 also sum_ref2 and sum_ref3
        sum_ref0 = _mm256_or_si256(sum_ref0, sum_ref1);
        sum_ref2 = _mm256_or_si256(sum_ref2, sum_ref3);

        // merge every 64 bit from each sum_ref-i
        sum_mlow  = _mm256_unpacklo_epi64(sum_ref0, sum_ref2);
        sum_mhigh = _mm256_unpackhi_epi64(sum_ref0, sum_ref2);

        // add the low 64 bit to the high 64 bit
        sum_mlow = _mm256_add_epi32(sum_mlow, sum_mhigh);

        // add the low 128 bit to the high 128 bit
        sum =
            _mm_add_epi32(_mm256_castsi256_si128(sum_mlow), _mm256_extractf128_si256(sum_mlow, 1));

        _mm_storeu_si128((__m128i *)(res), sum);
    }
    _mm256_zeroupper();
}

static INLINE void sad32x32x4d_avx2(const uint8_t *src, int src_stride, const uint8_t *const ref[4],
                                    int ref_stride, uint32_t res[4]) {
    __m256i        src_reg, ref0_reg, ref1_reg, ref2_reg, ref3_reg;
    __m256i        sum_ref0, sum_ref1, sum_ref2, sum_ref3;
    __m256i        sum_mlow, sum_mhigh;
    int            i;
    const uint8_t *ref0, *ref1, *ref2, *ref3;

    ref0     = ref[0];
    ref1     = ref[1];
    ref2     = ref[2];
    ref3     = ref[3];
    sum_ref0 = _mm256_set1_epi16(0);
    sum_ref1 = _mm256_set1_epi16(0);
    sum_ref2 = _mm256_set1_epi16(0);
    sum_ref3 = _mm256_set1_epi16(0);
    for (i = 0; i < 32; i++) {
        // load src and all refs
        src_reg  = _mm256_loadu_si256((const __m256i *)src);
        ref0_reg = _mm256_loadu_si256((const __m256i *)ref0);
        ref1_reg = _mm256_loadu_si256((const __m256i *)ref1);
        ref2_reg = _mm256_loadu_si256((const __m256i *)ref2);
        ref3_reg = _mm256_loadu_si256((const __m256i *)ref3);
        // sum of the absolute differences between every ref-i to src
        ref0_reg = _mm256_sad_epu8(ref0_reg, src_reg);
        ref1_reg = _mm256_sad_epu8(ref1_reg, src_reg);
        ref2_reg = _mm256_sad_epu8(ref2_reg, src_reg);
        ref3_reg = _mm256_sad_epu8(ref3_reg, src_reg);
        // sum every ref-i
        sum_ref0 = _mm256_add_epi32(sum_ref0, ref0_reg);
        sum_ref1 = _mm256_add_epi32(sum_ref1, ref1_reg);
        sum_ref2 = _mm256_add_epi32(sum_ref2, ref2_reg);
        sum_ref3 = _mm256_add_epi32(sum_ref3, ref3_reg);

        src += src_stride;
        ref0 += ref_stride;
        ref1 += ref_stride;
        ref2 += ref_stride;
        ref3 += ref_stride;
    }
    {
        __m128i sum;
        // in sum_ref-i the result is saved in the first 4 bytes
        // the other 4 bytes are zeroed.
        // sum_ref1 and sum_ref3 are shifted left by 4 bytes
        sum_ref1 = _mm256_slli_si256(sum_ref1, 4);
        sum_ref3 = _mm256_slli_si256(sum_ref3, 4);

        // merge sum_ref0 and sum_ref1 also sum_ref2 and sum_ref3
        sum_ref0 = _mm256_or_si256(sum_ref0, sum_ref1);
        sum_ref2 = _mm256_or_si256(sum_ref2, sum_ref3);

        // merge every 64 bit from each sum_ref-i
        sum_mlow  = _mm256_unpacklo_epi64(sum_ref0, sum_ref2);
        sum_mhigh = _mm256_unpackhi_epi64(sum_ref0, sum_ref2);

        // add the low 64 bit to the high 64 bit
        sum_mlow = _mm256_add_epi32(sum_mlow, sum_mhigh);

        // add the low 128 bit to the high 128 bit
        sum =
            _mm_add_epi32(_mm256_castsi256_si128(sum_mlow), _mm256_extractf128_si256(sum_mlow, 1));

        _mm_storeu_si128((__m128i *)(res), sum);
    }
    _mm256_zeroupper();
}

void eb_aom_sad32x32x4d_avx2(const uint8_t *src, int src_stride, const uint8_t *const ref[4],
                             int ref_stride, uint32_t res[4]) {
    sad32x32x4d_avx2(src, src_stride, ref, ref_stride, res);
}

void eb_aom_sad64x16x4d_avx2(const uint8_t *src, int src_stride, const uint8_t *const ref[4],
                             int ref_stride, uint32_t res[4]) {
    __m256i        src_reg, srcnext_reg, ref0_reg, ref0next_reg;
    __m256i        ref1_reg, ref1next_reg, ref2_reg, ref2next_reg;
    __m256i        ref3_reg, ref3next_reg;
    __m256i        sum_ref0, sum_ref1, sum_ref2, sum_ref3;
    __m256i        sum_mlow, sum_mhigh;
    int            i;
    const uint8_t *ref0, *ref1, *ref2, *ref3;

    ref0     = ref[0];
    ref1     = ref[1];
    ref2     = ref[2];
    ref3     = ref[3];
    sum_ref0 = _mm256_set1_epi16(0);
    sum_ref1 = _mm256_set1_epi16(0);
    sum_ref2 = _mm256_set1_epi16(0);
    sum_ref3 = _mm256_set1_epi16(0);
    for (i = 0; i < 16; i++) {
        // load 64 bytes from src and all refs
        src_reg      = _mm256_loadu_si256((const __m256i *)src);
        srcnext_reg  = _mm256_loadu_si256((const __m256i *)(src + 32));
        ref0_reg     = _mm256_loadu_si256((const __m256i *)ref0);
        ref0next_reg = _mm256_loadu_si256((const __m256i *)(ref0 + 32));
        ref1_reg     = _mm256_loadu_si256((const __m256i *)ref1);
        ref1next_reg = _mm256_loadu_si256((const __m256i *)(ref1 + 32));
        ref2_reg     = _mm256_loadu_si256((const __m256i *)ref2);
        ref2next_reg = _mm256_loadu_si256((const __m256i *)(ref2 + 32));
        ref3_reg     = _mm256_loadu_si256((const __m256i *)ref3);
        ref3next_reg = _mm256_loadu_si256((const __m256i *)(ref3 + 32));
        // sum of the absolute differences between every ref-i to src
        ref0_reg     = _mm256_sad_epu8(ref0_reg, src_reg);
        ref1_reg     = _mm256_sad_epu8(ref1_reg, src_reg);
        ref2_reg     = _mm256_sad_epu8(ref2_reg, src_reg);
        ref3_reg     = _mm256_sad_epu8(ref3_reg, src_reg);
        ref0next_reg = _mm256_sad_epu8(ref0next_reg, srcnext_reg);
        ref1next_reg = _mm256_sad_epu8(ref1next_reg, srcnext_reg);
        ref2next_reg = _mm256_sad_epu8(ref2next_reg, srcnext_reg);
        ref3next_reg = _mm256_sad_epu8(ref3next_reg, srcnext_reg);

        // sum every ref-i
        sum_ref0 = _mm256_add_epi32(sum_ref0, ref0_reg);
        sum_ref1 = _mm256_add_epi32(sum_ref1, ref1_reg);
        sum_ref2 = _mm256_add_epi32(sum_ref2, ref2_reg);
        sum_ref3 = _mm256_add_epi32(sum_ref3, ref3_reg);
        sum_ref0 = _mm256_add_epi32(sum_ref0, ref0next_reg);
        sum_ref1 = _mm256_add_epi32(sum_ref1, ref1next_reg);
        sum_ref2 = _mm256_add_epi32(sum_ref2, ref2next_reg);
        sum_ref3 = _mm256_add_epi32(sum_ref3, ref3next_reg);
        src += src_stride;
        ref0 += ref_stride;
        ref1 += ref_stride;
        ref2 += ref_stride;
        ref3 += ref_stride;
    }
    {
        __m128i sum;

        // in sum_ref-i the result is saved in the first 4 bytes
        // the other 4 bytes are zeroed.
        // sum_ref1 and sum_ref3 are shifted left by 4 bytes
        sum_ref1 = _mm256_slli_si256(sum_ref1, 4);
        sum_ref3 = _mm256_slli_si256(sum_ref3, 4);

        // merge sum_ref0 and sum_ref1 also sum_ref2 and sum_ref3
        sum_ref0 = _mm256_or_si256(sum_ref0, sum_ref1);
        sum_ref2 = _mm256_or_si256(sum_ref2, sum_ref3);

        // merge every 64 bit from each sum_ref-i
        sum_mlow  = _mm256_unpacklo_epi64(sum_ref0, sum_ref2);
        sum_mhigh = _mm256_unpackhi_epi64(sum_ref0, sum_ref2);

        // add the low 64 bit to the high 64 bit
        sum_mlow = _mm256_add_epi32(sum_mlow, sum_mhigh);

        // add the low 128 bit to the high 128 bit
        sum =
            _mm_add_epi32(_mm256_castsi256_si128(sum_mlow), _mm256_extractf128_si256(sum_mlow, 1));

        _mm_storeu_si128((__m128i *)(res), sum);
    }
    _mm256_zeroupper();
}

SIMD_INLINE void sad64x64x4d_avx2(const uint8_t *src, int src_stride, const uint8_t *const ref[4],
                                  int ref_stride, uint32_t res[4]) {
    __m256i        src_reg, srcnext_reg, ref0_reg, ref0next_reg;
    __m256i        ref1_reg, ref1next_reg, ref2_reg, ref2next_reg;
    __m256i        ref3_reg, ref3next_reg;
    __m256i        sum_ref0, sum_ref1, sum_ref2, sum_ref3;
    __m256i        sum_mlow, sum_mhigh;
    int            i;
    const uint8_t *ref0, *ref1, *ref2, *ref3;

    ref0     = ref[0];
    ref1     = ref[1];
    ref2     = ref[2];
    ref3     = ref[3];
    sum_ref0 = _mm256_set1_epi16(0);
    sum_ref1 = _mm256_set1_epi16(0);
    sum_ref2 = _mm256_set1_epi16(0);
    sum_ref3 = _mm256_set1_epi16(0);
    for (i = 0; i < 64; i++) {
        // load 64 bytes from src and all refs
        src_reg      = _mm256_loadu_si256((const __m256i *)src);
        srcnext_reg  = _mm256_loadu_si256((const __m256i *)(src + 32));
        ref0_reg     = _mm256_loadu_si256((const __m256i *)ref0);
        ref0next_reg = _mm256_loadu_si256((const __m256i *)(ref0 + 32));
        ref1_reg     = _mm256_loadu_si256((const __m256i *)ref1);
        ref1next_reg = _mm256_loadu_si256((const __m256i *)(ref1 + 32));
        ref2_reg     = _mm256_loadu_si256((const __m256i *)ref2);
        ref2next_reg = _mm256_loadu_si256((const __m256i *)(ref2 + 32));
        ref3_reg     = _mm256_loadu_si256((const __m256i *)ref3);
        ref3next_reg = _mm256_loadu_si256((const __m256i *)(ref3 + 32));
        // sum of the absolute differences between every ref-i to src
        ref0_reg     = _mm256_sad_epu8(ref0_reg, src_reg);
        ref1_reg     = _mm256_sad_epu8(ref1_reg, src_reg);
        ref2_reg     = _mm256_sad_epu8(ref2_reg, src_reg);
        ref3_reg     = _mm256_sad_epu8(ref3_reg, src_reg);
        ref0next_reg = _mm256_sad_epu8(ref0next_reg, srcnext_reg);
        ref1next_reg = _mm256_sad_epu8(ref1next_reg, srcnext_reg);
        ref2next_reg = _mm256_sad_epu8(ref2next_reg, srcnext_reg);
        ref3next_reg = _mm256_sad_epu8(ref3next_reg, srcnext_reg);

        // sum every ref-i
        sum_ref0 = _mm256_add_epi32(sum_ref0, ref0_reg);
        sum_ref1 = _mm256_add_epi32(sum_ref1, ref1_reg);
        sum_ref2 = _mm256_add_epi32(sum_ref2, ref2_reg);
        sum_ref3 = _mm256_add_epi32(sum_ref3, ref3_reg);
        sum_ref0 = _mm256_add_epi32(sum_ref0, ref0next_reg);
        sum_ref1 = _mm256_add_epi32(sum_ref1, ref1next_reg);
        sum_ref2 = _mm256_add_epi32(sum_ref2, ref2next_reg);
        sum_ref3 = _mm256_add_epi32(sum_ref3, ref3next_reg);
        src += src_stride;
        ref0 += ref_stride;
        ref1 += ref_stride;
        ref2 += ref_stride;
        ref3 += ref_stride;
    }
    {
        __m128i sum;

        // in sum_ref-i the result is saved in the first 4 bytes
        // the other 4 bytes are zeroed.
        // sum_ref1 and sum_ref3 are shifted left by 4 bytes
        sum_ref1 = _mm256_slli_si256(sum_ref1, 4);
        sum_ref3 = _mm256_slli_si256(sum_ref3, 4);

        // merge sum_ref0 and sum_ref1 also sum_ref2 and sum_ref3
        sum_ref0 = _mm256_or_si256(sum_ref0, sum_ref1);
        sum_ref2 = _mm256_or_si256(sum_ref2, sum_ref3);

        // merge every 64 bit from each sum_ref-i
        sum_mlow  = _mm256_unpacklo_epi64(sum_ref0, sum_ref2);
        sum_mhigh = _mm256_unpackhi_epi64(sum_ref0, sum_ref2);

        // add the low 64 bit to the high 64 bit
        sum_mlow = _mm256_add_epi32(sum_mlow, sum_mhigh);

        // add the low 128 bit to the high 128 bit
        sum =
            _mm_add_epi32(_mm256_castsi256_si128(sum_mlow), _mm256_extractf128_si256(sum_mlow, 1));

        _mm_storeu_si128((__m128i *)(res), sum);
    }
    _mm256_zeroupper();
}

void eb_aom_sad64x64x4d_avx2(const uint8_t *src, int src_stride, const uint8_t *const ref[4],
                             int ref_stride, uint32_t res[4]) {
    sad64x64x4d_avx2(src, src_stride, ref, ref_stride, res);
}

void eb_aom_sad32x64x4d_avx2(const uint8_t *src, int src_stride, const uint8_t *const ref[4],
                             int ref_stride, uint32_t res[4]) {
    const uint8_t *rf[4];
    uint32_t       sum0[4];
    uint32_t       sum1[4];

    rf[0] = ref[0];
    rf[1] = ref[1];
    rf[2] = ref[2];
    rf[3] = ref[3];
    sad32x32x4d_avx2(src, src_stride, rf, ref_stride, sum0);
    src += src_stride << 5;
    rf[0] += ref_stride << 5;
    rf[1] += ref_stride << 5;
    rf[2] += ref_stride << 5;
    rf[3] += ref_stride << 5;
    sad32x32x4d_avx2(src, src_stride, rf, ref_stride, sum1);
    res[0] = sum0[0] + sum1[0];
    res[1] = sum0[1] + sum1[1];
    res[2] = sum0[2] + sum1[2];
    res[3] = sum0[3] + sum1[3];
}

void eb_aom_sad64x32x4d_avx2(const uint8_t *src, int src_stride, const uint8_t *const ref[4],
                             int ref_stride, uint32_t res[4]) {
    const uint8_t *rf[4];
    uint32_t       sum0[4];
    uint32_t       sum1[4];
    unsigned int   half_width = 32;

    rf[0] = ref[0];
    rf[1] = ref[1];
    rf[2] = ref[2];
    rf[3] = ref[3];
    sad32x32x4d_avx2(src, src_stride, rf, ref_stride, sum0);
    src += half_width;
    rf[0] += half_width;
    rf[1] += half_width;
    rf[2] += half_width;
    rf[3] += half_width;
    sad32x32x4d_avx2(src, src_stride, rf, ref_stride, sum1);
    res[0] = sum0[0] + sum1[0];
    res[1] = sum0[1] + sum1[1];
    res[2] = sum0[2] + sum1[2];
    res[3] = sum0[3] + sum1[3];
}

static INLINE unsigned int sad32x32(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr,
                                    int ref_stride) {
    __m256i s1, s2, r1, r2;
    __m256i sum = _mm256_setzero_si256();
    __m128i sum_i128;
    int     i;

    for (i = 0; i < 16; ++i) {
        r1  = _mm256_loadu_si256((__m256i const *)ref_ptr);
        r2  = _mm256_loadu_si256((__m256i const *)(ref_ptr + ref_stride));
        s1  = _mm256_sad_epu8(r1, _mm256_loadu_si256((__m256i const *)src_ptr));
        s2  = _mm256_sad_epu8(r2, _mm256_loadu_si256((__m256i const *)(src_ptr + src_stride)));
        sum = _mm256_add_epi32(sum, _mm256_add_epi32(s1, s2));
        ref_ptr += ref_stride << 1;
        src_ptr += src_stride << 1;
    }

    sum      = _mm256_add_epi32(sum, _mm256_srli_si256(sum, 8));
    sum_i128 = _mm_add_epi32(_mm256_extracti128_si256(sum, 1), _mm256_castsi256_si128(sum));
    return _mm_cvtsi128_si32(sum_i128);
}

SIMD_INLINE unsigned int sad64x32(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr,
                                  int ref_stride) {
    unsigned int half_width = 32;
    uint32_t     sum        = sad32x32(src_ptr, src_stride, ref_ptr, ref_stride);
    src_ptr += half_width;
    ref_ptr += half_width;
    sum += sad32x32(src_ptr, src_stride, ref_ptr, ref_stride);
    return sum;
}

SIMD_INLINE unsigned int sad64x64(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr,
                                  int ref_stride) {
    uint32_t sum = sad64x32(src_ptr, src_stride, ref_ptr, ref_stride);
    src_ptr += src_stride << 5;
    ref_ptr += ref_stride << 5;
    sum += sad64x32(src_ptr, src_stride, ref_ptr, ref_stride);
    return sum;
}

unsigned int eb_aom_sad64x128_avx2(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr,
                                   int ref_stride) {
    uint32_t sum = sad64x64(src_ptr, src_stride, ref_ptr, ref_stride);
    src_ptr += src_stride << 6;
    ref_ptr += ref_stride << 6;
    sum += sad64x64(src_ptr, src_stride, ref_ptr, ref_stride);
    return sum;
}

SIMD_INLINE void sad64x64x4d(const uint8_t *src, int src_stride, const uint8_t *const ref[4],
                             int ref_stride, __m128i *res) {
    uint32_t sum[4];
    sad64x64x4d_avx2(src, src_stride, ref, ref_stride, sum);
    *res = _mm_loadu_si128((const __m128i *)sum);
}

void eb_aom_sad64x128x4d_avx2(const uint8_t *src, int src_stride, const uint8_t *const ref[4],
                              int ref_stride, uint32_t res[4]) {
    __m128i        sum0, sum1;
    const uint8_t *rf[4];

    rf[0] = ref[0];
    rf[1] = ref[1];
    rf[2] = ref[2];
    rf[3] = ref[3];
    sad64x64x4d(src, src_stride, rf, ref_stride, &sum0);
    src += src_stride << 6;
    rf[0] += ref_stride << 6;
    rf[1] += ref_stride << 6;
    rf[2] += ref_stride << 6;
    rf[3] += ref_stride << 6;
    sad64x64x4d(src, src_stride, rf, ref_stride, &sum1);
    sum0 = _mm_add_epi32(sum0, sum1);
    _mm_storeu_si128((__m128i *)res, sum0);
}

SIMD_INLINE void sad128x64x4d_avx2(const uint8_t *src, int src_stride, const uint8_t *const ref[4],
                                   int ref_stride, uint32_t res[4]) {
    __m128i        sum0, sum1;
    unsigned int   half_width = 64;
    const uint8_t *rf[4];

    rf[0] = ref[0];
    rf[1] = ref[1];
    rf[2] = ref[2];
    rf[3] = ref[3];
    sad64x64x4d(src, src_stride, rf, ref_stride, &sum0);
    src += half_width;
    rf[0] += half_width;
    rf[1] += half_width;
    rf[2] += half_width;
    rf[3] += half_width;
    sad64x64x4d(src, src_stride, rf, ref_stride, &sum1);
    sum0 = _mm_add_epi32(sum0, sum1);
    _mm_storeu_si128((__m128i *)res, sum0);
}

void eb_aom_sad128x64x4d_avx2(const uint8_t *src, int src_stride, const uint8_t *const ref[4],
                              int ref_stride, uint32_t res[4]) {
    sad128x64x4d_avx2(src, src_stride, ref, ref_stride, res);
}

void eb_aom_sad128x128x4d_avx2(const uint8_t *src, int src_stride, const uint8_t *const ref[4],
                               int ref_stride, uint32_t res[4]) {
    const uint8_t *rf[4];
    uint32_t       sum0[4];
    uint32_t       sum1[4];

    rf[0] = ref[0];
    rf[1] = ref[1];
    rf[2] = ref[2];
    rf[3] = ref[3];
    sad128x64x4d_avx2(src, src_stride, rf, ref_stride, sum0);
    src += src_stride << 6;
    rf[0] += ref_stride << 6;
    rf[1] += ref_stride << 6;
    rf[2] += ref_stride << 6;
    rf[3] += ref_stride << 6;
    sad128x64x4d_avx2(src, src_stride, rf, ref_stride, sum1);
    res[0] = sum0[0] + sum1[0];
    res[1] = sum0[1] + sum1[1];
    res[2] = sum0[2] + sum1[2];
    res[3] = sum0[3] + sum1[3];
}

void ext_all_sad_calculation_8x8_16x16_avx2(uint8_t *src, uint32_t src_stride, uint8_t *ref,
                                            uint32_t ref_stride, uint32_t mv,
                                            uint32_t *p_best_sad_8x8, uint32_t *p_best_sad_16x16,
                                            uint32_t *p_best_mv8x8, uint32_t *p_best_mv16x16,
                                            uint32_t p_eight_sad16x16[16][8],
                                            uint32_t p_eight_sad8x8[64][8]) {
    static const char offsets[16] = {0, 1, 4, 5, 2, 3, 6, 7, 8, 9, 12, 13, 10, 11, 14, 15};

    //---- 16x16 : 0, 1, 4, 5, 2, 3, 6, 7, 8, 9, 12, 13, 10, 11, 14, 15
    for (int y = 0; y < 4; y++) {
        for (int x = 0; x < 4; x++) {
            const uint32_t start_16x16_pos = offsets[4 * y + x];
            const uint32_t start_8x8_pos   = 4 * start_16x16_pos;
            const uint8_t *s               = src + 16 * y * src_stride + 16 * x;
            const uint8_t *r               = ref + 16 * y * ref_stride + 16 * x;
            __m256i        sad02           = _mm256_setzero_si256();
            __m256i        sad13           = _mm256_setzero_si256();

            for (int i = 0; i < 4; i++) {
                const __m128i src01 = _mm_loadu_si128((__m128i *)(s + 0 * src_stride));
                const __m128i src23 = _mm_loadu_si128((__m128i *)(s + 8 * src_stride));
                const __m128i ref0  = _mm_loadu_si128((__m128i *)(r + 0 * ref_stride + 0));
                const __m128i ref1  = _mm_loadu_si128((__m128i *)(r + 0 * ref_stride + 8));
                const __m128i ref2  = _mm_loadu_si128((__m128i *)(r + 8 * ref_stride + 0));
                const __m128i ref3  = _mm_loadu_si128((__m128i *)(r + 8 * ref_stride + 8));
                const __m256i src0123 =
                    _mm256_insertf128_si256(_mm256_castsi128_si256(src01), src23, 1);
                const __m256i ref02 =
                    _mm256_insertf128_si256(_mm256_castsi128_si256(ref0), ref2, 1);
                const __m256i ref13 =
                    _mm256_insertf128_si256(_mm256_castsi128_si256(ref1), ref3, 1);
                sad02 = _mm256_adds_epu16(sad02, _mm256_mpsadbw_epu8(ref02, src0123, 0)); // 000 000
                sad02 =
                    _mm256_adds_epu16(sad02, _mm256_mpsadbw_epu8(ref02, src0123, 45)); // 101 101
                sad13 =
                    _mm256_adds_epu16(sad13, _mm256_mpsadbw_epu8(ref13, src0123, 18)); // 010 010
                sad13 =
                    _mm256_adds_epu16(sad13, _mm256_mpsadbw_epu8(ref13, src0123, 63)); // 111 111
                s += 2 * src_stride;
                r += 2 * ref_stride;
            }

            sad02 = _mm256_slli_epi16(sad02, 1);
            sad13 = _mm256_slli_epi16(sad13, 1);

            const __m256i sad0202 = _mm256_permute4x64_epi64(sad02, 0xD8);
            const __m256i sad1313 = _mm256_permute4x64_epi64(sad13, 0xD8);
            const __m256i sad00   = _mm256_unpacklo_epi16(sad0202, _mm256_setzero_si256());
            const __m256i sad11   = _mm256_unpacklo_epi16(sad1313, _mm256_setzero_si256());
            const __m256i sad22   = _mm256_unpackhi_epi16(sad0202, _mm256_setzero_si256());
            const __m256i sad33   = _mm256_unpackhi_epi16(sad1313, _mm256_setzero_si256());

            _mm256_storeu_si256((__m256i *)(p_eight_sad8x8[0 + start_8x8_pos]), sad00);
            _mm256_storeu_si256((__m256i *)(p_eight_sad8x8[1 + start_8x8_pos]), sad11);
            _mm256_storeu_si256((__m256i *)(p_eight_sad8x8[2 + start_8x8_pos]), sad22);
            _mm256_storeu_si256((__m256i *)(p_eight_sad8x8[3 + start_8x8_pos]), sad33);

            const __m128i sad0 = _mm256_castsi256_si128(sad02);
            const __m128i sad1 = _mm256_castsi256_si128(sad13);
            const __m128i sad2 = _mm256_extracti128_si256(sad02, 1);
            const __m128i sad3 = _mm256_extracti128_si256(sad13, 1);

            const __m128i minpos0 = _mm_minpos_epu16(sad0);
            const __m128i minpos1 = _mm_minpos_epu16(sad1);
            const __m128i minpos2 = _mm_minpos_epu16(sad2);
            const __m128i minpos3 = _mm_minpos_epu16(sad3);

            const __m128i minpos01   = _mm_unpacklo_epi16(minpos0, minpos1);
            const __m128i minpos23   = _mm_unpacklo_epi16(minpos2, minpos3);
            const __m128i minpos0123 = _mm_unpacklo_epi32(minpos01, minpos23);
            const __m128i sad8x8     = _mm_unpacklo_epi16(minpos0123, _mm_setzero_si128());
            const __m128i pos0123    = _mm_unpackhi_epi16(minpos0123, _mm_setzero_si128());
            const __m128i pos8x8     = _mm_slli_epi32(pos0123, 2);

            __m128i best_sad8x8 = _mm_loadu_si128((__m128i *)(p_best_sad_8x8 + start_8x8_pos));
            const __m128i mask  = _mm_cmplt_epi32(sad8x8, best_sad8x8);
            best_sad8x8         = _mm_min_epi32(best_sad8x8, sad8x8);
            _mm_storeu_si128((__m128i *)(p_best_sad_8x8 + start_8x8_pos), best_sad8x8);

            __m128i       best_mv8x8 = _mm_loadu_si128((__m128i *)(p_best_mv8x8 + start_8x8_pos));
            const __m128i mvs        = _mm_set1_epi32(mv);
            const __m128i mv8x8      = _mm_add_epi16(mvs, pos8x8);
            best_mv8x8               = _mm_blendv_epi8(best_mv8x8, mv8x8, mask);
            _mm_storeu_si128((__m128i *)(p_best_mv8x8 + start_8x8_pos), best_mv8x8);

            const __m128i sum01       = _mm_add_epi16(sad0, sad1);
            const __m128i sum23       = _mm_add_epi16(sad2, sad3);
            const __m128i sad16x16_16 = _mm_add_epi16(sum01, sum23);
            const __m256i sad16x16_32 = _mm256_cvtepu16_epi32(sad16x16_16);
            _mm256_storeu_si256((__m256i *)(p_eight_sad16x16[start_16x16_pos]), sad16x16_32);

            const __m128i  minpos16x16 = _mm_minpos_epu16(sad16x16_16);
            const uint32_t min16x16    = _mm_extract_epi16(minpos16x16, 0);

            if (min16x16 < p_best_sad_16x16[start_16x16_pos]) {
                p_best_sad_16x16[start_16x16_pos] = min16x16;

                const __m128i pos               = _mm_srli_si128(minpos16x16, 2);
                const __m128i pos16x16          = _mm_slli_epi32(pos, 2);
                const __m128i mv16x16           = _mm_add_epi16(mvs, pos16x16);
                p_best_mv16x16[start_16x16_pos] = _mm_extract_epi32(mv16x16, 0);
            }
        }
    }
}

#define avx2_find_min_pos_init() \
    const __m256i idx_min_pos_prv = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7)

/* Horizontally search the minimum value in vector.
   Store the minimum index in best_id.
   To fast calculations, reduce available max value to 27 bits.
   First call avx2_find_min_pos_init in function
   */
#define avx2_find_min_pos(in, best_id)                                              \
    {                                                                               \
        __m256i x        = _mm256_slli_epi32(in, 3); /* x = x<<3 | idx */           \
        x                = _mm256_add_epi32(x, idx_min_pos_prv); /*x = [12345678]*/ \
        const __m256i y  = _mm256_permute2f128_si256(x, x, 1); /*y = [56781234]*/   \
        const __m256i m1 = _mm256_min_epu32(x, y); /*m1 = [12341234]*/              \
        const __m256i m2 = _mm256_permute4x64_epi64(m1, 5); /*m2 = [34123412]*/     \
        const __m256i m3 = _mm256_min_epu32(m1, m2); /*m3 = [12121212]*/            \
        const __m256i m4 = _mm256_shuffle_epi32(m3, 5); /*m4 = [21212121]*/         \
        const __m256i m  = _mm256_min_epu32(m3, m4); /*m5 = [11111111]*/            \
        /* (x<<3 | idx) & (0b000111) = idx */                                       \
        best_id = _mm256_extract_epi16(m, 0) & 0x07;                                \
    }

#if !SHUT_ME_NSQ_SEARCH
void ext_eigth_sad_calculation_nsq_avx2(
    uint32_t p_sad8x8[64][8], uint32_t p_sad16x16[16][8], uint32_t p_sad32x32[4][8],
    uint32_t *p_best_sad_64x32, uint32_t *p_best_mv64x32, uint32_t *p_best_sad_32x16,
    uint32_t *p_best_mv32x16, uint32_t *p_best_sad_16x8, uint32_t *p_best_mv16x8,
    uint32_t *p_best_sad_32x64, uint32_t *p_best_mv32x64, uint32_t *p_best_sad_16x32,
    uint32_t *p_best_mv16x32, uint32_t *p_best_sad_8x16, uint32_t *p_best_mv8x16,
    uint32_t *p_best_sad_32x8, uint32_t *p_best_mv32x8, uint32_t *p_best_sad_8x32,
    uint32_t *p_best_mv8x32, uint32_t *p_best_sad_64x16, uint32_t *p_best_mv64x16,
    uint32_t *p_best_sad_16x64, uint32_t *p_best_mv16x64, uint32_t mv) {
    avx2_find_min_pos_init();
    uint8_t search_index;
    DECLARE_ALIGNED(32, uint32_t, sad[8]);
    DECLARE_ALIGNED(32, uint32_t, sad_16x8[32][8]);
    DECLARE_ALIGNED(32, uint32_t, sad_8x16[32][8]);
    DECLARE_ALIGNED(32, uint32_t, sad_32x16[8][8]);
    DECLARE_ALIGNED(32, uint32_t, sad_16x32[8][8]);
    DECLARE_ALIGNED(32, uint32_t, computed_idx[8]);

    __m256i search_idx_avx2 = _mm256_setr_epi32(0, 4, 8, 12, 16, 20, 24, 28);
    __m256i mv_avx2         = _mm256_set1_epi32(mv);
    __m256i new_mv_avx2     = _mm256_add_epi32(search_idx_avx2, mv_avx2);
    new_mv_avx2             = _mm256_and_si256(new_mv_avx2, _mm256_set1_epi32(0xffff));
    *(__m256i *)computed_idx =
        _mm256_or_si256(new_mv_avx2, _mm256_and_si256(mv_avx2, _mm256_set1_epi32(0xffff0000)));

#define ADD_VECS(dst, a, b)                                                    \
    *(__m256i *)dst = _mm256_add_epi32(_mm256_loadu_si256((__m256i const *)a), \
                                       _mm256_loadu_si256((__m256i const *)b));

#define SEARCH_BEST(serch, best, best_mv)               \
    avx2_find_min_pos(*(__m256i *)serch, search_index); \
    if (serch[search_index] < best) {                   \
        best    = serch[search_index];                  \
        best_mv = computed_idx[search_index];           \
    }

#define SEARCH_CALC_BEST(a, b, best, best_mv) \
    ADD_VECS(sad, a, b);                      \
    SEARCH_BEST(sad, best, best_mv)

    // 32x16
    ADD_VECS(sad_32x16[0], p_sad16x16[0], p_sad16x16[1]);
    ADD_VECS(sad_32x16[1], p_sad16x16[2], p_sad16x16[3]);
    ADD_VECS(sad_32x16[2], p_sad16x16[4], p_sad16x16[5]);
    ADD_VECS(sad_32x16[3], p_sad16x16[6], p_sad16x16[7]);
    ADD_VECS(sad_32x16[4], p_sad16x16[8], p_sad16x16[9]);
    ADD_VECS(sad_32x16[5], p_sad16x16[10], p_sad16x16[11]);
    ADD_VECS(sad_32x16[6], p_sad16x16[12], p_sad16x16[13]);
    ADD_VECS(sad_32x16[7], p_sad16x16[14], p_sad16x16[15]);

    // 16x8
    ADD_VECS(sad_16x8[0], p_sad8x8[0], p_sad8x8[1]);
    ADD_VECS(sad_16x8[1], p_sad8x8[2], p_sad8x8[3]);
    ADD_VECS(sad_16x8[2], p_sad8x8[4], p_sad8x8[5]);
    ADD_VECS(sad_16x8[3], p_sad8x8[6], p_sad8x8[7]);
    ADD_VECS(sad_16x8[4], p_sad8x8[8], p_sad8x8[9]);
    ADD_VECS(sad_16x8[5], p_sad8x8[10], p_sad8x8[11]);
    ADD_VECS(sad_16x8[6], p_sad8x8[12], p_sad8x8[13]);
    ADD_VECS(sad_16x8[7], p_sad8x8[14], p_sad8x8[15]);
    ADD_VECS(sad_16x8[8], p_sad8x8[16], p_sad8x8[17]);
    ADD_VECS(sad_16x8[9], p_sad8x8[18], p_sad8x8[19]);
    ADD_VECS(sad_16x8[10], p_sad8x8[20], p_sad8x8[21]);
    ADD_VECS(sad_16x8[11], p_sad8x8[22], p_sad8x8[23]);
    ADD_VECS(sad_16x8[12], p_sad8x8[24], p_sad8x8[25]);
    ADD_VECS(sad_16x8[13], p_sad8x8[26], p_sad8x8[27]);
    ADD_VECS(sad_16x8[14], p_sad8x8[28], p_sad8x8[29]);
    ADD_VECS(sad_16x8[15], p_sad8x8[30], p_sad8x8[31]);
    ADD_VECS(sad_16x8[16], p_sad8x8[32], p_sad8x8[33]);
    ADD_VECS(sad_16x8[17], p_sad8x8[34], p_sad8x8[35]);
    ADD_VECS(sad_16x8[18], p_sad8x8[36], p_sad8x8[37]);
    ADD_VECS(sad_16x8[19], p_sad8x8[38], p_sad8x8[39]);
    ADD_VECS(sad_16x8[20], p_sad8x8[40], p_sad8x8[41]);
    ADD_VECS(sad_16x8[21], p_sad8x8[42], p_sad8x8[43]);
    ADD_VECS(sad_16x8[22], p_sad8x8[44], p_sad8x8[45]);
    ADD_VECS(sad_16x8[23], p_sad8x8[46], p_sad8x8[47]);
    ADD_VECS(sad_16x8[24], p_sad8x8[48], p_sad8x8[49]);
    ADD_VECS(sad_16x8[25], p_sad8x8[50], p_sad8x8[51]);
    ADD_VECS(sad_16x8[26], p_sad8x8[52], p_sad8x8[53]);
    ADD_VECS(sad_16x8[27], p_sad8x8[54], p_sad8x8[55]);
    ADD_VECS(sad_16x8[28], p_sad8x8[56], p_sad8x8[57]);
    ADD_VECS(sad_16x8[29], p_sad8x8[58], p_sad8x8[59]);
    ADD_VECS(sad_16x8[30], p_sad8x8[60], p_sad8x8[61]);
    ADD_VECS(sad_16x8[31], p_sad8x8[62], p_sad8x8[63]);
    // 16x32
    ADD_VECS(sad_16x32[0], p_sad16x16[0], p_sad16x16[2]);
    ADD_VECS(sad_16x32[1], p_sad16x16[1], p_sad16x16[3]);
    ADD_VECS(sad_16x32[2], p_sad16x16[4], p_sad16x16[6]);
    ADD_VECS(sad_16x32[3], p_sad16x16[5], p_sad16x16[7]);
    ADD_VECS(sad_16x32[4], p_sad16x16[8], p_sad16x16[10]);
    ADD_VECS(sad_16x32[5], p_sad16x16[9], p_sad16x16[11]);
    ADD_VECS(sad_16x32[6], p_sad16x16[12], p_sad16x16[14]);
    ADD_VECS(sad_16x32[7], p_sad16x16[13], p_sad16x16[15]);
    // 8x16
    ADD_VECS(sad_8x16[0], p_sad8x8[0], p_sad8x8[2]);
    ADD_VECS(sad_8x16[1], p_sad8x8[1], p_sad8x8[3]);
    ADD_VECS(sad_8x16[2], p_sad8x8[4], p_sad8x8[6]);
    ADD_VECS(sad_8x16[3], p_sad8x8[5], p_sad8x8[7]);
    ADD_VECS(sad_8x16[4], p_sad8x8[8], p_sad8x8[10]);
    ADD_VECS(sad_8x16[5], p_sad8x8[9], p_sad8x8[11]);
    ADD_VECS(sad_8x16[6], p_sad8x8[12], p_sad8x8[14]);
    ADD_VECS(sad_8x16[7], p_sad8x8[13], p_sad8x8[15]);
    ADD_VECS(sad_8x16[8], p_sad8x8[16], p_sad8x8[18]);
    ADD_VECS(sad_8x16[9], p_sad8x8[17], p_sad8x8[19]);
    ADD_VECS(sad_8x16[10], p_sad8x8[20], p_sad8x8[22]);
    ADD_VECS(sad_8x16[11], p_sad8x8[21], p_sad8x8[23]);
    ADD_VECS(sad_8x16[12], p_sad8x8[24], p_sad8x8[26]);
    ADD_VECS(sad_8x16[13], p_sad8x8[25], p_sad8x8[27]);
    ADD_VECS(sad_8x16[14], p_sad8x8[28], p_sad8x8[30]);
    ADD_VECS(sad_8x16[15], p_sad8x8[29], p_sad8x8[31]);
    ADD_VECS(sad_8x16[16], p_sad8x8[32], p_sad8x8[34]);
    ADD_VECS(sad_8x16[17], p_sad8x8[33], p_sad8x8[35]);
    ADD_VECS(sad_8x16[18], p_sad8x8[36], p_sad8x8[38]);
    ADD_VECS(sad_8x16[19], p_sad8x8[37], p_sad8x8[39]);
    ADD_VECS(sad_8x16[20], p_sad8x8[40], p_sad8x8[42]);
    ADD_VECS(sad_8x16[21], p_sad8x8[41], p_sad8x8[43]);
    ADD_VECS(sad_8x16[22], p_sad8x8[44], p_sad8x8[46]);
    ADD_VECS(sad_8x16[23], p_sad8x8[45], p_sad8x8[47]);
    ADD_VECS(sad_8x16[24], p_sad8x8[48], p_sad8x8[50]);
    ADD_VECS(sad_8x16[25], p_sad8x8[49], p_sad8x8[51]);
    ADD_VECS(sad_8x16[26], p_sad8x8[52], p_sad8x8[54]);
    ADD_VECS(sad_8x16[27], p_sad8x8[53], p_sad8x8[55]);
    ADD_VECS(sad_8x16[28], p_sad8x8[56], p_sad8x8[58]);
    ADD_VECS(sad_8x16[29], p_sad8x8[57], p_sad8x8[59]);
    ADD_VECS(sad_8x16[30], p_sad8x8[60], p_sad8x8[62]);
    ADD_VECS(sad_8x16[31], p_sad8x8[61], p_sad8x8[63]);

    // 32x16
    SEARCH_BEST(sad_32x16[0], p_best_sad_32x16[0], p_best_mv32x16[0]);
    SEARCH_BEST(sad_32x16[1], p_best_sad_32x16[1], p_best_mv32x16[1]);
    SEARCH_BEST(sad_32x16[2], p_best_sad_32x16[2], p_best_mv32x16[2]);
    SEARCH_BEST(sad_32x16[3], p_best_sad_32x16[3], p_best_mv32x16[3]);
    SEARCH_BEST(sad_32x16[4], p_best_sad_32x16[4], p_best_mv32x16[4]);
    SEARCH_BEST(sad_32x16[5], p_best_sad_32x16[5], p_best_mv32x16[5]);
    SEARCH_BEST(sad_32x16[6], p_best_sad_32x16[6], p_best_mv32x16[6]);
    SEARCH_BEST(sad_32x16[7], p_best_sad_32x16[7], p_best_mv32x16[7]);
    // 16x8
    SEARCH_BEST(sad_16x8[0], p_best_sad_16x8[0], p_best_mv16x8[0]);
    SEARCH_BEST(sad_16x8[1], p_best_sad_16x8[1], p_best_mv16x8[1]);
    SEARCH_BEST(sad_16x8[2], p_best_sad_16x8[2], p_best_mv16x8[2]);
    SEARCH_BEST(sad_16x8[3], p_best_sad_16x8[3], p_best_mv16x8[3]);
    SEARCH_BEST(sad_16x8[4], p_best_sad_16x8[4], p_best_mv16x8[4]);
    SEARCH_BEST(sad_16x8[5], p_best_sad_16x8[5], p_best_mv16x8[5]);
    SEARCH_BEST(sad_16x8[6], p_best_sad_16x8[6], p_best_mv16x8[6]);
    SEARCH_BEST(sad_16x8[7], p_best_sad_16x8[7], p_best_mv16x8[7]);
    SEARCH_BEST(sad_16x8[8], p_best_sad_16x8[8], p_best_mv16x8[8]);
    SEARCH_BEST(sad_16x8[9], p_best_sad_16x8[9], p_best_mv16x8[9]);
    SEARCH_BEST(sad_16x8[10], p_best_sad_16x8[10], p_best_mv16x8[10]);
    SEARCH_BEST(sad_16x8[11], p_best_sad_16x8[11], p_best_mv16x8[11]);
    SEARCH_BEST(sad_16x8[12], p_best_sad_16x8[12], p_best_mv16x8[12]);
    SEARCH_BEST(sad_16x8[13], p_best_sad_16x8[13], p_best_mv16x8[13]);
    SEARCH_BEST(sad_16x8[14], p_best_sad_16x8[14], p_best_mv16x8[14]);
    SEARCH_BEST(sad_16x8[15], p_best_sad_16x8[15], p_best_mv16x8[15]);
    SEARCH_BEST(sad_16x8[16], p_best_sad_16x8[16], p_best_mv16x8[16]);
    SEARCH_BEST(sad_16x8[17], p_best_sad_16x8[17], p_best_mv16x8[17]);
    SEARCH_BEST(sad_16x8[18], p_best_sad_16x8[18], p_best_mv16x8[18]);
    SEARCH_BEST(sad_16x8[19], p_best_sad_16x8[19], p_best_mv16x8[19]);
    SEARCH_BEST(sad_16x8[20], p_best_sad_16x8[20], p_best_mv16x8[20]);
    SEARCH_BEST(sad_16x8[21], p_best_sad_16x8[21], p_best_mv16x8[21]);
    SEARCH_BEST(sad_16x8[22], p_best_sad_16x8[22], p_best_mv16x8[22]);
    SEARCH_BEST(sad_16x8[23], p_best_sad_16x8[23], p_best_mv16x8[23]);
    SEARCH_BEST(sad_16x8[24], p_best_sad_16x8[24], p_best_mv16x8[24]);
    SEARCH_BEST(sad_16x8[25], p_best_sad_16x8[25], p_best_mv16x8[25]);
    SEARCH_BEST(sad_16x8[26], p_best_sad_16x8[26], p_best_mv16x8[26]);
    SEARCH_BEST(sad_16x8[27], p_best_sad_16x8[27], p_best_mv16x8[27]);
    SEARCH_BEST(sad_16x8[28], p_best_sad_16x8[28], p_best_mv16x8[28]);
    SEARCH_BEST(sad_16x8[29], p_best_sad_16x8[29], p_best_mv16x8[29]);
    SEARCH_BEST(sad_16x8[30], p_best_sad_16x8[30], p_best_mv16x8[30]);
    SEARCH_BEST(sad_16x8[31], p_best_sad_16x8[31], p_best_mv16x8[31]);
    // 16x32
    SEARCH_BEST(sad_16x32[0], p_best_sad_16x32[0], p_best_mv16x32[0]);
    SEARCH_BEST(sad_16x32[1], p_best_sad_16x32[1], p_best_mv16x32[1]);
    SEARCH_BEST(sad_16x32[2], p_best_sad_16x32[2], p_best_mv16x32[2]);
    SEARCH_BEST(sad_16x32[3], p_best_sad_16x32[3], p_best_mv16x32[3]);
    SEARCH_BEST(sad_16x32[4], p_best_sad_16x32[4], p_best_mv16x32[4]);
    SEARCH_BEST(sad_16x32[5], p_best_sad_16x32[5], p_best_mv16x32[5]);
    SEARCH_BEST(sad_16x32[6], p_best_sad_16x32[6], p_best_mv16x32[6]);
    SEARCH_BEST(sad_16x32[7], p_best_sad_16x32[7], p_best_mv16x32[7]);
    // 8x16
    SEARCH_BEST(sad_8x16[0], p_best_sad_8x16[0], p_best_mv8x16[0]);
    SEARCH_BEST(sad_8x16[1], p_best_sad_8x16[1], p_best_mv8x16[1]);
    SEARCH_BEST(sad_8x16[2], p_best_sad_8x16[2], p_best_mv8x16[2]);
    SEARCH_BEST(sad_8x16[3], p_best_sad_8x16[3], p_best_mv8x16[3]);
    SEARCH_BEST(sad_8x16[4], p_best_sad_8x16[4], p_best_mv8x16[4]);
    SEARCH_BEST(sad_8x16[5], p_best_sad_8x16[5], p_best_mv8x16[5]);
    SEARCH_BEST(sad_8x16[6], p_best_sad_8x16[6], p_best_mv8x16[6]);
    SEARCH_BEST(sad_8x16[7], p_best_sad_8x16[7], p_best_mv8x16[7]);
    SEARCH_BEST(sad_8x16[8], p_best_sad_8x16[8], p_best_mv8x16[8]);
    SEARCH_BEST(sad_8x16[9], p_best_sad_8x16[9], p_best_mv8x16[9]);
    SEARCH_BEST(sad_8x16[10], p_best_sad_8x16[10], p_best_mv8x16[10]);
    SEARCH_BEST(sad_8x16[11], p_best_sad_8x16[11], p_best_mv8x16[11]);
    SEARCH_BEST(sad_8x16[12], p_best_sad_8x16[12], p_best_mv8x16[12]);
    SEARCH_BEST(sad_8x16[13], p_best_sad_8x16[13], p_best_mv8x16[13]);
    SEARCH_BEST(sad_8x16[14], p_best_sad_8x16[14], p_best_mv8x16[14]);
    SEARCH_BEST(sad_8x16[15], p_best_sad_8x16[15], p_best_mv8x16[15]);
    SEARCH_BEST(sad_8x16[16], p_best_sad_8x16[16], p_best_mv8x16[16]);
    SEARCH_BEST(sad_8x16[17], p_best_sad_8x16[17], p_best_mv8x16[17]);
    SEARCH_BEST(sad_8x16[18], p_best_sad_8x16[18], p_best_mv8x16[18]);
    SEARCH_BEST(sad_8x16[19], p_best_sad_8x16[19], p_best_mv8x16[19]);
    SEARCH_BEST(sad_8x16[20], p_best_sad_8x16[20], p_best_mv8x16[20]);
    SEARCH_BEST(sad_8x16[21], p_best_sad_8x16[21], p_best_mv8x16[21]);
    SEARCH_BEST(sad_8x16[22], p_best_sad_8x16[22], p_best_mv8x16[22]);
    SEARCH_BEST(sad_8x16[23], p_best_sad_8x16[23], p_best_mv8x16[23]);
    SEARCH_BEST(sad_8x16[24], p_best_sad_8x16[24], p_best_mv8x16[24]);
    SEARCH_BEST(sad_8x16[25], p_best_sad_8x16[25], p_best_mv8x16[25]);
    SEARCH_BEST(sad_8x16[26], p_best_sad_8x16[26], p_best_mv8x16[26]);
    SEARCH_BEST(sad_8x16[27], p_best_sad_8x16[27], p_best_mv8x16[27]);
    SEARCH_BEST(sad_8x16[28], p_best_sad_8x16[28], p_best_mv8x16[28]);
    SEARCH_BEST(sad_8x16[29], p_best_sad_8x16[29], p_best_mv8x16[29]);
    SEARCH_BEST(sad_8x16[30], p_best_sad_8x16[30], p_best_mv8x16[30]);
    SEARCH_BEST(sad_8x16[31], p_best_sad_8x16[31], p_best_mv8x16[31]);

    SEARCH_CALC_BEST(p_sad32x32[0], p_sad32x32[1], p_best_sad_64x32[0], p_best_mv64x32[0]);
    SEARCH_CALC_BEST(p_sad32x32[2], p_sad32x32[3], p_best_sad_64x32[1], p_best_mv64x32[1]);
    SEARCH_CALC_BEST(sad_32x16[0], sad_32x16[2], p_best_sad_64x16[0], p_best_mv64x16[0]);
    SEARCH_CALC_BEST(sad_32x16[1], sad_32x16[3], p_best_sad_64x16[1], p_best_mv64x16[1]);
    SEARCH_CALC_BEST(sad_32x16[4], sad_32x16[6], p_best_sad_64x16[2], p_best_mv64x16[2]);
    SEARCH_CALC_BEST(sad_32x16[5], sad_32x16[7], p_best_sad_64x16[3], p_best_mv64x16[3]);
    SEARCH_CALC_BEST(p_sad32x32[0], p_sad32x32[2], p_best_sad_32x64[0], p_best_mv32x64[0]);
    SEARCH_CALC_BEST(p_sad32x32[1], p_sad32x32[3], p_best_sad_32x64[1], p_best_mv32x64[1]);
    SEARCH_CALC_BEST(sad_16x32[0], sad_16x32[4], p_best_sad_16x64[0], p_best_mv16x64[0]);
    SEARCH_CALC_BEST(sad_16x32[1], sad_16x32[5], p_best_sad_16x64[1], p_best_mv16x64[1]);
    SEARCH_CALC_BEST(sad_16x32[2], sad_16x32[6], p_best_sad_16x64[2], p_best_mv16x64[2]);
    SEARCH_CALC_BEST(sad_16x32[3], sad_16x32[7], p_best_sad_16x64[3], p_best_mv16x64[3]);
    SEARCH_CALC_BEST(sad_16x8[0], sad_16x8[2], p_best_sad_32x8[0], p_best_mv32x8[0]);
    SEARCH_CALC_BEST(sad_16x8[1], sad_16x8[3], p_best_sad_32x8[1], p_best_mv32x8[1]);
    SEARCH_CALC_BEST(sad_16x8[4], sad_16x8[6], p_best_sad_32x8[2], p_best_mv32x8[2]);
    SEARCH_CALC_BEST(sad_16x8[5], sad_16x8[7], p_best_sad_32x8[3], p_best_mv32x8[3]);
    SEARCH_CALC_BEST(sad_16x8[8], sad_16x8[10], p_best_sad_32x8[4], p_best_mv32x8[4]);
    SEARCH_CALC_BEST(sad_16x8[9], sad_16x8[11], p_best_sad_32x8[5], p_best_mv32x8[5]);
    SEARCH_CALC_BEST(sad_16x8[12], sad_16x8[14], p_best_sad_32x8[6], p_best_mv32x8[6]);
    SEARCH_CALC_BEST(sad_16x8[13], sad_16x8[15], p_best_sad_32x8[7], p_best_mv32x8[7]);
    SEARCH_CALC_BEST(sad_16x8[16], sad_16x8[18], p_best_sad_32x8[8], p_best_mv32x8[8]);
    SEARCH_CALC_BEST(sad_16x8[17], sad_16x8[19], p_best_sad_32x8[9], p_best_mv32x8[9]);
    SEARCH_CALC_BEST(sad_16x8[20], sad_16x8[22], p_best_sad_32x8[10], p_best_mv32x8[10]);
    SEARCH_CALC_BEST(sad_16x8[21], sad_16x8[23], p_best_sad_32x8[11], p_best_mv32x8[11]);
    SEARCH_CALC_BEST(sad_16x8[24], sad_16x8[26], p_best_sad_32x8[12], p_best_mv32x8[12]);
    SEARCH_CALC_BEST(sad_16x8[25], sad_16x8[27], p_best_sad_32x8[13], p_best_mv32x8[13]);
    SEARCH_CALC_BEST(sad_16x8[28], sad_16x8[30], p_best_sad_32x8[14], p_best_mv32x8[14]);
    SEARCH_CALC_BEST(sad_16x8[29], sad_16x8[31], p_best_sad_32x8[15], p_best_mv32x8[15]);
    SEARCH_CALC_BEST(sad_8x16[0], sad_8x16[4], p_best_sad_8x32[0], p_best_mv8x32[0]);
    SEARCH_CALC_BEST(sad_8x16[1], sad_8x16[5], p_best_sad_8x32[1], p_best_mv8x32[1]);
    SEARCH_CALC_BEST(sad_8x16[2], sad_8x16[6], p_best_sad_8x32[2], p_best_mv8x32[2]);
    SEARCH_CALC_BEST(sad_8x16[3], sad_8x16[7], p_best_sad_8x32[3], p_best_mv8x32[3]);
    SEARCH_CALC_BEST(sad_8x16[8], sad_8x16[12], p_best_sad_8x32[4], p_best_mv8x32[4]);
    SEARCH_CALC_BEST(sad_8x16[9], sad_8x16[13], p_best_sad_8x32[5], p_best_mv8x32[5]);
    SEARCH_CALC_BEST(sad_8x16[10], sad_8x16[14], p_best_sad_8x32[6], p_best_mv8x32[6]);
    SEARCH_CALC_BEST(sad_8x16[11], sad_8x16[15], p_best_sad_8x32[7], p_best_mv8x32[7]);
    SEARCH_CALC_BEST(sad_8x16[16], sad_8x16[20], p_best_sad_8x32[8], p_best_mv8x32[8]);
    SEARCH_CALC_BEST(sad_8x16[17], sad_8x16[21], p_best_sad_8x32[9], p_best_mv8x32[9]);
    SEARCH_CALC_BEST(sad_8x16[18], sad_8x16[22], p_best_sad_8x32[10], p_best_mv8x32[10]);
    SEARCH_CALC_BEST(sad_8x16[19], sad_8x16[23], p_best_sad_8x32[11], p_best_mv8x32[11]);
    SEARCH_CALC_BEST(sad_8x16[24], sad_8x16[28], p_best_sad_8x32[12], p_best_mv8x32[12]);
    SEARCH_CALC_BEST(sad_8x16[25], sad_8x16[29], p_best_sad_8x32[13], p_best_mv8x32[13]);
    SEARCH_CALC_BEST(sad_8x16[26], sad_8x16[30], p_best_sad_8x32[14], p_best_mv8x32[14]);
    SEARCH_CALC_BEST(sad_8x16[27], sad_8x16[31], p_best_sad_8x32[15], p_best_mv8x32[15]);

#undef SEARCH_CALC_BEST
#undef SEARCH_BEST
#undef ADD_VECS
}
#endif
void ext_eight_sad_calculation_32x32_64x64_avx2(uint32_t  p_sad16x16[16][8],
                                                uint32_t *p_best_sad_32x32,
                                                uint32_t *p_best_sad_64x64,
                                                uint32_t *p_best_mv32x32, uint32_t *p_best_mv64x64,
                                                uint32_t mv, uint32_t p_sad32x32[4][8]) {
    avx2_find_min_pos_init();
    uint32_t si_a, si_b, si_c, si_d, si_e;

    const __m256i tmp0 = _mm256_add_epi32(_mm256_loadu_si256((__m256i const *)p_sad16x16[0]),
                                          _mm256_loadu_si256((__m256i const *)p_sad16x16[1]));

    const __m256i tmp1 = _mm256_add_epi32(_mm256_loadu_si256((__m256i const *)p_sad16x16[2]),
                                          _mm256_loadu_si256((__m256i const *)p_sad16x16[3]));

    const __m256i sad32_a = _mm256_add_epi32(tmp0, tmp1);
    _mm256_storeu_si256((__m256i *)p_sad32x32[0], sad32_a);

    const __m256i tmp2 = _mm256_add_epi32(_mm256_loadu_si256((__m256i const *)p_sad16x16[4]),
                                          _mm256_loadu_si256((__m256i const *)p_sad16x16[5]));

    const __m256i tmp3 = _mm256_add_epi32(_mm256_loadu_si256((__m256i const *)p_sad16x16[6]),
                                          _mm256_loadu_si256((__m256i const *)p_sad16x16[7]));

    const __m256i sad32_b = _mm256_add_epi32(tmp2, tmp3);
    _mm256_storeu_si256((__m256i *)p_sad32x32[1], sad32_b);

    const __m256i tmp4 = _mm256_add_epi32(_mm256_loadu_si256((__m256i const *)p_sad16x16[8]),
                                          _mm256_loadu_si256((__m256i const *)p_sad16x16[9]));

    const __m256i tmp5 = _mm256_add_epi32(_mm256_loadu_si256((__m256i const *)p_sad16x16[10]),
                                          _mm256_loadu_si256((__m256i const *)p_sad16x16[11]));

    const __m256i sad32_c = _mm256_add_epi32(tmp4, tmp5);
    _mm256_storeu_si256((__m256i *)p_sad32x32[2], sad32_c);

    const __m256i tmp6 = _mm256_add_epi32(_mm256_loadu_si256((__m256i const *)p_sad16x16[12]),
                                          _mm256_loadu_si256((__m256i const *)p_sad16x16[13]));

    const __m256i tmp7 = _mm256_add_epi32(_mm256_loadu_si256((__m256i const *)p_sad16x16[14]),
                                          _mm256_loadu_si256((__m256i const *)p_sad16x16[15]));

    const __m256i sad32_d = _mm256_add_epi32(tmp6, tmp7);
    _mm256_storeu_si256((__m256i *)p_sad32x32[3], sad32_d);

    DECLARE_ALIGNED(32, uint32_t, sad64x64[8]);
    const __m256i tmp8     = _mm256_add_epi32(sad32_a, sad32_b);
    const __m256i tmp9     = _mm256_add_epi32(sad32_c, sad32_d);
    *((__m256i *)sad64x64) = _mm256_add_epi32(tmp8, tmp9);

    DECLARE_ALIGNED(32, uint32_t, computed_idx[8]);
    __m256i search_idx_avx2 = _mm256_setr_epi32(0, 4, 8, 12, 16, 20, 24, 28);
    __m256i mv_avx2         = _mm256_set1_epi32(mv);
    __m256i new_mv_avx2     = _mm256_add_epi32(search_idx_avx2, mv_avx2);
    new_mv_avx2             = _mm256_and_si256(new_mv_avx2, _mm256_set1_epi32(0xffff));
    *(__m256i *)computed_idx =
        _mm256_or_si256(new_mv_avx2, _mm256_and_si256(mv_avx2, _mm256_set1_epi32(0xffff0000)));

    avx2_find_min_pos(sad32_a, si_a);
    avx2_find_min_pos(sad32_b, si_b);
    avx2_find_min_pos(sad32_c, si_c);
    avx2_find_min_pos(sad32_d, si_d);
    avx2_find_min_pos(*(__m256i *)sad64x64, si_e);

    if (p_sad32x32[0][si_a] < p_best_sad_32x32[0]) {
        p_best_sad_32x32[0] = p_sad32x32[0][si_a];
        p_best_mv32x32[0]   = computed_idx[si_a];
    }

    if (p_sad32x32[1][si_b] < p_best_sad_32x32[1]) {
        p_best_sad_32x32[1] = p_sad32x32[1][si_b];
        p_best_mv32x32[1]   = computed_idx[si_b];
    }

    if (p_sad32x32[2][si_c] < p_best_sad_32x32[2]) {
        p_best_sad_32x32[2] = p_sad32x32[2][si_c];
        p_best_mv32x32[2]   = computed_idx[si_c];
    }

    if (p_sad32x32[3][si_d] < p_best_sad_32x32[3]) {
        p_best_sad_32x32[3] = p_sad32x32[3][si_d];
        p_best_mv32x32[3]   = computed_idx[si_d];
    }

    if (sad64x64[si_e] < p_best_sad_64x64[0]) {
        p_best_sad_64x64[0] = sad64x64[si_e];
        p_best_mv64x64[0]   = computed_idx[si_e];
    }
}

uint32_t nxm_sad_kernel_sub_sampled_helper_avx2(const uint8_t *src, uint32_t src_stride,
                                                const uint8_t *ref, uint32_t ref_stride,
                                                uint32_t height, uint32_t width) {
    uint32_t nxm_sad = 0;

    switch (width) {
    case 4:
        nxm_sad = eb_compute4x_m_sad_avx2_intrin(src, src_stride, ref, ref_stride, height, width);
        break;
    case 8:
        nxm_sad = eb_compute8x_m_sad_avx2_intrin(src, src_stride, ref, ref_stride, height, width);
        break;
    case 16:
        nxm_sad = eb_compute16x_m_sad_avx2_intrin(src, src_stride, ref, ref_stride, height, width);
        break;
    case 24:
        nxm_sad = eb_compute24x_m_sad_avx2_intrin(src, src_stride, ref, ref_stride, height, width);
        break;
    case 32:
        nxm_sad = eb_compute32x_m_sad_avx2_intrin(src, src_stride, ref, ref_stride, height, width);
        break;
    case 48:
        nxm_sad = eb_compute48x_m_sad_avx2_intrin(src, src_stride, ref, ref_stride, height, width);
        break;
    case 64:
        nxm_sad = eb_compute64x_m_sad_avx2_intrin(src, src_stride, ref, ref_stride, height, width);
        break;
    case 128:
        nxm_sad = eb_compute128x_m_sad_avx2_intrin(src, src_stride, ref, ref_stride, height, width);
        break;
    default: assert(0);
    }

    return nxm_sad;
};

uint32_t nxm_sad_kernel_helper_avx2(const uint8_t *src, uint32_t src_stride, const uint8_t *ref,
                                    uint32_t ref_stride, uint32_t height, uint32_t width) {
    uint32_t nxm_sad = 0;

    switch (width) {
    case 4:
        nxm_sad = eb_compute4x_m_sad_avx2_intrin(src, src_stride, ref, ref_stride, height, width);
        break;
    case 8:
        nxm_sad = eb_compute8x_m_sad_avx2_intrin(src, src_stride, ref, ref_stride, height, width);
        break;
    case 16:
        nxm_sad = eb_compute16x_m_sad_avx2_intrin(src, src_stride, ref, ref_stride, height, width);
        break;
    case 24:
        nxm_sad = eb_compute24x_m_sad_avx2_intrin(src, src_stride, ref, ref_stride, height, width);
        break;
    case 32:
        nxm_sad = eb_compute32x_m_sad_avx2_intrin(src, src_stride, ref, ref_stride, height, width);
        break;
    case 48:
        nxm_sad = eb_compute48x_m_sad_avx2_intrin(src, src_stride, ref, ref_stride, height, width);
        break;
    case 64:
        nxm_sad = eb_compute64x_m_sad_avx2_intrin(src, src_stride, ref, ref_stride, height, width);
        break;
    case 40:
    case 52: break; //void_func();
    default: assert(0);
    }

    return nxm_sad;
}

static INLINE uint32_t
sad_kernel_4xm_16bit_avx2(uint16_t *src, // input parameter, source samples Ptr
                          uint32_t  src_stride, // input parameter, source stride
                          uint16_t *ref, // input parameter, reference samples Ptr
                          uint32_t  ref_stride, // input parameter, reference stride
                          uint32_t  height, // input parameter, block height (M)
                          uint32_t  width) {
    (void)width;
    uint32_t      y;
    __m256i       _sad = _mm256_setzero_si256();
    const __m256i zero = _mm256_setzero_si256();
    assert(((height & 3) == 0) && (width == 4));

    for (y = 0; y < height; y += 4) {
        __m128i src01, src23, ref01, ref23;
        src01 = _mm_loadl_epi64((const __m128i *)src);
        src01 = _mm_castpd_si128(
            _mm_loadh_pd(_mm_castsi128_pd(src01), (double *)(const __m128i *)(src + src_stride)));
        src23              = _mm_loadl_epi64((__m128i *)(const __m128i *)(src + 2 * src_stride));
        src23              = _mm_castpd_si128(_mm_loadh_pd(_mm_castsi128_pd(src23),
                                              (double *)(const __m128i *)(src + 3 * src_stride)));
        const __m256i src0 = _mm256_setr_m128i(src01, src23);

        ref01 = _mm_loadl_epi64((const __m128i *)ref);
        ref01 = _mm_castpd_si128(
            _mm_loadh_pd(_mm_castsi128_pd(ref01), (double *)(const __m128i *)(ref + ref_stride)));
        ref23              = _mm_loadl_epi64((__m128i *)(const __m128i *)(ref + 2 * ref_stride));
        ref23              = _mm_castpd_si128(_mm_loadh_pd(_mm_castsi128_pd(ref23),
                                              (double *)(const __m128i *)(ref + 3 * ref_stride)));
        const __m256i ref0 = _mm256_setr_m128i(ref01, ref23);

        const __m256i min     = _mm256_min_epu16(src0, ref0);
        const __m256i max     = _mm256_max_epu16(src0, ref0);
        const __m256i sad_tmp = _mm256_subs_epu16(max, min);

        const __m256i diff_l = _mm256_unpacklo_epi16(sad_tmp, zero);
        const __m256i diff_h = _mm256_unpackhi_epi16(sad_tmp, zero);

        const __m256i diff = _mm256_add_epi32(diff_l, diff_h);
        _sad               = _mm256_add_epi32(_sad, diff);

        src += 4 * src_stride;
        ref += 4 * ref_stride;
    }

    return sad_final_8_val_avx2(_sad);
}

static INLINE uint32_t
sad_kernel_8xm_16bit_avx2(uint16_t *src, // input parameter, source samples Ptr
                          uint32_t  src_stride, // input parameter, source stride
                          uint16_t *ref, // input parameter, reference samples Ptr
                          uint32_t  ref_stride, // input parameter, reference stride
                          uint32_t  height, // input parameter, block height (M)
                          uint32_t  width) {
    (void)width;
    uint32_t      y;
    __m256i       _sad = _mm256_setzero_si256();
    const __m256i zero = _mm256_setzero_si256();
    assert(((height & 1) == 0) && (width == 8));

    for (y = 0; y < height; y += 2) {
        const __m256i src0 = _mm256_set_m128i(_mm_loadu_si128((const __m128i *)(src + src_stride)),
                                              _mm_loadu_si128((const __m128i *)src));
        const __m256i ref0 = _mm256_set_m128i(_mm_loadu_si128((const __m128i *)(ref + ref_stride)),
                                              _mm_loadu_si128((const __m128i *)ref));

        const __m256i min     = _mm256_min_epu16(src0, ref0);
        const __m256i max     = _mm256_max_epu16(src0, ref0);
        const __m256i sad_tmp = _mm256_subs_epu16(max, min);

        const __m256i diff_l = _mm256_unpacklo_epi16(sad_tmp, zero);
        const __m256i diff_h = _mm256_unpackhi_epi16(sad_tmp, zero);

        const __m256i diff = _mm256_add_epi32(diff_l, diff_h);
        _sad               = _mm256_add_epi32(_sad, diff);

        src += 2 * src_stride;
        ref += 2 * ref_stride;
    }

    return sad_final_8_val_avx2(_sad);
}

static INLINE uint32_t
sad_kernel_16xm_16bit_avx2(uint16_t *src, // input parameter, source samples Ptr
                           uint32_t  src_stride, // input parameter, source stride
                           uint16_t *ref, // input parameter, reference samples Ptr
                           uint32_t  ref_stride, // input parameter, reference stride
                           uint32_t  height, // input parameter, block height (M)
                           uint32_t  width) {
    uint32_t      x, y;
    __m256i       _sad   = _mm256_setzero_si256();
    const __m256i zero   = _mm256_setzero_si256();
    uint32_t      mul_16 = width >> 4;

    for (y = 0; y < height; y++) {
        for (x = 0; x < mul_16; ++x) {
            const __m256i src0    = _mm256_loadu_si256((const __m256i *)(src + x * 16));
            const __m256i ref0    = _mm256_loadu_si256((const __m256i *)(ref + x * 16));
            const __m256i min     = _mm256_min_epu16(src0, ref0);
            const __m256i max     = _mm256_max_epu16(src0, ref0);
            const __m256i sad_tmp = _mm256_subs_epu16(max, min);
            const __m256i diff_l  = _mm256_unpacklo_epi16(sad_tmp, zero);
            const __m256i diff_h  = _mm256_unpackhi_epi16(sad_tmp, zero);
            const __m256i diff    = _mm256_add_epi32(diff_l, diff_h);
            _sad                  = _mm256_add_epi32(_sad, diff);
        }
        src += src_stride;
        ref += ref_stride;
    }

    return sad_final_8_val_avx2(_sad);
}

uint32_t sad_16bit_kernel_avx2(uint16_t *src, // input parameter, source samples Ptr
                               uint32_t  src_stride, // input parameter, source stride
                               uint16_t *ref, // input parameter, reference samples Ptr
                               uint32_t  ref_stride, // input parameter, reference stride
                               uint32_t  height, // input parameter, block height (M)
                               uint32_t  width) // input parameter, block width (N))
{
    uint32_t sad          = 0;
    uint32_t width_offset = width & (~15);

    if (width_offset) {
        sad += sad_kernel_16xm_16bit_avx2(src, src_stride, ref, ref_stride, height, width);
    }
    if ((width - width_offset) > 4) {
        sad += sad_kernel_8xm_16bit_avx2(
            src + width_offset, src_stride, ref + width_offset, ref_stride, height, 8);
        width_offset += 8;
    }
    if (width - width_offset) {
        assert((width - width_offset) == 4);
        sad += sad_kernel_4xm_16bit_avx2(
            src + width_offset, src_stride, ref + width_offset, ref_stride, height, 4);
    }
    return sad;
}
