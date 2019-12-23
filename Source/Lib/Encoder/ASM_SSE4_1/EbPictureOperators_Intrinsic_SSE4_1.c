/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbPictureOperators_SSE4_1.h"
#include "smmintrin.h"

uint64_t compute8x8_satd_u8_sse4(uint8_t * src, // input parameter, diff samples Ptr
                                 uint64_t *dc_value, uint32_t src_stride) {
    uint64_t satd_block_8x8 = 0;
    int16_t  m2[8][8];

    uint32_t j;
    __m128i  s0, s1, s2, s3, s4, s5, s6, s7, s9, s10, s11, s12;
    __m128i  s8 = _mm_setzero_si128();
    __m128i  sum_01_neg, sum_01_pos, sum_23_neg, sum_23_pos, sum_45_neg, sum_45_pos, sum_67_neg,
        sum_67_pos;
    __m128i sum_0to3_pos, sum_4_to_7_pos, sum_0to3_neg, sum_4_to_7_neg, diff_0to3_pos,
        diff_4to7_pos, diff_0to3_neg, diff_4to7_neg;
    __m128i sum0, sum1, difference0, difference1;

    for (j = 0; j < 8; j += 2) {
        s0  = _mm_loadl_epi64((__m128i *)(src + (j * src_stride)));
        s10 = _mm_loadl_epi64((__m128i *)(src + ((j + 1) * src_stride)));
        s10 = _mm_unpacklo_epi8(s10, _mm_setzero_si128());
        s0  = _mm_unpacklo_epi8(s0, _mm_setzero_si128());

        sum0 = _mm_hadd_epi16(s0, s8);
        sum1 = _mm_hadd_epi16(s10, s8);

        difference0 = _mm_hsub_epi16(s0, s8);
        difference1 = _mm_hsub_epi16(s10, s8);

        // m2[j][0]
        // diff[jj] + diff[jj + 4] + diff[jj + 2] + diff[jj + 6] + diff[jj + 1] + diff[jj + 5] + diff[jj + 3] + diff[jj + 7]
        // diff[jj] + diff[jj + 1] + diff[jj + 2] + diff[jj + 3] + diff[jj + 4] + diff[jj + 5] + diff[jj + 6] + diff[jj + 7]
        s1           = _mm_hadd_epi16(sum0, sum1);
        s1           = _mm_hadd_epi16(s1, s8);
        m2[j][0]     = _mm_extract_epi16(s1, 0);
        m2[j + 1][0] = _mm_extract_epi16(s1, 2);

        //m2[j][1]
        //diff[jj] + diff[jj + 4] + diff[jj + 2] + diff[jj + 6] - diff[jj + 1] - diff[jj + 5] - diff[jj + 3] - diff[jj + 7]
        //diff[jj] - diff[jj + 1] + diff[jj + 2] - diff[jj + 3] + diff[jj + 4] - diff[jj + 5] + diff[jj + 6] - diff[jj + 7]
        //(diff[jj] - diff[jj + 1]) + (diff[jj + 2] - diff[jj + 3]) + (diff[jj + 4] - diff[jj + 5]) + (diff[jj + 6] - diff[jj + 7])
        s1           = _mm_hadd_epi16(difference0, difference1);
        s1           = _mm_hadd_epi16(s1, s8);
        m2[j][1]     = _mm_extract_epi16(s1, 0);
        m2[j + 1][1] = _mm_extract_epi16(s1, 2);

        //m2[j][2]
        //diff[jj] + diff[jj + 4] - diff[jj + 2] - diff[jj + 6] + diff[jj + 1] + diff[jj + 5] - diff[jj + 3] - diff[jj + 7]
        //diff[jj] + diff[jj + 1] - diff[jj + 2] - diff[jj + 3] + diff[jj + 4] + diff[jj + 5] - diff[jj + 6] - diff[jj + 7]
        s1           = _mm_hsub_epi16(sum0, sum1);
        s1           = _mm_hadd_epi16(s1, s8);
        m2[j][2]     = _mm_extract_epi16(s1, 0);
        m2[j + 1][2] = _mm_extract_epi16(s1, 2);

        //m2[j][3]
        //diff[jj] + diff[jj + 4] - diff[jj + 2] - diff[jj + 6] - diff[jj + 1] - diff[jj + 5] + diff[jj + 3] + diff[jj + 7]
        //diff[jj] - diff[jj + 1] - diff[jj + 2] + diff[jj + 3] + diff[jj + 4] - diff[jj + 5] - diff[jj + 6] + diff[jj + 7]
        //diff[jj] - diff[jj + 1] - diff[jj + 2] + diff[jj + 3] + diff[jj + 4] - diff[jj + 5] - diff[jj + 6] + diff[jj + 7]
        s1           = _mm_hsub_epi16(difference0, difference1);
        s1           = _mm_hadd_epi16(s1, s8);
        m2[j][3]     = _mm_extract_epi16(s1, 0);
        m2[j + 1][3] = _mm_extract_epi16(s1, 2);

        //m2[j][4]
        //diff[jj] - diff[jj + 4] + diff[jj + 2] - diff[jj + 6] + diff[jj + 1] - diff[jj + 5] + diff[jj + 3] - diff[jj + 7]
        //diff[jj] + diff[jj + 1] + diff[jj + 2] + diff[jj + 3] - diff[jj + 4] - diff[jj + 5] - diff[jj + 6] - diff[jj + 7]
        s1           = _mm_hadd_epi16(sum0, sum1);
        s1           = _mm_hsub_epi16(s1, s8);
        m2[j][4]     = _mm_extract_epi16(s1, 0);
        m2[j + 1][4] = _mm_extract_epi16(s1, 2);

        //m2[j][5]
        //m1[j][4] - m1[j][5]
        //diff[jj] - diff[jj + 1] + diff[jj + 2] - diff[jj + 3] - diff[jj + 4]  + diff[jj + 5] - diff[jj + 6] + diff[jj + 7]
        s1           = _mm_hadd_epi16(difference0, difference1);
        s1           = _mm_hsub_epi16(s1, s8);
        m2[j][5]     = _mm_extract_epi16(s1, 0);
        m2[j + 1][5] = _mm_extract_epi16(s1, 2);

        //m2[j][6]
        //diff[jj] - diff[jj + 4] - diff[jj + 2] + diff[jj + 6] + diff[jj + 1] - diff[jj + 5] - diff[jj + 3] + diff[jj + 7]
        //diff[jj] + diff[jj + 1] - diff[jj + 2] - diff[jj + 3] - diff[jj + 4] - diff[jj + 5] + diff[jj + 6] + diff[jj + 7]

        s1           = _mm_hsub_epi16(sum0, sum1);
        s1           = _mm_hsub_epi16(s1, s8);
        m2[j][6]     = _mm_extract_epi16(s1, 0);
        m2[j + 1][6] = _mm_extract_epi16(s1, 2);

        //m2[j][7]
        //diff[jj] - diff[jj + 4] - diff[jj + 2] + diff[jj + 6] - diff[jj + 1] + diff[jj + 5] + diff[jj + 3] - diff[jj + 7]
        //diff[jj] - diff[jj + 1] - diff[jj + 2] + diff[jj + 3] - diff[jj + 4] + diff[jj + 5] + diff[jj + 6] - diff[jj + 7]
        s1           = _mm_hsub_epi16(difference0, difference1);
        s1           = _mm_hsub_epi16(s1, s8);
        m2[j][7]     = _mm_extract_epi16(s1, 0);
        m2[j + 1][7] = _mm_extract_epi16(s1, 2);
    }

    // Vertical
    s0 = _mm_loadu_si128((__m128i *)(m2[0]));
    s1 = _mm_loadu_si128((__m128i *)(m2[1]));
    s2 = _mm_loadu_si128((__m128i *)(m2[2]));
    s3 = _mm_loadu_si128((__m128i *)(m2[3]));
    s4 = _mm_loadu_si128((__m128i *)(m2[4]));
    s5 = _mm_loadu_si128((__m128i *)(m2[5]));
    s6 = _mm_loadu_si128((__m128i *)(m2[6]));
    s7 = _mm_loadu_si128((__m128i *)(m2[7]));

    sum_01_pos = _mm_add_epi16(s0, s1);
    sum_23_pos = _mm_add_epi16(s2, s3);
    sum_45_pos = _mm_add_epi16(s4, s5);
    sum_67_pos = _mm_add_epi16(s6, s7);

    sum_01_neg = _mm_sub_epi16(s0, s1);
    sum_23_neg = _mm_sub_epi16(s2, s3);
    sum_45_neg = _mm_sub_epi16(s4, s5);
    sum_67_neg = _mm_sub_epi16(s6, s7);

    sum_0to3_pos   = _mm_add_epi16(sum_01_pos, sum_23_pos);
    sum_4_to_7_pos = _mm_add_epi16(sum_45_pos, sum_67_pos);
    diff_0to3_pos  = _mm_sub_epi16(sum_01_pos, sum_23_pos);
    diff_4to7_pos  = _mm_sub_epi16(sum_45_pos, sum_67_pos);

    sum_0to3_neg   = _mm_add_epi16(sum_01_neg, sum_23_neg);
    sum_4_to_7_neg = _mm_add_epi16(sum_45_neg, sum_67_neg);
    diff_0to3_neg  = _mm_sub_epi16(sum_01_neg, sum_23_neg);
    diff_4to7_neg  = _mm_sub_epi16(sum_45_neg, sum_67_neg);

    //m2[0][i] = m1[0][i] + m1[1][i]
    //m2[0][i] = m3[0][i] + m3[2][i] + m3[1][i] + m3[3][i]
    //m2[0][i] = m2[0][i] + m2[4][i] + m2[2][i] + m2[6][i] + m2[1][i] + m2[5][i] + m2[3][i] + m2[7][i]
    //m2[0][i] = m2[0][i] + m2[1][i] + m2[2][i] + m2[3][i] + m2[4][i] + m2[5][i] + m2[6][i] + m2[7][i]
    s9 = _mm_add_epi16(sum_0to3_pos, sum_4_to_7_pos);
    s9 = _mm_abs_epi16(s9);
    *dc_value += _mm_extract_epi16(s9, 0);

    s10 = _mm_unpacklo_epi16(s9, _mm_setzero_si128());
    s11 = _mm_unpackhi_epi16(s9, _mm_setzero_si128());
    s10 = _mm_add_epi32(s10, s11);
    s10 = _mm_hadd_epi32(s10, s8);
    s10 = _mm_hadd_epi32(s10, s8);

    //m2[1][i] = m1[0][i] - m1[1][i]
    //m2[1][i] = m3[0][i] + m3[2][i] -(m3[1][i] + m3[3][i])
    //m2[1][i] = m2[0][i] + m2[4][i] + m2[2][i] + m2[6][i] -(m2[1][i] + m2[5][i] + m2[3][i] + m2[7][i])
    //m2[1][i] = m2[0][i] - m2[1][i] + m2[2][i] - m2[3][i] + m2[4][i] - m2[5][i] + m2[6][i] - m2[7][i]
    s9  = _mm_add_epi16(sum_0to3_neg, sum_4_to_7_neg);
    s9  = _mm_abs_epi16(s9);
    s12 = _mm_unpacklo_epi16(s9, _mm_setzero_si128());
    s11 = _mm_unpackhi_epi16(s9, _mm_setzero_si128());
    s12 = _mm_add_epi32(s12, s11);
    s12 = _mm_hadd_epi32(s12, s8);
    s12 = _mm_hadd_epi32(s12, s8);
    s12 = _mm_add_epi32(s10, s12);

    //m2[2][i] = m1[2][i] + m1[3][i]
    //m2[2][i] = m3[0][i] - m3[2][i] + m3[1][i] - m3[3][i]
    //m2[2][i] = m2[0][i] + m2[4][i] - (m2[2][i] + m2[6][i]) + m2[1][i] + m2[5][i] - (m2[3][i] + m2[7][i])
    //m2[2][i] = m2[0][i] + m2[1][i] - m2[2][i] - m2[3][i] + m2[4][i] + m2[5][i] - m2[6][i] - m2[7][i]
    //m2[2][i] = m2[0][i] + m2[1][i] - (m2[2][i] + m2[3][i]) + m2[4][i] + m2[5][i] - (m2[6][i] + m2[7][i])
    s9  = _mm_add_epi16(diff_0to3_pos, diff_4to7_pos);
    s9  = _mm_abs_epi16(s9);
    s10 = _mm_unpacklo_epi16(s9, _mm_setzero_si128());
    s11 = _mm_unpackhi_epi16(s9, _mm_setzero_si128());
    s10 = _mm_add_epi32(s10, s11);
    s10 = _mm_hadd_epi32(s10, s8);
    s10 = _mm_hadd_epi32(s10, s8);
    s12 = _mm_add_epi32(s10, s12);

    //m2[3][i] = m1[2][i] - m1[3][i]
    //m2[3][i] = m3[0][i] - m3[2][i] - (m3[1][i] - m3[3][i])
    //m2[3][i] = m2[0][i] + m2[4][i] - (m2[2][i] + m2[6][i]) - (m2[1][i] + m2[5][i] - m2[3][i] - m2[7][i])
    //m2[3][i] = m2[0][i] - m2[1][i] - m2[2][i] + m2[3][i] + m2[4][i] - m2[5][i] - m2[6][i] + m2[7][i]
    //m2[3][i] = m2[0][i] - m2[1][i] - (m2[2][i] - m2[3][i]) + (m2[4][i] - m2[5][i]) - (m2[6][i] - m2[7][i])
    s9  = _mm_add_epi16(diff_0to3_neg, diff_4to7_neg);
    s9  = _mm_abs_epi16(s9);
    s10 = _mm_unpacklo_epi16(s9, _mm_setzero_si128());
    s11 = _mm_unpackhi_epi16(s9, _mm_setzero_si128());
    s10 = _mm_add_epi32(s10, s11);
    s10 = _mm_hadd_epi32(s10, s8);
    s10 = _mm_hadd_epi32(s10, s8);
    s12 = _mm_add_epi32(s10, s12);

    //m2[4][i] = m1[4][i] + m1[5][i]
    //m2[4][i] = m3[4][i] + m3[6][i] + m3[5][i] + m3[7][i]
    //m2[4][i] = m2[0][i] - m2[4][i] + m2[2][i] - m2[6][i] + m2[1][i] - m2[5][i] + m2[3][i] - m2[7][i]
    //m2[4][i] = m2[0][i] + m2[1][i] + m2[2][i] + m2[3][i] - m2[4][i] - m2[5][i] - m2[6][i] - m2[7][i]
    //m2[4][i] = m2[0][i] + m2[1][i] + m2[2][i] + m2[3][i] - ( (m2[4][i] + m2[5][i]) + (m2[6][i] + m2[7][i]) )
    s9  = _mm_sub_epi16(sum_0to3_pos, sum_4_to_7_pos);
    s9  = _mm_abs_epi16(s9);
    s10 = _mm_unpacklo_epi16(s9, _mm_setzero_si128());
    s11 = _mm_unpackhi_epi16(s9, _mm_setzero_si128());
    s10 = _mm_add_epi32(s10, s11);
    s10 = _mm_hadd_epi32(s10, s8);
    s10 = _mm_hadd_epi32(s10, s8);
    s12 = _mm_add_epi32(s10, s12);

    //m2[5][i] = m1[4][i] - m1[5][i]
    //m2[5][i] = m3[4][i] + m3[6][i] - (m3[5][i] + m3[7][i])
    //m2[5][i] = m2[0][i] - m2[4][i] + m2[2][i] - m2[6][i] - (m2[1][i] - m2[5][i] + m2[3][i] - m2[7][i])
    //m2[5][i] = m2[0][i] - m2[1][i] + m2[2][i] - m2[3][i] - m2[4][i] + m2[5][i] - m2[6][i] + m2[7][i]
    //m2[5][i] = m2[0][i] - m2[1][i] + (m2[2][i] - m2[3][i]) - ( (m2[4][i] - m2[5][i]) + (m2[6][i] - m2[7][i]) )
    s9  = _mm_sub_epi16(sum_0to3_neg, sum_4_to_7_neg);
    s9  = _mm_abs_epi16(s9);
    s10 = _mm_unpacklo_epi16(s9, _mm_setzero_si128());
    s11 = _mm_unpackhi_epi16(s9, _mm_setzero_si128());
    s10 = _mm_add_epi32(s10, s11);
    s10 = _mm_hadd_epi32(s10, s8);
    s10 = _mm_hadd_epi32(s10, s8);
    s12 = _mm_add_epi32(s10, s12);

    //m2[6][i] = m1[6][i] + m1[7][i]
    //m2[6][i] = m3[4][i] - m3[6][i] + m3[5][i] - m3[7][i]
    //m2[6][i] = m2[0][i] - m2[4][i] - (m2[2][i] - m2[6][i]) + m2[1][i] - m2[5][i] - (m2[3][i] - m2[7][i])
    //m2[6][i] = m2[0][i] + m2[1][i] - m2[2][i] - m2[3][i] - m2[4][i] - m2[5][i] + m2[6][i] + m2[7][i]
    //m2[6][i] = (m2[0][i] + m2[1][i]) - (m2[2][i] + m2[3][i]) - ( (m2[4][i] + m2[5][i]) - (m2[6][i] + m2[7][i]) )
    s9  = _mm_sub_epi16(diff_0to3_pos, diff_4to7_pos);
    s9  = _mm_abs_epi16(s9);
    s10 = _mm_unpacklo_epi16(s9, _mm_setzero_si128());
    s11 = _mm_unpackhi_epi16(s9, _mm_setzero_si128());
    s10 = _mm_add_epi32(s10, s11);
    s10 = _mm_hadd_epi32(s10, s8);
    s10 = _mm_hadd_epi32(s10, s8);
    s12 = _mm_add_epi32(s10, s12);

    //m2[7][i] = m1[6][i] - m1[7][i]
    //m2[7][i] = m3[4][i] - m3[6][i] - (m3[5][i] - m3[7][i])
    //m2[7][i] = m2[0][i] - m2[4][i] - (m2[2][i] - m2[6][i]) - ((m2[1][i] - m2[5][i]) - (m2[3][i] - m2[7][i]))
    //m2[7][i] = (m2[0][i] - m2[1][i]) - (m2[2][i] - m2[3][i]) - ( (m2[4][i] - m2[5][i]) - (m2[6][i] - m2[7][i]) )
    s9  = _mm_sub_epi16(diff_0to3_neg, diff_4to7_neg);
    s9  = _mm_abs_epi16(s9);
    s10 = _mm_unpacklo_epi16(s9, _mm_setzero_si128());
    s11 = _mm_unpackhi_epi16(s9, _mm_setzero_si128());
    s10 = _mm_add_epi32(s10, s11);
    s10 = _mm_hadd_epi32(s10, s8);
    s10 = _mm_hadd_epi32(s10, s8);
    s12 = _mm_add_epi32(s10, s12);

    satd_block_8x8 = (uint64_t)_mm_extract_epi32(s12, 0);

    satd_block_8x8 = ((satd_block_8x8 + 2) >> 2);

    return satd_block_8x8;
}
