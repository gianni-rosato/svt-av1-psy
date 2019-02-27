/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbDefinitions.h"
#include "emmintrin.h"
#include "EbMcp_SSE2.h"
#include "EbIntrinMacros_SSE2.h"
#include "EbIntraPrediction_SSE2.h"

#define OFFSET_17TO31   0
#define OFFSET_1TO15    (8+OFFSET_17TO31)
#define OFFSET_25TO32   (8+OFFSET_1TO15)
#define OFFSET_17TO24   (8+OFFSET_25TO32)
#define OFFSET_9TO16    (8+OFFSET_17TO24)
#define OFFSET_1TO8     (8+OFFSET_9TO16)
#define OFFSET_31TO24   (8+OFFSET_1TO8)
#define OFFSET_23TO16   (8+OFFSET_31TO24)
#define OFFSET_15TO8    (8+OFFSET_23TO16)
#define OFFSET_7TO0     (8+OFFSET_15TO8)
#define OFFSET_3TO0     (4+OFFSET_7TO0) // not a separate entry
#define OFFSET_0        (8+OFFSET_7TO0)
#define OFFSET_1        (8+OFFSET_0)
#define OFFSET_2        (8+OFFSET_1)
#define OFFSET_3        (8+OFFSET_2)
#define OFFSET_4        (8+OFFSET_3)
#define OFFSET_5        (8+OFFSET_4)
#define OFFSET_6        (8+OFFSET_5)
#define OFFSET_7        (8+OFFSET_6)
#define OFFSET_8        (8+OFFSET_7)
#define OFFSET_9        (8+OFFSET_8)
#define OFFSET_10       (8+OFFSET_9)
#define OFFSET_11       (8+OFFSET_10)
#define OFFSET_12       (8+OFFSET_11)
#define OFFSET_13       (8+OFFSET_12)
#define OFFSET_14       (8+OFFSET_13)
#define OFFSET_15       (8+OFFSET_14)
#define OFFSET_16       (8+OFFSET_15)
#define OFFSET_17       (8+OFFSET_16)
#define OFFSET_18       (8+OFFSET_17)
#define OFFSET_19       (8+OFFSET_18)
#define OFFSET_20       (8+OFFSET_19)
#define OFFSET_21       (8+OFFSET_20)
#define OFFSET_22       (8+OFFSET_21)
#define OFFSET_23       (8+OFFSET_22)
#define OFFSET_24       (8+OFFSET_23)
#define OFFSET_25       (8+OFFSET_24)
#define OFFSET_26       (8+OFFSET_25)
#define OFFSET_27       (8+OFFSET_26)
#define OFFSET_28       (8+OFFSET_27)
#define OFFSET_29       (8+OFFSET_28)
#define OFFSET_30       (8+OFFSET_29)
#define OFFSET_31       (8+OFFSET_30)
#define OFFSET_32       (8+OFFSET_31)

#define MACRO_VERTICAL_LUMA_4(A, B, C) \
    *(uint32_t*)prediction_ptr = _mm_cvtsi128_si32(_mm_or_si128(_mm_and_si128(A, B), C)); \
    A = _mm_srli_si128(A, 1); \
    *(uint32_t*)(prediction_ptr + pStride) = _mm_cvtsi128_si32(_mm_or_si128(_mm_and_si128(A, B), C)); \
    A = _mm_srli_si128(A, 1);

#define MACRO_HORIZONTAL_LUMA_32X16(A)\
    left0_14_even = _mm_and_si128(_mm_loadu_si128((__m128i *)(ref_samples+leftOffset+A)), skip_mask);\
    left0_14_even = _mm_packus_epi16(left0_14_even, left0_14_even);\
    left0_14_even = _mm_unpacklo_epi8(left0_14_even, left0_14_even);\
    left0_6_even = _mm_unpacklo_epi16(left0_14_even, left0_14_even);\
    left8_14_even = _mm_unpackhi_epi16(left0_14_even, left0_14_even);\
    left02 = _mm_unpacklo_epi32(left0_6_even, left0_6_even);\
    left46 = _mm_unpackhi_epi32(left0_6_even, left0_6_even);\
    left8_10 = _mm_unpacklo_epi32(left8_14_even, left8_14_even);\
    left12_14 = _mm_unpackhi_epi32(left8_14_even, left8_14_even);\
    left0 = _mm_unpacklo_epi64(left02, left02);\
    left2 = _mm_unpackhi_epi64(left02, left02);\
    left4 = _mm_unpacklo_epi64(left46, left46);\
    left6 = _mm_unpackhi_epi64(left46, left46);\
    left8 = _mm_unpacklo_epi64(left8_10, left8_10);\
    left10 = _mm_unpackhi_epi64(left8_10, left8_10);\
    left12 = _mm_unpacklo_epi64(left12_14, left12_14);\
    left14 = _mm_unpackhi_epi64(left12_14, left12_14);\
    _mm_storeu_si128((__m128i *)(prediction_ptr), left0);\
    _mm_storeu_si128((__m128i *)(prediction_ptr+16), left0);\
    _mm_storeu_si128((__m128i *)(prediction_ptr+pStride), left2);\
    _mm_storeu_si128((__m128i *)(prediction_ptr+pStride+16), left2);\
    _mm_storeu_si128((__m128i *)(prediction_ptr+2*pStride), left4);\
    _mm_storeu_si128((__m128i *)(prediction_ptr+2*pStride+16), left4);\
    _mm_storeu_si128((__m128i *)(prediction_ptr+3*pStride), left6);\
    _mm_storeu_si128((__m128i *)(prediction_ptr+3*pStride+16), left6);\
    prediction_ptr += (pStride << 2);\
    _mm_storeu_si128((__m128i *)(prediction_ptr), left8);\
    _mm_storeu_si128((__m128i *)(prediction_ptr+16), left8);\
    _mm_storeu_si128((__m128i *)(prediction_ptr+pStride), left10);\
    _mm_storeu_si128((__m128i *)(prediction_ptr+pStride+16), left10);\
    _mm_storeu_si128((__m128i *)(prediction_ptr+2*pStride), left12);\
    _mm_storeu_si128((__m128i *)(prediction_ptr+2*pStride+16), left12);\
    _mm_storeu_si128((__m128i *)(prediction_ptr+3*pStride), left14);\
    _mm_storeu_si128((__m128i *)(prediction_ptr+3*pStride+16), left14);


void intra_mode_vertical_luma_sse2_intrin(
    const uint32_t      size,                   //input parameter, denotes the size of the current PU
    uint8_t            *ref_samples,             //input parameter, pointer to the reference samples
    uint8_t            *prediction_ptr,          //output parameter, pointer to the prediction
    const uint32_t      prediction_buffer_stride, //input parameter, denotes the stride for the prediction ptr
    const EbBool     skip                    //skip one row
)
{
    uint32_t topLeftOffset = size << 1;
    uint32_t leftOffset = 0;
    uint32_t topOffset = (size << 1) + 1;
    __m128i xmm0;
    uint32_t pStride = prediction_buffer_stride;

    if (size != 32) {
        __m128i xmm_mask1, xmm_mask2, xmm_topLeft, xmm_topLeft_lo, xmm_topLeft_hi, xmm_mask_skip, xmm_top, xmm_left, xmm_left_lo, xmm_left_hi;

        xmm0 = _mm_setzero_si128();
        xmm_mask1 = _mm_slli_si128(_mm_set1_epi8((int8_t)0xFF), 1);
        xmm_mask2 = _mm_srli_si128(xmm_mask1, 15);

        xmm_topLeft = _mm_set_epi16((signed short)0xffff, (signed short)0xffff, (signed short)0xffff, (signed short)0xffff, (signed short)0xffff, (signed short)0xffff, (signed short)0xffff, *(uint16_t*)(ref_samples + topLeftOffset));
        xmm_topLeft = _mm_unpacklo_epi8(xmm_topLeft, xmm0);
        xmm_topLeft = _mm_unpacklo_epi16(xmm_topLeft, xmm_topLeft);
        xmm_topLeft = _mm_unpacklo_epi32(xmm_topLeft, xmm_topLeft);
        xmm_topLeft_hi = _mm_unpackhi_epi64(xmm_topLeft, xmm_topLeft);
        xmm_topLeft_lo = _mm_unpacklo_epi64(xmm_topLeft, xmm_topLeft);

        if (!skip) {

            if (size == 8) {

                xmm_left = _mm_packus_epi16(_mm_add_epi16(_mm_srai_epi16(_mm_sub_epi16(_mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(ref_samples + leftOffset)), xmm0), xmm_topLeft_lo), 1), xmm_topLeft_hi), xmm0);
                xmm_top = _mm_and_si128(_mm_loadl_epi64((__m128i *)(ref_samples + topOffset)), xmm_mask1);
                MACRO_VERTICAL_LUMA_8(xmm_left, xmm_mask2, xmm_top)
                    prediction_ptr += (pStride << 2);
                MACRO_VERTICAL_LUMA_8(xmm_left, xmm_mask2, xmm_top)
            }
            else if (size == 16) {
                xmm_left = _mm_loadu_si128((__m128i *)(ref_samples + leftOffset));
                xmm_left_lo = _mm_add_epi16(_mm_srai_epi16(_mm_sub_epi16(_mm_unpacklo_epi8(xmm_left, xmm0), xmm_topLeft_lo), 1), xmm_topLeft_hi);
                xmm_left_hi = _mm_add_epi16(_mm_srai_epi16(_mm_sub_epi16(_mm_unpackhi_epi8(xmm_left, xmm0), xmm_topLeft_lo), 1), xmm_topLeft_hi);
                xmm_left = _mm_packus_epi16(xmm_left_lo, xmm_left_hi);
                xmm_top = _mm_and_si128(_mm_loadu_si128((__m128i *)(ref_samples + topOffset)), xmm_mask1);
                MACRO_VERTICAL_LUMA_16(xmm_left, xmm_mask2, xmm_top)
                    prediction_ptr += (pStride << 2);
                MACRO_VERTICAL_LUMA_16(xmm_left, xmm_mask2, xmm_top)
                    prediction_ptr += (pStride << 2);
                MACRO_VERTICAL_LUMA_16(xmm_left, xmm_mask2, xmm_top)
                    prediction_ptr += (pStride << 2);
                MACRO_VERTICAL_LUMA_16(xmm_left, xmm_mask2, xmm_top)
            }
            else {
                xmm_left = _mm_packus_epi16(_mm_add_epi16(_mm_srai_epi16(_mm_sub_epi16(_mm_unpacklo_epi8(_mm_cvtsi32_si128(*(uint32_t*)(ref_samples + leftOffset)), xmm0), xmm_topLeft_lo), 1), xmm_topLeft_hi), xmm0);
                xmm_top = _mm_and_si128(_mm_loadl_epi64((__m128i *)(ref_samples + topOffset)), xmm_mask1);
                MACRO_VERTICAL_LUMA_4(xmm_left, xmm_mask2, xmm_top)
                    prediction_ptr += (pStride << 1);
                MACRO_VERTICAL_LUMA_4(xmm_left, xmm_mask2, xmm_top)
            }
        }
        else {
            pStride <<= 1;
            xmm_mask_skip = _mm_set1_epi16(0x00FF);

            if (size == 8) {
                xmm_left = _mm_packus_epi16(_mm_add_epi16(_mm_srai_epi16(_mm_sub_epi16(_mm_and_si128(_mm_loadl_epi64((__m128i *)(ref_samples + leftOffset)), xmm_mask_skip), xmm_topLeft_lo), 1), xmm_topLeft_hi), xmm0);
                xmm_top = _mm_and_si128(_mm_loadl_epi64((__m128i *)(ref_samples + topOffset)), xmm_mask1);
                MACRO_VERTICAL_LUMA_8(xmm_left, xmm_mask2, xmm_top)
            }
            else if (size == 16) {
                xmm_left = _mm_packus_epi16(_mm_add_epi16(_mm_srai_epi16(_mm_sub_epi16(_mm_and_si128(_mm_loadu_si128((__m128i *)(ref_samples + leftOffset)), xmm_mask_skip), xmm_topLeft_lo), 1), xmm_topLeft_hi), xmm0);
                xmm_top = _mm_and_si128(_mm_loadu_si128((__m128i *)(ref_samples + topOffset)), xmm_mask1);
                MACRO_VERTICAL_LUMA_16(xmm_left, xmm_mask2, xmm_top)
                    prediction_ptr += (pStride << 2);
                MACRO_VERTICAL_LUMA_16(xmm_left, xmm_mask2, xmm_top)
            }
            else {
                xmm_left = _mm_packus_epi16(_mm_add_epi16(_mm_srai_epi16(_mm_sub_epi16(_mm_and_si128(_mm_cvtsi32_si128(*(uint32_t*)(ref_samples + leftOffset)), xmm_mask_skip), xmm_topLeft_lo), 1), xmm_topLeft_hi), xmm0);
                xmm_top = _mm_and_si128(_mm_cvtsi32_si128(*(uint32_t*)(ref_samples + topOffset)), xmm_mask1);
                MACRO_VERTICAL_LUMA_4(xmm_left, xmm_mask2, xmm_top)
            }
        }
    }
    else {
        __m128i top0_15, top16_31;
        uint64_t size_to_write;
        uint32_t count;

        // Each 2 storeu calls stores 32 bytes. Hence each iteration stores 8 * 32 bytes.
        // Depending on skip, we need 4 or 2 iterations to store 32x32 bytes.
        size_to_write = 4 >> (skip ? 1 : 0);
        pStride = pStride << (skip ? 1 : 0);

        top0_15 = _mm_loadu_si128((__m128i *)(ref_samples + topOffset));
        top16_31 = _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 16));

        for (count = 0; count < size_to_write; count++) {
            _mm_storeu_si128((__m128i *)(prediction_ptr), top0_15);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 16), top16_31);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), top0_15);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride + 16), top16_31);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), top0_15);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride + 16), top16_31);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), top0_15);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride + 16), top16_31);
            prediction_ptr += (pStride << 2);
            _mm_storeu_si128((__m128i *)(prediction_ptr), top0_15);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 16), top16_31);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), top0_15);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride + 16), top16_31);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), top0_15);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride + 16), top16_31);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), top0_15);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride + 16), top16_31);
            prediction_ptr += (pStride << 2);
        }
    }

    return;
}

void intra_mode_vertical_chroma_sse2_intrin(
    const uint32_t      size,                   //input parameter, denotes the size of the current PU
    uint8_t            *ref_samples,             //input parameter, pointer to the reference samples
    uint8_t            *prediction_ptr,          //output parameter, pointer to the prediction
    const uint32_t      prediction_buffer_stride, //input parameter, denotes the stride for the prediction ptr
    const EbBool     skip)                    //skip one row
{
    uint32_t pStride = prediction_buffer_stride;
    uint32_t topOffset = (size << 1) + 1;

    if (!skip) {

        if (size == 16) {
            __m128i xmm0 = _mm_loadu_si128((__m128i *)(ref_samples + topOffset));
            _mm_storeu_si128((__m128i *)prediction_ptr, xmm0);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), xmm0);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), xmm0);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), xmm0);
            prediction_ptr = prediction_ptr + (pStride << 2);
            _mm_storeu_si128((__m128i *)prediction_ptr, xmm0);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), xmm0);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), xmm0);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), xmm0);
            prediction_ptr = prediction_ptr + (pStride << 2);
            _mm_storeu_si128((__m128i *)prediction_ptr, xmm0);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), xmm0);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), xmm0);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), xmm0);
            prediction_ptr = prediction_ptr + (pStride << 2);
            _mm_storeu_si128((__m128i *)prediction_ptr, xmm0);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), xmm0);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), xmm0);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), xmm0);
        }
        else if (size == 8) {
            __m128i xmm0 = _mm_loadl_epi64((__m128i *)(ref_samples + topOffset));
            _mm_storel_epi64((__m128i *)(prediction_ptr), xmm0);
            _mm_storel_epi64((__m128i *)(prediction_ptr + pStride), xmm0);
            _mm_storel_epi64((__m128i *)(prediction_ptr + 2 * pStride), xmm0);
            _mm_storel_epi64((__m128i *)(prediction_ptr + 3 * pStride), xmm0);
            prediction_ptr = prediction_ptr + (pStride << 2);
            _mm_storel_epi64((__m128i *)prediction_ptr, xmm0);
            _mm_storel_epi64((__m128i *)(prediction_ptr + pStride), xmm0);
            _mm_storel_epi64((__m128i *)(prediction_ptr + 2 * pStride), xmm0);
            _mm_storel_epi64((__m128i *)(prediction_ptr + 3 * pStride), xmm0);
        }
        else {
            uint32_t top = *(uint32_t*)(ref_samples + topOffset);
            *(uint32_t*)(prediction_ptr) = top;
            *(uint32_t*)(prediction_ptr + pStride) = top;
            *(uint32_t*)(prediction_ptr + 2 * pStride) = top;
            *(uint32_t*)(prediction_ptr + 3 * pStride) = top;
        }
    }
    else {
        pStride <<= 1;
        if (size == 16) {

            __m128i xmm0 = _mm_loadu_si128((__m128i *)(ref_samples + topOffset));
            _mm_storeu_si128((__m128i *)(prediction_ptr), xmm0);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), xmm0);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), xmm0);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), xmm0);
            prediction_ptr = prediction_ptr + (pStride << 2);
            _mm_storeu_si128((__m128i *)(prediction_ptr), xmm0);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), xmm0);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), xmm0);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), xmm0);
        }
        else if (size == 8) {

            __m128i xmm0 = _mm_loadl_epi64((__m128i *)(ref_samples + topOffset));
            _mm_storel_epi64((__m128i *)prediction_ptr, xmm0);
            _mm_storel_epi64((__m128i *)(prediction_ptr + pStride), xmm0);
            _mm_storel_epi64((__m128i *)(prediction_ptr + 2 * pStride), xmm0);
            _mm_storel_epi64((__m128i *)(prediction_ptr + 3 * pStride), xmm0);
        }
        else {
            uint32_t top = *(uint32_t*)(ref_samples + topOffset);
            *(uint32_t*)(prediction_ptr) = top;
            *(uint32_t*)(prediction_ptr + pStride) = top;
        }
    }
}

void intra_mode_horizontal_luma_sse2_intrin(
    const uint32_t      size,                   //input parameter, denotes the size of the current PU
    uint8_t            *ref_samples,             //input parameter, pointer to the reference samples
    uint8_t            *prediction_ptr,          //output parameter, pointer to the prediction
    const uint32_t      prediction_buffer_stride, //input parameter, denotes the stride for the prediction ptr
    const EbBool     skip)                    //skip one row
{
    uint32_t topOffset = (size << 1) + 1;
    uint32_t topLeftOffset = (size << 1);
    uint32_t leftOffset = 0;
    uint32_t pStride = prediction_buffer_stride;

    __m128i xmm_0, top, topLeft, filter_cmpnt0_7, unclipped, clipped;
    __m128i left0_15, left0_7, left8_15, left0_3, left4_7, left8_11, left12_15, left01, left23, left45, left67, left89, left10_11;
    __m128i left12_13, left14_15, left0, left1, left2, left3, left4, left5, left6, left7, left8, left9, left10, left11;
    __m128i left12, left13, left14, left15;

    if (size != 32) {

        topLeft = _mm_cvtsi32_si128(*(ref_samples + topLeftOffset));
        xmm_0 = _mm_setzero_si128();
        topLeft = _mm_unpacklo_epi8(topLeft, xmm_0);
        topLeft = _mm_unpacklo_epi16(topLeft, topLeft);
        topLeft = _mm_unpacklo_epi32(topLeft, topLeft);
        topLeft = _mm_unpacklo_epi64(topLeft, topLeft);

        if (!skip) {

            if (size == 16) {
                __m128i filter_cmpnt8_15;

                top = _mm_loadu_si128((__m128i *)(ref_samples + topOffset));

                filter_cmpnt0_7 = _mm_srai_epi16(_mm_sub_epi16(_mm_unpacklo_epi8(top, xmm_0), topLeft), 1);
                filter_cmpnt8_15 = _mm_srai_epi16(_mm_sub_epi16(_mm_unpackhi_epi8(top, xmm_0), topLeft), 1);

                left0_15 = _mm_loadu_si128((__m128i *)(ref_samples + leftOffset));
                left0_7 = _mm_unpacklo_epi8(left0_15, left0_15);
                left8_15 = _mm_unpackhi_epi8(left0_15, left0_15);

                left0_3 = _mm_unpacklo_epi16(left0_7, left0_7);
                left4_7 = _mm_unpackhi_epi16(left0_7, left0_7);
                left8_11 = _mm_unpacklo_epi16(left8_15, left8_15);
                left12_15 = _mm_unpackhi_epi16(left8_15, left8_15);

                left01 = _mm_unpacklo_epi32(left0_3, left0_3);
                left23 = _mm_unpackhi_epi32(left0_3, left0_3);
                left45 = _mm_unpacklo_epi32(left4_7, left4_7);
                left67 = _mm_unpackhi_epi32(left4_7, left4_7);
                left89 = _mm_unpacklo_epi32(left8_11, left8_11);
                left10_11 = _mm_unpackhi_epi32(left8_11, left8_11);
                left12_13 = _mm_unpacklo_epi32(left12_15, left12_15);
                left14_15 = _mm_unpackhi_epi32(left12_15, left12_15);

                left0 = _mm_unpacklo_epi8(left01, xmm_0);
                clipped = _mm_packus_epi16(_mm_add_epi16(left0, filter_cmpnt0_7), _mm_add_epi16(filter_cmpnt8_15, left0));

                left1 = _mm_unpackhi_epi64(left01, left01);
                left2 = _mm_unpacklo_epi64(left23, left23);
                left3 = _mm_unpackhi_epi64(left23, left23);
                left4 = _mm_unpacklo_epi64(left45, left45);
                left5 = _mm_unpackhi_epi64(left45, left45);
                left6 = _mm_unpacklo_epi64(left67, left67);
                left7 = _mm_unpackhi_epi64(left67, left67);
                left8 = _mm_unpacklo_epi64(left89, left89);
                left9 = _mm_unpackhi_epi64(left89, left89);
                left10 = _mm_unpacklo_epi64(left10_11, left10_11);
                left11 = _mm_unpackhi_epi64(left10_11, left10_11);
                left12 = _mm_unpacklo_epi64(left12_13, left12_13);
                left13 = _mm_unpackhi_epi64(left12_13, left12_13);
                left14 = _mm_unpacklo_epi64(left14_15, left14_15);
                left15 = _mm_unpackhi_epi64(left14_15, left14_15);

                _mm_storeu_si128((__m128i *)prediction_ptr, clipped);
                _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), left1);
                _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), left2);
                _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), left3);
                prediction_ptr += (pStride << 2);
                _mm_storeu_si128((__m128i *)prediction_ptr, left4);
                _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), left5);
                _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), left6);
                _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), left7);
                prediction_ptr += (pStride << 2);
                _mm_storeu_si128((__m128i *)prediction_ptr, left8);
                _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), left9);
                _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), left10);
                _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), left11);
                prediction_ptr += (pStride << 2);
                _mm_storeu_si128((__m128i *)prediction_ptr, left12);
                _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), left13);
                _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), left14);
                _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), left15);
            }
            else if (size == 8) {

                filter_cmpnt0_7 = _mm_srai_epi16(_mm_sub_epi16(_mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(ref_samples + topOffset)), xmm_0), topLeft), 1);

                left0_7 = _mm_loadl_epi64((__m128i *)(ref_samples + leftOffset));
                left0_7 = _mm_unpacklo_epi8(left0_7, left0_7);

                left0_3 = _mm_unpacklo_epi16(left0_7, left0_7);
                left4_7 = _mm_unpackhi_epi16(left0_7, left0_7);

                left01 = _mm_unpacklo_epi32(left0_3, left0_3);
                left23 = _mm_unpackhi_epi32(left0_3, left0_3);
                left45 = _mm_unpacklo_epi32(left4_7, left4_7);
                left67 = _mm_unpackhi_epi32(left4_7, left4_7);

                unclipped = _mm_add_epi16(_mm_unpacklo_epi8(left01, xmm_0), filter_cmpnt0_7);
                clipped = _mm_packus_epi16(unclipped, unclipped);

                _mm_storel_epi64((__m128i *)prediction_ptr, clipped);
                _mm_storel_epi64((__m128i *)(prediction_ptr + pStride), _mm_srli_si128(left01, 8));
                _mm_storel_epi64((__m128i *)(prediction_ptr + 2 * pStride), left23);
                _mm_storel_epi64((__m128i *)(prediction_ptr + 3 * pStride), _mm_srli_si128(left23, 8));
                prediction_ptr += (pStride << 2);
                _mm_storel_epi64((__m128i *)prediction_ptr, left45);
                _mm_storel_epi64((__m128i *)(prediction_ptr + pStride), _mm_srli_si128(left45, 8));
                _mm_storel_epi64((__m128i *)(prediction_ptr + 2 * pStride), left67);
                _mm_storel_epi64((__m128i *)(prediction_ptr + 3 * pStride), _mm_srli_si128(left67, 8));

            }
            else {
                __m128i filter_cmpnt0_3;

                filter_cmpnt0_3 = _mm_srai_epi16(_mm_sub_epi16(_mm_unpacklo_epi8(_mm_cvtsi32_si128(*(uint32_t*)(ref_samples + topOffset)), xmm_0), topLeft), 1);

                left0_3 = _mm_cvtsi32_si128(*(uint32_t*)(ref_samples + leftOffset));
                left0_3 = _mm_unpacklo_epi8(left0_3, left0_3);  //00112233
                left0_3 = _mm_unpacklo_epi16(left0_3, left0_3); //0000111122223333

                left01 = _mm_unpacklo_epi32(left0_3, left0_3);
                left23 = _mm_unpackhi_epi32(left0_3, left0_3);

                unclipped = _mm_add_epi16(_mm_unpacklo_epi8(left01, xmm_0), filter_cmpnt0_3);
                clipped = _mm_packus_epi16(unclipped, unclipped);

                *(uint32_t*)prediction_ptr = _mm_cvtsi128_si32(clipped);
                *(uint32_t*)(prediction_ptr + pStride) = _mm_cvtsi128_si32(_mm_srli_si128(left01, 8));
                *(uint32_t*)(prediction_ptr + 2 * pStride) = _mm_cvtsi128_si32(left23);
                *(uint32_t*)(prediction_ptr + 3 * pStride) = _mm_cvtsi128_si32(_mm_srli_si128(left23, 8));
            }
        }
        else {
            __m128i left0_14_even, left_0_6_even, left02, left46, skip_mask;

            skip_mask = _mm_set1_epi16(0x00FF);
            pStride <<= 1;

            if (size == 16) {
                __m128i filter_cmpnt8_15, left_8_14_even, left8_10, left12_14, left02_16wide;

                top = _mm_loadu_si128((__m128i *)(ref_samples + topOffset));
                filter_cmpnt0_7 = _mm_srai_epi16(_mm_sub_epi16(_mm_unpacklo_epi8(top, xmm_0), topLeft), 1);
                filter_cmpnt8_15 = _mm_srai_epi16(_mm_sub_epi16(_mm_unpackhi_epi8(top, xmm_0), topLeft), 1);

                left0_14_even = _mm_and_si128(_mm_loadu_si128((__m128i *)(ref_samples + leftOffset)), skip_mask);
                left0_14_even = _mm_packus_epi16(left0_14_even, left0_14_even);
                left0_14_even = _mm_unpacklo_epi8(left0_14_even, left0_14_even);

                left_0_6_even = _mm_unpacklo_epi16(left0_14_even, left0_14_even);
                left_8_14_even = _mm_unpackhi_epi16(left0_14_even, left0_14_even);

                left02 = _mm_unpacklo_epi32(left_0_6_even, left_0_6_even);
                left46 = _mm_unpackhi_epi32(left_0_6_even, left_0_6_even);
                left8_10 = _mm_unpacklo_epi32(left_8_14_even, left_8_14_even);
                left12_14 = _mm_unpackhi_epi32(left_8_14_even, left_8_14_even);

                left02_16wide = _mm_unpacklo_epi8(left02, xmm_0);
                clipped = _mm_packus_epi16(_mm_add_epi16(left02_16wide, filter_cmpnt0_7), _mm_add_epi16(filter_cmpnt8_15, left02_16wide));

                left2 = _mm_unpackhi_epi64(left02, left02);
                left4 = _mm_unpacklo_epi64(left46, left46);
                left6 = _mm_unpackhi_epi64(left46, left46);
                left8 = _mm_unpacklo_epi64(left8_10, left8_10);
                left10 = _mm_unpackhi_epi64(left8_10, left8_10);
                left12 = _mm_unpacklo_epi64(left12_14, left12_14);
                left14 = _mm_unpackhi_epi64(left12_14, left12_14);

                _mm_storeu_si128((__m128i *)(prediction_ptr), clipped);
                _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), left2);
                _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), left4);
                _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), left6);
                prediction_ptr += (pStride << 2);
                _mm_storeu_si128((__m128i *)(prediction_ptr), left8);
                _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), left10);
                _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), left12);
                _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), left14);

            }
            else {
                filter_cmpnt0_7 = _mm_srai_epi16(_mm_sub_epi16(_mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(ref_samples + topOffset)), xmm_0), topLeft), 1);

                left_0_6_even = _mm_and_si128(_mm_loadl_epi64((__m128i *)(ref_samples + leftOffset)), skip_mask);
                left_0_6_even = _mm_packus_epi16(left_0_6_even, left_0_6_even);
                left_0_6_even = _mm_unpacklo_epi8(left_0_6_even, left_0_6_even);
                left_0_6_even = _mm_unpacklo_epi16(left_0_6_even, left_0_6_even);

                left02 = _mm_unpacklo_epi32(left_0_6_even, left_0_6_even);
                left46 = _mm_unpackhi_epi32(left_0_6_even, left_0_6_even);

                unclipped = _mm_add_epi16(_mm_unpacklo_epi8(left02, xmm_0), filter_cmpnt0_7);
                clipped = _mm_packus_epi16(unclipped, unclipped);

                _mm_storel_epi64((__m128i *)prediction_ptr, clipped);
                _mm_storel_epi64((__m128i *)(prediction_ptr + pStride), _mm_srli_si128(left02, 8));
                _mm_storel_epi64((__m128i *)(prediction_ptr + 2 * pStride), left46);
                _mm_storel_epi64((__m128i *)(prediction_ptr + 3 * pStride), _mm_srli_si128(left46, 8));
            }
        }
    }
    else {
        if (!skip) {

            uint8_t count;

            for (count = 0; count < 2; ++count) {
                left0_15 = _mm_loadu_si128((__m128i *)(ref_samples + leftOffset));
                ref_samples += 16;

                left0_7 = _mm_unpacklo_epi8(left0_15, left0_15);
                left8_15 = _mm_unpackhi_epi8(left0_15, left0_15);

                left0_3 = _mm_unpacklo_epi16(left0_7, left0_7);
                left4_7 = _mm_unpackhi_epi16(left0_7, left0_7);
                left8_11 = _mm_unpacklo_epi16(left8_15, left8_15);
                left12_15 = _mm_unpackhi_epi16(left8_15, left8_15);

                left01 = _mm_unpacklo_epi32(left0_3, left0_3);
                left23 = _mm_unpackhi_epi32(left0_3, left0_3);
                left45 = _mm_unpacklo_epi32(left4_7, left4_7);
                left67 = _mm_unpackhi_epi32(left4_7, left4_7);
                left89 = _mm_unpacklo_epi32(left8_11, left8_11);
                left10_11 = _mm_unpackhi_epi32(left8_11, left8_11);
                left12_13 = _mm_unpacklo_epi32(left12_15, left12_15);
                left14_15 = _mm_unpackhi_epi32(left12_15, left12_15);

                left0 = _mm_unpacklo_epi64(left01, left01);
                left1 = _mm_unpackhi_epi64(left01, left01);
                left2 = _mm_unpacklo_epi64(left23, left23);
                left3 = _mm_unpackhi_epi64(left23, left23);
                left4 = _mm_unpacklo_epi64(left45, left45);
                left5 = _mm_unpackhi_epi64(left45, left45);
                left6 = _mm_unpacklo_epi64(left67, left67);
                left7 = _mm_unpackhi_epi64(left67, left67);
                left8 = _mm_unpacklo_epi64(left89, left89);
                left9 = _mm_unpackhi_epi64(left89, left89);
                left10 = _mm_unpacklo_epi64(left10_11, left10_11);
                left11 = _mm_unpackhi_epi64(left10_11, left10_11);
                left12 = _mm_unpacklo_epi64(left12_13, left12_13);
                left13 = _mm_unpackhi_epi64(left12_13, left12_13);
                left14 = _mm_unpacklo_epi64(left14_15, left14_15);
                left15 = _mm_unpackhi_epi64(left14_15, left14_15);

                _mm_storeu_si128((__m128i *)prediction_ptr, left0);
                _mm_storeu_si128((__m128i *)(prediction_ptr + 16), left0);
                _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), left1);
                _mm_storeu_si128((__m128i *)(prediction_ptr + pStride + 16), left1);
                _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), left2);
                _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride + 16), left2);
                _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), left3);
                _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride + 16), left3);
                prediction_ptr += (pStride << 2);
                _mm_storeu_si128((__m128i *)prediction_ptr, left4);
                _mm_storeu_si128((__m128i *)(prediction_ptr + 16), left4);
                _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), left5);
                _mm_storeu_si128((__m128i *)(prediction_ptr + pStride + 16), left5);
                _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), left6);
                _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride + 16), left6);
                _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), left7);
                _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride + 16), left7);
                prediction_ptr += (pStride << 2);
                _mm_storeu_si128((__m128i *)prediction_ptr, left8);
                _mm_storeu_si128((__m128i *)(prediction_ptr + 16), left8);
                _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), left9);
                _mm_storeu_si128((__m128i *)(prediction_ptr + pStride + 16), left9);
                _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), left10);
                _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride + 16), left10);
                _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), left11);
                _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride + 16), left11);
                prediction_ptr += (pStride << 2);
                _mm_storeu_si128((__m128i *)prediction_ptr, left12);
                _mm_storeu_si128((__m128i *)(prediction_ptr + 16), left12);
                _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), left13);
                _mm_storeu_si128((__m128i *)(prediction_ptr + pStride + 16), left13);
                _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), left14);
                _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride + 16), left14);
                _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), left15);
                _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride + 16), left15);
                prediction_ptr += (pStride << 2);
            }
        }
        else {
            __m128i left0_14_even, left0_6_even, left8_14_even, left8_10, left12_14, left02, left46, skip_mask;
            skip_mask = _mm_set1_epi16(0x00FF);
            pStride <<= 1;
            MACRO_HORIZONTAL_LUMA_32X16(0)
                prediction_ptr += (pStride << 2);
            MACRO_HORIZONTAL_LUMA_32X16(16)
        }
    }
}


void intra_mode_horizontal_chroma_sse2_intrin(
    const uint32_t      size,                   //input parameter, denotes the size of the current PU
    uint8_t            *ref_samples,             //input parameter, pointer to the reference samples
    uint8_t            *prediction_ptr,          //output parameter, pointer to the prediction
    const uint32_t      prediction_buffer_stride, //input parameter, denotes the stride for the prediction ptr
    const EbBool     skip)                    //skip one row
{
    uint32_t pStride = prediction_buffer_stride;
    uint32_t leftOffset = 0;

    if (!skip) {

        if (size == 16) {
            __m128i xmm0, xmm2, xmm4, xmm6, xmm8, xmm10, xmm12, xmm14;
            xmm0 = _mm_loadu_si128((__m128i *)(ref_samples + leftOffset));
            xmm8 = _mm_unpackhi_epi8(xmm0, xmm0);
            xmm0 = _mm_unpacklo_epi8(xmm0, xmm0);

            xmm4 = _mm_unpackhi_epi16(xmm0, xmm0);
            xmm0 = _mm_unpacklo_epi16(xmm0, xmm0);
            xmm12 = _mm_unpackhi_epi16(xmm8, xmm8);
            xmm8 = _mm_unpacklo_epi16(xmm8, xmm8);

            xmm2 = _mm_unpackhi_epi32(xmm0, xmm0);
            xmm0 = _mm_unpacklo_epi32(xmm0, xmm0);
            xmm6 = _mm_unpackhi_epi32(xmm4, xmm4);
            xmm4 = _mm_unpacklo_epi32(xmm4, xmm4);
            xmm10 = _mm_unpackhi_epi32(xmm8, xmm8);
            xmm8 = _mm_unpacklo_epi32(xmm8, xmm8);
            xmm14 = _mm_unpackhi_epi32(xmm12, xmm12);
            xmm12 = _mm_unpacklo_epi32(xmm12, xmm12);

            _mm_storeu_si128((__m128i *)prediction_ptr, _mm_unpacklo_epi64(xmm0, xmm0));
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_unpackhi_epi64(xmm0, xmm0));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_unpacklo_epi64(xmm2, xmm2));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_unpackhi_epi64(xmm2, xmm2));
            prediction_ptr += (pStride << 2);
            _mm_storeu_si128((__m128i *)prediction_ptr, _mm_unpacklo_epi64(xmm4, xmm4));
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_unpackhi_epi64(xmm4, xmm4));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_unpacklo_epi64(xmm6, xmm6));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_unpackhi_epi64(xmm6, xmm6));
            prediction_ptr += (pStride << 2);
            _mm_storeu_si128((__m128i *)prediction_ptr, _mm_unpacklo_epi64(xmm8, xmm8));
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_unpackhi_epi64(xmm8, xmm8));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_unpacklo_epi64(xmm10, xmm10));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_unpackhi_epi64(xmm10, xmm10));
            prediction_ptr += (pStride << 2);
            _mm_storeu_si128((__m128i *)prediction_ptr, _mm_unpacklo_epi64(xmm12, xmm12));
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_unpackhi_epi64(xmm12, xmm12));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_unpacklo_epi64(xmm14, xmm14));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_unpackhi_epi64(xmm14, xmm14));
        }
        else if (size == 8) {
            __m128i xmm0, xmm2, xmm4, xmm6;
            xmm0 = _mm_loadl_epi64((__m128i *)(ref_samples + leftOffset));
            xmm0 = _mm_unpacklo_epi8(xmm0, xmm0);
            xmm4 = _mm_unpackhi_epi16(xmm0, xmm0);
            xmm0 = _mm_unpacklo_epi16(xmm0, xmm0);

            xmm2 = _mm_unpackhi_epi32(xmm0, xmm0);
            xmm0 = _mm_unpacklo_epi32(xmm0, xmm0);
            xmm6 = _mm_unpackhi_epi32(xmm4, xmm4);
            xmm4 = _mm_unpacklo_epi32(xmm4, xmm4);

            _mm_storel_epi64((__m128i *)(prediction_ptr), xmm0);
            _mm_storel_epi64((__m128i *)(prediction_ptr + pStride), _mm_srli_si128(xmm0, 8));
            _mm_storel_epi64((__m128i *)(prediction_ptr + 2 * pStride), xmm2);
            _mm_storel_epi64((__m128i *)(prediction_ptr + 3 * pStride), _mm_srli_si128(xmm2, 8));
            prediction_ptr += (pStride << 2);
            _mm_storel_epi64((__m128i *)(prediction_ptr), xmm4);
            _mm_storel_epi64((__m128i *)(prediction_ptr + pStride), _mm_srli_si128(xmm4, 8));
            _mm_storel_epi64((__m128i *)(prediction_ptr + 2 * pStride), xmm6);
            _mm_storel_epi64((__m128i *)(prediction_ptr + 3 * pStride), _mm_srli_si128(xmm6, 8));
        }
        else {
            __m128i xmm0 = _mm_cvtsi32_si128(*(uint32_t*)(ref_samples + leftOffset));
            xmm0 = _mm_unpacklo_epi8(xmm0, xmm0);
            xmm0 = _mm_unpacklo_epi16(xmm0, xmm0);
            *(uint32_t*)prediction_ptr = _mm_cvtsi128_si32(xmm0);
            *(uint32_t*)(prediction_ptr + pStride) = _mm_cvtsi128_si32(_mm_srli_si128(xmm0, 4));
            *(uint32_t*)(prediction_ptr + 2 * pStride) = _mm_cvtsi128_si32(_mm_srli_si128(xmm0, 8));
            *(uint32_t*)(prediction_ptr + 3 * pStride) = _mm_cvtsi128_si32(_mm_srli_si128(xmm0, 12));
        }
    }
    else {
        pStride <<= 1;
        __m128i xmm15 = _mm_set1_epi16(0x00FF);
        if (size == 16) {
            __m128i xmm0, xmm2, xmm4, xmm6;

            xmm0 = _mm_and_si128(_mm_loadu_si128((__m128i *)(ref_samples + leftOffset)), xmm15);
            xmm0 = _mm_packus_epi16(xmm0, xmm0);
            xmm0 = _mm_unpacklo_epi8(xmm0, xmm0);
            xmm4 = _mm_unpackhi_epi16(xmm0, xmm0);
            xmm0 = _mm_unpacklo_epi16(xmm0, xmm0);

            xmm2 = _mm_unpackhi_epi32(xmm0, xmm0);
            xmm0 = _mm_unpacklo_epi32(xmm0, xmm0);
            xmm6 = _mm_unpackhi_epi32(xmm4, xmm4);
            xmm4 = _mm_unpacklo_epi32(xmm4, xmm4);

            _mm_storeu_si128((__m128i *)(prediction_ptr), _mm_unpacklo_epi64(xmm0, xmm0));
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_unpackhi_epi64(xmm0, xmm0));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_unpacklo_epi64(xmm2, xmm2));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_unpackhi_epi64(xmm2, xmm2));
            prediction_ptr += (pStride << 2);
            _mm_storeu_si128((__m128i *)(prediction_ptr), _mm_unpacklo_epi64(xmm4, xmm4));
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_unpackhi_epi64(xmm4, xmm4));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_unpacklo_epi64(xmm6, xmm6));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_unpackhi_epi64(xmm6, xmm6));
        }
        else if (size == 8) {

            __m128i xmm2, xmm0;
            xmm0 = _mm_and_si128(_mm_loadl_epi64((__m128i *)(ref_samples + leftOffset)), xmm15);
            xmm0 = _mm_packus_epi16(xmm0, xmm0);
            xmm0 = _mm_unpacklo_epi8(xmm0, xmm0);
            xmm0 = _mm_unpacklo_epi16(xmm0, xmm0);
            xmm2 = _mm_unpackhi_epi32(xmm0, xmm0);
            xmm0 = _mm_unpacklo_epi32(xmm0, xmm0);

            _mm_storel_epi64((__m128i *)(prediction_ptr), xmm0);
            _mm_storel_epi64((__m128i *)(prediction_ptr + pStride), _mm_srli_si128(xmm0, 8));
            _mm_storel_epi64((__m128i *)(prediction_ptr + 2 * pStride), xmm2);
            _mm_storel_epi64((__m128i *)(prediction_ptr + 3 * pStride), _mm_srli_si128(xmm2, 8));
        }
        else {
            __m128i xmm0;
            xmm0 = _mm_and_si128(_mm_cvtsi32_si128(*(uint32_t*)(ref_samples + leftOffset)), xmm15);
            xmm0 = _mm_packus_epi16(xmm0, xmm0);
            xmm0 = _mm_unpacklo_epi8(xmm0, xmm0);
            xmm0 = _mm_unpacklo_epi16(xmm0, xmm0);

            *(uint32_t*)(prediction_ptr) = _mm_cvtsi128_si32(xmm0);
            *(uint32_t*)(prediction_ptr + pStride) = _mm_cvtsi128_si32(_mm_srli_si128(xmm0, 4));
        }
    }
}

void intra_mode_dc_luma_sse2_intrin(
    const uint32_t      size,                       //input parameter, denotes the size of the current PU
    uint8_t            *ref_samples,                 //input parameter, pointer to the reference samples
    uint8_t            *prediction_ptr,              //output parameter, pointer to the prediction
    const uint32_t      prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
    const EbBool     skip)                       //skip one row
{
    __m128i xmm0 = _mm_setzero_si128();
    uint32_t pStride = prediction_buffer_stride;
    uint32_t topOffset = (size << 1) + 1;
    uint32_t leftOffset = 0;

    if (size != 32) {

        __m128i xmm_mask1 = _mm_slli_si128(_mm_set1_epi8((int8_t)0xFF), 1);
        __m128i xmm_mask2 = _mm_srli_si128(xmm_mask1, 15);
        __m128i xmm_C2 = _mm_set1_epi16(0x0002);

        if (!skip) {

            if (size == 16) {
                __m128i xmm_predictionDcValue, xmm_top, xmm_left, xmm_sum, xmm_predictionPtr_0;
                __m128i xmm_top_lo, xmm_top_hi, xmm_left_lo, xmm_left_hi, xmm_predictionDcValue_16, xmm_predictionDcValue_16_x2, xmm_predictionDcValue_16_x3;
                xmm_top = _mm_loadu_si128((__m128i *)(ref_samples + topOffset));
                xmm_left = _mm_loadu_si128((__m128i *)(ref_samples + leftOffset));
                xmm_top_lo = _mm_unpacklo_epi8(xmm_top, xmm0);
                xmm_top_hi = _mm_unpackhi_epi8(xmm_top, xmm0);
                xmm_left_lo = _mm_unpacklo_epi8(xmm_left, xmm0);
                xmm_left_hi = _mm_unpackhi_epi8(xmm_left, xmm0);

                xmm_sum = _mm_add_epi32(_mm_sad_epu8(xmm_top, xmm0), _mm_sad_epu8(xmm_left, xmm0));

                xmm_predictionDcValue = _mm_srli_epi32(_mm_add_epi32(_mm_add_epi32(_mm_srli_si128(xmm_sum, 8), xmm_sum), _mm_cvtsi32_si128(16)), 5);
                xmm_predictionDcValue = _mm_unpacklo_epi8(xmm_predictionDcValue, xmm_predictionDcValue);
                xmm_predictionDcValue = _mm_unpacklo_epi16(xmm_predictionDcValue, xmm_predictionDcValue);
                xmm_predictionDcValue = _mm_unpacklo_epi32(xmm_predictionDcValue, xmm_predictionDcValue);
                xmm_predictionDcValue = _mm_unpacklo_epi64(xmm_predictionDcValue, xmm_predictionDcValue);

                xmm_predictionDcValue_16 = _mm_srli_epi16(xmm_predictionDcValue, 8);
                xmm_predictionDcValue = _mm_and_si128(xmm_predictionDcValue, xmm_mask1);
                xmm_predictionDcValue_16_x2 = _mm_add_epi16(xmm_predictionDcValue_16, xmm_predictionDcValue_16);
                xmm_predictionDcValue_16_x3 = _mm_add_epi16(xmm_predictionDcValue_16_x2, xmm_predictionDcValue_16);

                xmm_top = _mm_packus_epi16(_mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(xmm_top_lo, xmm_predictionDcValue_16_x3), xmm_C2), 2),
                    _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(xmm_top_hi, xmm_predictionDcValue_16_x3), xmm_C2), 2));

                xmm_left = _mm_srli_si128(_mm_packus_epi16(_mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(xmm_left_lo, xmm_predictionDcValue_16_x3), xmm_C2), 2),
                    _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(xmm_left_hi, xmm_predictionDcValue_16_x3), xmm_C2), 2)), 1);

                xmm_predictionPtr_0 = _mm_or_si128(_mm_and_si128(_mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(_mm_add_epi16(xmm_top_lo, xmm_left_lo), xmm_predictionDcValue_16_x2), xmm_C2), 2), xmm_mask2), _mm_and_si128(xmm_top, xmm_mask1));
                _mm_storeu_si128((__m128i *)(prediction_ptr), xmm_predictionPtr_0);

                _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_or_si128(_mm_and_si128(xmm_left, xmm_mask2), xmm_predictionDcValue));
                xmm_left = _mm_srli_si128(xmm_left, 1);
                _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_or_si128(_mm_and_si128(xmm_left, xmm_mask2), xmm_predictionDcValue));
                xmm_left = _mm_srli_si128(xmm_left, 1);
                _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_or_si128(_mm_and_si128(xmm_left, xmm_mask2), xmm_predictionDcValue));
                xmm_left = _mm_srli_si128(xmm_left, 1);

                prediction_ptr += (pStride << 2);
                MACRO_VERTICAL_LUMA_16(xmm_left, xmm_mask2, xmm_predictionDcValue)
                    prediction_ptr += (pStride << 2);
                MACRO_VERTICAL_LUMA_16(xmm_left, xmm_mask2, xmm_predictionDcValue)
                    prediction_ptr += (pStride << 2);
                MACRO_VERTICAL_LUMA_16(xmm_left, xmm_mask2, xmm_predictionDcValue)
            }
            else if (size == 8) {

                __m128i xmm_left, xmm_top, xmm_top_lo, xmm_left_lo, xmm_predictionDcValue, xmm_predictionDcValue_16;
                __m128i xmm_predictionDcValue_16_x2, xmm_predictionDcValue_16_x3, xmm_predictionPtr_0;

                xmm_top = _mm_loadl_epi64((__m128i *)(ref_samples + topOffset));
                xmm_left = _mm_loadl_epi64((__m128i *)(ref_samples + leftOffset));

                xmm_top_lo = _mm_unpacklo_epi8(xmm_top, xmm0);
                xmm_left_lo = _mm_unpacklo_epi8(xmm_left, xmm0);

                xmm_predictionDcValue = _mm_srli_epi32(_mm_add_epi32(_mm_add_epi32(_mm_sad_epu8(xmm_top, xmm0), _mm_sad_epu8(xmm_left, xmm0)), _mm_cvtsi32_si128(8)), 4);
                xmm_predictionDcValue = _mm_unpacklo_epi8(xmm_predictionDcValue, xmm_predictionDcValue);
                xmm_predictionDcValue = _mm_unpacklo_epi16(xmm_predictionDcValue, xmm_predictionDcValue);
                xmm_predictionDcValue = _mm_unpacklo_epi32(xmm_predictionDcValue, xmm_predictionDcValue);
                xmm_predictionDcValue = _mm_unpacklo_epi64(xmm_predictionDcValue, xmm_predictionDcValue);

                xmm_predictionDcValue_16 = _mm_srli_epi16(xmm_predictionDcValue, 8);
                xmm_predictionDcValue = _mm_and_si128(xmm_predictionDcValue, xmm_mask1);

                xmm_predictionDcValue_16_x2 = _mm_add_epi16(xmm_predictionDcValue_16, xmm_predictionDcValue_16);
                xmm_predictionDcValue_16_x3 = _mm_add_epi16(xmm_predictionDcValue_16_x2, xmm_predictionDcValue_16);

                xmm_top = _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(xmm_top_lo, xmm_predictionDcValue_16_x3), xmm_C2), 2);
                xmm_left = _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(xmm_left_lo, xmm_predictionDcValue_16_x3), xmm_C2), 2);

                xmm_left = _mm_srli_si128(_mm_packus_epi16(xmm_left, xmm_left), 1);
                xmm_top = _mm_and_si128(_mm_packus_epi16(xmm_top, xmm_top), xmm_mask1);

                xmm_predictionPtr_0 = _mm_or_si128(_mm_and_si128(_mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(_mm_add_epi16(xmm_top_lo, xmm_left_lo), xmm_predictionDcValue_16_x2), xmm_C2), 2), xmm_mask2), xmm_top);
                _mm_storel_epi64((__m128i *)(prediction_ptr), xmm_predictionPtr_0);

                _mm_storel_epi64((__m128i *)(prediction_ptr + pStride), _mm_or_si128(_mm_and_si128(xmm_left, xmm_mask2), xmm_predictionDcValue));
                xmm_left = _mm_srli_si128(xmm_left, 1);
                _mm_storel_epi64((__m128i *)(prediction_ptr + 2 * pStride), _mm_or_si128(_mm_and_si128(xmm_left, xmm_mask2), xmm_predictionDcValue));
                xmm_left = _mm_srli_si128(xmm_left, 1);
                _mm_storel_epi64((__m128i *)(prediction_ptr + 3 * pStride), _mm_or_si128(_mm_and_si128(xmm_left, xmm_mask2), xmm_predictionDcValue));
                xmm_left = _mm_srli_si128(xmm_left, 1);

                prediction_ptr += (pStride << 2);
                MACRO_VERTICAL_LUMA_8(xmm_left, xmm_mask2, xmm_predictionDcValue)

            }
            else {
                __m128i xmm_left, xmm_top, xmm_top_lo, xmm_left_lo, xmm_predictionDcValue, xmm_predictionDcValue_16;
                __m128i xmm_predictionDcValue_16_x2, xmm_predictionDcValue_16_x3, xmm_predictionPtr_0;

                xmm_top = _mm_cvtsi32_si128(*(uint32_t*)(ref_samples + topOffset));
                xmm_left = _mm_cvtsi32_si128(*(uint32_t*)(ref_samples + leftOffset));

                xmm_top_lo = _mm_unpacklo_epi8(xmm_top, xmm0);
                xmm_left_lo = _mm_unpacklo_epi8(xmm_left, xmm0);
                xmm_predictionDcValue = _mm_srli_epi32(_mm_add_epi32(_mm_add_epi32(_mm_sad_epu8(xmm_top, xmm0), _mm_sad_epu8(xmm_left, xmm0)), _mm_cvtsi32_si128(4)), 3);
                xmm_predictionDcValue = _mm_unpacklo_epi8(xmm_predictionDcValue, xmm_predictionDcValue);
                xmm_predictionDcValue = _mm_unpacklo_epi16(xmm_predictionDcValue, xmm_predictionDcValue);
                xmm_predictionDcValue = _mm_unpacklo_epi32(xmm_predictionDcValue, xmm_predictionDcValue);

                xmm_predictionDcValue_16 = _mm_srli_epi16(xmm_predictionDcValue, 8);
                xmm_predictionDcValue = _mm_and_si128(xmm_predictionDcValue, xmm_mask1);

                xmm_predictionDcValue_16_x2 = _mm_add_epi16(xmm_predictionDcValue_16, xmm_predictionDcValue_16);
                xmm_predictionDcValue_16_x3 = _mm_add_epi16(xmm_predictionDcValue_16_x2, xmm_predictionDcValue_16);

                xmm_top = _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(xmm_top_lo, xmm_predictionDcValue_16_x3), xmm_C2), 2);
                xmm_left = _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(xmm_left_lo, xmm_predictionDcValue_16_x3), xmm_C2), 2);
                xmm_top = _mm_and_si128(_mm_packus_epi16(xmm_top, xmm_top), xmm_mask1);
                xmm_left = _mm_srli_si128(_mm_packus_epi16(xmm_left, xmm_left), 1);

                xmm_predictionPtr_0 = _mm_or_si128(_mm_and_si128(_mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(_mm_add_epi16(xmm_top_lo, xmm_left_lo), xmm_predictionDcValue_16_x2), xmm_C2), 2), xmm_mask2), xmm_top);
                *(uint32_t*)prediction_ptr = _mm_cvtsi128_si32(xmm_predictionPtr_0);

                *(uint32_t*)(prediction_ptr + pStride) = _mm_cvtsi128_si32(_mm_or_si128(_mm_and_si128(xmm_left, xmm_mask2), xmm_predictionDcValue));
                xmm_left = _mm_srli_si128(xmm_left, 1);

                *(uint32_t*)(prediction_ptr + 2 * pStride) = _mm_cvtsi128_si32(_mm_or_si128(_mm_and_si128(xmm_left, xmm_mask2), xmm_predictionDcValue));
                *(uint32_t*)(prediction_ptr + 3 * pStride) = _mm_cvtsi128_si32(_mm_or_si128(_mm_and_si128(_mm_srli_si128(xmm_left, 1), xmm_mask2), xmm_predictionDcValue));
            }
        }
        else {
            pStride <<= 1;

            __m128i xmm_skip_mask = _mm_set1_epi16(0x00FF);
            __m128i xmm_left, xmm_sum, xmm_top, xmm_top_lo, xmm_top_hi, xmm_left_skipped, xmm_predictionDcValue, xmm_predictionDcValue_16;
            __m128i xmm_predictionDcValue_16_x2, xmm_predictionDcValue_16_x3, xmm_predictionPtr_0;

            if (size == 16) {
                xmm_top = _mm_loadu_si128((__m128i *)(ref_samples + topOffset));
                xmm_left = _mm_loadu_si128((__m128i *)(ref_samples + leftOffset));

                xmm_top_lo = _mm_unpacklo_epi8(xmm_top, xmm0);
                xmm_top_hi = _mm_unpackhi_epi8(xmm_top, xmm0);

                xmm_left_skipped = _mm_and_si128(xmm_skip_mask, xmm_left);

                xmm_sum = _mm_add_epi32(_mm_sad_epu8(xmm_top, xmm0), _mm_sad_epu8(xmm_left, xmm0));

                xmm_predictionDcValue = _mm_srli_epi32(_mm_add_epi32(_mm_add_epi32(_mm_srli_si128(xmm_sum, 8), xmm_sum), _mm_cvtsi32_si128(16)), 5);
                xmm_predictionDcValue = _mm_unpacklo_epi8(xmm_predictionDcValue, xmm_predictionDcValue);
                xmm_predictionDcValue = _mm_unpacklo_epi16(xmm_predictionDcValue, xmm_predictionDcValue);
                xmm_predictionDcValue = _mm_unpacklo_epi32(xmm_predictionDcValue, xmm_predictionDcValue);
                xmm_predictionDcValue = _mm_unpacklo_epi64(xmm_predictionDcValue, xmm_predictionDcValue);

                xmm_predictionDcValue_16 = _mm_srli_epi16(xmm_predictionDcValue, 8);
                xmm_predictionDcValue = _mm_and_si128(xmm_predictionDcValue, xmm_mask1);

                xmm_predictionDcValue_16_x2 = _mm_add_epi16(xmm_predictionDcValue_16, xmm_predictionDcValue_16);
                xmm_predictionDcValue_16_x3 = _mm_add_epi16(xmm_predictionDcValue_16_x2, xmm_predictionDcValue_16);

                xmm_left = _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(xmm_left_skipped, xmm_predictionDcValue_16_x3), xmm_C2), 2);
                xmm_left = _mm_srli_si128(_mm_packus_epi16(xmm_left, xmm_left), 1);

                xmm_top = _mm_and_si128(_mm_packus_epi16(_mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(xmm_top_lo, xmm_predictionDcValue_16_x3), xmm_C2), 2),
                    _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(xmm_top_hi, xmm_predictionDcValue_16_x3), xmm_C2), 2)), xmm_mask1);

                xmm_predictionPtr_0 = _mm_or_si128(_mm_and_si128(_mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(_mm_add_epi16(xmm_top_lo, xmm_left_skipped), xmm_predictionDcValue_16_x2), xmm_C2), 2), xmm_mask2), xmm_top);
                _mm_storeu_si128((__m128i *)prediction_ptr, xmm_predictionPtr_0);

                _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_or_si128(_mm_and_si128(xmm_left, xmm_mask2), xmm_predictionDcValue));
                xmm_left = _mm_srli_si128(xmm_left, 1);

                _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_or_si128(_mm_and_si128(xmm_left, xmm_mask2), xmm_predictionDcValue));
                xmm_left = _mm_srli_si128(xmm_left, 1);

                _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_or_si128(_mm_and_si128(xmm_left, xmm_mask2), xmm_predictionDcValue));
                xmm_left = _mm_srli_si128(xmm_left, 1);
                prediction_ptr += (pStride << 2);
                MACRO_VERTICAL_LUMA_16(xmm_left, xmm_mask2, xmm_predictionDcValue)
            }
            else {

                xmm_top = _mm_loadl_epi64((__m128i *)(ref_samples + topOffset));
                xmm_left = _mm_loadl_epi64((__m128i *)(ref_samples + leftOffset));

                xmm_top_lo = _mm_unpacklo_epi8(xmm_top, xmm0);
                xmm_left_skipped = _mm_and_si128(xmm_skip_mask, xmm_left);

                xmm_predictionDcValue = _mm_srli_epi32(_mm_add_epi32(_mm_add_epi32(_mm_sad_epu8(xmm_top, xmm0), _mm_sad_epu8(xmm_left, xmm0)), _mm_cvtsi32_si128(8)), 4);
                xmm_predictionDcValue = _mm_unpacklo_epi8(xmm_predictionDcValue, xmm_predictionDcValue);
                xmm_predictionDcValue = _mm_unpacklo_epi16(xmm_predictionDcValue, xmm_predictionDcValue);
                xmm_predictionDcValue = _mm_unpacklo_epi32(xmm_predictionDcValue, xmm_predictionDcValue);
                xmm_predictionDcValue = _mm_unpacklo_epi64(xmm_predictionDcValue, xmm_predictionDcValue);

                xmm_predictionDcValue_16 = _mm_srli_epi16(xmm_predictionDcValue, 8);
                xmm_predictionDcValue = _mm_and_si128(xmm_predictionDcValue, xmm_mask1);

                xmm_predictionDcValue_16_x2 = _mm_add_epi16(xmm_predictionDcValue_16, xmm_predictionDcValue_16);
                xmm_predictionDcValue_16_x3 = _mm_add_epi16(xmm_predictionDcValue_16_x2, xmm_predictionDcValue_16);

                xmm_top = _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(xmm_top_lo, xmm_predictionDcValue_16_x3), xmm_C2), 2);
                xmm_left = _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(xmm_left_skipped, xmm_predictionDcValue_16_x3), xmm_C2), 2);
                xmm_predictionPtr_0 = _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(_mm_add_epi16(xmm_top_lo, xmm_left_skipped), xmm_predictionDcValue_16_x2), xmm_C2), 2);

                xmm_left = _mm_srli_si128(_mm_packus_epi16(xmm_left, xmm_left), 1);
                xmm_predictionPtr_0 = _mm_or_si128(_mm_and_si128(xmm_predictionPtr_0, xmm_mask2), _mm_and_si128(_mm_packus_epi16(xmm_top, xmm_top), xmm_mask1));
                _mm_storel_epi64((__m128i *)prediction_ptr, xmm_predictionPtr_0);

                _mm_storel_epi64((__m128i *)(prediction_ptr + pStride), _mm_or_si128(_mm_and_si128(xmm_left, xmm_mask2), xmm_predictionDcValue));
                xmm_left = _mm_srli_si128(xmm_left, 1);

                _mm_storel_epi64((__m128i *)(prediction_ptr + 2 * pStride), _mm_or_si128(_mm_and_si128(xmm_left, xmm_mask2), xmm_predictionDcValue));

                _mm_storel_epi64((__m128i *)(prediction_ptr + 3 * pStride), _mm_or_si128(_mm_and_si128(_mm_srli_si128(xmm_left, 1), xmm_mask2), xmm_predictionDcValue));

            }
        }
    }
    else {
        __m128i xmm_sum, xmm_predictionDcValue;

        xmm_sum = _mm_add_epi32(_mm_add_epi32(_mm_sad_epu8(_mm_loadu_si128((__m128i *)(ref_samples + topOffset)), xmm0),
            _mm_sad_epu8(_mm_loadu_si128((__m128i *)(ref_samples + topOffset + 16)), xmm0)),
            _mm_add_epi32(_mm_sad_epu8(_mm_loadu_si128((__m128i *)(ref_samples + leftOffset)), xmm0),
                _mm_sad_epu8(_mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 16)), xmm0)));

        xmm_predictionDcValue = _mm_srli_epi32(_mm_add_epi32(_mm_add_epi32(_mm_srli_si128(xmm_sum, 8), xmm_sum), _mm_cvtsi32_si128(32)), 6);
        xmm_predictionDcValue = _mm_unpacklo_epi8(xmm_predictionDcValue, xmm_predictionDcValue);
        xmm_predictionDcValue = _mm_unpacklo_epi16(xmm_predictionDcValue, xmm_predictionDcValue);
        xmm_predictionDcValue = _mm_unpacklo_epi32(xmm_predictionDcValue, xmm_predictionDcValue);
        xmm_predictionDcValue = _mm_unpacklo_epi64(xmm_predictionDcValue, xmm_predictionDcValue);

        uint32_t count, r1 = (4 >> (skip ? 1 : 0));
        pStride <<= (skip ? 1 : 0);

        for (count = 0; count < r1; ++count) {
            _mm_storeu_si128((__m128i *)prediction_ptr, xmm_predictionDcValue);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 16), xmm_predictionDcValue);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), xmm_predictionDcValue);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride + 16), xmm_predictionDcValue);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), xmm_predictionDcValue);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride + 16), xmm_predictionDcValue);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), xmm_predictionDcValue);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride + 16), xmm_predictionDcValue);
            prediction_ptr += (pStride << 2);
            _mm_storeu_si128((__m128i *)prediction_ptr, xmm_predictionDcValue);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 16), xmm_predictionDcValue);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), xmm_predictionDcValue);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride + 16), xmm_predictionDcValue);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), xmm_predictionDcValue);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride + 16), xmm_predictionDcValue);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), xmm_predictionDcValue);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride + 16), xmm_predictionDcValue);
            prediction_ptr += (pStride << 2);
        }
    }
}


void intra_mode_dc_chroma_sse2_intrin(
    const uint32_t      size,                       //input parameter, denotes the size of the current PU
    uint8_t            *ref_samples,                 //input parameter, pointer to the reference samples
    uint8_t            *prediction_ptr,              //output parameter, pointer to the prediction
    const uint32_t      prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
    const EbBool     skip)                       //skip one row
{
    __m128i xmm0 = _mm_setzero_si128();
    uint32_t pStride = prediction_buffer_stride;
    uint32_t topOffset = (size << 1) + 1;
    uint32_t leftOffset = 0;

    if (!skip) {

        if (size == 16) {
            __m128i sum, predictionDcValue;

            sum = _mm_add_epi32(_mm_sad_epu8(_mm_loadu_si128((__m128i *)(ref_samples + topOffset)), xmm0),
                _mm_sad_epu8(_mm_loadu_si128((__m128i *)(ref_samples + leftOffset)), xmm0));

            predictionDcValue = _mm_srli_epi32(_mm_add_epi32(_mm_add_epi32(_mm_srli_si128(sum, 8), sum), _mm_cvtsi32_si128(16)), 5);
            predictionDcValue = _mm_unpacklo_epi8(predictionDcValue, predictionDcValue);
            predictionDcValue = _mm_unpacklo_epi16(predictionDcValue, predictionDcValue);
            predictionDcValue = _mm_unpacklo_epi32(predictionDcValue, predictionDcValue);
            predictionDcValue = _mm_unpacklo_epi64(predictionDcValue, predictionDcValue);

            _mm_storeu_si128((__m128i *)(prediction_ptr), predictionDcValue);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), predictionDcValue);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), predictionDcValue);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), predictionDcValue);
            prediction_ptr += (pStride << 2);
            _mm_storeu_si128((__m128i *)(prediction_ptr), predictionDcValue);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), predictionDcValue);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), predictionDcValue);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), predictionDcValue);
            prediction_ptr += (pStride << 2);
            _mm_storeu_si128((__m128i *)(prediction_ptr), predictionDcValue);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), predictionDcValue);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), predictionDcValue);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), predictionDcValue);
            prediction_ptr += (pStride << 2);
            _mm_storeu_si128((__m128i *)(prediction_ptr), predictionDcValue);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), predictionDcValue);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), predictionDcValue);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), predictionDcValue);

        }
        else if (size == 8) {
            __m128i sum, predictionDcValue;

            sum = _mm_add_epi32(_mm_sad_epu8(_mm_loadl_epi64((__m128i *)(ref_samples + topOffset)), xmm0),
                _mm_sad_epu8(_mm_loadl_epi64((__m128i *)(ref_samples + leftOffset)), xmm0));

            predictionDcValue = _mm_srli_epi32(_mm_add_epi32(sum, _mm_cvtsi32_si128(8)), 4);
            predictionDcValue = _mm_unpacklo_epi8(predictionDcValue, predictionDcValue);
            predictionDcValue = _mm_unpacklo_epi16(predictionDcValue, predictionDcValue);
            predictionDcValue = _mm_unpacklo_epi32(predictionDcValue, predictionDcValue);

            _mm_storel_epi64((__m128i *)(prediction_ptr), predictionDcValue);
            _mm_storel_epi64((__m128i *)(prediction_ptr + pStride), predictionDcValue);
            _mm_storel_epi64((__m128i *)(prediction_ptr + 2 * pStride), predictionDcValue);
            _mm_storel_epi64((__m128i *)(prediction_ptr + 3 * pStride), predictionDcValue);
            prediction_ptr += (pStride << 2);
            _mm_storel_epi64((__m128i *)(prediction_ptr), predictionDcValue);
            _mm_storel_epi64((__m128i *)(prediction_ptr + pStride), predictionDcValue);
            _mm_storel_epi64((__m128i *)(prediction_ptr + 2 * pStride), predictionDcValue);
            _mm_storel_epi64((__m128i *)(prediction_ptr + 3 * pStride), predictionDcValue);
        }
        else {
            __m128i sum, predictionDcValue;

            sum = _mm_add_epi32(_mm_sad_epu8(_mm_cvtsi32_si128(*(uint32_t*)(ref_samples + topOffset)), xmm0),
                _mm_sad_epu8(_mm_cvtsi32_si128(*(uint32_t*)(ref_samples + leftOffset)), xmm0));

            predictionDcValue = _mm_srli_epi32(_mm_add_epi32(sum, _mm_cvtsi32_si128(4)), 3);
            predictionDcValue = _mm_unpacklo_epi8(predictionDcValue, predictionDcValue);
            predictionDcValue = _mm_unpacklo_epi16(predictionDcValue, predictionDcValue);

            *(uint32_t*)prediction_ptr = _mm_cvtsi128_si32(predictionDcValue);
            *(uint32_t*)(prediction_ptr + pStride) = _mm_cvtsi128_si32(predictionDcValue);
            *(uint32_t*)(prediction_ptr + 2 * pStride) = _mm_cvtsi128_si32(predictionDcValue);
            *(uint32_t*)(prediction_ptr + 3 * pStride) = _mm_cvtsi128_si32(predictionDcValue);
        }
    }
    else {

        pStride <<= 1;

        if (size == 16) {

            __m128i sum, predictionDcValue;

            sum = _mm_add_epi32(_mm_sad_epu8(_mm_loadu_si128((__m128i *)(ref_samples + topOffset)), xmm0),
                _mm_sad_epu8(_mm_loadu_si128((__m128i *)(ref_samples + leftOffset)), xmm0));

            predictionDcValue = _mm_srli_epi32(_mm_add_epi32(_mm_add_epi32(_mm_srli_si128(sum, 8), sum), _mm_cvtsi32_si128(16)), 5);
            predictionDcValue = _mm_unpacklo_epi8(predictionDcValue, predictionDcValue);
            predictionDcValue = _mm_unpacklo_epi16(predictionDcValue, predictionDcValue);
            predictionDcValue = _mm_unpacklo_epi32(predictionDcValue, predictionDcValue);
            predictionDcValue = _mm_unpacklo_epi64(predictionDcValue, predictionDcValue);

            _mm_storeu_si128((__m128i *)(prediction_ptr), predictionDcValue);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), predictionDcValue);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), predictionDcValue);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), predictionDcValue);
            prediction_ptr += (pStride << 2);
            _mm_storeu_si128((__m128i *)(prediction_ptr), predictionDcValue);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), predictionDcValue);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), predictionDcValue);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), predictionDcValue);
        }
        else if (size == 8) {
            __m128i sum, predictionDcValue;

            sum = _mm_add_epi32(_mm_sad_epu8(_mm_loadl_epi64((__m128i *)(ref_samples + topOffset)), xmm0),
                _mm_sad_epu8(_mm_loadl_epi64((__m128i *)(ref_samples + leftOffset)), xmm0));

            predictionDcValue = _mm_srli_epi32(_mm_add_epi32(sum, _mm_cvtsi32_si128(8)), 4);
            predictionDcValue = _mm_unpacklo_epi8(predictionDcValue, predictionDcValue);
            predictionDcValue = _mm_unpacklo_epi16(predictionDcValue, predictionDcValue);
            predictionDcValue = _mm_unpacklo_epi32(predictionDcValue, predictionDcValue);

            _mm_storel_epi64((__m128i *)(prediction_ptr), predictionDcValue);
            _mm_storel_epi64((__m128i *)(prediction_ptr + pStride), predictionDcValue);
            _mm_storel_epi64((__m128i *)(prediction_ptr + 2 * pStride), predictionDcValue);
            _mm_storel_epi64((__m128i *)(prediction_ptr + 3 * pStride), predictionDcValue);
        }
        else {
            __m128i sum, predictionDcValue;

            sum = _mm_add_epi32(_mm_sad_epu8(_mm_cvtsi32_si128(*(uint32_t*)(ref_samples + topOffset)), xmm0),
                _mm_sad_epu8(_mm_cvtsi32_si128(*(uint32_t*)(ref_samples + leftOffset)), xmm0));

            predictionDcValue = _mm_srli_epi32(_mm_add_epi32(sum, _mm_cvtsi32_si128(4)), 3);
            predictionDcValue = _mm_unpacklo_epi8(predictionDcValue, predictionDcValue);
            predictionDcValue = _mm_unpacklo_epi16(predictionDcValue, predictionDcValue);

            *(uint32_t*)prediction_ptr = _mm_cvtsi128_si32(predictionDcValue);
            *(uint32_t*)(prediction_ptr + pStride) = _mm_cvtsi128_si32(predictionDcValue);
        }
    }
}

void intra_mode_planar16bit_sse2_intrin(
    const uint32_t   size,                       //input parameter, denotes the size of the current PU
    uint16_t         *ref_samples,                 //input parameter, pointer to the reference samples
    uint16_t         *prediction_ptr,              //output parameter, pointer to the prediction
    const uint32_t   prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
    const EbBool  skip)                       //skip half rows
{
    uint32_t topOffset = (size << 1) + 1;
    uint32_t leftOffset = 0;
    uint32_t bottomLeftOffset = leftOffset + size;
    uint32_t topRightOffset = topOffset + size;
    uint32_t pStride = prediction_buffer_stride;
    uint32_t count, reverseCnt, coefOffset;
    __m128i  leftCoeff, pred, left, left16, top, topRight, topRightAddSize, bottomLeft, bottomLeftTotal, bottomLeftTotal16;
    EB_ALIGN(16) const int16_t * coeffArray = intra_prediction_const_sse2;

    if (size != 4) {
        __m128i xmm_TopRight, xmm_BottomLeft;
        __m128i pred_0, pred_1, pred_2, pred_3;

        xmm_TopRight = _mm_cvtsi32_si128((uint32_t)*(ref_samples + topRightOffset));
        xmm_TopRight = _mm_unpacklo_epi16(xmm_TopRight, xmm_TopRight);
        xmm_TopRight = _mm_unpacklo_epi32(xmm_TopRight, xmm_TopRight);
        xmm_TopRight = _mm_unpacklo_epi64(xmm_TopRight, xmm_TopRight);
        xmm_BottomLeft = _mm_cvtsi32_si128((uint32_t)*(ref_samples + bottomLeftOffset));
        xmm_BottomLeft = _mm_unpacklo_epi16(xmm_BottomLeft, xmm_BottomLeft);
        xmm_BottomLeft = _mm_unpacklo_epi32(xmm_BottomLeft, xmm_BottomLeft);
        xmm_BottomLeft = _mm_unpacklo_epi64(xmm_BottomLeft, xmm_BottomLeft);

        if (size == 32) {
            __m128i xmm_ref, xmm_ref16;
            __m128i xmm_topRightAddSize0, xmm_topRightAddSize1, xmm_topRightAddSize2, xmm_topRightAddSize3, xmm_top0, xmm_top1, xmm_top2, xmm_top3;

            coefOffset = OFFSET_31; // The coefficients will be size-1, size-2, ... so we start at 31, 30, ... and reverse towards 0

            xmm_topRightAddSize0 = _mm_add_epi16(_mm_mullo_epi16(xmm_TopRight, _mm_load_si128((__m128i *)(coeffArray + OFFSET_1TO8))), _mm_load_si128((__m128i *)(coeffArray + OFFSET_32)));
            xmm_topRightAddSize1 = _mm_add_epi16(_mm_mullo_epi16(xmm_TopRight, _mm_load_si128((__m128i *)(coeffArray + OFFSET_9TO16))), _mm_load_si128((__m128i *)(coeffArray + OFFSET_32)));
            xmm_topRightAddSize2 = _mm_add_epi16(_mm_mullo_epi16(xmm_TopRight, _mm_load_si128((__m128i *)(coeffArray + OFFSET_17TO24))), _mm_load_si128((__m128i *)(coeffArray + OFFSET_32)));
            xmm_topRightAddSize3 = _mm_add_epi16(_mm_mullo_epi16(xmm_TopRight, _mm_load_si128((__m128i *)(coeffArray + OFFSET_25TO32))), _mm_load_si128((__m128i *)(coeffArray + OFFSET_32)));

            xmm_top0 = _mm_loadu_si128((__m128i *)(ref_samples + topOffset));
            xmm_top1 = _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 8));
            xmm_top2 = _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 16));
            xmm_top3 = _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 24));
            ref_samples += leftOffset;

            if (!skip) {
                uint32_t count1;

                // The coefficients will be size-1, size-2, ... so we start at 31, 30, ... and reverse towards 0
                for (reverseCnt = 8; reverseCnt > 0; reverseCnt -= 2) {

                    bottomLeftTotal = _mm_mullo_epi16(xmm_BottomLeft, _mm_load_si128((__m128i *)(coeffArray + 4 * reverseCnt + OFFSET_1TO8 - 32)));

                    xmm_ref = _mm_loadu_si128((__m128i *)(ref_samples));
                    ref_samples += 8;

                    for (count1 = 0; count1 < 8; ++count1) {

                        xmm_ref16 = _mm_unpacklo_epi16(xmm_ref, xmm_ref);
                        xmm_ref16 = _mm_unpacklo_epi32(xmm_ref16, xmm_ref16);
                        xmm_ref16 = _mm_unpacklo_epi64(xmm_ref16, xmm_ref16);

                        bottomLeftTotal16 = _mm_unpacklo_epi16(bottomLeftTotal, bottomLeftTotal);
                        bottomLeftTotal16 = _mm_unpacklo_epi32(bottomLeftTotal16, bottomLeftTotal16);
                        bottomLeftTotal16 = _mm_unpacklo_epi64(bottomLeftTotal16, bottomLeftTotal16);

                        pred_0 = _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(_mm_add_epi16(_mm_mullo_epi16(xmm_ref16, _mm_load_si128((__m128i *)(coeffArray + OFFSET_31TO24))), xmm_topRightAddSize0), bottomLeftTotal16), _mm_mullo_epi16(xmm_top0, _mm_load_si128((__m128i *)(coeffArray + coefOffset)))), 6);
                        pred_1 = _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(_mm_add_epi16(_mm_mullo_epi16(xmm_ref16, _mm_load_si128((__m128i *)(coeffArray + OFFSET_23TO16))), xmm_topRightAddSize1), bottomLeftTotal16), _mm_mullo_epi16(xmm_top1, _mm_load_si128((__m128i *)(coeffArray + coefOffset)))), 6);
                        pred_2 = _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(_mm_add_epi16(_mm_mullo_epi16(xmm_ref16, _mm_load_si128((__m128i *)(coeffArray + OFFSET_15TO8))), xmm_topRightAddSize2), bottomLeftTotal16), _mm_mullo_epi16(xmm_top2, _mm_load_si128((__m128i *)(coeffArray + coefOffset)))), 6);
                        pred_3 = _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(_mm_add_epi16(_mm_mullo_epi16(xmm_ref16, _mm_load_si128((__m128i *)(coeffArray + OFFSET_7TO0))), xmm_topRightAddSize3), bottomLeftTotal16), _mm_mullo_epi16(xmm_top3, _mm_load_si128((__m128i *)(coeffArray + coefOffset)))), 6);

                        _mm_storeu_si128((__m128i *)(prediction_ptr), pred_0);
                        _mm_storeu_si128((__m128i *)(prediction_ptr + 8), pred_1);
                        _mm_storeu_si128((__m128i *)(prediction_ptr + 16), pred_2);
                        _mm_storeu_si128((__m128i *)(prediction_ptr + 24), pred_3);

                        // For next iteration
                        coefOffset -= 8;
                        prediction_ptr += pStride;
                        xmm_ref = _mm_srli_si128(xmm_ref, 2);
                        bottomLeftTotal = _mm_srli_si128(bottomLeftTotal, 2);
                    }
                }
            }
            else {
                for (reverseCnt = 4; reverseCnt > 0; reverseCnt -= 2) {

                    bottomLeftTotal = _mm_mullo_epi16(xmm_BottomLeft, _mm_load_si128((__m128i *)(coeffArray + 4 * reverseCnt + OFFSET_1TO15 - 16)));
                    xmm_ref = _mm_packs_epi32(_mm_srli_epi32(_mm_slli_epi32(_mm_loadu_si128((__m128i *)(ref_samples)), 16), 16),
                        _mm_srli_epi32(_mm_slli_epi32(_mm_loadu_si128((__m128i *)(ref_samples + 8)), 16), 16));
                    ref_samples += 16;

                    for (count = 0; count < 8; ++count) {

                        xmm_ref16 = _mm_unpacklo_epi16(xmm_ref, xmm_ref);
                        xmm_ref16 = _mm_unpacklo_epi32(xmm_ref16, xmm_ref16);
                        xmm_ref16 = _mm_unpacklo_epi64(xmm_ref16, xmm_ref16);

                        bottomLeftTotal16 = _mm_unpacklo_epi16(bottomLeftTotal, bottomLeftTotal);
                        bottomLeftTotal16 = _mm_unpacklo_epi32(bottomLeftTotal16, bottomLeftTotal16);
                        bottomLeftTotal16 = _mm_unpacklo_epi64(bottomLeftTotal16, bottomLeftTotal16);

                        pred_0 = _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(_mm_add_epi16(_mm_mullo_epi16(xmm_ref16, _mm_load_si128((__m128i *)(coeffArray + OFFSET_31TO24))), xmm_topRightAddSize0), bottomLeftTotal16), _mm_mullo_epi16(xmm_top0, _mm_load_si128((__m128i *)(coeffArray + coefOffset)))), 6);
                        pred_1 = _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(_mm_add_epi16(_mm_mullo_epi16(xmm_ref16, _mm_load_si128((__m128i *)(coeffArray + OFFSET_23TO16))), xmm_topRightAddSize1), bottomLeftTotal16), _mm_mullo_epi16(xmm_top1, _mm_load_si128((__m128i *)(coeffArray + coefOffset)))), 6);
                        pred_2 = _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(_mm_add_epi16(_mm_mullo_epi16(xmm_ref16, _mm_load_si128((__m128i *)(coeffArray + OFFSET_15TO8))), xmm_topRightAddSize2), bottomLeftTotal16), _mm_mullo_epi16(xmm_top2, _mm_load_si128((__m128i *)(coeffArray + coefOffset)))), 6);
                        pred_3 = _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(_mm_add_epi16(_mm_mullo_epi16(xmm_ref16, _mm_load_si128((__m128i *)(coeffArray + OFFSET_7TO0))), xmm_topRightAddSize3), bottomLeftTotal16), _mm_mullo_epi16(xmm_top3, _mm_load_si128((__m128i *)(coeffArray + coefOffset)))), 6);

                        _mm_storeu_si128((__m128i *)(prediction_ptr), pred_0);
                        _mm_storeu_si128((__m128i *)(prediction_ptr + 8), pred_1);
                        _mm_storeu_si128((__m128i *)(prediction_ptr + 16), pred_2);
                        _mm_storeu_si128((__m128i *)(prediction_ptr + 24), pred_3);

                        // For next iteration
                        coefOffset -= 16;
                        prediction_ptr += 2 * pStride;
                        bottomLeftTotal = _mm_srli_si128(bottomLeftTotal, 2);
                        xmm_ref = _mm_srli_si128(xmm_ref, 2);
                    }
                }
            }
        }
        else if (size == 16) {
            __m128i topRightTotal0, topRightTotal1, top0, top1, left1, bottomLeftTotal1, leftCoeff0, leftCoeff1;

            coefOffset = OFFSET_15; // The coefficients will be size-1, size-2, ... so we start at 15, 14, ... and reverse towards 0

            topRightTotal0 = _mm_add_epi16(_mm_mullo_epi16(xmm_TopRight, _mm_load_si128((__m128i *)(coeffArray + OFFSET_1TO8))), _mm_load_si128((__m128i *)(coeffArray + OFFSET_16)));
            topRightTotal1 = _mm_add_epi16(_mm_mullo_epi16(xmm_TopRight, _mm_load_si128((__m128i *)(coeffArray + OFFSET_9TO16))), _mm_load_si128((__m128i *)(coeffArray + OFFSET_16)));
            top0 = _mm_loadu_si128((__m128i *)(ref_samples + topOffset));
            top1 = _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 8));
            left = _mm_loadu_si128((__m128i *)(ref_samples + leftOffset));
            left1 = _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 8));
            leftCoeff0 = _mm_load_si128((__m128i *)(coeffArray + OFFSET_15TO8));
            leftCoeff1 = _mm_load_si128((__m128i *)(coeffArray + OFFSET_7TO0));

            if (!skip) {
                uint64_t count1;
                bottomLeftTotal = _mm_mullo_epi16(xmm_BottomLeft, _mm_load_si128((__m128i *)(coeffArray + OFFSET_1TO8)));
                bottomLeftTotal1 = _mm_mullo_epi16(xmm_BottomLeft, _mm_load_si128((__m128i *)(coeffArray + OFFSET_9TO16)));

                for (count = 0; count < 2; ++count) {
                    for (count1 = 0; count1 < 8; ++count1) {

                        left16 = _mm_unpacklo_epi16(left, left);
                        left16 = _mm_unpacklo_epi32(left16, left16);
                        left16 = _mm_unpacklo_epi64(left16, left16);

                        bottomLeftTotal16 = _mm_unpacklo_epi16(bottomLeftTotal, bottomLeftTotal);
                        bottomLeftTotal16 = _mm_unpacklo_epi32(bottomLeftTotal16, bottomLeftTotal16);
                        bottomLeftTotal16 = _mm_unpacklo_epi64(bottomLeftTotal16, bottomLeftTotal16);

                        pred_0 = _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(_mm_add_epi16(_mm_mullo_epi16(left16, leftCoeff0), topRightTotal0), bottomLeftTotal16), _mm_mullo_epi16(top0, _mm_load_si128((__m128i *)(coeffArray + coefOffset)))), 5);
                        pred_1 = _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(_mm_add_epi16(_mm_mullo_epi16(left16, leftCoeff1), topRightTotal1), bottomLeftTotal16), _mm_mullo_epi16(top1, _mm_load_si128((__m128i *)(coeffArray + coefOffset)))), 5);

                        _mm_storeu_si128((__m128i *)prediction_ptr, pred_0);
                        _mm_storeu_si128((__m128i *)(prediction_ptr + 8), pred_1);

                        // Next iteration
                        coefOffset -= 8;
                        prediction_ptr += pStride;
                        left = _mm_srli_si128(left, 2);
                        bottomLeftTotal = _mm_srli_si128(bottomLeftTotal, 2);
                    }
                    bottomLeftTotal = bottomLeftTotal1;
                    left = left1;
                }
            }
            else {
                bottomLeftTotal = _mm_mullo_epi16(xmm_BottomLeft, _mm_load_si128((__m128i *)(coeffArray + OFFSET_1TO15)));
                left = _mm_packs_epi32(_mm_srli_epi32(_mm_slli_epi32(left, 16), 16), _mm_srli_epi32(_mm_slli_epi32(left1, 16), 16));

                for (count = 0; count < 8; ++count) {

                    left16 = _mm_unpacklo_epi16(left, left);
                    left16 = _mm_unpacklo_epi32(left16, left16);
                    left16 = _mm_unpacklo_epi64(left16, left16);

                    bottomLeftTotal16 = _mm_unpacklo_epi16(bottomLeftTotal, bottomLeftTotal);
                    bottomLeftTotal16 = _mm_unpacklo_epi32(bottomLeftTotal16, bottomLeftTotal16);
                    bottomLeftTotal16 = _mm_unpacklo_epi64(bottomLeftTotal16, bottomLeftTotal16);

                    pred_0 = _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(_mm_add_epi16(_mm_mullo_epi16(left16, leftCoeff0), topRightTotal0), bottomLeftTotal16), _mm_mullo_epi16(top0, _mm_load_si128((__m128i *)(coeffArray + coefOffset)))), 5);
                    pred_1 = _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(_mm_add_epi16(_mm_mullo_epi16(left16, leftCoeff1), topRightTotal1), bottomLeftTotal16), _mm_mullo_epi16(top1, _mm_load_si128((__m128i *)(coeffArray + coefOffset)))), 5);

                    _mm_storeu_si128((__m128i *)prediction_ptr, pred_0);
                    _mm_storeu_si128((__m128i *)(prediction_ptr + 8), pred_1);

                    // Next iteration
                    coefOffset -= 16;
                    prediction_ptr += 2 * pStride;
                    bottomLeftTotal = _mm_srli_si128(bottomLeftTotal, 2);
                    left = _mm_srli_si128(left, 2);
                }
            }
        }
        else {

            topRightAddSize = _mm_add_epi16(_mm_mullo_epi16(xmm_TopRight, _mm_load_si128((__m128i *)(coeffArray + OFFSET_1TO8))), _mm_load_si128((__m128i *)(coeffArray + OFFSET_8)));
            top = _mm_loadu_si128((__m128i *)(ref_samples + topOffset));
            left = _mm_loadu_si128((__m128i *)(ref_samples + leftOffset));
            leftCoeff = _mm_load_si128((__m128i *)(coeffArray + OFFSET_7TO0));

            if (!skip) {

                bottomLeftTotal = _mm_mullo_epi16(xmm_BottomLeft, _mm_load_si128((__m128i *)(coeffArray + OFFSET_1TO8)));
                for (count = 0; count < 8; ++count) {

                    left16 = _mm_unpacklo_epi16(left, left);
                    left16 = _mm_unpacklo_epi32(left16, left16);
                    left16 = _mm_unpacklo_epi64(left16, left16);

                    bottomLeftTotal16 = _mm_unpacklo_epi16(bottomLeftTotal, bottomLeftTotal);
                    bottomLeftTotal16 = _mm_unpacklo_epi32(bottomLeftTotal16, bottomLeftTotal16);
                    bottomLeftTotal16 = _mm_unpacklo_epi64(bottomLeftTotal16, bottomLeftTotal16);

                    pred = _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(_mm_add_epi16(_mm_mullo_epi16(left16, leftCoeff), topRightAddSize), bottomLeftTotal16), _mm_mullo_epi16(top, _mm_load_si128((__m128i *)(coeffArray + OFFSET_7)))), 4);
                    _mm_storeu_si128((__m128i *)(prediction_ptr), pred);

                    // Next iteration
                    coeffArray -= 8;
                    prediction_ptr += pStride;
                    left = _mm_srli_si128(left, 2);
                    bottomLeftTotal = _mm_srli_si128(bottomLeftTotal, 2);
                }
            }
            else {
                bottomLeftTotal = _mm_mullo_epi16(xmm_BottomLeft, _mm_load_si128((__m128i *)(coeffArray + OFFSET_1TO15)));
                left = _mm_srli_epi32(_mm_slli_epi32(left, 16), 16);

                for (count = 0; count < 4; ++count) {

                    left16 = _mm_unpacklo_epi16(left, left);
                    left16 = _mm_unpacklo_epi32(left16, left16);
                    left16 = _mm_unpacklo_epi64(left16, left16);

                    bottomLeftTotal16 = _mm_unpacklo_epi16(bottomLeftTotal, bottomLeftTotal);
                    bottomLeftTotal16 = _mm_unpacklo_epi32(bottomLeftTotal16, bottomLeftTotal16);
                    bottomLeftTotal16 = _mm_unpacklo_epi64(bottomLeftTotal16, bottomLeftTotal16);

                    pred = _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(_mm_add_epi16(_mm_mullo_epi16(left16, leftCoeff), topRightAddSize), bottomLeftTotal16), _mm_mullo_epi16(top, _mm_load_si128((__m128i *)(coeffArray + OFFSET_7)))), 4);
                    _mm_storeu_si128((__m128i *)prediction_ptr, pred);

                    // Next iteration
                    coeffArray -= 16;
                    prediction_ptr += 2 * pStride;
                    left = _mm_srli_si128(left, 4);
                    bottomLeftTotal = _mm_srli_si128(bottomLeftTotal, 2);
                }
            }
        }
    }
    else {

        topRight = _mm_cvtsi32_si128((uint32_t)*(ref_samples + topRightOffset));
        topRight = _mm_unpacklo_epi16(topRight, topRight);
        topRight = _mm_unpacklo_epi32(topRight, topRight);

        bottomLeft = _mm_cvtsi32_si128((uint32_t)*(ref_samples + bottomLeftOffset));
        bottomLeft = _mm_unpacklo_epi16(bottomLeft, bottomLeft);
        bottomLeft = _mm_unpacklo_epi32(bottomLeft, bottomLeft);

        topRightAddSize = _mm_add_epi16(_mm_mullo_epi16(topRight, _mm_load_si128((__m128i *)(coeffArray + OFFSET_1TO8))), _mm_load_si128((__m128i *)(coeffArray + OFFSET_4)));

        top = _mm_loadl_epi64((__m128i *)(ref_samples + topOffset));
        left = _mm_loadl_epi64((__m128i *)(ref_samples + leftOffset));
        leftCoeff = _mm_loadl_epi64((__m128i *)(coeffArray + OFFSET_3TO0));

        if (!skip) {
            bottomLeftTotal = _mm_mullo_epi16(bottomLeft, _mm_load_si128((__m128i *)(coeffArray + OFFSET_1TO8)));

            for (count = 0; count < 4; ++count) {

                left16 = _mm_unpacklo_epi16(left, left);
                left16 = _mm_unpacklo_epi32(left16, left16);

                bottomLeftTotal16 = _mm_unpacklo_epi16(bottomLeftTotal, bottomLeftTotal);
                bottomLeftTotal16 = _mm_unpacklo_epi32(bottomLeftTotal16, bottomLeftTotal16);

                pred = _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(_mm_add_epi16(_mm_mullo_epi16(left16, leftCoeff), topRightAddSize), bottomLeftTotal16), _mm_mullo_epi16(top, _mm_load_si128((__m128i *)(coeffArray + OFFSET_3)))), 3);
                _mm_storel_epi64((__m128i *)(prediction_ptr), pred);

                // Next iteration
                coeffArray -= 8;
                prediction_ptr += pStride;
                bottomLeftTotal = _mm_srli_epi64(bottomLeftTotal, 16);
                left = _mm_srli_epi64(left, 16);
            }
        }
        else {
            bottomLeftTotal = _mm_mullo_epi16(bottomLeft, _mm_load_si128((__m128i *)(coeffArray + OFFSET_1TO15)));
            left = _mm_srli_epi32(_mm_slli_epi32(left, 16), 16);

            for (count = 0; count < 2; ++count) {

                left16 = _mm_unpacklo_epi16(left, left);
                left16 = _mm_unpacklo_epi32(left16, left16);

                bottomLeftTotal16 = _mm_unpacklo_epi16(bottomLeftTotal, bottomLeftTotal);
                bottomLeftTotal16 = _mm_unpacklo_epi32(bottomLeftTotal16, bottomLeftTotal16);

                pred = _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(_mm_add_epi16(_mm_mullo_epi16(left16, leftCoeff), topRightAddSize), bottomLeftTotal16), _mm_mullo_epi16(top, _mm_load_si128((__m128i *)(coeffArray + OFFSET_3)))), 3);
                _mm_storel_epi64((__m128i *)prediction_ptr, pred);

                // Next iteration
                coeffArray -= 16;
                prediction_ptr += 2 * pStride;
                left = _mm_srli_epi64(left, 32);
                bottomLeftTotal = _mm_srli_epi64(bottomLeftTotal, 16);
            }
        }
    }
}



void intra_mode_angular_2_sse2_intrin(
    const uint32_t      size,                       //input parameter, denotes the size of the current PU
    uint8_t            *ref_samples,                 //input parameter, pointer to the reference samples
    uint8_t            *prediction_ptr,              //output parameter, pointer to the prediction
    const uint32_t      prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
    const EbBool     skip)                       //skip one row
{
    uint32_t pStride = prediction_buffer_stride;
    uint32_t leftOffset = 0;

    if (!skip) {

        if (size == 32) {
            uint32_t count;
            for (count = 0; count < 8; ++count) {
                _mm_storeu_si128((__m128i *)prediction_ptr, _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 1)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + 16), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 17)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 2)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + pStride + 16), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 18)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 3)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride + 16), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 19)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 4)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride + 16), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 20)));
                ref_samples += 4;
                prediction_ptr += (pStride << 2);
            }
        }
        else if (size == 16) {
            _mm_storeu_si128((__m128i *)prediction_ptr, _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 1)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 2)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 3)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 4)));
            prediction_ptr += (pStride << 2);
            _mm_storeu_si128((__m128i *)prediction_ptr, _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 5)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 6)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 7)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 8)));
            prediction_ptr += (pStride << 2);
            _mm_storeu_si128((__m128i *)prediction_ptr, _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 9)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 10)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 11)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 12)));
            prediction_ptr += (pStride << 2);
            _mm_storeu_si128((__m128i *)prediction_ptr, _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 13)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 14)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 15)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 16)));
        }
        else if (size == 8) {
            _mm_storel_epi64((__m128i *)prediction_ptr, _mm_loadl_epi64((__m128i *)(ref_samples + leftOffset + 1)));
            _mm_storel_epi64((__m128i *)(prediction_ptr + pStride), _mm_loadl_epi64((__m128i *)(ref_samples + leftOffset + 2)));
            _mm_storel_epi64((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadl_epi64((__m128i *)(ref_samples + leftOffset + 3)));
            _mm_storel_epi64((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadl_epi64((__m128i *)(ref_samples + leftOffset + 4)));
            prediction_ptr += (pStride << 2);
            _mm_storel_epi64((__m128i *)prediction_ptr, _mm_loadl_epi64((__m128i *)(ref_samples + leftOffset + 5)));
            _mm_storel_epi64((__m128i *)(prediction_ptr + pStride), _mm_loadl_epi64((__m128i *)(ref_samples + leftOffset + 6)));
            _mm_storel_epi64((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadl_epi64((__m128i *)(ref_samples + leftOffset + 7)));
            _mm_storel_epi64((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadl_epi64((__m128i *)(ref_samples + leftOffset + 8)));
        }
        else {
            *(uint32_t *)prediction_ptr = *(uint32_t *)(ref_samples + leftOffset + 1);
            *(uint32_t *)(prediction_ptr + pStride) = *(uint32_t *)(ref_samples + leftOffset + 2);
            *(uint32_t *)(prediction_ptr + 2 * pStride) = *(uint32_t *)(ref_samples + leftOffset + 3);
            *(uint32_t *)(prediction_ptr + 3 * pStride) = *(uint32_t *)(ref_samples + leftOffset + 4);
        }
    }
    else {
        if (size != 4) {
            pStride <<= 1;

            if (size == 32) {
                uint32_t count;

                for (count = 0; count < 4; ++count) {

                    _mm_storeu_si128((__m128i *)prediction_ptr, _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 1)));
                    _mm_storeu_si128((__m128i *)(prediction_ptr + 16), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 17)));
                    _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 3)));
                    _mm_storeu_si128((__m128i *)(prediction_ptr + pStride + 16), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 19)));
                    _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 5)));
                    _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride + 16), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 21)));
                    _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 7)));
                    _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride + 16), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 23)));
                    ref_samples += 8;
                    prediction_ptr += (pStride << 2);
                }
            }
            else if (size == 16) {

                _mm_storeu_si128((__m128i *)prediction_ptr, _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 1)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 3)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 5)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 7)));
                prediction_ptr += (pStride << 2);
                _mm_storeu_si128((__m128i *)prediction_ptr, _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 9)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 11)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 13)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 15)));
            }
            else {
                _mm_storel_epi64((__m128i *)prediction_ptr, _mm_loadl_epi64((__m128i *)(ref_samples + leftOffset + 1)));
                _mm_storel_epi64((__m128i *)(prediction_ptr + pStride), _mm_loadl_epi64((__m128i *)(ref_samples + leftOffset + 3)));
                _mm_storel_epi64((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadl_epi64((__m128i *)(ref_samples + leftOffset + 5)));
                _mm_storel_epi64((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadl_epi64((__m128i *)(ref_samples + leftOffset + 7)));
            }
        }
        else {
            *(uint32_t*)prediction_ptr = *(uint32_t*)(ref_samples + leftOffset + 1);
            *(uint32_t*)(prediction_ptr + 2 * pStride) = *(uint32_t*)(ref_samples + leftOffset + 3);
        }
    }
}



void intra_mode_angular_34_sse2_intrin(
    const uint32_t      size,                       //input parameter, denotes the size of the current PU
    uint8_t            *ref_samples,                 //input parameter, pointer to the reference samples
    uint8_t            *prediction_ptr,              //output parameter, pointer to the prediction
    const uint32_t      prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
    const EbBool     skip)                       //skip one row
{
    uint32_t pStride = prediction_buffer_stride;
    uint32_t topOffset = ((size << 1) + 1);

    if (!skip) {

        if (size == 32) {

            uint32_t count;
            for (count = 0; count < 8; ++count) {

                _mm_storeu_si128((__m128i *)prediction_ptr, _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 1)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + 16), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 17)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 2)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + pStride + 16), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 18)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 3)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride + 16), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 19)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 4)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride + 16), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 20)));
                ref_samples += 4;
                prediction_ptr += (pStride << 2);
            }
        }
        else if (size == 16) {
            _mm_storeu_si128((__m128i *)prediction_ptr, _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 1)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 2)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 3)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 4)));
            prediction_ptr += (pStride << 2);
            _mm_storeu_si128((__m128i *)prediction_ptr, _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 5)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 6)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 7)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 8)));
            prediction_ptr += (pStride << 2);
            _mm_storeu_si128((__m128i *)prediction_ptr, _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 9)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 10)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 11)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 12)));
            prediction_ptr += (pStride << 2);
            _mm_storeu_si128((__m128i *)prediction_ptr, _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 13)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 14)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 15)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 16)));
        }
        else if (size == 8) {
            _mm_storel_epi64((__m128i *)prediction_ptr, _mm_loadl_epi64((__m128i *)(ref_samples + topOffset + 1)));
            _mm_storel_epi64((__m128i *)(prediction_ptr + pStride), _mm_loadl_epi64((__m128i *)(ref_samples + topOffset + 2)));
            _mm_storel_epi64((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadl_epi64((__m128i *)(ref_samples + topOffset + 3)));
            _mm_storel_epi64((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadl_epi64((__m128i *)(ref_samples + topOffset + 4)));
            prediction_ptr += (pStride << 2);
            _mm_storel_epi64((__m128i *)prediction_ptr, _mm_loadl_epi64((__m128i *)(ref_samples + topOffset + 5)));
            _mm_storel_epi64((__m128i *)(prediction_ptr + pStride), _mm_loadl_epi64((__m128i *)(ref_samples + topOffset + 6)));
            _mm_storel_epi64((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadl_epi64((__m128i *)(ref_samples + topOffset + 7)));
            _mm_storel_epi64((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadl_epi64((__m128i *)(ref_samples + topOffset + 8)));
        }
        else {
            *(uint32_t *)prediction_ptr = *(uint32_t *)(ref_samples + topOffset + 1);
            *(uint32_t *)(prediction_ptr + pStride) = *(uint32_t *)(ref_samples + topOffset + 2);
            *(uint32_t *)(prediction_ptr + 2 * pStride) = *(uint32_t *)(ref_samples + topOffset + 3);
            *(uint32_t *)(prediction_ptr + 3 * pStride) = *(uint32_t *)(ref_samples + topOffset + 4);
        }
    }
    else {
        if (size != 4) {
            pStride <<= 1;

            if (size == 32) {
                uint32_t count;

                for (count = 0; count < 4; ++count) {
                    _mm_storeu_si128((__m128i *)prediction_ptr, _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 1)));
                    _mm_storeu_si128((__m128i *)(prediction_ptr + 16), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 17)));
                    _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 3)));
                    _mm_storeu_si128((__m128i *)(prediction_ptr + pStride + 16), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 19)));
                    _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 5)));
                    _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride + 16), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 21)));
                    _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 7)));
                    _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride + 16), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 23)));
                    ref_samples += 8;
                    prediction_ptr += (pStride << 2);
                }
            }
            else if (size == 16) {

                _mm_storeu_si128((__m128i *)prediction_ptr, _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 1)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 3)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 5)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 7)));
                prediction_ptr += (pStride << 2);
                _mm_storeu_si128((__m128i *)prediction_ptr, _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 9)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 11)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 13)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 15)));
            }
            else {
                _mm_storel_epi64((__m128i *)prediction_ptr, _mm_loadl_epi64((__m128i *)(ref_samples + topOffset + 1)));
                _mm_storel_epi64((__m128i *)(prediction_ptr + pStride), _mm_loadl_epi64((__m128i *)(ref_samples + topOffset + 3)));
                _mm_storel_epi64((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadl_epi64((__m128i *)(ref_samples + topOffset + 5)));
                _mm_storel_epi64((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadl_epi64((__m128i *)(ref_samples + topOffset + 7)));
            }
        }
        else {
            *(uint32_t*)prediction_ptr = *(uint32_t*)(ref_samples + topOffset + 1);
            *(uint32_t*)(prediction_ptr + 2 * pStride) = *(uint32_t*)(ref_samples + topOffset + 3);
        }
    }
}



void intra_mode_angular_18_sse2_intrin(
    const uint32_t      size,                       //input parameter, denotes the size of the current PU
    uint8_t            *ref_samples,                 //input parameter, pointer to the reference samples
    uint8_t            *prediction_ptr,              //output parameter, pointer to the prediction
    const uint32_t      prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
    const EbBool     skip)                    //skip one row
{
    uint32_t pStride = prediction_buffer_stride;
    uint32_t topLeftOffset = (size << 1);

    if (!skip) {

        if (size == 32) {

            uint32_t count;
            for (count = 0; count < 8; ++count) {

                _mm_storeu_si128((__m128i *)prediction_ptr, _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + 16), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset + 16)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 1)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + pStride + 16), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset + 15)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 2)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride + 16), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset + 14)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 3)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride + 16), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset + 13)));
                ref_samples -= 4;
                prediction_ptr += (pStride << 2);
            }
        }
        else if (size == 16) {
            _mm_storeu_si128((__m128i *)prediction_ptr, _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 1)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 2)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 3)));
            prediction_ptr += (pStride << 2);
            _mm_storeu_si128((__m128i *)prediction_ptr, _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 4)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 5)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 6)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 7)));
            prediction_ptr += (pStride << 2);
            _mm_storeu_si128((__m128i *)prediction_ptr, _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 8)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 9)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 10)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 11)));
            prediction_ptr += (pStride << 2);
            _mm_storeu_si128((__m128i *)prediction_ptr, _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 12)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 13)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 14)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 15)));
        }
        else if (size == 8) {
            _mm_storel_epi64((__m128i *)prediction_ptr, _mm_loadl_epi64((__m128i *)(ref_samples + topLeftOffset)));
            _mm_storel_epi64((__m128i *)(prediction_ptr + pStride), _mm_loadl_epi64((__m128i *)(ref_samples + topLeftOffset - 1)));
            _mm_storel_epi64((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadl_epi64((__m128i *)(ref_samples + topLeftOffset - 2)));
            _mm_storel_epi64((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadl_epi64((__m128i *)(ref_samples + topLeftOffset - 3)));
            prediction_ptr += (pStride << 2);
            _mm_storel_epi64((__m128i *)prediction_ptr, _mm_loadl_epi64((__m128i *)(ref_samples + topLeftOffset - 4)));
            _mm_storel_epi64((__m128i *)(prediction_ptr + pStride), _mm_loadl_epi64((__m128i *)(ref_samples + topLeftOffset - 5)));
            _mm_storel_epi64((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadl_epi64((__m128i *)(ref_samples + topLeftOffset - 6)));
            _mm_storel_epi64((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadl_epi64((__m128i *)(ref_samples + topLeftOffset - 7)));
        }
        else {
            *(uint32_t *)prediction_ptr = *(uint32_t *)(ref_samples + topLeftOffset);
            *(uint32_t *)(prediction_ptr + pStride) = *(uint32_t *)(ref_samples + topLeftOffset - 1);
            *(uint32_t *)(prediction_ptr + 2 * pStride) = *(uint32_t *)(ref_samples + topLeftOffset - 2);
            *(uint32_t *)(prediction_ptr + 3 * pStride) = *(uint32_t *)(ref_samples + topLeftOffset - 3);
        }
    }
    else {
        if (size != 4) {
            pStride <<= 1;

            if (size == 32) {
                uint32_t count;

                for (count = 0; count < 4; ++count) {

                    _mm_storeu_si128((__m128i *)prediction_ptr, _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset)));
                    _mm_storeu_si128((__m128i *)(prediction_ptr + 16), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset + 16)));
                    _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 2)));
                    _mm_storeu_si128((__m128i *)(prediction_ptr + pStride + 16), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset + 14)));
                    _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 4)));
                    _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride + 16), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset + 12)));
                    _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 6)));
                    _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride + 16), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset + 10)));
                    ref_samples -= 8;
                    prediction_ptr += (pStride << 2);
                }
            }
            else if (size == 16) {
                _mm_storeu_si128((__m128i *)prediction_ptr, _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 2)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 4)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 6)));
                prediction_ptr += (pStride << 2);
                _mm_storeu_si128((__m128i *)prediction_ptr, _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 8)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 10)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 12)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 14)));
            }
            else {
                _mm_storel_epi64((__m128i *)prediction_ptr, _mm_loadl_epi64((__m128i *)(ref_samples + topLeftOffset)));
                _mm_storel_epi64((__m128i *)(prediction_ptr + pStride), _mm_loadl_epi64((__m128i *)(ref_samples + topLeftOffset - 2)));
                _mm_storel_epi64((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadl_epi64((__m128i *)(ref_samples + topLeftOffset - 4)));
                _mm_storel_epi64((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadl_epi64((__m128i *)(ref_samples + topLeftOffset - 6)));
            }
        }
        else {
            *(uint32_t*)prediction_ptr = *(uint32_t*)(ref_samples + topLeftOffset);
            *(uint32_t*)(prediction_ptr + 2 * pStride) = *(uint32_t*)(ref_samples + topLeftOffset - 2);
        }
    }
}

void intra_mode_planar_sse2_intrin(
    const uint32_t      size,                       //input parameter, denotes the size of the current PU
    uint8_t            *ref_samples,                 //input parameter, pointer to the reference samples
    uint8_t            *prediction_ptr,              //output parameter, pointer to the prediction
    const uint32_t      prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
    const EbBool     skip)                       //skip one row
{

    uint32_t topOffset = (size << 1) + 1;
    uint32_t leftOffset = 0;
    uint32_t topRightOffset = topOffset + size;
    uint32_t bottomLeftOffset = leftOffset + size;
    uint32_t count;
    __m128i xmm_TRP, xmm_BLP, xmm_LP, xmm_TO;
    __m128i xmm_TRP_tot, xmm_BLP_tot, xmm_LP_tot, xmm_TO_tot;
    __m128i xmm_LP_Coeff, xmm_LP_shift, xmm0 = _mm_setzero_si128(), xmm_sum, xmm_BLP_tot_shift;
    EB_ALIGN(16) const int16_t * IntraPredictionConst = intra_prediction_const_sse2;

    if (size == 4) {

        xmm_TRP = _mm_cvtsi32_si128((uint32_t)*(ref_samples + topRightOffset));
        xmm_TRP = _mm_unpacklo_epi16(xmm_TRP, xmm_TRP);
        xmm_TRP = _mm_unpacklo_epi32(xmm_TRP, xmm_TRP);

        xmm_BLP = _mm_cvtsi32_si128((uint32_t)*(ref_samples + bottomLeftOffset));
        xmm_BLP = _mm_unpacklo_epi16(xmm_BLP, xmm_BLP);
        xmm_BLP = _mm_unpacklo_epi32(xmm_BLP, xmm_BLP);

        xmm_TRP_tot = _mm_mullo_epi16(xmm_TRP, _mm_load_si128((__m128i *)(IntraPredictionConst + OFFSET_1TO8)));
        xmm_TRP_tot = _mm_add_epi16(xmm_TRP_tot, _mm_load_si128((__m128i *)(IntraPredictionConst + OFFSET_4)));

        xmm_TO = _mm_unpacklo_epi8(_mm_cvtsi32_si128(*(uint32_t*)(ref_samples + topOffset)), xmm0);
        xmm_LP = _mm_cvtsi32_si128(*(uint32_t*)(ref_samples + leftOffset));
        xmm_LP_Coeff = _mm_loadl_epi64((__m128i *)(IntraPredictionConst + OFFSET_3TO0));

        if (!skip) {
            xmm_BLP_tot = _mm_mullo_epi16(xmm_BLP, _mm_load_si128((__m128i *)(IntraPredictionConst + OFFSET_1TO8)));
            xmm_LP = _mm_unpacklo_epi8(xmm_LP, xmm0);
            for (count = 0; count < 4; ++count) {
                xmm_LP_shift = xmm_LP;
                xmm_BLP_tot_shift = xmm_BLP_tot;
                xmm_LP_shift = _mm_unpacklo_epi16(xmm_LP_shift, xmm_LP_shift);
                xmm_LP_shift = _mm_unpacklo_epi32(xmm_LP_shift, xmm_LP_shift);
                xmm_sum = _mm_add_epi16(xmm_TRP_tot, _mm_mullo_epi16(xmm_LP_shift, xmm_LP_Coeff));

                xmm_BLP_tot_shift = _mm_unpacklo_epi16(xmm_BLP_tot_shift, xmm_BLP_tot_shift);
                xmm_BLP_tot_shift = _mm_unpacklo_epi32(xmm_BLP_tot_shift, xmm_BLP_tot_shift);
                xmm_sum = _mm_add_epi16(xmm_sum, xmm_BLP_tot_shift);

                xmm_sum = _mm_srli_epi16(_mm_add_epi16(xmm_sum, _mm_mullo_epi16(xmm_TO, _mm_load_si128((__m128i *)(IntraPredictionConst + OFFSET_3)))), 3);
                xmm_sum = _mm_packus_epi16(xmm_sum, xmm_sum);
                *(uint32_t*)prediction_ptr = _mm_cvtsi128_si32(xmm_sum);
                prediction_ptr += prediction_buffer_stride;
                IntraPredictionConst -= 8;

                xmm_LP = _mm_srli_epi64(xmm_LP, 16);
                xmm_BLP_tot = _mm_srli_epi64(xmm_BLP_tot, 16);
            }
        }
        else {
            xmm_BLP_tot = _mm_mullo_epi16(xmm_BLP, _mm_load_si128((__m128i *)(IntraPredictionConst + OFFSET_1TO15)));
            xmm_LP = _mm_srli_epi16(_mm_slli_epi16(xmm_LP, 8), 8);

            for (count = 0; count < 2; ++count) {
                xmm_LP_shift = xmm_LP;
                xmm_LP_shift = _mm_unpacklo_epi16(xmm_LP_shift, xmm_LP_shift);
                xmm_LP_shift = _mm_unpacklo_epi32(xmm_LP_shift, xmm_LP_shift);
                xmm_sum = _mm_add_epi16(xmm_TRP_tot, _mm_mullo_epi16(xmm_LP_shift, xmm_LP_Coeff));

                xmm_BLP_tot_shift = xmm_BLP_tot;
                xmm_BLP_tot_shift = _mm_unpacklo_epi16(xmm_BLP_tot_shift, xmm_BLP_tot_shift);
                xmm_BLP_tot_shift = _mm_unpacklo_epi32(xmm_BLP_tot_shift, xmm_BLP_tot_shift);
                xmm_sum = _mm_add_epi16(xmm_sum, xmm_BLP_tot_shift);

                xmm_sum = _mm_srli_epi16(_mm_add_epi16(xmm_sum, _mm_mullo_epi16(xmm_TO, _mm_load_si128((__m128i *)(IntraPredictionConst + OFFSET_3)))), 3);
                xmm_sum = _mm_packus_epi16(xmm_sum, xmm_sum);
                *(uint32_t*)(prediction_ptr) = _mm_cvtsi128_si32(xmm_sum);

                xmm_LP = _mm_srli_epi64(xmm_LP, 16);
                xmm_BLP_tot = _mm_srli_epi64(xmm_BLP_tot, 16);
                prediction_ptr += (prediction_buffer_stride << 1);
                IntraPredictionConst -= 16;
            }
        }
    }
    else {
        xmm_TRP = _mm_cvtsi32_si128((uint32_t)*(ref_samples + topRightOffset));
        xmm_TRP = _mm_unpacklo_epi16(xmm_TRP, xmm_TRP);
        xmm_TRP = _mm_unpacklo_epi32(xmm_TRP, xmm_TRP);
        xmm_TRP = _mm_unpacklo_epi64(xmm_TRP, xmm_TRP);

        xmm_BLP = _mm_cvtsi32_si128((uint32_t)*(ref_samples + bottomLeftOffset));
        xmm_BLP = _mm_unpacklo_epi16(xmm_BLP, xmm_BLP);
        xmm_BLP = _mm_unpacklo_epi32(xmm_BLP, xmm_BLP);
        xmm_BLP = _mm_unpacklo_epi64(xmm_BLP, xmm_BLP);


        if (size == 32) {
            uint32_t count_in;
            __m128i xmm_TRP_tot_0, xmm_TRP_tot_1, xmm_TRP_tot_2, xmm_TRP_tot_3, xmm_TO_tot_0, xmm_TO_tot_1, xmm_TO_tot_2, xmm_TO_tot_3;
            __m128i xmm_TO_lo, xmm_TO_hi, xmm_sum0, xmm_sum1, xmm_sum2, xmm_sum3;
            uint32_t offset = OFFSET_31;

            xmm_TRP_tot_0 = _mm_add_epi16(_mm_mullo_epi16(xmm_TRP, _mm_load_si128((__m128i *)(IntraPredictionConst + OFFSET_1TO8))), _mm_load_si128((__m128i *)(IntraPredictionConst + OFFSET_32)));
            xmm_TRP_tot_1 = _mm_add_epi16(_mm_mullo_epi16(xmm_TRP, _mm_load_si128((__m128i *)(IntraPredictionConst + OFFSET_9TO16))), _mm_load_si128((__m128i *)(IntraPredictionConst + OFFSET_32)));
            xmm_TRP_tot_2 = _mm_add_epi16(_mm_mullo_epi16(xmm_TRP, _mm_load_si128((__m128i *)(IntraPredictionConst + OFFSET_17TO24))), _mm_load_si128((__m128i *)(IntraPredictionConst + OFFSET_32)));
            xmm_TRP_tot_3 = _mm_add_epi16(_mm_mullo_epi16(xmm_TRP, _mm_load_si128((__m128i *)(IntraPredictionConst + OFFSET_25TO32))), _mm_load_si128((__m128i *)(IntraPredictionConst + OFFSET_32)));
            xmm_TO_lo = _mm_loadu_si128((__m128i *)(ref_samples + topOffset));
            xmm_TO_hi = _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 16));
            ref_samples += leftOffset;

            if (!skip) {
                for (count = 0; count < 8; count += 2) {

                    xmm_BLP_tot = _mm_mullo_epi16(xmm_BLP, _mm_load_si128((__m128i *)(IntraPredictionConst - (count * 4) + OFFSET_1TO8)));
                    xmm_LP = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(ref_samples)), xmm0);
                    ref_samples += 8;

                    for (count_in = 0; count_in < 8; ++count_in) {
                        xmm_LP_shift = xmm_LP;
                        xmm_BLP_tot_shift = xmm_BLP_tot;

                        xmm_LP_shift = _mm_unpacklo_epi16(xmm_LP_shift, xmm_LP_shift);
                        xmm_LP_shift = _mm_unpacklo_epi32(xmm_LP_shift, xmm_LP_shift);
                        xmm_LP_shift = _mm_unpacklo_epi64(xmm_LP_shift, xmm_LP_shift);
                        xmm_sum0 = _mm_add_epi16(xmm_TRP_tot_0, _mm_mullo_epi16(xmm_LP_shift, _mm_load_si128((__m128i *)(IntraPredictionConst + OFFSET_31TO24))));
                        xmm_sum1 = _mm_add_epi16(xmm_TRP_tot_1, _mm_mullo_epi16(xmm_LP_shift, _mm_load_si128((__m128i *)(IntraPredictionConst + OFFSET_23TO16))));
                        xmm_sum2 = _mm_add_epi16(xmm_TRP_tot_2, _mm_mullo_epi16(xmm_LP_shift, _mm_load_si128((__m128i *)(IntraPredictionConst + OFFSET_15TO8))));
                        xmm_sum3 = _mm_add_epi16(xmm_TRP_tot_3, _mm_mullo_epi16(xmm_LP_shift, _mm_load_si128((__m128i *)(IntraPredictionConst + OFFSET_7TO0))));

                        xmm_BLP_tot_shift = _mm_unpacklo_epi16(xmm_BLP_tot_shift, xmm_BLP_tot_shift);
                        xmm_BLP_tot_shift = _mm_unpacklo_epi32(xmm_BLP_tot_shift, xmm_BLP_tot_shift);
                        xmm_BLP_tot_shift = _mm_unpacklo_epi64(xmm_BLP_tot_shift, xmm_BLP_tot_shift);
                        xmm_sum0 = _mm_add_epi16(xmm_sum0, xmm_BLP_tot_shift);
                        xmm_sum1 = _mm_add_epi16(xmm_sum1, xmm_BLP_tot_shift);
                        xmm_sum2 = _mm_add_epi16(xmm_sum2, xmm_BLP_tot_shift);
                        xmm_sum3 = _mm_add_epi16(xmm_sum3, xmm_BLP_tot_shift);

                        xmm_TO_tot_0 = _mm_mullo_epi16(_mm_unpacklo_epi8(xmm_TO_lo, xmm0), _mm_load_si128((__m128i *)(IntraPredictionConst + offset)));
                        xmm_TO_tot_2 = _mm_mullo_epi16(_mm_unpacklo_epi8(xmm_TO_hi, xmm0), _mm_load_si128((__m128i *)(IntraPredictionConst + offset)));
                        xmm_TO_tot_1 = _mm_mullo_epi16(_mm_unpackhi_epi8(xmm_TO_lo, xmm0), _mm_load_si128((__m128i *)(IntraPredictionConst + offset)));
                        xmm_TO_tot_3 = _mm_mullo_epi16(_mm_unpackhi_epi8(xmm_TO_hi, xmm0), _mm_load_si128((__m128i *)(IntraPredictionConst + offset)));
                        xmm_sum0 = _mm_srli_epi16(_mm_add_epi16(xmm_TO_tot_0, xmm_sum0), 6);
                        xmm_sum1 = _mm_srli_epi16(_mm_add_epi16(xmm_TO_tot_1, xmm_sum1), 6);
                        xmm_sum2 = _mm_srli_epi16(_mm_add_epi16(xmm_TO_tot_2, xmm_sum2), 6);
                        xmm_sum3 = _mm_srli_epi16(_mm_add_epi16(xmm_TO_tot_3, xmm_sum3), 6);

                        _mm_storeu_si128((__m128i *)(prediction_ptr), _mm_packus_epi16(xmm_sum0, xmm_sum1));
                        _mm_storeu_si128((__m128i *)(prediction_ptr + 16), _mm_packus_epi16(xmm_sum2, xmm_sum3));

                        offset -= 8;
                        xmm_LP = _mm_srli_si128(xmm_LP, 2);
                        xmm_BLP_tot = _mm_srli_si128(xmm_BLP_tot, 2);
                        prediction_ptr += prediction_buffer_stride;
                    }
                }
            }
            else {
                for (count = 0; count < 4; count += 2) {

                    xmm_BLP_tot = _mm_mullo_epi16(xmm_BLP, _mm_load_si128((__m128i *)(IntraPredictionConst - (count * 4) + OFFSET_1TO15)));
                    xmm_LP = _mm_srli_epi16(_mm_slli_epi16(_mm_loadu_si128((__m128i *)(ref_samples)), 8), 8);
                    ref_samples += 16;

                    for (count_in = 0; count_in < 8; ++count_in) {
                        xmm_LP_shift = xmm_LP;
                        xmm_BLP_tot_shift = xmm_BLP_tot;

                        xmm_LP_shift = _mm_unpacklo_epi16(xmm_LP_shift, xmm_LP_shift);
                        xmm_LP_shift = _mm_unpacklo_epi32(xmm_LP_shift, xmm_LP_shift);
                        xmm_LP_shift = _mm_unpacklo_epi64(xmm_LP_shift, xmm_LP_shift);
                        xmm_sum0 = _mm_add_epi16(xmm_TRP_tot_0, _mm_mullo_epi16(xmm_LP_shift, _mm_load_si128((__m128i *)(IntraPredictionConst + OFFSET_31TO24))));
                        xmm_sum1 = _mm_add_epi16(xmm_TRP_tot_1, _mm_mullo_epi16(xmm_LP_shift, _mm_load_si128((__m128i *)(IntraPredictionConst + OFFSET_23TO16))));
                        xmm_sum2 = _mm_add_epi16(xmm_TRP_tot_2, _mm_mullo_epi16(xmm_LP_shift, _mm_load_si128((__m128i *)(IntraPredictionConst + OFFSET_15TO8))));
                        xmm_sum3 = _mm_add_epi16(xmm_TRP_tot_3, _mm_mullo_epi16(xmm_LP_shift, _mm_load_si128((__m128i *)(IntraPredictionConst + OFFSET_7TO0))));

                        xmm_BLP_tot_shift = _mm_unpacklo_epi16(xmm_BLP_tot_shift, xmm_BLP_tot_shift);
                        xmm_BLP_tot_shift = _mm_unpacklo_epi32(xmm_BLP_tot_shift, xmm_BLP_tot_shift);
                        xmm_BLP_tot_shift = _mm_unpacklo_epi64(xmm_BLP_tot_shift, xmm_BLP_tot_shift);
                        xmm_sum0 = _mm_add_epi16(xmm_sum0, xmm_BLP_tot_shift);
                        xmm_sum1 = _mm_add_epi16(xmm_sum1, xmm_BLP_tot_shift);
                        xmm_sum2 = _mm_add_epi16(xmm_sum2, xmm_BLP_tot_shift);
                        xmm_sum3 = _mm_add_epi16(xmm_sum3, xmm_BLP_tot_shift);

                        xmm_TO_tot_0 = _mm_mullo_epi16(_mm_unpacklo_epi8(xmm_TO_lo, xmm0), _mm_load_si128((__m128i *)(IntraPredictionConst + offset)));
                        xmm_TO_tot_2 = _mm_mullo_epi16(_mm_unpacklo_epi8(xmm_TO_hi, xmm0), _mm_load_si128((__m128i *)(IntraPredictionConst + offset)));
                        xmm_TO_tot_1 = _mm_mullo_epi16(_mm_unpackhi_epi8(xmm_TO_lo, xmm0), _mm_load_si128((__m128i *)(IntraPredictionConst + offset)));
                        xmm_TO_tot_3 = _mm_mullo_epi16(_mm_unpackhi_epi8(xmm_TO_hi, xmm0), _mm_load_si128((__m128i *)(IntraPredictionConst + offset)));
                        xmm_sum0 = _mm_srli_epi16(_mm_add_epi16(xmm_TO_tot_0, xmm_sum0), 6);
                        xmm_sum1 = _mm_srli_epi16(_mm_add_epi16(xmm_TO_tot_1, xmm_sum1), 6);
                        xmm_sum2 = _mm_srli_epi16(_mm_add_epi16(xmm_TO_tot_2, xmm_sum2), 6);
                        xmm_sum3 = _mm_srli_epi16(_mm_add_epi16(xmm_TO_tot_3, xmm_sum3), 6);

                        _mm_storeu_si128((__m128i *)(prediction_ptr), _mm_packus_epi16(xmm_sum0, xmm_sum1));
                        _mm_storeu_si128((__m128i *)(prediction_ptr + 16), _mm_packus_epi16(xmm_sum2, xmm_sum3));

                        offset -= 16;
                        xmm_LP = _mm_srli_si128(xmm_LP, 2);
                        xmm_BLP_tot = _mm_srli_si128(xmm_BLP_tot, 2);
                        prediction_ptr += (prediction_buffer_stride << 1);
                    }
                }
            }
        }
        else if (size == 16) {
            uint32_t count_in;

            __m128i xmm_TRP_tot_lo, xmm_TRP_tot_hi, xmm_TO_lo, xmm_TO_hi, xmm_LP_Coeff_lo, xmm_LP_Coeff_hi;
            __m128i xmm_BLP_tot_lo, xmm_BLP_tot_hi, xmm_LP_lo, xmm_LP_hi, xmm_TO_tot_lo, xmm_TO_tot_hi;
            __m128i xmm_LP_tot_lo, xmm_LP_tot_hi, xmm_sum_lo, xmm_sum_hi;

            xmm_TRP_tot_lo = _mm_mullo_epi16(xmm_TRP, _mm_load_si128((__m128i *)(IntraPredictionConst + OFFSET_1TO8)));
            xmm_TRP_tot_hi = _mm_mullo_epi16(xmm_TRP, _mm_load_si128((__m128i *)(IntraPredictionConst + OFFSET_9TO16)));
            xmm_TRP_tot_lo = _mm_add_epi16(xmm_TRP_tot_lo, _mm_load_si128((__m128i *)(IntraPredictionConst + OFFSET_16)));
            xmm_TRP_tot_hi = _mm_add_epi16(xmm_TRP_tot_hi, _mm_load_si128((__m128i *)(IntraPredictionConst + OFFSET_16)));
            xmm_TO = _mm_loadu_si128((__m128i *)(ref_samples + topOffset));
            xmm_TO_lo = _mm_unpacklo_epi8(xmm_TO, xmm0);
            xmm_TO_hi = _mm_unpackhi_epi8(xmm_TO, xmm0);
            xmm_LP = _mm_loadu_si128((__m128i *)(ref_samples + leftOffset));
            xmm_LP_Coeff_lo = _mm_load_si128((__m128i *)(IntraPredictionConst + OFFSET_15TO8));
            xmm_LP_Coeff_hi = _mm_load_si128((__m128i *)(IntraPredictionConst + OFFSET_7TO0));

            if (!skip) {
                xmm_BLP_tot_lo = _mm_mullo_epi16(xmm_BLP, _mm_load_si128((__m128i *)(IntraPredictionConst + OFFSET_1TO8)));
                xmm_BLP_tot_hi = _mm_mullo_epi16(xmm_BLP, _mm_load_si128((__m128i *)(IntraPredictionConst + OFFSET_9TO16)));
                xmm_LP_lo = _mm_unpacklo_epi8(xmm_LP, xmm0);
                xmm_LP_hi = _mm_unpackhi_epi8(xmm_LP, xmm0);

                xmm_LP = xmm_LP_lo;
                xmm_BLP_tot = xmm_BLP_tot_lo;

                for (count = 0; count < 2; ++count) {
                    for (count_in = 0; count_in < 8; ++count_in) {
                        xmm_LP_shift = xmm_LP;
                        xmm_BLP_tot_shift = xmm_BLP_tot;

                        xmm_LP_shift = _mm_unpacklo_epi16(xmm_LP_shift, xmm_LP_shift);
                        xmm_LP_shift = _mm_unpacklo_epi32(xmm_LP_shift, xmm_LP_shift);
                        xmm_LP_shift = _mm_unpacklo_epi64(xmm_LP_shift, xmm_LP_shift);
                        xmm_LP_tot_lo = _mm_mullo_epi16(xmm_LP_shift, xmm_LP_Coeff_lo);
                        xmm_LP_tot_hi = _mm_mullo_epi16(xmm_LP_shift, xmm_LP_Coeff_hi);
                        xmm_sum_lo = _mm_add_epi16(xmm_LP_tot_lo, xmm_TRP_tot_lo);
                        xmm_sum_hi = _mm_add_epi16(xmm_LP_tot_hi, xmm_TRP_tot_hi);

                        xmm_BLP_tot_shift = _mm_unpacklo_epi16(xmm_BLP_tot_shift, xmm_BLP_tot_shift);
                        xmm_BLP_tot_shift = _mm_unpacklo_epi32(xmm_BLP_tot_shift, xmm_BLP_tot_shift);
                        xmm_BLP_tot_shift = _mm_unpacklo_epi64(xmm_BLP_tot_shift, xmm_BLP_tot_shift);
                        xmm_sum_lo = _mm_add_epi16(xmm_BLP_tot_shift, xmm_sum_lo);
                        xmm_sum_hi = _mm_add_epi16(xmm_BLP_tot_shift, xmm_sum_hi);

                        xmm_TO_tot_lo = _mm_mullo_epi16(xmm_TO_lo, _mm_load_si128((__m128i *)(IntraPredictionConst + OFFSET_15)));
                        xmm_TO_tot_hi = _mm_mullo_epi16(xmm_TO_hi, _mm_load_si128((__m128i *)(IntraPredictionConst + OFFSET_15)));
                        xmm_sum_lo = _mm_srli_epi16(_mm_add_epi16(xmm_sum_lo, xmm_TO_tot_lo), 5);
                        xmm_sum_hi = _mm_srli_epi16(_mm_add_epi16(xmm_sum_hi, xmm_TO_tot_hi), 5);

                        xmm_sum = _mm_packus_epi16(xmm_sum_lo, xmm_sum_hi);
                        _mm_storeu_si128((__m128i *)(prediction_ptr), xmm_sum);

                        xmm_LP = _mm_srli_si128(xmm_LP, 2);
                        xmm_BLP_tot = _mm_srli_si128(xmm_BLP_tot, 2);

                        IntraPredictionConst -= 8;
                        prediction_ptr += prediction_buffer_stride;
                    }
                    xmm_LP = xmm_LP_hi;
                    xmm_BLP_tot = xmm_BLP_tot_hi;
                }
            }
            else {
                xmm_BLP_tot = _mm_mullo_epi16(xmm_BLP, _mm_load_si128((__m128i *)(IntraPredictionConst + OFFSET_1TO15)));
                xmm_LP = _mm_srli_epi16(_mm_slli_epi16(xmm_LP, 8), 8);

                for (count = 0; count < 8; ++count) {
                    xmm_LP_shift = xmm_LP;
                    xmm_BLP_tot_shift = xmm_BLP_tot;

                    xmm_LP_shift = _mm_unpacklo_epi16(xmm_LP_shift, xmm_LP_shift);
                    xmm_LP_shift = _mm_unpacklo_epi32(xmm_LP_shift, xmm_LP_shift);
                    xmm_LP_shift = _mm_unpacklo_epi64(xmm_LP_shift, xmm_LP_shift);
                    xmm_LP_tot_lo = _mm_mullo_epi16(xmm_LP_shift, xmm_LP_Coeff_lo);
                    xmm_LP_tot_hi = _mm_mullo_epi16(xmm_LP_shift, xmm_LP_Coeff_hi);
                    xmm_sum_lo = _mm_add_epi16(xmm_LP_tot_lo, xmm_TRP_tot_lo);
                    xmm_sum_hi = _mm_add_epi16(xmm_LP_tot_hi, xmm_TRP_tot_hi);

                    xmm_BLP_tot_shift = _mm_unpacklo_epi16(xmm_BLP_tot_shift, xmm_BLP_tot_shift);
                    xmm_BLP_tot_shift = _mm_unpacklo_epi32(xmm_BLP_tot_shift, xmm_BLP_tot_shift);
                    xmm_BLP_tot_shift = _mm_unpacklo_epi64(xmm_BLP_tot_shift, xmm_BLP_tot_shift);
                    xmm_sum_lo = _mm_add_epi16(xmm_BLP_tot_shift, xmm_sum_lo);
                    xmm_sum_hi = _mm_add_epi16(xmm_BLP_tot_shift, xmm_sum_hi);

                    xmm_TO_tot_lo = _mm_mullo_epi16(xmm_TO_lo, _mm_load_si128((__m128i *)(IntraPredictionConst + OFFSET_15)));
                    xmm_TO_tot_hi = _mm_mullo_epi16(xmm_TO_hi, _mm_load_si128((__m128i *)(IntraPredictionConst + OFFSET_15)));
                    xmm_sum_lo = _mm_srli_epi16(_mm_add_epi16(xmm_sum_lo, xmm_TO_tot_lo), 5);
                    xmm_sum_hi = _mm_srli_epi16(_mm_add_epi16(xmm_sum_hi, xmm_TO_tot_hi), 5);
                    xmm_sum = _mm_packus_epi16(xmm_sum_lo, xmm_sum_hi);
                    _mm_storeu_si128((__m128i *)(prediction_ptr), xmm_sum);

                    xmm_LP = _mm_srli_si128(xmm_LP, 2);
                    xmm_BLP_tot = _mm_srli_si128(xmm_BLP_tot, 2);

                    IntraPredictionConst -= 16;
                    prediction_ptr += (prediction_buffer_stride << 1);
                }
            }
        }
        else {
            xmm_TRP_tot = _mm_mullo_epi16(xmm_TRP, _mm_load_si128((__m128i *)(IntraPredictionConst + OFFSET_1TO8)));
            xmm_TRP_tot = _mm_add_epi16(xmm_TRP_tot, _mm_load_si128((__m128i *)(IntraPredictionConst + OFFSET_8)));
            xmm_TO = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(ref_samples + topOffset)), xmm0);
            xmm_LP = _mm_loadl_epi64((__m128i *)(ref_samples + leftOffset));
            xmm_LP_Coeff = _mm_load_si128((__m128i *)(IntraPredictionConst + OFFSET_7TO0));

            if (!skip) {

                xmm_BLP_tot = _mm_mullo_epi16(xmm_BLP, _mm_load_si128((__m128i *)(IntraPredictionConst + OFFSET_1TO8)));
                xmm_LP = _mm_unpacklo_epi8(xmm_LP, xmm0);

                for (count = 0; count < 8; ++count) {
                    xmm_LP_shift = xmm_LP;
                    xmm_BLP_tot_shift = xmm_BLP_tot;

                    xmm_LP_shift = _mm_unpacklo_epi16(xmm_LP_shift, xmm_LP_shift);
                    xmm_LP_shift = _mm_unpacklo_epi32(xmm_LP_shift, xmm_LP_shift);
                    xmm_LP_shift = _mm_unpacklo_epi64(xmm_LP_shift, xmm_LP_shift);
                    xmm_LP_tot = _mm_mullo_epi16(xmm_LP_shift, xmm_LP_Coeff);
                    xmm_sum = _mm_add_epi16(xmm_LP_tot, xmm_TRP_tot);

                    xmm_BLP_tot_shift = _mm_unpacklo_epi16(xmm_BLP_tot_shift, xmm_BLP_tot_shift);
                    xmm_BLP_tot_shift = _mm_unpacklo_epi32(xmm_BLP_tot_shift, xmm_BLP_tot_shift);
                    xmm_BLP_tot_shift = _mm_unpacklo_epi64(xmm_BLP_tot_shift, xmm_BLP_tot_shift);
                    xmm_sum = _mm_add_epi16(xmm_BLP_tot_shift, xmm_sum);

                    xmm_TO_tot = _mm_mullo_epi16(xmm_TO, _mm_load_si128((__m128i *)(IntraPredictionConst + OFFSET_7)));
                    xmm_sum = _mm_srli_epi16(_mm_add_epi16(xmm_TO_tot, xmm_sum), 4);
                    xmm_sum = _mm_packus_epi16(xmm_sum, xmm_sum);
                    _mm_storel_epi64((__m128i *)(prediction_ptr), xmm_sum);

                    xmm_LP = _mm_srli_si128(xmm_LP, 2);
                    xmm_BLP_tot = _mm_srli_si128(xmm_BLP_tot, 2);

                    IntraPredictionConst -= 8;
                    prediction_ptr += prediction_buffer_stride;
                }
            }
            else {
                xmm_BLP_tot = _mm_mullo_epi16(xmm_BLP, _mm_load_si128((__m128i *)(IntraPredictionConst + OFFSET_1TO15)));
                xmm_LP = _mm_srli_epi16(_mm_slli_epi16(xmm_LP, 8), 8);

                for (count = 0; count < 4; ++count) {
                    xmm_LP_shift = xmm_LP;
                    xmm_BLP_tot_shift = xmm_BLP_tot;

                    xmm_LP_shift = _mm_unpacklo_epi16(xmm_LP_shift, xmm_LP_shift);
                    xmm_LP_shift = _mm_unpacklo_epi32(xmm_LP_shift, xmm_LP_shift);
                    xmm_LP_shift = _mm_unpacklo_epi64(xmm_LP_shift, xmm_LP_shift);
                    xmm_LP_tot = _mm_mullo_epi16(xmm_LP_shift, xmm_LP_Coeff);
                    xmm_sum = _mm_add_epi16(xmm_LP_tot, xmm_TRP_tot);

                    xmm_BLP_tot_shift = _mm_unpacklo_epi16(xmm_BLP_tot_shift, xmm_BLP_tot_shift);
                    xmm_BLP_tot_shift = _mm_unpacklo_epi32(xmm_BLP_tot_shift, xmm_BLP_tot_shift);
                    xmm_BLP_tot_shift = _mm_unpacklo_epi64(xmm_BLP_tot_shift, xmm_BLP_tot_shift);
                    xmm_sum = _mm_add_epi16(xmm_BLP_tot_shift, xmm_sum);

                    xmm_TO_tot = _mm_mullo_epi16(xmm_TO, _mm_load_si128((__m128i *)(IntraPredictionConst + OFFSET_7)));
                    xmm_sum = _mm_srli_epi16(_mm_add_epi16(xmm_TO_tot, xmm_sum), 4);
                    xmm_sum = _mm_packus_epi16(xmm_sum, xmm_sum);
                    _mm_storel_epi64((__m128i *)(prediction_ptr), xmm_sum);

                    xmm_LP = _mm_srli_si128(xmm_LP, 2);
                    xmm_BLP_tot = _mm_srli_si128(xmm_BLP_tot, 2);

                    IntraPredictionConst -= 16;
                    prediction_ptr += (prediction_buffer_stride << 1);

                }
            }
        }
    }
}
