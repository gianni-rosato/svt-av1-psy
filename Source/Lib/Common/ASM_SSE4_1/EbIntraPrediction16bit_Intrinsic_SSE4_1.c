/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbDefinitions.h"
#include "smmintrin.h"
#include "aom_dsp_rtcd.h"
EB_EXTERN void intra_mode_dc_luma16bit_sse4_1_intrin(
    const uint32_t   size,                       //input parameter, denotes the size of the current PU
    const uint16_t   *ref_samples,                 //input parameter, pointer to the reference samples
    uint16_t         *prediction_ptr,              //output parameter, pointer to the prediction
    const uint32_t   prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
    const EbBool  skip)                       //skip half rows
{
    uint32_t topOffset = (size << 1) + 1;
    uint32_t i;
    uint32_t pStride = skip ? (prediction_buffer_stride << 1) : prediction_buffer_stride;

    if (size == 4) {
        __m128i sum = _mm_setr_epi16(4, 0, 0, 0, 0, 0, 0, 0);
        __m128i leftT = _mm_loadl_epi64((__m128i *)ref_samples);
        __m128i topL = _mm_loadl_epi64((__m128i *)(ref_samples + topOffset));
        __m128i const2;
        __m128i temp0, temp1;
        sum = _mm_add_epi16(sum, topL);
        sum = _mm_add_epi16(sum, leftT);
        sum = _mm_hadd_epi16(sum, sum);
        sum = _mm_hadd_epi16(sum, sum);
        sum = _mm_srli_epi16(sum, 3);
        sum = _mm_unpacklo_epi16(sum, sum);
        sum = _mm_unpacklo_epi32(sum, sum);

        temp0 = _mm_add_epi16(sum, sum);      // 2*prediction_ptr[columnIndex]
        temp1 = _mm_add_epi16(leftT, topL);   // ref_samples[leftOffset] + ref_samples[topOffset]
        temp1 = _mm_add_epi16(temp1, temp0);  // ref_samples[leftOffset] + ref_samples[topOffset] + (prediction_ptr[0] << 1)
        temp0 = _mm_add_epi16(sum, temp0);    // 3*prediction_ptr[columnIndex]
        topL = _mm_add_epi16(topL, temp0);    // ref_samples[topOffset+columnIndex] + 3*prediction_ptr[columnIndex]
        leftT = _mm_add_epi16(leftT, temp0);  // ref_samples[leftOffset+row_index] + 3*prediction_ptr[writeIndex]
        const2 = _mm_setr_epi16(2, 2, 2, 2, 2, 2, 2, 2);
        topL = _mm_add_epi16(topL, const2);   // ref_samples[topOffset+columnIndex] + 3*prediction_ptr[columnIndex] + 2
        leftT = _mm_add_epi16(leftT, const2); // ref_samples[leftOffset+row_index] + 3*prediction_ptr[writeIndex] + 2
        temp1 = _mm_add_epi16(temp1, const2); // ref_samples[leftOffset] + ref_samples[topOffset] + (prediction_ptr[0] << 1) + 2
        topL = _mm_srli_epi16(topL, 2);       // (ref_samples[topOffset+columnIndex] + 3*prediction_ptr[columnIndex] + 2) >> 2
        leftT = _mm_srli_epi16(leftT, 2);     // (ref_samples[leftOffset+row_index] + 3*prediction_ptr[writeIndex] + 2) >> 2
        temp1 = _mm_srli_epi16(temp1, 2);     // (ref_samples[leftOffset] + ref_samples[topOffset] + (prediction_ptr[0] << 1) + 2) >> 2

        topL = _mm_blend_epi16(topL, temp1, 1); // prediction_ptr[0] = (uint16_t) ((ref_samples[leftOffset] + ref_samples[topOffset] + (prediction_ptr[0] << 1) + 2) >> 2);
        _mm_storel_epi64((__m128i *)prediction_ptr, topL); // prediction_ptr[columnIndex] = (uint16_t) ((ref_samples[topOffset+columnIndex] + 3*prediction_ptr[columnIndex] + 2) >> 2);

        if (skip) {
            leftT = _mm_srli_si128(leftT, 4);
            topL = _mm_blend_epi16(sum, leftT, 1);
            _mm_storel_epi64((__m128i *)(prediction_ptr + pStride), topL);
        }
        else {
            leftT = _mm_srli_si128(leftT, 2);
            topL = _mm_blend_epi16(sum, leftT, 1);
            _mm_storel_epi64((__m128i *)(prediction_ptr + pStride), topL);

            leftT = _mm_srli_si128(leftT, 2);
            topL = _mm_blend_epi16(sum, leftT, 1);
            prediction_ptr += 2 * pStride;
            _mm_storel_epi64((__m128i *)prediction_ptr, topL);
            leftT = _mm_srli_si128(leftT, 2);
            topL = _mm_blend_epi16(sum, leftT, 1);
            _mm_storel_epi64((__m128i *)(prediction_ptr + pStride), topL);
        }
    }
    else if (size == 8) {
        __m128i sum = _mm_setr_epi16(8, 0, 0, 0, 0, 0, 0, 0);
        __m128i leftT = _mm_loadu_si128((__m128i *)ref_samples);
        __m128i topL = _mm_loadu_si128((__m128i *)(ref_samples + topOffset));
        __m128i const2;
        __m128i temp0, temp1;
        sum = _mm_add_epi16(sum, topL);
        sum = _mm_add_epi16(sum, leftT);
        sum = _mm_hadd_epi16(sum, sum);
        sum = _mm_hadd_epi16(sum, sum);
        sum = _mm_hadd_epi16(sum, sum);
        sum = _mm_srli_epi16(sum, 4);
        sum = _mm_unpacklo_epi16(sum, sum);
        sum = _mm_unpacklo_epi32(sum, sum);
        sum = _mm_unpacklo_epi64(sum, sum);

        temp0 = _mm_add_epi16(sum, sum);      // 2*prediction_ptr[columnIndex]
        temp1 = _mm_add_epi16(leftT, topL);   // ref_samples[leftOffset] + ref_samples[topOffset]
        temp1 = _mm_add_epi16(temp1, temp0);  // ref_samples[leftOffset] + ref_samples[topOffset] + (prediction_ptr[0] << 1)
        temp0 = _mm_add_epi16(sum, temp0);    // 3*prediction_ptr[columnIndex]
        topL = _mm_add_epi16(topL, temp0);    // ref_samples[topOffset+columnIndex] + 3*prediction_ptr[columnIndex]
        leftT = _mm_add_epi16(leftT, temp0);  // ref_samples[leftOffset+row_index] + 3*prediction_ptr[writeIndex]
        const2 = _mm_setr_epi16(2, 2, 2, 2, 2, 2, 2, 2);
        topL = _mm_add_epi16(topL, const2);   // ref_samples[topOffset+columnIndex] + 3*prediction_ptr[columnIndex] + 2
        leftT = _mm_add_epi16(leftT, const2); // ref_samples[leftOffset+row_index] + 3*prediction_ptr[writeIndex] + 2
        temp1 = _mm_add_epi16(temp1, const2); // ref_samples[leftOffset] + ref_samples[topOffset] + (prediction_ptr[0] << 1) + 2
        topL = _mm_srli_epi16(topL, 2);       // (ref_samples[topOffset+columnIndex] + 3*prediction_ptr[columnIndex] + 2) >> 2
        leftT = _mm_srli_epi16(leftT, 2);     // (ref_samples[leftOffset+row_index] + 3*prediction_ptr[writeIndex] + 2) >> 2
        temp1 = _mm_srli_epi16(temp1, 2);     // (ref_samples[leftOffset] + ref_samples[topOffset] + (prediction_ptr[0] << 1) + 2) >> 2

        topL = _mm_blend_epi16(topL, temp1, 1); // prediction_ptr[0] = (uint16_t) ((ref_samples[leftOffset] + ref_samples[topOffset] + (prediction_ptr[0] << 1) + 2) >> 2);
        _mm_storeu_si128((__m128i *)prediction_ptr, topL); // prediction_ptr[columnIndex] = (uint16_t) ((ref_samples[topOffset+columnIndex] + 3*prediction_ptr[columnIndex] + 2) >> 2);

        if (skip) {
            leftT = _mm_srli_si128(leftT, 4);
            topL = _mm_blend_epi16(sum, leftT, 1);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), topL);

            leftT = _mm_srli_si128(leftT, 4);
            topL = _mm_blend_epi16(sum, leftT, 1);
            prediction_ptr += 2 * pStride;
            _mm_storeu_si128((__m128i *)prediction_ptr, topL);
            leftT = _mm_srli_si128(leftT, 4);
            topL = _mm_blend_epi16(sum, leftT, 1);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), topL);
        }
        else {
            leftT = _mm_srli_si128(leftT, 2);
            topL = _mm_blend_epi16(sum, leftT, 1);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), topL);

            leftT = _mm_srli_si128(leftT, 2);
            topL = _mm_blend_epi16(sum, leftT, 1);
            prediction_ptr += 2 * pStride;
            _mm_storeu_si128((__m128i *)prediction_ptr, topL);
            leftT = _mm_srli_si128(leftT, 2);
            topL = _mm_blend_epi16(sum, leftT, 1);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), topL);

            leftT = _mm_srli_si128(leftT, 2);
            topL = _mm_blend_epi16(sum, leftT, 1);
            prediction_ptr += 2 * pStride;
            _mm_storeu_si128((__m128i *)prediction_ptr, topL);
            leftT = _mm_srli_si128(leftT, 2);
            topL = _mm_blend_epi16(sum, leftT, 1);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), topL);

            leftT = _mm_srli_si128(leftT, 2);
            topL = _mm_blend_epi16(sum, leftT, 1);
            prediction_ptr += 2 * pStride;
            _mm_storeu_si128((__m128i *)prediction_ptr, topL);
            leftT = _mm_srli_si128(leftT, 2);
            topL = _mm_blend_epi16(sum, leftT, 1);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), topL);
        }
    }
    else if (size == 16) {
        __m128i sum = _mm_setr_epi16(16, 0, 0, 0, 0, 0, 0, 0);
        __m128i leftT = _mm_loadu_si128((__m128i *)ref_samples);
        __m128i leftB = _mm_loadu_si128((__m128i *)(ref_samples + 8));
        __m128i topL = _mm_loadu_si128((__m128i *)(ref_samples + topOffset));
        __m128i topR = _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 8));
        __m128i const2;
        __m128i temp0, temp1;
        sum = _mm_add_epi16(sum, topL);
        sum = _mm_add_epi16(sum, topR);
        sum = _mm_add_epi16(sum, leftT);
        sum = _mm_add_epi16(sum, leftB);
        sum = _mm_hadd_epi16(sum, sum);
        sum = _mm_hadd_epi16(sum, sum);
        sum = _mm_hadd_epi16(sum, sum);
        sum = _mm_srli_epi16(sum, 5);
        sum = _mm_unpacklo_epi16(sum, sum);
        sum = _mm_unpacklo_epi32(sum, sum);
        sum = _mm_unpacklo_epi64(sum, sum);

        temp0 = _mm_add_epi16(sum, sum);      // 2*prediction_ptr[columnIndex]
        temp1 = _mm_add_epi16(leftT, topL);   // ref_samples[leftOffset] + ref_samples[topOffset]
        temp1 = _mm_add_epi16(temp1, temp0);  // ref_samples[leftOffset] + ref_samples[topOffset] + (prediction_ptr[0] << 1)
        temp0 = _mm_add_epi16(sum, temp0);    // 3*prediction_ptr[columnIndex]
        topL = _mm_add_epi16(topL, temp0);    // ref_samples[topOffset+columnIndex] + 3*prediction_ptr[columnIndex]
        topR = _mm_add_epi16(topR, temp0);    // ref_samples[topOffset+columnIndex] + 3*prediction_ptr[columnIndex]
        leftT = _mm_add_epi16(leftT, temp0);  // ref_samples[leftOffset+row_index] + 3*prediction_ptr[writeIndex]
        leftB = _mm_add_epi16(leftB, temp0);  // ref_samples[leftOffset+row_index] + 3*prediction_ptr[writeIndex]
        const2 = _mm_setr_epi16(2, 2, 2, 2, 2, 2, 2, 2);
        topL = _mm_add_epi16(topL, const2);   // ref_samples[topOffset+columnIndex] + 3*prediction_ptr[columnIndex] + 2
        topR = _mm_add_epi16(topR, const2);   // ref_samples[topOffset+columnIndex] + 3*prediction_ptr[columnIndex] + 2
        leftT = _mm_add_epi16(leftT, const2); // ref_samples[leftOffset+row_index] + 3*prediction_ptr[writeIndex] + 2
        leftB = _mm_add_epi16(leftB, const2); // ref_samples[leftOffset+row_index] + 3*prediction_ptr[writeIndex] + 2
        temp1 = _mm_add_epi16(temp1, const2); // ref_samples[leftOffset] + ref_samples[topOffset] + (prediction_ptr[0] << 1) + 2
        topL = _mm_srli_epi16(topL, 2);       // (ref_samples[topOffset+columnIndex] + 3*prediction_ptr[columnIndex] + 2) >> 2
        topR = _mm_srli_epi16(topR, 2);       // (ref_samples[topOffset+columnIndex] + 3*prediction_ptr[columnIndex] + 2) >> 2
        leftT = _mm_srli_epi16(leftT, 2);     // (ref_samples[leftOffset+row_index] + 3*prediction_ptr[writeIndex] + 2) >> 2
        leftB = _mm_srli_epi16(leftB, 2);     // (ref_samples[leftOffset+row_index] + 3*prediction_ptr[writeIndex] + 2) >> 2
        temp1 = _mm_srli_epi16(temp1, 2);     // (ref_samples[leftOffset] + ref_samples[topOffset] + (prediction_ptr[0] << 1) + 2) >> 2

        topL = _mm_blend_epi16(topL, temp1, 1); // prediction_ptr[0] = (uint16_t) ((ref_samples[leftOffset] + ref_samples[topOffset] + (prediction_ptr[0] << 1) + 2) >> 2);
        _mm_storeu_si128((__m128i *)prediction_ptr, topL); // prediction_ptr[columnIndex] = (uint16_t) ((ref_samples[topOffset+columnIndex] + 3*prediction_ptr[columnIndex] + 2) >> 2);
        _mm_storeu_si128((__m128i *)(prediction_ptr + 8), topR);

        if (skip) {
            leftT = _mm_srli_si128(leftT, 4);
            topL = _mm_blend_epi16(sum, leftT, 1);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), topL);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride + 8), sum);

            leftT = _mm_srli_si128(leftT, 4);
            topL = _mm_blend_epi16(sum, leftT, 1);
            prediction_ptr += 2 * pStride;
            _mm_storeu_si128((__m128i *)prediction_ptr, topL);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 8), sum);
            leftT = _mm_srli_si128(leftT, 4);
            topL = _mm_blend_epi16(sum, leftT, 1);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), topL);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride + 8), sum);

            topL = _mm_blend_epi16(sum, leftB, 1);
            prediction_ptr += 2 * pStride;
            _mm_storeu_si128((__m128i *)prediction_ptr, topL);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 8), sum);
            leftB = _mm_srli_si128(leftB, 4);
            topL = _mm_blend_epi16(sum, leftB, 1);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), topL);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride + 8), sum);

            leftB = _mm_srli_si128(leftB, 4);
            topL = _mm_blend_epi16(sum, leftB, 1);
            prediction_ptr += 2 * pStride;
            _mm_storeu_si128((__m128i *)prediction_ptr, topL);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 8), sum);
            leftB = _mm_srli_si128(leftB, 4);
            topL = _mm_blend_epi16(sum, leftB, 1);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), topL);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride + 8), sum);
        }
        else {
            leftT = _mm_srli_si128(leftT, 2);
            topL = _mm_blend_epi16(sum, leftT, 1);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), topL);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride + 8), sum);

            leftT = _mm_srli_si128(leftT, 2);
            topL = _mm_blend_epi16(sum, leftT, 1);
            prediction_ptr += 2 * pStride;
            _mm_storeu_si128((__m128i *)prediction_ptr, topL);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 8), sum);
            leftT = _mm_srli_si128(leftT, 2);
            topL = _mm_blend_epi16(sum, leftT, 1);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), topL);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride + 8), sum);

            leftT = _mm_srli_si128(leftT, 2);
            topL = _mm_blend_epi16(sum, leftT, 1);
            prediction_ptr += 2 * pStride;
            _mm_storeu_si128((__m128i *)prediction_ptr, topL);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 8), sum);
            leftT = _mm_srli_si128(leftT, 2);
            topL = _mm_blend_epi16(sum, leftT, 1);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), topL);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride + 8), sum);

            leftT = _mm_srli_si128(leftT, 2);
            topL = _mm_blend_epi16(sum, leftT, 1);
            prediction_ptr += 2 * pStride;
            _mm_storeu_si128((__m128i *)prediction_ptr, topL);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 8), sum);
            leftT = _mm_srli_si128(leftT, 2);
            topL = _mm_blend_epi16(sum, leftT, 1);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), topL);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride + 8), sum);

            topL = _mm_blend_epi16(sum, leftB, 1);
            prediction_ptr += 2 * pStride;
            _mm_storeu_si128((__m128i *)prediction_ptr, topL);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 8), sum);
            leftB = _mm_srli_si128(leftB, 2);
            topL = _mm_blend_epi16(sum, leftB, 1);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), topL);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride + 8), sum);

            leftB = _mm_srli_si128(leftB, 2);
            topL = _mm_blend_epi16(sum, leftB, 1);
            prediction_ptr += 2 * pStride;
            _mm_storeu_si128((__m128i *)prediction_ptr, topL);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 8), sum);
            leftB = _mm_srli_si128(leftB, 2);
            topL = _mm_blend_epi16(sum, leftB, 1);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), topL);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride + 8), sum);

            leftB = _mm_srli_si128(leftB, 2);
            topL = _mm_blend_epi16(sum, leftB, 1);
            prediction_ptr += 2 * pStride;
            _mm_storeu_si128((__m128i *)prediction_ptr, topL);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 8), sum);
            leftB = _mm_srli_si128(leftB, 2);
            topL = _mm_blend_epi16(sum, leftB, 1);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), topL);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride + 8), sum);

            leftB = _mm_srli_si128(leftB, 2);
            topL = _mm_blend_epi16(sum, leftB, 1);
            prediction_ptr += 2 * pStride;
            _mm_storeu_si128((__m128i *)prediction_ptr, topL);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 8), sum);
            leftB = _mm_srli_si128(leftB, 2);
            topL = _mm_blend_epi16(sum, leftB, 1);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), topL);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride + 8), sum);
        }
    }
    else { /*if (size == 32) {*/
        __m128i sum = _mm_setr_epi16(32, 0, 0, 0, 0, 0, 0, 0);
        sum = _mm_add_epi16(sum, _mm_loadu_si128((__m128i *)ref_samples));
        sum = _mm_add_epi16(sum, _mm_loadu_si128((__m128i *)(ref_samples + 8)));
        sum = _mm_add_epi16(sum, _mm_loadu_si128((__m128i *)(ref_samples + 16)));
        sum = _mm_add_epi16(sum, _mm_loadu_si128((__m128i *)(ref_samples + 24)));
        sum = _mm_add_epi16(sum, _mm_loadu_si128((__m128i *)(ref_samples + topOffset)));
        sum = _mm_add_epi16(sum, _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 8)));
        sum = _mm_add_epi16(sum, _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 16)));
        sum = _mm_add_epi16(sum, _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 24)));
        sum = _mm_hadd_epi16(sum, sum);
        sum = _mm_hadd_epi16(sum, sum);
        sum = _mm_hadd_epi16(sum, sum);
        sum = _mm_srli_epi16(sum, 6);
        sum = _mm_unpacklo_epi16(sum, sum);
        sum = _mm_unpacklo_epi32(sum, sum);
        sum = _mm_unpacklo_epi64(sum, sum);
        for (i = 0; i < (uint32_t)(skip ? 8 : 16); i++) {
            _mm_storeu_si128((__m128i *)prediction_ptr, sum);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 8), sum);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 16), sum);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 24), sum);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), sum);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride + 8), sum);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride + 16), sum);
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride + 24), sum);
            prediction_ptr += 2 * pStride;
        }
    }
}

void av1_filter_intra_edge_sse4_1(uint8_t *p, int32_t sz, int32_t strength) {
    if (!strength) return;

    DECLARE_ALIGNED(16, static const int8_t, kern[3][16]) = {
        { 4, 8, 4, 0, 4, 8, 4, 0, 4, 8, 4, 0, 4, 8, 4, 0 },  // strength 1: 4,8,4
        { 5, 6, 5, 0, 5, 6, 5, 0, 5, 6, 5, 0, 5, 6, 5, 0 },  // strength 2: 5,6,5
        { 2, 4, 4, 4, 2, 0, 0, 0, 2, 4, 4, 4, 2, 0, 0, 0 }  // strength 3: 2,4,4,4,2
    };

    DECLARE_ALIGNED(16, static const int8_t, v_const[5][16]) = {
        { 0, 1, 2, 3, 1, 2, 3, 4, 2, 3, 4, 5, 3, 4, 5, 6 },
        { 4, 5, 6, 7, 5, 6, 7, 8, 6, 7, 8, 9, 7, 8, 9, 10 },
        { 0, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 8 },
        { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 },
    };

    // Extend the first and last samples to simplify the loop for the 5-tap case
    p[-1] = p[0];
    __m128i last = _mm_set1_epi8(p[sz - 1]);
    _mm_storeu_si128((__m128i *)&p[sz], last);

    // Adjust input pointer for filter support area
    uint8_t *in = (strength == 3) ? p - 1 : p;

    // Avoid modifying first sample
    uint8_t *out = p + 1;
    int32_t len = sz - 1;

    const int32_t use_3tap_filter = (strength < 3);

    if (use_3tap_filter) {
        __m128i coef0 = _mm_lddqu_si128((__m128i const *)kern[strength - 1]);
        __m128i shuf0 = _mm_lddqu_si128((__m128i const *)v_const[0]);
        __m128i shuf1 = _mm_lddqu_si128((__m128i const *)v_const[1]);
        __m128i iden = _mm_lddqu_si128((__m128i *)v_const[3]);
        __m128i in0 = _mm_lddqu_si128((__m128i *)in);
        while (len > 0) {
            int32_t n_out = (len < 8) ? len : 8;
            __m128i d0 = _mm_shuffle_epi8(in0, shuf0);
            __m128i d1 = _mm_shuffle_epi8(in0, shuf1);
            d0 = _mm_maddubs_epi16(d0, coef0);
            d1 = _mm_maddubs_epi16(d1, coef0);
            d0 = _mm_hadd_epi16(d0, d1);
            __m128i eight = _mm_set1_epi16(8);
            d0 = _mm_add_epi16(d0, eight);
            d0 = _mm_srai_epi16(d0, 4);
            d0 = _mm_packus_epi16(d0, d0);
            __m128i out0 = _mm_lddqu_si128((__m128i *)out);
            __m128i n0 = _mm_set1_epi8(n_out);
            __m128i mask = _mm_cmpgt_epi8(n0, iden);
            out0 = _mm_blendv_epi8(out0, d0, mask);
            _mm_storel_epi64((__m128i *)out, out0);
            __m128i in1 = _mm_lddqu_si128((__m128i *)(in + 16));
            in0 = _mm_alignr_epi8(in1, in0, 8);
            in += 8;
            out += 8;
            len -= n_out;
        }
    }
    else {  // 5-tap filter
        __m128i coef0 = _mm_lddqu_si128((__m128i const *)kern[strength - 1]);
        __m128i two = _mm_set1_epi8(2);
        __m128i shuf_a = _mm_lddqu_si128((__m128i const *)v_const[2]);
        __m128i shuf_b = _mm_add_epi8(shuf_a, two);
        __m128i shuf_c = _mm_add_epi8(shuf_b, two);
        __m128i shuf_d = _mm_add_epi8(shuf_c, two);
        __m128i iden = _mm_lddqu_si128((__m128i *)v_const[3]);
        __m128i in0 = _mm_lddqu_si128((__m128i *)in);
        while (len > 0) {
            int32_t n_out = (len < 8) ? len : 8;
            __m128i d0 = _mm_shuffle_epi8(in0, shuf_a);
            __m128i d1 = _mm_shuffle_epi8(in0, shuf_b);
            __m128i d2 = _mm_shuffle_epi8(in0, shuf_c);
            __m128i d3 = _mm_shuffle_epi8(in0, shuf_d);
            d0 = _mm_maddubs_epi16(d0, coef0);
            d1 = _mm_maddubs_epi16(d1, coef0);
            d2 = _mm_maddubs_epi16(d2, coef0);
            d3 = _mm_maddubs_epi16(d3, coef0);
            d0 = _mm_hadd_epi16(d0, d1);
            d2 = _mm_hadd_epi16(d2, d3);
            d0 = _mm_hadd_epi16(d0, d2);
            __m128i eight = _mm_set1_epi16(8);
            d0 = _mm_add_epi16(d0, eight);
            d0 = _mm_srai_epi16(d0, 4);
            d0 = _mm_packus_epi16(d0, d0);
            __m128i out0 = _mm_lddqu_si128((__m128i *)out);
            __m128i n0 = _mm_set1_epi8(n_out);
            __m128i mask = _mm_cmpgt_epi8(n0, iden);
            out0 = _mm_blendv_epi8(out0, d0, mask);
            _mm_storel_epi64((__m128i *)out, out0);
            __m128i in1 = _mm_lddqu_si128((__m128i *)(in + 16));
            in0 = _mm_alignr_epi8(in1, in0, 8);
            in += 8;
            out += 8;
            len -= n_out;
        }
    }
}

void av1_filter_intra_edge_high_sse4_1(uint16_t *p, int32_t sz, int32_t strength) {
    if (!strength) return;

    DECLARE_ALIGNED(16, static const int16_t, kern[3][8]) = {
        { 4, 8, 4, 8, 4, 8, 4, 8 },  // strength 1: 4,8,4
        { 5, 6, 5, 6, 5, 6, 5, 6 },  // strength 2: 5,6,5
        { 2, 4, 2, 4, 2, 4, 2, 4 }   // strength 3: 2,4,4,4,2
    };

    DECLARE_ALIGNED(16, static const int16_t,
    v_const[1][8]) = { { 0, 1, 2, 3, 4, 5, 6, 7 } };

    // Extend the first and last samples to simplify the loop for the 5-tap case
    p[-1] = p[0];
    __m128i last = _mm_set1_epi16(p[sz - 1]);
    _mm_storeu_si128((__m128i *)&p[sz], last);

    // Adjust input pointer for filter support area
    uint16_t *in = (strength == 3) ? p - 1 : p;

    // Avoid modifying first sample
    uint16_t *out = p + 1;
    int32_t len = sz - 1;

    const int32_t use_3tap_filter = (strength < 3);

    if (use_3tap_filter) {
        __m128i coef0 = _mm_lddqu_si128((__m128i const *)kern[strength - 1]);
        __m128i iden = _mm_lddqu_si128((__m128i *)v_const[0]);
        __m128i in0 = _mm_lddqu_si128((__m128i *)&in[0]);
        __m128i in8 = _mm_lddqu_si128((__m128i *)&in[8]);
        while (len > 0) {
            int32_t n_out = (len < 8) ? len : 8;
            __m128i in1 = _mm_alignr_epi8(in8, in0, 2);
            __m128i in2 = _mm_alignr_epi8(in8, in0, 4);
            __m128i in02 = _mm_add_epi16(in0, in2);
            __m128i d0 = _mm_unpacklo_epi16(in02, in1);
            __m128i d1 = _mm_unpackhi_epi16(in02, in1);
            d0 = _mm_mullo_epi16(d0, coef0);
            d1 = _mm_mullo_epi16(d1, coef0);
            d0 = _mm_hadd_epi16(d0, d1);
            __m128i eight = _mm_set1_epi16(8);
            d0 = _mm_add_epi16(d0, eight);
            d0 = _mm_srli_epi16(d0, 4);
            __m128i out0 = _mm_lddqu_si128((__m128i *)out);
            __m128i n0 = _mm_set1_epi16(n_out);
            __m128i mask = _mm_cmpgt_epi16(n0, iden);
            out0 = _mm_blendv_epi8(out0, d0, mask);
            _mm_storeu_si128((__m128i *)out, out0);
            in += 8;
            in0 = in8;
            in8 = _mm_lddqu_si128((__m128i *)&in[8]);
            out += 8;
            len -= n_out;
        }
    }
    else {  // 5-tap filter
        __m128i coef0 = _mm_lddqu_si128((__m128i const *)kern[strength - 1]);
        __m128i iden = _mm_lddqu_si128((__m128i *)v_const[0]);
        __m128i in0 = _mm_lddqu_si128((__m128i *)&in[0]);
        __m128i in8 = _mm_lddqu_si128((__m128i *)&in[8]);
        while (len > 0) {
            int32_t n_out = (len < 8) ? len : 8;
            __m128i in1 = _mm_alignr_epi8(in8, in0, 2);
            __m128i in2 = _mm_alignr_epi8(in8, in0, 4);
            __m128i in3 = _mm_alignr_epi8(in8, in0, 6);
            __m128i in4 = _mm_alignr_epi8(in8, in0, 8);
            __m128i in04 = _mm_add_epi16(in0, in4);
            __m128i in123 = _mm_add_epi16(in1, in2);
            in123 = _mm_add_epi16(in123, in3);
            __m128i d0 = _mm_unpacklo_epi16(in04, in123);
            __m128i d1 = _mm_unpackhi_epi16(in04, in123);
            d0 = _mm_mullo_epi16(d0, coef0);
            d1 = _mm_mullo_epi16(d1, coef0);
            d0 = _mm_hadd_epi16(d0, d1);
            __m128i eight = _mm_set1_epi16(8);
            d0 = _mm_add_epi16(d0, eight);
            d0 = _mm_srli_epi16(d0, 4);
            __m128i out0 = _mm_lddqu_si128((__m128i *)out);
            __m128i n0 = _mm_set1_epi16(n_out);
            __m128i mask = _mm_cmpgt_epi16(n0, iden);
            out0 = _mm_blendv_epi8(out0, d0, mask);
            _mm_storeu_si128((__m128i *)out, out0);
            in += 8;
            in0 = in8;
            in8 = _mm_lddqu_si128((__m128i *)&in[8]);
            out += 8;
            len -= n_out;
        }
    }
}

void av1_upsample_intra_edge_sse4_1(uint8_t *p, int32_t sz) {
    // interpolate half-sample positions
    assert(sz <= 24);

    DECLARE_ALIGNED(16, static const int8_t, kernel[1][16]) = {
        { -1, 9, 9, -1, -1, 9, 9, -1, -1, 9, 9, -1, -1, 9, 9, -1 }
    };

    DECLARE_ALIGNED(16, static const int8_t, v_const[2][16]) = {
        { 0, 1, 2, 3, 1, 2, 3, 4, 2, 3, 4, 5, 3, 4, 5, 6 },
        { 4, 5, 6, 7, 5, 6, 7, 8, 6, 7, 8, 9, 7, 8, 9, 10 }
    };

    // Extend first/last samples (upper-left p[-1], last p[sz-1])
    // to support 4-tap filter
    p[-2] = p[-1];
    p[sz] = p[sz - 1];

    uint8_t *in = &p[-2];
    uint8_t *out = &p[-2];

    int32_t n = sz + 1;  // Input length including upper-left sample

    __m128i in0 = _mm_lddqu_si128((__m128i *)&in[0]);
    __m128i in16 = _mm_lddqu_si128((__m128i *)&in[16]);

    __m128i coef0 = _mm_lddqu_si128((__m128i *)kernel[0]);
    __m128i shuf0 = _mm_lddqu_si128((__m128i *)v_const[0]);
    __m128i shuf1 = _mm_lddqu_si128((__m128i *)v_const[1]);

    while (n > 0) {
        __m128i in8 = _mm_alignr_epi8(in16, in0, 8);
        __m128i d0 = _mm_shuffle_epi8(in0, shuf0);
        __m128i d1 = _mm_shuffle_epi8(in0, shuf1);
        __m128i d2 = _mm_shuffle_epi8(in8, shuf0);
        __m128i d3 = _mm_shuffle_epi8(in8, shuf1);
        d0 = _mm_maddubs_epi16(d0, coef0);
        d1 = _mm_maddubs_epi16(d1, coef0);
        d2 = _mm_maddubs_epi16(d2, coef0);
        d3 = _mm_maddubs_epi16(d3, coef0);
        d0 = _mm_hadd_epi16(d0, d1);
        d2 = _mm_hadd_epi16(d2, d3);
        __m128i eight = _mm_set1_epi16(8);
        d0 = _mm_add_epi16(d0, eight);
        d2 = _mm_add_epi16(d2, eight);
        d0 = _mm_srai_epi16(d0, 4);
        d2 = _mm_srai_epi16(d2, 4);
        d0 = _mm_packus_epi16(d0, d2);
        __m128i in1 = _mm_alignr_epi8(in16, in0, 1);
        __m128i out0 = _mm_unpacklo_epi8(in1, d0);
        __m128i out1 = _mm_unpackhi_epi8(in1, d0);
        _mm_storeu_si128((__m128i *)&out[0], out0);
        _mm_storeu_si128((__m128i *)&out[16], out1);
        in0 = in16;
        in16 = _mm_setzero_si128();
        out += 32;
        n -= 16;
    }
}