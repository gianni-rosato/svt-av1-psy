/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbDefinitions.h"
#include "tmmintrin.h"
#include "EbIntraPrediction_SSSE3.h"

static void _mm_storeh_epi64(__m128i * p, __m128i x)
{
    _mm_storeh_pd((double *)p, _mm_castsi128_pd(x));
}

extern void intra_mode_angular_vertical_kernel_ssse3_intrin(
    uint32_t         size,                       //input parameter, denotes the size of the current PU
    uint8_t         *ref_samp_main,                //input parameter, pointer to the reference samples
    uint8_t         *prediction_ptr,              //output parameter, pointer to the prediction
    uint32_t   prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
    const EbBool  skip,
    int32_t   intra_pred_angle)
{
    uint32_t row_index;
    uint32_t height = size;
    int32_t deltaSum = intra_pred_angle;
    int32_t deltaInt;
    uint32_t deltaFract;
    __m128i top0, top1, top2, sum0, sum1, a0, a1;

    // --------- Reference Samples Structure ---------
    // ref_samp_main[-size+1] to ref_samp_main[-1] must be prepared (from bottom to top) for mode 19 to 25 (not required for mode 27 to 33)
    // ref_samp_main[0]      = TopLeft[0]
    // ref_samp_main[1]      = Top[0]
    // ref_samp_main[2]      = Top[1]
    // ...
    // ref_samp_main[size]   = Top[size-1]
    // ref_samp_main[size+1] = Top[size]     for mode 27 to 33 (not required for mode 19 to 25)
    // ...
    // ref_samp_main[2*size] = Top[2*size-1] for mode 27 to 33 (not required for mode 19 to 25)
    // -----------------------------------------------

    // Compute the prediction
    ref_samp_main += 1; // top0 sample
    if (skip) {
        height >>= 1;
        prediction_buffer_stride <<= 1;
        intra_pred_angle <<= 1;
    }

    if (size == 4) {
        for (row_index = 0; row_index < height; row_index += 2) {
            deltaInt = deltaSum >> 5;
            deltaFract = deltaSum & 31;
            a0 = _mm_unpacklo_epi8(_mm_cvtsi32_si128(32 - deltaFract), _mm_cvtsi32_si128(deltaFract));
            top0 = _mm_loadl_epi64((__m128i *)(ref_samp_main + deltaInt));
            deltaSum += intra_pred_angle;
            deltaInt = deltaSum >> 5;
            deltaFract = deltaSum & 31;
            a1 = _mm_unpacklo_epi8(_mm_cvtsi32_si128(32 - deltaFract), _mm_cvtsi32_si128(deltaFract));
            a0 = _mm_unpacklo_epi16(a0, a1);
            a0 = _mm_shuffle_epi8(a0, _mm_setr_epi8(0, 1, 0, 1, 0, 1, 0, 1, 2, 3, 2, 3, 2, 3, 2, 3));
            top0 = _mm_castps_si128(_mm_loadh_pi(_mm_castsi128_ps(top0), (__m64 *)(ref_samp_main + deltaInt)));
            top0 = _mm_shuffle_epi8(top0, _mm_setr_epi8(0, 1, 1, 2, 2, 3, 3, 4, 8, 9, 9, 10, 10, 11, 11, 12));
            sum0 = _mm_add_epi16(_mm_set1_epi16(16), _mm_maddubs_epi16(top0, a0));
            sum0 = _mm_srai_epi16(sum0, 5);
            sum0 = _mm_packus_epi16(sum0, sum0);
            *(uint32_t *)prediction_ptr = _mm_cvtsi128_si32(sum0);
            *(uint32_t *)(prediction_ptr + prediction_buffer_stride) = _mm_cvtsi128_si32(_mm_srli_si128(sum0, 4));
            prediction_ptr += 2 * prediction_buffer_stride;
            deltaSum += intra_pred_angle;
        }
    }
    else if (size == 8) {
        for (row_index = 0; row_index < height; row_index++) {
            deltaInt = deltaSum >> 5;
            deltaFract = deltaSum & 31;
            a0 = _mm_unpacklo_epi8(_mm_cvtsi32_si128(32 - deltaFract), _mm_cvtsi32_si128(deltaFract));
            a0 = _mm_shuffle_epi8(a0, _mm_setr_epi8(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1));
            top0 = _mm_loadu_si128((__m128i *)(ref_samp_main + deltaInt));
            top0 = _mm_shuffle_epi8(top0, _mm_setr_epi8(0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8));
            sum0 = _mm_add_epi16(_mm_set1_epi16(16), _mm_maddubs_epi16(top0, a0));
            sum0 = _mm_srai_epi16(sum0, 5);
            sum0 = _mm_packus_epi16(sum0, sum0);
            _mm_storel_epi64((__m128i *)prediction_ptr, sum0);
            prediction_ptr += prediction_buffer_stride;
            deltaSum += intra_pred_angle;
        }
    }
    else if (size == 16) {
        for (row_index = 0; row_index < height; row_index++) {
            deltaInt = deltaSum >> 5;
            deltaFract = deltaSum & 31;
            a0 = _mm_unpacklo_epi8(_mm_cvtsi32_si128(32 - deltaFract), _mm_cvtsi32_si128(deltaFract));
            a0 = _mm_shuffle_epi8(a0, _mm_setr_epi8(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1));
            top0 = _mm_loadu_si128((__m128i *)(ref_samp_main + deltaInt));
            top1 = _mm_loadu_si128((__m128i *)(ref_samp_main + deltaInt + 1));
            top2 = _mm_unpacklo_epi8(top0, top1);
            top0 = _mm_unpackhi_epi8(top0, top1);
            sum0 = _mm_add_epi16(_mm_set1_epi16(16), _mm_maddubs_epi16(top2, a0));
            sum1 = _mm_add_epi16(_mm_set1_epi16(16), _mm_maddubs_epi16(top0, a0));
            sum0 = _mm_srai_epi16(sum0, 5);
            sum1 = _mm_srai_epi16(sum1, 5);
            sum0 = _mm_packus_epi16(sum0, sum1);
            _mm_storeu_si128((__m128i *)prediction_ptr, sum0);
            prediction_ptr += prediction_buffer_stride;
            deltaSum += intra_pred_angle;
        }
    }
    else { // size == 32
        for (row_index = 0; row_index < height; row_index++) {
            deltaInt = deltaSum >> 5;
            deltaFract = deltaSum & 31;
            a0 = _mm_unpacklo_epi8(_mm_cvtsi32_si128(32 - deltaFract), _mm_cvtsi32_si128(deltaFract));
            a0 = _mm_shuffle_epi8(a0, _mm_setr_epi8(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1));
            top0 = _mm_loadu_si128((__m128i *)(ref_samp_main + deltaInt));
            top1 = _mm_loadu_si128((__m128i *)(ref_samp_main + deltaInt + 1));
            top2 = _mm_unpacklo_epi8(top0, top1);
            top0 = _mm_unpackhi_epi8(top0, top1);
            sum0 = _mm_add_epi16(_mm_set1_epi16(16), _mm_maddubs_epi16(top2, a0));
            sum1 = _mm_add_epi16(_mm_set1_epi16(16), _mm_maddubs_epi16(top0, a0));
            sum0 = _mm_srai_epi16(sum0, 5);
            sum1 = _mm_srai_epi16(sum1, 5);
            sum0 = _mm_packus_epi16(sum0, sum1);
            _mm_storeu_si128((__m128i *)prediction_ptr, sum0);
            top0 = _mm_loadu_si128((__m128i *)(ref_samp_main + deltaInt + 16));
            top1 = _mm_loadu_si128((__m128i *)(ref_samp_main + deltaInt + 17));
            top2 = _mm_unpacklo_epi8(top0, top1);
            top0 = _mm_unpackhi_epi8(top0, top1);
            sum0 = _mm_add_epi16(_mm_set1_epi16(16), _mm_maddubs_epi16(top2, a0));
            sum1 = _mm_add_epi16(_mm_set1_epi16(16), _mm_maddubs_epi16(top0, a0));
            sum0 = _mm_srai_epi16(sum0, 5);
            sum1 = _mm_srai_epi16(sum1, 5);
            sum0 = _mm_packus_epi16(sum0, sum1);
            _mm_storeu_si128((__m128i *)(prediction_ptr + 16), sum0);
            prediction_ptr += prediction_buffer_stride;
            deltaSum += intra_pred_angle;
        }
    }
}

extern void intra_mode_angular_horizontal_kernel_ssse3_intrin(
    uint32_t         size,                       //input parameter, denotes the size of the current PU
    uint8_t         *ref_samp_main,                //input parameter, pointer to the reference samples
    uint8_t         *prediction_ptr,              //output parameter, pointer to the prediction
    uint32_t         prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
    const EbBool  skip,
    int32_t         intra_pred_angle)
{
    uint32_t row_index, colIndex;
    int32_t deltaSum = 0;
    int32_t deltaInt;
    uint32_t deltaFract;
    __m128i top0, top1, top2, sum0, sum1, a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11;
    uint8_t temp_buf[32 * 32];
    uint8_t *p = temp_buf;

    // --------- Reference Samples Structure ---------
    // ref_samp_main[-size+1] to ref_samp_main[-1] must be prepared (from right to left) for mode 11 to 17 (not required for mode 3 to 9)
    // ref_samp_main[0]      = TopLeft[0]
    // ref_samp_main[1]      = Left[0]
    // ref_samp_main[2]      = Left[1]
    // ...
    // ref_samp_main[size]   = Left[size-1]
    // ref_samp_main[size+1] = Left[size]     for mode 3 to 9 (not required for mode 11 to 17)
    // ...
    // ref_samp_main[2*size] = Left[2*size-1] for mode 3 to 9 (not required for mode 11 to 17)
    // -----------------------------------------------

    // Compute the prediction
    ref_samp_main += 1; // left sample

    if (skip) {
        prediction_buffer_stride <<= 1;
        if (size == 4) {
            deltaSum += intra_pred_angle;
            deltaInt = deltaSum >> 5;
            deltaFract = deltaSum & 31;
            a0 = _mm_unpacklo_epi8(_mm_cvtsi32_si128(32 - deltaFract), _mm_cvtsi32_si128(deltaFract));
            top0 = _mm_cvtsi32_si128(*(uint32_t *)(ref_samp_main + deltaInt));

            deltaSum += intra_pred_angle;
            deltaInt = deltaSum >> 5;
            deltaFract = deltaSum & 31;
            a1 = _mm_unpacklo_epi8(_mm_cvtsi32_si128(32 - deltaFract), _mm_cvtsi32_si128(deltaFract));
            a0 = _mm_unpacklo_epi16(a0, a1);
            top0 = _mm_unpacklo_epi32(top0, _mm_cvtsi32_si128(*(uint32_t *)(ref_samp_main + deltaInt)));

            deltaSum += intra_pred_angle;
            deltaInt = deltaSum >> 5;
            deltaFract = deltaSum & 31;
            a2 = _mm_unpacklo_epi8(_mm_cvtsi32_si128(32 - deltaFract), _mm_cvtsi32_si128(deltaFract));
            a4 = _mm_cvtsi32_si128(*(uint32_t *)(ref_samp_main + deltaInt));

            deltaSum += intra_pred_angle;
            deltaInt = deltaSum >> 5;
            deltaFract = deltaSum & 31;
            a3 = _mm_unpacklo_epi8(_mm_cvtsi32_si128(32 - deltaFract), _mm_cvtsi32_si128(deltaFract));
            a2 = _mm_unpacklo_epi16(a2, a3);
            a0 = _mm_unpacklo_epi32(a0, a2);
            a0 = _mm_unpacklo_epi16(a0, a0);
            a4 = _mm_unpacklo_epi32(a4, _mm_cvtsi32_si128(*(uint32_t *)(ref_samp_main + deltaInt)));

            top0 = _mm_unpacklo_epi64(top0, a4);
            sum0 = _mm_add_epi16(_mm_set1_epi16(16), _mm_maddubs_epi16(top0, a0));
            sum0 = _mm_srai_epi16(sum0, 5);
            sum0 = _mm_packus_epi16(sum0, sum0);
            sum0 = _mm_shuffle_epi8(sum0, _mm_setr_epi8(0, 2, 4, 6, 1, 3, 5, 7, 2, 6, 10, 14, 3, 7, 11, 15));
            *(uint32_t *)prediction_ptr = _mm_cvtsi128_si32(sum0); sum0 = _mm_srli_si128(sum0, 4); prediction_ptr += prediction_buffer_stride;
            *(uint32_t *)prediction_ptr = _mm_cvtsi128_si32(sum0);
        }
        else if (size == 8) {
            for (colIndex = 0; colIndex < size; colIndex++) {
                deltaSum += intra_pred_angle;
                deltaInt = deltaSum >> 5;
                deltaFract = deltaSum & 31;
                a0 = _mm_unpacklo_epi8(_mm_cvtsi32_si128(32 - deltaFract), _mm_cvtsi32_si128(deltaFract));
                a0 = _mm_shuffle_epi8(a0, _mm_setr_epi8(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1));
                top0 = _mm_loadl_epi64((__m128i *)(ref_samp_main + deltaInt));
                sum0 = _mm_add_epi16(_mm_set1_epi16(16), _mm_maddubs_epi16(top0, a0));
                sum0 = _mm_srai_epi16(sum0, 5);
                sum0 = _mm_packus_epi16(sum0, sum0);
                *(uint32_t *)p = _mm_cvtsi128_si32(sum0);
                p += 4;
            }
            a0 = _mm_loadu_si128((__m128i *)temp_buf);
            a0 = _mm_shuffle_epi8(a0, _mm_setr_epi8(0, 4, 8, 12, 1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15));
            a1 = _mm_loadu_si128((__m128i *)(temp_buf + 16));
            a1 = _mm_shuffle_epi8(a1, _mm_setr_epi8(0, 4, 8, 12, 1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15));
            a2 = _mm_unpackhi_epi32(a0, a1);
            a0 = _mm_unpacklo_epi32(a0, a1);
            a1 = _mm_unpackhi_epi64(a0, a2);
            a0 = _mm_unpacklo_epi64(a0, a2);
            _mm_storel_epi64((__m128i *)prediction_ptr, a0); prediction_ptr += prediction_buffer_stride;
            _mm_storel_epi64((__m128i *)prediction_ptr, a1); prediction_ptr += prediction_buffer_stride;
            _mm_storeh_epi64((__m128i *)prediction_ptr, a0); prediction_ptr += prediction_buffer_stride;
            _mm_storeh_epi64((__m128i *)prediction_ptr, a1);
        }
        else if (size == 16) {
            for (colIndex = 0; colIndex < size; colIndex++) {
                deltaSum += intra_pred_angle;
                deltaInt = deltaSum >> 5;
                deltaFract = deltaSum & 31;
                a0 = _mm_unpacklo_epi8(_mm_cvtsi32_si128(32 - deltaFract), _mm_cvtsi32_si128(deltaFract));
                a0 = _mm_shuffle_epi8(a0, _mm_setr_epi8(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1));
                top0 = _mm_loadu_si128((__m128i *)(ref_samp_main + deltaInt));
                sum0 = _mm_add_epi16(_mm_set1_epi16(16), _mm_maddubs_epi16(top0, a0));
                sum0 = _mm_srai_epi16(sum0, 5);
                sum0 = _mm_packus_epi16(sum0, sum0);
                _mm_storel_epi64((__m128i *)p, sum0);
                p += 8;
            }
            p = temp_buf;
            a0 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x00)), _mm_loadl_epi64((__m128i *)(p + 0x08)));
            a1 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x10)), _mm_loadl_epi64((__m128i *)(p + 0x18)));
            a2 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x20)), _mm_loadl_epi64((__m128i *)(p + 0x28)));
            a3 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x30)), _mm_loadl_epi64((__m128i *)(p + 0x38)));
            a4 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x40)), _mm_loadl_epi64((__m128i *)(p + 0x48)));
            a5 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x50)), _mm_loadl_epi64((__m128i *)(p + 0x58)));
            a6 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x60)), _mm_loadl_epi64((__m128i *)(p + 0x68)));
            a7 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x70)), _mm_loadl_epi64((__m128i *)(p + 0x78)));

            a8 = _mm_unpackhi_epi16(a0, a1);
            a0 = _mm_unpacklo_epi16(a0, a1);
            a9 = _mm_unpackhi_epi16(a2, a3);
            a2 = _mm_unpacklo_epi16(a2, a3);
            a10 = _mm_unpackhi_epi16(a4, a5);
            a4 = _mm_unpacklo_epi16(a4, a5);
            a11 = _mm_unpackhi_epi16(a6, a7);
            a6 = _mm_unpacklo_epi16(a6, a7);

            a1 = _mm_unpackhi_epi32(a0, a2);
            a0 = _mm_unpacklo_epi32(a0, a2);
            a3 = _mm_unpackhi_epi32(a4, a6);
            a4 = _mm_unpacklo_epi32(a4, a6);
            a5 = _mm_unpackhi_epi32(a8, a9);
            a8 = _mm_unpacklo_epi32(a8, a9);
            a7 = _mm_unpackhi_epi32(a10, a11);
            a10 = _mm_unpacklo_epi32(a10, a11);

            a2 = _mm_unpackhi_epi64(a0, a4);
            a0 = _mm_unpacklo_epi64(a0, a4);
            a6 = _mm_unpackhi_epi64(a8, a10);
            a8 = _mm_unpacklo_epi64(a8, a10);
            a9 = _mm_unpackhi_epi64(a1, a3);
            a1 = _mm_unpacklo_epi64(a1, a3);
            a11 = _mm_unpackhi_epi64(a5, a7);
            a5 = _mm_unpacklo_epi64(a5, a7);

            _mm_storeu_si128((__m128i *)prediction_ptr, a0);  prediction_ptr += prediction_buffer_stride;
            _mm_storeu_si128((__m128i *)prediction_ptr, a2);  prediction_ptr += prediction_buffer_stride;
            _mm_storeu_si128((__m128i *)prediction_ptr, a1);  prediction_ptr += prediction_buffer_stride;
            _mm_storeu_si128((__m128i *)prediction_ptr, a9);  prediction_ptr += prediction_buffer_stride;
            _mm_storeu_si128((__m128i *)prediction_ptr, a8);  prediction_ptr += prediction_buffer_stride;
            _mm_storeu_si128((__m128i *)prediction_ptr, a6);  prediction_ptr += prediction_buffer_stride;
            _mm_storeu_si128((__m128i *)prediction_ptr, a5);  prediction_ptr += prediction_buffer_stride;
            _mm_storeu_si128((__m128i *)prediction_ptr, a11); prediction_ptr += prediction_buffer_stride;
        }
        else { // size == 32
            for (colIndex = 0; colIndex < size; colIndex++) {
                deltaSum += intra_pred_angle;
                deltaInt = deltaSum >> 5;
                deltaFract = deltaSum & 31;
                a0 = _mm_unpacklo_epi8(_mm_cvtsi32_si128(32 - deltaFract), _mm_cvtsi32_si128(deltaFract));
                a0 = _mm_shuffle_epi8(a0, _mm_setr_epi8(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1));
                top0 = _mm_loadu_si128((__m128i *)(ref_samp_main + deltaInt));
                top1 = _mm_loadu_si128((__m128i *)(ref_samp_main + deltaInt + 16));
                sum0 = _mm_add_epi16(_mm_set1_epi16(16), _mm_maddubs_epi16(top0, a0));
                sum1 = _mm_add_epi16(_mm_set1_epi16(16), _mm_maddubs_epi16(top1, a0));
                sum0 = _mm_srai_epi16(sum0, 5);
                sum1 = _mm_srai_epi16(sum1, 5);
                sum0 = _mm_packus_epi16(sum0, sum1);
                _mm_storeu_si128((__m128i *)p, sum0);
                p += 16;
            }
            p = temp_buf;
            for (colIndex = 0; colIndex < 2; colIndex++) {
                for (row_index = 0; row_index < 2; row_index++) {
                    a0 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x00)), _mm_loadl_epi64((__m128i *)(p + 0x10)));
                    a1 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x20)), _mm_loadl_epi64((__m128i *)(p + 0x30)));
                    a2 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x40)), _mm_loadl_epi64((__m128i *)(p + 0x50)));
                    a3 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x60)), _mm_loadl_epi64((__m128i *)(p + 0x70)));
                    a4 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x80)), _mm_loadl_epi64((__m128i *)(p + 0x90)));
                    a5 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0xA0)), _mm_loadl_epi64((__m128i *)(p + 0xB0)));
                    a6 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0xC0)), _mm_loadl_epi64((__m128i *)(p + 0xD0)));
                    a7 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0xE0)), _mm_loadl_epi64((__m128i *)(p + 0xF0)));

                    a8 = _mm_unpackhi_epi16(a0, a1);
                    a0 = _mm_unpacklo_epi16(a0, a1);
                    a9 = _mm_unpackhi_epi16(a2, a3);
                    a2 = _mm_unpacklo_epi16(a2, a3);
                    a10 = _mm_unpackhi_epi16(a4, a5);
                    a4 = _mm_unpacklo_epi16(a4, a5);
                    a11 = _mm_unpackhi_epi16(a6, a7);
                    a6 = _mm_unpacklo_epi16(a6, a7);

                    a1 = _mm_unpackhi_epi32(a0, a2);
                    a0 = _mm_unpacklo_epi32(a0, a2);
                    a3 = _mm_unpackhi_epi32(a4, a6);
                    a4 = _mm_unpacklo_epi32(a4, a6);
                    a5 = _mm_unpackhi_epi32(a8, a9);
                    a8 = _mm_unpacklo_epi32(a8, a9);
                    a7 = _mm_unpackhi_epi32(a10, a11);
                    a10 = _mm_unpacklo_epi32(a10, a11);

                    a2 = _mm_unpackhi_epi64(a0, a4);
                    a0 = _mm_unpacklo_epi64(a0, a4);
                    a6 = _mm_unpackhi_epi64(a8, a10);
                    a8 = _mm_unpacklo_epi64(a8, a10);
                    a9 = _mm_unpackhi_epi64(a1, a3);
                    a1 = _mm_unpacklo_epi64(a1, a3);
                    a11 = _mm_unpackhi_epi64(a5, a7);
                    a5 = _mm_unpacklo_epi64(a5, a7);

                    _mm_storeu_si128((__m128i *)prediction_ptr, a0);  prediction_ptr += prediction_buffer_stride;
                    _mm_storeu_si128((__m128i *)prediction_ptr, a2);  prediction_ptr += prediction_buffer_stride;
                    _mm_storeu_si128((__m128i *)prediction_ptr, a1);  prediction_ptr += prediction_buffer_stride;
                    _mm_storeu_si128((__m128i *)prediction_ptr, a9);  prediction_ptr += prediction_buffer_stride;
                    _mm_storeu_si128((__m128i *)prediction_ptr, a8);  prediction_ptr += prediction_buffer_stride;
                    _mm_storeu_si128((__m128i *)prediction_ptr, a6);  prediction_ptr += prediction_buffer_stride;
                    _mm_storeu_si128((__m128i *)prediction_ptr, a5);  prediction_ptr += prediction_buffer_stride;
                    _mm_storeu_si128((__m128i *)prediction_ptr, a11); prediction_ptr += prediction_buffer_stride;
                    p += 8;
                }
                p = temp_buf + 0x100;
                prediction_ptr -= 16 * prediction_buffer_stride - 16;
            }
        }
    }
    else {
        if (size == 4) {
            for (colIndex = 0; colIndex < size; colIndex += 2) {
                deltaSum += intra_pred_angle;
                deltaInt = deltaSum >> 5;
                deltaFract = deltaSum & 31;
                a0 = _mm_unpacklo_epi8(_mm_cvtsi32_si128(32 - deltaFract), _mm_cvtsi32_si128(deltaFract));
                top0 = _mm_loadl_epi64((__m128i *)(ref_samp_main + deltaInt));
                deltaSum += intra_pred_angle;
                deltaInt = deltaSum >> 5;
                deltaFract = deltaSum & 31;
                a1 = _mm_unpacklo_epi8(_mm_cvtsi32_si128(32 - deltaFract), _mm_cvtsi32_si128(deltaFract));
                a0 = _mm_unpacklo_epi16(a0, a1);
                a0 = _mm_shuffle_epi8(a0, _mm_setr_epi8(0, 1, 0, 1, 0, 1, 0, 1, 2, 3, 2, 3, 2, 3, 2, 3));
                top0 = _mm_castps_si128(_mm_loadh_pi(_mm_castsi128_ps(top0), (__m64 *)(ref_samp_main + deltaInt)));
                top0 = _mm_shuffle_epi8(top0, _mm_setr_epi8(0, 1, 1, 2, 2, 3, 3, 4, 8, 9, 9, 10, 10, 11, 11, 12));
                sum0 = _mm_add_epi16(_mm_set1_epi16(16), _mm_maddubs_epi16(top0, a0));
                sum0 = _mm_srai_epi16(sum0, 5);
                sum0 = _mm_packus_epi16(sum0, sum0);
                *(uint32_t *)p = _mm_cvtsi128_si32(sum0);
                *(uint32_t *)(p + 4) = _mm_cvtsi128_si32(_mm_srli_si128(sum0, 4));
                p += 8;
            }
            a0 = _mm_loadu_si128((__m128i *)temp_buf);
            a0 = _mm_shuffle_epi8(a0, _mm_setr_epi8(0, 4, 8, 12, 1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15));
            *(uint32_t *)prediction_ptr = _mm_cvtsi128_si32(a0); a0 = _mm_srli_si128(a0, 4); prediction_ptr += prediction_buffer_stride;
            *(uint32_t *)prediction_ptr = _mm_cvtsi128_si32(a0); a0 = _mm_srli_si128(a0, 4); prediction_ptr += prediction_buffer_stride;
            *(uint32_t *)prediction_ptr = _mm_cvtsi128_si32(a0); a0 = _mm_srli_si128(a0, 4); prediction_ptr += prediction_buffer_stride;
            *(uint32_t *)prediction_ptr = _mm_cvtsi128_si32(a0);
        }
        else if (size == 8) {
            for (colIndex = 0; colIndex < size; colIndex++) {
                deltaSum += intra_pred_angle;
                deltaInt = deltaSum >> 5;
                deltaFract = deltaSum & 31;
                a0 = _mm_unpacklo_epi8(_mm_cvtsi32_si128(32 - deltaFract), _mm_cvtsi32_si128(deltaFract));
                a0 = _mm_shuffle_epi8(a0, _mm_setr_epi8(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1));
                top0 = _mm_loadu_si128((__m128i *)(ref_samp_main + deltaInt));
                top0 = _mm_shuffle_epi8(top0, _mm_setr_epi8(0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8));
                sum0 = _mm_add_epi16(_mm_set1_epi16(16), _mm_maddubs_epi16(top0, a0));
                sum0 = _mm_srai_epi16(sum0, 5);
                sum0 = _mm_packus_epi16(sum0, sum0);
                _mm_storel_epi64((__m128i *)p, sum0);
                p += 8;
            }
            a0 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(temp_buf + 0x00)), _mm_loadl_epi64((__m128i *)(temp_buf + 0x08)));
            a1 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(temp_buf + 0x10)), _mm_loadl_epi64((__m128i *)(temp_buf + 0x18)));
            a2 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(temp_buf + 0x20)), _mm_loadl_epi64((__m128i *)(temp_buf + 0x28)));
            a3 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(temp_buf + 0x30)), _mm_loadl_epi64((__m128i *)(temp_buf + 0x38)));
            a4 = _mm_unpackhi_epi16(a0, a1);
            a0 = _mm_unpacklo_epi16(a0, a1);
            a5 = _mm_unpackhi_epi16(a2, a3);
            a2 = _mm_unpacklo_epi16(a2, a3);
            a1 = _mm_unpackhi_epi32(a0, a2);
            a0 = _mm_unpacklo_epi32(a0, a2);
            a3 = _mm_unpackhi_epi32(a4, a5);
            a4 = _mm_unpacklo_epi32(a4, a5);
            _mm_storel_epi64((__m128i *)prediction_ptr, a0); prediction_ptr += prediction_buffer_stride;
            _mm_storeh_epi64((__m128i *)prediction_ptr, a0); prediction_ptr += prediction_buffer_stride;
            _mm_storel_epi64((__m128i *)prediction_ptr, a1); prediction_ptr += prediction_buffer_stride;
            _mm_storeh_epi64((__m128i *)prediction_ptr, a1); prediction_ptr += prediction_buffer_stride;
            _mm_storel_epi64((__m128i *)prediction_ptr, a4); prediction_ptr += prediction_buffer_stride;
            _mm_storeh_epi64((__m128i *)prediction_ptr, a4); prediction_ptr += prediction_buffer_stride;
            _mm_storel_epi64((__m128i *)prediction_ptr, a3); prediction_ptr += prediction_buffer_stride;
            _mm_storeh_epi64((__m128i *)prediction_ptr, a3);
        }
        else if (size == 16) {
            for (colIndex = 0; colIndex < size; colIndex++) {
                deltaSum += intra_pred_angle;
                deltaInt = deltaSum >> 5;
                deltaFract = deltaSum & 31;
                a0 = _mm_unpacklo_epi8(_mm_cvtsi32_si128(32 - deltaFract), _mm_cvtsi32_si128(deltaFract));
                a0 = _mm_shuffle_epi8(a0, _mm_setr_epi8(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1));
                top0 = _mm_loadu_si128((__m128i *)(ref_samp_main + deltaInt));
                top1 = _mm_loadu_si128((__m128i *)(ref_samp_main + deltaInt + 1));
                top2 = _mm_unpacklo_epi8(top0, top1);
                top0 = _mm_unpackhi_epi8(top0, top1);
                sum0 = _mm_add_epi16(_mm_set1_epi16(16), _mm_maddubs_epi16(top2, a0));
                sum1 = _mm_add_epi16(_mm_set1_epi16(16), _mm_maddubs_epi16(top0, a0));
                sum0 = _mm_srai_epi16(sum0, 5);
                sum1 = _mm_srai_epi16(sum1, 5);
                sum0 = _mm_packus_epi16(sum0, sum1);
                _mm_storeu_si128((__m128i *)p, sum0);
                p += 16;
            }
            p = temp_buf;
            for (colIndex = 0; colIndex < 2; colIndex++) {
                a0 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x00)), _mm_loadl_epi64((__m128i *)(p + 0x10)));
                a1 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x20)), _mm_loadl_epi64((__m128i *)(p + 0x30)));
                a2 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x40)), _mm_loadl_epi64((__m128i *)(p + 0x50)));
                a3 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x60)), _mm_loadl_epi64((__m128i *)(p + 0x70)));
                a4 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x80)), _mm_loadl_epi64((__m128i *)(p + 0x90)));
                a5 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0xA0)), _mm_loadl_epi64((__m128i *)(p + 0xB0)));
                a6 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0xC0)), _mm_loadl_epi64((__m128i *)(p + 0xD0)));
                a7 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0xE0)), _mm_loadl_epi64((__m128i *)(p + 0xF0)));

                a8 = _mm_unpackhi_epi16(a0, a1);
                a0 = _mm_unpacklo_epi16(a0, a1);
                a9 = _mm_unpackhi_epi16(a2, a3);
                a2 = _mm_unpacklo_epi16(a2, a3);
                a10 = _mm_unpackhi_epi16(a4, a5);
                a4 = _mm_unpacklo_epi16(a4, a5);
                a11 = _mm_unpackhi_epi16(a6, a7);
                a6 = _mm_unpacklo_epi16(a6, a7);

                a1 = _mm_unpackhi_epi32(a0, a2);
                a0 = _mm_unpacklo_epi32(a0, a2);
                a3 = _mm_unpackhi_epi32(a4, a6);
                a4 = _mm_unpacklo_epi32(a4, a6);
                a5 = _mm_unpackhi_epi32(a8, a9);
                a8 = _mm_unpacklo_epi32(a8, a9);
                a7 = _mm_unpackhi_epi32(a10, a11);
                a10 = _mm_unpacklo_epi32(a10, a11);

                a2 = _mm_unpackhi_epi64(a0, a4);
                a0 = _mm_unpacklo_epi64(a0, a4);
                a6 = _mm_unpackhi_epi64(a8, a10);
                a8 = _mm_unpacklo_epi64(a8, a10);
                a9 = _mm_unpackhi_epi64(a1, a3);
                a1 = _mm_unpacklo_epi64(a1, a3);
                a11 = _mm_unpackhi_epi64(a5, a7);
                a5 = _mm_unpacklo_epi64(a5, a7);

                _mm_storeu_si128((__m128i *)prediction_ptr, a0);  prediction_ptr += prediction_buffer_stride;
                _mm_storeu_si128((__m128i *)prediction_ptr, a2);  prediction_ptr += prediction_buffer_stride;
                _mm_storeu_si128((__m128i *)prediction_ptr, a1);  prediction_ptr += prediction_buffer_stride;
                _mm_storeu_si128((__m128i *)prediction_ptr, a9);  prediction_ptr += prediction_buffer_stride;
                _mm_storeu_si128((__m128i *)prediction_ptr, a8);  prediction_ptr += prediction_buffer_stride;
                _mm_storeu_si128((__m128i *)prediction_ptr, a6);  prediction_ptr += prediction_buffer_stride;
                _mm_storeu_si128((__m128i *)prediction_ptr, a5);  prediction_ptr += prediction_buffer_stride;
                _mm_storeu_si128((__m128i *)prediction_ptr, a11); prediction_ptr += prediction_buffer_stride;
                p += 8;
            }
        }
        else { // size == 32
            for (colIndex = 0; colIndex < size; colIndex++) {
                deltaSum += intra_pred_angle;
                deltaInt = deltaSum >> 5;
                deltaFract = deltaSum & 31;
                a0 = _mm_unpacklo_epi8(_mm_cvtsi32_si128(32 - deltaFract), _mm_cvtsi32_si128(deltaFract));
                a0 = _mm_shuffle_epi8(a0, _mm_setr_epi8(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1));
                top0 = _mm_loadu_si128((__m128i *)(ref_samp_main + deltaInt));
                top1 = _mm_loadu_si128((__m128i *)(ref_samp_main + deltaInt + 1));
                top2 = _mm_unpacklo_epi8(top0, top1);
                top0 = _mm_unpackhi_epi8(top0, top1);
                sum0 = _mm_add_epi16(_mm_set1_epi16(16), _mm_maddubs_epi16(top2, a0));
                sum1 = _mm_add_epi16(_mm_set1_epi16(16), _mm_maddubs_epi16(top0, a0));
                sum0 = _mm_srai_epi16(sum0, 5);
                sum1 = _mm_srai_epi16(sum1, 5);
                sum0 = _mm_packus_epi16(sum0, sum1);
                _mm_storeu_si128((__m128i *)p, sum0);
                top0 = _mm_loadu_si128((__m128i *)(ref_samp_main + deltaInt + 16));
                top1 = _mm_loadu_si128((__m128i *)(ref_samp_main + deltaInt + 17));
                top2 = _mm_unpacklo_epi8(top0, top1);
                top0 = _mm_unpackhi_epi8(top0, top1);
                sum0 = _mm_add_epi16(_mm_set1_epi16(16), _mm_maddubs_epi16(top2, a0));
                sum1 = _mm_add_epi16(_mm_set1_epi16(16), _mm_maddubs_epi16(top0, a0));
                sum0 = _mm_srai_epi16(sum0, 5);
                sum1 = _mm_srai_epi16(sum1, 5);
                sum0 = _mm_packus_epi16(sum0, sum1);
                _mm_storeu_si128((__m128i *)(p + 16), sum0);
                p += 32;
            }
            p = temp_buf;
            for (colIndex = 0; colIndex < 2; colIndex++) {
                for (row_index = 0; row_index < 4; row_index++) {
                    a0 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x000)), _mm_loadl_epi64((__m128i *)(p + 0x020)));
                    a1 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x040)), _mm_loadl_epi64((__m128i *)(p + 0x060)));
                    a2 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x080)), _mm_loadl_epi64((__m128i *)(p + 0x0A0)));
                    a3 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x0C0)), _mm_loadl_epi64((__m128i *)(p + 0x0E0)));
                    a4 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x100)), _mm_loadl_epi64((__m128i *)(p + 0x120)));
                    a5 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x140)), _mm_loadl_epi64((__m128i *)(p + 0x160)));
                    a6 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x180)), _mm_loadl_epi64((__m128i *)(p + 0x1A0)));
                    a7 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x1C0)), _mm_loadl_epi64((__m128i *)(p + 0x1E0)));

                    a8 = _mm_unpackhi_epi16(a0, a1);
                    a0 = _mm_unpacklo_epi16(a0, a1);
                    a9 = _mm_unpackhi_epi16(a2, a3);
                    a2 = _mm_unpacklo_epi16(a2, a3);
                    a10 = _mm_unpackhi_epi16(a4, a5);
                    a4 = _mm_unpacklo_epi16(a4, a5);
                    a11 = _mm_unpackhi_epi16(a6, a7);
                    a6 = _mm_unpacklo_epi16(a6, a7);

                    a1 = _mm_unpackhi_epi32(a0, a2);
                    a0 = _mm_unpacklo_epi32(a0, a2);
                    a3 = _mm_unpackhi_epi32(a4, a6);
                    a4 = _mm_unpacklo_epi32(a4, a6);
                    a5 = _mm_unpackhi_epi32(a8, a9);
                    a8 = _mm_unpacklo_epi32(a8, a9);
                    a7 = _mm_unpackhi_epi32(a10, a11);
                    a10 = _mm_unpacklo_epi32(a10, a11);

                    a2 = _mm_unpackhi_epi64(a0, a4);
                    a0 = _mm_unpacklo_epi64(a0, a4);
                    a6 = _mm_unpackhi_epi64(a8, a10);
                    a8 = _mm_unpacklo_epi64(a8, a10);
                    a9 = _mm_unpackhi_epi64(a1, a3);
                    a1 = _mm_unpacklo_epi64(a1, a3);
                    a11 = _mm_unpackhi_epi64(a5, a7);
                    a5 = _mm_unpacklo_epi64(a5, a7);

                    _mm_storeu_si128((__m128i *)prediction_ptr, a0);  prediction_ptr += prediction_buffer_stride;
                    _mm_storeu_si128((__m128i *)prediction_ptr, a2);  prediction_ptr += prediction_buffer_stride;
                    _mm_storeu_si128((__m128i *)prediction_ptr, a1);  prediction_ptr += prediction_buffer_stride;
                    _mm_storeu_si128((__m128i *)prediction_ptr, a9);  prediction_ptr += prediction_buffer_stride;
                    _mm_storeu_si128((__m128i *)prediction_ptr, a8);  prediction_ptr += prediction_buffer_stride;
                    _mm_storeu_si128((__m128i *)prediction_ptr, a6);  prediction_ptr += prediction_buffer_stride;
                    _mm_storeu_si128((__m128i *)prediction_ptr, a5);  prediction_ptr += prediction_buffer_stride;
                    _mm_storeu_si128((__m128i *)prediction_ptr, a11); prediction_ptr += prediction_buffer_stride;
                    p += 8;
                }
                p = temp_buf + 0x200;
                prediction_ptr -= 32 * prediction_buffer_stride - 16;
            }
        }
    }
}

