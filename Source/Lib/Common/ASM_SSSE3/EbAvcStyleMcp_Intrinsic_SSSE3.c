/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbAvcStyleMcp_SSSE3.h"

#include "EbMcp_SSE2.h"
#include "EbDefinitions.h"
#include "EbAvcStyleMcp_SSE2.h"

#include "emmintrin.h"
#include "tmmintrin.h"

EB_EXTERN EB_ALIGN(16) const int8_t AvcStyleLumaIFCoeff8_SSSE3[] = {
    -1, 25, -1, 25, -1, 25, -1, 25, -1, 25, -1, 25, -1, 25, -1, 25,
     9, -1,  9, -1,  9, -1,  9, -1,  9, -1,  9, -1,  9, -1,  9, -1,
    -2, 18, -2, 18, -2, 18, -2, 18, -2, 18, -2, 18, -2, 18, -2, 18,
    18, -2, 18, -2, 18, -2, 18, -2, 18, -2, 18, -2, 18, -2, 18, -2,
    -1,  9, -1,  9, -1,  9, -1,  9, -1,  9, -1,  9, -1,  9, -1,  9,
    25, -1, 25, -1, 25, -1, 25, -1, 25, -1, 25, -1, 25, -1, 25, -1
};

void avc_style_luma_interpolation_filter_pose_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    EbByte               dst,
    uint32_t                dst_stride,
    uint32_t                pu_width,
    uint32_t                pu_height,
    EbByte               temp_buf,
    EbBool               skip,
    uint32_t                frac_pos)
{
    uint32_t tempBufSize = pu_width * pu_height;
    (void)frac_pos;
    avc_style_luma_interpolation_filter_horizontal_ssse3_intrin(ref_pic, src_stride, temp_buf, pu_width, pu_width, pu_height, 0, skip, 2);
    avc_style_luma_interpolation_filter_vertical_ssse3_intrin(ref_pic, src_stride, temp_buf + tempBufSize, pu_width, pu_width, pu_height, 0, skip, 2);
    picture_average_kernel_sse2_intrin(temp_buf, pu_width, temp_buf + tempBufSize, pu_width, dst, dst_stride, pu_width, pu_height);
}

void avc_style_luma_interpolation_filter_posf_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    EbByte               dst,
    uint32_t                dst_stride,
    uint32_t                pu_width,
    uint32_t                pu_height,
    EbByte               temp_buf,
    EbBool               skip,
    uint32_t                frac_pos)
{
    uint32_t tempBufSize = pu_width * pu_height;
    (void)frac_pos;
    avc_style_luma_interpolation_filter_horizontal_ssse3_intrin(ref_pic - src_stride, src_stride, temp_buf + tempBufSize, pu_width, pu_width, skip ? (2 * pu_height + 3) : (pu_height + 3), 0, EB_FALSE, 2);
    avc_style_luma_interpolation_filter_vertical_ssse3_intrin(temp_buf + tempBufSize + pu_width, pu_width, temp_buf, pu_width, pu_width, pu_height, 0, skip, 2);
    picture_average_kernel_sse2_intrin(temp_buf + tempBufSize + pu_width, skip ? 2 * pu_width : pu_width, temp_buf, pu_width, dst, dst_stride, pu_width, pu_height);
}

void avc_style_luma_interpolation_filter_posg_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    EbByte               dst,
    uint32_t                dst_stride,
    uint32_t                pu_width,
    uint32_t                pu_height,
    EbByte               temp_buf,
    EbBool               skip,
    uint32_t                frac_pos)
{
    uint32_t tempBufSize = pu_width * pu_height;
    (void)frac_pos;
    avc_style_luma_interpolation_filter_horizontal_ssse3_intrin(ref_pic, src_stride, temp_buf, pu_width, pu_width, pu_height, 0, skip, 2);
    avc_style_luma_interpolation_filter_vertical_ssse3_intrin(ref_pic + 1, src_stride, temp_buf + tempBufSize, pu_width, pu_width, pu_height, 0, skip, 2);
    picture_average_kernel_sse2_intrin(temp_buf, pu_width, temp_buf + tempBufSize, pu_width, dst, dst_stride, pu_width, pu_height);
}

void avc_style_luma_interpolation_filter_posi_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    EbByte               dst,
    uint32_t                dst_stride,
    uint32_t                pu_width,
    uint32_t                pu_height,
    EbByte               temp_buf,
    EbBool               skip,
    uint32_t                frac_pos)
{
    uint32_t tempBufSize = pu_width * pu_height;
    (void)frac_pos;
    avc_style_luma_interpolation_filter_vertical_ssse3_intrin(ref_pic, src_stride, temp_buf, pu_width, pu_width, pu_height, 0, skip, 2);
    avc_style_luma_interpolation_filter_posj_ssse3(ref_pic, src_stride, temp_buf + tempBufSize, pu_width, pu_width, pu_height, temp_buf + 2 * tempBufSize, skip, 2);
    picture_average_kernel_sse2_intrin(temp_buf, pu_width, temp_buf + tempBufSize, pu_width, dst, dst_stride, pu_width, pu_height);
}

void avc_style_luma_interpolation_filter_posj_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    EbByte               dst,
    uint32_t                dst_stride,
    uint32_t                pu_width,
    uint32_t                pu_height,
    EbByte               temp_buf,
    EbBool               skip,
    uint32_t                frac_pos)
{
    (void)frac_pos;
    if (skip)
        avc_style_luma_interpolation_filter_horizontal_ssse3_intrin(ref_pic - src_stride, src_stride, temp_buf, pu_width, pu_width, (pu_height + 3), 0, EB_FALSE, 2);
    else
        avc_style_luma_interpolation_filter_horizontal_ssse3_intrin(ref_pic - src_stride, src_stride, temp_buf, pu_width, pu_width, skip ? (2 * pu_height + 3) : (pu_height + 3), 0, EB_FALSE, 2);

    avc_style_luma_interpolation_filter_vertical_ssse3_intrin(temp_buf + pu_width, pu_width, dst, dst_stride, pu_width, pu_height, 0, skip, 2);
}

void avc_style_luma_interpolation_filter_posk_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    EbByte               dst,
    uint32_t                dst_stride,
    uint32_t                pu_width,
    uint32_t                pu_height,
    EbByte               temp_buf,
    EbBool               skip,
    uint32_t                frac_pos)
{
    uint32_t tempBufSize = pu_width * pu_height;
    (void)frac_pos;
    avc_style_luma_interpolation_filter_vertical_ssse3_intrin(ref_pic + 1, src_stride, temp_buf, pu_width, pu_width, pu_height, 0, skip, 2);
    avc_style_luma_interpolation_filter_posj_ssse3(ref_pic, src_stride, temp_buf + tempBufSize, pu_width, pu_width, pu_height, temp_buf + 2 * tempBufSize, skip, 2);
    picture_average_kernel_sse2_intrin(temp_buf, pu_width, temp_buf + tempBufSize, pu_width, dst, dst_stride, pu_width, pu_height);
}

void avc_style_luma_interpolation_filter_posp_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    EbByte               dst,
    uint32_t                dst_stride,
    uint32_t                pu_width,
    uint32_t                pu_height,
    EbByte               temp_buf,
    EbBool               skip,
    uint32_t                frac_pos)
{
    uint32_t tempBufSize = pu_width * pu_height;
    (void)frac_pos;
    avc_style_luma_interpolation_filter_vertical_ssse3_intrin(ref_pic, src_stride, temp_buf, pu_width, pu_width, pu_height, 0, skip, 2);
    avc_style_luma_interpolation_filter_horizontal_ssse3_intrin(ref_pic + src_stride, src_stride, temp_buf + tempBufSize, pu_width, pu_width, pu_height, 0, skip, 2);
    picture_average_kernel_sse2_intrin(temp_buf, pu_width, temp_buf + tempBufSize, pu_width, dst, dst_stride, pu_width, pu_height);
}

void avc_style_luma_interpolation_filter_posq_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    EbByte               dst,
    uint32_t                dst_stride,
    uint32_t                pu_width,
    uint32_t                pu_height,
    EbByte               temp_buf,
    EbBool               skip,
    uint32_t                frac_pos)
{
    uint32_t tempBufSize = pu_width * pu_height;
    (void)frac_pos;
    avc_style_luma_interpolation_filter_horizontal_ssse3_intrin(ref_pic - src_stride, src_stride, temp_buf + tempBufSize, pu_width, pu_width, skip ? (2 * pu_height + 3) : (pu_height + 3), 0, EB_FALSE, 2);
    avc_style_luma_interpolation_filter_vertical_ssse3_intrin(temp_buf + tempBufSize + pu_width, pu_width, temp_buf, pu_width, pu_width, pu_height, 0, skip, 2);
    picture_average_kernel_sse2_intrin(temp_buf + tempBufSize + 2 * pu_width, skip ? 2 * pu_width : pu_width, temp_buf, pu_width, dst, dst_stride, pu_width, pu_height);
}

void avc_style_luma_interpolation_filter_posr_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    EbByte               dst,
    uint32_t                dst_stride,
    uint32_t                pu_width,
    uint32_t                pu_height,
    EbByte               temp_buf,
    EbBool               skip,
    uint32_t                frac_pos)
{
    uint32_t tempBufSize = pu_width * pu_height;
    (void)frac_pos;
    avc_style_luma_interpolation_filter_vertical_ssse3_intrin(ref_pic + 1, src_stride, temp_buf, pu_width, pu_width, pu_height, 0, skip, 2);
    avc_style_luma_interpolation_filter_horizontal_ssse3_intrin(ref_pic + src_stride, src_stride, temp_buf + tempBufSize, pu_width, pu_width, pu_height, 0, skip, 2);
    picture_average_kernel_sse2_intrin(temp_buf, pu_width, temp_buf + tempBufSize, pu_width, dst, dst_stride, pu_width, pu_height);
}

void avc_style_luma_interpolation_filter_horizontal_ssse3_intrin(
    EbByte ref_pic,
    uint32_t src_stride,
    EbByte dst,
    uint32_t dst_stride,
    uint32_t pu_width,
    uint32_t pu_height,
    EbByte temp_buf,
    EbBool skip,
    uint32_t frac_pos)
{
    (void)temp_buf;
    __m128i IFOffset, IFCoeff_1_0, IFCoeff_3_2, sum_clip_U8;
    uint32_t width_cnt, height_cnt;
    uint32_t IFShift = 5;

    src_stride <<= skip;
    dst_stride <<= skip;
    pu_height >>= skip;
    frac_pos <<= 5;
    IFOffset = _mm_set1_epi16(0x0010);
    IFCoeff_1_0 = _mm_load_si128((__m128i *)(AvcStyleLumaIFCoeff8_SSSE3 + frac_pos - 32));
    IFCoeff_3_2 = _mm_load_si128((__m128i *)(AvcStyleLumaIFCoeff8_SSSE3 + frac_pos - 16));

    if (!(pu_width & 15)) { // 16x
        __m128i ref0, ref1, ref2, ref3, ref01_lo, ref01_hi, ref23_lo, ref23_hi, sum_lo, sum_hi;

        for (height_cnt = 0; height_cnt < pu_height; ++height_cnt) {
            for (width_cnt = 0; width_cnt < pu_width; width_cnt += 16) {
                ref0 = _mm_loadu_si128((__m128i *)(ref_pic + width_cnt - 1));
                ref1 = _mm_loadu_si128((__m128i *)(ref_pic + width_cnt));
                ref2 = _mm_loadu_si128((__m128i *)(ref_pic + width_cnt + 1));
                ref3 = _mm_loadu_si128((__m128i *)(ref_pic + width_cnt + 2));

                ref01_lo = _mm_unpacklo_epi8(ref0, ref1);
                ref01_hi = _mm_unpackhi_epi8(ref0, ref1);
                ref23_lo = _mm_unpacklo_epi8(ref2, ref3);
                ref23_hi = _mm_unpackhi_epi8(ref2, ref3);

                sum_lo = _mm_srai_epi16(_mm_add_epi16(_mm_add_epi16(_mm_maddubs_epi16(ref01_lo, IFCoeff_1_0), _mm_maddubs_epi16(ref23_lo, IFCoeff_3_2)), IFOffset), IFShift);
                sum_hi = _mm_srai_epi16(_mm_add_epi16(_mm_add_epi16(_mm_maddubs_epi16(ref01_hi, IFCoeff_1_0), _mm_maddubs_epi16(ref23_hi, IFCoeff_3_2)), IFOffset), IFShift);
                sum_clip_U8 = _mm_packus_epi16(sum_lo, sum_hi);
                _mm_storeu_si128((__m128i *)(dst + width_cnt), sum_clip_U8);
            }
            ref_pic += src_stride;
            dst += dst_stride;
        }
        //do the last row if sub-pred ON.
        if (skip) {
            ref_pic -= (src_stride >> 1);
            dst -= (dst_stride >> 1);
            for (width_cnt = 0; width_cnt < pu_width; width_cnt += 16) {
                ref0 = _mm_loadu_si128((__m128i *)(ref_pic + width_cnt - 1));
                ref1 = _mm_loadu_si128((__m128i *)(ref_pic + width_cnt));
                ref2 = _mm_loadu_si128((__m128i *)(ref_pic + width_cnt + 1));
                ref3 = _mm_loadu_si128((__m128i *)(ref_pic + width_cnt + 2));

                ref01_lo = _mm_unpacklo_epi8(ref0, ref1);
                ref01_hi = _mm_unpackhi_epi8(ref0, ref1);
                ref23_lo = _mm_unpacklo_epi8(ref2, ref3);
                ref23_hi = _mm_unpackhi_epi8(ref2, ref3);

                sum_lo = _mm_srai_epi16(_mm_add_epi16(_mm_add_epi16(_mm_maddubs_epi16(ref01_lo, IFCoeff_1_0), _mm_maddubs_epi16(ref23_lo, IFCoeff_3_2)), IFOffset), IFShift);
                sum_hi = _mm_srai_epi16(_mm_add_epi16(_mm_add_epi16(_mm_maddubs_epi16(ref01_hi, IFCoeff_1_0), _mm_maddubs_epi16(ref23_hi, IFCoeff_3_2)), IFOffset), IFShift);
                sum_clip_U8 = _mm_packus_epi16(sum_lo, sum_hi);
                _mm_storeu_si128((__m128i *)(dst + width_cnt), sum_clip_U8);
            }
        }
    }
    else { //8x
        __m128i  sum01, sum23, sum;

        for (height_cnt = 0; height_cnt < pu_height; ++height_cnt) {
            for (width_cnt = 0; width_cnt < pu_width; width_cnt += 8) {
                sum01 = _mm_maddubs_epi16(_mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(ref_pic + width_cnt - 1)),
                    _mm_loadl_epi64((__m128i *)(ref_pic + width_cnt))), IFCoeff_1_0);

                sum23 = _mm_maddubs_epi16(_mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(ref_pic + width_cnt + 1)),
                    _mm_loadl_epi64((__m128i *)(ref_pic + width_cnt + 2))), IFCoeff_3_2);

                sum = _mm_srai_epi16(_mm_add_epi16(_mm_add_epi16(sum01, sum23), IFOffset), IFShift);
                sum_clip_U8 = _mm_packus_epi16(sum, sum);

                _mm_storel_epi64((__m128i *)(dst + width_cnt), sum_clip_U8);
            }
            ref_pic += src_stride;
            dst += dst_stride;
        }

        //do the last row if sub-pred ON.
        if (skip) {
            ref_pic -= (src_stride >> 1);
            dst -= (dst_stride >> 1);
            for (width_cnt = 0; width_cnt < pu_width; width_cnt += 8) {
                sum01 = _mm_maddubs_epi16(_mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(ref_pic + width_cnt - 1)),
                    _mm_loadl_epi64((__m128i *)(ref_pic + width_cnt))), IFCoeff_1_0);

                sum23 = _mm_maddubs_epi16(_mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(ref_pic + width_cnt + 1)),
                    _mm_loadl_epi64((__m128i *)(ref_pic + width_cnt + 2))), IFCoeff_3_2);

                sum = _mm_srai_epi16(_mm_add_epi16(_mm_add_epi16(sum01, sum23), IFOffset), IFShift);
                sum_clip_U8 = _mm_packus_epi16(sum, sum);

                _mm_storel_epi64((__m128i *)(dst + width_cnt), sum_clip_U8);
            }
        }
    }
}

void avc_style_luma_interpolation_filter_vertical_ssse3_intrin(
    EbByte ref_pic,
    uint32_t src_stride,
    EbByte dst,
    uint32_t dst_stride,
    uint32_t pu_width,
    uint32_t pu_height,
    EbByte temp_buf,
    EbBool skip,
    uint32_t frac_pos)
{
    (void)temp_buf;
    __m128i IFOffset, IFCoeff_1_0, IFCoeff_3_2, sum_clip_U8;
    uint32_t width_cnt, height_cnt;
    uint32_t IFShift = 5;
    uint32_t srcStrideSkip = src_stride << (skip ? 1 : 0);
    EbByte refPicTemp, dstTemp;

    frac_pos <<= 5;
    ref_pic -= src_stride;
    IFOffset = _mm_set1_epi16(0x0010);
    IFCoeff_1_0 = _mm_load_si128((__m128i *)(AvcStyleLumaIFCoeff8_SSSE3 + frac_pos - 32));
    IFCoeff_3_2 = _mm_load_si128((__m128i *)(AvcStyleLumaIFCoeff8_SSSE3 + frac_pos - 16));
    dst_stride <<= skip;
    pu_height >>= skip;
    if (!(pu_width & 15)) { //16x

        __m128i sum_lo, sum_hi, ref0, refs, ref2s, ref3s;

        for (width_cnt = 0; width_cnt < pu_width; width_cnt += 16) {
            refPicTemp = ref_pic;
            dstTemp = dst;

            for (height_cnt = 0; height_cnt < pu_height; ++height_cnt) {
                ref0 = _mm_loadu_si128((__m128i *)(refPicTemp));
                refs = _mm_loadu_si128((__m128i *)(refPicTemp + src_stride));
                ref2s = _mm_loadu_si128((__m128i *)(refPicTemp + 2 * src_stride));
                ref3s = _mm_loadu_si128((__m128i *)(refPicTemp + 3 * src_stride));

                sum_lo = _mm_add_epi16(_mm_maddubs_epi16(_mm_unpacklo_epi8(ref0, refs), IFCoeff_1_0),
                    _mm_maddubs_epi16(_mm_unpacklo_epi8(ref2s, ref3s), IFCoeff_3_2));

                sum_hi = _mm_add_epi16(_mm_maddubs_epi16(_mm_unpackhi_epi8(ref0, refs), IFCoeff_1_0),
                    _mm_maddubs_epi16(_mm_unpackhi_epi8(ref2s, ref3s), IFCoeff_3_2));

                sum_lo = _mm_srai_epi16(_mm_add_epi16(sum_lo, IFOffset), IFShift);
                sum_hi = _mm_srai_epi16(_mm_add_epi16(sum_hi, IFOffset), IFShift);
                sum_clip_U8 = _mm_packus_epi16(sum_lo, sum_hi);
                _mm_storeu_si128((__m128i *)(dstTemp), sum_clip_U8);
                dstTemp += dst_stride;
                refPicTemp += srcStrideSkip;
            }
            //do the last row if sub-pred is ON.
            if (skip) {
                dstTemp -= (dst_stride >> 1);
                refPicTemp -= (srcStrideSkip >> 1);
                {
                    ref0 = _mm_loadu_si128((__m128i *)(refPicTemp));
                    refs = _mm_loadu_si128((__m128i *)(refPicTemp + src_stride));
                    ref2s = _mm_loadu_si128((__m128i *)(refPicTemp + 2 * src_stride));
                    ref3s = _mm_loadu_si128((__m128i *)(refPicTemp + 3 * src_stride));

                    sum_lo = _mm_add_epi16(_mm_maddubs_epi16(_mm_unpacklo_epi8(ref0, refs), IFCoeff_1_0),
                        _mm_maddubs_epi16(_mm_unpacklo_epi8(ref2s, ref3s), IFCoeff_3_2));

                    sum_hi = _mm_add_epi16(_mm_maddubs_epi16(_mm_unpackhi_epi8(ref0, refs), IFCoeff_1_0),
                        _mm_maddubs_epi16(_mm_unpackhi_epi8(ref2s, ref3s), IFCoeff_3_2));

                    sum_lo = _mm_srai_epi16(_mm_add_epi16(sum_lo, IFOffset), IFShift);
                    sum_hi = _mm_srai_epi16(_mm_add_epi16(sum_hi, IFOffset), IFShift);
                    sum_clip_U8 = _mm_packus_epi16(sum_lo, sum_hi);
                    _mm_storeu_si128((__m128i *)(dstTemp), sum_clip_U8);
                }
            }
            ref_pic += 16;
            dst += 16;
        }
    }
    else { //8x
        __m128i sum, sum01, sum23;

        for (width_cnt = 0; width_cnt < pu_width; width_cnt += 8) {
            refPicTemp = ref_pic;
            dstTemp = dst;

            for (height_cnt = 0; height_cnt < pu_height; ++height_cnt) {
                sum01 = _mm_maddubs_epi16(_mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(refPicTemp)),
                    _mm_loadl_epi64((__m128i *)(refPicTemp + src_stride))), IFCoeff_1_0);

                sum23 = _mm_maddubs_epi16(_mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(refPicTemp + 2 * src_stride)),
                    _mm_loadl_epi64((__m128i *)(refPicTemp + 3 * src_stride))), IFCoeff_3_2);

                sum = _mm_srai_epi16(_mm_add_epi16(_mm_add_epi16(sum01, sum23), IFOffset), IFShift);
                sum_clip_U8 = _mm_packus_epi16(sum, sum);
                _mm_storel_epi64((__m128i *)(dstTemp), sum_clip_U8);

                dstTemp += dst_stride;
                refPicTemp += srcStrideSkip;
            }
            //do the last row if sub-pred is ON.
            if (skip) {
                dstTemp -= (dst_stride >> 1);
                refPicTemp -= (srcStrideSkip >> 1);
                {
                    sum01 = _mm_maddubs_epi16(_mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(refPicTemp)),
                        _mm_loadl_epi64((__m128i *)(refPicTemp + src_stride))), IFCoeff_1_0);

                    sum23 = _mm_maddubs_epi16(_mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(refPicTemp + 2 * src_stride)),
                        _mm_loadl_epi64((__m128i *)(refPicTemp + 3 * src_stride))), IFCoeff_3_2);

                    sum = _mm_srai_epi16(_mm_add_epi16(_mm_add_epi16(sum01, sum23), IFOffset), IFShift);
                    sum_clip_U8 = _mm_packus_epi16(sum, sum);
                    _mm_storel_epi64((__m128i *)(dstTemp), sum_clip_U8);
                }
            }
            ref_pic += 8;
            dst += 8;
        }
    }
}

void avc_style_luma_interpolation_filter_helper_ssse3(
    EbByte ref_pic,
    uint32_t src_stride,
    EbByte dst,
    uint32_t dst_stride,
    uint32_t pu_width,
    uint32_t pu_height,
    EbByte temp_buf,
    EbBool skip,
    uint32_t frac_pos,
    uint8_t fractional_position)
{

    switch (fractional_position) {
    case 0:
        avc_style_copy_sse2(ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, skip, frac_pos); break;
    case 1:
        avc_style_luma_interpolation_filter_horizontal_ssse3_intrin(ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, skip, frac_pos); break;
    case 2:
        avc_style_luma_interpolation_filter_horizontal_ssse3_intrin(ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, skip, frac_pos); break;
    case 3:
        avc_style_luma_interpolation_filter_horizontal_ssse3_intrin(ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, skip, frac_pos); break;
    case 4:
        avc_style_luma_interpolation_filter_vertical_ssse3_intrin(ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, skip, frac_pos); break;
    case 5:
        avc_style_luma_interpolation_filter_pose_ssse3(ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, skip, frac_pos); break;
    case 6:
        avc_style_luma_interpolation_filter_posf_ssse3(ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, skip, frac_pos); break;
    case 7:
        avc_style_luma_interpolation_filter_posg_ssse3(ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, skip, frac_pos); break;
    case 8:
        avc_style_luma_interpolation_filter_vertical_ssse3_intrin(ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, skip, frac_pos); break;
    case 9:
        avc_style_luma_interpolation_filter_posi_ssse3(ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, skip, frac_pos); break;
    case 10:
        avc_style_luma_interpolation_filter_posj_ssse3(ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, skip, frac_pos); break;
    case 11:
        avc_style_luma_interpolation_filter_posk_ssse3(ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, skip, frac_pos); break;
    case 12:
        avc_style_luma_interpolation_filter_vertical_ssse3_intrin(ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, skip, frac_pos); break;
    case 13:
        avc_style_luma_interpolation_filter_posp_ssse3(ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, skip, frac_pos); break;
    case 14:
        avc_style_luma_interpolation_filter_posq_ssse3(ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, skip, frac_pos); break;
    case 15:
        avc_style_luma_interpolation_filter_posr_ssse3(ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, skip, frac_pos); break;
    default:
        assert(0);
    }
}
