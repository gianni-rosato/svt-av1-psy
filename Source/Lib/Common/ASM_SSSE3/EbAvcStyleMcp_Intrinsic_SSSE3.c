/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbAvcStyleMcp_SSSE3.h"

#include "EbMcp_SSE2.h"
#include "EbDefinitions.h"
#include "EbAvcStyleMcp_SSE2.h"
#include "common_dsp_rtcd.h"
#include <emmintrin.h>
#include <tmmintrin.h>
#if !REMOVE_UNUSED_CODE
EB_EXTERN EB_ALIGN(16) const int8_t avc_style_luma_if_coeff8_ssse3[] = {
    -1, 25, -1, 25, -1, 25, -1, 25, -1, 25, -1, 25, -1, 25, -1, 25, 9,  -1, 9,  -1, 9,  -1, 9,  -1,
    9,  -1, 9,  -1, 9,  -1, 9,  -1, -2, 18, -2, 18, -2, 18, -2, 18, -2, 18, -2, 18, -2, 18, -2, 18,
    18, -2, 18, -2, 18, -2, 18, -2, 18, -2, 18, -2, 18, -2, 18, -2, -1, 9,  -1, 9,  -1, 9,  -1, 9,
    -1, 9,  -1, 9,  -1, 9,  -1, 9,  25, -1, 25, -1, 25, -1, 25, -1, 25, -1, 25, -1, 25, -1, 25, -1};

void avc_style_luma_interpolation_filter_pose_ssse3(EbByte ref_pic, uint32_t src_stride, EbByte dst,
                                                    uint32_t dst_stride, uint32_t pu_width,
                                                    uint32_t pu_height, EbByte temp_buf,
                                                    uint32_t frac_pos) {
    uint32_t temp_buf_size = pu_width * pu_height;
    (void)frac_pos;
    avc_style_luma_interpolation_filter_horizontal_ssse3_intrin(
        ref_pic, src_stride, temp_buf, pu_width, pu_width, pu_height, 0, 2);
    avc_style_luma_interpolation_filter_vertical_ssse3_intrin(
        ref_pic, src_stride, temp_buf + temp_buf_size, pu_width, pu_width, pu_height, 0, 2);
    picture_average_kernel_sse2_intrin(temp_buf,
                                       pu_width,
                                       temp_buf + temp_buf_size,
                                       pu_width,
                                       dst,
                                       dst_stride,
                                       pu_width,
                                       pu_height);
}

void avc_style_luma_interpolation_filter_posf_ssse3(EbByte ref_pic, uint32_t src_stride, EbByte dst,
                                                    uint32_t dst_stride, uint32_t pu_width,
                                                    uint32_t pu_height, EbByte temp_buf,
                                                    uint32_t frac_pos) {
    uint32_t temp_buf_size = pu_width * pu_height;
    (void)frac_pos;
    avc_style_luma_interpolation_filter_horizontal_ssse3_intrin(
        ref_pic - src_stride,
        src_stride,
        temp_buf + temp_buf_size,
        pu_width,
        pu_width,
        (pu_height + 3),
        0,
        2);
    avc_style_luma_interpolation_filter_vertical_ssse3_intrin(temp_buf + temp_buf_size + pu_width,
                                                              pu_width,
                                                              temp_buf,
                                                              pu_width,
                                                              pu_width,
                                                              pu_height,
                                                              0,
                                                              2);
    picture_average_kernel_sse2_intrin(temp_buf + temp_buf_size + pu_width,
                                       pu_width,
                                       temp_buf,
                                       pu_width,
                                       dst,
                                       dst_stride,
                                       pu_width,
                                       pu_height);
}

void avc_style_luma_interpolation_filter_posg_ssse3(EbByte ref_pic, uint32_t src_stride, EbByte dst,
                                                    uint32_t dst_stride, uint32_t pu_width,
                                                    uint32_t pu_height, EbByte temp_buf,
                                                    uint32_t frac_pos) {
    uint32_t temp_buf_size = pu_width * pu_height;
    (void)frac_pos;
    avc_style_luma_interpolation_filter_horizontal_ssse3_intrin(
        ref_pic, src_stride, temp_buf, pu_width, pu_width, pu_height, 0, 2);
    avc_style_luma_interpolation_filter_vertical_ssse3_intrin(ref_pic + 1,
                                                              src_stride,
                                                              temp_buf + temp_buf_size,
                                                              pu_width,
                                                              pu_width,
                                                              pu_height,
                                                              0,
                                                              2);
    picture_average_kernel_sse2_intrin(temp_buf,
                                       pu_width,
                                       temp_buf + temp_buf_size,
                                       pu_width,
                                       dst,
                                       dst_stride,
                                       pu_width,
                                       pu_height);
}

void avc_style_luma_interpolation_filter_posi_ssse3(EbByte ref_pic, uint32_t src_stride, EbByte dst,
                                                    uint32_t dst_stride, uint32_t pu_width,
                                                    uint32_t pu_height, EbByte temp_buf,
                                                    uint32_t frac_pos) {
    uint32_t temp_buf_size = pu_width * pu_height;
    (void)frac_pos;
    avc_style_luma_interpolation_filter_vertical_ssse3_intrin(
        ref_pic, src_stride, temp_buf, pu_width, pu_width, pu_height, 0, 2);
    avc_style_luma_interpolation_filter_posj_ssse3(ref_pic,
                                                   src_stride,
                                                   temp_buf + temp_buf_size,
                                                   pu_width,
                                                   pu_width,
                                                   pu_height,
                                                   temp_buf + 2 * temp_buf_size,
                                                   2);
    picture_average_kernel_sse2_intrin(temp_buf,
                                       pu_width,
                                       temp_buf + temp_buf_size,
                                       pu_width,
                                       dst,
                                       dst_stride,
                                       pu_width,
                                       pu_height);
}

void avc_style_luma_interpolation_filter_posj_ssse3(EbByte ref_pic, uint32_t src_stride, EbByte dst,
                                                    uint32_t dst_stride, uint32_t pu_width,
                                                    uint32_t pu_height, EbByte temp_buf,
                                                    uint32_t frac_pos) {
    (void)frac_pos;
    avc_style_luma_interpolation_filter_horizontal_ssse3_intrin(
        ref_pic - src_stride,
        src_stride,
        temp_buf,
        pu_width,
        pu_width,
        (pu_height + 3),
        0,
        2);

    avc_style_luma_interpolation_filter_vertical_ssse3_intrin(
        temp_buf + pu_width, pu_width, dst, dst_stride, pu_width, pu_height, 0, 2);
}

void avc_style_luma_interpolation_filter_posk_ssse3(EbByte ref_pic, uint32_t src_stride, EbByte dst,
                                                    uint32_t dst_stride, uint32_t pu_width,
                                                    uint32_t pu_height, EbByte temp_buf,
                                                    uint32_t frac_pos) {
    uint32_t temp_buf_size = pu_width * pu_height;
    (void)frac_pos;
    avc_style_luma_interpolation_filter_vertical_ssse3_intrin(
        ref_pic + 1, src_stride, temp_buf, pu_width, pu_width, pu_height, 0, 2);
    avc_style_luma_interpolation_filter_posj_ssse3(ref_pic,
                                                   src_stride,
                                                   temp_buf + temp_buf_size,
                                                   pu_width,
                                                   pu_width,
                                                   pu_height,
                                                   temp_buf + 2 * temp_buf_size,
                                                   2);
    picture_average_kernel_sse2_intrin(temp_buf,
                                       pu_width,
                                       temp_buf + temp_buf_size,
                                       pu_width,
                                       dst,
                                       dst_stride,
                                       pu_width,
                                       pu_height);
}

void avc_style_luma_interpolation_filter_posp_ssse3(EbByte ref_pic, uint32_t src_stride, EbByte dst,
                                                    uint32_t dst_stride, uint32_t pu_width,
                                                    uint32_t pu_height, EbByte temp_buf,
                                                    uint32_t frac_pos) {
    uint32_t temp_buf_size = pu_width * pu_height;
    (void)frac_pos;
    avc_style_luma_interpolation_filter_vertical_ssse3_intrin(
        ref_pic, src_stride, temp_buf, pu_width, pu_width, pu_height, 0, 2);
    avc_style_luma_interpolation_filter_horizontal_ssse3_intrin(ref_pic + src_stride,
                                                                src_stride,
                                                                temp_buf + temp_buf_size,
                                                                pu_width,
                                                                pu_width,
                                                                pu_height,
                                                                0,
                                                                2);
    picture_average_kernel_sse2_intrin(temp_buf,
                                       pu_width,
                                       temp_buf + temp_buf_size,
                                       pu_width,
                                       dst,
                                       dst_stride,
                                       pu_width,
                                       pu_height);
}

void avc_style_luma_interpolation_filter_posq_ssse3(EbByte ref_pic, uint32_t src_stride, EbByte dst,
                                                    uint32_t dst_stride, uint32_t pu_width,
                                                    uint32_t pu_height, EbByte temp_buf,
                                                    uint32_t frac_pos) {
    uint32_t temp_buf_size = pu_width * pu_height;
    (void)frac_pos;
    avc_style_luma_interpolation_filter_horizontal_ssse3_intrin(
        ref_pic - src_stride,
        src_stride,
        temp_buf + temp_buf_size,
        pu_width,
        pu_width,
        (pu_height + 3),
        0,
        2);
    avc_style_luma_interpolation_filter_vertical_ssse3_intrin(temp_buf + temp_buf_size + pu_width,
                                                              pu_width,
                                                              temp_buf,
                                                              pu_width,
                                                              pu_width,
                                                              pu_height,
                                                              0,
                                                              2);
    picture_average_kernel_sse2_intrin(temp_buf + temp_buf_size + 2 * pu_width,
                                       pu_width,
                                       temp_buf,
                                       pu_width,
                                       dst,
                                       dst_stride,
                                       pu_width,
                                       pu_height);
}

void avc_style_luma_interpolation_filter_posr_ssse3(EbByte ref_pic, uint32_t src_stride, EbByte dst,
                                                    uint32_t dst_stride, uint32_t pu_width,
                                                    uint32_t pu_height, EbByte temp_buf,
                                                    uint32_t frac_pos) {
    uint32_t temp_buf_size = pu_width * pu_height;
    (void)frac_pos;
    avc_style_luma_interpolation_filter_vertical_ssse3_intrin(
        ref_pic + 1, src_stride, temp_buf, pu_width, pu_width, pu_height, 0, 2);
    avc_style_luma_interpolation_filter_horizontal_ssse3_intrin(ref_pic + src_stride,
                                                                src_stride,
                                                                temp_buf + temp_buf_size,
                                                                pu_width,
                                                                pu_width,
                                                                pu_height,
                                                                0,
                                                                2);
    picture_average_kernel_sse2_intrin(temp_buf,
                                       pu_width,
                                       temp_buf + temp_buf_size,
                                       pu_width,
                                       dst,
                                       dst_stride,
                                       pu_width,
                                       pu_height);
}

void avc_style_luma_interpolation_filter_horizontal_ssse3_intrin(
    EbByte ref_pic, uint32_t src_stride, EbByte dst, uint32_t dst_stride, uint32_t pu_width,
    uint32_t pu_height, EbByte temp_buf, uint32_t frac_pos) {
    (void)temp_buf;
    __m128i  if_offset, if_coeff_1_0, if_coeff_3_2, sum_clip_u8;
    uint32_t width_cnt, height_cnt;
    uint32_t if_shift = 5;

    frac_pos <<= 5;
    if_offset    = _mm_set1_epi16(0x0010);
    if_coeff_1_0 = _mm_loadu_si128((__m128i *)(avc_style_luma_if_coeff8_ssse3 + frac_pos - 32));
    if_coeff_3_2 = _mm_loadu_si128((__m128i *)(avc_style_luma_if_coeff8_ssse3 + frac_pos - 16));

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

                sum_lo = _mm_srai_epi16(
                    _mm_add_epi16(_mm_add_epi16(_mm_maddubs_epi16(ref01_lo, if_coeff_1_0),
                                                _mm_maddubs_epi16(ref23_lo, if_coeff_3_2)),
                                  if_offset),
                    if_shift);
                sum_hi = _mm_srai_epi16(
                    _mm_add_epi16(_mm_add_epi16(_mm_maddubs_epi16(ref01_hi, if_coeff_1_0),
                                                _mm_maddubs_epi16(ref23_hi, if_coeff_3_2)),
                                  if_offset),
                    if_shift);
                sum_clip_u8 = _mm_packus_epi16(sum_lo, sum_hi);
                _mm_storeu_si128((__m128i *)(dst + width_cnt), sum_clip_u8);
            }
            ref_pic += src_stride;
            dst += dst_stride;
        }
    } else { //8x
        __m128i sum01, sum23, sum;

        for (height_cnt = 0; height_cnt < pu_height; ++height_cnt) {
            for (width_cnt = 0; width_cnt < pu_width; width_cnt += 8) {
                sum01 = _mm_maddubs_epi16(
                    _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(ref_pic + width_cnt - 1)),
                                      _mm_loadl_epi64((__m128i *)(ref_pic + width_cnt))),
                    if_coeff_1_0);

                sum23 = _mm_maddubs_epi16(
                    _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(ref_pic + width_cnt + 1)),
                                      _mm_loadl_epi64((__m128i *)(ref_pic + width_cnt + 2))),
                    if_coeff_3_2);

                sum =
                    _mm_srai_epi16(_mm_add_epi16(_mm_add_epi16(sum01, sum23), if_offset), if_shift);
                sum_clip_u8 = _mm_packus_epi16(sum, sum);

                _mm_storel_epi64((__m128i *)(dst + width_cnt), sum_clip_u8);
            }
            ref_pic += src_stride;
            dst += dst_stride;
        }
    }
}

void avc_style_luma_interpolation_filter_vertical_ssse3_intrin(EbByte ref_pic, uint32_t src_stride,
                                                               EbByte dst, uint32_t dst_stride,
                                                               uint32_t pu_width,
                                                               uint32_t pu_height, EbByte temp_buf,
                                                               uint32_t frac_pos) {
    (void)temp_buf;
    __m128i  if_offset, if_coeff_1_0, if_coeff_3_2, sum_clip_u8;
    uint32_t width_cnt, height_cnt;
    uint32_t if_shift        = 5;

    EbByte   ref_pic_temp, dst_temp;

    frac_pos <<= 5;
    ref_pic -= src_stride;
    if_offset    = _mm_set1_epi16(0x0010);
    if_coeff_1_0 = _mm_loadu_si128((__m128i *)(avc_style_luma_if_coeff8_ssse3 + frac_pos - 32));
    if_coeff_3_2 = _mm_loadu_si128((__m128i *)(avc_style_luma_if_coeff8_ssse3 + frac_pos - 16));

    if (!(pu_width & 15)) { //16x

        __m128i sum_lo, sum_hi, ref0, refs, ref2s, ref3s;

        for (width_cnt = 0; width_cnt < pu_width; width_cnt += 16) {
            ref_pic_temp = ref_pic;
            dst_temp     = dst;

            for (height_cnt = 0; height_cnt < pu_height; ++height_cnt) {
                ref0  = _mm_loadu_si128((__m128i *)(ref_pic_temp));
                refs  = _mm_loadu_si128((__m128i *)(ref_pic_temp + src_stride));
                ref2s = _mm_loadu_si128((__m128i *)(ref_pic_temp + 2 * src_stride));
                ref3s = _mm_loadu_si128((__m128i *)(ref_pic_temp + 3 * src_stride));

                sum_lo =
                    _mm_add_epi16(_mm_maddubs_epi16(_mm_unpacklo_epi8(ref0, refs), if_coeff_1_0),
                                  _mm_maddubs_epi16(_mm_unpacklo_epi8(ref2s, ref3s), if_coeff_3_2));

                sum_hi =
                    _mm_add_epi16(_mm_maddubs_epi16(_mm_unpackhi_epi8(ref0, refs), if_coeff_1_0),
                                  _mm_maddubs_epi16(_mm_unpackhi_epi8(ref2s, ref3s), if_coeff_3_2));

                sum_lo      = _mm_srai_epi16(_mm_add_epi16(sum_lo, if_offset), if_shift);
                sum_hi      = _mm_srai_epi16(_mm_add_epi16(sum_hi, if_offset), if_shift);
                sum_clip_u8 = _mm_packus_epi16(sum_lo, sum_hi);
                _mm_storeu_si128((__m128i *)(dst_temp), sum_clip_u8);
                dst_temp += dst_stride;
                ref_pic_temp += src_stride;
            }
            ref_pic += 16;
            dst += 16;
        }
    } else { //8x
        __m128i sum, sum01, sum23;

        for (width_cnt = 0; width_cnt < pu_width; width_cnt += 8) {
            ref_pic_temp = ref_pic;
            dst_temp     = dst;

            for (height_cnt = 0; height_cnt < pu_height; ++height_cnt) {
                sum01 = _mm_maddubs_epi16(
                    _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(ref_pic_temp)),
                                      _mm_loadl_epi64((__m128i *)(ref_pic_temp + src_stride))),
                    if_coeff_1_0);

                sum23 = _mm_maddubs_epi16(
                    _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(ref_pic_temp + 2 * src_stride)),
                                      _mm_loadl_epi64((__m128i *)(ref_pic_temp + 3 * src_stride))),
                    if_coeff_3_2);

                sum =
                    _mm_srai_epi16(_mm_add_epi16(_mm_add_epi16(sum01, sum23), if_offset), if_shift);
                sum_clip_u8 = _mm_packus_epi16(sum, sum);
                _mm_storel_epi64((__m128i *)(dst_temp), sum_clip_u8);

                dst_temp += dst_stride;
                ref_pic_temp += src_stride;
            }
            ref_pic += 8;
            dst += 8;
        }
    }
}

void avc_style_luma_interpolation_filter_helper_ssse3(EbByte ref_pic, uint32_t src_stride,
                                                      EbByte dst, uint32_t dst_stride,
                                                      uint32_t pu_width, uint32_t pu_height,
                                                      EbByte temp_buf,
                                                      uint32_t frac_pos,
                                                      uint8_t  fractional_position) {
    switch (fractional_position) {
    case 0:
        avc_style_copy_sse2(
            ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos);
        break;
    case 1:
        avc_style_luma_interpolation_filter_horizontal_ssse3_intrin(
            ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos);
        break;
    case 2:
        avc_style_luma_interpolation_filter_horizontal_ssse3_intrin(
            ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos);
        break;
    case 3:
        avc_style_luma_interpolation_filter_horizontal_ssse3_intrin(
            ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos);
        break;
    case 4:
        avc_style_luma_interpolation_filter_vertical_ssse3_intrin(
            ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos);
        break;
    case 5:
        avc_style_luma_interpolation_filter_pose_ssse3(
            ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos);
        break;
    case 6:
        avc_style_luma_interpolation_filter_posf_ssse3(
            ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos);
        break;
    case 7:
        avc_style_luma_interpolation_filter_posg_ssse3(
            ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos);
        break;
    case 8:
        avc_style_luma_interpolation_filter_vertical_ssse3_intrin(
            ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos);
        break;
    case 9:
        avc_style_luma_interpolation_filter_posi_ssse3(
            ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos);
        break;
    case 10:
        avc_style_luma_interpolation_filter_posj_ssse3(
            ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos);
        break;
    case 11:
        avc_style_luma_interpolation_filter_posk_ssse3(
            ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos);
        break;
    case 12:
        avc_style_luma_interpolation_filter_vertical_ssse3_intrin(
            ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos);
        break;
    case 13:
        avc_style_luma_interpolation_filter_posp_ssse3(
            ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos);
        break;
    case 14:
        avc_style_luma_interpolation_filter_posq_ssse3(
            ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos);
        break;
    case 15:
        avc_style_luma_interpolation_filter_posr_ssse3(
            ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos);
        break;
    default: assert(0);
    }
}
#endif
