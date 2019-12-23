/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbAvcStyleMcp.h"
#include "EbPictureOperators.h"
#include "EbUtility.h"

static const int8_t avc_style_luma_if_coeff[4][4] = {
    {0, 0, 0, 0},
    {-1, 25, 9, -1},
    {-2, 18, 18, -2},
    {-1, 9, 25, -1},
};

void avc_style_copy_c(EbByte ref_pic, uint32_t src_stride, EbByte dst, uint32_t dst_stride,
                      uint32_t pu_width, uint32_t pu_height, EbByte temp_buf, uint32_t frac_pos) {
    (void)temp_buf;
    (void)frac_pos;
    picture_copy_kernel(ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, 1);
}

void avc_style_luma_interpolation_filter_horizontal_c(EbByte ref_pic, uint32_t src_stride,
                                                      EbByte dst, uint32_t dst_stride,
                                                      uint32_t pu_width, uint32_t pu_height,
                                                      EbByte temp_buf, uint32_t frac_pos) {
    const int8_t *if_coeff           = avc_style_luma_if_coeff[frac_pos];
    const int32_t if_init_pos_offset = -1;
    const uint8_t if_shift           = 5;
    const int16_t if_offset          = POW2(if_shift - 1);
    uint32_t      x, y;
    (void)temp_buf;

    ref_pic += if_init_pos_offset;
    for (y = 0; y < pu_height; y++) {
        for (x = 0; x < pu_width; x++) {
            dst[x] = (uint8_t)CLIP3(
                0,
                MAX_SAMPLE_VALUE,
                (ref_pic[x] * if_coeff[0] + ref_pic[x + 1] * if_coeff[1] +
                 ref_pic[x + 2] * if_coeff[2] + ref_pic[x + 3] * if_coeff[3] + if_offset) >>
                    if_shift);
        }
        ref_pic += src_stride;
        dst += dst_stride;
    }
}

void avc_style_luma_interpolation_filter_vertical_c(EbByte ref_pic, uint32_t src_stride, EbByte dst,
                                                    uint32_t dst_stride, uint32_t pu_width,
                                                    uint32_t pu_height, EbByte temp_buf,
                                                    uint32_t frac_pos) {
    const int8_t *if_coeff           = avc_style_luma_if_coeff[frac_pos];
    const int32_t if_stride          = src_stride;
    const int32_t if_init_pos_offset = -(int32_t)src_stride;
    const uint8_t if_shift           = 5;
    const int16_t if_offset          = POW2(if_shift - 1);
    uint32_t      x, y;
    (void)temp_buf;

    ref_pic += if_init_pos_offset;
    for (y = 0; y < pu_height; y++) {
        for (x = 0; x < pu_width; x++) {
            dst[x] =
                (uint8_t)CLIP3(0,
                               MAX_SAMPLE_VALUE,
                               (ref_pic[x] * if_coeff[0] + ref_pic[x + if_stride] * if_coeff[1] +
                                ref_pic[x + 2 * if_stride] * if_coeff[2] +
                                ref_pic[x + 3 * if_stride] * if_coeff[3] + if_offset) >>
                                   if_shift);
        }
        ref_pic += src_stride;
        dst += dst_stride;
    }
}

void avc_style_luma_interpolation_filter_pose_c(EbByte ref_pic, uint32_t src_stride, EbByte dst,
                                                uint32_t dst_stride, uint32_t pu_width,
                                                uint32_t pu_height, EbByte temp_buf,
                                                uint32_t frac_pos) {
    uint32_t temp_buf_size = pu_width * pu_height;
    (void)frac_pos;
    avc_style_luma_interpolation_filter_horizontal_c(
        ref_pic, src_stride, temp_buf, pu_width, pu_width, pu_height, 0, 2);
    avc_style_luma_interpolation_filter_vertical_c(
        ref_pic, src_stride, temp_buf + temp_buf_size, pu_width, pu_width, pu_height, 0, 2);
    picture_average_kernel_c(temp_buf,
                             pu_width,
                             temp_buf + temp_buf_size,
                             pu_width,
                             dst,
                             dst_stride,
                             pu_width,
                             pu_height);
}

void avc_style_luma_interpolation_filter_posf_c(EbByte ref_pic, uint32_t src_stride, EbByte dst,
                                                uint32_t dst_stride, uint32_t pu_width,
                                                uint32_t pu_height, EbByte temp_buf,
                                                uint32_t frac_pos) {
    uint32_t temp_buf_size = pu_width * pu_height;
    (void)frac_pos;
    avc_style_luma_interpolation_filter_horizontal_c(
        ref_pic, src_stride, temp_buf, pu_width, pu_width, pu_height, 0, 2);
    avc_style_luma_interpolation_filter_posj_c(ref_pic,
                                               src_stride,
                                               temp_buf + temp_buf_size,
                                               pu_width,
                                               pu_width,
                                               pu_height,
                                               temp_buf + 2 * temp_buf_size,
                                               2);
    picture_average_kernel_c(temp_buf,
                             pu_width,
                             temp_buf + temp_buf_size,
                             pu_width,
                             dst,
                             dst_stride,
                             pu_width,
                             pu_height);
}

void avc_style_luma_interpolation_filter_posg_c(EbByte ref_pic, uint32_t src_stride, EbByte dst,
                                                uint32_t dst_stride, uint32_t pu_width,
                                                uint32_t pu_height, EbByte temp_buf,
                                                uint32_t frac_pos) {
    uint32_t temp_buf_size = pu_width * pu_height;
    (void)frac_pos;
    avc_style_luma_interpolation_filter_horizontal_c(
        ref_pic, src_stride, temp_buf, pu_width, pu_width, pu_height, 0, 2);
    avc_style_luma_interpolation_filter_vertical_c(
        ref_pic + 1, src_stride, temp_buf + temp_buf_size, pu_width, pu_width, pu_height, 0, 2);
    picture_average_kernel_c(temp_buf,
                             pu_width,
                             temp_buf + temp_buf_size,
                             pu_width,
                             dst,
                             dst_stride,
                             pu_width,
                             pu_height);
}

void avc_style_luma_interpolation_filter_posi_c(EbByte ref_pic, uint32_t src_stride, EbByte dst,
                                                uint32_t dst_stride, uint32_t pu_width,
                                                uint32_t pu_height, EbByte temp_buf,
                                                uint32_t frac_pos) {
    uint32_t temp_buf_size = pu_width * pu_height;
    (void)frac_pos;
    avc_style_luma_interpolation_filter_vertical_c(
        ref_pic, src_stride, temp_buf, pu_width, pu_width, pu_height, 0, 2);
    avc_style_luma_interpolation_filter_posj_c(ref_pic,
                                               src_stride,
                                               temp_buf + temp_buf_size,
                                               pu_width,
                                               pu_width,
                                               pu_height,
                                               temp_buf + 2 * temp_buf_size,
                                               2);
    picture_average_kernel_c(temp_buf,
                             pu_width,
                             temp_buf + temp_buf_size,
                             pu_width,
                             dst,
                             dst_stride,
                             pu_width,
                             pu_height);
}

void avc_style_luma_interpolation_filter_posj_c(EbByte ref_pic, uint32_t src_stride, EbByte dst,
                                                uint32_t dst_stride, uint32_t pu_width,
                                                uint32_t pu_height, EbByte temp_buf,
                                                uint32_t frac_pos) {
    (void)frac_pos;
    avc_style_luma_interpolation_filter_horizontal_c(
        ref_pic - src_stride, src_stride, temp_buf, pu_width, pu_width, (pu_height + 4), 0, 2);
    avc_style_luma_interpolation_filter_vertical_c(
        temp_buf + pu_width, pu_width, dst, dst_stride, pu_width, pu_height, 0, 2);
}

void avc_style_luma_interpolation_filter_posk_c(EbByte ref_pic, uint32_t src_stride, EbByte dst,
                                                uint32_t dst_stride, uint32_t pu_width,
                                                uint32_t pu_height, EbByte temp_buf,
                                                uint32_t frac_pos) {
    uint32_t temp_buf_size = pu_width * pu_height;
    (void)frac_pos;
    avc_style_luma_interpolation_filter_vertical_c(
        ref_pic + 1, src_stride, temp_buf, pu_width, pu_width, pu_height, 0, 2);
    avc_style_luma_interpolation_filter_posj_c(ref_pic,
                                               src_stride,
                                               temp_buf + temp_buf_size,
                                               pu_width,
                                               pu_width,
                                               pu_height,
                                               temp_buf + 2 * temp_buf_size,
                                               2);
    picture_average_kernel_c(temp_buf,
                             pu_width,
                             temp_buf + temp_buf_size,
                             pu_width,
                             dst,
                             dst_stride,
                             pu_width,
                             pu_height);
}

void avc_style_luma_interpolation_filter_posp_c(EbByte ref_pic, uint32_t src_stride, EbByte dst,
                                                uint32_t dst_stride, uint32_t pu_width,
                                                uint32_t pu_height, EbByte temp_buf,
                                                uint32_t frac_pos) {
    uint32_t temp_buf_size = pu_width * pu_height;
    (void)frac_pos;
    avc_style_luma_interpolation_filter_vertical_c(
        ref_pic, src_stride, temp_buf, pu_width, pu_width, pu_height, 0, 2);
    avc_style_luma_interpolation_filter_horizontal_c(ref_pic + src_stride,
                                                     src_stride,
                                                     temp_buf + temp_buf_size,
                                                     pu_width,
                                                     pu_width,
                                                     pu_height,
                                                     0,
                                                     2);
    picture_average_kernel_c(temp_buf,
                             pu_width,
                             temp_buf + temp_buf_size,
                             pu_width,
                             dst,
                             dst_stride,
                             pu_width,
                             pu_height);
}

void avc_style_luma_interpolation_filter_posq_c(EbByte ref_pic, uint32_t src_stride, EbByte dst,
                                                uint32_t dst_stride, uint32_t pu_width,
                                                uint32_t pu_height, EbByte temp_buf,
                                                uint32_t frac_pos) {
    uint32_t temp_buf_size = pu_width * pu_height;
    (void)frac_pos;
    avc_style_luma_interpolation_filter_horizontal_c(
        ref_pic + src_stride, src_stride, temp_buf, pu_width, pu_width, pu_height, 0, 2);
    avc_style_luma_interpolation_filter_posj_c(ref_pic,
                                               src_stride,
                                               temp_buf + temp_buf_size,
                                               pu_width,
                                               pu_width,
                                               pu_height,
                                               temp_buf + 2 * temp_buf_size,
                                               2);
    picture_average_kernel_c(temp_buf,
                             pu_width,
                             temp_buf + temp_buf_size,
                             pu_width,
                             dst,
                             dst_stride,
                             pu_width,
                             pu_height);
}

void avc_style_luma_interpolation_filter_posr_c(EbByte ref_pic, uint32_t src_stride, EbByte dst,
                                                uint32_t dst_stride, uint32_t pu_width,
                                                uint32_t pu_height, EbByte temp_buf,
                                                uint32_t frac_pos) {
    uint32_t temp_buf_size = pu_width * pu_height;
    (void)frac_pos;
    avc_style_luma_interpolation_filter_vertical_c(
        ref_pic + 1, src_stride, temp_buf, pu_width, pu_width, pu_height, 0, 2);
    avc_style_luma_interpolation_filter_horizontal_c(ref_pic + src_stride,
                                                     src_stride,
                                                     temp_buf + temp_buf_size,
                                                     pu_width,
                                                     pu_width,
                                                     pu_height,
                                                     0,
                                                     2);
    picture_average_kernel_c(temp_buf,
                             pu_width,
                             temp_buf + temp_buf_size,
                             pu_width,
                             dst,
                             dst_stride,
                             pu_width,
                             pu_height);
}

void avc_style_luma_interpolation_filter_helper_c(EbByte ref_pic, uint32_t src_stride, EbByte dst,
                                                  uint32_t dst_stride, uint32_t pu_width,
                                                  uint32_t pu_height, EbByte temp_buf, EbBool skip,
                                                  uint32_t frac_pos, uint8_t fractional_position) {
    /* Code with 'skip' true are not used. Cleanup should remove 'skip' parameter. */
    /* frac_pos and fractional_position are redundant as well, cleanup should also unify the two*/
    (void)skip;
    assert(!skip);

    switch (fractional_position) {
    case 0:
        avc_style_copy_c(
            ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos);
        break;
    case 1:
        avc_style_luma_interpolation_filter_horizontal_c(
            ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos);
        break;
    case 2:
        avc_style_luma_interpolation_filter_horizontal_c(
            ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos);
        break;
    case 3:
        avc_style_luma_interpolation_filter_horizontal_c(
            ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos);
        break;
    case 4:
        avc_style_luma_interpolation_filter_vertical_c(
            ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos);
        break;
    case 5:
        avc_style_luma_interpolation_filter_pose_c(
            ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos);
        break;
    case 6:
        avc_style_luma_interpolation_filter_posf_c(
            ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos);
        break;
    case 7:
        avc_style_luma_interpolation_filter_posg_c(
            ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos);
        break;
    case 8:
        avc_style_luma_interpolation_filter_vertical_c(
            ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos);
        break;
    case 9:
        avc_style_luma_interpolation_filter_posi_c(
            ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos);
        break;
    case 10:
        avc_style_luma_interpolation_filter_posj_c(
            ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos);
        break;
    case 11:
        avc_style_luma_interpolation_filter_posk_c(
            ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos);
        break;
    case 12:
        avc_style_luma_interpolation_filter_vertical_c(
            ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos);
        break;
    case 13:
        avc_style_luma_interpolation_filter_posp_c(
            ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos);
        break;
    case 14:
        avc_style_luma_interpolation_filter_posq_c(
            ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos);
        break;
    case 15:
        avc_style_luma_interpolation_filter_posr_c(
            ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos);
        break;
    default: assert(0);
    }
}
