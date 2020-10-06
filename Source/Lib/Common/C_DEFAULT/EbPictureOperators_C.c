/*
* Copyright(c) 2019 Intel Corporation
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

#include "EbPictureOperators_C.h"
#include "EbUtility.h"
#include "common_dsp_rtcd.h"
/*********************************
* Picture Average
*********************************/
void svt_picture_average_kernel_c(EbByte src0, uint32_t src0_stride, EbByte src1,
                                  uint32_t src1_stride, EbByte dst, uint32_t dst_stride,
                                  uint32_t area_width, uint32_t area_height) {
    uint32_t x, y;

    for (y = 0; y < area_height; y++) {
        for (x = 0; x < area_width; x++) { dst[x] = (src0[x] + src1[x] + 1) >> 1; }
        src0 += src0_stride;
        src1 += src1_stride;
        dst += dst_stride;
    }
}

void svt_picture_average_kernel1_line_c(EbByte src0, EbByte src1, EbByte dst, uint32_t areaWidth) {
    uint32_t i;
    for (i = 0; i < areaWidth; i++) dst[i] = (src0[i] + src1[i] + 1) / 2;
}

/*********************************
* Picture Copy Kernel
*********************************/
void eb_memcpy_c(void  *dst_ptr, void  const*src_ptr, size_t size)
{
    memcpy(dst_ptr, src_ptr, size);
}
void picture_copy_kernel(EbByte src, uint32_t src_stride, EbByte dst, uint32_t dst_stride,
                         uint32_t area_width, uint32_t area_height,
                         uint32_t bytes_per_sample) //=1 always)
{
    uint32_t       sample_count       = 0;
    const uint32_t sample_total_count = area_width * area_height;
    const uint32_t copy_length        = area_width * bytes_per_sample;

    src_stride *= bytes_per_sample;
    dst_stride *= bytes_per_sample;

    while (sample_count < sample_total_count) {
        eb_memcpy_c(dst, src, copy_length);
        src += src_stride;
        dst += dst_stride;
        sample_count += area_width;
    }

    return;
}

// C equivalents

uint64_t svt_spatial_full_distortion_kernel_c(uint8_t *input, uint32_t input_offset,
                                              uint32_t input_stride, uint8_t *recon,
                                              int32_t recon_offset, uint32_t recon_stride,
                                              uint32_t area_width, uint32_t area_height) {
    uint64_t spatial_distortion = 0;
    input += input_offset;
    recon += recon_offset;

    for (uint32_t row_index = 0; row_index < area_height; ++row_index) {
        uint32_t column_index = 0;
        while (column_index < area_width) {
            spatial_distortion +=
                (int64_t)SQR((int64_t)(input[column_index]) - (recon[column_index]));
            ++column_index;
        }

        input += input_stride;
        recon += recon_stride;
    }
    return spatial_distortion;
}
