/*
* Copyright(c) 2024 Intel Corporation
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

#include "definitions.h"
#include <arm_neon.h>
#include <stdint.h>

void svt_enc_msb_un_pack2d_neon(uint16_t *in16_bit_buffer, uint32_t in_stride, uint8_t *out8_bit_buffer,
                                uint8_t *outn_bit_buffer, uint32_t out8_stride, uint32_t outn_stride, uint32_t width,
                                uint32_t height) {
    uint32_t   x, y;
    uint8x16_t temp_pixel0_u8, temp_pixel1_u8;

    const uint16x8_t xmm_3    = vdupq_n_u16(0x0003);
    const uint16x8_t xmm_00ff = vdupq_n_u16(0x00FF);

    if (width == 4) {
        for (y = 0; y < height; y += 2) {
            const uint16x8_t in_pixel0 = vcombine_u16(vld1_u16(in16_bit_buffer), vdup_n_u16(0));
            const uint16x8_t in_pixel1 = vcombine_u16(vld1_u16(in16_bit_buffer + in_stride), vdup_n_u16(0));

            if (outn_bit_buffer) {
                const uint16x8_t temp_pixel0 = vshlq_n_u16(vandq_u16(in_pixel0, xmm_3), 6);
                const uint16x8_t temp_pixel1 = vshlq_n_u16(vandq_u16(in_pixel1, xmm_3), 6);
                temp_pixel0_u8               = vcombine_u8(vqmovn_u16(temp_pixel0), vqmovn_u16(temp_pixel0));
                temp_pixel1_u8               = vcombine_u8(vqmovn_u16(temp_pixel1), vqmovn_u16(temp_pixel1));
            }

            const uint16x8_t in_pixel0_shft_r_2 = vandq_u16(vshrq_n_u16(in_pixel0, 2), xmm_00ff);
            const uint16x8_t in_pixel1_shft_r_2 = vandq_u16(vshrq_n_u16(in_pixel1, 2), xmm_00ff);

            const uint8x16_t in_pixel0_shft_r_2_u8 = vcombine_u8(vqmovn_u16(in_pixel0_shft_r_2),
                                                                 vqmovn_u16(in_pixel0_shft_r_2));
            const uint8x16_t in_pixel1_shft_r_2_u8 = vcombine_u8(vqmovn_u16(in_pixel1_shft_r_2),
                                                                 vqmovn_u16(in_pixel1_shft_r_2));

            if (outn_bit_buffer) {
                *(uint32_t *)outn_bit_buffer                 = vgetq_lane_u32(vreinterpretq_u32_u8(temp_pixel0_u8), 0);
                *(uint32_t *)(outn_bit_buffer + outn_stride) = vgetq_lane_u32(vreinterpretq_u32_u8(temp_pixel1_u8), 0);
            }
            *(uint32_t *)out8_bit_buffer = vgetq_lane_u32(vreinterpretq_u32_u8(in_pixel0_shft_r_2_u8), 0);
            *(uint32_t *)(out8_bit_buffer + out8_stride) = vgetq_lane_u32(vreinterpretq_u32_u8(in_pixel1_shft_r_2_u8),
                                                                          0);

            if (outn_bit_buffer) {
                outn_bit_buffer += 2 * outn_stride;
            }

            out8_bit_buffer += 2 * out8_stride;
            in16_bit_buffer += 2 * in_stride;
        }
    } else if (width == 8) {
        for (y = 0; y < height; y += 2) {
            const uint16x8_t in_pixel0 = vld1q_u16(in16_bit_buffer);
            const uint16x8_t in_pixel1 = vld1q_u16(in16_bit_buffer + in_stride);

            if (outn_bit_buffer) {
                const uint16x8_t temp_pixel0 = vshlq_n_u16(vandq_u16(in_pixel0, xmm_3), 6);
                const uint16x8_t temp_pixel1 = vshlq_n_u16(vandq_u16(in_pixel1, xmm_3), 6);
                temp_pixel0_u8               = vcombine_u8(vqmovn_u16(temp_pixel0), vqmovn_u16(temp_pixel0));
                temp_pixel1_u8               = vcombine_u8(vqmovn_u16(temp_pixel1), vqmovn_u16(temp_pixel1));
            }

            const uint16x8_t in_pixel0_shft_r_2 = vandq_u16(vshrq_n_u16(in_pixel0, 2), xmm_00ff);
            const uint16x8_t in_pixel1_shft_r_2 = vandq_u16(vshrq_n_u16(in_pixel1, 2), xmm_00ff);

            const uint8x16_t in_pixel0_shft_r_2_u8 = vcombine_u8(vqmovn_u16(in_pixel0_shft_r_2),
                                                                 vqmovn_u16(in_pixel0_shft_r_2));
            const uint8x16_t in_pixel1_shft_r_2_u8 = vcombine_u8(vqmovn_u16(in_pixel1_shft_r_2),
                                                                 vqmovn_u16(in_pixel1_shft_r_2));

            if (outn_bit_buffer) {
                vst1_u8(outn_bit_buffer, vget_low_u8(temp_pixel0_u8));
                vst1_u8(outn_bit_buffer + outn_stride, vget_low_u8(temp_pixel1_u8));
            }
            vst1_u8(out8_bit_buffer, vget_low_u8(in_pixel0_shft_r_2_u8));
            vst1_u8(out8_bit_buffer + out8_stride, vget_low_u8(in_pixel1_shft_r_2_u8));

            if (outn_bit_buffer) {
                outn_bit_buffer += 2 * outn_stride;
            }

            out8_bit_buffer += 2 * out8_stride;
            in16_bit_buffer += 2 * in_stride;
        }
    } else if (width == 16) {
        for (y = 0; y < height; y += 2) {
            const uint16x8_t in_pixel0 = vld1q_u16(in16_bit_buffer);
            const uint16x8_t in_pixel1 = vld1q_u16(in16_bit_buffer + 8);
            const uint16x8_t in_pixel2 = vld1q_u16(in16_bit_buffer + in_stride);
            const uint16x8_t in_pixel3 = vld1q_u16(in16_bit_buffer + in_stride + 8);

            if (outn_bit_buffer) {
                temp_pixel0_u8 = vcombine_u8(vqmovn_u16(vshlq_n_u16(vandq_u16(in_pixel0, xmm_3), 6)),
                                             vqmovn_u16(vshlq_n_u16(vandq_u16(in_pixel1, xmm_3), 6)));
                temp_pixel1_u8 = vcombine_u8(vqmovn_u16(vshlq_n_u16(vandq_u16(in_pixel2, xmm_3), 6)),
                                             vqmovn_u16(vshlq_n_u16(vandq_u16(in_pixel3, xmm_3), 6)));
            }

            const uint8x16_t in_pixel0_shft_r_2_u8 = vcombine_u8(
                vqmovn_u16(vandq_u16(vshrq_n_u16(in_pixel0, 2), xmm_00ff)),
                vqmovn_u16(vandq_u16(vshrq_n_u16(in_pixel1, 2), xmm_00ff)));
            const uint8x16_t in_pixel1_shft_r_2_u8 = vcombine_u8(
                vqmovn_u16(vandq_u16(vshrq_n_u16(in_pixel2, 2), xmm_00ff)),
                vqmovn_u16(vandq_u16(vshrq_n_u16(in_pixel3, 2), xmm_00ff)));

            if (outn_bit_buffer) {
                vst1q_u8(outn_bit_buffer, temp_pixel0_u8);
                vst1q_u8(outn_bit_buffer + outn_stride, temp_pixel1_u8);
            }
            vst1q_u8(out8_bit_buffer, in_pixel0_shft_r_2_u8);
            vst1q_u8(out8_bit_buffer + out8_stride, in_pixel1_shft_r_2_u8);

            if (outn_bit_buffer) {
                outn_bit_buffer += 2 * outn_stride;
            }

            out8_bit_buffer += 2 * out8_stride;
            in16_bit_buffer += 2 * in_stride;
        }
    } else if (width == 32) {
        uint8x16_t outn0_u8, outn1_u8, outn2_u8, outn3_u8, out8_0_u8, out8_1_u8, out8_2_u8, out8_3_u8;

        for (y = 0; y < height; y += 2) {
            const uint16x8_t in_pixel0 = vld1q_u16(in16_bit_buffer);
            const uint16x8_t in_pixel1 = vld1q_u16(in16_bit_buffer + 8);
            const uint16x8_t in_pixel2 = vld1q_u16(in16_bit_buffer + 16);
            const uint16x8_t in_pixel3 = vld1q_u16(in16_bit_buffer + 24);
            const uint16x8_t in_pixel4 = vld1q_u16(in16_bit_buffer + in_stride);
            const uint16x8_t in_pixel5 = vld1q_u16(in16_bit_buffer + in_stride + 8);
            const uint16x8_t in_pixel6 = vld1q_u16(in16_bit_buffer + in_stride + 16);
            const uint16x8_t in_pixel7 = vld1q_u16(in16_bit_buffer + in_stride + 24);

            if (outn_bit_buffer) {
                outn0_u8 = vcombine_u8(vqmovn_u16(vshlq_n_u16(vandq_u16(in_pixel0, xmm_3), 6)),
                                       vqmovn_u16(vshlq_n_u16(vandq_u16(in_pixel1, xmm_3), 6)));
                outn1_u8 = vcombine_u8(vqmovn_u16(vshlq_n_u16(vandq_u16(in_pixel2, xmm_3), 6)),
                                       vqmovn_u16(vshlq_n_u16(vandq_u16(in_pixel3, xmm_3), 6)));
                outn2_u8 = vcombine_u8(vqmovn_u16(vshlq_n_u16(vandq_u16(in_pixel4, xmm_3), 6)),
                                       vqmovn_u16(vshlq_n_u16(vandq_u16(in_pixel5, xmm_3), 6)));
                outn3_u8 = vcombine_u8(vqmovn_u16(vshlq_n_u16(vandq_u16(in_pixel6, xmm_3), 6)),
                                       vqmovn_u16(vshlq_n_u16(vandq_u16(in_pixel7, xmm_3), 6)));
            }

            out8_0_u8 = vcombine_u8(vqmovn_u16(vandq_u16(vshrq_n_u16(in_pixel0, 2), xmm_00ff)),
                                    vqmovn_u16(vandq_u16(vshrq_n_u16(in_pixel1, 2), xmm_00ff)));
            out8_1_u8 = vcombine_u8(vqmovn_u16(vandq_u16(vshrq_n_u16(in_pixel2, 2), xmm_00ff)),
                                    vqmovn_u16(vandq_u16(vshrq_n_u16(in_pixel3, 2), xmm_00ff)));
            out8_2_u8 = vcombine_u8(vqmovn_u16(vandq_u16(vshrq_n_u16(in_pixel4, 2), xmm_00ff)),
                                    vqmovn_u16(vandq_u16(vshrq_n_u16(in_pixel5, 2), xmm_00ff)));
            out8_3_u8 = vcombine_u8(vqmovn_u16(vandq_u16(vshrq_n_u16(in_pixel6, 2), xmm_00ff)),
                                    vqmovn_u16(vandq_u16(vshrq_n_u16(in_pixel7, 2), xmm_00ff)));

            if (outn_bit_buffer) {
                vst1q_u8(outn_bit_buffer, outn0_u8);
                vst1q_u8(outn_bit_buffer + 16, outn1_u8);
                vst1q_u8(outn_bit_buffer + outn_stride, outn2_u8);
                vst1q_u8(outn_bit_buffer + outn_stride + 16, outn3_u8);
            }

            vst1q_u8(out8_bit_buffer, out8_0_u8);
            vst1q_u8(out8_bit_buffer + 16, out8_1_u8);
            vst1q_u8(out8_bit_buffer + out8_stride, out8_2_u8);
            vst1q_u8(out8_bit_buffer + out8_stride + 16, out8_3_u8);

            if (outn_bit_buffer) {
                outn_bit_buffer += 2 * outn_stride;
            }

            out8_bit_buffer += 2 * out8_stride;
            in16_bit_buffer += 2 * in_stride;
        }
    } else if (width == 64) {
        uint8x16_t outn0_u8, outn1_u8, outn2_u8, outn3_u8, out8_0_u8, out8_1_u8, out8_2_u8, out8_3_u8;

        for (y = 0; y < height; ++y) {
            const uint16x8_t in_pixel0 = vld1q_u16(in16_bit_buffer);
            const uint16x8_t in_pixel1 = vld1q_u16(in16_bit_buffer + 8);
            const uint16x8_t in_pixel2 = vld1q_u16(in16_bit_buffer + 16);
            const uint16x8_t in_pixel3 = vld1q_u16(in16_bit_buffer + 24);
            const uint16x8_t in_pixel4 = vld1q_u16(in16_bit_buffer + 32);
            const uint16x8_t in_pixel5 = vld1q_u16(in16_bit_buffer + 40);
            const uint16x8_t in_pixel6 = vld1q_u16(in16_bit_buffer + 48);
            const uint16x8_t in_pixel7 = vld1q_u16(in16_bit_buffer + 56);

            if (outn_bit_buffer) {
                outn0_u8 = vcombine_u8(vqmovn_u16(vshlq_n_u16(vandq_u16(in_pixel0, xmm_3), 6)),
                                       vqmovn_u16(vshlq_n_u16(vandq_u16(in_pixel1, xmm_3), 6)));
                outn1_u8 = vcombine_u8(vqmovn_u16(vshlq_n_u16(vandq_u16(in_pixel2, xmm_3), 6)),
                                       vqmovn_u16(vshlq_n_u16(vandq_u16(in_pixel3, xmm_3), 6)));
                outn2_u8 = vcombine_u8(vqmovn_u16(vshlq_n_u16(vandq_u16(in_pixel4, xmm_3), 6)),
                                       vqmovn_u16(vshlq_n_u16(vandq_u16(in_pixel5, xmm_3), 6)));
                outn3_u8 = vcombine_u8(vqmovn_u16(vshlq_n_u16(vandq_u16(in_pixel6, xmm_3), 6)),
                                       vqmovn_u16(vshlq_n_u16(vandq_u16(in_pixel7, xmm_3), 6)));
            }

            out8_0_u8 = vcombine_u8(vqmovn_u16(vandq_u16(vshrq_n_u16(in_pixel0, 2), xmm_00ff)),
                                    vqmovn_u16(vandq_u16(vshrq_n_u16(in_pixel1, 2), xmm_00ff)));
            out8_1_u8 = vcombine_u8(vqmovn_u16(vandq_u16(vshrq_n_u16(in_pixel2, 2), xmm_00ff)),
                                    vqmovn_u16(vandq_u16(vshrq_n_u16(in_pixel3, 2), xmm_00ff)));
            out8_2_u8 = vcombine_u8(vqmovn_u16(vandq_u16(vshrq_n_u16(in_pixel4, 2), xmm_00ff)),
                                    vqmovn_u16(vandq_u16(vshrq_n_u16(in_pixel5, 2), xmm_00ff)));
            out8_3_u8 = vcombine_u8(vqmovn_u16(vandq_u16(vshrq_n_u16(in_pixel6, 2), xmm_00ff)),
                                    vqmovn_u16(vandq_u16(vshrq_n_u16(in_pixel7, 2), xmm_00ff)));

            if (outn_bit_buffer) {
                vst1q_u8(outn_bit_buffer, outn0_u8);
                vst1q_u8(outn_bit_buffer + 16, outn1_u8);
                vst1q_u8(outn_bit_buffer + 32, outn2_u8);
                vst1q_u8(outn_bit_buffer + 48, outn3_u8);
            }

            vst1q_u8(out8_bit_buffer, out8_0_u8);
            vst1q_u8(out8_bit_buffer + 16, out8_1_u8);
            vst1q_u8(out8_bit_buffer + 32, out8_2_u8);
            vst1q_u8(out8_bit_buffer + 48, out8_3_u8);

            if (outn_bit_buffer) {
                outn_bit_buffer += outn_stride;
            }

            out8_bit_buffer += out8_stride;
            in16_bit_buffer += in_stride;
        }
    } else {
        const uint32_t in_stride_diff    = (2 * in_stride) - width;
        const uint32_t out8_stride_diff  = (2 * out8_stride) - width;
        const uint32_t out_n_stride_diff = (2 * outn_stride) - width;

        const uint32_t in_stride_diff64    = in_stride - width;
        const uint32_t out8_stride_diff64  = out8_stride - width;
        const uint32_t out_n_stride_diff64 = outn_stride - width;

        if (!(width & 63)) {
            uint8x16_t outn0_u8, outn1_u8, outn2_u8, outn3_u8, out8_0_u8, out8_1_u8, out8_2_u8, out8_3_u8;

            for (x = 0; x < height; x += 1) {
                for (y = 0; y < width; y += 64) {
                    const uint16x8_t in_pixel0 = vld1q_u16(in16_bit_buffer);
                    const uint16x8_t in_pixel1 = vld1q_u16(in16_bit_buffer + 8);
                    const uint16x8_t in_pixel2 = vld1q_u16(in16_bit_buffer + 16);
                    const uint16x8_t in_pixel3 = vld1q_u16(in16_bit_buffer + 24);
                    const uint16x8_t in_pixel4 = vld1q_u16(in16_bit_buffer + 32);
                    const uint16x8_t in_pixel5 = vld1q_u16(in16_bit_buffer + 40);
                    const uint16x8_t in_pixel6 = vld1q_u16(in16_bit_buffer + 48);
                    const uint16x8_t in_pixel7 = vld1q_u16(in16_bit_buffer + 56);

                    if (outn_bit_buffer) {
                        outn0_u8 = vcombine_u8(vqmovn_u16(vshlq_n_u16(vandq_u16(in_pixel0, xmm_3), 6)),
                                               vqmovn_u16(vshlq_n_u16(vandq_u16(in_pixel1, xmm_3), 6)));
                        outn1_u8 = vcombine_u8(vqmovn_u16(vshlq_n_u16(vandq_u16(in_pixel2, xmm_3), 6)),
                                               vqmovn_u16(vshlq_n_u16(vandq_u16(in_pixel3, xmm_3), 6)));
                        outn2_u8 = vcombine_u8(vqmovn_u16(vshlq_n_u16(vandq_u16(in_pixel4, xmm_3), 6)),
                                               vqmovn_u16(vshlq_n_u16(vandq_u16(in_pixel5, xmm_3), 6)));
                        outn3_u8 = vcombine_u8(vqmovn_u16(vshlq_n_u16(vandq_u16(in_pixel6, xmm_3), 6)),
                                               vqmovn_u16(vshlq_n_u16(vandq_u16(in_pixel7, xmm_3), 6)));
                    }

                    out8_0_u8 = vcombine_u8(vqmovn_u16(vandq_u16(vshrq_n_u16(in_pixel0, 2), xmm_00ff)),
                                            vqmovn_u16(vandq_u16(vshrq_n_u16(in_pixel1, 2), xmm_00ff)));
                    out8_1_u8 = vcombine_u8(vqmovn_u16(vandq_u16(vshrq_n_u16(in_pixel2, 2), xmm_00ff)),
                                            vqmovn_u16(vandq_u16(vshrq_n_u16(in_pixel3, 2), xmm_00ff)));
                    out8_2_u8 = vcombine_u8(vqmovn_u16(vandq_u16(vshrq_n_u16(in_pixel4, 2), xmm_00ff)),
                                            vqmovn_u16(vandq_u16(vshrq_n_u16(in_pixel5, 2), xmm_00ff)));
                    out8_3_u8 = vcombine_u8(vqmovn_u16(vandq_u16(vshrq_n_u16(in_pixel6, 2), xmm_00ff)),
                                            vqmovn_u16(vandq_u16(vshrq_n_u16(in_pixel7, 2), xmm_00ff)));

                    if (outn_bit_buffer) {
                        vst1q_u8(outn_bit_buffer, outn0_u8);
                        vst1q_u8(outn_bit_buffer + 16, outn1_u8);
                        vst1q_u8(outn_bit_buffer + 32, outn2_u8);
                        vst1q_u8(outn_bit_buffer + 48, outn3_u8);
                    }

                    vst1q_u8(out8_bit_buffer, out8_0_u8);
                    vst1q_u8(out8_bit_buffer + 16, out8_1_u8);
                    vst1q_u8(out8_bit_buffer + 32, out8_2_u8);
                    vst1q_u8(out8_bit_buffer + 48, out8_3_u8);

                    if (outn_bit_buffer) {
                        outn_bit_buffer += 64;
                    }
                    out8_bit_buffer += 64;
                    in16_bit_buffer += 64;
                }
                in16_bit_buffer += in_stride_diff64;
                if (outn_bit_buffer) {
                    outn_bit_buffer += out_n_stride_diff64;
                }
                out8_bit_buffer += out8_stride_diff64;
            }
        } else if (!(width & 31)) {
            uint8x16_t outn0_u8, outn1_u8, outn2_u8, outn3_u8, out8_0_u8, out8_1_u8, out8_2_u8, out8_3_u8;

            for (x = 0; x < height; x += 2) {
                for (y = 0; y < width; y += 32) {
                    const uint16x8_t in_pixel0 = vld1q_u16(in16_bit_buffer);
                    const uint16x8_t in_pixel1 = vld1q_u16(in16_bit_buffer + 8);
                    const uint16x8_t in_pixel2 = vld1q_u16(in16_bit_buffer + 16);
                    const uint16x8_t in_pixel3 = vld1q_u16(in16_bit_buffer + 24);
                    const uint16x8_t in_pixel4 = vld1q_u16(in16_bit_buffer + in_stride);
                    const uint16x8_t in_pixel5 = vld1q_u16(in16_bit_buffer + in_stride + 8);
                    const uint16x8_t in_pixel6 = vld1q_u16(in16_bit_buffer + in_stride + 16);
                    const uint16x8_t in_pixel7 = vld1q_u16(in16_bit_buffer + in_stride + 24);

                    if (outn_bit_buffer) {
                        outn0_u8 = vcombine_u8(vqmovn_u16(vshlq_n_u16(vandq_u16(in_pixel0, xmm_3), 6)),
                                               vqmovn_u16(vshlq_n_u16(vandq_u16(in_pixel1, xmm_3), 6)));
                        outn1_u8 = vcombine_u8(vqmovn_u16(vshlq_n_u16(vandq_u16(in_pixel2, xmm_3), 6)),
                                               vqmovn_u16(vshlq_n_u16(vandq_u16(in_pixel3, xmm_3), 6)));
                        outn2_u8 = vcombine_u8(vqmovn_u16(vshlq_n_u16(vandq_u16(in_pixel4, xmm_3), 6)),
                                               vqmovn_u16(vshlq_n_u16(vandq_u16(in_pixel5, xmm_3), 6)));
                        outn3_u8 = vcombine_u8(vqmovn_u16(vshlq_n_u16(vandq_u16(in_pixel6, xmm_3), 6)),
                                               vqmovn_u16(vshlq_n_u16(vandq_u16(in_pixel7, xmm_3), 6)));
                    }

                    out8_0_u8 = vcombine_u8(vqmovn_u16(vandq_u16(vshrq_n_u16(in_pixel0, 2), xmm_00ff)),
                                            vqmovn_u16(vandq_u16(vshrq_n_u16(in_pixel1, 2), xmm_00ff)));
                    out8_1_u8 = vcombine_u8(vqmovn_u16(vandq_u16(vshrq_n_u16(in_pixel2, 2), xmm_00ff)),
                                            vqmovn_u16(vandq_u16(vshrq_n_u16(in_pixel3, 2), xmm_00ff)));
                    out8_2_u8 = vcombine_u8(vqmovn_u16(vandq_u16(vshrq_n_u16(in_pixel4, 2), xmm_00ff)),
                                            vqmovn_u16(vandq_u16(vshrq_n_u16(in_pixel5, 2), xmm_00ff)));
                    out8_3_u8 = vcombine_u8(vqmovn_u16(vandq_u16(vshrq_n_u16(in_pixel6, 2), xmm_00ff)),
                                            vqmovn_u16(vandq_u16(vshrq_n_u16(in_pixel7, 2), xmm_00ff)));

                    if (outn_bit_buffer) {
                        vst1q_u8(outn_bit_buffer, outn0_u8);
                        vst1q_u8(outn_bit_buffer + 16, outn1_u8);
                        vst1q_u8(outn_bit_buffer + outn_stride, outn2_u8);
                        vst1q_u8(outn_bit_buffer + outn_stride + 16, outn3_u8);
                    }

                    vst1q_u8(out8_bit_buffer, out8_0_u8);
                    vst1q_u8(out8_bit_buffer + 16, out8_1_u8);
                    vst1q_u8(out8_bit_buffer + out8_stride, out8_2_u8);
                    vst1q_u8(out8_bit_buffer + out8_stride + 16, out8_3_u8);

                    if (outn_bit_buffer) {
                        outn_bit_buffer += 32;
                    }
                    out8_bit_buffer += 32;
                    in16_bit_buffer += 32;
                }
                in16_bit_buffer += in_stride_diff;
                if (outn_bit_buffer) {
                    outn_bit_buffer += out_n_stride_diff;
                }
                out8_bit_buffer += out8_stride_diff;
            }
        } else if (!(width & 15)) {
            for (x = 0; x < height; x += 2) {
                for (y = 0; y < width; y += 16) {
                    const uint16x8_t in_pixel0 = vld1q_u16(in16_bit_buffer);
                    const uint16x8_t in_pixel1 = vld1q_u16(in16_bit_buffer + 8);
                    const uint16x8_t in_pixel2 = vld1q_u16(in16_bit_buffer + in_stride);
                    const uint16x8_t in_pixel3 = vld1q_u16(in16_bit_buffer + in_stride + 8);

                    if (outn_bit_buffer) {
                        temp_pixel0_u8 = vcombine_u8(vqmovn_u16(vshlq_n_u16(vandq_u16(in_pixel0, xmm_3), 6)),
                                                     vqmovn_u16(vshlq_n_u16(vandq_u16(in_pixel1, xmm_3), 6)));
                        temp_pixel1_u8 = vcombine_u8(vqmovn_u16(vshlq_n_u16(vandq_u16(in_pixel2, xmm_3), 6)),
                                                     vqmovn_u16(vshlq_n_u16(vandq_u16(in_pixel3, xmm_3), 6)));
                    }

                    const uint8x16_t in_pixel0_shft_r_2_u8 = vcombine_u8(
                        vqmovn_u16(vandq_u16(vshrq_n_u16(in_pixel0, 2), xmm_00ff)),
                        vqmovn_u16(vandq_u16(vshrq_n_u16(in_pixel1, 2), xmm_00ff)));
                    const uint8x16_t in_pixel1_shft_r_2_u8 = vcombine_u8(
                        vqmovn_u16(vandq_u16(vshrq_n_u16(in_pixel2, 2), xmm_00ff)),
                        vqmovn_u16(vandq_u16(vshrq_n_u16(in_pixel3, 2), xmm_00ff)));

                    if (outn_bit_buffer) {
                        vst1q_u8(outn_bit_buffer, temp_pixel0_u8);
                        vst1q_u8(outn_bit_buffer + outn_stride, temp_pixel1_u8);
                    }
                    vst1q_u8(out8_bit_buffer, in_pixel0_shft_r_2_u8);
                    vst1q_u8(out8_bit_buffer + out8_stride, in_pixel1_shft_r_2_u8);

                    if (outn_bit_buffer) {
                        outn_bit_buffer += 16;
                    }
                    out8_bit_buffer += 16;
                    in16_bit_buffer += 16;
                }
                in16_bit_buffer += in_stride_diff;
                if (outn_bit_buffer) {
                    outn_bit_buffer += out_n_stride_diff;
                }
                out8_bit_buffer += out8_stride_diff;
            }
        } else if (!(width & 7)) {
            for (x = 0; x < height; x += 2) {
                for (y = 0; y < width; y += 8) {
                    const uint16x8_t in_pixel0 = vld1q_u16(in16_bit_buffer);
                    const uint16x8_t in_pixel1 = vld1q_u16(in16_bit_buffer + in_stride);

                    if (outn_bit_buffer) {
                        const uint16x8_t temp_pixel0 = vshlq_n_u16(vandq_u16(in_pixel0, xmm_3), 6);
                        const uint16x8_t temp_pixel1 = vshlq_n_u16(vandq_u16(in_pixel1, xmm_3), 6);
                        temp_pixel0_u8               = vcombine_u8(vqmovn_u16(temp_pixel0), vqmovn_u16(temp_pixel0));
                        temp_pixel1_u8               = vcombine_u8(vqmovn_u16(temp_pixel1), vqmovn_u16(temp_pixel1));
                    }

                    const uint16x8_t in_pixel0_shft_r_2 = vandq_u16(vshrq_n_u16(in_pixel0, 2), xmm_00ff);
                    const uint16x8_t in_pixel1_shft_r_2 = vandq_u16(vshrq_n_u16(in_pixel1, 2), xmm_00ff);

                    const uint8x16_t in_pixel0_shft_r_2_u8 = vcombine_u8(vqmovn_u16(in_pixel0_shft_r_2),
                                                                         vqmovn_u16(in_pixel0_shft_r_2));
                    const uint8x16_t in_pixel1_shft_r_2_u8 = vcombine_u8(vqmovn_u16(in_pixel1_shft_r_2),
                                                                         vqmovn_u16(in_pixel1_shft_r_2));

                    if (outn_bit_buffer) {
                        vst1_u8(outn_bit_buffer, vget_low_u8(temp_pixel0_u8));
                        vst1_u8(outn_bit_buffer + outn_stride, vget_low_u8(temp_pixel1_u8));
                    }
                    vst1_u8(out8_bit_buffer, vget_low_u8(in_pixel0_shft_r_2_u8));
                    vst1_u8(out8_bit_buffer + out8_stride, vget_low_u8(in_pixel1_shft_r_2_u8));

                    if (outn_bit_buffer) {
                        outn_bit_buffer += 8;
                    }
                    out8_bit_buffer += 8;
                    in16_bit_buffer += 8;
                }
                in16_bit_buffer += in_stride_diff;
                if (outn_bit_buffer) {
                    outn_bit_buffer += out_n_stride_diff;
                }
                out8_bit_buffer += out8_stride_diff;
            }
        } else {
            for (x = 0; x < height; x += 2) {
                for (y = 0; y < width; y += 4) {
                    const uint16x8_t in_pixel0 = vcombine_u16(vld1_u16(in16_bit_buffer), vdup_n_u16(0));
                    const uint16x8_t in_pixel1 = vcombine_u16(vld1_u16(in16_bit_buffer + in_stride), vdup_n_u16(0));

                    if (outn_bit_buffer) {
                        const uint16x8_t temp_pixel0 = vshlq_n_u16(vandq_u16(in_pixel0, xmm_3), 6);
                        const uint16x8_t temp_pixel1 = vshlq_n_u16(vandq_u16(in_pixel1, xmm_3), 6);
                        temp_pixel0_u8               = vcombine_u8(vqmovn_u16(temp_pixel0), vqmovn_u16(temp_pixel0));
                        temp_pixel1_u8               = vcombine_u8(vqmovn_u16(temp_pixel1), vqmovn_u16(temp_pixel1));
                    }

                    const uint16x8_t in_pixel0_shft_r_2 = vandq_u16(vshrq_n_u16(in_pixel0, 2), xmm_00ff);
                    const uint16x8_t in_pixel1_shft_r_2 = vandq_u16(vshrq_n_u16(in_pixel1, 2), xmm_00ff);

                    const uint8x16_t in_pixel0_shft_r_2_u8 = vcombine_u8(vqmovn_u16(in_pixel0_shft_r_2),
                                                                         vqmovn_u16(in_pixel0_shft_r_2));
                    const uint8x16_t in_pixel1_shft_r_2_u8 = vcombine_u8(vqmovn_u16(in_pixel1_shft_r_2),
                                                                         vqmovn_u16(in_pixel1_shft_r_2));

                    if (outn_bit_buffer) {
                        *(uint32_t *)outn_bit_buffer = vgetq_lane_u32(vreinterpretq_u32_u8(temp_pixel0_u8), 0);
                        *(uint32_t *)(outn_bit_buffer + outn_stride) = vgetq_lane_u32(
                            vreinterpretq_u32_u8(temp_pixel1_u8), 0);
                    }
                    *(uint32_t *)out8_bit_buffer = vgetq_lane_u32(vreinterpretq_u32_u8(in_pixel0_shft_r_2_u8), 0);
                    *(uint32_t *)(out8_bit_buffer + out8_stride) = vgetq_lane_u32(
                        vreinterpretq_u32_u8(in_pixel1_shft_r_2_u8), 0);
                    if (outn_bit_buffer) {
                        outn_bit_buffer += 4;
                    }
                    out8_bit_buffer += 4;
                    in16_bit_buffer += 4;
                }
                in16_bit_buffer += in_stride_diff;
                if (outn_bit_buffer) {
                    outn_bit_buffer += out_n_stride_diff;
                }
                out8_bit_buffer += out8_stride_diff;
            }
        }
    }
}
