// clang-format off
/*
 * Copyright(c) 2019 Netflix, Inc.
 * Copyright(c) 2019 Intel Corporation
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 3-Clause Clear License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 3-Clause Clear License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at https://www.aomedia.org/license. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "EbTemporalFiltering.h"
#include "EbComputeSAD.h"
#include "EbMotionEstimation.h"
#include "EbMotionEstimationProcess.h"
#include "EbMotionEstimationContext.h"
#include "EbLambdaRateTables.h"
#include "EbPictureAnalysisProcess.h"
#include "EbModeDecisionProcess.h"
#include "EbMcp.h"
#include "av1me.h"
#ifdef ARCH_X86_64
#include <xmmintrin.h>
#endif
#include "EbObject.h"
#include "EbEncInterPrediction.h"
#include "EbLog.h"
#include <limits.h>
#include "EbPackUnPack_C.h"

#undef _MM_HINT_T2
#define _MM_HINT_T2 1

#include "EbPictureDecisionResults.h"
#include "EbUtility.h"

static const uint32_t subblock_xy_16x16[N_16X16_BLOCKS][2] = {{0, 0},
                                                              {0, 1},
                                                              {0, 2},
                                                              {0, 3},
                                                              {1, 0},
                                                              {1, 1},
                                                              {1, 2},
                                                              {1, 3},
                                                              {2, 0},
                                                              {2, 1},
                                                              {2, 2},
                                                              {2, 3},
                                                              {3, 0},
                                                              {3, 1},
                                                              {3, 2},
                                                              {3, 3}};
static const uint32_t idx_32x32_to_idx_16x16[4][4]         = {
    {0, 1, 4, 5}, {2, 3, 6, 7}, {8, 9, 12, 13}, {10, 11, 14, 15}};

int32_t get_frame_update_type(SequenceControlSet *scs_ptr, PictureParentControlSet *pcs_ptr);
int32_t svt_av1_compute_qdelta_fp(int32_t qstart_fp8, int32_t qtarget_fp8, EbBitDepth bit_depth);
int32_t svt_av1_compute_qdelta(double qstart, double qtarget, EbBitDepth bit_depth);
void generate_padding_compressed_10bit(EbByte src_pic, uint32_t src_stride,
                                       uint32_t original_src_width, uint32_t original_src_height,
                                       uint32_t padding_width, uint32_t padding_height);
void svt_c_unpack_compressed_10bit(const uint8_t *inn_bit_buffer, uint32_t inn_stride,
                                   uint8_t *in_compn_bit_buffer, uint32_t out_stride,
                                   uint32_t height);
#if DEBUG_SCALING
// save YUV to file - auxiliary function for debug
void save_YUV_to_file(char *filename, EbByte buffer_y, EbByte buffer_u, EbByte buffer_v,
                      uint16_t width, uint16_t height, uint16_t stride_y, uint16_t stride_u,
                      uint16_t stride_v, uint16_t origin_y, uint16_t origin_x, uint32_t ss_x,
                      uint32_t ss_y) {
    FILE *fid;

    // save current source picture to a YUV file
    FOPEN(fid, filename, "wb");

    if (!fid) {
        SVT_LOG("Unable to open file %s to write.\n", "temp_picture.yuv");
    } else {
        // the source picture saved in the enchanced_picture_ptr contains a border in x and y dimensions
        EbByte pic_point = buffer_y + (origin_y * stride_y) + origin_x;
        for (int h = 0; h < height; h++) {
            fwrite(pic_point, 1, (size_t)width, fid);
            pic_point = pic_point + stride_y;
        }
        pic_point = buffer_u + ((origin_y >> ss_y) * stride_u) + (origin_x >> ss_x);
        for (int h = 0; h < (height >> ss_y); h++) {
            fwrite(pic_point, 1, (size_t)width >> ss_x, fid);
            pic_point = pic_point + stride_u;
        }
        pic_point = buffer_v + ((origin_y >> ss_y) * stride_v) + (origin_x >> ss_x);
        for (int h = 0; h < (height >> ss_y); h++) {
            fwrite(pic_point, 1, (size_t)width >> ss_x, fid);
            pic_point = pic_point + stride_v;
        }
        fclose(fid);
    }
}

// save YUV to file - auxiliary function for debug
void save_YUV_to_file_highbd(char *filename, uint16_t *buffer_y, uint16_t *buffer_u,
                             uint16_t *buffer_v, uint16_t width, uint16_t height, uint16_t stride_y,
                             uint16_t stride_u, uint16_t stride_v, uint16_t origin_y,
                             uint16_t origin_x, uint32_t ss_x, uint32_t ss_y) {
    FILE *fid;

    // save current source picture to a YUV file
    FOPEN(fid, filename, "wb");

    if (!fid) {
        SVT_LOG("Unable to open file %s to write.\n", "temp_picture.yuv");
    } else {
        // the source picture saved in the enchanced_picture_ptr contains a border in x and y dimensions
        uint16_t *pic_point = buffer_y + (origin_y * stride_y) + origin_x;
        for (int h = 0; h < height; h++) {
            fwrite(pic_point, 2, (size_t)width, fid);
            pic_point = pic_point + stride_y;
        }
        pic_point = buffer_u + ((origin_y >> ss_y) * stride_u) + (origin_x >> ss_x);
        for (int h = 0; h < (height >> ss_y); h++) {
            fwrite(pic_point, 2, (size_t)width >> ss_x, fid);

            pic_point = pic_point + stride_u;
        }
        pic_point = buffer_v + ((origin_y >> ss_y) * stride_v) + (origin_x >> ss_x);
        for (int h = 0; h < (height >> ss_y); h++) {
            fwrite(pic_point, 2, (size_t)width >> ss_x, fid);
            pic_point = pic_point + stride_v;
        }
        fclose(fid);
    }
}
#endif
void pack_highbd_pic(const EbPictureBufferDesc *pic_ptr, uint16_t *buffer_16bit[3], uint32_t ss_x,
                     uint32_t ss_y, Bool include_padding) {
    uint16_t width  = pic_ptr->stride_y;
    uint16_t height = (uint16_t)(pic_ptr->origin_y * 2 + pic_ptr->height);

    assert_err(include_padding == 1, "not supporting OFF");

    uint32_t comp_stride_y = pic_ptr->stride_y / 4;

    compressed_pack_sb(pic_ptr->buffer_y,
                       pic_ptr->stride_y,
                       pic_ptr->buffer_bit_inc_y,
                       comp_stride_y,
                       buffer_16bit[C_Y],
                       pic_ptr->stride_y,
                       width,
                       height);

    uint32_t comp_stride_uv = pic_ptr->stride_cb / 4;
    if (buffer_16bit[C_U])
        compressed_pack_sb(pic_ptr->buffer_cb,
                           pic_ptr->stride_cb,
                           pic_ptr->buffer_bit_inc_cb,
                           comp_stride_uv,
                           buffer_16bit[C_U],
                           pic_ptr->stride_cb,
                           (width + ss_x) >> ss_x,
                           (height + ss_y) >> ss_y);
    if (buffer_16bit[C_V])
        compressed_pack_sb(pic_ptr->buffer_cr,
                           pic_ptr->stride_cr,
                           pic_ptr->buffer_bit_inc_cr,
                           comp_stride_uv,
                           buffer_16bit[C_V],
                           pic_ptr->stride_cr,
                           (width + ss_x) >> ss_x,
                           (height + ss_y) >> ss_y);
}

void unpack_highbd_pic(uint16_t *buffer_highbd[3], EbPictureBufferDesc *pic_ptr, uint32_t ss_x,
                       uint32_t ss_y, Bool include_padding) {
    uint16_t width  = pic_ptr->stride_y;
    uint16_t height = (uint16_t)(pic_ptr->origin_y * 2 + pic_ptr->height);

    assert_err(include_padding == 1, "not supporting OFF");

    uint32_t comp_stride_y  = pic_ptr->stride_y / 4;
    uint32_t comp_stride_uv = pic_ptr->stride_cb / 4;

    svt_unpack_and_2bcompress(
        buffer_highbd[C_Y],
        pic_ptr->stride_y,
        pic_ptr->buffer_y,
        pic_ptr->stride_y,
        pic_ptr->buffer_bit_inc_y,
        comp_stride_y,
        width,
        height);

    if (buffer_highbd[C_U])
        svt_unpack_and_2bcompress(
            buffer_highbd[C_U],
            pic_ptr->stride_cb,
            pic_ptr->buffer_cb,
            pic_ptr->stride_cb,
            pic_ptr->buffer_bit_inc_cb,
            comp_stride_uv,
            (width + ss_x) >> ss_x,
            (height + ss_y) >> ss_y);

    if (buffer_highbd[C_V])
        svt_unpack_and_2bcompress(
            buffer_highbd[C_V],
            pic_ptr->stride_cr,
            pic_ptr->buffer_cr,
            pic_ptr->stride_cr,
            pic_ptr->buffer_bit_inc_cr,
            comp_stride_uv,
            (width + ss_x) >> ss_x,
            (height + ss_y) >> ss_y);
}

void generate_padding_pic(EbPictureBufferDesc *pic_ptr, uint32_t ss_x, uint32_t ss_y,
                          Bool is_highbd) {
    if (!is_highbd) {
        generate_padding(pic_ptr->buffer_cb,
                         pic_ptr->stride_cb,
                         pic_ptr->width >> ss_x,
                         pic_ptr->height >> ss_y,
                         pic_ptr->origin_x >> ss_x,
                         pic_ptr->origin_y >> ss_y);

        generate_padding(pic_ptr->buffer_cr,
                         pic_ptr->stride_cr,
                         pic_ptr->width >> ss_x,
                         pic_ptr->height >> ss_y,
                         pic_ptr->origin_x >> ss_x,
                         pic_ptr->origin_y >> ss_y);
    } else {
        generate_padding(pic_ptr->buffer_cb,
                         pic_ptr->stride_cb,
                         pic_ptr->width >> ss_x,
                         pic_ptr->height >> ss_y,
                         pic_ptr->origin_x >> ss_x,
                         pic_ptr->origin_y >> ss_y);

        generate_padding(pic_ptr->buffer_cr,
                         pic_ptr->stride_cr,
                         pic_ptr->width >> ss_x,
                         pic_ptr->height >> ss_y,
                         pic_ptr->origin_x >> ss_x,
                         pic_ptr->origin_y >> ss_y);

        generate_padding_compressed_10bit(pic_ptr->buffer_bit_inc_cb,
                                          pic_ptr->stride_cr / 4,
                                          pic_ptr->width >> ss_x,
                                          pic_ptr->height >> ss_y,
                                          pic_ptr->origin_x >> ss_x,
                                          pic_ptr->origin_y >> ss_y);

        generate_padding_compressed_10bit(pic_ptr->buffer_bit_inc_cr,
                                          pic_ptr->stride_cr / 4,
                                          pic_ptr->width >> ss_x,
                                          pic_ptr->height >> ss_y,
                                          pic_ptr->origin_x >> ss_x,
                                          pic_ptr->origin_y >> ss_y);
    }
}
static void derive_tf_32x32_block_split_flag(MeContext *context_ptr) {
    int      subblock_errors[4];
    uint32_t idx_32x32   = context_ptr->idx_32x32;
    int      block_error = (int)context_ptr->tf_32x32_block_error[idx_32x32];

    // `block_error` is initialized as INT_MAX and will be overwritten after
    // motion search with reference frame, therefore INT_MAX can ONLY be accessed
    // by to-filter frame.
    if (block_error == INT_MAX) {
        context_ptr->tf_32x32_block_split_flag[idx_32x32] = 0;
    }

    int min_subblock_error = INT_MAX;
    int max_subblock_error = INT_MIN;
    int sum_subblock_error = 0;
    for (int i = 0; i < 4; ++i) {
        subblock_errors[i] = (int)context_ptr->tf_16x16_block_error[idx_32x32 * 4 + i];

        sum_subblock_error += subblock_errors[i];
        min_subblock_error = AOMMIN(min_subblock_error, subblock_errors[i]);
        max_subblock_error = AOMMAX(max_subblock_error, subblock_errors[i]);
    }
    if (block_error * 14 < sum_subblock_error * 16) { // No split.
        context_ptr->tf_32x32_block_split_flag[idx_32x32] = 0;
    } else { // Do split.
        context_ptr->tf_32x32_block_split_flag[idx_32x32] = 1;
    }
}
// Create and initialize all necessary ME context structures
static void create_me_context_and_picture_control(
    MotionEstimationContext_t *context_ptr, PictureParentControlSet *picture_control_set_ptr_frame,
    PictureParentControlSet *picture_control_set_ptr_central,
    EbPictureBufferDesc *input_picture_ptr_central, int blk_row, int blk_col, uint32_t ss_x,
    uint32_t ss_y) {
    // set reference picture for alt-refs
    context_ptr->me_context_ptr->alt_ref_reference_ptr =
        (EbPaReferenceObject *)
            picture_control_set_ptr_frame->pa_reference_picture_wrapper_ptr->object_ptr;
    context_ptr->me_context_ptr->me_type = ME_MCTF;

    // set the buffers with the original, quarter and sixteenth pixels version of the source frame
    EbPaReferenceObject *src_object = (EbPaReferenceObject *)picture_control_set_ptr_central
                                          ->pa_reference_picture_wrapper_ptr->object_ptr;
    EbPictureBufferDesc *padded_pic_ptr = src_object->input_padded_picture_ptr;
    // Set 1/4 and 1/16 ME reference buffer(s); filtered or decimated
    EbPictureBufferDesc *quarter_pic_ptr   = src_object->quarter_downsampled_picture_ptr;
    EbPictureBufferDesc *sixteenth_pic_ptr = src_object->sixteenth_downsampled_picture_ptr;
    // Parts from MotionEstimationKernel()
    uint32_t b64_origin_x = (uint32_t)(blk_col * BW);
    uint32_t b64_origin_y = (uint32_t)(blk_row * BH);

    // Load the SB from the input to the intermediate SB buffer
    int buffer_index = (input_picture_ptr_central->origin_y + b64_origin_y) *
        input_picture_ptr_central->stride_y +
        input_picture_ptr_central->origin_x + b64_origin_x;

    // set search method
    context_ptr->me_context_ptr->hme_search_method = FULL_SAD_SEARCH;
#ifdef ARCH_X86_64
    {
        uint8_t *src_ptr = &(padded_pic_ptr->buffer_y[buffer_index]);
        uint32_t b64_height = (input_picture_ptr_central->height - b64_origin_y) < BLOCK_SIZE_64
            ? input_picture_ptr_central->height - b64_origin_y
            : BLOCK_SIZE_64;
        //_MM_HINT_T0     //_MM_HINT_T1    //_MM_HINT_T2    //_MM_HINT_NTA
        uint32_t i;
        for (i = 0; i < b64_height; i++) {
            char const *p = (char const *)(src_ptr + i * padded_pic_ptr->stride_y);

            _mm_prefetch(p, _MM_HINT_T2);
        }
    }
#endif
    context_ptr->me_context_ptr->b64_src_ptr    = &(padded_pic_ptr->buffer_y[buffer_index]);
    context_ptr->me_context_ptr->b64_src_stride = padded_pic_ptr->stride_y;

    // Load the 1/4 decimated SB from the 1/4 decimated input to the 1/4 intermediate SB buffer
    buffer_index = (quarter_pic_ptr->origin_y + (b64_origin_y >> ss_y)) * quarter_pic_ptr->stride_y +
        quarter_pic_ptr->origin_x + (b64_origin_x >> ss_x);

    context_ptr->me_context_ptr->quarter_b64_buffer = &quarter_pic_ptr->buffer_y[buffer_index];
    context_ptr->me_context_ptr->quarter_b64_buffer_stride = quarter_pic_ptr->stride_y;

    // Load the 1/16 decimated SB from the 1/16 decimated input to the 1/16 intermediate SB buffer
    buffer_index = (sixteenth_pic_ptr->origin_y + (b64_origin_y >> 2)) * sixteenth_pic_ptr->stride_y +
        sixteenth_pic_ptr->origin_x + (b64_origin_x >> 2);

    context_ptr->me_context_ptr->sixteenth_b64_buffer = &sixteenth_pic_ptr->buffer_y[buffer_index];
    context_ptr->me_context_ptr->sixteenth_b64_buffer_stride = sixteenth_pic_ptr->stride_y;
}

static INLINE void calculate_squared_errors(const uint8_t *s, int s_stride, const uint8_t *p,
                                            int p_stride, uint16_t *diff_sse, unsigned int w,
                                            unsigned int h) {
    int          idx = 0;
    unsigned int i, j;

    for (i = 0; i < h; i++) {
        for (j = 0; j < w; j++) {
            const int16_t diff = s[i * s_stride + j] - p[i * p_stride + j];
            diff_sse[idx]      = (uint16_t)(diff * diff);
            idx++;
        }
    }
}

static INLINE void calculate_squared_errors_highbd(const uint16_t *s, int s_stride,
                                                   const uint16_t *p, int p_stride,
                                                   uint32_t *diff_sse, unsigned int w,
                                                   unsigned int h) {
    int          idx = 0;
    unsigned int i, j;

    for (i = 0; i < h; i++) {
        for (j = 0; j < w; j++) {
            const int32_t diff = s[i * s_stride + j] - p[i * p_stride + j];
            diff_sse[idx]      = (uint32_t)(diff * diff);
            idx++;
        }
    }
}

// Apply filtering to the central picture
void apply_filtering_central_c(MeContext *          context_ptr,
                               EbPictureBufferDesc *input_picture_ptr_central, EbByte *src,
                               uint32_t **accum, uint16_t **count, uint16_t blk_width,
                               uint16_t blk_height, uint32_t ss_x, uint32_t ss_y) {
    uint16_t blk_height_y  = blk_height;
    uint16_t blk_width_y   = blk_width;
    uint16_t blk_height_ch = blk_height >> ss_y;
    uint16_t blk_width_ch  = blk_width >> ss_x;
    uint16_t src_stride_y  = input_picture_ptr_central->stride_y;
    uint16_t src_stride_ch = src_stride_y >> ss_x;

    const int modifier = TF_PLANEWISE_FILTER_WEIGHT_SCALE;

    // Luma
    for (uint16_t k = 0, i = 0; i < blk_height_y; i++) {
        for (uint16_t j = 0; j < blk_width_y; j++) {
            accum[C_Y][k] = modifier * src[C_Y][i * src_stride_y + j];
            count[C_Y][k] = modifier;
            ++k;
        }
    }

    // Chroma
    if (context_ptr->tf_chroma)
        for (uint16_t k = 0, i = 0; i < blk_height_ch; i++) {
            for (uint16_t j = 0; j < blk_width_ch; j++) {
                accum[C_U][k] = modifier * src[C_U][i * src_stride_ch + j];
                count[C_U][k] = modifier;

                accum[C_V][k] = modifier * src[C_V][i * src_stride_ch + j];
                count[C_V][k] = modifier;
                ++k;
            }
        }
}

// Apply filtering to the central picture
void apply_filtering_central_highbd_c(MeContext *          context_ptr,
                                      EbPictureBufferDesc *input_picture_ptr_central,
                                      uint16_t **src_16bit, uint32_t **accum, uint16_t **count,
                                      uint16_t blk_width, uint16_t blk_height, uint32_t ss_x,
                                      uint32_t ss_y) {
    uint16_t blk_height_y  = blk_height;
    uint16_t blk_width_y   = blk_width;
    uint16_t blk_height_ch = blk_height >> ss_y;
    uint16_t blk_width_ch  = blk_width >> ss_x;
    uint16_t src_stride_y  = input_picture_ptr_central->stride_y;
    uint16_t src_stride_ch = src_stride_y >> ss_x;

    const int modifier = TF_PLANEWISE_FILTER_WEIGHT_SCALE;

    // Luma
    for (uint16_t k = 0, i = 0; i < blk_height_y; i++) {
        for (uint16_t j = 0; j < blk_width_y; j++) {
            accum[C_Y][k] = modifier * src_16bit[C_Y][i * src_stride_y + j];
            count[C_Y][k] = modifier;
            ++k;
        }
    }

    // Chroma
    if (context_ptr->tf_chroma)
        for (uint16_t k = 0, i = 0; i < blk_height_ch; i++) {
            for (uint16_t j = 0; j < blk_width_ch; j++) {
                accum[C_U][k] = modifier * src_16bit[C_U][i * src_stride_ch + j];
                count[C_U][k] = modifier;

                accum[C_V][k] = modifier * src_16bit[C_V][i * src_stride_ch + j];
                count[C_V][k] = modifier;
                ++k;
            }
        }
}

//log1p(x) for x in [-1..6], step 1/32 values in Fixed Points shift 16
static const int32_t log1p_tab_fp16[] = {
    (int32_t)-2147483647 - 1,
    (int32_t)-227130,
    (int32_t)-181704,
    (int32_t)-155131,
    (int32_t)-136278,
    (int32_t)-121654,
    (int32_t)-109705,
    (int32_t)-99603,
    (int32_t)-90852,
    (int32_t)-83133,
    (int32_t)-76228,
    (int32_t)-69982,
    (int32_t)-64279,
    (int32_t)-59033,
    (int32_t)-54177,
    (int32_t)-49655,
    (int32_t)-45426,
    (int32_t)-41452,
    (int32_t)-37707,
    (int32_t)-34163,
    (int32_t)-30802,
    (int32_t)-27604,
    (int32_t)-24555,
    (int32_t)-21642,
    (int32_t)-18853,
    (int32_t)-16178,
    (int32_t)-13607,
    (int32_t)-11134,
    (int32_t)-8751,
    (int32_t)-6451,
    (int32_t)-4229,
    (int32_t)-2080,
    (int32_t)0,
    (int32_t)2016,
    (int32_t)3973,
    (int32_t)5872,
    (int32_t)7719,
    (int32_t)9514,
    (int32_t)11262,
    (int32_t)12964,
    (int32_t)14623,
    (int32_t)16242,
    (int32_t)17821,
    (int32_t)19363,
    (int32_t)20870,
    (int32_t)22342,
    (int32_t)23783,
    (int32_t)25192,
    (int32_t)26572,
    (int32_t)27923,
    (int32_t)29247,
    (int32_t)30545,
    (int32_t)31818,
    (int32_t)33066,
    (int32_t)34291,
    (int32_t)35494,
    (int32_t)36674,
    (int32_t)37834,
    (int32_t)38974,
    (int32_t)40095,
    (int32_t)41196,
    (int32_t)42279,
    (int32_t)43345,
    (int32_t)44394,
    (int32_t)45426,
    (int32_t)46442,
    (int32_t)47442,
    (int32_t)48428,
    (int32_t)49399,
    (int32_t)50355,
    (int32_t)51298,
    (int32_t)52228,
    (int32_t)53145,
    (int32_t)54049,
    (int32_t)54940,
    (int32_t)55820,
    (int32_t)56688,
    (int32_t)57545,
    (int32_t)58390,
    (int32_t)59225,
    (int32_t)60050,
    (int32_t)60864,
    (int32_t)61668,
    (int32_t)62462,
    (int32_t)63247,
    (int32_t)64023,
    (int32_t)64789,
    (int32_t)65547,
    (int32_t)66296,
    (int32_t)67036,
    (int32_t)67769,
    (int32_t)68493,
    (int32_t)69209,
    (int32_t)69917,
    (int32_t)70618,
    (int32_t)71312,
    (int32_t)71998,
    (int32_t)72677,
    (int32_t)73349,
    (int32_t)74015,
    (int32_t)74673,
    (int32_t)75326,
    (int32_t)75971,
    (int32_t)76611,
    (int32_t)77244,
    (int32_t)77871,
    (int32_t)78492,
    (int32_t)79108,
    (int32_t)79717,
    (int32_t)80321,
    (int32_t)80920,
    (int32_t)81513,
    (int32_t)82101,
    (int32_t)82683,
    (int32_t)83261,
    (int32_t)83833,
    (int32_t)84400,
    (int32_t)84963,
    (int32_t)85521,
    (int32_t)86074,
    (int32_t)86622,
    (int32_t)87166,
    (int32_t)87705,
    (int32_t)88240,
    (int32_t)88771,
    (int32_t)89297,
    (int32_t)89820,
    (int32_t)90338,
    (int32_t)90852,
    (int32_t)91362,
    (int32_t)91868,
    (int32_t)92370,
    (int32_t)92868,
    (int32_t)93363,
    (int32_t)93854,
    (int32_t)94341,
    (int32_t)94825,
    (int32_t)95305,
    (int32_t)95782,
    (int32_t)96255,
    (int32_t)96725,
    (int32_t)97191,
    (int32_t)97654,
    (int32_t)98114,
    (int32_t)98571,
    (int32_t)99024,
    (int32_t)99475,
    (int32_t)99922,
    (int32_t)100366,
    (int32_t)100808,
    (int32_t)101246,
    (int32_t)101681,
    (int32_t)102114,
    (int32_t)102544,
    (int32_t)102971,
    (int32_t)103395,
    (int32_t)103816,
    (int32_t)104235,
    (int32_t)104651,
    (int32_t)105065,
    (int32_t)105476,
    (int32_t)105884,
    (int32_t)106290,
    (int32_t)106693,
    (int32_t)107094,
    (int32_t)107492,
    (int32_t)107888,
    (int32_t)108282,
    (int32_t)108673,
    (int32_t)109062,
    (int32_t)109449,
    (int32_t)109833,
    (int32_t)110215,
    (int32_t)110595,
    (int32_t)110973,
    (int32_t)111348,
    (int32_t)111722,
    (int32_t)112093,
    (int32_t)112462,
    (int32_t)112830,
    (int32_t)113195,
    (int32_t)113558,
    (int32_t)113919,
    (int32_t)114278,
    (int32_t)114635,
    (int32_t)114990,
    (int32_t)115344,
    (int32_t)115695,
    (int32_t)116044,
    (int32_t)116392,
    (int32_t)116738,
    (int32_t)117082,
    (int32_t)117424,
    (int32_t)117765,
    (int32_t)118103,
    (int32_t)118440,
    (int32_t)118776,
    (int32_t)119109,
    (int32_t)119441,
    (int32_t)119771,
    (int32_t)120100,
    (int32_t)120426,
    (int32_t)120752,
    (int32_t)121075,
    (int32_t)121397,
    (int32_t)121718,
    (int32_t)122037,
    (int32_t)122354,
    (int32_t)122670,
    (int32_t)122984,
    (int32_t)123297,
    (int32_t)123608,
    (int32_t)123918,
    (int32_t)124227,
    (int32_t)124534,
    (int32_t)124839,
    (int32_t)125143,
    (int32_t)125446,
    (int32_t)125747,
    (int32_t)126047,
    (int32_t)126346,
    (int32_t)126643,
    (int32_t)126939,
    (int32_t)127233,
    (int32_t)127527,
};

int32_t noise_log1p_fp16(int32_t noise_level_fp16) {
    int32_t base_fp16 = (65536 /*1:fp16*/ + noise_level_fp16);

    if (base_fp16 <= 0) {
        return (-2147483647 - 1) /*-MAX*/;
    } else if (base_fp16 < (458752) /*7:fp16*/) {
        //Aproximate value:
        int32_t id = base_fp16 >> 11; // //step 1/32 so multiple by 32 is reduce shift of 5
        FP_ASSERT(((size_t)id + 1) < sizeof(log1p_tab_fp16) / sizeof(log1p_tab_fp16[0]));
        int32_t rest     = base_fp16 & 0x7FF; //11 bits
        int32_t diff     = ((rest * (log1p_tab_fp16[id + 1] - log1p_tab_fp16[id])) >> 11);
        int32_t val_fp16 = log1p_tab_fp16[id] +
            diff; // + (rest*(log_tab_fp16[id + 1] - log_tab_fp16[id]))>>11;
        return val_fp16;
    } else {
        //approximation to line(fp16): y=1860*x + 116456
        FP_ASSERT((int64_t)(noise_level_fp16 >> 8) * 1860 < ((int64_t)1 << 31));
        return ((1860 * (noise_level_fp16 >> 8)) >> 8) + 116456;
    }
}

// T[X] =  exp(-(X)/16)  for x in [0..7], step 1/16 values in Fixed Points shift 16
static const int32_t expf_tab_fp16[] = {
    65536, 61565, 57835, 54331, 51039, 47947, 45042, 42313, 39749, 37341, 35078, 32953, 30957,
    29081, 27319, 25664, 24109, 22648, 21276, 19987, 18776, 17638, 16570, 15566, 14623, 13737,
    12904, 12122, 11388, 10698, 10050, 9441,  8869,  8331,  7827,  7352,  6907,  6488,  6095,
    5726,  5379,  5053,  4747,  4459,  4189,  3935,  3697,  3473,  3262,  3065,  2879,  2704,
    2541,  2387,  2242,  2106,  1979,  1859,  1746,  1640,  1541,  1447,  1360,  1277,  1200,
    1127,  1059,  995,   934,   878,   824,   774,   728,   683,   642,   603,   566,   532,
    500,   470,   441,   414,   389,   366,   343,   323,   303,   285,   267,   251,   236,
    222,   208,   195,   184,   172,   162,   152,   143,   134,   126,   118,   111,   104,
    98,    92,    86,    81,    76,    72,    67,    63,    59,    56,    52,    49,    46,
    43,    41,    38,    36,    34,    31,    30,    28,    26,    24,    23,    21};

/*value [i:0-15] (sqrt((float)i)*65536.0*/
static const uint32_t sqrt_array_fp16[16] = {0,
                                             65536,
                                             92681,
                                             113511,
                                             131072,
                                             146542,
                                             160529,
                                             173391,
                                             185363,
                                             196608,
                                             207243,
                                             217358,
                                             227023,
                                             236293,
                                             245213,
                                             253819};

/*Calc sqrt linear max error 10%*/
static uint32_t sqrt_fast(uint32_t x) {
    if (x > 15) {
        const int log2_half = svt_log2f(x) >> 1;
        const int mul2      = log2_half << 1;
        int       base      = x >> (mul2 - 2);
        assert(base < 16);
        return sqrt_array_fp16[base] >> (17 - log2_half);
    }
    return sqrt_array_fp16[x] >> 16;
}
//exp(-x) for x in [0..7]
double expf_tab[] = {1,        0.904837, 0.818731, 0.740818, 0.67032,  0.606531, 0.548812, 0.496585,
                     0.449329, 0.40657,  0.367879, 0.332871, 0.301194, 0.272532, 0.246597, 0.22313,
                     0.201896, 0.182683, 0.165299, 0.149569, 0.135335, 0.122456, 0.110803, 0.100259,
                     0.090718, 0.082085, 0.074274, 0.067206, 0.06081,  0.055023, 0.049787, 0.045049,
                     0.040762, 0.036883, 0.033373, 0.030197, 0.027324, 0.024724, 0.022371, 0.020242,
                     0.018316, 0.016573, 0.014996, 0.013569, 0.012277, 0.011109, 0.010052, 0.009095,
                     0.00823,  0.007447, 0.006738, 0.006097, 0.005517, 0.004992, 0.004517, 0.004087,
                     0.003698, 0.003346, 0.003028, 0.002739, 0.002479, 0.002243, 0.002029, 0.001836,
                     0.001662, 0.001503, 0.00136,  0.001231, 0.001114, 0.001008, 0.000912, 0.000825,
                     0.000747, 0.000676, 0.000611, 0.000553, 0.0005,   0.000453, 0.00041,  0.000371,
                     0.000335

};

static float exp_ps_c(float _x) {
    float   t, f, p, r;
    int32_t i, j;

    const float l2e = 1.442695041f; /* log2(e) */
    const float l2h = -6.93145752e-1f; /* -log(2)_hi */
    const float l2l = -1.42860677e-6f; /* -log(2)_lo */
    /* coefficients for core approximation to exp() in [-log(2)/2, log(2)/2] */
    const float c0 = 0.008301110f;
    const float c1 = 0.041906696f;
    const float c2 = 0.166674897f;
    const float c3 = 0.499990642f;
    const float c4 = 0.999999762f;
    const float c5 = 1.000000000f;

    /* exp(x) = 2^i * e^f; i = rint (log2(e) * x), f = x - log(2) * i */
    t = _x * l2e; /* t = log2(e) * x */
    r = rintf(t);

    f = r * l2h + _x;
    f = r * l2l + f;

    i = (int32_t)rint(t);

    /* p ~= exp (f), -log(2)/2 <= f <= log(2)/2 */
    p = c0; /* c0 */

    p = p * f + c1;
    p = p * f + c2;
    p = p * f + c3;
    p = p * f + c4;
    p = p * f + c5;

    /* exp(x) = 2^i * p */
    j = i << 23; /* i << 23 */

    int32_t d;
    memcpy(&d, &p, sizeof(d));
    d += j;
    memcpy(&r, &d, sizeof(r));
    return r;
}

/***************************************************************************************************
* Applies temporal filter plane by plane.
* Inputs:
*   y_src, u_src, v_src : Pointers to the frame to be filtered, which is used as
*                    reference to compute squared differece from the predictor.
*   block_width: Width of the block.
*   block_height: Height of the block
*   noise_levels: Pointer to the noise levels of the to-filter frame, estimated
*                 with each plane (in Y, U, V order).
*   y_pre, r_pre, v_pre: Pointers to the well-built predictors.
*   accum: Pointer to the pixel-wise accumulator for filtering.
*   count: Pointer to the pixel-wise counter fot filtering.
* Returns:
*   Nothing will be returned. But the content to which `accum` and `pred`
*   point will be modified.
***************************************************************************************************/

static uint32_t calc_scaled_diff_fp4(uint32_t distance_fp4, uint32_t distance_threshold_fp16,
                                     uint32_t block_error_fp8, uint32_t tf_decay_factor_fp16,
                                     uint32_t num_ref_pixels, uint64_t sum_square_diff) {
    //const double d_factor = AOMMAX(distance / distance_threshold, 1);
    uint32_t d_factor_fp12 = AOMMAX(
        (uint32_t)(((((uint64_t)distance_fp4) << 24) + (distance_threshold_fp16 >> 1)) /
                   (distance_threshold_fp16)),
        (1 << 12));

    /* Float calculation:
    //Combine window error and block error, and normalize it.
    double window_error = (double)sum_square_diff / num_ref_pixels;
    double combined_error = (TF_WINDOW_BLOCK_BALANCE_WEIGHT * window_error + block_error) /
                (TF_WINDOW_BLOCK_BALANCE_WEIGHT + 1);
    double scaled_diff = AOMMIN(combined_error * d_factor / context_ptr->tf_decay_factor[0], 7);
    equivalent of:
    double scaled_diff = AOMMIN( (TF_WINDOW_BLOCK_BALANCE_WEIGHT * (sum_square_diff / num_ref_pixels) + block_error) /
                (TF_WINDOW_BLOCK_BALANCE_WEIGHT + 1) * d_factor / context_ptr->tf_decay_factor[0], 7);
    equivalent of:
    double scaled_diff_B = (d_factor * block_error /(TF_WINDOW_BLOCK_BALANCE_WEIGHT + 1))/context_ptr->tf_decay_factor[0];
    double scaled_diff_A1_factor = d_factor * TF_WINDOW_BLOCK_BALANCE_WEIGHT / context_ptr->tf_decay_factor[0];
    double scaled_diff_A2_factor = 1 /((TF_WINDOW_BLOCK_BALANCE_WEIGHT + 1) * num_ref_pixels);
    //No div internal in loop:
    double scaled_diff = AOMMIN(sum_square_diff * scaled_diff_A2_factor * scaled_diff_A1_factor + scaled_diff_B, 7);
    */

    uint32_t scaled_diff_B_fp4 = (uint32_t)(
        (((uint64_t)d_factor_fp12) *
         ((block_error_fp8 + ((TF_WINDOW_BLOCK_BALANCE_WEIGHT + 1) / 2)) /
          (TF_WINDOW_BLOCK_BALANCE_WEIGHT + 1))) /
        ((tf_decay_factor_fp16)));

    uint32_t scaled_diff_A1_factor_fp14 = (uint32_t)(
        ((((uint64_t)d_factor_fp12) << 16) * TF_WINDOW_BLOCK_BALANCE_WEIGHT +
         (tf_decay_factor_fp16 >> 3)) /
        (tf_decay_factor_fp16 >> 2));
    uint32_t scaled_diff_A2_factor_fp10 =
        ((1 << 10) + ((TF_WINDOW_BLOCK_BALANCE_WEIGHT + 1) * num_ref_pixels) / 2) /
        ((TF_WINDOW_BLOCK_BALANCE_WEIGHT + 1) * num_ref_pixels);

    FP_ASSERT(((((uint64_t)sum_square_diff)) * scaled_diff_A2_factor_fp10) < ((uint64_t)1 << 31));
    uint32_t sum_square_part_fp10 = (((uint32_t)sum_square_diff) * scaled_diff_A2_factor_fp10);

    /*Dynamic multiply to keep precision for various range of values*/
    int move = 0;
    if (scaled_diff_A1_factor_fp14 > (1 << 20)) {
        move = 18;
    } else if (scaled_diff_A1_factor_fp14 > (1 << 16)) {
        move = 14;
    } else if (scaled_diff_A1_factor_fp14 > (1 << 12)) {
        move = 10;
    } else if (scaled_diff_A1_factor_fp14 > (1 << 8)) {
        move = 6;
    } else if (scaled_diff_A1_factor_fp14 > (1 << 4)) {
        move = 2;
    } else {
        move = 0;
    }

    FP_ASSERT((((uint64_t)sum_square_part_fp10) * ((uint64_t)(scaled_diff_A1_factor_fp14) >> move) +
               ((uint64_t)1 << (19 - move))) < ((uint64_t)1 << 31));
    uint32_t scaled_diff_A_fp4 = (uint32_t)(
        ((sum_square_part_fp10) * ((scaled_diff_A1_factor_fp14) >> move) + (1 << (19 - move))) >>
        (20 - move));

    FP_ASSERT(scaled_diff_A_fp4 < ((uint32_t)1 << 31));
    uint32_t scaled_diff_fp4 = AOMMIN(scaled_diff_A_fp4 + scaled_diff_B_fp4, 7 * 16);
    return scaled_diff_fp4;
}

void svt_av1_apply_temporal_filter_planewise_c(
    struct MeContext *context_ptr, const uint8_t *y_src, int y_src_stride, const uint8_t *y_pre,
    int y_pre_stride, const uint8_t *u_src, const uint8_t *v_src, int uv_src_stride,
    const uint8_t *u_pre, const uint8_t *v_pre, int uv_pre_stride, unsigned int block_width,
    unsigned int block_height, int ss_x, int ss_y,
    uint32_t *y_accum, uint16_t *y_count, uint32_t *u_accum, uint16_t *u_count, uint32_t *v_accum,
    uint16_t *v_count) {
    unsigned int       i, j, k, m;
    int                idx, idy;
    uint64_t           sum_square_diff;
    const unsigned int uv_block_width  = block_width >> ss_x;
    const unsigned int uv_block_height = block_height >> ss_y;
    DECLARE_ALIGNED(16, uint16_t, y_diff_se[BLK_PELS]);
    DECLARE_ALIGNED(16, uint16_t, u_diff_se[BLK_PELS]);
    DECLARE_ALIGNED(16, uint16_t, v_diff_se[BLK_PELS]);

    memset(y_diff_se, 0, BLK_PELS * sizeof(uint16_t));
    memset(u_diff_se, 0, BLK_PELS * sizeof(uint16_t));
    memset(v_diff_se, 0, BLK_PELS * sizeof(uint16_t));

    // Calculate squared differences for each pixel of the block (pred-orig)
    calculate_squared_errors(
        y_src, y_src_stride, y_pre, y_pre_stride, y_diff_se, block_width, block_height);
    if (context_ptr->tf_chroma) {
        calculate_squared_errors(
            u_src, uv_src_stride, u_pre, uv_pre_stride, u_diff_se, uv_block_width, uv_block_height);
        calculate_squared_errors(
            v_src, uv_src_stride, v_pre, uv_pre_stride, v_diff_se, uv_block_width, uv_block_height);
    }

    // Get window size for pixel-wise filtering.
    assert(TF_PLANEWISE_FILTER_WINDOW_LENGTH % 2 == 1);
    const int half_window = TF_PLANEWISE_FILTER_WINDOW_LENGTH >> 1;
    for (i = 0; i < block_height; i++) {
        for (j = 0; j < block_width; j++) {
            const int pixel_value = y_pre[i * y_pre_stride + j];
            // non-local mean approach
            uint32_t num_ref_pixels = 0;

            const int uv_r  = i >> ss_y;
            const int uv_c  = j >> ss_x;
            sum_square_diff = 0;
            for (idy = -half_window; idy <= half_window; ++idy) {
                for (idx = -half_window; idx <= half_window; ++idx) {
                    const int row = CLIP((int)i + idy, 0, (int)block_height - 1);
                    const int col = CLIP((int)j + idx, 0, (int)block_width - 1);
                    sum_square_diff += y_diff_se[row * (int)block_width + col];
                    ++num_ref_pixels;
                }
            }
            // Combine window error and block error, and normalize it.
            int      idx_32x32, subblock_idx;
            uint32_t block_error_fp8 = 0;
            double   combined_error, window_error, block_error = 0.0f;
            if (context_ptr->tf_ctrls.use_fixed_point) {
                subblock_idx = (i >= block_height / 2) * 2 + (j >= block_width / 2);
                idx_32x32    = context_ptr->tf_block_col + context_ptr->tf_block_row * 2;
                if (context_ptr->tf_32x32_block_split_flag[idx_32x32]) {
                    // 16x16
                    FP_ASSERT(context_ptr->tf_16x16_block_error[idx_32x32 * 4 + subblock_idx] <
                              ((uint64_t)1 << 31));
                    //block_error = (double)context_ptr->tf_16x16_block_error[idx_32x32 * 4 + subblock_idx] / 256;
                    block_error_fp8 = (uint32_t)(
                        context_ptr->tf_16x16_block_error[idx_32x32 * 4 + subblock_idx]);
                } else {
                    //32x32
                    FP_ASSERT(context_ptr->tf_32x32_block_error[idx_32x32] < ((uint64_t)1 << 30));
                    //block_error = (double)context_ptr->tf_32x32_block_error[idx_32x32] / 1024;
                    block_error_fp8 = (uint32_t)(context_ptr->tf_32x32_block_error[idx_32x32] >> 2);
                }
            } else {
                window_error = (double)sum_square_diff / num_ref_pixels;

                subblock_idx = (i >= block_height / 2) * 2 + (j >= block_width / 2);
                idx_32x32    = context_ptr->tf_block_col + context_ptr->tf_block_row * 2;
                if (context_ptr->tf_32x32_block_split_flag[idx_32x32])
                    // 16x16
                    block_error =
                        (double)context_ptr->tf_16x16_block_error[idx_32x32 * 4 + subblock_idx] /
                        256;
                else
                    //32x32
                    block_error = (double)context_ptr->tf_32x32_block_error[idx_32x32] / 1024;

                combined_error = (TF_WINDOW_BLOCK_BALANCE_WEIGHT * window_error + block_error) /
                    (TF_WINDOW_BLOCK_BALANCE_WEIGHT + 1);
            }
            MV mv;
            if (context_ptr->tf_32x32_block_split_flag[idx_32x32]) {
                // 16x16
                mv.col = context_ptr->tf_16x16_mv_x[idx_32x32 * 4 + subblock_idx];
                mv.row = context_ptr->tf_16x16_mv_y[idx_32x32 * 4 + subblock_idx];
            } else {
                //32x32
                mv.col = context_ptr->tf_32x32_mv_x[idx_32x32];
                mv.row = context_ptr->tf_32x32_mv_y[idx_32x32];
            }
            uint32_t adjusted_weight, distance_threshold_fp16 = 0;
            uint32_t distance_fp4 = 0;
            double   d_factor     = 0.0f;
            if (context_ptr->tf_ctrls.use_fixed_point) {
                //const float  distance           = sqrtf(powf(mv.row, 2) + powf(mv.col, 2));
                distance_fp4 = sqrt_fast(
                    ((uint32_t)((int32_t)mv.col * mv.col + (int32_t)mv.row * mv.row)) << 8);
                //const double distance_threshold = (double)AOMMAX(
                //    context_ptr->min_frame_size * TF_SEARCH_DISTANCE_THRESHOLD, 1);
                distance_threshold_fp16 = AOMMAX((context_ptr->min_frame_size << 16) / 10, 1 << 16);
                uint32_t scaled_diff_fp4 = calc_scaled_diff_fp4(
                    distance_fp4,
                    distance_threshold_fp16,
                    block_error_fp8,
                    context_ptr->tf_decay_factor_fp16[0],
                    num_ref_pixels,
                    sum_square_diff);
                // Compute filter weight.
                adjusted_weight = (expf_tab_fp16[scaled_diff_fp4] * TF_WEIGHT_SCALE) >> 16;
            } else {
                const float  distance           = sqrtf(powf(mv.row, 2) + powf(mv.col, 2));
                const double distance_threshold = (double)AOMMAX(
                    context_ptr->min_frame_size * TF_SEARCH_DISTANCE_THRESHOLD, 1);
                d_factor = AOMMAX(distance / distance_threshold, 1);

                // Compute filter weight.
                double scaled_diff = AOMMIN(
                    combined_error * d_factor / context_ptr->tf_decay_factor[0], 7);

                adjusted_weight = (int)(exp_ps_c((float)(-scaled_diff)) * TF_WEIGHT_SCALE);
            }
            k = i * y_pre_stride + j;
            y_count[k] += adjusted_weight;
            y_accum[k] += adjusted_weight * pixel_value;

            // Process chroma component
            if (context_ptr->tf_chroma)
                if (!(i & ss_y) && !(j & ss_x)) {
                    const int u_pixel_value = u_pre[uv_r * uv_pre_stride + uv_c];
                    const int v_pixel_value = v_pre[uv_r * uv_pre_stride + uv_c];
                    // non-local mean approach
                    num_ref_pixels             = 0;
                    uint64_t u_sum_square_diff = 0, v_sum_square_diff = 0;
                    sum_square_diff = 0;
                    // Filter U-plane and V-plane using Y-plane. This is because motion
                    // search is only done on Y-plane, so the information from Y-plane will
                    // be more accurate.
                    for (idy = 0; idy < (1 << ss_y); ++idy) {
                        for (idx = 0; idx < (1 << ss_x); ++idx) {
                            const int row = (int)i + idy;
                            const int col = (int)j + idx;
                            sum_square_diff += y_diff_se[row * (int)block_width + col];
                            ++num_ref_pixels;
                        }
                    }
                    u_sum_square_diff = sum_square_diff;
                    v_sum_square_diff = sum_square_diff;

                    for (idy = -half_window; idy <= half_window; ++idy) {
                        for (idx = -half_window; idx <= half_window; ++idx) {
                            const int row = CLIP((int)uv_r + idy, 0, (int)uv_block_height - 1);
                            const int col = CLIP((int)uv_c + idx, 0, (int)uv_block_width - 1);
                            u_sum_square_diff += u_diff_se[row * uv_block_width + col];
                            v_sum_square_diff += v_diff_se[row * uv_block_width + col];
                            ++num_ref_pixels;
                        }
                    }

                    m = (i >> ss_y) * uv_pre_stride + (j >> ss_x);
                    // Combine window error and block error, and normalize it.
                    if (context_ptr->tf_ctrls.use_fixed_point) {
                        uint32_t scaled_diff_fp4 = calc_scaled_diff_fp4(
                            distance_fp4,
                            distance_threshold_fp16,
                            block_error_fp8,
                            context_ptr->tf_decay_factor_fp16[1],
                            num_ref_pixels,
                            u_sum_square_diff);
                        // Compute filter weight.
                        adjusted_weight = (expf_tab_fp16[scaled_diff_fp4] * TF_WEIGHT_SCALE) >> 16;
                    } else {
                        window_error   = (double)u_sum_square_diff / num_ref_pixels;
                        combined_error = (TF_WINDOW_BLOCK_BALANCE_WEIGHT * window_error +
                                          block_error) /
                            (TF_WINDOW_BLOCK_BALANCE_WEIGHT + 1);

                        // Compute filter weight.
                        double scaled_diff = AOMMIN(
                            combined_error * d_factor / context_ptr->tf_decay_factor[1], 7);

                        adjusted_weight = (int)(exp_ps_c((float)(-scaled_diff)) * TF_WEIGHT_SCALE);
                    }
                    u_count[m] += adjusted_weight;
                    u_accum[m] += adjusted_weight * u_pixel_value;

                    if (context_ptr->tf_ctrls.use_fixed_point) {
                        uint32_t scaled_diff_fp4 = calc_scaled_diff_fp4(
                            distance_fp4,
                            distance_threshold_fp16,
                            block_error_fp8,
                            context_ptr->tf_decay_factor_fp16[2],
                            num_ref_pixels,
                            v_sum_square_diff);
                        // Compute filter weight.
                        adjusted_weight = (expf_tab_fp16[scaled_diff_fp4] * TF_WEIGHT_SCALE) >> 16;
                    } else {
                        // Combine window error and block error, and normalize it.
                        window_error   = (double)v_sum_square_diff / num_ref_pixels;
                        combined_error = (TF_WINDOW_BLOCK_BALANCE_WEIGHT * window_error +
                                          block_error) /
                            (TF_WINDOW_BLOCK_BALANCE_WEIGHT + 1);

                        // Compute filter weight.
                        double scaled_diff = AOMMIN(
                            combined_error * d_factor / context_ptr->tf_decay_factor[2], 7);

                        adjusted_weight = (int)(exp_ps_c((float)(-scaled_diff)) * TF_WEIGHT_SCALE);
                    }
                    v_count[m] += adjusted_weight;
                    v_accum[m] += adjusted_weight * v_pixel_value;
                }
        }
    }
}

/***************************************************************************************************
* Applies temporal filter plane by plane for hbd
* Inputs:
*   y_src, u_src, v_src : Pointers to the frame to be filtered, which is used as
*                    reference to compute squared differece from the predictor.
*   block_width: Width of the block.
*   block_height: Height of the block
*   noise_levels: Pointer to the noise levels of the to-filter frame, estimated
*                 with each plane (in Y, U, V order).
*   y_pre, r_pre, v_pre: Pointers to the well-built predictors.
*   accum: Pointer to the pixel-wise accumulator for filtering.
*   count: Pointer to the pixel-wise counter fot filtering.
* Returns:
*   Nothing will be returned. But the content to which `accum` and `pred`
*   point will be modified.
***************************************************************************************************/
void svt_av1_apply_temporal_filter_planewise_hbd_c(
    struct MeContext *context_ptr, const uint16_t *y_src, int y_src_stride, const uint16_t *y_pre,
    int y_pre_stride, const uint16_t *u_src, const uint16_t *v_src, int uv_src_stride,
    const uint16_t *u_pre, const uint16_t *v_pre, int uv_pre_stride, unsigned int block_width,
    unsigned int block_height, int ss_x, int ss_y,
    uint32_t *y_accum, uint16_t *y_count, uint32_t *u_accum, uint16_t *u_count, uint32_t *v_accum,
    uint16_t *v_count, uint32_t encoder_bit_depth) {
    unsigned int       i, j, k, m;
    int                idx, idy;
    uint64_t           sum_square_diff;
    const unsigned int uv_block_width  = block_width >> ss_x;
    const unsigned int uv_block_height = block_height >> ss_y;
    DECLARE_ALIGNED(16, uint32_t, y_diff_se[BLK_PELS]);
    DECLARE_ALIGNED(16, uint32_t, u_diff_se[BLK_PELS]);
    DECLARE_ALIGNED(16, uint32_t, v_diff_se[BLK_PELS]);

    memset(y_diff_se, 0, BLK_PELS * sizeof(uint32_t));
    memset(u_diff_se, 0, BLK_PELS * sizeof(uint32_t));
    memset(v_diff_se, 0, BLK_PELS * sizeof(uint32_t));

    // Calculate squared differences for each pixel of the block (pred-orig)
    calculate_squared_errors_highbd(
        y_src, y_src_stride, y_pre, y_pre_stride, y_diff_se, block_width, block_height);
    if (context_ptr->tf_chroma) {
        calculate_squared_errors_highbd(
            u_src, uv_src_stride, u_pre, uv_pre_stride, u_diff_se, uv_block_width, uv_block_height);
        calculate_squared_errors_highbd(
            v_src, uv_src_stride, v_pre, uv_pre_stride, v_diff_se, uv_block_width, uv_block_height);
    }
    // Get window size for pixel-wise filtering.
    assert(TF_PLANEWISE_FILTER_WINDOW_LENGTH % 2 == 1);
    const int half_window = TF_PLANEWISE_FILTER_WINDOW_LENGTH >> 1;
    for (i = 0; i < block_height; i++) {
        for (j = 0; j < block_width; j++) {
            const int pixel_value = y_pre[i * y_pre_stride + j];
            // non-local mean approach
            uint32_t num_ref_pixels = 0;

            const int uv_r  = i >> ss_y;
            const int uv_c  = j >> ss_x;
            sum_square_diff = 0;
            for (idy = -half_window; idy <= half_window; ++idy) {
                for (idx = -half_window; idx <= half_window; ++idx) {
                    const int row = CLIP((int)i + idy, 0, (int)block_height - 1);
                    const int col = CLIP((int)j + idx, 0, (int)block_width - 1);
                    sum_square_diff += y_diff_se[row * (int)block_width + col];
                    ++num_ref_pixels;
                }
            }
            // Combine window error and block error, and normalize it.
            // Scale down the difference for high bit depth input.
            sum_square_diff >>= ((encoder_bit_depth - 8) * 2);
            int      subblock_idx    = (i >= block_height / 2) * 2 + (j >= block_width / 2);
            int      idx_32x32       = context_ptr->tf_block_col + context_ptr->tf_block_row * 2;
            uint32_t block_error_fp8 = 0;
            double   combined_error, window_error, block_error = 0.0;
            if (context_ptr->tf_ctrls.use_fixed_point) {
                if (context_ptr->tf_32x32_block_split_flag[idx_32x32]) {
                    // 16x16
                    // Scale down the difference for high bit depth input.
                    FP_ASSERT(context_ptr->tf_16x16_block_error[idx_32x32 * 4 + subblock_idx] <
                              ((uint64_t)1 << 30));
                    //block_error =
                    //    (double)(context_ptr->tf_16x16_block_error[idx_32x32 * 4 + subblock_idx] >> 4) / 256;
                    block_error_fp8 = (uint32_t)(
                        context_ptr->tf_16x16_block_error[idx_32x32 * 4 + subblock_idx] >> 4);
                } else {
                    //32x32
                    // Scale down the difference for high bit depth input.
                    //block_error = (double)(context_ptr->tf_32x32_block_error[idx_32x32] >> 4) / 1024;
                    FP_ASSERT(context_ptr->tf_32x32_block_error[idx_32x32] < ((uint64_t)1 << 30));
                    block_error_fp8 = (uint32_t)(context_ptr->tf_32x32_block_error[idx_32x32] >> 6);
                }
            } else {
                window_error = (double)sum_square_diff / num_ref_pixels;
                if (context_ptr->tf_32x32_block_split_flag[idx_32x32])
                    // 16x16
                    // Scale down the difference for high bit depth input.
                    block_error =
                        (double)(context_ptr->tf_16x16_block_error[idx_32x32 * 4 + subblock_idx] >>
                                 4) /
                        256;
                else
                    //32x32
                    // Scale down the difference for high bit depth input.
                    block_error = (double)(context_ptr->tf_32x32_block_error[idx_32x32] >> 4) /
                        1024;

                combined_error = (TF_WINDOW_BLOCK_BALANCE_WEIGHT * window_error + block_error) /
                    (TF_WINDOW_BLOCK_BALANCE_WEIGHT + 1);
            }
            MV mv;
            if (context_ptr->tf_32x32_block_split_flag[idx_32x32]) {
                // 16x16
                mv.col = context_ptr->tf_16x16_mv_x[idx_32x32 * 4 + subblock_idx];
                mv.row = context_ptr->tf_16x16_mv_y[idx_32x32 * 4 + subblock_idx];
            } else {
                //32x32
                mv.col = context_ptr->tf_32x32_mv_x[idx_32x32];
                mv.row = context_ptr->tf_32x32_mv_y[idx_32x32];
            }
            uint32_t adjusted_weight, distance_threshold_fp16 = 0;
            uint32_t scaled_diff_fp4;
            uint32_t distance_fp4 = 0;
            double   scaled_diff, d_factor = 0.0;
            if (context_ptr->tf_ctrls.use_fixed_point) {
                distance_fp4 = sqrt_fast(
                    ((uint32_t)((int32_t)mv.col * mv.col + (int32_t)mv.row * mv.row)) << 8);
                distance_threshold_fp16 = AOMMAX((context_ptr->min_frame_size << 16) / 10, 1 << 16);
                scaled_diff_fp4         = calc_scaled_diff_fp4(distance_fp4,
                                                       distance_threshold_fp16,
                                                       block_error_fp8,
                                                       context_ptr->tf_decay_factor_fp16[0],
                                                       num_ref_pixels,
                                                       sum_square_diff);
                // Compute filter weight.
                adjusted_weight = (expf_tab_fp16[scaled_diff_fp4] * TF_WEIGHT_SCALE) >> 16;
            } else {
                const float  distance           = sqrtf(powf(mv.row, 2) + powf(mv.col, 2));
                const double distance_threshold = (double)AOMMAX(
                    context_ptr->min_frame_size * TF_SEARCH_DISTANCE_THRESHOLD, 1);
                d_factor = AOMMAX(distance / distance_threshold, 1);

                // Compute filter weight.
                scaled_diff = AOMMIN(combined_error * d_factor / context_ptr->tf_decay_factor[0],
                                     7);

                adjusted_weight = (int)(exp_ps_c((float)-scaled_diff) * TF_WEIGHT_SCALE);
            }
            k = i * y_pre_stride + j;
            y_count[k] += adjusted_weight;
            y_accum[k] += adjusted_weight * pixel_value;

            // Process chroma component
            if (context_ptr->tf_chroma)
                if (!(i & ss_y) && !(j & ss_x)) {
                    const int u_pixel_value = u_pre[uv_r * uv_pre_stride + uv_c];
                    const int v_pixel_value = v_pre[uv_r * uv_pre_stride + uv_c];
                    // non-local mean approach
                    num_ref_pixels             = 0;
                    uint64_t u_sum_square_diff = 0, v_sum_square_diff = 0;
                    sum_square_diff = 0;
                    // Filter U-plane and V-plane using Y-plane. This is because motion
                    // search is only done on Y-plane, so the information from Y-plane will
                    // be more accurate.
                    for (idy = 0; idy < (1 << ss_y); ++idy) {
                        for (idx = 0; idx < (1 << ss_x); ++idx) {
                            const int row = (int)i + idy;
                            const int col = (int)j + idx;
                            sum_square_diff += y_diff_se[row * (int)block_width + col];
                            ++num_ref_pixels;
                        }
                    }
                    u_sum_square_diff = sum_square_diff;
                    v_sum_square_diff = sum_square_diff;

                    for (idy = -half_window; idy <= half_window; ++idy) {
                        for (idx = -half_window; idx <= half_window; ++idx) {
                            const int row = CLIP((int)uv_r + idy, 0, (int)uv_block_height - 1);
                            const int col = CLIP((int)uv_c + idx, 0, (int)uv_block_width - 1);
                            u_sum_square_diff += u_diff_se[row * uv_block_width + col];
                            v_sum_square_diff += v_diff_se[row * uv_block_width + col];
                            ++num_ref_pixels;
                        }
                    }

                    m = (i >> ss_y) * uv_pre_stride + (j >> ss_x);
                    // Scale down the difference for high bit depth input.
                    u_sum_square_diff >>= ((encoder_bit_depth - 8) * 2);
                    v_sum_square_diff >>= ((encoder_bit_depth - 8) * 2);
                    // Combine window error and block error, and normalize it.
                    if (context_ptr->tf_ctrls.use_fixed_point) {
                        scaled_diff_fp4 = calc_scaled_diff_fp4(distance_fp4,
                                                               distance_threshold_fp16,
                                                               block_error_fp8,
                                                               context_ptr->tf_decay_factor_fp16[1],
                                                               num_ref_pixels,
                                                               u_sum_square_diff);
                        // Compute filter weight.
                        adjusted_weight = (expf_tab_fp16[scaled_diff_fp4] * TF_WEIGHT_SCALE) >> 16;
                    } else {
                        window_error   = (double)u_sum_square_diff / num_ref_pixels;
                        combined_error = (TF_WINDOW_BLOCK_BALANCE_WEIGHT * window_error +
                                          block_error) /
                            (TF_WINDOW_BLOCK_BALANCE_WEIGHT + 1);

                        // Compute filter weight.
                        scaled_diff = AOMMIN(
                            combined_error * d_factor / context_ptr->tf_decay_factor[1], 7);

                        adjusted_weight = (int)(exp_ps_c((float)-scaled_diff) * TF_WEIGHT_SCALE);
                    }
                    u_count[m] += adjusted_weight;
                    u_accum[m] += adjusted_weight * u_pixel_value;

                    if (context_ptr->tf_ctrls.use_fixed_point) {
                        scaled_diff_fp4 = calc_scaled_diff_fp4(distance_fp4,
                                                               distance_threshold_fp16,
                                                               block_error_fp8,
                                                               context_ptr->tf_decay_factor_fp16[2],
                                                               num_ref_pixels,
                                                               v_sum_square_diff);
                        // Compute filter weight.
                        adjusted_weight = (expf_tab_fp16[scaled_diff_fp4] * TF_WEIGHT_SCALE) >> 16;
                    } else {
                        // Combine window error and block error, and normalize it.
                        window_error   = (double)v_sum_square_diff / num_ref_pixels;
                        combined_error = (TF_WINDOW_BLOCK_BALANCE_WEIGHT * window_error +
                                          block_error) /
                            (TF_WINDOW_BLOCK_BALANCE_WEIGHT + 1);

                        // Compute filter weight.
                        scaled_diff = AOMMIN(
                            combined_error * d_factor / context_ptr->tf_decay_factor[2], 7);

                        adjusted_weight = (int)(exp_ps_c((float)-scaled_diff) * TF_WEIGHT_SCALE);
                    }
                    v_count[m] += adjusted_weight;
                    v_accum[m] += adjusted_weight * v_pixel_value;
                }
        }
    }
}
/* calculates SSE*/
uint32_t calculate_squared_errors_sum(const uint8_t *s, int s_stride, const uint8_t *p,
                                      int p_stride, unsigned int w, unsigned int h) {
    unsigned int i, j;
    uint32_t     sum = 0;
    for (i = 0; i < h; i++) {
        for (j = 0; j < w; j++) {
            const int16_t diff = s[i * s_stride + j] - p[i * p_stride + j];
            sum += diff * diff;
        }
    }
    return sum;
}

/* calculates SSE for 10bit*/
uint32_t calculate_squared_errors_sum_highbd(const uint16_t *s, int s_stride, const uint16_t *p,
                                             int p_stride,
                                             unsigned int w, unsigned int h, int shift_factor) {
    unsigned int i, j;
    uint32_t     sum = 0;
    for (i = 0; i < h; i++) {
        for (j = 0; j < w; j++) {
            const int32_t diff = s[i * s_stride + j] - p[i * p_stride + j];
            sum += diff * diff;
        }
    }
    return (sum >> shift_factor);
}
/*
apply fast TF filter
*/
void svt_av1_apply_temporal_filter_planewise_fast_c(struct MeContext *context_ptr,
                                                    const uint8_t *y_src, int y_src_stride,
                                                    const uint8_t *y_pre, int y_pre_stride,
                                                    unsigned int block_width,
                                                    unsigned int block_height, uint32_t *y_accum,
                                                    uint16_t *y_count) {
    //to add chroma if need be
    uint32_t avg_err = calculate_squared_errors_sum(
                           y_src, y_src_stride, y_pre, y_pre_stride, block_width, block_height) /
        (block_width * block_height);

    uint32_t adjusted_weight;
    if (context_ptr->tf_ctrls.use_fixed_point) {
        //16*avg_err/context_ptr->tf_decay_factor[0];
        uint32_t scaled_diff_fp4 = AOMMIN(
            (avg_err << 10) / (AOMMAX(context_ptr->tf_decay_factor_fp16[0] >> 10, 1)), 7 * 16);
        adjusted_weight = (expf_tab_fp16[scaled_diff_fp4] * TF_WEIGHT_SCALE) >> 16;
    } else {
        double scaled_diff = AOMMIN(avg_err / context_ptr->tf_decay_factor[0], 7);

        adjusted_weight = (int)(expf_tab[(int)(scaled_diff * 10)] * TF_WEIGHT_SCALE);
    }

    if (adjusted_weight) {
        unsigned int i, j, k;
        for (i = 0; i < block_height; i++) {
            for (j = 0; j < block_width; j++) {
                const int pixel_value = y_pre[i * y_pre_stride + j];
                k                     = i * y_pre_stride + j;
                y_count[k] += adjusted_weight;
                y_accum[k] += adjusted_weight * pixel_value;
            }
        }
    }
}
/*
apply fast TF filter for 10bit
*/
void svt_av1_apply_temporal_filter_planewise_fast_hbd_c(
    struct MeContext *context_ptr, const uint16_t *y_src, int y_src_stride, const uint16_t *y_pre,
    int y_pre_stride, unsigned int block_width,
    unsigned int block_height, uint32_t *y_accum, uint16_t *y_count, uint32_t encoder_bit_depth) {

    int shift_factor = ((encoder_bit_depth - 8) * 2);
    //to add chroma if need be
    uint32_t avg_err =
        calculate_squared_errors_sum_highbd(
            y_src, y_src_stride, y_pre, y_pre_stride, block_width, block_height, shift_factor) /
        (block_width * block_height);

    uint32_t adjusted_weight;
    if (context_ptr->tf_ctrls.use_fixed_point) {
        //16*avg_err/context_ptr->tf_decay_factor[0];
        uint32_t scaled_diff_fp4 = AOMMIN(
            (avg_err << 10) / (AOMMAX(context_ptr->tf_decay_factor_fp16[0] >> 10, 1)), 7 * 16);
        adjusted_weight = (expf_tab_fp16[scaled_diff_fp4] * TF_WEIGHT_SCALE) >> 16;
    } else {
        double scaled_diff = AOMMIN(avg_err / context_ptr->tf_decay_factor[0], 7);
        adjusted_weight = (int)(expf_tab[(int)(scaled_diff * 10)] * TF_WEIGHT_SCALE);
    }
    if (adjusted_weight) {
        unsigned int i, j, k;
        for (i = 0; i < block_height; i++) {
            for (j = 0; j < block_width; j++) {
                const int pixel_value = y_pre[i * y_pre_stride + j];
                k                     = i * y_pre_stride + j;
                y_count[k] += adjusted_weight;
                y_accum[k] += adjusted_weight * pixel_value;
            }
        }
    }
}

/***************************************************************************************************
* Applies temporal filter plane by plane.
* Inputs:
*   y_src, u_src, v_src : Pointers to the frame to be filtered, which is used as
*                    reference to compute squared differece from the predictor.
*   block_width: Width of the block.
*   block_height: Height of the block
*   noise_levels: Pointer to the noise levels of the to-filter frame, estimated
*                 with each plane (in Y, U, V order).
*   y_pre, r_pre, v_pre: Pointers to the well-built predictors.
*   accum: Pointer to the pixel-wise accumulator for filtering.
*   count: Pointer to the pixel-wise counter fot filtering.
* Returns:
*   Nothing will be returned. But the content to which `accum` and `pred`
*   point will be modified.
***************************************************************************************************/
static void svt_av1_apply_temporal_filter_planewise_medium_partial_c(
    struct MeContext *context_ptr, const uint8_t *y_src, int y_src_stride, const uint8_t *y_pre,
    int y_pre_stride, unsigned int block_width, unsigned int block_height, uint32_t *y_accum,
    uint16_t *y_count, const uint32_t tf_decay_factor_fp16, uint32_t luma_window_error_quad_fp8[4],
    int is_chroma) {
    unsigned int i, j, subblock_idx;

    // Decay factors for non-local mean approach.
    // Larger noise -> larger filtering weight.
    int32_t idx_32x32 = context_ptr->tf_block_col + context_ptr->tf_block_row * 2;

    //double distance_threshold = (double)AOMMAX(context_ptr->min_frame_size * TF_SEARCH_DISTANCE_THRESHOLD, 1);
    uint32_t distance_threshold_fp16 = AOMMAX((context_ptr->min_frame_size << 16) / 10, 1 << 16);

    //Calculation for every quarter
    uint32_t  d_factor_fp8[4];
    uint32_t  block_error_fp8[4];
    uint32_t  chroma_window_error_quad_fp8[4];
    uint32_t *window_error_quad_fp8 = is_chroma ? chroma_window_error_quad_fp8
                                                : luma_window_error_quad_fp8;

    if (context_ptr->tf_32x32_block_split_flag[idx_32x32]) {
        for (i = 0; i < 4; ++i) {
            int32_t col = context_ptr->tf_16x16_mv_x[idx_32x32 * 4 + i];
            int32_t row = context_ptr->tf_16x16_mv_y[idx_32x32 * 4 + i];
            //const float  distance = sqrtf((float)col*col + row*row);
            uint32_t distance_fp4 = sqrt_fast(((uint32_t)(col * col + row * row)) << 8);
            d_factor_fp8[i] = AOMMAX((distance_fp4 << 12) / (distance_threshold_fp16 >> 8), 1 << 8);
            FP_ASSERT(context_ptr->tf_16x16_block_error[idx_32x32 * 4 + i] < ((uint64_t)1 << 31));
            //block_error[i] = (double)context_ptr->tf_16x16_block_error[idx_32x32 * 4 + i] / 256;
            block_error_fp8[i] = (uint32_t)(context_ptr->tf_16x16_block_error[idx_32x32 * 4 + i]);
        }
    } else {
        int32_t col = context_ptr->tf_32x32_mv_x[idx_32x32];
        int32_t row = context_ptr->tf_32x32_mv_y[idx_32x32];

        uint32_t distance_fp4 = sqrt_fast(((uint32_t)(col * col + row * row)) << 8);
        //d_factor[0] = d_factor[1] = d_factor[2] = d_factor[3] = AOMMAX(distance / distance_threshold, 1);
        d_factor_fp8[0] = d_factor_fp8[1] = d_factor_fp8[2] = d_factor_fp8[3] = AOMMAX(
            (distance_fp4 << 12) / (distance_threshold_fp16 >> 8), 1 << 8);
        FP_ASSERT(context_ptr->tf_32x32_block_error[idx_32x32] < ((uint64_t)1 << 30));
        //block_error[0] = block_error[1] = block_error[2] = block_error[3] = (double)context_ptr->tf_32x32_block_error[idx_32x32] / 1024;
        block_error_fp8[0] = block_error_fp8[1] = block_error_fp8[2] = block_error_fp8[3] =
            (uint32_t)(context_ptr->tf_32x32_block_error[idx_32x32] >> 2);
    }
    const uint32_t bw_half = (block_width >> 1);
    const uint32_t bh_half = (block_height >> 1);
    uint32_t       sum;
    sum = calculate_squared_errors_sum(y_src, y_src_stride, y_pre, y_pre_stride, bw_half, bh_half);
    FP_ASSERT(sum <= (1 << (26)));
    window_error_quad_fp8[0] = (((sum << 4) / bw_half) << 4) / bh_half;

    sum = calculate_squared_errors_sum(
        y_src + bw_half, y_src_stride, y_pre + bw_half, y_pre_stride, bw_half, bh_half);
    FP_ASSERT(sum <= (1 << (26)));
    window_error_quad_fp8[1] = (((sum << 4) / bw_half) << 4) / bh_half;

    sum = calculate_squared_errors_sum(y_src + y_src_stride * (bh_half),
                                       y_src_stride,
                                       y_pre + y_pre_stride * (bh_half),
                                       y_pre_stride,
                                       bw_half,
                                       bh_half);
    FP_ASSERT(sum <= (1 << (26)));
    window_error_quad_fp8[2] = (((sum << 4) / bw_half) << 4) / bh_half;

    sum = calculate_squared_errors_sum(y_src + y_src_stride * (bh_half) + bw_half,
                                       y_src_stride,
                                       y_pre + y_pre_stride * (bh_half) + bw_half,
                                       y_pre_stride,
                                       bw_half,
                                       bh_half);
    FP_ASSERT(sum <= (1 << (26)));
    window_error_quad_fp8[3] = (((sum << 4) / bw_half) << 4) / bh_half;

    if (is_chroma) {
        for (i = 0; i < 4; ++i) {
            FP_ASSERT(((int64_t)window_error_quad_fp8[i] * 5 + luma_window_error_quad_fp8[i]) <
                      ((int64_t)1 << 31));
            window_error_quad_fp8[i] = (window_error_quad_fp8[i] * 5 +
                                        luma_window_error_quad_fp8[i]) /
                6;
        }
    }

    //Calculation for every quarter
    for (subblock_idx = 0; subblock_idx < 4; subblock_idx++) {
        FP_ASSERT(((int64_t)window_error_quad_fp8[subblock_idx] * TF_WINDOW_BLOCK_BALANCE_WEIGHT +
                   block_error_fp8[subblock_idx]) < ((int64_t)1 << 31));
        uint32_t combined_error_fp8 = (window_error_quad_fp8[subblock_idx] *
                                           TF_WINDOW_BLOCK_BALANCE_WEIGHT +
                                       block_error_fp8[subblock_idx]) /
            (TF_WINDOW_BLOCK_BALANCE_WEIGHT + 1);

        uint32_t avg_err_fp10 = ((combined_error_fp8 >> 3) * (d_factor_fp8[subblock_idx] >> 3));
        // int32_t avg_err_fp12 = ((int64_t)(combined_error_fp16) * d_factor_fp16[subblock_idx])>>20;
        FP_ASSERT((((int64_t)combined_error_fp8 >> 3) * (d_factor_fp8[subblock_idx] >> 3)) <
                  ((int64_t)1 << 31));
        //double scaled_diff = AOMMIN(combined_error * d_factor[subblock_idx] / (FP2FLOAT(tf_decay_factor_fp16)), 7);
        uint32_t scaled_diff16 = AOMMIN(
            /*((16*avg_err)<<8)*/ (avg_err_fp10) / AOMMAX((tf_decay_factor_fp16 >> 10), 1), 7 * 16);
        //int adjusted_weight = (int)(expf((float)(-scaled_diff)) * TF_WEIGHT_SCALE);
        uint32_t adjusted_weight = (expf_tab_fp16[scaled_diff16] * TF_WEIGHT_SCALE) >> 16;

        int x_offset = (subblock_idx % 2) * block_width / 2;
        int y_offset = (subblock_idx / 2) * block_height / 2;

        for (i = 0; i < block_height / 2; i++) {
            for (j = 0; j < block_width / 2; j++) {
                const int k           = (i + y_offset) * y_pre_stride + j + x_offset;
                const int pixel_value = y_pre[k];
                y_count[k] += adjusted_weight;
                y_accum[k] += adjusted_weight * pixel_value;
            }
        }
    }
}

void svt_av1_apply_temporal_filter_planewise_medium_c(
    struct MeContext *context_ptr, const uint8_t *y_src, int y_src_stride, const uint8_t *y_pre,
    int y_pre_stride, const uint8_t *u_src, const uint8_t *v_src, int uv_src_stride,
    const uint8_t *u_pre, const uint8_t *v_pre, int uv_pre_stride, unsigned int block_width,
    unsigned int block_height, int ss_x, int ss_y, uint32_t *y_accum, uint16_t *y_count,
    uint32_t *u_accum, uint16_t *u_count, uint32_t *v_accum, uint16_t *v_count) {
    uint32_t luma_window_error_quad_fp8[4];

    svt_av1_apply_temporal_filter_planewise_medium_partial_c(context_ptr,
                                                             y_src,
                                                             y_src_stride,
                                                             y_pre,
                                                             y_pre_stride,
                                                             (unsigned int)block_width,
                                                             (unsigned int)block_height,
                                                             y_accum,
                                                             y_count,
                                                             context_ptr->tf_decay_factor_fp16[C_Y],
                                                             luma_window_error_quad_fp8,
                                                             0);
    if (context_ptr->tf_chroma) {
        svt_av1_apply_temporal_filter_planewise_medium_partial_c(
            context_ptr,
            u_src,
            uv_src_stride,
            u_pre,
            uv_pre_stride,
            (unsigned int)block_width >> ss_x,
            (unsigned int)block_height >> ss_y,
            u_accum,
            u_count,
            context_ptr->tf_decay_factor_fp16[C_U],
            luma_window_error_quad_fp8,
            1);

        svt_av1_apply_temporal_filter_planewise_medium_partial_c(
            context_ptr,
            v_src,
            uv_src_stride,
            v_pre,
            uv_pre_stride,
            (unsigned int)block_width >> ss_x,
            (unsigned int)block_height >> ss_y,
            v_accum,
            v_count,
            context_ptr->tf_decay_factor_fp16[C_V],
            luma_window_error_quad_fp8,
            1);
    }
}

/***************************************************************************************************
* Applies temporal filter plane by plane for hbd
* Inputs:
*   y_src, u_src, v_src : Pointers to the frame to be filtered, which is used as
*                    reference to compute squared differece from the predictor.
*   block_width: Width of the block.
*   block_height: Height of the block
*   noise_levels: Pointer to the noise levels of the to-filter frame, estimated
*                 with each plane (in Y, U, V order).
*   y_pre, r_pre, v_pre: Pointers to the well-built predictors.
*   accum: Pointer to the pixel-wise accumulator for filtering.
*   count: Pointer to the pixel-wise counter fot filtering.
* Returns:
*   Nothing will be returned. But the content to which `accum` and `pred`
*   point will be modified.
***************************************************************************************************/
static void svt_av1_apply_temporal_filter_planewise_medium_hbd_partial_c(
    struct MeContext *context_ptr, const uint16_t *y_src, int y_src_stride, const uint16_t *y_pre,
    int y_pre_stride, unsigned int block_width, unsigned int block_height, uint32_t *y_accum,
    uint16_t *y_count, const uint32_t tf_decay_factor_fp16, uint32_t luma_window_error_quad_fp8[4],
    int is_chroma, uint32_t encoder_bit_depth) {
    unsigned int i, j, subblock_idx;
    // Decay factors for non-local mean approach.
    // Larger noise -> larger filtering weight.
    int idx_32x32    = context_ptr->tf_block_col + context_ptr->tf_block_row * 2;
    int shift_factor = ((encoder_bit_depth - 8) * 2);
    //const double distance_threshold = (double)AOMMAX(context_ptr->min_frame_size * TF_SEARCH_DISTANCE_THRESHOLD, 1);
    uint32_t distance_threshold_fp16 = AOMMAX((context_ptr->min_frame_size << 16) / 10, 1 << 16);

    //Calculation for every quarter
    uint32_t  d_factor_fp8[4];
    uint32_t  block_error_fp8[4];
    uint32_t  chroma_window_error_quad_fp8[4];
    uint32_t *window_error_quad_fp8 = is_chroma ? chroma_window_error_quad_fp8
                                                : luma_window_error_quad_fp8;

    if (context_ptr->tf_32x32_block_split_flag[idx_32x32]) {
        for (i = 0; i < 4; ++i) {
            int32_t col = context_ptr->tf_16x16_mv_x[idx_32x32 * 4 + i];
            int32_t row = context_ptr->tf_16x16_mv_y[idx_32x32 * 4 + i];
            //const float  distance = sqrtf((float)col*col + row*row);
            uint32_t distance_fp4 = sqrt_fast(((uint32_t)(col * col + row * row)) << 8);
            //d_factor[i] = AOMMAX(distance / distance_threshold, 1);
            d_factor_fp8[i] = AOMMAX((distance_fp4 << 12) / (distance_threshold_fp16 >> 8), 1 << 8);
            FP_ASSERT(context_ptr->tf_16x16_block_error[idx_32x32 * 4 + i] < ((uint64_t)1 << 35));
            //block_error[i] = (double)(context_ptr->tf_16x16_block_error[idx_32x32 * 4 + i] >>4) / 256;
            block_error_fp8[i] = (uint32_t)(context_ptr->tf_16x16_block_error[idx_32x32 * 4 + i] >>
                                            4);
        }
    } else {
        int32_t col = context_ptr->tf_32x32_mv_x[idx_32x32];
        int32_t row = context_ptr->tf_32x32_mv_y[idx_32x32];
        //const float  distance = sqrtf((float)col*col + row*row);
        uint32_t distance_fp4 = sqrt_fast(((uint32_t)(col * col + row * row)) << 8);
        //d_factor[i] = AOMMAX(distance / distance_threshold, 1);
        d_factor_fp8[0] = d_factor_fp8[1] = d_factor_fp8[2] = d_factor_fp8[3] = AOMMAX(
            (distance_fp4 << 12) / (distance_threshold_fp16 >> 8), 1 << 8);
        FP_ASSERT(context_ptr->tf_32x32_block_error[idx_32x32] < ((uint64_t)1 << 35));
        //= (double)(context_ptr->tf_32x32_block_error[idx_32x32]>> 4) / 1024;
        block_error_fp8[0] = block_error_fp8[1] = block_error_fp8[2] = block_error_fp8[3] =
            (uint32_t)(context_ptr->tf_32x32_block_error[idx_32x32] >> 6);
    }
    const uint32_t bw_half = (block_width >> 1);
    const uint32_t bh_half = (block_height >> 1);
    uint32_t       sum;

    sum = calculate_squared_errors_sum_highbd(
        y_src, y_src_stride, y_pre, y_pre_stride, bw_half, bh_half, shift_factor);
    FP_ASSERT(sum <= (1 << (26)));
    window_error_quad_fp8[0] = (((sum << 4) / bw_half) << 4) / bh_half;

    sum = calculate_squared_errors_sum_highbd(y_src + bw_half,
                                              y_src_stride,
                                              y_pre + bw_half,
                                              y_pre_stride,
                                              bw_half,
                                              bh_half,
                                              shift_factor);
    FP_ASSERT(sum <= (1 << (26)));
    window_error_quad_fp8[1] = (((sum << 4) / bw_half) << 4) / bh_half;
    sum                      = calculate_squared_errors_sum_highbd(y_src + y_src_stride * (bh_half),
                                              y_src_stride,
                                              y_pre + y_pre_stride * (bh_half),
                                              y_pre_stride,
                                              bw_half,
                                              bh_half,
                                              shift_factor);
    FP_ASSERT(sum <= (1 << (26)));
    window_error_quad_fp8[2] = (((sum << 4) / bw_half) << 4) / bh_half;
    sum = calculate_squared_errors_sum_highbd(y_src + y_src_stride * (bh_half) + bw_half,
                                              y_src_stride,
                                              y_pre + y_pre_stride * (bh_half) + bw_half,
                                              y_pre_stride,
                                              bw_half,
                                              bh_half,
                                              shift_factor);
    FP_ASSERT(sum <= (1 << (26)));
    window_error_quad_fp8[3] = (((sum << 4) / bw_half) << 4) / bh_half;

    if (is_chroma) {
        for (i = 0; i < 4; ++i) {
            FP_ASSERT(((int64_t)window_error_quad_fp8[i] * 5 + luma_window_error_quad_fp8[i]) <
                      ((int64_t)1 << 31));
            window_error_quad_fp8[i] = (window_error_quad_fp8[i] * 5 +
                                        luma_window_error_quad_fp8[i]) /
                6;
        }
    }

    //Calculation for every quarter
    for (subblock_idx = 0; subblock_idx < 4; subblock_idx++) {
        FP_ASSERT(((int64_t)window_error_quad_fp8[subblock_idx] * TF_WINDOW_BLOCK_BALANCE_WEIGHT +
                   block_error_fp8[subblock_idx]) < ((int64_t)1 << 31));
        uint32_t combined_error_fp8 = (window_error_quad_fp8[subblock_idx] *
                                           TF_WINDOW_BLOCK_BALANCE_WEIGHT +
                                       block_error_fp8[subblock_idx]) /
            (TF_WINDOW_BLOCK_BALANCE_WEIGHT + 1);

        uint32_t avg_err_fp10 = ((combined_error_fp8 >> 3) * (d_factor_fp8[subblock_idx] >> 3));
        // int32_t avg_err_fp12 = ((int64_t)(combined_error_fp16) * d_factor_fp16[subblock_idx])>>20;
        FP_ASSERT((((int64_t)combined_error_fp8 >> 3) * (d_factor_fp8[subblock_idx] >> 3)) <
                  ((int64_t)1 << 31));

        //double scaled_diff = AOMMIN(combined_error * d_factor[subblock_idx] / (FP2FLOAT(tf_decay_factor_fp16)), 7);
        uint32_t scaled_diff16 = AOMMIN(
            /*((16*avg_err)<<8)*/ (avg_err_fp10) / AOMMAX((tf_decay_factor_fp16 >> 10), 1), 7 * 16);
        //int adjusted_weight = (int)(expf((float)(-scaled_diff)) * TF_WEIGHT_SCALE);
        uint32_t adjusted_weight = (expf_tab_fp16[scaled_diff16] * TF_WEIGHT_SCALE) >> 16;

        int x_offset = (subblock_idx % 2) * block_width / 2;
        int y_offset = (subblock_idx / 2) * block_height / 2;

        for (i = 0; i < block_height / 2; i++) {
            for (j = 0; j < block_width / 2; j++) {
                const int k           = (i + y_offset) * y_pre_stride + j + x_offset;
                const int pixel_value = y_pre[k];
                y_count[k] += adjusted_weight;
                y_accum[k] += adjusted_weight * pixel_value;
            }
        }
    }
}

void svt_av1_apply_temporal_filter_planewise_medium_hbd_c(
    struct MeContext *context_ptr, const uint16_t *y_src, int y_src_stride, const uint16_t *y_pre,
    int y_pre_stride, const uint16_t *u_src, const uint16_t *v_src, int uv_src_stride,
    const uint16_t *u_pre, const uint16_t *v_pre, int uv_pre_stride, unsigned int block_width,
    unsigned int block_height, int ss_x, int ss_y, uint32_t *y_accum, uint16_t *y_count,
    uint32_t *u_accum, uint16_t *u_count, uint32_t *v_accum, uint16_t *v_count,
    uint32_t encoder_bit_depth) {
    uint32_t luma_window_error_quad_fp8[4];

    svt_av1_apply_temporal_filter_planewise_medium_hbd_partial_c(
        context_ptr,
        y_src,
        y_src_stride,
        y_pre,
        y_pre_stride,
        (unsigned int)block_width,
        (unsigned int)block_height,
        y_accum,
        y_count,
        context_ptr->tf_decay_factor_fp16[C_Y],
        luma_window_error_quad_fp8,
        0,
        encoder_bit_depth);

    if (context_ptr->tf_chroma) {
        svt_av1_apply_temporal_filter_planewise_medium_hbd_partial_c(
            context_ptr,
            u_src,
            uv_src_stride,
            u_pre,
            uv_pre_stride,
            (unsigned int)block_width >> ss_x,
            (unsigned int)block_height >> ss_y,
            u_accum,
            u_count,
            context_ptr->tf_decay_factor_fp16[C_U],
            luma_window_error_quad_fp8,
            1,
            encoder_bit_depth);

        svt_av1_apply_temporal_filter_planewise_medium_hbd_partial_c(
            context_ptr,
            v_src,
            uv_src_stride,
            v_pre,
            uv_pre_stride,
            (unsigned int)block_width >> ss_x,
            (unsigned int)block_height >> ss_y,
            v_accum,
            v_count,
            context_ptr->tf_decay_factor_fp16[C_V],
            luma_window_error_quad_fp8,
            1,
            encoder_bit_depth);
    }
}
/***************************************************************************************************
* Applies temporal filter for each block plane by plane. Passes the right inputs
* for 8 bit  and HBD path
* Inputs:
*   src : Pointer to the 8 bit frame to be filtered, which is used as
*                    reference to compute squared differece from the predictor.
*   src_16bit : Pointer to the hbd frame to be filtered, which is used as
*                    reference to compute squared differece from the predictor.
*   block_width: Width of the block.
*   block_height: Height of the block.
*   noise_levels: Pointer to the noise levels of the to-filter frame, estimated
*                 with each plane (in Y, U, V order).
*   pred: Pointers to the well-built 8 bit predictors.
*   pred_16bit: Pointers to the well-built hbd predictors.
*   accum: Pointer to the pixel-wise accumulator for filtering.
*   count: Pointer to the pixel-wise counter fot filtering.
* Returns:
*   Nothing will be returned. But the content to which `accum` and `pred`
*   point will be modified.
***************************************************************************************************/
static void apply_filtering_block_plane_wise(MeContext *context_ptr, int block_row, int block_col,
                                             EbByte *src, uint16_t **src_16bit, EbByte *pred,
                                             uint16_t **pred_16bit, uint32_t **accum,
                                             uint16_t **count, uint32_t *stride,
                                             uint32_t *stride_pred, int block_width,
                                             int block_height, uint32_t ss_x, uint32_t ss_y,
                                             uint32_t encoder_bit_depth) {
    int blk_h               = block_height;
    int blk_w               = block_width;
    int offset_src_buffer_Y = block_row * blk_h * stride[C_Y] + block_col * blk_w;
    int offset_src_buffer_U = block_row * (blk_h >> ss_y) * stride[C_U] +
        block_col * (blk_w >> ss_x);
    int offset_src_buffer_V = block_row * (blk_h >> ss_y) * stride[C_V] +
        block_col * (blk_w >> ss_x);

    int offset_block_buffer_Y = block_row * blk_h * stride_pred[C_Y] + block_col * blk_w;
    int offset_block_buffer_U = block_row * (blk_h >> ss_y) * stride_pred[C_U] +
        block_col * (blk_w >> ss_x);
    int offset_block_buffer_V = block_row * (blk_h >> ss_y) * stride_pred[C_V] +
        block_col * (blk_w >> ss_x);

    uint32_t *accum_ptr[COLOR_CHANNELS];
    uint16_t *count_ptr[COLOR_CHANNELS];

    accum_ptr[C_Y] = accum[C_Y] + offset_block_buffer_Y;
    accum_ptr[C_U] = accum[C_U] + offset_block_buffer_U;
    accum_ptr[C_V] = accum[C_V] + offset_block_buffer_V;

    count_ptr[C_Y] = count[C_Y] + offset_block_buffer_Y;
    count_ptr[C_U] = count[C_U] + offset_block_buffer_U;
    count_ptr[C_V] = count[C_V] + offset_block_buffer_V;

    if (encoder_bit_depth == 8) {
        uint8_t *src_ptr[COLOR_CHANNELS] = {
            src[C_Y] + offset_src_buffer_Y,
            src[C_U] + offset_src_buffer_U,
            src[C_V] + offset_src_buffer_V,
        };

        uint8_t *pred_ptr[COLOR_CHANNELS] = {
            pred[C_Y] + offset_block_buffer_Y,
            pred[C_U] + offset_block_buffer_U,
            pred[C_V] + offset_block_buffer_V,
        };

        if (context_ptr->tf_ctrls.use_fast_filter) {
            svt_av1_apply_temporal_filter_planewise_fast(context_ptr,
                                                         src_ptr[C_Y],
                                                         stride[C_Y],
                                                         pred_ptr[C_Y],
                                                         stride_pred[C_Y],
                                                         (unsigned int)block_width,
                                                         (unsigned int)block_height,
                                                         accum_ptr[C_Y],
                                                         count_ptr[C_Y]);
        } else {
            if (context_ptr->tf_ctrls.use_medium_filter) {
                svt_av1_apply_temporal_filter_planewise_medium(context_ptr,
                                                               src_ptr[C_Y],
                                                               stride[C_Y],
                                                               pred_ptr[C_Y],
                                                               stride_pred[C_Y],
                                                               src_ptr[C_U],
                                                               src_ptr[C_V],
                                                               stride[C_U],
                                                               pred_ptr[C_U],
                                                               pred_ptr[C_V],
                                                               stride_pred[C_U],
                                                               (unsigned int)block_width,
                                                               (unsigned int)block_height,
                                                               ss_x,
                                                               ss_y,
                                                               accum_ptr[C_Y],
                                                               count_ptr[C_Y],
                                                               accum_ptr[C_U],
                                                               count_ptr[C_U],
                                                               accum_ptr[C_V],
                                                               count_ptr[C_V]);
            } else
                svt_av1_apply_temporal_filter_planewise(context_ptr,
                                                        src_ptr[C_Y],
                                                        stride[C_Y],
                                                        pred_ptr[C_Y],
                                                        stride_pred[C_Y],
                                                        src_ptr[C_U],
                                                        src_ptr[C_V],
                                                        stride[C_U],
                                                        pred_ptr[C_U],
                                                        pred_ptr[C_V],
                                                        stride_pred[C_U],
                                                        (unsigned int)block_width,
                                                        (unsigned int)block_height,
                                                        ss_x,
                                                        ss_y,
                                                        accum_ptr[C_Y],
                                                        count_ptr[C_Y],
                                                        accum_ptr[C_U],
                                                        count_ptr[C_U],
                                                        accum_ptr[C_V],
                                                        count_ptr[C_V]);
        }
    } else {
        uint16_t *src_ptr_16bit[COLOR_CHANNELS] = {
            src_16bit[C_Y] + offset_src_buffer_Y,
            src_16bit[C_U] + offset_src_buffer_U,
            src_16bit[C_V] + offset_src_buffer_V,
        };

        uint16_t *pred_ptr_16bit[COLOR_CHANNELS] = {
            pred_16bit[C_Y] + offset_block_buffer_Y,
            pred_16bit[C_U] + offset_block_buffer_U,
            pred_16bit[C_V] + offset_block_buffer_V,
        };

        // Apply the temporal filtering strategy
        // TODO(any): avx2 version should also support high bit-depth.
        if (context_ptr->tf_ctrls.use_fast_filter) {
            svt_av1_apply_temporal_filter_planewise_fast_hbd(context_ptr,
                                                             src_ptr_16bit[C_Y],
                                                             stride[C_Y],
                                                             pred_ptr_16bit[C_Y],
                                                             stride_pred[C_Y],
                                                             (unsigned int)block_width,
                                                             (unsigned int)block_height,
                                                             accum_ptr[C_Y],
                                                             count_ptr[C_Y],
                                                             encoder_bit_depth);
        } else {
            if (context_ptr->tf_ctrls.use_medium_filter) {
                svt_av1_apply_temporal_filter_planewise_medium_hbd(context_ptr,
                                                                   src_ptr_16bit[C_Y],
                                                                   stride[C_Y],
                                                                   pred_ptr_16bit[C_Y],
                                                                   stride_pred[C_Y],
                                                                   src_ptr_16bit[C_U],
                                                                   src_ptr_16bit[C_V],
                                                                   stride[C_U],
                                                                   pred_ptr_16bit[C_U],
                                                                   pred_ptr_16bit[C_V],
                                                                   stride_pred[C_U],
                                                                   (unsigned int)block_width,
                                                                   (unsigned int)block_height,
                                                                   ss_x,
                                                                   ss_y,
                                                                   accum_ptr[C_Y],
                                                                   count_ptr[C_Y],
                                                                   accum_ptr[C_U],
                                                                   count_ptr[C_U],
                                                                   accum_ptr[C_V],
                                                                   count_ptr[C_V],
                                                                   encoder_bit_depth);
            } else
                svt_av1_apply_temporal_filter_planewise_hbd(context_ptr,
                                                            src_ptr_16bit[C_Y],
                                                            stride[C_Y],
                                                            pred_ptr_16bit[C_Y],
                                                            stride_pred[C_Y],
                                                            src_ptr_16bit[C_U],
                                                            src_ptr_16bit[C_V],
                                                            stride[C_U],
                                                            pred_ptr_16bit[C_U],
                                                            pred_ptr_16bit[C_V],
                                                            stride_pred[C_U],
                                                            (unsigned int)block_width,
                                                            (unsigned int)block_height,
                                                            ss_x,
                                                            ss_y,
                                                            accum_ptr[C_Y],
                                                            count_ptr[C_Y],
                                                            accum_ptr[C_U],
                                                            count_ptr[C_U],
                                                            accum_ptr[C_V],
                                                            count_ptr[C_V],
                                                            encoder_bit_depth);
        }
    }
}
uint32_t    get_mds_idx(uint32_t orgx, uint32_t orgy, uint32_t size, uint32_t use_128x128);
static void tf_16x16_sub_pel_search(PictureParentControlSet *pcs_ptr, MeContext *context_ptr,
                                    PictureParentControlSet *pcs_ref,
                                    EbPictureBufferDesc *pic_ptr_ref, EbByte *pred,
                                    uint16_t **pred_16bit, uint32_t *stride_pred, EbByte *src,
                                    uint16_t **src_16bit, uint32_t *stride_src,
                                    uint32_t sb_origin_x, uint32_t sb_origin_y, uint32_t ss_x,
                                    int encoder_bit_depth) {
    SequenceControlSet *scs_ptr = pcs_ptr->scs_ptr;
    InterpFilters interp_filters = av1_make_interp_filters(EIGHTTAP_REGULAR, EIGHTTAP_REGULAR);

    Bool is_highbd = (encoder_bit_depth == 8) ? (uint8_t)FALSE : (uint8_t)TRUE;

    BlkStruct   blk_ptr;
    MacroBlockD av1xd;
    blk_ptr.av1xd = &av1xd;
    MvUnit mv_unit;
    mv_unit.pred_direction = UNI_PRED_LIST_0;

    EbPictureBufferDesc reference_ptr;
    EbPictureBufferDesc prediction_ptr;

    UNUSED(ss_x);

    prediction_ptr.origin_x  = 0;
    prediction_ptr.origin_y  = 0;
    prediction_ptr.stride_y  = BW;
    prediction_ptr.stride_cb = (uint16_t)BW >> ss_x;
    prediction_ptr.stride_cr = (uint16_t)BW >> ss_x;

    if (!is_highbd) {
        assert(src[C_Y] != NULL);
        if (context_ptr->tf_chroma) {
            assert(src[C_U] != NULL);
            assert(src[C_V] != NULL);
        }
        prediction_ptr.buffer_y  = pred[C_Y];
        prediction_ptr.buffer_cb = pred[C_U];
        prediction_ptr.buffer_cr = pred[C_V];
    } else {
        assert(src_16bit[C_Y] != NULL);
        if (context_ptr->tf_chroma) {
            assert(src_16bit[C_U] != NULL);
            assert(src_16bit[C_V] != NULL);
        }
        prediction_ptr.buffer_y  = (uint8_t *)pred_16bit[C_Y];
        prediction_ptr.buffer_cb = (uint8_t *)pred_16bit[C_U];
        prediction_ptr.buffer_cr = (uint8_t *)pred_16bit[C_V];

        reference_ptr.buffer_y  = (uint8_t *)pcs_ref->altref_buffer_highbd[C_Y];
        reference_ptr.buffer_cb = (uint8_t *)pcs_ref->altref_buffer_highbd[C_U];
        reference_ptr.buffer_cr = (uint8_t *)pcs_ref->altref_buffer_highbd[C_V];
        reference_ptr.origin_x  = pic_ptr_ref->origin_x;
        reference_ptr.origin_y  = pic_ptr_ref->origin_y;
        reference_ptr.stride_y  = pic_ptr_ref->stride_y;
        reference_ptr.stride_cb = pic_ptr_ref->stride_cb;
        reference_ptr.stride_cr = pic_ptr_ref->stride_cr;
        reference_ptr.width     = pic_ptr_ref->width;
        reference_ptr.height    = pic_ptr_ref->height;
        reference_ptr.buffer_bit_inc_y  = NULL;
        reference_ptr.buffer_bit_inc_cb = NULL;
        reference_ptr.buffer_bit_inc_cr = NULL;

    }

    uint32_t bsize                             = 16;
    uint32_t idx_32x32                         = context_ptr->idx_32x32;
        for (uint32_t idx_16x16 = 0; idx_16x16 < 4; idx_16x16++) {
            uint32_t pu_index = idx_32x32_to_idx_16x16[idx_32x32][idx_16x16];

            uint32_t idx_y          = subblock_xy_16x16[pu_index][0];
            uint32_t idx_x          = subblock_xy_16x16[pu_index][1];
            uint16_t local_origin_x = idx_x * bsize;
            uint16_t local_origin_y = idx_y * bsize;
            uint16_t pu_origin_x    = sb_origin_x + local_origin_x;
            uint16_t pu_origin_y    = sb_origin_y + local_origin_y;
            int32_t mirow = pu_origin_y >> MI_SIZE_LOG2;
            int32_t micol = pu_origin_x >> MI_SIZE_LOG2;
            blk_ptr.mds_idx = get_mds_idx(local_origin_x,
                                          local_origin_y,
                                          bsize,
                                          pcs_ptr->scs_ptr->seq_header.sb_size == BLOCK_128X128);

            const int32_t bw                 = mi_size_wide[BLOCK_16X16];
            const int32_t bh                 = mi_size_high[BLOCK_16X16];
            blk_ptr.av1xd->mb_to_top_edge    = -(int32_t)((mirow * MI_SIZE) * 8);
            blk_ptr.av1xd->mb_to_bottom_edge = ((pcs_ptr->av1_cm->mi_rows - bw - mirow) * MI_SIZE) *
                8;
            blk_ptr.av1xd->mb_to_left_edge  = -(int32_t)((micol * MI_SIZE) * 8);
            blk_ptr.av1xd->mb_to_right_edge = ((pcs_ptr->av1_cm->mi_cols - bh - micol) * MI_SIZE) *
                8;

            uint32_t mv_index = tab16x16[pu_index];
            mv_unit.mv->x     = _MVXT(context_ptr->p_best_mv16x16[mv_index]);
            mv_unit.mv->y     = _MVYT(context_ptr->p_best_mv16x16[mv_index]);
            // AV1 MVs are always in 1/8th pel precision.
            mv_unit.mv->x = mv_unit.mv->x << 1;
            mv_unit.mv->y = mv_unit.mv->y << 1;

            context_ptr->tf_16x16_block_error[idx_32x32 * 4 + idx_16x16] = INT_MAX;
            signed short mv_x      = (_MVXT(context_ptr->p_best_mv16x16[mv_index])) << 1;
            signed short mv_y      = (_MVYT(context_ptr->p_best_mv16x16[mv_index])) << 1;
            signed short best_mv_x = mv_x;
            signed short best_mv_y = mv_y;

            if (!pcs_ptr->tf_ctrls.half_pel_mode && !pcs_ptr->tf_ctrls.quarter_pel_mode &&
                !pcs_ptr->tf_ctrls.eight_pel_mode) {
                mv_unit.mv->x = mv_x;
                mv_unit.mv->y = mv_y;

                av1_inter_prediction(scs_ptr,
                                     NULL, //pcs_ptr,
                                     (uint32_t)interp_filters,
                                     &blk_ptr,
                                     0, //ref_frame_type,
                                     &mv_unit,
                                     0, //use_intrabc,
                                     SIMPLE_TRANSLATION,
                                     0,
                                     0,
                                     1, //compound_idx not used
                                     NULL, // interinter_comp not used
                                     NULL,
                                     NULL,
                                     NULL,
                                     0,
                                     0,
                                     0,
                                     0,
                                     pu_origin_x,
                                     pu_origin_y,
                                     bsize,
                                     bsize,
                                     !is_highbd ? pic_ptr_ref : &reference_ptr,
                                     NULL, //ref_pic_list1,
                                     &prediction_ptr,
                                     local_origin_x,
                                     local_origin_y,
                                     PICTURE_BUFFER_DESC_LUMA_MASK,
                                     (uint8_t)encoder_bit_depth,
                                     0); // is_16bit_pipeline

                uint64_t distortion;
                if (!is_highbd) {
                    uint8_t *pred_y_ptr = pred[C_Y] + bsize * idx_y * stride_pred[C_Y] +
                        bsize * idx_x;
                    uint8_t *src_y_ptr = src[C_Y] + bsize * idx_y * stride_src[C_Y] + bsize * idx_x;

                    const AomVarianceFnPtr *fn_ptr = &mefn_ptr[BLOCK_16X16];

                    unsigned int sse;
                    distortion = fn_ptr->vf(
                        pred_y_ptr, stride_pred[C_Y], src_y_ptr, stride_src[C_Y], &sse);
                } else {
                    uint16_t *pred_y_ptr = pred_16bit[C_Y] + bsize * idx_y * stride_pred[C_Y] +
                        bsize * idx_x;
                    uint16_t *src_y_ptr = src_16bit[C_Y] + bsize * idx_y * stride_src[C_Y] +
                        bsize * idx_x;

                    unsigned int sse;
                    distortion = variance_highbd(
                        pred_y_ptr, stride_pred[C_Y], src_y_ptr, stride_src[C_Y], 16, 16, &sse);
                }
                if (distortion < context_ptr->tf_16x16_block_error[idx_32x32 * 4 + idx_16x16]) {
                    context_ptr->tf_16x16_block_error[idx_32x32 * 4 + idx_16x16] = distortion;
                    best_mv_x                                                    = mv_unit.mv->x;
                    best_mv_y                                                    = mv_unit.mv->y;
                }
            }
            // Perform 1/2 Pel MV Refinement
            for (signed short i = -4; i <= 4; i = i + 4) {
                for (signed short j = -4; j <= 4; j = j + 4) {
                    if (pcs_ptr->tf_ctrls.half_pel_mode == 2 && i != 0 && j != 0)
                        continue;
                    mv_unit.mv->x = mv_x + i;
                    mv_unit.mv->y = mv_y + j;

                    av1_inter_prediction(scs_ptr,
                                         NULL, //pcs_ptr,
                                         (uint32_t)interp_filters,
                                         &blk_ptr,
                                         0, //ref_frame_type,
                                         &mv_unit,
                                         0, //use_intrabc,
                                         SIMPLE_TRANSLATION,
                                         0,
                                         0,
                                         1, //compound_idx not used
                                         NULL, // interinter_comp not used
                                         NULL,
                                         NULL,
                                         NULL,
                                         0,
                                         0,
                                         0,
                                         0,
                                         pu_origin_x,
                                         pu_origin_y,
                                         bsize,
                                         bsize,
                                         !is_highbd ? pic_ptr_ref : &reference_ptr,
                                         NULL, //ref_pic_list1,
                                         &prediction_ptr,
                                         local_origin_x,
                                         local_origin_y,
                                         PICTURE_BUFFER_DESC_LUMA_MASK,
                                         (uint8_t)encoder_bit_depth,
                                         0); // is_16bit_pipeline

                    uint64_t distortion;
                    if (!is_highbd) {
                        uint8_t *pred_y_ptr = pred[C_Y] + bsize * idx_y * stride_pred[C_Y] +
                            bsize * idx_x;
                        uint8_t *src_y_ptr = src[C_Y] + bsize * idx_y * stride_src[C_Y] +
                            bsize * idx_x;

                        const AomVarianceFnPtr *fn_ptr = &mefn_ptr[BLOCK_16X16];

                        unsigned int sse;
                        distortion = fn_ptr->vf(
                            pred_y_ptr, stride_pred[C_Y], src_y_ptr, stride_src[C_Y], &sse);
                    } else {
                        uint16_t *pred_y_ptr = pred_16bit[C_Y] + bsize * idx_y * stride_pred[C_Y] +
                            bsize * idx_x;
                        uint16_t *src_y_ptr = src_16bit[C_Y] + bsize * idx_y * stride_src[C_Y] +
                            bsize * idx_x;

                        unsigned int sse;
                        distortion = variance_highbd(
                            pred_y_ptr, stride_pred[C_Y], src_y_ptr, stride_src[C_Y], 16, 16, &sse);
                    }
                    if (distortion < context_ptr->tf_16x16_block_error[idx_32x32 * 4 + idx_16x16]) {
                        context_ptr->tf_16x16_block_error[idx_32x32 * 4 + idx_16x16] = distortion;
                        best_mv_x = mv_unit.mv->x;
                        best_mv_y = mv_unit.mv->y;
                    }
                }
            }
            mv_x = best_mv_x;
            mv_y = best_mv_y;

            // Perform 1/4 Pel MV Refinement
            for (signed short i = -2; i <= 2; i = i + 2) {
                for (signed short j = -2; j <= 2; j = j + 2) {
                    if (pcs_ptr->tf_ctrls.quarter_pel_mode == 2 && i != 0 && j != 0)
                        continue;
                    mv_unit.mv->x = mv_x + i;
                    mv_unit.mv->y = mv_y + j;

                    av1_inter_prediction(scs_ptr,
                                         NULL, //pcs_ptr,
                                         (uint32_t)interp_filters,
                                         &blk_ptr,
                                         0, //ref_frame_type,
                                         &mv_unit,
                                         0, //use_intrabc,
                                         SIMPLE_TRANSLATION,
                                         0,
                                         0,
                                         1, //compound_idx not used
                                         NULL, // interinter_comp not used
                                         NULL,
                                         NULL,
                                         NULL,
                                         0,
                                         0,
                                         0,
                                         0,
                                         pu_origin_x,
                                         pu_origin_y,
                                         bsize,
                                         bsize,
                                         !is_highbd ? pic_ptr_ref : &reference_ptr,
                                         NULL, //ref_pic_list1,
                                         &prediction_ptr,
                                         local_origin_x,
                                         local_origin_y,
                                         PICTURE_BUFFER_DESC_LUMA_MASK,
                                         (uint8_t)encoder_bit_depth,
                                         0); // is_16bit_pipeline

                    uint64_t distortion;
                    if (!is_highbd) {
                        uint8_t *pred_y_ptr = pred[C_Y] + bsize * idx_y * stride_pred[C_Y] +
                            bsize * idx_x;
                        uint8_t *src_y_ptr = src[C_Y] + bsize * idx_y * stride_src[C_Y] +
                            bsize * idx_x;

                        const AomVarianceFnPtr *fn_ptr = &mefn_ptr[BLOCK_16X16];

                        unsigned int sse;
                        distortion = fn_ptr->vf(
                            pred_y_ptr, stride_pred[C_Y], src_y_ptr, stride_src[C_Y], &sse);
                    } else {
                        uint16_t *pred_y_ptr = pred_16bit[C_Y] + bsize * idx_y * stride_pred[C_Y] +
                            bsize * idx_x;
                        uint16_t *src_y_ptr = src_16bit[C_Y] + bsize * idx_y * stride_src[C_Y] +
                            bsize * idx_x;
                        ;

                        unsigned int sse;
                        distortion = variance_highbd(
                            pred_y_ptr, stride_pred[C_Y], src_y_ptr, stride_src[C_Y], 16, 16, &sse);
                    }
                    if (distortion < context_ptr->tf_16x16_block_error[idx_32x32 * 4 + idx_16x16]) {
                        context_ptr->tf_16x16_block_error[idx_32x32 * 4 + idx_16x16] = distortion;
                        best_mv_x = mv_unit.mv->x;
                        best_mv_y = mv_unit.mv->y;
                    }
                }
            }

            mv_x = best_mv_x;
            mv_y = best_mv_y;
            // Perform 1/8 Pel MV Refinement
            if (pcs_ptr->tf_ctrls.eight_pel_mode)
                for (signed short i = -1; i <= 1; i++) {
                    for (signed short j = -1; j <= 1; j++) {
                        if (pcs_ptr->tf_ctrls.eight_pel_mode == 2 && i != 0 && j != 0)
                            continue;
                        mv_unit.mv->x = mv_x + i;
                        mv_unit.mv->y = mv_y + j;

                        av1_inter_prediction(scs_ptr,
                                             NULL, //pcs_ptr,
                                             (uint32_t)interp_filters,
                                             &blk_ptr,
                                             0, //ref_frame_type,
                                             &mv_unit,
                                             0, //use_intrabc,
                                             SIMPLE_TRANSLATION,
                                             0,
                                             0,
                                             1, //compound_idx not used
                                             NULL, // interinter_comp not used
                                             NULL,
                                             NULL,
                                             NULL,
                                             0,
                                             0,
                                             0,
                                             0,
                                             pu_origin_x,
                                             pu_origin_y,
                                             bsize,
                                             bsize,
                                             !is_highbd ? pic_ptr_ref : &reference_ptr,
                                             NULL, //ref_pic_list1,
                                             &prediction_ptr,
                                             local_origin_x,
                                             local_origin_y,
                                             PICTURE_BUFFER_DESC_LUMA_MASK,
                                             (uint8_t)encoder_bit_depth,
                                             0); // is_16bit_pipeline

                        uint64_t distortion;
                        if (!is_highbd) {
                            uint8_t *pred_y_ptr = pred[C_Y] + bsize * idx_y * stride_pred[C_Y] +
                                bsize * idx_x;
                            uint8_t *src_y_ptr = src[C_Y] + bsize * idx_y * stride_src[C_Y] +
                                bsize * idx_x;

                            const AomVarianceFnPtr *fn_ptr = &mefn_ptr[BLOCK_16X16];

                            unsigned int sse;
                            distortion = fn_ptr->vf(
                                pred_y_ptr, stride_pred[C_Y], src_y_ptr, stride_src[C_Y], &sse);
                        } else {
                            uint16_t *pred_y_ptr = pred_16bit[C_Y] +
                                bsize * idx_y * stride_pred[C_Y] + bsize * idx_x;
                            uint16_t *src_y_ptr = src_16bit[C_Y] + bsize * idx_y * stride_src[C_Y] +
                                bsize * idx_x;

                            unsigned int sse;
                            distortion = variance_highbd(pred_y_ptr,
                                                         stride_pred[C_Y],
                                                         src_y_ptr,
                                                         stride_src[C_Y],
                                                         16,
                                                         16,
                                                         &sse);
                        }
                        if (distortion <
                            context_ptr->tf_16x16_block_error[idx_32x32 * 4 + idx_16x16]) {
                            context_ptr->tf_16x16_block_error[idx_32x32 * 4 + idx_16x16] =
                                distortion;
                            best_mv_x = mv_unit.mv->x;
                            best_mv_y = mv_unit.mv->y;
                        }
                    }
                }
            context_ptr->tf_16x16_mv_x[idx_32x32 * 4 + idx_16x16] = best_mv_x;
            context_ptr->tf_16x16_mv_y[idx_32x32 * 4 + idx_16x16] = best_mv_y;
        }
}

uint64_t svt_check_position_64x64(TF_SUBPEL_SEARCH_PARAMS  tf_sp_param,
                                  PictureParentControlSet *pcs_ptr, MeContext *context_ptr,
                                  BlkStruct *blk_ptr,
                                  MvUnit mv_unit, EbPictureBufferDesc *pic_ptr_ref,
                                  EbPictureBufferDesc prediction_ptr, EbByte *pred,
                                  uint16_t **pred_16bit, uint32_t *stride_pred, EbByte *src,
                                  uint16_t **src_16bit, uint32_t *stride_src,
                                  signed short *best_mv_x, signed short *best_mv_y) {
    if (tf_sp_param.subpel_pel_mode == 2 && tf_sp_param.xd != 0 && tf_sp_param.yd != 0)
        return UINT_MAX;
    // if previously checked position is good enough then quit
    if (context_ptr->tf_subpel_early_exit &&
        context_ptr->tf_64x64_block_error < (((tf_sp_param.bsize * tf_sp_param.bsize) << 2) << tf_sp_param.is_highbd))
        return UINT_MAX;
    SequenceControlSet *scs_ptr = pcs_ptr->scs_ptr;
    mv_unit.mv->x = tf_sp_param.mv_x + tf_sp_param.xd;
    mv_unit.mv->y = tf_sp_param.mv_y + tf_sp_param.yd;

    av1_simple_luma_unipred(
        scs_ptr,
        scs_ptr->sf_identity,
        tf_sp_param.interp_filters,
        blk_ptr,
        0, //ref_frame_type,
        &mv_unit,
        tf_sp_param.pu_origin_x,
        tf_sp_param.pu_origin_y,
        tf_sp_param.bsize,
        tf_sp_param.bsize,
        pic_ptr_ref,
        &prediction_ptr,
        tf_sp_param.local_origin_x,
        tf_sp_param.local_origin_y,
        (uint8_t)tf_sp_param.encoder_bit_depth,
        tf_sp_param.xd == 0 && tf_sp_param.yd == 0 ? tf_sp_param.subsampling_shift : 0);

    uint64_t distortion;
    if (!tf_sp_param.is_highbd) {
        uint8_t *pred_y_ptr = pred[C_Y] + tf_sp_param.bsize * tf_sp_param.idx_y * stride_pred[C_Y] +
            tf_sp_param.bsize * tf_sp_param.idx_x;
        uint8_t *src_y_ptr = src[C_Y] + tf_sp_param.bsize * tf_sp_param.idx_y * stride_src[C_Y] +
            tf_sp_param.bsize * tf_sp_param.idx_x;
        const AomVarianceFnPtr *fn_ptr = tf_sp_param.subsampling_shift ? &mefn_ptr[BLOCK_64X32]
                                                                       : &mefn_ptr[BLOCK_64X64];
        unsigned int            sse;
        distortion = fn_ptr->vf(pred_y_ptr,
                                stride_pred[C_Y] << tf_sp_param.subsampling_shift,
                                src_y_ptr,
                                stride_src[C_Y] << tf_sp_param.subsampling_shift,
                                &sse)
            << tf_sp_param.subsampling_shift;
    } else {
        uint16_t *pred_y_ptr = pred_16bit[C_Y] +
            tf_sp_param.bsize * tf_sp_param.idx_y * stride_pred[C_Y] +
            tf_sp_param.bsize * tf_sp_param.idx_x;
        uint16_t *src_y_ptr = src_16bit[C_Y] +
            tf_sp_param.bsize * tf_sp_param.idx_y * stride_src[C_Y] +
            tf_sp_param.bsize * tf_sp_param.idx_x;
        const AomVarianceFnPtr *fn_ptr = tf_sp_param.subsampling_shift ? &mefn_ptr[BLOCK_64X32]
                                                                       : &mefn_ptr[BLOCK_64X64];

        unsigned int sse;

        distortion = fn_ptr->vf_hbd_10(CONVERT_TO_BYTEPTR(pred_y_ptr),
                                       stride_pred[C_Y] << tf_sp_param.subsampling_shift,
                                       CONVERT_TO_BYTEPTR(src_y_ptr),
                                       stride_src[C_Y] << tf_sp_param.subsampling_shift,
                                       &sse)
            << tf_sp_param.subsampling_shift;
    }
    if (distortion < context_ptr->tf_64x64_block_error) {
        context_ptr->tf_64x64_block_error = distortion;
        *best_mv_x                        = mv_unit.mv->x;
        *best_mv_y                        = mv_unit.mv->y;
    }
    return distortion;
}
uint64_t svt_check_position(TF_SUBPEL_SEARCH_PARAMS tf_sp_param, PictureParentControlSet *pcs_ptr,
                            MeContext *context_ptr,
                            BlkStruct *blk_ptr,
                            MvUnit mv_unit, EbPictureBufferDesc *pic_ptr_ref,
                            EbPictureBufferDesc prediction_ptr, EbByte *pred, uint16_t **pred_16bit,
                            uint32_t *stride_pred, EbByte *src, uint16_t **src_16bit,
                            uint32_t *stride_src, signed short *best_mv_x,
                            signed short *best_mv_y) {
    if (tf_sp_param.subpel_pel_mode == 2 && tf_sp_param.xd != 0 && tf_sp_param.yd != 0)
        return UINT_MAX;
    // if previously checked position is good enough then quit
    if (context_ptr->tf_subpel_early_exit &&
        context_ptr->tf_32x32_block_error[context_ptr->idx_32x32] < (((tf_sp_param.bsize * tf_sp_param.bsize) >> 7) << tf_sp_param.is_highbd))
        return UINT_MAX;
    SequenceControlSet *scs_ptr = pcs_ptr->scs_ptr;
    mv_unit.mv->x               = tf_sp_param.mv_x + tf_sp_param.xd;
    mv_unit.mv->y               = tf_sp_param.mv_y + tf_sp_param.yd;

    av1_simple_luma_unipred(
        scs_ptr,
        scs_ptr->sf_identity,
        tf_sp_param.interp_filters,
        blk_ptr,
        0, //ref_frame_type,
        &mv_unit,
        tf_sp_param.pu_origin_x,
        tf_sp_param.pu_origin_y,
        tf_sp_param.bsize,
        tf_sp_param.bsize,
        pic_ptr_ref,
        &prediction_ptr,
        tf_sp_param.local_origin_x,
        tf_sp_param.local_origin_y,
        (uint8_t)tf_sp_param.encoder_bit_depth,
        tf_sp_param.xd == 0 && tf_sp_param.yd == 0 ? tf_sp_param.subsampling_shift : 0);

    uint64_t distortion;
    if (!tf_sp_param.is_highbd) {
        uint8_t *pred_y_ptr = pred[C_Y] + tf_sp_param.bsize * tf_sp_param.idx_y * stride_pred[C_Y] +
            tf_sp_param.bsize * tf_sp_param.idx_x;
        uint8_t *src_y_ptr = src[C_Y] + tf_sp_param.bsize * tf_sp_param.idx_y * stride_src[C_Y] +
            tf_sp_param.bsize * tf_sp_param.idx_x;
        const AomVarianceFnPtr *fn_ptr = tf_sp_param.subsampling_shift ? &mefn_ptr[BLOCK_32X16]
                                                                       : &mefn_ptr[BLOCK_32X32];
        unsigned int            sse;
        distortion = fn_ptr->vf(pred_y_ptr,
                                stride_pred[C_Y] << tf_sp_param.subsampling_shift,
                                src_y_ptr,
                                stride_src[C_Y] << tf_sp_param.subsampling_shift,
                                &sse)
            << tf_sp_param.subsampling_shift;
    } else {
        uint16_t *pred_y_ptr = pred_16bit[C_Y] +
            tf_sp_param.bsize * tf_sp_param.idx_y * stride_pred[C_Y] +
            tf_sp_param.bsize * tf_sp_param.idx_x;
        uint16_t *src_y_ptr = src_16bit[C_Y] +
            tf_sp_param.bsize * tf_sp_param.idx_y * stride_src[C_Y] +
            tf_sp_param.bsize * tf_sp_param.idx_x;
        const AomVarianceFnPtr *fn_ptr = tf_sp_param.subsampling_shift ? &mefn_ptr[BLOCK_32X16]
                                                                       : &mefn_ptr[BLOCK_32X32];

        unsigned int sse;

        distortion = fn_ptr->vf_hbd_10(CONVERT_TO_BYTEPTR(pred_y_ptr),
                                       stride_pred[C_Y] << tf_sp_param.subsampling_shift,
                                       CONVERT_TO_BYTEPTR(src_y_ptr),
                                       stride_src[C_Y] << tf_sp_param.subsampling_shift,
                                       &sse)
            << tf_sp_param.subsampling_shift;
    }
    if (distortion < context_ptr->tf_32x32_block_error[context_ptr->idx_32x32]) {
        context_ptr->tf_32x32_block_error[context_ptr->idx_32x32] = distortion;
        *best_mv_x                                                = mv_unit.mv->x;
        *best_mv_y                                                = mv_unit.mv->y;
    }
    return distortion;
}
static void tf_64x64_sub_pel_search(PictureParentControlSet *pcs_ptr, MeContext *context_ptr,
                                    PictureParentControlSet *pcs_ref,
                                    EbPictureBufferDesc *pic_ptr_ref, EbByte *pred,
                                    uint16_t **pred_16bit, uint32_t *stride_pred, EbByte *src,
                                    uint16_t **src_16bit, uint32_t *stride_src,
                                    uint32_t sb_origin_x, uint32_t sb_origin_y, uint32_t ss_x,
                                    int encoder_bit_depth) {
    InterpFilters interp_filters;
    if (context_ptr->tf_ctrls.use_2tap)
        interp_filters = av1_make_interp_filters(BILINEAR, BILINEAR);
    else
        interp_filters = av1_make_interp_filters(EIGHTTAP_REGULAR, EIGHTTAP_REGULAR);

    Bool is_highbd = (encoder_bit_depth == 8) ? (uint8_t)FALSE : (uint8_t)TRUE;
    BlkStruct blk_struct;
    MacroBlockD av1xd;
    blk_struct.av1xd = &av1xd;
    MvUnit mv_unit;
    mv_unit.pred_direction = UNI_PRED_LIST_0;
    EbPictureBufferDesc reference_ptr;
    EbPictureBufferDesc prediction_ptr;
    UNUSED(ss_x);
    prediction_ptr.origin_x  = 0;
    prediction_ptr.origin_y  = 0;
    prediction_ptr.stride_y  = BW;
    prediction_ptr.stride_cb = (uint16_t)BW >> ss_x;
    prediction_ptr.stride_cr = (uint16_t)BW >> ss_x;
    if (!is_highbd) {
        assert(src[C_Y] != NULL);
        if (context_ptr->tf_chroma) {
            assert(src[C_U] != NULL);
            assert(src[C_V] != NULL);
        }
        prediction_ptr.buffer_y  = pred[C_Y];
        prediction_ptr.buffer_cb = pred[C_U];
        prediction_ptr.buffer_cr = pred[C_V];
    } else {
        assert(src_16bit[C_Y] != NULL);
        if (context_ptr->tf_chroma) {
            assert(src_16bit[C_U] != NULL);
            assert(src_16bit[C_V] != NULL);
        }
        prediction_ptr.buffer_y  = (uint8_t *)pred_16bit[C_Y];
        prediction_ptr.buffer_cb = (uint8_t *)pred_16bit[C_U];
        prediction_ptr.buffer_cr = (uint8_t *)pred_16bit[C_V];
        reference_ptr.buffer_y   = (uint8_t *)pcs_ref->altref_buffer_highbd[C_Y];
        reference_ptr.buffer_cb  = (uint8_t *)pcs_ref->altref_buffer_highbd[C_U];
        reference_ptr.buffer_cr  = (uint8_t *)pcs_ref->altref_buffer_highbd[C_V];
        reference_ptr.origin_x   = pic_ptr_ref->origin_x;
        reference_ptr.origin_y   = pic_ptr_ref->origin_y;
        reference_ptr.stride_y   = pic_ptr_ref->stride_y;
        reference_ptr.stride_cb  = pic_ptr_ref->stride_cb;
        reference_ptr.stride_cr  = pic_ptr_ref->stride_cr;
        reference_ptr.width      = pic_ptr_ref->width;
        reference_ptr.height     = pic_ptr_ref->height;
        reference_ptr.buffer_bit_inc_y  = NULL;
        reference_ptr.buffer_bit_inc_cb = NULL;
        reference_ptr.buffer_bit_inc_cr = NULL;

    }
    uint32_t bsize = 64;

    uint16_t local_origin_x = 0;
    uint16_t local_origin_y = 0;
    uint16_t pu_origin_x    = sb_origin_x + local_origin_x;
    uint16_t pu_origin_y    = sb_origin_y + local_origin_y;
    int32_t mirow = pu_origin_y >> MI_SIZE_LOG2;
    int32_t micol = pu_origin_x >> MI_SIZE_LOG2;

    blk_struct.mds_idx                  = 0;
    const int32_t bw                    = mi_size_wide[BLOCK_64X64];
    const int32_t bh                    = mi_size_high[BLOCK_64X64];
    blk_struct.av1xd->mb_to_top_edge    = -(int32_t)((mirow * MI_SIZE) * 8);
    blk_struct.av1xd->mb_to_bottom_edge = ((pcs_ptr->av1_cm->mi_rows - bw - mirow) * MI_SIZE) * 8;
    blk_struct.av1xd->mb_to_left_edge   = -(int32_t)((micol * MI_SIZE) * 8);
    blk_struct.av1xd->mb_to_right_edge  = ((pcs_ptr->av1_cm->mi_cols - bh - micol) * MI_SIZE) * 8;
    context_ptr->tf_64x64_block_error   = INT_MAX;

    signed short mv_x = mv_unit.mv->x = (context_ptr->tf_use_pred_64x64_only_th == (uint8_t)~0)
        ? context_ptr->search_results[0][0].hme_sc_x << 3 : (_MVXT(context_ptr->p_best_mv64x64[0])) << 1;

    signed short mv_y = mv_unit.mv->y = (context_ptr->tf_use_pred_64x64_only_th == (uint8_t)~0)
        ? context_ptr->search_results[0][0].hme_sc_y << 3 : (_MVYT(context_ptr->p_best_mv64x64[0])) << 1;

    BlkStruct *blk_ptr = &blk_struct;
    signed short            best_mv_x = mv_x;
    signed short            best_mv_y = mv_y;
    TF_SUBPEL_SEARCH_PARAMS tf_sp_param;
    tf_sp_param.subsampling_shift = pcs_ptr->tf_ctrls.sub_sampling_shift;
    uint64_t center_err       = UINT_MAX;
    uint64_t min_error        = UINT_MAX;
    uint64_t left_bottom_err  = UINT_MAX;
    uint64_t left_top_err     = UINT_MAX;
    uint64_t top_right_err    = UINT_MAX;
    uint64_t bottom_right_err = UINT_MAX;

    uint8_t best_direction = 0; //horz
    // No refinement.
    if (!pcs_ptr->tf_ctrls.half_pel_mode && !pcs_ptr->tf_ctrls.quarter_pel_mode &&
        !pcs_ptr->tf_ctrls.eight_pel_mode) {
        tf_sp_param.subpel_pel_mode   = pcs_ptr->tf_ctrls.half_pel_mode;
        tf_sp_param.mv_x              = mv_x;
        tf_sp_param.mv_y              = mv_y;
        tf_sp_param.interp_filters    = (uint32_t)interp_filters;
        tf_sp_param.pu_origin_x       = pu_origin_x;
        tf_sp_param.pu_origin_y       = pu_origin_y;
        tf_sp_param.local_origin_x    = local_origin_x;
        tf_sp_param.local_origin_y    = local_origin_y;
        tf_sp_param.bsize             = bsize;
        tf_sp_param.is_highbd         = is_highbd;
        tf_sp_param.encoder_bit_depth = encoder_bit_depth;
        tf_sp_param.idx_x             = 0;
        tf_sp_param.idx_y             = 0;
        tf_sp_param.xd                = 0;
        tf_sp_param.yd                = 0;
        center_err                    = svt_check_position_64x64(tf_sp_param,
                                              pcs_ptr,
                                              context_ptr,
                                              blk_ptr,
                                              mv_unit,
                                              !is_highbd ? pic_ptr_ref : &reference_ptr,
                                              prediction_ptr,
                                              pred,
                                              pred_16bit,
                                              stride_pred,
                                              src,
                                              src_16bit,
                                              stride_src,
                                              &best_mv_x,
                                              &best_mv_y);
    }

    // Perform 1/2 Pel MV Refinement
    if (pcs_ptr->tf_ctrls.half_pel_mode) {
        tf_sp_param.subpel_pel_mode   = pcs_ptr->tf_ctrls.half_pel_mode;
        tf_sp_param.mv_x              = mv_x;
        tf_sp_param.mv_y              = mv_y;
        tf_sp_param.interp_filters    = (uint32_t)interp_filters;
        tf_sp_param.pu_origin_x       = pu_origin_x;
        tf_sp_param.pu_origin_y       = pu_origin_y;
        tf_sp_param.local_origin_x    = local_origin_x;
        tf_sp_param.local_origin_y    = local_origin_y;
        tf_sp_param.bsize             = bsize;
        tf_sp_param.is_highbd         = is_highbd;
        tf_sp_param.encoder_bit_depth = encoder_bit_depth;
        tf_sp_param.idx_x             = 0;
        tf_sp_param.idx_y             = 0;
        if (pcs_ptr->tf_ctrls.half_pel_mode == 1) {
            for (signed short i = -4; i <= 4; i = i + 4) {
                for (signed short j = -4; j <= 4; j = j + 4) {
                    tf_sp_param.xd   = i;
                    tf_sp_param.yd   = j;
                    uint64_t cur_err = svt_check_position_64x64(
                        tf_sp_param,
                        pcs_ptr,
                        context_ptr,
                        blk_ptr,
                        mv_unit,
                        !is_highbd ? pic_ptr_ref : &reference_ptr,
                        prediction_ptr,
                        pred,
                        pred_16bit,
                        stride_pred,
                        src,
                        src_16bit,
                        stride_src,
                        &best_mv_x,
                        &best_mv_y);
                    if (cur_err < min_error)
                        min_error = cur_err;
                }
            }
        } else {
            //check center
            tf_sp_param.xd = 0;
            tf_sp_param.yd = 0;
            center_err     = svt_check_position_64x64(tf_sp_param,
                                                  pcs_ptr,
                                                  context_ptr,
                                                  blk_ptr,
                                                  mv_unit,
                                                  !is_highbd ? pic_ptr_ref : &reference_ptr,
                                                  prediction_ptr,
                                                  pred,
                                                  pred_16bit,
                                                  stride_pred,
                                                  src,
                                                  src_16bit,
                                                  stride_src,
                                                  &best_mv_x,
                                                  &best_mv_y);

            //check left
            tf_sp_param.xd    = -4;
            tf_sp_param.yd    = 0;
            uint64_t left_err = svt_check_position_64x64(tf_sp_param,
                                                         pcs_ptr,
                                                         context_ptr,
                                                         blk_ptr,
                                                         mv_unit,
                                                         !is_highbd ? pic_ptr_ref : &reference_ptr,
                                                         prediction_ptr,
                                                         pred,
                                                         pred_16bit,
                                                         stride_pred,
                                                         src,
                                                         src_16bit,
                                                         stride_src,
                                                         &best_mv_x,
                                                         &best_mv_y);

            //check top
            tf_sp_param.xd   = 0;
            tf_sp_param.yd   = -4;
            uint64_t top_err = svt_check_position_64x64(tf_sp_param,
                                                        pcs_ptr,
                                                        context_ptr,
                                                        blk_ptr,
                                                        mv_unit,
                                                        !is_highbd ? pic_ptr_ref : &reference_ptr,
                                                        prediction_ptr,
                                                        pred,
                                                        pred_16bit,
                                                        stride_pred,
                                                        src,
                                                        src_16bit,
                                                        stride_src,
                                                        &best_mv_x,
                                                        &best_mv_y);

            //check right
            tf_sp_param.xd     = 4;
            tf_sp_param.yd     = 0;
            uint64_t right_err = svt_check_position_64x64(tf_sp_param,
                                                          pcs_ptr,
                                                          context_ptr,
                                                          blk_ptr,
                                                          mv_unit,
                                                          !is_highbd ? pic_ptr_ref : &reference_ptr,
                                                          prediction_ptr,
                                                          pred,
                                                          pred_16bit,
                                                          stride_pred,
                                                          src,
                                                          src_16bit,
                                                          stride_src,
                                                          &best_mv_x,
                                                          &best_mv_y);

            //check bottom
            tf_sp_param.xd      = 0;
            tf_sp_param.yd      = 4;
            uint64_t bottom_err = svt_check_position_64x64(
                tf_sp_param,
                pcs_ptr,
                context_ptr,
                blk_ptr,
                mv_unit,
                !is_highbd ? pic_ptr_ref : &reference_ptr,
                prediction_ptr,
                pred,
                pred_16bit,
                stride_pred,
                src,
                src_16bit,
                stride_src,
                &best_mv_x,
                &best_mv_y);

            min_error = MIN(left_err, MIN(top_err, MIN(right_err, bottom_err)));

            if (min_error == top_err || min_error == bottom_err)
                best_direction = 1; //vert

            if (min_error == left_err) {
                //check left_top
                tf_sp_param.xd = -4;
                tf_sp_param.yd = -4;
                left_top_err   = svt_check_position_64x64(tf_sp_param,
                                                        pcs_ptr,
                                                        context_ptr,
                                                        blk_ptr,
                                                        mv_unit,
                                                        !is_highbd ? pic_ptr_ref : &reference_ptr,
                                                        prediction_ptr,
                                                        pred,
                                                        pred_16bit,
                                                        stride_pred,
                                                        src,
                                                        src_16bit,
                                                        stride_src,
                                                        &best_mv_x,
                                                        &best_mv_y);
                //check left_bottom
                tf_sp_param.xd  = -4;
                tf_sp_param.yd  = 4;
                left_bottom_err = svt_check_position_64x64(
                    tf_sp_param,
                    pcs_ptr,
                    context_ptr,
                    blk_ptr,
                    mv_unit,
                    !is_highbd ? pic_ptr_ref : &reference_ptr,
                    prediction_ptr,
                    pred,
                    pred_16bit,
                    stride_pred,
                    src,
                    src_16bit,
                    stride_src,
                    &best_mv_x,
                    &best_mv_y);
            } else if (min_error == top_err) {
                //check left_top
                tf_sp_param.xd = -4;
                tf_sp_param.yd = -4;
                left_top_err   = svt_check_position_64x64(tf_sp_param,
                                                        pcs_ptr,
                                                        context_ptr,
                                                        blk_ptr,
                                                        mv_unit,
                                                        !is_highbd ? pic_ptr_ref : &reference_ptr,
                                                        prediction_ptr,
                                                        pred,
                                                        pred_16bit,
                                                        stride_pred,
                                                        src,
                                                        src_16bit,
                                                        stride_src,
                                                        &best_mv_x,
                                                        &best_mv_y);
                //check top_right
                tf_sp_param.xd = 4;
                tf_sp_param.yd = -4;
                top_right_err  = svt_check_position_64x64(tf_sp_param,
                                                         pcs_ptr,
                                                         context_ptr,
                                                         blk_ptr,
                                                         mv_unit,
                                                         !is_highbd ? pic_ptr_ref : &reference_ptr,
                                                         prediction_ptr,
                                                         pred,
                                                         pred_16bit,
                                                         stride_pred,
                                                         src,
                                                         src_16bit,
                                                         stride_src,
                                                         &best_mv_x,
                                                         &best_mv_y);
            } else if (min_error == right_err) {
                //check top_right
                tf_sp_param.xd = 4;
                tf_sp_param.yd = -4;
                top_right_err  = svt_check_position_64x64(tf_sp_param,
                                                         pcs_ptr,
                                                         context_ptr,
                                                         blk_ptr,
                                                         mv_unit,
                                                         !is_highbd ? pic_ptr_ref : &reference_ptr,
                                                         prediction_ptr,
                                                         pred,
                                                         pred_16bit,
                                                         stride_pred,
                                                         src,
                                                         src_16bit,
                                                         stride_src,
                                                         &best_mv_x,
                                                         &best_mv_y);
                //check bottom_right
                tf_sp_param.xd   = 4;
                tf_sp_param.yd   = 4;
                bottom_right_err = svt_check_position_64x64(
                    tf_sp_param,
                    pcs_ptr,
                    context_ptr,
                    blk_ptr,
                    mv_unit,
                    !is_highbd ? pic_ptr_ref : &reference_ptr,
                    prediction_ptr,
                    pred,
                    pred_16bit,
                    stride_pred,
                    src,
                    src_16bit,
                    stride_src,
                    &best_mv_x,
                    &best_mv_y);
            } else if (min_error == bottom_err) {
                //check bottom_left
                tf_sp_param.xd  = -4;
                tf_sp_param.yd  = 4;
                left_bottom_err = svt_check_position_64x64(
                    tf_sp_param,
                    pcs_ptr,
                    context_ptr,
                    blk_ptr,
                    mv_unit,
                    !is_highbd ? pic_ptr_ref : &reference_ptr,
                    prediction_ptr,
                    pred,
                    pred_16bit,
                    stride_pred,
                    src,
                    src_16bit,
                    stride_src,
                    &best_mv_x,
                    &best_mv_y);
                //check bottom_right
                tf_sp_param.xd   = 4;
                tf_sp_param.yd   = 4;
                bottom_right_err = svt_check_position_64x64(
                    tf_sp_param,
                    pcs_ptr,
                    context_ptr,
                    blk_ptr,
                    mv_unit,
                    !is_highbd ? pic_ptr_ref : &reference_ptr,
                    prediction_ptr,
                    pred,
                    pred_16bit,
                    stride_pred,
                    src,
                    src_16bit,
                    stride_src,
                    &best_mv_x,
                    &best_mv_y);
            }
        }
        min_error = MIN(
            min_error,
            MIN(left_top_err, MIN(left_bottom_err, MIN(top_right_err, bottom_right_err))));
    }

    mv_x = best_mv_x;
    mv_y = best_mv_y;
    // Perform 1/4 Pel MV Refinement
    if (pcs_ptr->tf_ctrls.quarter_pel_mode) {
        tf_sp_param.subpel_pel_mode   = pcs_ptr->tf_ctrls.quarter_pel_mode;
        tf_sp_param.mv_x              = mv_x;
        tf_sp_param.mv_y              = mv_y;
        tf_sp_param.interp_filters    = (uint32_t)interp_filters;
        tf_sp_param.pu_origin_x       = pu_origin_x;
        tf_sp_param.pu_origin_y       = pu_origin_y;
        tf_sp_param.local_origin_x    = local_origin_x;
        tf_sp_param.local_origin_y    = local_origin_y;
        tf_sp_param.bsize             = bsize;
        tf_sp_param.is_highbd         = is_highbd;
        tf_sp_param.encoder_bit_depth = encoder_bit_depth;
        tf_sp_param.idx_x             = 0;
        tf_sp_param.idx_y             = 0;
        if (pcs_ptr->tf_ctrls.quarter_pel_mode == 1) {
            for (signed short i = -2; i <= 2; i = i + 2) {
                for (signed short j = -2; j <= 2; j = j + 2) {
                    tf_sp_param.xd = i;
                    tf_sp_param.yd = j;
                    center_err     = svt_check_position_64x64(tf_sp_param,
                                                          pcs_ptr,
                                                          context_ptr,
                                                          blk_ptr,
                                                          mv_unit,
                                                          !is_highbd ? pic_ptr_ref : &reference_ptr,
                                                          prediction_ptr,
                                                          pred,
                                                          pred_16bit,
                                                          stride_pred,
                                                          src,
                                                          src_16bit,
                                                          stride_src,
                                                          &best_mv_x,
                                                          &best_mv_y);
                }
            }
        } else if (min_error < center_err) {
            //check center
            tf_sp_param.xd = 0;
            tf_sp_param.yd = 0;

            uint64_t left_err   = UINT_MAX;
            uint64_t right_err  = UINT_MAX;
            uint64_t top_err    = UINT_MAX;
            uint64_t bottom_err = UINT_MAX;
            if (best_direction == 0 || !context_ptr->tf_ctrls.avoid_2d_qpel) //horz
            {
                //check left
                tf_sp_param.xd = -2;
                tf_sp_param.yd = 0;
                left_err       = svt_check_position_64x64(tf_sp_param,
                                                    pcs_ptr,
                                                    context_ptr,
                                                    blk_ptr,
                                                    mv_unit,
                                                    !is_highbd ? pic_ptr_ref : &reference_ptr,
                                                    prediction_ptr,
                                                    pred,
                                                    pred_16bit,
                                                    stride_pred,
                                                    src,
                                                    src_16bit,
                                                    stride_src,
                                                    &best_mv_x,
                                                    &best_mv_y);

                //check right
                tf_sp_param.xd = 2;
                tf_sp_param.yd = 0;
                right_err      = svt_check_position_64x64(tf_sp_param,
                                                     pcs_ptr,
                                                     context_ptr,
                                                     blk_ptr,
                                                     mv_unit,
                                                     !is_highbd ? pic_ptr_ref : &reference_ptr,
                                                     prediction_ptr,
                                                     pred,
                                                     pred_16bit,
                                                     stride_pred,
                                                     src,
                                                     src_16bit,
                                                     stride_src,
                                                     &best_mv_x,
                                                     &best_mv_y);
            }

            if (best_direction == 1 || !context_ptr->tf_ctrls.avoid_2d_qpel) //vert
            {
                //check top
                tf_sp_param.xd = 0;
                tf_sp_param.yd = -2;
                top_err        = svt_check_position_64x64(tf_sp_param,
                                                   pcs_ptr,
                                                   context_ptr,
                                                   blk_ptr,
                                                   mv_unit,
                                                   !is_highbd ? pic_ptr_ref : &reference_ptr,
                                                   prediction_ptr,
                                                   pred,
                                                   pred_16bit,
                                                   stride_pred,
                                                   src,
                                                   src_16bit,
                                                   stride_src,
                                                   &best_mv_x,
                                                   &best_mv_y);

                //check bottom
                tf_sp_param.xd = 0;
                tf_sp_param.yd = 2;
                bottom_err     = svt_check_position_64x64(tf_sp_param,
                                                      pcs_ptr,
                                                      context_ptr,
                                                      blk_ptr,
                                                      mv_unit,
                                                      !is_highbd ? pic_ptr_ref : &reference_ptr,
                                                      prediction_ptr,
                                                      pred,
                                                      pred_16bit,
                                                      stride_pred,
                                                      src,
                                                      src_16bit,
                                                      stride_src,
                                                      &best_mv_x,
                                                      &best_mv_y);
            }

            min_error = MIN(left_err, MIN(top_err, MIN(right_err, bottom_err)));
            if (min_error == left_err) {
                //check left_top
                tf_sp_param.xd = -2;
                tf_sp_param.yd = -2;
                svt_check_position_64x64(tf_sp_param,
                                         pcs_ptr,
                                         context_ptr,
                                         blk_ptr,
                                         mv_unit,
                                         !is_highbd ? pic_ptr_ref : &reference_ptr,
                                         prediction_ptr,
                                         pred,
                                         pred_16bit,
                                         stride_pred,
                                         src,
                                         src_16bit,
                                         stride_src,
                                         &best_mv_x,
                                         &best_mv_y);
                //check left_bottom
                tf_sp_param.xd = -2;
                tf_sp_param.yd = 2;
                svt_check_position_64x64(tf_sp_param,
                                         pcs_ptr,
                                         context_ptr,
                                         blk_ptr,
                                         mv_unit,
                                         !is_highbd ? pic_ptr_ref : &reference_ptr,
                                         prediction_ptr,
                                         pred,
                                         pred_16bit,
                                         stride_pred,
                                         src,
                                         src_16bit,
                                         stride_src,
                                         &best_mv_x,
                                         &best_mv_y);
            } else if (min_error == top_err) {
                //check left_top
                tf_sp_param.xd = -2;
                tf_sp_param.yd = -2;
                svt_check_position_64x64(tf_sp_param,
                                         pcs_ptr,
                                         context_ptr,
                                         blk_ptr,
                                         mv_unit,
                                         !is_highbd ? pic_ptr_ref : &reference_ptr,
                                         prediction_ptr,
                                         pred,
                                         pred_16bit,
                                         stride_pred,
                                         src,
                                         src_16bit,
                                         stride_src,
                                         &best_mv_x,
                                         &best_mv_y);
                //check top_right
                tf_sp_param.xd = 2;
                tf_sp_param.yd = -2;
                svt_check_position_64x64(tf_sp_param,
                                         pcs_ptr,
                                         context_ptr,
                                         blk_ptr,
                                         mv_unit,
                                         !is_highbd ? pic_ptr_ref : &reference_ptr,
                                         prediction_ptr,
                                         pred,
                                         pred_16bit,
                                         stride_pred,
                                         src,
                                         src_16bit,
                                         stride_src,
                                         &best_mv_x,
                                         &best_mv_y);
            } else if (min_error == right_err) {
                //check top_right
                tf_sp_param.xd = 2;
                tf_sp_param.yd = -2;
                svt_check_position_64x64(tf_sp_param,
                                         pcs_ptr,
                                         context_ptr,
                                         blk_ptr,
                                         mv_unit,
                                         !is_highbd ? pic_ptr_ref : &reference_ptr,
                                         prediction_ptr,
                                         pred,
                                         pred_16bit,
                                         stride_pred,
                                         src,
                                         src_16bit,
                                         stride_src,
                                         &best_mv_x,
                                         &best_mv_y);
                //check bottom_right
                tf_sp_param.xd = 2;
                tf_sp_param.yd = 2;
                svt_check_position_64x64(tf_sp_param,
                                         pcs_ptr,
                                         context_ptr,
                                         blk_ptr,
                                         mv_unit,
                                         !is_highbd ? pic_ptr_ref : &reference_ptr,
                                         prediction_ptr,
                                         pred,
                                         pred_16bit,
                                         stride_pred,
                                         src,
                                         src_16bit,
                                         stride_src,
                                         &best_mv_x,
                                         &best_mv_y);
            } else if (min_error == bottom_err) {
                //check bottom_left
                tf_sp_param.xd = -2;
                tf_sp_param.yd = 2;
                svt_check_position_64x64(tf_sp_param,
                                         pcs_ptr,
                                         context_ptr,
                                         blk_ptr,
                                         mv_unit,
                                         !is_highbd ? pic_ptr_ref : &reference_ptr,
                                         prediction_ptr,
                                         pred,
                                         pred_16bit,
                                         stride_pred,
                                         src,
                                         src_16bit,
                                         stride_src,
                                         &best_mv_x,
                                         &best_mv_y);
                //check bottom_right
                tf_sp_param.xd = 2;
                tf_sp_param.yd = 2;
                svt_check_position_64x64(tf_sp_param,
                                         pcs_ptr,
                                         context_ptr,
                                         blk_ptr,
                                         mv_unit,
                                         !is_highbd ? pic_ptr_ref : &reference_ptr,
                                         prediction_ptr,
                                         pred,
                                         pred_16bit,
                                         stride_pred,
                                         src,
                                         src_16bit,
                                         stride_src,
                                         &best_mv_x,
                                         &best_mv_y);
            }
        }
    }

    mv_x = best_mv_x;
    mv_y = best_mv_y;
    // Perform 1/8 Pel MV Refinement
    if (pcs_ptr->tf_ctrls.eight_pel_mode) {
        tf_sp_param.subpel_pel_mode   = pcs_ptr->tf_ctrls.eight_pel_mode;
        tf_sp_param.mv_x              = mv_x;
        tf_sp_param.mv_y              = mv_y;
        tf_sp_param.interp_filters    = (uint32_t)interp_filters;
        tf_sp_param.pu_origin_x       = pu_origin_x;
        tf_sp_param.pu_origin_y       = pu_origin_y;
        tf_sp_param.local_origin_x    = local_origin_x;
        tf_sp_param.local_origin_y    = local_origin_y;
        tf_sp_param.bsize             = bsize;
        tf_sp_param.is_highbd         = is_highbd;
        tf_sp_param.encoder_bit_depth = encoder_bit_depth;
        tf_sp_param.idx_x             = 0;
        tf_sp_param.idx_y             = 0;
        for (signed short i = -1; i <= 1; i++) {
            for (signed short j = -1; j <= 1; j++) {
                tf_sp_param.xd = i;
                tf_sp_param.yd = j;
                svt_check_position_64x64(tf_sp_param,
                                         pcs_ptr,
                                         context_ptr,
                                         blk_ptr,
                                         mv_unit,
                                         !is_highbd ? pic_ptr_ref : &reference_ptr,
                                         prediction_ptr,
                                         pred,
                                         pred_16bit,
                                         stride_pred,
                                         src,
                                         src_16bit,
                                         stride_src,
                                         &best_mv_x,
                                         &best_mv_y);
            }
        }
    }

    context_ptr->tf_64x64_mv_x = best_mv_x;
    context_ptr->tf_64x64_mv_y = best_mv_y;
}
static void tf_32x32_sub_pel_search(PictureParentControlSet *pcs_ptr, MeContext *context_ptr,
                                    PictureParentControlSet *pcs_ref,
                                    EbPictureBufferDesc *pic_ptr_ref, EbByte *pred,
                                    uint16_t **pred_16bit, uint32_t *stride_pred, EbByte *src,
                                    uint16_t **src_16bit, uint32_t *stride_src,
                                    uint32_t sb_origin_x, uint32_t sb_origin_y, uint32_t ss_x,
                                    int encoder_bit_depth) {
    InterpFilters interp_filters;
    if (context_ptr->tf_ctrls.use_2tap)
        interp_filters = av1_make_interp_filters(BILINEAR, BILINEAR);
    else
        interp_filters = av1_make_interp_filters(EIGHTTAP_REGULAR, EIGHTTAP_REGULAR);

    Bool is_highbd = (encoder_bit_depth == 8) ? (uint8_t)FALSE : (uint8_t)TRUE;
    BlkStruct blk_struct;
    MacroBlockD av1xd;
    blk_struct.av1xd = &av1xd;
    MvUnit mv_unit;
    mv_unit.pred_direction = UNI_PRED_LIST_0;
    EbPictureBufferDesc reference_ptr;
    EbPictureBufferDesc prediction_ptr;
    UNUSED(ss_x);
    prediction_ptr.origin_x  = 0;
    prediction_ptr.origin_y  = 0;
    prediction_ptr.stride_y  = BW;
    prediction_ptr.stride_cb = (uint16_t)BW >> ss_x;
    prediction_ptr.stride_cr = (uint16_t)BW >> ss_x;
    if (!is_highbd) {
        assert(src[C_Y] != NULL);
        if (context_ptr->tf_chroma) {
            assert(src[C_U] != NULL);
            assert(src[C_V] != NULL);
        }
        prediction_ptr.buffer_y  = pred[C_Y];
        prediction_ptr.buffer_cb = pred[C_U];
        prediction_ptr.buffer_cr = pred[C_V];
    } else {
        assert(src_16bit[C_Y] != NULL);
        if (context_ptr->tf_chroma) {
            assert(src_16bit[C_U] != NULL);
            assert(src_16bit[C_V] != NULL);
        }
        prediction_ptr.buffer_y  = (uint8_t *)pred_16bit[C_Y];
        prediction_ptr.buffer_cb = (uint8_t *)pred_16bit[C_U];
        prediction_ptr.buffer_cr = (uint8_t *)pred_16bit[C_V];
        reference_ptr.buffer_y   = (uint8_t *)pcs_ref->altref_buffer_highbd[C_Y];
        reference_ptr.buffer_cb  = (uint8_t *)pcs_ref->altref_buffer_highbd[C_U];
        reference_ptr.buffer_cr  = (uint8_t *)pcs_ref->altref_buffer_highbd[C_V];
        reference_ptr.origin_x   = pic_ptr_ref->origin_x;
        reference_ptr.origin_y   = pic_ptr_ref->origin_y;
        reference_ptr.stride_y   = pic_ptr_ref->stride_y;
        reference_ptr.stride_cb  = pic_ptr_ref->stride_cb;
        reference_ptr.stride_cr  = pic_ptr_ref->stride_cr;
        reference_ptr.width      = pic_ptr_ref->width;
        reference_ptr.height     = pic_ptr_ref->height;
        reference_ptr.buffer_bit_inc_y  = NULL;
        reference_ptr.buffer_bit_inc_cb = NULL;
        reference_ptr.buffer_bit_inc_cr = NULL;

    }
    uint32_t bsize          = 32;
    uint32_t idx_32x32      = context_ptr->idx_32x32;
    uint32_t idx_x          = idx_32x32 & 0x1;
    uint32_t idx_y          = idx_32x32 >> 1;
    uint16_t local_origin_x = idx_x * bsize;
    uint16_t local_origin_y = idx_y * bsize;
    uint16_t pu_origin_x    = sb_origin_x + local_origin_x;
    uint16_t pu_origin_y    = sb_origin_y + local_origin_y;
    int32_t mirow = pu_origin_y >> MI_SIZE_LOG2;
    int32_t micol = pu_origin_x >> MI_SIZE_LOG2;
    blk_struct.mds_idx = get_mds_idx(local_origin_x,
                                     local_origin_y,
                                     bsize,
                                     pcs_ptr->scs_ptr->seq_header.sb_size == BLOCK_128X128);

    const int32_t bw                    = mi_size_wide[BLOCK_32X32];
    const int32_t bh                    = mi_size_high[BLOCK_32X32];
    blk_struct.av1xd->mb_to_top_edge    = -(int32_t)((mirow * MI_SIZE) * 8);
    blk_struct.av1xd->mb_to_bottom_edge = ((pcs_ptr->av1_cm->mi_rows - bw - mirow) * MI_SIZE) * 8;
    blk_struct.av1xd->mb_to_left_edge   = -(int32_t)((micol * MI_SIZE) * 8);
    blk_struct.av1xd->mb_to_right_edge  = ((pcs_ptr->av1_cm->mi_cols - bh - micol) * MI_SIZE) * 8;

    BlkStruct *blk_ptr = &blk_struct;
    uint32_t mv_index = idx_32x32;
    mv_unit.mv->x     = _MVXT(context_ptr->p_best_mv32x32[mv_index]);
    mv_unit.mv->y     = _MVYT(context_ptr->p_best_mv32x32[mv_index]);
    // AV1 MVs are always in 1/8th pel precision.
    mv_unit.mv->x                                = mv_unit.mv->x << 1;
    mv_unit.mv->y                                = mv_unit.mv->y << 1;
    context_ptr->tf_32x32_block_error[idx_32x32] = INT_MAX;
    signed short            mv_x      = (_MVXT(context_ptr->p_best_mv32x32[mv_index])) << 1;
    signed short            mv_y      = (_MVYT(context_ptr->p_best_mv32x32[mv_index])) << 1;
    signed short            best_mv_x = mv_x;
    signed short            best_mv_y = mv_y;
    TF_SUBPEL_SEARCH_PARAMS tf_sp_param;
    tf_sp_param.subsampling_shift = pcs_ptr->tf_ctrls.sub_sampling_shift;
    uint64_t center_err       = UINT_MAX;
    uint64_t min_error        = UINT_MAX;
    uint64_t left_bottom_err  = UINT_MAX;
    uint64_t left_top_err     = UINT_MAX;
    uint64_t top_right_err    = UINT_MAX;
    uint64_t bottom_right_err = UINT_MAX;

    uint8_t best_direction = 0; //horz
    // No refinement.
    if (!pcs_ptr->tf_ctrls.half_pel_mode && !pcs_ptr->tf_ctrls.quarter_pel_mode &&
        !pcs_ptr->tf_ctrls.eight_pel_mode) {
        tf_sp_param.subpel_pel_mode   = pcs_ptr->tf_ctrls.half_pel_mode;
        tf_sp_param.mv_x              = mv_x;
        tf_sp_param.mv_y              = mv_y;
        tf_sp_param.interp_filters    = (uint32_t)interp_filters;
        tf_sp_param.pu_origin_x       = pu_origin_x;
        tf_sp_param.pu_origin_y       = pu_origin_y;
        tf_sp_param.local_origin_x    = local_origin_x;
        tf_sp_param.local_origin_y    = local_origin_y;
        tf_sp_param.bsize             = bsize;
        tf_sp_param.is_highbd         = is_highbd;
        tf_sp_param.encoder_bit_depth = encoder_bit_depth;
        tf_sp_param.idx_x             = idx_x;
        tf_sp_param.idx_y             = idx_y;
        tf_sp_param.xd                = 0;
        tf_sp_param.yd                = 0;
        center_err                    = svt_check_position(tf_sp_param,
                                        pcs_ptr,
                                        context_ptr,
                                        blk_ptr,
                                        mv_unit,
                                        !is_highbd ? pic_ptr_ref : &reference_ptr,
                                        prediction_ptr,
                                        pred,
                                        pred_16bit,
                                        stride_pred,
                                        src,
                                        src_16bit,
                                        stride_src,
                                        &best_mv_x,
                                        &best_mv_y);
    }

    // Perform 1/2 Pel MV Refinement
    if (pcs_ptr->tf_ctrls.half_pel_mode) {
        tf_sp_param.subpel_pel_mode   = pcs_ptr->tf_ctrls.half_pel_mode;
        tf_sp_param.mv_x              = mv_x;
        tf_sp_param.mv_y              = mv_y;
        tf_sp_param.interp_filters    = (uint32_t)interp_filters;
        tf_sp_param.pu_origin_x       = pu_origin_x;
        tf_sp_param.pu_origin_y       = pu_origin_y;
        tf_sp_param.local_origin_x    = local_origin_x;
        tf_sp_param.local_origin_y    = local_origin_y;
        tf_sp_param.bsize             = bsize;
        tf_sp_param.is_highbd         = is_highbd;
        tf_sp_param.encoder_bit_depth = encoder_bit_depth;
        tf_sp_param.idx_x             = idx_x;
        tf_sp_param.idx_y             = idx_y;
        if (pcs_ptr->tf_ctrls.half_pel_mode == 1) {
            for (signed short i = -4; i <= 4; i = i + 4) {
                for (signed short j = -4; j <= 4; j = j + 4) {
                    tf_sp_param.xd   = i;
                    tf_sp_param.yd   = j;
                    uint64_t cur_err = svt_check_position(tf_sp_param,
                                                          pcs_ptr,
                                                          context_ptr,
                                                          blk_ptr,
                                                          mv_unit,
                                                          !is_highbd ? pic_ptr_ref : &reference_ptr,
                                                          prediction_ptr,
                                                          pred,
                                                          pred_16bit,
                                                          stride_pred,
                                                          src,
                                                          src_16bit,
                                                          stride_src,
                                                          &best_mv_x,
                                                          &best_mv_y);
                    if (cur_err < min_error)
                        min_error = cur_err;
                }
            }
        } else {
            //check center
            tf_sp_param.xd = 0;
            tf_sp_param.yd = 0;
            center_err     = svt_check_position(tf_sp_param,
                                            pcs_ptr,
                                            context_ptr,
                                            blk_ptr,
                                            mv_unit,
                                            !is_highbd ? pic_ptr_ref : &reference_ptr,
                                            prediction_ptr,
                                            pred,
                                            pred_16bit,
                                            stride_pred,
                                            src,
                                            src_16bit,
                                            stride_src,
                                            &best_mv_x,
                                            &best_mv_y);

            //check left
            tf_sp_param.xd    = -4;
            tf_sp_param.yd    = 0;
            uint64_t left_err = svt_check_position(tf_sp_param,
                                                   pcs_ptr,
                                                   context_ptr,
                                                   blk_ptr,
                                                   mv_unit,
                                                   !is_highbd ? pic_ptr_ref : &reference_ptr,
                                                   prediction_ptr,
                                                   pred,
                                                   pred_16bit,
                                                   stride_pred,
                                                   src,
                                                   src_16bit,
                                                   stride_src,
                                                   &best_mv_x,
                                                   &best_mv_y);

            //check top
            tf_sp_param.xd   = 0;
            tf_sp_param.yd   = -4;
            uint64_t top_err = svt_check_position(tf_sp_param,
                                                  pcs_ptr,
                                                  context_ptr,
                                                  blk_ptr,
                                                  mv_unit,
                                                  !is_highbd ? pic_ptr_ref : &reference_ptr,
                                                  prediction_ptr,
                                                  pred,
                                                  pred_16bit,
                                                  stride_pred,
                                                  src,
                                                  src_16bit,
                                                  stride_src,
                                                  &best_mv_x,
                                                  &best_mv_y);

            //check right
            tf_sp_param.xd     = 4;
            tf_sp_param.yd     = 0;
            uint64_t right_err = svt_check_position(tf_sp_param,
                                                    pcs_ptr,
                                                    context_ptr,
                                                    blk_ptr,
                                                    mv_unit,
                                                    !is_highbd ? pic_ptr_ref : &reference_ptr,
                                                    prediction_ptr,
                                                    pred,
                                                    pred_16bit,
                                                    stride_pred,
                                                    src,
                                                    src_16bit,
                                                    stride_src,
                                                    &best_mv_x,
                                                    &best_mv_y);

            //check bottom
            tf_sp_param.xd      = 0;
            tf_sp_param.yd      = 4;
            uint64_t bottom_err = svt_check_position(tf_sp_param,
                                                     pcs_ptr,
                                                     context_ptr,
                                                     blk_ptr,
                                                     mv_unit,
                                                     !is_highbd ? pic_ptr_ref : &reference_ptr,
                                                     prediction_ptr,
                                                     pred,
                                                     pred_16bit,
                                                     stride_pred,
                                                     src,
                                                     src_16bit,
                                                     stride_src,
                                                     &best_mv_x,
                                                     &best_mv_y);

            min_error = MIN(left_err, MIN(top_err, MIN(right_err, bottom_err)));

            if (min_error == top_err || min_error == bottom_err)
                best_direction = 1; //vert

            if (min_error == left_err) {
                //check left_top
                tf_sp_param.xd = -4;
                tf_sp_param.yd = -4;
                left_top_err   = svt_check_position(tf_sp_param,
                                                  pcs_ptr,
                                                  context_ptr,
                                                  blk_ptr,
                                                  mv_unit,
                                                  !is_highbd ? pic_ptr_ref : &reference_ptr,
                                                  prediction_ptr,
                                                  pred,
                                                  pred_16bit,
                                                  stride_pred,
                                                  src,
                                                  src_16bit,
                                                  stride_src,
                                                  &best_mv_x,
                                                  &best_mv_y);
                //check left_bottom
                tf_sp_param.xd  = -4;
                tf_sp_param.yd  = 4;
                left_bottom_err = svt_check_position(tf_sp_param,
                                                     pcs_ptr,
                                                     context_ptr,
                                                     blk_ptr,
                                                     mv_unit,
                                                     !is_highbd ? pic_ptr_ref : &reference_ptr,
                                                     prediction_ptr,
                                                     pred,
                                                     pred_16bit,
                                                     stride_pred,
                                                     src,
                                                     src_16bit,
                                                     stride_src,
                                                     &best_mv_x,
                                                     &best_mv_y);
            } else if (min_error == top_err) {
                //check left_top
                tf_sp_param.xd = -4;
                tf_sp_param.yd = -4;
                left_top_err   = svt_check_position(tf_sp_param,
                                                  pcs_ptr,
                                                  context_ptr,
                                                  blk_ptr,
                                                  mv_unit,
                                                  !is_highbd ? pic_ptr_ref : &reference_ptr,
                                                  prediction_ptr,
                                                  pred,
                                                  pred_16bit,
                                                  stride_pred,
                                                  src,
                                                  src_16bit,
                                                  stride_src,
                                                  &best_mv_x,
                                                  &best_mv_y);
                //check top_right
                tf_sp_param.xd = 4;
                tf_sp_param.yd = -4;
                top_right_err  = svt_check_position(tf_sp_param,
                                                   pcs_ptr,
                                                   context_ptr,
                                                   blk_ptr,
                                                   mv_unit,
                                                   !is_highbd ? pic_ptr_ref : &reference_ptr,
                                                   prediction_ptr,
                                                   pred,
                                                   pred_16bit,
                                                   stride_pred,
                                                   src,
                                                   src_16bit,
                                                   stride_src,
                                                   &best_mv_x,
                                                   &best_mv_y);
            } else if (min_error == right_err) {
                //check top_right
                tf_sp_param.xd = 4;
                tf_sp_param.yd = -4;
                top_right_err  = svt_check_position(tf_sp_param,
                                                   pcs_ptr,
                                                   context_ptr,
                                                   blk_ptr,
                                                   mv_unit,
                                                   !is_highbd ? pic_ptr_ref : &reference_ptr,
                                                   prediction_ptr,
                                                   pred,
                                                   pred_16bit,
                                                   stride_pred,
                                                   src,
                                                   src_16bit,
                                                   stride_src,
                                                   &best_mv_x,
                                                   &best_mv_y);
                //check bottom_right
                tf_sp_param.xd   = 4;
                tf_sp_param.yd   = 4;
                bottom_right_err = svt_check_position(tf_sp_param,
                                                      pcs_ptr,
                                                      context_ptr,
                                                      blk_ptr,
                                                      mv_unit,
                                                      !is_highbd ? pic_ptr_ref : &reference_ptr,
                                                      prediction_ptr,
                                                      pred,
                                                      pred_16bit,
                                                      stride_pred,
                                                      src,
                                                      src_16bit,
                                                      stride_src,
                                                      &best_mv_x,
                                                      &best_mv_y);
            } else if (min_error == bottom_err) {
                //check bottom_left
                tf_sp_param.xd  = -4;
                tf_sp_param.yd  = 4;
                left_bottom_err = svt_check_position(tf_sp_param,
                                                     pcs_ptr,
                                                     context_ptr,
                                                     blk_ptr,
                                                     mv_unit,
                                                     !is_highbd ? pic_ptr_ref : &reference_ptr,
                                                     prediction_ptr,
                                                     pred,
                                                     pred_16bit,
                                                     stride_pred,
                                                     src,
                                                     src_16bit,
                                                     stride_src,
                                                     &best_mv_x,
                                                     &best_mv_y);
                //check bottom_right
                tf_sp_param.xd   = 4;
                tf_sp_param.yd   = 4;
                bottom_right_err = svt_check_position(tf_sp_param,
                                                      pcs_ptr,
                                                      context_ptr,
                                                      blk_ptr,
                                                      mv_unit,
                                                      !is_highbd ? pic_ptr_ref : &reference_ptr,
                                                      prediction_ptr,
                                                      pred,
                                                      pred_16bit,
                                                      stride_pred,
                                                      src,
                                                      src_16bit,
                                                      stride_src,
                                                      &best_mv_x,
                                                      &best_mv_y);
            }
        }
        min_error = MIN(
            min_error,
            MIN(left_top_err, MIN(left_bottom_err, MIN(top_right_err, bottom_right_err))));
    }

    mv_x = best_mv_x;
    mv_y = best_mv_y;
    // Perform 1/4 Pel MV Refinement
    if (pcs_ptr->tf_ctrls.quarter_pel_mode) {
        tf_sp_param.subpel_pel_mode   = pcs_ptr->tf_ctrls.quarter_pel_mode;
        tf_sp_param.mv_x              = mv_x;
        tf_sp_param.mv_y              = mv_y;
        tf_sp_param.interp_filters    = (uint32_t)interp_filters;
        tf_sp_param.pu_origin_x       = pu_origin_x;
        tf_sp_param.pu_origin_y       = pu_origin_y;
        tf_sp_param.local_origin_x    = local_origin_x;
        tf_sp_param.local_origin_y    = local_origin_y;
        tf_sp_param.bsize             = bsize;
        tf_sp_param.is_highbd         = is_highbd;
        tf_sp_param.encoder_bit_depth = encoder_bit_depth;
        tf_sp_param.idx_x             = idx_x;
        tf_sp_param.idx_y             = idx_y;
        if (pcs_ptr->tf_ctrls.quarter_pel_mode == 1) {
            for (signed short i = -2; i <= 2; i = i + 2) {
                for (signed short j = -2; j <= 2; j = j + 2) {
                    tf_sp_param.xd = i;
                    tf_sp_param.yd = j;
                    center_err     = svt_check_position(tf_sp_param,
                                                    pcs_ptr,
                                                    context_ptr,
                                                    blk_ptr,
                                                    mv_unit,
                                                    !is_highbd ? pic_ptr_ref : &reference_ptr,
                                                    prediction_ptr,
                                                    pred,
                                                    pred_16bit,
                                                    stride_pred,
                                                    src,
                                                    src_16bit,
                                                    stride_src,
                                                    &best_mv_x,
                                                    &best_mv_y);
                }
            }
        } else if (min_error < center_err) {
            //check center
            tf_sp_param.xd = 0;
            tf_sp_param.yd = 0;

            uint64_t left_err   = UINT_MAX;
            uint64_t right_err  = UINT_MAX;
            uint64_t top_err    = UINT_MAX;
            uint64_t bottom_err = UINT_MAX;
            if (best_direction == 0 || !context_ptr->tf_ctrls.avoid_2d_qpel) //horz
            {
                //check left
                tf_sp_param.xd = -2;
                tf_sp_param.yd = 0;
                left_err       = svt_check_position(tf_sp_param,
                                              pcs_ptr,
                                              context_ptr,
                                              blk_ptr,
                                              mv_unit,
                                              !is_highbd ? pic_ptr_ref : &reference_ptr,
                                              prediction_ptr,
                                              pred,
                                              pred_16bit,
                                              stride_pred,
                                              src,
                                              src_16bit,
                                              stride_src,
                                              &best_mv_x,
                                              &best_mv_y);

                //check right
                tf_sp_param.xd = 2;
                tf_sp_param.yd = 0;
                right_err      = svt_check_position(tf_sp_param,
                                               pcs_ptr,
                                               context_ptr,
                                               blk_ptr,
                                               mv_unit,
                                               !is_highbd ? pic_ptr_ref : &reference_ptr,
                                               prediction_ptr,
                                               pred,
                                               pred_16bit,
                                               stride_pred,
                                               src,
                                               src_16bit,
                                               stride_src,
                                               &best_mv_x,
                                               &best_mv_y);
            }

            if (best_direction == 1 || !context_ptr->tf_ctrls.avoid_2d_qpel) //vert
            {
                //check top
                tf_sp_param.xd = 0;
                tf_sp_param.yd = -2;
                top_err        = svt_check_position(tf_sp_param,
                                             pcs_ptr,
                                             context_ptr,
                                             blk_ptr,
                                             mv_unit,
                                             !is_highbd ? pic_ptr_ref : &reference_ptr,
                                             prediction_ptr,
                                             pred,
                                             pred_16bit,
                                             stride_pred,
                                             src,
                                             src_16bit,
                                             stride_src,
                                             &best_mv_x,
                                             &best_mv_y);

                //check bottom
                tf_sp_param.xd = 0;
                tf_sp_param.yd = 2;
                bottom_err     = svt_check_position(tf_sp_param,
                                                pcs_ptr,
                                                context_ptr,
                                                blk_ptr,
                                                mv_unit,
                                                !is_highbd ? pic_ptr_ref : &reference_ptr,
                                                prediction_ptr,
                                                pred,
                                                pred_16bit,
                                                stride_pred,
                                                src,
                                                src_16bit,
                                                stride_src,
                                                &best_mv_x,
                                                &best_mv_y);
            }

            min_error = MIN(left_err, MIN(top_err, MIN(right_err, bottom_err)));
            if (min_error == left_err) {
                //check left_top
                tf_sp_param.xd = -2;
                tf_sp_param.yd = -2;
                svt_check_position(tf_sp_param,
                                   pcs_ptr,
                                   context_ptr,
                                   blk_ptr,
                                   mv_unit,
                                   !is_highbd ? pic_ptr_ref : &reference_ptr,
                                   prediction_ptr,
                                   pred,
                                   pred_16bit,
                                   stride_pred,
                                   src,
                                   src_16bit,
                                   stride_src,
                                   &best_mv_x,
                                   &best_mv_y);
                //check left_bottom
                tf_sp_param.xd = -2;
                tf_sp_param.yd = 2;
                svt_check_position(tf_sp_param,
                                   pcs_ptr,
                                   context_ptr,
                                   blk_ptr,
                                   mv_unit,
                                   !is_highbd ? pic_ptr_ref : &reference_ptr,
                                   prediction_ptr,
                                   pred,
                                   pred_16bit,
                                   stride_pred,
                                   src,
                                   src_16bit,
                                   stride_src,
                                   &best_mv_x,
                                   &best_mv_y);
            } else if (min_error == top_err) {
                //check left_top
                tf_sp_param.xd = -2;
                tf_sp_param.yd = -2;
                svt_check_position(tf_sp_param,
                                   pcs_ptr,
                                   context_ptr,
                                   blk_ptr,
                                   mv_unit,
                                   !is_highbd ? pic_ptr_ref : &reference_ptr,
                                   prediction_ptr,
                                   pred,
                                   pred_16bit,
                                   stride_pred,
                                   src,
                                   src_16bit,
                                   stride_src,
                                   &best_mv_x,
                                   &best_mv_y);
                //check top_right
                tf_sp_param.xd = 2;
                tf_sp_param.yd = -2;
                svt_check_position(tf_sp_param,
                                   pcs_ptr,
                                   context_ptr,
                                   blk_ptr,
                                   mv_unit,
                                   !is_highbd ? pic_ptr_ref : &reference_ptr,
                                   prediction_ptr,
                                   pred,
                                   pred_16bit,
                                   stride_pred,
                                   src,
                                   src_16bit,
                                   stride_src,
                                   &best_mv_x,
                                   &best_mv_y);
            } else if (min_error == right_err) {
                //check top_right
                tf_sp_param.xd = 2;
                tf_sp_param.yd = -2;
                svt_check_position(tf_sp_param,
                                   pcs_ptr,
                                   context_ptr,
                                   blk_ptr,
                                   mv_unit,
                                   !is_highbd ? pic_ptr_ref : &reference_ptr,
                                   prediction_ptr,
                                   pred,
                                   pred_16bit,
                                   stride_pred,
                                   src,
                                   src_16bit,
                                   stride_src,
                                   &best_mv_x,
                                   &best_mv_y);
                //check bottom_right
                tf_sp_param.xd = 2;
                tf_sp_param.yd = 2;
                svt_check_position(tf_sp_param,
                                   pcs_ptr,
                                   context_ptr,
                                   blk_ptr,
                                   mv_unit,
                                   !is_highbd ? pic_ptr_ref : &reference_ptr,
                                   prediction_ptr,
                                   pred,
                                   pred_16bit,
                                   stride_pred,
                                   src,
                                   src_16bit,
                                   stride_src,
                                   &best_mv_x,
                                   &best_mv_y);
            } else if (min_error == bottom_err) {
                //check bottom_left
                tf_sp_param.xd = -2;
                tf_sp_param.yd = 2;
                svt_check_position(tf_sp_param,
                                   pcs_ptr,
                                   context_ptr,
                                   blk_ptr,
                                   mv_unit,
                                   !is_highbd ? pic_ptr_ref : &reference_ptr,
                                   prediction_ptr,
                                   pred,
                                   pred_16bit,
                                   stride_pred,
                                   src,
                                   src_16bit,
                                   stride_src,
                                   &best_mv_x,
                                   &best_mv_y);
                //check bottom_right
                tf_sp_param.xd = 2;
                tf_sp_param.yd = 2;
                svt_check_position(tf_sp_param,
                                   pcs_ptr,
                                   context_ptr,
                                   blk_ptr,
                                   mv_unit,
                                   !is_highbd ? pic_ptr_ref : &reference_ptr,
                                   prediction_ptr,
                                   pred,
                                   pred_16bit,
                                   stride_pred,
                                   src,
                                   src_16bit,
                                   stride_src,
                                   &best_mv_x,
                                   &best_mv_y);
            }
        }
    }

    mv_x = best_mv_x;
    mv_y = best_mv_y;
    // Perform 1/8 Pel MV Refinement
    if (pcs_ptr->tf_ctrls.eight_pel_mode) {
        tf_sp_param.subpel_pel_mode   = pcs_ptr->tf_ctrls.eight_pel_mode;
        tf_sp_param.mv_x              = mv_x;
        tf_sp_param.mv_y              = mv_y;
        tf_sp_param.interp_filters    = (uint32_t)interp_filters;
        tf_sp_param.pu_origin_x       = pu_origin_x;
        tf_sp_param.pu_origin_y       = pu_origin_y;
        tf_sp_param.local_origin_x    = local_origin_x;
        tf_sp_param.local_origin_y    = local_origin_y;
        tf_sp_param.bsize             = bsize;
        tf_sp_param.is_highbd         = is_highbd;
        tf_sp_param.encoder_bit_depth = encoder_bit_depth;
        tf_sp_param.idx_x             = idx_x;
        tf_sp_param.idx_y             = idx_y;
        for (signed short i = -1; i <= 1; i++) {
            for (signed short j = -1; j <= 1; j++) {
                tf_sp_param.xd = i;
                tf_sp_param.yd = j;
                svt_check_position(tf_sp_param,
                                   pcs_ptr,
                                   context_ptr,
                                   blk_ptr,
                                   mv_unit,
                                   !is_highbd ? pic_ptr_ref : &reference_ptr,
                                   prediction_ptr,
                                   pred,
                                   pred_16bit,
                                   stride_pred,
                                   src,
                                   src_16bit,
                                   stride_src,
                                   &best_mv_x,
                                   &best_mv_y);
            }
        }
    }
    context_ptr->tf_32x32_mv_x[idx_32x32] = best_mv_x;
    context_ptr->tf_32x32_mv_y[idx_32x32] = best_mv_y;
}

static void tf_64x64_inter_prediction(PictureParentControlSet *pcs_ptr, MeContext *context_ptr,
                                      PictureParentControlSet *pcs_ref,
                                      EbPictureBufferDesc *pic_ptr_ref, EbByte *pred,
                                      uint16_t **pred_16bit, uint32_t sb_origin_x,
                                      uint32_t sb_origin_y, uint32_t ss_x, int encoder_bit_depth) {
    SequenceControlSet *scs_ptr = pcs_ptr->scs_ptr;
    const InterpFilters interp_filters = av1_make_interp_filters(MULTITAP_SHARP, MULTITAP_SHARP);

    Bool is_highbd = (encoder_bit_depth == 8) ? (uint8_t)FALSE : (uint8_t)TRUE;

    BlkStruct   blk_ptr;
    MacroBlockD av1xd;
    blk_ptr.av1xd = &av1xd;
    MvUnit mv_unit;
    mv_unit.pred_direction = UNI_PRED_LIST_0;

    EbPictureBufferDesc reference_ptr;
    EbPictureBufferDesc prediction_ptr;

    prediction_ptr.origin_x  = 0;
    prediction_ptr.origin_y  = 0;
    prediction_ptr.stride_y  = BW;
    prediction_ptr.stride_cb = (uint16_t)BW >> ss_x;
    prediction_ptr.stride_cr = (uint16_t)BW >> ss_x;

    if (!is_highbd) {
        prediction_ptr.buffer_y  = pred[C_Y];
        prediction_ptr.buffer_cb = pred[C_U];
        prediction_ptr.buffer_cr = pred[C_V];
    } else {
        prediction_ptr.buffer_y  = (uint8_t *)pred_16bit[C_Y];
        prediction_ptr.buffer_cb = (uint8_t *)pred_16bit[C_U];
        prediction_ptr.buffer_cr = (uint8_t *)pred_16bit[C_V];
        reference_ptr.buffer_y   = (uint8_t *)pcs_ref->altref_buffer_highbd[C_Y];
        reference_ptr.buffer_cb  = (uint8_t *)pcs_ref->altref_buffer_highbd[C_U];
        reference_ptr.buffer_cr  = (uint8_t *)pcs_ref->altref_buffer_highbd[C_V];
        reference_ptr.origin_x   = pic_ptr_ref->origin_x;
        reference_ptr.origin_y   = pic_ptr_ref->origin_y;
        reference_ptr.stride_y   = pic_ptr_ref->stride_y;
        reference_ptr.stride_cb  = pic_ptr_ref->stride_cb;
        reference_ptr.stride_cr  = pic_ptr_ref->stride_cr;
        reference_ptr.width      = pic_ptr_ref->width;
        reference_ptr.height     = pic_ptr_ref->height;
        reference_ptr.buffer_bit_inc_y  = NULL;
        reference_ptr.buffer_bit_inc_cb = NULL;
        reference_ptr.buffer_bit_inc_cr = NULL;

    }

    uint32_t bsize = 64;

    uint16_t local_origin_x = 0;
    uint16_t local_origin_y = 0;
    uint16_t pu_origin_x    = sb_origin_x + local_origin_x;
    uint16_t pu_origin_y    = sb_origin_y + local_origin_y;
    int32_t mirow = pu_origin_y >> MI_SIZE_LOG2;
    int32_t micol = pu_origin_x >> MI_SIZE_LOG2;
    blk_ptr.mds_idx = get_mds_idx(local_origin_x,
                                  local_origin_y,
                                  bsize,
                                  pcs_ptr->scs_ptr->seq_header.sb_size == BLOCK_128X128);

    const int32_t bw                 = mi_size_wide[BLOCK_32X32];
    const int32_t bh                 = mi_size_high[BLOCK_32X32];
    blk_ptr.av1xd->mb_to_top_edge    = -(int32_t)((mirow * MI_SIZE) * 8);
    blk_ptr.av1xd->mb_to_bottom_edge = ((pcs_ptr->av1_cm->mi_rows - bw - mirow) * MI_SIZE) * 8;
    blk_ptr.av1xd->mb_to_left_edge   = -(int32_t)((micol * MI_SIZE) * 8);
    blk_ptr.av1xd->mb_to_right_edge  = ((pcs_ptr->av1_cm->mi_cols - bh - micol) * MI_SIZE) * 8;

    // Perform final pass using the 1/8 MV
    // AV1 MVs are always in 1/8th pel precision.
    mv_unit.mv->x = context_ptr->tf_64x64_mv_x;
    mv_unit.mv->y = context_ptr->tf_64x64_mv_y;

    av1_inter_prediction(
        scs_ptr,
        NULL, //pcs_ptr,
        (uint32_t)interp_filters,
        &blk_ptr,
        0, //ref_frame_type,
        &mv_unit,
        0, //use_intrabc,
        SIMPLE_TRANSLATION,
        0,
        0,
        1, //compound_idx not used
        NULL, // interinter_comp not used
        NULL,
        NULL,
        NULL,
        0,
        0,
        0,
        0,
        pu_origin_x,
        pu_origin_y,
        bsize,
        bsize,
        !is_highbd ? pic_ptr_ref : &reference_ptr,
        NULL, //ref_pic_list1,
        &prediction_ptr,
        local_origin_x,
        local_origin_y,
        context_ptr->tf_chroma ? PICTURE_BUFFER_DESC_FULL_MASK : PICTURE_BUFFER_DESC_LUMA_MASK,
        (uint8_t)encoder_bit_depth,
        is_highbd);
}
static void tf_32x32_inter_prediction(PictureParentControlSet *pcs_ptr, MeContext *context_ptr,
                                      PictureParentControlSet *pcs_ref,
                                      EbPictureBufferDesc *pic_ptr_ref, EbByte *pred,
                                      uint16_t **pred_16bit, uint32_t sb_origin_x,
                                      uint32_t sb_origin_y, uint32_t ss_x, int encoder_bit_depth) {
    SequenceControlSet *scs_ptr = pcs_ptr->scs_ptr;
    const InterpFilters interp_filters = av1_make_interp_filters(MULTITAP_SHARP, MULTITAP_SHARP);

    Bool is_highbd = (encoder_bit_depth == 8) ? (uint8_t)FALSE : (uint8_t)TRUE;

    BlkStruct   blk_ptr;
    MacroBlockD av1xd;
    blk_ptr.av1xd = &av1xd;
    MvUnit mv_unit;
    mv_unit.pred_direction = UNI_PRED_LIST_0;

    EbPictureBufferDesc reference_ptr;
    EbPictureBufferDesc prediction_ptr;

    prediction_ptr.origin_x  = 0;
    prediction_ptr.origin_y  = 0;
    prediction_ptr.stride_y  = BW;
    prediction_ptr.stride_cb = (uint16_t)BW >> ss_x;
    prediction_ptr.stride_cr = (uint16_t)BW >> ss_x;

    if (!is_highbd) {
        prediction_ptr.buffer_y  = pred[C_Y];
        prediction_ptr.buffer_cb = pred[C_U];
        prediction_ptr.buffer_cr = pred[C_V];
    } else {
        prediction_ptr.buffer_y  = (uint8_t *)pred_16bit[C_Y];
        prediction_ptr.buffer_cb = (uint8_t *)pred_16bit[C_U];
        prediction_ptr.buffer_cr = (uint8_t *)pred_16bit[C_V];
        reference_ptr.buffer_y   = (uint8_t *)pcs_ref->altref_buffer_highbd[C_Y];
        reference_ptr.buffer_cb  = (uint8_t *)pcs_ref->altref_buffer_highbd[C_U];
        reference_ptr.buffer_cr  = (uint8_t *)pcs_ref->altref_buffer_highbd[C_V];
        reference_ptr.origin_x   = pic_ptr_ref->origin_x;
        reference_ptr.origin_y   = pic_ptr_ref->origin_y;
        reference_ptr.stride_y   = pic_ptr_ref->stride_y;
        reference_ptr.stride_cb  = pic_ptr_ref->stride_cb;
        reference_ptr.stride_cr  = pic_ptr_ref->stride_cr;
        reference_ptr.width      = pic_ptr_ref->width;
        reference_ptr.height     = pic_ptr_ref->height;
        reference_ptr.buffer_bit_inc_y  = NULL;
        reference_ptr.buffer_bit_inc_cb = NULL;
        reference_ptr.buffer_bit_inc_cr = NULL;

    }

    uint32_t idx_32x32 = context_ptr->idx_32x32;
    if (context_ptr->tf_32x32_block_split_flag[idx_32x32]) {
        uint32_t bsize = 16;

        for (uint32_t idx_16x16 = 0; idx_16x16 < 4; idx_16x16++) {
            uint32_t pu_index = idx_32x32_to_idx_16x16[idx_32x32][idx_16x16];

            uint32_t idx_y          = subblock_xy_16x16[pu_index][0];
            uint32_t idx_x          = subblock_xy_16x16[pu_index][1];
            uint16_t local_origin_x = idx_x * bsize;
            uint16_t local_origin_y = idx_y * bsize;
            uint16_t pu_origin_x    = sb_origin_x + local_origin_x;
            uint16_t pu_origin_y    = sb_origin_y + local_origin_y;
            int32_t mirow = pu_origin_y >> MI_SIZE_LOG2;
            int32_t micol = pu_origin_x >> MI_SIZE_LOG2;
            blk_ptr.mds_idx = get_mds_idx(local_origin_x,
                                          local_origin_y,
                                          bsize,
                                          pcs_ptr->scs_ptr->seq_header.sb_size == BLOCK_128X128);

            const int32_t bw                 = mi_size_wide[BLOCK_16X16];
            const int32_t bh                 = mi_size_high[BLOCK_16X16];
            blk_ptr.av1xd->mb_to_top_edge    = -(int32_t)((mirow * MI_SIZE) * 8);
            blk_ptr.av1xd->mb_to_bottom_edge = ((pcs_ptr->av1_cm->mi_rows - bw - mirow) * MI_SIZE) *
                8;
            blk_ptr.av1xd->mb_to_left_edge  = -(int32_t)((micol * MI_SIZE) * 8);
            blk_ptr.av1xd->mb_to_right_edge = ((pcs_ptr->av1_cm->mi_cols - bh - micol) * MI_SIZE) *
                8;
            // Perform final pass using the 1/8 MV
            //AV1 MVs are always in 1/8th pel precision.
            mv_unit.mv->x = context_ptr->tf_16x16_mv_x[idx_32x32 * 4 + idx_16x16];
            mv_unit.mv->y = context_ptr->tf_16x16_mv_y[idx_32x32 * 4 + idx_16x16];
            av1_inter_prediction(scs_ptr,
                                 NULL, //pcs_ptr,
                                 (uint32_t)interp_filters,
                                 &blk_ptr,
                                 0, //ref_frame_type,
                                 &mv_unit,
                                 0, //use_intrabc,
                                 SIMPLE_TRANSLATION,
                                 0,
                                 0,
                                 1, //compound_idx not used
                                 NULL, // interinter_comp not used
                                 NULL,
                                 NULL,
                                 NULL,
                                 0,
                                 0,
                                 0,
                                 0,
                                 pu_origin_x,
                                 pu_origin_y,
                                 bsize,
                                 bsize,
                                 !is_highbd ? pic_ptr_ref : &reference_ptr,
                                 NULL, //ref_pic_list1,
                                 &prediction_ptr,
                                 local_origin_x,
                                 local_origin_y,
                                 context_ptr->tf_chroma ? PICTURE_BUFFER_DESC_FULL_MASK
                                                        : PICTURE_BUFFER_DESC_LUMA_MASK,
                                 (uint8_t)encoder_bit_depth,
                                 0); // is_16bit_pipeline
        }
    } else {
        uint32_t bsize = 32;

        uint32_t idx_x = idx_32x32 & 0x1;
        uint32_t idx_y = idx_32x32 >> 1;

        uint16_t local_origin_x = idx_x * bsize;
        uint16_t local_origin_y = idx_y * bsize;
        uint16_t pu_origin_x    = sb_origin_x + local_origin_x;
        uint16_t pu_origin_y    = sb_origin_y + local_origin_y;
        int32_t mirow = pu_origin_y >> MI_SIZE_LOG2;
        int32_t micol = pu_origin_x >> MI_SIZE_LOG2;
        blk_ptr.mds_idx = get_mds_idx(local_origin_x,
                                      local_origin_y,
                                      bsize,
                                      pcs_ptr->scs_ptr->seq_header.sb_size == BLOCK_128X128);

        const int32_t bw                 = mi_size_wide[BLOCK_32X32];
        const int32_t bh                 = mi_size_high[BLOCK_32X32];
        blk_ptr.av1xd->mb_to_top_edge    = -(int32_t)((mirow * MI_SIZE) * 8);
        blk_ptr.av1xd->mb_to_bottom_edge = ((pcs_ptr->av1_cm->mi_rows - bw - mirow) * MI_SIZE) * 8;
        blk_ptr.av1xd->mb_to_left_edge   = -(int32_t)((micol * MI_SIZE) * 8);
        blk_ptr.av1xd->mb_to_right_edge  = ((pcs_ptr->av1_cm->mi_cols - bh - micol) * MI_SIZE) * 8;

        // Perform final pass using the 1/8 MV
        // AV1 MVs are always in 1/8th pel precision.
        mv_unit.mv->x = context_ptr->tf_32x32_mv_x[idx_32x32];
        mv_unit.mv->y = context_ptr->tf_32x32_mv_y[idx_32x32];

        av1_inter_prediction(
            scs_ptr,
            NULL, //pcs_ptr,
            (uint32_t)interp_filters,
            &blk_ptr,
            0, //ref_frame_type,
            &mv_unit,
            0, //use_intrabc,
            SIMPLE_TRANSLATION,
            0,
            0,
            1, //compound_idx not used
            NULL, // interinter_comp not used
            NULL,
            NULL,
            NULL,
            0,
            0,
            0,
            0,
            pu_origin_x,
            pu_origin_y,
            bsize,
            bsize,
            !is_highbd ? pic_ptr_ref : &reference_ptr,
            NULL, //ref_pic_list1,
            &prediction_ptr,
            local_origin_x,
            local_origin_y,
            context_ptr->tf_chroma ? PICTURE_BUFFER_DESC_FULL_MASK : PICTURE_BUFFER_DESC_LUMA_MASK,
            (uint8_t)encoder_bit_depth,
            0); // is_16bit_pipeline
    }
}

void get_final_filtered_pixels_c(MeContext *context_ptr, EbByte *src_center_ptr_start,
                                 uint16_t **altref_buffer_highbd_start, uint32_t **accum,
                                 uint16_t **count, const uint32_t *stride, int blk_y_src_offset,
                                 int blk_ch_src_offset, uint16_t blk_width_ch,
                                 uint16_t blk_height_ch, Bool is_highbd) {
    int i, j, k;

    if (!is_highbd) {
        // Process luma
        int pos = blk_y_src_offset;
        for (i = 0, k = 0; i < BH; i++) {
            for (j = 0; j < BW; j++, k++) {
                assert(OD_DIVU(
                    accum[C_Y][k] + (count[C_Y][k] >> 1), count[C_Y][k]) < 256);
                src_center_ptr_start[C_Y][pos] = (uint8_t)OD_DIVU(
                    accum[C_Y][k] + (count[C_Y][k] >> 1), count[C_Y][k]);
                pos++;
            }
            pos += stride[C_Y] - BW;
        }
        // Process chroma
        if (context_ptr->tf_chroma) {
            pos = blk_ch_src_offset;
            for (i = 0, k = 0; i < blk_height_ch; i++) {
                for (j = 0; j < blk_width_ch; j++, k++) {
                    assert(OD_DIVU(
                        accum[C_U][k] + (count[C_U][k] >> 1), count[C_U][k]) < 256);
                    src_center_ptr_start[C_U][pos] = (uint8_t)OD_DIVU(
                        accum[C_U][k] + (count[C_U][k] >> 1), count[C_U][k]);
                    assert(OD_DIVU(
                        accum[C_U][k] + (count[C_U][k] >> 1), count[C_U][k]) < 256);
                    src_center_ptr_start[C_V][pos] = (uint8_t)OD_DIVU(
                        accum[C_V][k] + (count[C_V][k] >> 1), count[C_V][k]);
                    pos++;
                }
                pos += stride[C_U] - blk_width_ch;
            }
        }
    } else {
        // Process luma
        int pos = blk_y_src_offset;
        for (i = 0, k = 0; i < BH; i++) {
            for (j = 0; j < BW; j++, k++) {
                altref_buffer_highbd_start[C_Y][pos] = (uint16_t)OD_DIVU(
                    accum[C_Y][k] + (count[C_Y][k] >> 1), count[C_Y][k]);
                pos++;
            }
            pos += stride[C_Y] - BW;
        }
        // Process chroma
        if (context_ptr->tf_chroma) {
            pos = blk_ch_src_offset;
            for (i = 0, k = 0; i < blk_height_ch; i++) {
                for (j = 0; j < blk_width_ch; j++, k++) {
                    altref_buffer_highbd_start[C_U][pos] = (uint16_t)OD_DIVU(
                        accum[C_U][k] + (count[C_U][k] >> 1), count[C_U][k]);
                    altref_buffer_highbd_start[C_V][pos] = (uint16_t)OD_DIVU(
                        accum[C_V][k] + (count[C_V][k] >> 1), count[C_V][k]);
                    pos++;
                }
                pos += stride[C_U] - blk_width_ch;
            }
        }
    }
}
/*
* Check whether to perform 64x64 pred only
*/
int8_t tf_use_64x64_pred(MeContext *context_ptr) {
    uint32_t dist_32x32 = 0;

    // 32x32
    for (unsigned i = 0; i < 4; i++) { dist_32x32 += context_ptr->p_best_sad_32x32[i]; }

    int64_t dev_64x64_to_32x32 = (int64_t)(((int64_t)MAX(context_ptr->p_best_sad_64x64[0], 1) -
                                            (int64_t)MAX(dist_32x32, 1)) *
                                           100) /
        (int64_t)MAX(dist_32x32, 1);
    if (dev_64x64_to_32x32 < context_ptr->tf_use_pred_64x64_only_th)
        return 1;
    else
        return 0;
}
// Produce the filtered alt-ref picture
// - core function
static EbErrorType produce_temporally_filtered_pic(
    PictureParentControlSet **list_picture_control_set_ptr,
    EbPictureBufferDesc **list_input_picture_ptr, uint8_t index_center,
    MotionEstimationContext_t *me_context_ptr,
    const int32_t *noise_levels_log1p_fp16,
    const double *noise_levels, int32_t segment_index, Bool is_highbd) {
    DECLARE_ALIGNED(16, uint32_t, accumulator[BLK_PELS * COLOR_CHANNELS]);
    DECLARE_ALIGNED(16, uint16_t, counter[BLK_PELS * COLOR_CHANNELS]);
    uint32_t *accum[COLOR_CHANNELS] = {
        accumulator, accumulator + BLK_PELS, accumulator + (BLK_PELS << 1)};
    uint16_t *count[COLOR_CHANNELS] = {counter, counter + BLK_PELS, counter + (BLK_PELS << 1)};

    EbByte    predictor       = {NULL};
    uint16_t *predictor_16bit = {NULL};
    PictureParentControlSet *picture_control_set_ptr_central =
        list_picture_control_set_ptr[index_center];
    SequenceControlSet *scs = picture_control_set_ptr_central->scs_ptr;
    EbPictureBufferDesc *input_picture_ptr_central = list_input_picture_ptr[index_center];
    MeContext *          context_ptr               = me_context_ptr->me_context_ptr;

    // Prep 8bit source if 8bit content or using 8bit for subpel
    if (!is_highbd || context_ptr->tf_ctrls.use_8bit_subpel)
        EB_MALLOC_ALIGNED_ARRAY(predictor, BLK_PELS * COLOR_CHANNELS);

    if (is_highbd)
        EB_MALLOC_ALIGNED_ARRAY(predictor_16bit, BLK_PELS * COLOR_CHANNELS);
    EbByte    pred[COLOR_CHANNELS] = {predictor, predictor + BLK_PELS, predictor + (BLK_PELS << 1)};
    uint16_t *pred_16bit[COLOR_CHANNELS] = {
        predictor_16bit, predictor_16bit + BLK_PELS, predictor_16bit + (BLK_PELS << 1)};
    int encoder_bit_depth = scs->static_config.encoder_bit_depth;

    // chroma subsampling
    uint32_t ss_x = scs->subsampling_x;
    uint32_t ss_y = scs->subsampling_y;
    uint16_t blk_width_ch  = (uint16_t)BW >> ss_x;
    uint16_t blk_height_ch = (uint16_t)BH >> ss_y;

    uint32_t blk_cols = (uint32_t)(input_picture_ptr_central->width + BW - 1) /
        BW; // I think only the part of the picture
    uint32_t blk_rows = (uint32_t)(input_picture_ptr_central->height + BH - 1) /
        BH; // that fits to the 32x32 blocks are actually filtered

    uint32_t stride[COLOR_CHANNELS]      = {input_picture_ptr_central->stride_y,
                                       input_picture_ptr_central->stride_cb,
                                       input_picture_ptr_central->stride_cr};
    uint32_t stride_pred[COLOR_CHANNELS] = {BW, blk_width_ch, blk_width_ch};
    uint32_t x_seg_idx;
    uint32_t y_seg_idx;
    uint32_t picture_width_in_b64  = blk_cols;
    uint32_t picture_height_in_b64 = blk_rows;
    SEGMENT_CONVERT_IDX_TO_XY(segment_index,
                              x_seg_idx,
                              y_seg_idx,
                              picture_control_set_ptr_central->tf_segments_column_count);
    uint32_t x_b64_start_idx = SEGMENT_START_IDX(
        x_seg_idx, picture_width_in_b64, picture_control_set_ptr_central->tf_segments_column_count);
    uint32_t x_b64_end_idx = SEGMENT_END_IDX(
        x_seg_idx, picture_width_in_b64, picture_control_set_ptr_central->tf_segments_column_count);
    uint32_t y_b64_start_idx = SEGMENT_START_IDX(
        y_seg_idx, picture_height_in_b64, picture_control_set_ptr_central->tf_segments_row_count);
    uint32_t y_b64_end_idx = SEGMENT_END_IDX(
        y_seg_idx, picture_height_in_b64, picture_control_set_ptr_central->tf_segments_row_count);

    // first position of the frame buffer according to the index center
    EbByte src_center_ptr_start[COLOR_CHANNELS] = {
        input_picture_ptr_central->buffer_y +
            input_picture_ptr_central->origin_y * input_picture_ptr_central->stride_y +
            input_picture_ptr_central->origin_x,
        input_picture_ptr_central->buffer_cb +
            (input_picture_ptr_central->origin_y >> ss_y) * input_picture_ptr_central->stride_cb +
            (input_picture_ptr_central->origin_x >> ss_x),
        input_picture_ptr_central->buffer_cr +
            (input_picture_ptr_central->origin_y >> ss_y) * input_picture_ptr_central->stride_cr +
            (input_picture_ptr_central->origin_x >> ss_x),
    };

    uint16_t *altref_buffer_highbd_start[COLOR_CHANNELS] = {
        picture_control_set_ptr_central->altref_buffer_highbd[C_Y] +
            input_picture_ptr_central->origin_y * input_picture_ptr_central->stride_y +
            input_picture_ptr_central->origin_x,
        picture_control_set_ptr_central->altref_buffer_highbd[C_U] +
            (input_picture_ptr_central->origin_y >> ss_y) *
                input_picture_ptr_central->stride_bit_inc_cb +
            (input_picture_ptr_central->origin_x >> ss_x),
        picture_control_set_ptr_central->altref_buffer_highbd[C_V] +
            (input_picture_ptr_central->origin_y >> ss_y) *
                input_picture_ptr_central->stride_bit_inc_cr +
            (input_picture_ptr_central->origin_x >> ss_x),
    };
    int decay_control;
    if (scs->vq_ctrls.sharpness_ctrls.tf && picture_control_set_ptr_central->is_noise_level && scs->calculate_variance && picture_control_set_ptr_central->pic_avg_variance < VQ_PIC_AVG_VARIANCE_TH) {
        decay_control = 1;
    }
    else {
        // Hyper-parameter for filter weight adjustment.
        decay_control = (context_ptr->tf_ctrls.use_fast_filter) ? 5
            : (context_ptr->tf_ctrls.use_medium_filter) ? 3
            : 4;
        // Decrease the filter strength for low QPs
        if (scs->static_config.qp <= ALT_REF_QP_THRESH)
            decay_control--;
    }
    // Adjust filtering based on q.
    // Larger q -> stronger filtering -> larger weight.
    // Smaller q -> weaker filtering -> smaller weight.

    // Fixed-QP offsets are use here since final picture QP(s) are not generated @ this early stage
    const int bit_depth = scs->static_config.encoder_bit_depth;
    int       active_best_quality = 0;
    int       active_worst_quality =
        quantizer_to_qindex[(uint8_t)scs->static_config.qp];
    int q;
    if (context_ptr->tf_ctrls.use_fixed_point || context_ptr->tf_ctrls.use_medium_filter) {
        FP_ASSERT(TF_FILTER_STRENGTH == 5);
        FP_ASSERT(TF_STRENGTH_THRESHOLD == 4);
        FP_ASSERT(TF_Q_DECAY_THRESHOLD == 20);

        int offset_idx;
        if (!picture_control_set_ptr_central->is_used_as_reference_flag)
            offset_idx = -1;
        else if (picture_control_set_ptr_central->idr_flag)
            offset_idx = 0;
        else
            offset_idx = MIN(picture_control_set_ptr_central->temporal_layer_index + 1,
                                FIXED_QP_OFFSET_COUNT - 1);

        // Fixed-QP offsets are use here since final picture QP(s) are not generated @ this early stage
        int32_t q_val_fp8 = svt_av1_convert_qindex_to_q_fp8(active_worst_quality, bit_depth);

        const int32_t q_val_target_fp8 = (offset_idx == -1)
            ? q_val_fp8
            : MAX(q_val_fp8 -
                        (q_val_fp8 *
                        percents[picture_control_set_ptr_central->hierarchical_levels <= 4]
                                [offset_idx] /
                        100),
                    0);

        const int32_t delta_qindex_f = svt_av1_compute_qdelta_fp(
            q_val_fp8, q_val_target_fp8, bit_depth);
        active_best_quality = (int32_t)(active_worst_quality + delta_qindex_f);
        q = active_best_quality;

        FP_ASSERT(q < (1 << 20));
        //double q_decay = pow((double)q / TF_Q_DECAY_THRESHOLD, 2);

        //TODO: FIX1 uint32_t???
        int32_t  q_treshold_fp8 = (q * ((int32_t)1 << 8)) / TF_Q_DECAY_THRESHOLD;
        uint32_t q_decay_fp8 = (q_treshold_fp8 * q_treshold_fp8) >> 8;

        //  q_decay_fp8 = ((q_decay_fp8) < (1) ? (1) : (q_decay_fp8) > (1<<8) ? (1<<8) : (q_decay_fp8));
        q_decay_fp8 = CLIP(q_decay_fp8, 1, 1 << 8);
        if (q >= TF_QINDEX_CUTOFF) {
            // Max q_factor is 255, therefore the upper bound of q_decay is 8.
            // We do not need a clip here.
            //q_decay = 0.5 * pow((double)q / 64, 2);
            FP_ASSERT(q < (1 << 15));
            q_decay_fp8 = (q * q) >> 5;
        }

        const int32_t const_0dot7_fp16 = 45875; //0.7
        /*Calculation of log and dceay_factor possible to move to estimate_noise() and calculate one time for GOP*/
        //decay_control * (0.7 + log1p(noise_levels[C_Y]))
        int32_t n_decay_fp10 = (decay_control * (const_0dot7_fp16 + noise_levels_log1p_fp16[C_Y])) /
            ((int32_t)1 << 6);
        //2 * n_decay * n_decay * q_decay * (s_decay always is 1);
        context_ptr->tf_decay_factor_fp16[C_Y] = (uint32_t)(
            (((((int64_t)n_decay_fp10) * ((int64_t)n_decay_fp10))) * q_decay_fp8) >> 11);

        if (context_ptr->tf_chroma) {
            n_decay_fp10 = (decay_control * (const_0dot7_fp16 + noise_levels_log1p_fp16[C_U])) /
                ((int32_t)1 << 6);
            context_ptr->tf_decay_factor_fp16[C_U] = (uint32_t)(
                (((((int64_t)n_decay_fp10) * ((int64_t)n_decay_fp10))) * q_decay_fp8) >> 11);

            n_decay_fp10 = (decay_control * (const_0dot7_fp16 + noise_levels_log1p_fp16[C_V])) /
                ((int32_t)1 << 6);
            context_ptr->tf_decay_factor_fp16[C_V] = (uint32_t)(
                (((((int64_t)n_decay_fp10) * ((int64_t)n_decay_fp10))) * q_decay_fp8) >> 11);
        }
    }
    else {
        double q_val = svt_av1_convert_qindex_to_q(active_worst_quality, bit_depth);
        int    offset_idx = -1;
        if (!picture_control_set_ptr_central->is_used_as_reference_flag)
            offset_idx = -1;
        else if (picture_control_set_ptr_central->idr_flag)
            offset_idx = 0;
        else
            offset_idx = MIN(picture_control_set_ptr_central->temporal_layer_index + 1,
                FIXED_QP_OFFSET_COUNT - 1);

        const double q_val_target = (offset_idx == -1)
            ? q_val
            : MAX(q_val -
                (q_val *
                    percents[picture_control_set_ptr_central->hierarchical_levels <= 4]
                    [offset_idx] /
                    100),
                0.0);

        const int32_t delta_qindex = svt_av1_compute_qdelta(q_val, q_val_target, bit_depth);

        active_best_quality = (int32_t)(active_worst_quality + delta_qindex);
        q = active_best_quality;
        clamp(q, active_best_quality, active_worst_quality);

        double q_decay = pow((double)q / TF_Q_DECAY_THRESHOLD, 2);
        q_decay = CLIP(q_decay, 1e-5, 1);
        if (q >= TF_QINDEX_CUTOFF) {
            // Max q_factor is 255, therefore the upper bound of q_decay is 8.
            // We do not need a clip here.
            q_decay = 0.5 * pow((double)q / 64, 2);
        }

        // Smaller strength -> smaller filtering weight.
        double s_decay = pow((double)TF_FILTER_STRENGTH / TF_STRENGTH_THRESHOLD, 2);
        s_decay = CLIP(s_decay, 1e-5, 1);
        double n_decay;
        n_decay = (double)decay_control * (0.7 + log1p(noise_levels[0]));
        context_ptr->tf_decay_factor[C_Y] = 2 * n_decay * n_decay * q_decay * s_decay;

        if (context_ptr->tf_chroma) {
            n_decay = (double)decay_control * (0.7 + log1p(noise_levels[1]));
            context_ptr->tf_decay_factor[C_U] = 2 * n_decay * n_decay * q_decay * s_decay;

            n_decay = (double)decay_control * (0.7 + log1p(noise_levels[2]));
            context_ptr->tf_decay_factor[C_V] = 2 * n_decay * n_decay * q_decay * s_decay;
        }
    }
    for (uint32_t blk_row = y_b64_start_idx; blk_row < y_b64_end_idx; blk_row++) {
        for (uint32_t blk_col = x_b64_start_idx; blk_col < x_b64_end_idx; blk_col++) {
            int blk_y_src_offset  = (blk_col * BW) + (blk_row * BH) * stride[C_Y];
            int blk_ch_src_offset = (blk_col * blk_width_ch) +
                (blk_row * blk_height_ch) * stride[C_U];

            // reset accumulator and count
            memset(accumulator, 0, BLK_PELS * COLOR_CHANNELS * sizeof(accumulator[0]));
            memset(counter, 0, BLK_PELS * COLOR_CHANNELS * sizeof(counter[0]));

            EbByte    src_center_ptr[COLOR_CHANNELS]           = {NULL};
            uint16_t *altref_buffer_highbd_ptr[COLOR_CHANNELS] = {NULL};
            // Prep 8bit source if 8bit content or using 8bit for subpel
            if (!is_highbd || context_ptr->tf_ctrls.use_8bit_subpel) {
                src_center_ptr[C_Y] = src_center_ptr_start[C_Y] + blk_y_src_offset;
                if (context_ptr->tf_chroma) {
                    src_center_ptr[C_U] = src_center_ptr_start[C_U] + blk_ch_src_offset;
                    src_center_ptr[C_V] = src_center_ptr_start[C_V] + blk_ch_src_offset;
                }
            }

            if (is_highbd) {
                altref_buffer_highbd_ptr[C_Y] = altref_buffer_highbd_start[C_Y] + blk_y_src_offset;
                if (context_ptr->tf_chroma) {
                    altref_buffer_highbd_ptr[C_U] = altref_buffer_highbd_start[C_U] +
                        blk_ch_src_offset;
                    altref_buffer_highbd_ptr[C_V] = altref_buffer_highbd_start[C_V] +
                        blk_ch_src_offset;
                }
            }

            if (!is_highbd)
                apply_filtering_central(context_ptr,
                                        input_picture_ptr_central,
                                        src_center_ptr,
                                        accum,
                                        count,
                                        BW,
                                        BH,
                                        ss_x,
                                        ss_y);
            else
                apply_filtering_central_highbd(context_ptr,
                                               input_picture_ptr_central,
                                               altref_buffer_highbd_ptr,
                                               accum,
                                               count,
                                               BW,
                                               BH,
                                               ss_x,
                                               ss_y);

            // for every frame to filter
            for (int frame_index = 0;
                 frame_index < (picture_control_set_ptr_central->past_altref_nframes +
                                picture_control_set_ptr_central->future_altref_nframes + 1);
                 frame_index++) {
                // ------------
                // Step 1: motion estimation + compensation
                // ------------
                me_context_ptr->me_context_ptr->tf_frame_index  = frame_index;
                me_context_ptr->me_context_ptr->tf_index_center = index_center;
                // if frame to process is the center frame
                if (frame_index == index_center) {
                    // skip MC (central frame)
                    if (!is_highbd) {
                        pic_copy_kernel_8bit(
                            src_center_ptr[C_Y], stride[C_Y], pred[C_Y], stride_pred[C_Y], BW, BH);
                        if (context_ptr->tf_chroma) {
                            pic_copy_kernel_8bit(src_center_ptr[C_U],
                                                 stride[C_U],
                                                 pred[C_U],
                                                 stride_pred[C_U],
                                                 blk_width_ch,
                                                 blk_height_ch);
                            pic_copy_kernel_8bit(src_center_ptr[C_V],
                                                 stride[C_V],
                                                 pred[C_V],
                                                 stride_pred[C_V],
                                                 blk_width_ch,
                                                 blk_height_ch);
                        }
                    } else {
                        pic_copy_kernel_16bit(altref_buffer_highbd_ptr[C_Y],
                                              stride[C_Y],
                                              pred_16bit[C_Y],
                                              stride_pred[C_Y],
                                              BW,
                                              BH);
                        if (context_ptr->tf_chroma) {
                            pic_copy_kernel_16bit(altref_buffer_highbd_ptr[C_U],
                                                  stride[C_U],
                                                  pred_16bit[C_U],
                                                  stride_pred[C_U],
                                                  blk_width_ch,
                                                  blk_height_ch);
                            pic_copy_kernel_16bit(altref_buffer_highbd_ptr[C_V],
                                                  stride[C_V],
                                                  pred_16bit[C_V],
                                                  stride_pred[C_V],
                                                  blk_width_ch,
                                                  blk_height_ch);
                        }
                    }

                } else {
                    // Initialize ME context
                    create_me_context_and_picture_control(
                        me_context_ptr,
                        list_picture_control_set_ptr[frame_index],
                        list_picture_control_set_ptr[index_center],
                        input_picture_ptr_central,
                        blk_row,
                        blk_col,
                        ss_x,
                        ss_y);
                    context_ptr->num_of_list_to_search = 1;
                    context_ptr->num_of_ref_pic_to_search[0] = 1;
                    context_ptr->num_of_ref_pic_to_search[1] = 0;
                    context_ptr->temporal_layer_index =
                        picture_control_set_ptr_central->temporal_layer_index;
                    context_ptr->is_used_as_reference_flag =
                        picture_control_set_ptr_central->is_used_as_reference_flag;

                    EbPaReferenceObject *reference_object = (EbPaReferenceObject *)
                                                                context_ptr->alt_ref_reference_ptr;
                    context_ptr->me_ds_ref_array[0][0].picture_ptr =
                        reference_object->input_padded_picture_ptr;
                    context_ptr->me_ds_ref_array[0][0].sixteenth_picture_ptr =
                        reference_object->sixteenth_downsampled_picture_ptr;
                    context_ptr->me_ds_ref_array[0][0].quarter_picture_ptr =
                        reference_object->quarter_downsampled_picture_ptr;
                    context_ptr->me_ds_ref_array[0][0].picture_number =
                        reference_object->picture_number;
                    context_ptr->tf_me_exit_th =
                        picture_control_set_ptr_central->tf_ctrls.me_exit_th;
                    ;
                    context_ptr->tf_use_pred_64x64_only_th =
                        picture_control_set_ptr_central->tf_ctrls.use_pred_64x64_only_th;
                    context_ptr->tf_subpel_early_exit =
                        picture_control_set_ptr_central->tf_ctrls.subpel_early_exit;
                    // Perform ME - context_ptr will store the outputs (MVs, buffers, etc)
                    // Block-based MC using open-loop HME + refinement
                    motion_estimation_b64(picture_control_set_ptr_central,
                        (uint32_t)blk_row * blk_cols + blk_col,
                        (uint32_t)blk_col * BW, // x block
                        (uint32_t)blk_row * BH, // y block
                        context_ptr,
                        input_picture_ptr_central); // source picture

                    if (context_ptr->tf_use_pred_64x64_only_th &&
                        (context_ptr->tf_use_pred_64x64_only_th == (uint8_t)~0 ||
                         tf_use_64x64_pred(context_ptr))) {
                        tf_64x64_sub_pel_search(
                            picture_control_set_ptr_central,
                            context_ptr,
                            list_picture_control_set_ptr[frame_index],
                            list_input_picture_ptr[frame_index],
                            pred,
                            pred_16bit,
                            stride_pred,
                            src_center_ptr,
                            altref_buffer_highbd_ptr,
                            stride,
                            (uint32_t)blk_col * BW,
                            (uint32_t)blk_row * BH,
                            ss_x,
                            (context_ptr->tf_ctrls.use_8bit_subpel) ? EB_EIGHT_BIT : encoder_bit_depth);

                        // Perform MC using the information acquired using the ME step
                        tf_64x64_inter_prediction(picture_control_set_ptr_central,
                                                  context_ptr,
                                                  list_picture_control_set_ptr[frame_index],
                                                  list_input_picture_ptr[frame_index],
                                                  pred,
                                                  pred_16bit,
                                                  (uint32_t)blk_col * BW,
                                                  (uint32_t)blk_row * BH,
                                                  ss_x,
                                                  encoder_bit_depth);

                    } else {
                        // 64x64 Sub-Pel search
                        tf_64x64_sub_pel_search(
                            picture_control_set_ptr_central,
                            context_ptr,
                            list_picture_control_set_ptr[frame_index],
                            list_input_picture_ptr[frame_index],
                            pred,
                            pred_16bit,
                            stride_pred,
                            src_center_ptr,
                            altref_buffer_highbd_ptr,
                            stride,
                            (uint32_t)blk_col* BW,
                            (uint32_t)blk_row* BH,
                            ss_x,
                            (context_ptr->tf_ctrls.use_8bit_subpel) ? EB_EIGHT_BIT : encoder_bit_depth);


                        // 32x32 Sub-Pel search
                        for (int block_row = 0; block_row < 2; block_row++) {
                            for (int block_col = 0; block_col < 2; block_col++) {

                                context_ptr->idx_32x32 = block_col + (block_row << 1);

                                tf_32x32_sub_pel_search(picture_control_set_ptr_central,
                                    context_ptr,
                                    list_picture_control_set_ptr[frame_index],
                                    list_input_picture_ptr[frame_index],
                                    pred,
                                    pred_16bit,
                                    stride_pred,
                                    src_center_ptr,
                                    altref_buffer_highbd_ptr,
                                    stride,
                                    (uint32_t)blk_col * BW,
                                    (uint32_t)blk_row * BH,
                                    ss_x,
                                    (context_ptr->tf_ctrls.use_8bit_subpel)
                                    ? EB_EIGHT_BIT
                                    : encoder_bit_depth);
                            }
                        }


                        uint64_t sum_32x32_block_error =
                            context_ptr->tf_32x32_block_error[0] +
                            context_ptr->tf_32x32_block_error[1] +
                            context_ptr->tf_32x32_block_error[2] +
                            context_ptr->tf_32x32_block_error[3];

                        if ((context_ptr->tf_64x64_block_error * 14 < sum_32x32_block_error * 16) &&
                            context_ptr->tf_64x64_block_error < (1 << 18)) {

                            tf_64x64_inter_prediction(picture_control_set_ptr_central,
                                context_ptr,
                                list_picture_control_set_ptr[frame_index],
                                list_input_picture_ptr[frame_index],
                                pred,
                                pred_16bit,
                                (uint32_t)blk_col * BW,
                                (uint32_t)blk_row * BH,
                                ss_x,
                                encoder_bit_depth);

                            context_ptr->tf_32x32_mv_x[0] = context_ptr->tf_64x64_mv_x;
                            context_ptr->tf_32x32_mv_y[0] = context_ptr->tf_64x64_mv_y;

                            context_ptr->tf_32x32_mv_x[1] = context_ptr->tf_64x64_mv_x;
                            context_ptr->tf_32x32_mv_y[1] = context_ptr->tf_64x64_mv_y;

                            context_ptr->tf_32x32_mv_x[2] = context_ptr->tf_64x64_mv_x;
                            context_ptr->tf_32x32_mv_y[2] = context_ptr->tf_64x64_mv_y;

                            context_ptr->tf_32x32_mv_x[3] = context_ptr->tf_64x64_mv_x;
                            context_ptr->tf_32x32_mv_y[3] = context_ptr->tf_64x64_mv_y;

                            context_ptr->tf_32x32_block_split_flag[0] = 0;
                            context_ptr->tf_32x32_block_split_flag[1] = 0;
                            context_ptr->tf_32x32_block_split_flag[2] = 0;
                            context_ptr->tf_32x32_block_split_flag[3] = 0;

                            // Update the 32x32 block-error
                            for (int block_row = 0; block_row < 2; block_row++) {
                                for (int block_col = 0; block_col < 2; block_col++) {

                                    uint32_t bsize = 32;
                                    uint64_t distortion;
                                    context_ptr->idx_32x32 = block_col + (block_row << 1);

                                    if (!is_highbd) {
                                        uint8_t* pred_y_ptr = pred[C_Y] + bsize * block_row * stride_pred[C_Y] +  bsize * block_col;
                                        uint8_t* src_y_ptr = src_center_ptr[C_Y] + bsize * block_row * stride[C_Y] + bsize * block_col;

                                        const AomVarianceFnPtr* fn_ptr = picture_control_set_ptr_central->tf_ctrls.sub_sampling_shift ? &mefn_ptr[BLOCK_32X16]
                                            : &mefn_ptr[BLOCK_32X32];
                                        unsigned int            sse;
                                        distortion = fn_ptr->vf(pred_y_ptr,
                                            stride_pred[C_Y] << picture_control_set_ptr_central->tf_ctrls.sub_sampling_shift,
                                            src_y_ptr,
                                            stride[C_Y] << picture_control_set_ptr_central->tf_ctrls.sub_sampling_shift,
                                            &sse)
                                            << picture_control_set_ptr_central->tf_ctrls.sub_sampling_shift;
                                    }
                                    else {
                                        uint16_t* pred_y_ptr = pred_16bit[C_Y] +
                                            bsize * block_row * stride_pred[C_Y] +
                                            bsize * block_col;
                                        uint16_t* src_y_ptr = altref_buffer_highbd_ptr[C_Y] +
                                            bsize * block_row * stride[C_Y] +
                                            bsize * block_col;
                                        const AomVarianceFnPtr* fn_ptr = picture_control_set_ptr_central->tf_ctrls.sub_sampling_shift ? &mefn_ptr[BLOCK_32X16]
                                            : &mefn_ptr[BLOCK_32X32];

                                        unsigned int sse;

                                        distortion = fn_ptr->vf_hbd_10(CONVERT_TO_BYTEPTR(pred_y_ptr),
                                            stride_pred[C_Y] << picture_control_set_ptr_central->tf_ctrls.sub_sampling_shift,
                                            CONVERT_TO_BYTEPTR(src_y_ptr),
                                            stride[C_Y] << picture_control_set_ptr_central->tf_ctrls.sub_sampling_shift,
                                            &sse)
                                            << picture_control_set_ptr_central->tf_ctrls.sub_sampling_shift;
                                    }
                                    context_ptr->tf_32x32_block_error[context_ptr->idx_32x32] = distortion;
                                }
                            }

                        }
                        else {
                            // 16x16 Sub-Pel search, and 32x32 partitioning
                            for (int block_row = 0; block_row < 2; block_row++) {
                                for (int block_col = 0; block_col < 2; block_col++) {

                                    context_ptr->idx_32x32 = block_col + (block_row << 1);
                                    if (context_ptr->tf_32x32_block_error[context_ptr->idx_32x32] < picture_control_set_ptr_central->tf_ctrls.pred_error_32x32_th) {
                                        context_ptr->tf_32x32_block_split_flag[context_ptr->idx_32x32] =
                                            0;
                                    } else {
                                        tf_16x16_sub_pel_search(picture_control_set_ptr_central,
                                            context_ptr,
                                            list_picture_control_set_ptr[frame_index],
                                            list_input_picture_ptr[frame_index],
                                            pred,
                                            pred_16bit,
                                            stride_pred,
                                            src_center_ptr,
                                            altref_buffer_highbd_ptr,
                                            stride,
                                            (uint32_t)blk_col * BW,
                                            (uint32_t)blk_row * BH,
                                            ss_x,
                                            (context_ptr->tf_ctrls.use_8bit_subpel)
                                            ? EB_EIGHT_BIT
                                            : encoder_bit_depth);

                                        // Derive tf_32x32_block_split_flag
                                        derive_tf_32x32_block_split_flag(context_ptr);
                                    }
                                        // Perform MC using the information acquired using the ME step
                                        tf_32x32_inter_prediction(picture_control_set_ptr_central,
                                            context_ptr,
                                            list_picture_control_set_ptr[frame_index],
                                            list_input_picture_ptr[frame_index],
                                            pred,
                                            pred_16bit,
                                            (uint32_t)blk_col * BW,
                                            (uint32_t)blk_row * BH,
                                            ss_x,
                                            encoder_bit_depth);
                                }
                            }
                        }
                    }

                    for (int block_row = 0; block_row < 2; block_row++) {
                        for (int block_col = 0; block_col < 2; block_col++) {
                            context_ptr->tf_block_col = block_col;
                            context_ptr->tf_block_row = block_row;
                            apply_filtering_block_plane_wise(context_ptr,
                                                             block_row,
                                                             block_col,
                                                             src_center_ptr,
                                                             altref_buffer_highbd_ptr,
                                                             pred,
                                                             pred_16bit,
                                                             accum,
                                                             count,
                                                             stride,
                                                             stride_pred,
                                                             BW >> 1,
                                                             BH >> 1,
                                                             ss_x,
                                                             ss_y,
                                                             encoder_bit_depth);
                        }
                    }
                }
            }

            // Normalize filter output to produce temporally filtered frame
            get_final_filtered_pixels(context_ptr,
                                      src_center_ptr_start,
                                      altref_buffer_highbd_start,
                                      accum,
                                      count,
                                      stride,
                                      blk_y_src_offset,
                                      blk_ch_src_offset,
                                      blk_width_ch,
                                      blk_height_ch,
                                      is_highbd);
        }
    }
    // Prep 8bit source if 8bit content or using 8bit for subpel
    if (!is_highbd || context_ptr->tf_ctrls.use_8bit_subpel)
        EB_FREE_ALIGNED_ARRAY(predictor);

    if (is_highbd)
        EB_FREE_ALIGNED_ARRAY(predictor_16bit);
    return EB_ErrorNone;
}

// This is an adaptation of the mehtod in the following paper:
// Shen-Chuan Tai, Shih-Ming Yang, "A fast method for image noise
// estimation using Laplacian operator and adaptive edge detection,"
// Proc. 3rd International Symposium on Communications, Control and
// Signal Processing, 2008, St Julians, Malta.
// Return noise estimate, or -1.0 if there was a failure
// function from libaom
// Standard bit depht input (=8 bits) to estimate the noise, I don't think there needs to be two methods for this
// Operates on the Y component only
int32_t svt_estimate_noise_fp16_c(const uint8_t *src, uint16_t width, uint16_t height, uint16_t stride_y) {
    int64_t sum = 0;
    int64_t num = 0;

    for (int i = 1; i < height - 1; ++i) {
        for (int j = 1; j < width - 1; ++j) {
            const int k = i * stride_y + j;
            // Sobel gradients
            const int g_x = (src[k - stride_y - 1] - src[k - stride_y + 1]) +
                (src[k + stride_y - 1] - src[k + stride_y + 1]) + 2 * (src[k - 1] - src[k + 1]);
            const int g_y = (src[k - stride_y - 1] - src[k + stride_y - 1]) +
                (src[k - stride_y + 1] - src[k + stride_y + 1]) +
                2 * (src[k - stride_y] - src[k + stride_y]);
            const int ga = abs(g_x) + abs(g_y);
            if (ga < EDGE_THRESHOLD) { // Do not consider edge pixels to estimate the noise
                // Find Laplacian
                const int v = 4 * src[k] -
                    2 * (src[k - 1] + src[k + 1] + src[k - stride_y] + src[k + stride_y]) +
                    (src[k - stride_y - 1] + src[k - stride_y + 1] + src[k + stride_y - 1] +
                     src[k + stride_y + 1]);
                sum += abs(v);
                ++num;
            }
        }
    }
    // If very few smooth pels, return -1 since the estimate is unreliable
    if (num < SMOOTH_THRESHOLD) {
        return -65536 /*-1:fp16*/;
    }

    FP_ASSERT((((int64_t)sum * SQRT_PI_BY_2_FP16) / (6 * num)) < ((int64_t)1 << 31));
    return (int32_t)((sum * SQRT_PI_BY_2_FP16) / (6 * num));
}

// Noise estimation for highbd
int32_t svt_estimate_noise_highbd_fp16_c(const uint16_t *src, int width, int height, int stride, int bd) {
    int64_t sum = 0;
    int64_t num = 0;

    for (int i = 1; i < height - 1; ++i) {
        for (int j = 1; j < width - 1; ++j) {
            const int k = i * stride + j;
            // Sobel gradients
            const int g_x = (src[k - stride - 1] - src[k - stride + 1]) +
                (src[k + stride - 1] - src[k + stride + 1]) + 2 * (src[k - 1] - src[k + 1]);
            const int g_y = (src[k - stride - 1] - src[k + stride - 1]) +
                (src[k - stride + 1] - src[k + stride + 1]) +
                2 * (src[k - stride] - src[k + stride]);
            const int ga = ROUND_POWER_OF_TWO(abs(g_x) + abs(g_y),
                                              bd - 8); // divide by 2^2 and round up
            if (ga < EDGE_THRESHOLD) { // Do not consider edge pixels to estimate the noise
                // Find Laplacian
                const int v = 4 * src[k] -
                    2 * (src[k - 1] + src[k + 1] + src[k - stride] + src[k + stride]) +
                    (src[k - stride - 1] + src[k - stride + 1] + src[k + stride - 1] +
                     src[k + stride + 1]);
                sum += ROUND_POWER_OF_TWO(abs(v), bd - 8);
                ++num;
            }
        }
    }
    // If very few smooth pels, return -1 since the estimate is unreliable
    if (num < SMOOTH_THRESHOLD) {
        return -65536 /*-1:fp16*/;
    }

    FP_ASSERT((((int64_t)sum * SQRT_PI_BY_2_FP16) / (6 * num)) < ((int64_t)1 << 31));
    return (int32_t)((sum * SQRT_PI_BY_2_FP16) / (6 * num));
}

double svt_estimate_noise_c(const uint8_t *src, uint16_t width, uint16_t height, uint16_t stride_y) {
    int64_t sum = 0;
    int64_t num = 0;

    for (int i = 1; i < height - 1; ++i) {
        for (int j = 1; j < width - 1; ++j) {
            const int k = i * stride_y + j;
            // Sobel gradients
            const int g_x = (src[k - stride_y - 1] - src[k - stride_y + 1]) +
                (src[k + stride_y - 1] - src[k + stride_y + 1]) + 2 * (src[k - 1] - src[k + 1]);
            const int g_y = (src[k - stride_y - 1] - src[k + stride_y - 1]) +
                (src[k - stride_y + 1] - src[k + stride_y + 1]) +
                2 * (src[k - stride_y] - src[k + stride_y]);
            const int ga = abs(g_x) + abs(g_y);
            if (ga < EDGE_THRESHOLD) { // Do not consider edge pixels to estimate the noise
                // Find Laplacian
                const int v = 4 * src[k] -
                    2 * (src[k - 1] + src[k + 1] + src[k - stride_y] + src[k + stride_y]) +
                    (src[k - stride_y - 1] + src[k - stride_y + 1] + src[k + stride_y - 1] +
                     src[k + stride_y + 1]);
                sum += abs(v);
                ++num;
            }
        }
    }
    // If very few smooth pels, return -1 since the estimate is unreliable
    if (num < SMOOTH_THRESHOLD)
        return -1.0;

    const double sigma = (double)sum / (6 * num) * SQRT_PI_BY_2;

    return sigma;
}

// Noise estimation for highbd
double svt_estimate_noise_highbd_c(const uint16_t *src, int width, int height, int stride, int bd) {
    int64_t sum = 0;
    int64_t num = 0;

    for (int i = 1; i < height - 1; ++i) {
        for (int j = 1; j < width - 1; ++j) {
            const int k = i * stride + j;
            // Sobel gradients
            const int g_x = (src[k - stride - 1] - src[k - stride + 1]) +
                (src[k + stride - 1] - src[k + stride + 1]) + 2 * (src[k - 1] - src[k + 1]);
            const int g_y = (src[k - stride - 1] - src[k + stride - 1]) +
                (src[k - stride + 1] - src[k + stride + 1]) +
                2 * (src[k - stride] - src[k + stride]);
            const int ga = ROUND_POWER_OF_TWO(abs(g_x) + abs(g_y),
                                              bd - 8); // divide by 2^2 and round up
            if (ga < EDGE_THRESHOLD) { // Do not consider edge pixels to estimate the noise
                // Find Laplacian
                const int v = 4 * src[k] -
                    2 * (src[k - 1] + src[k + 1] + src[k - stride] + src[k + stride]) +
                    (src[k - stride - 1] + src[k - stride + 1] + src[k + stride - 1] +
                     src[k + stride + 1]);
                sum += ROUND_POWER_OF_TWO(abs(v), bd - 8);
                ++num;
            }
        }
    }
    // If very few smooth pels, return -1 since the estimate is unreliable
    if (num < SMOOTH_THRESHOLD)
        return -1.0;

    const double sigma = (double)sum / (6 * num) * SQRT_PI_BY_2;
    return sigma;
}

void pad_and_decimate_filtered_pic(PictureParentControlSet *picture_control_set_ptr_central) {
    // reference structures (padded pictures + downsampled versions)
    SequenceControlSet *scs_ptr = picture_control_set_ptr_central->scs_ptr;
    EbPaReferenceObject *src_object = (EbPaReferenceObject *)picture_control_set_ptr_central
                                          ->pa_reference_picture_wrapper_ptr->object_ptr;
    EbPictureBufferDesc *input_picture_ptr = picture_control_set_ptr_central->enhanced_picture_ptr;

    // Refine the non-8 padding
    if (((input_picture_ptr->width - scs_ptr->pad_right) % 8 != 0) ||
        ((input_picture_ptr->height - scs_ptr->pad_bottom) % 8 != 0))
        pad_picture_to_multiple_of_min_blk_size_dimensions(scs_ptr, input_picture_ptr);

    //Generate padding first, then copy
    generate_padding(&(input_picture_ptr->buffer_y[C_Y]),
                     input_picture_ptr->stride_y,
                     input_picture_ptr->width,
                     input_picture_ptr->height,
                     input_picture_ptr->origin_x,
                     input_picture_ptr->origin_y);
    // Padding chroma after altref
    if (picture_control_set_ptr_central->tf_ctrls.do_chroma) {
        generate_padding(input_picture_ptr->buffer_cb,
                         input_picture_ptr->stride_cb,
                         input_picture_ptr->width >> scs_ptr->subsampling_x,
                         input_picture_ptr->height >> scs_ptr->subsampling_y,
                         input_picture_ptr->origin_x >> scs_ptr->subsampling_x,
                         input_picture_ptr->origin_y >> scs_ptr->subsampling_y);
        generate_padding(input_picture_ptr->buffer_cr,
                         input_picture_ptr->stride_cr,
                         input_picture_ptr->width >> scs_ptr->subsampling_x,
                         input_picture_ptr->height >> scs_ptr->subsampling_y,
                         input_picture_ptr->origin_x >> scs_ptr->subsampling_x,
                         input_picture_ptr->origin_y >> scs_ptr->subsampling_y);
    }

    // 1/4 & 1/16 input picture downsampling
    if (scs_ptr->down_sampling_method_me_search == ME_FILTERED_DOWNSAMPLED) {
        downsample_filtering_input_picture(picture_control_set_ptr_central,
                                           input_picture_ptr,
                                           src_object->quarter_downsampled_picture_ptr,
                                           src_object->sixteenth_downsampled_picture_ptr);
    } else {
        downsample_decimation_input_picture(picture_control_set_ptr_central,
                                            input_picture_ptr,
                                            src_object->quarter_downsampled_picture_ptr,
                                            src_object->sixteenth_downsampled_picture_ptr);
    }
}

// save original enchanced_picture_ptr buffer in a separate buffer (to be replaced by the temporally filtered pic)
static EbErrorType save_src_pic_buffers(PictureParentControlSet *picture_control_set_ptr_central,
                                        uint32_t ss_y, Bool is_highbd) {
    // save buffer from full size frame enhanced_unscaled_picture_ptr
    EbPictureBufferDesc *src_pic_ptr =
        picture_control_set_ptr_central->enhanced_unscaled_picture_ptr;
    assert(src_pic_ptr != NULL);
    // allocate memory for the copy of the original enhanced buffer
    EB_MALLOC_ARRAY(picture_control_set_ptr_central->save_source_picture_ptr[C_Y],
                    src_pic_ptr->luma_size);
    EB_MALLOC_ARRAY(picture_control_set_ptr_central->save_source_picture_ptr[C_U],
                    src_pic_ptr->chroma_size);
    EB_MALLOC_ARRAY(picture_control_set_ptr_central->save_source_picture_ptr[C_V],
                    src_pic_ptr->chroma_size);

    // if highbd, allocate memory for the copy of the original enhanced buffer - bit inc
    if (is_highbd) {
        EB_MALLOC_ARRAY(picture_control_set_ptr_central->save_source_picture_bit_inc_ptr[C_Y],
                        src_pic_ptr->luma_size);
        EB_MALLOC_ARRAY(picture_control_set_ptr_central->save_source_picture_bit_inc_ptr[C_U],
                        src_pic_ptr->chroma_size);
        EB_MALLOC_ARRAY(picture_control_set_ptr_central->save_source_picture_bit_inc_ptr[C_V],
                        src_pic_ptr->chroma_size);
    }
    picture_control_set_ptr_central->save_source_picture_width  = src_pic_ptr->width;
    picture_control_set_ptr_central->save_source_picture_height = src_pic_ptr->height;

    // copy buffers
    // Y
    uint32_t height_y  = (uint32_t)(src_pic_ptr->height + src_pic_ptr->origin_y +
                                   src_pic_ptr->origin_bot_y);
    uint32_t height_uv = (uint32_t)(
        (src_pic_ptr->height + src_pic_ptr->origin_y + src_pic_ptr->origin_bot_y) >> ss_y);

    assert(height_y * src_pic_ptr->stride_y == src_pic_ptr->luma_size);
    assert(height_uv * src_pic_ptr->stride_cb == src_pic_ptr->chroma_size);
    assert(height_uv * src_pic_ptr->stride_cr == src_pic_ptr->chroma_size);

    pic_copy_kernel_8bit(src_pic_ptr->buffer_y,
                         src_pic_ptr->stride_y,
                         picture_control_set_ptr_central->save_source_picture_ptr[C_Y],
                         src_pic_ptr->stride_y,
                         src_pic_ptr->stride_y,
                         height_y);

    pic_copy_kernel_8bit(src_pic_ptr->buffer_cb,
                         src_pic_ptr->stride_cb,
                         picture_control_set_ptr_central->save_source_picture_ptr[C_U],
                         src_pic_ptr->stride_cb,
                         src_pic_ptr->stride_cb,
                         height_uv);

    pic_copy_kernel_8bit(src_pic_ptr->buffer_cr,
                         src_pic_ptr->stride_cr,
                         picture_control_set_ptr_central->save_source_picture_ptr[C_V],
                         src_pic_ptr->stride_cr,
                         src_pic_ptr->stride_cr,
                         height_uv);

    if (is_highbd) {
        // if highbd, copy bit inc buffers
        // Y
        svt_c_unpack_compressed_10bit(
            picture_control_set_ptr_central->enhanced_picture_ptr->buffer_bit_inc_y,
            picture_control_set_ptr_central->enhanced_picture_ptr->stride_bit_inc_y / 4,
            picture_control_set_ptr_central->save_source_picture_bit_inc_ptr[C_Y],
            picture_control_set_ptr_central->enhanced_picture_ptr->stride_bit_inc_y,
            height_y);
        // U
        svt_c_unpack_compressed_10bit(
            picture_control_set_ptr_central->enhanced_picture_ptr->buffer_bit_inc_cb,
            picture_control_set_ptr_central->enhanced_picture_ptr->stride_bit_inc_cb / 4,
            picture_control_set_ptr_central->save_source_picture_bit_inc_ptr[C_U],
            picture_control_set_ptr_central->enhanced_picture_ptr->stride_bit_inc_cb,
            height_uv);
        // V
        svt_c_unpack_compressed_10bit(
            picture_control_set_ptr_central->enhanced_picture_ptr->buffer_bit_inc_cr,
            picture_control_set_ptr_central->enhanced_picture_ptr->stride_bit_inc_cr / 4,
            picture_control_set_ptr_central->save_source_picture_bit_inc_ptr[C_V],
            picture_control_set_ptr_central->enhanced_picture_ptr->stride_bit_inc_cr,
            height_uv);
    }

    return EB_ErrorNone;
}

EbErrorType svt_av1_init_temporal_filtering(
    PictureParentControlSet ** list_picture_control_set_ptr,
    PictureParentControlSet *  picture_control_set_ptr_central,
    MotionEstimationContext_t *me_context_ptr, int32_t segment_index) {
    uint8_t              index_center;
    EbPictureBufferDesc *central_picture_ptr;

    me_context_ptr->me_context_ptr->tf_chroma = picture_control_set_ptr_central->tf_ctrls.do_chroma;

    me_context_ptr->me_context_ptr->tf_ctrls = picture_control_set_ptr_central->tf_ctrls;

    me_context_ptr->me_context_ptr->tf_tot_horz_blks =
        me_context_ptr->me_context_ptr->tf_tot_vert_blks = 0;
    // index of the central source frame
    index_center = picture_control_set_ptr_central->past_altref_nframes;

    // if this assertion does not fail (as I think it should not, then remove picture_control_set_ptr_central from the input parameters of init_temporal_filtering())
    assert(list_picture_control_set_ptr[index_center] == picture_control_set_ptr_central);

    // source central frame picture buffer
    central_picture_ptr = picture_control_set_ptr_central->enhanced_picture_ptr;

    uint32_t encoder_bit_depth =
        picture_control_set_ptr_central->scs_ptr->static_config.encoder_bit_depth;
    Bool is_highbd = (encoder_bit_depth == 8) ? (uint8_t)FALSE : (uint8_t)TRUE;

    // chroma subsampling
    uint32_t ss_x = picture_control_set_ptr_central->scs_ptr->subsampling_x;
    uint32_t ss_y = picture_control_set_ptr_central->scs_ptr->subsampling_y;
    int32_t *noise_levels_log1p_fp16 = &(
        picture_control_set_ptr_central->noise_levels_log1p_fp16[0]);
    double *noise_levels = &(picture_control_set_ptr_central->noise_levels[0]);

    //only one performs any picture based prep
    svt_block_on_mutex(picture_control_set_ptr_central->temp_filt_mutex);
    if (picture_control_set_ptr_central->temp_filt_prep_done == 0) {
        picture_control_set_ptr_central->temp_filt_prep_done = 1;
        // Pad chroma reference samples - once only per picture
        for (int i = 0; i < (picture_control_set_ptr_central->past_altref_nframes +
                             picture_control_set_ptr_central->future_altref_nframes + 1);
             i++) {
            EbPictureBufferDesc *pic_ptr_ref =
                list_picture_control_set_ptr[i]->enhanced_picture_ptr;
            //10bit: for all the reference pictures do the packing once at the beggining.
            if (is_highbd && i != picture_control_set_ptr_central->past_altref_nframes) {
                EB_MALLOC_ARRAY(list_picture_control_set_ptr[i]->altref_buffer_highbd[C_Y],
                                central_picture_ptr->luma_size);
                if (me_context_ptr->me_context_ptr->tf_chroma) {
                    EB_MALLOC_ARRAY(list_picture_control_set_ptr[i]->altref_buffer_highbd[C_U],
                                    central_picture_ptr->chroma_size);
                    EB_MALLOC_ARRAY(list_picture_control_set_ptr[i]->altref_buffer_highbd[C_V],
                                    central_picture_ptr->chroma_size);
                }
                // pack byte buffers to 16 bit buffer
                pack_highbd_pic(pic_ptr_ref,
                                list_picture_control_set_ptr[i]->altref_buffer_highbd,
                                ss_x,
                                ss_y,
                                TRUE);
            }
        }

        picture_control_set_ptr_central->temporal_filtering_on =
            TRUE; // set temporal filtering flag ON for current picture

        // save original source picture (to be replaced by the temporally filtered pic)
        // if stat_report is enabled for PSNR computation
        // or if superres recode is enabled
        SUPERRES_MODE superres_mode =
            picture_control_set_ptr_central->scs_ptr->static_config.superres_mode;
        SUPERRES_AUTO_SEARCH_TYPE search_type =
            picture_control_set_ptr_central->scs_ptr->static_config.superres_auto_search_type;
        uint32_t frame_update_type = get_frame_update_type(picture_control_set_ptr_central->scs_ptr,
                                                           picture_control_set_ptr_central);
        Bool   superres_recode_enabled = (superres_mode == SUPERRES_AUTO) &&
            ((search_type == SUPERRES_AUTO_DUAL) ||
             (search_type == SUPERRES_AUTO_ALL)) // auto-dual or auto-all
            && ((frame_update_type == KF_UPDATE) ||
                (frame_update_type == ARF_UPDATE)); // recode only applies to key and arf
        if (picture_control_set_ptr_central->scs_ptr->static_config.stat_report ||
            superres_recode_enabled) {
            save_src_pic_buffers(picture_control_set_ptr_central, ss_y, is_highbd);
        }
    }
    svt_release_mutex(picture_control_set_ptr_central->temp_filt_mutex);
    me_context_ptr->me_context_ptr->min_frame_size = MIN(
        picture_control_set_ptr_central->aligned_height,
        picture_control_set_ptr_central->aligned_width);
    // index of the central source frame
    // index_center = picture_control_set_ptr_central->past_altref_nframes;
    // populate source frames picture buffer list
    EbPictureBufferDesc *list_input_picture_ptr[ALTREF_MAX_NFRAMES] = {NULL};
    for (int i = 0; i < (picture_control_set_ptr_central->past_altref_nframes +
                         picture_control_set_ptr_central->future_altref_nframes + 1);
         i++)
        list_input_picture_ptr[i] = list_picture_control_set_ptr[i]->enhanced_unscaled_picture_ptr;

    produce_temporally_filtered_pic(list_picture_control_set_ptr,
                                    list_input_picture_ptr,
                                    index_center,
                                    me_context_ptr,
                                    noise_levels_log1p_fp16,
                                    noise_levels,
                                    segment_index,
                                    is_highbd);

    svt_block_on_mutex(picture_control_set_ptr_central->temp_filt_mutex);
    picture_control_set_ptr_central->temp_filt_seg_acc++;

    picture_control_set_ptr_central->tf_tot_horz_blks +=
        me_context_ptr->me_context_ptr->tf_tot_horz_blks;
    picture_control_set_ptr_central->tf_tot_vert_blks +=
        me_context_ptr->me_context_ptr->tf_tot_vert_blks;

    if (picture_control_set_ptr_central->temp_filt_seg_acc ==
        picture_control_set_ptr_central->tf_segments_total_count) {
#if DEBUG_TF
        if (!is_highbd)
            save_YUV_to_file("filtered_picture.yuv",
                             central_picture_ptr->buffer_y,
                             central_picture_ptr->buffer_cb,
                             central_picture_ptr->buffer_cr,
                             central_picture_ptr->width,
                             central_picture_ptr->height,
                             central_picture_ptr->stride_y,
                             central_picture_ptr->stride_cb,
                             central_picture_ptr->stride_cr,
                             central_picture_ptr->origin_y,
                             central_picture_ptr->origin_x,
                             ss_x,
                             ss_y);
        else
            save_YUV_to_file_highbd("filtered_picture.yuv",
                                    picture_control_set_ptr_central->altref_buffer_highbd[C_Y],
                                    picture_control_set_ptr_central->altref_buffer_highbd[C_U],
                                    picture_control_set_ptr_central->altref_buffer_highbd[C_V],
                                    central_picture_ptr->width,
                                    central_picture_ptr->height,
                                    central_picture_ptr->stride_y,
                                    central_picture_ptr->stride_cb,
                                    central_picture_ptr->stride_cb,
                                    central_picture_ptr->origin_y,
                                    central_picture_ptr->origin_x,
                                    ss_x,
                                    ss_y);
#endif
        if (is_highbd) {
            unpack_highbd_pic(picture_control_set_ptr_central->altref_buffer_highbd,
                              central_picture_ptr,
                              ss_x,
                              ss_y,
                              TRUE);
            EB_FREE_ARRAY(picture_control_set_ptr_central->altref_buffer_highbd[C_Y]);
            if (me_context_ptr->me_context_ptr->tf_chroma) {
                EB_FREE_ARRAY(picture_control_set_ptr_central->altref_buffer_highbd[C_U]);
                EB_FREE_ARRAY(picture_control_set_ptr_central->altref_buffer_highbd[C_V]);
            }
            for (int i = 0; i < (picture_control_set_ptr_central->past_altref_nframes +
                                 picture_control_set_ptr_central->future_altref_nframes + 1);
                 i++) {
                if (i != picture_control_set_ptr_central->past_altref_nframes) {
                    EB_FREE_ARRAY(list_picture_control_set_ptr[i]->altref_buffer_highbd[C_Y]);
                    if (me_context_ptr->me_context_ptr->tf_chroma) {
                        EB_FREE_ARRAY(list_picture_control_set_ptr[i]->altref_buffer_highbd[C_U]);
                        EB_FREE_ARRAY(list_picture_control_set_ptr[i]->altref_buffer_highbd[C_V]);
                    }
                }
            }
        }

        // padding + decimation: even if highbd src, this is only performed on the 8 bit buffer (excluding the LSBs)
        pad_and_decimate_filtered_pic(picture_control_set_ptr_central);

        // signal that temp filt is done
        svt_post_semaphore(picture_control_set_ptr_central->temp_filt_done_semaphore);
    }

    svt_release_mutex(picture_control_set_ptr_central->temp_filt_mutex);

    return EB_ErrorNone;
}
// clang-format on
