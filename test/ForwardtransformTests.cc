/*
 * Copyright(c) 2019 Intel Corporation
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

#include "gtest/gtest.h"
#include "EbDefinitions.h"
#include "aom_dsp_rtcd.h"
#include "EbUnitTestUtility.h"
#include "EbUnitTest.h"

#ifndef NON_AVX512_SUPPORT

typedef void(*eb_av1_frwd_txfm_func)(int16_t *input, int32_t *coeff, uint32_t stride, TxType tx_type, uint8_t bitDepth);
eb_av1_frwd_txfm_func av1_frwd_txfm_func_ptr_array_base[9] = { eb_av1_fwd_txfm2d_16x16_avx2, eb_av1_fwd_txfm2d_32x32_avx2 , eb_av1_fwd_txfm2d_64x64_avx2 , eb_av1_fwd_txfm2d_16x64_avx2, eb_av1_fwd_txfm2d_64x16_avx2 , eb_av1_fwd_txfm2d_32x64_avx2 , eb_av1_fwd_txfm2d_64x32_avx2 , eb_av1_fwd_txfm2d_16x32_avx2 , eb_av1_fwd_txfm2d_32x16_avx2 };
eb_av1_frwd_txfm_func av1_frwd_txfm_func_ptr_array_opt[9] = { av1_fwd_txfm2d_16x16_avx512, av1_fwd_txfm2d_32x32_avx512, av1_fwd_txfm2d_64x64_avx512 , av1_fwd_txfm2d_16x64_avx512, av1_fwd_txfm2d_64x16_avx512 , av1_fwd_txfm2d_32x64_avx512 , av1_fwd_txfm2d_64x32_avx512 , av1_fwd_txfm2d_16x32_avx512 , av1_fwd_txfm2d_32x16_avx512 };
int tx_16[] = { 0 , 1 , 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
int tx_32[] = { 0 , 9 };
int tx_64[] = { 0 , 9 };
int bitDepth[] = { 8, 10, 12 };

static void init_data(int16_t **input, int16_t **input_opt, uint32_t input_stride) {
    TEST_ALLIGN_MALLOC(int16_t*, *input, sizeof(int16_t) * MAX_SB_SIZE * input_stride);
    TEST_ALLIGN_MALLOC(int16_t*, *input_opt, sizeof(int16_t) * MAX_SB_SIZE * input_stride);
    memset(*input, 0, MAX_SB_SIZE * input_stride);
    memset(*input_opt, 0, MAX_SB_SIZE * input_stride);
    eb_buf_random_s16(*input, MAX_SB_SIZE * input_stride);
    memcpy(*input_opt, *input, sizeof(**input) * MAX_SB_SIZE * input_stride);
}

static void uninit_data(int16_t *input, int16_t *input_opt) {
    TEST_ALLIGN_FREE(input);
    TEST_ALLIGN_FREE(input_opt);
}

static void uninit_coeff(int32_t *coeff, int32_t *coeff_opt) {
    TEST_ALLIGN_FREE(coeff);
    TEST_ALLIGN_FREE(coeff_opt);
}

static void init_coeff(int32_t **coeff, int32_t **coeff_opt, uint32_t *stride) {
    *stride = eb_create_random_aligned_stride(MAX_SB_SIZE, 64);
    TEST_ALLIGN_MALLOC(int32_t*, *coeff, sizeof(int32_t) * MAX_SB_SIZE * *stride);
    TEST_ALLIGN_MALLOC(int32_t*, *coeff_opt, sizeof(int32_t) * MAX_SB_SIZE * *stride);
    memset(*coeff, 0, MAX_SB_SIZE * *stride);
    memset(*coeff_opt, 0, MAX_SB_SIZE * *stride);

}

void compare_s32(int32_t *output_base, int32_t *output_opt, uint32_t stride, int height, int width) {
    for (int x = 0; x < height; x++) {
        for (int y = 0; y < width; y++) {
            EXPECT_EQ(output_base[y], output_opt[y]);
        }
        output_base += stride;
        output_opt += stride;
    }
}


TEST(ForwardTransformTest, av1_frwd_txfm_kernels)
{
    int16_t *input, *input_opt;
    int32_t* coeff, *coeff_opt;
    uint32_t stride;
    init_coeff(&coeff, &coeff_opt, &stride);
    ASSERT(eb_buf_compare_s32(coeff, coeff_opt, MAX_SB_SIZE * stride) == 1);
    for (int loop = 0; loop < 9; loop++) {         //Function Pairs
        for (int i = 0; i < 10; i++) {             //Number of Test Runs
            for (int x = 0; x < 2; x++) {//Bit Depth
                switch (loop) {
                case 0://16x16
                    for (int j = 0; j < 16; j++) {
                        init_data(&input, &input_opt, stride);
                        ASSERT(eb_buf_compare_s16(input, input_opt, MAX_SB_SIZE * stride) == 1);
                        av1_frwd_txfm_func_ptr_array_base[loop](input, coeff, stride, (TxType)tx_16[j], bitDepth[x]);
                        av1_frwd_txfm_func_ptr_array_opt[loop](input_opt, coeff_opt, stride, (TxType)tx_16[j], bitDepth[x]);
                        compare_s32(coeff, coeff_opt, stride, 16, 16);
                        uninit_data(input, input_opt);
                    }
                break;
                case 1://32x32
                    for (int j = 0; j < 2; j++) {
                        init_data(&input, &input_opt, stride);
                        ASSERT(eb_buf_compare_s16(input, input_opt, MAX_SB_SIZE * stride) == 1);
                        av1_frwd_txfm_func_ptr_array_base[loop](input, coeff, stride, (TxType)tx_32[j], bitDepth[x]);
                        av1_frwd_txfm_func_ptr_array_opt[loop](input_opt, coeff_opt, stride, (TxType)tx_32[j], bitDepth[x]);
                        compare_s32(coeff, coeff_opt, stride, 32, 32);
                        uninit_data(input, input_opt);
                    }
                break;
                case 2://64x64
                    for (int j = 0; j < 1; j++) {
                        init_data(&input, &input_opt, stride);
                        ASSERT(eb_buf_compare_s16(input, input_opt, MAX_SB_SIZE * stride) == 1);
                        av1_frwd_txfm_func_ptr_array_base[loop](input, coeff, stride, (TxType)tx_64[j], bitDepth[x]);
                        av1_frwd_txfm_func_ptr_array_opt[loop](input_opt, coeff_opt, stride, (TxType)tx_64[j], bitDepth[x]);
                        compare_s32(coeff, coeff_opt, stride, 32, 32); // top - left 32x32
                        uninit_data(input, input_opt);
                    }
                break;
                case 3://16x64
                    for (int j = 0; j < 1; j++) {
                        init_data(&input, &input_opt, stride);
                        ASSERT(eb_buf_compare_s16(input, input_opt, MAX_SB_SIZE * stride) == 1);
                        av1_frwd_txfm_func_ptr_array_base[loop](input, coeff, stride, (TxType)tx_64[j], bitDepth[x]);
                        av1_frwd_txfm_func_ptr_array_opt[loop](input_opt, coeff_opt, stride, (TxType)tx_64[j], bitDepth[x]);
                        compare_s32(coeff, coeff_opt, stride, 16, 32); // top - left 16x32
                        uninit_data(input, input_opt);
                    }
                    break;
                case 4://64x16
                    for (int j = 0; j < 1; j++) {
                        init_data(&input, &input_opt, stride);
                        ASSERT(eb_buf_compare_s16(input, input_opt, MAX_SB_SIZE * stride) == 1);
                        av1_frwd_txfm_func_ptr_array_base[loop](input, coeff, stride, (TxType)tx_64[j], bitDepth[x]);
                        av1_frwd_txfm_func_ptr_array_opt[loop](input_opt, coeff_opt, stride, (TxType)tx_64[j], bitDepth[x]);
                        compare_s32(coeff, coeff_opt, stride, 32, 16); // top - left 32x16
                        uninit_data(input, input_opt);
                    }
                    break;
                case 5://32x64
                    for (int j = 0; j < 1; j++) {
                        init_data(&input, &input_opt, stride);
                        ASSERT(eb_buf_compare_s16(input, input_opt, MAX_SB_SIZE * stride) == 1);
                        av1_frwd_txfm_func_ptr_array_base[loop](input, coeff, stride, (TxType)tx_64[j], bitDepth[x]);
                        av1_frwd_txfm_func_ptr_array_opt[loop](input_opt, coeff_opt, stride, (TxType)tx_64[j], bitDepth[x]);
                        compare_s32(coeff, coeff_opt, stride, 32, 32); // top - left 32x32
                        uninit_data(input, input_opt);
                    }
                    break;
                case 6://64x32
                    for (int j = 0; j < 1; j++) {
                        init_data(&input, &input_opt, stride);
                        ASSERT(eb_buf_compare_s16(input, input_opt, MAX_SB_SIZE * stride) == 1);
                        av1_frwd_txfm_func_ptr_array_base[loop](input, coeff, stride, (TxType)tx_64[j], bitDepth[x]);
                        av1_frwd_txfm_func_ptr_array_opt[loop](input_opt, coeff_opt, stride, (TxType)tx_64[j], bitDepth[x]);
                        compare_s32(coeff, coeff_opt, stride, 32, 32); // top - left 32x32
                        uninit_data(input, input_opt);
                    }
                    break;
                case 7://16x32
                    for (int j = 0; j < 2; j++) {
                        init_data(&input, &input_opt, stride);
                        ASSERT(eb_buf_compare_s16(input, input_opt, MAX_SB_SIZE * stride) == 1);
                        av1_frwd_txfm_func_ptr_array_base[loop](input, coeff, stride, (TxType)tx_32[j], bitDepth[x]);
                        av1_frwd_txfm_func_ptr_array_opt[loop](input_opt, coeff_opt, stride, (TxType)tx_32[j], bitDepth[x]);
                        compare_s32(coeff, coeff_opt, stride, 16, 32);
                        uninit_data(input, input_opt);
                    }
                    break;
                case 8://32x16
                    for (int j = 0; j < 2; j++) {
                        init_data(&input, &input_opt, stride);
                        ASSERT(eb_buf_compare_s16(input, input_opt, MAX_SB_SIZE * stride) == 1);
                        av1_frwd_txfm_func_ptr_array_base[loop](input, coeff, stride, (TxType)tx_32[j], bitDepth[x]);
                        av1_frwd_txfm_func_ptr_array_opt[loop](input_opt, coeff_opt, stride, (TxType)tx_32[j], bitDepth[x]);
                        compare_s32(coeff, coeff_opt, stride, 32, 16);
                        uninit_data(input, input_opt);
                    }
                    break;
                default:
                    ASSERT(0);
                }
            }
        }
    }
    uninit_coeff(coeff, coeff_opt);
}
#endif
