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

typedef void (*av1_inv_txfm_highbd_func)(const int32_t *coeff, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, int32_t bd);
typedef void (*av1_inv_txfm2d_highbd_rect_func)(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);
av1_inv_txfm_highbd_func av1_inv_txfm_highbd_func_ptr_array_base[3] = { eb_av1_inv_txfm2d_add_16x16_avx2 , eb_av1_inv_txfm2d_add_32x32_avx2 , eb_av1_inv_txfm2d_add_64x64_sse4_1 };
av1_inv_txfm_highbd_func av1_inv_txfm_highbd_func_ptr_array_opt[3] = { eb_av1_inv_txfm2d_add_16x16_avx512, eb_av1_inv_txfm2d_add_32x32_avx512 , eb_av1_inv_txfm2d_add_64x64_avx512 };
av1_inv_txfm2d_highbd_rect_func av1_inv_txfm_highbd_rect_func_ptr_array_base[6] = { eb_av1_inv_txfm2d_add_32x16_c , eb_av1_inv_txfm2d_add_16x32_c , eb_av1_inv_txfm2d_add_16x64_c , eb_av1_inv_txfm2d_add_32x64_c , eb_av1_inv_txfm2d_add_64x32_c , eb_av1_inv_txfm2d_add_64x16_c };
av1_inv_txfm2d_highbd_rect_func av1_inv_txfm_highbd_rect_func_ptr_array_opt[6] = { eb_av1_inv_txfm2d_add_32x16_avx512 , eb_av1_inv_txfm2d_add_16x32_avx512 , eb_av1_inv_txfm2d_add_16x64_avx512 , eb_av1_inv_txfm2d_add_32x64_avx512 , eb_av1_inv_txfm2d_add_64x32_avx512 , eb_av1_inv_txfm2d_add_64x16_avx512 };
int txsize_16[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8 ,9, 10, 11, 12, 13, 14, 15};
int rect_types[] = { 10 , 9 , 17 , 11 , 12 , 18};
int txsize_32[] = { 0 , 9};
int txsize_64[] = { 0};
int bd[] = { 10, 12};

static void init_data(int32_t **input, int32_t **input_opt, int32_t *input_stride) {
    *input_stride = eb_create_random_aligned_stride(MAX_SB_SIZE, 64);
    TEST_ALLIGN_MALLOC(int32_t*, *input, sizeof(int32_t) * MAX_SB_SIZE * *input_stride);
    TEST_ALLIGN_MALLOC(int32_t*, *input_opt, sizeof(int32_t) * MAX_SB_SIZE * *input_stride);
    memset(*input, 0, MAX_SB_SIZE * *input_stride);
    memset(*input_opt, 0, MAX_SB_SIZE * *input_stride);
    eb_buf_random_s32(*input, MAX_SB_SIZE * *input_stride);
    memcpy(*input_opt, *input, sizeof(**input) * MAX_SB_SIZE * *input_stride);
}

static void init_data_with_max(int32_t **input, int32_t **input_opt, int32_t *input_stride) {
    *input_stride = eb_create_random_aligned_stride(MAX_SB_SIZE, 64);
    TEST_ALLIGN_MALLOC(int32_t*, *input, sizeof(int32_t) * MAX_SB_SIZE * *input_stride);
    TEST_ALLIGN_MALLOC(int32_t*, *input_opt, sizeof(int32_t) * MAX_SB_SIZE * *input_stride);
    memset(*input, 0, MAX_SB_SIZE * *input_stride);
    memset(*input_opt, 0, MAX_SB_SIZE * *input_stride);
    eb_buf_random_s32_with_max(*input, MAX_SB_SIZE * *input_stride, 1023);
    memcpy(*input_opt, *input, sizeof(**input) * MAX_SB_SIZE * *input_stride);
}

static void uninit_data(int32_t *coeff, int32_t *coeff_opt, int32_t *stride) {
    (void)stride;
    TEST_ALLIGN_FREE(coeff);
    TEST_ALLIGN_FREE(coeff_opt);
}

static void uninit_output(uint16_t *output, uint16_t *output_opt) {
    TEST_ALLIGN_FREE(output);
    TEST_ALLIGN_FREE(output_opt);
}

static void init_output_r(uint16_t **output_r, uint16_t **output_opt_r, int32_t num) {
    TEST_ALLIGN_MALLOC(uint16_t*, *output_r, sizeof(uint16_t) * MAX_SB_SIZE * num);
    TEST_ALLIGN_MALLOC(uint16_t*, *output_opt_r, sizeof(uint16_t) * MAX_SB_SIZE * num);
    eb_buf_random_u16(*output_r, MAX_SB_SIZE * num);
    memcpy(*output_opt_r, *output_r, sizeof(uint16_t) * MAX_SB_SIZE * num);
}

static void init_output_w(uint16_t **output_w, uint16_t **output_opt_w, int32_t num) {
    TEST_ALLIGN_MALLOC(uint16_t*, *output_w, sizeof(uint16_t) * MAX_SB_SIZE * num);
    TEST_ALLIGN_MALLOC(uint16_t*, *output_opt_w, sizeof(uint16_t) * MAX_SB_SIZE * num);
    memset(*output_w, 0, sizeof(uint16_t) *MAX_SB_SIZE * num);
    memset(*output_opt_w, 0, sizeof(uint16_t) *MAX_SB_SIZE * num);
}

TEST(InverseTransformTest, av1_inv_txfm_2d_square_kernels)
{
    int32_t* coeff,* coeff_opt;
    uint16_t *output_r, *output_opt_r, *output_w, *output_opt_w;
    int32_t stride;
    for (int loop = 0; loop < 3; loop++) {          //Function Pairs
        for (int i = 0; i < 10; i++) {              //Number of Test Runs
            for (int x = 0; x < 2; x++) {           //Bit Depth
                switch (loop) {
                case 0://16x16
                    for (int j = 0; j < 16; j++) {
                        init_data(&coeff, &coeff_opt, &stride);
                        ASSERT(eb_buf_compare_s32(coeff, coeff_opt, MAX_SB_SIZE * stride) == 1);
                        init_output_r(&output_r, &output_opt_r, stride);
                        init_output_w(&output_w, &output_opt_w, stride);
                        av1_inv_txfm_highbd_func_ptr_array_base[loop](coeff, output_r, stride, output_w, stride, (TxType)txsize_16[j], bd[x]);
                        av1_inv_txfm_highbd_func_ptr_array_opt[loop](coeff, output_opt_r, stride, output_opt_w, stride, (TxType)txsize_16[j], bd[x]);
                        EXPECT_EQ(eb_buf_compare_u16(output_w, output_opt_w, MAX_SB_SIZE * stride), 1);
                        uninit_output(output_r, output_opt_r);
                        uninit_output(output_w, output_opt_w);
                        uninit_data(coeff, coeff_opt, &stride);
                    }
                break;
                case 1://32x32
                    for (int j = 0; j < 2; j++) {
                        init_data(&coeff, &coeff_opt, &stride);
                        ASSERT(eb_buf_compare_s32(coeff, coeff_opt, MAX_SB_SIZE * stride) == 1);
                        init_output_r(&output_r, &output_opt_r, stride);
                        init_output_w(&output_w, &output_opt_w, stride);
                        av1_inv_txfm_highbd_func_ptr_array_base[loop](coeff, output_r, stride, output_w, stride, (TxType)txsize_32[j], bd[x]);
                        av1_inv_txfm_highbd_func_ptr_array_opt[loop](coeff, output_opt_r, stride, output_opt_w, stride, (TxType)txsize_32[j], bd[x]);
                        EXPECT_EQ(eb_buf_compare_u16(output_w, output_opt_w, MAX_SB_SIZE * stride), 1);
                        uninit_output(output_r, output_opt_r);
                        uninit_output(output_w, output_opt_w);
                        uninit_data(coeff, coeff_opt, &stride);
                    }
                break;
                case 2://64x64
                    for (int j = 0; j < 1; j++) {
                        init_data(&coeff, &coeff_opt, &stride);
                        ASSERT(eb_buf_compare_s32(coeff, coeff_opt, MAX_SB_SIZE * stride) == 1);
                        init_output_r(&output_r, &output_opt_r, stride);
                        init_output_w(&output_w, &output_opt_w, stride);
                        av1_inv_txfm_highbd_func_ptr_array_base[loop](coeff, output_r, stride, output_w, stride, (TxType)txsize_64[j], bd[x]);
                        av1_inv_txfm_highbd_func_ptr_array_opt[loop](coeff, output_opt_r, stride, output_opt_w, stride, (TxType)txsize_64[j], bd[x]);
                        EXPECT_EQ(eb_buf_compare_u16(output_w, output_opt_w, MAX_SB_SIZE * stride), 1);
                        uninit_output(output_r, output_opt_r);
                        uninit_output(output_w, output_opt_w);
                        uninit_data(coeff, coeff_opt, &stride);
                    }
                break;

                default:
                    ASSERT(0);
                }
            }
        }
    }
}

TEST(InverseTransformTest, av1_inv_txfm_2d_rect_kernels)
{
    int32_t* coeff, *coeff_opt;
    uint16_t *output_r, *output_opt_r, *output_w, *output_opt_w;
    int32_t stride;
    for (int loop = 0; loop < 6; loop++) {          //Function Pairs
        for (int i = 0; i < 10; i++) {              //Number of Test Runs
            for (int x = 0; x < 2; x++) {           //Bit Depth
                switch (loop) {
                case 0://32x16
                    for (int j = 0; j < 2; j++) {
                        init_data_with_max(&coeff, &coeff_opt, &stride);
                        ASSERT(eb_buf_compare_s32(coeff, coeff_opt, MAX_SB_SIZE * stride) == 1);
                        init_output_r(&output_r, &output_opt_r, stride);
                        init_output_w(&output_w, &output_opt_w, stride);
                        av1_inv_txfm_highbd_rect_func_ptr_array_base[loop](coeff, output_r, stride, output_w, stride, (TxType)txsize_32[j], (TxSize)rect_types[loop], 0, bd[x]);
                        av1_inv_txfm_highbd_rect_func_ptr_array_opt[loop](coeff, output_opt_r, stride, output_opt_w, stride, (TxType)txsize_32[j], (TxSize)rect_types[loop], 0, bd[x]);
                        EXPECT_EQ(eb_buf_compare_u16(output_w, output_opt_w, MAX_SB_SIZE * stride), 1);
                        uninit_output(output_r, output_opt_r);
                        uninit_output(output_w, output_opt_w);
                        uninit_data(coeff, coeff_opt, &stride);
                    }
                    break;
                case 1://16x32
                    for (int j = 0; j < 2; j++) {
                        init_data_with_max(&coeff, &coeff_opt, &stride);
                        ASSERT(eb_buf_compare_s32(coeff, coeff_opt, MAX_SB_SIZE * stride) == 1);
                        init_output_r(&output_r, &output_opt_r, stride);
                        init_output_w(&output_w, &output_opt_w, stride);
                        av1_inv_txfm_highbd_rect_func_ptr_array_base[loop](coeff, output_r, stride, output_w, stride, (TxType)txsize_32[j], (TxSize)rect_types[loop], 0, bd[x]);
                        av1_inv_txfm_highbd_rect_func_ptr_array_opt[loop](coeff, output_opt_r, stride, output_opt_w, stride, (TxType)txsize_32[j], (TxSize)rect_types[loop], 0, bd[x]);
                        EXPECT_EQ(eb_buf_compare_u16(output_w, output_opt_w, MAX_SB_SIZE * stride), 1);
                        uninit_output(output_r, output_opt_r);
                        uninit_output(output_w, output_opt_w);
                        uninit_data(coeff, coeff_opt, &stride);
                    }
                    break;
                case 2://16x64
                    for (int j = 0; j < 1; j++) {
                        init_data_with_max(&coeff, &coeff_opt, &stride);
                        ASSERT(eb_buf_compare_s32(coeff, coeff_opt, MAX_SB_SIZE * stride) == 1);
                        init_output_r(&output_r, &output_opt_r, stride);
                        init_output_w(&output_w, &output_opt_w, stride);
                        av1_inv_txfm_highbd_rect_func_ptr_array_base[loop](coeff, output_r, stride, output_w, stride, (TxType)txsize_64[j], (TxSize)rect_types[loop], 0, bd[x]);
                        av1_inv_txfm_highbd_rect_func_ptr_array_opt[loop](coeff, output_opt_r, stride, output_opt_w, stride, (TxType)txsize_64[j], (TxSize)rect_types[loop], 0, bd[x]);
                        EXPECT_EQ(eb_buf_compare_u16(output_w, output_opt_w, MAX_SB_SIZE * stride), 1);
                        uninit_output(output_r, output_opt_r);
                        uninit_output(output_w, output_opt_w);
                        uninit_data(coeff, coeff_opt, &stride);
                    }
                    break;
                case 3://32x64
                    for (int j = 0; j < 1; j++) {
                        init_data_with_max(&coeff, &coeff_opt, &stride);
                        ASSERT(eb_buf_compare_s32(coeff, coeff_opt, MAX_SB_SIZE * stride) == 1);
                        init_output_r(&output_r, &output_opt_r, stride);
                        init_output_w(&output_w, &output_opt_w, stride);
                        av1_inv_txfm_highbd_rect_func_ptr_array_base[loop](coeff, output_r, stride, output_w, stride, (TxType)txsize_64[j], (TxSize)rect_types[loop], 0, bd[x]);
                        av1_inv_txfm_highbd_rect_func_ptr_array_opt[loop](coeff, output_opt_r, stride, output_opt_w, stride, (TxType)txsize_64[j], (TxSize)rect_types[loop], 0, bd[x]);
                        EXPECT_EQ(eb_buf_compare_u16(output_w, output_opt_w, MAX_SB_SIZE * stride), 1);
                        uninit_output(output_r, output_opt_r);
                        uninit_output(output_w, output_opt_w);
                        uninit_data(coeff, coeff_opt, &stride);
                    }
                    break;
                case 4://64x32
                    for (int j = 0; j < 1; j++) {
                        init_data_with_max(&coeff, &coeff_opt, &stride);
                        ASSERT(eb_buf_compare_s32(coeff, coeff_opt, MAX_SB_SIZE * stride) == 1);
                        init_output_r(&output_r, &output_opt_r, stride);
                        init_output_w(&output_w, &output_opt_w, stride);
                        av1_inv_txfm_highbd_rect_func_ptr_array_base[loop](coeff, output_r, stride, output_w, stride, (TxType)txsize_64[j], (TxSize)rect_types[loop], 0, bd[x]);
                        av1_inv_txfm_highbd_rect_func_ptr_array_opt[loop](coeff, output_opt_r, stride, output_opt_w, stride, (TxType)txsize_64[j], (TxSize)rect_types[loop], 0, bd[x]);
                        EXPECT_EQ(eb_buf_compare_u16(output_w, output_opt_w, MAX_SB_SIZE * stride), 1);
                        uninit_output(output_r, output_opt_r);
                        uninit_output(output_w, output_opt_w);
                        uninit_data(coeff, coeff_opt, &stride);
                    }
                    break;
                case 5://64x16
                    for (int j = 0; j < 1; j++) {
                        init_data_with_max(&coeff, &coeff_opt, &stride);
                        ASSERT(eb_buf_compare_s32(coeff, coeff_opt, MAX_SB_SIZE * stride) == 1);
                        init_output_r(&output_r, &output_opt_r, stride);
                        init_output_w(&output_w, &output_opt_w, stride);
                        av1_inv_txfm_highbd_rect_func_ptr_array_base[loop](coeff, output_r, stride, output_w, stride, (TxType)txsize_64[j], (TxSize)rect_types[loop], 0, bd[x]);
                        av1_inv_txfm_highbd_rect_func_ptr_array_opt[loop](coeff, output_opt_r, stride, output_opt_w, stride, (TxType)txsize_64[j], (TxSize)rect_types[loop], 0, bd[x]);
                        EXPECT_EQ(eb_buf_compare_u16(output_w, output_opt_w, MAX_SB_SIZE * stride), 1);
                        uninit_output(output_r, output_opt_r);
                        uninit_output(output_w, output_opt_w);
                        uninit_data(coeff, coeff_opt, &stride);
                    }
                    break;

                default:
                    ASSERT(0);
                }
            }
        }
    }
}
#endif
