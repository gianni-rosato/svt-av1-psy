#include "gtest/gtest.h"
#include "EbDefinitions.h"
#include "EbEncHandle.h"
#include "aom_dsp_rtcd.h"
#include "EbTransforms.h"
#include "EbUnitTestUtility.h"
#include <immintrin.h>

typedef void(*av1_frwd_txfm_func)(int16_t *input, int32_t *coeff, uint32_t stride, TxType tx_type, uint8_t bitDepth);
av1_frwd_txfm_func av1_frwd_txfm_func_ptr_array_base[9] = { av1_fwd_txfm2d_16x16_avx2, av1_fwd_txfm2d_32x32_avx2 , av1_fwd_txfm2d_64x64_avx2 , av1_fwd_txfm2d_16x64_avx2, av1_fwd_txfm2d_64x16_avx2 , av1_fwd_txfm2d_32x64_avx2 , av1_fwd_txfm2d_64x32_avx2 , av1_fwd_txfm2d_16x32_avx2 , av1_fwd_txfm2d_32x16_avx2 };
av1_frwd_txfm_func av1_frwd_txfm_func_ptr_array_opt[9] = { av1_fwd_txfm2d_16x16_avx512, av1_fwd_txfm2d_32x32_avx512, av1_fwd_txfm2d_64x64_avx512 , av1_fwd_txfm2d_16x64_avx512, av1_fwd_txfm2d_64x16_avx512 , av1_fwd_txfm2d_32x64_avx512 , av1_fwd_txfm2d_64x32_avx512 , av1_fwd_txfm2d_16x32_avx512 , av1_fwd_txfm2d_32x16_avx512 };
int tx_16[] = { 0 , 1 , 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
int tx_32[] = { 0 , 9 };
int tx_64[] = { 0 , 9 };
int bitDepth[] = { 8, 10, 12 };

static void init_data(int16_t **input, int16_t **input_opt, uint32_t input_stride) {
    *input = (int16_t*)malloc(sizeof(**input) * MAX_SB_SIZE * input_stride);
    *input_opt = (int16_t*)malloc(sizeof(**input_opt) * MAX_SB_SIZE * input_stride);
    memset(*input, 0, MAX_SB_SIZE * input_stride);
    memset(*input_opt, 0, MAX_SB_SIZE * input_stride);
    eb_buf_random_s16(*input, MAX_SB_SIZE * input_stride);
    memcpy(*input_opt, *input, sizeof(**input) * MAX_SB_SIZE * input_stride);
}

static void uninit_data(int16_t *input, int16_t *input_opt) {
    free(input);
    free(input_opt);
}

static void uninit_coeff(int32_t *coeff, int32_t *coeff_opt) {
    free(coeff);
    free(coeff_opt);
}

static void init_coeff(int32_t **coeff, int32_t **coeff_opt, uint32_t *stride) {
    *stride = eb_create_random_aligned_stride(MAX_SB_SIZE, 64);
    *coeff = (int32_t*)malloc(sizeof(**coeff) * MAX_SB_SIZE * *stride);
    *coeff_opt = (int32_t*)malloc(sizeof(**coeff_opt) * MAX_SB_SIZE * *stride);
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
            for (int x = 0; x < AOM_BITS_12; x++) {//Bit Depth
                switch (loop) {
                case 0://16x16
                    for (int j = 0; j < 16; j++) {
                        init_data(&input, &input_opt, stride);
                        ASSERT(eb_buf_compare_s16(input, input_opt, MAX_SB_SIZE * stride) == 1);
                        av1_frwd_txfm_func_ptr_array_base[loop](input,coeff,stride, (TxType)tx_16[j], bitDepth[x]);
                        av1_frwd_txfm_func_ptr_array_opt[loop](input_opt, coeff_opt, stride, (TxType)tx_16[j], bitDepth[x]);
                        EXPECT_EQ(eb_buf_compare_s32(coeff, coeff_opt, MAX_SB_SIZE * stride), 1);
                        uninit_data(input, input_opt);
                    }
                break;
                case 1://32x32
                    for (int j = 0; j < 2; j++) {
                        init_data(&input, &input_opt, stride);
                        ASSERT(eb_buf_compare_s16(input, input_opt, MAX_SB_SIZE * stride) == 1);
                        av1_frwd_txfm_func_ptr_array_base[loop](input, coeff, stride, (TxType)tx_32[j], bitDepth[x]);
                        av1_frwd_txfm_func_ptr_array_opt[loop](input_opt, coeff_opt, stride, (TxType)tx_32[j], bitDepth[x]);
                        EXPECT_EQ(eb_buf_compare_s32(coeff, coeff_opt, MAX_SB_SIZE * stride), 1);
                        uninit_data(input, input_opt);
                    }
                break;
                case 2://64x64
                    for (int j = 0; j < 1; j++) {
                        init_data(&input, &input_opt, stride);
                        ASSERT(eb_buf_compare_s16(input, input_opt, MAX_SB_SIZE * stride) == 1);
                        av1_frwd_txfm_func_ptr_array_base[loop](input, coeff, stride, (TxType)tx_64[j], bitDepth[x]);
                        av1_frwd_txfm_func_ptr_array_opt[loop](input_opt, coeff_opt, stride, (TxType)tx_64[j], bitDepth[x]);
                        EXPECT_EQ(eb_buf_compare_s32(coeff, coeff_opt, MAX_SB_SIZE * stride), 1);
                        uninit_data(input, input_opt);
                    }
                break;
                case 3://16x64
                    for (int j = 0; j < 1; j++) {
                        init_data(&input, &input_opt, stride);
                        ASSERT(eb_buf_compare_s16(input, input_opt, MAX_SB_SIZE * stride) == 1);
                        av1_frwd_txfm_func_ptr_array_base[loop](input, coeff, stride, (TxType)tx_64[j], bitDepth[x]);
                        av1_frwd_txfm_func_ptr_array_opt[loop](input_opt, coeff_opt, stride, (TxType)tx_64[j], bitDepth[x]);
                        EXPECT_EQ(eb_buf_compare_s32(coeff, coeff_opt, MAX_SB_SIZE * stride), 1);
                        uninit_data(input, input_opt);
                    }
                    break;
                case 4://64x16
                    for (int j = 0; j < 1; j++) {
                        init_data(&input, &input_opt, stride);
                        ASSERT(eb_buf_compare_s16(input, input_opt, MAX_SB_SIZE * stride) == 1);
                        av1_frwd_txfm_func_ptr_array_base[loop](input, coeff, stride, (TxType)tx_64[j], bitDepth[x]);
                        av1_frwd_txfm_func_ptr_array_opt[loop](input_opt, coeff_opt, stride, (TxType)tx_64[j], bitDepth[x]);
                        EXPECT_EQ(eb_buf_compare_s32(coeff, coeff_opt, MAX_SB_SIZE * stride), 1);
                        uninit_data(input, input_opt);
                    }
                    break;
                case 5://32x64
                    for (int j = 0; j < 1; j++) {
                        init_data(&input, &input_opt, stride);
                        ASSERT(eb_buf_compare_s16(input, input_opt, MAX_SB_SIZE * stride) == 1);
                        av1_frwd_txfm_func_ptr_array_base[loop](input, coeff, stride, (TxType)tx_64[j], bitDepth[x]);
                        av1_frwd_txfm_func_ptr_array_opt[loop](input_opt, coeff_opt, stride, (TxType)tx_64[j], bitDepth[x]);
                        EXPECT_EQ(eb_buf_compare_s32(coeff, coeff_opt, MAX_SB_SIZE * stride), 1);
                        uninit_data(input, input_opt);
                    }
                    break;
                case 6://64x32
                    for (int j = 0; j < 1; j++) {
                        init_data(&input, &input_opt, stride);
                        ASSERT(eb_buf_compare_s16(input, input_opt, MAX_SB_SIZE * stride) == 1);
                        av1_frwd_txfm_func_ptr_array_base[loop](input, coeff, stride, (TxType)tx_64[j], bitDepth[x]);
                        av1_frwd_txfm_func_ptr_array_opt[loop](input_opt, coeff_opt, stride, (TxType)tx_64[j], bitDepth[x]);
                        EXPECT_EQ(eb_buf_compare_s32(coeff, coeff_opt, MAX_SB_SIZE * stride), 1);
                        uninit_data(input, input_opt);
                    }
                    break;
                case 7://16x32
                    for (int j = 0; j < 2; j++) {
                        init_data(&input, &input_opt, stride);
                        ASSERT(eb_buf_compare_s16(input, input_opt, MAX_SB_SIZE * stride) == 1);
                        av1_frwd_txfm_func_ptr_array_base[loop](input, coeff, stride, (TxType)tx_32[j], bitDepth[x]);
                        av1_frwd_txfm_func_ptr_array_opt[loop](input_opt, coeff_opt, stride, (TxType)tx_32[j], bitDepth[x]);
                        EXPECT_EQ(eb_buf_compare_s32(coeff, coeff_opt, MAX_SB_SIZE * stride), 1);
                        uninit_data(input, input_opt);
                    }
                    break;
                case 8://32x16
                    for (int j = 0; j < 2; j++) {
                        init_data(&input, &input_opt, stride);
                        ASSERT(eb_buf_compare_s16(input, input_opt, MAX_SB_SIZE * stride) == 1);
                        av1_frwd_txfm_func_ptr_array_base[loop](input, coeff, stride, (TxType)tx_32[j], bitDepth[x]);
                        av1_frwd_txfm_func_ptr_array_opt[loop](input_opt, coeff_opt, stride, (TxType)tx_32[j], bitDepth[x]);
                        EXPECT_EQ(eb_buf_compare_s32(coeff, coeff_opt, MAX_SB_SIZE * stride), 1);
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
