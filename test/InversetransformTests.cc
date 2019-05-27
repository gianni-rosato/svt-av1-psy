#include "gtest/gtest.h"
#include "EbDefinitions.h"
#include "aom_dsp_rtcd.h"
#include "EbUnitTestUtility.h"

typedef void (*av1_inv_txfm_highbd_func)(const int32_t *coeff, uint16_t *output,int32_t stride, TxType tx_type, int32_t bd);
av1_inv_txfm_highbd_func av1_inv_txfm_highbd_func_ptr_array_base[3] = { eb_av1_inv_txfm2d_add_16x16_avx2 , eb_av1_inv_txfm2d_add_32x32_avx2 , eb_av1_inv_txfm2d_add_64x64_sse4_1 };
av1_inv_txfm_highbd_func av1_inv_txfm_highbd_func_ptr_array_opt[3] = { av1_inv_txfm2d_add_16x16_avx512, av1_inv_txfm2d_add_32x32_avx512 , av1_inv_txfm2d_add_64x64_avx512 };
int txsize_16[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8 ,9, 10, 11, 12, 13, 14, 15};
int txsize_32[] = { 0 , 9};
int txsize_64[] = { 0 };
int bd[] = {10, 12 };

static void init_data(int32_t **input, int32_t **input_opt, int32_t *input_stride) {
    *input_stride = eb_create_random_aligned_stride(MAX_SB_SIZE, 64);
    *input = (int32_t*)malloc(sizeof(**input) * MAX_SB_SIZE * *input_stride);
    *input_opt = (int32_t*)malloc(sizeof(**input_opt) * MAX_SB_SIZE * *input_stride);
    memset(*input, 0, MAX_SB_SIZE * *input_stride);
    memset(*input_opt, 0, MAX_SB_SIZE * *input_stride);
    eb_buf_random_s32(*input, MAX_SB_SIZE * *input_stride);
    memcpy(*input_opt, *input, sizeof(**input) * MAX_SB_SIZE * *input_stride);
}

static void uninit_data(int32_t *coeff, int32_t *coeff_opt,int32_t *stride) {
    free(coeff);
    free(coeff_opt);
}

static void uninit_output(uint16_t *output, uint16_t *output_opt) {
    free(output);
    free(output_opt);
}

static void init_output(uint16_t **output,uint16_t **output_opt,int32_t num) {
    *output = (uint16_t*)malloc(sizeof(**output) * MAX_SB_SIZE * num);
    *output_opt = (uint16_t*)malloc(sizeof(**output) * MAX_SB_SIZE * num);
    eb_buf_random_u16(*output, MAX_SB_SIZE * num);
    memcpy(*output_opt, *output, sizeof(**output) * MAX_SB_SIZE * num);
}

TEST(InverseTransformTest, av1_inv_txfm_2d_kernels)
{
    int32_t* coeff,* coeff_opt;
    uint16_t *output, *output_opt;
    int32_t stride;
    EbBool Result = EB_FALSE;
    for (int loop = 0; loop < 3; loop++) {          //Function Pairs
        for (int i = 0; i < 10; i++) {              //Number of Test Runs
            for (int x = 0; x < 3; x++) {           //Bit Depth
                switch (loop) {
                case 0://16x16
                    for (int j = 0; j < 16; j++) {
                        init_data(&coeff, &coeff_opt, &stride);
                        ASSERT(eb_buf_compare_s32(coeff, coeff_opt, MAX_SB_SIZE * stride) == 1);
                        init_output(&output, &output_opt, stride);
                        av1_inv_txfm_highbd_func_ptr_array_base[loop](coeff, output, stride, (TxType)txsize_16[j], bd[x]);
                        av1_inv_txfm_highbd_func_ptr_array_opt[loop](coeff, output_opt, stride, (TxType)txsize_16[j], bd[x]);
                        EXPECT_EQ(eb_buf_compare_u16(output, output_opt, MAX_SB_SIZE * stride), 1);
                        uninit_output(output, output_opt);
                        uninit_data(coeff, coeff_opt, &stride);
                    }
                break;
                case 1://32x32
                    for (int j = 0; j < 2; j++) {
                        init_data(&coeff, &coeff_opt, &stride);
                        ASSERT(eb_buf_compare_s32(coeff, coeff_opt, MAX_SB_SIZE * stride) == 1);
                        init_output(&output, &output_opt, stride);
                        av1_inv_txfm_highbd_func_ptr_array_base[loop](coeff, output, stride, (TxType)txsize_32[j], bd[x]);
                        av1_inv_txfm_highbd_func_ptr_array_opt[loop](coeff, output_opt, stride, (TxType)txsize_32[j], bd[x]);
                        EXPECT_EQ(eb_buf_compare_u16(output, output_opt,MAX_SB_SIZE * stride), 1);
                        uninit_output(output, output_opt);
                        uninit_data(coeff, coeff_opt, &stride);
                    }
                break;
                case 2://64x64
                    for (int j = 0; j < 1; j++) {
                        init_data(&coeff, &coeff_opt, &stride);
                        ASSERT(eb_buf_compare_s32(coeff, coeff_opt, MAX_SB_SIZE * stride) == 1);
                        init_output(&output, &output_opt, stride);
                        av1_inv_txfm_highbd_func_ptr_array_base[loop](coeff, output, stride, (TxType)txsize_64[j], bd[x]);
                        av1_inv_txfm_highbd_func_ptr_array_opt[loop](coeff, output_opt, stride, (TxType)txsize_64[j], bd[x]);
                        EXPECT_EQ(eb_buf_compare_u16(output, output_opt, MAX_SB_SIZE * stride), 1);
                        uninit_output(output, output_opt);
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
