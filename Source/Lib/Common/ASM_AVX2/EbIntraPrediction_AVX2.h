/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbIntraPrediction_AVX2_h
#define EbIntraPrediction_AVX2_h

#include "EbDefinitions.h"
#ifdef __cplusplus
extern "C" {
#endif
    // Weights are quadratic from '1' to '1 / block_size', scaled by
    // 2^sm_weight_log2_scale.
    static const int32_t sm_weight_log2_scale = 8;
    // max(block_size_wide[BLOCK_LARGEST], block_size_high[BLOCK_LARGEST])
#define MAX_BLOCK_DIM 64
/* clang-format off */
    static const uint8_t sm_weight_arrays[2 * MAX_BLOCK_DIM] = {
        // Unused, because we always offset by bs, which is at least 2.
        0, 0,
        // bs = 2
        255, 128,
        // bs = 4
        255, 149, 85, 64,
        // bs = 8
        255, 197, 146, 105, 73, 50, 37, 32,
        // bs = 16
        255, 225, 196, 170, 145, 123, 102, 84, 68, 54, 43, 33, 26, 20, 17, 16,
        // bs = 32
        255, 240, 225, 210, 196, 182, 169, 157, 145, 133, 122, 111, 101, 92, 83, 74,
        66, 59, 52, 45, 39, 34, 29, 25, 21, 17, 14, 12, 10, 9, 8, 8,
        // bs = 64
        255, 248, 240, 233, 225, 218, 210, 203, 196, 189, 182, 176, 169, 163, 156,
        150, 144, 138, 133, 127, 121, 116, 111, 106, 101, 96, 91, 86, 82, 77, 73, 69,
        65, 61, 57, 54, 50, 47, 44, 41, 38, 35, 32, 29, 27, 25, 22, 20, 18, 16, 15,
        13, 12, 10, 9, 8, 7, 6, 6, 5, 5, 4, 4, 4,
    };
    
#if !OIS_BASED_INTRA
    extern void intra_mode_planar_avx2_intrin(
        const uint32_t  size,                           //input parameter, denotes the size of the current PU
        uint8_t        *ref_samples,                    //input parameter, pointer to the reference samples
        uint8_t        *prediction_ptr,                 //output parameter, pointer to the prediction
        const uint32_t  prediction_buffer_stride,       //input parameter, denotes the stride for the prediction ptr
        const EbBool    skip);


    void intra_mode_planar_av1_avx2_intrin(
        const uint32_t  size,                           //input parameter, denotes the size of the current PU
        uint8_t        *ref_samples,                    //input parameter, pointer to the reference samples
        uint8_t        *dst,                            //output parameter, pointer to the prediction
        const uint32_t  prediction_buffer_stride,       //input parameter, denotes the stride for the prediction ptr
        const EbBool    skip);

    extern void intra_mode_angular_vertical_kernel_avx2_intrin(
        uint32_t      size,                             //input parameter, denotes the size of the current PU
        uint8_t      *ref_samp_main,                    //input parameter, pointer to the reference samples
        uint8_t      *prediction_ptr,                   //output parameter, pointer to the prediction
        uint32_t      prediction_buffer_stride,         //input parameter, denotes the stride for the prediction ptr
        const EbBool  skip,
        int32_t       intra_pred_angle);

    extern void intra_mode_angular_horizontal_kernel_avx2_intrin(
        uint32_t      size,                             //input parameter, denotes the size of the current PU
        uint8_t      *ref_samp_main,                    //input parameter, pointer to the reference samples
        uint8_t      *prediction_ptr,                   //output parameter, pointer to the prediction
        uint32_t      prediction_buffer_stride,         //input parameter, denotes the stride for the prediction ptr
        const EbBool  skip,
        int32_t       intra_pred_angle);

    extern void intra_mode_vertical_luma_avx2_intrin(
        const uint32_t  size,                           //input parameter, denotes the size of the current PU
        uint8_t        *ref_samples,                    //input parameter, pointer to the reference samples
        uint8_t        *prediction_ptr,                 //output parameter, pointer to the prediction
        const uint32_t  prediction_buffer_stride,       //input parameter, denotes the stride for the prediction ptr
        const EbBool    skip);                          //skip one row
   
    extern void intra_mode_dc_luma_avx2_intrin(
        const uint32_t  size,                           //input parameter, denotes the size of the current PU
        uint8_t        *ref_samples,                    //input parameter, pointer to the reference samples
        uint8_t        *prediction_ptr,                 //output parameter, pointer to the prediction
        const uint32_t  prediction_buffer_stride,       //input parameter, denotes the stride for the prediction ptr
        const EbBool    skip);


    extern void intra_mode_dc_32x32_av1_avx2_intrin(
        EbBool           is_left_availble,
        EbBool           is_above_availble,
        const uint32_t   size,                          //input parameter, denotes the size of the current PU
        uint8_t         *ref_samples,                   //input parameter, pointer to the reference samples
        uint8_t         *dst,                           //output parameter, pointer to the prediction
        const uint32_t   prediction_buffer_stride,      //input parameter, denotes the stride for the prediction ptr
        const EbBool     skip);                         //skip half rows

    void intra_mode_dc_64x64_av1_avx2_intrin(
        EbBool         is_left_availble,
        EbBool         is_above_availble,
        const uint32_t size,                            //input parameter, denotes the size of the current PU
        uint8_t       *ref_samples,                     //input parameter, pointer to the reference samples
        uint8_t       *dst,                             //output parameter, pointer to the prediction
        const uint32_t prediction_buffer_stride,        //input parameter, denotes the stride for the prediction ptr
        const EbBool   skip);

    extern void intra_mode_angular_2_avx2_intrin(
        const uint32_t  size,                           //input parameter, denotes the size of the current PU
        uint8_t        *ref_samples,                    //input parameter, pointer to the reference samples
        uint8_t        *prediction_ptr,                 //output parameter, pointer to the prediction
        const uint32_t  prediction_buffer_stride,       //input parameter, denotes the stride for the prediction ptr
        const EbBool    skip);

    extern void intra_mode_angular_18_avx2_intrin(
        const uint32_t  size,                           //input parameter, denotes the size of the current PU
        uint8_t        *ref_samples,                    //input parameter, pointer to the reference samples
        uint8_t        *prediction_ptr,                 //output parameter, pointer to the prediction
        const uint32_t  prediction_buffer_stride,       //input parameter, denotes the stride for the prediction ptr
        const EbBool    skip);

    extern void intra_mode_angular_34_avx2_intrin(
        const uint32_t  size,                           //input parameter, denotes the size of the current PU
        uint8_t        *ref_samples,                    //input parameter, pointer to the reference samples
        uint8_t        *prediction_ptr,                 //output parameter, pointer to the prediction
        const uint32_t  prediction_buffer_stride,       //input parameter, denotes the stride for the prediction ptr
        const EbBool    skip);

    uint32_t update_neighbor_dc_intra_pred_avx2_intrin(
        uint8_t  *y_intra_reference_array_reverse,
        uint16_t  input_height,
        uint16_t  stride_y,
        EbByte    buffer_y,
        uint16_t  origin_y,
        uint16_t  origin_x,
        uint32_t  src_origin_x,
        uint32_t  src_origin_y,
        uint32_t  block_size,
        EbAsm     asm_type);
#endif
    void intra_mode_angular_av1_z1_16bit_4x4_avx2(
        const uint32_t  size,
        uint16_t       *ref_samples,
        uint16_t       *dst,
        const uint32_t  prediction_buffer_stride,
        const EbBool    skip,
        uint16_t        dx,
        uint16_t        dy,
        uint16_t        bd);

    void intra_mode_angular_av1_z1_16bit_8x8_avx2(
        const uint32_t  size,
        uint16_t       *ref_samples,
        uint16_t       *dst,
        const uint32_t  prediction_buffer_stride,
        const EbBool    skip,
        uint16_t        dx,
        uint16_t        dy,
        uint16_t        bd);

    void intra_mode_angular_av1_z1_16bit_16x16_avx2(
        const uint32_t  size,
        uint16_t       *ref_samples,
        uint16_t       *dst,
        const uint32_t  prediction_buffer_stride,
        const EbBool    skip,
        uint16_t        dx,
        uint16_t        dy,
        uint16_t        bd);

    void intra_mode_angular_av1_z1_16bit_32x32_avx2(
        const uint32_t  size,
        uint16_t       *ref_samples,
        uint16_t       *dst,
        const uint32_t  prediction_buffer_stride,
        const EbBool    skip,
        uint16_t        dx,
        uint16_t        dy,
        uint16_t        bd);

    void intra_mode_angular_av1_z1_16bit_64x64_avx2(
        const uint32_t  size,
        uint16_t       *ref_samples,
        uint16_t       *dst,
        const uint32_t  prediction_buffer_stride,
        const EbBool    skip,
        uint16_t        dx,
        uint16_t        dy,
        uint16_t        bd);

    void intra_mode_angular_av1_z2_16bit_4x4_avx2(
        const uint32_t  size,
        uint16_t       *ref_samples,
        uint16_t       *dst,
        const uint32_t  prediction_buffer_stride,
        const EbBool    skip,
        uint16_t        dx,
        uint16_t        dy,
        uint16_t        bd);

    void intra_mode_angular_av1_z2_16bit_8x8_avx2(
        const uint32_t  size,
        uint16_t       *ref_samples,
        uint16_t       *dst,
        const uint32_t  prediction_buffer_stride,
        const EbBool    skip,
        uint16_t        dx,
        uint16_t        dy,
        uint16_t        bd);

    void intra_mode_angular_av1_z2_16bit_16x16_avx2(
        const uint32_t  size,
        uint16_t       *ref_samples,
        uint16_t       *dst,
        const uint32_t  prediction_buffer_stride,
        const EbBool    skip,
        uint16_t        dx,
        uint16_t        dy,
        uint16_t        bd);

    void intra_mode_angular_av1_z2_16bit_32x32_avx2(
        const uint32_t  size,
        uint16_t       *ref_samples,
        uint16_t       *dst,
        const uint32_t  prediction_buffer_stride,
        const EbBool    skip,
        uint16_t        dx,
        uint16_t        dy,
        uint16_t        bd);

    void intra_mode_angular_av1_z2_16bit_64x64_avx2(
        const uint32_t  size,
        uint16_t       *ref_samples,
        uint16_t       *dst,
        const uint32_t  prediction_buffer_stride,
        const EbBool    skip,
        uint16_t        dx,
        uint16_t        dy,
        uint16_t        bd);

    void intra_mode_angular_av1_z3_16bit_4x4_avx2(
        const uint32_t  size,
        uint16_t       *ref_samples,
        uint16_t       *dst,
        const uint32_t  prediction_buffer_stride,
        const EbBool    skip,
        uint16_t        dx,
        uint16_t        dy,
        uint16_t        bd);

    void intra_mode_angular_av1_z3_16bit_8x8_avx2(
        const uint32_t  size,
        uint16_t       *ref_samples,
        uint16_t       *dst,
        const uint32_t  prediction_buffer_stride,
        const EbBool    skip,
        uint16_t        dx,
        uint16_t        dy,
        uint16_t        bd);

    void intra_mode_angular_av1_z3_16bit_16x16_avx2(
        const uint32_t  size,
        uint16_t       *ref_samples,
        uint16_t       *dst,
        const uint32_t  prediction_buffer_stride,
        const EbBool    skip,
        uint16_t        dx,
        uint16_t        dy,
        uint16_t        bd);

    void intra_mode_angular_av1_z3_16bit_32x32_avx2(
        const uint32_t  size,
        uint16_t       *ref_samples,
        uint16_t       *dst,
        const uint32_t  prediction_buffer_stride,
        const EbBool    skip,
        uint16_t        dx,
        uint16_t        dy,
        uint16_t        bd);

    void intra_mode_angular_av1_z3_16bit_64x64_avx2(
        const uint32_t  size,
        uint16_t       *ref_samples,
        uint16_t       *dst,
        const uint32_t  prediction_buffer_stride,
        const EbBool    skip,
        uint16_t        dx,
        uint16_t        dy,
        uint16_t        bd);

#ifdef __cplusplus
}
#endif
#endif // EbIntraPrediction_AVX2_h
