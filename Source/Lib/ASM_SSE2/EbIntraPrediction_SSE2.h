/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbIntraPrediction_SSE2_h
#define EbIntraPrediction_SSE2_h

#include "EbDefinitions.h"

/*******************************************
* Function pointer Table
*******************************************/
#ifdef __cplusplus
extern "C" {
#endif

    extern void intra_mode_vertical_luma_sse2_intrin(
        const uint32_t  size,                         //input parameter, denotes the size of the current PU
        uint8_t        *ref_samples,                  //input parameter, pointer to the reference samples
        uint8_t        *prediction_ptr,               //output parameter, pointer to the prediction
        const uint32_t  prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
        const EbBool    skip);                        //skip one row
    
    extern void intra_mode_vertical_chroma_sse2_intrin(
        const uint32_t  size,                         //input parameter, denotes the size of the current PU
        uint8_t        *ref_samples,                  //input parameter, pointer to the reference samples
        uint8_t        *prediction_ptr,               //output parameter, pointer to the prediction
        const uint32_t  prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
        const EbBool    skip);                        //skip one row

    extern  void intra_mode_horizontal_luma_sse2_intrin(
        const uint32_t  size,                         //input parameter, denotes the size of the current PU
        uint8_t        *ref_samples,                  //input parameter, pointer to the reference samples
        uint8_t        *prediction_ptr,               //output parameter, pointer to the prediction
        const uint32_t  prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
        const EbBool    skip);                        //skip one row
    
    extern void intra_mode_horizontal_chroma_sse2_intrin(
        const uint32_t  size,                         //input parameter, denotes the size of the current PU
        uint8_t        *ref_samples,                  //input parameter, pointer to the reference samples
        uint8_t        *prediction_ptr,               //output parameter, pointer to the prediction
        const uint32_t  prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
        const EbBool    skip);                        //skip one row

    extern void intra_mode_dc_luma_sse2_intrin(
        const uint32_t  size,                         //input parameter, denotes the size of the current PU
        uint8_t        *ref_samples,                  //input parameter, pointer to the reference samples
        uint8_t        *prediction_ptr,               //output parameter, pointer to the prediction
        const uint32_t  prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
        const EbBool    skip);                        //skip one row

    extern void intra_mode_dc_16x16_av1_sse2_intrin(
        EbBool          is_left_availble,
        EbBool          is_above_availble,
        const uint32_t  size,                         //input parameter, denotes the size of the current PU
        uint8_t        *ref_samples,                  //input parameter, pointer to the reference samples
        uint8_t        *dst,                          //output parameter, pointer to the prediction
        const uint32_t  prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
        const EbBool    skip);                        //skip half rows

    extern void intra_mode_dc_8x8_av1_sse2_intrin(
        EbBool          is_left_availble,
        EbBool          is_above_availble,
        const uint32_t  size,                         //input parameter, denotes the size of the current PU
        uint8_t        *ref_samples,                  //input parameter, pointer to the reference samples
        uint8_t        *dst,                          //output parameter, pointer to the prediction
        const uint32_t  prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
        const EbBool    skip);                        //skip half rows

    extern void intra_mode_dc_4x4_av1_sse2_intrin(
        EbBool          is_left_availble,
        EbBool          is_above_availble,
        const uint32_t  size,                         //input parameter, denotes the size of the current PU
        uint8_t        *ref_samples,                  //input parameter, pointer to the reference samples
        uint8_t        *dst,                          //output parameter, pointer to the prediction
        const uint32_t  prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
        const EbBool    skip);                        //skip half rows

    void intra_mode_planar16bit_sse2_intrin(
        const uint32_t  size,                         //input parameter, denotes the size of the current PU
        uint16_t       *ref_samples,                  //input parameter, pointer to the reference samples
        uint16_t       *prediction_ptr,               //output parameter, pointer to the prediction
        const uint32_t  prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
        const EbBool    skip);                        //skip half rows

    extern void intra_mode_dc_chroma_sse2_intrin(
        const uint32_t  size,                         //input parameter, denotes the size of the current PU
        uint8_t        *ref_samples,                  //input parameter, pointer to the reference samples
        uint8_t        *prediction_ptr,               //output parameter, pointer to the prediction
        const uint32_t  prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
        const EbBool    skip);                        //skip one row


    extern void intra_mode_planar_sse2_intrin(
        const uint32_t  size,                         //input parameter, denotes the size of the current PU
        uint8_t        *ref_samples,                  //input parameter, pointer to the reference samples
        uint8_t        *prediction_ptr,               //output parameter, pointer to the prediction
        const uint32_t  prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
        const EbBool    skip);                        //skip one row

    extern  void intra_mode_angular_34_sse2_intrin(
        const uint32_t  size,                         //input parameter, denotes the size of the current PU
        uint8_t        *ref_samples,                  //input parameter, pointer to the reference samples
        uint8_t        *prediction_ptr,               //output parameter, pointer to the prediction
        const uint32_t  prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
        const EbBool    skip);                        //skip one row

    extern void intra_mode_angular_18_sse2_intrin(
        const uint32_t  size,                         //input parameter, denotes the size of the current PU
        uint8_t        *ref_samples,                  //input parameter, pointer to the reference samples
        uint8_t        *prediction_ptr,               //output parameter, pointer to the prediction
        const uint32_t  prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
        const EbBool    skip);                        //skip one row

    extern void intra_mode_angular_2_sse2_intrin(
        const uint32_t  size,                         //input parameter, denotes the size of the current PU
        uint8_t        *ref_samples,                  //input parameter, pointer to the reference samples
        uint8_t        *prediction_ptr,               //output parameter, pointer to the prediction
        const uint32_t  prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
        const EbBool    skip);                        //skip one row


#ifdef __cplusplus
}
#endif
#endif // EbIntraPrediction_SSE2_h
