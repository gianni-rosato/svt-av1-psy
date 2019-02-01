/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbIntraPrediction_SSE2_h
#define EbIntraPrediction_SSE2_h

#include "EbDefinitions.h"

/*******************************************
* Function Pointer Table
*******************************************/
#ifdef __cplusplus
extern "C" {
#endif

    extern void IntraModeVerticalLuma_SSE2_INTRIN(
        const uint32_t      size,                   //input parameter, denotes the size of the current PU
        uint8_t            *refSamples,             //input parameter, pointer to the reference samples
        uint8_t            *prediction_ptr,          //output parameter, pointer to the prediction
        const uint32_t      predictionBufferStride, //input parameter, denotes the stride for the prediction ptr
        const EbBool     skip                    //skip one row
    );



    extern void IntraModeVerticalChroma_SSE2_INTRIN(
        const uint32_t      size,                   //input parameter, denotes the size of the current PU
        uint8_t            *refSamples,             //input parameter, pointer to the reference samples
        uint8_t            *prediction_ptr,          //output parameter, pointer to the prediction
        const uint32_t      predictionBufferStride, //input parameter, denotes the stride for the prediction ptr
        const EbBool     skip                    //skip one row
    );





    extern  void IntraModeHorizontalLuma_SSE2_INTRIN(
        const uint32_t      size,                   //input parameter, denotes the size of the current PU
        uint8_t            *refSamples,             //input parameter, pointer to the reference samples
        uint8_t            *prediction_ptr,          //output parameter, pointer to the prediction
        const uint32_t      predictionBufferStride, //input parameter, denotes the stride for the prediction ptr
        const EbBool     skip                    //skip one row
    );



    extern void IntraModeHorizontalChroma_SSE2_INTRIN(
        const uint32_t      size,                   //input parameter, denotes the size of the current PU
        uint8_t            *refSamples,             //input parameter, pointer to the reference samples
        uint8_t            *prediction_ptr,          //output parameter, pointer to the prediction
        const uint32_t      predictionBufferStride, //input parameter, denotes the stride for the prediction ptr
        const EbBool     skip                    //skip one row
    );


    extern void IntraModeDCLuma_SSE2_INTRIN(
        const uint32_t      size,                       //input parameter, denotes the size of the current PU
        uint8_t            *refSamples,                 //input parameter, pointer to the reference samples
        uint8_t            *prediction_ptr,              //output parameter, pointer to the prediction
        const uint32_t      predictionBufferStride,     //input parameter, denotes the stride for the prediction ptr
        const EbBool     skip);               //skip one row


    extern void IntraModeDC_16x16_AV1_SSE2_INTRIN(
        EbBool                         is_left_availble,
        EbBool                         is_above_availble,
        const uint32_t   size,                       //input parameter, denotes the size of the current PU
        uint8_t         *refSamples,                 //input parameter, pointer to the reference samples
        uint8_t         *dst,              //output parameter, pointer to the prediction
        const uint32_t   predictionBufferStride,     //input parameter, denotes the stride for the prediction ptr
        const EbBool  skip);              //skip half rows

    extern void IntraModeDC_8x8_AV1_SSE2_INTRIN(
        EbBool                         is_left_availble,
        EbBool                         is_above_availble,
        const uint32_t   size,                       //input parameter, denotes the size of the current PU
        uint8_t         *refSamples,                 //input parameter, pointer to the reference samples
        uint8_t         *dst,              //output parameter, pointer to the prediction
        const uint32_t   predictionBufferStride,     //input parameter, denotes the stride for the prediction ptr
        const EbBool  skip);              //skip half rows

    extern void IntraModeDC_4x4_AV1_SSE2_INTRIN(
        EbBool                         is_left_availble,
        EbBool                         is_above_availble,
        const uint32_t   size,                       //input parameter, denotes the size of the current PU
        uint8_t         *refSamples,                 //input parameter, pointer to the reference samples
        uint8_t         *dst,              //output parameter, pointer to the prediction
        const uint32_t   predictionBufferStride,     //input parameter, denotes the stride for the prediction ptr
        const EbBool  skip);              //skip half rows

    void IntraModePlanar16bit_SSE2_INTRIN(
        const uint32_t   size,                       //input parameter, denotes the size of the current PU
        uint16_t         *refSamples,                 //input parameter, pointer to the reference samples
        uint16_t         *prediction_ptr,              //output parameter, pointer to the prediction
        const uint32_t   predictionBufferStride,     //input parameter, denotes the stride for the prediction ptr
        const EbBool  skip);                       //skip half rows

    extern void IntraModeDCChroma_SSE2_INTRIN(
        const uint32_t      size,                       //input parameter, denotes the size of the current PU
        uint8_t            *refSamples,                 //input parameter, pointer to the reference samples
        uint8_t            *prediction_ptr,              //output parameter, pointer to the prediction
        const uint32_t      predictionBufferStride,     //input parameter, denotes the stride for the prediction ptr
        const EbBool     skip);               //skip one row


    extern void IntraModePlanar_SSE2_INTRIN(
        const uint32_t      size,                       //input parameter, denotes the size of the current PU
        uint8_t            *refSamples,                 //input parameter, pointer to the reference samples
        uint8_t            *prediction_ptr,              //output parameter, pointer to the prediction
        const uint32_t      predictionBufferStride,     //input parameter, denotes the stride for the prediction ptr
        const EbBool     skip);                     //skip one row

    extern  void IntraModeAngular_34_SSE2_INTRIN(
        const uint32_t      size,                       //input parameter, denotes the size of the current PU
        uint8_t            *refSamples,                 //input parameter, pointer to the reference samples
        uint8_t            *prediction_ptr,              //output parameter, pointer to the prediction
        const uint32_t      predictionBufferStride,     //input parameter, denotes the stride for the prediction ptr
        const EbBool     skip);                     //skip one row



    extern void IntraModeAngular_18_SSE2_INTRIN(
        const uint32_t      size,                       //input parameter, denotes the size of the current PU
        uint8_t            *refSamples,                 //input parameter, pointer to the reference samples
        uint8_t            *prediction_ptr,              //output parameter, pointer to the prediction
        const uint32_t      predictionBufferStride,     //input parameter, denotes the stride for the prediction ptr
        const EbBool     skip);                   //skip one row



    extern void IntraModeAngular_2_SSE2_INTRIN(
        const uint32_t      size,                       //input parameter, denotes the size of the current PU
        uint8_t            *refSamples,                 //input parameter, pointer to the reference samples
        uint8_t            *prediction_ptr,              //output parameter, pointer to the prediction
        const uint32_t      predictionBufferStride,     //input parameter, denotes the stride for the prediction ptr
        const EbBool     skip);                      //skip one row








#ifdef __cplusplus
}
#endif
#endif // EbIntraPrediction_SSE2_h
