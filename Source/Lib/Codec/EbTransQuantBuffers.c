/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/


#include "EbDefinitions.h"
#include "EbTransQuantBuffers.h"
#include "EbPictureBufferDesc.h"


EbErrorType EbTransQuantBuffersCtor(
    EbTransQuantBuffers_t          *trans_quant_buffers_ptr)
{
    EbErrorType return_error = EB_ErrorNone;
    EbPictureBufferDescInitData_t transCoeffInitArray;
    transCoeffInitArray.maxWidth = SB_STRIDE_Y;
    transCoeffInitArray.maxHeight = SB_STRIDE_Y;
    transCoeffInitArray.bit_depth = EB_16BIT;
    transCoeffInitArray.bufferEnableMask = PICTURE_BUFFER_DESC_FULL_MASK;
    transCoeffInitArray.left_padding = 0;
    transCoeffInitArray.right_padding = 0;
    transCoeffInitArray.top_padding = 0;
    transCoeffInitArray.bot_padding = 0;
    transCoeffInitArray.splitMode = EB_FALSE;

    EbPictureBufferDescInitData_t ThirtyTwoBittransCoeffInitArray;
    ThirtyTwoBittransCoeffInitArray.maxWidth = SB_STRIDE_Y;
    ThirtyTwoBittransCoeffInitArray.maxHeight = SB_STRIDE_Y;
    ThirtyTwoBittransCoeffInitArray.bit_depth = EB_32BIT;
    ThirtyTwoBittransCoeffInitArray.bufferEnableMask = PICTURE_BUFFER_DESC_FULL_MASK;
    ThirtyTwoBittransCoeffInitArray.left_padding = 0;
    ThirtyTwoBittransCoeffInitArray.right_padding = 0;
    ThirtyTwoBittransCoeffInitArray.top_padding = 0;
    ThirtyTwoBittransCoeffInitArray.bot_padding = 0;
    ThirtyTwoBittransCoeffInitArray.splitMode = EB_FALSE;

    return_error = eb_picture_buffer_desc_ctor(
        (EbPtr*) &(trans_quant_buffers_ptr->tuTransCoeff2Nx2NPtr),
        (EbPtr)&ThirtyTwoBittransCoeffInitArray);
    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }

    return_error = eb_picture_buffer_desc_ctor(
        (EbPtr*) &(trans_quant_buffers_ptr->tuTransCoeffNxNPtr),
        (EbPtr)&ThirtyTwoBittransCoeffInitArray);
    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }


    return_error = eb_picture_buffer_desc_ctor(
        (EbPtr*) &(trans_quant_buffers_ptr->tuTransCoeffN2xN2Ptr),
        (EbPtr)&transCoeffInitArray);
    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }
    return_error = eb_picture_buffer_desc_ctor(
        (EbPtr*) &(trans_quant_buffers_ptr->tuQuantCoeffNxNPtr),
        (EbPtr)&transCoeffInitArray);
    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }

    return_error = eb_picture_buffer_desc_ctor(
        (EbPtr*) &(trans_quant_buffers_ptr->tuQuantCoeffN2xN2Ptr),
        (EbPtr)&transCoeffInitArray);
    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }

    return EB_ErrorNone;
}

