/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>
#include <string.h>

#include "EbPictureBufferDesc.h"
#include "EbReferenceObject.h"

void InitializeSamplesNeighboringReferencePicture16Bit(
    EbByte  reconSamplesBufferPtr,
    uint16_t   stride,
    uint16_t   reconWidth,
    uint16_t   reconHeight,
    uint16_t   left_padding,
    uint16_t   top_padding) {

    uint16_t  *reconSamplesPtr;
    uint16_t   sampleCount;

    // 1. Zero out the top row
    reconSamplesPtr = (uint16_t*)reconSamplesBufferPtr + (top_padding - 1) * stride + left_padding - 1;
    EB_MEMSET((uint8_t*)reconSamplesPtr, 0, sizeof(uint16_t)*(1 + reconWidth + 1));

    // 2. Zero out the bottom row
    reconSamplesPtr = (uint16_t*)reconSamplesBufferPtr + (top_padding + reconHeight) * stride + left_padding - 1;
    EB_MEMSET((uint8_t*)reconSamplesPtr, 0, sizeof(uint16_t)*(1 + reconWidth + 1));

    // 3. Zero out the left column
    reconSamplesPtr = (uint16_t*)reconSamplesBufferPtr + top_padding * stride + left_padding - 1;
    for (sampleCount = 0; sampleCount < reconHeight; sampleCount++) {
        reconSamplesPtr[sampleCount * stride] = 0;
    }

    // 4. Zero out the right column
    reconSamplesPtr = (uint16_t*)reconSamplesBufferPtr + top_padding * stride + left_padding + reconWidth;
    for (sampleCount = 0; sampleCount < reconHeight; sampleCount++) {
        reconSamplesPtr[sampleCount * stride] = 0;
    }
}

void InitializeSamplesNeighboringReferencePicture8Bit(
    EbByte  reconSamplesBufferPtr,
    uint16_t   stride,
    uint16_t   reconWidth,
    uint16_t   reconHeight,
    uint16_t   left_padding,
    uint16_t   top_padding) {

    uint8_t   *reconSamplesPtr;
    uint16_t   sampleCount;

    // 1. Zero out the top row
    reconSamplesPtr = reconSamplesBufferPtr + (top_padding - 1) * stride + left_padding - 1;
    EB_MEMSET(reconSamplesPtr, 0, sizeof(uint8_t)*(1 + reconWidth + 1));

    // 2. Zero out the bottom row
    reconSamplesPtr = reconSamplesBufferPtr + (top_padding + reconHeight) * stride + left_padding - 1;
    EB_MEMSET(reconSamplesPtr, 0, sizeof(uint8_t)*(1 + reconWidth + 1));

    // 3. Zero out the left column
    reconSamplesPtr = reconSamplesBufferPtr + top_padding * stride + left_padding - 1;
    for (sampleCount = 0; sampleCount < reconHeight; sampleCount++) {
        reconSamplesPtr[sampleCount * stride] = 0;
    }

    // 4. Zero out the right column
    reconSamplesPtr = reconSamplesBufferPtr + top_padding * stride + left_padding + reconWidth;
    for (sampleCount = 0; sampleCount < reconHeight; sampleCount++) {
        reconSamplesPtr[sampleCount * stride] = 0;
    }
}

void InitializeSamplesNeighboringReferencePicture(
    EbReferenceObject_t              *referenceObject,
    EbPictureBufferDescInitData_t    *pictureBufferDescInitDataPtr,
    EB_BITDEPTH                       bit_depth) {

    if (bit_depth == EB_10BIT) {

        InitializeSamplesNeighboringReferencePicture16Bit(
            referenceObject->referencePicture16bit->buffer_y,
            referenceObject->referencePicture16bit->stride_y,
            referenceObject->referencePicture16bit->width,
            referenceObject->referencePicture16bit->height,
            pictureBufferDescInitDataPtr->left_padding,
            pictureBufferDescInitDataPtr->top_padding);

        InitializeSamplesNeighboringReferencePicture16Bit(
            referenceObject->referencePicture16bit->bufferCb,
            referenceObject->referencePicture16bit->strideCb,
            referenceObject->referencePicture16bit->width >> 1,
            referenceObject->referencePicture16bit->height >> 1,
            pictureBufferDescInitDataPtr->left_padding >> 1,
            pictureBufferDescInitDataPtr->top_padding >> 1);

        InitializeSamplesNeighboringReferencePicture16Bit(
            referenceObject->referencePicture16bit->bufferCr,
            referenceObject->referencePicture16bit->strideCr,
            referenceObject->referencePicture16bit->width >> 1,
            referenceObject->referencePicture16bit->height >> 1,
            pictureBufferDescInitDataPtr->left_padding >> 1,
            pictureBufferDescInitDataPtr->top_padding >> 1);
    }
    else {

        InitializeSamplesNeighboringReferencePicture8Bit(
            referenceObject->referencePicture->buffer_y,
            referenceObject->referencePicture->stride_y,
            referenceObject->referencePicture->width,
            referenceObject->referencePicture->height,
            pictureBufferDescInitDataPtr->left_padding,
            pictureBufferDescInitDataPtr->top_padding);

        InitializeSamplesNeighboringReferencePicture8Bit(
            referenceObject->referencePicture->bufferCb,
            referenceObject->referencePicture->strideCb,
            referenceObject->referencePicture->width >> 1,
            referenceObject->referencePicture->height >> 1,
            pictureBufferDescInitDataPtr->left_padding >> 1,
            pictureBufferDescInitDataPtr->top_padding >> 1);

        InitializeSamplesNeighboringReferencePicture8Bit(
            referenceObject->referencePicture->bufferCr,
            referenceObject->referencePicture->strideCr,
            referenceObject->referencePicture->width >> 1,
            referenceObject->referencePicture->height >> 1,
            pictureBufferDescInitDataPtr->left_padding >> 1,
            pictureBufferDescInitDataPtr->top_padding >> 1);
    }
}


/*****************************************
 * eb_picture_buffer_desc_ctor
 *  Initializes the Buffer Descriptor's
 *  values that are fixed for the life of
 *  the descriptor.
 *****************************************/
EbErrorType eb_reference_object_ctor(
    EbPtr  *object_dbl_ptr,
    EbPtr   object_init_data_ptr)
{

    EbReferenceObject_t              *referenceObject;
    EbPictureBufferDescInitData_t    *pictureBufferDescInitDataPtr = (EbPictureBufferDescInitData_t*)object_init_data_ptr;
    EbPictureBufferDescInitData_t    pictureBufferDescInitData16BitPtr = *pictureBufferDescInitDataPtr;
    EbErrorType return_error = EB_ErrorNone;
    EB_MALLOC(EbReferenceObject_t*, referenceObject, sizeof(EbReferenceObject_t), EB_N_PTR);

    *object_dbl_ptr = (EbPtr)referenceObject;


    if (pictureBufferDescInitData16BitPtr.bit_depth == EB_10BIT) {

        return_error = eb_picture_buffer_desc_ctor(
            (EbPtr*)&(referenceObject->referencePicture16bit),
            (EbPtr)&pictureBufferDescInitData16BitPtr);

        InitializeSamplesNeighboringReferencePicture(
            referenceObject,
            &pictureBufferDescInitData16BitPtr,
            pictureBufferDescInitData16BitPtr.bit_depth);

    }
    else {

        return_error = eb_picture_buffer_desc_ctor(
            (EbPtr*)&(referenceObject->referencePicture),
            (EbPtr)pictureBufferDescInitDataPtr);

        InitializeSamplesNeighboringReferencePicture(
            referenceObject,
            pictureBufferDescInitDataPtr,
            pictureBufferDescInitData16BitPtr.bit_depth);
    }
    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }



    // Allocate SB based TMVP map
    EB_MALLOC(TmvpUnit_t *, referenceObject->tmvpMap, (sizeof(TmvpUnit_t) * (((pictureBufferDescInitDataPtr->maxWidth + (64 - 1)) >> 6) * ((pictureBufferDescInitDataPtr->maxHeight + (64 - 1)) >> 6))), EB_N_PTR);

    //RESTRICT THIS TO M4
    {
        EbPictureBufferDescInitData_t bufDesc;

        bufDesc.maxWidth = pictureBufferDescInitDataPtr->maxWidth;
        bufDesc.maxHeight = pictureBufferDescInitDataPtr->maxHeight;
        bufDesc.bit_depth = EB_8BIT;
        bufDesc.bufferEnableMask = PICTURE_BUFFER_DESC_FULL_MASK;
        bufDesc.left_padding = pictureBufferDescInitDataPtr->left_padding;
        bufDesc.right_padding = pictureBufferDescInitDataPtr->right_padding;
        bufDesc.top_padding = pictureBufferDescInitDataPtr->top_padding;
        bufDesc.bot_padding = pictureBufferDescInitDataPtr->bot_padding;
        bufDesc.splitMode = 0;


        return_error = eb_picture_buffer_desc_ctor((EbPtr*)&(referenceObject->refDenSrcPicture),
            (EbPtr)&bufDesc);
        if (return_error == EB_ErrorInsufficientResources)
            return EB_ErrorInsufficientResources;
    }

    memset(&referenceObject->film_grain_params, 0, sizeof(referenceObject->film_grain_params));

    return EB_ErrorNone;
}

/*****************************************
 * eb_pa_reference_object_ctor
 *  Initializes the Buffer Descriptor's
 *  values that are fixed for the life of
 *  the descriptor.
 *****************************************/
EbErrorType eb_pa_reference_object_ctor(
    EbPtr  *object_dbl_ptr,
    EbPtr   object_init_data_ptr)
{

    EbPaReferenceObject_t               *paReferenceObject;
    EbPictureBufferDescInitData_t       *pictureBufferDescInitDataPtr = (EbPictureBufferDescInitData_t*)object_init_data_ptr;
    EbErrorType return_error = EB_ErrorNone;
    EB_MALLOC(EbPaReferenceObject_t*, paReferenceObject, sizeof(EbPaReferenceObject_t), EB_N_PTR);
    *object_dbl_ptr = (EbPtr)paReferenceObject;

    // Reference picture constructor
    return_error = eb_picture_buffer_desc_ctor(
        (EbPtr*) &(paReferenceObject->inputPaddedPicturePtr),
        (EbPtr)pictureBufferDescInitDataPtr);
    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }

    // Quarter Decim reference picture constructor
    paReferenceObject->quarterDecimatedPicturePtr = (EbPictureBufferDesc_t*)EB_NULL;
    return_error = eb_picture_buffer_desc_ctor(
        (EbPtr*) &(paReferenceObject->quarterDecimatedPicturePtr),
        (EbPtr)(pictureBufferDescInitDataPtr + 1));
    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }

    // Sixteenth Decim reference picture constructor
    paReferenceObject->sixteenthDecimatedPicturePtr = (EbPictureBufferDesc_t*)EB_NULL;
    return_error = eb_picture_buffer_desc_ctor(
        (EbPtr*) &(paReferenceObject->sixteenthDecimatedPicturePtr),
        (EbPtr)(pictureBufferDescInitDataPtr + 2));
    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }

    return EB_ErrorNone;
}



