// clang-format off
/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>
#include <string.h>

#include "EbThreads.h"
#include "EbReferenceObject.h"
#include "EbPictureBufferDesc.h"


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
    for (sampleCount = 0; sampleCount < reconHeight; sampleCount++)
        reconSamplesPtr[sampleCount * stride] = 0;
    // 4. Zero out the right column
    reconSamplesPtr = (uint16_t*)reconSamplesBufferPtr + top_padding * stride + left_padding + reconWidth;
    for (sampleCount = 0; sampleCount < reconHeight; sampleCount++)
        reconSamplesPtr[sampleCount * stride] = 0;
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
    for (sampleCount = 0; sampleCount < reconHeight; sampleCount++)
        reconSamplesPtr[sampleCount * stride] = 0;
    // 4. Zero out the right column
    reconSamplesPtr = reconSamplesBufferPtr + top_padding * stride + left_padding + reconWidth;
    for (sampleCount = 0; sampleCount < reconHeight; sampleCount++)
        reconSamplesPtr[sampleCount * stride] = 0;
}

void InitializeSamplesNeighboringReferencePicture(
    EbReferenceObject              *referenceObject,
    EbPictureBufferDescInitData    *pictureBufferDescInitDataPtr,
    EbBitDepthEnum                       bit_depth) {
    if (bit_depth == EB_10BIT) {
        InitializeSamplesNeighboringReferencePicture16Bit(
            referenceObject->reference_picture16bit->buffer_y,
            referenceObject->reference_picture16bit->stride_y,
            referenceObject->reference_picture16bit->width,
            referenceObject->reference_picture16bit->height,
            pictureBufferDescInitDataPtr->left_padding,
            pictureBufferDescInitDataPtr->top_padding);

        InitializeSamplesNeighboringReferencePicture16Bit(
            referenceObject->reference_picture16bit->buffer_cb,
            referenceObject->reference_picture16bit->stride_cb,
            referenceObject->reference_picture16bit->width >> 1,
            referenceObject->reference_picture16bit->height >> 1,
            pictureBufferDescInitDataPtr->left_padding >> 1,
            pictureBufferDescInitDataPtr->top_padding >> 1);

        InitializeSamplesNeighboringReferencePicture16Bit(
            referenceObject->reference_picture16bit->buffer_cr,
            referenceObject->reference_picture16bit->stride_cr,
            referenceObject->reference_picture16bit->width >> 1,
            referenceObject->reference_picture16bit->height >> 1,
            pictureBufferDescInitDataPtr->left_padding >> 1,
            pictureBufferDescInitDataPtr->top_padding >> 1);
    }
    else {
        InitializeSamplesNeighboringReferencePicture8Bit(
            referenceObject->reference_picture->buffer_y,
            referenceObject->reference_picture->stride_y,
            referenceObject->reference_picture->width,
            referenceObject->reference_picture->height,
            pictureBufferDescInitDataPtr->left_padding,
            pictureBufferDescInitDataPtr->top_padding);

        InitializeSamplesNeighboringReferencePicture8Bit(
            referenceObject->reference_picture->buffer_cb,
            referenceObject->reference_picture->stride_cb,
            referenceObject->reference_picture->width >> 1,
            referenceObject->reference_picture->height >> 1,
            pictureBufferDescInitDataPtr->left_padding >> 1,
            pictureBufferDescInitDataPtr->top_padding >> 1);

        InitializeSamplesNeighboringReferencePicture8Bit(
            referenceObject->reference_picture->buffer_cr,
            referenceObject->reference_picture->stride_cr,
            referenceObject->reference_picture->width >> 1,
            referenceObject->reference_picture->height >> 1,
            pictureBufferDescInitDataPtr->left_padding >> 1,
            pictureBufferDescInitDataPtr->top_padding >> 1);
    }
}

static void eb_reference_object_dctor(EbPtr p)
{
    EbReferenceObject *obj = (EbReferenceObject*)p;
    EB_DELETE(obj->reference_picture16bit);
    EB_DELETE(obj->reference_picture);
    EB_FREE_ALIGNED_ARRAY(obj->mvs);
    EB_DESTROY_MUTEX(obj->referenced_area_mutex);
}


/*****************************************
 * eb_picture_buffer_desc_ctor
 *  Initializes the Buffer Descriptor's
 *  values that are fixed for the life of
 *  the descriptor.
 *****************************************/
EbErrorType eb_reference_object_ctor(
    EbReferenceObject  *referenceObject,
    EbPtr   object_init_data_ptr)
{

    EbPictureBufferDescInitData    *pictureBufferDescInitDataPtr = (EbPictureBufferDescInitData*)object_init_data_ptr;
    EbPictureBufferDescInitData    pictureBufferDescInitData16BitPtr = *pictureBufferDescInitDataPtr;

    referenceObject->dctor = eb_reference_object_dctor;
    //TODO:12bit
    if (pictureBufferDescInitData16BitPtr.bit_depth == EB_10BIT) {
        // Hsan: set split_mode to 0 to construct the packed reference buffer (used @ EP)
        pictureBufferDescInitData16BitPtr.split_mode = EB_FALSE;
        EB_NEW(
            referenceObject->reference_picture16bit,
            eb_picture_buffer_desc_ctor,
            (EbPtr)&pictureBufferDescInitData16BitPtr);

        InitializeSamplesNeighboringReferencePicture(
            referenceObject,
            &pictureBufferDescInitData16BitPtr,
            pictureBufferDescInitData16BitPtr.bit_depth);

        // Hsan: set split_mode to 1 to construct the unpacked reference buffer (used @ MD)
        pictureBufferDescInitData16BitPtr.split_mode = EB_TRUE;
        EB_NEW(
            referenceObject->reference_picture,
            eb_picture_buffer_desc_ctor,
            (EbPtr)&pictureBufferDescInitData16BitPtr);
    }
    else {
        // Hsan: set split_mode to 0 to as 8BIT input
        pictureBufferDescInitDataPtr->split_mode = EB_FALSE;
        EB_NEW(
            referenceObject->reference_picture,
            eb_picture_buffer_desc_ctor,
            (EbPtr)pictureBufferDescInitDataPtr);

        InitializeSamplesNeighboringReferencePicture(
            referenceObject,
            pictureBufferDescInitDataPtr,
            pictureBufferDescInitData16BitPtr.bit_depth);
    }
    if (pictureBufferDescInitDataPtr->mfmv)
    {
        //MFMV map is 8x8 based.
        uint32_t mi_rows = referenceObject->reference_picture->height >> MI_SIZE_LOG2;
        uint32_t mi_cols = referenceObject->reference_picture->width >> MI_SIZE_LOG2;
        const int mem_size = ((mi_rows + 1) >> 1) * ((mi_cols + 1) >> 1);
        EB_CALLOC_ALIGNED_ARRAY(referenceObject->mvs, mem_size);
    }
    memset(&referenceObject->film_grain_params, 0, sizeof(referenceObject->film_grain_params));
    EB_CREATE_MUTEX(referenceObject->referenced_area_mutex);
    return EB_ErrorNone;
}

EbErrorType eb_reference_object_creator(
    EbPtr  *object_dbl_ptr,
    EbPtr   object_init_data_ptr)
{
    EbReferenceObject* obj;

    *object_dbl_ptr = NULL;
    EB_NEW(obj, eb_reference_object_ctor, object_init_data_ptr);
    *object_dbl_ptr = obj;

    return EB_ErrorNone;
}

static void eb_pa_reference_object_dctor(EbPtr p)
{
    EbPaReferenceObject* obj = (EbPaReferenceObject*)p;
    EB_DELETE(obj->input_padded_picture_ptr);
    EB_DELETE(obj->quarter_decimated_picture_ptr);
    EB_DELETE(obj->sixteenth_decimated_picture_ptr);
    EB_DELETE(obj->quarter_filtered_picture_ptr);
    EB_DELETE(obj->sixteenth_filtered_picture_ptr);
}

/*****************************************
 * eb_pa_reference_object_ctor
 *  Initializes the Buffer Descriptor's
 *  values that are fixed for the life of
 *  the descriptor.
 *****************************************/
EbErrorType eb_pa_reference_object_ctor(
    EbPaReferenceObject  *paReferenceObject,
    EbPtr   object_init_data_ptr)
{
    EbPictureBufferDescInitData       *pictureBufferDescInitDataPtr = (EbPictureBufferDescInitData*)object_init_data_ptr;

    paReferenceObject->dctor = eb_pa_reference_object_dctor;

    // Reference picture constructor
    EB_NEW(
        paReferenceObject->input_padded_picture_ptr,
        eb_picture_buffer_desc_ctor,
        (EbPtr)pictureBufferDescInitDataPtr);
    // Quarter Decim reference picture constructor
    EB_NEW(
        paReferenceObject->quarter_decimated_picture_ptr,
        eb_picture_buffer_desc_ctor,
        (EbPtr)(pictureBufferDescInitDataPtr + 1));
    EB_NEW(
        paReferenceObject->sixteenth_decimated_picture_ptr,
        eb_picture_buffer_desc_ctor,
        (EbPtr)(pictureBufferDescInitDataPtr + 2));
    // Quarter Filtered reference picture constructor
    if ((pictureBufferDescInitDataPtr + 1)->down_sampled_filtered) {
        EB_NEW(
            paReferenceObject->quarter_filtered_picture_ptr,
            eb_picture_buffer_desc_ctor,
            (EbPtr)(pictureBufferDescInitDataPtr + 1));
    }
    // Sixteenth Filtered reference picture constructor
    if ((pictureBufferDescInitDataPtr + 2)->down_sampled_filtered) {
        EB_NEW(
            paReferenceObject->sixteenth_filtered_picture_ptr,
            eb_picture_buffer_desc_ctor,
            (EbPtr)(pictureBufferDescInitDataPtr + 2));
    }

    return EB_ErrorNone;
}

EbErrorType eb_pa_reference_object_creator(
    EbPtr  *object_dbl_ptr,
    EbPtr   object_init_data_ptr)
{
    EbPaReferenceObject* obj;

    *object_dbl_ptr = NULL;
    EB_NEW(obj, eb_pa_reference_object_ctor, object_init_data_ptr);
    *object_dbl_ptr = obj;

    return EB_ErrorNone;
}
// clang-format on
