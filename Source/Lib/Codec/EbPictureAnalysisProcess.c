/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>
#include <string.h>

#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"
#include "EbPictureControlSet.h"
#include "EbSequenceControlSet.h"
#include "EbPictureBufferDesc.h"

#include "EbResourceCoordinationResults.h"
#include "EbPictureAnalysisProcess.h"
#include "EbPictureAnalysisResults.h"
#include "EbMcp.h"
#include "EbMotionEstimation.h"
#include "EbReferenceObject.h"

#include "EbComputeMean.h"
#include "EbMeSadCalculation.h"
#include "EbComputeMean_SSE2.h"
#include "EbCombinedAveragingSAD_Intrinsic_AVX2.h"

#define VARIANCE_PRECISION        16
#define  LCU_LOW_VAR_TH                5
#define  PIC_LOW_VAR_PERCENTAGE_TH    60
#define    FLAT_MAX_VAR            50
#define FLAT_MAX_VAR_DECIM        (50-00)
#define    NOISE_MIN_LEVEL            70000//120000
#define NOISE_MIN_LEVEL_DECIM   (70000+000000)//(120000+000000)
#define    NOISE_MIN_LEVEL_M6_M7       120000
#define NOISE_MIN_LEVEL_DECIM_M6_M7    (120000+000000)
#define DENOISER_QP_TH            29
#define DENOISER_BITRATE_TH        14000000
#define SAMPLE_THRESHOLD_PRECENT_BORDER_LINE      15
#define SAMPLE_THRESHOLD_PRECENT_TWO_BORDER_LINES 10

/************************************************
* Picture Analysis Context Constructor
************************************************/
EbErrorType PictureAnalysisContextCtor(
    EbPictureBufferDescInitData_t * inputPictureBufferDescInitData,
    EbBool                         denoiseFlag,
    PictureAnalysisContext_t **context_dbl_ptr,
    EbFifo_t *resourceCoordinationResultsInputFifoPtr,
    EbFifo_t *pictureAnalysisResultsOutputFifoPtr)
{
    PictureAnalysisContext_t *context_ptr;
    EB_MALLOC(PictureAnalysisContext_t*, context_ptr, sizeof(PictureAnalysisContext_t), EB_N_PTR);
    *context_dbl_ptr = context_ptr;

    context_ptr->resourceCoordinationResultsInputFifoPtr = resourceCoordinationResultsInputFifoPtr;
    context_ptr->pictureAnalysisResultsOutputFifoPtr = pictureAnalysisResultsOutputFifoPtr;

    EbErrorType return_error = EB_ErrorNone;

    if (denoiseFlag == EB_TRUE) {

        //denoised
        return_error = EbPictureBufferDescCtor(
            (EbPtr*)&(context_ptr->denoisedPicturePtr),
            (EbPtr)inputPictureBufferDescInitData);

        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }

        //luma buffer could re-used to process chroma
        context_ptr->denoisedPicturePtr->bufferCb = context_ptr->denoisedPicturePtr->bufferY;
        context_ptr->denoisedPicturePtr->bufferCr = context_ptr->denoisedPicturePtr->bufferY + context_ptr->denoisedPicturePtr->chromaSize;

        // noise
        inputPictureBufferDescInitData->maxHeight = BLOCK_SIZE_64;
        inputPictureBufferDescInitData->bufferEnableMask = PICTURE_BUFFER_DESC_Y_FLAG;

        return_error = EbPictureBufferDescCtor(
            (EbPtr*)&(context_ptr->noisePicturePtr),
            (EbPtr)inputPictureBufferDescInitData);

        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }
    return EB_ErrorNone;


}

/************************************************
 * Picture Analysis Context Destructor
 ************************************************/

 /********************************************
  * Decimation2D
  *      decimates the input
  ********************************************/
void Decimation2D(
    uint8_t *  input_samples,      // input parameter, input samples Ptr
    uint32_t   inputStride,       // input parameter, input stride
    uint32_t   inputAreaWidth,    // input parameter, input area width
    uint32_t   inputAreaHeight,   // input parameter, input area height
    uint8_t *  decimSamples,      // output parameter, decimated samples Ptr
    uint32_t   decimStride,       // input parameter, output stride
    uint32_t   decimStep)        // input parameter, area height
{

    uint32_t horizontal_index;
    uint32_t vertical_index;


    for (vertical_index = 0; vertical_index < inputAreaHeight; vertical_index += decimStep) {
        for (horizontal_index = 0; horizontal_index < inputAreaWidth; horizontal_index += decimStep) {

            decimSamples[(horizontal_index >> (decimStep >> 1))] = input_samples[horizontal_index];

        }
        input_samples += (inputStride << (decimStep >> 1));
        decimSamples += decimStride;
    }

    return;
}

/********************************************
* CalculateHistogram
*      creates n-bins histogram for the input
********************************************/
void CalculateHistogram(
    uint8_t *  input_samples,      // input parameter, input samples Ptr
    uint32_t   inputAreaWidth,    // input parameter, input area width
    uint32_t   inputAreaHeight,   // input parameter, input area height
    uint32_t   stride,            // input parameter, input stride
    uint8_t    decimStep,         // input parameter, area height
    uint32_t  *histogram,            // output parameter, output histogram
    uint64_t  *sum)

{

    uint32_t horizontal_index;
    uint32_t vertical_index;
    *sum = 0;

    for (vertical_index = 0; vertical_index < inputAreaHeight; vertical_index += decimStep) {
        for (horizontal_index = 0; horizontal_index < inputAreaWidth; horizontal_index += decimStep) {
            ++(histogram[input_samples[horizontal_index]]);
            *sum += input_samples[horizontal_index];
        }
        input_samples += (stride << (decimStep >> 1));
    }

    return;
}


uint64_t ComputeVariance32x32(
    EbPictureBufferDesc_t       *inputPaddedPicturePtr,         // input parameter, Input Padded Picture
    uint32_t                       inputLumaOriginIndex,          // input parameter, SB index, used to point to source/reference samples
    uint64_t                        *variance8x8,
    EbAsm                         asm_type)
{

    uint32_t blockIndex;

    uint64_t meanOf8x8Blocks[16];
    uint64_t meanOf8x8SquaredValuesBlocks[16];

    uint64_t meanOf16x16Blocks[4];
    uint64_t meanOf16x16SquaredValuesBlocks[4];

    uint64_t meanOf32x32Blocks;
    uint64_t meanOf32x32SquaredValuesBlocks;
    /////////////////////////////////////////////
    // (0,0)
    blockIndex = inputLumaOriginIndex;

    meanOf8x8Blocks[0] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
    meanOf8x8SquaredValuesBlocks[0] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

    // (0,1)
    blockIndex = blockIndex + 8;
    meanOf8x8Blocks[1] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
    meanOf8x8SquaredValuesBlocks[1] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

    // (0,2)
    blockIndex = blockIndex + 8;
    meanOf8x8Blocks[2] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
    meanOf8x8SquaredValuesBlocks[2] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

    // (0,3)
    blockIndex = blockIndex + 8;
    meanOf8x8Blocks[3] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
    meanOf8x8SquaredValuesBlocks[3] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);



    // (1,0)
    blockIndex = inputLumaOriginIndex + (inputPaddedPicturePtr->strideY << 3);
    meanOf8x8Blocks[4] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
    meanOf8x8SquaredValuesBlocks[4] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

    // (1,1)
    blockIndex = blockIndex + 8;
    meanOf8x8Blocks[5] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
    meanOf8x8SquaredValuesBlocks[5] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

    // (1,2)
    blockIndex = blockIndex + 8;
    meanOf8x8Blocks[6] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
    meanOf8x8SquaredValuesBlocks[6] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

    // (1,3)
    blockIndex = blockIndex + 8;
    meanOf8x8Blocks[7] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
    meanOf8x8SquaredValuesBlocks[7] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);



    // (2,0)
    blockIndex = inputLumaOriginIndex + (inputPaddedPicturePtr->strideY << 4);
    meanOf8x8Blocks[8] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
    meanOf8x8SquaredValuesBlocks[8] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

    // (2,1)
    blockIndex = blockIndex + 8;
    meanOf8x8Blocks[9] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
    meanOf8x8SquaredValuesBlocks[9] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

    // (2,2)
    blockIndex = blockIndex + 8;
    meanOf8x8Blocks[10] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
    meanOf8x8SquaredValuesBlocks[10] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

    // (2,3)
    blockIndex = blockIndex + 8;
    meanOf8x8Blocks[11] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
    meanOf8x8SquaredValuesBlocks[11] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);



    // (3,0)
    blockIndex = inputLumaOriginIndex + (inputPaddedPicturePtr->strideY << 3) + (inputPaddedPicturePtr->strideY << 4);
    meanOf8x8Blocks[12] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
    meanOf8x8SquaredValuesBlocks[12] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

    // (3,1)
    blockIndex = blockIndex + 8;
    meanOf8x8Blocks[13] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
    meanOf8x8SquaredValuesBlocks[13] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

    // (3,2)
    blockIndex = blockIndex + 8;
    meanOf8x8Blocks[14] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
    meanOf8x8SquaredValuesBlocks[14] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

    // (3,3)
    blockIndex = blockIndex + 8;
    meanOf8x8Blocks[15] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
    meanOf8x8SquaredValuesBlocks[15] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);


    /////////////////////////////////////////////

    variance8x8[0] = meanOf8x8SquaredValuesBlocks[0] - (meanOf8x8Blocks[0] * meanOf8x8Blocks[0]);
    variance8x8[1] = meanOf8x8SquaredValuesBlocks[1] - (meanOf8x8Blocks[1] * meanOf8x8Blocks[1]);
    variance8x8[2] = meanOf8x8SquaredValuesBlocks[2] - (meanOf8x8Blocks[2] * meanOf8x8Blocks[2]);
    variance8x8[3] = meanOf8x8SquaredValuesBlocks[3] - (meanOf8x8Blocks[3] * meanOf8x8Blocks[3]);
    variance8x8[4] = meanOf8x8SquaredValuesBlocks[4] - (meanOf8x8Blocks[4] * meanOf8x8Blocks[4]);
    variance8x8[5] = meanOf8x8SquaredValuesBlocks[5] - (meanOf8x8Blocks[5] * meanOf8x8Blocks[5]);
    variance8x8[6] = meanOf8x8SquaredValuesBlocks[6] - (meanOf8x8Blocks[6] * meanOf8x8Blocks[6]);
    variance8x8[7] = meanOf8x8SquaredValuesBlocks[7] - (meanOf8x8Blocks[7] * meanOf8x8Blocks[7]);
    variance8x8[8] = meanOf8x8SquaredValuesBlocks[8] - (meanOf8x8Blocks[8] * meanOf8x8Blocks[8]);
    variance8x8[9] = meanOf8x8SquaredValuesBlocks[9] - (meanOf8x8Blocks[9] * meanOf8x8Blocks[9]);
    variance8x8[10] = meanOf8x8SquaredValuesBlocks[10] - (meanOf8x8Blocks[10] * meanOf8x8Blocks[10]);
    variance8x8[11] = meanOf8x8SquaredValuesBlocks[11] - (meanOf8x8Blocks[11] * meanOf8x8Blocks[11]);
    variance8x8[12] = meanOf8x8SquaredValuesBlocks[12] - (meanOf8x8Blocks[12] * meanOf8x8Blocks[12]);
    variance8x8[13] = meanOf8x8SquaredValuesBlocks[13] - (meanOf8x8Blocks[13] * meanOf8x8Blocks[13]);
    variance8x8[14] = meanOf8x8SquaredValuesBlocks[14] - (meanOf8x8Blocks[14] * meanOf8x8Blocks[14]);
    variance8x8[15] = meanOf8x8SquaredValuesBlocks[15] - (meanOf8x8Blocks[15] * meanOf8x8Blocks[15]);

    // 16x16
    meanOf16x16Blocks[0] = (meanOf8x8Blocks[0] + meanOf8x8Blocks[1] + meanOf8x8Blocks[8] + meanOf8x8Blocks[9]) >> 2;
    meanOf16x16Blocks[1] = (meanOf8x8Blocks[2] + meanOf8x8Blocks[3] + meanOf8x8Blocks[10] + meanOf8x8Blocks[11]) >> 2;
    meanOf16x16Blocks[2] = (meanOf8x8Blocks[4] + meanOf8x8Blocks[5] + meanOf8x8Blocks[12] + meanOf8x8Blocks[13]) >> 2;
    meanOf16x16Blocks[3] = (meanOf8x8Blocks[6] + meanOf8x8Blocks[7] + meanOf8x8Blocks[14] + meanOf8x8Blocks[15]) >> 2;


    meanOf16x16SquaredValuesBlocks[0] = (meanOf8x8SquaredValuesBlocks[0] + meanOf8x8SquaredValuesBlocks[1] + meanOf8x8SquaredValuesBlocks[8] + meanOf8x8SquaredValuesBlocks[9]) >> 2;
    meanOf16x16SquaredValuesBlocks[1] = (meanOf8x8SquaredValuesBlocks[2] + meanOf8x8SquaredValuesBlocks[3] + meanOf8x8SquaredValuesBlocks[10] + meanOf8x8SquaredValuesBlocks[11]) >> 2;
    meanOf16x16SquaredValuesBlocks[2] = (meanOf8x8SquaredValuesBlocks[4] + meanOf8x8SquaredValuesBlocks[5] + meanOf8x8SquaredValuesBlocks[12] + meanOf8x8SquaredValuesBlocks[13]) >> 2;
    meanOf16x16SquaredValuesBlocks[3] = (meanOf8x8SquaredValuesBlocks[6] + meanOf8x8SquaredValuesBlocks[7] + meanOf8x8SquaredValuesBlocks[14] + meanOf8x8SquaredValuesBlocks[15]) >> 2;

    // 32x32
    meanOf32x32Blocks = (meanOf16x16Blocks[0] + meanOf16x16Blocks[1] + meanOf16x16Blocks[2] + meanOf16x16Blocks[3]) >> 2;


    meanOf32x32SquaredValuesBlocks = (meanOf16x16SquaredValuesBlocks[0] + meanOf16x16SquaredValuesBlocks[1] + meanOf16x16SquaredValuesBlocks[2] + meanOf16x16SquaredValuesBlocks[3]) >> 2;


    return (meanOf32x32SquaredValuesBlocks - (meanOf32x32Blocks * meanOf32x32Blocks));
}

uint64_t ComputeVariance16x16(
    EbPictureBufferDesc_t       *inputPaddedPicturePtr,         // input parameter, Input Padded Picture
    uint32_t                       inputLumaOriginIndex,          // input parameter, SB index, used to point to source/reference samples
    uint64_t                        *variance8x8,
    EbAsm                         asm_type)
{

    uint32_t blockIndex;

    uint64_t meanOf8x8Blocks[4];
    uint64_t meanOf8x8SquaredValuesBlocks[4];

    uint64_t meanOf16x16Blocks;
    uint64_t meanOf16x16SquaredValuesBlocks;

    // (0,0)
    blockIndex = inputLumaOriginIndex;

    meanOf8x8Blocks[0] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
    meanOf8x8SquaredValuesBlocks[0] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

    // (0,1)
    blockIndex = blockIndex + 8;
    meanOf8x8Blocks[1] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
    meanOf8x8SquaredValuesBlocks[1] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

    // (1,0)
    blockIndex = inputLumaOriginIndex + (inputPaddedPicturePtr->strideY << 3);
    meanOf8x8Blocks[2] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
    meanOf8x8SquaredValuesBlocks[2] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

    // (1,1)
    blockIndex = blockIndex + 8;
    meanOf8x8Blocks[3] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
    meanOf8x8SquaredValuesBlocks[3] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

    variance8x8[0] = meanOf8x8SquaredValuesBlocks[0] - (meanOf8x8Blocks[0] * meanOf8x8Blocks[0]);
    variance8x8[1] = meanOf8x8SquaredValuesBlocks[1] - (meanOf8x8Blocks[1] * meanOf8x8Blocks[1]);
    variance8x8[2] = meanOf8x8SquaredValuesBlocks[2] - (meanOf8x8Blocks[2] * meanOf8x8Blocks[2]);
    variance8x8[3] = meanOf8x8SquaredValuesBlocks[3] - (meanOf8x8Blocks[3] * meanOf8x8Blocks[3]);

    // 16x16
    meanOf16x16Blocks = (meanOf8x8Blocks[0] + meanOf8x8Blocks[1] + meanOf8x8Blocks[2] + meanOf8x8Blocks[3]) >> 2;
    meanOf16x16SquaredValuesBlocks = (meanOf8x8SquaredValuesBlocks[0] + meanOf8x8SquaredValuesBlocks[1] + meanOf8x8SquaredValuesBlocks[2] + meanOf8x8SquaredValuesBlocks[3]) >> 2;

    return (meanOf16x16SquaredValuesBlocks - (meanOf16x16Blocks * meanOf16x16Blocks));
}

/*******************************************
ComputeVariance64x64
this function is exactly same as
PictureAnalysisComputeVarianceLcu excpet it
does not store data for every block,
just returns the 64x64 data point
*******************************************/
uint64_t ComputeVariance64x64(
    SequenceControlSet_t        *sequence_control_set_ptr,
    EbPictureBufferDesc_t       *inputPaddedPicturePtr,         // input parameter, Input Padded Picture
    uint32_t                       inputLumaOriginIndex,          // input parameter, SB index, used to point to source/reference samples
    uint64_t                        *variance32x32,
    EbAsm                         asm_type)
{


    uint32_t blockIndex;

    uint64_t meanOf8x8Blocks[64];
    uint64_t meanOf8x8SquaredValuesBlocks[64];

    uint64_t meanOf16x16Blocks[16];
    uint64_t meanOf16x16SquaredValuesBlocks[16];

    uint64_t meanOf32x32Blocks[4];
    uint64_t meanOf32x32SquaredValuesBlocks[4];

    uint64_t meanOf64x64Blocks;
    uint64_t meanOf64x64SquaredValuesBlocks;

    // (0,0)
    blockIndex = inputLumaOriginIndex;
    const uint16_t strideY = inputPaddedPicturePtr->strideY;
    if (sequence_control_set_ptr->block_mean_calc_prec == BLOCK_MEAN_PREC_FULL) {
        meanOf8x8Blocks[0] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[0] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (0,1)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[1] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[1] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (0,2)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[2] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[2] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (0,3)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[3] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[3] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (0,4)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[4] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[4] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (0,5)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[5] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[5] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (0,6)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[6] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[6] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (0,7)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[7] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[7] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (1,0)
        blockIndex = inputLumaOriginIndex + (inputPaddedPicturePtr->strideY << 3);
        meanOf8x8Blocks[8] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[8] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (1,1)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[9] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[9] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (1,2)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[10] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[10] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (1,3)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[11] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[11] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (1,4)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[12] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[12] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (1,5)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[13] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[13] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (1,6)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[14] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[14] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (1,7)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[15] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[15] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (2,0)
        blockIndex = inputLumaOriginIndex + (inputPaddedPicturePtr->strideY << 4);
        meanOf8x8Blocks[16] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[16] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (2,1)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[17] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[17] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (2,2)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[18] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[18] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (2,3)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[19] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[19] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        /// (2,4)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[20] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[20] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (2,5)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[21] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[21] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (2,6)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[22] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[22] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (2,7)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[23] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[23] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (3,0)
        blockIndex = inputLumaOriginIndex + (inputPaddedPicturePtr->strideY << 3) + (inputPaddedPicturePtr->strideY << 4);
        meanOf8x8Blocks[24] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[24] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (3,1)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[25] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[25] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (3,2)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[26] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[26] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (3,3)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[27] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[27] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (3,4)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[28] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[28] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (3,5)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[29] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[29] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (3,6)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[30] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[30] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (3,7)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[31] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[31] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (4,0)
        blockIndex = inputLumaOriginIndex + (inputPaddedPicturePtr->strideY << 5);
        meanOf8x8Blocks[32] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[32] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (4,1)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[33] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[33] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (4,2)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[34] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[34] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (4,3)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[35] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[35] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (4,4)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[36] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[36] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (4,5)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[37] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[37] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (4,6)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[38] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[38] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (4,7)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[39] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[39] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (5,0)
        blockIndex = inputLumaOriginIndex + (inputPaddedPicturePtr->strideY << 3) + (inputPaddedPicturePtr->strideY << 5);
        meanOf8x8Blocks[40] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[40] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (5,1)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[41] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[41] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (5,2)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[42] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[42] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (5,3)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[43] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[43] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (5,4)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[44] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[44] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (5,5)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[45] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[45] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (5,6)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[46] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[46] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (5,7)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[47] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[47] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (6,0)
        blockIndex = inputLumaOriginIndex + (inputPaddedPicturePtr->strideY << 4) + (inputPaddedPicturePtr->strideY << 5);
        meanOf8x8Blocks[48] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[48] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (6,1)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[49] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[49] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (6,2)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[50] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[50] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (6,3)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[51] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[51] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (6,4)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[52] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[52] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (6,5)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[53] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[53] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (6,6)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[54] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[54] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (6,7)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[55] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[55] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (7,0)
        blockIndex = inputLumaOriginIndex + (inputPaddedPicturePtr->strideY << 3) + (inputPaddedPicturePtr->strideY << 4) + (inputPaddedPicturePtr->strideY << 5);
        meanOf8x8Blocks[56] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[56] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (7,1)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[57] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[57] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (7,2)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[58] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[58] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (7,3)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[59] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[59] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (7,4)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[60] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[60] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (7,5)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[61] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[61] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (7,6)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[62] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[62] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (7,7)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[63] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[63] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

    }


    else {
        if (asm_type == ASM_AVX2) {

            ComputeIntermVarFour8x8_AVX2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY, &meanOf8x8Blocks[0], &meanOf8x8SquaredValuesBlocks[0]);

            // (0,1)
            blockIndex = blockIndex + 32;

            ComputeIntermVarFour8x8_AVX2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY, &meanOf8x8Blocks[4], &meanOf8x8SquaredValuesBlocks[4]);
            // (0,5)
            blockIndex = blockIndex + 24;

            // (1,0)
            blockIndex = inputLumaOriginIndex + (strideY << 3);

            ComputeIntermVarFour8x8_AVX2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY, &meanOf8x8Blocks[8], &meanOf8x8SquaredValuesBlocks[8]);

            // (1,1)
            blockIndex = blockIndex + 32;

            ComputeIntermVarFour8x8_AVX2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY, &meanOf8x8Blocks[12], &meanOf8x8SquaredValuesBlocks[12]);

            // (1,5)
            blockIndex = blockIndex + 24;

            // (2,0)
            blockIndex = inputLumaOriginIndex + (strideY << 4);

            ComputeIntermVarFour8x8_AVX2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY, &meanOf8x8Blocks[16], &meanOf8x8SquaredValuesBlocks[16]);

            // (2,1)
            blockIndex = blockIndex + 32;

            ComputeIntermVarFour8x8_AVX2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY, &meanOf8x8Blocks[20], &meanOf8x8SquaredValuesBlocks[20]);

            // (2,5)
            blockIndex = blockIndex + 24;

            // (3,0)
            blockIndex = inputLumaOriginIndex + (strideY << 3) + (strideY << 4);

            ComputeIntermVarFour8x8_AVX2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY, &meanOf8x8Blocks[24], &meanOf8x8SquaredValuesBlocks[24]);

            // (3,1)
            blockIndex = blockIndex + 32;

            ComputeIntermVarFour8x8_AVX2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY, &meanOf8x8Blocks[28], &meanOf8x8SquaredValuesBlocks[28]);

            // (3,5)
            blockIndex = blockIndex + 24;

            // (4,0)
            blockIndex = inputLumaOriginIndex + (strideY << 5);

            ComputeIntermVarFour8x8_AVX2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY, &meanOf8x8Blocks[32], &meanOf8x8SquaredValuesBlocks[32]);

            // (4,1)
            blockIndex = blockIndex + 32;

            ComputeIntermVarFour8x8_AVX2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY, &meanOf8x8Blocks[36], &meanOf8x8SquaredValuesBlocks[36]);

            // (4,5)
            blockIndex = blockIndex + 24;

            // (5,0)
            blockIndex = inputLumaOriginIndex + (strideY << 3) + (strideY << 5);

            ComputeIntermVarFour8x8_AVX2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY, &meanOf8x8Blocks[40], &meanOf8x8SquaredValuesBlocks[40]);

            // (5,1)
            blockIndex = blockIndex + 32;

            ComputeIntermVarFour8x8_AVX2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY, &meanOf8x8Blocks[44], &meanOf8x8SquaredValuesBlocks[44]);

            // (5,5)
            blockIndex = blockIndex + 24;

            // (6,0)
            blockIndex = inputLumaOriginIndex + (strideY << 4) + (strideY << 5);

            ComputeIntermVarFour8x8_AVX2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY, &meanOf8x8Blocks[48], &meanOf8x8SquaredValuesBlocks[48]);

            // (6,1)
            blockIndex = blockIndex + 32;

            ComputeIntermVarFour8x8_AVX2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY, &meanOf8x8Blocks[52], &meanOf8x8SquaredValuesBlocks[52]);

            // (6,5)
            blockIndex = blockIndex + 24;

            // (7,0)
            blockIndex = inputLumaOriginIndex + (strideY << 3) + (strideY << 4) + (strideY << 5);

            ComputeIntermVarFour8x8_AVX2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY, &meanOf8x8Blocks[56], &meanOf8x8SquaredValuesBlocks[56]);

            // (7,1)
            blockIndex = blockIndex + 32;

            ComputeIntermVarFour8x8_AVX2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY, &meanOf8x8Blocks[60], &meanOf8x8SquaredValuesBlocks[60]);


        }
        else {
            meanOf8x8Blocks[0] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[0] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (0,1)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[1] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[1] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (0,2)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[2] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[2] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (0,3)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[3] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[3] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (0,4)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[4] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[4] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (0,5)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[5] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[5] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (0,6)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[6] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[6] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (0,7)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[7] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[7] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (1,0)
            blockIndex = inputLumaOriginIndex + (strideY << 3);
            meanOf8x8Blocks[8] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[8] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (1,1)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[9] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[9] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (1,2)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[10] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[10] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (1,3)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[11] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[11] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (1,4)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[12] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[12] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (1,5)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[13] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[13] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (1,6)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[14] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[14] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (1,7)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[15] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[15] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (2,0)
            blockIndex = inputLumaOriginIndex + (strideY << 4);
            meanOf8x8Blocks[16] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[16] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (2,1)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[17] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[17] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (2,2)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[18] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[18] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (2,3)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[19] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[19] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            /// (2,4)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[20] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[20] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (2,5)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[21] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[21] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (2,6)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[22] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[22] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (2,7)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[23] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[23] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (3,0)
            blockIndex = inputLumaOriginIndex + (strideY << 3) + (strideY << 4);
            meanOf8x8Blocks[24] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[24] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (3,1)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[25] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[25] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (3,2)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[26] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[26] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (3,3)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[27] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[27] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (3,4)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[28] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[28] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (3,5)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[29] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[29] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (3,6)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[30] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[30] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (3,7)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[31] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[31] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (4,0)
            blockIndex = inputLumaOriginIndex + (strideY << 5);
            meanOf8x8Blocks[32] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[32] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (4,1)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[33] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[33] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (4,2)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[34] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[34] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (4,3)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[35] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[35] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (4,4)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[36] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[36] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (4,5)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[37] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[37] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (4,6)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[38] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[38] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (4,7)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[39] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[39] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (5,0)
            blockIndex = inputLumaOriginIndex + (strideY << 3) + (strideY << 5);
            meanOf8x8Blocks[40] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[40] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (5,1)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[41] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[41] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (5,2)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[42] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[42] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (5,3)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[43] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[43] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (5,4)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[44] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[44] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (5,5)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[45] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[45] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (5,6)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[46] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[46] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (5,7)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[47] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[47] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (6,0)
            blockIndex = inputLumaOriginIndex + (strideY << 4) + (strideY << 5);
            meanOf8x8Blocks[48] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[48] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (6,1)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[49] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[49] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (6,2)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[50] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[50] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (6,3)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[51] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[51] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (6,4)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[52] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[52] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (6,5)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[53] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[53] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (6,6)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[54] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[54] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (6,7)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[55] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[55] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (7,0)
            blockIndex = inputLumaOriginIndex + (strideY << 3) + (strideY << 4) + (strideY << 5);
            meanOf8x8Blocks[56] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[56] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (7,1)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[57] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[57] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (7,2)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[58] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[58] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (7,3)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[59] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[59] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (7,4)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[60] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[60] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (7,5)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[61] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[61] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (7,6)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[62] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[62] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (7,7)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[63] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[63] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);


        }
    }

    // 16x16
    meanOf16x16Blocks[0] = (meanOf8x8Blocks[0] + meanOf8x8Blocks[1] + meanOf8x8Blocks[8] + meanOf8x8Blocks[9]) >> 2;
    meanOf16x16Blocks[1] = (meanOf8x8Blocks[2] + meanOf8x8Blocks[3] + meanOf8x8Blocks[10] + meanOf8x8Blocks[11]) >> 2;
    meanOf16x16Blocks[2] = (meanOf8x8Blocks[4] + meanOf8x8Blocks[5] + meanOf8x8Blocks[12] + meanOf8x8Blocks[13]) >> 2;
    meanOf16x16Blocks[3] = (meanOf8x8Blocks[6] + meanOf8x8Blocks[7] + meanOf8x8Blocks[14] + meanOf8x8Blocks[15]) >> 2;

    meanOf16x16Blocks[4] = (meanOf8x8Blocks[16] + meanOf8x8Blocks[17] + meanOf8x8Blocks[24] + meanOf8x8Blocks[25]) >> 2;
    meanOf16x16Blocks[5] = (meanOf8x8Blocks[18] + meanOf8x8Blocks[19] + meanOf8x8Blocks[26] + meanOf8x8Blocks[27]) >> 2;
    meanOf16x16Blocks[6] = (meanOf8x8Blocks[20] + meanOf8x8Blocks[21] + meanOf8x8Blocks[28] + meanOf8x8Blocks[29]) >> 2;
    meanOf16x16Blocks[7] = (meanOf8x8Blocks[22] + meanOf8x8Blocks[23] + meanOf8x8Blocks[30] + meanOf8x8Blocks[31]) >> 2;

    meanOf16x16Blocks[8] = (meanOf8x8Blocks[32] + meanOf8x8Blocks[33] + meanOf8x8Blocks[40] + meanOf8x8Blocks[41]) >> 2;
    meanOf16x16Blocks[9] = (meanOf8x8Blocks[34] + meanOf8x8Blocks[35] + meanOf8x8Blocks[42] + meanOf8x8Blocks[43]) >> 2;
    meanOf16x16Blocks[10] = (meanOf8x8Blocks[36] + meanOf8x8Blocks[37] + meanOf8x8Blocks[44] + meanOf8x8Blocks[45]) >> 2;
    meanOf16x16Blocks[11] = (meanOf8x8Blocks[38] + meanOf8x8Blocks[39] + meanOf8x8Blocks[46] + meanOf8x8Blocks[47]) >> 2;

    meanOf16x16Blocks[12] = (meanOf8x8Blocks[48] + meanOf8x8Blocks[49] + meanOf8x8Blocks[56] + meanOf8x8Blocks[57]) >> 2;
    meanOf16x16Blocks[13] = (meanOf8x8Blocks[50] + meanOf8x8Blocks[51] + meanOf8x8Blocks[58] + meanOf8x8Blocks[59]) >> 2;
    meanOf16x16Blocks[14] = (meanOf8x8Blocks[52] + meanOf8x8Blocks[53] + meanOf8x8Blocks[60] + meanOf8x8Blocks[61]) >> 2;
    meanOf16x16Blocks[15] = (meanOf8x8Blocks[54] + meanOf8x8Blocks[55] + meanOf8x8Blocks[62] + meanOf8x8Blocks[63]) >> 2;

    meanOf16x16SquaredValuesBlocks[0] = (meanOf8x8SquaredValuesBlocks[0] + meanOf8x8SquaredValuesBlocks[1] + meanOf8x8SquaredValuesBlocks[8] + meanOf8x8SquaredValuesBlocks[9]) >> 2;
    meanOf16x16SquaredValuesBlocks[1] = (meanOf8x8SquaredValuesBlocks[2] + meanOf8x8SquaredValuesBlocks[3] + meanOf8x8SquaredValuesBlocks[10] + meanOf8x8SquaredValuesBlocks[11]) >> 2;
    meanOf16x16SquaredValuesBlocks[2] = (meanOf8x8SquaredValuesBlocks[4] + meanOf8x8SquaredValuesBlocks[5] + meanOf8x8SquaredValuesBlocks[12] + meanOf8x8SquaredValuesBlocks[13]) >> 2;
    meanOf16x16SquaredValuesBlocks[3] = (meanOf8x8SquaredValuesBlocks[6] + meanOf8x8SquaredValuesBlocks[7] + meanOf8x8SquaredValuesBlocks[14] + meanOf8x8SquaredValuesBlocks[15]) >> 2;

    meanOf16x16SquaredValuesBlocks[4] = (meanOf8x8SquaredValuesBlocks[16] + meanOf8x8SquaredValuesBlocks[17] + meanOf8x8SquaredValuesBlocks[24] + meanOf8x8SquaredValuesBlocks[25]) >> 2;
    meanOf16x16SquaredValuesBlocks[5] = (meanOf8x8SquaredValuesBlocks[18] + meanOf8x8SquaredValuesBlocks[19] + meanOf8x8SquaredValuesBlocks[26] + meanOf8x8SquaredValuesBlocks[27]) >> 2;
    meanOf16x16SquaredValuesBlocks[6] = (meanOf8x8SquaredValuesBlocks[20] + meanOf8x8SquaredValuesBlocks[21] + meanOf8x8SquaredValuesBlocks[28] + meanOf8x8SquaredValuesBlocks[29]) >> 2;
    meanOf16x16SquaredValuesBlocks[7] = (meanOf8x8SquaredValuesBlocks[22] + meanOf8x8SquaredValuesBlocks[23] + meanOf8x8SquaredValuesBlocks[30] + meanOf8x8SquaredValuesBlocks[31]) >> 2;

    meanOf16x16SquaredValuesBlocks[8] = (meanOf8x8SquaredValuesBlocks[32] + meanOf8x8SquaredValuesBlocks[33] + meanOf8x8SquaredValuesBlocks[40] + meanOf8x8SquaredValuesBlocks[41]) >> 2;
    meanOf16x16SquaredValuesBlocks[9] = (meanOf8x8SquaredValuesBlocks[34] + meanOf8x8SquaredValuesBlocks[35] + meanOf8x8SquaredValuesBlocks[42] + meanOf8x8SquaredValuesBlocks[43]) >> 2;
    meanOf16x16SquaredValuesBlocks[10] = (meanOf8x8SquaredValuesBlocks[36] + meanOf8x8SquaredValuesBlocks[37] + meanOf8x8SquaredValuesBlocks[44] + meanOf8x8SquaredValuesBlocks[45]) >> 2;
    meanOf16x16SquaredValuesBlocks[11] = (meanOf8x8SquaredValuesBlocks[38] + meanOf8x8SquaredValuesBlocks[39] + meanOf8x8SquaredValuesBlocks[46] + meanOf8x8SquaredValuesBlocks[47]) >> 2;

    meanOf16x16SquaredValuesBlocks[12] = (meanOf8x8SquaredValuesBlocks[48] + meanOf8x8SquaredValuesBlocks[49] + meanOf8x8SquaredValuesBlocks[56] + meanOf8x8SquaredValuesBlocks[57]) >> 2;
    meanOf16x16SquaredValuesBlocks[13] = (meanOf8x8SquaredValuesBlocks[50] + meanOf8x8SquaredValuesBlocks[51] + meanOf8x8SquaredValuesBlocks[58] + meanOf8x8SquaredValuesBlocks[59]) >> 2;
    meanOf16x16SquaredValuesBlocks[14] = (meanOf8x8SquaredValuesBlocks[52] + meanOf8x8SquaredValuesBlocks[53] + meanOf8x8SquaredValuesBlocks[60] + meanOf8x8SquaredValuesBlocks[61]) >> 2;
    meanOf16x16SquaredValuesBlocks[15] = (meanOf8x8SquaredValuesBlocks[54] + meanOf8x8SquaredValuesBlocks[55] + meanOf8x8SquaredValuesBlocks[62] + meanOf8x8SquaredValuesBlocks[63]) >> 2;

    // 32x32
    meanOf32x32Blocks[0] = (meanOf16x16Blocks[0] + meanOf16x16Blocks[1] + meanOf16x16Blocks[4] + meanOf16x16Blocks[5]) >> 2;
    meanOf32x32Blocks[1] = (meanOf16x16Blocks[2] + meanOf16x16Blocks[3] + meanOf16x16Blocks[6] + meanOf16x16Blocks[7]) >> 2;
    meanOf32x32Blocks[2] = (meanOf16x16Blocks[8] + meanOf16x16Blocks[9] + meanOf16x16Blocks[12] + meanOf16x16Blocks[13]) >> 2;
    meanOf32x32Blocks[3] = (meanOf16x16Blocks[10] + meanOf16x16Blocks[11] + meanOf16x16Blocks[14] + meanOf16x16Blocks[15]) >> 2;

    meanOf32x32SquaredValuesBlocks[0] = (meanOf16x16SquaredValuesBlocks[0] + meanOf16x16SquaredValuesBlocks[1] + meanOf16x16SquaredValuesBlocks[4] + meanOf16x16SquaredValuesBlocks[5]) >> 2;
    meanOf32x32SquaredValuesBlocks[1] = (meanOf16x16SquaredValuesBlocks[2] + meanOf16x16SquaredValuesBlocks[3] + meanOf16x16SquaredValuesBlocks[6] + meanOf16x16SquaredValuesBlocks[7]) >> 2;
    meanOf32x32SquaredValuesBlocks[2] = (meanOf16x16SquaredValuesBlocks[8] + meanOf16x16SquaredValuesBlocks[9] + meanOf16x16SquaredValuesBlocks[12] + meanOf16x16SquaredValuesBlocks[13]) >> 2;
    meanOf32x32SquaredValuesBlocks[3] = (meanOf16x16SquaredValuesBlocks[10] + meanOf16x16SquaredValuesBlocks[11] + meanOf16x16SquaredValuesBlocks[14] + meanOf16x16SquaredValuesBlocks[15]) >> 2;


    variance32x32[0] = meanOf32x32SquaredValuesBlocks[0] - (meanOf32x32Blocks[0] * meanOf32x32Blocks[0]);
    variance32x32[1] = meanOf32x32SquaredValuesBlocks[1] - (meanOf32x32Blocks[1] * meanOf32x32Blocks[1]);
    variance32x32[2] = meanOf32x32SquaredValuesBlocks[2] - (meanOf32x32Blocks[2] * meanOf32x32Blocks[2]);
    variance32x32[3] = meanOf32x32SquaredValuesBlocks[3] - (meanOf32x32Blocks[3] * meanOf32x32Blocks[3]);


    // 64x64
    meanOf64x64Blocks = (meanOf32x32Blocks[0] + meanOf32x32Blocks[1] + meanOf32x32Blocks[2] + meanOf32x32Blocks[3]) >> 2;
    meanOf64x64SquaredValuesBlocks = (meanOf32x32SquaredValuesBlocks[0] + meanOf32x32SquaredValuesBlocks[1] + meanOf32x32SquaredValuesBlocks[2] + meanOf32x32SquaredValuesBlocks[3]) >> 2;

    return (meanOf64x64SquaredValuesBlocks - (meanOf64x64Blocks * meanOf64x64Blocks));
}

uint8_t  getFilteredTypes(uint8_t  *ptr,
    uint32_t  stride,
    uint8_t   filterType)
{
    uint8_t *p = ptr - 1 - stride;

    uint32_t a = 0;

    if (filterType == 0) {

        //Luma
        a = (p[1] +
            p[0 + stride] + 4 * p[1 + stride] + p[2 + stride] +
            p[1 + 2 * stride]) / 8;

    }
    else if (filterType == 1) {
        a = (2 * p[1] +
            2 * p[0 + stride] + 4 * p[1 + stride] + 2 * p[2 + stride] +
            2 * p[1 + 2 * stride]);

        a = (((uint32_t)((a * 2730) >> 14) + 1) >> 1) & 0xFFFF;

        //fixed point version of a=a/12 to mimic x86 instruction _mm256_mulhrs_epi16;
        //a= (a*2730)>>15;
    }
    else if (filterType == 2) {


        a = (4 * p[1] +
            4 * p[0 + stride] + 4 * p[1 + stride] + 4 * p[2 + stride] +
            4 * p[1 + 2 * stride]) / 20;
    }
    else if (filterType == 3) {

        a = (1 * p[0] + 1 * p[1] + 1 * p[2] +
            1 * p[0 + stride] + 4 * p[1 + stride] + 1 * p[2 + stride] +
            1 * p[0 + 2 * stride] + 1 * p[1 + 2 * stride] + 1 * p[2 + 2 * stride]) / 12;


    }
    else if (filterType == 4) {

        //gaussian matrix(Chroma)
        a = (1 * p[0] + 2 * p[1] + 1 * p[2] +
            2 * p[0 + stride] + 4 * p[1 + stride] + 2 * p[2 + stride] +
            1 * p[0 + 2 * stride] + 2 * p[1 + 2 * stride] + 1 * p[2 + 2 * stride]) / 16;

    }
    else if (filterType == 5) {

        a = (2 * p[0] + 2 * p[1] + 2 * p[2] +
            2 * p[0 + stride] + 4 * p[1 + stride] + 2 * p[2 + stride] +
            2 * p[0 + 2 * stride] + 2 * p[1 + 2 * stride] + 2 * p[2 + 2 * stride]) / 20;

    }
    else if (filterType == 6) {

        a = (4 * p[0] + 4 * p[1] + 4 * p[2] +
            4 * p[0 + stride] + 4 * p[1 + stride] + 4 * p[2 + stride] +
            4 * p[0 + 2 * stride] + 4 * p[1 + 2 * stride] + 4 * p[2 + 2 * stride]) / 36;

    }


    return  (uint8_t)CLIP3EQ(0, 255, a);

}


/*******************************************
* noiseExtractLumaStrong
*  strong filter Luma.
*******************************************/
void noiseExtractLumaStrong(
    EbPictureBufferDesc_t       *inputPicturePtr,
    EbPictureBufferDesc_t       *denoisedPicturePtr,
    uint32_t                       sb_origin_y
    , uint32_t                       sb_origin_x
)
{
    uint32_t  ii, jj;
    uint32_t  picHeight, sb_height;
    uint32_t  picWidth;
    uint32_t  inputOriginIndex;
    uint32_t  inputOriginIndexPad;

    uint8_t *ptrIn;
    uint32_t strideIn;
    uint8_t *ptrDenoised;

    uint32_t strideOut;
    uint32_t idx = (sb_origin_x + BLOCK_SIZE_64 > inputPicturePtr->width) ? sb_origin_x : 0;

    //Luma
    {
        picHeight = inputPicturePtr->height;
        picWidth = inputPicturePtr->width;
        sb_height = MIN(BLOCK_SIZE_64, picHeight - sb_origin_y);

        strideIn = inputPicturePtr->strideY;
        inputOriginIndex = inputPicturePtr->origin_x + (inputPicturePtr->origin_y + sb_origin_y)* inputPicturePtr->strideY;
        ptrIn = &(inputPicturePtr->bufferY[inputOriginIndex]);

        inputOriginIndexPad = denoisedPicturePtr->origin_x + (denoisedPicturePtr->origin_y + sb_origin_y) * denoisedPicturePtr->strideY;
        strideOut = denoisedPicturePtr->strideY;
        ptrDenoised = &(denoisedPicturePtr->bufferY[inputOriginIndexPad]);

        for (jj = 0; jj < sb_height; jj++) {
            for (ii = idx; ii < picWidth; ii++) {

                if ((jj > 0 || sb_origin_y > 0) && (jj < sb_height - 1 || sb_origin_y + sb_height < picHeight) && ii > 0 && ii < picWidth - 1) {

                    ptrDenoised[ii + jj * strideOut] = getFilteredTypes(&ptrIn[ii + jj * strideIn], strideIn, 4);

                }
                else {
                    ptrDenoised[ii + jj * strideOut] = ptrIn[ii + jj * strideIn];

                }

            }
        }
    }

}
/*******************************************
* noiseExtractChromaStrong
*  strong filter chroma.
*******************************************/
void noiseExtractChromaStrong(
    EbPictureBufferDesc_t       *inputPicturePtr,
    EbPictureBufferDesc_t       *denoisedPicturePtr,
    uint32_t                       sb_origin_y
    , uint32_t                       sb_origin_x
)
{
    uint32_t  ii, jj;
    uint32_t  picHeight, sb_height;
    uint32_t  picWidth;
    uint32_t  inputOriginIndex;
    uint32_t  inputOriginIndexPad;

    uint8_t *ptrIn;
    uint32_t strideIn;
    uint8_t *ptrDenoised;

    uint32_t strideOut;
    uint32_t idx = (sb_origin_x + BLOCK_SIZE_64 > inputPicturePtr->width) ? sb_origin_x : 0;

    //Cb
    {
        picHeight = inputPicturePtr->height / 2;
        picWidth = inputPicturePtr->width / 2;
        sb_height = MIN(BLOCK_SIZE_64 / 2, picHeight - sb_origin_y);

        strideIn = inputPicturePtr->strideCb;
        inputOriginIndex = inputPicturePtr->origin_x / 2 + (inputPicturePtr->origin_y / 2 + sb_origin_y)  * inputPicturePtr->strideCb;
        ptrIn = &(inputPicturePtr->bufferCb[inputOriginIndex]);

        inputOriginIndexPad = denoisedPicturePtr->origin_x / 2 + (denoisedPicturePtr->origin_y / 2 + sb_origin_y)  * denoisedPicturePtr->strideCb;
        strideOut = denoisedPicturePtr->strideCb;
        ptrDenoised = &(denoisedPicturePtr->bufferCb[inputOriginIndexPad]);


        for (jj = 0; jj < sb_height; jj++) {
            for (ii = idx; ii < picWidth; ii++) {


                if ((jj > 0 || sb_origin_y > 0) && (jj < sb_height - 1 || (sb_origin_y + sb_height) < picHeight) && ii > 0 && ii < picWidth - 1) {
                    ptrDenoised[ii + jj * strideOut] = getFilteredTypes(&ptrIn[ii + jj * strideIn], strideIn, 6);
                }
                else {
                    ptrDenoised[ii + jj * strideOut] = ptrIn[ii + jj * strideIn];
                }

            }
        }
    }

    //Cr
    {
        picHeight = inputPicturePtr->height / 2;
        picWidth = inputPicturePtr->width / 2;
        sb_height = MIN(BLOCK_SIZE_64 / 2, picHeight - sb_origin_y);

        strideIn = inputPicturePtr->strideCr;
        inputOriginIndex = inputPicturePtr->origin_x / 2 + (inputPicturePtr->origin_y / 2 + sb_origin_y)  * inputPicturePtr->strideCr;
        ptrIn = &(inputPicturePtr->bufferCr[inputOriginIndex]);

        inputOriginIndexPad = denoisedPicturePtr->origin_x / 2 + (denoisedPicturePtr->origin_y / 2 + sb_origin_y)  * denoisedPicturePtr->strideCr;
        strideOut = denoisedPicturePtr->strideCr;
        ptrDenoised = &(denoisedPicturePtr->bufferCr[inputOriginIndexPad]);


        for (jj = 0; jj < sb_height; jj++) {
            for (ii = idx; ii < picWidth; ii++) {

                if ((jj > 0 || sb_origin_y > 0) && (jj < sb_height - 1 || (sb_origin_y + sb_height) < picHeight) && ii > 0 && ii < picWidth - 1) {
                    ptrDenoised[ii + jj * strideOut] = getFilteredTypes(&ptrIn[ii + jj * strideIn], strideIn, 6);
                }
                else {
                    ptrDenoised[ii + jj * strideOut] = ptrIn[ii + jj * strideIn];
                }

            }
        }
    }

}

/*******************************************
* noiseExtractChromaWeak
*  weak filter chroma.
*******************************************/
void noiseExtractChromaWeak(
    EbPictureBufferDesc_t       *inputPicturePtr,
    EbPictureBufferDesc_t       *denoisedPicturePtr,
    uint32_t                       sb_origin_y
    , uint32_t                       sb_origin_x
)
{
    uint32_t  ii, jj;
    uint32_t  picHeight, sb_height;
    uint32_t  picWidth;
    uint32_t  inputOriginIndex;
    uint32_t  inputOriginIndexPad;

    uint8_t *ptrIn;
    uint32_t strideIn;
    uint8_t *ptrDenoised;

    uint32_t strideOut;

    uint32_t idx = (sb_origin_x + BLOCK_SIZE_64 > inputPicturePtr->width) ? sb_origin_x : 0;

    //Cb
    {
        picHeight = inputPicturePtr->height / 2;
        picWidth = inputPicturePtr->width / 2;

        sb_height = MIN(BLOCK_SIZE_64 / 2, picHeight - sb_origin_y);

        strideIn = inputPicturePtr->strideCb;
        inputOriginIndex = inputPicturePtr->origin_x / 2 + (inputPicturePtr->origin_y / 2 + sb_origin_y)* inputPicturePtr->strideCb;
        ptrIn = &(inputPicturePtr->bufferCb[inputOriginIndex]);

        inputOriginIndexPad = denoisedPicturePtr->origin_x / 2 + (denoisedPicturePtr->origin_y / 2 + sb_origin_y)* denoisedPicturePtr->strideCb;
        strideOut = denoisedPicturePtr->strideCb;
        ptrDenoised = &(denoisedPicturePtr->bufferCb[inputOriginIndexPad]);


        for (jj = 0; jj < sb_height; jj++) {
            for (ii = idx; ii < picWidth; ii++) {


                if ((jj > 0 || sb_origin_y > 0) && (jj < sb_height - 1 || (sb_origin_y + sb_height) < picHeight) && ii > 0 && ii < picWidth - 1) {
                    ptrDenoised[ii + jj * strideOut] = getFilteredTypes(&ptrIn[ii + jj * strideIn], strideIn, 4);
                }
                else {
                    ptrDenoised[ii + jj * strideOut] = ptrIn[ii + jj * strideIn];
                }

            }
        }
    }

    //Cr
    {
        picHeight = inputPicturePtr->height / 2;
        picWidth = inputPicturePtr->width / 2;
        sb_height = MIN(BLOCK_SIZE_64 / 2, picHeight - sb_origin_y);

        strideIn = inputPicturePtr->strideCr;
        inputOriginIndex = inputPicturePtr->origin_x / 2 + (inputPicturePtr->origin_y / 2 + sb_origin_y)* inputPicturePtr->strideCr;
        ptrIn = &(inputPicturePtr->bufferCr[inputOriginIndex]);

        inputOriginIndexPad = denoisedPicturePtr->origin_x / 2 + (denoisedPicturePtr->origin_y / 2 + sb_origin_y)* denoisedPicturePtr->strideCr;
        strideOut = denoisedPicturePtr->strideCr;
        ptrDenoised = &(denoisedPicturePtr->bufferCr[inputOriginIndexPad]);


        for (jj = 0; jj < sb_height; jj++) {
            for (ii = idx; ii < picWidth; ii++) {

                if ((jj > 0 || sb_origin_y > 0) && (jj < sb_height - 1 || (sb_origin_y + sb_height) < picHeight) && ii > 0 && ii < picWidth - 1) {
                    ptrDenoised[ii + jj * strideOut] = getFilteredTypes(&ptrIn[ii + jj * strideIn], strideIn, 4);
                }
                else {
                    ptrDenoised[ii + jj * strideOut] = ptrIn[ii + jj * strideIn];
                }

            }
        }
    }

}

/*******************************************
* noiseExtractLumaWeak
*  weak filter Luma and store noise.
*******************************************/
void noiseExtractLumaWeak(
    EbPictureBufferDesc_t       *inputPicturePtr,
    EbPictureBufferDesc_t       *denoisedPicturePtr,
    EbPictureBufferDesc_t       *noisePicturePtr,
    uint32_t                       sb_origin_y
    , uint32_t                         sb_origin_x
)
{
    uint32_t  ii, jj;
    uint32_t  picHeight, sb_height;
    uint32_t  picWidth;
    uint32_t  inputOriginIndex;
    uint32_t  inputOriginIndexPad;
    uint32_t  noiseOriginIndex;

    uint8_t *ptrIn;
    uint32_t strideIn;
    uint8_t *ptrDenoised;

    uint8_t *ptrNoise;
    uint32_t strideOut;

    uint32_t idx = (sb_origin_x + BLOCK_SIZE_64 > inputPicturePtr->width) ? sb_origin_x : 0;

    //Luma
    {
        picHeight = inputPicturePtr->height;
        picWidth = inputPicturePtr->width;
        sb_height = MIN(BLOCK_SIZE_64, picHeight - sb_origin_y);

        strideIn = inputPicturePtr->strideY;
        inputOriginIndex = inputPicturePtr->origin_x + (inputPicturePtr->origin_y + sb_origin_y) * inputPicturePtr->strideY;
        ptrIn = &(inputPicturePtr->bufferY[inputOriginIndex]);

        inputOriginIndexPad = denoisedPicturePtr->origin_x + (denoisedPicturePtr->origin_y + sb_origin_y) * denoisedPicturePtr->strideY;
        strideOut = denoisedPicturePtr->strideY;
        ptrDenoised = &(denoisedPicturePtr->bufferY[inputOriginIndexPad]);

        noiseOriginIndex = noisePicturePtr->origin_x + noisePicturePtr->origin_y * noisePicturePtr->strideY;
        ptrNoise = &(noisePicturePtr->bufferY[noiseOriginIndex]);


        for (jj = 0; jj < sb_height; jj++) {
            for (ii = idx; ii < picWidth; ii++) {

                if ((jj > 0 || sb_origin_y > 0) && (jj < sb_height - 1 || sb_origin_y + sb_height < picHeight) && ii > 0 && ii < picWidth - 1) {

                    ptrDenoised[ii + jj * strideOut] = getFilteredTypes(&ptrIn[ii + jj * strideIn], strideIn, 0);
                    ptrNoise[ii + jj * strideOut] = CLIP3EQ(0, 255, ptrIn[ii + jj * strideIn] - ptrDenoised[ii + jj * strideOut]);

                }
                else {
                    ptrDenoised[ii + jj * strideOut] = ptrIn[ii + jj * strideIn];
                    ptrNoise[ii + jj * strideOut] = 0;
                }

            }
        }
    }

}

void noiseExtractLumaWeakLcu(
    EbPictureBufferDesc_t       *inputPicturePtr,
    EbPictureBufferDesc_t       *denoisedPicturePtr,
    EbPictureBufferDesc_t       *noisePicturePtr,
    uint32_t                       sb_origin_y
    , uint32_t                         sb_origin_x
)
{
    uint32_t  ii, jj;
    uint32_t  picHeight, sb_height;
    uint32_t  picWidth, sb_width;
    uint32_t  inputOriginIndex;
    uint32_t  inputOriginIndexPad;
    uint32_t  noiseOriginIndex;

    uint8_t *ptrIn;
    uint32_t strideIn;
    uint8_t *ptrDenoised;

    uint8_t *ptrNoise;
    uint32_t strideOut;

    uint32_t idx = (sb_origin_x + BLOCK_SIZE_64 > inputPicturePtr->width) ? sb_origin_x : 0;

    //Luma
    {
        picHeight = inputPicturePtr->height;
        picWidth = inputPicturePtr->width;
        sb_height = MIN(BLOCK_SIZE_64, picHeight - sb_origin_y);
        sb_width = MIN(BLOCK_SIZE_64, picWidth - sb_origin_x);

        strideIn = inputPicturePtr->strideY;
        inputOriginIndex = inputPicturePtr->origin_x + sb_origin_x + (inputPicturePtr->origin_y + sb_origin_y) * inputPicturePtr->strideY;
        ptrIn = &(inputPicturePtr->bufferY[inputOriginIndex]);

        inputOriginIndexPad = denoisedPicturePtr->origin_x + sb_origin_x + (denoisedPicturePtr->origin_y + sb_origin_y) * denoisedPicturePtr->strideY;
        strideOut = denoisedPicturePtr->strideY;
        ptrDenoised = &(denoisedPicturePtr->bufferY[inputOriginIndexPad]);

        noiseOriginIndex = noisePicturePtr->origin_x + sb_origin_x + noisePicturePtr->origin_y * noisePicturePtr->strideY;
        ptrNoise = &(noisePicturePtr->bufferY[noiseOriginIndex]);


        for (jj = 0; jj < sb_height; jj++) {
            for (ii = idx; ii < sb_width; ii++) {

                if ((jj > 0 || sb_origin_y > 0) && (jj < sb_height - 1 || sb_origin_y + sb_height < picHeight) && (ii > 0 || sb_origin_x > 0) && (ii + sb_origin_x) < picWidth - 1/* & ii < sb_width - 1*/) {

                    ptrDenoised[ii + jj * strideOut] = getFilteredTypes(&ptrIn[ii + jj * strideIn], strideIn, 0);
                    ptrNoise[ii + jj * strideOut] = CLIP3EQ(0, 255, ptrIn[ii + jj * strideIn] - ptrDenoised[ii + jj * strideOut]);

                }
                else {
                    ptrDenoised[ii + jj * strideOut] = ptrIn[ii + jj * strideIn];
                    ptrNoise[ii + jj * strideOut] = 0;
                }

            }
        }
    }

}

EbErrorType ZeroOutChromaBlockMean(
    PictureParentControlSet_t   *picture_control_set_ptr,          // input parameter, Picture Control Set Ptr
    uint32_t                       lcuCodingOrder                // input parameter, SB address
)
{

    EbErrorType return_error = EB_ErrorNone;
    // 16x16 mean
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_0] = 0;
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_1] = 0;
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_2] = 0;
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_3] = 0;
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_4] = 0;
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_5] = 0;
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_6] = 0;
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_7] = 0;
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_8] = 0;
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_9] = 0;
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_10] = 0;
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_11] = 0;
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_12] = 0;
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_13] = 0;
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_14] = 0;
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_15] = 0;

    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_0] = 0;
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_1] = 0;
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_2] = 0;
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_3] = 0;
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_4] = 0;
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_5] = 0;
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_6] = 0;
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_7] = 0;
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_8] = 0;
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_9] = 0;
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_10] = 0;
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_11] = 0;
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_12] = 0;
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_13] = 0;
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_14] = 0;
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_15] = 0;

    // 32x32 mean
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_32x32_0] = 0;
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_32x32_1] = 0;
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_32x32_2] = 0;
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_32x32_3] = 0;

    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_32x32_0] = 0;
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_32x32_1] = 0;
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_32x32_2] = 0;
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_32x32_3] = 0;

    // 64x64 mean
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_64x64] = 0;
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_64x64] = 0;

    return return_error;

}
/*******************************************
* ComputeChromaBlockMean
*   computes the chroma block mean for 64x64, 32x32 and 16x16 CUs inside the tree block
*******************************************/
EbErrorType ComputeChromaBlockMean(
    SequenceControlSet_t        *sequence_control_set_ptr,
    PictureParentControlSet_t   *picture_control_set_ptr,          // input parameter, Picture Control Set Ptr
    EbPictureBufferDesc_t       *inputPaddedPicturePtr,         // input parameter, Input Padded Picture
    uint32_t                       lcuCodingOrder,                // input parameter, SB address
    uint32_t                       inputCbOriginIndex,            // input parameter, SB index, used to point to source/reference samples
    uint32_t                       inputCrOriginIndex,            // input parameter, SB index, used to point to source/reference samples
    EbAsm                         asm_type)
{

    EbErrorType return_error = EB_ErrorNone;

    uint32_t cbBlockIndex, crBlockIndex;

    uint64_t cbMeanOf16x16Blocks[16];
    uint64_t crMeanOf16x16Blocks[16];

    uint64_t cbMeanOf32x32Blocks[4];
    uint64_t crMeanOf32x32Blocks[4];

    uint64_t cbMeanOf64x64Blocks;
    uint64_t crMeanOf64x64Blocks;


    // (0,0) 16x16 block
    cbBlockIndex = inputCbOriginIndex;
    crBlockIndex = inputCrOriginIndex;
    if (sequence_control_set_ptr->block_mean_calc_prec == BLOCK_MEAN_PREC_FULL) {
        cbMeanOf16x16Blocks[0] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferCb[cbBlockIndex]), inputPaddedPicturePtr->strideCb, 8, 8);
        crMeanOf16x16Blocks[0] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferCr[crBlockIndex]), inputPaddedPicturePtr->strideCr, 8, 8);

        // (0,1)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[1] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferCb[cbBlockIndex]), inputPaddedPicturePtr->strideCb, 8, 8);
        crMeanOf16x16Blocks[1] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferCr[crBlockIndex]), inputPaddedPicturePtr->strideCr, 8, 8);

        // (0,2)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[2] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferCb[cbBlockIndex]), inputPaddedPicturePtr->strideCb, 8, 8);
        crMeanOf16x16Blocks[2] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferCr[crBlockIndex]), inputPaddedPicturePtr->strideCr, 8, 8);

        // (0,3)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[3] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferCb[cbBlockIndex]), inputPaddedPicturePtr->strideCb, 8, 8);
        crMeanOf16x16Blocks[3] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferCr[crBlockIndex]), inputPaddedPicturePtr->strideCr, 8, 8);

        // (1,0)
        cbBlockIndex = inputCbOriginIndex + (inputPaddedPicturePtr->strideCb << 3);
        crBlockIndex = inputCrOriginIndex + (inputPaddedPicturePtr->strideCr << 3);
        cbMeanOf16x16Blocks[4] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferCb[cbBlockIndex]), inputPaddedPicturePtr->strideCb, 8, 8);
        crMeanOf16x16Blocks[4] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferCr[crBlockIndex]), inputPaddedPicturePtr->strideCr, 8, 8);

        // (1,1)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[5] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferCb[cbBlockIndex]), inputPaddedPicturePtr->strideCb, 8, 8);
        crMeanOf16x16Blocks[5] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferCr[crBlockIndex]), inputPaddedPicturePtr->strideCr, 8, 8);

        // (1,2)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[6] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferCb[cbBlockIndex]), inputPaddedPicturePtr->strideCb, 8, 8);
        crMeanOf16x16Blocks[6] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferCr[crBlockIndex]), inputPaddedPicturePtr->strideCr, 8, 8);

        // (1,3)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[7] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferCb[cbBlockIndex]), inputPaddedPicturePtr->strideCb, 8, 8);
        crMeanOf16x16Blocks[7] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferCr[crBlockIndex]), inputPaddedPicturePtr->strideCr, 8, 8);

        // (2,0)
        cbBlockIndex = inputCbOriginIndex + (inputPaddedPicturePtr->strideCb << 4);
        crBlockIndex = inputCrOriginIndex + (inputPaddedPicturePtr->strideCr << 4);
        cbMeanOf16x16Blocks[8] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferCb[cbBlockIndex]), inputPaddedPicturePtr->strideCb, 8, 8);
        crMeanOf16x16Blocks[8] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferCr[crBlockIndex]), inputPaddedPicturePtr->strideCr, 8, 8);

        // (2,1)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[9] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferCb[cbBlockIndex]), inputPaddedPicturePtr->strideCb, 8, 8);
        crMeanOf16x16Blocks[9] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferCr[crBlockIndex]), inputPaddedPicturePtr->strideCr, 8, 8);

        // (2,2)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[10] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferCb[cbBlockIndex]), inputPaddedPicturePtr->strideCb, 8, 8);
        crMeanOf16x16Blocks[10] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferCr[crBlockIndex]), inputPaddedPicturePtr->strideCr, 8, 8);

        // (2,3)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[11] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferCb[cbBlockIndex]), inputPaddedPicturePtr->strideCb, 8, 8);
        crMeanOf16x16Blocks[11] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferCr[crBlockIndex]), inputPaddedPicturePtr->strideCr, 8, 8);

        // (3,0)
        cbBlockIndex = inputCbOriginIndex + (inputPaddedPicturePtr->strideCb * 24);
        crBlockIndex = inputCrOriginIndex + (inputPaddedPicturePtr->strideCr * 24);
        cbMeanOf16x16Blocks[12] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferCb[cbBlockIndex]), inputPaddedPicturePtr->strideCb, 8, 8);
        crMeanOf16x16Blocks[12] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferCr[crBlockIndex]), inputPaddedPicturePtr->strideCr, 8, 8);

        // (3,1)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[13] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferCb[cbBlockIndex]), inputPaddedPicturePtr->strideCb, 8, 8);
        crMeanOf16x16Blocks[13] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferCr[crBlockIndex]), inputPaddedPicturePtr->strideCr, 8, 8);

        // (3,2)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[14] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferCb[cbBlockIndex]), inputPaddedPicturePtr->strideCb, 8, 8);
        crMeanOf16x16Blocks[14] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferCr[crBlockIndex]), inputPaddedPicturePtr->strideCr, 8, 8);

        // (3,3)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[15] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferCb[cbBlockIndex]), inputPaddedPicturePtr->strideCb, 8, 8);
        crMeanOf16x16Blocks[15] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferCr[crBlockIndex]), inputPaddedPicturePtr->strideCr, 8, 8);

    }
    else {
        (void)(asm_type);
        const uint16_t strideCb = inputPaddedPicturePtr->strideCb;
        const uint16_t strideCr = inputPaddedPicturePtr->strideCr;

        cbMeanOf16x16Blocks[0] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferCb[cbBlockIndex]), strideCb);
        crMeanOf16x16Blocks[0] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferCr[crBlockIndex]), strideCr);

        // (0,1)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[1] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferCb[cbBlockIndex]), strideCb);
        crMeanOf16x16Blocks[1] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferCr[crBlockIndex]), strideCr);

        // (0,2)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[2] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferCb[cbBlockIndex]), strideCb);
        crMeanOf16x16Blocks[2] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferCr[crBlockIndex]), strideCr);

        // (0,3)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[3] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferCb[cbBlockIndex]), strideCb);
        crMeanOf16x16Blocks[3] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferCr[crBlockIndex]), strideCr);

        // (1,0)
        cbBlockIndex = inputCbOriginIndex + (strideCb << 3);
        crBlockIndex = inputCrOriginIndex + (strideCr << 3);
        cbMeanOf16x16Blocks[4] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferCb[cbBlockIndex]), strideCb);
        crMeanOf16x16Blocks[4] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferCr[crBlockIndex]), strideCr);

        // (1,1)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[5] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferCb[cbBlockIndex]), strideCb);
        crMeanOf16x16Blocks[5] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferCr[crBlockIndex]), strideCr);

        // (1,2)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[6] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferCb[cbBlockIndex]), strideCb);
        crMeanOf16x16Blocks[6] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferCr[crBlockIndex]), strideCr);

        // (1,3)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[7] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferCb[cbBlockIndex]), strideCb);
        crMeanOf16x16Blocks[7] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferCr[crBlockIndex]), strideCr);

        // (2,0)
        cbBlockIndex = inputCbOriginIndex + (strideCb << 4);
        crBlockIndex = inputCrOriginIndex + (strideCr << 4);
        cbMeanOf16x16Blocks[8] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferCb[cbBlockIndex]), strideCb);
        crMeanOf16x16Blocks[8] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferCr[crBlockIndex]), strideCr);

        // (2,1)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[9] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferCb[cbBlockIndex]), strideCb);
        crMeanOf16x16Blocks[9] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferCr[crBlockIndex]), strideCr);

        // (2,2)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[10] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferCb[cbBlockIndex]), strideCb);
        crMeanOf16x16Blocks[10] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferCr[crBlockIndex]), strideCr);

        // (2,3)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[11] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferCb[cbBlockIndex]), strideCb);
        crMeanOf16x16Blocks[11] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferCr[crBlockIndex]), strideCr);

        // (3,0)
        cbBlockIndex = inputCbOriginIndex + (strideCb * 24);
        crBlockIndex = inputCrOriginIndex + (strideCr * 24);
        cbMeanOf16x16Blocks[12] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferCb[cbBlockIndex]), strideCb);
        crMeanOf16x16Blocks[12] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferCr[crBlockIndex]), strideCr);

        // (3,1)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[13] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferCb[cbBlockIndex]), strideCb);
        crMeanOf16x16Blocks[13] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferCr[crBlockIndex]), strideCr);

        // (3,2)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[14] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferCb[cbBlockIndex]), strideCb);
        crMeanOf16x16Blocks[14] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferCr[crBlockIndex]), strideCr);

        // (3,3)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[15] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferCb[cbBlockIndex]), strideCb);
        crMeanOf16x16Blocks[15] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferCr[crBlockIndex]), strideCr);
    }

    // 32x32
    cbMeanOf32x32Blocks[0] = (cbMeanOf16x16Blocks[0] + cbMeanOf16x16Blocks[1] + cbMeanOf16x16Blocks[4] + cbMeanOf16x16Blocks[5]) >> 2;
    crMeanOf32x32Blocks[0] = (crMeanOf16x16Blocks[0] + crMeanOf16x16Blocks[1] + crMeanOf16x16Blocks[4] + crMeanOf16x16Blocks[5]) >> 2;

    cbMeanOf32x32Blocks[1] = (cbMeanOf16x16Blocks[2] + cbMeanOf16x16Blocks[3] + cbMeanOf16x16Blocks[6] + cbMeanOf16x16Blocks[7]) >> 2;
    crMeanOf32x32Blocks[1] = (crMeanOf16x16Blocks[2] + crMeanOf16x16Blocks[3] + crMeanOf16x16Blocks[6] + crMeanOf16x16Blocks[7]) >> 2;


    cbMeanOf32x32Blocks[2] = (cbMeanOf16x16Blocks[8] + cbMeanOf16x16Blocks[9] + cbMeanOf16x16Blocks[12] + cbMeanOf16x16Blocks[13]) >> 2;
    crMeanOf32x32Blocks[2] = (crMeanOf16x16Blocks[8] + crMeanOf16x16Blocks[9] + crMeanOf16x16Blocks[12] + crMeanOf16x16Blocks[13]) >> 2;

    cbMeanOf32x32Blocks[3] = (cbMeanOf16x16Blocks[10] + cbMeanOf16x16Blocks[11] + cbMeanOf16x16Blocks[14] + cbMeanOf16x16Blocks[15]) >> 2;
    crMeanOf32x32Blocks[3] = (crMeanOf16x16Blocks[10] + crMeanOf16x16Blocks[11] + crMeanOf16x16Blocks[14] + crMeanOf16x16Blocks[15]) >> 2;

    // 64x64
    cbMeanOf64x64Blocks = (cbMeanOf32x32Blocks[0] + cbMeanOf32x32Blocks[1] + cbMeanOf32x32Blocks[3] + cbMeanOf32x32Blocks[3]) >> 2;
    crMeanOf64x64Blocks = (crMeanOf32x32Blocks[0] + crMeanOf32x32Blocks[1] + crMeanOf32x32Blocks[3] + crMeanOf32x32Blocks[3]) >> 2;
    // 16x16 mean
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_0] = (uint8_t)(cbMeanOf16x16Blocks[0] >> MEAN_PRECISION);
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_1] = (uint8_t)(cbMeanOf16x16Blocks[1] >> MEAN_PRECISION);
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_2] = (uint8_t)(cbMeanOf16x16Blocks[2] >> MEAN_PRECISION);
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_3] = (uint8_t)(cbMeanOf16x16Blocks[3] >> MEAN_PRECISION);
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_4] = (uint8_t)(cbMeanOf16x16Blocks[4] >> MEAN_PRECISION);
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_5] = (uint8_t)(cbMeanOf16x16Blocks[5] >> MEAN_PRECISION);
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_6] = (uint8_t)(cbMeanOf16x16Blocks[6] >> MEAN_PRECISION);
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_7] = (uint8_t)(cbMeanOf16x16Blocks[7] >> MEAN_PRECISION);
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_8] = (uint8_t)(cbMeanOf16x16Blocks[8] >> MEAN_PRECISION);
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_9] = (uint8_t)(cbMeanOf16x16Blocks[9] >> MEAN_PRECISION);
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_10] = (uint8_t)(cbMeanOf16x16Blocks[10] >> MEAN_PRECISION);
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_11] = (uint8_t)(cbMeanOf16x16Blocks[11] >> MEAN_PRECISION);
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_12] = (uint8_t)(cbMeanOf16x16Blocks[12] >> MEAN_PRECISION);
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_13] = (uint8_t)(cbMeanOf16x16Blocks[13] >> MEAN_PRECISION);
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_14] = (uint8_t)(cbMeanOf16x16Blocks[14] >> MEAN_PRECISION);
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_15] = (uint8_t)(cbMeanOf16x16Blocks[15] >> MEAN_PRECISION);

    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_0] = (uint8_t)(crMeanOf16x16Blocks[0] >> MEAN_PRECISION);
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_1] = (uint8_t)(crMeanOf16x16Blocks[1] >> MEAN_PRECISION);
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_2] = (uint8_t)(crMeanOf16x16Blocks[2] >> MEAN_PRECISION);
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_3] = (uint8_t)(crMeanOf16x16Blocks[3] >> MEAN_PRECISION);
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_4] = (uint8_t)(crMeanOf16x16Blocks[4] >> MEAN_PRECISION);
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_5] = (uint8_t)(crMeanOf16x16Blocks[5] >> MEAN_PRECISION);
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_6] = (uint8_t)(crMeanOf16x16Blocks[6] >> MEAN_PRECISION);
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_7] = (uint8_t)(crMeanOf16x16Blocks[7] >> MEAN_PRECISION);
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_8] = (uint8_t)(crMeanOf16x16Blocks[8] >> MEAN_PRECISION);
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_9] = (uint8_t)(crMeanOf16x16Blocks[9] >> MEAN_PRECISION);
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_10] = (uint8_t)(crMeanOf16x16Blocks[10] >> MEAN_PRECISION);
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_11] = (uint8_t)(crMeanOf16x16Blocks[11] >> MEAN_PRECISION);
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_12] = (uint8_t)(crMeanOf16x16Blocks[12] >> MEAN_PRECISION);
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_13] = (uint8_t)(crMeanOf16x16Blocks[13] >> MEAN_PRECISION);
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_14] = (uint8_t)(crMeanOf16x16Blocks[14] >> MEAN_PRECISION);
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_15] = (uint8_t)(crMeanOf16x16Blocks[15] >> MEAN_PRECISION);

    // 32x32 mean
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_32x32_0] = (uint8_t)(cbMeanOf32x32Blocks[0] >> MEAN_PRECISION);
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_32x32_1] = (uint8_t)(cbMeanOf32x32Blocks[1] >> MEAN_PRECISION);
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_32x32_2] = (uint8_t)(cbMeanOf32x32Blocks[2] >> MEAN_PRECISION);
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_32x32_3] = (uint8_t)(cbMeanOf32x32Blocks[3] >> MEAN_PRECISION);

    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_32x32_0] = (uint8_t)(crMeanOf32x32Blocks[0] >> MEAN_PRECISION);
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_32x32_1] = (uint8_t)(crMeanOf32x32Blocks[1] >> MEAN_PRECISION);
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_32x32_2] = (uint8_t)(crMeanOf32x32Blocks[2] >> MEAN_PRECISION);
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_32x32_3] = (uint8_t)(crMeanOf32x32Blocks[3] >> MEAN_PRECISION);

    // 64x64 mean
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_64x64] = (uint8_t)(cbMeanOf64x64Blocks >> MEAN_PRECISION);
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_64x64] = (uint8_t)(crMeanOf64x64Blocks >> MEAN_PRECISION);

    return return_error;
}


/*******************************************
* ComputeBlockMeanComputeVariance
*   computes the variance and the block mean of all CUs inside the tree block
*******************************************/
EbErrorType ComputeBlockMeanComputeVariance(
    SequenceControlSet_t        *sequence_control_set_ptr,
    PictureParentControlSet_t   *picture_control_set_ptr,          // input parameter, Picture Control Set Ptr
    EbPictureBufferDesc_t       *inputPaddedPicturePtr,         // input parameter, Input Padded Picture
    uint32_t                       sb_index,                // input parameter, SB address
    uint32_t                       inputLumaOriginIndex,          // input parameter, SB index, used to point to source/reference samples
    EbAsm                         asm_type)
{

    EbErrorType return_error = EB_ErrorNone;

    uint32_t blockIndex;

    uint64_t meanOf8x8Blocks[64];
    uint64_t meanOf8x8SquaredValuesBlocks[64];

    uint64_t meanOf16x16Blocks[16];
    uint64_t meanOf16x16SquaredValuesBlocks[16];

    uint64_t meanOf32x32Blocks[4];
    uint64_t meanOf32x32SquaredValuesBlocks[4];

    uint64_t meanOf64x64Blocks;
    uint64_t meanOf64x64SquaredValuesBlocks;

    // (0,0)
    blockIndex = inputLumaOriginIndex;
    if (sequence_control_set_ptr->block_mean_calc_prec == BLOCK_MEAN_PREC_FULL) {
        meanOf8x8Blocks[0] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[0] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (0,1)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[1] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[1] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (0,2)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[2] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[2] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (0,3)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[3] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[3] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (0,4)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[4] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[4] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (0,5)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[5] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[5] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (0,6)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[6] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[6] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (0,7)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[7] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[7] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (1,0)
        blockIndex = inputLumaOriginIndex + (inputPaddedPicturePtr->strideY << 3);
        meanOf8x8Blocks[8] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[8] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (1,1)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[9] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[9] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (1,2)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[10] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[10] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (1,3)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[11] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[11] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (1,4)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[12] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[12] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (1,5)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[13] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[13] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (1,6)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[14] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[14] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (1,7)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[15] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[15] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (2,0)
        blockIndex = inputLumaOriginIndex + (inputPaddedPicturePtr->strideY << 4);
        meanOf8x8Blocks[16] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[16] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (2,1)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[17] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[17] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (2,2)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[18] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[18] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (2,3)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[19] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[19] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        /// (2,4)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[20] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[20] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (2,5)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[21] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[21] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (2,6)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[22] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[22] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (2,7)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[23] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[23] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (3,0)
        blockIndex = inputLumaOriginIndex + (inputPaddedPicturePtr->strideY << 3) + (inputPaddedPicturePtr->strideY << 4);
        meanOf8x8Blocks[24] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[24] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (3,1)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[25] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[25] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (3,2)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[26] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[26] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (3,3)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[27] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[27] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (3,4)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[28] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[28] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (3,5)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[29] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[29] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (3,6)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[30] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[30] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (3,7)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[31] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[31] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (4,0)
        blockIndex = inputLumaOriginIndex + (inputPaddedPicturePtr->strideY << 5);
        meanOf8x8Blocks[32] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[32] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (4,1)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[33] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[33] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (4,2)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[34] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[34] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (4,3)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[35] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[35] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (4,4)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[36] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[36] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (4,5)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[37] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[37] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (4,6)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[38] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[38] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (4,7)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[39] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[39] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (5,0)
        blockIndex = inputLumaOriginIndex + (inputPaddedPicturePtr->strideY << 3) + (inputPaddedPicturePtr->strideY << 5);
        meanOf8x8Blocks[40] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[40] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (5,1)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[41] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[41] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (5,2)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[42] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[42] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (5,3)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[43] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[43] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (5,4)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[44] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[44] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (5,5)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[45] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[45] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (5,6)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[46] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[46] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (5,7)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[47] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[47] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (6,0)
        blockIndex = inputLumaOriginIndex + (inputPaddedPicturePtr->strideY << 4) + (inputPaddedPicturePtr->strideY << 5);
        meanOf8x8Blocks[48] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[48] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (6,1)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[49] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[49] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (6,2)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[50] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[50] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (6,3)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[51] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[51] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (6,4)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[52] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[52] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (6,5)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[53] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[53] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (6,6)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[54] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[54] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (6,7)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[55] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[55] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (7,0)
        blockIndex = inputLumaOriginIndex + (inputPaddedPicturePtr->strideY << 3) + (inputPaddedPicturePtr->strideY << 4) + (inputPaddedPicturePtr->strideY << 5);
        meanOf8x8Blocks[56] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[56] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (7,1)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[57] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[57] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (7,2)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[58] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[58] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (7,3)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[59] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[59] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (7,4)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[60] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[60] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (7,5)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[61] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[61] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (7,6)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[62] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[62] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);

        // (7,7)
        blockIndex = blockIndex + 8;
        meanOf8x8Blocks[63] = ComputeMeanFunc[0][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
        meanOf8x8SquaredValuesBlocks[63] = ComputeMeanFunc[1][asm_type](&(inputPaddedPicturePtr->bufferY[blockIndex]), inputPaddedPicturePtr->strideY, 8, 8);
    }
    else {
        const uint16_t strideY = inputPaddedPicturePtr->strideY;

        if (asm_type == ASM_AVX2) {

            ComputeIntermVarFour8x8_AVX2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY, &meanOf8x8Blocks[0], &meanOf8x8SquaredValuesBlocks[0]);

            // (0,1)
            blockIndex = blockIndex + 32;


            ComputeIntermVarFour8x8_AVX2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY, &meanOf8x8Blocks[4], &meanOf8x8SquaredValuesBlocks[4]);

            // (0,5)
            blockIndex = blockIndex + 24;

            // (1,0)
            blockIndex = inputLumaOriginIndex + (strideY << 3);

            ComputeIntermVarFour8x8_AVX2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY, &meanOf8x8Blocks[8], &meanOf8x8SquaredValuesBlocks[8]);

            // (1,1)
            blockIndex = blockIndex + 32;

            ComputeIntermVarFour8x8_AVX2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY, &meanOf8x8Blocks[12], &meanOf8x8SquaredValuesBlocks[12]);

            // (1,5)
            blockIndex = blockIndex + 24;

            // (2,0)
            blockIndex = inputLumaOriginIndex + (strideY << 4);

            ComputeIntermVarFour8x8_AVX2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY, &meanOf8x8Blocks[16], &meanOf8x8SquaredValuesBlocks[16]);

            // (2,1)
            blockIndex = blockIndex + 32;

            ComputeIntermVarFour8x8_AVX2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY, &meanOf8x8Blocks[20], &meanOf8x8SquaredValuesBlocks[20]);

            // (2,5)
            blockIndex = blockIndex + 24;

            // (3,0)
            blockIndex = inputLumaOriginIndex + (strideY << 3) + (strideY << 4);


            ComputeIntermVarFour8x8_AVX2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY, &meanOf8x8Blocks[24], &meanOf8x8SquaredValuesBlocks[24]);

            // (3,1)
            blockIndex = blockIndex + 32;

            ComputeIntermVarFour8x8_AVX2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY, &meanOf8x8Blocks[28], &meanOf8x8SquaredValuesBlocks[28]);

            // (3,5)
            blockIndex = blockIndex + 24;

            // (4,0)
            blockIndex = inputLumaOriginIndex + (strideY << 5);

            ComputeIntermVarFour8x8_AVX2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY, &meanOf8x8Blocks[32], &meanOf8x8SquaredValuesBlocks[32]);

            // (4,1)
            blockIndex = blockIndex + 32;

            ComputeIntermVarFour8x8_AVX2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY, &meanOf8x8Blocks[36], &meanOf8x8SquaredValuesBlocks[36]);

            // (4,5)
            blockIndex = blockIndex + 24;

            // (5,0)
            blockIndex = inputLumaOriginIndex + (strideY << 3) + (strideY << 5);

            ComputeIntermVarFour8x8_AVX2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY, &meanOf8x8Blocks[40], &meanOf8x8SquaredValuesBlocks[40]);

            // (5,1)
            blockIndex = blockIndex + 32;

            ComputeIntermVarFour8x8_AVX2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY, &meanOf8x8Blocks[44], &meanOf8x8SquaredValuesBlocks[44]);

            // (5,5)
            blockIndex = blockIndex + 24;

            // (6,0)
            blockIndex = inputLumaOriginIndex + (strideY << 4) + (strideY << 5);

            ComputeIntermVarFour8x8_AVX2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY, &meanOf8x8Blocks[48], &meanOf8x8SquaredValuesBlocks[48]);

            // (6,1)
            blockIndex = blockIndex + 32;

            ComputeIntermVarFour8x8_AVX2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY, &meanOf8x8Blocks[52], &meanOf8x8SquaredValuesBlocks[52]);

            // (6,5)
            blockIndex = blockIndex + 24;


            // (7,0)
            blockIndex = inputLumaOriginIndex + (strideY << 3) + (strideY << 4) + (strideY << 5);

            ComputeIntermVarFour8x8_AVX2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY, &meanOf8x8Blocks[56], &meanOf8x8SquaredValuesBlocks[56]);


            // (7,1)
            blockIndex = blockIndex + 32;

            ComputeIntermVarFour8x8_AVX2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY, &meanOf8x8Blocks[60], &meanOf8x8SquaredValuesBlocks[60]);


        }
        else {
            meanOf8x8Blocks[0] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[0] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (0,1)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[1] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[1] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (0,2)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[2] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[2] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (0,3)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[3] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[3] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (0,4)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[4] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[4] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (0,5)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[5] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[5] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (0,6)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[6] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[6] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (0,7)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[7] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[7] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (1,0)
            blockIndex = inputLumaOriginIndex + (strideY << 3);
            meanOf8x8Blocks[8] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[8] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (1,1)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[9] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[9] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (1,2)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[10] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[10] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (1,3)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[11] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[11] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (1,4)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[12] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[12] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (1,5)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[13] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[13] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (1,6)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[14] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[14] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (1,7)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[15] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[15] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (2,0)
            blockIndex = inputLumaOriginIndex + (strideY << 4);
            meanOf8x8Blocks[16] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[16] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (2,1)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[17] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[17] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (2,2)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[18] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[18] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (2,3)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[19] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[19] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            /// (2,4)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[20] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[20] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (2,5)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[21] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[21] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (2,6)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[22] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[22] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (2,7)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[23] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[23] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (3,0)
            blockIndex = inputLumaOriginIndex + (strideY << 3) + (strideY << 4);
            meanOf8x8Blocks[24] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[24] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (3,1)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[25] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[25] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (3,2)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[26] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[26] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (3,3)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[27] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[27] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (3,4)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[28] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[28] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (3,5)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[29] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[29] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (3,6)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[30] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[30] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (3,7)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[31] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[31] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (4,0)
            blockIndex = inputLumaOriginIndex + (strideY << 5);
            meanOf8x8Blocks[32] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[32] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (4,1)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[33] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[33] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (4,2)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[34] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[34] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (4,3)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[35] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[35] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (4,4)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[36] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[36] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (4,5)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[37] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[37] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (4,6)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[38] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[38] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (4,7)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[39] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[39] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (5,0)
            blockIndex = inputLumaOriginIndex + (strideY << 3) + (strideY << 5);
            meanOf8x8Blocks[40] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[40] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (5,1)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[41] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[41] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (5,2)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[42] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[42] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (5,3)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[43] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[43] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (5,4)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[44] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[44] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (5,5)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[45] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[45] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (5,6)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[46] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[46] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (5,7)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[47] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[47] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (6,0)
            blockIndex = inputLumaOriginIndex + (strideY << 4) + (strideY << 5);
            meanOf8x8Blocks[48] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[48] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (6,1)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[49] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[49] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (6,2)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[50] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[50] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (6,3)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[51] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[51] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (6,4)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[52] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[52] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (6,5)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[53] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[53] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (6,6)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[54] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[54] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (6,7)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[55] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[55] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (7,0)
            blockIndex = inputLumaOriginIndex + (strideY << 3) + (strideY << 4) + (strideY << 5);
            meanOf8x8Blocks[56] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[56] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (7,1)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[57] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[57] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (7,2)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[58] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[58] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (7,3)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[59] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[59] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (7,4)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[60] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[60] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (7,5)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[61] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[61] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (7,6)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[62] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[62] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);

            // (7,7)
            blockIndex = blockIndex + 8;
            meanOf8x8Blocks[63] = ComputeSubMean8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
            meanOf8x8SquaredValuesBlocks[63] = ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(&(inputPaddedPicturePtr->bufferY[blockIndex]), strideY);
        }
    }

    // 16x16
    meanOf16x16Blocks[0] = (meanOf8x8Blocks[0] + meanOf8x8Blocks[1] + meanOf8x8Blocks[8] + meanOf8x8Blocks[9]) >> 2;
    meanOf16x16Blocks[1] = (meanOf8x8Blocks[2] + meanOf8x8Blocks[3] + meanOf8x8Blocks[10] + meanOf8x8Blocks[11]) >> 2;
    meanOf16x16Blocks[2] = (meanOf8x8Blocks[4] + meanOf8x8Blocks[5] + meanOf8x8Blocks[12] + meanOf8x8Blocks[13]) >> 2;
    meanOf16x16Blocks[3] = (meanOf8x8Blocks[6] + meanOf8x8Blocks[7] + meanOf8x8Blocks[14] + meanOf8x8Blocks[15]) >> 2;

    meanOf16x16Blocks[4] = (meanOf8x8Blocks[16] + meanOf8x8Blocks[17] + meanOf8x8Blocks[24] + meanOf8x8Blocks[25]) >> 2;
    meanOf16x16Blocks[5] = (meanOf8x8Blocks[18] + meanOf8x8Blocks[19] + meanOf8x8Blocks[26] + meanOf8x8Blocks[27]) >> 2;
    meanOf16x16Blocks[6] = (meanOf8x8Blocks[20] + meanOf8x8Blocks[21] + meanOf8x8Blocks[28] + meanOf8x8Blocks[29]) >> 2;
    meanOf16x16Blocks[7] = (meanOf8x8Blocks[22] + meanOf8x8Blocks[23] + meanOf8x8Blocks[30] + meanOf8x8Blocks[31]) >> 2;

    meanOf16x16Blocks[8] = (meanOf8x8Blocks[32] + meanOf8x8Blocks[33] + meanOf8x8Blocks[40] + meanOf8x8Blocks[41]) >> 2;
    meanOf16x16Blocks[9] = (meanOf8x8Blocks[34] + meanOf8x8Blocks[35] + meanOf8x8Blocks[42] + meanOf8x8Blocks[43]) >> 2;
    meanOf16x16Blocks[10] = (meanOf8x8Blocks[36] + meanOf8x8Blocks[37] + meanOf8x8Blocks[44] + meanOf8x8Blocks[45]) >> 2;
    meanOf16x16Blocks[11] = (meanOf8x8Blocks[38] + meanOf8x8Blocks[39] + meanOf8x8Blocks[46] + meanOf8x8Blocks[47]) >> 2;

    meanOf16x16Blocks[12] = (meanOf8x8Blocks[48] + meanOf8x8Blocks[49] + meanOf8x8Blocks[56] + meanOf8x8Blocks[57]) >> 2;
    meanOf16x16Blocks[13] = (meanOf8x8Blocks[50] + meanOf8x8Blocks[51] + meanOf8x8Blocks[58] + meanOf8x8Blocks[59]) >> 2;
    meanOf16x16Blocks[14] = (meanOf8x8Blocks[52] + meanOf8x8Blocks[53] + meanOf8x8Blocks[60] + meanOf8x8Blocks[61]) >> 2;
    meanOf16x16Blocks[15] = (meanOf8x8Blocks[54] + meanOf8x8Blocks[55] + meanOf8x8Blocks[62] + meanOf8x8Blocks[63]) >> 2;

    meanOf16x16SquaredValuesBlocks[0] = (meanOf8x8SquaredValuesBlocks[0] + meanOf8x8SquaredValuesBlocks[1] + meanOf8x8SquaredValuesBlocks[8] + meanOf8x8SquaredValuesBlocks[9]) >> 2;
    meanOf16x16SquaredValuesBlocks[1] = (meanOf8x8SquaredValuesBlocks[2] + meanOf8x8SquaredValuesBlocks[3] + meanOf8x8SquaredValuesBlocks[10] + meanOf8x8SquaredValuesBlocks[11]) >> 2;
    meanOf16x16SquaredValuesBlocks[2] = (meanOf8x8SquaredValuesBlocks[4] + meanOf8x8SquaredValuesBlocks[5] + meanOf8x8SquaredValuesBlocks[12] + meanOf8x8SquaredValuesBlocks[13]) >> 2;
    meanOf16x16SquaredValuesBlocks[3] = (meanOf8x8SquaredValuesBlocks[6] + meanOf8x8SquaredValuesBlocks[7] + meanOf8x8SquaredValuesBlocks[14] + meanOf8x8SquaredValuesBlocks[15]) >> 2;

    meanOf16x16SquaredValuesBlocks[4] = (meanOf8x8SquaredValuesBlocks[16] + meanOf8x8SquaredValuesBlocks[17] + meanOf8x8SquaredValuesBlocks[24] + meanOf8x8SquaredValuesBlocks[25]) >> 2;
    meanOf16x16SquaredValuesBlocks[5] = (meanOf8x8SquaredValuesBlocks[18] + meanOf8x8SquaredValuesBlocks[19] + meanOf8x8SquaredValuesBlocks[26] + meanOf8x8SquaredValuesBlocks[27]) >> 2;
    meanOf16x16SquaredValuesBlocks[6] = (meanOf8x8SquaredValuesBlocks[20] + meanOf8x8SquaredValuesBlocks[21] + meanOf8x8SquaredValuesBlocks[28] + meanOf8x8SquaredValuesBlocks[29]) >> 2;
    meanOf16x16SquaredValuesBlocks[7] = (meanOf8x8SquaredValuesBlocks[22] + meanOf8x8SquaredValuesBlocks[23] + meanOf8x8SquaredValuesBlocks[30] + meanOf8x8SquaredValuesBlocks[31]) >> 2;

    meanOf16x16SquaredValuesBlocks[8] = (meanOf8x8SquaredValuesBlocks[32] + meanOf8x8SquaredValuesBlocks[33] + meanOf8x8SquaredValuesBlocks[40] + meanOf8x8SquaredValuesBlocks[41]) >> 2;
    meanOf16x16SquaredValuesBlocks[9] = (meanOf8x8SquaredValuesBlocks[34] + meanOf8x8SquaredValuesBlocks[35] + meanOf8x8SquaredValuesBlocks[42] + meanOf8x8SquaredValuesBlocks[43]) >> 2;
    meanOf16x16SquaredValuesBlocks[10] = (meanOf8x8SquaredValuesBlocks[36] + meanOf8x8SquaredValuesBlocks[37] + meanOf8x8SquaredValuesBlocks[44] + meanOf8x8SquaredValuesBlocks[45]) >> 2;
    meanOf16x16SquaredValuesBlocks[11] = (meanOf8x8SquaredValuesBlocks[38] + meanOf8x8SquaredValuesBlocks[39] + meanOf8x8SquaredValuesBlocks[46] + meanOf8x8SquaredValuesBlocks[47]) >> 2;

    meanOf16x16SquaredValuesBlocks[12] = (meanOf8x8SquaredValuesBlocks[48] + meanOf8x8SquaredValuesBlocks[49] + meanOf8x8SquaredValuesBlocks[56] + meanOf8x8SquaredValuesBlocks[57]) >> 2;
    meanOf16x16SquaredValuesBlocks[13] = (meanOf8x8SquaredValuesBlocks[50] + meanOf8x8SquaredValuesBlocks[51] + meanOf8x8SquaredValuesBlocks[58] + meanOf8x8SquaredValuesBlocks[59]) >> 2;
    meanOf16x16SquaredValuesBlocks[14] = (meanOf8x8SquaredValuesBlocks[52] + meanOf8x8SquaredValuesBlocks[53] + meanOf8x8SquaredValuesBlocks[60] + meanOf8x8SquaredValuesBlocks[61]) >> 2;
    meanOf16x16SquaredValuesBlocks[15] = (meanOf8x8SquaredValuesBlocks[54] + meanOf8x8SquaredValuesBlocks[55] + meanOf8x8SquaredValuesBlocks[62] + meanOf8x8SquaredValuesBlocks[63]) >> 2;

    // 32x32
    meanOf32x32Blocks[0] = (meanOf16x16Blocks[0] + meanOf16x16Blocks[1] + meanOf16x16Blocks[4] + meanOf16x16Blocks[5]) >> 2;
    meanOf32x32Blocks[1] = (meanOf16x16Blocks[2] + meanOf16x16Blocks[3] + meanOf16x16Blocks[6] + meanOf16x16Blocks[7]) >> 2;
    meanOf32x32Blocks[2] = (meanOf16x16Blocks[8] + meanOf16x16Blocks[9] + meanOf16x16Blocks[12] + meanOf16x16Blocks[13]) >> 2;
    meanOf32x32Blocks[3] = (meanOf16x16Blocks[10] + meanOf16x16Blocks[11] + meanOf16x16Blocks[14] + meanOf16x16Blocks[15]) >> 2;

    meanOf32x32SquaredValuesBlocks[0] = (meanOf16x16SquaredValuesBlocks[0] + meanOf16x16SquaredValuesBlocks[1] + meanOf16x16SquaredValuesBlocks[4] + meanOf16x16SquaredValuesBlocks[5]) >> 2;
    meanOf32x32SquaredValuesBlocks[1] = (meanOf16x16SquaredValuesBlocks[2] + meanOf16x16SquaredValuesBlocks[3] + meanOf16x16SquaredValuesBlocks[6] + meanOf16x16SquaredValuesBlocks[7]) >> 2;
    meanOf32x32SquaredValuesBlocks[2] = (meanOf16x16SquaredValuesBlocks[8] + meanOf16x16SquaredValuesBlocks[9] + meanOf16x16SquaredValuesBlocks[12] + meanOf16x16SquaredValuesBlocks[13]) >> 2;
    meanOf32x32SquaredValuesBlocks[3] = (meanOf16x16SquaredValuesBlocks[10] + meanOf16x16SquaredValuesBlocks[11] + meanOf16x16SquaredValuesBlocks[14] + meanOf16x16SquaredValuesBlocks[15]) >> 2;

    // 64x64
    meanOf64x64Blocks = (meanOf32x32Blocks[0] + meanOf32x32Blocks[1] + meanOf32x32Blocks[2] + meanOf32x32Blocks[3]) >> 2;
    meanOf64x64SquaredValuesBlocks = (meanOf32x32SquaredValuesBlocks[0] + meanOf32x32SquaredValuesBlocks[1] + meanOf32x32SquaredValuesBlocks[2] + meanOf32x32SquaredValuesBlocks[3]) >> 2;

    // 8x8 means
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_0] = (uint8_t)(meanOf8x8Blocks[0] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_1] = (uint8_t)(meanOf8x8Blocks[1] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_2] = (uint8_t)(meanOf8x8Blocks[2] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_3] = (uint8_t)(meanOf8x8Blocks[3] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_4] = (uint8_t)(meanOf8x8Blocks[4] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_5] = (uint8_t)(meanOf8x8Blocks[5] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_6] = (uint8_t)(meanOf8x8Blocks[6] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_7] = (uint8_t)(meanOf8x8Blocks[7] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_8] = (uint8_t)(meanOf8x8Blocks[8] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_9] = (uint8_t)(meanOf8x8Blocks[9] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_10] = (uint8_t)(meanOf8x8Blocks[10] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_11] = (uint8_t)(meanOf8x8Blocks[11] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_12] = (uint8_t)(meanOf8x8Blocks[12] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_13] = (uint8_t)(meanOf8x8Blocks[13] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_14] = (uint8_t)(meanOf8x8Blocks[14] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_15] = (uint8_t)(meanOf8x8Blocks[15] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_16] = (uint8_t)(meanOf8x8Blocks[16] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_17] = (uint8_t)(meanOf8x8Blocks[17] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_18] = (uint8_t)(meanOf8x8Blocks[18] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_19] = (uint8_t)(meanOf8x8Blocks[19] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_20] = (uint8_t)(meanOf8x8Blocks[20] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_21] = (uint8_t)(meanOf8x8Blocks[21] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_22] = (uint8_t)(meanOf8x8Blocks[22] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_23] = (uint8_t)(meanOf8x8Blocks[23] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_24] = (uint8_t)(meanOf8x8Blocks[24] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_25] = (uint8_t)(meanOf8x8Blocks[25] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_26] = (uint8_t)(meanOf8x8Blocks[26] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_27] = (uint8_t)(meanOf8x8Blocks[27] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_28] = (uint8_t)(meanOf8x8Blocks[28] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_29] = (uint8_t)(meanOf8x8Blocks[29] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_30] = (uint8_t)(meanOf8x8Blocks[30] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_31] = (uint8_t)(meanOf8x8Blocks[31] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_32] = (uint8_t)(meanOf8x8Blocks[32] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_33] = (uint8_t)(meanOf8x8Blocks[33] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_34] = (uint8_t)(meanOf8x8Blocks[34] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_35] = (uint8_t)(meanOf8x8Blocks[35] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_36] = (uint8_t)(meanOf8x8Blocks[36] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_37] = (uint8_t)(meanOf8x8Blocks[37] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_38] = (uint8_t)(meanOf8x8Blocks[38] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_39] = (uint8_t)(meanOf8x8Blocks[39] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_40] = (uint8_t)(meanOf8x8Blocks[40] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_41] = (uint8_t)(meanOf8x8Blocks[41] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_42] = (uint8_t)(meanOf8x8Blocks[42] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_43] = (uint8_t)(meanOf8x8Blocks[43] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_44] = (uint8_t)(meanOf8x8Blocks[44] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_45] = (uint8_t)(meanOf8x8Blocks[45] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_46] = (uint8_t)(meanOf8x8Blocks[46] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_47] = (uint8_t)(meanOf8x8Blocks[47] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_48] = (uint8_t)(meanOf8x8Blocks[48] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_49] = (uint8_t)(meanOf8x8Blocks[49] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_50] = (uint8_t)(meanOf8x8Blocks[50] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_51] = (uint8_t)(meanOf8x8Blocks[51] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_52] = (uint8_t)(meanOf8x8Blocks[52] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_53] = (uint8_t)(meanOf8x8Blocks[53] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_54] = (uint8_t)(meanOf8x8Blocks[54] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_55] = (uint8_t)(meanOf8x8Blocks[55] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_56] = (uint8_t)(meanOf8x8Blocks[56] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_57] = (uint8_t)(meanOf8x8Blocks[57] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_58] = (uint8_t)(meanOf8x8Blocks[58] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_59] = (uint8_t)(meanOf8x8Blocks[59] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_60] = (uint8_t)(meanOf8x8Blocks[60] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_61] = (uint8_t)(meanOf8x8Blocks[61] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_62] = (uint8_t)(meanOf8x8Blocks[62] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_8x8_63] = (uint8_t)(meanOf8x8Blocks[63] >> MEAN_PRECISION);

    // 16x16 mean
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_16x16_0] = (uint8_t)(meanOf16x16Blocks[0] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_16x16_1] = (uint8_t)(meanOf16x16Blocks[1] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_16x16_2] = (uint8_t)(meanOf16x16Blocks[2] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_16x16_3] = (uint8_t)(meanOf16x16Blocks[3] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_16x16_4] = (uint8_t)(meanOf16x16Blocks[4] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_16x16_5] = (uint8_t)(meanOf16x16Blocks[5] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_16x16_6] = (uint8_t)(meanOf16x16Blocks[6] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_16x16_7] = (uint8_t)(meanOf16x16Blocks[7] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_16x16_8] = (uint8_t)(meanOf16x16Blocks[8] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_16x16_9] = (uint8_t)(meanOf16x16Blocks[9] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_16x16_10] = (uint8_t)(meanOf16x16Blocks[10] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_16x16_11] = (uint8_t)(meanOf16x16Blocks[11] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_16x16_12] = (uint8_t)(meanOf16x16Blocks[12] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_16x16_13] = (uint8_t)(meanOf16x16Blocks[13] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_16x16_14] = (uint8_t)(meanOf16x16Blocks[14] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_16x16_15] = (uint8_t)(meanOf16x16Blocks[15] >> MEAN_PRECISION);

    // 32x32 mean
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_32x32_0] = (uint8_t)(meanOf32x32Blocks[0] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_32x32_1] = (uint8_t)(meanOf32x32Blocks[1] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_32x32_2] = (uint8_t)(meanOf32x32Blocks[2] >> MEAN_PRECISION);
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_32x32_3] = (uint8_t)(meanOf32x32Blocks[3] >> MEAN_PRECISION);

    // 64x64 mean
    picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_64x64] = (uint8_t)(meanOf64x64Blocks >> MEAN_PRECISION);

    // 8x8 variances
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_0] = (uint16_t)((meanOf8x8SquaredValuesBlocks[0] - (meanOf8x8Blocks[0] * meanOf8x8Blocks[0])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_1] = (uint16_t)((meanOf8x8SquaredValuesBlocks[1] - (meanOf8x8Blocks[1] * meanOf8x8Blocks[1])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_2] = (uint16_t)((meanOf8x8SquaredValuesBlocks[2] - (meanOf8x8Blocks[2] * meanOf8x8Blocks[2])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_3] = (uint16_t)((meanOf8x8SquaredValuesBlocks[3] - (meanOf8x8Blocks[3] * meanOf8x8Blocks[3])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_4] = (uint16_t)((meanOf8x8SquaredValuesBlocks[4] - (meanOf8x8Blocks[4] * meanOf8x8Blocks[4])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_5] = (uint16_t)((meanOf8x8SquaredValuesBlocks[5] - (meanOf8x8Blocks[5] * meanOf8x8Blocks[5])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_6] = (uint16_t)((meanOf8x8SquaredValuesBlocks[6] - (meanOf8x8Blocks[6] * meanOf8x8Blocks[6])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_7] = (uint16_t)((meanOf8x8SquaredValuesBlocks[7] - (meanOf8x8Blocks[7] * meanOf8x8Blocks[7])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_8] = (uint16_t)((meanOf8x8SquaredValuesBlocks[8] - (meanOf8x8Blocks[8] * meanOf8x8Blocks[8])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_9] = (uint16_t)((meanOf8x8SquaredValuesBlocks[9] - (meanOf8x8Blocks[9] * meanOf8x8Blocks[9])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_10] = (uint16_t)((meanOf8x8SquaredValuesBlocks[10] - (meanOf8x8Blocks[10] * meanOf8x8Blocks[10])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_11] = (uint16_t)((meanOf8x8SquaredValuesBlocks[11] - (meanOf8x8Blocks[11] * meanOf8x8Blocks[11])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_12] = (uint16_t)((meanOf8x8SquaredValuesBlocks[12] - (meanOf8x8Blocks[12] * meanOf8x8Blocks[12])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_13] = (uint16_t)((meanOf8x8SquaredValuesBlocks[13] - (meanOf8x8Blocks[13] * meanOf8x8Blocks[13])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_14] = (uint16_t)((meanOf8x8SquaredValuesBlocks[14] - (meanOf8x8Blocks[14] * meanOf8x8Blocks[14])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_15] = (uint16_t)((meanOf8x8SquaredValuesBlocks[15] - (meanOf8x8Blocks[15] * meanOf8x8Blocks[15])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_16] = (uint16_t)((meanOf8x8SquaredValuesBlocks[16] - (meanOf8x8Blocks[16] * meanOf8x8Blocks[16])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_17] = (uint16_t)((meanOf8x8SquaredValuesBlocks[17] - (meanOf8x8Blocks[17] * meanOf8x8Blocks[17])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_18] = (uint16_t)((meanOf8x8SquaredValuesBlocks[18] - (meanOf8x8Blocks[18] * meanOf8x8Blocks[18])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_19] = (uint16_t)((meanOf8x8SquaredValuesBlocks[19] - (meanOf8x8Blocks[19] * meanOf8x8Blocks[19])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_20] = (uint16_t)((meanOf8x8SquaredValuesBlocks[20] - (meanOf8x8Blocks[20] * meanOf8x8Blocks[20])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_21] = (uint16_t)((meanOf8x8SquaredValuesBlocks[21] - (meanOf8x8Blocks[21] * meanOf8x8Blocks[21])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_22] = (uint16_t)((meanOf8x8SquaredValuesBlocks[22] - (meanOf8x8Blocks[22] * meanOf8x8Blocks[22])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_23] = (uint16_t)((meanOf8x8SquaredValuesBlocks[23] - (meanOf8x8Blocks[23] * meanOf8x8Blocks[23])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_24] = (uint16_t)((meanOf8x8SquaredValuesBlocks[24] - (meanOf8x8Blocks[24] * meanOf8x8Blocks[24])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_25] = (uint16_t)((meanOf8x8SquaredValuesBlocks[25] - (meanOf8x8Blocks[25] * meanOf8x8Blocks[25])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_26] = (uint16_t)((meanOf8x8SquaredValuesBlocks[26] - (meanOf8x8Blocks[26] * meanOf8x8Blocks[26])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_27] = (uint16_t)((meanOf8x8SquaredValuesBlocks[27] - (meanOf8x8Blocks[27] * meanOf8x8Blocks[27])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_28] = (uint16_t)((meanOf8x8SquaredValuesBlocks[28] - (meanOf8x8Blocks[28] * meanOf8x8Blocks[28])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_29] = (uint16_t)((meanOf8x8SquaredValuesBlocks[29] - (meanOf8x8Blocks[29] * meanOf8x8Blocks[29])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_30] = (uint16_t)((meanOf8x8SquaredValuesBlocks[30] - (meanOf8x8Blocks[30] * meanOf8x8Blocks[30])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_31] = (uint16_t)((meanOf8x8SquaredValuesBlocks[31] - (meanOf8x8Blocks[31] * meanOf8x8Blocks[31])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_32] = (uint16_t)((meanOf8x8SquaredValuesBlocks[32] - (meanOf8x8Blocks[32] * meanOf8x8Blocks[32])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_33] = (uint16_t)((meanOf8x8SquaredValuesBlocks[33] - (meanOf8x8Blocks[33] * meanOf8x8Blocks[33])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_34] = (uint16_t)((meanOf8x8SquaredValuesBlocks[34] - (meanOf8x8Blocks[34] * meanOf8x8Blocks[34])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_35] = (uint16_t)((meanOf8x8SquaredValuesBlocks[35] - (meanOf8x8Blocks[35] * meanOf8x8Blocks[35])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_36] = (uint16_t)((meanOf8x8SquaredValuesBlocks[36] - (meanOf8x8Blocks[36] * meanOf8x8Blocks[36])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_37] = (uint16_t)((meanOf8x8SquaredValuesBlocks[37] - (meanOf8x8Blocks[37] * meanOf8x8Blocks[37])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_38] = (uint16_t)((meanOf8x8SquaredValuesBlocks[38] - (meanOf8x8Blocks[38] * meanOf8x8Blocks[38])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_39] = (uint16_t)((meanOf8x8SquaredValuesBlocks[39] - (meanOf8x8Blocks[39] * meanOf8x8Blocks[39])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_40] = (uint16_t)((meanOf8x8SquaredValuesBlocks[40] - (meanOf8x8Blocks[40] * meanOf8x8Blocks[40])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_41] = (uint16_t)((meanOf8x8SquaredValuesBlocks[41] - (meanOf8x8Blocks[41] * meanOf8x8Blocks[41])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_42] = (uint16_t)((meanOf8x8SquaredValuesBlocks[42] - (meanOf8x8Blocks[42] * meanOf8x8Blocks[42])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_43] = (uint16_t)((meanOf8x8SquaredValuesBlocks[43] - (meanOf8x8Blocks[43] * meanOf8x8Blocks[43])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_44] = (uint16_t)((meanOf8x8SquaredValuesBlocks[44] - (meanOf8x8Blocks[44] * meanOf8x8Blocks[44])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_45] = (uint16_t)((meanOf8x8SquaredValuesBlocks[45] - (meanOf8x8Blocks[45] * meanOf8x8Blocks[45])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_46] = (uint16_t)((meanOf8x8SquaredValuesBlocks[46] - (meanOf8x8Blocks[46] * meanOf8x8Blocks[46])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_47] = (uint16_t)((meanOf8x8SquaredValuesBlocks[47] - (meanOf8x8Blocks[47] * meanOf8x8Blocks[47])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_48] = (uint16_t)((meanOf8x8SquaredValuesBlocks[48] - (meanOf8x8Blocks[48] * meanOf8x8Blocks[48])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_49] = (uint16_t)((meanOf8x8SquaredValuesBlocks[49] - (meanOf8x8Blocks[49] * meanOf8x8Blocks[49])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_50] = (uint16_t)((meanOf8x8SquaredValuesBlocks[50] - (meanOf8x8Blocks[50] * meanOf8x8Blocks[50])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_51] = (uint16_t)((meanOf8x8SquaredValuesBlocks[51] - (meanOf8x8Blocks[51] * meanOf8x8Blocks[51])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_52] = (uint16_t)((meanOf8x8SquaredValuesBlocks[52] - (meanOf8x8Blocks[52] * meanOf8x8Blocks[52])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_53] = (uint16_t)((meanOf8x8SquaredValuesBlocks[53] - (meanOf8x8Blocks[53] * meanOf8x8Blocks[53])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_54] = (uint16_t)((meanOf8x8SquaredValuesBlocks[54] - (meanOf8x8Blocks[54] * meanOf8x8Blocks[54])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_55] = (uint16_t)((meanOf8x8SquaredValuesBlocks[55] - (meanOf8x8Blocks[55] * meanOf8x8Blocks[55])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_56] = (uint16_t)((meanOf8x8SquaredValuesBlocks[56] - (meanOf8x8Blocks[56] * meanOf8x8Blocks[56])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_57] = (uint16_t)((meanOf8x8SquaredValuesBlocks[57] - (meanOf8x8Blocks[57] * meanOf8x8Blocks[57])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_58] = (uint16_t)((meanOf8x8SquaredValuesBlocks[58] - (meanOf8x8Blocks[58] * meanOf8x8Blocks[58])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_59] = (uint16_t)((meanOf8x8SquaredValuesBlocks[59] - (meanOf8x8Blocks[59] * meanOf8x8Blocks[59])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_60] = (uint16_t)((meanOf8x8SquaredValuesBlocks[60] - (meanOf8x8Blocks[60] * meanOf8x8Blocks[60])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_61] = (uint16_t)((meanOf8x8SquaredValuesBlocks[61] - (meanOf8x8Blocks[61] * meanOf8x8Blocks[61])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_62] = (uint16_t)((meanOf8x8SquaredValuesBlocks[62] - (meanOf8x8Blocks[62] * meanOf8x8Blocks[62])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_63] = (uint16_t)((meanOf8x8SquaredValuesBlocks[63] - (meanOf8x8Blocks[63] * meanOf8x8Blocks[63])) >> VARIANCE_PRECISION);

    // 16x16 variances
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_0] = (uint16_t)((meanOf16x16SquaredValuesBlocks[0] - (meanOf16x16Blocks[0] * meanOf16x16Blocks[0])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_1] = (uint16_t)((meanOf16x16SquaredValuesBlocks[1] - (meanOf16x16Blocks[1] * meanOf16x16Blocks[1])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_2] = (uint16_t)((meanOf16x16SquaredValuesBlocks[2] - (meanOf16x16Blocks[2] * meanOf16x16Blocks[2])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_3] = (uint16_t)((meanOf16x16SquaredValuesBlocks[3] - (meanOf16x16Blocks[3] * meanOf16x16Blocks[3])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_4] = (uint16_t)((meanOf16x16SquaredValuesBlocks[4] - (meanOf16x16Blocks[4] * meanOf16x16Blocks[4])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_5] = (uint16_t)((meanOf16x16SquaredValuesBlocks[5] - (meanOf16x16Blocks[5] * meanOf16x16Blocks[5])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_6] = (uint16_t)((meanOf16x16SquaredValuesBlocks[6] - (meanOf16x16Blocks[6] * meanOf16x16Blocks[6])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_7] = (uint16_t)((meanOf16x16SquaredValuesBlocks[7] - (meanOf16x16Blocks[7] * meanOf16x16Blocks[7])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_8] = (uint16_t)((meanOf16x16SquaredValuesBlocks[8] - (meanOf16x16Blocks[8] * meanOf16x16Blocks[8])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_9] = (uint16_t)((meanOf16x16SquaredValuesBlocks[9] - (meanOf16x16Blocks[9] * meanOf16x16Blocks[9])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_10] = (uint16_t)((meanOf16x16SquaredValuesBlocks[10] - (meanOf16x16Blocks[10] * meanOf16x16Blocks[10])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_11] = (uint16_t)((meanOf16x16SquaredValuesBlocks[11] - (meanOf16x16Blocks[11] * meanOf16x16Blocks[11])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_12] = (uint16_t)((meanOf16x16SquaredValuesBlocks[12] - (meanOf16x16Blocks[12] * meanOf16x16Blocks[12])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_13] = (uint16_t)((meanOf16x16SquaredValuesBlocks[13] - (meanOf16x16Blocks[13] * meanOf16x16Blocks[13])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_14] = (uint16_t)((meanOf16x16SquaredValuesBlocks[14] - (meanOf16x16Blocks[14] * meanOf16x16Blocks[14])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_15] = (uint16_t)((meanOf16x16SquaredValuesBlocks[15] - (meanOf16x16Blocks[15] * meanOf16x16Blocks[15])) >> VARIANCE_PRECISION);

    // 32x32 variances
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_32x32_0] = (uint16_t)((meanOf32x32SquaredValuesBlocks[0] - (meanOf32x32Blocks[0] * meanOf32x32Blocks[0])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_32x32_1] = (uint16_t)((meanOf32x32SquaredValuesBlocks[1] - (meanOf32x32Blocks[1] * meanOf32x32Blocks[1])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_32x32_2] = (uint16_t)((meanOf32x32SquaredValuesBlocks[2] - (meanOf32x32Blocks[2] * meanOf32x32Blocks[2])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_32x32_3] = (uint16_t)((meanOf32x32SquaredValuesBlocks[3] - (meanOf32x32Blocks[3] * meanOf32x32Blocks[3])) >> VARIANCE_PRECISION);

    // 64x64 variance
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_64x64] = (uint16_t)((meanOf64x64SquaredValuesBlocks - (meanOf64x64Blocks * meanOf64x64Blocks)) >> VARIANCE_PRECISION);

    return return_error;
}

EbErrorType DenoiseInputPicture(
    PictureAnalysisContext_t    *context_ptr,
    SequenceControlSet_t        *sequence_control_set_ptr,
    PictureParentControlSet_t   *picture_control_set_ptr,
    uint32_t                       sb_total_count,
    EbPictureBufferDesc_t       *inputPicturePtr,
    EbPictureBufferDesc_t       *denoisedPicturePtr,
    uint32_t                         picture_width_in_sb,
    EbAsm                         asm_type)
{
    EbErrorType return_error = EB_ErrorNone;

    uint32_t         lcuCodingOrder;
    uint32_t       sb_origin_x;
    uint32_t       sb_origin_y;
    uint16_t       verticalIdx;
    //use denoised input if the source is extremly noisy
    if (picture_control_set_ptr->pic_noise_class >= PIC_NOISE_CLASS_4) {

        uint32_t inLumaOffSet = inputPicturePtr->origin_x + inputPicturePtr->origin_y      * inputPicturePtr->strideY;
        uint32_t inChromaOffSet = inputPicturePtr->origin_x / 2 + inputPicturePtr->origin_y / 2 * inputPicturePtr->strideCb;
        uint32_t denLumaOffSet = denoisedPicturePtr->origin_x + denoisedPicturePtr->origin_y   * denoisedPicturePtr->strideY;
        uint32_t denChromaOffSet = denoisedPicturePtr->origin_x / 2 + denoisedPicturePtr->origin_y / 2 * denoisedPicturePtr->strideCb;

        //filter Luma
        for (lcuCodingOrder = 0; lcuCodingOrder < sb_total_count; ++lcuCodingOrder) {

            sb_origin_x = (lcuCodingOrder % picture_width_in_sb) * sequence_control_set_ptr->sb_sz;
            sb_origin_y = (lcuCodingOrder / picture_width_in_sb) * sequence_control_set_ptr->sb_sz;

            if (sb_origin_x == 0)
                StrongLumaFilter_funcPtrArray[asm_type](
                    inputPicturePtr,
                    denoisedPicturePtr,
                    sb_origin_y,
                    sb_origin_x);

            if (sb_origin_x + BLOCK_SIZE_64 > inputPicturePtr->width)
            {
                noiseExtractLumaStrong(
                    inputPicturePtr,
                    denoisedPicturePtr,
                    sb_origin_y,
                    sb_origin_x);
            }

        }

        //copy Luma
        for (verticalIdx = 0; verticalIdx < inputPicturePtr->height; ++verticalIdx) {
            EB_MEMCPY(inputPicturePtr->bufferY + inLumaOffSet + verticalIdx * inputPicturePtr->strideY,
                denoisedPicturePtr->bufferY + denLumaOffSet + verticalIdx * denoisedPicturePtr->strideY,
                sizeof(uint8_t) * inputPicturePtr->width);
        }

        //copy chroma
        for (lcuCodingOrder = 0; lcuCodingOrder < sb_total_count; ++lcuCodingOrder) {

            sb_origin_x = (lcuCodingOrder % picture_width_in_sb) * sequence_control_set_ptr->sb_sz;
            sb_origin_y = (lcuCodingOrder / picture_width_in_sb) * sequence_control_set_ptr->sb_sz;

            if (sb_origin_x == 0)
                StrongChromaFilter_funcPtrArray[asm_type](
                    inputPicturePtr,
                    denoisedPicturePtr,
                    sb_origin_y / 2,
                    sb_origin_x / 2);

            if (sb_origin_x + BLOCK_SIZE_64 > inputPicturePtr->width)
            {
                noiseExtractChromaStrong(
                    inputPicturePtr,
                    denoisedPicturePtr,
                    sb_origin_y / 2,
                    sb_origin_x / 2);
            }

        }

        //copy chroma
        for (verticalIdx = 0; verticalIdx < inputPicturePtr->height / 2; ++verticalIdx) {

            EB_MEMCPY(inputPicturePtr->bufferCb + inChromaOffSet + verticalIdx * inputPicturePtr->strideCb,
                denoisedPicturePtr->bufferCb + denChromaOffSet + verticalIdx * denoisedPicturePtr->strideCb,
                sizeof(uint8_t) * inputPicturePtr->width / 2);

            EB_MEMCPY(inputPicturePtr->bufferCr + inChromaOffSet + verticalIdx * inputPicturePtr->strideCr,
                denoisedPicturePtr->bufferCr + denChromaOffSet + verticalIdx * denoisedPicturePtr->strideCr,
                sizeof(uint8_t) * inputPicturePtr->width / 2);
        }

    }
    else if (picture_control_set_ptr->pic_noise_class >= PIC_NOISE_CLASS_3_1) {

        uint32_t inLumaOffSet = inputPicturePtr->origin_x + inputPicturePtr->origin_y      * inputPicturePtr->strideY;
        uint32_t inChromaOffSet = inputPicturePtr->origin_x / 2 + inputPicturePtr->origin_y / 2 * inputPicturePtr->strideCb;
        uint32_t denLumaOffSet = denoisedPicturePtr->origin_x + denoisedPicturePtr->origin_y   * denoisedPicturePtr->strideY;
        uint32_t denChromaOffSet = denoisedPicturePtr->origin_x / 2 + denoisedPicturePtr->origin_y / 2 * denoisedPicturePtr->strideCb;


        for (verticalIdx = 0; verticalIdx < inputPicturePtr->height; ++verticalIdx) {
            EB_MEMCPY(inputPicturePtr->bufferY + inLumaOffSet + verticalIdx * inputPicturePtr->strideY,
                denoisedPicturePtr->bufferY + denLumaOffSet + verticalIdx * denoisedPicturePtr->strideY,
                sizeof(uint8_t) * inputPicturePtr->width);
        }

        //copy chroma
        for (lcuCodingOrder = 0; lcuCodingOrder < sb_total_count; ++lcuCodingOrder) {

            sb_origin_x = (lcuCodingOrder % picture_width_in_sb) * sequence_control_set_ptr->sb_sz;
            sb_origin_y = (lcuCodingOrder / picture_width_in_sb) * sequence_control_set_ptr->sb_sz;

            if (sb_origin_x == 0)
                WeakChromaFilter_funcPtrArray[asm_type](
                    inputPicturePtr,
                    denoisedPicturePtr,
                    sb_origin_y / 2,
                    sb_origin_x / 2);

            if (sb_origin_x + BLOCK_SIZE_64 > inputPicturePtr->width)
            {
                noiseExtractChromaWeak(
                    inputPicturePtr,
                    denoisedPicturePtr,
                    sb_origin_y / 2,
                    sb_origin_x / 2);
            }

        }



        for (verticalIdx = 0; verticalIdx < inputPicturePtr->height / 2; ++verticalIdx) {

            EB_MEMCPY(inputPicturePtr->bufferCb + inChromaOffSet + verticalIdx * inputPicturePtr->strideCb,
                denoisedPicturePtr->bufferCb + denChromaOffSet + verticalIdx * denoisedPicturePtr->strideCb,
                sizeof(uint8_t) * inputPicturePtr->width / 2);

            EB_MEMCPY(inputPicturePtr->bufferCr + inChromaOffSet + verticalIdx * inputPicturePtr->strideCr,
                denoisedPicturePtr->bufferCr + denChromaOffSet + verticalIdx * denoisedPicturePtr->strideCr,
                sizeof(uint8_t) * inputPicturePtr->width / 2);
        }

    }
    else if (context_ptr->picNoiseVarianceFloat >= 1.0) {
        //Luma : use filtered only for flatNoise LCUs
        for (lcuCodingOrder = 0; lcuCodingOrder < sb_total_count; ++lcuCodingOrder) {

            sb_origin_x = (lcuCodingOrder % picture_width_in_sb) * sequence_control_set_ptr->sb_sz;
            sb_origin_y = (lcuCodingOrder / picture_width_in_sb) * sequence_control_set_ptr->sb_sz;
            uint32_t  sb_height = MIN(BLOCK_SIZE_64, inputPicturePtr->height - sb_origin_y);
            uint32_t  sb_width = MIN(BLOCK_SIZE_64, inputPicturePtr->width - sb_origin_x);

            uint32_t inLumaOffSet = inputPicturePtr->origin_x + sb_origin_x + (inputPicturePtr->origin_y + sb_origin_y) * inputPicturePtr->strideY;
            uint32_t denLumaOffSet = denoisedPicturePtr->origin_x + sb_origin_x + (denoisedPicturePtr->origin_y + sb_origin_y) * denoisedPicturePtr->strideY;


            if (picture_control_set_ptr->sb_flat_noise_array[lcuCodingOrder] == 1) {


                for (verticalIdx = 0; verticalIdx < sb_height; ++verticalIdx) {

                    EB_MEMCPY(inputPicturePtr->bufferY + inLumaOffSet + verticalIdx * inputPicturePtr->strideY,
                        denoisedPicturePtr->bufferY + denLumaOffSet + verticalIdx * denoisedPicturePtr->strideY,
                        sizeof(uint8_t) * sb_width);

                }
            }
        }
    }

    return return_error;
}

EbErrorType DetectInputPictureNoise(
    PictureAnalysisContext_t    *context_ptr,
    SequenceControlSet_t        *sequence_control_set_ptr,
    PictureParentControlSet_t   *picture_control_set_ptr,
    uint32_t                       sb_total_count,
    EbPictureBufferDesc_t       *inputPicturePtr,
    EbPictureBufferDesc_t       *noisePicturePtr,
    EbPictureBufferDesc_t       *denoisedPicturePtr,
    uint32_t                     picture_width_in_sb,
    EbAsm                        asm_type)
{

    EbErrorType                 return_error = EB_ErrorNone;
    uint32_t                    lcuCodingOrder;

    uint64_t                    picNoiseVariance;

    uint32_t                    totLcuCount, noiseTh;

    uint32_t                    sb_origin_x;
    uint32_t                    sb_origin_y;
    uint32_t                    inputLumaOriginIndex;

    picNoiseVariance = 0;
    totLcuCount = 0;


    //Variance calc for noise picture
    for (lcuCodingOrder = 0; lcuCodingOrder < sb_total_count; ++lcuCodingOrder) {

        sb_origin_x = (lcuCodingOrder % picture_width_in_sb) * sequence_control_set_ptr->sb_sz;
        sb_origin_y = (lcuCodingOrder / picture_width_in_sb) * sequence_control_set_ptr->sb_sz;
        inputLumaOriginIndex = (noisePicturePtr->origin_y + sb_origin_y) * noisePicturePtr->strideY +
            noisePicturePtr->origin_x + sb_origin_x;


        uint32_t  noiseOriginIndex = noisePicturePtr->origin_x + sb_origin_x + noisePicturePtr->origin_y * noisePicturePtr->strideY;

        if (sb_origin_x == 0)
            WeakLumaFilter_funcPtrArray[asm_type](
                inputPicturePtr,
                denoisedPicturePtr,
                noisePicturePtr,
                sb_origin_y,
                sb_origin_x);

        if (sb_origin_x + BLOCK_SIZE_64 > inputPicturePtr->width)
        {
            noiseExtractLumaWeak(
                inputPicturePtr,
                denoisedPicturePtr,
                noisePicturePtr,
                sb_origin_y,
                sb_origin_x);
        }



        //do it only for complete 64x64 blocks
        if (sb_origin_x + 64 <= inputPicturePtr->width && sb_origin_y + 64 <= inputPicturePtr->height)
        {

            uint64_t noiseBlkVar32x32[4], denoiseBlkVar32x32[4];

            uint64_t noiseBlkVar = ComputeVariance64x64(
                sequence_control_set_ptr,
                noisePicturePtr,
                noiseOriginIndex,
                noiseBlkVar32x32,
                asm_type);

            uint64_t noiseBlkVarTh;
            uint64_t denBlkVarTh = FLAT_MAX_VAR;
#if ENCODER_MODE_CLEANUP
            if(0)
#else
            if (picture_control_set_ptr->enc_mode > ENC_M3)
#endif
                noiseBlkVarTh = NOISE_MIN_LEVEL;
            else
                noiseBlkVarTh = NOISE_MIN_LEVEL_M6_M7;

            picNoiseVariance += (noiseBlkVar >> 16);

            uint64_t denBlkVar = ComputeVariance64x64(
                sequence_control_set_ptr,
                denoisedPicturePtr,
                inputLumaOriginIndex,
                denoiseBlkVar32x32,
                asm_type) >> 16;

            if (denBlkVar < denBlkVarTh && noiseBlkVar > noiseBlkVarTh) {
                picture_control_set_ptr->sb_flat_noise_array[lcuCodingOrder] = 1;
            }

            totLcuCount++;
        }

    }

    context_ptr->picNoiseVarianceFloat = (double)picNoiseVariance / (double)totLcuCount;

    picNoiseVariance = picNoiseVariance / totLcuCount;


    //the variance of a 64x64 noise area tends to be bigger for small resolutions.
    if (sequence_control_set_ptr->luma_height <= 720)
        noiseTh = 25;
    else
        noiseTh = 0;

    if (picNoiseVariance >= 80 + noiseTh)
        picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_10;
    else if (picNoiseVariance >= 70 + noiseTh)
        picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_9;
    else if (picNoiseVariance >= 60 + noiseTh)
        picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_8;
    else if (picNoiseVariance >= 50 + noiseTh)
        picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_7;
    else if (picNoiseVariance >= 40 + noiseTh)
        picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_6;
    else if (picNoiseVariance >= 30 + noiseTh)
        picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_5;
    else if (picNoiseVariance >= 20 + noiseTh)
        picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_4;
    else if (picNoiseVariance >= 17 + noiseTh)
        picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_3_1;
    else if (picNoiseVariance >= 10 + noiseTh)
        picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_3;
    else if (picNoiseVariance >= 5 + noiseTh)
        picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_2;
    else
        picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_1;

    if (picture_control_set_ptr->pic_noise_class >= PIC_NOISE_CLASS_4)
        picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_3_1;

    return return_error;

}

static int32_t apply_denoise_2d(SequenceControlSet_t        *scsPtr,
    PictureParentControlSet_t   *pcsPtr,
    EbPictureBufferDesc_t *inputPicturePointer,
    EbAsm asm_type) {


    if (aom_denoise_and_model_run(pcsPtr->denoise_and_model, inputPicturePointer,
        &pcsPtr->film_grain_params,
        scsPtr->static_config.encoder_bit_depth > EB_8BIT, asm_type)) {
    }
    return 0;
}


EbErrorType denoise_estimate_film_grain(
    SequenceControlSet_t        *sequence_control_set_ptr,
    PictureParentControlSet_t   *picture_control_set_ptr,
    EbAsm asm_type)
{
    EbErrorType return_error = EB_ErrorNone;

    EbPictureBufferDesc_t    *inputPicturePtr = picture_control_set_ptr->enhanced_picture_ptr;
    picture_control_set_ptr->film_grain_params.apply_grain = 0;

    if (sequence_control_set_ptr->film_grain_denoise_strength) {
        if (apply_denoise_2d(sequence_control_set_ptr, picture_control_set_ptr, inputPicturePtr, asm_type) < 0)
            return 1;
    }

    sequence_control_set_ptr->film_grain_params_present |= picture_control_set_ptr->film_grain_params.apply_grain;

    return return_error;  //todo: add proper error handling
}


EbErrorType FullSampleDenoise(
    PictureAnalysisContext_t    *context_ptr,
    SequenceControlSet_t        *sequence_control_set_ptr,
    PictureParentControlSet_t   *picture_control_set_ptr,
    uint32_t                     sb_total_count,
    EbBool                       denoiseFlag,
    uint32_t                     picture_width_in_sb,
    EbAsm                        asm_type){

    EbErrorType return_error = EB_ErrorNone;

    uint32_t                     lcuCodingOrder;
    EbPictureBufferDesc_t    *inputPicturePtr = picture_control_set_ptr->enhanced_picture_ptr;
    EbPictureBufferDesc_t    *denoisedPicturePtr = context_ptr->denoisedPicturePtr;
    EbPictureBufferDesc_t    *noisePicturePtr = context_ptr->noisePicturePtr;

    //Reset the flat noise flag array to False for both RealTime/HighComplexity Modes
    for (lcuCodingOrder = 0; lcuCodingOrder < sb_total_count; ++lcuCodingOrder) {
        picture_control_set_ptr->sb_flat_noise_array[lcuCodingOrder] = 0;
    }

    picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_INV; //this init is for both REAL-TIME and BEST-QUALITY

    DetectInputPictureNoise(
        context_ptr,
        sequence_control_set_ptr,
        picture_control_set_ptr,
        sb_total_count,
        inputPicturePtr,
        noisePicturePtr,
        denoisedPicturePtr,
        picture_width_in_sb,
        asm_type);

    if (denoiseFlag == EB_TRUE)
    {

        DenoiseInputPicture(
            context_ptr,
            sequence_control_set_ptr,
            picture_control_set_ptr,
            sb_total_count,
            inputPicturePtr,
            denoisedPicturePtr,
            picture_width_in_sb,
            asm_type);
    }

    return return_error;

}

EbErrorType SubSampleFilterNoise(
    SequenceControlSet_t        *sequence_control_set_ptr,
    PictureParentControlSet_t   *picture_control_set_ptr,
    uint32_t                       sb_total_count,
    EbPictureBufferDesc_t       *inputPicturePtr,
    EbPictureBufferDesc_t       *noisePicturePtr,
    EbPictureBufferDesc_t       *denoisedPicturePtr,
    uint32_t                         picture_width_in_sb,
    EbAsm                         asm_type)
{
    EbErrorType return_error = EB_ErrorNone;

    uint32_t         lcuCodingOrder;
    uint32_t       sb_origin_x;
    uint32_t       sb_origin_y;
    uint16_t        verticalIdx;
    if (picture_control_set_ptr->pic_noise_class == PIC_NOISE_CLASS_3_1) {

        uint32_t inLumaOffSet = inputPicturePtr->origin_x + inputPicturePtr->origin_y      * inputPicturePtr->strideY;
        uint32_t inChromaOffSet = inputPicturePtr->origin_x / 2 + inputPicturePtr->origin_y / 2 * inputPicturePtr->strideCb;
        uint32_t denLumaOffSet = denoisedPicturePtr->origin_x + denoisedPicturePtr->origin_y   * denoisedPicturePtr->strideY;
        uint32_t denChromaOffSet = denoisedPicturePtr->origin_x / 2 + denoisedPicturePtr->origin_y / 2 * denoisedPicturePtr->strideCb;


        //filter Luma
        for (lcuCodingOrder = 0; lcuCodingOrder < sb_total_count; ++lcuCodingOrder) {

            sb_origin_x = (lcuCodingOrder % picture_width_in_sb) * sequence_control_set_ptr->sb_sz;
            sb_origin_y = (lcuCodingOrder / picture_width_in_sb) * sequence_control_set_ptr->sb_sz;

            if (sb_origin_x == 0)
                WeakLumaFilter_funcPtrArray[asm_type](
                    inputPicturePtr,
                    denoisedPicturePtr,
                    noisePicturePtr,
                    sb_origin_y,
                    sb_origin_x);

            if (sb_origin_x + BLOCK_SIZE_64 > inputPicturePtr->width)
            {
                noiseExtractLumaWeak(
                    inputPicturePtr,
                    denoisedPicturePtr,
                    noisePicturePtr,
                    sb_origin_y,
                    sb_origin_x);
            }
        }

        //copy luma
        for (verticalIdx = 0; verticalIdx < inputPicturePtr->height; ++verticalIdx) {
            EB_MEMCPY(inputPicturePtr->bufferY + inLumaOffSet + verticalIdx * inputPicturePtr->strideY,
                denoisedPicturePtr->bufferY + denLumaOffSet + verticalIdx * denoisedPicturePtr->strideY,
                sizeof(uint8_t) * inputPicturePtr->width);
        }

        //filter chroma
        for (lcuCodingOrder = 0; lcuCodingOrder < sb_total_count; ++lcuCodingOrder) {

            sb_origin_x = (lcuCodingOrder % picture_width_in_sb) * sequence_control_set_ptr->sb_sz;
            sb_origin_y = (lcuCodingOrder / picture_width_in_sb) * sequence_control_set_ptr->sb_sz;

            if (sb_origin_x == 0)
                WeakChromaFilter_funcPtrArray[asm_type](
                    inputPicturePtr,
                    denoisedPicturePtr,
                    sb_origin_y / 2,
                    sb_origin_x / 2);

            if (sb_origin_x + BLOCK_SIZE_64 > inputPicturePtr->width)
            {
                noiseExtractChromaWeak(
                    inputPicturePtr,
                    denoisedPicturePtr,
                    sb_origin_y / 2,
                    sb_origin_x / 2);
            }

        }

        //copy chroma
        for (verticalIdx = 0; verticalIdx < inputPicturePtr->height / 2; ++verticalIdx) {

            EB_MEMCPY(inputPicturePtr->bufferCb + inChromaOffSet + verticalIdx * inputPicturePtr->strideCb,
                denoisedPicturePtr->bufferCb + denChromaOffSet + verticalIdx * denoisedPicturePtr->strideCb,
                sizeof(uint8_t) * inputPicturePtr->width / 2);

            EB_MEMCPY(inputPicturePtr->bufferCr + inChromaOffSet + verticalIdx * inputPicturePtr->strideCr,
                denoisedPicturePtr->bufferCr + denChromaOffSet + verticalIdx * denoisedPicturePtr->strideCr,
                sizeof(uint8_t) * inputPicturePtr->width / 2);
        }

    }
    else if (picture_control_set_ptr->pic_noise_class == PIC_NOISE_CLASS_2) {

        uint32_t newTotFN = 0;

        //for each SB ,re check the FN information for only the FNdecim ones
        for (lcuCodingOrder = 0; lcuCodingOrder < sb_total_count; ++lcuCodingOrder) {

            sb_origin_x = (lcuCodingOrder % picture_width_in_sb) * sequence_control_set_ptr->sb_sz;
            sb_origin_y = (lcuCodingOrder / picture_width_in_sb) * sequence_control_set_ptr->sb_sz;
            uint32_t  inputLumaOriginIndex = noisePicturePtr->origin_x + sb_origin_x + (noisePicturePtr->origin_y + sb_origin_y) * noisePicturePtr->strideY;
            uint32_t  noiseOriginIndex = noisePicturePtr->origin_x + sb_origin_x + (noisePicturePtr->origin_y * noisePicturePtr->strideY);

            if (sb_origin_x + 64 <= inputPicturePtr->width && sb_origin_y + 64 <= inputPicturePtr->height && picture_control_set_ptr->sb_flat_noise_array[lcuCodingOrder] == 1)
            {

                WeakLumaFilterLcu_funcPtrArray[asm_type](
                    inputPicturePtr,
                    denoisedPicturePtr,
                    noisePicturePtr,
                    sb_origin_y,
                    sb_origin_x);

                if (sb_origin_x + BLOCK_SIZE_64 > inputPicturePtr->width)
                {
                    noiseExtractLumaWeakLcu(
                        inputPicturePtr,
                        denoisedPicturePtr,
                        noisePicturePtr,
                        sb_origin_y,
                        sb_origin_x);
                }

                uint64_t noiseBlkVar32x32[4], denoiseBlkVar32x32[4];
                uint64_t noiseBlkVar = ComputeVariance64x64(
                    sequence_control_set_ptr,
                    noisePicturePtr,
                    noiseOriginIndex,
                    noiseBlkVar32x32,
                    asm_type);
                uint64_t denBlkVar = ComputeVariance64x64(
                    sequence_control_set_ptr,
                    denoisedPicturePtr,
                    inputLumaOriginIndex,
                    denoiseBlkVar32x32,
                    asm_type) >> 16;

                uint64_t noiseBlkVarTh;
                uint64_t denBlkVarTh = FLAT_MAX_VAR;
#if ENCODER_MODE_CLEANUP
                if (0)
#else
                if (picture_control_set_ptr->enc_mode > ENC_M3)
#endif
                    noiseBlkVarTh = NOISE_MIN_LEVEL;
                else
                    noiseBlkVarTh = NOISE_MIN_LEVEL_M6_M7;

                if (denBlkVar<denBlkVarTh && noiseBlkVar> noiseBlkVarTh) {
                    picture_control_set_ptr->sb_flat_noise_array[lcuCodingOrder] = 1;
                    //printf("POC %i (%i,%i) denBlkVar: %i  noiseBlkVar :%i\n", picture_control_set_ptr->picture_number,sb_origin_x,sb_origin_y, denBlkVar, noiseBlkVar);
                    newTotFN++;

                }
                else {
                    picture_control_set_ptr->sb_flat_noise_array[lcuCodingOrder] = 0;
                }
            }
        }

        for (lcuCodingOrder = 0; lcuCodingOrder < sb_total_count; ++lcuCodingOrder) {

            sb_origin_x = (lcuCodingOrder % picture_width_in_sb) * sequence_control_set_ptr->sb_sz;
            sb_origin_y = (lcuCodingOrder / picture_width_in_sb) * sequence_control_set_ptr->sb_sz;

            if (sb_origin_x + 64 <= inputPicturePtr->width && sb_origin_y + 64 <= inputPicturePtr->height)
            {


                //use the denoised for FN LCUs
                if (picture_control_set_ptr->sb_flat_noise_array[lcuCodingOrder] == 1) {

                    uint32_t  sb_height = MIN(BLOCK_SIZE_64, inputPicturePtr->height - sb_origin_y);
                    uint32_t  sb_width = MIN(BLOCK_SIZE_64, inputPicturePtr->width - sb_origin_x);

                    uint32_t inLumaOffSet = inputPicturePtr->origin_x + sb_origin_x + (inputPicturePtr->origin_y + sb_origin_y) * inputPicturePtr->strideY;
                    uint32_t denLumaOffSet = denoisedPicturePtr->origin_x + sb_origin_x + (denoisedPicturePtr->origin_y + sb_origin_y) * denoisedPicturePtr->strideY;

                    for (verticalIdx = 0; verticalIdx < sb_height; ++verticalIdx) {

                        EB_MEMCPY(inputPicturePtr->bufferY + inLumaOffSet + verticalIdx * inputPicturePtr->strideY,
                            denoisedPicturePtr->bufferY + denLumaOffSet + verticalIdx * denoisedPicturePtr->strideY,
                            sizeof(uint8_t) * sb_width);

                    }
                }

            }

        }

    }
    return return_error;
}

EbErrorType QuarterSampleDetectNoise(
    PictureAnalysisContext_t    *context_ptr,
    PictureParentControlSet_t   *picture_control_set_ptr,
    EbPictureBufferDesc_t       *quarterDecimatedPicturePtr,
    EbPictureBufferDesc_t       *noisePicturePtr,
    EbPictureBufferDesc_t       *denoisedPicturePtr,
    uint32_t                         picture_width_in_sb,
    EbAsm                         asm_type)
{

    EbErrorType return_error = EB_ErrorNone;

    uint64_t                   picNoiseVariance;

    uint32_t                     totLcuCount, noiseTh;

    uint32_t                     blockIndex;

    picNoiseVariance = 0;
    totLcuCount = 0;


    uint16_t vert64x64Index;
    uint16_t horz64x64Index;
    uint32_t block64x64X;
    uint32_t block64x64Y;
    uint32_t vert32x32Index;
    uint32_t horz32x32Index;
    uint32_t block32x32X;
    uint32_t block32x32Y;
    uint32_t noiseOriginIndex;
    uint32_t lcuCodingOrder;

    // Loop over 64x64 blocks on the downsampled domain (each block would contain 16 LCUs on the full sampled domain)
    for (vert64x64Index = 0; vert64x64Index < (quarterDecimatedPicturePtr->height / 64); vert64x64Index++) {
        for (horz64x64Index = 0; horz64x64Index < (quarterDecimatedPicturePtr->width / 64); horz64x64Index++) {

            block64x64X = horz64x64Index * 64;
            block64x64Y = vert64x64Index * 64;

            if (block64x64X == 0)
                WeakLumaFilter_funcPtrArray[asm_type](
                    quarterDecimatedPicturePtr,
                    denoisedPicturePtr,
                    noisePicturePtr,
                    block64x64Y,
                    block64x64X);

            if (block64x64Y + BLOCK_SIZE_64 > quarterDecimatedPicturePtr->width)
            {
                noiseExtractLumaWeak(
                    quarterDecimatedPicturePtr,
                    denoisedPicturePtr,
                    noisePicturePtr,
                    block64x64Y,
                    block64x64X);
            }


            // Loop over 32x32 blocks (i.e, 64x64 blocks in full resolution)
            for (vert32x32Index = 0; vert32x32Index < 2; vert32x32Index++) {
                for (horz32x32Index = 0; horz32x32Index < 2; horz32x32Index++) {

                    block32x32X = block64x64X + horz32x32Index * 32;
                    block32x32Y = block64x64Y + vert32x32Index * 32;

                    //do it only for complete 32x32 blocks (i.e, complete 64x64 blocks in full resolution)
                    if ((block32x32X + 32 <= quarterDecimatedPicturePtr->width) && (block32x32Y + 32 <= quarterDecimatedPicturePtr->height))
                    {

                        lcuCodingOrder = ((vert64x64Index * 2) + vert32x32Index) * picture_width_in_sb + ((horz64x64Index * 2) + horz32x32Index);


                        uint64_t noiseBlkVar8x8[16], denoiseBlkVar8x8[16];

                        noiseOriginIndex = noisePicturePtr->origin_x + block32x32X + noisePicturePtr->origin_y * noisePicturePtr->strideY;

                        uint64_t noiseBlkVar = ComputeVariance32x32(
                            noisePicturePtr,
                            noiseOriginIndex,
                            noiseBlkVar8x8,
                            asm_type);


                        picNoiseVariance += (noiseBlkVar >> 16);

                        blockIndex = (noisePicturePtr->origin_y + block32x32Y) * noisePicturePtr->strideY + noisePicturePtr->origin_x + block32x32X;

                        uint64_t denBlkVar = ComputeVariance32x32(
                            denoisedPicturePtr,
                            blockIndex,
                            denoiseBlkVar8x8,
                            asm_type) >> 16;

                        uint64_t denBlkVarDecTh;
#if ENCODER_MODE_CLEANUP
                        if (0){
#else
                        if (picture_control_set_ptr->enc_mode > ENC_M3) {
#endif
                            denBlkVarDecTh = NOISE_MIN_LEVEL_DECIM;
                        }
                        else {
                            denBlkVarDecTh = NOISE_MIN_LEVEL_DECIM_M6_M7;
                        }
                        if (denBlkVar < FLAT_MAX_VAR_DECIM && noiseBlkVar> denBlkVarDecTh) {
                            picture_control_set_ptr->sb_flat_noise_array[lcuCodingOrder] = 1;
                        }

                        totLcuCount++;
                    }
                }
            }
        }
    }


    context_ptr->picNoiseVarianceFloat = (double)picNoiseVariance / (double)totLcuCount;

    picNoiseVariance = picNoiseVariance / totLcuCount;


    //the variance of a 64x64 noise area tends to be bigger for small resolutions.
    //if (sequence_control_set_ptr->luma_height <= 720)
    //    noiseTh = 25;
    //else if (sequence_control_set_ptr->luma_height <= 1080)
    //    noiseTh = 10;
    //else
    noiseTh = 0;

    //look for extreme noise or big enough flat noisy area to be denoised.
    if (picNoiseVariance > 60)
        picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_3_1; //Noise+Edge information is too big, so may be this is all noise (action: frame based denoising)
    else if (picNoiseVariance >= 10 + noiseTh)
        picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_3;   //Noise+Edge information is big enough, so there is no big enough flat noisy area (action : no denoising)
    else if (picNoiseVariance >= 5 + noiseTh)
        picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_2;   //Noise+Edge information is relatively small, so there might be a big enough flat noisy area(action : denoising only for FN blocks)
    else
        picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_1;   //Noise+Edge information is very small, so no noise nor edge area (action : no denoising)



    return return_error;

}



EbErrorType SubSampleDetectNoise(
    PictureAnalysisContext_t    *context_ptr,
    SequenceControlSet_t        *sequence_control_set_ptr,
    PictureParentControlSet_t   *picture_control_set_ptr,
    EbPictureBufferDesc_t       *sixteenthDecimatedPicturePtr,
    EbPictureBufferDesc_t       *noisePicturePtr,
    EbPictureBufferDesc_t       *denoisedPicturePtr,
    uint32_t                         picture_width_in_sb,
    EbAsm                         asm_type)
{

    EbErrorType return_error = EB_ErrorNone;

    uint64_t                   picNoiseVariance;

    uint32_t                     totLcuCount, noiseTh;

    uint32_t                     blockIndex;

    picNoiseVariance = 0;
    totLcuCount = 0;


    uint16_t vert64x64Index;
    uint16_t horz64x64Index;
    uint32_t block64x64X;
    uint32_t block64x64Y;
    uint32_t vert16x16Index;
    uint32_t horz16x16Index;
    uint32_t block16x16X;
    uint32_t block16x16Y;
    uint32_t noiseOriginIndex;
    uint32_t lcuCodingOrder;

    // Loop over 64x64 blocks on the downsampled domain (each block would contain 16 LCUs on the full sampled domain)
    for (vert64x64Index = 0; vert64x64Index < (sixteenthDecimatedPicturePtr->height / 64); vert64x64Index++) {
        for (horz64x64Index = 0; horz64x64Index < (sixteenthDecimatedPicturePtr->width / 64); horz64x64Index++) {

            block64x64X = horz64x64Index * 64;
            block64x64Y = vert64x64Index * 64;

            if (block64x64X == 0)
                WeakLumaFilter_funcPtrArray[asm_type](
                    sixteenthDecimatedPicturePtr,
                    denoisedPicturePtr,
                    noisePicturePtr,
                    block64x64Y,
                    block64x64X);

            if (block64x64Y + BLOCK_SIZE_64 > sixteenthDecimatedPicturePtr->width)
            {
                noiseExtractLumaWeak(
                    sixteenthDecimatedPicturePtr,
                    denoisedPicturePtr,
                    noisePicturePtr,
                    block64x64Y,
                    block64x64X);
            }


            // Loop over 16x16 blocks (i.e, 64x64 blocks in full resolution)
            for (vert16x16Index = 0; vert16x16Index < 4; vert16x16Index++) {
                for (horz16x16Index = 0; horz16x16Index < 4; horz16x16Index++) {

                    block16x16X = block64x64X + horz16x16Index * 16;
                    block16x16Y = block64x64Y + vert16x16Index * 16;

                    //do it only for complete 16x16 blocks (i.e, complete 64x64 blocks in full resolution)
                    if (block16x16X + 16 <= sixteenthDecimatedPicturePtr->width && block16x16Y + 16 <= sixteenthDecimatedPicturePtr->height)
                    {

                        lcuCodingOrder = ((vert64x64Index * 4) + vert16x16Index) * picture_width_in_sb + ((horz64x64Index * 4) + horz16x16Index);


                        uint64_t noiseBlkVar8x8[4], denoiseBlkVar8x8[4];

                        noiseOriginIndex = noisePicturePtr->origin_x + block16x16X + noisePicturePtr->origin_y * noisePicturePtr->strideY;

                        uint64_t noiseBlkVar = ComputeVariance16x16(
                            noisePicturePtr,
                            noiseOriginIndex,
                            noiseBlkVar8x8,
                            asm_type);


                        picNoiseVariance += (noiseBlkVar >> 16);

                        blockIndex = (noisePicturePtr->origin_y + block16x16Y) * noisePicturePtr->strideY + noisePicturePtr->origin_x + block16x16X;

                        uint64_t denBlkVar = ComputeVariance16x16(
                            denoisedPicturePtr,
                            blockIndex,
                            denoiseBlkVar8x8,
                            asm_type) >> 16;

                        uint64_t  noiseBlkVarDecTh;
                        uint64_t denBlkVarDecTh = FLAT_MAX_VAR_DECIM;
#if ENCODER_MODE_CLEANUP
                        if (0){
#else
                        if (picture_control_set_ptr->enc_mode > ENC_M3) {
#endif
                            noiseBlkVarDecTh = NOISE_MIN_LEVEL_DECIM;
                        }
                        else {
                            noiseBlkVarDecTh = NOISE_MIN_LEVEL_DECIM_M6_M7;
                        }

                        if (denBlkVar < denBlkVarDecTh && noiseBlkVar> noiseBlkVarDecTh) {
                            picture_control_set_ptr->sb_flat_noise_array[lcuCodingOrder] = 1;
                        }
                        totLcuCount++;
                    }
                }
            }
        }
    }


    context_ptr->picNoiseVarianceFloat = (double)picNoiseVariance / (double)totLcuCount;

    picNoiseVariance = picNoiseVariance / totLcuCount;


    //the variance of a 64x64 noise area tends to be bigger for small resolutions.
    if (sequence_control_set_ptr->luma_height <= 720)
        noiseTh = 25;
    else if (sequence_control_set_ptr->luma_height <= 1080)
        noiseTh = 10;
    else
        noiseTh = 0;

    //look for extreme noise or big enough flat noisy area to be denoised.
    if (picNoiseVariance >= 55 + noiseTh)
        picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_3_1; //Noise+Edge information is too big, so may be this is all noise (action: frame based denoising)
    else if (picNoiseVariance >= 10 + noiseTh)
        picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_3;   //Noise+Edge information is big enough, so there is no big enough flat noisy area (action : no denoising)
    else if (picNoiseVariance >= 5 + noiseTh)
        picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_2;   //Noise+Edge information is relatively small, so there might be a big enough flat noisy area(action : denoising only for FN blocks)
    else
        picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_1;   //Noise+Edge information is very small, so no noise nor edge area (action : no denoising)

    return return_error;

}

EbErrorType QuarterSampleDenoise(
    PictureAnalysisContext_t    *context_ptr,
    SequenceControlSet_t        *sequence_control_set_ptr,
    PictureParentControlSet_t   *picture_control_set_ptr,
    EbPictureBufferDesc_t        *quarterDecimatedPicturePtr,
    uint32_t                       sb_total_count,
    EbBool                      denoiseFlag,
    uint32_t                         picture_width_in_sb,
    EbAsm                         asm_type)
{

    EbErrorType return_error = EB_ErrorNone;

    uint32_t                     lcuCodingOrder;
    EbPictureBufferDesc_t    *inputPicturePtr = picture_control_set_ptr->enhanced_picture_ptr;
    EbPictureBufferDesc_t    *denoisedPicturePtr = context_ptr->denoisedPicturePtr;
    EbPictureBufferDesc_t    *noisePicturePtr = context_ptr->noisePicturePtr;

    //Reset the flat noise flag array to False for both RealTime/HighComplexity Modes
    for (lcuCodingOrder = 0; lcuCodingOrder < sb_total_count; ++lcuCodingOrder) {
        picture_control_set_ptr->sb_flat_noise_array[lcuCodingOrder] = 0;
    }

    picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_INV; //this init is for both REAL-TIME and BEST-QUALITY

    Decimation2D(
        &inputPicturePtr->bufferY[inputPicturePtr->origin_x + inputPicturePtr->origin_y * inputPicturePtr->strideY],
        inputPicturePtr->strideY,
        inputPicturePtr->width,
        inputPicturePtr->height,
        &quarterDecimatedPicturePtr->bufferY[quarterDecimatedPicturePtr->origin_x + (quarterDecimatedPicturePtr->origin_y * quarterDecimatedPicturePtr->strideY)],
        quarterDecimatedPicturePtr->strideY,
        2);


    QuarterSampleDetectNoise(
        context_ptr,
        picture_control_set_ptr,
        quarterDecimatedPicturePtr,
        noisePicturePtr,
        denoisedPicturePtr,
        picture_width_in_sb,
        asm_type);

    if (denoiseFlag == EB_TRUE) {

        // Turn OFF the de-noiser for Class 2 at QP=29 and lower (for Fixed_QP) and at the target rate of 14Mbps and higher (for RC=ON)
        if ((picture_control_set_ptr->pic_noise_class == PIC_NOISE_CLASS_3_1) ||
            ((picture_control_set_ptr->pic_noise_class == PIC_NOISE_CLASS_2) && ((sequence_control_set_ptr->static_config.rate_control_mode == 0 && sequence_control_set_ptr->qp > DENOISER_QP_TH) || (sequence_control_set_ptr->static_config.rate_control_mode != 0 && sequence_control_set_ptr->static_config.target_bit_rate < DENOISER_BITRATE_TH)))) {

            SubSampleFilterNoise(
                sequence_control_set_ptr,
                picture_control_set_ptr,
                sb_total_count,
                inputPicturePtr,
                noisePicturePtr,
                denoisedPicturePtr,
                picture_width_in_sb,
                asm_type);
        }
    }

    return return_error;

}


EbErrorType SubSampleDenoise(
    PictureAnalysisContext_t    *context_ptr,
    SequenceControlSet_t        *sequence_control_set_ptr,
    PictureParentControlSet_t   *picture_control_set_ptr,
    EbPictureBufferDesc_t        *sixteenthDecimatedPicturePtr,
    uint32_t                       sb_total_count,
    EbBool                      denoiseFlag,
    uint32_t                         picture_width_in_sb,
    EbAsm                         asm_type)
{

    EbErrorType return_error = EB_ErrorNone;

    uint32_t                     lcuCodingOrder;
    EbPictureBufferDesc_t    *inputPicturePtr = picture_control_set_ptr->enhanced_picture_ptr;
    EbPictureBufferDesc_t    *denoisedPicturePtr = context_ptr->denoisedPicturePtr;
    EbPictureBufferDesc_t    *noisePicturePtr = context_ptr->noisePicturePtr;

    //Reset the flat noise flag array to False for both RealTime/HighComplexity Modes
    for (lcuCodingOrder = 0; lcuCodingOrder < sb_total_count; ++lcuCodingOrder) {
        picture_control_set_ptr->sb_flat_noise_array[lcuCodingOrder] = 0;
    }

    picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_INV; //this init is for both REAL-TIME and BEST-QUALITY

    Decimation2D(
        &inputPicturePtr->bufferY[inputPicturePtr->origin_x + inputPicturePtr->origin_y * inputPicturePtr->strideY],
        inputPicturePtr->strideY,
        inputPicturePtr->width,
        inputPicturePtr->height,
        &sixteenthDecimatedPicturePtr->bufferY[sixteenthDecimatedPicturePtr->origin_x + (sixteenthDecimatedPicturePtr->origin_y * sixteenthDecimatedPicturePtr->strideY)],
        sixteenthDecimatedPicturePtr->strideY,
        4);

    SubSampleDetectNoise(
        context_ptr,
        sequence_control_set_ptr,
        picture_control_set_ptr,
        sixteenthDecimatedPicturePtr,
        noisePicturePtr,
        denoisedPicturePtr,
        picture_width_in_sb,
        asm_type);

    if (denoiseFlag == EB_TRUE) {

        // Turn OFF the de-noiser for Class 2 at QP=29 and lower (for Fixed_QP) and at the target rate of 14Mbps and higher (for RC=ON)
        if ((picture_control_set_ptr->pic_noise_class == PIC_NOISE_CLASS_3_1) ||
            ((picture_control_set_ptr->pic_noise_class == PIC_NOISE_CLASS_2) && ((sequence_control_set_ptr->static_config.rate_control_mode == 0 && sequence_control_set_ptr->qp > DENOISER_QP_TH) || (sequence_control_set_ptr->static_config.rate_control_mode != 0 && sequence_control_set_ptr->static_config.target_bit_rate < DENOISER_BITRATE_TH)))) {

            SubSampleFilterNoise(
                sequence_control_set_ptr,
                picture_control_set_ptr,
                sb_total_count,
                inputPicturePtr,
                noisePicturePtr,
                denoisedPicturePtr,
                picture_width_in_sb,
                asm_type);
        }
    }

    return return_error;

}


/************************************************
 * Set Picture Parameters based on input configuration
 ** Setting Number of regions per resolution
 ** Setting width and height for subpicture and when picture scan type is 1
 ************************************************/
void SetPictureParametersForStatisticsGathering(
    SequenceControlSet_t            *sequence_control_set_ptr
)
{

    sequence_control_set_ptr->picture_analysis_number_of_regions_per_width = HIGHER_THAN_CLASS_1_REGION_SPLIT_PER_WIDTH;
    sequence_control_set_ptr->picture_analysis_number_of_regions_per_height = HIGHER_THAN_CLASS_1_REGION_SPLIT_PER_HEIGHT;

    return;
}
/************************************************
 * Picture Pre Processing Operations *
 *** A function that groups all of the Pre proceesing
 * operations performed on the input picture
 *** Operations included at this point:
 ***** Borders preprocessing
 ***** Denoising
 ************************************************/
void PicturePreProcessingOperations(
    PictureParentControlSet_t       *picture_control_set_ptr,
    EbPictureBufferDesc_t           *inputPicturePtr,
    PictureAnalysisContext_t        *context_ptr,
    SequenceControlSet_t            *sequence_control_set_ptr,
    EbPictureBufferDesc_t           *quarterDecimatedPicturePtr,
    EbPictureBufferDesc_t           *sixteenthDecimatedPicturePtr,
    uint32_t                           sb_total_count,
    uint32_t                           picture_width_in_sb,
    EbAsm                           asm_type) {

    (void)inputPicturePtr;

    if (sequence_control_set_ptr->film_grain_denoise_strength) {

        denoise_estimate_film_grain(
            sequence_control_set_ptr,
            picture_control_set_ptr,
            asm_type);
    }
    // Use the decimated based detector for only for 4K.
    // When the decimated based detector is used and the enable noise flag is false, do not skip the detection since it is not costly
    // When the decimated based detector is used, turn OFF the de-noiser for Class 2 at QP=29 and lower (for Fixed_QP) and at the target rate of 14Mbps and higher (for RC=ON)
#if ENCODER_MODE_CLEANUP
    else if (0){
#else
    else if (sequence_control_set_ptr->input_resolution == INPUT_SIZE_4K_RANGE && picture_control_set_ptr->enc_mode > ENC_M3) {
#endif
        SubSampleDenoise(
            context_ptr,
            sequence_control_set_ptr,
            picture_control_set_ptr,
            sixteenthDecimatedPicturePtr,
            sb_total_count,
            sequence_control_set_ptr->static_config.enable_denoise_flag,
            picture_width_in_sb,
            asm_type);
    }
#if ENCODER_MODE_CLEANUP
    else if (0) {
#else
    else if (sequence_control_set_ptr->input_resolution == INPUT_SIZE_1080p_RANGE && picture_control_set_ptr->enc_mode > ENC_M3) {
#endif
        QuarterSampleDenoise(
            context_ptr,
            sequence_control_set_ptr,
            picture_control_set_ptr,
            quarterDecimatedPicturePtr,
            sb_total_count,
            sequence_control_set_ptr->static_config.enable_denoise_flag,
            picture_width_in_sb,
            asm_type);
    }
    else {

        FullSampleDenoise(
            context_ptr,
            sequence_control_set_ptr,
            picture_control_set_ptr,
            sb_total_count,
            sequence_control_set_ptr->static_config.enable_denoise_flag,
            picture_width_in_sb,
            asm_type);
    }
    return;

}

/**************************************************************
* Generate picture histogram bins for YUV pixel intensity *
* Calculation is done on a region based (Set previously, resolution dependent)
**************************************************************/
void SubSampleLumaGeneratePixelIntensityHistogramBins(
    SequenceControlSet_t            *sequence_control_set_ptr,
    PictureParentControlSet_t       *picture_control_set_ptr,
    EbPictureBufferDesc_t           *inputPicturePtr,
    uint64_t                          *sumAverageIntensityTotalRegionsLuma,
    EbAsm                           asm_type) {

    uint32_t                          regionWidth;
    uint32_t                          regionHeight;
    uint32_t                          regionWidthOffset;
    uint32_t                          regionHeightOffset;
    uint32_t                          regionInPictureWidthIndex;
    uint32_t                          regionInPictureHeightIndex;
    uint32_t                            histogramBin;
    uint64_t                          sum;

    regionWidth = inputPicturePtr->width / sequence_control_set_ptr->picture_analysis_number_of_regions_per_width;
    regionHeight = inputPicturePtr->height / sequence_control_set_ptr->picture_analysis_number_of_regions_per_height;

    // Loop over regions inside the picture
    for (regionInPictureWidthIndex = 0; regionInPictureWidthIndex < sequence_control_set_ptr->picture_analysis_number_of_regions_per_width; regionInPictureWidthIndex++) {  // loop over horizontal regions
        for (regionInPictureHeightIndex = 0; regionInPictureHeightIndex < sequence_control_set_ptr->picture_analysis_number_of_regions_per_height; regionInPictureHeightIndex++) { // loop over vertical regions


            // Initialize bins to 1
            InitializeBuffer_32bits_funcPtrArray[asm_type](picture_control_set_ptr->picture_histogram[regionInPictureWidthIndex][regionInPictureHeightIndex][0], 64, 0, 1);

            regionWidthOffset = (regionInPictureWidthIndex == sequence_control_set_ptr->picture_analysis_number_of_regions_per_width - 1) ?
                inputPicturePtr->width - (sequence_control_set_ptr->picture_analysis_number_of_regions_per_width * regionWidth) :
                0;

            regionHeightOffset = (regionInPictureHeightIndex == sequence_control_set_ptr->picture_analysis_number_of_regions_per_height - 1) ?
                inputPicturePtr->height - (sequence_control_set_ptr->picture_analysis_number_of_regions_per_height * regionHeight) :
                0;

            // Y Histogram
            CalculateHistogram(
                &inputPicturePtr->bufferY[(inputPicturePtr->origin_x + regionInPictureWidthIndex * regionWidth) + ((inputPicturePtr->origin_y + regionInPictureHeightIndex * regionHeight) * inputPicturePtr->strideY)],
                regionWidth + regionWidthOffset,
                regionHeight + regionHeightOffset,
                inputPicturePtr->strideY,
                1,
                picture_control_set_ptr->picture_histogram[regionInPictureWidthIndex][regionInPictureHeightIndex][0],
                &sum);

            picture_control_set_ptr->average_intensity_per_region[regionInPictureWidthIndex][regionInPictureHeightIndex][0] = (uint8_t)((sum + (((regionWidth + regionWidthOffset)*(regionHeight + regionHeightOffset)) >> 1)) / ((regionWidth + regionWidthOffset)*(regionHeight + regionHeightOffset)));
            (*sumAverageIntensityTotalRegionsLuma) += (sum << 4);
            for (histogramBin = 0; histogramBin < HISTOGRAM_NUMBER_OF_BINS; histogramBin++) { // Loop over the histogram bins
                picture_control_set_ptr->picture_histogram[regionInPictureWidthIndex][regionInPictureHeightIndex][0][histogramBin] =
                    picture_control_set_ptr->picture_histogram[regionInPictureWidthIndex][regionInPictureHeightIndex][0][histogramBin] << 4;
            }
        }
    }

    return;
}

void SubSampleChromaGeneratePixelIntensityHistogramBins(
    SequenceControlSet_t            *sequence_control_set_ptr,
    PictureParentControlSet_t       *picture_control_set_ptr,
    EbPictureBufferDesc_t           *inputPicturePtr,
    uint64_t                          *sumAverageIntensityTotalRegionsCb,
    uint64_t                          *sumAverageIntensityTotalRegionsCr,
    EbAsm                           asm_type) {

    uint64_t                          sum;
    uint32_t                          regionWidth;
    uint32_t                          regionHeight;
    uint32_t                          regionWidthOffset;
    uint32_t                          regionHeightOffset;
    uint32_t                          regionInPictureWidthIndex;
    uint32_t                          regionInPictureHeightIndex;

    uint16_t                          histogramBin;
    uint8_t                           decimStep = 4;

    regionWidth = inputPicturePtr->width / sequence_control_set_ptr->picture_analysis_number_of_regions_per_width;
    regionHeight = inputPicturePtr->height / sequence_control_set_ptr->picture_analysis_number_of_regions_per_height;

    // Loop over regions inside the picture
    for (regionInPictureWidthIndex = 0; regionInPictureWidthIndex < sequence_control_set_ptr->picture_analysis_number_of_regions_per_width; regionInPictureWidthIndex++) {  // loop over horizontal regions
        for (regionInPictureHeightIndex = 0; regionInPictureHeightIndex < sequence_control_set_ptr->picture_analysis_number_of_regions_per_height; regionInPictureHeightIndex++) { // loop over vertical regions


            // Initialize bins to 1
            InitializeBuffer_32bits_funcPtrArray[asm_type](picture_control_set_ptr->picture_histogram[regionInPictureWidthIndex][regionInPictureHeightIndex][1], 64, 0, 1);
            InitializeBuffer_32bits_funcPtrArray[asm_type](picture_control_set_ptr->picture_histogram[regionInPictureWidthIndex][regionInPictureHeightIndex][2], 64, 0, 1);

            regionWidthOffset = (regionInPictureWidthIndex == sequence_control_set_ptr->picture_analysis_number_of_regions_per_width - 1) ?
                inputPicturePtr->width - (sequence_control_set_ptr->picture_analysis_number_of_regions_per_width * regionWidth) :
                0;

            regionHeightOffset = (regionInPictureHeightIndex == sequence_control_set_ptr->picture_analysis_number_of_regions_per_height - 1) ?
                inputPicturePtr->height - (sequence_control_set_ptr->picture_analysis_number_of_regions_per_height * regionHeight) :
                0;


            // U Histogram
            CalculateHistogram(
                &inputPicturePtr->bufferCb[((inputPicturePtr->origin_x + regionInPictureWidthIndex * regionWidth) >> 1) + (((inputPicturePtr->origin_y + regionInPictureHeightIndex * regionHeight) >> 1) * inputPicturePtr->strideCb)],
                (regionWidth + regionWidthOffset) >> 1,
                (regionHeight + regionHeightOffset) >> 1,
                inputPicturePtr->strideCb,
                decimStep,
                picture_control_set_ptr->picture_histogram[regionInPictureWidthIndex][regionInPictureHeightIndex][1],
                &sum);

            sum = (sum << decimStep);
            *sumAverageIntensityTotalRegionsCb += sum;
            picture_control_set_ptr->average_intensity_per_region[regionInPictureWidthIndex][regionInPictureHeightIndex][1] = (uint8_t)((sum + (((regionWidth + regionWidthOffset)*(regionHeight + regionHeightOffset)) >> 3)) / (((regionWidth + regionWidthOffset)*(regionHeight + regionHeightOffset)) >> 2));

            for (histogramBin = 0; histogramBin < HISTOGRAM_NUMBER_OF_BINS; histogramBin++) { // Loop over the histogram bins
                picture_control_set_ptr->picture_histogram[regionInPictureWidthIndex][regionInPictureHeightIndex][1][histogramBin] =
                    picture_control_set_ptr->picture_histogram[regionInPictureWidthIndex][regionInPictureHeightIndex][1][histogramBin] << decimStep;
            }

            // V Histogram
            CalculateHistogram(
                &inputPicturePtr->bufferCr[((inputPicturePtr->origin_x + regionInPictureWidthIndex * regionWidth) >> 1) + (((inputPicturePtr->origin_y + regionInPictureHeightIndex * regionHeight) >> 1) * inputPicturePtr->strideCr)],
                (regionWidth + regionWidthOffset) >> 1,
                (regionHeight + regionHeightOffset) >> 1,
                inputPicturePtr->strideCr,
                decimStep,
                picture_control_set_ptr->picture_histogram[regionInPictureWidthIndex][regionInPictureHeightIndex][2],
                &sum);

            sum = (sum << decimStep);
            *sumAverageIntensityTotalRegionsCr += sum;
            picture_control_set_ptr->average_intensity_per_region[regionInPictureWidthIndex][regionInPictureHeightIndex][2] = (uint8_t)((sum + (((regionWidth + regionWidthOffset)*(regionHeight + regionHeightOffset)) >> 3)) / (((regionWidth + regionWidthOffset)*(regionHeight + regionHeightOffset)) >> 2));

            for (histogramBin = 0; histogramBin < HISTOGRAM_NUMBER_OF_BINS; histogramBin++) { // Loop over the histogram bins
                picture_control_set_ptr->picture_histogram[regionInPictureWidthIndex][regionInPictureHeightIndex][2][histogramBin] =
                    picture_control_set_ptr->picture_histogram[regionInPictureWidthIndex][regionInPictureHeightIndex][2][histogramBin] << decimStep;
            }
        }
    }
    return;

}

void EdgeDetectionMeanLumaChroma16x16(
    SequenceControlSet_t        *sequence_control_set_ptr,
    PictureParentControlSet_t   *picture_control_set_ptr,
    uint32_t                       totalLcuCount)
{

    uint32_t               sb_index;


    uint32_t maxGrad = 1;

    // The values are calculated for every 4th frame
    if ((picture_control_set_ptr->picture_number & 3) == 0) {
        for (sb_index = 0; sb_index < totalLcuCount; sb_index++) {

            SbStat_t *sb_stat_ptr = &picture_control_set_ptr->sb_stat_array[sb_index];

            EB_MEMSET(sb_stat_ptr, 0, sizeof(SbStat_t));
            SbParams_t     *sb_params = &sequence_control_set_ptr->sb_params_array[sb_index];
            if (sb_params->potential_logo_sb &&sb_params->is_complete_sb)

            {
                uint8_t *y_mean_ptr = picture_control_set_ptr->yMean[sb_index];
                uint8_t *cr_mean_ptr = picture_control_set_ptr->crMean[sb_index];
                uint8_t *cb_mean_ptr = picture_control_set_ptr->cbMean[sb_index];

                uint8_t rasterScanCuIndex;

                for (rasterScanCuIndex = RASTER_SCAN_CU_INDEX_16x16_0; rasterScanCuIndex <= RASTER_SCAN_CU_INDEX_16x16_15; rasterScanCuIndex++) {
                    uint8_t cu_index = rasterScanCuIndex - 5;
                    uint8_t x = cu_index & 3;
                    uint8_t y = (cu_index >> 2);
                    int32_t gradx = 0;
                    int32_t grady = 0;
                    int32_t nbcompx = 0;
                    int32_t nbcompy = 0;
                    if (x != 0)
                    {
                        gradx += ABS((int32_t)(y_mean_ptr[rasterScanCuIndex]) - (int32_t)(y_mean_ptr[rasterScanCuIndex - 1]));
                        gradx += ABS((int32_t)(cr_mean_ptr[rasterScanCuIndex]) - (int32_t)(cr_mean_ptr[rasterScanCuIndex - 1]));
                        gradx += ABS((int32_t)(cb_mean_ptr[rasterScanCuIndex]) - (int32_t)(cb_mean_ptr[rasterScanCuIndex - 1]));
                        nbcompx++;
                    }
                    if (x != 3)
                    {
                        gradx += ABS((int32_t)(y_mean_ptr[rasterScanCuIndex + 1]) - (int32_t)(y_mean_ptr[rasterScanCuIndex]));
                        gradx += ABS((int32_t)(cr_mean_ptr[rasterScanCuIndex + 1]) - (int32_t)(cr_mean_ptr[rasterScanCuIndex]));
                        gradx += ABS((int32_t)(cb_mean_ptr[rasterScanCuIndex + 1]) - (int32_t)(cb_mean_ptr[rasterScanCuIndex]));
                        nbcompx++;
                    }
                    gradx = gradx / nbcompx;


                    if (y != 0)
                    {
                        grady += ABS((int32_t)(y_mean_ptr[rasterScanCuIndex]) - (int32_t)(y_mean_ptr[rasterScanCuIndex - 4]));
                        grady += ABS((int32_t)(cr_mean_ptr[rasterScanCuIndex]) - (int32_t)(cr_mean_ptr[rasterScanCuIndex - 4]));
                        grady += ABS((int32_t)(cb_mean_ptr[rasterScanCuIndex]) - (int32_t)(cb_mean_ptr[rasterScanCuIndex - 4]));
                        nbcompy++;
                    }
                    if (y != 3)
                    {
                        grady += ABS((int32_t)(y_mean_ptr[rasterScanCuIndex + 4]) - (int32_t)(y_mean_ptr[rasterScanCuIndex]));
                        grady += ABS((int32_t)(cr_mean_ptr[rasterScanCuIndex + 4]) - (int32_t)(cr_mean_ptr[rasterScanCuIndex]));
                        grady += ABS((int32_t)(cb_mean_ptr[rasterScanCuIndex + 4]) - (int32_t)(cb_mean_ptr[rasterScanCuIndex]));

                        nbcompy++;
                    }

                    grady = grady / nbcompy;
                    sb_stat_ptr->cu_stat_array[rasterScanCuIndex].grad = (uint32_t)ABS(gradx) + ABS(grady);
                    if (sb_stat_ptr->cu_stat_array[rasterScanCuIndex].grad > maxGrad) {
                        maxGrad = sb_stat_ptr->cu_stat_array[rasterScanCuIndex].grad;
                    }
                }
            }
        }

        for (sb_index = 0; sb_index < totalLcuCount; sb_index++) {
            SbParams_t *sb_params = &sequence_control_set_ptr->sb_params_array[sb_index];
            if (sb_params->potential_logo_sb &&sb_params->is_complete_sb) {
                SbStat_t *sb_stat_ptr = &picture_control_set_ptr->sb_stat_array[sb_index];

                uint32_t rasterScanCuIndex;
                for (rasterScanCuIndex = RASTER_SCAN_CU_INDEX_16x16_0; rasterScanCuIndex <= RASTER_SCAN_CU_INDEX_16x16_15; rasterScanCuIndex++) {
                    sb_stat_ptr->cu_stat_array[rasterScanCuIndex].edge_cu = (uint16_t)MIN(((sb_stat_ptr->cu_stat_array[rasterScanCuIndex].grad * (255 * 3)) / maxGrad), 255) < 30 ? 0 : 1;
                }
            }
        }
    }
    else {
        for (sb_index = 0; sb_index < totalLcuCount; sb_index++) {

            SbStat_t *sb_stat_ptr = &picture_control_set_ptr->sb_stat_array[sb_index];

            EB_MEMSET(sb_stat_ptr, 0, sizeof(SbStat_t));
        }
    }
}

/******************************************************
* Edge map derivation
******************************************************/
void EdgeDetection(
    SequenceControlSet_t            *sequence_control_set_ptr,
    PictureParentControlSet_t       *picture_control_set_ptr)
{

    uint16_t  *variancePtr;
    uint32_t sb_total_count = picture_control_set_ptr->sb_total_count;
    uint64_t thrsldLevel0 = (picture_control_set_ptr->pic_avg_variance * 70) / 100;
    uint8_t  *meanPtr;
    uint32_t picture_width_in_sb = (sequence_control_set_ptr->luma_width + sequence_control_set_ptr->sb_sz - 1) / sequence_control_set_ptr->sb_sz;
    uint32_t picture_height_in_sb = (sequence_control_set_ptr->luma_height + sequence_control_set_ptr->sb_sz - 1) / sequence_control_set_ptr->sb_sz;
    uint32_t neighbourLcuIndex = 0;
    uint64_t similarityCount = 0;
    uint64_t similarityCount0 = 0;
    uint64_t similarityCount1 = 0;
    uint64_t similarityCount2 = 0;
    uint64_t similarityCount3 = 0;
    uint32_t sb_x = 0;
    uint32_t sb_y = 0;
    uint32_t sb_index;
    EbBool highVarianceLucFlag;

    uint32_t rasterScanCuIndex = 0;
    uint32_t numberOfEdgeLcu = 0;
    EbBool highIntensityLcuFlag;

    uint64_t neighbourLcuMean;
    int32_t i, j;

    uint8_t highIntensityTh = 180;
    uint8_t lowIntensityTh = 120;
    uint8_t highIntensityTh1 = 200;
    uint8_t veryLowIntensityTh = 20;

    for (sb_index = 0; sb_index < sb_total_count; ++sb_index) {

        sb_x = sb_index % picture_width_in_sb;
        sb_y = sb_index / picture_width_in_sb;

        EdgeLcuResults_t *edge_results_ptr = picture_control_set_ptr->edge_results_ptr;
        picture_control_set_ptr->edge_results_ptr[sb_index].edge_block_num = 0;
        picture_control_set_ptr->edge_results_ptr[sb_index].isolated_high_intensity_sb = 0;
        picture_control_set_ptr->sharp_edge_sb_flag[sb_index] = 0;

        if (sb_x > 0 && sb_x < (uint32_t)(picture_width_in_sb - 1) && sb_y >  0 && sb_y < (uint32_t)(picture_height_in_sb - 1)) {

            variancePtr = picture_control_set_ptr->variance[sb_index];
            meanPtr = picture_control_set_ptr->yMean[sb_index];


            similarityCount = 0;

            highVarianceLucFlag =
                (variancePtr[RASTER_SCAN_CU_INDEX_64x64] > thrsldLevel0) ? EB_TRUE : EB_FALSE;
            edge_results_ptr[sb_index].edge_block_num = highVarianceLucFlag;
            if (variancePtr[0] > highIntensityTh1) {
                uint8_t sharpEdge = 0;
                for (rasterScanCuIndex = RASTER_SCAN_CU_INDEX_16x16_0; rasterScanCuIndex <= RASTER_SCAN_CU_INDEX_16x16_15; rasterScanCuIndex++) {
                    sharpEdge = (variancePtr[rasterScanCuIndex] < veryLowIntensityTh) ? sharpEdge + 1 : sharpEdge;

                }
                if (sharpEdge > 4)
                {
                    picture_control_set_ptr->sharp_edge_sb_flag[sb_index] = 1;
                }
            }


            if (sb_x > 3 && sb_x < (uint32_t)(picture_width_in_sb - 4) && sb_y >  3 && sb_y < (uint32_t)(picture_height_in_sb - 4)) {

                highIntensityLcuFlag =
                    (meanPtr[RASTER_SCAN_CU_INDEX_64x64] > highIntensityTh) ? EB_TRUE : EB_FALSE;

                if (highIntensityLcuFlag) {

                    neighbourLcuIndex = sb_index - 1;
                    neighbourLcuMean = picture_control_set_ptr->yMean[neighbourLcuIndex][RASTER_SCAN_CU_INDEX_64x64];

                    similarityCount0 = (neighbourLcuMean < lowIntensityTh) ? 1 : 0;

                    neighbourLcuIndex = sb_index + 1;

                    neighbourLcuMean = picture_control_set_ptr->yMean[neighbourLcuIndex][RASTER_SCAN_CU_INDEX_64x64];
                    similarityCount1 = (neighbourLcuMean < lowIntensityTh) ? 1 : 0;

                    neighbourLcuIndex = sb_index - picture_width_in_sb;
                    neighbourLcuMean = picture_control_set_ptr->yMean[neighbourLcuIndex][RASTER_SCAN_CU_INDEX_64x64];
                    similarityCount2 = (neighbourLcuMean < lowIntensityTh) ? 1 : 0;

                    neighbourLcuIndex = sb_index + picture_width_in_sb;
                    neighbourLcuMean = picture_control_set_ptr->yMean[neighbourLcuIndex][RASTER_SCAN_CU_INDEX_64x64];
                    similarityCount3 = (neighbourLcuMean < lowIntensityTh) ? 1 : 0;

                    similarityCount = similarityCount0 + similarityCount1 + similarityCount2 + similarityCount3;

                    if (similarityCount > 0) {


                        for (i = -4; i < 5; i++) {
                            for (j = -4; j < 5; j++) {
                                neighbourLcuIndex = sb_index + (i * picture_width_in_sb) + j;
                                picture_control_set_ptr->edge_results_ptr[neighbourLcuIndex].isolated_high_intensity_sb = 1;
                            }
                        }
                    }
                }
            }


            if (highVarianceLucFlag) {
                numberOfEdgeLcu += edge_results_ptr[sb_index].edge_block_num;
            }
        }
    }
    return;
}
/******************************************************
* Calculate the variance of variance to determine Homogeneous regions. Note: Variance calculation should be on.
******************************************************/
void DetermineHomogeneousRegionInPicture(
    SequenceControlSet_t            *sequence_control_set_ptr,
    PictureParentControlSet_t       *picture_control_set_ptr)
{

    uint16_t  *variancePtr;
    uint32_t sb_index;

    uint32_t cuNum, cu_size, cuIndexOffset, cuH, cuW;
    uint64_t nullVarCnt = 0;
    uint64_t veryLowVarCnt = 0;
    uint64_t varLcuCnt = 0;
    uint32_t sb_total_count = picture_control_set_ptr->sb_total_count;

    for (sb_index = 0; sb_index < sb_total_count; ++sb_index) {
        uint64_t meanSqrVariance32x32Based[4] = { 0 }, meanVariance32x32Based[4] = { 0 };

        uint64_t meanSqrVariance64x64Based = 0, meanVariance64x64Based = 0;
        uint64_t varOfVar64x64Based = 0;

        SbParams_t sb_params = sequence_control_set_ptr->sb_params_array[sb_index];

        // Initialize
        picture_control_set_ptr->sb_homogeneous_area_array[sb_index] = EB_TRUE;

        variancePtr = picture_control_set_ptr->variance[sb_index];

        if (sb_params.is_complete_sb) {

            nullVarCnt += (variancePtr[ME_TIER_ZERO_PU_64x64] == 0) ? 1 : 0;

            varLcuCnt++;

            veryLowVarCnt += ((variancePtr[ME_TIER_ZERO_PU_64x64]) < LCU_LOW_VAR_TH) ? 1 : 0;
            cu_size = 8;
            cuIndexOffset = ME_TIER_ZERO_PU_8x8_0;
            cuNum = 64 / cu_size;

            //Variance of 8x8 blocks in a 32x32
            for (cuH = 0; cuH < (cuNum / 2); cuH++) {
                for (cuW = 0; cuW < (cuNum / 2); cuW++) {

                    meanSqrVariance32x32Based[0] += (variancePtr[cuIndexOffset + cuH * cuNum + cuW])*(variancePtr[cuIndexOffset + cuH * cuNum + cuW]);
                    meanVariance32x32Based[0] += (variancePtr[cuIndexOffset + cuH * cuNum + cuW]);

                    meanSqrVariance32x32Based[1] += (variancePtr[cuIndexOffset + cuH * cuNum + cuW + 4])*(variancePtr[cuIndexOffset + cuH * cuNum + cuW + 4]);
                    meanVariance32x32Based[1] += (variancePtr[cuIndexOffset + cuH * cuNum + cuW + 4]);

                    meanSqrVariance32x32Based[2] += (variancePtr[cuIndexOffset + (cuH + 4)*cuNum + cuW])*(variancePtr[cuIndexOffset + (cuH + 4)*cuNum + cuW]);
                    meanVariance32x32Based[2] += (variancePtr[cuIndexOffset + (cuH + 4)*cuNum + cuW]);

                    meanSqrVariance32x32Based[3] += (variancePtr[cuIndexOffset + (cuH + 4)*cuNum + cuW + 4])*(variancePtr[cuIndexOffset + (cuH + 4)*cuNum + cuW + 4]);
                    meanVariance32x32Based[3] += (variancePtr[cuIndexOffset + (cuH + 4)*cuNum + cuW + 4]);

                }
            }

            meanSqrVariance32x32Based[0] = meanSqrVariance32x32Based[0] >> 4;
            meanVariance32x32Based[0] = meanVariance32x32Based[0] >> 4;
            picture_control_set_ptr->var_of_var32x32_based_sb_array[sb_index][0] = meanSqrVariance32x32Based[0] - meanVariance32x32Based[0] * meanVariance32x32Based[0];

            meanSqrVariance32x32Based[1] = meanSqrVariance32x32Based[1] >> 4;
            meanVariance32x32Based[1] = meanVariance32x32Based[1] >> 4;
            picture_control_set_ptr->var_of_var32x32_based_sb_array[sb_index][1] = meanSqrVariance32x32Based[1] - meanVariance32x32Based[1] * meanVariance32x32Based[1];

            meanSqrVariance32x32Based[2] = meanSqrVariance32x32Based[2] >> 4;
            meanVariance32x32Based[2] = meanVariance32x32Based[2] >> 4;
            picture_control_set_ptr->var_of_var32x32_based_sb_array[sb_index][2] = meanSqrVariance32x32Based[2] - meanVariance32x32Based[2] * meanVariance32x32Based[2];

            meanSqrVariance32x32Based[3] = meanSqrVariance32x32Based[3] >> 4;
            meanVariance32x32Based[3] = meanVariance32x32Based[3] >> 4;
            picture_control_set_ptr->var_of_var32x32_based_sb_array[sb_index][3] = meanSqrVariance32x32Based[3] - meanVariance32x32Based[3] * meanVariance32x32Based[3];

            // Compute the 64x64 based variance of variance
            {
                uint32_t varIndex;
                // Loop over all 8x8s in a 64x64
                for (varIndex = ME_TIER_ZERO_PU_8x8_0; varIndex <= ME_TIER_ZERO_PU_8x8_63; varIndex++) {
                    meanSqrVariance64x64Based += variancePtr[varIndex] * variancePtr[varIndex];
                    meanVariance64x64Based += variancePtr[varIndex];
                }

                meanSqrVariance64x64Based = meanSqrVariance64x64Based >> 6;
                meanVariance64x64Based = meanVariance64x64Based >> 6;

                // Compute variance
                varOfVar64x64Based = meanSqrVariance64x64Based - meanVariance64x64Based * meanVariance64x64Based;

                // Turn off detail preservation if the varOfVar is greater than a threshold
                if (varOfVar64x64Based > VAR_BASED_DETAIL_PRESERVATION_SELECTOR_THRSLHD)
                {
                    picture_control_set_ptr->sb_homogeneous_area_array[sb_index] = EB_FALSE;
                }
            }

        }
        else {

            // Should be re-calculated and scaled properly
            picture_control_set_ptr->var_of_var32x32_based_sb_array[sb_index][0] = 0xFFFFFFFFFFFFFFFF;
            picture_control_set_ptr->var_of_var32x32_based_sb_array[sb_index][1] = 0xFFFFFFFFFFFFFFFF;
            picture_control_set_ptr->var_of_var32x32_based_sb_array[sb_index][2] = 0xFFFFFFFFFFFFFFFF;
            picture_control_set_ptr->var_of_var32x32_based_sb_array[sb_index][3] = 0xFFFFFFFFFFFFFFFF;
        }
    }
    picture_control_set_ptr->very_low_var_pic_flag = EB_FALSE;
    if (((veryLowVarCnt * 100) / varLcuCnt) > PIC_LOW_VAR_PERCENTAGE_TH) {
        picture_control_set_ptr->very_low_var_pic_flag = EB_TRUE;
    }

    picture_control_set_ptr->logo_pic_flag = EB_FALSE;
    if (((veryLowVarCnt * 100) / varLcuCnt) > 80) {
        picture_control_set_ptr->logo_pic_flag = EB_TRUE;
    }

    return;
}
/************************************************
 * ComputePictureSpatialStatistics
 ** Compute Block Variance
 ** Compute Picture Variance
 ** Compute Block Mean for all blocks in the picture
 ************************************************/
void ComputePictureSpatialStatistics(
    SequenceControlSet_t            *sequence_control_set_ptr,
    PictureParentControlSet_t       *picture_control_set_ptr,
    EbPictureBufferDesc_t           *inputPicturePtr,
    EbPictureBufferDesc_t           *inputPaddedPicturePtr,
    uint32_t                           sb_total_count,
    EbAsm                           asm_type)
{
    uint32_t sb_index;
    uint32_t sb_origin_x;        // to avoid using child PCS
    uint32_t sb_origin_y;
    uint32_t inputLumaOriginIndex;
    uint32_t inputCbOriginIndex;
    uint32_t inputCrOriginIndex;
    uint64_t picTotVariance;

    // Variance
    picTotVariance = 0;

    for (sb_index = 0; sb_index < picture_control_set_ptr->sb_total_count; ++sb_index) {
        SbParams_t   *sb_params = &sequence_control_set_ptr->sb_params_array[sb_index];

        sb_origin_x = sb_params->origin_x;
        sb_origin_y = sb_params->origin_y;
        inputLumaOriginIndex = (inputPaddedPicturePtr->origin_y + sb_origin_y) * inputPaddedPicturePtr->strideY +
            inputPaddedPicturePtr->origin_x + sb_origin_x;

        inputCbOriginIndex = ((inputPicturePtr->origin_y + sb_origin_y) >> 1) * inputPicturePtr->strideCb + ((inputPicturePtr->origin_x + sb_origin_x) >> 1);
        inputCrOriginIndex = ((inputPicturePtr->origin_y + sb_origin_y) >> 1) * inputPicturePtr->strideCr + ((inputPicturePtr->origin_x + sb_origin_x) >> 1);

        ComputeBlockMeanComputeVariance(
            sequence_control_set_ptr,
            picture_control_set_ptr,
            inputPaddedPicturePtr,
            sb_index,
            inputLumaOriginIndex,
            asm_type);

        if (sb_params->is_complete_sb) {

            ComputeChromaBlockMean(
                sequence_control_set_ptr,
                picture_control_set_ptr,
                inputPicturePtr,
                sb_index,
                inputCbOriginIndex,
                inputCrOriginIndex,
                asm_type);
        }
        else {
            ZeroOutChromaBlockMean(
                picture_control_set_ptr,
                sb_index);
        }

        picTotVariance += (picture_control_set_ptr->variance[sb_index][RASTER_SCAN_CU_INDEX_64x64]);
    }

    picture_control_set_ptr->pic_avg_variance = (uint16_t)(picTotVariance / sb_total_count);
    // Calculate the variance of variance to determine Homogeneous regions. Note: Variance calculation should be on.
    DetermineHomogeneousRegionInPicture(
        sequence_control_set_ptr,
        picture_control_set_ptr);

    EdgeDetectionMeanLumaChroma16x16(
        sequence_control_set_ptr,
        picture_control_set_ptr,
        sequence_control_set_ptr->sb_total_count);

    EdgeDetection(
        sequence_control_set_ptr,
        picture_control_set_ptr);


    return;
}

void CalculateInputAverageIntensity(
    SequenceControlSet_t            *sequence_control_set_ptr,
    PictureParentControlSet_t       *picture_control_set_ptr,
    EbPictureBufferDesc_t           *inputPicturePtr,
    uint64_t                           sumAverageIntensityTotalRegionsLuma,
    uint64_t                           sumAverageIntensityTotalRegionsCb,
    uint64_t                           sumAverageIntensityTotalRegionsCr,
    EbAsm                           asm_type)
{

    if (sequence_control_set_ptr->scd_mode == SCD_MODE_0) {
        uint16_t blockIndexInWidth;
        uint16_t blockIndexInHeight;
        uint64_t mean = 0;

        (void)(asm_type);
        const uint16_t strideY = inputPicturePtr->strideY;
        // Loop over 8x8 blocks and calculates the mean value
        if (sequence_control_set_ptr->block_mean_calc_prec == BLOCK_MEAN_PREC_FULL) {
            for (blockIndexInHeight = 0; blockIndexInHeight < inputPicturePtr->height >> 3; ++blockIndexInHeight) {
                for (blockIndexInWidth = 0; blockIndexInWidth < inputPicturePtr->width >> 3; ++blockIndexInWidth) {
                    mean += ComputeMeanFunc[0][asm_type](&(inputPicturePtr->bufferY[(blockIndexInWidth << 3) + (blockIndexInHeight << 3) * inputPicturePtr->strideY]), inputPicturePtr->strideY, 8, 8);
                }
            }
        }
        else {
            for (blockIndexInHeight = 0; blockIndexInHeight < inputPicturePtr->height >> 3; ++blockIndexInHeight) {
                for (blockIndexInWidth = 0; blockIndexInWidth < inputPicturePtr->width >> 3; ++blockIndexInWidth) {
                    mean += ComputeSubMean8x8_SSE2_INTRIN(&(inputPicturePtr->bufferY[(blockIndexInWidth << 3) + (blockIndexInHeight << 3) * strideY]), strideY);
                }
            }
        }
        mean = ((mean + ((inputPicturePtr->height* inputPicturePtr->width) >> 7)) / ((inputPicturePtr->height* inputPicturePtr->width) >> 6));
        mean = (mean + (1 << (MEAN_PRECISION - 1))) >> MEAN_PRECISION;
        picture_control_set_ptr->average_intensity[0] = (uint8_t)mean;
    }

    else {
        picture_control_set_ptr->average_intensity[0] = (uint8_t)((sumAverageIntensityTotalRegionsLuma + ((inputPicturePtr->width*inputPicturePtr->height) >> 1)) / (inputPicturePtr->width*inputPicturePtr->height));
        picture_control_set_ptr->average_intensity[1] = (uint8_t)((sumAverageIntensityTotalRegionsCb + ((inputPicturePtr->width*inputPicturePtr->height) >> 3)) / ((inputPicturePtr->width*inputPicturePtr->height) >> 2));
        picture_control_set_ptr->average_intensity[2] = (uint8_t)((sumAverageIntensityTotalRegionsCr + ((inputPicturePtr->width*inputPicturePtr->height) >> 3)) / ((inputPicturePtr->width*inputPicturePtr->height) >> 2));
    }

    return;
}

/************************************************
 * Gathering statistics per picture
 ** Calculating the pixel intensity histogram bins per picture needed for SCD
 ** Computing Picture Variance
 ************************************************/
void GatheringPictureStatistics(
    SequenceControlSet_t            *sequence_control_set_ptr,
    PictureParentControlSet_t       *picture_control_set_ptr,
    EbPictureBufferDesc_t           *inputPicturePtr,
    EbPictureBufferDesc_t           *inputPaddedPicturePtr,
    EbPictureBufferDesc_t            *sixteenthDecimatedPicturePtr,
    uint32_t                           sb_total_count,
    EbAsm                           asm_type)
{

    uint64_t                          sumAverageIntensityTotalRegionsLuma = 0;
    uint64_t                          sumAverageIntensityTotalRegionsCb = 0;
    uint64_t                          sumAverageIntensityTotalRegionsCr = 0;

    // Histogram bins
        // Use 1/16 Luma for Histogram generation
        // 1/16 input ready
    SubSampleLumaGeneratePixelIntensityHistogramBins(
        sequence_control_set_ptr,
        picture_control_set_ptr,
        sixteenthDecimatedPicturePtr,
        &sumAverageIntensityTotalRegionsLuma,
        asm_type);

    // Use 1/4 Chroma for Histogram generation
    // 1/4 input not ready => perform operation on the fly
    SubSampleChromaGeneratePixelIntensityHistogramBins(
        sequence_control_set_ptr,
        picture_control_set_ptr,
        inputPicturePtr,
        &sumAverageIntensityTotalRegionsCb,
        &sumAverageIntensityTotalRegionsCr,
        asm_type);
    //
    // Calculate the LUMA average intensity
    CalculateInputAverageIntensity(
        sequence_control_set_ptr,
        picture_control_set_ptr,
        inputPicturePtr,
        sumAverageIntensityTotalRegionsLuma,
        sumAverageIntensityTotalRegionsCb,
        sumAverageIntensityTotalRegionsCr,
        asm_type);


    ComputePictureSpatialStatistics(
        sequence_control_set_ptr,
        picture_control_set_ptr,
        inputPicturePtr,
        inputPaddedPicturePtr,
        sb_total_count,
        asm_type);

    return;
}
/************************************************
 * Pad Picture at the right and bottom sides
 ** To match a multiple of min CU size in width and height
 ************************************************/
void PadPictureToMultipleOfMinCuSizeDimensions(
    SequenceControlSet_t            *sequence_control_set_ptr,
    EbPictureBufferDesc_t           *inputPicturePtr)
{
    EbBool                          is16BitInput = (EbBool)(sequence_control_set_ptr->static_config.encoder_bit_depth > EB_8BIT);
    // Input Picture Padding
    pad_input_picture(
        &inputPicturePtr->bufferY[inputPicturePtr->origin_x + (inputPicturePtr->origin_y * inputPicturePtr->strideY)],
        inputPicturePtr->strideY,
        (inputPicturePtr->width - sequence_control_set_ptr->pad_right),
        (inputPicturePtr->height - sequence_control_set_ptr->pad_bottom),
        sequence_control_set_ptr->pad_right,
        sequence_control_set_ptr->pad_bottom);

    pad_input_picture(
        &inputPicturePtr->bufferCb[(inputPicturePtr->origin_x >> 1) + ((inputPicturePtr->origin_y >> 1) * inputPicturePtr->strideCb)],
        inputPicturePtr->strideCb,
        (inputPicturePtr->width - sequence_control_set_ptr->pad_right) >> 1,
        (inputPicturePtr->height - sequence_control_set_ptr->pad_bottom) >> 1,
        sequence_control_set_ptr->pad_right >> 1,
        sequence_control_set_ptr->pad_bottom >> 1);

    pad_input_picture(
        &inputPicturePtr->bufferCr[(inputPicturePtr->origin_x >> 1) + ((inputPicturePtr->origin_y >> 1) * inputPicturePtr->strideCr)],
        inputPicturePtr->strideCr,
        (inputPicturePtr->width - sequence_control_set_ptr->pad_right) >> 1,
        (inputPicturePtr->height - sequence_control_set_ptr->pad_bottom) >> 1,
        sequence_control_set_ptr->pad_right >> 1,
        sequence_control_set_ptr->pad_bottom >> 1);

    if (is16BitInput)
    {
        pad_input_picture(
            &inputPicturePtr->bufferBitIncY[inputPicturePtr->origin_x + (inputPicturePtr->origin_y * inputPicturePtr->strideBitIncY)],
            inputPicturePtr->strideBitIncY,
            (inputPicturePtr->width - sequence_control_set_ptr->pad_right),
            (inputPicturePtr->height - sequence_control_set_ptr->pad_bottom),
            sequence_control_set_ptr->pad_right,
            sequence_control_set_ptr->pad_bottom);

        pad_input_picture(
            &inputPicturePtr->bufferBitIncCb[(inputPicturePtr->origin_x >> 1) + ((inputPicturePtr->origin_y >> 1) * inputPicturePtr->strideBitIncCb)],
            inputPicturePtr->strideBitIncCb,
            (inputPicturePtr->width - sequence_control_set_ptr->pad_right) >> 1,
            (inputPicturePtr->height - sequence_control_set_ptr->pad_bottom) >> 1,
            sequence_control_set_ptr->pad_right >> 1,
            sequence_control_set_ptr->pad_bottom >> 1);

        pad_input_picture(
            &inputPicturePtr->bufferBitIncCr[(inputPicturePtr->origin_x >> 1) + ((inputPicturePtr->origin_y >> 1) * inputPicturePtr->strideBitIncCr)],
            inputPicturePtr->strideBitIncCr,
            (inputPicturePtr->width - sequence_control_set_ptr->pad_right) >> 1,
            (inputPicturePtr->height - sequence_control_set_ptr->pad_bottom) >> 1,
            sequence_control_set_ptr->pad_right >> 1,
            sequence_control_set_ptr->pad_bottom >> 1);

    }

    return;
}

/************************************************
 * Pad Picture at the right and bottom sides
 ** To complete border SB smaller than SB size
 ************************************************/
void PadPictureToMultipleOfLcuDimensions(
    EbPictureBufferDesc_t           *inputPaddedPicturePtr
)
{

    // Generate Padding
    generate_padding(
        &inputPaddedPicturePtr->bufferY[0],
        inputPaddedPicturePtr->strideY,
        inputPaddedPicturePtr->width,
        inputPaddedPicturePtr->height,
        inputPaddedPicturePtr->origin_x,
        inputPaddedPicturePtr->origin_y);

    return;
}

/************************************************
* 1/4 & 1/16 input picture decimation
************************************************/
void DecimateInputPicture(
    PictureParentControlSet_t       *picture_control_set_ptr,
    EbPictureBufferDesc_t           *inputPaddedPicturePtr,
    EbPictureBufferDesc_t           *quarterDecimatedPicturePtr,
    EbPictureBufferDesc_t           *sixteenthDecimatedPicturePtr) {


    // Decimate input picture for HME L0 and L1
    if (picture_control_set_ptr->enable_hme_flag) {

        if (picture_control_set_ptr->enable_hme_level1_flag) {
            Decimation2D(
                &inputPaddedPicturePtr->bufferY[inputPaddedPicturePtr->origin_x + inputPaddedPicturePtr->origin_y * inputPaddedPicturePtr->strideY],
                inputPaddedPicturePtr->strideY,
                inputPaddedPicturePtr->width,
                inputPaddedPicturePtr->height,
                &quarterDecimatedPicturePtr->bufferY[quarterDecimatedPicturePtr->origin_x + quarterDecimatedPicturePtr->origin_x*quarterDecimatedPicturePtr->strideY],
                quarterDecimatedPicturePtr->strideY,
                2);
            generate_padding(
                &quarterDecimatedPicturePtr->bufferY[0],
                quarterDecimatedPicturePtr->strideY,
                quarterDecimatedPicturePtr->width,
                quarterDecimatedPicturePtr->height,
                quarterDecimatedPicturePtr->origin_x,
                quarterDecimatedPicturePtr->origin_y);

        }

        if (picture_control_set_ptr->enable_hme_level0_flag) {

            // Sixteenth Input Picture Decimation
            Decimation2D(
                &inputPaddedPicturePtr->bufferY[inputPaddedPicturePtr->origin_x + inputPaddedPicturePtr->origin_y * inputPaddedPicturePtr->strideY],
                inputPaddedPicturePtr->strideY,
                inputPaddedPicturePtr->width,
                inputPaddedPicturePtr->height,
                &sixteenthDecimatedPicturePtr->bufferY[sixteenthDecimatedPicturePtr->origin_x + sixteenthDecimatedPicturePtr->origin_x*sixteenthDecimatedPicturePtr->strideY],
                sixteenthDecimatedPicturePtr->strideY,
                4);

            generate_padding(
                &sixteenthDecimatedPicturePtr->bufferY[0],
                sixteenthDecimatedPicturePtr->strideY,
                sixteenthDecimatedPicturePtr->width,
                sixteenthDecimatedPicturePtr->height,
                sixteenthDecimatedPicturePtr->origin_x,
                sixteenthDecimatedPicturePtr->origin_y);

        }
    }
}

/************************************************
 * Picture Analysis Kernel
 * The Picture Analysis Process pads & decimates the input pictures.
 * The Picture Analysis also includes creating an n-bin Histogram,
 * gathering picture 1st and 2nd moment statistics for each 8x8 block,
 * which are used to compute variance.
 * The Picture Analysis process is multithreaded, so pictures can be
 * processed out of order as long as all inputs are available.
 ************************************************/
void* PictureAnalysisKernel(void *input_ptr)
{
    PictureAnalysisContext_t        *context_ptr = (PictureAnalysisContext_t*)input_ptr;
    PictureParentControlSet_t       *picture_control_set_ptr;
    SequenceControlSet_t            *sequence_control_set_ptr;

    EbObjectWrapper_t               *inputResultsWrapperPtr;
    ResourceCoordinationResults_t   *inputResultsPtr;
    EbObjectWrapper_t               *outputResultsWrapperPtr;
    PictureAnalysisResults_t        *outputResultsPtr;
    EbPaReferenceObject_t           *paReferenceObject;

    EbPictureBufferDesc_t           *inputPaddedPicturePtr;
    EbPictureBufferDesc_t           *quarterDecimatedPicturePtr;
    EbPictureBufferDesc_t           *sixteenthDecimatedPicturePtr;
    EbPictureBufferDesc_t           *inputPicturePtr;

    // Variance
    uint32_t                          picture_width_in_sb;
    uint32_t                          pictureHeighInLcu;
    uint32_t                          sb_total_count;
    EbAsm                          asm_type;

    for (;;) {

        // Get Input Full Object
        EbGetFullObject(
            context_ptr->resourceCoordinationResultsInputFifoPtr,
            &inputResultsWrapperPtr);

        inputResultsPtr = (ResourceCoordinationResults_t*)inputResultsWrapperPtr->objectPtr;
        picture_control_set_ptr = (PictureParentControlSet_t*)inputResultsPtr->pictureControlSetWrapperPtr->objectPtr;
        sequence_control_set_ptr = (SequenceControlSet_t*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->objectPtr;
        inputPicturePtr = picture_control_set_ptr->enhanced_picture_ptr;

        paReferenceObject = (EbPaReferenceObject_t*)picture_control_set_ptr->pa_reference_picture_wrapper_ptr->objectPtr;
        inputPaddedPicturePtr = (EbPictureBufferDesc_t*)paReferenceObject->inputPaddedPicturePtr;
        quarterDecimatedPicturePtr = (EbPictureBufferDesc_t*)paReferenceObject->quarterDecimatedPicturePtr;
        sixteenthDecimatedPicturePtr = (EbPictureBufferDesc_t*)paReferenceObject->sixteenthDecimatedPicturePtr;

        // Variance
        picture_width_in_sb = (sequence_control_set_ptr->luma_width + sequence_control_set_ptr->sb_sz - 1) / sequence_control_set_ptr->sb_sz;
        pictureHeighInLcu = (sequence_control_set_ptr->luma_height + sequence_control_set_ptr->sb_sz - 1) / sequence_control_set_ptr->sb_sz;
        sb_total_count = picture_width_in_sb * pictureHeighInLcu;

        asm_type = sequence_control_set_ptr->encode_context_ptr->asm_type;

        // Set picture parameters to account for subpicture, picture scantype, and set regions by resolutions
        SetPictureParametersForStatisticsGathering(
            sequence_control_set_ptr);

        // Pad pictures to multiple min cu size
        PadPictureToMultipleOfMinCuSizeDimensions(
            sequence_control_set_ptr,
            inputPicturePtr);

        // Pre processing operations performed on the input picture
        PicturePreProcessingOperations(
            picture_control_set_ptr,
            inputPicturePtr,
            context_ptr,
            sequence_control_set_ptr,
            quarterDecimatedPicturePtr,
            sixteenthDecimatedPicturePtr,
            sb_total_count,
            picture_width_in_sb,
            asm_type);

        // Pad input picture to complete border LCUs
        PadPictureToMultipleOfLcuDimensions(
            inputPaddedPicturePtr);

        // 1/4 & 1/16 input picture decimation
        DecimateInputPicture(
            picture_control_set_ptr,
            inputPaddedPicturePtr,
            quarterDecimatedPicturePtr,
            sixteenthDecimatedPicturePtr);

        // Gathering statistics of input picture, including Variance Calculation, Histogram Bins
        GatheringPictureStatistics(
            sequence_control_set_ptr,
            picture_control_set_ptr,
            inputPicturePtr,
            inputPaddedPicturePtr,
            sixteenthDecimatedPicturePtr,
            sb_total_count,
            asm_type);


        // Hold the 64x64 variance and mean in the reference frame
        uint32_t sb_index;
        for (sb_index = 0; sb_index < picture_control_set_ptr->sb_total_count; ++sb_index) {
            paReferenceObject->variance[sb_index] = picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_64x64];
            paReferenceObject->yMean[sb_index] = picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_64x64];

        }

        // Get Empty Results Object
        EbGetEmptyObject(
            context_ptr->pictureAnalysisResultsOutputFifoPtr,
            &outputResultsWrapperPtr);

        outputResultsPtr = (PictureAnalysisResults_t*)outputResultsWrapperPtr->objectPtr;
        outputResultsPtr->pictureControlSetWrapperPtr = inputResultsPtr->pictureControlSetWrapperPtr;

        // Release the Input Results
        EbReleaseObject(inputResultsWrapperPtr);

        // Post the Full Results Object
        EbPostFullObject(outputResultsWrapperPtr);

    }
    return EB_NULL;
}
