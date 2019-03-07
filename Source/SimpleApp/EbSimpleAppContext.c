/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/
#include <stdlib.h>
#include "EbSimpleAppContext.h"
#define INPUT_SIZE_576p_TH                0x90000        // 0.58 Million
#define INPUT_SIZE_1080i_TH                0xB71B0        // 0.75 Million
#define INPUT_SIZE_1080p_TH                0x1AB3F0    // 1.75 Million
#define INPUT_SIZE_4K_TH                0x29F630    // 2.75 Million
#define EB_OUTPUTSTREAMBUFFERSIZE_MACRO(ResolutionSize)                ((ResolutionSize) < (INPUT_SIZE_1080i_TH) ? 0x1E8480 : (ResolutionSize) < (INPUT_SIZE_1080p_TH) ? 0x2DC6C0 : (ResolutionSize) < (INPUT_SIZE_4K_TH) ? 0x2DC6C0 : 0x2DC6C0  )
EbErrorType AllocateFrameBuffer(
    EbConfig_t        *config,
    uint8_t           *p_buffer)
{
    EbErrorType   return_error = EB_ErrorNone;
    const int32_t tenBitPackedMode = (config->encoderBitDepth > 8) && (config->compressedTenBitFormat == 0) ? 1 : 0;

    // Determine size of each plane
    const size_t luma8bitSize =
        config->inputPaddedWidth    *
        config->inputPaddedHeight   *
        (1 << tenBitPackedMode);

    const size_t chroma8bitSize = luma8bitSize >> 2;
    const size_t luma10bitSize = (config->encoderBitDepth > 8 && tenBitPackedMode == 0) ? luma8bitSize : 0;
    const size_t chroma10bitSize = (config->encoderBitDepth > 8 && tenBitPackedMode == 0) ? chroma8bitSize : 0;

    // Determine
    EbSvtIOFormat* inputPtr = (EbSvtIOFormat*)p_buffer;

    if (luma8bitSize) {
        inputPtr->luma = (uint8_t*)malloc(luma8bitSize);
        if (!inputPtr->luma) return EB_ErrorInsufficientResources;
    }
    else {
        inputPtr->luma = 0;
    }
    if (chroma8bitSize) {
        inputPtr->cb = (uint8_t*)malloc(chroma8bitSize);
        if (!inputPtr->cb) return EB_ErrorInsufficientResources;
    }
    else {
        inputPtr->cb = 0;
    }
    if (chroma8bitSize) {
        inputPtr->cr = (uint8_t*)malloc(chroma8bitSize);
        if (!inputPtr->cr) return EB_ErrorInsufficientResources;
    }
    else {
        inputPtr->cr = 0;
    }
    if (luma10bitSize) {
        inputPtr->lumaExt = (uint8_t*)malloc(luma10bitSize);
        if (!inputPtr->lumaExt) return EB_ErrorInsufficientResources;
    }
    else {
        inputPtr->lumaExt = 0;
    }
    if (chroma10bitSize) {
        inputPtr->cbExt = (uint8_t*)malloc(chroma10bitSize);
        if (!inputPtr->cbExt) return EB_ErrorInsufficientResources;
    }
    else {
        inputPtr->cbExt = 0;
    }
    if (chroma10bitSize) {
        inputPtr->crExt = (uint8_t*)malloc(chroma10bitSize);
        if (!inputPtr->crExt) return EB_ErrorInsufficientResources;
    }
    else {
        inputPtr->crExt = 0;
    }
    return return_error;
}
/***********************************
 * AppContext Constructor
 ***********************************/
EbErrorType EbAppContextCtor(
    EbAppContext_t *contextPtr,
    EbConfig_t     *config)
{
    EbErrorType   return_error = EB_ErrorInsufficientResources;

    // Input Buffer
    contextPtr->inputPictureBuffer = (EbBufferHeaderType*)malloc(sizeof(EbBufferHeaderType));
    if (!contextPtr->inputPictureBuffer) return return_error;

    contextPtr->inputPictureBuffer->p_buffer = (uint8_t*)malloc(sizeof(EbSvtIOFormat));
    if (!contextPtr->inputPictureBuffer->p_buffer) return return_error;

    contextPtr->inputPictureBuffer->size = sizeof(EbBufferHeaderType);
    contextPtr->inputPictureBuffer->p_app_private = NULL;
    contextPtr->inputPictureBuffer->pic_type = EB_AV1_INVALID_PICTURE;
    // Allocate frame buffer for the p_buffer
    AllocateFrameBuffer(
        config,
        contextPtr->inputPictureBuffer->p_buffer);

    // output buffer
    contextPtr->outputStreamBuffer = (EbBufferHeaderType*)malloc(sizeof(EbBufferHeaderType));
    if (!contextPtr->outputStreamBuffer) return return_error;

    contextPtr->outputStreamBuffer->p_buffer = (uint8_t*)malloc(EB_OUTPUTSTREAMBUFFERSIZE_MACRO(config->sourceWidth*config->sourceHeight));
    if (!contextPtr->outputStreamBuffer->p_buffer) return return_error;

    contextPtr->outputStreamBuffer->size = sizeof(EbBufferHeaderType);
    contextPtr->outputStreamBuffer->n_alloc_len = EB_OUTPUTSTREAMBUFFERSIZE_MACRO(config->sourceWidth*config->sourceHeight);
    contextPtr->outputStreamBuffer->p_app_private = NULL;
    contextPtr->outputStreamBuffer->pic_type = EB_AV1_INVALID_PICTURE;

    // recon buffer
    if (config->reconFile) {
        contextPtr->recon_buffer = (EbBufferHeaderType*)malloc(sizeof(EbBufferHeaderType));
        if (!contextPtr->recon_buffer) return return_error;
        const size_t lumaSize =
            config->inputPaddedWidth    *
            config->inputPaddedHeight;
        // both u and v
        const size_t chromaSize = lumaSize >> 1;
        const size_t tenBit = (config->encoderBitDepth > 8);
        const size_t frameSize = (lumaSize + chromaSize) << tenBit;

        // Initialize Header
        contextPtr->recon_buffer->size = sizeof(EbBufferHeaderType);

        contextPtr->recon_buffer->p_buffer = (uint8_t*)malloc(frameSize);
        if (!contextPtr->recon_buffer->p_buffer) return return_error;

        contextPtr->recon_buffer->n_alloc_len = (uint32_t)frameSize;
        contextPtr->recon_buffer->p_app_private = NULL;
    }
    else
        contextPtr->recon_buffer = NULL;
    return EB_ErrorNone;
}

/***********************************
 * AppContext Destructor
 ***********************************/
void EbAppContextDtor(
    EbAppContext_t *contextPtr)
{
    EbSvtIOFormat *inputPtr = (EbSvtIOFormat*)contextPtr->inputPictureBuffer->p_buffer;
    free(inputPtr->luma);
    free(inputPtr->cb);
    free(inputPtr->cr);
    free(inputPtr->lumaExt);
    free(inputPtr->cbExt);
    free(inputPtr->crExt);
    free(contextPtr->inputPictureBuffer->p_buffer);
    free(contextPtr->outputStreamBuffer->p_buffer);
    free(contextPtr->inputPictureBuffer);
    free(contextPtr->outputStreamBuffer);
    if(contextPtr->recon_buffer)
        free(contextPtr->recon_buffer);
}

/***********************************************
* Copy configuration parameters from
*  The config structure, to the
*  callback structure to send to the library
***********************************************/
EbErrorType CopyConfigurationParameters(
    EbConfig_t                *config,
    EbAppContext_t            *callbackData,
    uint32_t                 instanceIdx)
{
    EbErrorType   return_error = EB_ErrorNone;

    // Assign Instance index to the library
    callbackData->instanceIdx = (uint8_t)instanceIdx;

    // Initialize Port Activity Flags
    callbackData->ebEncParameters.source_width       = config->sourceWidth;
    callbackData->ebEncParameters.source_height      = config->sourceHeight;
    callbackData->ebEncParameters.encoder_bit_depth   = config->encoderBitDepth;
    //callbackData->ebEncParameters.code_vps_sps_pps     = 1;
    //callbackData->ebEncParameters.code_eos_nal        = 1;
    callbackData->ebEncParameters.recon_enabled      = config->reconFile ? 1 : 0;

    return return_error;

}

/***********************************
 * Initialize Core & Component
 ***********************************/
EbErrorType InitEncoder(
    EbConfig_t                *config,
    EbAppContext_t            *callbackData,
    uint32_t                 instanceIdx)
{
    EbErrorType        return_error = EB_ErrorNone;

    ///************************* LIBRARY INIT [START] *********************///
    // STEP 1: Call the library to construct a Component Handle
    return_error = eb_init_handle(&callbackData->svtEncoderHandle, callbackData, &callbackData->ebEncParameters);
    if (return_error != EB_ErrorNone) {return return_error;}

    // STEP 3: Copy all configuration parameters into the callback structure
    return_error = CopyConfigurationParameters(config,callbackData,instanceIdx);
    if (return_error != EB_ErrorNone) { return return_error; }

    // STEP 4: Send over all configuration parameters
    return_error = eb_svt_enc_set_parameter(callbackData->svtEncoderHandle,&callbackData->ebEncParameters);
    if (return_error != EB_ErrorNone) { return return_error; }

    // STEP 5: Init Encoder
    return_error = eb_init_encoder(callbackData->svtEncoderHandle);
    // Get ivf header
    if (config->bitstreamFile) {
        EbBufferHeaderType *outputStreamBuffer;
        return_error = eb_svt_enc_stream_header(callbackData->svtEncoderHandle, &outputStreamBuffer);
        if (return_error != EB_ErrorNone) {
            return return_error;
        }
        fwrite(outputStreamBuffer->p_buffer, 1, outputStreamBuffer->n_filled_len, config->bitstreamFile);
    }
    ///************************* LIBRARY INIT [END] *********************///
    return return_error;
}

/***********************************
 * Deinit Components
 ***********************************/
EbErrorType DeInitEncoder(
    EbAppContext_t *callbackDataPtr,
    uint32_t        instanceIndex)
{
    (void)instanceIndex;
    EbErrorType return_error = EB_ErrorNone;

    if (((EbComponentType*)(callbackDataPtr->svtEncoderHandle)) != NULL) {
            return_error = eb_deinit_encoder(callbackDataPtr->svtEncoderHandle);
    }

    // Destruct the buffer memory pool
    if (return_error != EB_ErrorNone) { return return_error; }

    // Destruct the component
    eb_deinit_handle(callbackDataPtr->svtEncoderHandle);

    return return_error;
}
