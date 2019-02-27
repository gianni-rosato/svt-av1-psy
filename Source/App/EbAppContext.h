/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbAppContext_h
#define EbAppContext_h

#include "EbApi.h"
#include "EbAppConfig.h"

/***************************************

 * App Callback data struct
 ***************************************/
typedef struct EbAppContext_s {
    void                               *cmdSemaphoreHandle;
    void                               *inputSemaphoreHandle;
    void                               *streamSemaphoreHandle;
    EbSvtAv1EncConfiguration              ebEncParameters;

    // Output Ports Active Flags
    APPPORTACTIVETYPE                   outputStreamPortActive;

    // Component Handle
    EbComponentType*                   svtEncoderHandle;

    // Buffer Pools
    EbBufferHeaderType                *inputBufferPool;
    EbBufferHeaderType                *streamBufferPool;
    EbBufferHeaderType                *recon_buffer;

    // Instance Index
    uint8_t                            instanceIdx;

} EbAppContext_t;


/********************************
 * External Function
 ********************************/
extern EbErrorType InitEncoder(EbConfig_t *config, EbAppContext_t *callbackData, uint32_t instanceIdx);
extern EbErrorType DeInitEncoder(EbAppContext_t *callbackDataPtr, uint32_t instanceIndex);

#endif // EbAppContext_h