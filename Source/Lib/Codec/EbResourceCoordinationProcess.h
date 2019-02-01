/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbResourceCoordination_h
#define EbResourceCoordination_h

#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"
#include "EbDefinitions.h"
#ifdef __cplusplus
extern "C" {
#endif
    /***************************************
     * Context
     ***************************************/
    typedef struct ResourceCoordinationContext_s
    {

        EbFifo_t                            *input_buffer_fifo_ptr;
        EbFifo_t                            *resourceCoordinationResultsOutputFifoPtr;
        EbFifo_t                           **pictureControlSetFifoPtrArray;
        EbSequenceControlSetInstance_t     **sequenceControlSetInstanceArray;
        EbObjectWrapper_t                  **sequenceControlSetActiveArray;
        EbFifo_t                            *sequenceControlSetEmptyFifoPtr;
        EbCallback_t                       **appCallbackPtrArray;
        
        // Compute Segments
        uint32_t                              *computeSegmentsTotalCountArray;
        uint32_t                               encodeInstancesTotalCount;

        // Picture Number Array
        uint64_t                              *pictureNumberArray;

        uint64_t                               averageEncMod;
        uint8_t                                prevEncMod;
        int8_t                                prevEncModeDelta;
        uint8_t                                prevChangeCond;

        int64_t                               previousModeChangeBuffer;
        int64_t                               previousModeChangeFrameIn;
        int64_t                               previousBufferCheck1;
        int64_t                               previousFrameInCheck1;
        int64_t                               previousFrameInCheck2;
        int64_t                               previousFrameInCheck3;


        uint64_t                               curSpeed; // speed x 1000
        uint64_t                               prevsTimeSeconds;
        uint64_t                               prevsTimeuSeconds;
        int64_t                               prevFrameOut;

        uint64_t                               firstInPicArrivedTimeSeconds;
        uint64_t                               firstInPicArrivedTimeuSeconds;
        EbBool                              startFlag;



    } ResourceCoordinationContext_t;

    /***************************************
     * Extern Function Declaration
     ***************************************/
    extern EbErrorType ResourceCoordinationContextCtor(
        ResourceCoordinationContext_t      **context_dbl_ptr,
        EbFifo_t                            *input_buffer_fifo_ptr,
        EbFifo_t                            *resourceCoordinationResultsOutputFifoPtr,
        EbFifo_t                           **pictureControlSetFifoPtrArray,
        EbSequenceControlSetInstance_t     **sequenceControlSetInstanceArray,
        EbFifo_t                            *sequenceControlSetEmptyFifoPtr,
        EbCallback_t                       **appCallbackPtrArray,
        uint32_t                            *computeSegmentsTotalCountArray,
        uint32_t                             encodeInstancesTotalCount);

    extern void* ResourceCoordinationKernel(void *input_ptr);
#ifdef __cplusplus
}
#endif
#endif // EbResourceCoordination_h