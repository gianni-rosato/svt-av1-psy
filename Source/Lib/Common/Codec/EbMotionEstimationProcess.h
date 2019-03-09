/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbEstimationProcess_h
#define EbEstimationProcess_h

#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"
#include "EbSequenceControlSet.h"
#include "EbPictureControlSet.h"
#include "EbMotionEstimationContext.h"

/**************************************
 * Context
 **************************************/
typedef struct MotionEstimationContext_s
{
    EbFifo_t                        *pictureDecisionResultsInputFifoPtr;
    EbFifo_t                        *motionEstimationResultsOutputFifoPtr;
    IntraReferenceSamplesOpenLoop_t *intra_ref_ptr;
    MeContext_t                     *me_context_ptr;

    uint8_t                       *indexTable0;
    uint8_t                       *indexTable1;

} MotionEstimationContext_t;

/***************************************
 * Extern Function Declaration
 ***************************************/
extern EbErrorType MotionEstimationContextCtor(
    MotionEstimationContext_t   **context_dbl_ptr,
    EbFifo_t                     *pictureDecisionResultsInputFifoPtr,
    EbFifo_t                     *motionEstimationResultsOutputFifoPtr);


extern void* MotionEstimationKernel(void *input_ptr);

#endif // EbMotionEstimationProcess_h