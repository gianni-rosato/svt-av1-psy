/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbInitialRateControl_h
#define EbInitialRateControl_h

#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"
#include "EbRateControlProcess.h"



/**************************************
 * Context
 **************************************/
typedef struct InitialRateControlContext_s
{
    EbFifo_t                    *motionEstimationResultsInputFifoPtr;
    EbFifo_t                    *initialrateControlResultsOutputFifoPtr;

} InitialRateControlContext_t;

/***************************************
 * Extern Function Declaration
 ***************************************/
extern EbErrorType InitialRateControlContextCtor(
    InitialRateControlContext_t **context_dbl_ptr,
    EbFifo_t                     *motionEstimationResultsInputFifoPtr,
    EbFifo_t                     *picture_demux_results_output_fifo_ptr);

extern void* InitialRateControlKernel(void *input_ptr);

#endif // EbInitialRateControl_h