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
typedef struct MotionEstimationContext
{
    EbFifo                        *picture_decision_results_input_fifo_ptr;
    EbFifo                        *motion_estimation_results_output_fifo_ptr;
    IntraReferenceSamplesOpenLoop *intra_ref_ptr;
    MeContext                     *me_context_ptr;

    uint8_t                       *indexTable0;
    uint8_t                       *indexTable1;
} MotionEstimationContext_t;

/***************************************
 * Extern Function Declaration
 ***************************************/
extern EbErrorType motion_estimation_context_ctor(
    MotionEstimationContext_t   **context_dbl_ptr,
    EbFifo                     *picture_decision_results_input_fifo_ptr,
    EbFifo                     *motion_estimation_results_output_fifo_ptr);


extern void* motion_estimation_kernel(void *input_ptr);

#endif // EbMotionEstimationProcess_h