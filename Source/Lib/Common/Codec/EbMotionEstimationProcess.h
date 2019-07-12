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
#include "EbObject.h"

/**************************************
 * Context
 **************************************/
typedef struct MotionEstimationContext
{
    EbDctor                       dctor;
    EbFifo                        *picture_decision_results_input_fifo_ptr;
    EbFifo                        *motion_estimation_results_output_fifo_ptr;
    MeContext                     *me_context_ptr;

    uint8_t                       *indexTable0;
    uint8_t                       *indexTable1;
} MotionEstimationContext_t;

/***************************************
 * Extern Function Declaration
 ***************************************/
extern EbErrorType motion_estimation_context_ctor(
    MotionEstimationContext_t  *context_ptr,
    EbFifo                     *picture_decision_results_input_fifo_ptr,
    EbFifo                     *motion_estimation_results_output_fifo_ptr,
    uint16_t                    max_input_luma_width,
    uint16_t                    max_input_luma_height,
    uint8_t                     nsq_present,
    uint8_t                     mrp_mode);

extern void* motion_estimation_kernel(void *input_ptr);

EbErrorType signal_derivation_me_kernel_oq(SequenceControlSet        *sequence_control_set_ptr,
                                           PictureParentControlSet   *picture_control_set_ptr,
                                           MotionEstimationContext_t *context_ptr);

#endif // EbMotionEstimationProcess_h
