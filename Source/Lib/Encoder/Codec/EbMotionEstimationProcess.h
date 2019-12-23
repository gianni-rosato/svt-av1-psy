/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbEstimationProcess_h
#define EbEstimationProcess_h

#include "EbDefinitions.h"
#include "EbSequenceControlSet.h"
#include "EbMotionEstimationContext.h"

/**************************************
 * Context
 **************************************/
typedef struct MotionEstimationContext {
    EbFifo *   picture_decision_results_input_fifo_ptr;
    EbFifo *   motion_estimation_results_output_fifo_ptr;
    MeContext *me_context_ptr;

    uint8_t *index_table0;
    uint8_t *index_table1;
} MotionEstimationContext_t;

/***************************************
 * Extern Function Declaration
 ***************************************/
EbErrorType motion_estimation_context_ctor(EbThreadContext *  thread_context_ptr,
                                           const EbEncHandle *enc_handle_ptr, int index);

extern void *motion_estimation_kernel(void *input_ptr);

EbErrorType signal_derivation_me_kernel_oq(SequenceControlSet *       scs_ptr,
                                           PictureParentControlSet *  pcs_ptr,
                                           MotionEstimationContext_t *context_ptr);

#endif // EbMotionEstimationProcess_h
