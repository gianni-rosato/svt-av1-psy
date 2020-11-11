/*
* Copyright(c) 2019 Intel Corporation
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
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

#if FEATURE_INL_ME
typedef struct InLoopMeContext {
    EbFifo *   input_fifo_ptr;
    EbFifo *   output_fifo_ptr;
    MeContext *me_context_ptr;

    uint8_t *index_table0;
    uint8_t *index_table1;
} InLoopMeContext;
#endif
/***************************************
 * Extern Function Declaration
 ***************************************/
EbErrorType motion_estimation_context_ctor(EbThreadContext *  thread_context_ptr,
                                           const EbEncHandle *enc_handle_ptr, int index);

extern void *motion_estimation_kernel(void *input_ptr);

#if FEATURE_INL_ME
EbErrorType ime_context_ctor(EbThreadContext *thread_context_ptr,
                             const EbEncHandle *enc_handle_ptr, int index);
extern void *inloop_me_kernel(void *input_ptr);
#endif
EbErrorType signal_derivation_me_kernel_oq(SequenceControlSet *       scs_ptr,
                                           PictureParentControlSet *  pcs_ptr,
                                           MotionEstimationContext_t *context_ptr);

#endif // EbMotionEstimationProcess_h
