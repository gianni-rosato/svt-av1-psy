/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbInitialRateControl_h
#define EbInitialRateControl_h

#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"
#include "EbRateControlProcess.h"

/***************************************
 * Extern Function Declaration
 ***************************************/
EbErrorType initial_rate_control_context_ctor(EbThreadContext *  thread_context_ptr,
                                              const EbEncHandle *enc_handle_ptr);

extern void *initial_rate_control_kernel(void *input_ptr);

#endif // EbInitialRateControl_h
