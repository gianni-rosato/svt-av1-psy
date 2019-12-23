/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbResourceCoordination_h
#define EbResourceCoordination_h

#include "EbDefinitions.h"
#ifdef __cplusplus
extern "C" {
#endif
/***************************************
     * Extern Function Declaration
     ***************************************/
EbErrorType resource_coordination_context_ctor(EbThreadContext* thread_context_ptr,
                                               EbEncHandle*     enc_handle_ptr);

extern void* resource_coordination_kernel(void* input_ptr);
#ifdef __cplusplus
}
#endif
#endif // EbResourceCoordination_h
