/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbRestProcess_h
#define EbRestProcess_h

#include "EbDefinitions.h"

/**************************************
 * Extern Function Declarations
 **************************************/
extern EbErrorType rest_context_ctor(EbThreadContext *  thread_context_ptr,
                                     const EbEncHandle *enc_handle_ptr, int index, int demux_index);

extern void *rest_kernel(void *input_ptr);

#endif
