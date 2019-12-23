/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbCdefProcess_h
#define EbCdefProcess_h

#include "EbSystemResourceManager.h"
#include "EbObject.h"

/**************************************
 * Extern Function Declarations
 **************************************/
extern EbErrorType cdef_context_ctor(EbThreadContext *  thread_context_ptr,
                                     const EbEncHandle *enc_handle_ptr, int index);

extern void *cdef_kernel(void *input_ptr);

#endif
