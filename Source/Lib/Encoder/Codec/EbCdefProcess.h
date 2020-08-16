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
