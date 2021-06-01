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

#ifndef EbRestProcess_h
#define EbRestProcess_h

#include "EbDefinitions.h"

/**************************************
 * Extern Function Declarations
 **************************************/
extern EbErrorType rest_context_ctor(EbThreadContext *  thread_context_ptr,
#if OPT_MEMORY_REST
                                     const EbEncHandle *enc_handle_ptr,EbPtr object_init_data_ptr,  int index,int demux_index);
#else
                                     const EbEncHandle *enc_handle_ptr, int index, int demux_index);
#endif

extern void *rest_kernel(void *input_ptr);

#endif
