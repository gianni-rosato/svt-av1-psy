/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

/*
* Copyright (c) 2016, Alliance for Open Media. All rights reserved
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at www.aomedia.org/license/software. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at www.aomedia.org/license/patent.
*/

#ifndef EbSourceBasedOperations_h
#define EbSourceBasedOperations_h

#include "EbDefinitions.h"

/***************************************
 * Extern Function Declaration
 ***************************************/
EbErrorType source_based_operations_context_ctor(EbThreadContext *  thread_context_ptr,
                                                 const EbEncHandle *enc_handle_ptr, int index);

extern void *source_based_operations_kernel(void *input_ptr);

#endif // EbSourceBasedOperations_h
