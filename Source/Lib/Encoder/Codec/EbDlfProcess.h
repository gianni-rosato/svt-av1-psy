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

#ifndef EbDlfProcess_h
#define EbDlfProcess_h

#include "EbSystemResourceManager.h"
#include "EbObject.h"
#include "EbPictureBufferDesc.h"
#include "EbSvtAv1Formats.h"

/**************************************
 * Dlf Context
 **************************************/
typedef struct DlfContext {
    EbFifo *             dlf_input_fifo_ptr;
    EbFifo *             dlf_output_fifo_ptr;
    EbPictureBufferDesc *temp_lf_recon_picture_ptr;
    EbPictureBufferDesc *temp_lf_recon_picture16bit_ptr;
} DlfContext;

/**************************************
 * Extern Function Declarations
 **************************************/
extern EbErrorType dlf_context_ctor(EbThreadContext *  thread_context_ptr,
                                    const EbEncHandle *enc_handle_ptr, int index);

extern void *dlf_kernel(void *input_ptr);

#endif // EbEntropyCodingProcess_h
