// clang-format off
/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
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
typedef struct DlfContext
{
    EbDctor              dctor;
    EbFifo              *dlf_input_fifo_ptr;
    EbFifo              *dlf_output_fifo_ptr;
    EbPictureBufferDesc *temp_lf_recon_picture_ptr;
    EbPictureBufferDesc *temp_lf_recon_picture16bit_ptr;
} DlfContext;

/**************************************
 * Extern Function Declarations
 **************************************/
extern EbErrorType dlf_context_ctor(
    DlfContext                   *context_ptr,
    EbFifo                       *dlf_input_fifo_ptr,
    EbFifo                       *dlf_output_fifo_ptr,
    EbBool                  is16bit,
    EbColorFormat           color_format,
    uint32_t                max_input_luma_width,
    uint32_t                max_input_luma_height
   );

extern void* dlf_kernel(void *input_ptr);

#endif // EbEntropyCodingProcess_h
// clang-format on
