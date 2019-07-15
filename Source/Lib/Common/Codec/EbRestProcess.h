/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbRestProcess_h
#define EbRestProcess_h

#include "EbDefinitions.h"

#include "EbSystemResourceManager.h"
#include "EbPictureBufferDesc.h"
#include "EbSequenceControlSet.h"
#include "EbUtility.h"
#include "EbPsnr.h"
#include "EbPictureControlSet.h"
#include "EbObject.h"

/**************************************
 * Rest Context
 **************************************/
typedef struct RestContext
{
    EbDctor                       dctor;
    EbFifo                       *rest_input_fifo_ptr;
    EbFifo                       *rest_output_fifo_ptr;
    EbFifo                       *picture_demux_fifo_ptr;

    EbPictureBufferDesc          *trial_frame_rst;

    EbPictureBufferDesc          *temp_lf_recon_picture_ptr;
    EbPictureBufferDesc          *temp_lf_recon_picture16bit_ptr;

    EbPictureBufferDesc           *org_rec_frame; // while doing the filtering recon gets updated uisng setup/restore processing_stripe_bounadaries
                                                    // many threads doing the above will result in race condition.
                                                    // each thread will hence have his own copy of recon to work on.
                                                    // later we can have a search version that does not need the exact right recon
    int32_t *rst_tmpbuf;
} RestContext;

/**************************************
 * Extern Function Declarations
 **************************************/
extern EbErrorType rest_context_ctor(
    RestContext                  *context_ptr,
    EbFifo                       *rest_input_fifo_ptr,
    EbFifo                       *rest_output_fifo_ptr,
    EbFifo                      *picture_demux_fifo_ptr,
    EbBool                  is16bit,
    EbColorFormat           color_format,
    uint32_t                max_input_luma_width,
    uint32_t                max_input_luma_height
   );

extern void* rest_kernel(void *input_ptr);

#endif
