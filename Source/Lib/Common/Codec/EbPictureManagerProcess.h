/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbPictureManager_h
#define EbPictureManager_h

#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"
#include "EbObject.h"
#ifdef __cplusplus
extern "C" {
#endif
    /***************************************
     * Context
     ***************************************/
    typedef struct PictureManagerContext
    {
        EbDctor  dctor;
        EbFifo  *picture_input_fifo_ptr;
        EbFifo  *picture_manager_output_fifo_ptr;
        EbFifo **picture_control_set_fifo_ptr_array;
    } PictureManagerContext;

    /***************************************
     * Extern Function Declaration
     ***************************************/
    extern EbErrorType picture_manager_context_ctor(
        PictureManagerContext  *context_ptr,
        EbFifo                 *pictureInputFifoPtr,
        EbFifo                 *picture_manager_output_fifo_ptr,
        EbFifo                **picture_control_set_fifo_ptr_array);

    extern void* picture_manager_kernel(void *input_ptr);

#ifdef __cplusplus
}
#endif
#endif // EbPictureManager_h
