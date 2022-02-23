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
typedef struct PictureManagerContext {
    EbDctor  dctor;
    EbFifo  *picture_input_fifo_ptr;
    EbFifo  *picture_manager_output_fifo_ptr;
    EbFifo  *picture_control_set_fifo_ptr;
    EbFifo  *recon_coef_fifo_ptr;
    uint64_t pmgr_dec_order;
} PictureManagerContext;
/***************************************
     * Extern Function Declaration
     ***************************************/
EbErrorType picture_manager_context_ctor(EbThreadContext   *thread_context_ptr,
                                         const EbEncHandle *enc_handle_ptr, int rate_control_index);

extern void *picture_manager_kernel(void *input_ptr);

#ifdef __cplusplus
}
#endif
#endif // EbPictureManager_h
