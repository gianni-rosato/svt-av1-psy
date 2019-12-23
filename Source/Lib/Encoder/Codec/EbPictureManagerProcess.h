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
     * Extern Function Declaration
     ***************************************/
EbErrorType picture_manager_context_ctor(EbThreadContext *  thread_context_ptr,
                                         const EbEncHandle *enc_handle_ptr, int rate_control_index);

extern void *picture_manager_kernel(void *input_ptr);

#ifdef __cplusplus
}
#endif
#endif // EbPictureManager_h
