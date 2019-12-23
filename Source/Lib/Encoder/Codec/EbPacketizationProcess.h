/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbPacketization_h
#define EbPacketization_h

#include "EbDefinitions.h"
#ifdef __cplusplus
extern "C" {
#endif

/**************************************
     * Extern Function Declarations
     **************************************/
EbErrorType packetization_context_ctor(EbThreadContext *  thread_context_ptr,
                                       const EbEncHandle *enc_handle_ptr, int rate_control_index,
                                       int demux_index);

extern void *packetization_kernel(void *input_ptr);
#ifdef __cplusplus
}
#endif
#endif // EbPacketization_h
