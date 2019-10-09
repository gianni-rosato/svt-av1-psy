/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbDecRestoration_h
#define EbDecRestoration_h

#ifdef __cplusplus
extern "C" {
#endif

#include "EbDecHandle.h"

void dec_av1_loop_restoration_filter_frame(EbDecHandle *dec_handle,
                                           int optimized_lr);
void dec_av1_loop_restoration_save_boundary_lines(EbDecHandle *dec_handle,
                                                  int after_cdef);

#ifdef __cplusplus
}
#endif

#endif
