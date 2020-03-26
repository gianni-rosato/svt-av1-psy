/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbPictureAnalysis_h
#define EbPictureAnalysis_h

#include "EbDefinitions.h"

#include "EbPictureControlSet.h"

/***************************************
 * Extern Function Declaration
 ***************************************/
EbErrorType picture_analysis_context_ctor(EbThreadContext *  thread_context_ptr,
                                          const EbEncHandle *enc_handle_ptr, int index);

extern void *picture_analysis_kernel(void *input_ptr);


void downsample_filtering_input_picture(PictureParentControlSet *pcs_ptr,
                                        EbPictureBufferDesc *    input_padded_picture_ptr,
                                        EbPictureBufferDesc *    quarter_picture_ptr,
                                        EbPictureBufferDesc *    sixteenth_picture_ptr);

#endif // EbPictureAnalysis_h
