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

#ifndef EbPictureAnalysis_h
#define EbPictureAnalysis_h

#include "EbDefinitions.h"

#include "EbPictureControlSet.h"
#if FEATURE_INL_ME
#include "EbSequenceControlSet.h"
#endif

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

#if FEATURE_INL_ME
void pad_input_pictures(SequenceControlSet *scs_ptr,
                               EbPictureBufferDesc *input_picture_ptr);
#endif

#endif // EbPictureAnalysis_h
