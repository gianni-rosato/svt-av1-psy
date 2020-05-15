/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbPictureDecision_h
#define EbPictureDecision_h

#include "EbDefinitions.h"
#include "EbPictureControlSet.h"
#include "EbSequenceControlSet.h"

/***************************************
 * Extern Function Declaration
 ***************************************/
EbErrorType picture_decision_context_ctor(EbThreadContext *  thread_context_ptr,
                                          const EbEncHandle *enc_handle_ptr);

extern void *picture_decision_kernel(void *input_ptr);

void downsample_decimation_input_picture(PictureParentControlSet *pcs_ptr,
                                         EbPictureBufferDesc *    inputPaddedPicturePtr,
                                         EbPictureBufferDesc *    quarterDecimatedPicturePtr,
                                         EbPictureBufferDesc *    sixteenthDecimatedPicturePtr);

void pad_picture_to_multiple_of_min_blk_size_dimensions(SequenceControlSet * scs_ptr,
                                                        EbPictureBufferDesc *input_picture_ptr);
void pad_picture_to_multiple_of_min_blk_size_dimensions_16bit(
    SequenceControlSet * scs_ptr, EbPictureBufferDesc *input_picture_ptr);
void picture_pre_processing_operations(PictureParentControlSet *pcs_ptr,
                                       SequenceControlSet *scs_ptr, uint32_t sb_total_count);
void pad_picture_to_multiple_of_sb_dimensions(EbPictureBufferDesc *input_padded_picture_ptr);

void gathering_picture_statistics(SequenceControlSet *scs_ptr, PictureParentControlSet *pcs_ptr,
                                  EbPictureBufferDesc *input_picture_ptr,
                                  EbPictureBufferDesc *input_padded_picture_ptr,
                                  EbPictureBufferDesc *sixteenth_decimated_picture_ptr,
                                  uint32_t             sb_total_count);

void down_sample_chroma(EbPictureBufferDesc *input_picture_ptr,
                        EbPictureBufferDesc *outputPicturePtr);

#endif // EbPictureDecision_h
