// clang-format off
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
EbErrorType picture_decision_context_ctor(
    EbThreadContext     *thread_context_ptr,
    const EbEncHandle   *enc_handle_ptr);

extern void* picture_decision_kernel(void *input_ptr);

void DownsampleDecimationInputPicture(
    PictureParentControlSet *picture_control_set_ptr,
    EbPictureBufferDesc     *inputPaddedPicturePtr,
    EbPictureBufferDesc     *quarterDecimatedPicturePtr,
    EbPictureBufferDesc     *sixteenthDecimatedPicturePtr);

void PadPictureToMultipleOfMinCuSizeDimensions(
        SequenceControlSet            *sequence_control_set_ptr,
        EbPictureBufferDesc           *input_picture_ptr);
void PicturePreProcessingOperations(
    PictureParentControlSet       *picture_control_set_ptr,
    SequenceControlSet            *sequence_control_set_ptr,
    uint32_t                       sb_total_count);
void PadPictureToMultipleOfLcuDimensions(
        EbPictureBufferDesc   *input_padded_picture_ptr);

void GatheringPictureStatistics(
        SequenceControlSet            *sequence_control_set_ptr,
        PictureParentControlSet       *picture_control_set_ptr,
        EbPictureBufferDesc           *input_picture_ptr,
        EbPictureBufferDesc           *input_padded_picture_ptr,
        EbPictureBufferDesc           *sixteenth_decimated_picture_ptr,
        uint32_t                      sb_total_count);

void DownSampleChroma(EbPictureBufferDesc* input_picture_ptr,
                      EbPictureBufferDesc* outputPicturePtr);

#endif // EbPictureDecision_h
// clang-format on
