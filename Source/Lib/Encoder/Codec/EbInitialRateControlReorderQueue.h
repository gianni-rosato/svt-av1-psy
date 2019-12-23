/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbInitialRateControlReorderQueue_h
#define EbInitialRateControlReorderQueue_h

#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"
#include "EbPictureControlSet.h"
#include "EbObject.h"
/************************************************
 * Initial Rate Control Reorder Queue Entry
 ************************************************/
typedef struct InitialRateControlReorderEntry {
    EbDctor          dctor;
    uint64_t         picture_number;
    EbObjectWrapper *parent_pcs_wrapper_ptr;
} InitialRateControlReorderEntry;

extern EbErrorType initial_rate_control_reorder_entry_ctor(
    InitialRateControlReorderEntry *entry_ptr, uint32_t picture_number);

/************************************************
 * High Level Rate Control Histogram Queue Entry
 ************************************************/
typedef struct HlRateControlHistogramEntry {
    EbDctor          dctor;
    uint64_t         picture_number;
    int16_t          life_count;
    EbBool           passed_to_hlrc;
    EbBool           is_coded;
    uint64_t         total_num_bits_coded;
    EbObjectWrapper *parent_pcs_wrapper_ptr;
    EbBool           end_of_sequence_flag;
    uint64_t         pred_bits_ref_qp[MAX_REF_QP_NUM];
    EB_SLICE         slice_type;
    uint32_t         temporal_layer_index;
    uint32_t         frames_in_sw;

    // Motion Estimation Distortion and OIS Historgram
    uint16_t *me_distortion_histogram;
    uint16_t *ois_distortion_histogram;
    uint32_t  full_sb_count;
} HlRateControlHistogramEntry;

extern EbErrorType hl_rate_control_histogram_entry_ctor(HlRateControlHistogramEntry *entry_ptr,
                                                        uint32_t picture_number);

#endif //EbInitialRateControlReorderQueue_h
