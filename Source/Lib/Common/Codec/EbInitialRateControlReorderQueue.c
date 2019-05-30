/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>
#include "EbInitialRateControlReorderQueue.h"

EbErrorType initial_rate_control_reorder_entry_ctor(
    InitialRateControlReorderEntry   **entry_dbl_ptr,
    uint32_t                          picture_number)
{
    EB_MALLOC(InitialRateControlReorderEntry*, *entry_dbl_ptr, sizeof(InitialRateControlReorderEntry), EB_N_PTR);

    (*entry_dbl_ptr)->picture_number = picture_number;
    (*entry_dbl_ptr)->parent_pcs_wrapper_ptr = (EbObjectWrapper *)EB_NULL;

    return EB_ErrorNone;
}

EbErrorType hl_rate_control_histogram_entry_ctor(
    HlRateControlHistogramEntry   **entry_dbl_ptr,
    uint32_t                          picture_number)
{
    EB_MALLOC(HlRateControlHistogramEntry*, *entry_dbl_ptr, sizeof(HlRateControlHistogramEntry), EB_N_PTR);

    (*entry_dbl_ptr)->picture_number = picture_number;
#if RC
    (*entry_dbl_ptr)->life_count = 0;
#else
    (*entry_dbl_ptr)->life_count = 0;
#endif
    (*entry_dbl_ptr)->parent_pcs_wrapper_ptr = (EbObjectWrapper *)EB_NULL;

    // ME and OIS Distortion Histograms
    EB_MALLOC(uint16_t*, (*entry_dbl_ptr)->me_distortion_histogram, sizeof(uint16_t) * NUMBER_OF_SAD_INTERVALS, EB_N_PTR);

    EB_MALLOC(uint16_t*, (*entry_dbl_ptr)->ois_distortion_histogram, sizeof(uint16_t) * NUMBER_OF_INTRA_SAD_INTERVALS, EB_N_PTR);

    return EB_ErrorNone;
}
