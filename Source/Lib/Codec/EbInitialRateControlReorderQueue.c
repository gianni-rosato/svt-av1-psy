/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>
#include "EbInitialRateControlReorderQueue.h"

EbErrorType InitialRateControlReorderEntryCtor(
    InitialRateControlReorderEntry_t   **entry_dbl_ptr,
    uint32_t                          picture_number)
{
    EB_MALLOC(InitialRateControlReorderEntry_t*, *entry_dbl_ptr, sizeof(InitialRateControlReorderEntry_t), EB_N_PTR);

    (*entry_dbl_ptr)->picture_number = picture_number;
    (*entry_dbl_ptr)->parentPcsWrapperPtr = (EbObjectWrapper_t *)EB_NULL;

    return EB_ErrorNone;
}


EbErrorType HlRateControlHistogramEntryCtor(
    HlRateControlHistogramEntry_t   **entry_dbl_ptr,
    uint32_t                          picture_number)
{
    EB_MALLOC(HlRateControlHistogramEntry_t*, *entry_dbl_ptr, sizeof(HlRateControlHistogramEntry_t), EB_N_PTR);

    (*entry_dbl_ptr)->picture_number = picture_number;
    (*entry_dbl_ptr)->lifeCount = 0;

    (*entry_dbl_ptr)->parentPcsWrapperPtr = (EbObjectWrapper_t *)EB_NULL;

    // ME and OIS Distortion Histograms
    EB_MALLOC(uint16_t*, (*entry_dbl_ptr)->me_distortion_histogram, sizeof(uint16_t) * NUMBER_OF_SAD_INTERVALS, EB_N_PTR);

    EB_MALLOC(uint16_t*, (*entry_dbl_ptr)->ois_distortion_histogram, sizeof(uint16_t) * NUMBER_OF_INTRA_SAD_INTERVALS, EB_N_PTR);

    return EB_ErrorNone;
}



