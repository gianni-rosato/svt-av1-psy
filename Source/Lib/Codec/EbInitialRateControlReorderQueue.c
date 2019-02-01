/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>
#include "EbInitialRateControlReorderQueue.h"

EbErrorType InitialRateControlReorderEntryCtor(
    InitialRateControlReorderEntry_t   **entryDblPtr,
    uint32_t                          picture_number)
{
    EB_MALLOC(InitialRateControlReorderEntry_t*, *entryDblPtr, sizeof(InitialRateControlReorderEntry_t), EB_N_PTR);

    (*entryDblPtr)->picture_number = picture_number;
    (*entryDblPtr)->parentPcsWrapperPtr = (EbObjectWrapper_t *)EB_NULL;

    return EB_ErrorNone;
}


EbErrorType HlRateControlHistogramEntryCtor(
    HlRateControlHistogramEntry_t   **entryDblPtr,
    uint32_t                          picture_number)
{
    EB_MALLOC(HlRateControlHistogramEntry_t*, *entryDblPtr, sizeof(HlRateControlHistogramEntry_t), EB_N_PTR);

    (*entryDblPtr)->picture_number = picture_number;
    (*entryDblPtr)->lifeCount = 0;

    (*entryDblPtr)->parentPcsWrapperPtr = (EbObjectWrapper_t *)EB_NULL;

    // ME and OIS Distortion Histograms
    EB_MALLOC(uint16_t*, (*entryDblPtr)->me_distortion_histogram, sizeof(uint16_t) * NUMBER_OF_SAD_INTERVALS, EB_N_PTR);

    EB_MALLOC(uint16_t*, (*entryDblPtr)->ois_distortion_histogram, sizeof(uint16_t) * NUMBER_OF_INTRA_SAD_INTERVALS, EB_N_PTR);

    return EB_ErrorNone;
}



