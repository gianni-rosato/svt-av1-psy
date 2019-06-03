/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbPictureAnalysisResults_h
#define EbPictureAnalysisResults_h

#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"

/**************************************
 * Process Results
 **************************************/
typedef struct PictureAnalysisResults {
    EbObjectWrapper *picture_control_set_wrapper_ptr;
} PictureAnalysisResults;

typedef struct PictureAnalysisResultInitData {
    int32_t junk;
} PictureAnalysisResultInitData;

/**************************************
 * Extern Function Declarations
 **************************************/
extern EbErrorType picture_analysis_result_ctor(
    EbPtr *object_dbl_ptr,
    EbPtr  object_init_data_ptr);

#endif //EbPictureAnalysisResults_h
