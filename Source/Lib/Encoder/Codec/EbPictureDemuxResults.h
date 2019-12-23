/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbPictureResults_h
#define EbPictureResults_h

#include "EbSystemResourceManager.h"
#include "EbObject.h"

/**************************************
 * Enums
 **************************************/
typedef enum EbPicType {
    EB_PIC_INVALID   = 0,
    EB_PIC_INPUT     = 1,
    EB_PIC_REFERENCE = 2,
    EB_PIC_FEEDBACK  = 3
} EbPicType;

/**************************************
 * Picture Demux Results
 **************************************/
typedef struct PictureDemuxResults {
    EbDctor   dctor;
    EbPicType picture_type;

    // Only valid for input pictures
    EbObjectWrapper *pcs_wrapper_ptr;

    // Only valid for reference pictures
    EbObjectWrapper *reference_picture_wrapper_ptr;
    EbObjectWrapper *scs_wrapper_ptr;
    uint64_t         picture_number;
} PictureDemuxResults;

typedef struct PictureResultInitData {
    int32_t junk;
} PictureResultInitData;

/**************************************
 * Extern Function Declarations
 **************************************/
extern EbErrorType picture_results_creator(EbPtr *object_dbl_ptr, EbPtr object_init_data_ptr);

#endif //EbPictureResults_h
