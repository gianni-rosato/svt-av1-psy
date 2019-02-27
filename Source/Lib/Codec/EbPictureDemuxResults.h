/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbPictureResults_h
#define EbPictureResults_h

#include "EbSystemResourceManager.h"

/**************************************
 * Enums
 **************************************/
typedef enum EB_PIC_TYPE {
    EB_PIC_INVALID = 0,
    EB_PIC_INPUT = 1,
    EB_PIC_REFERENCE = 2
} EB_PIC_TYPE;

/**************************************
 * Picture Demux Results
 **************************************/
typedef struct PictureDemuxResults_s
{
    EB_PIC_TYPE                    pictureType;

    // Only valid for input pictures
    EbObjectWrapper_t             *pictureControlSetWrapperPtr;

    // Only valid for reference pictures
    EbObjectWrapper_t             *reference_picture_wrapper_ptr;
    EbObjectWrapper_t             *sequence_control_set_wrapper_ptr;
    uint64_t                         picture_number;

} PictureDemuxResults_t;

typedef struct PictureResultInitData_s {
    int32_t junk;
} PictureResultInitData_t;

/**************************************
 * Extern Function Declarations
 **************************************/
extern EbErrorType picture_results_ctor(
    EbPtr *object_dbl_ptr,
    EbPtr  object_init_data_ptr);


#endif //EbPictureResults_h