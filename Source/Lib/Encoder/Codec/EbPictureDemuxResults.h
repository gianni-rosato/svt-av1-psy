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

#ifndef EbPictureResults_h
#define EbPictureResults_h

#include "EbSystemResourceManager.h"
#include "EbObject.h"

/**************************************
 * Enums
 **************************************/
typedef enum EbPicType {
    EB_PIC_INVALID        = 0,
    EB_PIC_INPUT          = 1,
    EB_PIC_REFERENCE      = 2,
    EB_PIC_FEEDBACK       = 3,
    EB_PIC_SUPERRES_INPUT = 4
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
    EbObjectWrapper           *reference_picture_wrapper_ptr;
    struct SequenceControlSet *scs_ptr;
    uint64_t                   picture_number;
    uint64_t                   decode_order;
} PictureDemuxResults;

typedef struct PictureResultInitData {
    int32_t junk;
} PictureResultInitData;

/**************************************
 * Extern Function Declarations
 **************************************/
extern EbErrorType picture_results_creator(EbPtr *object_dbl_ptr, EbPtr object_init_data_ptr);

typedef struct PictureManagerResults {
    EbDctor          dctor;
    EbObjectWrapper *pcs_wrapper_ptr;
    uint32_t         segment_index;
    uint8_t          task_type;
    uint8_t          tpl_ref_list0_count;
    uint8_t          tpl_ref_list1_count;
    uint8_t          temporal_layer_index;
    Bool             is_used_as_reference_flag;
} PictureManagerResults;

typedef struct PictureManagerResultInitData {
    int32_t junk;
} PictureManagerResultInitData;
#endif //EbPictureResults_h
