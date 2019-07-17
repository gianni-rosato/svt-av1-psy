/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbReferenceObject_h
#define EbReferenceObject_h

#include "EbDefinitions.h"
#include "EbDefinitions.h"
#include "EbAdaptiveMotionVectorPrediction.h"
#include "EbObject.h"

typedef struct EbReferenceObject
{
    EbDctor                      dctor;
    EbPictureBufferDesc          *reference_picture;
    EbPictureBufferDesc          *reference_picture16bit;
    uint64_t                        ref_poc;
#if ADD_DELTA_QP_SUPPORT
    uint16_t                        qp;
#else
    uint8_t                         qp;
#endif
    EB_SLICE                        slice_type;
    uint8_t                         intra_coded_area;//percentage of intra coded area 0-100%
    uint8_t                         intra_coded_area_sb[MAX_NUMBER_OF_TREEBLOCKS_PER_PICTURE];//percentage of intra coded area 0-100%
    uint32_t                        non_moving_index_array[MAX_NUMBER_OF_TREEBLOCKS_PER_PICTURE];//array to hold non-moving blocks in reference frames
    uint32_t                        picSampleValue[MAX_NUMBER_OF_REGIONS_IN_WIDTH][MAX_NUMBER_OF_REGIONS_IN_HEIGHT][3];// [Y U V];
    uint8_t                         tmp_layer_idx;
    EbBool                          is_scene_change;
    uint16_t                        pic_avg_variance;
    uint8_t                         average_intensity;
    aom_film_grain_t                film_grain_params; //Film grain parameters for a reference frame
    uint32_t                        cdef_frame_strength;
    int8_t                          sg_frame_ep;
#if ENABLE_CDF_UPDATE
    FRAME_CONTEXT                   frame_context;
    EbWarpedMotionParams            global_motion[TOTAL_REFS_PER_FRAME];
#endif
} EbReferenceObject;

typedef struct EbReferenceObjectDescInitData {
    EbPictureBufferDescInitData   reference_picture_desc_init_data;
} EbReferenceObjectDescInitData;

typedef struct EbPaReferenceObject
{
    EbDctor                      dctor;
    EbPictureBufferDesc          *input_padded_picture_ptr;
    EbPictureBufferDesc          *quarter_decimated_picture_ptr;
    EbPictureBufferDesc          *sixteenth_decimated_picture_ptr;
    EbPictureBufferDesc          *quarter_filtered_picture_ptr;
    EbPictureBufferDesc          *sixteenth_filtered_picture_ptr;
    uint16_t                      variance[MAX_NUMBER_OF_TREEBLOCKS_PER_PICTURE];
    uint8_t                       y_mean[MAX_NUMBER_OF_TREEBLOCKS_PER_PICTURE];
    EB_SLICE                      slice_type;
    uint32_t                      dependent_pictures_count; //number of pic using this reference frame

} EbPaReferenceObject;

typedef struct EbPaReferenceObjectDescInitData
{
    EbPictureBufferDescInitData   reference_picture_desc_init_data;
    EbPictureBufferDescInitData   quarter_picture_desc_init_data;
    EbPictureBufferDescInitData   sixteenth_picture_desc_init_data;
} EbPaReferenceObjectDescInitData;

/**************************************
 * Extern Function Declarations
 **************************************/
extern EbErrorType eb_reference_object_creator(
    EbPtr *object_dbl_ptr,
    EbPtr  object_init_data_ptr);

extern EbErrorType eb_pa_reference_object_creator(
    EbPtr *object_dbl_ptr,
    EbPtr  object_init_data_ptr);

#endif //EbReferenceObject_h
