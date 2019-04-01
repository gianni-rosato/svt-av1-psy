/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbReferenceObject_h
#define EbReferenceObject_h

#include "EbDefinitions.h"
#include "EbDefinitions.h"
#include "EbAdaptiveMotionVectorPrediction.h"

typedef struct EbReferenceObject 
{
    EbPictureBufferDesc_t          *reference_picture;
    EbPictureBufferDesc_t          *reference_picture16bit;
    EbPictureBufferDesc_t          *ref_den_src_picture;

    TmvpUnit_t                     *tmvp_map;
    EbBool                          tmvp_enable_flag;
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
    EbBool                          penalize_skipflag;
    uint8_t                         tmp_layer_idx;
    EbBool                          is_scene_change;
    uint16_t                        pic_avg_variance;
    uint8_t                         average_intensity;
    aom_film_grain_t                film_grain_params; //Film grain parameters for a reference frame
    uint32_t                        cdef_frame_strength;
    int8_t                          sg_frame_ep;
} EbReferenceObject;

typedef struct EbReferenceObjectDescInitData {
    EbPictureBufferDescInitData_t   reference_picture_desc_init_data;
} EbReferenceObjectDescInitData;

typedef struct EbPaReferenceObject 
{
    EbPictureBufferDesc_t          *input_padded_picture_ptr;
    EbPictureBufferDesc_t          *quarter_decimated_picture_ptr;
    EbPictureBufferDesc_t          *sixteenth_decimated_picture_ptr;
    uint16_t                        variance[MAX_NUMBER_OF_TREEBLOCKS_PER_PICTURE];
    uint8_t                         y_mean[MAX_NUMBER_OF_TREEBLOCKS_PER_PICTURE];
    EB_SLICE                        slice_type;
    uint32_t                        dependent_pictures_count; //number of pic using this reference frame
    PictureParentControlSet_t      *p_pcs_ptr;

} EbPaReferenceObject;

typedef struct EbPaReferenceObjectDescInitData 
{
    EbPictureBufferDescInitData_t   reference_picture_desc_init_data;
    EbPictureBufferDescInitData_t   quarter_picture_desc_init_data;
    EbPictureBufferDescInitData_t   sixteenth_picture_desc_init_data;
} EbPaReferenceObjectDescInitData;

/**************************************
 * Extern Function Declarations
 **************************************/
extern EbErrorType eb_reference_object_ctor(
    EbPtr *object_dbl_ptr,
    EbPtr  object_init_data_ptr);

extern EbErrorType eb_pa_reference_object_ctor(
    EbPtr *object_dbl_ptr,
    EbPtr  object_init_data_ptr);


#endif //EbReferenceObject_h