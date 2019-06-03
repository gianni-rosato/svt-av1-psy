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
    EbPictureBufferDesc          *reference_picture;
    EbPictureBufferDesc          *reference_picture16bit;
#if !OPT_LOSSLESS_1
    EbPictureBufferDesc          *ref_den_src_picture;

    TmvpUnit                     *tmvp_map;
    EbBool                          tmvp_enable_flag;
#endif
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
#if !DISABLE_OIS_USE
    EbBool                          penalize_skipflag;
#endif
    uint8_t                         tmp_layer_idx;
    EbBool                          is_scene_change;
    uint16_t                        pic_avg_variance;
    uint8_t                         average_intensity;
    aom_film_grain_t                film_grain_params; //Film grain parameters for a reference frame
    uint32_t                        cdef_frame_strength;
    int8_t                          sg_frame_ep;
} EbReferenceObject;

typedef struct EbReferenceObjectDescInitData {
    EbPictureBufferDescInitData   reference_picture_desc_init_data;
} EbReferenceObjectDescInitData;

typedef struct EbPaReferenceObject
{
    EbPictureBufferDesc          *input_padded_picture_ptr;
    EbPictureBufferDesc          *quarter_decimated_picture_ptr;
    EbPictureBufferDesc          *sixteenth_decimated_picture_ptr;
    uint16_t                        variance[MAX_NUMBER_OF_TREEBLOCKS_PER_PICTURE];
    uint8_t                         y_mean[MAX_NUMBER_OF_TREEBLOCKS_PER_PICTURE];
    EB_SLICE                        slice_type;
    uint32_t                        dependent_pictures_count; //number of pic using this reference frame
#if !BUG_FIX_PCS_LIVE_COUNT
    PictureParentControlSet      *p_pcs_ptr;
#endif

#if BUG_FIX_INPUT_LIVE_COUNT
    EbObjectWrapper              *input_picture_wrapper_ptr;
#endif
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
extern EbErrorType eb_reference_object_ctor(
    EbPtr *object_dbl_ptr,
    EbPtr  object_init_data_ptr);

extern EbErrorType eb_pa_reference_object_ctor(
    EbPtr *object_dbl_ptr,
    EbPtr  object_init_data_ptr);

#endif //EbReferenceObject_h
