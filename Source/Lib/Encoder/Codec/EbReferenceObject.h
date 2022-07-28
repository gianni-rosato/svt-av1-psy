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

#ifndef EbReferenceObject_h
#define EbReferenceObject_h

#include "EbDefinitions.h"
#include "EbObject.h"
#include "EbCabacContextModel.h"
#include "EbCodingUnit.h"
#include "EbSequenceControlSet.h"

typedef struct EbReferenceObject {
    EbDctor                     dctor;
    EbPictureBufferDesc        *reference_picture;
    EbPictureBufferDesc        *quarter_reference_picture;
    EbPictureBufferDesc        *sixteenth_reference_picture;
    EbDownScaledBufDescPtrArray ds_pics; // Pointer array for down scaled pictures
    EbPictureBufferDesc        *input_picture;
    EbPictureBufferDesc        *quarter_input_picture;
    EbPictureBufferDesc        *sixteenth_input_picture;
    EbPictureBufferDesc *downscaled_reference_picture[NUM_SR_SCALES + 1][NUM_RESIZE_SCALES + 1];
    uint64_t
        downscaled_picture_number[NUM_SR_SCALES + 1]
                                 [NUM_RESIZE_SCALES + 1]; // save the picture_number for each denom
    EbHandle  resize_mutex[NUM_SR_SCALES + 1][NUM_RESIZE_SCALES + 1];
    uint64_t  ref_poc;
    uint16_t  qp;
    SliceType slice_type;
    uint8_t   intra_coded_area; //percentage of intra coded area 0-100%
    uint8_t   skip_coded_area;

    uint8_t              tmp_layer_idx;
    Bool                 is_scene_change;
    uint16_t             pic_avg_variance;
    AomFilmGrain         film_grain_params; //Film grain parameters for a reference frame
    int8_t               sg_frame_ep;
    FRAME_CONTEXT        frame_context;
    EbWarpedMotionParams global_motion[TOTAL_REFS_PER_FRAME];
    MV_REF              *mvs;
    FrameType            frame_type;
    uint32_t             order_hint;
    uint32_t             ref_order_hint[7];
    double               r0;
    int32_t              filter_level[2];
    int32_t              filter_level_u;
    int32_t              filter_level_v;
    uint32_t             ref_cdef_strengths_num;
    uint8_t              ref_cdef_strengths[2][TOTAL_STRENGTHS];
    uint8_t             *sb_intra;
    uint8_t             *sb_skip;
    uint8_t             *sb_64x64_mvp;
    uint32_t            *sb_me_64x64_dist;
    uint32_t            *sb_me_8x8_cost_var;
    int32_t              mi_cols;
    int32_t              mi_rows;
    WienerUnitInfo *
        *unit_info; // per plane, per rest. unit; used for fwding wiener info to future frames
} EbReferenceObject;

typedef struct EbReferenceObjectDescInitData {
    EbPictureBufferDescInitData reference_picture_desc_init_data;
    int8_t                      hbd_mode_decision;
} EbReferenceObjectDescInitData;

typedef struct EbPaReferenceObject {
    EbDctor              dctor;
    EbPictureBufferDesc *input_padded_picture_ptr;
    EbPictureBufferDesc *quarter_downsampled_picture_ptr;
    EbPictureBufferDesc *sixteenth_downsampled_picture_ptr;
    // downscaled reference pointers
    // [super-res scales][resize scales]
    EbPictureBufferDesc
        *downscaled_input_padded_picture_ptr[NUM_SR_SCALES + 1][NUM_RESIZE_SCALES + 1];
    EbPictureBufferDesc
        *downscaled_quarter_downsampled_picture_ptr[NUM_SR_SCALES + 1][NUM_RESIZE_SCALES + 1];
    EbPictureBufferDesc
        *downscaled_sixteenth_downsampled_picture_ptr[NUM_SR_SCALES + 1][NUM_RESIZE_SCALES + 1];
    uint64_t
        downscaled_picture_number[NUM_SR_SCALES + 1]
                                 [NUM_RESIZE_SCALES + 1]; // save the picture_number for each denom
    EbHandle resize_mutex[NUM_SR_SCALES + 1][NUM_RESIZE_SCALES + 1];
    uint64_t picture_number;
    uint8_t  dummy_obj;
} EbPaReferenceObject;

typedef struct EbPaReferenceObjectDescInitData {
    EbPictureBufferDescInitData reference_picture_desc_init_data;
    EbPictureBufferDescInitData quarter_picture_desc_init_data;
    EbPictureBufferDescInitData sixteenth_picture_desc_init_data;
} EbPaReferenceObjectDescInitData;

/**************************************
 * Extern Function Declarations
 **************************************/
extern EbErrorType svt_reference_object_creator(EbPtr *object_dbl_ptr, EbPtr object_init_data_ptr);
extern EbErrorType svt_reference_object_reset(EbReferenceObject *obj, SequenceControlSet *scs_ptr);

extern EbErrorType svt_pa_reference_object_creator(EbPtr *object_dbl_ptr,
                                                   EbPtr  object_init_data_ptr);
void release_pa_reference_objects(SequenceControlSet *scs_ptr, PictureParentControlSet *pcs_ptr);

#endif //EbReferenceObject_h
