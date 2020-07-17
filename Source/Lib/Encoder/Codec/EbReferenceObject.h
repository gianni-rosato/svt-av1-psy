/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbReferenceObject_h
#define EbReferenceObject_h

#include "EbDefinitions.h"
#include "EbObject.h"
#include "EbCabacContextModel.h"
#include "EbCodingUnit.h"

typedef struct EbReferenceObject {
    EbDctor              dctor;
    EbPictureBufferDesc *reference_picture;
    EbPictureBufferDesc *reference_picture16bit;
    EbPictureBufferDesc *downscaled_reference_picture[NUM_SCALES];
    EbPictureBufferDesc *downscaled_reference_picture16bit[NUM_SCALES];
    uint64_t             ref_poc;
    uint16_t             qp;
    EB_SLICE             slice_type;
    uint8_t              intra_coded_area; //percentage of intra coded area 0-100%
    uint8_t              intra_coded_area_sb
        [MAX_NUMBER_OF_TREEBLOCKS_PER_PICTURE]; //percentage of intra coded area 0-100%
    uint32_t non_moving_index_array
        [MAX_NUMBER_OF_TREEBLOCKS_PER_PICTURE]; //array to hold non-moving blocks in reference frames
    uint8_t              tmp_layer_idx;
    EbBool               is_scene_change;
    uint16_t             pic_avg_variance;
    uint8_t              average_intensity;
    AomFilmGrain         film_grain_params; //Film grain parameters for a reference frame
    uint32_t             cdef_frame_strength;
    int8_t               sg_frame_ep;
    FRAME_CONTEXT        frame_context;
    EbWarpedMotionParams global_motion[TOTAL_REFS_PER_FRAME];
    MV_REF *             mvs;
    FrameType            frame_type;
    uint32_t             order_hint;
    uint32_t             ref_order_hint[7];
    StatStruct           stat_struct;
    EbHandle             referenced_area_mutex;
    uint64_t             referenced_area_avg;
#if TPL_1PASS_IMP
    double               r0;
#endif
#if ADAPTIVE_NSQ_CR
    uint32_t ref_part_cnt[NUMBER_OF_SHAPES-1][FB_NUM][SSEG_NUM];
#endif
#if ADAPTIVE_DEPTH_CR
#if SOFT_CYCLES_REDUCTION
    uint32_t ref_pred_depth_count[DEPTH_DELTA_NUM][NUMBER_OF_SHAPES-1];
#else
    uint32_t ref_pred_depth_count[DEPTH_DELTA_NUM];
#endif
#endif
#if ADAPTIVE_TXT_CR
    uint32_t ref_txt_cnt[TXT_DEPTH_DELTA_NUM][TX_TYPES];
#endif
    int32_t              mi_cols;
    int32_t              mi_rows;
} EbReferenceObject;

typedef struct EbReferenceObjectDescInitData {
    EbPictureBufferDescInitData reference_picture_desc_init_data;
#if MEM_OPT_10bit
#if CHANGE_HBD_MODE
    int8_t hbd_mode_decision;
#else
    uint8_t hbd_mode_decision;
#endif
#endif
} EbReferenceObjectDescInitData;

typedef struct EbPaReferenceObject {
    EbDctor              dctor;
    EbPictureBufferDesc *input_padded_picture_ptr;
    EbPictureBufferDesc *quarter_decimated_picture_ptr;
    EbPictureBufferDesc *sixteenth_decimated_picture_ptr;
    EbPictureBufferDesc *quarter_filtered_picture_ptr;
    EbPictureBufferDesc *sixteenth_filtered_picture_ptr;
    // downscaled reference pointers
    EbPictureBufferDesc *downscaled_input_padded_picture_ptr[NUM_SCALES];
    EbPictureBufferDesc *downscaled_quarter_decimated_picture_ptr[NUM_SCALES];
    EbPictureBufferDesc *downscaled_sixteenth_decimated_picture_ptr[NUM_SCALES];
    EbPictureBufferDesc *downscaled_quarter_filtered_picture_ptr[NUM_SCALES];
    EbPictureBufferDesc *downscaled_sixteenth_filtered_picture_ptr[NUM_SCALES];
    uint16_t             variance[MAX_NUMBER_OF_TREEBLOCKS_PER_PICTURE];
    uint8_t              y_mean[MAX_NUMBER_OF_TREEBLOCKS_PER_PICTURE];
    EB_SLICE             slice_type;
    uint32_t             dependent_pictures_count; //number of pic using this reference frame
} EbPaReferenceObject;

typedef struct EbPaReferenceObjectDescInitData {
    EbPictureBufferDescInitData reference_picture_desc_init_data;
    EbPictureBufferDescInitData quarter_picture_desc_init_data;
    EbPictureBufferDescInitData sixteenth_picture_desc_init_data;
} EbPaReferenceObjectDescInitData;

/**************************************
 * Extern Function Declarations
 **************************************/
extern EbErrorType eb_reference_object_creator(EbPtr *object_dbl_ptr, EbPtr object_init_data_ptr);

extern EbErrorType eb_pa_reference_object_creator(EbPtr *object_dbl_ptr,
                                                  EbPtr  object_init_data_ptr);

#endif //EbReferenceObject_h
