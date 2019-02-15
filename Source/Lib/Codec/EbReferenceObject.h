/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbReferenceObject_h
#define EbReferenceObject_h

#include "EbDefinitions.h"
#include "EbDefinitions.h"
#include "EbAdaptiveMotionVectorPrediction.h"

typedef struct EbReferenceObject_s {
    EbPictureBufferDesc_t          *referencePicture;
    EbPictureBufferDesc_t          *referencePicture16bit;
    EbPictureBufferDesc_t          *refDenSrcPicture;

    TmvpUnit_t                     *tmvpMap;
    EbBool                          tmvpEnableFlag;
    uint64_t                        refPOC;
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
    EbBool                          penalizeSkipflag;
    uint8_t                         tmpLayerIdx;
    EbBool                          isSceneChange;
    uint16_t                        pic_avg_variance;
    uint8_t                         average_intensity;
    aom_film_grain_t                film_grain_params; //Film grain parameters for a reference frame
#if FAST_CDEF
    uint32_t                        cdef_frame_strength;
#endif
#if FAST_SG
    int8_t                          sg_frame_ep;
#endif
} EbReferenceObject_t;

typedef struct EbReferenceObjectDescInitData_s {
    EbPictureBufferDescInitData_t   referencePictureDescInitData;
} EbReferenceObjectDescInitData_t;

typedef struct EbPaReferenceObject_s {
    EbPictureBufferDesc_t          *inputPaddedPicturePtr;
    EbPictureBufferDesc_t          *quarterDecimatedPicturePtr;
    EbPictureBufferDesc_t          *sixteenthDecimatedPicturePtr;
    uint16_t                        variance[MAX_NUMBER_OF_TREEBLOCKS_PER_PICTURE];
    uint8_t                         yMean[MAX_NUMBER_OF_TREEBLOCKS_PER_PICTURE];
    EB_SLICE                        slice_type;
    uint32_t                        dependentPicturesCount; //number of pic using this reference frame
    PictureParentControlSet_t      *pPcsPtr;

} EbPaReferenceObject_t;

typedef struct EbPaReferenceObjectDescInitData_s {
    EbPictureBufferDescInitData_t   referencePictureDescInitData;
    EbPictureBufferDescInitData_t   quarterPictureDescInitData;
    EbPictureBufferDescInitData_t   sixteenthPictureDescInitData;
} EbPaReferenceObjectDescInitData_t;

/**************************************
 * Extern Function Declarations
 **************************************/
extern EbErrorType EbReferenceObjectCtor(
    EbPtr *object_dbl_ptr,
    EbPtr object_init_data_ptr);

extern EbErrorType EbPaReferenceObjectCtor(
    EbPtr *object_dbl_ptr,
    EbPtr object_init_data_ptr);


#endif //EbReferenceObject_h