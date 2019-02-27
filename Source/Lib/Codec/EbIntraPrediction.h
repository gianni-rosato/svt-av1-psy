/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbIntraPrediction_h
#define EbIntraPrediction_h

#include "EbIntraPrediction_SSE2.h"
#include "EbIntraPrediction_SSSE3.h"
#include "EbIntraPrediction_SSE4_1.h"
#include "EbIntraPrediction_AVX2.h"

#include "EbDefinitions.h"
#include "EbUtility.h"
#include "EbPictureBufferDesc.h"
#include "EbPictureControlSet.h"
#include "EbCodingUnit.h"
#include "EbPictureControlSet.h"
#include "EbModeDecision.h"
#include "EbNeighborArrays.h"
#include "EbMotionEstimationProcess.h"

#ifdef __cplusplus
extern "C" {
#endif
#define MAX_PU_SIZE                            64

    struct ModeDecisionContext_s;

    typedef void(*intra_pred_fn_c)(uint8_t *dst, ptrdiff_t stride, int32_t w, int32_t h,
        const uint8_t *above, const uint8_t *left);
    typedef void(*intra_highbd_pred_fn_c)(uint16_t *dst, ptrdiff_t stride, int32_t w, int32_t h,
        const uint16_t *above, const uint16_t *left, int32_t bd);

    typedef void(*intra_pred_fn)(uint8_t *dst, ptrdiff_t stride,
        const uint8_t *above, const uint8_t *left);

    typedef void(*intra_high_pred_fn)(uint16_t *dst, ptrdiff_t stride,
        const uint16_t *above, const uint16_t *left,
        int32_t bd);

    typedef struct IntraReferenceSamples_s {

        uint8_t                  *y_intra_reference_array;
        uint8_t                  *cbIntraReferenceArray;
        uint8_t                  *crIntraReferenceArray;
        uint8_t                  *yIntraFilteredReferenceArray;

        uint8_t                  *y_intra_reference_array_reverse;
        uint8_t                  *yIntraFilteredReferenceArrayReverse;
        uint8_t                  *cbIntraReferenceArrayReverse;
        uint8_t                  *crIntraReferenceArrayReverse;

        // Scratch buffers used in the interpolaiton process
        uint8_t                   reference_above_line_y[(MAX_PU_SIZE << 2) + 1];
        uint8_t                   reference_left_line_y[(MAX_PU_SIZE << 2) + 1];
        EbBool                 above_ready_flag_y;
        EbBool                 left_ready_flag_y;

        uint8_t                   ReferenceAboveLineCb[(MAX_PU_SIZE << 2) + 2];
        uint8_t                   ReferenceLeftLineCb[(MAX_PU_SIZE << 2) + 2];
        EbBool                 AboveReadyFlagCb;
        EbBool                 LeftReadyFlagCb;

        uint8_t                   ReferenceAboveLineCr[(MAX_PU_SIZE << 2) + 2];
        uint8_t                   ReferenceLeftLineCr[(MAX_PU_SIZE << 2) + 2];
        EbBool                 AboveReadyFlagCr;
        EbBool                 LeftReadyFlagCr;

    } IntraReferenceSamples_t;

    typedef struct IntraReference16bitSamples_s {

        uint16_t                  *y_intra_reference_array;
        uint16_t                  *cbIntraReferenceArray;
        uint16_t                  *crIntraReferenceArray;
        uint16_t                  *yIntraFilteredReferenceArray;

        uint16_t                  *y_intra_reference_array_reverse;
        uint16_t                  *yIntraFilteredReferenceArrayReverse;
        uint16_t                  *cbIntraReferenceArrayReverse;
        uint16_t                  *crIntraReferenceArrayReverse;

        // Scratch buffers used in the interpolaiton process
        uint16_t                   reference_above_line_y[(MAX_PU_SIZE << 2) + 1];
        uint16_t                   reference_left_line_y[(MAX_PU_SIZE << 2) + 1];
        EbBool                  above_ready_flag_y;
        EbBool                  left_ready_flag_y;

        uint16_t                   ReferenceAboveLineCb[(MAX_PU_SIZE << 2) + 2];
        uint16_t                   ReferenceLeftLineCb[(MAX_PU_SIZE << 2) + 2];
        EbBool                  AboveReadyFlagCb;
        EbBool                  LeftReadyFlagCb;

        uint16_t                   ReferenceAboveLineCr[(MAX_PU_SIZE << 2) + 2];
        uint16_t                   ReferenceLeftLineCr[(MAX_PU_SIZE << 2) + 2];
        EbBool                  AboveReadyFlagCr;
        EbBool                  LeftReadyFlagCr;

    } IntraReference16bitSamples_t;

    extern EbErrorType IntraReferenceSamplesCtor(
        IntraReferenceSamples_t **context_dbl_ptr);



    extern EbErrorType IntraReference16bitSamplesCtor(
        IntraReference16bitSamples_t **context_dbl_ptr);



#define TOTAL_LUMA_MODES                   35
#define TOTAL_CHROMA_MODES                  5
#define TOTAL_INTRA_GROUPS                  5
#define INTRA_PLANAR_MODE                   0
#define INTRA_DC_MODE                       1
#define INTRA_HORIZONTAL_MODE              10
#define INTRA_VERTICAL_MODE                26
#define STRONG_INTRA_SMOOTHING_BLOCKSIZE   32
#define SMOOTHING_THRESHOLD                 8
#define SMOOTHING_THRESHOLD_10BIT          32

#if !QT_10BIT_SUPPORT

    extern EbErrorType GenerateIntraReferenceSamplesEncodePass(

        EbBool                         *is_left_availble,
        EbBool                         *is_above_availble,

        EbBool                     constrained_intra_flag,   //input parameter, indicates if constrained intra is switched on/off
        EbBool                     strongIntraSmoothingFlag,
        uint32_t                      origin_x,
        uint32_t                      origin_y,
        uint32_t                      size,
        uint32_t                      cu_depth,
        NeighborArrayUnit_t        *mode_type_neighbor_array,
        NeighborArrayUnit_t        *luma_recon_neighbor_array,
        NeighborArrayUnit_t        *cb_recon_neighbor_array,
        NeighborArrayUnit_t        *cr_recon_neighbor_array,
        void                       *refWrapperPtr,
        EbBool                     pictureLeftBoundary,
        EbBool                     pictureTopBoundary,
        EbBool                     pictureRightBoundary);
#endif




#if !QT_10BIT_SUPPORT

    extern EbErrorType GenerateIntraReference16bitSamplesEncodePass(
        EbBool                         *is_left_availble,
        EbBool                         *is_above_availble,
        EbBool                     constrained_intra_flag,   //input parameter, indicates if constrained intra is switched on/off
        EbBool                     strongIntraSmoothingFlag,
        uint32_t                      origin_x,
        uint32_t                      origin_y,
        uint32_t                      size,
        uint32_t                      cu_depth,
        NeighborArrayUnit_t        *mode_type_neighbor_array,
        NeighborArrayUnit_t        *luma_recon_neighbor_array,
        NeighborArrayUnit_t        *cb_recon_neighbor_array,
        NeighborArrayUnit_t        *cr_recon_neighbor_array,
        void                       *refWrapperPtr,
        EbBool                     pictureLeftBoundary,
        EbBool                     pictureTopBoundary,
        EbBool                     pictureRightBoundary);


    extern EbErrorType GenerateLumaIntraReference16bitSamplesEncodePass(
        EbBool                     *is_left_availble,
        EbBool                     *is_above_availble,
        EbBool                     constrained_intra_flag,   //input parameter, indicates if constrained intra is switched on/off
        EbBool                     strongIntraSmoothingFlag,
        uint32_t                      origin_x,
        uint32_t                      origin_y,
        uint32_t                      size,
        uint32_t                      sb_sz,
        uint32_t                      cu_depth,
        NeighborArrayUnit_t        *mode_type_neighbor_array,
        NeighborArrayUnit_t        *luma_recon_neighbor_array,
        NeighborArrayUnit_t        *cb_recon_neighbor_array,
        NeighborArrayUnit_t        *cr_recon_neighbor_array,
        void                       *refWrapperPtr,
        EbBool                     pictureLeftBoundary,
        EbBool                     pictureTopBoundary,
        EbBool                     pictureRightBoundary);


    extern EbErrorType GenerateChromaIntraReference16bitSamplesEncodePass(
        EbBool                     *is_left_availble,
        EbBool                     *is_above_availble,
        EbBool                     constrained_intra_flag,   //input parameter, indicates if constrained intra is switched on/off
        EbBool                     strongIntraSmoothingFlag,
        uint32_t                      origin_x,
        uint32_t                      origin_y,
        uint32_t                      size,
        uint32_t                      sb_sz,
        uint32_t                      cu_depth,
        NeighborArrayUnit_t        *mode_type_neighbor_array,
        NeighborArrayUnit_t        *luma_recon_neighbor_array,
        NeighborArrayUnit_t        *cb_recon_neighbor_array,
        NeighborArrayUnit_t        *cr_recon_neighbor_array,
        void                       *refWrapperPtr,
        EbBool                     pictureLeftBoundary,
        EbBool                     pictureTopBoundary,
        EbBool                     pictureRightBoundary);


    extern EbErrorType IntraPredictionCL(
        struct ModeDecisionContext_s           *context_ptr,
        uint32_t                                  component_mask,
        PictureControlSet_t                    *picture_control_set_ptr,
        ModeDecisionCandidateBuffer_t           *candidate_buffer_ptr,
        EbAsm                                  asm_type);
#endif

    extern EbErrorType AV1IntraPredictionCL(
        struct ModeDecisionContext_s           *context_ptr,
 #if !CHROMA_BLIND
        uint32_t                                  component_mask,
#endif
        PictureControlSet_t                    *picture_control_set_ptr,
        ModeDecisionCandidateBuffer_t           *candidate_buffer_ptr,
        EbAsm                                  asm_type);

#if !QT_10BIT_SUPPORT
    extern EbErrorType EncodePassIntraPrediction(

        uint8_t                         upsample_left,
        uint8_t                         upsample_above,
        uint8_t                          upsample_left_chroma,
        uint8_t                          upsample_above_chroma,

        EbBool                         is_left_availble,
        EbBool                         is_above_availble,
        void                                   *ref_samples,
        uint32_t                                  origin_x,
        uint32_t                                  origin_y,
        uint32_t                                  puSize,
        EbPictureBufferDesc_t                  *prediction_ptr,
        uint32_t                                  luma_mode,
        uint32_t                                  chroma_mode,
        int32_t                                  angle_delta,
        uint16_t                                  bitdepth,
        EbAsm                                  asm_type);
    extern EbErrorType EncodePassIntraPrediction16bit(

        uint8_t                         upsample_left,
        uint8_t                         upsample_above,
        uint8_t                          upsample_left_chroma,
        uint8_t                          upsample_above_chroma,

        EbBool                         is_left_availble,
        EbBool                         is_above_availble,
        void                                   *ref_samples,
        uint32_t                                  origin_x,
        uint32_t                                  origin_y,
        uint32_t                                  puSize,
        EbPictureBufferDesc_t                  *prediction_ptr,
        uint32_t                                  luma_mode,
        uint32_t                                  chroma_mode,
        int32_t                                  angle_delta,
        uint16_t                                  bitdepth,
        EbAsm                                  asm_type);

    extern EbErrorType EncodePassIntra4x4Prediction(
        uint8_t                         upsample_left,
        uint8_t                         upsample_above,
        uint8_t                          upsample_left_chroma,
        uint8_t                          upsample_above_chroma,

        EbBool                         is_left_availble,
        EbBool                         is_above_availble,
        IntraReferenceSamples_t                *referenceSamples,
        uint32_t                                  origin_x,
        uint32_t                                  origin_y,
        uint32_t                                  puSize,
        uint32_t                                  chromaPuSize,
        EbPictureBufferDesc_t                  *prediction_ptr,
        uint32_t                                  luma_mode,
        uint32_t                                  chroma_mode,
        uint32_t                                  component_mask,
        EbAsm                                  asm_type);

    extern EbErrorType EncodePassIntra4x4Prediction16bit(

        uint8_t                         upsample_left,
        uint8_t                         upsample_above,
        uint8_t                          upsample_left_chroma,
        uint8_t                          upsample_above_chroma,

        EbBool                         is_left_availble,
        EbBool                         is_above_availble,
        IntraReference16bitSamples_t           *referenceSamples,
        uint32_t                                  origin_x,
        uint32_t                                  origin_y,
        uint32_t                                  puSize,
        uint32_t                                  chromaPuSize,
        EbPictureBufferDesc_t                  *prediction_ptr,
        uint32_t                                  luma_mode,
        uint32_t                                  chroma_mode,
        uint32_t                                  component_mask,
        uint16_t                                  bitdepth,
        EbAsm                                  asm_type);

#endif

    static const uint32_t intraLumaModeNumber[] = {
        18,
        35,
        35,
        35,
        35 //4
    };

    static const uint32_t intraModeOrderMap[] = {
        34,    //  Planar
        9,     //  Vertical
        25,    //  Horizontal
        0,     //  DC
        1,     //  Intra mode 4
        5,     //  Intra mode 5
        13,    //  Intra mode 6
        17,    //  Intra mode 7
        21,    //  Intra mode 8
        29,    //  Intra mode 9
        33,    //  Intra mode 10
        3,     //  Intra mode 11
        7,     //  Intra mode 12
        11,    //  Intra mode 13
        15,    //  Intra mode 14
        19,    //  Intra mode 15
        23,    //  Intra mode 16
        27,    //  Intra mode 17
        31,    //  Intra mode 18
        2,     //  Intra mode 19
        4,     //  Intra mode 20
        6,     //  Intra mode 21
        8,     //  Intra mode 22
        10,    //  Intra mode 23
        12,    //  Intra mode 24
        14,    //  Intra mode 25
        16,    //  Intra mode 26
        18,    //  Intra mode 27
        20,    //  Intra mode 28
        22,    //  Intra mode 29
        24,    //  Intra mode 30
        26,    //  Intra mode 31
        28,    //  Intra mode 32
        30,    //  Intra mode 33
        32,    //  Intra mode 34
    };

    extern void intra_mode_angular_horizontal_kernel_ssse3_intrin(
        uint32_t            size,
        uint8_t            *ref_samp_main,
        uint8_t            *prediction_ptr,
        uint32_t            prediction_buffer_stride,
        const EbBool     skip,
        int32_t            intra_pred_angle);



    extern EbErrorType IntraOpenLoopReferenceSamplesCtor(
        IntraReferenceSamplesOpenLoop_t **context_dbl_ptr);
    extern void IntraOpenLoopReferenceSamplesDtor(
        IntraReferenceSamplesOpenLoop_t  *context_ptr);

    extern EbErrorType UpdateNeighborSamplesArrayOpenLoop(
        IntraReferenceSamplesOpenLoop_t *intra_ref_ptr,
        EbPictureBufferDesc_t           *input_ptr,
        uint32_t                           stride,
        uint32_t                           src_origin_x,
        uint32_t                           src_origin_y,
        uint32_t                           block_size);

    extern EbErrorType IntraPredictionOpenLoop(
        uint32_t                       cu_size,
        MotionEstimationContext_t   *context_ptr,
        uint32_t           openLoopIntraCandidate,
        EbAsm                       asm_type);

#if !QT_10BIT_SUPPORT
    extern EbErrorType Intra4x4IntraPredictionCL(
        uint32_t                                pu_index,
        uint32_t                                pu_origin_x,
        uint32_t                                pu_origin_y,
        uint32_t                                pu_width,
        uint32_t                                pu_height,
        uint32_t                                sb_sz,
        uint32_t                                component_mask,
        PictureControlSet_t                    *picture_control_set_ptr,
        ModeDecisionCandidateBuffer_t          *candidate_buffer_ptr,
        EbPtr                                  prediction_context_ptr,
        EbAsm                                   asm_type);
#endif

    /***************************************
    * Function Ptr Types
    ***************************************/
    typedef void(*EB_INTRA_NOANG_TYPE)(
        const uint32_t      size,
        uint8_t            *ref_samples,
        uint8_t            *prediction_ptr,
        const uint32_t      prediction_buffer_stride,
        const EbBool        skip);

    typedef void(*EB_INTRA_DC_AV1_TYPE)(
        EbBool        is_left_availble,
        EbBool        is_above_availble,
        const uint32_t   size,                       //input parameter, denotes the size of the current PU
        uint8_t         *ref_samples,                 //input parameter, pointer to the reference samples
        uint8_t         *dst,              //output parameter, pointer to the prediction
        const uint32_t   prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
        const EbBool  skip);                       //skip half rows

    typedef uint32_t(*EB_NEIGHBOR_DC_INTRA_TYPE)(
        MotionEstimationContext_t       *context_ptr,
        EbPictureBufferDesc_t           *input_ptr,
        uint32_t                           src_origin_x,
        uint32_t                           src_origin_y,
        uint32_t                           block_size,
        EbAsm                              asm_type);

    typedef void(*EB_INTRA_NOANG_16bit_TYPE)(
        const uint32_t   size,
        uint16_t         *ref_samples,
        uint16_t         *prediction_ptr,
        const uint32_t   prediction_buffer_stride,
        const EbBool  skip);


    typedef void(*EB_INTRA_ANG_Z1_Z2_Z3_16bit_TYPE)(
        const uint32_t   size,
        uint16_t         *ref_samples,
        uint16_t         *dst,
        const uint32_t   prediction_buffer_stride,
        const EbBool  skip,
        uint16_t          dx,
        uint16_t          dy,
        uint16_t          bd);


    typedef void(*EB_INTRA_ANG_TYPE)(
        uint32_t            size,
        uint8_t            *ref_samp_main,
        uint8_t            *prediction_ptr,
        uint32_t            prediction_buffer_stride,
        const EbBool     skip,
        int32_t            intra_pred_angle);

    typedef void(*EB_INTRA_ANG_16BIT_TYPE)(
        uint32_t          size,                       //input parameter, denotes the size of the current PU
        uint16_t         *ref_samp_main,                //input parameter, pointer to the reference samples
        uint16_t         *prediction_ptr,              //output parameter, pointer to the prediction
        uint32_t            prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
        const EbBool   skip,
        int32_t   intra_pred_angle);

    extern void IntraModePlanar(
        const uint32_t   size,                       //input parameter, denotes the size of the current PU
        uint8_t         *ref_samples,                 //input parameter, pointer to the reference samples
        uint8_t         *dst,              //output parameter, pointer to the prediction
        const uint32_t   prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
        const EbBool  skip);
#if !QT_10BIT_SUPPORT
    extern void highbd_smooth_v_predictor(
        const uint32_t   size,                       //input parameter, denotes the size of the current PU
        uint16_t         *ref_samples,                 //input parameter, pointer to the reference samples
        uint16_t         *dst,              //output parameter, pointer to the prediction
        const uint32_t   prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
        const EbBool  skip);
#endif
    extern void ebav1_smooth_v_predictor(
        const uint32_t   size,                       //input parameter, denotes the size of the current PU
        uint8_t         *ref_samples,                 //input parameter, pointer to the reference samples
        uint8_t         *dst,              //output parameter, pointer to the prediction
        const uint32_t   prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
        const EbBool  skip);
#if !QT_10BIT_SUPPORT
    extern void highbd_smooth_h_predictor(
        const uint32_t   size,                       //input parameter, denotes the size of the current PU
        uint16_t         *ref_samples,                 //input parameter, pointer to the reference samples
        uint16_t         *dst,              //output parameter, pointer to the prediction
        const uint32_t   prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
        const EbBool  skip);
#endif
    extern void ebav1_smooth_h_predictor(
        const uint32_t   size,                       //input parameter, denotes the size of the current PU
        uint8_t         *ref_samples,                 //input parameter, pointer to the reference samples
        uint8_t         *dst,              //output parameter, pointer to the prediction
        const uint32_t   prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
        const EbBool  skip);
#if !QT_10BIT_SUPPORT
    extern void highbd_dc_predictor(
        EbBool                         is_left_availble,
        EbBool                         is_above_availble,
        const uint32_t   size,                       //input parameter, denotes the size of the current PU
        uint8_t         *ref_samples,                 //input parameter, pointer to the reference samples
        uint8_t         *dst,              //output parameter, pointer to the prediction
        const uint32_t   prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
        const EbBool  skip);                     //skip half rows
#endif

    void IntraModeAngular_AV1_Z1_16bit(
        const uint32_t   size,                    //input parameter, denotes the size of the current PU
        uint16_t         *ref_samples,             //input parameter, pointer to the reference samples
        uint16_t         *dst,                    //output parameter, pointer to the prediction
        const uint32_t   prediction_buffer_stride,  //input parameter, denotes the stride for the prediction ptr
        const EbBool  skip,
        uint16_t          dx,                     //output parameter, pointer to the prediction
        uint16_t          dy,                      //output parameter, pointer to the prediction
        uint16_t          bd);

    void IntraModeAngular_AV1_Z2_16bit(
        const uint32_t   size,                       //input parameter, denotes the size of the current PU
        uint16_t         *ref_samples,                 //input parameter, pointer to the reference samples
        uint16_t         *dst,              //output parameter, pointer to the prediction
        const uint32_t   prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
        const EbBool  skip,
        uint16_t          dx,              //output parameter, pointer to the prediction
        uint16_t          dy,              //output parameter, pointer to the prediction
        uint16_t          bd);

    void IntraModeAngular_AV1_Z3_16bit(
        const uint32_t   size,                       //input parameter, denotes the size of the current PU
        uint16_t         *ref_samples,                 //input parameter, pointer to the reference samples
        uint16_t         *dst,              //output parameter, pointer to the prediction
        const uint32_t   prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
        const EbBool  skip,
        uint16_t          dx,              //output parameter, pointer to the prediction
        uint16_t          dy,              //output parameter, pointer to the prediction
        uint16_t          bd);


    /***************************************
    * Function Ptrs
    ***************************************/
    static EB_INTRA_NOANG_TYPE FUNC_TABLE IntraVerticalLuma_funcPtrArray[ASM_TYPE_TOTAL] = {
        // NON_AVX2
        intra_mode_vertical_luma_sse2_intrin,
        // AVX2
        intra_mode_vertical_luma_avx2_intrin,

    };


    static EB_INTRA_NOANG_TYPE FUNC_TABLE IntraVerticalChroma_funcPtrArray[ASM_TYPE_TOTAL] = {
        // NON_AVX2
        intra_mode_vertical_chroma_sse2_intrin,
        // AVX2
        intra_mode_vertical_chroma_sse2_intrin,
    };


    static EB_INTRA_NOANG_TYPE FUNC_TABLE IntraHorzLuma_funcPtrArray[ASM_TYPE_TOTAL] = {
        // NON_AVX2
        intra_mode_horizontal_luma_sse2_intrin,
        // AVX2
        intra_mode_horizontal_luma_sse2_intrin,
    };


    static EB_INTRA_NOANG_TYPE FUNC_TABLE IntraHorzChroma_funcPtrArray[ASM_TYPE_TOTAL] = {
        // NON_AVX2
        intra_mode_horizontal_chroma_sse2_intrin,
        // AVX2
        intra_mode_horizontal_chroma_sse2_intrin,
    };

#if !QT_10BIT_SUPPORT
    static EB_INTRA_DC_AV1_TYPE FUNC_TABLE IntraDC_Av1_funcPtrArray[9][ASM_TYPE_TOTAL] = {

        // 4x4
        {
            // NON_AVX2
            highbd_dc_predictor,
            // AVX2
            intra_mode_dc_4x4_av1_sse2_intrin,
        },
        // 8x8
        {
            // NON_AVX2
            highbd_dc_predictor,
            // AVX2
            intra_mode_dc_8x8_av1_sse2_intrin,
        },
        // 16x16
        {
            // NON_AVX2
            highbd_dc_predictor,
            // AVX2
            intra_mode_dc_16x16_av1_sse2_intrin,

        },
        // NxN
        {
            // NON_AVX2
            highbd_dc_predictor,
            // AVX2
            highbd_dc_predictor,

        },
        // 32x32
        {
            // NON_AVX2
            highbd_dc_predictor,
            // AVX2
            intra_mode_dc_32x32_av1_avx2_intrin,

        } ,
        // NxN
        {
            // NON_AVX2
            highbd_dc_predictor,
            // AVX2
            highbd_dc_predictor,

        },
        // NxN
        {
            // NON_AVX2
            highbd_dc_predictor,
            // AVX2
            highbd_dc_predictor,

        },
        // NxN
        {
            // NON_AVX2
            highbd_dc_predictor,
            // AVX2
            highbd_dc_predictor,

        },
        // 64x64
        {
            // NON_AVX2
            highbd_dc_predictor,
            // AVX2

            intra_mode_dc_64x64_av1_avx2_intrin,

        }

    };
#endif
    static EB_INTRA_NOANG_TYPE FUNC_TABLE IntraDCLuma_funcPtrArray[ASM_TYPE_TOTAL] = {
        // NON_AVX2
        intra_mode_dc_luma_sse2_intrin,
        // AVX2
        intra_mode_dc_luma_avx2_intrin,

    };

    uint32_t UpdateNeighborDcIntraPred(
        MotionEstimationContext_t       *context_ptr,
        EbPictureBufferDesc_t           *input_ptr,
        uint32_t                           src_origin_x,
        uint32_t                           src_origin_y,
        uint32_t                           block_size,
        EbAsm                             asm_type);

    static EB_INTRA_NOANG_16bit_TYPE FUNC_TABLE IntraDCLuma_16bit_funcPtrArray[ASM_TYPE_TOTAL] = {
        // NON_AVX2
        intra_mode_dc_luma16bit_sse4_1_intrin,
        // AVX2
        intra_mode_dc_luma16bit_sse4_1_intrin,
    };

    static EB_INTRA_NOANG_TYPE FUNC_TABLE IntraDCChroma_funcPtrArray[ASM_TYPE_TOTAL] = {
        // NON_AVX2
        intra_mode_dc_chroma_sse2_intrin,
        // AVX2
        intra_mode_dc_chroma_sse2_intrin,
    };


    static EB_INTRA_NOANG_TYPE FUNC_TABLE IntraPlanar_funcPtrArray[ASM_TYPE_TOTAL] = {
        // NON_AVX2
        intra_mode_planar_sse2_intrin,
        // AVX2
        intra_mode_planar_avx2_intrin,
    };

    void smooth_v_predictor_c(uint8_t *dst, ptrdiff_t stride, int32_t bw,
        int32_t bh, const uint8_t *above,
        const uint8_t *left);

    static EB_INTRA_NOANG_TYPE FUNC_TABLE IntraPlanar_Av1_funcPtrArray[ASM_TYPE_TOTAL] = {
        // NON_AVX2
        IntraModePlanar,
        // AVX2
        intra_mode_planar_av1_avx2_intrin,
    };
#if !QT_10BIT_SUPPORT
    static EB_INTRA_NOANG_16bit_TYPE FUNC_TABLE IntraSmoothV_16bit_Av1_funcPtrArray[ASM_TYPE_TOTAL] = {
        // NON_AVX2
        highbd_smooth_v_predictor,
        // AVX2
        highbd_smooth_v_predictor,
    };
#endif
    static EB_INTRA_NOANG_TYPE FUNC_TABLE IntraSmoothH_Av1_funcPtrArray[ASM_TYPE_TOTAL] = {
        // NON_AVX2
        ebav1_smooth_h_predictor,
        // AVX2
        ebav1_smooth_h_predictor,
    };
#if !QT_10BIT_SUPPORT
    static EB_INTRA_NOANG_16bit_TYPE FUNC_TABLE IntraSmoothH_16bit_Av1_funcPtrArray[ASM_TYPE_TOTAL] = {
        // NON_AVX2
        highbd_smooth_h_predictor,
        // AVX2
        highbd_smooth_h_predictor,
    };
#endif
    static EB_INTRA_NOANG_TYPE FUNC_TABLE IntraSmoothV_Av1_funcPtrArray[ASM_TYPE_TOTAL] = {
        // NON_AVX2
        ebav1_smooth_v_predictor,
        // AVX2
        ebav1_smooth_v_predictor,
    };

    static EB_INTRA_NOANG_16bit_TYPE FUNC_TABLE IntraPlanar_16bit_funcPtrArray[ASM_TYPE_TOTAL] = {
        // NON_AVX2
        intra_mode_planar16bit_sse2_intrin,
        // AVX2
        intra_mode_planar16bit_sse2_intrin,
    };

    static EB_INTRA_NOANG_TYPE FUNC_TABLE IntraAng34_funcPtrArray[ASM_TYPE_TOTAL] = {
        // NON_AVX2
        intra_mode_angular_34_sse2_intrin,
        // AVX2
        intra_mode_angular_34_avx2_intrin,
    };


    static EB_INTRA_NOANG_TYPE FUNC_TABLE IntraAng18_funcPtrArray[ASM_TYPE_TOTAL] = {
        // NON_AVX2
        intra_mode_angular_18_sse2_intrin,
        // AVX2
        intra_mode_angular_18_avx2_intrin,

    };


    static EB_INTRA_NOANG_TYPE FUNC_TABLE IntraAng2_funcPtrArray[ASM_TYPE_TOTAL] = {
        // NON_AVX2
        intra_mode_angular_2_sse2_intrin,
        // AVX2
        intra_mode_angular_2_avx2_intrin,
    };


    static EB_INTRA_ANG_TYPE FUNC_TABLE IntraAngVertical_funcPtrArray[ASM_TYPE_TOTAL] = {
        // NON_AVX2
        intra_mode_angular_vertical_kernel_ssse3_intrin,
        // AVX2
        intra_mode_angular_vertical_kernel_avx2_intrin,
    };


    static EB_INTRA_ANG_TYPE FUNC_TABLE IntraAngHorizontal_funcPtrArray[ASM_TYPE_TOTAL] = {
        // NON_AVX2
        intra_mode_angular_horizontal_kernel_ssse3_intrin,
        // AVX2
        intra_mode_angular_horizontal_kernel_avx2_intrin,
    };


    static EB_INTRA_ANG_Z1_Z2_Z3_16bit_TYPE FUNC_TABLE IntraModeAngular_AV1_Z1_16bit_funcPtrArray[9][ASM_TYPE_TOTAL] = {
        // 4x4
        {
            // NON_AVX2
            IntraModeAngular_AV1_Z1_16bit,
            // AVX2
            intra_mode_angular_av1_z1_16bit_4x4_avx2,
        },
        // 8x8
        {
            // NON_AVX2
            IntraModeAngular_AV1_Z1_16bit,
            // AVX2
            intra_mode_angular_av1_z1_16bit_8x8_avx2,
        },
        // 16x16
        {
            // NON_AVX2
            IntraModeAngular_AV1_Z1_16bit,
            // AVX2
            intra_mode_angular_av1_z1_16bit_16x16_avx2,
        },
        // NxN
        {
            // NON_AVX2
            IntraModeAngular_AV1_Z1_16bit,
            // AVX2
            IntraModeAngular_AV1_Z1_16bit,
        },
        // 32x32
        {
            // NON_AVX2
            IntraModeAngular_AV1_Z1_16bit,
            // AVX2
            intra_mode_angular_av1_z1_16bit_32x32_avx2,
        },
        // NxN
        {
            // NON_AVX2
            IntraModeAngular_AV1_Z1_16bit,
            // AVX2
            IntraModeAngular_AV1_Z1_16bit,
        },
        // NxN
        {
            // NON_AVX2
            IntraModeAngular_AV1_Z1_16bit,
            // AVX2
            IntraModeAngular_AV1_Z1_16bit,
        },
        // NxN
        {
            // NON_AVX2
            IntraModeAngular_AV1_Z1_16bit,
            // AVX2
            IntraModeAngular_AV1_Z1_16bit,
        },
        // 64x64
        {
            // NON_AVX2
            IntraModeAngular_AV1_Z1_16bit,
            // AVX2
            intra_mode_angular_av1_z1_16bit_64x64_avx2,
        }
    };
    static EB_INTRA_ANG_Z1_Z2_Z3_16bit_TYPE FUNC_TABLE IntraModeAngular_AV1_Z2_16bit_funcPtrArray[9][ASM_TYPE_TOTAL] = {
        // 4x4
        {
            // NON_AVX2
            IntraModeAngular_AV1_Z2_16bit,
            // AVX2
            intra_mode_angular_av1_z2_16bit_4x4_avx2,
        },
        // 8x8
        {
            // NON_AVX2
            IntraModeAngular_AV1_Z2_16bit,
            // AVX2
            intra_mode_angular_av1_z2_16bit_8x8_avx2,
        },
        // 16x16
        {
            // NON_AVX2
            IntraModeAngular_AV1_Z2_16bit,
            // AVX2
            intra_mode_angular_av1_z2_16bit_16x16_avx2,
        },
        // NxN
        {
            // NON_AVX2
            IntraModeAngular_AV1_Z2_16bit,
            // AVX2
            IntraModeAngular_AV1_Z2_16bit,
        },
        // 32x32
        {
            // NON_AVX2
            IntraModeAngular_AV1_Z2_16bit,
            // AVX2
            intra_mode_angular_av1_z2_16bit_32x32_avx2,
        },
        // NxN
        {
            // NON_AVX2
            IntraModeAngular_AV1_Z2_16bit,
            // AVX2
            IntraModeAngular_AV1_Z2_16bit,
        },
        // NxN
        {
            // NON_AVX2
            IntraModeAngular_AV1_Z2_16bit,
            // AVX2
            IntraModeAngular_AV1_Z2_16bit,
        },
        // NxN
        {
            // NON_AVX2
            IntraModeAngular_AV1_Z2_16bit,
            // AVX2
            IntraModeAngular_AV1_Z2_16bit,
        },
        // 64x64
        {
            // NON_AVX2
            IntraModeAngular_AV1_Z2_16bit,
            // AVX2
            intra_mode_angular_av1_z2_16bit_64x64_avx2,
        }
    };
    static EB_INTRA_ANG_Z1_Z2_Z3_16bit_TYPE FUNC_TABLE IntraModeAngular_AV1_Z3_16bit_funcPtrArray[9][ASM_TYPE_TOTAL] = {
        // 4x4
        {
            // NON_AVX2
            IntraModeAngular_AV1_Z3_16bit,
            // AVX2
            intra_mode_angular_av1_z3_16bit_4x4_avx2,
        },
        // 8x8
        {
            // NON_AVX2
            IntraModeAngular_AV1_Z3_16bit,
            // AVX2
            intra_mode_angular_av1_z3_16bit_8x8_avx2,
        },
        // 16x16
        {
            // NON_AVX2
            IntraModeAngular_AV1_Z3_16bit,
            // AVX2
            intra_mode_angular_av1_z3_16bit_16x16_avx2,
        },
        // NxN
        {
            // NON_AVX2
            IntraModeAngular_AV1_Z3_16bit,
            // AVX2
            IntraModeAngular_AV1_Z3_16bit,
        },
        // 32x32
        {
            // NON_AVX2
            IntraModeAngular_AV1_Z3_16bit,
            // AVX2
            intra_mode_angular_av1_z3_16bit_32x32_avx2,
        },
        // NxN
        {
            // NON_AVX2
            IntraModeAngular_AV1_Z3_16bit,
            // AVX2
            IntraModeAngular_AV1_Z3_16bit,
        },
        // NxN
        {
            // NON_AVX2
            IntraModeAngular_AV1_Z3_16bit,
            // AVX2
            IntraModeAngular_AV1_Z3_16bit,
        },
        // NxN
        {
            // NON_AVX2
            IntraModeAngular_AV1_Z3_16bit,
            // AVX2
            IntraModeAngular_AV1_Z3_16bit,
        },
        // 64x64
        {
            // NON_AVX2
            IntraModeAngular_AV1_Z3_16bit,
            // AVX2
            intra_mode_angular_av1_z3_16bit_64x64_avx2,
        }
    };

    extern void cfl_luma_subsampling_420_lbd_c(
        uint8_t *input, // AMIR-> Changed to 8 bit
        int32_t input_stride, int16_t *output_q3,
        int32_t width, int32_t height);
    extern void cfl_luma_subsampling_420_hbd_c(
        const uint16_t *input,
        int32_t input_stride, int16_t *output_q3,
        int32_t width, int32_t height);
    extern void subtract_average_c(
        int16_t *pred_buf_q3,
        int32_t width,
        int32_t height,
        int32_t round_offset,
        int32_t num_pel_log2);

    // Declare a size-specific wrapper for the size-generic function. The compiler
    // will inline the size generic function in here, the advantage is that the size
    // will be constant allowing for loop unrolling and other constant propagated
    // goodness.
    //#define CFL_SUB_AVG_X(arch, width, height, round_offset, num_pel_log2)
    //  void subtract_average_##width##x##height##_##arch(int16_t *pred_buf_q3) {
    //    subtract_average_##arch(pred_buf_q3, width, height, round_offset,
    //                            num_pel_log2);
    //  }
    //
    //// Declare size-specific wrappers for all valid CfL sizes.
    //#define CFL_SUB_AVG_FN(arch)
    //  CFL_SUB_AVG_X(arch, 4, 4, 8, 4)
    //  CFL_SUB_AVG_X(arch, 4, 8, 16, 5)
    //  CFL_SUB_AVG_X(arch, 4, 16, 32, 6)
    //  CFL_SUB_AVG_X(arch, 8, 4, 16, 5)
    //  CFL_SUB_AVG_X(arch, 8, 8, 32, 6)
    //  CFL_SUB_AVG_X(arch, 8, 16, 64, 7)
    //  CFL_SUB_AVG_X(arch, 8, 32, 128, 8)
    //  CFL_SUB_AVG_X(arch, 16, 4, 32, 6)
    //  CFL_SUB_AVG_X(arch, 16, 8, 64, 7)
    //  CFL_SUB_AVG_X(arch, 16, 16, 128, 8)
    //  CFL_SUB_AVG_X(arch, 16, 32, 256, 9)
    //  CFL_SUB_AVG_X(arch, 32, 8, 128, 8)
    //  CFL_SUB_AVG_X(arch, 32, 16, 256, 9)
    //  CFL_SUB_AVG_X(arch, 32, 32, 512, 10)
    //  cfl_subtract_average_fn get_subtract_average_fn_c(TxSize tx_size) {
    //    static const cfl_subtract_average_fn sub_avg[TX_SIZES_ALL] = {
    //      subtract_average_4x4_##arch,   /* 4x4 */
    //      subtract_average_8x8_##arch,   /* 8x8 */
    //      subtract_average_16x16_##arch, /* 16x16 */
    //      subtract_average_32x32_##arch, /* 32x32 */
    //      cfl_subtract_average_null,     /* 64x64 (invalid CFL size) */
    //      subtract_average_4x8_##arch,   /* 4x8 */
    //      subtract_average_8x4_##arch,   /* 8x4 */
    //      subtract_average_8x16_##arch,  /* 8x16 */
    //      subtract_average_16x8_##arch,  /* 16x8 */
    //      subtract_average_16x32_##arch, /* 16x32 */
    //      subtract_average_32x16_##arch, /* 32x16 */
    //      cfl_subtract_average_null,     /* 32x64 (invalid CFL size) */
    //      cfl_subtract_average_null,     /* 64x32 (invalid CFL size) */
    //      subtract_average_4x16_##arch,  /* 4x16 (invalid CFL size) */
    //      subtract_average_16x4_##arch,  /* 16x4 (invalid CFL size) */
    //      subtract_average_8x32_##arch,  /* 8x32 (invalid CFL size) */
    //      subtract_average_32x8_##arch,  /* 32x8 (invalid CFL size) */
    //      cfl_subtract_average_null,     /* 16x64 (invalid CFL size) */
    //      cfl_subtract_average_null,     /* 64x16 (invalid CFL size) */
    //    };
    //    /* Modulo TX_SIZES_ALL to ensure that an attacker won't be able to */
    //    /* index the function pointer array out of bounds. */
    //    return sub_avg[tx_size % TX_SIZES_ALL];
    //  }
    //
    //#define get_subtract_average_fn get_subtract_average_fn_c

    // Can we use CfL for the current block?
    //static INLINE CFL_ALLOWED_TYPE is_cfl_allowed(const MacroBlockD *xd) {
    //    const MbModeInfo *mbmi = xd->mi[0];
    //    const block_size bsize = mbmi->sb_type;
    //    assert(bsize < BlockSizeS_ALL);
    //    //if (0/*xd->lossless[mbmi->segment_id]*/) {
    //    //    // In lossless, CfL is available when the partition size is equal to the
    //    //    // transform size.
    //    //    const int32_t plane_bsize =
    //    //        get_plane_block_size(bsize, &xd->plane[AOM_PLANE_U]);
    //    //    return (CFL_ALLOWED_TYPE)(plane_bsize == BLOCK_4X4);
    //    //}
    //    // Spec: CfL is available to luma partitions lesser than or equal to 32x32
    //    return (CFL_ALLOWED_TYPE)(block_size_wide[bsize] <= 32 &&
    //        block_size_high[bsize] <= 32);
    //}
    /* Shift down with rounding for signed integers, for use when n >= 0 */

#define ROUND_POWER_OF_TWO_SIGNED(value, n)           \
  (((value) < 0) ? -ROUND_POWER_OF_TWO(-(value), (n)) \
                 : ROUND_POWER_OF_TWO((value), (n)))

    static INLINE int32_t get_scaled_luma_q0(int32_t alpha_q3, int16_t pred_buf_q3) {
        int32_t scaled_luma_q6 = alpha_q3 * pred_buf_q3;
        return ROUND_POWER_OF_TWO_SIGNED(scaled_luma_q6, 6);
    }

    //CFL_PREDICT_FN(c, lbd)



    void cfl_predict_lbd_c(
        const int16_t *pred_buf_q3,
        uint8_t *pred,// AMIR ADDED
        int32_t pred_stride,
        uint8_t *dst,// AMIR changed to 8 bit
        int32_t dst_stride,
        int32_t alpha_q3,
        int32_t bit_depth,
        int32_t width,
        int32_t height);

    void cfl_predict_hbd_c(
        const int16_t *pred_buf_q3,
        uint16_t *pred,// AMIR ADDED
        int32_t pred_stride,
        uint16_t *dst,// AMIR changed to 8 bit
        int32_t dst_stride,
        int32_t alpha_q3,
        int32_t bit_depth,
        int32_t width,
        int32_t height);


    static INLINE int32_t cfl_idx_to_alpha(int32_t alpha_idx, int32_t joint_sign,
        CFL_PRED_TYPE pred_type) {
        const int32_t alpha_sign = (pred_type == CFL_PRED_U) ? CFL_SIGN_U(joint_sign)
            : CFL_SIGN_V(joint_sign);
        if (alpha_sign == CFL_SIGN_ZERO) return 0;
        const int32_t abs_alpha_q3 =
            (pred_type == CFL_PRED_U) ? CFL_IDX_U(alpha_idx) : CFL_IDX_V(alpha_idx);
        return (alpha_sign == CFL_SIGN_POS) ? abs_alpha_q3 + 1 : -abs_alpha_q3 - 1;
    }





#ifdef __cplusplus
}
#endif
#endif // EbIntraPrediction_h
