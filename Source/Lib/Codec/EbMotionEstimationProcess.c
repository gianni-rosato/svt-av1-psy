/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>

#include "EbUtility.h"
#include "EbPictureControlSet.h"
#include "EbSequenceControlSet.h"
#include "EbPictureDecisionResults.h"
#include "EbMotionEstimationProcess.h"
#include "EbMotionEstimationResults.h"
#include "EbReferenceObject.h"
#include "EbMotionEstimation.h"
#include "EbIntraPrediction.h"
#include "EbLambdaRateTables.h"
#include "EbComputeSAD.h"

#include "emmintrin.h"

/* --32x32-
|00||01|
|02||03|
--------*/
/* ------16x16-----
|00||01||04||05|
|02||03||06||07|
|08||09||12||13|
|10||11||14||15|
----------------*/
/* ------8x8----------------------------
|00||01||04||05|     |16||17||20||21|
|02||03||06||07|     |18||19||22||23|
|08||09||12||13|     |24||25||28||29|
|10||11||14||15|     |26||27||30||31|

|32||33||36||37|     |48||49||52||53|
|34||35||38||39|     |50||51||54||55|
|40||41||44||45|     |56||57||60||61|
|42||43||46||47|     |58||59||62||63|
-------------------------------------*/
EbErrorType CheckZeroZeroCenter(
    PictureParentControlSet_t   *picture_control_set_ptr,
    EbPictureBufferDesc_t        *refPicPtr,
    MeContext_t                  *context_ptr,
    uint32_t                       sb_origin_x,
    uint32_t                       sb_origin_y,
    uint32_t                       sb_width,
    uint32_t                       sb_height,
    int16_t                       *xSearchCenter,
    int16_t                       *ySearchCenter,
    EbAsm                       asm_type);

#if ME_HME_OQ
/************************************************
 * Set ME/HME Params from Config
 ************************************************/
void* set_me_hme_params_from_config(
    SequenceControlSet_t        *sequence_control_set_ptr,
    MeContext_t                 *me_context_ptr)
{

    uint16_t hmeRegionIndex = 0;

    me_context_ptr->search_area_width = (uint8_t)sequence_control_set_ptr->static_config.search_area_width;
    me_context_ptr->search_area_height = (uint8_t)sequence_control_set_ptr->static_config.search_area_height;

    me_context_ptr->number_hme_search_region_in_width = (uint16_t)sequence_control_set_ptr->static_config.number_hme_search_region_in_width;
    me_context_ptr->number_hme_search_region_in_height = (uint16_t)sequence_control_set_ptr->static_config.number_hme_search_region_in_height;

    me_context_ptr->hme_level0_total_search_area_width = (uint16_t)sequence_control_set_ptr->static_config.hme_level0_total_search_area_width;
    me_context_ptr->hme_level0_total_search_area_height = (uint16_t)sequence_control_set_ptr->static_config.hme_level0_total_search_area_height;

    for (hmeRegionIndex = 0; hmeRegionIndex < me_context_ptr->number_hme_search_region_in_width; ++hmeRegionIndex) {
        me_context_ptr->hme_level0_search_area_in_width_array[hmeRegionIndex] = (uint16_t)sequence_control_set_ptr->static_config.hme_level0_search_area_in_width_array[hmeRegionIndex];
        me_context_ptr->hme_level1_search_area_in_width_array[hmeRegionIndex] = (uint16_t)sequence_control_set_ptr->static_config.hme_level1_search_area_in_width_array[hmeRegionIndex];
        me_context_ptr->hme_level2_search_area_in_width_array[hmeRegionIndex] = (uint16_t)sequence_control_set_ptr->static_config.hme_level2_search_area_in_width_array[hmeRegionIndex];
    }

    for (hmeRegionIndex = 0; hmeRegionIndex < me_context_ptr->number_hme_search_region_in_height; ++hmeRegionIndex) {
        me_context_ptr->hme_level0_search_area_in_height_array[hmeRegionIndex] = (uint16_t)sequence_control_set_ptr->static_config.hme_level0_search_area_in_height_array[hmeRegionIndex];
        me_context_ptr->hme_level1_search_area_in_height_array[hmeRegionIndex] = (uint16_t)sequence_control_set_ptr->static_config.hme_level1_search_area_in_height_array[hmeRegionIndex];
        me_context_ptr->hme_level2_search_area_in_height_array[hmeRegionIndex] = (uint16_t)sequence_control_set_ptr->static_config.hme_level2_search_area_in_height_array[hmeRegionIndex];
    }

    return EB_NULL;
}


/************************************************
 * Set ME/HME Params
 ************************************************/
void* set_me_hme_params_oq(
    MeContext_t                     *me_context_ptr,
    PictureParentControlSet_t       *picture_control_set_ptr,
    SequenceControlSet_t            *sequence_control_set_ptr,
    EbInputResolution                 input_resolution)
{
    (void)*picture_control_set_ptr;
    //uint8_t  hmeMeLevel = picture_control_set_ptr->enc_mode;
    uint8_t  hmeMeLevel =  picture_control_set_ptr->enc_mode; // OMK to be revised after new presets


// HME/ME default settings
    me_context_ptr->number_hme_search_region_in_width = 2;
    me_context_ptr->number_hme_search_region_in_height = 2;

    // HME Level0
    me_context_ptr->hme_level0_total_search_area_width = HmeLevel0TotalSearchAreaWidth[input_resolution][hmeMeLevel];
    me_context_ptr->hme_level0_total_search_area_height = HmeLevel0TotalSearchAreaHeight[input_resolution][hmeMeLevel];
    me_context_ptr->hme_level0_search_area_in_width_array[0] = HmeLevel0SearchAreaInWidthArrayRight[input_resolution][hmeMeLevel];
    me_context_ptr->hme_level0_search_area_in_width_array[1] = HmeLevel0SearchAreaInWidthArrayLeft[input_resolution][hmeMeLevel];
    me_context_ptr->hme_level0_search_area_in_height_array[0] = HmeLevel0SearchAreaInHeightArrayTop[input_resolution][hmeMeLevel];
    me_context_ptr->hme_level0_search_area_in_height_array[1] = HmeLevel0SearchAreaInHeightArrayBottom[input_resolution][hmeMeLevel];
    // HME Level1
    me_context_ptr->hme_level1_search_area_in_width_array[0] = HmeLevel1SearchAreaInWidthArrayRight[input_resolution][hmeMeLevel];
    me_context_ptr->hme_level1_search_area_in_width_array[1] = HmeLevel1SearchAreaInWidthArrayLeft[input_resolution][hmeMeLevel];
    me_context_ptr->hme_level1_search_area_in_height_array[0] = HmeLevel1SearchAreaInHeightArrayTop[input_resolution][hmeMeLevel];
    me_context_ptr->hme_level1_search_area_in_height_array[1] = HmeLevel1SearchAreaInHeightArrayBottom[input_resolution][hmeMeLevel];
    // HME Level2
    me_context_ptr->hme_level2_search_area_in_width_array[0] = HmeLevel2SearchAreaInWidthArrayRight[input_resolution][hmeMeLevel];
    me_context_ptr->hme_level2_search_area_in_width_array[1] = HmeLevel2SearchAreaInWidthArrayLeft[input_resolution][hmeMeLevel];
    me_context_ptr->hme_level2_search_area_in_height_array[0] = HmeLevel2SearchAreaInHeightArrayTop[input_resolution][hmeMeLevel];
    me_context_ptr->hme_level2_search_area_in_height_array[1] = HmeLevel2SearchAreaInHeightArrayBottom[input_resolution][hmeMeLevel];

    // ME
    me_context_ptr->search_area_width = SearchAreaWidth[input_resolution][hmeMeLevel];
    me_context_ptr->search_area_height = SearchAreaHeight[input_resolution][hmeMeLevel];


    // HME Level0 adjustment for low frame rate contents (frame rate <= 30)
    if (input_resolution == INPUT_SIZE_4K_RANGE) {
        if ((sequence_control_set_ptr->static_config.frame_rate >> 16) <= 30) {

            if (hmeMeLevel == ENC_M6) {
                me_context_ptr->hme_level0_total_search_area_width = MAX(96, me_context_ptr->hme_level0_total_search_area_width);
                me_context_ptr->hme_level0_total_search_area_height = MAX(64, me_context_ptr->hme_level0_total_search_area_height);
                me_context_ptr->hme_level0_search_area_in_width_array[0] = MAX(48, me_context_ptr->hme_level0_search_area_in_width_array[0]);
                me_context_ptr->hme_level0_search_area_in_width_array[1] = MAX(48, me_context_ptr->hme_level0_search_area_in_width_array[1]);
                me_context_ptr->hme_level0_search_area_in_height_array[0] = MAX(32, me_context_ptr->hme_level0_search_area_in_height_array[0]);
                me_context_ptr->hme_level0_search_area_in_height_array[1] = MAX(32, me_context_ptr->hme_level0_search_area_in_height_array[1]);

            }
        }
    }

    me_context_ptr->update_hme_search_center_flag = 1;

    if (input_resolution <= INPUT_SIZE_576p_RANGE_OR_LOWER) 
        me_context_ptr->update_hme_search_center_flag = 0;
    
    return EB_NULL;
};
/******************************************************
* Derive ME Settings for OQ
  Input   : encoder mode and tune
  Output  : ME Kernel signal(s)
******************************************************/
EbErrorType signal_derivation_me_kernel_oq(
    SequenceControlSet_t        *sequence_control_set_ptr,
    PictureParentControlSet_t   *picture_control_set_ptr,
    MotionEstimationContext_t   *context_ptr) {

    EbErrorType return_error = EB_ErrorNone;

    // Set ME/HME search regions
    if (sequence_control_set_ptr->static_config.use_default_me_hme) {
        set_me_hme_params_oq(
            context_ptr->me_context_ptr,
            picture_control_set_ptr,
            sequence_control_set_ptr,
            sequence_control_set_ptr->input_resolution);
    }
    else {
        set_me_hme_params_from_config(
            sequence_control_set_ptr,
            context_ptr->me_context_ptr);
    }

    return return_error;
};
#endif
/************************************************
 * Motion Analysis Context Constructor
 ************************************************/

EbErrorType MotionEstimationContextCtor(
    MotionEstimationContext_t   **context_dbl_ptr,
    EbFifo_t                     *pictureDecisionResultsInputFifoPtr,
    EbFifo_t                     *motionEstimationResultsOutputFifoPtr) {

    EbErrorType return_error = EB_ErrorNone;
    MotionEstimationContext_t *context_ptr;
    EB_MALLOC(MotionEstimationContext_t*, context_ptr, sizeof(MotionEstimationContext_t), EB_N_PTR);

    *context_dbl_ptr = context_ptr;

    context_ptr->pictureDecisionResultsInputFifoPtr = pictureDecisionResultsInputFifoPtr;
    context_ptr->motionEstimationResultsOutputFifoPtr = motionEstimationResultsOutputFifoPtr;
    return_error = IntraOpenLoopReferenceSamplesCtor(&context_ptr->intra_ref_ptr);
    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }
    return_error = MeContextCtor(&(context_ptr->me_context_ptr));
    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }

    return EB_ErrorNone;

}

/***************************************************************************************************
* ZZ Decimated SAD Computation
***************************************************************************************************/
EbErrorType ComputeDecimatedZzSad(
    MotionEstimationContext_t   *context_ptr,
    SequenceControlSet_t        *sequence_control_set_ptr,
    PictureParentControlSet_t   *picture_control_set_ptr,
    EbPictureBufferDesc_t       *sixteenthDecimatedPicturePtr,
    uint32_t                         xLcuStartIndex,
    uint32_t                         xLcuEndIndex,
    uint32_t                         yLcuStartIndex,
    uint32_t                         yLcuEndIndex) {

    EbErrorType return_error = EB_ErrorNone;

    EbAsm asm_type = sequence_control_set_ptr->encode_context_ptr->asm_type;

    PictureParentControlSet_t    *previous_picture_control_set_wrapper_ptr = ((PictureParentControlSet_t*)picture_control_set_ptr->previous_picture_control_set_wrapper_ptr->objectPtr);
    EbPictureBufferDesc_t        *previousInputPictureFull = previous_picture_control_set_wrapper_ptr->enhanced_picture_ptr;

    uint32_t sb_index;

    uint32_t sb_width;
    uint32_t sb_height;

    uint32_t decimatedLcuWidth;
    uint32_t decimatedLcuHeight;

    uint32_t sb_origin_x;
    uint32_t sb_origin_y;

    uint32_t blkDisplacementDecimated;
    uint32_t blkDisplacementFull;

    uint32_t decimatedLcuCollocatedSad;

    uint32_t xLcuIndex;
    uint32_t yLcuIndex;

    for (yLcuIndex = yLcuStartIndex; yLcuIndex < yLcuEndIndex; ++yLcuIndex) {
        for (xLcuIndex = xLcuStartIndex; xLcuIndex < xLcuEndIndex; ++xLcuIndex) {

            sb_index = xLcuIndex + yLcuIndex * sequence_control_set_ptr->picture_width_in_sb;
            SbParams_t *sb_params = &sequence_control_set_ptr->sb_params_array[sb_index];

            sb_width = sb_params->width;
            sb_height = sb_params->height;

            sb_origin_x = sb_params->origin_x;
            sb_origin_y = sb_params->origin_y;

            sb_width = sb_params->width;
            sb_height = sb_params->height;


            decimatedLcuWidth = sb_width >> 2;
            decimatedLcuHeight = sb_height >> 2;

            decimatedLcuCollocatedSad = 0;

            if (sb_params->is_complete_sb)
            {

                blkDisplacementDecimated = (sixteenthDecimatedPicturePtr->origin_y + (sb_origin_y >> 2)) * sixteenthDecimatedPicturePtr->strideY + sixteenthDecimatedPicturePtr->origin_x + (sb_origin_x >> 2);
                blkDisplacementFull = (previousInputPictureFull->origin_y + sb_origin_y)* previousInputPictureFull->strideY + (previousInputPictureFull->origin_x + sb_origin_x);


                // 1/16 collocated SB decimation
                Decimation2D(
                    &previousInputPictureFull->bufferY[blkDisplacementFull],
                    previousInputPictureFull->strideY,
                    BLOCK_SIZE_64,
                    BLOCK_SIZE_64,
                    context_ptr->me_context_ptr->sixteenth_sb_buffer,
                    context_ptr->me_context_ptr->sixteenth_sb_buffer_stride,
                    4);

                if (asm_type >= ASM_NON_AVX2 && asm_type < ASM_TYPE_TOTAL)
                    // ZZ SAD between 1/16 current & 1/16 collocated
                    decimatedLcuCollocatedSad = NxMSadKernel_funcPtrArray[asm_type][2](
                        &(sixteenthDecimatedPicturePtr->bufferY[blkDisplacementDecimated]),
                        sixteenthDecimatedPicturePtr->strideY,
                        context_ptr->me_context_ptr->sixteenth_sb_buffer,
                        context_ptr->me_context_ptr->sixteenth_sb_buffer_stride,
                        16, 16);

                // Background Enhancement Algorithm
                // Classification is important to:
                // 1. Avoid improving moving objects.
                // 2. Do not modulate when all the picture is background
                // 3. Do give different importance to different regions
                if (decimatedLcuCollocatedSad < BEA_CLASS_0_0_DEC_TH) {
                    previous_picture_control_set_wrapper_ptr->zz_cost_array[sb_index] = BEA_CLASS_0_ZZ_COST;
                }
                else if (decimatedLcuCollocatedSad < BEA_CLASS_0_DEC_TH) {
                    previous_picture_control_set_wrapper_ptr->zz_cost_array[sb_index] = BEA_CLASS_0_1_ZZ_COST;
                }
                else if (decimatedLcuCollocatedSad < BEA_CLASS_1_DEC_TH) {
                    previous_picture_control_set_wrapper_ptr->zz_cost_array[sb_index] = BEA_CLASS_1_ZZ_COST;
                }
                else if (decimatedLcuCollocatedSad < BEA_CLASS_2_DEC_TH) {
                    previous_picture_control_set_wrapper_ptr->zz_cost_array[sb_index] = BEA_CLASS_2_ZZ_COST;
                }
                else {
                    previous_picture_control_set_wrapper_ptr->zz_cost_array[sb_index] = BEA_CLASS_3_ZZ_COST;
                }


            }
            else {
                previous_picture_control_set_wrapper_ptr->zz_cost_array[sb_index] = INVALID_ZZ_COST;
                decimatedLcuCollocatedSad = (uint32_t)~0;
            }


            // Keep track of non moving LCUs for QP modulation
            if (decimatedLcuCollocatedSad < ((decimatedLcuWidth * decimatedLcuHeight) * 2)) {
                previous_picture_control_set_wrapper_ptr->non_moving_index_array[sb_index] = BEA_CLASS_0_ZZ_COST;
            }
            else if (decimatedLcuCollocatedSad < ((decimatedLcuWidth * decimatedLcuHeight) * 4)) {
                previous_picture_control_set_wrapper_ptr->non_moving_index_array[sb_index] = BEA_CLASS_1_ZZ_COST;
            }
            else if (decimatedLcuCollocatedSad < ((decimatedLcuWidth * decimatedLcuHeight) * 8)) {
                previous_picture_control_set_wrapper_ptr->non_moving_index_array[sb_index] = BEA_CLASS_2_ZZ_COST;
            }
            else { //if (decimatedLcuCollocatedSad < ((decimatedLcuWidth * decimatedLcuHeight) * 4)) {
                previous_picture_control_set_wrapper_ptr->non_moving_index_array[sb_index] = BEA_CLASS_3_ZZ_COST;
            }
        }
    }

    return return_error;
}




/************************************************
 * Motion Analysis Kernel
 * The Motion Analysis performs  Motion Estimation
 * This process has access to the current input picture as well as
 * the input pictures, which the current picture references according
 * to the prediction structure pattern.  The Motion Analysis process is multithreaded,
 * so pictures can be processed out of order as long as all inputs are available.
 ************************************************/
void* MotionEstimationKernel(void *input_ptr)
{
    MotionEstimationContext_t   *context_ptr = (MotionEstimationContext_t*)input_ptr;

    PictureParentControlSet_t   *picture_control_set_ptr;
    SequenceControlSet_t        *sequence_control_set_ptr;

    EbObjectWrapper_t           *inputResultsWrapperPtr;
    PictureDecisionResults_t    *inputResultsPtr;

    EbObjectWrapper_t           *outputResultsWrapperPtr;
    MotionEstimationResults_t   *outputResultsPtr;

    EbPictureBufferDesc_t       *inputPicturePtr;

    EbPictureBufferDesc_t       *inputPaddedPicturePtr;

    uint32_t                       bufferIndex;

    uint32_t                       sb_index;
    uint32_t                       xLcuIndex;
    uint32_t                       yLcuIndex;
    uint32_t                       picture_width_in_sb;
    uint32_t                       picture_height_in_sb;
    uint32_t                       sb_origin_x;
    uint32_t                       sb_origin_y;
    uint32_t                       sb_width;
    uint32_t                       sb_height;
    uint32_t                       lcuRow;



    EbPaReferenceObject_t       *paReferenceObject;
    EbPictureBufferDesc_t       *quarterDecimatedPicturePtr;
    EbPictureBufferDesc_t       *sixteenthDecimatedPicturePtr;

    // Segments
    uint32_t                      segment_index;
    uint32_t                      xSegmentIndex;
    uint32_t                      ySegmentIndex;
    uint32_t                      xLcuStartIndex;
    uint32_t                      xLcuEndIndex;
    uint32_t                      yLcuStartIndex;
    uint32_t                      yLcuEndIndex;

    uint32_t                      intra_sad_interval_index;

    EbAsm                      asm_type;
    MdRateEstimationContext_t   *md_rate_estimation_array;


    for (;;) {


        // Get Input Full Object
        EbGetFullObject(
            context_ptr->pictureDecisionResultsInputFifoPtr,
            &inputResultsWrapperPtr);

        inputResultsPtr = (PictureDecisionResults_t*)inputResultsWrapperPtr->objectPtr;
        picture_control_set_ptr = (PictureParentControlSet_t*)inputResultsPtr->pictureControlSetWrapperPtr->objectPtr;
        sequence_control_set_ptr = (SequenceControlSet_t*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->objectPtr;
        paReferenceObject = (EbPaReferenceObject_t*)picture_control_set_ptr->pa_reference_picture_wrapper_ptr->objectPtr;
        quarterDecimatedPicturePtr = (EbPictureBufferDesc_t*)paReferenceObject->quarterDecimatedPicturePtr;
        sixteenthDecimatedPicturePtr = (EbPictureBufferDesc_t*)paReferenceObject->sixteenthDecimatedPicturePtr;

        inputPaddedPicturePtr = (EbPictureBufferDesc_t*)paReferenceObject->inputPaddedPicturePtr;

        inputPicturePtr = picture_control_set_ptr->enhanced_picture_ptr;

        // Segments
        segment_index = inputResultsPtr->segment_index;
        picture_width_in_sb = (sequence_control_set_ptr->luma_width + sequence_control_set_ptr->sb_sz - 1) / sequence_control_set_ptr->sb_sz;
        picture_height_in_sb = (sequence_control_set_ptr->luma_height + sequence_control_set_ptr->sb_sz - 1) / sequence_control_set_ptr->sb_sz;
        SEGMENT_CONVERT_IDX_TO_XY(segment_index, xSegmentIndex, ySegmentIndex, picture_control_set_ptr->me_segments_column_count);
        xLcuStartIndex = SEGMENT_START_IDX(xSegmentIndex, picture_width_in_sb, picture_control_set_ptr->me_segments_column_count);
        xLcuEndIndex = SEGMENT_END_IDX(xSegmentIndex, picture_width_in_sb, picture_control_set_ptr->me_segments_column_count);
        yLcuStartIndex = SEGMENT_START_IDX(ySegmentIndex, picture_height_in_sb, picture_control_set_ptr->me_segments_row_count);
        yLcuEndIndex = SEGMENT_END_IDX(ySegmentIndex, picture_height_in_sb, picture_control_set_ptr->me_segments_row_count);
        asm_type = sequence_control_set_ptr->encode_context_ptr->asm_type;
        // Increment the MD Rate Estimation array pointer to point to the right address based on the QP and slice type
        md_rate_estimation_array = (MdRateEstimationContext_t*)sequence_control_set_ptr->encode_context_ptr->md_rate_estimation_array;
        md_rate_estimation_array += picture_control_set_ptr->slice_type * TOTAL_NUMBER_OF_QP_VALUES + picture_control_set_ptr->picture_qp;
        // Reset MD rate Estimation table to initial values by copying from md_rate_estimation_array
        EB_MEMCPY(&(context_ptr->me_context_ptr->mvd_bits_array[0]), &(md_rate_estimation_array->mvdBits[0]), sizeof(EB_BitFraction)*NUMBER_OF_MVD_CASES);
        ///context_ptr->me_context_ptr->lambda = lambdaModeDecisionLdSadQpScaling[picture_control_set_ptr->picture_qp];
#if ME_HME_OQ
   // ME Kernel Signal(s) derivation

        signal_derivation_me_kernel_oq(
            sequence_control_set_ptr,
            picture_control_set_ptr,
            context_ptr);
#endif
        // Lambda Assignement
        if (sequence_control_set_ptr->static_config.pred_structure == EB_PRED_RANDOM_ACCESS) {

            if (picture_control_set_ptr->temporal_layer_index == 0) {
                context_ptr->me_context_ptr->lambda = lambdaModeDecisionRaSad[picture_control_set_ptr->picture_qp];
            }
            else if (picture_control_set_ptr->temporal_layer_index < 3) {
                context_ptr->me_context_ptr->lambda = lambdaModeDecisionRaSadQpScalingL1[picture_control_set_ptr->picture_qp];
            }
            else {
                context_ptr->me_context_ptr->lambda = lambdaModeDecisionRaSadQpScalingL3[picture_control_set_ptr->picture_qp];
            }
        }
        else {
            if (picture_control_set_ptr->temporal_layer_index == 0) {
                context_ptr->me_context_ptr->lambda = lambdaModeDecisionLdSad[picture_control_set_ptr->picture_qp];
            }
            else {
                context_ptr->me_context_ptr->lambda = lambdaModeDecisionLdSadQpScaling[picture_control_set_ptr->picture_qp];
            }
        }

        // *** MOTION ESTIMATION CODE ***
        if (picture_control_set_ptr->slice_type != I_SLICE) {

            // SB Loop
            for (yLcuIndex = yLcuStartIndex; yLcuIndex < yLcuEndIndex; ++yLcuIndex) {
                for (xLcuIndex = xLcuStartIndex; xLcuIndex < xLcuEndIndex; ++xLcuIndex) {

                    sb_index = (uint16_t)(xLcuIndex + yLcuIndex * picture_width_in_sb);
                    sb_origin_x = xLcuIndex * sequence_control_set_ptr->sb_sz;
                    sb_origin_y = yLcuIndex * sequence_control_set_ptr->sb_sz;

                    sb_width = (sequence_control_set_ptr->luma_width - sb_origin_x) < BLOCK_SIZE_64 ? sequence_control_set_ptr->luma_width - sb_origin_x : BLOCK_SIZE_64;
                    sb_height = (sequence_control_set_ptr->luma_height - sb_origin_y) < BLOCK_SIZE_64 ? sequence_control_set_ptr->luma_height - sb_origin_y : BLOCK_SIZE_64;

                    // Load the SB from the input to the intermediate SB buffer
                    bufferIndex = (inputPicturePtr->origin_y + sb_origin_y) * inputPicturePtr->strideY + inputPicturePtr->origin_x + sb_origin_x;

                    context_ptr->me_context_ptr->hme_search_type = HME_RECTANGULAR;

                    for (lcuRow = 0; lcuRow < BLOCK_SIZE_64; lcuRow++) {
                        EB_MEMCPY((&(context_ptr->me_context_ptr->sb_buffer[lcuRow * BLOCK_SIZE_64])), (&(inputPicturePtr->bufferY[bufferIndex + lcuRow * inputPicturePtr->strideY])), BLOCK_SIZE_64 * sizeof(uint8_t));

                    }

                    {
                        uint8_t * srcPtr = &inputPaddedPicturePtr->bufferY[bufferIndex];

                        //_MM_HINT_T0     //_MM_HINT_T1    //_MM_HINT_T2//_MM_HINT_NTA
                        uint32_t i;
                        for (i = 0; i < sb_height; i++)
                        {
                            char const* p = (char const*)(srcPtr + i * inputPaddedPicturePtr->strideY);
                            _mm_prefetch(p, _MM_HINT_T2);
                        }
                    }


                    context_ptr->me_context_ptr->sb_src_ptr = &inputPaddedPicturePtr->bufferY[bufferIndex];
                    context_ptr->me_context_ptr->sb_src_stride = inputPaddedPicturePtr->strideY;


                    // Load the 1/4 decimated SB from the 1/4 decimated input to the 1/4 intermediate SB buffer
                    if (picture_control_set_ptr->enable_hme_level1_flag) {

                        bufferIndex = (quarterDecimatedPicturePtr->origin_y + (sb_origin_y >> 1)) * quarterDecimatedPicturePtr->strideY + quarterDecimatedPicturePtr->origin_x + (sb_origin_x >> 1);

                        for (lcuRow = 0; lcuRow < (sb_height >> 1); lcuRow++) {
                            EB_MEMCPY((&(context_ptr->me_context_ptr->quarter_sb_buffer[lcuRow * context_ptr->me_context_ptr->quarter_sb_buffer_stride])), (&(quarterDecimatedPicturePtr->bufferY[bufferIndex + lcuRow * quarterDecimatedPicturePtr->strideY])), (sb_width >> 1) * sizeof(uint8_t));

                        }
                    }

                    // Load the 1/16 decimated SB from the 1/16 decimated input to the 1/16 intermediate SB buffer
                    if (picture_control_set_ptr->enable_hme_level0_flag) {

                        bufferIndex = (sixteenthDecimatedPicturePtr->origin_y + (sb_origin_y >> 2)) * sixteenthDecimatedPicturePtr->strideY + sixteenthDecimatedPicturePtr->origin_x + (sb_origin_x >> 2);

                        {
                            uint8_t  *framePtr = &sixteenthDecimatedPicturePtr->bufferY[bufferIndex];
                            uint8_t  *localPtr = context_ptr->me_context_ptr->sixteenth_sb_buffer;

                            for (lcuRow = 0; lcuRow < (sb_height >> 2); lcuRow += 2) {
                                EB_MEMCPY(localPtr, framePtr, (sb_width >> 2) * sizeof(uint8_t));
                                localPtr += 16;
                                framePtr += sixteenthDecimatedPicturePtr->strideY << 1;
                            }
                        }
                    }

                    MotionEstimateLcu(
                        picture_control_set_ptr,
                        sb_index,
                        sb_origin_x,
                        sb_origin_y,
                        context_ptr->me_context_ptr,
                        inputPicturePtr);

                }
            }
        }

        // *** OPEN LOOP INTRA CANDIDATE SEARCH CODE ***
        {

            // SB Loop
            for (yLcuIndex = yLcuStartIndex; yLcuIndex < yLcuEndIndex; ++yLcuIndex) {
                for (xLcuIndex = xLcuStartIndex; xLcuIndex < xLcuEndIndex; ++xLcuIndex) {

                    sb_origin_x = xLcuIndex * sequence_control_set_ptr->sb_sz;
                    sb_origin_y = yLcuIndex * sequence_control_set_ptr->sb_sz;

                    sb_index = (uint16_t)(xLcuIndex + yLcuIndex * picture_width_in_sb);


                    OpenLoopIntraSearchLcu(
                        picture_control_set_ptr,
                        sb_index,
                        context_ptr,
                        inputPicturePtr,
                        asm_type);



                }
            }
        }

        // ZZ SADs Computation
        // 1 lookahead frame is needed to get valid (0,0) SAD
        if (sequence_control_set_ptr->static_config.look_ahead_distance != 0) {
            // when DG is ON, the ZZ SADs are computed @ the PD process
            {
                // ZZ SADs Computation using decimated picture
                if (picture_control_set_ptr->picture_number > 0) {

                    ComputeDecimatedZzSad(
                        context_ptr,
                        sequence_control_set_ptr,
                        picture_control_set_ptr,
                        sixteenthDecimatedPicturePtr,
                        xLcuStartIndex,
                        xLcuEndIndex,
                        yLcuStartIndex,
                        yLcuEndIndex);

                }
            }
        }


        // Calculate the ME Distortion and OIS Historgrams

        EbBlockOnMutex(picture_control_set_ptr->rc_distortion_histogram_mutex);

        if (sequence_control_set_ptr->static_config.rate_control_mode) {
            if (picture_control_set_ptr->slice_type != I_SLICE) {
                uint16_t sadIntervalIndex;
                for (yLcuIndex = yLcuStartIndex; yLcuIndex < yLcuEndIndex; ++yLcuIndex) {
                    for (xLcuIndex = xLcuStartIndex; xLcuIndex < xLcuEndIndex; ++xLcuIndex) {

                        sb_origin_x = xLcuIndex * sequence_control_set_ptr->sb_sz;
                        sb_origin_y = yLcuIndex * sequence_control_set_ptr->sb_sz;
                        sb_width = (sequence_control_set_ptr->luma_width - sb_origin_x) < BLOCK_SIZE_64 ? sequence_control_set_ptr->luma_width - sb_origin_x : BLOCK_SIZE_64;
                        sb_height = (sequence_control_set_ptr->luma_height - sb_origin_y) < BLOCK_SIZE_64 ? sequence_control_set_ptr->luma_height - sb_origin_y : BLOCK_SIZE_64;

                        sb_index = (uint16_t)(xLcuIndex + yLcuIndex * picture_width_in_sb);
                        picture_control_set_ptr->inter_sad_interval_index[sb_index] = 0;
                        picture_control_set_ptr->intra_sad_interval_index[sb_index] = 0;

                        if (sb_width == BLOCK_SIZE_64 && sb_height == BLOCK_SIZE_64) {


                            sadIntervalIndex = (uint16_t)(picture_control_set_ptr->rc_me_distortion[sb_index] >> (12 - SAD_PRECISION_INTERVAL));//change 12 to 2*log2(64)

                            // printf("%d\n", sadIntervalIndex);

                            sadIntervalIndex = (uint16_t)(sadIntervalIndex >> 2);
                            if (sadIntervalIndex > (NUMBER_OF_SAD_INTERVALS >> 1) - 1) {
                                uint16_t sadIntervalIndexTemp = sadIntervalIndex - ((NUMBER_OF_SAD_INTERVALS >> 1) - 1);

                                sadIntervalIndex = ((NUMBER_OF_SAD_INTERVALS >> 1) - 1) + (sadIntervalIndexTemp >> 3);

                            }
                            if (sadIntervalIndex >= NUMBER_OF_SAD_INTERVALS - 1)
                                sadIntervalIndex = NUMBER_OF_SAD_INTERVALS - 1;


                            picture_control_set_ptr->inter_sad_interval_index[sb_index] = sadIntervalIndex;

                            picture_control_set_ptr->me_distortion_histogram[sadIntervalIndex] ++;

                            uint32_t                       bestOisCuIndex = 0;

                            //DOUBLE CHECK THIS PIECE OF CODE
                            intra_sad_interval_index = (uint32_t)
                                (((picture_control_set_ptr->ois_cu32_cu16_results[sb_index]->sorted_ois_candidate[1][bestOisCuIndex].distortion +
                                    picture_control_set_ptr->ois_cu32_cu16_results[sb_index]->sorted_ois_candidate[2][bestOisCuIndex].distortion +
                                    picture_control_set_ptr->ois_cu32_cu16_results[sb_index]->sorted_ois_candidate[3][bestOisCuIndex].distortion +
                                    picture_control_set_ptr->ois_cu32_cu16_results[sb_index]->sorted_ois_candidate[4][bestOisCuIndex].distortion)) >> (12 - SAD_PRECISION_INTERVAL));//change 12 to 2*log2(64) ;

                            intra_sad_interval_index = (uint16_t)(intra_sad_interval_index >> 2);
                            if (intra_sad_interval_index > (NUMBER_OF_SAD_INTERVALS >> 1) - 1) {
                                uint32_t sadIntervalIndexTemp = intra_sad_interval_index - ((NUMBER_OF_SAD_INTERVALS >> 1) - 1);

                                intra_sad_interval_index = ((NUMBER_OF_SAD_INTERVALS >> 1) - 1) + (sadIntervalIndexTemp >> 3);

                            }
                            if (intra_sad_interval_index >= NUMBER_OF_SAD_INTERVALS - 1)
                                intra_sad_interval_index = NUMBER_OF_SAD_INTERVALS - 1;


                            picture_control_set_ptr->intra_sad_interval_index[sb_index] = intra_sad_interval_index;

                            picture_control_set_ptr->ois_distortion_histogram[intra_sad_interval_index] ++;

                            ++picture_control_set_ptr->full_sb_count;
                        }

                    }
                }
            }
            else {
                uint32_t                       bestOisCuIndex = 0;


                for (yLcuIndex = yLcuStartIndex; yLcuIndex < yLcuEndIndex; ++yLcuIndex) {
                    for (xLcuIndex = xLcuStartIndex; xLcuIndex < xLcuEndIndex; ++xLcuIndex) {
                        sb_origin_x = xLcuIndex * sequence_control_set_ptr->sb_sz;
                        sb_origin_y = yLcuIndex * sequence_control_set_ptr->sb_sz;
                        sb_width = (sequence_control_set_ptr->luma_width - sb_origin_x) < BLOCK_SIZE_64 ? sequence_control_set_ptr->luma_width - sb_origin_x : BLOCK_SIZE_64;
                        sb_height = (sequence_control_set_ptr->luma_height - sb_origin_y) < BLOCK_SIZE_64 ? sequence_control_set_ptr->luma_height - sb_origin_y : BLOCK_SIZE_64;

                        sb_index = (uint16_t)(xLcuIndex + yLcuIndex * picture_width_in_sb);

                        picture_control_set_ptr->inter_sad_interval_index[sb_index] = 0;
                        picture_control_set_ptr->intra_sad_interval_index[sb_index] = 0;

                        if (sb_width == BLOCK_SIZE_64 && sb_height == BLOCK_SIZE_64) {


                            //DOUBLE CHECK THIS PIECE OF CODE

                            intra_sad_interval_index = (uint32_t)
                                (((picture_control_set_ptr->ois_cu32_cu16_results[sb_index]->sorted_ois_candidate[1][bestOisCuIndex].distortion +
                                    picture_control_set_ptr->ois_cu32_cu16_results[sb_index]->sorted_ois_candidate[2][bestOisCuIndex].distortion +
                                    picture_control_set_ptr->ois_cu32_cu16_results[sb_index]->sorted_ois_candidate[3][bestOisCuIndex].distortion +
                                    picture_control_set_ptr->ois_cu32_cu16_results[sb_index]->sorted_ois_candidate[4][bestOisCuIndex].distortion)) >> (12 - SAD_PRECISION_INTERVAL));//change 12 to 2*log2(64) ;

                            intra_sad_interval_index = (uint16_t)(intra_sad_interval_index >> 2);
                            if (intra_sad_interval_index > (NUMBER_OF_SAD_INTERVALS >> 1) - 1) {
                                uint32_t sadIntervalIndexTemp = intra_sad_interval_index - ((NUMBER_OF_SAD_INTERVALS >> 1) - 1);

                                intra_sad_interval_index = ((NUMBER_OF_SAD_INTERVALS >> 1) - 1) + (sadIntervalIndexTemp >> 3);

                            }
                            if (intra_sad_interval_index >= NUMBER_OF_SAD_INTERVALS - 1)
                                intra_sad_interval_index = NUMBER_OF_SAD_INTERVALS - 1;

                            picture_control_set_ptr->intra_sad_interval_index[sb_index] = intra_sad_interval_index;

                            picture_control_set_ptr->ois_distortion_histogram[intra_sad_interval_index] ++;

                            ++picture_control_set_ptr->full_sb_count;
                        }

                    }
                }
            }
        }

        EbReleaseMutex(picture_control_set_ptr->rc_distortion_histogram_mutex);

        // Get Empty Results Object
        EbGetEmptyObject(
            context_ptr->motionEstimationResultsOutputFifoPtr,
            &outputResultsWrapperPtr);

        outputResultsPtr = (MotionEstimationResults_t*)outputResultsWrapperPtr->objectPtr;
        outputResultsPtr->pictureControlSetWrapperPtr = inputResultsPtr->pictureControlSetWrapperPtr;
        outputResultsPtr->segment_index = segment_index;

        // Release the Input Results
        EbReleaseObject(inputResultsWrapperPtr);

        // Post the Full Results Object
        EbPostFullObject(outputResultsWrapperPtr);
    }
    return EB_NULL;
}
