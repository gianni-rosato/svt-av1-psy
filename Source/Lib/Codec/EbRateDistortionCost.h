/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbRateDistortionCost_h
#define EbRateDistortionCost_h


/***************************************
 * Includes
 ***************************************/
#include "EbIntraPrediction.h"
#include "EbInterPrediction.h"
#include "EbLambdaRateTables.h"
#include "EbTransforms.h"
#include "EbModeDecisionProcess.h"
#include "EbEncDecProcess.h"
#include "EbEntropyCoding.h"

#ifdef __cplusplus
extern "C" {
#endif
    extern void CodingLoopContextGeneration(
        ModeDecisionContext_t      *context_ptr,
        CodingUnit_t            *cu_ptr,
        uint32_t                   cu_origin_x,
        uint32_t                   cu_origin_y,
        uint32_t                   sb_sz,

        NeighborArrayUnit_t     *skip_coeff_neighbor_array,
        NeighborArrayUnit_t     *luma_dc_sign_level_coeff_neighbor_array,
        NeighborArrayUnit_t     *cb_dc_sign_level_coeff_neighbor_array,
        NeighborArrayUnit_t     *cr_dc_sign_level_coeff_neighbor_array,
        NeighborArrayUnit_t     *inter_pred_dir_neighbor_array,
        NeighborArrayUnit_t     *ref_frame_type_neighbor_array,

        NeighborArrayUnit_t     *intraLumaNeighborArray,
        NeighborArrayUnit_t     *skip_flag_neighbor_array,
        NeighborArrayUnit_t     *mode_type_neighbor_array,
        NeighborArrayUnit_t     *leaf_depth_neighbor_array,
        NeighborArrayUnit_t       *leaf_partition_neighbor_array);

    extern EbErrorType Av1TuCalcCost(
        ModeDecisionCandidate_t *candidate_ptr,                        // input parameter, prediction result Ptr
        int16_t                   txbSkipCtx,
        uint32_t                   tu_index,                             // input parameter, TU index inside the CU
        uint32_t                   y_count_non_zero_coeffs,                 // input parameter, number of non zero Y quantized coefficients
        uint32_t                   cbCountNonZeroCoeffs,                // input parameter, number of non zero cb quantized coefficients
        uint32_t                   crCountNonZeroCoeffs,                // input parameter, number of non zero cr quantized coefficients
        uint64_t                   yTuDistortion[DIST_CALC_TOTAL],      // input parameter, Y distortion for both Normal and Cbf zero modes
        uint64_t                   cbTuDistortion[DIST_CALC_TOTAL],     // input parameter, Cb distortion for both Normal and Cbf zero modes
        uint64_t                   crTuDistortion[DIST_CALC_TOTAL],     // input parameter, Cr distortion for both Normal and Cbf zero modes
        COMPONENT_TYPE           componentType,
        uint64_t                  *yTuCoeffBits,                        // input parameter, Y quantized coefficients rate
        uint64_t                  *cbTuCoeffBits,                       // input parameter, Cb quantized coefficients rate
        uint64_t                  *crTuCoeffBits,                       // input parameter, Cr quantized coefficients rate
        TxSize                  txsize,
        uint64_t                   lambda);                              // input parameter, lambda for Luma


    extern EbErrorType TuCalcCost(
        uint32_t                   cu_size,
        ModeDecisionCandidate_t *candidate_ptr,
        uint32_t                   tu_index,
        uint32_t                   transform_size,
        uint32_t                   transform_chroma_size,
        uint32_t                   y_count_non_zero_coeffs,
        uint32_t                   cbCountNonZeroCoeffs,
        uint32_t                   crCountNonZeroCoeffs,
        uint64_t                   yTuDistortion[DIST_CALC_TOTAL],
        uint64_t                   cbTuDistortion[DIST_CALC_TOTAL],
        uint64_t                   crTuDistortion[DIST_CALC_TOTAL],
        uint32_t                   component_mask,
        uint64_t                  *yTuCoeffBits,
        uint64_t                  *cbTuCoeffBits,
        uint64_t                  *crTuCoeffBits,
        uint32_t                   qp,
        uint64_t                   lambda,
        uint64_t                   lambda_chroma);

    extern EbErrorType Av1TuCalcCostLuma(
        int16_t                   txbSkipCtx,
        ModeDecisionCandidate_t *candidate_ptr,                        // input parameter, prediction result Ptr
        uint32_t                   tu_index,                             // input parameter, TU index inside the CU
        TxSize                  txSize,
        uint32_t                   y_count_non_zero_coeffs,                 // input parameter, number of non zero Y quantized coefficients
        uint64_t                   yTuDistortion[DIST_CALC_TOTAL],      // input parameter, Y distortion for both Normal and Cbf zero modes
        uint64_t                  *yTuCoeffBits,                        // input parameter, Y quantized coefficients rate
        uint64_t                  *yFullCost,
        uint64_t                   lambda);                              // input parameter, lambda for Luma


    extern EbErrorType IntraLumaModeContext(
        CodingUnit_t                        *cu_ptr,
        uint32_t                                  lumaMode,
        int32_t                                 *predictionIndex);

    extern EbErrorType Intra2Nx2NFastCostIsliceOpt(
        struct ModeDecisionContext_s           *context_ptr,
        CodingUnit_t                          *cu_ptr,
        struct ModeDecisionCandidateBuffer_s   *candidate_buffer_ptr,
        uint32_t                                  qp,
        uint64_t                                  luma_distortion,
        uint64_t                                  chroma_distortion,
        uint64_t                                  lambda,
        PictureControlSet_t                    *picture_control_set_ptr);

    extern EbErrorType Intra2Nx2NFastCostIslice(
        CodingUnit_t                          *cu_ptr,
        struct ModeDecisionCandidateBuffer_s   *candidate_buffer_ptr,
        uint32_t                                  qp,
        uint64_t                                  luma_distortion,
        uint64_t                                  chroma_distortion,
        uint64_t                                  lambda,
        PictureControlSet_t                    *picture_control_set_ptr);

    extern EbErrorType Intra2Nx2NFastCostPsliceOpt(
        struct ModeDecisionContext_s           *context_ptr,
        CodingUnit_t                           *cu_ptr,
        struct ModeDecisionCandidateBuffer_s        *candidate_buffer_ptr,
        uint32_t                                  qp,
        uint64_t                                  luma_distortion,
        uint64_t                                  chroma_distortion,
        uint64_t                                  lambda,
        PictureControlSet_t                    *picture_control_set_ptr);





    EbErrorType IntraFullLumaCostIslice(
        CodingUnit_t                           *cu_ptr,
        uint32_t                                  cu_size,
        uint32_t                                  cu_size_log2,
        struct ModeDecisionCandidateBuffer_s   *candidate_buffer_ptr,
        uint64_t                                 *y_distortion,
        uint64_t                                  lambda,
        uint64_t                                 *y_coeff_bits,
        uint32_t                                  transform_size);

    EbErrorType IntraFullLumaCostPslice(
        CodingUnit_t                           *cu_ptr,
        uint32_t                                  cu_size,
        uint32_t                                  cu_size_log2,
        struct ModeDecisionCandidateBuffer_s   *candidate_buffer_ptr,
        uint64_t                                 *y_distortion,
        uint64_t                                  lambda,
        uint64_t                                 *y_coeff_bits,
        uint32_t                                  transform_size);

    extern EbErrorType InterFastCostPsliceOpt(
        struct ModeDecisionContext_s           *context_ptr,
        CodingUnit_t                        *cu_ptr,
        struct ModeDecisionCandidateBuffer_s   *candidate_buffer_ptr,
        uint32_t                                  qp,
        uint64_t                                  luma_distortion,
        uint64_t                                  chroma_distortion,
        uint64_t                                  lambda,
        PictureControlSet_t                    *picture_control_set_ptr);

    extern EbErrorType InterFastCostBsliceOpt(
        struct ModeDecisionContext_s           *context_ptr,
        CodingUnit_t                           *cu_ptr,
        struct ModeDecisionCandidateBuffer_s       *candidate_buffer_ptr,
        uint32_t                                  qp,
        uint64_t                                  luma_distortion,
        uint64_t                                  chroma_distortion,
        uint64_t                                  lambda,
        PictureControlSet_t                    *picture_control_set_ptr);


    EbErrorType InterFullLumaCost(
        CodingUnit_t                           *cu_ptr,
        uint32_t                                  cu_size,
        uint32_t                                  cu_size_log2,
        ModeDecisionCandidateBuffer_t          *candidate_buffer_ptr,
        uint64_t                                 *y_distortion,
        uint64_t                                  lambda,
        uint64_t                                 *y_coeff_bits,
        uint32_t                                  transform_size);


    extern EbErrorType  MergeSkipFullCost(

        LargestCodingUnit_t                    *sb_ptr,
        CodingUnit_t                           *cu_ptr,
        uint32_t                                  cu_size,
        uint32_t                                  cu_size_log2,
        ModeDecisionCandidateBuffer_t          *candidate_buffer_ptr,
        uint32_t                                  qp,
        uint64_t                                 *y_distortion,
        uint64_t                                 *cb_distortion,
        uint64_t                                 *cr_distortion,
        uint64_t                                  lambda,
        uint64_t                                  lambda_chroma,
        uint64_t                                 *y_coeff_bits,
        uint64_t                                 *cb_coeff_bits,
        uint64_t                                 *cr_coeff_bits,
        uint32_t                                  transform_size,
        uint32_t                                  transform_chroma_size,
        PictureControlSet_t                    *picture_control_set_ptr);

    extern EbErrorType SplitFlagRate(
        ModeDecisionContext_t               *context_ptr,
        CodingUnit_t                           *cu_ptr,
        uint32_t                                  split_flag,
        uint64_t                                 *splitRate,
        uint64_t                                  lambda,
        MdRateEstimationContext_t              *md_rate_estimation_ptr,
        uint32_t                                  tbMaxDepth);

#define RDDIV_BITS 7

#define RDCOST(RM, R, D)                                            \
  (ROUND_POWER_OF_TWO(((uint64_t)(R)) * (RM), AV1_PROB_COST_SHIFT) + \
   ((D) * (1 << RDDIV_BITS)))
    extern EbErrorType Av1SplitFlagRate(
        SequenceControlSet_t                  *sequence_control_set_ptr,
        ModeDecisionContext_t                  *context_ptr,
        CodingUnit_t                           *cu_ptr,
        uint32_t                                  leaf_index,
        PartitionType                          partitionType,
        uint64_t                                 *splitRate,
        uint64_t                                  lambda,
        MdRateEstimationContext_t              *md_rate_estimation_ptr,
        uint32_t                                  tbMaxDepth);

    extern EbErrorType Av1EncodeTuCalcCost(
        EncDecContext_t                        *context_ptr,
        uint32_t                                 *count_non_zero_coeffs,
        uint64_t                                  yTuDistortion[DIST_CALC_TOTAL],
        uint64_t                                 *yTuCoeffBits,
        uint32_t                                  component_mask
    );

    extern EbErrorType EncodeTuCalcCost(
        EncDecContext_t          *context_ptr,
        uint32_t                   *count_non_zero_coeffs,
        uint64_t                    yTuDistortion[DIST_CALC_TOTAL],
        uint64_t                   *yTuCoeffBits,
        uint32_t                    component_mask
    );



    extern EbErrorType Intra4x4FastCostIslice(
        ModeDecisionContext_t                  *context_ptr,
        uint32_t                                  pu_index,
        struct ModeDecisionCandidateBuffer_s   *candidate_buffer_ptr,
        uint64_t                                  luma_distortion,
        uint64_t                                  lambda);

    extern EbErrorType Intra4x4FastCostPslice(
        ModeDecisionContext_t                  *context_ptr,
        uint32_t                                  pu_index,
        struct ModeDecisionCandidateBuffer_s   *candidate_buffer_ptr,
        uint64_t                                  luma_distortion,
        uint64_t                                  lambda);

    extern EbErrorType Intra4x4FullCostIslice(
        struct ModeDecisionCandidateBuffer_s   *candidate_buffer_ptr,
        uint64_t                                 *y_distortion,
        uint64_t                                  lambda,
        uint64_t                                 *y_coeff_bits,
        uint32_t                                  transform_size);

    extern EbErrorType Intra4x4FullCostPslice(
        struct ModeDecisionCandidateBuffer_s   *candidate_buffer_ptr,
        uint64_t                                 *y_distortion,
        uint64_t                                  lambda,
        uint64_t                                 *y_coeff_bits,
        uint32_t                                  transform_size);


    extern EbErrorType Av1IntraFastCost(
        struct ModeDecisionContext_s            *context_ptr,
        CodingUnit_t                            *cu_ptr,
        ModeDecisionCandidateBuffer_t            *candidate_buffer_ptr,
        uint32_t                                  qp,
        uint64_t                                  luma_distortion,
        uint64_t                                  chroma_distortion,
        uint64_t                                  lambda,
        PictureControlSet_t                        *picture_control_set_ptr);

    extern EbErrorType Av1InterFastCost(
        struct ModeDecisionContext_s            *context_ptr,
        CodingUnit_t                            *cu_ptr,
        ModeDecisionCandidateBuffer_t            *candidate_buffer_ptr,
        uint32_t                                  qp,
        uint64_t                                  luma_distortion,
        uint64_t                                  chroma_distortion,
        uint64_t                                  lambda,
        PictureControlSet_t                        *picture_control_set_ptr);

    extern EbErrorType Av1IntraFullCost(
        PictureControlSet_t                    *picture_control_set_ptr,
        ModeDecisionContext_t                  *context_ptr,
        struct ModeDecisionCandidateBuffer_s   *candidate_buffer_ptr,
        CodingUnit_t                           *cu_ptr,
        uint64_t                                 *y_distortion,
        uint64_t                                 *cb_distortion,
        uint64_t                                 *cr_distortion,
        uint64_t                                  lambda,
        uint64_t                                 *y_coeff_bits,
        uint64_t                                 *cb_coeff_bits,
        uint64_t                                 *cr_coeff_bits,
        BlockSize                              bsize);



    extern EbErrorType Av1InterFullCost(
        PictureControlSet_t                    *picture_control_set_ptr,
        ModeDecisionContext_t                  *context_ptr,
        struct ModeDecisionCandidateBuffer_s   *candidate_buffer_ptr,
        CodingUnit_t                           *cu_ptr,
        uint64_t                                 *y_distortion,
        uint64_t                                 *cb_distortion,
        uint64_t                                 *cr_distortion,
        uint64_t                                  lambda,
        uint64_t                                 *y_coeff_bits,
        uint64_t                                 *cb_coeff_bits,
        uint64_t                                 *cr_coeff_bits,
        BlockSize                              bsize);





#ifdef __cplusplus
}
#endif
#endif //EbRateDistortionCost_h
