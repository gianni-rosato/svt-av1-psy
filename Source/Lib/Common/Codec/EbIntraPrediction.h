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

    struct ModeDecisionContext;

    typedef void(*IntraPredFnC)(uint8_t *dst, ptrdiff_t stride, int32_t w, int32_t h,
        const uint8_t *above, const uint8_t *left);
    typedef void(*IntraHighBdPredFnC)(uint16_t *dst, ptrdiff_t stride, int32_t w, int32_t h,
        const uint16_t *above, const uint16_t *left, int32_t bd);

    typedef void(*IntraPredFn)(uint8_t *dst, ptrdiff_t stride,
        const uint8_t *above, const uint8_t *left);

    typedef void(*IntraHighPredFn)(uint16_t *dst, ptrdiff_t stride,
        const uint16_t *above, const uint16_t *left,
        int32_t bd);

    typedef struct IntraReferenceSamples 
    {

        uint8_t                  *y_intra_reference_array;
        uint8_t                  *cb_intra_reference_array;
        uint8_t                  *cr_intra_reference_array;
        uint8_t                  *y_intra_filtered_reference_array;

        uint8_t                  *y_intra_reference_array_reverse;
        uint8_t                  *y_intra_filtered_reference_array_reverse;
        uint8_t                  *cb_intra_reference_array_reverse;
        uint8_t                  *cr_intra_reference_array_reverse;

        // Scratch buffers used in the interpolaiton process
        uint8_t                   reference_above_line_y[(MAX_PU_SIZE << 2) + 1];
        uint8_t                   reference_left_line_y[(MAX_PU_SIZE << 2) + 1];
        EbBool                 above_ready_flag_y;
        EbBool                 left_ready_flag_y;

        uint8_t                   reference_above_line_cb[(MAX_PU_SIZE << 2) + 2];
        uint8_t                   reference_left_line_cb[(MAX_PU_SIZE << 2) + 2];
        EbBool                 above_ready_flag_cb;
        EbBool                 left_ready_flag_cb;

        uint8_t                   reference_above_line_cr[(MAX_PU_SIZE << 2) + 2];
        uint8_t                   reference_left_line_cr[(MAX_PU_SIZE << 2) + 2];
        EbBool                 above_ready_flag_cr;
        EbBool                 left_ready_flag_cr;

    } IntraReferenceSamples;

    typedef struct IntraReference16bitSamples 
    {

        uint16_t                  *y_intra_reference_array;
        uint16_t                  *cb_intra_reference_array;
        uint16_t                  *cr_intra_reference_array;
        uint16_t                  *y_intra_filtered_reference_array;

        uint16_t                  *y_intra_reference_array_reverse;
        uint16_t                  *y_intra_filtered_reference_array_reverse;
        uint16_t                  *cb_intra_reference_array_reverse;
        uint16_t                  *cr_intra_reference_array_reverse;

        // Scratch buffers used in the interpolaiton process
        uint16_t                   reference_above_line_y[(MAX_PU_SIZE << 2) + 1];
        uint16_t                   reference_left_line_y[(MAX_PU_SIZE << 2) + 1];
        EbBool                  above_ready_flag_y;
        EbBool                  left_ready_flag_y;

        uint16_t                   reference_above_line_cb[(MAX_PU_SIZE << 2) + 2];
        uint16_t                   reference_left_line_cb[(MAX_PU_SIZE << 2) + 2];
        EbBool                  above_ready_flag_cb;
        EbBool                  left_ready_flag_cb;

        uint16_t                   reference_above_line_cr[(MAX_PU_SIZE << 2) + 2];
        uint16_t                   reference_left_line_cr[(MAX_PU_SIZE << 2) + 2];
        EbBool                  above_ready_flag_cr;
        EbBool                  left_ready_flag_cr;

    } IntraReference16bitSamples;

    extern EbErrorType intra_reference_samples_ctor(
        IntraReferenceSamples **context_dbl_ptr);

    extern EbErrorType intra_reference16bit_samples_ctor(
        IntraReference16bitSamples **context_dbl_ptr);

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


    extern EbErrorType av1_intra_prediction_cl(
        struct ModeDecisionContext           *context_ptr,
        PictureControlSet                    *picture_control_set_ptr,
        ModeDecisionCandidateBuffer           *candidate_buffer_ptr,
        EbAsm                                  asm_type);

    extern void intra_mode_angular_horizontal_kernel_ssse3_intrin(
        uint32_t            size,
        uint8_t            *ref_samp_main,
        uint8_t            *prediction_ptr,
        uint32_t            prediction_buffer_stride,
        const EbBool     skip,
        int32_t            intra_pred_angle);



    extern EbErrorType intra_open_loop_reference_samples_ctor(
        IntraReferenceSamplesOpenLoop **context_dbl_ptr);

    extern EbErrorType update_neighbor_samples_array_open_loop(
        uint8_t                           *above_ref,
        uint8_t                            *left_ref,
        EbPictureBufferDesc              *input_ptr,
        uint32_t                            stride,
        uint32_t                            srcOriginX,
        uint32_t                            srcOriginY,
        uint8_t                             bwidth,
        uint8_t                             bheight);
    extern EbErrorType intra_prediction_open_loop(
         int32_t  p_angle ,
        uint8_t                          ois_intra_mode,
        uint32_t                         srcOriginX,
        uint32_t                         srcOriginY,
        TxSize                          tx_size,
        uint8_t                         *above_row,
        uint8_t                         *left_col,
        MotionEstimationContext_t       *context_ptr);                  // input parameter, ME context

    extern void intra_mode_planar(
        const uint32_t   size,                       //input parameter, denotes the size of the current PU
        uint8_t         *ref_samples,                 //input parameter, pointer to the reference samples
        uint8_t         *dst,              //output parameter, pointer to the prediction
        const uint32_t   prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
        const EbBool  skip);
    extern void ebav1_smooth_v_predictor(
        const uint32_t   size,                       //input parameter, denotes the size of the current PU
        uint8_t         *ref_samples,                 //input parameter, pointer to the reference samples
        uint8_t         *dst,              //output parameter, pointer to the prediction
        const uint32_t   prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
        const EbBool  skip);
    extern void ebav1_smooth_h_predictor(
        const uint32_t   size,                       //input parameter, denotes the size of the current PU
        uint8_t         *ref_samples,                 //input parameter, pointer to the reference samples
        uint8_t         *dst,              //output parameter, pointer to the prediction
        const uint32_t   prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
        const EbBool  skip);

    void intra_mode_angular_av1_z1_16bit(
        const uint32_t   size,                    //input parameter, denotes the size of the current PU
        uint16_t         *ref_samples,             //input parameter, pointer to the reference samples
        uint16_t         *dst,                    //output parameter, pointer to the prediction
        const uint32_t   prediction_buffer_stride,  //input parameter, denotes the stride for the prediction ptr
        const EbBool  skip,
        uint16_t          dx,                     //output parameter, pointer to the prediction
        uint16_t          dy,                      //output parameter, pointer to the prediction
        uint16_t          bd);

    void intra_mode_angular_av1_z2_16bit(
        const uint32_t   size,                       //input parameter, denotes the size of the current PU
        uint16_t         *ref_samples,                 //input parameter, pointer to the reference samples
        uint16_t         *dst,              //output parameter, pointer to the prediction
        const uint32_t   prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
        const EbBool  skip,
        uint16_t          dx,              //output parameter, pointer to the prediction
        uint16_t          dy,              //output parameter, pointer to the prediction
        uint16_t          bd);

    void intra_mode_angular_av1_z3_16bit(
        const uint32_t   size,                       //input parameter, denotes the size of the current PU
        uint16_t         *ref_samples,                 //input parameter, pointer to the reference samples
        uint16_t         *dst,              //output parameter, pointer to the prediction
        const uint32_t   prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
        const EbBool  skip,
        uint16_t          dx,              //output parameter, pointer to the prediction
        uint16_t          dy,              //output parameter, pointer to the prediction
        uint16_t          bd);


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
        CflPredType pred_type) {
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
