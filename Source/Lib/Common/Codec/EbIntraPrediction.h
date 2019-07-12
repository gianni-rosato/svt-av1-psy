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
#include "EbObject.h"

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
        EbDctor                   dctor;
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
        EbDctor                    dctor;
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

/////####.... For recursive intra prediction.....#####///

#define FILTER_INTRA_SCALE_BITS 4
extern const int8_t av1_filter_intra_taps[FILTER_INTRA_MODES][8][8];

/////####.... To make functions common between EbIntraPrediction.c &
void *aom_memset16(void *dest, int32_t val, size_t length);

int32_t use_intra_edge_upsample(int32_t bs0, int32_t bs1, int32_t delta,
                                       int32_t type);

BlockSize scale_chroma_bsize(BlockSize bsize, int32_t subsampling_x,
                              int32_t subsampling_y);

int32_t intra_edge_filter_strength(int32_t bs0, int32_t bs1, int32_t delta, int32_t type);

 enum {
    NEED_LEFT = 1 << 1,
    NEED_ABOVE = 1 << 2,
    NEED_ABOVERIGHT = 1 << 3,
    NEED_ABOVELEFT = 1 << 4,
    NEED_BOTTOMLEFT = 1 << 5,
};

extern const uint8_t extend_modes[INTRA_MODES];

/* TODO: Need to harmonize with fun from EbAdaptiveMotionVectorPrediction.c */
int32_t intra_has_top_right(BlockSize   sb_size, BlockSize bsize, int32_t mi_row,
    int32_t mi_col, int32_t top_available, int32_t right_available,
    PartitionType partition, TxSize txsz, int32_t row_off,
    int32_t col_off, int32_t ss_x, int32_t ss_y);

extern int32_t intra_has_bottom_left(BlockSize sb_size, BlockSize bsize, int32_t mi_row,
    int32_t mi_col, int32_t bottom_available, int32_t left_available,
    PartitionType partition, TxSize txsz, int32_t row_off,
    int32_t col_off, int32_t ss_x, int32_t ss_y);

extern IntraPredFn pred[INTRA_MODES][TX_SIZES_ALL];
extern IntraPredFn dc_pred[2][2][TX_SIZES_ALL];

extern IntraHighPredFn pred_high[INTRA_MODES][TX_SIZES_ALL];
extern IntraHighPredFn dc_pred_high[2][2][TX_SIZES_ALL];

void dr_predictor(uint8_t *dst, ptrdiff_t stride, TxSize tx_size,
    const uint8_t *above, const uint8_t *left,
    int32_t upsample_above, int32_t upsample_left, int32_t angle);

void filter_intra_edge_corner(uint8_t *p_above, uint8_t *p_left);

void highbd_dr_predictor(uint16_t *dst, ptrdiff_t stride,
    TxSize tx_size, const uint16_t *above,
    const uint16_t *left, int32_t upsample_above,
    int32_t upsample_left, int32_t angle, int32_t bd);

void filter_intra_edge_corner_high(uint16_t *p_above, uint16_t *p_left);

void highbd_filter_intra_predictor(uint16_t *dst, ptrdiff_t stride,
                                          TxSize tx_size,
                                          const uint16_t *above,
                                          const uint16_t *left, int mode,
                                          int bd);

/////////..............................................//////////////////////////

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
        IntraReferenceSamplesOpenLoop *context_ptr);

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
        EbPictureBufferDesc           *input_ptr,
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

    static EB_INTRA_NOANG_TYPE FUNC_TABLE IntraSmoothH_Av1_funcPtrArray[ASM_TYPE_TOTAL] = {
        // NON_AVX2
        ebav1_smooth_h_predictor,
        // AVX2
        ebav1_smooth_h_predictor,
    };
    static EB_INTRA_NOANG_TYPE FUNC_TABLE IntraSmoothV_Av1_funcPtrArray[ASM_TYPE_TOTAL] = {
        // NON_AVX2
        ebav1_smooth_v_predictor,
        // AVX2
        ebav1_smooth_v_predictor,
    };

    static EB_INTRA_ANG_Z1_Z2_Z3_16bit_TYPE FUNC_TABLE IntraModeAngular_AV1_Z1_16bit_funcPtrArray[9][ASM_TYPE_TOTAL] = {
        // 4x4
        {
            // NON_AVX2
            intra_mode_angular_av1_z1_16bit,
            // AVX2
            intra_mode_angular_av1_z1_16bit_4x4_avx2,
        },
        // 8x8
        {
            // NON_AVX2
            intra_mode_angular_av1_z1_16bit,
            // AVX2
            intra_mode_angular_av1_z1_16bit_8x8_avx2,
        },
        // 16x16
        {
            // NON_AVX2
            intra_mode_angular_av1_z1_16bit,
            // AVX2
            intra_mode_angular_av1_z1_16bit_16x16_avx2,
        },
        // NxN
        {
            // NON_AVX2
            intra_mode_angular_av1_z1_16bit,
            // AVX2
            intra_mode_angular_av1_z1_16bit,
        },
        // 32x32
        {
            // NON_AVX2
            intra_mode_angular_av1_z1_16bit,
            // AVX2
            intra_mode_angular_av1_z1_16bit_32x32_avx2,
        },
        // NxN
        {
            // NON_AVX2
            intra_mode_angular_av1_z1_16bit,
            // AVX2
            intra_mode_angular_av1_z1_16bit,
        },
        // NxN
        {
            // NON_AVX2
            intra_mode_angular_av1_z1_16bit,
            // AVX2
            intra_mode_angular_av1_z1_16bit,
        },
        // NxN
        {
            // NON_AVX2
            intra_mode_angular_av1_z1_16bit,
            // AVX2
            intra_mode_angular_av1_z1_16bit,
        },
        // 64x64
        {
            // NON_AVX2
            intra_mode_angular_av1_z1_16bit,
            // AVX2
            intra_mode_angular_av1_z1_16bit_64x64_avx2,
        }
    };
    static EB_INTRA_ANG_Z1_Z2_Z3_16bit_TYPE FUNC_TABLE IntraModeAngular_AV1_Z2_16bit_funcPtrArray[9][ASM_TYPE_TOTAL] = {
        // 4x4
        {
            // NON_AVX2
            intra_mode_angular_av1_z2_16bit,
            // AVX2
            intra_mode_angular_av1_z2_16bit_4x4_avx2,
        },
        // 8x8
        {
            // NON_AVX2
            intra_mode_angular_av1_z2_16bit,
            // AVX2
            intra_mode_angular_av1_z2_16bit_8x8_avx2,
        },
        // 16x16
        {
            // NON_AVX2
            intra_mode_angular_av1_z2_16bit,
            // AVX2
            intra_mode_angular_av1_z2_16bit_16x16_avx2,
        },
        // NxN
        {
            // NON_AVX2
            intra_mode_angular_av1_z2_16bit,
            // AVX2
            intra_mode_angular_av1_z2_16bit,
        },
        // 32x32
        {
            // NON_AVX2
            intra_mode_angular_av1_z2_16bit,
            // AVX2
            intra_mode_angular_av1_z2_16bit_32x32_avx2,
        },
        // NxN
        {
            // NON_AVX2
            intra_mode_angular_av1_z2_16bit,
            // AVX2
            intra_mode_angular_av1_z2_16bit,
        },
        // NxN
        {
            // NON_AVX2
            intra_mode_angular_av1_z2_16bit,
            // AVX2
            intra_mode_angular_av1_z2_16bit,
        },
        // NxN
        {
            // NON_AVX2
            intra_mode_angular_av1_z2_16bit,
            // AVX2
            intra_mode_angular_av1_z2_16bit,
        },
        // 64x64
        {
            // NON_AVX2
            intra_mode_angular_av1_z2_16bit,
            // AVX2
            intra_mode_angular_av1_z2_16bit_64x64_avx2,
        }
    };
    static EB_INTRA_ANG_Z1_Z2_Z3_16bit_TYPE FUNC_TABLE IntraModeAngular_AV1_Z3_16bit_funcPtrArray[9][ASM_TYPE_TOTAL] = {
        // 4x4
        {
            // NON_AVX2
            intra_mode_angular_av1_z3_16bit,
            // AVX2
            intra_mode_angular_av1_z3_16bit_4x4_avx2,
        },
        // 8x8
        {
            // NON_AVX2
            intra_mode_angular_av1_z3_16bit,
            // AVX2
            intra_mode_angular_av1_z3_16bit_8x8_avx2,
        },
        // 16x16
        {
            // NON_AVX2
            intra_mode_angular_av1_z3_16bit,
            // AVX2
            intra_mode_angular_av1_z3_16bit_16x16_avx2,
        },
        // NxN
        {
            // NON_AVX2
            intra_mode_angular_av1_z3_16bit,
            // AVX2
            intra_mode_angular_av1_z3_16bit,
        },
        // 32x32
        {
            // NON_AVX2
            intra_mode_angular_av1_z3_16bit,
            // AVX2
            intra_mode_angular_av1_z3_16bit_32x32_avx2,
        },
        // NxN
        {
            // NON_AVX2
            intra_mode_angular_av1_z3_16bit,
            // AVX2
            intra_mode_angular_av1_z3_16bit,
        },
        // NxN
        {
            // NON_AVX2
            intra_mode_angular_av1_z3_16bit,
            // AVX2
            intra_mode_angular_av1_z3_16bit,
        },
        // NxN
        {
            // NON_AVX2
            intra_mode_angular_av1_z3_16bit,
            // AVX2
            intra_mode_angular_av1_z3_16bit,
        },
        // 64x64
        {
            // NON_AVX2
            intra_mode_angular_av1_z3_16bit,
            // AVX2
            intra_mode_angular_av1_z3_16bit_64x64_avx2,
        }
    };

typedef struct CflCtx {
        // Q3 reconstructed luma pixels (only Q2 is required, but Q3 is used to avoid
        // shifts)
        int16_t recon_buf_q3[CFL_BUF_SQUARE];

        // Height and width currently used in the CfL prediction buffer.
        int32_t buf_height, buf_width;

        int32_t are_parameters_computed;

        // Chroma subsampling
        int32_t subsampling_x, subsampling_y;
} CflCtx;

    extern void cfl_luma_subsampling_420_lbd_c(
        const uint8_t *input, // AMIR-> Changed to 8 bit
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

/* Function pointers return by CfL functions */

typedef void(*cfl_subsample_lbd_fn)(const uint8_t *input, int input_stride,
    int16_t *output_q3);

typedef void(*cfl_subsample_hbd_fn)(const uint16_t *input, int input_stride,
    int16_t *output_q3);

typedef void(*cfl_subtract_average_fn)(int16_t *dst);

typedef void(*cfl_predict_lbd_fn)(const int16_t *pred_buf_q3, uint8_t *pred, int32_t pred_stride,
    uint8_t *dst, int32_t dst_stride, int32_t alpha_q3, int32_t bit_depth, int32_t width, int32_t height);

typedef void(*cfl_predict_hbd_fn)(const int16_t *src, uint16_t *dst,
    int dst_stride, int alpha_q3, int bd);

#define cfl_get_luma_subsampling_420_hbd cfl_get_luma_subsampling_420_hbd_c
cfl_subsample_hbd_fn cfl_get_luma_subsampling_420_hbd_c(TxSize tx_size);

#define cfl_get_luma_subsampling_420_lbd cfl_get_luma_subsampling_420_lbd_c
cfl_subsample_lbd_fn cfl_get_luma_subsampling_420_lbd_c(TxSize tx_size);

#define cfl_get_luma_subsampling_422_hbd cfl_get_luma_subsampling_422_hbd_c
cfl_subsample_hbd_fn cfl_get_luma_subsampling_422_hbd_c(TxSize tx_size);

#define cfl_get_luma_subsampling_422_lbd cfl_get_luma_subsampling_422_lbd_c
cfl_subsample_lbd_fn cfl_get_luma_subsampling_422_lbd_c(TxSize tx_size);

#define cfl_get_luma_subsampling_444_hbd cfl_get_luma_subsampling_444_hbd_c
cfl_subsample_hbd_fn cfl_get_luma_subsampling_444_hbd_c(TxSize tx_size);

#define cfl_get_luma_subsampling_444_lbd cfl_get_luma_subsampling_444_lbd_c
cfl_subsample_lbd_fn cfl_get_luma_subsampling_444_lbd_c(TxSize tx_size);

cfl_subtract_average_fn get_subtract_average_fn_c(TxSize tx_size);
#define get_subtract_average_fn get_subtract_average_fn_c

   // Allows the CFL_SUBSAMPLE function to switch types depending on the bitdepth.
#define CFL_lbd_TYPE uint8_t *cfl_type
#define CFL_hbd_TYPE uint16_t *cfl_type

    // Declare a size-specific wrapper for the size-generic function. The compiler
    // will inline the size generic function in here, the advantage is that the size
    // will be constant allowing for loop unrolling and other constant propagated
    // goodness.
#define CFL_SUBSAMPLE(arch, sub, bd, width, height)                       \
      void subsample_##bd##_##sub##_##width##x##height##_##arch(              \
          const CFL_##bd##_TYPE, int input_stride, int16_t *output_q3) {     \
        cfl_luma_subsampling_##sub##_##bd##_##arch(cfl_type, input_stride,    \
                                                   output_q3, width, height); \
      }

    // Declare size-specific wrappers for all valid CfL sizes.
#define CFL_SUBSAMPLE_FUNCTIONS(arch, sub, bd)                            \
      CFL_SUBSAMPLE(arch, sub, bd, 4, 4)                                      \
      CFL_SUBSAMPLE(arch, sub, bd, 8, 8)                                      \
      CFL_SUBSAMPLE(arch, sub, bd, 16, 16)                                    \
      CFL_SUBSAMPLE(arch, sub, bd, 32, 32)                                    \
      CFL_SUBSAMPLE(arch, sub, bd, 4, 8)                                      \
      CFL_SUBSAMPLE(arch, sub, bd, 8, 4)                                      \
      CFL_SUBSAMPLE(arch, sub, bd, 8, 16)                                     \
      CFL_SUBSAMPLE(arch, sub, bd, 16, 8)                                     \
      CFL_SUBSAMPLE(arch, sub, bd, 16, 32)                                    \
      CFL_SUBSAMPLE(arch, sub, bd, 32, 16)                                    \
      CFL_SUBSAMPLE(arch, sub, bd, 4, 16)                                     \
      CFL_SUBSAMPLE(arch, sub, bd, 16, 4)                                     \
      CFL_SUBSAMPLE(arch, sub, bd, 8, 32)                                     \
      CFL_SUBSAMPLE(arch, sub, bd, 32, 8)                                     \
      cfl_subsample_##bd##_fn cfl_get_luma_subsampling_##sub##_##bd##_##arch( \
          TxSize tx_size) {                                                  \
        CFL_SUBSAMPLE_FUNCTION_ARRAY(arch, sub, bd)                           \
        return subfn_##sub[tx_size];                                          \
      }

    // Declare an architecture-specific array of function pointers for size-specific
    // wrappers.
#define CFL_SUBSAMPLE_FUNCTION_ARRAY(arch, sub, bd)                       \
             const cfl_subsample_##bd##_fn subfn_##sub[TX_SIZES_ALL] = {      \
        subsample_##bd##_##sub##_4x4_##arch,   /* 4x4 */                      \
        subsample_##bd##_##sub##_8x8_##arch,   /* 8x8 */                      \
        subsample_##bd##_##sub##_16x16_##arch, /* 16x16 */                    \
        subsample_##bd##_##sub##_32x32_##arch, /* 32x32 */                    \
        NULL,                                  /* 64x64 (invalid CFL size) */ \
        subsample_##bd##_##sub##_4x8_##arch,   /* 4x8 */                      \
        subsample_##bd##_##sub##_8x4_##arch,   /* 8x4 */                      \
        subsample_##bd##_##sub##_8x16_##arch,  /* 8x16 */                     \
        subsample_##bd##_##sub##_16x8_##arch,  /* 16x8 */                     \
        subsample_##bd##_##sub##_16x32_##arch, /* 16x32 */                    \
        subsample_##bd##_##sub##_32x16_##arch, /* 32x16 */                    \
        NULL,                                  /* 32x64 (invalid CFL size) */ \
        NULL,                                  /* 64x32 (invalid CFL size) */ \
        subsample_##bd##_##sub##_4x16_##arch,  /* 4x16  */                    \
        subsample_##bd##_##sub##_16x4_##arch,  /* 16x4  */                    \
        subsample_##bd##_##sub##_8x32_##arch,  /* 8x32  */                    \
        subsample_##bd##_##sub##_32x8_##arch,  /* 32x8  */                    \
        NULL,                                  /* 16x64 (invalid CFL size) */ \
        NULL,                                  /* 64x16 (invalid CFL size) */ \
      };

    // The RTCD script does not support passing in an array, so we wrap it in this
    // function.
#define CFL_GET_SUBSAMPLE_FUNCTION(arch)  \
      CFL_SUBSAMPLE_FUNCTIONS(arch, 420, lbd) \
      CFL_SUBSAMPLE_FUNCTIONS(arch, 422, lbd) \
      CFL_SUBSAMPLE_FUNCTIONS(arch, 444, lbd) \
      CFL_SUBSAMPLE_FUNCTIONS(arch, 420, hbd) \
      CFL_SUBSAMPLE_FUNCTIONS(arch, 422, hbd) \
      CFL_SUBSAMPLE_FUNCTIONS(arch, 444, hbd)

    // Declare a size-specific wrapper for the size-generic function. The compiler
    // will inline the size generic function in here, the advantage is that the size
    // will be constant allowing for loop unrolling and other constant propagated
    // goodness.
#define CFL_SUB_AVG_X(arch, width, height, round_offset, num_pel_log2)   \
  void subtract_average_##width##x##height##_##arch(int16_t *buf) {      \
    subtract_average_##arch(buf, width, height, round_offset,       \
                            num_pel_log2);                               \
  }

    // Declare size-specific wrappers for all valid CfL sizes.
#define CFL_SUB_AVG_FN(arch)                                                \
      CFL_SUB_AVG_X(arch, 4, 4, 8, 4)                                           \
      CFL_SUB_AVG_X(arch, 4, 8, 16, 5)                                          \
      CFL_SUB_AVG_X(arch, 4, 16, 32, 6)                                         \
      CFL_SUB_AVG_X(arch, 8, 4, 16, 5)                                          \
      CFL_SUB_AVG_X(arch, 8, 8, 32, 6)                                          \
      CFL_SUB_AVG_X(arch, 8, 16, 64, 7)                                         \
      CFL_SUB_AVG_X(arch, 8, 32, 128, 8)                                        \
      CFL_SUB_AVG_X(arch, 16, 4, 32, 6)                                         \
      CFL_SUB_AVG_X(arch, 16, 8, 64, 7)                                         \
      CFL_SUB_AVG_X(arch, 16, 16, 128, 8)                                       \
      CFL_SUB_AVG_X(arch, 16, 32, 256, 9)                                       \
      CFL_SUB_AVG_X(arch, 32, 8, 128, 8)                                        \
      CFL_SUB_AVG_X(arch, 32, 16, 256, 9)                                       \
      CFL_SUB_AVG_X(arch, 32, 32, 512, 10)                                      \
      cfl_subtract_average_fn get_subtract_average_fn_##arch(TxSize tx_size) { \
              const cfl_subtract_average_fn sub_avg[TX_SIZES_ALL] = {          \
          subtract_average_4x4_##arch,   /* 4x4 */                              \
          subtract_average_8x8_##arch,   /* 8x8 */                              \
          subtract_average_16x16_##arch, /* 16x16 */                            \
          subtract_average_32x32_##arch, /* 32x32 */                            \
          NULL,                          /* 64x64 (invalid CFL size) */         \
          subtract_average_4x8_##arch,   /* 4x8 */                              \
          subtract_average_8x4_##arch,   /* 8x4 */                              \
          subtract_average_8x16_##arch,  /* 8x16 */                             \
          subtract_average_16x8_##arch,  /* 16x8 */                             \
          subtract_average_16x32_##arch, /* 16x32 */                            \
          subtract_average_32x16_##arch, /* 32x16 */                            \
          NULL,                          /* 32x64 (invalid CFL size) */         \
          NULL,                          /* 64x32 (invalid CFL size) */         \
          subtract_average_4x16_##arch,  /* 4x16 (invalid CFL size) */          \
          subtract_average_16x4_##arch,  /* 16x4 (invalid CFL size) */          \
          subtract_average_8x32_##arch,  /* 8x32 (invalid CFL size) */          \
          subtract_average_32x8_##arch,  /* 32x8 (invalid CFL size) */          \
          NULL,                          /* 16x64 (invalid CFL size) */         \
          NULL,                          /* 64x16 (invalid CFL size) */         \
        };                                                                      \
        /* Modulo TX_SIZES_ALL to ensure that an attacker won't be able to */   \
        /* index the function pointer array out of bounds. */                   \
        return sub_avg[tx_size % TX_SIZES_ALL];                                 \
      }

#ifdef __cplusplus
}
#endif
#endif // EbIntraPrediction_h
