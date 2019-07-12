/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EBMCP_H
#define EBMCP_H

#include "EbMcp_SSE2.h"
#include "EbMcp_SSSE3.h"

#include "EbDefinitions.h"
#include "EbUtility.h"
#include "EbPictureBufferDesc.h"
#include "EbPictureControlSet.h"
#include "EbSequenceControlSet.h"
#include "EbMotionEstimationContext.h"
#include "EbObject.h"

#ifdef __cplusplus
extern "C" {
#endif
#define USE_PRE_COMPUTE             0

    typedef struct MotionCompensationPredictionContext
    {
        EbDctor                  dctor;
        EbByte                   avc_style_mcp_intermediate_result_buf0;                    // For short filter in MD
        EbByte                   avc_style_mcp_intermediate_result_buf1;                    // For short filter in MD
#if !USE_PRE_COMPUTE
        EbByte                   avc_style_mcp_two_d_interpolation_first_pass_filter_result_buf; // For short filter in MD
#endif
    } MotionCompensationPredictionContext;

    /** InterpolationFilter()
            is generally defined interpolation filter function.
            There is a whole group of these functions, each of which corresponds to a particular
            integer/fractional sample, and the function is indexed in a function pointer array
            in terms of the frac_pos_x and frac_pos_y.

        @param *ref_pic (8-bits input)
            ref_pic is the pointer to the reference picture data that was chosen by
            the integer pixel precision MV.
        @param src_stride (input)
        @param frac_pos_x (input)
            frac_pos_x is the horizontal fractional position of the predicted sample
        @param frac_pos_y (input)
            frac_pos_y is the veritcal fractional position of the predicted sample
        @param pu_width (input)
        @param pu_height (input)
        @param *dst (16-bits output)
            dst is the pointer to the destination where the prediction result will
            be stored.
        @param dst_stride (input)
        @param *first_pass_if_dst (16-bits input)
            first_pass_if_dst is the pointer to the buffer where the result of the first
            pass filtering of the 2D interpolation filter will be stored.
        @param is_last (input)
            is_last indicates if there is any further filtering (interpolation filtering)
            afterwards.
     */
    typedef void(*InterpolationFilterNew)(
        EbByte               ref_pic,               //8-bits input parameter, please refer to the detailed explanation above.
        uint32_t             src_stride,            //input parameter
        EbByte               dst,                  //output parameter, please refer to the detailed explanation above.
        uint32_t             dst_stride,            //input parameter
        uint32_t             pu_width,              //input parameter
        uint32_t             pu_height,             //input parameter
        int16_t             *first_pass_if_dst);      //input parameter, please refer to the detailed explanation above.

    typedef void(*InterpolationFilterOutRaw)(
        EbByte              ref_pic,               //8-bits input parameter, please refer to the detailed explanation above.
        uint32_t            src_stride,            //input parameter
        int16_t            *dst,                  //output parameter, please refer to the detailed explanation above.
        uint32_t            pu_width,              //input parameter
        uint32_t            pu_height,             //input parameter
        int16_t            *first_pass_if_dst);      //input parameter, please refer to the detailed explanation above.

    typedef void(*ChromaFilterNew)(
        EbByte             ref_pic,
        uint32_t           src_stride,
        EbByte             dst,
        uint32_t           dst_stride,
        uint32_t           pu_width,
        uint32_t           pu_height,
        int16_t           *first_pass_if_dst,
        uint32_t           frac_pos_x,
        uint32_t           frac_pos_y);

    typedef void(*ChromaFilterOutRaw)(
        EbByte             ref_pic,
        uint32_t           src_stride,
        int16_t           *dst,
        uint32_t           pu_width,
        uint32_t           pu_height,
        int16_t           *first_pass_if_dst,
        uint32_t           frac_pos_x,
        uint32_t           frac_pos_y);
    extern EbErrorType in_loop_me_context_ctor(
        SsMeContext                         *ss_mecontext);

    extern void generate_padding(
        EbByte              src_pic,
        uint32_t            src_stride,
        uint32_t            original_src_width,
        uint32_t            original_src_height,
        uint32_t            padding_width,
        uint32_t            padding_height);

    extern void generate_padding16_bit(
        EbByte              src_pic,
        uint32_t            src_stride,
        uint32_t            original_src_width,
        uint32_t            original_src_height,
        uint32_t            padding_width,
        uint32_t            padding_height);

    extern void pad_input_picture(
        EbByte              src_pic,
        uint32_t            src_stride,
        uint32_t            original_src_width,
        uint32_t            original_src_height,
        uint32_t            pad_right,
        uint32_t            pad_bottom);

    // Function Tables (Super-long, declared in EbMcpTables.c)
    extern const InterpolationFilterNew     uni_pred_luma_if_function_ptr_array_new[ASM_TYPE_TOTAL][16];
    extern const InterpolationFilterOutRaw  bi_pred_luma_if_function_ptr_array_new[ASM_TYPE_TOTAL][16];
    extern const ChromaFilterNew            uni_pred_chroma_if_function_ptr_array_new[ASM_TYPE_TOTAL][64];
    extern const ChromaFilterOutRaw         bi_pred_chroma_if_function_ptr_array_new[ASM_TYPE_TOTAL][64];

#ifdef __cplusplus
}
#endif
#endif // EBMCP_H
