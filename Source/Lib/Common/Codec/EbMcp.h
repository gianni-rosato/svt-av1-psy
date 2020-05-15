/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EBMCP_H
#define EBMCP_H


#ifdef __cplusplus
extern "C" {
#endif
#define USE_PRE_COMPUTE 0

typedef struct MotionCompensationPredictionContext {
    EbDctor dctor;
    EbByte  avc_style_mcp_intermediate_result_buf0; // For short filter in MD
    EbByte  avc_style_mcp_intermediate_result_buf1; // For short filter in MD
#if !USE_PRE_COMPUTE
    EbByte avc_style_mcp_two_d_interpolation_first_pass_filter_result_buf; // For short filter in MD
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
extern void generate_padding(EbByte src_pic, uint32_t src_stride, uint32_t original_src_width,
                             uint32_t original_src_height, uint32_t padding_width,
                             uint32_t padding_height);

extern void generate_padding16_bit(EbByte src_pic, uint32_t src_stride, uint32_t original_src_width,
                                   uint32_t original_src_height, uint32_t padding_width,
                                   uint32_t padding_height);

extern void pad_input_picture(EbByte src_pic, uint32_t src_stride, uint32_t original_src_width,
                              uint32_t original_src_height, uint32_t pad_right,
                              uint32_t pad_bottom);

extern void pad_input_picture_16bit(uint16_t* src_pic, uint32_t src_stride,
                                    uint32_t original_src_width, uint32_t original_src_height,
                                    uint32_t pad_right, uint32_t pad_bottom);

    void generate_padding_l(EbByte src_pic, uint32_t src_stride,
        uint32_t row_height, uint32_t padding_width);
    void generate_padding_r(EbByte src_pic, uint32_t src_stride,
        uint32_t row_width, uint32_t row_height, uint32_t padding_width);
    void generate_padding_t(EbByte src_pic, uint32_t src_stride,
        uint32_t row_width, uint32_t padding_height);
    void generate_padding_b(EbByte src_pic, uint32_t src_stride,
        uint32_t row_width, uint32_t row_height, uint32_t padding_height);

    void generate_padding_l_hbd(EbByte src_pic, uint32_t src_stride,
        uint32_t row_height, uint32_t padding_width);
    void generate_padding_r_hbd(EbByte src_pic, uint32_t src_stride,
        uint32_t row_width, uint32_t row_height, uint32_t padding_width);

#ifdef __cplusplus
}
#endif
#endif // EBMCP_H
