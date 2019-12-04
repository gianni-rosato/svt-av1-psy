/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbAvcStyleMcp.h"
#include "EbPictureOperators.h"
#include "aom_dsp_rtcd.h"

static const   uint8_t  frac_mapped_pos_tab_x[16] = { 0, 1, 2, 3,
    0, 1, 2, 3,
    0, 0, 2, 0,
    0, 1, 2, 3 };

static const   uint8_t  frac_mapped_pos_tab_y[16] = { 0, 0, 0, 0,
    1, 0, 0, 0,
    2, 2, 2, 2,
    3, 0, 0, 0 };

static const   uint8_t  integer_posoffset_tab_x[16] = { 0, 0, 0, 0,
    0, 0, 0, 0,
    0, 0, 0, 1,
    0, 0, 0, 0 };

static const   uint8_t  integer_posoffset_tab_y[16] = { 0, 0, 0, 0,
    0, 0, 0, 0,
    0, 0, 0, 0,
    0, 1, 1, 1 };

static const int8_t avc_style_luma_if_coeff[4][4] = {
    {0, 0, 0, 0},
    {-1, 25, 9, -1},
    {-2, 18, 18, -2},
    {-1, 9, 25, -1},
};

void avc_style_copy_c(
    EbByte                 ref_pic,
    uint32_t               src_stride,
    EbByte                 dst,
    uint32_t               dst_stride,
    uint32_t               pu_width,
    uint32_t               pu_height,
    EbByte                 temp_buf,
    uint32_t               frac_pos)
{
    (void)temp_buf;
    (void)frac_pos;
    picture_copy_kernel(ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, 1);
}

void estimate_uni_pred_interpolation_unpacked_avc_style(
    EbPictureBufferDesc   *ref_pic,
    uint32_t                 pos_x,
    uint32_t                 pos_y,
    uint32_t                 pu_width,
    uint32_t                 pu_height,
    EbPictureBufferDesc   *dst,
    uint32_t                 dst_luma_index,
    uint32_t                 dst_chroma_index,          //input parameter, please refer to the detailed explanation above.
    uint32_t                 component_mask,
    EbByte                  temp_buf,
    EbBool                   sub_sample_pred_flag)
{
    uint8_t    frac_posx, mapped_frac_posx;
    uint8_t    frac_posy, mapped_frac_posy;
    uint32_t   luma_stride;

    luma_stride = dst->stride_y;

    frac_posx = pos_x & 0x03;
    frac_posy = pos_y & 0x03;

    mapped_frac_posx = frac_posx;
    mapped_frac_posy = frac_posy;

    (void)component_mask;
    (void)dst_chroma_index;
    //doing the luma interpolation
    avc_style_luma_interpolation_filter(
        ref_pic->buffer_y + 2 + 2 * ref_pic->stride_y, ref_pic->stride_y,
        dst->buffer_y + dst_luma_index, luma_stride,
        pu_width, pu_height,
        temp_buf,
        sub_sample_pred_flag,
        mapped_frac_posx ? mapped_frac_posx : mapped_frac_posy,
        mapped_frac_posx + (mapped_frac_posy << 2));
}

/*******************************************************************************
* Requirement: pu_width      = 8, 16, 24, 32, 48 or 64
* Requirement: pu_height % 2 = 0
* Requirement: skip         = 0 or 1
* Requirement (x86 only): temp_buf % 16 = 0
* Requirement (x86 only): (dst->buffer_y  + dst_luma_index  ) % 16 = 0 when pu_width %16 = 0
* Requirement (x86 only): (dst->buffer_cb + dst_chroma_index) % 16 = 0 when pu_width %32 = 0
* Requirement (x86 only): (dst->buffer_cr + dst_chroma_index) % 16 = 0 when pu_width %32 = 0
* Requirement (x86 only): dst->stride_y   % 16 = 0 when pu_width %16 = 0
* Requirement (x86 only): dst->chromaStride % 16 = 0 when pu_width %32 = 0
*******************************************************************************/
void estimate_bi_pred_interpolation_unpacked_avc_style(
    EbPictureBufferDesc   *ref_pic_list0,
    EbPictureBufferDesc   *ref_pic_list1,
    uint32_t                 ref_list0_pos_x,
    uint32_t                 ref_list0_pos_y,
    uint32_t                 ref_list1_pos_x,
    uint32_t                 ref_list1_pos_y,
    uint32_t                 pu_width,
    uint32_t                 pu_height,
    EbPictureBufferDesc   *bi_dst,
    uint32_t                 dst_luma_index,
    uint32_t                 dst_chroma_index,
    uint32_t                 component_mask,
    EbByte                  ref_list0_temp_dst,
    EbByte                  ref_list1_temp_dst,
    EbByte                  first_pass_if_temp_dst,
    EbBool                   sub_sample_pred_flag)
{
    uint32_t   integ_pos_x;
    uint32_t   integ_pos_y;
    uint8_t    frac_posx, mapped_frac_posx;
    uint8_t    frac_posy, mapped_frac_posy;

    uint32_t   luma_stride, ref_luma_stride;

    luma_stride = bi_dst->stride_y;
    ref_luma_stride = ref_pic_list0->stride_y;

    (void)component_mask;
    (void)dst_chroma_index;
    uint8_t pad = 2;
    integ_pos_x = pad;
    integ_pos_y = pad;

    frac_posx = ref_list0_pos_x & 0x03;
    frac_posy = ref_list0_pos_y & 0x03;

    mapped_frac_posx = frac_posx;
    mapped_frac_posy = frac_posy;

    avc_style_luma_interpolation_filter(
        ref_pic_list0->buffer_y + integ_pos_x + integ_pos_y * ref_luma_stride, ref_luma_stride,
        ref_list0_temp_dst, pu_width,
        pu_width, pu_height,
        first_pass_if_temp_dst,
        sub_sample_pred_flag,
        mapped_frac_posx ? mapped_frac_posx : mapped_frac_posy,
        mapped_frac_posx + (mapped_frac_posy << 2));

    frac_posx = ref_list1_pos_x & 0x03;
    frac_posy = ref_list1_pos_y & 0x03;

    mapped_frac_posx = frac_posx;
    mapped_frac_posy = frac_posy;

    //doing the luma interpolation
    avc_style_luma_interpolation_filter(
        ref_pic_list1->buffer_y + integ_pos_x + integ_pos_y * ref_luma_stride, ref_luma_stride,
        ref_list1_temp_dst, pu_width,
        pu_width, pu_height,
        first_pass_if_temp_dst,
        sub_sample_pred_flag,
        mapped_frac_posx ? mapped_frac_posx : mapped_frac_posy,
        mapped_frac_posx + (mapped_frac_posy << 2));

    // bi-pred luma
    picture_average_kernel(ref_list0_temp_dst, pu_width << sub_sample_pred_flag, ref_list1_temp_dst, pu_width << sub_sample_pred_flag, bi_dst->buffer_y + dst_luma_index, luma_stride << sub_sample_pred_flag, pu_width, pu_height >> sub_sample_pred_flag);
    if (sub_sample_pred_flag)
        picture_average_kernel1_line(ref_list0_temp_dst + (pu_height - 1)*pu_width, ref_list1_temp_dst + (pu_height - 1)*pu_width, bi_dst->buffer_y + dst_luma_index + (pu_height - 1)*luma_stride, pu_width);
}

/*******************************************************************************
* Requirement: pu_width      = 8, 16, 24, 32, 48 or 64
* Requirement: pu_height % 2 = 0
* Requirement: skip         = 0 or 1
* Requirement (x86 only): temp_buf % 16 = 0
* Requirement (x86 only): (dst->buffer_y  + dst_luma_index  ) % 16 = 0 when pu_width %16 = 0
* Requirement (x86 only): (dst->buffer_cb + dst_chroma_index) % 16 = 0 when pu_width %32 = 0
* Requirement (x86 only): (dst->buffer_cr + dst_chroma_index) % 16 = 0 when pu_width %32 = 0
* Requirement (x86 only): dst->stride_y   % 16 = 0 when pu_width %16 = 0
* Requirement (x86 only): dst->chromaStride % 16 = 0 when pu_width %32 = 0
*******************************************************************************/
void estimate_uni_pred_interpolation_avc_luma(
    EbPictureBufferDesc   *ref_pic,
    uint32_t                 pos_x,
    uint32_t                 pos_y,
    uint32_t                 pu_width,
    uint32_t                 pu_height,
    EbPictureBufferDesc   *dst,
    uint32_t                 dst_luma_index,
    uint32_t                 dst_chroma_index,          //input parameter, please refer to the detailed explanation above.
    uint32_t                 component_mask,
    EbByte                  temp_buf,
    EbBool                   sub_sample_pred_flag)
{
    uint32_t   integ_pos_x;
    uint32_t   integ_pos_y;
    uint8_t    frac_posx, mapped_frac_posx;
    uint8_t    frac_posy, mapped_frac_posy;
    uint32_t   luma_stride;
    uint32_t   chroma_pu_width = pu_width >> 1;
    uint32_t   chroma_pu_height = pu_height >> 1;

    uint8_t    frac_pos;

    luma_stride = dst->stride_y;

    //luma
    //compute the luma fractional position
    integ_pos_x = (pos_x >> 2);
    integ_pos_y = (pos_y >> 2);

    frac_posx = pos_x & 0x03;
    frac_posy = pos_y & 0x03;

    //New Interpolation Mapping
    mapped_frac_posx = frac_posx;
    mapped_frac_posy = frac_posy;

    frac_pos = (mapped_frac_posx + (mapped_frac_posy << 2));
    mapped_frac_posx = frac_mapped_pos_tab_x[frac_pos];
    mapped_frac_posy = frac_mapped_pos_tab_y[frac_pos];
    integ_pos_x += integer_posoffset_tab_x[frac_pos];
    integ_pos_y += integer_posoffset_tab_y[frac_pos];
    if (component_mask & PICTURE_BUFFER_DESC_LUMA_MASK)
    {
        //doing the luma interpolation
        avc_style_luma_interpolation_filter(
            ref_pic->buffer_y + integ_pos_x + integ_pos_y * ref_pic->stride_y, ref_pic->stride_y,
            dst->buffer_y + dst_luma_index, luma_stride,
            pu_width, pu_height,
            temp_buf,
            sub_sample_pred_flag,
            mapped_frac_posx ? mapped_frac_posx : mapped_frac_posy,
            mapped_frac_posx + (mapped_frac_posy << 2));
    }
    //chroma
    if (component_mask & PICTURE_BUFFER_DESC_CHROMA_MASK) {
        //compute the chroma fractional position
        integ_pos_x = (pos_x >> 3);
        integ_pos_y = (pos_y >> 3);
        frac_posx = pos_x & 0x07;
        frac_posy = pos_y & 0x07;

        mapped_frac_posx = 0;
        if (frac_posx > 4)
            integ_pos_x++;
        mapped_frac_posy = 0;
        if (frac_posy > 4)
            integ_pos_y++;

        // Note: chroma_pu_width equals 4 is only supported in Intrinsic
       //       for integer positions ( mapped_frac_posx + (mapped_frac_posy << 3) equals 0 )
       //doing the chroma Cb interpolation
        avc_style_luma_interpolation_filter(
            ref_pic->buffer_cb + integ_pos_x + integ_pos_y * ref_pic->stride_cb,
            ref_pic->stride_cb,
            dst->buffer_cb + dst_chroma_index,
            dst->stride_cb,
            chroma_pu_width,
            chroma_pu_height,
            temp_buf,
            sub_sample_pred_flag,
            mapped_frac_posx ? mapped_frac_posx : mapped_frac_posy,
            mapped_frac_posx + (mapped_frac_posy << 3));

        //doing the chroma Cr interpolation
        avc_style_luma_interpolation_filter(
            ref_pic->buffer_cr + integ_pos_x + integ_pos_y * ref_pic->stride_cr,
            ref_pic->stride_cr,
            dst->buffer_cr + dst_chroma_index,
            dst->stride_cr,
            chroma_pu_width,
            chroma_pu_height,
            temp_buf,
            sub_sample_pred_flag,
            mapped_frac_posx ? mapped_frac_posx : mapped_frac_posy,
            mapped_frac_posx + (mapped_frac_posy << 3));
    }
}

/*******************************************************************************
 * Requirement: pu_width      = 8, 16, 24, 32, 48 or 64
 * Requirement: pu_height % 2 = 0
 * Requirement: skip         = 0 or 1
 * Requirement (x86 only): temp_buf % 16 = 0
 * Requirement (x86 only): (dst->buffer_y  + dst_luma_index  ) % 16 = 0 when pu_width %16 = 0
 * Requirement (x86 only): (dst->buffer_cb + dst_chroma_index) % 16 = 0 when pu_width %32 = 0
 * Requirement (x86 only): (dst->buffer_cr + dst_chroma_index) % 16 = 0 when pu_width %32 = 0
 * Requirement (x86 only): dst->stride_y   % 16 = 0 when pu_width %16 = 0
 * Requirement (x86 only): dst->chromaStride % 16 = 0 when pu_width %32 = 0
*******************************************************************************/
void estimate_bi_pred_interpolation_avc_luma(
    EbPictureBufferDesc   *ref_pic_list0,
    EbPictureBufferDesc   *ref_pic_list1,
    uint32_t                 ref_list0_pos_x,
    uint32_t                 ref_list0_pos_y,
    uint32_t                 ref_list1_pos_x,
    uint32_t                 ref_list1_pos_y,
    uint32_t                 pu_width,
    uint32_t                 pu_height,
    EbPictureBufferDesc   *bi_dst,
    uint32_t                 dst_luma_index,
    uint32_t                 dst_chroma_index,
    uint32_t                 component_mask,
    EbByte                  ref_list0_temp_dst,
    EbByte                  ref_list1_temp_dst,
    EbByte                  first_pass_if_temp_dst,
    EbBool                   sub_sample_pred_flag)
{
    uint32_t   integ_pos_x;
    uint32_t   integ_pos_y;
    uint8_t    frac_posx, mapped_frac_posx;
    uint8_t    frac_posy, mapped_frac_posy;
    uint8_t    frac_pos;

    uint32_t   luma_stride, ref_luma_stride;

    uint32_t   chroma_pu_width = pu_width >> 1;
    uint32_t   chroma_pu_height = pu_height >> 1;
    luma_stride = bi_dst->stride_y;
    ref_luma_stride = ref_pic_list0->stride_y;
    uint8_t shift = sub_sample_pred_flag ? 1 : 0;

    //Luma
    if (component_mask & PICTURE_BUFFER_DESC_LUMA_MASK) {
        //uni-prediction List0 luma
        //compute the luma fractional position
        integ_pos_x = (ref_list0_pos_x >> 2);
        integ_pos_y = (ref_list0_pos_y >> 2);

        frac_posx = ref_list0_pos_x & 0x03;
        frac_posy = ref_list0_pos_y & 0x03;

        mapped_frac_posx = frac_posx;
        mapped_frac_posy = frac_posy;

        frac_pos = (mapped_frac_posx + (mapped_frac_posy << 2));
        mapped_frac_posx = frac_mapped_pos_tab_x[frac_pos];
        mapped_frac_posy = frac_mapped_pos_tab_y[frac_pos];
        integ_pos_x += integer_posoffset_tab_x[frac_pos];
        integ_pos_y += integer_posoffset_tab_y[frac_pos];

        avc_style_luma_interpolation_filter(
            ref_pic_list0->buffer_y + integ_pos_x + integ_pos_y * ref_luma_stride, ref_luma_stride,
            ref_list0_temp_dst, pu_width,
            pu_width, pu_height,
            first_pass_if_temp_dst,
            sub_sample_pred_flag,
            mapped_frac_posx ? mapped_frac_posx : mapped_frac_posy,
            mapped_frac_posx + (mapped_frac_posy << 2));

        //uni-prediction List1 luma
        //compute the luma fractional position
        integ_pos_x = (ref_list1_pos_x >> 2);
        integ_pos_y = (ref_list1_pos_y >> 2);

        frac_posx = ref_list1_pos_x & 0x03;
        frac_posy = ref_list1_pos_y & 0x03;

        mapped_frac_posx = frac_posx;
        mapped_frac_posy = frac_posy;

        frac_pos = (mapped_frac_posx + (mapped_frac_posy << 2));
        mapped_frac_posx = frac_mapped_pos_tab_x[frac_pos];
        mapped_frac_posy = frac_mapped_pos_tab_y[frac_pos];
        integ_pos_x += integer_posoffset_tab_x[frac_pos];
        integ_pos_y += integer_posoffset_tab_y[frac_pos];
        //doing the luma interpolation
        avc_style_luma_interpolation_filter(
            ref_pic_list1->buffer_y + integ_pos_x + integ_pos_y * ref_luma_stride, ref_luma_stride,
            ref_list1_temp_dst, pu_width,
            pu_width, pu_height,
            first_pass_if_temp_dst,
            sub_sample_pred_flag,
            mapped_frac_posx ? mapped_frac_posx : mapped_frac_posy,
            mapped_frac_posx + (mapped_frac_posy << 2));

        // bi-pred luma
        picture_average_kernel(ref_list0_temp_dst, pu_width << sub_sample_pred_flag, ref_list1_temp_dst, pu_width << sub_sample_pred_flag, bi_dst->buffer_y + dst_luma_index, luma_stride << sub_sample_pred_flag, pu_width, pu_height >> sub_sample_pred_flag);
        if (sub_sample_pred_flag)
            picture_average_kernel1_line(ref_list0_temp_dst + (pu_height - 1)*pu_width, ref_list1_temp_dst + (pu_height - 1)*pu_width, bi_dst->buffer_y + dst_luma_index + (pu_height - 1)*luma_stride, pu_width);
    }

    //uni-prediction List0 chroma
    if (component_mask & PICTURE_BUFFER_DESC_CHROMA_MASK) {
        // bi-pred chroma  Cb
        // Note: chroma_pu_width equals 4 is only supported in Intrinsic
        //       for integer positions ( mapped_frac_posx + (mapped_frac_posy << 3) equals 0 )
        //doing the chroma Cb interpolation list 0
        //compute the chroma fractional position
        integ_pos_x = (ref_list0_pos_x >> 3);
        integ_pos_y = (ref_list0_pos_y >> 3);
        frac_posx = ref_list0_pos_x & 0x07;
        frac_posy = ref_list0_pos_y & 0x07;

        mapped_frac_posx = 0;
        if (frac_posx > 4)
            integ_pos_x++;
        mapped_frac_posy = 0;
        if (frac_posy > 4)
            integ_pos_y++;

        avc_style_luma_interpolation_filter(
            ref_pic_list0->buffer_cb + integ_pos_x + integ_pos_y * ref_pic_list0->stride_cb,
            ref_pic_list0->stride_cb,
            ref_list0_temp_dst,
            chroma_pu_width,
            chroma_pu_width,
            chroma_pu_height,
            first_pass_if_temp_dst,
            sub_sample_pred_flag,
            mapped_frac_posx ? mapped_frac_posx : mapped_frac_posy,
            mapped_frac_posx + (mapped_frac_posy << 3));

        //doing the chroma Cb interpolation list 1

        integ_pos_x = (ref_list1_pos_x >> 3);
        integ_pos_y = (ref_list1_pos_y >> 3);
        frac_posx = ref_list1_pos_x & 0x07;
        frac_posy = ref_list1_pos_y & 0x07;

        mapped_frac_posx = 0;
        if (frac_posx > 4)
            integ_pos_x++;
        mapped_frac_posy = 0;
        if (frac_posy > 4)
            integ_pos_y++;

        avc_style_luma_interpolation_filter(
            ref_pic_list1->buffer_cb + integ_pos_x + integ_pos_y * ref_pic_list1->stride_cb,
            ref_pic_list1->stride_cb,
            ref_list1_temp_dst,
            chroma_pu_width,
            chroma_pu_width,
            chroma_pu_height,
            first_pass_if_temp_dst,
            sub_sample_pred_flag,
            mapped_frac_posx ? mapped_frac_posx : mapped_frac_posy,
            mapped_frac_posx + (mapped_frac_posy << 3));

        // bi-pred Chroma Cb
        picture_average_kernel(
            ref_list0_temp_dst,
            chroma_pu_width << shift,
            ref_list1_temp_dst,
            chroma_pu_width << shift,
            bi_dst->buffer_cb + dst_chroma_index,
            bi_dst->stride_cb << shift,
            chroma_pu_width,
            chroma_pu_height >> shift
            );

        // bi-pred chroma  Cr
        // Note: chroma_pu_width equals 4 is only supported in Intrinsic
        //       for integer positions ( mapped_frac_posx + (mapped_frac_posy << 3) equals 0 )
        //doing the chroma Cb interpolation list 0
        //compute the chroma fractional position
        integ_pos_x = (ref_list0_pos_x >> 3);
        integ_pos_y = (ref_list0_pos_y >> 3);
        frac_posx = ref_list0_pos_x & 0x07;
        frac_posy = ref_list0_pos_y & 0x07;

        mapped_frac_posx = 0;
        if (frac_posx > 4)
            integ_pos_x++;
        mapped_frac_posy = 0;
        if (frac_posy > 4)
            integ_pos_y++;

        avc_style_luma_interpolation_filter(
            ref_pic_list0->buffer_cr + integ_pos_x + integ_pos_y * ref_pic_list0->stride_cr,
            ref_pic_list0->stride_cr,
            ref_list0_temp_dst,
            chroma_pu_width,
            chroma_pu_width,
            chroma_pu_height,
            first_pass_if_temp_dst,
            sub_sample_pred_flag,
            mapped_frac_posx ? mapped_frac_posx : mapped_frac_posy,
            mapped_frac_posx + (mapped_frac_posy << 3));

        //doing the chroma Cb interpolation list 1

        integ_pos_x = (ref_list1_pos_x >> 3);
        integ_pos_y = (ref_list1_pos_y >> 3);
        frac_posx = ref_list1_pos_x & 0x07;
        frac_posy = ref_list1_pos_y & 0x07;

        mapped_frac_posx = 0;
        if (frac_posx > 4)
            integ_pos_x++;
        mapped_frac_posy = 0;
        if (frac_posy > 4)
            integ_pos_y++;

        avc_style_luma_interpolation_filter(
            ref_pic_list1->buffer_cr + integ_pos_x + integ_pos_y * ref_pic_list1->stride_cr,
            ref_pic_list1->stride_cr,
            ref_list1_temp_dst,
            chroma_pu_width,
            chroma_pu_width,
            chroma_pu_height,
            first_pass_if_temp_dst,
            sub_sample_pred_flag,
            mapped_frac_posx ? mapped_frac_posx : mapped_frac_posy,
            mapped_frac_posx + (mapped_frac_posy << 3));

        // bi-pred Chroma Cr
        picture_average_kernel(
            ref_list0_temp_dst,
            chroma_pu_width << shift,
            ref_list1_temp_dst,
            chroma_pu_width << shift,
            bi_dst->buffer_cr + dst_chroma_index,
            bi_dst->stride_cr << shift,
            chroma_pu_width,
            chroma_pu_height >> shift
            );
    }
}

void estimate_uni_pred_interpolation_avc_lumaRef10Bit(
    EbPictureBufferDesc   *ref_frame_pic_list0,
    uint32_t                 pos_x,
    uint32_t                 pos_y,
    uint32_t                 pu_width,
    uint32_t                 pu_height,
    EbPictureBufferDesc   *dst,
    uint32_t                 dst_luma_index,
    uint32_t                 dst_chroma_index,          //input parameter, please refer to the detailed explanation above.
    uint32_t                 component_mask,
    EbByte                  temp_buf,
    EbBool                   sub_pred,
    EbBool                   sub_pred_chroma)
{
    uint32_t   chroma_pu_width = pu_width >> 1;
    uint32_t   chroma_pu_height = pu_height >> 1;

    //doing the luma interpolation
    if (component_mask & PICTURE_BUFFER_DESC_LUMA_MASK)
    {
        (void)temp_buf;
        uint32_t in_pos_x = (pos_x >> 2);
        uint32_t in_pos_y = (pos_y >> 2);
        uint16_t *ptr16 = (uint16_t *)ref_frame_pic_list0->buffer_y + in_pos_x + in_pos_y * ref_frame_pic_list0->stride_y;

        extract8_bitdata_safe_sub(
            ptr16,
            ref_frame_pic_list0->stride_y << sub_pred,
            dst->buffer_y + dst_luma_index,
            dst->stride_y << sub_pred,
            pu_width,
            pu_height >> sub_pred,
            sub_pred
        );
    }
    //chroma
    if (component_mask & PICTURE_BUFFER_DESC_CHROMA_MASK) {
        {
            sub_pred = sub_pred_chroma;
            uint32_t in_pos_x = (pos_x >> 3);
            uint32_t in_pos_y = (pos_y >> 3);
            uint16_t *ptr16 = (uint16_t *)ref_frame_pic_list0->buffer_cb + in_pos_x + in_pos_y * ref_frame_pic_list0->stride_cb;

            extract_8bit_data(
                ptr16,
                ref_frame_pic_list0->stride_cb << sub_pred,
                dst->buffer_cb + dst_chroma_index,
                dst->stride_cb << sub_pred,
                chroma_pu_width,
                chroma_pu_height >> sub_pred
            );

            ptr16 = (uint16_t *)ref_frame_pic_list0->buffer_cr + in_pos_x + in_pos_y * ref_frame_pic_list0->stride_cr;

            extract_8bit_data(
                ptr16,
                ref_frame_pic_list0->stride_cr << sub_pred,
                dst->buffer_cr + dst_chroma_index,
                dst->stride_cr << sub_pred,
                chroma_pu_width,
                chroma_pu_height >> sub_pred
            );
        }
    }
}

void estimate_uni_pred_interpolation_avc_chroma_ref10_bit(
    EbPictureBufferDesc   *ref_frame_pic_list0,
    uint32_t                 pos_x,
    uint32_t                 pos_y,
    uint32_t                 pu_width,
    uint32_t                 pu_height,
    EbPictureBufferDesc   *dst,
    uint32_t                 dst_luma_index,
    uint32_t                 dst_chroma_index,          //input parameter, please refer to the detailed explanation above.
    uint32_t                 component_mask,
    EbByte                  temp_buf,
    EbBool                   sub_pred)
{
    uint32_t   chroma_pu_width = pu_width >> 1;
    uint32_t   chroma_pu_height = pu_height >> 1;
    uint32_t   in_pos_x = (pos_x >> 3);
    uint32_t   in_pos_y = (pos_y >> 3);
    uint16_t  *ptr16 = (uint16_t *)ref_frame_pic_list0->buffer_cb + in_pos_x + in_pos_y * ref_frame_pic_list0->stride_cb;

    (void)temp_buf;
    (void)component_mask;
    (void)dst_luma_index;

    extract_8bit_data(
        ptr16,
        ref_frame_pic_list0->stride_cb << sub_pred,
        dst->buffer_cb + dst_chroma_index,
        dst->stride_cb << sub_pred,
        chroma_pu_width,
        chroma_pu_height >> sub_pred
    );

    ptr16 = (uint16_t *)ref_frame_pic_list0->buffer_cr + in_pos_x + in_pos_y * ref_frame_pic_list0->stride_cr;

    extract_8bit_data(
        ptr16,
        ref_frame_pic_list0->stride_cr << sub_pred,
        dst->buffer_cr + dst_chroma_index,
        dst->stride_cr << sub_pred,
        chroma_pu_width,
        chroma_pu_height >> sub_pred
    );
}
void estimate_bi_pred_interpolation_avc_chroma_ref10_bit(
    EbPictureBufferDesc   *ref_frame_pic_list0,
    EbPictureBufferDesc   *ref_frame_pic_list1,
    uint32_t                 ref_list0_pos_x,
    uint32_t                 ref_list0_pos_y,
    uint32_t                 ref_list1_pos_x,
    uint32_t                 ref_list1_pos_y,
    uint32_t                 pu_width,
    uint32_t                 pu_height,
    EbPictureBufferDesc   *bi_dst,
    uint32_t                 dst_luma_index,
    uint32_t                 dst_chroma_index,
    uint32_t                 component_mask,
    EbByte                  ref_list0_temp_dst,
    EbByte                  ref_list1_temp_dst,
    EbByte                  first_pass_if_temp_dst,
    EbBool                   sub_pred)
{
    uint32_t   chroma_pu_width = pu_width >> 1;
    uint32_t   chroma_pu_height = pu_height >> 1;

    (void)first_pass_if_temp_dst;
    (void)ref_list0_temp_dst;
    (void)ref_list1_temp_dst;
    (void)component_mask;
    (void)dst_luma_index;
    unpack_l0l1_avg(
        (uint16_t *)ref_frame_pic_list0->buffer_cb + (ref_list0_pos_x >> 3) + (ref_list0_pos_y >> 3)*ref_frame_pic_list0->stride_cb,
        ref_frame_pic_list0->stride_cb << sub_pred,
        (uint16_t *)ref_frame_pic_list1->buffer_cb + (ref_list1_pos_x >> 3) + (ref_list1_pos_y >> 3)*ref_frame_pic_list1->stride_cb,
        ref_frame_pic_list1->stride_cb << sub_pred,
        bi_dst->buffer_cb + dst_chroma_index,
        bi_dst->stride_cb << sub_pred,
        chroma_pu_width,
        chroma_pu_height >> sub_pred
    );

    unpack_l0l1_avg(
        (uint16_t *)ref_frame_pic_list0->buffer_cr + (ref_list0_pos_x >> 3) + (ref_list0_pos_y >> 3)*ref_frame_pic_list0->stride_cr,
        ref_frame_pic_list0->stride_cr << sub_pred,
        (uint16_t *)ref_frame_pic_list1->buffer_cr + (ref_list1_pos_x >> 3) + (ref_list1_pos_y >> 3)*ref_frame_pic_list1->stride_cr,
        ref_frame_pic_list1->stride_cr << sub_pred,
        bi_dst->buffer_cr + dst_chroma_index,
        bi_dst->stride_cr << sub_pred,
        chroma_pu_width,
        chroma_pu_height >> sub_pred
    );
}

void estimate_bi_pred_interpolation_avc_luma_ref10_bit(
    EbPictureBufferDesc   *ref_frame_pic_list0,
    EbPictureBufferDesc   *ref_frame_pic_list1,
    uint32_t                 ref_list0_pos_x,
    uint32_t                 ref_list0_pos_y,
    uint32_t                 ref_list1_pos_x,
    uint32_t                 ref_list1_pos_y,
    uint32_t                 pu_width,
    uint32_t                 pu_height,
    EbPictureBufferDesc   *bi_dst,
    uint32_t                 dst_luma_index,
    uint32_t                 dst_chroma_index,
    uint32_t                 component_mask,
    EbByte                  ref_list0_temp_dst,
    EbByte                  ref_list1_temp_dst,
    EbByte                  first_pass_if_temp_dst,
    EbBool                   sub_pred,
    EbBool                   sub_pred_chroma)
{
    uint32_t   chroma_pu_width = pu_width >> 1;
    uint32_t   chroma_pu_height = pu_height >> 1;

    //Luma
    if (component_mask & PICTURE_BUFFER_DESC_LUMA_MASK) {
        (void)first_pass_if_temp_dst;
        (void)ref_list0_temp_dst;
        (void)ref_list1_temp_dst;
        unpack_l0l1_avg_safe_sub(
            (uint16_t *)ref_frame_pic_list0->buffer_y + (ref_list0_pos_x >> 2) + (ref_list0_pos_y >> 2)*ref_frame_pic_list0->stride_y,
            ref_frame_pic_list0->stride_y << sub_pred,
            (uint16_t *)ref_frame_pic_list1->buffer_y + (ref_list1_pos_x >> 2) + (ref_list1_pos_y >> 2)*ref_frame_pic_list1->stride_y,
            ref_frame_pic_list1->stride_y << sub_pred,
            bi_dst->buffer_y + dst_luma_index,
            bi_dst->stride_y << sub_pred,
            pu_width,
            pu_height >> sub_pred,
            sub_pred
        );
    }

    //uni-prediction List0 chroma
    if (component_mask & PICTURE_BUFFER_DESC_CHROMA_MASK) {
        sub_pred = sub_pred_chroma;
        unpack_l0l1_avg(
            (uint16_t *)ref_frame_pic_list0->buffer_cb + (ref_list0_pos_x >> 3) + (ref_list0_pos_y >> 3)*ref_frame_pic_list0->stride_cb,
            ref_frame_pic_list0->stride_cb << sub_pred,
            (uint16_t *)ref_frame_pic_list1->buffer_cb + (ref_list1_pos_x >> 3) + (ref_list1_pos_y >> 3)*ref_frame_pic_list1->stride_cb,
            ref_frame_pic_list1->stride_cb << sub_pred,
            bi_dst->buffer_cb + dst_chroma_index,
            bi_dst->stride_cb << sub_pred,
            chroma_pu_width,
            chroma_pu_height >> sub_pred);

        unpack_l0l1_avg(
            (uint16_t *)ref_frame_pic_list0->buffer_cr + (ref_list0_pos_x >> 3) + (ref_list0_pos_y >> 3)*ref_frame_pic_list0->stride_cr,
            ref_frame_pic_list0->stride_cr << sub_pred,
            (uint16_t *)ref_frame_pic_list1->buffer_cr + (ref_list1_pos_x >> 3) + (ref_list1_pos_y >> 3)*ref_frame_pic_list1->stride_cr,
            ref_frame_pic_list1->stride_cr << sub_pred,
            bi_dst->buffer_cr + dst_chroma_index,
            bi_dst->stride_cr << sub_pred,
            chroma_pu_width,
            chroma_pu_height >> sub_pred);
    }
}

void uni_pred_i_free_ref8_bit(
    EbPictureBufferDesc   *ref_pic,
    uint32_t                 pos_x,
    uint32_t                 pos_y,
    uint32_t                 pu_width,
    uint32_t                 pu_height,
    EbPictureBufferDesc   *dst,
    uint32_t                 dst_luma_index,
    uint32_t                 dst_chroma_index,          //input parameter, please refer to the detailed explanation above.
    uint32_t                 component_mask,
    EbByte                  temp_buf,
    EbBool                   sub_sample_pred_flag,
    EbBool                   sub_sample_pred_flag_chroma)
{
    uint32_t   integ_pos_x;
    uint32_t   integ_pos_y;
    uint8_t    frac_posx, mapped_frac_posx;
    uint8_t    frac_posy, mapped_frac_posy;
    uint32_t   luma_stride;
    uint32_t   chroma_pu_width = pu_width >> 1;
    uint32_t   chroma_pu_height = pu_height >> 1;

    uint8_t    frac_pos;

    luma_stride = dst->stride_y;

    //luma
    //compute the luma fractional position
    integ_pos_x = (pos_x >> 2);
    integ_pos_y = (pos_y >> 2);

    frac_posx = pos_x & 0x03;
    frac_posy = pos_y & 0x03;

    //New Interpolation Mapping
    mapped_frac_posx = frac_posx;
    mapped_frac_posy = frac_posy;

    frac_pos = (mapped_frac_posx + (mapped_frac_posy << 2));
    mapped_frac_posx = frac_mapped_pos_tab_x[frac_pos];
    mapped_frac_posy = frac_mapped_pos_tab_y[frac_pos];
    integ_pos_x += integer_posoffset_tab_x[frac_pos];
    integ_pos_y += integer_posoffset_tab_y[frac_pos];
    if (component_mask & PICTURE_BUFFER_DESC_LUMA_MASK)
    {
        //doing the luma interpolation
        avc_style_luma_interpolation_filter(
            ref_pic->buffer_y + integ_pos_x + integ_pos_y * ref_pic->stride_y, ref_pic->stride_y,
            dst->buffer_y + dst_luma_index, luma_stride,
            pu_width, pu_height,
            temp_buf,
            sub_sample_pred_flag,
            mapped_frac_posx ? mapped_frac_posx : mapped_frac_posy,
            mapped_frac_posx + (mapped_frac_posy << 2));
    }
    //chroma
    if (component_mask & PICTURE_BUFFER_DESC_CHROMA_MASK) {
        //compute the chroma fractional position
        integ_pos_x = (pos_x >> 3);
        integ_pos_y = (pos_y >> 3);
        frac_posx = pos_x & 0x07;
        frac_posy = pos_y & 0x07;

        mapped_frac_posx = 0;
        if (frac_posx > 4)
            integ_pos_x++;
        mapped_frac_posy = 0;
        if (frac_posy > 4)
            integ_pos_y++;

        // Note: chroma_pu_width equals 4 is only supported in Intrinsic
       //       for integer positions ( mapped_frac_posx + (mapped_frac_posy << 3) equals 0 )
       //doing the chroma Cb interpolation
        avc_style_luma_interpolation_filter(
            ref_pic->buffer_cb + integ_pos_x + integ_pos_y * ref_pic->stride_cb,
            ref_pic->stride_cb,
            dst->buffer_cb + dst_chroma_index,
            dst->stride_cb,
            chroma_pu_width,
            chroma_pu_height,
            temp_buf,
            sub_sample_pred_flag_chroma,
            mapped_frac_posx ? mapped_frac_posx : mapped_frac_posy,
            mapped_frac_posx + (mapped_frac_posy << 3));

        //doing the chroma Cr interpolation
        avc_style_luma_interpolation_filter(
            ref_pic->buffer_cr + integ_pos_x + integ_pos_y * ref_pic->stride_cr,
            ref_pic->stride_cr,
            dst->buffer_cr + dst_chroma_index,
            dst->stride_cr,
            chroma_pu_width,
            chroma_pu_height,
            temp_buf,
            sub_sample_pred_flag_chroma,
            mapped_frac_posx ? mapped_frac_posx : mapped_frac_posy,
            mapped_frac_posx + (mapped_frac_posy << 3));
    }
}

void bi_pred_i_free_ref8_bit(
    EbPictureBufferDesc   *ref_pic_list0,
    EbPictureBufferDesc   *ref_pic_list1,
    uint32_t                 ref_list0_pos_x,
    uint32_t                 ref_list0_pos_y,
    uint32_t                 ref_list1_pos_x,
    uint32_t                 ref_list1_pos_y,
    uint32_t                 pu_width,
    uint32_t                 pu_height,
    EbPictureBufferDesc   *bi_dst,
    uint32_t                 dst_luma_index,
    uint32_t                 dst_chroma_index,
    uint32_t                 component_mask,
    EbByte                  ref_list0_temp_dst,
    EbByte                  ref_list1_temp_dst,
    EbByte                  first_pass_if_temp_dst,
    EbBool                   sub_sample_pred_flag,
    EbBool                   sub_sample_pred_flag_chroma) //needs to be connected
{
    uint32_t   integ_pos_x;
    uint32_t   integ_pos_y;
    uint8_t    frac_posx, mapped_frac_posx;
    uint8_t    frac_posy, mapped_frac_posy;
    uint8_t    frac_pos;

    uint32_t   luma_stride, ref_luma_stride;

    uint32_t   chroma_pu_width = pu_width >> 1;
    uint32_t   chroma_pu_height = pu_height >> 1;
    luma_stride = bi_dst->stride_y;
    ref_luma_stride = ref_pic_list0->stride_y;
    uint8_t shift = sub_sample_pred_flag ? 1 : 0;

    //Luma
    if (component_mask & PICTURE_BUFFER_DESC_LUMA_MASK) {
        //uni-prediction List0 luma
        //compute the luma fractional position
        integ_pos_x = (ref_list0_pos_x >> 2);
        integ_pos_y = (ref_list0_pos_y >> 2);

        frac_posx = ref_list0_pos_x & 0x03;
        frac_posy = ref_list0_pos_y & 0x03;

        mapped_frac_posx = frac_posx;
        mapped_frac_posy = frac_posy;

        frac_pos = (mapped_frac_posx + (mapped_frac_posy << 2));
        mapped_frac_posx = frac_mapped_pos_tab_x[frac_pos];
        mapped_frac_posy = frac_mapped_pos_tab_y[frac_pos];
        integ_pos_x += integer_posoffset_tab_x[frac_pos];
        integ_pos_y += integer_posoffset_tab_y[frac_pos];

        avc_style_luma_interpolation_filter(
            ref_pic_list0->buffer_y + integ_pos_x + integ_pos_y * ref_luma_stride, ref_luma_stride,
            ref_list0_temp_dst, pu_width,
            pu_width, pu_height,
            first_pass_if_temp_dst,
            sub_sample_pred_flag,
            mapped_frac_posx ? mapped_frac_posx : mapped_frac_posy,
            mapped_frac_posx + (mapped_frac_posy << 2));

        //uni-prediction List1 luma
        //compute the luma fractional position
        integ_pos_x = (ref_list1_pos_x >> 2);
        integ_pos_y = (ref_list1_pos_y >> 2);

        frac_posx = ref_list1_pos_x & 0x03;
        frac_posy = ref_list1_pos_y & 0x03;

        mapped_frac_posx = frac_posx;
        mapped_frac_posy = frac_posy;

        frac_pos = (mapped_frac_posx + (mapped_frac_posy << 2));
        mapped_frac_posx = frac_mapped_pos_tab_x[frac_pos];
        mapped_frac_posy = frac_mapped_pos_tab_y[frac_pos];
        integ_pos_x += integer_posoffset_tab_x[frac_pos];
        integ_pos_y += integer_posoffset_tab_y[frac_pos];

        //doing the luma interpolation
        avc_style_luma_interpolation_filter(
            ref_pic_list1->buffer_y + integ_pos_x + integ_pos_y * ref_luma_stride, ref_luma_stride,
            ref_list1_temp_dst, pu_width,
            pu_width, pu_height,
            first_pass_if_temp_dst,
            sub_sample_pred_flag,
            mapped_frac_posx ? mapped_frac_posx : mapped_frac_posy,
            mapped_frac_posx + (mapped_frac_posy << 2));

        // bi-pred luma
        picture_average_kernel(ref_list0_temp_dst, pu_width << sub_sample_pred_flag, ref_list1_temp_dst, pu_width << sub_sample_pred_flag, bi_dst->buffer_y + dst_luma_index, luma_stride << sub_sample_pred_flag, pu_width, pu_height >> sub_sample_pred_flag);
        if (sub_sample_pred_flag)
            picture_average_kernel1_line(ref_list0_temp_dst + (pu_height - 1)*pu_width, ref_list1_temp_dst + (pu_height - 1)*pu_width, bi_dst->buffer_y + dst_luma_index + (pu_height - 1)*luma_stride, pu_width);
    }

    //uni-prediction List0 chroma
    if (component_mask & PICTURE_BUFFER_DESC_CHROMA_MASK) {
        shift = sub_sample_pred_flag_chroma ? 1 : 0;

        // bi-pred chroma  Cb
        // Note: chroma_pu_width equals 4 is only supported in Intrinsic
        //       for integer positions ( mapped_frac_posx + (mapped_frac_posy << 3) equals 0 )
        //doing the chroma Cb interpolation list 0
        //compute the chroma fractional position
        integ_pos_x = (ref_list0_pos_x >> 3);
        integ_pos_y = (ref_list0_pos_y >> 3);
        frac_posx = ref_list0_pos_x & 0x07;
        frac_posy = ref_list0_pos_y & 0x07;

        mapped_frac_posx = 0;
        if (frac_posx > 4)
            integ_pos_x++;
        mapped_frac_posy = 0;
        if (frac_posy > 4)
            integ_pos_y++;

        avc_style_luma_interpolation_filter(
            ref_pic_list0->buffer_cb + integ_pos_x + integ_pos_y * ref_pic_list0->stride_cb,
            ref_pic_list0->stride_cb,
            ref_list0_temp_dst,
            chroma_pu_width,
            chroma_pu_width,
            chroma_pu_height,
            first_pass_if_temp_dst,
            sub_sample_pred_flag_chroma,
            mapped_frac_posx ? mapped_frac_posx : mapped_frac_posy,
            mapped_frac_posx + (mapped_frac_posy << 3));

        //doing the chroma Cb interpolation list 1
        integ_pos_x = (ref_list1_pos_x >> 3);
        integ_pos_y = (ref_list1_pos_y >> 3);
        frac_posx = ref_list1_pos_x & 0x07;
        frac_posy = ref_list1_pos_y & 0x07;

        mapped_frac_posx = 0;
        if (frac_posx > 4)
            integ_pos_x++;
        mapped_frac_posy = 0;
        if (frac_posy > 4)
            integ_pos_y++;

        avc_style_luma_interpolation_filter(
            ref_pic_list1->buffer_cb + integ_pos_x + integ_pos_y * ref_pic_list1->stride_cb,
            ref_pic_list1->stride_cb,
            ref_list1_temp_dst,
            chroma_pu_width,
            chroma_pu_width,
            chroma_pu_height,
            first_pass_if_temp_dst,
            sub_sample_pred_flag_chroma,
            mapped_frac_posx ? mapped_frac_posx : mapped_frac_posy,
            mapped_frac_posx + (mapped_frac_posy << 3));

        // bi-pred Chroma Cb
        picture_average_kernel(
            ref_list0_temp_dst,
            chroma_pu_width << shift,
            ref_list1_temp_dst,
            chroma_pu_width << shift,
            bi_dst->buffer_cb + dst_chroma_index,
            bi_dst->stride_cb << shift,
            chroma_pu_width,
            chroma_pu_height >> shift);

        // bi-pred chroma  Cr
        // Note: chroma_pu_width equals 4 is only supported in Intrinsic
        //       for integer positions ( mapped_frac_posx + (mapped_frac_posy << 3) equals 0 )
        //doing the chroma Cb interpolation list 0
        //compute the chroma fractional position
        integ_pos_x = (ref_list0_pos_x >> 3);
        integ_pos_y = (ref_list0_pos_y >> 3);
        frac_posx = ref_list0_pos_x & 0x07;
        frac_posy = ref_list0_pos_y & 0x07;

        mapped_frac_posx = 0;
        if (frac_posx > 4)
            integ_pos_x++;
        mapped_frac_posy = 0;
        if (frac_posy > 4)
            integ_pos_y++;

        avc_style_luma_interpolation_filter(
            ref_pic_list0->buffer_cr + integ_pos_x + integ_pos_y * ref_pic_list0->stride_cr,
            ref_pic_list0->stride_cr,
            ref_list0_temp_dst,
            chroma_pu_width,
            chroma_pu_width,
            chroma_pu_height,
            first_pass_if_temp_dst,
            sub_sample_pred_flag_chroma,
            mapped_frac_posx ? mapped_frac_posx : mapped_frac_posy,
            mapped_frac_posx + (mapped_frac_posy << 3));

        //doing the chroma Cb interpolation list 1
        integ_pos_x = (ref_list1_pos_x >> 3);
        integ_pos_y = (ref_list1_pos_y >> 3);
        frac_posx = ref_list1_pos_x & 0x07;
        frac_posy = ref_list1_pos_y & 0x07;

        mapped_frac_posx = 0;
        if (frac_posx > 4)
            integ_pos_x++;
        mapped_frac_posy = 0;
        if (frac_posy > 4)
            integ_pos_y++;

        avc_style_luma_interpolation_filter(
            ref_pic_list1->buffer_cr + integ_pos_x + integ_pos_y * ref_pic_list1->stride_cr,
            ref_pic_list1->stride_cr,
            ref_list1_temp_dst,
            chroma_pu_width,
            chroma_pu_width,
            chroma_pu_height,
            first_pass_if_temp_dst,
            sub_sample_pred_flag_chroma,
            mapped_frac_posx ? mapped_frac_posx : mapped_frac_posy,
            mapped_frac_posx + (mapped_frac_posy << 3));

        // bi-pred Chroma Cr
        picture_average_kernel(
            ref_list0_temp_dst,
            chroma_pu_width << shift,
            ref_list1_temp_dst,
            chroma_pu_width << shift,
            bi_dst->buffer_cr + dst_chroma_index,
            bi_dst->stride_cr << shift,
            chroma_pu_width,
            chroma_pu_height >> shift);
    }
}

void avc_style_luma_interpolation_filter_horizontal_c(
    EbByte                 ref_pic,
    uint32_t               src_stride,
    EbByte                 dst,
    uint32_t               dst_stride,
    uint32_t               pu_width,
    uint32_t               pu_height,
    EbByte                 temp_buf,
    uint32_t               frac_pos)
{
    const int8_t  *if_coeff = avc_style_luma_if_coeff[frac_pos];
    const int32_t  if_init_pos_offset = -1;
    const uint8_t  if_shift = 5;
    const int16_t  if_offset = POW2(if_shift - 1);
    uint32_t       x, y;
    (void)temp_buf;

    ref_pic += if_init_pos_offset;
    for (y = 0; y < pu_height; y++) {
        for (x = 0; x < pu_width; x++) {
            dst[x] = (uint8_t)CLIP3(0, MAX_SAMPLE_VALUE,
                                    (ref_pic[x] * if_coeff[0] +
                                     ref_pic[x + 1] * if_coeff[1] +
                                     ref_pic[x + 2] * if_coeff[2] +
                                     ref_pic[x + 3] * if_coeff[3] +
                                     if_offset) >> if_shift);
        }
        ref_pic += src_stride;
        dst += dst_stride;
    }
}

void avc_style_luma_interpolation_filter_vertical_c(
    EbByte                 ref_pic,
    uint32_t               src_stride,
    EbByte                 dst,
    uint32_t               dst_stride,
    uint32_t               pu_width,
    uint32_t               pu_height,
    EbByte                 temp_buf,
    uint32_t               frac_pos)
{
    const int8_t  *if_coeff = avc_style_luma_if_coeff[frac_pos];
    const int32_t  if_stride = src_stride;
    const int32_t  if_init_pos_offset = -(int32_t)src_stride;
    const uint8_t  if_shift = 5;
    const int16_t  if_offset = POW2(if_shift - 1);
    uint32_t       x, y;
    (void)temp_buf;

    ref_pic += if_init_pos_offset;
    for (y = 0; y < pu_height; y++) {
        for (x = 0; x < pu_width; x++) {
            dst[x] = (uint8_t)CLIP3(0, MAX_SAMPLE_VALUE,
                                    (ref_pic[x] * if_coeff[0] +
                                     ref_pic[x + if_stride] * if_coeff[1] +
                                     ref_pic[x + 2 * if_stride] * if_coeff[2] +
                                     ref_pic[x + 3 * if_stride] * if_coeff[3] +
                                     if_offset) >> if_shift);
        }
        ref_pic += src_stride;
        dst += dst_stride;
    }
}

void avc_style_luma_interpolation_filter_pose_c(
    EbByte                 ref_pic,
    uint32_t               src_stride,
    EbByte                 dst,
    uint32_t               dst_stride,
    uint32_t               pu_width,
    uint32_t               pu_height,
    EbByte                 temp_buf,
    uint32_t               frac_pos)
{
    uint32_t temp_buf_size = pu_width * pu_height;
    (void)frac_pos;
    avc_style_luma_interpolation_filter_horizontal_c(
        ref_pic, src_stride, temp_buf, pu_width, pu_width, pu_height, 0, 2);
    avc_style_luma_interpolation_filter_vertical_c(
        ref_pic, src_stride, temp_buf + temp_buf_size, pu_width, pu_width, pu_height, 0,
        2);
    picture_average_kernel_c(temp_buf, pu_width, temp_buf + temp_buf_size, pu_width, dst,
                         dst_stride, pu_width, pu_height);
}

void avc_style_luma_interpolation_filter_posf_c(
    EbByte                 ref_pic,
    uint32_t               src_stride,
    EbByte                 dst,
    uint32_t               dst_stride,
    uint32_t               pu_width,
    uint32_t               pu_height,
    EbByte                 temp_buf,
    uint32_t               frac_pos)
{
    uint32_t temp_buf_size = pu_width * pu_height;
    (void)frac_pos;
    avc_style_luma_interpolation_filter_horizontal_c(
        ref_pic, src_stride, temp_buf, pu_width, pu_width, pu_height, 0, 2);
    avc_style_luma_interpolation_filter_posj_c(
        ref_pic, src_stride, temp_buf + temp_buf_size, pu_width, pu_width, pu_height,
        temp_buf + 2 * temp_buf_size, 2);
    picture_average_kernel_c(temp_buf, pu_width, temp_buf + temp_buf_size, pu_width, dst,
                         dst_stride, pu_width, pu_height);
}

void avc_style_luma_interpolation_filter_posg_c(
    EbByte                 ref_pic,
    uint32_t               src_stride,
    EbByte                 dst,
    uint32_t               dst_stride,
    uint32_t               pu_width,
    uint32_t               pu_height,
    EbByte                 temp_buf,
    uint32_t               frac_pos)
{
    uint32_t temp_buf_size = pu_width * pu_height;
    (void)frac_pos;
    avc_style_luma_interpolation_filter_horizontal_c(
        ref_pic, src_stride, temp_buf, pu_width, pu_width, pu_height, 0, 2);
    avc_style_luma_interpolation_filter_vertical_c(
        ref_pic + 1, src_stride, temp_buf + temp_buf_size, pu_width, pu_width, pu_height,
        0, 2);
    picture_average_kernel_c(temp_buf, pu_width, temp_buf + temp_buf_size, pu_width, dst,
                         dst_stride, pu_width, pu_height);
}

void avc_style_luma_interpolation_filter_posi_c(
    EbByte                 ref_pic,
    uint32_t               src_stride,
    EbByte                 dst,
    uint32_t               dst_stride,
    uint32_t               pu_width,
    uint32_t               pu_height,
    EbByte                 temp_buf,
    uint32_t               frac_pos)
{
    uint32_t temp_buf_size = pu_width * pu_height;
    (void)frac_pos;
    avc_style_luma_interpolation_filter_vertical_c(ref_pic, src_stride, temp_buf, pu_width,
                                              pu_width, pu_height, 0, 2);
    avc_style_luma_interpolation_filter_posj_c(
        ref_pic, src_stride, temp_buf + temp_buf_size, pu_width, pu_width, pu_height,
        temp_buf + 2 * temp_buf_size, 2);
    picture_average_kernel_c(temp_buf, pu_width, temp_buf + temp_buf_size, pu_width, dst,
                         dst_stride, pu_width, pu_height);
}

void avc_style_luma_interpolation_filter_posj_c(
    EbByte                 ref_pic,
    uint32_t               src_stride,
    EbByte                 dst,
    uint32_t               dst_stride,
    uint32_t               pu_width,
    uint32_t               pu_height,
    EbByte                 temp_buf,
    uint32_t               frac_pos)
{
    (void)frac_pos;
    avc_style_luma_interpolation_filter_horizontal_c(
        ref_pic - src_stride, src_stride, temp_buf, pu_width, pu_width, (pu_height + 4),
        0, 2);
    avc_style_luma_interpolation_filter_vertical_c(temp_buf + pu_width, pu_width, dst,
                                              dst_stride, pu_width, pu_height, 0, 2);
}

void avc_style_luma_interpolation_filter_posk_c(
    EbByte                 ref_pic,
    uint32_t               src_stride,
    EbByte                 dst,
    uint32_t               dst_stride,
    uint32_t               pu_width,
    uint32_t               pu_height,
    EbByte                 temp_buf,
    uint32_t               frac_pos)
{
    uint32_t temp_buf_size = pu_width * pu_height;
    (void)frac_pos;
    avc_style_luma_interpolation_filter_vertical_c(ref_pic + 1, src_stride, temp_buf,
                                              pu_width, pu_width, pu_height, 0, 2);
    avc_style_luma_interpolation_filter_posj_c(
        ref_pic, src_stride, temp_buf + temp_buf_size, pu_width, pu_width, pu_height,
        temp_buf + 2 * temp_buf_size, 2);
    picture_average_kernel_c(temp_buf, pu_width, temp_buf + temp_buf_size, pu_width, dst,
                         dst_stride, pu_width, pu_height);
}

void avc_style_luma_interpolation_filter_posp_c(
    EbByte                 ref_pic,
    uint32_t               src_stride,
    EbByte                 dst,
    uint32_t               dst_stride,
    uint32_t               pu_width,
    uint32_t               pu_height,
    EbByte                 temp_buf,
    uint32_t               frac_pos)
{
    uint32_t temp_buf_size = pu_width * pu_height;
    (void)frac_pos;
    avc_style_luma_interpolation_filter_vertical_c(ref_pic, src_stride, temp_buf, pu_width,
                                              pu_width, pu_height, 0, 2);
    avc_style_luma_interpolation_filter_horizontal_c(
        ref_pic + src_stride, src_stride, temp_buf + temp_buf_size, pu_width, pu_width,
        pu_height, 0, 2);
    picture_average_kernel_c(temp_buf, pu_width, temp_buf + temp_buf_size, pu_width, dst,
                         dst_stride, pu_width, pu_height);
}

void avc_style_luma_interpolation_filter_posq_c(
    EbByte                 ref_pic,
    uint32_t               src_stride,
    EbByte                 dst,
    uint32_t               dst_stride,
    uint32_t               pu_width,
    uint32_t               pu_height,
    EbByte                 temp_buf,
    uint32_t               frac_pos)
{
    uint32_t temp_buf_size = pu_width * pu_height;
    (void)frac_pos;
    avc_style_luma_interpolation_filter_horizontal_c(
        ref_pic + src_stride, src_stride, temp_buf, pu_width, pu_width, pu_height, 0, 2);
    avc_style_luma_interpolation_filter_posj_c(
        ref_pic, src_stride, temp_buf + temp_buf_size, pu_width, pu_width, pu_height,
        temp_buf + 2 * temp_buf_size, 2);
    picture_average_kernel_c(temp_buf, pu_width, temp_buf + temp_buf_size, pu_width, dst,
                         dst_stride, pu_width, pu_height);
}

void avc_style_luma_interpolation_filter_posr_c(
    EbByte                 ref_pic,
    uint32_t               src_stride,
    EbByte                 dst,
    uint32_t               dst_stride,
    uint32_t               pu_width,
    uint32_t               pu_height,
    EbByte                 temp_buf,
    uint32_t               frac_pos)
{
    uint32_t temp_buf_size = pu_width * pu_height;
    (void)frac_pos;
    avc_style_luma_interpolation_filter_vertical_c(ref_pic + 1, src_stride, temp_buf,
                                              pu_width, pu_width, pu_height, 0, 2);
    avc_style_luma_interpolation_filter_horizontal_c(
        ref_pic + src_stride, src_stride, temp_buf + temp_buf_size, pu_width, pu_width,
        pu_height, 0, 2);
    picture_average_kernel_c(temp_buf, pu_width, temp_buf + temp_buf_size, pu_width, dst,
                         dst_stride, pu_width, pu_height);
}

void avc_style_luma_interpolation_filter_helper_c(
    EbByte ref_pic,
    uint32_t src_stride,
    EbByte dst,
    uint32_t dst_stride,
    uint32_t pu_width,
    uint32_t pu_height,
    EbByte temp_buf,
    EbBool skip,
    uint32_t frac_pos,
    uint8_t fractional_position)
{
    /* Code with 'skip' true are not used. Cleanup should remove 'skip' parameter. */
    (void)skip;
    assert(!skip);

    switch (fractional_position) {
    case 0:
        avc_style_copy_c(ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos); break;
    case 1:
        avc_style_luma_interpolation_filter_horizontal_c(ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos); break;
    case 2:
        avc_style_luma_interpolation_filter_horizontal_c(ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos); break;
    case 3:
        avc_style_luma_interpolation_filter_horizontal_c(ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos); break;
    case 4:
        avc_style_luma_interpolation_filter_vertical_c(ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos); break;
    case 5:
        avc_style_luma_interpolation_filter_pose_c(ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos); break;
    case 6:
        avc_style_luma_interpolation_filter_posf_c(ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos); break;
    case 7:
        avc_style_luma_interpolation_filter_posg_c(ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos); break;
    case 8:
        avc_style_luma_interpolation_filter_vertical_c(ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos); break;
    case 9:
        avc_style_luma_interpolation_filter_posi_c(ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos); break;
    case 10:
        avc_style_luma_interpolation_filter_posj_c(ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos); break;
    case 11:
        avc_style_luma_interpolation_filter_posk_c(ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos); break;
    case 12:
        avc_style_luma_interpolation_filter_vertical_c(ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos); break;
    case 13:
        avc_style_luma_interpolation_filter_posp_c(ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos); break;
    case 14:
        avc_style_luma_interpolation_filter_posq_c(ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos); break;
    case 15:
        avc_style_luma_interpolation_filter_posr_c(ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, temp_buf, frac_pos); break;
    default:
        assert(0);
    }
}
