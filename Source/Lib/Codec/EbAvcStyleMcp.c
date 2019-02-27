/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbAvcStyleMcp.h"
#include "EbPictureOperators.h"

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




void estimate_uni_pred_interpolation_unpacked_avc_style(
    EbPictureBufferDesc_t   *ref_pic,
    uint32_t                 pos_x,
    uint32_t                 pos_y,
    uint32_t                 pu_width,
    uint32_t                 pu_height,
    EbPictureBufferDesc_t   *dst,
    uint32_t                 dst_luma_index,
    uint32_t                 dst_chroma_index,          //input parameter, please refer to the detailed explanation above.
    uint32_t                 component_mask,
    EbByte                  temp_buf,
    EbBool                   sub_sample_pred_flag,
    EbAsm                   asm_type)
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
    avc_style_uni_pred_luma_if_function_ptr_array[asm_type][mapped_frac_posx + (mapped_frac_posy << 2)](
        ref_pic->buffer_y + 2 + 2 * ref_pic->stride_y, ref_pic->stride_y,
        dst->buffer_y + dst_luma_index, luma_stride,
        pu_width, pu_height,
        temp_buf,
        sub_sample_pred_flag,
        mapped_frac_posx ? mapped_frac_posx : mapped_frac_posy);
}

/*******************************************************************************
* Requirement: pu_width      = 8, 16, 24, 32, 48 or 64
* Requirement: pu_height % 2 = 0
* Requirement: skip         = 0 or 1
* Requirement (x86 only): temp_buf % 16 = 0
* Requirement (x86 only): (dst->buffer_y  + dst_luma_index  ) % 16 = 0 when pu_width %16 = 0
* Requirement (x86 only): (dst->bufferCb + dst_chroma_index) % 16 = 0 when pu_width %32 = 0
* Requirement (x86 only): (dst->bufferCr + dst_chroma_index) % 16 = 0 when pu_width %32 = 0
* Requirement (x86 only): dst->stride_y   % 16 = 0 when pu_width %16 = 0
* Requirement (x86 only): dst->chromaStride % 16 = 0 when pu_width %32 = 0
*******************************************************************************/
void estimate_bi_pred_interpolation_unpacked_avc_style(
    EbPictureBufferDesc_t   *ref_pic_list0,
    EbPictureBufferDesc_t   *ref_pic_list1,
    uint32_t                 ref_list0_pos_x,
    uint32_t                 ref_list0_pos_y,
    uint32_t                 ref_list1_pos_x,
    uint32_t                 ref_list1_pos_y,
    uint32_t                 pu_width,
    uint32_t                 pu_height,
    EbPictureBufferDesc_t   *bi_dst,
    uint32_t                 dst_luma_index,
    uint32_t                 dst_chroma_index,
    uint32_t                 component_mask,
    EbByte                  ref_list0_temp_dst,
    EbByte                  ref_list1_temp_dst,
    EbByte                  first_pass_if_temp_dst,
    EbBool                   sub_sample_pred_flag,
    EbAsm                    asm_type)
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


    avc_style_uni_pred_luma_if_function_ptr_array[asm_type][mapped_frac_posx + (mapped_frac_posy << 2)](
        ref_pic_list0->buffer_y + integ_pos_x + integ_pos_y * ref_luma_stride, ref_luma_stride,
        ref_list0_temp_dst, pu_width,
        pu_width, pu_height,
        first_pass_if_temp_dst,
        sub_sample_pred_flag,
        mapped_frac_posx ? mapped_frac_posx : mapped_frac_posy);

    frac_posx = ref_list1_pos_x & 0x03;
    frac_posy = ref_list1_pos_y & 0x03;

    mapped_frac_posx = frac_posx;
    mapped_frac_posy = frac_posy;

    //doing the luma interpolation
    avc_style_uni_pred_luma_if_function_ptr_array[asm_type][mapped_frac_posx + (mapped_frac_posy << 2)](
        ref_pic_list1->buffer_y + integ_pos_x + integ_pos_y * ref_luma_stride, ref_luma_stride,
        ref_list1_temp_dst, pu_width,
        pu_width, pu_height,
        first_pass_if_temp_dst,
        sub_sample_pred_flag,
        mapped_frac_posx ? mapped_frac_posx : mapped_frac_posy);

    // bi-pred luma
    picture_average_array[asm_type](ref_list0_temp_dst, pu_width << sub_sample_pred_flag, ref_list1_temp_dst, pu_width << sub_sample_pred_flag, bi_dst->buffer_y + dst_luma_index, luma_stride << sub_sample_pred_flag, pu_width, pu_height >> sub_sample_pred_flag);
    if (sub_sample_pred_flag) {
        picture_average1_line_array[asm_type](ref_list0_temp_dst + (pu_height - 1)*pu_width, ref_list1_temp_dst + (pu_height - 1)*pu_width, bi_dst->buffer_y + dst_luma_index + (pu_height - 1)*luma_stride, pu_width);
    }
}


/*******************************************************************************
* Requirement: pu_width      = 8, 16, 24, 32, 48 or 64
* Requirement: pu_height % 2 = 0
* Requirement: skip         = 0 or 1
* Requirement (x86 only): temp_buf % 16 = 0
* Requirement (x86 only): (dst->buffer_y  + dst_luma_index  ) % 16 = 0 when pu_width %16 = 0
* Requirement (x86 only): (dst->bufferCb + dst_chroma_index) % 16 = 0 when pu_width %32 = 0
* Requirement (x86 only): (dst->bufferCr + dst_chroma_index) % 16 = 0 when pu_width %32 = 0
* Requirement (x86 only): dst->stride_y   % 16 = 0 when pu_width %16 = 0
* Requirement (x86 only): dst->chromaStride % 16 = 0 when pu_width %32 = 0
*******************************************************************************/
void estimate_uni_pred_interpolation_avc_luma(
    EbPictureBufferDesc_t   *ref_pic,
    uint32_t                 pos_x,
    uint32_t                 pos_y,
    uint32_t                 pu_width,
    uint32_t                 pu_height,
    EbPictureBufferDesc_t   *dst,
    uint32_t                 dst_luma_index,
    uint32_t                 dst_chroma_index,          //input parameter, please refer to the detailed explanation above.
    uint32_t                 component_mask,
    EbByte                  temp_buf,
    EbBool                   sub_sample_pred_flag,
    EbAsm                    asm_type)
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
        avc_style_uni_pred_luma_if_function_ptr_array[asm_type][mapped_frac_posx + (mapped_frac_posy << 2)](
            ref_pic->buffer_y + integ_pos_x + integ_pos_y * ref_pic->stride_y, ref_pic->stride_y,
            dst->buffer_y + dst_luma_index, luma_stride,
            pu_width, pu_height,
            temp_buf,
            sub_sample_pred_flag,
            mapped_frac_posx ? mapped_frac_posx : mapped_frac_posy);

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
        avc_style_uni_pred_luma_if_function_ptr_array[asm_type][mapped_frac_posx + (mapped_frac_posy << 3)](
            ref_pic->bufferCb + integ_pos_x + integ_pos_y * ref_pic->strideCb,
            ref_pic->strideCb,
            dst->bufferCb + dst_chroma_index,
            dst->strideCb,
            chroma_pu_width,
            chroma_pu_height,
            temp_buf,
            sub_sample_pred_flag,
            mapped_frac_posx ? mapped_frac_posx : mapped_frac_posy);

        //doing the chroma Cr interpolation
        avc_style_uni_pred_luma_if_function_ptr_array[asm_type][mapped_frac_posx + (mapped_frac_posy << 3)](
            ref_pic->bufferCr + integ_pos_x + integ_pos_y * ref_pic->strideCr,
            ref_pic->strideCr,
            dst->bufferCr + dst_chroma_index,
            dst->strideCr,
            chroma_pu_width,
            chroma_pu_height,
            temp_buf,
            sub_sample_pred_flag,
            mapped_frac_posx ? mapped_frac_posx : mapped_frac_posy);

    }
}

/*******************************************************************************
 * Requirement: pu_width      = 8, 16, 24, 32, 48 or 64
 * Requirement: pu_height % 2 = 0
 * Requirement: skip         = 0 or 1
 * Requirement (x86 only): temp_buf % 16 = 0
 * Requirement (x86 only): (dst->buffer_y  + dst_luma_index  ) % 16 = 0 when pu_width %16 = 0
 * Requirement (x86 only): (dst->bufferCb + dst_chroma_index) % 16 = 0 when pu_width %32 = 0
 * Requirement (x86 only): (dst->bufferCr + dst_chroma_index) % 16 = 0 when pu_width %32 = 0
 * Requirement (x86 only): dst->stride_y   % 16 = 0 when pu_width %16 = 0
 * Requirement (x86 only): dst->chromaStride % 16 = 0 when pu_width %32 = 0
*******************************************************************************/
void estimate_bi_pred_interpolation_avc_luma(
    EbPictureBufferDesc_t   *ref_pic_list0,
    EbPictureBufferDesc_t   *ref_pic_list1,
    uint32_t                 ref_list0_pos_x,
    uint32_t                 ref_list0_pos_y,
    uint32_t                 ref_list1_pos_x,
    uint32_t                 ref_list1_pos_y,
    uint32_t                 pu_width,
    uint32_t                 pu_height,
    EbPictureBufferDesc_t   *bi_dst,
    uint32_t                 dst_luma_index,
    uint32_t                 dst_chroma_index,
    uint32_t                 component_mask,
    EbByte                  ref_list0_temp_dst,
    EbByte                  ref_list1_temp_dst,
    EbByte                  first_pass_if_temp_dst,
    EbBool                   sub_sample_pred_flag,
    EbAsm                    asm_type)
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

        avc_style_uni_pred_luma_if_function_ptr_array[asm_type][mapped_frac_posx + (mapped_frac_posy << 2)](
            ref_pic_list0->buffer_y + integ_pos_x + integ_pos_y * ref_luma_stride, ref_luma_stride,
            ref_list0_temp_dst, pu_width,
            pu_width, pu_height,
            first_pass_if_temp_dst,
            sub_sample_pred_flag,
            mapped_frac_posx ? mapped_frac_posx : mapped_frac_posy);

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
        avc_style_uni_pred_luma_if_function_ptr_array[asm_type][mapped_frac_posx + (mapped_frac_posy << 2)](
            ref_pic_list1->buffer_y + integ_pos_x + integ_pos_y * ref_luma_stride, ref_luma_stride,
            ref_list1_temp_dst, pu_width,
            pu_width, pu_height,
            first_pass_if_temp_dst,
            sub_sample_pred_flag,
            mapped_frac_posx ? mapped_frac_posx : mapped_frac_posy);


        // bi-pred luma
        picture_average_array[asm_type](ref_list0_temp_dst, pu_width << sub_sample_pred_flag, ref_list1_temp_dst, pu_width << sub_sample_pred_flag, bi_dst->buffer_y + dst_luma_index, luma_stride << sub_sample_pred_flag, pu_width, pu_height >> sub_sample_pred_flag);
        if (sub_sample_pred_flag) {
            picture_average1_line_array[asm_type](ref_list0_temp_dst + (pu_height - 1)*pu_width, ref_list1_temp_dst + (pu_height - 1)*pu_width, bi_dst->buffer_y + dst_luma_index + (pu_height - 1)*luma_stride, pu_width);
        }
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

        avc_style_uni_pred_luma_if_function_ptr_array[asm_type][mapped_frac_posx + (mapped_frac_posy << 3)](
            ref_pic_list0->bufferCb + integ_pos_x + integ_pos_y * ref_pic_list0->strideCb,
            ref_pic_list0->strideCb,
            ref_list0_temp_dst,
            chroma_pu_width,
            chroma_pu_width,
            chroma_pu_height,
            first_pass_if_temp_dst,
            sub_sample_pred_flag,
            mapped_frac_posx ? mapped_frac_posx : mapped_frac_posy);


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

        avc_style_uni_pred_luma_if_function_ptr_array[asm_type][mapped_frac_posx + (mapped_frac_posy << 3)](
            ref_pic_list1->bufferCb + integ_pos_x + integ_pos_y * ref_pic_list1->strideCb,
            ref_pic_list1->strideCb,
            ref_list1_temp_dst,
            chroma_pu_width,
            chroma_pu_width,
            chroma_pu_height,
            first_pass_if_temp_dst,
            sub_sample_pred_flag,
            mapped_frac_posx ? mapped_frac_posx : mapped_frac_posy);


        // bi-pred Chroma Cb
        picture_average_array[asm_type](
            ref_list0_temp_dst,
            chroma_pu_width << shift,
            ref_list1_temp_dst,
            chroma_pu_width << shift,
            bi_dst->bufferCb + dst_chroma_index,
            bi_dst->strideCb << shift,
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

        avc_style_uni_pred_luma_if_function_ptr_array[asm_type][mapped_frac_posx + (mapped_frac_posy << 3)](
            ref_pic_list0->bufferCr + integ_pos_x + integ_pos_y * ref_pic_list0->strideCr,
            ref_pic_list0->strideCr,
            ref_list0_temp_dst,
            chroma_pu_width,
            chroma_pu_width,
            chroma_pu_height,
            first_pass_if_temp_dst,
            sub_sample_pred_flag,
            mapped_frac_posx ? mapped_frac_posx : mapped_frac_posy);


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

        avc_style_uni_pred_luma_if_function_ptr_array[asm_type][mapped_frac_posx + (mapped_frac_posy << 3)](
            ref_pic_list1->bufferCr + integ_pos_x + integ_pos_y * ref_pic_list1->strideCr,
            ref_pic_list1->strideCr,
            ref_list1_temp_dst,
            chroma_pu_width,
            chroma_pu_width,
            chroma_pu_height,
            first_pass_if_temp_dst,
            sub_sample_pred_flag,
            mapped_frac_posx ? mapped_frac_posx : mapped_frac_posy);

        // bi-pred Chroma Cr
        picture_average_array[asm_type](
            ref_list0_temp_dst,
            chroma_pu_width << shift,
            ref_list1_temp_dst,
            chroma_pu_width << shift,
            bi_dst->bufferCr + dst_chroma_index,
            bi_dst->strideCr << shift,
            chroma_pu_width,
            chroma_pu_height >> shift
            );

    }
}


void estimate_uni_pred_interpolation_avc_lumaRef10Bit(
    EbPictureBufferDesc_t   *ref_frame_pic_list0,
    uint32_t                 pos_x,
    uint32_t                 pos_y,
    uint32_t                 pu_width,
    uint32_t                 pu_height,
    EbPictureBufferDesc_t   *dst,
    uint32_t                 dst_luma_index,
    uint32_t                 dst_chroma_index,          //input parameter, please refer to the detailed explanation above.
    uint32_t                 component_mask,
    EbByte                  temp_buf,
    EbBool                   sub_pred,
    EbBool                   sub_pred_chroma,
    EbAsm                    asm_type)
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
            sub_pred,
            asm_type
        );
    }
    //chroma
    if (component_mask & PICTURE_BUFFER_DESC_CHROMA_MASK) {
        {
            sub_pred = sub_pred_chroma;
            uint32_t in_pos_x = (pos_x >> 3);
            uint32_t in_pos_y = (pos_y >> 3);
            uint16_t *ptr16 = (uint16_t *)ref_frame_pic_list0->bufferCb + in_pos_x + in_pos_y * ref_frame_pic_list0->strideCb;

            extract_8bit_data(
                ptr16,
                ref_frame_pic_list0->strideCb << sub_pred,
                dst->bufferCb + dst_chroma_index,
                dst->strideCb << sub_pred,
                chroma_pu_width,
                chroma_pu_height >> sub_pred,
                asm_type
            );

            ptr16 = (uint16_t *)ref_frame_pic_list0->bufferCr + in_pos_x + in_pos_y * ref_frame_pic_list0->strideCr;

            extract_8bit_data(
                ptr16,
                ref_frame_pic_list0->strideCr << sub_pred,
                dst->bufferCr + dst_chroma_index,
                dst->strideCr << sub_pred,
                chroma_pu_width,
                chroma_pu_height >> sub_pred,
                asm_type
            );
        }
    }
}

void estimate_uni_pred_interpolation_avc_chroma_ref10_bit(
    EbPictureBufferDesc_t   *ref_frame_pic_list0,
    uint32_t                 pos_x,
    uint32_t                 pos_y,
    uint32_t                 pu_width,
    uint32_t                 pu_height,
    EbPictureBufferDesc_t   *dst,
    uint32_t                 dst_luma_index,
    uint32_t                 dst_chroma_index,          //input parameter, please refer to the detailed explanation above.
    uint32_t                 component_mask,
    EbByte                  temp_buf,
    EbBool                   sub_pred,
    EbAsm                    asm_type)
{
    uint32_t   chroma_pu_width = pu_width >> 1;
    uint32_t   chroma_pu_height = pu_height >> 1;
    uint32_t   in_pos_x = (pos_x >> 3);
    uint32_t   in_pos_y = (pos_y >> 3);
    uint16_t  *ptr16 = (uint16_t *)ref_frame_pic_list0->bufferCb + in_pos_x + in_pos_y * ref_frame_pic_list0->strideCb;

    (void)temp_buf;
    (void)component_mask;
    (void)dst_luma_index;

    extract_8bit_data(
        ptr16,
        ref_frame_pic_list0->strideCb << sub_pred,
        dst->bufferCb + dst_chroma_index,
        dst->strideCb << sub_pred,
        chroma_pu_width,
        chroma_pu_height >> sub_pred,
        asm_type
    );

    ptr16 = (uint16_t *)ref_frame_pic_list0->bufferCr + in_pos_x + in_pos_y * ref_frame_pic_list0->strideCr;

    extract_8bit_data(
        ptr16,
        ref_frame_pic_list0->strideCr << sub_pred,
        dst->bufferCr + dst_chroma_index,
        dst->strideCr << sub_pred,
        chroma_pu_width,
        chroma_pu_height >> sub_pred,
        asm_type
    );
}
void estimate_bi_pred_interpolation_avc_chroma_ref10_bit(
    EbPictureBufferDesc_t   *ref_frame_pic_list0,
    EbPictureBufferDesc_t   *ref_frame_pic_list1,
    uint32_t                 ref_list0_pos_x,
    uint32_t                 ref_list0_pos_y,
    uint32_t                 ref_list1_pos_x,
    uint32_t                 ref_list1_pos_y,
    uint32_t                 pu_width,
    uint32_t                 pu_height,
    EbPictureBufferDesc_t   *bi_dst,
    uint32_t                 dst_luma_index,
    uint32_t                 dst_chroma_index,
    uint32_t                 component_mask,
    EbByte                  ref_list0_temp_dst,
    EbByte                  ref_list1_temp_dst,
    EbByte                  first_pass_if_temp_dst,
    EbBool                   sub_pred,
    EbAsm                    asm_type)
{
    uint32_t   chroma_pu_width = pu_width >> 1;
    uint32_t   chroma_pu_height = pu_height >> 1;

    (void)first_pass_if_temp_dst;
    (void)ref_list0_temp_dst;
    (void)ref_list1_temp_dst;
    (void)component_mask;
    (void)dst_luma_index;
    unpack_l0l1_avg(
        (uint16_t *)ref_frame_pic_list0->bufferCb + (ref_list0_pos_x >> 3) + (ref_list0_pos_y >> 3)*ref_frame_pic_list0->strideCb,
        ref_frame_pic_list0->strideCb << sub_pred,
        (uint16_t *)ref_frame_pic_list1->bufferCb + (ref_list1_pos_x >> 3) + (ref_list1_pos_y >> 3)*ref_frame_pic_list1->strideCb,
        ref_frame_pic_list1->strideCb << sub_pred,
        bi_dst->bufferCb + dst_chroma_index,
        bi_dst->strideCb << sub_pred,
        chroma_pu_width,
        chroma_pu_height >> sub_pred,
        asm_type
    );

    unpack_l0l1_avg(
        (uint16_t *)ref_frame_pic_list0->bufferCr + (ref_list0_pos_x >> 3) + (ref_list0_pos_y >> 3)*ref_frame_pic_list0->strideCr,
        ref_frame_pic_list0->strideCr << sub_pred,
        (uint16_t *)ref_frame_pic_list1->bufferCr + (ref_list1_pos_x >> 3) + (ref_list1_pos_y >> 3)*ref_frame_pic_list1->strideCr,
        ref_frame_pic_list1->strideCr << sub_pred,
        bi_dst->bufferCr + dst_chroma_index,
        bi_dst->strideCr << sub_pred,
        chroma_pu_width,
        chroma_pu_height >> sub_pred,
        asm_type
    );
}

void estimate_bi_pred_interpolation_avc_luma_ref10_bit(
    EbPictureBufferDesc_t   *ref_frame_pic_list0,
    EbPictureBufferDesc_t   *ref_frame_pic_list1,
    uint32_t                 ref_list0_pos_x,
    uint32_t                 ref_list0_pos_y,
    uint32_t                 ref_list1_pos_x,
    uint32_t                 ref_list1_pos_y,
    uint32_t                 pu_width,
    uint32_t                 pu_height,
    EbPictureBufferDesc_t   *bi_dst,
    uint32_t                 dst_luma_index,
    uint32_t                 dst_chroma_index,
    uint32_t                 component_mask,
    EbByte                  ref_list0_temp_dst,
    EbByte                  ref_list1_temp_dst,
    EbByte                  first_pass_if_temp_dst,
    EbBool                   sub_pred,
    EbBool                   sub_pred_chroma,
    EbAsm                    asm_type)
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
            sub_pred,
            asm_type
        );
    }

    //uni-prediction List0 chroma
    if (component_mask & PICTURE_BUFFER_DESC_CHROMA_MASK) {
        sub_pred = sub_pred_chroma;
        unpack_l0l1_avg(
            (uint16_t *)ref_frame_pic_list0->bufferCb + (ref_list0_pos_x >> 3) + (ref_list0_pos_y >> 3)*ref_frame_pic_list0->strideCb,
            ref_frame_pic_list0->strideCb << sub_pred,
            (uint16_t *)ref_frame_pic_list1->bufferCb + (ref_list1_pos_x >> 3) + (ref_list1_pos_y >> 3)*ref_frame_pic_list1->strideCb,
            ref_frame_pic_list1->strideCb << sub_pred,
            bi_dst->bufferCb + dst_chroma_index,
            bi_dst->strideCb << sub_pred,
            chroma_pu_width,
            chroma_pu_height >> sub_pred,
            asm_type);

        unpack_l0l1_avg(
            (uint16_t *)ref_frame_pic_list0->bufferCr + (ref_list0_pos_x >> 3) + (ref_list0_pos_y >> 3)*ref_frame_pic_list0->strideCr,
            ref_frame_pic_list0->strideCr << sub_pred,
            (uint16_t *)ref_frame_pic_list1->bufferCr + (ref_list1_pos_x >> 3) + (ref_list1_pos_y >> 3)*ref_frame_pic_list1->strideCr,
            ref_frame_pic_list1->strideCr << sub_pred,
            bi_dst->bufferCr + dst_chroma_index,
            bi_dst->strideCr << sub_pred,
            chroma_pu_width,
            chroma_pu_height >> sub_pred,
            asm_type);
    }
}

void uni_pred_i_free_ref8_bit(
    EbPictureBufferDesc_t   *ref_pic,
    uint32_t                 pos_x,
    uint32_t                 pos_y,
    uint32_t                 pu_width,
    uint32_t                 pu_height,
    EbPictureBufferDesc_t   *dst,
    uint32_t                 dst_luma_index,
    uint32_t                 dst_chroma_index,          //input parameter, please refer to the detailed explanation above.
    uint32_t                 component_mask,
    EbByte                  temp_buf,
    EbBool                   sub_sample_pred_flag,
    EbBool                   sub_sample_pred_flag_chroma,
    EbAsm                    asm_type)
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
        avc_style_uni_pred_luma_if_function_ptr_array[asm_type][mapped_frac_posx + (mapped_frac_posy << 2)](
            ref_pic->buffer_y + integ_pos_x + integ_pos_y * ref_pic->stride_y, ref_pic->stride_y,
            dst->buffer_y + dst_luma_index, luma_stride,
            pu_width, pu_height,
            temp_buf,
            sub_sample_pred_flag,
            mapped_frac_posx ? mapped_frac_posx : mapped_frac_posy);
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
        avc_style_uni_pred_luma_if_function_ptr_array[asm_type][mapped_frac_posx + (mapped_frac_posy << 3)](
            ref_pic->bufferCb + integ_pos_x + integ_pos_y * ref_pic->strideCb,
            ref_pic->strideCb,
            dst->bufferCb + dst_chroma_index,
            dst->strideCb,
            chroma_pu_width,
            chroma_pu_height,
            temp_buf,
            sub_sample_pred_flag_chroma,
            mapped_frac_posx ? mapped_frac_posx : mapped_frac_posy);

        //doing the chroma Cr interpolation
        avc_style_uni_pred_luma_if_function_ptr_array[asm_type][mapped_frac_posx + (mapped_frac_posy << 3)](
            ref_pic->bufferCr + integ_pos_x + integ_pos_y * ref_pic->strideCr,
            ref_pic->strideCr,
            dst->bufferCr + dst_chroma_index,
            dst->strideCr,
            chroma_pu_width,
            chroma_pu_height,
            temp_buf,
            sub_sample_pred_flag_chroma,
            mapped_frac_posx ? mapped_frac_posx : mapped_frac_posy);
    }
}

void bi_pred_i_free_ref8_bit(
    EbPictureBufferDesc_t   *ref_pic_list0,
    EbPictureBufferDesc_t   *ref_pic_list1,
    uint32_t                 ref_list0_pos_x,
    uint32_t                 ref_list0_pos_y,
    uint32_t                 ref_list1_pos_x,
    uint32_t                 ref_list1_pos_y,
    uint32_t                 pu_width,
    uint32_t                 pu_height,
    EbPictureBufferDesc_t   *bi_dst,
    uint32_t                 dst_luma_index,
    uint32_t                 dst_chroma_index,
    uint32_t                 component_mask,
    EbByte                  ref_list0_temp_dst,
    EbByte                  ref_list1_temp_dst,
    EbByte                  first_pass_if_temp_dst,
    EbBool                   sub_sample_pred_flag,
    EbBool                   sub_sample_pred_flag_chroma,//needs to be connected
    EbAsm                    asm_type)
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

        avc_style_uni_pred_luma_if_function_ptr_array[asm_type][mapped_frac_posx + (mapped_frac_posy << 2)](
            ref_pic_list0->buffer_y + integ_pos_x + integ_pos_y * ref_luma_stride, ref_luma_stride,
            ref_list0_temp_dst, pu_width,
            pu_width, pu_height,
            first_pass_if_temp_dst,
            sub_sample_pred_flag,
            mapped_frac_posx ? mapped_frac_posx : mapped_frac_posy);

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
        avc_style_uni_pred_luma_if_function_ptr_array[asm_type][mapped_frac_posx + (mapped_frac_posy << 2)](
            ref_pic_list1->buffer_y + integ_pos_x + integ_pos_y * ref_luma_stride, ref_luma_stride,
            ref_list1_temp_dst, pu_width,
            pu_width, pu_height,
            first_pass_if_temp_dst,
            sub_sample_pred_flag,
            mapped_frac_posx ? mapped_frac_posx : mapped_frac_posy);

        // bi-pred luma
        picture_average_array[asm_type](ref_list0_temp_dst, pu_width << sub_sample_pred_flag, ref_list1_temp_dst, pu_width << sub_sample_pred_flag, bi_dst->buffer_y + dst_luma_index, luma_stride << sub_sample_pred_flag, pu_width, pu_height >> sub_sample_pred_flag);
        if (sub_sample_pred_flag) {
            picture_average1_line_array[asm_type](ref_list0_temp_dst + (pu_height - 1)*pu_width, ref_list1_temp_dst + (pu_height - 1)*pu_width, bi_dst->buffer_y + dst_luma_index + (pu_height - 1)*luma_stride, pu_width);
        }
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

        avc_style_uni_pred_luma_if_function_ptr_array[asm_type][mapped_frac_posx + (mapped_frac_posy << 3)](
            ref_pic_list0->bufferCb + integ_pos_x + integ_pos_y * ref_pic_list0->strideCb,
            ref_pic_list0->strideCb,
            ref_list0_temp_dst,
            chroma_pu_width,
            chroma_pu_width,
            chroma_pu_height,
            first_pass_if_temp_dst,
            sub_sample_pred_flag_chroma,
            mapped_frac_posx ? mapped_frac_posx : mapped_frac_posy);

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

        avc_style_uni_pred_luma_if_function_ptr_array[asm_type][mapped_frac_posx + (mapped_frac_posy << 3)](
            ref_pic_list1->bufferCb + integ_pos_x + integ_pos_y * ref_pic_list1->strideCb,
            ref_pic_list1->strideCb,
            ref_list1_temp_dst,
            chroma_pu_width,
            chroma_pu_width,
            chroma_pu_height,
            first_pass_if_temp_dst,
            sub_sample_pred_flag_chroma,
            mapped_frac_posx ? mapped_frac_posx : mapped_frac_posy);

        // bi-pred Chroma Cb
        picture_average_array[asm_type](
            ref_list0_temp_dst,
            chroma_pu_width << shift,
            ref_list1_temp_dst,
            chroma_pu_width << shift,
            bi_dst->bufferCb + dst_chroma_index,
            bi_dst->strideCb << shift,
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

        avc_style_uni_pred_luma_if_function_ptr_array[asm_type][mapped_frac_posx + (mapped_frac_posy << 3)](
            ref_pic_list0->bufferCr + integ_pos_x + integ_pos_y * ref_pic_list0->strideCr,
            ref_pic_list0->strideCr,
            ref_list0_temp_dst,
            chroma_pu_width,
            chroma_pu_width,
            chroma_pu_height,
            first_pass_if_temp_dst,
            sub_sample_pred_flag_chroma,
            mapped_frac_posx ? mapped_frac_posx : mapped_frac_posy);

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

        avc_style_uni_pred_luma_if_function_ptr_array[asm_type][mapped_frac_posx + (mapped_frac_posy << 3)](
            ref_pic_list1->bufferCr + integ_pos_x + integ_pos_y * ref_pic_list1->strideCr,
            ref_pic_list1->strideCr,
            ref_list1_temp_dst,
            chroma_pu_width,
            chroma_pu_width,
            chroma_pu_height,
            first_pass_if_temp_dst,
            sub_sample_pred_flag_chroma,
            mapped_frac_posx ? mapped_frac_posx : mapped_frac_posy);

        // bi-pred Chroma Cr
        picture_average_array[asm_type](
            ref_list0_temp_dst,
            chroma_pu_width << shift,
            ref_list1_temp_dst,
            chroma_pu_width << shift,
            bi_dst->bufferCr + dst_chroma_index,
            bi_dst->strideCr << shift,
            chroma_pu_width,
            chroma_pu_height >> shift);
    }
}
