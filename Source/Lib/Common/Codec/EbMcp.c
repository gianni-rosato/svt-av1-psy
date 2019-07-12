/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <string.h>

#include "EbMcp.h"
#include "EbDefinitions.h"
#include "EbPictureBufferDesc.h"
#include "EbPictureOperators.h"

#if (InternalBitDepthIncrement == 0)
#define ChromaOffset4 (1 << (Shift4 - 1))
#else
#define ChromaOffset4 Offset4
#endif
#if (InternalBitDepthIncrement == 0)
#define ChromaMinusOffset1 0
#else
#define ChromaMinusOffset1 MinusOffset1
#endif

void encode_uni_pred_interpolation(
    EbPictureBufferDesc *ref_pic,                  //input parameter, please refer to the detailed explanation above.
    uint32_t                 pos_x,                    //input parameter, please refer to the detailed explanation above.
    uint32_t                 pos_y,                    //input parameter, please refer to the detailed explanation above.
    uint32_t                 pu_width,                 //input parameter
    uint32_t                 pu_height,                //input parameter
    EbPictureBufferDesc *dst,                     //output parameter, please refer to the detailed explanation above.
    uint32_t                 dst_luma_index,            //input parameter, please refer to the detailed explanation above.
    uint32_t                 dst_chroma_index,          //input parameter, please refer to the detailed explanation above.
    int16_t                *tempBuf0,                //input parameter, please refer to the detailed explanation above.
    int16_t                *tempBuf1,                //input parameter, please refer to the detailed explanation above.
    EbAsm                 asm_type)
{
    uint32_t   integPosx;
    uint32_t   integPosy;
    uint8_t    frac_pos_x;
    uint8_t    frac_pos_y;
    uint32_t   chromaPuWidth = pu_width >> 1;
    uint32_t   chromaPuHeight = pu_height >> 1;

    (void)tempBuf1;

    //luma
    //compute the luma fractional position
    integPosx = (pos_x >> 2);
    integPosy = (pos_y >> 2);
    frac_pos_x = pos_x & 0x03;
    frac_pos_y = pos_y & 0x03;

    uni_pred_luma_if_function_ptr_array_new[asm_type][frac_pos_x + (frac_pos_y << 2)](
        ref_pic->buffer_y + integPosx + integPosy * ref_pic->stride_y,
        ref_pic->stride_y,
        dst->buffer_y + dst_luma_index,
        dst->stride_y,
        pu_width,
        pu_height,
        tempBuf0);

    //chroma
    //compute the chroma fractional position
    integPosx = (pos_x >> 3);
    integPosy = (pos_y >> 3);
    frac_pos_x = pos_x & 0x07;
    frac_pos_y = pos_y & 0x07;

    uni_pred_chroma_if_function_ptr_array_new[asm_type][frac_pos_x + (frac_pos_y << 3)](
        ref_pic->buffer_cb + integPosx + integPosy * ref_pic->stride_cb,
        ref_pic->stride_cb,
        dst->buffer_cb + dst_chroma_index,
        dst->stride_cb,
        chromaPuWidth,
        chromaPuHeight,
        tempBuf0,
        frac_pos_x,
        frac_pos_y);

    //doing the chroma Cr interpolation
    uni_pred_chroma_if_function_ptr_array_new[asm_type][frac_pos_x + (frac_pos_y << 3)](
        ref_pic->buffer_cr + integPosx + integPosy * ref_pic->stride_cr,
        ref_pic->stride_cr,
        dst->buffer_cr + dst_chroma_index,
        dst->stride_cr,
        chromaPuWidth,
        chromaPuHeight,
        tempBuf0,
        frac_pos_x,
        frac_pos_y);
}

/** generate_padding()
        is used to pad the target picture. The horizontal padding happens first and then the vertical padding.
 */
void generate_padding(
    EbByte  src_pic,                    //output paramter, pointer to the source picture to be padded.
    uint32_t   src_stride,                 //input paramter, the stride of the source picture to be padded.
    uint32_t   original_src_width,          //input paramter, the width of the source picture which excludes the padding.
    uint32_t   original_src_height,         //input paramter, the height of the source picture which excludes the padding.
    uint32_t   padding_width,              //input paramter, the padding width.
    uint32_t   padding_height)             //input paramter, the padding height.
{
    uint32_t   verticalIdx = original_src_height;
    EbByte  tempSrcPic0;
    EbByte  tempSrcPic1;
    EbByte  tempSrcPic2;
    EbByte  tempSrcPic3;

    tempSrcPic0 = src_pic + padding_width + padding_height * src_stride;
    while (verticalIdx)
    {
        // horizontal padding
        EB_MEMSET(tempSrcPic0 - padding_width, *tempSrcPic0, padding_width);
        EB_MEMSET(tempSrcPic0 + original_src_width, *(tempSrcPic0 + original_src_width - 1), padding_width);

        tempSrcPic0 += src_stride;
        --verticalIdx;
    }

    // vertical padding
    verticalIdx = padding_height;
    tempSrcPic0 = src_pic + padding_height * src_stride;
    tempSrcPic1 = src_pic + (padding_height + original_src_height - 1)*src_stride;
    tempSrcPic2 = tempSrcPic0;
    tempSrcPic3 = tempSrcPic1;
    while (verticalIdx)
    {
        // top part data copy
        tempSrcPic2 -= src_stride;
        EB_MEMCPY(tempSrcPic2, tempSrcPic0, sizeof(uint8_t)*src_stride);        // uint8_t to be modified
        // bottom part data copy
        tempSrcPic3 += src_stride;
        EB_MEMCPY(tempSrcPic3, tempSrcPic1, sizeof(uint8_t)*src_stride);        // uint8_t to be modified
        --verticalIdx;
    }

    return;
}
/** generate_padding16_bit()
is used to pad the target picture. The horizontal padding happens first and then the vertical padding.
*/
void generate_padding16_bit(
    EbByte  src_pic,                    //output paramter, pointer to the source picture to be padded.
    uint32_t   src_stride,                 //input paramter, the stride of the source picture to be padded.
    uint32_t   original_src_width,          //input paramter, the width of the source picture which excludes the padding.
    uint32_t   original_src_height,         //input paramter, the height of the source picture which excludes the padding.
    uint32_t   padding_width,              //input paramter, the padding width.
    uint32_t   padding_height)             //input paramter, the padding height.
{
    uint32_t   verticalIdx = original_src_height;
    EbByte  tempSrcPic0;
    EbByte  tempSrcPic1;
    EbByte  tempSrcPic2;
    EbByte  tempSrcPic3;

    tempSrcPic0 = src_pic + padding_width + padding_height * src_stride;
    while (verticalIdx)
    {
        // horizontal padding
        //EB_MEMSET(tempSrcPic0 - padding_width, tempSrcPic0, padding_width);
        memset16bit((uint16_t*)(tempSrcPic0 - padding_width), ((uint16_t*)(tempSrcPic0))[0], padding_width >> 1);
        memset16bit((uint16_t*)(tempSrcPic0 + original_src_width), ((uint16_t*)(tempSrcPic0 + original_src_width - 2/*1*/))[0], padding_width >> 1);

        tempSrcPic0 += src_stride;
        --verticalIdx;
    }

    // vertical padding
    verticalIdx = padding_height;
    tempSrcPic0 = src_pic + padding_height * src_stride;
    tempSrcPic1 = src_pic + (padding_height + original_src_height - 1)*src_stride;
    tempSrcPic2 = tempSrcPic0;
    tempSrcPic3 = tempSrcPic1;
    while (verticalIdx)
    {
        // top part data copy
        tempSrcPic2 -= src_stride;
        EB_MEMCPY(tempSrcPic2, tempSrcPic0, sizeof(uint8_t)*src_stride);        // uint8_t to be modified
        // bottom part data copy
        tempSrcPic3 += src_stride;
        EB_MEMCPY(tempSrcPic3, tempSrcPic1, sizeof(uint8_t)*src_stride);        // uint8_t to be modified
        --verticalIdx;
    }

    return;
}

/** pad_input_picture()
is used to pad the input picture in order to get . The horizontal padding happens first and then the vertical padding.
*/
void pad_input_picture(
    EbByte  src_pic,                //output paramter, pointer to the source picture to be padded.
    uint32_t   src_stride,             //input paramter, the stride of the source picture to be padded.
    uint32_t   original_src_width,      //input paramter, the width of the source picture which excludes the padding.
    uint32_t   original_src_height,     //input paramter, the height of the source picture which excludes the padding.
    uint32_t   pad_right,                //input paramter, the padding right.
    uint32_t   pad_bottom)             //input paramter, the padding bottom.
{
    uint32_t   verticalIdx;
    EbByte  tempSrcPic0;
    EbByte  tempSrcPic1;

    if (pad_right) {
        // Add padding @ the right
        verticalIdx = original_src_height;
        tempSrcPic0 = src_pic;

        while (verticalIdx)
        {
            EB_MEMSET(tempSrcPic0 + original_src_width, *(tempSrcPic0 + original_src_width - 1), pad_right);
            tempSrcPic0 += src_stride;
            --verticalIdx;
        }
    }

    if (pad_bottom) {
        // Add padding @ the bottom
        verticalIdx = pad_bottom;
        tempSrcPic0 = src_pic + (original_src_height - 1) * src_stride;
        tempSrcPic1 = tempSrcPic0;

        while (verticalIdx)
        {
            tempSrcPic1 += src_stride;
            EB_MEMCPY(tempSrcPic1, tempSrcPic0, sizeof(uint8_t)* (original_src_width + pad_right));
            --verticalIdx;
        }
    }

    return;
}
