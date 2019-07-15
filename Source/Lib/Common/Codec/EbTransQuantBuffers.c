/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbDefinitions.h"
#include "EbTransQuantBuffers.h"
#include "EbPictureBufferDesc.h"

void eb_trans_quant_buffers_dctor(EbPtr p)
{
    EbTransQuantBuffers* obj = (EbTransQuantBuffers*)p;
    EB_DELETE(obj->tu_trans_coeff2_nx2_n_ptr);
    EB_DELETE(obj->tu_trans_coeff_nxn_ptr);
    EB_DELETE(obj->tu_trans_coeff_n2x_n2_ptr);
    EB_DELETE(obj->tu_quant_coeff_nxn_ptr);
    EB_DELETE(obj->tu_quant_coeff_n2x_n2_ptr);
}

EbErrorType eb_trans_quant_buffers_ctor(
    EbTransQuantBuffers          *trans_quant_buffers_ptr)
{
    EbPictureBufferDescInitData transCoeffInitArray;

    trans_quant_buffers_ptr->dctor = eb_trans_quant_buffers_dctor;
    transCoeffInitArray.max_width = SB_STRIDE_Y;
    transCoeffInitArray.max_height = SB_STRIDE_Y;
    transCoeffInitArray.bit_depth = EB_16BIT;
    transCoeffInitArray.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
    transCoeffInitArray.color_format = EB_YUV420;
    transCoeffInitArray.left_padding = 0;
    transCoeffInitArray.right_padding = 0;
    transCoeffInitArray.top_padding = 0;
    transCoeffInitArray.bot_padding = 0;
    transCoeffInitArray.split_mode = EB_FALSE;

    EbPictureBufferDescInitData ThirtyTwoBittransCoeffInitArray;
    ThirtyTwoBittransCoeffInitArray.max_width = SB_STRIDE_Y;
    ThirtyTwoBittransCoeffInitArray.max_height = SB_STRIDE_Y;
    ThirtyTwoBittransCoeffInitArray.bit_depth = EB_32BIT;
    ThirtyTwoBittransCoeffInitArray.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
    ThirtyTwoBittransCoeffInitArray.color_format = EB_YUV420;
    ThirtyTwoBittransCoeffInitArray.left_padding = 0;
    ThirtyTwoBittransCoeffInitArray.right_padding = 0;
    ThirtyTwoBittransCoeffInitArray.top_padding = 0;
    ThirtyTwoBittransCoeffInitArray.bot_padding = 0;
    ThirtyTwoBittransCoeffInitArray.split_mode = EB_FALSE;

    EB_NEW(
        trans_quant_buffers_ptr->tu_trans_coeff2_nx2_n_ptr,
        eb_picture_buffer_desc_ctor,
        (EbPtr)&ThirtyTwoBittransCoeffInitArray);
    EB_NEW(
        trans_quant_buffers_ptr->tu_trans_coeff_nxn_ptr,
        eb_picture_buffer_desc_ctor,
        (EbPtr)&ThirtyTwoBittransCoeffInitArray);
    EB_NEW(
        trans_quant_buffers_ptr->tu_trans_coeff_n2x_n2_ptr,
        eb_picture_buffer_desc_ctor,
        (EbPtr)&transCoeffInitArray);
    EB_NEW(
        trans_quant_buffers_ptr->tu_quant_coeff_nxn_ptr,
        eb_picture_buffer_desc_ctor,
        (EbPtr)&transCoeffInitArray);
    EB_NEW(
        trans_quant_buffers_ptr->tu_quant_coeff_n2x_n2_ptr,
        eb_picture_buffer_desc_ctor,
        (EbPtr)&transCoeffInitArray);
    return EB_ErrorNone;
}
