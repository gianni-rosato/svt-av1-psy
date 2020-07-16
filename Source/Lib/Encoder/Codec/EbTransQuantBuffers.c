/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbTransQuantBuffers.h"

static void eb_trans_quant_buffers_dctor(EbPtr p) {
    EbTransQuantBuffers* obj = (EbTransQuantBuffers*)p;
    EB_DELETE(obj->txb_trans_coeff2_nx2_n_ptr);
    EB_DELETE(obj->txb_trans_coeff_nxn_ptr);
    EB_DELETE(obj->txb_trans_coeff_n2x_n2_ptr);
    EB_DELETE(obj->txb_quant_coeff_nxn_ptr);
    EB_DELETE(obj->txb_quant_coeff_n2x_n2_ptr);
}

#if SB64_MEM_OPT
EbErrorType eb_trans_quant_buffers_ctor(EbTransQuantBuffers* trans_quant_buffers_ptr, uint8_t sb_size) {
#else
EbErrorType eb_trans_quant_buffers_ctor(EbTransQuantBuffers* trans_quant_buffers_ptr) {
#endif
    EbPictureBufferDescInitData trans_coeff_init_array;

    trans_quant_buffers_ptr->dctor            = eb_trans_quant_buffers_dctor;
#if SB64_MEM_OPT
    trans_coeff_init_array.max_width          = sb_size;
    trans_coeff_init_array.max_height         = sb_size;
#else
    trans_coeff_init_array.max_width          = SB_STRIDE_Y;
    trans_coeff_init_array.max_height         = SB_STRIDE_Y;
#endif
    trans_coeff_init_array.bit_depth          = EB_16BIT;
    trans_coeff_init_array.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
    trans_coeff_init_array.color_format       = EB_YUV420;
    trans_coeff_init_array.left_padding       = 0;
    trans_coeff_init_array.right_padding      = 0;
    trans_coeff_init_array.top_padding        = 0;
    trans_coeff_init_array.bot_padding        = 0;
    trans_coeff_init_array.split_mode         = EB_FALSE;

    EbPictureBufferDescInitData trans_coeff_32bit_init_array;
#if SB64_MEM_OPT
    trans_coeff_32bit_init_array.max_width          = sb_size;
    trans_coeff_32bit_init_array.max_height         = sb_size;
#else
    trans_coeff_32bit_init_array.max_width          = SB_STRIDE_Y;
    trans_coeff_32bit_init_array.max_height         = SB_STRIDE_Y;
#endif
    trans_coeff_32bit_init_array.bit_depth          = EB_32BIT;
    trans_coeff_32bit_init_array.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
    trans_coeff_32bit_init_array.color_format       = EB_YUV420;
    trans_coeff_32bit_init_array.left_padding       = 0;
    trans_coeff_32bit_init_array.right_padding      = 0;
    trans_coeff_32bit_init_array.top_padding        = 0;
    trans_coeff_32bit_init_array.bot_padding        = 0;
    trans_coeff_32bit_init_array.split_mode         = EB_FALSE;

    EB_NEW(trans_quant_buffers_ptr->txb_trans_coeff2_nx2_n_ptr,
           eb_picture_buffer_desc_ctor,
           (EbPtr)&trans_coeff_32bit_init_array);
    EB_NEW(trans_quant_buffers_ptr->txb_trans_coeff_nxn_ptr,
           eb_picture_buffer_desc_ctor,
           (EbPtr)&trans_coeff_32bit_init_array);
    EB_NEW(trans_quant_buffers_ptr->txb_trans_coeff_n2x_n2_ptr,
           eb_picture_buffer_desc_ctor,
           (EbPtr)&trans_coeff_init_array);
    EB_NEW(trans_quant_buffers_ptr->txb_quant_coeff_nxn_ptr,
           eb_picture_buffer_desc_ctor,
           (EbPtr)&trans_coeff_init_array);
    EB_NEW(trans_quant_buffers_ptr->txb_quant_coeff_n2x_n2_ptr,
           eb_picture_buffer_desc_ctor,
           (EbPtr)&trans_coeff_init_array);
    return EB_ErrorNone;
}
