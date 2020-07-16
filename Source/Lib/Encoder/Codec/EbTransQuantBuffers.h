/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbTransQuantBuffers_h
#define EbTransQuantBuffers_h

#include "EbPictureBufferDesc.h"
#include "EbObject.h"

#ifdef __cplusplus
extern "C" {
#endif
typedef struct EbTransQuantBuffers {
    EbDctor              dctor;
    EbPictureBufferDesc *txb_trans_coeff2_nx2_n_ptr;
    EbPictureBufferDesc *txb_trans_coeff_nxn_ptr;
    EbPictureBufferDesc *txb_trans_coeff_n2x_n2_ptr;
    EbPictureBufferDesc *txb_quant_coeff_nxn_ptr;
    EbPictureBufferDesc *txb_quant_coeff_n2x_n2_ptr;
} EbTransQuantBuffers;

#if SB64_MEM_OPT
extern EbErrorType eb_trans_quant_buffers_ctor(EbTransQuantBuffers *trans_quant_buffers_ptr, uint8_t sb_size);
#else
extern EbErrorType eb_trans_quant_buffers_ctor(EbTransQuantBuffers *trans_quant_buffers_ptr);
#endif

#ifdef __cplusplus
}
#endif
#endif // EbTransQuantBuffers_h
