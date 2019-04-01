/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbTransQuantBuffers_h
#define EbTransQuantBuffers_h

#include "EbDefinitions.h"
#include "EbPictureBufferDesc.h"

#ifdef __cplusplus
extern "C" {
#endif
    typedef struct EbTransQuantBuffers
    {
        EbPictureBufferDesc_t         *tu_trans_coeff2_nx2_n_ptr;
        EbPictureBufferDesc_t         *tu_trans_coeff_nxn_ptr;
        EbPictureBufferDesc_t         *tu_trans_coeff_n2x_n2_ptr;
        EbPictureBufferDesc_t         *tu_quant_coeff_nxn_ptr;
        EbPictureBufferDesc_t         *tu_quant_coeff_n2x_n2_ptr;

    } EbTransQuantBuffers;


    extern EbErrorType eb_trans_quant_buffers_ctor(
        EbTransQuantBuffers            *trans_quant_buffers_ptr);


#ifdef __cplusplus
}
#endif
#endif // EbTransQuantBuffers_h