/*
* Copyright(c) 2019 Intel Corporation
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
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

extern EbErrorType svt_trans_quant_buffers_ctor(EbTransQuantBuffers *trans_quant_buffers_ptr, uint8_t sb_size);

#ifdef __cplusplus
}
#endif
#endif // EbTransQuantBuffers_h
