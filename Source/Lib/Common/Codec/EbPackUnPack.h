/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbPackUnPack_h
#define EbPackUnPack_h
#ifdef __cplusplus
extern "C" {
#endif

#include "EbPackUnPack_C.h"
#include "EbPackUnPack_SSE2.h"
#include "EbPackUnPack_AVX2.h"
#include "EbPictureOperators.h"

    typedef void(*EbEncPack2DType)(
        uint8_t     *in8_bit_buffer,
        uint32_t     in8_stride,
        uint8_t     *inn_bit_buffer,
        uint16_t    *out16_bit_buffer,
        uint32_t     inn_stride,
        uint32_t     out_stride,
        uint32_t     width,
        uint32_t     height);

    EbEncPack2DType pack2d_func_ptr_array_16_bit_src[2][ASM_TYPE_TOTAL] =
    {
        {
            // NON_AVX2
            eb_enc_msb_pack2_d,
            // AVX2
            eb_enc_msb_pack2_d,
        },
        {
            // NON_AVX2
            eb_enc_msb_pack2d_sse2_intrin,
            // AVX2
            eb_enc_msb_pack2d_avx2_intrin_al,//EB_ENC_msbPack2D_AVX2
        }
    };

    typedef void(*EbEncUnPack2DType)(
        uint16_t *in16_bit_buffer,
        uint32_t  in_stride,
        uint8_t  *out8_bit_buffer,
        uint8_t  *outn_bit_buffer,
        uint32_t  out8_stride,
        uint32_t  outn_stride,
        uint32_t  width,
        uint32_t  height);

    EbEncUnPack2DType un_pack2d_func_ptr_array_16_bit[2][ASM_TYPE_TOTAL] =
    {
        {
            // NON_AVX2
            eb_enc_msb_un_pack2_d,
            // AVX2
            eb_enc_msb_un_pack2_d,
        },
        {
            // NON_AVX2
            eb_enc_msb_un_pack2d_sse2_intrin,
            // AVX2
            eb_enc_msb_un_pack2d_sse2_intrin,
        }
    };

#ifdef __cplusplus
}
#endif
#endif // EbPackUnPack_h
