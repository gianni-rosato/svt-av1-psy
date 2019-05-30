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

    EbEncPack2DType compressed_pack_func_ptr_array[ASM_TYPE_TOTAL] =
    {
        // NON_AVX2
        compressed_packmsb,
        // AVX2
        compressed_packmsb_avx2_intrin,
    };

    typedef void(*CompPackType)(
        const uint8_t *inn_bit_buffer,
        uint32_t  inn_stride,
        uint8_t  *in_compn_bit_buffer,
        uint32_t  out_stride,
        uint8_t  *local_cache,
        uint32_t  width,
        uint32_t  height);

    CompPackType  convert_unpack_c_pack_func_ptr_array[ASM_TYPE_TOTAL] =
    {
        // NON_AVX2
        c_pack_c,
        // AVX2
        c_pack_avx2_intrin,
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

    typedef void(*EbEncUnpackAvgType)(
        uint16_t *ref16_l0,
        uint32_t  ref_l0_stride,
        uint16_t *ref16_l1,
        uint32_t  ref_l1_stride,
        uint8_t  *dst_ptr,
        uint32_t  dst_stride,
        uint32_t  width,
        uint32_t  height);

    EbEncUnpackAvgType un_pack_avg_func_ptr_array[ASM_TYPE_TOTAL] =
    {
        // NON_AVX2
        unpack_avg,
        // AVX2
        unpack_avg_avx2_intrin,//unpack_avg_sse2_intrin,
    };

    typedef void(*EbEncUnpackAvgSubType)(
        uint16_t *ref16_l0,
        uint32_t  ref_l0_stride,
        uint16_t *ref16_l1,
        uint32_t  ref_l1_stride,
        uint8_t  *dst_ptr,
        uint32_t  dst_stride,
        EbBool    sub_pred,
        uint32_t  width,
        uint32_t  height);

    EbEncUnpackAvgSubType unpack_avg_safe_sub_func_ptr_array[ASM_TYPE_TOTAL] =
    {
        // NON_AVX2
        unpack_avg_safe_sub,
        // AVX2  SafeSub
        unpack_avg_safe_sub_avx2_intrin,//unpack_avg_sse2_intrin,
    };

    typedef void(*EbEncUnPack8BitDataType)(
        uint16_t *in16_bit_buffer,
        uint32_t  in_stride,
        uint8_t  *out8_bit_buffer,
        uint32_t  out8_stride,
        uint32_t  width,
        uint32_t  height);

    EbEncUnPack8BitDataType unpack8_bit_func_ptr_array_16_bit[2][ASM_TYPE_TOTAL] =
    {
        {
           un_pack8_bit_data,
           un_pack8_bit_data,
        },
        {
            // NON_AVX2
            un_pack8_bit_data,
            // AVX2
            eb_enc_un_pack8_bit_data_avx2_intrin,
        }
    };

    typedef void(*EbEncUnPack8BitDataSubType)(
        uint16_t *in16_bit_buffer,
        uint32_t  in_stride,
        uint8_t  *out8_bit_buffer,
        uint32_t  out8_stride,
        uint32_t  width,
        uint32_t  height);

    EbEncUnPack8BitDataSubType unpack8_bit_safe_sub_func_ptr_array_16_bit[ASM_TYPE_TOTAL] =
    {
        // NON_AVX2
        un_pack8_bit_data,
        // AVX2
        eb_enc_un_pack8_bit_data_avx2_intrin,
    };

#ifdef __cplusplus
}
#endif
#endif // EbPackUnPack_h
