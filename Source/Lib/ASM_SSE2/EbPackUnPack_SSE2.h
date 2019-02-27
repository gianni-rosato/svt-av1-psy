/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbPackUnPack_asm_h
#define EbPackUnPack_asm_h

#include "EbDefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif

    void eb_enc_msb_pack2d_sse2_intrin(
        uint8_t  *in8_bit_buffer,
        uint32_t  in8_stride,
        uint8_t  *inn_bit_buffer,
        uint16_t *out16_bit_buffer,
        uint32_t  inn_stride,
        uint32_t  out_stride,
        uint32_t  width,
        uint32_t  height);

    void eb_enc_msb_un_pack2d_sse2_intrin(
        uint16_t *in16_bit_buffer,
        uint32_t  in_stride,
        uint8_t  *out8_bit_buffer,
        uint8_t  *outn_bit_buffer,
        uint32_t  out8_stride,
        uint32_t  outn_stride,
        uint32_t  width,
        uint32_t  height);

    void eb_enc_un_pack8_bit_data_sse2_intrin(
        uint16_t *in16_bit_buffer,
        uint32_t  in_stride,
        uint8_t  *out8_bit_buffer,
        uint32_t  out8_stride,
        uint32_t  width,
        uint32_t  height);

    void eb_enc_un_pack8_bit_data_safe_sub_sse2_intrin(
        uint16_t *in16_bit_buffer,
        uint32_t  in_stride,
        uint8_t  *out8_bit_buffer,
        uint32_t  out8_stride,
        uint32_t  width,
        uint32_t  height,
        EbBool    sub_pred);

    void unpack_avg_sse2_intrin(
        uint16_t *ref16_l0,
        uint32_t  ref_l0_stride,
        uint16_t *ref16_l1,
        uint32_t  ref_l1_stride,
        uint8_t  *dst_ptr,
        uint32_t  dst_stride,
        uint32_t  width,
        uint32_t  height);

#ifdef __cplusplus
}
#endif
#endif // EbPackUnPack_asm_h

