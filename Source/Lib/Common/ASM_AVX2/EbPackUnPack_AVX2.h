/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbDefinitions.h"

void eb_enc_un_pack8_bit_data_avx2_intrin(
    uint16_t *in_16bit_buffer,
    uint32_t  in_stride,
    uint8_t  *out_8bit_buffer,
    uint32_t  out_stride,
    uint32_t  width,
    uint32_t  height);
