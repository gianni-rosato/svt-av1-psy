/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/
#ifndef EbTransforms_AVX2_h
#define EbTransforms_AVX2_h

#include "EbDefinitions.h"
#ifdef __cplusplus
extern "C" {
#endif

    void mat_mult4x4_avx2_intrin(
        int16_t*        coeff,
        const uint32_t  coeff_stride,
        const uint16_t *masking_matrix,
        const uint32_t  masking_matrix_stride,  //Matrix size
        const uint32_t  compute_size,           //Computation area size
        const int32_t   offset,                 //(PMP_MAX >> 1)
        const int32_t   shift_num,              //PMP_PRECISION
        uint32_t*       nonzerocoeff);

    void mat_mult4x4_out_buff_avx2_intrin(
        int16_t*        coeff,
        const uint32_t  coeff_stride,
        int16_t*        coeff_out,
        const uint32_t  coeff_out_stride,
        const uint16_t *masking_matrix,
        const uint32_t  masking_matrix_stride,
        const uint32_t  compute_size,
        const int32_t   offset,
        const int32_t   shift_num,
        uint32_t*       nonzerocoeff);

    void mat_mult8x8_avx2_intrin(
        int16_t*        coeff,
        const uint32_t  coeff_stride,
        const uint16_t *masking_matrix,
        const uint32_t  masking_matrix_stride,  //Matrix size
        const uint32_t  compute_size,           //Computation area size
        const int32_t   offset,                 //(PMP_MAX >> 1)
        const int32_t   shift_num,              //PMP_PRECISION
        uint32_t*       nonzerocoeff);

    void mat_mult_nxn_avx2_intrin(
        int16_t*        coeff,
        const uint32_t  coeff_stride,
        const uint16_t *masking_matrix,
        const uint32_t  masking_matrix_stride,  //Matrix size
        const uint32_t  compute_size,           //Computation area size
        const int32_t   offset,                 //(PMP_MAX >> 1)
        const int32_t   shift_num,              //PMP_PRECISION
        uint32_t*       nonzerocoeff);

#ifdef __cplusplus
}
#endif
#endif // EbTransforms_AVX2_h
