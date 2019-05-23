/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbIntraPrediction_SSSE3_h
#define EbIntraPrediction_SSSE3_h

#include "EbDefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif

extern void intra_mode_angular_vertical_kernel_ssse3_intrin(
    uint32_t     size,
    uint8_t     *ref_samp_main,
    uint8_t     *prediction_ptr,
    uint32_t     prediction_buffer_stride,
    const EbBool skip,
    int32_t      intra_pred_angle);

extern void intra_mode_angular_horizontal_kernel_ssse3_intrin(
    uint32_t     size,
    uint8_t     *ref_samp_main,
    uint8_t     *prediction_ptr,
    uint32_t     prediction_buffer_stride,
    const EbBool skip,
    int32_t      intra_pred_angle);


#ifdef __cplusplus
}
#endif
#endif // EbIntraPrediction_SSSE3_h