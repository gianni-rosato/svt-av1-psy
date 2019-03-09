/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbDeblockingFilter_SSE2_h
#define EbDeblockingFilter_SSE2_h

#include "EbDefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif


    void aom_highbd_lpf_horizontal_14_sse2(
        uint16_t      *s, 
        int32_t        pitch, 
        const uint8_t *blimit, 
        const uint8_t *limit, 
        const uint8_t *thresh, 
        int32_t        bd);

    void aom_highbd_lpf_horizontal_14_dual_sse2(
        uint16_t      *s, 
        int32_t        pitch, 
        const uint8_t *blimit, 
        const uint8_t *limit, 
        const uint8_t *thresh, 
        int32_t        bd);

    void aom_highbd_lpf_horizontal_4_sse2(
        uint16_t      *s, 
        int32_t        pitch, 
        const uint8_t *blimit, 
        const uint8_t *limit, 
        const uint8_t *thresh, 
        int32_t        bd);

    void aom_highbd_lpf_horizontal_4_dual_sse2(
        uint16_t      *s, 
        int32_t        pitch, 
        const uint8_t *blimit0, 
        const uint8_t *limit0, 
        const uint8_t *thresh0, 
        const uint8_t *blimit1, 
        const uint8_t *limit1, 
        const uint8_t *thresh1,
        int32_t        bd);

    void aom_highbd_lpf_horizontal_6_sse2(
        uint16_t      *s, 
        int32_t        pitch, 
        const uint8_t *blimit, 
        const uint8_t *limit, 
        const uint8_t *thresh, 
        int32_t        bd);

    void aom_highbd_lpf_horizontal_8_sse2(
        uint16_t      *s, 
        int32_t        pitch, 
        const uint8_t *blimit, 
        const uint8_t *limit, 
        const uint8_t *thresh, 
        int32_t        bd);

    void aom_highbd_lpf_horizontal_8_dual_sse2(
        uint16_t      *s, 
        int32_t        pitch, 
        const uint8_t *blimit0, 
        const uint8_t *limit0, 
        const uint8_t *thresh0,
        const uint8_t *blimit1, 
        const uint8_t *limit1, 
        const uint8_t *thresh1, 
        int32_t        bd);

    void aom_highbd_lpf_vertical_14_sse2(
        uint16_t      *s, 
        int32_t        pitch, 
        const uint8_t *blimit,
        const uint8_t *limit, 
        const uint8_t *thresh, 
        int32_t        bd);

    void aom_highbd_lpf_vertical_14_dual_sse2(
        uint16_t      *s, 
        int32_t        pitch, 
        const uint8_t *blimit, 
        const uint8_t *limit, 
        const uint8_t *thresh, 
        int32_t        bd);

    void aom_highbd_lpf_vertical_4_sse2(
        uint16_t      *s, 
        int32_t        pitch, 
        const uint8_t *blimit, 
        const uint8_t *limit, 
        const uint8_t *thresh, 
        int32_t        bd);

    void aom_highbd_lpf_vertical_4_dual_sse2(
        uint16_t      *s, 
        int32_t        pitch,
        const uint8_t *blimit0, 
        const uint8_t *limit0, 
        const uint8_t *thresh0, 
        const uint8_t *blimit1, 
        const uint8_t *limit1, 
        const uint8_t *thresh1, 
        int32_t        bd);

    void aom_highbd_lpf_vertical_6_sse2(
        uint16_t      *s, 
        int32_t        pitch,
        const uint8_t *blimit,
        const uint8_t *limit,
        const uint8_t *thresh, 
        int32_t        bd);

    void aom_highbd_lpf_vertical_8_sse2(
        uint16_t      *s, 
        int32_t        pitch, 
        const uint8_t *blimit, 
        const uint8_t *limit, 
        const uint8_t *thresh, 
        int32_t        bd);

    void aom_highbd_lpf_vertical_8_dual_sse2(
        uint16_t      *s, 
        int32_t        pitch, 
        const uint8_t *blimit0, 
        const uint8_t *limit0, 
        const uint8_t *thresh0, 
        const uint8_t *blimit1, 
        const uint8_t *limit1, 
        const uint8_t *thresh1, 
        int32_t        bd);

    void aom_lpf_horizontal_14_sse2(
        uint8_t       *s, 
        int32_t        pitch, 
        const uint8_t *blimit, 
        const uint8_t *limit, 
        const uint8_t *thresh);

    void aom_lpf_horizontal_14_dual_sse2(
        uint8_t       *s, 
        int32_t        pitch, 
        const uint8_t *blimit, 
        const uint8_t *limit, 
        const uint8_t *thresh);

    void aom_lpf_horizontal_4_sse2(
        uint8_t       *s, 
        int32_t        pitch, 
        const uint8_t *blimit, 
        const uint8_t *limit, 
        const uint8_t *thresh);

    void aom_lpf_horizontal_6_sse2(
        uint8_t       *s, 
        int32_t        pitch, 
        const uint8_t *blimit, 
        const uint8_t *limit, 
        const uint8_t *thresh);

    void aom_lpf_horizontal_8_sse2(
        uint8_t       *s, 
        int32_t        pitch, 
        const uint8_t *blimit, 
        const uint8_t *limit, 
        const uint8_t *thresh);

    void aom_lpf_vertical_14_sse2(
        uint8_t       *s, 
        int32_t        pitch, 
        const uint8_t *blimit, 
        const uint8_t *limit, 
        const uint8_t *thresh);

    void aom_lpf_vertical_14_dual_sse2(
        uint8_t       *s, 
        int32_t        pitch, 
        const uint8_t *blimit, 
        const uint8_t *limit, 
        const uint8_t *thresh);

    void aom_lpf_vertical_4_sse2(
        uint8_t       *s, 
        int32_t        pitch, 
        const uint8_t *blimit, 
        const uint8_t *limit, 
        const uint8_t *thresh);

    void aom_lpf_vertical_6_sse2(
        uint8_t       *s, 
        int32_t        pitch, 
        const uint8_t *blimit, 
        const uint8_t *limit, 
        const uint8_t *thresh);

    void aom_lpf_vertical_8_sse2(
        uint8_t       *s, 
        int32_t        pitch, 
        const uint8_t *blimit,
        const uint8_t *limit, 
        const uint8_t *thresh);

#ifdef __cplusplus
}
#endif
#endif