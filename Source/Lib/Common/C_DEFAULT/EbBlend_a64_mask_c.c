/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

//#include "EbUtility.h"
#include "EbDefinitions.h"


void eb_aom_highbd_blend_a64_vmask_c(uint16_t *dst, uint32_t dst_stride, const uint16_t *src0,
                                     uint32_t src0_stride, const uint16_t *src1,
                                     uint32_t src1_stride, const uint8_t *mask, int w, int h,
                                     int bd) {
    (void)bd;
    int i, j;

    assert(IMPLIES(src0 == dst, src0_stride == dst_stride));
    assert(IMPLIES(src1 == dst, src1_stride == dst_stride));

    assert(h >= 1);
    assert(w >= 1);
    assert(IS_POWER_OF_TWO(h));
    assert(IS_POWER_OF_TWO(w));

    assert(bd == 8 || bd == 10 || bd == 12);

    for (i = 0; i < h; ++i) {
        const int m = mask[i];
        for (j = 0; j < w; ++j) {
            dst[i * dst_stride + j] =
                    AOM_BLEND_A64(m, src0[i * src0_stride + j], src1[i * src1_stride + j]);
        }
    }
}
