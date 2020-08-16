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

#ifndef EbMotionVectorUnit_h
#define EbMotionVectorUnit_h

#include "EbDefinitions.h"
#ifdef __cplusplus
extern "C" {
#endif
#pragma pack(push, 1)

typedef union Mv {
    struct {
        signed short x;
        signed short y;
    };
    uint32_t mv_union;
} Mv;

#pragma pack(pop)

#pragma pack(push, 1)
typedef struct Mvd {
    signed   mvd_x : 16;
    signed   mvd_y : 16;
    unsigned ref_idx : 1;
    unsigned : 7;
    unsigned pred_idx : 1;
    unsigned : 7;
} Mvd;
#pragma pack(pop)

typedef struct MvUnit {
    Mv      mv[MAX_NUM_OF_REF_PIC_LIST];
    uint8_t pred_direction;
} MvUnit;

#ifdef __cplusplus
}
#endif
#endif // EbMotionVectorUnit_h
