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

#ifndef EbTransformUnit_h
#define EbTransformUnit_h

#include "EbDefinitions.h"
#ifdef __cplusplus
extern "C" {
#endif
#define TRANSFORM_UNIT_MAX_COUNT 16
#pragma pack(push, 1)
typedef struct TransformUnit {
    uint16_t nz_coef_count[3];
    TxType   transform_type[PLANE_TYPES];
} TransformUnit;
#pragma pack(pop)
#ifdef __cplusplus
}
#endif
#endif // EbTransformUnit_h
