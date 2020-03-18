/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
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
