/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbPredictionUnit_h
#define EbPredictionUnit_h

#include "EbDefinitions.h"
#include "EbSyntaxElements.h"
#include "EbDefinitions.h"
#include "EbMotionVectorUnit.h"
#ifdef __cplusplus
extern "C" {
#endif
#pragma pack(push, 1)
typedef struct PredictionUnit {
    Mv mv[MAX_NUM_OF_REF_PIC_LIST]; // 16-bytes
    uint8_t inter_pred_direction_index;

    // Intra Mode
    int32_t  angle_delta[PLANE_TYPES];
    EbBool   is_directional_mode_flag;
    EbBool   is_directional_chroma_mode_flag;
    uint32_t intra_chroma_mode;

    // Inter Mode
    PredictionMode         inter_mode;
    EbBool                 is_compound;
    uint8_t  ref_frame_type;
    MotionMode           motion_mode;
    uint16_t             num_proj_ref;
    uint32_t             overlappable_neighbors[2];

    // Index of the alpha Cb and alpha Cr combination
    int32_t cfl_alpha_idx;
    // Joint sign of alpha Cb and alpha Cr
    int32_t cfl_alpha_signs;
} PredictionUnit;
#pragma pack(pop)

#ifdef __cplusplus
}
#endif
#endif //EbPredictionUnit_h
