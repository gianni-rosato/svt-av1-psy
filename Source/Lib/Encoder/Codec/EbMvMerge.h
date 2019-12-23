/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/
#ifndef EbMvMerge_h
#define EbMvMerge_h

#include "EbDefinitions.h"
#include "EbCodingUnit.h"
#include "EbPredictionUnit.h"
#include "EbNeighborArrays.h"
#ifdef __cplusplus
extern "C" {
#endif
/** This macro is used to compare if two PUs have the same mvs (same ref_pic_index, same MV_x and same MV_y)
        in a particular reference picture list.
    */
#define CHECK_MV_EQUIVALENT(pu1_pred_dir,                                                   \
                            pu1_ref_idx,                                                    \
                            pu1_mv_x,                                                       \
                            pu1_mv_y,                                                       \
                            pu2_pred_dir,                                                   \
                            pu2_ref_idx,                                                    \
                            pu2_mv_x,                                                       \
                            pu2_mv_y,                                                       \
                            ref_pic_list)                                                   \
    ((!(((pu1_pred_dir) + 1) & (1 << (ref_pic_list))) &&                                    \
      !(((pu2_pred_dir) + 1) & (1 << (ref_pic_list)))) ||                                   \
     ((((pu1_pred_dir) + 1) & (1 << (ref_pic_list))) &&                                     \
      (((pu2_pred_dir) + 1) & (1 << (ref_pic_list))) && ((pu1_ref_idx) == (pu2_ref_idx)) && \
      ((pu1_mv_x) == (pu2_mv_x)) && ((pu1_mv_y) == (pu2_mv_y))))

typedef struct MvMergeCandidate {
    Mv      mv[MAX_NUM_OF_REF_PIC_LIST];
    uint8_t prediction_direction;
} MvMergeCandidate;

#ifdef __cplusplus
}
#endif
#endif // EbMvMerge_h
