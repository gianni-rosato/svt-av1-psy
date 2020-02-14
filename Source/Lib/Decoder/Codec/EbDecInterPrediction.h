/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbDecInterPrediction_h
#define EbDecInterPrediction_h

#include "EbInterPrediction.h"
#include "EbDecProcessFrame.h"

#ifdef __cplusplus
extern "C" {
#endif

void svtav1_predict_inter_block(DecModCtxt *dec_mod_ctx, EbDecHandle *dec_hdl,
                                PartitionInfo *part_info, int32_t mi_row, int32_t mi_col,
                                int32_t num_planes);

void svtav1_predict_inter_block_plane(DecModCtxt *dec_mod_ctx, EbDecHandle *dec_hdl,
                                      PartitionInfo *part_info, int32_t plane,
                                      int32_t build_for_obmc, int32_t mi_x, int32_t mi_y, void *dst,
                                      int32_t dst_stride, int32_t some_use_intra,
                                      int32_t bit_depth);

#ifdef __cplusplus
}
#endif
#endif // EbDecInterPrediction_h
