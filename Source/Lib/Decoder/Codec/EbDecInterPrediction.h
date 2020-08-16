/*
* Copyright(c) 2019 Netflix, Inc.
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
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
