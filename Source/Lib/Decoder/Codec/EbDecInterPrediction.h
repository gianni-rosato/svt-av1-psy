/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbDecInterPrediction_h
#define EbDecInterPrediction_h

#ifdef __cplusplus
extern "C" {
#endif

void svtav1_predict_inter_block(
    EbDecHandle *dec_hdl, PartitionInfo_t *part_info,
    int32_t mi_row, int32_t mi_col, int32_t num_planes);

#ifdef __cplusplus
    }
#endif
#endif // EbDecInterPrediction_h
