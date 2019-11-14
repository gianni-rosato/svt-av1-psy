/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/
#ifndef EbDecInverseQuantize_h
#define EbDecInverseQuantize_h

int16_t eb_av1_dc_quant_Q3(int32_t qindex, int32_t delta, AomBitDepth bit_depth);
int16_t eb_av1_ac_quant_Q3(int32_t qindex, int32_t delta, AomBitDepth bit_depth);
void setup_segmentation_dequant(EbDecHandle *dec_handle_ptr, EbColorConfig *color_config);
void av1_inverse_qm_init(EbDecHandle *dec_handle_ptr);
void update_dequant(EbDecHandle *dec_handle, SBInfo *sb_info);
int32_t inverse_quantize(EbDecHandle * dec_handle, PartitionInfo_t *part, BlockModeInfo *mode,
    int32_t *level, int32_t *qcoeffs, TxType tx_type, TxSize tx_size, int plane);

#endif // EbDecInverseQuantize_h
