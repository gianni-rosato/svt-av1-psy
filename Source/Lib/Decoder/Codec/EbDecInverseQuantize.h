/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/
#ifndef EbDecInverseQuantize_h
#define EbDecInverseQuantize_h

#include "EbAv1Structs.h"
#include "EbDecProcessFrame.h"

int16_t eb_av1_dc_quant_q3(int32_t qindex, int32_t delta, AomBitDepth bit_depth);
int16_t eb_av1_ac_quant_q3(int32_t qindex, int32_t delta, AomBitDepth bit_depth);

void    setup_segmentation_dequant(DecModCtxt *dec_mod_ctxt);
void    av1_inverse_qm_init(DecModCtxt *dec_mod_ctxt, SeqHeader *seq_header);
void    update_dequant(DecModCtxt *dec_mod_ctxt, SBInfo *sb_info);
int32_t inverse_quantize(DecModCtxt *dec_mod_ctxt, PartitionInfo *part, BlockModeInfo *mode,
                         int32_t *level, int32_t *qcoeffs, TxType tx_type, TxSize tx_size,
                         int plane);

#endif // EbDecInverseQuantize_h
