/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbDecIntraPrediction_h
#define EbDecIntraPrediction_h

#ifdef __cplusplus
extern "C" {
#endif

// Do we need to save the luma pixels from the current block,
// for a possible future CfL prediction?

    CflAllowedType store_cfl_required(const EbColorConfig *cc,
                                      const PartitionInfo_t  *xd);

void svt_av1_predict_intra(DecModCtxt *dec_mod_ctxt, PartitionInfo_t *part_info,
        int32_t plane,
        TxSize tx_size, TileInfo *td,
        uint8_t *blk_recon_buf, int32_t recon_strd,
        EbBitDepthEnum bit_depth, int32_t blk_mi_col_off, int32_t blk_mi_row_off);

void cfl_store_tx(PartitionInfo_t *xd, CflCtx *cfl_ctx, int row, int col, TxSize tx_size,
    BlockSize  bsize, EbColorConfig *cc, void *pv_dst_buff,
    uint32_t dst_stride);

#ifdef __cplusplus
    }
#endif
#endif // EbDecIntraPrediction_h
