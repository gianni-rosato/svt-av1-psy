/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbDecProcessBlock_h
#define EbDecProcessBlock_h

#ifdef __cplusplus
extern "C" {
#endif

PartitionType get_partition(DecModCtxt *dec_mod_ctxt, FrameHeader *frame_header, uint32_t mi_row,
                            uint32_t mi_col, SBInfo *sb_info, BlockSize bsize);

// void decode_block();
void decode_block(DecModCtxt *dec_mod_ctxt, BlockModeInfo *mode_info, int32_t mi_row, int32_t mi_col,
                  BlockSize bsize, TileInfo *tile, SBInfo *sb_info /*, uint32_t *recon*/);

#ifdef __cplusplus
}
#endif
#endif // EbDecProcessBlock_h
