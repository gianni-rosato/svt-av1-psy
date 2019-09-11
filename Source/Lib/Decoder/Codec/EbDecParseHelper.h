/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbDecParseHelper_h
#define EbDecParseHelper_h

#include "EbObuParse.h"

#define ACCT_STR __func__

static const PredictionMode fimode_to_intradir[FILTER_INTRA_MODES] = {
  DC_PRED, V_PRED, H_PRED, D157_PRED, DC_PRED
};

typedef struct MvCount{
    uint8_t newmv_count;
    uint8_t num_mv_found[MODE_CTX_REF_FRAMES];
    uint8_t found_match;
    uint8_t found_above_match;
    uint8_t found_left_match;
}MvCount;

static INLINE int get_segdata(SegmentationParams *seg, int segment_id,
    SEG_LVL_FEATURES feature_id)
{
    return seg->feature_data[segment_id][feature_id];
}

static INLINE int is_interintra_allowed_bsize(const BlockSize bsize) {
    return (bsize >= BLOCK_8X8) && (bsize <= BLOCK_32X32);
}

static INLINE int is_interintra_allowed_mode(const PredictionMode mode) {
    return (mode >= SINGLE_INTER_MODE_START) && (mode < SINGLE_INTER_MODE_END);
}

static INLINE int is_interintra_allowed_ref(const MvReferenceFrame rf[2]) {
    return (rf[0] > INTRA_FRAME) && (rf[1] <= INTRA_FRAME);
}

static INLINE int is_interintra_allowed(const ModeInfo_t *mbmi) {
    return is_interintra_allowed_bsize(mbmi->sb_type) &&
        is_interintra_allowed_mode(mbmi->mode) &&
        is_interintra_allowed_ref(mbmi->ref_frame);
}

static INLINE int has_second_ref(const ModeInfo_t *mbmi) {
    return mbmi->ref_frame[1] > INTRA_FRAME;
}

static INLINE CflAllowedType is_cfl_allowed(PartitionInfo_t *xd,
    EbColorConfig* color_cfg, uint8_t *lossless_array)
{
    const ModeInfo_t *mbmi = xd->mi;
    const BlockSize bsize = mbmi->sb_type;
    assert(bsize < BlockSizeS_ALL);
    if (lossless_array[mbmi->segment_id]) {
        // In lossless, CfL is available when the partition size is equal to the
        // transform size.
        const int ssx = color_cfg->subsampling_x;
        const int ssy = color_cfg->subsampling_y;
        const int plane_bsize = get_plane_block_size(bsize, ssx, ssy);
        return (CflAllowedType)(plane_bsize == BLOCK_4X4);
    }
    // Spec: CfL is available to luma partitions lesser than or equal to 32x32
    return (CflAllowedType)(block_size_wide[bsize] <= 32 &&
        block_size_high[bsize] <= 32);
}

static INLINE int is_intrabc_block(const ModeInfo_t *mbmi) {
    return mbmi->use_intrabc;
}

/* TODO : Harmonize with is_inter_block*/
static INLINE int dec_is_inter_block(const ModeInfo_t *mbmi) {
    return is_intrabc_block(mbmi) || mbmi->ref_frame[0] > INTRA_FRAME;
}

static INLINE int allow_palette(int allow_screen_content_tools, BlockSize sb_type) {
    return allow_screen_content_tools && block_size_wide[sb_type] <= 64 &&
        block_size_high[sb_type] <= 64 && sb_type >= BLOCK_8X8;
}

static INLINE uint8_t *set_levels(uint8_t *const levels_buf, const int width) {
    return levels_buf + TX_PAD_TOP * (width + TX_PAD_HOR);
}

static INLINE int get_txb_bwl(TxSize tx_size) {
    tx_size = av1_get_adjusted_tx_size(tx_size);
    return tx_size_wide_log2[tx_size];
}

static INLINE int get_padded_idx(const int idx, const int bwl) {
    return idx + ((idx >> bwl) << TX_PAD_HOR_LOG2);
}

static INLINE int get_txb_wide(TxSize tx_size) {
    tx_size = av1_get_adjusted_tx_size(tx_size);
    return tx_size_wide[tx_size];
}

static INLINE int get_txb_high(TxSize tx_size) {
    tx_size = av1_get_adjusted_tx_size(tx_size);
    return tx_size_high[tx_size];
}

static INLINE int max_block_wide(PartitionInfo_t *part_info, int plane_bsize, int subx) {
    int max_blocks_wide = block_size_wide[plane_bsize];
    if (part_info->mb_to_right_edge < 0)
        max_blocks_wide += part_info->mb_to_right_edge >> (3 + subx);
    //Scale width in the transform block unit.
    return max_blocks_wide >> tx_size_wide_log2[0];
}

static INLINE int max_block_high(PartitionInfo_t *part_info, int plane_bsize, int suby) {
    int max_blocks_high = block_size_high[plane_bsize];
    if (part_info->mb_to_bottom_edge < 0)
        max_blocks_high += part_info->mb_to_bottom_edge >> (3 + suby);
    // Scale the height in the transform block unit.
    return max_blocks_high >> tx_size_high_log2[0];
}

TxSize read_selected_tx_size(PartitionInfo_t *xd, SvtReader *r,
    EbDecHandle *dec_handle);
PredictionMode read_intra_mode(SvtReader *r, AomCdfProb *cdf);
UvPredictionMode read_intra_mode_uv(FRAME_CONTEXT *ec_ctx, SvtReader *r,
    CflAllowedType cfl_allowed, PredictionMode y_mode);
IntMvDec gm_get_motion_vector(const GlobalMotionParams *gm, int allow_hp,
    BlockSize bsize, int mi_col, int mi_row, int is_integer);

void set_segment_id(EbDecHandle *dec_handle, int mi_offset,
    int x_mis, int y_mis, int segment_id);
void update_tx_context(ParseCtxt *parse_ctxt, PartitionInfo_t *pi,
    BlockSize bsize, TxSize txSize, int blk_row, int blk_col);

int neg_deinterleave(const int diff, int ref, int max);
int get_intra_inter_context(PartitionInfo_t *xd);
int dec_is_chroma_reference(int mi_row, int mi_col, BlockSize bsize,
    int subsampling_x, int subsampling_y);
int get_comp_reference_type_context(const PartitionInfo_t *xd);
int seg_feature_active(SegmentationParams *seg, int segment_id,
    SEG_LVL_FEATURES feature_id);
int find_warp_samples(EbDecHandle *dec_handle, PartitionInfo_t *pi,
    int mi_row, int mi_col, int *pts, int *pts_inref);
#endif  // EbDecParseHelper_h
