/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbDecParseHelper_h
#define EbDecParseHelper_h

#include "EbObuParse.h"
#include "EbDecParseFrame.h"
#include "EbUtility.h"

#define ACCT_STR __func__

#define ZERO_ARRAY(dest, n) memset(dest, 0, n * sizeof(*(dest)))

typedef struct MvCount {
    uint8_t newmv_count;
    uint8_t num_mv_found[MODE_CTX_REF_FRAMES];
    uint8_t found_above_match;
    uint8_t found_left_match;
} MvCount;

static INLINE CflAllowedType is_cfl_allowed(PartitionInfo *xd, EbColorConfig *color_cfg,
                                            uint8_t *lossless_array) {
    const BlockModeInfo *mbmi  = xd->mi;
    const BlockSize      bsize = mbmi->sb_type;
    assert(bsize < BlockSizeS_ALL);
    if (lossless_array[mbmi->segment_id]) {
        // In lossless, CfL is available when the partition size is equal to the
        // transform size.
        const int ssx         = color_cfg->subsampling_x;
        const int ssy         = color_cfg->subsampling_y;
        const int plane_bsize = get_plane_block_size(bsize, ssx, ssy);
        return (CflAllowedType)(plane_bsize == BLOCK_4X4);
    }
    // Spec: CfL is available to luma partitions lesser than or equal to 32x32
    return (CflAllowedType)(block_size_wide[bsize] <= 32 && block_size_high[bsize] <= 32);
}

//extern int is_inter_block(const BlockModeInfo *mbmi);

static INLINE int allow_palette(int allow_screen_content_tools, BlockSize sb_type) {
    return allow_screen_content_tools && block_size_wide[sb_type] <= 64 &&
           block_size_high[sb_type] <= 64 && sb_type >= BLOCK_8X8;
}

static INLINE int max_block_wide(PartitionInfo *part_info, int plane_bsize, int subx) {
    int max_blocks_wide = block_size_wide[plane_bsize];
    if (part_info->mb_to_right_edge < 0)
        max_blocks_wide += part_info->mb_to_right_edge >> (3 + subx);
    //Scale width in the transform block unit.
    return max_blocks_wide >> tx_size_wide_log2[0];
}

static INLINE int max_block_high(PartitionInfo *part_info, int plane_bsize, int suby) {
    int max_blocks_high = block_size_high[plane_bsize];
    if (part_info->mb_to_bottom_edge < 0)
        max_blocks_high += part_info->mb_to_bottom_edge >> (3 + suby);
    // Scale the height in the transform block unit.
    return max_blocks_high >> tx_size_high_log2[0];
}

TxSize           read_selected_tx_size(PartitionInfo *xd, ParseCtxt *parse_ctxt);
PredictionMode   read_intra_mode(SvtReader *r, AomCdfProb *cdf);
UvPredictionMode read_intra_mode_uv(FRAME_CONTEXT *ec_ctx, SvtReader *r, CflAllowedType cfl_allowed,
                                    PredictionMode y_mode);
IntMv gm_get_motion_vector(const GlobalMotionParams *gm, int allow_hp, BlockSize bsize, int mi_col,
                           int mi_row, int is_integer);

void set_segment_id(EbDecHandle *dec_handle, int mi_offset, int x_mis, int y_mis, int segment_id);
void update_tx_context(ParseCtxt *parse_ctxt, PartitionInfo *pi, BlockSize bsize, TxSize tx_size,
                       int blk_row, int blk_col);

int neg_deinterleave(const int diff, int ref, int max);
int get_intra_inter_context(PartitionInfo *xd);
int get_comp_reference_type_context(const PartitionInfo *xd);
int seg_feature_active(SegmentationParams *seg, int segment_id, SEG_LVL_FEATURES feature_id);
int find_warp_samples(EbDecHandle *dec_handle, TileInfo *tile, PartitionInfo *pi, int *pts,
                      int *pts_inref);
#endif // EbDecParseHelper_h
