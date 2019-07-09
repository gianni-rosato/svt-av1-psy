/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbDecParseHelper_h
#define EbDecParseHelper_h

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

int neg_deinterleave(const int diff, int ref, int max);
void set_segment_id(EbDecHandle *dec_handle, int mi_offset,
    int x_mis, int y_mis, int segment_id);
int bsize_to_max_depth(BlockSize bsize);
int get_tx_size_context(const PartitionInfo_t *xd, ParseCtxt *parse_ctx);
TxSize depth_to_tx_size(int depth, BlockSize bsize);
TxSize read_selected_tx_size(PartitionInfo_t *xd, SvtReader *r,
    EbDecHandle *dec_handle);
int dec_is_inter_block(const ModeInfo_t *mbmi);
int is_intrabc_block(const ModeInfo_t *mbmi);
int max_block_wide(PartitionInfo_t *part_info, int plane_bsize, int subx);
int max_block_high(PartitionInfo_t *part_info, int plane_bsize, int suby);
TxSize get_sqr_tx_size(int tx_dim);
int txfm_partition_context(TXFM_CONTEXT *above_ctx, TXFM_CONTEXT *left_ctx,
    BlockSize bsize, TxSize tx_size);
int get_segdata(SegmentationParams *seg, int segment_id, SEG_LVL_FEATURES feature_id);
int get_intra_inter_context(PartitionInfo_t *xd);
int use_angle_delta(BlockSize bsize);
PredictionMode read_intra_mode(SvtReader *r, AomCdfProb *cdf);
int dec_is_chroma_reference(int mi_row, int mi_col, BlockSize bsize,
    int subsampling_x, int subsampling_y);
UvPredictionMode read_intra_mode_uv(FRAME_CONTEXT *ec_ctx, SvtReader *r,
    CflAllowedType cfl_allowed, PredictionMode y_mode);
CflAllowedType is_cfl_allowed(PartitionInfo_t *xd, EbColorConfig* color_cfg,
    uint8_t *lossless_array);
int allow_palette(int allow_screen_content_tools, BlockSize sb_type);
int filter_intra_allowed_bsize(EbDecHandle *dec_handle, BlockSize bs);
int filter_intra_allowed(EbDecHandle *dec_handle, const ModeInfo_t *mbmi);
int allow_intrabc(const EbDecHandle *dec_handle);
PredictionMode dec_get_uv_mode(UvPredictionMode mode);
TxType intra_mode_to_tx_type(const ModeInfo_t *mbmi, PlaneType plane_type);
int has_second_ref(const ModeInfo_t *mbmi);
IntMv_dec gm_get_motion_vector(const GlobalMotionParams *gm, int allow_hp,
    BlockSize bsize, int mi_col, int mi_row, int is_integer);
int get_txb_wide(TxSize tx_size);
int get_txb_high(TxSize tx_size);
int get_lower_levels_ctx_eob(int bwl, int height, int scan_idx);
uint8_t *set_levels(uint8_t *const levels_buf, const int width);
int get_padded_idx(const int idx, const int bwl);
int get_txb_bwl(TxSize tx_size);
int get_comp_reference_type_context(const PartitionInfo_t *xd);
AomCdfProb *get_y_mode_cdf(FRAME_CONTEXT *tile_ctx, const ModeInfo_t *above_mi,
    const ModeInfo_t *left_mi);

int is_interintra_allowed_bsize(const BlockSize bsize);
int is_interintra_allowed_mode(const PredictionMode mode);
int is_interintra_allowed_ref(const MvReferenceFrame rf[2]);
int is_interintra_allowed(const ModeInfo_t *mbmi);
MotionMode dec_motion_mode_allowed();
void update_tx_context(ParseCtxt *parse_ctxt, PartitionInfo_t *pi,
    BlockSize bsize, TxSize txSize, int blk_row, int blk_col);
int seg_feature_active(SegmentationParams *seg, int segment_id,
    SEG_LVL_FEATURES feature_id);

int find_warp_samples(EbDecHandle *dec_handle, PartitionInfo_t *pi,
    int mi_row, int mi_col, int *pts, int *pts_inref);
#endif  // EbDecParseHelper_h
