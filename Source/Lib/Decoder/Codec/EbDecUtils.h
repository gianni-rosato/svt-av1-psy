/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbDecUtils_h
#define EbDecUtils_h

#ifdef __cplusplus
extern "C" {
#endif

#define WIENER_COEFFS 3

static INLINE int get_relative_dist(OrderHintInfo *ps_order_hint_info,
                                    int ref_hint, int order_hint)
{
    int diff, m;
    if (!ps_order_hint_info->enable_order_hint)
        return 0;
    diff = ref_hint - order_hint;
    m = 1 << (ps_order_hint_info->order_hint_bits - 1);
    diff = (diff & (m - 1)) - (diff & m);
    return diff;
}

EbErrorType check_add_tplmv_buf(EbDecHandle *dec_handle_ptr);

void derive_blk_pointers(EbPictureBufferDesc *recon_picture_buf, int32_t plane,
                         int32_t blk_col_px, int32_t blk_row_px,
                         void **pp_blk_recon_buf, int32_t *recon_stride,
                         int32_t sub_x, int32_t sub_y);

void pad_pic(EbPictureBufferDesc *recon_picture_buf, FrameHeader *frame_hdr);

int inverse_recenter(int r, int v);

static INLINE int is_interintra_allowed_bsize(const BlockSize bsize) {
    return (bsize >= BLOCK_8X8) && (bsize <= BLOCK_32X32);
}

static INLINE int is_interintra_allowed_mode(const PredictionMode mode) {
    return (mode >= SINGLE_INTER_MODE_START) && (mode < SINGLE_INTER_MODE_END);
}

static INLINE int is_interintra_allowed_ref(const MvReferenceFrame rf[2]) {
    return (rf[0] > INTRA_FRAME) && (rf[1] <= INTRA_FRAME);
}

static INLINE int is_interintra_allowed(const BlockModeInfo *mbmi) {
    return is_interintra_allowed_bsize(mbmi->sb_type) &&
        is_interintra_allowed_mode(mbmi->mode) &&
        is_interintra_allowed_ref(mbmi->ref_frame);
}

#if COMP_INTERINTRA
static INLINE int is_interintra_pred(const BlockModeInfo *mbmi) {
    return mbmi->ref_frame[0] > INTRA_FRAME &&
        mbmi->ref_frame[1] == INTRA_FRAME && is_interintra_allowed(mbmi);
}
#endif //comp_interintra

static INLINE int has_second_ref(const BlockModeInfo *mbmi) {
    return mbmi->ref_frame[1] > INTRA_FRAME;
}

#ifdef __cplusplus
}
#endif
#endif // EbDecUtils_h
