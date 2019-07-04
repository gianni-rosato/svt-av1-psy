/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbDecUtils_h
#define EbDecUtils_h

#ifdef __cplusplus
extern "C" {
#endif

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
                         void **pp_blk_recon_buf, int32_t *recon_strd,
                         int32_t sub_x, int32_t sub_y);

void pad_pic(EbPictureBufferDesc *recon_picture_buf);

#ifdef __cplusplus
}
#endif
#endif // EbDecUtils_h
