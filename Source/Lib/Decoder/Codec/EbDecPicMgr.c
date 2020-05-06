/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

// SUMMARY
//   Contains the Decoder Memory Init functions

/**************************************
 * Includes
 **************************************/
#include <stdlib.h>

#include "EbDefinitions.h"
#include "EbPictureBufferDesc.h"

#include "EbSvtAv1Dec.h"
#include "EbDecHandle.h"
#include "EbDecMemInit.h"
#include "EbDecUtils.h"

#include "EbDecPicMgr.h"

#define NUM_REF_FRAMES 8 // TODO: remove (reuse EbObuParse.h macro)

/**
*******************************************************************************
*
* @brief
*  Picture manager initializer
*
* @par Description:
*  Initialises the Picture manager structure
*
* @param[in] ps_pic_mgr
*  Pointer to the Picture manager structure
*
* @returns
*
* @remarks
*
*******************************************************************************
*/

EbErrorType dec_pic_mgr_init(EbDecHandle *dec_handle_ptr) {
    EbDecPicMgr **pps_pic_mgr = (EbDecPicMgr **)&dec_handle_ptr->pv_pic_mgr;
    uint32_t      mi_cols     = 2 * ((dec_handle_ptr->seq_header.max_frame_width + 7) >> 3);
    uint32_t      mi_rows     = 2 * ((dec_handle_ptr->seq_header.max_frame_height + 7) >> 3);
    int           size        = mi_cols * mi_rows;

    EbErrorType return_error = EB_ErrorNone;
    int32_t     i;

    EB_MALLOC_DEC(void *, *pps_pic_mgr, sizeof(EbDecPicMgr), EB_N_PTR);

    EbDecPicMgr *ps_pic_mgr = *pps_pic_mgr;

    for (i = 0; i < MAX_PIC_BUFS; i++) {
        ps_pic_mgr->as_dec_pic[i].ps_pic_buf = NULL;
        ps_pic_mgr->as_dec_pic[i].is_free    = 1;
        ps_pic_mgr->as_dec_pic[i].size       = 0;
        ps_pic_mgr->as_dec_pic[i].ref_count  = 0;
        ps_pic_mgr->as_dec_pic[i].mvs        = NULL;
        EB_MALLOC_DEC(
            uint8_t *, ps_pic_mgr->as_dec_pic[i].segment_maps, size * sizeof(uint8_t), EB_N_PTR);
        memset(ps_pic_mgr->as_dec_pic[i].segment_maps, 0, size);
    }

    ps_pic_mgr->num_pic_bufs = 0;

    return return_error;
}

static INLINE EbErrorType mvs_8x8_memory_alloc(TemporalMvRef **mvs, FrameHeader *frame_info) {
    const int frame_mvs_stride = ROUND_POWER_OF_TWO(frame_info->mi_cols, 1);
    const int frame_mvs_rows   = ROUND_POWER_OF_TWO(frame_info->mi_rows, 1);
    const int mvs_buff_size    = frame_mvs_stride * frame_mvs_rows;

    EB_MALLOC_DEC(TemporalMvRef *, *mvs, mvs_buff_size * sizeof(TemporalMvRef), EB_N_PTR);

    return EB_ErrorNone;
}

/**
*******************************************************************************
*
* @brief
*  Get current Picture buffer
*
* @par Description:
*  Give the free buffer from pool or dynamically allocate if nothing is free
*
* @param[in] ps_pic_mgr
*  Pointer to the Picture manager structure
*
* @returns
*
* @remarks
*
*******************************************************************************
*/

EbDecPicBuf *dec_pic_mgr_get_cur_pic(EbDecHandle *dec_handle_ptr) {
    EbDecPicMgr *ps_pic_mgr = (EbDecPicMgr *)dec_handle_ptr->pv_pic_mgr;
    SeqHeader   *seq_header = &dec_handle_ptr->seq_header;
    FrameHeader *frame_info = &dec_handle_ptr->frame_header;
    EbColorFormat color_format = seq_header->color_config.mono_chrome
        ? EB_YUV400
        : dec_handle_ptr->dec_config.max_color_format;
    int32_t      i;
    EbDecPicBuf *pic_buf = NULL;
    /* TODO: Add lock and unlock for MT */
    // Find a free buffer.
    for (i = 0; i < MAX_PIC_BUFS; i++) {
        if (ps_pic_mgr->as_dec_pic[i].is_free == 1) break;
    }

    if (i >= MAX_PIC_BUFS) return NULL;

    uint16_t       frame_width  = frame_info->frame_size.frame_width;
    uint16_t       frame_height = frame_info->frame_size.frame_height;
    EbColorConfig *cc           = &seq_header->color_config;
    size_t y_size = (frame_width + 2 * DEC_PAD_VALUE) * (frame_height + 2 * DEC_PAD_VALUE);
    size_t uv_size = cc->mono_chrome ? 0 :
                     (((frame_width + 2 * DEC_PAD_VALUE) >> cc->subsampling_x) *
                      ((frame_height + 2 * DEC_PAD_VALUE) >> cc->subsampling_y));
    size_t frame_size = y_size + uv_size;

    if (ps_pic_mgr->as_dec_pic[i].size < frame_size) {
        /* allocate the buffer. TODO: Should add free and allocate logic */

        EbPictureBufferDescInitData input_pic_buf_desc_init_data;
        // Init Picture Init data
        input_pic_buf_desc_init_data.max_width  = seq_header->max_frame_width;
        input_pic_buf_desc_init_data.max_height = seq_header->max_frame_height;
        input_pic_buf_desc_init_data.bit_depth  = (EbBitDepthEnum)cc->bit_depth;
        assert(IMPLIES(cc->mono_chrome, color_format == EB_YUV400));
        input_pic_buf_desc_init_data.color_format = cc->mono_chrome ? EB_YUV400 : color_format;
        input_pic_buf_desc_init_data.buffer_enable_mask =
            cc->mono_chrome ? PICTURE_BUFFER_DESC_LUMA_MASK : PICTURE_BUFFER_DESC_FULL_MASK;

        input_pic_buf_desc_init_data.left_padding  = DEC_PAD_VALUE;
        input_pic_buf_desc_init_data.right_padding = DEC_PAD_VALUE;
        input_pic_buf_desc_init_data.top_padding = DEC_PAD_VALUE;
        input_pic_buf_desc_init_data.bot_padding = DEC_PAD_VALUE;

        input_pic_buf_desc_init_data.split_mode = EB_FALSE;

        EbErrorType return_error = dec_eb_recon_picture_buffer_desc_ctor(
            (EbPtr *)&(ps_pic_mgr->as_dec_pic[i].ps_pic_buf),
            (EbPtr)&input_pic_buf_desc_init_data,
            dec_handle_ptr->is_16bit_pipeline);

        if (return_error != EB_ErrorNone) return NULL;

        ps_pic_mgr->as_dec_pic[i].size = frame_size;

        /* Memory for storing MV's at 8x8 lvl*/
        EbErrorType ret_err = mvs_8x8_memory_alloc(&ps_pic_mgr->as_dec_pic[i].mvs, frame_info);
        if (ret_err != EB_ErrorNone) return NULL;

        ps_pic_mgr->num_pic_bufs++;
    } else
        assert(ps_pic_mgr->as_dec_pic[i].ps_pic_buf != NULL);

    ps_pic_mgr->as_dec_pic[i].is_free   = 0;
    ps_pic_mgr->as_dec_pic[i].ref_count = 1;

    pic_buf = &ps_pic_mgr->as_dec_pic[i];

    return pic_buf;
}

static INLINE void dec_ref_count_and_rel(EbDecPicBuf *ps_pic_buf) {
    if (ps_pic_buf != NULL) {
        ps_pic_buf->ref_count--;
        assert(ps_pic_buf->ref_count >= 0);

        if (ps_pic_buf->ref_count == 0) ps_pic_buf->is_free = 1;
    }
}

/**
*******************************************************************************
*
* @brief
*  Update Reference Frames
*
* @par Description:
*  Update current frame to references. Remove any free frames.
*  Should be called at the end of frame parse/decode.
*
* @param[in] ps_pic_mgr
*  Pointer to the Picture manager structure
*
* @returns
*
* @remarks
*
*******************************************************************************
*/
void dec_pic_mgr_update_ref_pic(EbDecHandle *dec_handle_ptr, int32_t frame_decoded,
                                int32_t refresh_frame_flags) {
    int32_t ref_index = 0, mask;

    /* TODO: Add lock and unlock for MT */
    if (frame_decoded) {
        for (mask = refresh_frame_flags; mask; mask >>= 1) {
            dec_ref_count_and_rel(dec_handle_ptr->ref_frame_map[ref_index]);
            dec_handle_ptr->ref_frame_map[ref_index] =
                dec_handle_ptr->next_ref_frame_map[ref_index];
            dec_handle_ptr->next_ref_frame_map[ref_index] = NULL;
            ++ref_index;
        }

        for (; ref_index < REF_FRAMES; ++ref_index) {
            dec_ref_count_and_rel(dec_handle_ptr->ref_frame_map[ref_index]);
            dec_handle_ptr->ref_frame_map[ref_index] =
                dec_handle_ptr->next_ref_frame_map[ref_index];
            dec_handle_ptr->next_ref_frame_map[ref_index] = NULL;
        }

        if (dec_handle_ptr->frame_header.show_existing_frame) {
            //TODO: Add output Q logic
            //assert(0);
        } else
            dec_ref_count_and_rel(dec_handle_ptr->cur_pic_buf[0]);
    } else {
        // Nothing was decoded, so just drop this frame buffer
        dec_ref_count_and_rel(dec_handle_ptr->cur_pic_buf[0]);
    }

    /* Invalidate these references until the next frame starts. */
    for (ref_index = 0; ref_index < INTER_REFS_PER_FRAME; ref_index++) {
        dec_handle_ptr->remapped_ref_idx[ref_index] = INVALID_IDX;
    }

    for (int i = 0; i < NUM_REF_FRAMES; i++) {
        if ((dec_handle_ptr->frame_header.refresh_frame_flags >> i) & 1) {
            dec_handle_ptr->frame_header.ref_order_hint[i] =
                dec_handle_ptr->frame_header.order_hint;
            dec_handle_ptr->frame_header.ref_frame_id[i] =
                dec_handle_ptr->frame_header.current_frame_id;
        }
    }
}

// Generate next_ref_frame_map.
void generate_next_ref_frame_map(EbDecHandle *dec_handle_ptr) {
    /* TODO: Add lock and unlock for MT */

    // next_ref_frame_map holds references to frame buffers. After storing a
    // frame buffer index in next_ref_frame_map, we need to increase the
    // frame buffer's ref_count.
    int32_t ref_index = 0;
    for (int32_t mask = dec_handle_ptr->frame_header.refresh_frame_flags; mask; mask >>= 1) {
        if (mask & 1)
            dec_handle_ptr->next_ref_frame_map[ref_index] = dec_handle_ptr->cur_pic_buf[0];
        else
            dec_handle_ptr->next_ref_frame_map[ref_index] =
                dec_handle_ptr->ref_frame_map[ref_index];

        if (dec_handle_ptr->next_ref_frame_map[ref_index] != NULL)
            ++dec_handle_ptr->next_ref_frame_map[ref_index]->ref_count;
        ++ref_index;
    }

    for (; ref_index < REF_FRAMES; ++ref_index) {
        dec_handle_ptr->next_ref_frame_map[ref_index] = dec_handle_ptr->ref_frame_map[ref_index];
        if (dec_handle_ptr->next_ref_frame_map[ref_index] != NULL)
            ++dec_handle_ptr->next_ref_frame_map[ref_index]->ref_count;
    }
}

// These functions take a reference frame label between LAST_FRAME and
// EXTREF_FRAME inclusive.  Note that this is different to the indexing
// previously used by the frame_refs[] array.
static INLINE int32_t get_ref_frame_map_with_idx(EbDecHandle *          dec_handle_ptr,
                                                 const MvReferenceFrame ref_frame) {
    return (ref_frame >= LAST_FRAME && ref_frame <= REF_FRAMES)
               ? dec_handle_ptr->remapped_ref_idx[ref_frame - LAST_FRAME]
               : INVALID_IDX;
}

EbDecPicBuf *get_ref_frame_buf(EbDecHandle *dec_handle_ptr, const MvReferenceFrame ref_frame) {
    const int32_t map_idx = get_ref_frame_map_with_idx(dec_handle_ptr, ref_frame);
    return (map_idx != INVALID_IDX) ? dec_handle_ptr->ref_frame_map[map_idx] : NULL;
}

ScaleFactors *get_ref_scale_factors(EbDecHandle *dec_handle_ptr, const MvReferenceFrame ref_frame) {
    const int map_idx = get_ref_frame_map_with_idx(dec_handle_ptr, ref_frame);
    return (map_idx != INVALID_IDX) ? &dec_handle_ptr->ref_scale_factors[map_idx] : NULL;
}

EbDecPicBuf *get_primary_ref_frame_buf(EbDecHandle *dec_handle_ptr) {
    int primary_ref_frame = dec_handle_ptr->frame_header.primary_ref_frame;
    if (primary_ref_frame == PRIMARY_REF_NONE) return NULL;
    const int map_idx = get_ref_frame_map_with_idx(dec_handle_ptr, primary_ref_frame + 1);
    return (map_idx != INVALID_IDX) ? dec_handle_ptr->ref_frame_map[map_idx] : NULL;
}

/* Compares the sort_idx fields. If they are equal, then compares the map_idx
   fields to break the tie. This ensures a stable sort. */
static int compare_ref_frame_info(const void *arg_a, const void *arg_b) {
    const RefFrameInfo *info_a = (RefFrameInfo *)arg_a;
    const RefFrameInfo *info_b = (RefFrameInfo *)arg_b;

    const int sort_idx_diff = info_a->sort_idx - info_b->sort_idx;
    if (sort_idx_diff != 0) return sort_idx_diff;
    return info_a->map_idx - info_b->map_idx;
}

static void set_ref_frame_info(EbDecHandle *dec_handle_ptr, int frame_idx, RefFrameInfo *ref_info) {
    assert(frame_idx >= 0 && frame_idx < INTER_REFS_PER_FRAME);

    dec_handle_ptr->remapped_ref_idx[frame_idx] = ref_info->map_idx;
}

void svt_set_frame_refs(EbDecHandle *dec_handle_ptr, int32_t lst_map_idx, int32_t gld_map_idx) {
    int32_t lst_frame_sort_idx = -1;
    int32_t gld_frame_sort_idx = -1;

    assert(dec_handle_ptr->seq_header.order_hint_info.enable_order_hint);
    assert(dec_handle_ptr->seq_header.order_hint_info.order_hint_bits >= 0);

    const int32_t cur_order_hint = (int32_t)dec_handle_ptr->frame_header.order_hint;
    const int32_t cur_frame_sort_idx =
        1 << (dec_handle_ptr->seq_header.order_hint_info.order_hint_bits - 1);

    RefFrameInfo ref_frame_info[REF_FRAMES];
    int32_t      ref_flag_list[INTER_REFS_PER_FRAME] = {0, 0, 0, 0, 0, 0, 0};

    for (int32_t i = 0; i < REF_FRAMES; ++i) {
        const int32_t map_idx = i;

        ref_frame_info[i].map_idx  = map_idx;
        ref_frame_info[i].sort_idx = -1;

        EbDecPicBuf *const buf    = dec_handle_ptr->ref_frame_map[map_idx];
        ref_frame_info[i].pic_buf = buf;

        if (buf == NULL) continue;
        // If this assertion fails, there is a reference leak.
        assert(buf->ref_count > 0);
        if (buf->ref_count <= 0) continue;

        const int32_t offset = (int32_t)buf->order_hint;
        ref_frame_info[i].sort_idx =
            (offset == -1)
                ? -1
                : cur_frame_sort_idx +
                      get_relative_dist(
                          &dec_handle_ptr->seq_header.order_hint_info, offset, cur_order_hint);
        assert(ref_frame_info[i].sort_idx >= -1);

        if (map_idx == lst_map_idx) lst_frame_sort_idx = ref_frame_info[i].sort_idx;
        if (map_idx == gld_map_idx) gld_frame_sort_idx = ref_frame_info[i].sort_idx;
    }

    // Confirm both LAST_FRAME and GOLDEN_FRAME are valid forward reference
    // frames.
    if (lst_frame_sort_idx == -1 || lst_frame_sort_idx >= cur_frame_sort_idx)
        assert(0); //"Inter frame requests a look-ahead frame as LAST");
    if (gld_frame_sort_idx == -1 || gld_frame_sort_idx >= cur_frame_sort_idx)
        assert(0); //"Inter frame requests a look-ahead frame as GOLDEN");

    // Sort ref frames based on their frame_offset values.
    qsort(ref_frame_info, REF_FRAMES, sizeof(RefFrameInfo), compare_ref_frame_info);

    // Identify forward and backward reference frames.
    // Forward  reference: offset < order_hint
    // Backward reference: offset >= order_hint
    int32_t fwd_start_idx = 0, fwd_end_idx = REF_FRAMES - 1;

    for (int32_t i = 0; i < REF_FRAMES; i++) {
        if (ref_frame_info[i].sort_idx == -1) {
            fwd_start_idx++;
            continue;
        }

        if (ref_frame_info[i].sort_idx >= cur_frame_sort_idx) {
            fwd_end_idx = i - 1;
            break;
        }
    }

    int32_t bwd_start_idx = fwd_end_idx + 1;
    int32_t bwd_end_idx   = REF_FRAMES - 1;

    // === Backward Reference Frames ===

    // == ALTREF_FRAME ==
    if (bwd_start_idx <= bwd_end_idx) {
        set_ref_frame_info(dec_handle_ptr, ALTREF_FRAME - LAST_FRAME, &ref_frame_info[bwd_end_idx]);
        ref_flag_list[ALTREF_FRAME - LAST_FRAME] = 1;
        bwd_end_idx--;
    }

    // == BWDREF_FRAME ==
    if (bwd_start_idx <= bwd_end_idx) {
        set_ref_frame_info(
            dec_handle_ptr, BWDREF_FRAME - LAST_FRAME, &ref_frame_info[bwd_start_idx]);
        ref_flag_list[BWDREF_FRAME - LAST_FRAME] = 1;
        bwd_start_idx++;
    }

    // == ALTREF2_FRAME ==
    if (bwd_start_idx <= bwd_end_idx) {
        set_ref_frame_info(
            dec_handle_ptr, ALTREF2_FRAME - LAST_FRAME, &ref_frame_info[bwd_start_idx]);
        ref_flag_list[ALTREF2_FRAME - LAST_FRAME] = 1;
    }

    // === Forward Reference Frames ===
    for (int32_t i = fwd_start_idx; i <= fwd_end_idx; ++i) {
        // == LAST_FRAME ==
        if (ref_frame_info[i].map_idx == lst_map_idx) {
            set_ref_frame_info(dec_handle_ptr, LAST_FRAME - LAST_FRAME, &ref_frame_info[i]);
            ref_flag_list[LAST_FRAME - LAST_FRAME] = 1;
        }

        // == GOLDEN_FRAME ==
        if (ref_frame_info[i].map_idx == gld_map_idx) {
            set_ref_frame_info(dec_handle_ptr, GOLDEN_FRAME - LAST_FRAME, &ref_frame_info[i]);
            ref_flag_list[GOLDEN_FRAME - LAST_FRAME] = 1;
        }
    }

    assert(ref_flag_list[LAST_FRAME - LAST_FRAME] == 1 &&
           ref_flag_list[GOLDEN_FRAME - LAST_FRAME] == 1);

    // == LAST2_FRAME ==
    // == LAST3_FRAME ==
    // == BWDREF_FRAME ==
    // == ALTREF2_FRAME ==
    // == ALTREF_FRAME ==

    // Set up the reference frames in the anti-chronological order.
    static const MvReferenceFrame ref_frame_list[INTER_REFS_PER_FRAME - 2] = {
        LAST2_FRAME, LAST3_FRAME, BWDREF_FRAME, ALTREF2_FRAME, ALTREF_FRAME};

    int32_t ref_idx;
    for (ref_idx = 0; ref_idx < (INTER_REFS_PER_FRAME - 2); ref_idx++) {
        const MvReferenceFrame ref_frame = ref_frame_list[ref_idx];

        if (ref_flag_list[ref_frame - LAST_FRAME] == 1) continue;

        while (fwd_start_idx <= fwd_end_idx &&
               (ref_frame_info[fwd_end_idx].map_idx == lst_map_idx ||
                ref_frame_info[fwd_end_idx].map_idx == gld_map_idx)) {
            fwd_end_idx--;
        }
        if (fwd_start_idx > fwd_end_idx) break;

        set_ref_frame_info(dec_handle_ptr, ref_frame - LAST_FRAME, &ref_frame_info[fwd_end_idx]);
        ref_flag_list[ref_frame - LAST_FRAME] = 1;

        fwd_end_idx--;
    }

    // Assign all the remaining frame(s), if any, to the earliest reference frame.
    for (; ref_idx < (INTER_REFS_PER_FRAME - 2); ref_idx++) {
        const MvReferenceFrame ref_frame = ref_frame_list[ref_idx];
        if (ref_flag_list[ref_frame - LAST_FRAME] == 1) continue;
        set_ref_frame_info(dec_handle_ptr, ref_frame - LAST_FRAME, &ref_frame_info[fwd_start_idx]);
        ref_flag_list[ref_frame - LAST_FRAME] = 1;
    }

    for (int32_t i = 0; i < INTER_REFS_PER_FRAME; i++) assert(ref_flag_list[i] == 1);
}

void svt_setup_frame_buf_refs(EbDecHandle *dec_handle_ptr) {
    EbDecPicBuf *cur_pic_buf  = dec_handle_ptr->cur_pic_buf[0];
    FrameHeader *frame_header = &dec_handle_ptr->frame_header;

    cur_pic_buf->order_hint              = frame_header->order_hint;
    cur_pic_buf->frame_type              = frame_header->frame_type;
    cur_pic_buf->frame_width             = frame_header->frame_size.frame_width;
    cur_pic_buf->frame_height            = frame_header->frame_size.frame_height;
    cur_pic_buf->render_width            = frame_header->frame_size.render_width;
    cur_pic_buf->render_height           = frame_header->frame_size.render_height;
    cur_pic_buf->superres_upscaled_width = frame_header->frame_size.superres_upscaled_width;

    MvReferenceFrame ref_frame;
    for (ref_frame = LAST_FRAME; ref_frame <= ALTREF_FRAME; ++ref_frame) {
        const EbDecPicBuf *const buf = get_ref_frame_buf(dec_handle_ptr, ref_frame);
        if (buf != NULL)
            dec_handle_ptr->cur_pic_buf[0]->ref_order_hints[ref_frame - LAST_FRAME] =
                buf->order_hint;
    }
}
