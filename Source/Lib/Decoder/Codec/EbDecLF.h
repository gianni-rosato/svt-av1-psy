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

/*SUMMARY
Contains the Decoder Loop Filtering related functions*/

#ifndef EbDecLF_h
#define EbDecLF_h

extern const int32_t mode_lf_lut[];

typedef void(*svt_lbd_filter_tap_fn_t)(uint8_t *s,
    int32_t pitch, const uint8_t *blimit,
    const uint8_t *limit, const uint8_t *thresh);

typedef void(*svt_hbd_filter_tap_fn_t)(uint16_t *s,
    int32_t pitch, const uint8_t *blimit,
    const uint8_t *limit, const uint8_t *thresh, int32_t bd);

typedef struct LFBlockParamL {
    int8_t              skip;
    /*!< Specifies which segment is associated with the
    current intra block being decoded. */
    int8_t              segment_id;
    MvReferenceFrame    ref_frame_0;
    PredictionMode      mode;
    TxSize              tx_size_l;
    BlockSize           bsize;
} LFBlockParamL;

typedef struct LFBlockParamUV {
    TxSize              tx_size_uv;
} LFBlockParamUV;

typedef struct LFCtxt {
    LFBlockParamL *lf_block_luma;
    LFBlockParamUV *lf_block_uv;
    LoopFilterInfoN lf_info;
    int32_t         delta_lf_stride;
}LFCtxt;

void fill_4x4_param_luma(LFBlockParamL* lf_block_l,
    int32_t tu_x, int32_t tu_y, int32_t stride,
    TxSize tx_size, BlockModeInfo *mode_info);

void fill_4x4_param_uv(LFBlockParamUV* lf_block_uv, int32_t tu_x, int32_t tu_y,
    int32_t stride, TxSize tx_size, int32_t sub_x, int32_t sub_y);

void dec_av1_loop_filter_frame(
    EbDecHandle *dec_handle_ptr,
    EbPictureBufferDesc *recon_picture_buf, LFCtxt *lf_ctxt,
    int32_t plane_start, int32_t plane_end, int32_t is_mt,
    int enable_flag);

void set_lbd_lf_filter_tap_functions(void);
void set_hbd_lf_filter_tap_functions(void);

#endif  // EbDecLF_h
