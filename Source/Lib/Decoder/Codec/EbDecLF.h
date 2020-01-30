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

typedef void (*SvtLbdFilterTapFn)(uint8_t *s, int32_t pitch, const uint8_t *blimit,
                                  const uint8_t *limit, const uint8_t *thresh);

typedef void (*SvtHbdFilterTapFn)(uint16_t *s, int32_t pitch, const uint8_t *blimit,
                                  const uint8_t *limit, const uint8_t *thresh, int32_t bd);

typedef struct LfCtxt {
    TxSize *tx_size_l;
    TxSize *tx_size_uv;
    LoopFilterInfoN lf_info;
    int32_t         delta_lf_stride;
}LfCtxt;

void fill_4x4_lf_param(LfCtxt* lf_ctxt,
                       int32_t tu_x, int32_t tu_y,
                       int32_t stride,
                       TxSize tx_size,
                       int32_t sub_x, int32_t sub_y,
                       int plane);

void dec_av1_loop_filter_frame(EbDecHandle *dec_handle_ptr,
                               EbPictureBufferDesc *recon_picture_buf,
                               LfCtxt *lf_ctxt,
                               int32_t plane_start,
                               int32_t plane_end,
                               int32_t is_mt,
                               int enable_flag);

void set_lbd_lf_filter_tap_functions(void);
void set_hbd_lf_filter_tap_functions(void);

#endif // EbDecLF_h
