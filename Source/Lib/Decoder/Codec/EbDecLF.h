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

typedef enum EDGE_DIR { VERT_EDGE = 0, HORZ_EDGE = 1, NUM_EDGE_DIRS } EDGE_DIR;
static const int32_t mode_lf_lut[] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // INTRA_MODES
    1, 1, 0, 1,                             // INTER_MODES (GLOBALMV == 0)
    1, 1, 1, 1, 1, 1, 0, 1  // INTER_COMPOUND_MODES (GLOBAL_GLOBALMV == 0)
};

typedef struct AV1_DEBLOCKING_PARAMETERS {
    /*length of the filter applied to the outer edge*/
    uint32_t filter_length;
    /* deblocking limits*/
    const uint8_t *lim;
    const uint8_t *mblim;
    const uint8_t *hev_thr;
} AV1_DEBLOCKING_PARAMETERS;

static const SEG_LVL_FEATURES seg_lvl_lf_lut[MAX_MB_PLANE][2] = {
  { SEG_LVL_ALT_LF_Y_V, SEG_LVL_ALT_LF_Y_H },
  { SEG_LVL_ALT_LF_U, SEG_LVL_ALT_LF_U },
  { SEG_LVL_ALT_LF_V, SEG_LVL_ALT_LF_V }
};

static const int delta_lf_id_lut[MAX_MB_PLANE][2] = { { 0, 1 },
                                                      { 2, 2 },
                                                      { 3, 3 } };

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
}LFCtxt;

void fill_4x4_param_luma(LFBlockParamL* lf_block_l,
    int32_t tu_x, int32_t tu_y, int32_t stride,
    TxSize tx_size, ModeInfo_t *mode_info);

void fill_4x4_param_uv(LFBlockParamUV* lf_block_uv, int32_t tu_x, int32_t tu_y,
    int32_t stride, TxSize tx_size, int32_t sub_x, int32_t sub_y);

void dec_av1_loop_filter_frame(
    FrameHeader *frm_hdr, SeqHeader *seq_header,
    EbPictureBufferDesc *recon_picture_buf, LFCtxt *lf_ctxt,
    int32_t plane_start, int32_t plane_end);
