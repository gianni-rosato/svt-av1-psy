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

#ifndef EbDecHandle_h
#define EbDecHandle_h

#ifdef __cplusplus
extern "C" {
#endif

#include "EbDecStruct.h"
#include "EbDecBlock.h"
#include "EbDecProcess.h"
#include "EbCabacContextModel.h"
#include "Av1Common.h"
#include "EbThreads.h"

/* This value is set to 72 to make
   DEC_PAD_VALUE a multiple of 16. */
#define DYNIMIC_PAD_VALUE  72
#define DEC_PAD_VALUE    (DYNIMIC_PAD_VALUE + 8)

/* Maximum number of frames in parallel */
#define DEC_MAX_NUM_FRM_PRLL 1
/** Maximum picture buffers needed **/
#define MAX_PIC_BUFS (REF_FRAMES + 1 + DEC_MAX_NUM_FRM_PRLL)

/** Picture Structure **/
typedef struct EbDecPicBuf {
    uint8_t is_free;

    size_t size;

    /* Number of reference for this frame */
    uint8_t ref_count;

    uint32_t  order_hint;
    uint32_t  ref_order_hints[INTER_REFS_PER_FRAME];
    FrameType frame_type;

    /*!< Height of the frame in luma samples */
    uint16_t frame_width;
    /* Following 4 prms needed for frame_size_with_refs */
    /*!< Height of the frame in luma samples */
    uint16_t frame_height;
    /*!< Render width of the frame in luma samples */
    uint16_t render_width;
    /*!< Render height of the frame in luma samples */
    uint16_t render_height;
    /*!< Width of Upscaled SuperRes */
    uint16_t superres_upscaled_width;

    EbPictureBufferDesc *ps_pic_buf;

    FRAME_CONTEXT final_frm_ctx;

    GlobalMotionParams global_motion[REF_FRAMES];

    /* MV at 8x8 lvl */
    TemporalMvRef *mvs;

    /* seg map */
    uint8_t *          segment_maps;
    SegmentationParams seg_params;

    /* order hint */
    /* film grain */
    AomFilmGrain film_grain_params;

    // Inter frame reference frame delta for loop filter
    int8_t ref_deltas[REF_FRAMES];
    // 0 = ZERO_MV, MV
    int8_t mode_deltas[MAX_MODE_LF_DELTAS];
} EbDecPicBuf;

/* Frame level buffers */
typedef struct CurFrameBuf {
    SBInfo *sb_info;

    BlockModeInfo *mode_info;

    int32_t *coeff[MAX_MB_PLANE];

    TransformInfo_t *trans_info[MAX_MB_PLANE - 1];

    int8_t * cdef_strength;
    int32_t *delta_q;
    int32_t *delta_lf;

    // Loop Restoration Unit
    RestorationUnitInfo *lr_unit[MAX_MB_PLANE];

    /* Tile Map at SB level : TODO. Can be removed? */
    uint8_t *tile_map_sb;

    /*!< Global warp params of current frame */
    EbWarpedMotionParams global_motion_warp[REF_FRAMES];

    DecMtFrameData dec_mt_frame_data;

} CurFrameBuf;

#define FRAME_MI_MAP 1
/* Frame level buffers */
typedef struct FrameMiMap {
    /* SBInfo pointers for entire frame */
    SBInfo **pps_sb_info;
    /* ModeInfo offset wrt it's SB start */
    uint16_t *p_mi_offset;
    /*!< superblock size inlog2 unit */
    uint8_t sb_size_log2;

    int32_t mi_cols_algnsb;
    int32_t mi_rows_algnsb;
    int32_t sb_cols;
    int32_t sb_rows;
    /*  number of MI in SB width,
        is same as number of MI in SB height */
    int32_t num_mis_in_sb_wd;
} FrameMiMap;

/* Master Frame Buf containing all frame level bufs like ModeInfo
       for all the frames in parallel */
typedef struct MasterFrameBuf {
    CurFrameBuf cur_frame_bufs[DEC_MAX_NUM_FRM_PRLL];

    int32_t num_mis_in_sb;

    int32_t sb_cols;
    int32_t sb_rows;

    /* TODO : Should be moved to thread ctxt */
    FrameMiMap frame_mi_map;

    TemporalMvRef *tpl_mvs;
    int32_t        tpl_mvs_size;
    int8_t         ref_frame_side[REF_FRAMES];

} MasterFrameBuf;

/**************************************
 * Component Private Data
 **************************************/
typedef struct EbDecHandle {
    uint32_t size;
    uint32_t dec_cnt;

    /** Num frames in parallel */
    int32_t num_frms_prll;

    /** Flag to signal seq_header done */
    int32_t seq_header_done;

    /** Flag to signal decoder memory init is done */
    int32_t mem_init_done;

    /** Dec Configuration parameters */
    EbSvtAv1DecConfiguration dec_config;

    SeqHeader   seq_header;
    FrameHeader frame_header;

    uint8_t seen_frame_header;
    /* TODO: Move to frmae ctxt ? */
    uint8_t show_existing_frame;
    uint8_t show_frame;
    uint8_t showable_frame; // frame can be used as show existing frame in future

    // Thread Handles

    // Module Contexts
    void *pv_master_parse_ctxt;

    void *pv_dec_mod_ctxt;

    void *pv_lf_ctxt;

    void *pv_lr_ctxt;

    /** Pointer to Picture manager structure **/
    void *pv_pic_mgr;

    // * 'remapped_ref_idx[i - 1]' maps reference type 'i' (range: LAST_FRAME ...
    // EXTREF_FRAME) to a remapped index 'j' (in range: 0 ... REF_FRAMES - 1)
    // * Later, 'cm->ref_frame_map[j]' maps the remapped index 'j' to a pointer to
    // the reference counted buffer structure RefCntBuffer, taken from the buffer
    // pool cm->buffer_pool->frame_bufs.
    //
    // LAST_FRAME,                        ...,      EXTREF_FRAME
    //      |                                           |
    //      v                                           v
    // remapped_ref_idx[LAST_FRAME - 1],  ...,  remapped_ref_idx[EXTREF_FRAME - 1]
    //      |                                           |
    //      v                                           v
    // ref_frame_map[],                   ...,     ref_frame_map[]
    //
    // Note: INTRA_FRAME always refers to the current frame, so there's no need to
    // have a remapped index for the same.
    int32_t remapped_ref_idx[REF_FRAMES];

    struct ScaleFactors ref_scale_factors[REF_FRAMES];
    /*Scale of the current frame with respect to itself.*/
    struct ScaleFactors sf_identity;

    /* TODO:  Move ref_frame_map, remapped_ref_idx, cur_pic_buf and frame_header to a FrameStr! */
    EbDecPicBuf *ref_frame_map[REF_FRAMES];
    // Prepare ref_frame_map for the next frame.
    EbDecPicBuf *next_ref_frame_map[REF_FRAMES];

    /* For Reference frame loading process */
    EbDecPicBuf *prev_frame;
    /* TODO: Move to buffer pool. */
    EbDecPicBuf *cur_pic_buf[DEC_MAX_NUM_FRM_PRLL];

    // Callbacks

    //DPB + MV, ... buf

    /* Master Frame Buf containing all frame level bufs like ModeInfo
       for all the frames in parallel */
    MasterFrameBuf master_frame_buf;

    // Memory Map
    EbMemoryMapEntry *memory_map_init_address;
    EbMemoryMapEntry *memory_map;
    uint32_t          memory_map_index;
    uint64_t          total_lib_memory;
    struct Av1Common  cm;

    // Loop filter frame level flag
    uint8_t is_lf_enabled;

    // Thread Handles
    EbHandle *            decode_thread_handle_array;
    EbBool                start_thread_process;
    EbHandle              thread_semaphore;
    struct DecThreadCtxt *thread_ctxt_pa;

    EbBool is_16bit_pipeline; // internal bit-depth: when equals 1 internal bit-depth is 16bits regardless of the input bit-depth
} EbDecHandle;

/* Thread level context data */
typedef struct DecThreadCtxt {
    /* Unique ID for the thread */
    uint32_t thread_cnt;
    EbHandle thread_semaphore;
    /* Pointer to the decode handle */
    EbDecHandle *dec_handle_ptr;

    /* Pointer to the decode context */
    void *dec_mod_ctxt;

    /* Temporary block level scratch buffer to
       store LR output of [SB_Size x 64] block */
    uint8_t *dst;
} DecThreadCtxt;

#ifdef __cplusplus
}
#endif
#endif // EbEncHandle_h
