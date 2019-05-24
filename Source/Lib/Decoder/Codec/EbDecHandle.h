/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbDecHandle_h
#define EbDecHandle_h

#ifdef __cplusplus
extern "C" {
#endif

#include "EbDecStruct.h"
#include "EbDecBlock.h"

/* Maximum number of frames in parallel */
#define DEC_MAX_NUM_FRM_PRLL    1
/* Number of ref frame buffers needed */
#define DEC_MAX_REF_FRM_BUF   (REF_FRAMES + DEC_MAX_NUM_FRM_PRLL)

/* Frame level buffers */
typedef struct CurFrameBuf {

    SBInfo          *sb_info;

    ModeInfo_t      *mode_info;
    
    int32_t         *luma_coeff;
    int32_t         *chroma_coeff;

    TransformInfo_t *luma_trans_info;
    TransformInfo_t *chroma_trans_info;

    int8_t          *cdef_strength;
    int32_t         *delta_q;
    int32_t         *delta_lf;

    /* Tile Map at SB level : TODO. Can be removed? */
    uint8_t         *tile_map_sb;

} CurFrameBuf;

/* Frame level buffers */
typedef struct FrameMiMap {
    /* For cur SB> Allocated worst case 128x128 SB => 128/4 = 32.
      +1 for 1 top & left 4x4s */
    int16_t      cur_sb_mi_map[33][33];

    /* 2(for 4x4 chroma case) Top SB 4x4 row MI map */
    int16_t      *top_sbrow_mi_map;

    /*  number of MI in SB width, 
        is same as number of MI in SB height */
    int32_t     num_mis_in_sb_wd;
}FrameMiMap;

/* Master Frame Buf containing all frame level bufs like ModeInfo
       for all the frames in parallel */
typedef struct MasterFrameBuf {
    CurFrameBuf     cur_frame_bufs[DEC_MAX_NUM_FRM_PRLL];

    int32_t         num_mis_in_sb;
    //int32_t         mi_cols;
    //int32_t         mi_rows;

    int32_t         sb_cols;
    int32_t         sb_rows;

    /* TODO : Should be moved to thread ctxt */
    FrameMiMap      frame_mi_map;
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
    EbSvtAv1DecConfiguration    dec_config;

    SeqHeader   seq_header;
    FrameHeader frame_header;

    uint8_t seen_frame_header;
    uint8_t show_existing_frame;

    // Thread Handles

    // Module Contexts
    void   *pv_parse_ctxt;

    void   *pv_dec_mod_ctxt;
    

    // Callbacks

    //DPB + MV, ... buf

    /* Master Frame Buf containing all frame level bufs like ModeInfo
       for all the frames in parallel */
    MasterFrameBuf master_frame_buf;

    /* TODO: Move to buffer pool. */
    EbPictureBufferDesc   *recon_picture_buf[DEC_MAX_NUM_FRM_PRLL];

    // Memory Map
#if MEM_MAP_OPT
    EbMemoryMapEntry            *memory_map_init_address;
#endif
    EbMemoryMapEntry            *memory_map;
    uint32_t                     memory_map_index;
    uint64_t                     total_lib_memory;

}EbDecHandle;

#ifdef __cplusplus
    }
#endif
#endif // EbEncHandle_h