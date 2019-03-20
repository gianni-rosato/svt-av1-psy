/*
* Copyright(c) 2019 Intel Corporation
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

#include "EbDeblockingFilter_SSE2.h"

#include "EbPredictionUnit.h"
#include "EbNeighborArrays.h"
#include "EbEncDecProcess.h"
#if FILT_PROC
#include "EbDlfProcess.h"
#endif

#ifndef EbDeblockingFilter_h
#define EbDeblockingFilter_h
#ifdef __cplusplus
extern "C" {
#endif
#define BLK4X4_ADDR_TO_VERTICAL_EDGE_BS_ARRAY_IDX(blk4x4Addr)                   (((blk4x4Addr) & (MAX_LCU_SIZE_IN_4X4BLK - 1)) + (((blk4x4Addr) / MAX_LCU_SIZE_IN_4X4BLK) * MAX_LCU_SIZE_IN_4X4BLK))
#define BLK4X4_ADDR_TO_HORIZONTAL_EDGE_BS_ARRAY_IDX(blk4x4Addr)                 (((blk4x4Addr) & (MAX_LCU_SIZE_IN_4X4BLK - 1)) + (((blk4x4Addr) / MAX_LCU_SIZE_IN_4X4BLK) * MAX_LCU_SIZE_IN_4X4BLK))
#define GET_LUMA_4X4BLK_ADDR(lumaLcuWise4x4BlkPos_x, lumaLcuWise4x4BlkPos_y, logMaxLcuSizeIn4x4blk)           (((lumaLcuWise4x4BlkPos_x)>>2) + (((lumaLcuWise4x4BlkPos_y)>>2) << (logMaxLcuSizeIn4x4blk)))
#define GET_CHROMA_4X4BLK_ADDR(chromaLcuWise2x2BlkPos_x, chromaLcuWise2x2BlkPos_y, logMaxLcuSizeIn4x4blk)     (((chromaLcuWise2x2BlkPos_x)>>1) + (((chromaLcuWise2x2BlkPos_y)>>1) << (logMaxLcuSizeIn4x4blk)))
#define LUMA_SAMPLE_PIC_WISE_LOCATION_TO_QP_ARRAY_IDX(pos_x, pos_y, qp_array_stride)                            (((pos_x) >> LOG_MIN_BLOCK_SIZE) + ((pos_y) >> LOG_MIN_BLOCK_SIZE) * (qp_array_stride))
#define CHROMA_SAMPLE_PIC_WISE_LOCATION_TO_QP_ARRAY_IDX(pos_x, pos_y, qp_array_stride)                          ((2*(pos_x) >> LOG_MIN_BLOCK_SIZE) + (2*(pos_y) >> LOG_MIN_BLOCK_SIZE) * (qp_array_stride))
#define CHECK_MV_COMPONENT_EQUAL_OR_GREATER_THAN_4(pu1Ptr, pu2Ptr, pu1RefList, pu2RefList)                    (                         \
                            EB_ABS_DIFF((pu1Ptr)->mv[(pu1RefList)].x, (pu2Ptr)->mv[(pu2RefList)].x) >= 4 ||         \
                            EB_ABS_DIFF((pu1Ptr)->mv[(pu1RefList)].y, (pu2Ptr)->mv[(pu2RefList)].y) >= 4            \
                            )

    // Precision macros used in the mode decision
#define BIT_ESTIMATE_PRECISION              15
#define LAMBDA_PRECISION                    16
#define COST_PRECISION                      8
#define MD_SHIFT                            (BIT_ESTIMATE_PRECISION + LAMBDA_PRECISION - COST_PRECISION)
#define MD_OFFSET                           (1 << (MD_SHIFT-1))
#define VAR_QP                              1
#define MAX_QP_VALUE_PLUS_INTRA_TC_OFFSET   53
#define BETA_OFFSET_VALUE                   12 // range -12 to 12
#define TC_OFFSET_VALUE                     12//12 // range -12 to 12

#if AV1_LF
    typedef enum {
        // Try the full image with different values.
        LPF_PICK_FROM_FULL_IMAGE,
        // Try a small portion of the image with different values.
        LPF_PICK_FROM_SUBIMAGE,
        // Estimate the level based on quantizer and frame type
        LPF_PICK_FROM_Q,
        // Pick 0 to disable LPF if LPF was enabled last frame
        LPF_PICK_MINIMAL_LPF
    } LPF_PICK_METHOD;
#endif

    void SetQpArrayBasedOnCU(
        PictureControlSet_t *picture_control_set_ptr,          //input parameter
        uint32_t               cuPos_x,                       //input parameter, sample-based horizontal picture-wise locatin of the CU
        uint32_t               cuPos_y,                       //input parameter, sample-based vertical picture-wise locatin of the CU
        uint32_t               cuSizeInMinCuSize,             //input parameter
        uint32_t               cuQp);                         //input parameter, Qp of the CU

    extern void entropySetQpArrayBasedOnCU(
        PictureControlSet_t *picture_control_set_ptr,
        uint32_t               cuPos_x,
        uint32_t               cuPos_y,
        uint32_t               cuWidth,
        uint32_t               cuHeight,
        uint32_t               cuQp);


    /* assorted loopfilter functions which get used elsewhere */
    struct AV1Common;
    struct macroblockd;
    struct AV1LfSyncData;

    void av1_loop_filter_init(PictureControlSet_t *pcsPtr);

    void av1_loop_filter_frame_init(PictureControlSet_t *pcsPtr, int32_t plane_start,
        int32_t plane_end);


    void loop_filter_sb(
        EbPictureBufferDesc_t *frame_buffer,//reconpicture,
        //Yv12BufferConfig *frame_buffer,
        PictureControlSet_t *pcsPtr,
        MacroBlockD *xd, int32_t mi_row, int32_t mi_col,
        int32_t plane_start, int32_t plane_end,
        uint8_t LastCol);

    void av1_loop_filter_frame(
        EbPictureBufferDesc_t *frame_buffer,//reconpicture,
        //Yv12BufferConfig *frame_buffer,
        PictureControlSet_t *pcsPtr,
        /*MacroBlockD *xd,*/ int32_t plane_start, int32_t plane_end/*,
        int32_t partial_frame*/);

    void av1_pick_filter_level(
#if FILT_PROC
        DlfContext_t            *context_ptr,
#else
        EncDecContext_t         *context_ptr,
#endif
        EbPictureBufferDesc_t   *srcBuffer, // source input
        PictureControlSet_t     *pcsPtr,
        LPF_PICK_METHOD          method);


    void av1_filter_block_plane_vert(
        const PictureControlSet_t *const  pcsPtr,
        const MacroBlockD *const xd,
        const int32_t plane,
        const MacroblockdPlane *const plane_ptr,
        const uint32_t mi_row, const uint32_t mi_col);

    void av1_filter_block_plane_horz(
        const PictureControlSet_t *const  pcsPtr,
        const MacroBlockD *const xd, const int32_t plane,
        const MacroblockdPlane *const plane_ptr,
        const uint32_t mi_row, const uint32_t mi_col);

    typedef struct LoopFilterWorkerData {
        EbPictureBufferDesc_t *frame_buffer;//reconpicture,
        //Yv12BufferConfig *frame_buffer;
        PictureControlSet_t *pcsPtr;
        struct MacroblockdPlane planes[MAX_MB_PLANE];
        // TODO(Ranjit): When the filter functions are modified to use xd->lossless
        // add lossless as a member here.
        MacroBlockD *xd;
    } LFWorkerData;



    void aom_highbd_lpf_horizontal_14_sse2(uint16_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int32_t bd);
#define aom_highbd_lpf_horizontal_14 aom_highbd_lpf_horizontal_14_sse2

    void aom_highbd_lpf_horizontal_14_dual_sse2(uint16_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int32_t bd);
    void aom_highbd_lpf_horizontal_14_dual_avx2(uint16_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int32_t bd);
    //RTCD_EXTERN void(*aom_highbd_lpf_horizontal_14_dual)(uint16_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int32_t bd);

    void aom_highbd_lpf_horizontal_4_c(uint16_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int32_t bd);
    void aom_highbd_lpf_horizontal_4_sse2(uint16_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int32_t bd);
#define aom_highbd_lpf_horizontal_4 aom_highbd_lpf_horizontal_4_sse2

    void aom_highbd_lpf_horizontal_4_dual_sse2(uint16_t *s, int32_t pitch, const uint8_t *blimit0, const uint8_t *limit0, const uint8_t *thresh0, const uint8_t *blimit1, const uint8_t *limit1, const uint8_t *thresh1, int32_t bd);
    void aom_highbd_lpf_horizontal_4_dual_avx2(uint16_t *s, int32_t pitch, const uint8_t *blimit0, const uint8_t *limit0, const uint8_t *thresh0, const uint8_t *blimit1, const uint8_t *limit1, const uint8_t *thresh1, int32_t bd);
    //RTCD_EXTERN void(*aom_highbd_lpf_horizontal_4_dual)(uint16_t *s, int32_t pitch, const uint8_t *blimit0, const uint8_t *limit0, const uint8_t *thresh0, const uint8_t *blimit1, const uint8_t *limit1, const uint8_t *thresh1, int32_t bd);

    void aom_highbd_lpf_horizontal_6_sse2(uint16_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int32_t bd);
#define aom_highbd_lpf_horizontal_6 aom_highbd_lpf_horizontal_6_sse2

    void aom_highbd_lpf_horizontal_8_c(uint16_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int32_t bd);
    void aom_highbd_lpf_horizontal_8_sse2(uint16_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int32_t bd);
#define aom_highbd_lpf_horizontal_8 aom_highbd_lpf_horizontal_8_sse2

    void aom_highbd_lpf_horizontal_8_dual_sse2(uint16_t *s, int32_t pitch, const uint8_t *blimit0, const uint8_t *limit0, const uint8_t *thresh0, const uint8_t *blimit1, const uint8_t *limit1, const uint8_t *thresh1, int32_t bd);
    void aom_highbd_lpf_horizontal_8_dual_avx2(uint16_t *s, int32_t pitch, const uint8_t *blimit0, const uint8_t *limit0, const uint8_t *thresh0, const uint8_t *blimit1, const uint8_t *limit1, const uint8_t *thresh1, int32_t bd);
    //RTCD_EXTERN void(*aom_highbd_lpf_horizontal_8_dual)(uint16_t *s, int32_t pitch, const uint8_t *blimit0, const uint8_t *limit0, const uint8_t *thresh0, const uint8_t *blimit1, const uint8_t *limit1, const uint8_t *thresh1, int32_t bd);

    void aom_highbd_lpf_vertical_14_sse2(uint16_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int32_t bd);
#define aom_highbd_lpf_vertical_14 aom_highbd_lpf_vertical_14_sse2

    void aom_highbd_lpf_vertical_14_dual_sse2(uint16_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int32_t bd);
    void aom_highbd_lpf_vertical_14_dual_avx2(uint16_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int32_t bd);
    //RTCD_EXTERN void(*aom_highbd_lpf_vertical_14_dual)(uint16_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int32_t bd);

    void aom_highbd_lpf_vertical_4_c(uint16_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int32_t bd);
    void aom_highbd_lpf_vertical_4_sse2(uint16_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int32_t bd);
#define aom_highbd_lpf_vertical_4 aom_highbd_lpf_vertical_4_sse2

    void aom_highbd_lpf_vertical_4_dual_sse2(uint16_t *s, int32_t pitch, const uint8_t *blimit0, const uint8_t *limit0, const uint8_t *thresh0, const uint8_t *blimit1, const uint8_t *limit1, const uint8_t *thresh1, int32_t bd);
    void aom_highbd_lpf_vertical_4_dual_avx2(uint16_t *s, int32_t pitch, const uint8_t *blimit0, const uint8_t *limit0, const uint8_t *thresh0, const uint8_t *blimit1, const uint8_t *limit1, const uint8_t *thresh1, int32_t bd);
    //RTCD_EXTERN void(*aom_highbd_lpf_vertical_4_dual)(uint16_t *s, int32_t pitch, const uint8_t *blimit0, const uint8_t *limit0, const uint8_t *thresh0, const uint8_t *blimit1, const uint8_t *limit1, const uint8_t *thresh1, int32_t bd);

    void aom_highbd_lpf_vertical_6_sse2(uint16_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int32_t bd);
#define aom_highbd_lpf_vertical_6 aom_highbd_lpf_vertical_6_sse2

    void aom_highbd_lpf_vertical_8_c(uint16_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int32_t bd);
    void aom_highbd_lpf_vertical_8_sse2(uint16_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int32_t bd);
#define aom_highbd_lpf_vertical_8 aom_highbd_lpf_vertical_8_sse2

    void aom_highbd_lpf_vertical_8_dual_sse2(uint16_t *s, int32_t pitch, const uint8_t *blimit0, const uint8_t *limit0, const uint8_t *thresh0, const uint8_t *blimit1, const uint8_t *limit1, const uint8_t *thresh1, int32_t bd);
    void aom_highbd_lpf_vertical_8_dual_avx2(uint16_t *s, int32_t pitch, const uint8_t *blimit0, const uint8_t *limit0, const uint8_t *thresh0, const uint8_t *blimit1, const uint8_t *limit1, const uint8_t *thresh1, int32_t bd);
    //RTCD_EXTERN void(*aom_highbd_lpf_vertical_8_dual)(uint16_t *s, int32_t pitch, const uint8_t *blimit0, const uint8_t *limit0, const uint8_t *thresh0, const uint8_t *blimit1, const uint8_t *limit1, const uint8_t *thresh1, int32_t bd);


    void aom_lpf_horizontal_14_sse2(uint8_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh);
#define aom_lpf_horizontal_14 aom_lpf_horizontal_14_sse2

    void aom_lpf_horizontal_14_dual_sse2(uint8_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh);
#define aom_lpf_horizontal_14_dual aom_lpf_horizontal_14_dual_sse2

    void aom_lpf_horizontal_4_c(uint8_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh);
    void aom_lpf_horizontal_4_sse2(uint8_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh);
#define aom_lpf_horizontal_4 aom_lpf_horizontal_4_sse2

#define aom_lpf_horizontal_4_dual aom_lpf_horizontal_4_dual_sse2

    void aom_lpf_horizontal_6_c(uint8_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh);
    void aom_lpf_horizontal_6_sse2(uint8_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh);
#define aom_lpf_horizontal_6 aom_lpf_horizontal_6_sse2

    void aom_lpf_horizontal_8_c(uint8_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh);
    void aom_lpf_horizontal_8_sse2(uint8_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh);
#define aom_lpf_horizontal_8 aom_lpf_horizontal_8_sse2

    void aom_lpf_horizontal_8_dual_c(uint8_t *s, int32_t pitch, const uint8_t *blimit0, const uint8_t *limit0, const uint8_t *thresh0, const uint8_t *blimit1, const uint8_t *limit1, const uint8_t *thresh1);
#define aom_lpf_horizontal_8_dual aom_lpf_horizontal_8_dual_c

    void aom_lpf_vertical_14_sse2(uint8_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh);
#define aom_lpf_vertical_14 aom_lpf_vertical_14_sse2

    void aom_lpf_vertical_14_dual_c(uint8_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh);
    void aom_lpf_vertical_14_dual_sse2(uint8_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh);
#define aom_lpf_vertical_14_dual aom_lpf_vertical_14_dual_c

    void aom_lpf_vertical_4_c(uint8_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh);
    void aom_lpf_vertical_4_sse2(uint8_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh);
#define aom_lpf_vertical_4 aom_lpf_vertical_4_sse2

    void aom_lpf_vertical_4_dual_c(uint8_t *s, int32_t pitch, const uint8_t *blimit0, const uint8_t *limit0, const uint8_t *thresh0, const uint8_t *blimit1, const uint8_t *limit1, const uint8_t *thresh1);
#define aom_lpf_vertical_4_dual aom_lpf_vertical_4_dual_c

    void aom_lpf_vertical_6_sse2(uint8_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh);
#define aom_lpf_vertical_6 aom_lpf_vertical_6_sse2

    void aom_lpf_vertical_8_c(uint8_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh);
    void aom_lpf_vertical_8_sse2(uint8_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh);
#define aom_lpf_vertical_8 aom_lpf_vertical_8_sse2

    void aom_lpf_vertical_8_dual_c(uint8_t *s, int32_t pitch, const uint8_t *blimit0, const uint8_t *limit0, const uint8_t *thresh0, const uint8_t *blimit1, const uint8_t *limit1, const uint8_t *thresh1);
#define aom_lpf_vertical_8_dual aom_lpf_vertical_8_dual_c




#ifdef __cplusplus
}
#endif
#endif