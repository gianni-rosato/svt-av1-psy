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

//#include "EbDeblockingFilter_SSE2.h"

#include "EbPredictionUnit.h"
#include "EbNeighborArrays.h"
#include "EbEncDecProcess.h"
#include "EbDlfProcess.h"
#include "EbDeblockingCommon.h"
#include "common_dsp_rtcd.h"
#ifndef EbDeblockingFilter_h
#define EbDeblockingFilter_h
#ifdef __cplusplus
extern "C" {
#endif

typedef enum LpfPickMethod {
    // Try the full image with different values.
    LPF_PICK_FROM_FULL_IMAGE,
    // Try a small portion of the image with different values.
    LPF_PICK_FROM_SUBIMAGE,
    // Estimate the level based on quantizer and frame type
    LPF_PICK_FROM_Q,
    // Pick 0 to disable LPF if LPF was enabled last frame
    LPF_PICK_MINIMAL_LPF
} LpfPickMethod;

/* assorted LoopFilter functions which get used elsewhere */
struct AV1Common;
struct macroblockd;
struct AV1LfSyncData;

void eb_av1_loop_filter_init(PictureControlSet *pcs_ptr);

void loop_filter_sb(EbPictureBufferDesc *frame_buffer, //reconpicture,
                    //Yv12BufferConfig *frame_buffer,
                    PictureControlSet *pcs_ptr, MacroBlockD *xd, int32_t mi_row, int32_t mi_col,
                    int32_t plane_start, int32_t plane_end, uint8_t last_col);

void eb_av1_loop_filter_frame(
        EbPictureBufferDesc *frame_buffer,//reconpicture,
        //Yv12BufferConfig *frame_buffer,
        PictureControlSet *pcs_ptr,
        /*MacroBlockD *xd,*/ int32_t plane_start, int32_t plane_end/*,
        int32_t partial_frame*/);

void eb_av1_pick_filter_level(DlfContext *         context_ptr,
                              EbPictureBufferDesc *srcBuffer, // source input
                              PictureControlSet *pcs_ptr, LpfPickMethod method);

void eb_av1_filter_block_plane_vert(const PictureControlSet *const pcs_ptr,
                                    const MacroBlockD *const xd, const int32_t plane,
                                    const MacroblockdPlane *const plane_ptr, const uint32_t mi_row,
                                    const uint32_t mi_col);

void eb_av1_filter_block_plane_horz(const PictureControlSet *const pcs_ptr,
                                    const MacroBlockD *const xd, const int32_t plane,
                                    const MacroblockdPlane *const plane_ptr, const uint32_t mi_row,
                                    const uint32_t mi_col);

typedef struct LoopFilterWorkerData {
    EbPictureBufferDesc *   frame_buffer; //reconpicture,
    PictureControlSet *     pcs_ptr;
    struct MacroblockdPlane planes[MAX_MB_PLANE];
    MacroBlockD *xd;
} LFWorkerData;

#ifdef __cplusplus
}
#endif
#endif
