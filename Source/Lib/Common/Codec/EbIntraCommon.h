/*
* Copyright(c) 2019 Intel Corporation
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

#ifndef EbIntraCommon_h
#define EbIntraCommon_h


#ifdef __cplusplus
extern "C" {
#endif

typedef struct CflCtx {
    // Q3 reconstructed luma pixels (only Q2 is required, but Q3 is used to avoid
    // shifts)
    EB_ALIGN(64) int16_t recon_buf_q3[CFL_BUF_SQUARE];

    // Height and width currently used in the CfL prediction buffer.
    int32_t buf_height, buf_width;

    int32_t are_parameters_computed;

    // Chroma subsampling
    int32_t subsampling_x, subsampling_y;
} CflCtx;

#ifdef __cplusplus
}
#endif
#endif // EbIntraCommon_h
