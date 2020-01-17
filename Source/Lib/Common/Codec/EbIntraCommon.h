/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
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
