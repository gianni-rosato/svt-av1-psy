/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "EbEncCdef.h"
#include <stdint.h>
#include "aom_dsp_rtcd.h"
#include "EbLog.h"

extern int16_t svt_av1_ac_quant_q3(int32_t qindex, int32_t delta, AomBitDepth bit_depth);

#include "EbRateDistortionCost.h"

static INLINE uint64_t dist_8x8_16bit_c(const uint16_t *src, const uint16_t *dst,
                                        const int32_t dstride, const int32_t coeff_shift) {
    uint64_t svar   = 0;
    uint64_t dvar   = 0;
    uint64_t sum_s  = 0;
    uint64_t sum_d  = 0;
    uint64_t sum_s2 = 0;
    uint64_t sum_d2 = 0;
    uint64_t sum_sd = 0;
    int32_t  i, j;
    for (i = 0; i < 8; i++) {
        for (j = 0; j < 8; j++) {
            sum_s += src[8 * i + j];
            sum_d += dst[i * dstride + j];
            sum_s2 += src[8 * i + j] * src[8 * i + j];
            sum_d2 += dst[i * dstride + j] * dst[i * dstride + j];
            sum_sd += src[8 * i + j] * dst[i * dstride + j];
        }
    }
    /* Compute the variance -- the calculation cannot go negative. */
    svar = sum_s2 - ((sum_s * sum_s + 32) >> 6);
    dvar = sum_d2 - ((sum_d * sum_d + 32) >> 6);
    return (uint64_t)floor(.5 + (sum_d2 + sum_s2 - 2 * sum_sd) * .5 *
                                (svar + dvar + (400 << 2 * coeff_shift)) /
                                (sqrt((20000 << 4 * coeff_shift) + svar * (double)dvar)));
}

static INLINE uint64_t mse_8_16bit(const uint16_t *src, const uint16_t *dst, const int32_t dstride,
                                   const int32_t height) {
    uint64_t sum = 0;
    int32_t  i, j;
    for (i = 0; i < height; i++) {
        for (j = 0; j < 8; j++) {
            int32_t e = dst[i * dstride + j] - src[8 * i + j];
            sum += e * e;
        }
    }
    return sum;
}

static INLINE uint64_t mse_4_16bit_c(const uint16_t *src, const uint16_t *dst,
                                     const int32_t dstride, const int32_t height) {
    uint64_t sum = 0;
    int32_t  i, j;
    for (i = 0; i < height; i++) {
        for (j = 0; j < 4; j++) {
            int32_t e = dst[i * dstride + j] - src[4 * i + j];
            sum += e * e;
        }
    }
    return sum;
}

static INLINE uint64_t dist_8x8_8bit_c(const uint8_t *src, const uint8_t *dst,
                                       const int32_t dstride, const int32_t coeff_shift) {
    uint64_t svar   = 0;
    uint64_t dvar   = 0;
    uint64_t sum_s  = 0;
    uint64_t sum_d  = 0;
    uint64_t sum_s2 = 0;
    uint64_t sum_d2 = 0;
    uint64_t sum_sd = 0;
    int32_t  i, j;
    for (i = 0; i < 8; i++) {
        for (j = 0; j < 8; j++) {
            sum_s += src[8 * i + j];
            sum_d += dst[i * dstride + j];
            sum_s2 += src[8 * i + j] * src[8 * i + j];
            sum_d2 += dst[i * dstride + j] * dst[i * dstride + j];
            sum_sd += src[8 * i + j] * dst[i * dstride + j];
        }
    }
    /* Compute the variance -- the calculation cannot go negative. */
    svar = sum_s2 - ((sum_s * sum_s + 32) >> 6);
    dvar = sum_d2 - ((sum_d * sum_d + 32) >> 6);
    return (uint64_t)floor(.5 + (sum_d2 + sum_s2 - 2 * sum_sd) * .5 *
                                (svar + dvar + (400 << 2 * coeff_shift)) /
                                (sqrt((20000 << 4 * coeff_shift) + svar * (double)dvar)));
}

static INLINE uint64_t mse_8_8bit(const uint8_t *src, const uint8_t *dst, const int32_t dstride,
                                  const int32_t height) {
    uint64_t sum = 0;
    int32_t  i, j;
    for (i = 0; i < height; i++) {
        for (j = 0; j < 8; j++) {
            int32_t e = dst[i * dstride + j] - src[8 * i + j];
            sum += e * e;
        }
    }
    return sum;
}

static INLINE uint64_t mse_4_8bit_c(const uint8_t *src, const uint8_t *dst, const int32_t dstride,
                                    const int32_t height) {
    uint64_t sum = 0;
    int32_t  i, j;
    for (i = 0; i < height; i++) {
        for (j = 0; j < 4; j++) {
            int32_t e = dst[i * dstride + j] - src[4 * i + j];
            sum += e * e;
        }
    }
    return sum;
}


/* Compute MSE only on the blocks we filtered. */
uint64_t compute_cdef_dist_c(const uint16_t *dst, int32_t dstride, const uint16_t *src,
                             const CdefList *dlist, int32_t cdef_count, BlockSize bsize,
                             int32_t coeff_shift, int32_t pli) {
    uint64_t sum = 0;
    int32_t  bi, bx, by;
    if (bsize == BLOCK_8X8) {
        for (bi = 0; bi < cdef_count; bi++) {
            by = dlist[bi].by;
            bx = dlist[bi].bx;
            if (pli == 0) {
                sum += dist_8x8_16bit_c(&src[bi << (3 + 3)],
                                        &dst[(by << 3) * dstride + (bx << 3)],
                                        dstride,
                                        coeff_shift);
            } else
                sum += mse_8_16bit(
                        &src[bi << (3 + 3)], &dst[(by << 3) * dstride + (bx << 3)], dstride, 8);
        }
    } else if (bsize == BLOCK_4X8) {
        for (bi = 0; bi < cdef_count; bi++) {
            by = dlist[bi].by;
            bx = dlist[bi].bx;
            sum += mse_4_16bit_c(
                    &src[bi << (3 + 2)], &dst[(by << 3) * dstride + (bx << 2)], dstride, 8);
        }
    } else if (bsize == BLOCK_8X4) {
        for (bi = 0; bi < cdef_count; bi++) {
            by = dlist[bi].by;
            bx = dlist[bi].bx;
            sum +=
                    mse_8_16bit(&src[bi << (2 + 3)], &dst[(by << 2) * dstride + (bx << 3)], dstride, 4);
        }
    } else {
        assert(bsize == BLOCK_4X4);
        for (bi = 0; bi < cdef_count; bi++) {
            by = dlist[bi].by;
            bx = dlist[bi].bx;
            sum += mse_4_16bit_c(
                    &src[bi << (2 + 2)], &dst[(by << 2) * dstride + (bx << 2)], dstride, 4);
        }
    }
    return sum >> 2 * coeff_shift;
}

uint64_t compute_cdef_dist_8bit_c(const uint8_t *dst8, int32_t dstride, const uint8_t *src8,
                                  const CdefList *dlist, int32_t cdef_count, BlockSize bsize,
                                  int32_t coeff_shift, int32_t pli) {
    uint64_t sum = 0;
    int32_t  bi, bx, by;
    if (bsize == BLOCK_8X8) {
        for (bi = 0; bi < cdef_count; bi++) {
            by = dlist[bi].by;
            bx = dlist[bi].bx;
            if (pli == 0) {
                sum += dist_8x8_8bit_c(&src8[bi << (3 + 3)],
                                       &dst8[(by << 3) * dstride + (bx << 3)],
                                       dstride,
                                       coeff_shift);
            } else
                sum += mse_8_8bit(
                        &src8[bi << (3 + 3)], &dst8[(by << 3) * dstride + (bx << 3)], dstride, 8);
        }
    } else if (bsize == BLOCK_4X8) {
        for (bi = 0; bi < cdef_count; bi++) {
            by = dlist[bi].by;
            bx = dlist[bi].bx;
            sum += mse_4_8bit_c(
                    &src8[bi << (3 + 2)], &dst8[(by << 3) * dstride + (bx << 2)], dstride, 8);
        }
    } else if (bsize == BLOCK_8X4) {
        for (bi = 0; bi < cdef_count; bi++) {
            by = dlist[bi].by;
            bx = dlist[bi].bx;
            sum += mse_8_8bit(
                    &src8[bi << (2 + 3)], &dst8[(by << 2) * dstride + (bx << 3)], dstride, 4);
        }
    } else {
        assert(bsize == BLOCK_4X4);
        for (bi = 0; bi < cdef_count; bi++) {
            by = dlist[bi].by;
            bx = dlist[bi].bx;
            sum += mse_4_8bit_c(
                    &src8[bi << (2 + 2)], &dst8[(by << 2) * dstride + (bx << 2)], dstride, 4);
        }
    }
    return sum >> 2 * coeff_shift;
}


int32_t get_cdef_gi_step(int8_t cdef_level) {
    int32_t gi_step = cdef_level == 5
                        ? 1
                        : cdef_level == 4
                                ? 4
                                : cdef_level == 3 ? 8 : cdef_level == 2 ? 16 : 64;
    return gi_step;
}

int32_t svt_sb_all_skip(PictureControlSet *pcs_ptr, const Av1Common *const cm, int32_t mi_row,
                        int32_t mi_col) {
    int32_t maxc, maxr;
    int32_t skip = 1;
    maxc         = cm->mi_cols - mi_col;
    maxr         = cm->mi_rows - mi_row;

    maxr = AOMMIN(maxr, MI_SIZE_64X64);
    maxc = AOMMIN(maxc, MI_SIZE_64X64);

    for (int32_t r = 0; r < maxr; r++) {
        for (int32_t c = 0; c < maxc; c++) {
            skip = skip && pcs_ptr->mi_grid_base[(mi_row + r) * pcs_ptr->mi_stride + mi_col + c]
                               ->mbmi.block_mi.skip;
        }
    }
    return skip;
}

static int32_t is_8x8_block_skip(ModeInfo **grid, int32_t mi_row, int32_t mi_col,
                                 int32_t mi_stride) {
    int32_t is_skip = 1;
    for (int32_t r = 0; r < mi_size_high[BLOCK_8X8]; ++r)
        for (int32_t c = 0; c < mi_size_wide[BLOCK_8X8]; ++c)
            is_skip &= (int32_t)(grid[(mi_row + r) * mi_stride + (mi_col + c)]->mbmi.block_mi.skip);

    return is_skip;
}

int32_t svt_sb_compute_cdef_list(PictureControlSet *pcs_ptr, const Av1Common *const cm,
                                 int32_t mi_row, int32_t mi_col, CdefList *dlist, BlockSize bs) {
    //MbModeInfo **grid = cm->mi_grid_visible;
    ModeInfo **grid = pcs_ptr->mi_grid_base;

    int32_t maxc = cm->mi_cols - mi_col;
    int32_t maxr = cm->mi_rows - mi_row;

    if (bs == BLOCK_128X128 || bs == BLOCK_128X64)
        maxc = AOMMIN(maxc, MI_SIZE_128X128);
    else
        maxc = AOMMIN(maxc, MI_SIZE_64X64);
    if (bs == BLOCK_128X128 || bs == BLOCK_64X128)
        maxr = AOMMIN(maxr, MI_SIZE_128X128);
    else
        maxr = AOMMIN(maxr, MI_SIZE_64X64);

    const int32_t r_step  = mi_size_high[BLOCK_8X8];
    const int32_t c_step  = mi_size_wide[BLOCK_8X8];
    const int32_t r_shift = (r_step == 2);
    const int32_t c_shift = (c_step == 2);

    assert(r_step == 1 || r_step == 2);
    assert(c_step == 1 || c_step == 2);

    int32_t count = 0;

    for (int32_t r = 0; r < maxr; r += r_step) {
        for (int32_t c = 0; c < maxc; c += c_step) {
            if (!is_8x8_block_skip(grid, mi_row + r, mi_col + c, pcs_ptr->mi_stride)) {
                dlist[count].by   = (uint8_t)(r >> r_shift);
                dlist[count].bx   = (uint8_t)(c >> c_shift);
                dlist[count].skip = 0;
                count++;
            }
        }
    }
    return count;
}

void svt_av1_cdef_frame(EncDecContext *context_ptr, SequenceControlSet *scs_ptr,
                        PictureControlSet *pCs) {
    (void)context_ptr;

    struct PictureParentControlSet *ppcs    = pCs->parent_pcs_ptr;
    Av1Common *                     cm      = ppcs->av1_cm;
    FrameHeader *                   frm_hdr = &ppcs->frm_hdr;

    EbPictureBufferDesc *recon_picture_ptr;

    if (ppcs->is_used_as_reference_flag == EB_TRUE)
        recon_picture_ptr =
            ((EbReferenceObject *)pCs->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
                ->reference_picture;
    else
        recon_picture_ptr = pCs->recon_picture_ptr;

    EbByte recon_buffer_y =
        &((recon_picture_ptr->buffer_y)[recon_picture_ptr->origin_x +
                                        recon_picture_ptr->origin_y * recon_picture_ptr->stride_y]);
    EbByte recon_buffer_cb =
        &((recon_picture_ptr
               ->buffer_cb)[recon_picture_ptr->origin_x / 2 +
                            recon_picture_ptr->origin_y / 2 * recon_picture_ptr->stride_cb]);
    EbByte recon_buffer_cr =
        &((recon_picture_ptr
               ->buffer_cr)[recon_picture_ptr->origin_x / 2 +
                            recon_picture_ptr->origin_y / 2 * recon_picture_ptr->stride_cr]);

    const int32_t num_planes = av1_num_planes(&scs_ptr->seq_header.color_config);
    DECLARE_ALIGNED(16, uint16_t, src[CDEF_INBUF_SIZE]);
    uint16_t *linebuf[3];
    uint16_t *colbuf[3];
    CdefList  dlist[MI_SIZE_64X64 * MI_SIZE_64X64];
    uint8_t * row_cdef, *prev_row_cdef, *curr_row_cdef;
    int32_t   cdef_count;
    int32_t   dir[CDEF_NBLOCKS][CDEF_NBLOCKS] = {{0}};
    int32_t   var[CDEF_NBLOCKS][CDEF_NBLOCKS] = {{0}};
    int32_t   mi_wide_l2[3];
    int32_t   mi_high_l2[3];
    int32_t   xdec[3];
    int32_t   ydec[3];
    int32_t coeff_shift = AOMMAX(scs_ptr->static_config.encoder_bit_depth /*cm->bit_depth*/ - 8, 0);
    const int32_t nvfb  = (cm->mi_rows + MI_SIZE_64X64 - 1) / MI_SIZE_64X64;
    const int32_t nhfb  = (cm->mi_cols + MI_SIZE_64X64 - 1) / MI_SIZE_64X64;
    //svt_av1_setup_dst_planes(xd->plane, cm->seq_params.sb_size, frame, 0, 0, 0, num_planes);
    row_cdef = (uint8_t *)svt_aom_malloc(sizeof(*row_cdef) * (nhfb + 2) * 2);
    assert(row_cdef != NULL);
    memset(row_cdef, 1, sizeof(*row_cdef) * (nhfb + 2) * 2);
    prev_row_cdef = row_cdef + 1;
    curr_row_cdef = prev_row_cdef + nhfb + 2;
    for (int32_t pli = 0; pli < num_planes; pli++) {
        int32_t subsampling_x = (pli == 0) ? 0 : 1;
        int32_t subsampling_y = (pli == 0) ? 0 : 1;

        xdec[pli]       = subsampling_x; //CHKN xd->plane[pli].subsampling_x;
        ydec[pli]       = subsampling_y; //CHKN  xd->plane[pli].subsampling_y;
        mi_wide_l2[pli] = MI_SIZE_LOG2 - subsampling_x; //CHKN xd->plane[pli].subsampling_x;
        mi_high_l2[pli] = MI_SIZE_LOG2 - subsampling_y; //CHKN xd->plane[pli].subsampling_y;
    }

    const int32_t stride = (cm->mi_cols << MI_SIZE_LOG2) + 2 * CDEF_HBORDER;
    for (int32_t pli = 0; pli < num_planes; pli++) {
        linebuf[pli] = (uint16_t *)svt_aom_malloc(sizeof(*linebuf) * CDEF_VBORDER * stride);
        colbuf[pli]  = (uint16_t *)svt_aom_malloc(
            sizeof(*colbuf) * ((CDEF_BLOCKSIZE << mi_high_l2[pli]) + 2 * CDEF_VBORDER) *
            CDEF_HBORDER);
    }

    for (int32_t fbr = 0; fbr < nvfb; fbr++) {
        for (int32_t pli = 0; pli < num_planes; pli++) {
            const int32_t block_height = (MI_SIZE_64X64 << mi_high_l2[pli]) + 2 * CDEF_VBORDER;
            fill_rect(colbuf[pli], CDEF_HBORDER, block_height, CDEF_HBORDER, CDEF_VERY_LARGE);
        }

        int32_t cdef_left = 1;
        for (int32_t fbc = 0; fbc < nhfb; fbc++) {
            int32_t level, sec_strength;
            int32_t uv_level, uv_sec_strength;
            int32_t nhb, nvb;
            int32_t cstart     = 0;
            curr_row_cdef[fbc] = 0;

            //WAHT IS THIS  ?? CHKN -->for
            if (pCs->mi_grid_base[MI_SIZE_64X64 * fbr * cm->mi_stride + MI_SIZE_64X64 * fbc] ==
                    NULL ||
                pCs->mi_grid_base[MI_SIZE_64X64 * fbr * cm->mi_stride + MI_SIZE_64X64 * fbc]
                        ->mbmi.cdef_strength == -1) {
                cdef_left = 0;
                SVT_LOG("\n\n\nCDEF ERROR: Skipping Current FB\n\n\n");
                continue;
            }

            if (!cdef_left)
                cstart =
                    -CDEF_HBORDER; //CHKN if the left block has not been filtered, then we can use samples on the left as input.

            nhb = AOMMIN(MI_SIZE_64X64, cm->mi_cols - MI_SIZE_64X64 * fbc);
            nvb = AOMMIN(MI_SIZE_64X64, cm->mi_rows - MI_SIZE_64X64 * fbr);
            int32_t frame_top, frame_left, frame_bottom, frame_right;

            int32_t mi_row = MI_SIZE_64X64 * fbr;
            int32_t mi_col = MI_SIZE_64X64 * fbc;
            // for the current filter block, it's top left corner mi structure (mi_tl)
            // is first accessed to check whether the top and left boundaries are
            // frame boundaries. Then bottom-left and top-right mi structures are
            // accessed to check whether the bottom and right boundaries
            // (respectively) are frame boundaries.
            //
            // Note that we can't just check the bottom-right mi structure - eg. if
            // we're at the right-hand edge of the frame but not the bottom, then
            // the bottom-right mi is NULL but the bottom-left is not.
            frame_top  = (mi_row == 0) ? 1 : 0;
            frame_left = (mi_col == 0) ? 1 : 0;

            if (fbr != nvfb - 1)
                frame_bottom = (mi_row + MI_SIZE_64X64 == cm->mi_rows) ? 1 : 0;
            else
                frame_bottom = 1;

            if (fbc != nhfb - 1)
                frame_right = (mi_col + MI_SIZE_64X64 == cm->mi_cols) ? 1 : 0;
            else
                frame_right = 1;

            const int32_t mbmi_cdef_strength =
                pCs->mi_grid_base[MI_SIZE_64X64 * fbr * cm->mi_stride + MI_SIZE_64X64 * fbc]
                    ->mbmi.cdef_strength;
            level = frm_hdr->cdef_params.cdef_y_strength[mbmi_cdef_strength] / CDEF_SEC_STRENGTHS;
            sec_strength =
                frm_hdr->cdef_params.cdef_y_strength[mbmi_cdef_strength] % CDEF_SEC_STRENGTHS;
            sec_strength += sec_strength == 3;
            uv_level =
                frm_hdr->cdef_params.cdef_uv_strength[mbmi_cdef_strength] / CDEF_SEC_STRENGTHS;
            uv_sec_strength =
                frm_hdr->cdef_params.cdef_uv_strength[mbmi_cdef_strength] % CDEF_SEC_STRENGTHS;
            uv_sec_strength += uv_sec_strength == 3;
            if ((level == 0 && sec_strength == 0 && uv_level == 0 && uv_sec_strength == 0) ||
                (cdef_count = svt_sb_compute_cdef_list(
                     pCs, cm, fbr * MI_SIZE_64X64, fbc * MI_SIZE_64X64, dlist, BLOCK_64X64)) == 0) {
                cdef_left = 0;
                continue;
            }

            curr_row_cdef[fbc] = 1;
            for (int32_t pli = 0; pli < num_planes; pli++) {
                int32_t coffset;
                int32_t rend, cend;
                int32_t pri_damping = frm_hdr->cdef_params.cdef_damping;
                int32_t sec_damping = pri_damping;
                int32_t hsize       = nhb << mi_wide_l2[pli];
                int32_t vsize       = nvb << mi_high_l2[pli];

                if (pli) {
                    level        = uv_level;
                    sec_strength = uv_sec_strength;
                }

                if (fbc == nhfb - 1)
                    cend = hsize;
                else
                    cend = hsize + CDEF_HBORDER;

                if (fbr == nvfb - 1)
                    rend = vsize;
                else
                    rend = vsize + CDEF_VBORDER;

                coffset = fbc * MI_SIZE_64X64 << mi_wide_l2[pli];
                if (fbc == nhfb - 1) {
                    /* On the last superblock column, fill in the right border with
                       CDEF_VERY_LARGE to avoid filtering with the outside. */
                    fill_rect(&src[cend + CDEF_HBORDER],
                              CDEF_BSTRIDE,
                              rend + CDEF_VBORDER,
                              hsize + CDEF_HBORDER - cend,
                              CDEF_VERY_LARGE);
                }
                if (fbr == nvfb - 1) {
                    /* On the last superblock row, fill in the bottom border with
                       CDEF_VERY_LARGE to avoid filtering with the outside. */
                    fill_rect(&src[(rend + CDEF_VBORDER) * CDEF_BSTRIDE],
                              CDEF_BSTRIDE,
                              CDEF_VBORDER,
                              hsize + 2 * CDEF_HBORDER,
                              CDEF_VERY_LARGE);
                }

                uint8_t *rec_buff   = 0;
                uint32_t rec_stride = 0;

                switch (pli) {
                case 0:
                    rec_buff   = recon_buffer_y;
                    rec_stride = recon_picture_ptr->stride_y;
                    break;
                case 1:
                    rec_buff   = recon_buffer_cb;
                    rec_stride = recon_picture_ptr->stride_cb;

                    break;
                case 2:
                    rec_buff   = recon_buffer_cr;
                    rec_stride = recon_picture_ptr->stride_cr;
                    break;
                }

                /* Copy in the pixels we need from the current superblock for
                   deringing.*/
                copy_sb8_16(&src[CDEF_VBORDER * CDEF_BSTRIDE + CDEF_HBORDER + cstart],
                            CDEF_BSTRIDE,
                            rec_buff,
                            (MI_SIZE_64X64 << mi_high_l2[pli]) * fbr,
                            coffset + cstart,
                            rec_stride,
                            rend,
                            cend - cstart);
                if (!prev_row_cdef[fbc]) {
                    copy_sb8_16(&src[CDEF_HBORDER],
                                CDEF_BSTRIDE,
                                rec_buff,
                                (MI_SIZE_64X64 << mi_high_l2[pli]) * fbr - CDEF_VBORDER,
                                coffset,
                                rec_stride,
                                CDEF_VBORDER,
                                hsize);
                } else if (fbr > 0) {
                    copy_rect(&src[CDEF_HBORDER],
                              CDEF_BSTRIDE,
                              &linebuf[pli][coffset],
                              stride,
                              CDEF_VBORDER,
                              hsize);
                } else {
                    fill_rect(
                        &src[CDEF_HBORDER], CDEF_BSTRIDE, CDEF_VBORDER, hsize, CDEF_VERY_LARGE);
                }

                if (!prev_row_cdef[fbc - 1]) {
                    copy_sb8_16(src,
                                CDEF_BSTRIDE,
                                rec_buff,
                                (MI_SIZE_64X64 << mi_high_l2[pli]) * fbr - CDEF_VBORDER,
                                coffset - CDEF_HBORDER,
                                rec_stride,
                                CDEF_VBORDER,
                                CDEF_HBORDER);
                } else if (fbr > 0 && fbc > 0) {
                    copy_rect(src,
                              CDEF_BSTRIDE,
                              &linebuf[pli][coffset - CDEF_HBORDER],
                              stride,
                              CDEF_VBORDER,
                              CDEF_HBORDER);
                } else {
                    fill_rect(src, CDEF_BSTRIDE, CDEF_VBORDER, CDEF_HBORDER, CDEF_VERY_LARGE);
                }

                if (!prev_row_cdef[fbc + 1]) {
                    copy_sb8_16(&src[CDEF_HBORDER + (nhb << mi_wide_l2[pli])],
                                CDEF_BSTRIDE,
                                rec_buff,
                                (MI_SIZE_64X64 << mi_high_l2[pli]) * fbr - CDEF_VBORDER,
                                coffset + hsize,
                                rec_stride,
                                CDEF_VBORDER,
                                CDEF_HBORDER);
                } else if (fbr > 0 && fbc < nhfb - 1) {
                    copy_rect(&src[hsize + CDEF_HBORDER],
                              CDEF_BSTRIDE,
                              &linebuf[pli][coffset + hsize],
                              stride,
                              CDEF_VBORDER,
                              CDEF_HBORDER);
                } else {
                    fill_rect(&src[hsize + CDEF_HBORDER],
                              CDEF_BSTRIDE,
                              CDEF_VBORDER,
                              CDEF_HBORDER,
                              CDEF_VERY_LARGE);
                }

                if (cdef_left) {
                    /* If we deringed the superblock on the left then we need to copy in
                       saved pixels. */
                    copy_rect(src,
                              CDEF_BSTRIDE,
                              colbuf[pli],
                              CDEF_HBORDER,
                              rend + CDEF_VBORDER,
                              CDEF_HBORDER);
                }

                /* Saving pixels in case we need to dering the superblock on the
                    right. */
                if (fbc < nhfb - 1)
                    copy_rect(colbuf[pli],
                              CDEF_HBORDER,
                              src + hsize,
                              CDEF_BSTRIDE,
                              rend + CDEF_VBORDER,
                              CDEF_HBORDER);

                if (fbr < nvfb - 1)
                    copy_sb8_16(&linebuf[pli][coffset],
                                stride,
                                rec_buff,
                                (MI_SIZE_64X64 << mi_high_l2[pli]) * (fbr + 1) - CDEF_VBORDER,
                                coffset,
                                rec_stride,
                                CDEF_VBORDER,
                                hsize);

                if (frame_top) {
                    fill_rect(
                        src, CDEF_BSTRIDE, CDEF_VBORDER, hsize + 2 * CDEF_HBORDER, CDEF_VERY_LARGE);
                }
                if (frame_left) {
                    fill_rect(
                        src, CDEF_BSTRIDE, vsize + 2 * CDEF_VBORDER, CDEF_HBORDER, CDEF_VERY_LARGE);
                }
                if (frame_bottom) {
                    fill_rect(&src[(vsize + CDEF_VBORDER) * CDEF_BSTRIDE],
                              CDEF_BSTRIDE,
                              CDEF_VBORDER,
                              hsize + 2 * CDEF_HBORDER,
                              CDEF_VERY_LARGE);
                }
                if (frame_right) {
                    fill_rect(&src[hsize + CDEF_HBORDER],
                              CDEF_BSTRIDE,
                              vsize + 2 * CDEF_VBORDER,
                              CDEF_HBORDER,
                              CDEF_VERY_LARGE);
                }

                {
                    svt_cdef_filter_fb(
                        &rec_buff[rec_stride * (MI_SIZE_64X64 * fbr << mi_high_l2[pli]) +
                                  (fbc * MI_SIZE_64X64 << mi_wide_l2[pli])],
                        NULL,
                        rec_stride,
                        &src[CDEF_VBORDER * CDEF_BSTRIDE + CDEF_HBORDER],
                        xdec[pli],
                        ydec[pli],
                        dir,
                        NULL,
                        var,
                        pli,
                        dlist,
                        cdef_count,
                        level,
                        sec_strength,
                        pri_damping,
                        sec_damping,
                        coeff_shift);
                }
            }
            cdef_left = 1; //CHKN filtered data is written back directy to recFrame.
        }
        {
            uint8_t *tmp  = prev_row_cdef;
            prev_row_cdef = curr_row_cdef;
            curr_row_cdef = tmp;
        }
    }
    svt_aom_free(row_cdef);
    for (int32_t pli = 0; pli < num_planes; pli++) {
        svt_aom_free(linebuf[pli]);
        svt_aom_free(colbuf[pli]);
    }
}

void av1_cdef_frame16bit(EncDecContext *context_ptr, SequenceControlSet *scs_ptr,
                         PictureControlSet *pCs) {
    (void)context_ptr;
    struct PictureParentControlSet *ppcs    = pCs->parent_pcs_ptr;
    Av1Common *                     cm      = ppcs->av1_cm;
    FrameHeader *                   frm_hdr = &ppcs->frm_hdr;

    EbPictureBufferDesc *recon_picture_ptr;

    if (ppcs->is_used_as_reference_flag == EB_TRUE)
        recon_picture_ptr =
            ((EbReferenceObject *)pCs->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
                ->reference_picture16bit;

    else
        recon_picture_ptr = pCs->recon_picture16bit_ptr;

    uint16_t *recon_buffer_y =
        (uint16_t *)recon_picture_ptr->buffer_y +
        (recon_picture_ptr->origin_x + recon_picture_ptr->origin_y * recon_picture_ptr->stride_y);
    uint16_t *recon_buffer_cb = (uint16_t *)recon_picture_ptr->buffer_cb +
                                (recon_picture_ptr->origin_x / 2 +
                                 recon_picture_ptr->origin_y / 2 * recon_picture_ptr->stride_cb);
    uint16_t *recon_buffer_cr = (uint16_t *)recon_picture_ptr->buffer_cr +
                                (recon_picture_ptr->origin_x / 2 +
                                 recon_picture_ptr->origin_y / 2 * recon_picture_ptr->stride_cr);

    const int32_t num_planes = av1_num_planes(&scs_ptr->seq_header.color_config);
    DECLARE_ALIGNED(16, uint16_t, src[CDEF_INBUF_SIZE]);
    uint16_t *linebuf[3];
    uint16_t *colbuf[3];
    CdefList  dlist[MI_SIZE_64X64 * MI_SIZE_64X64];
    uint8_t * row_cdef, *prev_row_cdef, *curr_row_cdef;
    int32_t   cdef_count;
    int32_t   dir[CDEF_NBLOCKS][CDEF_NBLOCKS] = {{0}};
    int32_t   var[CDEF_NBLOCKS][CDEF_NBLOCKS] = {{0}};
    int32_t   mi_wide_l2[3];
    int32_t   mi_high_l2[3];
    int32_t   xdec[3];
    int32_t   ydec[3];
    int32_t coeff_shift = AOMMAX(scs_ptr->static_config.encoder_bit_depth /*cm->bit_depth*/ - 8, 0);
    const int32_t nvfb  = (cm->mi_rows + MI_SIZE_64X64 - 1) / MI_SIZE_64X64;
    const int32_t nhfb  = (cm->mi_cols + MI_SIZE_64X64 - 1) / MI_SIZE_64X64;
    row_cdef            = (uint8_t *)svt_aom_malloc(sizeof(*row_cdef) * (nhfb + 2) * 2);
    assert(row_cdef);
    memset(row_cdef, 1, sizeof(*row_cdef) * (nhfb + 2) * 2);
    prev_row_cdef = row_cdef + 1;
    curr_row_cdef = prev_row_cdef + nhfb + 2;
    for (int32_t pli = 0; pli < num_planes; pli++) {
        int32_t subsampling_x = (pli == 0) ? 0 : 1;
        int32_t subsampling_y = (pli == 0) ? 0 : 1;

        xdec[pli]       = subsampling_x; //CHKN xd->plane[pli].subsampling_x;
        ydec[pli]       = subsampling_y; //CHKN  xd->plane[pli].subsampling_y;
        mi_wide_l2[pli] = MI_SIZE_LOG2 - subsampling_x; //CHKN xd->plane[pli].subsampling_x;
        mi_high_l2[pli] = MI_SIZE_LOG2 - subsampling_y; //CHKN xd->plane[pli].subsampling_y;
    }

    const int32_t stride = (cm->mi_cols << MI_SIZE_LOG2) + 2 * CDEF_HBORDER;
    for (int32_t pli = 0; pli < num_planes; pli++) {
        linebuf[pli] = (uint16_t *)svt_aom_malloc(sizeof(*linebuf) * CDEF_VBORDER * stride);
        colbuf[pli]  = (uint16_t *)svt_aom_malloc(
            sizeof(*colbuf) * ((CDEF_BLOCKSIZE << mi_high_l2[pli]) + 2 * CDEF_VBORDER) *
            CDEF_HBORDER);
    }

    for (int32_t fbr = 0; fbr < nvfb; fbr++) {
        for (int32_t pli = 0; pli < num_planes; pli++) {
            const int32_t block_height = (MI_SIZE_64X64 << mi_high_l2[pli]) + 2 * CDEF_VBORDER;
            fill_rect(colbuf[pli], CDEF_HBORDER, block_height, CDEF_HBORDER, CDEF_VERY_LARGE);
        }

        int32_t cdef_left = 1;
        for (int32_t fbc = 0; fbc < nhfb; fbc++) {
            int32_t level, sec_strength;
            int32_t uv_level, uv_sec_strength;
            int32_t nhb, nvb;
            int32_t cstart     = 0;
            curr_row_cdef[fbc] = 0;

            //WAHT IS THIS  ?? CHKN -->for
            if (pCs->mi_grid_base[MI_SIZE_64X64 * fbr * cm->mi_stride + MI_SIZE_64X64 * fbc] ==
                    NULL ||
                pCs->mi_grid_base[MI_SIZE_64X64 * fbr * cm->mi_stride + MI_SIZE_64X64 * fbc]
                        ->mbmi.cdef_strength == -1) {
                cdef_left = 0;
                SVT_LOG("\n\n\nCDEF ERROR: Skipping Current FB\n\n\n");
                continue;
            }

            if (!cdef_left)
                cstart =
                    -CDEF_HBORDER; //CHKN if the left block has not been filtered, then we can use samples on the left as input.

            nhb = AOMMIN(MI_SIZE_64X64, cm->mi_cols - MI_SIZE_64X64 * fbc);
            nvb = AOMMIN(MI_SIZE_64X64, cm->mi_rows - MI_SIZE_64X64 * fbr);
            int32_t frame_top, frame_left, frame_bottom, frame_right;

            int32_t mi_row = MI_SIZE_64X64 * fbr;
            int32_t mi_col = MI_SIZE_64X64 * fbc;
            // for the current filter block, it's top left corner mi structure (mi_tl)
            // is first accessed to check whether the top and left boundaries are
            // frame boundaries. Then bottom-left and top-right mi structures are
            // accessed to check whether the bottom and right boundaries
            // (respectively) are frame boundaries.
            //
            // Note that we can't just check the bottom-right mi structure - eg. if
            // we're at the right-hand edge of the frame but not the bottom, then
            // the bottom-right mi is NULL but the bottom-left is not.
            frame_top  = (mi_row == 0) ? 1 : 0;
            frame_left = (mi_col == 0) ? 1 : 0;

            if (fbr != nvfb - 1)
                frame_bottom = (mi_row + MI_SIZE_64X64 == cm->mi_rows) ? 1 : 0;
            else
                frame_bottom = 1;

            if (fbc != nhfb - 1)
                frame_right = (mi_col + MI_SIZE_64X64 == cm->mi_cols) ? 1 : 0;
            else
                frame_right = 1;

            const int32_t mbmi_cdef_strength =
                pCs->mi_grid_base[MI_SIZE_64X64 * fbr * cm->mi_stride + MI_SIZE_64X64 * fbc]
                    ->mbmi.cdef_strength;
            level = frm_hdr->cdef_params.cdef_y_strength[mbmi_cdef_strength] / CDEF_SEC_STRENGTHS;
            sec_strength =
                frm_hdr->cdef_params.cdef_y_strength[mbmi_cdef_strength] % CDEF_SEC_STRENGTHS;
            sec_strength += sec_strength == 3;
            uv_level =
                frm_hdr->cdef_params.cdef_uv_strength[mbmi_cdef_strength] / CDEF_SEC_STRENGTHS;
            uv_sec_strength =
                frm_hdr->cdef_params.cdef_uv_strength[mbmi_cdef_strength] % CDEF_SEC_STRENGTHS;
            uv_sec_strength += uv_sec_strength == 3;
            if ((level == 0 && sec_strength == 0 && uv_level == 0 && uv_sec_strength == 0) ||
                (cdef_count = svt_sb_compute_cdef_list(
                     pCs, cm, fbr * MI_SIZE_64X64, fbc * MI_SIZE_64X64, dlist, BLOCK_64X64)) == 0) {
                cdef_left = 0;
                continue;
            }

            curr_row_cdef[fbc] = 1;
            for (int32_t pli = 0; pli < num_planes; pli++) {
                int32_t coffset;
                int32_t rend, cend;
                int32_t pri_damping = frm_hdr->cdef_params.cdef_damping;
                int32_t sec_damping = pri_damping;
                int32_t hsize       = nhb << mi_wide_l2[pli];
                int32_t vsize       = nvb << mi_high_l2[pli];

                if (pli) {
                    level        = uv_level;
                    sec_strength = uv_sec_strength;
                }

                if (fbc == nhfb - 1)
                    cend = hsize;
                else
                    cend = hsize + CDEF_HBORDER;

                if (fbr == nvfb - 1)
                    rend = vsize;
                else
                    rend = vsize + CDEF_VBORDER;

                coffset = fbc * MI_SIZE_64X64 << mi_wide_l2[pli];
                if (fbc == nhfb - 1) {
                    /* On the last superblock column, fill in the right border with
                    CDEF_VERY_LARGE to avoid filtering with the outside. */
                    fill_rect(&src[cend + CDEF_HBORDER],
                              CDEF_BSTRIDE,
                              rend + CDEF_VBORDER,
                              hsize + CDEF_HBORDER - cend,
                              CDEF_VERY_LARGE);
                }
                if (fbr == nvfb - 1) {
                    /* On the last superblock row, fill in the bottom border with
                    CDEF_VERY_LARGE to avoid filtering with the outside. */
                    fill_rect(&src[(rend + CDEF_VBORDER) * CDEF_BSTRIDE],
                              CDEF_BSTRIDE,
                              CDEF_VBORDER,
                              hsize + 2 * CDEF_HBORDER,
                              CDEF_VERY_LARGE);
                }

                uint16_t *rec_buff   = 0;
                uint32_t  rec_stride = 0;

                switch (pli) {
                case 0:
                    rec_buff   = recon_buffer_y;
                    rec_stride = recon_picture_ptr->stride_y;
                    break;
                case 1:
                    rec_buff   = recon_buffer_cb;
                    rec_stride = recon_picture_ptr->stride_cb;

                    break;
                case 2:
                    rec_buff   = recon_buffer_cr;
                    rec_stride = recon_picture_ptr->stride_cr;
                    break;
                }

                /* Copy in the pixels we need from the current superblock for
                deringing.*/

                copy_sb16_16(&src[CDEF_VBORDER * CDEF_BSTRIDE + CDEF_HBORDER + cstart],
                             CDEF_BSTRIDE,
                             rec_buff,
                             (MI_SIZE_64X64 << mi_high_l2[pli]) * fbr,
                             coffset + cstart,
                             rec_stride,
                             rend,
                             cend - cstart);

                if (!prev_row_cdef[fbc]) {
                    copy_sb16_16(&src[CDEF_HBORDER],
                                 CDEF_BSTRIDE,
                                 rec_buff,
                                 (MI_SIZE_64X64 << mi_high_l2[pli]) * fbr - CDEF_VBORDER,
                                 coffset,
                                 rec_stride,
                                 CDEF_VBORDER,
                                 hsize);
                } else if (fbr > 0) {
                    copy_rect(&src[CDEF_HBORDER],
                              CDEF_BSTRIDE,
                              &linebuf[pli][coffset],
                              stride,
                              CDEF_VBORDER,
                              hsize);
                } else {
                    fill_rect(
                        &src[CDEF_HBORDER], CDEF_BSTRIDE, CDEF_VBORDER, hsize, CDEF_VERY_LARGE);
                }

                if (!prev_row_cdef[fbc - 1]) {
                    copy_sb16_16(src,
                                 CDEF_BSTRIDE,
                                 rec_buff,
                                 (MI_SIZE_64X64 << mi_high_l2[pli]) * fbr - CDEF_VBORDER,
                                 coffset - CDEF_HBORDER,
                                 rec_stride,
                                 CDEF_VBORDER,
                                 CDEF_HBORDER);
                } else if (fbr > 0 && fbc > 0) {
                    copy_rect(src,
                              CDEF_BSTRIDE,
                              &linebuf[pli][coffset - CDEF_HBORDER],
                              stride,
                              CDEF_VBORDER,
                              CDEF_HBORDER);
                } else {
                    fill_rect(src, CDEF_BSTRIDE, CDEF_VBORDER, CDEF_HBORDER, CDEF_VERY_LARGE);
                }

                if (!prev_row_cdef[fbc + 1]) {
                    copy_sb16_16(&src[CDEF_HBORDER + (nhb << mi_wide_l2[pli])],
                                 CDEF_BSTRIDE,
                                 rec_buff,
                                 (MI_SIZE_64X64 << mi_high_l2[pli]) * fbr - CDEF_VBORDER,
                                 coffset + hsize,
                                 rec_stride,
                                 CDEF_VBORDER,
                                 CDEF_HBORDER);
                } else if (fbr > 0 && fbc < nhfb - 1) {
                    copy_rect(&src[hsize + CDEF_HBORDER],
                              CDEF_BSTRIDE,
                              &linebuf[pli][coffset + hsize],
                              stride,
                              CDEF_VBORDER,
                              CDEF_HBORDER);
                } else {
                    fill_rect(&src[hsize + CDEF_HBORDER],
                              CDEF_BSTRIDE,
                              CDEF_VBORDER,
                              CDEF_HBORDER,
                              CDEF_VERY_LARGE);
                }

                if (cdef_left) {
                    /* If we deringed the superblock on the left then we need to copy in
                    saved pixels. */
                    copy_rect(src,
                              CDEF_BSTRIDE,
                              colbuf[pli],
                              CDEF_HBORDER,
                              rend + CDEF_VBORDER,
                              CDEF_HBORDER);
                }

                /* Saving pixels in case we need to dering the superblock on the
                right. */
                if (fbc < nhfb - 1)
                    copy_rect(colbuf[pli],
                              CDEF_HBORDER,
                              src + hsize,
                              CDEF_BSTRIDE,
                              rend + CDEF_VBORDER,
                              CDEF_HBORDER);
                if (fbr < nvfb - 1)
                    copy_sb16_16(&linebuf[pli][coffset],
                                 stride,
                                 rec_buff,
                                 (MI_SIZE_64X64 << mi_high_l2[pli]) * (fbr + 1) - CDEF_VBORDER,
                                 coffset,
                                 rec_stride,
                                 CDEF_VBORDER,
                                 hsize);
                if (frame_top) {
                    fill_rect(
                        src, CDEF_BSTRIDE, CDEF_VBORDER, hsize + 2 * CDEF_HBORDER, CDEF_VERY_LARGE);
                }
                if (frame_left) {
                    fill_rect(
                        src, CDEF_BSTRIDE, vsize + 2 * CDEF_VBORDER, CDEF_HBORDER, CDEF_VERY_LARGE);
                }
                if (frame_bottom) {
                    fill_rect(&src[(vsize + CDEF_VBORDER) * CDEF_BSTRIDE],
                              CDEF_BSTRIDE,
                              CDEF_VBORDER,
                              hsize + 2 * CDEF_HBORDER,
                              CDEF_VERY_LARGE);
                }
                if (frame_right) {
                    fill_rect(&src[hsize + CDEF_HBORDER],
                              CDEF_BSTRIDE,
                              vsize + 2 * CDEF_VBORDER,
                              CDEF_HBORDER,
                              CDEF_VERY_LARGE);
                }

                svt_cdef_filter_fb(NULL,
                                   &rec_buff[rec_stride * (MI_SIZE_64X64 * fbr << mi_high_l2[pli]) +
                                             (fbc * MI_SIZE_64X64 << mi_wide_l2[pli])],
                                   rec_stride,
                                   &src[CDEF_VBORDER * CDEF_BSTRIDE + CDEF_HBORDER],
                                   xdec[pli],
                                   ydec[pli],
                                   dir,
                                   NULL,
                                   var,
                                   pli,
                                   dlist,
                                   cdef_count,
                                   level,
                                   sec_strength,
                                   pri_damping,
                                   sec_damping,
                                   coeff_shift);
            }
            cdef_left = 1; //CHKN filtered data is written back directy to recFrame.
        }
        {
            uint8_t *tmp  = prev_row_cdef;
            prev_row_cdef = curr_row_cdef;
            curr_row_cdef = tmp;
        }
    }
    svt_aom_free(row_cdef);
    for (int32_t pli = 0; pli < num_planes; pli++) {
        svt_aom_free(linebuf[pli]);
        svt_aom_free(colbuf[pli]);
    }
}

///-------search

/* Search for the best strength to add as an option, knowing we
already selected nb_strengths options. */
static uint64_t search_one(int32_t *lev, int32_t nb_strengths, uint64_t mse[][TOTAL_STRENGTHS],
                           int32_t sb_count, int32_t start_gi, int32_t end_gi) {
    uint64_t tot_mse[TOTAL_STRENGTHS];
    const int32_t total_strengths = end_gi;
    int32_t       i, j;
    uint64_t      best_tot_mse = (uint64_t)1 << 63;
    int32_t       best_id      = 0;
    memset(tot_mse, 0, sizeof(tot_mse));
    for (i = 0; i < sb_count; i++) {
        int32_t  gi;
        uint64_t best_mse = (uint64_t)1 << 63;
        /* Find best mse among already selected options. */
        for (gi = 0; gi < nb_strengths; gi++) {
            if (mse[i][lev[gi]] < best_mse) best_mse = mse[i][lev[gi]];
        }
        /* Find best mse when adding each possible new option. */

        for (j = start_gi; j < total_strengths; j++) {
            uint64_t best = best_mse;
            if (mse[i][j] < best) best = mse[i][j];
            tot_mse[j] += best;
        }
    }
    for (j = start_gi; j < total_strengths; j++) {
        if (tot_mse[j] < best_tot_mse) {
            best_tot_mse = tot_mse[j];
            best_id      = j;
        }
    }
    lev[nb_strengths] = best_id;
    return best_tot_mse;
}

/* Search for the best luma+chroma strength to add as an option, knowing we
already selected nb_strengths options. */
uint64_t svt_search_one_dual_c(int *lev0, int *lev1, int nb_strengths,
                               uint64_t (**mse)[TOTAL_STRENGTHS], int sb_count,
                               int start_gi, int end_gi) {
    uint64_t tot_mse[TOTAL_STRENGTHS][TOTAL_STRENGTHS];
    int32_t  i, j;
    uint64_t best_tot_mse = (uint64_t)1 << 63;
    int32_t  best_id0     = 0;
    int32_t  best_id1     = 0;
    const int32_t total_strengths = end_gi;
    memset(tot_mse, 0, sizeof(tot_mse));
    for (i = 0; i < sb_count; i++) {
        int32_t  gi;
        uint64_t best_mse = (uint64_t)1 << 63;
        /* Find best mse among already selected options. */
        for (gi = 0; gi < nb_strengths; gi++) {
            uint64_t curr = mse[0][i][lev0[gi]];
            curr += mse[1][i][lev1[gi]];
            if (curr < best_mse) best_mse = curr;
        }
        /* Find best mse when adding each possible new option. */
        for (j = start_gi; j < total_strengths; j++) {
            int32_t k;
            for (k = start_gi; k < total_strengths; k++) {
                uint64_t best = best_mse;
                uint64_t curr = mse[0][i][j];
                curr += mse[1][i][k];
                if (curr < best) best = curr;
                tot_mse[j][k] += best;
            }
        }
    }

    for (j = start_gi; j < total_strengths; j++) {
        int32_t k;
        for (k = start_gi; k < total_strengths; k++) {
            if (tot_mse[j][k] < best_tot_mse) {
                best_tot_mse = tot_mse[j][k];
                best_id0     = j;
                best_id1     = k;
            }
        }
    }
    lev0[nb_strengths] = best_id0;
    lev1[nb_strengths] = best_id1;
    return best_tot_mse;
}

/* Search for the set of strengths that minimizes mse. */
static uint64_t joint_strength_search(int32_t *best_lev, int32_t nb_strengths,
                                      uint64_t mse[][TOTAL_STRENGTHS], int32_t sb_count,
                                      int32_t start_gi, int32_t end_gi) {
    uint64_t best_tot_mse;
    int32_t  i;
    best_tot_mse = (uint64_t)1 << 63;
    /* Greedy search: add one strength options at a time. */
    for (i = 0; i < nb_strengths; i++)
        best_tot_mse = search_one(best_lev, i, mse, sb_count, start_gi, end_gi);
    /* Trying to refine the greedy search by reconsidering each
    already-selected option. */
    for (i = 0; i < 4 * nb_strengths; i++) {
        int32_t j;
        for (j = 0; j < nb_strengths - 1; j++) best_lev[j] = best_lev[j + 1];
        best_tot_mse =
            search_one(best_lev, nb_strengths - 1, mse, sb_count, start_gi, end_gi);
    }
    return best_tot_mse;
}

/* Search for the set of luma+chroma strengths that minimizes mse. */
static uint64_t joint_strength_search_dual(int32_t *best_lev0, int32_t *best_lev1,
                                           int32_t nb_strengths, uint64_t (**mse)[TOTAL_STRENGTHS],
                                           int32_t sb_count, int32_t start_gi,
                                           int32_t end_gi) {
    uint64_t best_tot_mse;
    int32_t  i;
    best_tot_mse = (uint64_t)1 << 63;
    /* Greedy search: add one strength options at a time. */
    for (i = 0; i < nb_strengths; i++)
        best_tot_mse = svt_search_one_dual(
            best_lev0, best_lev1, i, mse, sb_count, start_gi, end_gi);
    /* Trying to refine the greedy search by reconsidering each
    already-selected option. */
    for (i = 0; i < 4 * nb_strengths; i++) {
        int32_t j;
        for (j = 0; j < nb_strengths - 1; j++) {
            best_lev0[j] = best_lev0[j + 1];
            best_lev1[j] = best_lev1[j + 1];
        }
        best_tot_mse = svt_search_one_dual(
            best_lev0, best_lev1, nb_strengths - 1, mse, sb_count, start_gi, end_gi);
    }
    return best_tot_mse;
}

void finish_cdef_search(EncDecContext *context_ptr, PictureControlSet *pcs_ptr,
                        int32_t selected_strength_cnt[64]) {
    (void)context_ptr;
    struct PictureParentControlSet *ppcs    = pcs_ptr->parent_pcs_ptr;
    FrameHeader *                   frm_hdr = &ppcs->frm_hdr;
    Av1Common *                     cm      = ppcs->av1_cm;
    int32_t                         mi_rows = ppcs->av1_cm->mi_rows;
    int32_t                         mi_cols = ppcs->av1_cm->mi_cols;

    int32_t fbr, fbc;

    int32_t pli;

    uint64_t      best_tot_mse = (uint64_t)1 << 63;
    int32_t       sb_count;
    int32_t       nvfb              = (mi_rows + MI_SIZE_64X64 - 1) / MI_SIZE_64X64;
    int32_t       nhfb              = (mi_cols + MI_SIZE_64X64 - 1) / MI_SIZE_64X64;
    int32_t *     sb_index          = (int32_t *)malloc(nvfb * nhfb * sizeof(*sb_index));
    int32_t *     selected_strength = (int32_t *)malloc(nvfb * nhfb * sizeof(*sb_index));
    int32_t       best_frame_gi_cnt = 0;
    const int32_t total_strengths   = TOTAL_STRENGTHS;
    int32_t       gi_step;
    int32_t       mid_gi;
    int32_t       start_gi;
    int32_t       end_gi;

    assert(sb_index != NULL);
    assert(selected_strength != NULL);

    gi_step = get_cdef_gi_step(ppcs->cdef_level);

    mid_gi = ppcs->cdf_ref_frame_strength;
    start_gi = ppcs->use_ref_frame_cdef_strength && ppcs->cdef_level == 5
                    ? (AOMMAX(0, mid_gi - gi_step))
                    : 0;
    end_gi = ppcs->use_ref_frame_cdef_strength ? AOMMIN(total_strengths, mid_gi + gi_step)
                                               : ppcs->cdef_level == 5 ? 8 : total_strengths;

    uint64_t(*mse[2])[TOTAL_STRENGTHS];
    int32_t       pri_damping = 3 + (frm_hdr->quantization_params.base_q_idx >> 6);
    int32_t       i;
    int32_t       nb_strengths;
    int32_t       nb_strength_bits;
    uint64_t      lambda;
    const int32_t num_planes = 3; // av1_num_planes(cm);
    uint16_t qp_index = (uint8_t)pcs_ptr->parent_pcs_ptr->frm_hdr.quantization_params.base_q_idx;
    uint32_t fast_lambda, full_lambda = 0;
    (*av1_lambda_assignment_function_table[pcs_ptr->parent_pcs_ptr->pred_structure])(
        pcs_ptr,
        &fast_lambda,
        &full_lambda,
        (uint8_t)pcs_ptr->parent_pcs_ptr->enhanced_picture_ptr->bit_depth,
        qp_index,
        EB_FALSE);
    lambda = full_lambda;

    mse[0] = (uint64_t(*)[64])malloc(sizeof(**mse) * nvfb * nhfb);
    mse[1] = (uint64_t(*)[64])malloc(sizeof(**mse) * nvfb * nhfb);

    sb_count = 0;
    for (fbr = 0; fbr < nvfb; ++fbr) {
        for (fbc = 0; fbc < nhfb; ++fbc) {
            ModeInfo **mi =
                pcs_ptr->mi_grid_base + MI_SIZE_64X64 * fbr * cm->mi_stride + MI_SIZE_64X64 * fbc;
            const MbModeInfo *mbmi = &mi[0]->mbmi;

            if (((fbc & 1) && (mbmi->block_mi.sb_type == BLOCK_128X128 ||
                               mbmi->block_mi.sb_type == BLOCK_128X64)) ||
                ((fbr & 1) && (mbmi->block_mi.sb_type == BLOCK_128X128 ||
                               mbmi->block_mi.sb_type == BLOCK_64X128))) {
                continue;
            }

            // No filtering if the entire filter block is skipped
            if (svt_sb_all_skip(pcs_ptr, cm, fbr * MI_SIZE_64X64, fbc * MI_SIZE_64X64)) continue;

            for (pli = 0; pli < num_planes; pli++) {
                if (pli == 0)
                    svt_memcpy(mse[0][sb_count],
                               pcs_ptr->mse_seg[0][fbr * nhfb + fbc],
                               TOTAL_STRENGTHS * sizeof(uint64_t));
                if (pli == 2)
                    svt_memcpy(mse[1][sb_count],
                               pcs_ptr->mse_seg[1][fbr * nhfb + fbc],
                               TOTAL_STRENGTHS * sizeof(uint64_t));
                sb_index[sb_count] = MI_SIZE_64X64 * fbr * pcs_ptr->mi_stride + MI_SIZE_64X64 * fbc;
            }
            sb_count++;
        }
    }

    nb_strength_bits = 0;
    /* Search for different number of signalling bits. */
    for (i = 0; i <= 3; i++) {
        int32_t best_lev0[CDEF_MAX_STRENGTHS];
        int32_t best_lev1[CDEF_MAX_STRENGTHS] = {0};
        nb_strengths                          = 1 << i;
        uint64_t tot_mse                      = joint_strength_search_dual(
            best_lev0, best_lev1, nb_strengths, mse, sb_count, start_gi, end_gi);
        (void)joint_strength_search;
        /* Count superblock signalling cost. */
        const int total_bits =
            sb_count * i + nb_strengths * CDEF_STRENGTH_BITS * 2;
        const int      rate_cost = av1_cost_literal(total_bits);
        const uint64_t dist      = tot_mse * 16;
        tot_mse                  = RDCOST(lambda, rate_cost, dist);
        if (tot_mse < best_tot_mse) {
            best_tot_mse     = tot_mse;
            nb_strength_bits = i;
            for (int32_t j = 0; j < 1 << nb_strength_bits; j++) {
                frm_hdr->cdef_params.cdef_y_strength[j]  = best_lev0[j];
                frm_hdr->cdef_params.cdef_uv_strength[j] = best_lev1[j];
            }
        }
    }
    nb_strengths = 1 << nb_strength_bits;

    frm_hdr->cdef_params.cdef_bits = nb_strength_bits;
    ppcs->nb_cdef_strengths        = nb_strengths;
    for (i = 0; i < sb_count; i++) {
        int32_t  gi;
        int32_t  best_gi;
        uint64_t best_mse = (uint64_t)1 << 63;
        best_gi           = 0;
        for (gi = 0; gi < ppcs->nb_cdef_strengths; gi++) {
            uint64_t curr = mse[0][i][frm_hdr->cdef_params.cdef_y_strength[gi]];
            curr += mse[1][i][frm_hdr->cdef_params.cdef_uv_strength[gi]];
            if (curr < best_mse) {
                best_gi  = gi;
                best_mse = curr;
            }
        }
        selected_strength[i] = best_gi;
        selected_strength_cnt[best_gi]++;

        pcs_ptr->mi_grid_base[sb_index[i]]->mbmi.cdef_strength = (int8_t)best_gi;
        //in case the fb is within a block=128x128 or 128x64, or 64x128, then we genrate param only for the first 64x64.
        //since our mi map deos not have the multi pointer single data assignment, we need to duplicate data.
        BlockSize sb_type = pcs_ptr->mi_grid_base[sb_index[i]]->mbmi.block_mi.sb_type;

        switch (sb_type) {
        case BLOCK_128X128:
            pcs_ptr->mi_grid_base[sb_index[i] + MI_SIZE_64X64]->mbmi.cdef_strength =
                (int8_t)best_gi;
            pcs_ptr->mi_grid_base[sb_index[i] + MI_SIZE_64X64 * pcs_ptr->mi_stride]
                ->mbmi.cdef_strength = (int8_t)best_gi;
            pcs_ptr->mi_grid_base[sb_index[i] + MI_SIZE_64X64 * pcs_ptr->mi_stride + MI_SIZE_64X64]
                ->mbmi.cdef_strength = (int8_t)best_gi;
            break;
        case BLOCK_128X64:
            pcs_ptr->mi_grid_base[sb_index[i] + MI_SIZE_64X64]->mbmi.cdef_strength =
                (int8_t)best_gi;
            break;
        case BLOCK_64X128:
            pcs_ptr->mi_grid_base[sb_index[i] + MI_SIZE_64X64 * pcs_ptr->mi_stride]
                ->mbmi.cdef_strength = (int8_t)best_gi;
            break;
        default: break;
        }
    }

    //cdef_pri_damping & cdef_sec_damping consolidated to cdef_damping
    frm_hdr->cdef_params.cdef_damping = pri_damping;
    for (i = 0; i < total_strengths; i++)
        best_frame_gi_cnt += selected_strength_cnt[i] > best_frame_gi_cnt ? 1 : 0;
    ppcs->cdef_frame_strength = ((best_frame_gi_cnt + 4) / 4) * 4;

    free(mse[0]);
    free(mse[1]);
    free(sb_index);
    free(selected_strength);
}
