/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 3-Clause Clear License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 3-Clause Clear License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at https://www.aomedia.org/license. If the Alliance for Open
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
#include "EbRateDistortionCost.h"

void get_recon_pic(PictureControlSet *pcs_ptr, EbPictureBufferDesc **recon_ptr, Bool is_highbd);
static INLINE uint64_t dist_8xn_16bit_c(const uint16_t *src, const uint16_t *dst,
                                        const int32_t dstride, const int32_t coeff_shift,
                                        int8_t height, uint8_t subsampling_factor) {
    uint64_t svar   = 0;
    uint64_t dvar   = 0;
    uint64_t sum_s  = 0;
    uint64_t sum_d  = 0;
    uint64_t sum_s2 = 0;
    uint64_t sum_d2 = 0;
    uint64_t sum_sd = 0;
    int32_t  i, j;
    for (i = 0; i < height; i += subsampling_factor) {
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
    return (uint64_t)floor(.5 +
                           (sum_d2 + sum_s2 - 2 * sum_sd) * .5 *
                               (svar + dvar + (400 << 2 * coeff_shift)) /
                               (sqrt((20000 << 4 * coeff_shift) + svar * (double)dvar)));
}

static INLINE uint64_t mse_8xn_16bit_c(const uint16_t *src, const uint16_t *dst,
                                       const int32_t dstride, const int32_t height,
                                       uint8_t subsampling_factor) {
    uint64_t sum = 0;
    int32_t  i, j;
    for (i = 0; i < height; i += subsampling_factor) {
        for (j = 0; j < 8; j++) {
            int32_t e = dst[i * dstride + j] - src[8 * i + j];
            sum += e * e;
        }
    }
    return sum;
}

static INLINE uint64_t mse_4xn_16bit_c(const uint16_t *src, const uint16_t *dst,
                                       const int32_t dstride, const int32_t height,
                                       uint8_t subsampling_factor) {
    uint64_t sum = 0;
    int32_t  i, j;
    for (i = 0; i < height; i += subsampling_factor) {
        for (j = 0; j < 4; j++) {
            int32_t e = dst[i * dstride + j] - src[4 * i + j];
            sum += e * e;
        }
    }
    return sum;
}

static INLINE uint64_t dist_8xn_8bit_c(const uint8_t *src, const uint8_t *dst,
                                       const int32_t dstride, const int32_t coeff_shift,
                                       uint8_t height, uint8_t subsampling_factor) {
    uint64_t svar   = 0;
    uint64_t dvar   = 0;
    uint64_t sum_s  = 0;
    uint64_t sum_d  = 0;
    uint64_t sum_s2 = 0;
    uint64_t sum_d2 = 0;
    uint64_t sum_sd = 0;
    int32_t  i, j;
    for (i = 0; i < height; i += subsampling_factor) {
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
    return (uint64_t)floor(.5 +
                           (sum_d2 + sum_s2 - 2 * sum_sd) * .5 *
                               (svar + dvar + (400 << 2 * coeff_shift)) /
                               (sqrt((20000 << 4 * coeff_shift) + svar * (double)dvar)));
}
static INLINE uint64_t mse_8xn_8bit_c(const uint8_t *src, const uint8_t *dst, const int32_t dstride,
                                      const int32_t height, uint8_t subsampling_factor) {
    uint64_t sum = 0;
    int32_t  i, j;
    for (i = 0; i < height; i += subsampling_factor) {
        for (j = 0; j < 8; j++) {
            int32_t e = dst[i * dstride + j] - src[8 * i + j];
            sum += e * e;
        }
    }
    return sum;
}

static INLINE uint64_t mse_4xn_8bit_c(const uint8_t *src, const uint8_t *dst, const int32_t dstride,
                                      const int32_t height, uint8_t subsampling_factor) {
    uint64_t sum = 0;
    int32_t  i, j;
    for (i = 0; i < height; i += subsampling_factor) {
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
                             int32_t coeff_shift, int32_t pli, uint8_t subsampling_factor) {
    uint64_t sum = 0;
    int32_t  bi, bx, by;
    if (bsize == BLOCK_8X8) {
        for (bi = 0; bi < cdef_count; bi++) {
            by = dlist[bi].by;
            bx = dlist[bi].bx;
            if (pli == 0) {
                sum += dist_8xn_16bit_c(&src[bi << (3 + 3)],
                                        &dst[(by << 3) * dstride + (bx << 3)],
                                        dstride,
                                        coeff_shift,
                                        8,
                                        subsampling_factor);
            } else
                sum += mse_8xn_16bit_c(&src[bi << (3 + 3)],
                                       &dst[(by << 3) * dstride + (bx << 3)],
                                       dstride,
                                       8,
                                       subsampling_factor);
        }
    } else if (bsize == BLOCK_4X8) {
        for (bi = 0; bi < cdef_count; bi++) {
            by = dlist[bi].by;
            bx = dlist[bi].bx;
            sum += mse_4xn_16bit_c(&src[bi << (3 + 2)],
                                   &dst[(by << 3) * dstride + (bx << 2)],
                                   dstride,
                                   8,
                                   subsampling_factor);
        }
    } else if (bsize == BLOCK_8X4) {
        for (bi = 0; bi < cdef_count; bi++) {
            by = dlist[bi].by;
            bx = dlist[bi].bx;
            sum += mse_8xn_16bit_c(&src[bi << (2 + 3)],
                                   &dst[(by << 2) * dstride + (bx << 3)],
                                   dstride,
                                   4,
                                   subsampling_factor);
        }
    } else {
        assert(bsize == BLOCK_4X4);
        for (bi = 0; bi < cdef_count; bi++) {
            by = dlist[bi].by;
            bx = dlist[bi].bx;
            sum += mse_4xn_16bit_c(&src[bi << (2 + 2)],
                                   &dst[(by << 2) * dstride + (bx << 2)],
                                   dstride,
                                   4,
                                   subsampling_factor);
        }
    }
    return sum >> 2 * coeff_shift;
}

uint64_t compute_cdef_dist_8bit_c(const uint8_t *dst8, int32_t dstride, const uint8_t *src8,
                                  const CdefList *dlist, int32_t cdef_count, BlockSize bsize,
                                  int32_t coeff_shift, int32_t pli, uint8_t subsampling_factor) {
    uint64_t sum = 0;
    int32_t  bi, bx, by;
    if (bsize == BLOCK_8X8) {
        for (bi = 0; bi < cdef_count; bi++) {
            by = dlist[bi].by;
            bx = dlist[bi].bx;
            if (pli == 0) {
                sum += dist_8xn_8bit_c(&src8[bi << (3 + 3)],
                                       &dst8[(by << 3) * dstride + (bx << 3)],
                                       dstride,
                                       coeff_shift,
                                       8,
                                       subsampling_factor);
            } else
                sum += mse_8xn_8bit_c(&src8[bi << (3 + 3)],
                                      &dst8[(by << 3) * dstride + (bx << 3)],
                                      dstride,
                                      8,
                                      subsampling_factor);
        }
    } else if (bsize == BLOCK_4X8) {
        for (bi = 0; bi < cdef_count; bi++) {
            by = dlist[bi].by;
            bx = dlist[bi].bx;
            sum += mse_4xn_8bit_c(&src8[bi << (3 + 2)],
                                  &dst8[(by << 3) * dstride + (bx << 2)],
                                  dstride,
                                  8,
                                  subsampling_factor);
        }
    } else if (bsize == BLOCK_8X4) {
        for (bi = 0; bi < cdef_count; bi++) {
            by = dlist[bi].by;
            bx = dlist[bi].bx;
            sum += mse_8xn_8bit_c(&src8[bi << (2 + 3)],
                                  &dst8[(by << 2) * dstride + (bx << 3)],
                                  dstride,
                                  4,
                                  subsampling_factor);
        }
    } else {
        assert(bsize == BLOCK_4X4);
        for (bi = 0; bi < cdef_count; bi++) {
            by = dlist[bi].by;
            bx = dlist[bi].bx;
            sum += mse_4xn_8bit_c(&src8[bi << (2 + 2)],
                                  &dst8[(by << 2) * dstride + (bx << 2)],
                                  dstride,
                                  4,
                                  subsampling_factor);
        }
    }
    return sum >> 2 * coeff_shift;
}

int32_t svt_sb_all_skip(PictureControlSet *pcs_ptr, const Av1Common *const cm, int32_t mi_row,
                        int32_t mi_col) {
    int32_t maxc, maxr;
    maxc = cm->mi_cols - mi_col;
    maxr = cm->mi_rows - mi_row;

    maxr = AOMMIN(maxr, MI_SIZE_64X64);
    maxc = AOMMIN(maxc, MI_SIZE_64X64);

    for (int32_t r = 0; r < maxr; r++) {
        for (int32_t c = 0; c < maxc; c++) {
            if (!(pcs_ptr->mi_grid_base[(mi_row + r) * pcs_ptr->mi_stride + mi_col + c]
                      ->mbmi.block_mi.skip))
                return 0;
        }
    }
    return 1;
}

int32_t svt_sb_compute_cdef_list(PictureControlSet *pcs_ptr, const Av1Common *const cm,
                                 int32_t mi_row, int32_t mi_col, CdefList *dlist, BlockSize bs) {
    //MbModeInfo **grid = cm->mi_grid_visible;
    ModeInfo **grid      = pcs_ptr->mi_grid_base;
    int32_t    mi_stride = pcs_ptr->mi_stride;

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
            if (!grid[(mi_row + r) * mi_stride + (mi_col + c)]->mbmi.block_mi.skip ||
                !grid[(mi_row + r) * mi_stride + (mi_col + c + 1)]->mbmi.block_mi.skip ||
                !grid[(mi_row + r + 1) * mi_stride + (mi_col + c)]->mbmi.block_mi.skip ||
                !grid[(mi_row + r + 1) * mi_stride + (mi_col + c + 1)]->mbmi.block_mi.skip) {
                dlist[count].by = (uint8_t)(r >> r_shift);
                dlist[count].bx = (uint8_t)(c >> c_shift);
                count++;
            }
        }
    }
    return count;
}

/*
Loop over all 64x64 filter blocks and perform the CDEF filtering for each block, using
the filter strength pairs chosen in finish_cdef_search().
*/
void svt_av1_cdef_frame(SequenceControlSet *scs, PictureControlSet *pcs) {
    struct PictureParentControlSet *ppcs     = pcs->parent_pcs_ptr;
    Av1Common                      *cm       = ppcs->av1_cm;
    FrameHeader                    *frm_hdr  = &ppcs->frm_hdr;
    Bool                            is_16bit = scs->is_16bit_pipeline;

    EbPictureBufferDesc *recon_pic;
    get_recon_pic(pcs, &recon_pic, is_16bit);

    const uint32_t offset_y       = recon_pic->origin_x + recon_pic->origin_y * recon_pic->stride_y;
    EbByte         recon_buffer_y = recon_pic->buffer_y + (offset_y << is_16bit);
    const uint32_t offset_cb = (recon_pic->origin_x + recon_pic->origin_y * recon_pic->stride_cb) >>
        1;
    EbByte         recon_buffer_cb = recon_pic->buffer_cb + (offset_cb << is_16bit);
    const uint32_t offset_cr = (recon_pic->origin_x + recon_pic->origin_y * recon_pic->stride_cr) >>
        1;
    EbByte recon_buffer_cr = recon_pic->buffer_cr + (offset_cr << is_16bit);

    const int32_t num_planes = av1_num_planes(&scs->seq_header.color_config);
    DECLARE_ALIGNED(16, uint16_t, src[CDEF_INBUF_SIZE]);
    uint16_t      *linebuf[3];
    uint16_t      *colbuf[3];
    CdefList       dlist[MI_SIZE_64X64 * MI_SIZE_64X64];
    uint8_t       *row_cdef, *prev_row_cdef, *curr_row_cdef;
    int32_t        cdef_count;
    const uint32_t sb_size = scs->super_block_size;
    int32_t        mi_wide_l2[3];
    int32_t        mi_high_l2[3];
    int32_t        xdec[3];
    int32_t        ydec[3];
    int32_t        coeff_shift = AOMMAX(scs->static_config.encoder_bit_depth - 8, 0);
    const int32_t  nvfb        = (cm->mi_rows + MI_SIZE_64X64 - 1) / MI_SIZE_64X64;
    const int32_t  nhfb        = (cm->mi_cols + MI_SIZE_64X64 - 1) / MI_SIZE_64X64;

    row_cdef = (uint8_t *)svt_aom_malloc(sizeof(*row_cdef) * (nhfb + 2) * 2);
    assert(row_cdef != NULL);
    memset(row_cdef, 1, sizeof(*row_cdef) * (nhfb + 2) * 2);
    prev_row_cdef = row_cdef + 1;
    curr_row_cdef = prev_row_cdef + nhfb + 2;
    for (int32_t pli = 0; pli < num_planes; pli++) {
        int32_t subsampling_x = (pli == 0) ? 0 : 1;
        int32_t subsampling_y = (pli == 0) ? 0 : 1;
        xdec[pli]             = subsampling_x; //CHKN xd->plane[pli].subsampling_x;
        ydec[pli]             = subsampling_y; //CHKN  xd->plane[pli].subsampling_y;
        mi_wide_l2[pli]       = MI_SIZE_LOG2 - subsampling_x; //CHKN xd->plane[pli].subsampling_x;
        mi_high_l2[pli]       = MI_SIZE_LOG2 - subsampling_y; //CHKN xd->plane[pli].subsampling_y;
    }

    const int32_t stride = (cm->mi_cols << MI_SIZE_LOG2) + 2 * CDEF_HBORDER;
    for (int32_t pli = 0; pli < num_planes; pli++) {
        linebuf[pli] = (uint16_t *)svt_aom_malloc(sizeof(*linebuf) * CDEF_VBORDER * stride);
        colbuf[pli]  = (uint16_t *)svt_aom_malloc(
            sizeof(*colbuf) * ((CDEF_BLOCKSIZE << mi_high_l2[pli]) + 2 * CDEF_VBORDER) *
            CDEF_HBORDER);
    }

    for (int32_t fbr = 0; fbr < nvfb; fbr++) {
        int32_t cdef_left = 1;
        for (int32_t fbc = 0; fbc < nhfb; fbc++) {
            int32_t level, sec_strength;
            int32_t uv_level, uv_sec_strength;
            int32_t nhb, nvb;
            int32_t cstart     = 0;
            curr_row_cdef[fbc] = 0;
            assert(pcs->mi_grid_base[MI_SIZE_64X64 * fbr * cm->mi_stride + MI_SIZE_64X64 * fbc] !=
                       NULL &&
                   "CDEF ERROR: Skipping Current FB");
            assert(pcs->mi_grid_base[MI_SIZE_64X64 * fbr * cm->mi_stride + MI_SIZE_64X64 * fbc]
                           ->mbmi.cdef_strength != -1 &&
                   "CDEF ERROR: Skipping Current FB");
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

            // Find the index of the CDEF strength for the filter block
            const int32_t mbmi_cdef_strength =
                pcs->mi_grid_base[MI_SIZE_64X64 * fbr * cm->mi_stride + MI_SIZE_64X64 * fbc]
                    ->mbmi.cdef_strength;
            level = frm_hdr->cdef_params.cdef_y_strength[mbmi_cdef_strength] / CDEF_SEC_STRENGTHS;
            sec_strength = frm_hdr->cdef_params.cdef_y_strength[mbmi_cdef_strength] %
                CDEF_SEC_STRENGTHS;
            // Secondary luma strength takes values in {0, 1, 2, 4}. If sec_strength is equal to 3 from the step above, change it to 4.
            sec_strength += sec_strength == 3;
            // Set primary and secondary chroma strengths.
            uv_level = frm_hdr->cdef_params.cdef_uv_strength[mbmi_cdef_strength] /
                CDEF_SEC_STRENGTHS;
            uv_sec_strength = frm_hdr->cdef_params.cdef_uv_strength[mbmi_cdef_strength] %
                CDEF_SEC_STRENGTHS;
            // Secondary chroma strength takes values in {0, 1, 2, 4}. If sec_strength is equal to 3 from the step above, change it to 4.
            uv_sec_strength += uv_sec_strength == 3;
            if ((level == 0 && sec_strength == 0 && uv_level == 0 && uv_sec_strength == 0) ||
                (cdef_count = svt_sb_compute_cdef_list(
                     pcs, cm, fbr * MI_SIZE_64X64, fbc * MI_SIZE_64X64, dlist, BLOCK_64X64)) == 0) {
                cdef_left = 0;
                continue;
            }

            int dirinit = !(ppcs->cdef_ctrls.use_reference_cdef_fs);
            // When SB 128 is used, the search for certain blocks is skipped, so dir/var info is not generated
            // In those cases, must generate info here
            if (sb_size == 128) {
                const uint32_t    lc      = MI_SIZE_64X64 * fbc;
                const uint32_t    lr      = MI_SIZE_64X64 * fbr;
                ModeInfo        **mi      = pcs->mi_grid_base + lr * cm->mi_stride + lc;
                const MbModeInfo *mbmi    = &mi[0]->mbmi;
                const BlockSize   sb_type = mbmi->block_mi.sb_type;
                if (((fbc & 1) && (sb_type == BLOCK_128X128 || sb_type == BLOCK_128X64)) ||
                    ((fbr & 1) && (sb_type == BLOCK_128X128 || sb_type == BLOCK_64X128)))
                    dirinit = 0;
            }
            uint8_t(*dir)[CDEF_NBLOCKS][CDEF_NBLOCKS] = &pcs->cdef_dir_data[fbr * nhfb + fbc].dir;
            int32_t(*var)[CDEF_NBLOCKS][CDEF_NBLOCKS] = &pcs->cdef_dir_data[fbr * nhfb + fbc].var;
            curr_row_cdef[fbc]                        = 1;
            for (int32_t pli = 0; pli < num_planes; pli++) {
                int32_t coffset;
                int32_t rend, cend;
                int32_t pri_damping = frm_hdr->cdef_params.cdef_damping;
                int32_t sec_damping = pri_damping;
                int32_t hsize       = nhb << mi_wide_l2[pli];
                int32_t vsize       = nvb << mi_high_l2[pli];
                if (fbc == nhfb - 1)
                    cend = hsize;
                else
                    cend = hsize + CDEF_HBORDER;

                if (fbr == nvfb - 1)
                    rend = vsize;
                else
                    rend = vsize + CDEF_VBORDER;

                coffset             = fbc * MI_SIZE_64X64 << mi_wide_l2[pli];
                EbByte   rec_buff   = 0;
                uint32_t rec_stride = 0;

                switch (pli) {
                case 0:
                    rec_buff   = recon_buffer_y;
                    rec_stride = recon_pic->stride_y;
                    break;
                case 1:
                    rec_buff     = recon_buffer_cb;
                    rec_stride   = recon_pic->stride_cb;
                    level        = uv_level;
                    sec_strength = uv_sec_strength;
                    break;
                case 2:
                    rec_buff     = recon_buffer_cr;
                    rec_stride   = recon_pic->stride_cr;
                    level        = uv_level;
                    sec_strength = uv_sec_strength;
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
                            cend - cstart,
                            is_16bit);
                if (!prev_row_cdef[fbc]) {
                    copy_sb8_16(&src[CDEF_HBORDER],
                                CDEF_BSTRIDE,
                                rec_buff,
                                (MI_SIZE_64X64 << mi_high_l2[pli]) * fbr - CDEF_VBORDER,
                                coffset,
                                rec_stride,
                                CDEF_VBORDER,
                                hsize,
                                is_16bit);
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
                                CDEF_HBORDER,
                                is_16bit);
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
                                CDEF_HBORDER,
                                is_16bit);
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
                                hsize,
                                is_16bit);

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
                // if ppcs->cdef_ctrls.use_reference_cdef_fs is true, then search was not performed
                // Therefore, need to make sure dir and var are initialized
                if (level || sec_strength || !dirinit) {
                    svt_cdef_filter_fb(
                        is_16bit ? NULL
                                 : &rec_buff[rec_stride * (MI_SIZE_64X64 * fbr << mi_high_l2[pli]) +
                                             (fbc * MI_SIZE_64X64 << mi_wide_l2[pli])],
                        is_16bit
                            ? &((uint16_t *)rec_buff)[rec_stride *
                                                          (MI_SIZE_64X64 * fbr << mi_high_l2[pli]) +
                                                      (fbc * MI_SIZE_64X64 << mi_wide_l2[pli])]
                            : NULL,
                        rec_stride,
                        &src[CDEF_VBORDER * CDEF_BSTRIDE + CDEF_HBORDER],
                        xdec[pli],
                        ydec[pli],
                        *dir,
                        &dirinit,
                        *var,
                        pli,
                        dlist,
                        cdef_count,
                        level,
                        sec_strength,
                        pri_damping,
                        sec_damping,
                        coeff_shift,
                        1); // no subsampling
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

///-------search
/*
 * Search for the best luma+chroma strength to add as an option, knowing we
 * already selected nb_strengths options
 *
 * Params:
 *
 * lev0 : Array of indices of selected luma strengths.
 * lev1 : Array of indices of selected chroma strengths.
 * nb_strengths : Number of selected (Luma_strength, Chroma_strength) pairs.
 * mse : Array of luma and chroma filtering mse values.
 * sb_count : Number of filter blocks in the frame.
 * start_gi : starting strength index for the search of the additional strengths.
 * end_gi : End index for the for the search of the additional strengths.
*/
uint64_t svt_search_one_dual_c(int *lev0, int *lev1, int nb_strengths, uint64_t **mse[2],
                               int sb_count, int start_gi, int end_gi) {
    uint64_t      tot_mse[TOTAL_STRENGTHS][TOTAL_STRENGTHS];
    int32_t       i, j;
    uint64_t      best_tot_mse    = (uint64_t)1 << 63;
    int32_t       best_id0        = 0;
    int32_t       best_id1        = 0;
    const int32_t total_strengths = end_gi;
    memset(tot_mse, 0, sizeof(tot_mse));
    /* Loop over the filter blocks in the frame */
    for (i = 0; i < sb_count; i++) {
        int32_t  gi;
        uint64_t best_mse = (uint64_t)1 << 63;
        /* Loop over the already selected nb_strengths (Luma_strength,
           Chroma_strength) pairs, and find the pair that has the smallest mse
           (best_mse) for the current filter block.*/
        /* Find best mse among already selected options. */
        for (gi = 0; gi < nb_strengths; gi++) {
            uint64_t curr = mse[0][i][lev0[gi]];
            curr += mse[1][i][lev1[gi]];
            if (curr < best_mse)
                best_mse = curr;
        }
        /* Loop over the set of available (Luma_strength, Chroma_strength)
           pairs, identify any that provide an mse better than best_mse from the
           step above for the current filter block, and update any corresponding
           total mse (tot_mse[j][k]). */
        /* Find best mse when adding each possible new option. */
        for (j = start_gi; j < total_strengths; j++) {
            int32_t k;
            for (k = start_gi; k < total_strengths; k++) {
                uint64_t best = best_mse;
                uint64_t curr = mse[0][i][j];
                curr += mse[1][i][k];
                if (curr < best)
                    best = curr;
                tot_mse[j][k] += best;
            }
        }
    }
    /* Loop over the additionally searched (Luma_strength, Chroma_strength) pairs
       from the step above, and identify any such pair that provided the best mse for
       the whole frame. The identified pair would be added to the set of already selected pairs. */
    for (j = start_gi; j < total_strengths;
         j++) { // Loop over the additionally searched luma strengths
        int32_t k;
        for (k = start_gi; k < total_strengths;
             k++) { // Loop over the additionally searched chroma strengths
            if (tot_mse[j][k] < best_tot_mse) {
                best_tot_mse = tot_mse[j][k];
                best_id0     = j; // index for the best luma strength
                best_id1     = k; // index for the best chroma strength
            }
        }
    }
    lev0[nb_strengths] =
        best_id0; // Add the identified luma strength to the list of selected luma strengths
    lev1[nb_strengths] =
        best_id1; // Add the identified chroma strength to the list of selected chroma strengths
    return best_tot_mse;
}
/*
 * Search for the set of luma+chroma strengths that minimizes mse.
 *
 * Params:
 *
 * best_lev0 : Array of indices of selected luma strengths.
 * best_lev1 : Array of indices of selected chroma strengths.
 * nb_strengths : Number of selected (Luma_strength, Chroma_strength) pairs.
 * mse : Array of luma and chroma filtering mse values.
 * sb_count : Number of filter blocks in the frame.
 * start_gi : starting strength index for the search of the additional strengths.
 * end_gi : End index for the for the search of the additional strengths.
*/
static uint64_t joint_strength_search_dual(int32_t *best_lev0, int32_t *best_lev1,
                                           int32_t nb_strengths, uint64_t **mse[2],
                                           int32_t sb_count, int32_t start_gi, int32_t end_gi) {
    uint64_t best_tot_mse;
    int32_t  i;
    best_tot_mse = (uint64_t)1 << 63;
    /* Greedy search: add one strength options at a time.

    Determine nb_strengths (Luma_strength, Chroma_strength) pairs.
    The list of nb_strengths pairs is determined by adding one such pair at
    a time through the call to the function search_one_dual. When the
    function search_one_dual is called, the search accounts for the
    strength pairs that have already been added in the previous iteration of
    the loop below. The loop below returns in the end best_tot_mse
    representing the best filtering mse for the whole frame based on the
    selected list of best (Luma_strength, Chroma_strength) pairs.
    */
    for (i = 0; i < nb_strengths; i++)
        best_tot_mse = svt_search_one_dual(
            best_lev0, best_lev1, i, mse, sb_count, start_gi, end_gi);
    /* Performing further refinements on the search based on the results
    from the step above. Trying to refine the greedy search by reconsidering each
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
void finish_cdef_search(PictureControlSet *pcs_ptr) {
    struct PictureParentControlSet *ppcs    = pcs_ptr->parent_pcs_ptr;
    FrameHeader                    *frm_hdr = &ppcs->frm_hdr;
    Av1Common                      *cm      = ppcs->av1_cm;
    int32_t                         mi_rows = ppcs->av1_cm->mi_rows;
    int32_t                         mi_cols = ppcs->av1_cm->mi_cols;

    int32_t  fbr, fbc;
    uint64_t best_tot_mse = (uint64_t)1 << 63;
    int32_t  sb_count;
    int32_t  nvfb = (mi_rows + MI_SIZE_64X64 - 1) / MI_SIZE_64X64;
    int32_t  nhfb = (mi_cols + MI_SIZE_64X64 - 1) / MI_SIZE_64X64;
    //CDEF Settings
    CdefControls *cdef_ctrls                 = &pcs_ptr->parent_pcs_ptr->cdef_ctrls;
    const int     first_pass_fs_num          = cdef_ctrls->first_pass_fs_num;
    const int     default_second_pass_fs_num = cdef_ctrls->default_second_pass_fs_num;
    if (cdef_ctrls->use_reference_cdef_fs) {
        int32_t *sb_index = (int32_t *)malloc(nvfb * nhfb * sizeof(*sb_index));
        int32_t  best_gi  = 0;
        sb_count          = 0;
        assert(sb_index != NULL);
        for (fbr = 0; fbr < nvfb; ++fbr) {
            for (fbc = 0; fbc < nhfb; ++fbc) {
                ModeInfo **mi = pcs_ptr->mi_grid_base + MI_SIZE_64X64 * fbr * cm->mi_stride +
                    MI_SIZE_64X64 * fbc;
                const MbModeInfo *mbmi = &mi[0]->mbmi;

                if (((fbc & 1) &&
                     (mbmi->block_mi.sb_type == BLOCK_128X128 ||
                      mbmi->block_mi.sb_type == BLOCK_128X64)) ||
                    ((fbr & 1) &&
                     (mbmi->block_mi.sb_type == BLOCK_128X128 ||
                      mbmi->block_mi.sb_type == BLOCK_64X128))) {
                    continue;
                }
                // No filtering if the entire filter block is skipped
                if (svt_sb_all_skip(pcs_ptr, cm, fbr * MI_SIZE_64X64, fbc * MI_SIZE_64X64))
                    continue;
                sb_index[sb_count] = MI_SIZE_64X64 * fbr * pcs_ptr->mi_stride + MI_SIZE_64X64 * fbc;
                sb_count++;
            }
        }
        for (int32_t i = 0; i < sb_count; i++) {
            pcs_ptr->mi_grid_base[sb_index[i]]->mbmi.cdef_strength = (int8_t)best_gi;
            //in case the fb is within a block=128x128 or 128x64, or 64x128, then we genrate param only for the first 64x64.
            //since our mi map deos not have the multi pointer single data assignment, we need to duplicate data.
            BlockSize sb_type = pcs_ptr->mi_grid_base[sb_index[i]]->mbmi.block_mi.sb_type;
            switch (sb_type) {
            case BLOCK_128X128:
                pcs_ptr->mi_grid_base[sb_index[i] + MI_SIZE_64X64]->mbmi.cdef_strength = (int8_t)
                    best_gi;
                pcs_ptr->mi_grid_base[sb_index[i] + MI_SIZE_64X64 * pcs_ptr->mi_stride]
                    ->mbmi.cdef_strength = (int8_t)best_gi;
                pcs_ptr
                    ->mi_grid_base[sb_index[i] + MI_SIZE_64X64 * pcs_ptr->mi_stride + MI_SIZE_64X64]
                    ->mbmi.cdef_strength = (int8_t)best_gi;
                break;
            case BLOCK_128X64:
                pcs_ptr->mi_grid_base[sb_index[i] + MI_SIZE_64X64]->mbmi.cdef_strength = (int8_t)
                    best_gi;
                break;
            case BLOCK_64X128:
                pcs_ptr->mi_grid_base[sb_index[i] + MI_SIZE_64X64 * pcs_ptr->mi_stride]
                    ->mbmi.cdef_strength = (int8_t)best_gi;
                break;
            default: break;
            }
        }
        frm_hdr->cdef_params.cdef_bits = 0;
        ppcs->nb_cdef_strengths        = 1;
        //cdef_pri_damping & cdef_sec_damping consolidated to cdef_damping
        int32_t pri_damping               = 3 + (frm_hdr->quantization_params.base_q_idx >> 6);
        frm_hdr->cdef_params.cdef_damping = pri_damping;
        frm_hdr->cdef_params.cdef_y_strength[0]  = cdef_ctrls->pred_y_f;
        frm_hdr->cdef_params.cdef_uv_strength[0] = cdef_ctrls->pred_uv_f;
        free(sb_index);
        return;
    }
    int32_t *sb_index = (int32_t *)malloc(nvfb * nhfb * sizeof(*sb_index));
    // to keep track of the sb_address in units of SBs (not mi_size)
    int32_t *sb_addr  = (int32_t *)malloc(nvfb * nhfb * sizeof(*sb_index));
    int32_t  start_gi = 0;
    int32_t  end_gi   = first_pass_fs_num + default_second_pass_fs_num;
    assert(sb_index != NULL);
    uint64_t **mse[2];
    int32_t    i;
    int32_t    nb_strengths;
    int32_t    nb_strength_bits;
    uint64_t   lambda;
    uint32_t   fast_lambda, full_lambda = 0;
    (*av1_lambda_assignment_function_table[pcs_ptr->parent_pcs_ptr->pred_structure])(
        pcs_ptr,
        &fast_lambda,
        &full_lambda,
        (uint8_t)pcs_ptr->parent_pcs_ptr->enhanced_picture_ptr->bit_depth,
        pcs_ptr->parent_pcs_ptr->frm_hdr.quantization_params.base_q_idx,
        FALSE);
    lambda   = full_lambda;
    mse[0]   = (uint64_t **)malloc(sizeof(*mse) * nvfb * nhfb);
    mse[1]   = (uint64_t **)malloc(sizeof(*mse) * nvfb * nhfb);
    sb_count = 0;
    for (fbr = 0; fbr < nvfb; ++fbr) {
        for (fbc = 0; fbc < nhfb; ++fbc) {
            ModeInfo **mi = pcs_ptr->mi_grid_base + MI_SIZE_64X64 * fbr * cm->mi_stride +
                MI_SIZE_64X64 * fbc;
            const MbModeInfo *mbmi = &mi[0]->mbmi;

            if (((fbc & 1) &&
                 (mbmi->block_mi.sb_type == BLOCK_128X128 ||
                  mbmi->block_mi.sb_type == BLOCK_128X64)) ||
                ((fbr & 1) &&
                 (mbmi->block_mi.sb_type == BLOCK_128X128 ||
                  mbmi->block_mi.sb_type == BLOCK_64X128))) {
                continue;
            }

            // No filtering if the entire filter block is skipped
            if (pcs_ptr->skip_cdef_seg[fbr * nhfb + fbc])
                continue;
            // point to the MSE data
            mse[0][sb_count] = pcs_ptr->mse_seg[0][fbr * nhfb + fbc];
            mse[1][sb_count] = pcs_ptr->mse_seg[1][fbr * nhfb + fbc];

            sb_index[sb_count] = MI_SIZE_64X64 * fbr * pcs_ptr->mi_stride + MI_SIZE_64X64 * fbc;
            sb_addr[sb_count]  = fbr * nhfb + fbc;
            sb_count++;
        }
    }

    nb_strength_bits = 0;
    // Scale down the cost of the (0,0) filter strength to bias selection towards off.
    // When off, can save the cost of the application.
    if (cdef_ctrls->zero_fs_cost_bias) {
        const uint16_t factor = cdef_ctrls->zero_fs_cost_bias;
        for (i = 0; i < sb_count; i++) {
            mse[0][i][0] = (factor * mse[0][i][0]) >> 6;
            mse[1][i][0] = (factor * mse[1][i][0]) >> 6;
        }
    }
    /* Search for different number of signalling bits. */
    for (i = 0; i <= 3; i++) {
        int32_t best_lev0[CDEF_MAX_STRENGTHS] = {0};
        int32_t best_lev1[CDEF_MAX_STRENGTHS] = {0};
        nb_strengths                          = 1 << i;
        uint64_t tot_mse                      = joint_strength_search_dual(
            best_lev0, best_lev1, nb_strengths, mse, sb_count, start_gi, end_gi);
        /* Count superblock signalling cost. */
        const int      total_bits = sb_count * i + nb_strengths * CDEF_STRENGTH_BITS * 2;
        const int      rate_cost  = av1_cost_literal(total_bits);
        const uint64_t dist       = tot_mse * 16;
        tot_mse                   = RDCOST(lambda, rate_cost, dist);
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
        // skip this loop for SBs that are skipped in the search
        for (gi = 0; gi < ppcs->nb_cdef_strengths; gi++) {
            uint64_t curr = mse[0][i][frm_hdr->cdef_params.cdef_y_strength[gi]];
            curr += mse[1][i][frm_hdr->cdef_params.cdef_uv_strength[gi]];
            if (curr < best_mse) {
                best_gi  = gi;
                best_mse = curr;
            }
        }

        pcs_ptr->mi_grid_base[sb_index[i]]->mbmi.cdef_strength = (int8_t)best_gi;
        //in case the fb is within a block=128x128 or 128x64, or 64x128, then we genrate param only for the first 64x64.
        //since our mi map deos not have the multi pointer single data assignment, we need to duplicate data.
        BlockSize sb_type = pcs_ptr->mi_grid_base[sb_index[i]]->mbmi.block_mi.sb_type;

        switch (sb_type) {
        case BLOCK_128X128:
            pcs_ptr->mi_grid_base[sb_index[i] + MI_SIZE_64X64]->mbmi.cdef_strength = (int8_t)
                best_gi;
            pcs_ptr->mi_grid_base[sb_index[i] + MI_SIZE_64X64 * pcs_ptr->mi_stride]
                ->mbmi.cdef_strength = (int8_t)best_gi;
            pcs_ptr->mi_grid_base[sb_index[i] + MI_SIZE_64X64 * pcs_ptr->mi_stride + MI_SIZE_64X64]
                ->mbmi.cdef_strength = (int8_t)best_gi;
            break;
        case BLOCK_128X64:
            pcs_ptr->mi_grid_base[sb_index[i] + MI_SIZE_64X64]->mbmi.cdef_strength = (int8_t)
                best_gi;
            break;
        case BLOCK_64X128:
            pcs_ptr->mi_grid_base[sb_index[i] + MI_SIZE_64X64 * pcs_ptr->mi_stride]
                ->mbmi.cdef_strength = (int8_t)best_gi;
            break;
        default: break;
        }
    }
    int filter_map[TOTAL_STRENGTHS] = {0};
    for (i = 0; i < first_pass_fs_num; i++) filter_map[i] = cdef_ctrls->default_first_pass_fs[i];
    for (i = first_pass_fs_num; i < (first_pass_fs_num + default_second_pass_fs_num); i++)
        filter_map[i] = cdef_ctrls->default_second_pass_fs[i - first_pass_fs_num];

    for (i = 0; i < ppcs->nb_cdef_strengths; i++) {
        frm_hdr->cdef_params.cdef_y_strength[i] =
            filter_map[frm_hdr->cdef_params.cdef_y_strength[i]];
        frm_hdr->cdef_params.cdef_uv_strength[i] =
            filter_map[frm_hdr->cdef_params.cdef_uv_strength[i]];
    }
    //cdef_pri_damping & cdef_sec_damping consolidated to cdef_damping
    frm_hdr->cdef_params.cdef_damping = 3 + (frm_hdr->quantization_params.base_q_idx >> 6);
    free(mse[0]);
    free(mse[1]);
    free(sb_index);
    free(sb_addr);
}
