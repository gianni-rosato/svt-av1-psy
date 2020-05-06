/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

// SUMMARY
//   Contains the Decode Block related functions

/**************************************
 * Includes
 **************************************/

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "EbDefinitions.h"

#include "EbSvtAv1Dec.h"
#include "EbDecHandle.h"

#include "EbDecParseHelper.h"
#include "EbCommonUtils.h"

#include "EbDecInverseQuantize.h"
#include "EbDecProcessFrame.h"

#include "EbDecIntraPrediction.h"
#include "EbDecInterPrediction.h"

#include "EbDecNbr.h"
#include "EbDecUtils.h"

#include "EbInvTransforms.h"

#include "EbDecLF.h"
#include "EbDecPicMgr.h"

#include "EbWarpedMotion.h"

extern int select_samples(MV *mv, int *pts, int *pts_inref, int len, BlockSize bsize);

CflAllowedType store_cfl_required(const EbColorConfig *cc, PartitionInfo *xd,
                                  int32_t is_chroma_ref) {
    const BlockModeInfo *mbmi = xd->mi;

    if (cc->mono_chrome) return CFL_DISALLOWED;
    if (!is_chroma_ref) {
        // For non-chroma-reference blocks, we should always store the luma pixels,
        // in case the corresponding chroma-reference block uses CfL.
        // Note that this can only happen for block sizes which are <8 on
        // their shortest side, as otherwise they would be chroma reference
        // blocks.
        return CFL_ALLOWED;
    }

    // If this block has chroma information, we know whether we're
    // actually going to perform a CfL prediction
    return (CflAllowedType)(!is_inter_block(mbmi) && mbmi->uv_mode == UV_CFL_PRED);
}

/* decode partition */
PartitionType get_partition(DecModCtxt *dec_mod_ctxt, FrameHeader *frame_header, uint32_t mi_row,
                            uint32_t mi_col, SBInfo *sb_info, BlockSize bsize) {
    if (mi_row >= frame_header->mi_rows || mi_col >= frame_header->mi_cols)
        return PARTITION_INVALID;

    BlockModeInfo *mode_info =
        get_cur_mode_info(dec_mod_ctxt->dec_handle_ptr, mi_row, mi_col, sb_info);

    const BlockSize subsize = mode_info->sb_type;

    if (subsize == bsize) return PARTITION_NONE;

    const int bhigh  = mi_size_high[bsize];
    const int bwide  = mi_size_wide[bsize];
    const int sshigh = mi_size_high[subsize];
    const int sswide = mi_size_wide[subsize];

    if (bsize > BLOCK_8X8 && mi_row + bwide / 2 < frame_header->mi_rows &&
        mi_col + bhigh / 2 < frame_header->mi_cols) {
        // In this case, the block might be using an extended partition type.
        /* TODO: Fix the nbr access! */
        const BlockModeInfo *const mbmi_right =
            get_cur_mode_info(dec_mod_ctxt->dec_handle_ptr, mi_row, mi_col + (bwide / 2), sb_info);
        const BlockModeInfo *const mbmi_below =
            get_cur_mode_info(dec_mod_ctxt->dec_handle_ptr, mi_row + (bhigh / 2), mi_col, sb_info);

        if (sswide == bwide) {
            // Smaller height but same width. Is PARTITION_HORZ_4, PARTITION_HORZ or
            // PARTITION_HORZ_B. To distinguish the latter two, check if the lower
            // half was split.
            if (sshigh * 4 == bhigh) return PARTITION_HORZ_4;
            assert(sshigh * 2 == bhigh);

            if (mbmi_below->sb_type == subsize)
                return PARTITION_HORZ;
            else
                return PARTITION_HORZ_B;
        } else if (sshigh == bhigh) {
            // Smaller width but same height. Is PARTITION_VERT_4, PARTITION_VERT or
            // PARTITION_VERT_B. To distinguish the latter two, check if the right
            // half was split.
            if (sswide * 4 == bwide) return PARTITION_VERT_4;
            assert(sswide * 2 == bhigh);

            if (mbmi_right->sb_type == subsize)
                return PARTITION_VERT;
            else
                return PARTITION_VERT_B;
        } else {
            // Smaller width and smaller height. Might be PARTITION_SPLIT or could be
            // PARTITION_HORZ_A or PARTITION_VERT_A. If subsize isn't halved in both
            // dimensions, we immediately know this is a split (which will recurse to
            // get to subsize). Otherwise look down and to the right. With
            // PARTITION_VERT_A, the right block will have height bhigh; with
            // PARTITION_HORZ_A, the lower block with have width bwide. Otherwise
            // it's PARTITION_SPLIT.
            if (sswide * 2 != bwide || sshigh * 2 != bhigh) return PARTITION_SPLIT;

            if (mi_size_wide[mbmi_below->sb_type] == bwide) return PARTITION_HORZ_A;
            if (mi_size_high[mbmi_right->sb_type] == bhigh) return PARTITION_VERT_A;

            return PARTITION_SPLIT;
        }
    }

    const int vert_split = sswide < bwide;
    const int horz_split = sshigh < bhigh;
    const int split_idx  = (vert_split << 1) | horz_split;
    assert(split_idx != 0);

    static const PartitionType base_partitions[4] = {
        PARTITION_INVALID, PARTITION_HORZ, PARTITION_VERT, PARTITION_SPLIT};

    return base_partitions[split_idx];
}

void decode_block(DecModCtxt *dec_mod_ctxt, BlockModeInfo *mode_info, int32_t mi_row, int32_t mi_col,
                   BlockSize bsize, TileInfo *tile, SBInfo *sb_info) {
    EbDecHandle *        dec_handle        = (EbDecHandle *)dec_mod_ctxt->dec_handle_ptr;
    EbColorConfig *      color_config      = &dec_mod_ctxt->seq_header->color_config;
    EbPictureBufferDesc *recon_picture_buf = dec_handle->cur_pic_buf[0]->ps_pic_buf;
    uint32_t             mi_cols           = dec_mod_ctxt->frame_header->mi_cols;
    uint32_t             mi_rows           = dec_mod_ctxt->frame_header->mi_rows;

    int num_planes = av1_num_planes(color_config);

    bool inter_block = is_inter_block(mode_info);

    EbBool is16b = dec_handle->is_16bit_pipeline;
#if MODE_INFO_DBG
    assert(mode_info->mi_row == mi_row);
    assert(mode_info->mi_col == mi_col);
#endif
    int32_t bw4 = mi_size_wide[bsize];
    int32_t bh4 = mi_size_high[bsize];

    int hbd = (recon_picture_buf->bit_depth > EB_8BIT) || is16b;

    int sub_x, sub_y, n_coeffs;
    sub_x                 = color_config->subsampling_x;
    sub_y                 = color_config->subsampling_y;
    int32_t is_chroma_ref = is_chroma_reference(mi_row, mi_col, bsize, sub_x, sub_y);

    /* TODO: Can move to a common init fun for parse & decode */
    PartitionInfo part_info;
    part_info.mi          = mode_info;
    part_info.sb_info     = sb_info;
    part_info.mi_row      = mi_row;
    part_info.mi_col      = mi_col;
    part_info.pv_cfl_ctxt = &dec_mod_ctxt->cfl_ctx;
    part_info.mc_buf[0] = dec_mod_ctxt->mc_buf[0];
    part_info.mc_buf[1] = dec_mod_ctxt->mc_buf[1];
    part_info.is_chroma_ref = is_chroma_ref;

    /*!< Distance of MB away from frame edges in subpixels (1/8th pixel). */
    part_info.mb_to_top_edge    = -((mi_row * MI_SIZE) * 8);
    part_info.mb_to_bottom_edge = ((mi_rows - bh4 - mi_row) * MI_SIZE) * 8;
    part_info.mb_to_left_edge   = -((mi_col * MI_SIZE) * 8);
    part_info.mb_to_right_edge  = ((mi_cols - bw4 - mi_col) * MI_SIZE) * 8;

    /*!< Block Size width & height in pixels. */
    /* For Luma bock */
    part_info.wpx[0] = bw4 * MI_SIZE;
    part_info.hpx[0] = bh4 * MI_SIZE;

    /* For U plane chroma bock */
    part_info.wpx[1] = (AOMMAX(1, bw4 >> sub_x)) * MI_SIZE;
    part_info.hpx[1] = (AOMMAX(1, bh4 >> sub_y)) * MI_SIZE;

    /* For V plane chroma bock */
    part_info.wpx[2] = (AOMMAX(1, bw4 >> sub_x)) * MI_SIZE;
    part_info.hpx[2] = (AOMMAX(1, bh4 >> sub_y)) * MI_SIZE;

    /* TODO : tile->tile_rows boundary condn check is wrong */
    part_info.up_available          = (mi_row > tile->mi_row_start);
    part_info.left_available        = (mi_col > tile->mi_col_start);
    part_info.chroma_up_available   = part_info.up_available;
    part_info.chroma_left_available = part_info.left_available;
    if (is_chroma_ref) {
        if (bh4 == 1 && sub_y) part_info.chroma_up_available = (mi_row - 1) > tile->mi_row_start;
        if (bw4 == 1 && sub_x) part_info.chroma_left_available = (mi_col - 1) > tile->mi_col_start;
    } else {
        part_info.chroma_up_available   = 0;
        part_info.chroma_left_available = 0;
    }

    part_info.is_sec_rect = 0;
    if (bw4 < bh4) {
        if (!((mi_col + bw4) & (bh4 - 1))) part_info.is_sec_rect = 1;
    }

    if (bw4 > bh4)
        if (mi_row & (bw4 - 1)) part_info.is_sec_rect = 1;

    if (part_info.up_available)
        part_info.above_mbmi = get_top_mode_info(dec_handle, mi_row, mi_col, sb_info);
    else
        part_info.above_mbmi = NULL;
    if (part_info.left_available)
        part_info.left_mbmi = get_left_mode_info(dec_handle, mi_row, mi_col, sb_info);
    else
        part_info.left_mbmi = NULL;
    if (part_info.chroma_up_available) {
        part_info.chroma_above_mbmi =
            get_top_mode_info(dec_handle, (mi_row & (~sub_y)), (mi_col | sub_x), sb_info);
        // floored to nearest 4x4 based on sub subsampling x & y
    } else
        part_info.chroma_above_mbmi = NULL;
    if (part_info.chroma_left_available) {
        part_info.chroma_left_mbmi =
            get_left_mode_info(dec_handle, (mi_row | sub_y), (mi_col & (~sub_x)), sb_info);
        // floored to nearest 4x4 based on sub subsampling x & y
    } else
        part_info.chroma_left_mbmi = NULL;

    part_info.subsampling_x    = color_config->subsampling_x;
    part_info.subsampling_y    = color_config->subsampling_y;
    part_info.ps_global_motion = dec_handle->master_frame_buf.cur_frame_bufs[0].global_motion_warp;

    /* Wait until reference block's recon is complete for intrabc blocks */
    if (dec_handle->dec_config.threads > 1 && mode_info->use_intrabc) {
        assert(mode_info->ref_frame[1] == NONE_FRAME);
        const MV mv = mode_info->mv[0].as_mv;

        /* A negative mv.col will have ensured
        reference recon due to top sync */
        if (mv.col > 0) {
            /* mv.row is in (SUBPEL_BITS-1) Q format */
            int32_t ref_row = mi_row +
                (mv.row >> (SUBPEL_BITS-1 + MI_SIZE_LOG2));

            int32_t sb_size_log2 = dec_handle->seq_header.sb_size_log2;

            int32_t ref_sb_tile_row =
                (AOMMAX(ref_row - tile->mi_row_start, 0)) >>
                (sb_size_log2 - MI_SIZE_LOG2);
            int32_t cur_sb_tile_row = (mi_row - tile->mi_row_start) >>
                (sb_size_log2 - MI_SIZE_LOG2);

            if (ref_sb_tile_row != cur_sb_tile_row) {
                /* mv.col is in (SUBPEL_BITS-1) Q format */
                int32_t ref_col = mi_col + (mv.col >>
                    (SUBPEL_BITS-1 + MI_SIZE_LOG2));
                DecMtFrameData *dec_mt_frame_data = &dec_handle->
                    master_frame_buf.cur_frame_bufs[0].dec_mt_frame_data;

                int32_t ref_sb_tile_col = AOMMIN(ref_col, tile->mi_col_end) >>
                    (sb_size_log2 - MI_SIZE_LOG2);

                int32_t tiles_ctr = (tile->tile_row * dec_handle->
                    frame_header.tiles_info.tile_cols) + tile->tile_col;

                volatile int32_t *ref_sb_completed =
                    (volatile int32_t*) &dec_mt_frame_data->
                    parse_recon_tile_info_array[tiles_ctr].
                    sb_recon_completed_in_row[ref_sb_tile_row];
                while (*ref_sb_completed < ref_sb_tile_col + 1);
            }
        }
    }

    /* Derive warped params for local warp mode*/
    if (inter_block) {
        if (WARPED_CAUSAL == mode_info->motion_mode) {
            int32_t pts[SAMPLES_ARRAY_SIZE], pts_inref[SAMPLES_ARRAY_SIZE];
            int32_t nsamples = 0;
            int32_t apply_wm = 0;

            nsamples = find_warp_samples(dec_handle, tile, &part_info, pts, pts_inref);
            assert(nsamples > 0);

            MV mv                               = mode_info->mv[REF_LIST_0].as_mv;
            part_info.local_warp_params.wmtype  = DEFAULT_WMTYPE;
            part_info.local_warp_params.invalid = 0;

            if (nsamples > 1) nsamples = select_samples(&mv, pts, pts_inref, nsamples, bsize);

            part_info.num_samples = nsamples;

            apply_wm = !eb_find_projection(nsamples,
                                           pts,
                                           pts_inref,
                                           bsize,
                                           mv.row,
                                           mv.col,
                                           &part_info.local_warp_params,
                                           mi_row,
                                           mi_col);

            /* local warp mode should find valid projection */
            assert(apply_wm);
            part_info.local_warp_params.invalid = !apply_wm;
        }

        part_info.sf_identity = &dec_handle->sf_identity;
        for (int ref = 0; ref < 1 + has_second_ref(mode_info); ++ref) {
            const MvReferenceFrame frame = mode_info->ref_frame[ref];
            part_info.block_ref_sf[ref]  = get_ref_scale_factors(dec_handle, frame);
        }
    }

    if (inter_block)
        svtav1_predict_inter_block(
            dec_mod_ctxt, dec_handle, &part_info, mi_row, mi_col, num_planes);

    TxType           tx_type;
    int32_t *        coeffs;
    TransformInfo_t *trans_info = NULL;
    TxSize           tx_size;
    uint32_t         num_tu;

    const int max_blocks_wide = max_block_wide(&part_info, bsize, 0);
    const int max_blocks_high = max_block_high(&part_info, bsize, 0);

    uint8_t lossless       = dec_mod_ctxt->frame_header->lossless_array[part_info.mi->segment_id];
    int     lossless_block = (lossless && ((bsize >= BLOCK_64X64) && (bsize <= BLOCK_128X128)));
    int     num_chroma_tus = lossless_block
                             ? (max_blocks_wide * max_blocks_high) >>
                                   (color_config->subsampling_x + color_config->subsampling_y)
                             : mode_info->num_tus[AOM_PLANE_U];

    LfCtxt *               lf_ctxt     = (LfCtxt *)dec_handle->pv_lf_ctxt;
    int32_t                lf_stride   = dec_mod_ctxt->frame_header->mi_stride;

    for (int plane = 0; plane < num_planes; ++plane) {
        sub_x = (plane > 0) ? color_config->subsampling_x : 0;
        sub_y = (plane > 0) ? color_config->subsampling_y : 0;

        if (plane && !is_chroma_ref) continue;

        trans_info = (plane == 2)
                         ? (sb_info->sb_trans_info[plane - 1] +
                            mode_info->first_txb_offset[plane - 1] + num_chroma_tus)
                         : (sb_info->sb_trans_info[plane] + mode_info->first_txb_offset[plane]);

        if (lossless_block) {
            assert(trans_info->tx_size == TX_4X4);
            num_tu = (max_blocks_wide * max_blocks_high) >> (sub_x + sub_y);

        } else
            num_tu = mode_info->num_tus[!!plane];

        assert(num_tu != 0);
        void *  blk_recon_buf;
        int32_t recon_stride;

        derive_blk_pointers(recon_picture_buf,
                            plane,
                            (mi_col >> sub_x) << MI_SIZE_LOG2,
                            (mi_row >> sub_y) << MI_SIZE_LOG2,
                            &blk_recon_buf,
                            &recon_stride,
                            sub_x,
                            sub_y);

        for (uint32_t tu = 0; tu < num_tu; tu++) {
            void *  txb_recon_buf;
            int32_t txb_offset;

            tx_size = trans_info->tx_size;
            coeffs  = dec_mod_ctxt->cur_coeff[plane];

            txb_offset = (trans_info->txb_y_offset * recon_stride + trans_info->txb_x_offset)
                         << MI_SIZE_LOG2;
            txb_recon_buf = (void *)((uint8_t *)blk_recon_buf + (txb_offset << hbd));

            if (dec_handle->is_lf_enabled) {
                if(plane != 2)
                    fill_4x4_lf_param(lf_ctxt,
                                     (mi_col & (~sub_x)) +
                                     (trans_info->txb_x_offset << sub_x),
                                     (mi_row & (~sub_y)) +
                                     (trans_info->txb_y_offset << sub_y),
                                     lf_stride, tx_size, sub_x, sub_y, plane);
            }

            if (!inter_block)
                svt_av1_predict_intra(dec_mod_ctxt,
                                      &part_info,
                                      plane,
                                      tx_size,
                                      tile,
                                      txb_recon_buf,
                                      recon_stride,
                                      recon_picture_buf->bit_depth,
                                      trans_info->txb_x_offset,
                                      trans_info->txb_y_offset);

            n_coeffs = 0;

            if (!mode_info->skip && trans_info->cbf) {
                int32_t *qcoeffs = dec_mod_ctxt->iquant_cur_ptr;
                int32_t  iq_size = tx_size_wide[tx_size] * tx_size_high[tx_size];
                memset(dec_mod_ctxt->iquant_cur_ptr,
                       0,
                       iq_size * sizeof(*dec_mod_ctxt->iquant_cur_ptr));
                dec_mod_ctxt->iquant_cur_ptr = dec_mod_ctxt->iquant_cur_ptr + iq_size;
#if SVT_DEC_COEFF_DEBUG
                {
                    /* For debug purpose */
                    uint8_t *cur_coeff  = (uint8_t *)coeffs;
                    uint8_t  mi_row_der = cur_coeff[0];
                    uint8_t  mi_col_der = cur_coeff[1];

                    uint8_t cur_loc = (mi_row + trans_info->txb_y_offset) & 0xFF;
                    assert(mi_row_der == cur_loc);
                    cur_loc = (mi_col + trans_info->txb_x_offset) & 0xFF;
                    assert(mi_col_der == cur_loc);
                    (void)mi_row_der;
                    (void)mi_col_der;
                    (void)cur_loc;
                }
#endif
                tx_type = trans_info->txk_type;

                n_coeffs = inverse_quantize(
                    dec_mod_ctxt, &part_info, mode_info, coeffs, qcoeffs, tx_type, tx_size, plane);
                if (n_coeffs != 0) {
                    dec_mod_ctxt->cur_coeff[plane] += (n_coeffs + 1);

                    if (recon_picture_buf->bit_depth == EB_8BIT && !is16b)
                        av1_inv_transform_recon8bit(qcoeffs,
                                                    (uint8_t *)txb_recon_buf,
                                                    recon_stride,
                                                    (uint8_t *)txb_recon_buf,
                                                    recon_stride,
                                                    tx_size,
                                                    tx_type,
                                                    plane,
                                                    n_coeffs,
                                                    lossless);
                    else
                        av1_inv_transform_recon(qcoeffs,
                                                CONVERT_TO_BYTEPTR(txb_recon_buf),
                                                recon_stride,
                                                CONVERT_TO_BYTEPTR(txb_recon_buf),
                                                recon_stride,
                                                tx_size,
                                                recon_picture_buf->bit_depth,
                                                tx_type,
                                                plane,
                                                n_coeffs,
                                                lossless);
                }
            }

            /* Store Luma for CFL if required! */
            if (plane == AOM_PLANE_Y &&
                store_cfl_required(color_config, &part_info, is_chroma_ref)) {
                cfl_store_tx(&part_info,
                             &dec_mod_ctxt->cfl_ctx,
                             trans_info->txb_y_offset,
                             trans_info->txb_x_offset,
                             tx_size,
                             bsize,
                             color_config,
                             txb_recon_buf,
                             recon_stride,
                             is16b);
            }

            // increment transform pointer
            trans_info++;
        }
    }
    return;
}
