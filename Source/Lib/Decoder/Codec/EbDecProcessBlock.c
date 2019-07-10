/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

// SUMMARY
//   Contains the Decode Block related functions

/**************************************
 * Includes
 **************************************/

#include "EbDefinitions.h"

#include "EbSvtAv1Dec.h"
#include "EbDecHandle.h"

#include "EbObuParse.h"
#include "EbDecParseHelper.h"

#include "EbDecInverseQuantize.h"
#include "EbDecProcessFrame.h"

#include "EbDecIntraPrediction.h"
#include "EbDecInterPrediction.h"

#include "EbDecNbr.h"
#include "EbDecUtils.h"

#include "EbTransforms.h"

extern int select_samples(
    MV *mv,
    int *pts,
    int *pts_inref,
    int len,
    BlockSize bsize);

CflAllowedType store_cfl_required(const EbColorConfig *cc,
                                  const PartitionInfo_t  *xd)

{
    const ModeInfo_t *mbmi = xd->mi;

    if (cc->mono_chrome) return CFL_DISALLOWED;

    if (!xd->has_chroma) {
        // For non-chroma-reference blocks, we should always store the luma pixels,
        // in case the corresponding chroma-reference block uses CfL.
        // Note that this can only happen for block sizes which are <8 on
        // their shortest side, as otherwise they would be chroma reference
        // blocks.
        return CFL_ALLOWED;
    }

    // If this block has chroma information, we know whether we're
    // actually going to perform a CfL prediction
    return (CflAllowedType)(!dec_is_inter_block(mbmi) &&
        mbmi->uv_mode == UV_CFL_PRED);
}

 /* decode partition */
PartitionType get_partition(DecModCtxt *dec_mod_ctxt, FrameHeader *frame_header,
                            uint32_t mi_row, uint32_t mi_col, SBInfo *sb_info,
                            BlockSize bsize)
{
    if (mi_row >= frame_header->mi_rows || mi_col >= frame_header->mi_cols)
        return PARTITION_INVALID;

    ModeInfo_t *mode_info = get_cur_mode_info(dec_mod_ctxt->dec_handle_ptr, mi_row, mi_col, sb_info);

    const BlockSize subsize = mode_info->sb_type;

    if (subsize == bsize) return PARTITION_NONE;

    const int bhigh = mi_size_high[bsize];
    const int bwide = mi_size_wide[bsize];
    const int sshigh = mi_size_high[subsize];
    const int sswide = mi_size_wide[subsize];

    if (bsize > BLOCK_8X8 && mi_row + bwide / 2 < frame_header->mi_rows &&
        mi_col + bhigh / 2 < frame_header->mi_cols)
    {
        // In this case, the block might be using an extended partition type.
        /* TODO: Fix the nbr access! */
        const ModeInfo_t *const mbmi_right = get_cur_mode_info(dec_mod_ctxt->dec_handle_ptr,
            mi_row, mi_col + (bwide / 2), sb_info);
        const ModeInfo_t *const mbmi_below = get_cur_mode_info(dec_mod_ctxt->dec_handle_ptr,
            mi_row + (bhigh / 2), mi_col, sb_info);

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
        }
        else if (sshigh == bhigh) {
            // Smaller width but same height. Is PARTITION_VERT_4, PARTITION_VERT or
            // PARTITION_VERT_B. To distinguish the latter two, check if the right
            // half was split.
            if (sswide * 4 == bwide) return PARTITION_VERT_4;
            assert(sswide * 2 == bhigh);

            if (mbmi_right->sb_type == subsize)
                return PARTITION_VERT;
            else
                return PARTITION_VERT_B;
        }
        else {
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
    const int split_idx = (vert_split << 1) | horz_split;
    assert(split_idx != 0);

    static const PartitionType base_partitions[4] = {
      PARTITION_INVALID, PARTITION_HORZ, PARTITION_VERT, PARTITION_SPLIT
    };

    return base_partitions[split_idx];
}

void decode_block(DecModCtxt *dec_mod_ctxt, int32_t mi_row, int32_t mi_col,
    BlockSize bsize, TileInfo *tile, SBInfo *sb_info)
{
    EbDecHandle *dec_handle   = (EbDecHandle *)dec_mod_ctxt->dec_handle_ptr;
    EbColorConfig *color_config = &dec_handle->seq_header.color_config;
    EbPictureBufferDesc *recon_picture_buf = dec_handle->cur_pic_buf[0]->ps_pic_buf;
    uint32_t mi_cols = (&dec_handle->frame_header)->mi_cols;
    uint32_t mi_rows = (&dec_handle->frame_header)->mi_rows;

    int num_planes = av1_num_planes(color_config);

    ModeInfo_t *mode_info = get_cur_mode_info(dec_mod_ctxt->dec_handle_ptr, mi_row, mi_col, sb_info);
    bool inter_block = dec_is_inter_block(mode_info);
#if MODE_INFO_DBG
    assert(mode_info->mi_row == mi_row);
    assert(mode_info->mi_col == mi_col);
#endif
    int32_t bw4 = mi_size_wide[bsize];
    int32_t bh4 = mi_size_high[bsize];

    int sub_x, sub_y, n_coeffs;
    sub_x = color_config->subsampling_x;
    sub_y = color_config->subsampling_y;
    int is_chroma_reference = dec_is_chroma_reference(mi_row, mi_col, bsize,
        sub_x, sub_y);

    /* TODO: Can move to a common init fun for parse & decode */
    PartitionInfo_t part_info;
    part_info.mi = mode_info;
    part_info.sb_info = sb_info;
    part_info.mi_row = mi_row;
    part_info.mi_col = mi_col;
    part_info.has_chroma = (int8_t) is_chroma_reference;
    part_info.pv_cfl_ctxt = &dec_mod_ctxt->cfl_ctx;

    /*!< Distance of MB away from frame edges in subpixels (1/8th pixel). */
    part_info.mb_to_top_edge = -((mi_row * MI_SIZE) * 8);
    part_info.mb_to_bottom_edge = ((mi_rows - bh4 - mi_row) * MI_SIZE) * 8;
    part_info.mb_to_left_edge = -((mi_col * MI_SIZE) * 8);
    part_info.mb_to_right_edge = ((mi_cols - bw4 - mi_col) * MI_SIZE) * 8;

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
    part_info.up_available = (mi_row > tile->mi_row_start);
    part_info.left_available = (mi_col > tile->mi_col_start);
    part_info.chroma_up_available = part_info.up_available;
    part_info.chroma_left_available = part_info.left_available;

    if (part_info.has_chroma) {
        if (bh4 == 1 && sub_y)
            part_info.chroma_up_available = (mi_row - 1) > tile->mi_row_start;
        if (bw4 == 1 && sub_x)
            part_info.chroma_left_available = (mi_col - 1) > tile->mi_col_start;
    } else {
        part_info.chroma_up_available = 0;
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
        part_info.chroma_above_mbmi = get_top_mode_info
            (dec_handle, (mi_row & (~sub_x)), (mi_col | sub_y), sb_info); // floored to nearest 4x4 based on sub subsampling x & y
    }
    else
        part_info.chroma_above_mbmi = NULL;
    if (part_info.chroma_left_available) {
        part_info.chroma_left_mbmi = get_left_mode_info
            (dec_handle, (mi_row | sub_x), (mi_col & (~sub_y)), sb_info); // floored to nearest 4x4 based on sub subsampling x & y
    }
    else
        part_info.chroma_left_mbmi = NULL;

    part_info.subsampling_x = color_config->subsampling_x;
    part_info.subsampling_y = color_config->subsampling_y;
    part_info.ps_global_motion = dec_handle->master_frame_buf.cur_frame_bufs[0].global_motion_warp;

    /* Derive warped params for local warp mode*/
    if (WARPED_CAUSAL == mode_info->motion_mode && inter_block) {

        int32_t pts[SAMPLES_ARRAY_SIZE], pts_inref[SAMPLES_ARRAY_SIZE];
        int32_t nsamples = 0;
        int32_t apply_wm = 0;

        nsamples = find_warp_samples(dec_handle, &part_info, mi_row, mi_col, pts, pts_inref);
        assert(nsamples > 0);

        MV mv = mode_info->mv[REF_LIST_0].as_mv;
        part_info.local_warp_params.wmtype = DEFAULT_WMTYPE;
        part_info.local_warp_params.invalid = 0;

        if (nsamples > 1)
            nsamples = select_samples(&mv, pts, pts_inref, nsamples, bsize);

        part_info.num_samples = nsamples;

        apply_wm = !find_projection(
            nsamples,
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

    if (inter_block)
        svtav1_predict_inter_block(dec_handle, &part_info, mi_row, mi_col,
            num_planes);

    int32_t *qcoeffs = dec_mod_ctxt->sb_iquant_ptr;
    TxType tx_type;
    int32_t *coeffs;
    TransformInfo_t *trans_info = NULL;
    TxSize tx_size;
    int num_tu;

    const int max_blocks_wide = max_block_wide(&part_info, bsize, 0);
    const int max_blocks_high = max_block_high(&part_info, bsize, 0);

    int num_chroma_tus = (dec_handle->frame_header.lossless_array[part_info.mi->segment_id] &&
        ((bsize >= BLOCK_64X64) && (bsize <= BLOCK_128X128)) ) ?
        (max_blocks_wide * max_blocks_high) >>
        (color_config->subsampling_x + color_config->subsampling_y) : mode_info->num_chroma_tus;

    for (int plane = 0; plane < num_planes; ++plane) {
        sub_x = (plane > 0) ? color_config->subsampling_x : 0;
        sub_y = (plane > 0) ? color_config->subsampling_y : 0;

        if (!dec_is_chroma_reference(mi_row, mi_col, bsize, sub_x, sub_y))
            continue;

        trans_info = (plane == 0) ?
            (sb_info->sb_trans_info[plane] + mode_info->first_luma_tu_offset) :
            (plane == 1) ?
            (sb_info->sb_trans_info[plane] + mode_info->first_chroma_tu_offset) :
            (sb_info->sb_trans_info[plane - 1] + mode_info->first_chroma_tu_offset + num_chroma_tus);

        if (dec_handle->frame_header.lossless_array[part_info.mi->segment_id] &&
            (bsize >= BLOCK_64X64) && (bsize <= BLOCK_128X128) )
        {
            assert(trans_info->tx_size == TX_4X4);
            num_tu = (max_blocks_wide * max_blocks_high) >> (sub_x + sub_y);

        }
        else
            num_tu = plane ? mode_info->num_chroma_tus : mode_info->num_luma_tus;

        assert(num_tu != 0);

        for (int tu = 0; tu < num_tu; tu++)
        {
            void *blk_recon_buf;
            int32_t recon_strd;

            tx_size = trans_info->tx_size;
            coeffs = dec_mod_ctxt->cur_coeff[plane];

            derive_blk_pointers(recon_picture_buf, plane,
                ((mi_col >> sub_x) + trans_info->tu_x_offset)*MI_SIZE,
                ((mi_row >> sub_y) + trans_info->tu_y_offset)*MI_SIZE,
                &blk_recon_buf, &recon_strd, sub_x, sub_y);

            if (!inter_block)
                svt_av1_predict_intra(dec_mod_ctxt, &part_info, plane,
                    tx_size, dec_mod_ctxt->cur_tile_info, blk_recon_buf,
                    recon_strd, recon_picture_buf->bit_depth,
                    trans_info->tu_x_offset, trans_info->tu_y_offset);

            n_coeffs = 0;

            if (!mode_info->skip && trans_info->cbf)
            {
#if SVT_DEC_COEFF_DEBUG
                {
                    /* For debug purpose */
                    uint8_t    *cur_coeff = (uint8_t*)coeffs;
                    uint8_t mi_row_der = cur_coeff[0];
                    uint8_t mi_col_der = cur_coeff[1];

                    uint8_t  cur_loc = (mi_row + trans_info->tu_y_offset) & 0xFF;
                    assert(mi_row_der == cur_loc);
                    cur_loc = (mi_col + trans_info->tu_x_offset) & 0xFF;
                    assert(mi_col_der == cur_loc);
                    (void)mi_row_der; (void)mi_col_der; (void)cur_loc;
                }
#endif
                tx_type = trans_info->txk_type;

                n_coeffs = inverse_quantize(dec_handle, &part_info,
                    mode_info, coeffs, qcoeffs, tx_type, tx_size, plane);
                if (n_coeffs != 0) {
                    dec_mod_ctxt->cur_coeff[plane] += (n_coeffs + 1);

                    if (recon_picture_buf->bit_depth == EB_8BIT)
                        av1_inv_transform_recon8bit(qcoeffs,
                        (uint8_t *)blk_recon_buf,
                            recon_strd, tx_size, tx_type, plane, n_coeffs);
                    else
                        av1_inv_transform_recon(qcoeffs,
                            CONVERT_TO_BYTEPTR(blk_recon_buf), recon_strd,
                            tx_size, recon_picture_buf->bit_depth,
                            tx_type, plane, n_coeffs);
                }
            }

            /* Store Luma for CFL if required! */
            if (plane == AOM_PLANE_Y && store_cfl_required(
                color_config, &part_info))
            {
                cfl_store_tx(&part_info, &dec_mod_ctxt->cfl_ctx,
                    trans_info->tu_y_offset, trans_info->tu_x_offset, tx_size,
                    bsize, color_config, blk_recon_buf, recon_strd);
            }

            // increment transform pointer
            trans_info++;
        }
    }
    return;
}
