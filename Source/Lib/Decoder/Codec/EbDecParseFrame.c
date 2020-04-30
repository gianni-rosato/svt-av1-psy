/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

// SUMMARY
//   Contains the Parse related functions

/**************************************
 * Includes
 **************************************/

#include "EbDefinitions.h"
#include "EbPictureBufferDesc.h"

#include "EbDecHandle.h"
#include "EbDecBitReader.h"
#include "EbDecProcessFrame.h"
#include "EbDecUtils.h"

#include "EbDecParseFrame.h"
#include "EbDecParseHelper.h"

/* Inititalizes prms for current tile from Master TilesInfo ! */
void svt_tile_init(TileInfo *cur_tile_info, FrameHeader *frame_header, int32_t tile_row,
                   int32_t tile_col) {
    TilesInfo *tiles_info = &frame_header->tiles_info;

    /* tile_set_row */
    assert(tile_row < tiles_info->tile_rows);
    cur_tile_info->tile_row     = tile_row;
    cur_tile_info->mi_row_start = tiles_info->tile_row_start_mi[tile_row];
    cur_tile_info->mi_row_end =
        AOMMIN(tiles_info->tile_row_start_mi[tile_row + 1], frame_header->mi_rows);

    assert(cur_tile_info->mi_row_end > cur_tile_info->mi_row_start);

    /* tile_set_col */
    assert(tile_col < tiles_info->tile_cols);

    cur_tile_info->tile_col     = tile_col;
    cur_tile_info->mi_col_start = tiles_info->tile_col_start_mi[tile_col];
    cur_tile_info->mi_col_end =
        AOMMIN(tiles_info->tile_col_start_mi[tile_col + 1], frame_header->mi_cols);

    assert(cur_tile_info->mi_col_end > cur_tile_info->mi_col_start);
}

static int read_is_valid(const uint8_t *start, size_t len, const uint8_t *end) {
    return len != 0 && len <= (size_t)(end - start);
}

EbErrorType init_svt_reader(SvtReader *r, const uint8_t *data, const uint8_t *data_end,
                            const size_t read_size, uint8_t allow_update_cdf) {
    EbErrorType status = EB_ErrorNone;

    if (read_is_valid(data, read_size, data_end)) {
        if (0 == svt_reader_init(r, data, read_size))
            r->allow_update_cdf = allow_update_cdf;
        else
            status = EB_Corrupt_Frame;
    } else
        status = EB_Corrupt_Frame;
    return status;
}

void clear_above_context(ParseCtxt *parse_ctxt, int mi_col_start, int mi_col_end, int num_threads) {
    SeqHeader *seq_params = parse_ctxt->seq_header;
    int        num_planes = av1_num_planes(&seq_params->color_config);
    int        width      = mi_col_end - mi_col_start;
    /* ToDo: Bhavna : Can be optimized for ST */
    if (num_threads == 1) {
        int32_t num_mi_sb        = seq_params->sb_mi_size;
        int32_t sb_size_log2     = seq_params->sb_size_log2;
        int32_t sb_aligned_width = ALIGN_POWER_OF_TWO(seq_params->max_frame_width, sb_size_log2);
        int32_t sb_cols          = sb_aligned_width >> sb_size_log2;
        width                    = sb_cols * num_mi_sb;
    }

    int8_t num4_64x64 = mi_size_wide[BLOCK_64X64];

    for (int i = 0; i < num_planes; i++) {
        ZERO_ARRAY(parse_ctxt->parse_above_nbr4x4_ctxt->above_ctx[i], width);
        ZERO_ARRAY((parse_ctxt->parse_above_nbr4x4_ctxt->above_palette_colors[i]),
                   num4_64x64 * PALETTE_MAX_SIZE);
    }

    ZERO_ARRAY(parse_ctxt->parse_above_nbr4x4_ctxt->above_seg_pred_ctx, width);

    ZERO_ARRAY(parse_ctxt->parse_above_nbr4x4_ctxt->above_part_wd, width);

    ZERO_ARRAY(parse_ctxt->parse_above_nbr4x4_ctxt->above_comp_grp_idx, width);

    memset(parse_ctxt->parse_above_nbr4x4_ctxt->above_tx_wd,
           tx_size_wide[TX_SIZES_LARGEST],
           width * sizeof(uint8_t));
}

void clear_left_context(ParseCtxt *parse_ctxt) {
    SeqHeader *seq_params = parse_ctxt->seq_header;

    /* Maintained only for 1 left SB! */
    int     blk_cnt          = seq_params->sb_mi_size;
    int     num_planes       = av1_num_planes(&seq_params->color_config);
    int32_t num_4x4_neigh_sb = seq_params->sb_mi_size;

    /* TODO :  after Optimizing the allocation for Chroma fix here also */
    for (int i = 0; i < num_planes; i++)
        ZERO_ARRAY(parse_ctxt->parse_left_nbr4x4_ctxt->left_ctx[i], blk_cnt);

    ZERO_ARRAY(parse_ctxt->parse_left_nbr4x4_ctxt->left_seg_pred_ctx, blk_cnt);

    ZERO_ARRAY(parse_ctxt->parse_left_nbr4x4_ctxt->left_part_ht, blk_cnt);

    ZERO_ARRAY(parse_ctxt->parse_left_nbr4x4_ctxt->left_comp_grp_idx, blk_cnt);

    for (int i = 0; i < num_planes; i++) {
        ZERO_ARRAY((parse_ctxt->parse_left_nbr4x4_ctxt->left_palette_colors[i]),
                   num_4x4_neigh_sb * PALETTE_MAX_SIZE);
    }

    memset(parse_ctxt->parse_left_nbr4x4_ctxt->left_tx_ht,
           tx_size_high[TX_SIZES_LARGEST],
           blk_cnt * sizeof(parse_ctxt->parse_left_nbr4x4_ctxt->left_tx_ht[0]));
}

void clear_cdef(int8_t *sb_cdef_strength, int32_t cdef_factor) {
    memset(sb_cdef_strength, -1, cdef_factor * sizeof(*sb_cdef_strength));
}

void clear_loop_filter_delta(ParseCtxt *parse_ctx) {
    for (int lf_id = 0; lf_id < FRAME_LF_COUNT; ++lf_id) parse_ctx->delta_lf[lf_id] = 0;
}

EbErrorType start_parse_tile(EbDecHandle *dec_handle_ptr, ParseCtxt *parse_ctxt,
                             TilesInfo *tiles_info, int tile_num, int is_mt) {
    MasterParseCtxt *master_parse_ctxt = (MasterParseCtxt *)dec_handle_ptr->pv_master_parse_ctxt;
    FrameHeader *    frame_header      = &dec_handle_ptr->frame_header;
    ParseTileData *  parse_tile_data   = &master_parse_ctxt->parse_tile_data[tile_num];
    EbErrorType      status            = EB_ErrorNone;
    int              tile_row          = tile_num / tiles_info->tile_cols;
    int              tile_col          = tile_num % tiles_info->tile_cols;
    svt_tile_init(&parse_ctxt->cur_tile_info, frame_header, tile_row, tile_col);

    parse_ctxt->cur_q_ind = frame_header->quantization_params.base_q_idx;

    status = init_svt_reader(&parse_ctxt->r,
                             (const uint8_t *)parse_tile_data->data,
                             parse_tile_data->data_end,
                             parse_tile_data->tile_size,
                             !(parse_ctxt->frame_header->disable_cdf_update));
    if (status != EB_ErrorNone) return status;

    parse_ctxt->cur_tile_ctx = master_parse_ctxt->init_frm_ctx;

    /* Parse Tile */
    status =
        parse_tile(dec_handle_ptr, parse_ctxt, tiles_info, tile_num, tile_row, tile_col, is_mt);

    /* Save CDF */
    if (!frame_header->disable_frame_end_update_cdf &&
        (tile_num == tiles_info->context_update_tile_id)) {
        dec_handle_ptr->cur_pic_buf[0]->final_frm_ctx = parse_ctxt->cur_tile_ctx;
        eb_av1_reset_cdf_symbol_counters(&dec_handle_ptr->cur_pic_buf[0]->final_frm_ctx);
    }
    return status;
}

EbErrorType parse_tile(EbDecHandle *dec_handle_ptr, ParseCtxt *parse_ctx, TilesInfo *tile_info,
                       int tile_num, int32_t tile_row, int32_t tile_col, int32_t is_mt) {
    EbErrorType status = EB_ErrorNone;

    EbColorConfig *color_config = &dec_handle_ptr->seq_header.color_config;
    int            num_planes   = av1_num_planes(color_config);

    clear_above_context(parse_ctx,
                        tile_info->tile_col_start_mi[tile_col],
                        tile_info->tile_col_start_mi[tile_col + 1],
                        dec_handle_ptr->dec_config.threads);
    clear_loop_filter_delta(parse_ctx);

    /* Init ParseCtxt */
    RestorationUnitInfo *lr_unit[MAX_MB_PLANE];

    // Default initialization of Wiener and SGR Filter
    for (int p = 0; p < num_planes; ++p) {
        lr_unit[p] = &parse_ctx->ref_lr_unit[p];

        set_default_wiener(&lr_unit[p]->wiener_info);
        set_default_sgrproj(&lr_unit[p]->sgrproj_info);
    }

    // to-do access to wiener info that is currently part of PartitionInfo
    int32_t sb_row_tile_start = 0;
    if (is_mt) {
        DecMtParseReconTileInfo *parse_recon_tile_info =
            &dec_handle_ptr->master_frame_buf.cur_frame_bufs[0]
                 .dec_mt_frame_data.parse_recon_tile_info_array[tile_num];

        sb_row_tile_start = (parse_recon_tile_info->tile_info.mi_row_start << MI_SIZE_LOG2) >>
                            dec_handle_ptr->seq_header.sb_size_log2;
    }
    for (uint32_t mi_row = tile_info->tile_row_start_mi[tile_row];
         mi_row < tile_info->tile_row_start_mi[tile_row + 1];
         mi_row += dec_handle_ptr->seq_header.sb_mi_size) {
        int32_t sb_row = (mi_row << MI_SIZE_LOG2) >> dec_handle_ptr->seq_header.sb_size_log2;

        clear_left_context(parse_ctx);

        /*TODO: Move CFL to thread ctxt! We need to access DecModCtxt
          from parse_tile function . Add tile level cfl init. */
        if (!is_mt) {
            cfl_init(&((DecModCtxt *)dec_handle_ptr->pv_dec_mod_ctxt)->cfl_ctx, color_config);
        }

        for (uint32_t mi_col = tile_info->tile_col_start_mi[tile_col];
             mi_col < tile_info->tile_col_start_mi[tile_col + 1];
             mi_col += dec_handle_ptr->seq_header.sb_mi_size)

        {
            int32_t sb_col = (mi_col << MI_SIZE_LOG2) >> dec_handle_ptr->seq_header.sb_size_log2;
            uint8_t sx     = color_config->subsampling_x;
            uint8_t sy     = color_config->subsampling_y;

            //clear_block_decoded_flags(r, c, sbSize4)
            MasterFrameBuf *master_frame_buf = &dec_handle_ptr->master_frame_buf;
            CurFrameBuf *   frame_buf        = &master_frame_buf->cur_frame_bufs[0];
            int32_t         num_mis_in_sb    = master_frame_buf->num_mis_in_sb;

            SBInfo *sb_info = frame_buf->sb_info + (sb_row * master_frame_buf->sb_cols) + sb_col;
            *(master_frame_buf->frame_mi_map.pps_sb_info +
              sb_row * master_frame_buf->frame_mi_map.sb_cols + sb_col) = sb_info;
            sb_info->sb_mode_info                                       = frame_buf->mode_info +
                                    (sb_row * num_mis_in_sb * master_frame_buf->sb_cols) +
                                    sb_col * num_mis_in_sb;

            sb_info->sb_trans_info[AOM_PLANE_Y] =
                frame_buf->trans_info[AOM_PLANE_Y] +
                (sb_row * num_mis_in_sb * master_frame_buf->sb_cols) + sb_col * num_mis_in_sb;

            sb_info->sb_trans_info[AOM_PLANE_U] =
                frame_buf->trans_info[AOM_PLANE_U] +
                ((sb_row * num_mis_in_sb * master_frame_buf->sb_cols >> sy) +
                 (sb_col * num_mis_in_sb >> sx)) *
                    2;
            if (dec_handle_ptr->dec_config.threads == 1) {
                /*TODO : Change to macro */
                sb_info->sb_coeff[AOM_PLANE_Y] = frame_buf->coeff[AOM_PLANE_Y];
                sb_info->sb_coeff[AOM_PLANE_U] = frame_buf->coeff[AOM_PLANE_U];
                sb_info->sb_coeff[AOM_PLANE_V] = frame_buf->coeff[AOM_PLANE_V];
            }
            else {
                /*TODO : Change to macro */
                sb_info->sb_coeff[AOM_PLANE_Y] =
                    frame_buf->coeff[AOM_PLANE_Y] +
                    (sb_row * num_mis_in_sb * master_frame_buf->sb_cols * (16 + 1))
                        + sb_col * num_mis_in_sb* (16 + 1);
                /*TODO : Change to macro */
                sb_info->sb_coeff[AOM_PLANE_U] =
                    frame_buf->coeff[AOM_PLANE_U] +
                    (sb_row * num_mis_in_sb * master_frame_buf->sb_cols * (16 + 1) >> (sy + sx))
                        + (sb_col * num_mis_in_sb * (16 + 1) >> (sy + sx));
                sb_info->sb_coeff[AOM_PLANE_V] =
                    frame_buf->coeff[AOM_PLANE_V] +
                    (sb_row * num_mis_in_sb * master_frame_buf->sb_cols * (16 + 1) >> (sy + sx))
                        + (sb_col * num_mis_in_sb * (16 + 1) >> (sy + sx));
            }
            int cdef_factor = dec_handle_ptr->seq_header.use_128x128_superblock ? 4 : 1;
            sb_info->sb_cdef_strength =
                frame_buf->cdef_strength +
                (((sb_row * master_frame_buf->sb_cols) + sb_col) * cdef_factor);

            sb_info->sb_delta_lf =
                frame_buf->delta_lf +
                (FRAME_LF_COUNT * ((sb_row * master_frame_buf->sb_cols) + sb_col));

            sb_info->sb_delta_q =
                frame_buf->delta_q + (sb_row * master_frame_buf->sb_cols) + sb_col;

            clear_cdef(sb_info->sb_cdef_strength, cdef_factor);

            parse_ctx->first_txb_offset[AOM_PLANE_Y] = 0;
            parse_ctx->first_txb_offset[AOM_PLANE_U] = 0;
            parse_ctx->cur_mode_info                 = sb_info->sb_mode_info;
            parse_ctx->cur_mode_info_cnt             = 0;
            parse_ctx->sb_row_mi                     = mi_row;
            parse_ctx->sb_col_mi                     = mi_col;
            parse_ctx->cur_coeff_buf[AOM_PLANE_Y]    = sb_info->sb_coeff[AOM_PLANE_Y];
            parse_ctx->cur_coeff_buf[AOM_PLANE_U]    = sb_info->sb_coeff[AOM_PLANE_U];
            parse_ctx->cur_coeff_buf[AOM_PLANE_V]    = sb_info->sb_coeff[AOM_PLANE_V];
            parse_ctx->prev_blk_has_chroma           = 1; //default at start of frame / tile

            sb_info->num_block = 0;
            // Bit-stream parsing of the superblock
            parse_super_block(dec_handle_ptr, parse_ctx, mi_row, mi_col, sb_info);

            if (!is_mt) {
                /* Init DecModCtxt */
                DecModCtxt *dec_mod_ctxt = (DecModCtxt *)dec_handle_ptr->pv_dec_mod_ctxt;
                dec_mod_ctxt->cur_coeff[AOM_PLANE_Y] = sb_info->sb_coeff[AOM_PLANE_Y];
                dec_mod_ctxt->cur_coeff[AOM_PLANE_U] = sb_info->sb_coeff[AOM_PLANE_U];
                dec_mod_ctxt->cur_coeff[AOM_PLANE_V] = sb_info->sb_coeff[AOM_PLANE_V];

                dec_mod_ctxt->cur_tile_info = parse_ctx->cur_tile_info;

                /* TO DO : Will move later */
                // decoding of the superblock
                decode_super_block(dec_mod_ctxt, mi_row, mi_col, sb_info);
            }
        }
        if (is_mt) {
            DecMtFrameData *dec_mt_frame_data =
                &dec_handle_ptr->master_frame_buf.cur_frame_bufs[0]
                     .dec_mt_frame_data; //multi frame Parallel 0 -> idx
            assert(sb_row >= sb_row_tile_start);
            dec_mt_frame_data->parse_recon_tile_info_array[tile_num]
                .sb_recon_row_parsed[sb_row - sb_row_tile_start] = 1;
        }
    }

    return status;
}
