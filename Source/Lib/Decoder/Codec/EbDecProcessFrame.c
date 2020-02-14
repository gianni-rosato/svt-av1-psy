/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

// SUMMARY
//   Contains the Decode related functions

/**************************************
 * Includes
 **************************************/

#include "EbDefinitions.h"

#include "EbSvtAv1Dec.h"
#include "EbDecHandle.h"

#include "EbDecInverseQuantize.h"
#include "EbDecProcessFrame.h"
#include "EbDecProcessBlock.h"
#include "EbDecNbr.h"
#include "EbUtility.h"

/* decode partition */
static void decode_partition(DecModCtxt *dec_mod_ctxt,
                             uint32_t mi_row, uint32_t mi_col,
                             SBInfo *sb_info) {
    BlockSize     subsize;

    EbDecHandle *dec_handle = (EbDecHandle *)dec_mod_ctxt->dec_handle_ptr;

    if (mi_row >= dec_mod_ctxt->frame_header->mi_rows ||
        mi_col >= dec_mod_ctxt->frame_header->mi_cols)
        return;

    BlockModeInfo *mode_info = get_cur_mode_info(dec_handle,
                                                 mi_row, mi_col, sb_info);

    int n_blocks = sb_info->num_block;
    int sub_mi_row = 0;
    int sub_mi_col = 0;

    for (int i = 0; i < n_blocks; i++) {
        sub_mi_row = mode_info->mi_row_in_sb;
        sub_mi_col = mode_info->mi_col_in_sb;
        subsize = mode_info->sb_type;
        decode_block(dec_mod_ctxt,
                     mode_info,
                     mi_row + sub_mi_row,
                     mi_col + sub_mi_col,
                     subsize,
                     &dec_mod_ctxt->cur_tile_info,
                     sb_info);
        mode_info++;
    }
}

// decoding of the superblock
void decode_super_block(DecModCtxt *dec_mod_ctxt, uint32_t mi_row, uint32_t mi_col,
                        SBInfo *sb_info) {
    dec_mod_ctxt->iquant_cur_ptr = dec_mod_ctxt->sb_iquant_ptr;

    /* SB level dequant update */
    update_dequant(dec_mod_ctxt, sb_info);

    /* Decode partition */
    decode_partition(dec_mod_ctxt, mi_row, mi_col, sb_info);
}

EbErrorType decode_tile_row(DecModCtxt *dec_mod_ctxt, TilesInfo *tile_info,
                            DecMtParseReconTileInfo *parse_recon_tile_info_array, int32_t tile_col,
                            int32_t mi_row, int32_t sb_row) {
    EbErrorType       status                   = EB_ErrorNone;
    EbDecHandle *     dec_handle_ptr           = (EbDecHandle *)(dec_mod_ctxt->dec_handle_ptr);
    MasterFrameBuf *  master_frame_buf         = &dec_handle_ptr->master_frame_buf;
    CurFrameBuf *     frame_buf                = &master_frame_buf->cur_frame_bufs[0];
    volatile int32_t *sb_completed_in_prev_row = NULL;
    uint32_t *        sb_completed_in_row;
    int32_t           tile_wd_in_sb;
    int32_t           sb_mi_size_log2 = dec_mod_ctxt->seq_header->sb_size_log2 - MI_SIZE_LOG2;

    int32_t sb_row_tile_start =
        (parse_recon_tile_info_array->tile_info.mi_row_start << MI_SIZE_LOG2) >>
        dec_mod_ctxt->seq_header->sb_size_log2;

    int32_t sb_row_in_tile = sb_row - sb_row_tile_start;

    if (0 != sb_row_in_tile) {
        sb_completed_in_prev_row = (volatile int32_t *)&parse_recon_tile_info_array
                                       ->sb_recon_completed_in_row[sb_row_in_tile - 1];
    }

    sb_completed_in_row = &parse_recon_tile_info_array->sb_recon_completed_in_row[sb_row_in_tile];

    tile_wd_in_sb =
        (AOMMIN(tile_info->tile_col_start_mi[tile_col + 1], dec_handle_ptr->frame_header.mi_cols) +
         ((1 << sb_mi_size_log2) - 1)) >>
        sb_mi_size_log2;

    //tile_wd_in_sb = ( (tile_info->tile_col_start_mi[tile_col + 1] << MI_SIZE_LOG2) ) >>
    //   dec_mod_ctxt->seq_header->sb_size_log2;

    for (uint32_t mi_col = tile_info->tile_col_start_mi[tile_col];
         mi_col < tile_info->tile_col_start_mi[tile_col + 1];
         mi_col += dec_mod_ctxt->seq_header->sb_mi_size)

    {
        int32_t sb_col = (mi_col << MI_SIZE_LOG2) >> dec_mod_ctxt->seq_header->sb_size_log2;

        SBInfo *sb_info = frame_buf->sb_info + (sb_row * master_frame_buf->sb_cols) + sb_col;

        dec_mod_ctxt->cur_coeff[AOM_PLANE_Y] = sb_info->sb_coeff[AOM_PLANE_Y];
        dec_mod_ctxt->cur_coeff[AOM_PLANE_U] = sb_info->sb_coeff[AOM_PLANE_U];
        dec_mod_ctxt->cur_coeff[AOM_PLANE_V] = sb_info->sb_coeff[AOM_PLANE_V];
        /* Top-Right Sync*/
        if (sb_row_in_tile) {
            while (*sb_completed_in_prev_row < MIN((sb_col + 2), tile_wd_in_sb))
                ;
            //Sleep(5); /* ToDo : Change */
        }

        decode_super_block(dec_mod_ctxt, mi_row, mi_col, sb_info);
        *sb_completed_in_row = (uint32_t)(sb_col + 1);
    }

    DecMtFrameData *mt_frame_data = &frame_buf->dec_mt_frame_data;
    int             index         = mi_row / dec_mod_ctxt->seq_header->sb_mi_size;
    mt_frame_data->sb_recon_row_map[(index * tile_info->tile_cols) + tile_col] = 1;
    return status;
}
EbErrorType decode_tile(DecModCtxt *dec_mod_ctxt, TilesInfo *tile_info,
                        DecMtParseReconTileInfo *parse_recon_tile_info_array, int32_t tile_col) {
    EbErrorType status = EB_ErrorNone;

    while (1) {
        EbColorConfig *color_config;
        int32_t        sb_row, mi_row;
        int32_t        sb_row_in_tile = -1;

        int32_t sb_row_tile_start =
            (parse_recon_tile_info_array->tile_info.mi_row_start << MI_SIZE_LOG2) >>
            dec_mod_ctxt->seq_header->sb_size_log2;

        //lock mutex
        eb_block_on_mutex(parse_recon_tile_info_array->tile_sbrow_mutex);

        //pick up a row and increment the sb row counter
        if (parse_recon_tile_info_array->sb_row_to_process !=
            parse_recon_tile_info_array->tile_num_sb_rows) {
            sb_row_in_tile = parse_recon_tile_info_array->sb_row_to_process;
            parse_recon_tile_info_array->sb_row_to_process++;
        }

        //unlock mutex
        eb_release_mutex(parse_recon_tile_info_array->tile_sbrow_mutex);

        //wait for parse
        if (-1 != sb_row_in_tile) {
            volatile int32_t *sb_row_parsed = (volatile int32_t *)&parse_recon_tile_info_array
                                                  ->sb_recon_row_parsed[sb_row_in_tile];
            while (0 == *sb_row_parsed)
                ;

            sb_row = sb_row_in_tile + sb_row_tile_start;

            mi_row = (sb_row << dec_mod_ctxt->seq_header->sb_size_log2) >> MI_SIZE_LOG2;

            color_config = &dec_mod_ctxt->seq_header->color_config;
            cfl_init(&dec_mod_ctxt->cfl_ctx, color_config);

            //update the row started status
            parse_recon_tile_info_array->sb_recon_row_started[sb_row_in_tile] = 1;

            status = decode_tile_row(
                dec_mod_ctxt, tile_info, parse_recon_tile_info_array, tile_col, mi_row, sb_row);
        }

        /*if all sb rows have been picked up for processing then break the while loop */
        if (parse_recon_tile_info_array->sb_row_to_process ==
            parse_recon_tile_info_array->tile_num_sb_rows) {
            break;
        }
    }
    return status;
}

EbErrorType start_decode_tile(EbDecHandle *dec_handle_ptr, DecModCtxt *dec_mod_ctxt,
                              TilesInfo *tiles_info, int32_t tile_num) {
    EbErrorType     status = EB_ErrorNone;
    DecMtFrameData *dec_mt_frame_data =
        &dec_handle_ptr->master_frame_buf.cur_frame_bufs[0].dec_mt_frame_data;

    DecMtParseReconTileInfo *parse_recon_tile_info_array;

    dec_mod_ctxt->frame_header = &dec_handle_ptr->frame_header;
    dec_mod_ctxt->seq_header   = &dec_handle_ptr->seq_header;
    FrameHeader *frame_header  = &dec_handle_ptr->frame_header;
    int32_t      tile_row      = tile_num / tiles_info->tile_cols;
    int32_t      tile_col      = tile_num % tiles_info->tile_cols;
    svt_tile_init(&dec_mod_ctxt->cur_tile_info, frame_header, tile_row, tile_col);

    parse_recon_tile_info_array = &dec_mt_frame_data->parse_recon_tile_info_array[tile_num];
    status = decode_tile(dec_mod_ctxt, tiles_info, parse_recon_tile_info_array, tile_col);

    return status;
}
