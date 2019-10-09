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
#include "EbPictureBufferDesc.h"

#include "EbSvtAv1Dec.h"
#include "EbDecHandle.h"

#include "EbDecInverseQuantize.h"
#include "EbDecProcessFrame.h"
#include "EbDecProcessBlock.h"
#include "EbDecNbr.h"

/* decode partition */
static void decode_partition(DecModCtxt *dec_mod_ctxt,
                             uint32_t mi_row, uint32_t mi_col,
                             BlockSize bsize, SBInfo *sb_info)
{
    EbDecHandle *dec_handle = (EbDecHandle *)dec_mod_ctxt->dec_handle_ptr;
    BlockSize  subsize;
    PartitionType partition;

    uint8_t num4x4 = mi_size_wide[bsize];
    uint32_t half_block_4x4 =(uint32_t)num4x4 >> 1;
    uint32_t quarter_block_4x4 = half_block_4x4 >> 1;

    uint32_t has_rows = (mi_row + half_block_4x4) < dec_handle->frame_header.mi_rows;
    uint32_t has_cols = (mi_col + half_block_4x4) < dec_handle->frame_header.mi_cols;

    if (mi_row >= dec_handle->frame_header.mi_rows ||
        mi_col >= dec_handle->frame_header.mi_cols) return;

    partition = get_partition(dec_mod_ctxt, &dec_handle->frame_header,
                              mi_row, mi_col, sb_info, bsize);

    subsize = Partition_Subsize[partition][bsize];
    BlockSize splitSize = Partition_Subsize[PARTITION_SPLIT][bsize];

#define DECODE_BLOCK(db_r, db_c, db_subsize)                \
decode_block(dec_mod_ctxt, db_r, db_c, db_subsize,          \
    dec_mod_ctxt->cur_tile_info, sb_info)

#define DECODE_PARTITION(db_r, db_c, db_subsize)            \
decode_partition(dec_mod_ctxt, (db_r), (db_c),              \
                   (db_subsize), sb_info)

    switch ((int)partition) {
        case PARTITION_NONE:
            DECODE_BLOCK(mi_row, mi_col, subsize);
            break;
        case PARTITION_HORZ:
            DECODE_BLOCK(mi_row, mi_col, subsize);
            if (has_rows)
                DECODE_BLOCK(mi_row + half_block_4x4,
                             mi_col, subsize);
            break;
        case PARTITION_VERT:
            DECODE_BLOCK(mi_row, mi_col, subsize);
            if (has_cols)
                DECODE_BLOCK(mi_row, mi_col + half_block_4x4, subsize);
            break;
        case PARTITION_SPLIT:
            DECODE_PARTITION(mi_row, mi_col, subsize);
            DECODE_PARTITION(mi_row, mi_col + half_block_4x4, subsize);
            DECODE_PARTITION(mi_row + half_block_4x4, mi_col, subsize);
            DECODE_PARTITION(mi_row + half_block_4x4, mi_col + half_block_4x4, subsize);
            break;
        case PARTITION_HORZ_A:
            DECODE_BLOCK(mi_row, mi_col, splitSize);
            DECODE_BLOCK(mi_row, mi_col + half_block_4x4, splitSize);
            DECODE_BLOCK(mi_row + half_block_4x4, mi_col, subsize);
            break;
        case PARTITION_HORZ_B:
            DECODE_BLOCK(mi_row, mi_col, subsize);
            DECODE_BLOCK(mi_row + half_block_4x4, mi_col, splitSize);
            DECODE_BLOCK(mi_row + half_block_4x4, mi_col + half_block_4x4, splitSize);
            break;
        case PARTITION_VERT_A:
            DECODE_BLOCK(mi_row, mi_col, splitSize);
            DECODE_BLOCK(mi_row + half_block_4x4, mi_col, splitSize);
            DECODE_BLOCK(mi_row, mi_col + half_block_4x4, subsize);
            break;
        case PARTITION_VERT_B:
            DECODE_BLOCK(mi_row, mi_col, subsize);
            DECODE_BLOCK(mi_row, mi_col + half_block_4x4, splitSize);
            DECODE_BLOCK(mi_row + half_block_4x4, mi_col + half_block_4x4, splitSize);
            break;
        case PARTITION_HORZ_4:
            for (int i = 0; i < 4; ++i) {
                uint32_t this_mi_row = mi_row +  (i * quarter_block_4x4);
                if (i > 0 && this_mi_row >= dec_handle->frame_header.mi_rows) break;
                DECODE_BLOCK(this_mi_row, mi_col, subsize);
                }
            break;
        case PARTITION_VERT_4:
            for (int i = 0; i < 4; ++i) {
                uint32_t this_mi_col = mi_col + (i * quarter_block_4x4);
                if (i > 0 && this_mi_col >= dec_handle->frame_header.mi_cols) break;
                DECODE_BLOCK(mi_row, this_mi_col, subsize);
                }
            break;
        default: assert(0 && "Invalid partition type");
    }
}

// decoding of the superblock
void decode_super_block(DecModCtxt *dec_mod_ctxt,
                        uint32_t mi_row, uint32_t mi_col,
                        SBInfo *sb_info)
{
    EbDecHandle *dec_handle = (EbDecHandle *)dec_mod_ctxt->dec_handle_ptr;
    SeqHeader *seq = &dec_handle->seq_header;

    /* Pointer updates */
    bool do_memset = true;
    int left_available = (mi_col > (uint32_t)dec_mod_ctxt->cur_tile_info->mi_col_start);
    if (left_available) {
        BlockModeInfo *left_mode = get_left_mode_info(dec_handle, mi_row, mi_col, sb_info);
        if (left_mode->skip && left_mode->sb_type == seq->sb_size) {
            do_memset = false;
        }
    }

    if (do_memset) {
        EbColorConfig *color_config = &seq->color_config;
        int32_t sb_size_log2 = seq->sb_size_log2;

        int32_t y_size = (1 << sb_size_log2) * (1 << sb_size_log2);
        int32_t iq_size = y_size +
            (color_config->subsampling_x ? y_size >> 2 : y_size) +
            (color_config->subsampling_y ? y_size >> 2 : y_size);

        memset(dec_mod_ctxt->sb_iquant_ptr, 0, iq_size *
            sizeof(*dec_mod_ctxt->sb_iquant_ptr));
    }

    dec_mod_ctxt->iquant_cur_ptr = dec_mod_ctxt->sb_iquant_ptr;

    /* SB level dequant update */
    update_dequant(dec_handle, sb_info);

    /* Decode partition */
    decode_partition(dec_mod_ctxt, mi_row, mi_col,
                     dec_handle->seq_header.sb_size, sb_info);
}
