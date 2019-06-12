/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

// SUMMARY
//   Contains the Decoder Memory Init functions

/**************************************
 * Includes
 **************************************/
#include <stdlib.h>

#include "EbDefinitions.h"
#include "EbPictureBufferDesc.h"

#include "EbSvtAv1Dec.h"
#include "EbDecHandle.h"
#include "EbDecProcessFrame.h"

#include "EbObuParse.h"

#include "EbDecMemInit.h"
#include "EbDecInverseQuantize.h"

#include "EbDecPicMgr.h"

/*TODO: Remove and harmonize with encoder. Globals prevent harmonization now! */
/*****************************************
 * eb_recon_picture_buffer_desc_ctor
 *  Initializes the Buffer Descriptor's
 *  values that are fixed for the life of
 *  the descriptor.
 *****************************************/
EbErrorType dec_eb_recon_picture_buffer_desc_ctor(
    EbPtr  *object_dbl_ptr,
    EbPtr   object_init_data_ptr)
{
    EbPictureBufferDesc          *picture_buffer_desc_ptr;
    EbPictureBufferDescInitData  *pictureBufferDescInitDataPtr = (EbPictureBufferDescInitData*)object_init_data_ptr;

    uint32_t bytesPerPixel = (pictureBufferDescInitDataPtr->bit_depth == EB_8BIT) ? 1 : 2;

    EB_MALLOC_DEC(EbPictureBufferDesc*, picture_buffer_desc_ptr, sizeof(EbPictureBufferDesc), EB_N_PTR);

    // Allocate the PictureBufferDesc Object
    *object_dbl_ptr = (EbPtr)picture_buffer_desc_ptr;

    // Set the Picture Buffer Static variables
    picture_buffer_desc_ptr->max_width = pictureBufferDescInitDataPtr->max_width;
    picture_buffer_desc_ptr->max_height = pictureBufferDescInitDataPtr->max_height;
    picture_buffer_desc_ptr->width = pictureBufferDescInitDataPtr->max_width;
    picture_buffer_desc_ptr->height = pictureBufferDescInitDataPtr->max_height;
    picture_buffer_desc_ptr->bit_depth = pictureBufferDescInitDataPtr->bit_depth;
    picture_buffer_desc_ptr->color_format = pictureBufferDescInitDataPtr->color_format;
    picture_buffer_desc_ptr->stride_y = pictureBufferDescInitDataPtr->max_width + pictureBufferDescInitDataPtr->left_padding + pictureBufferDescInitDataPtr->right_padding;
    picture_buffer_desc_ptr->stride_cb = picture_buffer_desc_ptr->stride_cr = picture_buffer_desc_ptr->stride_y >> 1;
    picture_buffer_desc_ptr->origin_x = pictureBufferDescInitDataPtr->left_padding;
    picture_buffer_desc_ptr->origin_y = pictureBufferDescInitDataPtr->top_padding;

    picture_buffer_desc_ptr->luma_size = (pictureBufferDescInitDataPtr->max_width + pictureBufferDescInitDataPtr->left_padding + pictureBufferDescInitDataPtr->right_padding) *
        (pictureBufferDescInitDataPtr->max_height + pictureBufferDescInitDataPtr->top_padding + pictureBufferDescInitDataPtr->bot_padding);
    picture_buffer_desc_ptr->chroma_size = picture_buffer_desc_ptr->luma_size >> 2;
    picture_buffer_desc_ptr->packedFlag = EB_FALSE;

    picture_buffer_desc_ptr->stride_bit_inc_y = 0;
    picture_buffer_desc_ptr->stride_bit_inc_cb = 0;
    picture_buffer_desc_ptr->stride_bit_inc_cr = 0;

    // Allocate the Picture Buffers (luma & chroma)
    if (pictureBufferDescInitDataPtr->buffer_enable_mask & PICTURE_BUFFER_DESC_Y_FLAG) {
        EB_ALLIGN_MALLOC_DEC(EbByte, picture_buffer_desc_ptr->buffer_y, picture_buffer_desc_ptr->luma_size * bytesPerPixel, EB_A_PTR);
        memset(picture_buffer_desc_ptr->buffer_y, 0, picture_buffer_desc_ptr->luma_size      * bytesPerPixel);
    }
    else
        picture_buffer_desc_ptr->buffer_y = 0;
    if (pictureBufferDescInitDataPtr->buffer_enable_mask & PICTURE_BUFFER_DESC_Cb_FLAG) {
        EB_ALLIGN_MALLOC_DEC(EbByte, picture_buffer_desc_ptr->buffer_cb, picture_buffer_desc_ptr->chroma_size * bytesPerPixel, EB_A_PTR);
        memset(picture_buffer_desc_ptr->buffer_cb, 0, picture_buffer_desc_ptr->chroma_size      * bytesPerPixel);
    }
    else
        picture_buffer_desc_ptr->buffer_cb = 0;
    if (pictureBufferDescInitDataPtr->buffer_enable_mask & PICTURE_BUFFER_DESC_Cr_FLAG) {
        EB_ALLIGN_MALLOC_DEC(EbByte, picture_buffer_desc_ptr->buffer_cr, picture_buffer_desc_ptr->chroma_size * bytesPerPixel, EB_A_PTR);
        memset(picture_buffer_desc_ptr->buffer_cr, 0, picture_buffer_desc_ptr->chroma_size      * bytesPerPixel);
    }
    else
        picture_buffer_desc_ptr->buffer_cr = 0;
    return EB_ErrorNone;
}

/**********************************
* Master Frame Buf containing all frame level bufs like ModeInfo
for all the frames in parallel
**********************************/

static EbErrorType init_master_frame_ctxt(EbDecHandle  *dec_handle_ptr) {
    EbErrorType return_error = EB_ErrorNone;

    int32_t i, num_sb;
    CurFrameBuf *cur_frame_buf;
    MasterFrameBuf  *master_frame_buf = &dec_handle_ptr->master_frame_buf;
    SeqHeader   *seq_header = &dec_handle_ptr->seq_header;

    ///* 8x8 alignment for various tools like CDEF */
    //int32_t aligned_width   = ALIGN_POWER_OF_TWO(seq_header->max_frame_width, 3);
    //int32_t aligned_height  = ALIGN_POWER_OF_TWO(seq_header->max_frame_height, 3);
    /*int32_t mi_cols = aligned_width >> MI_SIZE_LOG2;
    int32_t mi_rows = aligned_height >> MI_SIZE_LOG2;*/

    int32_t sb_size_log2 = seq_header->sb_size_log2;
    int32_t sb_aligned_width = ALIGN_POWER_OF_TWO(seq_header->max_frame_width,
                                sb_size_log2);
    int32_t sb_aligned_height = ALIGN_POWER_OF_TWO(seq_header->max_frame_height,
                                sb_size_log2);
    int32_t sb_cols = sb_aligned_width >> sb_size_log2;
    int32_t sb_rows = sb_aligned_height >> sb_size_log2;

    num_sb = sb_cols * sb_rows;

    int32_t num_mis_in_sb = (1 << (sb_size_log2 - MI_SIZE_LOG2)) * (1 << (sb_size_log2 - MI_SIZE_LOG2));

    master_frame_buf->num_mis_in_sb = num_mis_in_sb;

    //master_frame_buf->mi_cols = mi_cols;
    //master_frame_buf->mi_rows = mi_rows;

    master_frame_buf->sb_cols = sb_cols;
    master_frame_buf->sb_rows = sb_rows;

    for (i = 0; i < dec_handle_ptr->num_frms_prll; i++) {
        cur_frame_buf = &master_frame_buf->cur_frame_bufs[i];

        /* SuperBlock str allocation at SB level */
        EB_MALLOC_DEC(SBInfo*, cur_frame_buf->sb_info,
            (num_sb * sizeof(SBInfo)), EB_N_PTR);

        /* ModeInfo str allocation at 4x4 level */
        EB_MALLOC_DEC(ModeInfo_t*, cur_frame_buf->mode_info,
                    (num_sb * num_mis_in_sb * sizeof(ModeInfo_t)), EB_N_PTR);

        /* TransformInfo str allocation at 4x4 level */
        EB_MALLOC_DEC(TransformInfo_t*, cur_frame_buf->trans_info[AOM_PLANE_Y],
            (num_sb * num_mis_in_sb * sizeof(TransformInfo_t)), EB_N_PTR);

        /* Coeff buf (1D compact) allocation for entire frame
            TODO: Should reduce this to save memory and
            dynammically allocate if needed */
            /*TODO : Change to macro */
            /* (16+1) : 1 for Length and 16 for all coeffs in 4x4 */
        EB_MALLOC_DEC(int32_t*, cur_frame_buf->coeff[AOM_PLANE_Y],
            (num_sb * num_mis_in_sb * sizeof(int32_t) * (16 + 1)), EB_N_PTR);

        if(seq_header->color_config.subsampling_x == 1 &&
           seq_header->color_config.subsampling_y == 1)
        {
            /*TODO : Change to macro */
            EB_MALLOC_DEC(TransformInfo_t*, cur_frame_buf->trans_info[AOM_PLANE_U],
                (num_sb * num_mis_in_sb * sizeof(TransformInfo_t) * 2), EB_N_PTR);

            /* Coeff buf (1D compact) allocation for entire frame
            TODO: Should reduce this to save memory and
            dynammically allocate if needed */
            /*TODO : Change to macro */
            /* (16+1) : 1 for Length and 16 for all coeffs in 4x4 */
            EB_MALLOC_DEC(int32_t*, cur_frame_buf->coeff[AOM_PLANE_U],
            (num_sb * num_mis_in_sb * sizeof(int32_t) * (16+1) >> 2), EB_N_PTR);
            EB_MALLOC_DEC(int32_t*, cur_frame_buf->coeff[AOM_PLANE_V],
            (num_sb * num_mis_in_sb * sizeof(int32_t) * (16+1) >> 2), EB_N_PTR);
        }
        else
            assert(0);

        /* delta_q allocation at SB level */
        EB_MALLOC_DEC(int32_t*, cur_frame_buf->delta_q,
            (num_sb * sizeof(int32_t)), EB_N_PTR);

        /* cdef_strength allocation at SB level */
        EB_MALLOC_DEC(int8_t*, cur_frame_buf->cdef_strength,
            (num_sb * (seq_header->use_128x128_superblock ? 4 : 1) *
            sizeof(int8_t)), EB_N_PTR);

        /* delta_lf allocation at SB level */
        EB_MALLOC_DEC(int32_t*, cur_frame_buf->delta_lf,
            (num_sb * sizeof(int32_t)), EB_N_PTR);

        /* tile map allocation at SB level */
        EB_MALLOC_DEC(uint8_t*, cur_frame_buf->tile_map_sb,
            (num_sb * sizeof(uint8_t)), EB_N_PTR);
    }
#if FRAME_MI_MAP
    FrameMiMap *frame_mi_map = &master_frame_buf->frame_mi_map;
    frame_mi_map->sb_cols = sb_cols;
    frame_mi_map->sb_rows = sb_rows;
    frame_mi_map->mi_cols_algnsb = sb_cols * (1 << (sb_size_log2 - MI_SIZE_LOG2));
    frame_mi_map->mi_rows_algnsb = sb_cols * (1 << (sb_size_log2 - MI_SIZE_LOG2));
    /* SBInfo pointers for entire frame */
    EB_MALLOC_DEC(SBInfo**, frame_mi_map->pps_sb_info,
        sb_rows * sb_cols * sizeof(SBInfo *), EB_N_PTR);
    /* ModeInfo offset wrt it's SB start for entire frame at 4x4 lvl */
    EB_MALLOC_DEC(uint16_t*, frame_mi_map->p_mi_offset, frame_mi_map->
    mi_rows_algnsb * frame_mi_map->mi_cols_algnsb * sizeof(uint16_t), EB_N_PTR);
    frame_mi_map->sb_size_log2 = sb_size_log2;
#else
    /* Top SB 4x4 row MI map */
    EB_MALLOC_DEC(int16_t*, frame_mi_map->top_sbrow_mi_map,
        (sb_cols * (1 << (sb_size_log2 - MI_SIZE_LOG2)) * sizeof(int16_t)), EB_N_PTR);
#endif
    frame_mi_map->num_mis_in_sb_wd = (1 << (sb_size_log2 - MI_SIZE_LOG2));
#if 0
    /* TODO: Recon Pic Buf. Should be generalized! */
    EbPictureBufferDescInitData input_picture_buffer_desc_init_data;
    // Init Picture Init data
    input_picture_buffer_desc_init_data.max_width = seq_header->max_frame_width;
    input_picture_buffer_desc_init_data.max_height = seq_header->max_frame_height;
    input_picture_buffer_desc_init_data.bit_depth  = (EbBitDepthEnum)seq_header->color_config.bit_depth;

    input_picture_buffer_desc_init_data.color_format    = dec_handle_ptr->
                                                dec_config.max_color_format;
    input_picture_buffer_desc_init_data.buffer_enable_mask =
                                                PICTURE_BUFFER_DESC_FULL_MASK;

    input_picture_buffer_desc_init_data.left_padding = PAD_VALUE;
    input_picture_buffer_desc_init_data.right_padding = PAD_VALUE;
    input_picture_buffer_desc_init_data.top_padding = PAD_VALUE;
    input_picture_buffer_desc_init_data.bot_padding = PAD_VALUE;

    input_picture_buffer_desc_init_data.split_mode = EB_FALSE;

    return_error = dec_eb_recon_picture_buffer_desc_ctor(
        (EbPtr*) &(dec_handle_ptr->recon_picture_buf[0]),
        (EbPtr)&input_picture_buffer_desc_init_data);
#endif

    master_frame_buf->tpl_mvs = NULL;
    master_frame_buf->tpl_mvs_size = 0;

    return return_error;
}

/*TODO: Move to module files */
static EbErrorType init_parse_context (EbDecHandle  *dec_handle_ptr) {
    EbErrorType return_error = EB_ErrorNone;

    EB_MALLOC_DEC(void *, dec_handle_ptr->pv_parse_ctxt, sizeof(ParseCtxt), EB_N_PTR);

    ParseCtxt *parse_ctx = (ParseCtxt*)dec_handle_ptr->pv_parse_ctxt;

    SeqHeader   *seq_header = &dec_handle_ptr->seq_header;

    int32_t num_mi_row, num_mi_col;
    //int32_t num_mi_frame;

    int32_t num_4x4_neigh_sb    = seq_header->sb_mi_size;
    int32_t sb_size_log2        = seq_header->sb_size_log2;
    int32_t sb_aligned_width    = ALIGN_POWER_OF_TWO(seq_header->max_frame_width,
                                    sb_size_log2);
    /*int32_t sb_aligned_height   = ALIGN_POWER_OF_TWO(seq_header->max_frame_height,
                                    sb_size_log2);*/
    int32_t sb_cols             = sb_aligned_width >> sb_size_log2;
    //int32_t sb_rows             = sb_aligned_height >> sb_size_log2;
    int8_t num_planes           = seq_header->color_config.mono_chrome ? 1 : MAX_MB_PLANE;

    ParseNbr4x4Ctxt *neigh_ctx = &parse_ctx->parse_nbr4x4_ctxt;

    num_mi_col = sb_cols * num_4x4_neigh_sb;
    num_mi_row = num_4x4_neigh_sb;
    //num_mi_frame = sb_cols * sb_rows * num_4x4_neigh_sb;

    EB_MALLOC_DEC(uint8_t*, neigh_ctx->above_tx_wd, num_mi_col * sizeof(uint8_t), EB_N_PTR);
    EB_MALLOC_DEC(uint8_t*, neigh_ctx->left_tx_ht, num_mi_row * sizeof(uint8_t), EB_N_PTR);

    EB_MALLOC_DEC(uint8_t*, neigh_ctx->above_part_wd, num_mi_col * sizeof(uint8_t), EB_N_PTR);
    EB_MALLOC_DEC(uint8_t*, neigh_ctx->left_part_ht, num_mi_row * sizeof(uint8_t), EB_N_PTR);
    /* TODO : Optimize the size for Chroma */
    for (int i = 0; i < num_planes; i++) {
        EB_MALLOC_DEC(uint8_t*, neigh_ctx->above_dc_ctx[i], num_mi_col * sizeof(uint8_t), EB_N_PTR);
        EB_MALLOC_DEC(uint8_t*, neigh_ctx->left_dc_ctx[i], num_mi_row * sizeof(uint8_t), EB_N_PTR);

        EB_MALLOC_DEC(uint8_t*, neigh_ctx->above_level_ctx[i], num_mi_col * sizeof(uint8_t), EB_N_PTR);
        EB_MALLOC_DEC(uint8_t*, neigh_ctx->left_level_ctx[i], num_mi_row * sizeof(uint8_t), EB_N_PTR);
    }

    EB_MALLOC_DEC(uint8_t*, neigh_ctx->above_seg_pred_ctx, num_mi_col * sizeof(uint8_t), EB_N_PTR);
    EB_MALLOC_DEC(uint8_t*, neigh_ctx->left_seg_pred_ctx, num_mi_row * sizeof(uint8_t), EB_N_PTR);

    return return_error;
}

/*TODO: Move to module files */
static EbErrorType init_dec_mod_ctxt(EbDecHandle  *dec_handle_ptr)
{
    EbErrorType return_error = EB_ErrorNone;

    EB_MALLOC_DEC(void *, dec_handle_ptr->pv_dec_mod_ctxt, sizeof(DecModCtxt), EB_N_PTR);

    DecModCtxt *dec_mod_ctxt = (DecModCtxt*)dec_handle_ptr->pv_dec_mod_ctxt;

    dec_mod_ctxt->dec_handle_ptr = (void *)dec_handle_ptr;

    int32_t sb_size_log2 = dec_handle_ptr->seq_header.sb_size_log2;
    EB_MALLOC_DEC(int32_t*, dec_mod_ctxt->sb_iquant_ptr, (1 << sb_size_log2) *
                  (1 << sb_size_log2) * sizeof(int32_t), EB_N_PTR);

    av1_inverse_qm_init(dec_handle_ptr);

    return return_error;
}

EbErrorType dec_mem_init(EbDecHandle  *dec_handle_ptr) {
    EbErrorType return_error = EB_ErrorNone;

    if (0 == dec_handle_ptr->seq_header_done)
        return EB_ErrorNone;

    /* init module ctxts */
    return_error |= dec_pic_mgr_init((EbDecPicMgr **)&dec_handle_ptr->pv_pic_mgr);

    return_error |= init_parse_context(dec_handle_ptr);

    return_error |= init_dec_mod_ctxt(dec_handle_ptr);

    /* init frame buffers */
    return_error |= init_master_frame_ctxt(dec_handle_ptr);

    /* Initialize the references to NULL */
    for (int i = 0; i < REF_FRAMES; i++) {
        dec_handle_ptr->ref_frame_map[i] = NULL;
        dec_handle_ptr->next_ref_frame_map[i] = NULL;
        dec_handle_ptr->remapped_ref_idx[i] = INVALID_IDX;
    }
    dec_handle_ptr->cur_pic_buf[0] = NULL;

    dec_handle_ptr->mem_init_done = 1;

    return return_error;
}
