// clang-format off
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
#include "EbDecParseFrame.h"

#include "EbDecMemInit.h"
#include "EbDecInverseQuantize.h"

#include "EbDecPicMgr.h"
#include "EbDecLF.h"

#include "EbUtility.h"

/*TODO: Remove and harmonize with encoder. Globals prevent harmonization now! */
/*****************************************
 * eb_recon_picture_buffer_desc_ctor
 *  Initializes the Buffer Descriptor's
 *  values that are fixed for the life of
 *  the descriptor.
 *****************************************/
EbErrorType dec_eb_recon_picture_buffer_desc_ctor(
    EbPtr  *object_dbl_ptr,
    EbPtr   object_init_data_ptr,
    EbBool is_16bit_pipeline /* can be removed as an extra argument once
                                EbPictureBufferDescInitData adds the support for this */
)
{
    EbPictureBufferDesc          *picture_buffer_desc_ptr;
    EbPictureBufferDescInitData  *picture_buffer_desc_init_data_ptr = (EbPictureBufferDescInitData*)object_init_data_ptr;

    EB_MALLOC_DEC(EbPictureBufferDesc*, picture_buffer_desc_ptr, sizeof(EbPictureBufferDesc), EB_N_PTR);

    uint32_t bytes_per_pixel = (picture_buffer_desc_init_data_ptr->bit_depth > EB_8BIT ||
        is_16bit_pipeline) ? 2 : 1;
    picture_buffer_desc_ptr->is_16bit_pipeline = is_16bit_pipeline;

    // Allocate the PictureBufferDesc Object
    *object_dbl_ptr = (EbPtr)picture_buffer_desc_ptr;

    // Set the Picture Buffer Static variables
    picture_buffer_desc_ptr->max_width = picture_buffer_desc_init_data_ptr->max_width;
    picture_buffer_desc_ptr->max_height = picture_buffer_desc_init_data_ptr->max_height;
    picture_buffer_desc_ptr->width = picture_buffer_desc_init_data_ptr->max_width;
    picture_buffer_desc_ptr->height = picture_buffer_desc_init_data_ptr->max_height;
    picture_buffer_desc_ptr->bit_depth = picture_buffer_desc_init_data_ptr->bit_depth;
    picture_buffer_desc_ptr->color_format = picture_buffer_desc_init_data_ptr->color_format;
    picture_buffer_desc_ptr->stride_y = picture_buffer_desc_init_data_ptr->max_width +
        picture_buffer_desc_init_data_ptr->left_padding + picture_buffer_desc_init_data_ptr->right_padding;
    uint32_t height_y = (picture_buffer_desc_init_data_ptr->max_height + picture_buffer_desc_init_data_ptr->top_padding
        + picture_buffer_desc_init_data_ptr->bot_padding);

    picture_buffer_desc_ptr->origin_x = picture_buffer_desc_init_data_ptr->left_padding;
    picture_buffer_desc_ptr->origin_y = picture_buffer_desc_init_data_ptr->top_padding;
    picture_buffer_desc_ptr->origin_bot_y = picture_buffer_desc_init_data_ptr->bot_padding;

    picture_buffer_desc_ptr->luma_size = (picture_buffer_desc_ptr->stride_y) * height_y;

    uint32_t stride_c = 0, height_c = 0;
    if (picture_buffer_desc_ptr->color_format == EB_YUV420) {// 420
        stride_c = (picture_buffer_desc_ptr->stride_y + 1) >> 1;
        height_c = (height_y + 1) >> 1;
    }
    else if (picture_buffer_desc_ptr->color_format == EB_YUV422) {// 422
        stride_c = (picture_buffer_desc_ptr->stride_y + 1) >> 1;
        height_c = height_y;
    }
    else if (picture_buffer_desc_ptr->color_format == EB_YUV444) {// 444
        stride_c = picture_buffer_desc_ptr->stride_y;
        height_c = height_y;
    }
    picture_buffer_desc_ptr->stride_cb = picture_buffer_desc_ptr->stride_cr = stride_c;
    picture_buffer_desc_ptr->chroma_size = (stride_c * height_c);

    picture_buffer_desc_ptr->packed_flag = EB_FALSE;

    picture_buffer_desc_ptr->stride_bit_inc_y = 0;
    picture_buffer_desc_ptr->stride_bit_inc_cb = 0;
    picture_buffer_desc_ptr->stride_bit_inc_cr = 0;

    // Allocate the Picture Buffers (luma & chroma)
    if (picture_buffer_desc_init_data_ptr->buffer_enable_mask & PICTURE_BUFFER_DESC_Y_FLAG) {
        EB_ALLIGN_MALLOC_DEC(EbByte, picture_buffer_desc_ptr->buffer_y, picture_buffer_desc_ptr->luma_size * bytes_per_pixel, EB_A_PTR);
        memset(picture_buffer_desc_ptr->buffer_y, 0, picture_buffer_desc_ptr->luma_size      * bytes_per_pixel);
    }
    else
        picture_buffer_desc_ptr->buffer_y = 0;
    if (picture_buffer_desc_init_data_ptr->buffer_enable_mask & PICTURE_BUFFER_DESC_Cb_FLAG) {
        EB_ALLIGN_MALLOC_DEC(EbByte, picture_buffer_desc_ptr->buffer_cb, picture_buffer_desc_ptr->chroma_size * bytes_per_pixel, EB_A_PTR);
        memset(picture_buffer_desc_ptr->buffer_cb, 0, picture_buffer_desc_ptr->chroma_size      * bytes_per_pixel);
    }
    else
        picture_buffer_desc_ptr->buffer_cb = 0;
    if (picture_buffer_desc_init_data_ptr->buffer_enable_mask & PICTURE_BUFFER_DESC_Cr_FLAG) {
        EB_ALLIGN_MALLOC_DEC(EbByte, picture_buffer_desc_ptr->buffer_cr, picture_buffer_desc_ptr->chroma_size * bytes_per_pixel, EB_A_PTR);
        memset(picture_buffer_desc_ptr->buffer_cr, 0, picture_buffer_desc_ptr->chroma_size      * bytes_per_pixel);
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

    EbBool is_st = dec_handle_ptr->dec_config.threads == 1 ? EB_TRUE : EB_FALSE;
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
        EB_MALLOC_DEC(BlockModeInfo*, cur_frame_buf->mode_info,
                    (num_sb * num_mis_in_sb * sizeof(BlockModeInfo)), EB_N_PTR);

        /* TransformInfo str allocation at 4x4 level
           TO-DO optimize memory based on the chroma subsampling.*/
        EB_MALLOC_DEC(TransformInfo_t*, cur_frame_buf->trans_info[AOM_PLANE_Y],
            (num_sb * num_mis_in_sb * sizeof(TransformInfo_t)), EB_N_PTR);

        /* Coeff buf (1D compact) allocation for entire frame
            TODO: Should reduce this to save memory and
            dynammically allocate if needed */
            /*TODO : Change to macro */
            /* (16+1) : 1 for Length and 16 for all coeffs in 4x4 */
        if (is_st) {
            /*Size of coeff buf reduced to sb_sizesss*/
            EB_MALLOC_DEC(int32_t*, cur_frame_buf->coeff[AOM_PLANE_Y],
                (num_mis_in_sb * sizeof(int32_t) * (16 + 1)), EB_N_PTR);
        }
        else {
            EB_MALLOC_DEC(int32_t*, cur_frame_buf->coeff[AOM_PLANE_Y],
                (num_sb * num_mis_in_sb * sizeof(int32_t) * (16 + 1)), EB_N_PTR);
        }

        /*TODO : Change to macro */
        EB_MALLOC_DEC(TransformInfo_t*, cur_frame_buf->trans_info[AOM_PLANE_U],
            (num_sb * num_mis_in_sb * sizeof(TransformInfo_t) * 2), EB_N_PTR);

        /* Coeff buf (1D compact) allocation for entire frame
        TODO: Should reduce this to save memory and
        dynammically allocate if needed */
        /*TODO : Change to macro */
        /* (16+1) : 1 for Length and 16 for all coeffs in 4x4 */
        if (seq_header->color_config.subsampling_x == 1 &&
            seq_header->color_config.subsampling_y == 1) // 420
        {
            if (is_st) {
                EB_MALLOC_DEC(int32_t*, cur_frame_buf->coeff[AOM_PLANE_U],
                    (num_mis_in_sb * sizeof(int32_t) * (16 + 1) >> 2), EB_N_PTR);
                EB_MALLOC_DEC(int32_t*, cur_frame_buf->coeff[AOM_PLANE_V],
                    (num_mis_in_sb * sizeof(int32_t) * (16 + 1) >> 2), EB_N_PTR);
            }
            else {
                EB_MALLOC_DEC(int32_t*,
                    cur_frame_buf->coeff[AOM_PLANE_U],
                    (num_sb * num_mis_in_sb * sizeof(int32_t) * (16 + 1) >> 2),
                    EB_N_PTR);
                EB_MALLOC_DEC(int32_t*,
                    cur_frame_buf->coeff[AOM_PLANE_V],
                    (num_sb * num_mis_in_sb * sizeof(int32_t) * (16 + 1) >> 2),
                    EB_N_PTR);
            }
        }
        else if (seq_header->color_config.subsampling_x == 1 &&
                 seq_header->color_config.subsampling_y == 0) // 422
        {
            if (is_st) {
                EB_MALLOC_DEC(int32_t*, cur_frame_buf->coeff[AOM_PLANE_U],
                    (num_mis_in_sb * sizeof(int32_t) * (16 + 1) >> 1), EB_N_PTR);
                EB_MALLOC_DEC(int32_t*, cur_frame_buf->coeff[AOM_PLANE_V],
                    (num_mis_in_sb * sizeof(int32_t) * (16 + 1) >> 1), EB_N_PTR);
            }
            else {
                EB_MALLOC_DEC(int32_t*,
                    cur_frame_buf->coeff[AOM_PLANE_U],
                    (num_sb * num_mis_in_sb * sizeof(int32_t) * (16 + 1) >> 1),
                    EB_N_PTR);
                EB_MALLOC_DEC(int32_t*,
                    cur_frame_buf->coeff[AOM_PLANE_V],
                    (num_sb * num_mis_in_sb * sizeof(int32_t) * (16 + 1) >> 1),
                    EB_N_PTR);
            }
        }
        else if (seq_header->color_config.subsampling_x == 0 &&
                 seq_header->color_config.subsampling_y == 0) // 444
        {
            if (is_st) {
                EB_MALLOC_DEC(int32_t*, cur_frame_buf->coeff[AOM_PLANE_U],
                    (num_mis_in_sb * sizeof(int32_t) * (16 + 1)), EB_N_PTR);
                EB_MALLOC_DEC(int32_t*, cur_frame_buf->coeff[AOM_PLANE_V],
                    (num_mis_in_sb * sizeof(int32_t) * (16 + 1)), EB_N_PTR);
            }
            else {
                EB_MALLOC_DEC(int32_t*,
                    cur_frame_buf->coeff[AOM_PLANE_U],
                    (num_sb * num_mis_in_sb * sizeof(int32_t) * (16 + 1)),
                    EB_N_PTR);
                EB_MALLOC_DEC(int32_t*,
                    cur_frame_buf->coeff[AOM_PLANE_V],
                    (num_sb * num_mis_in_sb * sizeof(int32_t) * (16 + 1)),
                    EB_N_PTR);
            }
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
        memset(cur_frame_buf->cdef_strength, -1, (num_sb *
            (seq_header->use_128x128_superblock ? 4 : 1) *
            sizeof(int8_t)));

        /* delta_lf allocation at SB level */
        EB_MALLOC_DEC(int32_t*, cur_frame_buf->delta_lf,
            (num_sb * FRAME_LF_COUNT * sizeof(int32_t)), EB_N_PTR);

        /* tile map allocation at SB level */
        EB_MALLOC_DEC(uint8_t*, cur_frame_buf->tile_map_sb,
            (num_sb * sizeof(uint8_t)), EB_N_PTR);

        // Allocating lr_unit based on SB_SIZE as worst case memory.
        // rest_unit_size cannot be less than SB_size.
        // if rest_unit_size > SB_size then holes are introduced in-between and
        // accessing will skip few SB in-between.
        // if rest_unit_size == SB_size then it's straight forward to access
        // every SB level loop restoration filter value.
        LrCtxt *lr_ctxt = (LrCtxt *)dec_handle_ptr->pv_lr_ctxt;
        for (int32_t plane = 0; plane <= AOM_PLANE_V; plane++) {
            EB_MALLOC_DEC(RestorationUnitInfo *, cur_frame_buf->lr_unit[plane],
                (num_sb * sizeof(RestorationUnitInfo)), EB_N_PTR);
            lr_ctxt->lr_unit[plane] = cur_frame_buf->lr_unit[plane];
            lr_ctxt->lr_stride[plane] = sb_cols;
        }
    }
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
    frame_mi_map->num_mis_in_sb_wd = (1 << (sb_size_log2 - MI_SIZE_LOG2));


    master_frame_buf->tpl_mvs = NULL;
    master_frame_buf->tpl_mvs_size = 0;

    return return_error;
}

/*TODO: Move to module files */
static EbErrorType init_parse_context (EbDecHandle  *dec_handle_ptr) {
    EbErrorType return_error = EB_ErrorNone;

    EB_MALLOC_DEC(void *, dec_handle_ptr->pv_master_parse_ctxt,
        sizeof(MasterParseCtxt), EB_N_PTR);
    MasterParseCtxt *master_parse_ctx =
        (MasterParseCtxt*)dec_handle_ptr->pv_master_parse_ctxt;

    master_parse_ctx->context_count = 0;
    master_parse_ctx->tile_parse_ctxt = NULL;

    master_parse_ctx->parse_above_nbr4x4_ctxt = NULL;
    master_parse_ctx->parse_left_nbr4x4_ctxt = NULL;

    master_parse_ctx->num_tiles = 0;
    master_parse_ctx->parse_tile_data = NULL;

    return return_error;
}

/*TODO: Move to module files */
EbErrorType init_dec_mod_ctxt(EbDecHandle  *dec_handle_ptr,
    void **pp_dec_mod_ctxt)
{
    EbErrorType return_error = EB_ErrorNone;
    SeqHeader *seq_header = &dec_handle_ptr->seq_header;
    EbColorConfig *color_config = &seq_header->color_config;

    EB_MALLOC_DEC(void *, *pp_dec_mod_ctxt, sizeof(DecModCtxt), EB_N_PTR);

    DecModCtxt *p_dec_mod_ctxt = (DecModCtxt*)*pp_dec_mod_ctxt;
    p_dec_mod_ctxt->dec_handle_ptr  = (void *)dec_handle_ptr;
    p_dec_mod_ctxt->seq_header      = &dec_handle_ptr->seq_header;
    p_dec_mod_ctxt->frame_header    = &dec_handle_ptr->frame_header;

    int32_t sb_size_log2 = seq_header->sb_size_log2;

    int32_t y_size = (1 << sb_size_log2) * (1 << sb_size_log2);
    int32_t iq_size = y_size +
        (color_config->subsampling_x ? y_size >> 2 : y_size) +
        (color_config->subsampling_y ? y_size >> 2 : y_size);

    EB_MALLOC_DEC(int32_t*, p_dec_mod_ctxt->sb_iquant_ptr,
        iq_size * sizeof(int32_t), EB_N_PTR);
    av1_inverse_qm_init(p_dec_mod_ctxt, seq_header);

    EbColorConfig *cc = &dec_handle_ptr->seq_header.color_config;
    uint32_t use_highbd = (cc->bit_depth > EB_8BIT ||
        dec_handle_ptr->is_16bit_pipeline);
    int32_t sb_size = 1 << sb_size_log2;
    uint16_t *hbd_mc_buf[2];
    for (int ref = 0; ref < 2; ref++) {

        //EB_MALLOC_DEC(uint8_t**, part_info->mc_buf[ref],
        //    sizeof(uint8_t*), EB_N_PTR);
        if (use_highbd) {
            EB_MALLOC_DEC(uint16_t*, hbd_mc_buf[ref],
                ((2 * sb_size) + (AOM_INTERP_EXTEND * 2))*
                ((2 * sb_size) + (AOM_INTERP_EXTEND * 2))*
                sizeof(uint16_t), EB_N_PTR);
            p_dec_mod_ctxt->mc_buf[ref] = (uint8_t *)hbd_mc_buf[ref];
        }
        else {
            EB_MALLOC_DEC(uint8_t*, p_dec_mod_ctxt->mc_buf[ref],
                ((2 * sb_size) + (AOM_INTERP_EXTEND * 2))*
                ((2 * sb_size) + (AOM_INTERP_EXTEND * 2))*
                sizeof(uint8_t), EB_N_PTR);
        }
    }

    return return_error;
}

/*mem init function for LF params*/
static EbErrorType init_lf_ctxt(EbDecHandle  *dec_handle_ptr) {

    EbErrorType return_error = EB_ErrorNone;

    SeqHeader *seq_header = &dec_handle_ptr->seq_header;
    /*Boundary checking of mi_row & mi_col are not done while populating,
    so more memory is allocated by alligning to sb_size */
    int32_t aligned_width   = ALIGN_POWER_OF_TWO(seq_header->max_frame_width,
        MAX_SB_SIZE_LOG2);
    int32_t aligned_height  = ALIGN_POWER_OF_TWO(seq_header->max_frame_height,
        MAX_SB_SIZE_LOG2);
    int32_t mi_cols = aligned_width >> MI_SIZE_LOG2;
    int32_t mi_rows = aligned_height >> MI_SIZE_LOG2;

    EB_MALLOC_DEC(void *, dec_handle_ptr->pv_lf_ctxt, sizeof(LfCtxt), EB_N_PTR);

    LfCtxt *lf_ctxt = (LfCtxt *)dec_handle_ptr->pv_lf_ctxt;
    /*Mem allocation for luma parmas 4x4 unit*/
    EB_MALLOC_DEC(TxSize *, lf_ctxt->tx_size_l,
        mi_rows * mi_cols * sizeof(TxSize), EB_N_PTR);

    /*Allocation of chroma params at 4x4 luma unit, can be optimized */
    EB_MALLOC_DEC(TxSize *, lf_ctxt->tx_size_uv,
        mi_rows * mi_cols * sizeof(TxSize), EB_N_PTR);

    return return_error;
}

static EbErrorType init_lr_ctxt(EbDecHandle  *dec_handle_ptr)
{
    EbErrorType return_error = EB_ErrorNone;
    EB_MALLOC_DEC(void *, dec_handle_ptr->pv_lr_ctxt, sizeof(LrCtxt), EB_N_PTR);

    LrCtxt *lr_ctxt = (LrCtxt*)dec_handle_ptr->pv_lr_ctxt;
    lr_ctxt->dec_handle_ptr = (void *)dec_handle_ptr;

    int32_t sb_size_h = block_size_high[dec_handle_ptr->seq_header.sb_size];
    uint32_t picture_height_in_sb = (dec_handle_ptr->seq_header.
        max_frame_height + sb_size_h - 1) / sb_size_h;
    EbBool is_mt = dec_handle_ptr->dec_config.threads > 1;

    picture_height_in_sb = (is_mt == 0) ? 1 : picture_height_in_sb;
    const int32_t num_planes = av1_num_planes(&dec_handle_ptr->seq_header.
        color_config);

    uint32_t num_instances = MIN(picture_height_in_sb,
        dec_handle_ptr->dec_config.threads);

    lr_ctxt->is_thread_min = EB_FALSE;
    if (num_instances == dec_handle_ptr->dec_config.threads)
        lr_ctxt->is_thread_min = EB_TRUE;
    EB_MALLOC_DEC(RestorationLineBuffers ***, lr_ctxt->rlbs,
        num_instances * sizeof(RestorationLineBuffers**), EB_N_PTR)
    EB_MALLOC_DEC(int32_t **, lr_ctxt->rst_tmpbuf,
         num_instances * RESTORATION_TMPBUF_SIZE, EB_N_PTR)
    for (uint32_t i = 0; i < num_instances; i++) {
        RestorationLineBuffers **p_rlbs;
        EB_MALLOC_DEC(RestorationLineBuffers**, lr_ctxt->rlbs[i],
            num_planes * sizeof(RestorationLineBuffers**), EB_N_PTR);
        p_rlbs = lr_ctxt->rlbs[i];
        for (int32_t pli = 0; pli < num_planes; pli++) {
            EB_MALLOC_DEC(RestorationLineBuffers *, p_rlbs[pli],
                sizeof(RestorationLineBuffers), EB_N_PTR)
        }
    }
    for (uint32_t i = 0; i < num_instances; i++) {
        EB_MALLOC_DEC(int32_t *, lr_ctxt->rst_tmpbuf[i],
            RESTORATION_TMPBUF_SIZE, EB_N_PTR)
    }

    int frame_width = dec_handle_ptr->seq_header.max_frame_width;
    int frame_height = dec_handle_ptr->seq_header.max_frame_height;
    int sub_x = dec_handle_ptr->seq_header.color_config.subsampling_x;

    // Allocate memory for Deblocked line buffer around stripe(64) boundary for a frame
    const int ext_h = RESTORATION_UNIT_OFFSET + frame_height;
    const int num_stripes = (ext_h + 63) / 64;
    int use_highbd = (dec_handle_ptr->seq_header.color_config.bit_depth > EB_8BIT ||
        dec_handle_ptr->is_16bit_pipeline);

    for (int plane = 0; plane < num_planes; plane++)
    {
        const int is_uv = plane > 0;
        const int ss_x = is_uv && sub_x;
        const int plane_w = ((frame_width + ss_x) >> ss_x) + 2 * RESTORATION_EXTRA_HORZ;
        const int stride = ALIGN_POWER_OF_TWO(plane_w, 5);
        const int buf_size = num_stripes * stride * RESTORATION_CTX_VERT << use_highbd;
        RestorationStripeBoundaries *boundaries = &lr_ctxt->boundaries[plane];

        EB_MALLOC_DEC(uint8_t *, boundaries->stripe_boundary_above,
                      buf_size * sizeof(uint8_t), EB_N_PTR);
        EB_MALLOC_DEC(uint8_t *, boundaries->stripe_boundary_below,
                      buf_size * sizeof(uint8_t), EB_N_PTR);
        boundaries->stripe_boundary_size = buf_size;
        boundaries->stripe_boundary_stride = stride;
    }
    EB_MALLOC_DEC(uint8_t *, lr_ctxt->dst, (MAX_SB_SIZE + 8) *
        RESTORATION_PROC_UNIT_SIZE * sizeof(uint8_t) << use_highbd, EB_N_PTR);
    return return_error;
}

EbErrorType dec_mem_init(EbDecHandle  *dec_handle_ptr) {
    EbErrorType return_error = EB_ErrorNone;

    if (0 == dec_handle_ptr->seq_header_done)
        return EB_ErrorNone;

    /* init module ctxts */
    return_error |= dec_pic_mgr_init(dec_handle_ptr);

    return_error |= init_parse_context(dec_handle_ptr);

    return_error |= init_dec_mod_ctxt(dec_handle_ptr,
                    &dec_handle_ptr->pv_dec_mod_ctxt);

    return_error |= init_lf_ctxt(dec_handle_ptr);

    return_error |= init_lr_ctxt(dec_handle_ptr);

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
// clang-format on
