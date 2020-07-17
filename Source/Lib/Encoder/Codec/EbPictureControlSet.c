/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

/*
* Copyright (c) 2016, Alliance for Open Media. All rights reserved
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at www.aomedia.org/license/software. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at www.aomedia.org/license/patent.
*/

#include <stdlib.h>

#include "EbDefinitions.h"
#include "EbPictureControlSet.h"
#include "EbSequenceControlSet.h"
#include "EbPictureBufferDesc.h"
#include "EbUtility.h"

void set_tile_info(PictureParentControlSet * pcs_ptr);

void *eb_aom_memalign(size_t align, size_t size);
void  eb_aom_free(void *memblk);
void *eb_aom_malloc(size_t size);

EbErrorType eb_av1_alloc_restoration_buffers(Av1Common *cm);

EbErrorType av1_hash_table_create(HashTable *p_hash_table);

static void set_restoration_unit_size(int32_t width, int32_t height, int32_t sx, int32_t sy,
                                      RestorationInfo *rst) {
    (void)width;
    (void)height;
    (void)sx;
    (void)sy;

    int32_t s = 0;

    if (width * height > 352 * 288)
        rst[0].restoration_unit_size = RESTORATION_UNITSIZE_MAX;
    else
        rst[0].restoration_unit_size = (RESTORATION_UNITSIZE_MAX >> 1);

    rst[1].restoration_unit_size = rst[0].restoration_unit_size >> s;
    rst[2].restoration_unit_size = rst[1].restoration_unit_size;
}

static void segmentation_map_dctor(EbPtr p) {
    SegmentationNeighborMap *obj = (SegmentationNeighborMap *)p;
    EB_FREE_ARRAY(obj->data);
}

static void eb_pcs_sb_structs_dctor(EbPtr p) {
    PictureParentControlSet *obj = (PictureParentControlSet *)p;
    EB_FREE_ARRAY(obj->sb_params_array);
    EB_FREE_ARRAY(obj->sb_geom);
}

EbErrorType segmentation_map_ctor(SegmentationNeighborMap *seg_neighbor_map, uint16_t pic_width,
                                  uint16_t pic_height) {
    uint32_t num_elements = (pic_width >> MI_SIZE_LOG2) * (pic_height >> MI_SIZE_LOG2);

    seg_neighbor_map->dctor = segmentation_map_dctor;

    seg_neighbor_map->map_size = num_elements;
    EB_CALLOC_ARRAY(seg_neighbor_map->data, num_elements);
    return EB_ErrorNone;
}

static void me_sb_results_dctor(EbPtr p) {
    MeSbResults *obj = (MeSbResults *)p;
#if  ME_MEM_OPT
    EB_FREE_ARRAY(obj->me_candidate_array);
    EB_FREE_ARRAY(obj->me_mv_array);
#else
    EB_FREE_ARRAY(obj->me_candidate);
    if (obj->me_mv_array) { EB_FREE_ARRAY(obj->me_mv_array[0]); }
    EB_FREE_ARRAY(obj->me_mv_array);
    EB_FREE_ARRAY(obj->me_candidate_array);
#endif
    EB_FREE_ARRAY(obj->total_me_candidate_index);
}
#if NSQ_REMOVAL_CODE_CLEAN_UP
#if REMOVE_MRP_MODE
EbErrorType me_sb_results_ctor(MeSbResults *obj_ptr) {
#else
EbErrorType me_sb_results_ctor(MeSbResults *obj_ptr, uint8_t mrp_mode,
                               uint32_t maxNumberOfMeCandidatesPerPU) {
#endif
#else
EbErrorType me_sb_results_ctor(MeSbResults *obj_ptr, uint32_t max_number_of_blks_per_sb,
                               uint8_t mrp_mode, uint32_t maxNumberOfMeCandidatesPerPU) {
#endif
    uint32_t pu_index;
#if  !REMOVE_MRP_MODE
    size_t count                      = ((mrp_mode == 0) ? ME_MV_MRP_MODE_0 : ME_MV_MRP_MODE_1);
#endif
    obj_ptr->dctor                    = me_sb_results_dctor;
#if !NSQ_REMOVAL_CODE_CLEAN_UP
    obj_ptr->max_number_of_pus_per_sb = max_number_of_blks_per_sb;
#endif
#if ME_MEM_OPT
#if NSQ_REMOVAL_CODE_CLEAN_UP
#if REMOVE_MRP_MODE
    EB_MALLOC_ARRAY(obj_ptr->me_mv_array, SQUARE_PU_COUNT * MAX_PA_ME_MV);
    EB_MALLOC_ARRAY(obj_ptr->me_candidate_array, SQUARE_PU_COUNT * MAX_PA_ME_CAND);
#else
    EB_MALLOC_ARRAY(obj_ptr->me_mv_array, SQUARE_PU_COUNT * count);
    EB_MALLOC_ARRAY(obj_ptr->me_candidate_array, SQUARE_PU_COUNT * maxNumberOfMeCandidatesPerPU);
#endif
#else
    EB_MALLOC_ARRAY(obj_ptr->me_mv_array, max_number_of_blks_per_sb * count);
    EB_MALLOC_ARRAY(obj_ptr->me_candidate_array,
        max_number_of_blks_per_sb * maxNumberOfMeCandidatesPerPU);
#endif
#else
    EB_MALLOC_ARRAY(obj_ptr->me_candidate, max_number_of_blks_per_sb);
    EB_MALLOC_ARRAY(obj_ptr->me_mv_array, max_number_of_blks_per_sb);
    EB_MALLOC_ARRAY(obj_ptr->me_candidate_array,
                    max_number_of_blks_per_sb * maxNumberOfMeCandidatesPerPU);
    EB_MALLOC_ARRAY(obj_ptr->me_mv_array[0], max_number_of_blks_per_sb * count);
#endif
#if NSQ_REMOVAL_CODE_CLEAN_UP
    for (pu_index = 0; pu_index < SQUARE_PU_COUNT; ++pu_index) {
#else
    for (pu_index = 0; pu_index < max_number_of_blks_per_sb; ++pu_index) {
#endif
#if  ME_MEM_OPT
#if REMOVE_MRP_MODE
        obj_ptr->me_candidate_array[pu_index*MAX_PA_ME_CAND + 0].ref_idx_l0 = 0;
        obj_ptr->me_candidate_array[pu_index*MAX_PA_ME_CAND + 0].ref_idx_l1 = 0;
        obj_ptr->me_candidate_array[pu_index*MAX_PA_ME_CAND + 1].ref_idx_l0 = 0;
        obj_ptr->me_candidate_array[pu_index*MAX_PA_ME_CAND + 1].ref_idx_l1 = 0;
        obj_ptr->me_candidate_array[pu_index*MAX_PA_ME_CAND + 2].ref_idx_l0 = 0;
        obj_ptr->me_candidate_array[pu_index*MAX_PA_ME_CAND + 2].ref_idx_l1 = 0;
        obj_ptr->me_candidate_array[pu_index*MAX_PA_ME_CAND + 0].direction = 0;
        obj_ptr->me_candidate_array[pu_index*MAX_PA_ME_CAND + 1].direction = 1;
        obj_ptr->me_candidate_array[pu_index*MAX_PA_ME_CAND + 2].direction = 2;
#else
        obj_ptr->me_candidate_array[pu_index*maxNumberOfMeCandidatesPerPU + 0].ref_idx_l0 = 0;
        obj_ptr->me_candidate_array[pu_index*maxNumberOfMeCandidatesPerPU + 0].ref_idx_l1 = 0;
        obj_ptr->me_candidate_array[pu_index*maxNumberOfMeCandidatesPerPU + 1].ref_idx_l0 = 0;
        obj_ptr->me_candidate_array[pu_index*maxNumberOfMeCandidatesPerPU + 1].ref_idx_l1 = 0;
        obj_ptr->me_candidate_array[pu_index*maxNumberOfMeCandidatesPerPU + 2].ref_idx_l0 = 0;
        obj_ptr->me_candidate_array[pu_index*maxNumberOfMeCandidatesPerPU + 2].ref_idx_l1 = 0;
        obj_ptr->me_candidate_array[pu_index*maxNumberOfMeCandidatesPerPU + 0].direction = 0;
        obj_ptr->me_candidate_array[pu_index*maxNumberOfMeCandidatesPerPU + 1].direction = 1;
        obj_ptr->me_candidate_array[pu_index*maxNumberOfMeCandidatesPerPU + 2].direction = 2;
#endif
#else
        obj_ptr->me_candidate[pu_index] =
            &obj_ptr->me_candidate_array[pu_index * maxNumberOfMeCandidatesPerPU];

        obj_ptr->me_candidate[pu_index][0].ref_idx_l0 = 0;
        obj_ptr->me_candidate[pu_index][0].ref_idx_l1 = 0;
        obj_ptr->me_candidate[pu_index][1].ref_idx_l0 = 0;
        obj_ptr->me_candidate[pu_index][1].ref_idx_l1 = 0;
        obj_ptr->me_candidate[pu_index][2].ref_idx_l0 = 0;
        obj_ptr->me_candidate[pu_index][2].ref_idx_l1 = 0;

        obj_ptr->me_candidate[pu_index][0].direction = 0;
        obj_ptr->me_candidate[pu_index][1].direction = 1;
        obj_ptr->me_candidate[pu_index][2].direction = 2;
        obj_ptr->me_mv_array[pu_index]               = obj_ptr->me_mv_array[0] + pu_index * count;
#endif
    }
#if NSQ_REMOVAL_CODE_CLEAN_UP
    EB_MALLOC_ARRAY(obj_ptr->total_me_candidate_index, SQUARE_PU_COUNT);
#else
    EB_MALLOC_ARRAY(obj_ptr->total_me_candidate_index, max_number_of_blks_per_sb);
#endif
    return EB_ErrorNone;
}

void picture_control_set_dctor(EbPtr p) {
    PictureControlSet *obj = (PictureControlSet *)p;
    uint16_t tile_cnt = obj->tile_row_count * obj->tile_column_count;
    uint8_t            depth;
    av1_hash_table_destroy(&obj->hash_table);
    EB_FREE_ALIGNED_ARRAY(obj->tpl_mvs);
    EB_FREE_ALIGNED(obj->rst_tmpbuf);
    EB_DELETE_PTR_ARRAY(obj->enc_dec_segment_ctrl, tile_cnt);
    EB_DELETE_PTR_ARRAY(obj->ep_intra_luma_mode_neighbor_array, tile_cnt);
    EB_DELETE_PTR_ARRAY(obj->ep_intra_chroma_mode_neighbor_array, tile_cnt);
    EB_DELETE_PTR_ARRAY(obj->ep_mv_neighbor_array, tile_cnt);
    EB_DELETE_PTR_ARRAY(obj->ep_skip_flag_neighbor_array, tile_cnt);
    EB_DELETE_PTR_ARRAY(obj->ep_mode_type_neighbor_array, tile_cnt);
    EB_DELETE_PTR_ARRAY(obj->ep_leaf_depth_neighbor_array, tile_cnt);
    EB_DELETE_PTR_ARRAY(obj->ep_luma_recon_neighbor_array, tile_cnt);
    EB_DELETE_PTR_ARRAY(obj->ep_cb_recon_neighbor_array, tile_cnt);
    EB_DELETE_PTR_ARRAY(obj->ep_cr_recon_neighbor_array, tile_cnt);
    EB_DELETE_PTR_ARRAY(obj->ep_luma_dc_sign_level_coeff_neighbor_array, tile_cnt);
    EB_DELETE_PTR_ARRAY(obj->ep_cb_dc_sign_level_coeff_neighbor_array, tile_cnt);
    EB_DELETE_PTR_ARRAY(obj->ep_cr_dc_sign_level_coeff_neighbor_array, tile_cnt);
    EB_DELETE_PTR_ARRAY(obj->mode_type_neighbor_array, tile_cnt);
    EB_DELETE_PTR_ARRAY(obj->partition_context_neighbor_array, tile_cnt);
    EB_DELETE_PTR_ARRAY(obj->skip_flag_neighbor_array, tile_cnt);
    EB_DELETE_PTR_ARRAY(obj->skip_coeff_neighbor_array, tile_cnt);
    EB_DELETE_PTR_ARRAY(obj->luma_dc_sign_level_coeff_neighbor_array, tile_cnt);
    EB_DELETE_PTR_ARRAY(obj->cr_dc_sign_level_coeff_neighbor_array, tile_cnt);
    EB_DELETE_PTR_ARRAY(obj->cb_dc_sign_level_coeff_neighbor_array, tile_cnt);
    EB_DELETE_PTR_ARRAY(obj->inter_pred_dir_neighbor_array, tile_cnt);
    EB_DELETE_PTR_ARRAY(obj->ref_frame_type_neighbor_array, tile_cnt);
    EB_DELETE_PTR_ARRAY(obj->intra_luma_mode_neighbor_array, tile_cnt);
    EB_DELETE_PTR_ARRAY(obj->txfm_context_array, tile_cnt);
    EB_DELETE_PTR_ARRAY(obj->segmentation_id_pred_array, tile_cnt);
    EB_DELETE(obj->segmentation_neighbor_map); //Jing, double check here
    EB_DELETE_PTR_ARRAY(obj->ep_luma_recon_neighbor_array16bit, tile_cnt);
    EB_DELETE_PTR_ARRAY(obj->ep_cb_recon_neighbor_array16bit, tile_cnt);
    EB_DELETE_PTR_ARRAY(obj->ep_cr_recon_neighbor_array16bit, tile_cnt);
    EB_DELETE_PTR_ARRAY(obj->interpolation_type_neighbor_array, tile_cnt);

    //EB_DELETE(obj->ep_partition_context_neighbor_array); //Jing: Double check here
    EB_DELETE_PTR_ARRAY(obj->ep_partition_context_neighbor_array, tile_cnt);

    for (depth = 0; depth < NEIGHBOR_ARRAY_TOTAL_COUNT; depth++) {
        EB_DELETE_PTR_ARRAY(obj->md_intra_luma_mode_neighbor_array[depth], tile_cnt);
        EB_DELETE_PTR_ARRAY(obj->md_intra_chroma_mode_neighbor_array[depth], tile_cnt);
        EB_DELETE_PTR_ARRAY(obj->md_mv_neighbor_array[depth], tile_cnt);
        EB_DELETE_PTR_ARRAY(obj->md_skip_flag_neighbor_array[depth], tile_cnt);
        EB_DELETE_PTR_ARRAY(obj->md_mode_type_neighbor_array[depth], tile_cnt);
        EB_DELETE_PTR_ARRAY(obj->md_leaf_depth_neighbor_array[depth], tile_cnt);
        EB_DELETE_PTR_ARRAY(obj->mdleaf_partition_neighbor_array[depth], tile_cnt);

        if (obj->hbd_mode_decision > EB_8_BIT_MD) {
            EB_DELETE_PTR_ARRAY(obj->md_luma_recon_neighbor_array16bit[depth], tile_cnt);
            EB_DELETE_PTR_ARRAY(obj->md_tx_depth_1_luma_recon_neighbor_array16bit[depth], tile_cnt);
            EB_DELETE_PTR_ARRAY(obj->md_tx_depth_2_luma_recon_neighbor_array16bit[depth], tile_cnt);
            EB_DELETE_PTR_ARRAY(obj->md_cb_recon_neighbor_array16bit[depth], tile_cnt);
            EB_DELETE_PTR_ARRAY(obj->md_cr_recon_neighbor_array16bit[depth], tile_cnt);
        }
        if (obj->hbd_mode_decision != EB_10_BIT_MD) {
            EB_DELETE_PTR_ARRAY(obj->md_luma_recon_neighbor_array[depth], tile_cnt);
            EB_DELETE_PTR_ARRAY(obj->md_tx_depth_1_luma_recon_neighbor_array[depth], tile_cnt);
            EB_DELETE_PTR_ARRAY(obj->md_tx_depth_2_luma_recon_neighbor_array[depth], tile_cnt);
            EB_DELETE_PTR_ARRAY(obj->md_cb_recon_neighbor_array[depth], tile_cnt);
            EB_DELETE_PTR_ARRAY(obj->md_cr_recon_neighbor_array[depth], tile_cnt);
        }

        EB_DELETE_PTR_ARRAY(obj->md_skip_coeff_neighbor_array[depth], tile_cnt);
        EB_DELETE_PTR_ARRAY(obj->md_luma_dc_sign_level_coeff_neighbor_array[depth], tile_cnt);
        EB_DELETE_PTR_ARRAY(obj->md_tx_depth_1_luma_dc_sign_level_coeff_neighbor_array[depth], tile_cnt);
        EB_DELETE_PTR_ARRAY(obj->md_cr_dc_sign_level_coeff_neighbor_array[depth], tile_cnt);
        EB_DELETE_PTR_ARRAY(obj->md_cb_dc_sign_level_coeff_neighbor_array[depth], tile_cnt);
        EB_DELETE_PTR_ARRAY(obj->md_txfm_context_array[depth], tile_cnt);
        EB_DELETE_PTR_ARRAY(obj->md_inter_pred_dir_neighbor_array[depth], tile_cnt);
        EB_DELETE_PTR_ARRAY(obj->md_ref_frame_type_neighbor_array[depth], tile_cnt);
        EB_DELETE_PTR_ARRAY(obj->md_interpolation_type_neighbor_array[depth], tile_cnt);
    }
    EB_DELETE_PTR_ARRAY(obj->sb_ptr_array, obj->sb_total_count_unscaled);
#if !MD_FRAME_CONTEXT_MEM_OPT
    EB_DELETE(obj->coeff_est_entropy_coder_ptr);
#endif
    EB_DELETE(obj->bitstream_ptr);
    EB_DELETE_PTR_ARRAY(obj->entropy_coding_info, tile_cnt);
#if !PCS_MEM_OPT
    EB_DELETE(obj->recon_picture32bit_ptr);
#endif
    EB_DELETE(obj->recon_picture16bit_ptr);
    EB_DELETE(obj->recon_picture_ptr);
    EB_DELETE(obj->film_grain_picture16bit_ptr);
    EB_DELETE(obj->film_grain_picture_ptr);
    EB_DELETE(obj->input_frame16bit);

    EB_FREE_ARRAY(obj->mse_seg[0]);
    EB_FREE_ARRAY(obj->mse_seg[1]);

    EB_FREE_ARRAY(obj->mi_grid_base);
    EB_FREE_ARRAY(obj->mip);
    EB_FREE_ARRAY(obj->md_rate_estimation_array);
    EB_FREE_ARRAY(obj->ec_ctx_array);
#if !REU_MEM_OPT
    EB_FREE_ARRAY(obj->rate_est_array);
#endif
#if !PAL_MEM_OPT
    if (obj->tile_tok[0][0]) EB_FREE_ARRAY(obj->tile_tok[0][0]);
#endif
#if !DEPTH_PART_CLEAN_UP
    EB_FREE_ARRAY(obj->mdc_sb_array);
#endif
    EB_DESTROY_MUTEX(obj->entropy_coding_pic_mutex);
    EB_DESTROY_MUTEX(obj->intra_mutex);
    EB_DESTROY_MUTEX(obj->cdef_search_mutex);
    EB_DESTROY_MUTEX(obj->rest_search_mutex);
}
// Token buffer is only used for palette tokens.
static INLINE unsigned int get_token_alloc(int mb_rows, int mb_cols, int sb_size_log2,
                                           const int num_planes) {
    // Calculate the maximum number of max superblocks in the image.
    const int shift          = sb_size_log2 - 4;
    const int sb_size        = 1 << sb_size_log2;
    const int sb_size_square = sb_size * sb_size;
    const int sb_rows        = ALIGN_POWER_OF_TWO(mb_rows, shift) >> shift;
    const int sb_cols        = ALIGN_POWER_OF_TWO(mb_cols, shift) >> shift;

    // One palette token for each pixel. There can be palettes on two planes.
    const int sb_palette_toks = AOMMIN(2, num_planes) * sb_size_square;

    return sb_rows * sb_cols * sb_palette_toks;
}

typedef struct InitData {
    NeighborArrayUnit **na_unit_dbl_ptr;
    uint32_t            max_picture_width;
    uint32_t            max_picture_height;
    uint32_t            unit_size;
    uint32_t            granularity_normal;
    uint32_t            granularity_top_left;
    uint32_t            type_mask;
} InitData;

#define DIM(array) (sizeof(array) / sizeof(array[0]))
EbErrorType create_neighbor_array_units(InitData *data, size_t count) {
    for (size_t i = 0; i < count; i++) {
        EB_NEW(*data[i].na_unit_dbl_ptr,
               neighbor_array_unit_ctor,
               data[i].max_picture_width,
               data[i].max_picture_height,
               data[i].unit_size,
               data[i].granularity_normal,
               data[i].granularity_top_left,
               data[i].type_mask);
    }
    return EB_ErrorNone;
}

#if PAL_MEM_OPT
EbErrorType  alloc_palette_tokens(SequenceControlSet * scs_ptr, PictureControlSet * child_pcs_ptr)
{
#if UPDATE_SC_DETECTION
    if (child_pcs_ptr->parent_pcs_ptr->frm_hdr.allow_screen_content_tools){
#else
    if (child_pcs_ptr->parent_pcs_ptr->sc_content_detected) {
#endif
        if (scs_ptr->static_config.screen_content_mode) {
            uint32_t     mi_cols = scs_ptr->max_input_luma_width >> MI_SIZE_LOG2;
            uint32_t     mi_rows = scs_ptr->max_input_luma_height >> MI_SIZE_LOG2;
            uint32_t     mb_cols = (mi_cols + 2) >> 2;
            uint32_t     mb_rows = (mi_rows + 2) >> 2;
            unsigned int tokens = get_token_alloc(mb_rows, mb_cols, MAX_SB_SIZE_LOG2, 2);
            EB_CALLOC_ARRAY(child_pcs_ptr->tile_tok[0][0], tokens);
        }
        else
            child_pcs_ptr->tile_tok[0][0] = NULL;
    }

    return EB_ErrorNone;
}
#endif
EbErrorType picture_control_set_ctor(PictureControlSet *object_ptr, EbPtr object_init_data_ptr) {
    PictureControlSetInitData *init_data_ptr = (PictureControlSetInitData *)object_init_data_ptr;

    EbPictureBufferDescInitData input_pic_buf_desc_init_data;
    EbPictureBufferDescInitData coeff_buffer_desc_init_data;

    // Max/Min CU Sizes
    const uint32_t max_blk_size = init_data_ptr->sb_size_pix;
    // SBs
    const uint16_t picture_sb_width = (uint16_t)(
        (init_data_ptr->picture_width + init_data_ptr->sb_sz - 1) / init_data_ptr->sb_sz);
    const uint16_t picture_sb_height = (uint16_t)(
        (init_data_ptr->picture_height + init_data_ptr->sb_sz - 1) / init_data_ptr->sb_sz);
    uint16_t    sb_index;
    uint16_t    sb_origin_x;
    uint16_t    sb_origin_y;
    EbErrorType return_error = EB_ErrorNone;

    EbBool         is_16bit      = init_data_ptr->bit_depth > 8 ? EB_TRUE : EB_FALSE;
    const uint16_t subsampling_x = (init_data_ptr->color_format == EB_YUV444 ? 1 : 2) - 1;
    const uint16_t subsampling_y = (init_data_ptr->color_format >= EB_YUV422 ? 1 : 2) - 1;

    uint32_t total_tile_cnt = init_data_ptr->tile_row_count * init_data_ptr->tile_column_count;
    uint32_t tile_idx = 0;
#if OUTPUT_MEM_OPT
    uint32_t output_buffer_size =
        (uint32_t)(EB_OUTPUTSTREAMBUFFERSIZE_MACRO(init_data_ptr->picture_width * init_data_ptr->picture_height));
#endif
    object_ptr->tile_row_count = init_data_ptr->tile_row_count;
    object_ptr->tile_column_count = init_data_ptr->tile_column_count;

    object_ptr->dctor = picture_control_set_dctor;

    // Init Picture Init data
    input_pic_buf_desc_init_data.max_width          = init_data_ptr->picture_width;
    input_pic_buf_desc_init_data.max_height         = init_data_ptr->picture_height;
    input_pic_buf_desc_init_data.bit_depth          = init_data_ptr->bit_depth;
    input_pic_buf_desc_init_data.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
    input_pic_buf_desc_init_data.color_format       = init_data_ptr->color_format;

    input_pic_buf_desc_init_data.left_padding  = PAD_VALUE;
    input_pic_buf_desc_init_data.right_padding = PAD_VALUE;
    input_pic_buf_desc_init_data.top_padding   = PAD_VALUE;
    input_pic_buf_desc_init_data.bot_padding   = PAD_VALUE;

    input_pic_buf_desc_init_data.split_mode = EB_FALSE;

    coeff_buffer_desc_init_data.max_width          = init_data_ptr->picture_width;
    coeff_buffer_desc_init_data.max_height         = init_data_ptr->picture_height;
    coeff_buffer_desc_init_data.bit_depth          = EB_16BIT;
    coeff_buffer_desc_init_data.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
    coeff_buffer_desc_init_data.color_format       = init_data_ptr->color_format;

    coeff_buffer_desc_init_data.left_padding  = PAD_VALUE;
    coeff_buffer_desc_init_data.right_padding = PAD_VALUE;
    coeff_buffer_desc_init_data.top_padding   = PAD_VALUE;
    coeff_buffer_desc_init_data.bot_padding   = PAD_VALUE;

    coeff_buffer_desc_init_data.split_mode = EB_FALSE;
    coeff_buffer_desc_init_data.is_16bit_pipeline = init_data_ptr->is_16bit_pipeline;

    object_ptr->scs_wrapper_ptr = (EbObjectWrapper *)NULL;

    object_ptr->recon_picture16bit_ptr = (EbPictureBufferDesc *)NULL;
    object_ptr->recon_picture_ptr      = (EbPictureBufferDesc *)NULL;
    object_ptr->color_format           = init_data_ptr->color_format;
#if !PCS_MEM_OPT
    EbPictureBufferDescInitData coeff_buffer_desc_32bit_init_data;
    coeff_buffer_desc_32bit_init_data.max_width          = init_data_ptr->picture_width;
    coeff_buffer_desc_32bit_init_data.max_height         = init_data_ptr->picture_height;
    coeff_buffer_desc_32bit_init_data.bit_depth          = EB_32BIT;
    coeff_buffer_desc_32bit_init_data.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
    coeff_buffer_desc_32bit_init_data.color_format       = init_data_ptr->color_format;
    coeff_buffer_desc_32bit_init_data.left_padding       = 0;
    coeff_buffer_desc_32bit_init_data.right_padding      = 0;
    coeff_buffer_desc_32bit_init_data.top_padding        = 0;
    coeff_buffer_desc_32bit_init_data.bot_padding        = 0;
    coeff_buffer_desc_32bit_init_data.split_mode         = EB_FALSE;

    object_ptr->recon_picture32bit_ptr = (EbPictureBufferDesc *)NULL;
    EB_NEW(object_ptr->recon_picture32bit_ptr,
           eb_recon_picture_buffer_desc_ctor,
           (EbPtr)&coeff_buffer_desc_32bit_init_data);
#endif
    // Reconstructed Picture Buffer
    if (is_16bit) {
        EB_NEW(object_ptr->recon_picture16bit_ptr,
               eb_recon_picture_buffer_desc_ctor,
               (EbPtr)&coeff_buffer_desc_init_data);
    } else {
        EB_NEW(object_ptr->recon_picture_ptr,
               eb_recon_picture_buffer_desc_ctor,
               (EbPtr)&input_pic_buf_desc_init_data);
    }
    if (init_data_ptr->is_16bit_pipeline && !is_16bit) {
        EB_NEW(object_ptr->recon_picture16bit_ptr,
            eb_recon_picture_buffer_desc_ctor,
            (EbPtr)&coeff_buffer_desc_init_data);
    }
    // Film Grain Picture Buffer
    if (init_data_ptr->film_grain_noise_level) {
        if (is_16bit) {
            EB_NEW(object_ptr->film_grain_picture16bit_ptr,
                   eb_recon_picture_buffer_desc_ctor,
                   (EbPtr)&coeff_buffer_desc_init_data);
        } else {
            EB_NEW(object_ptr->film_grain_picture_ptr,
                   eb_recon_picture_buffer_desc_ctor,
                   (EbPtr)&input_pic_buf_desc_init_data);
        }
    }

    if ((is_16bit) || (init_data_ptr->is_16bit_pipeline)) {
        EB_NEW(object_ptr->input_frame16bit,
               eb_picture_buffer_desc_ctor,
               (EbPtr)&coeff_buffer_desc_init_data);
    }
    // Entropy Coder
    EB_ALLOC_PTR_ARRAY(object_ptr->entropy_coding_info, total_tile_cnt);
    for (tile_idx = 0; tile_idx < total_tile_cnt; tile_idx++) {
#if OUTPUT_MEM_OPT
        EB_NEW(object_ptr->entropy_coding_info[tile_idx],
               entropy_tile_info_ctor,
               output_buffer_size / total_tile_cnt);
#else
        EB_NEW(object_ptr->entropy_coding_info[tile_idx],
               entropy_tile_info_ctor,
               SEGMENT_ENTROPY_BUFFER_SIZE / total_tile_cnt);
#endif
    }

    // Packetization process Bitstream
#if OUTPUT_MEM_OPT
    EB_NEW(object_ptr->bitstream_ptr, bitstream_ctor, output_buffer_size);
#else
    EB_NEW(object_ptr->bitstream_ptr, bitstream_ctor, PACKETIZATION_PROCESS_BUFFER_SIZE);
#endif

#if !MD_FRAME_CONTEXT_MEM_OPT
    // Rate estimation entropy coder
    EB_NEW(
        object_ptr->coeff_est_entropy_coder_ptr, entropy_coder_ctor, SEGMENT_ENTROPY_BUFFER_SIZE);
#endif
    // GOP
    object_ptr->picture_number       = 0;
    object_ptr->temporal_layer_index = 0;

    // SB Array
    object_ptr->sb_total_count = picture_sb_width * picture_sb_height;
    object_ptr->sb_total_count_unscaled = object_ptr->sb_total_count;
    EB_ALLOC_PTR_ARRAY(object_ptr->sb_ptr_array, object_ptr->sb_total_count);

    sb_origin_x = 0;
    sb_origin_y = 0;

    const uint16_t picture_sb_w =
        (uint16_t)((init_data_ptr->picture_width + init_data_ptr->sb_size_pix - 1) /
                   init_data_ptr->sb_size_pix);
    const uint16_t picture_sb_h =
        (uint16_t)((init_data_ptr->picture_height + init_data_ptr->sb_size_pix - 1) /
                   init_data_ptr->sb_size_pix);
    const uint16_t all_sb = picture_sb_w * picture_sb_h;

    object_ptr->sb_total_count_pix = all_sb;

    for (sb_index = 0; sb_index < all_sb; ++sb_index) {
        EB_NEW(object_ptr->sb_ptr_array[sb_index],
               largest_coding_unit_ctor,
               (uint8_t)init_data_ptr->sb_size_pix,
               (uint16_t)(sb_origin_x * max_blk_size),
               (uint16_t)(sb_origin_y * max_blk_size),
               (uint16_t)sb_index,
               object_ptr);
        // Increment the Order in coding order (Raster Scan Order)
        sb_origin_y = (sb_origin_x == picture_sb_w - 1) ? sb_origin_y + 1 : sb_origin_y;
        sb_origin_x = (sb_origin_x == picture_sb_w - 1) ? 0 : sb_origin_x + 1;
    }
    // MD Rate Estimation Array
    EB_MALLOC_ARRAY(object_ptr->md_rate_estimation_array, 1);
    memset(object_ptr->md_rate_estimation_array, 0, sizeof(MdRateEstimationContext));
    EB_MALLOC_ARRAY(object_ptr->ec_ctx_array, all_sb);
#if !REU_MEM_OPT
#if RATE_MEM_OPT
    if (init_data_ptr->serial_rate_est)
        EB_MALLOC_ARRAY(object_ptr->rate_est_array, 1);
    else
#endif
        EB_MALLOC_ARRAY(object_ptr->rate_est_array, all_sb);
#endif
#if ! PAL_MEM_OPT
    if (init_data_ptr->cfg_palette) {
        uint32_t     mi_cols = init_data_ptr->picture_width >> MI_SIZE_LOG2;
        uint32_t     mi_rows = init_data_ptr->picture_height >> MI_SIZE_LOG2;
        uint32_t     mb_cols = (mi_cols + 2) >> 2;
        uint32_t     mb_rows = (mi_rows + 2) >> 2;
        unsigned int tokens  = get_token_alloc(mb_rows, mb_cols, MAX_SB_SIZE_LOG2, 2);
        EB_CALLOC_ARRAY(object_ptr->tile_tok[0][0], tokens);
    } else
        object_ptr->tile_tok[0][0] = NULL;
#endif
#if !DEPTH_PART_CLEAN_UP
    // Mode Decision Control config
    EB_MALLOC_ARRAY(object_ptr->mdc_sb_array, object_ptr->sb_total_count);
#endif
#if CHANGE_HBD_MODE
    if (init_data_ptr->hbd_mode_decision == DEFAULT)
        object_ptr->hbd_mode_decision = init_data_ptr->hbd_mode_decision = 2;
    else
#endif
    object_ptr->hbd_mode_decision = init_data_ptr->hbd_mode_decision;
    // Mode Decision Neighbor Arrays
    uint8_t depth;
    for (depth = 0; depth < NEIGHBOR_ARRAY_TOTAL_COUNT; depth++) {
        EB_ALLOC_PTR_ARRAY(object_ptr->md_intra_luma_mode_neighbor_array[depth], total_tile_cnt);
        EB_ALLOC_PTR_ARRAY(object_ptr->md_intra_chroma_mode_neighbor_array[depth], total_tile_cnt);
        EB_ALLOC_PTR_ARRAY(object_ptr->md_mv_neighbor_array[depth], total_tile_cnt);
        EB_ALLOC_PTR_ARRAY(object_ptr->md_skip_flag_neighbor_array[depth], total_tile_cnt);
        EB_ALLOC_PTR_ARRAY(object_ptr->md_mode_type_neighbor_array[depth], total_tile_cnt);
        EB_ALLOC_PTR_ARRAY(object_ptr->md_leaf_depth_neighbor_array[depth], total_tile_cnt);
        EB_ALLOC_PTR_ARRAY(object_ptr->mdleaf_partition_neighbor_array[depth], total_tile_cnt);
        EB_ALLOC_PTR_ARRAY(object_ptr->md_skip_coeff_neighbor_array[depth], total_tile_cnt);
        EB_ALLOC_PTR_ARRAY(object_ptr->md_luma_dc_sign_level_coeff_neighbor_array[depth], total_tile_cnt);
        EB_ALLOC_PTR_ARRAY(object_ptr->md_tx_depth_1_luma_dc_sign_level_coeff_neighbor_array[depth], total_tile_cnt);
        EB_ALLOC_PTR_ARRAY(object_ptr->md_cr_dc_sign_level_coeff_neighbor_array[depth], total_tile_cnt);
        EB_ALLOC_PTR_ARRAY(object_ptr->md_cb_dc_sign_level_coeff_neighbor_array[depth], total_tile_cnt);
        EB_ALLOC_PTR_ARRAY(object_ptr->md_txfm_context_array[depth], total_tile_cnt);
        EB_ALLOC_PTR_ARRAY(object_ptr->md_inter_pred_dir_neighbor_array[depth], total_tile_cnt);
        EB_ALLOC_PTR_ARRAY(object_ptr->md_ref_frame_type_neighbor_array[depth], total_tile_cnt);
        if (init_data_ptr->hbd_mode_decision != EB_10_BIT_MD) {
            EB_ALLOC_PTR_ARRAY(object_ptr->md_luma_recon_neighbor_array[depth], total_tile_cnt);
            EB_ALLOC_PTR_ARRAY(object_ptr->md_tx_depth_1_luma_recon_neighbor_array[depth], total_tile_cnt);
            EB_ALLOC_PTR_ARRAY(object_ptr->md_tx_depth_2_luma_recon_neighbor_array[depth],
                               total_tile_cnt);
            EB_ALLOC_PTR_ARRAY(object_ptr->md_cb_recon_neighbor_array[depth], total_tile_cnt);
            EB_ALLOC_PTR_ARRAY(object_ptr->md_cr_recon_neighbor_array[depth], total_tile_cnt);
        }
        if (init_data_ptr->hbd_mode_decision > EB_8_BIT_MD) {
            EB_ALLOC_PTR_ARRAY(object_ptr->md_luma_recon_neighbor_array16bit[depth], total_tile_cnt);
            EB_ALLOC_PTR_ARRAY(object_ptr->md_tx_depth_1_luma_recon_neighbor_array16bit[depth], total_tile_cnt);
            EB_ALLOC_PTR_ARRAY(object_ptr->md_tx_depth_2_luma_recon_neighbor_array16bit[depth],
                               total_tile_cnt);
            EB_ALLOC_PTR_ARRAY(object_ptr->md_cb_recon_neighbor_array16bit[depth], total_tile_cnt);
            EB_ALLOC_PTR_ARRAY(object_ptr->md_cr_recon_neighbor_array16bit[depth], total_tile_cnt);
        }
        EB_ALLOC_PTR_ARRAY(object_ptr->md_interpolation_type_neighbor_array[depth], total_tile_cnt);
    }

    for (tile_idx = 0; tile_idx < total_tile_cnt; tile_idx++) {
        for (depth = 0; depth < NEIGHBOR_ARRAY_TOTAL_COUNT; depth++) {
            InitData data[] = {
                {&object_ptr->md_intra_luma_mode_neighbor_array[depth][tile_idx],
                    MAX_PICTURE_WIDTH_SIZE,
                    MAX_PICTURE_HEIGHT_SIZE,
                    sizeof(uint8_t),
                    PU_NEIGHBOR_ARRAY_GRANULARITY,
                    PU_NEIGHBOR_ARRAY_GRANULARITY,
                    NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK},
                {
                    &object_ptr->md_intra_chroma_mode_neighbor_array[depth][tile_idx],
                    MAX_PICTURE_WIDTH_SIZE >> subsampling_x,
                    MAX_PICTURE_HEIGHT_SIZE >> subsampling_y,
                    sizeof(uint8_t),
                    PU_NEIGHBOR_ARRAY_GRANULARITY,
                    PU_NEIGHBOR_ARRAY_GRANULARITY,
                    NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
                },
                {
                    &object_ptr->md_mv_neighbor_array[depth][tile_idx],
                    MAX_PICTURE_WIDTH_SIZE,
                    MAX_PICTURE_HEIGHT_SIZE,
                    sizeof(MvUnit),
                    PU_NEIGHBOR_ARRAY_GRANULARITY,
                    PU_NEIGHBOR_ARRAY_GRANULARITY,
                    NEIGHBOR_ARRAY_UNIT_FULL_MASK,
                },
                {
                    &object_ptr->md_skip_flag_neighbor_array[depth][tile_idx],
                    MAX_PICTURE_WIDTH_SIZE,
                    MAX_PICTURE_HEIGHT_SIZE,
                    sizeof(uint8_t),
                    PU_NEIGHBOR_ARRAY_GRANULARITY,
                    PU_NEIGHBOR_ARRAY_GRANULARITY,
                    NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
                },
                {
                    &object_ptr->md_mode_type_neighbor_array[depth][tile_idx],
                    MAX_PICTURE_WIDTH_SIZE,
                    MAX_PICTURE_HEIGHT_SIZE,
                    sizeof(uint8_t),
                    PU_NEIGHBOR_ARRAY_GRANULARITY,
                    PU_NEIGHBOR_ARRAY_GRANULARITY,
                    NEIGHBOR_ARRAY_UNIT_FULL_MASK,
                },
                {
                    &object_ptr->md_leaf_depth_neighbor_array[depth][tile_idx],
                    MAX_PICTURE_WIDTH_SIZE,
                    MAX_PICTURE_HEIGHT_SIZE,
                    sizeof(uint8_t),
                    PU_NEIGHBOR_ARRAY_GRANULARITY,
                    PU_NEIGHBOR_ARRAY_GRANULARITY,
                    NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
                },
                {
                    &object_ptr->mdleaf_partition_neighbor_array[depth][tile_idx],
                    MAX_PICTURE_WIDTH_SIZE,
                    MAX_PICTURE_HEIGHT_SIZE,
                    sizeof(struct PartitionContext),
                    PU_NEIGHBOR_ARRAY_GRANULARITY,
                    PU_NEIGHBOR_ARRAY_GRANULARITY,
                    NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
                },
                {
                    &object_ptr->md_skip_coeff_neighbor_array[depth][tile_idx],
                    MAX_PICTURE_WIDTH_SIZE,
                    MAX_PICTURE_HEIGHT_SIZE,
                    sizeof(uint8_t),
                    PU_NEIGHBOR_ARRAY_GRANULARITY,
                    PU_NEIGHBOR_ARRAY_GRANULARITY,
                    NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
                },
                // for each 4x4
                {
                    &object_ptr->md_luma_dc_sign_level_coeff_neighbor_array[depth][tile_idx],
                    MAX_PICTURE_WIDTH_SIZE,
                    MAX_PICTURE_HEIGHT_SIZE,
                    sizeof(uint8_t),
                    PU_NEIGHBOR_ARRAY_GRANULARITY,
                    PU_NEIGHBOR_ARRAY_GRANULARITY,
                    NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
                },
                // for each 4x4
                {
                    &object_ptr->md_tx_depth_1_luma_dc_sign_level_coeff_neighbor_array[depth][tile_idx],
                    MAX_PICTURE_WIDTH_SIZE,
                    MAX_PICTURE_HEIGHT_SIZE,
                    sizeof(uint8_t),
                    PU_NEIGHBOR_ARRAY_GRANULARITY,
                    PU_NEIGHBOR_ARRAY_GRANULARITY,
                    NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
                },
                // for each 4x4
                {
                    &object_ptr->md_cr_dc_sign_level_coeff_neighbor_array[depth][tile_idx],
                    MAX_PICTURE_WIDTH_SIZE,
                    MAX_PICTURE_HEIGHT_SIZE,
                    sizeof(uint8_t),
                    PU_NEIGHBOR_ARRAY_GRANULARITY,
                    PU_NEIGHBOR_ARRAY_GRANULARITY,
                    NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
                },
                // for each 4x4
                {
                    &object_ptr->md_cb_dc_sign_level_coeff_neighbor_array[depth][tile_idx],
                    MAX_PICTURE_WIDTH_SIZE,
                    MAX_PICTURE_HEIGHT_SIZE,
                    sizeof(uint8_t),
                    PU_NEIGHBOR_ARRAY_GRANULARITY,
                    PU_NEIGHBOR_ARRAY_GRANULARITY,
                    NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
                },
                {
                    &object_ptr->md_txfm_context_array[depth][tile_idx],
                    MAX_PICTURE_WIDTH_SIZE,
                    MAX_PICTURE_HEIGHT_SIZE,
                    sizeof(TXFM_CONTEXT),
                    PU_NEIGHBOR_ARRAY_GRANULARITY,
                    PU_NEIGHBOR_ARRAY_GRANULARITY,
                    NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
                },
                {
                    &object_ptr->md_inter_pred_dir_neighbor_array[depth][tile_idx],
                    MAX_PICTURE_WIDTH_SIZE,
                    MAX_PICTURE_HEIGHT_SIZE,
                    sizeof(uint8_t),
                    PU_NEIGHBOR_ARRAY_GRANULARITY,
                    PU_NEIGHBOR_ARRAY_GRANULARITY,
                    NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
                },
                {
                    &object_ptr->md_ref_frame_type_neighbor_array[depth][tile_idx],
                    MAX_PICTURE_WIDTH_SIZE,
                    MAX_PICTURE_HEIGHT_SIZE,
                    sizeof(uint8_t),
                    PU_NEIGHBOR_ARRAY_GRANULARITY,
                    PU_NEIGHBOR_ARRAY_GRANULARITY,
                    NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
                }};
            return_error = create_neighbor_array_units(data, DIM(data));
            if (return_error == EB_ErrorInsufficientResources) return EB_ErrorInsufficientResources;
            if (init_data_ptr->hbd_mode_decision != EB_10_BIT_MD) {
                InitData data[] = {

                    {
                        &object_ptr->md_luma_recon_neighbor_array[depth][tile_idx],
                        MAX_PICTURE_WIDTH_SIZE,
                        MAX_PICTURE_HEIGHT_SIZE,
                        sizeof(uint8_t),
                        SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                        SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                        NEIGHBOR_ARRAY_UNIT_FULL_MASK,
                    },
                    {
                        &object_ptr->md_tx_depth_1_luma_recon_neighbor_array[depth][tile_idx],
                        MAX_PICTURE_WIDTH_SIZE,
                        MAX_PICTURE_HEIGHT_SIZE,
                        sizeof(uint8_t),
                        SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                        SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                        NEIGHBOR_ARRAY_UNIT_FULL_MASK,
                    },
                    {
                        &object_ptr->md_tx_depth_2_luma_recon_neighbor_array[depth][tile_idx],
                        MAX_PICTURE_WIDTH_SIZE,
                        MAX_PICTURE_HEIGHT_SIZE,
                        sizeof(uint8_t),
                        SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                        SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                        NEIGHBOR_ARRAY_UNIT_FULL_MASK,
                    },
                    {
                        &object_ptr->md_cb_recon_neighbor_array[depth][tile_idx],
                        MAX_PICTURE_WIDTH_SIZE >> subsampling_x,
                        MAX_PICTURE_HEIGHT_SIZE >> subsampling_y,
                        sizeof(uint8_t),
                        SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                        SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                        NEIGHBOR_ARRAY_UNIT_FULL_MASK,
                    },
                    {
                        &object_ptr->md_cr_recon_neighbor_array[depth][tile_idx],
                        MAX_PICTURE_WIDTH_SIZE >> subsampling_x,
                        MAX_PICTURE_HEIGHT_SIZE >> subsampling_y,
                        sizeof(uint8_t),
                        SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                        SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                        NEIGHBOR_ARRAY_UNIT_FULL_MASK,
                    }

                };
                return_error = create_neighbor_array_units(data, DIM(data));
                if (return_error == EB_ErrorInsufficientResources) return EB_ErrorInsufficientResources;
            }
            if (init_data_ptr->hbd_mode_decision > EB_8_BIT_MD) {
                InitData data[] = {
                    {
                        &object_ptr->md_luma_recon_neighbor_array16bit[depth][tile_idx],
                        MAX_PICTURE_WIDTH_SIZE,
                        MAX_PICTURE_HEIGHT_SIZE,
                        sizeof(uint16_t),
                        SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                        SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                        NEIGHBOR_ARRAY_UNIT_FULL_MASK,
                    },
                    {
                        &object_ptr->md_tx_depth_1_luma_recon_neighbor_array16bit[depth][tile_idx],
                        MAX_PICTURE_WIDTH_SIZE,
                        MAX_PICTURE_HEIGHT_SIZE,
                        sizeof(uint16_t),
                        SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                        SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                        NEIGHBOR_ARRAY_UNIT_FULL_MASK,
                    },
                    {
                        &object_ptr->md_tx_depth_2_luma_recon_neighbor_array16bit[depth][tile_idx],
                        MAX_PICTURE_WIDTH_SIZE,
                        MAX_PICTURE_HEIGHT_SIZE,
                        sizeof(uint16_t),
                        SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                        SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                        NEIGHBOR_ARRAY_UNIT_FULL_MASK,
                    },
                    {
                        &object_ptr->md_cb_recon_neighbor_array16bit[depth][tile_idx],
                        MAX_PICTURE_WIDTH_SIZE >> subsampling_x,
                        MAX_PICTURE_HEIGHT_SIZE >> subsampling_y,
                        sizeof(uint16_t),
                        SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                        SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                        NEIGHBOR_ARRAY_UNIT_FULL_MASK,
                    },
                    {
                        &object_ptr->md_cr_recon_neighbor_array16bit[depth][tile_idx],
                        MAX_PICTURE_WIDTH_SIZE >> subsampling_x,
                        MAX_PICTURE_HEIGHT_SIZE >> subsampling_y,
                        sizeof(uint16_t),
                        SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                        SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                        NEIGHBOR_ARRAY_UNIT_FULL_MASK,
                    }
                };
                return_error    = create_neighbor_array_units(data, DIM(data));
                if (return_error == EB_ErrorInsufficientResources) return EB_ErrorInsufficientResources;
            }

            EB_NEW(object_ptr->md_interpolation_type_neighbor_array[depth][tile_idx],
                    neighbor_array_unit_ctor32,
                    MAX_PICTURE_WIDTH_SIZE,
                    MAX_PICTURE_HEIGHT_SIZE,
                    sizeof(uint32_t),
                    PU_NEIGHBOR_ARRAY_GRANULARITY,
                    PU_NEIGHBOR_ARRAY_GRANULARITY,
                    NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
        }
    }

    // EncDec Neighbor
    //EncDec
    EB_ALLOC_PTR_ARRAY(object_ptr->ep_intra_luma_mode_neighbor_array, total_tile_cnt);
    EB_ALLOC_PTR_ARRAY(object_ptr->ep_intra_chroma_mode_neighbor_array, total_tile_cnt);
    EB_ALLOC_PTR_ARRAY(object_ptr->ep_mv_neighbor_array, total_tile_cnt);
    EB_ALLOC_PTR_ARRAY(object_ptr->ep_skip_flag_neighbor_array, total_tile_cnt);
    EB_ALLOC_PTR_ARRAY(object_ptr->ep_mode_type_neighbor_array, total_tile_cnt);
    EB_ALLOC_PTR_ARRAY(object_ptr->ep_leaf_depth_neighbor_array, total_tile_cnt);
    EB_ALLOC_PTR_ARRAY(object_ptr->ep_luma_recon_neighbor_array, total_tile_cnt);
    EB_ALLOC_PTR_ARRAY(object_ptr->ep_cb_recon_neighbor_array, total_tile_cnt);
    EB_ALLOC_PTR_ARRAY(object_ptr->ep_cr_recon_neighbor_array, total_tile_cnt);
    EB_ALLOC_PTR_ARRAY(object_ptr->ep_luma_dc_sign_level_coeff_neighbor_array, total_tile_cnt);
    EB_ALLOC_PTR_ARRAY(object_ptr->ep_cb_dc_sign_level_coeff_neighbor_array, total_tile_cnt);
    EB_ALLOC_PTR_ARRAY(object_ptr->ep_cr_dc_sign_level_coeff_neighbor_array, total_tile_cnt);
    EB_ALLOC_PTR_ARRAY(object_ptr->ep_partition_context_neighbor_array, total_tile_cnt);

    //Entropy
    EB_ALLOC_PTR_ARRAY(object_ptr->mode_type_neighbor_array, total_tile_cnt);
    EB_ALLOC_PTR_ARRAY(object_ptr->partition_context_neighbor_array, total_tile_cnt);
    EB_ALLOC_PTR_ARRAY(object_ptr->skip_flag_neighbor_array, total_tile_cnt);
    EB_ALLOC_PTR_ARRAY(object_ptr->skip_coeff_neighbor_array, total_tile_cnt);
    EB_ALLOC_PTR_ARRAY(object_ptr->luma_dc_sign_level_coeff_neighbor_array, total_tile_cnt);
    EB_ALLOC_PTR_ARRAY(object_ptr->cr_dc_sign_level_coeff_neighbor_array, total_tile_cnt);
    EB_ALLOC_PTR_ARRAY(object_ptr->cb_dc_sign_level_coeff_neighbor_array, total_tile_cnt);
    EB_ALLOC_PTR_ARRAY(object_ptr->inter_pred_dir_neighbor_array, total_tile_cnt);
    EB_ALLOC_PTR_ARRAY(object_ptr->ref_frame_type_neighbor_array, total_tile_cnt);
    EB_ALLOC_PTR_ARRAY(object_ptr->intra_luma_mode_neighbor_array, total_tile_cnt);
    EB_ALLOC_PTR_ARRAY(object_ptr->txfm_context_array, total_tile_cnt);
    EB_ALLOC_PTR_ARRAY(object_ptr->segmentation_id_pred_array, total_tile_cnt);
    if ((is_16bit) || (init_data_ptr->is_16bit_pipeline)) {
        EB_ALLOC_PTR_ARRAY(object_ptr->ep_luma_recon_neighbor_array16bit, total_tile_cnt);
        EB_ALLOC_PTR_ARRAY(object_ptr->ep_cb_recon_neighbor_array16bit, total_tile_cnt);
        EB_ALLOC_PTR_ARRAY(object_ptr->ep_cr_recon_neighbor_array16bit, total_tile_cnt);
    }
    EB_ALLOC_PTR_ARRAY(object_ptr->interpolation_type_neighbor_array, total_tile_cnt);

    for (tile_idx = 0; tile_idx < total_tile_cnt; tile_idx++) {
        InitData data[] = {
            {
                &object_ptr->ep_intra_luma_mode_neighbor_array[tile_idx],
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            // Encode Pass Neighbor Arrays
            {
                &object_ptr->ep_intra_chroma_mode_neighbor_array[tile_idx],
                MAX_PICTURE_WIDTH_SIZE >> subsampling_x,
                MAX_PICTURE_HEIGHT_SIZE >> subsampling_y,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            {
                &object_ptr->ep_mv_neighbor_array[tile_idx],
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(MvUnit),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_FULL_MASK,
            },
            {
                &object_ptr->ep_skip_flag_neighbor_array[tile_idx],
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                CU_NEIGHBOR_ARRAY_GRANULARITY,
                CU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            {
                &object_ptr->ep_mode_type_neighbor_array[tile_idx],
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_FULL_MASK,
            },
            {
                &object_ptr->ep_leaf_depth_neighbor_array[tile_idx],
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            {
                &object_ptr->ep_luma_recon_neighbor_array[tile_idx],
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_FULL_MASK,
            },
            {
                &object_ptr->ep_cb_recon_neighbor_array[tile_idx],
                MAX_PICTURE_WIDTH_SIZE >> subsampling_x,
                MAX_PICTURE_HEIGHT_SIZE >> subsampling_y,
                sizeof(uint8_t),
                SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_FULL_MASK,
            },
            {
                &object_ptr->ep_cr_recon_neighbor_array[tile_idx],
                MAX_PICTURE_WIDTH_SIZE >> subsampling_x,
                MAX_PICTURE_HEIGHT_SIZE >> subsampling_y,
                sizeof(uint8_t),
                SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_FULL_MASK,
            },
            // for each 4x4
            {
                &object_ptr->ep_luma_dc_sign_level_coeff_neighbor_array[tile_idx],
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            // for each 4x4
            {
                &object_ptr->ep_cb_dc_sign_level_coeff_neighbor_array[tile_idx],
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            // for each 4x4
            {
                &object_ptr->ep_cr_dc_sign_level_coeff_neighbor_array[tile_idx],
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            // Encode pass partition neighbor array
            {
                &object_ptr->ep_partition_context_neighbor_array[tile_idx],
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(struct PartitionContext),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },

            // Entropy Coding Neighbor Arrays
            {
                &object_ptr->mode_type_neighbor_array[tile_idx],
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            {
                &object_ptr->partition_context_neighbor_array[tile_idx],
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(struct PartitionContext),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            {
                &object_ptr->skip_flag_neighbor_array[tile_idx],
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            {
                &object_ptr->skip_coeff_neighbor_array[tile_idx],
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            // for each 4x4
            {
                &object_ptr->luma_dc_sign_level_coeff_neighbor_array[tile_idx],
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            // for each 4x4
            {
                &object_ptr->cr_dc_sign_level_coeff_neighbor_array[tile_idx],
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            // for each 4x4
            {
                &object_ptr->cb_dc_sign_level_coeff_neighbor_array[tile_idx],
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            {
                &object_ptr->inter_pred_dir_neighbor_array[tile_idx],
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            {
                &object_ptr->ref_frame_type_neighbor_array[tile_idx],
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            {
                &object_ptr->intra_luma_mode_neighbor_array[tile_idx],
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            {
                &object_ptr->txfm_context_array[tile_idx],
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(TXFM_CONTEXT),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            {
                &object_ptr->segmentation_id_pred_array[tile_idx],
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_FULL_MASK,
            },
        };
        return_error = create_neighbor_array_units(data, DIM(data));
        if (return_error == EB_ErrorInsufficientResources) return EB_ErrorInsufficientResources;

        if ((is_16bit) || (init_data_ptr->is_16bit_pipeline)) {
            InitData data[] = {
                {
                    &object_ptr->ep_luma_recon_neighbor_array16bit[tile_idx],
                    MAX_PICTURE_WIDTH_SIZE,
                    MAX_PICTURE_HEIGHT_SIZE,
                    sizeof(uint16_t),
                    SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                    SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                    NEIGHBOR_ARRAY_UNIT_FULL_MASK,
                },
                {
                    &object_ptr->ep_cb_recon_neighbor_array16bit[tile_idx],
                    MAX_PICTURE_WIDTH_SIZE >> subsampling_x,
                    MAX_PICTURE_HEIGHT_SIZE >> subsampling_y,
                    sizeof(uint16_t),
                    SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                    SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                    NEIGHBOR_ARRAY_UNIT_FULL_MASK,
                },
                {
                    &object_ptr->ep_cr_recon_neighbor_array16bit[tile_idx],
                    MAX_PICTURE_WIDTH_SIZE >> subsampling_x,
                    MAX_PICTURE_HEIGHT_SIZE >> subsampling_y,
                    sizeof(uint16_t),
                    SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                    SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                    NEIGHBOR_ARRAY_UNIT_FULL_MASK,
                },
            };
            return_error = create_neighbor_array_units(data, DIM(data));
            if (return_error == EB_ErrorInsufficientResources) return EB_ErrorInsufficientResources;
        } else {
            object_ptr->ep_luma_recon_neighbor_array16bit = 0;
            object_ptr->ep_cb_recon_neighbor_array16bit   = 0;
            object_ptr->ep_cr_recon_neighbor_array16bit   = 0;
        }
        EB_NEW(object_ptr->interpolation_type_neighbor_array[tile_idx],
                neighbor_array_unit_ctor32,
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint32_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
    }

    //Segmentation neighbor arrays
    EB_NEW(object_ptr->segmentation_neighbor_map,
           segmentation_map_ctor,
           init_data_ptr->picture_width,
           init_data_ptr->picture_height);
    // Segments
    object_ptr->enc_dec_coded_sb_count = 0;

    EB_MALLOC_ARRAY(object_ptr->enc_dec_segment_ctrl, total_tile_cnt);

    for (tile_idx = 0; tile_idx < total_tile_cnt; tile_idx++) {
        EB_NEW(object_ptr->enc_dec_segment_ctrl[tile_idx],
                enc_dec_segments_ctor,
                init_data_ptr->enc_dec_segment_col,
                init_data_ptr->enc_dec_segment_row);
    }

    // Entropy Rows
    EB_CREATE_MUTEX(object_ptr->entropy_coding_pic_mutex);

    EB_CREATE_MUTEX(object_ptr->intra_mutex);

    EB_CREATE_MUTEX(object_ptr->cdef_search_mutex);

    //object_ptr->mse_seg[0] = (uint64_t(*)[64])eb_aom_malloc(sizeof(**object_ptr->mse_seg) *  picture_sb_width * picture_sb_height);
    // object_ptr->mse_seg[1] = (uint64_t(*)[64])eb_aom_malloc(sizeof(**object_ptr->mse_seg) *  picture_sb_width * picture_sb_height);

    EB_MALLOC_ARRAY(object_ptr->mse_seg[0], picture_sb_width * picture_sb_height);
    EB_MALLOC_ARRAY(object_ptr->mse_seg[1], picture_sb_width * picture_sb_height);

    EB_CREATE_MUTEX(object_ptr->rest_search_mutex);

    //the granularity is 4x4
    EB_MALLOC_ARRAY(object_ptr->mi_grid_base,
                    all_sb * (init_data_ptr->sb_size_pix >> MI_SIZE_LOG2) *
                        (init_data_ptr->sb_size_pix >> MI_SIZE_LOG2));

    EB_MALLOC_ARRAY(object_ptr->mip,
                    all_sb * (init_data_ptr->sb_size_pix >> MI_SIZE_LOG2) *
                        (init_data_ptr->sb_size_pix >> MI_SIZE_LOG2));

    memset(object_ptr->mip,
           0,
           sizeof(ModeInfo) * all_sb * (init_data_ptr->sb_size_pix >> MI_SIZE_LOG2) *
               (init_data_ptr->sb_size_pix >> MI_SIZE_LOG2));

    uint32_t mi_idx;
    for (mi_idx = 0; mi_idx < all_sb * (init_data_ptr->sb_size_pix >> MI_SIZE_LOG2) *
                                  (init_data_ptr->sb_size_pix >> MI_SIZE_LOG2);
         ++mi_idx)
        object_ptr->mi_grid_base[mi_idx] = object_ptr->mip + mi_idx;
    object_ptr->mi_stride = picture_sb_w * (init_data_ptr->sb_size_pix >> MI_SIZE_LOG2);
    if (init_data_ptr->mfmv) {
        //MFMV: map is 8x8 based.
        uint32_t  mi_rows  = init_data_ptr->picture_height >> MI_SIZE_LOG2;
        const int mem_size = ((mi_rows + MAX_MIB_SIZE) >> 1) * (object_ptr->mi_stride >> 1);

        EB_CALLOC_ALIGNED_ARRAY(object_ptr->tpl_mvs, mem_size);
    }
    object_ptr->hash_table.p_lookup_table = NULL;
    av1_hash_table_create(&object_ptr->hash_table);
    EB_MALLOC_ALIGNED(object_ptr->rst_tmpbuf, RESTORATION_TMPBUF_SIZE);
    return EB_ErrorNone;
}

EbErrorType picture_control_set_creator(EbPtr *object_dbl_ptr, EbPtr object_init_data_ptr) {
    PictureControlSet *obj;

    *object_dbl_ptr = NULL;
    EB_NEW(obj, picture_control_set_ctor, object_init_data_ptr);
    *object_dbl_ptr = obj;

    return EB_ErrorNone;
}

static void picture_parent_control_set_dctor(EbPtr p) {
    PictureParentControlSet *obj = (PictureParentControlSet *)p;
    uint32_t                 region_in_picture_width_index;
    uint32_t                 region_in_picture_height_index;

    EB_DELETE(obj->denoise_and_model);
#if !DECOUPLE_ME_RES
    EB_DELETE_PTR_ARRAY(obj->me_results, obj->sb_total_count_unscaled);
#endif
    if (obj->is_chroma_downsampled_picture_ptr_owner)
        EB_DELETE(obj->chroma_downsampled_picture_ptr);

    EB_FREE_2D(obj->variance);
    EB_FREE_2D(obj->y_mean);
    EB_FREE_2D(obj->cb_mean);
    EB_FREE_2D(obj->cr_mean);

    if (obj->picture_histogram) {
        for (region_in_picture_width_index = 0;
             region_in_picture_width_index < MAX_NUMBER_OF_REGIONS_IN_WIDTH;
             region_in_picture_width_index++) {
            if (obj->picture_histogram[region_in_picture_width_index]) {
                for (region_in_picture_height_index = 0;
                     region_in_picture_height_index < MAX_NUMBER_OF_REGIONS_IN_HEIGHT;
                     region_in_picture_height_index++) {
                    EB_FREE_2D(obj->picture_histogram[region_in_picture_width_index]
                                                     [region_in_picture_height_index]);
                }
            }
            EB_FREE_PTR_ARRAY(obj->picture_histogram[region_in_picture_width_index],
                              MAX_NUMBER_OF_REGIONS_IN_HEIGHT);
        }
        EB_FREE_PTR_ARRAY(obj->picture_histogram, MAX_NUMBER_OF_REGIONS_IN_WIDTH);
    }
#if !REMOVE_UNUSED_CODE_PH2

    EB_FREE_2D(obj->ois_sb_results);
#endif
#if TPL_LA
    if (obj->ois_mb_results)
        EB_FREE_2D(obj->ois_mb_results);
    if (obj->tpl_stats)
        EB_FREE_2D(obj->tpl_stats);
    if (obj->tpl_beta)
        EB_FREE_ARRAY(obj->tpl_beta);
#if TPL_LA_LAMBDA_SCALING
    if (obj->tpl_rdmult_scaling_factors)
        EB_FREE_ARRAY(obj->tpl_rdmult_scaling_factors);
    if (obj->tpl_sb_rdmult_scaling_factors)
        EB_FREE_ARRAY(obj->tpl_sb_rdmult_scaling_factors);
#endif
#endif
#if !REMOVE_UNUSED_CODE_PH2
    EB_FREE_2D(obj->ois_candicate);
#endif
    EB_FREE_ARRAY(obj->rc_me_distortion);
    // ME and OIS Distortion Histograms
    EB_FREE_ARRAY(obj->me_distortion_histogram);
    EB_FREE_ARRAY(obj->ois_distortion_histogram);
    EB_FREE_ARRAY(obj->intra_sad_interval_index);
    EB_FREE_ARRAY(obj->inter_sad_interval_index);
    // Non moving index array
    EB_FREE_ARRAY(obj->non_moving_index_array);
    // SB noise variance array
    EB_FREE_ARRAY(obj->sb_flat_noise_array);
    EB_FREE_ARRAY(obj->sb_depth_mode_array);

    if (obj->av1_cm) {
        const int32_t num_planes = 3; // av1_num_planes(cm);
        for (int32_t p = 0; p < num_planes; ++p) {
            RestorationInfo *            ri         = obj->av1_cm->rst_info + p;
            RestorationStripeBoundaries *boundaries = &ri->boundaries;
            EB_FREE_ARRAY(ri->unit_info);
            EB_FREE(boundaries->stripe_boundary_above);
            EB_FREE(boundaries->stripe_boundary_below);
        }
        EB_FREE_ARRAY(obj->av1_cm->frame_to_show);
        if (obj->av1_cm->rst_frame.buffer_alloc_sz) {
            EB_FREE_ARRAY(obj->av1_cm->rst_frame.buffer_alloc);
        }
        EB_FREE_ARRAY(obj->av1_cm);
    }
    EB_FREE_ARRAY(obj->rusi_picture[0]);
    EB_FREE_ARRAY(obj->rusi_picture[1]);
    EB_FREE_ARRAY(obj->rusi_picture[2]);

    EB_FREE_ARRAY(obj->av1x);

    EB_DESTROY_MUTEX(obj->rc_distortion_histogram_mutex);
    EB_DESTROY_SEMAPHORE(obj->temp_filt_done_semaphore);
    EB_DESTROY_MUTEX(obj->temp_filt_mutex);
    EB_DESTROY_MUTEX(obj->debug_mutex);
    EB_FREE_ARRAY(obj->tile_group_info);
    if(obj->frame_superres_enabled){
        eb_pcs_sb_structs_dctor(obj);
        EB_DELETE(obj->enhanced_picture_ptr);
    }
}
EbErrorType picture_parent_control_set_ctor(PictureParentControlSet *object_ptr,
                                            EbPtr                    object_init_data_ptr) {
    PictureControlSetInitData *init_data_ptr    = (PictureControlSetInitData *)object_init_data_ptr;
    EbErrorType                return_error     = EB_ErrorNone;
    const uint16_t             picture_sb_width = (uint16_t)(
        (init_data_ptr->picture_width + init_data_ptr->sb_sz - 1) / init_data_ptr->sb_sz);
    const uint16_t picture_sb_height = (uint16_t)(
        (init_data_ptr->picture_height + init_data_ptr->sb_sz - 1) / init_data_ptr->sb_sz);
    const uint16_t subsampling_x = (init_data_ptr->color_format == EB_YUV444 ? 1 : 2) - 1;
    const uint16_t subsampling_y = (init_data_ptr->color_format >= EB_YUV422 ? 1 : 2) - 1;
#if !REMOVE_UNUSED_CODE_PH2 || !DECOUPLE_ME_RES
    uint16_t       sb_index;
#endif
    uint32_t       region_in_picture_width_index;
    uint32_t       region_in_picture_height_index;

    object_ptr->dctor = picture_parent_control_set_dctor;

    object_ptr->scs_wrapper_ptr               = (EbObjectWrapper *)NULL;
    object_ptr->input_picture_wrapper_ptr     = (EbObjectWrapper *)NULL;
    object_ptr->reference_picture_wrapper_ptr = (EbObjectWrapper *)NULL;
    object_ptr->enhanced_picture_ptr          = (EbPictureBufferDesc *)NULL;
    object_ptr->enhanced_downscaled_picture_ptr = (EbPictureBufferDesc *)NULL;
    object_ptr->enhanced_unscaled_picture_ptr   = (EbPictureBufferDesc *)NULL;

    if (init_data_ptr->color_format >= EB_YUV422) {
        EbPictureBufferDescInitData input_pic_buf_desc_init_data;
        input_pic_buf_desc_init_data.max_width          = init_data_ptr->picture_width;
        input_pic_buf_desc_init_data.max_height         = init_data_ptr->picture_height;
        input_pic_buf_desc_init_data.bit_depth          = 8; //Should be 8bit
        input_pic_buf_desc_init_data.buffer_enable_mask = PICTURE_BUFFER_DESC_CHROMA_MASK;
        input_pic_buf_desc_init_data.left_padding       = init_data_ptr->left_padding;
        input_pic_buf_desc_init_data.right_padding      = init_data_ptr->right_padding;
        input_pic_buf_desc_init_data.top_padding        = init_data_ptr->top_padding;
        input_pic_buf_desc_init_data.bot_padding        = init_data_ptr->bot_padding;
        input_pic_buf_desc_init_data.color_format       = EB_YUV420; //set to 420 for MD
        input_pic_buf_desc_init_data.split_mode         = EB_FALSE;
        EB_NEW(object_ptr->chroma_downsampled_picture_ptr,
               eb_picture_buffer_desc_ctor,
               (EbPtr)&input_pic_buf_desc_init_data);
        object_ptr->is_chroma_downsampled_picture_ptr_owner = EB_TRUE;
    } else if (init_data_ptr->color_format == EB_YUV420) {
        object_ptr->chroma_downsampled_picture_ptr = NULL;
    } else
        return EB_ErrorBadParameter;
    // GOP
    object_ptr->pred_struct_index    = 0;
    object_ptr->picture_number       = 0;
    object_ptr->idr_flag             = EB_FALSE;
    object_ptr->temporal_layer_index = 0;
    object_ptr->total_num_bits       = 0;
    object_ptr->last_idr_picture     = 0;
    object_ptr->sb_total_count       = picture_sb_width * picture_sb_height;
    object_ptr->sb_total_count_unscaled = object_ptr->sb_total_count;

    object_ptr->data_ll_head_ptr         = (EbLinkedListNode *)NULL;
    object_ptr->app_out_data_ll_head_ptr = (EbLinkedListNode *)NULL;
    EB_MALLOC_2D(object_ptr->variance, object_ptr->sb_total_count, MAX_ME_PU_COUNT);
    EB_MALLOC_2D(object_ptr->y_mean, object_ptr->sb_total_count, MAX_ME_PU_COUNT);
    EB_MALLOC_2D(object_ptr->cb_mean, object_ptr->sb_total_count, 21);
    EB_MALLOC_2D(object_ptr->cr_mean, object_ptr->sb_total_count, 21);

    EB_ALLOC_PTR_ARRAY(object_ptr->picture_histogram, MAX_NUMBER_OF_REGIONS_IN_WIDTH);

    for (region_in_picture_width_index = 0;
         region_in_picture_width_index < MAX_NUMBER_OF_REGIONS_IN_WIDTH;
         region_in_picture_width_index++) { // loop over horizontal regions
        EB_ALLOC_PTR_ARRAY(object_ptr->picture_histogram[region_in_picture_width_index],
                           MAX_NUMBER_OF_REGIONS_IN_HEIGHT);
        for (region_in_picture_height_index = 0;
             region_in_picture_height_index < MAX_NUMBER_OF_REGIONS_IN_HEIGHT;
             region_in_picture_height_index++) {
            EB_MALLOC_2D(object_ptr->picture_histogram[region_in_picture_width_index]
                                                      [region_in_picture_height_index],
                         3,
                         HISTOGRAM_NUMBER_OF_BINS);
        }
    }
#if !REMOVE_UNUSED_CODE_PH2
    if (init_data_ptr->allocate_ois_struct) {
        EB_MALLOC_2D(object_ptr->ois_sb_results, object_ptr->sb_total_count, 1);
        EB_MALLOC_2D(
            object_ptr->ois_candicate, object_ptr->sb_total_count, MAX_OIS_CANDIDATES * CU_MAX_COUNT);

        for (sb_index = 0; sb_index < object_ptr->sb_total_count; ++sb_index) {
            uint32_t cu_idx;
            for (cu_idx = 0; cu_idx < CU_MAX_COUNT; ++cu_idx)
                object_ptr->ois_sb_results[sb_index]->ois_candidate_array[cu_idx] =
                &object_ptr->ois_candicate[sb_index][cu_idx * MAX_OIS_CANDIDATES];
        }
    }
    else {
        object_ptr->ois_sb_results = NULL;
    }
#endif
#if TPL_LA
    if(init_data_ptr->enable_tpl_la) {
        const uint16_t picture_width_in_mb  = (uint16_t)((init_data_ptr->picture_width + 15) / 16);
        const uint16_t picture_height_in_mb = (uint16_t)((init_data_ptr->picture_height + 15) / 16);
        object_ptr->r0 = 0;
        object_ptr->is_720p_or_larger = AOMMIN(init_data_ptr->picture_width, init_data_ptr->picture_height) >= 720;
        EB_MALLOC_2D(object_ptr->ois_mb_results, (uint32_t)(picture_width_in_mb * picture_height_in_mb), 1);
        EB_MALLOC_2D(object_ptr->tpl_stats, (uint32_t)((picture_width_in_mb << (1 - object_ptr->is_720p_or_larger)) * (picture_height_in_mb << (1 - object_ptr->is_720p_or_larger))), 1);
        EB_MALLOC_ARRAY(object_ptr->tpl_beta, object_ptr->sb_total_count);
#if TPL_LA_LAMBDA_SCALING
        EB_MALLOC_ARRAY(object_ptr->tpl_rdmult_scaling_factors, picture_width_in_mb * picture_height_in_mb);
        EB_MALLOC_ARRAY(object_ptr->tpl_sb_rdmult_scaling_factors, picture_width_in_mb * picture_height_in_mb);
#endif
    } else {
        object_ptr->r0 = 0;
        object_ptr->is_720p_or_larger = 0;
        object_ptr->ois_mb_results = NULL;
        object_ptr->tpl_stats = NULL;
        object_ptr->tpl_beta = NULL;
#if TPL_LA_LAMBDA_SCALING
        object_ptr->tpl_rdmult_scaling_factors = NULL;
        object_ptr->tpl_sb_rdmult_scaling_factors = NULL;
#endif
    }
#endif
#if !REMOVE_MRP_MODE
    object_ptr->max_number_of_candidates_per_block =
        (init_data_ptr->mrp_mode == 0)
            ? ME_RES_CAND_MRP_MODE_0
            : // [Single Ref = 7] + [BiDir = 12 = 3*4 ] + [UniDir = 4 = 3+1]
            ME_RES_CAND_MRP_MODE_1; // [BiDir = 1] + [UniDir = 2 = 1 + 1]
#endif
#if ! DECOUPLE_ME_RES
    EB_ALLOC_PTR_ARRAY(object_ptr->me_results, object_ptr->sb_total_count);

    for (sb_index = 0; sb_index < object_ptr->sb_total_count; ++sb_index) {
#if REMOVE_MRP_MODE
        EB_NEW(object_ptr->me_results[sb_index],
            me_sb_results_ctor);
#elif NSQ_REMOVAL_CODE_CLEAN_UP
        EB_NEW(object_ptr->me_results[sb_index],
            me_sb_results_ctor,
            init_data_ptr->mrp_mode,
            object_ptr->max_number_of_candidates_per_block);
#else
        EB_NEW(object_ptr->me_results[sb_index],
               me_sb_results_ctor,
               (init_data_ptr->nsq_present) ? MAX_ME_PU_COUNT : SQUARE_PU_COUNT,
               init_data_ptr->mrp_mode,
               object_ptr->max_number_of_candidates_per_block);
#endif
    }
#endif

    EB_MALLOC_ARRAY(object_ptr->rc_me_distortion, object_ptr->sb_total_count);
    // ME and OIS Distortion Histograms
    EB_MALLOC_ARRAY(object_ptr->me_distortion_histogram, NUMBER_OF_SAD_INTERVALS);
    EB_MALLOC_ARRAY(object_ptr->ois_distortion_histogram, NUMBER_OF_INTRA_SAD_INTERVALS);
    EB_MALLOC_ARRAY(object_ptr->intra_sad_interval_index, object_ptr->sb_total_count);
    EB_MALLOC_ARRAY(object_ptr->inter_sad_interval_index, object_ptr->sb_total_count);
    // Non moving index array
    EB_MALLOC_ARRAY(object_ptr->non_moving_index_array, object_ptr->sb_total_count);
    // SB noise variance array
    EB_MALLOC_ARRAY(object_ptr->sb_flat_noise_array, object_ptr->sb_total_count);
    EB_CREATE_MUTEX(object_ptr->rc_distortion_histogram_mutex);
    EB_MALLOC_ARRAY(object_ptr->sb_depth_mode_array, object_ptr->sb_total_count);
    EB_CREATE_SEMAPHORE(object_ptr->temp_filt_done_semaphore, 0, 1);
    EB_CREATE_MUTEX(object_ptr->temp_filt_mutex);
    EB_CREATE_MUTEX(object_ptr->debug_mutex);
    EB_MALLOC_ARRAY(object_ptr->av1_cm, 1);

    object_ptr->av1_cm->interp_filter = SWITCHABLE;

    object_ptr->av1_cm->mi_stride = picture_sb_width * (BLOCK_SIZE_64 / 4);

    EB_MALLOC_ARRAY(object_ptr->av1_cm->frame_to_show, 1);

    object_ptr->av1_cm->use_highbitdepth                  = ((init_data_ptr->bit_depth > 8) || (init_data_ptr->is_16bit_pipeline)) ? 1 : 0;
    object_ptr->av1_cm->bit_depth                         = init_data_ptr->bit_depth;
    object_ptr->av1_cm->color_format                      = init_data_ptr->color_format;
    object_ptr->av1_cm->subsampling_x                     = subsampling_x;
    object_ptr->av1_cm->subsampling_y                     = subsampling_y;
    object_ptr->av1_cm->frm_size.frame_width = init_data_ptr->picture_width- init_data_ptr->non_m8_pad_w;
    object_ptr->av1_cm->frm_size.frame_height = init_data_ptr->picture_height- init_data_ptr->non_m8_pad_h;
    object_ptr->av1_cm->frm_size.superres_upscaled_width = init_data_ptr->picture_width - init_data_ptr->non_m8_pad_w;;
    object_ptr->av1_cm->frm_size.superres_upscaled_height = init_data_ptr->picture_height - init_data_ptr->non_m8_pad_h;

    object_ptr->av1_cm->frm_size.superres_denominator     = SCALE_NUMERATOR;

    object_ptr->av1_cm->mi_cols = init_data_ptr->picture_width >> MI_SIZE_LOG2;
    object_ptr->av1_cm->mi_rows = init_data_ptr->picture_height >> MI_SIZE_LOG2;

    object_ptr->av1_cm->byte_alignment = 0;

    set_restoration_unit_size(init_data_ptr->picture_width,
                              init_data_ptr->picture_height,
                              1,
                              1,
                              object_ptr->av1_cm->rst_info);

    return_error = eb_av1_alloc_restoration_buffers(object_ptr->av1_cm);

    memset(&object_ptr->av1_cm->rst_frame, 0, sizeof(Yv12BufferConfig));

    int32_t ntiles[2];
    for (int32_t is_uv = 0; is_uv < 2; ++is_uv)
        ntiles[is_uv] =
            object_ptr->av1_cm->rst_info[is_uv].units_per_tile; //CHKN res_tiles_in_plane

    assert(ntiles[1] <= ntiles[0]);

    EB_CALLOC_ARRAY(object_ptr->rusi_picture[0], ntiles[0]);
    EB_CALLOC_ARRAY(object_ptr->rusi_picture[1], ntiles[1]);
    EB_CALLOC_ARRAY(object_ptr->rusi_picture[2], ntiles[1]);

    EB_MALLOC_ARRAY(object_ptr->av1x, 1);

    // Film grain noise model if film grain is applied
    if (init_data_ptr->film_grain_noise_level) {
        DenoiseAndModelInitData fg_init_data;
        fg_init_data.encoder_bit_depth    = init_data_ptr->bit_depth;
        fg_init_data.encoder_color_format = init_data_ptr->color_format;
        fg_init_data.noise_level          = init_data_ptr->film_grain_noise_level;
        fg_init_data.width                = init_data_ptr->picture_width;
        fg_init_data.height               = init_data_ptr->picture_height;
        fg_init_data.stride_y = init_data_ptr->picture_width + init_data_ptr->left_padding +
                                init_data_ptr->right_padding;
        fg_init_data.stride_cb = fg_init_data.stride_cr = fg_init_data.stride_y >> subsampling_x;

        EB_NEW(object_ptr->denoise_and_model, denoise_and_model_ctor, (EbPtr)&fg_init_data);
    }

    //Jing: need to know the tile split info at pcs initialize stage
    object_ptr->log2_tile_rows = init_data_ptr->log2_tile_rows;
    object_ptr->log2_tile_cols = init_data_ptr->log2_tile_cols;
    object_ptr->log2_sb_sz = init_data_ptr->log2_sb_sz;
    set_tile_info(object_ptr);
    EB_MALLOC_ARRAY(object_ptr->tile_group_info,
            (object_ptr->av1_cm->tiles_info.tile_rows * object_ptr->av1_cm->tiles_info.tile_cols));

    object_ptr->frame_superres_enabled = EB_FALSE;
    object_ptr->aligned_width = init_data_ptr->picture_width;
    object_ptr->aligned_height = init_data_ptr->picture_height;
    object_ptr->frame_width = init_data_ptr->picture_width;
    object_ptr->frame_height = init_data_ptr->picture_height;

    object_ptr->superres_denom = SCALE_NUMERATOR;

    return return_error;
}
#if DECOUPLE_ME_RES
static void me_dctor(EbPtr p) {
    MotionEstimationData *obj = (MotionEstimationData *)p;

    EB_DELETE_PTR_ARRAY(obj->me_results, obj->sb_total_count_unscaled);
}
EbErrorType me_ctor(MotionEstimationData *object_ptr,
    EbPtr                    object_init_data_ptr) {

    PictureControlSetInitData *init_data_ptr = (PictureControlSetInitData *)object_init_data_ptr;
    EbErrorType                return_error = EB_ErrorNone;
    const uint16_t             picture_sb_width = (uint16_t)(
        (init_data_ptr->picture_width + init_data_ptr->sb_sz - 1) / init_data_ptr->sb_sz);
    const uint16_t picture_sb_height = (uint16_t)(
        (init_data_ptr->picture_height + init_data_ptr->sb_sz - 1) / init_data_ptr->sb_sz);

    uint16_t       sb_index;
    object_ptr->dctor = me_dctor;
    uint32_t sb_total_count = picture_sb_width * picture_sb_height;
    object_ptr->sb_total_count_unscaled = sb_total_count;

    EB_ALLOC_PTR_ARRAY(object_ptr->me_results,  sb_total_count);

    for (sb_index = 0; sb_index < sb_total_count; ++sb_index) {
        EB_NEW(object_ptr->me_results[sb_index],
            me_sb_results_ctor);
    }

    return return_error;
}
#endif
EbErrorType sb_params_init_pcs(SequenceControlSet *scs_ptr,
                               PictureParentControlSet *pcs_ptr) {
    EbErrorType return_error = EB_ErrorNone;
    uint16_t sb_index;
    uint16_t raster_scan_blk_index;
    uint16_t encoding_width = pcs_ptr->aligned_width;
    uint16_t encoding_height = pcs_ptr->aligned_height;

    uint8_t picture_sb_width =
            (uint8_t)((encoding_width + scs_ptr->sb_sz - 1) / scs_ptr->sb_sz);
    uint8_t picture_sb_height =
            (uint8_t)((encoding_height + scs_ptr->sb_sz - 1) / scs_ptr->sb_sz);

    EB_MALLOC_ARRAY(pcs_ptr->sb_params_array, picture_sb_width * picture_sb_height);

    for (sb_index = 0; sb_index < picture_sb_width * picture_sb_height; ++sb_index) {
        pcs_ptr->sb_params_array[sb_index].horizontal_index =
                (uint8_t)(sb_index % picture_sb_width);
        pcs_ptr->sb_params_array[sb_index].vertical_index = (uint8_t)(sb_index / picture_sb_width);
        pcs_ptr->sb_params_array[sb_index].origin_x =
                pcs_ptr->sb_params_array[sb_index].horizontal_index * scs_ptr->sb_sz;
        pcs_ptr->sb_params_array[sb_index].origin_y =
                pcs_ptr->sb_params_array[sb_index].vertical_index * scs_ptr->sb_sz;

        pcs_ptr->sb_params_array[sb_index].width = (uint8_t)(
                ((encoding_width - pcs_ptr->sb_params_array[sb_index].origin_x) <
                 scs_ptr->sb_sz)
                ? encoding_width - pcs_ptr->sb_params_array[sb_index].origin_x
                : scs_ptr->sb_sz);

        pcs_ptr->sb_params_array[sb_index].height = (uint8_t)(
                ((encoding_height - pcs_ptr->sb_params_array[sb_index].origin_y) <
                 scs_ptr->sb_sz)
                ? encoding_height - pcs_ptr->sb_params_array[sb_index].origin_y
                : scs_ptr->sb_sz);

        pcs_ptr->sb_params_array[sb_index].is_complete_sb =
                (uint8_t)(((pcs_ptr->sb_params_array[sb_index].width == scs_ptr->sb_sz) &&
                           (pcs_ptr->sb_params_array[sb_index].height == scs_ptr->sb_sz))
                          ? 1
                          : 0);

        pcs_ptr->sb_params_array[sb_index].is_edge_sb =
                (pcs_ptr->sb_params_array[sb_index].origin_x < scs_ptr->sb_sz) ||
                (pcs_ptr->sb_params_array[sb_index].origin_y < scs_ptr->sb_sz) ||
                (pcs_ptr->sb_params_array[sb_index].origin_x >
                 encoding_width - scs_ptr->sb_sz) ||
                (pcs_ptr->sb_params_array[sb_index].origin_y >
                 encoding_height - scs_ptr->sb_sz)
                ? 1
                : 0;

        for (raster_scan_blk_index = RASTER_SCAN_CU_INDEX_64x64;
             raster_scan_blk_index <= RASTER_SCAN_CU_INDEX_8x8_63;
             raster_scan_blk_index++) {
            pcs_ptr->sb_params_array[sb_index].raster_scan_blk_validity[raster_scan_blk_index] =
                    ((pcs_ptr->sb_params_array[sb_index].origin_x +
                      raster_scan_blk_x[raster_scan_blk_index] +
                      raster_scan_blk_size[raster_scan_blk_index] >
                      encoding_width) ||
                     (pcs_ptr->sb_params_array[sb_index].origin_y +
                      raster_scan_blk_y[raster_scan_blk_index] +
                      raster_scan_blk_size[raster_scan_blk_index] >
                      encoding_height))
                    ? EB_FALSE
                    : EB_TRUE;
        }
    }

    return return_error;
}

EbErrorType sb_geom_init_pcs(SequenceControlSet *scs_ptr, PictureParentControlSet *pcs_ptr) {
    uint16_t sb_index;
    uint16_t md_scan_block_index;

    uint16_t encoding_width = pcs_ptr->aligned_width;
    uint16_t encoding_height = pcs_ptr->aligned_height;

    uint16_t picture_sb_width =
            (encoding_width + scs_ptr->sb_size_pix - 1) / scs_ptr->sb_size_pix;
    uint16_t picture_sb_height =
            (encoding_height + scs_ptr->sb_size_pix - 1) / scs_ptr->sb_size_pix;

    EB_MALLOC_ARRAY(pcs_ptr->sb_geom, picture_sb_width * picture_sb_height);

    for (sb_index = 0; sb_index < picture_sb_width * picture_sb_height; ++sb_index) {
        pcs_ptr->sb_geom[sb_index].horizontal_index = sb_index % picture_sb_width;
        pcs_ptr->sb_geom[sb_index].vertical_index   = sb_index / picture_sb_width;
        pcs_ptr->sb_geom[sb_index].origin_x =
                pcs_ptr->sb_geom[sb_index].horizontal_index * scs_ptr->sb_size_pix;
        pcs_ptr->sb_geom[sb_index].origin_y =
                pcs_ptr->sb_geom[sb_index].vertical_index * scs_ptr->sb_size_pix;

        pcs_ptr->sb_geom[sb_index].width = (uint8_t)(
                ((encoding_width - pcs_ptr->sb_geom[sb_index].origin_x) <
                 scs_ptr->sb_size_pix)
                ? encoding_width - pcs_ptr->sb_geom[sb_index].origin_x
                : scs_ptr->sb_size_pix);

        pcs_ptr->sb_geom[sb_index].height = (uint8_t)(
                ((encoding_height - pcs_ptr->sb_geom[sb_index].origin_y) <
                 scs_ptr->sb_size_pix)
                ? encoding_height - pcs_ptr->sb_geom[sb_index].origin_y
                : scs_ptr->sb_size_pix);

        pcs_ptr->sb_geom[sb_index].is_complete_sb =
                (uint8_t)(((pcs_ptr->sb_geom[sb_index].width == scs_ptr->sb_size_pix) &&
                           (pcs_ptr->sb_geom[sb_index].height == scs_ptr->sb_size_pix))
                          ? 1
                          : 0);

        uint16_t max_block_count = scs_ptr->max_block_cnt;

        for (md_scan_block_index = 0; md_scan_block_index < max_block_count;
             md_scan_block_index++) {
            const BlockGeom *blk_geom = get_blk_geom_mds(md_scan_block_index);
            if (scs_ptr->over_boundary_block_mode == 1) {
                pcs_ptr->sb_geom[sb_index].block_is_allowed[md_scan_block_index] =
                        ((pcs_ptr->sb_geom[sb_index].origin_x + blk_geom->origin_x +
                          blk_geom->bwidth / 2 <
                          encoding_width) &&
                         (pcs_ptr->sb_geom[sb_index].origin_y + blk_geom->origin_y +
                          blk_geom->bheight / 2 <
                          encoding_height))
                        ? EB_TRUE
                        : EB_FALSE;
#if !INCOMPLETE_SB_FIX
                // Temporary if the cropped width is not 4, 8, 16, 32, 64 and 128, the block is not allowed. To be removed after intrinsic functions for NxM spatial_full_distortion_kernel_func_ptr_array are added
                int32_t cropped_width =
                        MIN(blk_geom->bwidth,
                            encoding_width -
                            (pcs_ptr->sb_geom[sb_index].origin_x + blk_geom->origin_x));
                if (cropped_width != 4 && cropped_width != 8 && cropped_width != 16 &&
                    cropped_width != 32 && cropped_width != 64 && cropped_width != 128)
                    pcs_ptr->sb_geom[sb_index].block_is_allowed[md_scan_block_index] = EB_FALSE;
#endif

                if (blk_geom->shape != PART_N) blk_geom = get_blk_geom_mds(blk_geom->sqi_mds);
                pcs_ptr->sb_geom[sb_index].block_is_inside_md_scan[md_scan_block_index] =
                        ((pcs_ptr->sb_geom[sb_index].origin_x >= encoding_width) ||
                         (pcs_ptr->sb_geom[sb_index].origin_y >= encoding_height))
                        ? EB_FALSE
                        : EB_TRUE;
            } else {
                if (blk_geom->shape != PART_N) blk_geom = get_blk_geom_mds(blk_geom->sqi_mds);

                pcs_ptr->sb_geom[sb_index].block_is_allowed[md_scan_block_index] =
                        ((pcs_ptr->sb_geom[sb_index].origin_x + blk_geom->origin_x + blk_geom->bwidth >
                          encoding_width) ||
                         (pcs_ptr->sb_geom[sb_index].origin_y + blk_geom->origin_y + blk_geom->bheight >
                          encoding_height))
                        ? EB_FALSE
                        : EB_TRUE;

                pcs_ptr->sb_geom[sb_index].block_is_inside_md_scan[md_scan_block_index] =
                        ((pcs_ptr->sb_geom[sb_index].origin_x + blk_geom->origin_x + blk_geom->bwidth >
                          encoding_width) ||
                         (pcs_ptr->sb_geom[sb_index].origin_y + blk_geom->origin_y + blk_geom->bheight >
                          encoding_height))
                        ? EB_FALSE
                        : EB_TRUE;
            }
        }
    }

    return 0;
}

EbErrorType picture_parent_control_set_creator(EbPtr *object_dbl_ptr, EbPtr object_init_data_ptr) {
    PictureParentControlSet *obj;

    *object_dbl_ptr = NULL;
    EB_NEW(obj, picture_parent_control_set_ctor, object_init_data_ptr);
    *object_dbl_ptr = obj;

    return EB_ErrorNone;
}
#if DECOUPLE_ME_RES
EbErrorType me_creator(EbPtr *object_dbl_ptr, EbPtr object_init_data_ptr) {
    MotionEstimationData *obj;

    *object_dbl_ptr = NULL;
    EB_NEW(obj, me_ctor, object_init_data_ptr);
    *object_dbl_ptr = obj;

    return EB_ErrorNone;
}
#endif
