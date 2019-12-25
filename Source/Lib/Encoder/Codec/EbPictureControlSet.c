// clang-format off
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
#include "EbPictureBufferDesc.h"

void *eb_aom_memalign(size_t align, size_t size);
void eb_aom_free(void *memblk);
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

void segmentation_map_dctor(EbPtr p)
{
    SegmentationNeighborMap *obj = (SegmentationNeighborMap*)p;
    EB_FREE_ARRAY(obj->data);
}

EbErrorType segmentation_map_ctor(SegmentationNeighborMap *seg_neighbor_map,
                                uint16_t pic_width, uint16_t pic_height){

    uint32_t num_elements = (pic_width >> MI_SIZE_LOG2) * (pic_height >> MI_SIZE_LOG2);

    seg_neighbor_map->dctor = segmentation_map_dctor;

    seg_neighbor_map->map_size = num_elements;
    EB_CALLOC_ARRAY(seg_neighbor_map->data, num_elements);
    return  EB_ErrorNone;
}

static void me_sb_results_dctor(EbPtr p)
{
    MeLcuResults* obj = (MeLcuResults*)p;

    EB_FREE_ARRAY(obj->me_candidate);
    if (obj->me_mv_array) {
        EB_FREE_ARRAY(obj->me_mv_array[0]);
    }
    EB_FREE_ARRAY(obj->me_mv_array);
    EB_FREE_ARRAY(obj->me_candidate_array);
    EB_FREE_ARRAY(obj->total_me_candidate_index);

    EB_FREE_ARRAY(obj->me_nsq_0);
    EB_FREE_ARRAY(obj->me_nsq_1);
}

EbErrorType me_sb_results_ctor(
    MeLcuResults      *objectPtr,
    uint32_t           maxNumberOfPusPerLcu,
    uint8_t            mrp_mode,
    uint32_t           maxNumberOfMeCandidatesPerPU){
    uint32_t  puIndex;

    size_t count = ((mrp_mode == 0) ? ME_MV_MRP_MODE_0 : ME_MV_MRP_MODE_1);
    objectPtr->dctor = me_sb_results_dctor;
    objectPtr->max_number_of_pus_per_lcu = maxNumberOfPusPerLcu;

    EB_MALLOC_ARRAY(objectPtr->me_candidate, maxNumberOfPusPerLcu);
    EB_MALLOC_ARRAY(objectPtr->me_mv_array, maxNumberOfPusPerLcu);
    EB_MALLOC_ARRAY(objectPtr->me_candidate_array, maxNumberOfPusPerLcu * maxNumberOfMeCandidatesPerPU);
    EB_MALLOC_ARRAY(objectPtr->me_mv_array[0], maxNumberOfPusPerLcu * count);

    for (puIndex = 0; puIndex < maxNumberOfPusPerLcu; ++puIndex) {
        objectPtr->me_candidate[puIndex] = &objectPtr->me_candidate_array[puIndex * maxNumberOfMeCandidatesPerPU];

        objectPtr->me_candidate[puIndex][0].ref_idx_l0 = 0;
        objectPtr->me_candidate[puIndex][0].ref_idx_l1 = 0;
        objectPtr->me_candidate[puIndex][1].ref_idx_l0 = 0;
        objectPtr->me_candidate[puIndex][1].ref_idx_l1 = 0;
        objectPtr->me_candidate[puIndex][2].ref_idx_l0 = 0;
        objectPtr->me_candidate[puIndex][2].ref_idx_l1 = 0;

        objectPtr->me_candidate[puIndex][0].direction = 0;
        objectPtr->me_candidate[puIndex][1].direction = 1;
        objectPtr->me_candidate[puIndex][2].direction = 2;
        objectPtr->me_mv_array[puIndex] = objectPtr->me_mv_array[0] + puIndex * count;
    }
    EB_MALLOC_ARRAY(objectPtr->total_me_candidate_index, maxNumberOfPusPerLcu);

    EB_MALLOC_ARRAY(objectPtr->me_nsq_0, maxNumberOfPusPerLcu);
    EB_MALLOC_ARRAY(objectPtr->me_nsq_1, maxNumberOfPusPerLcu);

    //objectPtr->lcuDistortion = 0;
    return EB_ErrorNone;
}

void picture_control_set_dctor(EbPtr p)
{
    PictureControlSet* obj = (PictureControlSet*)p;
    uint8_t depth;
    av1_hash_table_destroy(&obj->hash_table);
    EB_FREE_ALIGNED_ARRAY(obj->tpl_mvs);
    EB_DELETE(obj->enc_dec_segment_ctrl);
    EB_DELETE(obj->ep_intra_luma_mode_neighbor_array);
    EB_DELETE(obj->ep_intra_chroma_mode_neighbor_array);
    EB_DELETE(obj->ep_mv_neighbor_array);
    EB_DELETE(obj->ep_skip_flag_neighbor_array);
    EB_DELETE(obj->ep_mode_type_neighbor_array);
    EB_DELETE(obj->ep_leaf_depth_neighbor_array);
    EB_DELETE(obj->ep_luma_recon_neighbor_array);
    EB_DELETE(obj->ep_cb_recon_neighbor_array);
    EB_DELETE(obj->ep_cr_recon_neighbor_array);
    EB_DELETE(obj->ep_luma_dc_sign_level_coeff_neighbor_array);
    EB_DELETE(obj->ep_cb_dc_sign_level_coeff_neighbor_array);
    EB_DELETE(obj->ep_cr_dc_sign_level_coeff_neighbor_array);
    EB_DELETE(obj->ep_partition_context_neighbor_array);
    EB_DELETE(obj->mode_type_neighbor_array);
    EB_DELETE(obj->partition_context_neighbor_array);
    EB_DELETE(obj->skip_flag_neighbor_array);
    EB_DELETE(obj->skip_coeff_neighbor_array);
    EB_DELETE(obj->luma_dc_sign_level_coeff_neighbor_array);
    EB_DELETE(obj->cr_dc_sign_level_coeff_neighbor_array);
    EB_DELETE(obj->cb_dc_sign_level_coeff_neighbor_array);
    EB_DELETE(obj->inter_pred_dir_neighbor_array);
    EB_DELETE(obj->ref_frame_type_neighbor_array);
    EB_DELETE(obj->intra_luma_mode_neighbor_array);
    EB_DELETE(obj->txfm_context_array);
    EB_DELETE(obj->segmentation_id_pred_array);
    EB_DELETE(obj->segmentation_neighbor_map);
    EB_DELETE(obj->ep_luma_recon_neighbor_array16bit);
    EB_DELETE(obj->ep_cb_recon_neighbor_array16bit);
    EB_DELETE(obj->ep_cr_recon_neighbor_array16bit);
    EB_DELETE(obj->interpolation_type_neighbor_array);

    for (depth = 0; depth < NEIGHBOR_ARRAY_TOTAL_COUNT; depth++) {
        EB_DELETE(obj->md_intra_luma_mode_neighbor_array[depth]);
        EB_DELETE(obj->md_intra_chroma_mode_neighbor_array[depth]);
        EB_DELETE(obj->md_mv_neighbor_array[depth]);
        EB_DELETE(obj->md_skip_flag_neighbor_array[depth]);
        EB_DELETE(obj->md_mode_type_neighbor_array[depth]);
        EB_DELETE(obj->md_leaf_depth_neighbor_array[depth]);
        EB_DELETE(obj->mdleaf_partition_neighbor_array[depth]);
        if (obj->hbd_mode_decision > EB_8_BIT_MD){
            EB_DELETE(obj->md_luma_recon_neighbor_array16bit[depth]);
            EB_DELETE(obj->md_tx_depth_1_luma_recon_neighbor_array16bit[depth]);
            EB_DELETE(obj->md_cb_recon_neighbor_array16bit[depth]);
            EB_DELETE(obj->md_cr_recon_neighbor_array16bit[depth]);
        }
        if (obj->hbd_mode_decision != EB_10_BIT_MD){
            EB_DELETE(obj->md_luma_recon_neighbor_array[depth]);
            EB_DELETE(obj->md_tx_depth_1_luma_recon_neighbor_array[depth]);
            EB_DELETE(obj->md_cb_recon_neighbor_array[depth]);
            EB_DELETE(obj->md_cr_recon_neighbor_array[depth]);
        }

        EB_DELETE(obj->md_skip_coeff_neighbor_array[depth]);
        EB_DELETE(obj->md_luma_dc_sign_level_coeff_neighbor_array[depth]);
        EB_DELETE(obj->md_tx_depth_1_luma_dc_sign_level_coeff_neighbor_array[depth]);
        EB_DELETE(obj->md_cr_dc_sign_level_coeff_neighbor_array[depth]);
        EB_DELETE(obj->md_cb_dc_sign_level_coeff_neighbor_array[depth]);
        EB_DELETE(obj->md_txfm_context_array[depth]);
        EB_DELETE(obj->md_inter_pred_dir_neighbor_array[depth]);
        EB_DELETE(obj->md_ref_frame_type_neighbor_array[depth]);
        EB_DELETE(obj->md_interpolation_type_neighbor_array[depth]);
    }
    EB_DELETE_PTR_ARRAY(obj->sb_ptr_array, obj->sb_total_count);
    EB_DELETE(obj->coeff_est_entropy_coder_ptr);
    EB_DELETE(obj->bitstream_ptr);
    EB_DELETE(obj->entropy_coder_ptr);
    EB_DELETE(obj->recon_picture32bit_ptr);
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
    EB_FREE_ARRAY(obj->rate_est_array);
    if(obj->tile_tok[0][0])
       EB_FREE_ARRAY(obj->tile_tok[0][0]);
    EB_FREE_ARRAY(obj->mdc_sb_array);
    EB_FREE_ARRAY(obj->qp_array);
    EB_DESTROY_MUTEX(obj->entropy_coding_mutex);
    EB_DESTROY_MUTEX(obj->intra_mutex);
    EB_DESTROY_MUTEX(obj->cdef_search_mutex);
    EB_DESTROY_MUTEX(obj->rest_search_mutex);

}
// Token buffer is only used for palette tokens.
static INLINE unsigned int get_token_alloc(int mb_rows, int mb_cols,
    int sb_size_log2,
    const int num_planes) {
    // Calculate the maximum number of max superblocks in the image.
    const int shift = sb_size_log2 - 4;
    const int sb_size = 1 << sb_size_log2;
    const int sb_size_square = sb_size * sb_size;
    const int sb_rows = ALIGN_POWER_OF_TWO(mb_rows, shift) >> shift;
    const int sb_cols = ALIGN_POWER_OF_TWO(mb_cols, shift) >> shift;

    // One palette token for each pixel. There can be palettes on two planes.
    const int sb_palette_toks = AOMMIN(2, num_planes) * sb_size_square;

    return sb_rows * sb_cols * sb_palette_toks;
}


typedef struct InitData {
    NeighborArrayUnit **na_unit_dbl_ptr;
    uint32_t   max_picture_width;
    uint32_t   max_picture_height;
    uint32_t   unit_size;
    uint32_t   granularity_normal;
    uint32_t   granularity_top_left;
    uint32_t   type_mask;
} InitData;

#define DIM(array) (sizeof(array) / sizeof(array[0]))
EbErrorType create_neighbor_array_units(InitData* data, size_t count)
{
    for (size_t i = 0; i < count; i++) {
        EB_NEW(
            *data[i].na_unit_dbl_ptr,
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

EbErrorType picture_control_set_ctor(
    PictureControlSet *object_ptr,
    EbPtr object_init_data_ptr)
{
    PictureControlSetInitData *initDataPtr = (PictureControlSetInitData*)object_init_data_ptr;

    EbPictureBufferDescInitData input_picture_buffer_desc_init_data;
    EbPictureBufferDescInitData coeffBufferDescInitData;

    // Max/Min CU Sizes
    const uint32_t maxCuSize = initDataPtr->sb_size_pix;
    // LCUs
    const uint16_t pictureLcuWidth = (uint16_t)((initDataPtr->picture_width + initDataPtr->sb_sz - 1) / initDataPtr->sb_sz);
    const uint16_t pictureLcuHeight = (uint16_t)((initDataPtr->picture_height + initDataPtr->sb_sz - 1) / initDataPtr->sb_sz);
    uint16_t sb_index;
    uint16_t sb_origin_x;
    uint16_t sb_origin_y;
    EbErrorType return_error = EB_ErrorNone;

    EbBool is16bit = initDataPtr->bit_depth > 8 ? EB_TRUE : EB_FALSE;
    const uint16_t subsampling_x = (initDataPtr->color_format == EB_YUV444 ? 1 : 2) - 1;
    const uint16_t subsampling_y = (initDataPtr->color_format >= EB_YUV422 ? 1 : 2) - 1;

    object_ptr->dctor = picture_control_set_dctor;

    // Init Picture Init data
    input_picture_buffer_desc_init_data.max_width = initDataPtr->picture_width;
    input_picture_buffer_desc_init_data.max_height = initDataPtr->picture_height;
    input_picture_buffer_desc_init_data.bit_depth = initDataPtr->bit_depth;
    input_picture_buffer_desc_init_data.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
    input_picture_buffer_desc_init_data.color_format = initDataPtr->color_format;

    input_picture_buffer_desc_init_data.left_padding = PAD_VALUE;
    input_picture_buffer_desc_init_data.right_padding = PAD_VALUE;
    input_picture_buffer_desc_init_data.top_padding = PAD_VALUE;
    input_picture_buffer_desc_init_data.bot_padding = PAD_VALUE;

    input_picture_buffer_desc_init_data.split_mode = EB_FALSE;

    coeffBufferDescInitData.max_width = initDataPtr->picture_width;
    coeffBufferDescInitData.max_height = initDataPtr->picture_height;
    coeffBufferDescInitData.bit_depth = EB_16BIT;
    coeffBufferDescInitData.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
    coeffBufferDescInitData.color_format = initDataPtr->color_format;

    coeffBufferDescInitData.left_padding = PAD_VALUE;
    coeffBufferDescInitData.right_padding = PAD_VALUE;
    coeffBufferDescInitData.top_padding = PAD_VALUE;
    coeffBufferDescInitData.bot_padding = PAD_VALUE;

    coeffBufferDescInitData.split_mode = EB_FALSE;

    object_ptr->sequence_control_set_wrapper_ptr = (EbObjectWrapper *)EB_NULL;

    object_ptr->recon_picture16bit_ptr = (EbPictureBufferDesc *)EB_NULL;
    object_ptr->recon_picture_ptr = (EbPictureBufferDesc *)EB_NULL;
    object_ptr->color_format = initDataPtr->color_format;

    EbPictureBufferDescInitData coeffBufferDes32bitInitData;
    coeffBufferDes32bitInitData.max_width = initDataPtr->picture_width;
    coeffBufferDes32bitInitData.max_height = initDataPtr->picture_height;
    coeffBufferDes32bitInitData.bit_depth = EB_32BIT;
    coeffBufferDes32bitInitData.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
    coeffBufferDes32bitInitData.color_format = initDataPtr->color_format;
    coeffBufferDes32bitInitData.left_padding = 0;
    coeffBufferDes32bitInitData.right_padding = 0;
    coeffBufferDes32bitInitData.top_padding = 0;
    coeffBufferDes32bitInitData.bot_padding = 0;
    coeffBufferDes32bitInitData.split_mode = EB_FALSE;

    object_ptr->recon_picture32bit_ptr = (EbPictureBufferDesc *)EB_NULL;
    EB_NEW(
        object_ptr->recon_picture32bit_ptr,
        eb_recon_picture_buffer_desc_ctor,
        (EbPtr)&coeffBufferDes32bitInitData);

    // Reconstructed Picture Buffer
    if (is16bit) {
        EB_NEW(
            object_ptr->recon_picture16bit_ptr,
            eb_recon_picture_buffer_desc_ctor,
            (EbPtr)&coeffBufferDescInitData);
    }
    else
    {
        EB_NEW(
            object_ptr->recon_picture_ptr,
            eb_recon_picture_buffer_desc_ctor,
            (EbPtr)&input_picture_buffer_desc_init_data);
    }
    // Film Grain Picture Buffer
    if (initDataPtr->film_grain_noise_level) {
        if (is16bit) {
            EB_NEW(
                object_ptr->film_grain_picture16bit_ptr,
                eb_recon_picture_buffer_desc_ctor,
                (EbPtr)&coeffBufferDescInitData);
        }
        else
        {
            EB_NEW(
                object_ptr->film_grain_picture_ptr,
                eb_recon_picture_buffer_desc_ctor,
                (EbPtr)&input_picture_buffer_desc_init_data);
        }
    }

    if (is16bit) {
        EB_NEW(
            object_ptr->input_frame16bit,
            eb_picture_buffer_desc_ctor,
            (EbPtr)&coeffBufferDescInitData);
    }
    // Entropy Coder
    EB_NEW(
        object_ptr->entropy_coder_ptr,
        entropy_coder_ctor,
        SEGMENT_ENTROPY_BUFFER_SIZE);

    // Packetization process Bitstream
    EB_NEW(
        object_ptr->bitstream_ptr,
        bitstream_ctor,
        PACKETIZATION_PROCESS_BUFFER_SIZE);

    if (return_error == EB_ErrorInsufficientResources)
        return EB_ErrorInsufficientResources;
    // Rate estimation entropy coder
    EB_NEW(
        object_ptr->coeff_est_entropy_coder_ptr,
        entropy_coder_ctor,
        SEGMENT_ENTROPY_BUFFER_SIZE);
    // GOP
    object_ptr->picture_number = 0;
    object_ptr->temporal_layer_index = 0;

    // SB Array
    object_ptr->sb_max_depth = (uint8_t)initDataPtr->max_depth;
    object_ptr->sb_total_count = pictureLcuWidth * pictureLcuHeight;
    EB_ALLOC_PTR_ARRAY(object_ptr->sb_ptr_array, object_ptr->sb_total_count);

    sb_origin_x = 0;
    sb_origin_y = 0;

    const uint16_t picture_sb_w   = (uint16_t)((initDataPtr->picture_width  + initDataPtr->sb_size_pix - 1) / initDataPtr->sb_size_pix);
    const uint16_t picture_sb_h   = (uint16_t)((initDataPtr->picture_height + initDataPtr->sb_size_pix - 1) / initDataPtr->sb_size_pix);
    const uint16_t all_sb = picture_sb_w * picture_sb_h;

    for (sb_index = 0; sb_index < all_sb; ++sb_index) {
        EB_NEW(
            object_ptr->sb_ptr_array[sb_index],
            largest_coding_unit_ctor,
            (uint8_t)initDataPtr->sb_size_pix,
            (uint16_t)(sb_origin_x * maxCuSize),
            (uint16_t)(sb_origin_y * maxCuSize),
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
    EB_MALLOC_ARRAY(object_ptr->rate_est_array, all_sb);

    if (initDataPtr->cfg_palette){
        uint32_t mi_cols = initDataPtr->picture_width >> MI_SIZE_LOG2;
        uint32_t mi_rows = initDataPtr->picture_height >> MI_SIZE_LOG2;
        uint32_t mb_cols = (mi_cols + 2) >> 2;
        uint32_t mb_rows = (mi_rows + 2) >> 2;
        unsigned int tokens =
            get_token_alloc(mb_rows, mb_cols, MAX_SB_SIZE_LOG2, 2);
        EB_CALLOC_ARRAY(object_ptr->tile_tok[0][0], tokens);
    }
    else
        object_ptr->tile_tok[0][0] = NULL;
    // Mode Decision Control config
    EB_MALLOC_ARRAY(object_ptr->mdc_sb_array, object_ptr->sb_total_count);
    object_ptr->qp_array_stride = (uint16_t)((initDataPtr->picture_width + MIN_BLOCK_SIZE - 1) / MIN_BLOCK_SIZE);
    object_ptr->qp_array_size = ((initDataPtr->picture_width + MIN_BLOCK_SIZE - 1) / MIN_BLOCK_SIZE) *
        ((initDataPtr->picture_height + MIN_BLOCK_SIZE - 1) / MIN_BLOCK_SIZE);

    // Allocate memory for qp array (used by DLF)
    EB_MALLOC_ARRAY(object_ptr->qp_array, object_ptr->qp_array_size);

    object_ptr->hbd_mode_decision = initDataPtr->hbd_mode_decision;
    // Mode Decision Neighbor Arrays
    uint8_t depth;
    for (depth = 0; depth < NEIGHBOR_ARRAY_TOTAL_COUNT; depth++) {
        InitData data[] = {
            {
                &object_ptr->md_intra_luma_mode_neighbor_array[depth],
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK
            },
            {
                &object_ptr->md_intra_chroma_mode_neighbor_array[depth],
                MAX_PICTURE_WIDTH_SIZE >> subsampling_x,
                MAX_PICTURE_HEIGHT_SIZE >> subsampling_y,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            {
                &object_ptr->md_mv_neighbor_array[depth],
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(MvUnit),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_FULL_MASK,
            },
            {
                &object_ptr->md_skip_flag_neighbor_array[depth],
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            {
                &object_ptr->md_mode_type_neighbor_array[depth],
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_FULL_MASK,
            },
            {
                &object_ptr->md_leaf_depth_neighbor_array[depth],
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            {
                &object_ptr->mdleaf_partition_neighbor_array[depth],
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(struct PartitionContext),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            {
                &object_ptr->md_skip_coeff_neighbor_array[depth],
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            // for each 4x4
            {
                &object_ptr->md_luma_dc_sign_level_coeff_neighbor_array[depth],
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            // for each 4x4
            {
                &object_ptr->md_tx_depth_1_luma_dc_sign_level_coeff_neighbor_array[depth],
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            // for each 4x4
            {
                &object_ptr->md_cr_dc_sign_level_coeff_neighbor_array[depth],
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            // for each 4x4
            {
                &object_ptr->md_cb_dc_sign_level_coeff_neighbor_array[depth],
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            {
                &object_ptr->md_txfm_context_array[depth],
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(TXFM_CONTEXT),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            {
                &object_ptr->md_inter_pred_dir_neighbor_array[depth],
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            {
                &object_ptr->md_ref_frame_type_neighbor_array[depth],
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            }
        };
        return_error = create_neighbor_array_units(data, DIM(data));
        if (return_error == EB_ErrorInsufficientResources)
            return EB_ErrorInsufficientResources;
        if (initDataPtr->hbd_mode_decision != EB_10_BIT_MD) {
            InitData data[] = {

                {
                    &object_ptr->md_luma_recon_neighbor_array[depth],
                    MAX_PICTURE_WIDTH_SIZE,
                    MAX_PICTURE_HEIGHT_SIZE,
                    sizeof(uint8_t),
                    SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                    SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                    NEIGHBOR_ARRAY_UNIT_FULL_MASK,
                },
                {
                    &object_ptr->md_tx_depth_1_luma_recon_neighbor_array[depth],
                    MAX_PICTURE_WIDTH_SIZE,
                    MAX_PICTURE_HEIGHT_SIZE,
                    sizeof(uint8_t),
                    SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                    SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                    NEIGHBOR_ARRAY_UNIT_FULL_MASK,
                },
                {
                    &object_ptr->md_cb_recon_neighbor_array[depth],
                    MAX_PICTURE_WIDTH_SIZE >> subsampling_x,
                    MAX_PICTURE_HEIGHT_SIZE >> subsampling_y,
                    sizeof(uint8_t),
                    SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                    SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                    NEIGHBOR_ARRAY_UNIT_FULL_MASK,
                },
                {
                    &object_ptr->md_cr_recon_neighbor_array[depth],
                    MAX_PICTURE_WIDTH_SIZE >> subsampling_x,
                    MAX_PICTURE_HEIGHT_SIZE >> subsampling_y,
                    sizeof(uint8_t),
                    SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                    SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                    NEIGHBOR_ARRAY_UNIT_FULL_MASK,
                }

            };
            return_error = create_neighbor_array_units(data, DIM(data));
            if (return_error == EB_ErrorInsufficientResources)
                return EB_ErrorInsufficientResources;
        }
        if (initDataPtr->hbd_mode_decision > EB_8_BIT_MD) {
            InitData data[] = {
                {
                    &object_ptr->md_luma_recon_neighbor_array16bit[depth],
                    MAX_PICTURE_WIDTH_SIZE,
                    MAX_PICTURE_HEIGHT_SIZE,
                    sizeof(uint16_t),
                    SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                    SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                    NEIGHBOR_ARRAY_UNIT_FULL_MASK,
                },
                {
                    &object_ptr->md_tx_depth_1_luma_recon_neighbor_array16bit[depth],
                    MAX_PICTURE_WIDTH_SIZE,
                    MAX_PICTURE_HEIGHT_SIZE,
                    sizeof(uint16_t),
                    SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                    SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                    NEIGHBOR_ARRAY_UNIT_FULL_MASK,
                },
                {
                    &object_ptr->md_cb_recon_neighbor_array16bit[depth],
                    MAX_PICTURE_WIDTH_SIZE >> subsampling_x,
                    MAX_PICTURE_HEIGHT_SIZE >> subsampling_y,
                    sizeof(uint16_t),
                    SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                    SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                    NEIGHBOR_ARRAY_UNIT_FULL_MASK,
                },
                {
                    &object_ptr->md_cr_recon_neighbor_array16bit[depth],
                    MAX_PICTURE_WIDTH_SIZE >> subsampling_x,
                    MAX_PICTURE_HEIGHT_SIZE >> subsampling_y,
                    sizeof(uint16_t),
                    SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                    SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                    NEIGHBOR_ARRAY_UNIT_FULL_MASK,
                }
            };
            return_error = create_neighbor_array_units(data, DIM(data));
            if (return_error == EB_ErrorInsufficientResources)
                return EB_ErrorInsufficientResources;
        }

        EB_NEW(
            object_ptr->md_interpolation_type_neighbor_array[depth],
            neighbor_array_unit_ctor32,
            MAX_PICTURE_WIDTH_SIZE,
            MAX_PICTURE_HEIGHT_SIZE,
            sizeof(uint32_t),
            PU_NEIGHBOR_ARRAY_GRANULARITY,
            PU_NEIGHBOR_ARRAY_GRANULARITY,
            NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
    }
    {
        InitData data[] = {
            {
                &object_ptr->ep_intra_luma_mode_neighbor_array,
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            // Encode Pass Neighbor Arrays
            {
                &object_ptr->ep_intra_chroma_mode_neighbor_array,
                MAX_PICTURE_WIDTH_SIZE >> subsampling_x,
                MAX_PICTURE_HEIGHT_SIZE >> subsampling_y,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            {
                &object_ptr->ep_mv_neighbor_array,
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(MvUnit),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_FULL_MASK,
            },
            {
                &object_ptr->ep_skip_flag_neighbor_array,
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                CU_NEIGHBOR_ARRAY_GRANULARITY,
                CU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            {
                &object_ptr->ep_mode_type_neighbor_array,
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_FULL_MASK,
            },
            {
                &object_ptr->ep_leaf_depth_neighbor_array,
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            {
                &object_ptr->ep_luma_recon_neighbor_array,
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_FULL_MASK,
            },
            {
                &object_ptr->ep_cb_recon_neighbor_array,
                MAX_PICTURE_WIDTH_SIZE >> subsampling_x,
                MAX_PICTURE_HEIGHT_SIZE >> subsampling_y,
                sizeof(uint8_t),
                SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_FULL_MASK,
            },
            {
                &object_ptr->ep_cr_recon_neighbor_array,
                MAX_PICTURE_WIDTH_SIZE >> subsampling_x,
                MAX_PICTURE_HEIGHT_SIZE >> subsampling_y,
                sizeof(uint8_t),
                SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_FULL_MASK,
            },
           // for each 4x4
            {
                &object_ptr->ep_luma_dc_sign_level_coeff_neighbor_array,
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            // for each 4x4
            {
                &object_ptr->ep_cb_dc_sign_level_coeff_neighbor_array,
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            // for each 4x4
            {
                &object_ptr->ep_cr_dc_sign_level_coeff_neighbor_array,
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            // Encode pass partition neighbor array
            {
                &object_ptr->ep_partition_context_neighbor_array,
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(struct PartitionContext),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },

            // Entropy Coding Neighbor Arrays
            {
                &object_ptr->mode_type_neighbor_array,
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            {
                &object_ptr->partition_context_neighbor_array,
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(struct PartitionContext),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            {
                &object_ptr->skip_flag_neighbor_array,
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            {
                &object_ptr->skip_coeff_neighbor_array,
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            // for each 4x4
            {
                &object_ptr->luma_dc_sign_level_coeff_neighbor_array,
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            // for each 4x4
            {
                &object_ptr->cr_dc_sign_level_coeff_neighbor_array,
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            // for each 4x4
            {
                &object_ptr->cb_dc_sign_level_coeff_neighbor_array,
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            {
                &object_ptr->inter_pred_dir_neighbor_array,
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            {
                &object_ptr->ref_frame_type_neighbor_array,
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            {
                &object_ptr->intra_luma_mode_neighbor_array,
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            {
                &object_ptr->txfm_context_array,
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(TXFM_CONTEXT),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK,
            },
            {
                &object_ptr->segmentation_id_pred_array,
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint8_t),
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                PU_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_FULL_MASK,
            },
        };
        return_error = create_neighbor_array_units(data, DIM(data));
        if (return_error == EB_ErrorInsufficientResources)
            return EB_ErrorInsufficientResources;
    }
    if (is16bit) {
        InitData data[] = {
            {
                &object_ptr->ep_luma_recon_neighbor_array16bit,
                MAX_PICTURE_WIDTH_SIZE,
                MAX_PICTURE_HEIGHT_SIZE,
                sizeof(uint16_t),
                SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_FULL_MASK,
            },
            {
                &object_ptr->ep_cb_recon_neighbor_array16bit,
                MAX_PICTURE_WIDTH_SIZE >> subsampling_x,
                MAX_PICTURE_HEIGHT_SIZE >> subsampling_y,
                sizeof(uint16_t),
                SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_FULL_MASK,
            },
            {
                &object_ptr->ep_cr_recon_neighbor_array16bit,
                MAX_PICTURE_WIDTH_SIZE >> subsampling_x,
                MAX_PICTURE_HEIGHT_SIZE >> subsampling_y,
                sizeof(uint16_t),
                SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
                NEIGHBOR_ARRAY_UNIT_FULL_MASK,
            },
        };
        return_error = create_neighbor_array_units(data, DIM(data));
        if (return_error == EB_ErrorInsufficientResources)
            return EB_ErrorInsufficientResources;
    }
    else {
        object_ptr->ep_luma_recon_neighbor_array16bit = 0;
        object_ptr->ep_cb_recon_neighbor_array16bit = 0;
        object_ptr->ep_cr_recon_neighbor_array16bit = 0;
    }
    EB_NEW(
        object_ptr->interpolation_type_neighbor_array,
        neighbor_array_unit_ctor32,
        MAX_PICTURE_WIDTH_SIZE,
        MAX_PICTURE_HEIGHT_SIZE,
        sizeof(uint32_t),
        PU_NEIGHBOR_ARRAY_GRANULARITY,
        PU_NEIGHBOR_ARRAY_GRANULARITY,
        NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

    //Segmentation neighbor arrays
    EB_NEW(
        object_ptr->segmentation_neighbor_map,
        segmentation_map_ctor,
        initDataPtr->picture_width, initDataPtr->picture_height);

    // Note - non-zero offsets are not supported (to be fixed later in DLF chroma filtering)
    object_ptr->cb_qp_offset = 0;
    object_ptr->cr_qp_offset = 0;

    object_ptr->slice_level_chroma_qp_flag = EB_TRUE;
    // slice level chroma QP offsets
    object_ptr->slice_cb_qp_offset = 0;
    object_ptr->slice_cr_qp_offset = 0;

    //object_ptr->total_num_bits = 0;

    // Error Resilience
    object_ptr->constrained_intra_flag = EB_FALSE;

    // Segments
    EB_NEW(
        object_ptr->enc_dec_segment_ctrl,
        enc_dec_segments_ctor,
        initDataPtr->enc_dec_segment_col,
        initDataPtr->enc_dec_segment_row);
    // Entropy Rows
    EB_CREATE_MUTEX(object_ptr->entropy_coding_mutex);

    EB_CREATE_MUTEX(object_ptr->intra_mutex);

    EB_CREATE_MUTEX(object_ptr->cdef_search_mutex);

    //object_ptr->mse_seg[0] = (uint64_t(*)[64])eb_aom_malloc(sizeof(**object_ptr->mse_seg) *  pictureLcuWidth * pictureLcuHeight);
   // object_ptr->mse_seg[1] = (uint64_t(*)[64])eb_aom_malloc(sizeof(**object_ptr->mse_seg) *  pictureLcuWidth * pictureLcuHeight);

    EB_MALLOC_ARRAY(object_ptr->mse_seg[0], pictureLcuWidth * pictureLcuHeight);
    EB_MALLOC_ARRAY(object_ptr->mse_seg[1], pictureLcuWidth * pictureLcuHeight);

    EB_CREATE_MUTEX(object_ptr->rest_search_mutex);

    //the granularity is 4x4
    EB_MALLOC_ARRAY(object_ptr->mi_grid_base, all_sb*(initDataPtr->sb_size_pix >> MI_SIZE_LOG2)*(initDataPtr->sb_size_pix >> MI_SIZE_LOG2));

    EB_MALLOC_ARRAY(object_ptr->mip, all_sb*(initDataPtr->sb_size_pix >> MI_SIZE_LOG2)*(initDataPtr->sb_size_pix >> MI_SIZE_LOG2));

    memset(object_ptr->mip, 0, sizeof(ModeInfo) * all_sb*(initDataPtr->sb_size_pix >> MI_SIZE_LOG2)*(initDataPtr->sb_size_pix >> MI_SIZE_LOG2));

    uint32_t miIdx;
    for (miIdx = 0; miIdx < all_sb*(initDataPtr->sb_size_pix >> MI_SIZE_LOG2)*(initDataPtr->sb_size_pix >> MI_SIZE_LOG2); ++miIdx)
        object_ptr->mi_grid_base[miIdx] = object_ptr->mip + miIdx;
    object_ptr->mi_stride = picture_sb_w * (initDataPtr->sb_size_pix >> MI_SIZE_LOG2);
    if (initDataPtr->mfmv)
    {
        //MFMV: map is 8x8 based.
        uint32_t mi_rows = initDataPtr->picture_height >> MI_SIZE_LOG2;
        const int mem_size = ((mi_rows + MAX_MIB_SIZE) >> 1) * (object_ptr->mi_stride >> 1);

        EB_CALLOC_ALIGNED_ARRAY(object_ptr->tpl_mvs, mem_size);
    }
    object_ptr->hash_table.p_lookup_table = NULL;
    av1_hash_table_create(&object_ptr->hash_table);
    return EB_ErrorNone;
}

EbErrorType picture_control_set_creator(
    EbPtr *object_dbl_ptr,
    EbPtr object_init_data_ptr)
{
    PictureControlSet* obj;

    *object_dbl_ptr = NULL;
    EB_NEW(obj, picture_control_set_ctor, object_init_data_ptr);
    *object_dbl_ptr = obj;

    return EB_ErrorNone;
}

static void picture_parent_control_set_dctor(EbPtr p)
{
    PictureParentControlSet *obj = (PictureParentControlSet*)p;
    uint32_t regionInPictureWidthIndex;
    uint32_t regionInPictureHeightIndex;

    EB_DELETE(obj->denoise_and_model);

    EB_DELETE_PTR_ARRAY(obj->me_results, obj->sb_total_count);
    if (obj->is_chroma_downsampled_picture_ptr_owner)
        EB_DELETE(obj->chroma_downsampled_picture_ptr);

    EB_FREE_2D(obj->variance);
    EB_FREE_2D(obj->y_mean);
    EB_FREE_2D(obj->cbMean);
    EB_FREE_2D(obj->crMean);

    if (obj->picture_histogram) {
        for (regionInPictureWidthIndex = 0; regionInPictureWidthIndex < MAX_NUMBER_OF_REGIONS_IN_WIDTH; regionInPictureWidthIndex++) {
            if (obj->picture_histogram[regionInPictureWidthIndex]) {
                for (regionInPictureHeightIndex = 0; regionInPictureHeightIndex < MAX_NUMBER_OF_REGIONS_IN_HEIGHT; regionInPictureHeightIndex++) {
                    EB_FREE_2D(obj->picture_histogram[regionInPictureWidthIndex][regionInPictureHeightIndex]);
                }
            }
            EB_FREE_PTR_ARRAY(obj->picture_histogram[regionInPictureWidthIndex], MAX_NUMBER_OF_REGIONS_IN_HEIGHT);
        }
        EB_FREE_PTR_ARRAY(obj->picture_histogram, MAX_NUMBER_OF_REGIONS_IN_WIDTH);
    }

    EB_FREE_2D(obj->ois_sb_results);
    EB_FREE_2D(obj->ois_candicate);
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
    EB_FREE_ARRAY(obj->edge_results_ptr);
    EB_FREE_ARRAY(obj->sharp_edge_sb_flag);
    EB_FREE_ARRAY(obj->sb_stat_array);
    EB_FREE_ARRAY(obj->sb_depth_mode_array);

    if (obj->av1_cm) {
        const int32_t num_planes = 3;// av1_num_planes(cm);
        for (int32_t p = 0; p < num_planes; ++p) {
            RestorationInfo* ri =obj->av1_cm->rst_info + p;
            RestorationStripeBoundaries *boundaries = &ri->boundaries;
            EB_FREE_ARRAY(ri->unit_info);
            EB_FREE(boundaries->stripe_boundary_above);
            EB_FREE(boundaries->stripe_boundary_below);
        }
        EB_FREE_ARRAY(obj->av1_cm->frame_to_show);
        EB_FREE_ALIGNED(obj->av1_cm->rst_tmpbuf);
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
}
EbErrorType picture_parent_control_set_ctor(
    PictureParentControlSet *object_ptr,
    EbPtr object_init_data_ptr)
{
    PictureControlSetInitData *initDataPtr = (PictureControlSetInitData*)object_init_data_ptr;
    EbErrorType return_error = EB_ErrorNone;
    const uint16_t pictureLcuWidth = (uint16_t)((initDataPtr->picture_width + initDataPtr->sb_sz - 1) / initDataPtr->sb_sz);
    const uint16_t pictureLcuHeight = (uint16_t)((initDataPtr->picture_height + initDataPtr->sb_sz - 1) / initDataPtr->sb_sz);
    const uint16_t subsampling_x = (initDataPtr->color_format == EB_YUV444 ? 1 : 2) - 1;
    const uint16_t subsampling_y = (initDataPtr->color_format >= EB_YUV422 ? 1 : 2) - 1;
    uint16_t sb_index;
    uint32_t regionInPictureWidthIndex;
    uint32_t regionInPictureHeightIndex;

    object_ptr->dctor = picture_parent_control_set_dctor;

    object_ptr->sequence_control_set_wrapper_ptr = (EbObjectWrapper *)EB_NULL;
    object_ptr->input_picture_wrapper_ptr = (EbObjectWrapper *)EB_NULL;
    object_ptr->reference_picture_wrapper_ptr = (EbObjectWrapper *)EB_NULL;
    object_ptr->enhanced_picture_ptr = (EbPictureBufferDesc *)EB_NULL;

    if (initDataPtr->color_format >= EB_YUV422) {
        EbPictureBufferDescInitData input_picture_buffer_desc_init_data;
        input_picture_buffer_desc_init_data.max_width     = initDataPtr->picture_width;
        input_picture_buffer_desc_init_data.max_height    = initDataPtr->picture_height;
        input_picture_buffer_desc_init_data.bit_depth = 8; //Should be 8bit
        input_picture_buffer_desc_init_data.buffer_enable_mask = PICTURE_BUFFER_DESC_CHROMA_MASK;
        input_picture_buffer_desc_init_data.left_padding = initDataPtr->left_padding;
        input_picture_buffer_desc_init_data.right_padding = initDataPtr->right_padding;
        input_picture_buffer_desc_init_data.top_padding = initDataPtr->top_padding;
        input_picture_buffer_desc_init_data.bot_padding = initDataPtr->bot_padding;
        input_picture_buffer_desc_init_data.color_format = EB_YUV420; //set to 420 for MD
        input_picture_buffer_desc_init_data.split_mode    = EB_FALSE;
        EB_NEW(
            object_ptr->chroma_downsampled_picture_ptr,
            eb_picture_buffer_desc_ctor,
            (EbPtr) &input_picture_buffer_desc_init_data);
        object_ptr->is_chroma_downsampled_picture_ptr_owner = EB_TRUE;
    } else if(initDataPtr->color_format == EB_YUV420) {
        object_ptr->chroma_downsampled_picture_ptr = NULL;
    } else
        return EB_ErrorBadParameter;
    // GOP
    object_ptr->pred_struct_index = 0;
    object_ptr->picture_number = 0;
    object_ptr->idr_flag = EB_FALSE;
    object_ptr->temporal_layer_index = 0;
    object_ptr->total_num_bits = 0;
    object_ptr->last_idr_picture = 0;
    object_ptr->sb_total_count = pictureLcuWidth * pictureLcuHeight;

    object_ptr->data_ll_head_ptr = (EbLinkedListNode *)EB_NULL;
    object_ptr->app_out_data_ll_head_ptr = (EbLinkedListNode *)EB_NULL;
    EB_MALLOC_2D(object_ptr->variance, object_ptr->sb_total_count, MAX_ME_PU_COUNT);
    EB_MALLOC_2D(object_ptr->y_mean, object_ptr->sb_total_count, MAX_ME_PU_COUNT);
    EB_MALLOC_2D(object_ptr->cbMean, object_ptr->sb_total_count, 21);
    EB_MALLOC_2D(object_ptr->crMean, object_ptr->sb_total_count, 21);


    EB_ALLOC_PTR_ARRAY(object_ptr->picture_histogram, MAX_NUMBER_OF_REGIONS_IN_WIDTH);

    for (regionInPictureWidthIndex = 0; regionInPictureWidthIndex < MAX_NUMBER_OF_REGIONS_IN_WIDTH; regionInPictureWidthIndex++) {  // loop over horizontal regions
        EB_ALLOC_PTR_ARRAY(object_ptr->picture_histogram[regionInPictureWidthIndex], MAX_NUMBER_OF_REGIONS_IN_HEIGHT);
        for (regionInPictureHeightIndex = 0; regionInPictureHeightIndex < MAX_NUMBER_OF_REGIONS_IN_HEIGHT; regionInPictureHeightIndex++) {
            EB_MALLOC_2D(object_ptr->picture_histogram[regionInPictureWidthIndex][regionInPictureHeightIndex], 3, HISTOGRAM_NUMBER_OF_BINS);
        }
    }

    EB_MALLOC_2D(object_ptr->ois_sb_results, object_ptr->sb_total_count, 1);
    EB_MALLOC_2D(object_ptr->ois_candicate, object_ptr->sb_total_count,  MAX_OIS_CANDIDATES * CU_MAX_COUNT);

    for (sb_index = 0; sb_index < object_ptr->sb_total_count; ++sb_index) {
        uint32_t cuIdx;
        for (cuIdx = 0; cuIdx < CU_MAX_COUNT; ++cuIdx)
            object_ptr->ois_sb_results[sb_index]->ois_candidate_array[cuIdx] = &object_ptr->ois_candicate[sb_index][cuIdx*MAX_OIS_CANDIDATES];
    }

    object_ptr->max_number_of_candidates_per_block = (initDataPtr->mrp_mode == 0) ?
        ME_RES_CAND_MRP_MODE_0 : // [Single Ref = 7] + [BiDir = 12 = 3*4 ] + [UniDir = 4 = 3+1]
        ME_RES_CAND_MRP_MODE_1 ; // [BiDir = 1] + [UniDir = 2 = 1 + 1]

    EB_ALLOC_PTR_ARRAY(object_ptr->me_results, object_ptr->sb_total_count);

    for (sb_index = 0; sb_index < object_ptr->sb_total_count; ++sb_index) {
        EB_NEW(
            object_ptr->me_results[sb_index],
            me_sb_results_ctor,
            (initDataPtr->nsq_present) ? MAX_ME_PU_COUNT : SQUARE_PU_COUNT,
            initDataPtr->mrp_mode,
            object_ptr->max_number_of_candidates_per_block);
    }

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
    EB_MALLOC_ARRAY(object_ptr->edge_results_ptr, object_ptr->sb_total_count);

    EB_MALLOC_ARRAY(object_ptr->sharp_edge_sb_flag, object_ptr->sb_total_count);
    EB_MALLOC_ARRAY(object_ptr->sb_stat_array, object_ptr->sb_total_count);
    EB_CREATE_MUTEX(object_ptr->rc_distortion_histogram_mutex);
    EB_MALLOC_ARRAY(object_ptr->sb_depth_mode_array, object_ptr->sb_total_count);
    EB_CREATE_SEMAPHORE(object_ptr->temp_filt_done_semaphore, 0, 1);
    EB_CREATE_MUTEX(object_ptr->temp_filt_mutex);
    EB_CREATE_MUTEX(object_ptr->debug_mutex);
    EB_MALLOC_ARRAY(object_ptr->av1_cm, 1);

    object_ptr->av1_cm->interp_filter = SWITCHABLE;

    object_ptr->av1_cm->mi_stride = pictureLcuWidth * (BLOCK_SIZE_64 / 4);

    object_ptr->av1_cm->p_pcs_ptr = object_ptr;

    EB_MALLOC_ARRAY(object_ptr->av1_cm->frame_to_show, 1);

    object_ptr->av1_cm->use_highbitdepth = (initDataPtr->bit_depth > 8 ? 1 : 0);
    object_ptr->av1_cm->bit_depth = initDataPtr->bit_depth;
    object_ptr->av1_cm->color_format = initDataPtr->color_format;
    object_ptr->av1_cm->subsampling_x = subsampling_x;
    object_ptr->av1_cm->subsampling_y = subsampling_y;
    object_ptr->av1_cm->frm_size.frame_width = initDataPtr->picture_width;
    object_ptr->av1_cm->frm_size.frame_height = initDataPtr->picture_height;
    object_ptr->av1_cm->frm_size.superres_upscaled_width = initDataPtr->picture_width;
    object_ptr->av1_cm->frm_size.superres_upscaled_height = initDataPtr->picture_height;

    object_ptr->av1_cm->mi_cols = initDataPtr->picture_width >> MI_SIZE_LOG2;
    object_ptr->av1_cm->mi_rows = initDataPtr->picture_height >> MI_SIZE_LOG2;

    object_ptr->av1_cm->byte_alignment = 0;

    set_restoration_unit_size(initDataPtr->picture_width, initDataPtr->picture_height, 1, 1, object_ptr->av1_cm->rst_info);

    return_error = eb_av1_alloc_restoration_buffers(object_ptr->av1_cm);

    memset(&object_ptr->av1_cm->rst_frame, 0, sizeof(Yv12BufferConfig));

    int32_t ntiles[2];
    for (int32_t is_uv = 0; is_uv < 2; ++is_uv)
        ntiles[is_uv] = object_ptr->av1_cm->rst_info[is_uv].units_per_tile; //CHKN res_tiles_in_plane

    assert(ntiles[1] <= ntiles[0]);

    EB_CALLOC_ARRAY(object_ptr->rusi_picture[0], ntiles[0]);
    EB_CALLOC_ARRAY(object_ptr->rusi_picture[1], ntiles[1]);
    EB_CALLOC_ARRAY(object_ptr->rusi_picture[2], ntiles[1]);

    EB_MALLOC_ARRAY(object_ptr->av1x, 1);

    // Film grain noise model if film grain is applied
    if (initDataPtr->film_grain_noise_level) {
        denoise_and_model_init_data_t fg_init_data;
        fg_init_data.encoder_bit_depth = initDataPtr->bit_depth;
        fg_init_data.encoder_color_format = initDataPtr->color_format;
        fg_init_data.noise_level = initDataPtr->film_grain_noise_level;
        fg_init_data.width = initDataPtr->picture_width;
        fg_init_data.height = initDataPtr->picture_height;
        fg_init_data.stride_y = initDataPtr->picture_width + initDataPtr->left_padding + initDataPtr->right_padding;
        fg_init_data.stride_cb = fg_init_data.stride_cr = fg_init_data.stride_y >> subsampling_x;

        EB_NEW(object_ptr->denoise_and_model, denoise_and_model_ctor,
            (EbPtr)&fg_init_data);
    }

    return return_error;
}

EbErrorType picture_parent_control_set_creator(
    EbPtr *object_dbl_ptr,
    EbPtr object_init_data_ptr)
{
    PictureParentControlSet* obj;

    *object_dbl_ptr = NULL;
    EB_NEW(obj, picture_parent_control_set_ctor, object_init_data_ptr);
    *object_dbl_ptr = obj;

    return EB_ErrorNone;
}
// clang-format on
