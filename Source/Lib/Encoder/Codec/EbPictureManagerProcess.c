/*
* Copyright(c) 2019 Intel Corporation
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

#include <stdlib.h>
#include <string.h>

#include "EbEncHandle.h"
#include "EbUtility.h"
#include "EbPictureControlSet.h"
#include "EbSequenceControlSet.h"
#include "EbPictureManagerProcess.h"
#include "EbReferenceObject.h"
#include "EbPictureDemuxResults.h"
#include "EbPictureManagerQueue.h"
#include "EbPredictionStructure.h"
#include "EbRateControlTasks.h"
#include "EbSvtAv1ErrorCodes.h"
#include "EbEntropyCoding.h"
#include "EbLog.h"

// Token buffer is only used for palette tokens.
static INLINE unsigned int get_token_alloc(int mb_rows, int mb_cols, int sb_size_log2, const int num_planes) {
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
void                    rtime_alloc_palette_tokens(SequenceControlSet *scs, PictureControlSet *child_pcs_ptr);
extern MvReferenceFrame svt_get_ref_frame_type(uint8_t list, uint8_t ref_idx);

void svt_aom_largest_coding_unit_dctor(EbPtr p);
void write_stat_to_file(SequenceControlSet *scs, StatStruct stat_struct, uint64_t ref_poc);

static void picture_manager_context_dctor(EbPtr p) {
    EbThreadContext       *thread_ctx = (EbThreadContext *)p;
    PictureManagerContext *obj        = (PictureManagerContext *)thread_ctx->priv;
    EB_FREE_ARRAY(obj);
}
/************************************************
 * Picture Manager Context Constructor
 ************************************************/
EbErrorType svt_aom_picture_manager_context_ctor(EbThreadContext *thread_ctx, const EbEncHandle *enc_handle_ptr,
                                                 int rate_control_index) {
    PictureManagerContext *context_ptr;
    EB_CALLOC_ARRAY(context_ptr, 1);
    thread_ctx->priv  = context_ptr;
    thread_ctx->dctor = picture_manager_context_dctor;

    context_ptr->picture_input_fifo_ptr = svt_system_resource_get_consumer_fifo(
        enc_handle_ptr->picture_demux_results_resource_ptr, 0);
    context_ptr->picture_manager_output_fifo_ptr = svt_system_resource_get_producer_fifo(
        enc_handle_ptr->rate_control_tasks_resource_ptr, rate_control_index);
    context_ptr->picture_control_set_fifo_ptr = svt_system_resource_get_producer_fifo(
        enc_handle_ptr->picture_control_set_pool_ptr_array[0], 0); //The Child PCS Pool here
    context_ptr->recon_coef_fifo_ptr = svt_system_resource_get_producer_fifo(enc_handle_ptr->enc_dec_pool_ptr_array[0],
                                                                             0); //The Child PCS Pool here

    context_ptr->consecutive_dec_order           = 0;
    context_ptr->started_pics_dec_order_head_idx = 0;
    context_ptr->started_pics_dec_order_tail_idx = 0;

    return EB_ErrorNone;
}

void svt_aom_copy_buffer_info(EbPictureBufferDesc *src_ptr, EbPictureBufferDesc *dst_ptr) {
    dst_ptr->width             = src_ptr->width;
    dst_ptr->height            = src_ptr->height;
    dst_ptr->max_width         = src_ptr->max_width;
    dst_ptr->max_height        = src_ptr->max_height;
    dst_ptr->stride_y          = src_ptr->stride_y;
    dst_ptr->stride_cb         = src_ptr->stride_cb;
    dst_ptr->stride_cr         = src_ptr->stride_cr;
    dst_ptr->org_x             = src_ptr->org_x;
    dst_ptr->origin_bot_y      = src_ptr->origin_bot_y;
    dst_ptr->org_y             = src_ptr->org_y;
    dst_ptr->stride_bit_inc_y  = src_ptr->stride_bit_inc_y;
    dst_ptr->stride_bit_inc_cb = src_ptr->stride_bit_inc_cb;
    dst_ptr->stride_bit_inc_cr = src_ptr->stride_bit_inc_cr;
    dst_ptr->luma_size         = src_ptr->luma_size;
    dst_ptr->chroma_size       = src_ptr->chroma_size;
}

void svt_aom_set_tile_info(PictureParentControlSet *pcs);
/*
  walks the ref picture list, and looks for a particular pic
  return NULL if not found.
*/
ReferenceQueueEntry *search_ref_in_ref_queue(EncodeContext *enc_ctx, uint64_t ref_poc) {
    ReferenceQueueEntry *ref_pic_entry = NULL;
    for (uint32_t i = 0; i < enc_ctx->ref_pic_list_length; i++) {
        ref_pic_entry = enc_ctx->ref_pic_list[i];
        if (ref_pic_entry->is_valid && ref_pic_entry->picture_number == ref_poc) {
            return ref_pic_entry;
        }
    }

    return NULL;
}

void svt_aom_init_enc_dec_segement(PictureParentControlSet *ppcs) {
    SequenceControlSet *scs                  = ppcs->scs;
    uint8_t             pic_width_in_sb      = (uint8_t)((ppcs->aligned_width + scs->sb_size - 1) / scs->sb_size);
    uint8_t             picture_height_in_sb = (uint8_t)((ppcs->aligned_height + scs->sb_size - 1) / scs->sb_size);
    svt_aom_set_tile_info(ppcs);
    int              sb_size_log2        = scs->seq_header.sb_size_log2;
    uint32_t         enc_dec_seg_col_cnt = scs->enc_dec_segment_col_count_array[ppcs->temporal_layer_index];
    uint32_t         enc_dec_seg_row_cnt = scs->enc_dec_segment_row_count_array[ppcs->temporal_layer_index];
    Av1Common *const cm                  = ppcs->av1_cm;
    const int        tile_cols           = ppcs->av1_cm->tiles_info.tile_cols;
    const int        tile_rows           = ppcs->av1_cm->tiles_info.tile_rows;
    uint8_t          tile_group_cols     = MIN(tile_cols, scs->tile_group_col_count_array[ppcs->temporal_layer_index]);
    uint8_t          tile_group_rows     = MIN(tile_rows, scs->tile_group_row_count_array[ppcs->temporal_layer_index]);

    if (tile_group_cols * tile_group_rows > 1) {
        enc_dec_seg_col_cnt = MIN(enc_dec_seg_col_cnt, (uint8_t)(pic_width_in_sb / tile_group_cols));
        enc_dec_seg_row_cnt = MIN(enc_dec_seg_row_cnt, (uint8_t)(picture_height_in_sb / tile_group_rows));
    }
    ppcs->tile_group_cols = tile_group_cols;
    ppcs->tile_group_rows = tile_group_rows;

    uint8_t tile_group_col_start_tile_idx[1024];
    uint8_t tile_group_row_start_tile_idx[1024];

    // Get the tile start index for tile group
    for (uint8_t c = 0; c <= tile_group_cols; c++) {
        tile_group_col_start_tile_idx[c] = c * tile_cols / tile_group_cols;
    }
    for (uint8_t r = 0; r <= tile_group_rows; r++) {
        tile_group_row_start_tile_idx[r] = r * tile_rows / tile_group_rows;
    }
    for (uint8_t r = 0; r < tile_group_rows; r++) {
        for (uint8_t c = 0; c < tile_group_cols; c++) {
            uint16_t tile_group_idx            = r * tile_group_cols + c;
            uint16_t top_left_tile_col_idx     = tile_group_col_start_tile_idx[c];
            uint16_t top_left_tile_row_idx     = tile_group_row_start_tile_idx[r];
            uint16_t bottom_right_tile_col_idx = tile_group_col_start_tile_idx[c + 1];
            uint16_t bottom_right_tile_row_idx = tile_group_row_start_tile_idx[r + 1];

            TileGroupInfo *tg_info_ptr = &ppcs->tile_group_info[tile_group_idx];

            tg_info_ptr->tile_group_tile_start_x = top_left_tile_col_idx;
            tg_info_ptr->tile_group_tile_end_x   = bottom_right_tile_col_idx;

            tg_info_ptr->tile_group_tile_start_y = top_left_tile_row_idx;
            tg_info_ptr->tile_group_tile_end_y   = bottom_right_tile_row_idx;

            tg_info_ptr->tile_group_sb_start_x = cm->tiles_info.tile_col_start_mi[top_left_tile_col_idx] >>
                sb_size_log2;
            tg_info_ptr->tile_group_sb_start_y = cm->tiles_info.tile_row_start_mi[top_left_tile_row_idx] >>
                sb_size_log2;

            // Get the SB end of the bottom right tile
            tg_info_ptr->tile_group_sb_end_x = (cm->tiles_info.tile_col_start_mi[bottom_right_tile_col_idx] >>
                                                sb_size_log2);
            tg_info_ptr->tile_group_sb_end_y = (cm->tiles_info.tile_row_start_mi[bottom_right_tile_row_idx] >>
                                                sb_size_log2);

            // Get the width/height of tile group in SB
            tg_info_ptr->tile_group_height_in_sb = tg_info_ptr->tile_group_sb_end_y -
                tg_info_ptr->tile_group_sb_start_y;
            tg_info_ptr->tile_group_width_in_sb = tg_info_ptr->tile_group_sb_end_x - tg_info_ptr->tile_group_sb_start_x;

            // Init segments within the tile group
            svt_aom_enc_dec_segments_init(ppcs->child_pcs->enc_dec_segment_ctrl[tile_group_idx],
                                          enc_dec_seg_col_cnt,
                                          enc_dec_seg_row_cnt,
                                          tg_info_ptr->tile_group_width_in_sb,
                                          tg_info_ptr->tile_group_height_in_sb);
            // Enable tile parallelism in Entropy Coding stage
            for (uint16_t s = top_left_tile_row_idx; s < bottom_right_tile_row_idx; s++) {
                for (uint16_t d = top_left_tile_col_idx; d < bottom_right_tile_col_idx; d++) {
                    uint16_t tileIdx                                            = s * tile_cols + d;
                    ppcs->child_pcs->ec_info[tileIdx]->entropy_coding_tile_done = FALSE;
                }
            }
            ppcs->child_pcs->entropy_coding_pic_reset_flag = TRUE;
        }
    }
}

void superres_setup_child_pcs(SequenceControlSet *entry_scs_ptr, PictureParentControlSet *entry_ppcs) {
    PictureControlSet *child_pcs = entry_ppcs->child_pcs;

    uint8_t pic_width_in_sb;
    uint8_t picture_height_in_sb;

    child_pcs->b64_total_count = entry_ppcs->b64_total_count;

    pic_width_in_sb      = (uint8_t)((entry_ppcs->aligned_width + entry_scs_ptr->sb_size - 1) / entry_scs_ptr->sb_size);
    picture_height_in_sb = (uint8_t)((entry_ppcs->aligned_height + entry_scs_ptr->sb_size - 1) /
                                     entry_scs_ptr->sb_size);

    child_pcs->sb_total_count = pic_width_in_sb * picture_height_in_sb;

    // if (entry_ppcs->frame_superres_enabled)
    {
        // Modify sb_prt_array in child pcs
        uint16_t sb_index;
        uint16_t sb_origin_x = 0;
        uint16_t sb_origin_y = 0;
        for (sb_index = 0; sb_index < child_pcs->sb_total_count; ++sb_index) {
            svt_aom_largest_coding_unit_dctor(child_pcs->sb_ptr_array[sb_index]);
            svt_aom_largest_coding_unit_ctor(child_pcs->sb_ptr_array[sb_index],
                                             (uint8_t)entry_scs_ptr->sb_size,
                                             (uint16_t)(sb_origin_x * entry_scs_ptr->sb_size),
                                             (uint16_t)(sb_origin_y * entry_scs_ptr->sb_size),
                                             (uint16_t)sb_index,
                                             child_pcs->enc_mode,
                                             entry_scs_ptr->max_block_cnt,
                                             child_pcs);
            // Increment the Order in coding order (Raster Scan Order)
            sb_origin_y = (sb_origin_x == pic_width_in_sb - 1) ? sb_origin_y + 1 : sb_origin_y;
            sb_origin_x = (sb_origin_x == pic_width_in_sb - 1) ? 0 : sb_origin_x + 1;
        }
    }

    // Update pcs->mi_stride
    child_pcs->mi_stride = pic_width_in_sb * (entry_scs_ptr->sb_size >> MI_SIZE_LOG2);

    // init segment since picture scaled
    svt_aom_init_enc_dec_segement(entry_ppcs);

    // Tile Loop
    int              sb_size_log2 = entry_scs_ptr->seq_header.sb_size_log2;
    Av1Common *const cm           = entry_ppcs->av1_cm;
    const int        tile_cols    = entry_ppcs->av1_cm->tiles_info.tile_cols;
    const int        tile_rows    = entry_ppcs->av1_cm->tiles_info.tile_rows;
    uint32_t         x_sb_index, y_sb_index;
    uint16_t         tile_row, tile_col;
    TileInfo         tile_info;
    for (tile_row = 0; tile_row < tile_rows; tile_row++) {
        svt_av1_tile_set_row(&tile_info, &cm->tiles_info, cm->mi_rows, tile_row);

        for (tile_col = 0; tile_col < tile_cols; tile_col++) {
            svt_av1_tile_set_col(&tile_info, &cm->tiles_info, cm->mi_cols, tile_col);
            tile_info.tile_rs_index = tile_col + tile_row * tile_cols;

            for ((y_sb_index = cm->tiles_info.tile_row_start_mi[tile_row] >> sb_size_log2);
                 (y_sb_index < ((uint32_t)cm->tiles_info.tile_row_start_mi[tile_row + 1] >> sb_size_log2));
                 y_sb_index++) {
                for ((x_sb_index = cm->tiles_info.tile_col_start_mi[tile_col] >> sb_size_log2);
                     (x_sb_index < ((uint32_t)cm->tiles_info.tile_col_start_mi[tile_col + 1] >> sb_size_log2));
                     x_sb_index++) {
                    int sb_index = (uint16_t)(x_sb_index + y_sb_index * pic_width_in_sb);
                    child_pcs->sb_ptr_array[sb_index]->tile_info = tile_info;
                }
            }
        }
    }

    child_pcs->enc_dec_coded_sb_count = 0;
}

/* Picture Manager Kernel */

/***************************************************************************************************
 *
 * @brief
 *  The Picture Manager Process performs the function of managing both the Input Picture buffers and
 *  the Reference Picture buffers and subdividing the Input Picture into Tiles.
 *
 * @par Description:
 *  Both the Input Picture and Reference Picture buffers particular management depends on
 *  the GoP structure already implemented in the Picture decision. Also note that the
 *  Picture Manager sets up the RPS for Entropy Coding as well.
 *
 * @param[in] Input Picture
 *  Input Picture Data
 *
 * @param[in] Reference Picture
 *  Reference Picture Data
 *
 * @param[out] Picture Control Set
 *  Picture Control Set with fully available Reference List
 *
 ***************************************************************************************************/
void *svt_aom_picture_manager_kernel(void *input_ptr) {
    EbThreadContext       *thread_ctx  = (EbThreadContext *)input_ptr;
    PictureManagerContext *context_ptr = (PictureManagerContext *)thread_ctx->priv;

    EbObjectWrapper *enc_dec_wrapper;
    EncDecSet       *enc_dec_ptr;

    EbObjectWrapper         *child_pcs_wrapper;
    PictureControlSet       *child_pcs;
    PictureParentControlSet *pcs;
    SequenceControlSet      *scs;
    EncodeContext           *enc_ctx;

    EbObjectWrapper     *input_pic_demux_wrapper;
    PictureDemuxResults *input_pic_demux;

    Bool availability_flag;

    InputQueueEntry         *input_entry;
    uint32_t                 input_queue_index;
    ReferenceQueueEntry     *ref_entry = NULL;
    PictureParentControlSet *entry_ppcs;
    SequenceControlSet      *entry_scs_ptr;

    // Initialization
    uint16_t pic_width_in_sb;
    uint16_t picture_height_in_sb;
    uint64_t decode_order = 0;

    for (;;) {
        // Get Input Full Object
        EB_GET_FULL_OBJECT(context_ptr->picture_input_fifo_ptr, &input_pic_demux_wrapper);

        input_pic_demux = (PictureDemuxResults *)input_pic_demux_wrapper->object_ptr;

        // *Note - This should be overhauled and/or replaced when we
        //   need hierarchical support.

        switch (input_pic_demux->picture_type) {
        case EB_PIC_SUPERRES_INPUT: {
            pcs = (PictureParentControlSet *)input_pic_demux->pcs_wrapper->object_ptr;
            scs = pcs->scs;

            assert(scs->static_config.superres_mode == SUPERRES_QTHRESH ||
                   scs->static_config.superres_mode == SUPERRES_AUTO);

            // setup child pcs to reflect superres config. E.g. sb count, sb orig, tile info, etc.
            superres_setup_child_pcs(scs, pcs);

            EbObjectWrapper *out_results_wrapper;
            // Get Empty Results Object
            svt_get_empty_object(context_ptr->picture_manager_output_fifo_ptr, &out_results_wrapper);
            RateControlTasks *rc_tasks = (RateControlTasks *)out_results_wrapper->object_ptr;
            rc_tasks->pcs_wrapper      = pcs->child_pcs->c_pcs_wrapper_ptr;
            rc_tasks->task_type        = RC_INPUT_SUPERRES_RECODE;
            // Post the Full Results Object
            svt_post_full_object(out_results_wrapper);

            pcs     = (PictureParentControlSet *)NULL;
            enc_ctx = (EncodeContext *)NULL;
            break;
        }
        case EB_PIC_INPUT:

            pcs     = (PictureParentControlSet *)input_pic_demux->pcs_wrapper->object_ptr;
            scs     = pcs->scs;
            enc_ctx = scs->enc_ctx;

            // SVT_LOG("\nPicture Manager Process @ %d \n ", pcs->picture_number);
            CHECK_REPORT_ERROR(
                (((enc_ctx->input_picture_queue_head_index != enc_ctx->input_picture_queue_tail_index) ||
                  (enc_ctx->input_picture_queue[enc_ctx->input_picture_queue_head_index]->input_object_ptr == NULL))),
                enc_ctx->app_callback_ptr,
                EB_ENC_PM_ERROR4);

            // Place Picture in input queue
            input_entry                   = enc_ctx->input_picture_queue[enc_ctx->input_picture_queue_tail_index];
            input_entry->input_object_ptr = input_pic_demux->pcs_wrapper;
            enc_ctx->input_picture_queue_tail_index = (enc_ctx->input_picture_queue_tail_index ==
                                                       INPUT_QUEUE_MAX_DEPTH - 1)
                ? 0
                : enc_ctx->input_picture_queue_tail_index + 1;

            // Overlay pics should be NREF
            if (pcs->is_ref) {
#if OPT_LD_LATENCY2
                svt_block_on_mutex(enc_ctx->ref_pic_list_mutex);
#endif
                for (uint32_t i = 0; i < enc_ctx->ref_pic_list_length; i++) {
                    ref_entry = enc_ctx->ref_pic_list[i];
                    if (!ref_entry->is_valid) {
                        ref_entry->is_valid = 1;
                        break;
                    }
                    // Check if the reference list is full
                    CHECK_REPORT_ERROR(
                        i != enc_ctx->ref_pic_list_length - 1, enc_ctx->app_callback_ptr, EB_ENC_PM_ERROR5);
                }
                ref_entry->picture_number                = pcs->picture_number;
                ref_entry->reference_object_ptr          = (EbObjectWrapper *)NULL;
                ref_entry->release_enable                = TRUE;
                ref_entry->reference_available           = FALSE;
                ref_entry->slice_type                    = pcs->slice_type;
                ref_entry->temporal_layer_index          = pcs->temporal_layer_index;
                ref_entry->frame_context_updated         = FALSE;
                ref_entry->is_alt_ref                    = pcs->is_alt_ref;
                ref_entry->feedback_arrived              = FALSE;
                ref_entry->is_ref                        = pcs->is_ref;
                ref_entry->decode_order                  = pcs->decode_order;
                ref_entry->refresh_frame_mask            = pcs->av1_ref_signal.refresh_frame_mask;
                ref_entry->dec_order_of_last_ref         = pcs->is_ref ? UINT64_MAX : 0;
                ref_entry->frame_end_cdf_update_required = pcs->frame_end_cdf_update_mode;

                CHECK_REPORT_ERROR(
                    (pcs->pred_struct_ptr->pred_struct_period * REF_LIST_MAX_DEPTH < MAX_ELAPSED_IDR_COUNT),
                    enc_ctx->app_callback_ptr,
                    EB_ENC_PM_ERROR6);
#if OPT_LD_LATENCY2
                svt_release_mutex(enc_ctx->ref_pic_list_mutex);
#endif
            }
            break;

        case EB_PIC_REFERENCE:

            scs     = input_pic_demux->scs;
            enc_ctx = scs->enc_ctx;
            ((EbReferenceObject *)input_pic_demux->ref_pic_wrapper->object_ptr)->ds_pics.picture_number =
                input_pic_demux->picture_number;
            // Find the Reference in the Reference List
#if OPT_LD_LATENCY2
            svt_block_on_mutex(enc_ctx->ref_pic_list_mutex);
#endif
            for (uint32_t i = 0; i < enc_ctx->ref_pic_list_length; i++) {
                ref_entry = enc_ctx->ref_pic_list[i];
                if (ref_entry->is_valid && ref_entry->picture_number == input_pic_demux->picture_number) {
                    // Assign the reference object if there is a match
                    ref_entry->reference_object_ptr = input_pic_demux->ref_pic_wrapper;

                    // Set the reference availability
                    ref_entry->reference_available = TRUE;
                    break;
                }
                // Check if the reference list is full
                CHECK_REPORT_ERROR(i != enc_ctx->ref_pic_list_length - 1, enc_ctx->app_callback_ptr, EB_ENC_PM_ERROR5);
            }
            CHECK_REPORT_ERROR((ref_entry->picture_number == input_pic_demux->picture_number),
                               enc_ctx->app_callback_ptr,
                               EB_ENC_PM_ERROR8);
#if OPT_LD_LATENCY2
            svt_release_mutex(enc_ctx->ref_pic_list_mutex);
#endif
            break;
        case EB_PIC_FEEDBACK:
            scs     = input_pic_demux->scs;
            enc_ctx = scs->enc_ctx;

            // Find the Reference in the Reference Queue
#if OPT_LD_LATENCY2
            svt_block_on_mutex(enc_ctx->ref_pic_list_mutex);
#endif
            for (uint32_t i = 0; i < enc_ctx->ref_pic_list_length; i++) {
                ref_entry = enc_ctx->ref_pic_list[i];
                if (ref_entry->is_valid && ref_entry->picture_number == input_pic_demux->picture_number) {
                    // Set the feedback arrived
                    ref_entry->feedback_arrived      = TRUE;
                    ref_entry->frame_context_updated = TRUE;
                    break;
                }
                // Sometimes the reference may not be in the queue (e.g. if a delayed-I causes a ref
                // pic to not be needed, it may be dropped from the queue before feedback arrives).
            }
#if OPT_LD_LATENCY2
            svt_release_mutex(enc_ctx->ref_pic_list_mutex);
#endif
            // Update the last decode order
            if (input_pic_demux->decode_order == decode_order)
                decode_order++;
            break;
        default:
            scs     = input_pic_demux->scs;
            enc_ctx = scs->enc_ctx;

            CHECK_REPORT_ERROR_NC(enc_ctx->app_callback_ptr, EB_ENC_PM_ERROR9);

            pcs     = (PictureParentControlSet *)NULL;
            enc_ctx = (EncodeContext *)NULL;

            break;
        }

        // ***********************************
        //  Common Code
        // *************************************

        // Walk the input queue and start all ready pictures.  Mark entry as null after started.
        // Increment the head as you go.
        if (enc_ctx != (EncodeContext *)NULL) {
            input_queue_index = enc_ctx->input_picture_queue_head_index;
            while (input_queue_index != enc_ctx->input_picture_queue_tail_index) {
                input_entry = enc_ctx->input_picture_queue[input_queue_index];

                if (input_entry->input_object_ptr != NULL) {
                    entry_ppcs        = (PictureParentControlSet *)input_entry->input_object_ptr->object_ptr;
                    entry_scs_ptr     = entry_ppcs->scs;
                    availability_flag = TRUE;
                    if (entry_ppcs->decode_order != decode_order && (scs->enable_dec_order))
                        availability_flag = FALSE;

                    // pic mgr starts pictures in dec order (no need to wait for feedback)
                    if (entry_scs_ptr->enable_pic_mgr_dec_order)
                        if (entry_ppcs->picture_number > 0 &&
                            entry_ppcs->decode_order != context_ptr->pmgr_dec_order + 1)
                            availability_flag = FALSE;
                    // Only start pictures that can be reached with the given number of reference
                    // buffers. PA may have more ref buffers than PM, so can send pictures faster,
                    // but PM can't start those pictures until it has sufficient buffers to reach
                    // it.
                    if (entry_ppcs->decode_order >
                        (context_ptr->consecutive_dec_order + scs->reference_picture_buffer_init_count - REF_FRAMES))
                        availability_flag = FALSE;
                    if ((entry_ppcs->slice_type == P_SLICE) || (entry_ppcs->slice_type == B_SLICE)) {
                        uint8_t max_ref_count = (entry_ppcs->slice_type == B_SLICE) ? ALT + 1
                                                                                    : BWD; // no list1 refs for P_SLICE
                        for (REF_FRAME_MINUS1 ref = LAST; ref < max_ref_count; ref++) {
                            // hardcode the reference for the overlay frame
                            uint64_t ref_poc = entry_ppcs->is_overlay ? entry_ppcs->picture_number
                                                                      : entry_ppcs->av1_ref_signal.ref_poc_array[ref];

                            uint8_t list_idx = get_list_idx(ref + 1);
                            uint8_t ref_idx  = get_ref_frame_idx(ref + 1);
                            assert(IMPLIES(entry_ppcs->is_overlay, list_idx == 0));
                            if ((list_idx == 0 && ref_idx >= entry_ppcs->ref_list0_count) ||
                                (list_idx == 1 && ref_idx >= entry_ppcs->ref_list1_count))
                                continue;
#if OPT_LD_LATENCY2
                            svt_block_on_mutex(enc_ctx->ref_pic_list_mutex);
#endif
                            ref_entry = search_ref_in_ref_queue(enc_ctx, ref_poc);

                            if (ref_entry != NULL) {
                                availability_flag = (availability_flag == FALSE) ? FALSE
                                                                                 : // Don't update if already False
                                    (scs->static_config.rate_control_mode && entry_ppcs->slice_type != I_SLICE &&
                                     entry_ppcs->temporal_layer_index == 0 && !ref_entry->feedback_arrived &&
                                     !enc_ctx->terminating_sequence_flag_received)
                                    ? FALSE
                                    : (entry_ppcs->frame_end_cdf_update_mode && !ref_entry->frame_context_updated)
                                    ? FALSE
                                    : (ref_entry->reference_available) ? TRUE
                                                                       : // The Reference has been completed
                                    FALSE; // The Reference has not been completed
                            } else {
                                availability_flag = FALSE;
                            }
#if OPT_LD_LATENCY2
                            svt_release_mutex(enc_ctx->ref_pic_list_mutex);
#endif
                        }
                    }

                    if (availability_flag == TRUE) {
                        if (entry_ppcs->is_ref) {
                            EbObjectWrapper *ref_pic_wrapper;
                            // Get Empty Reference Picture Object
                            svt_get_empty_object(scs->enc_ctx->reference_picture_pool_fifo_ptr, &ref_pic_wrapper);
                            entry_ppcs->ref_pic_wrapper = ref_pic_wrapper;
                            // reset reference object in case of its members are altered by superres
                            // tool
                            EbReferenceObject *ref = (EbReferenceObject *)ref_pic_wrapper->object_ptr;
                            svt_reference_object_reset(ref, scs);
                            // Give the new Reference a nominal live_count of 1
                            svt_object_inc_live_count(entry_ppcs->ref_pic_wrapper, 1);
#if SRM_REPORT
                            pcs->ref_pic_wrapper->pic_number = pcs->picture_number;
#endif
                        } else {
                            entry_ppcs->ref_pic_wrapper = NULL;
                        }
                        // Get New  Empty recon-coef from recon-coef  Pool
                        svt_get_empty_object(context_ptr->recon_coef_fifo_ptr, &enc_dec_wrapper);
                        // Child PCS is released by Packetization
                        svt_object_inc_live_count(enc_dec_wrapper, 1);
                        enc_dec_ptr                  = (EncDecSet *)enc_dec_wrapper->object_ptr;
                        enc_dec_ptr->enc_dec_wrapper = enc_dec_wrapper;

                        // 1.Link The Child PCS to its Parent
                        enc_dec_ptr->ppcs_wrapper = input_entry->input_object_ptr;
                        enc_dec_ptr->ppcs         = entry_ppcs;

                        enc_dec_ptr->ppcs->enc_dec_ptr = enc_dec_ptr;
                        // Get New  Empty Child PCS from PCS Pool
                        svt_get_empty_object(context_ptr->picture_control_set_fifo_ptr, &child_pcs_wrapper);

                        // Child PCS is released by Packetization
                        svt_object_inc_live_count(child_pcs_wrapper, 1);

                        child_pcs = (PictureControlSet *)child_pcs_wrapper->object_ptr;

                        child_pcs->c_pcs_wrapper_ptr = child_pcs_wrapper;

                        // 1.Link The Child PCS to its Parent
                        child_pcs->ppcs_wrapper = input_entry->input_object_ptr;
                        child_pcs->ppcs         = entry_ppcs;

                        child_pcs->ppcs->child_pcs = child_pcs;
                        // 1b Link The Child PCS to av1_cm to be used by Restoration
                        child_pcs->ppcs->av1_cm->child_pcs = child_pcs;

                        // 2. Have some common information between  ChildPCS and ParentPCS.
                        child_pcs->scs                  = entry_ppcs->scs;
                        child_pcs->picture_qp           = entry_ppcs->picture_qp;
                        child_pcs->picture_number       = entry_ppcs->picture_number;
                        child_pcs->slice_type           = entry_ppcs->slice_type;
                        child_pcs->temporal_layer_index = entry_ppcs->temporal_layer_index;

                        child_pcs->ppcs->total_num_bits = 0;
                        child_pcs->ppcs->picture_qp     = entry_ppcs->picture_qp;
                        child_pcs->enc_mode             = entry_ppcs->enc_mode;
                        child_pcs->b64_total_count      = entry_ppcs->b64_total_count;

                        child_pcs->enc_dec_coded_sb_count = 0;
                        child_pcs->hbd_md                 = entry_ppcs->hbd_md;
                        context_ptr->pmgr_dec_order       = child_pcs->ppcs->decode_order;

                        // Update the consecutive decode order count, if this picture is the next
                        // picture in decode order. Otherwise, add the picture to the
                        // started_pics_dec_order list so the consecutive decode order count can be
                        // properly updated later.
                        if (entry_ppcs->decode_order == context_ptr->consecutive_dec_order + 1) {
                            context_ptr->consecutive_dec_order++;

                            if (context_ptr->started_pics_dec_order_head_idx !=
                                context_ptr->started_pics_dec_order_tail_idx) {
                                for (int idx = context_ptr->started_pics_dec_order_head_idx;
                                     idx != context_ptr->started_pics_dec_order_tail_idx;) {
                                    if (context_ptr->started_pics_dec_order[idx] ==
                                        context_ptr->consecutive_dec_order + 1) {
                                        context_ptr->consecutive_dec_order++;
                                        idx = context_ptr->started_pics_dec_order_head_idx;
                                    } else {
                                        idx = (idx == REFERENCE_QUEUE_MAX_DEPTH - 1) ? 0 : idx + 1;
                                    }
                                }

                                /* clang-format off */
                                while (context_ptr->started_pics_dec_order_head_idx != context_ptr->started_pics_dec_order_tail_idx &&
                                       context_ptr->started_pics_dec_order[context_ptr->started_pics_dec_order_head_idx] <= context_ptr->consecutive_dec_order) {
                                    context_ptr->started_pics_dec_order_head_idx =
                                        (context_ptr->started_pics_dec_order_head_idx == REFERENCE_QUEUE_MAX_DEPTH - 1)
                                        ? 0
                                        : context_ptr->started_pics_dec_order_head_idx + 1;
                                }
                                /* clang-format on */
                            }
                        } else if (entry_ppcs->decode_order > 0) {
                            context_ptr->started_pics_dec_order[context_ptr->started_pics_dec_order_tail_idx] =
                                entry_ppcs->decode_order;
                            context_ptr->started_pics_dec_order_tail_idx =
                                (context_ptr->started_pics_dec_order_tail_idx == REFERENCE_QUEUE_MAX_DEPTH - 1)
                                ? 0
                                : context_ptr->started_pics_dec_order_tail_idx + 1;
                        }
                        // 3.make all  init for ChildPCS
                        pic_width_in_sb = (entry_ppcs->aligned_width + entry_scs_ptr->sb_size - 1) /
                            entry_scs_ptr->sb_size;
                        picture_height_in_sb = (entry_ppcs->aligned_height + entry_scs_ptr->sb_size - 1) /
                            entry_scs_ptr->sb_size;

                        svt_aom_init_enc_dec_segement(entry_ppcs);

                        int                             sb_size_log2 = entry_scs_ptr->seq_header.sb_size_log2;
                        struct PictureParentControlSet *ppcs         = child_pcs->ppcs;
                        const int                       tile_cols    = ppcs->av1_cm->tiles_info.tile_cols;
                        const int                       tile_rows    = ppcs->av1_cm->tiles_info.tile_rows;
                        Av1Common *const                cm           = ppcs->av1_cm;
                        uint16_t                        tile_row, tile_col;
                        uint32_t                        x_sb_index, y_sb_index;
                        TileInfo                        tile_info;

                        child_pcs->sb_total_count = pic_width_in_sb * picture_height_in_sb;

                        // force re-ctor sb_ptr since child_pcs_ptrs are reused, and sb_ptr could be
                        // altered by superres tool when coding previous pictures
                        if (scs->static_config.superres_mode > SUPERRES_NONE ||
                            scs->static_config.resize_mode > RESIZE_NONE) {
                            // Modify sb_prt_array in child pcs
                            uint16_t sb_index;
                            uint16_t sb_origin_x = 0;
                            uint16_t sb_origin_y = 0;
                            for (sb_index = 0; sb_index < child_pcs->sb_total_count; ++sb_index) {
                                svt_aom_largest_coding_unit_dctor(child_pcs->sb_ptr_array[sb_index]);
                                svt_aom_largest_coding_unit_ctor(child_pcs->sb_ptr_array[sb_index],
                                                                 (uint8_t)scs->sb_size,
                                                                 (uint16_t)(sb_origin_x * scs->sb_size),
                                                                 (uint16_t)(sb_origin_y * scs->sb_size),
                                                                 (uint16_t)sb_index,
                                                                 child_pcs->enc_mode,
                                                                 scs->max_block_cnt,
                                                                 child_pcs);
                                // Increment the Order in coding order (Raster Scan Order)
                                sb_origin_y = (sb_origin_x == pic_width_in_sb - 1) ? sb_origin_y + 1 : sb_origin_y;
                                sb_origin_x = (sb_origin_x == pic_width_in_sb - 1) ? 0 : sb_origin_x + 1;
                            }

                            // reset input_frame16bit to align with enhanced_pic
                            if (scs->static_config.encoder_bit_depth > EB_EIGHT_BIT) {
                                svt_aom_copy_buffer_info(entry_ppcs->enhanced_pic, child_pcs->input_frame16bit);
                            }
                        }

                        // Update pcs->mi_stride
                        child_pcs->mi_stride = pic_width_in_sb * (scs->sb_size >> MI_SIZE_LOG2);

                        // copy buffer info from the downsampled picture to the input frame 16 bit
                        // buffer
                        if ((entry_ppcs->frame_superres_enabled || entry_ppcs->frame_resize_enabled) &&
                            scs->static_config.encoder_bit_depth > EB_EIGHT_BIT) {
                            svt_aom_copy_buffer_info(entry_ppcs->enhanced_downscaled_pic, child_pcs->input_frame16bit);
                        }

                        //                        child_pcs->ppcs->av1_cm->pcs = child_pcs;
                        // Palette
                        rtime_alloc_palette_tokens(scs, child_pcs);
                        TOKENEXTRA  *pre_tok  = child_pcs->tile_tok[0][0];
                        unsigned int tile_tok = 0;
                        // Tile Loop
                        for (tile_row = 0; tile_row < tile_rows; tile_row++) {
                            svt_av1_tile_set_row(&tile_info, &cm->tiles_info, cm->mi_rows, tile_row);

                            for (tile_col = 0; tile_col < tile_cols; tile_col++) {
                                svt_av1_tile_set_col(&tile_info, &cm->tiles_info, cm->mi_cols, tile_col);
                                tile_info.tile_rs_index = tile_col + tile_row * tile_cols;

                                // Palette
                                if (child_pcs->tile_tok[0][0]) {
                                    child_pcs->tile_tok[tile_row][tile_col] = pre_tok + tile_tok;
                                    pre_tok                                 = child_pcs->tile_tok[tile_row][tile_col];
                                    int tile_mb_rows = (tile_info.mi_row_end - tile_info.mi_row_start + 2) >> 2;
                                    int tile_mb_cols = (tile_info.mi_col_end - tile_info.mi_col_start + 2) >> 2;
                                    tile_tok = get_token_alloc(tile_mb_rows, tile_mb_cols, (sb_size_log2 + 2), 2);
                                }
                                for ((y_sb_index = cm->tiles_info.tile_row_start_mi[tile_row] >> sb_size_log2);
                                     (y_sb_index <
                                      ((uint32_t)cm->tiles_info.tile_row_start_mi[tile_row + 1] >> sb_size_log2));
                                     y_sb_index++) {
                                    for ((x_sb_index = cm->tiles_info.tile_col_start_mi[tile_col] >> sb_size_log2);
                                         (x_sb_index <
                                          ((uint32_t)cm->tiles_info.tile_col_start_mi[tile_col + 1] >> sb_size_log2));
                                         x_sb_index++) {
                                        int sb_index = (uint16_t)(x_sb_index + y_sb_index * pic_width_in_sb);
                                        child_pcs->sb_ptr_array[sb_index]->tile_info = tile_info;
                                    }
                                }
                            }
                        }
                        cm->mi_stride = child_pcs->mi_stride;
                        // Reset the Reference Lists
                        EB_MEMSET(child_pcs->ref_pic_ptr_array[REF_LIST_0],
                                  0,
                                  REF_LIST_MAX_DEPTH * sizeof(EbObjectWrapper *));
                        EB_MEMSET(child_pcs->ref_pic_ptr_array[REF_LIST_1],
                                  0,
                                  REF_LIST_MAX_DEPTH * sizeof(EbObjectWrapper *));
                        EB_MEMSET(child_pcs->ref_pic_qp_array[REF_LIST_0], 0, REF_LIST_MAX_DEPTH * sizeof(uint8_t));
                        EB_MEMSET(child_pcs->ref_pic_qp_array[REF_LIST_1], 0, REF_LIST_MAX_DEPTH * sizeof(uint8_t));
                        EB_MEMSET(
                            child_pcs->ref_slice_type_array[REF_LIST_0], 0, REF_LIST_MAX_DEPTH * sizeof(SliceType));
                        EB_MEMSET(
                            child_pcs->ref_slice_type_array[REF_LIST_1], 0, REF_LIST_MAX_DEPTH * sizeof(SliceType));
                        EB_MEMSET(child_pcs->ref_pic_r0[REF_LIST_0], 0, REF_LIST_MAX_DEPTH * sizeof(double));
                        EB_MEMSET(child_pcs->ref_pic_r0[REF_LIST_1], 0, REF_LIST_MAX_DEPTH * sizeof(double));
                        int8_t ref_index = 0;
                        if ((entry_ppcs->slice_type == P_SLICE) || (entry_ppcs->slice_type == B_SLICE)) {
                            int8_t  max_temporal_index = -1;
                            uint8_t max_ref_count      = (entry_ppcs->slice_type == B_SLICE)
                                     ? ALT + 1
                                     : BWD; // no list1 refs for P_SLICE
                            for (REF_FRAME_MINUS1 ref = LAST; ref < max_ref_count; ref++) {
                                // hardcode the reference for the overlay frame
                                uint64_t ref_poc = entry_ppcs->is_overlay
                                    ? entry_ppcs->picture_number
                                    : entry_ppcs->av1_ref_signal.ref_poc_array[ref];

                                uint8_t list_idx = get_list_idx(ref + 1);
                                uint8_t ref_idx  = get_ref_frame_idx(ref + 1);
                                assert(IMPLIES(entry_ppcs->is_overlay, list_idx == 0));
                                if ((list_idx == 0 && ref_idx >= entry_ppcs->ref_list0_count) ||
                                    (list_idx == 1 && ref_idx >= entry_ppcs->ref_list1_count))
                                    continue;

                                ref_entry = search_ref_in_ref_queue(enc_ctx, ref_poc);
                                assert(ref_entry != 0);
                                CHECK_REPORT_ERROR((ref_entry), enc_ctx->app_callback_ptr, EB_ENC_PM_ERROR10);

                                /* clang-format off */
                                if (entry_ppcs->frame_end_cdf_update_mode) {
                                    child_pcs->ref_frame_context[svt_get_ref_frame_type(list_idx, ref_idx) - LAST_FRAME] =
                                        ((EbReferenceObject*)ref_entry->reference_object_ptr->object_ptr)->frame_context;
                                    if (max_temporal_index < (int8_t)ref_entry->temporal_layer_index &&
                                        (int8_t)ref_entry->temporal_layer_index <= child_pcs->temporal_layer_index) {
                                        max_temporal_index = (int8_t)ref_entry->temporal_layer_index;
                                        ref_index = svt_get_ref_frame_type(list_idx, ref_idx) - LAST_FRAME;
                                        for (int frame = LAST_FRAME; frame <= ALTREF_FRAME; ++frame) {
                                            EbReferenceObject *ref_obj =
                                                (EbReferenceObject *)ref_entry->reference_object_ptr->object_ptr;

                                            child_pcs->ref_global_motion[frame] =
                                                ref_obj->slice_type != I_SLICE
                                                ? ref_obj->global_motion[frame]
                                                : default_warp_params;
                                        }
                                    }
                                }
                                /* clang-format on */
                                // Set the Reference Object
                                child_pcs->ref_pic_ptr_array[list_idx][ref_idx] = ref_entry->reference_object_ptr;
                                child_pcs->ref_pic_qp_array[list_idx][ref_idx] =
                                    (uint8_t)((EbReferenceObject *)ref_entry->reference_object_ptr->object_ptr)->qp;
                                child_pcs->ref_slice_type_array[list_idx][ref_idx] =
                                    ((EbReferenceObject *)ref_entry->reference_object_ptr->object_ptr)->slice_type;
                                child_pcs->ref_pic_r0[list_idx][ref_idx] =
                                    ((EbReferenceObject *)ref_entry->reference_object_ptr->object_ptr)->r0;
                                // Increment the Reference's liveCount by the number of tiles in the
                                // input picture
                                svt_object_inc_live_count(ref_entry->reference_object_ptr, 1);

#if DEBUG_SFRAME
                                if (scs->static_config.pass != ENC_FIRST_PASS) {
                                    fprintf(stderr,
                                            "\nframe %d, layer %d, ref-list-0 count %u, ref "
                                            "frame %d, ref frame remain dep count %d\n",
                                            (int)child_pcs->picture_number,
                                            entry_ppcs->temporal_layer_index,
                                            entry_ppcs->ref_list0_count,
                                            (int)child_pcs->ppcs->ref_pic_poc_array[REF_LIST_0][ref_idx],
                                            ref_entry->dependent_count);
                                }
#endif
                            }

                            /* clang-format off */
                            //fill the non used spots to be used in MFMV.
                            for (uint8_t ref_idx = entry_ppcs->ref_list0_count; ref_idx < 4; ++ref_idx)
                                child_pcs->ref_pic_ptr_array[REF_LIST_0][ref_idx] =
                                    child_pcs->ref_pic_ptr_array[REF_LIST_0][0];

                            if (entry_ppcs->ref_list1_count == 0) {
                                for (uint8_t ref_idx = entry_ppcs->ref_list1_count; ref_idx < 3; ++ref_idx) {
                                    child_pcs->ref_pic_ptr_array[REF_LIST_1][ref_idx] =
                                        child_pcs->ref_pic_ptr_array[REF_LIST_0][0];
                                    child_pcs->ref_pic_qp_array[REF_LIST_1][ref_idx] =
                                        child_pcs->ref_pic_qp_array[REF_LIST_0][0];
                                    child_pcs->ref_slice_type_array[REF_LIST_1][ref_idx] =
                                        child_pcs->ref_slice_type_array[REF_LIST_0][0];
                                    child_pcs->ref_pic_r0[REF_LIST_1][ref_idx] =
                                        child_pcs->ref_pic_r0[REF_LIST_0][0];
                                }
                            }

                            if (entry_ppcs->slice_type == B_SLICE) {
                                //fill the non used spots to be used in MFMV.
                                if (entry_ppcs->ref_list1_count) {
                                    for (uint8_t ref_idx = entry_ppcs->ref_list1_count; ref_idx < 3; ++ref_idx)
                                        child_pcs->ref_pic_ptr_array[REF_LIST_1][ref_idx] =
                                            child_pcs->ref_pic_ptr_array[REF_LIST_1][0];
                                }
                            }
                            /* clang-format on */
                        }

                        if (entry_ppcs->frame_end_cdf_update_mode) {
                            if (entry_ppcs->slice_type != I_SLICE && entry_ppcs->frm_hdr.frame_type != S_FRAME)
                                child_pcs->ppcs->frm_hdr.primary_ref_frame = ref_index;
                            else
                                child_pcs->ppcs->frm_hdr.primary_ref_frame = PRIMARY_REF_NONE;
                            child_pcs->ppcs->refresh_frame_context = REFRESH_FRAME_CONTEXT_BACKWARD;

                        } else {
                            child_pcs->ppcs->frm_hdr.primary_ref_frame = PRIMARY_REF_NONE;
                            // Tells each frame whether to forward their data to the next frames;
                            // never disabled so that the feature can be on in higher layers, while
                            // off in low layers.
                            child_pcs->ppcs->refresh_frame_context = REFRESH_FRAME_CONTEXT_BACKWARD;
                        }
                        EbObjectWrapper *out_results_wrapper;
                        // Get Empty Results Object
                        svt_get_empty_object(context_ptr->picture_manager_output_fifo_ptr, &out_results_wrapper);
                        RateControlTasks *rc_tasks = (RateControlTasks *)out_results_wrapper->object_ptr;
                        rc_tasks->pcs_wrapper      = child_pcs->c_pcs_wrapper_ptr;
                        rc_tasks->task_type        = RC_INPUT;

                        // printf("picMgr sending:%x \n", rc_tasks->pcs_wrapper);

                        // Post the Full Results Object
                        svt_post_full_object(out_results_wrapper);
                        // Remove the Input Entry from the Input Queue
                        input_entry->input_object_ptr = (EbObjectWrapper *)NULL;
#if OPT_LD_LATENCY2
                        svt_block_on_mutex(enc_ctx->ref_pic_list_mutex);
#endif
                        for (uint32_t i = 0; i < enc_ctx->ref_pic_list_length; i++) {
                            ref_entry = enc_ctx->ref_pic_list[i];
                            if (ref_entry->is_valid) {
                                for (uint8_t idx = 0; idx < entry_ppcs->released_pics_count; idx++) {
                                    if (ref_entry->decode_order == entry_ppcs->released_pics[idx]) {
                                        ref_entry->dec_order_of_last_ref = entry_ppcs->decode_order;
                                    }
                                }
                            }
                        }
#if OPT_LD_LATENCY2
                        svt_release_mutex(enc_ctx->ref_pic_list_mutex);
#endif
                    }
                }

                // Increment the head_index if the head is null
                enc_ctx->input_picture_queue_head_index =
                    (enc_ctx->input_picture_queue[enc_ctx->input_picture_queue_head_index]->input_object_ptr)
                    ? enc_ctx->input_picture_queue_head_index
                    : (enc_ctx->input_picture_queue_head_index == INPUT_QUEUE_MAX_DEPTH - 1)
                    ? 0
                    : enc_ctx->input_picture_queue_head_index + 1;

                // Increment the input_queue_index Iterator
                input_queue_index = (input_queue_index == INPUT_QUEUE_MAX_DEPTH - 1) ? 0 : input_queue_index + 1;
            }
#if OPT_LD_LATENCY2
            svt_block_on_mutex(enc_ctx->ref_pic_list_mutex);
#endif
            for (uint32_t i = 0; i < enc_ctx->ref_pic_list_length; i++) {
                ref_entry = enc_ctx->ref_pic_list[i];
                if (ref_entry->is_valid) {
                    // Remove the entry & release the reference if there are no remaining references
                    if (ref_entry->dec_order_of_last_ref <= context_ptr->consecutive_dec_order &&
                        ref_entry->release_enable && ref_entry->reference_available &&
                        ref_entry->reference_object_ptr &&
                        (!ref_entry->frame_end_cdf_update_required || ref_entry->frame_context_updated)) {
                        // Release the nominal live_count value
                        svt_release_object(ref_entry->reference_object_ptr);
                        ref_entry->reference_object_ptr  = (EbObjectWrapper *)NULL;
                        ref_entry->reference_available   = FALSE;
                        ref_entry->is_ref                = FALSE;
                        ref_entry->is_valid              = false;
                        ref_entry->frame_context_updated = FALSE;
                        ref_entry->feedback_arrived      = FALSE;
                        svt_post_semaphore(scs->ref_buffer_available_semaphore);
                    }
                }
            }
#if OPT_LD_LATENCY2
            svt_release_mutex(enc_ctx->ref_pic_list_mutex);
#endif
        }

        // Release the Input Picture Demux Results
        svt_release_object(input_pic_demux_wrapper);
    }
    return NULL;
}
