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
void rtime_alloc_palette_tokens(SequenceControlSet *scs_ptr, PictureControlSet *child_pcs_ptr);
extern MvReferenceFrame svt_get_ref_frame_type(uint8_t list, uint8_t ref_idx);
/************************************************
 * Defines
 ************************************************/
#define POC_CIRCULAR_ADD(base, offset) (((base) + (offset)))

void largest_coding_unit_dctor(EbPtr p);

/************************************************
  * Configure Picture edges
  ************************************************/
static void configure_picture_edges(SequenceControlSet *scs_ptr, PictureControlSet *ppsPtr) {
    // Tiles Initialisation
    const uint16_t pic_width_in_sb = (ppsPtr->parent_pcs_ptr->aligned_width + scs_ptr->sb_size -
                                      1) /
        scs_ptr->sb_size;
    const uint16_t picture_height_in_sb = (ppsPtr->parent_pcs_ptr->aligned_height +
                                           scs_ptr->sb_size - 1) /
        scs_ptr->sb_size;
    unsigned x_sb_index, y_sb_index, sb_index;

    // SB-loops
    for (y_sb_index = 0; y_sb_index < picture_height_in_sb; ++y_sb_index) {
        for (x_sb_index = 0; x_sb_index < pic_width_in_sb; ++x_sb_index) {
            sb_index = (uint16_t)(x_sb_index + y_sb_index * pic_width_in_sb);
            ppsPtr->sb_ptr_array[sb_index]->picture_left_edge_flag = (x_sb_index == 0) ? TRUE
                                                                                       : FALSE;
            ppsPtr->sb_ptr_array[sb_index]->picture_top_edge_flag  = (y_sb_index == 0) ? TRUE
                                                                                       : FALSE;
            ppsPtr->sb_ptr_array[sb_index]->picture_right_edge_flag =
                (x_sb_index == (unsigned)(pic_width_in_sb - 1)) ? TRUE : FALSE;
        }
    }

    return;
}
void write_stat_to_file(SequenceControlSet *scs_ptr, StatStruct stat_struct, uint64_t ref_poc);

static void picture_manager_context_dctor(EbPtr p) {
    EbThreadContext       *thread_context_ptr = (EbThreadContext *)p;
    PictureManagerContext *obj                = (PictureManagerContext *)thread_context_ptr->priv;
    EB_FREE_ARRAY(obj);
}
/************************************************
 * Picture Manager Context Constructor
 ************************************************/
EbErrorType picture_manager_context_ctor(EbThreadContext   *thread_context_ptr,
                                         const EbEncHandle *enc_handle_ptr,
                                         int                rate_control_index) {
    PictureManagerContext *context_ptr;
    EB_CALLOC_ARRAY(context_ptr, 1);
    thread_context_ptr->priv  = context_ptr;
    thread_context_ptr->dctor = picture_manager_context_dctor;

    context_ptr->picture_input_fifo_ptr = svt_system_resource_get_consumer_fifo(
        enc_handle_ptr->picture_demux_results_resource_ptr, 0);
    context_ptr->picture_manager_output_fifo_ptr = svt_system_resource_get_producer_fifo(
        enc_handle_ptr->rate_control_tasks_resource_ptr, rate_control_index);
    context_ptr->picture_control_set_fifo_ptr = svt_system_resource_get_producer_fifo(
        enc_handle_ptr->picture_control_set_pool_ptr_array[0], 0); //The Child PCS Pool here
    context_ptr->recon_coef_fifo_ptr = svt_system_resource_get_producer_fifo(
        enc_handle_ptr->enc_dec_pool_ptr_array[0], 0); //The Child PCS Pool here

    return EB_ErrorNone;
}

void copy_buffer_info(EbPictureBufferDesc *src_ptr, EbPictureBufferDesc *dst_ptr) {
    dst_ptr->width             = src_ptr->width;
    dst_ptr->height            = src_ptr->height;
    dst_ptr->max_width         = src_ptr->max_width;
    dst_ptr->max_height        = src_ptr->max_height;
    dst_ptr->stride_y          = src_ptr->stride_y;
    dst_ptr->stride_cb         = src_ptr->stride_cb;
    dst_ptr->stride_cr         = src_ptr->stride_cr;
    dst_ptr->origin_x          = src_ptr->origin_x;
    dst_ptr->origin_bot_y      = src_ptr->origin_bot_y;
    dst_ptr->origin_y          = src_ptr->origin_y;
    dst_ptr->stride_bit_inc_y  = src_ptr->stride_bit_inc_y;
    dst_ptr->stride_bit_inc_cb = src_ptr->stride_bit_inc_cb;
    dst_ptr->stride_bit_inc_cr = src_ptr->stride_bit_inc_cr;
    dst_ptr->luma_size         = src_ptr->luma_size;
    dst_ptr->chroma_size       = src_ptr->chroma_size;
}

void set_tile_info(PictureParentControlSet *pcs_ptr);
/*
  walk the ref queue, and looks for a particular pic
  return NULL if not found.
*/
ReferenceQueueEntry *search_ref_in_ref_queue(EncodeContext *encode_context_ptr, uint64_t ref_poc) {
    ReferenceQueueEntry *ref_entry_ptr = NULL;
    uint32_t             ref_queue_i   = encode_context_ptr->reference_picture_queue_head_index;
    // Find the Reference in the Reference Queue
    do {
        ref_entry_ptr = encode_context_ptr->reference_picture_queue[ref_queue_i];
        if (ref_entry_ptr->picture_number == ref_poc)
            return ref_entry_ptr;

        // Increment the reference_queue_index Iterator
        ref_queue_i = (ref_queue_i == REFERENCE_QUEUE_MAX_DEPTH - 1) ? 0 : ref_queue_i + 1;

    } while (ref_queue_i != encode_context_ptr->reference_picture_queue_tail_index);

    return NULL;
}
/*
  update dependent count of any existing pic in the ref queue.
  if the picture is not in the ref queue, subesequent calls
  will attempt to do the cleaning
*/
void clean_pictures_in_ref_queue(EncodeContext *ctx) {
    PicQueueEntry      **queue         = ctx->dep_cnt_picture_queue;
    ReferenceQueueEntry *ref_entry_ptr = NULL;

    //all pictures needing the clean-up are stored in dep cnt queue
    //go through all of elements, search it in the ref queue, and clean
    //it if found
    uint32_t q_idx = ctx->dep_q_head;
    while (q_idx != ctx->dep_q_tail) {
        if (queue[q_idx]->is_done == 0) {
            ref_entry_ptr = search_ref_in_ref_queue(ctx, queue[q_idx]->pic_num);
            if (ref_entry_ptr != NULL) {
                int new_dep_cnt = (int)ref_entry_ptr->dependent_count + queue[q_idx]->dep_cnt_diff;
                if (new_dep_cnt < 0)
                    SVT_ERROR(" PicMgr: Negative dep count\n");
                //printf(" search %I64i  ====Pic-MGR-Dep-Cnt-reduce found: %I64i  %i --> %i\n",
                //    queue[q_idx]->pic_num, ref_entry_ptr->picture_number, ref_entry_ptr->dependent_count, new_dep_cnt);
                ref_entry_ptr->dependent_count = new_dep_cnt;
                queue[q_idx]->is_done          = 1;
            }
        }

        // Increment the head_index if the head is done
        if (queue[ctx->dep_q_head]->is_done)
            ctx->dep_q_head = (ctx->dep_q_head == REFERENCE_QUEUE_MAX_DEPTH - 1)
                ? 0
                : ctx->dep_q_head + 1;

        //Increment the queue_index Iterator
        q_idx = (q_idx == REFERENCE_QUEUE_MAX_DEPTH - 1) ? 0 : q_idx + 1;
    }
}

/*
  Pic Mgr maintains a queue of pic to clear dep cnt
  if a picture comes from PD, it copies its clean-up list to the queue
*/
void copy_dep_cnt_cleaning_list(EncodeContext *ctx, PictureParentControlSet *pcs) {
    //triggering picture hands list to PicMgr, which will do the clearing
    PicQueueEntry **queue = ctx->dep_cnt_picture_queue;

    for (uint32_t pic_i = 0; pic_i < pcs->other_updated_links_cnt; pic_i++) {
        //place pic to be cleared in queue
        //the spot should be empty
        if (queue[ctx->dep_q_tail]->is_done == 0)
            SVT_ERROR("\n pic %I64i  ERR  DepCount Queue is done\n",
                      queue[ctx->dep_q_tail]->pic_num);

        queue[ctx->dep_q_tail]->pic_num      = pcs->updated_links_arr[pic_i].pic_num;
        queue[ctx->dep_q_tail]->dep_cnt_diff = pcs->updated_links_arr[pic_i].dep_cnt_diff;
        queue[ctx->dep_q_tail]->is_done      = 0;

        //advance the queue
        ctx->dep_q_tail = (ctx->dep_q_tail == REFERENCE_QUEUE_MAX_DEPTH - 1) ? 0
                                                                             : ctx->dep_q_tail + 1;
    }
}

void init_enc_dec_segement(PictureParentControlSet *parentpicture_control_set_ptr) {
    SequenceControlSet *scs_ptr         = parentpicture_control_set_ptr->scs_ptr;
    uint8_t             pic_width_in_sb = (uint8_t)((parentpicture_control_set_ptr->aligned_width +
                                         scs_ptr->sb_size - 1) /
                                        scs_ptr->sb_size);
    uint8_t picture_height_in_sb        = (uint8_t)((parentpicture_control_set_ptr->aligned_height +
                                              scs_ptr->sb_size - 1) /
                                             scs_ptr->sb_size);
    set_tile_info(parentpicture_control_set_ptr);
    int      sb_size_log2 = scs_ptr->seq_header.sb_size_log2;
    uint32_t enc_dec_seg_col_cnt =
        scs_ptr
            ->enc_dec_segment_col_count_array[parentpicture_control_set_ptr->temporal_layer_index];
    uint32_t enc_dec_seg_row_cnt =
        scs_ptr
            ->enc_dec_segment_row_count_array[parentpicture_control_set_ptr->temporal_layer_index];
    Av1Common *const cm              = parentpicture_control_set_ptr->av1_cm;
    const int        tile_cols       = parentpicture_control_set_ptr->av1_cm->tiles_info.tile_cols;
    const int        tile_rows       = parentpicture_control_set_ptr->av1_cm->tiles_info.tile_rows;
    uint8_t          tile_group_cols = MIN(
        tile_cols,
        scs_ptr->tile_group_col_count_array[parentpicture_control_set_ptr->temporal_layer_index]);
    uint8_t tile_group_rows = MIN(
        tile_rows,
        scs_ptr->tile_group_row_count_array[parentpicture_control_set_ptr->temporal_layer_index]);

    if (tile_group_cols * tile_group_rows > 1) {
        enc_dec_seg_col_cnt = MIN(enc_dec_seg_col_cnt,
                                  (uint8_t)(pic_width_in_sb / tile_group_cols));
        enc_dec_seg_row_cnt = MIN(enc_dec_seg_row_cnt,
                                  (uint8_t)(picture_height_in_sb / tile_group_rows));
    }
    parentpicture_control_set_ptr->tile_group_cols = tile_group_cols;
    parentpicture_control_set_ptr->tile_group_rows = tile_group_rows;

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

            TileGroupInfo *tg_info_ptr =
                &parentpicture_control_set_ptr->tile_group_info[tile_group_idx];

            tg_info_ptr->tile_group_tile_start_x = top_left_tile_col_idx;
            tg_info_ptr->tile_group_tile_end_x   = bottom_right_tile_col_idx;

            tg_info_ptr->tile_group_tile_start_y = top_left_tile_row_idx;
            tg_info_ptr->tile_group_tile_end_y   = bottom_right_tile_row_idx;

            tg_info_ptr->tile_group_sb_start_x =
                cm->tiles_info.tile_col_start_mi[top_left_tile_col_idx] >> sb_size_log2;
            tg_info_ptr->tile_group_sb_start_y =
                cm->tiles_info.tile_row_start_mi[top_left_tile_row_idx] >> sb_size_log2;

            // Get the SB end of the bottom right tile
            tg_info_ptr->tile_group_sb_end_x =
                (cm->tiles_info.tile_col_start_mi[bottom_right_tile_col_idx] >> sb_size_log2);
            tg_info_ptr->tile_group_sb_end_y =
                (cm->tiles_info.tile_row_start_mi[bottom_right_tile_row_idx] >> sb_size_log2);

            // Get the width/height of tile group in SB
            tg_info_ptr->tile_group_height_in_sb = tg_info_ptr->tile_group_sb_end_y -
                tg_info_ptr->tile_group_sb_start_y;
            tg_info_ptr->tile_group_width_in_sb = tg_info_ptr->tile_group_sb_end_x -
                tg_info_ptr->tile_group_sb_start_x;

            // Init segments within the tile group
            enc_dec_segments_init(
                parentpicture_control_set_ptr->child_pcs->enc_dec_segment_ctrl[tile_group_idx],
                enc_dec_seg_col_cnt,
                enc_dec_seg_row_cnt,
                tg_info_ptr->tile_group_width_in_sb,
                tg_info_ptr->tile_group_height_in_sb);
            // Enable tile parallelism in Entropy Coding stage
            for (uint16_t s = top_left_tile_row_idx; s < bottom_right_tile_row_idx; s++) {
                for (uint16_t d = top_left_tile_col_idx; d < bottom_right_tile_col_idx; d++) {
                    uint16_t tileIdx = s * tile_cols + d;
                    parentpicture_control_set_ptr->child_pcs->entropy_coding_info[tileIdx]
                        ->entropy_coding_tile_done = FALSE;
                }
            }
            parentpicture_control_set_ptr->child_pcs->entropy_coding_pic_reset_flag = TRUE;
        }
    }
}

void superres_setup_child_pcs(SequenceControlSet      *entry_scs_ptr,
                              PictureParentControlSet *entry_pcs_ptr) {
    PictureControlSet *child_pcs_ptr = entry_pcs_ptr->child_pcs;

    uint8_t pic_width_in_sb;
    uint8_t picture_height_in_sb;

    child_pcs_ptr->b64_total_count = entry_pcs_ptr->b64_total_count;

    pic_width_in_sb      = (uint8_t)((entry_pcs_ptr->aligned_width + entry_scs_ptr->sb_size - 1) /
                                entry_scs_ptr->sb_size);
    picture_height_in_sb = (uint8_t)((entry_pcs_ptr->aligned_height + entry_scs_ptr->sb_size - 1) /
                                     entry_scs_ptr->sb_size);

    child_pcs_ptr->sb_total_count = pic_width_in_sb * picture_height_in_sb;

    //if (entry_pcs_ptr->frame_superres_enabled)
    {
        // Modify sb_prt_array in child pcs
        uint16_t sb_index;
        uint16_t sb_origin_x = 0;
        uint16_t sb_origin_y = 0;
        for (sb_index = 0; sb_index < child_pcs_ptr->sb_total_count; ++sb_index) {
            largest_coding_unit_dctor(child_pcs_ptr->sb_ptr_array[sb_index]);
            largest_coding_unit_ctor(child_pcs_ptr->sb_ptr_array[sb_index],
                                     (uint8_t)entry_scs_ptr->sb_size,
                                     (uint16_t)(sb_origin_x * entry_scs_ptr->sb_size),
                                     (uint16_t)(sb_origin_y * entry_scs_ptr->sb_size),
                                     (uint16_t)sb_index,
                                     child_pcs_ptr->enc_mode,
                                     entry_scs_ptr->max_block_cnt,
                                     child_pcs_ptr);
            // Increment the Order in coding order (Raster Scan Order)
            sb_origin_y = (sb_origin_x == pic_width_in_sb - 1) ? sb_origin_y + 1 : sb_origin_y;
            sb_origin_x = (sb_origin_x == pic_width_in_sb - 1) ? 0 : sb_origin_x + 1;
        }
    }

    // Update pcs_ptr->mi_stride
    child_pcs_ptr->mi_stride = pic_width_in_sb * (entry_scs_ptr->sb_size >> MI_SIZE_LOG2);

    // init segment since picture scaled
    init_enc_dec_segement(entry_pcs_ptr);

    //Tile Loop
    int              sb_size_log2 = entry_scs_ptr->seq_header.sb_size_log2;
    Av1Common *const cm           = entry_pcs_ptr->av1_cm;
    const int        tile_cols    = entry_pcs_ptr->av1_cm->tiles_info.tile_cols;
    const int        tile_rows    = entry_pcs_ptr->av1_cm->tiles_info.tile_rows;
    uint32_t         x_sb_index, y_sb_index;
    uint16_t         tile_row, tile_col;
    TileInfo         tile_info;
    for (tile_row = 0; tile_row < tile_rows; tile_row++) {
        svt_av1_tile_set_row(&tile_info, &cm->tiles_info, cm->mi_rows, tile_row);

        for (tile_col = 0; tile_col < tile_cols; tile_col++) {
            svt_av1_tile_set_col(&tile_info, &cm->tiles_info, cm->mi_cols, tile_col);
            tile_info.tile_rs_index = tile_col + tile_row * tile_cols;

            for ((y_sb_index = cm->tiles_info.tile_row_start_mi[tile_row] >> sb_size_log2);
                 (y_sb_index <
                  ((uint32_t)cm->tiles_info.tile_row_start_mi[tile_row + 1] >> sb_size_log2));
                 y_sb_index++) {
                for ((x_sb_index = cm->tiles_info.tile_col_start_mi[tile_col] >> sb_size_log2);
                     (x_sb_index <
                      ((uint32_t)cm->tiles_info.tile_col_start_mi[tile_col + 1] >> sb_size_log2));
                     x_sb_index++) {
                    int sb_index = (uint16_t)(x_sb_index + y_sb_index * pic_width_in_sb);
                    child_pcs_ptr->sb_ptr_array[sb_index]->tile_info = tile_info;
                }
            }
        }
    }

    child_pcs_ptr->enc_dec_coded_sb_count = 0;
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
void *picture_manager_kernel(void *input_ptr) {
    EbThreadContext       *thread_context_ptr = (EbThreadContext *)input_ptr;
    PictureManagerContext *context_ptr        = (PictureManagerContext *)thread_context_ptr->priv;

    EbObjectWrapper *enc_dec_wrapper_ptr;
    EncDecSet       *enc_dec_ptr;

    EbObjectWrapper         *child_pcs_wrapper_ptr;
    PictureControlSet       *child_pcs_ptr;
    PictureParentControlSet *pcs_ptr;
    SequenceControlSet      *scs_ptr;
    EncodeContext           *encode_context_ptr;

    EbObjectWrapper     *input_picture_demux_wrapper_ptr;
    PictureDemuxResults *input_picture_demux_ptr;

    Bool availability_flag;

    PredictionStructureEntry *pred_position_ptr;
    InputQueueEntry          *input_entry_ptr;
    uint32_t                  input_queue_index;
    ReferenceQueueEntry      *reference_entry_ptr;
    uint32_t                  reference_queue_index;
    uint64_t                  ref_poc;
    uint32_t                  dep_idx;
    PictureParentControlSet  *entry_pcs_ptr;
    SequenceControlSet       *entry_scs_ptr;

    // Initialization
    uint16_t pic_width_in_sb;
    uint16_t picture_height_in_sb;
    uint64_t decode_order = 0;
    // Debug
    uint32_t loop_count = 0;

    for (;;) {
        // Get Input Full Object
        EB_GET_FULL_OBJECT(context_ptr->picture_input_fifo_ptr, &input_picture_demux_wrapper_ptr);

        input_picture_demux_ptr = (PictureDemuxResults *)
                                      input_picture_demux_wrapper_ptr->object_ptr;

        // *Note - This should be overhauled and/or replaced when we
        //   need hierarchical support.
        loop_count++;

        switch (input_picture_demux_ptr->picture_type) {
        case EB_PIC_SUPERRES_INPUT: {
            pcs_ptr = (PictureParentControlSet *)
                          input_picture_demux_ptr->pcs_wrapper_ptr->object_ptr;
            scs_ptr = pcs_ptr->scs_ptr;

            assert(scs_ptr->static_config.superres_mode == SUPERRES_QTHRESH ||
                   scs_ptr->static_config.superres_mode == SUPERRES_AUTO);

            // setup child pcs to reflect superres config. E.g. sb count, sb orig, tile info, etc.
            superres_setup_child_pcs(scs_ptr, pcs_ptr);

            EbObjectWrapper *out_results_wrapper_ptr;
            // Get Empty Results Object
            svt_get_empty_object(context_ptr->picture_manager_output_fifo_ptr,
                                 &out_results_wrapper_ptr);
            RateControlTasks *rate_control_tasks_ptr = (RateControlTasks *)
                                                           out_results_wrapper_ptr->object_ptr;
            rate_control_tasks_ptr->pcs_wrapper_ptr = pcs_ptr->child_pcs->c_pcs_wrapper_ptr;
            rate_control_tasks_ptr->task_type       = RC_INPUT_SUPERRES_RECODE;
            // Post the Full Results Object
            svt_post_full_object(out_results_wrapper_ptr);

            pcs_ptr            = (PictureParentControlSet *)NULL;
            encode_context_ptr = (EncodeContext *)NULL;
            break;
        }
        case EB_PIC_INPUT:

            pcs_ptr = (PictureParentControlSet *)
                          input_picture_demux_ptr->pcs_wrapper_ptr->object_ptr;
            scs_ptr            = pcs_ptr->scs_ptr;
            encode_context_ptr = scs_ptr->encode_context_ptr;

            //SVT_LOG("\nPicture Manager Process @ %d \n ", pcs_ptr->picture_number);
            pred_position_ptr =
                pcs_ptr->pred_struct_ptr->pred_struct_entry_ptr_array[pcs_ptr->pred_struct_index];

            // overlay dep is counted by alt-ref, does not copy to cleaning list
            if (!pcs_ptr->is_overlay)
                copy_dep_cnt_cleaning_list(encode_context_ptr, pcs_ptr);

            clean_pictures_in_ref_queue(encode_context_ptr);

            CHECK_REPORT_ERROR(
                (((encode_context_ptr->input_picture_queue_head_index !=
                   encode_context_ptr->input_picture_queue_tail_index) ||
                  (encode_context_ptr
                       ->input_picture_queue[encode_context_ptr->input_picture_queue_head_index]
                       ->input_object_ptr == NULL))),
                encode_context_ptr->app_callback_ptr,
                EB_ENC_PM_ERROR4);

            // Place Picture in input queue
            input_entry_ptr =
                encode_context_ptr
                    ->input_picture_queue[encode_context_ptr->input_picture_queue_tail_index];
            input_entry_ptr->input_object_ptr = input_picture_demux_ptr->pcs_wrapper_ptr;
            encode_context_ptr->input_picture_queue_tail_index =
                (encode_context_ptr->input_picture_queue_tail_index == INPUT_QUEUE_MAX_DEPTH - 1)
                ? 0
                : encode_context_ptr->input_picture_queue_tail_index + 1;

            // Copy the reference lists into the inputEntry and
            // set the Reference Counts Based on Temporal Layer and how many frames are active
            input_entry_ptr->list0_ptr = &pred_position_ptr->ref_list0;
            input_entry_ptr->list1_ptr = &pred_position_ptr->ref_list1;
            if (!pcs_ptr->is_overlay) {
                // Check if the ReferencePictureQueue is full.
                CHECK_REPORT_ERROR(
                    (((encode_context_ptr->reference_picture_queue_head_index !=
                       encode_context_ptr->reference_picture_queue_tail_index) ||
                      (encode_context_ptr
                           ->reference_picture_queue[encode_context_ptr
                                                         ->reference_picture_queue_head_index]
                           ->reference_object_ptr == NULL))),
                    encode_context_ptr->app_callback_ptr,
                    EB_ENC_PM_ERROR5);

                // Create Reference Queue Entry even if picture will not be referenced
                reference_entry_ptr = encode_context_ptr->reference_picture_queue
                                          [encode_context_ptr->reference_picture_queue_tail_index];
                reference_entry_ptr->picture_number       = pcs_ptr->picture_number;
                reference_entry_ptr->reference_object_ptr = (EbObjectWrapper *)NULL;
                reference_entry_ptr->ref_wraper =
                    pcs_ptr
                        ->reference_picture_wrapper_ptr; //at this time only the src data is valid
                reference_entry_ptr->release_enable            = TRUE;
                reference_entry_ptr->reference_available       = FALSE;
                reference_entry_ptr->slice_type                = pcs_ptr->slice_type;
                reference_entry_ptr->temporal_layer_index      = pcs_ptr->temporal_layer_index;
                reference_entry_ptr->frame_context_updated     = FALSE;
                reference_entry_ptr->is_alt_ref                = pcs_ptr->is_alt_ref;
                reference_entry_ptr->feedback_arrived          = FALSE;
                reference_entry_ptr->is_used_as_reference_flag = pcs_ptr->is_used_as_reference_flag;
                encode_context_ptr->reference_picture_queue_tail_index =
                    (encode_context_ptr->reference_picture_queue_tail_index ==
                     REFERENCE_QUEUE_MAX_DEPTH - 1)
                    ? 0
                    : encode_context_ptr->reference_picture_queue_tail_index + 1;

                // Copy the Dependent Lists
                // *Note - we are removing any leading picture dependencies for now
                reference_entry_ptr->list0.list_count = 0;
                for (dep_idx = 0; dep_idx < pred_position_ptr->dep_list0.list_count; ++dep_idx) {
                    if (pred_position_ptr->dep_list0.list[dep_idx] >= 0)
                        reference_entry_ptr->list0.list[reference_entry_ptr->list0.list_count++] =
                            pred_position_ptr->dep_list0.list[dep_idx];
                }

                reference_entry_ptr->list1.list_count = pred_position_ptr->dep_list1.list_count;
                for (dep_idx = 0; dep_idx < pred_position_ptr->dep_list1.list_count; ++dep_idx)
                    reference_entry_ptr->list1.list[dep_idx] =
                        pred_position_ptr->dep_list1.list[dep_idx];
                reference_entry_ptr->dep_list0_count = (pcs_ptr->is_alt_ref)
                    ? reference_entry_ptr->list0.list_count + 1
                    : reference_entry_ptr->list0.list_count;

                reference_entry_ptr->dep_list1_count = reference_entry_ptr->list1.list_count;
                reference_entry_ptr->dependent_count = reference_entry_ptr->dep_list0_count +
                    reference_entry_ptr->dep_list1_count;
                //remove dependency to alredy knowmn broken links from PD
                reference_entry_ptr->dependent_count =
                    (int32_t)reference_entry_ptr->dependent_count + pcs_ptr->self_updated_links;

                CHECK_REPORT_ERROR(
                    (pcs_ptr->pred_struct_ptr->pred_struct_period * REF_LIST_MAX_DEPTH <
                     MAX_ELAPSED_IDR_COUNT),
                    encode_context_ptr->app_callback_ptr,
                    EB_ENC_PM_ERROR6);
            }
            break;

        case EB_PIC_REFERENCE:

            scs_ptr            = input_picture_demux_ptr->scs_ptr;
            encode_context_ptr = scs_ptr->encode_context_ptr;
            clean_pictures_in_ref_queue(scs_ptr->encode_context_ptr);
            ((EbReferenceObject *)
                 input_picture_demux_ptr->reference_picture_wrapper_ptr->object_ptr)
                ->ds_pics.picture_number = input_picture_demux_ptr->picture_number;
            // Check if Reference Queue is full
            CHECK_REPORT_ERROR((encode_context_ptr->reference_picture_queue_head_index !=
                                encode_context_ptr->reference_picture_queue_tail_index),
                               encode_context_ptr->app_callback_ptr,
                               EB_ENC_PM_ERROR7);

            reference_queue_index = encode_context_ptr->reference_picture_queue_head_index;
            // Find the Reference in the Reference Queue
            do {
                reference_entry_ptr =
                    encode_context_ptr->reference_picture_queue[reference_queue_index];

                if (reference_entry_ptr->picture_number ==
                    input_picture_demux_ptr->picture_number) {
                    // Assign the reference object if there is a match
                    reference_entry_ptr->reference_object_ptr =
                        input_picture_demux_ptr->reference_picture_wrapper_ptr;

                    // Set the reference availability
                    reference_entry_ptr->reference_available = TRUE;
                }

                // Increment the reference_queue_index Iterator
                reference_queue_index = (reference_queue_index == REFERENCE_QUEUE_MAX_DEPTH - 1)
                    ? 0
                    : reference_queue_index + 1;
            } while (
                (reference_queue_index != encode_context_ptr->reference_picture_queue_tail_index) &&
                (reference_entry_ptr->picture_number != input_picture_demux_ptr->picture_number));

            CHECK_REPORT_ERROR(
                (reference_entry_ptr->picture_number == input_picture_demux_ptr->picture_number),
                encode_context_ptr->app_callback_ptr,
                EB_ENC_PM_ERROR8);
            break;
        case EB_PIC_FEEDBACK:
            scs_ptr            = input_picture_demux_ptr->scs_ptr;
            encode_context_ptr = scs_ptr->encode_context_ptr;

            clean_pictures_in_ref_queue(scs_ptr->encode_context_ptr);
            reference_queue_index = encode_context_ptr->reference_picture_queue_head_index;
            // Find the Reference in the Reference Queue
            do {
                reference_entry_ptr =
                    encode_context_ptr->reference_picture_queue[reference_queue_index];
                if (reference_entry_ptr->picture_number ==
                    input_picture_demux_ptr->picture_number) {
                    // Set the feedback arrived
                    reference_entry_ptr->feedback_arrived      = TRUE;
                    reference_entry_ptr->frame_context_updated = TRUE;
                }
                // Increment the reference_queue_index Iterator
                reference_queue_index = (reference_queue_index == REFERENCE_QUEUE_MAX_DEPTH - 1)
                    ? 0
                    : reference_queue_index + 1;

            } while (
                (reference_queue_index != encode_context_ptr->reference_picture_queue_tail_index) &&
                (reference_entry_ptr->picture_number != input_picture_demux_ptr->picture_number));
            // Update the last decode order
            if (input_picture_demux_ptr->decode_order == decode_order)
                decode_order++;
            break;
        default:
            scs_ptr            = input_picture_demux_ptr->scs_ptr;
            encode_context_ptr = scs_ptr->encode_context_ptr;

            CHECK_REPORT_ERROR_NC(encode_context_ptr->app_callback_ptr, EB_ENC_PM_ERROR9);

            pcs_ptr            = (PictureParentControlSet *)NULL;
            encode_context_ptr = (EncodeContext *)NULL;

            break;
        }

        // ***********************************
        //  Common Code
        // *************************************

        // Walk the input queue and start all ready pictures.  Mark entry as null after started.  Increment the head as you go.
        if (encode_context_ptr != (EncodeContext *)NULL) {
            input_queue_index = encode_context_ptr->input_picture_queue_head_index;
            while (input_queue_index != encode_context_ptr->input_picture_queue_tail_index) {
                input_entry_ptr = encode_context_ptr->input_picture_queue[input_queue_index];

                if (input_entry_ptr->input_object_ptr != NULL) {
                    entry_pcs_ptr = (PictureParentControlSet *)
                                        input_entry_ptr->input_object_ptr->object_ptr;
                    entry_scs_ptr     = entry_pcs_ptr->scs_ptr;
                    availability_flag = TRUE;
                    if (entry_pcs_ptr->decode_order != decode_order && (scs_ptr->enable_dec_order))
                        availability_flag = FALSE;

                    //pic mgr starts pictures in dec order (no need to wait for feedback)
                    if (entry_scs_ptr->enable_pic_mgr_dec_order)
                        if (entry_pcs_ptr->picture_number > 0 &&
                            entry_pcs_ptr->decode_order != context_ptr->pmgr_dec_order + 1)
                            availability_flag = FALSE;
                    // Check RefList0 Availability
                    for (uint8_t ref_idx = 0; ref_idx < entry_pcs_ptr->ref_list0_count; ++ref_idx) {
                        //if (entry_pcs_ptr->ref_list0_count)  // NM: to double check.
                        {
                            if (entry_pcs_ptr->is_overlay)
                                // hardcode the reference for the overlay frame
                                ref_poc = entry_pcs_ptr->picture_number;
                            else
                                ref_poc = POC_CIRCULAR_ADD(
                                    (int64_t)entry_pcs_ptr->picture_number,
                                    -input_entry_ptr->list0_ptr->reference_list[ref_idx]);

                            reference_entry_ptr = search_ref_in_ref_queue(encode_context_ptr,
                                                                          ref_poc);

                            if (reference_entry_ptr != NULL) {
                                availability_flag = (availability_flag == FALSE)
                                    ? FALSE
                                    : // Don't update if already False
                                    (scs_ptr->static_config.rate_control_mode &&
                                     entry_pcs_ptr->slice_type != I_SLICE &&
                                     entry_pcs_ptr->temporal_layer_index == 0 &&
                                     !reference_entry_ptr->feedback_arrived &&
                                     !encode_context_ptr->terminating_sequence_flag_received)
                                    ? FALSE
                                    : (entry_pcs_ptr->frame_end_cdf_update_mode &&
                                       !reference_entry_ptr->frame_context_updated)
                                    ? FALSE
                                    : (reference_entry_ptr->reference_available)
                                    ? TRUE
                                    : // The Reference has been completed
                                    FALSE; // The Reference has not been completed
                            } else {
                                availability_flag = FALSE;
                            }
                        }
                    }
                    // Check RefList1 Availability
                    if (entry_pcs_ptr->slice_type == B_SLICE) {
                        for (uint8_t ref_idx = 0; ref_idx < entry_pcs_ptr->ref_list1_count;
                             ++ref_idx) {
                            // if (entry_pcs_ptr->ref_list1_count) // NM: To double check
                            {
                                // If Reference is valid (non-zero), update the availability
                                if (input_entry_ptr->list1_ptr->reference_list[ref_idx] !=
                                    (int32_t)INVALID_POC) {
                                    ref_poc = POC_CIRCULAR_ADD(
                                        (int64_t)entry_pcs_ptr->picture_number,
                                        -input_entry_ptr->list1_ptr->reference_list[ref_idx]);

                                    reference_entry_ptr = search_ref_in_ref_queue(
                                        encode_context_ptr, ref_poc);

                                    if (reference_entry_ptr != NULL) {
                                        availability_flag = (availability_flag == FALSE)
                                            ? FALSE
                                            : // Don't update if already False
                                            (scs_ptr->static_config.rate_control_mode &&
                                             entry_pcs_ptr->slice_type != I_SLICE &&
                                             entry_pcs_ptr->temporal_layer_index == 0 &&
                                             !reference_entry_ptr->feedback_arrived &&
                                             !encode_context_ptr
                                                  ->terminating_sequence_flag_received)
                                            ? FALSE
                                            : (entry_pcs_ptr->frame_end_cdf_update_mode &&
                                               !reference_entry_ptr->frame_context_updated)
                                            ? FALSE
                                            : (reference_entry_ptr->reference_available)
                                            ? TRUE
                                            : // The Reference has been completed
                                            FALSE; // The Reference has not been completed

                                    } else {
                                        availability_flag = FALSE;
                                    }
                                }
                            }
                        }
                    }

                    if (availability_flag == TRUE) {
                        // Get New  Empty recon-coef from recon-coef  Pool
                        svt_get_empty_object(context_ptr->recon_coef_fifo_ptr,
                                             &enc_dec_wrapper_ptr);
                        // Child PCS is released by Packetization
                        svt_object_inc_live_count(enc_dec_wrapper_ptr, 1);
                        enc_dec_ptr = (EncDecSet *)enc_dec_wrapper_ptr->object_ptr;
                        enc_dec_ptr->enc_dec_wrapper_ptr = enc_dec_wrapper_ptr;

                        //1.Link The Child PCS to its Parent
                        enc_dec_ptr->picture_parent_control_set_wrapper_ptr =
                            input_entry_ptr->input_object_ptr;
                        enc_dec_ptr->parent_pcs_ptr = entry_pcs_ptr;

                        enc_dec_ptr->parent_pcs_ptr->enc_dec_ptr = enc_dec_ptr;
                        // Get New  Empty Child PCS from PCS Pool
                        svt_get_empty_object(context_ptr->picture_control_set_fifo_ptr,
                                             &child_pcs_wrapper_ptr);

                        // Child PCS is released by Packetization
                        svt_object_inc_live_count(child_pcs_wrapper_ptr, 1);

                        child_pcs_ptr = (PictureControlSet *)child_pcs_wrapper_ptr->object_ptr;

                        child_pcs_ptr->c_pcs_wrapper_ptr = child_pcs_wrapper_ptr;

                        //1.Link The Child PCS to its Parent
                        child_pcs_ptr->picture_parent_control_set_wrapper_ptr =
                            input_entry_ptr->input_object_ptr;
                        child_pcs_ptr->parent_pcs_ptr = entry_pcs_ptr;

                        child_pcs_ptr->parent_pcs_ptr->child_pcs = child_pcs_ptr;
                        //1b Link The Child PCS to av1_cm to be used by Restoration
                        child_pcs_ptr->parent_pcs_ptr->av1_cm->child_pcs = child_pcs_ptr;

                        //2. Have some common information between  ChildPCS and ParentPCS.
                        child_pcs_ptr->scs_ptr              = entry_pcs_ptr->scs_ptr;
                        child_pcs_ptr->picture_qp           = entry_pcs_ptr->picture_qp;
                        child_pcs_ptr->picture_number       = entry_pcs_ptr->picture_number;
                        child_pcs_ptr->slice_type           = entry_pcs_ptr->slice_type;
                        child_pcs_ptr->temporal_layer_index = entry_pcs_ptr->temporal_layer_index;

                        child_pcs_ptr->parent_pcs_ptr->total_num_bits = 0;
                        child_pcs_ptr->parent_pcs_ptr->picture_qp     = entry_pcs_ptr->picture_qp;
                        child_pcs_ptr->enc_mode                       = entry_pcs_ptr->enc_mode;
                        child_pcs_ptr->b64_total_count = entry_pcs_ptr->b64_total_count;

                        child_pcs_ptr->enc_dec_coded_sb_count = 0;
                        child_pcs_ptr->parent_pcs_ptr->av1_cm->rst_tmpbuf =
                            child_pcs_ptr->rst_tmpbuf;

                        child_pcs_ptr->hbd_mode_decision = entry_pcs_ptr->hbd_mode_decision;
                        context_ptr->pmgr_dec_order = child_pcs_ptr->parent_pcs_ptr->decode_order;
                        //3.make all  init for ChildPCS

                        pic_width_in_sb = (entry_pcs_ptr->aligned_width + entry_scs_ptr->sb_size -
                                           1) /
                            entry_scs_ptr->sb_size;
                        picture_height_in_sb = (entry_pcs_ptr->aligned_height +
                                                entry_scs_ptr->sb_size - 1) /
                            entry_scs_ptr->sb_size;

                        init_enc_dec_segement(entry_pcs_ptr);

                        int sb_size_log2 = entry_scs_ptr->seq_header.sb_size_log2;
                        struct PictureParentControlSet *ppcs_ptr = child_pcs_ptr->parent_pcs_ptr;
                        const int        tile_cols = ppcs_ptr->av1_cm->tiles_info.tile_cols;
                        const int        tile_rows = ppcs_ptr->av1_cm->tiles_info.tile_rows;
                        Av1Common *const cm        = ppcs_ptr->av1_cm;
                        uint16_t         tile_row, tile_col;
                        uint32_t         x_sb_index, y_sb_index;
                        TileInfo         tile_info;

                        child_pcs_ptr->sb_total_count = pic_width_in_sb * picture_height_in_sb;

                        // force re-ctor sb_ptr since child_pcs_ptrs are reused, and sb_ptr could be altered by superres tool when coding previous pictures
                        if (scs_ptr->static_config.superres_mode > SUPERRES_NONE ||
                            scs_ptr->static_config.resize_mode > RESIZE_NONE) {
                            // Modify sb_prt_array in child pcs
                            uint16_t sb_index;
                            uint16_t sb_origin_x = 0;
                            uint16_t sb_origin_y = 0;
                            for (sb_index = 0; sb_index < child_pcs_ptr->sb_total_count;
                                 ++sb_index) {
                                largest_coding_unit_dctor(child_pcs_ptr->sb_ptr_array[sb_index]);
                                largest_coding_unit_ctor(child_pcs_ptr->sb_ptr_array[sb_index],
                                                         (uint8_t)scs_ptr->sb_size,
                                                         (uint16_t)(sb_origin_x * scs_ptr->sb_size),
                                                         (uint16_t)(sb_origin_y * scs_ptr->sb_size),
                                                         (uint16_t)sb_index,
                                                         child_pcs_ptr->enc_mode,
                                                         scs_ptr->max_block_cnt,
                                                         child_pcs_ptr);
                                // Increment the Order in coding order (Raster Scan Order)
                                sb_origin_y = (sb_origin_x == pic_width_in_sb - 1) ? sb_origin_y + 1
                                                                                   : sb_origin_y;
                                sb_origin_x = (sb_origin_x == pic_width_in_sb - 1)
                                    ? 0
                                    : sb_origin_x + 1;
                            }

                            // reset input_frame16bit to align with enhanced_picture_ptr
                            if (scs_ptr->static_config.encoder_bit_depth > EB_EIGHT_BIT) {
                                copy_buffer_info(entry_pcs_ptr->enhanced_picture_ptr,
                                                 child_pcs_ptr->input_frame16bit);
                            }
                        }

                        // Update pcs_ptr->mi_stride
                        child_pcs_ptr->mi_stride = pic_width_in_sb *
                            (scs_ptr->sb_size >> MI_SIZE_LOG2);

                        // copy buffer info from the downsampled picture to the input frame 16 bit buffer
                        if ((entry_pcs_ptr->frame_superres_enabled ||
                             entry_pcs_ptr->frame_resize_enabled) &&
                            scs_ptr->static_config.encoder_bit_depth > EB_EIGHT_BIT) {
                            copy_buffer_info(entry_pcs_ptr->enhanced_downscaled_picture_ptr,
                                             child_pcs_ptr->input_frame16bit);
                        }

                        //                        child_pcs_ptr->parent_pcs_ptr->av1_cm->pcs_ptr = child_pcs_ptr;
                        // Palette
                        rtime_alloc_palette_tokens(scs_ptr, child_pcs_ptr);
                        TOKENEXTRA  *pre_tok  = child_pcs_ptr->tile_tok[0][0];
                        unsigned int tile_tok = 0;
                        //Tile Loop
                        for (tile_row = 0; tile_row < tile_rows; tile_row++) {
                            svt_av1_tile_set_row(
                                &tile_info, &cm->tiles_info, cm->mi_rows, tile_row);

                            for (tile_col = 0; tile_col < tile_cols; tile_col++) {
                                svt_av1_tile_set_col(
                                    &tile_info, &cm->tiles_info, cm->mi_cols, tile_col);
                                tile_info.tile_rs_index = tile_col + tile_row * tile_cols;

                                // Palette
                                if (child_pcs_ptr->tile_tok[0][0]) {
                                    child_pcs_ptr->tile_tok[tile_row][tile_col] = pre_tok +
                                        tile_tok;
                                    pre_tok          = child_pcs_ptr->tile_tok[tile_row][tile_col];
                                    int tile_mb_rows = (tile_info.mi_row_end -
                                                        tile_info.mi_row_start + 2) >>
                                        2;
                                    int tile_mb_cols = (tile_info.mi_col_end -
                                                        tile_info.mi_col_start + 2) >>
                                        2;
                                    tile_tok = get_token_alloc(
                                        tile_mb_rows, tile_mb_cols, (sb_size_log2 + 2), 2);
                                }
                                for ((y_sb_index = cm->tiles_info.tile_row_start_mi[tile_row] >>
                                          sb_size_log2);
                                     (y_sb_index <
                                      ((uint32_t)cm->tiles_info.tile_row_start_mi[tile_row + 1] >>
                                       sb_size_log2));
                                     y_sb_index++) {
                                    for ((x_sb_index = cm->tiles_info.tile_col_start_mi[tile_col] >>
                                              sb_size_log2);
                                         (x_sb_index < ((uint32_t)cm->tiles_info
                                                            .tile_col_start_mi[tile_col + 1] >>
                                                        sb_size_log2));
                                         x_sb_index++) {
                                        int sb_index = (uint16_t)(x_sb_index +
                                                                  y_sb_index * pic_width_in_sb);
                                        child_pcs_ptr->sb_ptr_array[sb_index]->tile_info =
                                            tile_info;
                                    }
                                }
                            }
                        }
                        cm->mi_stride = child_pcs_ptr->mi_stride;
                        // Picture edges
                        configure_picture_edges(entry_scs_ptr, child_pcs_ptr);
                        // Error resilience related
                        child_pcs_ptr->colocated_pu_ref_list = REF_LIST_0; // to be modified
                        // Reset the Reference Lists
                        EB_MEMSET(child_pcs_ptr->ref_pic_ptr_array[REF_LIST_0],
                                  0,
                                  REF_LIST_MAX_DEPTH * sizeof(EbObjectWrapper *));
                        EB_MEMSET(child_pcs_ptr->ref_pic_ptr_array[REF_LIST_1],
                                  0,
                                  REF_LIST_MAX_DEPTH * sizeof(EbObjectWrapper *));
                        EB_MEMSET(child_pcs_ptr->ref_pic_qp_array[REF_LIST_0],
                                  0,
                                  REF_LIST_MAX_DEPTH * sizeof(uint8_t));
                        EB_MEMSET(child_pcs_ptr->ref_pic_qp_array[REF_LIST_1],
                                  0,
                                  REF_LIST_MAX_DEPTH * sizeof(uint8_t));
                        EB_MEMSET(child_pcs_ptr->ref_slice_type_array[REF_LIST_0],
                                  0,
                                  REF_LIST_MAX_DEPTH * sizeof(SliceType));
                        EB_MEMSET(child_pcs_ptr->ref_slice_type_array[REF_LIST_1],
                                  0,
                                  REF_LIST_MAX_DEPTH * sizeof(SliceType));
                        EB_MEMSET(child_pcs_ptr->ref_pic_r0[REF_LIST_0],
                                  0,
                                  REF_LIST_MAX_DEPTH * sizeof(double));
                        EB_MEMSET(child_pcs_ptr->ref_pic_r0[REF_LIST_1],
                                  0,
                                  REF_LIST_MAX_DEPTH * sizeof(double));
                        int8_t max_temporal_index = -1, ref_index = 0;
                        // Configure List0
                        if ((entry_pcs_ptr->slice_type == P_SLICE) ||
                            (entry_pcs_ptr->slice_type == B_SLICE)) {
                            for (uint8_t ref_idx = 0; ref_idx < entry_pcs_ptr->ref_list0_count;
                                 ++ref_idx) {
                                if (entry_pcs_ptr->ref_list0_count) {
                                    if (entry_pcs_ptr->is_overlay)
                                        ref_poc = entry_pcs_ptr->picture_number;
                                    else
                                        ref_poc = POC_CIRCULAR_ADD(
                                            (int64_t)entry_pcs_ptr->picture_number,
                                            -input_entry_ptr->list0_ptr->reference_list[ref_idx]);

                                    reference_entry_ptr = search_ref_in_ref_queue(
                                        encode_context_ptr, ref_poc);
                                    assert(reference_entry_ptr != 0);
                                    CHECK_REPORT_ERROR((reference_entry_ptr),
                                                       encode_context_ptr->app_callback_ptr,
                                                       EB_ENC_PM_ERROR10);
                                    if (entry_pcs_ptr->frame_end_cdf_update_mode) {
                                        child_pcs_ptr->ref_frame_context[svt_get_ref_frame_type(
                                                                             REF_LIST_0, ref_idx) -
                                                                         LAST_FRAME] =
                                            ((EbReferenceObject *)reference_entry_ptr
                                                 ->reference_object_ptr->object_ptr)
                                                ->frame_context;
                                        if (max_temporal_index <
                                                (int8_t)reference_entry_ptr->temporal_layer_index &&
                                            (int8_t)reference_entry_ptr->temporal_layer_index <=
                                                child_pcs_ptr->temporal_layer_index) {
                                            max_temporal_index =
                                                (int8_t)reference_entry_ptr->temporal_layer_index;
                                            ref_index = svt_get_ref_frame_type(REF_LIST_0,
                                                                               ref_idx) -
                                                LAST_FRAME;
                                            for (int frame = LAST_FRAME; frame <= ALTREF_FRAME;
                                                 ++frame) {
                                                EbReferenceObject *ref =
                                                    (EbReferenceObject *)reference_entry_ptr
                                                        ->reference_object_ptr->object_ptr;

                                                child_pcs_ptr->ref_global_motion[frame] =
                                                    ref->slice_type != I_SLICE
                                                    ? ref->global_motion[frame]
                                                    : default_warp_params;
                                            }
                                        }
                                    }
                                    // Set the Reference Object
                                    child_pcs_ptr->ref_pic_ptr_array[REF_LIST_0][ref_idx] =
                                        reference_entry_ptr->reference_object_ptr;
                                    child_pcs_ptr->ref_pic_qp_array[REF_LIST_0][ref_idx] =
                                        (uint8_t)((EbReferenceObject *)reference_entry_ptr
                                                      ->reference_object_ptr->object_ptr)
                                            ->qp;
                                    child_pcs_ptr->ref_slice_type_array[REF_LIST_0][ref_idx] =
                                        ((EbReferenceObject *)
                                             reference_entry_ptr->reference_object_ptr->object_ptr)
                                            ->slice_type;
                                    child_pcs_ptr->ref_pic_r0[REF_LIST_0][ref_idx] =
                                        ((EbReferenceObject *)
                                             reference_entry_ptr->reference_object_ptr->object_ptr)
                                            ->r0;
                                    // Increment the Reference's liveCount by the number of tiles in the input picture
                                    svt_object_inc_live_count(
                                        reference_entry_ptr->reference_object_ptr, 1);

                                    // Decrement the Reference's dependent_count Count
                                    --reference_entry_ptr->dependent_count;

#if DEBUG_SFRAME
                                    if (scs_ptr->static_config.pass != ENC_FIRST_PASS) {
                                        fprintf(stderr,
                                                "\nframe %d, layer %d, ref-list-0 count %u, ref "
                                                "frame %d, ref frame remain dep count %d\n",
                                                (int)child_pcs_ptr->picture_number,
                                                entry_pcs_ptr->temporal_layer_index,
                                                entry_pcs_ptr->ref_list0_count,
                                                (int)child_pcs_ptr->parent_pcs_ptr
                                                    ->ref_pic_poc_array[REF_LIST_0][ref_idx],
                                                reference_entry_ptr->dependent_count);
                                    }
#endif

                                    CHECK_REPORT_ERROR(
                                        (reference_entry_ptr->dependent_count != ~0u),
                                        encode_context_ptr->app_callback_ptr,
                                        EB_ENC_PM_ERROR1);
                                }
                            }
                            //fill the non used spots to be used in MFMV.
                            for (uint8_t ref_idx = entry_pcs_ptr->ref_list0_count; ref_idx < 4;
                                 ++ref_idx)
                                child_pcs_ptr->ref_pic_ptr_array[REF_LIST_0][ref_idx] =
                                    child_pcs_ptr->ref_pic_ptr_array[REF_LIST_0][0];

                            if (entry_pcs_ptr->ref_list1_count == 0) {
                                for (uint8_t ref_idx = entry_pcs_ptr->ref_list1_count; ref_idx < 3;
                                     ++ref_idx)
                                    child_pcs_ptr->ref_pic_ptr_array[REF_LIST_1][ref_idx] =
                                        child_pcs_ptr->ref_pic_ptr_array[REF_LIST_0][0];
                            }
                        }

                        // Configure List1
                        if (entry_pcs_ptr->slice_type == B_SLICE) {
                            for (uint8_t ref_idx = 0; ref_idx < entry_pcs_ptr->ref_list1_count;
                                 ++ref_idx) {
                                if (entry_pcs_ptr->ref_list1_count) {
                                    ref_poc = POC_CIRCULAR_ADD(
                                        (int64_t)entry_pcs_ptr->picture_number,
                                        -input_entry_ptr->list1_ptr->reference_list[ref_idx]);

                                    reference_entry_ptr = search_ref_in_ref_queue(
                                        encode_context_ptr, ref_poc);

                                    assert(reference_entry_ptr != 0);
                                    CHECK_REPORT_ERROR((reference_entry_ptr),
                                                       encode_context_ptr->app_callback_ptr,
                                                       EB_ENC_PM_ERROR10);
                                    if (entry_pcs_ptr->frame_end_cdf_update_mode) {
                                        child_pcs_ptr->ref_frame_context[svt_get_ref_frame_type(
                                                                             REF_LIST_1, ref_idx) -
                                                                         LAST_FRAME] =
                                            ((EbReferenceObject *)reference_entry_ptr
                                                 ->reference_object_ptr->object_ptr)
                                                ->frame_context;
                                        if (max_temporal_index <
                                                (int8_t)reference_entry_ptr->temporal_layer_index &&
                                            reference_entry_ptr->slice_type != I_SLICE &&
                                            (int8_t)reference_entry_ptr->temporal_layer_index <=
                                                child_pcs_ptr->temporal_layer_index) {
                                            max_temporal_index =
                                                (int8_t)reference_entry_ptr->temporal_layer_index;
                                            ref_index = svt_get_ref_frame_type(REF_LIST_1,
                                                                               ref_idx) -
                                                LAST_FRAME;
                                            for (int frame = LAST_FRAME; frame <= ALTREF_FRAME;
                                                 ++frame) {
                                                EbReferenceObject *ref =
                                                    (EbReferenceObject *)reference_entry_ptr
                                                        ->reference_object_ptr->object_ptr;

                                                child_pcs_ptr->ref_global_motion[frame] =
                                                    ref->slice_type != I_SLICE
                                                    ? ref->global_motion[frame]
                                                    : default_warp_params;
                                            }
                                        }
                                    }
                                    // Set the Reference Object
                                    child_pcs_ptr->ref_pic_ptr_array[REF_LIST_1][ref_idx] =
                                        reference_entry_ptr->reference_object_ptr;
                                    child_pcs_ptr->ref_pic_qp_array[REF_LIST_1][ref_idx] =
                                        (uint8_t)((EbReferenceObject *)reference_entry_ptr
                                                      ->reference_object_ptr->object_ptr)
                                            ->qp;
                                    child_pcs_ptr->ref_slice_type_array[REF_LIST_1][ref_idx] =
                                        ((EbReferenceObject *)
                                             reference_entry_ptr->reference_object_ptr->object_ptr)
                                            ->slice_type;
                                    child_pcs_ptr->ref_pic_r0[REF_LIST_1][ref_idx] =
                                        ((EbReferenceObject *)
                                             reference_entry_ptr->reference_object_ptr->object_ptr)
                                            ->r0;
                                    // Increment the Reference's liveCount by the number of tiles in the input picture
                                    svt_object_inc_live_count(
                                        reference_entry_ptr->reference_object_ptr, 1);

                                    // Decrement the Reference's dependent_count Count
                                    --reference_entry_ptr->dependent_count;

#if DEBUG_SFRAME
                                    if (scs_ptr->static_config.pass != ENC_FIRST_PASS) {
                                        fprintf(stderr,
                                                "\nframe %d, layer %d, ref-list-1 count %u, ref "
                                                "frame %d, ref frame remain dep count %d\n",
                                                (int)child_pcs_ptr->picture_number,
                                                entry_pcs_ptr->temporal_layer_index,
                                                entry_pcs_ptr->ref_list1_count,
                                                (int)child_pcs_ptr->parent_pcs_ptr
                                                    ->ref_pic_poc_array[REF_LIST_1][ref_idx],
                                                reference_entry_ptr->dependent_count);
                                    }
#endif

                                    CHECK_REPORT_ERROR(
                                        (reference_entry_ptr->dependent_count != ~0u),
                                        encode_context_ptr->app_callback_ptr,
                                        EB_ENC_PM_ERROR1);
                                }
                            }
                            //fill the non used spots to be used in MFMV.
                            if (entry_pcs_ptr->ref_list1_count) {
                                for (uint8_t ref_idx = entry_pcs_ptr->ref_list1_count; ref_idx < 3;
                                     ++ref_idx)
                                    child_pcs_ptr->ref_pic_ptr_array[REF_LIST_1][ref_idx] =
                                        child_pcs_ptr->ref_pic_ptr_array[REF_LIST_1][0];
                            }
                        }

                        if (entry_pcs_ptr->frame_end_cdf_update_mode) {
                            if (entry_pcs_ptr->slice_type != I_SLICE &&
                                entry_pcs_ptr->frm_hdr.frame_type != S_FRAME)
                                child_pcs_ptr->parent_pcs_ptr->frm_hdr.primary_ref_frame =
                                    ref_index;
                            else
                                child_pcs_ptr->parent_pcs_ptr->frm_hdr.primary_ref_frame =
                                    PRIMARY_REF_NONE;
                            child_pcs_ptr->parent_pcs_ptr->refresh_frame_context =
                                REFRESH_FRAME_CONTEXT_BACKWARD;

                        } else {
                            child_pcs_ptr->parent_pcs_ptr->frm_hdr.primary_ref_frame =
                                PRIMARY_REF_NONE;
                            // Tells each frame whether to forward their data to the next frames;
                            // never disabled so that the feature can be on in higher layers, while off
                            // in low layers.
                            child_pcs_ptr->parent_pcs_ptr->refresh_frame_context =
                                REFRESH_FRAME_CONTEXT_BACKWARD;
                        }
                        EbObjectWrapper *out_results_wrapper_ptr;
                        // Get Empty Results Object
                        svt_get_empty_object(context_ptr->picture_manager_output_fifo_ptr,
                                             &out_results_wrapper_ptr);
                        RateControlTasks *rate_control_tasks_ptr =
                            (RateControlTasks *)out_results_wrapper_ptr->object_ptr;
                        rate_control_tasks_ptr->pcs_wrapper_ptr = child_pcs_ptr->c_pcs_wrapper_ptr;
                        rate_control_tasks_ptr->task_type       = RC_INPUT;

                        // printf("picMgr sending:%x \n", rate_control_tasks_ptr->pcs_wrapper_ptr);

                        // Post the Full Results Object
                        svt_post_full_object(out_results_wrapper_ptr);
                        // Remove the Input Entry from the Input Queue
                        input_entry_ptr->input_object_ptr = (EbObjectWrapper *)NULL;
                    }
                }

                // Increment the head_index if the head is null
                encode_context_ptr->input_picture_queue_head_index =
                    (encode_context_ptr
                         ->input_picture_queue[encode_context_ptr->input_picture_queue_head_index]
                         ->input_object_ptr)
                    ? encode_context_ptr->input_picture_queue_head_index
                    : (encode_context_ptr->input_picture_queue_head_index ==
                       INPUT_QUEUE_MAX_DEPTH - 1)
                    ? 0
                    : encode_context_ptr->input_picture_queue_head_index + 1;

                // Increment the input_queue_index Iterator
                input_queue_index = (input_queue_index == INPUT_QUEUE_MAX_DEPTH - 1)
                    ? 0
                    : input_queue_index + 1;
            }

            // Walk the reference queue and remove entries that have been completely referenced.
            reference_queue_index = encode_context_ptr->reference_picture_queue_head_index;
            while (reference_queue_index !=
                   encode_context_ptr->reference_picture_queue_tail_index) {
                reference_entry_ptr =
                    encode_context_ptr->reference_picture_queue[reference_queue_index];

                // Remove the entry & release the reference if there are no remaining references
                if ((reference_entry_ptr->dependent_count == 0) &&
                    (reference_entry_ptr->reference_available) &&
                    (reference_entry_ptr->release_enable) &&
                    (reference_entry_ptr->reference_object_ptr)) {
                    // Release the nominal live_count value
                    svt_release_object(reference_entry_ptr->reference_object_ptr);
                    reference_entry_ptr->reference_object_ptr      = (EbObjectWrapper *)NULL;
                    reference_entry_ptr->reference_available       = FALSE;
                    reference_entry_ptr->is_used_as_reference_flag = FALSE;
                }

                // Increment the head_index if the head is empty
                encode_context_ptr->reference_picture_queue_head_index =
                    (encode_context_ptr
                         ->reference_picture_queue[encode_context_ptr
                                                       ->reference_picture_queue_head_index]
                         ->release_enable == FALSE)
                    ? encode_context_ptr->reference_picture_queue_head_index
                    : (encode_context_ptr
                               ->reference_picture_queue[encode_context_ptr
                                                             ->reference_picture_queue_head_index]
                               ->reference_available == FALSE &&
                       encode_context_ptr
                               ->reference_picture_queue[encode_context_ptr
                                                             ->reference_picture_queue_head_index]
                               ->is_used_as_reference_flag == TRUE)
                    ? encode_context_ptr->reference_picture_queue_head_index
                    : (encode_context_ptr
                           ->reference_picture_queue[encode_context_ptr
                                                         ->reference_picture_queue_head_index]
                           ->dependent_count > 0)
                    ? encode_context_ptr->reference_picture_queue_head_index
                    : (encode_context_ptr->reference_picture_queue_head_index ==
                       REFERENCE_QUEUE_MAX_DEPTH - 1)
                    ? 0
                    : encode_context_ptr->reference_picture_queue_head_index + 1;
                // Increment the reference_queue_index Iterator
                reference_queue_index = (reference_queue_index == REFERENCE_QUEUE_MAX_DEPTH - 1)
                    ? 0
                    : reference_queue_index + 1;
            }
        }

        // Release the Input Picture Demux Results
        svt_release_object(input_picture_demux_wrapper_ptr);
    }
    return NULL;
}
