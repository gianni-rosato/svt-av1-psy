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

#include "EbSequenceControlSet.h"
#include "EbUtility.h"

static void free_scale_evts(SvtAv1FrameScaleEvts *evts) {
    EB_FREE_ARRAY(evts->resize_denoms);
    EB_FREE_ARRAY(evts->resize_kf_denoms);
    EB_FREE_ARRAY(evts->start_frame_nums);
    evts->evt_num = 0;
}

static void svt_sequence_control_set_dctor(EbPtr p) {
    SequenceControlSet *obj = (SequenceControlSet *)p;
    if (!obj)
        return;
    EB_FREE_ARRAY(obj->b64_geom);
    EB_FREE_ARRAY(obj->sb_geom);
    free_scale_evts(&obj->static_config.frame_scale_evts);
}
/**************************************************************************************************
    General notes on how Sequence Control Sets (SCS) are used.

    SequenceControlSetInstance
        is the primary copy that interacts with the API in real-time.  When a
        change happens, the changeFlag is signaled so that appropriate action can
        be taken.  There is one scsInstance per stream/encode instance.  The scsInstance
        owns the encodeContext

    encodeContext
        has context type variables (i.e. non-config) that keep track of global parameters.

    SequenceControlSets
        general SCSs are controled by a system resource manager.  They are kept completely
        separate from the instances.  In general there is one active SCS at a time.  When the
        changeFlag is signaled, the old active SCS is no longer used for new input pictures.
        A fresh copy of the scsInstance is made to a new SCS, which becomes the active SCS.  The
        old SCS will eventually be released back into the SCS pool when its current pictures are
        finished encoding.

    Motivations
        The whole reason for this structure is due to the nature of the pipeline.  We have to
        take great care not to have pipeline mismanagement.  Once an object enters use in the
        pipeline, it cannot be changed on the fly or you will have pipeline coherency problems.

    **** Currently, real-time updates to the SCS are not supported.  Therefore, each instance
    has a single SCS (from the SequenceControlSetInstance) that is used for encoding the entire
    stream.  At the resource coordination kernel, a pointer to the SequenceControlSetInstance
    is saved in the PCS, that is not managed by an SRM.
 ***************************************************************************************************/
EbErrorType svt_sequence_control_set_ctor(SequenceControlSet *scs, EbPtr object_init_data_ptr) {
    UNUSED(object_init_data_ptr);
    scs->dctor = svt_sequence_control_set_dctor;

    // Allocation will happen in resource-coordination
    scs->b64_geom = NULL;

    scs->mvrate_set                   = 0;
    scs->bits_for_picture_order_count = 16;
    scs->film_grain_random_seed       = 7391;

    // Initialize certain sequence header variables here for write_sequence_header(),
    // which may be called before the first picture hits resource coordination thread
    // (e.g. when ffmpeg is used it may be called first to construct mkv/mp4 container headers).
    // Whenever possible, it is recommended to initialize all sequence header info here
    // instead of in resource coordination.
    scs->seq_header.frame_width_bits              = 16;
    scs->seq_header.frame_height_bits             = 16;
    scs->seq_header.frame_id_numbers_present_flag = 0;
    scs->seq_header.frame_id_length               = FRAME_ID_LENGTH;
    scs->seq_header.delta_frame_id_length         = DELTA_FRAME_ID_LENGTH;

    // 0 - disable dual interpolation filter
    // 1 - enable vertical and horiz filter selection
    scs->seq_header.enable_dual_filter = 0;

    // 0 - force off
    // 1 - force on
    // 2 - adaptive
    scs->seq_header.seq_force_screen_content_tools = 2;

    // 0 - Not to force. MV can be in 1/4 or 1/8
    // 1 - force to integer
    // 2 - adaptive
    scs->seq_header.seq_force_integer_mv = 2;

    scs->seq_header.order_hint_info.enable_ref_frame_mvs = 1;
    scs->seq_header.order_hint_info.enable_order_hint    = 1;
    scs->seq_header.order_hint_info.order_hint_bits      = 7;

    return EB_ErrorNone;
}
extern EbErrorType svt_aom_derive_input_resolution(EbInputResolution *input_resolution, uint32_t inputSize) {
    EbErrorType return_error = EB_ErrorNone;
    if (inputSize < INPUT_SIZE_240p_TH)
        *input_resolution = INPUT_SIZE_240p_RANGE;
    else if (inputSize < INPUT_SIZE_360p_TH)
        *input_resolution = INPUT_SIZE_360p_RANGE;
    else if (inputSize < INPUT_SIZE_480p_TH)
        *input_resolution = INPUT_SIZE_480p_RANGE;
    else if (inputSize < INPUT_SIZE_720p_TH)
        *input_resolution = INPUT_SIZE_720p_RANGE;
    else if (inputSize < INPUT_SIZE_1080p_TH)
        *input_resolution = INPUT_SIZE_1080p_RANGE;
    else if (inputSize < INPUT_SIZE_4K_TH)
        *input_resolution = INPUT_SIZE_4K_RANGE;
    else
        *input_resolution = INPUT_SIZE_8K_RANGE;

    return return_error;
}

static void svt_sequence_control_set_instance_dctor(EbPtr p) {
    EbSequenceControlSetInstance *obj = (EbSequenceControlSetInstance *)p;
    EB_DELETE(obj->enc_ctx);
    EB_DESTROY_SEMAPHORE(obj->scs->ref_buffer_available_semaphore);
    EB_DELETE(obj->scs);
}

EbErrorType svt_sequence_control_set_instance_ctor(EbSequenceControlSetInstance *object_ptr) {
    object_ptr->dctor = svt_sequence_control_set_instance_dctor;

    EB_NEW(object_ptr->enc_ctx, svt_aom_encode_context_ctor, NULL);
    EB_NEW(object_ptr->scs, svt_sequence_control_set_ctor, NULL);
    object_ptr->scs->enc_ctx = object_ptr->enc_ctx;

    return EB_ErrorNone;
}

extern EbErrorType svt_aom_b64_geom_init(SequenceControlSet *scs) {
    EbErrorType return_error = EB_ErrorNone;
    uint16_t    b64_idx;
    uint16_t    raster_scan_blk_index;
    uint8_t     b64_size = scs->b64_size;

    uint16_t picture_b64_width  = (scs->max_input_luma_width + b64_size - 1) / b64_size;
    uint16_t picture_b64_height = (scs->max_input_luma_height + b64_size - 1) / b64_size;
    //free old one;
    EB_FREE_ARRAY(scs->b64_geom);
    EB_MALLOC_ARRAY(scs->b64_geom, picture_b64_width * picture_b64_height);

    for (b64_idx = 0; b64_idx < picture_b64_width * picture_b64_height; ++b64_idx) {
        B64Geom *b64_geom          = &scs->b64_geom[b64_idx];
        b64_geom->horizontal_index = (uint8_t)(b64_idx % picture_b64_width);
        b64_geom->vertical_index   = (uint8_t)(b64_idx / picture_b64_width);
        b64_geom->org_x            = b64_geom->horizontal_index * b64_size;
        b64_geom->org_y            = b64_geom->vertical_index * b64_size;

        b64_geom->width = (uint8_t)(((scs->max_input_luma_width - b64_geom->org_x) < b64_size)
                                        ? scs->max_input_luma_width - b64_geom->org_x
                                        : b64_size);

        b64_geom->height = (uint8_t)(((scs->max_input_luma_height - b64_geom->org_y) < b64_size)
                                         ? scs->max_input_luma_height - b64_geom->org_y
                                         : b64_size);

        b64_geom->is_complete_b64 = (uint8_t)(((b64_geom->width == b64_size) && (b64_geom->height == b64_size)) ? 1
                                                                                                                : 0);

        b64_geom->is_edge_sb = (b64_geom->org_x < b64_size) || (b64_geom->org_y < b64_size) ||
                (b64_geom->org_x > scs->max_input_luma_width - b64_size) ||
                (b64_geom->org_y > scs->max_input_luma_height - b64_size)
            ? 1
            : 0;

        for (raster_scan_blk_index = RASTER_SCAN_CU_INDEX_64x64; raster_scan_blk_index <= RASTER_SCAN_CU_INDEX_8x8_63;
             raster_scan_blk_index++) {
            b64_geom->raster_scan_blk_validity[raster_scan_blk_index] =
                ((b64_geom->org_x + raster_scan_blk_x[raster_scan_blk_index] +
                      raster_scan_blk_size[raster_scan_blk_index] >
                  scs->max_input_luma_width) ||
                 (b64_geom->org_y + raster_scan_blk_y[raster_scan_blk_index] +
                      raster_scan_blk_size[raster_scan_blk_index] >
                  scs->max_input_luma_height))
                ? FALSE
                : TRUE;
        }
    }

    scs->pic_width_in_b64  = picture_b64_width;
    scs->pic_height_in_b64 = picture_b64_height;
    scs->b64_total_count   = picture_b64_width * picture_b64_height;

    return return_error;
}

EbErrorType rtime_alloc_sb_geom(SequenceControlSet *scs, uint32_t size) {
    EB_MALLOC_ARRAY(scs->sb_geom, size);
    return EB_ErrorNone;
}
EbErrorType svt_aom_sb_geom_init(SequenceControlSet *scs) {
    uint16_t sb_index;
    uint16_t md_scan_block_index;
    uint16_t picture_sb_width  = (scs->max_input_luma_width + scs->sb_size - 1) / scs->sb_size;
    uint16_t picture_sb_height = (scs->max_input_luma_height + scs->sb_size - 1) / scs->sb_size;

    EB_FREE_ARRAY(scs->sb_geom);
    rtime_alloc_sb_geom(scs, picture_sb_width * picture_sb_height);

    for (sb_index = 0; sb_index < picture_sb_width * picture_sb_height; ++sb_index) {
        scs->sb_geom[sb_index].horizontal_index = sb_index % picture_sb_width;
        scs->sb_geom[sb_index].vertical_index   = sb_index / picture_sb_width;
        scs->sb_geom[sb_index].org_x            = scs->sb_geom[sb_index].horizontal_index * scs->sb_size;
        scs->sb_geom[sb_index].org_y            = scs->sb_geom[sb_index].vertical_index * scs->sb_size;

        scs->sb_geom[sb_index].width = (uint8_t)(((scs->max_input_luma_width - scs->sb_geom[sb_index].org_x) <
                                                  scs->sb_size)
                                                     ? scs->max_input_luma_width - scs->sb_geom[sb_index].org_x
                                                     : scs->sb_size);

        scs->sb_geom[sb_index].height = (uint8_t)(((scs->max_input_luma_height - scs->sb_geom[sb_index].org_y) <
                                                   scs->sb_size)
                                                      ? scs->max_input_luma_height - scs->sb_geom[sb_index].org_y
                                                      : scs->sb_size);

        scs->sb_geom[sb_index].is_complete_sb = (uint8_t)(((scs->sb_geom[sb_index].width == scs->sb_size) &&
                                                           (scs->sb_geom[sb_index].height == scs->sb_size))
                                                              ? 1
                                                              : 0);

        uint16_t max_block_count = scs->max_block_cnt;

        for (md_scan_block_index = 0; md_scan_block_index < max_block_count; md_scan_block_index++) {
            const BlockGeom *blk_geom = get_blk_geom_mds(md_scan_block_index);
            if (scs->over_boundary_block_mode == 1) {
                const BlockGeom *sq_blk_geom = get_blk_geom_mds(blk_geom->sqi_mds);
                uint8_t has_rows = (scs->sb_geom[sb_index].org_y + sq_blk_geom->org_y + sq_blk_geom->bheight / 2 <
                                    scs->max_input_luma_height);
                uint8_t has_cols = (scs->sb_geom[sb_index].org_x + sq_blk_geom->org_x + sq_blk_geom->bwidth / 2 <
                                    scs->max_input_luma_width);

                // See AV1 spec section 5.11.4 for allowable blocks
                if (has_rows && has_cols &&
                    (scs->sb_geom[sb_index].org_y + blk_geom->org_y < scs->max_input_luma_height) &&
                    (scs->sb_geom[sb_index].org_x + blk_geom->org_x < scs->max_input_luma_width)) {
                    scs->sb_geom[sb_index].block_is_allowed[md_scan_block_index] = 1;
                } else if (blk_geom->shape == PART_H && has_cols &&
                           (scs->sb_geom[sb_index].org_y + blk_geom->org_y < scs->max_input_luma_height) &&
                           (scs->sb_geom[sb_index].org_x + blk_geom->org_x < scs->max_input_luma_width)) {
                    scs->sb_geom[sb_index].block_is_allowed[md_scan_block_index] = 1;
                } else if (blk_geom->shape == PART_V && has_rows &&
                           (scs->sb_geom[sb_index].org_y + blk_geom->org_y < scs->max_input_luma_height) &&
                           (scs->sb_geom[sb_index].org_x + blk_geom->org_x < scs->max_input_luma_width)) {
                    scs->sb_geom[sb_index].block_is_allowed[md_scan_block_index] = 1;
                } else {
                    scs->sb_geom[sb_index].block_is_allowed[md_scan_block_index] = 0;
                }

            } else {
                if (blk_geom->shape != PART_N)
                    blk_geom = get_blk_geom_mds(blk_geom->sqi_mds);

                scs->sb_geom[sb_index].block_is_allowed[md_scan_block_index] =
                    ((scs->sb_geom[sb_index].org_x + blk_geom->org_x + blk_geom->bwidth > scs->max_input_luma_width) ||
                     (scs->sb_geom[sb_index].org_y + blk_geom->org_y + blk_geom->bheight > scs->max_input_luma_height))
                    ? FALSE
                    : TRUE;
            }
        }
    }

    scs->sb_total_count = picture_sb_width * picture_sb_height;

    return 0;
}
