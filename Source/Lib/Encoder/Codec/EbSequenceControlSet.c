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

static void svt_sequence_control_set_dctor(EbPtr p) {
    SequenceControlSet *obj = (SequenceControlSet *)p;
    EB_FREE_ARRAY(obj->sb_params_array);
    EB_FREE_ARRAY(obj->sb_geom);
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
    scs->sb_params_array = NULL;

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
extern EbErrorType derive_input_resolution(EbInputResolution *input_resolution,
                                           uint32_t           inputSize) {
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
    EB_DELETE(obj->encode_context_ptr);
    EB_DELETE(obj->scs_ptr);
}

EbErrorType svt_sequence_control_set_instance_ctor(EbSequenceControlSetInstance *object_ptr) {
    object_ptr->dctor = svt_sequence_control_set_instance_dctor;

    EB_NEW(object_ptr->encode_context_ptr, encode_context_ctor, NULL);
    EB_NEW(object_ptr->scs_ptr, svt_sequence_control_set_ctor, NULL);
    object_ptr->scs_ptr->encode_context_ptr = object_ptr->encode_context_ptr;

    return EB_ErrorNone;
}

extern EbErrorType sb_params_init(SequenceControlSet *scs_ptr) {
    EbErrorType return_error = EB_ErrorNone;
    uint16_t    sb_index;
    uint16_t    raster_scan_blk_index;

    uint16_t picture_sb_width = (scs_ptr->max_input_luma_width + scs_ptr->sb_sz - 1) /
        scs_ptr->sb_sz;
    uint16_t picture_sb_height = (scs_ptr->max_input_luma_height + scs_ptr->sb_sz - 1) /
        scs_ptr->sb_sz;
    //free old one;
    EB_FREE_ARRAY(scs_ptr->sb_params_array);

    EB_MALLOC_ARRAY(scs_ptr->sb_params_array, picture_sb_width * picture_sb_height);

    for (sb_index = 0; sb_index < picture_sb_width * picture_sb_height; ++sb_index) {
        scs_ptr->sb_params_array[sb_index].horizontal_index = (uint8_t)(sb_index %
                                                                        picture_sb_width);
        scs_ptr->sb_params_array[sb_index].vertical_index = (uint8_t)(sb_index / picture_sb_width);
        scs_ptr->sb_params_array[sb_index].origin_x =
            scs_ptr->sb_params_array[sb_index].horizontal_index * scs_ptr->sb_sz;
        scs_ptr->sb_params_array[sb_index].origin_y =
            scs_ptr->sb_params_array[sb_index].vertical_index * scs_ptr->sb_sz;

        scs_ptr->sb_params_array[sb_index].width =
            (uint8_t)(((scs_ptr->max_input_luma_width -
                        scs_ptr->sb_params_array[sb_index].origin_x) < scs_ptr->sb_sz)
                          ? scs_ptr->max_input_luma_width -
                              scs_ptr->sb_params_array[sb_index].origin_x
                          : scs_ptr->sb_sz);

        scs_ptr->sb_params_array[sb_index].height =
            (uint8_t)(((scs_ptr->max_input_luma_height -
                        scs_ptr->sb_params_array[sb_index].origin_y) < scs_ptr->sb_sz)
                          ? scs_ptr->max_input_luma_height -
                              scs_ptr->sb_params_array[sb_index].origin_y
                          : scs_ptr->sb_sz);

        scs_ptr->sb_params_array[sb_index].is_complete_sb =
            (uint8_t)(((scs_ptr->sb_params_array[sb_index].width == scs_ptr->sb_sz) &&
                       (scs_ptr->sb_params_array[sb_index].height == scs_ptr->sb_sz))
                          ? 1
                          : 0);

        scs_ptr->sb_params_array[sb_index].is_edge_sb =
            (scs_ptr->sb_params_array[sb_index].origin_x < scs_ptr->sb_sz) ||
                (scs_ptr->sb_params_array[sb_index].origin_y < scs_ptr->sb_sz) ||
                (scs_ptr->sb_params_array[sb_index].origin_x >
                 scs_ptr->max_input_luma_width - scs_ptr->sb_sz) ||
                (scs_ptr->sb_params_array[sb_index].origin_y >
                 scs_ptr->max_input_luma_height - scs_ptr->sb_sz)
            ? 1
            : 0;

        for (raster_scan_blk_index = RASTER_SCAN_CU_INDEX_64x64;
             raster_scan_blk_index <= RASTER_SCAN_CU_INDEX_8x8_63;
             raster_scan_blk_index++) {
            scs_ptr->sb_params_array[sb_index].raster_scan_blk_validity[raster_scan_blk_index] =
                ((scs_ptr->sb_params_array[sb_index].origin_x +
                      raster_scan_blk_x[raster_scan_blk_index] +
                      raster_scan_blk_size[raster_scan_blk_index] >
                  scs_ptr->max_input_luma_width) ||
                 (scs_ptr->sb_params_array[sb_index].origin_y +
                      raster_scan_blk_y[raster_scan_blk_index] +
                      raster_scan_blk_size[raster_scan_blk_index] >
                  scs_ptr->max_input_luma_height))
                ? FALSE
                : TRUE;
        }
    }

    scs_ptr->pic_width_in_sb      = picture_sb_width;
    scs_ptr->picture_height_in_sb = picture_sb_height;
    scs_ptr->sb_total_count       = picture_sb_width * picture_sb_height;

    return return_error;
}

EbErrorType rtime_alloc_sb_geom(SequenceControlSet *scs_ptr, uint32_t size) {
    EB_MALLOC_ARRAY(scs_ptr->sb_geom, size);
    return EB_ErrorNone;
}
EbErrorType sb_geom_init(SequenceControlSet *scs_ptr) {
    uint16_t sb_index;
    uint16_t md_scan_block_index;
    uint16_t picture_sb_width = (scs_ptr->max_input_luma_width + scs_ptr->sb_size_pix - 1) /
        scs_ptr->sb_size_pix;
    uint16_t picture_sb_height = (scs_ptr->max_input_luma_height + scs_ptr->sb_size_pix - 1) /
        scs_ptr->sb_size_pix;

    EB_FREE_ARRAY(scs_ptr->sb_geom);
    rtime_alloc_sb_geom(scs_ptr, picture_sb_width * picture_sb_height);

    for (sb_index = 0; sb_index < picture_sb_width * picture_sb_height; ++sb_index) {
        scs_ptr->sb_geom[sb_index].horizontal_index = sb_index % picture_sb_width;
        scs_ptr->sb_geom[sb_index].vertical_index   = sb_index / picture_sb_width;
        scs_ptr->sb_geom[sb_index].origin_x         = scs_ptr->sb_geom[sb_index].horizontal_index *
            scs_ptr->sb_size_pix;
        scs_ptr->sb_geom[sb_index].origin_y = scs_ptr->sb_geom[sb_index].vertical_index *
            scs_ptr->sb_size_pix;

        scs_ptr->sb_geom[sb_index].width =
            (uint8_t)(((scs_ptr->max_input_luma_width - scs_ptr->sb_geom[sb_index].origin_x) <
                       scs_ptr->sb_size_pix)
                          ? scs_ptr->max_input_luma_width - scs_ptr->sb_geom[sb_index].origin_x
                          : scs_ptr->sb_size_pix);

        scs_ptr->sb_geom[sb_index].height =
            (uint8_t)(((scs_ptr->max_input_luma_height - scs_ptr->sb_geom[sb_index].origin_y) <
                       scs_ptr->sb_size_pix)
                          ? scs_ptr->max_input_luma_height - scs_ptr->sb_geom[sb_index].origin_y
                          : scs_ptr->sb_size_pix);

        scs_ptr->sb_geom[sb_index].is_complete_sb =
            (uint8_t)(((scs_ptr->sb_geom[sb_index].width == scs_ptr->sb_size_pix) &&
                       (scs_ptr->sb_geom[sb_index].height == scs_ptr->sb_size_pix))
                          ? 1
                          : 0);

        uint16_t max_block_count = scs_ptr->max_block_cnt;

        for (md_scan_block_index = 0; md_scan_block_index < max_block_count;
             md_scan_block_index++) {
            const BlockGeom *blk_geom = get_blk_geom_mds(md_scan_block_index);
            if (scs_ptr->over_boundary_block_mode == 1) {
                scs_ptr->sb_geom[sb_index].block_is_allowed[md_scan_block_index] =
                    ((scs_ptr->sb_geom[sb_index].origin_x + blk_geom->origin_x +
                          blk_geom->bwidth / 2 <
                      scs_ptr->max_input_luma_width) &&
                     (scs_ptr->sb_geom[sb_index].origin_y + blk_geom->origin_y +
                          blk_geom->bheight / 2 <
                      scs_ptr->max_input_luma_height))
                    ? TRUE
                    : FALSE;

                if (blk_geom->shape != PART_N)
                    blk_geom = get_blk_geom_mds(blk_geom->sqi_mds);
                scs_ptr->sb_geom[sb_index].block_is_inside_md_scan[md_scan_block_index] =
                    ((scs_ptr->sb_geom[sb_index].origin_x >= scs_ptr->max_input_luma_width) ||
                     (scs_ptr->sb_geom[sb_index].origin_y >= scs_ptr->max_input_luma_height))
                    ? FALSE
                    : TRUE;
            } else {
                if (blk_geom->shape != PART_N)
                    blk_geom = get_blk_geom_mds(blk_geom->sqi_mds);

                scs_ptr->sb_geom[sb_index].block_is_allowed[md_scan_block_index] =
                    ((scs_ptr->sb_geom[sb_index].origin_x + blk_geom->origin_x + blk_geom->bwidth >
                      scs_ptr->max_input_luma_width) ||
                     (scs_ptr->sb_geom[sb_index].origin_y + blk_geom->origin_y + blk_geom->bheight >
                      scs_ptr->max_input_luma_height))
                    ? FALSE
                    : TRUE;

                scs_ptr->sb_geom[sb_index].block_is_inside_md_scan[md_scan_block_index] =
                    ((scs_ptr->sb_geom[sb_index].origin_x + blk_geom->origin_x + blk_geom->bwidth >
                      scs_ptr->max_input_luma_width) ||
                     (scs_ptr->sb_geom[sb_index].origin_y + blk_geom->origin_y + blk_geom->bheight >
                      scs_ptr->max_input_luma_height))
                    ? FALSE
                    : TRUE;
            }
        }
    }

    scs_ptr->sb_tot_cnt = picture_sb_width * picture_sb_height;

    return 0;
}
