/*
* Copyright(c) 2019 Netflix, Inc.
* Copyright (c) 2016, Alliance for Open Media. All rights reserved
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

#include <stdlib.h>
#include <limits.h>

#include "EbDefinitions.h"
#include "EbUtility.h"
#include "EbDecBlock.h"
#include "EbDecHandle.h"
#include "EbObuParse.h"
#include "EbDecMemInit.h"
#include "EbDecPicMgr.h"
#include "EbDecRestoration.h"
#include "EbDecParseObuUtil.h"
#include "EbDecParseFrame.h"
#include "EbInterPrediction.h"
#include "EbWarpedMotion.h"

/*TODO : Should be removed */
#include "EbDecInverseQuantize.h"
#include "EbDecProcessFrame.h"

#include "EbDecUtils.h"
#include "EbDecLF.h"

#include "EbDecCdef.h"
#include "EbLog.h"

void dec_av1_loop_filter_frame_mt(EbDecHandle *        dec_handle_ptr,
                                  EbPictureBufferDesc *recon_picture_buf,
                                  LfCtxt *lf_ctxt,
                                  int32_t plane_start,
                                  int32_t plane_end,
                                  DecThreadCtxt *thread_ctxt);

EbErrorType dec_system_resource_init(EbDecHandle *dec_handle_ptr, TilesInfo *tiles_info);

/* Scan through the Tiles to find Bitstream offsets */
void svt_av1_scan_tiles(EbDecHandle *dec_handle_ptr, TilesInfo *tiles_info, ObuHeader *obu_header,
                        Bitstrm *bs, uint32_t tg_start, uint32_t tg_end);

void svt_av1_queue_parse_jobs(EbDecHandle *dec_handle_ptr, TilesInfo *tiles_info);
void parse_frame_tiles(EbDecHandle *dec_handle_ptr, DecThreadCtxt *thread_ctxt);
void decode_frame_tiles(EbDecHandle *dec_handle_ptr, DecThreadCtxt *thread_ctxt);
void svt_av1_queue_lf_jobs(EbDecHandle *dec_handle_ptr);
void svt_av1_queue_cdef_jobs(EbDecHandle *dec_handle_ptr);
void svt_cdef_frame_mt(EbDecHandle *dec_handle_ptr, DecThreadCtxt *thread_ctxt);

void svt_av1_queue_lr_jobs(EbDecHandle *dec_handle_ptr);
void dec_av1_loop_restoration_filter_frame_mt(EbDecHandle *dec_handle,
                                              DecThreadCtxt *thread_ctxt);

#define CONFIG_MAX_DECODE_PROFILE 2

void dec_init_intra_predictors_12b_internal(void);

int remap_lr_type[4] = {RESTORE_NONE, RESTORE_SWITCHABLE, RESTORE_WIENER, RESTORE_SGRPROJ};

void av1_superres_upscale(Av1Common *cm, FrameHeader *frm_hdr, SeqHeader *seq_hdr,
                          EbPictureBufferDesc *recon_picture_src, int enable_flag);

/* Checks that the remaining bits start with a 1 and ends with 0s.
 * It consumes an additional byte, if already byte aligned before the check. */
int av1_check_trailing_bits(Bitstrm *bs) {
    // bit_offset is set to 0 (mod 8) when the reader is already byte aligned
    int bits_before_alignment = 8 - bs->bit_ofst % 8;
    int trailing              = dec_get_bits(bs, bits_before_alignment);
    if (trailing != (1 << (bits_before_alignment - 1))) return EB_Corrupt_Frame;
    return 0;
}

int byte_alignment(Bitstrm *bs) {
    while (bs->bit_ofst & 7) {
        if (dec_get_bits(bs, 1)) return EB_Corrupt_Frame;
    }
    return 0;
}

void compute_image_size(SeqHeader *seq_header, FrameHeader *frm) {
    frm->mi_cols = 2 * ((frm->frame_size.frame_width + 7) >> 3);
    frm->mi_rows = 2 * ((frm->frame_size.frame_height + 7) >> 3);
    frm->mi_stride =
        (ALIGN_POWER_OF_TWO(seq_header->max_frame_width, MAX_SB_SIZE_LOG2)) >> MI_SIZE_LOG2;
}

// Returns 1 when OBU type is valid, and 0 otherwise.
static int is_valid_obu_type(int obu_type) {
    int valid_type = 0;
    switch (obu_type) {
    case OBU_SEQUENCE_HEADER:
    case OBU_TEMPORAL_DELIMITER:
    case OBU_FRAME_HEADER:
    case OBU_TILE_GROUP:
    case OBU_METADATA:
    case OBU_FRAME:
    case OBU_REDUNDANT_FRAME_HEADER:
        //case OBU_TILE_LIST:
    case OBU_PADDING: valid_type = 1; break;
    default: break;
    }
    return valid_type;
}
// Read Operating point parameters
void read_operating_params_info(Bitstrm *bs, EbOperatingParametersInfo *op_info,
                                DecoderModelInfo *model_info, int index) {
    if (index > MAX_NUM_OPERATING_POINTS) return; // EB_DecUnsupportedBitstream;
    op_info->decoder_buffer_delay = dec_get_bits(bs, model_info->buffer_delay_length_minus_1 + 1);
    PRINT("decoder_buffer_delay", op_info->decoder_buffer_delay);
    op_info->encoder_buffer_delay = dec_get_bits(bs, model_info->buffer_delay_length_minus_1 + 1);
    PRINT("encoder_buffer_delay", op_info->encoder_buffer_delay);
    op_info->low_delay_mode_flag = dec_get_bits(bs, 1);
    PRINT("low_delay_mode_flag", op_info->low_delay_mode_flag);
}

// Read Timing information
void read_timing_info(Bitstrm *bs, EbTimingInfo *timing_info) {
    timing_info->num_units_in_display_tick = dec_get_bits(bs, 32);
    PRINT("num_units_in_display_tick", timing_info->num_units_in_display_tick);
    timing_info->time_scale = dec_get_bits(bs, 32);
    PRINT("time_scale", timing_info->time_scale);
    if (!timing_info->num_units_in_display_tick || !timing_info->time_scale)
        return; // EB_DecUnsupportedBitstream;
    timing_info->equal_picture_interval = dec_get_bits(bs, 1);
    PRINT("equal_picture_interval", timing_info->equal_picture_interval);
    if (timing_info->equal_picture_interval) {
        timing_info->num_ticks_per_picture = dec_get_bits_uvlc(bs) + 1;
        PRINT("num_ticks_per_picture", timing_info->num_ticks_per_picture);
        if ((timing_info->num_ticks_per_picture) == UINT32_MAX)
            return; // EB_DecUnsupportedBitstream;
    }
}

// Read Decoder model information
void read_decoder_model_info(Bitstrm *bs, DecoderModelInfo *model_info) {
    model_info->buffer_delay_length_minus_1 = dec_get_bits(bs, 5);
    PRINT("buffer_delay_length_minus_1", model_info->buffer_delay_length_minus_1);
    model_info->num_units_in_decoding_tick = dec_get_bits(bs, 32);
    PRINT("num_units_in_decoding_tick", model_info->num_units_in_decoding_tick);
    model_info->buffer_removal_time_length_minus_1 = dec_get_bits(bs, 5);
    PRINT("buffer_removal_time_length_minus_1",
        model_info->buffer_removal_time_length_minus_1);
    model_info->frame_presentation_time_length_minus_1 = dec_get_bits(bs, 5);
    PRINT("frame_presentation_time_length_minus_1",
          model_info->frame_presentation_time_length_minus_1);
}

// Read bit depth
void read_bit_depth(Bitstrm *bs, EbColorConfig *color_info, SeqHeader *seq_header) {
    uint8_t high_bitdepth = dec_get_bits(bs, 1);
    PRINT("high_bitdepth", high_bitdepth);
    if (seq_header->seq_profile == PROFESSIONAL_PROFILE && high_bitdepth) {
        uint8_t twelve_bit = dec_get_bits(bs, 1);
        PRINT("twelve_bit", twelve_bit);
        color_info->bit_depth = twelve_bit ? AOM_BITS_12 : AOM_BITS_10;
    } else if (seq_header->seq_profile <= PROFESSIONAL_PROFILE)
        color_info->bit_depth = high_bitdepth ? AOM_BITS_10 : AOM_BITS_8;
    else
        return; // EB_DecUnsupportedBitstream;
}

// Read Color configuration
void read_color_config(Bitstrm *bs, EbColorConfig *color_info, SeqHeader *seq_header) {
    read_bit_depth(bs, color_info, seq_header);
    color_info->mono_chrome = (seq_header->seq_profile != HIGH_PROFILE) ? dec_get_bits(bs, 1) : 0;
    PRINT("mono_chrome", color_info->mono_chrome);
    color_info->color_description_present_flag = dec_get_bits(bs, 1);
    PRINT("color_description_present_flag", color_info->color_description_present_flag);
    if (color_info->color_description_present_flag) {
        color_info->color_primaries = dec_get_bits(bs, 8);
        PRINT("color_primaries", color_info->color_primaries);
        color_info->transfer_characteristics = dec_get_bits(bs, 8);
        PRINT("transfer_characteristics", color_info->transfer_characteristics);
        color_info->matrix_coefficients = dec_get_bits(bs, 8);
        PRINT("matrix_coefficients", color_info->matrix_coefficients);
    } else {
        color_info->color_primaries          = EB_CICP_CP_UNSPECIFIED;
        color_info->transfer_characteristics = EB_CICP_TC_UNSPECIFIED;
        color_info->matrix_coefficients      = EB_CICP_MC_UNSPECIFIED;
    }
    if (color_info->mono_chrome) {
        color_info->color_range = dec_get_bits(bs, 1);
        PRINT("color_range", color_info->color_range);
        color_info->subsampling_y = color_info->subsampling_x = 1;
        color_info->chroma_sample_position                    = EB_CSP_UNKNOWN;
        color_info->separate_uv_delta_q                       = 0;
        return;
    } else if (color_info->color_primaries == EB_CICP_CP_BT_709 &&
               color_info->transfer_characteristics == EB_CICP_TC_SRGB &&
               color_info->matrix_coefficients == EB_CICP_MC_IDENTITY) {
        color_info->subsampling_x = color_info->subsampling_y = 0;
        color_info->color_range                               = 1; // assume full color-range
        if (!(seq_header->seq_profile == HIGH_PROFILE ||
              (seq_header->seq_profile == PROFESSIONAL_PROFILE &&
               color_info->bit_depth == AOM_BITS_12)))
            return; // EB_DecUnsupportedBitstream;
    } else {
        color_info->color_range = dec_get_bits(bs, 1);
        PRINT("color_range", color_info->color_range);
        if (seq_header->seq_profile == MAIN_PROFILE)
            color_info->subsampling_x = color_info->subsampling_y = 1; // 420 only
        else if (seq_header->seq_profile == HIGH_PROFILE)
            color_info->subsampling_x = color_info->subsampling_y = 0; // 444 only
        else {
            if (color_info->bit_depth == AOM_BITS_12) {
                color_info->subsampling_x = dec_get_bits(bs, 1);
                PRINT("subsampling_x", color_info->subsampling_x);
                if (color_info->subsampling_x) {
                    color_info->subsampling_y = dec_get_bits(bs, 1);
                    PRINT("subsampling_y", color_info->subsampling_y);
                } else
                    color_info->subsampling_y = 0;
            } else { // 422
                color_info->subsampling_x = 1;
                color_info->subsampling_y = 0;
            }
        }
        if (color_info->matrix_coefficients == EB_CICP_MC_IDENTITY &&
            (color_info->subsampling_x || color_info->subsampling_y))
            return; // EB_DecUnsupportedBitstream;
        if (color_info->subsampling_x && color_info->subsampling_y)
            color_info->chroma_sample_position = dec_get_bits(bs, 2);
        PRINT("chroma_sample_position", color_info->chroma_sample_position);
    }
    color_info->separate_uv_delta_q = dec_get_bits(bs, 1);
    PRINT("separate_uv_delta_q", color_info->separate_uv_delta_q);
}

void read_temporal_delimitor_obu(uint8_t *seen_frame_header) { *seen_frame_header = 0; }

// Read Sequence header
EbErrorType read_sequence_header_obu(Bitstrm *bs, SeqHeader *seq_header) {
    EbErrorType status;

    seq_header->seq_profile = (EbAv1SeqProfile)dec_get_bits(bs, 3);
    PRINT("seq_profile", seq_header->seq_profile);
    if (seq_header->seq_profile > CONFIG_MAX_DECODE_PROFILE) return EB_Corrupt_Frame;

    seq_header->still_picture = dec_get_bits(bs, 1);
    PRINT("still_picture", seq_header->still_picture);
    seq_header->reduced_still_picture_header = dec_get_bits(bs, 1);
    PRINT("reduced_still_picture_header", seq_header->reduced_still_picture_header);

    // Video must have reduced_still_picture_header = 0
    if (!seq_header->still_picture && seq_header->reduced_still_picture_header)
        return EB_DecUnsupportedBitstream;

    if (seq_header->reduced_still_picture_header) {
        seq_header->timing_info.timing_info_present    = 0;
        seq_header->decoder_model_info_present_flag    = 0;
        seq_header->initial_display_delay_present_flag = 0;
        seq_header->operating_points_cnt_minus_1       = 0;
        seq_header->operating_point[0].op_idc          = 0;
        seq_header->operating_point[0].seq_level_idx   = dec_get_bits(bs, LEVEL_BITS);
        PRINT("seq_level_idx", seq_header->operating_point[0].seq_level_idx);
        if (!is_valid_seq_level_idx(seq_header->operating_point->seq_level_idx))
            return EB_Corrupt_Frame;
        seq_header->operating_point[0].seq_tier                                  = 0;
        seq_header->operating_point[0].decoder_model_present_for_this_op         = 0;
        seq_header->operating_point[0].initial_display_delay_present_for_this_op = 0;
    } else {
        seq_header->timing_info.timing_info_present = dec_get_bits(bs, 1);
        PRINT("timing_info_present_flag", seq_header->timing_info.timing_info_present);
        if (seq_header->timing_info.timing_info_present) {
            read_timing_info(bs, &seq_header->timing_info);
            seq_header->decoder_model_info_present_flag = dec_get_bits(bs, 1);
            PRINT("decoder_model_info_present_flag", seq_header->decoder_model_info_present_flag);
            if (seq_header->decoder_model_info_present_flag)
                read_decoder_model_info(bs, &seq_header->decoder_model_info);
        } else
            seq_header->decoder_model_info_present_flag = 0;
    }

    seq_header->initial_display_delay_present_flag = dec_get_bits(bs, 1);
    PRINT("initial_display_delay_present_flag", seq_header->initial_display_delay_present_flag);
    seq_header->operating_points_cnt_minus_1 = dec_get_bits(bs, OP_POINTS_CNT_MINUS_1_BITS);
    PRINT("operating_points_cnt_minus_1", seq_header->operating_points_cnt_minus_1);
    for (int i = 0; i <= seq_header->operating_points_cnt_minus_1; i++) {
        seq_header->operating_point[i].op_idc = dec_get_bits(bs, OP_POINTS_IDC_BITS);
        PRINT("operating_point_idc", seq_header->operating_point[i].op_idc);
        seq_header->operating_point[i].seq_level_idx = dec_get_bits(bs, LEVEL_BITS);
        PRINT("seq_level_idx", seq_header->operating_point[i].seq_level_idx);
        if (!is_valid_seq_level_idx(seq_header->operating_point[i].seq_level_idx))
            return EB_Corrupt_Frame;
        if (seq_header->operating_point[i].seq_level_idx > 7) {
            seq_header->operating_point[i].seq_tier = dec_get_bits(bs, 1);
            PRINT("seq_tier", seq_header->operating_point[i].seq_tier);
        } else
            seq_header->operating_point[i].seq_tier = 0;

        if (seq_header->decoder_model_info_present_flag) {
            seq_header->operating_point[i].decoder_model_present_for_this_op = dec_get_bits(bs, 1);
            PRINT("decoder_model_present_for_this_op",
                  seq_header->operating_point[i].decoder_model_present_for_this_op);
            if (seq_header->operating_point[i].decoder_model_present_for_this_op)
                read_operating_params_info(
                    bs,
                    &seq_header->operating_point[i].operating_parameters_info,
                    &seq_header->decoder_model_info,
                    i);
        } else
            seq_header->operating_point[i].decoder_model_present_for_this_op = 0;

        if (seq_header->initial_display_delay_present_flag) {
            seq_header->operating_point[i].initial_display_delay_present_for_this_op =
                dec_get_bits(bs, 1);
            PRINT("initial_display_delay_present_for_this_op",
                  seq_header->operating_point[i].initial_display_delay_present_for_this_op);
            if (seq_header->operating_point[i].initial_display_delay_present_for_this_op)
                seq_header->operating_point[i].initial_display_delay = dec_get_bits(bs, 4) + 1;
            PRINT("initial_display_delay",
                  seq_header->operating_point[i].initial_display_delay);
        }
    }

    // TODO: operating_point = choose_operating_point( )

    seq_header->frame_width_bits = dec_get_bits(bs, 4) + 1;
    PRINT("frame_width_bits", seq_header->frame_width_bits);
    seq_header->frame_height_bits = dec_get_bits(bs, 4) + 1;
    PRINT("frame_height_bits", seq_header->frame_height_bits);
    seq_header->max_frame_width = dec_get_bits(bs, (seq_header->frame_width_bits)) + 1;
    PRINT("max_frame_width", seq_header->max_frame_width);
    seq_header->max_frame_height = dec_get_bits(bs, (seq_header->frame_height_bits)) + 1;
    PRINT("max_frame_height", seq_header->max_frame_height);

    if (seq_header->reduced_still_picture_header)
        seq_header->frame_id_numbers_present_flag = 0;
    else {
        seq_header->frame_id_numbers_present_flag = dec_get_bits(bs, 1);
        PRINT("frame_id_numbers_present_flag", seq_header->frame_id_numbers_present_flag);
    }
    if (seq_header->frame_id_numbers_present_flag) {
        seq_header->delta_frame_id_length = dec_get_bits(bs, 4) + 2;
        PRINT("delta_frame_id_length", seq_header->delta_frame_id_length);
        seq_header->frame_id_length = dec_get_bits(bs, 3) + 1;
        PRINT("frame_id_length", seq_header->frame_id_length + seq_header->delta_frame_id_length);
        if (seq_header->frame_id_length - 1 > 16) return EB_Corrupt_Frame;
    }

    seq_header->use_128x128_superblock = dec_get_bits(bs, 1);
    seq_header->sb_size      = seq_header->use_128x128_superblock ? BLOCK_128X128 : BLOCK_64X64;
    seq_header->sb_mi_size   = seq_header->use_128x128_superblock ? 32 : 16;
    seq_header->sb_size_log2 = seq_header->use_128x128_superblock ? 7 : 6;
    PRINT("use_128x128_superblock", seq_header->use_128x128_superblock);
    seq_header->filter_intra_level = dec_get_bits(bs, 1);
    PRINT("filter_intra_level", seq_header->filter_intra_level);
    seq_header->enable_intra_edge_filter = dec_get_bits(bs, 1);
    PRINT("enable_intra_edge_filter", seq_header->enable_intra_edge_filter);

    if (seq_header->reduced_still_picture_header) {
        seq_header->enable_interintra_compound           = 0;
        seq_header->enable_masked_compound               = 0;
        seq_header->enable_warped_motion                 = 0;
        seq_header->enable_dual_filter                   = 0;
        seq_header->order_hint_info.enable_jnt_comp      = 0;
        seq_header->order_hint_info.enable_ref_frame_mvs = 0;
        seq_header->seq_force_screen_content_tools       = 2;
        seq_header->seq_force_integer_mv                 = 2;
        seq_header->order_hint_info.order_hint_bits      = 0;
    } else {
        seq_header->enable_interintra_compound = dec_get_bits(bs, 1);
        PRINT("enable_interintra_compound", seq_header->enable_interintra_compound);
        seq_header->enable_masked_compound = dec_get_bits(bs, 1);
        PRINT("enable_masked_compound", seq_header->enable_masked_compound);
        seq_header->enable_warped_motion = dec_get_bits(bs, 1);
        PRINT("enable_warped_motion", seq_header->enable_warped_motion);
        seq_header->enable_dual_filter = dec_get_bits(bs, 1);
        PRINT("enable_dual_filter", seq_header->enable_dual_filter);
        seq_header->order_hint_info.enable_order_hint = dec_get_bits(bs, 1);
        PRINT("enable_order_hint", seq_header->order_hint_info.enable_order_hint);
        if (seq_header->order_hint_info.enable_order_hint) {
            seq_header->order_hint_info.enable_jnt_comp = dec_get_bits(bs, 1);
            PRINT("enable_jnt_comp", seq_header->order_hint_info.enable_jnt_comp);
            seq_header->order_hint_info.enable_ref_frame_mvs = dec_get_bits(bs, 1);
            PRINT("enable_ref_frame_mvs", seq_header->order_hint_info.enable_ref_frame_mvs);
        } else {
            seq_header->order_hint_info.enable_jnt_comp      = 0;
            seq_header->order_hint_info.enable_ref_frame_mvs = 0;
        }
        if (dec_get_bits(bs, 1)) {
            PRINT_NAME("seq_choose_screen_content_tools");
            seq_header->seq_force_screen_content_tools = SELECT_SCREEN_CONTENT_TOOLS;
        } else
            seq_header->seq_force_screen_content_tools = dec_get_bits(bs, 1);
        PRINT("seq_force_screen_content_tools", seq_header->seq_force_screen_content_tools);
        if (seq_header->seq_force_screen_content_tools > 0) {
            if (dec_get_bits(bs, 1)) {
                PRINT_NAME("seq_choose_screen_content_tools");
                seq_header->seq_force_integer_mv = SELECT_INTEGER_MV;
            } else
                seq_header->seq_force_integer_mv = dec_get_bits(bs, 1);
            PRINT("seq_force_integer_mv", seq_header->seq_force_integer_mv);
        } else
            seq_header->seq_force_integer_mv = SELECT_INTEGER_MV;

        if (seq_header->order_hint_info.enable_order_hint) {
            seq_header->order_hint_info.order_hint_bits = dec_get_bits(bs, 3) + 1;
            PRINT("order_hint_bits", seq_header->order_hint_info.order_hint_bits - 1);
        } else
            seq_header->order_hint_info.order_hint_bits = 0;
    }
    seq_header->enable_superres = dec_get_bits(bs, 1);
    PRINT("enable_superres", seq_header->enable_superres);
    seq_header->cdef_level = dec_get_bits(bs, 1);
    PRINT("cdef_level", seq_header->cdef_level);
    seq_header->enable_restoration = dec_get_bits(bs, 1);
    PRINT("enable_restoration", seq_header->enable_restoration);

    read_color_config(bs, &seq_header->color_config, seq_header);
    seq_header->film_grain_params_present = dec_get_bits(bs, 1);
    PRINT("film_grain_params_present", seq_header->film_grain_params_present);
    status = av1_check_trailing_bits(bs);
    if (status != EB_ErrorNone) return status;
    return EB_ErrorNone;
}

// Read OBU header
EbErrorType read_obu_header(Bitstrm *bs, ObuHeader *header) {
    PRINT_NL;
    if (!bs || !header) return EB_ErrorBadParameter;

    header->size = 1;

    if (dec_get_bits(bs, 1) != 0) {
        // obu_forbidden_bit must be set to 0.
        return EB_Corrupt_Frame;
    }
    PRINT_NAME("obu_forbidden_bit");
    header->obu_type = (ObuType)dec_get_bits(bs, 4);
    PRINT("obu_type", header->obu_type);
    if (!is_valid_obu_type(header->obu_type)) return EB_Corrupt_Frame;

    header->obu_extension_flag = dec_get_bits(bs, 1);
    PRINT("obu_extension_flag", header->obu_extension_flag);
    header->obu_has_size_field = dec_get_bits(bs, 1);
    PRINT("obu_has_size_field", header->obu_has_size_field);

    if (dec_get_bits(bs, 1) != 0) {
        // obu_reserved_1bit must be set to 0
        return EB_Corrupt_Frame;
    }
    PRINT_NAME("obu_reserved_1bit");

    if (header->obu_extension_flag) {
        header->size += 1;

        header->temporal_id = dec_get_bits(bs, 3);
        PRINT("temporal_id", header->temporal_id);
        header->spatial_id = dec_get_bits(bs, 2);
        PRINT("spatial_id", header->spatial_id);

        if (dec_get_bits(bs, 3) != 0) {
            // extension_header_reserved_3bits must be set to 0
            return EB_Corrupt_Frame;
        }
        PRINT_NAME("obu_extension_header_reserved_3bits");
    } else
        header->temporal_id = header->spatial_id = 0;

    return EB_ErrorNone;
}

// Read OBU size
EbErrorType read_obu_size(Bitstrm *bs, size_t bytes_available, size_t *const obu_size,
                          size_t *const length_field_size) {
    size_t u_obu_size = 0;
    dec_get_bits_leb128(bs, bytes_available, &u_obu_size, length_field_size);

    if (u_obu_size > UINT32_MAX) return EB_Corrupt_Frame;
    *obu_size = u_obu_size;
    PRINT("obu_size", *obu_size);
    return EB_ErrorNone;
}

/** Reads OBU header and size */
EbErrorType read_obu_header_size(Bitstrm *bs, ObuHeader *header, size_t size,
                                 size_t *const length_size) {
    EbErrorType status;

    status = read_obu_header(bs, header);
    if (status != EB_ErrorNone) return status;

    if (header->obu_has_size_field) {
        status = read_obu_size(bs, size, &header->payload_size, length_size);
        if (status != EB_ErrorNone) return status;
    }

    return EB_ErrorNone;
}

void temporal_point_info(Bitstrm *bs, DecoderModelInfo *model_info, FrameHeader *frame_info) {
    int n                               = model_info->frame_presentation_time_length_minus_1 + 1;
    frame_info->frame_presentation_time = dec_get_bits(bs, n);
}

void superres_params(Bitstrm *bs, SeqHeader *seq_header, FrameHeader *frame_info) {
    int use_superres = seq_header->enable_superres ? dec_get_bits(bs, 1) : 0;
    PRINT_NAME("use_superres");

    if (use_superres)
        frame_info->frame_size.superres_denominator = dec_get_bits(bs, SUPERRES_SCALE_BITS) + SUPERRES_SCALE_DENOMINATOR_MIN;
    else
        frame_info->frame_size.superres_denominator = SCALE_NUMERATOR;
    frame_info->frame_size.superres_upscaled_width = frame_info->frame_size.frame_width;
    frame_info->frame_size.frame_width =
        (frame_info->frame_size.superres_upscaled_width * SCALE_NUMERATOR +
         (frame_info->frame_size.superres_denominator / 2)) /
        frame_info->frame_size.superres_denominator;

    if (frame_info->frame_size.superres_denominator != SCALE_NUMERATOR) {
        /* We need to ensure the constraint in "Appendix A" of the spec:
        * FrameWidth is greater than or equal to 16
        * FrameHeight is greater than or equal to 16
        For this, we clamp the downscaled dimension to at least 16. One
        exception: if original dimension itself was < 16, then we keep the
        downscaled dimension to be same as the original, to ensure that resizin
        is valid.*/
        const int min_w = MIN(16, frame_info->frame_size.superres_upscaled_width);
        frame_info->frame_size.frame_width = MAX(min_w, frame_info->frame_size.frame_width);
    }

    PRINT_FRAME("superres_scale_denominator", frame_info->frame_size.superres_denominator)
}

// Read frame size
static void read_frame_size(Bitstrm *bs, SeqHeader *seq_header, FrameHeader *frame_info,
                            int frame_size_override_flag) {
    if (frame_size_override_flag) {
        frame_info->frame_size.frame_width  = dec_get_bits(bs, seq_header->frame_width_bits) + 1;
        frame_info->frame_size.frame_height = dec_get_bits(bs, seq_header->frame_height_bits) + 1;
    } else {
        frame_info->frame_size.frame_width  = seq_header->max_frame_width;
        frame_info->frame_size.frame_height = seq_header->max_frame_height;
    }
    PRINT_FRAME("frame_width", frame_info->frame_size.frame_width);
    PRINT_FRAME("frame_height", frame_info->frame_size.frame_height);
    superres_params(bs, seq_header, frame_info);
    compute_image_size(seq_header, frame_info);

    assert((frame_info->frame_size.frame_width) <= seq_header->max_frame_width);
    assert((frame_info->frame_size.frame_height) <= seq_header->max_frame_height);
}

void read_render_size(Bitstrm *bs, FrameHeader *frame_info) {
    uint8_t render_and_frame_size_different;
    render_and_frame_size_different = dec_get_bits(bs, 1);
    PRINT_NAME("render_and_frame_size_different");
    if (render_and_frame_size_different == 1) {
        frame_info->frame_size.render_width  = dec_get_bits(bs, 16) + 1;
        frame_info->frame_size.render_height = dec_get_bits(bs, 16) + 1;
    } else {
        frame_info->frame_size.render_width  = frame_info->frame_size.superres_upscaled_width;
        frame_info->frame_size.render_height = frame_info->frame_size.frame_height;
    }
    PRINT_FRAME("render_width", frame_info->frame_size.render_width);
    PRINT_FRAME("render_height", frame_info->frame_size.render_height);
}

static void frame_size_with_refs(Bitstrm *bs, EbDecHandle *dec_handle_ptr,
                                 int frame_size_override_flag) {
    SeqHeader *  seq_header = &dec_handle_ptr->seq_header;
    FrameHeader *frame_info = &dec_handle_ptr->frame_header;

    int found_ref;
    for (int i = 0; i < REFS_PER_FRAME; i++) {
        found_ref = dec_get_bits(bs, 1);
        PRINT_FRAME("found_ref", found_ref);
        if (found_ref == 1) {
            FrameSize *  frame_size = &frame_info->frame_size;
            EbDecPicBuf *ref_buf    = get_ref_frame_buf(dec_handle_ptr, (i + 1));
            assert(ref_buf != NULL);

            frame_size->superres_upscaled_width = ref_buf->superres_upscaled_width;
            frame_size->frame_width             = ref_buf->superres_upscaled_width;
            frame_size->frame_height            = ref_buf->frame_height;
            frame_size->render_width            = ref_buf->render_width;
            frame_size->render_height           = ref_buf->render_height;
            break;
        }
    }
    if (found_ref == 0) {
        read_frame_size(bs, seq_header, frame_info, frame_size_override_flag);
        read_render_size(bs, frame_info);
    } else {
        superres_params(bs, seq_header, frame_info);
        compute_image_size(seq_header, frame_info);
    }
}

// Read IP filter parameters
void read_interpolation_filter(Bitstrm *bs, FrameHeader *frame_info) {
    uint8_t is_filter_switchable;
    is_filter_switchable = dec_get_bits(bs, 1);
    PRINT_NAME("is_filter_switchable");
    if (is_filter_switchable == 1)
        frame_info->interpolation_filter = SWITCHABLE;
    else
        frame_info->interpolation_filter = dec_get_bits(bs, 2);
    PRINT_FRAME("interpolation_filter", frame_info->interpolation_filter);
}

// Read Tile information
void read_tile_info(Bitstrm *bs, TilesInfo *tile_info, SeqHeader *seq_header,
                    FrameHeader *frame_info) {
    int      start_sb, i;
    int      sb_cols = seq_header->use_128x128_superblock ? ((frame_info->mi_cols + 31) >> 5)
                                                     : ((frame_info->mi_cols + 15) >> 4);
    int sb_rows = seq_header->use_128x128_superblock ? ((frame_info->mi_rows + 31) >> 5)
                                                     : ((frame_info->mi_rows + 15) >> 4);
    int sb_shift         = seq_header->use_128x128_superblock ? 5 : 4;
    int sb_size          = sb_shift + 2;
    int max_tile_area_sb = MAX_TILE_AREA >> (2 * sb_size);

    tile_info->max_tile_width_sb  = MAX_TILE_WIDTH >> sb_size;
    tile_info->max_tile_height_sb = (MAX_TILE_AREA / MAX_TILE_WIDTH) >> sb_size;
    tile_info->min_log2_tile_cols = tile_log2(tile_info->max_tile_width_sb, sb_cols);
    tile_info->max_log2_tile_cols = tile_log2(1, MIN(sb_cols, MAX_TILE_COLS));
    tile_info->max_log2_tile_rows = tile_log2(1, MIN(sb_rows, MAX_TILE_ROWS));
    tile_info->min_log2_tiles =
        MAX(tile_info->min_log2_tile_cols, tile_log2(max_tile_area_sb, sb_rows * sb_cols));
    tile_info->uniform_tile_spacing_flag = dec_get_bits(bs, 1);
    PRINT_FRAME("uniform_tile_spacing_flag", tile_info->uniform_tile_spacing_flag);
    if (tile_info->uniform_tile_spacing_flag) {
        tile_info->tile_cols_log2 = tile_info->min_log2_tile_cols;
        while (tile_info->tile_cols_log2 < tile_info->max_log2_tile_cols) {
            PRINT_NAME("increment_tile_cols_log2");
            if (dec_get_bits(bs, 1) == 1)
                tile_info->tile_cols_log2++;
            else
                break;
        }
        int tile_width_sb =
            (sb_cols + (1 << tile_info->tile_cols_log2) - 1) >> tile_info->tile_cols_log2;
        assert(tile_width_sb <= tile_info->max_tile_width_sb); // Bitstream conformance
        i = 0;
        for (start_sb = 0; start_sb < sb_cols; start_sb += tile_width_sb) {
            tile_info->tile_col_start_mi[i] = start_sb << sb_shift;
            i += 1;
        }
        tile_info->tile_col_start_mi[i] = frame_info->mi_cols;
        tile_info->tile_cols            = i;

        tile_info->min_log2_tile_rows =
            MAX(tile_info->min_log2_tiles - tile_info->tile_cols_log2, 0);
        tile_info->tile_rows_log2 = tile_info->min_log2_tile_rows;
        while (tile_info->tile_rows_log2 < tile_info->max_log2_tile_rows) {
            PRINT_NAME("Some read")
            if (dec_get_bits(bs, 1) == 1)
                tile_info->tile_rows_log2++;
            else
                break;
        }
        int tile_height_sb =
            (sb_rows + (1 << tile_info->tile_rows_log2) - 1) >> tile_info->tile_rows_log2;
        assert(tile_height_sb <= tile_info->max_tile_height_sb); // Bitstream conformance
        i = 0;
        for (start_sb = 0; start_sb < sb_rows; start_sb += tile_height_sb) {
            tile_info->tile_row_start_mi[i] = start_sb << sb_shift;
            i += 1;
        }
        tile_info->tile_row_start_mi[i] = frame_info->mi_rows;
        tile_info->tile_rows            = i;
    } else {
        int widest_tile_sb = 0;
        start_sb           = 0;
        for (i = 0; start_sb < sb_cols; i++) {
            tile_info->tile_col_start_mi[i] = start_sb << sb_shift;
            int max_width                   = MIN(sb_cols - start_sb, tile_info->max_tile_width_sb);
            uint32_t width_in_sbs_minus_1        = dec_get_bits_ns(bs, max_width);
            PRINT("width_in_sbs_minus_1", width_in_sbs_minus_1)
            int size_sb    = width_in_sbs_minus_1 + 1;
            widest_tile_sb = MAX(size_sb, widest_tile_sb);
            start_sb += size_sb;
        }
        assert(start_sb == sb_cols); // Bitstream conformance

        tile_info->tile_col_start_mi[i] = frame_info->mi_cols;
        tile_info->tile_cols            = i;
        tile_info->tile_cols_log2       = tile_log2(1, tile_info->tile_cols);
        if (tile_info->min_log2_tiles > 0)
            max_tile_area_sb = (sb_rows * sb_cols) >> (tile_info->min_log2_tiles + 1);
        else
            max_tile_area_sb = sb_rows * sb_cols;
        assert(widest_tile_sb > 0);
        tile_info->max_tile_height_sb = MAX(max_tile_area_sb / widest_tile_sb, 1);

        start_sb = 0;
        for (i = 0; start_sb < sb_rows; i++) {
            tile_info->tile_row_start_mi[i] = start_sb << sb_shift;
            int max_height            = MIN(sb_rows - start_sb, tile_info->max_tile_height_sb);
            uint32_t height_in_sbs_minus_1  = dec_get_bits_ns(bs, max_height);
            PRINT("height_in_sbs_minus_1", height_in_sbs_minus_1)
            start_sb += height_in_sbs_minus_1 + 1;
        }
        assert(start_sb == sb_rows); // Bitstream conformance

        tile_info->tile_row_start_mi[i] = frame_info->mi_rows;
        tile_info->tile_rows            = i;
        tile_info->tile_rows_log2       = tile_log2(1, tile_info->tile_rows);
    }

    // Bitstream conformance
    assert(tile_info->tile_cols <= MAX_TILE_ROWS);
    assert(tile_info->tile_rows <= MAX_TILE_COLS);

    if (tile_info->tile_cols_log2 > 0 || tile_info->tile_rows_log2 > 0) {
        tile_info->context_update_tile_id =
            dec_get_bits(bs, tile_info->tile_rows_log2 + tile_info->tile_cols_log2);
        PRINT("context_update_tile_id", tile_info->context_update_tile_id)
        tile_info->tile_size_bytes = dec_get_bits(bs, 2) + 1;
        PRINT("tile_size_bytes", tile_info->tile_size_bytes)
    } else
        tile_info->context_update_tile_id = 0;
    assert(tile_info->context_update_tile_id < (tile_info->tile_cols * tile_info->tile_rows));
}

uint8_t read_delta_q(Bitstrm *bs) {
    uint8_t delta_q;
    if (dec_get_bits(bs, 1))
        delta_q = dec_get_bits_su(bs, 7);
    else
        delta_q = 0;
    return delta_q;
}

void read_frame_delta_q_params(Bitstrm *bs, FrameHeader *frame_info) {
    frame_info->delta_q_params.delta_q_res     = 0;
    frame_info->delta_q_params.delta_q_present = 0;
    if (frame_info->quantization_params.base_q_idx > 0) {
        frame_info->delta_q_params.delta_q_present = dec_get_bits(bs, 1);
    }
    if (frame_info->delta_q_params.delta_q_present) {
        frame_info->delta_q_params.delta_q_res = dec_get_bits(bs, 2);
        PRINT_FRAME("delta_q_res", 1 << frame_info->delta_q_params.delta_q_res);
    }
    PRINT_FRAME("delta_q_present", frame_info->delta_q_params.delta_q_present);
}

void read_frame_delta_lf_params(Bitstrm *bs, FrameHeader *frame_info) {
    frame_info->delta_lf_params.delta_lf_present = 0;
    frame_info->delta_lf_params.delta_lf_res     = 0;
    frame_info->delta_lf_params.delta_lf_multi   = 0;
    if (frame_info->delta_q_params.delta_q_present) {
        if (!frame_info->allow_intrabc) {
            frame_info->delta_lf_params.delta_lf_present = dec_get_bits(bs, 1);
            PRINT_FRAME("delta_lf_present_flag", frame_info->delta_lf_params.delta_lf_present);
        }
        if (frame_info->delta_lf_params.delta_lf_present) {
            frame_info->delta_lf_params.delta_lf_res   = dec_get_bits(bs, 2);
            frame_info->delta_lf_params.delta_lf_multi = dec_get_bits(bs, 1);
            PRINT_FRAME("delta_lf_res", frame_info->delta_lf_params.delta_lf_res);
            PRINT_FRAME("delta_lf_multi", frame_info->delta_lf_params.delta_lf_multi);
        }
    }
}

void read_quantization_params(Bitstrm *bs, QuantizationParams *quant_params,
                              EbColorConfig *color_info, int num_planes) {
    quant_params->base_q_idx = dec_get_bits(bs, 8);
    PRINT_FRAME("base_q_idx", quant_params->base_q_idx);
    quant_params->delta_q_dc[AOM_PLANE_Y] = read_delta_q(bs);
    PRINT_FRAME("delta_q_y_dc", quant_params->delta_q_dc[0]);
    quant_params->delta_q_ac[AOM_PLANE_Y] = 0;
    if (num_planes > 1) {
        uint8_t diff_uv_delta = 0;
        if (color_info->separate_uv_delta_q) {
            diff_uv_delta = dec_get_bits(bs, 1);
            PRINT_FRAME("diff_uv_delta", diff_uv_delta);
        }
        quant_params->delta_q_dc[AOM_PLANE_U] = read_delta_q(bs); //U Plane
        quant_params->delta_q_ac[AOM_PLANE_U] = read_delta_q(bs); //U Plane
        if (diff_uv_delta) {
            quant_params->delta_q_dc[AOM_PLANE_V] = read_delta_q(bs); //V Plane
            quant_params->delta_q_ac[AOM_PLANE_V] = read_delta_q(bs); //V Plane
        } else {
            quant_params->delta_q_dc[AOM_PLANE_V] = quant_params->delta_q_dc[AOM_PLANE_U];
            quant_params->delta_q_ac[AOM_PLANE_V] = quant_params->delta_q_ac[AOM_PLANE_U];
        }
    } else {
        quant_params->delta_q_dc[AOM_PLANE_U] = 0;
        quant_params->delta_q_ac[AOM_PLANE_U] = 0;
        quant_params->delta_q_dc[AOM_PLANE_V] = 0;
        quant_params->delta_q_ac[AOM_PLANE_V] = 0;
    }
    PRINT_FRAME("u_dc_delta_q", quant_params->delta_q_dc[AOM_PLANE_U]);
    PRINT_FRAME("u_ac_delta_q", quant_params->delta_q_ac[AOM_PLANE_U]);
    PRINT_FRAME("v_dc_delta_q", quant_params->delta_q_dc[AOM_PLANE_V]);
    PRINT_FRAME("v_ac_delta_q", quant_params->delta_q_ac[AOM_PLANE_V]);
    quant_params->using_qmatrix = dec_get_bits(bs, 1);
    PRINT_FRAME("using_qmatrix", quant_params->using_qmatrix);
    if (quant_params->using_qmatrix) {
        quant_params->qm[AOM_PLANE_Y] = dec_get_bits(bs, 4);
        quant_params->qm[AOM_PLANE_U] = dec_get_bits(bs, 4);
        if (!color_info->separate_uv_delta_q)
            quant_params->qm[AOM_PLANE_V] = quant_params->qm[AOM_PLANE_U];
        else
            quant_params->qm[AOM_PLANE_V] = dec_get_bits(bs, 4);
    } else {
        quant_params->qm[AOM_PLANE_Y] = 0;
        quant_params->qm[AOM_PLANE_U] = 0;
        quant_params->qm[AOM_PLANE_V] = 0;
    }
    PRINT_FRAME("qm_y", quant_params->qm[AOM_PLANE_Y]);
    PRINT_FRAME("qm_u", quant_params->qm[AOM_PLANE_U]);
    PRINT_FRAME("qm_v", quant_params->qm[AOM_PLANE_V]);
}

static INLINE void segfeatures_copy(SegmentationParams *dst, SegmentationParams *src) {
    int i, j;
    for (i = 0; i < MAX_SEGMENTS; i++) {
        for (j = 0; j < SEG_LVL_MAX; j++) {
            dst->feature_data[i][j]    = src->feature_data[i][j];
            dst->feature_enabled[i][j] = src->feature_enabled[i][j];
        }
    }
    dst->seg_id_pre_skip    = src->seg_id_pre_skip;
    dst->last_active_seg_id = src->last_active_seg_id;
}

void read_segmentation_params(Bitstrm *bs, EbDecHandle *dec_handle_ptr, FrameHeader *frame_info) {
    SegmentationParams *seg_params = &frame_info->segmentation_params;
    EbDecPicBuf *       cur_buf    = dec_handle_ptr->cur_pic_buf[0];
    EbDecPicBuf *       prev_buf = dec_handle_ptr->prev_frame;

    seg_params->segmentation_enabled = dec_get_bits(bs, 1);
    PRINT_FRAME("segmentation_enabled", seg_params->segmentation_enabled);
    if (!seg_params->segmentation_enabled) {
        if (cur_buf->segment_maps)
            memset(cur_buf->segment_maps, 0, (frame_info->mi_rows * frame_info->mi_cols));

        memset(seg_params, 0, sizeof(*seg_params));
        segfeatures_copy(&cur_buf->seg_params, seg_params);
        return;
    }
    {
        EbDecPicBuf *prev_frame   = dec_handle_ptr->prev_frame;
        dec_handle_ptr->cm.last_frame_seg_map = NULL;
        if (prev_frame) {
            uint32_t prev_mi_cols = 2 * ((prev_frame->frame_width + 7) >> 3);
            uint32_t prev_mi_rows = 2 * ((prev_frame->frame_height + 7) >> 3);
            if (seg_params->segmentation_enabled      &&
                (frame_info->mi_rows == prev_mi_rows) &&
                (frame_info->mi_cols == prev_mi_cols))
            {
                dec_handle_ptr->cm.last_frame_seg_map = prev_frame->segment_maps;
            }
        }
    }
    if (frame_info->primary_ref_frame == PRIMARY_REF_NONE) {
        seg_params->segmentation_update_map      = 1;
        seg_params->segmentation_temporal_update = 0;
        seg_params->segmentation_update_data     = 1;
    } else {
        seg_params->segmentation_update_map = dec_get_bits(bs, 1);
        seg_params->segmentation_temporal_update = seg_params->segmentation_update_map
            ? dec_get_bits(bs, 1)
            : 0;
        seg_params->segmentation_update_data = dec_get_bits(bs, 1);
    }
    PRINT_FRAME("segmentation_update_map", seg_params->segmentation_update_map);
    PRINT_FRAME("segmentation_temporal_update", seg_params->segmentation_temporal_update);
    PRINT_FRAME("segmentation_update_data", seg_params->segmentation_update_data);
    if (seg_params->segmentation_update_data == 1) {
        for (int i = 0; i < MAX_SEGMENTS; i++) {
            for (int j = 0; j < SEG_LVL_MAX; j++) {
                int feature_enabled = dec_get_bits(bs, 1);
                PRINT_FRAME("feature_enabled", feature_enabled);
                seg_params->feature_enabled[i][j] = feature_enabled;
                int clipped_value                 = 0;
                if (feature_enabled == 1) {
                    const int bits_to_read = segmentation_feature_bits[j];
                    const int limit        = segmentation_feature_max[j];
                    if (segmentation_feature_signed[j] == 1) {
                        int feature_value = dec_get_bits_su(bs, 1 + bits_to_read);
                        clipped_value = CLIP3(-limit, limit, feature_value);
                    } else {
                        int feature_value = dec_get_bits(bs, bits_to_read);
                        clipped_value = CLIP3(0, limit, feature_value);
                    }
                    PRINT_FRAME("data", clipped_value)
                }
                if (clipped_value < 0) {
                    assert(segmentation_feature_signed[j]);
                    assert(-clipped_value <= segmentation_feature_max[j]);
                } else
                    assert(clipped_value <= segmentation_feature_max[j]);
                seg_params->feature_data[i][j] = clipped_value;
            }
        }
    } else if (prev_buf) {
        segfeatures_copy(seg_params, &prev_buf->seg_params);
    }
    segfeatures_copy(&cur_buf->seg_params, seg_params);

    seg_params->last_active_seg_id = 0;
    seg_params->seg_id_pre_skip    = 0;
    for (int i = 0; i < MAX_SEGMENTS; i++) {
        for (int j = 0; j < SEG_LVL_MAX; j++) {
            if (seg_params->feature_enabled[i][j]) {
                seg_params->last_active_seg_id = i;
                if (j >= SEG_LVL_REF_FRAME) seg_params->seg_id_pre_skip = 1;
            }
        }
    }
}

static void av1_set_default_ref_and_mode_deltas(int8_t *ref_deltas, int8_t *mode_deltas) {
    assert(ref_deltas != NULL);
    assert(mode_deltas != NULL);

    ref_deltas[INTRA_FRAME] = 1;
    ref_deltas[LAST_FRAME] = 0;
    ref_deltas[LAST2_FRAME] = 0;
    ref_deltas[LAST3_FRAME] = 0;
    ref_deltas[BWDREF_FRAME] = 0;
    ref_deltas[GOLDEN_FRAME] = -1;
    ref_deltas[ALTREF_FRAME] = -1;
    ref_deltas[ALTREF2_FRAME] = -1;

    mode_deltas[0] = 0;
    mode_deltas[1] = 0;
}

void read_loop_filter_params(Bitstrm *bs, EbDecHandle *dec_handle, int num_planes) {
    FrameHeader *frame_info = &dec_handle->frame_header;
    struct LoopFilter *lf = &frame_info->loop_filter_params;

    if (frame_info->coded_lossless || frame_info->allow_intrabc) {
        lf->filter_level[0] = 0;
        lf->filter_level[1] = 0;

        av1_set_default_ref_and_mode_deltas(dec_handle->cur_pic_buf[0]->ref_deltas,
                                            dec_handle->cur_pic_buf[0]->mode_deltas);
        return;
    }

    if (dec_handle->prev_frame) {
        // write deltas to frame buffer
        svt_memcpy(lf->ref_deltas, dec_handle->prev_frame->ref_deltas, REF_FRAMES);
        svt_memcpy(lf->mode_deltas, dec_handle->prev_frame->mode_deltas, MAX_MODE_LF_DELTAS);
    }
    else {
        av1_set_default_ref_and_mode_deltas(lf->ref_deltas, lf->mode_deltas);
    }

    lf->filter_level[0] = dec_get_bits(bs, 6);
    lf->filter_level[1] = dec_get_bits(bs, 6);
    PRINT_FRAME("loop_filter_level[0]", lf->filter_level[0]);
    PRINT_FRAME("loop_filter_level[1]", lf->filter_level[1]);
    if (num_planes > 1) {
        if (lf->filter_level[0] || lf->filter_level[1]) {
            lf->filter_level_u = dec_get_bits(bs, 6);
            lf->filter_level_v = dec_get_bits(bs, 6);
            PRINT_FRAME("loop_filter_level[2]", lf->filter_level_u);
            PRINT_FRAME("loop_filter_level[3]", lf->filter_level_v);
        }
    }
    lf->sharpness_level = dec_get_bits(bs, 3);
    lf->mode_ref_delta_enabled = dec_get_bits(bs, 1);
    PRINT_FRAME("loop_filter_sharpness", lf->sharpness_level);
    PRINT_FRAME("loop_filter_delta_enabled", lf->mode_ref_delta_enabled);

    if (lf->mode_ref_delta_enabled == 1) {
        lf->mode_ref_delta_update = dec_get_bits(bs, 1);
        PRINT_FRAME("loop_filter_delta_update", lf->mode_ref_delta_update);

        if (lf->mode_ref_delta_update == 1) {
            for (int i = 0; i < TOTAL_REFS_PER_FRAME; i++) {
                PRINT_NAME("Some read");
                if (dec_get_bits(bs, 1) == 1) {
                    lf->ref_deltas[i] = dec_get_bits_su(bs, 1 + 6);
                    PRINT_FRAME("lf_ref_deltas[i]", lf->ref_deltas[i]);
                }
            }
            for (int i = 0; i < 2; i++) {
                PRINT_NAME("Some read");
                if (dec_get_bits(bs, 1) == 1) {
                    lf->mode_deltas[i] = dec_get_bits_su(bs, 1 + 6);
                    PRINT_FRAME("lf_mode_deltas[i]", lf->mode_deltas[i]);
                }
            }
        }
    }

    /*write deltas to prev_frame buffer*/
    svt_memcpy(dec_handle->cur_pic_buf[0]->ref_deltas, lf->ref_deltas, REF_FRAMES);
    svt_memcpy(dec_handle->cur_pic_buf[0]->mode_deltas, lf->mode_deltas, MAX_MODE_LF_DELTAS);
}

void read_tx_mode(Bitstrm *bs, FrameHeader *frame_info) {
    if (frame_info->coded_lossless == 1)
        frame_info->tx_mode = ONLY_4X4;
    else {
        if (dec_get_bits(bs, 1))
            frame_info->tx_mode = TX_MODE_SELECT;
        else
            frame_info->tx_mode = TX_MODE_LARGEST;
    }
    PRINT_FRAME("tx_mode", frame_info->tx_mode);
}

void read_lr_params(Bitstrm *bs, FrameHeader *frame_info, SeqHeader *seq_header, int num_planes) {
    int uses_lr, uses_chroma_lr;

    if (frame_info->coded_lossless || frame_info->allow_intrabc ||
        !seq_header->enable_restoration) {
        frame_info->lr_params[0].frame_restoration_type = RESTORE_NONE;
        frame_info->lr_params[1].frame_restoration_type = RESTORE_NONE;
        frame_info->lr_params[2].frame_restoration_type = RESTORE_NONE;
        return;
    }
    uses_lr        = 0;
    uses_chroma_lr = 0;
    for (int i = 0; i < num_planes; i++) {
        int lr_type = dec_get_bits(bs, 2);
        PRINT_NAME("lr_type")
        frame_info->lr_params[i].frame_restoration_type = remap_lr_type[lr_type];
        PRINT_FRAME("frame_restoration_type", frame_info->lr_params[i].frame_restoration_type);
        if (frame_info->lr_params[i].frame_restoration_type != RESTORE_NONE) {
            uses_lr = 1;
            if (i > 0) uses_chroma_lr = 1;
        }
    }
    if (uses_lr) {
        int lr_unit_shift = dec_get_bits(bs, 1);
        if (seq_header->use_128x128_superblock) {
            lr_unit_shift++;
        } else if (lr_unit_shift) {
            int lr_unit_extra_shift = dec_get_bits(bs, 1);
            PRINT_FRAME("lr_unit_extra_shift", lr_unit_extra_shift);
            lr_unit_shift += lr_unit_extra_shift;
        }
        frame_info->lr_params[0].loop_restoration_size =
            (RESTORATION_TILESIZE_MAX >> (2 - lr_unit_shift));
        frame_info->lr_params[0].lr_size_log2 = RESTORATION_UNIT_OFFSET - (2 - lr_unit_shift);
        PRINT_FRAME("restoration_unit_size", frame_info->lr_params[0].loop_restoration_size);
        int lr_uv_shift = seq_header->color_config.subsampling_x &&
                seq_header->color_config.subsampling_y && uses_chroma_lr
            ? dec_get_bits(bs, 1)
            : 0;

        frame_info->lr_params[1].loop_restoration_size =
            frame_info->lr_params[0].loop_restoration_size >> lr_uv_shift;
        frame_info->lr_params[2].loop_restoration_size =
            frame_info->lr_params[0].loop_restoration_size >> lr_uv_shift;
        frame_info->lr_params[1].lr_size_log2 = frame_info->lr_params[0].lr_size_log2 - lr_uv_shift;
        frame_info->lr_params[2].lr_size_log2 = frame_info->lr_params[1].lr_size_log2;
    } else {
        frame_info->lr_params[0].loop_restoration_size = RESTORATION_TILESIZE_MAX;
        frame_info->lr_params[1].loop_restoration_size = RESTORATION_TILESIZE_MAX;
        frame_info->lr_params[2].loop_restoration_size = RESTORATION_TILESIZE_MAX;
        frame_info->lr_params[0].lr_size_log2          = RESTORATION_UNIT_OFFSET;
        frame_info->lr_params[1].lr_size_log2          = RESTORATION_UNIT_OFFSET;
        frame_info->lr_params[2].lr_size_log2          = RESTORATION_UNIT_OFFSET;
    }
    PRINT_FRAME("cm->rst_info[1].restoration_unit_size",
                frame_info->lr_params[1].loop_restoration_size);
}

void read_frame_cdef_params(Bitstrm *bs, FrameHeader *frame_info, SeqHeader *seq_header,
                            int num_planes) {
    int i;
    if (frame_info->coded_lossless || frame_info->allow_intrabc || !seq_header->cdef_level) {
        frame_info->cdef_params.cdef_bits           = 0;
        frame_info->cdef_params.cdef_y_strength[0]  = 0;
        frame_info->cdef_params.cdef_y_strength[4]  = 0;
        frame_info->cdef_params.cdef_uv_strength[0] = 0;
        frame_info->cdef_params.cdef_uv_strength[4] = 0;
        frame_info->cdef_params.cdef_damping        = 3;
        return;
    }
    frame_info->cdef_params.cdef_damping = dec_get_bits(bs, 2) + 3;
    frame_info->cdef_params.cdef_bits    = dec_get_bits(bs, 2);
    PRINT_FRAME("cdef_damping", frame_info->cdef_params.cdef_damping);
    PRINT_FRAME("cdef_bits", frame_info->cdef_params.cdef_bits);
    for (i = 0; i < (1 << frame_info->cdef_params.cdef_bits); i++) {
        frame_info->cdef_params.cdef_y_strength[i] = dec_get_bits(bs, 6);
        PRINT_FRAME("Primary Y cdef", frame_info->cdef_params.cdef_y_strength[i]);

        if (num_planes > 1) {
            frame_info->cdef_params.cdef_uv_strength[i] = dec_get_bits(bs, 6);
            PRINT_FRAME("Primary UV cdef", frame_info->cdef_params.cdef_uv_strength[i]);
        }
    }
}

int decode_subexp(Bitstrm *bs, int numSyms) {
    int i = 0, mk = 0, k = 3;

    while (1) {
        int b2 = i ? k + i - 1 : k;
        int a = 1 << b2;
        if (numSyms <= mk + 3 * a) {
            PRINT_NAME("subexp_final_bits");
            return dec_get_bits_ns(bs, numSyms - mk) + mk;
        } else {
            PRINT_NAME("subexp_more_bits");
            if (dec_get_bits(bs, 1)) {
                i++;
                mk += a;
            } else {
                PRINT_NAME("subexp_bits");
                return dec_get_bits(bs, b2) + mk;
            }
        }
    }
    return 0;
}

int decode_unsigned_subexp_with_ref(Bitstrm *bs, int mx, int r) {
    int v = decode_subexp(bs, mx);
    if ((r << 1) <= mx)
        return inverse_recenter(r, v);
    else
        return mx - 1 - inverse_recenter(mx - 1 - r, v);
}

int decode_signed_subexp_with_ref(Bitstrm *bs, int low, int high, int r) {
    int x = decode_unsigned_subexp_with_ref(bs, high - low, r - low);
    return x + low;
}

void read_global_param(Bitstrm *bs, EbDecHandle *dec_handle, TransformationType type, int ref_idx,
                       int idx, FrameHeader *frame_info) {
    GlobalMotionParams
        prev_gm_params[ALTREF_FRAME +
                       1]; // Need to initialize in setup_past_independence() section: 6.8.2
    for (int ref = LAST_FRAME; ref <= ALTREF_FRAME; ref++)
        for (int i = 0; i <= 5; i++)
            prev_gm_params[ref].gm_params[i] = ((i % 3 == 2) ? 1 << WARPEDMODEL_PREC_BITS : 0);

    int abs_bits  = GM_ABS_ALPHA_BITS, prec_diff, round, sub, mx, r;
    int prec_bits = GM_ALPHA_PREC_BITS;
    if (idx < 2) {
        if (type == TRANSLATION) {
            abs_bits  = GM_ABS_TRANS_ONLY_BITS - !frame_info->allow_high_precision_mv;
            prec_bits = GM_TRANS_ONLY_PREC_BITS - !frame_info->allow_high_precision_mv;
        } else {
            abs_bits  = GM_ABS_TRANS_BITS;
            prec_bits = GM_TRANS_PREC_BITS;
        }
    }

    prec_diff = WARPEDMODEL_PREC_BITS - prec_bits;
    round     = (idx % 3) == 2 ? (1 << WARPEDMODEL_PREC_BITS) : 0;
    sub       = (idx % 3) == 2 ? (1 << prec_bits) : 0;
    mx        = (1 << abs_bits);

    EbDecPicBuf *cur_buf  = dec_handle->cur_pic_buf[0];
    EbDecPicBuf *prev_buf = dec_handle->prev_frame;

    GlobalMotionParams *gm_params = prev_buf != NULL ? prev_buf->global_motion : prev_gm_params;
    r                             = (gm_params[ref_idx].gm_params[idx] >> prec_diff) - sub;
    cur_buf->global_motion[ref_idx].gm_params[idx] =
        (decode_signed_subexp_with_ref(bs, -mx, mx + 1, r) << prec_diff) + round;
}

void read_global_motion_params(Bitstrm *bs, EbDecHandle *dec_handle, FrameHeader *frame_info,
                               int frame_is_intra) {
    int                ref, i;
    TransformationType type;
    EbDecPicBuf *      cur_buf = dec_handle->cur_pic_buf[0];
    for (ref = LAST_FRAME; ref <= ALTREF_FRAME; ref++) {
        cur_buf->global_motion[ref].gm_type = IDENTITY;
        for (i = 0; i < 6; i++) {
            cur_buf->global_motion[ref].gm_params[i] =
                ((i % 3 == 2) ? 1 << WARPEDMODEL_PREC_BITS : 0);
        }
    }
    if (frame_is_intra) return;
    for (ref = LAST_FRAME; ref <= ALTREF_FRAME; ref++) {
        PRINT_NAME("Some read");
        if (dec_get_bits(bs, 1)) {
            PRINT_NAME("Some read");
            if (dec_get_bits(bs, 1))
                type = ROTZOOM;
            else
                type = dec_get_bits(bs, 1) ? TRANSLATION : AFFINE;
        } else
            type = IDENTITY;
        PRINT_FRAME("Transform_type", type);

        cur_buf->global_motion[ref].gm_type = type;

        if (type >= ROTZOOM) {
            read_global_param(bs, dec_handle, type, ref, 2, frame_info);
            read_global_param(bs, dec_handle, type, ref, 3, frame_info);
        }
        if (type >= AFFINE) {
            read_global_param(bs, dec_handle, type, ref, 4, frame_info);
            read_global_param(bs, dec_handle, type, ref, 5, frame_info);
        } else {
            cur_buf->global_motion[ref].gm_params[4] = -cur_buf->global_motion[ref].gm_params[3];
            cur_buf->global_motion[ref].gm_params[5] = cur_buf->global_motion[ref].gm_params[2];
        }
        if (type >= TRANSLATION) {
            read_global_param(bs, dec_handle, type, ref, 0, frame_info);
            read_global_param(bs, dec_handle, type, ref, 1, frame_info);
        }

        /* TODO: Can we remove one of the type? */
        /* Convert to EbWarpedMotionParams type */
        {
            EbWarpedMotionParams *wm_global =
                &dec_handle->master_frame_buf.cur_frame_bufs[0].global_motion_warp[ref];
            wm_global->wmtype = cur_buf->global_motion[ref].gm_type;
            svt_memcpy(wm_global->wmmat,
                   cur_buf->global_motion[ref].gm_params,
                   sizeof(cur_buf->global_motion[ref].gm_params));
            int return_val = svt_get_shear_params(wm_global);
            assert(1 == return_val);
            (void)return_val;
        }
    }
}

uint8_t read_frame_reference_mode(Bitstrm *bs, int frame_is_intra) {
    if (frame_is_intra)
        return SINGLE_REFERENCE;
    else
        return dec_get_bits(bs, 1);
}

// Read skip mode paramters
void read_skip_mode_params(Bitstrm *bs, FrameHeader *frame_info, int frame_is_intra,
                           SeqHeader *seq_header, int reference_select) {
    if (frame_is_intra || !reference_select || !seq_header->order_hint_info.enable_order_hint)
        frame_info->skip_mode_params.skip_mode_allowed = 0;
    else {
        int forward_idx = -1, backward_idx = -1, forward_hint = -1, backward_hint = INT_MAX,
            second_forward_hint = -1;
        //frame_info->skip_mode_params.skip_mode_allowed = 1;
        for (int i = 0; i < REFS_PER_FRAME; i++) {
            int ref_hint = frame_info->ref_order_hint[frame_info->ref_frame_idx[i]];
            if (get_relative_dist(&seq_header->order_hint_info, ref_hint, frame_info->order_hint) <
                0) {
                if (forward_idx < 0 ||
                    get_relative_dist(&seq_header->order_hint_info, ref_hint, forward_hint) > 0) {
                    forward_idx  = i;
                    forward_hint = ref_hint;
                }
            } else if (get_relative_dist(
                           &seq_header->order_hint_info, ref_hint, frame_info->order_hint) > 0) {
                if (backward_idx < 0 ||
                    get_relative_dist(&seq_header->order_hint_info, ref_hint, backward_hint) < 0) {
                    backward_idx  = i;
                    backward_hint = ref_hint;
                }
            }
        }
        if (forward_idx < 0)
            frame_info->skip_mode_params.skip_mode_allowed = 0;
        else if (backward_idx >= 0) {
            frame_info->skip_mode_params.skip_mode_allowed = 1;
            frame_info->skip_mode_params.ref_frame_idx_0 =
                LAST_FRAME + MIN(forward_idx, backward_idx);
            frame_info->skip_mode_params.ref_frame_idx_1 =
                LAST_FRAME + MAX(forward_idx, backward_idx);
        } else {
            int second_forward_idx = -1;
            for (int i = 0; i < REFS_PER_FRAME; i++) {
                int ref_hint = frame_info->ref_order_hint[frame_info->ref_frame_idx[i]];
                if (get_relative_dist(&seq_header->order_hint_info, ref_hint, forward_hint) < 0) {
                    if (second_forward_idx < 0 ||
                        get_relative_dist(
                            &seq_header->order_hint_info, ref_hint, second_forward_hint) > 0) {
                        second_forward_idx  = i;
                        second_forward_hint = ref_hint;
                    }
                }
            }
            if (second_forward_idx < 0) {
                frame_info->skip_mode_params.skip_mode_allowed = 0;
            } else {
                frame_info->skip_mode_params.skip_mode_allowed = 1;
                frame_info->skip_mode_params.ref_frame_idx_0 =
                    LAST_FRAME + MIN(forward_idx, second_forward_idx);
                frame_info->skip_mode_params.ref_frame_idx_1 =
                    LAST_FRAME + MAX(forward_idx, second_forward_idx);
            }
        }
    }

    if (frame_info->skip_mode_params.skip_mode_allowed)
        frame_info->skip_mode_params.skip_mode_flag = dec_get_bits(bs, 1);
    else
        frame_info->skip_mode_params.skip_mode_flag = 0;
    PRINT_FRAME("skip_mode_present", frame_info->skip_mode_params.skip_mode_flag);
}

void load_grain_params(EbDecHandle *dec_handle_ptr, AomFilmGrain *grain_params,
                       int film_grain_params_ref_idx) {
    EbDecPicBuf *ref_buf = dec_handle_ptr->ref_frame_map[film_grain_params_ref_idx];
    assert(ref_buf != NULL);
    *grain_params = ref_buf->film_grain_params;
}

// Read film grain parameters
void read_film_grain_params(EbDecHandle *dec_handle, Bitstrm *bs, AomFilmGrain *grain_params) {
    SeqHeader *  seq_header = &dec_handle->seq_header;
    FrameHeader *frame_info = &dec_handle->frame_header;
    int          i, num_pos_luma, num_pos_chroma;

    if (!seq_header->film_grain_params_present ||
        (!frame_info->show_frame && !frame_info->showable_frame)) {
        memset(grain_params, 0, sizeof(*grain_params));
        return;
    }
    grain_params->apply_grain = dec_get_bits(bs, 1);
    PRINT_FRAME("apply_grain", grain_params->apply_grain);

    if (!grain_params->apply_grain) {
        memset(grain_params, 0, sizeof(*grain_params));
        return;
    }

    grain_params->random_seed = dec_get_bits(bs, 16);
    PRINT_FRAME("grain_seed", grain_params->random_seed);
    if (frame_info->frame_type == INTER_FRAME)
        grain_params->update_parameters = dec_get_bits(bs, 1);
    else
        grain_params->update_parameters = 1;
    PRINT_FRAME("update_parameters", grain_params->update_parameters);
    if (!grain_params->update_parameters) {
        int film_grain_params_ref_idx = dec_get_bits(bs, 3);
        PRINT_FRAME("film_grain_params_ref_idx", film_grain_params_ref_idx);
        uint16_t temp_grain_seed = grain_params->random_seed;
        load_grain_params(dec_handle, grain_params, film_grain_params_ref_idx);
        grain_params->random_seed = temp_grain_seed;
        return;
    }
    grain_params->num_y_points = dec_get_bits(bs, 4);
    assert(grain_params->num_y_points <= 14);
    PRINT_FRAME("num_y_points", grain_params->num_y_points);
    for (i = 0; i < grain_params->num_y_points; i++) {
        grain_params->scaling_points_y[i][0] = dec_get_bits(bs, 8);
        grain_params->scaling_points_y[i][1] = dec_get_bits(bs, 8);
        if (i > 0)
            assert(grain_params->scaling_points_y[i][0] > grain_params->scaling_points_y[i - 1][0]);
        PRINT_FRAME("scaling_points_y[i][0]", grain_params->scaling_points_y[i][0]);
        PRINT_FRAME("scaling_points_y[i][1]", grain_params->scaling_points_y[i][1]);
    }
    if (seq_header->color_config.mono_chrome)
        grain_params->chroma_scaling_from_luma = 0;
    else
        grain_params->chroma_scaling_from_luma = dec_get_bits(bs, 1);
    PRINT_FRAME("chroma_scaling_from_luma", grain_params->chroma_scaling_from_luma);

    if (seq_header->color_config.mono_chrome || grain_params->chroma_scaling_from_luma ||
        ((seq_header->color_config.subsampling_y == 1) &&
         (seq_header->color_config.subsampling_x == 1) && grain_params->num_y_points == 0)) {
        grain_params->num_cb_points = 0;
        grain_params->num_cr_points = 0;
    } else {
        grain_params->num_cb_points = dec_get_bits(bs, 4);
        PRINT_FRAME("num_cb_points", grain_params->num_cb_points);
        assert(grain_params->num_cb_points <= 10);
        for (i = 0; i < grain_params->num_cb_points; i++) {
            grain_params->scaling_points_cb[i][0] = dec_get_bits(bs, 8);
            grain_params->scaling_points_cb[i][1] = dec_get_bits(bs, 8);
            PRINT_FRAME("scaling_points_cb[i][0]", grain_params->scaling_points_cb[i][0]);
            PRINT_FRAME("scaling_points_cb[i][1]", grain_params->scaling_points_cb[i][1]);
            if (i > 0)
                assert(grain_params->scaling_points_cb[i][0] >
                       grain_params->scaling_points_cb[i - 1][0]);
        }
        grain_params->num_cr_points = dec_get_bits(bs, 4);
        PRINT_FRAME("num_cr_points", grain_params->num_cr_points);
        assert(grain_params->num_cr_points <= 14);
        for (i = 0; i < grain_params->num_cr_points; i++) {
            grain_params->scaling_points_cr[i][0] = dec_get_bits(bs, 8);
            grain_params->scaling_points_cr[i][1] = dec_get_bits(bs, 8);
            PRINT_FRAME("scaling_points_cr[i][0]", grain_params->scaling_points_cr[i][0]);
            PRINT_FRAME("scaling_points_cr[i][1]", grain_params->scaling_points_cr[i][1]);
            if (i > 0)
                assert(grain_params->scaling_points_cr[i][0] >
                       grain_params->scaling_points_cr[i - 1][0]);
        }
    }

    if ((seq_header->color_config.subsampling_x == 1) &&
        (seq_header->color_config.subsampling_y == 1) &&
        (((grain_params->num_cb_points == 0) && (grain_params->num_cr_points != 0)) ||
         ((grain_params->num_cb_points != 0) && (grain_params->num_cr_points == 0))))
        return; // EB_DecUnsupportedBitstream;

    grain_params->scaling_shift = dec_get_bits(bs, 2) + 8;
    grain_params->ar_coeff_lag  = dec_get_bits(bs, 2);
    PRINT_FRAME("scaling_shift", grain_params->grain_scale_shift);
    PRINT_FRAME("ar_coeff_lag", grain_params->ar_coeff_lag);

    num_pos_luma = 2 * grain_params->ar_coeff_lag * (grain_params->ar_coeff_lag + 1);
    if (grain_params->num_y_points) {
        num_pos_chroma = num_pos_luma + 1;
        for (i = 0; i < num_pos_luma; i++) {
            grain_params->ar_coeffs_y[i] = dec_get_bits(bs, 8) - 128;
            PRINT_FRAME("ar_coeffs_y[i]", grain_params->ar_coeffs_y[i]);
        }
    } else
        num_pos_chroma = num_pos_luma;
    if (grain_params->chroma_scaling_from_luma || grain_params->num_cb_points) {
        for (i = 0; i < num_pos_chroma; i++) {
            grain_params->ar_coeffs_cb[i] = dec_get_bits(bs, 8) - 128;
            PRINT_FRAME("ar_coeffs_cb[i]", grain_params->ar_coeffs_cb[i]);
        }
    }
    if (grain_params->chroma_scaling_from_luma || grain_params->num_cr_points) {
        for (i = 0; i < num_pos_chroma; i++) {
            grain_params->ar_coeffs_cr[i] = dec_get_bits(bs, 8) - 128;
            PRINT_FRAME("ar_coeffs_cr[i]", grain_params->ar_coeffs_cr[i]);
        }
    }
    grain_params->ar_coeff_shift    = dec_get_bits(bs, 2) + 6;
    grain_params->grain_scale_shift = dec_get_bits(bs, 2);
    PRINT_FRAME("ar_coeff_shift", grain_params->ar_coeff_shift);
    PRINT_FRAME("grain_scale_shift", grain_params->grain_scale_shift);
    if (grain_params->num_cb_points) {
        grain_params->cb_mult      = dec_get_bits(bs, 8);
        grain_params->cb_luma_mult = dec_get_bits(bs, 8);
        grain_params->cb_offset    = dec_get_bits(bs, 9);
        PRINT_FRAME("cb_mult", grain_params->cb_mult);
        PRINT_FRAME("cb_luma_mult", grain_params->cb_luma_mult);
        PRINT_FRAME("cb_offset", grain_params->cb_offset);
    }
    if (grain_params->num_cr_points) {
        grain_params->cr_mult      = dec_get_bits(bs, 8);
        grain_params->cr_luma_mult = dec_get_bits(bs, 8);
        grain_params->cr_offset    = dec_get_bits(bs, 9);
        PRINT_FRAME("cr_mult", grain_params->cr_mult);
        PRINT_FRAME("cr_luma_mult", grain_params->cr_luma_mult);
        PRINT_FRAME("cr_offset", grain_params->cr_offset);
    }
    grain_params->overlap_flag             = dec_get_bits(bs, 1);
    grain_params->clip_to_restricted_range = dec_get_bits(bs, 1);
    PRINT_FRAME("overlap_flag", grain_params->overlap_flag);
    PRINT_FRAME("clip_to_restricted_range", grain_params->clip_to_restricted_range);
}

int seg_feature_active_idx(SegmentationParams *seg_params, int segment_id,
                           SEG_LVL_FEATURES feature_id) {
    return seg_params->segmentation_enabled &&
           (seg_params->feature_enabled[segment_id][feature_id]);
}

int get_qindex(SegmentationParams *seg_params, int segment_id, int base_q_idx) {
    if (seg_feature_active_idx(seg_params, segment_id, SEG_LVL_ALT_Q)) {
        int data    = seg_params->feature_data[segment_id][SEG_LVL_ALT_Q];
        int q_index = base_q_idx + data;
        return clamp(q_index, 0, MAXQ);
    } else
        return base_q_idx;
}

EbErrorType reset_parse_ctx(FRAME_CONTEXT *frm_ctx, uint8_t base_qp) {
    EbErrorType return_error = EB_ErrorNone;

    svt_av1_default_coef_probs(frm_ctx, base_qp);
    init_mode_probs(frm_ctx);

    return return_error;
}

void setup_frame_sign_bias(EbDecHandle *dec_handle) {
    MvReferenceFrame ref_frame;
    for (ref_frame = LAST_FRAME; ref_frame <= ALTREF_FRAME; ++ref_frame) {
        const EbDecPicBuf *const buf = get_ref_frame_buf(dec_handle, ref_frame);
        if (dec_handle->seq_header.order_hint_info.enable_order_hint && buf != NULL) {
            const int ref_order_hint = buf->order_hint;
            dec_handle->frame_header.ref_frame_sign_bias[ref_frame] =
                (get_relative_dist(&dec_handle->seq_header.order_hint_info,
                                   ref_order_hint,
                                   (int)dec_handle->cur_pic_buf[0]->order_hint) <= 0)
                    ? 0
                    : 1;
        } else {
            dec_handle->frame_header.ref_frame_sign_bias[ref_frame] = 0;
        }
    }
}

void setup_past_independence(EbDecHandle *dec_handle_ptr, FrameHeader *frame_info) {
    int                ref, i, j;
    EbDecPicBuf *      cur_buf = dec_handle_ptr->cur_pic_buf[0];
    SegmentationParams seg     = dec_handle_ptr->frame_header.segmentation_params;

    SeqHeader *seq_header = &dec_handle_ptr->seq_header;
    int        size       = (seq_header->max_frame_width * seq_header->max_frame_height) >> 4;
    if (cur_buf->segment_maps) memset(cur_buf->segment_maps, 0, size);

    for (i = 0; i < MAX_SEGMENTS; i++)
        for (j = 0; j < SEG_LVL_MAX; j++) {
            seg.feature_data[i][j]    = 0;
            seg.feature_enabled[i][j] = 0;
        }
    UNUSED(seg);

    for (ref = LAST_FRAME; ref <= ALTREF_FRAME; ref++)
        cur_buf->global_motion[ref].gm_type = IDENTITY;

    frame_info->loop_filter_params.mode_ref_delta_enabled    = 1;
    frame_info->loop_filter_params.ref_deltas[INTRA_FRAME]   = 1;
    frame_info->loop_filter_params.ref_deltas[LAST_FRAME]    = 0;
    frame_info->loop_filter_params.ref_deltas[LAST2_FRAME]   = 0;
    frame_info->loop_filter_params.ref_deltas[LAST3_FRAME]   = 0;
    frame_info->loop_filter_params.ref_deltas[BWDREF_FRAME]  = 0;
    frame_info->loop_filter_params.ref_deltas[GOLDEN_FRAME]  = -1;
    frame_info->loop_filter_params.ref_deltas[ALTREF_FRAME]  = -1;
    frame_info->loop_filter_params.ref_deltas[ALTREF2_FRAME] = -1;

    frame_info->loop_filter_params.mode_deltas[0] = 0;
    frame_info->loop_filter_params.mode_deltas[1] = 0;
}

static INLINE EbErrorType reallocate_parse_context_memory(EbDecHandle *    dec_handle_ptr,
                                                          MasterParseCtxt *master_parse_ctx,
                                                          int              num_instances) {
    SeqHeader *seq_header = &dec_handle_ptr->seq_header;
    int32_t    num_mi_frame;

    master_parse_ctx->context_count = num_instances;

    int32_t num_mi_sb        = seq_header->sb_mi_size;
    int32_t sb_size_log2     = seq_header->sb_size_log2;
    int32_t sb_aligned_width = ALIGN_POWER_OF_TWO(seq_header->max_frame_width, sb_size_log2);
    int32_t sb_cols          = sb_aligned_width >> sb_size_log2;
    int8_t  num_planes       = seq_header->color_config.mono_chrome ? 1 : MAX_MB_PLANE;
    num_mi_frame             = sb_cols * num_mi_sb;
    int num_mi_64x64         = mi_size_wide[BLOCK_64X64];

    TilesInfo tiles_info = dec_handle_ptr->frame_header.tiles_info;
    int       num_tiles  = tiles_info.tile_cols * tiles_info.tile_rows;
    int32_t   num_ctx    = num_instances == 1 ? 1 : num_tiles;
    if (num_instances == 1) master_parse_ctx->context_count = num_tiles;

    /* TO-DO this memory will be freed at the end of decode.
       Can be optimized by reallocating the memory when
       the number of tiles changes within a sequence. */
    EB_MALLOC_DEC(
        ParseCtxt *, master_parse_ctx->tile_parse_ctxt, sizeof(ParseCtxt) * num_ctx, EB_N_PTR);

    EB_MALLOC_DEC(ParseAboveNbr4x4Ctxt *,
                  master_parse_ctx->parse_above_nbr4x4_ctxt,
                  sizeof(ParseAboveNbr4x4Ctxt) * num_ctx,
                  EB_N_PTR);
    EB_MALLOC_DEC(ParseLeftNbr4x4Ctxt *,
                  master_parse_ctx->parse_left_nbr4x4_ctxt,
                  sizeof(ParseLeftNbr4x4Ctxt) * num_ctx,
                  EB_N_PTR);
    int total_rows = num_instances == 1 ? 1 : tiles_info.tile_rows;
    int total_cols = num_instances == 1 ? 1 : tiles_info.tile_cols;
    for (int row = 0; row < total_rows; row++) {
        for (int col = 0; col < total_cols; col++) {
            int     instance = (row * total_cols) + col;
            int32_t num_mi_tile =
                tiles_info.tile_col_start_mi[col + 1] - tiles_info.tile_col_start_mi[col];
            int32_t num_mi_wide = num_instances == 1 ? num_mi_frame : num_mi_tile;
            num_mi_wide         = ALIGN_POWER_OF_TWO(num_mi_wide, sb_size_log2 - MI_SIZE_LOG2);
            ParseAboveNbr4x4Ctxt *above_ctx = &master_parse_ctx->parse_above_nbr4x4_ctxt[instance];
            ParseLeftNbr4x4Ctxt * left_ctx  = &master_parse_ctx->parse_left_nbr4x4_ctxt[instance];
            EB_MALLOC_DEC(
                uint8_t *, above_ctx->above_tx_wd, num_mi_wide * sizeof(uint8_t), EB_N_PTR);
            EB_MALLOC_DEC(
                uint8_t *, above_ctx->above_part_wd, num_mi_wide * sizeof(uint8_t), EB_N_PTR);
            EB_MALLOC_DEC(uint8_t *, left_ctx->left_tx_ht, num_mi_sb * sizeof(uint8_t), EB_N_PTR);
            EB_MALLOC_DEC(uint8_t *, left_ctx->left_part_ht, num_mi_sb * sizeof(uint8_t), EB_N_PTR);
            /* TODO : Optimize the size for Chroma */
            for (int i = 0; i < num_planes; i++) {
                EB_MALLOC_DEC(
                    uint8_t *, above_ctx->above_ctx[i], num_mi_wide * sizeof(uint8_t), EB_N_PTR);
                EB_MALLOC_DEC(uint16_t *,
                              above_ctx->above_palette_colors[i],
                              num_mi_64x64 * PALETTE_MAX_SIZE * sizeof(uint16_t),
                              EB_N_PTR);

                EB_MALLOC_DEC(
                    uint8_t *, left_ctx->left_ctx[i], num_mi_sb * sizeof(uint8_t), EB_N_PTR);
                EB_MALLOC_DEC(uint16_t *,
                              left_ctx->left_palette_colors[i],
                              num_mi_sb * PALETTE_MAX_SIZE * sizeof(uint16_t),
                              EB_N_PTR);
            }
            EB_MALLOC_DEC(
                int8_t *, above_ctx->above_comp_grp_idx, num_mi_wide * sizeof(int8_t), EB_N_PTR);
            EB_MALLOC_DEC(
                uint8_t *, above_ctx->above_seg_pred_ctx, num_mi_wide * sizeof(uint8_t), EB_N_PTR);
            EB_MALLOC_DEC(
                int8_t *, left_ctx->left_comp_grp_idx, num_mi_sb * sizeof(int8_t), EB_N_PTR);
            EB_MALLOC_DEC(
                uint8_t *, left_ctx->left_seg_pred_ctx, num_mi_sb * sizeof(uint8_t), EB_N_PTR);
        }
    }
    return EB_ErrorNone;
}

static INLINE EbErrorType reallocate_parse_tile_data(MasterParseCtxt *master_parse_ctx,
                                                     int              num_tiles) {
    master_parse_ctx->num_tiles = num_tiles;
    /* TO-DO this memory will be freed at the end of decode.
       Can be optimized by reallocating the memory when
       the number of tiles changes within a sequence. */
    EB_MALLOC_DEC(ParseTileData *,
                  master_parse_ctx->parse_tile_data,
                  sizeof(ParseTileData) * num_tiles,
                  EB_N_PTR);
    return EB_ErrorNone;
}

void set_prev_frame_info(EbDecHandle *dec_handle_ptr) {
    DecMtFrameData *dec_mt_frame_data =
        &dec_handle_ptr->master_frame_buf.cur_frame_bufs[0].dec_mt_frame_data;
    TilesInfo tiles_info = dec_handle_ptr->frame_header.tiles_info;
    dec_mt_frame_data->prev_frame_info.prev_max_frame_width =
        dec_handle_ptr->frame_header.frame_size.frame_width;
    dec_mt_frame_data->prev_frame_info.prev_max_frame_height =
        dec_handle_ptr->frame_header.frame_size.frame_height;
    dec_mt_frame_data->prev_frame_info.frame_header_read = EB_TRUE;
    dec_mt_frame_data->prev_frame_info.prev_sb_size =
        dec_handle_ptr->seq_header.sb_size;
    svt_memcpy(&dec_mt_frame_data->prev_frame_info.prev_tiles_info,
        &tiles_info, sizeof(TilesInfo));
}

static void realloc_parse_memory(EbDecHandle *dec_handle_ptr) {
    MasterParseCtxt *master_parse_ctx =
        (MasterParseCtxt *)dec_handle_ptr->pv_master_parse_ctxt;
    TilesInfo tiles_info = dec_handle_ptr->frame_header.tiles_info;
    int num_tiles = tiles_info.tile_cols * tiles_info.tile_rows;
    int num_instances = MIN((int32_t)dec_handle_ptr->dec_config.threads,
        num_tiles);
    if (dec_handle_ptr->dec_config.threads == 1) {
        /* For single thread case, allocate memory for one
           frame row above and one sb column for the left context. */
        reallocate_parse_context_memory(dec_handle_ptr,
            master_parse_ctx, 1);
    }
    else {
        reallocate_parse_context_memory(dec_handle_ptr,
            master_parse_ctx, num_instances);
    }
    if (num_tiles != master_parse_ctx->num_tiles)
        reallocate_parse_tile_data(master_parse_ctx, num_tiles);
}

static void check_mt_support(EbDecHandle *dec_handle_ptr) {
    TilesInfo       tiles_info = dec_handle_ptr->frame_header.tiles_info;

    DecMtFrameData *dec_mt_frame_data =
        &dec_handle_ptr->master_frame_buf.cur_frame_bufs[0].dec_mt_frame_data;

    if (dec_mt_frame_data->prev_frame_info.frame_header_read != EB_TRUE) {
        set_prev_frame_info(dec_handle_ptr);
        return;
    }

    EbBool do_realloc = EB_FALSE;
    if (dec_mt_frame_data->prev_frame_info.prev_max_frame_width !=
        dec_handle_ptr->frame_header.frame_size.frame_width ||
        dec_mt_frame_data->prev_frame_info.prev_max_frame_height !=
        dec_handle_ptr->frame_header.frame_size.frame_height)
    {
        do_realloc = EB_TRUE;
    }

    if (dec_mt_frame_data->prev_frame_info.prev_sb_size !=
        dec_handle_ptr->seq_header.sb_size)
    {
        do_realloc = EB_TRUE;
    }

    if (dec_mt_frame_data->prev_frame_info.prev_tiles_info.tile_cols !=
        tiles_info.tile_cols ||
        dec_mt_frame_data->prev_frame_info.prev_tiles_info.tile_rows !=
        tiles_info.tile_rows)
    {
        do_realloc = EB_TRUE;
    }

    for (int i = 0; i <= tiles_info.tile_cols; i++) {
        if (dec_mt_frame_data->prev_frame_info.
            prev_tiles_info.tile_col_start_mi[i] !=
            tiles_info.tile_col_start_mi[i])
        {
            do_realloc = EB_TRUE;
            break;
        }
    }

    for (int i = 0; i <= tiles_info.tile_rows; i++) {
        if (dec_mt_frame_data->prev_frame_info.prev_tiles_info.
            tile_row_start_mi[i] !=
            tiles_info.tile_row_start_mi[i])
        {
            do_realloc = EB_TRUE;
            break;
        }
    }

    if (do_realloc) {
        EbMemoryMapEntry *memory_entry = svt_dec_memory_map;
        EbMemoryMapEntry *previous_entry = NULL;
        if (memory_entry != memory_map_end_address) {
            while ((EbMemoryMapEntry *)memory_entry->prev_entry !=
                memory_map_end_address)
            {
                memory_entry = (EbMemoryMapEntry *)memory_entry->prev_entry;
            }
            previous_entry = memory_entry;
            memory_entry = (EbMemoryMapEntry *)memory_entry->prev_entry;
        }
        do {
            switch (memory_entry->ptr_type) {
            case EB_N_PTR: free(memory_entry->ptr); break;
            case EB_A_PTR:
#ifdef _WIN32
                _aligned_free(memory_entry->ptr);
#else
                free(memory_entry->ptr);
#endif
                break;
            case EB_SEMAPHORE: svt_destroy_semaphore(memory_entry->ptr); break;
            case EB_THREAD: svt_destroy_thread(memory_entry->ptr); break;
            case EB_MUTEX: svt_destroy_mutex(memory_entry->ptr); break;
            default: break;
            }
            EbMemoryMapEntry *tmp_memory_entry = memory_entry;
            memory_entry = (EbMemoryMapEntry *)tmp_memory_entry->prev_entry;
            free(tmp_memory_entry);
        } while (memory_entry != memory_map_start_address && memory_entry);
        if (previous_entry != NULL)
            previous_entry->prev_entry = memory_map_start_address;
        else
            svt_dec_memory_map = memory_map_start_address;
        dec_system_resource_init(dec_handle_ptr, &tiles_info);
        set_prev_frame_info(dec_handle_ptr);
        realloc_parse_memory(dec_handle_ptr);
    }
}

void read_uncompressed_header(Bitstrm *bs, EbDecHandle *dec_handle_ptr, ObuHeader *obu_header,
                              int num_planes) {
    SeqHeader *  seq_header = &dec_handle_ptr->seq_header;
    FrameHeader *frame_info = &dec_handle_ptr->frame_header;
    int          id_len = 0, all_frames, frame_is_intra = 0, frame_size_override_flag = 0;
    uint32_t     prev_frame_id = 0;

    if (seq_header->frame_id_numbers_present_flag) {
        id_len = seq_header->frame_id_length - 1 + seq_header->delta_frame_id_length - 2 + 3;
        assert(id_len <= 16);
    }
    all_frames = (1 << NUM_REF_FRAMES) - 1;
    if (seq_header->reduced_still_picture_header) {
        frame_info->show_existing_frame  = 0;
        frame_info->frame_type           = KEY_FRAME;
        frame_is_intra                   = 1;
        frame_info->show_frame           = 1;
        frame_info->showable_frame       = 0;
        frame_info->error_resilient_mode = 1;
    } else {
        frame_info->show_existing_frame = dec_get_bits(bs, 1);
        PRINT_FRAME("show_existing_frame", frame_info->show_existing_frame);
        if (frame_info->show_existing_frame) {
            int frame_to_show_map_idx = dec_get_bits(bs, 3);
            PRINT_FRAME("frame_to_show_map_idx", frame_to_show_map_idx);
            if (seq_header->decoder_model_info_present_flag &&
                !seq_header->timing_info.equal_picture_interval)
                temporal_point_info(bs, &seq_header->decoder_model_info, frame_info);
            frame_info->refresh_frame_flags = 0;
            if (seq_header->frame_id_numbers_present_flag) {
                uint32_t display_frame_id = dec_get_bits(bs, id_len);
                PRINT_FRAME("display_frame_id", display_frame_id);
                if (display_frame_id != frame_info->ref_frame_idx[frame_to_show_map_idx] &&
                    frame_info->ref_valid[frame_to_show_map_idx] == 1)
                    return; // EB_Corrupt_Frame;
            }

            dec_handle_ptr->cur_pic_buf[0] = dec_handle_ptr->ref_frame_map[frame_to_show_map_idx];
            frame_info->frame_type         = dec_handle_ptr->cur_pic_buf[0]->frame_type;

            if (frame_info->frame_type == KEY_FRAME) {
                frame_info->refresh_frame_flags = all_frames;
                frame_info->showable_frame      = 0;
            }

            if (seq_header->film_grain_params_present)
                load_grain_params(
                    dec_handle_ptr, &frame_info->film_grain_params, frame_to_show_map_idx);

            generate_next_ref_frame_map(dec_handle_ptr);

            frame_info->show_frame = 1;
            dec_handle_ptr->cur_pic_buf[0]->film_grain_params =
                dec_handle_ptr->frame_header.film_grain_params;
            dec_handle_ptr->show_existing_frame = frame_info->show_existing_frame;
            dec_handle_ptr->show_frame          = frame_info->show_frame;
            dec_handle_ptr->showable_frame      = frame_info->showable_frame;
            return;
        }

        frame_info->frame_type = dec_get_bits(bs, 2);

        frame_is_intra =
            (frame_info->frame_type == INTRA_ONLY_FRAME || frame_info->frame_type == KEY_FRAME);
        frame_info->show_frame = dec_get_bits(bs, 1);
        if (frame_info->show_frame && seq_header->decoder_model_info_present_flag &&
            !seq_header->timing_info.equal_picture_interval)
            temporal_point_info(bs, &seq_header->decoder_model_info, frame_info);
        if (frame_info->show_frame)
            frame_info->showable_frame = (frame_info->frame_type != KEY_FRAME);
        else
            frame_info->showable_frame = dec_get_bits(bs, 1);
        if (frame_info->frame_type == S_FRAME ||
            (frame_info->frame_type == KEY_FRAME && frame_info->show_frame))
            frame_info->error_resilient_mode = 1;
        else
            frame_info->error_resilient_mode = dec_get_bits(bs, 1);
    }
    PRINT_FRAME("frame_type", frame_info->frame_type);
    PRINT_FRAME("show_frame", frame_info->show_frame);
    PRINT_FRAME("showable_frame", frame_info->showable_frame);
    PRINT_FRAME("error_resilient_mode", frame_info->error_resilient_mode);
    if (frame_info->frame_type == KEY_FRAME && frame_info->show_frame) {
        for (int i = 0; i < NUM_REF_FRAMES; i++) {
            frame_info->ref_valid[i] = 0;
            // TODO: Need to differentate RefOrderHint and ref_order_hint
            frame_info->ref_order_hint[i] = 0;
        }
    }
    frame_info->disable_cdf_update = dec_get_bits(bs, 1);
    PRINT_FRAME("disable_cdf_update", frame_info->disable_cdf_update);
    if (seq_header->seq_force_screen_content_tools == SELECT_SCREEN_CONTENT_TOOLS)
        frame_info->allow_screen_content_tools = dec_get_bits(bs, 1);
    else
        frame_info->allow_screen_content_tools = seq_header->seq_force_screen_content_tools;
    PRINT_FRAME("allow_screen_content_tools", frame_info->allow_screen_content_tools);

    if (frame_info->allow_screen_content_tools) {
        if (seq_header->seq_force_integer_mv == SELECT_INTEGER_MV)
            frame_info->force_integer_mv = dec_get_bits(bs, 1);
        else
            frame_info->force_integer_mv = seq_header->seq_force_integer_mv;
    } else
        frame_info->force_integer_mv = 0;
    PRINT_FRAME("force_integer_mv", frame_info->force_integer_mv);

    if (frame_is_intra) frame_info->force_integer_mv = 1;
    int have_prev_frame_id = /*!pbi->decoding_first_frame && */
        !(frame_info->frame_type == KEY_FRAME && frame_info->show_frame);
    if (have_prev_frame_id) prev_frame_id = frame_info->current_frame_id;

    if (seq_header->frame_id_numbers_present_flag) {
        // int PrevFrameID = frame_info->current_frame_id;
        frame_info->current_frame_id = dec_get_bits(bs, id_len);
        PRINT_FRAME("current_frame_id", frame_info->current_frame_id);

        if (have_prev_frame_id) {
            int diff_frame_id = frame_info->current_frame_id > prev_frame_id
                ? frame_info->current_frame_id - prev_frame_id
                : (1 << id_len) + frame_info->current_frame_id - prev_frame_id;
            // Bitstream conformance
            if (frame_info->current_frame_id == prev_frame_id || diff_frame_id >= 1 << (id_len - 1))
                return; // EB_Corrupt_Frame;
        }

        //mark_ref_frames( id_len )
        uint32_t diff_len = seq_header->delta_frame_id_length;
        for (int i = 0; i < REF_FRAMES; i++) {
            if (frame_info->current_frame_id > (uint32_t)(1 << diff_len)) {
                if (frame_info->ref_frame_idx[i] > frame_info->current_frame_id ||
                    frame_info->ref_frame_idx[i] > (frame_info->current_frame_id - (1 - diff_len)))
                    frame_info->ref_valid[i] = 0;
            } else if (frame_info->ref_frame_idx[i] > frame_info->current_frame_id &&
                    frame_info->ref_frame_idx[i] <
                        (uint32_t)((1 << id_len) + frame_info->current_frame_id - (1 << diff_len)))
                    frame_info->ref_valid[i] = 0;
        }
    } else
        frame_info->current_frame_id = 0;
    if (frame_info->frame_type == S_FRAME)
        frame_size_override_flag = 1;
    else if (seq_header->reduced_still_picture_header)
        frame_size_override_flag = 0;
    else
        frame_size_override_flag = dec_get_bits(bs, 1);
    PRINT_FRAME("frame_size_override_flag", frame_size_override_flag);
    frame_info->order_hint = dec_get_bits(bs, seq_header->order_hint_info.order_hint_bits);
    PRINT_FRAME("order_hint", frame_info->order_hint);

    if (frame_is_intra || frame_info->error_resilient_mode)
        frame_info->primary_ref_frame = PRIMARY_REF_NONE;
    else {
        frame_info->primary_ref_frame = dec_get_bits(bs, PRIMARY_REF_BITS);
        PRINT_FRAME("primary_ref_frame", frame_info->primary_ref_frame)
    }
    if (seq_header->decoder_model_info_present_flag) {
        frame_info->buffer_removal_time_present_flag = dec_get_bits(bs, 1);
        PRINT_FRAME("buffer_removal_time_present_flag",
                    frame_info->buffer_removal_time_present_flag);
        if (frame_info->buffer_removal_time_present_flag) {
            for (int op_num = 0; op_num <= seq_header->operating_points_cnt_minus_1; op_num++) {
                if (seq_header->operating_point[op_num].decoder_model_present_for_this_op) {
                    uint16_t op_pt_idc = seq_header->operating_point[op_num].op_idc;
                    int      in_temporal_layer = (op_pt_idc >> obu_header->temporal_id) & 1;
                    int      in_spatial_layer  = (op_pt_idc >> (obu_header->spatial_id + 8)) & 1;
                    if (op_pt_idc == 0 || (in_temporal_layer && in_spatial_layer))
                        frame_info->buffer_removal_time[op_num] = dec_get_bits(
                            bs,
                            seq_header->decoder_model_info.buffer_removal_time_length_minus_1 + 1);
                } else
                    frame_info->buffer_removal_time[op_num] = 0;
                PRINT_FRAME("buffer_removal_time[op_num]", frame_info->buffer_removal_time[op_num]);
            }
        }
    }

    frame_info->allow_high_precision_mv = 0;
    frame_info->use_ref_frame_mvs       = 0;
    frame_info->allow_intrabc           = 0;
    if (frame_info->frame_type == S_FRAME ||
        (frame_info->frame_type == KEY_FRAME && frame_info->show_frame))
        frame_info->refresh_frame_flags = 0xFF;
    else
        frame_info->refresh_frame_flags = dec_get_bits(bs, 8);

    if (frame_info->frame_type == INTRA_ONLY_FRAME) assert(frame_info->refresh_frame_flags != 0xFF);

    PRINT_FRAME("refresh_frame_flags", frame_info->refresh_frame_flags);
    if (!frame_is_intra || (frame_info->refresh_frame_flags != 0xFF)) {
        if (frame_info->error_resilient_mode && seq_header->order_hint_info.enable_order_hint) {
            for (int i = 0; i < NUM_REF_FRAMES; i++) {
                int ref_order_hint = dec_get_bits(bs, seq_header->order_hint_info.order_hint_bits);
                PRINT_FRAME("ref_order_hint[i]", ref_order_hint);

                if (ref_order_hint != (int)frame_info->ref_order_hint[i])
                    frame_info->ref_valid[i] = 0;
            }
        }
    }
    if (frame_is_intra) {
        read_frame_size(bs, seq_header, frame_info, frame_size_override_flag);
        read_render_size(bs, frame_info);
        if (frame_info->allow_screen_content_tools && frame_info->frame_size.render_width) {
            if (frame_info->allow_screen_content_tools &&
                frame_info->frame_size.frame_width ==
                    frame_info->frame_size.superres_upscaled_width) {
                frame_info->allow_intrabc = dec_get_bits(bs, 1);
                PRINT_FRAME("allow_intrabc", frame_info->allow_intrabc);
            }
        }
        dec_handle_ptr->prev_frame = NULL;
    } else {
        int          frame_refs_short_signaling;
        if (!seq_header->order_hint_info.enable_order_hint)
            frame_refs_short_signaling = 0;
        else {
            frame_refs_short_signaling = dec_get_bits(bs, 1);
            PRINT_FRAME("frame_refs_short_signaling", frame_refs_short_signaling);
            if (frame_refs_short_signaling) {
                int last_frame_idx = dec_get_bits(bs, 3);
                int gold_frame_idx = dec_get_bits(bs, 3);
                PRINT_FRAME("last_frame_idx", last_frame_idx);
                PRINT_FRAME("gold_frame_idx", gold_frame_idx);
                svt_set_frame_refs(dec_handle_ptr, last_frame_idx, gold_frame_idx);
            }
        }

        //int DeltaFrameId;
        for (int i = 0; i < INTER_REFS_PER_FRAME; i++) {
            if (!frame_refs_short_signaling) {
                frame_info->ref_frame_idx[i] = dec_get_bits(bs, 3);
                PRINT_FRAME("ref_frame_idx", frame_info->ref_frame_idx[i]);
                dec_handle_ptr->remapped_ref_idx[i] = frame_info->ref_frame_idx[i];
            }
            int ref_frm_id = dec_handle_ptr->remapped_ref_idx[i];

            frame_info->ref_frame_sign_bias[LAST_FRAME + i] = 0;

            if (seq_header->frame_id_numbers_present_flag) {
                int delta_frame_id_length_minus_1 = dec_get_bits(bs,
                                                                 seq_header->delta_frame_id_length);
                PRINT_FRAME("delta_frame_id_length_minus_1", delta_frame_id_length_minus_1);
                uint32_t expected_frame_id = ((frame_info->current_frame_id + (1 << id_len) -
                                               (delta_frame_id_length_minus_1 + 1)) %
                                              (1 << id_len));
                if (expected_frame_id != frame_info->ref_frame_id[ref_frm_id]) {
                    assert(0);
                    return; // EB_Corrupt_Frame;
                }
            }
        }

        if (frame_size_override_flag && !frame_info->error_resilient_mode)
            frame_size_with_refs(bs, dec_handle_ptr, frame_size_override_flag);
        else {
            read_frame_size(bs, seq_header, frame_info, frame_size_override_flag);
            read_render_size(bs, frame_info);
        }
        if (frame_info->force_integer_mv)
            frame_info->allow_high_precision_mv = 0;
        else {
            frame_info->allow_high_precision_mv = dec_get_bits(bs, 1);
            PRINT_FRAME("allow_high_precision_mv", frame_info->allow_high_precision_mv);
        }
        read_interpolation_filter(bs, frame_info);
        frame_info->is_motion_mode_switchable = dec_get_bits(bs, 1);
        PRINT_FRAME("is_motion_mode_switchable", frame_info->is_motion_mode_switchable);
        dec_handle_ptr->prev_frame = get_primary_ref_frame_buf(dec_handle_ptr);
        if (frame_info->primary_ref_frame != PRIMARY_REF_NONE &&
            dec_handle_ptr->prev_frame == NULL)
        {
            SVT_LOG("Reference frame containing this frame's initial "
                "frame context is unavailable.");
            assert(0);
        }
        if (frame_info->error_resilient_mode || !seq_header->order_hint_info.enable_ref_frame_mvs)
            frame_info->use_ref_frame_mvs = 0;
        else
            frame_info->use_ref_frame_mvs = dec_get_bits(bs, 1);
        PRINT_FRAME("use_ref_frame_mvs", frame_info->use_ref_frame_mvs);

        for (int i = LAST_FRAME; i <= ALTREF_FRAME; ++i) {
            const EbDecPicBuf *const   ref_buf           = get_ref_frame_buf(dec_handle_ptr, i);
            struct ScaleFactors *const ref_scale_factors = get_ref_scale_factors(dec_handle_ptr, i);

            svt_av1_setup_scale_factors_for_frame(ref_scale_factors,
                                                  ref_buf->superres_upscaled_width,
                                                  ref_buf->frame_height,
                                                  frame_info->frame_size.frame_width,
                                                  frame_info->frame_size.frame_height);

            if ((!av1_is_valid_scale(ref_scale_factors))) {
                SVT_LOG("\n Reference frame has invalid dimensions \n");
                assert(0);
            }
        }
    }

    if (seq_header->color_config.subsampling_x == 1 && seq_header->color_config.subsampling_y == 1)
        dec_handle_ptr->dec_config.max_color_format = EB_YUV420;
    else if (seq_header->color_config.subsampling_x == 1 &&
             seq_header->color_config.subsampling_y == 0)
        dec_handle_ptr->dec_config.max_color_format = EB_YUV422;
    else if (seq_header->color_config.subsampling_x == 0 &&
             seq_header->color_config.subsampling_y == 0)
        dec_handle_ptr->dec_config.max_color_format = EB_YUV444;

    dec_handle_ptr->cur_pic_buf[0] =
        dec_pic_mgr_get_cur_pic(dec_handle_ptr);

    svt_setup_frame_buf_refs(dec_handle_ptr);
    /*Temporal MVs allocation */
    check_add_tplmv_buf(dec_handle_ptr);

    setup_frame_sign_bias(dec_handle_ptr);

    if (seq_header->reduced_still_picture_header || frame_info->disable_cdf_update)
        frame_info->disable_frame_end_update_cdf = 1;
    else {
        frame_info->disable_frame_end_update_cdf = dec_get_bits(bs, 1);
        PRINT_FRAME("disable_frame_end_update_cdf",
                    frame_info->disable_frame_end_update_cdf ? REFRESH_FRAME_CONTEXT_DISABLED
                                                             : REFRESH_FRAME_CONTEXT_BACKWARD);
    }

    if (frame_info->primary_ref_frame == PRIMARY_REF_NONE)
        setup_past_independence(dec_handle_ptr, frame_info);

    generate_next_ref_frame_map(dec_handle_ptr);

    read_tile_info(bs, &frame_info->tiles_info, seq_header, frame_info);
    read_quantization_params(
        bs, &frame_info->quantization_params, &seq_header->color_config, num_planes);
    read_segmentation_params(bs, dec_handle_ptr, frame_info);
    read_frame_delta_q_params(bs, frame_info);
    read_frame_delta_lf_params(bs, frame_info);
    setup_segmentation_dequant((DecModCtxt *)dec_handle_ptr->pv_dec_mod_ctxt);

    MasterParseCtxt *master_parse_ctx = (MasterParseCtxt *)dec_handle_ptr->pv_master_parse_ctxt;
    if (frame_info->primary_ref_frame == PRIMARY_REF_NONE)
        reset_parse_ctx(&master_parse_ctx->init_frm_ctx,
                        frame_info->quantization_params.base_q_idx);
    else
        /* Load CDF */
        master_parse_ctx->init_frm_ctx = dec_handle_ptr->prev_frame->final_frm_ctx;

    TilesInfo tiles_info = dec_handle_ptr->frame_header.tiles_info;

    if (dec_handle_ptr->dec_config.threads > 1) {
        /* Call System Resource Init only once */
        if (EB_FALSE == dec_handle_ptr->start_thread_process) {
            dec_system_resource_init(dec_handle_ptr, &tiles_info);
            dec_handle_ptr->start_thread_process = EB_TRUE;
        }
        check_mt_support(dec_handle_ptr);
    }

    int       num_tiles = tiles_info.tile_cols * tiles_info.tile_rows;
    int       num_instances = num_tiles;

    if(dec_handle_ptr->dec_config.threads != 1)
        num_instances = MIN((int32_t)dec_handle_ptr->dec_config.threads, num_tiles);

    if (num_instances != master_parse_ctx->context_count)
        realloc_parse_memory(dec_handle_ptr);

    frame_info->coded_lossless = 1;
    for (int i = 0; i < MAX_SEGMENTS; ++i) {
        int qindex = get_qindex(
            &frame_info->segmentation_params, i, frame_info->quantization_params.base_q_idx);
        frame_info->quantization_params.qindex[i] = qindex;
        frame_info->lossless_array[i] =
            qindex == 0 && frame_info->quantization_params.delta_q_dc[AOM_PLANE_Y] == 0 &&
            frame_info->quantization_params.delta_q_ac[AOM_PLANE_U] == 0 &&
            frame_info->quantization_params.delta_q_dc[AOM_PLANE_U] == 0 &&
            frame_info->quantization_params.delta_q_ac[AOM_PLANE_V] == 0 &&
            frame_info->quantization_params.delta_q_dc[AOM_PLANE_V] == 0;
        if (!frame_info->lossless_array[i]) frame_info->coded_lossless = 0;
        if (frame_info->quantization_params.using_qmatrix) {
            if (frame_info->lossless_array[i]) {
                frame_info->segmentation_params.seg_qm_level[0][i] = 15;
                frame_info->segmentation_params.seg_qm_level[1][i] = 15;
                frame_info->segmentation_params.seg_qm_level[2][i] = 15;
            } else {
                frame_info->segmentation_params.seg_qm_level[0][i] =
                    frame_info->quantization_params.qm[AOM_PLANE_Y];
                frame_info->segmentation_params.seg_qm_level[1][i] =
                    frame_info->quantization_params.qm[AOM_PLANE_U];
                frame_info->segmentation_params.seg_qm_level[2][i] =
                    frame_info->quantization_params.qm[AOM_PLANE_V];
            }
        }
    }

    if (frame_info->coded_lossless == 1) assert(frame_info->delta_q_params.delta_q_present == 0);

    frame_info->all_lossless =
        frame_info->coded_lossless &&
        (frame_info->frame_size.frame_width == frame_info->frame_size.superres_upscaled_width);
    read_loop_filter_params(bs, dec_handle_ptr, num_planes);
    read_frame_cdef_params(bs, frame_info, seq_header, num_planes);
    read_lr_params(bs, frame_info, seq_header, num_planes);
    read_tx_mode(bs, frame_info);

    frame_info->reference_mode =
        read_frame_reference_mode(bs, frame_is_intra) ? REFERENCE_MODE_SELECT : SINGLE_REFERENCE;
    PRINT_FRAME("reference_mode",
                frame_info->reference_mode ? REFERENCE_MODE_SELECT : SINGLE_REFERENCE);
    read_skip_mode_params(bs, frame_info, frame_is_intra, seq_header, frame_info->reference_mode);

    if (frame_is_intra || frame_info->error_resilient_mode || !seq_header->enable_warped_motion)
        frame_info->allow_warped_motion = 0;
    else
        frame_info->allow_warped_motion = dec_get_bits(bs, 1);
    frame_info->reduced_tx_set = dec_get_bits(bs, 1);
    PRINT_FRAME("allow_warped_motion", frame_info->allow_warped_motion);
    PRINT_FRAME("reduced_tx_set", frame_info->reduced_tx_set);
    read_global_motion_params(bs, dec_handle_ptr, frame_info, frame_is_intra);

    read_film_grain_params(dec_handle_ptr, bs, &frame_info->film_grain_params);

    dec_handle_ptr->cur_pic_buf[0]->film_grain_params =
        dec_handle_ptr->frame_header.film_grain_params;

    dec_handle_ptr->show_existing_frame = frame_info->show_existing_frame;
    dec_handle_ptr->show_frame          = frame_info->show_frame;
    dec_handle_ptr->showable_frame      = frame_info->showable_frame;

    /* TODO: Should be moved to caller */
    if (dec_handle_ptr->dec_config.threads == 1) {
        if (!frame_info->show_existing_frame)
            svt_setup_motion_field(dec_handle_ptr, NULL);
    }
}

EbErrorType read_frame_header_obu(Bitstrm *bs, EbDecHandle *dec_handle_ptr, ObuHeader *obu_header,
                                  int trailing_bit) {
    EbErrorType status = EB_ErrorNone;

    int      num_planes = av1_num_planes(&dec_handle_ptr->seq_header.color_config);
    uint32_t start_position, end_position, header_bytes;

    start_position = get_position(bs);
    read_uncompressed_header(bs, dec_handle_ptr, obu_header, num_planes);

    if (allow_intrabc(dec_handle_ptr)) {
        svt_av1_setup_scale_factors_for_frame(&dec_handle_ptr->sf_identity,
                                              dec_handle_ptr->cur_pic_buf[0]->frame_width,
                                              dec_handle_ptr->cur_pic_buf[0]->frame_height,
                                              dec_handle_ptr->cur_pic_buf[0]->frame_width,
                                              dec_handle_ptr->cur_pic_buf[0]->frame_height);
    }

    if (trailing_bit) {
        status = av1_check_trailing_bits(bs);
        if (status != EB_ErrorNone) return status;
    }

    byte_alignment(bs);

    end_position = get_position(bs);
    header_bytes = (end_position - start_position) / 8;
    obu_header->payload_size -= header_bytes;

    return status;
}

// Read Tile group information
EbErrorType read_tile_group_obu(Bitstrm *bs, EbDecHandle *dec_handle_ptr, TilesInfo *tiles_info,
                                ObuHeader *obu_header, int *is_last_tg) {
    EbErrorType status = EB_ErrorNone;

    MasterParseCtxt *master_parse_ctxt = (MasterParseCtxt *)dec_handle_ptr->pv_master_parse_ctxt;

    FrameHeader *frame_header = &dec_handle_ptr->frame_header;

    DecMtFrameData *dec_mt_frame_data =
        &dec_handle_ptr->master_frame_buf.cur_frame_bufs[0].dec_mt_frame_data;

    int      num_tiles, tg_start, tg_end, tile_start_and_end_present_flag = 0;
    uint32_t start_position, end_position, header_bytes;
    num_tiles = tiles_info->tile_cols * tiles_info->tile_rows;

    int32_t sb_size_log2 = dec_handle_ptr->seq_header.sb_size_log2;
    int32_t sb_aligned_width =
        ALIGN_POWER_OF_TWO(frame_header->frame_size.frame_width, sb_size_log2);
    int32_t sb_aligned_height =
        ALIGN_POWER_OF_TWO(frame_header->frame_size.frame_height, sb_size_log2);
    int32_t sb_cols            = sb_aligned_width >> sb_size_log2;
    int32_t sb_rows            = sb_aligned_height >> sb_size_log2;
    dec_mt_frame_data->sb_cols = sb_cols;
    dec_mt_frame_data->sb_rows = sb_rows;

    start_position = get_position(bs);
    if (num_tiles > 1) {
        tile_start_and_end_present_flag = dec_get_bits(bs, 1);
        PRINT_FRAME("tile_start_and_end_present_flag", tile_start_and_end_present_flag);
    }

    if (obu_header->obu_type == OBU_FRAME) assert(tile_start_and_end_present_flag == 0);
    if (num_tiles == 1 || !tile_start_and_end_present_flag) {
        tg_start = 0;
        tg_end   = num_tiles - 1;
    } else {
        uint8_t tile_bits = tiles_info->tile_cols_log2 + tiles_info->tile_rows_log2;
        tg_start  = dec_get_bits(bs, tile_bits);
        tg_end    = dec_get_bits(bs, tile_bits);
    }
    assert(tg_end >= tg_start);
    PRINT_FRAME("tg_start", tg_start);
    PRINT_FRAME("tg_end", tg_end);

    *is_last_tg = ((tg_end + 1) == num_tiles);

    byte_alignment(bs);
    end_position = get_position(bs);
    header_bytes = (end_position - start_position) / 8;
    obu_header->payload_size -= header_bytes;

    dec_handle_ptr->cm.mi_cols       = dec_handle_ptr->frame_header.mi_cols;
    dec_handle_ptr->cm.mi_rows       = dec_handle_ptr->frame_header.mi_rows;
    dec_handle_ptr->cm.mi_stride     = dec_handle_ptr->frame_header.mi_stride;
    dec_handle_ptr->cm.bit_depth     = dec_handle_ptr->seq_header.color_config.bit_depth;
    dec_handle_ptr->cm.subsampling_x = dec_handle_ptr->seq_header.color_config.subsampling_x;
    dec_handle_ptr->cm.subsampling_y = dec_handle_ptr->seq_header.color_config.subsampling_y;
    dec_handle_ptr->cm.frm_size      = dec_handle_ptr->frame_header.frame_size;
    dec_handle_ptr->cm.tiles_info    = dec_handle_ptr->frame_header.tiles_info;
    dec_handle_ptr->is_lf_enabled =
        (!dec_handle_ptr->frame_header.allow_intrabc &&
         (dec_handle_ptr->frame_header.loop_filter_params.filter_level[0] ||
          dec_handle_ptr->frame_header.loop_filter_params.filter_level[1]));

    uint32_t num_threads = dec_handle_ptr->dec_config.threads;
    int      is_mt       = num_threads != 1;

    /* PPF flags derivation */
    EbBool no_ibc = !dec_handle_ptr->frame_header.allow_intrabc;
    /* LF */
    EbBool do_lf_flag =
        no_ibc && (dec_handle_ptr->frame_header.loop_filter_params.filter_level[0] ||
            dec_handle_ptr->frame_header.loop_filter_params.filter_level[1]);
    /* CDEF */
    EbBool do_cdef = no_ibc && (!frame_header->coded_lossless &&
        (frame_header->cdef_params.cdef_bits ||
            frame_header->cdef_params.cdef_y_strength[0] ||
            frame_header->cdef_params.cdef_uv_strength[0]));

    EbBool do_upscale = no_ibc &&
        !av1_superres_unscaled(&dec_handle_ptr->frame_header.frame_size);
    /* LR */
    //EbBool opt_lr = !do_cdef && !do_upscale;
    LrParams *lr_param = dec_handle_ptr->frame_header.lr_params;
    EbBool    do_lr = no_ibc &&
        (lr_param[AOM_PLANE_Y].frame_restoration_type != RESTORE_NONE ||
        lr_param[AOM_PLANE_U].frame_restoration_type != RESTORE_NONE ||
        lr_param[AOM_PLANE_V].frame_restoration_type != RESTORE_NONE);

    /* Set Parse Jobs */
    if (is_mt) {
        svt_av1_scan_tiles(dec_handle_ptr, tiles_info, obu_header, bs, tg_start, tg_end);
        if ((tg_end + 1) != num_tiles) return 0;
        {
            int32_t tiles_ctr;

            for (tiles_ctr = 0; tiles_ctr < num_tiles; tiles_ctr++) {
                uint32_t *sb_recon_completed_in_row, *sb_recon_row_started;
                uint32_t *sb_recon_row_parsed;
                uint32_t  tile_num_sb_rows;

                sb_recon_row_parsed =
                    dec_mt_frame_data->parse_recon_tile_info_array[tiles_ctr].sb_recon_row_parsed;
                sb_recon_completed_in_row =
                    dec_mt_frame_data->parse_recon_tile_info_array[tiles_ctr]
                        .sb_recon_completed_in_row;
                sb_recon_row_started =
                    dec_mt_frame_data->parse_recon_tile_info_array[tiles_ctr].sb_recon_row_started;
                tile_num_sb_rows =
                    dec_mt_frame_data->parse_recon_tile_info_array[tiles_ctr].tile_num_sb_rows;

                dec_mt_frame_data->parse_recon_tile_info_array[tiles_ctr].sb_row_to_process = 0;

                memset(sb_recon_row_parsed, 0, tile_num_sb_rows * sizeof(uint32_t));
                memset(sb_recon_completed_in_row, 0, tile_num_sb_rows * sizeof(uint32_t));
                memset(sb_recon_row_started, 0, tile_num_sb_rows * sizeof(uint32_t));
            }
        }

        const int mvs_rows    = (dec_handle_ptr->frame_header.mi_rows + 1) >> 1; //8x8 unit level
        const int sb_mvs_rows = (mvs_rows + 7) >> 3; //64x64 unit level
        dec_mt_frame_data->motion_proj_info.num_motion_proj_rows       = sb_mvs_rows;
        dec_mt_frame_data->motion_proj_info.motion_proj_row_to_process = 0;
        dec_mt_frame_data->motion_proj_info.motion_proj_init_done      = EB_FALSE;
        dec_mt_frame_data->num_threads_header                          = 0;

        svt_block_on_mutex(dec_mt_frame_data->temp_mutex);
        dec_mt_frame_data->start_motion_proj = EB_TRUE;
        svt_release_mutex(dec_mt_frame_data->temp_mutex);
        svt_post_semaphore(dec_handle_ptr->thread_semaphore);
        for (uint32_t lib_thrd = 0; lib_thrd < num_threads - 1; lib_thrd++)
            svt_post_semaphore(dec_handle_ptr->thread_ctxt_pa[lib_thrd].thread_semaphore);

        svt_setup_motion_field(dec_handle_ptr, NULL);

        svt_av1_queue_parse_jobs(dec_handle_ptr, tiles_info);

        svt_block_on_mutex(dec_mt_frame_data->temp_mutex);
        dec_mt_frame_data->start_parse_frame = EB_TRUE;

        dec_mt_frame_data->num_threads_cdefed = 0;
        dec_mt_frame_data->num_threads_lred   = 0;

        svt_release_mutex(dec_mt_frame_data->temp_mutex);
        svt_post_semaphore(dec_handle_ptr->thread_semaphore);
        for (uint32_t lib_thrd = 0; lib_thrd < num_threads - 1; lib_thrd++)
            svt_post_semaphore(dec_handle_ptr->thread_ctxt_pa[lib_thrd].thread_semaphore);

        svt_av1_queue_lf_jobs(dec_handle_ptr);
        svt_av1_queue_cdef_jobs(dec_handle_ptr);
        svt_block_on_mutex(dec_mt_frame_data->temp_mutex);

        dec_mt_frame_data->start_lf_frame = EB_TRUE;
        /*ToDo : Post outside mutex lock */
        svt_post_semaphore(dec_handle_ptr->thread_semaphore);
        for (uint32_t lib_thrd = 0; lib_thrd < num_threads - 1; lib_thrd++)
            svt_post_semaphore(dec_handle_ptr->thread_ctxt_pa[lib_thrd].thread_semaphore);
        dec_mt_frame_data->start_cdef_frame = EB_TRUE;
        svt_post_semaphore(dec_handle_ptr->thread_semaphore);
        for (uint32_t lib_thrd = 0; lib_thrd < num_threads - 1; lib_thrd++)
            svt_post_semaphore(dec_handle_ptr->thread_ctxt_pa[lib_thrd].thread_semaphore);
        svt_release_mutex(dec_mt_frame_data->temp_mutex);

        if(!do_upscale) svt_av1_queue_lr_jobs(dec_handle_ptr);

        parse_frame_tiles(dec_handle_ptr, 0);

        decode_frame_tiles(dec_handle_ptr, NULL);
    } else {
        //TO-DO assign to appropriate tile_parse_ctxt
        ParseCtxt *parse_ctxt               = &master_parse_ctxt->tile_parse_ctxt[0];
        parse_ctxt->seq_header              = &dec_handle_ptr->seq_header;
        parse_ctxt->frame_header            = &dec_handle_ptr->frame_header;
        parse_ctxt->parse_above_nbr4x4_ctxt = &master_parse_ctxt->parse_above_nbr4x4_ctxt[0];
        parse_ctxt->parse_left_nbr4x4_ctxt  = &master_parse_ctxt->parse_left_nbr4x4_ctxt[0];

        for (int tile_num = tg_start; tile_num <= tg_end; tile_num++) {
            size_t tile_size;
            if (tile_num == tg_end)
                tile_size = obu_header->payload_size;
            else {
                tile_size = dec_get_bits_le(bs, tiles_info->tile_size_bytes) + 1;
                obu_header->payload_size -= (tiles_info->tile_size_bytes + tile_size);
            }

            ParseTileData *parse_tile_data      = master_parse_ctxt->parse_tile_data;
            parse_tile_data[tile_num].data      = get_bitsteam_buf(bs);
            parse_tile_data[tile_num].data_end  = bs->buf_max;
            parse_tile_data[tile_num].tile_size = tile_size;

            start_parse_tile(dec_handle_ptr, parse_ctxt, tiles_info, tile_num, is_mt);
            dec_bits_init(bs, (get_bitsteam_buf(bs) + tile_size), obu_header->payload_size);
        }
    }

    if ((tg_end + 1) != num_tiles) return 0;

    if (is_mt) {
        dec_av1_loop_filter_frame_mt(dec_handle_ptr,
                                     dec_handle_ptr->cur_pic_buf[0]->ps_pic_buf,
                                     dec_handle_ptr->pv_lf_ctxt,
                                     AOM_PLANE_Y,
                                     MAX_MB_PLANE,
                                     NULL);
    } else {
        dec_av1_loop_filter_frame(dec_handle_ptr,
                                  dec_handle_ptr->cur_pic_buf[0]->ps_pic_buf,
                                  dec_handle_ptr->pv_lf_ctxt,
                                  AOM_PLANE_Y,
                                  MAX_MB_PLANE,
                                  is_mt,
                                  do_lf_flag);
    }

    if (!is_mt && do_lr) dec_av1_loop_restoration_save_boundary_lines(dec_handle_ptr, 0);

    if (is_mt) {
        svt_cdef_frame_mt(dec_handle_ptr, NULL);
    } else
        svt_cdef_frame(dec_handle_ptr, do_cdef);

    av1_superres_upscale(&dec_handle_ptr->cm,
                         &dec_handle_ptr->frame_header,
                         &dec_handle_ptr->seq_header,
                         dec_handle_ptr->cur_pic_buf[0]->ps_pic_buf,
                         do_upscale);

    if (do_upscale)
        dec_handle_ptr->cm.frm_size.frame_width =
            dec_handle_ptr->frame_header.frame_size.frame_width;

    if (do_lr && (!is_mt || do_upscale))
        dec_av1_loop_restoration_save_boundary_lines(dec_handle_ptr, 1);

    if (is_mt) {
        if (do_upscale) svt_av1_queue_lr_jobs(dec_handle_ptr);
        dec_handle_ptr->master_frame_buf.cur_frame_bufs[0].dec_mt_frame_data.start_lr_frame =
            EB_TRUE;
        svt_post_semaphore(dec_handle_ptr->thread_semaphore);
        for (uint32_t lib_thrd = 0; lib_thrd < num_threads - 1; lib_thrd++)
            svt_post_semaphore(dec_handle_ptr->thread_ctxt_pa[lib_thrd].thread_semaphore);
        dec_av1_loop_restoration_filter_frame_mt(dec_handle_ptr, NULL);
    } else
        dec_av1_loop_restoration_filter_frame(dec_handle_ptr, 0, /*opt_lr*/ do_lr);

    /* Save CDF */
    if (frame_header->disable_frame_end_update_cdf)
        dec_handle_ptr->cur_pic_buf[0]->final_frm_ctx = master_parse_ctxt->init_frm_ctx;

    if (!is_mt) { pad_pic(dec_handle_ptr); }

    return status;
}

// Decode all OBUs in a Frame
EbErrorType decode_multiple_obu(EbDecHandle *dec_handle_ptr, uint8_t **data, size_t data_size,
                                uint32_t is_annexb) {
    Bitstrm   bs;
    EbErrorType status = EB_ErrorNone;
    ObuHeader   obu_header;
    int         frame_decoding_finished = 0;

#if ENABLE_ENTROPY_TRACE
    enable_dump = 1;
#if FRAME_LEVEL_TRACE
    if (enable_dump) {
        char str[1000];
        sprintf(str, "SVT_fr_%d.txt", dec_handle_ptr->dec_cnt);
        if (temp_fp == NULL) temp_fp = fopen(str, "w");
    }
#else
    if (temp_fp == NULL) temp_fp = fopen("SVT.txt", "w");
#endif
#endif

    while (!frame_decoding_finished) {
        size_t payload_size = 0, length_size = 0;

        /* Decoder memory init if not done */
        if (0 == dec_handle_ptr->mem_init_done && 1 == dec_handle_ptr->seq_header_done)
            status = dec_mem_init(dec_handle_ptr);
        if (status != EB_ErrorNone) return status;

        dec_bits_init(&bs, *data, data_size);

        if (is_annexb) {
            // read the size of OBU
            status = read_obu_size(&bs, data_size, &obu_header.payload_size, &length_size);
            if (status != EB_ErrorNone) return status;

            *data += length_size;
            data_size -= length_size;
            length_size = 0;
        }

        status = read_obu_header_size(&bs, &obu_header, data_size, &length_size);
        if (status != EB_ErrorNone) return status;

        if (is_annexb) obu_header.payload_size -= obu_header.size;

        payload_size = obu_header.payload_size;

        *data += (obu_header.size + length_size);
        data_size -= (obu_header.size + length_size);

        if (data_size < payload_size) return EB_Corrupt_Frame;

        dec_bits_init(&bs, *data, payload_size);

        switch (obu_header.obu_type) {
        case OBU_TEMPORAL_DELIMITER:
            PRINT_NAME("**************OBU_TEMPORAL_DELIMITER*******************");
            read_temporal_delimitor_obu(&dec_handle_ptr->seen_frame_header);
            break;

        case OBU_SEQUENCE_HEADER: {
            PRINT_NAME("**************OBU_SEQUENCE_HEADER*******************")
            BlockSize prev_sb_size          = dec_handle_ptr->seq_header.sb_size;
            uint16_t  prev_max_frame_width  = dec_handle_ptr->seq_header.max_frame_width;
            uint16_t  prev_max_frame_height = dec_handle_ptr->seq_header.max_frame_height;

            status = read_sequence_header_obu(&bs, &dec_handle_ptr->seq_header);
            if (status != EB_ErrorNone) return status;
            if (dec_handle_ptr->seq_header.color_config.bit_depth == EB_TWELVE_BIT)
                dec_init_intra_predictors_12b_internal();
            dec_handle_ptr->seq_header_done = 1;
            if (prev_sb_size != dec_handle_ptr->seq_header.sb_size ||
                prev_max_frame_width != dec_handle_ptr->seq_header.max_frame_width ||
                prev_max_frame_height != dec_handle_ptr->seq_header.max_frame_height) {
                dec_handle_ptr->mem_init_done = 0;
            }
            break;
        }
        case OBU_FRAME_HEADER:
        case OBU_REDUNDANT_FRAME_HEADER:
        case OBU_FRAME:
            if (obu_header.obu_type == OBU_FRAME) {
                PRINT_NAME("**************OBU_FRAME*******************");
                dec_handle_ptr->show_existing_frame = 0;
            } else if (obu_header.obu_type == OBU_FRAME_HEADER) {
                PRINT_NAME("**************OBU_FRAME_HEADER*******************");
                assert(dec_handle_ptr->seen_frame_header == 0);
            } else {
                PRINT_NAME("**************OBU_REDUNDANT_FRAME_HEADER*******************");
                assert(dec_handle_ptr->seen_frame_header == 1);
            }

            if (!dec_handle_ptr->seen_frame_header) {
                dec_handle_ptr->seen_frame_header = 1;
                status                            = read_frame_header_obu(
                    &bs, dec_handle_ptr, &obu_header, obu_header.obu_type != OBU_FRAME);
            }
            /*else {
                 For OBU_REDUNDANT_FRAME_HEADER, previous frame_header is taken from dec_handle_ptr->frame_header
                //frame_header_copy(); TODO()
            }*/

            if (obu_header.obu_type != OBU_FRAME) break; // For OBU_TILE_GROUP comes under OBU_FRAME
            goto TITLE_GROUP;

        case OBU_TILE_GROUP:
        TITLE_GROUP:
            PRINT_NAME("**************OBU_TILE_GROUP*******************");
            if (!dec_handle_ptr->seen_frame_header) return EB_Corrupt_Frame;
            status = read_tile_group_obu(&bs,
                                         dec_handle_ptr,
                                         &dec_handle_ptr->frame_header.tiles_info,
                                         &obu_header,
                                         &frame_decoding_finished);
            if (status != EB_ErrorNone) return status;
            if (frame_decoding_finished) dec_handle_ptr->seen_frame_header = 0;
            break;

        default: PRINT_NAME("**************UNKNOWN OBU*******************"); break;
        }

        *data += payload_size;
        data_size -= payload_size;
        if (!data_size) frame_decoding_finished = 1;
    }

#if ENABLE_ENTROPY_TRACE
#if FRAME_LEVEL_TRACE
    if (enable_dump) {
        fclose(temp_fp);
        temp_fp = NULL;
    }
#endif
#endif

    return status;
}

EB_API EbErrorType svt_get_sequence_info(const uint8_t *obu_data, size_t size,
                                         SeqHeader *sequence_info) {
    if (obu_data == NULL || size == 0 || sequence_info == NULL) return EB_ErrorBadParameter;
    const uint8_t *frame_buf = obu_data;
    size_t         frame_sz  = size;
    EbErrorType    status    = EB_ErrorNone;
    do {
        Bitstrm bs;
        dec_bits_init(&bs, frame_buf, frame_sz);

        ObuHeader ou;
        memset(&ou, 0, sizeof(ou));
        size_t length_size = 0;
        status             = read_obu_header_size(&bs, &ou, frame_sz, &length_size);
        if (status != EB_ErrorNone) return status;

        frame_buf += ou.size + length_size;
        frame_sz -= (uint32_t)(ou.size + length_size);

        if (ou.obu_type == OBU_SEQUENCE_HEADER) {
            // check the ou type and parse sequence header
            status = read_sequence_header_obu(&bs, sequence_info);
            if (status == EB_ErrorNone) return status;
        }

        frame_buf += ou.payload_size;
        frame_sz -= ou.payload_size;
    } while (status == EB_ErrorNone && frame_sz > 0);
    return EB_ErrorUndefined;
}
