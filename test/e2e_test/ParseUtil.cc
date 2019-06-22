/*
 * Copyright(c) 2019 Netflix, Inc.
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

/******************************************************************************
 * @file ParseUtil.cc
 *
 * @brief Impelmentation of sequence header parser.
 *
 ******************************************************************************/

#include "EbDefinitions.h"
#include "ParseUtil.h"
#include "gtest/gtest.h"

namespace svt_av1_e2e_tools {

typedef struct aom_timing {
    uint32_t num_units_in_display_tick;
    uint32_t time_scale;
    int equal_picture_interval;
    uint32_t num_ticks_per_picture;
} aom_timing_info_t;

typedef struct aom_dec_model_info {
    uint32_t num_units_in_decoding_tick;
    int encoder_decoder_buffer_delay_length;
    int buffer_removal_delay_length;
    int frame_presentation_delay_length;
} aom_dec_model_info_t;

typedef struct aom_dec_model_op_parameters {
    int decoder_model_present_for_this_op;
    int64_t bitrate;
    int64_t buffer_size;
    int decoder_buffer_delay;
    int encoder_buffer_delay;
    int low_delay_mode_flag;
    int initial_display_delay_present_for_this_op;
    int initial_display_delay;
} aom_dec_model_op_parameters_t;

typedef struct aom_op_timing_info_t {
    int64_t buffer_removal_delay;
} aom_op_timing_info_t;

typedef struct SequenceHeader {
    int error;  // indicate what kind of error on parsing the headers
    int ready;  // 1 - valid sequence header;
                // 0 - invalid sequence header

    int still_picture;                 // Video is a single frame still picture
    int reduced_still_picture_header;  // Use reduced header for still picture

    /* timing info */
    int timing_info_present_flag;
    int decoder_model_info_present_flag;
    aom_timing_info_t timing_info;
    aom_dec_model_info_t buffer_model;

    /* operating points */
    aom_dec_model_op_parameters_t op_params[MAX_NUM_OPERATING_POINTS + 1];
    aom_op_timing_info_t op_frame_timing[MAX_NUM_OPERATING_POINTS + 1];
    int operating_points_cnt_minus_1;
    int operating_point_idc[MAX_NUM_OPERATING_POINTS];
    int initial_display_delay_present_flag;
    uint8_t tier[MAX_NUM_OPERATING_POINTS];  // seq_tier in the spec. One bit: 0

    /* profile and levels */
    BitstreamProfile profile;
    uint8_t seq_level_idx;
    unsigned int number_temporal_layers;
    unsigned int number_spatial_layers;
    BitstreamLevel level[MAX_NUM_OPERATING_POINTS];

    /* resolution and superblock size */
    int frame_width_bits;
    int frame_height_bits;
    int max_frame_width;
    int max_frame_height;
    BlockSize sb_size;  // Size of the superblock used for this frame
    int mib_size;       // Size of the superblock in units of MI blocks
    int mib_size_log2;  // Log 2 of above.

    /* frame id */
    int frame_id_numbers_present_flag;
    int frame_id_length;
    int delta_frame_id_length;

    /* coding tools */
    int order_hint_bits;
    int force_screen_content_tools;  // 0 - force off
                                     // 1 - force on
                                     // 2 - adaptive
    int force_integer_mv;          // 0 - Not to force. MV can be in 1/4 or 1/8
                                   // 1 - force to integer
                                   // 2 - adaptive
    int enable_filter_intra;       // enables/disables filterintra
    int enable_intra_edge_filter;  // enables/disables corner/edge/upsampling
    int enable_interintra_compound;  // enables/disables interintra_compound
    int enable_masked_compound;      // enables/disables masked compound
    int enable_dual_filter;          // 0 - disable dual interpolation filter
                                     // 1 - enable vert/horiz filter selection
    int enable_order_hint;     // 0 - disable order hint, and related tools
                               // jnt_comp, ref_frame_mvs, frame_sign_bias
                               // if 0, enable_jnt_comp and
                               // enable_ref_frame_mvs must be set zs 0.
    int enable_jnt_comp;       // 0 - disable joint compound modes
                               // 1 - enable it
    int enable_ref_frame_mvs;  // 0 - disable ref frame mvs
                               // 1 - enable it
    int enable_warped_motion;  // 0 - disable warped motion for sequence
                               // 1 - enable it for the sequence
    int enable_superres;  // 0 - Disable superres for the sequence, and disable
                          //     transmitting per-frame superres enabled flag.
                          // 1 - Enable superres for the sequence, and also
                          //     enable per-frame flag to denote if superres is
                          //     enabled for that frame.
    int enable_cdef;      // To turn on/off CDEF
    int enable_restoration;  // To turn on/off loop restoration
    int film_grain_params_present;

    /* color config */
    int monochrome;  // Monochorme video
    int bit_depth;
    int color_range;
    int subsampling_x;
    int subsampling_y;
    int separate_uv_delta_q;
    aom_color_primaries_t color_primaries;
    aom_transfer_characteristics_t transfer_characteristics;
    aom_matrix_coefficients_t matrix_coefficients;
    aom_chroma_sample_position_t chroma_sample_position;
} SequenceHeader;

typedef void (*aom_rb_error_handler)(void *data);

struct aom_read_bit_buffer {
    const uint8_t *bit_buffer;
    const uint8_t *bit_buffer_end;
    uint32_t bit_offset;

    void *error_handler_data;
    aom_rb_error_handler error_handler;
};

// ivf reader
#define IVF_FRAME_HDR_SZ (4 + 8) /* 4 byte size + 8 byte timestamp */
#define IVF_FILE_HDR_SZ 32
#define RAW_FRAME_HDR_SZ sizeof(uint32_t)

struct AvxRational {
    int numerator;
    int denominator;
};

struct AvxInputContext {
    const char *filename;
    FILE *file;
    uint32_t width;
    uint32_t height;
    uint32_t fourcc;
    struct AvxRational framerate;
};

typedef int64_t aom_codec_pts_t;

typedef enum ATTRIBUTE_PACKED {
    OBU_SEQUENCE_HEADER = 1,
    OBU_TEMPORAL_DELIMITER = 2,
    OBU_FRAME_HEADER = 3,
    OBU_TILE_GROUP = 4,
    OBU_METADATA = 5,
    OBU_FRAME = 6,
    OBU_REDUNDANT_FRAME_HEADER = 7,
    OBU_PADDING = 15,
} obuType;

/* obu reader */
typedef struct {
    size_t size;  // Size (1 or 2 bytes) of the OBU header (including the
                  // optional OBU extension header) in the bitstream.
    obuType type;
    int has_size_field;
    int has_extension;
    // The following fields come from the OBU extension header and therefore are
    // only used if has_extension is true.
    int temporal_layer_id;
    int spatial_layer_id;
} ObuHeader;

#define OBU_HEADER_SIZE 1
#define OBU_EXTENSION_SIZE 1
#define OBU_MAX_LENGTH_FIELD_SIZE 8

// bit reader
static size_t aom_rb_bytes_read(struct aom_read_bit_buffer *rb) {
    return (rb->bit_offset + 7) >> 3;
}

static int aom_rb_read_bit(struct aom_read_bit_buffer *rb) {
    const uint32_t off = rb->bit_offset;
    const uint32_t p = off >> 3;
    const int q = 7 - (int)(off & 0x7);
    if (rb->bit_buffer + p < rb->bit_buffer_end) {
        const int bit = (rb->bit_buffer[p] >> q) & 1;
        rb->bit_offset = off + 1;
        return bit;
    } else {
        if (rb->error_handler)
            rb->error_handler(rb->error_handler_data);
        return 0;
    }
}

static int aom_rb_read_literal(struct aom_read_bit_buffer *rb, int bits) {
    int value = 0, bit;
    for (bit = bits - 1; bit >= 0; bit--)
        value |= aom_rb_read_bit(rb) << bit;
    return value;
}

static uint32_t aom_rb_read_unsigned_literal(struct aom_read_bit_buffer *rb,
                                             int bits) {
    uint32_t value = 0;
    int bit;
    for (bit = bits - 1; bit >= 0; bit--)
        value |= (uint32_t)aom_rb_read_bit(rb) << bit;
    return value;
}

static int aom_rb_read_inv_signed_literal(struct aom_read_bit_buffer *rb,
                                          int bits) {
    const int nbits = sizeof(unsigned) * 8 - bits - 1;
    const unsigned value = (unsigned)aom_rb_read_literal(rb, bits + 1) << nbits;
    return ((int)value) >> nbits;
}

static uint32_t aom_rb_read_uvlc(struct aom_read_bit_buffer *rb) {
    int leading_zeros = 0;
    while (!aom_rb_read_bit(rb))
        ++leading_zeros;
    // Maximum 32 bits.
    if (leading_zeros >= 32)
        return UINT32_MAX;
    const uint32_t base = (1u << leading_zeros) - 1;
    const uint32_t value = aom_rb_read_literal(rb, leading_zeros);
    return base + value;
}

static const size_t kMaximumLeb128Size = 8;
static const uint8_t kLeb128ByteMask = 0x7f;  // Binary: 01111111

// uleb decoder
int aom_uleb_decode(const uint8_t *buffer, size_t available, uint64_t *value,
                    size_t *length) {
    if (buffer && value) {
        *value = 0;
        for (size_t i = 0; i < kMaximumLeb128Size && i < available; ++i) {
            const uint8_t decoded_byte = *(buffer + i) & kLeb128ByteMask;
            *value |= ((uint64_t)decoded_byte) << (i * 7);
            if ((*(buffer + i) >> 7) == 0) {
                if (length) {
                    *length = i + 1;
                }

                // Fail on values larger than 32-bits to ensure consistent
                // behavior on 32 and 64 bit targets: value is typically used to
                // determine buffer allocation size.
                if (*value > UINT32_MAX)
                    return -1;

                return 0;
            }
        }
    }

    // If we get here, either the buffer/value pointers were invalid,
    // or we ran over the available space
    return -1;
}

// ivf reader
static const char *IVF_SIGNATURE = "DKIF";
static void fix_framerate(int *num, int *den) {
    if (*den <= 0 || *den >= 1000000000 || *num <= 0 || *num >= 1000) {
        // framerate seems to be invalid, just default to 30fps.
        *num = 30;
        *den = 1;
    }
}

static unsigned int mem_get_le16(const void *vmem) {
    const uint8_t *mem = (const uint8_t *)vmem;
    return mem[1] << 8 | mem[0];
}

static unsigned int mem_get_le32(const void *vmem) {
    const uint8_t *mem = (const uint8_t *)vmem;
    return mem[3] << 24 | mem[2] << 16 | mem[1] << 8 | mem[0];
}

int file_is_ivf(struct AvxInputContext *input_ctx) {
    char raw_hdr[32];
    int is_ivf = 0;

    if (fread(raw_hdr, 1, 32, input_ctx->file) == 32) {
        if (memcmp(IVF_SIGNATURE, raw_hdr, 4) == 0) {
            is_ivf = 1;

            if (mem_get_le16(raw_hdr + 4) != 0) {
                fprintf(stderr,
                        "Error: Unrecognized IVF version! This file may not"
                        " decode properly.");
            }

            input_ctx->fourcc = mem_get_le32(raw_hdr + 8);
            input_ctx->width = mem_get_le16(raw_hdr + 12);
            input_ctx->height = mem_get_le16(raw_hdr + 14);
            input_ctx->framerate.numerator = mem_get_le32(raw_hdr + 16);
            input_ctx->framerate.denominator = mem_get_le32(raw_hdr + 20);
            fix_framerate(&input_ctx->framerate.numerator,
                          &input_ctx->framerate.denominator);
        }
    }

    if (!is_ivf) {
        rewind(input_ctx->file);
    }
#if 0
  if (!is_ivf) {
    rewind(input_ctx->file);
    input_ctx->detect.buf_read = 0;
  } else {
    input_ctx->detect.position = 4;
  }
#endif
    return is_ivf;
}

int ivf_read_frame(FILE *infile, uint8_t **buffer, size_t *bytes_read,
                   size_t *buffer_size, aom_codec_pts_t *pts) {
    char raw_header[IVF_FRAME_HDR_SZ] = {0};
    size_t frame_size = 0;

    if (fread(raw_header, IVF_FRAME_HDR_SZ, 1, infile) != 1) {
        if (!feof(infile)) {
            printf("Failed to read frame size");
            return -1;
        }
    } else {
        frame_size = mem_get_le32(raw_header);

        if (frame_size > 256 * 1024 * 1024) {
            printf("Read invalid frame size (%u)", (unsigned int)frame_size);
            frame_size = 0;
        }

        if (frame_size > *buffer_size) {
            uint8_t *new_buffer =
                reinterpret_cast<uint8_t *>(realloc(*buffer, 2 * frame_size));

            if (new_buffer) {
                *buffer = new_buffer;
                *buffer_size = 2 * frame_size;
            } else {
                printf("Failed to allocate compressed data buffer");
                frame_size = 0;
            }
        }

        if (pts) {
            *pts = mem_get_le32(&raw_header[4]);
            *pts += ((aom_codec_pts_t)mem_get_le32(&raw_header[8]) << 32);
        }
    }

    if (!feof(infile)) {
        if (fread(*buffer, 1, frame_size, infile) != frame_size) {
            printf("Failed to read full frame");
            return 1;
        }

        *bytes_read = frame_size;
        return 0;
    }

    return 1;
}

// obu reader
// Returns 1 when OBU type is valid, and 0 otherwise.
static int valid_obu_type(int obu_type) {
    int valid_type = 0;
    switch (obu_type) {
    case OBU_SEQUENCE_HEADER:
    case OBU_TEMPORAL_DELIMITER:
    case OBU_FRAME_HEADER:
    case OBU_TILE_GROUP:
    case OBU_METADATA:
    case OBU_FRAME:
    case OBU_REDUNDANT_FRAME_HEADER:
        /*case OBU_TILE_LIST:*/
    case OBU_PADDING: valid_type = 1; break;
    default: break;
    }
    return valid_type;
}

AomCodecErr read_obu_header(struct aom_read_bit_buffer *rb, int is_annexb,
                            ObuHeader *header) {
    if (!rb || !header)
        return AOM_CODEC_INVALID_PARAM;

    const ptrdiff_t bit_buffer_byte_length =
        rb->bit_buffer_end - rb->bit_buffer;
    if (bit_buffer_byte_length < 1)
        return AOM_CODEC_CORRUPT_FRAME;

    header->size = 1;

    if (aom_rb_read_bit(rb) != 0) {
        // Forbidden bit. Must not be set.
        return AOM_CODEC_CORRUPT_FRAME;
    }

    header->type = (obuType)aom_rb_read_literal(rb, 4);

    if (!valid_obu_type(header->type))
        return AOM_CODEC_CORRUPT_FRAME;

    header->has_extension = aom_rb_read_bit(rb);
    header->has_size_field = aom_rb_read_bit(rb);

    if (!header->has_size_field && !is_annexb) {
        // section 5 obu streams must have obu_size field set.
        return AOM_CODEC_UNSUP_BITSTREAM;
    }

    if (aom_rb_read_bit(rb) != 0) {
        // obu_reserved_1bit must be set to 0.
        return AOM_CODEC_CORRUPT_FRAME;
    }

    if (header->has_extension) {
        if (bit_buffer_byte_length == 1)
            return AOM_CODEC_CORRUPT_FRAME;

        header->size += 1;
        header->temporal_layer_id = aom_rb_read_literal(rb, 3);
        header->spatial_layer_id = aom_rb_read_literal(rb, 2);
        if (aom_rb_read_literal(rb, 3) != 0) {
            // extension_header_reserved_3bits must be set to 0.
            return AOM_CODEC_CORRUPT_FRAME;
        }
    }

    return AOM_CODEC_OK;
}

// sequence header parser
static int read_bitstream_level(BitstreamLevel *bl,
                                struct aom_read_bit_buffer *rb) {
    const uint8_t seq_level_idx = aom_rb_read_literal(rb, LEVEL_BITS);
    if (!is_valid_seq_level_idx(seq_level_idx))
        return 0;
    bl->major = (seq_level_idx >> LEVEL_MINOR_BITS) + LEVEL_MAJOR_MIN;
    bl->minor = seq_level_idx & ((1 << LEVEL_MINOR_BITS) - 1);
    return 1;
}
/* Max Bitrates for levels of Main Tier in kbps. Bitrate in main_kbps [31] */
/* is a dummy value. The decoder model is not applicable for level 31. */
static int32_t main_kbps[1 << LEVEL_BITS] = {
    1500,  3000,   0,      0,      6000,  10000, 0,     0,
    12000, 20000,  0,      0,      30000, 40000, 60000, 60000,
    60000, 100000, 160000, 160000, 0,     0,     0,     0,
    0,     0,      0,      0,      0,     0,     0,     (1 << 26)};

/* Max Bitrates for levels of High Tier in kbps. Bitrate in high_kbps [31] */
/* is a dummy value. The decoder model is not applicable for level 31. */
static int32_t high_kbps[1 << LEVEL_BITS] = {
    0,      0,      0,      0,      0,      0,      0,      0,
    30000,  50000,  0,      0,      100000, 160000, 240000, 240000,
    240000, 480000, 800000, 800000, 0,      0,      0,      0,
    0,      0,      0,      0,      0,      0,      0,      (1 << 26)};

/* BitrateProfileFactor */
static int bitrate_profile_factor[1 << PROFILE_BITS] = {1, 2, 3, 0, 0, 0, 0, 0};

static int64_t max_level_bitrate(BitstreamProfile seq_profile,
                                 int seq_level_idx, int seq_tier) {
    int64_t bitrate;

    if (seq_tier) {
        bitrate =
            high_kbps[seq_level_idx] * bitrate_profile_factor[seq_profile];
    } else {
        bitrate =
            main_kbps[seq_level_idx] * bitrate_profile_factor[seq_profile];
    }

    return bitrate * 1000;
}

// Reads the high_bitdepth and twelve_bit fields in color_config() and sets
// cm->bit_depth based on the values of those fields and cm->profile. Reports
// errors by calling rb->error_handler() or aom_internal_error().
static void av1_read_bitdepth(SequenceHeader *cm,
                              struct aom_read_bit_buffer *rb) {
    const int high_bitdepth = aom_rb_read_bit(rb);
    if (cm->profile == PROFILE_2 && high_bitdepth) {
        const int twelve_bit = aom_rb_read_bit(rb);
        cm->bit_depth = twelve_bit ? AOM_BITS_12 : AOM_BITS_10;
    } else if (cm->profile <= PROFILE_2) {
        cm->bit_depth = high_bitdepth ? AOM_BITS_10 : AOM_BITS_8;
    } else {
        cm->error = AOM_CODEC_UNSUP_BITSTREAM;
        printf("Unsupported profile/bit-depth combination");
    }
}

static AomCodecErr aom_get_num_layers_from_operating_point_idc(
    int operating_point_idc, unsigned int *number_spatial_layers,
    unsigned int *number_temporal_layers) {
    // derive number of spatial/temporal layers from operating_point_idc

    if (!number_spatial_layers || !number_temporal_layers)
        return AOM_CODEC_INVALID_PARAM;

    if (operating_point_idc == 0) {
        *number_temporal_layers = 1;
        *number_spatial_layers = 1;
    } else {
        *number_spatial_layers = 0;
        *number_temporal_layers = 0;
        for (int j = 0; j < MAX_NUM_SPATIAL_LAYERS; j++) {
            *number_spatial_layers +=
                (operating_point_idc >> (j + MAX_NUM_TEMPORAL_LAYERS)) & 0x1;
        }
        for (int j = 0; j < MAX_NUM_TEMPORAL_LAYERS; j++) {
            *number_temporal_layers += (operating_point_idc >> j) & 0x1;
        }
    }

    return AOM_CODEC_OK;
}

static void av1_read_timing_info_header(SequenceHeader *cm,
                                        struct aom_read_bit_buffer *rb) {
    cm->timing_info.num_units_in_display_tick = aom_rb_read_unsigned_literal(
        rb, 32);  // Number of units in a display tick
    cm->timing_info.time_scale =
        aom_rb_read_unsigned_literal(rb, 32);  // Time scale
    if (cm->timing_info.num_units_in_display_tick == 0 ||
        cm->timing_info.time_scale == 0) {
        printf(
            "num_units_in_display_tick and time_scale must be greater than 0.");
        cm->error = AOM_CODEC_UNSUP_BITSTREAM;
    }
    cm->timing_info.equal_picture_interval =
        aom_rb_read_bit(rb);  // Equal picture interval bit
    if (cm->timing_info.equal_picture_interval) {
        cm->timing_info.num_ticks_per_picture =
            aom_rb_read_uvlc(rb) + 1;  // ticks per picture
        if (cm->timing_info.num_ticks_per_picture == 0) {
            cm->error = AOM_CODEC_UNSUP_BITSTREAM;
            printf("num_ticks_per_picture_minus_1 cannot be (1 << 32) - 1.");
        }
    }
}

static void av1_read_decoder_model_info(SequenceHeader *cm,
                                        struct aom_read_bit_buffer *rb) {
    cm->buffer_model.encoder_decoder_buffer_delay_length =
        aom_rb_read_literal(rb, 5) + 1;
    cm->buffer_model.num_units_in_decoding_tick = aom_rb_read_unsigned_literal(
        rb, 32);  // Number of units in a decoding tick
    cm->buffer_model.buffer_removal_delay_length =
        aom_rb_read_literal(rb, 5) + 1;
    cm->buffer_model.frame_presentation_delay_length =
        aom_rb_read_literal(rb, 5) + 1;
}

static void av1_read_op_parameters_info(SequenceHeader *const cm,
                                        struct aom_read_bit_buffer *rb,
                                        int op_num) {
    // The cm->op_params array has MAX_NUM_OPERATING_POINTS + 1 elements.
    if (op_num > MAX_NUM_OPERATING_POINTS) {
        cm->error = AOM_CODEC_UNSUP_BITSTREAM;
        printf("AV1 does not support %d decoder model operating points",
               op_num + 1);
    }

    cm->op_params[op_num].decoder_buffer_delay = aom_rb_read_literal(
        rb, cm->buffer_model.encoder_decoder_buffer_delay_length);

    cm->op_params[op_num].encoder_buffer_delay = aom_rb_read_literal(
        rb, cm->buffer_model.encoder_decoder_buffer_delay_length);

    cm->op_params[op_num].low_delay_mode_flag = aom_rb_read_bit(rb);
}

static void av1_read_color_config(struct aom_read_bit_buffer *rb,
                                  SequenceHeader *seq_params) {
    av1_read_bitdepth(seq_params, rb);

    // monochrome bit (not needed for PROFILE_1)
    const int is_monochrome =
        seq_params->profile != PROFILE_1 ? aom_rb_read_bit(rb) : 0;
    seq_params->monochrome = is_monochrome;
    int color_description_present_flag = aom_rb_read_bit(rb);
    if (color_description_present_flag) {
        seq_params->color_primaries =
            (aom_color_primaries_t)aom_rb_read_literal(rb, 8);
        seq_params->transfer_characteristics =
            (aom_transfer_characteristics_t)aom_rb_read_literal(rb, 8);
        seq_params->matrix_coefficients =
            (aom_matrix_coefficients_t)aom_rb_read_literal(rb, 8);
    } else {
        seq_params->color_primaries = AOM_CICP_CP_UNSPECIFIED;
        seq_params->transfer_characteristics = AOM_CICP_TC_UNSPECIFIED;
        seq_params->matrix_coefficients = AOM_CICP_MC_UNSPECIFIED;
    }
    if (is_monochrome) {
        // [16,235] (including xvycc) vs [0,255] range
        seq_params->color_range = aom_rb_read_bit(rb);
        seq_params->subsampling_y = seq_params->subsampling_x = 1;
        seq_params->chroma_sample_position = AOM_CSP_UNKNOWN;
        seq_params->separate_uv_delta_q = 0;
        return;
    }
    if (seq_params->color_primaries == AOM_CICP_CP_BT_709 &&
        seq_params->transfer_characteristics == AOM_CICP_TC_SRGB &&
        seq_params->matrix_coefficients ==
            AOM_CICP_MC_IDENTITY) {  // it would be better
                                     // to remove this
                                     // dependency too
        seq_params->subsampling_y = seq_params->subsampling_x = 0;
        seq_params->color_range = 1;  // assume full color-range
        if (!(seq_params->profile == PROFILE_1 ||
              (seq_params->profile == PROFILE_2 &&
               seq_params->bit_depth == AOM_BITS_12))) {
            seq_params->error = AOM_CODEC_UNSUP_BITSTREAM,
            printf("sRGB colorspace not compatible with specified profile");
        }
    } else {
        // [16,235] (including xvycc) vs [0,255] range
        seq_params->color_range = aom_rb_read_bit(rb);
        if (seq_params->profile == PROFILE_0) {
            // 420 only
            seq_params->subsampling_x = seq_params->subsampling_y = 1;
        } else if (seq_params->profile == PROFILE_1) {
            // 444 only
            seq_params->subsampling_x = seq_params->subsampling_y = 0;
        } else {
            assert(seq_params->profile == PROFILE_2);
            if (seq_params->bit_depth == AOM_BITS_12) {
                seq_params->subsampling_x = aom_rb_read_bit(rb);
                if (seq_params->subsampling_x)
                    seq_params->subsampling_y =
                        aom_rb_read_bit(rb);  // 422 or 420
                else
                    seq_params->subsampling_y = 0;  // 444
            } else {
                // 422
                seq_params->subsampling_x = 1;
                seq_params->subsampling_y = 0;
            }
        }
        if (seq_params->matrix_coefficients == AOM_CICP_MC_IDENTITY &&
            (seq_params->subsampling_x || seq_params->subsampling_y)) {
            seq_params->error = AOM_CODEC_UNSUP_BITSTREAM;
            printf(
                "Identity CICP Matrix incompatible with non 4:4:4 color "
                "sampling");
        }
        if (seq_params->subsampling_x && seq_params->subsampling_y) {
            seq_params->chroma_sample_position =
                (aom_chroma_sample_position_t)aom_rb_read_literal(rb, 2);
        }
    }
    seq_params->separate_uv_delta_q = aom_rb_read_bit(rb);
}

// Checks that the remaining bits start with a 1 and ends with 0s.
// It consumes an additional byte, if already byte aligned before the check.
static int av1_check_trailing_bits(SequenceHeader *seq_params,
                                   struct aom_read_bit_buffer *rb) {
    // bit_offset is set to 0 (mod 8) when the reader is already byte aligned
    int bits_before_alignment = 8 - rb->bit_offset % 8;
    int trailing = aom_rb_read_literal(rb, bits_before_alignment);
    if (trailing != (1 << (bits_before_alignment - 1))) {
        seq_params->error = AOM_CODEC_CORRUPT_FRAME;
        return -1;
    }
    return 0;
}

static void av1_read_sequence_header(struct aom_read_bit_buffer *rb,
                                     SequenceHeader *seq_params) {
    const int num_bits_width = aom_rb_read_literal(rb, 4) + 1;
    const int num_bits_height = aom_rb_read_literal(rb, 4) + 1;
    const int max_frame_width = aom_rb_read_literal(rb, num_bits_width) + 1;
    const int max_frame_height = aom_rb_read_literal(rb, num_bits_height) + 1;

    seq_params->frame_width_bits = num_bits_width;
    seq_params->frame_height_bits = num_bits_height;
    seq_params->max_frame_width = max_frame_width;
    seq_params->max_frame_height = max_frame_height;

    if (seq_params->reduced_still_picture_header) {
        seq_params->frame_id_numbers_present_flag = 0;
    } else {
        seq_params->frame_id_numbers_present_flag = aom_rb_read_bit(rb);
    }
    if (seq_params->frame_id_numbers_present_flag) {
        // We must always have delta_frame_id_length < frame_id_length,
        // in order for a frame to be referenced with a unique delta.
        // Avoid wasting bits by using a coding that enforces this restriction.
        seq_params->delta_frame_id_length = aom_rb_read_literal(rb, 4) + 2;
        seq_params->frame_id_length =
            aom_rb_read_literal(rb, 3) + seq_params->delta_frame_id_length + 1;
        if (seq_params->frame_id_length > 16)
            seq_params->error = AOM_CODEC_CORRUPT_FRAME;
        printf("Invalid frame_id_length");
    }

    seq_params->sb_size = aom_rb_read_bit(rb) ? BLOCK_128X128 : BLOCK_64X64;
    seq_params->mib_size = mi_size_wide[seq_params->sb_size];
    seq_params->mib_size_log2 = mi_size_wide_log2[seq_params->sb_size];

    seq_params->enable_filter_intra = aom_rb_read_bit(rb);
    seq_params->enable_intra_edge_filter = aom_rb_read_bit(rb);

    if (seq_params->reduced_still_picture_header) {
        seq_params->enable_interintra_compound = 0;
        seq_params->enable_masked_compound = 0;
        seq_params->enable_warped_motion = 0;
        seq_params->enable_dual_filter = 0;
        seq_params->enable_order_hint = 0;
        seq_params->enable_jnt_comp = 0;
        seq_params->enable_ref_frame_mvs = 0;
        seq_params->force_screen_content_tools =
            2;                             // SELECT_SCREEN_CONTENT_TOOLS
        seq_params->force_integer_mv = 2;  // SELECT_INTEGER_MV
        seq_params->order_hint_bits = 0;
    } else {
        seq_params->enable_interintra_compound = aom_rb_read_bit(rb);
        seq_params->enable_masked_compound = aom_rb_read_bit(rb);
        seq_params->enable_warped_motion = aom_rb_read_bit(rb);
        seq_params->enable_dual_filter = aom_rb_read_bit(rb);

        seq_params->enable_order_hint = aom_rb_read_bit(rb);
        seq_params->enable_jnt_comp =
            seq_params->enable_order_hint ? aom_rb_read_bit(rb) : 0;
        seq_params->enable_ref_frame_mvs =
            seq_params->enable_order_hint ? aom_rb_read_bit(rb) : 0;

        if (aom_rb_read_bit(rb)) {
            seq_params->force_screen_content_tools =
                2;  // SELECT_SCREEN_CONTENT_TOOLS
        } else {
            seq_params->force_screen_content_tools = aom_rb_read_bit(rb);
        }

        if (seq_params->force_screen_content_tools > 0) {
            if (aom_rb_read_bit(rb)) {
                seq_params->force_integer_mv = 2;  // SELECT_INTEGER_MV
            } else {
                seq_params->force_integer_mv = aom_rb_read_bit(rb);
            }
        } else {
            seq_params->force_integer_mv = 2;  // SELECT_INTEGER_MV
        }
        seq_params->order_hint_bits =
            seq_params->enable_order_hint ? aom_rb_read_literal(rb, 3) + 1 : 0;
    }

    seq_params->enable_superres = aom_rb_read_bit(rb);
    seq_params->enable_cdef = aom_rb_read_bit(rb);
    seq_params->enable_restoration = aom_rb_read_bit(rb);
}

static INLINE uint8_t major_minor_to_seq_level_idx(BitstreamLevel bl) {
    assert(bl.major >= LEVEL_MAJOR_MIN && bl.major <= LEVEL_MAJOR_MAX);
    // assert(bl.minor >= LEVEL_MINOR_MIN && bl.minor <= LEVEL_MINOR_MAX);
    return ((bl.major - LEVEL_MAJOR_MIN) << LEVEL_MINOR_BITS) + bl.minor;
}

// parse sequence header;
// On success, sets pbi->sequence_header_ready to 1 and returns the number of
// bytes read from 'rb'.
// On failure, sets seq_params.error and returns 0.
uint32_t read_sequence_header_obu(SequenceHeader *seq_params,
                                  struct aom_read_bit_buffer *rb) {
    const uint32_t saved_bit_offset = rb->bit_offset;

    seq_params->profile =
        (BitstreamProfile)aom_rb_read_literal(rb, PROFILE_BITS);
    if (seq_params->profile > PROFILE_2) {
        seq_params->error = AOM_CODEC_UNSUP_BITSTREAM;
        return 0;
    }

    // Still picture or not
    seq_params->still_picture = aom_rb_read_bit(rb);
    seq_params->reduced_still_picture_header = aom_rb_read_bit(rb);
    // Video must have reduced_still_picture_hdr = 0
    if (!seq_params->still_picture &&
        seq_params->reduced_still_picture_header) {
        seq_params->error = AOM_CODEC_UNSUP_BITSTREAM;
        return 0;
    }

    if (seq_params->reduced_still_picture_header) {
        seq_params->timing_info_present_flag = 0;
        seq_params->decoder_model_info_present_flag = 0;
        seq_params->initial_display_delay_present_flag = 0;
        seq_params->operating_points_cnt_minus_1 = 0;
        seq_params->operating_point_idc[0] = 0;
        if (!read_bitstream_level(&seq_params->level[0], rb)) {
            seq_params->error = AOM_CODEC_UNSUP_BITSTREAM;
            return 0;
        }
        seq_params->tier[0] = 0;
        seq_params->op_params[0].decoder_model_present_for_this_op = 0;
        seq_params->op_params[0].initial_display_delay_present_for_this_op = 0;
    } else {
        // read timing info
        seq_params->timing_info_present_flag =
            aom_rb_read_bit(rb);  // timing_info_present_flag
        if (seq_params->timing_info_present_flag) {
            av1_read_timing_info_header(seq_params, rb);

            seq_params->decoder_model_info_present_flag = aom_rb_read_bit(rb);
            if (seq_params->decoder_model_info_present_flag)
                av1_read_decoder_model_info(seq_params, rb);
        } else {
            seq_params->decoder_model_info_present_flag = 0;
        }

        // read operating points
        seq_params->initial_display_delay_present_flag = aom_rb_read_bit(rb);
        seq_params->operating_points_cnt_minus_1 =
            aom_rb_read_literal(rb, OP_POINTS_CNT_MINUS_1_BITS);
        for (int i = 0; i < seq_params->operating_points_cnt_minus_1 + 1; i++) {
            seq_params->operating_point_idc[i] =
                aom_rb_read_literal(rb, OP_POINTS_IDC_BITS);
            if (!read_bitstream_level(&seq_params->level[i], rb)) {
                seq_params->error = AOM_CODEC_UNSUP_BITSTREAM;
                return 0;
            }
            // This is the seq_level_idx[i] > 7 check in the spec. seq_level_idx
            // 7 is equivalent to level 3.3.
            if (seq_params->level[i].major > 3)
                seq_params->tier[i] = aom_rb_read_bit(rb);
            else
                seq_params->tier[i] = 0;
            if (seq_params->decoder_model_info_present_flag) {
                seq_params->op_params[i].decoder_model_present_for_this_op =
                    aom_rb_read_bit(rb);
                if (seq_params->op_params[i].decoder_model_present_for_this_op)
                    av1_read_op_parameters_info(seq_params, rb, i);
            } else {
                seq_params->op_params[i].decoder_model_present_for_this_op = 0;
            }
            if (seq_params->timing_info_present_flag &&
                (seq_params->timing_info.equal_picture_interval ||
                 seq_params->op_params[i].decoder_model_present_for_this_op)) {
                seq_params->op_params[i].bitrate = max_level_bitrate(
                    seq_params->profile,
                    major_minor_to_seq_level_idx(seq_params->level[i]),
                    seq_params->tier[i]);
                // Level with seq_level_idx = 31 returns a high "dummy" bitrate
                // to pass the check
                if (seq_params->op_params[i].bitrate == 0)
                    seq_params->error = AOM_CODEC_UNSUP_BITSTREAM,
                    printf(
                        "AV1 does not support this combination of "
                        "profile, level, and tier.");
                // Buffer size in bits/s is bitrate in bits/s * 1 s
                seq_params->op_params[i].buffer_size =
                    seq_params->op_params[i].bitrate;
            }
            if (seq_params->timing_info_present_flag &&
                seq_params->timing_info.equal_picture_interval &&
                !seq_params->op_params[i].decoder_model_present_for_this_op) {
                // When the decoder_model_parameters are not sent for this op,
                // set the default ones that can be used with the resource
                // availability mode
                seq_params->op_params[i].decoder_buffer_delay = 70000;
                seq_params->op_params[i].encoder_buffer_delay = 20000;
                seq_params->op_params[i].low_delay_mode_flag = 0;
            }

            if (seq_params->initial_display_delay_present_flag) {
                seq_params->op_params[i]
                    .initial_display_delay_present_for_this_op =
                    aom_rb_read_bit(rb);
                if (seq_params->op_params[i]
                        .initial_display_delay_present_for_this_op) {
                    seq_params->op_params[i].initial_display_delay =
                        aom_rb_read_literal(rb, 4) + 1;
                    if (seq_params->op_params[i].initial_display_delay > 10)

                        seq_params->error = AOM_CODEC_UNSUP_BITSTREAM,
                        printf(
                            "AV1 does not support more than 10 decoded frames "
                            "delay");
                } else {
                    seq_params->op_params[i].initial_display_delay = 10;
                }
            } else {
                seq_params->op_params[i]
                    .initial_display_delay_present_for_this_op = 0;
                seq_params->op_params[i].initial_display_delay = 10;
            }
        }
    }

    // This decoder supports all levels.  Choose operating point provided by
    // external means
    int operating_point = 0;  // TODO: choose operating points
    if (operating_point < 0 ||
        operating_point > seq_params->operating_points_cnt_minus_1)
        operating_point = 0;
    int current_operating_point =
        seq_params->operating_point_idc[operating_point];
    if (aom_get_num_layers_from_operating_point_idc(
            current_operating_point,
            &seq_params->number_spatial_layers,
            &seq_params->number_temporal_layers) != AOM_CODEC_OK) {
        seq_params->error = AOM_CODEC_ERROR;
        return 0;
    }

    av1_read_sequence_header(rb, seq_params);

    av1_read_color_config(rb, seq_params);
    seq_params->film_grain_params_present = aom_rb_read_bit(rb);

    if (av1_check_trailing_bits(seq_params, rb) != 0) {
        // cm->error.error_code is already set.
        return 0;
    }

    seq_params->ready = 1;
    return ((rb->bit_offset - saved_bit_offset + 7) >> 3);
}

int parse_sequence_header_from_file(const char *ivf_file) {
    FILE *f = fopen(ivf_file, "rb");
    struct AvxInputContext input_ctx = {ivf_file, f, 0, 0, 0, {0, 0}};
    if (!file_is_ivf(&input_ctx)) {
        printf("File is NOT valid ivf\n");
        return -1;
    }

    uint8_t *stream_buf = reinterpret_cast<uint8_t *>(malloc(1024));
    size_t buf_sz = 1024;
    aom_codec_pts_t pts;
    size_t frame_sz = 0;
    int frame_cnt = 0;
    // keep reading ivf frames
    while (ivf_read_frame(f, &stream_buf, &frame_sz, &buf_sz, &pts) == 0) {
        printf("ivf frame count: %d\n", frame_cnt++);
        uint8_t *frame_buf = stream_buf;
        AomCodecErr err = AOM_CODEC_OK;
        do {
            // one ivf frame may contain multiple obus
            ObuHeader ou = {0};
            struct aom_read_bit_buffer rb = {
                frame_buf, frame_buf + frame_sz, 0, NULL, NULL};

            err = read_obu_header(&rb, 0, &ou);
            int header_size = ou.has_extension ? 2 : 1;
            frame_buf += header_size;
            frame_sz -= header_size;

            if (ou.has_size_field) {
                uint64_t u64_payload_length = 0;
                int value_len = 0;

                // read payload length
                for (int len = 0; len < OBU_MAX_LENGTH_FIELD_SIZE; ++len) {
                    if ((frame_buf[len] >> 7) == 0) {
                        ++len;
                        value_len = len;
                        break;
                    }
                }

                aom_uleb_decode(
                    frame_buf, value_len, &u64_payload_length, NULL);
                frame_buf += value_len;
                frame_sz -= value_len;
                printf("OBU type: %d, payload length: %u\n",
                       ou.type,
                       (uint32_t)u64_payload_length);

                // check the ou type and parse sequence header
                if (ou.type == OBU_SEQUENCE_HEADER) {
                    struct aom_read_bit_buffer rb = {
                        frame_buf, frame_buf + frame_sz, 0, NULL, NULL};
                    SequenceHeader sqs_headers = {0};
                    if (read_sequence_header_obu(&sqs_headers, &rb) == 0) {
                        printf("read seqence header fail\n");
                    }
                }

                frame_buf += u64_payload_length;
                frame_sz -= u64_payload_length;
            }
        } while (err == 0 && frame_sz > 0);
    }

    free(stream_buf);
    fclose(f);
    return 0;
}

void SequenceHeaderParser::input_obu_data(const uint8_t *obu_data,
                                          const uint32_t size,
                                          RefDecoder::StreamInfo *stream_info) {
    const uint8_t *frame_buf = obu_data;
    uint32_t frame_sz = size;
    AomCodecErr err = AOM_CODEC_OK;
    do {
        struct aom_read_bit_buffer rb = {
            frame_buf, frame_buf + frame_sz, 0, NULL, NULL};
        ObuHeader ou;
        AomCodecErr err = read_obu_header(&rb, 0, &ou);
        if (ou.has_size_field) {
            uint64_t u64_payload_length = 0;
            int header_size = ou.has_extension ? 2 : 1;
            int value_len = 0;

            frame_buf += header_size;
            frame_sz -= header_size;
            for (int len = 0; len < OBU_MAX_LENGTH_FIELD_SIZE; ++len) {
                if ((frame_buf[len] >> 7) == 0) {
                    ++len;
                    value_len = len;
                    break;
                }
            }
            aom_uleb_decode(frame_buf, value_len, &u64_payload_length, NULL);

            frame_buf += value_len;
            frame_sz -= value_len;
            if (ou.type == OBU_SEQUENCE_HEADER) {
                // check the ou type and parse sequence header
                struct aom_read_bit_buffer rb = {
                    frame_buf, frame_buf + frame_sz, 0, NULL, NULL};
                SequenceHeader sqs_headers = {0};
                ASSERT_NE(read_sequence_header_obu(&sqs_headers, &rb), 0)
                    << "read seqence header fail";
                profile_ = sqs_headers.profile;
                switch (sqs_headers.sb_size) {
                case BLOCK_64X64: sb_size_ = 64; break;
                case BLOCK_128X128: sb_size_ = 128; break;
                default:
                    ASSERT_TRUE(false)
                        << "super block size is invalid in sequence header!";
                    break;
                }
                printf("SPS header: profile(%u), sb_size(%u)\n",
                       profile_,
                       sb_size_);

                // update stream info
                stream_info->profile = sqs_headers.profile;
                stream_info->bit_depth = sqs_headers.bit_depth;
                stream_info->monochrome = sqs_headers.monochrome;
                stream_info->sb_size = sb_size_;
                stream_info->force_integer_mv = sqs_headers.force_integer_mv;
                stream_info->enable_filter_intra =
                    sqs_headers.enable_filter_intra;
                stream_info->enable_intra_edge_filter =
                    sqs_headers.enable_intra_edge_filter;
                stream_info->enable_masked_compound =
                    sqs_headers.enable_masked_compound;
                stream_info->enable_dual_filter =
                    sqs_headers.enable_dual_filter;
                stream_info->enable_jnt_comp = sqs_headers.enable_jnt_comp;
                stream_info->enable_ref_frame_mvs =
                    sqs_headers.enable_ref_frame_mvs;
                stream_info->enable_warped_motion =
                    sqs_headers.enable_warped_motion;
                stream_info->enable_cdef = sqs_headers.enable_cdef;
                stream_info->enable_restoration =
                    sqs_headers.enable_restoration;
            }
            frame_buf += u64_payload_length;
            frame_sz -= u64_payload_length;
        }
    } while (err == 0 && frame_sz > 0);
}

}  // namespace svt_av1_e2e_tools
