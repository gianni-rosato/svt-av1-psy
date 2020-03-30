/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <assert.h>

#include "EbFileUtils.h"

const char *ivf_signature = "DKIF";

static const size_t  k_maximum_leb_128_size = 8;
static const uint8_t k_leb_128byte_mask     = 0x7f; // Binary: 01111111

static unsigned int mem_get_le16(const void *vmem) {
    unsigned int   val;
    const uint8_t *mem = (const uint8_t *)vmem;

    val = mem[1] << 8;
    val |= mem[0];
    return val;
}

static unsigned int mem_get_le32(const void *vmem) {
    unsigned int   val;
    const uint8_t *mem = (const uint8_t *)vmem;

    val = ((unsigned int)mem[3]) << 24;
    val |= mem[2] << 16;
    val |= mem[1] << 8;
    val |= mem[0];
    return val;
}

static void fix_framerate(int *num, int *den) {
    if (*den <= 0 || *den >= 1000000000 || *num <= 0 || *num >= 1000) {
        // framerate seems to be invalid, just default to 30fps.
        *num = 30;
        *den = 1;
    }
}

typedef struct ReadBitBuffer {
    const uint8_t *bit_buffer;
    const uint8_t *bit_buffer_end;
    uint32_t       bit_offset;
} ReadBitBuffer;

int uleb_decode(const uint8_t *buffer, size_t available, uint64_t *value, size_t *length) {
    if (buffer && value) {
        *value = 0;
        for (size_t i = 0; i < k_maximum_leb_128_size && i < available; ++i) {
            const uint8_t decoded_byte = *(buffer + i) & k_leb_128byte_mask;
            *value |= ((uint64_t)decoded_byte) << (i * 7);
            if ((*(buffer + i) >> 7) == 0) {
                if (length) { *length = i + 1; }

                // Fail on values larger than 32-bits to ensure consistent
                // behavior on 32 and 64 bit targets: value is typically
                // used to determine buffer allocation size.
                if (*value > UINT32_MAX) return -1;

                return 0;
            }
        }
    }

    // If we get here, either the buffer/value pointers were invalid,
    // or we ran over the available space
    return -1;
}

// Reads unsigned LEB128 integer and returns 0 upon successful read and decode.
// Stores raw bytes in 'value_buffer', length of the number in 'value_length',
// and decoded value in 'value'.
static int obudec_read_leb128(FILE *f, uint8_t *value_buffer, size_t *value_length,
                              uint64_t *value) {
    if (!f || !value_buffer || !value_length || !value) return -1;
    size_t len;
    for (len = 0; len < OBU_MAX_LENGTH_FIELD_SIZE; ++len) {
        const size_t num_read = fread(&value_buffer[len], 1, 1, f);
        if (num_read == 0) {
            if (len == 0 && feof(f)) {
                *value_length = 0;
                return 0;
            }
            // Ran out of data before completing read of value.
            return -1;
        }
        if ((value_buffer[len] >> 7) == 0) {
            ++len;
            *value_length = len;
            break;
        }
    }

    return uleb_decode(value_buffer, len, value, NULL);
}

int rb_read_bit(ReadBitBuffer *rb) {
    const uint32_t off = rb->bit_offset;
    const uint32_t p   = off >> 3;
    const int      q   = 7 - (int)(off & 0x7);
    if (rb->bit_buffer + p < rb->bit_buffer_end) {
        const int bit  = (rb->bit_buffer[p] >> q) & 1;
        rb->bit_offset = off + 1;
        return bit;
    } else
        return 0;
}

int rb_read_literal(ReadBitBuffer *rb, int bits) {
    assert(bits <= 31);
    int value = 0, bit;
    for (bit = bits - 1; bit >= 0; bit--) value |= rb_read_bit(rb) << bit;
    return value;
}

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
    case OBU_TILE_LIST:
    case OBU_PADDING: valid_type = 1; break;
    default: break;
    }
    return valid_type;
}

// Parses OBU header and stores values in 'header'.
static int read_obu_header(ReadBitBuffer *rb, uint32_t is_annexb, ObuHeader *header) {
    if (!rb || !header) return -1;

    const ptrdiff_t bit_buffer_byte_length = rb->bit_buffer_end - rb->bit_buffer;
    if (bit_buffer_byte_length < 1) return -1;

    header->size = 1;

    if (rb_read_bit(rb) != 0) {
        // Forbidden bit. Must not be set.
        return -1;
    }

    header->type = (OBU_TYPE)rb_read_literal(rb, 4);

    if (!valid_obu_type(header->type)) return -1;

    header->has_extension  = rb_read_bit(rb);
    header->has_size_field = rb_read_bit(rb);

    if (!header->has_size_field && !is_annexb) {
        // section 5 obu streams must have obu_size field set.
        return -1;
    }

    if (rb_read_bit(rb) != 0) {
        // obu_reserved_1bit must be set to 0.
        return -1;
    }

    if (header->has_extension) {
        if (bit_buffer_byte_length == 1) return -1;

        header->size += 1;
        header->temporal_layer_id = rb_read_literal(rb, 3);
        header->spatial_layer_id  = rb_read_literal(rb, 2);
        if (rb_read_literal(rb, 3) != 0) {
            // extension_header_reserved_3bits must be set to 0.
            return -1;
        }
    }

    return 0;
}

int svt_read_obu_header(uint8_t *buffer, size_t buffer_length, size_t *consumed, ObuHeader *header,
                        uint32_t is_annexb) {
    if (buffer_length < 1 || !consumed || !header) return -1;

    ReadBitBuffer rb           = {buffer, buffer + buffer_length, 0};
    int           parse_result = read_obu_header(&rb, is_annexb, header);
    if (parse_result == 0) *consumed = header->size;
    return parse_result;
}

// Reads OBU header from 'f'. The 'buffer_capacity' passed in must be large
// enough to store an OBU header with extension (2 bytes). Raw OBU data is
// written to 'obu_data', parsed OBU header values are written to 'obu_header',
// and total bytes read from file are written to 'bytes_read'. Returns 0 for
// success, and non-zero on failure. When end of file is reached, the return
// value is 0 and the 'bytes_read' value is set to 0.
static int obudec_read_obu_header(FILE *f, size_t buffer_capacity, uint32_t is_annexb,
                                  uint8_t *obu_data, ObuHeader *obu_header, size_t *bytes_read) {
    if (!f || buffer_capacity < (OBU_HEADER_SIZE + OBU_EXTENSION_SIZE) || !obu_data ||
        !obu_header || !bytes_read) {
        return -1;
    }
    *bytes_read = fread(obu_data, 1, 1, f);

    if (feof(f) && *bytes_read == 0) {
        return 0;
    } else if (*bytes_read != 1) {
        fprintf(stderr, "obudec: Failure reading OBU header.\n");
        return -1;
    }

    const int has_extension = (obu_data[0] >> 2) & 0x1;
    if (has_extension) {
        if (fread(&obu_data[1], 1, 1, f) != 1) {
            fprintf(stderr, "obudec: Failure reading OBU extension.");
            return -1;
        }
        ++*bytes_read;
    }

    size_t obu_bytes_parsed = 0;
    svt_read_obu_header(obu_data, *bytes_read, &obu_bytes_parsed, obu_header, is_annexb);

    return 0;
}

static int obudec_read_obu_header_and_size(FILE *f, size_t buffer_capacity, uint32_t is_annexb,
                                           uint8_t *buffer, size_t *bytes_read,
                                           size_t *payload_length, ObuHeader *obu_header) {
    const size_t k_minimum_buffer_size = OBU_MAX_HEADER_SIZE;
    if (!f || !buffer || !bytes_read || !payload_length || !obu_header ||
        buffer_capacity < k_minimum_buffer_size) {
        return -1;
    }

    size_t   leb128_length_obu     = 0;
    size_t   leb128_length_payload = 0;
    uint64_t obu_size              = 0;
    if (is_annexb) {
        if (obudec_read_leb128(f, &buffer[0], &leb128_length_obu, &obu_size) != 0) {
            fprintf(stderr, "obudec: Failure reading OBU size length.\n");
            return -1;
        } else if (leb128_length_obu == 0) {
            *payload_length = 0;
            return 0;
        }
        if (obu_size > UINT32_MAX) {
            fprintf(stderr, "obudec: OBU payload length too large.\n");
            return -1;
        }
    }

    size_t header_size = 0;
    if (obudec_read_obu_header(f,
                               buffer_capacity - leb128_length_obu,
                               is_annexb,
                               buffer + leb128_length_obu,
                               obu_header,
                               &header_size) != 0) {
        return -1;
    } else if (header_size == 0) {
        *payload_length = 0;
        return 0;
    }

    if (!obu_header->has_size_field) {
        assert(is_annexb);
        if (obu_size < header_size) {
            fprintf(stderr, "obudec: OBU size is too small.\n");
            return -1;
        }
        *payload_length = (size_t)obu_size - header_size;
    } else {
        uint64_t u64_payload_length = 0;
        if (obudec_read_leb128(f,
                               &buffer[leb128_length_obu + header_size],
                               &leb128_length_payload,
                               &u64_payload_length) != 0) {
            fprintf(stderr, "obudec: Failure reading OBU payload length.\n");
            return -1;
        }
        if (u64_payload_length > UINT32_MAX) {
            fprintf(stderr, "obudec: OBU payload length too large.\n");
            return -1;
        }

        *payload_length = (size_t)u64_payload_length;
    }

    *bytes_read = leb128_length_obu + header_size + leb128_length_payload;
    return 0;
}

// Reads OBU payload from 'f' and returns 0 for success when all payload bytes
// are read from the file. Payload data is written to 'obu_data', and actual
// bytes read added to 'bytes_read'.
static int obudec_read_obu_payload(FILE *f, size_t payload_length, uint8_t *obu_data,
                                   size_t *bytes_read) {
    if (!f || payload_length == 0 || !obu_data || !bytes_read) return -1;

    if (fread(obu_data, 1, payload_length, f) != payload_length) {
        fprintf(stderr, "obudec: Failure reading OBU payload.\n");
        return -1;
    }

    *bytes_read += payload_length;
    return 0;
}

int file_is_obu(CliInput *cli, ObuDecInputContext *obu_ctx) {
    if (!obu_ctx || !cli) return 0;

    uint8_t        detect_buf[OBU_DETECTION_SIZE] = {0};
    const uint32_t is_annexb                      = obu_ctx->is_annexb;
    FILE *         f                              = cli->in_file;
    size_t         payload_length                 = 0;
    ObuHeader      obu_header;
    memset(&obu_header, 0, sizeof(obu_header));
    size_t   length_of_unit_size  = 0;
    size_t   annexb_header_length = 0;
    uint64_t unit_size            = 0;

    if (is_annexb) {
        // read the size of first temporal unit
        if (obudec_read_leb128(f, &detect_buf[0], &length_of_unit_size, &unit_size) != 0) {
            fprintf(stderr, "obudec: Failure reading temporal unit header\n");
            return 0;
        }

        // read the size of first frame unit
        if (obudec_read_leb128(
                f, &detect_buf[length_of_unit_size], &annexb_header_length, &unit_size) != 0) {
            fprintf(stderr, "obudec: Failure reading frame unit header\n");
            return 0;
        }
        annexb_header_length += length_of_unit_size;
    }

    size_t bytes_read = 0;
    if (obudec_read_obu_header_and_size(f,
                                        OBU_DETECTION_SIZE - annexb_header_length,
                                        is_annexb,
                                        &detect_buf[annexb_header_length],
                                        &bytes_read,
                                        &payload_length,
                                        &obu_header) != 0) {
        fprintf(stderr, "obudec: Failure reading first OBU.\n");
        rewind(f);
        return 0;
    }

    if (is_annexb) { bytes_read += annexb_header_length; }

    if (obu_header.type != OBU_TEMPORAL_DELIMITER && obu_header.type != OBU_SEQUENCE_HEADER) {
        return 0;
    }

    if (obu_header.has_size_field) {
        if (obu_header.type == OBU_TEMPORAL_DELIMITER && payload_length != 0) {
            fprintf(stderr, "obudec: Invalid OBU_TEMPORAL_DELIMITER payload length (non-zero).");
            rewind(f);
            return 0;
        }
    } else if (!is_annexb) {
        fprintf(stderr, "obudec: OBU size fields required, cannot decode input.\n");
        rewind(f);
        return 0;
    }

    // Appears that input is valid Section 5 AV1 stream.
    obu_ctx->buffer = (uint8_t *)malloc(OBU_BUFFER_SIZE);
    if (!obu_ctx->buffer) {
        fprintf(stderr, "Out of memory.\n");
        rewind(f);
        return 0;
    }
    obu_ctx->buffer_capacity = OBU_BUFFER_SIZE;

    memcpy(obu_ctx->buffer, &detect_buf[0], bytes_read);
    obu_ctx->bytes_buffered = bytes_read;
    // If the first OBU is a SEQUENCE_HEADER, then it will have a payload.
    // We need to read this in so that our buffer only contains complete OBUs.
    if (payload_length > 0) {
        if (payload_length > (obu_ctx->buffer_capacity - bytes_read)) {
            fprintf(stderr, "obudec: First OBU's payload is too large\n");
            rewind(f);
            return 0;
        }

        size_t    payload_bytes = 0;
        const int status        = obudec_read_obu_payload(
            f, payload_length, &obu_ctx->buffer[bytes_read], &payload_bytes);
        if (status < 0) {
            rewind(f);
            return 0;
        }
        obu_ctx->bytes_buffered += payload_bytes;
    }

    /* This is because to avoid to many conditions while reading
    frame by frame information in TU's */
    if (is_annexb) {
        rewind(f);
        obu_ctx->bytes_buffered = 0;
    }
    return 1;
}

static int obudec_grow_buffer(size_t growth_amount, uint8_t **obu_buffer,
                              size_t *obu_buffer_capacity) {
    if (!*obu_buffer || !obu_buffer_capacity || growth_amount == 0) { return -1; }

    const size_t capacity = *obu_buffer_capacity;
    if (SIZE_MAX - growth_amount < capacity) {
        fprintf(stderr, "obudec: cannot grow buffer, capacity will roll over.\n");
        return -1;
    }

    const size_t new_capacity = capacity + growth_amount;

    uint8_t *new_buffer = (uint8_t *)realloc(*obu_buffer, new_capacity);
    if (!new_buffer) {
        fprintf(stderr, "obudec: Failed to allocate compressed data buffer.\n");
        return -1;
    }

    *obu_buffer          = new_buffer;
    *obu_buffer_capacity = new_capacity;
    return 0;
}

static int obudec_read_one_obu(FILE *f, uint8_t **obu_buffer, size_t obu_bytes_buffered,
                               size_t *obu_buffer_capacity, size_t *obu_length,
                               ObuHeader *obu_header, uint32_t is_annexb) {
    if (!f || !(*obu_buffer) || !obu_buffer_capacity || !obu_length || !obu_header) { return -1; }

    size_t bytes_read                = 0;
    size_t obu_payload_length        = 0;
    size_t available_buffer_capacity = *obu_buffer_capacity - obu_bytes_buffered;

    if (available_buffer_capacity < OBU_MAX_HEADER_SIZE) {
        if (obudec_grow_buffer(DECAPP_MAX(*obu_buffer_capacity, OBU_MAX_HEADER_SIZE),
                               obu_buffer,
                               obu_buffer_capacity) != 0) {
            *obu_length = bytes_read;
            return -1;
        }
        available_buffer_capacity += DECAPP_MAX(*obu_buffer_capacity, OBU_MAX_HEADER_SIZE);
    }

    const int status = obudec_read_obu_header_and_size(f,
                                                       available_buffer_capacity,
                                                       is_annexb,
                                                       *obu_buffer + obu_bytes_buffered,
                                                       &bytes_read,
                                                       &obu_payload_length,
                                                       obu_header);
    if (status < 0) return status;

    if (obu_payload_length > SIZE_MAX - bytes_read) return -1;

    if (obu_payload_length > 256 * 1024 * 1024) {
        fprintf(stderr, "obudec: Read invalid OBU size (%u)\n", (unsigned int)obu_payload_length);
        *obu_length = bytes_read + obu_payload_length;
        return -1;
    }

    if (bytes_read + obu_payload_length > available_buffer_capacity &&
        obudec_grow_buffer(DECAPP_MAX(*obu_buffer_capacity, obu_payload_length),
                           obu_buffer,
                           obu_buffer_capacity) != 0) {
        *obu_length = bytes_read + obu_payload_length;
        return -1;
    }

    if (obu_payload_length > 0 &&
        obudec_read_obu_payload(
            f, obu_payload_length, *obu_buffer + obu_bytes_buffered + bytes_read, &bytes_read) !=
            0) {
        return -1;
    }

    *obu_length = bytes_read;
    return 0;
}

int obudec_read_temporal_unit(DecInputContext *input, uint8_t **buffer, size_t *bytes_read,
                              size_t *buffer_size) {
    CliInput *          cli     = input->cli_ctx;
    FILE *              f       = cli->in_file;
    ObuDecInputContext *obu_ctx = input->obu_ctx;
    if (!f) return 0;

    *buffer_size = 0;
    *bytes_read  = 0;

    if (feof(f)) { return 0; }

    size_t  txb_size = 0, fr_size = 0;
    size_t  obu_size                            = 0;
    size_t  length_of_temporal_unit_size        = 0;
    size_t  length_of_frame_unit_size           = 0;
    uint8_t frheader[OBU_MAX_LENGTH_FIELD_SIZE] = {0};

    if (obu_ctx->is_annexb) {
        uint64_t size = 0;

        assert(obu_ctx->bytes_buffered == 0);

        if (!obu_ctx->rem_txb_size) {
            if (obudec_read_leb128(f, &frheader[0], &length_of_temporal_unit_size, &size) != 0) {
                fprintf(stderr, "obudec: Failure reading temporal unit header\n");
                return 0;
            }
            if (size == 0 && feof(f)) { return 0; }
            /*Stores only tu size ie excluding tu header*/
            obu_ctx->rem_txb_size = size;
        }

        if (size > UINT32_MAX || size + length_of_temporal_unit_size > UINT32_MAX) {
            fprintf(stderr, "obudec: TU too large.\n");
            return 0;
        }

        if (obudec_read_leb128(f, &frheader[0], &length_of_frame_unit_size, &size) != 0) {
            fprintf(stderr, "obudec: Failure reading frame header\n");
            return 0;
        }
        if (size == 0 || feof(f)) { return 0; }

        fr_size  = (size_t)size;
        txb_size = fr_size;
    } else {
        while (1) {
            ObuHeader obu_header;
            memset(&obu_header, 0, sizeof(obu_header));

            if (obudec_read_one_obu(f,
                                    &obu_ctx->buffer,
                                    obu_ctx->bytes_buffered,
                                    &obu_ctx->buffer_capacity,
                                    &obu_size,
                                    &obu_header,
                                    0) != 0) {
                fprintf(stderr, "obudec: read_one_obu failed in TU loop\n");
                return 0;
            }

            if (obu_header.type == OBU_TEMPORAL_DELIMITER || obu_size == 0) {
                txb_size = obu_ctx->bytes_buffered;
                break;
            } else {
                obu_ctx->bytes_buffered += obu_size;
            }
        }
    }

    uint8_t *new_buffer = (uint8_t *)realloc(*buffer, txb_size);
    if (!new_buffer) {
        free(*buffer);
        fprintf(stderr, "obudec: Out of memory.\n");
        return 0;
    }
    *buffer      = new_buffer;
    *bytes_read  = txb_size;
    *buffer_size = txb_size;

    if (!obu_ctx->is_annexb) {
        memcpy(*buffer, obu_ctx->buffer, txb_size);

        // At this point, (obu_ctx->buffer + obu_ctx->bytes_buffered + obu_size)
        // points to the end of the buffer.
        memmove(obu_ctx->buffer, obu_ctx->buffer + obu_ctx->bytes_buffered, obu_size);
        obu_ctx->bytes_buffered = obu_size;
    } else {
        if (!feof(f)) {
            if (fread(*buffer, 1, fr_size, f) != fr_size) {
                fprintf(stderr, "obudec: Failed to read full temporal unit\n");
                return 0;
            }
            obu_ctx->rem_txb_size -= (fr_size + length_of_frame_unit_size);
        }
    }
    return 1;
}

int file_is_ivf(CliInput *cli) {
    char raw_hdr[32];
    int  is_ivf = 0;

    if (fread(raw_hdr, 1, 32, cli->in_file) == 32) {
        if (memcmp(ivf_signature, raw_hdr, 4) == 0) {
            is_ivf = 1;

            if (mem_get_le16(raw_hdr + 4) != 0) {
                fprintf(stderr,
                        "Error: Unrecognized IVF version! This file may not"
                        " decode properly.");
            }

            cli->fourcc                = mem_get_le32(raw_hdr + 8);
            cli->width                 = mem_get_le16(raw_hdr + 12);
            cli->height                = mem_get_le16(raw_hdr + 14);
            cli->framerate.numerator   = mem_get_le32(raw_hdr + 16);
            cli->framerate.denominator = mem_get_le32(raw_hdr + 20);
            fix_framerate(&cli->framerate.numerator, &cli->framerate.denominator);
        }
    }

    if (!is_ivf) {
        rewind(cli->in_file);
        cli->detect.buf_read = 0;
    } else
        cli->detect.position = 4;
    return is_ivf;
}

int read_ivf_frame(FILE *infile, uint8_t **buffer, size_t *bytes_read, size_t *buffer_size,
                   int64_t *pts) {
    char   raw_header[IVF_FRAME_HDR_SZ] = {0};
    size_t frame_size                   = 0;

    if (fread(raw_header, IVF_FRAME_HDR_SZ, 1, infile) != 1) {
        if (!feof(infile)) fprintf(stderr, "Failed to read frame size. \n");
    } else {
        frame_size = mem_get_le32(raw_header);

        if (frame_size > 256 * 1024 * 1024) {
            fprintf(stderr, "Read invalid frame size (%u) \n", (unsigned int)frame_size);
            frame_size = 0;
        }

        if (frame_size > *buffer_size) {
            uint8_t *new_buffer = (uint8_t *)realloc(*buffer, 2 * frame_size);

            if (new_buffer) {
                *buffer      = new_buffer;
                *buffer_size = 2 * frame_size;
            } else {
                fprintf(stderr, "Failed to allocate compressed data buffer. \n");
                frame_size = 0;
            }
        }

        if (pts) {
            *pts = mem_get_le32(&raw_header[4]);
            *pts += ((int64_t)mem_get_le32(&raw_header[8]) << 32);
        }
    }

    if (!feof(infile)) {
        if (fread(*buffer, 1, frame_size, infile) != frame_size) {
            fprintf(stderr, "Failed to read full frame. \n");
            return 0;
        }
        *bytes_read = frame_size;
        return 1;
    }
    return 0;
}
