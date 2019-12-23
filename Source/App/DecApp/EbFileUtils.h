/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbFileUtils_h
#define EbFileUtils_h

#include <stdio.h>

#include "EbSvtAv1Dec.h"

#define OBU_BUFFER_SIZE (500 * 1024)

#define OBU_HEADER_SIZE 1
#define OBU_EXTENSION_SIZE 1
#define OBU_MAX_LENGTH_FIELD_SIZE 8

#define OBU_MAX_HEADER_SIZE (OBU_HEADER_SIZE + OBU_EXTENSION_SIZE + 2 * OBU_MAX_LENGTH_FIELD_SIZE)

#define OBU_DETECTION_SIZE (OBU_HEADER_SIZE + OBU_EXTENSION_SIZE + 4 * OBU_MAX_LENGTH_FIELD_SIZE)

#define IVF_FRAME_HDR_SZ (4 + 8) /* 4 byte size + 8 byte timestamp */

#define DECAPP_MIN(x, y) (((x) < (y)) ? (x) : (y))
#define DECAPP_MAX(x, y) (((x) > (y)) ? (x) : (y))

static const char *const csp_names[] = {"400", "420", "422", "444", 0};

/**********************************
 * Input image properties
 **********************************/

enum VideoFileType { FILE_TYPE_OBU, FILE_TYPE_RAW, FILE_TYPE_IVF, FILE_TYPE_Y4M, FILE_TYPE_WEBM };

struct FileTypeDetectionBuffer {
    char   buf[4];
    size_t buf_read;
    size_t position;
};

struct Rational {
    int numerator;
    int denominator;
};

typedef struct CliInput {
    const char *                   in_filename;
    const char *                   out_filename;
    FILE *                         in_file;
    FILE *                         out_file;
    uint32_t                       width;
    uint32_t                       height;
    uint32_t                       fourcc;
    enum VideoFileType             in_file_type;
    struct FileTypeDetectionBuffer detect;
    struct Rational                pixel_aspect_ratio;
    struct Rational                framerate;
    EbColorFormat                  fmt;
    uint32_t                       enable_md5;
    uint32_t                       fps_frm;
    uint32_t                       fps_summary;
    uint32_t                       skip_film_grain;
} CliInput;

typedef struct ObuDecInputContext {
    uint8_t *buffer;
    size_t   buffer_capacity;
    size_t   bytes_buffered;
    uint32_t is_annexb;
    uint64_t rem_txb_size;
} ObuDecInputContext;

typedef struct DecInputContext {
    CliInput *          cli_ctx;
    ObuDecInputContext *obu_ctx;
} DecInputContext;

/*!\brief OBU types. */
typedef enum ATTRIBUTE_PACKED {
    OBU_SEQUENCE_HEADER        = 1,
    OBU_TEMPORAL_DELIMITER     = 2,
    OBU_FRAME_HEADER           = 3,
    OBU_TILE_GROUP             = 4,
    OBU_METADATA               = 5,
    OBU_FRAME                  = 6,
    OBU_REDUNDANT_FRAME_HEADER = 7,
    OBU_TILE_LIST              = 8,
    OBU_PADDING                = 15,
} OBU_TYPE;

typedef struct {
    size_t size; // Size (1 or 2 bytes) of the OBU header (including the
        // optional OBU extension header) in the Bitstream.
    OBU_TYPE type;
    int      has_size_field;
    int      has_extension;
    // The following fields come from the OBU extension header and therefore are
    // only used if has_extension is true.
    int temporal_layer_id;
    int spatial_layer_id;
} ObuHeader;

int file_is_obu(CliInput *cli, ObuDecInputContext *obu_ctx);
int obudec_read_temporal_unit(DecInputContext *input, uint8_t **buffer, size_t *bytes_read,
                              size_t *buffer_size);

int file_is_ivf(CliInput *cli);
int read_ivf_frame(FILE *infile, uint8_t **buffer, size_t *bytes_read, size_t *buffer_size,
                   int64_t *pts);

#endif
