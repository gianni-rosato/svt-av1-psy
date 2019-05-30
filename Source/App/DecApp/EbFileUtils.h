/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdio.h>

#include "EbSvtAv1Dec.h"

#define IVF_FRAME_HDR_SZ (4 + 8) /* 4 byte size + 8 byte timestamp */

static const char * const csp_names[] = { "400", "420", "422", "444", 0 };

/**********************************
 * Input image properties
 **********************************/

enum VideoFileType {
    FILE_TYPE_OBU,
    FILE_TYPE_RAW,
    FILE_TYPE_IVF,
    FILE_TYPE_Y4M,
    FILE_TYPE_WEBM
};

struct FileTypeDetectionBuffer {
    char buf[4];
    size_t buf_read;
    size_t position;
};

struct Rational {
    int numerator;
    int denominator;
};

typedef struct CLInput{
    const char *inFilename;
    const char *outFilename;
    FILE *inFile;
    FILE *outFile;
    uint32_t width;
    uint32_t height;
    uint32_t fourcc;
    enum VideoFileType inFileType;
    struct FileTypeDetectionBuffer detect;
    struct Rational pixel_aspect_ratio;
    struct Rational framerate;
    EbColorFormat fmt;
    EbBitDepth bit_depth;
    uint32_t   enable_md5;
}CLInput;

int file_is_ivf(CLInput *cli);
int read_ivf_frame(FILE *infile, uint8_t **buffer, size_t *bytes_read,
    size_t *buffer_size, int64_t *pts);
