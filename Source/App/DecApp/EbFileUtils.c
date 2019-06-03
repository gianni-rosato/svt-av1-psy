/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>
#include <string.h>

#include "EbFileUtils.h"

const char *IVF_SIGNATURE = "DKIF";

static unsigned int mem_get_le16(const void *vmem) {
    unsigned int val;
    const uint8_t *mem = (const uint8_t *)vmem;

    val = mem[1] << 8;
    val |= mem[0];
    return val;
}

static unsigned int mem_get_le32(const void *vmem) {
    unsigned int val;
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

int file_is_ivf(CLInput *cli) {
    char raw_hdr[32];
    int is_ivf = 0;

    if (fread(raw_hdr, 1, 32, cli->inFile) == 32) {
        if (memcmp(IVF_SIGNATURE, raw_hdr, 4) == 0) {
            is_ivf = 1;

            if (mem_get_le16(raw_hdr + 4) != 0) {
                fprintf(stderr,
                    "Error: Unrecognized IVF version! This file may not"
                    " decode properly.");
            }

            cli->fourcc = mem_get_le32(raw_hdr + 8);
            cli->width = mem_get_le16(raw_hdr + 12);
            cli->height = mem_get_le16(raw_hdr + 14);
            cli->framerate.numerator = mem_get_le32(raw_hdr + 16);
            cli->framerate.denominator = mem_get_le32(raw_hdr + 20);
            fix_framerate(&cli->framerate.numerator,
                &cli->framerate.denominator);
        }
    }

    if (!is_ivf) {
        rewind(cli->inFile);
        cli->detect.buf_read = 0;
    }
    else
        cli->detect.position = 4;
    return is_ivf;
}

int read_ivf_frame(FILE *infile, uint8_t **buffer, size_t *bytes_read,
    size_t *buffer_size, int64_t *pts)
{
    char raw_header[IVF_FRAME_HDR_SZ] = { 0 };
    size_t frame_size = 0;

    if (fread(raw_header, IVF_FRAME_HDR_SZ, 1, infile) != 1) {
        if (!feof(infile))
            printf("Failed to read frame size. \n");
    }
    else {
        frame_size = mem_get_le32(raw_header);

        if (frame_size > 256 * 1024 * 1024) {
            printf("Read invalid frame size (%u) \n", (unsigned int)frame_size);
            frame_size = 0;
        }

        if (frame_size > *buffer_size) {
            uint8_t *new_buffer = (uint8_t *)realloc(*buffer, 2 * frame_size);

            if (new_buffer) {
                *buffer = new_buffer;
                *buffer_size = 2 * frame_size;
            }
            else {
                printf("Failed to allocate compressed data buffer. \n");
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
            printf("Failed to read full frame. \n");
            return 0;
        }
        *bytes_read = frame_size;
        return 1;
    }
    return 0;
}
