/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbAppString.h"
#include "EbAppInputy4m.h"
#define YFM_HEADER_MAX 80
#define YUV4MPEG2_IND_SIZE 9
#define PRINT_HEADER 0
#define CHROMA_MAX 4

/* copy a string until a specified character or a new line is found */
char *copy_until_char_or_newline(char *src, char *dst, char chr) {
    rsize_t count    = 0;
    char *  src_init = src;

    while (*src != chr && *src != '\n') {
        src++;
        count++;
    }

    EB_STRNCPY(dst, YFM_HEADER_MAX, src_init, count);

    return src;
}

/* reads the y4m header and parses the input parameters */
int32_t read_y4m_header(EbConfig *cfg) {
    FILE *   ptr_in;
    char     buffer[YFM_HEADER_MAX];
    char *   fresult, *tokstart, *tokend, format_str[YFM_HEADER_MAX];
    uint32_t bitdepth = 8, width = 0, height = 0, fr_n = 0, fr_d = 0, aspect_n, aspect_d;
    char     chroma[CHROMA_MAX] = "420", scan_type = 'p';

    /* pointer to the input file */
    ptr_in = cfg->input_file;

    /* get first line after YUV4MPEG2 */
    fresult = fgets(buffer, sizeof(buffer), ptr_in);
    if (fresult == NULL) return EB_ErrorBadParameter;

    /* print header */
    if (PRINT_HEADER) {
        fprintf(stderr, "y4m header:");
        fputs(buffer, stdout);
    }

    /* read header parameters */
    tokstart = &(buffer[0]);

    while (*tokstart != '\0') {
        if (*tokstart == 0x20) {
            tokstart++;
            continue;
        }

        switch (*tokstart++) {
        case 'W': /* width, required. */
            width = (uint32_t)strtol(tokstart, &tokend, 10);
            if (PRINT_HEADER) fprintf(stderr, "width = %d\n", width);
            tokstart = tokend;
            break;
        case 'H': /* height, required. */
            height = (uint32_t)strtol(tokstart, &tokend, 10);
            if (PRINT_HEADER) fprintf(stderr, "height = %d\n", height);
            tokstart = tokend;
            break;
        case 'I': /* scan type, not required, default: 'p' */
            switch (*tokstart++) {
            case 'p': scan_type = 'p'; break;
            case 't': scan_type = 't'; break;
            case 'b': scan_type = 'b'; break;
            case '?':
            default:
                fprintf(cfg->error_log_file, "interlace type not supported\n");
                return EB_ErrorBadParameter;
            }
            if (PRINT_HEADER) fprintf(stderr, "scan_type = %c\n", scan_type);
            break;
        case 'C': /* color space, not required: default "420" */
            tokstart = copy_until_char_or_newline(tokstart, format_str, 0x20);
            if (EB_STRCMP("420mpeg2", format_str) == 0) {
                EB_STRCPY(chroma, CHROMA_MAX, "420");
                // chroma left
                bitdepth = 8;
            } else if (EB_STRCMP("420paldv", format_str) == 0) {
                EB_STRCPY(chroma, CHROMA_MAX, "420");
                // chroma top-left
                bitdepth = 8;
            } else if (EB_STRCMP("420jpeg", format_str) == 0) {
                EB_STRCPY(chroma, CHROMA_MAX, "420");
                // chroma center
                bitdepth = 8;
            } else if (EB_STRCMP("420p16", format_str) == 0) {
                EB_STRCPY(chroma, CHROMA_MAX, "420");
                bitdepth = 16;
            } else if (EB_STRCMP("422p16", format_str) == 0) {
                EB_STRCPY(chroma, CHROMA_MAX, "422");
                bitdepth = 16;
            } else if (EB_STRCMP("444p16", format_str) == 0) {
                EB_STRCPY(chroma, CHROMA_MAX, "444");
                bitdepth = 16;
            } else if (EB_STRCMP("420p14", format_str) == 0) {
                EB_STRCPY(chroma, CHROMA_MAX, "420");
                bitdepth = 14;
            } else if (EB_STRCMP("422p14", format_str) == 0) {
                EB_STRCPY(chroma, CHROMA_MAX, "422");
                bitdepth = 14;
            } else if (EB_STRCMP("444p14", format_str) == 0) {
                EB_STRCPY(chroma, CHROMA_MAX, "444");
                bitdepth = 14;
            } else if (EB_STRCMP("420p12", format_str) == 0) {
                EB_STRCPY(chroma, CHROMA_MAX, "420");
                bitdepth = 12;
            } else if (EB_STRCMP("422p12", format_str) == 0) {
                EB_STRCPY(chroma, CHROMA_MAX, "422");
                bitdepth = 12;
            } else if (EB_STRCMP("444p12", format_str) == 0) {
                EB_STRCPY(chroma, CHROMA_MAX, "444");
                bitdepth = 12;
            } else if (EB_STRCMP("420p10", format_str) == 0) {
                EB_STRCPY(chroma, CHROMA_MAX, "420");
                bitdepth = 10;
            } else if (EB_STRCMP("422p10", format_str) == 0) {
                EB_STRCPY(chroma, CHROMA_MAX, "422");
                bitdepth = 10;
            } else if (EB_STRCMP("444p10", format_str) == 0) {
                EB_STRCPY(chroma, CHROMA_MAX, "444");
                bitdepth = 10;
            } else if (EB_STRCMP("420p9", format_str) == 0) {
                EB_STRCPY(chroma, CHROMA_MAX, "420");
                bitdepth = 9;
            } else if (EB_STRCMP("422p9", format_str) == 0) {
                EB_STRCPY(chroma, CHROMA_MAX, "422");
                bitdepth = 9;
            } else if (EB_STRCMP("444p9", format_str) == 0) {
                EB_STRCPY(chroma, CHROMA_MAX, "444");
                bitdepth = 9;
            } else if (EB_STRCMP("420", format_str) == 0) {
                EB_STRCPY(chroma, CHROMA_MAX, "420");
                bitdepth = 8;
            } else if (EB_STRCMP("411", format_str) == 0) {
                EB_STRCPY(chroma, CHROMA_MAX, "411");
                bitdepth = 8;
            } else if (EB_STRCMP("422", format_str) == 0) {
                EB_STRCPY(chroma, CHROMA_MAX, "422");
                bitdepth = 8;
            } else if (EB_STRCMP("444", format_str) == 0) {
                EB_STRCPY(chroma, CHROMA_MAX, "444");
                bitdepth = 8;
            } else if (EB_STRCMP("mono16", format_str) == 0) {
                EB_STRCPY(chroma, CHROMA_MAX, "400");
                bitdepth = 16;
            } else if (EB_STRCMP("mono12", format_str) == 0) {
                EB_STRCPY(chroma, CHROMA_MAX, "400");
                bitdepth = 12;
            } else if (EB_STRCMP("mono10", format_str) == 0) {
                EB_STRCPY(chroma, CHROMA_MAX, "400");
                bitdepth = 10;
            } else if (EB_STRCMP("mono9", format_str) == 0) {
                EB_STRCPY(chroma, CHROMA_MAX, "400");
                bitdepth = 9;
            } else if (EB_STRCMP("mono", format_str) == 0) {
                EB_STRCPY(chroma, CHROMA_MAX, "400");
                bitdepth = 8;
            } else {
                fprintf(cfg->error_log_file, "chroma format not supported\n");
                return EB_ErrorBadParameter;
            }
            if (PRINT_HEADER) fprintf(stderr, "chroma = %s, bitdepth = %d\n", chroma, bitdepth);
            break;
        case 'F': /* frame rate, required */
            tokstart = copy_until_char_or_newline(tokstart, format_str, ':');
            fr_n     = (uint32_t)strtol(format_str, (char **)NULL, 10);
            tokstart++;
            tokstart = copy_until_char_or_newline(tokstart, format_str, 0x20);
            fr_d     = (uint32_t)strtol(format_str, (char **)NULL, 10);
            if (PRINT_HEADER) {
                fprintf(stderr, "framerate_n = %d\n", fr_n);
                fprintf(stderr, "framerate_d = %d\n", fr_d);
            }
            break;
        case 'A': /* aspect ratio, not required */
            tokstart = copy_until_char_or_newline(tokstart, format_str, ':');
            aspect_n = (uint32_t)strtol(format_str, (char **)NULL, 10);
            tokstart++;
            tokstart = copy_until_char_or_newline(tokstart, format_str, 0x20);
            aspect_d = (uint32_t)strtol(format_str, (char **)NULL, 10);
            if (PRINT_HEADER) {
                fprintf(stderr, "aspect_n = %d\n", aspect_n);
                fprintf(stderr, "aspect_d = %d\n", aspect_d);
            }
            break;
        default:
            /* Unknown section: skip it */
            while (*tokstart != 0x20 && *tokstart != '\0') tokstart++;
            break;
        }
    }

    /* Check that we did not try to parse further the end of the header string */
    assert(fresult + strlen(fresult) == tokstart);

    /*check if required parameters were read*/
    if (width == 0) {
        fprintf(cfg->error_log_file, "width not found in y4m header\n");
        return EB_ErrorBadParameter;
    }
    if (height == 0) {
        fprintf(cfg->error_log_file, "height not found in y4m header\n");
        return EB_ErrorBadParameter;
    }
    if (fr_n == 0 || fr_d == 0) {
        fprintf(cfg->error_log_file, "frame rate not found in y4m header\n");
        return EB_ErrorBadParameter;
    }

    /* Assign parameters to cfg */
    cfg->source_width           = width;
    cfg->source_height          = height;
    cfg->frame_rate_numerator   = fr_n;
    cfg->frame_rate_denominator = fr_d;
    cfg->frame_rate             = fr_n / fr_d;
    cfg->encoder_bit_depth      = bitdepth;
    /* TODO: when implemented, need to set input bit depth
        (instead of the encoder bit depth) and chroma format */

    return EB_ErrorNone;
}

/* read next line which contains the "FRAME" delimiter */
int32_t read_y4m_frame_delimiter(EbConfig *cfg) {
    unsigned char buffer_y4m_header[10];
    char *        fresult;

    fresult = fgets((char *)buffer_y4m_header, sizeof(buffer_y4m_header), cfg->input_file);

    if (fresult == NULL) {
        assert(feof(cfg->input_file));
        return EB_ErrorNone;
    }

    if (EB_STRCMP((const char *)buffer_y4m_header, "FRAME\n") !=0
        && EB_STRCMP((const char *)buffer_y4m_header, "FRAME\r\n") != 0) {
        fprintf(cfg->error_log_file, "Failed to read proper y4m frame delimeter. Read broken.\n");
        return EB_ErrorBadParameter;
    }

    return EB_ErrorNone;
}

/* check if the input file is in YUV4MPEG2 (y4m) format */
EbBool check_if_y4m(EbConfig *cfg) {
    size_t len;
    char   buf[YUV4MPEG2_IND_SIZE + 1];
    len = fread(buf, YUV4MPEG2_IND_SIZE, 1, cfg->input_file);
    if (len != 1) return EB_FALSE;

    if ((cfg->input_file != stdin) && (!cfg->input_file_is_fifo)) {
        fseek(cfg->input_file, 0, SEEK_SET);
    } else {
        memcpy(cfg->y4m_buf, buf, YUV4MPEG2_IND_SIZE);
    }

    buf[YUV4MPEG2_IND_SIZE] = 0;
    return (EB_STRCMP(buf, "YUV4MPEG2") == 0);
}
