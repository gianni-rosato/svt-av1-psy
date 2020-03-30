/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

// Command line argument parsing

/***************************************
 * Includes
 ***************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "EbDecParamParser.h"

static int parse_name(const char *value, const char *const *names) {
    for (int i = 0; names[i]; i++)
        if (!strcmp(value, names[i])) return i;
    return -1;
}

static void set_skip_frame(const char *value, EbSvtAv1DecConfiguration *cfg) {
    cfg->skip_frames = strtoul(value, NULL, 0);
};
static void set_limit_frame(const char *value, EbSvtAv1DecConfiguration *cfg) {
    cfg->frames_to_be_decoded = strtoul(value, NULL, 0);
};
static void set_bit_depth(const char *value, EbSvtAv1DecConfiguration *cfg) {
    cfg->max_bit_depth = strtoul(value, NULL, 0);
};
static void set_decoder_16bit_pipeline(const char *value, EbSvtAv1DecConfiguration *cfg) {
    cfg->is_16bit_pipeline = (EbBool)strtoul(value, NULL, 0);
    if (cfg->is_16bit_pipeline != 1 && cfg->is_16bit_pipeline != 0) {
        fprintf(stderr,
            "Warning : Invalid value for is_16bit_pipeline, setting value to 0. \n");
        cfg->is_16bit_pipeline = 0;
    }
};
static void set_pic_width(const char *value, EbSvtAv1DecConfiguration *cfg) {
    cfg->max_picture_width = strtoul(value, NULL, 0);
};
static void set_pic_height(const char *value, EbSvtAv1DecConfiguration *cfg) {
    cfg->max_picture_height = strtoul(value, NULL, 0);
};
static void set_colour_space(const char *value, EbSvtAv1DecConfiguration *cfg) {
    cfg->max_color_format = parse_name(value, csp_names);
};
static void set_num_thread(const char *value, EbSvtAv1DecConfiguration *cfg) {
    cfg->threads = strtoul(value, NULL, 0);
};
static void set_num_pframes(const char *value, EbSvtAv1DecConfiguration *cfg) {
    cfg->num_p_frames = strtoul(value, NULL, 0);
    if (cfg->num_p_frames != 1) {
        fprintf(
            stderr,
            "Warning : Multi frame parallelism not supported. Setting parallel frames to 1. \n");
        cfg->num_p_frames = 1;
    }
};

/**********************************
  * Config Entry Array
  **********************************/
ConfigEntry config_entry[] = {
    // Decoder settings
    {SKIP_FRAME_TOKEN, "SkipFrame", 1, set_skip_frame},
    {LIMIT_FRAME_TOKEN, "LimitFrame", 1, set_limit_frame},
    // Picture properties
    {BIT_DEPTH_TOKEN, "InputBitDepth", 1, set_bit_depth},
    {DECODER_16BIT_PIPELINE, "Decoder16BitPipeline", 1, set_decoder_16bit_pipeline},
    {PIC_WIDTH_TOKEN, "PictureWidth", 1, set_pic_width},
    {PIC_HEIGHT_TOKEN, "PictureHeight", 1, set_pic_height},
    {COLOUR_SPACE_TOKEN, "InputColourSpace", 1, set_colour_space},
    {THREADS_TOKEN, "ThreadCount", 1, set_num_thread},
    {FRAME_PLL_TOKEN, "PllFrameCount", 1, set_num_pframes},
    // Termination
    {NULL, NULL, 0, NULL}};

static void show_help() {
#define H0 printf
    H0(" Options : \n");
    H0( " -help                     Show usage options and exit \n");
    H0( " -i <arg>                  Input file name \n");
    H0( " -o <arg>                  Output file name \n");
    H0( " -skip <arg>               Skip the first n input frames \n");
    H0( " -limit <arg>              Stop decoding after n frames \n");
    H0( " -bit-depth <arg>          Input bitdepth. [8, 10] \n");
    H0( " -w <arg>                  Input picture width \n");
    H0( " -h <arg>                  Input picture height \n");
    H0( " -colour-space <arg>       Input picture colour space. [400, 420, 422, 444]\n");
    H0( " -threads <arg>            Number of threads to be launched \n");
    H0( " -parallel-frames <arg>    Number of frames to be processed in parallel \n");
    H0( " -md5                      MD5 support flag \n");
    H0( " -fps-frm                  Show fps after each frame decoded\n");
    H0( " -fps-summary              Show fps summary");
    H0( " -skip-film-grain          Disable Film Grain");
    H0( " -16bit-pipeline           Enable 16b pipeline. [1 - enable, 0 - disable]");

    exit(1);
}

EbErrorType read_command_line(int32_t argc, char *const argv[], EbSvtAv1DecConfiguration *configs,
                              CliInput *cli, ObuDecInputContext *obu_ctx) {
    char *  cmd_copy[MAX_NUM_TOKENS]       = {NULL};
    char *  config_strings[MAX_NUM_TOKENS] = {NULL};
    int32_t cmd_token_cnt                  = 0;
    int     token_index                    = 0;

    cli->skip_film_grain = configs->skip_film_grain;

    for (token_index = 1; token_index < argc; token_index++, cmd_token_cnt++) {
        if (argv[token_index][0] == '-') {
            cmd_copy[cmd_token_cnt] = argv[token_index];
            if (argv[token_index + 1] != NULL && (argv[token_index + 1][0] != '-'))
                config_strings[cmd_token_cnt] = argv[++token_index];
        } else {
            fprintf(stderr, " Invalid CLI: %s \n", argv[token_index]);
            return EB_ErrorBadParameter;
        }
    }

    token_index = 0;
    // Parse command line for tokens
    while (token_index < cmd_token_cnt) {
        if (cmd_copy[token_index] != NULL) {
            if (EB_STRCMP(cmd_copy[token_index], INPUT_FILE_TOKEN) == 0) {
                FILE *fin = NULL;
                FOPEN(fin, config_strings[token_index], "rb");
                if (!fin) {
                    fprintf(stderr, "Invalid input file \n");
                    return EB_ErrorBadParameter;
                } else {
                    cli->in_file     = fin;
                    cli->in_filename = config_strings[token_index];
                }
            } else if (EB_STRCMP(cmd_copy[token_index], OUTPUT_FILE_TOKEN) == 0) {
                FILE *fout = NULL;
                FOPEN(fout, config_strings[token_index], "wb");
                if (!fout) {
                    fprintf(stderr, "Invalid output file \n");
                    return EB_ErrorBadParameter;
                } else {
                    cli->out_file     = fout;
                    cli->out_filename = config_strings[token_index];
                }
            } else if (EB_STRCMP(cmd_copy[token_index], MD5_SUPPORT_TOKEN) == 0)
                cli->enable_md5 = 1;
            else if (EB_STRCMP(cmd_copy[token_index], FPS_FRM_TOKEN) == 0)
                cli->fps_frm = 1;
            else if (EB_STRCMP(cmd_copy[token_index], FPS_SUMMARY_TOKEN) == 0)
                cli->fps_summary = 1;
            else if (EB_STRCMP(cmd_copy[token_index], FILM_GRAIN_TOKEN) == 0)
                cli->skip_film_grain = 1;
            else if (EB_STRCMP(cmd_copy[token_index], ANNEX_B_TOKEN) == 0)
                obu_ctx->is_annexb = 1;
            else if (EB_STRCMP(cmd_copy[token_index], HELP_TOKEN) == 0)
                show_help();
            else {
                int temp_ind = 0;
                int cli_read = 0;
                while (config_entry[temp_ind].name != NULL) {
                    if (EB_STRCMP(cmd_copy[token_index], config_entry[temp_ind].token) == 0) {
                        if (config_strings[token_index] == NULL) {
                            if (config_entry[temp_ind].value_required == 1) {
                                fprintf(stderr, "Invalid CLI option: %s \n", cmd_copy[token_index]);
                                return EB_ErrorBadParameter;
                            } else
                                config_strings[token_index] = "1";
                        }
                        (*config_entry[temp_ind].scf)(config_strings[token_index], configs);
                        cli_read = 1;
                    }
                    temp_ind++;
                }
                if (!cli_read) {
                    fprintf(stderr, "Invalid CLI option: %s \n", cmd_copy[token_index]);
                    return EB_ErrorBadParameter;
                }
            }
        } else
            break;
        token_index++;
    }

    if (!cli->in_file) {
        fprintf(stderr, "Input file not specified. \n");
        show_help();
        return EB_ErrorBadParameter;
    }

    cli->fmt = configs->max_color_format;
    configs->skip_film_grain = cli->skip_film_grain;

    if (file_is_ivf(cli)) {
        cli->in_file_type = FILE_TYPE_IVF;
        assert(0 == obu_ctx->is_annexb);
    } else if (file_is_obu(cli, obu_ctx))
        cli->in_file_type = FILE_TYPE_OBU;
    else {
        fprintf(stderr, "Unsupported input file format. \n");
        return EB_ErrorBadParameter;
    }
    configs->max_picture_height = cli->height;
    configs->max_picture_width = cli->width;
    return EB_ErrorNone;
}
