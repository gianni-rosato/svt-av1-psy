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

#include "EbDecParamParser.h"

static int parse_name(const char* value, const char* const* names) {
    for (int i = 0; names[i]; i++)
        if (!strcmp(value, names[i]))
            return i;
    return -1;
}

static void set_skip_frame(const char *value, EbSvtAv1DecConfiguration *cfg) { cfg->skip_frames = strtoul(value, NULL, 0); };
static void set_limit_frame(const char *value, EbSvtAv1DecConfiguration *cfg) { cfg->frames_to_be_decoded = strtoul(value, NULL, 0); };
static void set_bit_depth(const char *value, EbSvtAv1DecConfiguration *cfg) { cfg->max_bit_depth = strtoul(value, NULL, 0); };
static void set_pic_width(const char *value, EbSvtAv1DecConfiguration *cfg) { cfg->max_picture_width = strtoul(value, NULL, 0); };
static void set_pic_height(const char *value, EbSvtAv1DecConfiguration *cfg) { cfg->max_picture_height = strtoul(value, NULL, 0); };
static void set_colour_space(const char *value, EbSvtAv1DecConfiguration *cfg) { cfg->max_color_format = parse_name(value, csp_names); };

 /**********************************
  * Config Entry Array
  **********************************/
ConfigEntry config_entry[] = {
    // Decoder settings
    { SKIP_FRAME_TOKEN, "SkipFrame", 1, set_skip_frame },
    { LIMIT_FRAME_TOKEN, "LimitFrame", 1, set_limit_frame },
    // Picture properties
    { BIT_DEPTH_TOKEN,"InputBitDepth", 1, set_bit_depth },
    { PIC_WIDTH_TOKEN, "PictureWidth", 1, set_pic_width},
    { PIC_HEIGHT_TOKEN, "PictureHeight", 1, set_pic_height},
    { COLOUR_SPACE_TOKEN,"InputColourSpace", 1, set_colour_space},
    // Termination
    { NULL, NULL, 0, NULL}
};

static void showHelp()
{
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
    H0( " -md5                      MD5 support flag \n");

    exit(1);
}

EbErrorType read_command_line(int32_t argc, char *const argv[],
                              EbSvtAv1DecConfiguration *configs,
                              CLInput *cli)
{
    char    *cmd_copy[MAX_NUM_TOKENS] = { NULL };
    char    *config_strings[MAX_NUM_TOKENS] = { NULL };
    int32_t cmd_token_cnt = 0;
    int token_index = 0;

    for (token_index = 1; token_index < argc; token_index++, cmd_token_cnt++) {
        if (argv[token_index][0] == '-') {
            cmd_copy[cmd_token_cnt] = argv[token_index];
            if (argv[token_index + 1] != NULL && (argv[token_index + 1][0] != '-'))
                config_strings[cmd_token_cnt] = argv[++token_index];
        }
        else {
            printf(" Invalid CLI: %s \n", argv[token_index]);
            return EB_ErrorBadParameter;
        }
    }

    token_index = 0;
    // Parse command line for tokens
    while (token_index < cmd_token_cnt) {
        if (cmd_copy[token_index] != NULL) {
            if (EB_STRCMP(cmd_copy[token_index], INPUT_FILE_TOKEN) == 0) {
                FILE * fin = NULL;
                FOPEN(fin, config_strings[token_index], "rb");
                if (!fin) {
                    printf("Invalid input file \n");
                    return EB_ErrorBadParameter;
                }
                else {
                    cli->inFile = fin;
                    cli->inFilename = config_strings[token_index];
                    if (file_is_ivf(cli))
                        cli->inFileType = FILE_TYPE_IVF;
                    else {
                        printf("Unsupported input file format. \n");
                        return EB_ErrorBadParameter;
                    }
                }
            }
            else if (EB_STRCMP(cmd_copy[token_index], OUTPUT_FILE_TOKEN) == 0) {
                FILE * fout = NULL;
                FOPEN(fout, config_strings[token_index], "wb");
                if (!fout) {
                    printf("Invalid output file \n");
                    return EB_ErrorBadParameter;
                }
                else {
                    cli->outFile = fout;
                    cli->outFilename = config_strings[token_index];
                }
            }
            else if (EB_STRCMP(cmd_copy[token_index], MD5_SUPPORT_TOKEN) == 0)
                cli->enable_md5 = 1;
            else if (EB_STRCMP(cmd_copy[token_index], HELP_TOKEN) == 0)
                showHelp();
            else {
                int temp_ind = 0;
                int cli_read = 0;
                while (config_entry[temp_ind].name != NULL)
                {
                    if (EB_STRCMP(cmd_copy[token_index], config_entry[temp_ind].token) == 0) {
                        if (config_strings[token_index] == NULL) {
                            if (config_entry[temp_ind].valueRequired == 1) {
                                printf("Invalid CLI option: %s \n", cmd_copy[token_index]);
                                return EB_ErrorBadParameter;
                            }
                            else
                                config_strings[token_index] = "1";
                        }
                        (*config_entry[temp_ind].scf)(config_strings[token_index], configs);
                        cli_read = 1;
                    }
                    temp_ind++;
                }
                if (!cli_read) {
                    printf("Invalid CLI option: %s \n", cmd_copy[token_index]);
                    return EB_ErrorBadParameter;
                }
            }
        }
        else
            break;
        token_index++;
    }

    if (!cli->inFile || !cli->outFile) {
        printf("Input/output file not specified. \n");
        showHelp();
        return EB_ErrorBadParameter;
    }

    cli->bit_depth = configs->max_bit_depth;
    cli->fmt = configs->max_color_format;

    if(cli->height != configs->max_picture_height)
        configs->max_picture_height = cli->height;
    if (cli->width != configs->max_picture_width)
        configs->max_picture_width = cli->width;

    return EB_ErrorNone;
}
