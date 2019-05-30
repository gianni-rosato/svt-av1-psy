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

#include "EbSvtAv1Dec.h"
#include "EbFileUtils.h"

#ifdef _MSC_VER
#define FOPEN(f,s,m) fopen_s(&f,s,m)
#else
#define FOPEN(f,s,m) f=fopen(s,m)
#endif

 /**********************************
  * CLI options
  **********************************/
#define HELP_TOKEN                      "-help"
#define COMMAND_LINE_MAX_SIZE           2048
#define INPUT_FILE_TOKEN                "-i"
#define OUTPUT_FILE_TOKEN               "-o"
#define SKIP_FRAME_TOKEN                "-skip"
#define LIMIT_FRAME_TOKEN               "-limit"
#define BIT_DEPTH_TOKEN                 "-bit-depth"
#define PIC_WIDTH_TOKEN                 "-w"
#define PIC_HEIGHT_TOKEN                "-h"
#define COLOUR_SPACE_TOKEN              "-colour-space"
#define MD5_SUPPORT_TOKEN               "-md5"
#define MAX_NUM_TOKENS 200

#define EB_STRCMP(target,token) \
    strcmp(target,token)

/**********************************
 * Config Entry Struct
 **********************************/
typedef struct ConfigEntry {
    const char *token;
    const char *name;
    EbBool valueRequired;
    void(*scf)(const char *, EbSvtAv1DecConfiguration *);
} ConfigEntry;

EbErrorType read_command_line(int32_t argc, char *const argv[], EbSvtAv1DecConfiguration *configs, CLInput *cli);
