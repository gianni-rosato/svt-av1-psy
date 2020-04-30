/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/
#include "EbLog.h"
//for getenv and fopen on windows
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

static SvtLogLevel g_log_level;
static FILE*       g_log_file;

static void svt_log_set_level(SvtLogLevel level) { g_log_level = level; }

static void svt_log_set_log_file(const char* file) {
    if (file) g_log_file = fopen(file, "w+");
}

void svt_log_init() {
    const char* log   = getenv("SVT_LOG");
    SvtLogLevel level = SVT_LOG_INFO;
    if (log) level = (SvtLogLevel)atoi(log);
    svt_log_set_level(level);

    if (!g_log_file) {
        const char* file = getenv("SVT_LOG_FILE");
        svt_log_set_log_file(file);
    }
}

static const char* log_level_str(SvtLogLevel level) {
    switch (level) {
    case SVT_LOG_FATAL: return "fatal";
    case SVT_LOG_ERROR: return "error";
    case SVT_LOG_WARN: return "warn";
    case SVT_LOG_INFO: return "info";
    case SVT_LOG_DEBUG: return "debug";
    default: return "unknown";
    }
}

void svt_log(SvtLogLevel level, const char* tag, const char* format, ...) {
    if (level > g_log_level) return;

    if (!g_log_file) g_log_file = stderr;

    if (tag) fprintf(g_log_file, "%s[%s]: ", tag, log_level_str(level));

    va_list args;
    va_start(args, format);
    vfprintf(g_log_file, format, args);
    va_end(args);
}
