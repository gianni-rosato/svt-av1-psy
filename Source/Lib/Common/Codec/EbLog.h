/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/
#ifndef EbLog_h
#define EbLog_h

#ifndef LOG_TAG
#define LOG_TAG "Svt"
#endif

typedef enum {
    SVT_LOG_ALL   = -1,
    SVT_LOG_FATAL = 0,
    SVT_LOG_ERROR = 1,
    SVT_LOG_WARN  = 2,
    SVT_LOG_INFO  = 3,
    SVT_LOG_DEBUG = 4,
} SvtLogLevel;

//define this to turn off all library log
//#define SVT_LOG_QUIET
#ifndef SVT_LOG_QUIET

//SVT_LOG will not output the prefix. you can contorl the output style.
#define SVT_LOG(format, ...) svt_log(SVT_LOG_ALL, NULL, format, ##__VA_ARGS__)

#define SVT_DEBUG(format, ...) svt_log(SVT_LOG_DEBUG, LOG_TAG, format, ##__VA_ARGS__)
#define SVT_INFO(format, ...) svt_log(SVT_LOG_INFO, LOG_TAG, format, ##__VA_ARGS__)
#define SVT_WARN(format, ...) svt_log(SVT_LOG_WARN, LOG_TAG, format, ##__VA_ARGS__)
#define SVT_ERROR(format, ...) svt_log(SVT_LOG_ERROR, LOG_TAG, format, ##__VA_ARGS__)
#define SVT_FATAL(format, ...) svt_log(SVT_LOG_FATAL, LOG_TAG, format, ##__VA_ARGS__)

#else

#define SVT_LOG(format, ...) \
    do {                     \
    } while (0)
#define SVT_DEBUG(format, ...) \
    do {                       \
    } while (0)
#define SVT_INFO(format, ...) \
    do {                      \
    } while (0)
#define SVT_WARN(format, ...) \
    do {                      \
    } while (0)
#define SVT_ERROR(format, ...) \
    do {                       \
    } while (0)
#define SVT_FATAL(format, ...) \
    do {                       \
    } while (0)

#endif //SVT_LOG_QUIET

void svt_log_init();
void svt_log(SvtLogLevel level, const char* tag, const char* format, ...);

#endif //EbLog_h
