/*
 * Copyright(c) 2019 Intel Corporation
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

// EbUnitTestUtility.c
//  -Unit test utility
//      -Randomize Data
//      -Compare Data
//      -Log Data
//      -Compute Overall Elapsed Millisecond

/***************************************
 * Includes
 ***************************************/
#include <time.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include "EbUnitTest.h"
#include "EbUnitTestUtility.h"
//#include "EbTime.h"

/***************************************
 * Randomize Data
 ***************************************/
void svt_buf_random_void(void *const buf, const uint32_t sizeBuf) {
    uint8_t *const buffer = (uint8_t *)buf;

    for (uint32_t i = 0; i < sizeBuf; i++)
        buffer[i] = (uint8_t)(rand() % 256);
}

void svt_buf_random_s8(int8_t *const buf, const uint32_t sizeBuf) {
    svt_buf_random_void(buf, sizeBuf);
}

void svt_buf_random_s8_with_max(int8_t *const buf, const uint32_t sizeBuf,
                                const uint8_t max) {
    for (uint32_t i = 0; i < sizeBuf; i++) {
        const uint32_t val = rand() % (max + 1);
        const uint32_t sign = rand() & 1;
        buf[i] = (int8_t)(sign ? val : 0 - val);
    }
}

void svt_buf_random_u8(uint8_t *const buf, const uint32_t sizeBuf) {
    svt_buf_random_void(buf, sizeBuf);
}

void svt_buf_random_u8_to_0_or_255(uint8_t *const buf, const uint32_t sizeBuf) {
    for (uint32_t i = 0; i < sizeBuf; i++)
        buf[i] = (rand() > (RAND_MAX >> 1)) ? 255 : 0;
}

void svt_buf_random_u8_to_255(uint8_t *const buf, const uint32_t sizeBuf) {
    memset(buf, 255, sizeBuf);
}

void svt_buf_random_u8_to_large(uint8_t *const buf, const uint32_t sizeBuf) {
    for (uint32_t i = 0; i < sizeBuf; i++)
        buf[i] = 255 - (rand() % 10);
}

void svt_buf_random_u8_to_near_value(uint8_t *const buf, const uint32_t sizeBuf,
                                     const uint8_t val, const uint32_t range) {
    for (uint32_t i = 0; i < sizeBuf; i++) {
        int32_t v = val;
        v += (rand() % range);
        v -= (rand() % range);
        if (v > 255)
            v = 255;
        if (v < 0)
            v = 0;
        buf[i] = (uint8_t)v;
    }
}

void svt_buf_random_u8_to_small(uint8_t *const buf, const uint32_t sizeBuf) {
    for (uint32_t i = 0; i < sizeBuf; i++)
        buf[i] = (rand() % 10);
}

void svt_buf_random_u8_to_small_or_large(uint8_t *const buf,
                                         const uint32_t sizeBuf) {
    for (uint32_t i = 0; i < sizeBuf; i++)
        buf[i] =
            (rand() > (RAND_MAX >> 1)) ? (rand() % 10) : (255 - (rand() % 10));
}

void svt_buf_random_u8_with_max(uint8_t *const buf, const uint32_t sizeBuf,
                                const uint8_t max) {
    for (uint32_t i = 0; i < sizeBuf; i++)
        buf[i] = (uint8_t)((uint32_t)rand() % (max + 1));
}

void svt_buf_random_s16(int16_t *const buf, const uint32_t sizeBuf) {
    svt_buf_random_void(buf, sizeof(*buf) * sizeBuf);
}

void svt_buf_random_s16_to_bd(int16_t *const buf, const uint32_t sizeBuf,
                              const uint32_t bd) {
    const int16_t max = (1 << bd) - 1;

    for (uint32_t i = 0; i < sizeBuf; i++)
        buf[i] = (rand() & 1) ? max : -max;
}

void svt_buf_random_s16_with_bd(int16_t *const buf, const uint32_t sizeBuf,
                                const uint32_t bd) {
    assert(bd >= 8);

    for (uint32_t i = 0; i < sizeBuf; i++) {
        const uint32_t val = rand() % (1 << bd);
        const uint32_t sign = rand() & 1;
        buf[i] = (int16_t)(sign ? val : 0 - val);
    }
}

void svt_buf_random_s16_with_max(int16_t *const buf, const uint32_t sizeBuf,
                                 const uint16_t max) {
    for (uint32_t i = 0; i < sizeBuf; i++) {
        const uint32_t val = rand() % (max + 1);
        const uint32_t sign = rand() & 1;
        buf[i] = (int16_t)(sign ? val : 0 - val);
    }
}

void svt_buf_random_u16(uint16_t *const buf, const uint32_t sizeBuf) {
    svt_buf_random_void(buf, sizeof(*buf) * sizeBuf);
}

void svt_buf_random_u16_to_0_or_1023(uint16_t *const buf,
                                     const uint32_t sizeBuf) {
    for (uint32_t i = 0; i < sizeBuf; i++)
        buf[i] = (rand() > (RAND_MAX >> 1)) ? 1023 : 0;
}

void svt_buf_random_u16_to_0_or_bd(uint16_t *const buf, const uint32_t sizeBuf,
                                   const uint32_t bd) {
    const uint16_t max = (uint16_t)((1 << bd) - 1);

    for (uint32_t i = 0; i < sizeBuf; i++)
        buf[i] = (rand() > (RAND_MAX >> 1)) ? max : 0;
}

void svt_buf_random_u16_to_bd(uint16_t *const buf, const uint32_t sizeBuf,
                              const uint32_t bd) {
    const uint16_t max = (uint16_t)((1 << bd) - 1);

    for (uint32_t i = 0; i < sizeBuf; i++)
        buf[i] = max;
}

void svt_buf_random_u16_with_bd(uint16_t *const buf, const uint32_t sizeBuf,
                                const uint32_t bd) {
    assert(bd >= 8);

    for (uint32_t i = 0; i < sizeBuf; i++)
        buf[i] = (uint16_t)((uint32_t)rand() % (1 << bd));
}

void svt_buf_random_u16_with_max(uint16_t *const buf, const uint32_t sizeBuf,
                                 const uint32_t max) {
    for (uint32_t i = 0; i < sizeBuf; i++)
        buf[i] = (uint16_t)((uint32_t)rand() % (max + 1));
}

void svt_buf_random_s32(int32_t *const buf, const uint32_t sizeBuf) {
    svt_buf_random_void(buf, sizeof(*buf) * sizeBuf);
}

void svt_buf_random_s32_with_max(int32_t *const buf, const uint32_t sizeBuf,
                                 const int32_t max) {
    for (uint32_t i = 0; i < sizeBuf; i++)
        buf[i] = (rand() & 1) ? max : -max;
}

void svt_buf_random_u32(uint32_t *const buf, const uint32_t sizeBuf) {
    svt_buf_random_void(buf, sizeof(*buf) * sizeBuf);
}

void svt_buf_random_u32_with_max(uint32_t *const buf, const uint32_t sizeBuf,
                                 const uint32_t max) {
    svt_buf_random_void(buf, sizeof(*buf) * sizeBuf);

    for (uint32_t i = 0; i < sizeBuf; i++)
        buf[i] %= (max + 1);
}

void svt_buf_random_s64(int64_t *const buf, const uint32_t sizeBuf) {
    svt_buf_random_void(buf, sizeof(*buf) * sizeBuf);
}

void svt_buf_random_u64(uint64_t *const buf, const uint32_t sizeBuf) {
    svt_buf_random_void(buf, sizeof(*buf) * sizeBuf);
}

void svt_buf_random_ptrdiff_with_max(ptrdiff_t *const buf,
                                     const uint32_t sizeBuf,
                                     const ptrdiff_t max) {
    svt_buf_random_void(buf, sizeof(*buf) * sizeBuf);

    for (uint32_t i = 0; i < sizeBuf; i++)
        buf[i] %= (max + 1);
}

uint32_t svt_create_random_aligned_stride(const uint32_t width,
                                          const uint32_t align) {
    uint32_t stride;

    svt_buf_random_u32_with_max(&stride, 1, 100);
    stride += width + align;
    stride -= stride % align;
    return stride;
}

/***************************************
 * Compare Data
 ***************************************/
#define BUF_COMPARE_FUNCTION(name, type)                         \
    EbBool name(const type *const buf1,                          \
                const type *const buf2,                          \
                const size_t bufSize) {                          \
        EbBool result = EB_TRUE;                                 \
                                                                 \
        for (uint32_t i = 0; i < bufSize; i++) {                 \
            if (buf1[i] != buf2[i]) {                            \
                printf("\nbuf1[%3d] = 0x%8x\tbuf2[%3d] = 0x%8x", \
                       i,                                        \
                       buf1[i],                                  \
                       i,                                        \
                       buf2[i]);                                 \
                result = EB_FALSE;                               \
            }                                                    \
        }                                                        \
                                                                 \
        if (!result)                                             \
            printf("\n\n");                                      \
                                                                 \
        return result;                                           \
    }

BUF_COMPARE_FUNCTION(svt_buf_compare_u8, uint8_t)
BUF_COMPARE_FUNCTION(svt_buf_compare_s8, int8_t)
BUF_COMPARE_FUNCTION(svt_buf_compare_u16, uint16_t)
BUF_COMPARE_FUNCTION(svt_buf_compare_s16, int16_t)
BUF_COMPARE_FUNCTION(svt_buf_compare_u32, uint32_t)
BUF_COMPARE_FUNCTION(svt_buf_compare_s32, int32_t)

#define IMAGE_COMPARE_FUNCTION(name, type)                                 \
    EbBool name(const type *buf1,                                          \
                const type *buf2,                                          \
                const size_t width,                                        \
                const size_t height,                                       \
                const size_t stride) {                                     \
        EbBool result = EB_TRUE;                                           \
                                                                           \
        buf1 -= stride;                                                    \
        buf2 -= stride;                                                    \
                                                                           \
        for (int32_t i = -1; i < (int32_t)height + 1;                      \
             i++, buf1 += stride, buf2 += stride) {                        \
            for (uint32_t j = 0; j < stride; j++) {                        \
                if (buf1[j] != buf2[j]) {                                  \
                    printf("\nbuf1[%3d][%3d] = 0x%16" PRIx64               \
                           "\tbuf2[%3d][%3d] = 0x%16" PRIx64,              \
                           i,                                              \
                           j,                                              \
                           (uint64_t)buf1[j],                              \
                           i,                                              \
                           j,                                              \
                           (uint64_t)buf2[j]);                             \
                    if ((i < 0) || (i >= (int32_t)height) || (j >= width)) \
                        printf(" (outside image)");                        \
                    result = EB_FALSE;                                     \
                }                                                          \
            }                                                              \
        }                                                                  \
                                                                           \
        if (!result)                                                       \
            printf("\n\n");                                                \
        return result;                                                     \
    }

IMAGE_COMPARE_FUNCTION(svt_image_compare_u8, uint8_t)
IMAGE_COMPARE_FUNCTION(svt_image_compare_s8, int8_t)
IMAGE_COMPARE_FUNCTION(svt_image_compare_u16, uint16_t)
IMAGE_COMPARE_FUNCTION(svt_image_compare_s16, int16_t)
IMAGE_COMPARE_FUNCTION(svt_image_compare_u32, uint32_t)
IMAGE_COMPARE_FUNCTION(svt_image_compare_s32, int32_t)
IMAGE_COMPARE_FUNCTION(svt_image_compare_u64, uint64_t)

/***************************************
 * Log Data
 ***************************************/
void svt_unit_test_log_s8(const char *const nameBuf, const int8_t *const buf,
                          const uint32_t sizeBuf) {
    printf("%16s = ", nameBuf);

    if (sizeBuf == 1)
        printf("%d\n", buf[0]);
    else {
        for (uint32_t i = 0; i < sizeBuf; i++) {
            printf("%4d,", buf[i]);
            if (!((i + 1) % 32))
                printf("\n  ");
        }
        printf("\n");
    }
}

void svt_unit_test_log_u8(const char *const nameBuf, const uint8_t *const buf,
                          const uint32_t sizeBuf) {
    printf("%16s = ", nameBuf);

    if (sizeBuf == 1)
        printf("%u\n", buf[0]);
    else {
        for (uint32_t i = 0; i < sizeBuf; i++) {
            printf("%4u,", buf[i]);
            if (!((i + 1) % 32))
                printf("\n  ");
        }
        printf("\n");
    }
}

void svt_unit_test_log_s16(const char *const nameBuf, const int16_t *const buf,
                           const uint32_t sizeBuf) {
    printf("%16s = ", nameBuf);

    if (sizeBuf == 1)
        printf("%d\n", buf[0]);
    else {
        for (uint32_t i = 0; i < sizeBuf; i++) {
            printf("%10d,", buf[i]);
            if (!((i + 1) % 32))
                printf("\n  ");
        }
        printf("\n");
    }
}

void svt_unit_test_log_u16(const char *const nameBuf, const uint16_t *const buf,
                           const uint32_t sizeBuf) {
    printf("%16s = ", nameBuf);

    if (sizeBuf == 1)
        printf("%d\n", buf[0]);
    else {
        for (uint32_t i = 0; i < sizeBuf; i++) {
            printf("%10u,", buf[i]);
            if (!((i + 1) % 32))
                printf("\n  ");
        }
        printf("\n");
    }
}

void svt_unit_test_log_s32(const char *const nameBuf, const int32_t *const buf,
                           const uint32_t sizeBuf) {
    printf("%16s = ", nameBuf);

    if (sizeBuf == 1)
        printf("%d\n", buf[0]);
    else {
        for (uint32_t i = 0; i < sizeBuf; i++) {
            printf("%10d,", buf[i]);
            if (!((i + 1) % 32))
                printf("\n  ");
        }
        printf("\n");
    }
}

void svt_unit_test_log_u32(const char *const nameBuf, const uint32_t *const buf,
                           const uint32_t sizeBuf) {
    printf("%16s = ", nameBuf);

    if (sizeBuf == 1)
        printf("%u\n", buf[0]);
    else {
        for (uint32_t i = 0; i < sizeBuf; i++) {
            printf("%10u,", buf[i]);
            if (!((i + 1) % 32))
                printf("\n  ");
        }
        printf("\n");
    }
}

void svt_unit_test_log_s64(const char *const nameBuf, const int64_t *const buf,
                           const uint32_t sizeBuf) {
    printf("%16s = ", nameBuf);

    if (sizeBuf == 1)
        printf("%" PRIx64 "\n", buf[0]);
    else {
        for (uint32_t i = 0; i < sizeBuf; i++) {
            printf("%10" PRIx64 ",", buf[i]);
            if (!((i + 1) % 32))
                printf("\n  ");
        }
        printf("\n");
    }
}

void svt_unit_test_log_u64(const char *const nameBuf, const uint64_t *const buf,
                           const uint32_t sizeBuf) {
    printf("%16s = ", nameBuf);

    if (sizeBuf == 1)
        printf("%" PRIu64 "\n", buf[0]);
    else {
        for (uint32_t i = 0; i < sizeBuf; i++) {
            printf("%10" PRIu64 ",", buf[i]);
            if (!((i + 1) % 32))
                printf("\n  ");
        }
        printf("\n");
    }
}

void svt_unit_test_log_ptrdiff(const char *const nameBuf,
                               const ptrdiff_t *const buf,
                               const uint32_t sizeBuf) {
    printf("%16s = ", nameBuf);

    if (sizeBuf == 1)
        printf("%td\n", buf[0]);
    else {
        for (uint32_t i = 0; i < sizeBuf; i++) {
            printf("%10td,", buf[i]);
            if (!((i + 1) % 32))
                printf("\n  ");
        }
        printf("\n");
    }
}

void svt_unit_test_log_bool(const char *const nameBuf, const EbBool *const buf,
                            const uint32_t sizeBuf) {
    printf("%16s = ", nameBuf);

    if (sizeBuf == 1)
        printf("%u\n", buf[0]);
    else {
        for (uint32_t i = 0; i < sizeBuf; i++) {
            printf("%10u,", buf[i]);
            if (!((i + 1) % 32))
                printf("\n  ");
        }
        printf("\n");
    }
}

void svt_unit_test_log_double(const char *const nameBuf, const double *const buf,
                              const uint32_t sizeBuf) {
    printf("%16s = ", nameBuf);

    if (sizeBuf == 1)
        printf("%f\n", buf[0]);
    else {
        for (uint32_t i = 0; i < sizeBuf; i++) {
            printf("%10f,", buf[i]);
            if (!((i + 1) % 32))
                printf("\n  ");
        }
        printf("\n");
    }
}

#define LOG_IMAGE_FUNCTION(name, type)                \
    void name(const char *const nameBuf,              \
              const type *const buf,                  \
              const uint32_t width,                   \
              const uint32_t height,                  \
              const ptrdiff_t stride) {               \
        printf("%16s = ", nameBuf);                   \
                                                      \
        for (uint32_t i = 0; i < height; i++) {       \
            for (uint32_t j = 0; j < width; j++)      \
                printf("%10d,", buf[i * stride + j]); \
            printf("\n  ");                           \
        }                                             \
                                                      \
        printf("\n");                                 \
    }

LOG_IMAGE_FUNCTION(svt_unit_test_log_image_u8, uint8_t)
LOG_IMAGE_FUNCTION(svt_unit_test_log_image_s16, int16_t)
LOG_IMAGE_FUNCTION(svt_unit_test_log_image_u16, uint16_t)
LOG_IMAGE_FUNCTION(svt_unit_test_log_image_s32, int32_t)
LOG_IMAGE_FUNCTION(svt_unit_test_log_image_u32, uint32_t)
