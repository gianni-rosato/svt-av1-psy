/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

/*
* Copyright (c) 2016, Alliance for Open Media. All rights reserved
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at www.aomedia.org/license/software. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at www.aomedia.org/license/patent.
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#if defined(_WIN32)
/*
  * Win32 specific includes
  */
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#else
/*
  * POSIX specific includes
  */
#include <sys/time.h>

/* timersub is not provided by msys at this time. */
#ifndef timersub
#define timersub(a, b, result)                           \
    do {                                                 \
        (result)->tv_sec  = (a)->tv_sec - (b)->tv_sec;   \
        (result)->tv_usec = (a)->tv_usec - (b)->tv_usec; \
        if ((result)->tv_usec < 0) {                     \
            --(result)->tv_sec;                          \
            (result)->tv_usec += 1000000;                \
        }                                                \
    } while (0)
#endif
#endif

struct EbDecTimer {
#if defined(_WIN32)
    LARGE_INTEGER begin, end;
#else
    struct timeval begin, end;
#endif
};

void    dec_timer_start(struct EbDecTimer *t);
void    dec_timer_mark(struct EbDecTimer *t);
int64_t dec_timer_elapsed(struct EbDecTimer *t);
