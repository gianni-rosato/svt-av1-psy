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

#include "EbDecTime.h"

void dec_timer_start(struct EbDecTimer *t) {
#if defined(_WIN32)
    QueryPerformanceCounter(&t->begin);
#else
    gettimeofday(&t->begin, NULL);
#endif
}

void dec_timer_mark(struct EbDecTimer *t) {
#if defined(_WIN32)
    QueryPerformanceCounter(&t->end);
#else
    gettimeofday(&t->end, NULL);
#endif
}

int64_t dec_timer_elapsed(struct EbDecTimer *t) {
#if defined(_WIN32)
    LARGE_INTEGER freq, diff;

    diff.QuadPart = t->end.QuadPart - t->begin.QuadPart;

    QueryPerformanceFrequency(&freq);
    return diff.QuadPart * 1000000 / freq.QuadPart;
#else
    struct timeval diff;

    timersub(&t->end, &t->begin, &diff);
    return ((int64_t)diff.tv_sec) * 1000000 + diff.tv_usec;
#endif
}
