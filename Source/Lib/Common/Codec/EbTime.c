/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef __USE_POSIX199309
#define __USE_POSIX199309
#endif

#include <time.h>

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

#include "EbTime.h"

#if defined(__GNUC__) && !defined(__clang__) && !defined(__ICC__)
__attribute__((optimize("unroll-loops")))
#endif

void eb_start_time(uint64_t *start_seconds, uint64_t *start_u_seconds) {
#ifdef _WIN32
    *start_seconds = (uint64_t)clock();
    (void)(*start_u_seconds);
#else
    struct timeval start;
    gettimeofday(&start, NULL);
    *start_seconds   = start.tv_sec;
    *start_u_seconds = start.tv_usec;
#endif
}

void eb_finish_time(uint64_t *finish_seconds, uint64_t *finish_u_seconds) {
#ifdef _WIN32
    *finish_seconds = (uint64_t)clock();
    (void)(*finish_u_seconds);
#else
    struct timeval finish;
    gettimeofday(&finish, NULL);
    *finish_seconds   = finish.tv_sec;
    *finish_u_seconds = finish.tv_usec;
#endif
}

void eb_compute_overall_elapsed_time(uint64_t start_seconds, uint64_t start_u_seconds,
                                     uint64_t finish_seconds, uint64_t finish_u_seconds,
                                     double *duration) {
#ifdef _WIN32
    //double  duration;
    *duration = (double)(finish_seconds - start_seconds) / CLOCKS_PER_SEC;
    (void)(start_u_seconds);
    (void)(finish_u_seconds);
#else
    long mtime, seconds, useconds;
    seconds   = finish_seconds - start_seconds;
    useconds  = finish_u_seconds - start_u_seconds;
    mtime     = ((seconds)*1000 + useconds / 1000.0) + 0.5;
    *duration = (double)mtime / 1000;
#endif
}

void eb_compute_overall_elapsed_time_ms(uint64_t start_seconds, uint64_t start_u_seconds,
                                        uint64_t finish_seconds, uint64_t finish_u_seconds,
                                        double *duration) {
#ifdef _WIN32
    //double  duration;
    *duration = (double)(finish_seconds - start_seconds);
    (void)(start_u_seconds);
    (void)(finish_u_seconds);
#else
    long mtime, seconds, useconds;
    seconds   = finish_seconds - start_seconds;
    useconds  = finish_u_seconds - start_u_seconds;
    mtime     = ((seconds)*1000 + useconds / 1000.0) + 0.5;
    *duration = (double)mtime;
#endif
}

void eb_sleep_ms(uint64_t milli_seconds) {
    if (milli_seconds) {
#ifdef _WIN32
        Sleep((DWORD)milli_seconds);
#else
        struct timespec req, rem;
        req.tv_sec = (int32_t)(milli_seconds / 1000);
        milli_seconds -= req.tv_sec * 1000;
        req.tv_nsec = milli_seconds * 1000000UL;
        nanosleep(&req, &rem);
#endif
    }
}

