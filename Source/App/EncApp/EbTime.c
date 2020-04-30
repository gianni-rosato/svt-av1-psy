/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdint.h>

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

void start_time(uint64_t *start_seconds, uint64_t *start_u_seconds) {
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

void finish_time(uint64_t *finish_seconds, uint64_t *finish_u_seconds) {
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

void compute_overall_elapsed_time(uint64_t start_seconds, uint64_t start_u_seconds,
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

static void sleep_ms(uint64_t milli_seconds) {
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

void injector(uint64_t processed_frame_count, uint32_t injector_frame_rate) {
#ifdef _WIN32
    static LARGE_INTEGER start_count; // this is the start time
    static LARGE_INTEGER counter_freq; // performance counter frequency
    LARGE_INTEGER        now_count; // this is the current time
#else
    uint64_t        currentTimesSeconds  = 0;
    uint64_t        currentTimesuSeconds = 0;
    static uint64_t startTimesSeconds;
    static uint64_t startTimesuSeconds;
#endif

    double injector_interval =
        (double)(1 << 16) /
        (double)injector_frame_rate; // 1.0 / injector frame rate (in this case, 1.0/encodRate)
    double         elapsed_time;
    double         predicted_time;
    int32_t        buffer_frames = 1; // How far ahead of time should we let it get
    int32_t        millisec_ahead;
    static int32_t first_time = 0;

    if (first_time == 0) {
        first_time = 1;

#ifdef _WIN32
        QueryPerformanceFrequency(&counter_freq);
        QueryPerformanceCounter(&start_count);
#else
        start_time((uint64_t *)&startTimesSeconds, (uint64_t *)&startTimesuSeconds);
#endif
    } else {
#ifdef _WIN32
        QueryPerformanceCounter(&now_count);
        elapsed_time =
            (double)(now_count.QuadPart - start_count.QuadPart) / (double)counter_freq.QuadPart;
#else
        finish_time((uint64_t *)&currentTimesSeconds, (uint64_t *)&currentTimesuSeconds);
        compute_overall_elapsed_time(startTimesSeconds,
                                     startTimesuSeconds,
                                     currentTimesSeconds,
                                     currentTimesuSeconds,
                                     &elapsed_time);
#endif

        predicted_time = (processed_frame_count - buffer_frames) * injector_interval;
        millisec_ahead = (int32_t)(1000 * (predicted_time - elapsed_time));
        if (millisec_ahead > 0) {
            //  timeBeginPeriod(1);
            sleep_ms(millisec_ahead);
            //  timeEndPeriod (1);
        }
    }
}
