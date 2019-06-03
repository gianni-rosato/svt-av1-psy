/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdint.h>

#ifdef _WIN32
//#if  (WIN_ENCODER_TIMING || WIN_DECODER_TIMING)
#include <time.h>
#include <stdio.h>
#include <windows.h>
//#endif

#elif defined(__linux__) || defined(__APPLE__)
#include <stdio.h>
#include <stdlib.h>
//#if   (LINUX_ENCODER_TIMING || LINUX_DECODER_TIMING)
#include <sys/time.h>
#include <time.h>
//#endif

#else
#error OS/Platform not supported.
#endif

#include "EbTime.h"

#if defined(__linux__)
#ifndef __clang__
__attribute__((optimize("unroll-loops")))
#endif
#endif

void EbStartTime(uint64_t *Startseconds, uint64_t *Startuseconds) {
#if defined(__linux__) || defined(__APPLE__) //(LINUX_ENCODER_TIMING || LINUX_DECODER_TIMING)
    struct timeval start;
    gettimeofday(&start, NULL);
    *Startseconds = start.tv_sec;
    *Startuseconds = start.tv_usec;
#elif _WIN32 //(WIN_ENCODER_TIMING || WIN_DECODER_TIMING)
    *Startseconds = (uint64_t)clock();
    (void)(*Startuseconds);
#else
    (void)(*Startuseconds);
    (void)(*Startseconds);
#endif
}

void EbFinishTime(uint64_t *Finishseconds, uint64_t *Finishuseconds) {
#if defined(__linux__) || defined(__APPLE__) //(LINUX_ENCODER_TIMING || LINUX_DECODER_TIMING)
    struct timeval finish;
    gettimeofday(&finish, NULL);
    *Finishseconds = finish.tv_sec;
    *Finishuseconds = finish.tv_usec;
#elif _WIN32 //(WIN_ENCODER_TIMING || WIN_DECODER_TIMING)
    *Finishseconds = (uint64_t)clock();
    (void)(*Finishuseconds);
#else
    (void)(*Finishuseconds);
    (void)(*Finishseconds);
#endif
}

void EbComputeOverallElapsedTime(uint64_t Startseconds, uint64_t Startuseconds, uint64_t Finishseconds, uint64_t Finishuseconds, double *duration)
{
#if defined(__linux__) || defined(__APPLE__) //(LINUX_ENCODER_TIMING || LINUX_DECODER_TIMING)
    long   mtime, seconds, useconds;
    seconds = Finishseconds - Startseconds;
    useconds = Finishuseconds - Startuseconds;
    mtime = ((seconds) * 1000 + useconds / 1000.0) + 0.5;
    *duration = (double)mtime / 1000;
    //printf("\nElapsed time: %3.3ld seconds\n", mtime/1000);
#elif _WIN32 //(WIN_ENCODER_TIMING || WIN_DECODER_TIMING)
    //double  duration;
    *duration = (double)(Finishseconds - Startseconds) / CLOCKS_PER_SEC;
    //printf("\nElapsed time: %3.3f seconds\n", *duration);
    (void)(Startuseconds);
    (void)(Finishuseconds);
#else
    (void)(Startuseconds);
    (void)(Startseconds);
    (void)(Finishuseconds);
    (void)(Finishseconds);

#endif
}

void EbComputeOverallElapsedTimeMs(uint64_t Startseconds, uint64_t Startuseconds, uint64_t Finishseconds, uint64_t Finishuseconds, double *duration)
{
#if defined(__linux__) || defined(__APPLE__) //(LINUX_ENCODER_TIMING || LINUX_DECODER_TIMING)
    long   mtime, seconds, useconds;
    seconds = Finishseconds - Startseconds;
    useconds = Finishuseconds - Startuseconds;
    mtime = ((seconds) * 1000 + useconds / 1000.0) + 0.5;
    *duration = (double)mtime;
    //printf("\nElapsed time: %3.3ld seconds\n", mtime/1000);
#elif _WIN32 //(WIN_ENCODER_TIMING || WIN_DECODER_TIMING)
    //double  duration;
    *duration = (double)(Finishseconds - Startseconds);
    //printf("\nElapsed time: %3.3f seconds\n", *duration);
    (void)(Startuseconds);
    (void)(Finishuseconds);
#else
    (void)(Startuseconds);
    (void)(Startseconds);
    (void)(Finishuseconds);
    (void)(Finishseconds);

#endif
}

static void EbSleepMs(uint64_t milliSeconds)
{
    if(milliSeconds) {
#if defined(__linux__) || defined(__APPLE__)
        struct timespec req,rem;
        req.tv_sec=(int32_t)(milliSeconds/1000);
        milliSeconds -= req.tv_sec * 1000;
        req.tv_nsec = milliSeconds * 1000000UL;
        nanosleep(&req,&rem);
#elif _WIN32
        Sleep((DWORD) milliSeconds);
#else
#error OS Not supported
#endif
    }
}

void EbInjector(uint64_t processed_frame_count,
              uint32_t injector_frame_rate)
{
#if defined(__linux__) || defined(__APPLE__)
    uint64_t                  currentTimesSeconds = 0;
    uint64_t                  currentTimesuSeconds = 0;
    static uint64_t           startTimesSeconds;
    static uint64_t           startTimesuSeconds;
#elif _WIN32
    static LARGE_INTEGER    startCount;               // this is the start time
    static LARGE_INTEGER    counterFreq;              // performance counter frequency
    LARGE_INTEGER           nowCount;                 // this is the current time
#else
#error OS Not supported
#endif

    double                 injectorInterval  = (double)(1<<16)/(double)injector_frame_rate;     // 1.0 / injector frame rate (in this case, 1.0/encodRate)
    double                  elapsedTime;
    double                  predictedTime;
    int32_t                     bufferFrames = 1;         // How far ahead of time should we let it get
    int32_t                     milliSecAhead;
    static int32_t              firstTime = 0;

    if (firstTime == 0)
    {
        firstTime = 1;

#if defined(__linux__) || defined(__APPLE__)
        EbStartTime((uint64_t*)&startTimesSeconds, (uint64_t*)&startTimesuSeconds);
#elif _WIN32
        QueryPerformanceFrequency(&counterFreq);
        QueryPerformanceCounter(&startCount);
#endif
    }
    else
    {
#if defined(__linux__) || defined(__APPLE__)
        EbFinishTime((uint64_t*)&currentTimesSeconds, (uint64_t*)&currentTimesuSeconds);
        EbComputeOverallElapsedTime(startTimesSeconds, startTimesuSeconds, currentTimesSeconds, currentTimesuSeconds, &elapsedTime);
#elif _WIN32
        QueryPerformanceCounter(&nowCount);
        elapsedTime = (double)(nowCount.QuadPart - startCount.QuadPart) / (double)counterFreq.QuadPart;
#endif

        predictedTime = (processed_frame_count - bufferFrames) * injectorInterval;
        milliSecAhead = (int32_t)(1000 * (predictedTime - elapsedTime ));
        if (milliSecAhead>0)
        {
            //  timeBeginPeriod(1);
            EbSleepMs(milliSecAhead);
            //  timeEndPeriod (1);
        }
    }
}
