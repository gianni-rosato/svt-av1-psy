/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbSvtAv1Time_h
#define EbSvtAv1Time_h

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

#include "EbSvtAv1.h"

#define NANOSECS_PER_SEC ((uint32_t)(1000000000L))

EB_API void EbStartTime(uint64_t *Startseconds, uint64_t *Startuseconds);
EB_API void EbFinishTime(uint64_t *Finishseconds, uint64_t *Finishuseconds);
EB_API void EbComputeOverallElapsedTime(uint64_t Startseconds, uint64_t Startuseconds, uint64_t Finishseconds, uint64_t Finishuseconds, double *duration);
EB_API void EbComputeOverallElapsedTimeMs(uint64_t Startseconds, uint64_t Startuseconds, uint64_t Finishseconds, uint64_t Finishuseconds, double *duration);
EB_API void EbSleep(uint64_t milliSeconds);
EB_API void EbInjector(uint64_t processedFrameCount, uint32_t injector_frame_rate);

#ifdef __cplusplus
}
#endif // __cplusplus

#endif // EbSvtAv1Time_h
/* File EOF */
