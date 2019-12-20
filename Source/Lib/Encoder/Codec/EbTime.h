// clang-format off
/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbTime_h
#define EbTime_h

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

#define NANOSECS_PER_SEC ((uint32_t)(1000000000L))

void EbStartTime(uint64_t *Startseconds, uint64_t *Startuseconds);
void EbFinishTime(uint64_t *Finishseconds, uint64_t *Finishuseconds);
void EbComputeOverallElapsedTime(uint64_t Startseconds, uint64_t Startuseconds, uint64_t Finishseconds, uint64_t Finishuseconds, double *duration);
void EbComputeOverallElapsedTimeMs(uint64_t Startseconds, uint64_t Startuseconds, uint64_t Finishseconds, uint64_t Finishuseconds, double *duration);
void EbInjector(uint64_t processedFrameCount, uint32_t injector_frame_rate);
void EbSleepMs(uint64_t milliSeconds);

#ifdef __cplusplus
}
#endif // __cplusplus

#endif // EbTime_h
/* File EOF */
// clang-format on
