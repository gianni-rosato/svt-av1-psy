/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbThreads_h
#define EbThreads_h

#include "EbDefinitions.h"

#ifdef _MSC_VER
#include <Windows.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif
    // Create wrapper functions that hide thread calls,
    // semaphores, mutex, etc. These wrappers also hide
    // platform specific implementations of these objects.

    /**************************************
     * Threads
     **************************************/
    extern EbHandle EbCreateThread(
        void *threadFunction(void *),
        void *threadContext);
    extern EbErrorType EbStartThread(
        EbHandle threadHandle);
    extern EbErrorType EbStopThread(
        EbHandle threadHandle);
    extern EbErrorType EbDestroyThread(
        EbHandle threadHandle);

    /**************************************
     * Semaphores
     **************************************/
    extern EbHandle EbCreateSemaphore(
        uint32_t initialCount,
        uint32_t maxCount);
    extern EbErrorType EbPostSemaphore(
        EbHandle semaphoreHandle);
    extern EbErrorType EbBlockOnSemaphore(
        EbHandle semaphoreHandle);
    extern EbErrorType EbDestroySemaphore(
        EbHandle semaphoreHandle);
    /**************************************
     * Mutex
     **************************************/
    extern EbHandle EbCreateMutex(
        void);
    extern EbErrorType EbReleaseMutex(
        EbHandle mutexHandle);
    extern EbErrorType EbBlockOnMutex(
        EbHandle mutexHandle);
    extern EbErrorType EbBlockOnMutexTimeout(
        EbHandle mutexHandle,
        uint32_t timeout);
    extern EbErrorType EbDestroyMutex(
        EbHandle mutexHandle);

    extern    EbMemoryMapEntry        *memoryMap;               // library Memory table
    extern    uint32_t                  *memoryMapIndex;          // library memory index
    extern    uint64_t                  *totalLibMemory;          // library Memory malloc'd

#ifdef _MSC_VER
    extern    GROUP_AFFINITY           groupAffinity;
    extern    uint8_t                    numGroups;
    extern    EbBool                  alternateGroups;
#define EB_CREATETHREAD(type, pointer, nElements, pointerClass, threadFunction, threadContext) \
    pointer = EbCreateThread(threadFunction, threadContext); \
    if (pointer == (type)EB_NULL) { \
        return EB_ErrorInsufficientResources; \
    } \
    else { \
        memoryMap[*(memoryMapIndex)].ptrType = pointerClass; \
        memoryMap[(*(memoryMapIndex))++].ptr = pointer; \
        if (nElements % 8 == 0) { \
            *totalLibMemory += (nElements); \
        } \
        else { \
            *totalLibMemory += ((nElements) + (8 - ((nElements) % 8))); \
        } \
        if (numGroups == 2 && alternateGroups){ \
            groupAffinity.Group = 1 - groupAffinity.Group; \
            SetThreadGroupAffinity(pointer,&groupAffinity,NULL); \
        } \
        else if (numGroups == 2 && !alternateGroups){ \
            SetThreadGroupAffinity(pointer,&groupAffinity,NULL); \
        } \
    } \
    if (*(memoryMapIndex) >= MAX_NUM_PTR) { \
        return EB_ErrorInsufficientResources; \
    } \
    libThreadCount++;
#else
#define EB_CREATETHREAD(type, pointer, nElements, pointerClass, threadFunction, threadContext) \
    pointer = EbCreateThread(threadFunction, threadContext); \
    if (pointer == (type)EB_NULL) { \
        return EB_ErrorInsufficientResources; \
    } \
    else { \
        memoryMap[*(memoryMapIndex)].ptrType = pointerClass; \
        memoryMap[(*(memoryMapIndex))++].ptr = pointer; \
        if (nElements % 8 == 0) { \
            *totalLibMemory += (nElements); \
        } \
        else { \
            *totalLibMemory += ((nElements) + (8 - ((nElements) % 8))); \
        } \
    } \
    if (*(memoryMapIndex) >= MAX_NUM_PTR) { \
        return EB_ErrorInsufficientResources; \
    } \
    libThreadCount++;
#endif


#ifdef __cplusplus
}
#endif
#endif // EbThreads_h
