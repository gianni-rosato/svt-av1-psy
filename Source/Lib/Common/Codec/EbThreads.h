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

#ifndef EbThreads_h
#define EbThreads_h

#include "EbDefinitions.h"

#ifdef _WIN32
#include <windows.h>
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
extern EbHandle svt_create_thread(void *thread_function(void *), void *thread_context);

extern EbErrorType svt_start_thread(EbHandle thread_handle);

extern EbErrorType svt_stop_thread(EbHandle thread_handle);

extern EbErrorType svt_destroy_thread(EbHandle thread_handle);

/**************************************
     * Semaphores
     **************************************/
extern EbHandle svt_create_semaphore(uint32_t initial_count, uint32_t max_count);

extern EbErrorType svt_post_semaphore(EbHandle semaphore_handle);

extern EbErrorType svt_block_on_semaphore(EbHandle semaphore_handle);

extern EbErrorType svt_destroy_semaphore(EbHandle semaphore_handle);

/**************************************
     * Mutex
     **************************************/
extern EbHandle svt_create_mutex(void);
extern EbErrorType svt_release_mutex(EbHandle mutex_handle);
extern EbErrorType svt_block_on_mutex(EbHandle mutex_handle);
extern EbErrorType svt_destroy_mutex(EbHandle mutex_handle);
extern EbMemoryMapEntry *memory_map; // library Memory table
extern uint32_t *        memory_map_index; // library memory index
extern uint64_t *        total_lib_memory; // library Memory malloc'd
#ifdef _WIN32

#define EB_CREATE_THREAD(pointer, thread_function, thread_context)   \
    do {                                                             \
        pointer = svt_create_thread(thread_function, thread_context);\
        EB_ADD_MEM(pointer, 1, EB_THREAD);                           \
        if (num_groups == 1)                                         \
            SetThreadAffinityMask(pointer, group_affinity.Mask);     \
        else if (num_groups == 2 && alternate_groups) {              \
            group_affinity.Group = 1 - group_affinity.Group;         \
            SetThreadGroupAffinity(pointer, &group_affinity, NULL);  \
        } else if (num_groups == 2 && !alternate_groups)             \
            SetThreadGroupAffinity(pointer, &group_affinity, NULL);  \
    } while (0)

#elif defined(__linux__)
#ifndef __USE_GNU
#define __USE_GNU
#endif
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <sched.h>
#include <pthread.h>
#define EB_CREATE_THREAD(pointer, thread_function, thread_context)                           \
    do {                                                                                     \
        pointer = svt_create_thread(thread_function, thread_context);                        \
        EB_ADD_MEM(pointer, 1, EB_THREAD);                                                   \
        pthread_setaffinity_np(*((pthread_t *)pointer), sizeof(cpu_set_t), &group_affinity); \
    } while (0)
#else
#define EB_CREATE_THREAD(pointer, thread_function, thread_context)   \
    do {                                                             \
        pointer = svt_create_thread(thread_function, thread_context);\
        EB_ADD_MEM(pointer, 1, EB_THREAD);                           \
    } while (0)
#endif

#define EB_DESTROY_THREAD(pointer)                   \
    do {                                             \
        if (pointer) {                               \
            svt_destroy_thread(pointer);             \
            EB_REMOVE_MEM_ENTRY(pointer, EB_THREAD); \
            pointer = NULL;                          \
        }                                            \
    } while (0);

#define EB_CREATE_THREAD_ARRAY(pa, count, thread_function, thread_contexts) \
    do {                                                                    \
        EB_ALLOC_PTR_ARRAY(pa, count);                                      \
        for (uint32_t i = 0; i < count; i++)                                \
            EB_CREATE_THREAD(pa[i], thread_function, thread_contexts[i]);   \
    } while (0)

#define EB_DESTROY_THREAD_ARRAY(pa, count)                                 \
    do {                                                                   \
        if (pa) {                                                          \
            for (uint32_t i = 0; i < count; i++) EB_DESTROY_THREAD(pa[i]); \
            EB_FREE_PTR_ARRAY(pa, count);                                  \
        }                                                                  \
    } while (0)

#ifdef __cplusplus
}
#endif
#endif // EbThreads_h
