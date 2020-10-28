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
#ifndef EbMalloc_h
#define EbMalloc_h
#include <stdlib.h>
#include <stdio.h>

#include "EbSvtAv1Enc.h"
#include "EbDefinitions.h"

#ifndef NDEBUG
#define DEBUG_MEMORY_USAGE
#endif

void svt_print_alloc_fail(const char* file, int line);

#ifdef DEBUG_MEMORY_USAGE
void svt_print_memory_usage(void);
void svt_increase_component_count(void);
void svt_decrease_component_count(void);
void svt_add_mem_entry(void* ptr, EbPtrType type, size_t count, const char* file, uint32_t line);
void svt_remove_mem_entry(void* ptr, EbPtrType type);

#define EB_ADD_MEM_ENTRY(p, type, count) svt_add_mem_entry(p, type, count, __FILE__, __LINE__)

#define EB_REMOVE_MEM_ENTRY(p, type) svt_remove_mem_entry(p, type);

#else
#define svt_print_memory_usage() \
    do {                         \
    } while (0)
#define svt_increase_component_count() \
    do {                               \
    } while (0)
#define svt_decrease_component_count() \
    do {                               \
    } while (0)
#define EB_ADD_MEM_ENTRY(p, type, count) \
    do {                                 \
    } while (0)
#define EB_REMOVE_MEM_ENTRY(p, type) \
    do {                             \
    } while (0)

#endif //DEBUG_MEMORY_USAGE

#define EB_NO_THROW_ADD_MEM(p, size, type)                                               \
    do {                                                                                 \
        if (!p)                                                                          \
            svt_print_alloc_fail(__FILE__, __LINE__);                                    \
        else                                                                             \
            EB_ADD_MEM_ENTRY(p, type, size);                                             \
    } while (0)

#define EB_CHECK_MEM(p)                           \
    do {                                          \
        if (!p)                                   \
            return EB_ErrorInsufficientResources; \
    } while (0)

#define EB_ADD_MEM(p, size, type)           \
    do {                                    \
        EB_NO_THROW_ADD_MEM(p, size, type); \
        EB_CHECK_MEM(p);                    \
    } while (0)

#define EB_NO_THROW_MALLOC(pointer, size)                \
    do {                                                 \
        void* malloced_p = malloc(size);                 \
        EB_NO_THROW_ADD_MEM(malloced_p, size, EB_N_PTR); \
        pointer = malloced_p;                            \
    } while (0)

#define EB_MALLOC(pointer, size)           \
    do {                                   \
        EB_NO_THROW_MALLOC(pointer, size); \
        EB_CHECK_MEM(pointer);             \
    } while (0)

#define EB_NO_THROW_CALLOC(pointer, count, size)             \
    do {                                                     \
        pointer = calloc(count, size);                       \
        EB_NO_THROW_ADD_MEM(pointer, count* size, EB_C_PTR); \
    } while (0)

#define EB_CALLOC(pointer, count, size)           \
    do {                                          \
        EB_NO_THROW_CALLOC(pointer, count, size); \
        EB_CHECK_MEM(pointer);                    \
    } while (0)

#define EB_FREE(pointer)                        \
    do {                                        \
        EB_REMOVE_MEM_ENTRY(pointer, EB_N_PTR); \
        free(pointer);                          \
        pointer = NULL;                         \
    } while (0)

#define EB_MALLOC_ARRAY(pa, count) \
    do { EB_MALLOC(pa, sizeof(*(pa)) * (count)); } while (0)

#define EB_REALLOC_ARRAY(pa, count)             \
    do {                                        \
        size_t size = sizeof(*(pa)) * (count);  \
        void* p = realloc(pa, size);            \
        if (p) {                                \
            EB_REMOVE_MEM_ENTRY(pa, EB_N_PTR);  \
        }                                       \
        EB_ADD_MEM(p, size, EB_N_PTR);          \
        pa = p;                                 \
    } while (0)

#define EB_CALLOC_ARRAY(pa, count) \
    do { EB_CALLOC(pa, count, sizeof(*(pa))); } while (0)

#define EB_FREE_ARRAY(pa) EB_FREE(pa);

#define EB_ALLOC_PTR_ARRAY(pa, count) \
    do { EB_CALLOC(pa, count, sizeof(*(pa))); } while (0)

#define EB_FREE_PTR_ARRAY(pa, count)                           \
    do {                                                       \
        if (pa) {                                              \
            for (size_t i = 0; i < count; i++) EB_FREE(pa[i]); \
            EB_FREE(pa);                                       \
        }                                                      \
    } while (0)

#define EB_MALLOC_2D(p2d, width, height)                                     \
    do {                                                                     \
        EB_MALLOC_ARRAY(p2d, width);                                         \
        EB_MALLOC_ARRAY(p2d[0], (width) * (height));                         \
        for (size_t w = 1; w < (width); w++) p2d[w] = p2d[0] + w * (height); \
    } while (0)

#define EB_CALLOC_2D(p2d, width, height)                                     \
    do {                                                                     \
        EB_MALLOC_ARRAY(p2d, width);                                         \
        EB_CALLOC_ARRAY(p2d[0], (width) * (height));                         \
        for (size_t w = 1; w < (width); w++) p2d[w] = p2d[0] + w * (height); \
    } while (0)

#define EB_FREE_2D(p2d)            \
    do {                           \
        if (p2d)                   \
            EB_FREE_ARRAY(p2d[0]); \
        EB_FREE_ARRAY(p2d);        \
    } while (0)

#ifdef _WIN32
#define EB_MALLOC_ALIGNED(pointer, size)          \
    do {                                          \
        pointer = _aligned_malloc(size, ALVALUE); \
        EB_ADD_MEM(pointer, size, EB_A_PTR);      \
    } while (0)

#define EB_FREE_ALIGNED(pointer)                \
    do {                                        \
        EB_REMOVE_MEM_ENTRY(pointer, EB_A_PTR); \
        _aligned_free(pointer);                 \
        pointer = NULL;                         \
    } while (0)
#else
#define EB_MALLOC_ALIGNED(pointer, size)                            \
    do {                                                            \
        if (posix_memalign((void**)&(pointer), ALVALUE, size) != 0) \
            return EB_ErrorInsufficientResources;                   \
        EB_ADD_MEM(pointer, size, EB_A_PTR);                        \
    } while (0)

#define EB_FREE_ALIGNED(pointer)                \
    do {                                        \
        EB_REMOVE_MEM_ENTRY(pointer, EB_A_PTR); \
        free(pointer);                          \
        pointer = NULL;                         \
    } while (0)
#endif

#define EB_MALLOC_ALIGNED_ARRAY(pa, count) EB_MALLOC_ALIGNED(pa, sizeof(*(pa)) * (count))

#define EB_CALLOC_ALIGNED_ARRAY(pa, count)              \
    do {                                                \
        EB_MALLOC_ALIGNED(pa, sizeof(*(pa)) * (count)); \
        memset(pa, 0, sizeof(*(pa)) * (count));         \
    } while (0)

#define EB_FREE_ALIGNED_ARRAY(pa) EB_FREE_ALIGNED(pa)

#endif //EbMalloc_h
