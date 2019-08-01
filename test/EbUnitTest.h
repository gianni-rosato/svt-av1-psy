/*
 * Copyright(c) 2019 Intel Corporation
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

#ifndef EbUnitTest_h
#define EbUnitTest_h

#include "EbDefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif

static const int EB_UNIT_TEST_NUM = 10;
static const int EB_UNIT_TEST_MAX_BD = 12;

extern void *eb_unit_test_buf_alligned;    // 256-byte aligned buffer
extern void *eb_unit_test_buf_unalligned;  // unaligned buffer
extern const char eb_unit_test_result_str[2][25];

#ifdef __cplusplus
}
#endif

#ifdef _WIN32
#define TEST_ALLIGN_MALLOC(type, pointer, n_elements) \
pointer = (type) _aligned_malloc(n_elements, ALVALUE); \

#define TEST_ALLIGN_FREE(pointer) _aligned_free(pointer);

#else
#define TEST_ALLIGN_MALLOC(type, pointer, n_elements) \
posix_memalign((void**)(&(pointer)), ALVALUE, n_elements); \

#define TEST_ALLIGN_FREE(pointer) free(pointer);

#endif

#endif // EbUnitTest_h
