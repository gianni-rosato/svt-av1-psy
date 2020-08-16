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
