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

#ifndef EbEncDecResults_h
#define EbEncDecResults_h

#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"
#include "EbObject.h"
#ifdef __cplusplus
extern "C" {
#endif
/**************************************
     * Process Results
     **************************************/
typedef struct EncDecResults {
    EbDctor          dctor;
    EbObjectWrapper *pcs_wrapper_ptr;
    uint32_t         completed_sb_row_index_start;
    uint32_t         completed_sb_row_count;
} EncDecResults;

typedef struct DlfResults {
    EbDctor          dctor;
    EbObjectWrapper *pcs_wrapper_ptr;
    uint32_t         segment_index;
} DlfResults;

typedef struct CdefResults {
    EbDctor          dctor;
    EbObjectWrapper *pcs_wrapper_ptr;
    uint32_t         segment_index;
} CdefResults;

typedef struct RestResults {
    EbDctor          dctor;
    EbObjectWrapper *pcs_wrapper_ptr;
    uint32_t         completed_sb_row_index_start;
    uint32_t         completed_sb_row_count;
    uint16_t         tile_index;
} RestResults;

typedef struct EncDecResultsInitData {
    uint32_t junk;
} EncDecResultsInitData;

/**************************************
     * Extern Function Declarations
     **************************************/
extern EbErrorType enc_dec_results_creator(EbPtr *object_dbl_ptr, EbPtr object_init_data_ptr);

#ifdef __cplusplus
}
#endif
#endif // EbEncDecResults_h
