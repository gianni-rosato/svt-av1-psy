/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
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
