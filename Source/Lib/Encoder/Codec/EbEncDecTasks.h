/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbEncDecTasks_h
#define EbEncDecTasks_h

#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"
#include "EbObject.h"
#ifdef __cplusplus
extern "C" {
#endif
#define ENCDEC_TASKS_MDC_INPUT 0
#define ENCDEC_TASKS_ENCDEC_INPUT 1
#define ENCDEC_TASKS_CONTINUE 2

/**************************************
     * Process Results
     **************************************/
typedef struct EncDecTasks {
    EbDctor          dctor;
    EbObjectWrapper *pcs_wrapper_ptr;
    uint32_t         input_type;
    int16_t          enc_dec_segment_row;
    uint16_t         tile_group_index;
} EncDecTasks;

typedef struct EncDecTasksInitData {
    unsigned enc_dec_segment_row_count;
} EncDecTasksInitData;

/**************************************
     * Extern Function Declarations
     **************************************/
extern EbErrorType enc_dec_tasks_creator(EbPtr *object_dbl_ptr, EbPtr object_init_data_ptr);

#ifdef __cplusplus
}
#endif
#endif // EbEncDecTasks_h
