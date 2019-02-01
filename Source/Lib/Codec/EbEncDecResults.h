/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbEncDecResults_h
#define EbEncDecResults_h

#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"
#ifdef __cplusplus
extern "C" {
#endif
    /**************************************
     * Process Results
     **************************************/
    typedef struct EncDecResults_s
    {
        EbObjectWrapper_t      *pictureControlSetWrapperPtr;
        uint32_t                  completedLcuRowIndexStart;
        uint32_t                  completedLcuRowCount;

    } EncDecResults_t;

    typedef struct EncDecResultsInitData_s
    {
        uint32_t         junk;
    } EncDecResultsInitData_t;

    /**************************************
     * Extern Function Declarations
     **************************************/
    extern EbErrorType EncDecResultsCtor(
        EbPtr *object_dbl_ptr,
        EbPtr object_init_data_ptr);


#ifdef __cplusplus
}
#endif
#endif // EbEncDecResults_h