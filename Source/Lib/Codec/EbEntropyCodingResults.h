/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbEntropyCodingResults_h
#define EbEntropyCodingResults_h

#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"
#ifdef __cplusplus
extern "C" {
#endif
    /**************************************
     * Process Results
     **************************************/
    typedef struct
    {
        EbObjectWrapper_t      *pictureControlSetWrapperPtr;

    } EntropyCodingResults_t;

    typedef struct
    {
        uint32_t         junk;
    } EntropyCodingResultsInitData_t;

    /**************************************
     * Extern Function Declarations
     **************************************/
    extern EbErrorType EntropyCodingResultsCtor(
        EbPtr *object_dbl_ptr,
        EbPtr object_init_data_ptr);


#ifdef __cplusplus
}
#endif
#endif // EbEntropyCodingResults_h