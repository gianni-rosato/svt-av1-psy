/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbEntropyCodingUtil_h
#define EbEntropyCodingUtil_h

#include "EbBitstreamUnit.h"

/**************************************
* Defines
**************************************/
#ifdef __cplusplus
extern "C" {
#endif

/**************************************
    * Data Structures
    **************************************/
typedef struct BacEncContext {
    OutputBitstreamUnit* m_pc_t_com_bit_if;
} BacEncContext;

typedef struct CabacEncodeContext {
    BacEncContext bac_enc_context;
} CabacEncodeContext;

#ifdef __cplusplus
}
#endif
#endif //EbEntropyCodingUtil_h
