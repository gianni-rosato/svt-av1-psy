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
