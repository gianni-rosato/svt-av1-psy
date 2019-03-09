/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "EbAppConfig.h"

int32_t readY4mHeader(EbConfig_t *cfg);

int32_t readY4mFrameDelimiter(EbConfig_t *cfg);

EbBool checkIfY4m(EbConfig_t *cfg);
