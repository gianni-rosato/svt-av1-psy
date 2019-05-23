/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "EbAppConfig.h"

int32_t read_y4m_header(EbConfig *cfg);

int32_t read_y4m_frame_delimiter(EbConfig *cfg);

EbBool check_if_y4m(EbConfig *cfg);
