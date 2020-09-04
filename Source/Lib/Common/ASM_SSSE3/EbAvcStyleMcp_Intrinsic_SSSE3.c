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

#include "EbAvcStyleMcp_SSSE3.h"

#include "EbMcp_SSE2.h"
#include "EbDefinitions.h"
#include "EbAvcStyleMcp_SSE2.h"
#include "common_dsp_rtcd.h"
#include <emmintrin.h>
#include <tmmintrin.h>
