/*
 * Copyright (c) 2018, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#ifndef AOM_AV1_ENCODER_DWT_H_
#define AOM_AV1_ENCODER_DWT_H_

//#include "av1/common/common.h"
//#include "av1/common/enums.h"
#include "EbDefinitions.h"
#define DWT_MAX_LENGTH 64
// Note:
// tran_low_t  is the datatype used for final transform coefficients.
// tran_high_t is the datatype used for intermediate transform stages.
typedef int64_t tran_high_t;
typedef int32_t tran_low_t;
void av1_fdwt8x8(tran_low_t *input, tran_low_t *output, int stride);
void svt_av1_fdwt8x8_uint8_input_c(uint8_t *input, tran_low_t *output, int stride,
                                   int hbd);
#endif  // AOM_AV1_ENCODER_DWT_H_
