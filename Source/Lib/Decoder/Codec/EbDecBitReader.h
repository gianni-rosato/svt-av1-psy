/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

/*
* Copyright (c) 2016, Alliance for Open Media. All rights reserved
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at www.aomedia.org/license/software. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at www.aomedia.org/license/patent.
*/

#ifndef EbDecBitReader_h
#define EbDecBitReader_h

#include "EbDecBitstreamUnit.h"

#ifdef __cplusplus
extern "C" {
#endif

/********************************************************************************************************************************/
/********************************************************************************************************************************/
/********************************************************************************************************************************/
/********************************************************************************************************************************/
// bitreader.h from AOM

#define ACCT_STR_PARAM
#define ACCT_STR_ARG(s)

#define svt_read(r, prob, ACCT_STR_NAME) \
  aom_read_(r, prob ACCT_STR_ARG(ACCT_STR_NAME))
#define svt_read_bit(r, ACCT_STR_NAME) \
  aom_read_bit_(r ACCT_STR_ARG(ACCT_STR_NAME))
#define svt_read_literal(r, bits, ACCT_STR_NAME) \
  aom_read_literal_(r, bits ACCT_STR_ARG(ACCT_STR_NAME))
#define svt_read_cdf(r, cdf, nsymbs, ACCT_STR_NAME) \
  aom_read_cdf_(r, cdf, nsymbs ACCT_STR_ARG(ACCT_STR_NAME))
#define svt_read_symbol(r, cdf, nsymbs, ACCT_STR_NAME) \
  aom_read_symbol_(r, cdf, nsymbs ACCT_STR_ARG(ACCT_STR_NAME))

typedef DaalaReader_t SvtReader;

int svt_reader_init(SvtReader   *r,
    const uint8_t   *buffer,
    size_t          size);

const uint8_t *aom_reader_find_begin(SvtReader *r);

const uint8_t *aom_reader_find_end(SvtReader *r);

int aom_read_(SvtReader *r, int prob ACCT_STR_PARAM);

int aom_read_bit_(SvtReader *r ACCT_STR_PARAM);

int aom_read_literal_(SvtReader *r, int bits ACCT_STR_PARAM);

int aom_read_cdf_(SvtReader *r,
    const AomCdfProb *cdf,
    int                 nsymbs ACCT_STR_PARAM);

int aom_read_symbol_(SvtReader   *r,
    AomCdfProb *cdf,
    int          nsymbs ACCT_STR_PARAM);

#ifdef __cplusplus
}
#endif
#endif // EbDecBitReader_h
