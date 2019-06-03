/*
* Copyright(c) 2019 Netflix, Inc.
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

// SUMMARY
//   Contains the Entropy coding functions

/**************************************
 * Includes
 **************************************/

#include "EbDefinitions.h"
#include "EbDecBitstreamUnit.h"
#include "EbDecBitReader.h"

/********************************************************************************************************************************/
/********************************************************************************************************************************/
/********************************************************************************************************************************/
/********************************************************************************************************************************/
// bitreader.h from AOM

int svt_reader_init(SvtReader   *r,
                                const uint8_t   *buffer,
                                size_t          size)
{
  return aom_daala_reader_init(r, buffer, (int)size);
}

const uint8_t *aom_reader_find_begin(SvtReader *r) {
  return aom_daala_reader_find_begin(r);
}

const uint8_t *aom_reader_find_end(SvtReader *r) {
  return aom_daala_reader_find_end(r);
}

int aom_read_(SvtReader *r, int prob ACCT_STR_PARAM) {
  int ret;
  ret = aom_daala_read(r, prob);
  return ret;
}

int aom_read_bit_(SvtReader *r ACCT_STR_PARAM) {
  int ret;
  ret = svt_read(r, 128, NULL);  // aom_prob_half
  return ret;
}

int aom_read_literal_(SvtReader *r, int bits ACCT_STR_PARAM) {
  int literal = 0, bit;

  for (bit = bits - 1; bit >= 0; bit--) literal |= svt_read_bit(r, NULL) << bit;
  return literal;
}

int aom_read_cdf_(SvtReader *r,
                        const AomCdfProb *cdf,
                        int                 nsymbs ACCT_STR_PARAM)
{
  int ret;
  ret = daala_read_symbol(r, cdf, nsymbs);

  return ret;
}

int aom_read_symbol_(SvtReader   *r,
                                   AomCdfProb *cdf,
                                   int          nsymbs ACCT_STR_PARAM)
{
  int ret;
  ret = svt_read_cdf(r, cdf, nsymbs, ACCT_STR_NAME);
  if (r->allow_update_cdf) dec_update_cdf(cdf, ret, nsymbs);
  return ret;
}
