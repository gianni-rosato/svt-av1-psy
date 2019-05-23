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

#ifndef SVT_AV1_BITSTREAM_READER_MOCK_H_
#define SVT_AV1_BITSTREAM_READER_MOCK_H_
#include <stdint.h>

#define ACCT_STR_PARAM
#define ACCT_STR_ARG(s)

#define aom_read(r, prob, ACCT_STR_NAME) \
    aom_read_(r, prob ACCT_STR_ARG(ACCT_STR_NAME))
#define aom_read_bit(r, ACCT_STR_NAME) \
    aom_read_bit_(r ACCT_STR_ARG(ACCT_STR_NAME))
#define aom_read_tree(r, tree, probs, ACCT_STR_NAME) \
    aom_read_tree_(r, tree, probs ACCT_STR_ARG(ACCT_STR_NAME))
#define aom_read_literal(r, bits, ACCT_STR_NAME) \
    aom_read_literal_(r, bits ACCT_STR_ARG(ACCT_STR_NAME))
#define aom_read_cdf(r, cdf, nsymbs, ACCT_STR_NAME) \
    aom_read_cdf_(r, cdf, nsymbs ACCT_STR_ARG(ACCT_STR_NAME))
#define aom_read_symbol(r, cdf, nsymbs, ACCT_STR_NAME) \
    aom_read_symbol_(r, cdf, nsymbs ACCT_STR_ARG(ACCT_STR_NAME))

#ifdef __cplusplus
extern "C" {
#endif

#include "EbBitstreamUnit.h"
#include "EbDefinitions.h"

// entdec.h
typedef struct od_ec_dec od_ec_dec;

#if defined(OD_ACCOUNTING) && OD_ACCOUNTING
#define OD_ACC_STR , char *acc_str
#define od_ec_dec_bits(dec, ftb, str) od_ec_dec_bits_(dec, ftb, str)
#else
#define OD_ACC_STR
#define od_ec_dec_bits(dec, ftb, str) od_ec_dec_bits_(dec, ftb)
#endif

/*The entropy decoder context.*/
struct od_ec_dec {
    /*The start of the current input buffer.*/
    const unsigned char *buf;
    /*An offset used to keep track of tell after reaching the end of the stream.
      This is constant throughout most of the decoding process, but becomes
       important once we hit the end of the buffer and stop incrementing bptr
       (and instead pretend cnt has lots of bits).*/
    int32_t tell_offs;
    /*The end of the current input buffer.*/
    const unsigned char *end;
    /*The read pointer for the entropy-coded bits.*/
    const unsigned char *bptr;
    /*The difference between the high end of the current range, (low + rng), and
       the coded value, minus 1.
      This stores up to OD_EC_WINDOW_SIZE bits of that difference, but the
       decoder only uses the top 16 bits of the window to decode the next
      symbol. As we shift up during renormalization, if we don't have enough
      bits left in the window to fill the top 16, we'll read in more bits of the
      coded value.*/
    od_ec_window dif;
    /*The number of values in the current range.*/
    uint16_t rng;
    /*The number of bits of data in the current value.*/
    int16_t cnt;
};

/*See entdec.c for further documentation.*/

void od_ec_dec_init(od_ec_dec *dec, const unsigned char *buf, uint32_t storage)
    OD_ARG_NONNULL(1) OD_ARG_NONNULL(2);

OD_WARN_UNUSED_RESULT int od_ec_decode_bool_q15(od_ec_dec *dec, unsigned f)
    OD_ARG_NONNULL(1);
OD_WARN_UNUSED_RESULT int od_ec_decode_cdf_q15(od_ec_dec *dec,
                                               const uint16_t *cdf, int nsyms)
    OD_ARG_NONNULL(1) OD_ARG_NONNULL(2);

OD_WARN_UNUSED_RESULT uint32_t od_ec_dec_bits_(od_ec_dec *dec, unsigned ftb)
    OD_ARG_NONNULL(1);

OD_WARN_UNUSED_RESULT int od_ec_dec_tell(const od_ec_dec *dec)
    OD_ARG_NONNULL(1);
OD_WARN_UNUSED_RESULT uint32_t od_ec_dec_tell_frac(const od_ec_dec *dec)
    OD_ARG_NONNULL(1);

// daalaboolreader.h
struct daala_reader {
    const uint8_t *buffer;
    const uint8_t *buffer_end;
    od_ec_dec ec;
    uint8_t allow_update_cdf;
};

typedef struct daala_reader daala_reader;

int aom_daala_reader_init(daala_reader *r, const uint8_t *buffer, int size);
const uint8_t *aom_daala_reader_find_begin(daala_reader *r);
const uint8_t *aom_daala_reader_find_end(daala_reader *r);
uint32_t aom_daala_reader_tell(const daala_reader *r);
uint32_t aom_daala_reader_tell_frac(const daala_reader *r);
// Returns true if the reader has tried to decode more data from the buffer
// than was actually provided.
int aom_daala_reader_has_overflowed(const daala_reader *r);

static INLINE int aom_daala_read(daala_reader *r, int prob) {
    int bit;
    int p = (0x7FFFFF - (prob << 15) + prob) >> 8;

    bit = od_ec_decode_bool_q15(&r->ec, p);
    return bit;
}

static INLINE int daala_read_symbol(daala_reader *r, const AomCdfProb *cdf,
                                    int nsymbs) {
    int symb;
    assert(cdf != nullptr);
    symb = od_ec_decode_cdf_q15(&r->ec, cdf, nsymbs);

    return symb;
}

/* bitreader.h */
typedef struct daala_reader aom_reader;

static INLINE int aom_reader_init(aom_reader *r, const uint8_t *buffer,
                                  size_t size) {
    return aom_daala_reader_init(r, buffer, (int)size);
}

static INLINE const uint8_t *aom_reader_find_begin(aom_reader *r) {
    return aom_daala_reader_find_begin(r);
}

static INLINE const uint8_t *aom_reader_find_end(aom_reader *r) {
    return aom_daala_reader_find_end(r);
}

// Returns true if the bit reader has tried to decode more data from the buffer
// than was actually provided.
static INLINE int aom_reader_has_overflowed(const aom_reader *r) {
    return aom_daala_reader_has_overflowed(r);
}

// Returns the position in the bit reader in bits.
static INLINE uint32_t aom_reader_tell(const aom_reader *r) {
    return aom_daala_reader_tell(r);
}

// Returns the position in the bit reader in 1/8th bits.
static INLINE uint32_t aom_reader_tell_frac(const aom_reader *r) {
    return aom_daala_reader_tell_frac(r);
}

static INLINE int aom_read_(aom_reader *r, int prob) {
    int ret;
    ret = aom_daala_read(r, prob);
    return ret;
}

static INLINE int aom_read_bit_(aom_reader *r) {
    int ret;
    ret = aom_read(r, 128, nullptr);  // aom_prob_half
    return ret;
}

static INLINE int aom_read_literal_(aom_reader *r, int bits) {
    int literal = 0, bit;

    for (bit = bits - 1; bit >= 0; bit--)
        literal |= aom_read_bit(r, nullptr) << bit;
    return literal;
}

static INLINE int aom_read_cdf_(aom_reader *r, const AomCdfProb *cdf,
                                int nsymbs ACCT_STR_PARAM) {
    int ret;
    ret = daala_read_symbol(r, cdf, nsymbs);
    return ret;
}

static INLINE int aom_read_symbol_(aom_reader *r, AomCdfProb *cdf,
                                   int nsymbs ACCT_STR_PARAM) {
    int ret;
    ret = aom_read_cdf(r, cdf, nsymbs, ACCT_STR_NAME);
    if (r->allow_update_cdf)
        update_cdf(cdf, ret, nsymbs);
    return ret;
}

#ifdef __cplusplus
}  // extern "C"
#endif

#endif
