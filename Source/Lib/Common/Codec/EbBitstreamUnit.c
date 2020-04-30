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

#include <stdlib.h>

#include "EbBitstreamUnit.h"
#include "EbDefinitions.h"
#include "EbUtility.h"

#if OD_MEASURE_EC_OVERHEAD
#include <stdio.h>
#endif

static void output_bitstream_unit_dctor(EbPtr p) {
    OutputBitstreamUnit *obj = (OutputBitstreamUnit *)p;
    EB_FREE_ARRAY(obj->buffer_begin_av1);
}

/**********************************
 * Constructor
 **********************************/
EbErrorType output_bitstream_unit_ctor(OutputBitstreamUnit *bitstream_ptr, uint32_t buffer_size) {
    bitstream_ptr->dctor = output_bitstream_unit_dctor;
    if (buffer_size) {
        bitstream_ptr->size = buffer_size / sizeof(uint32_t);
        EB_MALLOC_ARRAY(bitstream_ptr->buffer_begin_av1, bitstream_ptr->size);
        bitstream_ptr->buffer_av1 = bitstream_ptr->buffer_begin_av1;
    } else {
        bitstream_ptr->size             = 0;
        bitstream_ptr->buffer_begin_av1 = 0;
        bitstream_ptr->buffer_av1       = 0;
    }

    return EB_ErrorNone;
}

/**********************************
 * Reset Bitstream
 **********************************/
EbErrorType output_bitstream_reset(OutputBitstreamUnit *bitstream_ptr) {
    EbErrorType return_error = EB_ErrorNone;

    // Reset the write ptr to the beginning of the buffer
    bitstream_ptr->buffer_av1 = bitstream_ptr->buffer_begin_av1;

    return return_error;
}

/********************************************************************************************************************************/
/********************************************************************************************************************************/
/********************************************************************************************************************************/
/********************************************************************************************************************************/
// daalaboolwriter.c
void eb_aom_daala_start_encode(DaalaWriter *br, uint8_t *source) {
    br->buffer = source;
    br->pos    = 0;
    eb_od_ec_enc_init(&br->ec, 62025);
}

int32_t eb_aom_daala_stop_encode(DaalaWriter *br) {
    int32_t  nb_bits;
    uint32_t daala_bytes;
    uint8_t *daala_data;
    daala_data = eb_od_ec_enc_done(&br->ec, &daala_bytes);
    nb_bits    = eb_od_ec_enc_tell(&br->ec);
    memcpy(br->buffer, daala_data, daala_bytes);
    br->pos = daala_bytes;
    eb_od_ec_enc_clear(&br->ec);
    return nb_bits;
}

/*A range encoder.
See entdec.c and the references for implementation details \cite{Mar79,MNW98}.

@INPROCEEDINGS{Mar79,
author="Martin, G.N.N.",
title="Range encoding: an algorithm for removing redundancy from a digitised
message",
booktitle="Video \& Data Recording Conference",
year=1979,
address="Southampton",
month=Jul,
URL="http://www.compressconsult.com/rangecoder/rngcod.pdf.gz"
}
@ARTICLE{MNW98,
author="Alistair Moffat and Radford Neal and Ian H. Witten",
title="Arithmetic Coding Revisited",
journal="{ACM} Transactions on Information Systems",
year=1998,
volume=16,
number=3,
pages="256--294",
month=Jul,
URL="http://researchcommons.waikato.ac.nz/Bitstream/handle/10289/78/content.pdf"
}*/

/*Takes updated low and range values, renormalizes them so that
32768 <= rng < 65536 (flushing bytes from low to the pre-carry buffer if
necessary), and stores them back in the encoder context.
low: The new value of low.
rng: The new value of the range.*/
static void od_ec_enc_normalize(OdEcEnc *enc, OdEcWindow low, unsigned rng) {
    int32_t d;
    int32_t c;
    int32_t s;
    c = enc->cnt;
    assert(rng <= 65535U);
    d = 16 - OD_ILOG_NZ(rng);
    s = c + d;
    /*TODO: Right now we flush every time we have at least one byte available.
    Instead we should use an OdEcWindow and flush right before we're about to
    shift bits off the end of the window.
    For a 32-bit window this is about the same amount of work, but for a 64-bit
    window it should be a fair win.*/
    if (s >= 0) {
        uint16_t *buf;
        uint32_t  storage;
        uint32_t  offs;
        unsigned  m;
        buf     = enc->precarry_buf;
        storage = enc->precarry_storage;
        offs    = enc->offs;
        if (offs + 2 > storage) {
            storage = 2 * storage + 2;
            buf = realloc(enc->precarry_buf, sizeof(*buf) * storage);
            if (!buf) {
                enc->error = -1;
                enc->offs  = 0;
                return;
            }
            enc->precarry_buf     = buf;
            enc->precarry_storage = storage;
        }
        c += 16;
        m = (1 << c) - 1;
        if (s >= 8) {
            assert(offs < storage);
            buf[offs++] = (uint16_t)(low >> c);
            low &= m;
            c -= 8;
            m >>= 8;
        }
        assert(offs < storage);
        buf[offs++] = (uint16_t)(low >> c);
        s           = c + d - 24;
        low &= m;
        enc->offs = offs;
    }
    enc->low = low << d;
    enc->rng = (int16_t)(rng << d);
    enc->cnt = (int16_t)s;
}

/*Initializes the encoder.
size: The initial size of the buffer, in bytes.*/
void eb_od_ec_enc_init(OdEcEnc *enc, uint32_t size) {
    eb_od_ec_enc_reset(enc);
    enc->buf     = (uint8_t *)malloc(sizeof(*enc->buf) * size);
    enc->storage = size;
    if (size > 0 && enc->buf == NULL) {
        enc->storage = 0;
        enc->error   = -1;
    }
    enc->precarry_buf     = (uint16_t *)malloc(sizeof(*enc->precarry_buf) * size);
    enc->precarry_storage = size;
    if (size > 0 && enc->precarry_buf == NULL) {
        enc->precarry_storage = 0;
        enc->error            = -1;
    }
}

/*Reinitializes the encoder.*/
void eb_od_ec_enc_reset(OdEcEnc *enc) {
    enc->offs = 0;
    enc->low  = 0;
    enc->rng  = 0x8000;
    /*This is initialized to -9 so that it crosses zero after we've accumulated
    one byte + one carry bit.*/
    enc->cnt   = -9;
    enc->error = 0;
#if OD_MEASURE_EC_OVERHEAD
    enc->entropy    = 0;
    enc->nb_symbols = 0;
#endif
}

/*Frees the buffers used by the encoder.*/
void eb_od_ec_enc_clear(OdEcEnc *enc) {
    free(enc->precarry_buf);
    free(enc->buf);
}

/*Encodes a symbol given its frequency in Q15.
fl: CDF_PROB_TOP minus the cumulative frequency of all symbols that come
before the
one to be encoded.
fh: CDF_PROB_TOP minus the cumulative frequency of all symbols up to and
including
the one to be encoded.*/
static void od_ec_encode_q15(OdEcEnc *enc, unsigned fl, unsigned fh, int32_t s, int32_t nsyms) {
    OdEcWindow l;
    unsigned   r;
    l = enc->low;
    r = enc->rng;
    assert(32768U <= r);
    assert(fh <= fl);
    assert(fl <= 32768U);
    assert(7 - EC_PROB_SHIFT - CDF_SHIFT >= 0);
    const int32_t N = nsyms - 1;
    if (fl < CDF_PROB_TOP) {
        unsigned u;
        unsigned v;
        u = ((r >> 8) * (uint32_t)(fl >> EC_PROB_SHIFT) >> (7 - EC_PROB_SHIFT - CDF_SHIFT)) +
            EC_MIN_PROB * (N - (s - 1));
        v = ((r >> 8) * (uint32_t)(fh >> EC_PROB_SHIFT) >> (7 - EC_PROB_SHIFT - CDF_SHIFT)) +
            EC_MIN_PROB * (N - (s + 0));
        l += r - u;
        r = u - v;
    } else {
        r -= ((r >> 8) * (uint32_t)(fh >> EC_PROB_SHIFT) >> (7 - EC_PROB_SHIFT - CDF_SHIFT)) +
             EC_MIN_PROB * (N - (s + 0));
    }
    od_ec_enc_normalize(enc, l, r);
#if OD_MEASURE_EC_OVERHEAD
    enc->entropy -= OD_LOG2((double)(OD_ICDF(fh) - OD_ICDF(fl)) / CDF_PROB_TOP.);
    enc->nb_symbols++;
#endif
}

/*Encode a single binary value.
val: The value to encode (0 or 1).
f: The probability that the val is one, scaled by 32768.*/
void eb_od_ec_encode_bool_q15(OdEcEnc *enc, int32_t val, unsigned f) {
    OdEcWindow l;
    unsigned   r;
    unsigned   v;
    assert(0 < f);
    assert(f < 32768U);
    l = enc->low;
    r = enc->rng;
    assert(32768U <= r);
    v = ((r >> 8) * (uint32_t)(f >> EC_PROB_SHIFT) >> (7 - EC_PROB_SHIFT));
    v += EC_MIN_PROB;
    if (val) l += r - v;
    r = val ? v : r - v;
    od_ec_enc_normalize(enc, l, r);
#if OD_MEASURE_EC_OVERHEAD
    enc->entropy -= OD_LOG2((double)(val ? f : (32768 - f)) / 32768.);
    enc->nb_symbols++;
#endif
}

/*Encodes a symbol given a cumulative distribution function (CDF) table in Q15.
s: The index of the symbol to encode.
icdf: 32768 minus the CDF, such that symbol s falls in the range
[s > 0 ? (32768 - icdf[s - 1]) : 0, 32768 - icdf[s]).
The values must be monotonically decreasing, and icdf[nsyms - 1] must
be 0.
nsyms: The number of symbols in the alphabet.
This should be at most 16.*/
void eb_od_ec_encode_cdf_q15(OdEcEnc *enc, int32_t s, const uint16_t *icdf, int32_t nsyms) {
    (void)nsyms;
    assert(s >= 0);
    assert(s < nsyms);
    assert(icdf[nsyms - 1] == OD_ICDF(CDF_PROB_TOP));
    od_ec_encode_q15(enc, s > 0 ? icdf[s - 1] : OD_ICDF(0), icdf[s], s, nsyms);
}

uint8_t *eb_od_ec_enc_done(OdEcEnc *enc, uint32_t *nbytes) {
    uint8_t *  out;
    uint32_t   storage;
    uint16_t * buf;
    uint32_t   offs;
    OdEcWindow m;
    OdEcWindow e;
    OdEcWindow l;
    int32_t    c;
    int32_t    s;
    if (enc->error) return NULL;
#if OD_MEASURE_EC_OVERHEAD
    {
        uint32_t tell;
        /* Don't count the 1 bit we lose to raw bits as overhead. */
        tell = eb_od_ec_enc_tell(enc) - 1;
        SVT_ERROR("overhead: %f%%\n", 100 * (tell - enc->entropy) / enc->entropy);
        SVT_ERROR("efficiency: %f bits/symbol\n", (double)tell / enc->nb_symbols);
    }
#endif
    /*We output the minimum number of bits that ensures that the symbols encoded
    thus far will be decoded correctly regardless of the bits that follow.*/
    l = enc->low;
    c = enc->cnt;
    s = 10;
    m = 0x3FFF;
    e = ((l + m) & ~m) | (m + 1);
    s += c;
    offs = enc->offs;
    buf  = enc->precarry_buf;
    if (s > 0) {
        unsigned n;
        storage = enc->precarry_storage;
        if (offs + ((s + 7) >> 3) > storage) {
            storage = storage * 2 + ((s + 7) >> 3);
            buf = realloc(enc->precarry_buf, sizeof(*buf) * storage);
            if (!buf) {
                enc->error = -1;
                return NULL;
            }
            enc->precarry_buf     = buf;
            enc->precarry_storage = storage;
        }
        n = (1 << (c + 16)) - 1;
        do {
            assert(offs < storage);
            buf[offs++] = (uint16_t)(e >> (c + 16));
            e &= n;
            s -= 8;
            c -= 8;
            n >>= 8;
        } while (s > 0);
    }
    /*Make sure there's enough room for the entropy-coded bits.*/
    out     = enc->buf;
    storage = enc->storage;
    c       = OD_MAXI((s + 7) >> 3, 0);
    if (offs + c > storage) {
        storage = offs + c;
        out = realloc(enc->buf, sizeof(*buf) * storage);
        if (!out) {
            enc->error = -1;
            return NULL;
        }
        enc->buf     = out;
        enc->storage = storage;
    }
    *nbytes = offs;
    /*Perform carry propagation.*/
    assert(offs <= storage);
    out = out + storage - offs;
    c   = 0;
    while (offs > 0) {
        offs--;
        c         = buf[offs] + c;
        out[offs] = (uint8_t)c;
        c >>= 8;
    }
    /*Note: Unless there's an allocation error, if you keep encoding into the
    current buffer and call this function again later, everything will work
    just fine (you won't get a new packet out, but you will get a single
    buffer with the new data appended to the old).
    However, this function is O(N) where N is the amount of data coded so far,
    so calling it more than once for a given packet is a bad idea.*/
    return out;
}

/*Returns the number of bits "used" by the encoded symbols so far.
This same number can be computed in either the encoder or the decoder, and is
suitable for making coding decisions.
Warning: The value returned by this function can decrease compared to an
earlier call, even after encoding more data, if there is an encoding error
(i.e., a failure to allocate enough space for the output buffer).
Return: The number of bits.
This will always be slightly larger than the exact value (e.g., all
rounding error is in the positive direction).*/
int32_t eb_od_ec_enc_tell(const OdEcEnc *enc) {
    /*The 10 here counteracts the offset of -9 baked into cnt, and adds 1 extra
    bit, which we reserve for terminating the stream.*/
    return (enc->cnt + 10) + enc->offs * 8;
}
/********************************************************************************************************************************/
/********************************************************************************************************************************/
/********************************************************************************************************************************/
