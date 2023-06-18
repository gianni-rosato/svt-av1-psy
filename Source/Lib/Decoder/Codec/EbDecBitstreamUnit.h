/*
* Copyright(c) 2019 Intel Corporation
* Copyright (c) 2016, Alliance for Open Media. All rights reserved
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

#ifndef EbDecBitstreamUnit_h
#define EbDecBitstreamUnit_h

#include "EbCabacContextModel.h"
#include "EbBitstreamUnit.h"
//Added this EbBitstreamUnit.h because OdEcWindow is defined in it, but
//we also defining it, so it leads to warning,  so i commented our defination & added EbBitstreamUnit.h file.

#ifdef __cplusplus
extern "C" {
#endif

#define CONFIG_BITSTREAM_DEBUG 0

/********************************************************************************************************************************/
/********************************************************************************************************************************/
/********************************************************************************************************************************/
/********************************************************************************************************************************/
#if (CHAR_BIT != 8)
#undef CHAR_BIT
#define CHAR_BIT 8 /* number of bits in a char */
#endif
// entcode.h from AOM

#define EC_PROB_SHIFT 6
#define EC_MIN_PROB 4 // must be <= (1<<EC_PROB_SHIFT)/16

/*OPT: OdEcWindow must be at least 32 bits, but if you have fast arithmetic
   on a larger type, you can speed up the decoder by using it here.*/
//typedef uint32_t OdEcWindow;

/*The size in bits of OdEcWindow.*/
//#define OD_EC_WINDOW_SIZE ((int)sizeof(OdEcWindow) * CHAR_BIT)

/********************************************************************************************************************************/
/********************************************************************************************************************************/
/********************************************************************************************************************************/
/********************************************************************************************************************************/
// prob.h from AOM
//typedef uint16_t aom_cdf_prob;
#define ENABLE_ENTROPY_TRACE 0
#define EXTRA_DUMP 0
#if ENABLE_ENTROPY_TRACE
#define ENTROPY_TRACE_FILE_BASED 1
#define FRAME_LEVEL_TRACE 1
#include <stdio.h>
extern FILE *temp_fp;
extern int   enable_dump;
#endif

#define CDF_PROB_BITS 15
#define CDF_PROB_TOP (1 << CDF_PROB_BITS)
#define CDF_INIT_TOP 32768
#define CDF_SHIFT (15 - CDF_PROB_BITS)
/*The value stored in an iCDF is CDF_PROB_TOP minus the actual cumulative
  probability (an "inverse" CDF).
  This function converts from one representation to the other (and is its own
  inverse).*/
#define AOM_ICDF(x) (CDF_PROB_TOP - (x))

static INLINE void dec_update_cdf(AomCdfProb *cdf, int8_t val, int nsymbs) {
    int rate;
    int i, tmp;

    static const int nsymbs2speed[17] = {0, 0, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
    assert(nsymbs < 17);
    rate = 3 + (cdf[nsymbs] > 15) + (cdf[nsymbs] > 31) + nsymbs2speed[nsymbs]; // + get_msb(nsymbs);
    tmp  = AOM_ICDF(0);

    // Single loop (faster)
    for (i = 0; i < nsymbs - 1; ++i) {
        tmp = (i == val) ? 0 : tmp;
        if (tmp < cdf[i]) {
            cdf[i] -= (AomCdfProb)((cdf[i] - tmp) >> rate);
        } else
            cdf[i] += (AomCdfProb)((tmp - cdf[i]) >> rate);
    }
    cdf[nsymbs] += (cdf[nsymbs] < 32);
}

/********************************************************************************************************************************/
/********************************************************************************************************************************/
/********************************************************************************************************************************/
/********************************************************************************************************************************/
// entdec.h from AOM

/*The entropy decoder context.*/
typedef struct OdEcDec {
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
    decoder only uses the top 16 bits of the window to decode the next symbol.
    As we shift up during renormalization, if we don't have enough bits left in
    the window to fill the top 16, we'll read in more bits of the coded
    value.*/
    OdEcWindow dif;
    /*The number of values in the current range.*/
    uint16_t rng;
    /*The number of bits of data in the current value.*/
    int16_t cnt;
} OdEcDec;

static void od_ec_dec_refill(OdEcDec *dec);

/*Takes updated dif and range values, renormalizes them so that
   32768 <= rng < 65536 (reading more bytes from the stream into dif if
   necessary), and stores them back in the decoder context.
  dif: The new value of dif.
  rng: The new value of the range.
  ret: The value to return.
  Return: ret.
          This allows the compiler to jump to this function via a tail-call.*/
static int od_ec_dec_normalize(OdEcDec *dec, OdEcWindow dif, unsigned rng, int ret) {
    int d;
    assert(rng <= 65535U);
    /*The number of leading zeros in the 16-bit binary representation of rng.*/
    d = 16 - OD_ILOG_NZ(rng);
    /*d bits in dec->dif are consumed.*/
    dec->cnt -= d;
    /*This is equivalent to shifting in 1's instead of 0's.*/
    dec->dif = ((dif + 1) << d) - 1;
    dec->rng = rng << d;
    if (dec->cnt < 0)
        od_ec_dec_refill(dec);
    return ret;
}

/*Decode a single binary value.
  f: The probability that the bit is one, scaled by 32768.
  Return: The value decoded (0 or 1).*/
static int od_ec_decode_bool_q15(OdEcDec *dec, unsigned f) {
    OdEcWindow dif;
    OdEcWindow vw;
    unsigned   r;
    unsigned   r_new;
    unsigned   v;
    int        ret;
    assert(0 < f);
    assert(f < 32768U);
    dif = dec->dif;
    r   = dec->rng;
    assert(dif >> (OD_EC_WINDOW_SIZE - 16) < r);
    assert(32768U <= r);
    v = ((r >> 8) * (uint32_t)(f >> EC_PROB_SHIFT) >> (7 - EC_PROB_SHIFT));
    v += EC_MIN_PROB;
    vw    = (OdEcWindow)v << (OD_EC_WINDOW_SIZE - 16);
    ret   = 1;
    r_new = v;
    if (dif >= vw) {
        r_new = r - v;
        dif -= vw;
        ret = 0;
    }
    return od_ec_dec_normalize(dec, dif, r_new, ret);
}

/*Decodes a symbol given an inverse cumulative distribution function (CDF)
   table in Q15.
  icdf: CDF_PROB_TOP minus the CDF, such that symbol s falls in the range
         [s > 0 ? (CDF_PROB_TOP - icdf[s - 1]) : 0, CDF_PROB_TOP - icdf[s]).
        The values must be monotonically non-increasing, and icdf[nsyms - 1]
         must be 0.
  nsyms: The number of symbols in the alphabet.
         This should be at most 16.
  Return: The decoded symbol s.*/
static int od_ec_decode_cdf_q15(OdEcDec *dec, const uint16_t *icdf, int nsyms) {
    OdEcWindow dif;
    unsigned   r;
    unsigned   c;
    unsigned   u;
    unsigned   v;
    int        ret;
    (void)nsyms;
    dif         = dec->dif;
    r           = dec->rng;
    const int N = nsyms - 1;

    assert(dif >> (OD_EC_WINDOW_SIZE - 16) < r);
    assert(icdf[nsyms - 1] == OD_ICDF(CDF_PROB_TOP));
    assert(32768U <= r);
    assert(7 - EC_PROB_SHIFT - CDF_SHIFT >= 0);
    c   = (unsigned)(dif >> (OD_EC_WINDOW_SIZE - 16));
    v   = r;
    ret = -1;
    do {
        u = v;
        v = ((r >> 8) * (uint32_t)(icdf[++ret] >> EC_PROB_SHIFT) >> (7 - EC_PROB_SHIFT - CDF_SHIFT));
        v += EC_MIN_PROB * (N - ret);
    } while (c < v);
    assert(v < u);
    assert(u <= r);
    r = u - v;
    dif -= (OdEcWindow)v << (OD_EC_WINDOW_SIZE - 16);
    return od_ec_dec_normalize(dec, dif, r, ret);
}

/********************************************************************************************************************************/
/********************************************************************************************************************************/
/********************************************************************************************************************************/
/********************************************************************************************************************************/
// daalaboolreader.h from AOM

/*This is meant to be a large, positive constant that can still be efficiently
   loaded as an immediate (on platforms like ARM, for example).
  Even relatively modest values like 100 would work fine.*/
#define OD_EC_LOTS_OF_BITS (0x4000)

/*The return value of od_ec_dec_tell does not change across an od_ec_dec_refill
   call.*/
static void od_ec_dec_refill(OdEcDec *dec) {
    int                  s;
    OdEcWindow           dif;
    int16_t              cnt;
    const unsigned char *bptr;
    const unsigned char *end;
    dif  = dec->dif;
    cnt  = dec->cnt;
    bptr = dec->bptr;
    end  = dec->end;
    s    = OD_EC_WINDOW_SIZE - 9 - (cnt + 15);
    for (; s >= 0 && bptr < end; s -= 8, bptr++) {
        /*Each time a byte is inserted into the window (dif), bptr advances and cnt
       is incremented by 8, so the total number of consumed bits (the return
       value of od_ec_dec_tell) does not change.*/
        assert(s <= OD_EC_WINDOW_SIZE - 8);
        dif ^= (OdEcWindow)bptr[0] << s;
        cnt += 8;
    }
    if (bptr >= end) {
        /*We've reached the end of the buffer. It is perfectly valid for us to need
       to fill the window with additional bits past the end of the buffer (and
       this happens in normal operation). These bits should all just be taken
       as zero. But we cannot increment bptr past 'end' (this is undefined
       behavior), so we start to increment dec->tell_offs. We also don't want
       to keep testing bptr against 'end', so we set cnt to OD_EC_LOTS_OF_BITS
       and adjust dec->tell_offs so that the total number of unconsumed bits in
       the window (dec->cnt - dec->tell_offs) does not change. This effectively
       puts lots of zero bits into the window, and means we won't try to refill
       it from the buffer for a very long time (at which point we'll put lots
       of zero bits into the window again).*/
        dec->tell_offs += OD_EC_LOTS_OF_BITS - cnt;
        cnt = OD_EC_LOTS_OF_BITS;
    }
    dec->dif  = dif;
    dec->cnt  = cnt;
    dec->bptr = bptr;
}

/*Initializes the decoder.
  buf: The input buffer to use.
  storage: The size in bytes of the input buffer.*/
static void od_ec_dec_init(OdEcDec *dec, const unsigned char *buf, uint32_t storage) {
    dec->buf       = buf;
    dec->tell_offs = 10 - (OD_EC_WINDOW_SIZE - 8);
    dec->end       = buf + storage;
    dec->bptr      = buf;
    dec->dif       = ((OdEcWindow)1 << (OD_EC_WINDOW_SIZE - 1)) - 1;
    dec->rng       = 0x8000;
    dec->cnt       = -15;
    od_ec_dec_refill(dec);
}

typedef struct DaalaReader {
    const uint8_t *buffer;
    const uint8_t *buffer_end;
    OdEcDec        ec;
    uint8_t        allow_update_cdf;
} DaalaReader_t;

static inline int svt_aom_daala_reader_init(DaalaReader_t *r, const uint8_t *buffer, int size) {
    if (size && !buffer)
        return 1;
    r->buffer_end = buffer + size;
    r->buffer     = buffer;
    od_ec_dec_init(&r->ec, buffer, size);
#if CONFIG_ACCOUNTING
    r->accounting = NULL;
#endif
    return 0;
}

static inline const uint8_t *svt_aom_daala_reader_find_begin(DaalaReader_t *r) { return r->buffer; }

static inline const uint8_t *svt_aom_daala_reader_find_end(DaalaReader_t *r) { return r->buffer_end; }

static INLINE int aom_daala_read(DaalaReader_t *r, int prob) {
    int bit;
    int p = (0x7FFFFF - (prob << 15) + prob) >> 8;

    bit = od_ec_decode_bool_q15(&r->ec, p);
#if CONFIG_BITSTREAM_DEBUG
    {
        int        i;
        int        ref_bit, ref_nsymbs;
        AomCdfProb ref_cdf[16];
        const int  queue_r   = bitstream_queue_get_read();
        const int  frame_idx = bitstream_queue_get_frame_read();
        bitstream_queue_pop(&ref_bit, ref_cdf, &ref_nsymbs);
        if (ref_nsymbs != 2) {
            SVT_ERROR(
                "\n *** [bit] nsymbs error, frame_idx_r %d nsymbs %d ref_nsymbs "
                "%d queue_r %d\n",
                frame_idx,
                2,
                ref_nsymbs,
                queue_r);
            assert(0);
        }
        if ((ref_nsymbs != 2) || (ref_cdf[0] != (AomCdfProb)p) || (ref_cdf[1] != 32767)) {
            SVT_ERROR(
                "\n *** [bit] cdf error, frame_idx_r %d cdf {%d, %d} ref_cdf {%d", frame_idx, p, 32767, ref_cdf[0]);
            for (i = 1; i < ref_nsymbs; ++i) SVT_ERROR(", %d", ref_cdf[i]);
            SVT_ERROR("} queue_r %d\n", queue_r);
            assert(0);
        }
        if (bit != ref_bit) {
            SVT_ERROR(
                "\n *** [bit] symb error, frame_idx_r %d symb %d ref_symb %d "
                "queue_r %d\n",
                frame_idx,
                bit,
                ref_bit,
                queue_r);
            assert(0);
        }
    }
#endif
#if ENABLE_ENTROPY_TRACE
#if ENTROPY_TRACE_FILE_BASED
    if (enable_dump) {
        assert(temp_fp);
        fprintf(temp_fp, "\n *** p %d \t", p);
        fprintf(temp_fp, "symb : %d \t", bit);
        fflush(temp_fp);
    }
#else
    if (enable_dump) {
        SVT_LOG("\n *** p %d \t", p);
        SVT_LOG("symb : %d \t", bit);
        fflush(stdout);
    }
#endif
#endif
    return bit;
}

static INLINE int daala_read_symbol(DaalaReader_t *r, const AomCdfProb *cdf, int nsymbs) {
    int symb;
    assert(cdf != NULL);
    symb = od_ec_decode_cdf_q15(&r->ec, cdf, nsymbs);
#if CONFIG_BITSTREAM_DEBUG
    {
        int        i;
        int        cdf_error = 0;
        int        ref_symb, ref_nsymbs;
        AomCdfProb ref_cdf[16];
        const int  queue_r   = bitstream_queue_get_read();
        const int  frame_idx = bitstream_queue_get_frame_read();
        bitstream_queue_pop(&ref_symb, ref_cdf, &ref_nsymbs);
        if (nsymbs != ref_nsymbs) {
            SVT_ERROR(
                "\n *** nsymbs error, frame_idx_r %d nsymbs %d ref_nsymbs %d "
                "queue_r %d\n",
                frame_idx,
                nsymbs,
                ref_nsymbs,
                queue_r);
            cdf_error = 0;
            assert(0);
        } else {
            for (i = 0; i < nsymbs; ++i)
                if (cdf[i] != ref_cdf[i])
                    cdf_error = 1;
        }
        if (cdf_error) {
            SVT_ERROR("\n *** cdf error, frame_idx_r %d cdf {%d", frame_idx, cdf[0]);
            for (i = 1; i < nsymbs; ++i) SVT_ERROR(", %d", cdf[i]);
            SVT_ERROR("} ref_cdf {%d", ref_cdf[0]);
            for (i = 1; i < ref_nsymbs; ++i) SVT_ERROR(", %d", ref_cdf[i]);
            SVT_ERROR("} queue_r %d\n", queue_r);
            assert(0);
        }
        if (symb != ref_symb) {
            fprintf(stderr,
                    "\n *** symb error, frame_idx_r %d symb %d ref_symb %d queue_r %d\n",
                    frame_idx,
                    symb,
                    ref_symb,
                    queue_r);
            assert(0);
        }
    }
#endif
#if ENABLE_ENTROPY_TRACE
#if ENTROPY_TRACE_FILE_BASED
    if (enable_dump) {
        fprintf(temp_fp, "\n *** nsymbs %d \t", nsymbs);
        for (int i = 0; i < nsymbs; ++i) fprintf(temp_fp, "cdf[%d] : %d \t", i, cdf[i]);
        fprintf(temp_fp, "symb : %d \t", symb);
        fflush(temp_fp);
    }
#else
    if (enable_dump) {
        SVT_LOG("\n *** nsymbs %d \t", nsymbs);
        for (int i = 0; i < nsymbs; ++i) SVT_LOG("cdf[%d] : %d \t", i, cdf[i]);
        SVT_LOG("symb : %d \t", symb);
        fflush(stdout);
    }
#endif
#endif
    return symb;
}

#ifdef __cplusplus
}
#endif
#endif // EbDecBitstreamUnit_h
