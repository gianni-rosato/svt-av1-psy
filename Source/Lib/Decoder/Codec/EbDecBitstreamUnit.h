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

int od_ec_decode_bool_q15(OdEcDec *dec, unsigned f);
int od_ec_decode_cdf_q15(OdEcDec *dec, const uint16_t *cdf, int nsyms);

/********************************************************************************************************************************/
/********************************************************************************************************************************/
/********************************************************************************************************************************/
/********************************************************************************************************************************/
// daalaboolreader.h from AOM

typedef struct DaalaReader {
    const uint8_t *buffer;
    const uint8_t *buffer_end;
    OdEcDec        ec;
    uint8_t        allow_update_cdf;
} DaalaReader_t;

int aom_daala_reader_init(DaalaReader_t *r, const uint8_t *buffer, int size);

const uint8_t *aom_daala_reader_find_begin(DaalaReader_t *r);

const uint8_t *aom_daala_reader_find_end(DaalaReader_t *r);

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
            SVT_ERROR("\n *** [bit] cdf error, frame_idx_r %d cdf {%d, %d} ref_cdf {%d",
                      frame_idx,
                      p,
                      32767,
                      ref_cdf[0]);
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
                if (cdf[i] != ref_cdf[i]) cdf_error = 1;
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
