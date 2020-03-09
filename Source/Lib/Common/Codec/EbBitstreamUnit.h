/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbBitstreamUnit_h
#define EbBitstreamUnit_h

#include "EbObject.h"
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif

#ifdef _MSC_VER
#if defined(_M_X64) || defined(_M_IX86)
#include <intrin.h>
#define USE_MSC_INTRINSICS
#endif
#endif

// Bistream Slice Buffer Size
#define EB_BITSTREAM_SLICE_BUFFER_SIZE 0x300000
#define SLICE_HEADER_COUNT 256

/**********************************
 * Bitstream Unit Types
 **********************************/
typedef struct OutputBitstreamUnit {
    EbDctor  dctor;
    uint32_t size; // allocated buffer size
    uint8_t *buffer_begin_av1; // the byte buffer
    uint8_t *buffer_av1; // the byte buffer
} OutputBitstreamUnit;

/**********************************
     * Extern Function Declarations
     **********************************/
extern EbErrorType output_bitstream_unit_ctor(OutputBitstreamUnit *bitstream_ptr,
                                              uint32_t             buffer_size);

extern EbErrorType output_bitstream_reset(OutputBitstreamUnit *bitstream_ptr);

/********************************************************************************************************************************/
/********************************************************************************************************************************/
/********************************************************************************************************************************/
#include "EbCabacContextModel.h"
/********************************************************************************************************************************/
// bitops.h
// These versions of get_msb() are only valid when n != 0 because all
// of the optimized versions are undefined when n == 0:
// https://gcc.gnu.org/onlinedocs/gcc/Other-Builtins.html

// use GNU builtins where available.
#if defined(__GNUC__) && ((__GNUC__ == 3 && __GNUC_MINOR__ >= 4) || __GNUC__ >= 4)
static INLINE int32_t get_msb(uint32_t n) {
    assert(n != 0);
    return 31 ^ __builtin_clz(n);
}
#elif defined(USE_MSC_INTRINSICS)
#pragma intrinsic(_BitScanReverse)

static INLINE int32_t get_msb(uint32_t n) {
    unsigned long first_set_bit;
    assert(n != 0);
    _BitScanReverse(&first_set_bit, n);
    return first_set_bit;
}
#undef USE_MSC_INTRINSICS
#else
// Returns (int32_t)floor(log2(n)). n must be > 0.
/*static*/ INLINE int32_t get_msb(uint32_t n) {
    int32_t  log   = 0;
    uint32_t value = n;
    int32_t  i;

    assert(n != 0);

    for (i = 4; i >= 0; --i) {
        const int32_t  shift = (1 << i);
        const uint32_t x     = value >> shift;
        if (x != 0) {
            value = x;
            log += shift;
        }
    }
    return log;
}
#endif
/********************************************************************************************************************************/
//odintrin.h

#define OD_DIVU_DMAX (1024)

extern uint32_t od_divu_small_consts[OD_DIVU_DMAX][2];

/*Translate unsigned division by small divisors into multiplications.*/
#define OD_DIVU_SMALL(_x, _d)                                                                    \
    ((uint32_t)(                                                                                 \
         (od_divu_small_consts[(_d)-1][0] * (uint64_t)(_x) + od_divu_small_consts[(_d)-1][1]) >> \
         32) >>                                                                                  \
     (OD_ILOG_NZ(_d) - 1))

#define OD_DIVU(_x, _d) (((_d) < OD_DIVU_DMAX) ? (OD_DIVU_SMALL((_x), (_d))) : ((_x) / (_d)))

#define OD_MINI MIN
#define OD_MAXI MAX
#define OD_CLAMPI(min, val, max) (OD_MAXI(min, OD_MINI(val, max)))

#define OD_CLZ0 (1)
#define OD_CLZ(x) (-get_msb(x))
#define OD_ILOG_NZ(x) (OD_CLZ0 - OD_CLZ(x))

/*Enable special features for gcc and compatible compilers.*/
#if defined(__GNUC__) && defined(__GNUC_MINOR__) && defined(__GNUC_PATCHLEVEL__)
#define OD_GNUC_PREREQ(maj, min, pat)                                  \
    ((__GNUC__ << 16) + (__GNUC_MINOR__ << 8) + __GNUC_PATCHLEVEL__ >= \
     ((maj) << 16) + ((min) << 8) + pat) // NOLINT
#else
#define OD_GNUC_PREREQ(maj, min, pat) (0)
#endif

#if OD_GNUC_PREREQ(3, 4, 0)
#define OD_WARN_UNUSED_RESULT __attribute__((__warn_unused_result__))
#else
#define OD_WARN_UNUSED_RESULT
#endif

#if OD_GNUC_PREREQ(3, 4, 0)
#define OD_ARG_NONNULL(x) __attribute__((__nonnull__(x)))
#else
#define OD_ARG_NONNULL(x)
#endif

/** Copy n elements of memory from src to dst. The 0* term provides
compile-time type checking  */
#if !defined(OVERRIDE_OD_COPY)
#define OD_COPY(dst, src, n) (memcpy((dst), (src), sizeof(*(dst)) * (n) + 0 * ((dst) - (src))))
#endif

/** Copy n elements of memory from src to dst, allowing overlapping regions.
The 0* term provides compile-time type checking */
#if !defined(OVERRIDE_OD_MOVE)
#define OD_MOVE(dst, src, n) (memmove((dst), (src), sizeof(*(dst)) * (n) + 0 * ((dst) - (src))))
#endif

/*All of these macros should expect floats as arguments.*/
#define OD_SIGNMASK(a) (-((a) < 0))
#define OD_FLIPSIGNI(a, b) (((a) + OD_SIGNMASK(b)) ^ OD_SIGNMASK(b))

/********************************************************************************************************************************/
//entcode.h
#define EC_PROB_SHIFT 6
#define EC_MIN_PROB 4 // must be <= (1<<EC_PROB_SHIFT)/16

/*OPT: OdEcWindow must be at least 32 bits, but if you have fast arithmetic
on a larger type, you can speed up the decoder by using it here.*/
typedef uint32_t OdEcWindow;

#define OD_EC_WINDOW_SIZE ((int32_t)sizeof(OdEcWindow) * CHAR_BIT)

/*The resolution of fractional-precision bit usage measurements, i.e.,
    3 => 1/8th bits.*/
#define OD_BITRES (3)

#define OD_ICDF AOM_ICDF

/********************************************************************************************************************************/
//entenc.h
typedef struct OdEcEnc OdEcEnc;

#define OD_MEASURE_EC_OVERHEAD (0)

/*The entropy encoder context.*/
struct OdEcEnc {
    /*Buffered output.
        This contains only the raw bits until the final call to eb_od_ec_enc_done(),
        where all the arithmetic-coded data gets prepended to it.*/
    uint8_t *buf;
    /*The size of the buffer.*/
    uint32_t storage;
    /*The offset at which the last byte containing raw bits was written.*/

    /*A buffer for output bytes with their associated carry flags.*/
    uint16_t *precarry_buf;
    /*The size of the pre-carry buffer.*/
    uint32_t precarry_storage;
    /*The offset at which the next entropy-coded byte will be written.*/
    uint32_t offs;
    /*The low end of the current range.*/
    OdEcWindow low;
    /*The number of values in the current range.*/
    uint16_t rng;
    /*The number of bits of data in the current value.*/
    int16_t cnt;
    /*Nonzero if an error occurred.*/
    int32_t error;
#if OD_MEASURE_EC_OVERHEAD
    double  entropy;
    int32_t nb_symbols;
#endif
};

/*See entenc.c for further documentation.*/

void eb_od_ec_enc_init(OdEcEnc *enc, uint32_t size) OD_ARG_NONNULL(1);
void eb_od_ec_enc_reset(OdEcEnc *enc) OD_ARG_NONNULL(1);
void eb_od_ec_enc_clear(OdEcEnc *enc) OD_ARG_NONNULL(1);

void eb_od_ec_encode_bool_q15(OdEcEnc *enc, int32_t val, unsigned f_q15) OD_ARG_NONNULL(1);
void eb_od_ec_encode_cdf_q15(OdEcEnc *enc, int32_t s, const uint16_t *cdf, int32_t nsyms)
    OD_ARG_NONNULL(1) OD_ARG_NONNULL(3);

void od_ec_enc_bits(OdEcEnc *enc, uint32_t fl, unsigned ftb) OD_ARG_NONNULL(1);

OD_WARN_UNUSED_RESULT uint8_t *eb_od_ec_enc_done(OdEcEnc *enc, uint32_t *nbytes) OD_ARG_NONNULL(1)
    OD_ARG_NONNULL(2);

OD_WARN_UNUSED_RESULT int32_t eb_od_ec_enc_tell(const OdEcEnc *enc) OD_ARG_NONNULL(1);

/********************************************************************************************************************************/
//daalaboolwriter.h
struct DaalaWriter {
    uint32_t pos;
    uint8_t *buffer;
    OdEcEnc  ec;
    uint8_t  allow_update_cdf;
};

typedef struct DaalaWriter DaalaWriter;

void    eb_aom_daala_start_encode(DaalaWriter *w, uint8_t *buffer);
int32_t eb_aom_daala_stop_encode(DaalaWriter *w);

static INLINE void aom_daala_write(DaalaWriter *w, int32_t bit, int32_t prob) {
    int32_t p = (0x7FFFFF - (prob << 15) + prob) >> 8;
#if CONFIG_BITSTREAM_DEBUG
    AomCdfProb cdf[2] = {(AomCdfProb)p, 32767};
    bitstream_queue_push(bit, cdf, 2);
#endif
    eb_od_ec_encode_bool_q15(&w->ec, bit, p);
}

static INLINE void daala_write_symbol(DaalaWriter *w, int32_t symb, const AomCdfProb *cdf,
                                      int32_t nsymbs) {
#if CONFIG_BITSTREAM_DEBUG
    bitstream_queue_push(symb, cdf, nsymbs);
#endif
    eb_od_ec_encode_cdf_q15(&w->ec, symb, cdf, nsymbs);
}

/********************************************************************************************************************************/
// bitwriter.h
typedef struct DaalaWriter AomWriter;

static INLINE void aom_start_encode(AomWriter *bc, uint8_t *buffer) {
    eb_aom_daala_start_encode(bc, buffer);
}

static INLINE int32_t aom_stop_encode(AomWriter *bc) { return eb_aom_daala_stop_encode(bc); }

static INLINE void aom_write(AomWriter *br, int32_t bit, int32_t probability) {
    aom_daala_write(br, bit, probability);
}

static INLINE void aom_write_bit(AomWriter *w, int32_t bit) {
    aom_write(w, bit, 128); // aom_prob_half
}

static INLINE void aom_write_literal(AomWriter *w, int32_t data, int32_t bits) {
    int32_t bit;

    for (bit = bits - 1; bit >= 0; bit--) aom_write_bit(w, 1 & (data >> bit));
}

static INLINE void aom_write_cdf(AomWriter *w, int32_t symb, const AomCdfProb *cdf,
                                 int32_t nsymbs) {
    daala_write_symbol(w, symb, cdf, nsymbs);
}

static INLINE void aom_write_symbol(AomWriter *w, int32_t symb, AomCdfProb *cdf, int32_t nsymbs) {
    aom_write_cdf(w, symb, cdf, nsymbs);
    if (w->allow_update_cdf) update_cdf(cdf, symb, nsymbs);
}

/********************************************************************************************************************************/
/********************************************************************************************************************************/
#ifdef __cplusplus
}
#endif

#endif // EbBitstreamUnit_h
