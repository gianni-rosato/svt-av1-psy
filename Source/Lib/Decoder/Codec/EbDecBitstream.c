/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbDefinitions.h"
#include "EbBitstreamUnit.h"

#include "EbDecBitstream.h"

/* Function used for bitstream structure initialization
Assumes data is aligned to 4 bytes. If not aligned  then all bitstream
accesses will be unaligned and hence  costlier. Since this is codec memory that
holds emulation prevented data, assumption of aligned to 4 bytes is valid */
void dec_bits_init(bitstrm_t *bs, const uint8_t *data, size_t numbytes) {
    uint32_t cur_word;
    uint32_t nxt_word;
    uint32_t temp;
    uint32_t *buf;
    buf = (uint32_t *)data;
    temp = *buf++;
    cur_word = TO_BIG_ENDIAN(temp);
    temp = *buf++;
    nxt_word = TO_BIG_ENDIAN(temp);
    bs->bit_ofst = 0;
    bs->buf_base = (uint8_t *)data;
    bs->buf = buf;
    bs->cur_word = cur_word;
    bs->nxt_word = nxt_word;
    bs->buf_max = (uint8_t *)data + numbytes + 8;
    return;
}

/* Reads next numbits number of bits from the bitstream  this updates the
bitstream offset and consumes the bits. Section: 4.10.2 -> f(n) */
uint32_t dec_get_bits(bitstrm_t *bs, uint32_t numbits) {
    uint32_t bits_read;
    GET_BITS(bits_read,
        bs->buf,
        bs->bit_ofst,
        bs->cur_word,
        bs->nxt_word,
        numbits);
    return bits_read;
}

/* Get unsigned integer represented by a variable number of little-endian bytes */
void dec_get_bits_leb128(bitstrm_t *bs, size_t available, size_t *value,
    size_t *length)
{
    (void)available;
    *value = 0;
    *length = 0;
    int i;
    uint32_t leb128_byte = 0;

    for (i = 0; i < 8; i++) {
        leb128_byte = dec_get_bits(bs, 8);
        *value |= (((uint64_t)leb128_byte & 0x7f) << (i * 7));
        *length += 1;
        if (!(leb128_byte & 0x80))
            break;
    }
}

/* Get variable length unsigned n-bit number appearing directly in the bitstream */
uint32_t dec_get_bits_uvlc(bitstrm_t *bs) {
    int leading_zeros = 0;
    while (leading_zeros < 32 && !dec_get_bits(bs, 1)) ++leading_zeros;
    // Maximum 32 bits.
    if (leading_zeros == 32) return UINT32_MAX;
    const uint32_t base = (1u << leading_zeros) - 1;
    const uint32_t value = dec_get_bits(bs, leading_zeros);
    return base + value;
}

/* Unsigned encoded integer with maximum number of values n */
uint32_t dec_get_bits_ns(bitstrm_t *bs, uint32_t n) {
    if (n <= 1) return 0;
    int w = get_msb(n) + 1; //w = FloorLog2(n) + 1
    int m = (1 << w) - n;
    int v = dec_get_bits(bs, w - 1);
    if (v < m)
        return v;
    return (v << 1) - m + dec_get_bits(bs, 1);
}

/* Signed integer converted from an n bits unsigned integer in the bitstream */
int32_t dec_get_bits_su(bitstrm_t *bs, uint32_t n) {
    int value = dec_get_bits(bs, n);
    int signMask = 1 << (n - 1);
    if (value & signMask)
        value = value - 2 * signMask;
    return value;
}

/* Unsigned little-endian n-byte number appearing directly in the bitstream */
uint32_t dec_get_bits_le(bitstrm_t *bs, uint32_t n) {
    uint32_t t = 0, byte;
        for (uint32_t i = 0; i < n; i++) {
            byte = dec_get_bits(bs, 8);
            t += (byte << (i * 8));
        }
        return t;
}

uint32_t get_position(bitstrm_t *bs) {
    return (uint32_t)( (( ((uint8_t *)bs->buf) - bs->buf_base) * 8) - WORD_SIZE/*nxt_word*/
            - (WORD_SIZE - bs->bit_ofst)/*cur_word*/ );
}

uint8_t * get_bitsteam_buf(bitstrm_t *bs) {
    uint8_t *bitsteam_buf = (uint8_t *)bs->buf;
    bitsteam_buf -= ((WORD_SIZE/*nxt_word*/ >> 3) + ((WORD_SIZE - bs->bit_ofst)/*cur_word*/ >> 3));

    assert(bitsteam_buf == (bs->buf_base + (get_position(bs) >> 3) ) );

    return bitsteam_buf;
}
