/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/
#ifndef EbDecBitstream_h
#define EbDecBitstream_h

// Defines the maximum number of bits in a bitstream word
#define WORD_SIZE         32

// Twice the WORD_SIZE
#define DBL_WORD_SIZE     (2 * (WORD_SIZE))

#define SHL(x,y) (((y) < 32) ? ((x) << (y)) : 0)
#define SHR(x,y) (((y) < 32) ? ((x) >> (y)) : 0)

#define TO_BIG_ENDIAN(x)       ((x << 24))                |   \
                            ((x & 0x0000ff00) << 8)    |   \
                            ((x & 0x00ff0000) >> 8)    |   \
                            ((uint32_t)x >> 24);

static const uint8_t kLeb128ByteMask = 0x7f; // Binary: 01111111

/*
 * Bitstream structure
 */
typedef struct
{
    /* Bitstream buffer base pointer */
    uint8_t *buf_base;

    /* Bitstream bit offset in current word. Value between 0 and 31 */
    uint32_t bit_ofst;

    /* Current bitstream buffer pointer */
    uint32_t *buf;

    /* Current word */
    uint32_t cur_word;

    /* Next word */
    uint32_t nxt_word;

    /* Max address for bitstream */
    uint8_t *buf_max;
} bitstrm_t;

// Get m_cnt number of bits and update bffer pointers and offset.
#define GET_BITS(bits, m_pu4_buf, bit_ofst,                   \
                          cur_word,nxt_word, m_cnt)           \
{                                                             \
    bits = (cur_word << bit_ofst)                             \
                             >> (WORD_SIZE - m_cnt);          \
    bit_ofst += m_cnt;                                        \
    if(bit_ofst > WORD_SIZE)                                  \
    {                                                         \
        bits |= SHR(nxt_word,                                 \
                     (DBL_WORD_SIZE - bit_ofst));             \
    }                                                         \
                                                              \
    if( bit_ofst >=   WORD_SIZE )                             \
    {                                                         \
        uint32_t pu4_word_tmp;                                \
        cur_word  = nxt_word;                                 \
        /* Getting the next word */                           \
        pu4_word_tmp = *(m_pu4_buf++);                        \
                                                              \
        bit_ofst -= WORD_SIZE;                                \
        /* Swapping little endian to big endian conversion*/  \
        nxt_word  = TO_BIG_ENDIAN(pu4_word_tmp);              \
    }                                                         \
}

void dec_bits_init(bitstrm_t *bs, const uint8_t *data, size_t u4_numbytes);

uint32_t dec_get_bits_uvlc(bitstrm_t *bs);
uint32_t dec_get_bits(bitstrm_t *bs, uint32_t numbits);
void dec_get_bits_leb128(bitstrm_t *bs, size_t available, size_t *value,
                    size_t *length);
uint32_t dec_get_bits_ns(bitstrm_t *bs, uint32_t n);
int32_t dec_get_bits_su(bitstrm_t *bs, uint32_t n);
uint32_t dec_get_bits_le(bitstrm_t *bs, uint32_t n);

uint32_t get_position(bitstrm_t *bs);
uint8_t * get_bitsteam_buf(bitstrm_t *bs);

#endif // EbDecBitstream_h
