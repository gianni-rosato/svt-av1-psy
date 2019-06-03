/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

/*
 * This is the header file for the MD5 message-digest algorithm.
 * The algorithm is due to Ron Rivest.  This code was
 * written by Colin Plumb in 1993, no copyright is claimed.
 * This code is in the public domain; do with it what you wish.
 *
 * Equivalent code is available from RSA Data Security, Inc.
 * This code has been tested against that, and is equivalent,
 * except that you don't need to include two pages of legalese
 * with every copy.
 *
 * To compute the message digest of a chunk of bytes, declare an
 * MD5Context structure, pass it to MD5Init, call MD5Update as
 * needed on buffers full of bytes, and then call MD5Final, which
 * will fill a supplied 16-byte array with the digest.
 *
 * Changed so as no longer to depend on Colin Plumb's `usual.h'
 * header definitions
 *  - Ian Jackson <ian@chiark.greenend.org.uk>.
 * Still in the public domain.
*/

typedef struct MD5Context {
    unsigned int buf[4];
    unsigned int bytes[2];
    unsigned int in[16];
}MD5Context;

void md5_init(MD5Context *context);
void md5_update(MD5Context *context, unsigned char const *buf, unsigned int len);
void md5_final(unsigned char digest[16], MD5Context *context);
void md5_transform(unsigned int buf[4], unsigned int const in[16]);

void print_md5(unsigned char digest[16]);
void write_md5(EbBufferHeaderType *recon_buffer, CLInput *cli, MD5Context *md5);
