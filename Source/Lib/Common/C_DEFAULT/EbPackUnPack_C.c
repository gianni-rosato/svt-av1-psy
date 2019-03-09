/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbPackUnPack_C.h"
#include "EbDefinitions.h"

/************************************************
* pack 8 and 2 bit 2D data into 10 bit data
************************************************/
void eb_enc_msb_pack2_d(
    uint8_t     *in8_bit_buffer,
    uint32_t     in8_stride,
    uint8_t     *inn_bit_buffer,
    uint16_t    *out16_bit_buffer,
    uint32_t     inn_stride,
    uint32_t     out_stride,
    uint32_t     width,
    uint32_t     height)
{
    uint64_t   j, k;
    uint16_t   outPixel;
    uint8_t    nBitPixel;

    //SIMD hint: use _mm_unpacklo_epi8 +_mm_unpackhi_epi8 to do the concatenation

    for (j = 0; j < height; j++)
    {
        for (k = 0; k < width; k++)
        {
            outPixel = (in8_bit_buffer[k + j * in8_stride]) << 2;
            nBitPixel = (inn_bit_buffer[k + j * inn_stride] >> 6) & 3;

            out16_bit_buffer[k + j * out_stride] = outPixel | nBitPixel;
        }
    }
}

/************************************************
* pack 8 and 2 bit 2D data into 10 bit data
2bit data storage : 4 2bit-pixels in one byte
************************************************/
void compressed_packmsb(
    uint8_t     *in8_bit_buffer,
    uint32_t     in8_stride,
    uint8_t     *inn_bit_buffer,
    uint16_t    *out16_bit_buffer,
    uint32_t     inn_stride,
    uint32_t     out_stride,
    uint32_t     width,
    uint32_t     height)
{
    uint64_t   row, kIdx;
    uint16_t   outPixel;
    uint8_t    nBitPixel;
    uint8_t   four2bitPels;

    for (row = 0; row < height; row++)
    {
        for (kIdx = 0; kIdx < width / 4; kIdx++)
        {

            four2bitPels = inn_bit_buffer[kIdx + row * inn_stride];

            nBitPixel = (four2bitPels >> 6) & 3;

            outPixel = in8_bit_buffer[kIdx * 4 + 0 + row * in8_stride] << 2;
            out16_bit_buffer[kIdx * 4 + 0 + row * out_stride] = outPixel | nBitPixel;

            nBitPixel = (four2bitPels >> 4) & 3;
            outPixel = in8_bit_buffer[kIdx * 4 + 1 + row * in8_stride] << 2;
            out16_bit_buffer[kIdx * 4 + 1 + row * out_stride] = outPixel | nBitPixel;

            nBitPixel = (four2bitPels >> 2) & 3;
            outPixel = in8_bit_buffer[kIdx * 4 + 2 + row * in8_stride] << 2;
            out16_bit_buffer[kIdx * 4 + 2 + row * out_stride] = outPixel | nBitPixel;

            nBitPixel = (four2bitPels >> 0) & 3;
            outPixel = in8_bit_buffer[kIdx * 4 + 3 + row * in8_stride] << 2;
            out16_bit_buffer[kIdx * 4 + 3 + row * out_stride] = outPixel | nBitPixel;


        }
    }
}

/************************************************
* convert unpacked nbit (n=2) data to compressedPAcked
2bit data storage : 4 2bit-pixels in one byte
************************************************/
void c_pack_c(
    const uint8_t     *inn_bit_buffer,
    uint32_t     inn_stride,
    uint8_t     *in_compn_bit_buffer,
    uint32_t     out_stride,
    uint8_t    *local_cache,
    uint32_t     width,
    uint32_t     height)
{
    uint32_t row_index, colIndex;
    (void)local_cache;

    for (row_index = 0; row_index < height; row_index++)
    {
        for (colIndex = 0; colIndex < width; colIndex += 4)
        {
            uint32_t i = colIndex + row_index * inn_stride;

            uint8_t compressedUnpackedPixel = 0;
            compressedUnpackedPixel = compressedUnpackedPixel | ((inn_bit_buffer[i + 0] >> 0) & 0xC0);//1100.0000
            compressedUnpackedPixel = compressedUnpackedPixel | ((inn_bit_buffer[i + 1] >> 2) & 0x30);//0011.0000
            compressedUnpackedPixel = compressedUnpackedPixel | ((inn_bit_buffer[i + 2] >> 4) & 0x0C);//0000.1100
            compressedUnpackedPixel = compressedUnpackedPixel | ((inn_bit_buffer[i + 3] >> 6) & 0x03);//0000.0011

            uint32_t j = colIndex / 4 + row_index * out_stride;
            in_compn_bit_buffer[j] = compressedUnpackedPixel;
        }
    }

}


/************************************************
* unpack 10 bit data into  8 and 2 bit 2D data
************************************************/
void eb_enc_msb_un_pack2_d(
    uint16_t      *in16_bit_buffer,
    uint32_t       in_stride,
    uint8_t       *out8_bit_buffer,
    uint8_t       *outn_bit_buffer,
    uint32_t       out8_stride,
    uint32_t       outn_stride,
    uint32_t       width,
    uint32_t       height)
{
    uint64_t   j, k;
    uint16_t   inPixel;
    uint8_t    tmpPixel;
    for (j = 0; j < height; j++)
    {
        for (k = 0; k < width; k++)
        {
            inPixel = in16_bit_buffer[k + j * in_stride];
            out8_bit_buffer[k + j * out8_stride] = (uint8_t)(inPixel >> 2);
            tmpPixel = (uint8_t)(inPixel << 6);
            outn_bit_buffer[k + j * outn_stride] = tmpPixel;
        }
    }

}
void un_pack8_bit_data(
    uint16_t      *in16_bit_buffer,
    uint32_t       in_stride,
    uint8_t       *out8_bit_buffer,
    uint32_t       out8_stride,
    uint32_t       width,
    uint32_t       height)
{
    uint64_t   j, k;
    uint16_t   inPixel;
    //uint8_t    tmpPixel;
    for (j = 0; j < height; j++)
    {
        for (k = 0; k < width; k++)
        {
            inPixel = in16_bit_buffer[k + j * in_stride];
            out8_bit_buffer[k + j * out8_stride] = (uint8_t)(inPixel >> 2);
            //tmpPixel = (uint8_t)(inPixel << 6);
            //outn_bit_buffer[k + j*outn_stride] = tmpPixel;
        }
    }

}
void unpack_avg(
    uint16_t *ref16_l0,
    uint32_t  ref_l0_stride,
    uint16_t *ref16_l1,
    uint32_t  ref_l1_stride,
    uint8_t  *dst_ptr,
    uint32_t  dst_stride,
    uint32_t  width,
    uint32_t  height)
{

    uint64_t   j, k;
    uint8_t   inPixelL0, inPixelL1;

    for (j = 0; j < height; j++)
    {
        for (k = 0; k < width; k++)
        {
            inPixelL0 = (uint8_t)(ref16_l0[k + j * ref_l0_stride] >> 2);
            inPixelL1 = (uint8_t)(ref16_l1[k + j * ref_l1_stride] >> 2);
            dst_ptr[k + j * dst_stride] = (inPixelL0 + inPixelL1 + 1) >> 1;

        }
    }


}

void unpack_avg_safe_sub(
    uint16_t *ref16_l0,
    uint32_t  ref_l0_stride,
    uint16_t *ref16_l1,
    uint32_t  ref_l1_stride,
    uint8_t  *dst_ptr,
    uint32_t  dst_stride,
    EbBool      sub_pred,
    uint32_t  width,
    uint32_t  height)
{

    uint64_t   j, k;
    uint8_t   inPixelL0, inPixelL1;

    for (j = 0; j < height; j++)
    {
        for (k = 0; k < width; k++)
        {
            inPixelL0 = (uint8_t)(ref16_l0[k + j * ref_l0_stride] >> 2);
            inPixelL1 = (uint8_t)(ref16_l1[k + j * ref_l1_stride] >> 2);
            dst_ptr[k + j * dst_stride] = (inPixelL0 + inPixelL1 + 1) >> 1;

        }
    }

    if (sub_pred) {
        //Last row
        j = height * 2 - 1;
        for (k = 0; k < width; k++)
        {
            inPixelL0 = (uint8_t)(ref16_l0[k + j * ref_l0_stride / 2] >> 2);
            inPixelL1 = (uint8_t)(ref16_l1[k + j * ref_l1_stride / 2] >> 2);
            dst_ptr[k + j * dst_stride / 2] = (inPixelL0 + inPixelL1 + 1) >> 1;

        }

    }

}
