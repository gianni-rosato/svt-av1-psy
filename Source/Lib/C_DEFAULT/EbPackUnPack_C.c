/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbPackUnPack_C.h"
#include "EbDefinitions.h"

/************************************************
* pack 8 and 2 bit 2D data into 10 bit data
************************************************/
void EB_ENC_msbPack2D(
    uint8_t     *in8BitBuffer,
    uint32_t     in8Stride,
    uint8_t     *innBitBuffer,
    uint16_t    *out16BitBuffer,
    uint32_t     innStride,
    uint32_t     outStride,
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
            outPixel = (in8BitBuffer[k + j * in8Stride]) << 2;
            nBitPixel = (innBitBuffer[k + j * innStride] >> 6) & 3;

            out16BitBuffer[k + j * outStride] = outPixel | nBitPixel;
        }
    }
}

/************************************************
* pack 8 and 2 bit 2D data into 10 bit data
2bit data storage : 4 2bit-pixels in one byte
************************************************/
void CompressedPackmsb(
    uint8_t     *in8BitBuffer,
    uint32_t     in8Stride,
    uint8_t     *innBitBuffer,
    uint16_t    *out16BitBuffer,
    uint32_t     innStride,
    uint32_t     outStride,
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

            four2bitPels = innBitBuffer[kIdx + row * innStride];

            nBitPixel = (four2bitPels >> 6) & 3;

            outPixel = in8BitBuffer[kIdx * 4 + 0 + row * in8Stride] << 2;
            out16BitBuffer[kIdx * 4 + 0 + row * outStride] = outPixel | nBitPixel;

            nBitPixel = (four2bitPels >> 4) & 3;
            outPixel = in8BitBuffer[kIdx * 4 + 1 + row * in8Stride] << 2;
            out16BitBuffer[kIdx * 4 + 1 + row * outStride] = outPixel | nBitPixel;

            nBitPixel = (four2bitPels >> 2) & 3;
            outPixel = in8BitBuffer[kIdx * 4 + 2 + row * in8Stride] << 2;
            out16BitBuffer[kIdx * 4 + 2 + row * outStride] = outPixel | nBitPixel;

            nBitPixel = (four2bitPels >> 0) & 3;
            outPixel = in8BitBuffer[kIdx * 4 + 3 + row * in8Stride] << 2;
            out16BitBuffer[kIdx * 4 + 3 + row * outStride] = outPixel | nBitPixel;


        }
    }
}

/************************************************
* convert unpacked nbit (n=2) data to compressedPAcked
2bit data storage : 4 2bit-pixels in one byte
************************************************/
void CPack_C(
    const uint8_t     *innBitBuffer,
    uint32_t     innStride,
    uint8_t     *inCompnBitBuffer,
    uint32_t     outStride,
    uint8_t    *localCache,
    uint32_t     width,
    uint32_t     height)
{
    uint32_t rowIndex, colIndex;
    (void)localCache;

    for (rowIndex = 0; rowIndex < height; rowIndex++)
    {
        for (colIndex = 0; colIndex < width; colIndex += 4)
        {
            uint32_t i = colIndex + rowIndex * innStride;

            uint8_t compressedUnpackedPixel = 0;
            compressedUnpackedPixel = compressedUnpackedPixel | ((innBitBuffer[i + 0] >> 0) & 0xC0);//1100.0000
            compressedUnpackedPixel = compressedUnpackedPixel | ((innBitBuffer[i + 1] >> 2) & 0x30);//0011.0000
            compressedUnpackedPixel = compressedUnpackedPixel | ((innBitBuffer[i + 2] >> 4) & 0x0C);//0000.1100
            compressedUnpackedPixel = compressedUnpackedPixel | ((innBitBuffer[i + 3] >> 6) & 0x03);//0000.0011

            uint32_t j = colIndex / 4 + rowIndex * outStride;
            inCompnBitBuffer[j] = compressedUnpackedPixel;
        }
    }

}


/************************************************
* unpack 10 bit data into  8 and 2 bit 2D data
************************************************/
void EB_ENC_msbUnPack2D(
    uint16_t      *in16BitBuffer,
    uint32_t       inStride,
    uint8_t       *out8BitBuffer,
    uint8_t       *outnBitBuffer,
    uint32_t       out8Stride,
    uint32_t       outnStride,
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
            inPixel = in16BitBuffer[k + j * inStride];
            out8BitBuffer[k + j * out8Stride] = (uint8_t)(inPixel >> 2);
            tmpPixel = (uint8_t)(inPixel << 6);
            outnBitBuffer[k + j * outnStride] = tmpPixel;
        }
    }

}
void UnPack8BitData(
    uint16_t      *in16BitBuffer,
    uint32_t       inStride,
    uint8_t       *out8BitBuffer,
    uint32_t       out8Stride,
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
            inPixel = in16BitBuffer[k + j * inStride];
            out8BitBuffer[k + j * out8Stride] = (uint8_t)(inPixel >> 2);
            //tmpPixel = (uint8_t)(inPixel << 6);
            //outnBitBuffer[k + j*outnStride] = tmpPixel;
        }
    }

}
void UnpackAvg(
    uint16_t *ref16L0,
    uint32_t  refL0Stride,
    uint16_t *ref16L1,
    uint32_t  refL1Stride,
    uint8_t  *dstPtr,
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
            inPixelL0 = (uint8_t)(ref16L0[k + j * refL0Stride] >> 2);
            inPixelL1 = (uint8_t)(ref16L1[k + j * refL1Stride] >> 2);
            dstPtr[k + j * dst_stride] = (inPixelL0 + inPixelL1 + 1) >> 1;

        }
    }


}

void UnpackAvgSafeSub(
    uint16_t *ref16L0,
    uint32_t  refL0Stride,
    uint16_t *ref16L1,
    uint32_t  refL1Stride,
    uint8_t  *dstPtr,
    uint32_t  dst_stride,
    EbBool      subPred,
    uint32_t  width,
    uint32_t  height)
{

    uint64_t   j, k;
    uint8_t   inPixelL0, inPixelL1;

    for (j = 0; j < height; j++)
    {
        for (k = 0; k < width; k++)
        {
            inPixelL0 = (uint8_t)(ref16L0[k + j * refL0Stride] >> 2);
            inPixelL1 = (uint8_t)(ref16L1[k + j * refL1Stride] >> 2);
            dstPtr[k + j * dst_stride] = (inPixelL0 + inPixelL1 + 1) >> 1;

        }
    }

    if (subPred) {
        //Last row
        j = height * 2 - 1;
        for (k = 0; k < width; k++)
        {
            inPixelL0 = (uint8_t)(ref16L0[k + j * refL0Stride / 2] >> 2);
            inPixelL1 = (uint8_t)(ref16L1[k + j * refL1Stride / 2] >> 2);
            dstPtr[k + j * dst_stride / 2] = (inPixelL0 + inPixelL1 + 1) >> 1;

        }

    }

}
