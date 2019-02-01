/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbPictureOperators_C.h"
#include "EbUtility.h"

/*********************************
* Picture Copy Kernel
*********************************/
void PictureCopyKernel(
    EbByte                  src,
    uint32_t                   src_stride,
    EbByte                  dst,
    uint32_t                   dst_stride,
    uint32_t                   areaWidth,
    uint32_t                   areaHeight,
    uint32_t                   bytesPerSample)  //=1 always)
{
    uint32_t sampleCount = 0;
    const uint32_t sampleTotalCount = areaWidth * areaHeight;
    const uint32_t copyLength = areaWidth * bytesPerSample;

    src_stride *= bytesPerSample;
    dst_stride *= bytesPerSample;

    while (sampleCount < sampleTotalCount) {
        EB_MEMCPY(dst, src, copyLength);
        src += src_stride;
        dst += dst_stride;
        sampleCount += areaWidth;
    }

    return;
}

// C equivalents

uint64_t SpatialFullDistortionKernel(
    uint8_t   *input,
    uint32_t   inputStride,
    uint8_t   *recon,
    uint32_t   reconStride,
    uint32_t   areaWidth,
    uint32_t   areaHeight)
{
    uint32_t  columnIndex;
    uint32_t  rowIndex = 0;

    uint64_t  spatialDistortion = 0;

    while (rowIndex < areaHeight) {

        columnIndex = 0;
        while (columnIndex < areaWidth) {
            spatialDistortion += (int64_t)SQR((int64_t)(input[columnIndex]) - (recon[columnIndex]));
            ++columnIndex;
        }

        input += inputStride;
        recon += reconStride;
        ++rowIndex;
    }
    return spatialDistortion;
}





