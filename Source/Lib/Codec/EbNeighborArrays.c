/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>
#include <string.h>

#include "EbNeighborArrays.h"
#include "EbUtility.h"
#include "EbPictureOperators.h"

#define UNUSED(x) (void)(x)
/*************************************************
 * Neighbor Array Unit Ctor
 *************************************************/
EbErrorType NeighborArrayUnitCtor32(
    NeighborArrayUnit32_t **naUnitDblPtr,
    uint32_t   maxPictureWidth,
    uint32_t   maxPictureHeight,
    uint32_t   unitSize,
    uint32_t   granularityNormal,
    uint32_t   granularityTopLeft,
    uint32_t   typeMask)
{
    NeighborArrayUnit32_t *naUnitPtr;
    EB_MALLOC(NeighborArrayUnit32_t*, naUnitPtr, sizeof(NeighborArrayUnit32_t), EB_N_PTR);

    *naUnitDblPtr = naUnitPtr;
    naUnitPtr->unitSize = (uint8_t)(unitSize);
    naUnitPtr->granularityNormal = (uint8_t)(granularityNormal);
    naUnitPtr->granularityNormalLog2 = (uint8_t)(Log2f(naUnitPtr->granularityNormal));
    naUnitPtr->granularityTopLeft = (uint8_t)(granularityTopLeft);
    naUnitPtr->granularityTopLeftLog2 = (uint8_t)(Log2f(naUnitPtr->granularityTopLeft));
    naUnitPtr->leftArraySize = (uint16_t)((typeMask & NEIGHBOR_ARRAY_UNIT_LEFT_MASK) ? maxPictureHeight >> naUnitPtr->granularityNormalLog2 : 0);
    naUnitPtr->topArraySize = (uint16_t)((typeMask & NEIGHBOR_ARRAY_UNIT_TOP_MASK) ? maxPictureWidth >> naUnitPtr->granularityNormalLog2 : 0);
    naUnitPtr->topLeftArraySize = (uint16_t)((typeMask & NEIGHBOR_ARRAY_UNIT_TOPLEFT_MASK) ? (maxPictureWidth + maxPictureHeight) >> naUnitPtr->granularityTopLeftLog2 : 0);

    if (naUnitPtr->leftArraySize) {
        EB_MALLOC(uint32_t*, naUnitPtr->leftArray, naUnitPtr->unitSize * naUnitPtr->leftArraySize, EB_N_PTR);
    }
    else {
        naUnitPtr->leftArray = (uint32_t*)EB_NULL;
    }

    if (naUnitPtr->topArraySize) {
        EB_MALLOC(uint32_t*, naUnitPtr->topArray, naUnitPtr->unitSize * naUnitPtr->topArraySize, EB_N_PTR);
    }
    else {
        naUnitPtr->topArray = (uint32_t*)EB_NULL;
    }

    if (naUnitPtr->topLeftArraySize) {
        EB_MALLOC(uint32_t*, naUnitPtr->topLeftArray, naUnitPtr->unitSize * naUnitPtr->topLeftArraySize, EB_N_PTR);
    }
    else {
        naUnitPtr->topLeftArray = (uint32_t*)EB_NULL;
    }

    return EB_ErrorNone;
}

EbErrorType NeighborArrayUnitCtor(
    NeighborArrayUnit_t **naUnitDblPtr,
    uint32_t   maxPictureWidth,
    uint32_t   maxPictureHeight,
    uint32_t   unitSize,
    uint32_t   granularityNormal,
    uint32_t   granularityTopLeft,
    uint32_t   typeMask)
{
    NeighborArrayUnit_t *naUnitPtr;
    EB_MALLOC(NeighborArrayUnit_t*, naUnitPtr, sizeof(NeighborArrayUnit_t), EB_N_PTR);

    *naUnitDblPtr = naUnitPtr;
    naUnitPtr->unitSize = (uint8_t)(unitSize);
    naUnitPtr->granularityNormal = (uint8_t)(granularityNormal);
    naUnitPtr->granularityNormalLog2 = (uint8_t)(Log2f(naUnitPtr->granularityNormal));
    naUnitPtr->granularityTopLeft = (uint8_t)(granularityTopLeft);
    naUnitPtr->granularityTopLeftLog2 = (uint8_t)(Log2f(naUnitPtr->granularityTopLeft));
    naUnitPtr->leftArraySize = (uint16_t)((typeMask & NEIGHBOR_ARRAY_UNIT_LEFT_MASK) ? maxPictureHeight >> naUnitPtr->granularityNormalLog2 : 0);
    naUnitPtr->topArraySize = (uint16_t)((typeMask & NEIGHBOR_ARRAY_UNIT_TOP_MASK) ? maxPictureWidth >> naUnitPtr->granularityNormalLog2 : 0);
    naUnitPtr->topLeftArraySize = (uint16_t)((typeMask & NEIGHBOR_ARRAY_UNIT_TOPLEFT_MASK) ? (maxPictureWidth + maxPictureHeight) >> naUnitPtr->granularityTopLeftLog2 : 0);

    if (naUnitPtr->leftArraySize) {
        EB_MALLOC(uint8_t*, naUnitPtr->leftArray, naUnitPtr->unitSize * naUnitPtr->leftArraySize, EB_N_PTR);
    }
    else {
        naUnitPtr->leftArray = (uint8_t*)EB_NULL;
    }

    if (naUnitPtr->topArraySize) {
        EB_MALLOC(uint8_t*, naUnitPtr->topArray, naUnitPtr->unitSize * naUnitPtr->topArraySize, EB_N_PTR);
    }
    else {
        naUnitPtr->topArray = (uint8_t*)EB_NULL;
    }

    if (naUnitPtr->topLeftArraySize) {
        EB_MALLOC(uint8_t*, naUnitPtr->topLeftArray, naUnitPtr->unitSize * naUnitPtr->topLeftArraySize, EB_N_PTR);
    }
    else {
        naUnitPtr->topLeftArray = (uint8_t*)EB_NULL;
    }

    return EB_ErrorNone;
}


/*************************************************
 * Neighbor Array Unit Reset
 *************************************************/

void NeighborArrayUnitReset32(NeighborArrayUnit32_t *naUnitPtr)
{
    if (naUnitPtr->leftArray) {
        EB_MEMSET(naUnitPtr->leftArray, ~0, naUnitPtr->unitSize * naUnitPtr->leftArraySize);
    }

    if (naUnitPtr->topArray) {
        EB_MEMSET(naUnitPtr->topArray, ~0, naUnitPtr->unitSize * naUnitPtr->topArraySize);
    }

    if (naUnitPtr->topLeftArray) {
        EB_MEMSET(naUnitPtr->topLeftArray, ~0, naUnitPtr->unitSize * naUnitPtr->topLeftArraySize);
    }

    return;
}
void NeighborArrayUnitReset(NeighborArrayUnit_t *naUnitPtr)
{
    if (naUnitPtr->leftArray) {
        EB_MEMSET(naUnitPtr->leftArray, ~0, naUnitPtr->unitSize * naUnitPtr->leftArraySize);
    }

    if (naUnitPtr->topArray) {
        EB_MEMSET(naUnitPtr->topArray, ~0, naUnitPtr->unitSize * naUnitPtr->topArraySize);
    }

    if (naUnitPtr->topLeftArray) {
        EB_MEMSET(naUnitPtr->topLeftArray, ~0, naUnitPtr->unitSize * naUnitPtr->topLeftArraySize);
    }

    return;
}


/*************************************************
 * Neighbor Array Unit Get Left Index
 *************************************************/
uint32_t GetNeighborArrayUnitLeftIndex32(
    NeighborArrayUnit32_t *naUnitPtr,
    uint32_t               locY)
{
    return (locY >> naUnitPtr->granularityNormalLog2);
}

uint32_t GetNeighborArrayUnitLeftIndex(
    NeighborArrayUnit_t *naUnitPtr,
    uint32_t               locY)
{
    return (locY >> naUnitPtr->granularityNormalLog2);
}

/*************************************************
 * Neighbor Array Unit Get Top Index
 *************************************************/
uint32_t GetNeighborArrayUnitTopIndex32(
    NeighborArrayUnit32_t *naUnitPtr,
    uint32_t               locX)
{
    return (locX >> naUnitPtr->granularityNormalLog2);
}

uint32_t GetNeighborArrayUnitTopIndex(
    NeighborArrayUnit_t *naUnitPtr,
    uint32_t               locX)
{
    return (locX >> naUnitPtr->granularityNormalLog2);
}

/*************************************************
 * Neighbor Array Unit Get Top Index
 *************************************************/
uint32_t GetNeighborArrayUnitTopLeftIndex32(
    NeighborArrayUnit32_t *naUnitPtr,
    int32_t               locX,
    int32_t               locY)
{
    return naUnitPtr->leftArraySize + (locX >> naUnitPtr->granularityTopLeftLog2) - (locY >> naUnitPtr->granularityTopLeftLog2);
}

uint32_t GetNeighborArrayUnitTopLeftIndex(
    NeighborArrayUnit_t *naUnitPtr,
    int32_t               locX,
    int32_t               locY)
{
    return naUnitPtr->leftArraySize + (locX >> naUnitPtr->granularityTopLeftLog2) - (locY >> naUnitPtr->granularityTopLeftLog2);
}

void update_recon_neighbor_array(
    NeighborArrayUnit_t *naUnitPtr,
    uint8_t               *srcPtrTop,
    uint8_t               *srcPtrLeft,
    uint32_t               picOriginX,
    uint32_t               picOriginY,
    uint32_t               blockWidth,
    uint32_t               blockHeight)
{

    uint8_t  *dstPtr;

    dstPtr = naUnitPtr->topArray +
        GetNeighborArrayUnitTopIndex(
            naUnitPtr,
            picOriginX) * naUnitPtr->unitSize;
#if OPT_MEMCPY
    EB_MEMCPY(dstPtr, srcPtrTop, blockWidth);
#else
    memcpy(dstPtr, srcPtrTop, blockWidth);
#endif


    dstPtr = naUnitPtr->leftArray +
        GetNeighborArrayUnitLeftIndex(
            naUnitPtr,
            picOriginY) * naUnitPtr->unitSize;
#if OPT_MEMCPY
    EB_MEMCPY(dstPtr, srcPtrLeft, blockHeight);
#else
    memcpy(dstPtr, srcPtrLeft, blockHeight);
#endif


    //naUnitPtr->topLeftArray[ (MAX_PICTURE_HEIGHT_SIZE>>is_chroma) + picOriginX - picOriginY] = srcPtr2[blockHeight-1];

     /*
        //   Top-left Neighbor Array
        //
        //    4-5--6--7--------------
        //    3 \      \
        //    2  \      \
        //    1   \      \
        //    |\   xxxxxx7
        //    | \  x     6
        //    |  \ x     5
        //    |   \1x2x3x4
        //    |
        //
        //  The top-left neighbor array is updated with the reversed samples
        //    from the right column and bottom row of the source block
        //
        // Index = origin_x - origin_y
        */

    uint32_t idx;

    uint8_t  *readPtr;

    int32_t dstStep;
    int32_t readStep;
    uint32_t count;


    readPtr = srcPtrTop;//+ ((blockHeight - 1) * stride);

      // Copy bottom row
    dstPtr =
        //    topLeftArray_chkn+
        naUnitPtr->topLeftArray +
        GetNeighborArrayUnitTopLeftIndex(
            naUnitPtr,
            picOriginX,
            picOriginY + (blockHeight - 1)) * naUnitPtr->unitSize;

    EB_MEMCPY(dstPtr, readPtr, blockWidth);

    // Reset readPtr to the right-column
    readPtr = srcPtrLeft;// + (blockWidth - 1);

    // Copy right column
    dstPtr =
        //  topLeftArray_chkn+
        naUnitPtr->topLeftArray +
        GetNeighborArrayUnitTopLeftIndex(
            naUnitPtr,
            picOriginX + (blockWidth - 1),
            picOriginY) * naUnitPtr->unitSize;

    dstStep = -1;
    readStep = 1;//stride;
    count = blockHeight;

    for (idx = 0; idx < count; ++idx) {

        *dstPtr = *readPtr;

        dstPtr += dstStep;
        readPtr += readStep;
    }


    return;
}

/*************************************************
 * Neighbor Array Sample Update
 *************************************************/
void NeighborArrayUnitSampleWrite(
    NeighborArrayUnit_t *naUnitPtr,
    uint8_t               *srcPtr,
    uint32_t               stride,
    uint32_t               srcOriginX,
    uint32_t               srcOriginY,
    uint32_t               picOriginX,
    uint32_t               picOriginY,
    uint32_t               blockWidth,
    uint32_t               blockHeight,
    uint32_t               neighborArrayTypeMask)
{
    uint32_t idx;
    uint8_t  *dstPtr;
    uint8_t  *readPtr;

    int32_t dstStep;
    int32_t readStep;
    uint32_t count;

    // Adjust the Source ptr to start at the origin of the block being updated.
    srcPtr += ((srcOriginY * stride) + srcOriginX) * naUnitPtr->unitSize;

    if (neighborArrayTypeMask & NEIGHBOR_ARRAY_UNIT_TOP_MASK) {

        //
        //     ----------12345678---------------------  Top Neighbor Array
        //                ^    ^
        //                |    |
        //                |    |
        //               xxxxxxxx
        //               x      x
        //               x      x
        //               12345678
        //
        //  The top neighbor array is updated with the samples from the
        //    bottom row of the source block
        //
        //  Index = origin_x

        // Adjust readPtr to the bottom-row
        readPtr = srcPtr + ((blockHeight - 1) * stride);

        dstPtr = naUnitPtr->topArray +
            GetNeighborArrayUnitTopIndex(
                naUnitPtr,
                picOriginX) * naUnitPtr->unitSize;

        dstStep = naUnitPtr->unitSize;
        readStep = naUnitPtr->unitSize;
        count = blockWidth;

        for (idx = 0; idx < count; ++idx) {

            *dstPtr = *readPtr;

            dstPtr += dstStep;
            readPtr += readStep;
        }

    }

    if (neighborArrayTypeMask & NEIGHBOR_ARRAY_UNIT_LEFT_MASK) {

        //   Left Neighbor Array
        //
        //    |
        //    |
        //    1         xxxxxxx1
        //    2  <----  x      2
        //    3  <----  x      3
        //    4         xxxxxxx4
        //    |
        //    |
        //
        //  The left neighbor array is updated with the samples from the
        //    right column of the source block
        //
        //  Index = origin_y

        // Adjust readPtr to the right-column
        readPtr = srcPtr + (blockWidth - 1);

        dstPtr = naUnitPtr->leftArray +
            GetNeighborArrayUnitLeftIndex(
                naUnitPtr,
                picOriginY) * naUnitPtr->unitSize;

        dstStep = 1;
        readStep = stride;
        count = blockHeight;

        for (idx = 0; idx < count; ++idx) {

            *dstPtr = *readPtr;

            dstPtr += dstStep;
            readPtr += readStep;
        }

    }

    if (neighborArrayTypeMask & NEIGHBOR_ARRAY_UNIT_TOPLEFT_MASK) {

        /*
        //   Top-left Neighbor Array
        //
        //    4-5--6--7--------------
        //    3 \      \
        //    2  \      \
        //    1   \      \
        //    |\   xxxxxx7
        //    | \  x     6
        //    |  \ x     5
        //    |   \1x2x3x4
        //    |
        //
        //  The top-left neighbor array is updated with the reversed samples
        //    from the right column and bottom row of the source block
        //
        // Index = origin_x - origin_y
        */

        // Adjust readPtr to the bottom-row
        readPtr = srcPtr + ((blockHeight - 1) * stride);

        // Copy bottom row
        dstPtr =
            naUnitPtr->topLeftArray +
            GetNeighborArrayUnitTopLeftIndex(
                naUnitPtr,
                picOriginX,
                picOriginY + (blockHeight - 1)) * naUnitPtr->unitSize;

        EB_MEMCPY(dstPtr, readPtr, blockWidth);

        // Reset readPtr to the right-column
        readPtr = srcPtr + (blockWidth - 1);

        // Copy right column
        dstPtr =
            naUnitPtr->topLeftArray +
            GetNeighborArrayUnitTopLeftIndex(
                naUnitPtr,
                picOriginX + (blockWidth - 1),
                picOriginY) * naUnitPtr->unitSize;

        dstStep = -1;
        readStep = stride;
        count = blockHeight;

        for (idx = 0; idx < count; ++idx) {

            *dstPtr = *readPtr;

            dstPtr += dstStep;
            readPtr += readStep;
        }

    }

    return;
}

/*************************************************
 * Neighbor Array Sample Update for 16 bit case
 *************************************************/
void NeighborArrayUnit16bitSampleWrite(
    NeighborArrayUnit_t *naUnitPtr,
    uint16_t               *srcPtr,
    uint32_t               stride,
    uint32_t               srcOriginX,
    uint32_t               srcOriginY,
    uint32_t               picOriginX,
    uint32_t               picOriginY,
    uint32_t               blockWidth,
    uint32_t               blockHeight,
    uint32_t               neighborArrayTypeMask)
{
    uint32_t idx;
    uint16_t  *dstPtr;
    uint16_t  *readPtr;

    int32_t dstStep;
    int32_t readStep;
    uint32_t count;

    // Adjust the Source ptr to start at the origin of the block being updated.
    srcPtr += ((srcOriginY * stride) + srcOriginX)/*CHKN  * naUnitPtr->unitSize*/;

    if (neighborArrayTypeMask & NEIGHBOR_ARRAY_UNIT_TOP_MASK) {

        //
        //     ----------12345678---------------------  Top Neighbor Array
        //                ^    ^
        //                |    |
        //                |    |
        //               xxxxxxxx
        //               x      x
        //               x      x
        //               12345678
        //
        //  The top neighbor array is updated with the samples from the
        //    bottom row of the source block
        //
        //  Index = origin_x

        // Adjust readPtr to the bottom-row
        readPtr = srcPtr + ((blockHeight - 1) * stride);

        dstPtr = (uint16_t*)(naUnitPtr->topArray) +
            GetNeighborArrayUnitTopIndex(
                naUnitPtr,
                picOriginX);//CHKN * naUnitPtr->unitSize;

        dstStep = naUnitPtr->unitSize;
        readStep = naUnitPtr->unitSize;
        count = blockWidth;

        for (idx = 0; idx < count; ++idx) {

            *dstPtr = *readPtr;

            dstPtr += 1;
            readPtr += 1;
        }

    }

    if (neighborArrayTypeMask & NEIGHBOR_ARRAY_UNIT_LEFT_MASK) {

        //   Left Neighbor Array
        //
        //    |
        //    |
        //    1         xxxxxxx1
        //    2  <----  x      2
        //    3  <----  x      3
        //    4         xxxxxxx4
        //    |
        //    |
        //
        //  The left neighbor array is updated with the samples from the
        //    right column of the source block
        //
        //  Index = origin_y

        // Adjust readPtr to the right-column
        readPtr = srcPtr + (blockWidth - 1);

        dstPtr = (uint16_t*)(naUnitPtr->leftArray) +
            GetNeighborArrayUnitLeftIndex(
                naUnitPtr,
                picOriginY);//CHKN * naUnitPtr->unitSize;

        dstStep = 1;
        readStep = stride;
        count = blockHeight;

        for (idx = 0; idx < count; ++idx) {

            *dstPtr = *readPtr;

            dstPtr += dstStep;
            readPtr += readStep;
        }

    }

    if (neighborArrayTypeMask & NEIGHBOR_ARRAY_UNIT_TOPLEFT_MASK) {

        /*
        //   Top-left Neighbor Array
        //
        //    4-5--6--7--------------
        //    3 \      \
        //    2  \      \
        //    1   \      \
        //    |\   xxxxxx7
        //    | \  x     6
        //    |  \ x     5
        //    |   \1x2x3x4
        //    |
        //
        //  The top-left neighbor array is updated with the reversed samples
        //    from the right column and bottom row of the source block
        //
        // Index = origin_x - origin_y
        */

        // Adjust readPtr to the bottom-row
        readPtr = srcPtr + ((blockHeight - 1) * stride);

        // Copy bottom row
        dstPtr =
            (uint16_t*)(naUnitPtr->topLeftArray) +
            GetNeighborArrayUnitTopLeftIndex(
                naUnitPtr,
                picOriginX,
                picOriginY + (blockHeight - 1));

        dstStep = 1;
        readStep = 1;
        count = blockWidth;

        for (idx = 0; idx < count; ++idx) {

            *dstPtr = *readPtr;

            dstPtr += dstStep;
            readPtr += readStep;
        }

        // Reset readPtr to the right-column
        readPtr = srcPtr + (blockWidth - 1);

        // Copy right column
        dstPtr =
            (uint16_t*)(naUnitPtr->topLeftArray) +
            GetNeighborArrayUnitTopLeftIndex(
                naUnitPtr,
                picOriginX + (blockWidth - 1),
                picOriginY);//CHKN  * naUnitPtr->unitSize;

        dstStep = -1;
        readStep = stride;
        count = blockHeight;

        for (idx = 0; idx < count; ++idx) {

            *dstPtr = *readPtr;

            dstPtr += dstStep;
            readPtr += readStep;
        }

    }

    return;
}
/*************************************************
 * Neighbor Array Unit Mode Write
 *************************************************/
void NeighborArrayUnitModeWrite32(
    NeighborArrayUnit32_t *naUnitPtr,
    uint32_t               value,
    uint32_t               origin_x,
    uint32_t               origin_y,
    uint32_t               blockWidth,
    uint32_t               blockHeight,
    uint32_t               neighborArrayTypeMask)
{
    uint32_t idx;
    uint32_t  *dstPtr;

    uint32_t count;
    uint32_t naOffset;
    uint32_t naUnitSize;

    naUnitSize = 1;//naUnitPtr->unitSize;

    if (neighborArrayTypeMask & NEIGHBOR_ARRAY_UNIT_TOP_MASK) {

        //
        //     ----------12345678---------------------  Top Neighbor Array
        //                ^    ^
        //                |    |
        //                |    |
        //               xxxxxxxx
        //               x      x
        //               x      x
        //               12345678
        //
        //  The top neighbor array is updated with the samples from the
        //    bottom row of the source block
        //
        //  Index = origin_x

        naOffset = GetNeighborArrayUnitTopIndex32(
            naUnitPtr,
            origin_x);

        dstPtr = naUnitPtr->topArray +
            naOffset * naUnitSize;

        count = blockWidth >> naUnitPtr->granularityNormalLog2;

        for (idx = 0; idx < count; ++idx) {

            memset32bit(dstPtr, value, naUnitSize);

            dstPtr += naUnitSize;
        }
    }

    if (neighborArrayTypeMask & NEIGHBOR_ARRAY_UNIT_LEFT_MASK) {

        //   Left Neighbor Array
        //
        //    |
        //    |
        //    1         xxxxxxx1
        //    2  <----  x      2
        //    3  <----  x      3
        //    4         xxxxxxx4
        //    |
        //    |
        //
        //  The left neighbor array is updated with the samples from the
        //    right column of the source block
        //
        //  Index = origin_y

        naOffset = GetNeighborArrayUnitLeftIndex32(
            naUnitPtr,
            origin_y);

        dstPtr = naUnitPtr->leftArray +
            naOffset * naUnitSize;

        count = blockHeight >> naUnitPtr->granularityNormalLog2;

        for (idx = 0; idx < count; ++idx) {

            memset32bit(dstPtr, value, naUnitSize);

            dstPtr += naUnitSize;
        }
    }

    if (neighborArrayTypeMask & NEIGHBOR_ARRAY_UNIT_TOPLEFT_MASK) {

        /*
        //   Top-left Neighbor Array
        //
        //    4-5--6--7------------
        //    3 \      \
        //    2  \      \
        //    1   \      \
        //    |\   xxxxxx7
        //    | \  x     6
        //    |  \ x     5
        //    |   \1x2x3x4
        //    |
        //
        //  The top-left neighbor array is updated with the reversed samples
        //    from the right column and bottom row of the source block
        //
        // Index = origin_x - origin_y
        */

        naOffset = GetNeighborArrayUnitTopLeftIndex32(
            naUnitPtr,
            origin_x,
            origin_y + (blockHeight - 1));

        // Copy bottom-row + right-column
        // *Note - start from the bottom-left corner
        dstPtr = naUnitPtr->topLeftArray +
            naOffset * naUnitSize;

        count = ((blockWidth + blockHeight) >> naUnitPtr->granularityTopLeftLog2) - 1;

        for (idx = 0; idx < count; ++idx) {

            memset32bit(dstPtr, value, naUnitSize);

            dstPtr += naUnitSize;
        }
    }

    return;
}

void NeighborArrayUnitModeWrite(
    NeighborArrayUnit_t *naUnitPtr,
    uint8_t               *value,
    uint32_t               origin_x,
    uint32_t               origin_y,
    uint32_t               blockWidth,
    uint32_t               blockHeight,
    uint32_t               neighborArrayTypeMask)
{
    uint32_t idx;
    uint8_t  *dstPtr;

    uint32_t count;
    uint32_t naOffset;
    uint32_t naUnitSize;

    naUnitSize = naUnitPtr->unitSize;

    if (neighborArrayTypeMask & NEIGHBOR_ARRAY_UNIT_TOP_MASK) {

        //
        //     ----------12345678---------------------  Top Neighbor Array
        //                ^    ^
        //                |    |
        //                |    |
        //               xxxxxxxx
        //               x      x
        //               x      x
        //               12345678
        //
        //  The top neighbor array is updated with the samples from the
        //    bottom row of the source block
        //
        //  Index = origin_x

        naOffset = GetNeighborArrayUnitTopIndex(
            naUnitPtr,
            origin_x);

        dstPtr = naUnitPtr->topArray +
            naOffset * naUnitSize;

        count = blockWidth >> naUnitPtr->granularityNormalLog2;

        for (idx = 0; idx < count; ++idx) {

            EB_MEMCPY(dstPtr, value, naUnitSize);

            dstPtr += naUnitSize;
        }
    }

    if (neighborArrayTypeMask & NEIGHBOR_ARRAY_UNIT_LEFT_MASK) {

        //   Left Neighbor Array
        //
        //    |
        //    |
        //    1         xxxxxxx1
        //    2  <----  x      2
        //    3  <----  x      3
        //    4         xxxxxxx4
        //    |
        //    |
        //
        //  The left neighbor array is updated with the samples from the
        //    right column of the source block
        //
        //  Index = origin_y

        naOffset = GetNeighborArrayUnitLeftIndex(
            naUnitPtr,
            origin_y);

        dstPtr = naUnitPtr->leftArray +
            naOffset * naUnitSize;

        count = blockHeight >> naUnitPtr->granularityNormalLog2;

        for (idx = 0; idx < count; ++idx) {

            EB_MEMCPY(dstPtr, value, naUnitSize);

            dstPtr += naUnitSize;
        }
    }

    if (neighborArrayTypeMask & NEIGHBOR_ARRAY_UNIT_TOPLEFT_MASK) {

        /*
        //   Top-left Neighbor Array
        //
        //    4-5--6--7------------
        //    3 \      \
        //    2  \      \
        //    1   \      \
        //    |\   xxxxxx7
        //    | \  x     6
        //    |  \ x     5
        //    |   \1x2x3x4
        //    |
        //
        //  The top-left neighbor array is updated with the reversed samples
        //    from the right column and bottom row of the source block
        //
        // Index = origin_x - origin_y
        */

        naOffset = GetNeighborArrayUnitTopLeftIndex(
            naUnitPtr,
            origin_x,
            origin_y + (blockHeight - 1));

        // Copy bottom-row + right-column
        // *Note - start from the bottom-left corner
        dstPtr = naUnitPtr->topLeftArray +
            naOffset * naUnitSize;

        count = ((blockWidth + blockHeight) >> naUnitPtr->granularityTopLeftLog2) - 1;

        for (idx = 0; idx < count; ++idx) {

            EB_MEMCPY(dstPtr, value, naUnitSize);

            dstPtr += naUnitSize;
        }
    }

    return;
}

void copy_neigh_arr(
    NeighborArrayUnit_t   *na_src,
    NeighborArrayUnit_t   *na_dst,
    uint32_t               origin_x,
    uint32_t               origin_y,
    uint32_t               bw,
    uint32_t               bh,
    uint32_t                 neighborArrayTypeMask)
{

    uint32_t idx;
    uint8_t  *dstPtr, *srcPtr;

    uint32_t count;
    uint32_t naOffset;
    uint32_t naUnitSize;

    UNUSED(idx);
    naUnitSize = na_src->unitSize;

    if (neighborArrayTypeMask & NEIGHBOR_ARRAY_UNIT_TOP_MASK) {

        naOffset = GetNeighborArrayUnitTopIndex(na_src, origin_x);
        srcPtr = na_src->topArray + naOffset * naUnitSize;
        dstPtr = na_dst->topArray + naOffset * naUnitSize;
        count = bw >> na_src->granularityNormalLog2;

        EB_MEMCPY(dstPtr, srcPtr, naUnitSize*count);

    }

    if (neighborArrayTypeMask & NEIGHBOR_ARRAY_UNIT_LEFT_MASK) {

        naOffset = GetNeighborArrayUnitLeftIndex(na_src, origin_y);
        srcPtr = na_src->leftArray + naOffset * naUnitSize;
        dstPtr = na_dst->leftArray + naOffset * naUnitSize;
        count = bh >> na_src->granularityNormalLog2;

        EB_MEMCPY(dstPtr, srcPtr, naUnitSize*count);
    }

    if (neighborArrayTypeMask & NEIGHBOR_ARRAY_UNIT_TOPLEFT_MASK) {

        /*
        //   Top-left Neighbor Array
        //
        //    4-5--6--7------------
        //    3 \      \
        //    2  \      \
        //    1   \      \
        //    |\   xxxxxx7
        //    | \  x     6
        //    |  \ x     5
        //    |   \1x2x3x4
        //    |
        //
        //  The top-left neighbor array is updated with the reversed samples
        //    from the right column and bottom row of the source block
        //
        // Index = origin_x - origin_y
        */

        naOffset = GetNeighborArrayUnitTopLeftIndex(
            na_src,
            origin_x,
            origin_y + (bh - 1));

        // Copy bottom-row + right-column
        // *Note - start from the bottom-left corner
        srcPtr = na_src->topLeftArray + naOffset * naUnitSize;
        dstPtr = na_dst->topLeftArray + naOffset * naUnitSize;

        count = ((bw + bh) >> na_src->granularityTopLeftLog2) - 1;

        EB_MEMCPY(dstPtr, srcPtr, naUnitSize*count);
    }

    return;
}

void copy_neigh_arr_32(
    NeighborArrayUnit32_t   *na_src,
    NeighborArrayUnit32_t   *na_dst,
    uint32_t               origin_x,
    uint32_t               origin_y,
    uint32_t               bw,
    uint32_t               bh,
    uint32_t                 neighborArrayTypeMask)
{

    uint32_t idx;
    uint32_t  *dstPtr, *srcPtr;

    uint32_t count;
    uint32_t naOffset;
    uint32_t naUnitSize;

    UNUSED(idx);

    naUnitSize = na_src->unitSize;

    if (neighborArrayTypeMask & NEIGHBOR_ARRAY_UNIT_TOP_MASK) {

        naOffset = GetNeighborArrayUnitTopIndex32(na_src, origin_x);
        srcPtr = na_src->topArray + naOffset;
        dstPtr = na_dst->topArray + naOffset;
        count = bw >> na_src->granularityNormalLog2;

        EB_MEMCPY(dstPtr, srcPtr, naUnitSize*count);

    }

    if (neighborArrayTypeMask & NEIGHBOR_ARRAY_UNIT_LEFT_MASK) {

        naOffset = GetNeighborArrayUnitLeftIndex32(na_src, origin_y);
        srcPtr = na_src->leftArray + naOffset;
        dstPtr = na_dst->leftArray + naOffset;
        count = bh >> na_src->granularityNormalLog2;

        EB_MEMCPY(dstPtr, srcPtr, naUnitSize*count);
    }
    if (neighborArrayTypeMask & NEIGHBOR_ARRAY_UNIT_TOPLEFT_MASK) {



        /*
        //   Top-left Neighbor Array
        //
        //    4-5--6--7------------
        //    3 \      \
        //    2  \      \
        //    1   \      \
        //    |\   xxxxxx7
        //    | \  x     6
        //    |  \ x     5
        //    |   \1x2x3x4
        //    |
        //
        //  The top-left neighbor array is updated with the reversed samples
        //    from the right column and bottom row of the source block
        //
        // Index = origin_x - origin_y
        */

        naOffset = GetNeighborArrayUnitTopLeftIndex32(
            na_src,
            origin_x,
            origin_y + (bh - 1));

        // Copy bottom-row + right-column
        // *Note - start from the bottom-left corner
        srcPtr = na_src->topLeftArray + naOffset;
        dstPtr = na_dst->topLeftArray + naOffset;

        count = ((bw + bh) >> na_src->granularityTopLeftLog2) - 1;

        EB_MEMCPY(dstPtr, srcPtr, naUnitSize*count);
    }
    return;
}

/*************************************************
 * Neighbor Array Unit Mode Write
 *************************************************/
void NeighborArrayUnitMvWrite(
    NeighborArrayUnit_t *naUnitPtr,
    uint8_t               *value,
    uint32_t               origin_x,
    uint32_t               origin_y,
    uint32_t               blockSize)
{
    uint32_t idx;
    uint8_t  *dstPtr;
    uint8_t   *naUnittopArray;
    uint8_t   *naUnitleftArray;
    uint8_t   *naUnittopLeftArray;

    uint32_t count;
    uint32_t naOffset;
    uint32_t naUnitSize;

    naUnitSize = naUnitPtr->unitSize;
    naUnittopArray = naUnitPtr->topArray;
    naUnitleftArray = naUnitPtr->leftArray;
    naUnittopLeftArray = naUnitPtr->topLeftArray;


    //
    //     ----------12345678---------------------  Top Neighbor Array
    //                ^    ^
    //                |    |
    //                |    |
    //               xxxxxxxx
    //               x      x
    //               x      x
    //               12345678
    //
    //  The top neighbor array is updated with the samples from the
    //    bottom row of the source block
    //
    //  Index = origin_x

    naOffset = origin_x >> 2;

    dstPtr = naUnittopArray +
        naOffset * naUnitSize;

    //dstStep = naUnitSize;
    count = blockSize >> 2;

    for (idx = 0; idx < count; ++idx) {

        EB_MEMCPY(dstPtr, value, naUnitSize);

        dstPtr += naUnitSize;
    }


    //   Left Neighbor Array
    //
    //    |
    //    |
    //    1         xxxxxxx1
    //    2  <----  x      2
    //    3  <----  x      3
    //    4         xxxxxxx4
    //    |
    //    |
    //
    //  The left neighbor array is updated with the samples from the
    //    right column of the source block
    //
    //  Index = origin_y

    naOffset = origin_y >> 2;

    dstPtr = naUnitleftArray +
        naOffset * naUnitSize;


    for (idx = 0; idx < count; ++idx) {

        EB_MEMCPY(dstPtr, value, naUnitSize);

        dstPtr += naUnitSize;
    }


    /*
    //   Top-left Neighbor Array
    //
    //    4-5--6--7------------
    //    3 \      \
    //    2  \      \
    //    1   \      \
    //    |\   xxxxxx7
    //    | \  x     6
    //    |  \ x     5
    //    |   \1x2x3x4
    //    |
    //
    //  The top-left neighbor array is updated with the reversed samples
    //    from the right column and bottom row of the source block
    //
    // Index = origin_x - origin_y
    */

    naOffset = GetNeighborArrayUnitTopLeftIndex(
        naUnitPtr,
        origin_x,
        origin_y + (blockSize - 1));

    // Copy bottom-row + right-column
    // *Note - start from the bottom-left corner
    dstPtr = naUnittopLeftArray +
        naOffset * naUnitSize;

    count = ((blockSize + blockSize) >> 2) - 1;

    for (idx = 0; idx < count; ++idx) {

        EB_MEMCPY(dstPtr, value, naUnitSize);

        dstPtr += naUnitSize;
    }

    return;
}
