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
EbErrorType neighbor_array_unit_ctor32(
    NeighborArrayUnit32_t **na_unit_dbl_ptr,
    uint32_t   max_picture_width,
    uint32_t   max_picture_height,
    uint32_t   unit_size,
    uint32_t   granularity_normal,
    uint32_t   granularity_top_left,
    uint32_t   type_mask)
{
    NeighborArrayUnit32_t *na_unit_ptr;
    EB_MALLOC(NeighborArrayUnit32_t*, na_unit_ptr, sizeof(NeighborArrayUnit32_t), EB_N_PTR);

    *na_unit_dbl_ptr = na_unit_ptr;
    na_unit_ptr->unit_size = (uint8_t)(unit_size);
    na_unit_ptr->granularity_normal = (uint8_t)(granularity_normal);
    na_unit_ptr->granularityNormalLog2 = (uint8_t)(Log2f(na_unit_ptr->granularity_normal));
    na_unit_ptr->granularity_top_left = (uint8_t)(granularity_top_left);
    na_unit_ptr->granularityTopLeftLog2 = (uint8_t)(Log2f(na_unit_ptr->granularity_top_left));
    na_unit_ptr->leftArraySize = (uint16_t)((type_mask & NEIGHBOR_ARRAY_UNIT_LEFT_MASK) ? max_picture_height >> na_unit_ptr->granularityNormalLog2 : 0);
    na_unit_ptr->topArraySize = (uint16_t)((type_mask & NEIGHBOR_ARRAY_UNIT_TOP_MASK) ? max_picture_width >> na_unit_ptr->granularityNormalLog2 : 0);
    na_unit_ptr->topLeftArraySize = (uint16_t)((type_mask & NEIGHBOR_ARRAY_UNIT_TOPLEFT_MASK) ? (max_picture_width + max_picture_height) >> na_unit_ptr->granularityTopLeftLog2 : 0);

    if (na_unit_ptr->leftArraySize) {
        EB_MALLOC(uint32_t*, na_unit_ptr->leftArray, na_unit_ptr->unit_size * na_unit_ptr->leftArraySize, EB_N_PTR);
    }
    else {
        na_unit_ptr->leftArray = (uint32_t*)EB_NULL;
    }

    if (na_unit_ptr->topArraySize) {
        EB_MALLOC(uint32_t*, na_unit_ptr->topArray, na_unit_ptr->unit_size * na_unit_ptr->topArraySize, EB_N_PTR);
    }
    else {
        na_unit_ptr->topArray = (uint32_t*)EB_NULL;
    }

    if (na_unit_ptr->topLeftArraySize) {
        EB_MALLOC(uint32_t*, na_unit_ptr->topLeftArray, na_unit_ptr->unit_size * na_unit_ptr->topLeftArraySize, EB_N_PTR);
    }
    else {
        na_unit_ptr->topLeftArray = (uint32_t*)EB_NULL;
    }

    return EB_ErrorNone;
}

EbErrorType neighbor_array_unit_ctor(
    NeighborArrayUnit_t **na_unit_dbl_ptr,
    uint32_t   max_picture_width,
    uint32_t   max_picture_height,
    uint32_t   unit_size,
    uint32_t   granularity_normal,
    uint32_t   granularity_top_left,
    uint32_t   type_mask)
{
    NeighborArrayUnit_t *na_unit_ptr;
    EB_MALLOC(NeighborArrayUnit_t*, na_unit_ptr, sizeof(NeighborArrayUnit_t), EB_N_PTR);

    *na_unit_dbl_ptr = na_unit_ptr;
    na_unit_ptr->unit_size = (uint8_t)(unit_size);
    na_unit_ptr->granularity_normal = (uint8_t)(granularity_normal);
    na_unit_ptr->granularityNormalLog2 = (uint8_t)(Log2f(na_unit_ptr->granularity_normal));
    na_unit_ptr->granularity_top_left = (uint8_t)(granularity_top_left);
    na_unit_ptr->granularityTopLeftLog2 = (uint8_t)(Log2f(na_unit_ptr->granularity_top_left));
    na_unit_ptr->leftArraySize = (uint16_t)((type_mask & NEIGHBOR_ARRAY_UNIT_LEFT_MASK) ? max_picture_height >> na_unit_ptr->granularityNormalLog2 : 0);
    na_unit_ptr->topArraySize = (uint16_t)((type_mask & NEIGHBOR_ARRAY_UNIT_TOP_MASK) ? max_picture_width >> na_unit_ptr->granularityNormalLog2 : 0);
    na_unit_ptr->topLeftArraySize = (uint16_t)((type_mask & NEIGHBOR_ARRAY_UNIT_TOPLEFT_MASK) ? (max_picture_width + max_picture_height) >> na_unit_ptr->granularityTopLeftLog2 : 0);

    if (na_unit_ptr->leftArraySize) {
        EB_MALLOC(uint8_t*, na_unit_ptr->leftArray, na_unit_ptr->unit_size * na_unit_ptr->leftArraySize, EB_N_PTR);
    }
    else {
        na_unit_ptr->leftArray = (uint8_t*)EB_NULL;
    }

    if (na_unit_ptr->topArraySize) {
        EB_MALLOC(uint8_t*, na_unit_ptr->topArray, na_unit_ptr->unit_size * na_unit_ptr->topArraySize, EB_N_PTR);
    }
    else {
        na_unit_ptr->topArray = (uint8_t*)EB_NULL;
    }

    if (na_unit_ptr->topLeftArraySize) {
        EB_MALLOC(uint8_t*, na_unit_ptr->topLeftArray, na_unit_ptr->unit_size * na_unit_ptr->topLeftArraySize, EB_N_PTR);
    }
    else {
        na_unit_ptr->topLeftArray = (uint8_t*)EB_NULL;
    }

    return EB_ErrorNone;
}


/*************************************************
 * Neighbor Array Unit Reset
 *************************************************/

void neighbor_array_unit_reset32(NeighborArrayUnit32_t *na_unit_ptr)
{
    if (na_unit_ptr->leftArray) {
        EB_MEMSET(na_unit_ptr->leftArray, ~0, na_unit_ptr->unit_size * na_unit_ptr->leftArraySize);
    }

    if (na_unit_ptr->topArray) {
        EB_MEMSET(na_unit_ptr->topArray, ~0, na_unit_ptr->unit_size * na_unit_ptr->topArraySize);
    }

    if (na_unit_ptr->topLeftArray) {
        EB_MEMSET(na_unit_ptr->topLeftArray, ~0, na_unit_ptr->unit_size * na_unit_ptr->topLeftArraySize);
    }

    return;
}
void neighbor_array_unit_reset(NeighborArrayUnit_t *na_unit_ptr)
{
    if (na_unit_ptr->leftArray) {
        EB_MEMSET(na_unit_ptr->leftArray, ~0, na_unit_ptr->unit_size * na_unit_ptr->leftArraySize);
    }

    if (na_unit_ptr->topArray) {
        EB_MEMSET(na_unit_ptr->topArray, ~0, na_unit_ptr->unit_size * na_unit_ptr->topArraySize);
    }

    if (na_unit_ptr->topLeftArray) {
        EB_MEMSET(na_unit_ptr->topLeftArray, ~0, na_unit_ptr->unit_size * na_unit_ptr->topLeftArraySize);
    }

    return;
}


/*************************************************
 * Neighbor Array Unit Get Top Index
 *************************************************/
uint32_t GetNeighborArrayUnitTopLeftIndex32(
    NeighborArrayUnit32_t *na_unit_ptr,
    int32_t               loc_x,
    int32_t               loc_y)
{
    return na_unit_ptr->leftArraySize + (loc_x >> na_unit_ptr->granularityTopLeftLog2) - (loc_y >> na_unit_ptr->granularityTopLeftLog2);
}

uint32_t get_neighbor_array_unit_top_left_index(
    NeighborArrayUnit_t *na_unit_ptr,
    int32_t               loc_x,
    int32_t               loc_y)
{
    return na_unit_ptr->leftArraySize + (loc_x >> na_unit_ptr->granularityTopLeftLog2) - (loc_y >> na_unit_ptr->granularityTopLeftLog2);
}

void update_recon_neighbor_array(
    NeighborArrayUnit_t *na_unit_ptr,
    uint8_t               *src_ptr_top,
    uint8_t               *src_ptr_left,
    uint32_t               pic_origin_x,
    uint32_t               pic_origin_y,
    uint32_t               block_width,
    uint32_t               block_height)
{

    uint8_t  *dst_ptr;

    dst_ptr = na_unit_ptr->topArray +
        get_neighbor_array_unit_top_index(
            na_unit_ptr,
            pic_origin_x) * na_unit_ptr->unit_size;
    EB_MEMCPY(dst_ptr, src_ptr_top, block_width);


    dst_ptr = na_unit_ptr->leftArray +
        get_neighbor_array_unit_left_index(
            na_unit_ptr,
            pic_origin_y) * na_unit_ptr->unit_size;
    EB_MEMCPY(dst_ptr, src_ptr_left, block_height);


    //na_unit_ptr->topLeftArray[ (MAX_PICTURE_HEIGHT_SIZE>>is_chroma) + pic_origin_x - pic_origin_y] = srcPtr2[block_height-1];

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


    readPtr = src_ptr_top;//+ ((block_height - 1) * stride);

      // Copy bottom row
    dst_ptr =
        //    topLeftArray_chkn+
        na_unit_ptr->topLeftArray +
        get_neighbor_array_unit_top_left_index(
            na_unit_ptr,
            pic_origin_x,
            pic_origin_y + (block_height - 1)) * na_unit_ptr->unit_size;

    EB_MEMCPY(dst_ptr, readPtr, block_width);

    // Reset readPtr to the right-column
    readPtr = src_ptr_left;// + (block_width - 1);

    // Copy right column
    dst_ptr =
        //  topLeftArray_chkn+
        na_unit_ptr->topLeftArray +
        get_neighbor_array_unit_top_left_index(
            na_unit_ptr,
            pic_origin_x + (block_width - 1),
            pic_origin_y) * na_unit_ptr->unit_size;

    dstStep = -1;
    readStep = 1;//stride;
    count = block_height;

    for (idx = 0; idx < count; ++idx) {

        *dst_ptr = *readPtr;

        dst_ptr += dstStep;
        readPtr += readStep;
    }


    return;
}

/*************************************************
 * Neighbor Array Sample Update
 *************************************************/
void neighbor_array_unit_sample_write(
    NeighborArrayUnit_t *na_unit_ptr,
    uint8_t               *src_ptr,
    uint32_t               stride,
    uint32_t               src_origin_x,
    uint32_t               src_origin_y,
    uint32_t               pic_origin_x,
    uint32_t               pic_origin_y,
    uint32_t               block_width,
    uint32_t               block_height,
    uint32_t               neighbor_array_type_mask)
{
    uint32_t idx;
    uint8_t  *dst_ptr;
    uint8_t  *readPtr;

    int32_t dstStep;
    int32_t readStep;
    uint32_t count;

    // Adjust the Source ptr to start at the origin of the block being updated.
    src_ptr += ((src_origin_y * stride) + src_origin_x) * na_unit_ptr->unit_size;

    if (neighbor_array_type_mask & NEIGHBOR_ARRAY_UNIT_TOP_MASK) {

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
        readPtr = src_ptr + ((block_height - 1) * stride);

        dst_ptr = na_unit_ptr->topArray +
            get_neighbor_array_unit_top_index(
                na_unit_ptr,
                pic_origin_x) * na_unit_ptr->unit_size;

        dstStep = na_unit_ptr->unit_size;
        readStep = na_unit_ptr->unit_size;
        count = block_width;

        for (idx = 0; idx < count; ++idx) {

            *dst_ptr = *readPtr;

            dst_ptr += dstStep;
            readPtr += readStep;
        }

    }

    if (neighbor_array_type_mask & NEIGHBOR_ARRAY_UNIT_LEFT_MASK) {

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
        readPtr = src_ptr + (block_width - 1);

        dst_ptr = na_unit_ptr->leftArray +
            get_neighbor_array_unit_left_index(
                na_unit_ptr,
                pic_origin_y) * na_unit_ptr->unit_size;

        dstStep = 1;
        readStep = stride;
        count = block_height;

        for (idx = 0; idx < count; ++idx) {

            *dst_ptr = *readPtr;

            dst_ptr += dstStep;
            readPtr += readStep;
        }

    }

    if (neighbor_array_type_mask & NEIGHBOR_ARRAY_UNIT_TOPLEFT_MASK) {

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
        readPtr = src_ptr + ((block_height - 1) * stride);

        // Copy bottom row
        dst_ptr =
            na_unit_ptr->topLeftArray +
            get_neighbor_array_unit_top_left_index(
                na_unit_ptr,
                pic_origin_x,
                pic_origin_y + (block_height - 1)) * na_unit_ptr->unit_size;

        EB_MEMCPY(dst_ptr, readPtr, block_width);

        // Reset readPtr to the right-column
        readPtr = src_ptr + (block_width - 1);

        // Copy right column
        dst_ptr =
            na_unit_ptr->topLeftArray +
            get_neighbor_array_unit_top_left_index(
                na_unit_ptr,
                pic_origin_x + (block_width - 1),
                pic_origin_y) * na_unit_ptr->unit_size;

        dstStep = -1;
        readStep = stride;
        count = block_height;

        for (idx = 0; idx < count; ++idx) {

            *dst_ptr = *readPtr;

            dst_ptr += dstStep;
            readPtr += readStep;
        }

    }

    return;
}

/*************************************************
 * Neighbor Array Sample Update for 16 bit case
 *************************************************/
void neighbor_array_unit16bit_sample_write(
    NeighborArrayUnit_t *na_unit_ptr,
    uint16_t               *src_ptr,
    uint32_t               stride,
    uint32_t               src_origin_x,
    uint32_t               src_origin_y,
    uint32_t               pic_origin_x,
    uint32_t               pic_origin_y,
    uint32_t               block_width,
    uint32_t               block_height,
    uint32_t               neighbor_array_type_mask)
{
    uint32_t idx;
    uint16_t  *dst_ptr;
    uint16_t  *readPtr;

    int32_t dstStep;
    int32_t readStep;
    uint32_t count;

    // Adjust the Source ptr to start at the origin of the block being updated.
    src_ptr += ((src_origin_y * stride) + src_origin_x)/*CHKN  * na_unit_ptr->unit_size*/;

    if (neighbor_array_type_mask & NEIGHBOR_ARRAY_UNIT_TOP_MASK) {

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
        readPtr = src_ptr + ((block_height - 1) * stride);

        dst_ptr = (uint16_t*)(na_unit_ptr->topArray) +
            get_neighbor_array_unit_top_index(
                na_unit_ptr,
                pic_origin_x);//CHKN * na_unit_ptr->unit_size;

        dstStep = na_unit_ptr->unit_size;
        readStep = na_unit_ptr->unit_size;
        count = block_width;

        for (idx = 0; idx < count; ++idx) {

            *dst_ptr = *readPtr;

            dst_ptr += 1;
            readPtr += 1;
        }

    }

    if (neighbor_array_type_mask & NEIGHBOR_ARRAY_UNIT_LEFT_MASK) {

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
        readPtr = src_ptr + (block_width - 1);

        dst_ptr = (uint16_t*)(na_unit_ptr->leftArray) +
            get_neighbor_array_unit_left_index(
                na_unit_ptr,
                pic_origin_y);//CHKN * na_unit_ptr->unit_size;

        dstStep = 1;
        readStep = stride;
        count = block_height;

        for (idx = 0; idx < count; ++idx) {

            *dst_ptr = *readPtr;

            dst_ptr += dstStep;
            readPtr += readStep;
        }

    }

    if (neighbor_array_type_mask & NEIGHBOR_ARRAY_UNIT_TOPLEFT_MASK) {

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
        readPtr = src_ptr + ((block_height - 1) * stride);

        // Copy bottom row
        dst_ptr =
            (uint16_t*)(na_unit_ptr->topLeftArray) +
            get_neighbor_array_unit_top_left_index(
                na_unit_ptr,
                pic_origin_x,
                pic_origin_y + (block_height - 1));

        dstStep = 1;
        readStep = 1;
        count = block_width;

        for (idx = 0; idx < count; ++idx) {

            *dst_ptr = *readPtr;

            dst_ptr += dstStep;
            readPtr += readStep;
        }

        // Reset readPtr to the right-column
        readPtr = src_ptr + (block_width - 1);

        // Copy right column
        dst_ptr =
            (uint16_t*)(na_unit_ptr->topLeftArray) +
            get_neighbor_array_unit_top_left_index(
                na_unit_ptr,
                pic_origin_x + (block_width - 1),
                pic_origin_y);//CHKN  * na_unit_ptr->unit_size;

        dstStep = -1;
        readStep = stride;
        count = block_height;

        for (idx = 0; idx < count; ++idx) {

            *dst_ptr = *readPtr;

            dst_ptr += dstStep;
            readPtr += readStep;
        }

    }

    return;
}
/*************************************************
 * Neighbor Array Unit Mode Write
 *************************************************/
void neighbor_array_unit_mode_write32(
    NeighborArrayUnit32_t *na_unit_ptr,
    uint32_t               value,
    uint32_t               origin_x,
    uint32_t               origin_y,
    uint32_t               block_width,
    uint32_t               block_height,
    uint32_t               neighbor_array_type_mask)
{
    uint32_t idx;
    uint32_t  *dst_ptr;

    uint32_t count;
    uint32_t naOffset;
    uint32_t naUnitSize;

    naUnitSize = 1;//na_unit_ptr->unit_size;

    if (neighbor_array_type_mask & NEIGHBOR_ARRAY_UNIT_TOP_MASK) {

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

        naOffset = get_neighbor_array_unit_top_index32(
            na_unit_ptr,
            origin_x);

        dst_ptr = na_unit_ptr->topArray +
            naOffset * naUnitSize;

        count = block_width >> na_unit_ptr->granularityNormalLog2;

        for (idx = 0; idx < count; ++idx) {

            memset32bit(dst_ptr, value, naUnitSize);

            dst_ptr += naUnitSize;
        }
    }

    if (neighbor_array_type_mask & NEIGHBOR_ARRAY_UNIT_LEFT_MASK) {

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

        naOffset = get_neighbor_array_unit_left_index32(
            na_unit_ptr,
            origin_y);

        dst_ptr = na_unit_ptr->leftArray +
            naOffset * naUnitSize;

        count = block_height >> na_unit_ptr->granularityNormalLog2;

        for (idx = 0; idx < count; ++idx) {

            memset32bit(dst_ptr, value, naUnitSize);

            dst_ptr += naUnitSize;
        }
    }

    if (neighbor_array_type_mask & NEIGHBOR_ARRAY_UNIT_TOPLEFT_MASK) {

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
            na_unit_ptr,
            origin_x,
            origin_y + (block_height - 1));

        // Copy bottom-row + right-column
        // *Note - start from the bottom-left corner
        dst_ptr = na_unit_ptr->topLeftArray +
            naOffset * naUnitSize;

        count = ((block_width + block_height) >> na_unit_ptr->granularityTopLeftLog2) - 1;

        for (idx = 0; idx < count; ++idx) {

            memset32bit(dst_ptr, value, naUnitSize);

            dst_ptr += naUnitSize;
        }
    }

    return;
}

void neighbor_array_unit_mode_write(
    NeighborArrayUnit_t *na_unit_ptr,
    uint8_t               *value,
    uint32_t               origin_x,
    uint32_t               origin_y,
    uint32_t               block_width,
    uint32_t               block_height,
    uint32_t               neighbor_array_type_mask)
{
    uint32_t idx, j;
    uint8_t  *dst_ptr;

    uint32_t count;
    uint32_t naOffset;
    uint32_t naUnitSize;

    naUnitSize = na_unit_ptr->unit_size;

    if (neighbor_array_type_mask & NEIGHBOR_ARRAY_UNIT_TOP_MASK) {

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

        naOffset = get_neighbor_array_unit_top_index(
            na_unit_ptr,
            origin_x);

        dst_ptr = na_unit_ptr->topArray +
            naOffset * naUnitSize;

        count = block_width >> na_unit_ptr->granularityNormalLog2;

        for (idx = 0; idx < count; ++idx) {
            /* memcpy less that 10 bytes*/
            for (j = 0; j < naUnitSize; ++j)
                dst_ptr[j] = value[j];

            dst_ptr += naUnitSize;
        }
    }

    if (neighbor_array_type_mask & NEIGHBOR_ARRAY_UNIT_LEFT_MASK) {

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

        naOffset = get_neighbor_array_unit_left_index(
            na_unit_ptr,
            origin_y);

        dst_ptr = na_unit_ptr->leftArray +
            naOffset * naUnitSize;

        count = block_height >> na_unit_ptr->granularityNormalLog2;

        for (idx = 0; idx < count; ++idx) {
            /* memcpy less that 10 bytes*/
            for (j = 0; j < naUnitSize; ++j)
                dst_ptr[j] = value[j];

            dst_ptr += naUnitSize;
        }
    }

    if (neighbor_array_type_mask & NEIGHBOR_ARRAY_UNIT_TOPLEFT_MASK) {

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

        naOffset = get_neighbor_array_unit_top_left_index(
            na_unit_ptr,
            origin_x,
            origin_y + (block_height - 1));

        // Copy bottom-row + right-column
        // *Note - start from the bottom-left corner
        dst_ptr = na_unit_ptr->topLeftArray +
            naOffset * naUnitSize;

        count = ((block_width + block_height) >> na_unit_ptr->granularityTopLeftLog2) - 1;

        for (idx = 0; idx < count; ++idx) {
            /* memcpy less that 10 bytes*/
            for (j = 0; j < naUnitSize; ++j)
                dst_ptr[j] = value[j];

            dst_ptr += naUnitSize;
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
    uint32_t                 neighbor_array_type_mask)
{

    uint32_t idx;
    uint8_t  *dst_ptr, *src_ptr;

    uint32_t count;
    uint32_t naOffset;
    uint32_t naUnitSize;

    UNUSED(idx);
    naUnitSize = na_src->unit_size;

    if (neighbor_array_type_mask & NEIGHBOR_ARRAY_UNIT_TOP_MASK) {

        naOffset = get_neighbor_array_unit_top_index(na_src, origin_x);
        src_ptr = na_src->topArray + naOffset * naUnitSize;
        dst_ptr = na_dst->topArray + naOffset * naUnitSize;
        count = bw >> na_src->granularityNormalLog2;

        EB_MEMCPY(dst_ptr, src_ptr, naUnitSize*count);

    }

    if (neighbor_array_type_mask & NEIGHBOR_ARRAY_UNIT_LEFT_MASK) {

        naOffset = get_neighbor_array_unit_left_index(na_src, origin_y);
        src_ptr = na_src->leftArray + naOffset * naUnitSize;
        dst_ptr = na_dst->leftArray + naOffset * naUnitSize;
        count = bh >> na_src->granularityNormalLog2;

        EB_MEMCPY(dst_ptr, src_ptr, naUnitSize*count);
    }

    if (neighbor_array_type_mask & NEIGHBOR_ARRAY_UNIT_TOPLEFT_MASK) {

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

        naOffset = get_neighbor_array_unit_top_left_index(
            na_src,
            origin_x,
            origin_y + (bh - 1));

        // Copy bottom-row + right-column
        // *Note - start from the bottom-left corner
        src_ptr = na_src->topLeftArray + naOffset * naUnitSize;
        dst_ptr = na_dst->topLeftArray + naOffset * naUnitSize;

        count = ((bw + bh) >> na_src->granularityTopLeftLog2) - 1;

        EB_MEMCPY(dst_ptr, src_ptr, naUnitSize*count);
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
    uint32_t                 neighbor_array_type_mask)
{

    uint32_t idx;
    uint32_t  *dst_ptr, *src_ptr;

    uint32_t count;
    uint32_t naOffset;
    uint32_t naUnitSize;

    UNUSED(idx);

    naUnitSize = na_src->unit_size;

    if (neighbor_array_type_mask & NEIGHBOR_ARRAY_UNIT_TOP_MASK) {

        naOffset = get_neighbor_array_unit_top_index32(na_src, origin_x);
        src_ptr = na_src->topArray + naOffset;
        dst_ptr = na_dst->topArray + naOffset;
        count = bw >> na_src->granularityNormalLog2;

        EB_MEMCPY(dst_ptr, src_ptr, naUnitSize*count);

    }

    if (neighbor_array_type_mask & NEIGHBOR_ARRAY_UNIT_LEFT_MASK) {

        naOffset = get_neighbor_array_unit_left_index32(na_src, origin_y);
        src_ptr = na_src->leftArray + naOffset;
        dst_ptr = na_dst->leftArray + naOffset;
        count = bh >> na_src->granularityNormalLog2;

        EB_MEMCPY(dst_ptr, src_ptr, naUnitSize*count);
    }
    if (neighbor_array_type_mask & NEIGHBOR_ARRAY_UNIT_TOPLEFT_MASK) {



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
        src_ptr = na_src->topLeftArray + naOffset;
        dst_ptr = na_dst->topLeftArray + naOffset;

        count = ((bw + bh) >> na_src->granularityTopLeftLog2) - 1;

        EB_MEMCPY(dst_ptr, src_ptr, naUnitSize*count);
    }
    return;
}

/*************************************************
 * Neighbor Array Unit Mode Write
 *************************************************/
void neighbor_array_unit_mv_write(
    NeighborArrayUnit_t *na_unit_ptr,
    uint8_t               *value,
    uint32_t               origin_x,
    uint32_t               origin_y,
    uint32_t               block_size)
{
    uint32_t idx;
    uint8_t  *dst_ptr;
    uint8_t   *naUnittopArray;
    uint8_t   *naUnitleftArray;
    uint8_t   *naUnittopLeftArray;

    uint32_t count;
    uint32_t naOffset;
    uint32_t naUnitSize;

    naUnitSize = na_unit_ptr->unit_size;
    naUnittopArray = na_unit_ptr->topArray;
    naUnitleftArray = na_unit_ptr->leftArray;
    naUnittopLeftArray = na_unit_ptr->topLeftArray;


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

    dst_ptr = naUnittopArray +
        naOffset * naUnitSize;

    //dstStep = naUnitSize;
    count = block_size >> 2;

    for (idx = 0; idx < count; ++idx) {

        EB_MEMCPY(dst_ptr, value, naUnitSize);

        dst_ptr += naUnitSize;
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

    dst_ptr = naUnitleftArray +
        naOffset * naUnitSize;


    for (idx = 0; idx < count; ++idx) {

        EB_MEMCPY(dst_ptr, value, naUnitSize);

        dst_ptr += naUnitSize;
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

    naOffset = get_neighbor_array_unit_top_left_index(
        na_unit_ptr,
        origin_x,
        origin_y + (block_size - 1));

    // Copy bottom-row + right-column
    // *Note - start from the bottom-left corner
    dst_ptr = naUnittopLeftArray +
        naOffset * naUnitSize;

    count = ((block_size + block_size) >> 2) - 1;

    for (idx = 0; idx < count; ++idx) {

        EB_MEMCPY(dst_ptr, value, naUnitSize);

        dst_ptr += naUnitSize;
    }

    return;
}
