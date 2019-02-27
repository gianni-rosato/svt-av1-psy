/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbNeighborArrays_h
#define EbNeighborArrays_h

#include "EbDefinitions.h"
#include "EbSyntaxElements.h"
#include "EbMotionVectorUnit.h"
#ifdef __cplusplus
extern "C" {
#endif
    // Neighbor Array Granulairity
#define LCU_NEIGHBOR_ARRAY_GRANULARITY                  64
#define CU_NEIGHBOR_ARRAY_GRANULARITY                   8
#define PU_NEIGHBOR_ARRAY_GRANULARITY                   4
#define TU_NEIGHBOR_ARRAY_GRANULARITY                   4
#define SAMPLE_NEIGHBOR_ARRAY_GRANULARITY               1

    typedef enum NEIGHBOR_ARRAY_TYPE
    {
        NEIGHBOR_ARRAY_LEFT = 0,
        NEIGHBOR_ARRAY_TOP = 1,
        NEIGHBOR_ARRAY_TOPLEFT = 2,
        NEIGHBOR_ARRAY_INVALID = ~0
    } NEIGHBOR_ARRAY_TYPE;

#define NEIGHBOR_ARRAY_UNIT_LEFT_MASK                   (1 << NEIGHBOR_ARRAY_LEFT)
#define NEIGHBOR_ARRAY_UNIT_TOP_MASK                    (1 << NEIGHBOR_ARRAY_TOP)
#define NEIGHBOR_ARRAY_UNIT_TOPLEFT_MASK                (1 << NEIGHBOR_ARRAY_TOPLEFT)

#define NEIGHBOR_ARRAY_UNIT_FULL_MASK                   (NEIGHBOR_ARRAY_UNIT_LEFT_MASK | NEIGHBOR_ARRAY_UNIT_TOP_MASK | NEIGHBOR_ARRAY_UNIT_TOPLEFT_MASK)
#define NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK      (NEIGHBOR_ARRAY_UNIT_LEFT_MASK | NEIGHBOR_ARRAY_UNIT_TOP_MASK)

    typedef struct NeighborArrayUnit_s
    {
        uint8_t   *leftArray;
        uint8_t   *topArray;
        uint8_t   *topLeftArray;
        uint16_t   leftArraySize;
        uint16_t   topArraySize;
        uint16_t   topLeftArraySize;
        uint8_t    unit_size;
        uint8_t    granularity_normal;
        uint8_t    granularityNormalLog2;
        uint8_t    granularity_top_left;
        uint8_t    granularityTopLeftLog2;

    } NeighborArrayUnit_t;

    typedef struct NeighborArrayUnit32_s
    {
        uint32_t   *leftArray;
        uint32_t   *topArray;
        uint32_t   *topLeftArray;
        uint16_t   leftArraySize;
        uint16_t   topArraySize;
        uint16_t   topLeftArraySize;
        uint8_t    unit_size;
        uint8_t    granularity_normal;
        uint8_t    granularityNormalLog2;
        uint8_t    granularity_top_left;
        uint8_t    granularityTopLeftLog2;

    } NeighborArrayUnit32_t;

    extern EbErrorType neighbor_array_unit_ctor32(
        NeighborArrayUnit32_t **na_unit_dbl_ptr,
        uint32_t   max_picture_width,
        uint32_t   max_picture_height,
        uint32_t   unit_size,
        uint32_t   granularity_normal,
        uint32_t   granularity_top_left,
        uint32_t   type_mask);


    extern EbErrorType neighbor_array_unit_ctor(
        NeighborArrayUnit_t **na_unit_dbl_ptr,
        uint32_t   max_picture_width,
        uint32_t   max_picture_height,
        uint32_t   unit_size,
        uint32_t   granularity_normal,
        uint32_t   granularity_top_left,
        uint32_t   type_mask);

    extern void neighbor_array_unit_reset(NeighborArrayUnit_t *na_unit_ptr);

    extern void neighbor_array_unit_reset32(NeighborArrayUnit32_t *na_unit_ptr);

    extern uint32_t get_neighbor_array_unit_left_index32(
        NeighborArrayUnit32_t *na_unit_ptr,
        uint32_t               loc_y);

    extern uint32_t get_neighbor_array_unit_top_index32(
        NeighborArrayUnit32_t *na_unit_ptr,
        uint32_t               loc_x);

    extern uint32_t get_neighbor_array_unit_left_index(
        NeighborArrayUnit_t *na_unit_ptr,
        uint32_t             loc_y);

    extern uint32_t get_neighbor_array_unit_top_index(
        NeighborArrayUnit_t *na_unit_ptr,
        uint32_t             loc_x);

    extern uint32_t get_neighbor_array_unit_top_left_index(
        NeighborArrayUnit_t *na_unit_ptr,
        int32_t              loc_x,
        int32_t              loc_y);

    extern void neighbor_array_unit_sample_write(
        NeighborArrayUnit_t *na_unit_ptr,
        uint8_t             *src_ptr,
        uint32_t             stride,
        uint32_t             src_origin_x,
        uint32_t             src_origin_y,
        uint32_t             pic_origin_x,
        uint32_t             pic_origin_y,
        uint32_t             block_width,
        uint32_t             block_height,
        uint32_t             neighbor_array_type_mask);

    void update_recon_neighbor_array(
        NeighborArrayUnit_t *na_unit_ptr,
        uint8_t             *src_ptr_top,
        uint8_t             *src_ptr_left,
        uint32_t             pic_origin_x,
        uint32_t             pic_origin_y,
        uint32_t             block_width,
        uint32_t             block_height);


    void copy_neigh_arr(
        NeighborArrayUnit_t *na_src,
        NeighborArrayUnit_t *na_dst,
        uint32_t             origin_x,
        uint32_t             origin_y,
        uint32_t             bw,
        uint32_t             bh,
        uint32_t             neighbor_array_type_mask);

    void copy_neigh_arr_32(
        NeighborArrayUnit32_t *na_src,
        NeighborArrayUnit32_t *na_dst,
        uint32_t               origin_x,
        uint32_t               origin_y,
        uint32_t               bw,
        uint32_t               bh,
        uint32_t               neighbor_array_type_mask);

    extern void neighbor_array_unit16bit_sample_write(
        NeighborArrayUnit_t *na_unit_ptr,
        uint16_t            *src_ptr,
        uint32_t             stride,
        uint32_t             src_origin_x,
        uint32_t             src_origin_y,
        uint32_t             pic_origin_x,
        uint32_t             pic_origin_y,
        uint32_t             block_width,
        uint32_t             block_height,
        uint32_t             neighbor_array_type_mask);

    extern void neighbor_array_unit_mode_write32(
        NeighborArrayUnit32_t *na_unit_ptr,
        uint32_t               value,
        uint32_t               origin_x,
        uint32_t               origin_y,
        uint32_t               block_width,
        uint32_t               block_height,
        uint32_t               neighbor_array_type_mask);

    extern void neighbor_array_unit_mode_write(
        NeighborArrayUnit_t *na_unit_ptr,
        uint8_t             *value,
        uint32_t             origin_x,
        uint32_t             origin_y,
        uint32_t             block_width,
        uint32_t             block_height,
        uint32_t             neighbor_array_type_mask);

    extern void neighbor_array_unit_mv_write(
        NeighborArrayUnit_t *na_unit_ptr,
        uint8_t             *value,
        uint32_t             origin_x,
        uint32_t             origin_y,
        uint32_t             block_size);
#ifdef __cplusplus
}
#endif
#endif //EbNeighborArrays_h
