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
        uint8_t    unitSize;
        uint8_t    granularityNormal;
        uint8_t    granularityNormalLog2;
        uint8_t    granularityTopLeft;
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
        uint8_t    unitSize;
        uint8_t    granularityNormal;
        uint8_t    granularityNormalLog2;
        uint8_t    granularityTopLeft;
        uint8_t    granularityTopLeftLog2;

    } NeighborArrayUnit32_t;

    extern EbErrorType NeighborArrayUnitCtor32(
        NeighborArrayUnit32_t **naUnitDblPtr,
        uint32_t   maxPictureWidth,
        uint32_t   maxPictureHeight,
        uint32_t   unitSize,
        uint32_t   granularityNormal,
        uint32_t   granularityTopLeft,
        uint32_t   typeMask);


    extern EbErrorType NeighborArrayUnitCtor(
        NeighborArrayUnit_t **naUnitDblPtr,
        uint32_t   maxPictureWidth,
        uint32_t   maxPictureHeight,
        uint32_t   unitSize,
        uint32_t   granularityNormal,
        uint32_t   granularityTopLeft,
        uint32_t   typeMask);

    extern void NeighborArrayUnitDtor(NeighborArrayUnit_t  *naUnitPtr);

    extern void NeighborArrayUnitReset(NeighborArrayUnit_t *naUnitPtr);

    extern void NeighborArrayUnitReset32(NeighborArrayUnit32_t *naUnitPtr);

    extern uint32_t GetNeighborArrayUnitLeftIndex32(
        NeighborArrayUnit32_t *naUnitPtr,
        uint32_t               locY);

    extern uint32_t GetNeighborArrayUnitTopIndex32(
        NeighborArrayUnit32_t *naUnitPtr,
        uint32_t               locX);

    extern uint32_t GetNeighborArrayUnitLeftIndex(
        NeighborArrayUnit_t *naUnitPtr,
        uint32_t               locY);

    extern uint32_t GetNeighborArrayUnitTopIndex(
        NeighborArrayUnit_t *naUnitPtr,
        uint32_t               locX);

    extern uint32_t GetNeighborArrayUnitTopLeftIndex(
        NeighborArrayUnit_t *naUnitPtr,
        int32_t               locX,
        int32_t               locY);

    extern void NeighborArrayUnitSampleWrite(
        NeighborArrayUnit_t *naUnitPtr,
        uint8_t               *srcPtr,
        uint32_t               stride,
        uint32_t               srcOriginX,
        uint32_t               srcOriginY,
        uint32_t               picOriginX,
        uint32_t               picOriginY,
        uint32_t               blockWidth,
        uint32_t               blockHeight,
        uint32_t               neighborArrayTypeMask);

    void update_recon_neighbor_array(
        NeighborArrayUnit_t *naUnitPtr,
        uint8_t               *srcPtrTop,
        uint8_t               *srcPtrLeft,
        uint32_t               picOriginX,
        uint32_t               picOriginY,
        uint32_t               blockWidth,
        uint32_t               blockHeight);


    void copy_neigh_arr(
        NeighborArrayUnit_t   *na_src,
        NeighborArrayUnit_t   *na_dst,
        uint32_t               origin_x,
        uint32_t               origin_y,
        uint32_t               bw,
        uint32_t               bh,
        uint32_t                 neighborArrayTypeMask);
    void copy_neigh_arr_32(
        NeighborArrayUnit32_t   *na_src,
        NeighborArrayUnit32_t   *na_dst,
        uint32_t               origin_x,
        uint32_t               origin_y,
        uint32_t               bw,
        uint32_t               bh,
        uint32_t                 neighborArrayTypeMask);


    extern void NeighborArrayUnit16bitSampleWrite(
        NeighborArrayUnit_t *naUnitPtr,
        uint16_t               *srcPtr,
        uint32_t               stride,
        uint32_t               srcOriginX,
        uint32_t               srcOriginY,
        uint32_t               picOriginX,
        uint32_t               picOriginY,
        uint32_t               blockWidth,
        uint32_t               blockHeight,
        uint32_t               neighborArrayTypeMask);

    extern void NeighborArrayUnitModeWrite32(
        NeighborArrayUnit32_t *naUnitPtr,
        uint32_t               value,
        uint32_t               origin_x,
        uint32_t               origin_y,
        uint32_t               blockWidth,
        uint32_t               blockHeight,
        uint32_t               neighborArrayTypeMask);

    extern void NeighborArrayUnitModeWrite(
        NeighborArrayUnit_t *naUnitPtr,
        uint8_t               *value,
        uint32_t               origin_x,
        uint32_t               origin_y,
        uint32_t               blockWidth,
        uint32_t               blockHeight,
        uint32_t               neighborArrayTypeMask);

    extern void NeighborArrayUnitMvWrite(
        NeighborArrayUnit_t *naUnitPtr,
        uint8_t               *value,
        uint32_t               origin_x,
        uint32_t               origin_y,
        uint32_t               blockSize);
#ifdef __cplusplus
}
#endif
#endif //EbNeighborArrays_h
