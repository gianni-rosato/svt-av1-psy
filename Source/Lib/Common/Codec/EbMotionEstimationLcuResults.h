/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbMotionEstimationLcuResults_h
#define EbMotionEstimationLcuResults_h

#include "EbDefinitions.h"
#include "EbDefinitions.h"
#ifdef __cplusplus
extern "C" {
#endif

#define MAX_ME_PU_COUNT         209  // Sum of all the possible partitions which have both deminsions greater than 4.

#define MAX_SS_ME_PU_COUNT       (849 * 4 + 5) // Sum of all the possible partitions which have both deminsions greater or equal to 4.

    // i.e. no 4x4, 8x4, or 4x8 partitions
#define SQUARE_PU_COUNT          85
#define MAX_ME_CANDIDATE_PER_PU   3

    typedef struct MeCandidate 
    {

        union {
            struct {
                signed short     x_mv_l0;  //Note: Do not change the order of these fields
                signed short     y_mv_l0;
                signed short     x_mv_l1;
                signed short     y_mv_l1;
            };
            uint64_t mvs;
        };

        unsigned    distortion : 32;     // 20-bits holds maximum SAD of 64x64 PU

        unsigned    direction : 8;      // 0: uni-pred L0, 1: uni-pred L1, 2: bi-pred

    } MeCandidate;

    // move this to a new file with ctor & dtor
    typedef struct MeLcuResults 
    {
        uint32_t          lcu_distortion;
        uint8_t          *total_me_candidate_index;
        int16_t          x_mv_hme_search_center[MAX_NUM_OF_REF_PIC_LIST][MAX_REF_IDX];
        int16_t          y_mv_hme_search_center[MAX_NUM_OF_REF_PIC_LIST][MAX_REF_IDX];
        MeCandidate **me_candidate;
        MeCandidate  *me_candidate_array;

    } MeLcuResults;




    typedef struct  DistDir 
    {
        unsigned    distortion : 32; //20bits are enough
        unsigned    direction : 2;
    } DistDir;


    typedef struct MeCuResults 
    {
        union {
            struct {
                signed short     x_mv_l0;
                signed short     y_mv_l0;
                signed short     x_mv_l1;
                signed short     y_mv_l1;
            };
            uint64_t mvs;
        };

        DistDir    distortion_direction[3];

        uint8_t        total_me_candidate_index;
        uint8_t       me_nsq[2]; // 2 Number of reference lists
    } MeCuResults;

#ifdef __cplusplus
}
#endif
#endif // EbMotionEstimationLcuResults_h
