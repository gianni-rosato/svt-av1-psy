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

    typedef struct MeCandidate_s {

        union {
            struct {
                signed short     xMvL0;  //Note: Do not change the order of these fields
                signed short     yMvL0;
                signed short     xMvL1;
                signed short     yMvL1;
            };
            uint64_t MVs;
        };

        unsigned    distortion : 32;     // 20-bits holds maximum SAD of 64x64 PU

        unsigned    direction : 8;      // 0: uni-pred L0, 1: uni-pred L1, 2: bi-pred

    } MeCandidate_t;

    // move this to a new file with ctor & dtor
    typedef struct MeLcuResults_s {

        uint32_t          lcuDistortion;
        uint8_t          *totalMeCandidateIndex;
        int16_t          xMvHmeSearchCenter[MAX_NUM_OF_REF_PIC_LIST][MAX_REF_IDX];
        int16_t          yMvHmeSearchCenter[MAX_NUM_OF_REF_PIC_LIST][MAX_REF_IDX];
        MeCandidate_t **me_candidate;
        MeCandidate_t  *meCandidateArray;

    } MeLcuResults_t;




    typedef struct  DistDir_s {
        unsigned    distortion : 32; //20bits are enough
        unsigned    direction : 2;
    } DistDir_t;


    typedef struct MeCuResults_s {
        union {
            struct {
                signed short     xMvL0;
                signed short     yMvL0;
                signed short     xMvL1;
                signed short     yMvL1;
            };
            uint64_t MVs;
        };

        DistDir_t    distortionDirection[3];

        uint8_t        totalMeCandidateIndex;
#if NSQ_OPTIMASATION
        uint8_t       me_nsq[2]; // 2 Number of reference lists
#endif
    } MeCuResults_t;

#ifdef __cplusplus
}
#endif
#endif // EbMotionEstimationLcuResults_h
