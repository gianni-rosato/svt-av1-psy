/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbMotionEstimationLcuResults_h
#define EbMotionEstimationLcuResults_h

#include "EbDefinitions.h"
#include "EbObject.h"
#ifdef __cplusplus
extern "C" {
#endif

#define MAX_ME_PU_COUNT \
    209 // Sum of all the possible partitions which have both deminsions greater than 4.

#define ME_RES_CAND_MRP_MODE_0 23 // [Single Ref = 7] + [BiDir = 12 = 3*4 ] + [UniDir = 4 = 3+1]
#define ME_MV_MRP_MODE_0 7 // [7 = 4+3]

#define ME_RES_CAND_MRP_MODE_1 10 // [Single Ref = 4] + [UniDir = 4 = 2*2] + [UniDir = 2 = 1+1]
#define ME_MV_MRP_MODE_1 4 // [4 = 2+2]

#define MAX_SS_ME_PU_COUNT \
    (849 * 4 +             \
     5) // Sum of all the possible partitions which have both deminsions greater or equal to 4.

// i.e. no 4x4, 8x4, or 4x8 partitions
#define SQUARE_PU_COUNT 85
typedef struct MeCandidate {
    uint8_t direction : 2;
    uint8_t ref_idx_l0 : 2; // allows for up to 4 references
    uint8_t ref_idx_l1 : 2;
    uint8_t ref0_list : 1;
    uint8_t ref1_list : 1;
} MeCandidate;
typedef struct MvCandidate {
    signed short x_mv;
    signed short y_mv;
} MvCandidate;
// move this to a new file with ctor & dtor
typedef struct MeSbResults {
    EbDctor       dctor;
    uint32_t      sb_distortion;
    uint8_t *     total_me_candidate_index;
#if ME_MEM_OPT
    MvCandidate *me_mv_array;
    MeCandidate *me_candidate_array;
#else
    MeCandidate **me_candidate;
    MeCandidate * me_candidate_array;
    MvCandidate **me_mv_array;
#endif
    // [PU][LAST, LAST2, LAST3, GOLD, BWD, ALT2, ALT] if MRP Mode 0,
    // [PU][LAST, LAST2, BWD, ALT2] if MRP Mode 1,
#if !NSQ_REMOVAL_CODE_CLEAN_UP
    uint32_t max_number_of_pus_per_sb;
#endif
#if INTER_COMP_REDESIGN
    uint8_t do_comp[MAX_NUM_OF_REF_PIC_LIST][REF_LIST_MAX_DEPTH];
#endif
} MeSbResults;
#ifdef __cplusplus
}
#endif
#endif // EbMotionEstimationLcuResults_h
