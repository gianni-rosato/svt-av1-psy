/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbModeDecisionConfiguration_h
#define EbModeDecisionConfiguration_h

#include "EbModeDecisionConfigurationProcess.h"
#include "EbSequenceControlSet.h"
#include "EbPictureControlSet.h"
#ifdef __cplusplus
extern "C" {
#endif

/*******************************************
 * EarlyModeDecisionLcu 
 *   predicts candidates (LCU)
 *******************************************/
extern EbErrorType EarlyModeDecisionLcu(
    SequenceControlSet                   *sequence_control_set_ptr,
    PictureControlSet_t                    *picture_control_set_ptr,
    LargestCodingUnit_t                    *sb_ptr,
    uint32_t                                  sb_index,
    ModeDecisionConfigurationContext_t     *context_ptr);


/*******************************************
* DeriveDeltaQPForEachLeafLcu
*   Derive Lcu For Each Leaf (LCU)
*******************************************/
extern EbErrorType DeriveDeltaQPForEachLeafLcu(
    SequenceControlSet                   *sequence_control_set_ptr,
    PictureControlSet_t                    *picture_control_set_ptr,
    LargestCodingUnit_t                    *sb_ptr,
    uint32_t                                  sb_index,
    int32_t                                  intra_min_distance,
    int32_t                                  intra_max_distance,
    int32_t                                  inter_min_distance,
    int32_t                                  inter_max_distance,
    ModeDecisionConfigurationContext_t     *context_ptr);

void QpmDeriveDeltaQpMapWeights(
    ModeDecisionConfigurationContext_t    *context_ptr,
    PictureControlSet_t                  *picture_control_set_ptr);

extern uint8_t DeriveContouringClass(
    PictureParentControlSet_t   *parentPcsPtr,
    uint16_t                       sb_index,
    uint8_t                        leaf_index);  
/**************************************
* Function Ptrs Definitions
**************************************/
typedef EbErrorType(*EB_MDC_FUNC)(
    MdcpLocalCodingUnit_t                   *localCuArray,
    uint32_t                                   cu_index,
    uint32_t                                   depth,
    EbBool                                 *mdcPrediction64);

#define Pred        0x01
#define Predp1      0x02
#define Predp2      0x04
#define Predp3      0x08
#define Predm1      0x10
#define Predm2      0x20
#define Predm3      0x40
#define ALL64       0x0F
#define ALL32       0x17
#define ALL16       0x33
#define ALL8        0x71
#define AllD        0x80

EB_ALIGN(16) static const uint8_t ndp_level_0[4] = {Pred + Predp1 + Predp2, Pred + Predp1, Pred + Predp1, Pred + Predm1};
EB_ALIGN(16) static const uint8_t ndp_level_1[4] = {Pred + Predp1         , Pred + Predp1, Pred + Predp1, Pred + Predm1 };

#ifdef __cplusplus
}
#endif
#endif // EbModeDecisionConfiguration_h