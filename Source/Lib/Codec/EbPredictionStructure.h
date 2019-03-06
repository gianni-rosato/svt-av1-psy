/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbPredictionStructure_h
#define EbPredictionStructure_h

#include "EbDefinitions.h"
#include "EbDefinitions.h"
#ifdef __cplusplus
extern "C" {
#endif
    /************************************************
     * Defines
     ************************************************/
#define FLAT_PREDICTION_STRUCTURE_PERIOD                                    1
#define TWO_LEVEL_HIERARCHICAL_PREDICTION_STRUCTURE_PERIOD                  2
#define THREE_LEVEL_HIERARCHICAL_PREDICTION_STRUCTURE_PERIOD                4
#define FOUR_LEVEL_HIERARCHICAL_PREDICTION_STRUCTURE_PERIOD                 8
#define MAX_PREDICTION_STRUCTURE_PERIOD                                     64

     /************************************************
      * RPS defines
      ************************************************/
#define RPS_UNUSED                                                          ~0
#define MAX_NUM_OF_NEGATIVE_REF_PICS                                        5
#define MAX_NUM_OF_POSITIVE_REF_PICS                                        5
#define MAX_NUM_OF_REF_PICS_TOTAL                                           (MAX_NUM_OF_POSITIVE_REF_PICS + MAX_NUM_OF_NEGATIVE_REF_PICS)

      /************************************************
       * Reference List
       *
       *   referenceList - Contains the deltaPOCs of
       *    the pictures referenced by the current
       *    picture.
       ************************************************/
    typedef struct ReferenceList_s
    {
        int32_t                              referenceList;
        uint32_t                              referenceListCount;

    } ReferenceList_t;

    /************************************************
     * Dependent List
     *
     *   listCount - Contains count of how
     *     deep into list should be used
     *     depending on how many references are
     *     being used in the prediction structure.
     *
     *   list - Contains the deltaPOCs of
     *     pictures that reference the current picture.
     *     The dependent list pictures must be grouped
     *     by the referenceCount group in ascending
     *     order.  The grouping is not display order!
     ************************************************/
    typedef struct DependentList_s
    {
        int32_t                             *list;
        uint32_t                              listCount;

    } DependentList_t;

    /************************************************
     * Prediction Structure Config Entry
     *   Contains the basic reference lists and
     *   configurations for each Prediction Structure
     *   Config Entry.
     ************************************************/
    typedef struct PredictionStructureConfigEntry_s {
        uint32_t                              temporal_layer_index;
        uint32_t                              decode_order;
        int32_t                              refList0;
        int32_t                              refList1;
    } PredictionStructureConfigEntry_t;

    /************************************************
     * Prediction Structure Config
     *   Contains a collection of basic control data
     *   for the basic prediction structure.
     ************************************************/
    typedef struct PredictionStructureConfig_s {
        uint32_t                              entryCount;
        PredictionStructureConfigEntry_t   *entryArray;
    } PredictionStructureConfig_t;

    /************************************************
     * Prediction Structure Entry
     *   Contains the reference and dependent lists
     *   for a particular picture in the Prediction
     *   Structure.
     ************************************************/
    typedef struct PredictionStructureEntry_s {
        ReferenceList_t                     refList0;
        ReferenceList_t                     refList1;
        DependentList_t                     depList0;
        DependentList_t                     depList1;
        uint32_t                              temporal_layer_index;
        uint32_t                              decode_order;
        EbBool                             isReferenced;

        // High-level RPS
        EbBool                             shortTermRpsInSpsFlag;
        uint32_t                              shortTermRpsInSpsIndex;
        EbBool                             interRpsPredictionFlag;
        EbBool                             longTermRpsPresentFlag;
        uint32_t                              gopPositionLeastSignificantBits;

        // Predicted Short-Term RPS
        uint32_t                              deltaRpsIndexMinus1;
        uint32_t                              absoluteDeltaRpsMinus1;
        uint32_t                              deltaRpsSign;
        EbBool                             usedByCurrPicFlag[MAX_NUM_OF_REF_PICS_TOTAL];
        EbBool                             usedByFuturePicFlag[MAX_NUM_OF_REF_PICS_TOTAL];

        // Non-Predicted Short-Term RPS
        uint32_t                              negativeRefPicsTotalCount;
        uint32_t                              positiveRefPicsTotalCount;
        uint32_t                              deltaNegativeGopPosMinus1[MAX_NUM_OF_NEGATIVE_REF_PICS];
        uint32_t                              deltaPositiveGopPosMinus1[MAX_NUM_OF_POSITIVE_REF_PICS];
        EbBool                             usedByNegativeCurrPicFlag[MAX_NUM_OF_NEGATIVE_REF_PICS];
        EbBool                             usedByPositiveCurrPicFlag[MAX_NUM_OF_POSITIVE_REF_PICS];

        // Long-Term RPS
        uint32_t                              longTermRefPicsTotalCount;
        uint32_t                              deltaGopPoslsb[MAX_NUM_OF_REF_PICS_TOTAL];
        EbBool                             deltaGopPosMsbPresentFlag[MAX_NUM_OF_REF_PICS_TOTAL];
        uint32_t                              deltaGopPosMsbMinus1[MAX_NUM_OF_REF_PICS_TOTAL];
        EbBool                             usedByLtCurrPicFlagArray[MAX_NUM_OF_REF_PICS_TOTAL];

        // List Construction
        EbBool                             refPicsOverrideTotalCountFlag;
        int32_t                              refPicsList0TotalCountMinus1;
        int32_t                              refPicsList1TotalCountMinus1;
        EbBool                             listsModificationPresentFlag;
        EbBool                             restrictedRefPicListsFlag;      // Same list enable flag (if set,
                                                                            //   it implies all slices of the
                                                                            //   same type in the same picture
                                                                            //   have identical lists)

        // List Modification
        // *Note - This should probably be moved to the slice header since its a dynamic control - JMJ Jan 2, 2013
        EbBool                             list0ModificationFlag;
        EbBool                             list1ModificationFlag;
        uint32_t                              list0ModIndex[MAX_NUM_OF_REF_PICS_TOTAL];
        uint32_t                              list1ModIndex[MAX_NUM_OF_REF_PICS_TOTAL];

        // Lists Combination (STUB)

    } PredictionStructureEntry_t;

    /************************************************
     * Prediction Structure
     *   Contains a collection of control and RPS
     *   data types for an entire Prediction Structure
     ************************************************/
    typedef struct PredictionStructure_s {

        uint32_t                              predStructEntryCount;
        PredictionStructureEntry_t        **predStructEntryPtrArray;
        EbPred                             predType;
        uint32_t                              temporalLayerCount;
        uint32_t                              predStructPeriod;
        uint32_t                              maximumExtent;

        // Section Indices
        uint32_t                              leadingPicIndex;
        uint32_t                              initPicIndex;
        uint32_t                              steadyStateIndex;

        // RPS Related Entries
        EbBool                             restrictedRefPicListsEnableFlag;
        EbBool                             listsModificationEnableFlag;
        EbBool                             longTermEnableFlag;
        uint32_t                              defaultRefPicsList0TotalCountMinus1;
        uint32_t                              defaultRefPicsList1TotalCountMinus1;

    } PredictionStructure_t;

    /************************************************
     * Prediction Structure Group
     *   Contains the control structures for all
     *   supported prediction structures.
     ************************************************/
    typedef struct PredictionStructureGroup_s {
        PredictionStructure_t             **predictionStructurePtrArray;
        uint32_t                              predictionStructureCount;
    } PredictionStructureGroup_t;

    /************************************************
     * Declarations
     ************************************************/
    extern EbErrorType PredictionStructureGroupCtor(
        PredictionStructureGroup_t   **predictionStructureGroupDblPtr,
        uint32_t                        base_layer_switch_mode);

    extern PredictionStructure_t* GetPredictionStructure(
        PredictionStructureGroup_t    *prediction_structure_group_ptr,
        EbPred                        pred_structure,
        uint32_t                         numberOfReferences,
        uint32_t                         levelsOfHierarchy);




    typedef struct Av1RpsNode_s {
        uint8_t refreshFrameMask;
        uint8_t refDpbIndex[7];//LAST-LAST2-LAST3-GOLDEN-BWD-ALT2-ALT
    } Av1RpsNode_t;



#ifdef __cplusplus
}
#endif
#endif // EbPredictionStructure_h