/*
* Copyright(c) 2019 Intel Corporation
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

#ifndef EbPredictionStructure_h
#define EbPredictionStructure_h

#include "EbDefinitions.h"
#include "EbObject.h"
#ifdef __cplusplus
extern "C" {
#endif
/************************************************
     * Defines
     ************************************************/
#define FLAT_PREDICTION_STRUCTURE_PERIOD 1
#define TWO_LEVEL_HIERARCHICAL_PREDICTION_STRUCTURE_PERIOD 2
#define THREE_LEVEL_HIERARCHICAL_PREDICTION_STRUCTURE_PERIOD 4
#define FOUR_LEVEL_HIERARCHICAL_PREDICTION_STRUCTURE_PERIOD 8
#define MAX_PREDICTION_STRUCTURE_PERIOD 64

/************************************************
      * RPS defines
      ************************************************/
#define RPS_UNUSED ~0
#define MAX_NUM_OF_NEGATIVE_REF_PICS 5
#define MAX_NUM_OF_POSITIVE_REF_PICS 5
#define MAX_NUM_OF_REF_PICS_TOTAL (MAX_NUM_OF_POSITIVE_REF_PICS + MAX_NUM_OF_NEGATIVE_REF_PICS)

/************************************************
       * Reference List
       *
       *   reference_list - Contains the deltaPOCs of
       *    the pictures referenced by the current
       *    picture.
       ************************************************/
typedef struct ReferenceList {
    int32_t *reference_list;
    uint32_t reference_list_count;
} ReferenceList;

/************************************************
     * Dependent List
     *
     *   list_count - Contains count of how
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
typedef struct DependentList {
    int32_t *list;
    uint32_t list_count;
} DependentList;

/************************************************
     * Prediction Structure Config
     *   Contains a collection of basic control data
     *   for the basic prediction structure.
     ************************************************/
typedef struct PredictionStructureConfig {
    uint32_t                        entry_count;
    PredictionStructureConfigEntry *entry_array;
} PredictionStructureConfig;

/************************************************
     * Prediction Structure Entry
     *   Contains the reference and dependent lists
     *   for a particular picture in the Prediction
     *   Structure.
     ************************************************/
typedef struct PredictionStructureEntry {
    ReferenceList ref_list0;
    ReferenceList ref_list1;
    DependentList dep_list0;
    DependentList dep_list1;
    uint32_t      temporal_layer_index;
    uint32_t      decode_order;
    EbBool        is_referenced;
} PredictionStructureEntry;

/************************************************
     * Prediction Structure
     *   Contains a collection of control and RPS
     *   data types for an entire Prediction Structure
     ************************************************/
typedef struct PredictionStructure {
    EbDctor                    dctor;
    uint32_t                   pred_struct_entry_count;
    PredictionStructureEntry **pred_struct_entry_ptr_array;
    EbPred                     pred_type;
    uint32_t                   temporal_layer_count;
    uint32_t                   pred_struct_period;
    uint32_t                   maximum_extent;

    // Section Indices
    uint32_t leading_pic_index;
    uint32_t init_pic_index;
    uint32_t steady_state_index;
} PredictionStructure;

/************************************************
     * Prediction Structure Group
     *   Contains the control structures for all
     *   supported prediction structures.
     ************************************************/
typedef struct PredictionStructureGroup {
    EbDctor               dctor;
    PredictionStructure **prediction_structure_ptr_array;
    uint32_t              prediction_structure_count;
    void                 *priv; /* private member*/
    uint8_t               ref_count_used;
} PredictionStructureGroup;

/************************************************
     * Declarations
     ************************************************/
//EbErrorType prediction_structure_group_ctor(
//   PredictionStructureGroup* pred_struct_group_ptr,
//   struct SequenceControlSet* scs_ptr);
extern PredictionStructure *get_prediction_structure(
    PredictionStructureGroup *prediction_structure_group_ptr, EbPred pred_structure,
    uint32_t number_of_references, uint32_t levels_of_hierarchy);
typedef enum {
    LAST  = 0,
    LAST2 = 1,
    LAST3 = 2,
    GOLD  = 3,
    BWD   = 4,
    ALT2  = 5,
    ALT   = 6
} REF_FRAME_MINUS1;
typedef struct Av1RpsNode {
    uint8_t  refresh_frame_mask;
    uint8_t  ref_dpb_index[7]; //LAST-LAST2-LAST3-GOLDEN-BWD-ALT2-ALT
    uint64_t ref_poc_array[7]; //decoder based ref poc array //LAST-LAST2-LAST3-GOLDEN-BWD-ALT2-ALT
} Av1RpsNode;

#ifdef __cplusplus
}
#endif
#endif // EbPredictionStructure_h
