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
#if !OPT_RPS_CONSTR_3
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
#endif
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
#if !OPT_RPS_CONSTR_3
    ReferenceList ref_list0;
    ReferenceList ref_list1;
    uint32_t      dep_list0_count;
    uint32_t      dep_list1_count;
#endif
    uint32_t temporal_layer_index;
    uint32_t decode_order;
#if !OPT_RPS_CONSTR_3
    Bool is_referenced;
#endif
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
    SvtAv1PredStructure        pred_type;
#if !OPT_RPS_CONSTR_3
    uint32_t temporal_layer_count;
#endif
    uint32_t pred_struct_period;
#if !OPT_RPS_CONSTR_3
    uint32_t maximum_extent;
#endif
    // Section Indices
    uint32_t init_pic_index;
#if !OPT_RPS_CONSTR_3
    uint32_t steady_state_index;
#endif
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
#if !OPT_RPS_CONSTR_3
    uint8_t ref_count_used;
#endif
} PredictionStructureGroup;

/************************************************
     * Declarations
     ************************************************/
#if !OPT_RPS_CONSTR_3
//EbErrorType prediction_structure_group_ctor(
//   PredictionStructureGroup* pred_struct_group_ptr,
//   struct SequenceControlSet* scs_ptr);
#endif
#if CLN_REMOVE_REF_CNT
extern PredictionStructure *get_prediction_structure(
    PredictionStructureGroup *prediction_structure_group_ptr, SvtAv1PredStructure pred_structure,
    uint32_t levels_of_hierarchy);
#else
extern PredictionStructure *get_prediction_structure(
    PredictionStructureGroup *prediction_structure_group_ptr, SvtAv1PredStructure pred_structure,
    uint32_t number_of_references, uint32_t levels_of_hierarchy);
#endif
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
