/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
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

    // High-level RPS
    EbBool   short_term_rps_in_sps_flag;
    uint32_t short_term_rps_in_sps_index;
    EbBool   inter_rps_prediction_flag;
    EbBool   long_term_rps_present_flag;
    uint32_t gop_position_least_significant_bits;

    // Predicted Short-Term RPS
    uint32_t delta_rps_index_minus1;
    uint32_t absolute_delta_rps_minus1;
    uint32_t delta_rps_sign;
    EbBool   used_by_curr_pic_flag[MAX_NUM_OF_REF_PICS_TOTAL];
    EbBool   used_by_future_pic_flag[MAX_NUM_OF_REF_PICS_TOTAL];
    // Non-Predicted Short-Term RPS
    uint32_t negative_ref_pics_total_count;
    uint32_t positive_ref_pics_total_count;
    uint32_t delta_negative_gop_pos_minus1[MAX_NUM_OF_NEGATIVE_REF_PICS];
    uint32_t delta_positive_gop_pos_minus1[MAX_NUM_OF_POSITIVE_REF_PICS];
    EbBool   used_by_negative_curr_pic_flag[MAX_NUM_OF_NEGATIVE_REF_PICS];
    EbBool   used_by_positive_curr_pic_flag[MAX_NUM_OF_POSITIVE_REF_PICS];

    // Long-Term RPS
    uint32_t long_term_ref_pics_total_count;
    uint32_t delta_gop_poslsb[MAX_NUM_OF_REF_PICS_TOTAL];
    EbBool   delta_gop_pos_msb_present_flag[MAX_NUM_OF_REF_PICS_TOTAL];
    uint32_t delta_gop_pos_msb_minus1[MAX_NUM_OF_REF_PICS_TOTAL];
    EbBool   used_by_lt_curr_pic_flag_array[MAX_NUM_OF_REF_PICS_TOTAL];
    // List Construction
    EbBool  ref_pics_override_total_count_flag;
    int32_t ref_pics_list0_total_count_minus1;
    int32_t ref_pics_list1_total_count_minus1;
    EbBool  lists_modification_present_flag;
    EbBool  restricted_ref_pic_lists_flag; // Same list enable flag (if set,
        //   it implies all slices of the
        //   same type in the same picture
        //   have identical lists)

    // List Modification
    // *Note - This should probably be moved to the slice header since its a dynamic control - JMJ Jan 2, 2013
    EbBool   list0_modification_flag;
    EbBool   list1_modification_flag;
    uint32_t list0_mod_index[MAX_NUM_OF_REF_PICS_TOTAL];
    uint32_t list1_mod_index[MAX_NUM_OF_REF_PICS_TOTAL];

    // Lists Combination (STUB)
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

    // RPS Related Entries
    EbBool   restricted_ref_pic_lists_enable_flag;
    EbBool   lists_modification_enable_flag;
    EbBool   long_term_enable_flag;
    uint32_t default_ref_pics_list0_total_count_minus1;
    uint32_t default_ref_pics_list1_total_count_minus1;

    //private use only
    int32_t * decode_order_table;
    uint32_t *display_order_table;
    EbBool *  timeline_map;
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
} PredictionStructureGroup;

/************************************************
     * Declarations
     ************************************************/
extern EbErrorType prediction_structure_group_ctor(PredictionStructureGroup *pred_struct_group_ptr,
#if !MAR12_M8_ADOPTIONS || (M8_MRP && !UPGRADE_M6_M7_M8)
                                                   uint8_t                   enc_mode,
#endif
                                                   EbSvtAv1EncConfiguration *config);

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
