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
    uint32_t temporal_layer_index;
    uint32_t decode_order;
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
    uint32_t                   pred_struct_period;
    // Section Indices
    uint32_t init_pic_index;
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
} PredictionStructureGroup;

/************************************************
     * Declarations
     ************************************************/
extern PredictionStructure *svt_aom_get_prediction_structure(PredictionStructureGroup *prediction_structure_group_ptr,
                                                             SvtAv1PredStructure       pred_structure,
                                                             uint32_t                  levels_of_hierarchy);
typedef enum { LAST = 0, LAST2 = 1, LAST3 = 2, GOLD = 3, BWD = 4, ALT2 = 5, ALT = 6 } REF_FRAME_MINUS1;
typedef struct Av1RpsNode {
    uint8_t  refresh_frame_mask;
    uint8_t  ref_dpb_index[7]; //LAST-LAST2-LAST3-GOLDEN-BWD-ALT2-ALT
    uint64_t ref_poc_array[7]; //decoder based ref poc array //LAST-LAST2-LAST3-GOLDEN-BWD-ALT2-ALT
} Av1RpsNode;

#ifdef __cplusplus
}
#endif
#endif // EbPredictionStructure_h
