/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>
#include <string.h>

#include "EbDefinitions.h"
#include "EbPredictionStructure.h"
#include "EbUtility.h"

/**********************************************************
 * Macros
 **********************************************************/
#define PRED_STRUCT_INDEX(hierarchicalLevelCount, predType, refCount) (((hierarchicalLevelCount) * EB_PRED_TOTAL_COUNT + (predType)) * REF_LIST_MAX_DEPTH + (refCount))

 //#define DEP_INDEX(predIdx, entryIndex, entryTotalCount) ((((predIdx) - ((int32_t) entryIndex)) % (entryTotalCount)))
#define DEP_INDEX(predIdx, entryIndex, entryTotalCount) (((((int32_t) entryIndex) - (predIdx)) % (entryTotalCount)))

/**********************************************************
 * Instructions for how to create a Predicion Structure
 *
 * Overview:
 *   The prediction structure consists of a collection
 *   of Prediction Structure Entires, which themselves
 *   consist of reference and dependent lists.  The
 *   reference lists are exactly like those found in
 *   the standard and can be clipped in order to reduce
 *   the number of references.
 *
 *   Dependent lists are the corollary to reference lists,
 *   the describe how a particular picture is referenced.
 *   Dependent lists can also be clipped at predefined
 *   junctions (i.e. the list_count array) in order
 *   to reduce the number of references.  Note that the
 *   dependent deltaPOCs must be grouped together in order
 *   of ascending referencePicture in order for the Dependent
 *   List clip to work properly.
 *
 *   All control and RPS information is derived from
 *   these lists.  The lists for a structure are defined
 *   for both P & B-picture variants.  In the case of
 *   P-pictures, only Lists 0 are used.
 *
 *   Negative deltaPOCs are for backward-referencing pictures
 *   in display order and positive deltaPOCs are for
 *   forward-referencing pictures.
 *
 *   Please note that there is no assigned coding order,
 *   the PictureManager will start pictures as soon as
 *   their references become available.
 *
 *   Any prediction structure is possible; however, we are
 *     restricting usage to the following controls:
 *     # Hierarchical Levels
 *     # Number of References
 *     # B-pictures enabled
 *     # Intra Refresh Period
 *
 *  To Get Low Delay P, only use List 0
 *  To Get Low Delay B, replace List 1 with List 0
 *  To Get Random Access, use the preduction structure as is
 **********************************************************/

 /************************************************
  * Flat
  *
  *  I-B-B-B-B-B-B-B-B
  *
  * Display & Coding Order:
  *  0 1 2 3 4 5 6 7 8
  *
  ************************************************/
static PredictionStructureConfigEntry flat_pred_struct[] = {
    {
        0,               // GOP Index 0 - Temporal Layer
        0,               // GOP Index 0 - Decode Order
        {1, 2, 3, 4},    // GOP Index 0 - Ref List 0
        {1, 2, 3, 4}     // GOP Index 0 - Ref List 1
    }
};

/************************************************
 * Random Access - Two-Level Hierarchical
 *
 *    b   b   b   b      Temporal Layer 1
 *   / \ / \ / \ / \
 *  I---B---B---B---B    Temporal Layer 0
 *
 * Display Order:
 *  0 1 2 3 4 5 6 7 8
 *
 * Coding Order:
 *  0 2 1 4 3 6 5 8 7
 ************************************************/
static PredictionStructureConfigEntry two_level_hierarchical_pred_struct[] = {
    {
        0,                // GOP Index 0 - Temporal Layer
        0,                // GOP Index 0 - Decode Order
        {2, 4, 6, 8},     // GOP Index 0 - Ref List 0
        {2, 4, 6, 8}      // GOP Index 0 - Ref List 1
    },
    {
        1,                // GOP Index 1 - Temporal Layer
        1,                // GOP Index 1 - Decode Order
        { 1, 3, 5, 7},    // GOP Index 1 - Ref List 0
        {-1, 1, 3, 5}     // GOP Index 1 - Ref List 1
    }
};

/************************************************
 * Three-Level Hierarchical
 *
 *      b   b       b   b       b   b        Temporal Layer 2
 *     / \ / \     / \ / \     / \ / \
 *    /   B   \   /   B   \   /   B   \      Temporal Layer 1
 *   /   / \   \ /   / \   \ /   / \   \
 *  I-----------B-----------B-----------B    Temporal Layer 0
 *
 * Display Order:
 *  0   1 2 3   4   5 6 7   8   9 1 1   1
 *                                0 1   2
 *
 * Coding Order:
 *  0   3 2 4   1   7 6 8   5   1 1 1   9
 *                              1 0 2
 ************************************************/
static PredictionStructureConfigEntry three_level_hierarchical_pred_struct[] = {
    {
        0,                  // GOP Index 0 - Temporal Layer
        0,                  // GOP Index 0 - Decode Order
        {4, 8, 12, 16},     // GOP Index 0 - Ref List 0
        {4, 8, 12, 16}      // GOP Index 0 - Ref List 1
    },
    {
        2,                  // GOP Index 1 - Temporal Layer
        2,                  // GOP Index 1 - Decode Order
        { 1,  3, 5, 0},     // GOP Index 1 - Ref List 0
        {-1, -3, 1, 3}      // GOP Index 1 - Ref List 1
    },
    {
        1,                  // GOP Index 2 - Temporal Layer
        1,                  // GOP Index 2 - Decode Order
        { 2, 6, 10, 14},    // GOP Index 2 - Ref List 0
        {-2, 2,  6, 10}     // GOP Index 2 - Ref List 1
    },
    {
        2,                  // GOP Index 3 - Temporal Layer
        3,                  // GOP Index 3 - Decode Order
        { 1, 3, 5, 7},      // GOP Index 3 - Ref List 0
        {-1, 1, 3, 5}       // GOP Index 3 - Ref List 1
    }
};

/************************************************************************************************************
 * Four-Level Hierarchical
 *
 *
 *          b     b           b     b               b     b           b     b           Temporal Layer 3
 *         / \   / \         / \   / \             / \   / \         / \   / \
 *        /   \ /   \       /   \ /   \           /   \ /   \       /   \ /   \
 *       /     B     \     /     B     \         /     B     \     /     B     \        Temporal Layer 2
 *      /     / \     \   /     / \     \       /     / \     \   /     / \     \
 *     /     /   \     \ /     /   \     \     /     /   \     \ /     /   \     \
 *    /     /     ------B------     \     \   /     /     ------B------     \     \     Temporal Layer 1
 *   /     /           / \           \     \ /     /           / \           \     \
 *  I---------------------------------------B---------------------------------------B   Temporal Layer 0
 *
 * Display Order:
 *  0       1  2  3     4     5  6  7       8       9  1  1     1     1  1  1       1
 *                                                     0  1     2     3  4  5       6
 *
 * Coding Order:
 *  0       4  3  5     2     7  6  8       1       1  1  1     1     1  1  1       9
 *                                                  2  1  3     0     5  4  6
 *
 ***********************************************************************************************************/
PredictionStructureConfigEntry four_level_hierarchical_pred_struct[] = {
    {
        0,                      // GOP Index 0 - Temporal Layer
        0,                      // GOP Index 0 - Decode Order
        { 8, 0, 0, 0},          // GOP Index 0 - Ref List 0
        { 8, 0, 0, 0}           // GOP Index 0 - Ref List 1
    },
    {
        3,                      // GOP Index 1 - Temporal Layer
        3,                      // GOP Index 1 - Decode Order
#if PRED_CHANGE
        { 1, 3, 5, 8},          // GOP Index 1 - Ref List 0
#else
        { 1, 3, 5, 9},          // GOP Index 1 - Ref List 0
#endif
        {-1,-3,-7, 0}           // GOP Index 1 - Ref List 1
    },
    {
        2,                      // GOP Index 2 - Temporal Layer
        2,                      // GOP Index 2 - Decode Order
        { 2, 4, 6, 10},         // GOP Index 2 - Ref List 0
        {-2,-6, 0,  0}          // GOP Index 2 - Ref List 1
    },
    {
        3,                      // GOP Index 3 - Temporal Layer
        4,                      // GOP Index 3 - Decode Order
#if PRED_CHANGE_MOD
        { 1, 3, 2, 5},          // GOP Index 3 - Ref List 0
#elif PRED_CHANGE
        { 1, 2, 3, 5},          // GOP Index 3 - Ref List 0
#else
        { 1, 3, 5, 7},          // GOP Index 3 - Ref List 0
#endif
        {-1, -5, 0, 0}          // GOP Index 3 - Ref List 1
    },
    {
        1,                      // GOP Index 4 - Temporal Layer
        1,                      // GOP Index 4 - Decode Order
        { 4, 8, 12, 0},         // GOP Index 4 - Ref List 0
        {-4, 0,  0, 0}          // GOP Index 4 - Ref List 1
    },
    {
        3,                      // GOP Index 5 - Temporal Layer
        6,                      // GOP Index 5 - Decode Order
#if PRED_CHANGE
        { 1, 3, 5, 4},          // GOP Index 5 - Ref List 0
#else
        { 1, 3, 5, 9},          // GOP Index 5 - Ref List 0
#endif
        {-1, -3, 0, 0}          // GOP Index 5 - Ref List 1
    },
    {
        2,                      // GOP Index 6 - Temporal Layer
        5,                      // GOP Index 6 - Decode Order
        { 2, 4, 6, 10},         // GOP Index 6 - Ref List 0
        {-2, 0, 0,  0}          // GOP Index 6 - Ref List 1
    },
    {
        3,                      // GOP Index 7 - Temporal Layer
        7,                      // GOP Index 7 - Decode Order
#if PRED_CHANGE
        { 1, 3, 5, 6},          // GOP Index 7 - Ref List 0
#else
        { 1, 3, 5, 7},          // GOP Index 7 - Ref List 0
#endif
        {-1, 0, 0, 0}           // GOP Index 7 - Ref List 1
    }
};

/***********************************************************************************************************
 * Five-Level Level Hierarchical
 *
 *           b     b           b     b               b     b           b     b              Temporal Layer 4
 *          / \   / \         / \   / \             / \   / \         / \   / \
 *         /   \ /   \       /   \ /   \           /   \ /   \       /   \ /   \
 *        /     B     \     /     B     \         /     B     \     /     B     \           Temporal Layer 3
 *       /     / \     \   /     / \     \       /     / \     \   /     / \     \
 *      /     /   \     \ /     /   \     \     /     /   \     \ /     /   \     \
 *     /     /     ------B------     \     \   /     /     ------B------     \     \        Temporal Layer 2
 *    /     /           / \           \     \ /     /           / \           \     \
 *   /     /           /   \-----------------B------------------   \           \     \      Temporal Layer 1
 *  /     /           /                     / \                     \           \     \
 * I-----------------------------------------------------------------------------------B    Temporal Layer 0
 *
 * Display Order:
 *  0        1  2  3     4     5  6  7       8       9  1  1     1     1  1  1         1
 *                                                      0  1     2     3  4  5         6
 *
 * Coding Order:
 *  0        5  4  6     3     8  7  9       2       1  1  1     1     1  1  1         1
 *                                                   2  1  3     0     5  4  6
 *
 ***********************************************************************************************************/
PredictionStructureConfigEntry five_level_hierarchical_pred_struct[] = {

    {
        0,                      // GOP Index 0 - Temporal Layer
        0,                      // GOP Index 0 - Decode Order
        {16, 48, 0, 0},         // GOP Index 0 - Ref List 0
        {16, 32, 0, 0}          // GOP Index 0 - Ref List 1     //we need keep 16 as first entry in List1, this will ensure that for poc=16 there is a valid ref frame in List1.
    },
    {
        4,                      // GOP Index 1 - Temporal Layer
        4,                      // GOP Index 1 - Decode Order
#if PRED_CHANGE_MOD
        { 1, 9, 8, 17},          // GOP Index 1 - Ref List 0
#elif PRED_CHANGE_5L
        { 1, 8, 9, 17},          // GOP Index 1 - Ref List 0
#else
        { 1, 5, 9, 17},          // GOP Index 1 - Ref List 0
#endif
        { -1, -3, -7, 0}         // GOP Index 1 - Ref List 1
    },
    {
        3,                      // GOP Index 2 - Temporal Layer
        3,                      // GOP Index 2 - Decode Order
#if PRED_CHANGE_5L
        { 2, 4, 10, 18},        // GOP Index 2 - Ref List 0
#else
        { 2, 4, 6, 10},         // GOP Index 2 - Ref List 0
#endif
        { -2, -6, -14, 0}       // GOP Index 2 - Ref List 1
    },
    {
        4,                      // GOP Index 3 - Temporal Layer
        5,                      // GOP Index 3 - Decode Order
#if PRED_CHANGE_MOD
        { 1, 3, 2, 11},         // GOP Index 3 - Ref List 0
#elif PRED_CHANGE_5L
        { 1, 2, 3, 11},         // GOP Index 3 - Ref List 0
#else
        { 1, 3, 7, 11},         // GOP Index 3 - Ref List 0
#endif
        { -1, -5, -13, 0}       // GOP Index 3 - Ref List 1
    },
    {
        2,                      // GOP Index 4 - Temporal Layer
        2,                      // GOP Index 4 - Decode Order
        {  4, 8, 12, 20},       // GOP Index 4 - Ref List 0
        { -4, -12, 0, 0}        // GOP Index 4 - Ref List 1
    },
    {
        4,                      // GOP Index 5 - Temporal Layer
        7,                      // GOP Index 5 - Decode Order
#if PRED_CHANGE_MOD
        { 1, 5, 4, 13},         // GOP Index 5 - Ref List 0
#elif PRED_CHANGE_5L
        { 1, 4, 5, 13},         // GOP Index 5 - Ref List 0
#else
        { 1, 5, 9, 13},         // GOP Index 5 - Ref List 0
#endif
        { -1, -3, -11, 0}       // GOP Index 5 - Ref List 1
    },
    {
        3,                      // GOP Index 6 - Temporal Layer
        6,                      // GOP Index 6 - Decode Order
#if PRED_CHANGE_5L
        { 2, 4, 6, 14},         // GOP Index 6 - Ref List 0
#else
        { 2, 4, 6, 10},         // GOP Index 6 - Ref List 0
#endif
        { -2, -10, 0, 0}        // GOP Index 6 - Ref List 1
    },
    {
        4,                      // GOP Index 7 - Temporal Layer
        8,                      // GOP Index 7 - Decode Order
#if PRED_CHANGE_MOD
        { 1, 3, 6, 7},          // GOP Index 7 - Ref List 0
#elif PRED_CHANGE_5L
        { 1, 3, 6, 7},          // GOP Index 7 - Ref List 0
#else
        { 1, 3, 7, 11},         // GOP Index 7 - Ref List 0
#endif
        { -1, -9, 0, 0}         // GOP Index 7 - Ref List 1
    },
    {
        1,                      // GOP Index 8 - Temporal Layer
        1,                      // GOP Index 8 - Decode Order
        {  8, 16, 24, 0},       // GOP Index 8 - Ref List 0
        { -8, 0, 0, 0}          // GOP Index 8 - Ref List 1
    },
    {
        4,                      // GOP Index 9 - Temporal Layer
        11,                     // GOP Index 9 - Decode Order
#if PRED_CHANGE_MOD
        { 1, 9, 8, 17},         // GOP Index 9 - Ref List 0
#elif PRED_CHANGE_5L
        { 1, 8, 9, 17},         // GOP Index 9 - Ref List 0
#else
        { 1, 5, 9, 17},         // GOP Index 9 - Ref List 0
#endif
        { -1, -3, -7, 0}        // GOP Index 9 - Ref List 1
    },
    {
        3,                      // GOP Index 10 - Temporal Layer
        10,                     // GOP Index 10 - Decode Order
#if PRED_CHANGE_5L
        { 2, 4, 10, 18},        // GOP Index 10 - Ref List 0
#else
        { 2, 4, 6, 10},         // GOP Index 10 - Ref List 0
#endif
        { -2, -6, 0, 0}         // GOP Index 10 - Ref List 1
    },
    {
        4,                      // GOP Index 11 - Temporal Layer
        12,                     // GOP Index 11 - Decode Order
#if PRED_CHANGE_MOD
        { 1, 3, 2, 11},         // GOP Index 11 - Ref List 0
#elif PRED_CHANGE_5L
        { 1, 2, 3, 11},         // GOP Index 11 - Ref List 0
#else
        { 1, 3, 7, 11},         // GOP Index 11 - Ref List 0
#endif
        { -1, -5, 0, 0}         // GOP Index 11 - Ref List 1
    },
    {
        2,                      // GOP Index 12 - Temporal Layer
        9,                      // GOP Index 12 - Decode Order
        { 4, 8, 12, 0},         // GOP Index 12 - Ref List 0
        { -4, 0, 0, 0}          // GOP Index 12 - Ref List 1
    },
    {
        4,                      // GOP Index 13 - Temporal Layer
        14,                     // GOP Index 13 - Decode Order
#if PRED_CHANGE_MOD
        { 1, 5, 4, 13},         // GOP Index 13 - Ref List 0
#elif PRED_CHANGE_5L
        { 1, 4, 5, 13},         // GOP Index 13 - Ref List 0
#else
        { 1, 5, 9, 13},         // GOP Index 13 - Ref List 0
#endif
        { -1, -3, 0, 0}         // GOP Index 13 - Ref List 1
    },
    {
        3,                      // GOP Index 14 - Temporal Layer
        13,                     // GOP Index 14 - Decode Order
        { 2, 4, 6, 14},         // GOP Index 14 - Ref List 0
        { -2, 0, 0, 0}          // GOP Index 14 - Ref List 1

    },
    {
        4,                      // GOP Index 15 - Temporal Layer
        15,                     // GOP Index 15 - Decode Order
#if PRED_CHANGE_MOD
        { 1, 3, 6, 7},          // GOP Index 15 - Ref List 0
#elif PRED_CHANGE_5L
        { 1, 3, 6, 7},          // GOP Index 15 - Ref List 0
#else
        { 1, 3, 7, 11},         // GOP Index 15 - Ref List 0
#endif
        { -1, 0, 0, 0}          // GOP Index 15 - Ref List 1
    }
};

/**********************************************************************************************************************************************************************************************************************
 * Six-Level Level Hierarchical
 *
 *
 *              b     b           b     b               b     b           b     b                   b     b           b     b               b     b           b     b               Temporal Layer 5
 *             / \   / \         / \   / \             / \   / \         / \   / \                 / \   / \         / \   / \             / \   / \         / \   / \
 *            /   \ /   \       /   \ /   \           /   \ /   \       /   \ /   \               /   \ /   \       /   \ /   \           /   \ /   \       /   \ /   \
 *           /     B     \     /     B     \         /     B     \     /     B     \             /     B     \     /     B     \         /     B     \     /     B     \            Temporal Layer 4
 *          /     / \     \   /     / \     \       /     / \     \   /     / \     \           /     / \     \   /     / \     \       /     / \     \   /     / \     \
 *         /     /   \     \ /     /   \     \     /     /   \     \ /     /   \     \         /     /   \     \ /     /   \     \     /     /   \     \ /     /   \     \
 *        /     /     ------B------     \     \   /     /     ------B------     \     \       /     /     ------B------     \     \   /     /     ------B------     \     \         Temporal Layer 3
 *       /     /           / \           \     \ /     /           / \           \     \     /     /           / \           \     \ /     /           / \           \     \
 *      /     /           /   \-----------------B------------------   \           \     \   /     /           /   \-----------------B------------------   \           \     \       Temporal Layer 2
 *     /     /           /                     / \                     \           \     \ /     /           /                     / \                     \           \     \
 *    /     /           /                     /   \---------------------------------------B---------------------------------------/   \                     \           \     \     Temporal Layer 1
 *   /     /           /                     /                                           / \                                           \                     \           \     \
 *  I---------------------------------------------------------------------------------------------------------------------------------------------------------------------------B   Temporal Layer 0
 *
 * Display Order:
 *  0           1  2  3     4     5  6  7       8       9  1  1     1     1  1  1         1         1  1  1     2     2  2  2       2       2  2  2     2     2  3  3           3
 *                                                         0  1     2     3  4  5         6         7  8  9     0     1  2  3       4       5  6  7     8     9  0  1           2
 *
 * Coding Order:
 *  0           6  5  7     4     9  8  1       3       1  1  1     1     1  1  1         2         2  2  2     1     2  2  2       1       2  2  2     2     3  3  3           1
 *                                      0               3  2  4     1     6  5  7                   1  0  2     9     4  3  5       8       8  7  9     6     1  0  2
 *
 **********************************************************************************************************************************************************************************************************************/
static PredictionStructureConfigEntry six_level_hierarchical_pred_struct[] = {
    {
        0,                  // GOP Index 0 - Temporal Layer
        0,                  // GOP Index 0 - Decode Order
        {32,  0, 0, 0},     // GOP Index 0 - Ref List 0
        {32,  0, 0, 0}      // GOP Index 0 - Ref List 1
    },
    {
        5,                  // GOP Index 1 - Temporal Layer
        5,                  // GOP Index 1 - Decode Order
        { 1,  0,  0,   0},  // GOP Index 1 - Ref List 0
        {-1, -3, -7, -15}   // GOP Index 1 - Ref List 1
    },
    {
        4,                  // GOP Index 2 - Temporal Layer
        4,                  // GOP Index 2 - Decode Order
        {2,   0,   0,   0}, // GOP Index 2 - Ref List 0
        {-2, -6, -14, -30}  // GOP Index 2 - Ref List 1
    },
    {
        5,                  // GOP Index 3 - Temporal Layer
        6,                  // GOP Index 3 - Decode Order
        { 1,  3,   0,   0}, // GOP Index 3 - Ref List 0
        {-1, -5, -13,   0}  // GOP Index 3 - Ref List 1
    },
    {
        3,                  // GOP Index 4 - Temporal Layer
        3,                  // GOP Index 4 - Decode Order
        { 4,   0,   0,  0}, // GOP Index 4 - Ref List 0
        {-4, -12, -28,  4}  // GOP Index 4 - Ref List 1
    },
    {
        5,                  // GOP Index 5 - Temporal Layer
        8,                  // GOP Index 5 - Decode Order
        { 1,  5,   0, 0},   // GOP Index 5 - Ref List 0
        {-1, -3, -11, 0}    // GOP Index 5 - Ref List 1
    },
    {
        4,                  // GOP Index 6 - Temporal Layer
        7,                  // GOP Index 6 - Decode Order
        { 2,   6,   0, 0},  // GOP Index 6 - Ref List 0
        {-2, -10, -26, 2}   // GOP Index 6 - Ref List 1
    },
    {
        5,                  // GOP Index 7 - Temporal Layer
        9,                  // GOP Index 7 - Decode Order
        { 1,  3,   0, 0},   // GOP Index 7 - Ref List 0
        {-1, -9, -25, 1}    // GOP Index 7 - Ref List 1
    },
    {
        2,                  // GOP Index 8 - Temporal Layer
        2,                  // GOP Index 8 - Decode Order
        { 8,   0, 0, 0},    // GOP Index 8 - Ref List 0
        {-8, -24, 8, 0}     // GOP Index 8 - Ref List 1
    },
    {
        5,                  // GOP Index 9 - Temporal Layer
        12,                 // GOP Index 9 - Decode Order
        { 1,  9,  0, 0},    // GOP Index 9 - Ref List 0
        {-1, -3, -7, 0}     // GOP Index 9 - Ref List 1
    },
    {
        4,                  // GOP Index 10 - Temporal Layer
        11,                 // GOP Index 10 - Decode Order
        {2,  10,   0,  0},  // GOP Index 10 - Ref List 0
        {-2, -6, -22,  2}   // GOP Index 10 - Ref List 1
    },
    {
        5,                  // GOP Index 11 - Temporal Layer
        13,                 // GOP Index 11 - Decode Order
        { 1,  3,   0, 0},   // GOP Index 11 - Ref List 0
        {-1, -5, -21, 1}    // GOP Index 11 - Ref List 1
    },
    {
        3,                  // GOP Index 12 - Temporal Layer
        10,                 // GOP Index 12 - Decode Order
        { 4,  12,  0,  0},  // GOP Index 12 - Ref List 0
        {-4, -20,  4, 12}   // GOP Index 12 - Ref List 1
    },
    {
        5,                  // GOP Index 13 - Temporal Layer
        15,                 // GOP Index 13 - Decode Order
        {  1,  5,   0, 0},  // GOP Index 13 - Ref List 0
        { -1, -3, -19, 1}   // GOP Index 13 - Ref List 1
    },
    {
        4,                  // GOP Index 14 - Temporal Layer
        14,                 // GOP Index 14 - Decode Order
        { 2,   6, 14,  0},  // GOP Index 14 - Ref List 0
        {-2, -18,  2,  6}   // GOP Index 14 - Ref List 1
    },
    {
        5,                  // GOP Index 15 - Temporal Layer
        16,                 // GOP Index 15 - Decode Order
        { 1,   3, 7,  0},   // GOP Index 15 - Ref List 0
        {-1, -17, 1,  3}    // GOP Index 15 - Ref List 1
    },
    {
        1,                  // GOP Index 16 - Temporal Layer
        1,                  // GOP Index 16 - Decode Order
        { 16,  0, 0, 0},    // GOP Index 16 - Ref List 0
        {-16, 16, 0, 0}     // GOP Index 16 - Ref List 1
    },
    {
        5,                  // GOP Index 17 - Temporal Layer
        20,                 // GOP Index 17 - Decode Order
        { 1, 17,  0,  0},   // GOP Index 17 - Ref List 0
        {-1, -3, -7,  0}    // GOP Index 17 - Ref List 1
    },
    {
        4,                  // GOP Index 18 - Temporal Layer
        19,                 // GOP Index 18 - Decode Order
        { 2, 18,   0,  0},  // GOP Index 18 - Ref List 0
        {-2, -6, -14,  2}   // GOP Index 18 - Ref List 1
    },
    {
        5,                  // GOP Index 19 - Temporal Layer
        21,                 // GOP Index 19 - Decode Order
        { 1,  3,   0, 0},   // GOP Index 19 - Ref List 0
        {-1, -5, -13, 1}    // GOP Index 19 - Ref List 1
    },
    {
        3,                  // GOP Index 20 - Temporal Layer
        18,                 // GOP Index 20 - Decode Order
        { 4,  20, 0,  0},   // GOP Index 20 - Ref List 0
        {-4, -12, 4, 20}    // GOP Index 20 - Ref List 1
    },
    {
        5,                  // GOP Index 21 - Temporal Layer
        23,                 // GOP Index 21 - Decode Order
        { 1,  5,   0, 0},   // GOP Index 21 - Ref List 0
        {-1, -3, -11, 1}    // GOP Index 21 - Ref List 1
    },
    {
        4,                  // GOP Index 22 - Temporal Layer
        22,                 // GOP Index 22 - Decode Order
        { 2,   6, 22, 0},   // GOP Index 22 - Ref List 0
        {-2, -10,  2, 6}    // GOP Index 22 - Ref List 1
    },
    {
        5,                  // GOP Index 23 - Temporal Layer
        24,                 // GOP Index 23 - Decode Order
        { 1,  3, 7,  0},    // GOP Index 23 - Ref List 0
        {-1, -9, 1,  3}     // GOP Index 23 - Ref List 1
    },
    {
        2,                  // GOP Index 24 - Temporal Layer
        17,                 // GOP Index 24 - Decode Order
        { 8, 24,  0, 0},    // GOP Index 24 - Ref List 0
        {-8,  8, 24, 0}     // GOP Index 24 - Ref List 1
    },
    {
        5,                  // GOP Index 25 - Temporal Layer
        27,                 // GOP Index 25 - Decode Order
        { 1,  9,  0, 0},    // GOP Index 25 - Ref List 0
        {-1, -3, -7, 1}     // GOP Index 25 - Ref List 1
    },
    {
        4,                  // GOP Index 26 - Temporal Layer
        26,                 // GOP Index 26 - Decode Order
        { 2, 10, 26,  0},   // GOP Index 26 - Ref List 0
        {-2, -6,  2, 10}    // GOP Index 26 - Ref List 1
    },
    {
        5,                  // GOP Index 27 - Temporal Layer
        28,                 // GOP Index 27 - Decode Order
        { 1,  3, 11,  0},   // GOP Index 27 - Ref List 0
        {-1, -5,  1,  3}    // GOP Index 27 - Ref List 1
    },
    {
        3,                  // GOP Index 28 - Temporal Layer
        25,                 // GOP Index 28 - Decode Order
        { 4, 12, 28,  0},   // GOP Index 28 - Ref List 0
        {-4,  4, 12, 28}    // GOP Index 28 - Ref List 1
    },
    {
        5,                  // GOP Index 29 - Temporal Layer
        30,                 // GOP Index 29 - Decode Order
        { 1,  5, 13,  0},   // GOP Index 29 - Ref List 0
        {-1, -3,  1,  5}    // GOP Index 29 - Ref List 1
    },
    {
        4,                  // GOP Index 30 - Temporal Layer
        29,                 // GOP Index 30 - Decode Order
        { 2, 6, 14,  0},    // GOP Index 30 - Ref List 0
        {-2, 2,  6, 14}     // GOP Index 30 - Ref List 1
    },
    {
        5,                  // GOP Index 31 - Temporal Layer
        31,                 // GOP Index 31 - Decode Order
        { 1, 3, 5, 15},     // GOP Index 31 - Ref List 0
        {-1, 1, 3,  5}      // GOP Index 31 - Ref List 1
    }
};

/************************************************
 * Prediction Structure Config Array
 ************************************************/
static const PredictionStructureConfig PredictionStructureConfigArray[] = {
    {1,     flat_pred_struct},
    {2,     two_level_hierarchical_pred_struct},
    {4,     three_level_hierarchical_pred_struct},
    {8,     four_level_hierarchical_pred_struct},
    {16,    five_level_hierarchical_pred_struct},
    {32,    six_level_hierarchical_pred_struct},
    {0,     (PredictionStructureConfigEntry*)EB_NULL} // Terminating Code, must always come last!
};

/************************************************
 * Get Prediction Structure
 ************************************************/
PredictionStructure* get_prediction_structure(
    PredictionStructureGroup     *predictionStructureGroupPtr,
    EbPred                         predStructure,
    uint32_t                          numberOfReferences,
    uint32_t                          levelsOfHierarchy)
{
    PredictionStructure *predStructPtr;
    uint32_t predStructIndex;

    // Convert numberOfReferences to an index
    --numberOfReferences;

    // Determine the Index value
    predStructIndex = PRED_STRUCT_INDEX(levelsOfHierarchy, (uint32_t)predStructure, numberOfReferences);

    predStructPtr = predictionStructureGroupPtr->prediction_structure_ptr_array[predStructIndex];

    return predStructPtr;
}

static void PredictionStructureDctor(EbPtr p)
{
    PredictionStructure *obj = (PredictionStructure*)p;
    PredictionStructureEntry** pe = obj->pred_struct_entry_ptr_array;
    uint32_t count = obj->pred_struct_entry_count;
    if (pe) {
        for (uint32_t i = 0; i < count; i++) {
            EB_FREE_ARRAY(pe[i]->ref_list0.reference_list);
            EB_FREE_ARRAY(pe[i]->ref_list1.reference_list);
            EB_FREE_ARRAY(pe[i]->dep_list0.list);
            EB_FREE_ARRAY(pe[i]->dep_list1.list);
        }
        EB_FREE_2D(obj->pred_struct_entry_ptr_array);
    }
    EB_FREE_ARRAY(obj->decodeOrderTable);
    EB_FREE_ARRAY(obj->displayOrderTable);
    EB_FREE_ARRAY(obj->timelineMap);
}


/********************************************************************************************
 * Prediction Structure Ctor
 *
 * GOP Type:
 *   For Low Delay P, eliminate Ref List 0
 *   Fow Low Delay B, copy Ref List 0 into Ref List 1
 *   For Random Access, leave config as is
 *
 * numberOfReferences:
 *   Clip the Ref Lists
 *
 *  Summary:
 *
 *  The Pred Struct Ctor constructs the Reference Lists, Dependent Lists, and RPS for each
 *    valid prediction structure position. The full prediction structure is composed of four
 *    sections:
 *    a. Leading Pictures
 *    b. Initialization Pictures
 *    c. Steady-state Pictures
 *    d. Trailing Pictures
 *
 *  By definition, the Prediction Structure Config describes the Steady-state Picture
 *    Set. From the PS Config, the other sections are determined by following a simple
 *    set of construction rules. These rules are:
 *    -Leading Pictures use only List 1 of the Steady-state for forward-prediction
 *    -Init Pictures don't violate CRA mechanics
 *    -Steady-state Pictures come directly from the PS Config
 *    -Following pictures use only List 0 of the Steady-state for rear-prediction
 *
 *  In general terms, Leading and Trailing pictures are useful when trying to reduce
 *    the number of base-layer pictures in the presense of scene changes.  Trailing
 *    pictures are also useful for terminating sequences.  Init pictures are needed
 *    when using multiple references that expand outside of a Prediction Structure.
 *    Steady-state pictures are the normal use cases.
 *
 *  Leading and Trailing Pictures are not applicable to Low Delay prediction structures.
 *
 *  Below are a set of example PS diagrams
 *
 *  Low-delay P, Flat, 2 reference:
 *
 *                    I---P---P
 *
 *  Display Order     0   1   2
 *
 *  Sections:
 *    Let PredStructSize = N
 *    Leading Pictures:     [null]  Size: 0
 *    Init Pictures:        [0-1]   Size: Ceil(MaxReference, N) - N + 1
 *    Stead-state Pictures: [2]     Size: N
 *    Trailing Pictures:    [null]  Size: 0
 *    ------------------------------------------
 *      Total Size: Ceil(MaxReference, N) + 1
 *
 *  Low-delay B, 3-level, 2 references:
 *
 *                         b   b     b   b
 *                        /   /     /   /
 *                       /   B     /   B
 *                      /   /     /   /
 *                     I---------B-----------B
 *
 *  Display Order      0   1 2 3 4   5 6 7   8
 *
 *  Sections:
 *    Let PredStructSize = N
 *    Leading Pictures:     [null]  Size: 0
 *    Init Pictures:        [1-4]   Size: Ceil(MaxReference, N) - N + 1
 *    Stead-state Pictures: [5-8]   Size: N
 *    Trailing Pictures:    [null]  Size: 0
 *    ------------------------------------------
 *      Total Size: Ceil(MaxReference, N) + 1
 *
 *  Random Access, 3-level structure with 3 references:
 *
 *                   p   p       b   b       b   b       b   b       p   p
 *                    \   \     / \ / \     / \ / \     / \ / \     /   /
 *                     P   \   /   B   \   /   B   \   /   B   \   /   P
 *                      \   \ /   / \   \ /   / \   \ /   / \   \ /   /
 *                       ----I-----------B-----------B-----------B----
 *  Display Order:   0 1 2   3   4 5 6   7   8 9 1   1   1 1 1   1   1 1 1
 *                                               0   1   2 3 4   5   6 7 8
 *
 *  Decode Order:    2 1 3   0   6 5 7   4   1 9 1   8   1 1 1   1   1 1 1
 *                                           0   1       4 3 5   2   6 7 8
 *
 *  Sections:
 *    Let PredStructSize = N
 *    Leading Pictures:      [0-2]   Size: N - 1
 *    Init Pictures:         [3-11]  Size: Ceil(MaxReference, N) - N + 1
 *    Steady-state Pictures: [12-15] Size: N
 *    Trailing Pictures:     [16-18] Size: N - 1
 *    ------------------------------------------
 *      Total Size: 2*N + Ceil(MaxReference, N) - 1
 *
 *  Encoding Order:
 *                   -------->----------->----------->-----------|------->
 *                                                   |           |
 *                                                   |----<------|
 *
 *
 *  Timeline:
 *
 *  The timeline is a tool that is used to determine for how long a
 *    picture should be preserved for future reference. The concept of
 *    future reference is equivalently defined as dependence.  The RPS
 *    mechanism works by signaling for each picture in the DPB whether
 *    it is used for direct reference, kept for future reference, or
 *    discarded.  The timeline merely provides a means of determing
 *    the reference picture states at each prediction structure position.
 *    Its also important to note that all signaling should be done relative
 *    to decode order, not display order. Display order is irrelevant except
 *    for signaling the POC.
 *
 *  Timeline Example: 3-Level Hierarchical with Leading and Trailing Pictures, 3 References
 *
 *                   p   p       b   b       b   b       b   b       p   p   Temporal Layer 2
 *                    \   \     / \ / \     / \ / \     / \ / \     /   /
 *                     P   \   /   B   \   /   B   \   /   B   \   /   P     Temporal Layer 1
 *                      \   \ /   / \   \ /   / \   \ /   / \   \ /   /
 *                       ----I-----------B-----------B-----------B----       Temporal Layer 0
 *
 *  Decode Order:    2 1 3   0   6 5 7   4   1 9 1   8   1 1 1   1   1 1 1
 *                                           0   1       4 3 5   2   6 7 8
 *
 *  Display Order:   0 1 2   3   4 5 6   7   8 9 1   1   1 1 1   1   1 1 1
 *                                               0   1   2 3 4   5   6 7 8
 *            X --->
 *
 *                                1 1 1 1 1 1 1 1 1   DECODE ORDER
 *            0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8
 *
 *         |----------------------------------------------------------
 *         |
 *  Y D  0 |  \ x---x-x-x-x---x-x-x---x-x-x
 *  | E  1 |    \ x
 *  | C  2 |      \
 *  | O  3 |        \
 *  v D  4 |          \ x---x-x-x-x---x-x-x---x-x-x
 *    E  5 |            \ x-x
 *       6 |              \
 *    O  7 |                \
 *    R  8 |                  \ x---x-x-x-x---x-x-x
 *    D  9 |                    \ x-x
 *    E 10 |                      \
 *    R 11 |                        \
 *      12 |                          \ x---x-x-x-x
 *      13 |                            \ x-x
 *      14 |                              \
 *      15 |                                \
 *      16 |                                  \
 *      17 |                                    \ x
 *      18 |                                      \
 *
 *  Interpreting the timeline:
 *
 *  The most important detail to keep in mind is that all signaling
 *    is done in Decode Order space. The symbols mean the following:
 *    'x' directly referenced picture
 *    '-' picture kept for future reference
 *    ' ' not referenced, inferred discard
 *    '\' eqivalent to ' ', deliminiter that nothing can be to the left of
 *
 *  The basic steps for constructing the timeline are to increment through
 *    each position in the prediction structure (Y-direction on the timeline)
 *    and mark the appropriate state: directly referenced, kept for future reference,
 *    or discarded.  As shown, base-layer pictures are referenced much more
 *    frequently than by the other layers.
 *
 *  The RPS is constructed by looking at each 'x' position in the timeline and
 *    signaling each 'y' reference as depicted in the timeline. DPB analysis is
 *    fairly straigtforward - the total number of directly-referenced and
 *    kept-for-future-reference pictures should not exceed the DPB size.
 *
 *  The RPS Ctor code follows these construction steps.
 ******************************************************************************************/
static EbErrorType PredictionStructureCtor(
    PredictionStructure              *predictionStructurePtr,
    const PredictionStructureConfig  *predictionStructureConfigPtr,
    EbPred                       predType,
    uint32_t                        numberOfReferences)
{
    uint32_t                  entryIndex;
    uint32_t                  configEntryIndex;
    uint32_t                  refIndex;

    // Section Variables
    uint32_t                  leadingPicCount;
    uint32_t                  initPicCount;
    uint32_t                  steadyStatePicCount;

    predictionStructurePtr->dctor = PredictionStructureDctor;

    predictionStructurePtr->pred_type = predType;

    // Set the Pred Struct Period
    predictionStructurePtr->pred_struct_period = predictionStructureConfigPtr->entry_count;

    //----------------------------------------
    // Find the Pred Struct Entry Count
    //   There are four sections of the pred struct:
    //     -Leading Pictures        Size: N-1
    //     -Init Pictures           Size: Ceil(MaxReference, N) - N + 1
    //     -Steady-state Pictures   Size: N
    //     -Trailing Pictures       Size: N-1
    //----------------------------------------

    //----------------------------------------
    // Determine the Prediction Structure Size
    //   First, start by determining
    //   Ceil(MaxReference, N)
    //----------------------------------------
    {
        int32_t maxRef = MIN_SIGNED_VALUE;
        for (configEntryIndex = 0, entryIndex = predictionStructureConfigPtr->entry_count - 1; configEntryIndex < predictionStructureConfigPtr->entry_count; ++configEntryIndex) {
            // Increment through Reference List 0
            refIndex = 0;
            while (refIndex < numberOfReferences && predictionStructureConfigPtr->entry_array[configEntryIndex].ref_list0[refIndex] != 0) {
                //maxRef = MAX(predictionStructureConfigPtr->entry_array[configEntryIndex].ref_list0[refIndex], maxRef);
                maxRef = MAX((int32_t)(predictionStructureConfigPtr->entry_count - entryIndex - 1) + predictionStructureConfigPtr->entry_array[configEntryIndex].ref_list0[refIndex], maxRef);
                ++refIndex;
            }

            // Increment through Reference List 1 (Random Access only)
            if (predType == EB_PRED_RANDOM_ACCESS) {
                refIndex = 0;
                while (refIndex < numberOfReferences && predictionStructureConfigPtr->entry_array[configEntryIndex].ref_list1[refIndex] != 0) {
                    //maxRef = MAX(predictionStructureConfigPtr->entry_array[configEntryIndex].ref_list1[refIndex], maxRef);
                    maxRef = MAX((int32_t)(predictionStructureConfigPtr->entry_count - entryIndex - 1) + predictionStructureConfigPtr->entry_array[configEntryIndex].ref_list1[refIndex], maxRef);
                    ++refIndex;
                }
            }

            // Increment entryIndex
            entryIndex = (entryIndex == predictionStructureConfigPtr->entry_count - 1) ? 0 : entryIndex + 1;
        }

        // Perform the Ceil(MaxReference, N) operation
        predictionStructurePtr->maximum_extent = CEILING(maxRef, predictionStructurePtr->pred_struct_period);

        // Set the Section Sizes
        leadingPicCount = (predType == EB_PRED_RANDOM_ACCESS) ?     // No leading pictures in low-delay configurations
            predictionStructurePtr->pred_struct_period - 1 :
            0;
        initPicCount = predictionStructurePtr->maximum_extent - predictionStructurePtr->pred_struct_period + 1;
        steadyStatePicCount = predictionStructurePtr->pred_struct_period;
        //trailingPicCount        = (predType == EB_PRED_RANDOM_ACCESS) ?     // No trailing pictures in low-delay configurations
        //    predictionStructurePtr->pred_struct_period - 1:
        //    0;

        // Set the total Entry Count
        predictionStructurePtr->pred_struct_entry_count =
            leadingPicCount +
            initPicCount +
            steadyStatePicCount;

        // Set the Section Indices
        predictionStructurePtr->leading_pic_index = 0;
        predictionStructurePtr->init_pic_index = predictionStructurePtr->leading_pic_index + leadingPicCount;
        predictionStructurePtr->steady_state_index = predictionStructurePtr->init_pic_index + initPicCount;
    }

    // Allocate the entry array
    EB_CALLOC_2D(predictionStructurePtr->pred_struct_entry_ptr_array, predictionStructurePtr->pred_struct_entry_count, 1);

    // Find the Max Temporal Layer Index
    predictionStructurePtr->temporal_layer_count = 0;
    for (configEntryIndex = 0; configEntryIndex < predictionStructureConfigPtr->entry_count; ++configEntryIndex)
        predictionStructurePtr->temporal_layer_count = MAX(predictionStructureConfigPtr->entry_array[configEntryIndex].temporal_layer_index, predictionStructurePtr->temporal_layer_count);
    // Increment the Zero-indexed temporal layer index to get the total count
    ++predictionStructurePtr->temporal_layer_count;

    //----------------------------------------
    // Construct Leading Pictures
    //   -Use only Ref List1 from the Config
    //   -Note the Config starts from the 2nd position to construct the leading pictures
    //----------------------------------------
    {
        for (entryIndex = 0, configEntryIndex = 1; entryIndex < leadingPicCount; ++entryIndex, ++configEntryIndex) {
            // Find the Size of the Config's Reference List 1
            refIndex = 0;
            while (refIndex < numberOfReferences && predictionStructureConfigPtr->entry_array[configEntryIndex].ref_list1[refIndex] != 0)
                ++refIndex;
            // Set Leading Picture's Reference List 0 Count {Config List1 => LeadingPic List 0}
            predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list_count = refIndex;

            // Allocate the Leading Picture Reference List 0
            if (predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list_count) {
                EB_MALLOC_ARRAY(predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list, predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list_count);
            }
            else
                predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list = (int32_t*)EB_NULL;
            // Copy Config List1 => LeadingPic Reference List 0
            for (refIndex = 0; refIndex < predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list_count; ++refIndex)
                predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list[refIndex] = predictionStructureConfigPtr->entry_array[configEntryIndex].ref_list1[refIndex];
            // Null out List 1
            predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list_count = 0;
            predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list = (int32_t*)EB_NULL;

            // Set the Temporal Layer Index
            predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->temporal_layer_index = predictionStructureConfigPtr->entry_array[configEntryIndex].temporal_layer_index;

            // Set the Decode Order
            predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->decode_order = (predType == EB_PRED_RANDOM_ACCESS) ?
                predictionStructureConfigPtr->entry_array[configEntryIndex].decode_order :
                entryIndex;
        }
    }

    //----------------------------------------
    // Construct Init Pictures
    //   -Use only references from Ref List0 & Ref List1 from the Config that don't violate CRA mechanics
    //   -The Config Index cycles through continuously
    //----------------------------------------
    {
        uint32_t terminatingEntryIndex = entryIndex + initPicCount;
        int32_t pocValue;

        for (configEntryIndex = 0, pocValue = 0; entryIndex < terminatingEntryIndex; ++entryIndex, ++pocValue) {
            // REFERENCE LIST 0

            // Find the Size of the Config's Reference List 0
            refIndex = 0;
            while (
                refIndex < numberOfReferences &&
                predictionStructureConfigPtr->entry_array[configEntryIndex].ref_list0[refIndex] != 0 &&
                pocValue - predictionStructureConfigPtr->entry_array[configEntryIndex].ref_list0[refIndex] >= 0)  // Stop when we violate the CRA (i.e. reference past it)
            {
                ++refIndex;
            }

            // Set Leading Picture's Reference List 0 Count
            predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list_count = refIndex;

            // Allocate the Leading Picture Reference List 0
            if (predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list_count) {
                EB_MALLOC_ARRAY(predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list, predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list_count);
            }
            else
                predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list = (int32_t*)EB_NULL;
            // Copy Reference List 0
            for (refIndex = 0; refIndex < predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list_count; ++refIndex)
                predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list[refIndex] = predictionStructureConfigPtr->entry_array[configEntryIndex].ref_list0[refIndex];
            // REFERENCE LIST 1
            switch (predType) {
            case EB_PRED_LOW_DELAY_P:

                // Null out List 1
                predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list_count = 0;
                predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list = (int32_t*)EB_NULL;

                break;

            case EB_PRED_LOW_DELAY_B:

                // Copy List 0 => List 1

                // Set Leading Picture's Reference List 1 Count
                predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list_count = predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list_count;

                // Allocate the Leading Picture Reference List 1
                if (predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list_count) {
                    EB_MALLOC_ARRAY(predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list, predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list_count);
                }
                else
                    predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list = (int32_t*)EB_NULL;
                // Copy Reference List 1
                for (refIndex = 0; refIndex < predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list_count; ++refIndex)
                    predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list[refIndex] = predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list[refIndex];
                break;

            case EB_PRED_RANDOM_ACCESS:

                // Find the Size of the Config's Reference List 1
                refIndex = 0;
                while (
                    refIndex < numberOfReferences &&
                    predictionStructureConfigPtr->entry_array[configEntryIndex].ref_list1[refIndex] != 0 &&
                    pocValue - predictionStructureConfigPtr->entry_array[configEntryIndex].ref_list1[refIndex] >= 0) // Stop when we violate the CRA (i.e. reference past it)
                {
                    ++refIndex;
                }

                // Set Leading Picture's Reference List 1 Count
                predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list_count = refIndex;

                // Allocate the Leading Picture Reference List 1
                if (predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list_count) {
                    EB_MALLOC_ARRAY(predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list, predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list_count);
                }
                else
                    predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list = (int32_t*)EB_NULL;
                // Copy Reference List 1
                for (refIndex = 0; refIndex < predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list_count; ++refIndex)
                    predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list[refIndex] = predictionStructureConfigPtr->entry_array[configEntryIndex].ref_list1[refIndex];
                break;

            default:
                break;
            }

            // Set the Temporal Layer Index
            predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->temporal_layer_index = predictionStructureConfigPtr->entry_array[configEntryIndex].temporal_layer_index;

            // Set the Decode Order
            predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->decode_order = (predType == EB_PRED_RANDOM_ACCESS) ?
                predictionStructureConfigPtr->entry_array[configEntryIndex].decode_order :
                entryIndex;

            // Rollover the Config Index
            configEntryIndex = (configEntryIndex == predictionStructureConfigPtr->entry_count - 1) ?
                0 :
                configEntryIndex + 1;
        }
    }

    //----------------------------------------
    // Construct Steady-state Pictures
    //   -Copy directly from the Config
    //----------------------------------------
    {
        uint32_t terminatingEntryIndex = entryIndex + steadyStatePicCount;

        for (/*configEntryIndex = 0*/; entryIndex < terminatingEntryIndex; ++entryIndex/*, ++configEntryIndex*/) {
            // Find the Size of Reference List 0
            refIndex = 0;
            while (refIndex < numberOfReferences && predictionStructureConfigPtr->entry_array[configEntryIndex].ref_list0[refIndex] != 0)
                ++refIndex;
            // Set Reference List 0 Count
            predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list_count = refIndex;

            // Allocate Reference List 0
            EB_MALLOC_ARRAY(predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list, predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list_count);
            // Copy Reference List 0
            for (refIndex = 0; refIndex < predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list_count; ++refIndex)
                predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list[refIndex] = predictionStructureConfigPtr->entry_array[configEntryIndex].ref_list0[refIndex];
            // REFERENCE LIST 1
            switch (predType) {
            case EB_PRED_LOW_DELAY_P:

                // Null out List 1
                predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list_count = 0;
                predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list = (int32_t*)EB_NULL;

                break;

            case EB_PRED_LOW_DELAY_B:

                // Copy List 0 => List 1

                // Set Leading Picture's Reference List 1 Count
                predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list_count = predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list_count;

                // Allocate the Leading Picture Reference List 1
                if (predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list_count) {
                    EB_MALLOC_ARRAY(predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list,  predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list_count);
                }
                else
                    predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list = (int32_t*)EB_NULL;
                // Copy Reference List 1
                for (refIndex = 0; refIndex < predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list_count; ++refIndex)
                    predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list[refIndex] = predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list[refIndex];
                break;

            case EB_PRED_RANDOM_ACCESS:

                // Find the Size of the Config's Reference List 1
                refIndex = 0;
                while (refIndex < numberOfReferences && predictionStructureConfigPtr->entry_array[configEntryIndex].ref_list1[refIndex] != 0)
                    ++refIndex;
                // Set Leading Picture's Reference List 1 Count
                predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list_count = refIndex;

                // Allocate the Leading Picture Reference List 1
                if (predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list_count) {
                    EB_MALLOC_ARRAY(predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list, predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list_count);
                }
                else
                    predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list = (int32_t*)EB_NULL;
                // Copy Reference List 1
                for (refIndex = 0; refIndex < predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list_count; ++refIndex)
                    predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list[refIndex] = predictionStructureConfigPtr->entry_array[configEntryIndex].ref_list1[refIndex];
                break;

            default:
                break;
            }

            // Set the Temporal Layer Index
            predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->temporal_layer_index = predictionStructureConfigPtr->entry_array[configEntryIndex].temporal_layer_index;

            // Set the Decode Order
            predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->decode_order = (predType == EB_PRED_RANDOM_ACCESS) ?
                predictionStructureConfigPtr->entry_array[configEntryIndex].decode_order :
                entryIndex;

            // Rollover the Config Index
            configEntryIndex = (configEntryIndex == predictionStructureConfigPtr->entry_count - 1) ?
                0 :
                configEntryIndex + 1;
        }
    }

    //----------------------------------------
    // Construct Trailing Pictures
    //   -Use only Ref List0 from the Config
    //----------------------------------------
    //{
    //    uint32_t terminatingEntryIndex = entryIndex + trailingPicCount;
    //
    //    for(configEntryIndex = 0; entryIndex < terminatingEntryIndex; ++entryIndex, ++configEntryIndex) {
    //
    //        // Set Reference List 0 Count
    //        // *Note - only 1 reference is used for trailing pictures.  If you have frequent CRAs and the Pred Struct
    //        //   has many references, you can run into edge conditions.
    //        predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list_count = 1;
    //
    //        // Allocate Reference List 0
    //        predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list = (int32_t*) malloc(sizeof(int32_t) * predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list_count);
    //
    //        // Copy Reference List 0
    //        for(refIndex = 0; refIndex < predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list_count; ++refIndex) {
    //            predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list[refIndex] = predictionStructureConfigPtr->entry_array[configEntryIndex].ref_list0[refIndex];
    //        }
    //
    //        // Null out List 1
    //        predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list_count = 0;
    //        predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list = (int32_t*) EB_NULL;
    //
    //        // Set the Temporal Layer Index
    //        predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->temporal_layer_index = predictionStructureConfigPtr->entry_array[configEntryIndex].temporal_layer_index;
    //
    //        // Set the Decode Order
    //        predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->decode_order = (predType == EB_PRED_RANDOM_ACCESS) ?
    //            predictionStructureConfigPtr->entry_array[configEntryIndex].decode_order :
    //            entryIndex;
    //    }
    //}

    //----------------------------------------
    // CONSTRUCT DEPENDENT LIST 0
    //----------------------------------------

    {
        int64_t  depIndex;
        uint64_t  pictureNumber;

        // First, determine the Dependent List Size for each Entry by incrementing the dependent list length
        {
            // Go through a single pass of the Leading Pictures and Init pictures
            for (pictureNumber = 0, entryIndex = 0; pictureNumber < predictionStructurePtr->steady_state_index; ++pictureNumber) {
                // Go through each Reference picture and accumulate counts
                for (refIndex = 0; refIndex < predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list_count; ++refIndex) {
                    depIndex = pictureNumber - predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list[refIndex];

                    if (depIndex >= 0 && depIndex < (int32_t)(predictionStructurePtr->steady_state_index + predictionStructurePtr->pred_struct_period))
                        ++predictionStructurePtr->pred_struct_entry_ptr_array[depIndex]->dep_list0.list_count;
                }

                // Increment the entryIndex
                ++entryIndex;
            }

            // Go through an entire maximum extent pass for the Steady-state pictures
            for (entryIndex = predictionStructurePtr->steady_state_index; pictureNumber <= predictionStructurePtr->steady_state_index + 2 * predictionStructurePtr->maximum_extent; ++pictureNumber) {
                // Go through each Reference picture and accumulate counts
                for (refIndex = 0; refIndex < predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list_count; ++refIndex) {
                    depIndex = pictureNumber - predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list[refIndex];

                    if (depIndex >= 0 && depIndex < (int32_t)(predictionStructurePtr->steady_state_index + predictionStructurePtr->pred_struct_period))
                        ++predictionStructurePtr->pred_struct_entry_ptr_array[depIndex]->dep_list0.list_count;
                }

                // Rollover the entryIndex each time it reaches the end of the steady state index
                entryIndex = (entryIndex == predictionStructurePtr->pred_struct_entry_count - 1) ?
                    predictionStructurePtr->steady_state_index :
                    entryIndex + 1;
            }
        }

        // Second, allocate memory for each dependent list of each Entry
        for (entryIndex = 0; entryIndex < predictionStructurePtr->pred_struct_entry_count; ++entryIndex) {
            // If the dependent list count is non-zero, allocate the list, else the list is NULL.
            if (predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->dep_list0.list_count > 0) {
                EB_MALLOC_ARRAY(predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->dep_list0.list, predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->dep_list0.list_count);
            }
            else
                predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->dep_list0.list = (int32_t*)EB_NULL;
        }

        // Third, reset the Dependent List Length (they are re-derived)
        for (entryIndex = 0; entryIndex < predictionStructurePtr->pred_struct_entry_count; ++entryIndex)
            predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->dep_list0.list_count = 0;
        // Fourth, run through each Reference List entry again and populate the Dependent Lists and Dep List Counts
        {
            // Go through a single pass of the Leading Pictures and Init pictures
            for (pictureNumber = 0, entryIndex = 0; pictureNumber < predictionStructurePtr->steady_state_index; ++pictureNumber) {
                // Go through each Reference picture and accumulate counts
                for (refIndex = 0; refIndex < predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list_count; ++refIndex) {
                    depIndex = pictureNumber - predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list[refIndex];

                    if (depIndex >= 0 && depIndex < (int32_t)(predictionStructurePtr->steady_state_index + predictionStructurePtr->pred_struct_period)) {
                        predictionStructurePtr->pred_struct_entry_ptr_array[depIndex]->dep_list0.list[predictionStructurePtr->pred_struct_entry_ptr_array[depIndex]->dep_list0.list_count++] =
                            predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list[refIndex];
                    }
                }

                // Increment the entryIndex
                ++entryIndex;
            }

            // Go through an entire maximum extent pass for the Steady-state pictures
            for (entryIndex = predictionStructurePtr->steady_state_index; pictureNumber <= predictionStructurePtr->steady_state_index + 2 * predictionStructurePtr->maximum_extent; ++pictureNumber) {
                // Go through each Reference picture and accumulate counts
                for (refIndex = 0; refIndex < predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list_count; ++refIndex) {
                    depIndex = pictureNumber - predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list[refIndex];

                    // Assign the Reference to the Dep List and Increment the Dep List Count
                    if (depIndex >= 0 && depIndex < (int32_t)(predictionStructurePtr->steady_state_index + predictionStructurePtr->pred_struct_period)) {
                        predictionStructurePtr->pred_struct_entry_ptr_array[depIndex]->dep_list0.list[predictionStructurePtr->pred_struct_entry_ptr_array[depIndex]->dep_list0.list_count++] =
                            predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list[refIndex];
                    }
                }

                // Rollover the entryIndex each time it reaches the end of the steady state index
                entryIndex = (entryIndex == predictionStructurePtr->pred_struct_entry_count - 1) ?
                    predictionStructurePtr->steady_state_index :
                    entryIndex + 1;
            }
        }
    }

    //----------------------------------------
    // CONSTRUCT DEPENDENT LIST 1
    //----------------------------------------

    {
        int32_t  depIndex;
        uint32_t  pictureNumber;

        // First, determine the Dependent List Size for each Entry by incrementing the dependent list length
        {
            // Go through a single pass of the Leading Pictures and Init pictures
            for (pictureNumber = 0, entryIndex = 0; pictureNumber < predictionStructurePtr->steady_state_index; ++pictureNumber) {
                // Go through each Reference picture and accumulate counts
                for (refIndex = 0; refIndex < predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list_count; ++refIndex) {
                    depIndex = pictureNumber - predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list[refIndex];

                    if (depIndex >= 0 && depIndex < (int32_t)(predictionStructurePtr->steady_state_index + predictionStructurePtr->pred_struct_period))
                        ++predictionStructurePtr->pred_struct_entry_ptr_array[depIndex]->dep_list1.list_count;
                }

                // Increment the entryIndex
                ++entryIndex;
            }

            // Go through an entire maximum extent pass for the Steady-state pictures
            for (entryIndex = predictionStructurePtr->steady_state_index; pictureNumber <= predictionStructurePtr->steady_state_index + 2 * predictionStructurePtr->maximum_extent; ++pictureNumber) {
                // Go through each Reference picture and accumulate counts
                for (refIndex = 0; refIndex < predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list_count; ++refIndex) {
                    depIndex = pictureNumber - predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list[refIndex];

                    if (depIndex >= 0 && depIndex < (int32_t)(predictionStructurePtr->steady_state_index + predictionStructurePtr->pred_struct_period))
                        ++predictionStructurePtr->pred_struct_entry_ptr_array[depIndex]->dep_list1.list_count;
                }

                // Rollover the entryIndex each time it reaches the end of the steady state index
                entryIndex = (entryIndex == predictionStructurePtr->pred_struct_entry_count - 1) ?
                    predictionStructurePtr->steady_state_index :
                    entryIndex + 1;
            }
        }

        // Second, allocate memory for each dependent list of each Entry
        for (entryIndex = 0; entryIndex < predictionStructurePtr->pred_struct_entry_count; ++entryIndex) {
            // If the dependent list count is non-zero, allocate the list, else the list is NULL.
            if (predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->dep_list1.list_count > 0) {
                EB_MALLOC_ARRAY(predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->dep_list1.list, predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->dep_list1.list_count);
            }
            else
                predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->dep_list1.list = (int32_t*)EB_NULL;
        }

        // Third, reset the Dependent List Length (they are re-derived)
        for (entryIndex = 0; entryIndex < predictionStructurePtr->pred_struct_entry_count; ++entryIndex)
            predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->dep_list1.list_count = 0;
        // Fourth, run through each Reference List entry again and populate the Dependent Lists and Dep List Counts
        {
            // Go through a single pass of the Leading Pictures and Init pictures
            for (pictureNumber = 0, entryIndex = 0; pictureNumber < predictionStructurePtr->steady_state_index; ++pictureNumber) {
                // Go through each Reference picture and accumulate counts
                for (refIndex = 0; refIndex < predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list_count; ++refIndex) {
                    depIndex = pictureNumber - predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list[refIndex];

                    if (depIndex >= 0 && depIndex < (int32_t)(predictionStructurePtr->steady_state_index + predictionStructurePtr->pred_struct_period)) {
                        predictionStructurePtr->pred_struct_entry_ptr_array[depIndex]->dep_list1.list[predictionStructurePtr->pred_struct_entry_ptr_array[depIndex]->dep_list1.list_count++] =
                            predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list[refIndex];
                    }
                }

                // Increment the entryIndex
                ++entryIndex;
            }

            // Go through an entire maximum extent pass for the Steady-state pictures
            for (entryIndex = predictionStructurePtr->steady_state_index; pictureNumber <= predictionStructurePtr->steady_state_index + 2 * predictionStructurePtr->maximum_extent; ++pictureNumber) {
                // Go through each Reference picture and accumulate counts
                for (refIndex = 0; refIndex < predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list_count; ++refIndex) {
                    depIndex = pictureNumber - predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list[refIndex];

                    // Assign the Reference to the Dep List and Increment the Dep List Count
                    if (depIndex >= 0 && depIndex < (int32_t)(predictionStructurePtr->steady_state_index + predictionStructurePtr->pred_struct_period)) {
                        predictionStructurePtr->pred_struct_entry_ptr_array[depIndex]->dep_list1.list[predictionStructurePtr->pred_struct_entry_ptr_array[depIndex]->dep_list1.list_count++] =
                            predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list[refIndex];
                    }
                }

                // Rollover the entryIndex each time it reaches the end of the steady state index
                entryIndex = (entryIndex == predictionStructurePtr->pred_struct_entry_count - 1) ?
                    predictionStructurePtr->steady_state_index :
                    entryIndex + 1;
            }
        }
    }

    // Set is_referenced for each entry
    for (entryIndex = 0; entryIndex < predictionStructurePtr->pred_struct_entry_count; ++entryIndex) {
        predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->is_referenced =
            ((predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->dep_list0.list_count > 0) ||
            (predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->dep_list1.list_count > 0)) ?
            EB_TRUE :
            EB_FALSE;
    }

    //----------------------------------------
    // CONSTRUCT THE RPSes
    //----------------------------------------
    {
        // Counts & Indicies
        uint32_t      refIndex;
        uint32_t      depIndex;
        uint32_t      entryIndex;
        uint32_t      currentPocIndex;
        uint32_t      refPocIndex;

        uint32_t      decodeOrderTableSize;
        int32_t     *decodeOrderTable;
        uint32_t     *displayOrderTable;
        uint32_t      gopNumber;
        uint32_t      baseNumber;

        // Timeline Map Variables
        EbBool    *timelineMap;
        uint32_t      timelineSize;

        int32_t      depListMax;
        int32_t      depListMin;

        int32_t      lifetimeStart;
        int32_t      lifetimeSpan;

        int32_t      deltaPoc;
        int32_t      prevDeltaPoc;
        EbBool     pocInreference_list0;
        EbBool     pocInreference_list1;
        EbBool     pocInTimeline;

        int32_t      adjustedDepIndex;

        // Allocate & Initialize the Timeline map
        timelineSize = predictionStructurePtr->pred_struct_entry_count;
        decodeOrderTableSize = CEILING(predictionStructurePtr->pred_struct_entry_count + predictionStructurePtr->maximum_extent, predictionStructurePtr->pred_struct_entry_count);
        EB_CALLOC_ARRAY(predictionStructurePtr->timelineMap, SQR(timelineSize));
        timelineMap = predictionStructurePtr->timelineMap;

        // Construct the Decode & Display Order
        EB_MALLOC_ARRAY(predictionStructurePtr->decodeOrderTable, decodeOrderTableSize);
        decodeOrderTable = predictionStructurePtr->decodeOrderTable;

        EB_MALLOC_ARRAY(predictionStructurePtr->displayOrderTable, decodeOrderTableSize);
        displayOrderTable = predictionStructurePtr->displayOrderTable;

        for (currentPocIndex = 0, entryIndex = 0; currentPocIndex < decodeOrderTableSize; ++currentPocIndex) {
            // Set the Decode Order
            gopNumber = (currentPocIndex / predictionStructurePtr->pred_struct_period);
            baseNumber = gopNumber * predictionStructurePtr->pred_struct_period;

            if (predType == EB_PRED_RANDOM_ACCESS)
                decodeOrderTable[currentPocIndex] = baseNumber + predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->decode_order;
            else
                decodeOrderTable[currentPocIndex] = currentPocIndex;
            displayOrderTable[decodeOrderTable[currentPocIndex]] = currentPocIndex;

            // Increment the entryIndex
            entryIndex = (entryIndex == predictionStructurePtr->pred_struct_entry_count - 1) ?
                predictionStructurePtr->pred_struct_entry_count - predictionStructurePtr->pred_struct_period :
                entryIndex + 1;
        }

        // Construct the timeline map from the dependency lists
        for (refPocIndex = 0, entryIndex = 0; refPocIndex < timelineSize; ++refPocIndex) {
            // Initialize Max to most negative signed value and Min to most positive signed value
            depListMax = MIN_SIGNED_VALUE;
            depListMin = MAX_SIGNED_VALUE;

            // Find depListMax and depListMin for the entryIndex in the prediction structure for dep_list0
            for (depIndex = 0; depIndex < predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->dep_list0.list_count; ++depIndex) {
                adjustedDepIndex = predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->dep_list0.list[depIndex] + (int32_t)refPocIndex;

                //if(adjustedDepIndex >= 0 && adjustedDepIndex < (int32_t) timelineSize) {
                if (adjustedDepIndex >= 0) {
                    // Update Max
                    depListMax = MAX(decodeOrderTable[adjustedDepIndex], depListMax);

                    // Update Min
                    depListMin = MIN(decodeOrderTable[adjustedDepIndex], depListMin);
                }
            }

            // Continue search for depListMax and depListMin for the entryIndex in the prediction structure for dep_list1,
            //   the lists are combined in the RPS logic
            for (depIndex = 0; depIndex < predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->dep_list1.list_count; ++depIndex) {
                adjustedDepIndex = predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->dep_list1.list[depIndex] + (int32_t)refPocIndex;

                //if(adjustedDepIndex >= 0 && adjustedDepIndex < (int32_t) timelineSize)  {
                if (adjustedDepIndex >= 0) {
                    // Update Max
                    depListMax = MAX(decodeOrderTable[adjustedDepIndex], depListMax);

                    // Update Min
                    depListMin = MIN(decodeOrderTable[adjustedDepIndex], depListMin);
                }
            }

            // If the Dependent Lists are empty, ensure that no RPS signaling is set
            if ((predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->dep_list0.list_count > 0) ||
                (predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->dep_list1.list_count > 0)) {
                // Determine lifetimeStart and lifetimeSpan - its important to note that out-of-range references are
                //   clipped/eliminated to not violate IDR/CRA referencing rules
                lifetimeStart = depListMin;

                if (lifetimeStart < (int32_t)timelineSize) {
                    lifetimeStart = CLIP3(0, (int32_t)(timelineSize - 1), lifetimeStart);

                    lifetimeSpan = depListMax - depListMin + 1;
                    lifetimeSpan = CLIP3(0, (int32_t)timelineSize - lifetimeStart, lifetimeSpan);

                    // Set the timelineMap
                    for (currentPocIndex = (uint32_t)lifetimeStart; currentPocIndex < (uint32_t)(lifetimeStart + lifetimeSpan); ++currentPocIndex)
                        timelineMap[refPocIndex*timelineSize + displayOrderTable[currentPocIndex]] = EB_TRUE;
                }
            }

            // Increment the entryIndex
            entryIndex = (entryIndex == predictionStructurePtr->pred_struct_entry_count - 1) ?
                predictionStructurePtr->pred_struct_entry_count - predictionStructurePtr->pred_struct_period :
                entryIndex + 1;
        }

        //--------------------------------------------------------
        // Create the RPS for Prediction Structure Entry
        //--------------------------------------------------------

        // *Note- many of the below Syntax Elements are signaled
        //    in the Slice Header and not in the RPS.  These syntax
        //    elements can be configured during runtime to manipulate
        //    existing RPS structures.  E.g. a reference list could
        //    be shortened...

        // Initialize the RPS Group
        predictionStructurePtr->restricted_ref_pic_lists_enable_flag = EB_TRUE;
        predictionStructurePtr->lists_modification_enable_flag = EB_FALSE;
        predictionStructurePtr->long_term_enable_flag = EB_FALSE;
        predictionStructurePtr->default_ref_pics_list0_total_count_minus1 = 0;
        predictionStructurePtr->default_ref_pics_list1_total_count_minus1 = 0;

        // For each RPS Index
        for (entryIndex = 0; entryIndex < predictionStructurePtr->pred_struct_entry_count; ++entryIndex) {
            // Determine the Current POC Index
            currentPocIndex = entryIndex;

            // Initialize the RPS
            predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->short_term_rps_in_sps_flag = EB_TRUE;
            predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->short_term_rps_in_sps_index = entryIndex;
            predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->inter_rps_prediction_flag = EB_FALSE;
            predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->long_term_rps_present_flag = EB_FALSE;
            predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->list0_modification_flag = EB_FALSE;
            predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->negative_ref_pics_total_count = 0;
            predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->positive_ref_pics_total_count = 0;
            predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_pics_list0_total_count_minus1 = ~0;
            predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_pics_list1_total_count_minus1 = ~0;

            // Create the Negative List
            prevDeltaPoc = 0;
            for (refPocIndex = currentPocIndex - 1; (int32_t)refPocIndex >= 0; --refPocIndex) {
                // Find the deltaPoc value
                deltaPoc = (int32_t)currentPocIndex - (int32_t)refPocIndex;

                // Check to see if the deltaPoc is in Reference List 0
                pocInreference_list0 = EB_FALSE;
                for (refIndex = 0; (refIndex < predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list_count) && (pocInreference_list0 == EB_FALSE); ++refIndex) {
                    // Reference List 0
                    if ((predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list[refIndex] != 0) &&
                        (predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list[refIndex] == deltaPoc))
                    {
                        pocInreference_list0 = EB_TRUE;
                        ++predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_pics_list0_total_count_minus1;
                    }
                }

                // Check to see if the deltaPoc is in Reference List 1
                pocInreference_list1 = EB_FALSE;
                for (refIndex = 0; (refIndex < predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list_count) && (pocInreference_list1 == EB_FALSE); ++refIndex) {
                    // Reference List 1
                    if ((predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list[refIndex] != 0) &&
                        (predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list[refIndex] == deltaPoc))
                    {
                        pocInreference_list1 = EB_TRUE;
                        ++predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_pics_list1_total_count_minus1;
                    }
                }

                // Check to see if the refPocIndex is in the timeline
                pocInTimeline = timelineMap[refPocIndex*timelineSize + currentPocIndex];

                // If the deltaPoc is in the timeline
                if (pocInTimeline == EB_TRUE) {
                    predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->used_by_negative_curr_pic_flag[predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->negative_ref_pics_total_count] = (pocInreference_list0 == EB_TRUE || pocInreference_list1 == EB_TRUE) ? EB_TRUE : EB_FALSE;
                    predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->delta_negative_gop_pos_minus1[predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->negative_ref_pics_total_count++] = deltaPoc - 1 - prevDeltaPoc;
                    prevDeltaPoc = deltaPoc;
                }
            }

            // Create the Positive List
            prevDeltaPoc = 0;
            for (refPocIndex = currentPocIndex + 1; refPocIndex < timelineSize; ++refPocIndex) {
                // Find the deltaPoc value
                deltaPoc = (int32_t)currentPocIndex - (int32_t)refPocIndex;

                // Check to see if the deltaPoc is in Ref List 0
                pocInreference_list0 = EB_FALSE;
                for (refIndex = 0; (refIndex < predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list_count) && (pocInreference_list0 == EB_FALSE); ++refIndex) {
                    // Reference List 0
                    if ((predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list[refIndex] != 0) &&
                        (predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list0.reference_list[refIndex] == deltaPoc))
                    {
                        pocInreference_list0 = EB_TRUE;
                        ++predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_pics_list0_total_count_minus1;
                    }
                }

                // Check to see if the deltaPoc is in Ref List 1
                pocInreference_list1 = EB_FALSE;
                for (refIndex = 0; (refIndex < predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list_count) && (pocInreference_list1 == EB_FALSE); ++refIndex) {
                    if ((predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list[refIndex] != 0) &&
                        (predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_list1.reference_list[refIndex] == deltaPoc))
                    {
                        pocInreference_list1 = EB_TRUE;
                        ++predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_pics_list1_total_count_minus1;
                    }
                }

                // Check to see if the Y-position is in the timeline
                pocInTimeline = timelineMap[refPocIndex*timelineSize + currentPocIndex];

                // If the Y-position is in the time lime
                if (pocInTimeline == EB_TRUE) {
                    predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->used_by_positive_curr_pic_flag[predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->positive_ref_pics_total_count] = (pocInreference_list0 == EB_TRUE || pocInreference_list1 == EB_TRUE) ? EB_TRUE : EB_FALSE;
                    predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->delta_positive_gop_pos_minus1[predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->positive_ref_pics_total_count++] = -deltaPoc - 1 - prevDeltaPoc;
                    prevDeltaPoc = -deltaPoc;
                }
            }

            // Adjust Reference Counts if list is empty
            predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_pics_list0_total_count_minus1 =
                (predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_pics_list0_total_count_minus1 == ~0) ?
                0 :
                predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_pics_list0_total_count_minus1;

            predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_pics_list1_total_count_minus1 =
                (predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_pics_list1_total_count_minus1 == ~0) ?
                0 :
                predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_pics_list1_total_count_minus1;

            // Set ref_pics_override_total_count_flag to TRUE if Reflist_count is different than the default
            predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_pics_override_total_count_flag =
                ((predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_pics_list0_total_count_minus1 != (int32_t)predictionStructurePtr->default_ref_pics_list0_total_count_minus1) ||
                (predictionStructurePtr->pred_struct_entry_ptr_array[entryIndex]->ref_pics_list1_total_count_minus1 != (int32_t)predictionStructurePtr->default_ref_pics_list1_total_count_minus1)) ?
                EB_TRUE :
                EB_FALSE;
        }
    }

    return EB_ErrorNone;
}

static void prediction_structure_group_dctor(EbPtr p)
{
    PredictionStructureGroup *obj = (PredictionStructureGroup*)p;
    EB_DELETE_PTR_ARRAY(obj->prediction_structure_ptr_array, obj->prediction_structure_count);
}
/*************************************************
 * Prediction Structure Group Ctor
 *
 * Summary: Converts the Prediction Structure Config
 *   into the usable Prediction Structure with RPS and
 *   Dependent List control.
 *
 * From each config, several prediction structures
 *   are created. These include:
 *   -Variable Number of References
 *      # [1 - 4]
 *   -Temporal Layers
 *      # [1 - 6]
 *   -GOP Type
 *      # Low Delay P
 *      # Low Delay B
 *      # Random Access
 *
 *************************************************/

EbErrorType prediction_structure_group_ctor(
    PredictionStructureGroup   *predictionStructureGroupPtr,
    uint8_t          enc_mode,
    uint32_t                         baseLayerSwitchMode)
{
    uint32_t          predStructIndex = 0;
    uint32_t          refIdx;
    uint32_t          hierarchicalLevelIdx;
    uint32_t          predTypeIdx;
    uint32_t          numberOfReferences;

    predictionStructureGroupPtr->dctor = prediction_structure_group_dctor;

    if (enc_mode > ENC_M0) {
        for (int gop_i = 1; gop_i < 8; ++gop_i) {
            for (int i = 1; i < 4; ++i) {
                four_level_hierarchical_pred_struct[gop_i].ref_list0[i] = 0;
                four_level_hierarchical_pred_struct[gop_i].ref_list1[i] = 0;
            }
        }

        for (int gop_i = 1; gop_i < 16; ++gop_i) {
            for (int i = 1; i < 4; ++i) {
                five_level_hierarchical_pred_struct[gop_i].ref_list0[i] = 0;
                five_level_hierarchical_pred_struct[gop_i].ref_list1[i] = 0;
            }
        }
    }

    // Count the number of Prediction Structures
    while ((PredictionStructureConfigArray[predStructIndex].entry_array != 0) && (PredictionStructureConfigArray[predStructIndex].entry_count != 0)) {
        // Get Random Access + P for temporal ID 0
        if (PredictionStructureConfigArray[predStructIndex].entry_array->temporal_layer_index == 0 && baseLayerSwitchMode) {
            for (refIdx = 0; refIdx < REF_LIST_MAX_DEPTH; ++refIdx)
                PredictionStructureConfigArray[predStructIndex].entry_array->ref_list1[refIdx] = 0;
        }
        ++predStructIndex;
    }

    predictionStructureGroupPtr->prediction_structure_count = MAX_TEMPORAL_LAYERS * EB_PRED_TOTAL_COUNT * REF_LIST_MAX_DEPTH;
    EB_ALLOC_PTR_ARRAY(predictionStructureGroupPtr->prediction_structure_ptr_array, predictionStructureGroupPtr->prediction_structure_count);
    for (hierarchicalLevelIdx = 0; hierarchicalLevelIdx < MAX_TEMPORAL_LAYERS; ++hierarchicalLevelIdx) {
        for (predTypeIdx = 0; predTypeIdx < EB_PRED_TOTAL_COUNT; ++predTypeIdx) {
            for (refIdx = 0; refIdx < REF_LIST_MAX_DEPTH; ++refIdx) {
                predStructIndex = PRED_STRUCT_INDEX(hierarchicalLevelIdx, predTypeIdx, refIdx);
                numberOfReferences = refIdx + 1;

                EB_NEW(
                    predictionStructureGroupPtr->prediction_structure_ptr_array[predStructIndex],
                    PredictionStructureCtor,
                    &(PredictionStructureConfigArray[hierarchicalLevelIdx]),
                    (EbPred)predTypeIdx,
                    numberOfReferences);
            }
        }
    }

    return EB_ErrorNone;
}
