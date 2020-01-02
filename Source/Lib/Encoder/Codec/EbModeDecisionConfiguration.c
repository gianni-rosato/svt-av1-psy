// clang-format off
/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbModeDecisionConfiguration.h"
#include "EbRateDistortionCost.h"
#include "EbUtility.h"
#include "EbModeDecisionProcess.h"
#include "EbDefinitions.h"

#include "EbLog.h"
/********************************************
 * Constants
 ********************************************/
int pa_to_ep_block_index[85] = {
    0    ,
    25   ,
    50   ,
    75   ,    84   ,    93   ,    102  ,
    111  ,
    136  ,    145  ,    154  ,    163  ,
    172  ,
    197  ,    206  ,    215  ,    224  ,
    233  ,
    258  ,    267  ,    276  ,    285  ,
    294  ,
    319  ,
    344  ,    353  ,    362  ,    371  ,
    380  ,
    405  ,    414  ,    423  ,    432  ,
    441  ,
    466  ,    475  ,    484  ,   493  ,
    502  ,
    527  ,    536  ,    545  ,    554  ,
    563  ,
    588  ,
    613  ,    622  ,    631  ,    640  ,
    649  ,
    674  ,    683  ,    692  ,    701  ,
    710  ,
    735  ,    744  ,    753  ,    762  ,
    771  ,
    796  ,    805  ,    814  ,    823  ,
    832  ,
    857  ,
    882  ,    891  ,    900  ,    909  ,
    918  ,
    943  ,    952  ,    961  ,    970  ,
    979  ,
    1004 ,    1013 ,    1022 ,    1031 ,
    1040 ,
    1065 ,    1074 ,    1083 ,    1092
};
#define ADD_CU_STOP_SPLIT             0   // Take into account & Stop Splitting
#define ADD_CU_CONTINUE_SPLIT         1   // Take into account & Continue Splitting
#define DO_NOT_ADD_CU_CONTINUE_SPLIT  2   // Do not take into account & Continue Splitting

#define DEPTH_64                      0   // Depth corresponding to the CU size
#define DEPTH_32                      1   // Depth corresponding to the CU size
#define DEPTH_16                      2   // Depth corresponding to the CU size
#define DEPTH_8                       3   // Depth corresponding to the CU size

static const uint8_t parent_blk_index[85] =
{
    0,
    0, 0, 0, 1, 2, 3, 5, 0, 1, 2, 3, 10, 0, 1, 2, 3, 15, 0, 1, 2, 3,
    21, 0, 0, 1, 2, 3, 5, 0, 1, 2, 3, 10, 0, 1, 2, 3, 15, 0, 1, 2, 3,
    42, 0, 0, 1, 2, 3, 5, 0, 1, 2, 3, 10, 0, 1, 2, 3, 15, 0, 1, 2, 3,
    36, 0, 0, 1, 2, 3, 5, 0, 1, 2, 3, 10, 0, 1, 2, 3, 15, 0, 1, 2, 3,
};

const uint8_t incremental_count[85] = {
    //64x64
    0,
    //32x32
    4, 4,
    4, 4,
    //16x16
    0, 0, 0, 0,
    0, 4, 0, 4,
    0, 0, 0, 0,
    0, 4, 0, 4,
    //8x8
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 4, 0, 0, 0, 4,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 4, 0, 0, 0, 4
};


/*******************************************
mdcSetDepth : set depth to be tested
*******************************************/
#define REFINEMENT_P        0x01
#define REFINEMENT_Pp1      0x02
#define REFINEMENT_Pp2      0x04
#define REFINEMENT_Pp3      0x08
#define REFINEMENT_Pm1      0x10
#define REFINEMENT_Pm2      0x20
#define REFINEMENT_Pm3      0x40

EbErrorType mdc_refinement(
    MdcpLocalBlkStruct                   *local_blk_array,
    uint32_t                                  blk_index,
    uint32_t                                  depth,
    uint8_t                                   refinement_level,
    uint8_t                                   lowest_level)
{
    EbErrorType return_error = EB_ErrorNone;

    if (refinement_level & REFINEMENT_P) {
        if (lowest_level == REFINEMENT_P)
            local_blk_array[blk_index].stop_split = EB_TRUE;
    }
    else
        local_blk_array[blk_index].selected_cu = EB_FALSE;
    if (refinement_level & REFINEMENT_Pp1) {
        if (depth < 3 && blk_index < 81) {
            local_blk_array[blk_index + 1].selected_cu = EB_TRUE;
            local_blk_array[blk_index + 1 + depth_offset[depth + 1]].selected_cu = EB_TRUE;
            local_blk_array[blk_index + 1 + 2 * depth_offset[depth + 1]].selected_cu = EB_TRUE;
            local_blk_array[blk_index + 1 + 3 * depth_offset[depth + 1]].selected_cu = EB_TRUE;
        }
        if (lowest_level == REFINEMENT_Pp1) {
            if (depth < 3 && blk_index < 81) {
                local_blk_array[blk_index + 1].stop_split = EB_TRUE;
                local_blk_array[blk_index + 1 + depth_offset[depth + 1]].stop_split = EB_TRUE;
                local_blk_array[blk_index + 1 + 2 * depth_offset[depth + 1]].stop_split = EB_TRUE;
                local_blk_array[blk_index + 1 + 3 * depth_offset[depth + 1]].stop_split = EB_TRUE;
            }
        }
    }

    if (refinement_level & REFINEMENT_Pp2) {
        if (depth < 2 && blk_index < 65) {
            local_blk_array[blk_index + 1 + 1].selected_cu = EB_TRUE;
            local_blk_array[blk_index + 1 + 1 + depth_offset[depth + 2]].selected_cu = EB_TRUE;
            local_blk_array[blk_index + 1 + 1 + 2 * depth_offset[depth + 2]].selected_cu = EB_TRUE;
            local_blk_array[blk_index + 1 + 1 + 3 * depth_offset[depth + 2]].selected_cu = EB_TRUE;

            local_blk_array[blk_index + 1 + depth_offset[depth + 1] + 1].selected_cu = EB_TRUE;
            local_blk_array[blk_index + 1 + depth_offset[depth + 1] + 1 + depth_offset[depth + 2]].selected_cu = EB_TRUE;
            local_blk_array[blk_index + 1 + depth_offset[depth + 1] + 1 + 2 * depth_offset[depth + 2]].selected_cu = EB_TRUE;
            local_blk_array[blk_index + 1 + depth_offset[depth + 1] + 1 + 3 * depth_offset[depth + 2]].selected_cu = EB_TRUE;

            local_blk_array[blk_index + 1 + 2 * depth_offset[depth + 1] + 1].selected_cu = EB_TRUE;
            local_blk_array[blk_index + 1 + 2 * depth_offset[depth + 1] + 1 + depth_offset[depth + 2]].selected_cu = EB_TRUE;
            local_blk_array[blk_index + 1 + 2 * depth_offset[depth + 1] + 1 + 2 * depth_offset[depth + 2]].selected_cu = EB_TRUE;
            local_blk_array[blk_index + 1 + 2 * depth_offset[depth + 1] + 1 + 3 * depth_offset[depth + 2]].selected_cu = EB_TRUE;

            local_blk_array[blk_index + 1 + 3 * depth_offset[depth + 1] + 1].selected_cu = EB_TRUE;
            local_blk_array[blk_index + 1 + 3 * depth_offset[depth + 1] + 1 + depth_offset[depth + 2]].selected_cu = EB_TRUE;
            local_blk_array[blk_index + 1 + 3 * depth_offset[depth + 1] + 1 + 2 * depth_offset[depth + 2]].selected_cu = EB_TRUE;
            local_blk_array[blk_index + 1 + 3 * depth_offset[depth + 1] + 1 + 3 * depth_offset[depth + 2]].selected_cu = EB_TRUE;
        }
        if (lowest_level == REFINEMENT_Pp2) {
            if (depth < 2 && blk_index < 65) {
                local_blk_array[blk_index + 1 + 1].stop_split = EB_TRUE;
                local_blk_array[blk_index + 1 + 1 + depth_offset[depth + 2]].stop_split = EB_TRUE;
                local_blk_array[blk_index + 1 + 1 + 2 * depth_offset[depth + 2]].stop_split = EB_TRUE;
                local_blk_array[blk_index + 1 + 1 + 3 * depth_offset[depth + 2]].stop_split = EB_TRUE;

                local_blk_array[blk_index + 1 + depth_offset[depth + 1] + 1].stop_split = EB_TRUE;
                local_blk_array[blk_index + 1 + depth_offset[depth + 1] + 1 + depth_offset[depth + 2]].stop_split = EB_TRUE;
                local_blk_array[blk_index + 1 + depth_offset[depth + 1] + 1 + 2 * depth_offset[depth + 2]].stop_split = EB_TRUE;
                local_blk_array[blk_index + 1 + depth_offset[depth + 1] + 1 + 3 * depth_offset[depth + 2]].stop_split = EB_TRUE;

                local_blk_array[blk_index + 1 + 2 * depth_offset[depth + 1] + 1].stop_split = EB_TRUE;
                local_blk_array[blk_index + 1 + 2 * depth_offset[depth + 1] + 1 + depth_offset[depth + 2]].stop_split = EB_TRUE;
                local_blk_array[blk_index + 1 + 2 * depth_offset[depth + 1] + 1 + 2 * depth_offset[depth + 2]].stop_split = EB_TRUE;
                local_blk_array[blk_index + 1 + 2 * depth_offset[depth + 1] + 1 + 3 * depth_offset[depth + 2]].stop_split = EB_TRUE;

                local_blk_array[blk_index + 1 + 3 * depth_offset[depth + 1] + 1].stop_split = EB_TRUE;
                local_blk_array[blk_index + 1 + 3 * depth_offset[depth + 1] + 1 + depth_offset[depth + 2]].stop_split = EB_TRUE;
                local_blk_array[blk_index + 1 + 3 * depth_offset[depth + 1] + 1 + 2 * depth_offset[depth + 2]].stop_split = EB_TRUE;
                local_blk_array[blk_index + 1 + 3 * depth_offset[depth + 1] + 1 + 3 * depth_offset[depth + 2]].stop_split = EB_TRUE;
            }
        }
    }

    if (refinement_level & REFINEMENT_Pp3) {
        uint8_t in_loop;
        uint8_t out_loop;
        uint8_t blk_index = 2;
        if (depth == 0) {
            for (out_loop = 0; out_loop < 16; ++out_loop) {
                for (in_loop = 0; in_loop < 4; ++in_loop)
                    local_blk_array[++blk_index].selected_cu = EB_TRUE;
                blk_index += blk_index == 21 ? 2 : blk_index == 42 ? 2 : blk_index == 63 ? 2 : 1;
            }
            if (lowest_level == REFINEMENT_Pp3) {
                blk_index = 2;
                for (out_loop = 0; out_loop < 16; ++out_loop) {
                    for (in_loop = 0; in_loop < 4; ++in_loop)
                        local_blk_array[++blk_index].stop_split = EB_TRUE;
                    blk_index += blk_index == 21 ? 2 : blk_index == 42 ? 2 : blk_index == 63 ? 2 : 1;
                }
            }
        }
    }

    if (refinement_level & REFINEMENT_Pm1) {
        if (depth > 0)
            local_blk_array[blk_index - 1 - parent_blk_index[blk_index]].selected_cu = EB_TRUE;
        if (lowest_level == REFINEMENT_Pm1) {
            if (depth > 0)
                local_blk_array[blk_index - 1 - parent_blk_index[blk_index]].stop_split = EB_TRUE;
        }
    }

    if (refinement_level & REFINEMENT_Pm2) {
        if (depth == 2)
            local_blk_array[0].selected_cu = EB_TRUE;
        if (depth == 3) {
            local_blk_array[1].selected_cu = EB_TRUE;
            local_blk_array[22].selected_cu = EB_TRUE;
            local_blk_array[43].selected_cu = EB_TRUE;
            local_blk_array[64].selected_cu = EB_TRUE;
        }
        if (lowest_level == REFINEMENT_Pm2) {
            if (depth == 2)
                local_blk_array[0].stop_split = EB_TRUE;
            if (depth == 3) {
                local_blk_array[1].stop_split = EB_TRUE;
                local_blk_array[22].stop_split = EB_TRUE;
                local_blk_array[43].stop_split = EB_TRUE;
                local_blk_array[64].stop_split = EB_TRUE;
            }
        }
    }

    if (refinement_level & REFINEMENT_Pm3) {
        if (depth == 3)
            local_blk_array[0].selected_cu = EB_TRUE;
        if (lowest_level == REFINEMENT_Pm2) {
            if (depth == 3)
                local_blk_array[0].stop_split = EB_TRUE;
        }
    }

    return return_error;
}

void refinement_prediction_loop(
    SequenceControlSet                   *scs_ptr,
    PictureControlSet                    *pcs_ptr,
    uint32_t                              sb_index,
    ModeDecisionConfigurationContext     *context_ptr)
{
    MdcpLocalBlkStruct    *local_blk_array         = context_ptr->local_blk_array;
    SbParams               *sb_params            = &scs_ptr->sb_params_array[sb_index];
    uint32_t                  blk_index             = 0;
    while (blk_index < CU_MAX_COUNT)
    {
        if (sb_params->raster_scan_blk_validity[md_scan_to_raster_scan[blk_index]] && (local_blk_array[blk_index].early_split_flag == EB_FALSE))
        {
            local_blk_array[blk_index].selected_cu = EB_TRUE;
            uint32_t depth = get_coded_blk_stats(blk_index)->depth;
            uint8_t refinement_level;
            {
                if (pcs_ptr->parent_pcs_ptr->pic_depth_mode == PIC_SB_SWITCH_DEPTH_MODE && pcs_ptr->parent_pcs_ptr->sb_depth_mode_array[sb_index] == SB_PRED_OPEN_LOOP_DEPTH_MODE)
                    refinement_level = Pred;
                else

                    if (pcs_ptr->parent_pcs_ptr->pic_depth_mode == PIC_SB_SWITCH_DEPTH_MODE && pcs_ptr->parent_pcs_ptr->sb_depth_mode_array[sb_index] == SB_FAST_OPEN_LOOP_DEPTH_MODE)
                        refinement_level = ndp_level_1[depth];
                    else  { // SB_OPEN_LOOP_DEPTH_MODE
                        refinement_level = ndp_level_0[depth];
                    }

                if (pcs_ptr->parent_pcs_ptr->cu8x8_mode == CU_8x8_MODE_1) {
                    refinement_level = ((refinement_level & REFINEMENT_Pp1) && depth == 2) ? refinement_level - REFINEMENT_Pp1 :
                        ((refinement_level & REFINEMENT_Pp2) && depth == 1) ? refinement_level - REFINEMENT_Pp2 :
                        ((refinement_level & REFINEMENT_Pp3) && depth == 0) ? refinement_level - REFINEMENT_Pp3 : refinement_level;
                }

                uint8_t lowest_level = 0x00;

                lowest_level = (refinement_level & REFINEMENT_Pp3) ? REFINEMENT_Pp3 : (refinement_level & REFINEMENT_Pp2) ? REFINEMENT_Pp2 : (refinement_level & REFINEMENT_Pp1) ? REFINEMENT_Pp1 :
                    (refinement_level & REFINEMENT_P) ? REFINEMENT_P :
                    (refinement_level & REFINEMENT_Pm1) ? REFINEMENT_Pm1 : (refinement_level & REFINEMENT_Pm2) ? REFINEMENT_Pm2 : (refinement_level & REFINEMENT_Pm3) ? REFINEMENT_Pm3 : 0x00;

                mdc_refinement(
                    &(*context_ptr->local_blk_array),
                    blk_index,
                    depth,
                    refinement_level,
                    lowest_level);
            }

            blk_index += depth_offset[depth];
        }
        else
            blk_index++;
    } // End while 1 CU Loop
}

void forward_blk_to_mode_decision(
    SequenceControlSet                   *scs_ptr,
    PictureControlSet                    *pcs_ptr,
    uint32_t                                  sb_index,
    ModeDecisionConfigurationContext     *context_ptr
)
{
    uint8_t                   blk_index = 0;
    uint32_t                  blk_class = DO_NOT_ADD_CU_CONTINUE_SPLIT;
    EbBool                 split_flag = EB_TRUE;
    MdcSbData           *results_ptr = &pcs_ptr->mdc_sb_array[sb_index];
    SbParams            *sb_params = &scs_ptr->sb_params_array[sb_index];
    MdcpLocalBlkStruct  *local_blk_array = context_ptr->local_blk_array;
    EB_SLICE                slice_type = pcs_ptr->slice_type;

    // CU Loop
    const CodedBlockStats *blk_stats_ptr = get_coded_blk_stats(0);

    results_ptr->leaf_count = 0;
    uint8_t   enable_blk_4x4 = 0;
    blk_index = 0;

    while (blk_index < CU_MAX_COUNT)
    {
        split_flag = EB_TRUE;
        if (sb_params->raster_scan_blk_validity[md_scan_to_raster_scan[blk_index]])
        {
            blk_stats_ptr = get_coded_blk_stats(blk_index);

            switch (blk_stats_ptr->depth) {
            case 0:
            case 1:
            case 2:

                blk_class = DO_NOT_ADD_CU_CONTINUE_SPLIT;

                if (slice_type == I_SLICE) {
                    blk_class = local_blk_array[blk_index].selected_cu == EB_TRUE ? ADD_CU_CONTINUE_SPLIT : blk_class;
                    blk_class = local_blk_array[blk_index].stop_split == EB_TRUE ? ADD_CU_STOP_SPLIT : blk_class;
                }
                else {
                    blk_class = local_blk_array[blk_index].selected_cu == EB_TRUE ? ADD_CU_CONTINUE_SPLIT : blk_class;
                    blk_class = local_blk_array[blk_index].stop_split == EB_TRUE ? ADD_CU_STOP_SPLIT : blk_class;
                }

                // Take into account MAX CU size & MAX intra size (from the API)
                blk_class = (blk_stats_ptr->size > scs_ptr->max_blk_size || (slice_type == I_SLICE && blk_stats_ptr->size > scs_ptr->max_intra_size)) ?
                    DO_NOT_ADD_CU_CONTINUE_SPLIT :
                    blk_class;

                // Take into account MIN CU size & Min intra size(from the API)
                blk_class = (blk_stats_ptr->size == scs_ptr->min_blk_size || (slice_type == I_SLICE && blk_stats_ptr->size == scs_ptr->min_intra_size)) ?
                    ADD_CU_STOP_SPLIT :
                    blk_class;

                switch (blk_class) {
                case ADD_CU_STOP_SPLIT:
                    // Stop
                    results_ptr->leaf_data_array[results_ptr->leaf_count].leaf_index = blk_index;
                    results_ptr->leaf_data_array[results_ptr->leaf_count].mds_idx = pa_to_ep_block_index[blk_index];
                    results_ptr->leaf_data_array[results_ptr->leaf_count].tot_d1_blocks = 1;
                    results_ptr->leaf_data_array[results_ptr->leaf_count++].split_flag = split_flag = EB_FALSE;

                    break;

                case ADD_CU_CONTINUE_SPLIT:
                    // Go Down + consider the current CU as candidate
                    results_ptr->leaf_data_array[results_ptr->leaf_count].leaf_index = blk_index;
                    results_ptr->leaf_data_array[results_ptr->leaf_count].mds_idx = pa_to_ep_block_index[blk_index];
                    results_ptr->leaf_data_array[results_ptr->leaf_count].tot_d1_blocks = 1;
                    results_ptr->leaf_data_array[results_ptr->leaf_count++].split_flag = split_flag = EB_TRUE;

                    break;

                case DO_NOT_ADD_CU_CONTINUE_SPLIT:
                    // Go Down + do not consider the current CU as candidate
                    split_flag = EB_TRUE;

                    break;

                default:
                    results_ptr->leaf_data_array[results_ptr->leaf_count].leaf_index = blk_index;
                    results_ptr->leaf_data_array[results_ptr->leaf_count].mds_idx = pa_to_ep_block_index[blk_index];
                    results_ptr->leaf_data_array[results_ptr->leaf_count].tot_d1_blocks = 1;
                    results_ptr->leaf_data_array[results_ptr->leaf_count++].split_flag = split_flag = EB_TRUE;

                    break;
                }

                break;
            case 3:

                results_ptr->leaf_data_array[results_ptr->leaf_count].leaf_index = blk_index;
                results_ptr->leaf_data_array[results_ptr->leaf_count].mds_idx = pa_to_ep_block_index[blk_index];
                results_ptr->leaf_data_array[results_ptr->leaf_count].tot_d1_blocks = 1;

                if (enable_blk_4x4) {
                    results_ptr->leaf_data_array[results_ptr->leaf_count++].split_flag = split_flag = EB_TRUE;

                    int first_4_index = pa_to_ep_block_index[blk_index] + d1_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_stats_ptr->depth];
                    for (int i = 0; i < 4; ++i) {
                        results_ptr->leaf_data_array[results_ptr->leaf_count].leaf_index = blk_index;

                        results_ptr->leaf_data_array[results_ptr->leaf_count].mds_idx = first_4_index + i;
                        results_ptr->leaf_data_array[results_ptr->leaf_count].tot_d1_blocks = 1;

                        results_ptr->leaf_data_array[results_ptr->leaf_count++].split_flag = split_flag = EB_FALSE;
                    }
                }else
                    results_ptr->leaf_data_array[results_ptr->leaf_count++].split_flag = split_flag = EB_FALSE;

                break;

            default:
                results_ptr->leaf_data_array[results_ptr->leaf_count].leaf_index = blk_index;
                results_ptr->leaf_data_array[results_ptr->leaf_count].mds_idx = pa_to_ep_block_index[blk_index];
                results_ptr->leaf_data_array[results_ptr->leaf_count].tot_d1_blocks = 1;
                results_ptr->leaf_data_array[results_ptr->leaf_count++].split_flag = split_flag = EB_TRUE;
                break;
            }
        }

        blk_index += (split_flag == EB_TRUE) ? 1 : depth_offset[blk_stats_ptr->depth];
    } // End CU Loop
}

void mdc_inter_depth_decision(
    ModeDecisionConfigurationContext     *context_ptr,
    uint32_t                                 origin_x,
    uint32_t                                 origin_y,
    uint32_t                                 end_depth,
    uint32_t                                 blk_index)
{
    uint32_t               left_blk_index;
    uint32_t               top_blk_index;
    uint32_t               top_left_blk_index;
    uint32_t               depth_0_cand_blk_idx;
    uint32_t               depth_1_cand_blk_idx = blk_index;
    uint32_t               depth_2_cand_blk_idx = blk_index;
    uint64_t               depth_n_rate = 0;
    uint64_t               depth_n_plus_1_rate = 0;
    uint64_t               depth_n_cost = 0;
    uint64_t               depth_n_plus_1_cost = 0;
    MdcpLocalBlkStruct *local_blk_array = context_ptr->local_blk_array;
    /*** Stage 0: Inter depth decision: depth 2 vs depth 3 ***/
    // Walks to the last coded 8x8 block for merging
    uint8_t  group_of8x8_blocks_count = context_ptr->group_of8x8_blocks_count;
    uint8_t  group_of16x16_blocks_count = context_ptr->group_of16x16_blocks_count;
    if ((GROUP_OF_4_8x8_BLOCKS(origin_x, origin_y))) {
        group_of8x8_blocks_count++;

        // From the last coded cu index, get the indices of the left, top, and top left cus
        left_blk_index = blk_index - DEPTH_THREE_STEP;
        top_blk_index = left_blk_index - DEPTH_THREE_STEP;
        top_left_blk_index = top_blk_index - DEPTH_THREE_STEP;

        // From the top left index, get the index of the candidate pu for merging
        depth_2_cand_blk_idx = top_left_blk_index - 1;

        // Compute depth N cost
        local_blk_array[depth_2_cand_blk_idx].split_context = 0;
        depth_n_cost = (local_blk_array[depth_2_cand_blk_idx]).early_cost + depth_n_rate;

        if (end_depth < 3) {
            (local_blk_array[depth_2_cand_blk_idx]).early_split_flag = EB_FALSE;
            (local_blk_array[depth_2_cand_blk_idx]).early_cost = depth_n_cost;
        }
        else {
            depth_n_plus_1_cost = (local_blk_array[blk_index]).early_cost + (local_blk_array[left_blk_index]).early_cost + (local_blk_array[top_blk_index]).early_cost + (local_blk_array[top_left_blk_index]).early_cost + depth_n_plus_1_rate;

            if (depth_n_cost <= depth_n_plus_1_cost) {
                // If the cost is low enough to warrant not spliting further:
                // 1. set the split flag of the candidate pu for merging to false
                // 2. update the last pu index
                (local_blk_array[depth_2_cand_blk_idx]).early_split_flag = EB_FALSE;
                (local_blk_array[depth_2_cand_blk_idx]).early_cost = depth_n_cost;
            }
            else {
                // If the cost is not low enough:
                // update the cost of the candidate pu for merging
                // this update is required for the next inter depth decision
                (&local_blk_array[depth_2_cand_blk_idx])->early_cost = depth_n_plus_1_cost;
            }
        }
    }

    // Walks to the last coded 16x16 block for merging
    if (GROUP_OF_4_16x16_BLOCKS(get_coded_blk_stats(depth_2_cand_blk_idx)->origin_x, get_coded_blk_stats(depth_2_cand_blk_idx)->origin_y) &&
        (group_of8x8_blocks_count == 4)) {
        group_of8x8_blocks_count = 0;
        group_of16x16_blocks_count++;

        // From the last coded pu index, get the indices of the left, top, and top left pus
        left_blk_index = depth_2_cand_blk_idx - DEPTH_TWO_STEP;
        top_blk_index = left_blk_index - DEPTH_TWO_STEP;
        top_left_blk_index = top_blk_index - DEPTH_TWO_STEP;

        // From the top left index, get the index of the candidate pu for merging
        depth_1_cand_blk_idx = top_left_blk_index - 1;

        if (get_coded_blk_stats(depth_1_cand_blk_idx)->depth == 1) {
            depth_n_cost = local_blk_array[depth_1_cand_blk_idx].early_cost + depth_n_rate;
            if (end_depth < 2) {
                local_blk_array[depth_1_cand_blk_idx].early_split_flag = EB_FALSE;
                local_blk_array[depth_1_cand_blk_idx].early_cost = depth_n_cost;
            }
            else {
                // Compute depth N+1 cost
                depth_n_plus_1_cost = local_blk_array[depth_2_cand_blk_idx].early_cost +
                    local_blk_array[left_blk_index].early_cost +
                    local_blk_array[top_blk_index].early_cost +
                    local_blk_array[top_left_blk_index].early_cost +
                    depth_n_plus_1_rate;

                // Inter depth comparison: depth 1 vs depth 2
                if (depth_n_cost <= depth_n_plus_1_cost) {
                    // If the cost is low enough to warrant not spliting further:
                    // 1. set the split flag of the candidate pu for merging to false
                    // 2. update the last pu index
                    local_blk_array[depth_1_cand_blk_idx].early_split_flag = EB_FALSE;
                    local_blk_array[depth_1_cand_blk_idx].early_cost = depth_n_cost;
                }
                else {
                    // If the cost is not low enough:
                    // update the cost of the candidate pu for merging
                    // this update is required for the next inter depth decision
                    local_blk_array[depth_1_cand_blk_idx].early_cost = depth_n_plus_1_cost;
                }
            }
        }
    }

    // Stage 2: Inter depth decision: depth 0 vs depth 1

    // Walks to the last coded 32x32 block for merging
    // Stage 2 isn't performed in I slices since the abcense of 64x64 candidates
    if (GROUP_OF_4_32x32_BLOCKS(get_coded_blk_stats(depth_1_cand_blk_idx)->origin_x, get_coded_blk_stats(depth_1_cand_blk_idx)->origin_y) &&
        (group_of16x16_blocks_count == 4)) {
        group_of16x16_blocks_count = 0;

        // From the last coded pu index, get the indices of the left, top, and top left pus
        left_blk_index = depth_1_cand_blk_idx - DEPTH_ONE_STEP;
        top_blk_index = left_blk_index - DEPTH_ONE_STEP;
        top_left_blk_index = top_blk_index - DEPTH_ONE_STEP;

        // From the top left index, get the index of the candidate pu for merging
        depth_0_cand_blk_idx = top_left_blk_index - 1;

        if (get_coded_blk_stats(depth_0_cand_blk_idx)->depth == 0) {
            // Compute depth N cost
            depth_n_cost = (&local_blk_array[depth_0_cand_blk_idx])->early_cost + depth_n_rate;
            if (end_depth < 1)
                (&local_blk_array[depth_0_cand_blk_idx])->early_split_flag = EB_FALSE;
            else {
                // Compute depth N+1 cost
                depth_n_plus_1_cost = local_blk_array[depth_1_cand_blk_idx].early_cost +
                    local_blk_array[left_blk_index].early_cost +
                    local_blk_array[top_blk_index].early_cost +
                    local_blk_array[top_left_blk_index].early_cost +
                    depth_n_plus_1_rate;

                // Inter depth comparison: depth 0 vs depth 1
                if (depth_n_cost <= depth_n_plus_1_cost) {
                    // If the cost is low enough to warrant not spliting further:
                    // 1. set the split flag of the candidate pu for merging to false
                    // 2. update the last pu index
                    (&local_blk_array[depth_0_cand_blk_idx])->early_split_flag = EB_FALSE;
                }
            }
        }
    }

    context_ptr->group_of8x8_blocks_count = group_of8x8_blocks_count;
    context_ptr->group_of16x16_blocks_count = group_of16x16_blocks_count;
}

void prediction_partition_loop(
    SequenceControlSet                   *scs_ptr,
    PictureControlSet                    *pcs_ptr,
    uint32_t                                sb_index,
    uint32_t                                tb_origin_x,
    uint32_t                                tb_origin_y,
    uint32_t                                start_depth,
    uint32_t                                end_depth,
    ModeDecisionConfigurationContext     *context_ptr){
    MdcpLocalBlkStruct *local_blk_array = context_ptr->local_blk_array;
    MdcpLocalBlkStruct   *blk_ptr;

    SbParams *sb_params = &scs_ptr->sb_params_array[sb_index];
    uint32_t      raster_scan_blk_index;
    uint32_t      blk_index = 0;
    uint32_t      start_index = 0;

    (void)tb_origin_x;
    (void)tb_origin_y;

    const CodedBlockStats *blk_stats_ptr;

    FrameHeader *frm_hdr = &pcs_ptr->parent_pcs_ptr->frm_hdr;

    for (blk_index = start_index; blk_index < CU_MAX_COUNT; ++blk_index)

    {
        local_blk_array[blk_index].selected_cu = EB_FALSE;
        local_blk_array[blk_index].stop_split = EB_FALSE;

        blk_ptr = &local_blk_array[blk_index];
        raster_scan_blk_index = md_scan_to_raster_scan[blk_index];
        if (sb_params->raster_scan_blk_validity[raster_scan_blk_index])
        {
            uint32_t depth;
            blk_stats_ptr = get_coded_blk_stats(blk_index);

            depth = blk_stats_ptr->depth;
            blk_ptr->early_split_flag = (depth < end_depth) ? EB_TRUE : EB_FALSE;

            if (depth >= start_depth && depth <= end_depth) {
                //reset the flags here:   all CU splitFalg=TRUE. default: we always split. interDepthDecision will select where  to stop splitting(ie setting the flag to False)

                if (pcs_ptr->slice_type != I_SLICE) {
                    const MeSbResults *me_results = pcs_ptr->parent_pcs_ptr->me_results[sb_index];
                    const MeCandidate *me_block_results = me_results->me_candidate[raster_scan_blk_index];
                    uint8_t total_me_cnt = me_results->total_me_candidate_index[raster_scan_blk_index];
                    uint8_t me_index = 0;
                    for (uint8_t me_candidate_index = 0; me_candidate_index < total_me_cnt; me_candidate_index++) {
                        const MeCandidate *me_block_results_ptr = &me_block_results[me_candidate_index];
                        if (frm_hdr->reference_mode == SINGLE_REFERENCE) {
                            if (me_block_results_ptr->direction == 0) {
                                me_index = me_candidate_index;
                                break;
                            }
                        }
                        else {
                            if (me_block_results_ptr->direction == 2) {
                                me_index = me_candidate_index;
                                break;
                            }
                        }
                    }

                    //const MeCandidate_t *me_results = &me_block_results[me_index];

                    // Initialize the mdc candidate (only av1 rate estimation inputs)
                    // Hsan: mode, direction, .. could be modified toward better early inter depth decision (e.g. NEARESTMV instead of NEWMV)
                    context_ptr->mdc_candidate_ptr->md_rate_estimation_ptr = context_ptr->md_rate_estimation_ptr;
                    context_ptr->mdc_candidate_ptr->type = INTER_MODE;
                    context_ptr->mdc_candidate_ptr->merge_flag = EB_FALSE;
                    context_ptr->mdc_candidate_ptr->prediction_direction[0] = (pcs_ptr->parent_pcs_ptr->temporal_layer_index == 0) ?
                        UNI_PRED_LIST_0 :
                        me_block_results[me_index].direction;
                    // Hsan: what's the best mode for rate simulation
                    context_ptr->mdc_candidate_ptr->inter_mode = NEARESTMV;
                    context_ptr->mdc_candidate_ptr->pred_mode = NEARESTMV;
                    context_ptr->mdc_candidate_ptr->motion_mode = SIMPLE_TRANSLATION;
                    context_ptr->mdc_candidate_ptr->is_new_mv = 1;
                    context_ptr->mdc_candidate_ptr->is_zero_mv = 0;
                    context_ptr->mdc_candidate_ptr->drl_index = 0;
                    context_ptr->mdc_candidate_ptr->motion_vector_xl0 = me_results->me_mv_array[raster_scan_blk_index][0].x_mv << 1;
                    context_ptr->mdc_candidate_ptr->motion_vector_yl0 = me_results->me_mv_array[raster_scan_blk_index][0].y_mv << 1;
                    context_ptr->mdc_candidate_ptr->motion_vector_xl1 = me_results->me_mv_array[raster_scan_blk_index][((scs_ptr->mrp_mode == 0) ? 4 : 2)].x_mv << 1;
                    context_ptr->mdc_candidate_ptr->motion_vector_yl1 = me_results->me_mv_array[raster_scan_blk_index][((scs_ptr->mrp_mode == 0) ? 4 : 2)].y_mv << 1;
                    context_ptr->mdc_candidate_ptr->ref_mv_index = 0;
                    context_ptr->mdc_candidate_ptr->pred_mv_weight = 0;
                    if (context_ptr->mdc_candidate_ptr->prediction_direction[0] == BI_PRED) {
                        context_ptr->mdc_candidate_ptr->ref_frame_type = LAST_BWD_FRAME;
                        context_ptr->mdc_candidate_ptr->is_compound = 1;
                    }
                    else if (context_ptr->mdc_candidate_ptr->prediction_direction[0] == UNI_PRED_LIST_0) {
                        context_ptr->mdc_candidate_ptr->ref_frame_type = LAST_FRAME;
                        context_ptr->mdc_candidate_ptr->is_compound = 0;
                    }
                    else { // context_ptr->mdc_candidate_ptr->prediction_direction[0]
                        context_ptr->mdc_candidate_ptr->ref_frame_type = BWDREF_FRAME;
                        context_ptr->mdc_candidate_ptr->is_compound = 0;
                    }
                    context_ptr->mdc_candidate_ptr->motion_vector_pred_x[REF_LIST_0] = 0;
                    context_ptr->mdc_candidate_ptr->motion_vector_pred_y[REF_LIST_0] = 0;
                    // Initialize the ref mv
                    memset(context_ptr->mdc_ref_mv_stack,0,sizeof(CandidateMv));
                    context_ptr->blk_geom = get_blk_geom_mds(pa_to_ep_block_index[blk_index]);
                    // Initialize mdc cu (only av1 rate estimation inputs)
                    context_ptr->mdc_blk_ptr->is_inter_ctx = 0;
                    context_ptr->mdc_blk_ptr->skip_flag_context = 0;
                    context_ptr->mdc_blk_ptr->inter_mode_ctx[context_ptr->mdc_candidate_ptr->ref_frame_type] = 0;
                    context_ptr->mdc_blk_ptr->reference_mode_context = 0;
                    context_ptr->mdc_blk_ptr->compoud_reference_type_context = 0;
                    av1_zero(context_ptr->mdc_blk_ptr->av1xd->neighbors_ref_counts); // Hsan: neighbor not generated @ open loop partitioning => assumes always (0,0)

                    // Fast Cost Calc
                    blk_ptr->early_cost = av1_inter_fast_cost(
                        context_ptr->mdc_blk_ptr,
                        context_ptr->mdc_candidate_ptr,
                        context_ptr->qp,
                        me_block_results[me_index].distortion,
                        (uint64_t) 0,
                        context_ptr->lambda,
                        0,
                        pcs_ptr,
                        context_ptr->mdc_ref_mv_stack,
                        context_ptr->blk_geom,
                        (tb_origin_y + context_ptr->blk_geom->origin_y) >> MI_SIZE_LOG2,
                        (tb_origin_x + context_ptr->blk_geom->origin_x) >> MI_SIZE_LOG2,
                        0,
                        0,
                        0,
                        DC_PRED,        // Hsan: neighbor not generated @ open loop partitioning
                        DC_PRED);       // Hsan: neighbor not generated @ open loop partitioning
                }

                if (end_depth == 2)
                    context_ptr->group_of8x8_blocks_count = depth == 2 ? incremental_count[raster_scan_blk_index] : 0;
                if (end_depth == 1)
                    context_ptr->group_of16x16_blocks_count = depth == 1 ? incremental_count[raster_scan_blk_index] : 0;
                mdc_inter_depth_decision(
                    context_ptr,
                    blk_stats_ptr->origin_x,
                    blk_stats_ptr->origin_y,
                    end_depth,
                    blk_index);
            }
            else
                blk_ptr->early_cost = ~0u;
        }
    }// End CU Loop
}

EbErrorType early_mode_decision_sb(
    SequenceControlSet                   *scs_ptr,
    PictureControlSet                    *pcs_ptr,
    SuperBlock                           *sb_ptr,
    uint32_t                                  sb_index,
    ModeDecisionConfigurationContext     *context_ptr){
    EbErrorType    return_error = EB_ErrorNone;
    uint32_t       tb_origin_x    = sb_ptr->origin_x;
    uint32_t       tb_origin_y    = sb_ptr->origin_y;

    uint32_t      start_depth = DEPTH_64;

    uint32_t      end_depth =  DEPTH_8 ;
    context_ptr->group_of8x8_blocks_count = 0;
    context_ptr->group_of16x16_blocks_count = 0;

    prediction_partition_loop(
        scs_ptr,
        pcs_ptr,
        sb_index,
        tb_origin_x,
        tb_origin_y,
        start_depth,
        end_depth,
        context_ptr
    );

    refinement_prediction_loop(
        scs_ptr,
        pcs_ptr,
        sb_index,
        context_ptr);

    forward_blk_to_mode_decision(
        scs_ptr,
        pcs_ptr,
        sb_index,
        context_ptr);

    return return_error;
}

// clang-format on
