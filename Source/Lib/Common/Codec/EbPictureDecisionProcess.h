/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbPictureDecision_h
#define EbPictureDecision_h

#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"
#if ALTREF_FILTERING_SUPPORT
#include "EbPictureControlSet.h"
#endif

/**************************************
 * Context
 **************************************/
typedef struct PictureDecisionContext
{
    EbFifo       *picture_analysis_results_input_fifo_ptr;
    EbFifo       *picture_decision_results_output_fifo_ptr;
                 
    uint64_t      last_solid_color_frame_poc;
                 
    EbBool        reset_running_avg;
                 
    uint32_t    **ahd_running_avg_cb;
    uint32_t    **ahd_running_avg_cr;
    uint32_t    **ahd_running_avg;
    EbBool        is_scene_change_detected;

    // Dynamic GOP
    uint32_t        totalRegionActivityCost[MAX_NUMBER_OF_REGIONS_IN_WIDTH][MAX_NUMBER_OF_REGIONS_IN_HEIGHT];

    uint32_t      total_number_of_mini_gops;
                  
    uint32_t      mini_gop_start_index[MINI_GOP_WINDOW_MAX_COUNT];
    uint32_t      mini_gop_end_index[MINI_GOP_WINDOW_MAX_COUNT];
    uint32_t      mini_gop_length[MINI_GOP_WINDOW_MAX_COUNT];
    uint32_t      mini_gop_intra_count[MINI_GOP_WINDOW_MAX_COUNT];
    uint32_t      mini_gop_idr_count[MINI_GOP_WINDOW_MAX_COUNT];
    uint32_t      mini_gop_hierarchical_levels[MINI_GOP_WINDOW_MAX_COUNT];
    EbBool        mini_gop_activity_array[MINI_GOP_MAX_COUNT];
    uint32_t        mini_gop_region_activity_cost_array[MINI_GOP_MAX_COUNT][MAX_NUMBER_OF_REGIONS_IN_WIDTH][MAX_NUMBER_OF_REGIONS_IN_HEIGHT];

    uint32_t        mini_gop_group_faded_in_pictures_count[MINI_GOP_MAX_COUNT];
    uint32_t        mini_gop_group_faded_out_pictures_count[MINI_GOP_MAX_COUNT];
#if NEW_RPS
    uint8_t     lay0_toggle; //3 way toggle 0->1->2
    uint8_t     lay1_toggle; //2 way toggle 0->1
    uint8_t     lay2_toggle; //2 way toggle 0->1
#endif
    EbBool        mini_gop_toggle;    //mini GOP toggling since last Key Frame  K-0-1-0-1-0-K-0-1-0-1-K-0-1.....
#if BASE_LAYER_REF
    uint64_t      last_islice_picture_number;
#endif
    uint8_t       last_i_picture_sc_detection;
#if REF_ORDER
    uint64_t         key_poc;
#endif
} PictureDecisionContext;

/***************************************
 * Extern Function Declaration
 ***************************************/
extern EbErrorType picture_decision_context_ctor(
    PictureDecisionContext **context_dbl_ptr,
    EbFifo                  *picture_analysis_results_input_fifo_ptr,
    EbFifo                  *picture_decision_results_output_fifo_ptr);


extern void* picture_decision_kernel(void *input_ptr);

#if ALTREF_FILTERING_SUPPORT
void DecimateInputPicture(PictureParentControlSet       *picture_control_set_ptr,
                          EbPictureBufferDesc           *inputPaddedPicturePtr,
                          EbPictureBufferDesc           *quarterDecimatedPicturePtr,
                          EbPictureBufferDesc           *sixteenthDecimatedPicturePtr);
#endif

#endif // EbPictureDecision_h