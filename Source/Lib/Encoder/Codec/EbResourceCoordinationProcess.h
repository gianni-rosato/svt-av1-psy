/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbResourceCoordination_h
#define EbResourceCoordination_h

#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"
#include "EbObject.h"
#ifdef __cplusplus
extern "C" {
#endif
    /***************************************
     * Context
     ***************************************/
    typedef struct ResourceCoordinationContext
    {
        EbDctor                               dctor;
        EbFifo                                *input_buffer_fifo_ptr;
        EbFifo                                *resource_coordination_results_output_fifo_ptr;
        EbFifo                               **picture_control_set_fifo_ptr_array;
        EbSequenceControlSetInstance         **sequence_control_set_instance_array;
        EbObjectWrapper                      **sequenceControlSetActiveArray;
        EbFifo                                *sequence_control_set_empty_fifo_ptr;
        EbCallback                           **app_callback_ptr_array;

        // Compute Segments
        uint32_t                               compute_segments_total_count_array;
        uint32_t                               encode_instances_total_count;

        // Picture Number Array
        uint64_t                              *picture_number_array;

        uint64_t                               average_enc_mod;
        uint8_t                                prev_enc_mod;
        int8_t                                 prev_enc_mode_delta;
        uint8_t                                prev_change_cond;

        int64_t                                previous_mode_change_buffer;
        int64_t                                previous_mode_change_frame_in;
        int64_t                                previous_buffer_check1;
        int64_t                                previous_frame_in_check1;
        int64_t                                previous_frame_in_check2;
        int64_t                                previous_frame_in_check3;

        uint64_t                               cur_speed; // speed x 1000
        uint64_t                               prevs_time_seconds;
        uint64_t                               prevs_timeu_seconds;
        int64_t                                prev_frame_out;

        uint64_t                               first_in_pic_arrived_time_seconds;
        uint64_t                               first_in_pic_arrived_timeu_seconds;
        EbBool                                 start_flag;
    } ResourceCoordinationContext;

    /***************************************
     * Extern Function Declaration
     ***************************************/
    extern EbErrorType resource_coordination_context_ctor(
        ResourceCoordinationContext  *context_ptr,
        EbFifo                        *input_buffer_fifo_ptr,
        EbFifo                        *resource_coordination_results_output_fifo_ptr,
        EbFifo                       **picture_control_set_fifo_ptr_array,
        EbSequenceControlSetInstance **sequence_control_set_instance_array,
        EbFifo                        *sequence_control_set_empty_fifo_ptr,
        EbCallback                   **app_callback_ptr_array,
        uint32_t                       compute_segments_total_count_array,
        uint32_t                       encode_instances_total_count);

    extern void* resource_coordination_kernel(void *input_ptr);
#ifdef __cplusplus
}
#endif
#endif // EbResourceCoordination_h
