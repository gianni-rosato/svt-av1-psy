/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbResourceCoordination_h
#define EbResourceCoordination_h

#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"
#include "EbDefinitions.h"
#ifdef __cplusplus
extern "C" {
#endif
    /***************************************
     * Context
     ***************************************/
    typedef struct ResourceCoordinationContext_s
    {

        EbFifo_t                            *input_buffer_fifo_ptr;
        EbFifo_t                            *resource_coordination_results_output_fifo_ptr;
        EbFifo_t                           **picture_control_set_fifo_ptr_array;
        EbSequenceControlSetInstance_t     **sequence_control_set_instance_array;
        EbObjectWrapper_t                  **sequenceControlSetActiveArray;
        EbFifo_t                            *sequence_control_set_empty_fifo_ptr;
        EbCallback_t                       **app_callback_ptr_array;
        
        // Compute Segments
        uint32_t                              *compute_segments_total_count_array;
        uint32_t                               encode_instances_total_count;

        // Picture Number Array
        uint64_t                              *pictureNumberArray;

        uint64_t                               averageEncMod;
        uint8_t                                prevEncMod;
        int8_t                                prevEncModeDelta;
        uint8_t                                prevChangeCond;

        int64_t                               previousModeChangeBuffer;
        int64_t                               previousModeChangeFrameIn;
        int64_t                               previousBufferCheck1;
        int64_t                               previousFrameInCheck1;
        int64_t                               previousFrameInCheck2;
        int64_t                               previousFrameInCheck3;


        uint64_t                               curSpeed; // speed x 1000
        uint64_t                               prevsTimeSeconds;
        uint64_t                               prevsTimeuSeconds;
        int64_t                               prevFrameOut;

        uint64_t                               firstInPicArrivedTimeSeconds;
        uint64_t                               firstInPicArrivedTimeuSeconds;
        EbBool                              startFlag;



    } ResourceCoordinationContext_t;

    /***************************************
     * Extern Function Declaration
     ***************************************/
    extern EbErrorType resource_coordination_context_ctor(
        ResourceCoordinationContext_t  **context_dbl_ptr,
        EbFifo_t                        *input_buffer_fifo_ptr,
        EbFifo_t                        *resource_coordination_results_output_fifo_ptr,
        EbFifo_t                       **picture_control_set_fifo_ptr_array,
        EbSequenceControlSetInstance_t **sequence_control_set_instance_array,
        EbFifo_t                        *sequence_control_set_empty_fifo_ptr,
        EbCallback_t                   **app_callback_ptr_array,
        uint32_t                        *compute_segments_total_count_array,
        uint32_t                         encode_instances_total_count);

    extern void* resource_coordination_kernel(void *input_ptr);
#ifdef __cplusplus
}
#endif
#endif // EbResourceCoordination_h