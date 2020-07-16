/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbEncHandle_h
#define EbEncHandle_h

#include "EbDefinitions.h"
#include "EbSvtAv1Enc.h"
#include "EbPictureBufferDesc.h"
#include "EbSystemResourceManager.h"
#include "EbSequenceControlSet.h"
#include "EbObject.h"

struct _EbThreadContext {
    EbDctor dctor;
    EbPtr   priv;
};

/**************************************
 * Component Private Data
 **************************************/
struct _EbEncHandle {
    EbDctor dctor;
    // Encode Instances & Compute Segments
    uint32_t encode_instance_total_count;
    uint32_t compute_segments_total_count_array;

    // Config Set Counts
    uint32_t scs_pool_total_count;

    // Full Results Count
    uint32_t pcs_pool_total_count;

    // Picture Buffer Count
    uint32_t ref_pic_pool_total_count;

    // Config Set Pool & Active Array
    EbSystemResource *             scs_pool_ptr; // sequence_control_set_pool
    EbSequenceControlSetInstance **scs_instance_array;

    // Full Results
    EbSystemResource **picture_control_set_pool_ptr_array;

    //ParentControlSet
    EbSystemResource **picture_parent_control_set_pool_ptr_array;
#if DECOUPLE_ME_RES
    EbSystemResource **me_pool_ptr_array;
#endif
    // Picture Buffers
    EbSystemResource **reference_picture_pool_ptr_array;
    EbSystemResource **pa_reference_picture_pool_ptr_array;

    // Overlay input picture
    EbSystemResource **overlay_input_picture_pool_ptr_array;

    // Thread Handles
    EbHandle  resource_coordination_thread_handle;
    EbHandle *picture_analysis_thread_handle_array;
    EbHandle  picture_decision_thread_handle;
    EbHandle *motion_estimation_thread_handle_array;
    EbHandle  initial_rate_control_thread_handle;
    EbHandle *source_based_operations_thread_handle_array;
    EbHandle  picture_manager_thread_handle;
    EbHandle  rate_control_thread_handle;
    EbHandle *mode_decision_configuration_thread_handle_array;
    EbHandle *enc_dec_thread_handle_array;
    EbHandle *entropy_coding_thread_handle_array;
    EbHandle *dlf_thread_handle_array;
    EbHandle *cdef_thread_handle_array;
    EbHandle *rest_thread_handle_array;

    EbHandle packetization_thread_handle;

    // Contexts
    EbThreadContext * resource_coordination_context_ptr;
    EbThreadContext **picture_analysis_context_ptr_array;
    EbThreadContext * picture_decision_context_ptr;
    EbThreadContext **motion_estimation_context_ptr_array;
    EbThreadContext * initial_rate_control_context_ptr;
    EbThreadContext **source_based_operations_context_ptr_array;
    EbThreadContext * picture_manager_context_ptr;
    EbThreadContext * rate_control_context_ptr;
    EbThreadContext **mode_decision_configuration_context_ptr_array;
    EbThreadContext **enc_dec_context_ptr_array;
    EbThreadContext **entropy_coding_context_ptr_array;
    EbThreadContext **dlf_context_ptr_array;
    EbThreadContext **cdef_context_ptr_array;
    EbThreadContext **rest_context_ptr_array;
    EbThreadContext * packetization_context_ptr;

    // System Resource Managers
    EbSystemResource * input_buffer_resource_ptr;
    EbSystemResource **output_stream_buffer_resource_ptr_array;
    EbSystemResource **output_recon_buffer_resource_ptr_array;
    EbSystemResource **output_statistics_buffer_resource_ptr_array;
    EbSystemResource * resource_coordination_results_resource_ptr;
    EbSystemResource * picture_analysis_results_resource_ptr;
    EbSystemResource * picture_decision_results_resource_ptr;
    EbSystemResource * motion_estimation_results_resource_ptr;
    EbSystemResource * initial_rate_control_results_resource_ptr;
    EbSystemResource * picture_demux_results_resource_ptr;
    EbSystemResource * rate_control_tasks_resource_ptr;
    EbSystemResource * rate_control_results_resource_ptr;
    EbSystemResource * enc_dec_tasks_resource_ptr;
    EbSystemResource * enc_dec_results_resource_ptr;
    EbSystemResource * entropy_coding_results_resource_ptr;
    EbSystemResource * dlf_results_resource_ptr;
    EbSystemResource * cdef_results_resource_ptr;
    EbSystemResource * rest_results_resource_ptr;

    // Callbacks
    EbCallback **app_callback_ptr_array;

    EbFifo *input_buffer_producer_fifo_ptr;
    EbFifo *output_stream_buffer_consumer_fifo_ptr;
    EbFifo *output_recon_buffer_consumer_fifo_ptr;
};

#endif // EbEncHandle_h
