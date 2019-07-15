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
#include "EbResourceCoordinationProcess.h"
#include "EbPictureAnalysisProcess.h"
#include "EbResourceCoordinationResults.h"
#include "EbPictureDemuxResults.h"
#include "EbRateControlResults.h"
#include "EbPictureDecisionProcess.h"
#include "EbMotionEstimationProcess.h"
#include "EbInitialRateControlProcess.h"
#include "EbSourceBasedOperationsProcess.h"
#include "EbPictureManagerProcess.h"
#include "EbRateControlProcess.h"
#include "EbModeDecisionConfigurationProcess.h"
#include "EbEncDecProcess.h"
#include "EbDlfProcess.h"
#include "EbCdefProcess.h"
#include "EbRestProcess.h"
#include "EbEntropyCodingProcess.h"
#include "EbPacketizationProcess.h"
#include "EbObject.h"

/**************************************
 * Component Private Data
 **************************************/
typedef struct EbEncHandle
{
    EbDctor                                   dctor;
    // Encode Instances & Compute Segments
    uint32_t                                  encode_instance_total_count;
    uint32_t                                  compute_segments_total_count_array;

    // Config Set Counts
    uint32_t                                  sequence_control_set_pool_total_count;

    // Full Results Count
    uint32_t                                  picture_control_set_pool_total_count;

    // Picture Buffer Count
    uint32_t                                  reference_picture_pool_total_count;

    // Config Set Pool & Active Array
    EbSystemResource                         *sequence_control_set_pool_ptr;
    EbFifo                                  **sequence_control_set_pool_producer_fifo_ptr_array;
    EbSequenceControlSetInstance            **sequence_control_set_instance_array;

    // Full Results
    EbSystemResource                    **picture_control_set_pool_ptr_array;
    EbFifo                             ***picture_control_set_pool_producer_fifo_ptr_dbl_array;

    //ParentControlSet
    EbSystemResource                    **picture_parent_control_set_pool_ptr_array;
    EbFifo                             ***picture_parent_control_set_pool_producer_fifo_ptr_dbl_array;

    // Picture Buffers
    EbSystemResource                    **reference_picture_pool_ptr_array;
    EbSystemResource                    **pa_reference_picture_pool_ptr_array;

    // Overlay input picture
    EbSystemResource                    **overlay_input_picture_pool_ptr_array;
    EbFifo                             ***overlay_input_picture_pool_producer_fifo_ptr_dbl_array;
    // Picture Buffer Producer Fifos
    EbFifo                             ***reference_picture_pool_producer_fifo_ptr_dbl_array;
    EbFifo                             ***pa_reference_picture_pool_producer_fifo_ptr_dbl_array;

    // Thread Handles
    EbHandle                               resource_coordination_thread_handle;
    EbHandle                              *picture_analysis_thread_handle_array;
    EbHandle                               picture_decision_thread_handle;
    EbHandle                              *motion_estimation_thread_handle_array;
    EbHandle                               initial_rate_control_thread_handle;
    EbHandle                              *source_based_operations_thread_handle_array;
    EbHandle                               picture_manager_thread_handle;
    EbHandle                               rate_control_thread_handle;
    EbHandle                              *mode_decision_configuration_thread_handle_array;
    EbHandle                              *enc_dec_thread_handle_array;
    EbHandle                              *entropy_coding_thread_handle_array;
    EbHandle                              *dlf_thread_handle_array;
    EbHandle                              *cdef_thread_handle_array;
    EbHandle                              *rest_thread_handle_array;

    EbHandle                               packetization_thread_handle;

    // Contexts
    ResourceCoordinationContext            *resource_coordination_context_ptr;
    PictureAnalysisContext                 **picture_analysis_context_ptr_array;
    PictureDecisionContext                *picture_decision_context_ptr;
    MotionEstimationContext_t             **motion_estimation_context_ptr_array;
    InitialRateControlContext             *initial_rate_control_context_ptr;
    SourceBasedOperationsContext         **source_based_operations_context_ptr_array;
    PictureManagerContext                 *picture_manager_context_ptr;
    RateControlContext                    *rate_control_context_ptr;
    ModeDecisionConfigurationContext     **mode_decision_configuration_context_ptr_array;
    EncDecContext                        **enc_dec_context_ptr_array;
    EntropyCodingContext                 **entropy_coding_context_ptr_array;
    DlfContext                           **dlf_context_ptr_array;
    CdefContext_t                        **cdef_context_ptr_array;
    RestContext                          **rest_context_ptr_array;
    PacketizationContext                  *packetization_context_ptr;

    // System Resource Managers
    EbSystemResource                     *input_buffer_resource_ptr;
    EbSystemResource                    **output_stream_buffer_resource_ptr_array;
    EbSystemResource                    **output_recon_buffer_resource_ptr_array;
    EbSystemResource                    **output_statistics_buffer_resource_ptr_array;
    EbSystemResource                     *resource_coordination_results_resource_ptr;
    EbSystemResource                     *picture_analysis_results_resource_ptr;
    EbSystemResource                     *picture_decision_results_resource_ptr;
    EbSystemResource                     *motion_estimation_results_resource_ptr;
    EbSystemResource                     *initial_rate_control_results_resource_ptr;
    EbSystemResource                     *picture_demux_results_resource_ptr;
    EbSystemResource                     *rate_control_tasks_resource_ptr;
    EbSystemResource                     *rate_control_results_resource_ptr;
    EbSystemResource                     *enc_dec_tasks_resource_ptr;
    EbSystemResource                     *enc_dec_results_resource_ptr;
    EbSystemResource                     *entropy_coding_results_resource_ptr;
    EbSystemResource                     *dlf_results_resource_ptr;
    EbSystemResource                     *cdef_results_resource_ptr;
    EbSystemResource                     *rest_results_resource_ptr;

    // Inter-Process Producer Fifos
    EbFifo                              **input_buffer_producer_fifo_ptr_array;
    EbFifo                             ***output_stream_buffer_producer_fifo_ptr_dbl_array;
    EbFifo                             ***output_recon_buffer_producer_fifo_ptr_dbl_array;
    EbFifo                             ***output_statistics_buffer_producer_fifo_ptr_dbl_array;
    EbFifo                              **resource_coordination_results_producer_fifo_ptr_array;
    EbFifo                              **picture_analysis_results_producer_fifo_ptr_array;
    EbFifo                              **picture_decision_results_producer_fifo_ptr_array;
    EbFifo                              **motion_estimation_results_producer_fifo_ptr_array;
    EbFifo                              **initial_rate_control_results_producer_fifo_ptr_array;
    EbFifo                              **picture_demux_results_producer_fifo_ptr_array;
    EbFifo                              **picture_manager_results_producer_fifo_ptr_array;
    EbFifo                              **rate_control_tasks_producer_fifo_ptr_array;
    EbFifo                              **rate_control_results_producer_fifo_ptr_array;
    EbFifo                              **enc_dec_tasks_producer_fifo_ptr_array;
    EbFifo                              **enc_dec_results_producer_fifo_ptr_array;
    EbFifo                              **entropy_coding_results_producer_fifo_ptr_array;
    EbFifo                              **dlf_results_producer_fifo_ptr_array;
    EbFifo                              **cdef_results_producer_fifo_ptr_array;
    EbFifo                              **rest_results_producer_fifo_ptr_array;

    // Inter-Process Consumer Fifos
    EbFifo                              **input_buffer_consumer_fifo_ptr_array;
    EbFifo                             ***output_stream_buffer_consumer_fifo_ptr_dbl_array;
    EbFifo                             ***output_recon_buffer_consumer_fifo_ptr_dbl_array;
    EbFifo                             ***output_statistics_buffer_consumer_fifo_ptr_dbl_array;
    EbFifo                              **resource_coordination_results_consumer_fifo_ptr_array;
    EbFifo                              **picture_analysis_results_consumer_fifo_ptr_array;
    EbFifo                              **picture_decision_results_consumer_fifo_ptr_array;
    EbFifo                              **motion_estimation_results_consumer_fifo_ptr_array;
    EbFifo                              **initial_rate_control_results_consumer_fifo_ptr_array;
    EbFifo                              **picture_demux_results_consumer_fifo_ptr_array;
    EbFifo                              **rate_control_tasks_consumer_fifo_ptr_array;
    EbFifo                              **rate_control_results_consumer_fifo_ptr_array;
    EbFifo                              **enc_dec_tasks_consumer_fifo_ptr_array;
    EbFifo                              **enc_dec_results_consumer_fifo_ptr_array;
    EbFifo                              **entropy_coding_results_consumer_fifo_ptr_array;
    EbFifo                              **dlf_results_consumer_fifo_ptr_array;
    EbFifo                              **cdef_results_consumer_fifo_ptr_array;
    EbFifo                              **rest_results_consumer_fifo_ptr_array;

    // Callbacks
    EbCallback                          **app_callback_ptr_array;

} EbEncHandle;

#endif // EbEncHandle_h
