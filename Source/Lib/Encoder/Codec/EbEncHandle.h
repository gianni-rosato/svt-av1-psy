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

#include "EbResourceCoordinationResults.h"
#include "EbPictureDemuxResults.h"
#include "EbRateControlResults.h"

/**************************************
 * Component Private Data
 **************************************/
typedef struct EbEncHandle_s
{
    // Encode Instances & Compute Segments
    uint32_t                                  encodeInstanceTotalCount;
    uint32_t                                 *compute_segments_total_count_array;

    // Config Set Counts
    uint32_t                                  sequenceControlSetPoolTotalCount;

    // Full Results Count
    uint32_t                                  pictureControlSetPoolTotalCount;

    // Picture Buffer Counts
    uint32_t                                  reconPicturePoolTotalCount;
    uint32_t                                  coeffPicturePoolTotalCount;
    uint32_t                                  referencePicturePoolTotalCount;
    uint32_t                                  paReferencePicturPooleBufferInitCount;

    // Config Set Pool & Active Array
    EbSystemResource                     *sequenceControlSetPoolPtr;
    EbFifo                              **sequenceControlSetPoolProducerFifoPtrArray;
    EbSequenceControlSetInstance        **sequence_control_set_instance_array;

    // Full Results
    EbSystemResource                    **pictureControlSetPoolPtrArray;
    EbFifo                             ***pictureControlSetPoolProducerFifoPtrDblArray;

    //ParentControlSet
    EbSystemResource                    **pictureParentControlSetPoolPtrArray;
    EbFifo                             ***pictureParentControlSetPoolProducerFifoPtrDblArray;

    // Picture Buffers
    EbSystemResource                    **referencePicturePoolPtrArray;
    EbSystemResource                    **paReferencePicturePoolPtrArray;

    // Picture Buffer Producer Fifos
    EbFifo                             ***referencePicturePoolProducerFifoPtrDblArray;
    EbFifo                             ***paReferencePicturePoolProducerFifoPtrDblArray;

    // Thread Handles
    EbHandle                               resourceCoordinationThreadHandle;
    EbHandle                               pictureEnhancementThreadHandle;
    EbHandle                              *pictureAnalysisThreadHandleArray;
    EbHandle                               pictureDecisionThreadHandle;
    EbHandle                              *motionEstimationThreadHandleArray;
    EbHandle                               initialRateControlThreadHandle;
    EbHandle                              *sourceBasedOperationsThreadHandleArray;
    EbHandle                               pictureManagerThreadHandle;
    EbHandle                               rateControlThreadHandle;
    EbHandle                              *modeDecisionConfigurationThreadHandleArray;
    EbHandle                              *encDecThreadHandleArray;
    EbHandle                              *entropyCodingThreadHandleArray;
    EbHandle                              *dlfThreadHandleArray;
    EbHandle                              *cdefThreadHandleArray;
    EbHandle                              *restThreadHandleArray;

    EbHandle                               packetizationThreadHandle;

    // Contexts
    EbPtr                                  resourceCoordinationContextPtr;
    EbPtr                                  pictureEnhancementContextPtr;
    EbPtr                                 *pictureAnalysisContextPtrArray;
    EbPtr                                  pictureDecisionContextPtr;
    EbPtr                                 *motionEstimationContextPtrArray;
    EbPtr                                  initialRateControlContextPtr;
    EbPtr                                 *sourceBasedOperationsContextPtrArray;
    EbPtr                                  pictureManagerContextPtr;
    EbPtr                                  rateControlContextPtr;
    EbPtr                                 *modeDecisionConfigurationContextPtrArray;
    EbPtr                                 *encDecContextPtrArray;
    EbPtr                                 *entropyCodingContextPtrArray;
    EbPtr                                 *dlfContextPtrArray;
    EbPtr                                 *cdefContextPtrArray;
    EbPtr                                 *restContextPtrArray;
    EbPtr                                  packetizationContextPtr;

    // System Resource Managers
    EbSystemResource                     *input_buffer_resource_ptr;
    EbSystemResource                    **output_stream_buffer_resource_ptr_array;
    EbSystemResource                    **output_recon_buffer_resource_ptr_array;
    EbSystemResource                    **output_statistics_buffer_resource_ptr_array;
    EbSystemResource                     *resourceCoordinationResultsResourcePtr;
    EbSystemResource                     *pictureAnalysisResultsResourcePtr;
    EbSystemResource                     *pictureDecisionResultsResourcePtr;
    EbSystemResource                     *motionEstimationResultsResourcePtr;
    EbSystemResource                     *initialRateControlResultsResourcePtr;
    EbSystemResource                     *pictureDemuxResultsResourcePtr;
    EbSystemResource                     *rateControlTasksResourcePtr;
    EbSystemResource                     *rateControlResultsResourcePtr;
    EbSystemResource                     *encDecTasksResourcePtr;
    EbSystemResource                     *encDecResultsResourcePtr;
    EbSystemResource                     *entropyCodingResultsResourcePtr;
    EbSystemResource                     *dlfResultsResourcePtr;
    EbSystemResource                     *cdefResultsResourcePtr;
    EbSystemResource                     *restResultsResourcePtr;

    // Inter-Process Producer Fifos
    EbFifo                              **input_buffer_producer_fifo_ptr_array;
    EbFifo                             ***output_stream_buffer_producer_fifo_ptr_dbl_array;
    EbFifo                             ***output_recon_buffer_producer_fifo_ptr_dbl_array;
    EbFifo                             ***output_statistics_buffer_producer_fifo_ptr_dbl_array;
    EbFifo                              **resourceCoordinationResultsProducerFifoPtrArray;
    EbFifo                              **pictureAnalysisResultsProducerFifoPtrArray;
    EbFifo                              **pictureDecisionResultsProducerFifoPtrArray;
    EbFifo                              **motionEstimationResultsProducerFifoPtrArray;
    EbFifo                              **initialRateControlResultsProducerFifoPtrArray;
    EbFifo                              **pictureDemuxResultsProducerFifoPtrArray;
    EbFifo                              **pictureManagerResultsProducerFifoPtrArray;
    EbFifo                              **rateControlTasksProducerFifoPtrArray;
    EbFifo                              **rateControlResultsProducerFifoPtrArray;
    EbFifo                              **encDecTasksProducerFifoPtrArray;
    EbFifo                              **encDecResultsProducerFifoPtrArray;
    EbFifo                              **entropyCodingResultsProducerFifoPtrArray;
    EbFifo                              **dlfResultsProducerFifoPtrArray;
    EbFifo                              **cdefResultsProducerFifoPtrArray;
    EbFifo                              **restResultsProducerFifoPtrArray;


    // Inter-Process Consumer Fifos
    EbFifo                              **input_buffer_consumer_fifo_ptr_array;
    EbFifo                             ***output_stream_buffer_consumer_fifo_ptr_dbl_array;
    EbFifo                             ***output_recon_buffer_consumer_fifo_ptr_dbl_array;
    EbFifo                             ***output_statistics_buffer_consumer_fifo_ptr_dbl_array;
    EbFifo                              **resourceCoordinationResultsConsumerFifoPtrArray;
    EbFifo                              **pictureAnalysisResultsConsumerFifoPtrArray;
    EbFifo                              **pictureDecisionResultsConsumerFifoPtrArray;
    EbFifo                              **motionEstimationResultsConsumerFifoPtrArray;
    EbFifo                              **initialRateControlResultsConsumerFifoPtrArray;
    EbFifo                              **pictureDemuxResultsConsumerFifoPtrArray;
    EbFifo                              **rateControlTasksConsumerFifoPtrArray;
    EbFifo                              **rateControlResultsConsumerFifoPtrArray;
    EbFifo                              **encDecTasksConsumerFifoPtrArray;
    EbFifo                              **encDecResultsConsumerFifoPtrArray;
    EbFifo                              **entropyCodingResultsConsumerFifoPtrArray;
    EbFifo                              **dlfResultsConsumerFifoPtrArray;
    EbFifo                              **cdefResultsConsumerFifoPtrArray;
    EbFifo                              **restResultsConsumerFifoPtrArray;

    // Callbacks
    EbCallback_t                          **app_callback_ptr_array;

    // Memory Map
    EbMemoryMapEntry                       *memory_map;
    uint32_t                                memory_map_index;
    uint64_t                                total_lib_memory;

} EbEncHandle_t;


#endif // EbEncHandle_h