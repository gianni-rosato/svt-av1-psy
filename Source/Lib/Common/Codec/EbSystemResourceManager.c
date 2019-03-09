/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>

#include "EbSystemResourceManager.h"

/**************************************
 * EbFifoCtor
 **************************************/
static EbErrorType EbFifoCtor(
    EbFifo_t           *fifoPtr,
    uint32_t              initial_count,
    uint32_t              max_count,
    EbObjectWrapper_t  *firstWrapperPtr,
    EbObjectWrapper_t  *lastWrapperPtr,
    EbMuxingQueue_t    *queuePtr)
{
    // Create Counting Semaphore
    EB_CREATESEMAPHORE(EbHandle, fifoPtr->countingSemaphore, sizeof(EbHandle), EB_SEMAPHORE, initial_count, max_count);

    // Create Buffer Pool Mutex
    EB_CREATEMUTEX(EbHandle, fifoPtr->lockoutMutex, sizeof(EbHandle), EB_MUTEX);

    // Initialize Fifo First & Last ptrs
    fifoPtr->firstPtr = firstWrapperPtr;
    fifoPtr->lastPtr = lastWrapperPtr;

    // Copy the Muxing Queue ptr this Fifo belongs to
    fifoPtr->queuePtr = queuePtr;

    return EB_ErrorNone;
}


/**************************************
 * EbFifoPushBack
 **************************************/
static EbErrorType EbFifoPushBack(
    EbFifo_t            *fifoPtr,
    EbObjectWrapper_t   *wrapper_ptr)
{
    EbErrorType return_error = EB_ErrorNone;

    // If FIFO is empty
    if (fifoPtr->firstPtr == (EbObjectWrapper_t*)EB_NULL) {
        fifoPtr->firstPtr = wrapper_ptr;
        fifoPtr->lastPtr = wrapper_ptr;
    }
    else {
        fifoPtr->lastPtr->nextPtr = wrapper_ptr;
        fifoPtr->lastPtr = wrapper_ptr;
    }

    fifoPtr->lastPtr->nextPtr = (EbObjectWrapper_t*)EB_NULL;

    return return_error;
}

/**************************************
 * EbFifoPopFront
 **************************************/
static EbErrorType EbFifoPopFront(
    EbFifo_t            *fifoPtr,
    EbObjectWrapper_t  **wrapper_ptr)
{
    EbErrorType return_error = EB_ErrorNone;

    // Set wrapper_ptr to head of BufferPool
    *wrapper_ptr = fifoPtr->firstPtr;

    // Update tail of BufferPool if the BufferPool is now empty
    fifoPtr->lastPtr = (fifoPtr->firstPtr == fifoPtr->lastPtr) ? (EbObjectWrapper_t*)EB_NULL : fifoPtr->lastPtr;

    // Update head of BufferPool
    fifoPtr->firstPtr = fifoPtr->firstPtr->nextPtr;

    return return_error;
}

/**************************************
 * EbCircularBufferCtor
 **************************************/
static EbErrorType EbCircularBufferCtor(
    EbCircularBuffer_t  **buffer_dbl_ptr,
    uint32_t                buffer_total_count)
{
    uint32_t bufferIndex;
    EbCircularBuffer_t *bufferPtr;

    EB_MALLOC(EbCircularBuffer_t*, bufferPtr, sizeof(EbCircularBuffer_t), EB_N_PTR);

    *buffer_dbl_ptr = bufferPtr;

    bufferPtr->buffer_total_count = buffer_total_count;

    EB_MALLOC(EbPtr*, bufferPtr->arrayPtr, sizeof(EbPtr) * bufferPtr->buffer_total_count, EB_N_PTR);

    for (bufferIndex = 0; bufferIndex < bufferPtr->buffer_total_count; ++bufferIndex) {
        bufferPtr->arrayPtr[bufferIndex] = EB_NULL;
    }

    bufferPtr->headIndex = 0;
    bufferPtr->tailIndex = 0;

    bufferPtr->currentCount = 0;

    return EB_ErrorNone;
}



/**************************************
 * EbCircularBufferEmptyCheck
 **************************************/
static EbBool EbCircularBufferEmptyCheck(
    EbCircularBuffer_t   *bufferPtr)
{
    return ((bufferPtr->headIndex == bufferPtr->tailIndex) && (bufferPtr->arrayPtr[bufferPtr->headIndex] == EB_NULL)) ? EB_TRUE : EB_FALSE;
}

/**************************************
 * EbCircularBufferPopFront
 **************************************/
static EbErrorType EbCircularBufferPopFront(
    EbCircularBuffer_t   *bufferPtr,
    EbPtr               *object_ptr)
{
    EbErrorType return_error = EB_ErrorNone;

    // Copy the head of the buffer into the object_ptr
    *object_ptr = bufferPtr->arrayPtr[bufferPtr->headIndex];
    bufferPtr->arrayPtr[bufferPtr->headIndex] = EB_NULL;

    // Increment the head & check for rollover
    bufferPtr->headIndex = (bufferPtr->headIndex == bufferPtr->buffer_total_count - 1) ? 0 : bufferPtr->headIndex + 1;

    // Decrement the Current Count
    --bufferPtr->currentCount;

    return return_error;
}

/**************************************
 * EbCircularBufferPushBack
 **************************************/
static EbErrorType EbCircularBufferPushBack(
    EbCircularBuffer_t   *bufferPtr,
    EbPtr                object_ptr)
{
    EbErrorType return_error = EB_ErrorNone;

    // Copy the pointer into the array
    bufferPtr->arrayPtr[bufferPtr->tailIndex] = object_ptr;

    // Increment the tail & check for rollover
    bufferPtr->tailIndex = (bufferPtr->tailIndex == bufferPtr->buffer_total_count - 1) ? 0 : bufferPtr->tailIndex + 1;

    // Increment the Current Count
    ++bufferPtr->currentCount;

    return return_error;
}

/**************************************
 * EbCircularBufferPushFront
 **************************************/
static EbErrorType EbCircularBufferPushFront(
    EbCircularBuffer_t   *bufferPtr,
    EbPtr                object_ptr)
{
    EbErrorType return_error = EB_ErrorNone;

    // Decrement the headIndex
    bufferPtr->headIndex = (bufferPtr->headIndex == 0) ? bufferPtr->buffer_total_count - 1 : bufferPtr->headIndex - 1;

    // Copy the pointer into the array
    bufferPtr->arrayPtr[bufferPtr->headIndex] = object_ptr;

    // Increment the Current Count
    ++bufferPtr->currentCount;

    return return_error;
}

/**************************************
 * EbMuxingQueueCtor
 **************************************/
static EbErrorType EbMuxingQueueCtor(
    EbMuxingQueue_t   **queueDblPtr,
    uint32_t              object_total_count,
    uint32_t              processTotalCount,
    EbFifo_t         ***processFifoPtrArrayPtr)
{
    EbMuxingQueue_t *queuePtr;
    uint32_t processIndex;
    EbErrorType     return_error = EB_ErrorNone;

    EB_MALLOC(EbMuxingQueue_t *, queuePtr, sizeof(EbMuxingQueue_t), EB_N_PTR);
    *queueDblPtr = queuePtr;

    queuePtr->processTotalCount = processTotalCount;

    // Lockout Mutex
    EB_CREATEMUTEX(EbHandle, queuePtr->lockoutMutex, sizeof(EbHandle), EB_MUTEX);

    // Construct Object Circular Buffer
    return_error = EbCircularBufferCtor(
        &queuePtr->objectQueue,
        object_total_count);
    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }
    // Construct Process Circular Buffer
    return_error = EbCircularBufferCtor(
        &queuePtr->processQueue,
        queuePtr->processTotalCount);
    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }
    // Construct the Process Fifos
    EB_MALLOC(EbFifo_t**, queuePtr->processFifoPtrArray, sizeof(EbFifo_t*) * queuePtr->processTotalCount, EB_N_PTR);

    for (processIndex = 0; processIndex < queuePtr->processTotalCount; ++processIndex) {
        EB_MALLOC(EbFifo_t*, queuePtr->processFifoPtrArray[processIndex], sizeof(EbFifo_t) * queuePtr->processTotalCount, EB_N_PTR);
        return_error = EbFifoCtor(
            queuePtr->processFifoPtrArray[processIndex],
            0,
            object_total_count,
            (EbObjectWrapper_t *)EB_NULL,
            (EbObjectWrapper_t *)EB_NULL,
            queuePtr);
        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }

    *processFifoPtrArrayPtr = queuePtr->processFifoPtrArray;

    return return_error;
}

/**************************************
 * EbMuxingQueueAssignation
 **************************************/
static EbErrorType EbMuxingQueueAssignation(
    EbMuxingQueue_t *queuePtr)
{
    EbErrorType return_error = EB_ErrorNone;
    EbFifo_t *processFifoPtr;
    EbObjectWrapper_t *wrapper_ptr;

    // while loop
    while ((EbCircularBufferEmptyCheck(queuePtr->objectQueue) == EB_FALSE) &&
        (EbCircularBufferEmptyCheck(queuePtr->processQueue) == EB_FALSE)) {
        // Get the next process
        EbCircularBufferPopFront(
            queuePtr->processQueue,
            (void **)&processFifoPtr);

        // Get the next object
        EbCircularBufferPopFront(
            queuePtr->objectQueue,
            (void **)&wrapper_ptr);

        // Block on the Process Fifo's Mutex
        eb_block_on_mutex(processFifoPtr->lockoutMutex);

        // Put the object on the fifo
        EbFifoPushBack(
            processFifoPtr,
            wrapper_ptr);

        // Release the Process Fifo's Mutex
        eb_release_mutex(processFifoPtr->lockoutMutex);

        // Post the semaphore
        eb_post_semaphore(processFifoPtr->countingSemaphore);
    }

    return return_error;
}

/**************************************
 * EbMuxingQueueObjectPushBack
 **************************************/
static EbErrorType EbMuxingQueueObjectPushBack(
    EbMuxingQueue_t    *queuePtr,
    EbObjectWrapper_t  *object_ptr)
{
    EbErrorType return_error = EB_ErrorNone;

    EbCircularBufferPushBack(
        queuePtr->objectQueue,
        object_ptr);

    EbMuxingQueueAssignation(queuePtr);

    return return_error;
}

/**************************************
* EbMuxingQueueObjectPushFront
**************************************/
static EbErrorType EbMuxingQueueObjectPushFront(
    EbMuxingQueue_t    *queuePtr,
    EbObjectWrapper_t  *object_ptr)
{
    EbErrorType return_error = EB_ErrorNone;

    EbCircularBufferPushFront(
        queuePtr->objectQueue,
        object_ptr);

    EbMuxingQueueAssignation(queuePtr);

    return return_error;
}

/*********************************************************************
 * eb_object_release_enable
 *   Enables the releaseEnable member of EbObjectWrapper.  Used by
 *   certain objects (e.g. SequenceControlSet) to control whether
 *   EbObjectWrappers are allowed to be released or not.
 *
 *   resource_ptr
 *      pointer to the SystemResource that manages the EbObjectWrapper.
 *      The emptyFifo's lockoutMutex is used to write protect the
 *      modification of the EbObjectWrapper.
 *
 *   wrapper_ptr
 *      pointer to the EbObjectWrapper to be modified.
 *********************************************************************/
EbErrorType eb_object_release_enable(
    EbObjectWrapper_t   *wrapper_ptr)
{
    EbErrorType return_error = EB_ErrorNone;

    eb_block_on_mutex(wrapper_ptr->systemResourcePtr->emptyQueue->lockoutMutex);

    wrapper_ptr->releaseEnable = EB_TRUE;

    eb_release_mutex(wrapper_ptr->systemResourcePtr->emptyQueue->lockoutMutex);

    return return_error;
}

/*********************************************************************
 * eb_object_release_disable
 *   Disables the releaseEnable member of EbObjectWrapper.  Used by
 *   certain objects (e.g. SequenceControlSet) to control whether
 *   EbObjectWrappers are allowed to be released or not.
 *
 *   resource_ptr
 *      pointer to the SystemResource that manages the EbObjectWrapper.
 *      The emptyFifo's lockoutMutex is used to write protect the
 *      modification of the EbObjectWrapper.
 *
 *   wrapper_ptr
 *      pointer to the EbObjectWrapper to be modified.
 *********************************************************************/
EbErrorType eb_object_release_disable(
    EbObjectWrapper_t   *wrapper_ptr)
{
    EbErrorType return_error = EB_ErrorNone;

    eb_block_on_mutex(wrapper_ptr->systemResourcePtr->emptyQueue->lockoutMutex);

    wrapper_ptr->releaseEnable = EB_FALSE;

    eb_release_mutex(wrapper_ptr->systemResourcePtr->emptyQueue->lockoutMutex);

    return return_error;
}

/*********************************************************************
 * eb_object_inc_live_count
 *   Increments the liveCount member of EbObjectWrapper.  Used by
 *   certain objects (e.g. SequenceControlSet) to count the number of active
 *   pointers of a EbObjectWrapper in pipeline at any point in time.
 *
 *   resource_ptr
 *      pointer to the SystemResource that manages the EbObjectWrapper.
 *      The emptyFifo's lockoutMutex is used to write protect the
 *      modification of the EbObjectWrapper.
 *
 *   wrapper_ptr
 *      pointer to the EbObjectWrapper to be modified.
 *********************************************************************/
EbErrorType eb_object_inc_live_count(
    EbObjectWrapper_t   *wrapper_ptr,
    uint32_t               increment_number)
{
    EbErrorType return_error = EB_ErrorNone;

    eb_block_on_mutex(wrapper_ptr->systemResourcePtr->emptyQueue->lockoutMutex);

    wrapper_ptr->liveCount += increment_number;

    eb_release_mutex(wrapper_ptr->systemResourcePtr->emptyQueue->lockoutMutex);

    return return_error;
}

/*********************************************************************
 * eb_system_resource_ctor
 *   Constructor for EbSystemResource.  Fully constructs all members
 *   of EbSystemResource including the object with the passed
 *   object_ctor function.
 *
 *   resource_ptr
 *     pointer that will contain the SystemResource to be constructed.
 *
 *   object_total_count
 *     Number of objects to be managed by the SystemResource.
 *
 *   full_fifo_enabled
 *     Bool that describes if the SystemResource is to have an output
 *     fifo.  An outputFifo is not used by certain objects (e.g.
 *     SequenceControlSet).
 *
 *   object_ctor
 *     Function pointer to the constructor of the object managed by
 *     SystemResource referenced by resource_ptr. No object level
 *     construction is performed if object_ctor is NULL.
 *
 *   object_init_data_ptr

 *     pointer to data block to be used during the construction of
 *     the object. object_init_data_ptr is passed to object_ctor when
 *     object_ctor is called.
 *********************************************************************/
EbErrorType eb_system_resource_ctor(
    EbSystemResource_t **resource_dbl_ptr,
    uint32_t               object_total_count,
    uint32_t               producer_process_total_count,
    uint32_t               consumer_process_total_count,
    EbFifo_t          ***producer_fifo_ptr_array_ptr,
    EbFifo_t          ***consumer_fifo_ptr_array_ptr,
    EbBool              full_fifo_enabled,
    EB_CTOR              object_ctor,
    EbPtr               object_init_data_ptr)
{
    uint32_t wrapperIndex;
    EbErrorType return_error = EB_ErrorNone;
    // Allocate the System Resource
    EbSystemResource_t *resource_ptr;

    EB_MALLOC(EbSystemResource_t*, resource_ptr, sizeof(EbSystemResource_t), EB_N_PTR);
    *resource_dbl_ptr = resource_ptr;

    resource_ptr->object_total_count = object_total_count;

    // Allocate array for wrapper pointers
    EB_MALLOC(EbObjectWrapper_t**, resource_ptr->wrapperPtrPool, sizeof(EbObjectWrapper_t*) * resource_ptr->object_total_count, EB_N_PTR);

    // Initialize each wrapper
    for (wrapperIndex = 0; wrapperIndex < resource_ptr->object_total_count; ++wrapperIndex) {
        EB_MALLOC(EbObjectWrapper_t*, resource_ptr->wrapperPtrPool[wrapperIndex], sizeof(EbObjectWrapper_t), EB_N_PTR);
        resource_ptr->wrapperPtrPool[wrapperIndex]->liveCount = 0;
        resource_ptr->wrapperPtrPool[wrapperIndex]->releaseEnable = EB_TRUE;
        resource_ptr->wrapperPtrPool[wrapperIndex]->systemResourcePtr = resource_ptr;

        // Call the Constructor for each element
        if (object_ctor) {
            return_error = object_ctor(
                &resource_ptr->wrapperPtrPool[wrapperIndex]->object_ptr,
                object_init_data_ptr);
            if (return_error == EB_ErrorInsufficientResources) {
                return EB_ErrorInsufficientResources;
            }
        }
    }

    // Initialize the Empty Queue
    return_error = EbMuxingQueueCtor(
        &resource_ptr->emptyQueue,
        resource_ptr->object_total_count,
        producer_process_total_count,
        producer_fifo_ptr_array_ptr);
    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }
    // Fill the Empty Fifo with every ObjectWrapper
    for (wrapperIndex = 0; wrapperIndex < resource_ptr->object_total_count; ++wrapperIndex) {
        EbMuxingQueueObjectPushBack(
            resource_ptr->emptyQueue,
            resource_ptr->wrapperPtrPool[wrapperIndex]);
    }

    // Initialize the Full Queue
    if (full_fifo_enabled == EB_TRUE) {
        return_error = EbMuxingQueueCtor(
            &resource_ptr->fullQueue,
            resource_ptr->object_total_count,
            consumer_process_total_count,
            consumer_fifo_ptr_array_ptr);
        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }
    else {
        resource_ptr->fullQueue = (EbMuxingQueue_t *)EB_NULL;
        consumer_fifo_ptr_array_ptr = (EbFifo_t ***)EB_NULL;
    }

    return return_error;
}



/*********************************************************************
 * EbSystemResourceReleaseProcess
 *********************************************************************/
static EbErrorType EbReleaseProcess(
    EbFifo_t   *processFifoPtr)
{
    EbErrorType return_error = EB_ErrorNone;

    eb_block_on_mutex(processFifoPtr->queuePtr->lockoutMutex);

    EbCircularBufferPushFront(
        processFifoPtr->queuePtr->processQueue,
        processFifoPtr);

    EbMuxingQueueAssignation(processFifoPtr->queuePtr);

    eb_release_mutex(processFifoPtr->queuePtr->lockoutMutex);

    return return_error;
}

/*********************************************************************
 * EbSystemResourcePostObject
 *   Queues a full EbObjectWrapper to the SystemResource. This
 *   function posts the SystemResource fullFifo countingSemaphore.
 *   This function is write protected by the SystemResource fullFifo
 *   lockoutMutex.
 *
 *   resource_ptr
 *      pointer to the SystemResource that the EbObjectWrapper is
 *      posted to.
 *
 *   wrapper_ptr
 *      pointer to EbObjectWrapper to be posted.
 *********************************************************************/
EbErrorType eb_post_full_object(
    EbObjectWrapper_t   *object_ptr)
{
    EbErrorType return_error = EB_ErrorNone;

    eb_block_on_mutex(object_ptr->systemResourcePtr->fullQueue->lockoutMutex);

    EbMuxingQueueObjectPushBack(
        object_ptr->systemResourcePtr->fullQueue,
        object_ptr);

    eb_release_mutex(object_ptr->systemResourcePtr->fullQueue->lockoutMutex);

    return return_error;
}

/*********************************************************************
 * EbSystemResourceReleaseObject
 *   Queues an empty EbObjectWrapper to the SystemResource. This
 *   function posts the SystemResource emptyFifo countingSemaphore.
 *   This function is write protected by the SystemResource emptyFifo
 *   lockoutMutex.
 *
 *   object_ptr
 *      pointer to EbObjectWrapper to be released.
 *********************************************************************/
EbErrorType eb_release_object(
    EbObjectWrapper_t   *object_ptr)
{
    EbErrorType return_error = EB_ErrorNone;

    eb_block_on_mutex(object_ptr->systemResourcePtr->emptyQueue->lockoutMutex);

    // Decrement liveCount
    object_ptr->liveCount = (object_ptr->liveCount == 0) ? object_ptr->liveCount : object_ptr->liveCount - 1;

    if ((object_ptr->releaseEnable == EB_TRUE) && (object_ptr->liveCount == 0)) {

        // Set liveCount to EB_ObjectWrapperReleasedValue
        object_ptr->liveCount = EB_ObjectWrapperReleasedValue;

        EbMuxingQueueObjectPushFront(
            object_ptr->systemResourcePtr->emptyQueue,
            object_ptr);

    }

    eb_release_mutex(object_ptr->systemResourcePtr->emptyQueue->lockoutMutex);

    return return_error;
}

/*********************************************************************
 * EbSystemResourceGetEmptyObject
 *   Dequeues an empty EbObjectWrapper from the SystemResource.  This
 *   function blocks on the SystemResource emptyFifo countingSemaphore.
 *   This function is write protected by the SystemResource emptyFifo
 *   lockoutMutex.
 *
 *   resource_ptr
 *      pointer to the SystemResource that provides the empty
 *      EbObjectWrapper.
 *
 *   wrapper_dbl_ptr
 *      Double pointer used to pass the pointer to the empty
 *      EbObjectWrapper pointer.
 *********************************************************************/
EbErrorType eb_get_empty_object(
    EbFifo_t   *empty_fifo_ptr,
    EbObjectWrapper_t **wrapper_dbl_ptr)
{
    EbErrorType return_error = EB_ErrorNone;

    // Queue the Fifo requesting the empty fifo
    EbReleaseProcess(empty_fifo_ptr);

    // Block on the counting Semaphore until an empty buffer is available
    eb_block_on_semaphore(empty_fifo_ptr->countingSemaphore);

    // Acquire lockout Mutex
    eb_block_on_mutex(empty_fifo_ptr->lockoutMutex);

    // Get the empty object
    EbFifoPopFront(
        empty_fifo_ptr,
        wrapper_dbl_ptr);

    // Reset the wrapper's liveCount
    (*wrapper_dbl_ptr)->liveCount = 0;

    // Object release enable
    (*wrapper_dbl_ptr)->releaseEnable = EB_TRUE;

    // Release Mutex
    eb_release_mutex(empty_fifo_ptr->lockoutMutex);

    return return_error;
}

/*********************************************************************
 * EbSystemResourceGetFullObject
 *   Dequeues an full EbObjectWrapper from the SystemResource. This
 *   function blocks on the SystemResource fullFifo countingSemaphore.
 *   This function is write protected by the SystemResource fullFifo
 *   lockoutMutex.
 *
 *   resource_ptr
 *      pointer to the SystemResource that provides the full
 *      EbObjectWrapper.
 *
 *   wrapper_dbl_ptr
 *      Double pointer used to pass the pointer to the full
 *      EbObjectWrapper pointer.
 *********************************************************************/
EbErrorType eb_get_full_object(
    EbFifo_t   *full_fifo_ptr,
    EbObjectWrapper_t **wrapper_dbl_ptr)
{
    EbErrorType return_error = EB_ErrorNone;

    // Queue the Fifo requesting the full fifo
    EbReleaseProcess(full_fifo_ptr);

    // Block on the counting Semaphore until an empty buffer is available
    eb_block_on_semaphore(full_fifo_ptr->countingSemaphore);

    // Acquire lockout Mutex
    eb_block_on_mutex(full_fifo_ptr->lockoutMutex);

    EbFifoPopFront(
        full_fifo_ptr,
        wrapper_dbl_ptr);

    // Release Mutex
    eb_release_mutex(full_fifo_ptr->lockoutMutex);

    return return_error;
}

/**************************************
* EbFifoPopFront
**************************************/
static EbBool EbFifoPeakFront(
    EbFifo_t            *fifoPtr)
{

    // Set wrapper_ptr to head of BufferPool
    if (fifoPtr->firstPtr == (EbObjectWrapper_t*)EB_NULL) {
        return EB_TRUE;
    }
    else {
        return EB_FALSE;
    }
}


EbErrorType eb_get_full_object_non_blocking(
    EbFifo_t   *full_fifo_ptr,
    EbObjectWrapper_t **wrapper_dbl_ptr)
{
    EbErrorType return_error = EB_ErrorNone;
    EbBool      fifoEmpty;
    // Queue the Fifo requesting the full fifo
    EbReleaseProcess(full_fifo_ptr);

    // Acquire lockout Mutex
    eb_block_on_mutex(full_fifo_ptr->lockoutMutex);

    fifoEmpty = EbFifoPeakFront(
        full_fifo_ptr);

    // Release Mutex
    eb_release_mutex(full_fifo_ptr->lockoutMutex);

    if (fifoEmpty == EB_FALSE)
        eb_get_full_object(
            full_fifo_ptr,
            wrapper_dbl_ptr);
    else
        *wrapper_dbl_ptr = (EbObjectWrapper_t*)EB_NULL;

    return return_error;
}