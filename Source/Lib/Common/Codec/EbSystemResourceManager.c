// clang-format off
/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>

#include "EbSystemResourceManager.h"
#include "EbDefinitions.h"
#include "EbThreads.h"

void EbFifoDctor(EbPtr p)
{
    EbFifo *obj = (EbFifo*)p;
    EB_DESTROY_SEMAPHORE(obj->counting_semaphore);
    EB_DESTROY_MUTEX(obj->lockout_mutex);
}
/**************************************
 * EbFifoCtor
 **************************************/
static EbErrorType EbFifoCtor(
    EbFifo           *fifoPtr,
    uint32_t              initial_count,
    uint32_t              max_count,
    EbObjectWrapper  *firstWrapperPtr,
    EbObjectWrapper  *lastWrapperPtr,
    EbMuxingQueue    *queue_ptr)
{
    fifoPtr->dctor = EbFifoDctor;
    // Create Counting Semaphore
    EB_CREATE_SEMAPHORE(fifoPtr->counting_semaphore, initial_count, max_count);

    // Create Buffer Pool Mutex
    EB_CREATE_MUTEX(fifoPtr->lockout_mutex);

    // Initialize Fifo First & Last ptrs
    fifoPtr->first_ptr = firstWrapperPtr;
    fifoPtr->last_ptr = lastWrapperPtr;

    // Copy the Muxing Queue ptr this Fifo belongs to
    fifoPtr->queue_ptr = queue_ptr;

    return EB_ErrorNone;
}

/**************************************
 * EbFifoPushBack
 **************************************/
static EbErrorType EbFifoPushBack(
    EbFifo            *fifoPtr,
    EbObjectWrapper   *wrapper_ptr)
{
    EbErrorType return_error = EB_ErrorNone;

    // If FIFO is empty
    if (fifoPtr->first_ptr == (EbObjectWrapper*)EB_NULL) {
        fifoPtr->first_ptr = wrapper_ptr;
        fifoPtr->last_ptr = wrapper_ptr;
    }
    else {
        fifoPtr->last_ptr->next_ptr = wrapper_ptr;
        fifoPtr->last_ptr = wrapper_ptr;
    }

    fifoPtr->last_ptr->next_ptr = (EbObjectWrapper*)EB_NULL;

    return return_error;
}

/**************************************
 * EbFifoPopFront
 **************************************/
static EbErrorType EbFifoPopFront(
    EbFifo            *fifoPtr,
    EbObjectWrapper  **wrapper_ptr)
{
    EbErrorType return_error = EB_ErrorNone;

    // Set wrapper_ptr to head of BufferPool
    *wrapper_ptr = fifoPtr->first_ptr;

    // Update tail of BufferPool if the BufferPool is now empty
    fifoPtr->last_ptr = (fifoPtr->first_ptr == fifoPtr->last_ptr) ? (EbObjectWrapper*)EB_NULL : fifoPtr->last_ptr;

    // Update head of BufferPool
    fifoPtr->first_ptr = fifoPtr->first_ptr->next_ptr;

    return return_error;
}

void EbCircularBufferDctor(EbPtr p)
{
    EbCircularBuffer* obj = (EbCircularBuffer*)p;
    EB_FREE(obj->array_ptr);
}

/**************************************
 * EbCircularBufferCtor
 **************************************/
static EbErrorType EbCircularBufferCtor(
    EbCircularBuffer  *bufferPtr,
    uint32_t                buffer_total_count)
{
    bufferPtr->dctor = EbCircularBufferDctor;

    bufferPtr->buffer_total_count = buffer_total_count;

    EB_CALLOC(bufferPtr->array_ptr, bufferPtr->buffer_total_count, sizeof(EbPtr));

    return EB_ErrorNone;
}

/**************************************
 * EbCircularBufferEmptyCheck
 **************************************/
static EbBool EbCircularBufferEmptyCheck(
    EbCircularBuffer   *bufferPtr)
{
    return ((bufferPtr->head_index == bufferPtr->tail_index) && (bufferPtr->array_ptr[bufferPtr->head_index] == EB_NULL)) ? EB_TRUE : EB_FALSE;
}

/**************************************
 * EbCircularBufferPopFront
 **************************************/
static EbErrorType EbCircularBufferPopFront(
    EbCircularBuffer   *bufferPtr,
    EbPtr               *object_ptr)
{
    EbErrorType return_error = EB_ErrorNone;

    // Copy the head of the buffer into the object_ptr
    *object_ptr = bufferPtr->array_ptr[bufferPtr->head_index];
    bufferPtr->array_ptr[bufferPtr->head_index] = EB_NULL;

    // Increment the head & check for rollover
    bufferPtr->head_index = (bufferPtr->head_index == bufferPtr->buffer_total_count - 1) ? 0 : bufferPtr->head_index + 1;

    // Decrement the Current Count
    --bufferPtr->current_count;

    return return_error;
}

/**************************************
 * EbCircularBufferPushBack
 **************************************/
static EbErrorType EbCircularBufferPushBack(
    EbCircularBuffer   *bufferPtr,
    EbPtr                object_ptr)
{
    EbErrorType return_error = EB_ErrorNone;

    // Copy the pointer into the array
    bufferPtr->array_ptr[bufferPtr->tail_index] = object_ptr;

    // Increment the tail & check for rollover
    bufferPtr->tail_index = (bufferPtr->tail_index == bufferPtr->buffer_total_count - 1) ? 0 : bufferPtr->tail_index + 1;

    // Increment the Current Count
    ++bufferPtr->current_count;

    return return_error;
}

/**************************************
 * EbCircularBufferPushFront
 **************************************/
static EbErrorType EbCircularBufferPushFront(
    EbCircularBuffer   *bufferPtr,
    EbPtr                object_ptr)
{
    EbErrorType return_error = EB_ErrorNone;

    // Decrement the head_index
    bufferPtr->head_index = (bufferPtr->head_index == 0) ? bufferPtr->buffer_total_count - 1 : bufferPtr->head_index - 1;

    // Copy the pointer into the array
    bufferPtr->array_ptr[bufferPtr->head_index] = object_ptr;

    // Increment the Current Count
    ++bufferPtr->current_count;

    return return_error;
}

void EbMuxingQueueDctor(EbPtr p)
{
    EbMuxingQueue* obj = (EbMuxingQueue*)p;
    EB_DELETE_PTR_ARRAY(obj->process_fifo_ptr_array, obj->process_total_count);
    EB_DELETE(obj->object_queue);
    EB_DELETE(obj->process_queue);
    EB_DESTROY_MUTEX(obj->lockout_mutex);
}

/**************************************
 * EbMuxingQueueCtor
 **************************************/
static EbErrorType EbMuxingQueueCtor(
    EbMuxingQueue        *queue_ptr,
    uint32_t              object_total_count,
    uint32_t              process_total_count)
{
    uint32_t processIndex;
    EbErrorType     return_error = EB_ErrorNone;

    queue_ptr->dctor = EbMuxingQueueDctor;
    queue_ptr->process_total_count = process_total_count;

    // Lockout Mutex
    EB_CREATE_MUTEX(queue_ptr->lockout_mutex);

    // Construct Object Circular Buffer
    EB_NEW(
        queue_ptr->object_queue,
        EbCircularBufferCtor,
        object_total_count);
    // Construct Process Circular Buffer
    EB_NEW(
        queue_ptr->process_queue,
        EbCircularBufferCtor,
        queue_ptr->process_total_count);
    // Construct the Process Fifos
    EB_ALLOC_PTR_ARRAY(queue_ptr->process_fifo_ptr_array, queue_ptr->process_total_count);

    for (processIndex = 0; processIndex < queue_ptr->process_total_count; ++processIndex) {
        EB_NEW(
            queue_ptr->process_fifo_ptr_array[processIndex],
            EbFifoCtor,
            0,
            object_total_count,
            (EbObjectWrapper *)EB_NULL,
            (EbObjectWrapper *)EB_NULL,
            queue_ptr);
    }

    return return_error;
}

/**************************************
 * EbMuxingQueueAssignation
 **************************************/
static EbErrorType EbMuxingQueueAssignation(
    EbMuxingQueue *queue_ptr)
{
    EbErrorType return_error = EB_ErrorNone;
    EbFifo *processFifoPtr;
    EbObjectWrapper *wrapper_ptr;

    // while loop
    while ((EbCircularBufferEmptyCheck(queue_ptr->object_queue) == EB_FALSE) &&
        (EbCircularBufferEmptyCheck(queue_ptr->process_queue) == EB_FALSE)) {
        // Get the next process
        EbCircularBufferPopFront(
            queue_ptr->process_queue,
            (void **)&processFifoPtr);

        // Get the next object
        EbCircularBufferPopFront(
            queue_ptr->object_queue,
            (void **)&wrapper_ptr);

        // Block on the Process Fifo's Mutex
        eb_block_on_mutex(processFifoPtr->lockout_mutex);

        // Put the object on the fifo
        EbFifoPushBack(
            processFifoPtr,
            wrapper_ptr);

        // Release the Process Fifo's Mutex
        eb_release_mutex(processFifoPtr->lockout_mutex);

        // Post the semaphore
        eb_post_semaphore(processFifoPtr->counting_semaphore);
    }

    return return_error;
}

static EbErrorType EbMuxingQueueAssignationNonBlocking(
    EbMuxingQueue *queue_ptr)
{
    EbErrorType return_error = EB_ErrorNone;
    EbFifo *processFifoPtr;
    EbObjectWrapper *wrapper_ptr;

    // while loop
    while ((EbCircularBufferEmptyCheck(queue_ptr->object_queue) == EB_FALSE) &&
        (EbCircularBufferEmptyCheck(queue_ptr->process_queue) == EB_FALSE))
    {
        // Get the next process
        EbCircularBufferPopFront(
            queue_ptr->process_queue,
            (void **)&processFifoPtr);

        // Get the next object
        EbCircularBufferPopFront(
            queue_ptr->object_queue,
            (void **)&wrapper_ptr);

        // Put the object on the fifo
        EbFifoPushBack(
            processFifoPtr,
            wrapper_ptr);

        // Post the semaphore
        eb_post_semaphore(processFifoPtr->counting_semaphore);
    }

    return return_error;
}

/**************************************
 * EbMuxingQueueObjectPushBack
 **************************************/
static EbErrorType EbMuxingQueueObjectPushBack(
    EbMuxingQueue    *queue_ptr,
    EbObjectWrapper  *object_ptr)
{
    EbErrorType return_error = EB_ErrorNone;

    EbCircularBufferPushBack(
        queue_ptr->object_queue,
        object_ptr);

    EbMuxingQueueAssignation(queue_ptr);

    return return_error;
}

/**************************************
* EbMuxingQueueObjectPushFront
**************************************/
static EbErrorType EbMuxingQueueObjectPushFront(
    EbMuxingQueue    *queue_ptr,
    EbObjectWrapper  *object_ptr)
{
    EbErrorType return_error = EB_ErrorNone;

    EbCircularBufferPushFront(
        queue_ptr->object_queue,
        object_ptr);

    EbMuxingQueueAssignation(queue_ptr);

    return return_error;
}

static EbFifo* eb_muxing_queue_get_fifo(
    EbMuxingQueue        *queue_ptr,
    uint32_t index)
{
    assert(queue_ptr->process_fifo_ptr_array && (queue_ptr->process_total_count > index));
    return queue_ptr->process_fifo_ptr_array[index];
}


/*********************************************************************
 * eb_object_release_enable
 *   Enables the release_enable member of EbObjectWrapper.  Used by
 *   certain objects (e.g. SequenceControlSet) to control whether
 *   EbObjectWrappers are allowed to be released or not.
 *
 *   resource_ptr
 *      pointer to the SystemResource that manages the EbObjectWrapper.
 *      The emptyFifo's lockout_mutex is used to write protect the
 *      modification of the EbObjectWrapper.
 *
 *   wrapper_ptr
 *      pointer to the EbObjectWrapper to be modified.
 *********************************************************************/
EbErrorType eb_object_release_enable(
    EbObjectWrapper   *wrapper_ptr)
{
    EbErrorType return_error = EB_ErrorNone;

    eb_block_on_mutex(wrapper_ptr->system_resource_ptr->empty_queue->lockout_mutex);

    wrapper_ptr->release_enable = EB_TRUE;

    eb_release_mutex(wrapper_ptr->system_resource_ptr->empty_queue->lockout_mutex);

    return return_error;
}

/*********************************************************************
 * eb_object_release_disable
 *   Disables the release_enable member of EbObjectWrapper.  Used by
 *   certain objects (e.g. SequenceControlSet) to control whether
 *   EbObjectWrappers are allowed to be released or not.
 *
 *   resource_ptr
 *      pointer to the SystemResource that manages the EbObjectWrapper.
 *      The emptyFifo's lockout_mutex is used to write protect the
 *      modification of the EbObjectWrapper.
 *
 *   wrapper_ptr
 *      pointer to the EbObjectWrapper to be modified.
 *********************************************************************/
EbErrorType eb_object_release_disable(
    EbObjectWrapper   *wrapper_ptr)
{
    EbErrorType return_error = EB_ErrorNone;

    eb_block_on_mutex(wrapper_ptr->system_resource_ptr->empty_queue->lockout_mutex);

    wrapper_ptr->release_enable = EB_FALSE;

    eb_release_mutex(wrapper_ptr->system_resource_ptr->empty_queue->lockout_mutex);

    return return_error;
}

/*********************************************************************
 * eb_object_inc_live_count
 *   Increments the live_count member of EbObjectWrapper.  Used by
 *   certain objects (e.g. SequenceControlSet) to count the number of active
 *   pointers of a EbObjectWrapper in pipeline at any point in time.
 *
 *   resource_ptr
 *      pointer to the SystemResource that manages the EbObjectWrapper.
 *      The emptyFifo's lockout_mutex is used to write protect the
 *      modification of the EbObjectWrapper.
 *
 *   wrapper_ptr
 *      pointer to the EbObjectWrapper to be modified.
 *********************************************************************/
EbErrorType eb_object_inc_live_count(
    EbObjectWrapper   *wrapper_ptr,
    uint32_t               increment_number)
{
    EbErrorType return_error = EB_ErrorNone;

    eb_block_on_mutex(wrapper_ptr->system_resource_ptr->empty_queue->lockout_mutex);

    wrapper_ptr->live_count += increment_number;

    eb_release_mutex(wrapper_ptr->system_resource_ptr->empty_queue->lockout_mutex);

    return return_error;
}

//ugly hack
typedef struct DctorAble
{
    EbDctor dctor;
} DctorAble;

void eb_object_wrapper_dctor(EbPtr p)
{
    EbObjectWrapper* wrapper = (EbObjectWrapper*)p;
    if (wrapper->object_destroyer) {
        //customized destoryer
        if (wrapper->object_ptr)
            wrapper->object_destroyer(wrapper->object_ptr);
    }
    else {
        //hack....
        DctorAble* obj= (DctorAble*)wrapper->object_ptr;
        EB_DELETE(obj);
    }
}

static EbErrorType eb_object_wrapper_ctor(EbObjectWrapper* wrapper,
    EbSystemResource    *resource,
    EbCreator           object_creator,
    EbPtr               object_init_data_ptr,
    EbDctor             object_destroyer)
{
    EbErrorType ret;

    wrapper->dctor = eb_object_wrapper_dctor;
    ret = object_creator(&wrapper->object_ptr, object_init_data_ptr);
    if (ret != EB_ErrorNone)
        return ret;
    wrapper->release_enable = EB_TRUE;
    wrapper->system_resource_ptr = resource;
    wrapper->object_destroyer = object_destroyer;
    return EB_ErrorNone;
}

static void eb_system_resource_dctor(EbPtr p)
{
    EbSystemResource* obj = (EbSystemResource*)p;
    EB_DELETE(obj->full_queue);
    EB_DELETE(obj->empty_queue);
    EB_DELETE_PTR_ARRAY(obj->wrapper_ptr_pool, obj->object_total_count);
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
 *   object_ctor
 *     Function pointer to the constructor of the object managed by
 *     SystemResource referenced by resource_ptr. No object level
 *     construction is performed if object_ctor is NULL.
 *
 *   object_init_data_ptr

 *     pointer to data block to be used during the construction of
 *     the object. object_init_data_ptr is passed to object_ctor when
 *     object_ctor is called.
 *   object_destroyer
 *     object destroyer, will call dctor if this is null
 *********************************************************************/
EbErrorType eb_system_resource_ctor(
    EbSystemResource *resource_ptr,
    uint32_t               object_total_count,
    uint32_t               producer_process_total_count,
    uint32_t               consumer_process_total_count,
    EbCreator           object_creator,
    EbPtr               object_init_data_ptr,
    EbDctor             object_destroyer)
{
    uint32_t wrapperIndex;
    EbErrorType return_error = EB_ErrorNone;
    resource_ptr->dctor = eb_system_resource_dctor;

    resource_ptr->object_total_count = object_total_count;

    // Allocate array for wrapper pointers
    EB_ALLOC_PTR_ARRAY(resource_ptr->wrapper_ptr_pool, resource_ptr->object_total_count);

    // Initialize each wrapper
    for (wrapperIndex = 0; wrapperIndex < resource_ptr->object_total_count; ++wrapperIndex) {
        EB_NEW(resource_ptr->wrapper_ptr_pool[wrapperIndex], eb_object_wrapper_ctor, resource_ptr,
            object_creator, object_init_data_ptr, object_destroyer);
    }

    // Initialize the Empty Queue
    EB_NEW(
        resource_ptr->empty_queue,
        EbMuxingQueueCtor,
        resource_ptr->object_total_count,
        producer_process_total_count);
    // Fill the Empty Fifo with every ObjectWrapper
    for (wrapperIndex = 0; wrapperIndex < resource_ptr->object_total_count; ++wrapperIndex) {
        EbMuxingQueueObjectPushBack(
            resource_ptr->empty_queue,
            resource_ptr->wrapper_ptr_pool[wrapperIndex]);
    }

    // Initialize the Full Queue
    if (consumer_process_total_count) {
        EB_NEW(
            resource_ptr->full_queue,
            EbMuxingQueueCtor,
            resource_ptr->object_total_count,
            consumer_process_total_count);
    }
    else {
        resource_ptr->full_queue = (EbMuxingQueue *)EB_NULL;
    }

    return return_error;
}

EbFifo* eb_system_resource_get_producer_fifo(const EbSystemResource *resource_ptr, uint32_t index)
{
    return eb_muxing_queue_get_fifo(resource_ptr->empty_queue, index);
}

EbFifo* eb_system_resource_get_consumer_fifo(const EbSystemResource *resource_ptr, uint32_t index)
{
    return eb_muxing_queue_get_fifo(resource_ptr->full_queue, index);
}


/*********************************************************************
 * EbSystemResourceReleaseProcess
 *********************************************************************/
static EbErrorType EbReleaseProcess(
    EbFifo   *processFifoPtr)
{
    EbErrorType return_error = EB_ErrorNone;

    eb_block_on_mutex(processFifoPtr->queue_ptr->lockout_mutex);

    EbCircularBufferPushFront(
        processFifoPtr->queue_ptr->process_queue,
        processFifoPtr);

    EbMuxingQueueAssignation(processFifoPtr->queue_ptr);

    eb_release_mutex(processFifoPtr->queue_ptr->lockout_mutex);

    return return_error;
}

static EbErrorType EbReleaseProcessNonBlocking (
    EbFifo   *processFifoPtr)
{
    EbErrorType return_error = EB_ErrorNone;

    eb_block_on_mutex(processFifoPtr->queue_ptr->lockout_mutex);

    EbCircularBufferPushFront(
        processFifoPtr->queue_ptr->process_queue,
        processFifoPtr);

    EbMuxingQueueAssignationNonBlocking(processFifoPtr->queue_ptr);

    eb_release_mutex(processFifoPtr->queue_ptr->lockout_mutex);

    return return_error;
}

/*********************************************************************
 * EbSystemResourcePostObject
 *   Queues a full EbObjectWrapper to the SystemResource. This
 *   function posts the SystemResource fullFifo counting_semaphore.
 *   This function is write protected by the SystemResource fullFifo
 *   lockout_mutex.
 *
 *   resource_ptr
 *      pointer to the SystemResource that the EbObjectWrapper is
 *      posted to.
 *
 *   wrapper_ptr
 *      pointer to EbObjectWrapper to be posted.
 *********************************************************************/
EbErrorType eb_post_full_object(
    EbObjectWrapper   *object_ptr)
{
    EbErrorType return_error = EB_ErrorNone;

    eb_block_on_mutex(object_ptr->system_resource_ptr->full_queue->lockout_mutex);

    EbMuxingQueueObjectPushBack(
        object_ptr->system_resource_ptr->full_queue,
        object_ptr);

    eb_release_mutex(object_ptr->system_resource_ptr->full_queue->lockout_mutex);

    return return_error;
}

/*********************************************************************
 * EbSystemResourceReleaseObject
 *   Queues an empty EbObjectWrapper to the SystemResource. This
 *   function posts the SystemResource emptyFifo counting_semaphore.
 *   This function is write protected by the SystemResource emptyFifo
 *   lockout_mutex.
 *
 *   object_ptr
 *      pointer to EbObjectWrapper to be released.
 *********************************************************************/
EbErrorType eb_release_object(
    EbObjectWrapper   *object_ptr)
{
    EbErrorType return_error = EB_ErrorNone;

    eb_block_on_mutex(object_ptr->system_resource_ptr->empty_queue->lockout_mutex);

    // Decrement live_count
    object_ptr->live_count = (object_ptr->live_count == 0) ? object_ptr->live_count : object_ptr->live_count - 1;

    if ((object_ptr->release_enable == EB_TRUE) && (object_ptr->live_count == 0)) {
        // Set live_count to EB_ObjectWrapperReleasedValue
        object_ptr->live_count = EB_ObjectWrapperReleasedValue;

        EbMuxingQueueObjectPushFront(
            object_ptr->system_resource_ptr->empty_queue,
            object_ptr);
    }

    eb_release_mutex(object_ptr->system_resource_ptr->empty_queue->lockout_mutex);

    return return_error;
}

/*********************************************************************
 * EbSystemResourceGetEmptyObject
 *   Dequeues an empty EbObjectWrapper from the SystemResource.  This
 *   function blocks on the SystemResource emptyFifo counting_semaphore.
 *   This function is write protected by the SystemResource emptyFifo
 *   lockout_mutex.
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
    EbFifo   *empty_fifo_ptr,
    EbObjectWrapper **wrapper_dbl_ptr)
{
    EbErrorType return_error = EB_ErrorNone;

    // Queue the Fifo requesting the empty fifo
    EbReleaseProcess(empty_fifo_ptr);

    // Block on the counting Semaphore until an empty buffer is available
    eb_block_on_semaphore(empty_fifo_ptr->counting_semaphore);

    // Acquire lockout Mutex
    eb_block_on_mutex(empty_fifo_ptr->lockout_mutex);

    // Get the empty object
    EbFifoPopFront(
        empty_fifo_ptr,
        wrapper_dbl_ptr);

    // Reset the wrapper's live_count
    (*wrapper_dbl_ptr)->live_count = 0;

    // Object release enable
    (*wrapper_dbl_ptr)->release_enable = EB_TRUE;

    // Release Mutex
    eb_release_mutex(empty_fifo_ptr->lockout_mutex);

    return return_error;
}

/*********************************************************************
 * EbSystemResourceGetFullObject
 *   Dequeues an full EbObjectWrapper from the SystemResource. This
 *   function blocks on the SystemResource fullFifo counting_semaphore.
 *   This function is write protected by the SystemResource fullFifo
 *   lockout_mutex.
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
    EbFifo   *full_fifo_ptr,
    EbObjectWrapper **wrapper_dbl_ptr)
{
    EbErrorType return_error = EB_ErrorNone;

    // Queue the Fifo requesting the full fifo
    EbReleaseProcess(full_fifo_ptr);

    // Block on the counting Semaphore until an empty buffer is available
    eb_block_on_semaphore(full_fifo_ptr->counting_semaphore);

    // Acquire lockout Mutex
    eb_block_on_mutex(full_fifo_ptr->lockout_mutex);

    EbFifoPopFront(
        full_fifo_ptr,
        wrapper_dbl_ptr);

    // Release Mutex
    eb_release_mutex(full_fifo_ptr->lockout_mutex);

    return return_error;
}

EbErrorType eb_get_non_blocking_full_object(
    EbFifo   *full_fifo_ptr,
    EbObjectWrapper **wrapper_dbl_ptr)
{
    EbErrorType return_error = EB_ErrorNone;

    // Queue the Fifo requesting the full fifo
    EbReleaseProcessNonBlocking(full_fifo_ptr);

    if (full_fifo_ptr->first_ptr == NULL)
        *wrapper_dbl_ptr = NULL;
    else
        EbFifoPopFront(
            full_fifo_ptr,
            wrapper_dbl_ptr);

    return return_error;
}


/**************************************
* EbFifoPopFront
**************************************/
static EbBool EbFifoPeakFront(
    EbFifo            *fifoPtr)
{
    // Set wrapper_ptr to head of BufferPool
    if (fifoPtr->first_ptr == (EbObjectWrapper*)EB_NULL)
        return EB_TRUE;
    else
        return EB_FALSE;
}

/* NonBlocking Get Object Modified for Faster Row Level Jobs of Decoder */
EbErrorType eb_dec_get_full_object_non_blocking(
    EbFifo   *full_fifo_ptr,
    EbObjectWrapper **wrapper_dbl_ptr)
{
    EbErrorType return_error = EB_ErrorNone;
    EbBool      fifoEmpty;
    // Queue the Fifo requesting the full fifo
    EbReleaseProcessNonBlocking(full_fifo_ptr);

    // Acquire lockout Mutex
    eb_block_on_mutex(full_fifo_ptr->lockout_mutex);

    fifoEmpty = EbFifoPeakFront(
        full_fifo_ptr);

    // Release Mutex
    eb_release_mutex(full_fifo_ptr->lockout_mutex);

    if (fifoEmpty == EB_FALSE)
        eb_get_non_blocking_full_object(
            full_fifo_ptr,
            wrapper_dbl_ptr);
    else
        *wrapper_dbl_ptr = (EbObjectWrapper*)EB_NULL;

    return return_error;
}

EbErrorType eb_get_full_object_non_blocking(
    EbFifo   *full_fifo_ptr,
    EbObjectWrapper **wrapper_dbl_ptr)
{
    EbErrorType return_error = EB_ErrorNone;
    EbBool      fifoEmpty;
    // Queue the Fifo requesting the full fifo
    EbReleaseProcess(full_fifo_ptr);

    // Acquire lockout Mutex
    eb_block_on_mutex(full_fifo_ptr->lockout_mutex);

    fifoEmpty = EbFifoPeakFront(
        full_fifo_ptr);

    // Release Mutex
    eb_release_mutex(full_fifo_ptr->lockout_mutex);

    if (fifoEmpty == EB_FALSE)
        eb_get_full_object(
            full_fifo_ptr,
            wrapper_dbl_ptr);
    else
        *wrapper_dbl_ptr = (EbObjectWrapper*)EB_NULL;

    return return_error;
}
// clang-format on
