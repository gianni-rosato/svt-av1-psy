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
    uint32_t              initialCount,
    uint32_t              maxCount,
    EbObjectWrapper_t  *firstWrapperPtr,
    EbObjectWrapper_t  *lastWrapperPtr,
    EbMuxingQueue_t    *queuePtr)
{
    // Create Counting Semaphore
    EB_CREATESEMAPHORE(EbHandle, fifoPtr->countingSemaphore, sizeof(EbHandle), EB_SEMAPHORE, initialCount, maxCount);

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
    EbPtr               *objectPtr)
{
    EbErrorType return_error = EB_ErrorNone;

    // Copy the head of the buffer into the objectPtr
    *objectPtr = bufferPtr->arrayPtr[bufferPtr->headIndex];
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
    EbPtr                objectPtr)
{
    EbErrorType return_error = EB_ErrorNone;

    // Copy the pointer into the array
    bufferPtr->arrayPtr[bufferPtr->tailIndex] = objectPtr;

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
    EbPtr                objectPtr)
{
    EbErrorType return_error = EB_ErrorNone;

    // Decrement the headIndex
    bufferPtr->headIndex = (bufferPtr->headIndex == 0) ? bufferPtr->buffer_total_count - 1 : bufferPtr->headIndex - 1;

    // Copy the pointer into the array
    bufferPtr->arrayPtr[bufferPtr->headIndex] = objectPtr;

    // Increment the Current Count
    ++bufferPtr->currentCount;

    return return_error;
}

/**************************************
 * EbMuxingQueueCtor
 **************************************/
static EbErrorType EbMuxingQueueCtor(
    EbMuxingQueue_t   **queueDblPtr,
    uint32_t              objectTotalCount,
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
        objectTotalCount);
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
            objectTotalCount,
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
        EbBlockOnMutex(processFifoPtr->lockoutMutex);

        // Put the object on the fifo
        EbFifoPushBack(
            processFifoPtr,
            wrapper_ptr);

        // Release the Process Fifo's Mutex
        EbReleaseMutex(processFifoPtr->lockoutMutex);

        // Post the semaphore
        EbPostSemaphore(processFifoPtr->countingSemaphore);
    }

    return return_error;
}

/**************************************
 * EbMuxingQueueObjectPushBack
 **************************************/
static EbErrorType EbMuxingQueueObjectPushBack(
    EbMuxingQueue_t    *queuePtr,
    EbObjectWrapper_t  *objectPtr)
{
    EbErrorType return_error = EB_ErrorNone;

    EbCircularBufferPushBack(
        queuePtr->objectQueue,
        objectPtr);

    EbMuxingQueueAssignation(queuePtr);

    return return_error;
}

/**************************************
* EbMuxingQueueObjectPushFront
**************************************/
static EbErrorType EbMuxingQueueObjectPushFront(
    EbMuxingQueue_t    *queuePtr,
    EbObjectWrapper_t  *objectPtr)
{
    EbErrorType return_error = EB_ErrorNone;

    EbCircularBufferPushFront(
        queuePtr->objectQueue,
        objectPtr);

    EbMuxingQueueAssignation(queuePtr);

    return return_error;
}

/*********************************************************************
 * EbObjectReleaseEnable
 *   Enables the releaseEnable member of EbObjectWrapper.  Used by
 *   certain objects (e.g. SequenceControlSet) to control whether
 *   EbObjectWrappers are allowed to be released or not.
 *
 *   resourcePtr
 *      Pointer to the SystemResource that manages the EbObjectWrapper.
 *      The emptyFifo's lockoutMutex is used to write protect the
 *      modification of the EbObjectWrapper.
 *
 *   wrapper_ptr
 *      Pointer to the EbObjectWrapper to be modified.
 *********************************************************************/
EbErrorType EbObjectReleaseEnable(
    EbObjectWrapper_t   *wrapper_ptr)
{
    EbErrorType return_error = EB_ErrorNone;

    EbBlockOnMutex(wrapper_ptr->systemResourcePtr->emptyQueue->lockoutMutex);

    wrapper_ptr->releaseEnable = EB_TRUE;

    EbReleaseMutex(wrapper_ptr->systemResourcePtr->emptyQueue->lockoutMutex);

    return return_error;
}

/*********************************************************************
 * EbObjectReleaseDisable
 *   Disables the releaseEnable member of EbObjectWrapper.  Used by
 *   certain objects (e.g. SequenceControlSet) to control whether
 *   EbObjectWrappers are allowed to be released or not.
 *
 *   resourcePtr
 *      Pointer to the SystemResource that manages the EbObjectWrapper.
 *      The emptyFifo's lockoutMutex is used to write protect the
 *      modification of the EbObjectWrapper.
 *
 *   wrapper_ptr
 *      Pointer to the EbObjectWrapper to be modified.
 *********************************************************************/
EbErrorType EbObjectReleaseDisable(
    EbObjectWrapper_t   *wrapper_ptr)
{
    EbErrorType return_error = EB_ErrorNone;

    EbBlockOnMutex(wrapper_ptr->systemResourcePtr->emptyQueue->lockoutMutex);

    wrapper_ptr->releaseEnable = EB_FALSE;

    EbReleaseMutex(wrapper_ptr->systemResourcePtr->emptyQueue->lockoutMutex);

    return return_error;
}

/*********************************************************************
 * EbObjectIncLiveCount
 *   Increments the liveCount member of EbObjectWrapper.  Used by
 *   certain objects (e.g. SequenceControlSet) to count the number of active
 *   pointers of a EbObjectWrapper in pipeline at any point in time.
 *
 *   resourcePtr
 *      Pointer to the SystemResource that manages the EbObjectWrapper.
 *      The emptyFifo's lockoutMutex is used to write protect the
 *      modification of the EbObjectWrapper.
 *
 *   wrapper_ptr
 *      Pointer to the EbObjectWrapper to be modified.
 *********************************************************************/
EbErrorType EbObjectIncLiveCount(
    EbObjectWrapper_t   *wrapper_ptr,
    uint32_t               incrementNumber)
{
    EbErrorType return_error = EB_ErrorNone;

    EbBlockOnMutex(wrapper_ptr->systemResourcePtr->emptyQueue->lockoutMutex);

    wrapper_ptr->liveCount += incrementNumber;

    EbReleaseMutex(wrapper_ptr->systemResourcePtr->emptyQueue->lockoutMutex);

    return return_error;
}

/*********************************************************************
 * EbSystemResourceCtor
 *   Constructor for EbSystemResource.  Fully constructs all members
 *   of EbSystemResource including the object with the passed
 *   ObjectCtor function.
 *
 *   resourcePtr
 *     Pointer that will contain the SystemResource to be constructed.
 *
 *   objectTotalCount
 *     Number of objects to be managed by the SystemResource.
 *
 *   fullFifoEnabled
 *     Bool that describes if the SystemResource is to have an output
 *     fifo.  An outputFifo is not used by certain objects (e.g.
 *     SequenceControlSet).
 *
 *   ObjectCtor
 *     Function pointer to the constructor of the object managed by
 *     SystemResource referenced by resourcePtr. No object level
 *     construction is performed if ObjectCtor is NULL.
 *
 *   object_init_data_ptr

 *     Pointer to data block to be used during the construction of
 *     the object. object_init_data_ptr is passed to ObjectCtor when
 *     ObjectCtor is called.
 *********************************************************************/
EbErrorType EbSystemResourceCtor(
    EbSystemResource_t **resourceDblPtr,
    uint32_t               objectTotalCount,
    uint32_t               producerProcessTotalCount,
    uint32_t               consumerProcessTotalCount,
    EbFifo_t          ***producerFifoPtrArrayPtr,
    EbFifo_t          ***consumerFifoPtrArrayPtr,
    EbBool              fullFifoEnabled,
    EB_CTOR              ObjectCtor,
    EbPtr               object_init_data_ptr)
{
    uint32_t wrapperIndex;
    EbErrorType return_error = EB_ErrorNone;
    // Allocate the System Resource
    EbSystemResource_t *resourcePtr;

    EB_MALLOC(EbSystemResource_t*, resourcePtr, sizeof(EbSystemResource_t), EB_N_PTR);
    *resourceDblPtr = resourcePtr;

    resourcePtr->objectTotalCount = objectTotalCount;

    // Allocate array for wrapper pointers
    EB_MALLOC(EbObjectWrapper_t**, resourcePtr->wrapperPtrPool, sizeof(EbObjectWrapper_t*) * resourcePtr->objectTotalCount, EB_N_PTR);

    // Initialize each wrapper
    for (wrapperIndex = 0; wrapperIndex < resourcePtr->objectTotalCount; ++wrapperIndex) {
        EB_MALLOC(EbObjectWrapper_t*, resourcePtr->wrapperPtrPool[wrapperIndex], sizeof(EbObjectWrapper_t), EB_N_PTR);
        resourcePtr->wrapperPtrPool[wrapperIndex]->liveCount = 0;
        resourcePtr->wrapperPtrPool[wrapperIndex]->releaseEnable = EB_TRUE;
        resourcePtr->wrapperPtrPool[wrapperIndex]->systemResourcePtr = resourcePtr;

        // Call the Constructor for each element
        if (ObjectCtor) {
            return_error = ObjectCtor(
                &resourcePtr->wrapperPtrPool[wrapperIndex]->objectPtr,
                object_init_data_ptr);
            if (return_error == EB_ErrorInsufficientResources) {
                return EB_ErrorInsufficientResources;
            }
        }
    }

    // Initialize the Empty Queue
    return_error = EbMuxingQueueCtor(
        &resourcePtr->emptyQueue,
        resourcePtr->objectTotalCount,
        producerProcessTotalCount,
        producerFifoPtrArrayPtr);
    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }
    // Fill the Empty Fifo with every ObjectWrapper
    for (wrapperIndex = 0; wrapperIndex < resourcePtr->objectTotalCount; ++wrapperIndex) {
        EbMuxingQueueObjectPushBack(
            resourcePtr->emptyQueue,
            resourcePtr->wrapperPtrPool[wrapperIndex]);
    }

    // Initialize the Full Queue
    if (fullFifoEnabled == EB_TRUE) {
        return_error = EbMuxingQueueCtor(
            &resourcePtr->fullQueue,
            resourcePtr->objectTotalCount,
            consumerProcessTotalCount,
            consumerFifoPtrArrayPtr);
        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }
    else {
        resourcePtr->fullQueue = (EbMuxingQueue_t *)EB_NULL;
        consumerFifoPtrArrayPtr = (EbFifo_t ***)EB_NULL;
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

    EbBlockOnMutex(processFifoPtr->queuePtr->lockoutMutex);

    EbCircularBufferPushFront(
        processFifoPtr->queuePtr->processQueue,
        processFifoPtr);

    EbMuxingQueueAssignation(processFifoPtr->queuePtr);

    EbReleaseMutex(processFifoPtr->queuePtr->lockoutMutex);

    return return_error;
}

/*********************************************************************
 * EbSystemResourcePostObject
 *   Queues a full EbObjectWrapper to the SystemResource. This
 *   function posts the SystemResource fullFifo countingSemaphore.
 *   This function is write protected by the SystemResource fullFifo
 *   lockoutMutex.
 *
 *   resourcePtr
 *      Pointer to the SystemResource that the EbObjectWrapper is
 *      posted to.
 *
 *   wrapper_ptr
 *      Pointer to EbObjectWrapper to be posted.
 *********************************************************************/
EbErrorType EbPostFullObject(
    EbObjectWrapper_t   *objectPtr)
{
    EbErrorType return_error = EB_ErrorNone;

    EbBlockOnMutex(objectPtr->systemResourcePtr->fullQueue->lockoutMutex);

    EbMuxingQueueObjectPushBack(
        objectPtr->systemResourcePtr->fullQueue,
        objectPtr);

    EbReleaseMutex(objectPtr->systemResourcePtr->fullQueue->lockoutMutex);

    return return_error;
}

/*********************************************************************
 * EbSystemResourceReleaseObject
 *   Queues an empty EbObjectWrapper to the SystemResource. This
 *   function posts the SystemResource emptyFifo countingSemaphore.
 *   This function is write protected by the SystemResource emptyFifo
 *   lockoutMutex.
 *
 *   objectPtr
 *      Pointer to EbObjectWrapper to be released.
 *********************************************************************/
EbErrorType EbReleaseObject(
    EbObjectWrapper_t   *objectPtr)
{
    EbErrorType return_error = EB_ErrorNone;

    EbBlockOnMutex(objectPtr->systemResourcePtr->emptyQueue->lockoutMutex);

    // Decrement liveCount
    objectPtr->liveCount = (objectPtr->liveCount == 0) ? objectPtr->liveCount : objectPtr->liveCount - 1;

    if ((objectPtr->releaseEnable == EB_TRUE) && (objectPtr->liveCount == 0)) {

        // Set liveCount to EB_ObjectWrapperReleasedValue
        objectPtr->liveCount = EB_ObjectWrapperReleasedValue;

        EbMuxingQueueObjectPushFront(
            objectPtr->systemResourcePtr->emptyQueue,
            objectPtr);

    }

    EbReleaseMutex(objectPtr->systemResourcePtr->emptyQueue->lockoutMutex);

    return return_error;
}

/*********************************************************************
 * EbSystemResourceGetEmptyObject
 *   Dequeues an empty EbObjectWrapper from the SystemResource.  This
 *   function blocks on the SystemResource emptyFifo countingSemaphore.
 *   This function is write protected by the SystemResource emptyFifo
 *   lockoutMutex.
 *
 *   resourcePtr
 *      Pointer to the SystemResource that provides the empty
 *      EbObjectWrapper.
 *
 *   wrapperDblPtr
 *      Double pointer used to pass the pointer to the empty
 *      EbObjectWrapper pointer.
 *********************************************************************/
EbErrorType EbGetEmptyObject(
    EbFifo_t   *emptyFifoPtr,
    EbObjectWrapper_t **wrapperDblPtr)
{
    EbErrorType return_error = EB_ErrorNone;

    // Queue the Fifo requesting the empty fifo
    EbReleaseProcess(emptyFifoPtr);

    // Block on the counting Semaphore until an empty buffer is available
    EbBlockOnSemaphore(emptyFifoPtr->countingSemaphore);

    // Acquire lockout Mutex
    EbBlockOnMutex(emptyFifoPtr->lockoutMutex);

    // Get the empty object
    EbFifoPopFront(
        emptyFifoPtr,
        wrapperDblPtr);

    // Reset the wrapper's liveCount
    (*wrapperDblPtr)->liveCount = 0;

    // Object release enable
    (*wrapperDblPtr)->releaseEnable = EB_TRUE;

    // Release Mutex
    EbReleaseMutex(emptyFifoPtr->lockoutMutex);

    return return_error;
}

/*********************************************************************
 * EbSystemResourceGetFullObject
 *   Dequeues an full EbObjectWrapper from the SystemResource. This
 *   function blocks on the SystemResource fullFifo countingSemaphore.
 *   This function is write protected by the SystemResource fullFifo
 *   lockoutMutex.
 *
 *   resourcePtr
 *      Pointer to the SystemResource that provides the full
 *      EbObjectWrapper.
 *
 *   wrapperDblPtr
 *      Double pointer used to pass the pointer to the full
 *      EbObjectWrapper pointer.
 *********************************************************************/
EbErrorType EbGetFullObject(
    EbFifo_t   *fullFifoPtr,
    EbObjectWrapper_t **wrapperDblPtr)
{
    EbErrorType return_error = EB_ErrorNone;

    // Queue the Fifo requesting the full fifo
    EbReleaseProcess(fullFifoPtr);

    // Block on the counting Semaphore until an empty buffer is available
    EbBlockOnSemaphore(fullFifoPtr->countingSemaphore);

    // Acquire lockout Mutex
    EbBlockOnMutex(fullFifoPtr->lockoutMutex);

    EbFifoPopFront(
        fullFifoPtr,
        wrapperDblPtr);

    // Release Mutex
    EbReleaseMutex(fullFifoPtr->lockoutMutex);

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


EbErrorType EbGetFullObjectNonBlocking(
    EbFifo_t   *fullFifoPtr,
    EbObjectWrapper_t **wrapperDblPtr)
{
    EbErrorType return_error = EB_ErrorNone;
    EbBool      fifoEmpty;
    // Queue the Fifo requesting the full fifo
    EbReleaseProcess(fullFifoPtr);

    // Acquire lockout Mutex
    EbBlockOnMutex(fullFifoPtr->lockoutMutex);

    fifoEmpty = EbFifoPeakFront(
        fullFifoPtr);

    // Release Mutex
    EbReleaseMutex(fullFifoPtr->lockoutMutex);

    if (fifoEmpty == EB_FALSE)
        EbGetFullObject(
            fullFifoPtr,
            wrapperDblPtr);
    else
        *wrapperDblPtr = (EbObjectWrapper_t*)EB_NULL;

    return return_error;
}