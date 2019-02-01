/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbSystemResource_h
#define EbSystemResource_h

#include "EbDefinitions.h"
#include "EbThreads.h"
#ifdef __cplusplus
extern "C" {
#endif
    /*********************************
     * Defines
     *********************************/
#define EB_ObjectWrapperReleasedValue   ~0u

     /*********************************************************************
      * Object Wrapper
      *   Provides state information for each type of object in the
      *   encoder system (i.e. SequenceControlSet, PictureBufferDesc,
      *   ProcessResults, and GracefulDegradation)
      *********************************************************************/
    typedef struct EbObjectWrapper_s {
        // objectPtr - pointer to the object being managed.
        void                     *objectPtr;

        // liveCount - a count of the number of pictures actively being
        //   encoded in the pipeline at any given time.  Modification
        //   of this value by any process must be protected by a mutex.
        uint32_t                    liveCount;

        // releaseEnable - a flag that enables the release of
        //   EbObjectWrapper for reuse in the encoding of subsequent
        //   pictures in the encoder pipeline.
        EbBool                   releaseEnable;

        // systemResourcePtr - a pointer to the SystemResourceManager
        //   that the object belongs to.
        struct EbSystemResource_s *systemResourcePtr;

        // nextPtr - a pointer to a different EbObjectWrapper.  Used
        //   only in the implemenation of a single-linked Fifo.
        struct EbObjectWrapper_s *nextPtr;

    } EbObjectWrapper_t;

    /*********************************************************************
     * Fifo
     *   Defines a static (i.e. no dynamic memory allocation) single
     *   linked-list, constant time fifo implmentation. The fifo uses
     *   the EbObjectWrapper member nextPtr to create the linked-list.
     *   The Fifo also contains a countingSemaphore for OS thread-blocking
     *   and dynamic EbObjectWrapper counting.
     *********************************************************************/
    typedef struct EbFifo_s {
        // countingSemaphore - used for OS thread-blocking & dynamically
        //   counting the number of EbObjectWrappers currently in the
        //   EbFifo.
        EbHandle countingSemaphore;

        // lockoutMutex - used to prevent more than one thread from
        //   modifying EbFifo simultaneously.
        EbHandle lockoutMutex;

        // firstPtr - pointer the the head of the Fifo
        EbObjectWrapper_t *firstPtr;

        // lastPtr - pointer to the tail of the Fifo
        EbObjectWrapper_t *lastPtr;

        // queuePtr - pointer to MuxingQueue that the EbFifo is
        //   associated with.
        struct EbMuxingQueue_s *queuePtr;

    } EbFifo_t;

    /*********************************************************************
     * CircularBuffer
     *********************************************************************/
    typedef struct EbCircularBuffer_s {
        EbPtr *arrayPtr;
        uint32_t  headIndex;
        uint32_t  tailIndex;
        uint32_t  buffer_total_count;
        uint32_t  currentCount;

    } EbCircularBuffer_t;

    /*********************************************************************
     * MuxingQueue
     *********************************************************************/
    typedef struct EbMuxingQueue_s {
        EbHandle           lockoutMutex;
        EbCircularBuffer_t *objectQueue;
        EbCircularBuffer_t *processQueue;
        uint32_t              processTotalCount;
        EbFifo_t          **processFifoPtrArray;

    } EbMuxingQueue_t;

    /*********************************************************************
     * SystemResource
     *   Defines a complete solution for managing objects in the encoder
     *   system (i.e. SequenceControlSet, PictureBufferDesc, ProcessResults, and
     *   GracefulDegradation).  The objectTotalCount and wrapperPtrPool are
     *   only used to construct and destruct the SystemResource.  The
     *   fullFifo provides downstream pipeline data flow control.  The
     *   emptyFifo provides upstream pipeline backpressure flow control.
     *********************************************************************/
    typedef struct EbSystemResource_s {
        // objectTotalCount - A count of the number of objects contained in the
        //   System Resoruce.
        uint32_t              objectTotalCount;

        // wrapperPtrPool - An array of pointers to the EbObjectWrappers used
        //   to construct and destruct the SystemResource.
        EbObjectWrapper_t **wrapperPtrPool;

        // The empty FIFO contains a queue of empty buffers
        //EbFifo_t           *emptyFifo;
        EbMuxingQueue_t     *emptyQueue;

        // The full FIFO contains a queue of completed buffers
        //EbFifo_t           *fullFifo;
        EbMuxingQueue_t     *fullQueue;

    } EbSystemResource_t;

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
    extern EbErrorType EbObjectReleaseEnable(
        EbObjectWrapper_t   *wrapper_ptr);

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
    extern EbErrorType EbObjectReleaseDisable(
        EbObjectWrapper_t   *wrapper_ptr);

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
     *
     *   incrementNumber
     *      The number to increment the live count by.
     *********************************************************************/
    extern EbErrorType EbObjectIncLiveCount(
        EbObjectWrapper_t   *wrapper_ptr,
        uint32_t               incrementNumber);

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
    extern EbErrorType EbSystemResourceCtor(
        EbSystemResource_t **resourceDblPtr,
        uint32_t               objectTotalCount,
        uint32_t               producerProcessTotalCount,
        uint32_t               consumerProcessTotalCount,
        EbFifo_t          ***producerFifoPtrArrayPtr,
        EbFifo_t          ***consumerFifoPtrArrayPtr,
        EbBool              fullFifoEnabled,
        EB_CTOR              ObjectCtor,
        EbPtr               object_init_data_ptr);

    /*********************************************************************
     * EbSystemResourceDtor
     *   Destructor for EbSystemResource.  Fully destructs all members
     *   of EbSystemResource including the object with the passed
     *   ObjectDtor function.
     *
     *   resourcePtr
     *     Pointer to the SystemResource to be destructed.
     *
     *   ObjectDtor
     *     Function pointer to the destructor of the object managed by
     *     SystemResource referenced by resourcePtr. No object level
     *     destruction is performed if ObjectDtor is NULL.
     *********************************************************************/
    extern void EbSystemResourceDtor(
        EbSystemResource_t  *resourcePtr,
        EB_DTOR              ObjectDtor);

    /*********************************************************************
     * EbSystemResourceGetEmptyObject
     *   Dequeues an empty EbObjectWrapper from the SystemResource.  The
     *   new EbObjectWrapper will be populated with the contents of the
     *   wrapperCopyPtr if wrapperCopyPtr is not NULL. This function blocks
     *   on the SystemResource emptyFifo countingSemaphore. This function
     *   is write protected by the SystemResource emptyFifo lockoutMutex.
     *
     *   resourcePtr
     *      Pointer to the SystemResource that provides the empty
     *      EbObjectWrapper.
     *
     *   wrapperDblPtr
     *      Double pointer used to pass the pointer to the empty
     *      EbObjectWrapper pointer.
     *********************************************************************/
    extern EbErrorType EbGetEmptyObject(
        EbFifo_t   *emptyFifoPtr,
        EbObjectWrapper_t **wrapperDblPtr);

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
    extern EbErrorType EbPostFullObject(
        EbObjectWrapper_t   *objectPtr);

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
    extern EbErrorType EbGetFullObject(
        EbFifo_t   *fullFifoPtr,
        EbObjectWrapper_t **wrapperDblPtr);

    extern EbErrorType EbGetFullObjectNonBlocking(
        EbFifo_t   *fullFifoPtr,
        EbObjectWrapper_t **wrapperDblPtr);

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
    extern EbErrorType EbReleaseObject(
        EbObjectWrapper_t   *objectPtr);
#ifdef __cplusplus
}
#endif
#endif //EbSystemResource_h
