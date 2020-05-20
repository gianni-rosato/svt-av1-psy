/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

// Summary:
// EbThreads contains wrappers functions that hide
// platform specific objects such as threads, semaphores,
// and mutexs.  The goal is to eliminiate platform #define
// in the code.

#if defined(__has_feature)
#  if __has_feature(thread_sanitizer)
#    define EB_THREAD_SANITIZER_ENABLED
#  endif
#endif

/****************************************
 * Universal Includes
 ****************************************/
#include <stdlib.h>
#include "EbThreads.h"
#include "EbLog.h"
/****************************************
  * Win32 Includes
  ****************************************/
#ifdef _WIN32
#include <windows.h>
#else
#include <stdio.h>
#include <errno.h>
#include <fcntl.h>
#include <pthread.h>
#include <semaphore.h>
#include <unistd.h>
#endif // _WIN32
#ifdef __APPLE__
#include <dispatch/dispatch.h>
#endif
#if PRINTF_TIME
#include <time.h>
#ifdef _WIN32
void printfTime(const char *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    SVT_LOG("  [%i ms]\t", ((int32_t)clock()));
    vprintf(fmt, args);
    va_end(args);
}
#endif
#endif

/****************************************
 * eb_create_thread
 ****************************************/
EbHandle eb_create_thread(void *thread_function(void *), void *thread_context) {
    EbHandle thread_handle = NULL;

#ifdef _WIN32

    thread_handle = (EbHandle)CreateThread(
        NULL, // default security attributes
        0, // default stack size
        (LPTHREAD_START_ROUTINE)thread_function, // function to be tied to the new thread
        thread_context, // context to be tied to the new thread
        0, // thread active when created
        NULL); // new thread ID

#else

    int ret;
    pthread_t *th;
    pthread_attr_t attr;

    th = malloc(sizeof(pthread_t));
    if (th == NULL)
        return NULL;

#ifndef EB_THREAD_SANITIZER_ENABLED
    pthread_attr_init(&attr);
    pthread_attr_setschedpolicy(&attr, SCHED_FIFO);
    pthread_attr_setinheritsched(&attr, PTHREAD_EXPLICIT_SCHED);

    struct sched_param param = {.sched_priority = 99};
    pthread_attr_setschedparam(&attr, &param);

    ret = pthread_create(th, &attr, thread_function, thread_context);
    pthread_attr_destroy(&attr);

    if (ret == EPERM) {
        // When creating the thread failed because setting scheduling
        // parameters failed, retry creating the thread without them.
        ret = pthread_create(th, NULL, thread_function, thread_context);
    }
#else
    // When running with thread sanitizer, we are not running as root
    // so the above priority change will always fail, which will cause
    // issues with the thread sanitizer.
    // See https://github.com/google/sanitizers/issues/1088
    // Therefore we never try this, with the thread sanitizer
    // and just create a normal thread here.
    ret = pthread_create(th, NULL, thread_function, thread_context);
#endif

    if (ret != 0) {
        free(th);
        return NULL;
    }
    thread_handle = th;
#endif // _WIN32

    return thread_handle;
}

///****************************************
// * eb_start_thread
// ****************************************/
//EbErrorType eb_start_thread(
//    EbHandle thread_handle)
//{
//    EbErrorType error_return = EB_ErrorNone;
//
//    /* Note JMJ 9/6/2011
//        The thread Pause/Resume functionality is being removed.  The main reason is that
//        POSIX Threads (aka pthreads) does not support this functionality.  The destructor
//        and deinit code is safe as along as when EbDestropyThread is called on a thread,
//        the thread is immediately destroyed and its stack cleared.
//
//        The Encoder Start/Stop functionality, which previously used the thread Pause/Resume
//        functions could be implemented with mutex checks either at the head of the pipeline,
//        or throughout the code if a more responsive Pause is needed.
//    */
//
//#ifdef _WIN32
//    //error_return = ResumeThread((HANDLE) thread_handle) ? EB_ErrorThreadUnresponsive : EB_ErrorNone;
//#else
//#endif // _WIN32
//
//    error_return = (thread_handle) ? EB_ErrorNone : EB_ErrorNullThread;
//
//    return error_return;
//}
//
///****************************************
// * eb_stop_thread
// ****************************************/
//EbErrorType eb_stop_thread(
//    EbHandle thread_handle)
//{
//    EbErrorType error_return = EB_ErrorNone;
//
//#ifdef _WIN32
//    //error_return = SuspendThread((HANDLE) thread_handle) ? EB_ErrorThreadUnresponsive : EB_ErrorNone;
//#else
//#endif // _WIN32
//
//    error_return = (thread_handle) ? EB_ErrorNone : EB_ErrorNullThread;
//
//    return error_return;
//}
//
/****************************************
 * eb_destroy_thread
 ****************************************/
EbErrorType eb_destroy_thread(EbHandle thread_handle) {
    EbErrorType error_return = EB_ErrorNone;

#ifdef _WIN32
    WaitForSingleObject(thread_handle, INFINITE);
    error_return = CloseHandle(thread_handle) ? EB_ErrorNone : EB_ErrorDestroyThreadFailed;
#else
    pthread_join(*((pthread_t *)thread_handle), NULL);
    free(thread_handle);
#endif // _WIN32

    return error_return;
}

/***************************************
 * eb_create_semaphore
 ***************************************/
EbHandle eb_create_semaphore(uint32_t initial_count, uint32_t max_count)
{
    EbHandle semaphore_handle;

#if defined(_WIN32)
    semaphore_handle = (EbHandle)CreateSemaphore(NULL, // default security attributes
                                                 initial_count, // initial semaphore count
                                                 max_count, // maximum semaphore count
                                                 NULL); // semaphore is not named
#elif defined(__APPLE__)
    UNUSED(max_count);
    semaphore_handle = (EbHandle)dispatch_semaphore_create(initial_count);
#else
    UNUSED(max_count);

    semaphore_handle = (sem_t *)malloc(sizeof(sem_t));
    sem_init((sem_t *)semaphore_handle, // semaphore handle
             0, // shared semaphore (not local)
             initial_count); // initial count
#endif

    return semaphore_handle;
}

/***************************************
 * eb_post_semaphore
 ***************************************/
EbErrorType eb_post_semaphore(EbHandle semaphore_handle)
{
    EbErrorType return_error;

#ifdef _WIN32
    return_error = !ReleaseSemaphore(semaphore_handle, // semaphore handle
                                     1, // amount to increment the semaphore
                                     NULL) // pointer to previous count (optional)
                       ? EB_ErrorSemaphoreUnresponsive
                       : EB_ErrorNone;
#elif defined(__APPLE__)
    dispatch_semaphore_signal((dispatch_semaphore_t)semaphore_handle);
    return_error = EB_ErrorNone;
#else
    return_error =
        sem_post((sem_t *)semaphore_handle) ? EB_ErrorSemaphoreUnresponsive : EB_ErrorNone;
#endif

    return return_error;
}

/***************************************
 * eb_block_on_semaphore
 ***************************************/
EbErrorType eb_block_on_semaphore(EbHandle semaphore_handle)
{
    EbErrorType return_error;

#ifdef _WIN32
    return_error = WaitForSingleObject((HANDLE)semaphore_handle, INFINITE)
                       ? EB_ErrorSemaphoreUnresponsive
                       : EB_ErrorNone;
#elif defined(__APPLE__)
    return_error =
        dispatch_semaphore_wait((dispatch_semaphore_t)semaphore_handle, DISPATCH_TIME_FOREVER)
            ? EB_ErrorSemaphoreUnresponsive
            : EB_ErrorNone;
#else
    return_error =
        sem_wait((sem_t *)semaphore_handle) ? EB_ErrorSemaphoreUnresponsive : EB_ErrorNone;
#endif

    return return_error;
}

/***************************************
 * eb_destroy_semaphore
 ***************************************/
EbErrorType eb_destroy_semaphore(EbHandle semaphore_handle)
{
    EbErrorType return_error;

#ifdef _WIN32
    return_error =
        !CloseHandle((HANDLE)semaphore_handle) ? EB_ErrorDestroySemaphoreFailed : EB_ErrorNone;
#elif defined(__APPLE__)
    dispatch_release((dispatch_semaphore_t)semaphore_handle);
    return_error = EB_ErrorNone;
#else
    return_error =
        sem_destroy((sem_t *)semaphore_handle) ? EB_ErrorDestroySemaphoreFailed : EB_ErrorNone;
    free(semaphore_handle);
#endif

    return return_error;
}
/***************************************
 * eb_create_mutex
 ***************************************/
EbHandle eb_create_mutex(void)
{
    EbHandle mutex_handle;

#ifdef _WIN32
    mutex_handle = (EbHandle)CreateMutex(NULL, // default security attributes
                                         FALSE, // FALSE := not initially owned
                                         NULL); // mutex is not named

#else

    mutex_handle = (EbHandle)malloc(sizeof(pthread_mutex_t));

    if (mutex_handle != NULL) {
        pthread_mutex_init((pthread_mutex_t *)mutex_handle,
                           NULL); // default attributes
    }
#endif

    return mutex_handle;
}

/***************************************
 * EbPostMutex
 ***************************************/
EbErrorType eb_release_mutex(EbHandle mutex_handle)
{
    EbErrorType return_error;

#ifdef _WIN32
    return_error = !ReleaseMutex((HANDLE)mutex_handle) ? EB_ErrorCreateMutexFailed : EB_ErrorNone;
#else
    return_error = pthread_mutex_unlock((pthread_mutex_t *)mutex_handle) ? EB_ErrorCreateMutexFailed
                                                                         : EB_ErrorNone;
#endif

    return return_error;
}

/***************************************
 * eb_block_on_mutex
 ***************************************/
EbErrorType eb_block_on_mutex(EbHandle mutex_handle)
{
    EbErrorType return_error;

#ifdef _WIN32
    return_error = WaitForSingleObject((HANDLE)mutex_handle, INFINITE) ? EB_ErrorMutexUnresponsive
                                                                       : EB_ErrorNone;
#else
    return_error = pthread_mutex_lock((pthread_mutex_t *)mutex_handle) ? EB_ErrorMutexUnresponsive
                                                                       : EB_ErrorNone;
#endif

    return return_error;
}

/***************************************
 * eb_destroy_mutex
 ***************************************/
EbErrorType eb_destroy_mutex(EbHandle mutex_handle)
{
    EbErrorType return_error;

#ifdef _WIN32
    return_error = CloseHandle((HANDLE)mutex_handle) ? EB_ErrorDestroyMutexFailed : EB_ErrorNone;
#else
    return_error = pthread_mutex_destroy((pthread_mutex_t *)mutex_handle)
                       ? EB_ErrorDestroyMutexFailed
                       : EB_ErrorNone;
    free(mutex_handle);
#endif

    return return_error;
}
