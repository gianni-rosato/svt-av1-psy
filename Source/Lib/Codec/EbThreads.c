/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

// Summary:
// EbThreads contains wrappers functions that hide
// platform specific objects such as threads, semaphores,
// and mutexs.  The goal is to eliminiate platform #define
// in the code.

/****************************************
 * Universal Includes
 ****************************************/
#include <stdlib.h>
#include "EbDefinitions.h"
#include "EbThreads.h"
 /****************************************
  * Win32 Includes
  ****************************************/
#ifdef _WIN32
#include <windows.h>
#elif defined(__linux__) || defined(__APPLE__)
#include <stdio.h>
#include <errno.h>
#include <fcntl.h>
#include <pthread.h>
#include <semaphore.h>
#include <time.h>
#include <unistd.h>
#else
#error OS/Platform not supported.
#endif // _WIN32
#if PRINTF_TIME
#ifdef _WIN32
#include <time.h>
void printfTime(const char *fmt, ...)
{
    va_list args;
    va_start(args, fmt);
    printf("  [%i ms]\t", ((int32_t)clock()));
    vprintf(fmt, args);
    va_end(args);
}
#endif
#endif

/****************************************
 * eb_create_thread
 ****************************************/
EbHandle eb_create_thread(
    void *thread_function(void *),
    void *thread_context)
{
    EbHandle thread_handle = NULL;

#ifdef _WIN32

    thread_handle = (EbHandle)CreateThread(
        NULL,                           // default security attributes
        0,                              // default stack size
        (LPTHREAD_START_ROUTINE)thread_function, // function to be tied to the new thread
        thread_context,                  // context to be tied to the new thread
        0,                              // thread active when created
        NULL);                          // new thread ID

#elif defined(__linux__) || defined(__APPLE__)

    pthread_attr_t attr;
    struct sched_param param = {
        .sched_priority = 99
    };
    pthread_attr_init(&attr);
    pthread_attr_setschedpolicy(&attr, SCHED_FIFO);
    pthread_attr_setschedparam(&attr, &param);

    pthread_attr_setinheritsched(&attr, PTHREAD_EXPLICIT_SCHED);

    thread_handle = (pthread_t*)malloc(sizeof(pthread_t));
    int32_t ret = pthread_create(
        (pthread_t*)thread_handle,      // Thread handle
        &attr,                       // attributes
        thread_function,                 // function to be run by new thread
        thread_context);

    if (ret != 0)
        if (ret == EPERM) {

            pthread_cancel(*((pthread_t*)thread_handle));
            pthread_join(*((pthread_t*)thread_handle), NULL);
            free(thread_handle);

            thread_handle = (pthread_t*)malloc(sizeof(pthread_t));

            pthread_create(
                (pthread_t*)thread_handle,      // Thread handle
                (const pthread_attr_t*)EB_NULL,                        // attributes
                thread_function,                 // function to be run by new thread
                thread_context);
        }

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
//#elif defined(__linux__) || defined(__APPLE__)
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
//#elif defined(__linux__) || defined(__APPLE__)
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
EbErrorType eb_destroy_thread(
    EbHandle thread_handle)
{
    EbErrorType error_return = EB_ErrorNone;

#ifdef _WIN32
    error_return = TerminateThread((HANDLE)thread_handle, 0) ? EB_ErrorDestroyThreadFailed : EB_ErrorNone;
#elif defined(__linux__) || defined(__APPLE__)
    error_return = pthread_cancel(*((pthread_t*)thread_handle)) ? EB_ErrorDestroyThreadFailed : EB_ErrorNone;
    pthread_join(*((pthread_t*)thread_handle), NULL);
    free(thread_handle);
#endif // _WIN32

    return error_return;
}
#if defined(__APPLE__)
static int32_t semaphore_id(void)
{
    static unsigned id = 0;
    return id++;
}
#endif

/***************************************
 * eb_create_semaphore
 ***************************************/
EbHandle eb_create_semaphore(uint32_t initial_count, uint32_t max_count)
{
#ifdef _WIN32
    EbHandle semaphore_handle = NULL;
    (void)max_count;

    semaphore_handle = (EbHandle)CreateSemaphore(
        NULL,                           // default security attributes
        initial_count,                   // initial semaphore count
        max_count,                       // maximum semaphore count
        NULL);                          // semaphore is not named
    return semaphore_handle;

#elif defined(__linux__) 
    EbHandle semaphore_handle = NULL;
    (void)max_count;

    semaphore_handle = (sem_t*)malloc(sizeof(sem_t));
    sem_init(
        (sem_t*)semaphore_handle,       // semaphore handle
        0,                              // shared semaphore (not local)
        initial_count);                  // initial count
    return semaphore_handle;

#elif defined(__APPLE__)
    UNUSED(max_count);
    char name[15];
    sprintf(name, "/sem_%05d_%03d", getpid(), semaphore_id());
    sem_t *s = sem_open(name, O_CREAT | O_EXCL, 0644, initial_count);
    if (s == SEM_FAILED) {
        fprintf(stderr, "errno: %d\n", errno);
        return NULL;
    }
    return s;
#endif // _WIN32

}

/***************************************
 * eb_post_semaphore
 ***************************************/
EbErrorType eb_post_semaphore(
    EbHandle semaphore_handle)
{
    EbErrorType return_error = EB_ErrorNone;

#ifdef _WIN32
    return_error = ReleaseSemaphore(
        semaphore_handle,    // semaphore handle
        1,                  // amount to increment the semaphore
        NULL)               // pointer to previous count (optional)
        ? EB_ErrorSemaphoreUnresponsive : EB_ErrorNone;
#elif defined(__linux__) || defined(__APPLE__)
    return_error = sem_post((sem_t*)semaphore_handle) ? EB_ErrorSemaphoreUnresponsive : EB_ErrorNone;
#endif // _WIN32

    return return_error;
}

/***************************************
 * eb_block_on_semaphore
 ***************************************/
EbErrorType eb_block_on_semaphore(
    EbHandle semaphore_handle)
{
    EbErrorType return_error = EB_ErrorNone;

#ifdef _WIN32
    return_error = WaitForSingleObject((HANDLE)semaphore_handle, INFINITE) ? EB_ErrorSemaphoreUnresponsive : EB_ErrorNone;
#elif defined(__linux__) || defined(__APPLE__)
    return_error = sem_wait((sem_t*)semaphore_handle) ? EB_ErrorSemaphoreUnresponsive : EB_ErrorNone;
#endif // _WIN32

    return return_error;
}

/***************************************
 * eb_destroy_semaphore
 ***************************************/
EbErrorType eb_destroy_semaphore(EbHandle semaphore_handle)
{
    EbErrorType return_error = EB_ErrorNone;

#ifdef _WIN32
    return_error = CloseHandle((HANDLE)semaphore_handle) ? EB_ErrorDestroySemaphoreFailed : EB_ErrorNone;
#elif defined(__linux__) 
    return_error = sem_destroy((sem_t*)semaphore_handle) ? EB_ErrorDestroySemaphoreFailed : EB_ErrorNone;
    free(semaphore_handle);
#elif defined(__APPLE__)
    return_error = sem_close(semaphore_handle);
#endif // _WIN32

    return return_error;
}
/***************************************
 * eb_create_mutex
 ***************************************/
EbHandle eb_create_mutex(
    void)
{
    EbHandle mutex_handle = NULL;

#ifdef _WIN32
    mutex_handle = (EbHandle)CreateMutex(
        NULL,                   // default security attributes
        FALSE,                  // FALSE := not initially owned
        NULL);                  // mutex is not named

#elif defined(__linux__) || defined(__APPLE__)

    mutex_handle = (EbHandle)malloc(sizeof(pthread_mutex_t));

    pthread_mutex_init(
        (pthread_mutex_t*)mutex_handle,
        NULL);                  // default attributes

#endif // _WIN32

    return mutex_handle;
}

/***************************************
 * EbPostMutex
 ***************************************/
EbErrorType eb_release_mutex(
    EbHandle mutex_handle)
{
    EbErrorType return_error = EB_ErrorNone;

#ifdef _WIN32
    return_error = ReleaseMutex((HANDLE)mutex_handle) ? EB_ErrorCreateMutexFailed : EB_ErrorNone;
#elif defined(__linux__) || defined(__APPLE__)
    return_error = pthread_mutex_unlock((pthread_mutex_t*)mutex_handle) ? EB_ErrorCreateMutexFailed : EB_ErrorNone;
#endif // _WIN32

    return return_error;
}

/***************************************
 * eb_block_on_mutex
 ***************************************/
EbErrorType eb_block_on_mutex(
    EbHandle mutex_handle)
{
    EbErrorType return_error = EB_ErrorNone;

#ifdef _WIN32
    return_error = WaitForSingleObject((HANDLE)mutex_handle, INFINITE) ? EB_ErrorMutexUnresponsive : EB_ErrorNone;
#elif defined(__linux__) || defined(__APPLE__)
    return_error = pthread_mutex_lock((pthread_mutex_t*)mutex_handle) ? EB_ErrorMutexUnresponsive : EB_ErrorNone;
#endif // _WIN32

    return return_error;
}

/***************************************
 * eb_block_on_mutex_timeout
 ***************************************/
EbErrorType eb_block_on_mutex_timeout(
    EbHandle mutex_handle,
    uint32_t    timeout)
{
    EbErrorType return_error = EB_ErrorNone;

#ifdef _WIN32
    WaitForSingleObject((HANDLE)mutex_handle, timeout);
#elif defined(__linux__) || defined(__APPLE__)
    return_error = pthread_mutex_lock((pthread_mutex_t*)mutex_handle) ? EB_ErrorMutexUnresponsive : EB_ErrorNone;
    (void)timeout;
#endif // _WIN32

    return return_error;
}

/***************************************
 * eb_destroy_mutex
 ***************************************/
EbErrorType eb_destroy_mutex(
    EbHandle mutex_handle)
{
    EbErrorType return_error = EB_ErrorNone;

#ifdef _WIN32
    return_error = CloseHandle((HANDLE)mutex_handle) ? EB_ErrorDestroyMutexFailed : EB_ErrorNone;
#elif defined(__linux__) || defined(__APPLE__)
    return_error = pthread_mutex_destroy((pthread_mutex_t*)mutex_handle) ? EB_ErrorDestroyMutexFailed : EB_ErrorNone;
    free(mutex_handle);
#endif // _WIN32

    return return_error;
}
