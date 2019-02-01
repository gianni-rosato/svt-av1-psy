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
 * EbCreateThread
 ****************************************/
EbHandle EbCreateThread(
    void *threadFunction(void *),
    void *threadContext)
{
    EbHandle threadHandle = NULL;

#ifdef _WIN32

    threadHandle = (EbHandle)CreateThread(
        NULL,                           // default security attributes
        0,                              // default stack size
        (LPTHREAD_START_ROUTINE)threadFunction, // function to be tied to the new thread
        threadContext,                  // context to be tied to the new thread
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

    threadHandle = (pthread_t*)malloc(sizeof(pthread_t));
    int32_t ret = pthread_create(
        (pthread_t*)threadHandle,      // Thread handle
        &attr,                       // attributes
        threadFunction,                 // function to be run by new thread
        threadContext);

    if (ret != 0)
        if (ret == EPERM) {

            pthread_cancel(*((pthread_t*)threadHandle));
            pthread_join(*((pthread_t*)threadHandle), NULL);
            free(threadHandle);

            threadHandle = (pthread_t*)malloc(sizeof(pthread_t));

            pthread_create(
                (pthread_t*)threadHandle,      // Thread handle
                (const pthread_attr_t*)EB_NULL,                        // attributes
                threadFunction,                 // function to be run by new thread
                threadContext);
        }

#endif // _WIN32

    return threadHandle;
}

///****************************************
// * EbStartThread
// ****************************************/
//EbErrorType EbStartThread(
//    EbHandle threadHandle)
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
//    //error_return = ResumeThread((HANDLE) threadHandle) ? EB_ErrorThreadUnresponsive : EB_ErrorNone;
//#elif defined(__linux__) || defined(__APPLE__)
//#endif // _WIN32
//
//    error_return = (threadHandle) ? EB_ErrorNone : EB_ErrorNullThread;
//
//    return error_return;
//}
//
///****************************************
// * EbStopThread
// ****************************************/
//EbErrorType EbStopThread(
//    EbHandle threadHandle)
//{
//    EbErrorType error_return = EB_ErrorNone;
//
//#ifdef _WIN32
//    //error_return = SuspendThread((HANDLE) threadHandle) ? EB_ErrorThreadUnresponsive : EB_ErrorNone;
//#elif defined(__linux__) || defined(__APPLE__)
//#endif // _WIN32
//
//    error_return = (threadHandle) ? EB_ErrorNone : EB_ErrorNullThread;
//
//    return error_return;
//}
//
/****************************************
 * EbDestroyThread
 ****************************************/
EbErrorType EbDestroyThread(
    EbHandle threadHandle)
{
    EbErrorType error_return = EB_ErrorNone;

#ifdef _WIN32
    error_return = TerminateThread((HANDLE)threadHandle, 0) ? EB_ErrorDestroyThreadFailed : EB_ErrorNone;
#elif defined(__linux__) || defined(__APPLE__)
    error_return = pthread_cancel(*((pthread_t*)threadHandle)) ? EB_ErrorDestroyThreadFailed : EB_ErrorNone;
    pthread_join(*((pthread_t*)threadHandle), NULL);
    free(threadHandle);
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
 * EbCreateSemaphore
 ***************************************/
EbHandle EbCreateSemaphore(uint32_t initialCount, uint32_t maxCount)
{
#ifdef _WIN32
    EbHandle semaphoreHandle = NULL;
    (void)maxCount;

    semaphoreHandle = (EbHandle)CreateSemaphore(
        NULL,                           // default security attributes
        initialCount,                   // initial semaphore count
        maxCount,                       // maximum semaphore count
        NULL);                          // semaphore is not named
    return semaphoreHandle;

#elif defined(__linux__) 
    EbHandle semaphoreHandle = NULL;
    (void)maxCount;

    semaphoreHandle = (sem_t*)malloc(sizeof(sem_t));
    sem_init(
        (sem_t*)semaphoreHandle,       // semaphore handle
        0,                              // shared semaphore (not local)
        initialCount);                  // initial count
    return semaphoreHandle;

#elif defined(__APPLE__)
    UNUSED(maxCount);
    char name[15];
    sprintf(name, "/sem_%05d_%03d", getpid(), semaphore_id());
    sem_t *s = sem_open(name, O_CREAT | O_EXCL, 0644, initialCount);
    if (s == SEM_FAILED) {
        fprintf(stderr, "errno: %d\n", errno);
        return NULL;
    }
    return s;
#endif // _WIN32

}

/***************************************
 * EbPostSemaphore
 ***************************************/
EbErrorType EbPostSemaphore(
    EbHandle semaphoreHandle)
{
    EbErrorType return_error = EB_ErrorNone;

#ifdef _WIN32
    return_error = ReleaseSemaphore(
        semaphoreHandle,    // semaphore handle
        1,                  // amount to increment the semaphore
        NULL)               // pointer to previous count (optional)
        ? EB_ErrorSemaphoreUnresponsive : EB_ErrorNone;
#elif defined(__linux__) || defined(__APPLE__)
    return_error = sem_post((sem_t*)semaphoreHandle) ? EB_ErrorSemaphoreUnresponsive : EB_ErrorNone;
#endif // _WIN32

    return return_error;
}

/***************************************
 * EbBlockOnSemaphore
 ***************************************/
EbErrorType EbBlockOnSemaphore(
    EbHandle semaphoreHandle)
{
    EbErrorType return_error = EB_ErrorNone;

#ifdef _WIN32
    return_error = WaitForSingleObject((HANDLE)semaphoreHandle, INFINITE) ? EB_ErrorSemaphoreUnresponsive : EB_ErrorNone;
#elif defined(__linux__) || defined(__APPLE__)
    return_error = sem_wait((sem_t*)semaphoreHandle) ? EB_ErrorSemaphoreUnresponsive : EB_ErrorNone;
#endif // _WIN32

    return return_error;
}

/***************************************
 * EbDestroySemaphore
 ***************************************/
EbErrorType EbDestroySemaphore(EbHandle semaphoreHandle)
{
    EbErrorType return_error = EB_ErrorNone;

#ifdef _WIN32
    return_error = CloseHandle((HANDLE)semaphoreHandle) ? EB_ErrorDestroySemaphoreFailed : EB_ErrorNone;
#elif defined(__linux__) 
    return_error = sem_destroy((sem_t*)semaphoreHandle) ? EB_ErrorDestroySemaphoreFailed : EB_ErrorNone;
    free(semaphoreHandle);
#elif defined(__APPLE__)
    return_error = sem_close(semaphoreHandle);
#endif // _WIN32

    return return_error;
}
/***************************************
 * EbCreateMutex
 ***************************************/
EbHandle EbCreateMutex(
    void)
{
    EbHandle mutexHandle = NULL;

#ifdef _WIN32
    mutexHandle = (EbHandle)CreateMutex(
        NULL,                   // default security attributes
        FALSE,                  // FALSE := not initially owned
        NULL);                  // mutex is not named

#elif defined(__linux__) || defined(__APPLE__)

    mutexHandle = (EbHandle)malloc(sizeof(pthread_mutex_t));

    pthread_mutex_init(
        (pthread_mutex_t*)mutexHandle,
        NULL);                  // default attributes

#endif // _WIN32

    return mutexHandle;
}

/***************************************
 * EbPostMutex
 ***************************************/
EbErrorType EbReleaseMutex(
    EbHandle mutexHandle)
{
    EbErrorType return_error = EB_ErrorNone;

#ifdef _WIN32
    return_error = ReleaseMutex((HANDLE)mutexHandle) ? EB_ErrorCreateMutexFailed : EB_ErrorNone;
#elif defined(__linux__) || defined(__APPLE__)
    return_error = pthread_mutex_unlock((pthread_mutex_t*)mutexHandle) ? EB_ErrorCreateMutexFailed : EB_ErrorNone;
#endif // _WIN32

    return return_error;
}

/***************************************
 * EbBlockOnMutex
 ***************************************/
EbErrorType EbBlockOnMutex(
    EbHandle mutexHandle)
{
    EbErrorType return_error = EB_ErrorNone;

#ifdef _WIN32
    return_error = WaitForSingleObject((HANDLE)mutexHandle, INFINITE) ? EB_ErrorMutexUnresponsive : EB_ErrorNone;
#elif defined(__linux__) || defined(__APPLE__)
    return_error = pthread_mutex_lock((pthread_mutex_t*)mutexHandle) ? EB_ErrorMutexUnresponsive : EB_ErrorNone;
#endif // _WIN32

    return return_error;
}

/***************************************
 * EbBlockOnMutexTimeout
 ***************************************/
EbErrorType EbBlockOnMutexTimeout(
    EbHandle mutexHandle,
    uint32_t    timeout)
{
    EbErrorType return_error = EB_ErrorNone;

#ifdef _WIN32
    WaitForSingleObject((HANDLE)mutexHandle, timeout);
#elif defined(__linux__) || defined(__APPLE__)
    return_error = pthread_mutex_lock((pthread_mutex_t*)mutexHandle) ? EB_ErrorMutexUnresponsive : EB_ErrorNone;
    (void)timeout;
#endif // _WIN32

    return return_error;
}

/***************************************
 * EbDestroyMutex
 ***************************************/
EbErrorType EbDestroyMutex(
    EbHandle mutexHandle)
{
    EbErrorType return_error = EB_ErrorNone;

#ifdef _WIN32
    return_error = CloseHandle((HANDLE)mutexHandle) ? EB_ErrorDestroyMutexFailed : EB_ErrorNone;
#elif defined(__linux__) || defined(__APPLE__)
    return_error = pthread_mutex_destroy((pthread_mutex_t*)mutexHandle) ? EB_ErrorDestroyMutexFailed : EB_ErrorNone;
    free(mutexHandle);
#endif // _WIN32

    return return_error;
}
