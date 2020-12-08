#include <stdbool.h>
#include <stddef.h>

#include <cpuinfo.h>
#include <cpuinfo/internal-api.h>
#include <cpuinfo/log.h>

#ifdef __linux__
    #include <unistd.h>
    #include <sys/syscall.h>
    #if !defined(__NR_getcpu)
        #include <asm-generic/unistd.h>
    #endif
#endif

bool cpuinfo_is_initialized = false;
