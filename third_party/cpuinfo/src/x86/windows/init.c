#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include <cpuinfo.h>
#include <x86/api.h>
#include <cpuinfo/internal-api.h>
#include <cpuinfo/log.h>

#include <windows.h>

#ifdef __GNUC__
  #define CPUINFO_ALLOCA __builtin_alloca
#else
  #define CPUINFO_ALLOCA _alloca
#endif



BOOL CALLBACK cpuinfo_x86_windows_init(PINIT_ONCE init_once, PVOID parameter, PVOID* context) {
    struct cpuinfo_x86_processor x86_processor;
    ZeroMemory(&x86_processor, sizeof(x86_processor));
    cpuinfo_x86_init_processor(&x86_processor);

    return TRUE;
}
