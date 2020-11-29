#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include <cpuinfo.h>
#include <x86/api.h>
#include <cpuinfo/internal-api.h>
#include <cpuinfo/log.h>



void cpuinfo_x86_linux_init(void) {
	struct cpuinfo_x86_processor x86_processor;
	memset(&x86_processor, 0, sizeof(x86_processor));
	cpuinfo_x86_init_processor(&x86_processor);
}
