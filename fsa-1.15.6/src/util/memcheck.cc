
/**
 * \file memcheck.cc
 * This file is part of FSA, a sequence alignment algorithm.
 * \author Source code in this file was written by Colin Dewey.
 */

#include "config.h"
#include "memcheck.h"

#ifdef HAVE_STRUCT_SYSINFO_TOTALRAM
// Memory check for Linux systems
#include <sys/sysinfo.h>
int total_ram() {
	struct sysinfo info;
	if (sysinfo(&info) == 0) {
		unsigned long totalram = info.totalram;
#ifdef HAVE_STRUCT_SYSINFO_MEM_UNIT
		// info.totalram is in units of info.mem_unit bytes
		totalram *= info.mem_unit;
#endif
		// Return total RAM in megabytes
		return totalram >> 20;
	} else {
		return -1;
	}
}
#elif HAVE_DECL_SYSCTL && HAVE_DECL_CTL_HW && HAVE_DECL_HW_PHYSMEM
// Memory check for BSD/Darwin systems
#include <sys/sysctl.h>
int total_ram() {
	// Management Information Base (MIB) codes for physical memory
	static int mib[] = {CTL_HW, HW_PHYSMEM};

	// Variable for result of physical memory query and its size
	size_t physmem;
	size_t len = sizeof(physmem);

	if (sysctl(mib, sizeof(mib) / sizeof(mib[0]), &physmem, &len, NULL, 0) == 0
		and len == sizeof(physmem)) {
		// sysctl returns total RAM in bytes
		return physmem >> 20;
	} else {
		return -1;
	}
}
#else
int total_ram() {
	return -1;
}
#endif
