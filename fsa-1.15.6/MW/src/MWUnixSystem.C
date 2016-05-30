
#include "MWSystem.h"

#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>

#include <sys/resource.h> // getrusage
#include <sys/time.h> // gettimeofday

int
MWSystem::getpid() {
	return ::getpid();
}

void
MWSystem::gethostname(char *h, int len) {
	::gethostname(h, len);
}

double 
timeval_to_double( struct timeval t )
{
	return (double) t.tv_sec + ( (double) t.tv_usec * (double) 1e-6 );
}

double
MWSystem::gettimeofday() {
	struct timeval t;
	::gettimeofday(&t, NULL);
	return timeval_to_double(t);
}

double
MWSystem::getcputime() {
	struct rusage rt;
	::getrusage ( RUSAGE_SELF, &rt );
	double d =  timeval_to_double ( rt.ru_utime ) + timeval_to_double ( rt.ru_stime );
	::getrusage ( RUSAGE_CHILDREN, &rt );
	return d +  timeval_to_double ( rt.ru_utime ) + timeval_to_double ( rt.ru_stime );
}

void
MWSystem::sleep(int t) {
	::sleep(t);
}

void
MWSystem::mkdir(char *c) {
	::mkdir(c, 0755);
}
