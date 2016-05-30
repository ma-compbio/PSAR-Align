/***************************Copyright-DO-NOT-REMOVE-THIS-LINE**
  *
  * Condor Software Copyright Notice
  * Copyright (C) 1990-2004, Condor Team, Computer Sciences Department,
  * University of Wisconsin-Madison, WI.
  *
  * This source code is covered by the Condor Public License, which can
  * be found in the accompanying LICENSE.TXT file, or online at
  * www.condorproject.org.
  *
  * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  * AND THE UNIVERSITY OF WISCONSIN-MADISON "AS IS" AND ANY EXPRESS OR
  * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  * WARRANTIES OF MERCHANTABILITY, OF SATISFACTORY QUALITY, AND FITNESS
  * FOR A PARTICULAR PURPOSE OR USE ARE DISCLAIMED. THE COPYRIGHT
  * HOLDERS AND CONTRIBUTORS AND THE UNIVERSITY OF WISCONSIN-MADISON
  * MAKE NO MAKE NO REPRESENTATION THAT THE SOFTWARE, MODIFICATIONS,
  * ENHANCEMENTS OR DERIVATIVE WORKS THEREOF, WILL NOT INFRINGE ANY
  * PATENT, COPYRIGHT, TRADEMARK, TRADE SECRET OR OTHER PROPRIETARY
  * RIGHT.
  *
  ****************************Copyright-DO-NOT-REMOVE-THIS-LINE**/
/* MWprintf.C - The place for the MWPrintf stuff. */

#include "MW.h"
#include <stdarg.h>
#include <time.h>
#include <string.h>
#include <errno.h>

#ifndef WINDOWS
#include <sys/time.h>
#endif

int MWprintf_level = 50;
static int fopen_count = 0;

FILE* Open ( char *filename, char *mode )
{
    FILE *fp;

    fp = fopen ( filename, mode );
    if ( fp ) { 
      fopen_count++;
    }
    else{
      //MWprintf ( 10, "Could not open file %s as errno is %d\n", filename, errno );
    }
    return fp;
}

void Close ( FILE *fp )
{
    int retval;
	fflush(fp);
    retval = fclose ( fp );
    if ( retval != 0 )
    	MWprintf ( 10, "Could not close the file %x as errno is %d\n", fp, errno );
    else
    	fopen_count--;
}

int 
get_MWprintf_level () {
   return MWprintf_level;
}


int set_MWprintf_level ( int level ) {
	int foo = MWprintf_level;
	if ( (level<0) || (level>99) ) {
		MWprintf( 10, "Bad arg \"%d\" in set_MWprintf_level().\n", level );
	} else {
		MWprintf_level = level;
	}
	return foo;
}

void MWprintf ( int level, char *fmt, ... ) {

	static int printTime = TRUE;

	if ( level > MWprintf_level ) {
		return;
	}

#ifndef WINDOWS
	if ( printTime ) {
		struct timeval tv;
		::gettimeofday(&tv, NULL);
		struct tm *t = localtime(&tv.tv_sec);
		printf( "%d:%02d:%02d.%03d ",  t->tm_hour, t->tm_min, t->tm_sec, (int) (tv.tv_usec / 1000));
	}
#endif

	va_list ap;
	va_start( ap, fmt );
	vprintf( fmt, ap );
	va_end( ap );
	fflush( stdout );

	if ( fmt[strlen(fmt)-1] == '\n' ) 
		printTime = TRUE;
	else
		printTime = FALSE;
}
