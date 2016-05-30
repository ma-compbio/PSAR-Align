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
#ifndef MW_H
#define MW_H


#include "MWprintf.h"

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>

/// FALSE is defined as 0
#ifndef FALSE
#define FALSE 0
#endif

/// TRUE is defined as 1
#ifndef TRUE
#define TRUE 1
#endif

/// UNDEFINED is defined as -1
#ifndef UNDEFINED
#define UNDEFINED -1
#endif

/** An upper limit on the number of workers -- so we can allocated arrays
    to keep statistics */
const int MW_MAX_WORKERS = 8192;

/**@name Introduction.

   MW is a class library that can be a useful tool for building
   opportunistic, fault-tolerant applications for high throughput
   computing.  

   In order to build an application, there are three classes
   that the user {\rm must} rederive.

   \begin{itemize}
   \item \Ref{MWDriver}
   \item \Ref{MWWorker}
   \item \Ref{MWTask}
   \end{itemize}

   The documentation of these classes includes a description 
   of the pure virtual methods that must be implemented for 
   a user's particular application.

   Using the MW library allows the user to focus on the application
   specific implementation at hand.  All details related to
   being fault tolerant and opportunistic are implemented in the
   MW library.

   Also included is a small, naive, example of how to create
   an application with the MW class library.  The three classes

   \begin{itemize}
   \item \Ref{Driver_Fib}
   \item \Ref{Worker_Fib}
   \item \Ref{Task_Fib}
   \end{itemize}
   
   are concrete classes derived from MW's abstract classes.  
   Using these classes, a simple program that makes a lot of 
   fibonacci sequences is created.

 */

/**
   A list of the message tags that the Master-Worker application
   will send and receive.  See the switch statement in master_mainloop 
   for more information.
 */

enum MWmessages{
  HOSTADD = 33,
  HOSTDELETE,   
  HOSTSUSPEND, 
  HOSTRESUME, 
  TASKEXIT,   
  DO_THIS_WORK, 
  RESULTS,      
  INIT,         
  INIT_REPLY,   
  BENCH_RESULTS,
  KILL_YOURSELF,
  CHECKSUM_ERROR,
  RE_INIT,
  REFRESH,
  SUBRESULTS,
  /* added for async message: */
  UPDATE_FROM_WORKER, /* worker send to driver a update message */
  UPDATE_FROM_DRIVER,  /* driver send to worker a new message */
  NO_MESSAGE,
  STOP_WORK 
};

/** Possible return values from many virtual functions */
enum MWReturn {
		/// Normal return
	OK, 
		/// We want to exit, not an error.
	QUIT, 
		/// We want to exit immediately; a bad error ocurred
	ABORT,
		/// by jae
		/// Some job is still running
	PENDING
};

typedef enum { MWNORMAL, MWNWS, MWNUMTASKTYPES } MWTaskType;

#endif
