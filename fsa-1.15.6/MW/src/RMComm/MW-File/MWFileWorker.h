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
#ifndef MWFILEWORKER_H
#define MWFILEWORKER_H

enum FileWorkerState
{ 
    FILE_FREE,
    FILE_SUBMITTED,
    FILE_EXECUTE,	// Started executing for the first time.
    FILE_RUNNING, 
    FILE_KILLED, 
    FILE_SUSPENDED,
    FILE_RESUMED,
    FILE_IDLE,
    FILE_TRANSIT,
};


struct FileWorker
{
    int id;

    // What is the message number that I have to look next for.
    int counter;

    // What is the message that the worker is looking for.
    int worker_counter;

    // What is the condor_id of the worker.
    int condorID;

    // What is the proc id.
    int condorprocID;

    FileWorkerState state;

    int arch;

    int event_no;

    /** What was last served */
    int served;

	int exec_class;
};

#endif
