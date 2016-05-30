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
#ifndef MWWORKER_H
#define MWWORKER_H

#include <stdlib.h>
#include <string.h>
#include "MWTask.h"

/**  This is the worker class that performs the tasks in the
    opportunistic condor environment.  It is an oppressed worker class
    in that in simply executes the tasks given to it and reports the
    results back to the master.  {\tt :-)}

    Capitalist stooges who wish to create an application
    must derive a class from this class, and implement the following
    two methods

    \begin{itemize}
    \item unpack_init_data()
    \item execute_task()
    \end{itemize}

    @see MWDriver
    @see MWWorker
    @author Mike Yoder, modified by Jeff Linderoth and Jean-Pierre Goux */

class MWWorker
{

public:
    
  /// Default constructor
  MWWorker();  

  /// Default Destructor
  virtual ~MWWorker();

  /// Giddyap!
  void go( int argc, char *argv[] );

  MWReturn do_benchmark_cmd( );

  MWReturn worker_mainloop_oneshot ();


	  /// Our RM / Comm class.  Used here only for communication.
  static MWRMComm * RMC;

 protected:

  /// The task ID of the master - used for sending messages.
  int master;

  /// The task instance that a worker will use for packing/unpacking 
  /// information from the master
  MWTask *workingTask;
  MWTask *controlTask;

  /// The name of the machine the worker is running on.
  char mach_name[64];

  /// The name of the master machine.
  char master_mach_name[64];

  /** 
      Here we might in the future pack some useful information about
      the specific machine on which we're running.  Right now,
      all workers are equal, and we pass only the hostname.

      There must be a "matching" unpack_worker_initinfo() in
      the MWDriver class.
  */
  virtual void pack_worker_initinfo() {};


//virtual MWTask* gimme_a_task() { return ((MWTask *) 0xdeadbeef); }

  /**
     This unpacks the initial data that is sent to the worker
     once the master knows that he has started.

     There must be a "matching" pack_worker_init_data() in
     the MWDriver class derived for your application.
     
   */
  virtual MWReturn unpack_init_data() = 0;

  /** This function performs the actions that happen
      once the Worker finds out there is a task to do.
      You will need to cast the MWTask * to a pointer of the Task type
      derived for your application.  For example

      \begin{verbatim}
      Task_Fib *dt = dynamic_cast<Task_Fib *> ( t );
      assert( dt );     
      \end{verbatim}    

   */
  virtual void execute_task( MWTask * ) = 0;

  /**
     This function performs the action that happens once the worker
	 receives a subtask. Subtask id starts at 0 and ends at n-1 subtasks.
	
	Deprecated - changing subtask model for container of tasks model
   */
  //virtual void execute_subtask(MWTask *, int) {}

  /**
     The number of subtasks per worker is defined here.
	 If this function is not implemented, then data streaming/subtasks
	 is turned off and Master-Worker executes normally.

	    Deprecated - changing subtask model for container of tasks model
   */
  //virtual void set_num_subtask(int *) {}

  /** Run a benchmark, given an MWTASK.  The default implementation 
   *  is to call \begin{verbatim} execute_task(t) \end{verbatim}
   *  and return 1.0/(Task Time) as a benchmark of how fast the machine
   *  is
   */
  virtual double benchmark ( MWTask *t );

  /**
    If you have some driver data that you would like to use on the
    worker in order to execute the task, you should unpack it here.
   */
  virtual void unpack_driver_task_data( void ) {};
  
  /**
     This is run before the worker_mainloop().  It prints a message 
     that the worker has been spawned, sends the master an INIT
     message, and awaits a reply from the master, upon which 
     it calls  unpack_init_data().
    
  */
  MWReturn greet_master(); 

  /**
     This sits in a loop in which it asks for work from the master
     and does the work.  Maybe we should name this class GradStudent.
     {\tt :-)}
   */
  void worker_mainloop();

  /// Die!!!!!!
  void suicide();

  private:

public:
    /* send a update message to the master, without stopping current execution
	 * assume that the user has already called RMC->initsend, and RMC->pack */
  	int send_update_message();

	/* Non-blocking check for DRIVER_UPDATE message, calls unpack_update */
	void check_update_message();

	/* Unpack data sent by send_update in master. Defined by user and called by check_update_message */
	virtual void unpack_update(){};

	int workerid; // id of worker.

};

#endif


