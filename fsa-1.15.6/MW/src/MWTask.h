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
#ifndef MWTASK_H
#define MWTASK_H

#include "MW.h"
#include "MWWorkerID.h"
#include "MWGroup.h"
#include <MWRMComm.h>
#include <stdio.h>

//I'm moving this into MW.h so that task container can compile
//typedef enum { MWNORMAL, MWNWS, MWNUMTASKTYPES } MWTaskType;

/** 
    This the class that represents a unit of work. 

    The task consits of two main components.  The "work" to be done,
    and the "result" of the completed task.
    In order to create an application, the user must specify methods
    for packing and unpacking both the "work" and "result" portions
    of the task.
    The user is not responsible for initializing the send
    or for packing the task's number    

    When the task is being serviced by a worker, it contains a link
    to the ID of that instance of the worker.


    @see MWDriver
    @see MWWorker
    @author Mike Yoder, modified by Jeff Linderoth and Jean-Pierre Goux
*/

class MWTask
{
 public:   

  /// Default constructor
  MWTask();
  
  /// Default Destructor
  virtual ~MWTask();

  /// The task's number
  int number;
int randomstop;
int numsubtask;
  /// The task's type.
  MWTaskType taskType;

	  /** @name Time and usage data */
	  //@{
	  /** The time (wall clock) that it takes to run the 'execute_task'
		  function on the worker's side. */
  double working_time;
	  /** The amount of user+system time taken by this task, measured
		  from start to finish of the 'execute_task' function. */
  double cpu_time;

	  //@}

	  /**@name Packing and Unpacking 

		 The user must pack and unpack the contents of the class
		 so that it can be sent to the worker and back.  These 
		 functions will make RMC->pack() and RMC->unpack() calls.
   */
  //@{

  /// Pack the work portion of the task
  virtual void pack_work( void ) = 0;

  /// Unpack the work portion of the task
  virtual void unpack_work( void ) = 0;

  /// Pack the result portion of the task
  virtual void pack_results( void ) = 0;
  
  /// Unpack the result portion of the task
  virtual void unpack_results( void ) = 0;

  /// Pack the subresult portion of the task
  virtual void pack_subresults( int ) {}
  
  /// Unpack the subresult portion of the task
  virtual void unpack_subresults( int ) {}


  /// Dump this task to the screen
  virtual void printself( int level = 60 );

	  //@}

  /**@name Checkpointing Utilities

     These two functions are used when checkpointing.  Simply put, 
     they write/read the state of the task to this file pointer 
     provided. */

      //@{

      /** Write the state of this task out to a FILE* */
  virtual void write_ckpt_info( FILE *fp ) {};
      /** Read the state of this task from a FILE* (overwriting
          any existing state */
  virtual void read_ckpt_info( FILE *fp ) {};

  void write_group_info ( FILE *fp );

  void read_group_info ( FILE *fp );


      //@}

  void initGroups ( int num );
  /// add task to a workclass/group	
  void addGroup ( int num );
  /// remove a task from a workclass/group
  void deleteGroup ( int num );
  bool doesBelong ( int num );
  MWGroup *getGroup ( );


  /**@name List management 

     The task also has pointers that aid for managing the list
     of tasks to be done
   */

  //@{
  /// A pointer to the worker ID executing this task (NULL if none)
  MWWorkerID *worker;

  /*
  /// A pointer to the next task in the list (NULL if none)
  MWTask *next;
  */
  //@}

  MWGroup *group;

  static MWRMComm * RMC;

	//by jae
	void completedTask(double wall_time, double cpu_time);

#ifdef MEASURE
public:
  // double _measure_start_time;	// the time called give_task
  // double _measure_MP_master_time;	// master send_data + recv_result
  double _measure_MP_worker_time;	// worker recv_data + send_result
  double _measure_MP_worker_cpu_time;	// CPU time (worker recv_data + send_result)
  // double _measure_Exe_wall_time;	// worker execute_task wall time = working_time
  // double _measure_Exe_cpu_time;	// worker execute_task CPU time = cpu_time
  // double _measure_Life_time;   	// from creation to deletion
#endif // MEASURE
  
};

#endif










