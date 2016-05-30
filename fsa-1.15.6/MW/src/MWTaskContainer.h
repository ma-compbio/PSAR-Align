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

#ifndef MWTASKCONTAINER_H
#define MWTASKCONTAINER_H

#include "MW.h"
#include "MWList.h"
#include "MWTask.h"
//#include <MWRMComm.h> MWTask.h already had this header

//forward declarations necessary for compilation
class MWTask;
class MWWorkerID;
//enum MWTaskType;

/**
	MWTaskContainer allows MWDriver to send multiple tasks to a single worker.
	The exact number of tasks is hardcoded in send_task_to_worker() in 
	MWDriver.C
	The worker then works on each task and sends an update to MWDriver upon the 
	completion of each task; act_on_completed_subtask(MWTask*) is called.
	When all tasks are done, the results of all the tasks are sent back to 
	MWDriver and act_on_completed_task(MWTask*) is called for each completed 
	task. 
	MWDriver then sends a new MWTaskContainer to the Worker.

	The following functions in this class have to be defined in a derived class:
	- pack_work()
	- unpack_work()
	- pack_subresults()
	- unpack_subresults()
	- pack_results()
	- unpack_results()

	The following functions in MWDriver.C have to be defined:
	- MWDriver::act_on_completed_subtask(MWTask*)
	
	Notes:	
	- None of the functions in MWTask are called directly by MW. They are now
	called by MWTaskContainer.
*/
class MWTaskContainer
{
public:
	MWTaskType taskType;
	int number; // have to define numbering system for task containers, let's make it the number of the first task for now
	MWList<MWTask> *tasks;
	MWWorkerID *worker;
	MWRMComm *RMC; //pointer to RMC instance in MWTask.h

	MWTaskContainer();
	virtual ~MWTaskContainer();
	void addTask(MWTask *t);
	void removeAll();
	MWTask* First();
	bool AfterEnd();
	MWTask* Current();
	MWTask* Next();

	/// Pack all the tasks in this container
	virtual void pack_work();

	/// Unpack all the tasks in this container
	virtual void unpack_work();

	/// Pack all the results in this container
	virtual void pack_results();

	/// Unpack all the results in this container
	virtual void unpack_results();

	/// Called by MWWorker at the completion of each task(not task container). This could be left empty if we only want to let MWDriver know that a task was completed
	virtual void pack_subresults(int tasknum);

	/// Called by MWDriver at the completion of each task. @param tasknum number of the task that was completed
	virtual void unpack_subresults(int tasknum);

	virtual void printself(int level);

	int FirstNum();

	int LastNum();

	/// Number of tasks in this task container
	int Number();

private:
	double Side;
};

#endif
