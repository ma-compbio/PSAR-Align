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
#ifndef MWDRIVER_H
#define MWDRIVER_H

#include "MWTask.h"
#include "MWWorker.h"
#include "MWStats.h"
#include "MWList.h"
#include "MW.h"
#include <MWRMComm.h>
#include <pthread.h>

/// This is the "key" by which the task list can be managed.
#define MWKey double
//typedef double MWKey;

/// This is if you wish to have XML output.
//#define XML_OUTPUT
// (Note -- it should be included in 'configure' output later)

/** The ways in which tasks may be added to the list.
*/
enum MWTaskAdditionMode {
  /// Tasks will be added at the end of the list
  ADD_AT_END,
  /// Tasks will be added at the beginning
  ADD_AT_BEGIN,
  /// Tasks will be added based on their key (low keys before high keys)
  ADD_BY_KEY 
};

/** The ways in which tasks may be removed from the list.
*/
enum MWTaskRetrievalMode {
  /// Task at head of list will be returned.
  GET_FROM_BEGIN,
  /// Task with lowest key will be retrieved
  GET_FROM_KEY 
};

enum MWREFRESH_TYPE
{
	MW_ALL,
	MW_THIS,
	MW_NONE
};

/** The suspension policy to use - What do we do when it happens? */
enum MWSuspensionPolicy {
		/// Normally do nothing unless there are idle workers.
	DEFAULT,
		/** Always reassign the task; move it to the front of the 
			todo list */
	REASSIGN 
};

/** Tasks are always assigned to he first "IDLE" machine.  By ordering
    the machine list, we can implement a number of scheduling policies.
    Insert your own favorite policy here.
*/
enum MWMachineOrderingPolicy {

  /// The machines are ordered simply as they become available
  NO_ORDER,
  
  /** Machines are ordered by the "benchmark" result.  Larger benchmark results
      go first, so if using time to complete a task as the benchmark, you
      should return 1/time as the benchmark value.
  */
  BY_USER_BENCHMARK,

  /// Machines are ordered by KFLOPS reported by Condor.
  BY_KFLOPS
  
};

/* TASK CLUSTERING / WORK CYCLE */

/** This class is responsible for managing an application in an 
    opportunistic environment.  The goal is to be completely fault - 
    tolerant, dealing with all possiblities of host (worker) problems.  
    To do this, the MWDriver class manages a set of tasks and a set 
    of workers.  It monitors messages about hosts coming up and 
    going down, and assigns tasks appropriately.

	This class is built upon some sort of resource management and 
	message passing lower layer.  Previously, it was built directly 
	on top of Condor - PVM, but the interface to that has been 
	abstracted away so that it can use any facility that provides 
	for resource management and message passing.  See the abstract
	MWRMComm class for details of this lower layer.  When interfacing
	with this level, you'll have use the RMC object that's a static
	member of the MWDriver, MWTask, and MWWorker class.  

    To implement an application, a user must derive a class from this
    base class and implement the following methods:

    - get_userinfo()
    - setup_initial_tasks()
    - pack_worker_init_data()
    - act_on_completed_task()
  

	For a higher level of control regarding distribution of tasks to workers, 
	tasks and workers may be enrolled to workclasses/groups by calling the
	following methods:

	- MWWorkerID::addGroup()
	- MWWorkerID::deleteGroup()
	- MWTask::addGroup()
	- MWTask::deleteGroup()

	These methods can be called by an application using workclasses:
	
	- workClasses_set()  (Required)
	- workClasses_get()
	- workClasses_gettasks()
	- workClasses_getworkers()

	To enroll workers in workclasses, there are two options:
	- act_on_starting_worker() can be implemented to enroll workers to 
	workclasses based on MWWorkerID info at the start of execution
	- pack_worker_init_data() can be implemented to enroll workers. When used 
	in conjunction with refreshWorkers(), can change workclass of workers in the 	middle of execution.	

    Similar application dependent methods must be implemented
    for the "Task" of work to be done and the "Worker" who performs
    the tasks.
    
    @see MWTask
    @see MWWorker
    @see MWRMComm
    @author Mike Yoder, Jeff Linderoth, Jean-Pierre Goux, Sanjeev Kulkarni
*/


class MWDriver {

	friend class MWTask;

	public:
	
		/// Default constructor
		MWDriver();
	
		/** Destructor - walks through lists of tasks & workers and
		deletes them. */	
		virtual ~MWDriver();
	
		/** This method runs the entire fault-tolerant
		application in the condor environment.  What is *really* does
		is call setup_master(), then master(), then printresults(), 
		and then ends.  See the other functions	for details. 
		*/
		void go( int argc, char *argv[] );
	
		/** This version of go simply calls go(0, NULL).*/
		void go() { go ( 0, NULL ); };
  
		/**
		Prints the Results. Applications may re-implement this
		to print their application specific results.
		It is meant to be over-ridden.
		*/
		virtual void printresults();


		/** A static instance of our Resource Management / Communication
		class.  It's a member of this class because that way derived
		classes can use it easily; it's static because there should
		only be one instance EVER.  The instance of RMC in the MWTask
		class is actually a pointer to this one... */
		static MWRMComm * RMC;


	protected:
	
		/**@name A. Pure Virtual Methods
		These are the methods from the MWDriver class that a user 
		{\bf must} reimplement in order to have to create an application. 
		*/
	
		//@{
		/** This function is called to read in all information
		specific to a user's application and do any initialization on
		this information.
		*/
		virtual MWReturn get_userinfo( int argc, char *argv[] )=0;
  
		/** This function must return a number n > 0 of pointers
		to Tasks to "jump start" the application.
			
		The MWTasks pointed to should be of the task type derived
		for your application
		*/
		virtual MWReturn setup_initial_tasks( int *n, MWTask ***task ) = 0;
	
		/** 
		This function performs actions that happen
		once the Driver receives notification of a completed task.  
		You will need to cast the MWTask * to a pointer of the Task type 
		derived for your application.  For example

				\begin{verbatim}
				My_Task *dt = dynamic_cast<My_Task *> ( t );
				assert( dt );     
				\end{verbatim}
		*/
		virtual MWReturn act_on_completed_task( MWTask * ) = 0;

		/**
			This is called once the Worker is done with the TaskContainer.
			- Put the logic for this in MWDriver
		*/
		//virtual MWReturn act_on_completed_task( MWTaskContainer * ) = 0;

		/**
		This function should be implemented by the application
		to assign the workClass number to the worker if it is doing
		intelligent work scheduling.
		*/
		virtual MWReturn act_on_starting_worker ( MWWorkerID *w );


		/**
		A common theme of Master-Worker applications is that there is 
		a base amount of "initial" data defining the problem, and then 
		just incremental data defining "Tasks" to be done by the Workers.

		This one packs all the user's initial data.  It is unpacked 
		int the worker class, in unpack_init_data().
		*/
		virtual MWReturn pack_worker_init_data( void ) = 0;
	
		/**
     This one unpacks the "initial" information sent to the driver
     once the worker initializes. 
     
     Potential "initial" information that might be useful is...
     \begin{itemize}
     \item Information on the worker characteristics  etc...
     \item Information on the bandwith between MWDriver and worker
     \end{itemize}
		   
     These sorts of things could be useful in building some 
     scheduling intelligence into the driver.
     
  */
  virtual void unpack_worker_initinfo( MWWorkerID *w ) {};
	
  /**
     OK, This one is not pure virtual either, but if you have some 
     "driver" data that is conceptually part of the task and you wish
     not to replicate the data in each task, you can pack it in a
     message buffer by implementing this function.  If you do this, 
     you must implement a matching unpack_worker_task_data()
     function.
  */
	
  virtual void pack_driver_task_data( void ) {};
	
        /**@name Data streaming/subtasks
        These are the methods that must be implemented to enable data streaming/
		subtasks.
		- act_on_completed_subtask()
		- MWWorker::execute_subtask
		- MWWorker::set_num_subtask
		- MWTask::pack_subresults
		- MWTask::unpack_subresults

		A subtask is identified by a MWTask and a subtask id. The subtask id 
		starts at 0 and ends at n-1 subtasks.

		When subtasks are enabled, the Master sends the task to the Worker as 
	 	usual. But instead of calling execute_task once, execute_subtask is
		called once per subtask id (starting at 0). At the end of each call, the
		result will be returned to the Master. After all subtasks have completed,
		the Master is notified that the Task has completed.	
        */

		/**
   		This function performs an action that happens once the 
		Driver receives notification of a completed subtask. 
   		*/
        virtual MWReturn act_on_completed_subtask(MWTask *) { return ABORT; }
  //@}
	
	
  /**@name B. Task List Management
		   
     These functions are to manage the list of Tasks.  MW provides
     default useful functionality for managing the list of tasks.
		   
  */
	
  //@{

  /// Set number of workclasses
  void workClasses_set ( int num );

  /// Get the number of workclasses
  int workClasses_get ( );

  /// Stop the worker working on this task container and send it a new task container. This should only be called in act_on_completed_task(MWTask*)
  int stop_work();

  /// get number of workers in the specified work class
  int workClasses_getworkers ( int num );
  /// get number of tasks in the specified work class
  int workClasses_gettasks ( int num );
  /// Pack new init data for all workers in group i
  int refreshWorkers ( int i, MWREFRESH_TYPE );

	
  /// Add a task to the list
  void addTask( MWTask * );
	
  /** Add a bunch of tasks to the list.  You do this by making
      an array of pointers to MWTasks and giving that array to 
      this function.  The MWDriver will take over memory 
      management for the MWTasks, but not for the array of 
      pointers, so don't forget to delete [] it! */
  void addTasks( int, MWTask ** );

  /** This will add a list of tasks that are sorted by key.  Efficiency can 
      be greatly improved by using this function */  
  void addSortedTasks( int n, MWTask **add_tasks );

private:
  /// This is a helper function for addSortedTasks().
  void addTaskByKey( MWTask *add_task );

public:
  /// Sets the function that MWDriver users to get the "key" for a task
  void set_task_key_function( MWKey (*)( MWTask * ) );
	
  /// Set the mode you wish for task addition.
  int set_task_add_mode( MWTaskAdditionMode );
	
  /// Set the mode you wish for task retrieval.  
  int set_task_retrieve_mode( MWTaskRetrievalMode );
	
  /// This sorts the task list by the key that is set
  int sort_task_list( void );
	
  /** This deletes all tasks in the task list with a key worse than 
      the one specified */
  int delete_tasks_worse_than( MWKey );
	
  /// returns the number of tasks on the todo list.
  int get_number_tasks();

  int getNumWorkers ( );
	
  /// returns the number of running tasks.
  int get_number_running_tasks();

  /// (Mostly for debugging) -- Prints the task keys in the todo list
  int print_task_keys( void );


  //@}
  void reassignSuspendedTask ( MWList<MWTask> * running, MWWorkerID * w );
  void reassignIdleTask ( MWList<MWWorkerID> * workers, MWWorkerID * w );
  int matchTask ( void *arg1, void *arg2 );
  MWWorkerID* numberworker ( void *arg1, void *arg2 );

	
  /** @name C. Worker Policy Management */
  //@{
  /** Set the policy to use when suspending.  Currently 
      this can be either DEFAULT or REASSIGN */
  void set_suspension_policy( MWSuspensionPolicy );

  ///   Sets the machine ordering policy.
  int set_machine_ordering_policy( MWMachineOrderingPolicy );

  //@}

  /** @name.  task timeout policy.  

      MW provides a mechanism for performing tasks on workers that are 
      potentially "lost".  If the RMComm fails to notify MW of a worker
      going away in a timely fashion, the state of the computing platform
      and MW's vision of its state may become out of synch.  In order to
      make sure that all tasks are done in a timely fashion, the user may set
      a time limit after which a task running on a "lost" worker 
      may be rescheduled.
   */

  //@{

private:	

  /** If false : workers never timeout and can potentially work forever on a task
      If true : workers time out after worker_timeout_limit seconds */	
  bool worker_timeout;

  /** Limit of seconds after which workers are considered time out and 
      tasks are re-assigned */
  double worker_timeout_limit;

  /** frequency at which we check if there are timed out workers */
  int worker_timeout_check_frequency;

  /** based on the time out frequency, next timeout check time*/
  int next_worker_timeout_check;

  /** Go through the list of timed out WORKING workers and reschedule tasks */
  void reassign_tasks_timedout_workers();

public:

  /** Sets the timeout_limit and turn worker_timeout to 1 */
  void set_worker_timeout_limit(double timeout_limit, int timeout_frequency);

  //@}


  /**@name D. Event Handling Methods
		   
     In the case that the user wants to take specific actions
     when notified of processors going away, these methods
     may be reimplemented.  Care must be taken when
     reimplementing these, or else things may get messed up.
     
     Probably a better solution in the long run is to provide 
     users hooks into these functions or something. 
		   
     Basic default functionality that updates the known
     status of our virtual machine is provided. 
		   
  */
	
  //@{
	
  /** Here, we get back the benchmarking
      results, which tell us something about the worker we've got.
      Also, we could get some sort of error back from the worker
      at this stage, in which case we remove it. */

  virtual MWReturn handle_benchmark( MWWorkerID *w );

  /** This is what gets called when a host goes away.  We figure out
      who died, remove that worker from our records, remove its task
      from the running queue (if it was running one) and put that
      task back on the todo list. */
	
  virtual void handle_hostdel();
	
  /** Implements a suspension policy.  Currently either DEFAULT or
      REASSIGN, depending on how suspensionPolicy is set. */
  virtual void handle_hostsuspend();
	
  /** Here's where you go when a host gets resumed.  Usually, 
      you do nothing...but it's nice to know...*/
  virtual void handle_hostresume();
  
  /** We do basically the same thing as handle_hostdel().  One might 
      {\em think} that we could restart something on that host; 
      in practice, however -- especially with the Condor-PVM RMComm
      implementation -- it means that the host has gone down, too.
      We put that host's task back on the todo list.
  */
  virtual void handle_taskexit();

  /** Routine to handle when the communication layer says that a
      checksum error happened. If the underlying Communitor
      gives a reliably reliable communication then this messge
      need not be generated. But for some Communicators like
      MW-File we may need some thing like this.
  */
  virtual void handle_checksum ();

  //@}
	
  /** @name E. Checkpoint Handling Functions
			
      These are logical checkpoint handling functions.  They are
      virtual, and are *entirely* application-specific.  In them, the
      user must save the "state" of the application to permanent
      storage (disk).  To do this, you need to:
			
      \begin{itemize}
      \item Implement the methods write_master_state() and
      read_master_state() in your derived MWDriver app.
      \item Implement the methods write_ckpt_info() and 
      read_ckpt_info() in your derived MWTask class.
      \end{itemize}
			
      Then MWDriver does the rest for you.  When checkpoint() is
      called (see below) it opens up a known filename for writing.
      It passes the file pointer of that file to write_master_state(), 
      which dumps the "state" of the master to that fp.  Here 
      "sate" includes all the variables, info, etc of YOUR
      CLASS THAT WAS DERIVED FROM MWDRIVER.  All state in
      MWDriver.C is taken care of (there's not much).  Next, 
      checkpoint will walk down the running queue and the todo
      queue and call each member's write_ckpt_info().  
      
      Upon restart, MWDriver will detect the presence of the 
      checkpoint file and restart from it.  It calls 
      read_master_state(), which is the inverse of 
      write_master_state().  Then, for each task in the 
      checkpoint file, it creates a new MWTask, calls 
      read_ckpt_info() on it, and adds it to the todo queue.
			
      We start from there and proceed as normal.
			
      One can set the "frequency" that checkpoint files will be 
      written (using set_checkpoint_frequency()).  The default
      frequency is zero - no checkpointing.  When the frequency is
      set to n, every nth time that act_on_completed_task gets 
      called, we checkpoint immediately afterwards.  If your
      application involves "work steps", you probably will want to 
      leave the frequency at zero and call checkpoint yourself
      at the end of a work step.
			
  */
  //@{
	
  /** This function writes the current state of the job to disk.  
      See the section header to see how it does this.
      @see MWTask
  */
  void checkpoint();
	
  
  /** This function does the inverse of checkpoint.  
      It opens the checkpoint file, calls read_master_state(), 
      then, for each task class in the file, creates a MWTask, 
      calls read_ckpt_info on it, and adds that class to the
      todo list. */
  void restart_from_ckpt();
	
	
  /** This function sets the frequency with with checkpoints are
      done.  It returns the former frequency value.  The default
      frequency is zero (no checkpoints).  If the frequency is n, 
      then a checkpoint will occur after the nth call to 
      act_on_completed_task().  A good place to set this is in
      get_userinfo().
      @param freq The frequency to set checkpoints to.
      @return The former frequency value.
  */

  int set_checkpoint_frequency( int freq );
	
  /** Set a time-based frequency for checkpoints.  The time units
      are in seconds.  A value of 0 "turns off" time-based 
      checkpointing.  Time-based checkpointing cannot be "turned 
      on" unless the checkpoint_frequency is set to 0.  A good
      place to do this is in get_userinfo().
      @param secs Checkpoint every "secs" seconds
      @return The former time frequency value.
  */
  int set_checkpoint_time( int secs );
	
  /** Here you write out all 'state' of the driver to fp.
      @param fp A file pointer that has been opened for writing. */
  virtual void write_master_state( FILE *fp ) {};
	
  /** Here, you read in the 'state' of the driver from fp.  Note
      that this is the reverse of write_master_state().
      @param fp A file pointer that has been opened for reading. */
  virtual void read_master_state( FILE *fp ) {};

  
  /** Swap the unused part of the (sorted and indexed) TODO list 
      onto disk (file name: TODO_tasks.id), so that the TODO tasks
      with smaller key values are kept in memory, and others are kept
      in a file. When there are too many tasks in memory, we swap most
      of them out; when there are too few tasks in memory, we swap in
      more. The MW application will be done only when all TODO tasks 
      are done. When checkpointing, it only saves TODO tasks in memory, 
      and the pointers (the file name) to the swap file. 
      
      To help swapping in tasks, we will record the number of tasks
      swapped in at the beginning of the file, and skip then when read 
      to swap file. The num_to_skip is initially set to 0.
      
      The swap file for TODO tasks is fully-sorted, and looks like:
      
        num_to_skip
	trash_task1
	trash_task2
	...
	...
	task1
	task2
	...
	...
  	task[num_tasks]
   **/
   
  /** When swapping out, we walk through the sorted TODO list, keep the first
      num_in_mem tasks in memory (copy into a newly created list). For the 
      rest of the list, merge it in order with existing tasks in the file 
      TODO_tasks, and write them into a newly created swap file. Tasks
      in swap file whose key values are larger than max_key are removed. When 
      done, delete the old TODO list and swap file, switch to the new ones.
      
      If the old TODO list is not too large (<max_in_mem), won't swap out.
      Return value indicates whether the swapping is successful or not. 
   **/

int MIN_IN_MEM;	   /* == 256            */
int NUM_IN_MEM_LO; /* == MIN_IN_MEM*8   */
int NUM_IN_MEM_HI; /* == MIN_IN_MEM*100 */
int MAX_IN_MEM;    /* == MIN_IN_MEM*800 */
  
public:
  bool swap_out_todo_tasks(int num_in_mem = 2048, 
		  int max_in_mem = 204800, 
		  double max_key = DBL_MAX);

  /** When swapping in, we first read first num_in_mem (or all tasks if not so
      many tasks in swap file) tasks from the swap file, create tasks and 
      append to a newly created list. Then we sorted insert the existing tasks 
      in memory into the new list. When done, switch the list and file.

      If the old TODO list is not too small (>num_in_mem), won't swap in.
      Return value indicates whether the swapping is successful or not. 
   **/
  bool swap_in_todo_tasks(int min_in_mem = 256, 
		  int num_in_mem = 25600);

  /** See whether we still have TODO tasks */
  bool is_TODO_empty();

  /* Read and write tasks */
  MWTask * read_mem_task(MWList<MWTask> *tasks);
  MWTask * read_file_task(FILE *f);
  void write_task(FILE *f, MWTask *t);

private:
  bool has_task_swapped;  /* inited as false */
  double max_task_key; 	  /* inited as -DBL_MAX, 
			     tasks with larger key values should be removed */
  
public:
  /** It's really annoying that the user has to do this, but
      they do.  The thing is, we have to make a new task of the
      user's derived type when we read in the checkpoint file.
      
      If your application coredumps when trying to restart from
      a checkpoint, it might be becasue you haven't implemented this function.

      \begin{verbatim}
      MWTask* gimme_a_task() {
      return new <your derived task class>;
      }
      \end{verbatim}
  */
  virtual MWTask *gimme_a_task() = 0;
  //@}

  /**   @name Benchmarking
	We now have a user-defined benchmarking phase.  
	The user can "register" a task that is sent to each worker upon startup.
	This way, the user knows which machines are fastest, and MW can perform
	can automatic "normalization" of the equivalent CPU time.
  */
  //@{

  /** register the task that will be used for benchmarking. */
  void register_benchmark_task( MWTask *t ) { bench_task = t; };

  /** get the benchmark task */
  MWTask * get_benchmark_task() { return bench_task; };

  //@}

private:
	
  /** @name Main Internal Handling Routines */
  //@{
  
  /** This method is called before master_mainloop() is.  It does 
      some setup, including calling the get_userinfo() and 
      create_initial_tasks() methods.  It then figures out how 
      many machines it has and starts worker processes on them.
      @param argc The argc from the command line
      @param argv The argv from the command line
      @return This is the from the user's get_userinfo() routine.
      If get_userinfo() returns OK, then the return value is from 
      the user's setup_initial_tasks() function.
  */

  MWReturn master_setup( int argc, char *argv[] );  
	
  /** This is the main controlling routine of the master.  It sits
      in a loop that accepts a message and then (in a big switch 
      statement) calls routines to deal with that message.  This loop
      ends when there are no jobs on either the running or todo queues.
      It is probably best to see the switch staement yourself to see
      which routines are called when. */

  MWReturn master_mainloop();      
	
  /**
     unpacks the initial worker information, and sends the
     application startup information (by calling pure virtual
     {\tt pack_worker_init_data()}
     
     The return value is taken as the return value from the user's 
     {\tt pack_worker_init_data()} function.
     
  */

  MWReturn worker_init( MWWorkerID *w );
	
  /** 
      This routine sets up the list of initial tasks to do on the 
      todo list.  In calls the pure virtual function 
      {\tt setup_initial_tasks()}.
      @return Is taken from the return value of
      {\tt setup_initial_tasks()}.
  */
  MWReturn create_initial_tasks();  
	
  /**
     Act on a "completed task" message from a worker.
     Calls pure virtual function {\tt act_on_completed_task()}.
     @return Is from the return value of {\tt act_on_completed_task()}.
  */
  MWReturn handle_worker_results( MWWorkerID *w );

	
  /** We grab the next task off the todo list, make and send a work
      message, and send it to a worker.  That worker is marked as 
      "working" and has its runningtask pointer set to that task.  The
      worker pointer in the task is set to that worker.  The task
      is then placed on the running queue. */
  void send_task_to_worker( MWWorkerID *w );

  //by jae
  void lockSendTaskToWorker( MWWorkerID *w ); 
	
  /**  After each result message is processed, we try to match up
       tasks with workers.  (New tasks might have been added to the list
       during processing of a message).  Don't send a task to
       "nosend", since he just reported in.
  */
  void rematch_tasks_to_workers( MWWorkerID *nosend );
	
  /** A wrapper around the lower level's hostaddlogic.  Handles
      things like counting machines and deleting surplus */
  void call_hostaddlogic();

  /// Kill all the workers
  void kill_workers();
	
  /** This is called in both handle_hostdelete and handle_taskexit.
      It removes the host from our records and cleans up 
      relevent pointers with the task it's running. */
 
 void hostPostmortem ( MWWorkerID *w );



  /** The control panel that controls the execution of the
      independent mode. */
  void ControlPanel ( );	

  MWReturn master_mainloop_oneshot ( int buf_id, int sending_host);

  int outstanding_spawn;
  
  //@}
	
  /**@name Internal Task List Routines 
		   
     These methods and data are responsible for managing the 
     list of tasks to be done
  */
	
  //@{
  
  /// This puts a (generally failed) task at the beginning of the list
  void pushTask( MWTask * );
	
  /// Get a Task.  
  MWTask *getNextTask( MWGroup *grp );
	
  /// This puts a task at the end of the running list
  void putOnRunQ( MWTask *t );
	
  /// Removes a task from the queue of things to run
  MWTask * rmFromRunQ ( int jobnum );
	
  /// Print the tasks in the list of tasks to do
  void printRunQ();
	
  /** Add one task to the todo list; do NOT set the 'number' of
      the task - useful in restarting from a checkpoint */
  void ckpt_addTask( MWTask * );
	
  /// returns the worker this task is assigned to, NULL if none.
  MWWorkerID * task_assigned( MWTask *t );
	
  /// Returns true if "t" is still in the todo list
  bool task_in_todo_list( MWTask *t );
	
  /** A pointer to a (user written) function that takes an MWTask
      and returns the "key" for this task.  The user is allowed to 
      change the "key" by simply changing the function
  */
protected:
  MWKey (*task_key)( MWTask * );

	  void give_work ( MWWorkerID *w, MWTask *t );

private:
	
  /// Where should tasks be added to the list?
  MWTaskAdditionMode addmode;
	
  /// Where should tasks by retrived from the list
  MWTaskRetrievalMode getmode;

  /** A pointer to the function that returns the "key" by which machines are ranked.
      Right now, we offer only some (hopefully useful) default functions that are 
      set through the machine_ordering_policy
   */

  MWKey (*worker_key)( MWWorkerID * );
	

  MWMachineOrderingPolicy machine_ordering_policy;

  /** MWDriver keeps a unique identifier for each task -- 
      here's the counter */
  int task_counter;
	
  /// Is the list sorted by the current key function
  bool listsorted;
	
  /// The head of the list of tasks to do
  MWList<MWTask> * todo;

public:
  /// This is Jeff's nasty addition so that he can get access
  ///  to the tasks on the master
  MWTask *get_todo_head();

  //by jae
  void increasePendingNum(void);
  void decreasePendingNum(void);
  int getPendingNum(void) { return m_pending_num; }

  void getGlobalLock(void);
  void getGlobalUnLock(void);

private:
	//by jae
  int m_pending_num;
  pthread_mutex_t m_pending_mutex;
  pthread_cond_t m_pending_cond;
	
  /// The head of the list of tasks that are actually running
  MWList<MWTask> * running;

  /// The reassigned task list
  MWList<int> * reassigned_tasks;
  /// The reassigned and done task list
  MWList<int> * reassigned_tasks_done;
  /// All running tasks are reassigned_done
  bool normal_tasks_all_done;

  /// Insert a task with an id of tid into a int* list
  inline void insert_into_set(MWList<int> * list, int tid) {
	  int i, num = list->number();
	  int *t;
	  bool found = false;
	  
	  for (i=0; i<num; i++) {
		t = (int*)list->Remove();
		if (*t==tid)
			found = true;
		list->Append(t);
	  }
	  if (!found) {
	  	  t = new int;
		  *t = tid;
		  list->Append(t);
	  }
  }
  /// See if a task is in a given int* list
  inline bool in_set(MWList<int> * list, int tid) {
	  int i, num = list->number();
	  int *t;
	  
	  for (i=0; i<num; i++) {
		t = (int*)list->Remove();
		list->Append(t);
		if (*t==tid) 
			return true;
	  }
	  return false;
  } 
  //@}
  
  /**@name Worker management methods
		   
     These methods act on the list of workers
     (or specifically) ID's of workers, that the
     driver knows about.
     
  */
	
  //@{
	
  /// Adds a worker to the list of avaiable workers
  void addWorker ( MWWorkerID *w );
  
  /// Looks up information about a worker given its task ID
  MWWorkerID * lookupWorker( int tid );
	
  /// Removes a worker from the list of available workers
  MWWorkerID * rmWorker ( int tid );

  /// This function removes worker from the list, removes it and deletes
  /// the structure.
  void worker_last_rites ( MWWorkerID *w );

public:	
  /// Prints the available workers
  void printWorkers();

  /** Another terrible addition so that Jeff can print out the worker list
      in his own format */
  MWWorkerID *get_workers_head();

private:

  /// The head of the list of workers.
  MWList<MWWorkerID> * workers;

  /// Here's where we store what should happen on a suspension...
  MWSuspensionPolicy suspensionPolicy;

  /// Based on the ordering policy, place w in the worker list appropriately
  void sort_worker_list();
public:
  /// Counts the existing workers
  int numWorkers();
	
  /// Counts the number of workers in the given arch class
  int numWorkers( int arch );

  /// Counts the number of workers in the given state
  int numWorkersInState( int ThisState );

  /// Returns the value (only) of the best key in the Todo list 
  MWKey return_best_todo_keyval( void );

  /// Returns the best value (only) of the best key in the Running list.
  MWKey return_best_running_keyval( void );  

private:
	
  //@}
	
  /** @name Checkpoint internal helpers... */
  //@{

  /** How often to checkpoint?  Task frequency based. */
  int checkpoint_frequency;
	
  /** How often to checkpoint?  Time based. */
  int checkpoint_time_freq;
  /** Time to do next checkpoint...valid when using time-based 
      checkpointing. */

  long next_ckpt_time;
	
  /** The number of tasks acted upon up to now.  Used with 
      checkpoint_frequency */
  //int num_completed_tasks;
	
  /** The benchmark task. */
  MWTask *bench_task;
	
  /** The name of the checkpoint file */
  const char *ckpt_filename; // JGS
  
  /** When submitting multiple jobs from the same directory, the checkpoint files all the 
   *  the same name if we use the fixed ckpt_filename. So we made some changes:
   *  (1) Allow the user set their checkpoint file name by calling 
   *  		set_ckpt_filename(const char*) before finishing get_user_info()
   *  (2) If the user doesn't specify the ckpt_filename, then use
   *  		$(master_exec_name)-ckpt.$(Process)
   *  (3) When reading checkpoints during execution, MW will just read the file. When the
   *      application restarts, it will find the file with the same name as user specified,
   *      or as the default name (choose the first one when there are multiple files). 
   */
  // char *ckpt_filename;
  
	
  //@}

  /** The instance of the stats class that takes workers and later
      prints out relevant stats... */
  // Changed into MWStatistics * to make insure happy (then I can new and delete it).
  MWStatistics* stats;

  /// This returns the sum of the bench results for the currently working machines

  double get_instant_pool_perf();

  /// This returns the hostname
  char *getHostName ();
  


#if defined( XML_OUTPUT )
  
  /**@name XML and Status Methods.

     This function is called by the CORBA layer to get the
     XML status of the MWDriver.
  */

  //@{

public:
  ///
  void  write_XML_status();


  /**
     If you want to display information about status of some
     results variables of your solver, you have to dump a string in
     ASCII, HTML or XML format out of the following method. 
     The iMW interface will be in charge of displaying this information
     on the user's browser.
   */

  virtual char* get_XML_results_status(void );

private:

  ///
  char* get_XML_status();

  ///
  char* get_XML_job_information();

  ///
  char* get_XML_problem_description();

  ///
  char* get_XML_interface_remote_files();
  
  ///
  char* get_XML_resources_status();

  ///
  const char *xml_filename; // JGS

  ///
  const char *xml_menus_filename;

  ///
  const char *xml_jobinfo_filename;

  ///
  const char *xml_pbdescrib_filename;

  /// Set the current machine information
  void get_machine_info();

  /// Returns a pointer to the machine's Arch
  char *get_Arch(){ return Arch; }
  
  /// Returns a pointer to the machine's OpSys
  char *get_OpSys(){ return OpSys; }

  /// Returns a pointer to the machine's IPAddress
  char *get_IPAddress(){ return IPAddress; }

  ///
  double get_CondorLoadAvg(){ return CondorLoadAvg; }

  ///
  double get_LoadAvg(){ return LoadAvg; }

  ///
  int get_Memory(){return Memory;}

  ///
  int get_Cpus(){return Cpus;}

  ///
  int get_VirtualMemory(){return VirtualMemory;}

  ///
  int get_Disk(){return Disk;}

  ///
  int get_KFlops(){return KFlops;}

  ///
  int get_Mips(){return Mips;}

  ///
  char Arch[64];

  ///
  char OpSys[64];

  ///
  char IPAddress[64];

  ///
  double CondorLoadAvg;

  ///
  double LoadAvg;

  ///
  int Memory;

  ///
  int Cpus;

  ///
  int VirtualMemory;

  ///
  int Disk;

  ///
  int KFlops;

  ///
  int Mips;

  /// Utility functions used by get_machine info
  int check_for_int_val(char* name, char* key, char* value);
  double check_for_float_val(char* name, char* key, char* value);
  int check_for_string_val(char* name, char* key, char* value);

#endif

  /// The name of the machine the worker is running on.
  char mach_name[64];

  /// for measuring network connectivity
  double defaultTimeInterval;
  int defaultIterations;
  int defaultPortNo;

#ifdef MEASURE
  // Gather performance measurement information at runtime
public:
  // read/write measurement files
  void _measure_read_options(const char* opt_fname);
  void _measure_dump_header();
  void _measure_dump_records();
private:
  // Do the measurement
  void _measure();
  // Reset the numbers
  void _measure_reset();

private:
  // Names for the measurement option and record files
  // We want to make the record file unique across different runs, so better
  //  add the ProcID part. But how to get the $(Cluster).$(Process) part? 
  //  Here is Peter's method: use _ProcID_$(Cluster).$(Process) in the 
  //  argument line in the master's submit file, and Condor will substitute 
  //  them into the real ProcID before it executes the master executable, 
  //  so that I can read this from argv. 
  // Suggestion from Peter: don't make this with fixed argc position.
  char _measure_opt_file_name[_POSIX_PATH_MAX]; 	// = "_measure_opt";
  char _measure_rec_fname_prefix[_POSIX_PATH_MAX]; 	//  = "_measure_rec";
  char _measure_ProcID_prefix[_POSIX_PATH_MAX]; 	//  = "_ProcID.";
  char _measure_rec_file_name[_POSIX_PATH_MAX];
  
  // Time intervals to read options and dump measurement record (append)
  int  _measure_read_opt_interval; 	//  = 300;
  int  _measure_dump_rec_interval; 	//  = 120;
  time_t  _measure_last_read_opt_time;
  time_t  _measure_last_dump_rec_time;
  
  // Recent queue statistics
  // int  _measure_num_worker;		// = workers->number();
  // int  _measure_num_task_TODO;	// = todo->number();
  int  _measure_num_task_Done;
  
  // Recent master performance
  double _measure_master_recv_time;	// time used for mainloop's recv()
  double _measure_master_recv_cpu_time;	// CPU time used for mainloop's recv()
  double _measure_master_wall_time;	// time elapsed
  double _measure_master_cpu_time;	// CPU time (both user and system)
  double _measure_master_act_on_completed_task_time;	// Time in user code
  double _measure_master_act_on_completed_task_cpu_time;// CPU time in user code

  // Performance of recent finished tasks
  double _measure_task_RoundTrip_time;	// time from sending task to deleting it
  double _measure_task_MP_master_time;	// time used when master packing/sending 
  double _measure_task_MP_worker_time;	// time used when worker unpacking/recving
  double _measure_task_MP_worker_cpu_time; // CPU time used when worker unpacking/recving
  double _measure_task_Exe_wall_time;	// time elapsed when worker executing task
  double _measure_task_Exe_cpu_time;	// CPU time when worker executing task

  // Below are adaptation logic 
  int _measure_use_adaptation;
  // To keep track of the master and worker utilization
  double _measure_master_wait_rate;
  int _measure_master_busy_times;
  int _measure_master_wait_times;
  MWWorkerID *_measure_current_worker;
  int _measure_target_num_workers;
  
  // Adapt the number of workers when necessary: 
  //  assumes that it's called in the mainloop (can access w that just returns a RESULT)
  //  also assumes that it can access the most recent measurement numbers
  void _measure_adapt(); 
  void _measure_remove_worker();
#endif

public:
	/* Send an asynchronous UPDATE_FROM_DRIVER message to the worker
		-2 to send to all workers
		-1 to send to current worker
		>0 to send to a specific worker
		User should override pack_update if some data needs to be sent with update.
	*/
	void prepare_update(int mode);

	/* user packs update data here. Used in conjunction with prepare_update */
	virtual void pack_update(){};	

	/* user tells Driver to send a DRIVER_UPDATE message to the current worker 
	 * from who a task result or async message has been send */
	int send_update_message();

	/* user tells Driver to send a DRIVER_UPDATE message to a specific worker */
	int send_update_message_to(int worker);
};

///  This is the function allowing the workers to be sorted by KFlops
MWKey kflops( MWWorkerID * );

/// This is the function allowing the workers to be sorted by benchmarkresult
MWKey benchmark_result( MWWorkerID * );

#endif


