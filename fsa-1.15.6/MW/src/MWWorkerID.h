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
#ifndef MWWORKERID_H
#define MWWORKERID_H
#include "MW.h"
#include "MWGroup.h"
#include <string.h>

#define NO_VAL -1

/// A forward reference; needed to define the MWWorkerID class.
class MWTask;

/**
   This is a list of the states in which a worker can be.
 */
enum MWworker_states { 
    /// Some initialization; not ready for work
    INITIALIZING, 
    /// Waiting for the master to send work
    IDLE, 
    /// Benchmarking, or doing application initialization on the initial data
    BENCHMARKING,
    /// Working actively on its task
    WORKING, 
    /// The machine has been suspended
    SUSPENDED,
    /// This worker has finished life.
    EXITED
};

/** This array is used to hold strings that represent the
    MWworker_states.  They are handy for printing out a worker - now
    you'll see "WORKING" instead of "2" in printself().
    @see MWworker_states */
extern const char* MWworker_statenames[];


/**
   This class keeps an identification of the worker application
   useful for the driver class.  It also keeps statistical information
   that is useful to the stats class.  This information will be used at
   the end of a run.
 */
class MWWorkerID  {

        // allow access to worker's data...makes my life easier.
    friend class MWStatistics;

public:

    /// Default constructor
    MWWorkerID();
    
    /// Default destructor
    virtual ~MWWorkerID();
    
		/// Return primary id of this worker
	int get_id1() { return id1; };

		/// Return secondard id of this worker
	int get_id2() { return id2; };

		/** Return the virtual id.  This is an int starting at zero and
			increasing by one with each worker.  When a worker dies, this
			number can be reused by a worker that joins us. */
	int get_vid() { return virtual_id; };

	/// Set primary id
	void set_id1( int i ) { id1 = i; };

	/// Set secondary id
	void set_id2( int i ) { id2 = i; };

    /// Set the machine's name
    void set_machine_name( char * );

    /// Get the current machine information
    void get_machine_info();

    /// Returns a pointer to the machine's name
    char *machine_name();
    
	/// Set this worker's arch class
	void set_arch ( int a ) { arch = a; };

	/// Get this worker's arch class
	int get_arch() { return arch; };

	/// Sets the exec_class
	void set_exec_class ( int num ) { exec_class = num; };

	/// Gets the exec_class
	int get_exec_class ( ) { return exec_class; };

	/// Sets the executable
	void set_executable ( char *exec ) { executable_name = exec; };

	/// returns the executable name
	char* get_executable ( ) { return executable_name; };

		/// Mark this worker so that when the task completes it will exit
	void mark_for_removal() { doomed = TRUE; };

		/// returns TRUE if mark_for_removal was called previously
	int is_doomed() { return doomed; };

		/// 
	void set_bench_result( double bres );

		///
	double get_bench_result() { return bench_result; };

		/// 
	double getNetworkConnectivity ( double &t );
	void setNetworkConnectivity ( double );


    /// The task running on this worker.  NULL if there isn't one.
    MWTask *runningtask;
    
	/*
    /// The next worker ID in the list of available workers.  NULL otherwise
    MWWorkerID *next;
	*/

	/// The state of the worker.
	enum MWworker_states currentState();

	/// Set the state of the worker
	void setState ( MWworker_states st ) { state = st; };
    
		/** Print out a description of this worker.
			@param level The debug level to use */
    virtual void printself( int level = 40 ); 

    /**@name Statistics Collection Calls

       These methods are called when events happen to workers.
    */

    //@{

        /// Returns true if this worker is idle.
    bool isIdle() { return (state==IDLE); }

        /// Returns true if this worker is suspended.
    bool isSusp() { return (state==SUSPENDED); }

    /** This should be called when the master becomes aware of the 
        existence of this worker. */
    void started();
    /// Called when the worker is doing the benchmarking task
    void benchmark();
    /// Called when the worker has done the benchmarking task
    void benchmarkOver();
    /// Called when this worker gets a task.
    void gottask( int tasknum );
    /// Called when this worker just finished the task it was working on.
    void completedtask( double wall_time = 0.0, double cpu_time = 0.0 );
    /// Called when the worker becomes suspended
    void suspended();
    /// Called when the worker resumes.
    void resumed();
    /// Send flowers...this is called when the worker dies.
    void ended();

    //@}

	void initGroups ( int num );
	/// Add worker to a workclass/group
	void addGroup ( int num );
	/// Remove worker from a workclass/group
	void deleteGroup ( int num );
	bool doesBelong ( int num );
	MWGroup *getGroup ( );

    /**@name Checkpointing Call

	   Yes, each instance of a MWWorkerID needs to be able to checkpoint
	   itself.  Why?  It has to store statistics information on itself.
	   It never has to read them in, though, because stats that we write
	   out here wind up being read in by the stats class when we restart.
    */

    //@{

		/// Return the relevant stats info for this worker.

	void ckpt_stats( double *up, 
			 double *working, 
			 double *susp, 
			 double *cpu, 
			 double *norm, 
			 double* s_bench, 
			 int* n_bench );

	//@}
	
		/** @name Worker Information
			This is now collected from condor_status.  It's the work of 
			one of JP's students... 
			
			Note that this Arch is different from arch.  This one is 
			what condor claims its arch is to the outside world. 
		*/
		//@{

  ///
    char Arch[64];

  ///
    char OpSys[64];

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


  ///
    double get_total_time(){ return total_time;};
  ///
    double get_last_event(){return  last_event;};
  ///
    double get_total_suspended(){return total_suspended;};
  ///
    double get_total_working(){return total_working;};

 private:

    int check_for_int_val(char* name, char* key, char* value);
    double check_for_float_val(char* name, char* key, char* value);
    int check_for_string_val(char* name, char* key, char* value);

		//@}

		/** @name Private Data... */
		//@{
	
		/// The machine name of the worker
    char mach_name[64];
	
		/// A "primary" identifier of this machine.  (Was pvmd_tid)
    int id1;
    
		/// A "secondary" identifier of this machine.  (Was user_tid)
    int id2;
    
		/** A "virtual" number for this worker.  Assigned starting at 
			0 and working upwards by one.  Also, when a worker dies, 
			this number can get taken over by another worker */
	int virtual_id;

		/// This worker's arch class
	int arch;   // (should be just a 0, 1, 2, etc...)

	/// This is the exec_class
	int exec_class;

	/// This is the executable_name
	char *executable_name;

		/// TRUE if marked for removal, FALSE otherwise	
	int doomed;  

		/// The results of the benchmarking process.
	double bench_result;

		/// The state of this worker.
    enum MWworker_states state;

		/** @name Time information 
		    Note that all time information is really stored as a double.
		    We use gettimeofday() internally, and convert the result to 
		    a double.  That way, we don't have to mess around with the
		    struct timeval... */
		//@{
		/// The time that this worker started.
	double start_time;
	
		/// The total time that this worker ran for.
	double total_time;

		/// The time of the last 'event'.
	double last_event;
	
		/// The time spent in the suspended state
	double total_suspended;

		/// The time spent working
	double total_working;

		/// The cpu usage while working
	double cpu_while_working;

                /// The benchmarked (weighted) cpu working time
        double benchmarked_cpu_working_time;

               /// The benchmarked (normalized) cpu working time
	double normalized_cpu_working_time;

               /// The sum of the benchmark values for that vid
        double sum_benchmark;

               /// The number of the benchmark values for that vid
        int num_benchmark;
  

		//@}
		//@}

	bool isBenchMark;
	bool isBenchMarkAvailable;

	/** @name Virtual Id Helpers */
	//@{
	/// vids[i] is 1 if virtual id i is taken, 0 if not.
	static int *vids; 

	// So here is the fix to MWStats bugs: I moved the vids[MW_MAX_WORKERS] init logic 
	// into the init_vids() function. 
	static int vids_inited;
	void init_vids();
	
	/// Set virtual id
	void set_vid( int i );
	
	/// Returns the lowest available virtual id; also sets it as used.
	int get_next_vid();
	
	/// Returns a virtual id to the pool
	void release_vid( int vid );
	//@}

	/** @name A Special Printf
		Used to make output for John Bent's ulv program... */
	
	//@{
	/** The log printf - use it to print ONE log event.  It will 
		automagically add the "..." delimiter for you. */
	void lprintf ( int vid, double now, char *fmt, ... );
	/// The name of the ulv log file
	static char *ulv_filename;
	/// The log file fp.
	static FILE *lfp;
	//@}

	/// Variables regarding the nw bandwidth.
	double networkLatency;
	double networkLatency_lastMeasuredTime;
	MWGroup *group;
};


#endif
