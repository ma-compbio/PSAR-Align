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
/* MWWorkerID.C
*/
#include "MW.h"
#include "MWWorkerID.h"
#include "MWTask.h"
#include "MWSystem.h"

#include <memory.h>
#include <string.h>
#include <time.h>
#include <errno.h>
#include <stdio.h>
#include <stdarg.h>



int *MWWorkerID::vids = new int[MW_MAX_WORKERS];  
int MWWorkerID::vids_inited = 0;

char *MWWorkerID::ulv_filename = "/tmp/mw-ulv";
FILE *MWWorkerID::lfp = NULL;
extern int *MWworkClassWorkers;
extern int MWworkClasses;

/** This static array is used to hold strings that represent the
    MWworker_states.  They are handy for printing out a worker - now
    you'll see "WORKING" instead of "2" in printself().
    @see MWworker_states */

const char *MWworker_statenames[] = {
    "INITIALIZING", 
    "IDLE", 
    "BENCHMARKING/INITIAL DATA HANDLING",
    "WORKING", 
    "SUSPENDED",
    "EXITED"
};


void 
MWWorkerID::init_vids() {
	if (vids_inited == 0) {
		for ( int i = 0; i < MW_MAX_WORKERS; i++ )
			vids[i] = 0;
		vids_inited = 1;
	}
}

MWWorkerID::MWWorkerID() 
{
    init_vids();
    memset ( mach_name, 0, 64 ); 

    id1 = -1;
    id2 = -1;
	virtual_id = -1;
	doomed = FALSE;
	// We use a negative number to signify that no benchmark has been set
	bench_result = -1.0;

    runningtask = NULL;
	/*
    next = NULL;
	*/

    state = INITIALIZING;
    isBenchMark = FALSE;
    isBenchMarkAvailable = FALSE;

	total_suspended = 0.0;
	total_working = 0.0;
	start_time = 0.0;
	total_time = 0.0;
	last_event = 0.0;
	cpu_while_working = 0.0;


	normalized_cpu_working_time = 0.0;
	sum_benchmark = 0.0;
	num_benchmark = 0;

	
    CondorLoadAvg = NO_VAL;
    LoadAvg = NO_VAL;
    Memory = NO_VAL;
    Cpus = NO_VAL;
    VirtualMemory = NO_VAL;
    Disk = NO_VAL;
    KFlops = NO_VAL;
    Mips = NO_VAL;

    memset(OpSys, 0, sizeof(OpSys)); 
    memset(Arch, 0 , sizeof(Arch)); 

	set_vid ( get_next_vid() );
    double currentTime = MWSystem::gettimeofday();
    networkLatency = 0.0;
    networkLatency_lastMeasuredTime = currentTime;

	group = NULL;
	executable_name = NULL;
}

MWWorkerID::~MWWorkerID()
{
	delete group;
}


enum MWworker_states MWWorkerID::currentState()
{
	return state;
}

void MWWorkerID::printself( int level ) 
{

    MWprintf ( level, "id1: %d id2: %d, vid: %d, name: \"%s\"\n",
			   id1, id2, virtual_id, mach_name ? mach_name : "" );
    MWprintf ( level, "Memory: %d, KFlops: %d, Mips: %d\n", Memory, KFlops, Mips );
    MWprintf ( level, "arch: \"%d\" state: %s\n", 
			   arch, MWworker_statenames[state] );
	
    if ( runningtask != NULL ) {
        MWprintf ( level, " -> Running task: " );
        runningtask->printself( level );
    }
    else 
        MWprintf ( level, " -> There's no task on this worker.\n" );
    
    /* Print Here other Worker specific Info */    
    
}    

void MWWorkerID::set_machine_name( char *name ) {
    strncpy( mach_name, name, 64 );
}

char * MWWorkerID::machine_name() {
    return mach_name;
}

void MWWorkerID::started() {
        // this is called when a worker first starts up.
	if( state != INITIALIZING && state != BENCHMARKING ) {
		MWprintf( 10, "Danger Will Robinson!  started() worker whose "
				  "state != INITIALIZING && state != BENCHMARKING\n" );
		printself(10);
    }
    
    // The following line is commented out because the state should be INITIALIZING !? - jichuan
    // state = IDLE;

    //XXX Jeff's hack.  To get parallel efficiency "right" on
    //  MWfiles (with worker timeouts), we don't want to count
    //  the machines that just sit there in INITIALIZING state
    //  until they time out.  Thus, we do not start the counter here, but
    //  only when we send a benchmark.

    lprintf ( virtual_id, start_time, "Started on host %s\n", mach_name );
}


void MWWorkerID::benchmark() 
{
	// this is called when a worker is sent the benchmark task
	if ( state != INITIALIZING )
	{
		MWprintf( 10, "Danger Will Robinson! benchmarking a worker "
				"state != INITIALIZING\n" );
		printself(10);
	}
	state = BENCHMARKING;
	isBenchMark = TRUE;

}

void MWWorkerID::benchmarkOver()
{
	if ( state != BENCHMARKING )
	{
		MWprintf( 10, "Danger Will Robinson! benchmarking over of "
				" a worker whose state != BENCHMARKING\n" );
		printself(10);
	}
	isBenchMark = FALSE;

    //XXX Jeff's hack.  To get parallel efficiency "right" on
    //  MWfiles (with worker timeouts), we don't want to count
    //  the machines that just sit there in INITIALIZING state
    //  until they time out.  Thus, we start the counter here

	start_time = MWSystem::gettimeofday();
	state = IDLE;

}

void MWWorkerID::gottask( int tasknum ) {
        // called when we begin a task.
	state = WORKING;
	last_event = MWSystem::gettimeofday();

	lprintf ( virtual_id, last_event, "Assigned task %d\n", tasknum );
}

void MWWorkerID::completedtask( double wall_time, double cpu_time ) {
        // called when the master is told that we finished this task.
    	state = IDLE;
		last_event = MWSystem::gettimeofday();

	total_working += wall_time;
	cpu_while_working += cpu_time;

	if (isBenchMarkAvailable){
	  if ( (cpu_time * bench_result) < -1.0e-6 ) {
	        MWprintf ( 10, "PROBLEM !!! cpu_time = %lf ; bench = %lf\n", cpu_time, bench_result);
	  } else {
	  	normalized_cpu_working_time += cpu_time * bench_result;
	  }
	}
	
	lprintf ( virtual_id, last_event, "Completed task. (tot: %12.4f)\n"
		  "\tWall time: %12.4f\n"
		  "\tCpu time:  %12.4f\n",
		  total_working, wall_time, cpu_time );

}

void MWWorkerID::suspended() {
        // Called when this worker has a task and gets suspended.
	last_event = MWSystem::gettimeofday();
    state = SUSPENDED;

	lprintf ( virtual_id, last_event, "Suspend.\n" );
}

void MWWorkerID::resumed() {

	double now = MWSystem::gettimeofday();

  // called when we resume and have a task.
  if ( state != SUSPENDED ) {
    MWprintf ( 10, "Got resume while not suspended!\n" );
  }
  else
    {
	if ( runningtask ) {
		total_suspended += now - last_event;
		state = WORKING;
	} else if ( isBenchMark == TRUE ) {
		state = BENCHMARKING;
	} else {
		state = IDLE;
	}

    }
  last_event = now;

  lprintf ( virtual_id, now, "Resume.\n" );
}

void MWWorkerID::set_bench_result( double bres )
{
  if( bres < 1.0e-8 ) {
    MWprintf( 10, "Benchmark result must be a *positive* number\nSetting to 0\n" );
    bench_result = 0.0;
  }
  else {
    bench_result = bres;
    sum_benchmark += bres;
    num_benchmark++;

    isBenchMarkAvailable = TRUE;
  }
}

void MWWorkerID::ended() {
        // called when the master is aware of the worker's death.
	double now = MWSystem::gettimeofday();

	if ( state == WORKING ) { 
/* We don't want to do this!  If we're here, it means that we're running
   a task, but we're done - therefore whatever we were doing wasn't
   important and we shouldn't count it towards our work.  (It's not
   part of "GoodPut", in other words.) */
//		total_working += now - last_event; 
	} else if ( state == SUSPENDED ) {
/* But what the heck, we'll count suspended */
		total_suspended += now - last_event;
	}
	state = EXITED;


	if( start_time > 1.0 ) {
	  total_time = now - start_time;
	}
	else {
	  total_time = 0.0;
	}

	lprintf ( virtual_id, now, "Finished.\n" );

	release_vid ( get_vid() );

	if ( group )
	{
		for ( int i = 0; i < MWworkClasses; i++ )
			if ( doesBelong ( i ) )
				deleteGroup ( i );

	}
}

void MWWorkerID::ckpt_stats( double *up, 
			     double *working, 
			     double *susp, 
			     double *cpu, 
			     double *norm, 
			     double *s_bench, 
			     int *n_bench ) {

		// print duration, total_working, total_suspended 

	double now = MWSystem::gettimeofday();

	if ( start_time > 1.0 ) { 
		*up = now - start_time;  
	} else {
		*up = 0.0;   // possible if not yet started().
	}

	*working = total_working;
	*susp = total_suspended;
	*cpu = cpu_while_working;
	*norm = normalized_cpu_working_time;

	*s_bench = sum_benchmark;
	*n_bench = num_benchmark;

/* Don't do any of this - only count goodput... 
	if ( state == WORKING ) {
		*work += now - last_event;
		*cpu += now - last_event;  // punt... 
	} else if ( state == SUSPENDED ) {
		*susp += now - last_event; // we can lose some cpu time here.  Oops. 
	}
*/
}



/* get the lowest vid from 'vids' */
int
MWWorkerID::get_next_vid() {
	int i;
	for ( i=0 ; i< MW_MAX_WORKERS ; i++ ) {
		if ( MWWorkerID::vids[i] == 0 ) {
			MWWorkerID::vids[i] = 1;
			return i;
		}
	}
	MWprintf ( 10, "ERROR: Ran out of virtual IDs in MWWorkerID!\n" );
	return -1;
}

/* Set virtual id */
void 
MWWorkerID::set_vid( int i ) { 
	virtual_id = i; 
}

/* Returns a virtual id to the pool */
void 
MWWorkerID::release_vid( int vid ) { 
	MWWorkerID::vids[vid] = 0; 
}	

void
MWWorkerID::lprintf ( int vid, double now, char *fmt, ... ) {

		/* let's return now so we don't inflict this code on others. */
	return;

	if ( !lfp ) {
		if ( (lfp=fopen( ulv_filename, "w" )) == NULL ) {
			MWprintf ( 10, "Failed to open %s in MWWorkerID::lprintf, "
					   "errno %d\n", ulv_filename, errno );
			return;
		} else {
			fprintf ( lfp, "#MW\n" );
		}
	}

	fprintf ( lfp, "(%03d) ", vid );

	struct tm	*tm;
	time_t time = (time_t) now;
	tm = localtime( (time_t *)&time );
	fprintf( lfp, "%d/%d %02d:%02d:%02d.%03d ",
			 tm->tm_mon + 1, tm->tm_mday,
			 tm->tm_hour, tm->tm_min, tm->tm_sec, 
			 ((int)(((double) ((double) now -((int) now ))) * 1000))
	);

	va_list ap;
	va_start( ap, fmt );
	vfprintf( lfp, fmt, ap );
	va_end( ap );

	if ( fmt[strlen(fmt)-1] == '\n' ) {
		fprintf ( lfp, "...\n" );
	} else {
		fprintf ( lfp, "\n...\n" );
	}

	fflush ( lfp );
}

/* The following (until the end of the file) was written by one of 
   Jean-Pierre's students.  It isn't exactly efficient, dumping
   a condor_status output to a file and reading from that.  I understand
   that it does work, however.  -MEY (8-18-99) */

int MWWorkerID::check_for_int_val(char* name, char* key, char* value) {
	if (strcmp(name, key) == 0)
		return atoi(value);
	else return NO_VAL;
}

double MWWorkerID::check_for_float_val(char* name, char* key, char* value) {
	if (strcmp(name, key) == 0)
		return atof(value);
	else return NO_VAL;  
}

int  MWWorkerID::check_for_string_val(char* name, char* key, char* value) {
	if (strcmp(name, key) == 0)
		return 0;
	else return NO_VAL;  
}

void MWWorkerID::get_machine_info() {

	FILE* inputfile;
	char key[200];
	char value[1300];
	char raw_line[1500];
	char* equal_pos;
	int found;
	char temp_str[256];
	
	char zero_string[2];
	
#ifdef INDEPENDENT
		return;
#endif

	memset(zero_string, 0 , sizeof(zero_string));
	
	memset(key, '\0', sizeof(key));
	memset(value, '\0', sizeof(value));

	// This is a terrible kludge, but I'd rather do this in a script.
	// We're going to move this into the RMC class anyway, so it's 
	//  not so bad to fake it for now.

	const bool jeff_hack = false;;

	if( jeff_hack ) {
	// note to jeff_hack: change this script to output to stdout, not filename
	  sprintf( temp_str, "./create_condor_status_file %s", mach_name /*, filename */);
	}
	else {
	  
	  memset(temp_str, '\0', sizeof(temp_str));
	  sprintf(temp_str, "%s/bin/condor_status -l %s", 
		  CONDOR_DIR, mach_name);
	}

	
#ifdef WINDOWS
	return;

#else
	if (( inputfile = popen(temp_str, "r")) == 0) {
	  MWprintf( 10, "Error occurred during attempt to get "
		    "condor_status for %s.\n",  mach_name );
	  return;
	}
	
	while (fgets(raw_line, 1500, inputfile) != 0) {
		found = 0;
		equal_pos = strchr(raw_line, '=');
		
		if (equal_pos != NULL) {
			strncpy(key, raw_line, equal_pos - (raw_line+1));
			strcpy(value, equal_pos+2);
			
			if (CondorLoadAvg == NO_VAL && !found) {
				CondorLoadAvg = check_for_float_val("CondorLoadAvg", 
													key, value);
				if (CondorLoadAvg != NO_VAL)
					found = 1;
			}
			
			if (LoadAvg == NO_VAL && !found) {
				LoadAvg = check_for_float_val("LoadAvg", key, value);
				if (LoadAvg != NO_VAL)
					found = 1;
			}
			
			if (Memory == NO_VAL && !found) {
				Memory = check_for_int_val("Memory", key, value);
				if (Memory != NO_VAL)
					found = 1;
			}
			
			if (Cpus == NO_VAL && !found) {
				Cpus = check_for_int_val("Cpus", key, value);
				if (Cpus != NO_VAL)
					found = 1;
			}
			
			if (VirtualMemory == NO_VAL && !found) {
				VirtualMemory = check_for_int_val("VirtualMemory", key, value);
				if (VirtualMemory != NO_VAL)
					found = 1;
			}
			
			if (Disk == NO_VAL && !found) {
				Disk = check_for_int_val("Disk", key, value);
				if (Disk != NO_VAL)
					found = 1;
			}

			if (KFlops == NO_VAL && !found) {
				KFlops = check_for_int_val("KFlops", key, value);
				if (KFlops != NO_VAL)
					found = 1;
			}

			if (Mips == NO_VAL && !found) {
				Mips = check_for_int_val("Mips", key, value);
				if (Mips != NO_VAL)
					found = 1;
			}

			if ( (strncmp(Arch, zero_string, 1) == 0) && !found) {	     
				if (check_for_string_val("Arch", key, value) == 0){
					strncpy( Arch, value, strlen(value)-1 );
				}
				
				if (strncmp(Arch, zero_string, 1) != 0)
					found = 1;
			}
			
			if ( (strncmp(OpSys, zero_string, 1) == 0) && !found) {	     
				if (check_for_string_val("OpSys", key, value) == 0){
					strncpy( OpSys, value, strlen(value)-1 );
				}
				
				if (strncmp(OpSys, zero_string, 1) != 0)
					found = 1;
			}

			memset(key, '\0', sizeof(key));
			memset(value, '\0', sizeof(value));
			
		}
		
    }

    pclose ( inputfile);
	
	MWprintf(90,"CURRENT MACHINE  : %s \n", mach_name);
	MWprintf(90,"Architecture : %s \n", Arch);
	MWprintf(90,"Operating System : %s \n", OpSys);
	
	MWprintf(90,"CondorLoadAvg : %f\n", CondorLoadAvg);
	MWprintf(90,"LoadAvg : %f\n", LoadAvg);
	MWprintf(90,"Memory : %d\n", Memory);
	MWprintf(90,"Cpus : %d\n", Cpus);
	MWprintf(90,"VirtualMemory : %d\n", VirtualMemory);
	MWprintf(90,"Disk : %d\n", Disk);
	MWprintf(90,"KFlops : %d\n", KFlops);
	MWprintf(90,"Mips : %d\n", Mips);
#endif
}

void
MWWorkerID::setNetworkConnectivity ( double latency )
{
	networkLatency_lastMeasuredTime = MWSystem::gettimeofday();
	networkLatency = latency;
	MWprintf ( 10, "The network latency is %f\n", networkLatency);
}

double
MWWorkerID::getNetworkConnectivity ( double &retVal )
{
	retVal = networkLatency_lastMeasuredTime;
	return networkLatency;
}

void
MWWorkerID::initGroups ( int num )
{
	group = new MWGroup ( num );
}

void
MWWorkerID::addGroup ( int num )
{
	MWworkClassWorkers[num]++;
	group->join ( num );
    return;
}

void
MWWorkerID::deleteGroup ( int num )
{
	MWworkClassWorkers[num]--;
	group->leave ( num );
    return;
}

bool
MWWorkerID::doesBelong ( int num )
{
	return group->belong ( num );
}

MWGroup*
MWWorkerID::getGroup ( )
{
	return group;
}
