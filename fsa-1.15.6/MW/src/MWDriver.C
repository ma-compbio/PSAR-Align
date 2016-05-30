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

#include <cassert>
#include <cerrno>
#include <cstdarg>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
//by jae
#include <pthread.h>
#include <sys/time.h>

#include "MW.h"
#include "MWDriver.h"
#include "MWWorkerID.h"
#include "MWWorker.h"
#include "MWTask.h"
#include "MWGroup.h"
#include "MWSystem.h"

#ifndef NWSENABLED
#include "MWNWSTask.h"
#endif

/* 8/28/00
   Qun Chen wrote signal handling code that Michael likes.
   We add it here */

#include <setjmp.h>

void statsgather ( MWList<MWWorkerID>* workers, MWStatistics * st );
MWTask * first_task_in_group(MWList<MWTask> * tasks, MWGroup * grp);

/// Indication of how many workers in each group.
int MWworkClasses = 0;
int *MWworkClassWorkers;
int *MWworkClassTasks;
MWWorkerID *MWcurrentWorker;


//by jae
pthread_mutex_t global_mw_mutex = PTHREAD_MUTEX_INITIALIZER;
int num_completed_tasks = 0;
// Signal code not used now
// jmp_buf env; 
// static void signal_handler(int sig){
  
//   switch ( sig ){
//   case SIGINT:       /* process for interrupt */
//     longjmp(env,sig);
//     /* break never reached */
//   case SIGCONT:     /* process for alarm */
//     longjmp(env,sig);
//     /* break never reached */
//   default:        exit(sig);
//   }
  
//   return;
// }
/* end Qun */

MWDriver::MWDriver() 
{
	stats = new MWStatistics();
	workers = new MWList<MWWorkerID>("workers");
	todo = new MWList<MWTask>("TODO_tasks");
	running = new MWList<MWTask>("running_tasks");

		/* Task_Swap */
	has_task_swapped = false;
	max_task_key = DBL_MAX;
	
		/* XXX DEBUG value for swap for ATR */
	MIN_IN_MEM = 256;    			/* == 256            */
	MAX_IN_MEM = MIN_IN_MEM * 800;   	/* == MIN_IN_MEM*800 */
	NUM_IN_MEM_LO = MIN_IN_MEM * 8; 	/* == MIN_IN_MEM*8   */
	NUM_IN_MEM_HI = MIN_IN_MEM * 100; 	/* == MIN_IN_MEM*100 */

		// XXX reassigned
	reassigned_tasks = new MWList<int>();
	reassigned_tasks_done = new MWList<int>();
	normal_tasks_all_done = false;
    
	task_counter = 0;
    
	checkpoint_frequency = 0;
	checkpoint_time_freq = 0;
	next_ckpt_time = 0;
	bench_task = NULL;

		// These are for list management
	addmode = ADD_AT_END;
	getmode = GET_FROM_BEGIN;
	machine_ordering_policy = NO_ORDER;
	listsorted = false;
	task_key = NULL;
	worker_key = NULL;

	ckpt_filename = "checkpoint";

	suspensionPolicy = DEFAULT;
  
	worker_timeout = false;
	worker_timeout_limit = 0.0;
	worker_timeout_check_frequency = 0;
	next_worker_timeout_check = 0;

	defaultTimeInterval = 10;
	defaultIterations = 5;
	defaultPortNo = 7;

	MWworkClasses = 0;
	MWworkClassWorkers = NULL;
	MWworkClassTasks = NULL;
	MWcurrentWorker = NULL;

#if defined( XML_OUTPUT )
  
	xml_menus_filename = "menus";
	xml_jobinfo_filename = "job_info";
	xml_pbdescrib_filename = "pbdescrib";
	xml_filename = "/u/m/e/metaneos/public/html/iMW/status.xml";

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
	memset(mach_name, 0 , sizeof(mach_name)); 

	get_machine_info();
#endif    

	normal_tasks_all_done = false;

#ifdef MEASURE
	strcpy(_measure_opt_file_name, "_measure_opt");
	strcpy(_measure_rec_fname_prefix, "_measure_rec");
	strcpy(_measure_ProcID_prefix, "_ProcID.");
	_measure_read_opt_interval = 1200;
	_measure_dump_rec_interval = 10;
	_measure_use_adaptation = 0;
	_measure_read_options(_measure_opt_file_name);
#endif //MEASURE

	m_pending_num = 0;
	pthread_mutex_init(&m_pending_mutex, NULL);
	pthread_cond_init(&m_pending_cond, NULL);
}


MWDriver::~MWDriver() 
{
	MWTask *t;
	MWWorkerID *w;
	
	while ( ( t = (MWTask *)todo->Remove() ) != NULL )
		delete t;
	delete todo;

	while ( ( t = (MWTask *)running->Remove() ) != NULL )
		delete t;
	delete running;

	while ( ( w = (MWWorkerID *)workers->Remove() ) != NULL )
		delete w;
	delete workers;

	if ( bench_task ) 
		{
			delete bench_task;
		}

	if ( MWworkClassWorkers )
		delete [] MWworkClassWorkers;
	if ( MWworkClassTasks )
		delete [] MWworkClassTasks;
	if ( stats )
		delete stats;
}

void 
MWDriver::go( int argc, char *argv[] ) 
{
	MWReturn ustat;

		// do some setup when we start:
	
	ustat = master_setup( argc, argv );

	if ( ustat != OK ) 
		{
			MWprintf(10, "master_setup did not return OK, exiting\n");
			RMC->exit(1);
			return; 
		}
    
		// the main processing loop:
	ustat = master_mainloop ();
	
	if( ustat == OK || ustat == QUIT || ustat == PENDING ) 
		{
				// We're done - now print the results.
				
			printresults();

				//JTL  Since kill_workers can take a long time --we like
				//   to collect the statistics beforehand.  Therefore, 
				//   we gather all the stats from the currently workign workers,
				//   then call makestats here (before kill_workers).
		
			statsgather(workers, stats);
				// tell the stats class to do its magic
			stats->makestats();
		}else{
			MWprintf(10,"no print results because ustat is %d  %d\n",ustat,HOSTRESUME);
		}

	kill_workers();

		// remove checkpoint file.
	//remove(ckpt_filename);

		/* Does not return. */
	RMC->exit(0);
}

void
statsgather ( MWList<MWWorkerID>* workers, MWStatistics * st )
{
	MWWorkerID * w = workers->First();
	while ( workers->AfterEnd() == false ) {
		w = workers->Current();
		st->gather(w);
		workers->Next();
	}
}

MWReturn 
MWDriver::master_setup( int argc, char *argv[] ) 
{
	MWReturn ustat = OK;
	int myid, master_id;
	int num_curr_hosts, i;
	MWWorkerID **tempWorkers = NULL;

	MWprintf ( 70, "MWDriver is pid %ld.\n", MWSystem::getpid() );
    
	RMC->setup( argc, argv, &myid, &master_id );

#ifdef MEASURE
	{
			// First set the file name for the measurement records
		int hasProcID = 0;
		for (i=0; i<argc; i++) {
			if ( strstr(argv[i], _measure_ProcID_prefix) != NULL ) {
					// Use the %(ProcID)
				hasProcID = 1;
				sprintf(_measure_rec_file_name, "%s%s", _measure_rec_fname_prefix, argv[i]);
			}
		}
		if (hasProcID == 0) {
				// No ProcID info from the arguments, so use getpid()
			sprintf(_measure_rec_file_name, "%s.%d", _measure_rec_fname_prefix, getpid());
		}
		MWprintf(31, "MEASURE|Setting measurement record fname = %s\n", _measure_rec_file_name);
	}
		
		// Now start counting until next time read-opt-file/dump_rec-file
	_measure_last_read_opt_time = _measure_last_dump_rec_time = time(0);
		// Write the header lines into the record file
	_measure_dump_header();

		// Init/Reset the measurement data
	_measure_reset();
#endif // MEASURE
    
		// See if checkpoint file exists
	FILE *chp_file = NULL;
	if (  (chp_file = fopen(ckpt_filename, "r")) == NULL)
		{
				// checkpoint file not found - start from scratch:
			MWprintf ( 50, "Starting from the beginning.\n" );
			MWprintf ( 50, "argc=%d, argv[0]=%s\n", argc, argv[0] );

			ustat = get_userinfo( argc, argv );
		
			if( ustat == OK ) 
				{
					if ( RMC->MW_exec_class_num_workers == NULL )
						{
							MWprintf ( 10, "In Get User Info set_exec_classes not called\n");
							MWprintf ( 10, "Aborting\n");
							assert(0);
						}
					ustat = create_initial_tasks();
						/* Memory for workers allocated in here: */
					RMC->init_beginning_workers ( &num_curr_hosts, &tempWorkers );
				} 
			else 
				{
					return ustat;
				}
		}
	else 
		{  // restart from that checkpoint file.
			fclose(chp_file);
			MWprintf ( 50, "Restarting from a checkpoint.\n" );

#if 0
				// Yet another Hack for Michael
			get_userinfo( argc, argv );
#endif
			restart_from_ckpt();
				/* Memory for workers allocated in here: */
			RMC->restart_beginning_workers ( &num_curr_hosts, &tempWorkers, RE_INIT );
		}

	MWprintf ( 70, "Good to go.\n" );

	if ( tempWorkers ) 
		{ /* there are indeed workers... */

			for ( i=0 ; i<num_curr_hosts ; i++ ) 
				{
					if ( tempWorkers[i]->get_id1() == -1 ) 
						{
								/* There was a problem with this worker! */
							tempWorkers[i]->started();
							tempWorkers[i]->ended();
							stats->gather( tempWorkers[i] );
						} 
					else 
						{
								/* Success! */
								/* Memory now taken care of by list... */
							addWorker ( tempWorkers[i] );
						}
				}
				/* don't forget to remove the array of ptrs .. */
			delete [] tempWorkers;
		}

	return OK;
}

MWReturn 
MWDriver::create_initial_tasks() 
{    
	int orig_ntasks;
	MWTask **orig_tasks;
    
	MWReturn ustat = setup_initial_tasks( &orig_ntasks, &orig_tasks );
	if ( MWworkClassWorkers == NULL )
		workClasses_set ( 1 );
	if( ustat == OK ) 
		{
			addTasks( orig_ntasks, orig_tasks );    
		}
	else 
		{
			MWprintf( 10, "setup_initial_tasks() returned %d\n", ustat );
		}
	delete [] orig_tasks;

	return ustat;
}


MWReturn 
MWDriver::master_mainloop() 
{

    call_hostaddlogic();
#ifdef INDEPENDENT

    ControlPanel ( );
    return OK;
#else

	MWReturn ustat = OK;

	while(1) {

		bool need_continue = false;

		// while there is at least one job in the todo or running queue:
		bool has_something = !is_TODO_empty() || 
			( !running->isEmpty() && (ustat == OK || ustat == PENDING) 
			  && !normal_tasks_all_done );

		if( has_something ) {

			// while there is at least one job in the todo or running queue:
			int buf_id = 0;
			int sending_host = 0;
			ustat = master_mainloop_oneshot(buf_id, sending_host);
			if (ustat != OK && ustat != PENDING ) {
				break;
			}

		}else {

			struct timespec ts;
			struct timeval tp;

			pthread_mutex_lock(&m_pending_mutex);
			while( m_pending_num > 0 ) {

				// Get current Time
				gettimeofday(&tp, NULL);

				// Convert from timeval to timespec
				ts.tv_sec = tp.tv_sec;
				ts.tv_nsec = tp.tv_usec * 1000;
				ts.tv_sec += 5; // pthread_cond_wait will wake up every 5 secs

				pthread_mutex_lock(&global_mw_mutex);
				if( !is_TODO_empty() ) {
					need_continue = true;
					pthread_mutex_unlock(&global_mw_mutex);
					break;
				}
				pthread_mutex_unlock(&global_mw_mutex);

				pthread_cond_timedwait(&m_pending_cond, 
						&m_pending_mutex, &ts);

				pthread_mutex_lock(&global_mw_mutex);
				if( !is_TODO_empty() ) {
					need_continue = true;
					pthread_mutex_unlock(&global_mw_mutex);
					break;
				}
				pthread_mutex_unlock(&global_mw_mutex);
			}
			pthread_mutex_unlock(&m_pending_mutex);

			if( need_continue ) {
				continue;
			}

			// Now, we have no todo, running, and pending jobs
			break;
		}
	}

	return ustat;
#endif
}

MWReturn
MWDriver::master_mainloop_oneshot(int buf_id, int sending_host)
{
	MWWorkerID *w = NULL;
	MWReturn ustat = OK;
	int len, tag, info;

	int tasknum=0;
	int subtasknum=0;
	MWTask *t;
#ifdef MEASURE
	double _task_recv_time = 0.0;
	double _task_recv_cpu_time = 0.0;
#endif


	MWprintf(90, "num_TODO = %d, num_run = %d, num_done = %d\n", 
			 todo->Number(), running->Number(), num_completed_tasks );

		/*	if (todo->Number() > MAX_IN_MEM)
			swap_out_todo_tasks(NUM_IN_MEM_LO, MAX_IN_MEM);
			else if (todo->Number() < MIN_IN_MEM)
			swap_in_todo_tasks(MIN_IN_MEM, NUM_IN_MEM_HI);
		*/	
	if( !todo->isEmpty() ) 
		MWprintf( 99, " CONTINUE -- todo list has at least task number: %d\n", todo->number()); 
	else if( !running->isEmpty() ) 
		MWprintf( 60, " CONTINUE -- running list has at least task number: %d\n", running->number() ); 

#ifdef MEASURE	
		// _measure: accumulate the _measure_master_recv_time value
	double _recv_start_time = _measure_cur_wall_time();
	double _recv_start_cpu_time = _measure_cur_cpu_time();
#endif // MEASURE

	// get any message from anybody
	//		MWList<void> * recv_buf_list = RMC->recv_buffers();

	// BUFFER: move to the next buffer
	//		RMC->recv_all(-1, -1);

	//		if (recv_buf_list == NULL || recv_buf_list->number()==0) {
	buf_id = RMC->recv( -1, -1 );
	//			MWprintf(91, "BUFFER: Have to use the blocking recv(), got buffer %d.\n", buf_id);
	//		} else {
	//			int* bid =  (int*) recv_buf_list->Remove();
	//			buf_id = *bid;
	//			recv_buf_list->Prepend(bid);
	//		}

#ifdef MEASURE	
	_measure_master_recv_time += _measure_cur_wall_time() - _recv_start_time;
	_measure_master_recv_cpu_time += _measure_cur_cpu_time() - _recv_start_cpu_time;
#endif // MEASURE

	if ( buf_id < 0 ) 
		{
			MWprintf ( 10, "ERROR, Problem with RMC->recv!\n" );
#ifdef INDEPENDENT
			return QUIT;
#else
				// GGT -- was continue
				//continue;
				
			//by jae	
			rematch_tasks_to_workers(NULL);
			return OK;
#endif
		}

	info = RMC->bufinfo ( buf_id, &len, &tag, &sending_host );

		/* 8/28/00
		   Qun Chen's signal handling code for Michael Ferris
		*/
/*

int returned_from_longjump, processing = 1;
	
if ((returned_from_longjump = setjmp(env)) != 0)
switch (returned_from_longjump)     {
case SIGINT:
MWprintf( 10, "!!!! Get user interruption.\n" );
checkpoint();
	//printresults();
	return QUIT;
	case SIGCONT: //Our program normally won't have this signal
	MWprintf( 10, "FATCOP was signaled to checkpoint...\n" );
	checkpoint();
	MWprintf( 10, "finish checkpointing ...\n" );
	break; 
	}
	(void) signal(SIGINT, signal_handler);
	(void) signal(SIGCONT, signal_handler);
		// alarm( 10 );
		*/
	switch ( tag ) 
		{
		case INIT: 
#ifdef MEASURE
			if (_measure_target_num_workers < workers->number()) {
				MWprintf ( 31, "MEASURE|Got INIT message, but we need to remove workers ...\n");
				RMC->initsend ( );
				RMC->send ( sending_host, KILL_YOURSELF );
				break;
			}
#endif // MEASURE
			if ( (w = lookupWorker ( sending_host )) == NULL ) 
				{
					MWprintf ( 10, "We got an INIT message from worker %08x "
							   "who we don't have records for.\n", sending_host);
					RMC->initsend ( );
					RMC->send ( sending_host, KILL_YOURSELF );
					break;
				}
			ustat = worker_init( w );
			break;

		case BENCH_RESULTS:
#ifdef MEASURE
			if (_measure_target_num_workers < workers->number()) {
				MWprintf ( 31, "MEASURE|Got BENCH_RESULTS message, but we need to remove workers ...\n");
				RMC->initsend ( );
				RMC->send ( sending_host, KILL_YOURSELF );
				break;
			}
#endif // MEASURE
		    if ( (w = lookupWorker ( sending_host )) == NULL ) 
				{
					MWprintf ( 10, "We got an BENCH_RESULTS message from (%08x), "
							   "who we don't have records for.\n", sending_host);
					RMC->initsend ( );
					RMC->send ( sending_host, KILL_YOURSELF );
					break;
				}
	
			ustat = handle_benchmark( w );
		    if ( ustat == OK ) 
				{
					printWorkers();
					//send_task_to_worker( w );		
					lockSendTaskToWorker( w );
				} 
		    else 
				{
					ustat = OK;  // we don't want to shut down...
				}
		    break;
		    
		case SUBRESULTS:
            if ( (w = lookupWorker ( sending_host )) == NULL )
				{                    
					MWprintf ( 10, "We got a RESULTS message from worker %08x who we don't have records for.\n",sending_host);
					RMC->initsend ();
					RMC->send ( sending_host, KILL_YOURSELF );
					break;
				}
				// deal with subresults
		    RMC->unpack( &tasknum, 1, 1 );
			MWprintf(51,"Received SUBRESULT message for tasknum %d\n",tasknum);
		    RMC->unpack( &subtasknum, 1,1);
		    t = rmFromRunQ( tasknum );
			if(t == NULL) MWprintf(51,"huh, task is not in runQ\n");
		    w->runningtask->unpack_subresults(tasknum);
				// computational steering
/*
  if(subtasknum==t->randomstop)
  {
  RMC->initsend();
  RMC->pack(&tasknum, 1, 1);
  RMC->send(sending_host, STOP_WORK);
  }
*/
		    act_on_completed_subtask(t);
				//act_on_completed_task(t);
		    putOnRunQ(t); // return task to run queue to complete the other tasks - with task container model, this task is done
			break;
		  
		case RESULTS:

			MWprintf ( 10, "We got a RESULTS message from worker %08x.\n", sending_host);

		    if ( (w = lookupWorker ( sending_host )) == NULL ) 
				{
	    		    MWprintf ( 10, "We got a RESULTS message from worker %08x who we don't have records for.\n", 
							   sending_host);
					RMC->initsend ();
					RMC->send ( sending_host, KILL_YOURSELF );
					break;
				}

		    ustat = handle_worker_results ( w );

		    if ( w->is_doomed() ) 
				{
						/* This worker has been marked for removal */
					worker_last_rites ( w );
				} 
		    else 
				{
					if ( ustat == OK ) 
						{
							//send_task_to_worker( w );		
							lockSendTaskToWorker( w );
						}
					else if( ustat != PENDING )
						{
							RMC->initsend ( );
							RMC->send ( w->get_id1(), KILL_YOURSELF );
							worker_last_rites ( w );
							call_hostaddlogic ( );
						}
				}

#ifdef MEASURE
    		    // Summarize the measurement results
			_measure_current_worker = w;
		    _measure();
#endif // MEASURE	
		    break;

		case HOSTADD:
			{
				MWWorkerID *w = new MWWorkerID;

				int r=0;
				r = RMC->start_worker( w );
				if ( r >= 0 )
					addWorker( w );
				else
					{
						delete w;
					}
		
				MWprintf(50, "Got a HOSTADD message.\n");
				call_hostaddlogic();
				break;
			}

		case HOSTDELETE:
			handle_hostdel();
			break;

		case TASKEXIT:
			handle_taskexit();
			break;

		case CHECKSUM_ERROR:
		    handle_checksum();
		    break;

		case HOSTSUSPEND:
			handle_hostsuspend();
			break;

		case HOSTRESUME:
			handle_hostresume();
			break;

		default:
			MWprintf ( 10, "Unknown message %d.  What the heck?\n", tag );

		} // switch
		
		/* Jeff's addition.  At this point, we need to do some assigning of
		   tasks.  If many tasks are added all in a batch, then many 
		   "SEND_WORK" commands are ignored, and the workers sit
		   idle.  Checking again here whether or not there are workers
		   waiting (or suspended) is a good idea.
		*/
	rematch_tasks_to_workers( w );

#if defined( SIMULATE_FAILURE )
		// Fail with probability p
	const double p = 0.0025;
	
	if( drand48() < p ) 
		{
			MWprintf ( 10, "Simulating FAILURE!\n" );
			kill_workers();
			RMC->exit(1);
		}
#endif	
	
		// If there is a worker_timeout_limit check all workers if there are
		// timed-out workers in working state

	if (worker_timeout)
		{
			time_t now = time(0);

			if (next_worker_timeout_check  <= now)
				{
					next_worker_timeout_check = now + worker_timeout_check_frequency; 
					MWprintf ( 80, "About to check timedout workers\n"); 
					reassign_tasks_timedout_workers();
				}
		}

	if( ustat != OK ) 
		{
			MWprintf( 10, "The user signalled %d to stop execution.\n", ustat );
		}
	
	return ustat;
}

int 
MWDriver::stop_work()
{
    int who;
    RMC->who ( &who );
	MWWorkerID *w = lookupWorker ( who );
	RMC->initsend();
	RMC->send ( w->get_id1(), STOP_WORK );
	return 0;
}

void
MWDriver::printresults ( )
{
	MWprintf ( 0, "In MWDriver packed %lf and unpacked %lf\n", RMC->get_bytes_packed(), 
			   RMC->get_bytes_unpacked() );
	MWprintf ( 0, "The MWDriver is done\n");
}

int
MWDriver::refreshWorkers ( int i, MWREFRESH_TYPE type )
{
		// send only to idle workers
	MWWorkerID *so = MWcurrentWorker;
	int j;

	if ( type == MW_NONE )
		return 0;
	if ( type == MW_THIS )
		{
			for ( j = 0; j < MWworkClasses; j++ )
				{
					if ( MWcurrentWorker->doesBelong ( j ) )
						MWcurrentWorker->deleteGroup ( j );
				}
			RMC->initsend ( );
			pack_worker_init_data ( );
			RMC->send ( MWcurrentWorker->get_id1(), REFRESH );
			return 0;
		}

	int num = 0;
	int max = workers->number();
	MWWorkerID *w;

	while ( num < max )
		{
			w = (MWWorkerID *)workers->Remove();
			if ( w != MWcurrentWorker && w->currentState() != IDLE )
				workers->Append ( w );
			else if ( !w->doesBelong ( i ) )
				workers->Append ( w );
			else
				{
						// You have to remove them from the list.
					for ( j = 0; j < MWworkClasses; j++ )
						{
							if ( w->doesBelong ( j ) )
								w->deleteGroup ( j );
						}
					RMC->initsend();
					MWcurrentWorker = w;
					pack_worker_init_data ( );
					RMC->send ( w->get_id1(), REFRESH );
					workers->Append ( w );
				}
			num++;
		}
	MWcurrentWorker = so;
	return 0;
}

void
MWDriver::worker_last_rites ( MWWorkerID *w )
{
	w->ended();
	stats->gather(w);
	rmWorker ( w->get_id1() );
	RMC->removeWorker ( w );
	RMC->MW_exec_class_num_workers[w->get_exec_class()]--;
	
		/* Commented out by Jichuan, to avoid calling deleteGroup() twice */
		/*
		  for ( int i = 0; i < MWworkClasses; i++ )
		  {
		  if ( w->doesBelong ( i ) )
		  w->deleteGroup ( i );
		  }
		*/
}

MWReturn 
MWDriver::worker_init( MWWorkerID *w ) 
{
	char workerhostname[64];
	MWReturn ustat = OK;

		/* We pack the hostname on the worker side... */
	RMC->unpack( workerhostname );
    
	MWprintf ( 40, "The MWDriver recognizes the existence of worker \"%s\".\n",
			   workerhostname );

		// 9/5/00 Jeff changes this at Jim Basney's behest.  
		//          Now, MWPvmRC will return the proper hostname, so we
		//          check to see if the name has already been set before
		//          setting it here
	
	w->set_machine_name( workerhostname );
	w->started();

	unpack_worker_initinfo( w );

	w->get_machine_info();

	int stat = RMC->initsend();
	if( stat < 0 ) 
		{
			MWprintf( 10, "Error in initializing send in pack_worker_init_data\n");
			worker_last_rites ( w );
			return OK;
		}

	if ( RMC->pack ( getHostName() ) != 0 )
		{
			MWprintf ( 10, "packing of hostName returned !OK.\n");
			worker_last_rites ( w );
			return ustat;
		}
	MWcurrentWorker = w;
	if ( (ustat = pack_worker_init_data()) != OK ) 
		{
			MWprintf ( 10, "pack_worker_init_data returned !OK.\n" );
			w->printself(10);
			worker_last_rites ( w );
			return ustat;
		}

	MWTask *b = get_benchmark_task();
	int benchmarking = b ? TRUE : FALSE;
	
	RMC->pack( &benchmarking, 1 );
	if ( b ) 
		{
			b->pack_work();
		}

	MWprintf ( 10, "Master sending an INIT_REPLY\n");
	RMC->send( w->get_id1(), INIT_REPLY );
		// Note the fact that the worker is doing a benchmark. 
	w->benchmark();
	return ustat;
}

MWReturn
MWDriver::handle_benchmark ( MWWorkerID *w ) 
{
	int okay = TRUE;
	double bres = 0.0;

	RMC->unpack( &okay, 1 );
	if ( okay == -1 ) 
		{
				/* Very bad! */
			MWprintf ( 10, "A worker failed to initialize!\n" );
			w->printself(10);
			worker_last_rites ( w );
			return QUIT;
		} 

	RMC->unpack( &bres, 1 );
	w->set_bench_result( bres );

	MWprintf( 10, "Received bench result %lf from %s (%08x)\n", 
			  bres, w->machine_name(), w->get_id1() );

	stats->update_best_bench_results( bres );

	w->benchmarkOver();

		//JTL -- Now you need to order the workerlist
	sort_worker_list();

	return OK;
}

void 
MWDriver::kill_workers() 
{

	MWprintf( 40, "Killing workers:\n");

    MWWorkerID *w = (MWWorkerID *)workers->Remove();
    
    while ( w ) 
		{
			RMC->initsend ( );
			RMC->send ( w->get_id1(), KILL_YOURSELF );
			w->ended();
			stats->gather(w);
			RMC->MW_exec_class_num_workers[w->get_exec_class()]--;
			for ( int i = 0; i < MWworkClasses; i++ )
				{
					if ( w->doesBelong ( i ) )
						w->deleteGroup ( i );
				}
			RMC->removeWorker( w );
		
			w = (MWWorkerID *)workers->Remove();
		}
}

void
MWDriver::give_work ( MWWorkerID *w, MWTask *t)
{
	RMC->initsend();
	int itmp = t->taskType;
	RMC->pack ( &itmp, 1, 1 );
	RMC->pack( &(t->number), 1, 1 );
	pack_driver_task_data();
	t->pack_work();
	RMC->send( w->get_id1(), DO_THIS_WORK );
}

void
MWDriver::lockSendTaskToWorker( MWWorkerID *w ) 
{
	pthread_mutex_lock(&global_mw_mutex);
	send_task_to_worker(w);
	pthread_mutex_unlock(&global_mw_mutex);
}

void 
MWDriver::send_task_to_worker( MWWorkerID *w ) 
{

		/* Jeff added this!  Due to the "rematch" function,
		 * it can be the case that a worker who asks for work might
		 * already have some work, in which case we don't want to give
		 * him another task.  (If he goes down then only one of the tasks
		 * will be rescheduled.
		 */

	MWTask *t = NULL;
	if( w->currentState() != IDLE ) 
		{
			return;
		}
		/* End of Jeff's kludgy fix */

	double lastMeasured;
	double nowTime;
	w->getNetworkConnectivity ( lastMeasured );
	nowTime = MWSystem::gettimeofday();

#ifndef NWSENABLED
/*
  if ( nowTime > lastMeasured + 60 * 120 )
  {
  MWprintf ( 10, "Creating NWSTask: nowTime = %f, lastMeasured = %f\n", nowTime, lastMeasured);
  t = new MWNWSTask ( task_counter++, defaultTimeInterval, defaultIterations, defaultPortNo, getHostName() );
	  // t->printself ( 10 );
	  MWprintf ( 10, "Sending MWNWSTask %f %d\n", defaultTimeInterval, defaultIterations);
	  }
	  else 
*/
		// in the future, there will be some policy to determine how many tasks go into one container
	int tasks_per_container = 1;
	for(int i = 0; i < tasks_per_container; i++)
		{
			t = getNextTask( w->getGroup() );
			if( t != NULL) {
				w->gottask(t->number);
				putOnRunQ(t);
			}
		}
#else
	int tasks_per_container = 1;
    for(int i = 0; i < tasks_per_container; i++)
		{
			t = getNextTask( w->getGroup() );
			if( t == NULL)
				break;
			tc->addTask(t);
			w->gottask(t->number);
			putOnRunQ(t);
		}
#endif

	if ( t != 0)
		{
			give_work ( w, t);
			w->runningtask = t;
			t->worker = w;

			MWprintf ( 10, "Sent following task to %s (%08x).\n",
					   w->machine_name(), w->get_id1() );
			MWprintf(10,"%d\n",t->number);
	
				//MWprintf ( 10, "Sent task (%d-%d) to %s  (%08x).\n",
				//       tc->FirstNum(), tc->LastNum(), w->machine_name(), w->get_id1() );
				//remember to delete tc when you are done
		}	
/*
  if ( t != NULL ) 
  {
  give_work ( w, tc );
  w->gottask( t->number );
  w->runningtask = t;
  t->worker = w;
  putOnRunQ( t );
        
  MWprintf ( 10, "Sent task %d to %s  (%08x).\n", 
  t->number, w->machine_name(), w->get_id1() );

  }
  else 
  {
  this added to make things work better at program's end:
  The end result will be that a suspended worker will have its
  task adopted.  We first look for a task with a suspended worker.
  We call it t.  The suspended worker's 'runningtask' pointer 
  will still point to t, but t->worker will point to this
  worker (w).  Now there will be two workers working on the 
  same task.  The first one to finish reports the results, and
  the other worker's results get discarded when (if) they come in.
		

  if ( !running->isEmpty() ) 
  {
	  // here there are no tasks on the todo list, but there are 
		  // tasks in the run queue.  See if any of those tasks are
			  // being held by a suspended worker.  If so, adopt that task.
			
				  // will reimplement later
					  //reassignSuspendedTask(running, w);
					  }

					  XXX Jichuan's hack to avoid the master and workers busy working on 
					  * NWSTasks only, this happens if there are many workers starting from 
					  * different time which happen to make this cycle of NWSTasks overlap 
					  * with the next. So we will check to see if all the tasks in runQ are
					  * NWSTasks, and stop the application if so! 
					  *
		
					  XXX A fix to solve the waiting for already finished task problem. 
					  * The problem is that because a task can be assigned to more than one workers, 
					  * it can happen that the some of the workers already return the result, but
					  * either the master is waiting for the other workers to finish (if the task is
					  * the last one), or the master will reassign the task to more workers. 
					  * 
					  * To fix it, we use two marker lists to record the set (no dup) of tasks that 
					  * can were reassigned and already finished: List "reassigned_tasks" and List
					  * "reassigned_tasks_done". A task id is inserted to "reassigned_tasks" when 
					  * a its worker is suspended or timed-out, and inserted to "reassigned_tasks_done"
					  * when a result is reported back. The difference of two set can be problematic
					  * tasks.
					  * 
					  * If the TODO list is empty, and all running task is in the "reassigned_tasks_done"
					  * set, then QUIT; if a task is suspended or timed-out, and it's in the 
					  * "reassigned_tasks_done" set, it will not be reassigned.
					  *
		
					  bool all_NWS = true, all_reassigned = true, all_done = true;
		
						  // MWprintf ( 60, "Now that no task in TODO list, let's see if all tasks in the runQ are NWSTaks.\n");
						  MWTask *tmp_t = running->First();

						  while (running->AfterEnd() == false) {
						  tmp_t = running->Current();
						  if (tmp_t->taskType != MWNWS )
						  all_NWS = false;
						  if (!in_set(reassigned_tasks, tmp_t->number))
						  all_reassigned = false;
						  else if (!in_set(reassigned_tasks_done, tmp_t->number))
						  all_done = false;
						  running->Next();
						  }

						  if (all_NWS) {
						  MWprintf( 10, "Since all running tasks are NWSTask, and no TODO task, we are done! \n");
						  normal_tasks_all_done = true;
						  } else if (all_reassigned && all_done) {
						  MWprintf( 10, "Since all running tasks are reassigned and all done, and no TODO task, we are done! \n");
						  normal_tasks_all_done = true;
						  }
						  }
*/
}

void
MWDriver::reassignSuspendedTask ( MWList<MWTask> * running, MWWorkerID * w )
{
	MWGroup *grp1, *grp2;
	grp2 = w->getGroup();

	MWTask *t = running->First();
	while (running->AfterEnd() == false) {
		t = running->Current();
		grp1 = t->getGroup();
	
		if ( grp1->doesOverlap ( grp2 ) && t->worker->isSusp() ) 
			{
		       	MWprintf ( 60, "Re-assigning task %d to %s (%08x).\n", 
						   t->number, w->machine_name(), w->get_id1() );
		       	t->worker->printself( 60 ); // temp
		       	w->gottask( t->number );
		       	t->worker = w;
		       	w->runningtask = t;
					// note that old wkr->runningtask points to t also.
				give_work ( w, t );
		       	break;
	       	}

		running->Next();
	}
}

void 
MWDriver::handle_hostdel () 
{
	
		// Here we find out a host has gone down.
		// this tells us who died:
	int who;
	RMC->who ( &who );

	MWWorkerID *w = rmWorker( who );
	if ( w == NULL ) 
		{
			MWprintf ( 40, "Got HOSTDELETE of worker %08x, but it isn't in the worker list.\n", who );
		}
	else 
		{
			MWprintf ( 40, "Got HOSTDELETE!  \" %s \" (%08x) went down.\n", 
					   w->machine_name(), who );   

			hostPostmortem ( w );
			MWprintf ( 40, "There are now %d workers\n", workers->number() );

				/* We may have to ask for some hosts at this time... */
			delete w;
			call_hostaddlogic();
		}
}

void 
MWDriver::handle_taskexit () 
{

		/*   I have found that we usually get this message just before a HOSTDELETE, 
			 so I'm not going to do anything like remove that host...
		*/

	int who;
	RMC->who ( &who );
		// search workers for that id2.
	MWWorkerID *w = rmWorker ( who );
	if ( w == NULL ) 
		{
			MWprintf ( 40, "Got TASKEXIT of worker %08x, but it isn't in "
					   "the worker list.\n", who );
			return;
		}

	MWprintf ( 40, "Got TASKEXIT from machine %s (%08x).\n", 
			   w->machine_name(), w->get_id1());

	hostPostmortem( w );
	RMC->removeWorker ( w );
		/* We may have to ask for some hosts at this time... */
	call_hostaddlogic();
}

MWReturn
MWDriver::act_on_starting_worker ( MWWorkerID *wid )
{
	return OK;
}

void
MWDriver::handle_checksum ( )
{
		// Some worker returned us a checksum error message
		// Actually there may be some bootstrapping problem here
	int who;
	RMC->who ( &who );

	MWReturn ustat = OK;
	MWWorkerID *w = lookupWorker ( who );

	MWprintf( 10, "Worker \" %s \" (%08x) returned a checksum error!\n", 
			  w->machine_name(), who );

	if ( w == NULL )
		{
			MWprintf ( 10, "Got CHECKSUM_ERROR from worker %08x, but it isn't in "
					   "the worker list.\n", who );
			return;
		}
	else
		{
			switch ( w->currentState ( ) ) 
				{
				case BENCHMARKING:
					{
						int stat = RMC->initsend();
						if( stat < 0 ) 
							{
								MWprintf( 10, "Error in initializing send in pack_worker_init_data\n");
							}
						MWcurrentWorker = w;
						if ( (ustat = pack_worker_init_data()) != OK ) 
							{
								MWprintf ( 10, "pack_worker_init_data returned !OK.\n" );
								w->printself(10);
								w->started();
								w->ended();
								RMC->MW_exec_class_num_workers[w->get_exec_class()]--;
								rmWorker ( w->get_id1() );
								RMC->removeWorker ( w );
								return;
							}

						MWTask *b = get_benchmark_task();
						int benchmarking = b ? TRUE : FALSE;

						RMC->pack( &benchmarking, 1 );
						if ( b ) 
							{
								b->pack_work();
							}

						RMC->send( w->get_id1(), INIT_REPLY );
							// Note the fact that the worker is doing a benchmark.
						w->benchmark();
						return;
						break;
					}
				case WORKING:
					{
						if ( w->runningtask )
							{
									// It was running a task
									// See if this task is there in the done list
									// implement this later
									//	if ( w->runningtask->taskType == MWNORMAL )
									//		pushTask ( w->runningtask );

								//send_task_to_worker( w );		
								lockSendTaskToWorker( w );
								return;
							}
						else
							{
									// How is this possible that a worker is not working on something
								MWprintf ( 10, "Worker was not working on any task but still gave a CHECKSUM_ERROR\n");
								return;
							}
						break;
					}
			
				default:
					MWprintf(30, "MWDriver::handle_checksum() - other worker state.\n");
				}
		}
}

void 
MWDriver::handle_hostsuspend () 
{
	int who;
	RMC->who ( &who );
	MWprintf(40,"Got HOSTSUSPEND\n");

	MWWorkerID *w = lookupWorker ( who );
	if ( w != NULL ) 
		{
			w->suspended();
			MWprintf ( 40, "Got HOSTSUSPEND!  \"%s\" is suspended.\n", w->machine_name() );   

			int task;
			if ( w->runningtask && w->runningtask->taskType == MWNORMAL )
				{
					task = w->runningtask->number;
					MWprintf ( 60, " --> It was running task %d.\n", task );

					switch ( suspensionPolicy ) 
						{
						case REASSIGN: 
							{
									// remove task from running list; push on todo.

// 						 First see if it's *on* the running list.  If it's 
// 						not, then it's either back on the todo list or
//						its already done.  If this is the case, ignore. 

								MWTask *rmTask = rmFromRunQ( task );
								if ( !rmTask ) 
									{
										MWprintf ( 60, "Not Re-Assigning task %d.\n", task );
										break;
									}
								MWTask *rt = w->runningtask;
								if ( rt != rmTask ) 
									{
										MWprintf ( 10, "REASSIGN: something big-time screwy.\n");
									}

// 						 Now we put that task on the todo list and rematch
// 						the tasks to workers, in case any are waiting around.

								if (!in_set(reassigned_tasks_done, rt->number)) {
									insert_into_set(reassigned_tasks, rt->number);
									MWprintf(60,"REASSIGN: Will reassign task %d.\n", rt->number);
									pushTask( rt );
										// just in case...
									rematch_tasks_to_workers ( w );
								}
								break;
							}
						case DEFAULT: 
							{
									// we only take action if the todo list is empty.
								if ( todo->isEmpty() ) 
									{
											// here there are no more tasks left to do, we want
											// to re-allocate this job.  First, find idle worker:
											//will reimplement later
											//reassignIdleTask(workers, w);
									}
							}
						} // switch
				}
			else if ( w->runningtask && w->runningtask->taskType != MWNORMAL )
				{
					MWprintf ( 60, " --> It was running a control task.\n");
				}
			else 
				{
					MWprintf ( 60, " --> It was not running a task.\n" );
				}
    	}
	else
		{
			MWprintf ( 40, "Got HOSTSUSPEND of worker %08x, but it isn't in "
					   "the worker list.\n", who );
		}
}

void
MWDriver::reassignIdleTask (MWList<MWWorkerID> * workers, MWWorkerID * w)
{
	MWGroup *grp1, *grp2; 
	grp2 = w->runningtask->getGroup();

	MWWorkerID *o = workers->First();
	while (workers->AfterEnd() == false) {
		o = workers->Current();
		grp1 = o->getGroup();
		
		if ( o->isIdle() && grp1->doesOverlap ( grp2 ) ) 
			{
		       	MWprintf ( 60, "Reassigning task %d. to %s\n", 
						   w->runningtask->number, o->machine_name() );
		       	w->runningtask->worker = o;
		       	o->runningtask = w->runningtask;
		       	o->gottask( o->runningtask->number );

				give_work ( o, o->runningtask );
				break;
	       	}
		o = workers->Next();
	}
}

void 
MWDriver::handle_hostresume () 
{

		/* Do nothing; applaud politely */

	int who;
	RMC->who ( &who );

	MWWorkerID *w = lookupWorker ( who );
	if ( w != NULL ) 
		{
			w->resumed();
			MWprintf ( 40, "Got HOSTRESUME.  \"%s\" has resumed.\n", 
					   w->machine_name() );
			if ( w->runningtask != NULL ) 
				{
					MWprintf ( 60, " --> It's still running task %d.\n", 
							   w->runningtask->number );
				}
			else 
				{
					MWprintf ( 60, " --> It's not running a task.\n" );
				}
		}
	else 
		{
			MWprintf ( 40, "Got HOSTRESUME of worker %08x, but it isn't in "
					   "the worker list.\n", who );
		}
}

void 
MWDriver::hostPostmortem ( MWWorkerID *w ) 
{
	MWTask *t = 0;
	w->ended();
	RMC->MW_exec_class_num_workers[w->get_exec_class()]--;
	for ( int i = 0; i < MWworkClasses; i++ )
		{
			if ( w->doesBelong ( i ) )
				w->deleteGroup ( i );
		}
	t = w->runningtask;	

	if ( t == NULL ) 
		{
			MWprintf ( 60, " %s (%08x) was not running a task.\n", 
					   w->machine_name(), w->get_id1() );
		}
	else if ( t->taskType == MWNORMAL )
		{
			MWprintf ( 60, " %s (%08x) was running task %d.\n", 
					   w->machine_name(), w->get_id1(), t->number );
				// remove that task from the running queue...iff that
				// task hasn't ALREADY been reassigned.  

			w->runningtask = NULL;
			MWWorkerID *o;
			if ( !(o = task_assigned( t )) ) 
				{
					MWTask *check = rmFromRunQ( t->number );
					if ( t == check ) 
						{
							t->worker = NULL;
							rmFromRunQ(t->number);
							pushTask( t );
						}
				}		
			else 
				{
					t->worker = o;
					MWprintf ( 60, " Task %d already has been reassigned\n", 
							   t->number );
				}
		}
	else
		{
				// It is a control task
				// XXX Maybe also need to delete the task from RunQ, Insure complained about this.
			MWTask *check = rmFromRunQ( t->number );
			if ( t == check )	// GGT gcc complains about t used before set here, I think this is a bug
				MWprintf(21, "The control task IS in the RunQ, and should haven't been deleted!");
			delete t;
		}
	stats->gather( w );
}
		
MWReturn 
MWDriver::handle_worker_results( MWWorkerID *w ) 
{
	MWReturn ustat = OK;

	int tasknum;  
	double wall_time = 0;
	double cpu_time = 0;
    
	RMC->unpack( &tasknum, 1, 1 );

#ifdef MEASURE
		// get the _task_recv_time 
	double _task_recv_time = 0.0;
	double _task_recv_cpu_time = 0.0;
	RMC->unpack( &_task_recv_time, 1, 1);
	RMC->unpack( &_task_recv_cpu_time, 1, 1);
#endif // MEASURE
	
	RMC->unpack( &wall_time, 1, 1 );
	RMC->unpack( &cpu_time, 1, 1 );

		// We once received CPU times that were NAN from a worker,
		// and this messes EVERYTHING up

	if( cpu_time != cpu_time ) 
		{
			MWprintf( 10, "Error.  cpu_time of %lf, assigning to wall_time of %lf\n",
					  cpu_time, wall_time );
		}
	
		// Jeff's CPU time sanity check!!!!
	if( cpu_time > 0.1 && (cpu_time > 1.1 * wall_time || cpu_time < 0.01 * wall_time ))
		{
			MWprintf( 10, "Warning.  Cpu time =%lf,  wall time = %lf seems weird.\n", 
					  cpu_time, wall_time );
			MWprintf( 10, "Assigning cpu time to be %lf\n", wall_time );
			cpu_time = wall_time;
		}
		//	total_time = now - start_time;	


	MWprintf ( 60, "Unpacked task %d result from %s (%08x).\n", 
			   tasknum, w->machine_name(), w->get_id1() );

	MWprintf ( 81, "It took %5.2f secs wall clock and %5.2f secs cpu time\n", 
			   wall_time, cpu_time );

	int mwtasknum = w->runningtask ? w->runningtask->number : -2;
	if( tasknum != mwtasknum ) 
		{
			MWprintf( 10, "What the $%!^$!^!!!.  MW thought that %s (%08x) was running %d, not %d\n", 
					  w->machine_name(), w->get_id1(), mwtasknum, tasknum );
			return QUIT;
		}

	if ( w->isSusp() ) 
		{
			MWprintf ( 60, "Got results from a suspended worker!\n" );
		}

	MWTask *t = rmFromRunQ( tasknum );

		/* The task *could* have been suspended, REASSIGNed, and put
		   on the todo list.  Let's just check the first entry on the
		   todo list and see if it's there.  If so, remove it and handle. */
	if ( suspensionPolicy == REASSIGN )
		{
			MWTask * tmp_t = todo->First();
			while (todo->AfterEnd() == false) {
				tmp_t = todo->Current();
				if (tmp_t->number == t->number)
					todo->RemoveCurrent();
				else todo->Next();
			}
		}

		// in_TODO inited to be 0
	int in_TODO = 0;
	
	if ( t == NULL ) {		
			// XXX Check to see if the task is in TODO list
			// It is possible that when re-assigning a task for many times, 
			// the task can be temporarily in TODO list, but a worker at the 
			// same time can return the result for this task. In this case, 
			// we will just process the result (t = tmp_t) and remove the task from TODO 
			// list. 
			// 		-- Reported by Jeff
			// 		-- temp fix by Jichuan
	
		if (in_set(reassigned_tasks_done, tasknum))
			MWprintf ( 60, "Task %d already done, no need to process it.\n", tasknum);
		else {
			if (!in_set(reassigned_tasks, tasknum)) 
				MWprintf ( 30, "Something wrong with REASSIGN, task %d.\n", tasknum);
			
			MWTask * tmp_t = todo->First();
			while (todo->AfterEnd() == false) {
				tmp_t = todo->Current();
				if (tmp_t->number == tasknum) { /* found */
					todo->RemoveCurrent();
					in_TODO = 1;
					t = tmp_t;
					MWprintf(60, "Done task %d found in TODO, delete from it.\n", tasknum);
				} else todo->Next();
			}
		}
	}
		
	if ( ( !in_TODO ) && ( t == NULL ) ) {
		
		w->completedtask( wall_time, cpu_time );

			// XXX Delete the task if there is no other workers are running it. 
			// It's possible that the same task is assigned to many different
			// workers (everytime a worker suspend, it's running task may be 
			// reassigned to another worker). 
			// 
			// 	if ( w->runningtask ) 
			//
			// The original statement assumes that the number of workers running
			// task t is maximally 2, so when the second worker report result, 
			// MW delete the task. If a new task is created using the same region
			// of memory, the runningtask pointer is corrupted. 

		if (w->runningtask ) {
			if (task_assigned(w->runningtask) == NULL) {
				delete w->runningtask;
			}
			w->runningtask = NULL;
		}
	}
	else 
		{
				// tell t the results
			t->unpack_results();
				/* give these to the task class; now the user can use them. */
			t->working_time = wall_time;
			t->cpu_time = cpu_time;
			stats->update_best_task_results( cpu_time);

#ifdef MEASURE
				// set the task's _measure_MP_worker_time
			t->_measure_MP_worker_time = _task_recv_time;
			t->_measure_MP_worker_cpu_time = _task_recv_cpu_time;
#endif // MEASURE
	
			if ( t->taskType == MWNORMAL )
				{
					MWcurrentWorker = w;
#ifdef MEASURE
						// _measure: accumulate the _measure_master_act_on_completed_task_time (and cpu)values
					double act_start_time = _measure_cur_wall_time();
					double act_start_cpu_time = _measure_cur_cpu_time();
#endif // MEASURE

					ustat = act_on_completed_task( t );
#ifdef MEASURE
					_measure_master_act_on_completed_task_time += _measure_cur_wall_time() - act_start_time;
					_measure_master_act_on_completed_task_cpu_time += _measure_cur_cpu_time() - act_start_cpu_time;
					_measure_num_task_Done ++;
#endif // MEASURE
			
					if (in_set(reassigned_tasks, t->number)) {
						insert_into_set(reassigned_tasks_done, t->number);
						MWprintf(60, "REASSIGN: Task %d once reassigned and now finished.\n", t->number);
					}
			
					num_completed_tasks++;

					if ( (checkpoint_frequency > 0) && 
						 (num_completed_tasks % checkpoint_frequency) == 0 ) 
						{
							checkpoint();
						} 
					else 
						{
							if ( checkpoint_time_freq > 0 ) 
								{
									time_t now = time(0);
									MWprintf ( 91, "Trying ckpt-time, %d %d\n", 
											   next_ckpt_time, now );
									if ( next_ckpt_time <= now ) 
										{
											next_ckpt_time = now + checkpoint_time_freq;
											MWprintf ( 80, "About to...next_ckpt_time = %d\n", 
													   next_ckpt_time );
											checkpoint();
										}
								}
						}
				}
			else
				{
#ifndef NWSENABLED
					if ( t->taskType == MWNWS )
						{
							if ( w ) 
								w->setNetworkConnectivity ( ((MWNWSTask *)t)->median );
							MWprintf ( 10, "Got results of NWSTask\n");
						}
#endif
				}

			if ( t->worker == NULL ) 
				{
					MWprintf ( 40, "Task had no worker...huh?\n" );
				}

			if ( w != t->worker ) 
				{
					MWprintf ( 90, "This task not done by 'primary' worker, "
							   "but that's ok.\n" );
				}
			else 
				{
					MWprintf ( 90, "Task from 'primary' worker %s.\n", 
							   w->machine_name() );
				}

			//by jae
			if( ustat != PENDING ) {	
				w->runningtask = NULL;
				w->completedtask( wall_time, cpu_time );
				t->worker = NULL;

				if( ! task_assigned( t ) ) 
				{
					MWprintf( 80, "Deleting task %d\n", t->number );

					for ( int tempi = 0; tempi < MWworkClasses; tempi++ )
						{
							if ( t->doesBelong ( tempi ) )
								t->deleteGroup ( tempi );
						}

#ifdef MEASURE
						// When a task is done, accumulate the following values:
						// 	_measure_num_task_Done, task-related summary values
						// _measure_task_RoundTrip_time += _measure_cur_wall_time() - t->_measure_start_time;
					_measure_task_Exe_wall_time += t->working_time;
					_measure_task_Exe_cpu_time += t->cpu_time;
						// _measure_task_MP_master_time += t->_measure_MP_master_time;
					_measure_task_MP_worker_time += t->_measure_MP_worker_time;
					_measure_task_MP_worker_cpu_time += t->_measure_MP_worker_cpu_time;
#endif // MEASURE
					delete t;
				}
			}
		}
	return ustat;
}

void 
MWDriver::rematch_tasks_to_workers( MWWorkerID *nosend )
{
	pthread_mutex_lock(&global_mw_mutex);

		// For each worker that is "waiting for work", you can safely(?)
		// now do a send_task_to_worker.  Don't send it to the id "nosend"
		// (He just reported in).  
	MWWorkerID *w;
	int ns = nosend ? nosend->get_id1() : -1;
    
	w = workers->First();
	while (workers->AfterEnd() == false) {
		w = workers->Current();
		if ( (w->currentState() == IDLE) && ( (w->get_id1() != ns) || (w->get_id2() != ns) ) ) {
			MWprintf ( 90, "In rematch_tasks_to_workers(), trying send...\n" );
			send_task_to_worker(w);		
		}
		workers->Next();
	}	

	pthread_mutex_unlock(&global_mw_mutex);
}

MWTask * 
MWDriver::getNextTask( MWGroup *grp )  
{
	MWTask *t = NULL;
		// return next task to do.  NULL if none.
	if ( todo->isEmpty() ) 
		return NULL;
	
	switch( getmode ) {
		
	case GET_FROM_BEGIN: 
		{
			t = first_task_in_group(todo, grp);
			break;
		}
	case GET_FROM_KEY: 
		{
			if( task_key == NULL ) 
				{
					MWprintf ( 10, " GET_FROM_KEY retrieval mode, "
							   "but no key function set.\n"); 
					MWprintf ( 10, " Returning task at head of list\n");
					t = first_task_in_group(todo, grp);
				}
			else 
				{
					if ( listsorted )
						t = first_task_in_group(todo, grp);
					else
						{
							MWKey retval;

							if( task_key == NULL )
								{
									MWprintf( 10, " return_bst_keyval:  no task key "
											  "function defined!\n");      
									assert(0);
								}
							else
								{
									retval = -DBL_MAX;
									t = NULL;
									MWListElement<MWTask> * pos = NULL;
									todo->ToHead();
									while (todo->AfterEnd() == false) {
										MWTask * tt = todo->Current();
										MWKey tmp = (*task_key) (tt);
										MWGroup * g = tt->getGroup();
										if ( (retval < tmp) && g->doesOverlap ( grp ) ) {
											pos = todo->CurPosition();
											retval = tmp;
											t = tt;
										}
										todo->Next();
									}
									todo->Remove(pos);
								}
						}
				}
			break;
		}
	default: 
		{
			MWprintf ( 10, "Unknown getmode %d\n", getmode );
			assert( 0 );
		}
	} // switch
	return t;
}

MWTask * 
first_task_in_group ( MWList<MWTask> * tasks, MWGroup * grp)
{
	MWTask *t = tasks->First();
	MWGroup *grp2;
	while (tasks->AfterEnd() == false) {
		t = tasks->Current();
		grp2 = t->getGroup();
		if ( grp->doesOverlap ( grp2 ) ) {
			tasks->RemoveCurrent();
			return t;
		}
		else tasks->Next();
	}
	return NULL;
}

void 
MWDriver::pushTask( MWTask *push_task ) 
{
	switch( addmode ) 
		{
		case ADD_AT_BEGIN: 
			{
				todo->Prepend ( push_task );
				break;
			}
		case ADD_AT_END: 
			{
				todo->Append ( push_task );
				break;
			}
		case ADD_BY_KEY: 
			{
					// This will insert at the first key WORSE than the 
					// given one
				if( task_key == NULL ) 
					{
						MWprintf ( 10, "ERROR!  Adding by key, but no key "
								   "function defined!\n" );
						assert(0);
					}

				if ((*task_key)(push_task)<=max_task_key) {
					todo->SortedInsert ( push_task, (*task_key)(push_task) );
				}
				else {
					MWprintf(10, "Won't insert task whose key value %f > max_task_key %f.\n",
							 (*task_key)(push_task), max_task_key);
				}
				break;
			}
		default: 
			{
				MWprintf ( 10, "What the heck kinda addmode is %d?\n", addmode );
				assert( 0 );
			}
		} // switch
}

void 
MWDriver::addTasks( int n, MWTask **add_tasks )
{
	if ( n <= 0 ) 
		{
			MWprintf ( 10, "Please add a positive number of tasks!\n" );
			RMC->exit(1);
		}
    
	for ( int i = 0; i < n; i++ )
		{
			add_tasks[i]->number = task_counter++;
			if ( MWworkClasses <= 1 )
				{
					add_tasks[i]->initGroups ( MWworkClasses );
					add_tasks[i]->addGroup ( 0 );
				}

			pushTask ( add_tasks[i] );
		}
}

void 
MWDriver::addTask( MWTask *add_task )  
{
	add_task->number = task_counter++;
	if ( MWworkClasses <= 1 )
		{
			add_task->initGroups ( MWworkClasses );
			add_task->addGroup ( 0 );
		}
	pushTask ( add_task );
}

void 
MWDriver::ckpt_addTask( MWTask *add_task )  
{
	pushTask ( add_task );
}

void 
MWDriver::set_task_key_function( MWKey (*t)( MWTask * ) ) 
{
	if ( ( t != task_key ) && (task_key == NULL) ) {
		sort_task_list();
		listsorted = true;
	} else if ( t != task_key )
		listsorted = false;
		
	task_key = t;
}

int 
MWDriver::set_task_add_mode( MWTaskAdditionMode mode )
{
	int errcode = 0;
	if( mode < ADD_AT_END || mode > ADD_BY_KEY ) {
		errcode = UNDEFINED;
		MWprintf( 10, " Bad Add Task Mode.  Switching to ADD_AT_END\n");
		addmode = ADD_AT_END;
    }
	else {
		addmode = mode;
	}
	return errcode;
}

int 
MWDriver::set_task_retrieve_mode( MWTaskRetrievalMode mode )
{
	int errcode = 0;
	if( mode < GET_FROM_BEGIN || mode > GET_FROM_KEY ) {
		errcode = UNDEFINED;
		MWprintf( 10, " Bad Retrive Task Mode.  "
				  "Switching to GET_FROM_BEGIN\n" );
		getmode = GET_FROM_BEGIN;
    }
	else {
		getmode = mode;
		if( ( mode == GET_FROM_KEY ) && ( task_key == NULL ) ) {
			MWprintf( 10, " Warning.  Call set_task_key() function before ");
			MWprintf( 10, " retrieving task\n");
		}
    }
	
	return errcode;
}

int 
MWDriver::set_machine_ordering_policy( MWMachineOrderingPolicy mode ) 
{
	int errcode = 0;
	if( mode < NO_ORDER || mode > BY_KFLOPS ) {
		errcode = UNDEFINED;
		MWprintf( 10, "Bad machine ordering policy %d.  Switching to NO_ORDER\n", 
				  (int) mode );
		machine_ordering_policy = NO_ORDER;
		worker_key = NULL;
	}
	else {
		MWprintf( 30, "Set machine ordering policy to %d\n", (int) mode );
		machine_ordering_policy = mode;
		if( mode == BY_USER_BENCHMARK ) {
			worker_key = &benchmark_result;
		}
		else {
			worker_key = &kflops;
		}

	}
	return errcode;
}


void 
MWDriver::set_worker_timeout_limit( double timeout_limit, int timeout_freq ) 
{
	if (timeout_limit > 0)
		{
			worker_timeout = true;
			worker_timeout_limit = timeout_limit;

			worker_timeout_check_frequency = timeout_freq ;
			next_worker_timeout_check = time(0) + timeout_freq ;

			MWprintf( 30, "Set worker timeout limit to %lf\n", timeout_limit);
			MWprintf( 30, "Set worker timeout frequency to %d\n", timeout_freq);    

		}
	else 
		{
			MWprintf(10, "Urgh, Timeout limit has to be > 0\n");
		}
}

void MWDriver::reassign_tasks_timedout_workers()
{
	MWprintf ( 10, "Enter : reassign_tasks_timedout_workers()\n");
	MWWorkerID *tempw;
	int total = 0;

	if (worker_timeout_limit > 0) 
		{
			double now = MWSystem::gettimeofday();
			MWprintf( 60, "Now: %lf\n", now );

			while( total < workers->number() ) 
				{
					total++;
					tempw = (MWWorkerID *)workers->Remove();
					if (tempw->currentState() == WORKING) 
						{
							MWprintf( 50, "Last event from %s (%08x) : %lf\n", tempw->machine_name(), 
									  tempw->get_id1(), tempw->get_last_event() );

							if ( (now - tempw->get_last_event() ) > worker_timeout_limit )
								{
										// this worker is timed out
									if (tempw->runningtask)
										{
											int task = tempw->runningtask->number;

											MWTask *rmTask = rmFromRunQ( task );
											if ( rmTask ) 
												{
													MWprintf ( 10, "Worker %08x is timed-out\n",tempw->get_id1());
													MWprintf ( 10, "Task %d rescheduled \n",rmTask->number);
												}
											else 
												{
													MWprintf ( 60, "Task not found in running list; assumed done or in todo list.\n" );
												}
						
											if ( !rmTask || rmTask->taskType != MWNORMAL ) 
												{
													MWprintf ( 60, "No need to reassign task %d.\n", task );
												}
											else 
												{
													MWTask *rt = tempw->runningtask;
													if ( rt != rmTask ) 
														{
															MWprintf ( 10, "REASSIGN: something big-time screwy.\n");
														}
														// We don't need to do that.
														// tempw->runningtask = NULL;

														// Now we put that task on the todo list and rematch
														// the tasks to workers, in case any are waiting around.
													if (!in_set(reassigned_tasks_done, rt->number)) {
														insert_into_set(reassigned_tasks, rt->number);
														MWprintf(60, "REASSIGN: Will assign task %d.\n", rt->number);
														pushTask( rt );
															// just in case...
														rematch_tasks_to_workers ( tempw );
													}
												}
										}
								}
						}
					workers->Append ( tempw );
				}

/*
  This is Jeff's code.  We go through the *running* list 
  too to see if there is a task there that needs to be reassigned.
  For example -- if there is a task in the running list
  not assigned to anyone.  This can happen.  For example...
*/		

			total = 0;
			MWTask *tempr;

			while( total < running->number() ) 
				{
					total++;
					tempr = (MWTask *)running->Remove();
					if ( tempr->taskType != MWNORMAL )
						{
							running->Append ( tempr );
							continue;
						}

						// Find the worker running this task.
					bool foundworker = false;
					int subtotal = 0;
					MWWorkerID *tempw = NULL;
			
						// XXX // while( subtotal < workers->number() ) 
					while( (foundworker == false) && (subtotal < workers->number()) ) 
						{
							subtotal++;
							tempw = (MWWorkerID *)workers->Remove();
							if( tempw->runningtask ) 
								{
									MWprintf(99, "worker=%08x, task=%d, worker->runningTask=%d\n", 
											 tempw->get_id1(), tempr->number, tempw->runningtask->number);
									if( tempw->runningtask->number == tempr->number ) 
										{
											foundworker = true;
											MWprintf(99, "Found worker!\n");
										}
								}

								// XXX Should be the bug that Steve reported - jichuan
								// running->Append ( tempw );
							workers->Append ( tempw );
						}
			
					if( foundworker ) 
						{
							assert( tempw != NULL );
							MWprintf( 60, "Worker id %08x  has task %d\n", tempw->get_id1(), 
									  tempr->number );
							running->Append ( tempr );
						}
					else
						{
							MWprintf( 10, "Woah!!!  Task number %08x lost his worker.  Rescheduling\n", 
									  tempr->number );
								// First remove him from the run queue
							MWTask *check = rmFromRunQ( tempr->number);
							if( check != tempr ) 
								{
									MWprintf( 10, "There is a HUGE logic error in the code.  Aborting\n" );
									assert( 0 );
								}
							pushTask( tempr );
							break;
						}
				}
		}
}

int 
MWDriver::print_task_keys() 
{
	int retcode = 0;
	
	if( task_key == NULL ) 
		{
			retcode = UNDEFINED;
			MWprintf( 10, " No key function assigned to tasks -- Can't print\n");
		}
	else 
		{
			MWprintf( 10, "Task Keys:\n");
			MWTask *t;

				/* Won't print the keys if todo->number() too large */
			if (todo->number() <= 10000) {
				t = todo->First();
				while (todo->AfterEnd() == false) {
					t = todo->Current();
					MWprintf(10, "   %f\t(%d)\n", (*task_key)(t), t->number);
					todo->Next();
				}
			} else {
				MWprintf(10, "Too many task keys (%d), skip the printing.\n", 
						 todo->number());
			}
		}
	return retcode;
}

int 
MWDriver::delete_tasks_worse_than( MWKey key ) 
{
	if( task_key == NULL ) 
		{
			MWprintf( 10, " delete_task_worse_than:  no task key "
					  "function defined!\n");
			return UNDEFINED;
		}
	
	MWTask *t = todo->First();

	while( todo->AfterEnd() == false) {
		t = todo->Current();
		if( (*task_key)( t ) > key ) {
			todo->RemoveCurrent();
			delete t;
		} else todo->Next();
	}

		//Jeff says -- DO NOT save this state!  What if the user changes the
		//     task key?
		/* Update max_task_key */
#if 0
	if (key < max_task_key)
		max_task_key = key;
#endif
	return 0;
}

void 
MWDriver::putOnRunQ( MWTask *t ) 
{

		// Here's the deal:  running points to the head of the runQ, 
		// runningend points at the last element in the runQ.
		// we insert tasks at the end.
	running->Append ( t );
}

MWTask * 
MWDriver::rmFromRunQ( int jobnum ) 
{
	MWTask *t;

	running->First();
	while (running->AfterEnd() == false) {
		t = running->Current();
		if (t->number == jobnum) {
			running->RemoveCurrent();
			return t;
		} else running->Next();
	}
	return NULL;
}

MWWorkerID * 
MWDriver::task_assigned( MWTask *t ) 
{
	MWWorkerID * w;
	
	workers->ToHead();
	while (workers->AfterEnd() == false) {
		w = workers->Current();
		if ( (w->runningtask) && (w->runningtask->number == t->number) )
			return w;
		else workers->Next();
	}
	return NULL;
}

bool 
MWDriver::task_in_todo_list( MWTask *t )
{
		// Before rescheduling a task on a TASKEXIT message, we need 
		// to see that the task is still in the todo list
		// Tasks NOT in the todo list have already been finished,
		// hence they should not be reassigned and this function should
		// return false.

	MWTask * tt;
	todo->ToHead();

	while (todo->AfterEnd() == false) {
		tt = todo->Current();
		if (tt->number == t->number)
			return true;
		else todo->Next();
	}
	return false;
}  

int
MWDriver::matchTask ( void *arg1, void *arg2 )
{
	MWTask *t1 = (MWTask *)arg1;
	MWTask *t2 = (MWTask *)arg2;
	if ( t1 == t2 )
		return 1;
	return 0;
}

void 
MWDriver::printRunQ() 
{
	MWprintf ( 60, "PrintRunQ start:\n" );
	MWTask *t = running->First();
   	while ( running->AfterEnd() == false )
		{
			t = running->Current();
			t->printself();
			running->Next();
		}
	MWprintf ( 60, "PrintRunQ end.\n\n" );
}

void 
MWDriver::addWorker( MWWorkerID *w ) 
{
		// put a worker on the workers list.  Add to front.
		//  At this point, no benchmark information or machine information
		//  is known, so it is impossible to put the new machine in its proper
		//  place

	w->initGroups ( MWworkClasses );
	MWcurrentWorker = w;
	act_on_starting_worker ( w );
	if ( MWworkClasses <= 1 )
		w->addGroup ( 0 );
	RMC->MW_exec_class_num_workers[w->get_exec_class()]++;

	workers->Append ( w );
}


// This could be more efficient, since we would need only insert items
//  in order to ensure that things are sorted.  But this gives us
//  more flexibility for later, and we will likely rewrite these classes anyway

void 
MWDriver::sort_worker_list()
{
	if( machine_ordering_policy == NO_ORDER || worker_key == NULL  )
		return;

	MWWorkerID *w;
	MWList<MWWorkerID> * newList = new MWList<MWWorkerID>;
	while ( !workers->isEmpty() )
		{
			w = (MWWorkerID *)workers->Remove();
			newList->SortedInsert ( w, (*worker_key)( w ) );
		}
	delete workers;
	workers = newList;
}

MWWorkerID* 
MWDriver::lookupWorker ( int id )
{
	MWWorkerID * w = workers->First();
	while (workers->AfterEnd() == false) {
		w = workers->Current();
		if ( (w->get_id1() == id) || (w->get_id2() == id) )
			return w;
		else workers->Next();
	}
	return NULL;
}

void
MWDriver::call_hostaddlogic() 
{
		/* This is a wrapper around the lower level's hostaddlogic(), 
		   which will return 0 as a basic OK value, and a positive number
		   of hosts to delete if it thinks we should delete hosts.  This
		   can happen if the user changes the target number of workers 
		   to a lower value. */

	int *num_workers = new int[RMC->get_num_exec_classes()];
	int j, numtodelete;
	int total = 0;
	MWWorkerID *w; // , *wx, *wn;
	for ( j=0 ; j < RMC->get_num_exec_classes() ; j++ ) {
		num_workers[j] = numWorkers( j );
	}
	
	numtodelete = RMC->hostaddlogic( num_workers );
	delete [] num_workers;
	
	
		/* Make sure we don't count already doomed machines: */
	if ( numtodelete > 0 ) 
		{
			while ( total < workers->number() )
				{
					total++;
					w = (MWWorkerID *)workers->Remove();
					if ( w->is_doomed() ) numtodelete--;
					workers->Append ( w );
				}
		}

	if ( numtodelete > 0 ) 
		{
			MWprintf ( 40, "Deleting %d hosts.\n", numtodelete );
			MWprintf ( 40, "Ignore the HOSTDELETE or TASKEXIT messages.\n" );

				/* Walk thru list and kill idle workers. */
			total = 0;
			while ( total < workers->number() ) 
				{
					total++;
					w = (MWWorkerID *)workers->Remove();
					if ( numtodelete > 0 && w->isIdle() ) 
						{
							MWprintf ( 80, "Idle worker:\n" );
							worker_last_rites ( w );
							numtodelete--;
						} 
					else 
						{
							workers->Append ( w );
						}
				}
			if ( numtodelete <= 0 ) return;
				 
				/* Now walk thru workers and remove suspended machines */
			total = 0;
			while ( total < workers->number() ) 
				{
					total++;
					w = (MWWorkerID *)workers->Remove();
					if ( numtodelete > 0 && w->isSusp() ) 
						{
							MWprintf ( 80, "Suspended worker:\n" );
							w->printself ( 80 );
				
							hostPostmortem( w ); /* this really does do what we want! */
							RMC->removeWorker( w );
							numtodelete--;
						} 
					else 
						{
							workers->Append ( w );
						}
				}
			if ( numtodelete <= 0 ) return;

				/* At this point, we could do something really cool like 
				   sort the workers based on some attributes.... */

				/* Next, walk thru list and mark working workers for removal */
			total = 0;
			while ( total < workers->number() ) 
				{
					total++;
					w = (MWWorkerID *)workers->Remove();
					if ( numtodelete > 0 && !(w->is_doomed()) ) 
						{
							w->mark_for_removal();
							MWprintf ( 40, "Host %s marked for removal.\n", 
									   w->machine_name() );
							numtodelete--;
						}
					workers->Append ( w );
				}
		}

}

int 
MWDriver::numWorkers() 
{
	return workers->number();
}

int 
MWDriver::numWorkers( int ex_cl ) 
{
	int count = 0;
	int total = 0;
	MWWorkerID *w;
	while ( total < workers->number() ) 
		{
			total++;
			w = (MWWorkerID *)workers->Remove();
			if ( w->get_exec_class() == ex_cl ) 
				{
					count++;
				}
			workers->Append ( w );
		}
	return count;
}


int 
MWDriver::numWorkersInState( int ThisState ) 
{
		// We handle the case WORKING separately such that
		// we don't count the workers with no runningtask

	int total = 0;
	int count = 0;   
	MWWorkerID *w;
  
	while ( total < workers->number() )
		{
			total++;
			w = (MWWorkerID *)workers->Remove();
			if (ThisState == WORKING)
				{
					if ( (w->currentState() == ThisState) && (w->runningtask != NULL) ) 
						{
							count++;
						}
				}
			else
				{
					if ( w->currentState() == ThisState ) 
						{
							count++;
						}
				}
			workers->Append ( w );
		}
	return count;
}

void 
MWDriver::printWorkers() 
{
	MWprintf ( 10, "---- A list of Workers follows: ----------------\n" );

	MWWorkerID *w = workers->First();

	while ( workers->AfterEnd() == false ) 
		{
			w = workers->Current();
			w->printself( 10 );
			MWprintf( 10, "\n" );
			workers->Next();
		}

	MWprintf ( 10, "---- End worker list --- %d workers ------------\n", workers->number() );
}

MWKey 
MWDriver::return_best_todo_keyval()
{
	MWKey retval = DBL_MAX;

	if( task_key == NULL ) {
		MWprintf( 10, " return_bst_keyval:  no task key "
				  "function defined!\n");      
	}
	else {
		if ( todo->isEmpty() )
			return retval;
		if ( listsorted ) {
				//XXX Jeff -- fixed the bug, if you take something off the list that has a key, 
				//  You'd better damn well put it back with the same key (not 0.0)
				//  but I don't know if it is the best way to remove
				//  and put the key back on the list.  I know I didn' write it this way
			MWTask *t = (MWTask *)todo->Remove();
			double val = (*task_key)(t);
			todo->Prepend(t, val);
			return (*task_key)( t );
		}
		else {
			MWTask * t = todo->First();
			while ( todo->AfterEnd() == false ) {
				t = todo->Current();
				MWKey tmp = (*task_key) ( t );
				if ( retval > tmp )
					retval = tmp;
				todo->Next();
			}
			return retval;
		}
	}
	return retval;
}

MWKey MWDriver::return_best_running_keyval()
{
	MWKey retval = DBL_MAX;

	if( task_key == NULL )
		{
			MWprintf( 10, " return_bst_keyval:  no task key "
					  "function defined!\n");      
    	}
	else
		{
			if ( running->isEmpty() ) return retval;
			MWTask * t = running->First();
			while ( running->AfterEnd() == false )
				{
					t = running->Current();
					if ( t->taskType == MWNORMAL )
						{
							MWKey tmp = (*task_key)( t );
							if( tmp < retval )
								retval = tmp;
						}
					running->Next();
				}
			return retval;
		}
	return retval;
}



MWWorkerID * 
MWDriver::rmWorker ( int id ) 
{
        // search for a worker with the given id and remove it.
		// The id can refer to either the primary or secondary id.
	
	MWWorkerID * w = workers->First();
	while (workers->AfterEnd() == false) {
		w = workers->Current();
		if ( (w->get_id1() == id) || (w->get_id2() == id) ) {
			workers->RemoveCurrent();
			return w;
		}
		else workers->Next();
	}
	return NULL;
}

int
MWDriver::set_checkpoint_frequency( int freq ) {
	if ( checkpoint_time_freq != 0 ) {
		MWprintf ( 10, "Cannot set_checkpoint_frequency while time-based "
				   "frequency is not zero!\n" );
		return checkpoint_frequency;
	}
    int old = checkpoint_frequency;
    checkpoint_frequency = freq;
    return old;
}

int 
MWDriver::set_checkpoint_time ( int secs ) {
	if ( checkpoint_frequency != 0 ) {
		MWprintf ( 10, "Cannot set_checkpoint_time while task-based "
				   "frequency is not zero!\n" );
		return checkpoint_time_freq;
	}
	int old = checkpoint_time_freq;
	checkpoint_time_freq = secs;
	next_ckpt_time = time(0) + secs;
	return old;
}

void 
MWDriver::set_suspension_policy( MWSuspensionPolicy policy ) {
	if ( (policy!=DEFAULT) && (policy!=REASSIGN) ) {
		MWprintf ( 10, "Bad suspension policy %d.  Using DEFAULT.\n" );
		suspensionPolicy = DEFAULT;
		return;
	}
	suspensionPolicy = policy;
}

void
MWDriver::checkpoint() 
{
    FILE *cfp;  // Checkpoint File Pointer
		/* We're going to write the checkpoint to a temporary file, 
		   then move the file to "ckpt_filename" */

		/* We're not going to use tempnam(), because it sometimes 
		   gives us a filename in /var/tmp, which is bad for rename() */

	char *tempName = "mw_tmp_ckpt_file";

    if ( ( cfp=fopen( tempName, "w" ) ) == NULL ) {
        MWprintf ( 10, "Failed to open %s for writing! errno %d\n", 
				   tempName, errno );
        return;
    }

    MWprintf ( 50, "Beginning checkpoint()\n" );

        // header
    fprintf ( cfp, "MWDriver Checkpoint File.  %ld\n", time(0) );
    
        // some internal MWDriver state:
    fprintf ( cfp, "%d %d %d %d\n%d %d %d %d %d %lf %d\n", task_counter, 
			  checkpoint_frequency, checkpoint_time_freq, 
			  num_completed_tasks,
			  (int)addmode, (int)getmode, (int)suspensionPolicy,
			  (int) machine_ordering_policy,
	          worker_timeout,worker_timeout_limit, worker_timeout_check_frequency);

	if ( bench_task ) {
		fprintf ( cfp, "1 " );
		bench_task->write_ckpt_info ( cfp );
	} else {
		fprintf ( cfp, "0\n" );
	}

	MWprintf ( 80, "Wrote internal MWDriver state.\n" );

	RMC->write_checkpoint( cfp );

	MWprintf ( 80, "Wrote RMC state.\n" );

		/* Yes, I really am passing it the list of workers... */
	stats->write_checkpoint ( cfp, workers );
	
	MWprintf ( 80, "Wrote the Worker stats.\n" );

		// Write the user master's state:
	write_master_state( cfp );

	MWprintf ( 80, "Wrote the user master's state.\n" );

		// write the number of work classes
	fprintf ( cfp, "%d\n", MWworkClasses );

	int tempnum = 0;
	MWTask * ttt = running->First();
	while ( running->AfterEnd() == false )
		{
			ttt = running->Current();
			if ( ttt->taskType == MWNORMAL )
				tempnum++;
			running->Next();
		}
	int num_tasks = tempnum + todo->number();
    
		// Tasks separator:
	fprintf ( cfp, "Tasks: %d \n", num_tasks );

		/* Write ordered todo before unsorted running tasks */
	MWTask *t;
	
	fprintf ( cfp, "TODO_num: %d \n", todo->number() );
	t = todo->First();
	while ( todo->AfterEnd() == false ) 
		{
			t = todo->Current();
			fprintf ( cfp, "%d ", t->number );
			t->write_ckpt_info( cfp );
			t->write_group_info ( cfp );
			todo->Next();
		}

	fprintf( cfp, "running_num: %d \n", running->number() ); 
	t = running->First();
	while ( running->AfterEnd() == false) 
		{
			t = running->Current();
			if ( t->taskType == MWNORMAL )
				{
					fprintf ( cfp, "%d ", t->number );
					t->write_ckpt_info( cfp );
					t->write_group_info ( cfp );
				}
			running->Next();
		}

	MWprintf ( 80, "Wrote task list.\n" );

	fclose ( cfp );

		/* Now we rename our temp file to the real thing. Atomicity! */
	if ( rename( tempName, ckpt_filename ) == -1 ) {
		MWprintf ( 10, "rename( %s, %s ) failed! errno %d.\n", 
				   tempName, ckpt_filename, errno );
	}

	MWprintf ( 50, "Done checkpointing.\n" );
}

void 
MWDriver::restart_from_ckpt() 
{

	int i;
	
	FILE *cfp;
	if ( ( cfp=fopen( ckpt_filename, "r" ) ) == NULL ) {
		MWprintf ( 10, "Failed to open %s for reading! errno %d\n", 
				   ckpt_filename, errno );
		return;	
	}	

	char buf[128];
	time_t then;
	fscanf( cfp, "%s %s %s %ld", buf, buf, buf, &then );
    
	MWprintf ( 10, "This checkpoint made on %s\n", ctime( &then ) );

	fscanf( cfp, "%d %d %d %d", &task_counter, 
			&checkpoint_frequency, &checkpoint_time_freq,
			&num_completed_tasks );

	fscanf( cfp, "%d %d %d %d", 
			(int*) &addmode, (int*) &getmode, (int*) &suspensionPolicy,
			(int*) &machine_ordering_policy );

	int worker_timeout_int;
	fscanf( cfp, "%d %lf %d", &worker_timeout_int, 
		&worker_timeout_limit, &worker_timeout_check_frequency);
	worker_timeout_int ? worker_timeout = true : worker_timeout = false;

	int hasbench = 0;
	fscanf( cfp, "%d ", &hasbench );
	if ( hasbench ) {
		bench_task = gimme_a_task();
		bench_task->read_ckpt_info( cfp );
		MWprintf ( 10, "Got benchmark task:\n" );
		bench_task->printself( 10 );
	} 

	MWprintf ( 10, "Read internal MW info.\n" );

	RMC->read_checkpoint( cfp );

	MWprintf ( 10, "Read RMC state.\n" );

		// read old stats 
	stats->read_checkpoint( cfp );
	
	MWprintf ( 10, "Read defunct workers' stats.\n" );

	read_master_state( cfp );

	MWprintf ( 10, "Read the user master's state.\n" );

	fscanf ( cfp, "%d", &MWworkClasses );
	MWworkClassWorkers = new int[MWworkClasses];
	MWworkClassTasks = new int[MWworkClasses];
	MWcurrentWorker = NULL;
	for ( i = 0; i < MWworkClasses; i++ )
		{
			MWworkClassWorkers[i] = 0;
			MWworkClassTasks[i] = 0;
		}

	int num_tasks;
	fscanf ( cfp, "%s %d", buf, &num_tasks );

	if ( strcmp( buf, "Tasks:" ) ) {
		MWprintf ( 10, "Problem in reading Tasks separator. buf = %s\n", buf );
		fclose ( cfp );
		return;
	}

	MWTask *t;
	
	int num_todo; 
	fscanf( cfp, "%s %d", buf, &num_todo);
	if ( strcmp( buf, "TODO_num:" ) ) {
		MWprintf ( 10, "Problem in reading TODO_num separator. buf = %s\n", buf );
		fclose( cfp );
		return;
	}

	if (task_key != NULL) {
		/* Add all tasks by key, if user has set it. */
		set_task_add_mode( ADD_BY_KEY );
	}
	else {
		set_task_add_mode( ADD_AT_END );
	}
	
	for ( i=0 ; i<num_todo ; i++ ) 	{
		t = gimme_a_task();   // return a derived task class...
		fscanf( cfp, "%d ", &(t->number) );
		t->read_ckpt_info( cfp );
		t->read_group_info ( cfp );
		t->printself( 90 );
		pushTask(t);
	}

	int num_running;
	fscanf( cfp, "%s %d", buf, &num_running);
	if ( strcmp( buf, "running_num:" ) ) {
		MWprintf ( 10, "Problem in reading running_num separator. buf = %s\n", buf );
		fclose( cfp );
		return;
	}

	for ( i=0 ; i<num_running ; i++ ) 
		{
			t = gimme_a_task();   // return a derived task class...
			fscanf( cfp, "%d ", &(t->number) );
			t->read_ckpt_info( cfp );
			t->read_group_info ( cfp );
			t->printself( 90 );
			pushTask(t);
		}

	MWprintf ( 10, "Read the task list.\n" );

		// This is debugging code for Michael

#if 1
	if( task_key != NULL ) {
		MWprintf( 10, "Sort the TODO task list...\n" );
		sort_task_list();
		MWprintf( 10, "Printing keys of user's task list...\n" );
		print_task_keys();
	}
#endif
	fclose( cfp );
    
	FILE *swapfile;
	if ( (swapfile = fopen(todo->Name(), "r")) == NULL ) {
		MWprintf(90, "Task_Swap:: No task swap file exists, has_task_swapped = false.\n");
		has_task_swapped = false;
	} else {
		fclose(swapfile);
	}
	
		/* Task_Swap testing code */
	if (todo->Number() > MAX_IN_MEM) {
		swap_out_todo_tasks();
		checkpoint();
	}
	
	MWprintf ( 10, "We are done restarting.\n" );
}

/* Read/write task */
MWTask*
MWDriver::read_mem_task(MWList<MWTask> *tasks)
{
	if (tasks->AfterEnd())
		return NULL;

	MWTask * t = tasks->RemoveCurrent();
	todo->Next();
	return t;
}

MWTask*
MWDriver::read_file_task(FILE *f)
{
	char strbuf[128];

	fscanf(f, "%s ", strbuf);
	if (strcmp(strbuf, "~~Task~~")) {
		MWprintf(10, "Problem in reading ~~Task~~ separator, buf=%s\n", strbuf);
		return NULL;
	}

	MWTask *t = gimme_a_task();
	fscanf(f, "%d ", &(t->number));
	t->read_ckpt_info(f);
	t->read_group_info(f);
	MWprintf(90, "Read task %d\n", t->number);
	return t;
}

void
MWDriver::write_task(FILE *f, MWTask *t)
{
	fprintf(f, "%s \n", "~~Task~~");
	fprintf(f, "%d \n", t->number);
	t->write_ckpt_info(f);
	t->write_group_info(f);
	MWprintf(90, "Written task %d\n", t->number);
}

/* Remember to update (1) how to check whether we are done;
 * (2) give each list a "name", so that they can be dumped into 
 * a file; (3) change remove_tasks_worse_than() to update the 
 * max_task_key value; (4) add logic to trigger the swapping */

/** TODO Tasks Swapping **/
bool 
MWDriver::swap_out_todo_tasks(int num_in_mem, int max_in_mem, double max_key)
{
	char strbuf[128];
	MWTask *mem_t = NULL, *file_t = NULL;
	double mem_key = -DBL_MAX, file_key = -DBL_MAX;
	FILE *swap, *new_swap;
	int num_swapped, num_to_skip, num_file_task=0, num_total;
	
	if ( (task_key == NULL)||(todo->isSorted()==false) ) {
		MWprintf(10, "Task_Swap:: No task_key func defined for the list.\n");
		MWprintf(10, "Task_Swap:: I won't swap these tasks are not sorted.\n");
		return false;
	}

	if (todo->number() < max_in_mem)
		return true;
	if (todo->number() < num_in_mem) /* just in case :=) */
		return true;
	
		/* Skip the first num_in_mem tasks */
	mem_t = todo->First();
	for (int i=0; i<num_in_mem; i++)
		todo->Next();

	FILE *swapfile;
	if ( (swapfile = fopen(todo->Name(), "r")) == NULL) {
		MWprintf(10, "Task_Swap:: No existing swap file, dump memory directly into the file.\n");
		
			/* Create the swap file for the first time */
		swap = Open(todo->Name(), "a");
		if (swap == NULL) {
			MWprintf(10, "Task_Swap:: Can't create swap file.\n");
			return false;
		}
		
			/* As will rebuild index, we ClearIndex() first */
		todo->ClearIndex();
		
			/* First write the num_to_skip to 0 */
		fprintf(swap, "num_to_skip %10d\n", 0);
		fprintf(swap, "num_of_file_task %10d\n", todo->number()-num_in_mem);
		
			/* Remove tasks until we only have num_in_mem in list */
		num_swapped = 0;
		while ( (mem_t=read_mem_task(todo)) != NULL ) {
			write_task(swap, mem_t);
			delete mem_t;
			num_swapped ++;
		}
		fflush(swap);
		Close(swap);
		
			/* rebuild index for todo */
		todo->BuildIndex();

			/* Remember that we have tasks swapped out */
		has_task_swapped = true;

		MWprintf(10, "Task_Swap:: Done swapping, swapped out %d tasks, keep %d tasks in memory.\n",
				 num_swapped, todo->number());
	} else {
		
		MWprintf(10, "Task_Swap:: Will merge tasks in memory with those in file.\n");

		fclose(swapfile);
			/* As will rebuild index, we ClearIndex() first */
		todo->ClearIndex();
		
			/* Create new temp file */
		new_swap = Open("tmp_TODO_swap", "w");
		if (new_swap == NULL) {
			MWprintf(10, "Task_Swap:: Can't create tmp swap file.\n");
			return false;
		}

		fprintf(new_swap, "num_to_skip %10d\n", 0);
		fprintf(new_swap, "num_of_file_task %10d\n", 0);
		
			/* Open the swap file, read num_to_skip */ 
		swap = Open(todo->Name(), "r");
		if (swap == NULL) {
			MWprintf(10, "Task_Swap:: Can't open swap file to read.\n");
			return false;
		}
		
		fscanf(swap, "%s %10d", strbuf, &num_to_skip);
		if ( strcmp(strbuf, "num_to_skip") ) {
			MWprintf(10, "Task_Swap:: problem in reading num_to_skip marker.\n");
			Close(swap);
			return false;
		}

			/* Read num_of_file_task */
		fscanf(swap, "%s %10d", strbuf, &num_file_task);
		if ( strcmp(strbuf, "num_of_file_task") ) {
			MWprintf(10, "Task_Swap:: problem in reading num_of_file_task marker.\n");
			Close(swap);
			return false;
		}

			/* Skip tasks at the beginning */
		for (int i=0; i<num_to_skip; i++) {
			file_t = read_file_task(swap);
			if (file_t != NULL)
				delete file_t;
		}

		if ( (mem_t = read_mem_task(todo)) != NULL ) {
			mem_key = (*task_key)(mem_t);
		}

		if ( (file_t = read_file_task(swap)) != NULL ) {
			file_key = (*task_key)(file_t);
			num_file_task --;
		}
		
			/* Begin merging the rest of the file and the rest of the list */
		num_total = 0;
		while ( mem_t && file_t ) /* while both valid */ {
			if (mem_key <= file_key) { 
					/* dump current mem_t */
				if (mem_key < max_task_key) {
					write_task(new_swap, mem_t);
					num_total ++;
				}
				delete mem_t;
				if ( (mem_t = read_mem_task(todo)) != NULL ) {
					mem_key = (*task_key)(mem_t);
				}
			} else { 
					/* dump current file_t */
				if (file_key < max_task_key) {
					write_task(new_swap, file_t);
					num_total ++;
				}
				delete file_t;
				
				if (num_file_task > 0) {
					if ( (file_t = read_file_task(swap) ) != NULL ) {
						file_key = (*task_key)(file_t);
						num_file_task --;
					} 
				} else file_t = NULL;
			}
		}

			/* dump the remaining */
		if ((mem_t != NULL)) {
			write_task(new_swap, mem_t);
			num_total ++;
			delete mem_t;
			
			while ( (mem_t = read_mem_task(todo)) != NULL ) {
				mem_key = (*task_key)(mem_t);
				if (mem_key < max_task_key) {
					write_task(new_swap, mem_t);
					num_total ++;
				}
				delete mem_t;
			}
		} else if (file_t != NULL) {
			write_task(new_swap, file_t);
			num_total ++;
			delete file_t;
			
			while ( num_file_task > 0) {
				if ( (file_t = read_file_task(swap)) != NULL ) {
					file_key = (*task_key)(file_t);
					if (file_key < max_task_key) {
						write_task(new_swap, file_t);
						num_total ++;
					}
					delete file_t;
				}
				num_file_task --;
			}
		} else {
			MWprintf(10, "Task_Swap:: Write task error, either mem_t or file_t should be NULL.\n");
		}
		
		MWprintf(10, " num_file_task now is %d\n", num_file_task);
		Close(swap);
			
		fflush(new_swap);
		Close(new_swap);

			/* Update header info */
		new_swap = Open("tmp_TODO_swap", "r+");
		if (new_swap == NULL) {
			MWprintf(10, "Task_Swap:: Can't open swap file to read.\n");
			return false;
		}
		
		fprintf(new_swap, "num_to_skip %10d\n", 0);
		fprintf(new_swap, "num_of_file_task %10d\n", num_total);
		fflush(new_swap);
		Close(new_swap);

			/* Swith the files */
		char cmd[1024];
		sprintf(cmd, "mv %s old_TODO_tasks", todo->Name());
		if (system(cmd) < 0) {
			MWprintf(10, "Task_Swap:: Can't rename the old swap file, errno = %d.\n", errno);
			return false;
		}

		sprintf(cmd, "mv -f tmp_TODO_swap %s", todo->Name());
		if (system(cmd) < 0) {
			MWprintf(10, "Task_Swap:: Can't rename the new swap file, errno = %d.\n", errno);
			return false;
		}

		remove("old_TODO_tasks");

			/* rebuild index for todo */
		todo->BuildIndex();
		
			/* Remember that we have tasks swapped out */
		has_task_swapped = true;

			/* Print some trash */
		MWprintf(10, "Task_Swap:: Done swapping, swapped out %d tasks, keep %d tasks in memory.\n",
				 num_total, todo->number());
	}

		/* Keep a snapshot */
	checkpoint();
	return true;
}

bool 
MWDriver::swap_in_todo_tasks(int min_in_mem, int num_in_mem)
{
	char strbuf[128];
	MWTask *file_t = NULL, * mem_t = NULL;
	MWList<MWTask> * todo_new, * todo_old;
	FILE *swap;
	int num_swapped, num_to_skip, num_file_task = 0;
	
	if ( (todo->number() > num_in_mem) )
		return true;
	
	FILE *swapfp;
	if ( (swapfp = fopen(todo->Name(), "r")) == NULL) {
		MWprintf(90, "Task_Swap:: No task swap file exists, will not swap in.\n");
		has_task_swapped = false;
		return true;
	}
	
	fclose(swapfp);
		/* Now need to bring more task on disk into mem */
		/* Open the swap file, read num_to_skip */ 
	swap = Open(todo->Name(), "r+");
	if (swap == NULL) {
		MWprintf(10, "Task_Swap:: Can't open swap file to read.\n");
		return false;
	}

	fscanf(swap, "%s %10d", strbuf, &num_to_skip);
	if ( strcmp(strbuf, "num_to_skip") ) {
		MWprintf(10, "Task_Swap:: problem in reading num_to_skip marker.\n");
		Close(swap);
		return false;
	}

		/* Read num_of_file_task */
	fscanf(swap, "%s %10d", strbuf, &num_file_task);
	if ( strcmp(strbuf, "num_of_file_task") ) {
		MWprintf(10, "Task_Swap:: problem in reading num_of_file_task marker.\n");
		Close(swap);
		return false;
	}

		/* Skip tasks at the beginning */
	file_t = gimme_a_task();
	if (file_t == NULL) {
		MWprintf(10, "Task_Swap:: Can't allocate a task.\n");
		return false;
	}

		/* Skip some tasks */
	for (int i=0; i<num_to_skip; i++) {
		file_t = read_file_task(swap);
		if (file_t != NULL)
			delete file_t;
	}

		/* Now begin read into mem */
	todo_new = new MWList<MWTask>("TODO_tasks");
	if (todo_new == NULL) {
		MWprintf(10, "Task_Swap:: Can't allocate a new task list.\n");
		Close(swap);
		return false;
	}
	
	num_swapped = (num_in_mem < num_file_task) ? num_in_mem : num_file_task;
	for (int ii=0; ii<num_swapped; ii++) 
		if ( (file_t = read_file_task(swap)) != NULL )
			todo_new->SortedInsert(file_t, (*task_key)(file_t));
	
		/* Update header info */
	if (num_swapped == num_file_task) { /* No task on disk then */
		Close(swap);
		remove(todo->Name());
		has_task_swapped = false;
	} else {
		if (fseek(swap, 0L, SEEK_SET) == -1) {
			MWprintf(10, "Task_Swap:: Can update header_info in the swap file, errno = %d.\n", errno);
			return false;
		}

		fprintf(swap, "num_to_skip %10d\n", num_to_skip + num_swapped);
		fprintf(swap, "num_of_file_task %10d\n", num_file_task-num_swapped);
		fflush(swap);

			/* Close new swap file and switch */
		Close(swap);
	}

		/* Sorted insert old tasks into new todo list */
	todo->First();
	while ( (mem_t = read_mem_task(todo))!=NULL ) {
		todo_new->SortedInsert(mem_t, (*task_key)(mem_t));
	}

	todo_old = todo;
	delete todo_old;
	todo = todo_new;
	
	MWprintf(10, "Task_Swap:: done swap in, moved %d tasks into memory. num_todo = %d.\n", 
			 num_swapped, todo->Number());

		/* Keep a snapshot */
	checkpoint();
	return true;
}

bool 
MWDriver::is_TODO_empty()
{
	if (todo->number() > 0)
		return false;
	
	if (has_task_swapped)
		return false;

	return true;
}
	
int	/* Added as in $GAMS */
MWDriver::sort_task_list ( )
{
		/* check and only sort the list when listsorted is FALSE */
	if ( (task_key == NULL) || (listsorted) )
		return 0;
	
	MWList<MWTask> * newList = new MWList<MWTask>;
	while ( todo->number() > 0 )
        {
			MWTask *t = (MWTask *)todo->Remove();
			newList->SortedInsert ( t, (*task_key)(t) );
        }
	delete todo; 
        
	todo = newList; 
	
		/* reset listsorted as TRUE */
	listsorted = true;
	return 0;
}


int 
MWDriver::get_number_tasks() 
{
	return todo->number();
}

int MWDriver::get_number_running_tasks()
{
	return running->number();
}

void
MWDriver::workClasses_set ( int num )
{
	if ( num < 1 )
		{
			MWprintf ( 10, "Use of workClasses_set with <= 1 classes not allowed\n");
			assert ( 0 );
		}
	MWworkClasses = num;
	MWworkClassWorkers = new int[MWworkClasses];
	MWworkClassTasks = new int[MWworkClasses];

	for ( int i = 0; i < MWworkClasses; i++ )
		{
			MWworkClassWorkers[i] = 0;
			MWworkClassTasks[i] = 0;
		}
}

int
MWDriver::workClasses_get ( )
{
	return MWworkClasses;
}

int
MWDriver::workClasses_getworkers ( int num )
{
	return MWworkClassWorkers[num];
}

int
MWDriver::workClasses_gettasks ( int num )
{
	return MWworkClassTasks[num];
}


/* 8/28/00

Qun added Jeff's idea of adding a SORTED list of tasks.
This can greatly improve the efficiency in many applications.
   
XXX I haven't debugged this, and probably there needs to be more
checking that the list is actually sorted.  
*/

/* added by Qun: for insert a list of sorted tasks */

void 
MWDriver::addSortedTasks( int n, MWTask **add_tasks )
{
	if ( n <= 0 ) {
		MWprintf ( 10, "Please add a positive number of tasks!\n" );
		RMC->exit(1);
	}

	for ( int i = 0; i < n; i++ )
		{
			add_tasks[i]->number = task_counter++;
			if ( MWworkClasses <= 1 )
				{
					add_tasks[i]->initGroups ( MWworkClasses );
					add_tasks[i]->addGroup ( 0 );
				}
			pushTask ( add_tasks[i] );
		}
}


void 
MWDriver::addTaskByKey( MWTask *add_task )  
{
	add_task->number = task_counter++;
	if ( MWworkClasses <= 1 )
		{
			add_task->initGroups ( MWworkClasses );
			add_task->addGroup ( 0 );
		}
	if ( task_key == NULL )
		{
			MWprintf ( 10, "ERROR!  Adding by key, but no key "
					   "function defined!\n" );
			assert(0);
		}
	todo->SortedInsert ( add_task, (*task_key)(add_task) );
}

//end the changes done  by Qun

extern MWWorker *gimme_a_worker ();

void
MWDriver::ControlPanel ( )
{

#ifdef INDEPENDENT
    MWReturn ustat = OK;
    RMC->hostadd ();
    master_mainloop_oneshot ( 0, 2 );
    MWWorker *worker = gimme_a_worker();
    worker->go( 1, NULL );
    master_mainloop_oneshot ( 0, 2 );
    worker->do_benchmark_cmd( );
    ustat = master_mainloop_oneshot ( 0, 2 );
    
		// while ( ( (todo != NULL) || (running != NULL) ) && ( ustat == OK ) )
    while ( ( (!todo->isEmpty()) || (!running->isEmpty()) ) && ( ustat == OK ) )
		{
			ustat = worker->worker_mainloop_oneshot ( );
			ustat = master_mainloop_oneshot ( 0, 2 );
		}
    
    MWprintf(71, "\n\n\n ***** Almost done ***** \n\n\n\n");

		// found LEAK_SCOPE! need to deallocate worker
    delete worker;
#endif
    return;
}


MWKey
kflops( MWWorkerID *w )
{
	return (MWKey) w->KFlops;
}

MWKey
benchmark_result( MWWorkerID *w )
{
	return (MWKey) w->get_bench_result();
}

// These are functions that might be useful, but I am not sure if they should
//  be included in the final release.

double 
MWDriver::get_instant_pool_perf( ) 
{
	double total_perf = 0;
	int total = 0;
	MWWorkerID *w;

	while ( total < workers->number() ) 
		{
			total++;
			w = (MWWorkerID *)workers->Remove();
			if ( (w->currentState() == WORKING) && (w->runningtask != NULL) ) 
				{
					total_perf += w->get_bench_result();
				}
			workers->Append ( w );
		}
	return total_perf;
}

MWTask*
MWDriver::get_todo_head ( )
{
	MWTask *t = (MWTask *)todo->Remove();
	todo->Prepend ( t );
	return t;
}

MWWorkerID*
MWDriver::get_workers_head ( )
{
	MWWorkerID *w = (MWWorkerID *)workers->Remove();
	workers->Prepend ( w );
	return w;
}


#if defined( XML_OUTPUT )

void MWDriver::write_XML_status()
{
	ofstream xmlfile("/u/m/e/metaneos/public/html/iMW/status.xml",
					 ios::trunc|ios::out);


	if( ! xmlfile ) {
		cerr << "Cannot open 'status.xml data file' file!\n";
	}

		//system("/bin/rm -f ~/public/html/iMW/status.xml");
		//system("/bin/mv -f ./status.xml ~/public/html/iMW/");   

		// system("/bin/more status.xml");

	char *temp_ptr ;
	temp_ptr = get_XML_status();

	xmlfile << temp_ptr;
	delete temp_ptr;

		//  xmlfile << ends;

	xmlfile.close();

		//system("/bin/rm -f ~/public/html/iMW/status.xml");
		//system("/bin/cp -f status.xml /u/m/e/metaneos/public/html/iMW/");
		// system("/bin/cp -f status.xml ~metaneos/public/html/iMW/");
}

#include <strstream.h>

char* MWDriver::get_XML_status(){

	ifstream menus(xml_menus_filename);

	ostrstream xmlstr;

	xmlstr << "<?xml version=\"1.0\"?>" << endl;
	xmlstr << "<?xml:stylesheet type=\"text/xsl\" href=\"menus-QAP.xsl\"?>" << endl;
	xmlstr << "<INTERFACE TYPE=\"iMW-QAP\">" << endl;

		// xmlstr << menus;

	char ch;
	while (menus.get(ch)) xmlstr.put(ch);

		// menus >> xmlstr ;

	xmlstr << "<MWOutput TYPE=\"QAP\">" << endl;

	char *temp_ptr ;

	temp_ptr = get_XML_job_information();

	xmlstr << temp_ptr;
	delete temp_ptr;

	temp_ptr = get_XML_problem_description();

	xmlstr << temp_ptr;
	delete temp_ptr;

	temp_ptr = get_XML_interface_remote_files();

	xmlstr << temp_ptr;
	delete temp_ptr;

	temp_ptr = get_XML_resources_status();

	xmlstr << temp_ptr;
	delete temp_ptr;

	temp_ptr = get_XML_results_status();

	xmlstr << temp_ptr;
	delete temp_ptr;

	xmlstr << "</MWOutput>" << endl;
	xmlstr << "</INTERFACE>" << endl;

	xmlstr << ends;



	menus.close();


	return xmlstr.str();

}

char* MWDriver::get_XML_job_information()
{

	ifstream jobinfo(xml_jobinfo_filename);

	ostrstream xmlstr;

	char ch;
	while (jobinfo.get(ch)) xmlstr.put(ch);

	xmlstr << ends;

	jobinfo.close();
 

	return xmlstr.str();

};

char* MWDriver::get_XML_problem_description(){

	ifstream pbdescrib(xml_pbdescrib_filename);

	ostrstream xmlstr;

	char ch;
	while (pbdescrib.get(ch)) xmlstr.put(ch);

	xmlstr << ends;

	pbdescrib.close();


	return xmlstr.str();

};

char* MWDriver::get_XML_interface_remote_files(){

	ostrstream xmlstr;

	xmlstr << "<InterfaceRemoteFiles>" << endl;

		// dump here content of file

	xmlstr << "</InterfaceRemoteFiles>" << endl;

	xmlstr << ends;

	return xmlstr.str();

}


char* MWDriver::get_XML_resources_status()
{

	int total;
	int i;
	ostrstream xmlstr;

		// Begin XML string

	xmlstr << "<ResourcesStatus>" << endl;

	double average_bench;
	double equivalent_bench;
	double min_bench;
	double max_bench;
	double av_present_workers;
	double av_nonsusp_workers;
	double av_active_workers;
	double equi_pool_performance;
	double equi_run_time;
	double parallel_performance;
	double wall_time;

	stats->get_stats(&average_bench,
					 &equivalent_bench,
					 &min_bench,
					 &max_bench,
					 &av_present_workers,
					 &av_nonsusp_workers,
					 &av_active_workers,
					 &equi_pool_performance,
					 &equi_run_time,
					 &parallel_performance,
					 &wall_time,
					 workers
					 );

		// MWStats Information

	xmlstr << "<MWStats>" << endl;

	xmlstr << "<WallClockTime>" << wall_time << "</WallClockTime>" << endl;
	xmlstr << "<NormalizedTotalCPU>" << equi_run_time << "</NormalizedTotalCPU>" << endl;
	xmlstr << "<ParallelEff>" <<  parallel_performance << "</ParallelEff>" << endl;

	xmlstr << "<BenchInfo>" << endl;
	xmlstr << "<InstantPoolPerf>" << get_instant_pool_perf() << "</InstantPoolPerf>" << endl;
	xmlstr << "<EquivalentPoolPerf>" << equi_pool_performance << "</EquivalentPoolPerf>" << endl;
	xmlstr << "<AverageBench>" << average_bench << "</AverageBench>" << endl;
	xmlstr << "<EquivalentBench>" << equivalent_bench << "</EquivalentBench>" << endl;
	xmlstr << "<MinBench>" << min_bench << "</MinBench>" << endl;
	xmlstr << "<MaxBench>" << max_bench << "</MaxBench>" << endl;
	xmlstr << "</BenchInfo>" << endl;

	xmlstr << "<AverageWorkersStats>" << endl;
	xmlstr << "<AveragePresentWorkers>" << av_present_workers << "</AveragePresentWorkers>" << endl;
	xmlstr << "<AverageNonSuspWorkers>" << av_nonsusp_workers << "</AverageNonSuspWorkers>" << endl;
	xmlstr << "<AverageActiveWorkers>"  << av_nonsusp_workers << "</AverageActiveWorkers>" << endl;
	xmlstr << "</AverageWorkersStats>" << endl;


	xmlstr << "</MWStats>" << endl;

		// Master Information

	xmlstr << "<Master>" << endl;

	xmlstr << "<MasterPhysicalProperties>" << endl;
	xmlstr << "<Name>" << mach_name << "</Name>" << endl;
		//     xmlstr << "<IPAddress>" << get_IPAddress() << "</IPAddress>" << endl;
	xmlstr << "<OpSys>" << get_OpSys() << "</OpSys>" << endl;
	xmlstr << "<Arch>" << get_Arch() <<"</Arch>" << endl;
	xmlstr << "<Memory>" << get_Memory() << "</Memory>" << endl;
	xmlstr << "<VirtualMemory>" << get_VirtualMemory() << "</VirtualMemory>" <<  endl; 
	xmlstr << "<DiskSpace>" << get_Disk() << "</DiskSpace>" << endl;
	xmlstr << "<KFlops>" << get_KFlops() << "</KFlops>" << endl;
	xmlstr << "<Mips>" << get_Mips() << "</Mips>" << endl;
	xmlstr << "<CPUs>" << get_Cpus() << "</CPUs>" << endl;
	xmlstr << "<NWSinfos/>" << endl;
	xmlstr << "</MasterPhysicalProperties>" << endl;
                       


	xmlstr << "<MasterUsageProperties>" << endl;
	xmlstr << "<StartTime>July 1st 1999, 18:21:12s GMT</StartTime>"  << endl;
	xmlstr << "</MasterUsageProperties>" << endl;

	xmlstr << "</Master>" << endl;

		// Worker Information

	xmlstr << "<Workers>" << endl;

	int numworkers = numWorkersInState( INITIALIZING ) +  numWorkersInState( BENCHMARKING ) + numWorkersInState( IDLE ) + numWorkersInState( WORKING ) ;

	xmlstr << "<WorkersNumber>" << numworkers << "</WorkersNumber>" << endl;

	xmlstr << "<WorkersStats>" << endl;

	xmlstr << "<WorkersInitializing>" << numWorkersInState( INITIALIZING ) << "</WorkersInitializing>" << endl;
	xmlstr << "<WorkersBenchMarking>" << numWorkersInState( BENCHMARKING ) << "</WorkersBenchMarking>" << endl;
	xmlstr << "<WorkersWaiting>" << numWorkersInState( IDLE ) << "</WorkersWaiting>" << endl;
	xmlstr << "<WorkersWorking>" << numWorkersInState( WORKING ) << "</WorkersWorking>" << endl;
	xmlstr << "<WorkersSuspended>" << numWorkersInState( SUSPENDED ) << "</WorkersSuspended>" << endl;
	xmlstr << "<WorkersDone>" << numWorkersInState( EXITED ) << "</WorkersDone>" << endl;

	xmlstr << "</WorkersStats>" << endl;

	xmlstr << "<WorkersList>" << endl;

	total = 0;
	MWWorkerID *w;

	while ( total < workers->number() ) 
		{
			total++;
			w = (MWWorkerID *)workers->Remove();

			if ((w->currentState() != SUSPENDED) || (w->currentState() != EXITED))
				{
		         
					xmlstr << "<Worker>" << endl;

					xmlstr << "<WorkerPhysicalProperties>" << endl;

					xmlstr << "<Name>" << w->machine_name() << "</Name>" << endl;
						//          xmlstr << "<IPAddress>" << w->get_IPAddress() << "</IPAddress>" << endl; 
					xmlstr << "<Status>" << MWworker_statenames[w->currentState()] << "</Status>" << endl; 
					xmlstr << "<OpSys>" << w->OpSys << "</OpSys>" << endl; 
					xmlstr << "<Arch>" <<  w->Arch << "</Arch>" << endl; 
					xmlstr << "<Bandwidth>" << "N/A" << "</Bandwidth>" << endl; 
					xmlstr << "<Latency>" << "N/A" << "</Latency>" << endl; 
					xmlstr << "<Memory>" << w->Memory << "</Memory>" <<  endl; 
					xmlstr << "<VirtualMemory>" << w->VirtualMemory << "</VirtualMemory>" <<  endl; 
					xmlstr << "<DiskSpace>" << w->Disk << "</DiskSpace>" << endl; 
					xmlstr << "<BenchResult>" << w->get_bench_result() <<"</BenchResult>" << endl;
					xmlstr << "<KFlops>" << w->KFlops <<"</KFlops>" << endl;
					xmlstr << "<Mips>" << w->Mips <<"</Mips>" << endl;   
					xmlstr << "<CPUs>" << w->Cpus <<"</CPUs>" << endl;
					xmlstr << "<NWSinfos>" << "N/A" << "</NWSinfos>" << endl; 

					xmlstr << "</WorkerPhysicalProperties>" << endl;

					xmlstr << "<WorkerUsageProperties>" << endl;

					xmlstr << "<TotalTime>" << "N/A" << "</TotalTime>" << endl;
					xmlstr << "<TotalWorking>" << w->get_total_working() << "</TotalWorking>" << endl;
					xmlstr << "<TotalSuspended>" << w->get_total_suspended() << "</TotalSuspended>" << endl;

					xmlstr << "</WorkerUsageProperties>" << endl;

					xmlstr << "</Worker>" << endl;
				}
			workers->Append ( w );
		}

	xmlstr << "</WorkersList>" << endl;

	xmlstr << "</Workers>" << endl;

		// Global Worker Statistics

	xmlstr << "<GlobalStats>" << endl;

      

	xmlstr << "</GlobalStats>" << endl;

		// Task Pool Information


	xmlstr << "<TaskPoolInfos>" << endl;

	int nrt = get_number_running_tasks();
	int tt = get_number_tasks();
	int tdt = tt - nrt;

	xmlstr << "<TotalTasks>" << tt << "</TotalTasks>" << endl;
	xmlstr << "<TodoTasks>" << tdt << "</TodoTasks>" << endl;
	xmlstr << "<RunningTasks>" << nrt << "</RunningTasks>" << endl;
	xmlstr << "<NumberCompletedTasks>" << num_completed_tasks << "</NumberCompletedTasks>" << endl;
		//xmlstr << "<MaxNumberTasks>" << max_number_tasks << "</MaxNumberTasks>" << endl;
	xmlstr << "</TaskPoolInfos>" << endl;

		// Memory Info

		// Secondary storage info


		// End  XML string

	xmlstr << "</ResourcesStatus>" << endl;
	xmlstr << ends ;

	return xmlstr.str();

}

char* MWDriver::get_XML_results_status(){

	ostrstream xmlstr;

	xmlstr << "<ResultsStatus>" << endl;
	xmlstr << "</ResultsStatus>" << endl;

	xmlstr << ends;

	return xmlstr.str();

}

/* The following (until the end of the file) was written by one of 
   Jean-Pierre's students.  It isn't exactly efficient, dumping
   a condor_status output to a file and reading from that.  I understand
   that it does work, however.  -MEY (8-18-99) */

/* UPDATE:  Note that it doesn't work for flocked machines - GGT 9/26/06 */

int MWDriver::check_for_int_val(char* name, char* key, char* value) {
	if (strcmp(name, key) == 0)
		return atoi(value);
	else return NO_VAL;
}

double MWDriver::check_for_float_val(char* name, char* key, char* value) {
	if (strcmp(name, key) == 0)
		return atof(value);
	else return NO_VAL;  
}

int  MWDriver::check_for_string_val(char* name, char* key, char* value) {
	if (strcmp(name, key) == 0)
		return 0;
	else return NO_VAL;  
}


void MWDriver::get_machine_info()
{
	FILE* inputfile;
	char filename[50];
	char key[200];
	char value[1300];
	char raw_line[1500];
	char* equal_pos;
	int found;
	char temp_str[256];

	char zero_string[2];

	memset(zero_string, 0 , sizeof(zero_string));

	memset(filename, '\0', sizeof(filename));
	memset(key, '\0', sizeof(key));
	memset(value, '\0', sizeof(value));

	strcpy(filename, "/tmp/metaneos_file2");

	memset(temp_str, '\0', sizeof(temp_str));
	sprintf(temp_str, "%s/bin/condor_status -l %s > %s", CONDOR_DIR, mach_name, filename);

	if (system(temp_str) < 0)
		{

			MWprintf( 10, "Error occurred during attempt to get condor_status for %s.\n",  mach_name );
			return;
		}

	if ((inputfile = fopen(filename, "r")) == 0)
		{
			MWprintf( 10, "Failed to open condor_status file!\n");
			return;
		}
	else
		{
			MWprintf( 90, "Successfully opened condor_status file.\n");      
		}

	while (fgets(raw_line, 1500, inputfile) != 0)
		{
			found = 0;
			equal_pos = strchr(raw_line, '=');

			if (equal_pos != NULL)
				{
					strncpy(key, raw_line, equal_pos - (raw_line+1));
					strcpy(value, equal_pos+2);

					if (CondorLoadAvg == NO_VAL && !found)
						{
							CondorLoadAvg = check_for_float_val("CondorLoadAvg", key, value);
							if (CondorLoadAvg != NO_VAL)
								found = 1;
						}

					if (LoadAvg == NO_VAL && !found)
						{
							LoadAvg = check_for_float_val("LoadAvg", key, value);
							if (LoadAvg != NO_VAL)
								found = 1;
						}

					if (Memory == NO_VAL && !found)
						{
							Memory = check_for_int_val("Memory", key, value);
							if (Memory != NO_VAL)
								found = 1;
						}

					if (Cpus == NO_VAL && !found)
						{
							Cpus = check_for_int_val("Cpus", key, value);
							if (Cpus != NO_VAL)
								found = 1;
						}

					if (VirtualMemory == NO_VAL && !found)
						{
							VirtualMemory = check_for_int_val("VirtualMemory", key, value);
							if (VirtualMemory != NO_VAL)
								found = 1;
						}

					if (Disk == NO_VAL && !found)
						{
							Disk = check_for_int_val("Disk", key, value);
							if (Disk != NO_VAL)
								found = 1;
						}

					if (KFlops == NO_VAL && !found)
						{
							KFlops = check_for_int_val("KFlops", key, value);
							if (KFlops != NO_VAL)
								found = 1;
						}

					if (Mips == NO_VAL && !found)
						{
							Mips = check_for_int_val("Mips", key, value);
							if (Mips != NO_VAL)
								found = 1;
						}

	  
					if ( (strncmp(Arch, zero_string, 1) == 0) && !found)
						{	     
							if (check_for_string_val("Arch", key, value) == 0){
								strncpy( Arch, value, sizeof(Arch) );	       
							}
	     
							if (strncmp(Arch, zero_string, 1) != 0)
								found = 1;
						}

					if ( (strncmp(OpSys, zero_string, 1) == 0) && !found)
						{	     
							if (check_for_string_val("OpSys", key, value) == 0){
								strncpy( OpSys, value, sizeof(OpSys) );	       
							}
	     
							if (strncmp(OpSys, zero_string, 1) != 0)
								found = 1;
						}

					if ( (strncmp(IPAddress, zero_string, 1) == 0) && !found)
						{	     
							if (check_for_string_val("StartdIpAddr", key, value) == 0){
								strncpy(IPAddress , value, sizeof(IPAddress) );	       
							}
	     
							if (strncmp(IPAddress, zero_string, 1) != 0)
								found = 1;
						}


					memset(key, '\0', sizeof(key));
					memset(value, '\0', sizeof(value));

				}

		}

	MWprintf(90,"CURRENT MACHINE  : %s \n", mach_name);
	MWprintf(90,"Architecture : %s \n", Arch);
	MWprintf(90,"Operating System : %s \n", OpSys);
	MWprintf(90,"IP address : %s \n", IPAddress);

	MWprintf(90,"CondorLoadAvg : %f\n", CondorLoadAvg);
	MWprintf(90,"LoadAvg : %f\n", LoadAvg);
	MWprintf(90,"Memory : %d\n", Memory);
	MWprintf(90,"Cpus : %d\n", Cpus);
	MWprintf(90,"VirtualMemory : %d\n", VirtualMemory);
	MWprintf(90,"Disk : %d\n", Disk);
	MWprintf(90,"KFlops : %d\n", KFlops);
	MWprintf(90,"Mips : %d\n", Mips);

	fclose( inputfile );

	if (remove(filename) != 0)
		{
			MWprintf(10,"Condor status file NOT removed!\n");
		}
	else
		{
			MWprintf(90,"Condor status file removed.\n");
		}


}
#endif


static char hostname[255];

char* 	
MWDriver::getHostName ()
{
    MWSystem::gethostname(hostname, 255);
	return hostname;
}

int
MWDriver::getNumWorkers()
{
	return workers->number();
}

#ifdef MEASURE
void 
MWDriver::_measure()
{
	time_t now_time = time(0);
	if ( now_time - _measure_last_read_opt_time >= _measure_read_opt_interval ) {
			// Need to read the option file again
		_measure_read_options(_measure_opt_file_name);
		_measure_last_read_opt_time = now_time;
	} else if ( now_time - _measure_last_dump_rec_time >= _measure_dump_rec_interval ) {
		if (_measure_num_task_Done == 0)
			return;

			// Need to dump measurement result again
			// First collect the information
		statsgather(workers, stats);
			
		_measure_master_wall_time += _measure_cur_wall_time();
		_measure_master_cpu_time += _measure_cur_cpu_time();
		if (_measure_master_cpu_time < 0.0001)
			_measure_master_cpu_time = 0.0;
		if (_measure_master_recv_cpu_time < 0.0001)
			_measure_master_recv_cpu_time = 0.0;
		if (_measure_master_act_on_completed_task_cpu_time < 0.0001)
			_measure_master_act_on_completed_task_cpu_time = 0.0;
		if (_measure_task_MP_worker_cpu_time < 0.0001)
			_measure_task_MP_worker_cpu_time = 0.0;

	    	// Dump record to file
		_measure_dump_records();

		if (_measure_use_adaptation == 1) {
				// Adapt the number of workers when necessary
			_measure_adapt();
		}

			// Reset the measurement numbers for the next run
		_measure_reset();
	}
}

void MWDriver::_measure_reset()
{
		// Reset the measurements for the next interval
	MWprintf(91, "MEASURE: reset all measurement related values to default.\n");
	_measure_num_task_Done = 0;
	_measure_master_wall_time = 0.0 - _measure_cur_wall_time();
	_measure_master_cpu_time = 0.0 - _measure_cur_cpu_time();
	_measure_master_recv_time = 0.0;
	_measure_master_recv_cpu_time = 0.0;
	_measure_master_act_on_completed_task_time = 0.0;
	_measure_master_act_on_completed_task_cpu_time = 0.0;
	_measure_task_RoundTrip_time = 0.0;
	_measure_task_MP_master_time = 0.0;
	_measure_task_MP_worker_time = 0.0;
	_measure_task_MP_worker_cpu_time = 0.0;
	_measure_task_Exe_wall_time = 0.0;
	_measure_task_Exe_cpu_time = 0.0;
	_measure_target_num_workers = 99999;
}

void
MWDriver::_measure_dump_header() 
{
		// Create the file and write the first several lines		
	FILE *rec;
	if ( (rec = fopen(_measure_rec_file_name, "a")) == NULL ) {
		MWprintf(31, "MEASURE|Failed to open %s for measurement records! errno = %d\n", 
				 _measure_rec_file_name, errno);
	} else {
		fprintf(rec, "Start_time = %s\n", ctime(&_measure_last_read_opt_time));
		fprintf(rec, "NOW_TIME LastDump #Wkr #RUN #TODO #Done(sum) *M_wall(-CPU-) *M_recv(-CPU-) *W_recv(-CPU-) *M_exec(-CPU-) *W_exec(-CPU-) #MSGS\n");
		fclose(rec);
	}
}

void 
MWDriver::_measure_dump_records()
{
	time_t now_time = time(0);
	
	FILE *rec;
	if ( (rec = fopen(_measure_rec_file_name, "a")) == NULL ) {
		MWprintf(31, "MEASURE|Failed to open %s for measurement records! errno = %d \n", _measure_rec_file_name, errno);
		return;
	} 
	
		// NOW_TIME LastDump 
	char *now = ctime(&now_time);
	now[19] = '\0';
	fprintf(rec, "%s ", &now[11]);
	char *last = ctime(&_measure_last_dump_rec_time);
	last[19] = '\0';
	fprintf(rec, "%s ", &last[11]);
		// Then the numbers
        // #Wkr #RUN #TODO #Done(sum) *M_wall(-CPU-) *M_recv(-CPU-) *W_recv(-CPU-) *M_exec(-CPU-) *W_exec(-CPU-) *T_RTrp
	if (_measure_num_task_Done) {
		fprintf(rec, "%3d %3d %6d %3d(%6d) %7.2f(%5.1f) %7.2f(%5.1f) %7.2f(%5.1f) %7.2f(%5.1f) %7.2f(%5.1f)\n",
				(workers->number()>999)?999:workers->number(),		// %3d #Wkr
				(running->number()>999)?999:running->number(), 		// %3d #RUN
				(todo->number()>999999)?999999:todo->number(), 		// %6d #TODO
				_measure_num_task_Done,					// %3d #Done
				(num_completed_tasks>999999)?999999:num_completed_tasks,// %6d (sum)
			
				_measure_master_wall_time, 				// %7.2f *M_wall
				_measure_master_cpu_time,				// %5.1f (-CPU-)
				_measure_master_recv_time,				// %7.2f *M_recv
				_measure_master_recv_cpu_time,				// %5.1f (-CPU-)
				_measure_task_MP_worker_time/_measure_num_task_Done,	// %7.2f *W_recv 
				_measure_task_MP_worker_cpu_time/_measure_num_task_Done,// %5.1f (-CPU-)
				_measure_master_act_on_completed_task_time, 		// %7.2f *M_exec
				_measure_master_act_on_completed_task_cpu_time,		// %5.1f (-CPU-)
				_measure_task_Exe_wall_time/_measure_num_task_Done, 	// %7.2f *T_exec
				_measure_task_Exe_cpu_time/_measure_num_task_Done 	// %5.1f (-CPU-)
					// 	_measure_task_RoundTrip_time/_measure_num_task_Done 	// %7.2f *T_RTrp
				);
		MWList<void> * recv_buf_list = RMC->recv_buffers();
		if (recv_buf_list) {
			fprintf(rec, "BUFFER: size = %4d\n", RMC->recv_buffers()->number() );
		}
	} else {
		fprintf(rec, "%3d %3d %6d %3d(%6d) %7.2f(%5.1f) %7.2f(%5.1f) %7.2f(%5.1f) %7.2f(%5.1f) %7.2f(%5.1f)\n",
				(workers->number()>999)?999:workers->number(),		// %3d #Wkr
				(running->number()>999)?999:running->number(), 		// %3d #RUN
				(todo->number()>999999)?999999:todo->number(), 		// %6d #TODO
				_measure_num_task_Done,					// %3d #Done
				(num_completed_tasks>999999)?999999:num_completed_tasks,// %6d (sum)
			
				_measure_master_wall_time, 				// %7.2f *M_wall
				_measure_master_cpu_time,				// %5.1f (-CPU-)
				_measure_master_recv_time,				// %7.2f *M_recv
				_measure_master_recv_cpu_time,				// %5.1f (-CPU-)
				0.0,							// %7.2f *W_recv 
				0.0,							// %5.1f (-CPU-)
				_measure_master_act_on_completed_task_time, 		// %7.2f *M_exec
				_measure_master_act_on_completed_task_cpu_time,		// %5.1f (-CPU-)
				0.0,							// %7.2f *T_exec
				0.0							// %5.1f (-CPU-)
					// 0.0  						// %7.2f *T_RTrp
				);
		MWList<void> * recv_buf_list = RMC->recv_buffers();
		if (recv_buf_list) {
			fprintf(rec, "BUFFER: size = %4d\n", RMC->recv_buffers()->number() );
		}
	}
		
	fclose(rec);
	_measure_last_dump_rec_time = now_time;

		/*
		  if (num_completed_tasks > 20000) {
		  MWprintf(31, "Jichuan's hack - exit when running too long!");
		  assert(0);
		  exit(1);
		  }
		*/
	
}

void 
MWDriver::_measure_read_options(const char* opt_fname)
{
	FILE *opt;
	if ( (opt = fopen(opt_fname, "r")) == NULL ) {
		MWprintf (31, "MEASURE|%s not found, uses default.\n", opt_fname);
		return;
	}
	
	char str[512];
	char *ret = NULL;
	while (fgets(str, 512, opt)) {
		if ((str[0] == '#')||(str[0] == ' ')) // comments
			continue;
		ret = strchr(str, '\n');
		if (ret) *ret = '\0';
		if (strstr(str, "dump_record_interval") != NULL) {
			_measure_dump_rec_interval = atoi( &(str[strlen("dump_record_interval")+1]) );
			MWprintf(31, "MEASURE|Setting dump_record_interval = %d\n", _measure_dump_rec_interval);
		} else if (strstr(str, "read_option_interval") != NULL) {
			_measure_read_opt_interval = atoi( &(str[strlen("read_option_interval")+1]) );
			MWprintf(31, "MEASURE|Setting read_option_interval = %d\n", _measure_read_opt_interval);
		} else if (strstr(str, "use_adaptation") != NULL) {
			_measure_use_adaptation = 1;
			MWprintf(21, "MEASURE|Setting use_adaptation to yes. \n");
		}

	}
	fclose(opt);
}

void 
MWDriver::_measure_remove_worker()
{
	MWWorkerID * w = _measure_current_worker;
	MWprintf(31, "MEASURE|ADAPT|Removing worker id1 = %08x \n", w->get_id1());
	RMC->initsend ( );
	RMC->send ( w->get_id1(), KILL_YOURSELF );
	worker_last_rites ( w );
}

void 
MWDriver::_measure_adapt()
{
	if (_measure_num_task_Done == 0)  {
		return;
	}
	
	if (_measure_target_num_workers < workers->number()) {
			// need to remove more workers, still adapting ...
		_measure_remove_worker();		
		return;
	}
		
	_measure_master_wait_rate = _measure_master_recv_time/_measure_master_wall_time;
	if (_measure_master_wait_rate  < 0.001) { // then the workers MIGHT BE WAITING
		_measure_master_busy_times ++ ; 
	} else if (_measure_master_wait_rate > 0.50) { // then the master is WAITING
		_measure_master_wait_times ++;
	}
				
	if (_measure_master_busy_times > workers->number()/10 ) {
			// reduce target worker numbers
		_measure_master_busy_times = 0;
		
		_measure_target_num_workers = (workers->number()>12)? workers->number()-12 : 6;
		RMC->set_target_num_workers(_measure_target_num_workers);
		
		_measure_remove_worker();
		
		call_hostaddlogic();
		MWprintf(31, "MEASURE|ADAPT|decrease worker number to %d\n", _measure_target_num_workers);
	}
	
	if (_measure_master_wait_times > 10) {
			// ask for more workers
		_measure_master_wait_times = 0;
		
		if (RMC->get_num_exec_classes() == 0) {
			int target = RMC->get_target_num_workers(-1);
			if (workers->number() == target) {
				RMC->set_target_num_workers(target + 6);
				MWprintf(31, "MEASURE|ADAPT|increase worker number to %d\n", target+6);
			}
		}
	}
}
#endif // MEASURE

void 
MWDriver::prepare_update(int mode)
{
	MWWorkerID *w = workers->First();
	switch(mode)
		{
		case -2:
			while ( workers->AfterEnd() == false )
				{
					if( w->currentState() != BENCHMARKING && w->currentState() != INITIALIZING )
						{
							w = workers->Current();
							RMC->initsend();
							pack_update();
							RMC->send(w->get_id1(), UPDATE_FROM_DRIVER);
						}
					workers->Next();
				}
			break;
		case -1:
			RMC->initsend();
			pack_update();
			send_update_message();
			break;
		default:
			RMC->initsend();
			pack_update();
			send_update_message_to(mode);
		}
}

/* send a UPDATE_FROM_DRIVER message to MWcurrentWorker, assuming the user has 
 * already called RMC->initsend(), and RMC->pack() for the read data. */
int 
MWDriver::send_update_message()
{
	if (MWcurrentWorker == NULL)
		return -1;

	return RMC->send(MWcurrentWorker->get_id1(), UPDATE_FROM_DRIVER);
}

/* send a UPDATE_FROM_DRIVER message to a specific worker, assuming the user has 
 * already called RMC->initsend(), and RMC->pack() for the read data. */
int
MWDriver::send_update_message_to(int worker)
{
	return RMC->send(worker, UPDATE_FROM_DRIVER);
}

/*
  Local Variables:
  mode: c++
  eval: (setq c-basic-offset 4)
  eval: (setq c-comment-only-line-offset 4)
  eval: (setq c-indent-level 4)
  eval: (setq c-brace-imaginary-offset 0)
  eval: (setq c-brace-offset 0)
  eval: (setq c-argdecl-indent 0)
  eval: (setq c-label-offset -4)
  eval: (setq c-continued-statement-offset 4)
  eval: (setq c-continued-brace-offset -4)
  eval: (setq c-tab-always-indent nil)
  eval: (setq tab-width 4)
  End:
*/

void
MWDriver::increasePendingNum(void)
{
  	pthread_mutex_lock(&m_pending_mutex);
  	m_pending_num++;
  	pthread_mutex_unlock(&m_pending_mutex);
}
void
MWDriver::decreasePendingNum(void)
{
  	pthread_mutex_lock(&m_pending_mutex);
  	m_pending_num--;
	pthread_cond_signal(&m_pending_cond);
  	pthread_mutex_unlock(&m_pending_mutex);
}

void 
MWDriver::getGlobalLock(void)
{
	pthread_mutex_lock(&global_mw_mutex);
}

void 
MWDriver::getGlobalUnLock(void)
{
	pthread_mutex_unlock(&global_mw_mutex);
}
