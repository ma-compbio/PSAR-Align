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
#include "MW.h"
#include "MWWorker.h"
#include "MWSystem.h"
#include "RMComm/MWRMComm.h"

MWWorker::MWWorker() 
{
	master = UNDEFINED;
	workingTask = 0; // This must be initialized by subclass
}

MWWorker::~MWWorker() 
{
}

double
MWWorker::benchmark(MWTask *t)
{
    if (!t) return 1.0;
        
    double begTime = MWSystem::gettimeofday();
    execute_task( t );
    double endTime = MWSystem::gettimeofday();
    return 1.0/(endTime-begTime);

}

void MWWorker::go( int argc, char *argv[] ) 
{
	MWprintf ( 10, "About to call setup\n");

	if (workingTask == 0) {
		MWprintf( 1, "Fatal Error - Worker did not set field workingTask to valid value\n");
		exit(-1);
	}

	RMC->setup( argc, argv, &workerid, &master );
	MWprintf ( 10, "Worker %x started.\n", workerid );
	
	 greet_master();
#ifndef INDEPENDENT
	worker_mainloop();
#endif
}

MWReturn MWWorker::greet_master() {
	MWSystem::gethostname ( mach_name, 64 );

	MWprintf ( 10, "Worker started on machine %s.\npid = %d\n", mach_name, MWSystem::getpid() );
	
		/* Pack and send to the master all these information
		   concerning the host and the worker specificities  */
	RMC->initsend();
	RMC->pack( mach_name );
	pack_worker_initinfo();
	
	int status = RMC->send( master, INIT );
	MWprintf ( 10, "Sent the master %x an INIT message.\n", master );
	
	if ( status < 0 ) {
		MWprintf ( 10, "Had a problem sending my name to master.  Exiting.\n");
		RMC->exit(1);
	}

	
#ifdef INDEPENDENT
    return OK;
#else
    return do_benchmark_cmd();
#endif
}

MWReturn MWWorker::do_benchmark_cmd ( )
{
	int status;
	int len, tag, tid;
	
	// wait for the setup info from the master 
	int buf_id = RMC->recv( master, -1 );  
	if( buf_id < 0 ) {
		MWprintf ( 10, "Had a problem receiving INIT_REPLY.  Exiting.\n" );
		RMC->exit( INIT_REPLY_FAILURE );
	}
	status = RMC->bufinfo ( buf_id, &len, &tag, &tid );
	MWprintf ( 10, "Got Something from the master in reply to INIT %d\n", tag);
	
	MWReturn ustat = OK;
		// unpack initial data to set up the worker state

	switch ( tag )
	{
		case INIT_REPLY:
		{

			if ( RMC->unpack ( master_mach_name ) != 0 )
			{
				int err = -1;
				MWprintf ( 10, "Error unpacking master hostname. \n");
				RMC->initsend ( );
				RMC->pack ( &err, 1 );
				RMC->send ( master, BENCH_RESULTS );
				return ustat;
			}

			if ( (ustat = unpack_init_data()) != OK ) {
				int err = -1;
				MWprintf ( 10, "Error unpacking initial data.\n" );
				RMC->initsend();
				RMC->pack( &err, 1 );
				RMC->send( master, BENCH_RESULTS );
				return ustat;
			}

			int bench_tf = FALSE;
			RMC->unpack( &bench_tf, 1 );

			if ( bench_tf ) {
				MWprintf ( 10, "Recvd INIT_REPLY, now benchmarking.\n" );
				workingTask->unpack_work();
				double bench_result = benchmark( workingTask );

				MWprintf ( 40, "Benchmark completed....%f\n", bench_result );
				int zero = 0;
				RMC->initsend();
				RMC->pack( &zero, 1 );  // zero means that unpack_init_data is OK.
				RMC->pack( &bench_result, 1 );
			} else {
				MWprintf ( 10, "Recvd INIT_REPLY, no benchmark.\n" );
				double z = 0.0;
				int zero = 0;
				RMC->initsend();
				RMC->pack( &zero, 1 );  // zero means that unpack_init_data is OK.
				RMC->pack( &z, 1 );
			}
			MWprintf ( 10, "Worker Sending BENCH_RESULTS\n");
			RMC->send( master, BENCH_RESULTS );

			return ustat;
			break;
		}

		case CHECKSUM_ERROR:
		{
			MWprintf ( 10, "Got a checksum error\n");
			RMC->exit( CHECKSUM_ERROR_EXIT );
		}
	}
	return OK;
}


void MWWorker::worker_mainloop() {
	
	for (;;) {
      worker_mainloop_oneshot();
    }
}
		
MWReturn MWWorker::worker_mainloop_oneshot () 
{
	int status = 0, len = 0, tag = 0, tid = 0;
    double wall_time = 0.0;
    double cpu_time = 0.0;
	
	int curSubTask = 0;

#ifdef MEASURE
		double _recv_start_time = _measure_cur_wall_time();
		double _recv_start_cpu_time = _measure_cur_cpu_time();
#endif // MEASURE	
	
		// get message from master. block if task done
        int buf_id;
		buf_id = RMC->recv ( master, -1 );
		if( buf_id < 0 ) {
		  MWprintf( 10, "Could not receive message from master, ret = %d.  Exiting\n", buf_id );
		  RMC->exit( buf_id );
		}

#ifdef MEASURE
		double _task_recv_time = _measure_cur_wall_time() - _recv_start_time;
		double _task_recv_cpu_time = _measure_cur_cpu_time() - _recv_start_cpu_time;
#endif // MEASURE

		status = -2; len = -2; tag = -2; tid = -2;		
		status = RMC->bufinfo ( buf_id, &len, &tag, &tid );

		switch ( tag ) {
			
		case RE_INIT: {   /* This can happen:  the lower level can tell us 
							 that the master has gone down and come back 
							 up, and we hve to re-initialize ourself. */
			greet_master();
            do_benchmark_cmd();
			break;
		}

		case REFRESH:
		{
			unpack_init_data ( );
			break;
		}

		case DO_THIS_WORK: {
			MWprintf(51,"received DO_THIS_WORK signal\n");
/*
time_t now = time(0);
while(time(0)-now<30){}	
*/
			int num = UNDEFINED;
			MWTaskType thisTaskType;
			//moved out of switch statement
			//double wall_time = 0.0;
			//double cpu_time = 0.0;
			int tstat;
			int mytemp;

			MWTask *curTask = NULL;
			
			tstat = RMC->unpack ( &mytemp, 1, 1);

			if ( tstat != 0 )
			{
				MWprintf ( 10, "Error: The receive buffer not unpacked on %d\n",
					mach_name );
				fflush ( stdout );
				RMC->exit ( UNPACK_FAILURE );
			}
			thisTaskType = (MWTaskType)mytemp;

			switch ( thisTaskType )
			{
				case MWNORMAL:
				  //curTask = workingTaskContainer->Current();
					break;
				case MWNWS:
#ifndef NWSENABLED
/*
					MWprintf(31, "Got a MWNWS task. \n");
					controlTask = getNWSTask ( );
					curTask = controlTask;
*/
#endif
					break;
					
				default:
					MWprintf(30, "MWWorker::worker_mainloop_ind() - other task type\n");
			}
			curTask = workingTask;

			curTask->taskType = thisTaskType;
			tstat = RMC->unpack( &num, 1, 1 );
			if( tstat != 0 ) {
			  MWprintf( 10, "Error.  The receive buffer not unpacked on %s\n", 
				    mach_name );
			  fflush( stdout );
			  RMC->exit( UNPACK_FAILURE );
			}
			curTask->number = num; 

			unpack_driver_task_data();
			curTask->unpack_work();
			
			MWprintf(51,"Finished unpacking work, now executing...\n");			
			/* Set our stopwatch.  :-) */
			wall_time = - MWSystem::gettimeofday();
			cpu_time = - MWSystem::getcputime();
		      
				/* do it! */
			switch ( thisTaskType )
			{
				case MWNORMAL:
						execute_task(curTask);
						// Record times for it...
                		wall_time += MWSystem::gettimeofday();
                		cpu_time +=  MWSystem::getcputime();
/*
				      RMC->initsend();
				      RMC->pack( &(curTask->number), 1, 1 );   

						  // GGT Not sure what to do here..

					  RMC->pack(&curSubTask, 1, 1);
				      curTask->pack_subresults(curSubTask);
					  
				      status = RMC->send(master, SUBRESULTS);
				      if ( status < 0 ){
					MWprintf ( 10, "Bummer!  Could not send SUBresults of task %d\n", num );
					MWprintf( 10, "Exiting worker!" );
					RMC->exit( FAIL_MASTER_SEND );
				      }
*/
            // if all subtasks done, do task and send RESULT message to master
            // worker then goes back into blocking recv. This gets called if there is only 1 task in the container
			{
                // Record times for it...
            	//wall_time += MWSystem::getcputime();
            	//cpu_time += MWSystem::gettimeofday();

                RMC->initsend();
				RMC->pack( &num, 1, 1);
#ifdef MEASURE
                RMC->pack(&_task_recv_time, 1, 1);
                RMC->pack(&_task_recv_cpu_time, 1, 1);
#endif // END MEASURE
                RMC->pack( &wall_time, 1, 1 );
                RMC->pack( &cpu_time, 1, 1 );
                curTask->pack_results();
                status = RMC->send(master, RESULTS);
                MWprintf(51,"sent result\n");
			}

				      break;
				case MWNWS:
#ifndef NWSENABLED
/*
					curTask->printself ( 10 );
					execute_nws ( curTask );
					curTask->printself ( 10 );
*/
#endif
					break;
				default:
					MWprintf ( 10, "Unidentified %d \n", thisTaskType);
					exit(1);
			}

/*
			if(curTask->numsubtask == -1) // partial results disabled
			{	
				// Record times for it... 
			wall_time += _measure_cur_wall_time();
			cpu_time += _measure_cur_cpu_time();

				// Now send... 
			RMC->initsend();
			//RMC->pack( &num, 1, 1 );
#ifdef MEASURE
			RMC->pack(&_task_recv_time, 1, 1);
			RMC->pack(&_task_recv_cpu_time, 1, 1);
#endif // END MEASURE
			RMC->pack( &wall_time, 1, 1 );
			RMC->pack( &cpu_time, 1, 1 );
			curTask->pack_results();
			
			status = RMC->send(master, RESULTS);

			if ( status < 0 ){
			  MWprintf ( 10, "Bummer!  Could not send results of task %d\n", num ); 
			  MWprintf( 10, "Exiting worker!" ); 
			  RMC->exit( FAIL_MASTER_SEND );
			}
                        
		
			MWprintf ( 40, "%s sent results of job %d.\n", 
					   mach_name, curTask->number );
			}
*/
			switch ( thisTaskType )
			{
				case MWNORMAL:
					break;
				case MWNWS:
#ifndef NWSENABLED
					delete curTask;
					controlTask = NULL;
#endif
					break;

				default:
					MWprintf(30, "MWWorker::worker_mainloop_ind() - other task type\n");
			}
			break;
			
		}
        case NO_MESSAGE:{
            MWTask *curTask = workingTask;
			int num = curTask->number;

            /* Set our stopwatch. again  :-) */
            wall_time -= MWSystem::gettimeofday();
            cpu_time -= MWSystem::getcputime();
	
			execute_task(curTask);

			// Record times for it...
			wall_time += MWSystem::gettimeofday();
			cpu_time += MWSystem::getcputime();

            curSubTask++;

            if(curSubTask==curTask->numsubtask)
            {   MWprintf(50,"sending result\n");
                curSubTask = 0;
                execute_task( curTask );
                RMC->initsend();
                RMC->pack( &num, 1, 1 );
                #ifdef MEASURE
                    RMC->pack(&_task_recv_time, 1, 1);
                    RMC->pack(&_task_recv_cpu_time, 1, 1);
                #endif // END MEASURE
                RMC->pack( &wall_time, 1, 1 );
                RMC->pack( &cpu_time, 1, 1 );
                curTask->pack_results();
                status = RMC->send(master, RESULTS);
                MWprintf(50,"sent result\n");
            }  

            break;
		}
		case UPDATE_FROM_DRIVER:
		{
			MWprintf(51,"Received update while not working on task\n");
			break;
		}       
		case STOP_WORK: {
// GGT This shouldn't happen anymore -- ignore it

// 		  int tasknum;
// 		  RMC->unpack(&tasknum, 1, 1);
// 		  // stop our current task only if the signal is for current task (not a previous one) and task is still running
// 		  if(tasknum == curTask->number)
// 		    {
// 		      //curSubTask = 0;
// 		      //execute_task( workingTaskContainer );
// 		      RMC->initsend();
// 		      RMC->pack( &workingTaskContainer->number, 1, 1 );
// #ifdef MEASURE
// 		      RMC->pack(&_task_recv_time, 1, 1);
// 		      RMC->pack(&_task_recv_cpu_time, 1, 1);
// #endif // END MEASURE
// 		      RMC->pack( &wall_time, 1, 1 );
// 		      RMC->pack( &cpu_time, 1, 1 );
// 		      workingTaskContainer->pack_results();
// 		      status = RMC->send(master, RESULTS);
// 		      MWprintf(50,"sent result\n");
//		    }

		  break;
		}
		case KILL_YOURSELF: {
			suicide();
		}
		case CHECKSUM_ERROR:
		{
			MWprintf ( 10, "Got a checksum error\n");
			RMC->exit( CHECKSUM_ERROR_EXIT );
		}
		default: {
			MWprintf ( 10, "Received strange command %d.\n", tag );
			RMC->exit( UNKNOWN_COMMAND );
		}
		} // switch

	return OK;

}

/* We've received orders to kill ourself; we're not needed anymore.
   Fall on own sword. */

void MWWorker::suicide () 
{   
	MWprintf ( 10, "\"Goodbye, cruel world...\" says %s\n", mach_name );
	
	RMC->exit(0);
}


/* For steering */
int 
MWWorker::send_update_message()
{
	return RMC->send(master, UPDATE_FROM_WORKER);	
}

void
MWWorker::check_update_message()
{
	int bid = 0 ;
	int info, len, tag, sending_host;

	RMC->nrecv(-1, UPDATE_FROM_DRIVER);
	info = RMC->bufinfo(bid, &len, &tag, &sending_host);
	if (tag == UPDATE_FROM_DRIVER)
	{
		unpack_update();
	}
	else if(tag != NO_MESSAGE)
		MWprintf(10, "check_update_message: \
			got a message which is not UPDATE_FROM_DRIVER.\n");
/* moving recv_all functionality into recv
	MWList<void> * recv_buf_list = RMC->recv_buffers();
	
	RMC->recv_all(-1, UPDATE_FROM_DRIVER);
       	len = recv_buf_list->number();
	recv_buf_list->First();	
	while (recv_buf_list->AfterEnd() == false) {
		bid = (int*)recv_buf_list->Current();
		info = RMC->bufinfo(*bid, &len, &tag, &sending_host);
		if (tag != UPDATE_FROM_DRIVER) {
			MWprintf(10, "check_update_message: \
				got a message which is not UPDATE_FROM_DRIVER.\n");
			recv_buf_list->Remove();
		}
		recv_buf_list->Next();
	}
	return RMC->recv_buffers()->number();
*/
}
// Local Variables:
// mode: c++
// eval: (c-set-style "K&R")
// eval: (setq fill-column 79)
// eval: (setq c-basic-offset 2)
// eval: (setq tab-width 4)
// eval: (c-set-offset 'arglist-close 0)
// eval: (setq indent-tabs-mode nil)
// End:
