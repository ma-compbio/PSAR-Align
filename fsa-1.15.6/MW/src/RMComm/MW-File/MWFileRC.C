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
#include "MWFileTypes.h"
#include "MWFileSend.h"
#include "MWFileRC.h"
#include "MWFileRCSymbol.h"
#include "MWFileError.h"

#include <assert.h>
#include <time.h>

#ifdef FILE_MASTER
#include <user_log.c++.h>
#endif

#include <ctype.h>
#include <errno.h>

#include <MWSystem.h>
#include <MWDriver.h>
#include <MWWorkerID.h>
#include <MWTask.h>
#include <MWWorker.h>

#ifdef USE_CHIRP
#include "chirp_client.h"
#endif

extern int *MW_exec_class_num_workers;

// These static declarations are now necessary in the RMComm implementation
#ifdef FILE_MASTER
MWRMComm * MWDriver::RMC = new MWFileRC( TRUE, 0 );
MWRMComm * MWTask::RMC = MWDriver::RMC;
MWRMComm * MWWorker::RMC = NULL;
#else
MWRMComm * MWWorker::RMC = new MWFileRC ( FALSE, 0 );
MWRMComm * MWTask::RMC = MWWorker::RMC;
MWRMComm * MWDriver::RMC = NULL;
#endif

/// We need to do something better with this
const int want_checkpointing = 0;

/** If target_num_hosts is larger than this, we cut it back, since
    we need to allocate arrays of reasonable size.
*/
const int REAL_TARGET_NUM_WORKERS = 2048;

///
const unsigned int SLEEP_DELAY = 5;

int workerEvents[FileNumEvents] = { HOSTDELETE, HOSTADD, HOSTDELETE, HOSTSUSPEND, HOSTRESUME, CHECKSUM_ERROR };

int min( int *arr, int len ) 
{
                // return index of min element in array of length len.
    int min = 0;
    for ( int i=1 ; i<len ; i++ ) 
    {
        if ( arr[i] < arr[min] ) 
	{
             min = i;
        }
    }
    return min;
}

/* // To make gcc -Wall happy
static void Truncate ( char *file )
{
    if ( truncate ( file, 0 ) != 0 );
    	//MWprintf ( 10, "Could not truncate file %s. Errno is %d\n", file, errno );
}
*/

static void HardDeleteFile ( char *file )
{
	if ( !file ) return;
	remove(file);
}

MWFileRC::MWFileRC ( bool val, int id )
{
	isMaster = val;
	FileRCID = id;
	sendList = NULL;
	recvList = NULL;
}

void
MWFileRC::exit ( int retval )
{

#ifndef FILE_MASTER
    ::exit( retval );
#else

//    struct stat statbuf;
//    char cmd[_POSIX_PATH_MAX];

    for ( int i = 0; i < target_num_workers; i++ )
    {
	if ( &fileWorkers[i] && fileWorkers[i].state != FILE_FREE )
	    killWorker ( i );
    }

#if 1
    // Jeff doesn't want to remove the directories for now...

    /* We have to unlink the directories */
 /*   
	if ( stat ( output_directory, &statbuf ) == 0 )
    {
    	strcpy ( cmd, "/bin/rm -rf " );
    	strcat ( cmd, output_directory );
    	if ( system ( cmd ) < 0 )
	    MWprintf ( 10, "Cannot unlink %s directory\n", output_directory );
    }

    if ( stat ( input_directory, &statbuf ) == 0 )
    {
    	strcpy ( cmd, "/bin/rm -rf " );
    	strcat ( cmd, input_directory );
    	if ( system ( cmd ) < 0 )
	    MWprintf ( 10, "Cannot unlink %s directory\n", input_directory );
    }

    if ( stat ( control_directory, &statbuf ) == 0 )
    {
    	strcpy ( cmd, "/bin/rm -rf " );
    	strcat ( cmd, control_directory );
	if ( system ( cmd ) < 0 ) 
	    MWprintf ( 10, "Cannot unlink %s directory\n", control_directory );
    }
*/
#endif

    ::exit( retval );
#endif
}


int
MWFileRC::setup( int argc, char *argv[], int *mytid, int *master_tid )
{
    SIMUL_SEND = 100;
    worker_timeout = 80*60;

    if ( isMaster == TRUE )
    {

	strcpy ( output_directory, "worker_output" );
	strcpy ( input_directory, "worker_input" );
	strcpy ( control_directory, "submit_files" );
	strcpy ( moment_worker_file, "moment_worker_file" );
	cyclePosition = 0;
	turnNo = 0;
	*mytid = FileRCID;
	*master_tid = FileRCID;
	CHECKLOG_FREQ = 100;
	subId = 0;
	hostadd_reqs = NULL;

	MWSystem::mkdir ( output_directory);
	MWSystem::mkdir ( input_directory);
    MWSystem::mkdir ( control_directory);
    }
    else
    {
#ifdef USE_CHIRP
	strcpy ( output_directory, "worker_output" );
	strcpy ( input_directory, "worker_input" );
	strcpy ( control_directory, "submit_files" );
	strcpy ( moment_worker_file, "moment_worker_file" );
	
	MWSystem::mkdir( output_directory);
	MWSystem::mkdir( input_directory);
	MWSystem::mkdir( control_directory);
#endif

	// Slave FileRC instance. Only now the argc and argv make any sense
	FileRCID = atoi ( argv[1] );
	expected_number = atoi ( argv[2] );
	master_expected_number = atoi(argv[3] );
	strcpy ( output_directory, argv[4] );
	strcpy ( input_directory, argv[5] );
		strcpy ( control_directory, "submit_files" );
		*master_tid = atoi ( argv[6] );
		*mytid = FileRCID;
		}

		sendList = recvList = NULL;
		return 0;
	}


	//---------------------------------------------------------------------
	// MWFileRC::init_beginning_workers():
	//	Called every time a master gets up either at the init time or 
	// aftermath a brutal crash. It is here that we do a resuscicate and
	// fill our worker lists.
	// We start off by generating all the needed workers.
	//---------------------------------------------------------------------
	int
	MWFileRC::init_beginning_workers ( int *nworkers, MWWorkerID ***workers )
	{
		int j = 0;
		int i;

		*nworkers = j;

		if ( !hostadd_reqs )
		{
			// It is not restarting from checkpoint
			MWprintf ( 10, "In MWFileRC::init_beginning_workers()\n");
			hostadd_reqs = new int[exec_classes];
		
			for ( i = 0 ; i < exec_classes ; i++ ) 
			{
				hostadd_reqs[i] = 0;
			}

			fileWorkers = new struct FileWorker[target_num_workers];

			if ( !fileWorkers )
			{
				MWprintf ( 10, "Cannot allocate memory\n");
				return -1;
			}
			for ( i = 0; i < target_num_workers; i++ )
			{
				fileWorkers[i].state = FILE_FREE;
				fileWorkers[i].served = -1;
				fileWorkers[i].event_no = 0;
				fileWorkers[i].id = -1;
				fileWorkers[i].counter = 0;
				fileWorkers[i].worker_counter = 0;
				fileWorkers[i].arch = -1;
				fileWorkers[i].condorID = -1;
				fileWorkers[i].condorprocID = -1;
			}
		}

		if ( j == 0 ) return 0;

		return 0;
	}

int
MWFileRC::restart_beginning_workers ( int *nworkers, MWWorkerID ***workers, MWmessages msg )
{
	int j = 0;
	int i;
	FILE *temp;
	char worker_waitfile[_POSIX_PATH_MAX];
	char master_waitfile[_POSIX_PATH_MAX];

	for ( i = 0; i < target_num_workers; i++ )
	{
		if ( fileWorkers[i].state == FILE_RUNNING || fileWorkers[i].state == FILE_SUSPENDED || 
					fileWorkers[i].state == FILE_RESUMED || fileWorkers[i].state == FILE_EXECUTE )
		{
			j++;
		}
	}

	*nworkers = j;
	if ( j == 0 ) return 0;

	(*workers) = new MWWorkerID*[j];

	j = 0;
	for ( i = 0; i < target_num_workers; i++ )
	{
		if ( fileWorkers[i].state == FILE_RUNNING || fileWorkers[i].state == FILE_SUSPENDED || 
					fileWorkers[i].state == FILE_RESUMED || fileWorkers[i].state == FILE_EXECUTE )
		{
			(*workers)[j] = new MWWorkerID;
			(*workers)[j]->set_id1 ( fileWorkers[i].id );
			(*workers)[j]->set_id2 ( fileWorkers[i].id );
			(*workers)[j]->set_exec_class ( fileWorkers[i].exec_class );
			j++;

			sprintf ( worker_waitfile, "./%s/worker_waitfile.%d", input_directory, fileWorkers[i].id );
			sprintf ( master_waitfile, "./%s/master_waitfile.%d", output_directory, fileWorkers[i].id );
			temp = Open ( worker_waitfile, "a" );
			if ( temp == NULL )
			{
				MWprintf ( 10, "Oops couldn't open the worker_waitfile for syncing\n");
				continue;
			}
			fprintf( temp, "%c %d %c", SYNC, fileWorkers[i].worker_counter + SIMUL_SEND, ESYNC );
			fileWorkers[i].worker_counter += SIMUL_SEND;
			fflush( temp );
			Close( temp );

			initsend ( );
			send ( i, msg );
			fileWorkers[i].counter = 0;
			fileWorkers[i].worker_counter = 0;
			HardDeleteFile ( master_waitfile );
		}
	}
	MWprintf (10, "Created in init %d workers\n", j );

	return 0;
}


	//---------------------------------------------------------------------------
	// MWFileRC::start_worker ()
	//	In responce to a HOSTADD message. You know in whomRecv who has started.
	//---------------------------------------------------------------------------
	int
	MWFileRC::start_worker ( MWWorkerID *w )
	{
		if ( !w )
		{
			MWprintf ( 10, "WorkerID w cannot be null in start_worker\n");
			return -1;
		}

		w->set_id1 ( fileWorkers[whomRecv].id );
		w->set_id2 ( fileWorkers[whomRecv].id );
		w->set_arch ( fileWorkers[whomRecv].arch );
		w->set_exec_class ( fileWorkers[whomRecv].exec_class );

		hostadd_reqs[fileWorkers[whomRecv].exec_class]--; 
		return 0;
	}

	int MWFileRC::removeWorker ( MWWorkerID *w )
	{
		if ( !w ) 
		{
		MWprintf (10, "Worker ID cannot be NULL in removeWorker\n" );
		return -1;
		}

		char buf[_POSIX_PATH_MAX];

		int condorID = fileWorkers[w->get_id1()].condorID;
		int procID = fileWorkers[w->get_id1()].condorprocID;
		MWprintf ( 10, "Removing workers %d.%d of %d\n", condorID, procID, w->get_id1() );
		sprintf ( buf, "%s/bin/condor_rm %d.%d", CONDOR_DIR, condorID, procID );

		if ( system ( buf ) < 0 )
		{
			MWprintf ( 10, "Could not condor_rm the worker\n" );
		}

		if ( fileWorkers[w->get_id1()].state == FILE_SUBMITTED )
		{
	//		submitted_num_workers--;

			hostadd_reqs[w->get_exec_class()]--;
	/*
		//XXX HOSTADDDEBUG
		MWprintf( 20, "Subtracting from %d hostadd_reqs.  Value %d\n", w->get_arch(),
			  hostadd_reqs[w->get_arch()] );
			for ( int i = 0; i < 2 * HOSTINC; i++ )
				if ( hostaddind_reqs[w->get_arch()][i] == w->get_id1() ) hostaddind_reqs[w->get_arch()][i] = -1;
	*/
		}
	/*
		else
			current_num_workers--;
	*/
		fileWorkers[w->get_id1()].state = FILE_FREE;
		write_RMstate ( NULL );
		return 0;
	}


//---------------------------------------------------------------------
// MWFileRC::do_spawn(  int workerid ):
//	Spawns the worker.
//---------------------------------------------------------------------
int
MWFileRC::do_spawn ( int nWorkers, int ex_cl )
{
	int cID = -1;
	int i;
	int ids[nWorkers];
	char sub_file[_POSIX_PATH_MAX];
	char exe[_POSIX_PATH_MAX];
	char logfile[_POSIX_PATH_MAX];
	char master_waitfile[_POSIX_PATH_MAX], worker_waitfile[_POSIX_PATH_MAX];
	FILE *ptr;
	char temp[_POSIX_PATH_MAX];
	char requirements[num_executables * _POSIX_PATH_MAX];

	if ( nWorkers <= 0 ) return 0;

	sprintf( sub_file, "%s/submit_file.%d", control_directory, subId++ );
	sprintf( exe, "%s/bin/condor_submit %s", CONDOR_DIR, sub_file );

	FILE *f = Open( sub_file, "w" );
	if ( f == NULL )
	{
		MWprintf ( 10, "Couldn't open the submitfile for writing\n");
		return -1;
	}
    	
	int index = 0;
	for ( i = 0; i < target_num_workers && index < nWorkers; i++ ) 
	{
		if ( fileWorkers[i].state == FILE_FREE ) 
		{
			ids[index++] = i;
		}
	}

	if ( index < nWorkers ) 
	{
		MWprintf ( 10, "In MW-File hostaddlogic asking for more workers than target_num_workers\n");
		MWprintf ( 10, "Had asked %d workers but I am adding only %d\n", nWorkers, index );
		nWorkers = index;
	}

#ifdef USE_CHIRP
	fprintf ( f, "Universe = Vanilla\n");
	fprintf( f, "+WantIOProxy=True\n");
	fprintf( f, "should_transfer_files = Yes\n"); // These two just needed to move to the execute dir
	fprintf( f, "when_to_transfer_output = ON_EXIT\n");
#else
	fprintf ( f, "Universe = Standard\n");
#endif
	fprintf( f, "Executable = mw_exec%d.$$(Opsys).$$(Arch)\n", ex_cl );

	bool tempfirst = TRUE;
	for ( i = 0; i < num_executables; i++ )
	{
		if ( worker_executables[i]->exec_class == ex_cl )
		{
			if ( tempfirst == TRUE )
			{
				sprintf ( requirements, "( ( %s ) ", arch_class_attributes[worker_executables[i]->arch_class] );
				tempfirst = FALSE;
			}
			else
			{
				sprintf ( requirements, "%s || ( %s ) ", requirements, arch_class_attributes[worker_executables[i]->arch_class] );
			}
		}
	}
	strcat ( requirements, " )" );
	for ( i = 0; i < nWorkers; i++ ) 
	{
		sprintf( logfile, "%s/log_file.%d", control_directory, ids[i] );
		sprintf( master_waitfile, "./%s/master_waitfile.%d", 
							output_directory, ids[i] );
		sprintf( worker_waitfile, "./%s/worker_waitfile.%d", 
							input_directory, ids[i] );
    
		fprintf( f, "arguments = %d %d %d %s %s -1\n", ids[i], 0, 0,
							output_directory, input_directory ); 
		fprintf( f, "log = %s/log_file.%d\n", control_directory, ids[i] );

		fprintf( f, "Output = %s/output_file.%d.$(Cluster).$(Process)\n", control_directory, ids[i] );
		fprintf( f, "Error = %s/error_file.%d.$(Cluster).$(Process)\n", control_directory, ids[i] );

		fprintf( f, "Requirements = %s\n ", requirements );
		fprintf ( f, "getenv = True\n" );

		// Jeff changed this

		FILE *userf = Open( "worker_attributes", "r" );
		if( userf == NULL ) 
		{
			MWprintf( 10, "No worker_attributes file, assuming Condor defaults\n" );
		}
		else 
		{
			int c;
			while( ( c = getc( userf ) ) != EOF ) 
			{
				putc( c, f );
			}
			Close( userf );
		}

		fprintf( f, "Queue\n" );
    
		HardDeleteFile( logfile );
		HardDeleteFile( master_waitfile );
		HardDeleteFile( worker_waitfile );

	}
	Close( f );

	MWprintf (10, "About to call condor_submit %s\n", exe);
	if ( (ptr = popen(exe, "r") ) != NULL ) 
	{
		while ( fscanf(ptr, "%s", temp ) >= 0 ) 
		{

			if ( strcmp ( temp, "cluster" ) == 0 ) 
			{

				// This is for Jeff's debugging with JP
				//MWprintf( 10, "Got cluster\n" );

				fscanf( ptr, "%s", temp );

				// This is for Jeff's debugging with JP
				//MWprintf( 10, "Read %s from the pipe.\n", temp );
	
				cID = atoi ( temp );
				pclose(ptr);
				ptr = NULL;
				MWprintf ( 10, "Spawned to cluster %d\n", cID );
				break;
			}
		}
		if ( ptr ) pclose ( ptr ); 
		for ( i = 0; i < nWorkers; i++ ) 
		{
			int tempp = ids[i];
			fileWorkers[tempp].condorID = cID;
			fileWorkers[tempp].condorprocID = i;
			fileWorkers[tempp].arch = -1;
			fileWorkers[tempp].counter = 0;
			fileWorkers[tempp].worker_counter = 0;
			fileWorkers[tempp].id = tempp;
			fileWorkers[tempp].state = FILE_SUBMITTED;
			fileWorkers[tempp].event_no = 1;
			fileWorkers[tempp].exec_class = ex_cl;
/*
			for ( int j = 0; j < 2 * HOSTINC; j++ ) {
				if ( hostaddind_reqs[archid][j] == -1 ) { 
				hostaddind_reqs[archid][j] = tempp; 
				break; 
				}
			}
*/
		}
	}
	else 
	{
		MWprintf ( 10, "Couldn't popen in FileRC\n");
		return -1;
	}

	write_RMstate ( NULL );
	return nWorkers;
}

int
MWFileRC::initsend ( int useless )
{
	/* This prepares for a series of packs before we do the actual send */
	if ( !sendList )
		sendList = new MWList<void> ;

	while ( sendList->number() > 0 )
	{
		struct FileSendBuffer *buf = (struct FileSendBuffer *)sendList->Remove();
		delete buf;
	}

	checksum = 0;
	return 0;
}

void
MWFileRC::who ( int *wh )
{
	unpack ( wh, 1, 1 );
}


//---------------------------------------------------------------------------
//
// 
//---------------------------------------------------------------------------
int
MWFileRC::send ( int to_whom, int msgtag )
{
	char filename[_POSIX_PATH_MAX];
	char control_file[_POSIX_PATH_MAX];
	FILE *fp;
	if ( isMaster == TRUE )
	{
		// In the master mode 
		sprintf ( filename, "%s/worker_input.%d.%d", input_directory, to_whom,
					fileWorkers[to_whom].worker_counter );
		sprintf ( control_file, "%s/worker_waitfile.%d", input_directory, 
					to_whom );
	}
	else
	{
		// We are sending something to the master. So to_whom actually has
		// no meaning.
		sprintf ( filename, "%s/worker_output.%d.%d", output_directory,
					FileRCID, master_expected_number );
		sprintf ( control_file, "%s/master_waitfile.%d", output_directory,
					FileRCID );
	}

	fp = Open ( filename, "w" );
	if ( fp == NULL )
	{
		MWprintf ( 10, "Couldn't open file %s for writing work.  Bailing out!\n", filename );

		assert( filename && 0 );
		return -1;
	}

	fprintf ( fp, "%d %lld ", msgtag, checksum );

	while ( sendList->IsEmpty() != TRUE )
	{
		struct FileSendBuffer *buf = ( struct FileSendBuffer *)sendList->Remove();
		switch ( buf->type )
		{
			case INT:
				fprintf ( fp, "%d %d ", INT, *(int *)buf->data );
				break;
			case CHAR:
				fprintf ( fp, "%d %c ", CHAR, *(char *)buf->data );
				break;
			case LONG:
				fprintf ( fp, "%d %ld ", LONG, *(long *)buf->data );
				break;
			case FLOAT:
				fprintf ( fp, "%d %f ", FLOAT, *(float *)buf->data );
				break;
			case DOUBLE:
				fprintf ( fp, "%d %f ", DOUBLE, *(double*)buf->data );
				break;
			case UNSIGNED_INT:
				fprintf ( fp, "%d %o ", UNSIGNED_INT, *(unsigned int *)buf->data );
				break;
			case SHORT:
				fprintf ( fp, "%d %hd ", SHORT, *(short *)buf->data );
				break;
			case UNSIGNED_SHORT:
				fprintf ( fp, "%d %ho ", UNSIGNED_SHORT, *(unsigned short *)buf->data );
				break;
			case UNSIGNED_LONG:
				fprintf ( fp, "%d %lo ", UNSIGNED_LONG, *(unsigned long *)buf->data );
				break;
			case STRING:
				fprintf ( fp, "%d %d %s ", STRING, buf->size - 1, 
						(char *)buf->data );
				break;
			default:
				MWprintf ( 10, "Couldnot decipher a data type while sending\n");
				MWprintf ( 10, "Got %d\n", buf->type );
				Close ( fp );
				return -1;
		}

		delete buf;
	}

	fflush ( fp );
	Close ( fp );

	if ( isMaster || ( !isMaster && master_expected_number == 0 ) )
	{
		fp = Open ( control_file, "a" );
	}
	else
	{
		// You are a slave and expected_number is not zero.
		fp = Open ( control_file, "a" );
		if (fp == NULL)
		{
			MWprintf ( 10, "In slave my %s doesn't exist\n");
			HardDeleteFile ( filename );
			return 0;
		}
	}
	if ( fp == NULL )
	{
		MWprintf ( 10, "Could not open the waitfile %s.  Send will fail!\n", control_file );

		assert( control_file && 0 );
		return -1;
	}

	if ( isMaster == TRUE )
	{
		fprintf ( fp, "%d %c", fileWorkers[to_whom].worker_counter, OK );
		fileWorkers[to_whom].worker_counter++;
	}
	else
	{
		fprintf ( fp, "%d %c",  master_expected_number, OK );
		master_expected_number++;
	}
	fflush ( fp );
	Close ( fp );

#ifdef USE_CHIRP
	if (!isMaster) {
		sendFileToMaster(filename);
		sendFileToMaster(control_file);
	}
#endif
	return 0;
}


int
MWFileRC::recv ( int from_whom, int msgtag )
{
	if ( !recvList )
		recvList = new MWList<void> ;
	while ( recvList->number() > 0 )
	{
		struct FileSendBuffer *buf = (struct FileSendBuffer *)recvList->Remove();
		delete buf;
	}

	if ( isMaster == TRUE )
		return master_recv ( from_whom, msgtag );
	else
		return worker_recv ( from_whom, msgtag );
}

int
MWFileRC::nrecv ( int from_whom, int msgtag )
{
  if ( !recvList )
    recvList = new MWList<void> ;
   
  while ( recvList->number() > 0 )
    {
      struct FileSendBuffer *buf = (struct FileSendBuffer *)recvList->Remove();
      delete buf;
    }

    //from_whom has no meaning as it is invariably from the master.
    static time_t last_comm = time(0);
    int calculated_number = 0;
    char c;
    char filename[_POSIX_PATH_MAX];
    FILE *op = NULL;
    bool toFlush = FALSE;
    time_t starttime;
    time_t waitedtime;
    starttime = waitedtime = time(0);

    sprintf ( filename, "%s/worker_waitfile.%d", input_directory, FileRCID );
    op = Open ( filename, "r" );
    while ( 1 )
      {
      waitedtime = time(0);
      /*
      if(waitedtime - starttime > SLEEP_DELAY*2)
      {
        Close(op);
        msgTag = NO_MESSAGE;
    return 0;
      }
      */
    if ( worker_timeout > 0 && last_comm + worker_timeout < time(0) )
      {
        MWprintf ( 10, "Worker Timed out\n");
        if ( op ) Close ( op );
        return WORKER_TIMEOUT;
      }
    if ( !op )
      {
        MWSystem::sleep( SLEEP_DELAY );
        op = Open ( filename, "r" );
        //MWprintf(10,"Cannot open file, try again\n");
        continue;
      }
    if ( fscanf( op, "%c", &c ) < 0 )
      {// end of file, return from non-blocking receive
            Close(op);
            msgTag = NO_MESSAGE;
        return 0;
      /*
        Close( op );
        sleep( SLEEP_DELAY );
        op = Open ( filename, "r" );
        MWprintf(10,"opened file no data, try again\n"); // wenhan
        continue;
        */
      }
    if ( c == OK )
      {
        if ( calculated_number != expected_number )
          {
        if ( calculated_number > expected_number )
          {
            MWprintf(10, "Fault in master resuscicate logic,exiting\n");
            Close ( op );
            return UNABLE_TO_WAKE;
          }
        calculated_number = 0;
        continue;
          }
             else
        {
                MWprintf ( 10, "Got the work in slave %d\n", expected_number );
                expected_number++;
                Close(op);
                whomRecv = 0;
                last_comm = time(0);
                int hh = handle_work( msgtag );
                if ( toFlush )
                  {
                expected_number = 0;
                master_expected_number = 0;
                HardDeleteFile ( filename );
                toFlush = FALSE;
                  }
                return hh;
              }
      }
    else if ( c == ' ' )
      {
      }
    else if ( c == SYNC)
      {
        // Master wants it to ignore some things.
        calculated_number = 0;
      }
    else if ( c == ESYNC )
      {
        expected_number = calculated_number;
        calculated_number = 0;
        toFlush = TRUE;
      }
    else if ( isdigit ( c ) )
      {
        calculated_number = calculated_number * 10 + c - '0';
      }
                else
          {
            Close ( op );
            return WAITFILE_PROTOCOL_ERROR;
          }
      }
}


int MWFileRC::master_recv ( int from_whom, int msgtag )
{
	// msgtag has no meaning as we have to recv all sorts of control info.
    int i;
    bool change = FALSE;
    while ( 1 )
    {
    	if ( turnNo == CHECKLOG_FREQ )
    	{
	    CheckLogFilesRunning();
	    turnNo = 0;
    	}
    	for ( i = cyclePosition; i < target_num_workers; i++ )
    	{
	    cyclePosition++;
	    switch ( fileWorkers[i].state )
	    {
	    	case FILE_RUNNING:
		    change = IsComplete ( i );
		    if ( change == TRUE )
		    {
    			MWprintf ( 10, "Master received something from worker %d\n", i );
			return handle_finished_worker ( i );
	    }
		    break;

		case FILE_EXECUTE:
		    if ( fileWorkers[i].served != FILE_RUNNING )
		    {
    		    	MWprintf ( 10, "Master found that worker is Started %d\n", i );
		    	return handle_executing_worker (i );
		    }
		    break;

	    	case FILE_KILLED:
		    if ( fileWorkers[i].served != FILE_KILLED )
		    {
		    	MWprintf ( 10, "Master found that worker was Killed %d\n", i);
		    	return handle_killed_worker(i);
		    }
		    break;
		
		case FILE_SUSPENDED:
		    if ( fileWorkers[i].served != FILE_SUSPENDED )
		    {
		    	MWprintf ( 10, "Master found that worker was Suspended %d \n", i);
		    	return handle_suspended_worker(i);
		    }
		    break;

		case FILE_RESUMED:
		    if ( fileWorkers[i].served != FILE_RESUMED )
		    {
		    	MWprintf ( 10, "Master found that worker was Resumed %d\n", i);
		    	return handle_resumed_worker(i);
		    }
		    break;

		case FILE_TRANSIT:
		    if ( fileWorkers[i].served != FILE_TRANSIT )
		    {
		    	MWprintf ( 10, "Master found that worker transited %d \n", i );
			return handle_transited_worker(i);
		    }

		default:
		    break;
	    }

	}

	cyclePosition = 0;
	turnNo++;
    }
}


int
MWFileRC::worker_recv ( int from_whom, int msgtag )
{
	//from_whom has no meaning as it is invariably from the master.
	static time_t last_comm = time(0);
	int calculated_number = 0;
	char c;
	char filename[_POSIX_PATH_MAX];
	FILE *op = NULL;
	bool toFlush = FALSE;


	sprintf ( filename, "%s/worker_waitfile.%d", input_directory, FileRCID );
#ifdef USE_CHIRP
	getFileFromMaster(filename);
#endif
	op = Open ( filename, "r" );
	while ( 1 )
	{
		if ( worker_timeout > 0 && last_comm + worker_timeout < time(0) )
		{
			MWprintf ( 10, "Worker Timed out\n");
			if ( op ) Close ( op );
			return WORKER_TIMEOUT;
		}
		if ( !op )
		{
			MWSystem::sleep( SLEEP_DELAY );
#ifdef USE_CHIRP
			getFileFromMaster(filename);
#endif
			op = Open ( filename, "r" );
			continue;
		}
		if ( fscanf( op, "%c", &c ) < 0 )
		{
			Close( op );
			MWSystem::sleep( SLEEP_DELAY );
#ifdef USE_CHIRP
			getFileFromMaster(filename);
#endif
			op = Open ( filename, "r" );
			continue;
		}
		if ( c == OK )
		{
			if ( calculated_number != expected_number )
			{
				if ( calculated_number > expected_number )
				{
					MWprintf(10, "Fault in master resuscicate logic,exiting\n");
					Close ( op );
					return UNABLE_TO_WAKE;
				}
				calculated_number = 0;
				continue;
			}
			else
			{
				MWprintf ( 10, "Got the work in slave %d\n", expected_number );
				expected_number++;
				Close(op);
				whomRecv = 0;
				last_comm = time(0);
				int hh = handle_work( msgtag );
				if ( toFlush )
				{
					expected_number = 0;
					master_expected_number = 0;
					HardDeleteFile ( filename );
					toFlush = FALSE;
				}
				return hh;
			}
		}
		else if ( c == ' ' )
		{
		}
		else if ( c == SYNC)
		{
			// Master wants it to ignore some things.
			calculated_number = 0;
		}
		else if ( c == ESYNC )
		{
			expected_number = calculated_number;
			calculated_number = 0;
			toFlush = TRUE;
		}
		else if ( isdigit ( c ) )
		{
			calculated_number = calculated_number * 10 + c - '0';
		}
		else
		{
			Close ( op );
			return WAITFILE_PROTOCOL_ERROR;
		}
	}
}


int
MWFileRC::pack ( const char *bytes, int nitem, int stride )
{
    int i;

    if ( nitem < 0 ) return -1;
    if ( nitem == 0 ) return 0;

    if ( sendList == NULL )
    {
	MWprintf ( 10, "Didn't do an init send before packing something\n");
	return -1;
    }

    for ( i = 0; i < nitem; i++ )
    {
	FileSendBuffer *buf = 
		new FileSendBuffer ( (void *)&bytes[i * stride], CHAR, sizeof(char) );
	sendList->Append ( (void *)buf );
	checksum += CHAR;
    }

    return 0;
}


int
MWFileRC::pack ( const float *f, int nitem, int stride )
{
    int i;

    if ( nitem < 0 ) return -1;
    if ( nitem == 0 ) return 0;

    if ( sendList == NULL )
    {
	MWprintf ( 10, "Didn't do an init send before packing something\n");
	return -1;
    }

    for ( i = 0; i < nitem; i++ )
    {
	FileSendBuffer *buf = 
		new FileSendBuffer ( (void *)&(f[i * stride]), FLOAT, sizeof(float));
	sendList->Append ( (void *)buf );
	checksum += FLOAT;
    }
    return 0;
}


int
MWFileRC::pack ( const double *d, int nitem, int stride )
{
    int i;

    if ( nitem < 0 ) return -1;
    if ( nitem == 0 ) return 0;

    if ( sendList == NULL )
    {
	MWprintf ( 10, "Didn't do an init send before packing something\n");
	return -1;
    }

    for ( i = 0; i < nitem; i++ )
    {
	FileSendBuffer *buf = 
		new FileSendBuffer ( (void *)&(d[i * stride]), DOUBLE, sizeof(double) );
	sendList->Append ( (void *)buf );
	checksum += DOUBLE;
    }
    return 0;
}


int
MWFileRC::pack ( const int *j, int nitem, int stride )
{
    int i;

    if ( nitem < 0 ) return -1;
    if ( nitem == 0 ) return 0;
    
    if ( sendList == NULL )
    {
	MWprintf ( 10, "Didn't do an init send before packing something\n");
	return -1;
    }
    for ( i = 0; i < nitem; i++ )
    {
	FileSendBuffer *buf = 
		new FileSendBuffer ( (void *)&(j[i * stride]), INT, sizeof(int) );
	sendList->Append ( (void *)buf );
	checksum += INT;
    }
    return 0;
}


int
MWFileRC::pack ( const unsigned int *ui, int nitem, int stride )
{
    int i;

    if ( nitem < 0 ) return -1;
    if ( nitem == 0 ) return 0;

    if ( sendList == NULL )
    {
	MWprintf ( 10, "Didn't do an init send before packing something\n");
	return -1;
    }

    for ( i = 0; i < nitem; i++ )
    {
	FileSendBuffer *buf = new FileSendBuffer ( (void *)&(ui[i * stride]), UNSIGNED_INT, sizeof(unsigned int) );
	sendList->Append ( (void *)buf );
	checksum += UNSIGNED_INT;
    }
    return 0;
}


int
MWFileRC::pack ( const short *sh, int nitem, int stride )
{
    int i;

    if ( nitem < 0 ) return -1;
    if ( nitem == 0 ) return 0;

    if ( sendList == NULL )
    {
	MWprintf ( 10, "Didn't do an init send before packing something\n");
	return -1;
    }

    for ( i = 0; i < nitem; i++ )
    {
	FileSendBuffer *buf = 
		new FileSendBuffer ( (void *)&(sh[i * stride]), SHORT, sizeof(short) );
	sendList->Append ( (void *)buf );
	checksum += SHORT;
    }
    return 0;
}


int
MWFileRC::pack ( const unsigned short *ush, int nitem, int stride )
{
    int i;

    if ( nitem < 0 ) return -1;
    if ( nitem == 0 ) return 0;

    if ( sendList == NULL )
    {
	MWprintf ( 10, "Didn't do an init send before packing something\n");
	return -1;
    }

    for ( i = 0; i < nitem; i++ )
    {
	FileSendBuffer *buf = 
		new FileSendBuffer ( (void *)&(ush[i * stride]), UNSIGNED_SHORT,
				sizeof(unsigned short) );
	sendList->Append ( (void *)buf );
	checksum += UNSIGNED_SHORT;
    }
    return 0;
}


int
MWFileRC::pack ( const long *l, int nitem, int stride )
{
    int i;

    if ( nitem < 0 ) return -1;
    if ( nitem == 0 ) return 0;

    if ( sendList == NULL )
    {
	MWprintf ( 10, "Didn't do an init send before packing something\n");
	return -1;
    }

    for ( i = 0; i < nitem; i++ )
    {
	FileSendBuffer *buf = 
		new FileSendBuffer ( (void *)&(l[i * stride]), LONG, sizeof(long) );
	sendList->Append ( (void *)buf );
	checksum += LONG;
    }
    return 0;
}


int
MWFileRC::pack ( const unsigned long *ul, int nitem, int stride )
{
    int i;

    if ( nitem < 0 ) return -1;
    if ( nitem == 0 ) return 0;

    if ( sendList == NULL )
    {
	MWprintf ( 10, "Didn't do an init send before packing something\n");
	return -1;
    }

    for ( i = 0; i < nitem; i++ )
    {
	FileSendBuffer *buf = 
		new FileSendBuffer ( (void *)&(ul[i * stride]), UNSIGNED_LONG, 
						sizeof(unsigned long) );
	sendList->Append ( (void *)buf );
	checksum += UNSIGNED_LONG;
    }
    return 0;
}


int
MWFileRC::pack ( const char *str )
{
    if ( str == NULL ) return -1;
    if ( strlen(str) == 0 ) return 0;

    if ( sendList == NULL )
    {
	MWprintf ( 10, "Didn't do an init send before packing something\n");
	return -1;
    }

    FileSendBuffer *buf = 
		new FileSendBuffer ( (void *)str, STRING, 
					strlen(str) + 1 );
    sendList->Append ( (void *)buf );
    checksum += STRING;

    return 0;
}


int
MWFileRC::unpack ( char *bytes, int nitems, int stride )
{
    int i;

    if ( nitems <= 0 ) return 0;

    if ( recvList == NULL || recvList->IsEmpty ( ) == TRUE )
    {
	MWprintf ( 10, "Trying to unpack before receiving \n" );
	return -1;
    }

    for ( i = 0; i < nitems; i++ )
    {
	FileSendBuffer *buf = (FileSendBuffer *)recvList->Remove();
	if ( buf == NULL || buf->type != CHAR )
	{
	    MWprintf ( 20, "Recieve buffer contained something else \n");
	    recvList->Prepend ( (void *)buf );
	    return -1;
	}

	memcpy ( &bytes[i * stride], buf->data, sizeof(char) );
	
	delete buf;
    }

    return 0;
}


int
MWFileRC::unpack ( float *f, int nitems, int stride )
{
    int i;

    if ( nitems <= 0 ) return -1;

    if ( recvList == NULL || recvList->IsEmpty ( ) == TRUE )
    {
	MWprintf ( 10, "Trying to unpack before receiving \n" );
	return -1;
    }

    for ( i = 0; i < nitems; i++ )
    {
	FileSendBuffer *buf = (FileSendBuffer *)recvList->Remove();
	if ( buf == NULL || buf->type != FLOAT )
	{
	    MWprintf ( 20, "Recieve buffer contained something else \n");
	    recvList->Prepend ( (void *)buf );
	    return -1;
	}

	memcpy ( &f[i * stride], buf->data, sizeof(float) );
	
	delete buf;
    }

    return 0;
}


int
MWFileRC::unpack ( double *d, int nitems, int stride )
{
    int i;

    if ( nitems <= 0 ) return -1;

    if ( recvList == NULL || recvList->IsEmpty ( ) == TRUE )
    {
	MWprintf ( 10, "Trying to unpack before receiving \n" );
	return -1;
    }

    for ( i = 0; i < nitems; i++ )
    {
	FileSendBuffer *buf = (FileSendBuffer *)recvList->Remove();
	if ( buf == NULL || buf->type != DOUBLE )
	{
	    MWprintf ( 20, "Recieve buffer contained something else \n");
	    recvList->Prepend ( (void *)buf );
	    return -1;
	}

	memcpy ( &d[i * stride], buf->data, sizeof(double) );
	
	delete buf;
    }

    return 0;
}


int
MWFileRC::unpack ( int *j, int nitems, int stride )
{
    int i;

    if ( nitems <= 0 ) return -1;

    if ( recvList == NULL || recvList->IsEmpty ( ) == TRUE )
    {
	MWprintf ( 10, "Trying to unpack before receiving \n" );
	return -1;
    }

    for ( i = 0; i < nitems; i++ )
    {
	FileSendBuffer *buf = (FileSendBuffer *)recvList->Remove();
	if ( buf == NULL || buf->type != INT )
	{
	    MWprintf ( 20, "Recieve buffer contained something else \n");
	    recvList->Prepend ( (void *)buf );
	    return -1;
	}
	memcpy ( &j[i * stride], buf->data, sizeof(int) );
	delete buf;
    }
    return 0;
}


int
MWFileRC::unpack ( unsigned int *ui, int nitems, int stride )
{
    int i;

    if ( nitems <= 0 ) return -1;

    if ( recvList == NULL || recvList->IsEmpty ( ) == TRUE )
    {
	MWprintf ( 10, "Trying to unpack before receiving \n" );
	return -1;
    }

    for ( i = 0; i < nitems; i++ )
    {
	FileSendBuffer *buf = (FileSendBuffer *)recvList->Remove();
	if ( buf == NULL || buf->type != UNSIGNED_INT )
	{
	    MWprintf ( 20, "Recieve buffer contained something else \n");
	    recvList->Prepend ( (void *)buf );
	    return -1;
	}

	memcpy ( &ui[i * stride], buf->data, sizeof(unsigned int) );
	
	delete buf;
    }
    return 0;
}


int
MWFileRC::unpack ( short *sh, int nitems, int stride )
{
    int i;

    if ( nitems <= 0 ) return -1;

    if ( recvList == NULL || recvList->IsEmpty ( ) == TRUE )
    {
	MWprintf ( 10, "Trying to unpack before receiving \n" );
	return -1;
    }

    for ( i = 0; i < nitems; i++ )
    {
	FileSendBuffer *buf = (FileSendBuffer *)recvList->Remove();
	if ( buf == NULL || buf->type != SHORT )
	{
	    MWprintf ( 20, "Recieve buffer contained something else \n");
	    recvList->Prepend ( (void *)buf );
	    return -1;
	}

	memcpy ( &sh[i * stride], buf->data, sizeof(short) );
	
	delete buf;
    }

    return 0;
}

int
MWFileRC::unpack ( unsigned short *ush, int nitems, int stride )
{
    int i;

    if ( nitems <= 0 ) return -1;

    if ( recvList == NULL || recvList->IsEmpty ( ) == TRUE )
    {
	MWprintf ( 10, "Trying to unpack before receiving \n" );
	return -1;
    }

    for ( i = 0; i < nitems; i++ )
    {
	FileSendBuffer *buf = (FileSendBuffer *)recvList->Remove();
	if ( buf == NULL || buf->type != UNSIGNED_SHORT )
	{
	    MWprintf ( 20, "Recieve buffer contained something else \n");
	    recvList->Prepend ( (void *)buf );
	    return -1;
	}

	memcpy ( &ush[i * stride], buf->data, sizeof(unsigned short) );
	
	delete buf;
    }

    return 0;
}


int
MWFileRC::unpack ( long *l, int nitems, int stride )
{
    int i;

    if ( nitems <= 0 ) return -1;

    if ( recvList == NULL || recvList->IsEmpty ( ) == TRUE )
    {
	MWprintf ( 10, "Trying to unpack before receiving \n" );
	return -1;
    }

    for ( i = 0; i < nitems; i++ )
    {
	FileSendBuffer *buf = (FileSendBuffer *)recvList->Remove();
	if ( buf == NULL || buf->type != LONG )
	{
	    MWprintf ( 20, "Recieve buffer contained something else \n");
	    recvList->Prepend ( (void *)buf );
	    return -1;
	}

	memcpy ( &l[i * stride], buf->data, sizeof(long) );
	
	delete buf;
    }

    return 0;
}


int
MWFileRC::unpack ( unsigned long *ul, int nitems, int stride )
{
    int i;

    if ( nitems <= 0 ) return -1;

    if ( recvList == NULL || recvList->IsEmpty ( ) == TRUE )
    {
	MWprintf ( 10, "Trying to unpack before receiving \n" );
	return -1;
    }

    for ( i = 0; i < nitems; i++ )
    {
	FileSendBuffer *buf = (FileSendBuffer *)recvList->Remove();
	if ( buf == NULL || buf->type != UNSIGNED_LONG )
	{
	    MWprintf ( 20, "Recieve buffer contained something else \n");
	    recvList->Prepend ( (void *)buf );
	    return -1;
	}

	memcpy ( &ul[i * stride], buf->data, sizeof(unsigned long) );
	
	delete buf;
    }

    return 0;
}


int
MWFileRC::unpack ( char *str )
{
    if ( recvList == NULL || recvList->IsEmpty ( ) == TRUE )
    {
	MWprintf ( 10, "Trying to unpack before receiving \n" );
	return -1;
    }

    FileSendBuffer *buf = (FileSendBuffer *)recvList->Remove();
    if ( buf == NULL || buf->type != STRING )
    {
    	MWprintf ( 20, "Recieve buffer contained something else %d\n", buf->type);
    	recvList->Prepend ( (void *)buf );
    	return -1;
    }

    memcpy ( str, buf->data, sizeof(char) * buf->size );

    delete buf;

    return 0;
}


bool 
MWFileRC::IsComplete ( int i )
{
    char c;
    char f_name[_POSIX_PATH_MAX];

    int number_calculated = 0;
    int num = fileWorkers[i].counter;

    sprintf( f_name, "./%s/master_waitfile.%d", output_directory, i );

    FILE *fp = Open ( f_name, "r" );
    if ( fp == NULL )
    {
	return FALSE;
    }

    while ( fscanf(fp, "%c", &c ) > 0 )
    {
	if ( c == OK )
	{
	    if ( number_calculated == num )
	    {
		Close(fp);
		return TRUE;
	    }
	    else
		number_calculated = 0;
	}
	else if ( c == WRONGINITFILE )
	{
		// The process found a wrong init file.
		// It probably has killed itself. Thus do nothing.
	}
	else if ( c == SYNC )
	{
	    number_calculated = 0;
	}
	else if ( isdigit(c) > 0 )
	{
	    number_calculated = number_calculated * 10 + c - '0';
	}
    }

    Close(fp);
    return FALSE;               
}


int
MWFileRC::handle_finished_worker ( int i )
{
    FILE *fp;
    char filename[_POSIX_PATH_MAX];
    int typeno;

    fileWorkers[i].served = FILE_RUNNING;
    sprintf ( filename, "%s/worker_output.%d.%d", output_directory, fileWorkers[i].id, fileWorkers[i].counter );   
    fp = Open ( filename, "r" );
    if ( fp == NULL )
    {
	MWprintf ( 10, "Couldn't open the worker output file %s for reading: Errno is %d\n", filename, errno );
	return -1;
    }
    fscanf ( fp, "%d %lld ", &msgTag, &checksum );

    while ( fscanf ( fp, "%d ", &typeno ) >= 0 )
    {
	FileSendBuffer *buf;

	switch ( typeno )
	{
	    case INT:
		int it1;
		fscanf ( fp, "%d ", &it1 );
		buf = new FileSendBuffer ( (void *)&it1, INT, sizeof(int) );
		if (buf == NULL)
			MWprintf(30, "MWDEBUG: FileSendBuffer construction failed\n");
		recvList->Append ( (void *)buf );
		checksum -= INT;
		break;

	    case CHAR:
		char it2;
		fscanf ( fp, "%c ", &it2 );
		buf = new FileSendBuffer ( (void *)&it2, CHAR, sizeof(char) );
		recvList->Append ( (void *)buf );
		checksum -= CHAR;
		break;

	    case LONG:
		long it3;
		fscanf ( fp, "%ld ", &it3 );
		buf = new FileSendBuffer ( (void *)&it3, LONG, sizeof(long) );
		recvList->Append ( (void *)buf );
		checksum -= LONG;
		break;

	    case FLOAT:
		float it4;
		fscanf ( fp, "%f ", &it4 );
		buf = new FileSendBuffer ( (void *)&it4, FLOAT, sizeof(float) );
		recvList->Append ( (void *)buf );
		checksum -= FLOAT;
		break;

	    case DOUBLE:
		double it5;
		fscanf ( fp, "%lf ", &it5 );
		buf = new FileSendBuffer ( (void *)&it5, DOUBLE, sizeof(double));
		if (buf == NULL)
			MWprintf(30, "MWDEBUG: FileSendBuffer construction failed\n");
		recvList->Append ( (void *)buf );
		checksum -= DOUBLE;
		break;

	    case UNSIGNED_INT:
		unsigned int it6;
		fscanf ( fp, "%o ", &it6 );
		buf = new FileSendBuffer ( (void *)&it6, UNSIGNED_INT, 
						sizeof(unsigned int));
		recvList->Append ( (void *)buf );
		checksum -= UNSIGNED_INT;
		break;

	    case SHORT:
		short it7;
		fscanf ( fp, "%hd ", &it7 );
		buf = new FileSendBuffer ( (void *)&it7, SHORT, sizeof(short));
		recvList->Append ( (void *)buf );
		checksum -= SHORT;
		break;

	    case UNSIGNED_SHORT:
		unsigned short it8;
		fscanf ( fp, "%ho ", &it8 );
		buf = new FileSendBuffer ( (void *)&it8, UNSIGNED_SHORT, 
						sizeof(unsigned short));
		recvList->Append ( (void *)buf );
		checksum -= UNSIGNED_SHORT;
		break;

	    case UNSIGNED_LONG:
		unsigned long it9;
		fscanf ( fp, "%lo ", &it9 );
		buf = new FileSendBuffer ( (void *)&it9, UNSIGNED_LONG, 
						sizeof(unsigned long));
		recvList->Append ( (void *)buf );
		checksum -= UNSIGNED_LONG;
		break;

	    case STRING:
		{
		int sz, z;
		fscanf ( fp, "%d ", &sz );
		char *it10 = new char[sz+1];
		for ( z = 0; z < sz; z++ )
		    fscanf ( fp, "%c", &it10[z] );
		it10[z] = '\0';
		buf = new FileSendBuffer ( (void *)it10, STRING, 
					strlen(it10) + 1);
		recvList->Append ( (void *)buf );
		delete []it10;
		checksum -= STRING;
		break;
		}

	    default:
		MWprintf (10, "Couldn't decipher the data type in receiving\n");
		MWprintf (10, "Got %d\n", typeno );
		Close ( fp );
		return -1;
		break;
	}
    }
    if ( checksum != 0 )
    {
    	// Problem 
	MWprintf ( 10, "The checksum is not zero in recv\n");
	while ( recvList->IsEmpty() != TRUE )
	{
		FileSendBuffer *buff = (FileSendBuffer *)recvList->Remove();
		delete buff;
	}

	FileSendBuffer *newbuf = new FileSendBuffer ( &fileWorkers[i].id, INT,
				sizeof(int) );
	recvList->Append ( (void *)newbuf );

	msgTag = workerEvents[FileChecksumError];
    }
    
    whomRecv = i;
    Close ( fp );
    fileWorkers[i].counter++;
    HardDeleteFile ( filename );
    return 0;
}


int
MWFileRC::bufinfo ( int buf_id, int *len, int *tag, int *sending_host )
{
    *sending_host = whomRecv;
    *tag = msgTag;
    *sending_host = whomRecv;
    return 0;
}

int
MWFileRC::handle_killed_worker ( int i )
{
	// We got the info that a worker has been killed.
	// Inform the upper layer that the worker has disappeared.
    fileWorkers[i].served = FILE_KILLED;

    FileSendBuffer *buf = new FileSendBuffer ( &fileWorkers[i].id, INT,
				sizeof(int) );
    recvList->Append ( (void *)buf );

    whomRecv = i;
    msgTag = workerEvents[FileTaskExit];

    // Jeff debugging 
    MWprintf( 10, "State of killed worker was %d\n",  fileWorkers[i].state );

    if ( fileWorkers[i].state == FILE_SUBMITTED )
    {

      //XXX HOSTADDDEBUG

    	hostadd_reqs[fileWorkers[i].exec_class]--;
/*
	MWprintf( 20, "Substracting from %d hostadd_reqs.  Value %d\n", fileWorkers[i].arch,
		  hostadd_reqs[fileWorkers[i].arch] );

	submitted_num_workers--;
    	for ( int j = 0; j < 2 * HOSTINC; j++ )
    	    if ( hostaddind_reqs[fileWorkers[i].arch][j] == i ) hostaddind_reqs[fileWorkers[i].arch][j] = -1;
*/
    }
/*
    else
    	current_num_workers--;
*/

    fileWorkers[i].state = FILE_FREE;
    return 0;
}

int
MWFileRC::handle_transited_worker ( int i )
{
    // int j;

    MWprintf( 10, "Transit event was from worker %d, whose previous state was %d\n", 
	      i,  fileWorkers[whomRecv].served );

    fileWorkers[i].served = FILE_TRANSIT;
    FileSendBuffer *buf = new FileSendBuffer ( &fileWorkers[i].id, INT,
				sizeof(int) );
    recvList->Append ( (void *)buf );

    whomRecv = i;
    msgTag = workerEvents[FileTaskExit];

    if ( fileWorkers[i].state == FILE_SUBMITTED )
    {
    	MWprintf ( 10, "How can a transit event happen when I am in submit state\n");
		hostadd_reqs[fileWorkers[i].exec_class]--;
    }

	char exe[_POSIX_PATH_MAX];
	sprintf ( exe, "%s/bin/condor_rm %d.%d", CONDOR_DIR, fileWorkers[i].condorID, fileWorkers[i].condorprocID );
	system ( exe );
	fileWorkers[i].state = FILE_FREE;

    return 0;
}


int
MWFileRC::handle_suspended_worker ( int i )
{
	// We got the info from the log file that the worker has been suspended.
	// Inform the upper layer that this is the case.
	// No need to change the state of fileWorkers as that has been
	// modified by the CheckLogFiles.

    fileWorkers[i].served = FILE_SUSPENDED;
    FileSendBuffer *buf = new FileSendBuffer ( &fileWorkers[i].id, INT,
				sizeof(int) );
    recvList->Append ( (void *)buf );

    whomRecv = i;
    msgTag = workerEvents[FileTaskSuspend];

    return 0;
}



int
MWFileRC::handle_resumed_worker ( int i )
{
	// We got the info from the log file that the worker has been resumed.
	// Inform the upper layer that this is the case.
	// No need to change the state of fileWorkers as that has been
	// modified by the CheckLogFiles.

    fileWorkers[i].state = FILE_RUNNING;
    fileWorkers[i].served = FILE_RESUMED;
    FileSendBuffer *buf = new FileSendBuffer ( &fileWorkers[i].id, INT,
				sizeof(int) );
    recvList->Append ( (void *)buf );

    whomRecv = i;
    msgTag = workerEvents[FileTaskResume];

    return 0;
}

int
MWFileRC::handle_executing_worker ( int i )
{
    // We got the info that the worker i started executin.

    if ( !(&fileWorkers[i]) )
    {
	MWprintf( 10, "How can the corresponding entry in execute worker be null\n");
	return -1;
    }

    fileWorkers[i].served = FILE_EXECUTE;
    fileWorkers[i].state = FILE_RUNNING;
    FileSendBuffer *buf = new FileSendBuffer ( &fileWorkers[i].id, INT,
				sizeof(int) );
    recvList->Append ( (void *)buf );

    whomRecv = fileWorkers[i].id;
    msgTag = workerEvents[FileHostAdd];

    return 0;
}


int
MWFileRC::handle_work ( int msgtag )
{
	// Called by the worker recv. 
	// Read the worker_input file and prepare the recvList.
    FILE *fp;
    char filename[_POSIX_PATH_MAX];
    int typeno;

    sprintf ( filename, "%s/worker_input.%d.%d", input_directory, FileRCID,
				expected_number - 1 );   
#ifdef USE_CHIRP
	getFileFromMaster(filename);
#endif
    fp = Open ( filename, "r" );
    if ( fp == NULL )
    {
	MWprintf ( 10, "Couldn't open the worker input file %s for reading: errno = %d\n", filename, errno );
	return CANNOT_OPEN_INPUT_FILE;
    }

    if ( fscanf ( fp, "%d %lld ", &msgTag, &checksum ) <= 0 ) 
    {
    	Close ( fp );
    	return SCANF_ERROR_ON_INPUT_FILE;
    }

    if ( msgtag != -1 && msgTag != msgtag ) 
    {
	MWprintf( 10, "FileRC wrong message sequence got %d \n", msgTag);
	Close ( fp );
	return MESSAGE_SEQUENCE_ERROR;
    }

    while ( fscanf ( fp, "%d ", &typeno ) > 0 )
    {
	FileSendBuffer *buf;

	switch ( typeno )
	{
	    case INT:
		int it1;
		fscanf ( fp, "%d ", &it1 );
		buf = new FileSendBuffer ( (void *)&it1, INT, sizeof(int) );
		recvList->Append ( (void *)buf );
		checksum -= INT;
		break;

	    case CHAR:
		char it2;
		fscanf ( fp, "%c ", &it2 );
		buf = new FileSendBuffer ( (void *)&it2, CHAR, sizeof(char) );
		recvList->Append ( (void *)buf );
		checksum -= CHAR;
		break;

	    case LONG:
		long it3;
		fscanf ( fp, "%ld ", &it3 );
		buf = new FileSendBuffer ( (void *)&it3, LONG, sizeof(long) );
		recvList->Append ( (void *)buf );
		checksum -= LONG;
		break;

	    case FLOAT:
		float it4;
		fscanf ( fp, "%f ", &it4 );
		buf = new FileSendBuffer ( (void *)&it4, FLOAT, sizeof(float) );
		recvList->Append ( (void *)buf );
		checksum -= FLOAT;
		break;

	    case DOUBLE:
		double it5;
		fscanf ( fp, "%lf ", &it5 );
		buf = new FileSendBuffer ( (void *)&it5, DOUBLE, sizeof(double));
		recvList->Append ( (void *)buf );
		checksum -= DOUBLE;
		break;

	    case UNSIGNED_INT:
		unsigned int it6;
		fscanf ( fp, "%o ", &it6 );
		buf = new FileSendBuffer ( (void *)&it6, UNSIGNED_INT, 
						sizeof(unsigned int));
		recvList->Append ( (void *)buf );
		checksum -= UNSIGNED_INT;
		break;

	    case SHORT:
		short it7;
		fscanf ( fp, "%hd ", &it7 );
		buf = new FileSendBuffer ( (void *)&it7, SHORT, sizeof(short));
		recvList->Append ( (void *)buf );
		checksum -= SHORT;
		break;

	    case UNSIGNED_SHORT:
		unsigned short it8;
		fscanf ( fp, "%ho ", &it8 );
		buf = new FileSendBuffer ( (void *)&it8, UNSIGNED_SHORT, 
						sizeof(unsigned short));
		recvList->Append ( (void *)buf );
		checksum -= UNSIGNED_SHORT;
		break;

	    case UNSIGNED_LONG:
		unsigned long it9;
		fscanf ( fp, "%lo ", &it9 );
		buf = new FileSendBuffer ( (void *)&it9, UNSIGNED_LONG, 
						sizeof(unsigned long));
		recvList->Append ( (void *)buf );
		checksum -= UNSIGNED_LONG;
		break;

	    case STRING:
		{
		int sz;
		fscanf ( fp, "%d ", &sz );
		char *it10 = new char[sz + 1];
		it10[sz] = '\0';
		for ( int z = 0; z < sz; z++ )
		    fscanf ( fp, "%c", &it10[z] );
		buf = new FileSendBuffer ( (void *)it10, STRING, 
					strlen(it10) + 1);
		recvList->Append ( (void *)buf );
		delete []it10;
		checksum -= STRING;
		break;
		}

	    default:
		MWprintf (10, "Couldn't decipher the data type in receiving\n");
		MWprintf (10, "Got %d\n", typeno );
		Close ( fp );
		return UNKNOWN_DATA_TYPE;
	}
    }
    if ( checksum != 0 )
    {
    	// Problem 
	MWprintf ( 10, "The checksum is not zero in recv\n");
	while ( recvList->IsEmpty() != TRUE )
	{
		FileSendBuffer *buff = (FileSendBuffer *)recvList->Remove();
		delete buff;
	}

	FileSendBuffer *newbuf = new FileSendBuffer ( &FileRCID, INT,
				sizeof(int) );
	recvList->Append ( (void *)newbuf );

	msgTag = CHECKSUM_ERROR;
    }
    
    Close ( fp );
    HardDeleteFile ( filename );
    return 0;
}


//---------------------------------------------------------------------------
// 	MWFileRC::CheckLogFiles():
//	Sophisticated check log function which never fails to detect any event
// in the log.
//---------------------------------------------------------------------------
void 
MWFileRC::CheckLogFilesRunning ( )
{
#ifdef FILE_MASTER
    int i;
    char f_name[_POSIX_PATH_MAX];

    for ( i = 0; i < target_num_workers; i++ ) {
				//Discuss about this thing.
			if ( fileWorkers[i].state == FILE_FREE ) continue;
			int num_event = 0;

			ReadUserLog log;
			ULogEvent *event;
			ULogEventOutcome outcome;

			sprintf( f_name, "%s/log_file.%d", control_directory, i );

			FILE *f;
			if ( (f = fopen( f_name, "r" )) == NULL )
				{
					fileWorkers[i].state = FILE_IDLE;
					fileWorkers[i].event_no = 0;
					continue;
				}
            fclose(f);

			log.initialize ( f_name );

				// Skip over events we've seen before.
			while( (outcome = log.readEvent(event)) == ULOG_OK)	{
				num_event++;
				if ( num_event >= fileWorkers[i].event_no)	{
					break;
				}
				delete event;
			}

				// If we hit EOF in this userlog, no new events
			if (outcome != ULOG_OK) {
				delete event;
				continue; 
			}

				// Make sure this event belongs to this job
			while ( event->cluster != fileWorkers[i].condorID || event->proc != fileWorkers[i].condorprocID )	{
					MWprintf ( 10, "How can an event belonging to a different cluster be printed here\n");
					num_event++;
					fileWorkers[i].event_no++;
					delete event;
					outcome = log.readEvent(event);
					if (outcome != ULOG_OK) {
						break;
					}
			}

			if (outcome != ULOG_OK) {
				delete event;
				continue;
			}

			switch (event->eventNumber) {

				case ULOG_JOB_ABORTED:	
				case ULOG_JOB_TERMINATED:
				case ULOG_EXECUTABLE_ERROR:
					{
						fileWorkers[i].state = FILE_KILLED;
						break;
					}

				case ULOG_SHADOW_EXCEPTION:
					{
						fileWorkers[i].state = FILE_TRANSIT;
						break;
					}
				case ULOG_CHECKPOINTED:
					{
							//fileWorkers[i].state = FILE_SUSPENDED;
						break;
					}
				case ULOG_SUBMIT:
					{
						MWprintf ( 10, "Got a submit event\n"); 
						fileWorkers[i].state = FILE_SUBMITTED;
						break;
					}
				case ULOG_JOB_EVICTED:
					{
						if ( want_checkpointing == 1 )
							{      
								JobEvictedEvent *evictEvent = (JobEvictedEvent *)event;
					
								if ( evictEvent->checkpointed == TRUE ) 
									{
										fileWorkers[i].state = FILE_SUSPENDED;
										MWprintf ( 10, "Got an evict with checkpoint event from %d\n",i );
									}
								else 
									{
										fileWorkers[i].state = FILE_TRANSIT;
									}
							}
						else
							{
								fileWorkers[i].state = FILE_TRANSIT;
							}

						break;
					}	
				case ULOG_EXECUTE:
					{
						if ( fileWorkers[i].state == FILE_SUSPENDED )
							{
								MWprintf (10,"In LogFiles gotResumed after suspend %d\n",i);
								fileWorkers[i].state = FILE_RESUMED;
							}
						else if ( fileWorkers[i].state == FILE_SUBMITTED )
							fileWorkers[i].state = FILE_EXECUTE;
						else {
							MWprintf ( 10, "Unexpected state transition to execute from %d in submitfile %s\n", fileWorkers[i].state, f_name );
							MWprintf ( 10, "Forcing this worker to TRANSIT state\n");
							fileWorkers[i].state = FILE_TRANSIT;
						}
						break;
					}

				default:
					{
						fileWorkers[i].event_no++;
						break;
					}
			}
			fileWorkers[i].event_no++;
			delete event;
	}
#endif
}


//---------------------------------------------------------------------------
// MWFileRC::CheckFilesResuscicate ():
//	Called when recovering from a crash.
//---------------------------------------------------------------------------
void 
MWFileRC::CheckLogFilesResuscicate ( char *f_name, struct FileWorker &fw )
{
#ifdef FILE_MASTER

	ReadUserLog log;
	ULogEvent *event;
	ULogEventOutcome outcome;

	FILE *f;
	if ( (f = fopen ( f_name, "r" )) == NULL )
	{
		fw.state = FILE_FREE;
		fw.event_no = 0;
		return;
	}
	fclose(f);

	log.initialize ( f_name );

	outcome = log.readEvent( event );
	for ( ;  outcome == ULOG_OK ;outcome = log.readEvent( event ) )
	{
		if ( event->eventNumber == ULOG_JOB_ABORTED )
		{
			fw.state = FILE_KILLED;
			fw.event_no++;
			delete event;
			return;
		}
		else if ( event->eventNumber == ULOG_SHADOW_EXCEPTION )
		{
			char exe[128];
			sprintf ( exe, "%s/bin/condor_rm %d.%d", CONDOR_DIR, fw.condorID, fw.condorprocID);
			system ( exe );
			fw.state = FILE_KILLED;
			fw.event_no++;
			delete event;
			return;
		}
		else if ( event->eventNumber == ULOG_JOB_TERMINATED )
		{
			fw.state = FILE_KILLED;
			fw.event_no++;
			delete event;
			return;
		}
		else if ( event->eventNumber == ULOG_EXECUTABLE_ERROR )
		{
			fw.state = FILE_KILLED;
			fw.event_no++;
			delete event;
			return;
		}
		else if ( event->eventNumber == ULOG_CHECKPOINTED )
		{
			fw.event_no++;
			delete event;
		}
		else if ( event->eventNumber == ULOG_SUBMIT )
		{
			fw.state = FILE_SUBMITTED;
			fw.event_no++;
			delete event;
		}
		else if ( event->eventNumber == ULOG_JOB_EVICTED )
		{
			if ( want_checkpointing == 1 )
	   		{
				JobEvictedEvent *evictEvent = (JobEvictedEvent *)event;
				if ( evictEvent->checkpointed == TRUE ) 
				{
					fw.state = FILE_SUSPENDED;
					fw.event_no++;
					delete event;
				}
				else 
				{
					char exe[128];
					sprintf ( exe, "%s/bin/condor_rm %d.%d", CONDOR_DIR, fw.condorID,
										fw.condorprocID);
					system ( exe );
					fw.state = FILE_KILLED;
					delete event;
					break;
				}
			}
			else 	      
			{
				char exe[128];
				sprintf ( exe, "%s/bin/condor_rm %d.%d", CONDOR_DIR, fw.condorID,
									fw.condorprocID);
				system ( exe );
				fw.state = FILE_KILLED;
				fw.event_no++;
				delete event;
			}
		}
		else if ( event->eventNumber == ULOG_EXECUTE )
		{
			fw.state = FILE_RUNNING;
			fw.event_no++;
			delete event;
		}
		else
		{
			fw.event_no++;
			delete event;
		}
	}

	if ( fw.state == FILE_KILLED )
	{
		fw.state = FILE_FREE;
	}

#endif
}


//---------------------------------------------------------------------------
// MWFileRC::resuscicate ( )
// Sanjeev : This is THE function of the whole program. It is this routine that
// gives MWFileRC that robustness that one can only dream right now.
// The real Sanjeev : Robustness? You must be Just kidding!!
//---------------------------------------------------------------------------
void
MWFileRC::resuscicate ()
{
	int  i;
	char c[2]; c[1] = '\0';
	char worker_state_file[_POSIX_PATH_MAX];
	char master_waitfile[_POSIX_PATH_MAX], worker_waitfile[_POSIX_PATH_MAX];
	for ( i = 0; i < target_num_workers; i++ )
	{
		if ( fileWorkers[i].id < 0 ) continue;

		sprintf ( worker_state_file, "%s/log_file.%d", control_directory, fileWorkers[i].id );
		sprintf ( master_waitfile, "./%s/master_waitfile.%d", output_directory, fileWorkers[i].id );
		sprintf ( worker_waitfile, "./%s/worker_waitfile.%d", input_directory, fileWorkers[i].id );

		fileWorkers[i].worker_counter = GetWorkerCounter ( worker_waitfile );
		fileWorkers[i].counter = GetCounter ( master_waitfile );
		MWprintf ( 10, "The counters of %d are %d and %d and %d\n", i, fileWorkers[i].counter, 
										fileWorkers[i].worker_counter, fileWorkers[i].condorID );

		CheckLogFilesResuscicate( worker_state_file, fileWorkers[i] );
	}
}

//---------------------------------------------------------------------------
//Master::GetCondorId():
//      Read the open log file and get the condor cluster of the log file.
//---------------------------------------------------------------------------
void
MWFileRC::GetCondorId ( char *lgfile, int *cID, int *pID )
{
#ifdef FILE_MASTER

    ReadUserLog log;
    ULogEvent *event;
    ULogEventOutcome outcome;

    log.initialize ( lgfile );
    outcome = log.readEvent( event );
    while ( outcome == ULOG_OK && event->eventNumber != ULOG_SUBMIT )
    {
	outcome = log.readEvent ( event );
    }

    if ( outcome != ULOG_OK ) 
    {
	*cID = -1;
	*pID = -1;
	return;
    }

    *cID = event->cluster;
    *pID = event->proc;
    return;
#endif
    return;
}


void
MWFileRC::sort_exec_class_ratio ( int *temp )
{
	int i, j;
	double ii, jj;
	for ( i = 0; i < exec_classes; i++ )
	{
		temp[i] = i;
	}
	for ( i = 0; i < exec_classes; i++ )
	{
		for ( j = i + 1; j < exec_classes; j++ )
		{
			ii = ((double)MW_exec_class_num_workers[temp[i]]) / exec_class_target_num_workers[temp[i]];
			jj = ((double)MW_exec_class_num_workers[temp[j]]) / exec_class_target_num_workers[temp[j]];
			if ( jj > ii )
			{
				int temp = j;
				j = i;
				i = temp;
			}
		}
	}
	return;
}

//---------------------------------------------------------------------------
// hostaddlogic():
//	If the number of workers are less than the target number of workers 
// just spawns the remaining. current_num_workers holds the information about the
// the current number of workers.
//---------------------------------------------------------------------------
int
MWFileRC::hostaddlogic ( int *num_workers )
{
	/* This awfully complex function determines the conditions
	under which it being called, and then calls condor_submit
	asking for the correct number of hosts.

	- num_workers - an array of size exec_classes that contains
	the number of workers in each arch class.
	- hostadd_reqs: The number of outstanding HOSTADDs out there
	- target_num_workers: Set by the user...
	*/

	MWprintf ( 10, "In hostaddlogic \n" );
	if ( !hostadd_reqs ) 
	{
		/* if this doesn't exist yet, we won't do anything */
		return 0; 
	}       

	/* number of hosts to ask for at a time. The idea is that we'll
	have double this outstanding at a time - and this amounts
	to 12 for 1, 2, or 3 arch classes. */
	int i;          
	int *sorted_order = new int [exec_classes];
	int *num_execs = new int[exec_classes];
	int cur;
	int req;
	sort_exec_class_ratio ( sorted_order );
	int status;
                        
	MWprintf ( 60, "hal: target: %d.\n", target_num_workers );

	int HOSTINC = hostinc_;

	for ( i = 0; i < exec_classes; i++ )
	{
		cur = sorted_order[i];
		if ( MW_exec_class_num_workers[cur] + hostadd_reqs[cur] >= exec_class_target_num_workers[cur] ) continue;
		if ( hostadd_reqs[cur] > 0 ) continue;
		if ( exec_class_target_num_workers[cur] - ( MW_exec_class_num_workers[cur] + hostadd_reqs[cur] ) > HOSTINC )
			req = HOSTINC;
		else
		       req = exec_class_target_num_workers[cur] - ( MW_exec_class_num_workers[cur] + hostadd_reqs[cur] );
		status = do_spawn ( req, cur );
		if ( status > 0 )
		{
			MWprintf ( 10, "Added %d workers of exec_class %d\n", req, i );
			hostadd_reqs[cur] += status;
		}
		else
		{
			MWprintf ( 10, "ERROR!! Failed to submit workers of exec_class %d\n", i );
		}
	}

	write_RMstate ( );

	delete [] sorted_order;
	delete [] num_execs;
	return 0;
}


//---------------------------------------------------------------------------
// GetWorkerCounter ():
//	Gets the next message to be sent to the worker.
//---------------------------------------------------------------------------
int 
MWFileRC::GetWorkerCounter ( char *file )
{
    int cal = 0;
    int task_num = 0; 
    char c[2]; c[1] = '\0';
    bool end = FALSE;
    bool first = TRUE;

    if ( !file ) return 0;

    FILE *fp = Open ( file, "r" );
    if ( fp == NULL ) return 0;

    while ( fscanf(fp, "%c", &c[0] ) >= 0 )
    {
	if ( isdigit ( c[0] ) > 0 )
	{
	    cal = cal * 10 + atoi(c);
	    end = FALSE;
	}
	else if ( isspace(c[0]) > 0 )
	{
	}
	else if ( c[0] == OK )
	{
	    task_num = cal;
	    cal = 0;
	    end = TRUE;
	    first = FALSE;
	}
	else if ( c[0] == SYNC || c[0] == ESYNC )
	{
	    cal = 0;
	    end = TRUE;
	}
    }
    Close ( fp );
    if ( first == TRUE ) return 0;
    if ( end == TRUE ) 
    {
	return task_num + 1;
    }
    else 
    {
	return task_num;
    }
}



//---------------------------------------------------------------------------
// GetCounter ():
//	Gets the next message to be received from the worker.
//---------------------------------------------------------------------------
int 
MWFileRC::GetCounter ( char *file )
{

	char c[2]; c[1] = '\0';
	int num_cal = 0;
	int task_num = 0;
	bool end = TRUE;

	FILE *fp = Open ( file, "r" );
	if ( !fp ) return 0;

	while ( fscanf(fp, "%c", &c[0] ) >= 0 )
    {
		if ( isdigit ( c[0] ) > 0 )
		{
			num_cal = num_cal * 10 + atoi(c);
			end = FALSE;
		}
		else if ( c[0] == OK )
		{
			task_num = num_cal;
			num_cal = 0;
			end = TRUE;
		}
	}
	Close ( fp );
	
	if ( end ) return task_num + 1;
	else return task_num;
}


void
MWFileRC::killWorker ( int i )
{
    char exe[_POSIX_PATH_MAX];
    MWprintf ( 10, "Killing condor job %d.%d\n", fileWorkers[i].condorID,
    						fileWorkers[i].condorprocID );
    sprintf ( exe, "%s/bin/condor_rm %d.%d", CONDOR_DIR, fileWorkers[i].condorID,
						fileWorkers[i].condorprocID);
    system ( exe );
}


void
MWFileRC::InitStructures ( )
{
	int i;

	fileWorkers = new struct FileWorker[target_num_workers];

	for ( i = 0; i < target_num_workers; i++ )
	{
		fileWorkers[i].state = FILE_FREE;
		fileWorkers[i].served = -1;
		fileWorkers[i].event_no = 0;
		fileWorkers[i].id = -1;
	}

	hostadd_reqs = new int[exec_classes];
	for ( i = 0; i < exec_classes; i++ )
		hostadd_reqs[i] = 0;
}


// Here you have to resuscicate and fill in the appropriate structures
// By This time target_num_workers has been filled up.
// We should not touch target_num_workers. All things go into current_num_workers.
int
MWFileRC::read_RMstate ( FILE *fp )
{
	int i, temp;
	int temp1, temp2, temp3, temp4;
	char target_file[_POSIX_PATH_MAX];

	sprintf ( target_file, "%s/%s", control_directory, moment_worker_file );

	
	FILE *numwf = Open ( target_file, "r" );
	if ( numwf == NULL )
	{
		MWprintf ( 10, "Coudn't open the moment_worker_file for reading\n");
		InitStructures ( );
		return 0;
	}

	if ( fscanf ( numwf, "%d %d", &temp, &subId ) < 0 )
	{
		MWprintf (10, "Was unable to read the moment worker file\n");
		InitStructures ( );
		Close ( numwf );
		return 0;
	}

	if ( temp != target_num_workers )
	{
		MWprintf ( 10, "Something wrong with the target_num_workers\n");
		if ( temp > target_num_workers )
			target_num_workers = temp;
	}

	InitStructures ( );

	while ( fscanf ( numwf, "%d %d %d %d ", &temp1, &temp2, &temp3, &temp4 ) > 0 )
	{
		if ( temp1 > target_num_workers )
		{
			MWprintf ( 10, "ERROR!! Some Big Time Screwup, how can id be greater than target_num_workers\n");
			::exit(1);
		}
		fileWorkers[temp1].id = temp1;
		fileWorkers[temp1].condorID = temp2;
		fileWorkers[temp1].condorprocID = temp3;
		fileWorkers[temp1].exec_class = temp4;
	}

    Close ( numwf );

	resuscicate();

	for ( i = 0; i < target_num_workers; i++ )
	{
		if ( fileWorkers[i].state == FILE_SUBMITTED )
			hostadd_reqs[fileWorkers[i].exec_class]++;
	}

    return 0;
}


int
MWFileRC::write_RMstate ( FILE *filp )
{
	int i;
	FILE *fp;
	char target_file[_POSIX_PATH_MAX];
	char cmd[4 * _POSIX_PATH_MAX];
    
	sprintf ( target_file, "%s/%s", control_directory, moment_worker_file );
	
	fp = Open ( "temp_File", "w" );
	if  ( fp == NULL )
	{
		MWprintf ( 10, "Cannot open the file temp_File for writing\n");
		return -1;
	}

//	fprintf ( fp, "%d %d ", max_num_workers, subId );
	fprintf ( fp, "%d %d ", target_num_workers, subId );

	for ( i = 0; i < exec_classes; i++ )
	{
//		fprintf ( fp, "%d ", hostadd_reqs[i] );
		for ( int j = 0; j < target_num_workers; j++ )
		{
			if ( fileWorkers[j].state != FILE_FREE && fileWorkers[j].exec_class == i )
				fprintf ( fp, "%d %d %d %d ", j, fileWorkers[j].condorID, fileWorkers[j].condorprocID, fileWorkers[j].exec_class );
		}
/*
		for ( int j = 0; j < 2 * HOSTINC; j++ )
			if ( hostaddind_reqs[i][j] != -1 ) 
				fprintf ( fp, "%d ", hostaddind_reqs[i][j] );
*/
	}

	Close ( fp );

	sprintf ( cmd, "mv temp_File %s", target_file );
//	if ( rename ( "temp_File", target_file ) < 0 )
	if ( system ( cmd ) < 0 )
	{
		MWprintf (10, "Could not move the new RM_status file %d\n", errno);
		return -1;
	}

	return 0;
}

char*
MWFileRC::process_executable_name ( char *exec, int ex_cl, int ar_cl )
{
	int i;
	char newstring[_POSIX_PATH_MAX];
	char arch[100], opsys[100];
	char *index;
	char *req = arch_class_attributes[ar_cl];
	char *newexec = new char[_POSIX_PATH_MAX];
	char exe[4 * _POSIX_PATH_MAX];
	strcpy ( newstring, req );
	sprintf ( newexec, "mw_exec%d", ex_cl );

	for ( unsigned i = 0; i < strlen ( newstring ); i++ )
	{
		int k = newstring[i];
		if( ! ( ( k == 'u' || k == 'x' ) && newstring[i-1] == '4' ) ) 
		{
			k = toupper(k);
		}
		newstring[i] = k;
	}

    i = 0;
    index = strstr ( newstring, "OPSYS" );
    index += 5;
    while ( *index == ' ' || *index == '=' || *index == '"' || *index == '\\' ) index++;
    while ( *index != ')' && *index != '"' && *index != '\\' ) opsys[i++] = *index++;
    opsys[i] = '\0';

    i = 0;
    index = strstr ( newstring, "ARCH" );
    index += 4;
    while ( *index == ' ' || *index == '=' || *index == '"' || *index == '\\' ) index++;
    while ( *index != ')' && *index != '"' && *index != '\\' ) arch[i++] = *index++;
    arch[i] = '\0';

    sprintf ( newexec, "%s.%s.%s", newexec, opsys, arch );
    remove ( newexec );

    MWprintf ( 10, "Making a link from %s to %s\n", exec, newexec );
#ifdef WINDOWS
    sprintf ( exe, "copy %s %s", exec, newexec );
#else
    sprintf ( exe, "ln -s %s %s", exec, newexec );
#endif
    system ( exe );

    return newexec;
}

//BUFFER
int 
MWFileRC::recv_all( int from_whom, int msgtag )
{
	MWprintf(91, "Not implemented yet for MW-File.\n");
	return 0;
}

int 
MWFileRC::setrbuf( int bid )
{
	MWprintf(91, "Not implemented yet for MW-File.\n");
	return -1;
}

int 
MWFileRC::freebuf( int bid )
{
	MWprintf(91, "Not implemented yet for MW-File.\n");
	return -1;
}

MWList<void> * 
MWFileRC::recv_buffers()
{
	MWprintf(91, "Not implemented yet for MW-File.\n");
	return NULL;
}

int 
MWFileRC::next_buf()
{
	MWprintf(91, "Not implemented yet for MW-File.\n");
	return -1;
}
	


int
MWFileRC::sendFileToMaster(char *filename) {
#ifdef USE_CHIRP
#ifndef FILE_MASTER
	struct chirp_client *client;

		// We connect each time, so that we don't have thousands
		// of idle connections hanging around the master
	MWprintf(60, "Sending %s to master\n", filename);
	client = chirp_client_connect_default();
	if (!client) {
		MWprintf(10, "MWFileRC:sendFileToMaster cannot chirp_connect to master, exitting\n");
		exit(-1);
	}

	FILE *rfd = ::fopen(filename, "r");
	if (rfd == NULL) {
		chirp_client_disconnect(client);
		MWprintf(10, "MWFileRC:sendFileToMaster can't open filename %s:%s\n", filename, strerror(errno));
		return -1;
	}

	int wfd = chirp_client_open(client, filename, "cwat", 0777);
	if (wfd < 0) {
		::fclose(rfd);
		chirp_client_disconnect(client);
		MWprintf(10, "MWFfileRC:sendFileToMaster can't chirp_client_open %s:%d\n", filename, wfd);
		return -1;
	}

	char buf[1024];

	int num_read = 0;
	do {
		num_read = ::fread(buf, 1, 1024, rfd);
		if (num_read < 0) {
			fclose(rfd);
			chirp_client_close(client, wfd);
			chirp_client_disconnect(client);
			MWprintf(10, "MWFileRC:sendFileToMaster local read error on %s\n", filename);
			return -1;
		}
		int num_written = chirp_client_write(client, wfd, buf, num_read);
		if (num_written != num_read) {
			fclose(rfd);
			chirp_client_close(client, wfd);
			chirp_client_disconnect(client);
			MWprintf(10, "MWFileRC:sendFileToMaster couldn't chirp_write as much as we read\n");
			return -1;
		}

	} while (num_read == 1024);

	::fclose(rfd);
	chirp_client_close(client, wfd);
	chirp_client_disconnect(client);
		
#endif
#endif
	return 0;
}

int
MWFileRC::getFileFromMaster(char *filename) {
#ifdef USE_CHIRP
#ifndef FILE_MASTER
	struct chirp_client *client = 0;

	MWprintf(60, "Fetching %s from master\n", filename);

		// We connect each time, so that we don't have thousands
		// of idle connections hanging around the master

	client = chirp_client_connect_default();
	if (!client) {
		MWprintf(10, "MWFileRC:getFileFromMaster cannot chirp_connect to master, waiting\n");
		return -1;
	}

	int rfd = chirp_client_open(client, filename, "r", 0);
	if (rfd < 0) {
		MWprintf(10, "MWFileRC:getFileFromMaster can't chirp_client_open %s:%d\n", filename, rfd);
		chirp_client_disconnect(client);
		return -1;
	}

	FILE *wfd = ::fopen(filename, "w+");
	if (wfd == NULL) {
		MWprintf(10, "MWFileRC:getFileFromMaster can't open filename %s: %s\n", filename, strerror(errno));
		chirp_client_close(client, rfd);
		chirp_client_disconnect(client);
		return -1;
	}


	char buf[1024];

	int num_read = 0;
	do {
		num_read = chirp_client_read(client, rfd, buf, 1024);
		if (num_read < 0) {
			MWprintf(10, "MWFileRC:getFileFromMaster couldn't chirp_read\n");
			::fclose(wfd);
			chirp_client_close(client, rfd);
			chirp_client_disconnect(client);
			return -1;
		}

		int num_written = fwrite(buf, 1, num_read, wfd);
		if (num_written < 0) {
			MWprintf(10, "MWFileRC:GetFileFromMaster local read error on %s\n", filename);
			::fclose(wfd);
			chirp_client_close(client, rfd);
			chirp_client_disconnect(client);
			return -1;
		}

	} while (num_read == 1024);

	::fclose(wfd);
	chirp_client_close(client, rfd);
	chirp_client_disconnect(client);
	return 0;
#endif
#endif
	return 0;

}

