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
#include "MWCondorPvmRC.h"
#include <MW.h>
#include <MWDriver.h>
#include <MWTask.h>
#include <MWWorker.h>

#ifdef PVM_MASTER
MWRMComm * MWDriver::RMC = new MWPvmRC;
MWRMComm * MWTask::RMC = MWDriver::RMC;
MWRMComm * MWWorker::RMC = NULL;
#else
MWRMComm * MWWorker::RMC = new MWPvmRC;
MWRMComm * MWTask::RMC = MWWorker::RMC;
MWRMComm * MWDriver::RMC = NULL;
#endif

extern int *MW_exec_class_num_workers;

MWPvmRC::MWPvmRC ( )
{ 
	hostadd_reqs = NULL;
	is_master = FALSE;   // assume the worst...
	submit_list = new MWList<void>;
}

MWPvmRC::~MWPvmRC ( )
{
	if ( hostadd_reqs )
	{
		for ( int i = 0; i < exec_classes; i++ )
			delete [] hostadd_reqs[i];
		delete [] hostadd_reqs;
	}
	if ( submit_list )
	{
		while ( submit_list->number() > 0 )
			delete ((int *)submit_list->Remove());
		delete submit_list;
	}
	
	if (recv_buf_list) {
		while ( recv_buf_list->number() > 0 )
			delete ((int *)recv_buf_list->Remove());
		delete recv_buf_list;
	}
}

/* Implementation of comm / rm class for pvm. */
void 
MWPvmRC::exit( int exitval ) 
{
	MWprintf ( 50, "Before pvm_exit\n" );
	pvm_catchout( 0 );
	pvm_exit();
	MWprintf ( 50, "After pvm_exit\n" );
	::exit(exitval);
}

int  
MWPvmRC::setup( int argc, char *argv[], int *id, int *master_id ) 
{
	
	MWprintf ( 10, "In MWPvmRC::setup()\n");

	*id = pvm_mytid();
	MWprintf(10, "MyPid is %d(%x)\n", *id, *id );

	// BUFFER - init the recv_buf_list
	recv_buf_list = new MWList<void>;
	
	if ( ( *master_id = pvm_parent() ) == PvmNoParent ) {
	   	MWprintf ( 40, "I have no parent, therefore I'm the master.\n" );
	   	is_master = TRUE;
	   	*master_id = 0;
   	} else {
	   	is_master = FALSE;
		return 0;
	}
    
	pvm_catchout( stderr );
		/* We tell PVM to tell us when hosts are added.  We will be sent
		   a message with tag HOSTADD.  The -1 turns the notification on.
		*/
	pvm_notify ( PvmHostAdd, HOSTADD, -1, NULL );

	return 0;
}

int 
MWPvmRC::config( int *nhosts, int *narches, MWWorkerID ***workers ) 
{
	if ( !is_master ) {
		MWprintf ( 10, "Slaves not allowed in this privaledged call!\n" );
		return -1;
	}

	struct pvmhostinfo *hi= NULL;
	int r;
	
	MWprintf ( 70, "In MWPvmRC::config()\n" );

	if ( *workers ) {
		MWprintf ( 10, "workers should be NULL when called.\n" );
		return -1;
	}

	if ( (r = pvm_config ( nhosts, narches, &hi ) ) < 0 ) {
		MWprintf ( 10, "Return value of %d from pvm_config!\n", r );
		return -1;
	}
	
	(*workers) = new MWWorkerID*[*nhosts];
	for ( int i=0 ; i<*nhosts ; i++ ) {
		(*workers)[i] = new MWWorkerID;
		(*workers)[i]->set_arch ( atoi ( hi[i].hi_arch ) );
			/* Do domething creative with speed later! */
	}

	return 0;
} 

int
MWPvmRC::start_worker ( MWWorkerID *w ) 
{
	struct condorpvm_submit_element *elem = NULL;
	// int ex_cl;
	int retval;

	if ( !is_master ) 
	{
		MWprintf ( 10, "Slaves not allowed in this privaledged call!\n" );
		return -1;
	}

	MWprintf ( 50, "In MWPvmRC::start_worker()\n" );
	if ( !w ) 
	{
		MWprintf ( 10, "w cannot be NULL!\n" );
		return -1;
	}

		/* First, we will unpack some goodies that pvm gave us: */
    int val;
    char name[64], arch[64];
    int status, dsig;

		/* The following is packed by the multishadow along with the 
		   HOSTADD message.  We unpack it here to put into the MWWorkerID */
    unpack ( &val, 1, 1 );
    unpack ( name );
    unpack ( arch );
    unpack ( &status, 1, 1 );
    unpack ( &dsig, 1, 1 );  // new in 3.4.1
    val &= PVM_MASK;   // remove that annoying high bit
    MWprintf ( 40, "Received a host: " );
    MWprintf ( 40, "Host %x added (%s), arch: %s, status:%d\n", 
			   val, name, arch, status );

	w->set_arch( atoi(arch) );
//	w->set_machine_name( name );
	elem = find_exec_class ( atoi(arch) );
	if ( !elem )
	{
		MWprintf ( 10, "Strange!! No pending request found for arch %d\n", atoi(arch) );
		w->set_id1(-1);
		return -1;
	}
	else
	{
		w->set_exec_class ( elem->exec_class );
		w->set_executable ( elem->exec );
		retval = do_spawn( w );
	}

		/* We can do this now because the request will go down regardless
		   of whether or not the spawn succeeds... */
	hostadd_reqs[elem->exec_class][w->get_arch()]--;
	delete elem;
	return retval;
}

int
MWPvmRC::init_beginning_workers( int *nworkers, MWWorkerID ***workers ) 
{
	if ( !is_master ) 
	{
		MWprintf ( 10, "Slaves not allowed in this privaledged call!\n" );
		return -1;
	}

	int i, j, narches;
	MWprintf ( 50, "In MWPvmRC::init_beginning_workers\n" );

		/* config allocates memory for workers */
	if ( config( nworkers, &narches, workers ) < 0 ) 
	{
		return -1;
	}
	
		/* num_arches is a member of MWCommRC */
	if ( (num_arches != -1) && (num_arches != narches) ) 
	{
			/* We were told by someone that we had the wrong number of 
			   arch classes.  Complain here. */
		MWprintf ( 10, "num_arches was set to %d, pvm just told us %d!\n", 
				   num_arches, narches );
		return -1;
	}
	num_arches = narches;

    MWprintf ( 40, "Number of hosts at start: %d.\n", *nworkers );
    MWprintf ( 40, "Number of arch. classes:  %d.\n", num_arches );

	hostadd_reqs = new int*[exec_classes];
	for ( i = 0; i < exec_classes; i++ )
	{
		hostadd_reqs[i] = new int[num_arches];
		for ( j = 0 ; j < num_arches ; j++ ) 
		{
			hostadd_reqs[i][j] = 0;
		}
	}

	for ( i=0 ; i<(*nworkers) ; i++ ) 
	{
		do_spawn( (*workers)[i] );
	}
	
	return 0;
}

int
MWPvmRC::restart_beginning_workers ( int *nworkers, MWWorkerID ***workers, MWmessages msg )
{
	return init_beginning_workers ( nworkers, workers );
}

int
MWPvmRC::do_spawn ( MWWorkerID *w ) 
{
	int spawned_tid;
	int exec_class;
	char arch[4];
	char *exec;
	int archnum = w->get_arch();
	sprintf ( arch, "%d", archnum );
	int i;

	if ( w->get_executable() == NULL )
	{
		int tempi[exec_classes];

		for ( i = 0; i < exec_classes; i++ )
		{
			tempi[i] = 0;
		}

		for ( i = 0; i < num_executables; i++ )
		{
			if ( worker_executables[i]->arch_class == w->get_arch() )
			tempi[worker_executables[i]->exec_class] = 1;
		}
	
		exec_class = choose_exec_class ( tempi );
		if ( exec_class < 0 )
		{
			MWprintf ( 10, "ERROR!! No suitable exec_class chosen for arch_class %d\n", w->get_arch() );
			return -1;
		}
		exec = find_executable ( exec_class, w->get_arch() );
		w->set_executable ( exec );
		w->set_exec_class ( exec_class );
	}
	else
	{
		exec_class = w->get_exec_class();
		exec = w->get_executable();
	}

	MWprintf ( 40, "About to call pvm_spawn (%s)...\n", exec );
	
#if ! defined( STANDALONE )
    int status = pvm_spawn( exec, NULL, 2, arch, 1, &spawned_tid );
#else    
    int status = pvm_spawn( exec, NULL, 1, "wheel", 1, &spawned_tid );
#endif
	
	if ( status == 1 ) 
	{
        MWprintf ( 40, "Success spawning as tid %x.\n", spawned_tid );
		w->set_id1( ( spawned_tid & PVM_MASK ) );
		w->set_id2( ( pvm_tidtohost( spawned_tid ) & PVM_MASK ) );
		setup_notifies ( spawned_tid );
		return spawned_tid;
	}
	else
	{
        MWprintf ( 40, "Huh? Error: status:%d, tid:%d\n", 
				   status, spawned_tid );
		MWprintf ( 40, "There was a problem spawning on %x.\n", w->get_id1() );

			/* We make sure that id2 is -1, signifying error! */
		// w->set_id1( -1 );	// FIXME ??? Isn't it set_id2(-1) ?
		w->set_id2(-1);
		return -1;
    }
}

int
MWPvmRC::removeWorker ( MWWorkerID *w ) 
{
	char *machine = w->machine_name();
	int info;
	MWprintf ( 70, "Removing worker %s.\n", machine );
	int retval = pvm_delhosts( &machine, 1, &info );
	MWprintf(31, "pvm_delhosts returned %d and infos[0] = %d\n", retval, info);
	if ( retval < 0 ) {
		MWprintf ( 10, "pvm_delhosts returns an error of %d!\n", retval );
	} 
	if ( info < 0 ) {
		MWprintf ( 10, "pvm_delhosts of %s has error %d.\n", machine, info );
	}
	delete w;
	return retval;
}

int
MWPvmRC::hostaddlogic( int *num_workers ) 
{
		/* This awfully complex function determines the conditions 
		   under which it being called, and then calls pvm_addhosts()
		   asking for the correct number of hosts.
		   
		   - num_workers - an array of size num_arches that contains 
		      the number of workers in each arch class.
		   - hostadd_reqs: The number of outstanding HOSTADDs out there
		   - target_num_workers: Set by the user...
		*/

		/* If we're in STANDALONE mode, we never ask for hosts,
		   so we basically ignore this whole function: */
#if !defined (STANDALONE)

	if ( !hostadd_reqs ) {
			/* if this doesn't exist yet, we won't do anything */
		return 0;
	}

		/* number of hosts to ask for at a time. The idea is that we'll
		   have double this outstanding at a time - and this amounts 
		   to 12 for 1, 2, or 3 arch classes. */

	// Jeff changed this -- we always ask for six extra for each arch class...
	int HOSTINC = hostinc_;
	

	int i, j;
	int *sorted_order = new int[exec_classes];
	sort_exec_class_ratio ( sorted_order );
	int *howmany = new int[num_arches];  // how many to ask for
	int *current_reqs_arches = new int[num_arches];
	int *current_reqs_execs = new int[num_arches];
	int *num_execs = new int[exec_classes];
	// int total_workers = 0;
	// int total_reqs = 0;
	// int total_howmany = 0;
	
	
	// DEBUGGING
	struct pvmhostinfo *hostp;
	int info, k, nhost, narch;

	MWprintf ( 60, "hal: target: %d.\n", target_num_workers );

	for ( i = 0; i < num_arches; i++ )
	{
		current_reqs_arches[i] = 0;
		howmany[i] = 0;
		for ( j = 0; j < exec_classes; j++ )
		{
			current_reqs_arches[i] += hostadd_reqs[j][i];
		}
	}

	for ( i = 0; i < exec_classes; i++ )
	{
		current_reqs_execs[i] = 0;
		num_execs[i] = 0;
		for ( j = 0; j < num_arches; j++ )
		{
			current_reqs_execs[i] += hostadd_reqs[i][j];
		}
	}

	for ( i = 0; i < num_executables; i++ )
	{
		num_execs[ worker_executables[i]->exec_class ]++;
	}

	for ( i = 0 ; i < exec_classes ; i++ ) 
	{
		MWprintf ( 60, "hal: %d workers in class %d\n", num_workers[i], i );
		MWprintf ( 60, "   -- outstanding reqs: %d\n", current_reqs_arches[i] );
		int cur = sorted_order[i];
		int req;

		// DEBUGGING

	        info = pvm_config( &nhost, &narch, &hostp );
		if (info!=0)
			MWprintf(31, "DEBUG: pvm_config() failed, returns %d.\n", info);
		else {
			MWprintf(31, "DEBUG: pvm_config() returned, nhost = %d.\n", nhost);
			if ( (exec_classes == 1) && (nhost > MW_exec_class_num_workers[cur]) ) {
				MWprintf(81, "DEBUG: CondorPvm and MW thinks differently about the number of workers.\n");
				MWprintf(81, "DEBUG: nhost(MW) = %d, nhost(CondorPvm) = %d. Below is the worker list: \n", 
					        MW_exec_class_num_workers[cur], nhost);
				for (k = 0; k < nhost; k++)
					MWprintf(91, "HOST_NAME = %s\n", hostp[k].hi_name);
			}
		}
		
		
		if ( MW_exec_class_num_workers[cur] + current_reqs_execs[cur] >= exec_class_target_num_workers[cur] ) continue;
		if ( current_reqs_execs[cur] > num_execs[cur] * HOSTINC ) continue;
		if ( exec_class_target_num_workers[cur] - ( MW_exec_class_num_workers[cur] + current_reqs_execs[cur] ) > 2 * HOSTINC )
			req = 2 * HOSTINC;
		else
			req = exec_class_target_num_workers[cur] - ( MW_exec_class_num_workers[cur] + current_reqs_execs[cur] );

		for ( j = 0; j < num_arches; j++ )
		{
			int given;
			char *temp;
			if ( !(temp = exists_executable ( cur, j )) ) continue;
			if ( current_reqs_arches[j] > HOSTINC ) continue;
			if ( current_reqs_arches[j] + req > 2 * HOSTINC )
			{
				given = 2 * HOSTINC - current_reqs_arches[j];
			}
			else
			{
				given = req;
			}

			req -= given;
			current_reqs_execs[cur] += given;
			hostadd_reqs[cur][j] += given;
			current_reqs_arches[j] += given;
			howmany[j] += given;

			for ( int k = 0; k < given; k++ )
			{
				struct condorpvm_submit_element *news = new struct condorpvm_submit_element;
				news->arch_class = j;
				news->exec_class = cur;
				news->exec = temp;
				submit_list->Append ( (void *)news );
			}
		}
	}

	/* We're FINALLY ready to do the deed! */
	for ( i = 0 ; i < num_arches ; i++ ) 
	{
		if ( howmany[i] > 0 ) 
			ask_for_host( howmany[i], i );
	}

	delete [] howmany;
	delete [] current_reqs_execs;
	delete [] current_reqs_arches;
	delete [] sorted_order;

#endif /* ifndef STANDALONE */
	return 0;
}

int
MWPvmRC::setup_notifies ( int task_tid ) 
{
		/* Takes a *task* tid, and sets up all notifies at once. */

	MWprintf ( 40, "Setting up notifies for tid %x.\n", task_tid );
	int pvmd_tid = pvm_tidtohost( task_tid );
	
    pvm_notify ( PvmHostDelete, HOSTDELETE, 1, &pvmd_tid );
	pvm_notify ( PvmTaskExit, TASKEXIT, 1, &task_tid );
#if !defined( STANDALONE )
    pvm_notify ( PvmHostSuspend, HOSTSUSPEND, 1, &pvmd_tid );
	pvm_notify ( PvmHostResume,  HOSTRESUME, 1, &pvmd_tid );
#endif

	return 0;
}

int 
MWPvmRC::ask_for_host( int howmany, int arch ) 
{
	int status;
	MWprintf ( 10, "Asking for %d host(s) with arch %d\n", howmany, arch );
	
	char **arches = new char*[howmany];
	for ( int i=0 ; i<howmany ; i++ ) {
		arches[i] = new char[4];
		sprintf( arches[i], "%d", arch );
	}
	int *infos = new int[howmany];
	
	status = pvm_addhosts ( arches, howmany, infos );
	
	for ( int j=0 ; j<howmany ; j++ ) {
		delete [] arches[j];
	}
	delete [] arches;
	delete [] infos;

	if ( status < 0 ) {
		MWprintf ( 10, "pvm_addhosts returns %d, may not be able to "
				   "aquire requested resources...\n", status );
	} else {
		MWprintf ( 90, "pvm_addhosts called with status = %d\n", status );
	}
    
	return status; // JGS
}	

struct condorpvm_submit_element*
MWPvmRC::find_exec_class ( int arch )
{
	struct condorpvm_submit_element *elem;
       
	elem = (struct condorpvm_submit_element *)submit_list->First();
	while ( submit_list->AfterEnd() == false )
	{
		elem = (struct condorpvm_submit_element *)submit_list->Current();
		if ( elem->arch_class == arch ) {
			submit_list->RemoveCurrent();
			return elem;
		}
		else submit_list->Next ();
	}
	return NULL;
}

int
MWPvmRC::choose_exec_class ( int *tempi )
{
	for ( int i = 0; i < exec_classes; i++ )
	{
		if ( tempi[i] == 1 )
			return i;
	}
	return -1;
}

char*
MWPvmRC::find_executable ( int ex_cl, int ar_cl )
{
	for ( int i = 0; i < num_executables; i++ )
	{
		if ( worker_executables[i]->arch_class == ar_cl && 
				worker_executables[i]->exec_class == ex_cl )
			return worker_executables[i]->executable_name;
	}
	return NULL;
}

void
MWPvmRC::sort_exec_class_ratio ( int *temp )
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

char*
MWPvmRC::exists_executable ( int ex_cl, int ar_cl )
{
	for ( int i = 0; i < num_executables; i++ )
	{
		if ( worker_executables[i]->arch_class == ar_cl &&
				worker_executables[i]->exec_class == ex_cl )
			return worker_executables[i]->executable_name;
	}
	return NULL;
}

void 
MWPvmRC::conf() {

    int num_curr_hosts, num_curr_arch;

    struct pvmhostinfo *hi;

        // Foolish pvm_config allocs its own memory...
    pvm_config ( &num_curr_hosts, &num_curr_arch, &hi );

    MWprintf ( 40, "* Pvm says that there are %d machines and %d arches.\n", 
             num_curr_hosts, num_curr_arch );
    
    MWprintf ( 40, "* These machines have the following configurations:\n" );
    MWprintf ( 40, "* tid\t\tname\t\t\tarch\tspeed\n" );
    MWprintf ( 40, "* ---\t\t----\t\t\t----\t-----\n" );

    for ( int i=0 ; i<num_curr_hosts ; i++ ) {
        MWprintf ( 40, "* %x\t%s\t%s\t%d\n", hi[i].hi_tid, hi[i].hi_name, 
                 hi[i].hi_arch, hi[i].hi_speed );
    }

    MWprintf ( 40, "\n" );
}

int
MWPvmRC::initsend ( int encoding ) 
{
	return pvm_initsend ( encoding ); 
}


int
MWPvmRC::send ( int to_whom, int msgtag ) 
{
	return pvm_send ( to_whom, msgtag );
}

int
MWPvmRC::nrecv( int from_whom, int msgtag)
{
    if ( (!is_master) )
        return pvm_nrecv( from_whom, msgtag );

    int *tmp;
    int *buf_id = NULL;

    if ( !(recv_buf_list->IsEmpty()) ) {
        /* Assume that the first message in the list is already processed */
        tmp = (int*) recv_buf_list->Remove();
        MWprintf(91, "Deleted buffer %d\n", *tmp);
        freebuf(*tmp);
    }

    if (recv_buf_list->IsEmpty()) {
        buf_id = new int;
        (*buf_id) = pvm_recv ( from_whom, msgtag );
        if ( buf_id > 0 )  /* attach the buf_id  */
            recv_buf_list->Append((void*)buf_id);
    } else {
        buf_id = (int*)recv_buf_list->First();
    }

    return *buf_id;
}

/* old blocking recv() */
int
MWPvmRC::recv ( int from_whom, int msgtag ) 
{
	recv_all(-1, -1);	
	int *buf_id = NULL;

	if (recv_buf_list == NULL || recv_buf_list->number()==0)
	{
	int *tmp;
	
	if ( !(recv_buf_list->IsEmpty()) ) {
		/* Assume that the first message in the list is already processed */
		tmp = (int*) recv_buf_list->Remove();
		MWprintf(91, "Deleted buffer %d\n", *tmp);
        freebuf(*tmp);
	}
	
	if (recv_buf_list->IsEmpty()) {
		buf_id = new int;
		(*buf_id) = pvm_recv ( from_whom, msgtag );
		if ( buf_id > 0 )  /* attach the buf_id  */
			recv_buf_list->Append((void*)buf_id);
	} else {
		buf_id = (int*)recv_buf_list->First();
	}
	}
	else
	{
		int* bid =  (int*) recv_buf_list->Remove();
		buf_id = bid;
		recv_buf_list->Prepend(bid);
	}

	return *buf_id;
}

/* BUFFER: returns the new buf_id if success, -1 otherwise. */
int 
MWPvmRC::next_buf()
{
	int setr_ret;
	int *bid1, *bid2;

	MWprintf(91, "BUFFER: Entered MWPvmRC::next_buf()\n");

	if ( recv_buf_list->number()<2 )
		return -1;
	
	bid1 = (int*) recv_buf_list->Remove();
	bid2 = (int*) recv_buf_list->Remove();

	if (bid2) {
		setr_ret = pvm_setrbuf(*bid2);
		while ( setr_ret == PvmBadParam || setr_ret == PvmNoSuchBuf ) {
			MWprintf(91, "BUFFER: pvm_setrbuf() failed, can't move to the next buffer. Pvm returns %d.\n", setr_ret);
				
			freebuf(*bid1);
			delete bid1;
			bid1 = bid2;
			bid2 = (int*) recv_buf_list->Remove();
			
			if (bid2)
				setr_ret = pvm_setrbuf(*bid2);
			else {
				MWprintf(91, "BUFFER: In MWPvmRC::next_buf(), the next buffer points to NULL.\n");
				return -1;
			}
		}

		// Prepend the current active buffer
		recv_buf_list->Prepend(bid2);
	}

	return (*bid2);
}

/* BUFFER: returns 1 if switched, returns 0 if no need to switch */
int 
MWPvmRC::switch_buf_back()
{
	int *buf_id;
	
	// if the newly received buffer is the only buffer in recv_buf_list, then return.
	if (recv_buf_list->number()<1)
		return 0;
	
	// otherwise, need to set the active buffer back to the head buffer of recv_buf_list
	buf_id = (int*) recv_buf_list->Remove();
	
	int setr_ret = pvm_setrbuf(*buf_id);
	if (setr_ret == PvmBadParam || setr_ret == PvmNoSuchBuf)
		MWprintf(61, "BUFFER: Can't set active receive buffer to %d.\n", *buf_id);
	else MWprintf(91, "BUFFER: Set active receive buffer back from %d to %d.\n", setr_ret, *buf_id);
	recv_buf_list->Prepend(buf_id);

	return 1;
}

/* BUFFER: receive from pvm and pack into the buffer list */
int	
MWPvmRC::recv_all( int from_whom, int msgtag )
{
	int *buf_id = new int[1];
	int *tmp, len = 0;
	
	if ( !(recv_buf_list->IsEmpty()) ) {
		/* Assume that the first message in the list is already processed */
		tmp = (int*) recv_buf_list->Remove();
		MWprintf(91, "Deleted buffer %d\n", *tmp);
        freebuf(*tmp);
	}
	
	(*buf_id) = pvm_nrecv( from_whom, msgtag );
	
	if ( *buf_id <= 0 ) {
		MWprintf(91, "BUFFER: recv_all() got nothing, Pvm returns %d.\n", *buf_id);
		if (recv_buf_list->number()) 
			switch_buf_back();
		return *buf_id;
	}

	while ( *buf_id > 0 ) {
		MWprintf(91, "BUFFER: got a new buffer %d.\n", *buf_id);
		recv_buf_list->Append(buf_id);
		buf_id = new int [1];
		len ++;
		
		pvm_setrbuf(0);
		(*buf_id) = pvm_nrecv( from_whom, msgtag );
	}
	delete buf_id;

	MWprintf(31, "BUFFER: recv_all() got %d new messages, recv_buf_list length is %d.\n", len, recv_buf_list->number());

	// Return the number of newly received buffers
	return (switch_buf_back()) ? len : 1;
}

/* BUFFER: pvm_setrbuf wrapper */
int
MWPvmRC::setrbuf( int bid )
{
	return pvm_setrbuf( bid );
}

// BUFFER: pvm_freebuf wrapper */
int 
MWPvmRC::freebuf( int bid )
{
	int ret = pvm_freebuf( bid );
	
	if ( ret == PvmNoSuchBuf || ret == PvmBadParam )
		MWprintf(61, "BUFFER: pvm_freebuf(%d) failed, returns %d.\n", bid, ret);
	
	return ret;
}

/* BUFFER: pvm_bufinfo wrapper */
int
MWPvmRC::bufinfo ( int buf_id, int *len, int *tag, int *from ) 
{
	if ( buf_id == 0 ) // might want to check if 0 is a valid id
    {
        *tag = NO_MESSAGE;
        return 0;
    }

	int retval = pvm_bufinfo ( buf_id, len, tag, from );
	/* to clear out the annoying high bit! */
	*from &= PVM_MASK;
	return retval;
}

// BUFFER
MWList<void>*
MWPvmRC::recv_buffers() 
{
	return recv_buf_list;
}

void
MWPvmRC::who ( int *wh )
{
	unpack ( wh, 1, 1 );
	*wh &= PVM_MASK;
}

int
MWPvmRC::pack ( const char *bytes,         int nitem, int stride ) 
{
	bytes_packed_ += nitem / stride;
	return pvm_pkbyte ( const_cast<char *>(bytes), nitem, stride );
}

int
MWPvmRC::pack ( const float *f,            int nitem, int stride ) 
{
	bytes_packed_ += nitem / stride * sizeof(float);
	return pvm_pkfloat ( const_cast<float *>(f), nitem, stride );
}

int
MWPvmRC::pack ( const double *d,           int nitem, int stride ) 
{
	bytes_packed_ += nitem / stride * sizeof(double);
	return pvm_pkdouble ( const_cast<double *>(d), nitem, stride );
}

int
MWPvmRC::pack ( const int *i,              int nitem, int stride ) 
{
	bytes_packed_ += nitem / stride * sizeof(int);
	return pvm_pkint ( const_cast<int *>(i), nitem, stride );
}

int
MWPvmRC::pack ( const unsigned int *ui,    int nitem, int stride ) 
{
	bytes_packed_ += nitem / stride * sizeof(unsigned int);
	return pvm_pkuint ( const_cast<unsigned int *>(ui), nitem, stride );
}

int
MWPvmRC::pack ( const short *sh,           int nitem, int stride ) 
{
	bytes_packed_ += nitem / stride * sizeof(short);
	return pvm_pkshort ( const_cast<short *>(sh), nitem, stride );
}

int
MWPvmRC::pack ( const unsigned short *ush, int nitem, int stride ) 
{
	bytes_packed_ += nitem / stride * sizeof(unsigned short);
	return pvm_pkushort ( const_cast<unsigned short *>(ush), nitem, stride );
}

int
MWPvmRC::pack ( const long *l,             int nitem, int stride ) 
{
	bytes_packed_ += nitem / stride * sizeof(long);
	return pvm_pklong ( const_cast<long *>(l), nitem, stride );
}

int
MWPvmRC::pack ( const unsigned long *ul,   int nitem, int stride ) 
{
	bytes_packed_ += nitem / stride * sizeof(unsigned long);
	return pvm_pkulong ( const_cast<unsigned long *>(ul), nitem, stride );
}

int
MWPvmRC::pack ( const char *str ) 
{
	bytes_packed_ += strlen(str);
	return pvm_pkstr ( const_cast<char *>(str) );
}

// BUFFER	
int
MWPvmRC::unpack ( char *bytes,         int nitem, int stride ) 
{
	int ret, setr_ret;
	
	// receive from the current active buffer
	ret = pvm_upkbyte ( bytes, nitem, stride );
	
	if ( ret==0 ) // Success!
		return 0;
	
	if ( recv_buf_list == NULL )
	       return ret;
	
	if ( ( ret == PvmNoData ) && (recv_buf_list->number()>=2) ) {
		setr_ret = next_buf();
		
		if (setr_ret!=-1) { // Switch to the next valid buffer succeeded!
			ret = pvm_upkbyte ( bytes, nitem, stride );
				
			if ( ret!=0 ) {
				MWprintf(61, "BUFFER: The first two buffers returns PvmNoData, the second can't be unpacked %d.\n", ret);
				return ret;
			} else {
				bytes_unpacked_ += nitem / stride * sizeof(char);
				return 0;
			}
		} else return setr_ret;

	} else  return ret;
}

int
MWPvmRC::unpack ( float *f,            int nitem, int stride ) 
{
	int ret, setr_ret;
	
	// receive from the current active buffer
	ret = pvm_upkfloat ( f, nitem, stride );
	
	if ( ret==0 ) // Success!
		return 0;
	
	if ( recv_buf_list == NULL )
	       return ret;
	
	if ( ( ret == PvmNoData ) && (recv_buf_list->number()>=2) ) {
		setr_ret = next_buf();
		
		if (setr_ret!=-1) { // Switch to the next valid buffer succeeded!
			ret = pvm_upkfloat ( f, nitem, stride );
				
			if ( ret!=0 ) {
				MWprintf(61, "BUFFER: The first two buffers returns PvmNoData, the second can't be unpacked %d.\n", ret);
				return ret;
			} else {
				bytes_unpacked_ += nitem / stride * sizeof(float);
				return 0;
			}
		} else return setr_ret;
	} else  return ret;
}

int
MWPvmRC::unpack ( double *d,           int nitem, int stride ) 
{
	int ret, setr_ret;
	
	// receive from the current active buffer
	ret = pvm_upkdouble ( d, nitem, stride );
	
	if ( ret==0 ) // Success!
		return 0;
	
	if ( recv_buf_list == NULL )
	       return ret;
	
	if ( ( ret == PvmNoData ) && (recv_buf_list->number()>=2) ) {
		setr_ret = next_buf();
		
		if (setr_ret!=-1) { // Switch to the next valid buffer succeeded!
			ret = pvm_upkdouble ( d, nitem, stride );
				
			if ( ret!=0 ) {
				MWprintf(61, "BUFFER: The first two buffers returns PvmNoData, the second can't be unpacked %d.\n", ret);
				return ret;
			} else {
				bytes_unpacked_ += nitem / stride * sizeof(double);
				return 0;
			}
		} else return setr_ret;
	} else  return ret;
}

int
MWPvmRC::unpack ( int *i,              int nitem, int stride ) 
{
	int ret, setr_ret;
	
	// receive from the current active buffer
	ret = pvm_upkint ( i, nitem, stride );
	
	if ( ret==0 ) // Success!
		return 0;
	
	if ( recv_buf_list == NULL )
	       return ret;
	
	if ( ( ret == PvmNoData ) && (recv_buf_list->number()>=2) ) {
		setr_ret = next_buf();
		
		if (setr_ret!=-1) { // Switch to the next valid buffer succeeded!
			ret = pvm_upkint ( i, nitem, stride );
				
			if ( ret!=0 ) {
				MWprintf(61, "BUFFER: The first two buffers returns PvmNoData, the second can't be unpacked %d.\n", ret);
				return ret;
			} else {
				bytes_unpacked_ += nitem / stride * sizeof(int);
				return 0;
			}
		} else return setr_ret;
	} else  return ret;
}

int
MWPvmRC::unpack ( unsigned int *ui,    int nitem, int stride ) 
{
	int ret, setr_ret;
	
	// receive from the current active buffer
	ret = pvm_upkuint( ui, nitem, stride );
	
	if ( ret==0 ) // Success!
		return 0;
	
	if ( recv_buf_list == NULL )
	       return ret;
	
	if ( ( ret == PvmNoData ) && (recv_buf_list->number()>=2) ) {
		setr_ret = next_buf();
		
		if (setr_ret!=-1) { // Switch to the next valid buffer succeeded!
			ret = pvm_upkuint ( ui, nitem, stride );
				
			if ( ret!=0 ) {
				MWprintf(61, "BUFFER: The first two buffers returns PvmNoData, the second can't be unpacked %d.\n", ret);
				return ret;
			} else {
				bytes_unpacked_ += nitem / stride * sizeof(unsigned int);
				return 0;
			}
		} else return setr_ret;
	} else  return ret;
}

int
MWPvmRC::unpack ( short *sh,           int nitem, int stride ) 
{
	int ret, setr_ret;
	
	// receive from the current active buffer
	ret = pvm_upkshort ( sh, nitem, stride );
	
	if ( ret==0 ) // Success!
		return 0;
	
	if ( recv_buf_list == NULL )
	       return ret;
	
	if ( ( ret == PvmNoData ) && (recv_buf_list->number()>=2) ) {
		setr_ret = next_buf();
		
		if (setr_ret!=-1) { // Switch to the next valid buffer succeeded!
			ret = pvm_upkshort ( sh, nitem, stride );
				
			if ( ret!=0 ) {
				MWprintf(61, "BUFFER: The first two buffers returns PvmNoData, the second can't be unpacked %d.\n", ret);
				return ret;
			} else {
				bytes_unpacked_ += nitem / stride * sizeof(short);
				return 0;
			}
		} else return setr_ret;
	} else  return ret;
}

int
MWPvmRC::unpack ( unsigned short *ush, int nitem, int stride ) 
{
	int ret, setr_ret;
	
	// receive from the current active buffer
	ret = pvm_upkushort ( ush, nitem, stride );
	
	if ( ret==0 ) // Success!
		return 0;
	
	if ( recv_buf_list == NULL )
	       return ret;
	
	if ( ( ret == PvmNoData ) && (recv_buf_list->number()>=2) ) {
		setr_ret = next_buf();
		
		if (setr_ret!=-1) { // Switch to the next valid buffer succeeded!
			ret = pvm_upkushort ( ush, nitem, stride );
				
			if ( ret!=0 ) {
				MWprintf(61, "BUFFER: The first two buffers returns PvmNoData, the second can't be unpacked %d.\n", ret);
				return ret;
			} else {
				bytes_unpacked_ += nitem / stride * sizeof(unsigned short);
				return 0;
			}
		} else return setr_ret;
	} else  return ret;
}

int
MWPvmRC::unpack ( long *l,             int nitem, int stride ) 
{
	int ret, setr_ret;
	
	// receive from the current active buffer
	ret = pvm_upklong ( l, nitem, stride );
	
	if ( ret==0 ) // Success!
		return 0;
	
	if ( recv_buf_list == NULL )
	       return ret;
	
	if ( ( ret == PvmNoData ) && (recv_buf_list->number()>=2) ) {
		setr_ret = next_buf();
		
		if (setr_ret!=-1) { // Switch to the next valid buffer succeeded!
			ret = pvm_upklong ( l, nitem, stride );
				
			if ( ret!=0 ) {
				MWprintf(61, "BUFFER: The first two buffers returns PvmNoData, the second can't be unpacked %d.\n", ret);
				return ret;
			} else {
				bytes_unpacked_ += nitem / stride * sizeof(long);
				return 0;
			}
		} else return setr_ret;
	} else  return ret;
}

int
MWPvmRC::unpack ( unsigned long *ul,   int nitem, int stride ) 
{
	int ret, setr_ret;
	
	// receive from the current active buffer
	ret = pvm_upkulong ( ul, nitem, stride );
	
	if ( ret==0 ) // Success!
		return 0;
	
	if ( recv_buf_list == NULL )
	       return ret;
	
	if ( ( ret == PvmNoData ) && (recv_buf_list->number()>=2) ) {
		setr_ret = next_buf();
		
		if (setr_ret!=-1) { // Switch to the next valid buffer succeeded!
			ret = pvm_upkulong ( ul, nitem, stride );
				
			if ( ret!=0 ) {
				MWprintf(61, "BUFFER: The first two buffers returns PvmNoData, the second can't be unpacked %d.\n", ret);
				return ret;
			} else {
				bytes_unpacked_ += nitem / stride * sizeof(unsigned long);
				return 0;
			}
		} else return setr_ret;
	} else  return ret;
}

int
MWPvmRC::unpack ( char *str ) 
{
	int ret, setr_ret;
	
	// receive from the current active buffer
	ret = pvm_upkstr ( str );
	
	if ( ret==0 ) // Success!
		return 0;
	
	if ( recv_buf_list == NULL )
	       return ret;
	
	if ( ( ret == PvmNoData ) && (recv_buf_list->number()>=2) ) {
		setr_ret = next_buf();
		
		if (setr_ret!=-1) { // Switch to the next valid buffer succeeded!
			ret = pvm_upkstr ( str );
				
			if ( ret!=0 ) {
				MWprintf(61, "BUFFER: The first two buffers returns PvmNoData, the second can't be unpacked %d.\n", ret);
				return ret;
			} else {
				bytes_unpacked_ += strlen(str);
				return 0;
			}
		} else return setr_ret;
	} else  return ret;
}

int 
MWPvmRC::read_RMstate ( FILE *fp = NULL )
{
    // No Pvm specific functions
    return 0;
}


int 
MWPvmRC::write_RMstate ( FILE *fp = NULL )
{
    // No Pvm specific functions
    return 0;
}

#ifdef NWSENABLED
typedef int (*callBackFcnDefn)(void* arg, int status);
void setRestartFunction ( 
				callBackFcnDefn /*cbRestart*/, 
				void* 			/*cbRestartArg*/, 
				callBackFcnDefn /*cbCkpt*/, 
				void* 			/*cbCkptArg*/ )
{
}

#endif
