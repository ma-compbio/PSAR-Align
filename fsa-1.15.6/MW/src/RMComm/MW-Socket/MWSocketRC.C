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

#ifdef WINDOWS

#define FD_SETSIZE 1024
#include <winsock2.h>
#define close closesocket
#define popen _popen
#define pclose _pclose
#define sleep _sleep
#define EINTR WSAEINTR
#define errno WSAGetLastError()
#else 

#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <netdb.h>
#include <unistd.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <signal.h>
#include <netinet/in.h>

#endif

#include <ctype.h>

#ifdef __svr4__
#include <strings.h>
#endif

#include "MWSocketRC.h"
#include <MWDriver.h>
#include <MWWorkerID.h>
#include <MWTask.h>
#include <MWWorker.h>
#include <MWSystem.h>

#define MAX_TCP_BUFFER (64*1024)

#ifdef SOCKET_MASTER
MWRMComm * MWDriver::RMC = new MWSocketRC( TRUE, 0 );
MWRMComm * MWTask::RMC = MWDriver::RMC;
MWRMComm * MWWorker::RMC = NULL;
#else
MWRMComm * MWWorker::RMC = new MWSocketRC ( FALSE, 0 );
MWRMComm * MWTask::RMC = MWWorker::RMC;
MWRMComm * MWDriver::RMC = NULL;
#endif

extern int *MW_exec_class_num_workers;

#ifdef USE_POLL
#include <sys/poll.h>
#endif

/* 
 * Instead of ifdef'ing this for platforms that don't have it
 * we have our own, so that we are always compiling it, and
 * always know that it is correct.
 */

char *mw_inet_ntop( int family, const void *addrptr, char *strptr, size_t len)
{
	const u_char *p = ( const u_char *) addrptr;

	if ( family == AF_INET )
	{
		char temp[16];
		sprintf( temp, "%d.%d.%d.%d", p[0], p[1], p[2], p[3]);
		strcpy( strptr, temp);
		return(strptr);
	}
	return( NULL);
}

int cyclePosition = 0;
bool toBlock = TRUE;
static bool timedOut;
static int  timeOutSocket = -1;

#ifndef WINDOWS
void alarmed(int signo) {
      timedOut = true;
      close(timeOutSocket);
      write(1, "Alarm timed out\n", 17);
}

void brokenPipeHandler ( int sig )
{
	MWprintf ( 10, "Got the signal %d\n", sig );
	signal ( SIGPIPE, brokenPipeHandler );
	return;
}
#endif

MWSocketRC::MWSocketRC ( bool master, int id )
{
	isMaster = master;
	if ( isMaster )
		socketId = 8997;
	else
		socketId = id;
	strcpy ( control_directory, "submit_files" );
	strcpy ( worker_number_file, "worker_number_file" );

#ifdef WINDOWS
	WORD version = MAKEWORD(2,2);
	WSADATA d;
	WSAStartup(version, &d);
#endif
	
	//by jae	
	pre_send_done = false;
	m_poll_timeout = 5; // Default value is 5 sec
}

MWSocketRC::~MWSocketRC ( )
{
	delete []hostadd_reqs;
	delete []socketWorkers;
}

int
MWSocketRC::setup ( int argc, char *argv[], int *mytid, int *mastertid )
{
	struct hostent *hptr;

	char **ptr;
	char str[16];

	subId = 0;

	if ( isMaster == TRUE )
	{
		turnNo = 0;
		*mytid = socketId;
		*mastertid = socketId;
		hostadd_reqs = NULL;
		masterId = socketId;
		masterSocket = 8997;

		char myname[256];
		MWSystem::gethostname(myname, 256);

		if ( (hptr = gethostbyname(myname)) == NULL)
		{
			MWprintf ( 10, "Cannot do a gethostbyname: errno %d\n", errno );
			::exit(1);
		}

		ptr = hptr->h_addr_list;
		strcpy ( masterAddress, ((char *)mw_inet_ntop(hptr->h_addrtype, *ptr, str, sizeof(str))) );

		socketComm = creatAndBind ( masterSocket );

		if ( listen( socketComm, 50 ) < 0 )
		{
			MWprintf ( 10, "Listen failed: errno %d\n", errno );
			close ( socketComm );
			::exit(1);
		}

		FD_ZERO ( &rfds );
		FD_SET ( socketComm, &rfds );

		for ( int i = 0; i < 4096; i++ )
		{
			reverseMap[i] = -1;
		}

		// If the directory doesn't exist make it, ignoring the error
		// if it does exist

		char buf[256];
		sprintf(buf, "mkdir %s", control_directory);
		system(buf);

		// Windows doesn't raise signals on broken pipe writes
#ifndef WINDOWS
		signal ( SIGPIPE, brokenPipeHandler );
#endif
	}
	else
	{
		struct sockaddr_in myAddr;

		*mastertid = atoi ( argv[2] );
		socketId = atoi ( argv[1] );
		masterSocket = atoi ( argv[3] );
		*mytid = socketId;
		strcpy ( masterAddress, argv[4] );
		socketComm = creatAndBind ( 0 );

		memset(&myAddr, 0, sizeof(myAddr));
		myAddr.sin_family = AF_INET;
		(myAddr.sin_addr).s_addr = inet_addr(masterAddress);
		myAddr.sin_port = htons(masterSocket);

		int retryCount = 0;

tryAgain:
		if ( connect ( socketComm, (struct sockaddr *)&myAddr, sizeof(myAddr) ) < 0 )
		{
			MWprintf ( 10, "Connect failed with the master: errno %d\n", errno );
			if (retryCount++ < 24) {
				sleep(5);
				goto tryAgain;
			}
			::exit(1);
		}
		MWprintf ( 10, "Connect succeeded\n");

		char buf[sizeof(int)];
		marshall_int ( socketId, buf );
		if ( finalSend ( socketComm, buf, sizeof(int), 0 ) < 0 )
		{
			MWprintf ( 10, "Initial Send failed to the master: errno %d\n", errno );
			::exit(1);
		}
	}

	return 0;
}

void
MWSocketRC::exit ( int retval )
{
	if ( isMaster )
	{
		killWorkers ( );
		closeSockets ( );
	}else {
		// by jae
		::exit ( retval );
	}
}

int
MWSocketRC::init_beginning_workers ( int *nworkers, MWWorkerID ***workers )
{
	int i;

	MWprintf ( 10, "In MWSocketRC::init_beginning_workers()\n");
	hostadd_reqs = new int[num_arches];

	for ( i = 0; i < num_arches; i++ )
	{
		hostadd_reqs[i] = 0;
	}


	socketWorkers = new SocketWorker[target_num_workers];
	for ( i = 0; i < target_num_workers; i++ )
	{
		socketWorkers[i].state = MWSocket_FREE;
		socketWorkers[i].socket = -1;
	}
	*nworkers = 0;

	return 0;
}

int
MWSocketRC::restart_beginning_workers ( int *nworkers, MWWorkerID ***workers, MWmessages msg )
{
	return init_beginning_workers ( nworkers, workers );
}

void
MWSocketRC::who ( int *wh )
{
	unpack ( wh, 1, 1 );
}

int
MWSocketRC::start_worker ( MWWorkerID *w )
{
	w->set_id1 ( whomRecv );
	w->set_id2 ( whomRecv );
	w->set_arch ( socketWorkers[whomRecv].arch );
	w->set_exec_class ( socketWorkers[whomRecv].exec_class );

	hostadd_reqs[socketWorkers[whomRecv].exec_class]--;

	return 0;
}

int
MWSocketRC::removeWorker ( MWWorkerID *w )
{
	// char buf[_POSIX_PATH_MAX];

	if ( !w ) 
	{
		MWprintf (10, "Worker ID cannot be NULL in removeWorker\n" );
		return -1;
	}

	return killWorker ( w->get_id1() );
}

int
MWSocketRC::hostaddlogic ( int *numworkers )
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
	int cur;
	int req;
	sort_exec_class_ratio ( sorted_order );
	int status;
                        
	MWprintf ( 60, "hal: target: %d\n", target_num_workers );

	int HOSTINC = hostinc_;
	for ( i = 0; i < exec_classes; i++ )
	{
		cur = sorted_order[i];
		MWprintf(60, "hal:  cur: %d, num_workers: %d, hostadd_reqs: %d, target: %d\n", cur, 
			 MW_exec_class_num_workers[cur], hostadd_reqs[cur], 
			 exec_class_target_num_workers[cur]);
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

	delete [] sorted_order;
	return 0;
}

void
MWSocketRC::sort_exec_class_ratio ( int *temp )
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


int
MWSocketRC::read_RMstate ( FILE *myfp )
{
	char temp[_POSIX_PATH_MAX];
	int cluster, subcluster;
	FILE *fp;

	sprintf ( temp, "%s/%s", control_directory, worker_number_file );
	fp = fopen ( temp, "r" );
	if ( fp )
	{
		while ( fscanf ( fp, "%d %d", &cluster, &subcluster ) > 0 )
		{
			char condor_cmd[ 3 * _POSIX_PATH_MAX];
			sprintf ( condor_cmd, "condor_rm %d.%d", cluster, subcluster );
			if ( system ( condor_cmd ) < 0 )
			{
				MWprintf ( 10, "Couldn't remove the job %d.%d\n", cluster, subcluster );
			}
		}
		fclose ( fp );
	}
	return 0;
}

int
MWSocketRC::write_RMstate ( FILE *myfp )
{
	FILE *fp;
	char temp[_POSIX_PATH_MAX];
	char cmd[3 * _POSIX_PATH_MAX];

	sprintf ( temp, "%s/%s", control_directory, worker_number_file );
#ifdef WINDOWS
	sprintf ( cmd, "copy temp %s\\%s", control_directory, worker_number_file );
	printf("Before rename cmd is %s\n", cmd);
#else
	sprintf ( cmd, "/bin/mv temp %s/%s", control_directory, worker_number_file );
#endif
	fp = fopen ( "temp", "w" );
	if ( !fp )
	{
		MWprintf ( 10, "DANGER!!, MW was not able to write into file temp\n" );
		return -1;
	}

	for ( int i = 0; i < target_num_workers; i++ )
	{
		if ( socketWorkers[i].state == MWSocket_FREE )
			continue;

		fprintf ( fp, "%d %d\n", socketWorkers[i].cId, socketWorkers[i].subcId );
	}

	fclose ( fp );

	if ( system ( cmd ) < 0 )
	{
		MWprintf ( 10, "Couldn't mv file temp into the actual file\n");
		return -1;
	}

	return 0;
}

char*
MWSocketRC::process_executable_name ( char *exec, int ex_cl, int ar_cl )
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


	MWprintf(20, "MWSocketRC::process_executable_name %s %d %d\n", exec, ex_cl, ar_cl);
	// Use unsigned int to make gcc -Wall happy

	// Change everything to capitals, except in cases...
	// 1) u, x newstring[i-1] == '4' ?
	// 2) k == 'u' || k == 'x' and newstring[i+1] == '8' (handles arch==x86_64)
	for ( unsigned ii = 0; ii < strlen ( newstring ); ii++ )
	{
		int k = newstring[ii];
		if( ! ( ( k == 'u' || k == 'x' ) && 
			( (newstring[ii-1] == '4') || (newstring[ii+1] == '8') ) ) )
		{
			k = toupper(k);
		}
		newstring[ii] = k;
	}

    i = 0;
    index = strstr ( newstring, "OPSYS" );
	
	if (index == NULL) {
		MWprintf(10, "ERROR: OPSYS not specified in add_executable requirements string\n");
		MWprintf(10, "newstring is %s\n", newstring);
		exit(-1);
	}

    index += 5;
    while ( *index == ' ' || *index == '=' || *index == '"' || *index == '\\' ) index++;
    while ( *index && *index != ')' && *index != '"' && *index != '\\' ) opsys[i++] = *index++;
    opsys[i] = '\0';

    i = 0;
    index = strstr ( newstring, "ARCH" );
	if (index == NULL) {
		MWprintf(10, "ERROR: ARCH not specified in add_executable requirements string\n");
		exit(-1);
	}
    index += 4;
    while ( *index == ' ' || *index == '=' || *index == '"' || *index == '\\' ) index++;
    while ( *index && *index != ')' && *index != '"' && *index != '\\' ) arch[i++] = *index++;
    arch[i] = '\0';

    sprintf ( newexec, "%s.%s.%s.exe", newexec, opsys, arch );
    unlink ( newexec );

    MWprintf ( 10, "Making a link from %s to %s\n", exec, newexec );
#if defined(WINDOWS) || defined(CYGWIN)
    sprintf ( exe, "copy %s %s", exec, newexec );
#else
    sprintf ( exe, "ln -s %s %s", exec, newexec );
#endif
    system ( exe );

    return newexec;
}

int
MWSocketRC::initsend ( int useless )
{
	messageSize = 0;
	pre_send_done = false;
	return 0;
}

int
MWSocketRC::send ( int toWhom, int msgtag )
{
	char buffer[3*sizeof(int)];
	int retVal;
	int sendingSocket;

	if ( isMaster )
		sendingSocket = socketWorkers[toWhom].socket;
	else {
		// by jae
		// Worker
		if( pre_send_done ) {
			return 0;
		}
		sendingSocket = socketComm;
	}

	marshall_int ( msgtag, buffer );
	marshall_int ( messageSize, buffer + sizeof(msgtag) );
	marshall_int ( socketId, buffer + sizeof(msgtag) + sizeof(messageSize) );

	retVal = finalSend ( sendingSocket, buffer, 3 * sizeof(int), 0 );
	if ( retVal < 0 ) 
		return retVal;
	retVal = finalSend ( sendingSocket, sendBuffer, messageSize, 0 );

	return retVal;
}

//by jae
int
MWSocketRC::pre_send (int msgtag )
{
	if( pre_send_done ) {
		return 0;
	}

	if ( isMaster ) {
		MWprintf ( 10, "ERROR: pre_send is not implemented for Master\n");
		return -1;
	}

	int retVal = send(0, msgtag);
	pre_send_done = true;

	return retVal;
}

int
MWSocketRC::finalSend ( int s, char *buf, size_t len, int flags )
{
#ifdef __svr4__
	ssize_t retVal;
#else
	int retVal; 
#endif


	do
	{
		retVal = ::send ( s, buf, len, flags );
	} while ( (retVal == -1) && (errno == EINTR ));

	if (retVal == -1) {
		MWprintf(10, "FinalSend::send returned -1 errno %d\n", errno);
		return -1;
	}

	if ( retVal < (int)len ) // Make gcc -Wall happy
	{
		return finalSend(s, buf + retVal, len - retVal, flags);
	}

	return retVal;
}

int
MWSocketRC::recv ( int fromWhom, int msgtag )
{
	recvPointer = 0;
	if ( isMaster )
		return masterrecv ( fromWhom, msgtag );
	else
		return workerrecv ( fromWhom, msgtag );
}

int
MWSocketRC::masterrecv ( int fromWhom, int msgtag )
{
	char buffer[ 3 * sizeof(int) ];

	// if fromWhom == -1, receive from anyone
    if ( fromWhom <= 0 ) {
		return masterrecvany();
	}

	// Otherwise, from one particular someone
	{
		// to recv from someone.
		int retVal = finalRecv ( socketWorkers[fromWhom].socket, buffer, 3 * sizeof(int), 0 );
		if ( retVal < 0 )
			return retVal;

		while ( msgtag > 0 && unmarshall_int ( buffer ) != msgtag )
		{
			retVal = finalRecv ( socketWorkers[fromWhom].socket, recvBuffer,
							unmarshall_int ( &buffer[sizeof(int)] ), 0 );
			if ( retVal < 0 )
				return retVal;

			retVal = finalRecv ( socketWorkers[fromWhom].socket, buffer, 3 * sizeof(int), 0 );
			if ( retVal < 0 )
				return retVal;
		}

		retVal = finalRecv ( socketWorkers[fromWhom].socket, recvBuffer,
							unmarshall_int ( &buffer[sizeof(int)] ), 0 );

		recvSize = retVal;
		msgTag = unmarshall_int ( buffer );
		whomRecv = fromWhom;
		return retVal;
	}
}

#ifdef USE_POLL

int
MWSocketRC::masterrecvany()
{
	struct pollfd pollfds[target_num_workers + 1];

	int index = 0;
	for (int i = 0; i < target_num_workers; i++) {
		if (socketWorkers[i].socket > 0) {
			if( socketWorkers[i].state == MWSocket_THREAD ) {
				// skip sockets attacehd to threads
				continue;
			}
			pollfds[index].fd = socketWorkers[i].socket;
			pollfds[index].events = POLLIN; 
			index++;
		}
	}

	pollfds[index].fd = socketComm;
	pollfds[index].events = POLLIN; 

	index++;

	int r;
	do {
			errno = 0;
			//jae
			//r = poll(pollfds, index, -1);
			r = poll(pollfds, index, m_poll_timeout*1000);
	} while (errno == EINTR);

	if (r == -1) {
		MWprintf(10, "poll returned -1, errno = %d\n", errno);
		return -1;
	}

	//by jae
	if( r == 0 ) {
		// timeout
		MWprintf(10, "poll timeout\n");
		return -1;
	}

	// For all the non-accept'ing sockets
	for (int j = 0; j < index - 1; j++) {
		if (pollfds[j].revents & POLLHUP) {
			// It doesn't look like linux returns this on a socket, but
			// just in case it does...
			return handle_killed_worker(reverseMap[pollfds[j].fd]);
		}
		if (pollfds[j].revents & POLLIN) {
			return handle_finished_worker(reverseMap[pollfds[j].fd]);
		}

		if (pollfds[j].revents != 0) {
			MWprintf(10, "Unexpected return from poll for fd %d: %x\n", pollfds[j].fd, pollfds[j].revents);
		}
	}

	if (pollfds[index - 1].revents) {
		return handle_starting_worker();
	}


	MWprintf(10, "USE_POLL::masterrecvany shouldn't get here\n");
	return -1;

}
#endif

#ifndef USE_POLL
int 
MWSocketRC::masterrecvany()
{

		int retVal;
		int i;
		FD_ZERO ( &rfds );
		FD_SET ( socketComm, &rfds );

		for ( i = 0; i < FD_SETSIZE; i++ )
		{
			int worker_id = reverseMap[i];
			if ( worker_id >= 0 )
			{
				if( socketWorkers[worker_id].state == MWSocket_THREAD ) {
					// skip sockets attacehd to threads
					continue;
				}
				FD_SET ( i, &rfds );
			}
		}


		struct timeval timeout;
		timeout.tv_usec = 0;
		timeout.tv_sec = m_poll_timeout;

		do
		{
			retVal = select ( FD_SETSIZE, &rfds, NULL, NULL, &timeout );
		} while ((retVal == -1) && ( errno == EINTR ));

		//by jae
		if( retVal == 0 ) {
			// timeout
			MWprintf(10, "select timeout\n");
			return -1;
		}

		if ( retVal < 0 )
		{
			MWprintf ( 10, "Select returned -1, errno = %d\n", errno );
			//exit(1);
			// don't need to exit if select returns zero, just remove the worker that went down
            fd_set temprfds;
            FD_ZERO ( &rfds );
            for ( i = 0; i < FD_SETSIZE; i++ )
            {
                int worker_id = reverseMap[i];
                if ( worker_id >= 0 )
                {
					if( socketWorkers[worker_id].state == MWSocket_THREAD ) {
						// skip sockets attacehd to threads
						continue;
					}

                    struct timeval zero;
                    zero.tv_sec=0;
                    zero.tv_usec=0;
                    FD_ZERO ( &temprfds );
                    FD_SET ( i, &temprfds );
                    int selRet;
                    do
                    {
                        selRet = select ( FD_SETSIZE, &temprfds, NULL, NULL, &zero );
                    } while ( errno == EINTR );
                    if( selRet < 0 )
                    {
                        int len = handle_killed_worker ( reverseMap[i] );
                        reverseMap[i] = -1;
                        return len;
                    }
                }
			}
		}
		if ( retVal == 0 )
		{
			MWprintf ( 10, "Unexpected 0 return from select, exitting\n" );
			exit(1);
		}

		return checkSockets ( );
}

#endif
int
MWSocketRC::workerrecv ( int fromWhom, int msgtag )
{
	// int temp;
	int retVal;
	char buffer [ 3 * sizeof(int) ];

	retVal = finalRecv ( socketComm, buffer, 3 * sizeof(int), 0 );
	if ( retVal < 0 )
		return retVal;

	while ( msgtag > 0 && unmarshall_int( buffer ) != msgtag )
	{
		retVal = finalRecv ( socketComm, recvBuffer, unmarshall_int ( &buffer[sizeof(int)] ), 0 );
		if ( retVal < 0 )
			return retVal;

		retVal = finalRecv ( socketComm, buffer, 3 * sizeof(int), 0 );
		if ( retVal < 0 )
			return retVal;
	}

	retVal = finalRecv ( socketComm, recvBuffer, unmarshall_int( &buffer[sizeof(int)] ), 0 );

	recvSize = retVal;
	msgTag = unmarshall_int ( buffer );
	whomRecv = fromWhom;
	return retVal;
}

int
MWSocketRC::nrecv ( int fromWhom, int msgtag )
{
    recvPointer = 0;
    if( isMaster )
	{
		MWprintf(10, "Master should not call this\n");
		exit(-1);
	}
	int retVal;
	char buffer [ 3 * sizeof(int) ];

	retVal = nfinalRecv ( socketComm, buffer, 3 * sizeof(int), 0 );
	if ( retVal == -2 ) // nothing in the socket, return from call
	{
		recvSize = 0;//retVal;
		msgTag = NO_MESSAGE;
		whomRecv = fromWhom;
        return 0;
	}
	while ( msgtag > 0 && unmarshall_int( buffer ) != msgtag )
	{
		retVal = finalRecv ( socketComm, recvBuffer, unmarshall_int ( &buffer[sizeof(int)] ), 0 );
		if ( retVal < 0 )
			return retVal;

		retVal = finalRecv ( socketComm, buffer, 3 * sizeof(int), 0 );
		if ( retVal < 0 )
			return retVal;
	}
	retVal = finalRecv ( socketComm, recvBuffer, unmarshall_int( &buffer[sizeof(int)] ), 0 );
	recvSize = retVal;
	msgTag = unmarshall_int ( buffer );
	whomRecv = fromWhom;
	return retVal;
}

int
MWSocketRC::nfinalRecv ( int s, void *buf, size_t len, int flags )
{
	int retVal;

	fd_set sfd;
	FD_ZERO( &sfd );
	FD_SET( s, &sfd );
	struct timeval zero;
	zero.tv_sec=0; zero.tv_usec=0;
	int selRet = select(s+1, &sfd, NULL, NULL, &zero);
	if( selRet == 0 )
		return -2;

	do
	{
		retVal = ::recv ( s, (char *)buf, len, MSG_PEEK );
/*  
if(retVal == -1 && errno == EAGAIN){
            fcntl(s, F_SETFL, fdflags);
            MWprintf(51,"Nothing in socket, returning\n");
            return -1;
            }
*/
	} while ( errno == EINTR || retVal < (int)len ); // Make -Wall happy

	do
	{
        retVal = ::recv ( s, (char *)buf, len, flags );
	} while ((retVal == -1) && ( errno == EINTR ));

	if ( retVal < (int)len ) // Make gcc -Wall happy
	{
		MWprintf ( 10, "Couldn't properly do a recv from socket %d: retVal = %d, errno = %d\n", s, retVal, errno );
		return -1;
	}

	return retVal;
}

int
MWSocketRC::finalRecv ( int s, char *b, size_t len, int flags )
{
	int retVal; 

	char *buf = (char *) b;

	int total = 0;
	retVal = 0;

	timedOut = false;
	timeOutSocket = s;

	do
	{

			// First time through, these are no-ops
		len -= retVal;
		buf += retVal;

		retVal = ::recv ( s, buf, len, flags );

		if( (timedOut == true) || (retVal < 0 && errno != EINTR)) {
			MWprintf(10, "recv returned -1: errno is %d\n", errno);
			return -1;
		}

        if ( retVal == 0 ) return 0; // Socket closed

		if (errno == EINTR) {
			retVal = 0;
		}
		total += retVal;
		
	} while (( errno == EINTR ) || (retVal < (int)len)) ;

	if ( retVal < (int)len ) // Make gcc -Wall happy
	{
		return 0;
	}

	return total;
}

int
MWSocketRC::creatAndBind ( int port )
{
	int retVal;
	int bufSize = MAX_TCP_BUFFER;
	struct sockaddr_in myAddr;

	retVal = socket ( AF_INET, SOCK_STREAM, 0 );
	if ( retVal < 0 )
	{
		MWprintf ( 10, "Couldn't open a socket errno %d\n", errno );
		::exit(1);
	}
	if ( setsockopt( retVal, SOL_SOCKET, SO_SNDBUF, (char *)&bufSize, sizeof(int)) != 0)
	{
		MWprintf ( 10, "Couldn't set the send buffer size: errno %d\n", errno );
		close ( retVal );
		::exit(1);
	}
	if ( setsockopt( retVal, SOL_SOCKET, SO_RCVBUF, (char *)&bufSize, sizeof(int)) != 0)
	{
		MWprintf ( 10, "Couldn't set the recv buffer size: errno %d\n", errno );
		close ( retVal );
		::exit(1);
	}
	memset(&myAddr, 0, sizeof(myAddr));
	myAddr.sin_family = AF_INET;
	(myAddr.sin_addr).s_addr = htonl(INADDR_ANY);
	myAddr.sin_port = htons(port);

	int range = 50; // range of ports to scan
	int range_counter = 0;
	while( bind( retVal, (struct sockaddr *) &myAddr, sizeof(myAddr) ) < 0 && range_counter < range)
	{
		//MWprintf ( 10, "Couldn't bind the socket to port: errno %d\n", errno );
		port++;
		myAddr.sin_port = htons(port);
		range_counter++;
		//close ( retVal );
		//::exit(1);
	}
	if ( range_counter >= range )
	{
		MWprintf ( 10, "Couldn't bind the socket to port: errno %d\n", errno );
		close ( retVal );
		::exit(1);
	}
	else
		MWprintf ( 10, "Socket bound to port: %d\n", port );
	if ( isMaster == TRUE )
		masterSocket = port;

	return retVal;
}

int
MWSocketRC::checkSockets ( )
{
	char buf;
	int len = 1;
	// First we check the read list
	for ( int i = 0; i < FD_SETSIZE; i++ )
	{
		if ( i == socketComm )
			continue;

		if ( FD_ISSET ( i, &rfds ) )
		{
			cyclePosition = i;
			toBlock = FALSE;
			FD_CLR ( i, &rfds );

			//by jae
			int worker_id = reverseMap[i];
			if( worker_id >= 0 ) {
				if( socketWorkers[worker_id].state == MWSocket_THREAD ) {
					continue;
				}
			}

			int retVal;
			do
			{
				retVal = ::recv ( i, &buf, len, MSG_PEEK );
			} while ( (retVal < 0) && (errno == EINTR ));
			if ( errno != 0 )
			{
				len = handle_killed_worker ( reverseMap[i] );
				reverseMap[i] = -1;
				return len;
			}
			else if ( retVal <= 0 )
			{
				len = handle_killed_worker ( reverseMap[i] );
                reverseMap[i] = -1;
                return len;
			}
			else
			{
				return handle_finished_worker ( reverseMap[i] );
			}
		}
	}

	// Finally let us examine the socketComm.
	if ( FD_ISSET ( socketComm, &rfds ) )
	{
		// someone is trying to connect to us.
		FD_CLR ( socketComm, &rfds );
		return handle_starting_worker ( );
	}
	cyclePosition = 0;
	toBlock = TRUE;
	return -1;
}

int
MWSocketRC::handle_starting_worker ( )
{
	int temp;
	int remote;
	char buf[sizeof(int)];

	temp = accept ( socketComm, 0, 0 ); 
	if ( temp < 0 )
	{
		MWprintf ( 10, "Accept failed: errno %d\n", errno );
		return -1;
	}

#ifndef WINDOWS
	signal(SIGALRM, alarmed);
	alarm(5);
#endif
	if ( finalRecv ( temp, buf, sizeof(int), 0 ) < 0 )
	{
		MWprintf ( 10, "First recv from the new guy failed: errno %d\n", errno );
#ifndef WINDOWS
	alarm(0);
#endif
		close ( temp );
		return -1;
	}
#ifndef WINDOWS
	alarm(0);
#endif
	remote = unmarshall_int ( buf );

	if ( !&socketWorkers[remote] )
	{
		MWprintf ( 10, "How can the socket worker structure be null from a worker that is starting\n");
		return -1;
	}

    if ( socketWorkers[remote].socket != -1 )
    {  
		close(temp);

		// What's happened here is that the worker has died and been restarted by
		// Condor, but the socket connection hasn't gone away.  This can happen!
		// Since we don't have user-level keep-alives, if a machine's kernel panics,
		// or loses power, it is possible for Condor to realize the machine has gone
		// away before we do.  Since we can only return one value to the higher level,
		// we close the new socket (which will make the new worker restart again), 
		// and mark the old worker dead to the upper level. 

        int len = handle_killed_worker ( remote );
        if ( reverseMap[socketWorkers[remote].socket] == remote ) // if mapping of socket to worker not already taken by another worker, invalidate
            reverseMap[socketWorkers[remote].socket] = -1;
        MWprintf(10, "worker %d tried to come back from eviction, killed it\n",remote);
        return len;
   }

	socketWorkers[remote].state = MWSocket_EXECUTING;
	socketWorkers[remote].socket = temp;
	if ( socketWorkers[remote].socket < 0 )
	{
		MWprintf ( 10, "Accept failed: errno %d\n", errno );
		return -1;
	}
	FD_SET ( socketWorkers[remote].socket, &rfds );
	reverseMap[socketWorkers[remote].socket] = remote;

	whomRecv = remote;
	msgTag = MWSocketMessageTags[MWSocketHostAdd];
	recvSize = 0;
	return 0;
}

int
MWSocketRC::handle_finished_worker ( int i )
{
	// int temp;
	int retVal;
	char buffer [ 3 * sizeof(int) ];

	//by jae
	if( i < 0 ) {
		return -1;
	}

	if ( !&socketWorkers[i] )
	{
		MWprintf ( 10, "How can the socket worker structure be null from a worker that is executing\n");
		return -1;
	}

	//by jae
	if( socketWorkers[i].state == MWSocket_THREAD ) {
		return -1;
	}

#ifndef WINDOWS
	signal(SIGALRM, alarmed);
	alarm(5);
#endif
	retVal = finalRecv ( socketWorkers[i].socket, buffer, 3 * sizeof(int), 0 );
	if ( retVal < 0 ) {
#ifndef WINDOWS
	alarm(0);
#endif
		return retVal;
	}
#ifndef WINDOWS
	alarm(0);
#endif

	if (retVal == 0) {
		MWprintf(10, "read from worker returned 0, killing off disconnected worker\n");
		return handle_killed_worker(i);
	}
	retVal = finalRecv ( socketWorkers[i].socket, recvBuffer, unmarshall_int ( &buffer[sizeof(int)] ), 0 );

	recvSize = retVal;
	msgTag = unmarshall_int ( buffer ) ;
	whomRecv = i;
	return retVal;
}

int
MWSocketRC::handle_killed_worker ( int i )
{
	MWprintf ( 10, "A worker is killed %d\n", i );
	//by jae
	if( i < 0 ) {
		return -1;
	}
	if ( !&socketWorkers[i] )
	{
		MWprintf ( 10, "How can the socket worker structure be null from a worker that is executing\n");
		return -1;
	}

	marshall_int ( i, recvBuffer );
	recvSize = sizeof(int);
	msgTag = MWSocketMessageTags[MWSocketTaskExit];
	whomRecv = i;

	killWorker(i);
	return 0;
}

int
MWSocketRC::bufinfo ( int buf_id, int *len, int *tag, int *sending_host )
{
	*sending_host = whomRecv;
	*tag = msgTag;
	*len = recvSize;
	return 0;
}

int
MWSocketRC::do_spawn ( int nWorkers, int ex_cl )
{
	int cID = -1;
	int i;

	char sub_file[_POSIX_PATH_MAX];
	char sub_file1[_POSIX_PATH_MAX];
	char temp[65536];
	char exe[_POSIX_PATH_MAX];
	FILE *ptr;
	char requirements[10 * _POSIX_PATH_MAX]; // 10 is max num executables

	if ( nWorkers <= 0 ) return 0;

	int *ids = new int[nWorkers];

	sprintf( sub_file1, "submit_file.%d", subId );
	sprintf( sub_file, "./%s/submit_file.%d", control_directory, subId++ );
	sprintf( exe, "condor_submit ./%s/%s", control_directory, sub_file1 );


	FILE *f = Open( sub_file, "w" );
	if ( f == NULL )
	{
		MWprintf ( 10, "Couldn't open the submitfile for writing\n");
		return -1;
	}
    	
	int index = 0;
	for ( i = 0; i < target_num_workers && index < nWorkers; i++ ) 
	{
		if ( socketWorkers[i].state == MWSocket_FREE ) 
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

	fprintf ( f, "Universe = Vanilla\n");
	fprintf( f, "Executable = mw_exec%d.$$(Opsys).$$(Arch).exe\n", ex_cl );

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
	//strcat ( requirements, "&& (Machine==\"quasar-16.biostat.wisc.edu\")");
	strcat ( requirements, " )" );
	for ( i = 0; i < nWorkers; i++ ) 
	{
		fprintf( f, "arguments = %d %d %d %s\n", ids[i], masterId, masterSocket, masterAddress );
		fprintf( f, "log = %s/log_file\n", control_directory );

		fprintf( f, "Output = output_file.%d\n", ids[i] );
		fprintf( f, "Error = error_file.%d\n", ids[i] );
						       
		fprintf( f, "Requirements = %s\n ", requirements );

		fprintf( f, "should_transfer_files = Yes\n"); 
		fprintf( f, "when_to_transfer_output = ON_EXIT\n");

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


			// send worker MWprintf back immediately, don't wait
			// for worker exit.  May want to turn this off for
			// better performance, but this is really handy for debugging

		// This causes windows clients to hang, so leave off by default
		//fprintf( f, "stream_output = true\n");
		//fprintf( f, "stream_error = true\n");

			// Have Condor give us workers sorted by MIPS
		fprintf( f, "rank = Mips\n");
		
			// If a job starts, then exits before connecting to the
			// socket master, we want condor to restart the job.
			// otherwise, the master doen't know the job has left
			// the queue, and will wait indefinitely for it to start
			// only need this for Socket
		fprintf( f, "on_exit_remove = false\n");
		fprintf( f, "Queue\n" );
	}
	fclose( f );

	MWprintf (10, "About to call condor_submit %s\n", exe);
	ptr = popen ( exe, "r" );
	if ( ptr != NULL ) 
	{
		while ( fscanf(ptr, "%s", temp ) >= 0 ) 
		{

			if ( strcmp ( temp, "cluster" ) == 0 ) 
			{

				fscanf( ptr, "%s", temp );
				cID = atoi ( temp );
				MWprintf ( 10, "Spawned to cluster %d\n", cID );
				break;
			}
		}
		if ( ptr ) pclose ( ptr ); 
		for ( i = 0; i < nWorkers; i++ ) 
		{
			int tempp = ids[i];
			socketWorkers[tempp].state = MWSocket_SUBMITTED;
			socketWorkers[tempp].cId = cID;
			socketWorkers[tempp].subcId = i;
			socketWorkers[tempp].arch = -1;
			socketWorkers[tempp].socket = -1;
			socketWorkers[tempp].exec_class = ex_cl;
		}
	}
	else 
	{
		MWprintf ( 10, "Couldn't popen in FileRC\n");
		return -1;
	}

	delete []ids;
	write_RMstate ( NULL );
	return nWorkers;
}

int
MWSocketRC::pack ( const char *bytes, int nitem, int stride )
{
	if ( messageSize + nitem * sizeof(char) > MWSOCKET_MAX_MSG_SIZE )
	{
		MWprintf ( 10, "ERROR: Maximum Message size exceeded, (%d > %d), cannot pack\n",
			   messageSize + nitem * sizeof(char), MWSOCKET_MAX_MSG_SIZE);
		return -1;
	}

	bytes_packed_ += nitem / stride;
	for ( int i = 0; i < nitem; i++ )
	{
		memcpy ( &sendBuffer[messageSize], &(bytes [ i * stride ]), sizeof(char) );
		messageSize += sizeof(char);
	}

	return 0;
}

int
MWSocketRC::pack ( const float *f, int nitem, int stride )
{
	if ( messageSize + nitem * sizeof(float) > MWSOCKET_MAX_MSG_SIZE )
	{
		MWprintf ( 10, "ERROR: Maximum Message size exceeded, (%d > %d), cannot pack\n",
			   messageSize + nitem * sizeof(float), MWSOCKET_MAX_MSG_SIZE);
		return -1;
	}

	bytes_packed_ += nitem / stride * sizeof(float);
	for ( int i = 0; i < nitem; i++ )
	{
		marshall_float ( f[ i * stride ], &sendBuffer[messageSize] );
		messageSize += sizeof(float);
	}

	return 0;
}

int
MWSocketRC::pack ( const double *d, int nitem, int stride )
{
	if ( messageSize + nitem * sizeof(double) > MWSOCKET_MAX_MSG_SIZE )
	{
		MWprintf ( 10, "ERROR: Maximum Message size exceeded, (%d > %d), cannot pack\n",
			   messageSize + nitem * sizeof(double), MWSOCKET_MAX_MSG_SIZE);
		return -1;
	}

	bytes_packed_ += nitem / stride * sizeof(double);
	for ( int i = 0; i < nitem; i++ )
	{
		marshall_double ( d [ i * stride ], &sendBuffer[messageSize] );
		messageSize += sizeof(double);
	}

	return 0;
}

int
MWSocketRC::pack ( const int *ii, int nitem, int stride )
{
	if ( messageSize + nitem * sizeof(int) > MWSOCKET_MAX_MSG_SIZE )
	{
		MWprintf ( 10, "ERROR: Maximum Message size exceeded, (%d > %d), cannot pack\n",
			   messageSize + nitem * sizeof(int), MWSOCKET_MAX_MSG_SIZE);
		return -1;
	}

	bytes_packed_ += nitem / stride * sizeof(int);
	for ( int i = 0; i < nitem; i++ )
	{
		marshall_int ( ii [ i * stride ],  &sendBuffer[messageSize] );
		messageSize += sizeof(int);
	}

	return 0;
}

int
MWSocketRC::pack ( const unsigned int *ui, int nitem, int stride )
{
	if ( messageSize + nitem * sizeof(unsigned int) > MWSOCKET_MAX_MSG_SIZE )
	{
		MWprintf ( 10, "ERROR: Maximum Message size exceeded, (%d > %d), cannot pack\n",
			   messageSize + nitem * sizeof(unsigned int), MWSOCKET_MAX_MSG_SIZE);
		return -1;
	}


	bytes_packed_ += nitem / stride * sizeof(unsigned int);
	for ( int i = 0; i < nitem; i++ )
	{
		marshall_uint ( ui [ i * stride ], &sendBuffer[messageSize] );
		messageSize += sizeof(unsigned int);
	}

	return 0;
}

int
MWSocketRC::pack ( const short *sh, int nitem, int stride )
{
	if ( messageSize + nitem * sizeof(short) > MWSOCKET_MAX_MSG_SIZE )
	{
		MWprintf ( 10, "ERROR: Maximum Message size exceeded, (%d > %d), cannot pack\n",
			   messageSize + nitem * sizeof(short), MWSOCKET_MAX_MSG_SIZE);
		return -1;
	}

	bytes_packed_ += nitem / stride * sizeof(short);
	for ( int i = 0; i < nitem; i++ )
	{
		marshall_short ( sh [ i * stride ], &sendBuffer[messageSize] );
		messageSize += sizeof(short);
	}

	return 0;
}

int
MWSocketRC::pack ( const unsigned short *us, int nitem, int stride )
{
	if ( messageSize + nitem * sizeof(unsigned short) > MWSOCKET_MAX_MSG_SIZE )
	{
		MWprintf ( 10, "ERROR: Maximum Message size exceeded, (%d > %d), cannot pack\n",
			   messageSize + nitem * sizeof(unsigned short), MWSOCKET_MAX_MSG_SIZE);
		return -1;
	}

	bytes_packed_ += nitem / stride * sizeof(unsigned short);
	for ( int i = 0; i < nitem; i++ )
	{
		marshall_ushort ( us [ i * stride ], &sendBuffer[messageSize] );
		messageSize += sizeof(unsigned short);
	}

	return 0;
}

int
MWSocketRC::pack ( const long *l, int nitem, int stride )
{
   
	if ( messageSize + nitem * sizeof(long) > MWSOCKET_MAX_MSG_SIZE )
	{
		MWprintf ( 10, "ERROR: Maximum Message size exceeded, (%d > %d), cannot pack\n",
			   messageSize + nitem * sizeof(long), MWSOCKET_MAX_MSG_SIZE);
		return -1;
	}

	bytes_packed_ += nitem / stride * sizeof(long);
	for ( int i = 0; i < nitem; i++ )
	{
		marshall_long ( l [ i * stride ], &sendBuffer[messageSize] );
		messageSize += sizeof(long);
	}

	return 0;
}

int
MWSocketRC::pack ( const unsigned long *ul, int nitem, int stride )
{
	if ( messageSize + nitem * sizeof(unsigned long) > MWSOCKET_MAX_MSG_SIZE )
	{
		MWprintf ( 10, "ERROR: Maximum Message size exceeded, (%d > %d), cannot pack\n",
			   messageSize + nitem * sizeof(unsigned long), MWSOCKET_MAX_MSG_SIZE);
		return -1;
	}

	bytes_packed_ += nitem / stride * sizeof(unsigned long);
	for ( int i = 0; i < nitem; i++ )
	{
		marshall_ulong ( ul [ i * stride ], &sendBuffer[messageSize] );
		messageSize += sizeof(unsigned long);
	}

	return 0;
}

int
MWSocketRC::pack ( const char *str )
{
	int len = strlen ( str );

	if ( messageSize + sizeof(int) * len > MWSOCKET_MAX_MSG_SIZE )
	{
		MWprintf ( 10, "ERROR: Maximum Message size exceeded, (%d > %d), cannot pack\n",
			   messageSize + len * sizeof(int), MWSOCKET_MAX_MSG_SIZE);
		return -1;
	}

	bytes_packed_ += strlen(str);
	marshall_int ( len, &sendBuffer[messageSize] );
	messageSize += sizeof(int);
	memcpy ( &sendBuffer[messageSize], str, len );
	messageSize += len;

	return 0;
}

// by jae
// Send raw data directly
int MWSocketRC::raw_pack( char* buf, int len ) 
{
	// Currently, this function works only for Worker
	if( isMaster ) {
		MWprintf ( 10, "ERROR: raw_pack is not implemented for Master\n");
		return -1;
	}

	return finalSend ( socketComm, buf, len, 0 );
}

void MWSocketRC::thread_start(int worker_id) 
{
	socketWorkers[worker_id].state = MWSocket_THREAD;
}

void MWSocketRC::thread_stop(int worker_id) 
{
	if( socketWorkers[worker_id].state == MWSocket_THREAD ) {
		socketWorkers[worker_id].state = MWSocket_EXECUTING;
	}
}

/// Read raw data from Worker ID
int MWSocketRC::raw_unpack( int from_whom, char* buf, int len ) 
{
	if( !isMaster ) {
		MWprintf ( 10, "ERROR: raw_unpack is not implemented for Worker\n");
		return -1;
	}
	if ( from_whom < 0 ) {
		MWprintf ( 10, "ERROR: raw_unpack has invalid Worker ID\n");
		return -1;
	}

	return finalRecv ( socketWorkers[from_whom].socket, buf, len, 0 );
}

int
MWSocketRC::unpack ( char *bytes, int nitem, int stride )
{
	if ( int(recvPointer + nitem * sizeof(char)) > recvSize )
	{
		MWprintf ( 10, "ERROR: You can't unpack more than what's received\n");
		return -1;
	}

	bytes_unpacked_ += nitem / stride * sizeof(char);
	for ( int i = 0; i < nitem; i++ )
	{
		memcpy ( &(bytes [ i * stride ]), &recvBuffer[recvPointer], sizeof(char) );
		recvPointer += sizeof(char);
	}

	return 0;
}

int
MWSocketRC::unpack ( float *f, int nitem, int stride )
{
	if ( int(recvPointer + nitem * sizeof(float)) > recvSize )
	{
		MWprintf ( 10, "ERROR: You can't unpack more than what's received\n");
		return -1;
	}

	bytes_unpacked_ += nitem / stride * sizeof(float);
	for ( int i = 0; i < nitem; i++ )
	{
		f [ i * stride ] = unmarshall_float ( &recvBuffer[recvPointer] );
		recvPointer += sizeof(float);
	}

	return 0;
}

int
MWSocketRC::unpack ( double *d, int nitem, int stride )
{
	if ( int (recvPointer + nitem * sizeof(double)) > recvSize )
	{
		MWprintf ( 10, "ERROR: You can't unpack more than what's received\n");
		return -1;
	}

	bytes_unpacked_ += nitem / stride * sizeof(double);
	for ( int i = 0; i < nitem; i++ )
	{
		 d [ i * stride ] = unmarshall_double ( &recvBuffer[recvPointer] );
		recvPointer += sizeof(double);
	}

	return 0;
}

int
MWSocketRC::unpack ( int *ii, int nitem, int stride )
{
	if ( int (recvPointer + nitem * sizeof(int)) > recvSize ) // Make gcc -Wall happy
	{
		MWprintf ( 10, "ERROR: You can't unpack more than what's received\n");
		return -1;
	}

	bytes_unpacked_ += nitem / stride * sizeof(int);
	for ( int i = 0; i < nitem; i++ )
	{
		ii [ i * stride ] = unmarshall_int ( &recvBuffer[recvPointer] );
		recvPointer += sizeof(int);
	}

	return 0;
}

int
MWSocketRC::unpack ( unsigned int *ui, int nitem, int stride )
{
	if ( int (recvPointer + nitem * sizeof(unsigned int)) > recvSize ) // Make gcc -Wall happy
	{
		MWprintf ( 10, "ERROR: You can't unpack more than what's received\n");
		return -1;
	}

	bytes_unpacked_ += nitem / stride * sizeof(unsigned int);
	for ( int i = 0; i < nitem; i++ )
	{
		ui [ i * stride ] = unmarshall_uint ( &recvBuffer[recvPointer] );
		recvPointer += sizeof(unsigned int);
	}

	return 0;
}

int
MWSocketRC::unpack ( short *sh, int nitem, int stride )
{
	if ( int (recvPointer + nitem * sizeof(short)) > recvSize )
	{
		MWprintf ( 10, "ERROR: You can't unpack more than what's received\n");
		return -1;
	}

	bytes_unpacked_ += nitem / stride * sizeof(short);
	for ( int i = 0; i < nitem; i++ )
	{
		sh [ i * stride ] = unmarshall_short ( &recvBuffer[recvPointer] );
		recvPointer += sizeof(short);
	}

	return 0;
}

int
MWSocketRC::unpack ( unsigned short *us, int nitem, int stride )
{
	if ( int (recvPointer + nitem * sizeof(unsigned short)) > recvSize )
	{
		MWprintf ( 10, "ERROR: You can't unpack more than what's received\n");
		return -1;
	}

	bytes_unpacked_ += nitem / stride * sizeof(unsigned short);
	for ( int i = 0; i < nitem; i++ )
	{
		us [ i * stride ] = unmarshall_ushort ( &recvBuffer[recvPointer] );
		recvPointer += sizeof(unsigned short);
	}

	return 0;
}

int
MWSocketRC::unpack ( long *l, int nitem, int stride )
{
	if ( int (recvPointer + nitem * sizeof(long)) > recvSize )
	{
		MWprintf ( 10, "ERROR: You can't unpack more than what's received\n");
		return -1;
	}

	bytes_unpacked_ += nitem / stride * sizeof(long);
	for ( int i = 0; i < nitem; i++ )
	{
		l [ i * stride ] = unmarshall_long ( &recvBuffer[recvPointer] );
		recvPointer += sizeof(long);
	}

	return 0;
}

int
MWSocketRC::unpack ( unsigned long *ul, int nitem, int stride )
{
	if ( int (recvPointer + nitem * sizeof(unsigned long)) > recvSize )
	{
		MWprintf ( 10, "ERROR: You can't unpack more than what's received\n");
		return -1;
	}

	bytes_unpacked_ += nitem / stride * sizeof(unsigned long);
	for ( int i = 0; i < nitem; i++ )
	{
		ul [ i * stride ] = unmarshall_ulong ( &recvBuffer[recvPointer] );
		recvPointer += sizeof(unsigned long);
	}

	return 0;
}

int
MWSocketRC::unpack ( char *str )
{
	int len;

	if ( int (recvPointer + sizeof(int)) > recvSize )
	{
		MWprintf ( 10, "ERROR: You can't unpack more than what's received\n");
		return -1;
	}

	len = unmarshall_int ( &recvBuffer[recvPointer] );
	recvPointer += sizeof(int);

	if ( recvPointer + len > recvSize )
	{
		MWprintf ( 10, "ERROR: You can't unpack more than what's received\n");
		return -1;
	}

	memcpy ( str, &recvBuffer[recvPointer], len );
	str[len] = '\0';
	recvPointer += len;

	bytes_unpacked_ += strlen(str);

	return 0;
}

int
MWSocketRC::killWorkers ( )
{
	for ( int i = 0; i < target_num_workers; i++ )
	{
		if ( socketWorkers[i].state == MWSocket_SUBMITTED || socketWorkers[i].state == MWSocket_EXECUTING || socketWorkers[i].state == MWSocket_THREAD) 
		{
			char exe[ 8 * _POSIX_PATH_MAX];

			sprintf ( exe, "condor_rm %d.%d", socketWorkers[i].cId, socketWorkers[i].subcId );
			if ( system ( exe ) < 0 )
			{
				MWprintf ( 10, "ERROR: Couldnt kill the worker %d, condor job %d.%d\n", i, socketWorkers[i].cId, socketWorkers[i].subcId );
			}
		}
	}

	return 0;
}

void
MWSocketRC::closeSockets ( )
{
	close ( socketComm );
	return;
}

int
MWSocketRC::killWorker ( int i )
{
	//by jae
	if( i < 0 ) {
		return -1;
	}

	close ( socketWorkers[i].socket );
	hostadd_reqs[socketWorkers[i].exec_class]++;

    socketWorkers[i].state = MWSocket_SUBMITTED;
	socketWorkers[i].socket = -1;
    write_RMstate ( NULL );
	return 0;
}

inline void
MWSocketRC::marshall_int ( int n, char *buf )
{
	unsigned int k = htonl ( n );
	memcpy ( buf, &k, sizeof(int) );
	return;
}

inline void
MWSocketRC::marshall_uint ( unsigned int n, char *buf )
{
	unsigned int k = htonl ( n );
	memcpy ( buf, &k, sizeof(unsigned int) );
	return;
}

inline void
MWSocketRC::marshall_short ( short n, char *buf )
{
	unsigned short k = htons ( n );
	memcpy ( buf, &k, sizeof(short) );
	return;
}

inline void
MWSocketRC::marshall_ushort ( unsigned short n, char *buf )
{
	unsigned short k = htons ( n );
	memcpy ( buf, &k, sizeof(unsigned short) );
	return;
}

inline void
MWSocketRC::marshall_long ( long n, char *buf )
{
	long k = htonl ( n );
	memcpy ( buf, &k, sizeof(long) );
	return;
}

inline void
MWSocketRC::marshall_ulong ( unsigned long n, char *buf )
{
	unsigned long k = htonl ( n );
	memcpy ( buf, &k, sizeof(unsigned long) );
	return;
}

inline void
MWSocketRC::marshall_float ( float n, char *buf )
{
	int k = *(int *)&n; 
	marshall_int ( k, buf );
	return;
}

inline void
MWSocketRC::marshall_double ( double n, char *buf )
{

		// Union utopia is the best flavor of babcock ice cream
	union utopia {
		double d;
		struct {
			int low;
			int high;
		} ints;
		char   bytes[8];
	};

	union utopia u;
	u.d = n;


	if ( htonl ( 3 ) == 3 )
	{
		marshall_int ( u.ints.low,  buf );
		marshall_int ( u.ints.high, buf + sizeof(int) );
	}
	else
	{
		marshall_int ( u.ints.high, buf );
		marshall_int ( u.ints.low, buf + sizeof(int) );
	}
	return;
}

inline int
MWSocketRC::unmarshall_int ( char *buf )
{
	int k;
	memcpy ( &k, buf, sizeof(int) );
	return ntohl ( k );
}

inline unsigned int
MWSocketRC::unmarshall_uint ( char *buf )
{
	unsigned int k;
	memcpy ( &k, buf, sizeof(unsigned int) );
	return ntohl ( k );
}

inline short
MWSocketRC::unmarshall_short ( char *buf )
{
	short k;
	memcpy ( &k, buf, sizeof(short) );
	return ntohs ( k );
}

inline unsigned short
MWSocketRC::unmarshall_ushort ( char *buf )
{
	unsigned short k;
	memcpy ( &k, buf, sizeof(unsigned short) );
	return ntohs ( k );
}

inline long
MWSocketRC::unmarshall_long ( char *buf )
{
	long k;
	memcpy ( &k, buf, sizeof(long) );
	return ntohl ( k );
}

inline unsigned long
MWSocketRC::unmarshall_ulong ( char *buf )
{
	unsigned long k;
	memcpy ( &k, buf, sizeof(unsigned long) );
	return ntohl ( k );
}

inline float
MWSocketRC::unmarshall_float ( char *buf )
{
	int k = unmarshall_int ( buf );
	float f = *(float *)&k;
	return f;
}

inline double
MWSocketRC::unmarshall_double ( char *buf )
{
	double d;
	int f,s;

	f = unmarshall_int ( buf );
	s = unmarshall_int ( buf + sizeof(int) );
	if ( ntohl ( 3 ) == 3 )
	{
		memcpy ( &d, &f, sizeof(int) );
		memcpy ( ((char *)&d) + sizeof(int), &s, sizeof(int) );
	}
	else
	{
		memcpy ( &d, &s, sizeof(int) );
		memcpy ( ((char *)&d) + sizeof(int), &f, sizeof(int) );
	}
	return d;
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
//BUFFER
int 
MWSocketRC::recv_all( int from_whom, int msgtag )
{
	MWprintf(91, "Not implemented yet for MW-Socket.\n");
	return 0;
}

int 
MWSocketRC::setrbuf( int bid )
{
	MWprintf(91, "Not implemented yet for MW-Socket.\n");
	return -1;
}

int 
MWSocketRC::freebuf( int bid )
{
	MWprintf(91, "Not implemented yet for MW-Socket.\n");
	return -1;
}

MWList<void>* 
MWSocketRC::recv_buffers()
{
	MWprintf(91, "Not implemented yet for MW-Socket.\n");
	return NULL;
}

int 
MWSocketRC::next_buf()
{
	MWprintf(91, "Not implemented yet for MW-Socket.\n");
	return -1;
}
