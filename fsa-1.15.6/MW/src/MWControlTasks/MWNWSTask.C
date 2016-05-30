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
#include "MWNWSTask.h"
#include "MWSystem.h"

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <sys/socket.h>
#include <signal.h>
#include <unistd.h>
#include <errno.h>
#include <strings.h>


MWNWSTask::MWNWSTask ( int nu, double interval, int length, int port, char *addr )
{
	taskType = MWNWS;
	timeInterval = interval;
	iterations = length;
	min = max = median = -1;
	portNo = port;
	number = nu;
	strcpy ( machineAddress, addr );
}

MWNWSTask::MWNWSTask ( )
{
	taskType = MWNWS;
	min = max = median = -1;
	portNo = 7;
	number = -1;
}

MWNWSTask::~MWNWSTask()
{
}

void 
MWNWSTask::pack_work( void )
{
	RMC->pack ( &timeInterval, 1, 1 );
	RMC->pack ( &iterations, 1, 1 );
	RMC->pack ( &portNo, 1, 1 );
	RMC->pack ( machineAddress );
}

void 
MWNWSTask::unpack_work( void )
{
	RMC->unpack ( &timeInterval, 1, 1 );
	RMC->unpack ( &iterations, 1, 1 );
	RMC->unpack ( &portNo, 1, 1 );
	RMC->unpack ( machineAddress );
}

void 
MWNWSTask::pack_results( void )
{
	RMC->pack ( &max, 1, 1 );
	RMC->pack ( &min, 1, 1 );
	RMC->pack ( &median, 1, 1 );
}
    
void 
MWNWSTask::unpack_results( void )
{
	double temp;
	MWprintf ( 10, "NWS Unpacking results\n");
	RMC->unpack ( &temp, 1, 1 );
	MWprintf ( 10, "%f\n", temp );
	max = temp;
	RMC->unpack ( &min, 1, 1 );
	RMC->unpack ( &median, 1, 1 );
}

void 
MWNWSTask::pack_subresults( int )
{
}
    
void 
MWNWSTask::unpack_subresults( int )
{
}

void 
MWNWSTask::printself(int level)
{
	MWprintf ( level, "Interval:%f Iterations:%d PortNum:%d \n", timeInterval, 
								iterations, portNo );
	if ( max > 0 )
		MWprintf ( level, "results:- min:%f max:%f median:%f\n",
						min, max, median );
}

void 
MWNWSTask::write_ckpt_info( FILE *fp )
{
	fprintf ( fp, "%f %d ", timeInterval, iterations );
	if ( max >= 0 )
		fprintf ( fp, "%f %f %f", max, min, median );
	else
		fprintf ( fp, "-1 " );
}

void 
MWNWSTask::read_ckpt_info( FILE *fp )
{
	fscanf ( fp, "%lf %d ", &timeInterval, &iterations );
	fscanf ( fp, "%lf ", &max );
	if ( max >= 0 )
		fscanf ( fp, "%lf %lf ", &min, &median );
}

void
NWSControlTask ( MWTask *task )
{
	int listenFd;
	struct sockaddr_in servaddr;
	struct sockaddr_in cliaddr;
	MWNWSTask *nwsTask = (MWNWSTask *)task;
	char sendData = 'a';
	char recvData;
	double startTime, endTime;
	double tempStart, tempEnd;
	int cc;
	int iter = 0;
	int temp;

	MWprintf (10, "In NWSControlTask\n");
	
	listenFd = socket( AF_INET, SOCK_DGRAM, 0);
	if ( listenFd < 0 )
	{
		MWprintf( 10, "TCP socket creation error. %d\n", errno);
		::exit(1);
	}

	bzero( &cliaddr, sizeof(cliaddr));
	cliaddr.sin_family = AF_INET;
	cliaddr.sin_addr.s_addr = htonl(INADDR_ANY);
	cliaddr.sin_port = htons(0);


	if ( bind( listenFd, (struct sockaddr *) &cliaddr, sizeof(cliaddr)) < 0)
	{
		MWprintf( 10, "TCP socket bind error %d.\n", errno);
		::exit(2);
	}
  
	bzero( &servaddr, sizeof(servaddr));
	servaddr.sin_family = AF_INET;
	servaddr.sin_addr.s_addr = inet_addr( nwsTask->machineAddress);
	servaddr.sin_port = htons(nwsTask->portNo);

	startTime = MWSystem::gettimeofday();
	nwsTask->min = -1;
	nwsTask->max = -1;
	nwsTask->median = -1;
	do
	{
		tempStart = MWSystem::gettimeofday();
		cc = sendto ( listenFd, &sendData, 1, 0, (struct sockaddr *)&servaddr, sizeof(servaddr) );
		if ( cc < 0 )
		{
			MWprintf ( 10, "Could not sendto %d\n", errno);
			// XXX Instead of 
			// 	::exit(4);
			// which will cause a worker goes down (because the master's portNo is not open.
			// 
			// We probably just need to ignore the error. 
			
		}

		temp = sizeof(servaddr);	// XXX Added because error when checkpointing
		cc = recvfrom ( listenFd, &recvData, 1, 0, (struct sockaddr *)&servaddr,
						(socklen_t *)&temp);
		if ( cc < 0 )
		{
			MWprintf ( 10, "Could not recvfrom %d\n", errno);
			// XXX Instead of 
			// 	::exit(5);
			// which will cause a worker goes down (because the master's portNo is not open.
			// 
			// We probably just need to ignore the error. 
		}
		tempEnd = MWSystem::gettimeofday();
		iter++;

		if ( nwsTask->min < 0 || tempEnd - tempStart < nwsTask->min )
			nwsTask->min = tempEnd - tempStart;
		if ( nwsTask->max < 0 || tempEnd - tempStart > nwsTask->max )
			nwsTask->max = tempEnd - tempStart;
	} while ( iter < nwsTask->iterations );

	endTime = MWSystem::gettimeofday();
	nwsTask->median = ( endTime - startTime ) / nwsTask->iterations;

	close(listenFd);
}

MWTask*
NWSNewTask ( )
{
	return new MWNWSTask ( );
}
