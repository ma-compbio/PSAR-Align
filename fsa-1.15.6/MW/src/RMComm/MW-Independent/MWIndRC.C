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
#include "MWIndRC.h"
#include <MW.h>
#include <MWDriver.h>
#include <MWTask.h>
#include <MWWorker.h>

MWRMComm * GlobalRMComm = new MWIndRC ( );
MWRMComm * MWTask::RMC = GlobalRMComm;
MWRMComm * MWDriver::RMC = GlobalRMComm;
MWRMComm * MWWorker::RMC = GlobalRMComm;

MWIndRC::MWIndRC() 
{
    target_num_workers = -1;
    num_arches = -1;
    // memset ( worker_requirements, 0, (16 * 1024 ));
    // num_arches = 1;
}
		/// Destructor...
MWIndRC::~MWIndRC()
{
	/// Nothing to delete
}

void 
MWIndRC::hostadd ( )
{
    last_tag = HOSTADD;
    return;
}

void 
MWIndRC::exit( int exitval ) 
{ 
    // ::exit ( exitval ); 
}

int 
MWIndRC::setup ( int argc, char *argv[], int *my_id, int *master_id ) 
{
    *my_id = 2;
    *master_id = 0;
    return 0; 
};


int 
MWIndRC::config( int *nhosts, int *narches, MWWorkerID ***w ) 
{ 
	// Nothing to configure either
    return 0; 
}
	
int 
MWIndRC::start_worker ( MWWorkerID *w ) 
{ 
    if (!w) return -1; 
    w->set_id1 ( 2 );
    w->set_arch ( 0 );
    w->set_exec_class(0); 	// XXX should have done this!
    return 0;
};

int 
MWIndRC::init_beginning_workers ( int *nworkers, MWWorkerID ***workers ) 
{
    return 0;
};

int
MWIndRC::restart_beginning_workers( int *nworkers, MWWorkerID ***workers, MWmessages msg ){
	return init_beginning_workers( nworkers, workers);
}

int 
MWIndRC::removeWorker( MWWorkerID *w ) 
{
    return 0;
};

int 
MWIndRC::read_RMstate( FILE *fp )
{
    return 0;
}

int 
MWIndRC::write_RMstate ( FILE *fp )
{ 
    return 0;
}

int 
MWIndRC::hostaddlogic( int *num_workers )
{
    return 0; 
}

int 
MWIndRC::initsend ( int encoding )
{ 
    ind_sendbuf.num_sent = 0;
    return 0;
}

int 
MWIndRC::send ( int to_whom, int msgtag ) 
{ 
    last_tag = msgtag;
	sent_data = true; 
    return 0;
}

int 
MWIndRC::recv ( int from_whom, int msgtag ) 
{ 
    memcpy ( &ind_recvbuf.buf, &ind_sendbuf.buf, INDP_BUF_SIZE * sizeof(char) ); 
    ind_recvbuf.num_sent = 0;
	sent_data = false;
    return 0;
}

int
MWIndRC::nrecv ( int from_whom, int msgtag )
{
	if ( !sent_data )
	{
		last_tag = NO_MESSAGE;
		sent_data = false;
		return 0;
	}
	else	
	{
		recv ( from_whom, msgtag );
	}
	return 0;
}

int 
MWIndRC::bufinfo ( int buf_id, int *len, int *tag, int *from ) 
{ 
    *tag = last_tag; 
    return 0;
}

void
MWIndRC::who ( int *wh )
{
        unpack ( wh, 1, 1 );
}


int 
MWIndRC::pack ( const char *bytes,         int nitem, int stride )
{
    int num_sent = ind_sendbuf.num_sent;
    for ( int i = 0; i < nitem; i++ )
    {
	// ASSERT ( num_sent + sizeof ( char ) <= 1024 * 64 );
	memcpy ( &(ind_sendbuf.buf[num_sent]), &bytes[ i * stride ], sizeof(char) );
	ind_sendbuf.num_sent += sizeof(char);
	num_sent += sizeof(char);
    }
    return 0;
}

int 
MWIndRC::pack ( const float *f,            int nitem, int stride )
{
    int num_sent = ind_sendbuf.num_sent;
    for ( int i = 0; i < nitem; i++ )
    {
	// ASSERT ( num_sent + sizeof ( float ) <= 1024 * 64 );
	memcpy ( &(ind_sendbuf.buf[num_sent]), &f[ i * stride ], sizeof(float) );
	ind_sendbuf.num_sent += sizeof(float);
	num_sent += sizeof(float);
    }
    return 0;
}

int 
MWIndRC::pack ( const double *d,           int nitem, int stride )
{
    int num_sent = ind_sendbuf.num_sent;
    for ( int i = 0; i < nitem; i++ )
    {
	// ASSERT ( num_sent + sizeof ( double ) <= 1024 * 64 );
	memcpy ( &(ind_sendbuf.buf[num_sent]), &d[ i * stride ], sizeof(double) );
	ind_sendbuf.num_sent += sizeof(double);
	num_sent += sizeof(double);
    }
    return 0;
}

int 
MWIndRC::pack ( const int *in,              int nitem, int stride )
{
    int num_sent = ind_sendbuf.num_sent;
    for ( int i = 0; i < nitem; i++ )
    {
	// ASSERT ( num_sent + sizeof ( int ) <= 1024 * 64 );
	memcpy ( &(ind_sendbuf.buf[num_sent]), &in[ i * stride ], sizeof(int) );
	ind_sendbuf.num_sent += sizeof(int);
	num_sent += sizeof(int);
    }
    return 0;
}
		
int 
MWIndRC::pack ( const unsigned int *ui,    int nitem, int stride )
{
    int num_sent = ind_sendbuf.num_sent;
    for ( int i = 0; i < nitem; i++ )
    {
	// ASSERT ( num_sent + sizeof ( unsigned int ) <= 1024 * 64 );
	memcpy ( &(ind_sendbuf.buf[num_sent]), &ui[ i * stride ], sizeof(unsigned int) );
	ind_sendbuf.num_sent += sizeof(unsigned int);
	num_sent += sizeof(unsigned int);
    }
    return 0;
}

int 
MWIndRC::pack ( const short *sh,           int nitem, int stride )
{
    int num_sent = ind_sendbuf.num_sent;
    for ( int i = 0; i < nitem; i++ )
    {
	// ASSERT ( num_sent + sizeof ( short ) <= 1024 * 64 );
	memcpy ( &(ind_sendbuf.buf[num_sent]), &sh[ i * stride ], sizeof(short) );
	ind_sendbuf.num_sent += sizeof(short);
	num_sent += sizeof(short);
    }
    return 0;
}

int 
MWIndRC::pack ( const unsigned short *ush, int nitem, int stride )
{
    int num_sent = ind_sendbuf.num_sent;
    for ( int i = 0; i < nitem; i++ )
    {
	// ASSERT ( num_sent + sizeof ( unsigned short ) <= 1024 * 64 );
	memcpy ( &(ind_sendbuf.buf[num_sent]), &ush[ i * stride ], sizeof(unsigned short) );
	ind_sendbuf.num_sent += sizeof(unsigned short);
	num_sent += sizeof(unsigned short);
    }
    return 0;
}

int 
MWIndRC::pack ( const long *l,             int nitem, int stride )
{
    int num_sent = ind_sendbuf.num_sent;
    for ( int i = 0; i < nitem; i++ )
    {
	// ASSERT ( num_sent + sizeof ( long ) <= 1024 * 64 );
	memcpy ( &(ind_sendbuf.buf[num_sent]), &l[ i * stride ], sizeof(long) );
	ind_sendbuf.num_sent += sizeof(long);
	num_sent += sizeof(long);
    }
    return 0;
}

int 
MWIndRC::pack ( const unsigned long *ul,   int nitem, int stride )
{
    int num_sent = ind_sendbuf.num_sent;
    for ( int i = 0; i < nitem; i++ )
    {
	// ASSERT ( num_sent + sizeof ( unsigned long ) <= 1024 * 64 );
	memcpy ( &(ind_sendbuf.buf[num_sent]), &ul[ i * stride ], sizeof(unsigned long) );
	ind_sendbuf.num_sent += sizeof(unsigned long);
	num_sent += sizeof(unsigned long);
    }
    return 0;
}

int 
MWIndRC::pack ( const char *string )
{
    int num_sent = ind_sendbuf.num_sent;
    // ASSERT ( num_sent + sizeof(int) + sizeof( char ) * strlen (string) <= 1024 * 64 );
    int len = strlen(string);
    memcpy ( &(ind_sendbuf.buf[num_sent]), &len, sizeof(int) );
    ind_sendbuf.num_sent += sizeof(int);
    num_sent += sizeof(int);
    memcpy ( &(ind_sendbuf.buf[num_sent]), string, len * sizeof(char) );
    ind_sendbuf.num_sent += sizeof(char) * len;
    num_sent += sizeof(char) * len;
    return 0;
}

int 
MWIndRC::unpack ( char *bytes,         int nitem, int stride )
{
    int num_sent = ind_recvbuf.num_sent;
    for ( int i = 0; i < nitem; i++ )
    {
	char temp;
	memcpy ( &temp,&(ind_recvbuf.buf[num_sent]), sizeof(char) );
	ind_recvbuf.num_sent += sizeof(char);
	num_sent += sizeof(char);
	bytes[i * stride] = temp;
    }
    return 0;
}

int 
MWIndRC::unpack ( float *f,            int nitem, int stride )
{
    int num_sent = ind_recvbuf.num_sent;
    for ( int i = 0; i < nitem; i++ )
    {
	float temp;
	memcpy ( &temp, &(ind_recvbuf.buf[num_sent]), sizeof(float) );
	ind_recvbuf.num_sent += sizeof(float);
	num_sent += sizeof(float);
	f[i * stride] = temp;
    }
    return 0;
}

int 
MWIndRC::unpack ( double *d,           int nitem, int stride )
{
    int num_sent = ind_recvbuf.num_sent;
    for ( int i = 0; i < nitem; i++ )
    {
	double temp;
	memcpy ( &temp, &(ind_recvbuf.buf[num_sent]), sizeof(double) );
	ind_recvbuf.num_sent += sizeof(double);
	num_sent += sizeof(double);
	d[i * stride] = temp;
    }
    return 0;
}

int 
MWIndRC::unpack ( int *i,              int nitem, int stride )
{
    int num_sent = ind_recvbuf.num_sent;
    for ( int ii = 0; ii < nitem; ii++ )
    {
	int temp;
	memcpy ( &temp, &(ind_recvbuf.buf[num_sent]), sizeof(int) );
	ind_recvbuf.num_sent += sizeof(int);
	num_sent += sizeof(int);
	i[ii * stride] = temp;
    }
    return 0;
}

int 
MWIndRC::unpack ( unsigned int *ui,    int nitem, int stride )
{
    int num_sent = ind_recvbuf.num_sent;
    for ( int i = 0; i < nitem; i++ )
    {
	unsigned int temp;
	memcpy ( &temp, &(ind_recvbuf.buf[num_sent]), sizeof(unsigned int) );
	ind_recvbuf.num_sent += sizeof(unsigned int);
	num_sent += sizeof(unsigned int);
	ui[i * stride] = temp;
    }
    return 0;
}

int 
MWIndRC::unpack ( short *sh,           int nitem, int stride )
{
    int num_sent = ind_recvbuf.num_sent;
    for ( int i = 0; i < nitem; i++ )
    {
	short temp;
	memcpy ( &temp, &(ind_recvbuf.buf[num_sent]), sizeof(short) );
	ind_recvbuf.num_sent += sizeof(short);
	num_sent += sizeof(short);
	sh[i * stride] = temp;
    }
    return 0;
}

int 
MWIndRC::unpack ( unsigned short *ush, int nitem, int stride )
{
    int num_sent = ind_recvbuf.num_sent;
    for ( int i = 0; i < nitem; i++ )
    {
	unsigned short temp;
	memcpy ( &temp, &(ind_recvbuf.buf[num_sent]), sizeof(unsigned short) );
	ind_recvbuf.num_sent += sizeof(unsigned short);
	num_sent += sizeof(unsigned short);
	ush[i * stride] = temp;
    }
    return 0;
}

int 
MWIndRC::unpack ( long *l,             int nitem, int stride )
{
    int num_sent = ind_recvbuf.num_sent;
    for ( int i = 0; i < nitem; i++ )
    {
	long temp;
	memcpy ( &temp, &(ind_recvbuf.buf[num_sent]), sizeof(long) );
	ind_recvbuf.num_sent += sizeof(long);
	num_sent += sizeof(long);
	l[i * stride] = temp;
    }
    return 0;
}

int 
MWIndRC::unpack ( unsigned long *ul,   int nitem, int stride )
{
    int num_sent = ind_recvbuf.num_sent;
    for ( int i = 0; i < nitem; i++ )
    {
	unsigned long temp;
	memcpy ( &temp, &(ind_recvbuf.buf[num_sent]), sizeof(unsigned long) );
	ind_recvbuf.num_sent += sizeof(unsigned long);
	num_sent += sizeof(unsigned long);
	ul[i * stride] = temp;
    }
    return 0;
}

int 
MWIndRC::unpack ( char *string )
{
    int num_sent = ind_recvbuf.num_sent;
    int len;
    memcpy ( &len, &(ind_recvbuf.buf[num_sent]), sizeof(int) );
    ind_recvbuf.num_sent += sizeof(int);
    num_sent += sizeof(int);
    memcpy ( string, &(ind_recvbuf.buf[num_sent]), len * sizeof(char) );
    ind_recvbuf.num_sent += sizeof(char) * len;
    num_sent += sizeof(char) * len;
    string[len] = '\0';
    return 0;
}
//BUFFER
int 
MWIndRC::recv_all( int from_whom, int msgtag )
{
	MWprintf(91, "Not implemented yet for MW-Ind.\n");
	return 0;
}

int 
MWIndRC::setrbuf( int bid )
{
	MWprintf(91, "Not implemented yet for MW-Ind.\n");
	return -1;
}

int 
MWIndRC::freebuf( int bid )
{
	MWprintf(91, "Not implemented yet for MW-Ind.\n");
	return -1;
}

MWList<void> * 
MWIndRC::recv_buffers()
{
	MWprintf(91, "Not implemented yet for MW-Ind.\n");
	return NULL;
}

int 
MWIndRC::next_buf()
{
	MWprintf(91, "Not implemented yet for MW-Ind.\n");
	return -1;
}
	
