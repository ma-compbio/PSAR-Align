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
#include "MWStaticMPIRC.h"
#include <MW.h>
#include <MWDriver.h>
#include <MWTask.h>
#include <MWWorker.h>
#include <stdio.h>

MWStaticMPIRC * MWGlobalMPIRMComm = new MWStaticMPIRC();
MWRMComm * MWDriver::RMC = MWGlobalMPIRMComm;
MWRMComm * MWTask::RMC = MWGlobalMPIRMComm;
MWRMComm * MWWorker::RMC = MWGlobalMPIRMComm;

extern int *MW_exec_class_num_workers;
//#define S_PRINT

MWStaticMPIRC::MWStaticMPIRC ( )
{
	
	is_master = FALSE;

#ifdef S_PRINT
	MWprintf(10, "**** In MWStaticMPIRC -- constructor ****\n");
#endif
}	

MWStaticMPIRC::~MWStaticMPIRC ( )
{
	
}	
	
void MWStaticMPIRC::init_env(int argc, char* argv[])
{

	MPI_Init(&argc, &argv);
#ifdef S_PRINT
	MWprintf(10, "**** In MWStaticMPIRC -- init_env ****\n");
#endif
}

bool MWStaticMPIRC::is_master_proc()
{
	int myId;
   	MPI_Comm_rank(MPI_COMM_WORLD, &myId);
#ifdef S_PRINT
	MWprintf(10, "****PID = %d In MWStaticMPIRC -- is_master_proc ****\n", myId);
#endif
  	if (myId == 0)
		return true;
 	else
		return false;
}
/* Implementation of comm / rm class for mpi_finalize. */
void MWStaticMPIRC::exit( int exitval ) 
{
	MWprintf ( 50, "Before MPI_Finalize\n" );
	MPI_Finalize();
	MWprintf ( 50, "After MPI_Finalize\n" );
	::exit(exitval);
	
}	


int  
MWStaticMPIRC::setup( int argc, char *argv[], int *id, int *master_id ) 
{
	MWprintf ( 10, "In MWStaticMPIRC::setup()\n");
	//Check the process id
        MPI_Comm_rank(MPI_COMM_WORLD, id);
#ifdef S_PRINT
	MWprintf(10, "**** PID = %d In MWStaticMPIRC -- setup ****\n", id);
#endif
	
	sendBuffer = new Msg_buf();
	recvBuffer = new Msg_buf();
	//In MPI the master has pid=0
	
	*master_id = 0;

	if ( *id == 0 ) {
	   	MWprintf ( 40, "I have no parent, therefore I'm the master.\n" );
	   	is_master = TRUE;
	   	*master_id = 0;
		return 0;
   	} else {
	   	is_master = FALSE;
		return 0;
	
	}
    
	
}

int 
MWStaticMPIRC::config( int *nhosts, int *narches, MWWorkerID ***workers ) 
{
	if ( !is_master ) {
		MWprintf ( 10, "Slaves not allowed in this privaledged call!\n" );
		return -1;
	}
	
	*narches = 1;
	MWprintf ( 70, "In MWStaticMPIRC::config()\n" );
#ifdef S_PRINT
	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MWprintf(10, "**** PID = %d In MWStaticMPIRC -- config() ****\n", id);
#endif

	if ( *workers ) {
		MWprintf ( 10, "workers should be NULL when called.\n" );
		return -1;
	}
	//Numero di processi compreso il master
	int numProc;
	MPI_Comm_size(MPI_COMM_WORLD, &numProc);
	*nhosts = numProc - 1; 
	MWprintf ( 70, "In MWStaticMPIRC::config(). Number of processes are %d \n", *nhosts );
	
	
	//The number of workers is number of processes-master (nhost -1 ) 
	(*workers) = new MWWorkerID*[*nhosts];
	for ( int i=0 ; i<*nhosts ; i++ ) {
		(*workers)[i] = new MWWorkerID;
		(*workers)[i]->set_arch ( 0 );
			/* Do domething creative with speed later! */
	}

	return 0;
} 

int
MWStaticMPIRC::start_worker ( MWWorkerID *w ) 
{
	
	if ( !is_master ) 
	{
		MWprintf ( 10, "Slaves not allowed in this privaledged call!\n" );
		return -1;
	}
#ifdef S_PRINT
	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MWprintf(10, "**** PID = %d In MWStaticMPIRC -- start_worker ****\n", id);
#endif

	MWprintf ( 50, "In MWStaticMPIRC::start_worker()\n" );
	if ( !w ) 
	{
		MWprintf ( 10, "w cannot be NULL!\n" );
		return -1;
	}

		return 0;
}

int
MWStaticMPIRC::init_beginning_workers( int *nworkers, MWWorkerID ***workers ) 
{
	if ( !is_master ) 
	{
		MWprintf ( 10, "Slaves not allowed in this privaledged call!\n" );
		return -1;
	}
#ifdef S_PRINT
	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MWprintf(10, "**** PID = %d In MWStaticMPIRC -- init_beginning workers() ****\n", id);
#endif

	int i, j, narches;
	MWprintf ( 50, "In MWStaticMPIRC::init_beginning_workers\n" );

		/* config allocates memory for workers */
	if ( config( nworkers, &narches, workers ) < 0 ) 
	{
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
		//do_spawn( (*workers)[i] );
		(*workers)[i]->set_id1(i+1);
		(*workers)[i]->set_id2(0);
		(*workers)[i]->set_exec_class(0);
		(*workers)[i]->set_executable(0);
		(*workers)[i]->set_arch(0);
	}
	
	return 0;
}


int
MWStaticMPIRC::restart_beginning_workers ( int *nworkers, MWWorkerID ***workers, MWmessages msg )
{
#ifdef S_PRINT
	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MWprintf(10, "**** PID = %d In MWStaticMPIRC -- restart_beginning_workers ****\n", id);
#endif
	return init_beginning_workers ( nworkers, workers );
}


int
MWStaticMPIRC::initsend ( int encoding ) 
{
#ifdef S_PRINT
	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MWprintf(10, "**** PID = %d In MWStaticMPIRC -- initsend ****\n", id);
#endif
	sendBuffer->reset();
	sendPosition = 0;
	return 0; 
}


int
MWStaticMPIRC::send ( int to_whom, int msgtag ) 
{  
	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
#ifdef S_PRINT
	MWprintf(10, "***** PID %d : Send to %d mstgag %d size = %d ***** \n", id, to_whom, msgtag, sendBuffer->size);
#endif
   	return MPI_Send(sendBuffer->data, sendBuffer->size, MPI_PACKED, to_whom, msgtag, MPI_COMM_WORLD);
	
}

int
MWStaticMPIRC::nrecv( int from_whom, int msgtag)
{
    recvPosition = 0;	
    if (is_master)  
      return masterrecv(from_whom, msgtag);	    
    else
      return iworkerrecv(from_whom, msgtag);   
}
    
   
 
int MWStaticMPIRC::iworkerrecv(int from_whom, int msgtag)
{
   int flag, bytes;
   MPI_Status status;
   int id;
   MPI_Comm_rank(MPI_COMM_WORLD, &id);
#ifdef S_PRINT
   MWprintf(10, "**** PID %d : iworkerrecv from = %d msgtag = %d ****\n", id, from_whom, msgtag);	
#endif
   MPI_Iprobe(from_whom, msgtag, MPI_COMM_WORLD, &flag, &status); 
   if (flag <= 0) {
    return flag;
   }
   else {
     recvBuffer->reset(); 	
     MPI_Get_count( &status, MPI_PACKED, &bytes );
     recvBuffer->resize(bytes); 
     recvBuffer->sender = status.MPI_SOURCE;
     recvBuffer->tag = status.MPI_TAG;
     recvBuffer->size = bytes;   
     return MPI_Recv(recvBuffer->data, recvBuffer->size, MPI_PACKED, recvBuffer->sender, recvBuffer->tag, MPI_COMM_WORLD,&status);
     
    }

   

}

int MWStaticMPIRC::masterrecv (int from_whom, int msgtag) {
 
    
    int *buf_id = NULL;
    int bytes;
    MPI_Status status;
    buf_id = new int;
    int id;
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
#ifdef S_PRINT
    MWprintf(10, "**** PID %d : imasterrecv from = %d msgtag = %d ****\n", id, from_whom, msgtag);		 	
#endif
    (*buf_id) = MPI_Probe(from_whom, msgtag, MPI_COMM_WORLD,&status); 
    if (buf_id > 0) {
    	MPI_Get_count( &status, MPI_PACKED, &bytes );
        recvBuffer->reset();
	recvBuffer->resize(bytes);
        recvBuffer->sender = status.MPI_SOURCE;
        recvBuffer->tag = status.MPI_TAG;
        recvBuffer->size = bytes;   
        MPI_Recv(recvBuffer->data, recvBuffer->size, MPI_PACKED, recvBuffer->sender, recvBuffer->tag, MPI_COMM_WORLD,&status);
        return *buf_id;
     } 
     else
      return 0;
 
} 

 int
MWStaticMPIRC::recv( int from_whom, int msgtag)
{
	recvPosition = 0;
	if (from_whom == -1)
	   from_whom = MPI_ANY_SOURCE;
	if (msgtag == -1)
	   msgtag = MPI_ANY_TAG;	
	
	int *buf_id = new int;
	int retVal;
    	MPI_Status status;
   	int bytes, myId;
     	MPI_Comm_rank(MPI_COMM_WORLD, &myId);
#ifdef S_PRINT

	MWprintf(40, "**** PID %d : MPI_Probe in recv myId = %d from_whom = %d and tag = % d ****\n", myId, myId, from_whom, msgtag);
#endif
	retVal = MPI_Probe(from_whom, msgtag, MPI_COMM_WORLD,&status); 
        if (retVal >= 0) 
        {
        	MPI_Get_count( &status, MPI_PACKED, &bytes );
		recvBuffer->reset();
		recvBuffer->resize(bytes);
		recvBuffer->sender = status.MPI_SOURCE;
		recvBuffer->tag = status.MPI_TAG;
		recvBuffer->size = bytes; 
		MPI_Recv(recvBuffer->data, recvBuffer->size, MPI_PACKED, recvBuffer->sender, recvBuffer->tag, MPI_COMM_WORLD,&status);
		*buf_id = 1;
	 }
	 else
	 {
	  *buf_id = 0;
	 }
	
	return *buf_id;

}

int MWStaticMPIRC::bufinfo ( int buf_id, int *len, int *tag, int *from ) 
{
	
        *len = recvBuffer->size;
        *tag = recvBuffer->tag;
	*from = recvBuffer->sender;
	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
#ifdef S_PRINT

	MWprintf(40, "**** PID %d : MWStaticMPIRC::bufinfo: len %d tag %d from %d ****\n",id, *len, *tag, *from);
#endif
	return 0;

}
int 
MWStaticMPIRC::recv_all( int from_whom, int msgtag )
{

	MWprintf(91, "Not implemented yet recv_all for MW-StaticMPI.\n");
	return 0;
}



int 
MWStaticMPIRC::setrbuf( int bid )
{
	MWprintf(91, "Not implemented yet setrbuf for MW-StaticMPI.\n");
	return -1;
}

int 
MWStaticMPIRC::freebuf( int bid )
{
	MWprintf(91, "Not implemented yet freebuf for MW-StaticMPI.\n");
	return -1;
}

MWList<void>* 
MWStaticMPIRC::recv_buffers()
{
	MWprintf(91, "Not implemented yet recv_buffers for MW-StaticMPI.\n");
	return NULL;
}

int 
MWStaticMPIRC::next_buf()
{
	MWprintf(91, "Not implemented yet next_buf for MW-StaticMPI.\n");
	return -1;
}

int
MWStaticMPIRC::pack (const char *bytes, int nitem, int stride ) 
{
	bytes_packed_ += nitem / stride;
	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	sendBuffer->resize(nitem*sizeof bytes);
	MPI_Pack(const_cast<char *> (bytes), (nitem*stride), MPI_CHAR, 
		 sendBuffer->data, sendBuffer->max_size, 
		 &sendPosition, MPI_COMM_WORLD);
	sendBuffer->size = sendPosition;
#ifdef S_PRINT

	MWprintf(10, "**** PID %d : Pack char bytes = %s nitem = %d sendPosition = %d stride = %d ****\n", id, bytes, nitem, sendPosition, stride);
#endif
	return  0;
}

int
MWStaticMPIRC::pack (const float *f, int nitem, int stride  )
{
	bytes_packed_ += nitem / stride * sizeof(float);
	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	sendBuffer->resize(nitem*sizeof f);
	MPI_Pack(const_cast<float *> (f), (nitem*stride), MPI_FLOAT, 
		 sendBuffer->data, sendBuffer->max_size, 
		 &sendPosition, MPI_COMM_WORLD);
	sendBuffer->size = sendPosition;

#ifdef S_PRINT	
	MWprintf(10, "**** PID %d : Pack float data = %f nitem = %d sendPosition = %d stride = %d ****\n", id, *f, nitem, sendPosition, stride);
#endif

 	return 0; 
}


int 
MWStaticMPIRC::pack(const double *d, int nitem, int stride)
{ 
	bytes_packed_ += nitem / stride * sizeof(double);
	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	sendBuffer->resize(nitem*sizeof d);
	MPI_Pack(const_cast<double *> (d), (nitem*stride), MPI_DOUBLE, sendBuffer->data, sendBuffer->max_size, &sendPosition, MPI_COMM_WORLD);
	sendBuffer->size = sendPosition;

#ifdef S_PRINT
	MWprintf(10, "**** PID %d : Pack double sendPosition = %d nitem*stride = %d size = %d stride = %d ****\n", id,  sendPosition, nitem*stride, sendBuffer->size, stride);
#endif
	return 0; 
}


int 
MWStaticMPIRC::pack(const int *i, int nitem, int stride)
{ 
	bytes_packed_ += nitem / stride * sizeof(int);
	int id, dim;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	sendBuffer->resize(nitem*sizeof i);
        MPI_Pack(const_cast<int *> (i), (nitem*stride), MPI_INT, sendBuffer->data, 
		 sendBuffer->max_size, &sendPosition, MPI_COMM_WORLD);
	sendBuffer->size = sendPosition;

#ifdef S_PRINT
	MWprintf(10, "**** PID %d : Pack int data = %d nitem = %d sendPosition = %d stride = %d ****\n", id, *i, nitem, sendPosition, stride); 
#endif

	return 0; 
}

	
int 
MWStaticMPIRC::pack(const unsigned int *ui, int nitem, int stride)
{ 
	bytes_packed_ += nitem / stride * sizeof(unsigned int);
	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Pack(const_cast<unsigned int *> (ui), (nitem*stride), MPI_UNSIGNED, 
		 sendBuffer->data, sendBuffer->max_size, &sendPosition, MPI_COMM_WORLD);
	sendBuffer->size = sendPosition;
#ifdef S_PRINT

	MWprintf(10, "**** PID %d : Pack ui data = %ui nitem = %d sendPosition = %d stride = %d ****\n", id, *ui, nitem, sendPosition, stride); 
#endif
	return 0; 
};

int 
MWStaticMPIRC::pack(const short *sh, int nitem, int stride )
{ 	
	bytes_packed_ += nitem / stride * sizeof(short);
	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	sendBuffer->resize(nitem*sizeof sh);
	MPI_Pack(const_cast<short *> (sh), (nitem*stride), MPI_SHORT, 
		 sendBuffer->data, sendBuffer->max_size, &sendPosition, MPI_COMM_WORLD);
	sendBuffer->size = sendPosition;

#ifdef S_PRINT
	MWprintf(10, "**** PID %d : Pack sh data = %sh nitem = %d sendPosition = %d stride = %d ****\n", id, *sh, nitem, sendPosition, stride); 
#endif
	return 0; 
};

int 
MWStaticMPIRC::pack (const unsigned short *ush, int nitem, int stride  )
{ 
	bytes_packed_ += nitem / stride * sizeof(unsigned short);
	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	sendBuffer->resize(nitem*sizeof ush);
	MPI_Pack(const_cast<unsigned short *> (ush), (nitem*stride), MPI_UNSIGNED_SHORT, 
		 sendBuffer->data, sendBuffer->max_size, &sendPosition, MPI_COMM_WORLD);
	sendBuffer->size = sendPosition;

#ifdef S_PRINT
	MWprintf(10, "**** PID %d : Pack ush data = %ush nitem = %d sendPosition = %d stride = %d ****\n", id, *ush, nitem, sendPosition, stride ); 
#endif
	return 0; 
};

int 
MWStaticMPIRC::pack(const long *l, int nitem, int stride )
{ 
	bytes_packed_ += nitem / stride * sizeof(long);
	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	sendBuffer->resize(nitem*sizeof l);
	MPI_Pack(const_cast<long *> (l), (nitem*stride), MPI_LONG, sendBuffer->data, 
		 sendBuffer->max_size, &sendPosition, MPI_COMM_WORLD);
	sendBuffer->size = sendPosition;
#ifdef S_PRINT

	MWprintf(10, "**** PID %d : Pack l data = %l nitem = %d sendPosition = %d stride = %d ****\n", id, *l, nitem, sendPosition, stride); 
#endif
	return 0; 
};


int 
MWStaticMPIRC::pack(const unsigned long *ul,   int nitem, int stride )
{ 
	bytes_packed_ += nitem / stride * sizeof(unsigned long);
	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	sendBuffer->resize(nitem*sizeof ul);
	MPI_Pack(const_cast<unsigned long *> (ul), (nitem*stride), MPI_UNSIGNED_LONG, 
		 sendBuffer->data, sendBuffer->max_size, &sendPosition, MPI_COMM_WORLD);
	sendBuffer->size = sendPosition;

#ifdef S_PRINT
	MWprintf(10, "**** PID %d : Pack ul data = %ul nitem = %d sendPosition = %d ****\n", id, *ul, nitem, sendPosition); 
#endif
	return 0; 
};

int MWStaticMPIRC::pack(const char *str)
{ 
	bytes_packed_ += strlen(str);
	int nitem = strlen(str);
	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	sendBuffer->resize(nitem);

	MPI_Pack(&nitem, 1, MPI_INT, sendBuffer->data, sendBuffer->max_size, &sendPosition, MPI_COMM_WORLD);
	MPI_Pack(const_cast<char *> (str), nitem, MPI_CHAR, sendBuffer->data, 
		 sendBuffer->max_size, &sendPosition, MPI_COMM_WORLD);
	sendBuffer->size = sendPosition;

#ifdef S_PRINT	
	MWprintf(10, "**** PID %d : Pack str data = %s nitem = %d sendPosition = %d ****\n", id, str, nitem, sendPosition);
#endif
	return 0; 
};
	
int MWStaticMPIRC::unpack( char *bytes, int nitem, int stride )
{ 
	
	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Unpack(recvBuffer->data, recvBuffer->size, &recvPosition, bytes, nitem, MPI_CHAR, MPI_COMM_WORLD);
#ifdef S_PRINT

	MWprintf(10, "**** PID %d : Unpack bytes data = %s nitem = %d recvPosition = %d ****\n", id, bytes, (nitem*stride), recvPosition);
#endif
	bytes_unpacked_ += nitem / stride * sizeof(char);
	return 0; 
};

	/// float
int MWStaticMPIRC::unpack( float *f,            int nitem, int stride )
{ 
	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Unpack(recvBuffer->data, recvBuffer->size, &recvPosition , &f, (nitem*stride), MPI_FLOAT, MPI_COMM_WORLD);
#ifdef S_PRINT

	MWprintf(10, "**** PID %d : Unpack float data = %f nitem = %d recvPosition = %d ****\n", id, *f, nitem, recvPosition);
#endif
	bytes_unpacked_ += nitem / stride * sizeof(float);
	return 0; 
};

		/// double
int MWStaticMPIRC::unpack( double *d, int nitem, int stride )
{ 
	int id, dim;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	dim = (sizeof d)*nitem;
	MPI_Unpack(recvBuffer->data, recvBuffer->size, &recvPosition , d, (nitem*stride), MPI_DOUBLE, MPI_COMM_WORLD);
#ifdef S_PRINT

	MWprintf(10, "**** PID %d : Unpack double recvPosition = %d ****\n", id, recvPosition);
#endif
	bytes_unpacked_ += nitem / stride * sizeof(double);
	return 0; 
};

		/// int
int MWStaticMPIRC::unpack( int *i, int nitem, int stride )
{ 
	int id, dim = (sizeof i)*nitem;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Unpack(recvBuffer->data, recvBuffer->size, &recvPosition ,i , (nitem*stride), MPI_INT, MPI_COMM_WORLD);
#ifdef S_PRINT

	MWprintf(10, "**** PID %d : Unpack int data = %i nitem = %d recvPosition = %d ****\n", id, *i, nitem, recvPosition);
#endif
	bytes_unpacked_ += nitem / stride * sizeof(int);
	return 0; 
};

		/// unsigned int
int MWStaticMPIRC::unpack( unsigned int *ui,    int nitem, int stride )
{ 
	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Unpack(recvBuffer->data, recvBuffer->size, &recvPosition, &ui, (nitem*stride), MPI_UNSIGNED, MPI_COMM_WORLD);
#ifdef S_PRINT

	MWprintf(10, "**** PID %d : Unpack ui data = %ui nitem = %d recvPosition = %d ****\n", id, *ui, nitem, recvPosition);
#endif
	bytes_unpacked_ += nitem / stride * sizeof(unsigned int);
	return 0; 
};

		/// short
int MWStaticMPIRC::unpack( short *sh,           int nitem, int stride )
{ 
	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Unpack(recvBuffer->data, recvBuffer->size, &recvPosition , &sh, (nitem*stride), MPI_SHORT, MPI_COMM_WORLD);
#ifdef S_PRINT

	MWprintf(10, "**** PID %d : Unpack sh data = %sh nitem = %d recvPosition = %d ****\n", id, *sh, nitem, recvPosition);
#endif
	bytes_unpacked_ += nitem / stride * sizeof(short);
	return 0; 
};

		/// unsigned short
int MWStaticMPIRC::unpack( unsigned short *ush, int nitem, int stride )
{ 
	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Unpack(recvBuffer->data, recvBuffer->size, &recvPosition, &ush , (nitem*stride), MPI_UNSIGNED_SHORT, MPI_COMM_WORLD);
#ifdef S_PRINT

	MWprintf(10, "**** PID %d : Unpack ush data = %ush nitem = %d recvPosition = %d ****\n", id, *ush, nitem, recvPosition);
#endif
	bytes_unpacked_ += nitem / stride * sizeof(unsigned short);
	return 0; 
};

		/// long
int MWStaticMPIRC::unpack( long *l,             int nitem, int stride )
{ 
	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Unpack(recvBuffer->data, recvBuffer->size, &recvPosition, &l , (nitem*stride), MPI_LONG, MPI_COMM_WORLD);
#ifdef S_PRINT 

	MWprintf(10, "**** PID %d : Unpack l data = %l nitem = %d recvPosition = %d ****\n", id, *l, nitem, recvPosition);
#endif
	bytes_unpacked_ += nitem / stride * sizeof(long);
	return 0; 
};

		/// unsigned long
int MWStaticMPIRC::unpack( unsigned long *ul,   int nitem, int stride )
{ 
	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Unpack(recvBuffer->data, recvBuffer->size, &recvPosition, &ul, (nitem*stride), MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#ifdef S_PRINT 
	MWprintf(10, "**** PID %d : Unpack ul data = %ul nitem = %d recvPosition = %d ****\n", id, *ul, nitem, recvPosition);
#endif
	bytes_unpacked_ += nitem / stride * sizeof(unsigned long);
	return 0; 
};

		/// Unpack a NULL-terminated string
int MWStaticMPIRC::unpack( char *str )
{ 
 	char* data, nitem;
	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	//Spacchetto prima la dimensione della stringa e poi la stringa stessa
	MPI_Unpack(recvBuffer->data, recvBuffer->size, &recvPosition , &nitem , 1 ,MPI_INT, MPI_COMM_WORLD);	
	MPI_Unpack(recvBuffer->data, recvBuffer->size, &recvPosition ,str , nitem ,MPI_CHAR, MPI_COMM_WORLD);
	str[nitem] = '\0'; 	
#ifdef S_PRINT 
	MWprintf(10, "**** PID %d : Unpack str data = %s nitem = %d recvPosition = %d ****\n", id, str, nitem, recvPosition);
#endif
	bytes_unpacked_ += strlen(str);
	return 0; 
};

int 
MWStaticMPIRC::write_RMstate ( FILE *fp = NULL )
{
    // No MPI specific functions
    return 0;
}
