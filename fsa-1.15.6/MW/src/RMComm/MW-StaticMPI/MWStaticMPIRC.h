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
#ifndef MWSTATICMPIRC_H
#define MWSTATICMPIRC_H

#define MPICH_SKIP_MPICXX
#define KB64 1048576

#include "MWRMComm.h"
#include "MWList.h"

extern "C" {
#include "mpi.h"
}
//#include "/usr/local/mpich-g2/include/mpi.h"

class Msg_buf  
{
public:
	Msg_buf(int max_size = KB64)
	{ 
          data = new char[max_size]; 
          size = 0;
          pos = 0;
	  tag = 0;
        };

	
	/// Destructor
	virtual ~Msg_buf( ) { delete[] data; };

	/// Set everything to default
	void reset(int max_size = KB64)
	{
	  tag = 0;
	  max_size = max_size;
	  size = 0;
	  pos = 0;
	  id = -1;
	  delete[] data;
	  data = new char[max_size];
	}
	//@}

	/// Resize the buffer if it is not ehough place to add <code> add_size </code>.
	void resize(int add_size)
	{
	  if( size + add_size <= max_size)
	    return;
	  max_size = size + add_size + KB64;
  	  char* tmp = new char[max_size];
	  memcpy(tmp, data, size);
          delete[] data;
	  data = tmp;
	}
	
	/// Message tag.
	int tag;

	/// Pointer to the data
	char* data;


	/// Allocated size.
	int max_size;

	/// Current message size
	int size;

	/// Current position in the buffer to read from the buffer
	int pos;
	
	/// ID of the buffer
	int id;
	
	/// Sender pid.
	int sender;
};




class MWStaticMPIRC: public MWRMComm 
{

	

public:

		/** Constructor.  Sets data to -1's. */
	MWStaticMPIRC();

		/// Destructor...
	~MWStaticMPIRC();

		/** @name A. Resource Management Functions

		Here are all the methods you could want to have for managing 
		a set of machines.  See each method for details...
		*/
		//@{

		/** System shutdown.  Does not return.  
		    @param exitval The value to call ::exit() with 
		*/
		 void exit( int exitval );

		/** Initialization of the master process.  This will be called
			one time only when the master starts up.  
			@return 0 on success, -1 on failure
		*/
		int setup( int argc, char* argv[], int *my_id, int *master_id ) ;

		void init_env(int argc, char* argv[]);
		bool is_master_proc();

		/** This returns the state of the "virtual machine" - aka the 
			set of workers under the master's control.  
			@param nhosts The number of hosts
			@param narches The number of architechture classes
			@param w A pointer to a pointer to an array of pointers.  
			This should be NULL when called; config() will allocate
			the memory using new and return it to you.  Don't forget
			to delete not only the elements in that array, but 
			also the array itself.  The array will have nhosts
			elements, and they will each represent a worker 
			machine.
			@return 0 on success, -1 on error.
		*/
		int config( int *nhosts, int *narches, MWWorkerID ***w ) ;
	
		/** Start a worker on a machine that has been given to you.
			This is really only important if the process of starting
			a worker is two-stage.  For instance, in pvm, you
			first get a machine.  THEN you have to spawn an 
			executable on that machine.  For something like MW-files, 
			once you get a worker it already has the executable 
			started on it.  In that case this will basically be
			a no-op.  
			@param w A pointer to a MWWorkerID.  This must point
			to allocated memory.  It is filled with info like the
			machine name, id2, and the arch.
			@return id2 on success, -1 on failure
		*/
		int start_worker( MWWorkerID *w ) ;

		
	
		/** Start up the workers that are given to us at the beginning. 
			See base class comments for more details.
			@param nworkers The number of workers at start
			@param workers A pointer to an array of pointers to 
			       MWWorkerID classes.  This call will new() memory 
				   for the MWWorkerIDs.  Also, if (*w)[n]->id2 is
				   -1, that means that the spawn failed for 
				   worker number n.
		*/
		int init_beginning_workers ( int *nworkers, MWWorkerID ***workers );

		/** Called at the time of restart */
		int restart_beginning_workers ( int *nworkers, MWWorkerID ***worker, MWmessages msg );

		/** Remove a worker.  Basically, call pvm_delhosts(). */
		int removeWorker( MWWorkerID *w ){ return 0; };

		/** Read some state. A null function */
		int read_RMstate ( FILE *fp ){ return 0; };

		/** Write some state. A null function */
		int write_RMstate ( FILE *fp );
		
		
		//@}

	private:
		/** @name Other Functions to be implemented */

		//@{

		/** Figure out wether or not to ask for hosts, and how many... 
			@param num_workers An array of size num_arches that contains
			the number of workers in each arch class.
		 */
		int hostaddlogic( int *num_workers ){ return 0; };

		/** Do the pvm_spawn and associated stuff */
		int do_spawn( MWWorkerID *w ){ return 0; };

		/** Set up the proper notifies for a task & host tid */
		int setup_notifies ( int task_tid ){ return 0; };

		/** This function says to pvm: "I would like another machine, 
			please".
			@param howmany The number of machines to request
			@param archnum The arch class number to ask for. */
		int ask_for_host( int howmany, int archnum ){ return 0; };

		/** Shows the pvm virtual machine according to pvm.  Used
			for debugging purposes{ return 0; }; not normally called. */
		void conf(){ };
	
		/** Helper for hostaddlogic().  Returns the index of the min
			element in array, which is of length len */
		int min ( int *array, int len ){ return 0; };

		/** Funzioni di ricezione del worker e del master */

 		int iworkerrecv(int from_whom, int msgtag);
		int masterrecv(int from_whom, int msgtag);
		
		/** The number of outstanding requests for a particular 
			arch. */
		//int *hostadd_reqs{ return 0; };

		/** The requests per exec_classes */
		int **hostadd_reqs;

		MWList<void> * submit_list;

		//@}

		/** Am I the master process or not? */
		int is_master;

	public:


		/** @name C. Communication Routines
			
			These are essentially thin wrappers of PVM calls.
		*/
		//@{

		///
		int initsend ( int encoding = 0);

		///
		int send ( int to_whom, int msgtag );

		///
		int recv ( int from_whom, int msgtag );

        	/// non-blocking version of receive
       		int nrecv (int fromWhom, int msgtag);

		/** Provide info on the message just received */
		int bufinfo ( int buf_id, int *len, int *tag, int *from );

		/** Tells the affected party */
		void who ( int *id ){  };

		/// Pack some bytes
		int pack ( const char *bytes,         int nitem, int stride = 1 );

		/// float
		int pack ( const float *f,            int nitem, int stride = 1 );

		/// double
		int pack ( const double *d,           int nitem, int stride = 1 );

		/// int
		int pack ( const int *i,              int nitem, int stride = 1 );

		/// unsigned int
		int pack ( const unsigned int *ui,    int nitem, int stride = 1 );

		/// short
		int pack ( const short *sh,           int nitem, int stride = 1 );

		/// unsigned short
		int pack ( const unsigned short *ush, int nitem, int stride = 1 );

		/// long
		int pack ( const long *l,             int nitem, int stride = 1 );

		/// unsigned long
		int pack ( const unsigned long *ul,   int nitem, int stride = 1 );

		/// Pack a NULL-terminated string
		int pack ( const char *str );
	
		/// Unpack some bytes
		int unpack ( char *bytes,         int nitem, int stride = 1 );

		/// float
		int unpack ( float *f,            int nitem, int stride = 1 );

		/// double
		int unpack ( double *d,           int nitem, int stride = 1 );

		/// int
		int unpack ( int *i,              int nitem, int stride = 1 );

		/// unsigned int
		int unpack ( unsigned int *ui,    int nitem, int stride = 1 );

		/// short
		int unpack ( short *sh,           int nitem, int stride = 1 );

		/// unsigned short
		int unpack ( unsigned short *ush, int nitem, int stride = 1 );

		/// long
		int unpack ( long *l,             int nitem, int stride = 1 );

		/// unsigned long
		int unpack ( unsigned long *ul,   int nitem, int stride = 1 );

		/// Unpack a NULL-terminated string
		int unpack ( char *str );
		
		
		//@}

	private:
		/** @name C. Some Misc Helper functions
		*/

		//@{
	    
        /** Number of sent messages */  	
            int num_send_msgs;

	/** Number of received messages */
            int num_recv_msgs;
	

 	/** The current message to send */
            Msg_buf *sendBuffer;

        /** For the sake of transmission */
	    Msg_buf *recvBuffer;
	  	
	    int sendPosition;
	    int recvPosition;
		//@}
		/** The non-blocking recv functions */
	public:
		int	recv_all( int from_whom, int msgtag );	// returns the number of new buffers
		int 	setrbuf(int bid);		// switch the active receive buffer
							// return the old active buf_id
		int	freebuf(int bid);		// free the receive buffer.
		MWList<void>* recv_buffers();			// return the recv_buf_list
		int 	next_buf();			// advance the active buffer to the next valid buffer 
							// in the recv_buf_list.
		int 	switch_buf_back();		// Switch the active receive buffer to the buffer
							// on the head of the recv_buf_list;
};

extern MWStaticMPIRC *MWGlobalMPIRMComm;
#endif
