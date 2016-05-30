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
#ifndef MWINDRC_H
#define MWINDRC_H
#include "../MWRMComm.h"

/** This class is an derived from the Resource Management (RM) and
	Communication (Comm) layer. In MW-Ind, there is just one worker and one master
	and both are on the same machine. And not just that both are same process!!
	Send and Recv are merely memcpy. This class I deem is useful for debugging purposes.
	The api etc all will remain the same when using this layer. All the user has to do
	is to compile it with the INDEPENDENT flag and use this as the RMComm layer. See into
	the fibonacci example for more details on the makefile. And then just launch the 
	application as a process. That's it.
*/

class MWIndRC : public MWRMComm 
{
    public:
		/** Constructor.  Sets data to -1's. */
	MWIndRC();
		/// Destructor...
	~MWIndRC();

	void hostadd ( );

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
	int setup( int argc, char* argv[], int *my_id, int *master_id);

		/** What is there to configure? Nothing.
		*/
	int config( int *nhosts, int *narches, MWWorkerID ***w );
	
		/** Just init the workerID to a worker. This is an idempotent function.
		*/
	int start_worker( MWWorkerID *w );

		/** This routine is used to start up multiple workers at 
			the beginning of a run.  In MW-Ind there is no initial worker.
			Simply return 0.
		*/
	int init_beginning_workers( int *nworkers, MWWorkerID ***workers );

		/** Called at the time of restart */
	int restart_beginning_workers ( int *nworkers, MWWorkerID ***workers, MWmessages msg );
			

		/** Remove a worker from the virtual machine.  MW-Ind has 
			no concept of host management. Hence all are null.
		 */
	int removeWorker( MWWorkerID *w );
		//@}

		/** @name B. Checkpointing Functions */
		//@{

		/** Do nothing */
	int read_RMstate( FILE *fp );

	int write_RMstate ( FILE *fp );
		//@}

 protected: 
	
		/** @name D. Host Management Members */
		//@{
		/** This will figure out if we need to ask for more hosts
			or remove hosts.  In MW-Ind version there is just one worker.
			So this function is NULL.
		*/
	int hostaddlogic( int *num_workers );

		//@}

#define INDP_BUF_SIZE 1024*1024*64
	
	struct SendBuf {
		int num_sent;
		char buf[INDP_BUF_SIZE];
	};
	struct SendBuf ind_sendbuf, ind_recvbuf;
	int last_tag;
	bool sent_data;

 public:

		/** @name The Communication Routines
			These communication primitives use memcpy as the communication
			mechanism for the transfer of data.

		*/
	
		//@{
		
		/**     A Function called to know which worker had an exception event */
	void who ( int *wh );

		/** Initialize a buffer for sending some data.
			@param encoding Defined by each application.  0 = default */
	int initsend ( int encoding = 0 );
		/** Send the data that has been packed. 
			@param to_whom An identifier for the recipient
			@param msgtag A 'tag' to identify that type of message */
	int send ( int to_whom, int msgtag );

		/** Receive some data that has been packed.  Just a memcpy
			@param from_whom From a specific source; -1 is from all
			@param msgtag With a certain tag; -1 is all. */
	int recv ( int from_whom, int msgtag );

	int nrecv ( int from_whom, int msgtag );	

		/** Provide info on the message just received */
	int bufinfo ( int buf_id, int *len, int *tag, int *from );

		/** @name Pack Functions
			
			In the following pack() functions, there are some common themes.
			First, each stuffs some data into a buffer to be sent.  The
			nitem parameter is just a count of the number of items.  The 
			stride parameter specifies *which* items to pack.  1 implies
			all, 2 would be every 2nd item, 3 is every 3rd item, etc. 

			The return value is user defined.  It should be standardized, 
			but I'll do that later.
		*/
		//@{

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
	int pack ( const char *string );

		//@}

		/** @name Unpack Functions
			
			These unpack functions unpack data packed with the pack() 
			functions.  See the pack() functions for more details.

		*/
		//@{

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
	int unpack ( char *string );
		//@}
		//@}

	/** The non-blocking recv() functions */
	public:
		int	recv_all( int from_whom, int msgtag );	// returns the number of new buffers
		int 	setrbuf(int bid);		// switch the active receive buffer
							// return the old active buf_id
		int	freebuf(int bid);		// free the receive buffer.
		MWList<void> * recv_buffers();			// return the recv_buf_list
		int 	next_buf();			// advance the active buffer to the next valid buffer 
							// in the recv_buf_list.

};
#endif

