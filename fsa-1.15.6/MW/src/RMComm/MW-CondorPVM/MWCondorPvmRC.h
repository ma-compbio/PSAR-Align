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
#ifndef MWPVMRC_H
#define MWPVMRC_H

#include "../MWRMComm.h"
#include "MWList.h"
extern "C" {
#include "pvm3.h"
}

/// A mask for getting rid of that annoying high set bit.
#define PVM_MASK 0x7fffffff

struct condorpvm_submit_element
{
	int arch_class;
	int exec_class;
	char *exec;
};

/** 
	A Resource Management and Communication class that uses PVM 
	for underlying support of inter-process communication and 
	resource management.

	We *could* do cool things here...like keep track of amount of
	data sent, to whom, etc.  We could also MWprintf ( 99, "..." )
	everything that gets sent for heaps and heaps of debugging.

*/

class MWPvmRC : public MWRMComm 
{

	public:

		/// Constructor...
		MWPvmRC();

		/// Destructor...
		~MWPvmRC();

		/** @name A. Resource Management Routines

			Here we implement the pure virtual functions found in 
			ur parent class, MWRMComm.
		*/
		//@{

		/** Shut down.  Calls pvm_exit(), then ::exit(exitval). */
		void exit( int exitval );

		/** Initialization.  Does pvm_catchout(), pvm_parent(), 
			and pvm_notifies(). */
		int setup( int argc, char *argv[], int *my_id, int *master_id );

		/** Does a pvm_config, and stuffs the returned hostinfo
		    struct information into a lot of MWWorkerID's.  See
		    the parent class for return details.  */
		int config( int *nhosts, int *narches, MWWorkerID ***workers );

		/** Basically do a pvm_spawn.  We also put some information
			into w upon return.
			@return The tid spawned.  Negative number on error.
		*/
		int start_worker ( MWWorkerID *w );
	
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
		int restart_beginning_workers ( int *nworkers, MWWorkerID ***workers, MWmessages msg );

		/** Remove a worker.  Basically, call pvm_delhosts(). */
		int removeWorker( MWWorkerID *w );

		/** Read some state. A null function */
		int read_RMstate ( FILE *fp );

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
		int hostaddlogic( int *num_workers );

		/** Do the pvm_spawn and associated stuff */
		int do_spawn( MWWorkerID *w );

		/** Set up the proper notifies for a task & host tid */
		int setup_notifies ( int task_tid );

		/** This function says to pvm: "I would like another machine, 
			please".
			@param howmany The number of machines to request
			@param archnum The arch class number to ask for. */
		int ask_for_host( int howmany, int archnum );

		/** Shows the pvm virtual machine according to pvm.  Used
			for debugging purposes; not normally called. */
		void conf();
	
		/** Helper for hostaddlogic().  Returns the index of the min
			element in array, which is of length len */
		int min ( int *array, int len );

		/** The number of outstanding requests for a particular 
			arch. */
		//int *hostadd_reqs;

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
		int initsend ( int encoding = PvmDataDefault);

		///
		int send ( int to_whom, int msgtag );

		///
		int recv ( int from_whom, int msgtag );

        /// non-blocking version of receive
        int nrecv (int fromWhom, int msgtag);

		/** Provide info on the message just received */
		int bufinfo ( int buf_id, int *len, int *tag, int *from );

		/** Tells the affected party */
		void who ( int *id );

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

		/** Find a suitable executable_class for a particular arch_class.
		This function is useful to pick up the executable */
		struct condorpvm_submit_element* find_exec_class ( int arch );

		/** In the beginning, find a exec_class to an already
		acquired arch_class */
		int choose_exec_class ( int *tempi );

		/** Find the appropriate executable name that belongs to exec_class ex_cl and
		arch_class ar_cl */
		char* find_executable ( int ex_cl, int ar_cl );

		/** Sort the array temp based on the ratio of current/desired
		number of workers of each exec_class */
		void sort_exec_class_ratio ( int *temp );

		/** Find if there is an executable that belongs to exec_class ex_cl and
		arch_class ar_cl and if it exists return the executable name */
		char* exists_executable ( int ex_cl, int ar_cl );

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

#endif
