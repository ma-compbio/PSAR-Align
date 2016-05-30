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
#ifndef MWRMCOMM_H
#define MWRMCOMM_H
#include <MWWorkerID.h>
#include "MWList.h"

/** This class is an abstract layer over a Resource Management (RM) and
	Communication (Comm) layer.  Hence the RMComm name.  It is designed
	to abstract away the differences between various systems.  

	The member functions in this class fall into two categories:  
	Resource Management and Message Passing (Communication).  The "master"
	in a master-worker application will use both; the worker will only
	use the communication facilities.
*/

enum MWRMCommErrors {
	WORKER_TIMEOUT = -50,
	UNABLE_TO_WAKE = -49,
	CANNOT_OPEN_INPUT_FILE = -48,
	SCANF_ERROR_ON_INPUT_FILE = -47,
	MESSAGE_SEQUENCE_ERROR = -46,
	UNKNOWN_DATA_TYPE = -45,
	WAITFILE_PROTOCOL_ERROR = -44,
	FAIL_MASTER_SEND = -43,
	CHECKSUM_ERROR_EXIT = -42,
	UNKNOWN_COMMAND = -41,
	INIT_REPLY_FAILURE = -40,
	UNPACK_FAILURE = -39
};

struct RMC_executable
{
	int arch_class;
	int exec_class;
	char *executable;
	char *executable_name;
	char *attributes;
};

class MWRMComm {

	friend class MWDriver;

public:

		/** Constructor.  Sets data to -1's. */
	MWRMComm();

		/// Destructor...
	virtual ~MWRMComm();

		/** @name A. Resource Management Functions

		Here are all the methods you could want to have for managing 
		a set of machines.  See each method for details...
		*/
		//@{

		/** System shutdown.  Does not return.  
		    @param exitval The value to call ::exit() with 
		*/
	virtual void exit( int exitval );

		/** Initialization of the master process.  This will be called
			one time only when the master starts up.  
			@return 0 on success, -1 on failure
		*/
	virtual int setup( int argc, char* argv[], int *my_id, int *master_id ) = 0;

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
	virtual int config( int *nhosts, int *narches, MWWorkerID ***w ) { return 0; };
	
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
	virtual int start_worker( MWWorkerID *w ) = 0;

		/** This routine is used to start up multiple workers at 
			the beginning of a run.  It should only be called one
			time.  It basically does a config() to find out what
			machines are available, and then starts a worker on 
			each of them.  You may want to check the implementations
			for details...
			@param nworkers The number of workers at the start (returned!)
			@param workers A pointer to a pointer to an array of
			MWWorkerID pointers.  The memory management for this 
			is the same as it is for the config() call - it should
			be null when called and will point to allocated memory
			when it comes back.
			@return 0 on success, -1 on failure
		*/
	virtual int init_beginning_workers( int *nworkers, MWWorkerID ***workers ) = 0;

		/** Remove a worker from the virtual machine.  This call 
			will delete w, so don't reference that memory again!
			@param w The MWWorkerID of the machine to remove.
			@return 0 on success, a negative number on failure
		*/
	virtual int removeWorker( MWWorkerID *w ) = 0;
		//@}

		/** @name B. Checkpointing Functions */
		//@{
		/** Write out internal state to a FILE *
		 */
	virtual int write_checkpoint( FILE * fp );

		/** Read in restart information from a FILE * 
		 */
	virtual int read_checkpoint( FILE * fp );

		/** Some Low level specific read/write functions */
	virtual int read_RMstate( FILE *fp ); 
		///
	virtual int write_RMstate ( FILE *fp );

		/** A function to restart the workers. This is called at restart from
			a ckpt file and it is meant to re-init all the workers.
		*/
	virtual int restart_beginning_workers ( int *num, MWWorkerID ***tempWorkers, MWmessages msg ) = 0;
		//@}

		/** @name C. Executable Management Interface functions */
		//@{

		/** Set the number of executable classes */
	void set_num_exec_classes ( int num );

		/** Returns the above number */
	int get_num_exec_classes ( );

		/** Set the number of arch classes */
	void set_num_arch_classes( int n );

		/** set the arch attribute */
	void set_arch_class_attributes ( int arch_class, const char *attr );

		/** Return the number of arch classes */
	int get_num_arch_classes();

		/** Set the number of executables */
	void set_num_executables ( int num );

	int get_num_executables();

		/** Set the name of the binary and the requirements string 
			for an arch class.  Technically, the requirements string 
			is not needed for MWPvmRC - it is defined in the submit
			file.  It *is* needed in the MWFileRC, however, for job
			submission.
			@param arch_class This specifies which arch class the above
			requirements will apply to.
			@param exec_name The name of the worker executable for 
			this arch class.
			@param requirements A string containing the "requirements" 
			attribute of a given arch class of workers.  This will 
			be something that can be used in a submit file, like 
			"(arch == \"INTEL\" && opsys == \"SOLARIS26\")"
		*/
	void add_executable ( int exec_class, int arch_class, char *exec_name, 
						  char *requirements ); 

		// The new, better way to do the above
	void add_executable(const char *exec_name, const char *requirements);

		/** If the RM software needs to "process" the executable name in
			some way, this is done here.
		*/
	virtual char* process_executable_name ( char* exec_name, int ex_cl, int num_arc );

		/** Set whether or not you would like worker checkpointing 
			(if the CommRM implementation has the capability)
		*/
	virtual void set_worker_checkpointing( bool wc );

protected:

		/** The number of exec classes */
	int exec_classes;

		/** The number of different arch classes. */
	int num_arches;

		/** The arch attributes */
	char **arch_class_attributes;

		/** The number of executables */
	int num_executables;
	int tempnum_executables;

		/** An array containing the {\bf full} pathnames of the executables.
			Element 0 is for arch "0", element 1 is for arch "1", etc.
			Usually read in get get_userinfo().
		*/
	struct RMC_executable **worker_executables;

		/// Would you like the workers to be checkpointed
	bool worker_checkpointing;

		/* moved here from MWDriver.C, to separate RMC and MW more */
	int *MW_exec_class_num_workers;

		//@}

public: 
	
		/** @name D. Host Management Members */
		//@{
		/** Set a "target" number of workers across all arches.  This 
			is useful if you don't care how many you get of each arch... 
			@param num_workers The target number of workers 
		*/
	void set_target_num_workers( int num_workers );

	void set_target_num_workers( int exec_class, int num_workers );
		
		// XXX Let the MWDriver/user know our current target number of workers
	int get_target_num_workers(int exec_class = -1);

		/** The RMComm will requests this many workers at a time */
	void set_worker_increment(int newinc);
	
	int get_worker_increment() const;

protected:

		/** This will figure out if we need to ask for more hosts
			or remove hosts.  It is called whenever a host is added
			or removed from the system, or set_target_num_workers()
			is called.
			@param num_workers A pointer to an array of length
			num_arches that contains the number of workers for 
			each arch class.
			@return If we have more workers than we need, we return a
			positive number as the "excess" that can be deleted.  
		*/
	virtual int hostaddlogic( int *num_workers ) = 0;

		/** The desired number of workers */
	int target_num_workers;

		/** The desired number of workers, exec_class wise */
	int *exec_class_target_num_workers;

		/** The current */
	int hostinc_;

		//@}
  
public:

		/** @name The Communication Routines

		These message passing routines are very closely modeled on
		PVM's message passing facility.  They are, however, pretty
		generic in that any message passing interface should be 
		able to implement them.
		*/
	
		//@{

		/** Initialize a buffer for sending some data.
			@param encoding Defined by each application.  0 = default */
	virtual int initsend ( int encoding = 0 ) = 0;
		/** Send the data that has been packed. 
			@param to_whom An identifier for the recipient
			@param msgtag A 'tag' to identify that type of message */
	virtual int send ( int to_whom, int msgtag ) = 0;

	//by jae
	virtual int pre_send ( int msgtag) { return -1; }
	virtual void setPollTimeOut(int secs) {}

		/** Receive some data that has been packed.  Should make this 
			more PVM-independent; will do this sometime.
			@param from_whom From a specific source; -1 is from all
			@param msgtag With a certain tag; -1 is all. */
	virtual int recv ( int from_whom, int msgtag ) = 0;

        /** non-blocking receive. This is the same as recv except it 
			returns the NO_MESSAGE tag instead of blocking if there is
			no data in the buffer */
	virtual int nrecv( int from_whom, int msgtag){return ABORT;};

		/** Provide info on the message just received */
	virtual int bufinfo ( int buf_id, int *len, int *tag, int *from ) = 0;

		/** For some system events like HOSTDELETE, TASKEXIT, etc, this will tell
			who was affected */
	virtual void who ( int *wh ) = 0;

		/** Needed only for MW-Independent */
	virtual void hostadd ( ) { };

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
	virtual int pack ( const char *bytes,         int nitem, int stride = 1 ) = 0;
		/// float
	virtual int pack ( const float *f,            int nitem, int stride = 1 ) = 0;
		/// double
	virtual int pack ( const double *d,           int nitem, int stride = 1 ) = 0;
		/// int
	virtual int pack ( const int *i,              int nitem, int stride = 1 ) = 0;
		/// unsigned int
	virtual int pack ( const unsigned int *ui,    int nitem, int stride = 1 ) = 0;
		/// short
	virtual int pack ( const short *sh,           int nitem, int stride = 1 ) = 0;
		/// unsigned short
	virtual int pack ( const unsigned short *ush, int nitem, int stride = 1 ) = 0;
		/// long
	virtual int pack ( const long *l,             int nitem, int stride = 1 ) = 0;
		/// unsigned long
	virtual int pack ( const unsigned long *ul,   int nitem, int stride = 1 ) = 0;
		/// Pack a NULL-terminated string
	virtual int pack ( const char *string ) = 0;

		/// by jae
	virtual void thread_start(int worker_id) {};
	virtual void thread_stop(int worker_id) {};
		/// Send raw data directly
	virtual int raw_pack( char* buf, int len ) { return 0; };

		/// Read raw data from Worker ID
	virtual int raw_unpack( int from_whom, char* buf, int len ) { return 0; };

		//@}

		/** @name Unpack Functions
			
		These unpack functions unpack data packed with the pack() 
		functions.  See the pack() functions for more details.

		*/
		//@{

		/// Unpack some bytes
	virtual int unpack ( char *bytes,         int nitem, int stride = 1 ) = 0;
		/// float
	virtual int unpack ( float *f,            int nitem, int stride = 1 ) = 0;
		/// double
	virtual int unpack ( double *d,           int nitem, int stride = 1 ) = 0;
		/// int
	virtual int unpack ( int *i,              int nitem, int stride = 1 ) = 0;
		/// unsigned int
	virtual int unpack ( unsigned int *ui,    int nitem, int stride = 1 ) = 0;
		/// short
	virtual int unpack ( short *sh,           int nitem, int stride = 1 ) = 0;
		/// unsigned short
	virtual int unpack ( unsigned short *ush, int nitem, int stride = 1 ) = 0;
		/// long
	virtual int unpack ( long *l,             int nitem, int stride = 1 ) = 0;
		/// unsigned long
	virtual int unpack ( unsigned long *ul,   int nitem, int stride = 1 ) = 0;
		/// Unpack a NULL-terminated string
	virtual int unpack ( char *string ) = 0;

		//@}
		//@}

		///
	double get_bytes_packed() const { return bytes_packed_; }

		///
	double get_bytes_unpacked() const { return bytes_unpacked_; }	
	
		/** When investigating the master-side bottleneck problem, we realized 
		 *  that sometimes the master can't handle all the messages sent from the 
		 *  many masters and the RMComm layer. The messages are accumulated in the 
		 *  RMComm layer's buffer (especially in the case of using MW-CondorPvm), 
		 *  but the master can't see them, not to mention process them or delete 
		 *  useless messages. This makes the master unable to accormodate the chanegs
		 *  in the environment, also makes most of the workers waiting for new tasks.

		 *  One way to deal with this situation is to pull all the messages in the 
		 *  underlying message passing layer and let the master buffer them by 
		 *  itself, to separate the message receiving and message processing steps. 
		 *
		 *  Thus we will introduce a FIFO list of buf_id (list of *int), each of the
		 *  buf_id points to a message buffer in PVM (and likewise in the other 
		 *  RMComms). The original message processing steps are: 
		 *  	RMC->recv()
		 *  	RMC->bufinfo()
		 *  	RMC->unpack(), RMC->unpack()
		 *  And now the process will be the same, but implementation is different:
		 *  	RMC->recv_all()	pull all messages out (in a non-blocking way) and 
		 *  			put each message's buf_id into the FIFO list.
		 *  	RMC->recv()	blocking recv	
		 *  	RMC->bufinfo()	examine the buf_info of the buffer with given buf_id
		 *  	RMC->setrbuf()	will be used to switch the active receive buffer
		 *  	RMC->freebuf()  will be used to free useless receive buffer
		 *  	RMC->unpack()
		 *  	
		 *  All recv will via the buffer_list now!
		 */
protected:	
		/* both the master and worker can have it! */
	MWList<void> *recv_buf_list;	

		///	Keeps the number of bytes packed
	double bytes_packed_;

		/// Holds the number of bytes unpacked
	double bytes_unpacked_;

public:
		/* returns the number of new buffers */
	virtual int recv_all( int from_whom, int msgtag ) = 0;	
		/* switch the active receive buffer, return the old active buf_id */
	virtual int setrbuf(int bid) = 0;		
		/* free the receive buffer. */
	virtual int freebuf(int bid) = 0;		
		/* advance the active buffer to the next valid buffer in the recv_buf_list. */
	virtual int next_buf() = 0;			
		/* return the recv_buf_list */
	virtual MWList<void>* recv_buffers() = 0;	
 private:
	void verify_file_exists(const char *);
};

#endif

/*
  Local Variables:
  mode: c++
  eval: (setq c-basic-offset 4)
  eval: (setq c-comment-only-line-offset 4)
  eval: (setq c-indent-level 4)
  eval: (setq c-brace-imaginary-offset 0)
  eval: (setq c-brace-offset 0)
  eval: (setq c-argdecl-indent 0)
  eval: (setq c-label-offset -4)
  eval: (setq c-continued-statement-offset 4)
  eval: (setq c-continued-brace-offset -4)
  eval: (setq c-tab-always-indent nil)
  eval: (setq tab-width 4)
  End:
*/
