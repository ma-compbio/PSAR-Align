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

#ifndef MWSOCKETRC_H
#define MWSOCKETRC_H

#include "../MWRMComm.h"
#include <stdio.h>

#ifdef WINDOWS
#define _POSIX_PATH_MAX 1024
#endif

typedef enum { MWSocket_FREE, MWSocket_SUBMITTED, MWSocket_EXECUTING, MWSocket_THREAD } SocketWorkerState;
#define MWSOCKET_MAX_MSG_SIZE (16 * 1024 * 1024)

#define MWSocketHostAdd 0
#define MWSocketHostDelete 1
#define MWSocketTaskExit 2
#define MWSocketTaskSuspend 3
#define MWSocketTaskResume 4
#define MWSocketChecksumError 5

static int MWSocketMessageTags[] = 
		{	HOSTADD,
			HOSTDELETE,
			HOSTDELETE,
			HOSTSUSPEND,
			HOSTRESUME,
			CHECKSUM_ERROR
		};

typedef struct
{
	SocketWorkerState state;
	int arch;
	int cId;
	int subcId;
	int socket;
	int exec_class;
}SocketWorker;

/** 
	A Resource Management and Communication class that uses Condor 
	for underlying support of resource managament. Some crude inter-process 
	communication is provided using the userlog feature of Condor.
	resource management.

*/

class MWSocketRC : public MWRMComm
{

  public:

	///  Constructor 
	MWSocketRC( bool maste, int id );

	///  Destructor 
	~MWSocketRC();

 /** @name A. Resource Management Routines

	 Here we implement the pure virtual functions found in 
	 ur parent class, MWRMComm.
 */
 //@{

	/**  Initialises. Depending on whether it is master or worker instance
	     it initializes all the internal variables. */
	int setup ( int argc, char *argv[], int *mytid, int *mastertid );

	/**  Shutdown. Kills all the workers if it is master */
	void exit ( int exitval );

	/**  Initialize workers if already some have started up */
	int init_beginning_workers ( int *nworkers, MWWorkerID ***workers );

	/** Restart function */
	int restart_beginning_workers ( int *nworkers, MWWorkerID ***workers, MWmessages msg );

	/**  This function is actually a misonomer. It DOES NOT spawn a new
	     worker. Rather it just inits the structure that is passed on to 
	     it */
	int start_worker ( MWWorkerID *w );

	/**  This function removes a existing worker */
	int removeWorker ( MWWorkerID *w );

	/**  Figure out whether or not to generate a new worker depending on
	     whether new requests have been made */
	int hostaddlogic ( int *w );

	/** A routine for reading in the MW-File state at the time 
	 *  of checkpointing.
	 */
	int read_RMstate ( FILE *fp = NULL );

	/** A routine for writing in the MW-File state at the time 
	 *  of checkpointing.
	 */
	int write_RMstate ( FILE *fp = NULL );

	char* process_executable_name ( char *execut, int, int );

		//@}

	/** @name B. Communication Routines

		Unlike MWPvmRC, the communication routines are non-trivial 
		because Condor provides no inter-process comminucation.
		Thus we use files for communication. So a send is essentially
		a file write operation and a recv is a file-read operation.
		We maintain 2 lists:- The sendList and recvList for taking 
		care of what is to be written/read to/from the files. As in
		pvm a user beings by calling initsend which creates a new list.
		Calls to pack insert into the list what is being packed. And
		finally a send writes the entire thing into a file 
		identified by the destination. Corresponding things happpen
		in recv.
	*/
	//@{

	/** A Function called to know which worker had an exception event */
	void who ( int *wh );

	/**  Initialize the send buffer */
	int initsend ( int useless = 0 );

	/**  Send function */
	int send ( int toWhom, int msgtag );

	// by jae
	/** Send function for Jae */
	int pre_send ( int msgtag );
	void setPollTimeOut(int secs) { m_poll_timeout = secs; }

	/**  Recv function */
	int  recv ( int fromWhom, int msgtag );

	/**  Non-blocking version of recv */
	int nrecv (int fromWhom, int msgtag);

	/**  Get some info about the recv buffer */
	int bufinfo ( int buf_id, int *len, int *tag, int *sending_host );

	/// pack some bytes
	int pack ( const char *bytes, int nitem, int stride = 1 );

	/// float
	int pack ( const float *f, int nitem, int stride = 1 );

	/// double
	int pack ( const double *d, int nitem, int stride = 1 );

	/// int
	int pack ( const int *i, int nitem, int stride = 1 );

	/// unsigned int
	int pack ( const unsigned int *ui, int nitem, int stride = 1 );

	/// short
	int pack ( const short *sh, int nitem, int stride = 1 );

	/// unsigned short
	int pack ( const unsigned short *ush, int nitem, int stride = 1 );

	/// long
	int pack ( const long *l, int nitem, int stride = 1 );

	/// unsigned long
	int pack ( const unsigned long *ul, int nitem, int stride = 1 );

	/// string
	int pack ( const char  *str );

	/// by jae
	/***************************/
	void thread_start(int worker_id);
	void thread_stop(int worker_id);

	/// Send raw data directly
	int raw_pack( char* buf, int len );

	/// Read raw data from Worker ID
	int raw_unpack( int from_whom, char* buf, int len );
	/***************************/

	/// Unpack some bytes
	int unpack ( char *bytes, int nitem, int stride = 1 );

	/// float
	int unpack ( float *f, int nitem, int stride = 1 );

	///double
	int unpack ( double *d, int nitem, int stride = 1 );

	/// int
	int unpack ( int *i, int nitem, int stride = 1 );

	/// unsigned int
	int unpack ( unsigned int *ui, int nitem, int stride = 1 );

	/// short
	int unpack ( short *sh, int nitem, int stride = 1 );

	/// unsigned short
	int unpack ( unsigned short *ush, int nitem, int stride = 1 );

	/// long
	int unpack ( long *l, int nitem, int stride = 1 );

	/// unsigned long
	int unpack ( unsigned long *ul, int nitem, int stride = 1 );

	/// string
	int unpack ( char  *str );

		//@}
  private:

	/**	Some private functions */
	/** Some marshalling functions */
	/// int
	inline void marshall_int ( int n, char *buf );

	/// unsigned int
	inline void marshall_uint ( unsigned int n, char *buf );

	/// short
	inline void marshall_short ( short n, char *buf );

	/// unsigned short
	inline void marshall_ushort ( unsigned short n, char *buf );

	/// long
	inline void marshall_long ( long n, char *buf );

	/// unsigned long
	inline void marshall_ulong ( unsigned long n, char *buf );

	/// float
	inline void marshall_float ( float n, char *buf );

	/// double
	inline void marshall_double ( double n, char *buf );

	/** Some unmarhsalling routines */
	/// int
	inline int unmarshall_int ( char *buf );

	/// unsigned int
	inline unsigned int unmarshall_uint ( char *buf );

	/// short
	inline short unmarshall_short ( char *buf );

	/// unsigned short
	inline unsigned short unmarshall_ushort ( char *buf );

	/// long
	inline long unmarshall_long ( char *buf );

	/// unsigned long
	inline unsigned long unmarshall_ulong ( char *buf );

	/// float
	inline float unmarshall_float ( char *buf );

	/// double
	inline double unmarshall_double ( char *buf );

	/**	Handle a message from the worker */
	int handle_finished_worker ( int i );

	/**	Is Called when a task is dead */
	int handle_killed_worker ( int i );

	/** 	Called when the worker starts executing first */
	int handle_starting_worker( );

	/**	Create a worker */
	int do_spawn ( int numworkers, int arch );

	/** Kill all the workers */
	int killWorkers ( );

	/**	Kill a worker */
	int killWorker ( int i );

	/** Final closing of sockets */
	void closeSockets ( );

	/** A primitive Send function */
	int finalSend ( int s, char *buf, size_t len, int flags );

	/** A primitive Recv function */
	int finalRecv ( int s, char *buf, size_t len, int flags );

	int nfinalRecv ( int s, void *buf, size_t len, int flags );

	/**	Master receive */
	int masterrecv ( int fromWhom, int msgtag );

	/**	Worker receive */
	int workerrecv ( int fromWhom, int msgtag );

	/** Master recive from anyone */
	int masterrecvany();

	/** A generic function to create and bind a socket to a address */
	int creatAndBind ( int port );

	/** The function that checks the sockets */
	int checkSockets ( );

	/** A function to sort the exec_classes */
	void sort_exec_class_ratio ( int *temp );

	/**	The bool indicates the mode of the RC instance 
	 *  a true value means that it is a master and a false means that 
	 * it is a worker.
	 *  We need to make it a tristate as somtimes an RC can be both
	 *  a master and a worker.
	 */
	bool isMaster;

	/**	The value is the id that the worker gets. For master it is
	 * of no use.
	 */
	int socketId;

	/** The id of the master */
	int masterId;

	/** The address of the master */
	char masterAddress[_POSIX_PATH_MAX];
	/** The value of the master Socket */
	int masterSocket;

	/** The actual Socket */
	int socketComm;

	char control_directory[_POSIX_PATH_MAX];

	char worker_number_file[_POSIX_PATH_MAX];

	/**  	The number of workers that were present. This is an internal
	 *   variable that will be used for copying the entire thing when 
	 *   the target number changes. This is the number that will be maintained close to 
	 *   target_num_workers.
	 */
	int current_num_workers;

	/**  The number of submitted workers */
	int submitted_num_workers;

	/**  An array of the number of workers */
	SocketWorker *socketWorkers;

	/** A reverse map in the master from socket to socketWorkers */
	int reverseMap[4096];

	/**  This keeps the track of how many cycles were made.
	 */
	int turnNo;

	/**  The tag of the message that just came in	*/
	int msgTag;

	/**  The message came from whom 	*/
	int whomRecv;

	/**  A variable to keep track of submit files */
	int subId;

	/**  A variable array keeping track of how many have been requested */
	int *hostadd_reqs;

	/**  A variable array keeping track of how many workers of each arch */
	int *num_workers;

	/**  The worker timeout in minutes */
	long worker_timeout;

	/** For the sake of transmission */
	char sendBuffer[MWSOCKET_MAX_MSG_SIZE];
	/** The current message size */
	int messageSize;

	/** For the sake of reception */
	char recvBuffer[MWSOCKET_MAX_MSG_SIZE];
	/** The received message size */
	int recvSize;
	int recvPointer;
	fd_set rfds;

	//by jae
	bool	pre_send_done;
	int 	m_poll_timeout;

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

