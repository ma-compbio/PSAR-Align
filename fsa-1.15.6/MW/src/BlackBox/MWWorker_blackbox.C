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
#include "MWWorker_blackbox.h"
#include "MWTask_blackbox.h"

#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <fcntl.h>

/* init */
MWWorker_blackbox::MWWorker_blackbox() 
{
  workingTask = new MWTask_blackbox;
}

/* destruct */
MWWorker_blackbox::~MWWorker_blackbox() 
{
    delete workingTask;
}

/* Do benchmark and return result (usually the time to task t), t is supposed
 * to be a benchmark task.  In this app, it just send back a PI. */
double
MWWorker_blackbox::benchmark( MWTask *t ) 
{
        MWTask_blackbox *bb = dynamic_cast<MWTask_blackbox *> ( t );
		bb->printself(30);
        return 1.0;
}

/* unpack the init data from the driver */
MWReturn 
MWWorker_blackbox::unpack_init_data( void ) 
{
		// Get the executable
	list<string> l;
	MWTask_blackbox::recv_files(l);
	MWprintf(10, "Worker_blackbox got %d files\n", l.size());
	if (l.size() < 1) {
	  MWprintf(1, "Warning!?!  Where is your executable!?!?!\n");
	}
	else {
	   list<string>::const_iterator i = l.begin();
	   const char *c = i->c_str();	
	   executable_ = *i;
	   MWprintf(20, "marking worker executable %s as 0755\n", c);
	   ::chmod(c, 0755);  // and make it executable
	}
	return OK;
}

/* Execute each task */
void 
MWWorker_blackbox::execute_task( MWTask *t ) 
{
	int retVal;
	
	MWprintf(30, "Enter Worker_blackbox::execute_task\n");
	MWTask_blackbox *bb = dynamic_cast<MWTask_blackbox *> ( t );
	
	MWprintf(40, "     task is \n");
	bb->printself(40);

	const list<string>& args = bb->getArgs();
	int argc = args.size() + 2;

	char **argv = (char **)malloc( sizeof(char *) * argc);
	
	int index = 0;
	argv[index] = (char *) executable_.c_str();
	index++;

	for (list<string>::const_iterator i = args.begin(); i != args.end(); i++) {
		argv[index] = (char *) (*i).c_str();
		MWprintf(50, "execute_task:: argv[%d] is %s\n", index, argv[index]);
		index++;
	}
	argv[index] = 0;

	int pid = fork();
	if (pid == 0) {
			// Child

			// make stdout to to the file named "stdout"
		close(1);
		close(2);
		
		open("stdout", O_CREAT | O_TRUNC);
		open("stderr", O_CREAT | O_TRUNC);

		execve(argv[0], argv, 0);
		MWprintf(10, "execute_task failed to exec %s errno = %d (%s)\n",	
				 argv[0], errno, strerror(errno));

		exit(-1);
	} 

	::waitpid(pid, &retVal, 0);
	bb->setReturnVal(retVal);
	MWprintf(30, "Leave Worker_blackbox::execute_task retVal was %d\n", retVal);
}

MWTask* 
MWWorker_blackbox::gimme_a_task()
{
	return new MWTask_blackbox;
}

/* Just return a newly created application worker object */
MWWorker*
gimme_a_worker ()
{
       	return new MWWorker_blackbox;
}
