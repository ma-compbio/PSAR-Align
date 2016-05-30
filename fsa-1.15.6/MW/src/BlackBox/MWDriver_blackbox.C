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
/* These methods will be implemented to reflect the application behavior */

#include "MW.h"
#include "MWDriver_blackbox.h"
#include "MWWorker_blackbox.h"
#include <unistd.h>

/* Statstics for how many bytes are packed (sent) and unpacked (received) */
extern double RMCOMM_bytes_packed;
extern double RMCOMM_bytes_unpacked;

/* initialization */
MWDriver_blackbox::MWDriver_blackbox() 
{
}

/* destruction */
MWDriver_blackbox::~MWDriver_blackbox() 
{
}

/* */
MWReturn 
MWDriver_blackbox::pack_worker_init_data( void ) 
{
	string exe = get_executable();
	list<string> l;
	l.push_back(exe);

	for (set<string>::const_iterator it = stageFiles_.begin(); 
		 it != stageFiles_.end(); ++it) {
		l.push_back(*it);
	}

	MWTask_blackbox::send_files(l);

	

	return OK;
}

/* Print out the result when MW is done. MW assume that the application 
 * is keeping track of the results :-) */
void 
MWDriver_blackbox::printresults() 
{
	MWprintf ( 10, "BlackBox Driver Complete\n");
}

/* Return a new application task object */
MWTask*
MWDriver_blackbox::gimme_a_task() 
{
	return new MWTask_blackbox;
}
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
  eval: (setq tab-width 4)
  End:
*/
