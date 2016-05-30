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
/* The main() of worker executable. Simply instantiate a Worker_blackbox 
 * class and go()!  */

#include "MW.h"
#include "MWWorker_blackbox.h"

int main(int argc, char *argv[]) 
{
		/* init a worker object */
	MWWorker_blackbox graduate_student;

		/* How much information you want the workers to print */
	set_MWprintf_level ( 75 );
	MWprintf ( 10, "A worker is starting.\n" );

		/* Go ! */ 
	graduate_student.go(argc, argv);
	return 0;
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
