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
#ifndef _blackbox_DRIVER_H
#define _blackbox_DRIVER_H

#include "MWDriver.h"
#include "MWTask_blackbox.h"

#include <set>

/** Blackbox Driver subclass derived from MWDriver */

class MWDriver_blackbox : public MWDriver 
{
public:
	MWDriver_blackbox();
	~MWDriver_blackbox();

	virtual MWReturn get_userinfo( int argc, char *argv[] ) = 0;
	virtual MWReturn setup_initial_tasks( int *, MWTask *** ) = 0;
	virtual MWReturn act_on_completed_task( MWTask * ) = 0;

		// Made final here -- you pass the executable -- worker can't pack initial data.
		//  Rather, they are allowed set set the "stage files"
	MWReturn pack_worker_init_data( void );
	virtual void printresults();
	virtual void write_master_state(FILE *fp) {}
	virtual void read_master_state(FILE *fp) {}

	void add_staged_file(const string &s) { stageFiles_.insert(s); }
	string get_executable() const { return theExecutable_; }
	void set_executable(const string &s) { theExecutable_ = s; }
	MWTask* gimme_a_task();

private: 
	string theExecutable_;
	std::set<string> stageFiles_;
	
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
