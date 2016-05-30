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
#ifndef _blackbox_TASK_H
#define _blackbox_TASK_H

#include "MWTask.h"

using namespace std;

#include <string>
#include <list>

class MWTask_blackbox : public MWTask 
{
public:
	/* constructors */
  MWTask_blackbox();
  MWTask_blackbox(const list<string> &args,
		  const list<string> &input_files,
		  const list<string> &output_files);
  
  /* destructor */
  virtual ~MWTask_blackbox();

	const std::list<string> & getArgs() const { return args_; }
	const std::list<string> & getInputFiles() const { return input_files_; }
	const std::list<string> & getOutputFiles() const { return output_files_; }

	int getReturnVal() const { return retVal_; }
	void setReturnVal(int rv){ retVal_ = rv; }
  
  /* App is required to implement the following functions. */
  void pack_work( void );
  void unpack_work( void );
  void pack_results( void );
  void unpack_results( void );
  
  /* The following functions have default implementation. */
  void printself( int level = 70 );
  void printList(int level, char *prefix, list<string>& l);
  void write_ckpt_info( FILE *fp );
  void read_ckpt_info( FILE *fp );

  /* The application specific information goes here */
  static void send_files(list<string> &files);
  static void recv_files(list<string> &files);
	
  static void saveStringList(FILE *fp, list<string> &sl);
  static void restoreStringList(FILE *fp, list<string> &sl);

protected:  
  list<string> args_;
  list<string> input_files_;
  list<string> output_files_;
  int retVal_;
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
