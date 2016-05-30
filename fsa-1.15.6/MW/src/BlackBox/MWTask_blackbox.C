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
#include "MWTask_blackbox.h"
#include "MW.h"
#include <fstream>

/* init */
MWTask_blackbox::MWTask_blackbox() 
{
}


MWTask_blackbox::MWTask_blackbox(const list<string> &args,
								 const list<string> &input_files,
								 const list<string> &output_files)
{
	this->args_ = args;
	this->input_files_ = input_files;
	this->output_files_ = output_files;
}


/* destruction */
MWTask_blackbox::~MWTask_blackbox() {
}

/* print the task to stdout */
void
MWTask_blackbox::printself( int level ) 
{
	printList( level, "args: ", args_);
	printList( level, "inputs: ", input_files_);
	printList( level, "outputs: ", output_files_);
}

void
MWTask_blackbox::printList(int level, char *prefix, list<string>& l) {
	MWprintf(level, "%s", prefix);
	for (list<string>::iterator it = l.begin();	
		 it != l.end();
		 it++) {
		MWprintf(level, "%s ", (*it).c_str());
	}
	MWprintf (level, "\n");
}

/* The driver packs the input data via RMC, the data which will be sent to a worker. */
void 
MWTask_blackbox::pack_work( void ) 
{
	int size = args_.size();

	RMC->pack(&size, 1);
	MWprintf(30, "MWTask::pack_work packed arg size of %d\n", size);
	for (list<string>::iterator it = args_.begin();
		 it != args_.end();
		 it++) {
		RMC->pack((char *)(*it).c_str());
	}
	MWprintf(30, "MWTask::pack_work packing input files\n");
	send_files(input_files_);

	int num_output_files = output_files_.size();
	RMC->pack(&num_output_files, 1);
	for (list<string>::iterator it = output_files_.begin();
		 it != output_files_.end();
		 it++) {
		RMC->pack((char *)(*it).c_str());
	}
}

/* The worker unpacks input data via RMC, need to allocate space for data */
void 
MWTask_blackbox::unpack_work( void ) 
{
	args_.clear();
	int size;
	RMC->unpack(&size, 1);
	MWprintf(30,"MWTask::unpack_work unpacked arg size of %d\n", size);
	while (size--) {
		char str[8096];
		RMC->unpack(str);
		args_.push_back(str);
	}
	MWprintf(30, "MWTask::pack_work unpacking input files\n");
	recv_files(input_files_);

	int num_output_files = 0;
	output_files_.clear();
	RMC->unpack(&num_output_files, 1);
	while (num_output_files--) {
		char str[8096];
		RMC->unpack(str);
		output_files_.push_back(str);
	}
}

/* The worker packs result data via RMC, the result will be sent back to driver */
void 
MWTask_blackbox::pack_results( void ) 
{
	RMC->pack(&retVal_, 1);
	send_files(output_files_);
	remove("stdout");
	remove("stderr");
}

/* The driver unpacks result data via RMC */
void 
MWTask_blackbox::unpack_results( void ) 
{
	RMC->unpack(&retVal_, 1);
	recv_files(output_files_);
}

/* write checkpoint info per task, for each task haven't been finished */
void 
MWTask_blackbox::write_ckpt_info( FILE *fp ) 
{
  //list<string> args_;
  //list<string> input_files_;
  //list<string> output_files_;
  //int retVal_;

  fprintf(fp, "%d\n", retVal_);
  saveStringList(fp, args_);
  saveStringList(fp, input_files_);
  saveStringList(fp, output_files_);
}

/* Read checkpoint info, in the order written into the file */
void 
MWTask_blackbox::read_ckpt_info( FILE *fp ) 
{
  fscanf(fp, "%d", &retVal_);
  restoreStringList(fp, args_);
  restoreStringList(fp, input_files_);
  restoreStringList(fp, output_files_);
}


void 
MWTask_blackbox::send_files(list<string> &files) {
	int num_files = files.size();
	RMC->pack(&num_files, 1);
	for (list<string>::iterator i = files.begin(); i != files.end(); i++) {
		const char *file_name = i->c_str();

		fstream file(file_name);
		file.seekg(0, ios::end);
		size_t file_size = file.tellg();

			// If file doesn't exist, force file_size to 0
		if (file_size == (size_t) -1) {
			file_size = 0;
		}

		char *b = new char[file_size];
		file.seekg(0, ios::beg);
		file.read(b, file_size);

		RMC->pack((char *)file_name); 
		RMC->pack(&file_size, 1, 1);
		RMC->pack(b, file_size, 1);

		delete [] b;
		MWprintf(30, "packed file named %s of length %d\n", i->c_str(), file_size);
	}
}

void 
MWTask_blackbox::recv_files(list<string> &files) {
	int num_files;
	RMC->unpack(&num_files, 1);
	for ( int i = 0; i < num_files; i++) {
		char file_name[4096];
		int file_size;

		RMC->unpack(file_name);
		RMC->unpack(&file_size, 1);
		MWprintf(30, "unpacked file named %s of size %d\n", file_name, file_size);
		
		ofstream file(file_name);
		char *b = new char[file_size];
		RMC->unpack(b, file_size, 1);
		file.write(b, file_size);
		delete [] b;

		files.push_back(string(file_name));
	}
}

void
MWTask_blackbox::saveStringList(FILE *fp, list<string> &sl) {
	fprintf(fp, "%d\n", sl.size());
	for (list<string>::iterator i = sl.begin(); i != sl.end(); i++) {
		const char *name = i->c_str();
		fprintf(fp, "%s\n", name);
	}
}

void
MWTask_blackbox::restoreStringList(FILE *fp, list<string> &sl) {
	int size = 0;
	fscanf(fp, "%d", &size);
	while (size--) {
		char name[4096];
		fscanf(fp, "%s", name);
		sl.push_back(string(name));
	}
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
