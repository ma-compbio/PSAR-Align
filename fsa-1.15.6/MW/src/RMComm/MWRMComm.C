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
#include "MWRMComm.h"

#ifdef WINDOWS
#define _POSIX_PATH_MAX 1024
#endif

MWRMComm::MWRMComm()
{
	exec_classes = 0;
	num_arches = 0;
	arch_class_attributes = NULL;
	num_executables = tempnum_executables = 0;
	hostinc_ = 6;

	worker_executables = NULL;
	worker_checkpointing = false;	

	target_num_workers = 0;
	exec_class_target_num_workers = NULL;
	MW_exec_class_num_workers = NULL;
}


MWRMComm::~MWRMComm()
{
	int i;

	if ( arch_class_attributes )
		{
			for ( i = 0; i < num_arches; i++ )
				if ( arch_class_attributes[i] )
					delete [] arch_class_attributes[i];
		}

	if ( worker_executables )
		{
			for ( i = 0; i < num_executables; i++ )
				if ( worker_executables[i] )
					{
						delete [] worker_executables[i]->executable;
						delete [] worker_executables[i]->executable_name;
						delete [] worker_executables[i]->attributes;
						delete worker_executables[i];
					}
		}

	if ( exec_class_target_num_workers )
		delete [] exec_class_target_num_workers;

	if ( MW_exec_class_num_workers )
		delete [] MW_exec_class_num_workers;
}

void MWRMComm::exit( int exitval )
{
	::exit ( exitval );
}

int MWRMComm::write_checkpoint( FILE *fp )
{
	int i;

	MWprintf ( 80, "Writing checkpoint for MWPvmRC layer.\n" );

	fprintf ( fp, "%d \n", exec_classes );
	fprintf ( fp, "%d \n", num_arches );
	for ( i = 0 ; i<num_arches ; i++ ) 
		{
			if ( arch_class_attributes && arch_class_attributes[i] )
				fprintf ( fp, "1 %d %s\n", strlen(arch_class_attributes[i]), arch_class_attributes[i] );
			else
				fprintf ( fp, "0 " );
		}

	fprintf ( fp, "%d \n", num_executables );
	for ( i = 0; i < num_executables; i++ )
		{
			fprintf ( fp, "%d %d %s %s\n", worker_executables[i]->arch_class, 
					  worker_executables[i]->exec_class, worker_executables[i]->executable, 
					  worker_executables[i]->executable_name );
			if ( worker_executables[i]->attributes != NULL )
				fprintf ( fp, "1 %d %s\n", strlen(worker_executables[i]->attributes), 
						  worker_executables[i]->attributes );
			else
				fprintf ( fp, "0 " );
		}

	for ( i = 0; i < exec_classes; i++ )
		{
			fprintf ( fp, "%d ", exec_class_target_num_workers[i] );
		}
	fprintf ( fp, "\n");

	write_RMstate ( fp );
	return 0;
}

int MWRMComm::read_checkpoint( FILE *fp )
{
	int temp;
	int i;
	MWprintf ( 50, "Reading checkpoint in MWPvmRC layer.\n" );
	fscanf ( fp, "%d %d", &exec_classes, &num_arches );
	MWprintf ( 50, "Target num workers: %d    Num arches %d.\n", 
			   target_num_workers, num_arches );
	arch_class_attributes = new char*[num_arches];
	target_num_workers = 0;
	for ( i = 0; i < num_arches; i++ )
		{
			fscanf ( fp, "%d ", &temp );
			if ( temp == 1 )
				{
					fscanf ( fp, "%d ", &temp );
					arch_class_attributes[i] = new char[temp+1];
					fgets ( arch_class_attributes[i], temp + 1, fp );
				}
			else
				arch_class_attributes[i] = NULL;
		}

	fscanf ( fp, "%d ", &num_executables );
	MWprintf ( 50, "%d Worker executables:\n", num_executables );
	worker_executables = new struct RMC_executable*[num_executables];
	for ( i = 0 ; i < num_executables; i++ ) 
		{
			worker_executables[i] = new struct RMC_executable;
			worker_executables[i]->executable = new char[_POSIX_PATH_MAX];
			worker_executables[i]->executable_name = new char[_POSIX_PATH_MAX];
			fscanf ( fp, "%d %d %s %s", &worker_executables[i]->arch_class,
					 &worker_executables[i]->exec_class, worker_executables[i]->executable,
					 worker_executables[i]->executable_name );
			fscanf ( fp, "%d ", &temp );
			if ( temp == 1 )
				{
					fscanf ( fp, "%d ", &temp );
					worker_executables[i]->attributes = new char[temp+1];
					fgets ( worker_executables[i]->attributes, temp + 1, fp );
				}
			else
				worker_executables[i]->attributes = NULL;
		}

	MW_exec_class_num_workers = new int[exec_classes];
	exec_class_target_num_workers = new int[exec_classes];
	for ( i = 0; i < exec_classes; i++ )
		{
			fscanf ( fp, "%d ", &exec_class_target_num_workers[i] );
			target_num_workers += exec_class_target_num_workers[i];
			MW_exec_class_num_workers[i] = 0;
		}

	read_RMstate ( fp );
	return 0;
}

int MWRMComm::read_RMstate( FILE *fp ) 
{ 
	return 0;
}

int MWRMComm::write_RMstate( FILE *fp )
{
	return 0;
}


// An exec class is a partition of the total number
// of running workers, and can be of Any architecture 
void
MWRMComm::set_num_exec_classes ( int num )
{
	exec_classes = num;

	if (exec_class_target_num_workers) {
		delete [] exec_class_target_num_workers;
		delete [] MW_exec_class_num_workers;
	}

	exec_class_target_num_workers = new int[num];
	MW_exec_class_num_workers = new int[num];
	for ( int i = 0; i < num; i++ )
		{
			exec_class_target_num_workers[i] = 0;
			MW_exec_class_num_workers[i] = 0;
		}
	
}

int
MWRMComm::get_num_exec_classes ( )
{
	return exec_classes;
}

// arch_classes are the different opsys/arch combos
// we will be running with
void
MWRMComm::set_num_arch_classes ( int num )
{
	char **tmp_arch_class_attributes = new char*[num];

	MWprintf(80, "MWRMComm::set_num_arch_classes to %d\n", num);
	for (int i = 0 ; i < num; i++) {
		tmp_arch_class_attributes[i] = NULL;
	}

	if (arch_class_attributes) {
		for (int i = 0; i < num_arches; i++) {
			tmp_arch_class_attributes[i] = arch_class_attributes[i];
		}
		delete [] arch_class_attributes;
	}

	num_arches = num;
	arch_class_attributes = tmp_arch_class_attributes;
}

void
MWRMComm::set_arch_class_attributes ( int arch_class, const char *attr )
{
	MWprintf(80, "MWRMComm::set_arch_class_attributes for arch class %d to %s\n", arch_class, attr);
	if ( arch_class >= num_arches ) {
		MWprintf(10, "Warning: arch_class %d >= num_arches (%d)\n", arch_class, num_arches);
		exit(-1);
		return;
	}
	arch_class_attributes[arch_class] = new char[strlen(attr)+1];
	strcpy ( arch_class_attributes[arch_class], attr );
}

int
MWRMComm::get_num_arch_classes ( )
{
	return num_arches;
}

void
MWRMComm::set_num_executables ( int num )
{
	int i = 0; 

	struct RMC_executable **tmp_executables = new RMC_executable*[num];

	if (worker_executables != NULL) {
		for ( i = 0; i < num_executables; i++ ) {
			tmp_executables[i] = worker_executables[i];
				//Jeff changed this.
			worker_executables[i] = NULL;
		}
		
		delete [] worker_executables;
	}

	worker_executables = tmp_executables;
	num_executables = num;
}

int
MWRMComm::get_num_executables() {
	return num_executables;
}

void MWRMComm::add_executable( int exec_class, int arch_class, char *exec_name, 
							   char *requirements )
{
	if ( !exec_name )
		return;

	if( arch_class >= num_arches ) 
		{
			MWprintf( 10, "set_worker_attributes(): incrementing num_arches to %d\n", 
					  arch_class + 1 );
			num_arches = arch_class + 1;
		}

	MWprintf(31, "tempnum_executables = %d\n", tempnum_executables);
	worker_executables[tempnum_executables] = new struct RMC_executable;
	worker_executables[tempnum_executables]->arch_class = arch_class;
	worker_executables[tempnum_executables]->exec_class = exec_class;
	worker_executables[tempnum_executables]->executable = new char [ strlen(exec_name) + 1 ];
	strcpy ( worker_executables[tempnum_executables]->executable, exec_name );
	worker_executables[tempnum_executables]->executable_name = 
		process_executable_name ( worker_executables[tempnum_executables]->executable, 
								  exec_class, arch_class );
	if ( !requirements )
		{
			worker_executables[tempnum_executables]->attributes = NULL;
		}
	else
		{
			worker_executables[tempnum_executables]->attributes = new char [ strlen(requirements) + 1 ];
			strcpy ( worker_executables[tempnum_executables]->attributes, requirements );
		}
	tempnum_executables++;
}

// This is the new, easy way to add executables.  No more worrying about getting
// calling order right, or setting up the wrong number of arrays.
// We assume that this will be called a couple of times at startup, so 
// efficiency isn't a big concern.

void
MWRMComm::add_executable(const char *exec_name, const char *requirements)
{
	verify_file_exists(exec_name);
		// exec class is the notion that you want to partition
		// the number of executables into two or more groups, and
		// try to maintain a different number in each group.

		// Most codes don't need this, so always force it to one,
		// no matter how many executables we have
	set_num_exec_classes(1);

		// make a new arch_class, and give it this requirements
	int num_arches = get_num_arch_classes();
	if (num_arches == -1) {
		num_arches = 0;
	}
	set_num_arch_classes(1 + num_arches);
	set_arch_class_attributes(get_num_arch_classes() - 1, requirements);

		// make a new executable
	set_num_executables(1 + get_num_executables());

	RMC_executable *e = new struct RMC_executable;
	e->exec_class = 0;
	e->arch_class = get_num_arch_classes() - 1;
	e->executable   = strdup(exec_name);
	e->attributes = ""; // Where is this used?
	e->executable_name = 
		process_executable_name ( e->executable, 
								  e->exec_class, e->arch_class );
	worker_executables[get_num_executables() - 1] = e;
}

char* MWRMComm::process_executable_name( char *exec_name, int ex_cl, int ar_cl )
{
	char *newone = new char[strlen(exec_name) + 1];
	strcpy ( newone, exec_name );
	return newone;
}

// XXX Newly added get function so that we can know the target number of workers
//   and change the number for runtime adaptation

int 
MWRMComm::get_target_num_workers ( int exec_class )
{
	if (exec_class == -1)
		return target_num_workers;
	
	if ( (exec_class > exec_classes) || (exec_class < -1) ) 
		return -1;

	return exec_class_target_num_workers[exec_class];
}

void 
MWRMComm::set_worker_increment(int newinc)
{
	if (newinc <= 0) {
		MWprintf(1, "Unreasonable value: %d for set_worker_increment().  Ignoring.\n", newinc); 
		return;
	}
	hostinc_ = newinc < target_num_workers ? newinc : target_num_workers;
}

int 
MWRMComm::get_worker_increment() const
{
	return hostinc_; 
}

void
MWRMComm::set_target_num_workers ( int num_workers )
{
	set_target_num_workers ( -1, num_workers );
}

void
MWRMComm::set_target_num_workers ( int exec_class, int num_workers )
{
	if ( exec_classes <= 0 )
		{
			exec_classes = 1;
			exec_class_target_num_workers = new int[exec_classes];
			exec_class_target_num_workers[0] = 0;
		}

	if ( exec_class < 0 )
		{
			exec_class_target_num_workers[0] = num_workers;
			target_num_workers = exec_class_target_num_workers[0];
			return;
		}

	target_num_workers -= exec_class_target_num_workers[exec_class];
	exec_class_target_num_workers[exec_class] = num_workers;
	target_num_workers += exec_class_target_num_workers[exec_class];
	return;
}

void MWRMComm::set_worker_checkpointing( bool wc )
{
	if( wc == true )
		MWprintf( 10, "Warning!  Worker checkpointing not available in this CommRM implementation\n" );
	worker_checkpointing = false;
}

void 
MWRMComm::verify_file_exists(const char *exec_name) {
	FILE *f = fopen(exec_name, "r");
	if (f == NULL) {
		MWprintf(10, "ERROR: add_executable with %s, which doesn't exist or isn't readable\n", exec_name);
	} else {
		fclose(f);
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
  eval: (setq c-tab-always-indent nil)
  eval: (setq tab-width 4)
  End:
*/
