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
/* MWTask.C

   The implementation of the MWTask class

*/

#include "MWTask.h"
#include "MWDriver.h"
#include <stdio.h>

extern int MWworkClasses;
extern int *MWworkClassTasks;

MWTask::MWTask() {
	number = -1;
	numsubtask = -1;
	worker = NULL;
	taskType = MWNORMAL;
	group = new MWGroup ( MWworkClasses );
}

MWTask::~MWTask() {
  if ( worker )
    if ( worker->runningtask )
      worker->runningtask = NULL;

  delete group;
}

void MWTask::printself( int level ) 
{    
	MWprintf ( level, "  Task %d\n", number);
}

void
MWTask::initGroups ( int num )
{
	// found LEAK_ASSIGN (group probably already allocated!)
	// group = new MWGroup ( num );
	if (group != NULL)
		delete group;
	group = new MWGroup( num );
}

void
MWTask::addGroup ( int num )
{
	if ( doesBelong ( num ) )
		return;
	MWworkClassTasks[num]++;
	group->join ( num );
}

void
MWTask::deleteGroup ( int num )
{
	if ( !doesBelong ( num ) )
		return;
	MWworkClassTasks[num]--;
	group->leave ( num );
}

bool
MWTask::doesBelong ( int num )
{
	return group->belong ( num );
}

MWGroup*
MWTask::getGroup ( )
{
	return group;
}

void
MWTask::write_group_info ( FILE *fp )
{
	group->write_checkpoint ( fp );
}

void
MWTask::read_group_info ( FILE *fp )
{
	// found LEAK_ASSIGN (group probably already allocated!)
        // group = new MWGroup ( num );
	if (group != NULL)
		delete group;
	group = new MWGroup ( MWworkClasses );
	group->read_checkpoint ( fp );
	for ( int i = 0; i < MWworkClasses; i++ )
	{
		if ( doesBelong ( i ) )
			addGroup ( i );
	}
}

void 
MWTask::completedTask(double wall_time, double cpu_time)
{
	worker->runningtask = NULL;
	worker->completedtask( wall_time, cpu_time);
	worker = NULL;

	for ( int tempi = 0; tempi < MWworkClasses; tempi++ ) {
		deleteGroup ( tempi );
	}
}
