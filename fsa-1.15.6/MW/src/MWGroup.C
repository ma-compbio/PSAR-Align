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

#include "MWGroup.h"

MWGroup::MWGroup ( int grps )
{
	maxGroups = grps;
	group = new bool[maxGroups];
	for ( int i = 0; i < maxGroups; i++ )
		group[i] = FALSE;
}

MWGroup::~MWGroup ( )
{
	delete []group;
}

void
MWGroup::join ( int num )
{
	if ( num >= 0 && num < maxGroups )
		group[num] = TRUE;
}

void
MWGroup::leave ( int num )
{
	if ( num >= 0 && num < maxGroups )
		group[num] = FALSE;
}

bool
MWGroup::belong ( int num )
{
	if ( num >= 0 && num < maxGroups )
		return group[num];

	return FALSE;
}

bool
MWGroup::doesOverlap ( MWGroup *grp )
{
	for ( int i = 0; i < maxGroups; i++ )
	{
		if ( belong ( i ) && grp->belong ( i ) )
			return TRUE;
	}
	return FALSE;
}

void
MWGroup::write_checkpoint ( FILE *fp )
{
	int num = 0;
	int i;
	for ( i = 0; i < maxGroups; i++ )
		if ( belong ( i ) ) num++;
	fprintf ( fp, "%d ", num );

	for ( i = 0; i < maxGroups; i++ )
		if ( belong ( i ) )
			fprintf ( fp, "%d ", i );
	fprintf ( fp, "\n" );
}

void
MWGroup::read_checkpoint ( FILE *fp )
{
	int i, temp;
	int gp;
	fscanf ( fp, "%d ", &temp );
	for ( i = 0; i < temp; i++ )
	{
		fscanf ( fp, "%d ", &gp );
		join ( gp );
	}
}
