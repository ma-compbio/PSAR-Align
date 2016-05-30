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
#ifndef MWNWSTASK_H
#define MWNWSTASK_H

#include <stdio.h>
#include "MWTask.h"

/** The task class specific to NWS weather monitoring
*/

class MWNWSTask : public MWTask {

 public:
    
    /// Default Constructor
    MWNWSTask();

    /** 
    */
    MWNWSTask( int, double, int, int, char* );

    MWNWSTask( const MWNWSTask& );

    /// Default Destructor
    ~MWNWSTask();

    MWNWSTask& operator = ( const MWNWSTask& );

    /**@name Implemented methods
       
       These are the task methods that must be implemented 
       in order to create an application.
    */
    //@{
    /// Pack the work for this task into the PVM buffer
    void pack_work( void );
    
    /// Unpack the work for this task from the PVM buffer
    void unpack_work( void );
    
    /// Pack the results from this task into the PVM buffer
    void pack_results( void );
    
    /// Unpack the results from this task into the PVM buffer
    void unpack_results( void );
    //@}

	// not used in this class
    void pack_subresults( int );
    void unpack_subresults( int );   

        /// dump self to screen:
    void printself( int level = 60 );

    /**@name Checkpointing Implementation 
       These members used when checkpointing. */
    //@{
        /// Write state
    void write_ckpt_info( FILE *fp );
        /// Read state
    void read_ckpt_info( FILE *fp );
    //@}

    /**@name Input Parameters.
       These parameters are sent to the other side */
    //@{
    	/// the number of measurements.
    int iterations;
	/// the interval between measurements
    double timeInterval;
    	/// the port to connect to.
    int portNo;
    	/// the machine name.
    char machineAddress[1024];
    //@}

    /**@name Output Parameters.
       These parameters are returned by the workers */
    //@{
    	/// the maximum latency observed.
    double max;
    	/// the minimum latency observed.
    double min ;
    	/// the median latency observed.
    double median ;
    //@}
};

#endif
