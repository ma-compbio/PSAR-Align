/*
 *  This file is part of FSA, a sequence alignment algorithm.
 *  Source code in this file was written by Jaeyoung Do.
 */

#ifndef DB_ADAPTER_INCLUDED
#define DB_ADAPTER_INCLUDED

#include "seq/sequence.h"
#include "util/hash_fcn.h"
#include "manager/db_misc.h"
#ifdef HAVE_POSTGRES
#include "manager/db_postgres.h"
#endif

namespace fsa {

  /// Class DB_adapter
  /*
   * This is a bridge class between the Manager class and Database classes. 
   * The purpose of this class is to support many different databases. 
   * Note that currently only PostgresSQL is supported. 
   *
   * The organization of db tables and schemas
   *
   * fsa_root +- schema_1 ------+-- sequence_info                 // for input sequences
   *          +- schema_2       +-- parameter_info                // for parameters of the input sequences  
   *          +- ...            +-- h1_0, h1_1, ...               // heap tables for candidate edges
   *                            +-- ncells1_0, ncells1_1, ...     // num cells tables
   *                            +-- s1_0, s1_1, ...               // sparse matrix tables for pairwise posterior probabilities
   *
   * Note 1) schema_x (i.e., x is the sequence schema id) matches the unique set of sequences 
   * Note 2) hy_z, ncellsy_z, sy_z: y is the unique parameter table id, and 0 <= z <= # of workers that have been used 
   *                                when generating the data
   */

  class DB_adapter
  {
  public:
    /// Constructor
    DB_adapter();

    /// Distructor
    ~DB_adapter();

    /// attempt to establish a connection to the given database host
    /*
     * With given information, try to make a connection to database server.
     * Note that to establish a connection, at least db_hostname (or db_hostaddr) and db_name
     * must be provided.
     */
    bool connect_db (const char *db_hostname, const char* db_hostaddr, const char *db_name, 
		     const int db_port, const char *db_user, const char *db_password);

    /// disconnect the connection
    void disconnect_db (); 

    /// look up data to be used when the sequence annealing.
    /*
     * seq_db_internal contains information about the input sequence strings. 
     * db_opts contains the options that have been used when running fsa
     */
    bool look_up_data (const Sequence_database &seq_db_internal, const DB_opts &db_opts);

    /// initialize the database 
    /* 
     * If necessary, create tables and schemas, and insert tuples to the tables.
     */
    bool init_database(const Sequence_database &seq_db_internal, const DB_opts &db_opts); 

    ///  set up various ids
    /*
     * <seqs_schema_id, params_table_id> is the unique pair of ids denoting 
     * specific sequences with the set of options that has been used when 
     * the data was generated. 
     */
    bool set_up_ids (const int seqs_schema_id, const int params_table_id, const int worker_id);

#ifdef HAVE_POSTGRES

  private:
    int             m_seqs_schema_id;       /// sequence schema id
    int             m_params_table_id;      /// parameter table id
    DB_postgres*    m_db_connection;        /// db connection object

  public:

    /// get the db connection
    DB_postgres* get_connection ();

    /// true if the database server is running
    bool is_db_running ();

    /// is data available? true if (is_seqs_available() && is_params_available())
    bool is_data_available ();

    /// is input sequences available?
    bool is_seqs_available ();

    /// does the data in database correspond with the input parameters?
    bool is_params_available ();
#endif

  };

}

#endif
