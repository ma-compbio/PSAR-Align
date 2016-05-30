/*
 *  This file is part of FSA, a sequence alignment algorithm.
 *  Source code in this file was written by Jaeyoung Do.
 */

#include "manager/db_adapter.h"

using namespace fsa;
 
DB_adapter::DB_adapter() {
#ifdef HAVE_POSTGRES
	m_seqs_schema_id  = DB_NOT_AVAILABLE;
	m_params_table_id = DB_NOT_AVAILABLE;

	m_db_connection = NULL;
#endif
}

DB_adapter::~DB_adapter() {
	// disconnect the connection
	disconnect_db();
}

void DB_adapter::disconnect_db () {
#ifdef HAVE_POSTGERS
	if (!m_db_connection) {
		// if there's previous connection garbage, delete it and create new one 
		m_db_connection->disconnect_db();
		delete m_db_connection;
	}
	
	// set various ids as "UNDETERMINED" ids
	m_seqs_schema_id = DB_NOT_AVAILABLE;
	m_params_table_id = DB_NOT_AVAILABLE;
	m_db_connection = NULL;

#endif

}

bool DB_adapter::connect_db (const char *db_hostname, const char *db_hostaddr, const char* db_name, 
		const int db_port, const char *db_user, const char *db_password) {
#ifndef HAVE_POSTGRES
	return DB_BAD;
#else
	m_db_connection = new DB_postgres();

	// can we establish a connection to the database server?
	bool is_db_available = m_db_connection->connect_db(db_hostname, db_hostaddr, db_name, db_port, db_user, db_password);

	if ( !is_db_available) {
		// in the case that the database we are willing to connect is not running (or does not response)
		delete m_db_connection;
		m_db_connection = NULL;
	}

	return is_db_available;
#endif
}

bool DB_adapter::set_up_ids (const int seqs_schema_id, const int params_table_id, const int worker_id) {
#ifndef HAVE_POSTGRES
	return DB_BAD;
#else
	// set sequence schema id
	m_db_connection->set_seqs_schema_id (seqs_schema_id);

	// set parameter table id
	m_db_connection->set_params_table_id (params_table_id);

	// set worker id
	m_db_connection->set_worker_id (worker_id);

	m_seqs_schema_id = seqs_schema_id;
	m_params_table_id = params_table_id;

	return DB_OK;
#endif
}

bool DB_adapter::init_database(const Sequence_database &seq_db_internal, const DB_opts &db_opts) {

#ifndef HAVE_POSTGRES
	return DB_BAD;
#else
	if ( is_data_available() )  {
		// in the case of the input sequences that has been previously considered
		m_db_connection->update_fsa_table ();
		m_db_connection->update_params_table (db_opts.num_parallelized_jobs);
	}
	else {
		// if the input sequences are available - but possibly there's no relevant parameters
		if ( is_seqs_available() ) 
			m_db_connection->update_fsa_table ();
		else {
		  std::string whole_seqs    = "";
			float seqs_avr_length = 0.0;
			uint32_t hash_key     = 0;
			int num_seqs          = seq_db_internal.size();

			// generate a hash key based on the input sequences
			for (int i=0; i<num_seqs; i++) {
			  whole_seqs += seq_db_internal.get_seq (i).seq;
			  seqs_avr_length += (float) seq_db_internal.get_seq (i).length();
			}
			seqs_avr_length /= (float) num_seqs; 
			seqs_avr_length = (int) seqs_avr_length;
			hash_key = Hash_functions::hsieh_hash(whole_seqs.c_str());
			
			// create fsa main table "public.fsa_root"
			m_db_connection->create_fsa_table ();
			
			// add a row to the fsa main table for the input sequences
			m_db_connection->insert_fsa_schema (hash_key, num_seqs, seqs_avr_length);

			// set the id
			m_seqs_schema_id = m_db_connection->get_seqs_schema_id (hash_key, num_seqs, seqs_avr_length);
			
			if (m_seqs_schema_id != DB_NOT_AVAILABLE) 
				m_db_connection->set_seqs_schema_id (m_seqs_schema_id);

			// create a new schema for the input sequences with the hash key
			m_db_connection->create_seqs_schema ();
			
			// create a parameter table in the schema
			m_db_connection->create_params_table ();
			
			// create seqs table - for information of the input sequences
			m_db_connection->create_seqs_table ();
		
			// insert input sequences into the seqs table
			m_db_connection->insert_seqs_table (seq_db_internal);
		}

		// insert parameter info. into the parameter table
		m_db_connection->insert_params_table (db_opts);

		// set the parameter id
		m_params_table_id = m_db_connection->get_params_table_id (db_opts);
		if (m_params_table_id != DB_NOT_AVAILABLE) 
			m_db_connection->set_params_table_id (m_params_table_id);
	}

	for (int i=0; i<db_opts.num_parallelized_jobs; i++) {

		// create sparse matrix tables to store pairwise posterior probabilities
		if ( m_db_connection->create_sparse_matrix_table (i) == DB_BAD ) {
			m_db_connection->delete_sparse_matrix_table (i); 
			m_db_connection->drop_sparse_matrix_table_index (i);
		}

		// create num cells tables
		// num cells table : how many entries are in each sequence pair
		// Note that seq1 = 0 and seq2 = 0 indicates the number of edges
		if ( m_db_connection->create_num_cells_table (i) == DB_BAD ) 
			m_db_connection->delete_num_cells_table (i);

		// create heap tables to store candidate edges of the null alignment
		if ( m_db_connection->create_heap_table(i) == DB_BAD ) {
			m_db_connection->delete_heap_table (i);
		}

		// create merged heap table
		// merged heap table: when doing the sequence annealing, the merged heap table acts like the priority queue
		if ( m_db_connection->create_merged_heap_table() == DB_BAD ) {
			m_db_connection->delete_merged_heap_table();
			m_db_connection->drop_merged_heap_table_index ();
		}
	}
	
	return DB_OK;
#endif

}

bool DB_adapter::look_up_data (const Sequence_database &seq_db_internal, const DB_opts &db_opts) {
#ifndef HAVE_POSTGRES
	return DB_BAD;
#else
	// to look up data, the database must be running 
	if (!is_db_running()) 
		return DB_NO;

	std::string whole_seqs    = "";
	float seqs_avr_length = 0.0;
	uint32_t hash_key     = 0;

	// generate the hash key based on input sequences

	for (int i=0; i<seq_db_internal.size(); i++) {
	  whole_seqs += seq_db_internal.get_seq (i).seq;
	  seqs_avr_length += (float) seq_db_internal.get_seq (i).length();
	}

	seqs_avr_length /= (float) seq_db_internal.size(); 
	seqs_avr_length = (int) seqs_avr_length;

	hash_key = Hash_functions::hsieh_hash(whole_seqs.c_str());

	// set the sequence table id and the parameter id as looking up the appropriate data based 
	// on the input sequences (which can be referenced by seq_db_internal) and the options 
	// (which can be gotten from db_opts) that have been used when running fsa.  

	if ( (m_seqs_schema_id = m_db_connection->get_seqs_schema_id (hash_key, seq_db_internal.size(), seqs_avr_length) ) 
			!= DB_NOT_AVAILABLE) {
		m_db_connection->set_seqs_schema_id (m_seqs_schema_id);

		if ( (m_params_table_id = m_db_connection->get_params_table_id (db_opts) ) 
				!= DB_NOT_AVAILABLE) 
			m_db_connection->set_params_table_id (m_params_table_id);
	}

	// see if there are tables for input sequences and input parameters of the sequences. 	
	if ( m_seqs_schema_id != DB_NOT_AVAILABLE && m_params_table_id != DB_NOT_AVAILABLE)
		return DB_OK;
	else
		return DB_BAD;
#endif

}

#ifdef HAVE_POSTGRES
/// get the connection
DB_postgres* DB_adapter::get_connection () {
	return m_db_connection;
}

/// is the database server is running?
bool DB_adapter::is_db_running () 
{
	if ( m_db_connection )
		return DB_YES;
	else 
		return DB_NO;
}

/// data is available only if there're tables related to input sequences and parameters of the input sequences.
bool DB_adapter::is_data_available () 
{
	if ( is_seqs_available () && is_params_available () )
		return DB_YES;
	else 
		return DB_NO;
}

/// is seqs available?
bool DB_adapter::is_seqs_available () 
{
	if ( m_seqs_schema_id != DB_NOT_AVAILABLE )
		return DB_YES;
	else
		return DB_NO;
}

/// is parameters available?
bool DB_adapter::is_params_available () 
{
	if (m_params_table_id != DB_NOT_AVAILABLE )
		return DB_YES;
	else
		return DB_NO;
}
#endif
