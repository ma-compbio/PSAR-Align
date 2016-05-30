
/**
 * \file gapcleaner.cc
 * This file is part of FSA, a sequence alignment algorithm.
 * \author Source code in this file was written by Robert Bradley.
 */

#include <getopt.h>

#include "config.h"
#include "annealing/alignment_DAG.h"

using namespace fsa;

static std::string program_name = "gapcleaner";

// short options letters
const char* const short_options = "vh";

// long options
static struct option long_options[] = {
  {"version", no_argument, NULL, 'v'},
  {"help", no_argument, NULL, 'h'},

  /* required NULL termination of array */
  {NULL, 0, NULL, 0}
};

static void print_version (std::ostream& o) {
  o << program_name << " from " << PACKAGE_STRING << endl;
}

static void print_usage (std::ostream& o) {
  print_version (o);
  o << endl
    << "Usage: " << program_name << " <multi-FASTA or Stockholm alignment>" << endl
    << endl
    << "Find the most-parsimonious ordering of indels." << endl
    << "Finds the minimal chain decomposition of the POSET of the alignment;" << endl
    << "this corresponds to minimizing the number of gap-openings across the sequences." << endl
    << endl
    << "Input protein alignment must be in multi-FASTA or Stockholm format." << endl
    << endl;
}


int main (int argc, char** argv) {

  std::string alignfilename;

  // parse options
  int c;
  int option_index = 0; // getopt_long stores option index here
  while (1) {

    // get next option
    c = getopt_long (argc, argv, short_options,
		     long_options, &option_index);

    if (c == -1)
      break;

    switch (c) {

    case 'v': // version message
      print_version (cout);
      return 0;

    case 'h': // help message
      print_usage (cout);
      return 0;

    case '?': // invalid option
      print_usage (cerr);
      return 1;

    default:  // unexpected
      abort();

    }
  }

  // stuff
  if (optind + 1 != argc) {
    print_usage (cerr);
    return 1;
  }

  alignfilename = std::string (argv[optind++]);

  assert (optind == argc);

  // now do stuff!
  Sequence_database seq_db;
  Stockholm stockholm (seq_db);

  bool is_mfa = false;
  if (Alignment::detect_mfa (alignfilename))
    is_mfa = true;
  if (is_mfa)
    stockholm.read_mfa (alignfilename);
  else
    stockholm.read_stockholm (alignfilename);
    
  // perform greedy DFS decomposition
  Alignment_DAG dag (seq_db, stockholm);
  dag.dfs_topological_sort();

  // show results
  if (is_mfa)
    dag.get_stockholm().write_mfa (cout);
  else
    dag.get_stockholm().write_stockholm (cout);

  return 0;
}
