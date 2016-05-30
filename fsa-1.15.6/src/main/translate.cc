
/**
 * \file translate.cc
 * This file is part of FSA.
 * \author Source code in this file was written by Robert Bradley.
 */

#include <getopt.h>

#include "config.h"
#include "seq/sequence.h"
#include "seq/alignment.h"

using namespace fsa;

static std::string program_name = "translate";

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
    << "Usage: " << program_name << " <FASTA file>" << endl
    << endl
    << "Translate input nucleotide sequences into protein sequence." << endl
    << endl
    << "Input nucleotide sequences must be in FASTA format." << endl
    << "Uses the first reading frame and drops incomplete codons at ends of sequences." << endl
    << endl;
}

int main (int argc, char** argv) {

  std::string filename;

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

  filename = std::string (argv[optind++]);

  assert (optind == argc);

  // now do stuff!
  Sequence_database seq_db;
  if (!Sequence::detect_fasta (filename)) {
    cerr << "ERROR: Input nucleotide sequences must be in FASTA format.";
    return 1;
  }
  seq_db.read_fasta (filename, Alignment::is_gap_char,
		     false, // strip_leading_chr
		     true); // verbose

  // translate sequences & display  
  seq_db.translate().write_fasta (cout);

  return 0;
    
}
