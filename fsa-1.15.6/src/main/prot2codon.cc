
/**
 * \file prot2codon.cc
 * This file is part of FSA.
 * \author Source code in this file was written by Robert Bradley.
 */

#include <getopt.h>

#include "config.h"
#include "seq/sequence.h"
#include "seq/alignment.h"

using namespace fsa;

static std::string program_name = "prot2codon";

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
    << "Usage: " << program_name << " <multi-FASTA or Stockholm alignment> <FASTA file>" << endl
    << endl
    << "Find the codon alignment corresponding to the given protein alignment." << endl
    << endl
    << "Input protein alignment must be in multi-FASTA or Stockholm format." << endl
    << "Input nucleotide sequences must be in FASTA format." << endl
    << endl;
}

int main (int argc, char** argv) {

  std::string alignfilename;
  std::string seqfilename;

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
  if (optind + 2 != argc) {
    print_usage (cerr);
    return 1;
  }

  alignfilename = std::string (argv[optind++]);
  seqfilename = std::string (argv[optind++]);

  assert (optind == argc);

  // now do stuff!

  // read in protein alignment
  Sequence_database seq_db_aa;
  Stockholm stockholm_aa (seq_db_aa);
  bool is_mfa = false;
  if (Alignment::detect_mfa (alignfilename))
    is_mfa = true;
  if (is_mfa)
    stockholm_aa.read_mfa (alignfilename);
  else
    stockholm_aa.read_stockholm (alignfilename);
  
  // check correct alphabet
  Protein_alphabet prot = Protein_alphabet();
  if (!seq_db_aa.matches_alphabet (prot, 0.75)) { // use loose threshold for calling it protein
    cerr << "ERROR: Input alignment doesn't seem to be protein sequence; I'm bailing." << endl;
    return 1;
  }

  // read in nucleotide fasta file
  if (!Sequence::detect_fasta (seqfilename)) {
    cerr << "ERROR: Input nucleotide sequences must be in FASTA format.";
    return 1;
  }
  Sequence_database seq_db_codon;
  seq_db_codon.read_fasta (seqfilename, Alignment::is_gap_char,
			   false, // strip_leading_chr
			   true); // verbose

  // translate sequences & display  
  Stockholm stockholm_codon = stockholm_aa.get_codon_from_aa_alignment (seq_db_codon);
  if (is_mfa)
    stockholm_codon.write_mfa (cout);
  else
    stockholm_codon.write_stockholm (cout);

  return 0;
    
}
