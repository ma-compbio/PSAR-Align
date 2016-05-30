
/**
 * \file slice_fasta.cc
 * This file is part of FSA.
 * \author Source code in this file was written by Robert Bradley.
 */

#include <getopt.h>

#include "config.h"
#include "seq/sequence.h"
#include "seq/alignment.h"

using namespace fsa;

static std::string program_name = "slice_fasta";

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
    << "Usage: " << program_name << " <FASTA file> <sequence name> <start> <end> <strand>" << endl
    << endl
    << "Slice a subsequence from an input FASTA file." << endl
    << "Assumes 1-based, fully-closed coordinates." << endl
    << "If <strand> is omitted, then the + strand is assumed." << endl
    << "If the - strand is requested, then the subsequence is pulled out and then reverse-complemented." << endl
    << endl;
}

int main (int argc, char** argv) {

  std::string filename;
  std::string seqname;
  unsigned start;
  unsigned end;
  char strand = '+';

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
  if ((optind + 4 != argc) && (optind + 5 != argc)) {
    print_usage (cerr);
    return 1;
  }

  filename = std::string (argv[optind++]);
  seqname = std::string (argv[optind++]);
  Util::strip_leading_chr (seqname);
  start = static_cast<unsigned> (atoi (argv[optind++]));
  end = static_cast<unsigned> (atoi (argv[optind++]));
  // optionally get strand
  if (optind + 1 == argc)
    strand = argv[optind++][0];

  // check that arguments are sane
  if (strand != '+' && strand != '-') {
    print_usage (cerr);
    return 1;
  }

  // check that coordinates are sane
  if (end < start) {
    cerr << "ERROR: Please enter valid coordinates." << endl;
    return 1;
  }
  if (start == 0 || end == 0) {
    cerr << "ERROR: Please enter valid coordinates (remember that they are 1-based)." << endl;
    return 1;
  }

  assert (optind == argc);

  // now do stuff!
  Sequence_database seq_db;
  if (!Sequence::detect_fasta (filename)) {
    cerr << "ERROR: Input nucleotide sequences must be in FASTA format.";
    return 1;
  }
  seq_db.read_fasta (filename, Alignment::is_gap_char,
		     true,  // strip_leading_chr
		     true); // verbose

  // get requested subsequence
  if (!seq_db.exists_seq (seqname)) {
    cerr << "ERROR: No sequence named '" << seqname << "'." << endl
	 << "Available sequences are:" << endl
	 << Util::join (seq_db.get_sequence_list(), "\t") << endl;
    return 1;
  }
  const Sequence& sequence = seq_db.get_seq (seqname);

  // convert to 0-based coordinates
  --start;
  --end;

  // check coordinates sane
  if (start >= sequence.length()) {
    cerr << "ERROR: Requested coordinates are out of bounds." << endl;
    exit (1);
  }

  // take subsequence
  Sequence* subseq = sequence.subsequence (start, end);

  // reverse-complement if relevant
  const DNA_alphabet dna_alphabet;
  if (strand == '-') {
    if (seq_db.matches_alphabet (dna_alphabet))
      subseq->revcomp (dna_alphabet);
    else {
      cerr << "ERROR: Sequences in FASTA file don't match DNA alphabet, so I can't take the reverse-complement." << endl;
      return 1;
    }
  }

  // display  
  subseq->write_fasta (cout);

  // clean up
  delete subseq;

  return 0;
    
}
