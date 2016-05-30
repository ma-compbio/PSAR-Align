
/**
 * \file slice_fasta.cc
 * This file is part of FSA.
 * \author Source code in this file was written by Robert Bradley.
 */

#include <getopt.h>

#include "config.h"
#include "seq/gff.h"
#include "seq/sequence.h"
#include "seq/alignment.h"

using namespace fsa;

static std::string program_name = "slice_fasta";

static std::string type = "";

// short options letters
const char* const short_options = "vht:";

// long options
static struct option long_options[] = {
  {"version", no_argument, NULL, 'v'},
  {"help", no_argument, NULL, 'h'},

  {"type", required_argument, NULL, 't'},

  /* required NULL termination of array */
  {NULL, 0, NULL, 0}
};

static void print_version (std::ostream& o) {
  o << program_name << " from " << PACKAGE_STRING << endl;
}

static void print_usage (std::ostream& o) {
  print_version (o);
  o << endl
    << "Usage: " << program_name << " [options] <FASTA file> <GFF file>" << endl
    << endl
    << "Options:" << endl
    << "    -t, --type <string>         only look at features of particular type" << endl
    << endl
    << "Slice subsequences from an input FASTA file." << endl
    << "Assumes 1-based, fully-closed coordinates." << endl
    << "If no strand information is available, then the + strand is assumed." << endl
    << "If the - strand is requested, then the subsequence is pulled out and then reverse-complemented." << endl
    << endl
    << "Assumes a DNA alphabet." << endl
    << endl;
}

int main (int argc, char** argv) {

  std::string filename;
  std::string gff_filename;

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

    case 'h': // help message
      print_usage (cout);
      return 0;

    case 'v': // version message
      print_version (cout);
      return 0;

    case 't': // type
      type = std::string (optarg);
      break;

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

  filename = std::string (argv[optind++]);
  gff_filename = std::string (argv[optind++]);

  assert (optind == argc);

  // read in sequence
  Sequence_database seq_db;
  if (!Sequence::detect_fasta (filename)) {
    cerr << "ERROR: Input nucleotide sequences must be in FASTA format.";
    return 1;
  }
  seq_db.read_fasta (filename, Alignment::is_gap_char,
		     true,  // strip_leading_chr
		     true); // verbose

  // read GFF file
  GFF_database gff_db;
  gff_db.from_file (gff_filename,
		    true); // strip_leading_chr

  // initialize alphabet for reverse-complementing
  const DNA_alphabet dna_alphabet;

  // iterate through GFF
  for (std::vector<GFF>::const_iterator gff = gff_db.begin(); gff != gff_db.end(); ++gff) {
    
    // if requested, only look at features of a particular type
    if (type.length() && type != gff->type)
      continue;

    // confirm that we have sequence data for this chromosome
    if (!seq_db.exists_seq (gff->seqid)) {
      cerr << "WARNING: No sequence data found for requested chromosome '" << gff->seqid << "'." << endl;
      continue;
    }

    // get sequence data
    const Sequence& sequence = seq_db.get_seq (gff->seqid);

    // convert to 0-based coordinates
    const size_t start = gff->start - 1;
    const size_t end = gff->end - 1;

    // check coordinates sane
    if (gff->start - 1 >= sequence.length()) {
      cerr << "ERROR: Requested coordinates are out of bounds:" << endl
	   << *gff << endl;
      exit (1);
    }

    // take subsequence
    Sequence* subseq = sequence.subsequence (start, end);

    // replace chromosome name with GFF entry
    subseq->name = gff->to_string();

    // reverse-complement if relevant
    if (gff->strand == '-')
      subseq->revcomp (dna_alphabet);

    // display  
    subseq->write_fasta (cout);

    // clean up
    delete subseq;

  }

  return 0;
    
}
