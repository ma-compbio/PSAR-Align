
/**
 * \file slice_mercator_alignment.cc
 * This file is part of FSA.
 * \author Source code in this file was written by Robert Bradley.
 */

#include <getopt.h>
#include <fstream>

#include "config.h"
#include "seq/alignment.h"
#include "seq/mercator.h"

using namespace fsa;

static std::string program_name = "slice_mercator_alignment";

static std::string datadir = ".";
static std::string mapdir = ".";
static std::string aligndir = ".";
static bool stockholm = false;
static bool lazy = false;
static bool truncate_ok = false;
static bool zerobased = false;
static bool halfopen = false;

// short options letters
const char* const short_options = "vhD:M:A:LUs0o";

// long options
static struct option long_options[] = {
  {"version", no_argument, NULL, 'v'},
  {"help", no_argument, NULL, 'h'},

  {"data", required_argument, NULL, 'D'},
  {"map", required_argument, NULL, 'M'},
  {"align", required_argument, NULL, 'A'},
  {"lazy", no_argument, NULL, 'L'},
  {"truncate", no_argument, NULL, 'U'},

  {"stockholm", no_argument, NULL, 's'},

  {"zerobased", no_argument, NULL, '0'},
  {"halfopen", no_argument, NULL, 'o'},

  /* required NULL termination of array */
  {NULL, 0, NULL, 0}
};

static void print_version (std::ostream& o) {
  o << program_name << " from " << PACKAGE_STRING << endl;
}

static void print_usage (std::ostream& o) {
  print_version (o);
  o << endl
    << "Usage: " << program_name << " [options] <genome> <chromosome> <start> <end> <strand>" << endl
    << endl
    << "Extracts the corresponding subalignment from a Mercator multiple alignment." << endl
    << endl
    << "Options:" << endl
    << "    -h, --help                  show this message" << endl
    << endl
    << "    -D, --data <directory>      path to map, genome and alignment files" << endl
    << "    -M, --map <directory>       path to map and genome files" << endl
    << "    -A, --align <directory>     path to alignment files" << endl
    << "    -L, --lazy                  warn, rather than die, if the subalignment can't be obtained" << endl
    << "    -U, --truncate              truncate unmappable sequence (rather than skipping) and show truncated subalignment" << endl
    << endl
    << "    -s, --stockholm             use and display Stockholm-format alignments with conservation statistics (default is multi-FASTA)" << endl
    << endl
    << "    -0, --zerobased             coordinates are 0-based (default is 1-based)" << endl
    << "    -o, --halfopen              end coordinate is open, i.e., [start, end)" << endl
    << endl
    << "Assumes that coordinates are 1-based and fully-closed," << endl
    << "therefore representing the interval [start, end]." << endl
    << endl
    << "If requested, unmappable sequence will be truncated to the mappable portion;" << endl
    << "note that the truncation will favor the beginning of the requested sequence." << endl
    << endl
    << "If the requested sequence is on the - strand, then the corresponding" << endl
    << "subalignment will be reverse-complemented." << endl
    << endl;
}

int main (int argc, char** argv) {

  std::string genome;
  std::string chromosome;
  unsigned start;
  unsigned end;
  char strand;

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

    case 'D': // data directory
      datadir = std::string (optarg);
      break;

    case 'M': // map directory
      mapdir = std::string (optarg);
      break;

    case 'A': // align directory
      aligndir = std::string (optarg);
      break;

    case 'L': // lazy
      lazy = true;
      break;

    case 'U': // truncate
      truncate_ok = true;
      break;

    case 's': // stockholm
      stockholm = true;
      break;

    case '0': // zerobased
      zerobased = true;
      break;

    case 'o': // halfopen
      halfopen = true;
      break;

    case '?': // invalid option
      print_usage (cerr);
      return 1;

    default:  // unexpected
      abort();

    }
  }

  // get the rest of the command-line arguments
  if (optind + 5 != argc) {
    print_usage (cerr);
    return 1;
  }

  genome = std::string (argv[optind++]);
  chromosome = std::string (argv[optind++]);
  start = static_cast<unsigned> (atoi (argv[optind++]));
  end = static_cast<unsigned> (atoi (argv[optind++]));
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
  if (!zerobased && start == 0) {
    cerr << "ERROR:Please enter a valid start coordinate (remember that they are 1-based by default)." << endl;
    return 1;
  }

  assert (optind == argc);


  // now do stuff!
  Sequence_database seq_db;

  // map coordinates to make them 0-based and fully-closed
  if (!zerobased) {
    --start; --end;
  }
  if (halfopen) {
    --end;
  }

  // turn off sync with stdio to try to increase speed of input/output to standard streams
  std::ios::sync_with_stdio (false);

  // take alignment slice
  Mercator_alignment mercator_alignment (mapdir != "." ? mapdir : datadir, aligndir != "." ? aligndir : datadir);
  Sequence_database seq_db_subalign;
  Stockholm* slice = mercator_alignment.slice (seq_db_subalign,
					       genome, chromosome,
					       strand,
					       start, end,
					       lazy, truncate_ok,
					       stockholm);   // annotate with homology information if stockholm

  // display
  if (stockholm) {
    slice->annotate_with_statistics();
    slice->write_stockholm (cout);
  }
  else
    slice->write_mfa (cout);

  // clean up
  delete slice;

  return 0;
    
}
