
/**
 * \file isect_mercator_alignment_gff.cc
 * This file is part of FSA.
 * \author Source code in this file was written by Robert Bradley.
 */

#include <getopt.h>
#include <fstream>

#include "config.h"
#include "seq/gff.h"
#include "seq/alignment.h"
#include "seq/mercator.h"

using namespace fsa;

static std::string program_name = "isect_mercator_alignment_gff";

static std::string type = "";

static std::string datadir = ".";
static std::string mapdir = ".";
static std::string aligndir = ".";
static bool stockholm = false;
static bool lazy = false;
static bool truncate_ok = false;

static bool verbose = false;

// short options letters
const char* const short_options = "hvt:D:M:A:LUse";

// long options
static struct option long_options[] = {
  {"help", no_argument, NULL, 'h'},
  {"version", no_argument, NULL, 'v'},

  {"type", required_argument, NULL, 't'},

  {"data", required_argument, NULL, 'D'},
  {"map", required_argument, NULL, 'M'},
  {"align", required_argument, NULL, 'A'},
  {"lazy", no_argument, NULL, 'L'},
  {"truncate", no_argument, NULL, 'U'},

  {"stockholm", no_argument, NULL, 's'},

  {"verbose", no_argument, NULL, 'e'},

  /* required NULL termination of array */
  {NULL, 0, NULL, 0}
};

static void print_version (std::ostream& o) {
  o << program_name << " from " << PACKAGE_STRING << endl;
}

static void print_usage (std::ostream& o) {
  print_version (o);
  o << endl
    << "Usage: " << program_name << " [options] <genome> <GFF file>" << endl
    << endl
    << "Extracts subalignments from a Mercator multiple alignment for the features in the GFF file." << endl
    << endl
    << "Options:" << endl
    << "    -h, --help                  show this message" << endl
    << "    -v, --version               show version information" << endl
    << endl
    << "    -t, --type <string>         only look at features of particular type" << endl
    << endl
    << "    -D, --data <directory>      path to map, genome and alignment files" << endl
    << "    -M, --map <directory>       path to map and genome files" << endl
    << "    -A, --align <directory>     path to alignment files" << endl
    << "    -L, --lazy                  warn, rather than die, if the subalignment can't be obtained" << endl
    << "    -U, --truncate              truncate unmappable sequence (rather than skipping) and show truncated subalignment" << endl
    << endl
    << "    -s, --stockholm             use and display Stockholm-format alignments with conservation statistics (default is multi-FASTA)" << endl
    << endl
    << "    -e, --verbose               report progress" << endl
    << endl
    << "PLEASE NOTE: While this program is reasonably fast if the GFF is properly" << endl
    << "ordered by chromosome and the start and end coordinates of features, it" << endl
    << "will be *very slow* if the GFF is not sorted."
    << endl
    << "Assumes that the \"seqid\" or \"name\" field (the first field) of the GFF entries" << endl
    << "holds the chromosome name." << endl
    << endl
    << "Note that the GFF specification defines coordinates to be 1-based" << endl
    << "and fully-closed, therefore representing the interval [start, end]." << endl
    << "Conformance to this specification is assumed internally." << endl
    << endl
    << "If requested, unmappable sequence will be truncated to the mappable portion;" << endl
    << "note that the truncation will favor the beginning of the requested sequence." << endl
    << endl
    << "If a GFF feature is on the - strand, then the corresponding" << endl
    << "subalignment will be reverse-complemented."
    << endl;
}

int main (int argc, char** argv) {

  std::string genome;
  std::string gfffile;

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

    case 'e': // verbose
      verbose = true;
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

  genome = std::string (argv[optind++]);
  gfffile = std::string (argv[optind++]);

  assert (optind == argc);


  // now do stuff!

  // read GFF file
  GFF_database gff_db;
  gff_db.from_file (gfffile,
		    true); // strip_leading_chr

  // create Mercator map
  Mercator_alignment mercator_alignment (mapdir != "." ? mapdir : datadir, aligndir != "." ? aligndir : datadir);

  // initialize persistent sequence & alignment information for current bin
  Sequence_database seq_db_bin;
  Stockholm stockholm_bin (seq_db_bin);
  unsigned bin = 0; // initialize to dummy 0 value

  // turn off sync with stdio to try to increase speed of input/output to standard streams
  std::ios::sync_with_stdio (false);

  // iterate through GFF
  for (std::vector<GFF>::const_iterator gff = gff_db.begin(); gff != gff_db.end(); ++gff) {
    
    // if requested, only look at features of a particular type
    if (type.length() && type != gff->type)
      continue;

    // hold subalignment sequence information
    Sequence_database seq_db_subalign;

    // take slice
    Stockholm* slice = mercator_alignment.slice (seq_db_subalign,
						 bin,
						 stockholm_bin,
						 genome, gff->seqid,
						 gff->strand,
						 gff->start - 1, gff->end - 1, // convert to 0-based coordinates
						 lazy, truncate_ok,
						 stockholm, // annotate with homology information if stockholm
						 verbose);

    // skip empty subalignments
    // (indicating that no mapping was found)
    if (slice->columns() == 0)
      continue;

    // annotate entire alignment with #=GF GFF line
    slice->add_gf_annot (Stockholm::gff_annotation, gff->to_string());

    // display
    if (stockholm) {
      slice->annotate_with_statistics();
      slice->write_stockholm (cout);
    }
    else
      slice->write_mfa (cout);

    // clean up
    delete slice;

  }

  return 0;
    
}
