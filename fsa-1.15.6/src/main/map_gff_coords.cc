
/**
 * \file map_gff_coords.cc
 * This file is part of FSA.
 * \author Source code in this file was written by Robert Bradley.
 */

#include <getopt.h>

#include "config.h"
#include "seq/gff.h"
#include "seq/alignment.h"
#include "seq/mercator.h"

using namespace fsa;

static std::string program_name = "map_gff_coords";

static std::string type = "";

static std::string datadir = ".";
static std::string mapdir = ".";
static std::string aligndir = ".";
static bool lazy = false;
static bool truncate_ok = false;
static bool force_entry = false;

static bool verbose = false;

// short options letters
const char* const short_options = "vht:D:M:A:LUf:e";

// long options
static struct option long_options[] = {
  {"version", no_argument, NULL, 'v'},
  {"help", no_argument, NULL, 'h'},

  {"type", required_argument, NULL, 't'},

  {"data", required_argument, NULL, 'D'},
  {"map", required_argument, NULL, 'M'},
  {"align", required_argument, NULL, 'A'},
  {"lazy", no_argument, NULL, 'L'},
  {"truncate", no_argument, NULL, 'U'},

  {"force-entry", no_argument, NULL, 'f'},

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
    << "Usage: " << program_name << " [options] <source genome> <source GFF file> <target genome>" << endl
    << endl
    << "Map coordinates of GFF features from one genome to another using a Mercator multiple alignment." << endl
    << endl
    << "Options:" << endl
    << "    -h, --help                  show this message" << endl
    << endl
    << "    -t, --type <string>         only look at features of particular type" << endl
    << endl
    << "    -D, --data <directory>      path to map, genome and alignment files" << endl
    << "    -M, --map <directory>       path to map and genome files" << endl
    << "    -A, --align <directory>     path to alignment files" << endl
    << "    -L, --lazy                  warn, rather than die, if the subalignment can't be obtained" << endl
    << "    -U, --truncate              truncate unmappable sequence (rather than skipping) and show truncated subalignment" << endl
    << endl
    << "    -f, --force-entry           if a feature can't be mapped, then add an empty entry to the GFF file (rather than skipping it entirely); implies --lazy" << endl
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
    << "If a GFF feature is on the + strand for the source genome but" << endl
    << "the corresponding homologous region in the target genome is on the - strand," << endl
    << "then the mapped GFF feature will be reported as on the - strand." << endl
    << endl;
}

int main (int argc, char** argv) {

  std::string genome_source;
  std::string gfffile;
  std::string genome_target;

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

    case 'f': // force-entry
      force_entry = true;
      lazy = true;
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
  if (optind + 3 != argc) {
    print_usage (cerr);
    return 1;
  }

  genome_source = std::string (argv[optind++]);
  gfffile = std::string (argv[optind++]);
  genome_target = std::string (argv[optind++]);

  assert (optind == argc);


  // now do stuff!

  // read GFF file
  GFF_database gff_db_from;
  gff_db_from.from_file (gfffile);

  // create Mercator map
  Mercator_alignment mercator_alignment (mapdir != "." ? mapdir : datadir, aligndir != "." ? aligndir : datadir);

  // map the GFF!
  GFF_database gff_db_to = mercator_alignment.map_gff_database (genome_source, genome_target,
								gff_db_from,
								lazy, truncate_ok,
								verbose,
								force_entry);

  gff_db_to.write (cout);

  return 0;
    
}
