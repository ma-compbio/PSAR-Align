
/**
 * \file map_coords.cc
 * This file is part of FSA.
 * \author Source code in this file was written by Robert Bradley.
 */

#include <getopt.h>

#include "config.h"
#include "seq/alignment.h"
#include "seq/mercator.h"

using namespace fsa;

static std::string program_name = "map_coords";

static std::string datadir = ".";
static std::string mapdir = ".";
static std::string aligndir = ".";
static bool lazy = false;
static bool truncate_ok = false;

// short options letters
const char* const short_options = "vhD:M:A:LU";

// long options
static struct option long_options[] = {
  {"version", no_argument, NULL, 'v'},
  {"help", no_argument, NULL, 'h'},

  {"data", required_argument, NULL, 'D'},
  {"map", required_argument, NULL, 'M'},
  {"align", required_argument, NULL, 'A'},
  {"lazy", no_argument, NULL, 'L'},
  {"truncate", no_argument, NULL, 'U'},

  /* required NULL termination of array */
  {NULL, 0, NULL, 0}
};

static void print_version (std::ostream& o) {
  o << program_name << " from " << PACKAGE_STRING << endl;
}

static void print_usage (std::ostream& o) {
  print_version (o);
  o << endl
    << "Usage: " << program_name << " [options] <source genome> <chromosome> <start> <end> <strand> <target genome>" << endl
    << endl
    << "Map coordinates from one genome to another using a Mercator multiple alignment." << endl
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
    << "Assumes that coordinates are 1-based and fully-closed," << endl
    << "therefore representing the interval [start, end]." << endl
    << endl
    << "If requested, unmappable sequence will be truncated to the mappable portion;" << endl
    << "note that the truncation will favor the beginning of the requested sequence." << endl
    << endl;
}

int main (int argc, char** argv) {

  std::string genome_source;
  std::string chromosome_source;
  unsigned start_source, end_source;
  char strand_source;
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

    case '?': // invalid option
      print_usage (cerr);
      return 1;

    default:  // unexpected
      abort();

    }
  }

  // stuff
  if (optind + 6 != argc) {
    print_usage (cerr);
    return 1;
  }

  genome_source = std::string (argv[optind++]);
  chromosome_source = std::string (argv[optind++]);
  start_source = static_cast<size_t> (atoi (argv[optind++]));
  end_source = static_cast<size_t> (atoi (argv[optind++]));
  strand_source = argv[optind++][0];
  genome_target = std::string (argv[optind++]);

  assert (optind == argc);


  // now do stuff!

  // create Mercator map
  Mercator_alignment mercator_alignment (mapdir != "." ? mapdir : datadir, aligndir != "." ? aligndir : datadir);

  // initialize persistent sequence & alignment information for current bin
  Sequence_database seq_db_bin;
  Stockholm stockholm_bin (seq_db_bin);
  unsigned bin = 0; // initialize to dummy 0 value

  // use the alignment to find the homologous interval
  Genomic_interval region_target = mercator_alignment.find_homologous_interval (bin,
										stockholm_bin,
										genome_source, chromosome_source,
										start_source - 1, end_source - 1, // convert to 0-based coordinates
										genome_target,
										lazy, truncate_ok);


  // check that it was mapped successfully
  // (Mercator_alignment::find_homologous_interval returns an empty value
  // if it wasn't or if homologous subsequence is empty)
  // NB no need to enforce lazy here, b/c find_homologous_interval does this for us
  if (region_target.genome == "") {
    cerr << "WARNING: No homologous sequence found (all gaps)." << endl;
    return 0;
  }

  // get the Mercator_interval for the source
  const size_t idx_source = mercator_alignment.find_mercator_interval (genome_source, chromosome_source,
								       start_source - 1, end_source - 1, // convert to 0-based coordinates
								       !truncate_ok);
  assert (idx_source < mercator_alignment.size());
  const Mercator_interval& interval_source = mercator_alignment.get_interval (idx_source);

  // use the source interval to determine whether the source feature was 
  // reverse-complemented w.r.t. the strand in the Mercator_interval
  // if it was, then the target feature will also be reverse-complemented
  // w.r.t. the strand in the Mercator_interval
  const bool is_rc = (interval_source.strand != strand_source);

  std::string chromosome_target = region_target.chromosome;
  size_t start_target = region_target.start + 1; // convert back to 1-based coordinates
  size_t end_target = region_target.end + 1;
  char strand_target;
  if (strand_source == GFF::unknown_strand)   // if the source feature had no annotated strand, then do the same for the target feature
      strand_target = GFF::unknown_strand;
  else {                           // otherwise set the target strand appropriately using the homology mapping implied by the alignment
    strand_target = region_target.strand;
    if (is_rc)
      Sequence::complement_strand (strand_target);
  }

  cout << genome_target << '\t' << chromosome_target << '\t'
       << start_target << '\t' << end_target << '\t'
       << strand_target << endl;

  return 0;
    
}
