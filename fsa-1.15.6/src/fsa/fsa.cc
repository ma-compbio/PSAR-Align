
/**
 * \file fsa.cc
 * This file is part of FSA, a sequence alignment algorithm.
 * \author Source code in this file was written by Robert Bradley.
 * Jaeyoung Do wrote the parallelization and database code.
 */

#include "math/mst.h"
#include "seq/similarity_matrix.h"
#include "fsa/fsa.h"

using namespace fsa;

FSA::FSA (int argc, char** argv) :
  opts (argc, argv),
  alphabet (DNA_alphabet()),
  tree_weights (Tree_weights()), // default constructor sets all weights to 1.0
  w_worker (false)
{ }

void FSA::init_opts() {

  opts.newline();
  INIT_CONSTRUCTED_OPTS_LIST (opts, -1, "[options] <sequence file(s)>", "Distance-based alignment of DNA, RNA and proteins.");

  opts.print_title ("Output options");
  opts.add ("-stockholm", write_stockholm = false, "output Stockholm alignments (default is multi-FASTA format)", false);
  opts.add ("-gui", output_for_gui = false, "record alignment & statistical model for interactive Java GUI", false);
  opts.add ("-write-params", write_params = false, "write learned emission distributions (substitution matrices) to disk", false);
  opts.add ("-write-posteriors", write_dotplots = false, "write learned pairwise posterior alignment probability matrices to disk", false);
  // --write-divergences option disabled for now (9/14/08)
  //  opts.add ("-write-divergences", write_divergences = false, "write pairwise divergences using log-det transform on learned substitution matrices to disk", false);
  write_divergences = false;
  opts.newline();

  opts.print_title ("Parallelization options");
#if defined(HAVE_CONDOR)
  opts.add ("-parallelize", num_parallelized_jobs = 0, "collect pairwise posterior probabilities and candidate edges simultaneously", false);
  opts.add ("-noannealing", noannealing = false, "disable annealing", false);
  opts.print ("      Note that --parallelize is only relevant when database connection options are given;"); opts.newline();
  opts.print ("      similarly, --noannealing is only relevant when --parallelize is set."); opts.newline();
#else
  opts.print ("(Parallelization not available; please see the manual for more information.)");
  opts.newline();
  num_parallelized_jobs = 0;
  noannealing = false;
#endif
  opts.newline();

  opts.print_title ("Database options");
#if defined(HAVE_POSTGRES)
  opts.add ("-db-maxram", db_max_ram = 0, "maximum RAM to use when the database mode (in megabytes)", true); 
  opts.add ("-hostname", db_hostname = "", "database server host name", false);
  opts.add ("-hostaddr", db_hostaddr = "", "database server host IP address", false);
  opts.add ("-dbname", db_name = "", "database name", false);
  opts.add ("-port", db_port = 0, "database server port", false);
  opts.add ("-user", db_user = "", "database username", false);
  opts.add ("-password", db_password = "", "database password", false);
  opts.add ("-noposteriors", read_posteriors_from_db = false, "do not compute posteriors; instead read them from database", false);
  opts.print ("      Note that to access the database server, you must provide at least --hostname and --dbname.");
  opts.newline();
#else
  opts.print ("(Database not available; please see the manual for more information.)");
  opts.newline();
#endif
  opts.newline();

  opts.print_title ("Pair HMM model options");
  opts.add ("-nucprot", nucprot = false, "align input nucleotide sequences (must all be nucleotide) in protein space", false);
  opts.add ("-indel2", is_indel2 = true, "(default) use two sets of indel states in Pair HMM (use --noindel2 to use 1 set only)", false);
  opts.add ("-gapopen1", gap_open1 = 0., "initial gap-open probability (for set 1 of indel states)", false);       // the command-line values, if passed, are used in init_indel()
  opts.add ("-gapextend1", gap_extend1 = 0., "initial gap-extend probability (for set 1 of indel states)", false);
  opts.add ("-gapopen2", gap_open2 = 0., "initial gap-open probability (for set 2 of indel states)", false);
  opts.add ("-gapextend2", gap_extend2 = 0., "initial gap-extend probability (for set 2 of indel states)", false);
  // --raggedends option disabled for now (8/6/08)
  //  opts.add ("-raggedends", ragged_ends = false, "favor padding indels at starts and ends of sequences", true);
  ragged_ends = false;
  opts.add ("-model", model = MODEL_DEFAULT, "initial substitution model: 0 = Jukes-Cantor, 1 = Tamura-Nei / BLOSUM62-like (proteins)", true);
  opts.add ("-time", time = TIME_DEFAULT, "Jukes-Cantor/Tamura-Nei evolutionary time parameter", true);
  opts.add ("-alphar", alpha_R = ALPHA_R_DEFAULT, "Tamura-Nei rate alpha_R (transition: purine)", true);
  opts.add ("-alphay", alpha_Y = ALPHA_Y_DEFAULT, "Tamura-Nei rate alpha_Y (transition: pyrimidine)", true);
  opts.add ("-beta", beta = BETA_DEFAULT, "Tamura-Nei rate beta (transversion)", true);
  opts.add ("-load-probs", load_probs_file = "", "load pairwise posterior probabilities from a file rather than performing inference with Pair HMM", false);
  opts.newline();

  opts.print_title ("Parameter estimation options");
  opts.add ("-learngap", learn_gap, "estimate indel probabilities for each pair of sequences (--nolearngap to disable)", true);
  opts.add ("-learnemit-bypair", learn_emit_bypair, "(default for DNA and RNA) estimate emission probabilities for each pair of sequences (--nolearnemit-bypair to disable)", false);
  opts.add ("-learnemit-all", learn_emit_all, "(default for proteins) estimate emission probabilities averaged over all sequences (--nolearnemit-all to disable)", false);
  opts.add ("-nolearn", nolearn = false, "disable ALL parameter learning (use ProbCons defaults)", false);
  opts.add ("-regularize", regularize = true, "regularize learned emission and gap probabilities with Dirichlet prior", true);
  opts.add ("-regularization-gapscale", regularization_transition_scale, "scaling factor for transition prior", false);
  opts.add ("-regularization-emitscale", regularization_emission_scale, "scaling factor for emission Dirichlet prior", false);
  opts.add ("-mininc", em_min_inc = EM_MIN_INC_DEFAULT, "minimum fractional increase in log-likelihood per round of EM", true);
  opts.add ("-maxrounds", em_max_iter = EM_MAX_ITER_DEFAULT, "maximum number of iterations of EM", true);
  opts.add ("-mingapdata", min_gap_training_data, "minimum amount of sequence data (# of aligned pairs of characters) for training gap probs", false);
  opts.add ("-minemitdata", min_emit_training_data, "minimum amount of sequence data (# of aligned pairs of characters) for training emission probs", false);
  opts.newline();

  opts.print_title ("Multiple alignment options: sequence annealing");
  opts.add ("-refinement", num_refinement_steps = -1, "number of iterative refinement steps (default is unlimited; 0 for none)", false);
  opts.add ("-maxsn", maxsn = false, "maximum sensitivity (instead of highest accuracy)", false);
  // hack to use #define directive
  std::string opt_string_hack = "gap factor; 0 for highest sensitivity (the internal effective minimum is " + Util::to_string (DEFAULT_POSTERIOR_PROB_CUTOFF) + "); >1 for higher specificity";
  opts.add ("-gapfactor", gap_factor = 1, opt_string_hack.c_str(), true);
  opts.add ("-dynamicweights", enable_dynamic_weights = true, "enable dynamic edge re-weighting", true);
  opts.add ("-treeweights", tree_weights_file = "", "weights for sequence pairs based on a tree", false);
  opts.add ("-require-homology", require_homology = false, "require that there be some detectable homology between all input sequences", true);
  opts.newline();

  opts.print_title ("Alignment speedup options: many sequences");
  opts.add ("-fast", fast_defaults = false, "fast alignment: use 5 * Erdos-Renyi threshold percent of sequence pairs for alignment and 2 * for learning", false);
  opts.add ("-refalign", refalign = false, "alignment to a reference sequence only (reference must be first sequence in file)", false);
  opts.add ("-mst-min", num_minimum_spanning_trees = 3, "build --mst-min minimum spanning trees on input sequences for pairwise comparisons", true);
  opts.add ("-mst-max", num_maximum_spanning_trees = 0, "build --mst-max maximum spanning trees on input sequences for pairwise comparisons", true);
  opts.add ("-mst-palm", num_minimum_spanning_palm_trees = 0, "build --mst-palm minimum spanning palm trees on input sequences for pairwise comparisons", true);
  opts.add ("-degree", degree = 0, "use --degree number of pairwise comparisons between closest sequences", true);
  opts.add ("-kmer", kmer_length = 0, "length of k-mers to use when determining sequence similarity", false);
  opts.add ("-alignment-fraction", fraction_alignment_pairs = 1., "randomized fraction of all (n choose 2) pairs of sequences to consider during alignment inference", true);
  opts.add ("-alignment-number", num_alignment_pairs = 0, "total number of (randomized) pairs of sequences to consider during alignment inference", false);
  opts.newline();

  opts.print_title ("Alignment speedup options: long sequences (MUMmer)");
//#ifdef MUMMER_EXEC
  opts.add ("-anchored", anchored = false, "use anchoring (--noanchored to disable)", false);
#ifdef MUMMER_EXEC
  opts.add ("-translated", use_translated = false, "perform anchoring in protein space", true);
  opts.add ("-minlen", anchor_minlen, "minimum length of exact matches for anchoring", false);
  opts.add ("-maxjoinlen", anchor_max_join_length = 2, "maximum ungapped separation of parallel adjacent anchors to join", true);
  // hack to use #define directive
  opt_string_hack.clear(); opt_string_hack = "leave hardmasked sequence >" + Util::to_string (MIN_HARDMASK_LENGTH) + " nt unaligned instead of randomizing it (default for long DNA)";
  opts.add ("-hardmasked", hardmasked = false, opt_string_hack.c_str(), false);
#else
  anchored = false;
  use_translated = false;
  hardmasked = false;
  opts.print ("(MUMmer not available; please see the manual for more information.)");
  opts.newline();
#endif
  opts.newline();

  opts.print_title ("Alignment speedup options: long sequences (exonerate)");
#ifdef EXONERATE_EXEC
  opts.add ("-exonerate", exonerate = false, "call exonerate to get anchors (implies --anchored)", false);
  opts.add ("-minscore", exonerate_minscore = EXONERATE_MINSCORE_DEFAULT, "minimum score of alignments found by exonerate", true);
  opts.add ("-softmasked", softmasked = false, "input sequences are softmasked", false);
#else
  exonerate = false;
  softmasked = false;
  opts.print ("(exonerate not available; please see the manual for more information.)");
  opts.newline();
#endif
  opts.newline();

  opts.print_title ("Alignment speedup options: long sequences (Mercator)");
  opts.add ("-mercator", mercator_constraint_file = "", "input Mercator constraints", false);
  opts.newline();

  opts.print_title ("Memory savings");
  opts.add ("-maxram", max_ram = (total_ram() > 0 ? static_cast<int> (0.8 * total_ram()) : -1), "maximum RAM to use (in megabytes)", true); // use 80% of max RAM by default
  opts.add ("-bandwidth", bandwidth = 0, "banding (default is no banding)", false);
  opts.add ("-minprob", minprob = DEFAULT_POSTERIOR_PROB_CUTOFF, "minimum posterior probability to store", true);
  opts.newline();

  opts.newline();
  opts.print ("Input sequence file(s) must be in FASTA format."); opts.newline();
  opts.newline();
  opts.print ("FSA attempts to automatically figure out appropriate settings;"); opts.newline();
  opts.print ("you can override its automated choices with the above options."); opts.newline();
  opts.newline();

  opts.print ("Please contact the FSA team at fsa@math.berkeley.edu with any questions or comments.");

}

void FSA::input_data() {

  // get sequence filenames
  const std::vector<sstring>& seq_filename = opts.get_not_opts();

  // check for at least one sequence file
  if (seq_filename.size() == 0) {
    cerr << opts.short_help() << endl
	 << "ERROR: Please specify at least one sequence file." << endl;
    exit (1);
  }
  // read all files into a single Sequence_database
  for (std::vector<sstring>::const_iterator sf = seq_filename.begin(); sf != seq_filename.end(); ++sf)
    seq_db.read_fasta (*sf, &Alignment::is_gap_char,
		       false, // strip_leading_chr
		       true); // verbose

  // calculate some statistics
  num_seqs = seq_db.size();
  num_seq_pairs = num_seqs * (num_seqs - 1) / 2;
  if (num_seqs < 2)
    THROWEXPR (opts.short_help() << "You must include at least 2 sequences to be aligned.");
  CTAG(9,FSA) << "Read " << num_seqs << " sequences." << endl;

  // figure out which alphabet we're using
  if (seq_db.matches_alphabet (DNA_alphabet()) || seq_db.matches_alphabet (RNA_alphabet())) {
    is_dna = true;
    is_protein = false;
    alphabet = DNA_alphabet();
    // store alphabet as string
    // NB can't use e.g. Alphabet::chars_uc because not alphabetical sorting (which we require).
    alphabet_string = DNA_ALPHABET_STRING;
  }
  else if (seq_db.matches_alphabet (Protein_alphabet())) {
    is_dna = false;
    is_protein = true;
    alphabet = Protein_alphabet();
    alphabet_string = PROTEIN_ALPHABET_STRING;
  }
  else {
    THROWEXPR ("ERROR: This doesn't seem to be nucleotide or amino acid sequence; I'm bailing.");
  }

  // get gui filename
  gui_prefix = *(seq_filename.begin());

}

void FSA::set_up_defaults() {

  // set regularization scales to defaults
  regularization_transition_scale = is_dna ? REGULARIZATION_NUC_TRANSITION_SCALE_DEFAULT : REGULARIZATION_AA_TRANSITION_SCALE_DEFAULT;
  regularization_emission_scale = is_dna ? REGULARIZATION_NUC_EMISSION_SCALE_DEFAULT : REGULARIZATION_AA_EMISSION_SCALE_DEFAULT;

  // set training data cutoffs to defaults
  min_emit_training_data = is_dna ? MIN_NUC_EMISSION_TRAINING_DATA : MIN_AA_EMISSION_TRAINING_DATA;
  min_gap_training_data = is_dna ? MIN_NUC_TRANSITION_TRAINING_DATA : MIN_AA_TRANSITION_TRAINING_DATA;

  // set up anchoring defaults
  anchor_minlen = is_dna ? ANCHOR_NUC_MINLEN_DEFAULT : ANCHOR_AA_MINLEN_DEFAULT;

  const size_t meanlen = static_cast<size_t> (std::floor (seq_db.meanlength()));

  // model
  const bool dna_defaults = (is_dna && (meanlen <= LONG_DNA_DEFAULT)) ? true : false;
  const bool longdna_defaults = (is_dna && (meanlen > LONG_DNA_DEFAULT)) ? true : false;
  const bool protein_defaults = (is_protein) ? true : false;

  // It's best to set the minimal number of options here:
  // Only set things which are /different/ from the initialization
  // in init_opts.
  //  - explicitly turn on anchoring if desired, but don't turn it off
  // For example, I used to explicitly set
  //  model = MODEL_DEFAULT
  // but it's best to just leave this job to init_opts.

  // small DNA
  if (dna_defaults) {
    learn_gap = true;
    learn_emit_bypair = true;
    learn_emit_all = false;
    hardmasked = true;
  }

  // same as dna_defaults, but with anchoring
  if (longdna_defaults) {
    learn_gap = true;
    learn_emit_bypair = true;
    learn_emit_all = false;
    anchored = true;
    hardmasked = true;
  }

  // protein:
  if (protein_defaults) {
    learn_gap = (num_seq_pairs * meanlen > static_cast<size_t> (min_gap_training_data)) ? true : false;
    learn_emit_bypair = false;
    learn_emit_all = (num_seq_pairs * meanlen > static_cast<size_t> (min_emit_training_data)) ? true : false;
  }

  // warning about anchoring (MUMmer)
#ifndef MUMMER_EXEC
  if (meanlen > 5000) {
    CTAG(9,FSA) << endl << "WARNING: MUMmer is not available, so we are running in unanchored mode; inference may be slow or even fail." << endl << endl;
    anchored = false;
  }
#endif

  // set up default word length to use for sequence similarity comparisons
  kmer_length = Sequence_kmer_counts::choose_minimum_word_length (seq_db,
								  alphabet);

}

void FSA::parse_opts() {

  // parse command-line options
  opts.parse_or_die();


  /***********************************************
   * Output options
   ***********************************************/

  // no point in writing divergences unless learning by pair
  if (!learn_emit_bypair && write_divergences)
    THROWEXPR ("ERROR: Writing divergences doesn't make sense unless --learnemit-bypair is set.");


  /***********************************************
   * Parallelization options
   ***********************************************/

  // ensure valid number of parallelized jobs
  if ((static_cast<size_t> (num_parallelized_jobs) > num_seq_pairs) || (static_cast<size_t> (num_parallelized_jobs) > MAX_NUM_PARALLELIZED_JOBS)) {
    CTAG(9,FSA) << "WARNING: Reducing --parallelize " << num_parallelized_jobs << " to ";
    if (MAX_NUM_PARALLELIZED_JOBS < num_seq_pairs) {
      CL << MAX_NUM_PARALLELIZED_JOBS <<", the maximum number of parallelized jobs" << endl;
      num_parallelized_jobs = MAX_NUM_PARALLELIZED_JOBS;
    } else {
      CTAG(9,FSA) << num_seq_pairs  <<", the maximum possible number of unique sequence pairs." << endl;
      num_parallelized_jobs = num_seq_pairs;
    }
  }
  else if (num_parallelized_jobs != 0 && num_parallelized_jobs < MIN_NUM_PARALLELIZED_JOBS) {
    CTAG(9,FSA) << "WARNING: Increasing --parallelize " << num_parallelized_jobs << " to " << MIN_NUM_PARALLELIZED_JOBS <<", the minimum number of parallelized jobs" << endl;
    num_parallelized_jobs = MIN_NUM_PARALLELIZED_JOBS;
  }
  parallelizing = (num_parallelized_jobs > 0) ? true : false;

  /***********************************************
   * Database options
   ***********************************************/

  write_db = !read_posteriors_from_db && (db_name.length() > 0);

  
  /***********************************************
   * Pair HMM model options
   ***********************************************/

  // handle nucprot option
  if (nucprot) {

    // sanity checks
    if (!is_dna)
      THROWEXPR ("ERROR: Input sequences must be nucleotide in order to use --nucprot option.");
    if (output_for_gui)
      THROWEXPR ("ERROR: You currently cannot view GUI output for --nucprot alignments.  Sorry!");
    if (mercator_constraint_file != "")
      THROWEXPR ("ERROR: You cannot specify Mercator constraints for --nucprot alignments.");
    if (hardmasked) {
      CTAG(9,FSA) << "WARNING: Hardmasking is not allowed for --nucprot alignments; I'm disabling hardmasking." << endl;
      hardmasked = false;
    }
    if (anchored) {
      CTAG(9,FSA) << "WARNING: Anchoring is not allowed for --nucprot alignments; I'm disabling anchoring." << endl;
      anchored = false;
    }

    // fix the alphabet as appropriate
    is_dna = false;
    is_protein = true;
    alphabet = Protein_alphabet();
    alphabet_string = PROTEIN_ALPHABET_STRING;

    // warn if learning strategy isn't acceptable.
    if (learn_emit_bypair) {
      CTAG(9,FSA) << "WARNING: You cannot use --learnemit-bypair with --nucprot; I'm disabling --learnemit-bypair and enabling --learnemit-all." << endl;
      learn_emit_bypair = false;
    }

    // now set default options all over again (but for proteins)
    // NB this prevent users from changing some options on the command line
    set_up_defaults();

  }

  // sane parameters to Tamura-Nei
  if ((time < 0) || (alpha_R < 0) || (alpha_Y < 0) || (beta < 0))
    THROWEXPR ("ERROR: Jukes-Cantor and Tamura-Nei parameters --time, --alphar, --alphay and --beta must be positive.");


  /***********************************************
   * Parameter estimation options
   ***********************************************/

  // don't allow too many iterations of EM
  if (em_max_iter > EM_MAX_ITER_DEFAULT) {
    CL << "Reducing max rounds of EM to " << EM_MAX_ITER_DEFAULT << endl;
    em_max_iter = EM_MAX_ITER_DEFAULT;
  }

  // enforce no learning if so requested
  // initialize parameters to defaults which don't depend on the data:
  // ProbCons defaults for RNA and BLOSUM62 for protein
  if (nolearn) {
    learn_emit_all = learn_emit_bypair = learn_gap = false;
    model = 5;
  }

  // check learning strategy is consistent
  if (learn_emit_all && learn_emit_bypair)
    THROWEXPR ("ERROR: You cannot learn both across all sequences and by pairs of sequences.  Please tell me which to ignore with --nolearnemit-all or --nolearnemit-bypair.");

  // check that regularization is sane
  if ((regularization_transition_scale < 0) || (regularization_emission_scale < 0))
    THROWEXPR ("ERROR: Regularization scales must be >= 0.");

  // can't learn by all if anchored
  if (learn_emit_all && anchored)
    THROWEXPR ("ERROR: Only --learnemit-bypair is allowed for anchoring; try using --nolearnemit-all or --nolearn.");


  /***********************************************
   * Multiple alignment options: sequence annealing
   ***********************************************/

  // refinement
  if (num_refinement_steps == -1)
    num_refinement_steps = 99999999;

  // warning about low gap factor
  if (gap_factor < minprob)
    CL << "A gap factor < " << minprob << " is not meaningful; " << minprob << " will be used internally instead." << endl;

  // maximum sensitivity option
  if (maxsn)
    gap_factor = 0;
  
  // use tgf weighting always
  // (originally the maxsn option triggered maxstep weighting,
  // but empirically I found that tgf weighting gave essentially identical performance
  // + is truly steepest-ascent => use tgf weighting, even for maxsn mode)
  use_tgf = true;


  /***********************************************
   * Alignment speedup options: many sequences
   ***********************************************/

  // recommendation to use --fast
  if (!fast_defaults && (num_seqs > 100))
    CTAG(9,FSA) << endl << "WARNING: It is highly recommended that you invoke the --fast option when aligning many sequences.  If you do not, then inference may take a very long time or even fail." << endl << endl;

  // k-mer length must be > 0
  // (catch the case of k = 0, which occurs automatically if the median sequence length is 0)
  if (kmer_length <= 0) {
    if (seq_db.median_length() != 0)
      THROWEXPR ("ERROR: k-mer word length must be > 0.");
  }

  // degree must be >= 0
  if (degree < 0)
    THROWEXPR ("ERROR: The requested --degree setting must be in 0, ..., (N - 1).");
  // and less than N
  if (degree > static_cast<int> (seq_db.size() - 1))
    degree = seq_db.size() - 1;

  // can't have both refalign and fast
  if (refalign) {
    if (fast_defaults)
      THROWEXPR ("ERROR: You can't select both --fast and --refalign; please choose one or the other.");
    num_minimum_spanning_trees = 0;
    num_maximum_spanning_trees = 0;
    num_minimum_spanning_palm_trees = 0;
    degree = 0;
  }

  // fast options:
  //  Scaled by 5 (so we do ~90% of all (10 choose 2) pairs for an alignment of 15 sequences)
  if (fast_defaults)
    fraction_alignment_pairs = 5 * Sequence_pair_selector::erdos_renyi_p_cutoff (num_seqs);

  // if alignment to a reference, use minimal number of sequence pairs
  if (refalign)
    num_alignment_pairs = num_seqs - 1;

  // convert from fraction to absolute number of pairs for inference
  // if num_alignment_pairs not specified, then initialize with 
  //  fraction_alignment_pairs if specified
  //  all pairs if not
  if (num_alignment_pairs < DOUBLE_TINY)
    num_alignment_pairs = static_cast<int> (std::floor (num_seq_pairs * fraction_alignment_pairs));

  // make sure that we can generate a spanning tree to get a complete alignment
  if (static_cast<size_t> (num_alignment_pairs) < (num_seqs - 1)) {
    CTAG(9,FSA) << "WARNING: Increasing --alignment-number to (N - 1) = " << num_seqs - 1 << " to ensure a complete alignment." << endl;
    num_alignment_pairs = num_seqs - 1;
  }

  // make sure that a reasonable number of alignment pairs are requested
  if (static_cast<size_t> (num_alignment_pairs) > num_seq_pairs) {
    // only warn if --fast wasn't invoked
    if (!fast_defaults)
      CTAG(9,FSA) << "WARNING: Reducing --alignment-number " << num_alignment_pairs << " to " << num_seq_pairs << ", the maximum possible number of unique sequence pairs." << endl;
    num_alignment_pairs = num_seq_pairs;
  }

  // now store fractions (used for logging)
  fraction_alignment_pairs = static_cast<double> (num_alignment_pairs) / num_seq_pairs;



  /***********************************************
   * Alignment speedup options: long sequences (MUMmer)
   ***********************************************/

  // Exonerate implies anchoring
  if (exonerate)
    anchored = true;

  // check sane
#ifndef MUMMER_EXEC
  if (anchored)
    THROWEXPR ("ERROR: FSA was compiled without MUMmer support.  Please recompile with MUMmer if you wish to use anchoring.");
#endif

  // don't try to translate sequence if we're already in protein space
  if (is_protein && use_translated)
    use_translated = false;

  // translated implies anchoring
  if (use_translated)
    anchored = true;

  // set up proper minlen if doing translated anchoring
  if (use_translated)
    anchor_minlen = ANCHOR_AA_MINLEN_DEFAULT;
  
  // no hardmasking for protein sequence
  if (!is_dna && hardmasked)
    THROWEXPR ("ERROR: Hardmasking is only available for nucleotide sequence.");

  /***********************************************
   * Alignment speedup options: long sequences (exonerate)
   ***********************************************/

  // check sane
#ifndef EXONERATE_EXEC
  if (exonerate)
    THROWEXPR ("ERROR: FSA was compiled without exonerate support.  Please recompile with exonerate if you wish to use exonerate for anchoring.");
#endif

  // recommend softmasking
  if (exonerate && !softmasked)
    CTAG(9,FSA) << endl << "WARNING: It is HIGHLY RECOMMENDED that you softmask sequence and use the --softmasked option when calling exonerate." << endl << endl;

  // warning about anchoring (exonerate)
  if ((seq_db.meanlength() > 50000) && !exonerate) {
    CTAG(9,FSA) << endl << "WARNING: Unless you are aligning well-conserved sequences, it is highly recommend that you invoke the exonerate program (as well as MUMmer) when aligning very long sequences." << endl << endl;
  }

  /***********************************************
   * Alignment speedup options: long sequences (Mercator)
   ***********************************************/


  /***********************************************
   * Memory savings
   ***********************************************/

  // check that bandwidth is sane
  if (bandwidth < 0)
    THROWEXPR ("ERROR: Please specify a positive banding width.");

}

void FSA::init_indel_params (Params& params) {

  params.is_indel2 = is_indel2;
  params.to_end = DOUBLE_TINY;

  // if nucleotide
  if (is_dna) {
    if (!is_indel2) {
      params.gap_open1 = 0.02;
      params.gap_extend1 = 0.3;
      params.gap_open2 = 0;
      params.gap_extend2 = 0;
    } else {
      params.gap_open1 = 0.012; // have mixture components still sum to gap_open prob of 0.02 as above
      params.gap_extend1 = 0.4;
      params.gap_open2 = 0.008;
      params.gap_extend2 = 0.9;
    }
  }

  // if amino acid (ProbCons defaults for these -- taken from probcons/Defaults.h)
  else if (is_protein) {
    if (!is_indel2) {
      params.gap_open1 = 0.01993141696f;
      params.gap_extend1 = 0.7943345308f;
      params.gap_open2 = 0;
      params.gap_extend2 = 0;
    } else {
      params.gap_open1 = 0.0119511066f;
      params.gap_extend1 = 0.3965826333f;
      params.gap_open2 = 0.008008334786f;
      params.gap_extend2 = 0.8988758326f;
    }
  }

  // unreachable code
  else {
    THROWEXPR ("ERROR: These don't seem to be nucleotide or amino acid sequences!  This message should never occur.");
  }

  // initialize indel params to command-line values if passed
  if (gap_open1 > DOUBLE_TINY)
    params.gap_open1 = gap_open1;
  if (gap_extend1 > DOUBLE_TINY)
    params.gap_extend1 = gap_extend1;
  if (gap_open2 > DOUBLE_TINY)
    params.gap_open2 = gap_open2;
  if (gap_extend2 > DOUBLE_TINY)
    params.gap_extend2 = gap_extend2;

  // log seed params
  if (CTAGGING(-1,FSAEM)) {
    CL << "Seed indel parameters:" << endl;
    params.show_transition_matrix (CL);
  }

}


void FSA::init_subst_params (Params& params, const std::vector<double>& char_freq, const unsigned model /* = MODEL_DEFAULT */, const double time /* = TIME_DEFAULT */) {

  // time parameter
  params.time = time;

  // alphabet_string
  params.alphabet_string = alphabet_string;

  // allocate memory
  params.single_dist.resize (alphabet_string.length());
  params.pair_dist.resize (alphabet_string.length());
  for (size_t i = 0; i < alphabet_string.length(); i++)
    params.pair_dist[i].resize (alphabet_string.length());

  // init nuc models
  if (is_dna) {

    // ProbConsRNA defaults (taken from probconsRNA/Defaults.h)
    const double probcons_single[4] = {
      0.2270790040f, 0.2422080040f, 0.2839320004f, 0.2464679927f };
    const double probcons_pair[4][4] = {
      { 0.1487240046f, 0.0184142999f, 0.0361397006f, 0.0238473993f },
      { 0.0184142999f, 0.1583919972f, 0.0275536999f, 0.0389291011f },
      { 0.0361397006f, 0.0275536999f, 0.1979320049f, 0.0244289003f },
      { 0.0238473993f, 0.0389291011f, 0.0244289003f, 0.1557479948f } };


    // now set up the substitution model
    switch (model) {

    case 0: // Jukes-Cantor model
      for (size_t i = 0; i < 4; i++) {
	params.single_dist[i] = 1 / 4.0;
	for (size_t j = 0; j < 4; j++) {
	  if (i == j) {
	    params.pair_dist[i][j] = params.single_dist[i] * (1/4.0 + (3/4.0) * std::exp (-4.0 * time / 3.0));
	  } else {
	    params.pair_dist[i][j] = params.single_dist[i] * (1/4.0 - (1/4.0) * std::exp (-4.0 * time / 3.0));
	  }
	}
      }
      break;

    case 1: // Tamura-Nei model
      // empirical char freqs
      params.single_dist.assign (char_freq.begin(), char_freq.end());
      for (size_t i = 0; i < 4; i++) {
	for (size_t j = 0; j < 4; j++) {
	  const double alphai = (i == 0 || i == 2) ? alpha_R : alpha_Y;
	  const double pj = (j == 0 || j == 2) ? (char_freq[0] + char_freq[2]) : (char_freq[1] + char_freq[3]); // piR or piY, as appropriate for j
	  assert (alphai > 0 && pj > 0);
	  // if no change
	  if (i == j) {
	    params.pair_dist[i][j] = std::exp (-(alphai + beta) * time)
	      + std::exp (-beta * time) * (1 - std::exp (-alphai * time)) * (char_freq[j] / pj)
	      + (1 - std::exp (-beta * time)) * char_freq[j];
	    params.pair_dist[i][j] *= char_freq[i]; // convert from subst. prob. to joint prob.
	  }
	  // else if transition
	  //	  else if (std::abs (static_cast<int> (static_cast<int> (i) - j)) == 2) {
	  else if (abs (static_cast<int> (static_cast<int> (i) - j)) == 2) {
	    params.pair_dist[i][j] = 0.
	      + std::exp (-beta * time) * (1 - std::exp (-alphai * time)) * (char_freq[j] / pj)
	      + (1 - std::exp (-beta * time)) * char_freq[j];
	    params.pair_dist[i][j] *= char_freq[i];
	  }
	  // else transversion
	  else {
	    params.pair_dist[i][j] = 0.
	      + 0.
	      + (1 - std::exp (-beta * time)) * char_freq[j];
	    params.pair_dist[i][j] *= char_freq[i];
	  }
	}
      }
      break;

    case 5: // ProbCons
      for (size_t i = 0; i < 4; i++) {
	params.single_dist[i] = probcons_single[i];
	for (size_t j = 0; j < 4; j++) {
	  params.pair_dist[i][j] = probcons_pair[i][j];
	}
      }
      break;

    default:
      THROWEXPR ("ERROR: Model " << model << " does not exist.  Please choose a valid model." << endl);
    }

  }

  // proteins
  else {

    // BLOSUM62 matrix (taken from probcons/Defaults.h and converted per my AA sorting)
    double BLOSUM62_single[20] = {
      0.07831005, 0.02189704, 0.05130349, 0.05615771, 0.04463228, 0.07783433, 0.02601093, 0.06511648, 0.05877077, 0.09716489, 0.02438117, 0.04433257, 0.03940142, 0.03585766, 0.05246024, 0.05849916, 0.05115306, 0.07343426, 0.01203523, 0.03124726 };
    double BLOSUM62_pair[20][20] = {
      { 0.02373072, 0.00145515, 0.00223549, 0.00332218, 0.00165004, 0.00597898, 0.00114353, 0.00318853, 0.00331693, 0.00449576, 0.00148878, 0.00210228, 0.00230618, 0.00219102, 0.00244502, 0.00631752, 0.00389995, 0.00533241, 0.00039119, 0.0013184 },
      { 0.00145515, 0.0101347, 0.00036798, 0.00037956, 0.00052274, 0.00071206, 0.00026421, 0.0009404, 0.00046951, 0.00138494, 0.00037421, 0.00042479, 0.00034766, 0.00032102, 0.00044701, 0.00094867, 0.00073798, 0.00119152, 0.00010666, 0.00036626 },
      { 0.00223549, 0.00036798, 0.01911178, 0.004968, 0.00069041, 0.00235249, 0.00097077, 0.00105355, 0.00252518, 0.00161966, 0.00047808, 0.0035354, 0.00125381, 0.00176784, 0.00161657, 0.00285226, 0.00180488, 0.00127915, 0.00016015, 0.00066005 },
      { 0.00332218, 0.00037956, 0.004968, 0.01676565, 0.00078814, 0.0021486, 0.00131767, 0.00124207, 0.0042842, 0.00222063, 0.00076105, 0.00224738, 0.0015155, 0.00345128, 0.00268865, 0.00293898, 0.0021676, 0.00178697, 0.00023815, 0.00092548 },
      { 0.00165004, 0.00052274, 0.00069041, 0.00078814, 0.01661038, 0.00115204, 0.00072545, 0.00279948, 0.00087222, 0.00533369, 0.00116111, 0.00084658, 0.00060701, 0.00059248, 0.00090768, 0.00119036, 0.00107595, 0.00256159, 0.00085751, 0.00368739 },
      { 0.00597898, 0.00071206, 0.00235249, 0.0021486, 0.00115204, 0.04062876, 0.00103704, 0.0014252, 0.00259311, 0.00212853, 0.00066504, 0.00288882, 0.00155601, 0.00142432, 0.00194865, 0.00381962, 0.00214841, 0.00194579, 0.00038786, 0.00089301 },
      { 0.00114353, 0.00026421, 0.00097077, 0.00131767, 0.00072545, 0.00103704, 0.00867996, 0.00059716, 0.00121376, 0.00111754, 0.00042237, 0.00141205, 0.00049078, 0.00113901, 0.00132105, 0.00116422, 0.00077747, 0.00071553, 0.00019097, 0.00131038 },
      { 0.00318853, 0.0009404, 0.00105355, 0.00124207, 0.00279948, 0.0014252, 0.00059716, 0.01778263, 0.00157852, 0.01071834, 0.00224097, 0.00104273, 0.00103767, 0.00100883, 0.00138145, 0.00173565, 0.00248968, 0.01117956, 0.00039549, 0.00127857 },
      { 0.00331693, 0.00046951, 0.00252518, 0.0042842, 0.00087222, 0.00259311, 0.00121376, 0.00157852, 0.01612228, 0.00259626, 0.0009612, 0.0025731, 0.00154836, 0.00312308, 0.0059565, 0.00312633, 0.00250862, 0.00210897, 0.00028448, 0.00100817 },
      { 0.00449576, 0.00138494, 0.00161966, 0.00222063, 0.00533369, 0.00212853, 0.00111754, 0.01071834, 0.00259626, 0.03583921, 0.00461939, 0.00160275, 0.0015731, 0.00180553, 0.00246811, 0.00250962, 0.00302273, 0.0091446, 0.00076736, 0.00219713 },
      { 0.00148878, 0.00037421, 0.00047808, 0.00076105, 0.00116111, 0.00066504, 0.00042237, 0.00224097, 0.0009612, 0.00461939, 0.00409522, 0.00063401, 0.00046718, 0.00075546, 0.00076734, 0.00087787, 0.00093371, 0.00197461, 0.00016253, 0.00054105 },
      { 0.00210228, 0.00042479, 0.0035354, 0.00224738, 0.00084658, 0.00288882, 0.00141205, 0.00104273, 0.0025731, 0.00160275, 0.00063401, 0.01281864, 0.00100282, 0.00158223, 0.00207782, 0.00301397, 0.00220144, 0.00136609, 0.00021006, 0.0007496 },
      { 0.00230618, 0.00034766, 0.00125381, 0.0015155, 0.00060701, 0.00155601, 0.00049078, 0.00103767, 0.00154836, 0.0015731, 0.00046718, 0.00100282, 0.01846071, 0.00090111, 0.00106268, 0.00180037, 0.00147982, 0.00135781, 0.00015674, 0.00047608 },
      { 0.00219102, 0.00032102, 0.00176784, 0.00345128, 0.00059248, 0.00142432, 0.00113901, 0.00100883, 0.00312308, 0.00180553, 0.00075546, 0.00158223, 0.00090111, 0.00756604, 0.00253532, 0.00191155, 0.00154526, 0.00132844, 0.00020592, 0.00070192 },
      { 0.00244502, 0.00044701, 0.00161657, 0.00268865, 0.00090768, 0.00194865, 0.00132105, 0.00138145, 0.0059565, 0.00246811, 0.00076734, 0.00207782, 0.00106268, 0.00253532, 0.01775118, 0.0022454, 0.00186053, 0.00169359, 0.00029139, 0.0009943 },
      { 0.00631752, 0.00094867, 0.00285226, 0.00293898, 0.00119036, 0.00381962, 0.00116422, 0.00173565, 0.00312633, 0.00250962, 0.00087787, 0.00301397, 0.00180037, 0.00191155, 0.0022454, 0.01346609, 0.00487295, 0.00241601, 0.00026525, 0.00102648 },
      { 0.00389995, 0.00073798, 0.00180488, 0.0021676, 0.00107595, 0.00214841, 0.00077747, 0.00248968, 0.00250862, 0.00302273, 0.00093371, 0.00220144, 0.00147982, 0.00154526, 0.00186053, 0.00487295, 0.01299436, 0.00343452, 0.00024961, 0.00094759 },
      { 0.00533241, 0.00119152, 0.00127915, 0.00178697, 0.00256159, 0.00194579, 0.00071553, 0.01117956, 0.00210897, 0.0091446, 0.00197461, 0.00136609, 0.00135781, 0.00132844, 0.00169359, 0.00241601, 0.00343452, 0.02075171, 0.00038538, 0.00148001 },
      { 0.00039119, 0.00010666, 0.00016015, 0.00023815, 0.00085751, 0.00038786, 0.00019097, 0.00039549, 0.00028448, 0.00076736, 0.00016253, 0.00021006, 0.00015674, 0.00020592, 0.00029139, 0.00026525, 0.00024961, 0.00038538, 0.00563625, 0.00069226 },
      { 0.0013184, 0.00036626, 0.00066005, 0.00092548, 0.00368739, 0.00089301, 0.00131038, 0.00127857, 0.00100817, 0.00219713, 0.00054105, 0.0007496, 0.00047608, 0.00070192, 0.0009943, 0.00102648, 0.00094759, 0.00148001, 0.00069226, 0.00999315 } };


    switch (model) {

    case 1: // BLOSUM62, but converted to empirical target aa frequencies
      params.single_dist.assign (char_freq.begin(), char_freq.end());
      for (size_t i = 0; i < 20; i++) {
	for (size_t j = 0; j < 20; j++) {
	  params.pair_dist[i][j] = BLOSUM62_pair[i][j] * (char_freq[i] * char_freq[j]) / (BLOSUM62_single[i] * BLOSUM62_single[j]);
	}
      }
      params.normalize();
      break;

    case 5: // BLOSUM62
      for (size_t i = 0; i < 20; i++) {
	params.single_dist[i] = BLOSUM62_single[i];
	for (size_t j = 0; j < 20; j++) {
	  params.pair_dist[i][j] = BLOSUM62_pair[i][j];
	}
      }
      break;

    default:
      THROWEXPR ("ERROR: Model " << model << " does not exist.  Please choose a valid model." << endl);
    }

  }

  params.assert_normalized();

  // log seed params
  if (CTAGGING(-1,FSAEM)) {
    CL << "Seed emit parameters:" << endl;
    params.show_emission (CL);
  }

  params.assert_normalized();

}

void FSA::init_pseudocounts (Params& pseudo, const Params& seed, double transition_scale, double emission_scale) {

  if (CTAGGING(2,FSA))
    CL << "Initializing pseudocounts." << endl;

  // get seed distributions
  pseudo.copy_all (seed);

  // now scale to get counts
  pseudo.gap_open1 *= transition_scale;
  pseudo.gap_open2 *= transition_scale;
  pseudo.gap_extend1 *= transition_scale;
  pseudo.gap_extend2 *= transition_scale;
  for (size_t i = 0; i < pseudo.pair_dist.size(); ++i) {
    pseudo.single_dist[i] *= emission_scale;
    for (size_t j = 0; j < pseudo.pair_dist.size(); ++j)
      pseudo.pair_dist[i][j] *= emission_scale;
  }

}

bool FSA::train_params (Params& params, const Sequence_database& seq_db_train, const Sequence_pairs& subset,
			bool left_match, bool right_match, bool ragged_ends,
			const Params& pseudocounts, bool learn_emit) {

  // enforce max_ram cutoff
  for (Sequence_pairs::const_iterator pair = subset.begin(); pair != subset.end(); ++pair) {
    const size_t i = pair->first;
    const size_t j = pair->second;
    const size_t needed = estimate_ram_needed (seq_db_train.get_seq (i).seq.length(), seq_db_train.get_seq (j).seq.length());
    if ((max_ram > 0) && (static_cast<int> (needed) > max_ram)) {
      CTAG(9,FSAEM) << "WARNING: Skipping inference step because estimated memory usage (" << needed << " MB) "
		    << "is greater than specified limit (" << max_ram << " MB)." << endl
		    << "Offending sequences are '" << seq_db_train.get_seq (i).name << "' and '" << seq_db_train.get_seq (j).name << "'." << endl;
      return false;
    }
  }

  // without banding
  if (bandwidth == 0) {
    if (is_dna && !is_indel2)
      return Model::train_params_engine<NucleotideWithoutBanding> (params, seq_db_train, subset,
								   left_match, right_match, ragged_ends,
								   learn_gap, learn_emit, regularize, pseudocounts, em_max_iter, em_min_inc);
    else if (is_dna && is_indel2)
      return Model::train_params_engine<NucleotideIndel2WithoutBanding> (params, seq_db_train, subset,
									 left_match, right_match, ragged_ends,
									 learn_gap, learn_emit, regularize, pseudocounts, em_max_iter, em_min_inc);
    else if (!is_dna && !is_indel2)
      return Model::train_params_engine<AminoAcidWithoutBanding> (params, seq_db_train, subset,
								  left_match, right_match, ragged_ends,
								  learn_gap, learn_emit, regularize, pseudocounts, em_max_iter, em_min_inc);
    else if (!is_dna && is_indel2)
      return Model::train_params_engine<AminoAcidIndel2WithoutBanding> (params, seq_db_train, subset,
									left_match, right_match, ragged_ends,
									learn_gap, learn_emit, regularize, pseudocounts, em_max_iter, em_min_inc);
    else
      THROWEXPR ("ERROR: Unreachable code.");
  }
  // with banding
  else {
    if (is_dna && !is_indel2)
      return Model::train_params_engine<NucleotideWithBanding> (params, seq_db_train, subset,
								left_match, right_match, ragged_ends,
								learn_gap, learn_emit, regularize, pseudocounts, em_max_iter, em_min_inc);
    else if (is_dna && is_indel2)
      return Model::train_params_engine<NucleotideIndel2WithBanding> (params, seq_db_train, subset,
								      left_match, right_match, ragged_ends,
								      learn_gap, learn_emit, regularize, pseudocounts, em_max_iter, em_min_inc);
    else if (!is_dna && !is_indel2)
      return Model::train_params_engine<AminoAcidWithBanding> (params, seq_db_train, subset,
							       left_match, right_match, ragged_ends,
							       learn_gap, learn_emit, regularize, pseudocounts, em_max_iter, em_min_inc);
    else if (!is_dna && is_indel2)
      return Model::train_params_engine<AminoAcidIndel2WithBanding> (params, seq_db_train, subset,
								     left_match, right_match, ragged_ends,
								     learn_gap, learn_emit, regularize, pseudocounts, em_max_iter, em_min_inc);
    else
      THROWEXPR ("ERROR: Unreachable code.");
  }

}

Dotplot FSA::get_pairwise_dotplot (Params& params, const Sequence& xseq, const Sequence& yseq,
				   bool left_match, bool right_match, bool ragged_ends) {

  // enforce max_ram cutoff
  size_t needed = estimate_ram_needed (xseq.length(), yseq.length());
  if ((max_ram > 0) && (static_cast<int> (needed) > max_ram)) {
    CTAG(9,FSAEM) << "WARNING: Skipping inference step because estimated memory usage (" << needed << " MB) "
		  << "is greater than specified limit (" << max_ram << " MB)." << endl
		  << "Offending sequences are '" << xseq.name << "' and '" << yseq.name << "'." << endl;
    return Dotplot (xseq.length(), yseq.length());
  }

  // without banding
  if (bandwidth == 0) {
    if (is_dna && !is_indel2)
      return Model::get_pairwise_dotplot_engine<NucleotideWithoutBanding> (params, xseq, yseq,
									   left_match, right_match, ragged_ends);
    else if (is_dna && is_indel2)
      return Model::get_pairwise_dotplot_engine<NucleotideIndel2WithoutBanding> (params, xseq, yseq,
										 left_match, right_match, ragged_ends);
    else if (!is_dna && !is_indel2)
      return Model::get_pairwise_dotplot_engine<AminoAcidWithoutBanding> (params, xseq, yseq,
									  left_match, right_match, ragged_ends);
    else if (!is_dna && is_indel2)
      return Model::get_pairwise_dotplot_engine<AminoAcidIndel2WithoutBanding> (params, xseq, yseq,
										left_match, right_match, ragged_ends);
    else
      THROWEXPR ("ERROR: Unreachable code.");
  }
  // with banding
  else {
    if (is_dna && !is_indel2)
      return Model::get_pairwise_dotplot_engine<NucleotideWithBanding> (params, xseq, yseq,
									left_match, right_match, ragged_ends);
    else if (is_dna && is_indel2)
      return Model::get_pairwise_dotplot_engine<NucleotideIndel2WithBanding> (params, xseq, yseq,
									      left_match, right_match, ragged_ends);
    else if (!is_dna && !is_indel2)
      return Model::get_pairwise_dotplot_engine<AminoAcidWithBanding> (params, xseq, yseq,
								       left_match, right_match, ragged_ends);
    else if (!is_dna && is_indel2)
      return Model::get_pairwise_dotplot_engine<AminoAcidIndel2WithBanding> (params, xseq, yseq,
									     left_match, right_match, ragged_ends);
    else
      THROWEXPR ("ERROR: Unreachable code.");
  }

}

Post_probs FSA::get_pairwise_post_probs (Params& params, const Sequence& xseq, const Sequence& yseq,
					 bool left_match, bool right_match, bool ragged_ends) {

  // enforce max_ram cutoff
  const size_t needed = estimate_ram_needed (xseq.length(), yseq.length());
  if ((max_ram > 0) && (needed > static_cast<size_t> (max_ram))) {
    CTAG(9,FSAEM) << "WARNING: Skipping inference step because estimated memory usage (" << needed << " MB) "
		  << "is greater than specified limit (" << max_ram << " MB)." << endl
		  << "Offending sequences are '" << xseq.name << "' and '" << yseq.name << "'." << endl;
    return Post_probs();
  }

  // without banding
  if (bandwidth == 0) {
    if (is_dna && !is_indel2)
      return Model::get_pairwise_post_probs_engine<NucleotideWithoutBanding> (params, xseq, yseq,
									      minprob,
									      left_match, right_match, ragged_ends);
    else if (is_dna && is_indel2)
      return Model::get_pairwise_post_probs_engine<NucleotideIndel2WithoutBanding> (params, xseq, yseq,
										    minprob,
										    left_match, right_match, ragged_ends);  
    else if (!is_dna && !is_indel2)
      return Model::get_pairwise_post_probs_engine<AminoAcidWithoutBanding> (params, xseq, yseq,
									     minprob,
									     left_match, right_match, ragged_ends);
    else if (!is_dna && is_indel2)
      return Model::get_pairwise_post_probs_engine<AminoAcidIndel2WithoutBanding> (params, xseq, yseq,
										   minprob,
										   left_match, right_match, ragged_ends);
    else
      THROWEXPR ("ERROR: Unreachable code.");
  }
  // with banding
  else {
    if (is_dna && !is_indel2)
      return Model::get_pairwise_post_probs_engine<NucleotideWithBanding> (params, xseq, yseq,
									   minprob,
									   left_match, right_match, ragged_ends);
    else if (is_dna && is_indel2)
      return Model::get_pairwise_post_probs_engine<NucleotideIndel2WithBanding> (params, xseq, yseq,
										 minprob,
										 left_match, right_match, ragged_ends);
    else if (!is_dna && !is_indel2)
      return Model::get_pairwise_post_probs_engine<AminoAcidWithBanding> (params, xseq, yseq,
									  minprob,
									  left_match, right_match, ragged_ends);
    else if (!is_dna && is_indel2)
      return Model::get_pairwise_post_probs_engine<AminoAcidIndel2WithBanding> (params, xseq, yseq,
										minprob,
										left_match, right_match, ragged_ends);
    else
      THROWEXPR ("ERROR: Unreachable code.");
  }

}

void FSA::show_divergences (std::ostream& o) const {

  // first compute a matrix of divergences (-1 for no information)
  std::vector<std::vector<double> > distmat (num_seqs, std::vector<double> (num_seqs, -1));
  for (size_t cnt = 0; cnt < alignment_seq_pairs.size(); ++cnt) {
    const size_t i = alignment_seq_pairs[cnt].first;
    const size_t j = alignment_seq_pairs[cnt].second;
    distmat[i][j] = alignment_params[cnt].branch_length (0);
    distmat[j][i] = alignment_params[cnt].branch_length (1);
  }

  // now display: first get width
  size_t width = 0;
  for (size_t i = 0; i < seq_db_internal.size(); ++i)
    width = max (width, seq_db_internal.get_seq (i).name.length());
  width += 2;

  // state label row
  o.width (width); o << "";
  for (size_t i = 0; i < num_seqs; ++i) {
    o.width (width);
    o << seq_db_internal.get_seq (i).name;
  }
  o << endl;
  // divergences
  for (size_t i = 0; i < num_seqs; ++i) {
    o.width (width);
    o << seq_db_internal.get_seq (i).name;
    for (size_t j = 0; j < num_seqs; ++j) {
      o.width (width);
      if (i == j)
	o << 0.;
      else
	o << std::setprecision (PRECISION_DEFAULT) << distmat[i][j];
    }
    o << endl;
  }

}

size_t FSA::estimate_ram_needed (size_t xlen, size_t ylen) const {

  size_t needed = static_cast<size_t> (std::floor (((xlen + 1.0) / (1.0 * MEGABYTE)) * (ylen + 1) * static_cast<double> (sizeof (bfloat))));
  needed *= is_indel2 ? 5 : 3; // hardcode in size of Pair HMM state space (not counting state and end)
  needed *= 2;                 // factor of 2 for when we store both forward and backward DP matrices (for posterior decoding)

  return needed;
}

void FSA::assemble_sequence_data() {

  // assemble clean sequence data
  // (hardmasked sequence stripped if so requested and degenerate characters randomized)
  for (size_t i = 0; i < num_seqs; ++i) {

    // get original sequence data
    Sequence& sequence = seq_db.get_seq (i);

    // if hardmasked, keep track of hardmasking
    if (hardmasked)
      sequence.init_hardmasking (MIN_HARDMASK_LENGTH, DNA_alphabet::is_hardmask_char);

    // if hardmasking, store stripped sequence
    std::string orig_seq = hardmasked ? sequence.get_stripped_sequence().seq : sequence.seq;

    // remove degenerate characters
    std::string nondegen_seq;
    if (nucprot)    // handle case of nucprot (store the translated sequence)
      nondegen_seq = Protein_alphabet().get_nondegen ((Translated_sequence (Sequence (sequence.name, orig_seq))).get_forward (0).seq);
    else
      nondegen_seq = alphabet.get_nondegen (orig_seq);

    // store
    Sequence* nondegen_sequence = new Sequence (sequence.name, nondegen_seq);
    seq_db_internal.add_seq (nondegen_sequence);
  }

}

void FSA::choose_seq_pairs() {

  // initialize Sequence_pair_selector
  Sequence_pair_selector sequence_pair_selector (seq_db_internal,
						 alphabet,
						 kmer_length);

  // refalign first
  if (refalign) {

    // first sequence lies at the center if --refalign is select
    sequence_pair_selector.choose_palm_tree (seq_db_internal.get_seq (0).name);

  }

  // if all are requested, then go ahead and select all pairs
  else if (static_cast<size_t> (num_alignment_pairs) == num_seq_pairs) {

      sequence_pair_selector.choose_all();

  }

  // else select sequence pairs
  // order as:
  //  - minimum spanning tree
  //  - maximum spanning tree
  //  - minimum spanning palm tree
  //  - degree per sequence
  //  - randomly-chosen pairs
  else {

    // first add sequence pairs in a structured manner
    if (num_minimum_spanning_trees > 0)
      sequence_pair_selector.choose_minimum_spanning_tree (num_minimum_spanning_trees);
    if (num_maximum_spanning_trees > 0)
      sequence_pair_selector.choose_maximum_spanning_tree (num_maximum_spanning_trees);
    if (num_minimum_spanning_palm_trees > 0)
      sequence_pair_selector.choose_minimum_spanning_palm_tree (num_minimum_spanning_palm_trees);
    if (degree > 0)
      sequence_pair_selector.choose_kmer_similarity (degree);

    // then add randomly until we hit the requested num_alignment_pairs target
    const size_t num_selected = sequence_pair_selector.num_selected();
    if (num_selected < static_cast<size_t> (num_alignment_pairs))
      sequence_pair_selector.choose_random (num_alignment_pairs - num_selected);

  }
  
  // now actually store the pairs
  Sequence_pairs selected_sequence_pairs = sequence_pair_selector.get_chosen_sequence_pairs();
  assert (selected_sequence_pairs.size() == static_cast<size_t> (num_alignment_pairs));
  size_t num_added_alignment_pairs = 0;
  for (Sequence_pairs::const_iterator sequence_pair = selected_sequence_pairs.begin(); sequence_pair != selected_sequence_pairs.end(); ++sequence_pair) {

    const size_t i = sequence_pair->first;
    const size_t j = sequence_pair->second;
    assert (i < j);

    // now actually store the sequence pair
    // if a worker, do we store this seq pair?
    if (w_worker) {
      if (is_valid_worker_seq_pair (++num_added_alignment_pairs, i, j))
	alignment_seq_pairs.push_back (*sequence_pair);
    }

    // if not parallelizing, just store it
    else
      alignment_seq_pairs.push_back (*sequence_pair);

  }

}

bool FSA::is_valid_worker_seq_pair (const int cnt, const int i, const int j) const {

  int start_pos = (int) (w_start_seq_pair.first * (2 * num_seqs - w_start_seq_pair.first - 1) / 2) + (w_start_seq_pair.second - w_start_seq_pair.first);
  int input_pos = (int) (i * (2 * num_seqs - i - 1) / 2) + (j - i);

  if (cnt != -1) {
	if (w_prev_seq_pairs_sum <= cnt && cnt < (w_prev_seq_pairs_sum + w_num_seq_pairs))
	  return true;
  } else {
    if (start_pos <= input_pos && input_pos < (start_pos + w_num_seq_pairs))
      return true;
  }
  return false;
}

void FSA::init_for_mw_worker (const std::pair<int, int> start_seq_pair, const int prev_seq_pairs_sum, const int num_seq_pairs) {
  w_start_seq_pair = start_seq_pair;
  w_prev_seq_pairs_sum = prev_seq_pairs_sum;
  w_num_seq_pairs = num_seq_pairs;
  w_worker = true;
}

Post_probs FSA::perform_pairwise_inference (Params& params, const Sequence& xseq, const Sequence& yseq,
					    const bool left_match, const bool right_match, const bool ragged_ends,
					    const Params& pseudocounts) {

  const size_t xlen = xseq.length();
  const size_t ylen = yseq.length();

  Sequence_database seq_db_pair;
  seq_db_pair.add_seq (xseq);
  seq_db_pair.add_seq (yseq);
  Sequence_pairs both;
  both.push_back (std::make_pair (0, 1));

  // train params and get alignment posteriors:
  // only train if there is enough sequence data
  if (learn_gap || learn_emit_bypair) {
    // 2 cases where we learn:
    // 1) learn emit and possibly gap:
    //     only learn emit if there's enough sequence data (which implies sufficient data for learning gap)
    // 2) learn only gap
    //     only learn gap if there's enough sequence data
    if (learn_emit_bypair && ((0.5 * (xlen + ylen)) > min_emit_training_data))
      train_params (params, seq_db_pair, both,
		    left_match, right_match, ragged_ends,
		    pseudocounts, learn_emit_bypair);
    else if (learn_gap && ((0.5 * (xlen + ylen)) > min_gap_training_data))
      train_params (params, seq_db_pair, both,
		    left_match, right_match, ragged_ends,
		    pseudocounts, false);
  }

  // perform inference
  // note that we don't need to enforce lexical ordering here (it's guaranteed by get_pairwise_post_probs())
  return get_pairwise_post_probs (params, xseq, yseq,
				  left_match, right_match, ragged_ends);

}

Post_probs FSA::perform_anchored_pairwise_inference (const Params& params_seed, const Sequence& xseq, const Sequence& yseq,
						     const Params& pseudocounts,
						     const Anchors& anchors) {

  const size_t xlen = xseq.length();
  const size_t ylen = yseq.length();

  // hold posterior information for this sequence pair
  Post_probs post_probs;

  // parameters for anchored subsequence pairs
  // (modified for each pair)
  Params params;

  // loop through the anchors for this sequence pair:
  // look at the subsequence from the previous to the current anchor
  // This all relies on the list of anchors being properly sorted, but we guarantee this
  // during our anchor resolution.
  unsigned xprevright = 0; // right-hand X coordinate of previous anchor (closed interval)
  unsigned yprevright = 0; // set to fake value before sequence start
  for (size_t a = 0; a <= anchors.size(); ++a) { // we're going over the bound of Anchors anchors here
                                                 // in order to catch the last anchored subsequence pair

    // specify homology assumptions beyond (to left and right of) current subsequences
    // note that this gives appropriate settings for the case of no anchors as well
    bool left_match = true;          // for interior anchored pairs, homology is assumed to both left and right
    bool right_match = true;
    if (a == 0)                      // if first anchored pair, then no homology assumed to left
      left_match = false;
    if (a == anchors.size())         // if last anchored pair, then no homology assumed to right
      right_match = false;

    // to do:
    // make concrete decision about whether to use {left,right}_match information
    // for now don't force homology
    left_match = right_match = false;

    // calculate the minimum required subsequence length
    // (left_match forces Start -> Match and right_match forces Match -> End)
    const size_t minsublen = static_cast<size_t> (left_match) + static_cast<size_t> (right_match);

    // get coordinates for current anchor, making sure to catch end case,
    // for which we have no right-anchoring information
    unsigned xcurrleft, xcurrright;  // current anchor is the closed interval [xcurrleft, xcurrright]
    unsigned ycurrleft, ycurrright;
    if (a == anchors.size()) {       // if last anchored pair,
      xcurrleft = xlen;              // create a fake anchor which lies just beyond the sequence boundaries (0-based indexing)
      xcurrright = xlen - 1;         // (fake value)
      ycurrleft = ylen;
      ycurrright = ylen - 1;
    } else {
      xcurrleft = anchors[a].xcoords.start;
      xcurrright = anchors[a].xcoords.end;
      ycurrleft = anchors[a].ycoords.start;
      ycurrright = anchors[a].ycoords.end;
    }

    // catch the case of the last real anchor (a == anchors.size() - 1) ending at a sequence boundary
    // (no need to perform inference since overhanging sequence will be left unaligned regardless)
    if (a == anchors.size()
	&& ((xprevright == xlen) || (yprevright == ylen)))
      continue;

    // assert that anchors are non-overlapping (at a pairwise level)
    if (a != 0) {
      if (xprevright >= xcurrleft)
	THROWEXPR ("ERROR: Overlapping anchors: xprevright = " << xprevright << "; xcurrleft = " << xcurrleft << endl
		   << anchors[a]);
      if (yprevright >= ycurrleft)
	THROWEXPR ("ERROR: Overlapping anchors: yprevright = " << yprevright << "; ycurrleft = " << ycurrleft << endl
		   << anchors[a]);
    }

    // create subsequences (xsubleft + xsublen) and (ysubleft + ysublen)

    unsigned xsubleft, ysubleft; // left coordinate of closed interval
    unsigned xsublen, ysublen;   // length of closed interval

    // if first anchored pair, then make sure to get first character in subseq
    xsubleft = xprevright + (a != 0 ? 1 : 0); // this catches the case of an anchor beginning at the coordinate 0,
    ysubleft = yprevright + (a != 0 ? 1 : 0); // since we initialize these to 0 before looping over anchors (so we don't want to increment the boundary)

    // cover case of abutting (adjacent) anchors
    xsubleft = (xsubleft > xcurrleft) ? xcurrleft : xsubleft;
    ysubleft = (ysubleft > ycurrleft) ? ycurrleft : ysubleft;

    xsublen = xcurrleft - xsubleft; // we don't include the anchors themselves
    ysublen = ycurrleft - ysubleft; // the subseq [xsubleft + xsublen] is a closed interval (this is important to remember!)

    // subsequence names and subsequences themselves
    // cast to int to cover the case of 0-length subseqs at the beginning,
    // in which case the second coordinate is negative
    const Sequence xsubseq (xseq.name + "_[" + Util::to_string (xsubleft) + "-" + Util::to_string (static_cast<int> (xsubleft + xsublen) - 1) + "]", xseq.seq.substr (xsubleft, xsublen));
    const Sequence ysubseq (yseq.name + "_[" + Util::to_string (ysubleft) + "-" + Util::to_string (static_cast<int> (ysubleft + ysublen) - 1) + "]", yseq.seq.substr (ysubleft, ysublen));

    Sequence_database subseq_pair;
    subseq_pair.add_seq (xsubseq);
    subseq_pair.add_seq (ysubseq);
    Sequence_pairs both;
    both.push_back (std::make_pair (0, 1));

    // cover case of empty subsequences:
    //  (occur when anchors begin or end at sequence boundaries or are adjacent)
    // store anchor indices and update posteriors,
    // then continue to next anchor
    if ((xsublen <= minsublen) || (ysublen <= minsublen)) {
      // update posteriors to record current anchor info if not a fake anchor
      if (a < anchors.size()) {
	for (unsigned x = 0; x <= xcurrright - xcurrleft; ++x) // closed interval!
	  post_probs.push_back (Post_prob (xcurrleft + x, ycurrleft + x, anchors[a].score_to_prob())); // force high probability for anchor
      }
      // increment the anchor indices
      xprevright = xcurrright;
      yprevright = ycurrright;
      // log (but only if at least one sequence has length > 0)
      if (xsublen + ysublen > 0) {
	const double anchors_done = 100.0 * (a + 1) / (anchors.size() + 1);
	if (CTAGGING(6,FSA))
	  CTAG(6,FSA) << "Processed anchored subsequence pair '" << xsubseq.name << "' and '" << ysubseq.name << "'; "
		      << std::floor (anchors_done + 0.5) << "% (" << a + 1 << "/" << anchors.size() + 1 << ") complete."
		      << endl;
      }
      // skip inference step
      continue;
    }

    // re-initialize params for this anchored pair of subseqs
    params.copy_all (params_seed);

    // train params and get alignment posteriors
    // only train if there is enough sequence data
    if (learn_gap || learn_emit_bypair) {
      // 2 cases where we learn:
      // 1) learn emit and possibly gap:
      //     only learn emit if there's enough sequence data (which implies sufficient data for learning gap)
      // 2) learn only gap
      //     only learn gap if there's enough sequence data
      if (learn_emit_bypair && ((0.5 * (xsublen + ysublen)) > min_emit_training_data))
	train_params (params, subseq_pair, both,
		      left_match, right_match, ragged_ends,
		      pseudocounts, learn_emit_bypair);
      else if (learn_gap && ((0.5 * (xsublen + ysublen)) > min_gap_training_data))
	train_params (params, subseq_pair, both,
		      left_match, right_match, ragged_ends,
		      pseudocounts, false);
    }

    // perform inference
    Post_probs subseq_post_probs = get_pairwise_post_probs (params, xsubseq, ysubseq,
							    left_match, right_match, ragged_ends);

    // store this subsequence information in the posterior maps for the entire sequence pair:
    // first store the posteriors calculated for these subsequences
    for (Post_probs::iterator iter = subseq_post_probs.begin(); iter != subseq_post_probs.end(); ++iter) {
      const unsigned subx = iter->x;    // coordinates within the subseqs
      const unsigned suby = iter->y;
      const double& p = iter->prob;
      const unsigned x = xsubleft + subx; // coordinates for the entire sequences
      const unsigned y = ysubleft + suby;
      assert (x < xlen);       // (0-based coordinates)
      assert (y < ylen);
      post_probs.push_back (Post_prob (x, y, p));
    }

    // then fix the current anchors
    //  (unless we're at the end case, for which we've created a fake anchor beyond the sequence boundaries)
    if (a < anchors.size()) {
      for (unsigned x = 0; x <= xcurrright - xcurrleft; ++x) // closed interval!
	post_probs.push_back (Post_prob (xcurrleft + x, ycurrleft + x, anchors[a].score_to_prob())); // force high probability for anchor
    }

    // log progress through anchors
    const unsigned anchors_done = static_cast<unsigned> (std::floor ((100.0 * (a+1) / (anchors.size() + 1)) + 0.5));
    if (CTAGGING(6,FSA))
      CTAG(6,FSA) << "Processed anchored subsequence pair '" << xsubseq.name << "' and '" << ysubseq.name << "'; "
		   << anchors_done << "% (" << a+1 << "/" << anchors.size() + 1 << ") complete."
		   << endl;

    // increment the anchor indices
    xprevright = xcurrright;
    yprevright = ycurrright;

  }

  // Enforce lexical ordering of post_probs:
  // Unfortunately, storing anchor posteriors before subsequence posteriors isn't sufficient to preserve
  // lexical ordering; overlapping anchors will still cause problems.  We therefore force a sort
  // before constructing the SparseMatrix.  This costs O (n log n) time, but is unavoidable
  // unless we deal more cleverly with overlapping anchors.
  std::sort (post_probs.begin(), post_probs.end());

  return post_probs;

}

void FSA::load_probabilities_from_file (const std::string& filename,
					std::vector<std::vector<SparseMatrix*> >& sparse_matrices) {

  CTAG(8,FSA) << "Reading posterior probabilities file '" << filename << "'." << endl;

  // format is:
  // ; match posteriors
  // (0, 0) ~ (1, 0) => 0.999999
  // ; gap posteriors
  // (0, 0) ~ (1, -1) => 0.0001
  Regexp re_probs ("^[ \t]*\\([ \t]*([0-9]+)[ \t]*,[ \t]*([0-9]+)[ \t]*\\)[ \t]*~[ \t]*\\([ \t]*([0-9]+)[ \t]*,[ \t]*([0-9]+)[ \t]*\\)[ \t]*=>[ \t]*(.*)$");

  std::vector<std::vector<Post_probs> > post_probs_storage (num_seqs, std::vector<Post_probs> (num_seqs));

  std::ifstream filestream;
  filestream.open (filename.c_str(), std::ios::in);
  if (!filestream.is_open())
    THROWEXPR ("ERROR: Couldn't open file '" << filename << "' for reading.")

  std::string line;
  while (!filestream.eof()) {

    getline (filestream, line);
    Util::chomp (line);
    if (!line.length()) { continue; }

    // are we at a match line?
    if (re_probs.Match (line.c_str())) {

      // check that format is correct
      if (re_probs.SubStrings() != 5) {
	CTAG(8,FSA) << "WARNING: Couldn't parse line: " << line << endl;
	continue;
      }

      // pull out data
      const size_t i = static_cast<size_t> (atoi (re_probs[1].c_str()));
      const unsigned ii = static_cast<unsigned> (atoi (re_probs[2].c_str()));
      const size_t j = static_cast<size_t> (atoi (re_probs[3].c_str()));
      const unsigned jj = static_cast<unsigned> (atoi (re_probs[4].c_str()));
      const float prob = atof (re_probs[5].c_str());

      // check sane
      if (i >= num_seqs || j >= num_seqs)
	THROWEXPR ("ERROR: Invalid sequence indices i = " << i << ", j = " << j << ".");
      if (ii >= seq_db.get_seq (i).length())
	THROWEXPR ("ERROR: Sequence position out of bounds: (" << i << ", " << ii << ").");
      if (jj >= seq_db.get_seq (j).length())
	THROWEXPR ("ERROR: Sequence position out of bounds: (" << j << ", " << jj << ").");

      // store
      post_probs_storage[i][j].push_back (Post_prob (ii, jj, prob));

    }

  }
  filestream.close();

  // initialize SparseMatrix objects
  for (size_t i = 0; i < num_seqs; ++i) {
    for (size_t j = 0; j < num_seqs; ++j) {
      Post_probs& post_probs = post_probs_storage[i][j];
      // ensure we don't attempt to store empty data (potentially overwriting good data)
      if (i == j || !post_probs.size())
	continue;
      // must sort before initializing SparseMatrix
      std::sort (post_probs.begin(), post_probs.end());
      sparse_matrices[i][j] = new SparseMatrix (i, j,
						seq_db.get_seq (i).length(), seq_db.get_seq (j).length(),
						post_probs);
      // clean up as we go to save memory
      post_probs.clear();
      // get transpose too
      sparse_matrices[j][i] = sparse_matrices[i][j]->ComputeTranspose();

      // check that we've actually created a matrix!
      if (!sparse_matrices[i][j]->size())
	THROWEXPR ("ERROR: Failed to create sparse matrix of posterior probabilities.");

    }
  }

}

void FSA::build_multiple_alignment() {

  // read in weights for sequence pairs if information present
  if (tree_weights_file != "")
    tree_weights.from_file (seq_db, tree_weights_file);

  // set up model topology
  bool left_match = false;  // do not require homology beyond sequence boundaries
  bool right_match = false;

  // initialize parameter seed
  Params params_seed;
  const std::vector<double> char_freq_total (Model::count_char_freq (seq_db_internal, this->alphabet_string, is_dna));
  params_seed.bandwidth = bandwidth;
  init_indel_params (params_seed);
  init_subst_params (params_seed, char_freq_total, model, time);

  // initialize pseudocounts (for regularization) according to the params_seed distribution
  Params pseudocounts;
  init_pseudocounts (pseudocounts, params_seed, regularization_transition_scale, regularization_emission_scale);

  // if requested, learn params_seed over alignment_seq_pairs
  if (learn_emit_all) {
    train_params (params_seed, seq_db_internal, alignment_seq_pairs,
		  left_match, right_match, ragged_ends,
		  pseudocounts, learn_emit_all);
  }

  // store sparse matrices for annealing: store only sparse matrices for memory efficiency
  std::vector<std::vector<SparseMatrix*> > sparse_matrices (num_seqs, std::vector<SparseMatrix*> (num_seqs, reinterpret_cast<SparseMatrix*> (NULL)));

  // option list to be used for a database connection
  DB_opts db_opts;
  db_opts.copy_opts (this);

  // create database manager 
  Manager manager (seq_db_internal, db_opts);

  if (parallelizing) 
    manager.mw_master_run (opts.init_argc, opts.init_argv, params_seed, pseudocounts);
  else if (write_db)
    manager.mw_single_worker_run (params_seed, pseudocounts);

  if (noannealing)
    return;

  // if available, get pairwise posterior probabilities from the database
  if (manager.is_sparse_matrices_available()) {
#if defined(HAVE_CONDOR) || defined(HAVE_POSTGRES)
    CTAG(9,FSA) << "Getting pairwise posterior probabilities from database." << endl;
    manager.get_sparse_matrices (sparse_matrices);
    CTAG(8,FSA) << "Got pairwise posterior probabilities." << endl;
#endif
  }

  // if requested, load pairwise probabilities from a file
  // rather than re-estimating them with DP
  else if (load_probs_file.length()) {

    load_probabilities_from_file (load_probs_file,
				  sparse_matrices);

  }

  // else perform alignment as usual
  else {

    // initialize params for each alignment_seq_pair
    alignment_params.resize (alignment_seq_pairs.size());
    for (size_t cnt = 0; cnt < alignment_seq_pairs.size(); ++cnt)
      alignment_params[cnt].copy_all (params_seed);

    // now do all-pairs comparison
    CTAG(9,FSA) << "Collecting pairwise posterior probabilities." << endl;

    // loop through sequence database
    for (size_t cnt = 0; cnt < alignment_seq_pairs.size(); ++cnt) {

      // get the sequence indices
      const size_t i = alignment_seq_pairs[cnt].first;
      const size_t j = alignment_seq_pairs[cnt].second;

      // initialize all the sequence data
      const Sequence& xseq = seq_db_internal.get_seq (i);
      const Sequence& yseq = seq_db_internal.get_seq (j);

      // if there is an empty sequence, then skip inference
      if (!xseq.length() || !yseq.length())
	continue;

      // perform pairwise inference
      Post_probs post_probs = perform_pairwise_inference (alignment_params[cnt], xseq, yseq,
							  left_match, right_match, ragged_ends,
							  pseudocounts);

      // we need the original sequences in case we've been hardmasking
      const Sequence& xseq_orig = seq_db.get_seq (i);
      const Sequence& yseq_orig = seq_db.get_seq (j);

      // if hardmasking, map coords back to original sequence
      if (hardmasked) {
	for (Post_probs::iterator p = post_probs.begin(); p != post_probs.end(); ++p) {
	  p->x = xseq_orig.map_stripped_to_orig (p->x);
	  p->y = yseq_orig.map_stripped_to_orig (p->y);
	}
      }

      // create sparse matrix for sequence pair
      sparse_matrices[i][j] = new SparseMatrix (i, j,
						xseq_orig.length(), yseq_orig.length(),
						post_probs);
      sparse_matrices[j][i] = sparse_matrices[i][j]->ComputeTranspose(); // pre-compute the transpose for speed

      // check for success (whether we hit max_ram limit)
      // throw error if so
      if (!sparse_matrices[i][j]->size() && xseq_orig.length() && yseq_orig.length()) {
	CTAG(9,FSA) << "WARNING: Unable to detect any homology between sequences '" << xseq.name << "' and '" << yseq.name << "'."
		    << " This may be due to RAM constraints; check --maxram or try using anchoring." << endl;
	if (require_homology)
	  THROWEXPR ("ERROR: Alignment failure.");
      }

      // log progress through sequence pairs
      const unsigned percent_done = static_cast<unsigned> (std::floor ((100.0 * (cnt+1) / alignment_seq_pairs.size()) + 0.5));
      if (CTAGGING(7,FSA))
	CTAG(7,FSA) << "Processed sequence pair '" << xseq.name << "' and '" << yseq.name << "'; "
		    << percent_done << "% (" << cnt+1 << "/" << alignment_seq_pairs.size() << ") complete."
		    << endl;

    } // end loop over sequences

    // log
    CTAG(8,FSA) << "Processed a total of " << std::floor ((100.0 * alignment_seq_pairs.size() / num_seq_pairs) + 0.5) << "% ("
		<< alignment_seq_pairs.size() << "/" << num_seq_pairs << ") of all sequence pairs."
		<< endl;

  }

  // write trained params and/or dotplots for all seq pairs if requested
  if (write_params || write_dotplots) {
    CTAG(8,FSA) << "Writing learned parameters and/or dotplots to disk." << endl;
    std::string filename;

    for (size_t cnt = 0; cnt < alignment_seq_pairs.size(); ++cnt) {
      const Sequence& xseq = seq_db_internal.get_seq (alignment_seq_pairs[cnt].first);
      const Sequence& yseq = seq_db_internal.get_seq (alignment_seq_pairs[cnt].second);

      if (write_params) {
	// emission in human-readable format
	filename.clear();
	filename = std::string ("params.emission.learned")
	  + "." + xseq.name + "-" + yseq.name;
	alignment_params[cnt].write_emission (filename);

	// emission in .dat format
	filename.clear();
	filename = std::string ("params.emission.learned")
	  + "." + xseq.name + "-" + yseq.name
	  + ".dat";
	alignment_params[cnt].write_emission_dat (filename);

	// transition in human-readable format
	filename.clear();
	filename = std::string ("params.transition.learned")
	  + "." + xseq.name + "-" + yseq.name;
	alignment_params[cnt].write_transition (filename);
      }

      if (write_dotplots) {
	filename.clear();
	filename = std::string ("posterior.final.")
	  + xseq.name + '-' + yseq.name;
	get_pairwise_dotplot (alignment_params[cnt], xseq, yseq,
			      left_match, right_match, ragged_ends).write_dotplot (filename);
      }

    }
  }

  // write pairwise distance estimates
  if (write_divergences) {
    CTAG(8,FSA) << "Writing distance matrix to disk." << endl;
    const std::string filename = "divergences";
    std::ofstream file (filename.c_str());
    if (!file)
      THROWEXPR ("ERROR: Couldn't create file with name '" << filename << "'.");
    show_divergences (file);
    file.close();
  }

  // now do sequence annealing
  // handle case of nucprot
  // who knows why, but I got weird memory errors when I didn't explicitly create
  // a reference like this...
  const Sequence_database& aa_db = nucprot
    ? seq_db.translate()
    : Sequence_database();

  Alignment_DAG dag (nucprot ? aa_db : seq_db);
  dag.anneal (sparse_matrices, tree_weights,
	      manager,
	      use_tgf,
	      gap_factor, enable_dynamic_weights, 0,  // edge_weight_threshold = 0
	      num_refinement_steps,
	      output_for_gui, gui_prefix);
  dag.dfs_topological_sort();

  // we're done!
  Stockholm stock = nucprot
    ? dag.get_stockholm (sparse_matrices, tree_weights).get_codon_from_aa_alignment (seq_db)
    : dag.get_stockholm (sparse_matrices, tree_weights);

  if (write_stockholm)
    stock.write_stockholm (cout);
  else
    stock.write_mfa (cout);

}

void FSA::build_anchored_multiple_alignment() {

  // read in weights for sequence pairs if information present
  if (tree_weights_file != "")
    tree_weights.from_file (seq_db, tree_weights_file);

  // initialize parameter seed
  Params params_seed;
  const std::vector<double> char_freq_total (Model::count_char_freq (seq_db_internal, this->alphabet_string, is_dna));
  params_seed.bandwidth = bandwidth;
  init_indel_params (params_seed);
  init_subst_params (params_seed, char_freq_total, model, time);

  // initialize pseudocounts (for regularization) according to the params_seed distribution
  Params pseudocounts;
  init_pseudocounts (pseudocounts, params_seed, regularization_transition_scale, regularization_emission_scale);

  // store sparse matrices for annealing: store only sparse matrices for memory efficiency
  std::vector<std::vector<SparseMatrix*> > sparse_matrices (num_seqs, std::vector<SparseMatrix*> (num_seqs, reinterpret_cast<SparseMatrix*> (NULL)));

  // option list to be used for a database connection
  DB_opts db_opts;
  db_opts.copy_opts (this);

  // create database manager 
  Manager manager (seq_db_internal, db_opts);

  if (parallelizing)
    manager.mw_master_run (opts.init_argc, opts.init_argv, params_seed, pseudocounts);
  else if (write_db)
    manager.mw_single_worker_run (params_seed, pseudocounts);

  if (noannealing)
    return;

  // if available, get pairwise posterior probabilities from the database
  if (manager.is_sparse_matrices_available()) {
#if defined(HAVE_CONDOR) || defined(HAVE_POSTGRES)
    CTAG(9,FSA) << "Getting pairwise posterior probabilities from database." << endl;
    manager.get_sparse_matrices (sparse_matrices);
    CTAG(8,FSA) << "Got pairwise posterior probabilities." << endl;
#endif
  }

  // if requested, load pairwise probabilities from a file
  // rather than re-estimating them with DP
  else if (load_probs_file.length()) {

    load_probabilities_from_file (load_probs_file,
				  sparse_matrices);

  }

  // else perform alignment as usual
  else {

    // no need to bother with alignment_params
    // (different parameters are learned for each anchored subsequence pair)

    // handle case of nucprot
    const Sequence_database& aa_db = nucprot
      ? seq_db.translate()
      : Sequence_database();

    // get resolved anchors for the sequence pairs which we're going use for alignment
    Anchor_resolver anchor_resolver = nucprot
      ? Anchor_resolver (aa_db, seq_db_internal, alignment_seq_pairs)
      : Anchor_resolver (seq_db, seq_db_internal, alignment_seq_pairs);

    // add Mercator constraints if present
    if (mercator_constraint_file != "")
      anchor_resolver.add_mercator_constraints (mercator_constraint_file);

    // get resolved anchors
    CTAG(9,FSA) << "Getting anchors for " << std::floor ((100.0 * alignment_seq_pairs.size() / num_seq_pairs) + 0.5)
		<< "% (" << alignment_seq_pairs.size() << "/" << num_seq_pairs << ") of all sequence pairs."
		<< endl;

    // now actually get resolved anchors
    std::vector<Anchors> resolved_anchors_list = anchor_resolver.get_resolved_anchors (tree_weights,
										       anchor_minlen, anchor_max_join_length, use_translated,
										       exonerate, exonerate_minscore, softmasked,
										       hardmasked,
										       num_refinement_steps,
										       output_for_gui, gui_prefix);

    // now align the subsequences between the anchors:
    // do all-pairs comparison

    // log
    CTAG(9,FSA) << "Collecting pairwise posterior probabilities for each anchored subsequence pair." << endl;

    // loop through sequence database
    for (size_t cnt = 0; cnt < alignment_seq_pairs.size(); ++cnt) {

      // initialize all the sequence data
      const size_t i = alignment_seq_pairs[cnt].first;
      const size_t j = alignment_seq_pairs[cnt].second;

      const Sequence& xseq = seq_db_internal.get_seq (i);
      const Sequence& yseq = seq_db_internal.get_seq (j);

      // if there is an empty sequence, then skip inference
      if (!xseq.length() || !yseq.length())
	continue;

      // perform anchored pairwise inference
      Post_probs post_probs = perform_anchored_pairwise_inference (params_seed, xseq, yseq,
								   pseudocounts,
								   resolved_anchors_list[cnt]);

      // we need the original sequences in case we've been hardmasking
      const Sequence& xseq_orig = seq_db.get_seq (i);
      const Sequence& yseq_orig = seq_db.get_seq (j);

      // if hardmasking, map coords back to original sequence
      if (hardmasked) {
	for (Post_probs::iterator p = post_probs.begin(); p != post_probs.end(); ++p) {
	  p->x = xseq_orig.map_stripped_to_orig (p->x);
	  p->y = yseq_orig.map_stripped_to_orig (p->y);
	}
      }

      // create sparse matrix for sequence pair
      sparse_matrices[i][j] = new SparseMatrix (i, j,
						xseq_orig.length(), yseq_orig.length(), post_probs);
      sparse_matrices[j][i] = sparse_matrices[i][j]->ComputeTranspose();

      // check for success (whether we hit max_ram limit)
      // throw error if so
      if (!sparse_matrices[i][j]->size() && xseq_orig.length() && yseq_orig.length()) {
	CTAG(9,FSA) << "WARNING: Unable to detect any homology between sequences '" << xseq.name << "' and '" << yseq.name << "'."
		    << " This may be due to RAM constraints; check --maxram or try using anchoring." << endl;
	if (require_homology)
	  THROWEXPR ("ERROR: Alignment failure.");
      }

      // log progress through sequence pairs
      const unsigned percent_done = static_cast<unsigned> (std::floor ((100.0 * (cnt+1) / alignment_seq_pairs.size()) + 0.5));
      if (CTAGGING(7,FSA))
	CTAG(7,FSA) << "Processed sequence pair '" << xseq.name << "' and '" << yseq.name << "'; "
		    << percent_done << "% (" << cnt+1 << "/" << alignment_seq_pairs.size() << ") complete."
		    << endl;

    } // end loop over sequences

    // log
    CTAG(8,FSA) << "Processed a total of " << std::floor ((100.0 * alignment_seq_pairs.size() / num_seq_pairs) + 0.5) << "% ("
		<< alignment_seq_pairs.size() << "/" << num_seq_pairs << ") of all sequence pairs."
		<< endl;
  }

  // now do sequence annealing
  // handle case of nucprot
  // who knows why, but I got weird memory errors when I didn't explicitly create
  // a reference like this...
  const Sequence_database& aa_db = nucprot
    ? seq_db.translate()
    : Sequence_database();

  Alignment_DAG dag (nucprot ? aa_db : seq_db);
  dag.anneal (sparse_matrices, tree_weights,
	      manager,
	      use_tgf,
	      gap_factor, enable_dynamic_weights, 0,  // edge_weight_threshold = 0
	      num_refinement_steps,
	      output_for_gui, gui_prefix);
  dag.dfs_topological_sort();

  // we're done!
  Stockholm stock = nucprot
    ? dag.get_stockholm (sparse_matrices, tree_weights).get_codon_from_aa_alignment (seq_db)
    : dag.get_stockholm (sparse_matrices, tree_weights);

  if (write_stockholm)
    stock.write_stockholm (cout);
  else
    stock.write_mfa (cout);

}

int FSA::run() {

  try {
    // It's important that we input the data /before/ parsing the command-line options.
    // This allows us to set preset values for DNA, RNA and proteins and then let
    // them be overriden by the command-line options.
    init_opts();
    input_data();
    set_up_defaults();
    parse_opts();
    assemble_sequence_data();
    choose_seq_pairs();
    if (anchored)
      build_anchored_multiple_alignment();
    else
      build_multiple_alignment();
  }

  catch (const Dart_exception& e) {
    CLOGERR << e.what();
    THROWEXPR ("ERROR: Exception thrown.");
  }

  return 0;
}
