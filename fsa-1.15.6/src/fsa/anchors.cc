
/**
 * \file anchors.cc
 * This file is part of FSA, a sequence alignment algorithm.
 * \author Source code in this file was written by Robert Bradley.
 */

#include <cstdio>
#include <iostream>
#include <sys/types.h>
#include <sys/wait.h>
#include <algorithm>

#include "seq/alignment.h"
#include "fsa/anchors.h"

#define FILENAME_PREFIX "fsa_tmp_seq"
#define TMP_DIR "/tmp"

#define GUI_FILENAME_ANCHORS_SUFFIX ".anchors"

using namespace fsa;

// Note that unless specified otherwise, we use a 0-based coordinate system.
// The SparseMatrix code uses a 1-based system, but its constructor expects
// 0-based coordinates as input and the conversion happens behind the scenes.

// constructor
Anchor::Anchor (Interval xcoords, Interval ycoords)
  : xcoords (xcoords), ycoords (ycoords),
    __score (-1.),
    __external_scoring (false),
    __immutable (false) {

  update();

}

// constructor
Anchor::Anchor (Interval xcoords, Interval ycoords, double score)
  : xcoords (xcoords), ycoords (ycoords),
    __score (score),
    __external_scoring (false),
    __immutable (false) {

  update();

}

void Anchor::update() {

  // ensure valid coordinates
  assert ((xcoords.start >= 0) && (xcoords.end >= xcoords.start));
  assert ((ycoords.start >= 0) && (ycoords.end >= ycoords.start));
  assert ((xcoords.end - xcoords.start) == (ycoords.end - ycoords.start));

  x = static_cast<unsigned> (std::floor (0.5 * (xcoords.end + xcoords.start)));
  y = static_cast<unsigned> (std::floor (0.5 * (ycoords.end + ycoords.start)));

  assert (x - y == xcoords.end - ycoords.end);

  length = xcoords.end - xcoords.start + 1;
  assert (length > 0);

}

double Anchor::p_value (const std::string& xseq, const std::string& yseq) const {
  return p_value (&xseq, &yseq);
}

double Anchor::p_value (const std::string* xseq, const std::string* yseq) const {

  // extract anchored subseqs
  const std::string& xsubseq = xseq->substr (xcoords.start, length);
  const std::string& ysubseq = yseq->substr (ycoords.start, length);

  // calculate hamming distance (# of mismatches)
  size_t h = 0;
  for (size_t i = 0; i < length; ++i) {
    if (xsubseq[i] != ysubseq[i])
      ++h;
  }

  // get p-value (normalize to [0,1])
  const double pval =  1.0 - static_cast<double> (length) / static_cast<double> (std::max (xseq->length(), yseq->length()));
  assert ((pval >= 0.) && (pval <= 1.0));

  return pval;
}

void Anchors::show (std::ostream& o) const {
  for (std::vector<Anchor>::const_iterator anchor = begin(); anchor != end(); ++anchor)
    o << *anchor << endl;
}

void Anchors::write_gui_output (std::ostream& o, const size_t xidx, const size_t yidx) const {
  for (std::vector<Anchor>::const_iterator anchor = this->begin(); anchor != this->end(); ++anchor)
    o << "(" << xidx << " ~ " << yidx << ") == " << *anchor << endl;
}

Anchors Anchors::convert_constraints_to_anchors (const Constraints& constraints, const Sequence& xseq_orig, const Sequence& yseq_orig, const bool hardmasked) {

  Anchors anchors;
  for (std::vector<Constraint>::const_iterator constraint = constraints.begin(); constraint != constraints.end(); ++constraint) {

    // left end
    Anchor left (Interval (constraint->xcoords.start, constraint->xcoords.start),
		 Interval (constraint->ycoords.start, constraint->ycoords.start));
    left.set_immutable();

    // right end
    Anchor right (Interval (constraint->xcoords.end, constraint->xcoords.end),
		  Interval (constraint->ycoords.end, constraint->ycoords.end));
    right.set_immutable();

    // move anchors if necessary to get out of hardmasked region
    bool is_anchor_entirely_hardmasked = false;
    if (hardmasked) {

      // if necessary, increment left boundary along the sequences until it doesn't fall in a hardmasked region
      bool shifted = false; // for logging
      while (xseq_orig.is_pos_hardmasked (left.x) || yseq_orig.is_pos_hardmasked (left.y)) {
	shifted = true;
	left.xcoords.start = left.xcoords.end = (left.x + 1);
	left.ycoords.start = left.ycoords.end = (left.y + 1);
	left.update();
	// check whether the increment boundary is valid
	// (falls within sequence limits and not past right-hand boundary)
	if ((left.x >= xseq_orig.length()) || (left.x >= right.x)) {
	  is_anchor_entirely_hardmasked = true;
	  break;
	}
      }

      if (shifted && !is_anchor_entirely_hardmasked)
	CTAG(6,ANCHORING) << "WARNING: Shifted Mercator constraint boundary which fell in hardmasked region; new left-hand boundary is: " << left << endl;

      // if necessary, decrement right boundary along the sequences until doesn't fall in a hardmasked region
      shifted = false;
      while (xseq_orig.is_pos_hardmasked (right.x) || yseq_orig.is_pos_hardmasked (right.y)) {
	shifted = true;
	right.xcoords.start = right.xcoords.end = (right.x - 1);
	right.ycoords.start = right.ycoords.end = (right.y - 1);
	right.update();
	if ((right.x >= xseq_orig.length()) || (left.x >= right.x)) {
	  is_anchor_entirely_hardmasked = true;
	  break;
	}
      }

      if (shifted && !is_anchor_entirely_hardmasked)
	CTAG(6,ANCHORING) << "WARNING: Shifted Mercator constraint boundary which fell in hardmasked region; new right-hand boundary is: " << right << endl;

    }

    // if anchor is entirely hardmasked in one sequence, then ignore it
    if (is_anchor_entirely_hardmasked) {
      CTAG(6,ANCHORING) << "WARNING: Mercator constraint fell entirely within a hardmasked region; ignoring it entirely." << endl;
      continue;
    }

    // store anchors
    anchors.store (left);
    anchors.store (right);

  }

  // if hardmasking, map Mercator coordinates to stripped sequence
  if (hardmasked) {
    for (std::vector<Anchor>::iterator anchor = anchors.begin(); anchor != anchors.end(); ++anchor) {
      anchor->xcoords.start = xseq_orig.map_orig_to_stripped (anchor->xcoords.start);
      anchor->xcoords.end = xseq_orig.map_orig_to_stripped (anchor->xcoords.end);
      anchor->ycoords.start = yseq_orig.map_orig_to_stripped (anchor->ycoords.start);
      anchor->ycoords.end = yseq_orig.map_orig_to_stripped (anchor->ycoords.end);
      anchor->update();
      // sanity check
      if ((anchor->xcoords.end == xseq_orig.length()) || (anchor->ycoords.end == yseq_orig.length())) {
	THROWEXPR ("ERROR: Mercator constraint falls in a hardmasked interval:" << endl << *anchor);
      }
    }
  }

  return anchors;

}

const Post_probs Anchors::convert_to_post_probs() const {

  // store each Anchor as an entry
  Post_probs post_probs;
  for (std::vector<Anchor>::const_iterator anchor = begin(); anchor != end(); ++anchor)
    post_probs.push_back (Post_prob (anchor->x, anchor->y, anchor->get_score()));

  // now enforce lexical ordering on centroids
  // (required by SparseMatrix constructor)
  // forgetting this step caused me BIG PAIN!
  std::sort (post_probs.begin(), post_probs.end());

  return post_probs;
}

bool Anchors::create_spanning_map() {

  __spanning_map.clear();

  bool nondegen = true;
  for (size_t i = 0; i < this->size(); ++i)
    nondegen = __spanning_map.insert (std::make_pair (std::make_pair ((*this)[i].x, (*this)[i].y), i)).second;

  return nondegen;
}

size_t Anchors::exists_spanning_anchor (const unsigned x, const unsigned y) const {

  if (__spanning_map.find (std::make_pair (x, y)) != __spanning_map.end())
    return (*__spanning_map.find (std::make_pair (x, y))).second;

  return size();

}

const Anchor Anchors::get_spanning_anchor (const unsigned x, const unsigned y) {

  const size_t idx = exists_spanning_anchor (x, y);
  assert (idx < this->size());
  assert (((*this)[idx].xcoords.start <= x) && ((*this)[idx].xcoords.end >= x));
  assert (((*this)[idx].ycoords.start <= y) && ((*this)[idx].ycoords.end >= y));
  return (*this)[idx];

}

bool Anchors::resolve_parallel (const size_t max_join_distance /* = 0 */, const bool concatenate_immutable /* = false */) {

  CTAG (6,ANCHORING ANCHORING_VERBOSE) << "Resolving parallel anchors." << endl;

  // catch case of no anchors (very important!)
  if (size() == 0)
    return true;

  // lexical sort in x coordinate (essential)
  // we get this automatically from MUMmer if not using translated anchors,
  // but must do it by hand if using translated anchors
  std::sort (begin(), end(), Anchor::lexical_comparison_x);

  // first deal with overlapping degenerate parallel anchors
  Anchors resolved (this->xseq, this->yseq);
  std::vector<bool> processed (size(), false);
  for (size_t i = 0; i < this->size(); ++i) {

    // if we've already dealt with this anchor, ignore it and go to the next
    if (processed[i])
      continue;
    // else mark as processed
    else
      processed[i] = true;   // this actually isn't necessary but it makes me happy

    Anchor& curr = (*this)[i];

    // don't concatenate immutable anchors unless so requested; leave them disjoint
    // this is crucial when concatenating parallel anchors prior to performing anchor annealing:
    // if two original immutable anchors abut and one is concatenated so
    // as to give a centroid overlapping the centroid of the other immutable anchor,
    // then we will get negative gap posterior errors
    // => very very bad
    // So: just store current (immutable) anchor
    if (curr.is_immutable() && !concatenate_immutable) {
      resolved.store (curr);
      continue;
    }

    // right-hand x boundary of the current anchor
    unsigned currright = curr.xcoords.end;

    // look for overlapping parallel anchors
    for (size_t j = i + 1; j < this->size(); ++j) { // i + 1 because lexical sorting

      const Anchor& overlap = (*this)[j];

      // don't concatenate immutable anchors unless so requested; leave them disjoint
      if (overlap.is_immutable() && !concatenate_immutable)
	continue;

      // does this anchor overlap or is it within max_join_distance of the previous?
      unsigned overlapleft = overlap.xcoords.start;
      unsigned overlapright = overlap.xcoords.end;
      if (static_cast<int> (overlapleft - currright) > static_cast<int> (max_join_distance)) // if not, then stop;
	break;                                         // lexical ordering means there can be no other overlaps or acceptably adjacent anchors

      // is this anchor parallel with the current anchor?
      if (!curr.is_parallel (overlap))  // if not, then go on to the next
	continue;
      else
	processed[j] = true;            // else mark as processed and deal with it

      // if this overlapping anchor is completely degenerate (contained within anchor curr),
      // then just ignore it (conveys no new information) 
      // (although handle special cases of immutable anchors)
      if (overlapright <= currright) {

	// handle case of immutable anchor
	// deal with it by marking the containing anchor as immutable
	// and ignoring (implicitly dropping) the original immutable anchor
	if (overlap.is_immutable() && concatenate_immutable) {
	  curr.set_immutable();
	  // log
	  if (CTAGGING(-1,ANCHORING ANCHORING_VERBOSE)) {
	    CL << "Marking anchor containing immutable anchor as immutable itself:" << endl
	       << curr << endl;
	  }
	}

	// handle case of external_scoring anchor
	// deal with it by marking the containing anchor as external_scoring
	// and ignoring (implicitly dropping) the original external_scoring anchor
	if (overlap.is_external_scoring()) {
	  curr.set_external_scoring();
	  curr.set_score (overlap.get_score()); // won't affect immutable anchors
	  // log
	  if (CTAGGING(-1,ANCHORING ANCHORING_VERBOSE)) {
	    CL << "Marking anchor containing externally-scored anchor as externally-scored itself:" << endl
	       << curr << endl;
	  }
	}

	// else ignore (implicitly drop) the degenerate anchor
	else {
	  // log
	  if (CTAGGING(-1,ANCHORING ANCHORING_VERBOSE)) {
	    CL << "Dropping degenerate first anchor in favor of second:" << endl
	       << overlap << endl
	       << curr << endl;
	  }
	}

	continue;
      }

      // else increment boundaries of current anchor as appropriate
      // (although handle special cases of immutable anchors)
      else {

	// log
	if (CTAGGING(-1,ANCHORING ANCHORING_VERBOSE)) {
	  CL << "Concatenating these two parallel anchors:" << endl
	     << curr << endl
	     << overlap << endl;
	}

	// increment coordinates for current anchor
	curr.xcoords.end = overlap.xcoords.end;
	curr.ycoords.end = overlap.ycoords.end;
	currright = curr.xcoords.end;

	// handle case of immutable anchor
	// deal with it by marking the current anchor as immutable
	if (overlap.is_immutable() && concatenate_immutable) {
	  curr.set_immutable();
	  // log
	  if (CTAGGING(-1,ANCHORING ANCHORING_VERBOSE)) {
	    CL << "Marking anchor containing immutable anchor as immutable itself:" << endl
	       << curr << endl;
	  }
	}

	// handle case of external_scoring anchor
	// deal with it by marking the current anchor as external_scoring
	if (overlap.is_external_scoring()) {
	  curr.set_external_scoring();
	  curr.set_score (overlap.get_score()); // won't affect immutable anchors
	  // log
	  if (CTAGGING(-1,ANCHORING ANCHORING_VERBOSE)) {
	    CL << "Marking anchor containing externally-scored anchor as externally-scored itself:" << endl
	       << curr << endl;
	  }
	}

      }

    }

    // update and store current anchor
    curr.update();
    resolved.store (curr);
    resolved.update_score (resolved.size() - 1);

  }

  // store the results
  __anchors.assign (resolved.begin(), resolved.end());

  // enforce immutable unique
  // (we can't call create_spanning_map directly before enforce_immutable_unique
  // because we need to ensure that the immutable anchors are nondegenerate first)
  enforce_immutable_unique();

  // now re-create the spanning_map
  return create_spanning_map();

}

void Anchors::enforce_immutable_unique() {

  std::set<unsigned> immutable_x;
  std::set<unsigned> immutable_y;

  // get immutable coordinates
  const Anchors immutable_anchors;
  for (std::vector<Anchor>::const_iterator anchor = begin(); anchor != end(); ++anchor) {
    if (anchor->is_immutable()) {
      immutable_x.insert (anchor->x);
      immutable_y.insert (anchor->y);
    }
  }

  // if there are no immutable anchors, then we're done!
  if (!immutable_x.size())
    return;
  
  Anchors passed;
  for (std::vector<Anchor>::const_iterator anchor = begin(); anchor != end(); ++anchor) {

    // if immutable, store and continue    
    if (anchor->is_immutable()) {
      passed.store (*anchor);
      continue;
    }

    // else see if it overlaps an immutable anchor
    else {

      // drop it if overlaps with immutable
      if ((immutable_x.find (anchor->x) != immutable_x.end()) || (immutable_y.find (anchor->y) != immutable_y.end())) {
	if (CTAGGING(-1,ANCHORING ANCHORING_VERBOSE)) {
	  CL << "Dropping anchor which overlaps an immutable anchor:" << endl
	     << *anchor << endl;
	}
	continue;
      }

      // else keep it
      else {
	passed.store (*anchor);
      }

    }

  }

  __anchors.assign (passed.begin(), passed.end());

}

void Anchors::remove_overlaps (const size_t minlen /* = ANCHOR_NUC_MINLEN_DEFAULT */) {

  // the x and y coordinates are handled separately in the hope of greater efficiency
  // (since pruning x means pruning y as well => fewer possiblities for overlapping anchors)

  // sort by x, prune by x
  std::sort (begin(), end(), Anchor::lexical_comparison_x);
  remove_overlaps (minlen, 0);

  // sort by x, prune by y
  std::sort (begin(), end(), Anchor::lexical_comparison_x);
  remove_overlaps (minlen, 1);

  // sort by y, prune by x
  std::sort (begin(), end(), Anchor::lexical_comparison_y);
  remove_overlaps (minlen, 0);

  // sort by y, prune by y
  std::sort (begin(), end(), Anchor::lexical_comparison_y);
  remove_overlaps (minlen, 1);

}

void Anchors::remove_overlaps (const size_t minlen, const unsigned which) {

  // loop over all anchors:
  // for each anchor i
  //   find all overlapping anchors
  //   for each overlapping anchor j
  //      if anchor j scores below the current anchor i, then prune the left-hand boundary of j
  //      if anchor j scores above the current anchor i, then prune the right-hand boundary of i; break
  //       (break because the lexical sorting means that anchor j will be the lefthand-most overlapping anchor)

  // We require both lexical ordering and consistency of anchors.
  // Lexical ordering guarantees that the anchors are sorted by increasing left-hand x coordinate.

  std::vector<bool> is_dead (size(), false);
  for (size_t i = 0; i < size(); ++i) {

    // skip if dead
    if (is_dead[i])
      continue;

    Anchor& curr = (*this)[i]; // must be a reference so that we aren't just modifying a copy!

    // right-hand x boundary of the current anchor
    unsigned xright = (which == 0) ? curr.xcoords.end : curr.ycoords.end;

     // assemble list of overlapping anchors
    std::list<size_t> xoverlaps;
    unsigned xleftmost = xright;  // leftmost boundary of overlaps
    for (size_t j = i + 1; j < size(); ++j) { // i + 1 because lexical sorting

      // ignore dead anchors
      if (is_dead[j])
	continue;

      // does this anchor overlap with the previous?
      unsigned jxleft = (which == 0) ? (*this)[j].xcoords.start : (*this)[j].ycoords.start;
      if (jxleft > xright) // if it doesn't overlap, then stop;
	break;             // lexical ordering means there can be no other overlaps

      // store and increment boundaries
      xoverlaps.push_back (j);
      if (jxleft < xleftmost)
	xleftmost = jxleft;

    }

    // if there are overlaps, make sure that coordinate boundaries are sane
    if (xoverlaps.size())
      assert (xleftmost <= xright);
    // if no overlaps, go to the next anchor
    else
      continue;

    // go through the list of overlapping coordinates:
    //  at each overlapping anchor, choose in favor of either current anchor i or overlapping anchor *overlap
    //  prune anchors as appropriate
    assert (xright >= xleftmost);
    for (unsigned s = 0; s <= xright - xleftmost; ++s) { // s is the distance "within" anchor i
      unsigned x = xleftmost + s;                        // x is the actual x coordinate
      
      assert (x <= ((which == 0) ? curr.xcoords.end : curr.ycoords.end));

      std::list<size_t>::iterator overlapidx = xoverlaps.begin(); // use STL erase-remove idiom
      while (overlapidx != xoverlaps.end()) { // *overlapidx is index of overlapping anchor in (*this)

	Anchor& overlap = (*this)[*overlapidx]; // must be a reference so that we aren't just modifying a copy!

	// does this anchor overlap at coordinate x?
	// if not, then go on to the next
	if (((which == 0) && overlap.xcoords.start > x) || (!(which == 0) && overlap.ycoords.start > x)) {
	  ++overlapidx;
	  continue;
	}

	// log before
	if (CTAGGING(4,ANCHORING ANCHORING_VERBOSE)) {
	  if (which == 0)
	    CL << "Overlapping anchors (x = " << x << "):" << endl
	       << curr << endl
	       << overlap << endl;
	  else
	    CL << "Overlapping anchors (y = " << x << "):" << endl
	       << curr << endl
	       << overlap << endl;

	}

	// if anchor overlap is a new best, 
	// then prune the right-hand boundary of current anchor i here
	// and end loop over coordinates (because we're done with anchor i)
	if (overlap.get_score() > curr.get_score()) {

	  // calculate distance to prune anchor
#ifndef NDEBUG
	  if (which == 0)
	    assert (static_cast<int> (curr.xcoords.end) - (x - 1) >= 0);
	  else
	    assert (static_cast<int> (curr.ycoords.end) - (x - 1) >= 0);
#endif // NDEBUG
	  const unsigned delta = ((which == 0) ? curr.xcoords.end : curr.ycoords.end) - (x - 1);

	  // catch case of current anchor being shorter than the distance to prune
	  // or new length shorter than minlen
	  // (note that this always includes the case of length-1 (transitive) anchors)
	  // if so, then mark it as dead
	  if ((curr.length <= delta) || ((curr.length - delta) < minlen)) {
	    assert (!curr.is_immutable() || overlap.is_immutable()); // sanity check: either the current anchor must not be immutable,
	    is_dead[i] = true;                                       // or the better overlapping anchor must be
	    if (CTAGGING(4,ANCHORING ANCHORING_VERBOSE))
	      CL << "Pruned by removing first anchor." << endl;
	  }
	  // else prune the current anchor *iter
	  else {
	    curr.xcoords.end -= delta;
	    curr.ycoords.end -= delta;
	    curr.update();
	    update_score (i);
	    // log after
	    if (CTAGGING(4,ANCHORING ANCHORING_VERBOSE)) {
	      CL << "Pruned as:" << endl
		 << curr << endl
		 << overlap << endl;
	    }
	  }
	  // if we found a better overlapping anchor, then stop the loop over overlaps
	  // as well as coordinates
	  // (because we've already pruned and are done with the current anchor)
	  goto FOUND_BETTER;

	}

	// otherwise overlapping anchor overlap isn't as good,
	// so prune its left-hand boundary and drop it from the list of overlaps
	else {

	  // calculate distance to prune anchor
#ifndef NDEBUG
	  if (which == 0)
	    assert (static_cast<int> (curr.xcoords.end + 1) - overlap.xcoords.start);
	  else
	    assert (static_cast<int> (curr.ycoords.end + 1) - overlap.ycoords.start);
#endif // NDEBUG
	  const unsigned delta = (which == 0) ? ((curr.xcoords.end + 1) - overlap.xcoords.start) : ((curr.ycoords.end + 1) - overlap.ycoords.start);

	  // catch case of *overlap being shorter than the distance to prune
	  // or new length shorter than minlen
	  // (note that this always includes the case of length-1 (transitive) anchors)
	  // if so, then erase *overlap entirely
	  if ((overlap.length <= delta) || ((overlap.length - delta) < minlen)) {
	    assert (!overlap.is_immutable() || curr.is_immutable()); // sanity check: either the current anchor must not be immutable,
	    is_dead[*overlapidx] = true;                             // or the better overlapping anchor must be
	    if (CTAGGING(4,ANCHORING ANCHORING_VERBOSE))
	      CL << "Pruned by removing second anchor." << endl;
	  }
	  // else prune overlap
	  else {
	    overlap.ycoords.start += delta;
	    overlap.xcoords.start += delta;
	    overlap.update();
	    update_score (*overlapidx);
	    // log after
	    if (CTAGGING(4,ANCHORING ANCHORING_VERBOSE)) {
	      CL << "Pruned as:" << endl
		 << curr << endl
		 << overlap << endl;
	    }
	  }
	  // we're done with anchor *overlap; drop it from the list
	  overlapidx = xoverlaps.erase (overlapidx);

	}

	// Note that we don't need to explicitly increment overlapidx here at all:
	// If *overlap is better than the curren anchor *iter, then we break out of
	// the overlapidx loop with an explicit goto FOUND_BETTER.
	// If *overlap is worse, then we prune it, remove it from the list of overlaps
	// and reset overlapidx when using erase.

      } // end loop over overlapping anchors

    } // end loop over coordinates

  FOUND_BETTER:
    ; // hack to get around requirement of a statement post-label

  }

  // assemble and store list of pruned anchors
  Anchors pruned (this->xseq, this->yseq);
  for (size_t i = 0; i < size(); ++i) {
    if (!is_dead[i])
      pruned.store ((*this)[i]);
  }

  // store list of pruned anchors
  __anchors.assign (pruned.begin(), pruned.end());

}


// to do: fix this for proper scoring
// pull out the anchored subseqs, get hamming distance (mismatches)
// and use BLAST theory / forward approximation
void Anchors::update_score (const size_t idx) {

  assert ((idx >= 0) && (idx < size()));

  Anchor& anchor = (*this)[idx];

  // compute anchor score as (1 - p_value)
  if (!anchor.is_external_scoring()) {
    anchor.set_score (1.0 - anchor.p_value (this->xseq, this->yseq));
    assert ((anchor.get_score() >= 0.0) && (anchor.get_score() <= 1.0));
  }

}

void Anchors::impose_normalized_probability_distribution() {

  // sum of scores of anchors with particular coordinates
  std::vector<double> sum_scorex (xseq->length(), 0.);
  std::vector<double> sum_scorey (yseq->length(), 0.);

  // max of scores of anchors with particular coordinates
  std::vector<double> max_scorex (xseq->length(), 0.);
  std::vector<double> max_scorey (yseq->length(), 0.);

  // accumulate counts for normalization
  for (size_t i = 0; i < size(); ++i) {

    // update score
    update_score (i);

    const Anchor& anchor = (*this)[i];

    // check sane
    assert ((anchor.get_score() >= 0.0) && (anchor.get_score() <= 1.0));

    // accumulate counts for normalization
    sum_scorex[anchor.x] += anchor.get_score();
    sum_scorey[anchor.y] += anchor.get_score();

    max_scorex[anchor.x] = max (anchor.get_score(), max_scorex[anchor.x]);
    max_scorey[anchor.y] = max (anchor.get_score(), max_scorey[anchor.y]);

  }

  for (size_t i = 0; i < size(); ++i) {

    Anchor& anchor = (*this)[i];

    if (anchor.is_immutable())
      continue;

    // sanity check for finiteness
    assert ((sum_scorex[anchor.x] > 0) && (sum_scorey[anchor.y] > 0));

    // reweight scores
    anchor.set_score (max_scorex[anchor.x] * (anchor.get_score() / sum_scorex[anchor.x]));
    anchor.set_score (max_scorey[anchor.y] * (anchor.get_score() / sum_scorey[anchor.y]));

    // check sane
    assert ((anchor.get_score() >= 0.0) && (anchor.get_score() <= 1.0));

  }

}

Anchors Cabbage_adapter::get_candidate_anchors() const {

  CTAG (6,ANCHORING ANCHORING_VERBOSE) << "Getting candidate anchors for sequence pair '" << xseq.name << "' and '" << yseq.name << "'." << endl;

  Anchors anchors (xseq.seq, yseq.seq);
  
  // (90,110)  -- (190,210) => 100 -- 200; 0.1
  // (140,160) -- (170,190) => 150 -- 180; 0.2
  // (240,260) -- (290,310) => 250 -- 300; 0.3
  // Because we use a 0-based coordinate system,
  // these anchor the first "different" nucleotides starting a new sequence
  // (rather than the last "same" ending the old sequence).

  Anchor one (Interval (90, 110), Interval (190, 210), .1);
  Anchor two (Interval (140, 160), Interval (170, 190), .2);
  Anchor three (Interval (240, 260), Interval (290, 310), .3);

  anchors.store (one);
  anchors.store (two);
  anchors.store (three);

  if (!anchors.create_spanning_map())
    THROWEXPR ("ERROR: Degenerate anchors.");

  return anchors;
}

// Wrapper for call_exonerate.
Anchors Exonerate_adapter::get_candidate_anchors (const int minscore /* = EXONERATE_MINSCORE_DEFAULT */) const {

  CTAG (8,ANCHORING ANCHORING_VERBOSE) << "Getting candidate anchors (with exonerate) for sequence pair '" << xseq.name << "' and '" << yseq.name << "'." << endl;

  Anchors anchors = call_exonerate (xseq.seq, yseq.seq, minscore);
  return anchors;

}

Anchors Exonerate_adapter::call_exonerate (const std::string& xseq, const std::string& yseq, const int minscore /* = EXONERATE_MINSCORE_DEFAULT */) const {

  // check that exonerate executable is properly defined;
  // if not, then die
#ifndef EXONERATE_EXEC
  THROWEXPR ("ERROR: exonerate executable isn't defined; try recompiling and specifying the exonerate path.");
  // else call exonerate!
#else

  // get pid of current process for unique filenames
  pid_t curr_pid = getpid();

  // write sequence files to scratch directory
  std::string xfile = std::string (TMP_DIR) + "/" + FILENAME_PREFIX + "_" + Util::to_string (curr_pid) + "_x.fa";
  std::string yfile = std::string (TMP_DIR) + "/" + FILENAME_PREFIX + "_" + Util::to_string (curr_pid) + "_y.fa";
  // X
  std::ofstream filestream;
  filestream.open (xfile.c_str());
  if (!filestream.is_open())
      THROWEXPR ("ERROR: Couldn't create file with name '" << xfile << "'.");
  filestream << ">X" << endl;
  filestream << xseq << endl;
  filestream.close();
  // Y
  filestream.open (yfile.c_str());
  if (!filestream.is_open())
      THROWEXPR ("ERROR: Couldn't create file with name '" << yfile << "'.");
  filestream << ">Y" << endl;
  filestream << yseq << endl;
  filestream.close();

  // create exonerate executable string
  // prepare for use with execv()
  // (which requires null-terminated string)
  std::string exec;
  exec = std::string (EXONERATE_EXEC)
    + " " + xfile + " " + yfile
    + " --querytype dna --targettype dna "
    + " --score " + Util::to_string (minscore) + " --gappedextension false "
    + " --saturatethreshold 5 "
    + " --bigseq true "
    + " --showvulgar false --showalignment false ";
  if (__softmasked)
    exec += " --softmaskquery true --softmasktarget true ";
  if (__use_translated)
    exec += " --model ungapped:trans --proteinsubmat blosum62 ";
  else
    exec += " --model ungapped ";

  std::vector<std::string> exonerate_args = Util::split (exec);
  exonerate_args.push_back (static_cast<std::string> ("--ryo")); // this hack is necessary to prevent splitting the --ryo args
  exonerate_args.push_back (static_cast<std::string> ("%qab %qae %qS %tab %tae %tS %s\\n")); // note that this is NOT enclosed in double quotes
                                                                                         // these are only necessary for command-line usage
                                                                                         // (to prevent the shell from splitting the argument)

  char* exonerate_argv [exonerate_args.size() + 1];
  for (size_t i = 0; i < exonerate_args.size(); ++i)
    exonerate_argv[i] = static_cast<char*> (const_cast<char*> (exonerate_args[i].c_str()));
  exonerate_argv[exonerate_args.size()] = static_cast<char*> (0);
  if (CTAGGING(6,ANCHORING ANCHORING_VERBOSE)) {
    std::string exec_line;
    for (std::vector<std::string>::iterator iter = exonerate_args.begin(); iter != exonerate_args.end(); ++iter)
      exec_line += ' ' + *iter;
    CTAG (6,ANCHORING ANCHORING_VERBOSE) << "Calling exonerate as '" << exec_line << "'" << endl;
  }

  // run exonerate
  pid_t exonerate_pid;

  // clear cin
  std::cin.clear();

  // tmp file to hold exonerate output
  std::string output_file = std::string (TMP_DIR) + "/" + FILENAME_PREFIX + "_" + Util::to_string (curr_pid) + "_output";

  // fork and call exonerate
  int exonerate_status = 0;
  if ((exonerate_pid = fork()) < 0) {  // a negative PID indicates failure
    THROWEXPR ("ERROR: Couldn't start fork; exiting.");
  }
  else if (exonerate_pid == 0) {       // a 0 PID indicates the child process
    freopen (output_file.c_str(), "w", stdout);    // capture stdout (send to file output_file) to store exonerate results
    freopen ("/dev/null", "a", stderr);            // capture stderr (send to trash) so that it doesn't mess up our logging
    if (execv (exonerate_argv[0], exonerate_argv) == -1) // replace child fork with a new process
      THROWEXPR ("ERROR: Couldn't run exonerate (execv problem).");
    fclose (stdout);
    fclose (stderr);
  }
  else {                               // a positive PID indicates the parent process
    waitpid (exonerate_pid, &exonerate_status, 0);  // wait for child process to end
    if (exonerate_status)                           // did it work?
      THROWEXPR ("ERROR: Running exonerate failed with exit code " << exonerate_status << ".");
  }

  // read output_file and parse results to get anchors
  // exonerate output format is:
  //    564000 563704 - 243683 243386 - 168
  //    xstart xend xstrand ystart yend ystrand raw_score
  // NB exonerate uses 0-based, half-open coordinates [start, end);
  // we convert from this format to our 0-based, closed interval coordinates
  Anchors anchors (xseq, yseq);
  std::ifstream output_stream (output_file.c_str(), std::ifstream::in);
  std::string line;
  if (CTAGGING(-1,ANCHORING ANCHORING_VERBOSE_VERBOSE))
    CL << "exonerate output:" << endl;

  while (!output_stream.eof()) {

    getline (output_stream, line);
    if (line.length() == 0) // skip empty lines
      continue;

    // log
    if (CTAGGING(-1,ANCHORING ANCHORING_VERBOSE_VERBOSE))
      CL << line;
    std::string buffer;
    std::stringstream ss (line);
    std::vector<std::string> tokens; // vector to hold whitespace-separated tokens (words)
    while (ss >> buffer)
      tokens.push_back (buffer);

    // skip exonerate's progress messages
    if ((tokens.size() > 0) && (tokens[0] == "Command") || (tokens[0] == "Hostname:") || (tokens[0] == "Message:") || (tokens[0] == "--"))
      continue;

    // now parse tokens into anchors
    if (tokens.size() == 7) {

      // parse output
      const unsigned xstart = atoi (tokens[0].c_str());   // 0-based, half-open coordinates [start, end)
      const unsigned xend = atoi (tokens[1].c_str()) - 1;
      assert (tokens[2].length() == 1);        // assert single character indicating strand
      const char xstrand = (tokens[2])[0];
      const unsigned ystart = atoi (tokens[3].c_str());
      const unsigned yend = atoi (tokens[4].c_str()) - 1;
      assert (tokens[5].length() == 1);        // assert single character indicating strand
      const char ystrand = (tokens[5])[0];
      const double raw_score = static_cast<double> (atoi (tokens[6].c_str()));

      // only consider collinear hits on forward strand
      if ((xstrand != '+') || (ystrand != '+'))
	continue;

      // check sane
      assert (xstart <= xend);
      assert (ystart <= yend);
      assert (raw_score >= minscore - DOUBLE_TINY);

      // create anchor for ungapped aligned intervalval
      anchors.store (Anchor (Interval (xstart, xend),
			     Interval (ystart, yend),
			     raw_score));

      // mark the anchor as scored externally (no score updating later on)
      anchors.back().set_external_scoring();

    } else if (tokens.size() == 0) {
      continue;
    } else {
      if (CTAGGING(4,ANCHORING ANCHORING_VERBOSE))
	CL << "WARNING: I don't know how to parse exonerate output line: " << line;
      continue;
    }
  }

  output_stream.close();

  // clean up temporary files
  if (remove (xfile.c_str()) != 0)
    CTAG (9,ANCHORING ANCHORING_VERBOSE) << "WARNING: Couldn't delete temporary file '" << xfile << "'." << endl;
  if (remove (yfile.c_str()) != 0)
    CTAG (9,ANCHORING ANCHORING_VERBOSE) << "WARNING: Couldn't delete temporary file '" << yfile << "'." << endl;
  if (remove (output_file.c_str()) != 0)
    CTAG (9,ANCHORING ANCHORING_VERBOSE) << "WARNING: Couldn't delete temporary file '" << output_file << "'." << endl;

  return anchors;

#endif /* EXONERATE_EXEC */

}


// Wrapper for call_mummer.
Anchors Mummer_adapter::get_candidate_anchors (size_t minlen /* = ANCHOR_NUC_MINLEN_DEFAULT */) const {

  CTAG (8,ANCHORING ANCHORING_VERBOSE) << "Getting candidate anchors (with MUMmer) for sequence pair '" << xseq.name << "' and '" << yseq.name << "'." << endl;

  Anchors anchors (xseq.seq, yseq.seq);

  // anchoring on translated sequences in protein space
  if (__use_translated) {

    // get translated sequences
    const Translated_sequence tr_x (xseq);
    const Translated_sequence tr_y (yseq);
    
    // forward strand for all frame combinations
    for (size_t fx = 0; fx < 3; ++fx) {
      for (size_t fy = 0; fy < 3; ++fy) {

	// get anchors
	Anchors tr_anchors = call_mummer (tr_x.get_forward (fx).seq, tr_y.get_forward (fy).seq, minlen);

	// map anchor coordinates back to nucleotide sequence
	for (std::vector<Anchor>::iterator anchor = tr_anchors.begin(); anchor != tr_anchors.end(); ++anchor) {
	  anchor->xcoords = tr_x.map_interval_to_orig (true, fx, anchor->xcoords);
	  anchor->ycoords = tr_y.map_interval_to_orig (true, fy, anchor->ycoords);
	  anchor->update();
	  // sanity check on mapping
	  assert (anchor->xcoords.end < xseq.length());
	  assert (anchor->ycoords.end < yseq.length());
	  // store anchor
	  anchors.store (*anchor);
	  anchors.update_score (anchors.size() - 1);
	}

      }
    }

//    // reverse strand for all frame combinations
//    for (size_t fx = 0; fx < 3; ++fx) {
//      for (size_t fy = 0; fy < 3; ++fy) {
//
//	// get anchors
//	Anchors tr_anchors = call_mummer (tr_x.reverse[fx], tr_y.reverse[fy], minlen);
//
//	// map anchor coordinates back to nucleotide sequence
//	for (std::vector<Anchor>::iterator anchor = tr_anchors.begin(); anchor != tr_anchors.end(); ++anchor) {
//	  anchor->xcoords = tr_x.map_interval_to_orig (false, fx, anchor->xcoords);
//	  anchor->ycoords = tr_y.map_interval_to_orig (false, fy, anchor->ycoords);
//	  anchor->update();
//	  // store anchor
//	  anchors.store (*anchor);
//	  anchors.update_score (anchors.size() - 1);
//	}
//
//      }
//    }
// reverse strand doesn't convey much more information
//   -- RKB 10/11/08

  }

  // anchoring on nucleotide sequence
  else {
    anchors = call_mummer (xseq.seq, yseq.seq, minlen);
  }

  return anchors;

}

// I originally wrote this function using pipes to redirect the child process (mummer) stdout
// to the parent process stdin.  For some inscrutable reason, it only worked the first
// time it was called; subsequent calls had nothing in cin.
// This was very weird; I was saving and restoring the original stdin and stdout, but
// it happened anyways.
// Hence the current way, where I redirect stdout of the child process to a file.
Anchors Mummer_adapter::call_mummer (const std::string& xseq, const std::string& yseq, size_t minlen /* = ANCHOR_NUC_MINLEN_DEFAULT */) const {

  // check that MUMmer executable is properly defined;
  // if not, then die
#ifndef MUMMER_EXEC
  THROWEXPR ("ERROR: MUMmer executable isn't defined; try recompiling and specifying the MUMmer path.");
  // else call MUMmer!
#else

  // get pid of current process for unique filenames
  const pid_t curr_pid = getpid();

  // write sequence files to scratch directory
  std::string xfile = std::string (TMP_DIR) + "/" + FILENAME_PREFIX + "_" + Util::to_string (curr_pid) + "_x.fa";
  std::string yfile = std::string (TMP_DIR) + "/" + FILENAME_PREFIX + "_" + Util::to_string (curr_pid) + "_y.fa";
  // X
  std::ofstream filestream;
  filestream.open (xfile.c_str());
  if (!filestream.is_open())
      THROWEXPR ("ERROR: Couldn't create file with name '" << xfile << "'.");
  filestream << ">X" << endl;
  filestream << xseq << endl;
  filestream.close();
  // Y
  filestream.open (yfile.c_str());
  if (!filestream.is_open())
      THROWEXPR ("ERROR: Couldn't create file with name '" << yfile << "'.");
  filestream << ">Y" << endl;
  filestream << yseq << endl;
  filestream.close();

  // create mummer executable string
  // prepare for use with execv()
  // (which requires null-terminated string)
  std::string exec = std::string (MUMMER_EXEC)
    + " -mum -l " + Util::to_string (minlen)
    + " " + xfile + " " + yfile;
  const std::vector<std::string> mummer_args = Util::split (exec);
  char* mummer_argv [mummer_args.size() + 1];
  for (size_t i = 0; i < mummer_args.size(); ++i)
    mummer_argv[i] = static_cast<char*> (const_cast<char*> (mummer_args[i].c_str()));
  mummer_argv[mummer_args.size()] = static_cast<char*> (0);
  if (CTAGGING(6,ANCHORING ANCHORING_VERBOSE))
    CTAG (6,ANCHORING ANCHORING_VERBOSE) << "Calling MUMmer as '" << exec << "'" << endl;
  
  // run mummer
  pid_t mummer_pid;

  // clear cin
  std::cin.clear();

  // tmp file to hold mummer output
  std::string output_file = std::string (TMP_DIR) + "/" + FILENAME_PREFIX + "_" + Util::to_string (curr_pid) + "_output";

  // fork and call mummer
  int mummer_status = 0;
  if ((mummer_pid = fork()) < 0) {  // a negative PID indicates failure
    THROWEXPR ("ERROR: Couldn't start fork; exiting.");
  }
  else if (mummer_pid == 0) {       // a 0 PID indicates the child process
    freopen (output_file.c_str(), "w", stdout);    // capture stdout (send to file output_file) to store mummer results
    freopen ("/dev/null", "a", stderr);            // capture stderr (send to trash) so that it doesn't mess up our logging
    if (execv (mummer_argv[0], mummer_argv) == -1) // replace child fork with a new process
      THROWEXPR ("ERROR: Couldn't run MUMmer (execv problem).");
    fclose (stdout);
    fclose (stderr);
  }
  else {                            // a positive PID indicates the parent process
    waitpid (mummer_pid, &mummer_status, 0);  // wait for child process to end
    if (mummer_status)                        // did it work?
      THROWEXPR ("ERROR: Running MUMmer failed with exit code " << mummer_status << ".");
  }

  // read output_file and parse results to get anchors
  // mummer output format is:
  // >seqname
  //             xpos                         ypos                      match_length
  // (position in reference sequence /t position in query sequence /t length of match)
  // Note that we convert from 1-based (MUMmer's output) to 0-based coordinates here
  Anchors anchors (xseq, yseq);
  std::ifstream output_stream (output_file.c_str(), std::ifstream::in);
  std::string line;
  if (CTAGGING(-1,ANCHORING ANCHORING_VERBOSE_VERBOSE))
    CL << "MUMmer output:" << endl;
  while (!output_stream.eof()) {
    getline (output_stream, line);
    // log
    if (CTAGGING(-1,ANCHORING ANCHORING_VERBOSE_VERBOSE))
      CL << line;
    std::string buffer;
    // skip comments and sequence names
    if ((line.length() > 0) && ((line[0] == '#') || (line[0] == '>')))
      continue;
    std::stringstream ss (line);
    std::vector<std::string> tokens; // vector to hold whitespace-separated tokens (words)
    while (ss >> buffer)
      tokens.push_back (buffer);

    // now parse tokens into anchors
    if (tokens.size() == 3) {
      unsigned xstart = atoi (tokens[0].c_str()) - 1; // convert from 1-based to 0-based coordinates
      unsigned ystart = atoi (tokens[1].c_str()) - 1;
      size_t length = atoi (tokens[2].c_str());
      assert (length >= minlen);
      anchors.store (Anchor (Interval (xstart, xstart + length - 1),
			     Interval (ystart, ystart + length - 1))); // -1 because closed interval
      anchors.update_score (anchors.size() - 1);
    } else if (tokens.size() == 0) {
      continue;
    } else {
      if (CTAGGING(4,ANCHORING ANCHORING_VERBOSE))
	CL << "WARNING: I don't know how to parse MUMmer output line: " << line;
      continue;
    }
  }

  output_stream.close();
  
  // clean up temporary files
  if (remove (xfile.c_str()) != 0)
    CTAG (9,ANCHORING ANCHORING_VERBOSE) << "WARNING: Couldn't delete temporary file '" << xfile << "'." << endl;
  if (remove (yfile.c_str()) != 0)
    CTAG (9,ANCHORING ANCHORING_VERBOSE) << "WARNING: Couldn't delete temporary file '" << yfile << "'." << endl;
  if (remove (output_file.c_str()) != 0)
    CTAG (9,ANCHORING ANCHORING_VERBOSE) << "WARNING: Couldn't delete temporary file '" << output_file << "'." << endl;

  return anchors;

#endif /* MUMMER_EXEC */

}

void Anchor_resolver::add_mercator_constraints (const std::string& filename) {

  // read Mercator file
  __constraints_set.read_mercator (filename);
 
}

const std::vector<Anchors> Anchor_resolver::get_resolved_anchors (const size_t minlen /* = ANCHOR_NUC_MINLEN_DEFAULT */, const size_t max_join_distance /* = 0 */, bool const use_translated /* = false */,
								  const bool use_exonerate /* = false */, const int minscore /* = EXONERATE_MINSCORE_DEFAULT */, const bool softmasked /* = false */,
								  const bool hardmasked /* = true */,
								  const size_t num_refinement_steps /* = 0 */,
								  const bool output_for_gui /* = false */, const std::string gui_prefix /* = "" */) const {

  // initialize dummy Tree_weights (returns all weights as 1.0)
  Tree_weights tree_weights;
  return get_resolved_anchors (tree_weights,
			       minlen, max_join_distance, use_translated,
			       use_exonerate, minscore, softmasked,
			       hardmasked,
			       num_refinement_steps,
			       output_for_gui, gui_prefix);

}

const std::vector<Anchors> Anchor_resolver::get_resolved_anchors (const Tree_weights& tree_weights,
								  const size_t minlen /* = ANCHOR_NUC_MINLEN_DEFAULT */, const size_t max_join_distance /* = 0 */, const bool use_translated /* false */,
								  const bool use_exonerate /* = false */, const int minscore /* = EXONERATE_MINSCORE_DEFAULT */, const bool softmasked /* = false */,
								  const bool hardmasked /* = true */,
								  const size_t num_refinement_steps /* = 0 */,
								  const bool output_for_gui /* = false */, std::string gui_prefix /* = "" */) const {

  // store all candidate anchors for later recovery
  std::vector<Anchors> candidate_anchors_list (seq_pairs.size());

  // This proceeds as:
  // - get candidate anchors for all sequence pairs we're considering
  // - create and store Post_probs objects for each sequence pair
  // - create "reduced" sequences, where a reduced sequence is obtained from
  //    the original sequence by taking only positions corresponding to the
  //    centroids of candidate anchors
  // - perform anchor annealing on the reduced sequences
  // - map the resolved anchors back to the original sequence coordinates

  // information to build reduced sequences:
  //  set of positions within the sequence which correspond to candidate anchors
  std::vector<std::set<unsigned> > anchored_positions (seq_db.size());

  // get list of candidate anchors
  // list of Post_probs for each sequence pair
  std::vector<Post_probs> post_probs_list (seq_pairs.size());
  for (size_t cnt = 0; cnt < seq_pairs.size(); ++cnt) {

    // get sequence indices for this sequence pair
    const size_t i = seq_pairs[cnt].first;
    const size_t j = seq_pairs[cnt].second;

    const Sequence& xseq = seq_db_internal.get_seq (i);
    const Sequence& yseq = seq_db_internal.get_seq (j);

    // if one of the seqs is empty, don't try to find anchors
    // (MUMmer doesn't handle empty sequences gracefully)
    if (!xseq.length() || !yseq.length())
      continue;

    // get candidate anchors
    const Mummer_adapter mummer_adapter (xseq, yseq, use_translated);            // find in protein space if requested
    candidate_anchors_list[cnt] = mummer_adapter.get_candidate_anchors (minlen);

    // if present, add Mercator constraint information as anchors
    if (__constraints_set.size()) {
      const Sequence& xseq_orig = seq_db.get_seq (i);
      const Sequence& yseq_orig = seq_db.get_seq (j);
      const Anchors constraints_anchors = Anchors::convert_constraints_to_anchors (__constraints_set.get_constraints (i, j), xseq_orig, yseq_orig, hardmasked);
      // store constraints as anchors
      candidate_anchors_list[cnt].store (constraints_anchors);
    }

    // get exonerate anchors if so requested
    if (use_exonerate) {
      const Exonerate_adapter exonerate_adapter (xseq, yseq, use_translated, softmasked);
      const Anchors exonerate_anchors = exonerate_adapter.get_candidate_anchors (minscore);
      candidate_anchors_list[cnt].store (exonerate_anchors);
    }

    // enforce immutable unique
    // important that we do this before any other manipulations
    // (otherwise, e.g., resolve_parallel might fail)
    candidate_anchors_list[cnt].enforce_immutable_unique();

    // resolve parallel anchors:
    // remove degenerate anchors, merge overlapping anchors and concatenate adjacent anchors if requested
    // essential if using translated anchors
    if (!candidate_anchors_list[cnt].resolve_parallel (max_join_distance, false)) // don't concatenate immutable anchors
      THROWEXPR ("ERROR: Couldn't resolve parallel anchors.");

    // record the centroids of the candidate anchors to later build the reduced sequences
    for (std::vector<Anchor>::const_iterator anchor = candidate_anchors_list[cnt].begin(); anchor != candidate_anchors_list[cnt].end(); ++anchor) {
      anchored_positions[i].insert (anchor->x);
      anchored_positions[j].insert (anchor->y);
    }

  }

  // now properly normalize the scores of externally-scored anchors:
  // first get the maximum score for an externally-scored anchor
  // over the anchors for /all/ sequence pairs
  double max_external_score = 0;
  for (size_t cnt = 0; cnt < seq_pairs.size(); ++cnt) {
    for (std::vector<Anchor>::const_iterator anchor = candidate_anchors_list[cnt].begin(); anchor != candidate_anchors_list[cnt].end(); ++anchor) {
      if (anchor->is_external_scoring() && (anchor->get_score() > max_external_score))
	max_external_score = anchor->get_score();
    }
  }

  // then use this max score as a normalization constant
  for (size_t cnt = 0; cnt < seq_pairs.size(); ++cnt) {
    for (std::vector<Anchor>::iterator anchor = candidate_anchors_list[cnt].begin(); anchor != candidate_anchors_list[cnt].end(); ++anchor) {
      if (anchor->is_external_scoring())
	anchor->set_score ((anchor->get_score() / max_external_score) - DOUBLE_TINY); // (ensure all scores < 1.0)
      // check sane after normalization to [0,1]
      assert ((anchor->get_score() >= 0.0) && (anchor->get_score() <= 1.0));
    }
  }

  // now finish pre-processing the anchors prior to anchor annealing and log them
  for (size_t cnt = 0; cnt < seq_pairs.size(); ++cnt) {

    // convert score distribution over anchors to a normalized probability distribution
    candidate_anchors_list[cnt].impose_normalized_probability_distribution();

    // convert candidate anchors to Post_probs format
    post_probs_list[cnt] = candidate_anchors_list[cnt].convert_to_post_probs();

    // get sequence indices for this sequence pair
    const size_t i = seq_pairs[cnt].first;
    const size_t j = seq_pairs[cnt].second;

    const Sequence& xseq = seq_db_internal.get_seq (i);
    const Sequence& yseq = seq_db_internal.get_seq (j);

    // log candidate anchors
    if (CTAGGING(4,ANCHORING ANCHORING_VERBOSE)) {
      CL << "Sequences '" << xseq.name << "' and '" << yseq.name << "' have candidate anchors:" << endl;
      candidate_anchors_list[cnt].show (CL);
    }

    // log
    size_t num_anchors = candidate_anchors_list[cnt].size();
    CTAG (7,ANCHORING ANCHORING_VERBOSE) << "Found " << num_anchors << " candidate anchors for sequence pair '" << xseq.name << "' and '" << yseq.name << "'." << endl;

  }

  // create reduced_seq_db to hold actual reduced sequences
  Sequence_database reduced_seq_db;
  for (size_t i = 0; i < seq_db.size(); ++i) {
    const Sequence& sequence = seq_db.get_seq (i);
    std::string seq (anchored_positions[i].size(), 'R'); // R for Rob!
    reduced_seq_db.add_seq (Sequence (sequence.name, seq));
  }

  // create map from positions in "reduced" sequences to positions in actual sequences
  // note that the anchor positions in anchored_positions[i] are already sorted because it's a set
  std::vector<std::map<unsigned, unsigned> > reduced_position_map (seq_db.size());            // map reduced -> original
  std::vector<std::map<unsigned, unsigned> > reduced_position_reverse_map (seq_db.size());    // map original -> reduced
  if (CTAGGING(-1,ANCHORING_VERBOSE))
    CL << "reduced_position_map:" << endl;
  for (size_t i = 0; i < seq_db.size(); ++i) {
    unsigned ii = 0;
    for (std::set<unsigned>::iterator x = anchored_positions[i].begin(); x != anchored_positions[i].end(); ++x, ++ii) { // ii is the position within the reduced sequence
      reduced_position_map[i].insert (std::make_pair (ii, *x));                                                                // x is the position within the original sequence
      reduced_position_reverse_map[i].insert (std::make_pair (*x, ii));
      if (CTAGGING(-1,ANCHORING_VERBOSE))
	CL << "(" << i << ", " << *x << ") == " << ii << endl;
    }
  }

  // store sparse matrices on reduced sequences for anchor annealing
  std::vector<std::vector<SparseMatrix*> > sparse_matrices (seq_db.size(), std::vector<SparseMatrix*> (seq_db.size(), reinterpret_cast<SparseMatrix*> (NULL)));
  for (size_t cnt = 0; cnt < seq_pairs.size(); ++cnt) {

    // get sequence indices for this sequence pair
    const size_t i = seq_pairs[cnt].first;
    const size_t j = seq_pairs[cnt].second;

    // map the original sequence positions to the new reduced ones
    for (Post_probs::iterator post = post_probs_list[cnt].begin(); post != post_probs_list[cnt].end(); ++post) {
      (*post).x = reduced_position_reverse_map[i][(*post).x];
      (*post).y = reduced_position_reverse_map[j][(*post).y];
    }

    // if reduced seqs are length 0, then leave null
    if ((anchored_positions[i].size() == 0) || (anchored_positions[j].size() == 0))
      continue;

    // else create corresponding sparse matrix and store transpose
    sparse_matrices[i][j] = new SparseMatrix (i, j,
					      anchored_positions[i].size(), anchored_positions[j].size(), post_probs_list[cnt]);
    sparse_matrices[j][i] = sparse_matrices[i][j]->ComputeTranspose();

  }

  // do database stuff
  // create empty database manager 
  Manager manager;

  // now resolve anchors
  CTAG (9,ANCHORING ANCHORING_VERBOSE) << "Performing anchor annealing to resolve anchors." << endl;

  // do anchor annealing
  // It's important to perform use gap-factor 0;
  // otherwise the result will depend heavily on how we scale the scores for the HSPs.
  // This way will give the maximal set (well, a greedy version thereof) of consistent anchors.
  Alignment_DAG dag (reduced_seq_db);
  dag.anneal (sparse_matrices, tree_weights,
	      manager,
	      false,                    // use_tgf = false
	      0, true, 0,               // gap_factor = 0, enable_dynamic_weights = true, edge_weight_threshold = 0
	      num_refinement_steps,
	      false, "");               // output_for_gui = false, gui_prefix = ""
  dag.dfs_topological_sort();  // this isn't really necessary but it does force removal of dead columns, which is nice

  // convert the annealed alignment into a set of pairwise anchors
  // initialize them with sequence lengths for score calculation
  std::vector<Anchors> resolved_anchors_list;
  for (size_t cnt = 0; cnt < seq_pairs.size(); ++cnt)
    resolved_anchors_list.push_back (Anchors (seq_db_internal.get_seq (seq_pairs[cnt].first).seq, seq_db_internal.get_seq (seq_pairs[cnt].second).seq));

  // we need a mapping from sequence pairs to the corresponding index in seq_pairs
  // -1 means not present; otherwise is the index 
  std::map<std::pair<size_t, size_t>, size_t> seq_pair_lookup;
  for (size_t cnt = 0; cnt < seq_pairs.size(); ++cnt)
    seq_pair_lookup.insert (std::make_pair (std::make_pair (seq_pairs[cnt].first, seq_pairs[cnt].second), cnt));

  // loop over columns of the alignment,
  // interpreting aligned characters as constituting a resolved anchor

  // don't store an anchor more than once:
  // candidate_anchors_added[cnt][a] = true indicates that candidate_anchors_list[cnt][a]
  // has already been added as a resolved anchor
  std::vector<std::vector<bool> > candidate_anchors_added (candidate_anchors_list.size(), std::vector<bool> ());
  for (size_t cnt = 0; cnt < candidate_anchors_list.size(); ++cnt)
    candidate_anchors_added[cnt].assign (candidate_anchors_list[cnt].size(), false);

  // now loop over columns of the DAG to look for aligned positions
  const std::vector<Column*>& columns = dag.get_columns();
  for (std::vector<Column*>::const_iterator col = columns.begin(); col != columns.end(); ++col) {
    
    // there shouldn't be any dead columns here
    assert (!(*col)->is_dead());

    const Seq_pos_map seq_pos_map = (*col)->get_seq_pos_map();

    // if no aligned positions, continue on to the next column
    if (seq_pos_map.size() < 2)
      continue;

    // create a pairwise anchor for each pair of aligned characters
    for (Seq_pos_map::const_iterator seq_pos1 = seq_pos_map.begin(); seq_pos1 != seq_pos_map.end(); ++seq_pos1) {
      for (Seq_pos_map::const_iterator seq_pos2 = seq_pos_map.begin(); seq_pos2 != seq_pos_map.end(); ++seq_pos2) {

	// sequence indices
	const size_t i = seq_pos1->first;
	const size_t j = seq_pos2->first;
	assert ((i < seq_db.size()) && (j < seq_db.size()));

	// require i < j to prevent duplicates
	if (i >= j)
	  continue;

	// make sure that the sequence pair (i, j) is in the list seq_pairs of sequence pairs to process;
	// if not, then continue on to the next pair
	if (seq_pair_lookup.find (std::make_pair (i, j)) == seq_pair_lookup.end())
	  continue;
	const size_t cnt = seq_pair_lookup[std::make_pair (i, j)];

	const unsigned x = (reduced_position_map[i])[(*seq_pos1).second]; // map from positions within the reduced sequences to the original sequences
	const unsigned y = (reduced_position_map[j])[(*seq_pos2).second]; //     (*seq_pos1).second   ->     reduced_position_map[j][(*seq_pos1).second]
	assert ((x < seq_db.get_seq (i).length()) && (y < seq_db.get_seq (j).length()));

	// here "transitive anchoring" comes into play:
	// If it's a "transitive anchor" then we won't have information stored for it
	// (ie no entry in candidate_anchors[cnt]), but we don't want to discard the transitive information either.
	// We therefore should do a lookup to see if a candidate anchor (x,y) exists for seq. pair cnt;
	// if not, then just create an anchor of length one with midpoints (x,y) and store it.
	const size_t candidate_index = candidate_anchors_list[cnt].exists_spanning_anchor (x, y);
	if (candidate_index != candidate_anchors_list[cnt].size()) {                                // if there is a spanning candidate anchor
	  if (!candidate_anchors_added[cnt][candidate_index]) {    // if it hasn't already been added, then add it
	    resolved_anchors_list[cnt].store (candidate_anchors_list[cnt].get_spanning_anchor (x, y));
	    candidate_anchors_added[cnt][candidate_index] = true;
	  } else                                                   // else we've already added it, so just continue
	    continue;
	} else {                                                   // if there is no spanning candidate anchor
	  resolved_anchors_list[cnt].store (Anchor (Interval (x, x), Interval (y, y)));
	  resolved_anchors_list[cnt].update_score (resolved_anchors_list[cnt].size() - 1);
	}

      }
    }

  }

  // log resolved anchors
  if (CTAGGING(4,ANCHORING ANCHORING_VERBOSE)) {
    for (size_t cnt = 0; cnt < resolved_anchors_list.size(); ++cnt) {
      CL << "Resolved anchors after anchor annealing for '" << seq_db_internal.get_seq (seq_pairs[cnt].first).name << "' and '" << seq_db_internal.get_seq (seq_pairs[cnt].second).name << "':" << endl;
      resolved_anchors_list[cnt].show (CL);
    }
  }

  // perform various administrative tasks:
  // - resolve adjacent parallel anchors
  // - enforce immutable anchors unique
  // - prune overlapping anchors
  // because we have already resolved candidate parallel anchors before anchor annealing,
  // this should be relevant only for catching adjacent parallel anchors which arise due to transitive anchoring
  // also sort lexically
  for (size_t cnt = 0; cnt < resolved_anchors_list.size(); ++cnt) {
    resolved_anchors_list[cnt].enforce_immutable_unique();
    if (!resolved_anchors_list[cnt].resolve_parallel (max_join_distance, true)) // concatenate immutable anchors
      THROWEXPR ("ERROR: Couldn't resolve parallel anchors.");
    resolved_anchors_list[cnt].remove_overlaps (use_translated ? 3 * minlen : minlen); // remove pruned anchors which are shorter than minlen
  }

  // debugging
  if (CTAGGING(-1,ANCHORING ANCHORING_VERBOSE_VERBOSE)) {
    dag.get_stockholm().write_stockholm (CL);
  }

  // log resolved anchors
  if (CTAGGING(4,ANCHORING ANCHORING_VERBOSE)) {
    for (size_t cnt = 0; cnt < resolved_anchors_list.size(); ++cnt) {
      CL << "Resolved anchors for '" << seq_db_internal.get_seq (seq_pairs[cnt].first).name << "' and '" << seq_db_internal.get_seq (seq_pairs[cnt].second).name << "':" << endl;
      resolved_anchors_list[cnt].show (CL);
    }
  }

  // write GUI output if requested
  if (output_for_gui) {
    std::string gui_filename = std::string (gui_prefix) + GUI_FILENAME_ANCHORS_SUFFIX;
    std::ofstream gui_file;
    gui_file.open (gui_filename.c_str());
    if (!gui_file.is_open())
      THROWEXPR ("ERROR: Couldn't create file with name '" << gui_filename << "'.");
    for (size_t cnt = 0; cnt < resolved_anchors_list.size(); ++cnt) {
      gui_file << "; Resolved anchors after concatenation for '" << seq_db_internal.get_seq (seq_pairs[cnt].first).name << "' and '" << seq_db_internal.get_seq (seq_pairs[cnt].second).name << "'." << endl
	       << "; Format is " << endl
	       << ";   (0 ~ 1) == [1,6] ~ [3,8] => 0.9" << endl
	       << ";   meaning that sequences 0 and 1 have an ungapped anchor " << endl
	       << ";   which spans [1,6] in sequence 0 and [3,8] in sequence 1" << endl
	       << ";   and has a score of 0.9." << endl
	       << "; sequence is 0-based and position is 0-based" << endl
	       << endl;
      resolved_anchors_list[cnt].write_gui_output (gui_file, seq_pairs[cnt].first, seq_pairs[cnt].second);
      gui_file << endl;
    }
    gui_file.close();
    CTAG (9,ANCHORING ANCHORING_VERBOSE) << "Created GUI output file '" << gui_filename << "'" << endl;
  }

  // impose normalized probability distribution
  for (size_t cnt = 0; cnt < resolved_anchors_list.size(); ++cnt)
    resolved_anchors_list[cnt].impose_normalized_probability_distribution();

  return resolved_anchors_list;
}
