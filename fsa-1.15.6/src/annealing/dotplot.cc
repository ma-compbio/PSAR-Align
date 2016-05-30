
/**
 * \file dotplot.cc
 * This file is part of FSA, a sequence alignment algorithm.
 * \author Source code in this file was written by Ian Holmes, Lars Barquist and Robert Bradley.
 */

#include "dotplot.h"

using namespace fsa;

Dotplot::Dotplot (const Sequence& x, const Sequence& y)
  : array2d<double> (x.length(), y.length(), 0.), xseq (x.name), yseq (y.name)
{ }

void Dotplot::write_dotplot (const std::string& filename, const double cutoff /* = PROB_DISPLAY_CUTOFF */) const
{
  std::ofstream file (filename.c_str());
  if (!file)
    THROWEXPR ("ERROR: Couldn't create dotplot file with name '" << filename << "'.");
  // print horizontal sequence axis labels
  file << ".";
  for (unsigned x = 0; x < xseq.size(); ++x)
    file << ' ' << xseq[x];
  file << '\n';
  // print rows
  for (unsigned y = 0; y < yseq.size(); ++y) {
    file << yseq[y];
    for (unsigned x = 0; x < xseq.size(); ++x) {
      const double& p = (*this) (x, y);
      file << ' ' << (p < cutoff ? '0' : p);
    }
    file << '\n';
  }
	
  file.close();
}


void Dotplot::write_dotplot (const std::string& prefix, const std::string& seqname, const double cutoff /* = PROB_DISPLAY_CUTOFF */) const
{
  std::string filename;
  filename += prefix + '-' + seqname;
  write_dotplot (filename, cutoff);
}

void Dotplot::write_dotplot (const std::string& prefix, const std::string& xseqname, const std::string& yseqname, const double cutoff /* = PROB_DISPLAY_CUTOFF */) const
{
  std::string filename;
  filename += prefix + '-' + xseqname + '-' + yseqname;
  write_dotplot (filename, cutoff);
}

