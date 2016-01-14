#include "VariantBamWalker.h"
#include "htslib/khash.h"

void VariantBamWalker::writeVariantBam() 
{

#ifndef __APPLE__
  // start the timer
  clock_gettime(CLOCK_MONOTONIC, &start);
#endif

  SnowTools::BamRead r;
  bool rule;

  bool COV_A = true;
  int32_t buffer_size = 10000;
  SnowTools::BamReadVector buffer;

  // check if the BAM is sorted by looking at the header
  std::string hh = std::string(header()->text);
  bool sorted = hh.find("SO:coord") != std::string::npos;

  if (!sorted && max_cov > 0) {
    std::cerr << "ERROR: BAM file does not appear to be sorted (no SO:coordinate) found in header." << std::endl;
    std::cerr << "       Sorted BAMs are required for coverage-based rules (max/min coverage)." << std::endl;
    exit(EXIT_FAILURE);
  }

  while (GetNextRead(r, rule)) {

      // prepare for case of long reads
      buffer_size = std::max((int32_t)r.Length() * 5, buffer_size);

      //std::cerr << "...r " << r << " rule " << rule << std::endl;
      TrackSeenRead(r);

      // add coverage
      if (max_cov != 0) {
	cov_a.addRead(r);
	cov_b.addRead(r);
      }

      // read is valid
      if (rule) {

	if (max_cov == 0 && fop) {// if we specified an output file, write it
	  WriteAlignment(r);
	  ++rc_main.keep;
	} else if (fop) {
	  buffer.push_back(r);

	  // clear buffer
	  // pass back and forth between cov_a and cov_b.
	  if (buffer.size()) {

	    // error if BAM not sorted
	    if (buffer[0].Position() - buffer.back().Position() > 0 && buffer[0].ChrID() == buffer.back().ChrID())
	      {
		std::cerr << "ERROR: BAM file is not sorted. " << std::endl;
		std::cerr << " ------ Found read:  " << buffer[0] << std::endl;
		std::cerr << " ------ before read: " << buffer.back() << std::endl;
		std::cerr << " ------ BAM must be sorted if using the -m flag for max coverage. Exiting" << std::endl;
		exit(EXIT_FAILURE);
	      }

	    if ( (buffer.back().Position() - buffer[0].Position() > buffer_size) || buffer.back().ChrID() != buffer[0].ChrID()) {
	      COV_A ? subSampleWrite(buffer, cov_a) : subSampleWrite(buffer, cov_b);
	      COV_A ? cov_a.clear() : cov_b.clear();
	      COV_A = !COV_A;
	      buffer.clear();
	    }
	  }
	}
      }
      
      if (++rc_main.total % 1000000 == 0 && m_verbose)
	printMessage(r);

    }

  // clear last buffer
  if (buffer.size()) {
    if ( (buffer.back().Position() - buffer[0].Position() > buffer_size) || buffer.back().ChrID() != buffer[0].ChrID()) {
      COV_A ? subSampleWrite(buffer, cov_a) : subSampleWrite(buffer, cov_b);
      COV_A ? cov_a.clear() : cov_b.clear();
      COV_A = !COV_A;
      buffer.clear();
    }
  }

  if (m_verbose)
    printMessage(r);
}

void VariantBamWalker::subSampleWrite(SnowTools::BamReadVector& buff, const SnowTools::STCoverage& cov) {

  for (auto& r : buff)
    {
      double this_cov1 = cov.getCoverageAtPosition(r.ChrID(), r.Position());
      double this_cov2 = cov.getCoverageAtPosition(r.ChrID(), r.PositionEnd());
      //double this_cov3 = cov.getCoverageAtPosition(r.ChrID(), r.Position());
      //double this_cov4 = cov.getCoverageAtPosition(r.ChrID(), r.PositionEnd());
      double this_cov = std::max(this_cov1, this_cov2);
      double sample_rate = 1; // dummy, always set if max_coverage > 0
      if (this_cov > 0) 
	sample_rate = 1 - (this_cov - max_cov) / this_cov; // if cov->inf, sample_rate -> 0. if cov -> max_cov, sample_rate -> 1
      
      // this read should be randomly sampled, cov is too high
      if (this_cov > max_cov && max_cov > 0) 
	{
	  uint32_t k = __ac_Wang_hash(__ac_X31_hash_string(r.Qname().c_str()) ^ m_seed);
	  if ((double)(k&0xffffff) / 0x1000000 <= sample_rate) { // passed the random filter
	    WriteAlignment(r);
	    ++rc_main.keep;
	  }
	}
      // only take if reaches minimum coverage
      else if (this_cov < -max_cov) // max_cov = -10 
	{
	  std::cerr << "not writing because this cov is " << this_cov << " and min cov is " << (-max_cov) << std::endl;
	}
      else // didn't have a coverage problems
	{
	  ++rc_main.keep;
	  WriteAlignment(r);
	}
      
    }


}

void VariantBamWalker::TrackSeenRead(SnowTools::BamRead &r)
{
  m_stats.addRead(r);
}

void VariantBamWalker::printMessage(const SnowTools::BamRead &r) const 
{

  if (rc_main.total <= 0) {
    std::cerr << "NO READS FOUND" << std::endl;
    return;
  }

  char buffer[120];
  std::string posstring = SnowTools::AddCommas<int>(r.Position());
  std::sprintf (buffer, "Reading read %11s at %2s:%-11s. Kept %10s (%2d%%) [running count]",  
		rc_main.totalString().c_str(), SnowTools::GenomicRegion::chrToString(r.ChrID()).c_str(), posstring.c_str(),  
		rc_main.keepString().c_str(), rc_main.percent());
  //buffer[99] = '\0';
  std::string outr(buffer);
  //std::printf ("%s | ",buffer);
  std::cerr << outr << SnowTools::displayRuntime(start) << std::endl;;
  
}
