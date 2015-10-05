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
  int32_t buffer_size = 500;
  SnowTools::BamReadVector buffer;

  while (GetNextRead(r, rule))
    {

      // prepare for case of long reads
      buffer_size = std::max((int32_t)r.Length() * 5, buffer_size);

      //std::cerr << "...r " << r << " rule " << rule << std::endl;
      TrackSeenRead(r);

      // add coverage
      if (max_cov > 0) {
	cov_a.addRead(r);
	cov_b.addRead(r);
      }

      // read is valid
      if (rule) {

	if (max_cov <= 0 && fop) {// if we specified an output file, write it
	  WriteAlignment(r);
	  ++rc_main.keep;
	} else if (fop) {
	  buffer.push_back(r);

	  // clear buffer
	  // pass back and forth between cov_a and cov_b.
	  if (buffer.size()) {
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
      double this_cov = std::max(this_cov1, this_cov2);
      double sample_rate = 1; // dummy, always set if max_coverage > 0
      if (this_cov > 0) 
	sample_rate = 1 - (this_cov - max_cov) / this_cov; // if cov->inf, sample_rate -> 0. if cov -> max_cov, sample_rate -> 1
      
      // this read should be randomly sampled, cov is too high
      if (this_cov > max_cov) 
	{
	  uint32_t k = __ac_Wang_hash(__ac_X31_hash_string(r.Qname().c_str()) ^ m_seed);
	  if ((double)(k&0xffffff) / 0x1000000 <= sample_rate) { // passed the random filter
	    WriteAlignment(r);
	    ++rc_main.keep;
	  }
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
