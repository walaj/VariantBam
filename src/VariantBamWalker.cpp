#include "VariantBamWalker.h"
#include "htslib/khash.h"

void VariantBamWalker::writeVariantBam() {

#ifndef __APPLE__
  // start the timer
  //clock_gettime(CLOCK_MONOTONIC, &start);
#endif

  SeqLib::BamRecord r;

  bool COV_A = true;
  int32_t buffer_size = 10000;
  SeqLib::BamRecordVector buffer;

  // check if the BAM is sorted by looking at the header
  std::string hh = Header().AsString(); //std::string(header()->text);
  bool sorted = hh.find("SO:coord") != std::string::npos;

  if (!sorted && max_cov > 0) {
    std::cerr << "ERROR: BAM file does not appear to be sorted (no SO:coordinate) found in header." << std::endl;
    std::cerr << "       Sorted BAMs are required for coverage-based rules (max/min coverage)." << std::endl;
    exit(EXIT_FAILURE);
  }

  // check that regions are sufficient size
  for (auto& k : m_region)
    if (k.Width() < 1000)
      k.Pad(1000);

  while (GetNextRecord(r)) {

    int s, e;
    if (phred  > 0) {
      std::string seq = r.Sequence();
      r.QualityTrimmedSequence(phred, s, e);
      int new_len = e - s;
      if (e != -1 && new_len < r.Length() && new_len > 0 && new_len - s >= 0 && s + new_len <= r.Length())
	r.AddZTag("GV", seq.substr(s, new_len));
    }

    bool rule = m_mr.isValid(r);
    
    // prepare for case of long reads
    buffer_size = std::max((int32_t)r.Length() * 5, buffer_size);
    
    TrackSeenRead(r);
    
    // add coverage
    if (max_cov != 0) {
      cov_a.addRead(r, 0, false);
      cov_b.addRead(r, 0, false);
    }
    
    // read is valid
    if (rule) {
      
      if (max_cov == 0 && m_writer.IsOpen()) { // if we specified an output file, write it
	write_record(r);
      } else if (m_writer.IsOpen()) {
	buffer.push_back(r);
	
	// clear buffer
	// pass back and forth between cov_a and cov_b.
	if (buffer.size()) {
	  
	  // error if BAM not sorted
	  if (buffer[0].Position() - buffer.back().Position() > 0 && buffer[0].ChrID() == buffer.back().ChrID()) {
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
      } else if (!m_writer.IsOpen()) { // we are not outputting anything
	++rc_main.keep;
      }
      
    } else if (m_mark_qc_fail) { // fails, but we should mark it and write
      r.SetQCFail(true);
      write_record(r);
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

  if (r.isEmpty()) {
    std::cerr << "NO READS RETRIEVED FROM THESE REGIONS" << std::endl;
    return;
  }

  if (m_verbose)
    printMessage(r);

}

void VariantBamWalker::subSampleWrite(SeqLib::BamRecordVector& buff, const STCoverage& cov) {

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
	    write_record(r);
	  } else if (m_mark_qc_fail) {
	    r.SetQCFail(true);
	    write_record(r);
	  }
	}
      // only take if reaches minimum coverage
      else if (this_cov < -max_cov) { // max_cov = -10 
	if (m_mark_qc_fail) {
	  r.SetQCFail(true);
          write_record(r);
	} else {
	  //std::cerr << "not writing because this cov is " << this_cov << " and min cov is " << (-max_cov) << std::endl;
	}
      } else {
        write_record(r);
      }
      
    }
  
  
}

void VariantBamWalker::TrackSeenRead(SeqLib::BamRecord &r)
{
  m_stats.addRead(r);
}

void VariantBamWalker::printMessage(const SeqLib::BamRecord &r) const 
{

  char buffer[90];

  std::string posstring = SeqLib::AddCommas<int>(r.Position());
  std::string chrname;
  if (r.ChrID() == -1)
    chrname = "Unmapped";
  else {
    try { 
      chrname = Header().IDtoName(r.ChrID()); 
    } catch(...) { 
      chrname = "CHR_NAME_FAIL"; 
    } 
  }

  if (rc_main.total <= 0) {
    std::sprintf(buffer, "NO READS FOUND at %2s:%-11s",chrname.c_str(), posstring.c_str());
    std::cerr << std::string(buffer) << std::endl;
    return;
  }

  std::sprintf (buffer, "Read %11s at %2s:%-11s. Kept %10s (%2d%%) -- ",  
		rc_main.totalString().c_str(), chrname.c_str(), posstring.c_str(),  
		rc_main.keepString().c_str(), rc_main.percent());
  std::cerr << std::string(buffer) << std::endl;
}

void VariantBamWalker::write_record(SeqLib::BamRecord& r) {

  if (m_write_trimmed) {
    r.SetSequence(r.QualitySequence());
    r.SetQualities(std::string(), 0);
  }

  // strip tags
  if (m_strip_all_tags)
    r.RemoveAllTags();
  else if (m_tags_to_strip.size()) {
    for (const auto& i : m_tags_to_strip)
      r.RemoveTag(i.c_str());
  }

  // write it
  m_writer.WriteRecord(r);

  ++rc_main.keep;
}
