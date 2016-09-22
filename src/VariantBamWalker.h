#ifndef VARIANT_VARIANT_BAM_WALKER_H__
#define VARIANT_VARIANT_BAM_WALKER_H__

#include "SeqLib/BamReader.h"
#include "SeqLib/BamWriter.h"
#include "SeqLib/ReadFilter.h"
#include "BamStats.h"
//#include "SnowTools/BamRead.h"
#include "STCoverage.h"

class VariantBamWalker: public SeqLib::BamReader
{
 public:

   VariantBamWalker() {}

  void writeVariantBam();
  
  void TrackSeenRead(SeqLib::BamRecord &r);
  
  void printMessage(const SeqLib::BamRecord &r) const;
  
  BamStats m_stats;

  bool m_write_trimmed = false; // output the phred trimmed instead of orig sequence

  STCoverage cov_a;
  STCoverage cov_b;

  int m_seed = 0;
  
  int max_cov = 0;

  void subSampleWrite(SeqLib::BamRecordVector& buff, const STCoverage& cov);

  int phred = -1;

  ReadCount rc_main;

  bool m_verbose = false;
  
  SeqLib::Filter::ReadFilterCollection m_mr;

  SeqLib::BamWriter m_writer;

  bool m_strip_all_tags = false;

  std::vector<std::string> m_tags_to_strip;

 private:

  void write_record(SeqLib::BamRecord& r);

};
#endif
