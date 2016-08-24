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

  VariantBamWalker(const std::string in) : SeqLib::BamReader(in) {}

  void writeVariantBam();
  
  void TrackSeenRead(SeqLib::BamRecord &r);
  
  void printMessage(const SeqLib::BamRecord &r) const;
  
  BamStats m_stats;

  STCoverage cov_a;
  STCoverage cov_b;

  int m_seed = 0;
  
  int max_cov = 0;

  void subSampleWrite(SeqLib::BamRecordVector& buff, const STCoverage& cov);

  ReadCount rc_main;

  bool m_verbose = false;
  
  SeqLib::ReadFilterCollection m_mr;

  SeqLib::BamWriter m_writer;

  bool m_strip_all_tags = false;

  std::vector<std::string> m_tags_to_strip;

};
#endif
