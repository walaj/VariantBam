#ifndef VARIANT_VARIANT_BAM_WALKER_H__
#define VARIANT_VARIANT_BAM_WALKER_H__

#include "SnowTools/BamWalker.h"
#include "SnowTools/BamStats.h"
#include "SnowTools/BamRead.h"
#include "SnowTools/STCoverage.h"

class VariantBamWalker: public SnowTools::BamWalker
{
 public:

  VariantBamWalker() {}

  VariantBamWalker(const std::string in) : SnowTools::BamWalker(in) {}

  void writeVariantBam();
  
  void TrackSeenRead(SnowTools::BamRead &r);
  
  void printMessage(const SnowTools::BamRead &r) const;
  
  SnowTools::BamStats m_stats;

  SnowTools::STCoverage cov_a;
  SnowTools::STCoverage cov_b;

  int m_seed = 0;
  
  int max_cov = 0;

  void subSampleWrite(SnowTools::BamReadVector& buff, const SnowTools::STCoverage& cov);

  SnowTools::ReadCount rc_main;

};
#endif
