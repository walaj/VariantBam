#ifndef SNOWTOOLS_BAMSTATS_H__
#define SNOWTOOLS_BAMSTATS_H__

#include <unordered_map>
#include <cstdint>
#include <iostream>

#include "Histogram.h"
#include "SeqLib/BamRecord.h"

/** Small class to store a counter to measure BamWalker progress.
 *
 * Currently only stores number of reads seen / kept. 
 */
struct ReadCount {

  uint64_t keep = 0;
  uint64_t total = 0;
  
  /** Return the percent of total reads kept
   */
  int percent () const {
    int perc  = SeqLib::percentCalc<uint64_t>(keep, total); 
    return perc;
  }

  /** Return the total reads visited as a comma-formatted string
   */
  std::string totalString() const {
    return SeqLib::AddCommas<uint64_t>(total);
  }

  /** Return the kept reads as a comma-formatted string
   */
  std::string keepString() const {
    return SeqLib::AddCommas<uint64_t>(keep);
  }

};

  /** Store information pertaining to a given read group *
   *
   * This class will collect statistics on number of: read, supplementary reads, unmapped reads, qcfail reads, duplicate reads.
   * It will also create Histogram objects to store counts of: mapq, nm, isize, clip, mean phred score, length
   */
class BamReadGroup {

  friend class BamStats;
  
 public:

  /** Construct an empty BamReadGroup */
  BamReadGroup() {}

  /** Construct an empty BamReadGroup for the specified read group
   * @param name Name of the read group
   */
  BamReadGroup(const std::string& name);

  /** Display basic information about this read group
   */
  friend std::ostream& operator<<(std::ostream& out, const BamReadGroup& rg);

  /** Add a BamRecord to this read group */
  void addRead(SeqLib::BamRecord &r);

 private:

  size_t reads;
  size_t supp;
  size_t unmap;  
  size_t qcfail;
  size_t duplicate;
  size_t mate_unmap;

  Histogram mapq;
  Histogram nm;
  Histogram isize;
  Histogram clip;
  Histogram phred;
  Histogram len;

  std::string m_name;

};

/** Class to store statistics on a BAM file.
 *
 * BamStats currently stores a map of BamReadGroup objects. Bam statistics
 * are collected then on a read-group basis, but can be output in aggregate. See
 * BamReadGroup for description of relevant BAM statistics.
 */
class BamStats
{

 public:
  
  /** Loop through the BamReadGroup objections and print them */
  friend std::ostream& operator<<(std::ostream& out, const BamStats& qc);

  /** Add a read by finding which read group it belongs to and calling the 
   * addRead function for that BamReadGroup.
   */
  void addRead(SeqLib::BamRecord &r);

  std::unordered_map<std::string, BamReadGroup> m_group_map;

};

#endif
