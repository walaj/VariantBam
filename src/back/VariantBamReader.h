#ifndef VARIANT_BAM_READER_H
#define VARIANT_BAM_READER_H

#include "MiniRules.h"
#include "GenomicRegion.h"
#include "BamQC.h"
#include <time.h>
#include "reads.h"
#include "SnowUtils.h"
#include "assert.h"

using namespace std;

// Phred score transformations
inline int char2phred(char b) {
  uint8_t v = b;
  assert(v >= 33);
  return v - 33;
}

/////////////// 
// Hold read counts
//////////////
struct ReadCount {

  int keep = 0;
  int total = 0;
  
  int percent () const {
    int perc  = SnowUtils::percentCalc<int>(keep, total); 
    return perc;
  }

  string totalString() const {
    return SnowUtils::AddCommas<int>(total);
  }

  string keepString() const {
    return SnowUtils::AddCommas<int>(keep);
  }

};


class VariantBamReader {

 public:
  VariantBamReader() { }
  ~VariantBamReader() {
#ifdef HAVE_BAMTOOLS
    m_writer->Close();
    m_reader->Close();
    delete m_writer;
    delete m_reader;
#endif

#ifdef HAVE_HTSLIB
    bgzf_close(fp);
    bam_hdr_destroy(br);
    hts_itr_destroy(hts_itr);
    hts_idx_destroy(idx);
    sam_close(fop);
#endif
  }

  void saveAlignment(Read &r);

  unordered_map<string, bool> m_hash;

  bool twopass = false;

  struct timespec start;
  
  VariantBamReader(string inbam, string outbam, MiniRulesCollection* mr, int verbose);

  void dumpBuffer(ReadVec &buff, ReadVec &store, int mapq);

  //static unsigned getClipCount(BamAlignment a);
  static void qualityTrimRead(int qualTrim, string &seq, string &qual);

  static int32_t qualityTrimRead(int qualTrim, int32_t &startpoint, Read &r);

  //bool writeVariantBam(BamQC &qc, bool qc_only);
  //bool writeVariantBam(BamQC &qc, BamAlignmentVector &bav);
  //bool writeVariantBam(BamQC &qc, bam1_v &bav);
  bool writeVariantBam(BamQC &qc, ReadVec &bav);

  void printRuleCounts(unordered_map<string, size_t> &rm) const;
  
  // set which part of the bam to read
  bool setBamRegion(GenomicRegion gp);

  // create the index file for the output bam
  void MakeIndex();

  // print to stdout
  void printMessage(const ReadCount &rc_main, const Read &r) const;

  void writeVariantBamFromHash();

  void rewind() { 
    #ifdef HAVE_BAMTOOLS
    m_reader->Rewind(); 
    #endif
  }

 private:
  
  string m_bam;
  string m_out;
#ifdef HAVE_BAMTOOLS
  BamTools::BamReader * m_reader;
  BamTools::BamWriter * m_writer;
#endif
  GenomicRegion m_region;
  MiniRulesCollection * m_mr;
  int m_verbose;

#ifdef HAVE_HTSLIB
  // hts
  BGZF * fp = 0;
  hts_idx_t * idx = 0; // hts_idx_load(bamfile.c_str(), HTS_FMT_BAI);
  hts_itr_t * hts_itr = 0; // sam_itr_queryi(idx, 3, 60000, 80000);
  bam_hdr_t * br = 0;

  samFile* fop = 0;
  //fp = sam_open(fn, mode);
#endif

};


#endif 


