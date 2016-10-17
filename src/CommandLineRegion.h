#ifndef VBAM_COMMAND_LINE_REGION_H__
#define VBAM_COMMAND_LINE_REGION_H__

#include "SeqLib/ReadFilter.h"

struct CommandLineRegion {
  
CommandLineRegion(const std::string& mf, int t) : f(mf), type(t), pad(0), i_flag(0), e_flag(0) {}

  std::string f; // file
  int type; // mate linked, excluder, etc
  int pad;
  uint32_t i_flag; // inclusive flags
  uint32_t e_flag; // exclusive flags

  uint32_t any_i_flag; // inclusive any flags (eg if read has any bit of these, keep)
  uint32_t any_e_flag; // exclusive any flags (eg if read has any bit of these, fail)

  int len = 0;
  int mapq = 0;
  int nbases = INT_MAX;
  int phred = 0;
  int clip = 0;
  int ins = 0;
  int del = 0;
  std::string rg;
  std::string motif;

  bool all() const { 
    return !len && !mapq && !nbases && !phred && rg.empty() && !i_flag && !e_flag && !any_i_flag && !any_e_flag; 
  }

};

static SeqLib::Filter::ReadFilter BuildReadFilterFromCommandLineRegion(const CommandLineRegion& c, const SeqLib::BamHeader& hdr) {

  SeqLib::Filter::ReadFilter r;
  //r.m_region_file = c.f;
  
  // set a whole genome ALL rule
  if (c.type < 0) {
    //m_grv.clear();
    //id = "WG";
  } else {
    // set the genomic region this rule applies to
    SeqLib::GRC regr(c.f, hdr);
    regr.Pad(c.pad);
    r.setRegions(regr);
    //debug setRegionFromFile(c.f, hdr);
  }
  
  // add the abstract rule
  SeqLib::Filter::AbstractRule ar;
  
  // set the flag
  if (c.i_flag || c.e_flag) {
    ar.fr.setAllOnFlag(c.i_flag);
    ar.fr.setAllOffFlag(c.e_flag);
  }
  
  // set the other fields
  if (c.len)
    ar.len = SeqLib::Filter::Range(c.len, INT_MAX, false);
  if (c.nbases != INT_MAX)
    ar.nbases = SeqLib::Filter::Range(0, c.nbases, false);
  if (c.phred)
    ar.phred = SeqLib::Filter::Range(c.phred, INT_MAX, false);
  if (c.mapq)
    ar.mapq = SeqLib::Filter::Range(c.mapq, INT_MAX, false);
  if (c.clip)
    ar.clip = SeqLib::Filter::Range(c.clip, INT_MAX, false);
  if (c.del)
    ar.del = SeqLib::Filter::Range(c.del, INT_MAX, false);
  if (c.ins)
    ar.ins = SeqLib::Filter::Range(c.ins, INT_MAX, false);
  
  // set the id
  //ar.id = id + "_CMD_RULE";
  
  // add a motif rule
  if (!c.motif.empty())
    ar.addMotifRule(c.motif, false);
  
  // add read group rule
  ar.SetReadGroup(c.rg);
  
  r.AddRule(ar);
  
  // set the properties of the region
  if (c.type >= 0) {
    switch(c.type) {
    case MINIRULES_MATE_LINKED:
      r.SetMateLinked(true);
      r.SetExcluder(false);
      //m_applies_to_mate = true;
      //excluder = false;
      break;
    case MINIRULES_MATE_LINKED_EXCLUDE:
      r.SetMateLinked(true);
      r.SetExcluder(true);
      //m_applies_to_mate = true;
      //excluder = true;
      break;
    case MINIRULES_REGION:
      r.SetMateLinked(false);
      r.SetExcluder(false);
      //m_applies_to_mate = false;
      //excluder = false;
      break;
    case MINIRULES_REGION_EXCLUDE:
      r.SetMateLinked(false);
      r.SetExcluder(true);
      //m_applies_to_mate = false;
      //excluder = true;
      break;
    default:
      std::cerr << "Unexpected type in ReadFilter. Exiting" << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  return r;
}

#endif
