#ifndef MINI_RULES_H
#define MINI_RULES_H

/* Define a set of rules for creating a variant bam. The syntax is:
   all@!isize:[0,800],mapq:[0,60]
   region@REGION_FILE
   rule1@isize:[0,800],mapq:[0,60]
   rule2@!isize[0,800]:mapq[0,60],:ardclip:supplementary:duplicate:qcfail
   rule3@
   
   A file of NA indicates that the rule should be applied genome-wide.
   The ordering of the lines sets the hierarchical rule. For instance, a rule on line 2 will be applied 
   before a rule on line 3 for all regions that are the union of regions in level 3 and below.
   
   e.g. Level 3 region file has region chr1   100   1000
        Level 2 region file has region chr1   150   1200
	The union of these will produce a new region chr1   100   1200, with level 2
*/

#include <string>
#include <vector>
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "GenomicRegion.h"

using namespace BamTools;
using namespace std;

// hold a range of valid numeric values (e.g. isize). 
// can optionally invert the range to make rule the complement of the range
struct Range {

  Range() : min(0), max(0), inverted(false), pattern("") {}
  Range(int mn, int mx, int in, string p) : min(mn), max(mx), inverted(in), pattern(p) {}
  ~Range() {}

  int min;
  int max;
  bool inverted;
  string pattern;
  
  bool isValid(int val) {
    if (!inverted)
      return (val >= min && val <= max);
    else
      return (val < min || val > max);
  }

  void parseRuleLine(string line);

  friend ostream& operator<<(ostream &out, const Range &r);
};

// a container to hold boolean rules based mostly on alignment flag
struct FlagRule {
  
  FlagRule() : duplicate(true), supp(true), qcfail(true), hardclip(true), 
               fwd_strand(true), rev_strand(true), mate_fwd_strand(true), mate_rev_strand(true),
	       unmapped(true), unmapped_mate(true), mapped(true), mapped_mate(true)
               {}

  void parseRuleLine(string line);
  
  bool isValid(BamAlignment &a);

  bool duplicate;
  bool supp;
  bool qcfail;
  bool hardclip;
  bool fwd_strand;
  bool rev_strand;
  bool mate_fwd_strand;
  bool mate_rev_strand;
  bool unmapped;
  bool unmapped_mate;
  bool mapped;
  bool mapped_mate;

  friend ostream& operator<<(ostream &out, const FlagRule &fr);
};

//
class AbstractRule {

 public:

  string name = "";
  Range isize = {-1, -1, true, "isize"}; // include all
  Range mapq =  {-1, -1, true, "mapq"}; 
  Range len =   {-1, -1, true, "length"};
  Range clip =  {-1, -1, true, "clip"};
  Range phred = {-1, -1, true, "phred"};

  bool keep_all = false;
  bool keep_none = false;

  // set to true if you want a read to belong to the region if its mate does
  //bool mate = false; 

  FlagRule fr;

  bool isValid(BamAlignment &a);
  
  void parseRuleLine(string line);

  friend ostream& operator<<(ostream &out, const AbstractRule &fr);
};

class MiniRulesCollection;

class MiniRules {
  
  friend class MiniRulesCollection;

  public:
  MiniRules() {}
  ~MiniRules() {}
  //MiniRules(const MiniRules * mr); // transfer defaults from one to another
    
  bool isValid(BamAlignment &a);
   
  void setIntervalTreeMap(string file);

  bool isOverlapping(BamAlignment &a);

  friend ostream& operator<<(ostream& out, const MiniRules &mr);


 private:

  bool m_whole_genome = false;
  GenomicRegionVector m_grv;
  GenomicIntervalTreeMap m_tree;
  string m_region_file;
  int m_level = -1;
  int m_width = 0;

  int pad = 0; // how much should we pad the region?

  vector<AbstractRule> m_abstract_rules;

  // rule applies to mate too
  bool m_applies_to_mate = false;

  // count the total number of valid reads
  int m_count = 0;
  
};

// a hierarchy of mini rules to operate on
class MiniRulesCollection {

 public: 
  MiniRulesCollection() {}
  ~MiniRulesCollection() {}
  MiniRulesCollection(string file);

  int isValid(BamAlignment &a);
  
  friend ostream& operator<<(ostream& out, const MiniRulesCollection &mr);
  
  void sendToBed(string file);

 private:

  vector<MiniRules*> m_rules;
  
};


#endif
