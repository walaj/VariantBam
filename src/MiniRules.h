#ifndef MINI_RULES_H
#define MINI_RULES_H

/* Define a set of rules for creating a variant bam. The syntax is:
   REGION_FILE   rule1:value1    rule1:value1
  
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

class MiniRulesCollection;

class MiniRules {
  
  friend class MiniRulesCollection;

  public:
  MiniRules() {}
  ~MiniRules() {}
  MiniRules(const MiniRules * mr); // transfer defaults from one to another
    
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

  // full include and full exclude
  bool m_full_include = false;
  bool m_full_exclude = false;

  // rule applies to mate too
  bool m_applies_to_mate = false;

  // numeric properities of reads to filter on
  int m_isize = 800;
  int m_mapq = 0;
  int m_qual = 0;
  int m_minclip = 5;
  int m_length = 50;
  int m_nmlim = 20;

  // boolean properites of reads
  bool m_unmap = true;
  bool m_supp = true;
  bool m_hardclip = true;
  bool m_duplicate = true;
  bool m_failqc = true;
  bool m_exclude_n = false;

  //bool m_hardclip_DISC = true;
  //bool m_duplicate_DISC = true;
  //bool m_failqc_DISC = true;
  //bool m_exclude_n_DISC = false;
  //bool m_supp_DISC = false;
  int  m_mapq_DISC = 0;

  // if true, ALL discordant reads are kept, regardless of filters
  //bool m_force_disckeep = false; 

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
