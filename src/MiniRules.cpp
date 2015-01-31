#include "MiniRules.h"
#include "VariantBamReader.h"
#include <regex>

using namespace std;
using namespace BamTools;

bool MiniRules::isValid(BamAlignment &a) {

  if (m_full_include)
    return true;
  if (m_full_exclude)
    return false;

  // check for unmapped read or mate
  bool unmap = !a.IsMapped() || !a.IsMateMapped();
  unmap = unmap && m_unmap;
  
  // check if is discordant
  bool discr = (abs(a.InsertSize) > m_isize) || (a.RefID != a.MateRefID) || unmap; 

  // check that read orientation is as expected
  if (!discr) {
    bool FR_f = !a.IsReverseStrand() && (a.Position < a.MatePosition) && (a.RefID == a.MateRefID) &&  a.IsMateReverseStrand();
    bool FR_r =  a.IsReverseStrand() && (a.Position > a.MatePosition) && (a.RefID == a.MateRefID) && !a.IsMateReverseStrand();
    bool FR = FR_f || FR_r;
    discr = discr || !FR;
  }

  // check for hard clips
  bool hardclip_pass = true;
  if (!m_hardclip) {// check that we want to chunk hard clip
    for (auto cig : a.CigarData)
      if (cig.Type == 'H')
	hardclip_pass = false;
  }

  // filter the discordant-specific params
  discr = discr && m_isize > 0;
  discr = discr && (m_supp || a.IsPrimaryAlignment()); // check its primary
  discr = discr && (m_failqc || !a.IsFailedQC()); // check that its QC is ok
  discr = discr && (m_duplicate || !a.IsDuplicate());
  discr = discr && a.MapQuality >= m_mapq_DISC;
  discr = discr && (hardclip_pass || m_hardclip);

  if (discr)
    return true;

  // check for mapq pass
  bool mapq_pass = a.MapQuality >= m_mapq;
    
  // check for flag quality
  bool flag = (!a.IsFailedQC() || !m_failqc) && (!a.IsDuplicate() && !m_duplicate);

  // check for supplementatry
  bool supp  = a.IsPrimaryAlignment() || m_supp;

  // check for clipping
  unsigned clipnum = VariantBamReader::getClipCount(a);
  
  // trim for quality
  string trimmed_bases = a.QueryBases;
  string trimmed_quals = a.Qualities;
  VariantBamReader::qualityTrimRead(m_qual, trimmed_bases, trimmed_quals); 
  int new_clipnum = max(0, static_cast<int>(clipnum - (a.Length - trimmed_bases.length())));

  // check for QUALITY CLIP
  bool clip  = (new_clipnum >= m_minclip) || (m_minclip == 0);
  
  // check for QUALITY LENGTH
  int newlen = trimmed_bases.length();
  bool len = newlen >= m_length;

  // check to remove reads with N
  bool npass = a.QueryBases.find("N") == string::npos;
  npass = npass || !m_exclude_n;
  
  // check NM tag
  uint32_t nm;
  if (a.GetTag("NM", nm)) {} else { nm = 0; }
  int intnm = nm;
  bool nmpass = intnm <= m_nmlim;

  // accumulate all the information to make a judgement
  bool qual_read = len && nmpass && npass && hardclip_pass && mapq_pass && flag && supp; 
  bool save_read = clip && qual_read;

  // update the read counter
  if (save_read)
    m_count++;

  return save_read;
}

// check whether a BamAlignment (or optionally it's mate) is overlapping the regions
// contained in these rules
bool MiniRules::isOverlapping(BamAlignment &a) {

  // if this is a whole genome rule, it overlaps
  if (m_whole_genome)
    return true;

  // TODO fix a.MatePosition + a.Length is using wrong length

  // check whether a read (or maybe its mate) hits a rule
  GenomicIntervalVector grv;
  if (m_tree.count(a.RefID) == 1) // check that we have a tree for this chr
    m_tree[a.RefID].findOverlapping(a.Position, a.Position + a.Length, grv);
  if (m_tree.count(a.MateRefID) == 1 && m_applies_to_mate) // check that we have a tree for this chr
    m_tree[a.MateRefID].findOverlapping (a.MatePosition, a.MatePosition + a.Length, grv);
  return grv.size() > 0;
  
}

// checks which rule a read applies to (using the hiearchy stored in m_rules).
// if a read does not satisfy a rule it is excluded.
int MiniRulesCollection::isValid(BamAlignment &a) {

  size_t which_rule = 0;
  
  // find out which rule it is a part of
  // lower number rules dominate
  for (auto it : m_rules) {
    if (it->isOverlapping(a))
      if (it->isValid(a))
	break;
    which_rule++;
  }
  
  // isn't in a rule or it never satisfied one. Remove
  if (which_rule >= m_rules.size())
    return 0; 

  int level = which_rule + 1;// rules start at 1
  return level; 
  
}

// convert a region BED file into an interval tree map
void MiniRules::setIntervalTreeMap(string file) {
  
  m_region_file = file;
  GenomicRegionVector grv = GenomicRegion::regionFileToGRV(file, 0); // 0 is pad
  m_grv = GenomicRegion::mergeOverlappingIntervals(grv); // 0 is pad
  sort(m_grv.begin(), m_grv.end());

  // set the width
  for (auto it : m_grv)
    m_width += it.width();
 
  size_t grv_size = m_grv.size();
  if (grv_size == 0) {
    cerr << "Warning: No regions dected in file: " << file << endl;
    return;
  }

  m_tree = GenomicRegion::createTreeMap(m_grv);
  return;
}

// constructor to make a MiniRulesCollection from a rules file.
// This will reduce each individual BED file and make the 
// GenomicIntervalTreeMap
MiniRulesCollection::MiniRulesCollection(string file) {

  // parse the rules file
  vector<string> region_files;
  ifstream iss_rules(file.c_str());
  if (!iss_rules) {
    cerr << "Rules file " << file << " is not able to be opened" << endl;
    exit(EXIT_FAILURE);
  }
  
  // loop through the rules file and grab the rules
  string line;
  int level = 1;
  MiniRules * mr_all = new MiniRules();
  while(getline(iss_rules, line, '\n')) {

    if (line.length() == 0)
      break;
    //exclude commentes
    if (line.at(0) != '#') {

    MiniRules * mr = new MiniRules(mr_all);
    int counter = 0;    
    istringstream iss_line(line);
    string val;
    bool isall = false;

    // loop through the rules
    while(getline(iss_line, val, ',')) {

      // load the region file
      if (counter == 0) {
	if (val == "NA" || val == "WG") 
	  mr->m_whole_genome = true;
	else if (val == "ALL")
	  isall = true;
	else 
	  mr->setIntervalTreeMap(val);
      } 

      // check for rules
      if (counter > 0) {
	regex mapqreg("mapq:(.*)");  
	regex isizereg("isize:(.*)");  
	regex clipreg("clip:(.*)");  
	regex mapqDISCreg("mapq-discordant:(.*)");
	regex phredreg("phred-qual:(.*)");

	smatch match;
	if (regex_search(val, match, mapqreg))
	  try { mr->m_mapq = stoi(match[1].str()); } catch (...) { cerr << "Error trying to convert MAPQ to num. Setting to 0. Val: " << val << " match " << match[1].str() << endl; mr->m_mapq = 0; }
	else if (regex_search(val, match, isizereg))
	  try { mr->m_isize = stoi(match[1].str()); } catch (...) { cerr << "Error trying to convert ISIZE to num. Setting to -1. Val: " << val << " match " << match[1].str() << endl; mr->m_isize = 0; }
	else if (regex_search(val, match, clipreg))
	  try { mr->m_minclip = stoi(match[1].str()); } catch (...) { cerr << "Error trying to convert CLIP to num. Setting to 0"; mr->m_minclip = 0; }
	else if (regex_search(val, match, mapqDISCreg))
	  try { mr->m_mapq_DISC = stoi(match[1].str()); } catch (...) { cerr << "Error trying to convert mapqDISC to num. Setting to 0"; mr->m_mapq_DISC = 0; }
	else if (regex_search(val, match, phredreg))
	  try { mr->m_qual = stoi(match[1].str()); } catch (...) { cerr << "Error trying to convert PHRED to num. Setting to 0"; mr->m_qual = 0; }
	else if (val == "all")
	  mr->m_full_include = true;
	else if (val == "none")
	  mr->m_full_exclude = true;
	else if (val == "mate")
	  mr->m_applies_to_mate = true;
	else if (val == "keep-unmapped-and-mapped-mate")
	  mr->m_unmap = true;
	else if (val == "remove-hardclip")
	  mr->m_hardclip = false;
	else if (val == "remove-supplementary")
	  mr->m_supp = false;
	else if (val == "remove-duplicate")
	  mr->m_duplicate = false;
	else if (val == "remove-failqc")
	  mr->m_failqc = false;
	else if (val == "remove-reads-with-N")
	  mr->m_exclude_n = true;
	/*else if (val == "remove-hardclip-DISCORDANT")
	  mr->m_hardclip_DISC = false;
	else if (val == "remove-supplementary-DISCORDANT")
	  mr->m_supp = false;
	else if (val == "remove-duplicate-DISCORDANT")
	  mr->m_duplicate = false;
	else if (val == "remove-failqc-DISCORDANT")
	  mr->m_failqc = false;
	*/
      }

      counter++;
    } // end , parse
    
    // if it's ALL, then set as default but dont make rule
    if (isall) {
      mr_all = mr;
    } else {
      mr->m_level = level++;
      m_rules.push_back(mr);
    }
    } //end comment check
  } // end \n parse
  delete mr_all;
}

// print the MiniRulesCollectoin
ostream& operator<<(ostream &out, const MiniRulesCollection &mr) {

  out << "---- MiniRulesCollection ----" << endl;
  for (auto it : mr.m_rules)
    out << (*it);

  return out;

}

// print a MiniRules information
ostream& operator<<(ostream &out, const MiniRules &mr) {
   
  string file_print = mr.m_whole_genome ? "WHOLE GENOME" : VarUtils::getFileName(mr.m_region_file);
  out << "LEVEL " << mr.m_level << " Rules applying to: " << file_print << endl;
  cout << "--Include Mate: " << (mr.m_applies_to_mate ? "ON" : "OFF") << endl;
  if (mr.m_full_include) {
    out << "   Including all reads and mates"  << endl;
  } else if (mr.m_full_exclude) {
    out << "   Excluding all reads and mates"  << endl;
  } else {
    out << "Non-discordant read filters: " << endl;
    out << "   Min MAPQ:             " << mr.m_mapq << endl;
    out << "   Min CLIP:             " << mr.m_minclip << endl;
    out << "   Phred trim:           " << mr.m_qual << endl;
    out << "   Remove Hardclip:      " << (mr.m_hardclip ? "OFF" : "ON") << endl;
    out << "   Remove Duplicate:     " << (mr.m_duplicate ? "OFF" : "ON") << endl;
    out << "   Remove QCFail:        " << (mr.m_failqc ? "OFF" : "ON") << endl;
    out << "   Remove Supplementary: " << (mr.m_supp ? "OFF" : "ON") << endl;
    out << "   Min ISIZE:    " << mr.m_isize << endl;
    out << "Discordant read filters: " << endl;
    out << "   Min MAPQ:             " << mr.m_mapq_DISC << endl;
    out << "   Min ISIZE:            " << mr.m_isize << endl;
    //out << "   Remove Hardclip:      " << (mr.m_hardclip_DISC ? "OFF" : "ON") << endl;
    //out << "   Remove Duplicate:     " << (mr.m_duplicate_DISC ? "OFF" : "ON") << endl;
    //out << "   Remove QCFail:        " << (mr.m_failqc_DISC ? "OFF" : "ON") << endl;
    //out << "   Remove Supplementary: " << (mr.m_supp_DISC ? "OFF" : "ON") << endl;
  }
  out << "-- Size of rules region: " << (mr.m_whole_genome ? "WHOLE GENOME" : VarUtils::AddCommas<int>(mr.m_width)) << endl;
  return out;
}

// make a MiniRules that inherits flags from other MiniRules
MiniRules::MiniRules(const MiniRules * mr) {

  // full include and full exclude
  m_full_include = mr->m_full_include;
  m_full_exclude = mr->m_full_exclude;

  // rule applies to mate too
  m_applies_to_mate = mr->m_applies_to_mate;

  // numeric properities of reads to filter on
  m_isize = mr->m_isize;
  m_mapq = mr->m_mapq;
  m_qual = mr->m_qual;
  m_minclip = mr->m_minclip;
  m_length = mr->m_length;
  m_nmlim = mr->m_nmlim;

  m_unmap = mr->m_unmap;
  m_supp  = mr->m_supp; 
  m_hardclip = mr->m_hardclip; 
  m_duplicate = mr->m_duplicate; 
  m_failqc = mr->m_failqc; 
  m_exclude_n = mr->m_exclude_n; 

  m_mapq_DISC = mr->m_mapq_DISC;


}

// merge all of the intervals into one and send to a bed file
void MiniRulesCollection::sendToBed(string file) {

  ofstream out(file);
  if (!out) {
    cerr << "Cannot write BED file: " << file << endl;
    return;
  }
  out.close();

  // make a composite
  GenomicRegionVector comp;
  for (auto it : m_rules)
    comp.insert(comp.begin(), it->m_grv.begin(), it->m_grv.end()); 
  
  // merge it down
  GenomicRegionVector merged = GenomicRegion::mergeOverlappingIntervals(comp);

  // send to BED file
  GenomicRegion::sendToBed(merged, file);
  return;
}
