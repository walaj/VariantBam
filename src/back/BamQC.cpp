#include "BamQC.h"

#define MAX_MAPQ 60
#define MAX_NM 100
#define MAX_ISIZE 2000
#define MAX_CLIP 100
#define MAX_AS 100
#define MAX_XP 100
#define MAX_LEN 250
#define MAX_PHRED 60

template <typename T> void printQCVec(ostream &out, const vector<T> &vec) {
  for (auto it = vec.begin(); it != vec.end(); it++)
    out << "," << *it;
  return;
}

// instantiate a new read group
BamQCReadGroup::BamQCReadGroup() {

  mapq = vector<size_t>(MAX_MAPQ + 1,0); // 60 elems of 0
  nm   = vector<size_t>(MAX_NM + 1,0); // 101 elems of 0
  isize= vector<size_t>(MAX_ISIZE + 1,0); // (everything above 2000 is inter)
  clip = vector<size_t>(MAX_CLIP + 1,0); 
  as   = vector<size_t>(MAX_AS + 1,0);
  xp   = vector<size_t>(MAX_XP + 1,0);
  len  = vector<size_t>(MAX_LEN + 1,0);
  phred= vector<size_t>(MAX_PHRED+1, 0);
}

// make the output
std::ostream& operator<<(std::ostream& out, const BamQC& qc) {

  for (auto it = qc.map.begin(); it != qc.map.end(); it++)
    out << "READGROUP:" << it->first << endl << it->second << endl;
  return out;
}

// make the output
std::ostream& operator<<(std::ostream& out, const BamQCReadGroup& rg) {

  out << "total," << rg.num_reads << endl;
  out << "unmap," << rg.unmap << endl;
  out << "qcfail," << rg.qcfail << endl;
  out << "duplicate," << rg.duplicate << endl;
  out << "supplementary," << rg.supp << endl;

  out << "mapq";
  printQCVec<size_t>(out, rg.mapq);
  out << endl;

  out << "nm";
  printQCVec<size_t>(out, rg.nm);
  out << endl;

  out << "isize";
  printQCVec<size_t>(out, rg.isize);
  out << endl;

  out << "as";
  printQCVec<size_t>(out, rg.as);
  out << endl;
  
  out << "xp";
  printQCVec<size_t>(out, rg.xp);
  out << endl;

  out << "clip";
  printQCVec<size_t>(out, rg.clip);
  out << endl;

  out << "len";
  printQCVec<size_t>(out, rg.len);
  out << endl;

  out << "phred";
  printQCVec<size_t>(out, rg.phred);

  return out;
}

// add an addional alignment
void BamQC::addRead(Read &r) {

  try {

      string rgroup;
      r_get_Z_tag(r, "RG", rgroup);
      
      int this_isize = r_isize(r);
      this_isize = (r_mid(r) != r_id(r) || this_isize > 2000) ? 2000 : this_isize;
      
      // get clip num
      unsigned clipnum = 0;
      r_get_clip(r, clipnum);
      //int clipnum = VariantBamReader::getClipCount(a);
      
      // get the mean phred quality
      //size_t i = 0;
      //int phred = 0;
      //while(i < a.Qualities.length()) {
      //  phred += char2phred(a.Qualities[i]);
//	i++;
 //     }
   //   if (a.Qualities.length() > 0)
	//phred = static_cast<int>(floor(static_cast<float>(phred) / a.Qualities.length()));

      // get the NM tag
      int32_t nm;
      r_get_int32_tag(r, "NM", nm);
      //if (a.GetTag("NM", nm)) {} else { nm = 0; }
      int nmr = nm;

      int mapqr = r_mapq(r);

      // discordant
      bool FR_f = !r_is_rev(r) && (r_pos(r) < r_mpos(r)) && (r_id(r) == r_mid(r)) && r_is_mrev(r) && r_is_pmapped(r);
      bool FR_r = r_is_rev(r) && (r_pos(r) > r_mpos(r)) && (r_id(r) == r_mid(r)) && !r_is_mrev(r) && r_is_pmapped(r);
      bool FR = FR_f || FR_r;
      if (r_isize(r) > 0 && this_isize <= 2000 && r_is_pmapped(r) && FR ) // only count "proper" reads 
	map[rgroup].isize[this_isize]++;
      
      // all the rest
      //map[rgroup].xp[xp]++;
      map[rgroup].len[min(r_length(r), MAX_LEN)]++;
      //map[rgroup].as[as]++;
      map[rgroup].clip[min((int)clipnum,MAX_CLIP)]++;
      //map[rgroup].phred[min(phred, MAX_PHRED)]++;
      map[rgroup].num_reads++;
      map[rgroup].mapq[min(mapqr, MAX_MAPQ)]++;
      map[rgroup].nm[min(nmr, MAX_NM)]++;
      
      if (!r_is_mapped(r))
	map[rgroup].unmap++;
      if (r_is_qc_fail(r)) 
	map[rgroup].qcfail++;
      if (r_is_dup(r))
	map[rgroup].duplicate++;
      if (!r_is_primary(r))
	map[rgroup].supp++;
      
  } catch (...) {
    cerr << "Failed at adding to QC" << endl;
    //cerr << "Readgroup " << "NM " << nm << " mapq " << a.MapQuality << 
    //" xp " << xp << 
    //" len " << a.Length <<
    //" as " << as << 
    //" phred " << phred << endl;
  }
  
}

