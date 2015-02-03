#include "VariantBamReader.h"

using namespace std;
using namespace BamTools;

// Trim the sequence by removing low quality bases from either end
void VariantBamReader::qualityTrimRead(int qualTrim, std::string &seq, std::string &qual) {

    assert(seq.size() == qual.size());

    int endpoint = -1; //seq.length();
    int startpoint = 0;
    int i = 0; 
 
    // get the start point (loop forward)
    while(i < (int)seq.length()) {
        int ps = char2phred(qual[i]);
        if (ps >= qualTrim) {
          startpoint = i;
          break;
	}
	i++;
    }

    // get the end point (loop backwards)
    i = seq.length() - 1;
    while(i >= 0) {
        int ps = char2phred(qual[i]);
        if (ps >= qualTrim) {
          endpoint = i + 1; // endpoint is one past edge
          break;
	}
	i--;
    }
    // check that they aren't all bad
    if (startpoint == 0 && endpoint == -1) {
      seq = "";
      qual = "";
      return;
    }

    // Clip the read
    seq =   seq.substr(startpoint, endpoint);
    qual = qual.substr(startpoint, endpoint);

    return;

}

// obtain the clipped count
unsigned VariantBamReader::getClipCount(BamAlignment a) {
  
  std::vector<int> clipSize;
  std::vector<int> readPos;
  std::vector<int> genPos;
  a.GetSoftClips(clipSize, readPos, genPos, false);
  
  // get the clip number
  unsigned clipnum = 0;
  for(std::vector<int>::iterator j=clipSize.begin(); j != clipSize.end(); ++j)
    clipnum += *j;
  return clipnum;
}


//bool VariantBamReader::writeVariantBam(BamQC &qc, bool qc_only) {
bool VariantBamReader::writeVariantBam() {

  int keep_counter = 0;
  int total = 0;
  int keep_counter_MAIN = 0;
  int total_MAIN = 0;

  int mapq0_keep_counter = 0;
  int discordant_keep_counter = 0;
  int n_keep_counter = 0;
  int clipped_keep_counter = 0;

  int mapq0_counter = 0;
  int discordant_counter = 0;
  int n_counter = 0;
  int clipped_counter = 0;

  int pileup = 0;
  
  BamTools::BamAlignment a;

  BamAlignmentVector bam_buffer;
  vector<int> mapq_buffer;

  while (m_reader->GetNextAlignmentCore(a)) {

    total++;
    total_MAIN++;
        
    if (m_verbose > 0 && total % 500000 == 0) {

      char buffer[100];
      int perc  = VarUtils::percentCalc<int>(keep_counter_MAIN, total_MAIN); 
      string posstring = VarUtils::AddCommas<int>(a.Position);
      sprintf (buffer, "Reading read %11s at position %2s:%-11s. Kept %11s (%2d%%) [running count across whole BAM]",  VarUtils::AddCommas<int>(total_MAIN).c_str(), GenomicRegion::chrToString(a.RefID).c_str(), posstring.c_str(),  VarUtils::AddCommas<int>(keep_counter_MAIN).c_str(), perc);
      printf ("%s\n",buffer);
      char buffer2[100];
      sprintf(buffer2, "   Filter (%% of kept)  -- Reads with N (%2d%%), Mapq0 (%2d%%), Discordant (%2d%%), Clipped (%2d%%)", 
	      VarUtils::percentCalc<int>(n_keep_counter, keep_counter), 
	      VarUtils::percentCalc<int>(mapq0_keep_counter, keep_counter), 
	      VarUtils::percentCalc<int>(discordant_keep_counter, keep_counter), 
	      VarUtils::percentCalc<int>(clipped_keep_counter, keep_counter));
      char buffer3[100];
      sprintf(buffer3, "   Filter (%% of total) -- Reads with N (%2d%%), Mapq0 (%2d%%), Discordant (%2d%%), Clipped (%2d%%)", 
	      VarUtils::percentCalc<int>(n_counter, total), 
	      VarUtils::percentCalc<int>(mapq0_counter, total), 
	      VarUtils::percentCalc<int>(discordant_counter, total), 
	      VarUtils::percentCalc<int>(clipped_counter, total));
      if (m_verbose > 1) {
	printf("%s\n", buffer2);
	printf("%s\n", buffer3);
      }

      // zero the counters
      mapq0_keep_counter = 0;
      discordant_keep_counter = 0;
      n_keep_counter = 0;
      clipped_keep_counter = 0;
      
      mapq0_counter = 0;
      discordant_counter = 0;
      n_counter = 0;
      clipped_counter = 0;
      keep_counter = 0;
      total = 0;

      // kill if seen 50m reads, and it's looking bad
      int perclimit = 50;
      if (perc >= perclimit && total > 25000000) { 
	cerr << "This is a a really bad BAM after checking out 25m+ reads. Killing job. Percent weird reads: " << perc << " is above limit of " << perclimit << endl;
	cerr << "Reading in region" << m_region << endl;
	exit(EXIT_FAILURE);
      }
    }
    string rule_pass = m_mr->isValid(a);

    // build the qc
    /*
    try {
      string rgroup;
      if (!a.GetTag("RG",rgroup))
	cerr << "Failed to read rgroup" << endl;

      int this_isize = a.InsertSize;
      this_isize = (a.MateRefID != a.RefID || this_isize > 2000) ? 2000 : this_isize;

      assert(a.MapQuality <= 60);
      //assert(clipnum <= 101);
      //assert(as <= 101);
      //assert(xp <= 101);
      assert(a.Length <= 101 && a.Length  >= 0);
      //assert(phred <= 60 && phred  >= 0);
      //assert(nm <= 101);

      //qc.map[rgroup].nm[nm]++;
      qc.map[rgroup].mapq[a.MapQuality]++;
      //if (a.InsertSize > 0 && a.IsPaired() && (FR_f || FR_r) ) // only count "proper" reads
	//qc.map[rgroup].isize[this_isize]++;
      //qc.map[rgroup].xp[xp]++;
      //qc.map[rgroup].len[a.Length]++;
      //qc.map[rgroup].as[as]++;
      //qc.map[rgroup].clip[clipnum]++;
      //qc.map[rgroup].phred[phred]++;
      qc.map[rgroup].num_reads++;
      if (!a.IsMapped())
	qc.map[rgroup].unmap++;
      if (a.IsFailedQC()) 
	qc.map[rgroup].qcfail++;
      if (a.IsDuplicate())
	qc.map[rgroup].duplicate++;
      if (!a.IsPrimaryAlignment())
	qc.map[rgroup].supp++;
    } catch (...) {
      cerr << "Failed at adding to QC" << endl;
      //cerr << "Readgroup " << "NM " << nm << " mapq " << a.MapQuality << " xp " << xp << " len " << a.Length <<
      //	" as " << as << " phred " << phred << endl;
    }
    */

    if ( rule_pass != "" /*&& !qc_only*/ ) {

      mapq0_keep_counter += (a.MapQuality == 0 ) ? 1 : 0; 

      // keep track of pile
      if (a.MapQuality == 0) 
	pileup++;

      // add a tag to say which rule it pass
      a.AddTag("RL","Z",rule_pass);

      bam_buffer.push_back(a);
      keep_counter++;
      keep_counter_MAIN++;
      
      size_t buffer_lim = 100;
      // deal with bam buff
      if (bam_buffer.size() >= buffer_lim/* && !in_full_region*/) {
	// check if bad region
	int buf_width = bam_buffer.back().Position - bam_buffer[0].Position;
	if (pileup >= buffer_lim * 0.8 && buf_width <= 40) {
	  for (auto it = bam_buffer.begin(); it != bam_buffer.end(); it++) 
	    if (it->MapQuality > 0)
	      m_writer->SaveAlignment(*it);
	  if (m_verbose > 2)
	    cout << "Detected mapq 0 pileup of " << pileup << " at " << a.RefID+1 << ":" << bam_buffer[0].Position << "-" << bam_buffer.back().Position << endl;
	} 
	// it's OK or its in full region
	else if (bam_buffer.size() >= buffer_lim) {
	  for (auto it = bam_buffer.begin(); it != bam_buffer.end(); it++) 
	    m_writer->SaveAlignment(*it);
	}

	bam_buffer.clear();
	pileup = 0;

      } // end buffer check
      
    } // end save read checking

  } // end read while loop

  // write the final buffer
  for (auto it = bam_buffer.begin(); it != bam_buffer.end(); it++)
    m_writer->SaveAlignment(*it);
  
  // print the final message
  if (m_verbose > 0) {
    char buffer[100];
    int perc  = VarUtils::percentCalc<int>(keep_counter_MAIN, total_MAIN); 
    string posstring = VarUtils::AddCommas<int>(a.Position);
    sprintf (buffer, "Finished region at %20s. Kept %11s (%2d%%) [running count across whole BAM]",  m_region.toString().c_str(), VarUtils::AddCommas<int>(keep_counter_MAIN).c_str(), perc);
    printf ("%s\n",buffer);
    char buffer2[100];
    sprintf(buffer2, "   Filter (%% of kept)  -- Reads with N (%2d%%), Mapq0 (%2d%%), Discordant (%2d%%), Clipped (%2d%%)", 
	    VarUtils::percentCalc<int>(n_keep_counter, keep_counter), 
	    VarUtils::percentCalc<int>(mapq0_keep_counter, keep_counter), 
	    VarUtils::percentCalc<int>(discordant_keep_counter, keep_counter), 
	    VarUtils::percentCalc<int>(clipped_keep_counter, keep_counter));
    char buffer3[100];
    sprintf(buffer3, "   Filter (%% of total) -- Reads with N (%2d%%), Mapq0 (%2d%%), Discordant (%2d%%), Clipped (%2d%%)", 
	    VarUtils::percentCalc<int>(n_counter, total), 
	    VarUtils::percentCalc<int>(mapq0_counter, total), 
	    VarUtils::percentCalc<int>(discordant_counter, total), 
	    VarUtils::percentCalc<int>(clipped_counter, total));
    if (m_verbose > 1) {
      printf("%s\n", buffer2);
      printf("%s\n", buffer3);
    }
  }
  
  
  return true;

}

bool VariantBamReader::setBamRegion(GenomicRegion gp) {
  m_region = gp;

  // set the region
  if (!m_reader->SetRegion(m_region.chr, m_region.pos1, m_region.chr, m_region.pos2)) {
    std::cerr << "Error: Failed to set region: " << gp << endl; 
    exit(EXIT_FAILURE);
  }

  return true;
}

// closes the BamWriter and makes an index file
void VariantBamReader::MakeIndex() {

  m_writer->Close();
  
  // open the file 
  BamReader reader;
  if (!reader.Open(m_out)) {
    cerr << "Error: Could not open the output BAM to create index " << m_out << endl;
    exit(EXIT_FAILURE);
  }

  // create the index
  if (!reader.CreateIndex()) {
    cerr << "Error: Could not create the output BAM index for " << m_out << endl;
    exit(EXIT_FAILURE);
  }

  reader.Close();

}

// make a new object and put the reader and writer on the heap.
// this is also opens the files for reading/writing and checks
// that they are readable/writable
VariantBamReader::VariantBamReader(string in, string out, MiniRulesCollection *mr, int verbose) {

  m_mr = mr;
  m_bam = in;
  m_out = out;
  m_verbose = verbose;

  // open the reader
  m_reader = new BamReader();
  if (!m_reader->Open(in)) {
    cerr << "Error: Cannot open " << in << " for reading" << endl;
    exit(EXIT_FAILURE);
  }

  // get the index
  if (!m_reader->LocateIndex()) {
    cerr << "Error: Cannot locate index file for " << in << endl;
    exit(EXIT_FAILURE);
  }

  // open the writer
  m_writer = new BamWriter();
  if (!m_writer->Open(out, m_reader->GetHeaderText(), m_reader->GetReferenceData())) {
    cerr << "Error: Cannot open BAM for writing " << out << endl;
    exit(EXIT_FAILURE);
  }


}

