#include "VariantBamReader.h"

using namespace std;

// Trim the sequence by removing low quality bases from either end
int32_t VariantBamReader::qualityTrimRead(int qualTrim, int32_t &startpoint, Read &r) {

    int endpoint = -1; //seq.length();
    startpoint = 0;
    int i = 0; 

#ifdef HAVE_HTSLIB    
    uint8_t * qual = bam_get_qual(r.get());
#endif

    // get the start point (loop forward)
#ifdef HAVE_HTSLIB
    while(i < r->core.l_qseq) {
      int ps = qual[i];
#endif
#ifdef HAVE_BAMTOOLS
      int thislen = r->Qualities.length();
      while(i < thislen) {
      int ps = char2phred(r->Qualities[i]);
#endif
      if (ps >= qualTrim) {
          startpoint = i;
          break;
	}
	i++;
    }

    // get the end point (loop backwards)
#ifdef HAVE_HTSLIB
    i = r->core.l_qseq - 1; //seq.length() - 1;
#endif
#ifdef HAVE_BAMTOOLS
    i = r->Qualities.length() - 1; //core.l_qseq - 1; //seq.length() - 1;
#endif
    while(i >= 0) {
#ifdef HAVE_HTSLIB
      int ps = qual[i];
#endif
#ifdef HAVE_BAMTOOLS
      int ps = char2phred(r->Qualities[i]);
#endif
        if (ps >= qualTrim) { //ps >= qualTrim) {
	  endpoint = i + 1; // endpoint is one past edge
          break;
	}
	i--;
    }

    // check that they aren't all bad
    if (startpoint == 0 && endpoint == -1) {
      //trimmed_seq = 0; //trimmed_seq = "";
      //seq = "";
      //qual = "";
      return 0;
    }

    return (endpoint - startpoint);

}

// set the bam region
bool VariantBamReader::setBamRegion(GenomicRegion gp) {
  m_region = gp;

#ifdef HAVE_BAMTOOLS
  // set the region
  if (!m_reader->SetRegion(m_region.chr, m_region.pos1, m_region.chr, m_region.pos2)) {
    std::cerr << "Error: Failed to set region: " << gp << endl; 
    exit(EXIT_FAILURE);
  }
#endif 

#ifdef HAVE_HTSLIB
  //HTS set region
  if (!idx)
    idx = hts_idx_load(m_bam.c_str(), HTS_FMT_BAI);
  hts_itr = sam_itr_queryi(idx, gp.chr, gp.pos1, gp.pos2);
  if (!hts_itr) {
    std::cerr << "Error: Failed to set region: " << gp << endl; 
    exit(EXIT_FAILURE);
  }
#endif

  return true;
}

// closes the BamWriter and makes an index file
void VariantBamReader::MakeIndex() {

  /*
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
  */
}

// make a new object and put the reader and writer on the heap.
// this is also opens the files for reading/writing and checks
// that they are readable/writable
VariantBamReader::VariantBamReader(string in, string out, MiniRulesCollection *mr, int verbose) {

  m_mr = mr;
  m_bam = in;
  m_out = out;
  m_verbose = verbose;

#ifdef HAVE_BAMTOOLS
  // open the reader
  m_reader = new BamTools::BamReader();
  if (!m_reader->Open(in)) {
    cerr << "Error: Cannot open " << in << " for reading" << endl;
    exit(EXIT_FAILURE);
  }
#endif

#ifdef HAVE_HTSLIB
  // HTS open the reader
  const char rflag = 'r';
  fp = bgzf_open(in.c_str(), &rflag); 

  if (!fp) {
    cerr << "Error using HTS reader on opening " << in << endl;
    exit(EXIT_FAILURE);
  }
  br = bam_hdr_read(fp);

  if (!br) {
    cerr << "Error using HTS reader on opening " << in << endl;
    exit(EXIT_FAILURE);
  }
#endif
  // HTS region
  ///////////////
  //hts_idx_t * idx = hts_idx_load(in.c_str(), HTS_FMT_BAI);                                                                                                                                                                
  //hts_itr_t * hts_itr = sam_itr_queryi(idx, 3, 50000, 60001);                                                                                                                                                       
  //void* dum;                                                                                                                                                                                                        
                                                                                                                                                                                                                    
  //bam1_t * r = bam_init1();                                                                                                                                                                                         
  //int countr = 0;                                                                                                                                                                                                   
  //while (hts_itr_next(fp, hts_itr, r, dum) > 0) {                                                                                                                                                                   
  //  std::cout << r->core.pos << " " << r->core.tid << std::endl;                                                                                                                                                    
  //}                                                    
  /////////////////////////

  // get the index
#ifdef HAVE_BAMTOOLS
  if (!m_reader->LocateIndex()) {
    
    // try finding it manually
    string bai = in;
    if (!m_reader->OpenIndex(bai + ".bai")) {
      bai = SnowUtils::scrubString(bai, ".bam");
      bai += ".bai";
      if (!m_reader->OpenIndex(bai)) {
	cerr << "Error: Cannot locate index file for " << in << endl;
	exit(EXIT_FAILURE);
      }
    }

  }

  // open the writer
  m_writer = new BamTools::BamWriter();
  if (!m_writer->Open(out, m_reader->GetHeaderText(), m_reader->GetReferenceData())) {
    cerr << "Error: Cannot open BAM for writing " << out << endl;
    exit(EXIT_FAILURE);
  }
#endif

#ifdef HAVE_HTSLIB
  // hts open the writer
  fop = sam_open(out.c_str(), "wb");
  if (!fop) {
    cerr << "Error: Cannot open BAM for writing " << out << endl;
    exit(EXIT_FAILURE);
  }

  // hts write the header

  sam_hdr_write(fop, br);
#endif

}


void VariantBamReader::printMessage(const ReadCount &rc_main, const Read &r) const {

  char buffer[100];
  string posstring = SnowUtils::AddCommas<int>(r_pos(r));
  sprintf (buffer, "Reading read %11s at position %2s:%-11s. Kept %11s (%2d%%) [running count across whole BAM]",  
	   rc_main.totalString().c_str(), GenomicRegion::chrToString(r_id(r)).c_str(), posstring.c_str(),  
	   rc_main.keepString().c_str(), rc_main.percent());
  
  printf ("%s | ",buffer);
  SnowUtils::displayRuntime(start);
  cout << endl;
  
}

void VariantBamReader::writeVariantBamFromHash() {

#ifdef HAVE_BAMTOOLS
  BamTools::BamAlignment a;
  twopass = false;

  ReadCount rc_main;
  /*
  while (m_reader->GetNextAlignment(a)) {

    rc_main.total++;
    if (m_hash.count(a.Name)) {
      rc_main.keep++;
      saveAlignment(a);
    }

    if (m_verbose > 0 && rc_main.total % 500000 == 0) 
      printMessage(rc_main, a);
    
  }
  */
#else
  /*
  void* dum;
  int count = 0, keep_count = 0;
  for (;;) {
    bam1_t * b = bam_init1();
    if (hts_itr == 0) { // whole genome
      if (bam_read1(fp, b) < 0)
	break;
    } else { // region
      if (hts_itr_next(fp, hts_itr, b, dum) <= 0)
	break;
    }
    
    if (m_hash.count(string(bam_get_qname(b)))) {
      keep_count++;
      
      // remove tags if need
      uint8_t * p;
      if ( (p = bam_aux_get(b, "OQ")) ) bam_aux_del(b, p);
      if ( (p = bam_aux_get(b, "R2")) ) bam_aux_del(b, p);
      if ( (p = bam_aux_get(b, "Q2")) ) bam_aux_del(b, p);
      
      sam_write1(fop, br, b); 
    }
    
    if (++count % 2000000 == 0) {
      char buffer[100];
      string posstring = SnowUtils::AddCommas<int>(b->core.pos);
      sprintf (buffer, "Writing read at position %2s:%-11s from hash. Kept %11d of %11d",  
	       GenomicRegion::chrToString(b->core.tid).c_str(), posstring.c_str(),  
	       keep_count, count);
      
      printf ("%s | ",buffer);
      SnowUtils::displayRuntime(start);
      cout << endl;
    }

    bam_destroy1(b); // its written, so remove from heap

  }
  */
#endif
  

}

void VariantBamReader::printRuleCounts(unordered_map<string, size_t> &rm) const {

  size_t total = 0;
  for (auto& i : rm)
    total += i.second;
  for (auto& i : rm) {
    cout << "  " << i.first << ":" << i.second << "(" << SnowUtils::percentCalc<size_t>(i.second, total) << "%)" << endl;
  }
  
}

// save alignment
void VariantBamReader::saveAlignment(Read &r) {

  if (!twopass) {

    // remove tags if need
    r_remove_tag(r, "R2");
    r_remove_tag(r, "OQ");
    r_remove_tag(r, "Q2");

#ifdef HAVE_HTSLIB
    sam_write1(fop, br, r.get()); 
#endif

#ifdef HAVE_BAMTOOLS
    m_writer->SaveAlignment(*(r.get()));
#endif
    //bam_destroy1(b); // its written, so remove from heap

  } else {
    //char * s = bam_get_qname(b);
    //m_hash[string(s)] = true;
    //bam_destroy1(b);
  }

}





bool VariantBamReader::writeVariantBam(BamQC &qc, ReadVec &bav) {

  ReadCount rc_main;
  ReadCount rc_this;

  int pileup = 0;
  
  unordered_map<string, size_t> rule_count;

  ReadVec bam_buffer;

  for (;;) { 
    
    Read r;
    GET_READ(r);

    /*    bam1_t * b = bam_init1();
    if (hts_itr == 0) { // whole genome
      if (bam_read1(fp, b) < 0)
	break;
    } else { // region
      if (hts_itr_next(fp, hts_itr, b, dum) <= 0)
	  break;
	  }*/
    
    // pre-process it
    rc_this.total++;
    rc_main.total++;
    
    // decide what to do
    // check if read passes rules. 
    string rule_pass = m_mr->isValid(r);
    
    // build the qc
    if (qc.use)
      qc.addRead(r);
    
    if ( rule_pass != "" /*&& !qc_only*/ ) {
      
      auto ff = rule_count.find(rule_pass);
      if (ff != rule_count.end())
	ff->second++;
      else
	rule_count[rule_pass] = 1;
      
      // keep track of pile
      if (r_mapq(r) == 0) 
	pileup++;
      
      r_add_Z_tag(r, "RL", rule_pass);
      //bam_aux_append(b, "RL", 'Z', rule_pass.length()+1, (uint8_t*)rule_pass.c_str()); // need +1. Dunno why
      
      //bam_buffer.push_back(bam_dup1(b));
      bam_buffer.push_back(r);
      
      rc_this.keep++;
      rc_main.keep++;
      
      size_t buffer_lim = 100;
      
      // deal with bam buff
      if (bam_buffer.size() >= buffer_lim) {
	
	// check if bad region
	int buf_width = r_pos(bam_buffer.back()) - r_pos(bam_buffer[0]);
	//int buf_width = bam_buffer.back()->core.pos - bam_buffer[0]->core.pos;
	if (pileup >= buffer_lim * 0.8 && buf_width <= 40) { // TODO option to turn off
	  dumpBuffer(bam_buffer, bav, 0);
	  if (m_verbose > 2)
	    cout << "Detected mapq 0 pileup of " << pileup << " at " << (r_id(r)+1) << ":" << r_pos(bam_buffer[0]) << "-" << r_pos(bam_buffer.back()) << endl;
	}
	// it's OK or its in full region
	else if (bam_buffer.size() >= buffer_lim) {
	  dumpBuffer(bam_buffer, bav, -1);
	}
	
	pileup = 0;
	
      } // end buffer check
    } // end save read checking
    
      // print the message every 500k reads and check 
      // if we should continue or quit
    if (m_verbose > 0 && rc_main.total % 2000000 == 0) {
      
      printMessage(rc_main, r);
      
      if (m_verbose > 1)
	printRuleCounts(rule_count);
      
      // zero the counters
      rc_this = ReadCount();
      
      // kill if seen 25m reads, and it's looking bad
      int perclimit = 50;
      if (rc_main.percent() >= perclimit && rc_main.total > 25000000 && false) { 
	cerr << "This is a a really bad BAM after checking out 25m+ reads. Killing job. Percent weird reads: " << rc_main.percent() << " is above limit of " << perclimit << endl;
	cerr << "Reading in region" << m_region << endl;
	exit(EXIT_FAILURE);
      }
    }
    
    //bam_destroy1(b);
    //b = bam_init1();
    //b = bam_init1();
    
  } // end read while loop
  
    // print the final message
  if (m_verbose > 0)
    printMessage(rc_main, bam_buffer.back());
  
  // write the final buffer
  dumpBuffer(bam_buffer, bav, -1);
  
  
  //bam_destroy1(b); // destroy the one we created to be used for next read (which never came)
  return true;
}
  
  
void VariantBamReader::dumpBuffer(ReadVec &buff, ReadVec &store, int mapq) {

  for (auto& i : buff) { // = bam_buffer.begin(); it != bam_buffer.end(); it++) {
    if (r_mapq(i) > mapq) {
      if (true)
	saveAlignment(i);
      else 
	store.push_back(i);
    } //else {
      //bam_destroy1(i);
    //}
  }
  
  buff.clear();
}
