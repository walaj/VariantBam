#include <string>
#include <vector>
#include <getopt.h>
#include <iostream>
#include "GenomicRegion.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "BamQC.h"
#include "VarUtils.h"
#include "MiniRules.h"
#include "VariantBamReader.h"
#include "gzstream.h"
#include "reads.h"

//#define HAVE_BAMTOOLS

//
//bamt=/broad/software/free/Linux/redhat_5_x86_64/pkgs/pezmaster31_bamtools-6708a21
//./configure --with-bamtools=$bamt

using namespace std;
using namespace BamTools;

static const char *VARIANT_BAM_USAGE_MESSAGE =
"Usage: varbam -i <input.bam> -o <output.bam> [OPTIONS] \n\n"
"  Description: Process a BAM file for use with rearrangement variant callers by removing proper pairs and bad regions\n"
"\n"
" General options\n"
"  -v, --verbose                        Select verbosity level (0-4). Default: 1 \n"
"  -h, --help                           Display this help and exit\n"
"  -i, --input-bam                      BAM file to preprocess\n"
"  -o, --output-bam                     File to write final BAM to\n"
"  -q, --qc-only                        Loop through the BAM, but only to make the QC file\n"
"  -g, --region                         Regions (e.g. myvcf.vcf or WG for whole genome) or newline seperated subsequence file.  Applied in same order as -r for multiple\n"
"  -r, --rules                          A script for the rules. If specified multiple times, will be applied in same order as -g\n"
"  -f, --rules-script                   A file (script) for the rules and regions/sequences. If specified, ignores -r and -g flags.\n"
  //"  -k, --proc-regions-file              csv file of regions to proess reads from, of format chr,pos1,pos2 eg 11,2232423,2235000. Default: Whole genome\n"
"  -t, --two-pass                       If set, will pass through BAM second time and retrieve reads AND MATES of reads satisfying a rule. Currently unbounded memory usage...\n"
  //" Read Filters\n"
  //"  -w, --min-mapq                       Minimum mapping quality a read must have to be included. Default 0\n"
  //"  -c, --min-clip                       Minimum number of bases clipped bases to considered clipped. Default: 5\n"
  //"  -n  --nm-limit                       Skip reads that have more than nm mismatches (NM tag). Default: 15\n"
  //"  -p  --perc-limit                     Exit with failure if more than this percent of reads are weird (Must see over 50m). Default: 30\n"
  //"  -s, --isize                          Insert size threshold to consider discordant. Default: 800\n"
  //"  -m, --min-readlength                 Removes reads < this cutoff. If quality score trimming is implemented, evaluates after trimming. Default: 50\n"
  //"  -z, --min-phred                      When measuring clipping and read length, dont count bases with less than this quality. Set to 0 to turn off. Default: 4\n"
  //"  -e, --exclude-with-n                 Remove reads that have N in their sequence. Default: OFF\n"
  //"      --no-pileup-check                Don't perform the low mapq pileup check. Default is to exclude reads with a pileup of > 100 mapq0 reads within 200 bp.\n"
  //"      --skip-centromeres               Don't processes centromers at all. Default OFF\n"
" Optional Input\n"
"\n";

namespace opt {

  static string bam;
  static string out;
  //static int isize = 800;
  static size_t verbose = 1;
  //static int pad = 500;
  static string mutect_callstats = "";
  static string mutect2_regions = "";
  
  static string rules = "";
  static string rules_file = "";
  static int mapq = 0;
  static bool skip_cent = false;
  static string proc_regions = "";

  static bool twopass = false;

  //static int perc_limit = 30;

  //static bool no_pileup_check = false;
  //static bool qc_only = false;

  static vector<string> rules_vec;
  static vector<string> region_vec;
}

static const char* shortopts = "hv:qejb:i:o:w:n:p:s:r:m:z:c:k:u:x:f:g:tc:";
static const struct option longopts[] = {
  { "help",                       no_argument, NULL, 'h' },
  { "verbose",                    required_argument, NULL, 'v' },
  { "input-bam",                  required_argument, NULL, 'i' },
  { "output-bam",                 required_argument, NULL, 'o' },
  { "two-pass",                   required_argument, NULL, 't' },
  //{ "min-mapq",                 required_argument, NULL, 'w' },
  //{ "min-clip",                 required_argument, NULL, 'c' },
  //{ "nm-limit",                 required_argument, NULL, 'n' },
  //{ "perc-limit",                 required_argument, NULL, 'p' },
  //{ "isize",                 required_argument, NULL, 's' },
  { "qc-only",                 no_argument, NULL, 'q' },
  { "rules-file",                 required_argument, NULL, 'k' },
  { "rules",                 required_argument, NULL, 'r' },
  { "region",                 required_argument, NULL, 'g' },
  { "region-with-mates",                   required_argument, NULL, 'c' },
  { "proc-regions-file",                 required_argument, NULL, 'k' },
  //  { "min-length",                 required_argument, NULL, 'm' },
  // { "min-phred",                 required_argument, NULL, 'z' },
  //{ "exclude-with-n",                 no_argument, NULL, 'e' },
  { "no-pileup-check",                 no_argument, NULL, 'j' },
  { "skip-centromeres",                 no_argument, NULL, 'b' },
  { "mutect-callstats",               required_argument, NULL, 'u' },
  { "mutect2-regions",               required_argument, NULL, 'x' },
  { NULL, 0, NULL, 0 }
};

static struct timespec start;

// forward declare
void parseVarOptions(int argc, char** argv);

int main(int argc, char** argv) {

  // start the timer
  clock_gettime(CLOCK_MONOTONIC, &start);

  parseVarOptions(argc, argv);

  /*
  // convert the region / rules vecs into rules
  if (!VarUtils::existTest(opt::rules_file) && opt::rules_vec.size() == 0) {
    cerr << "Rules / regions not read in. Rules file not found: " << opt::rules_file << endl;;
    exit(EXIT_FAILURE);
  }

  // add the "region" and "rules"
  for(size_t i = 0; i < min(opt::rules_vec.size(), opt::region_vec.size()); i++) {
    if (i == 0)
      opt::rules = "region@" + opt::region_vec[0] + "%" + opt::rules_vec[0];
    else
      opt::rules += "%region@" + opt::region_vec[i] + "%" + opt::rules_vec[i];
      }*/

  if (opt::verbose > 0) {
    cout << "Input BAM:  " << opt::bam << endl;
    cout << "Output BAM: " << opt::out << endl;
    cout << "Input rules and regions: " << opt::rules_file << endl;
    cout << "Input proc regions file: " << opt::proc_regions << endl;
    cout << "TWO-PASS solution?:      " << (opt::twopass ? "ON" : "OFF") << endl;
#ifdef HAVE_BAMTOOLS
    cout << "RUNNING AS BAMTOOLS" << endl;
#else
    cout << "RUNNING AS HTS" << endl;
#endif
  }

  // make the mini rules collection from the rules file
  // this also calls function to parse the BED files
  MiniRulesCollection * mr;
  //if (opt::rules.length() > 0)
  //  mr = new MiniRulesCollection(opt::rules);
  //else
  mr = new MiniRulesCollection(opt::rules_file);
  // check that it's valid
  if (mr->size() == 0) {
    cerr << "No rules or regions specified. Provide via a script and pass with -r flag." << endl;
    exit(EXIT_FAILURE);
  }

  assert(mr->numRules() > 0);

  if (opt::verbose > 0)
    cout << (*mr);
 
    // make a bed file
  if (opt::verbose > 0)
    cout << "...sending merged regions to BED file" << endl;
  mr->sendToBed("merged_rules.bed");

  // parse the proc region file
  bool runWholeGenome = mr->hasWholeGenome();
  GenomicRegionVector grv_proc_regions;
  if (opt::proc_regions.length() > 0)
    grv_proc_regions = GenomicRegion::regionFileToGRV(opt::proc_regions, 0);
  if (opt::proc_regions.length() > 0)
    assert(grv_proc_regions.size() > 0);
  
  // run just the regions in the mask file
  if (grv_proc_regions.size() > 0) {
    sort(grv_proc_regions.begin(), grv_proc_regions.end());      
    if (opt::verbose > 0) {
      if (opt::verbose > 1)
	for (auto it = grv_proc_regions.begin(); it != grv_proc_regions.end(); it++)
	  cout << "   Process-region: " << *it << endl;
    }
  } 
  // run everything, but skip the centromerss
  else if (opt::skip_cent && runWholeGenome) {
    grv_proc_regions = GenomicRegion::non_centromeres;
    if (opt::verbose > 0)
      cout << "Processing whole genome (minus centromeres)" << endl;
  } 
  // run everything including centromeres
  else if (runWholeGenome){
    //grv_proc_regions = GenomicRegion::getWholeGenome();
    //if (opt::verbose > 0)
    cout << "Processing whole genome (including centromeres)" << endl;
  }
  // if there is no proc_region file, take the mini collection as a region set
  else if (grv_proc_regions.size() == 0 && !runWholeGenome) {
    grv_proc_regions = mr->sendToGrv();
  }
  // check that there are some regions to run
  else  {
    cerr << "No regions to run. Did you specify a file as a region?" << endl;
    cerr << "Defaulting to whole genome" << endl;
    exit(EXIT_FAILURE);
  }

  // make sure the global mask is sorted
  if (grv_proc_regions.size() > 0)
    sort(grv_proc_regions.begin(), grv_proc_regions.end());      

  BamQC qc; 

  VariantBamReader sv_reader(opt::bam, opt::out, mr, opt::verbose);
  sv_reader.start = start;

  // dummy vector to store reads. Won't be used because writer is != NULL
  ReadVec bav;
 
  sv_reader.twopass = opt::twopass;

  if (runWholeGenome) {
    // first pass
    //sv_reader.writeVariantBam(qc, bav);
    sv_reader.writeVariantBam(qc, bav);
  }

  // loop through the process regions and extract files
  for (auto it = grv_proc_regions.begin(); it != grv_proc_regions.end(); it++) {

    sv_reader.setBamRegion(*it);

    if (opt::verbose > 0)
      cout << "Running region: "  << (*it) << endl;

    //sv_reader.writeVariantBam(qc, bav);

    sv_reader.writeVariantBam(qc, bav);

  }

  sv_reader.rewind();

  // second pass
  if (sv_reader.twopass && runWholeGenome) {

    sv_reader.writeVariantBamFromHash();

  } else if (sv_reader.twopass) {
    
    // loop through the process regions and extract files
    for (auto it = grv_proc_regions.begin(); it != grv_proc_regions.end(); it++) {
      sv_reader.setBamRegion(*it);
      if (opt::verbose > 0)
	cout << "Running region: "  << (*it) << endl;
      sv_reader.writeVariantBamFromHash();
    }

  }

  // clean up
  delete mr;

  // index it
#ifdef HAVE_BAMTOOLS
  sv_reader.MakeIndex();
#endif

  // write out the stats
  ofstream stats_out;
  stats_out.open("qcreport.txt");
  stats_out << qc;

  if (opt::verbose > 0) {
    VarUtils::displayRuntime(start);
    cout << endl;
  }
  
  return 0;
}

void parseVarOptions(int argc, char** argv) {

  bool die = false;

  if (argc < 2) 
    die = true;

  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'h': die = true; break;
    case 'v': arg >> opt::verbose; break;
    case 't': opt::twopass = true; break;
    case 'i': arg >> opt::bam; break;
    case 'o': arg >> opt::out; break;
    case 'w': arg >> opt::mapq; break;
    case 'g': 
      {
	string tmp;
	arg >> tmp;
	if (tmp.length() == 0)
	  break;
	if (opt::rules_file.length())
	  opt::rules_file += "%";
	opt::rules_file += "region@" + tmp;
      }
      break;
    case 'c':  // call stats hack
      {
	string tmp;
	arg >> tmp;
	if (tmp.length() == 0)
	  break;
	if (opt::rules_file.length())
	  opt::rules_file += "%";
	opt::rules_file += "region@" + tmp + ";mate";
      }
      break;

    case 'r': 
      {
	string tmp;
	arg >> tmp;
	if (tmp.length() == 0)
	  break;
	if (opt::rules_file.length())
	  opt::rules_file += "%";
	else 
	  opt::rules_file = "region@WG%"; // need to specify a region
	opt::rules_file += tmp;
      }
      break;
    case 'k': arg >> opt::proc_regions; break;
    case 'f': 
      {
	string file;
	arg >> file;
	
	string line;
	igzstream iss(file.c_str());
	while(getline(iss,line, '\n')) {
	  if (opt::rules_file.length() && line.length())
	    opt::rules_file += "%";
	  if (line.length())
	    opt::rules_file += line;
	}
      }
      break;
      //case 'm': arg >> opt::min_length; break;
      //case 'z': arg >> opt::min_phred; break;
      //case 'c': arg >> opt::min_clip; break;
      //case 'e': opt::exclude_n = true; break;
    case 'b': opt::skip_cent = true; break;
    }
  }

  if (opt::bam == "")
    die = true;
  if (opt::out == "")
    die = true;

  // dont stop the run for bad bams for quality checking only
  //opt::perc_limit = opt::qc_only ? 101 : opt::perc_limit;

  // something went wrong, kill
  if (die) {
    cout << "\n" << VARIANT_BAM_USAGE_MESSAGE;
    exit(1);
  }

}


