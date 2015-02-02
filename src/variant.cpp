#include <string>
#include <vector>
#include <getopt.h>
#include <iostream>
#include "GenomicRegion.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
//#include "BamQC.h"
#include "VarUtils.h"
#include "MiniRules.h"
#include "VariantBamReader.h"

using namespace std;
using namespace BamTools;

static unordered_map<string, bool> valid;

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
"  -r, --rules-file                     Rules file\n"
"  -k, --proc-regions-file              csv file of regions to proess reads from, of format chr,pos1,pos2 eg 11,2232423,2235000. Default: Whole genome\n"
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
  

  static int mapq = 0;
  static bool skip_cent = false;
  static string rules_file = "";
  static string proc_regions = "";

  //static int perc_limit = 30;

  //static bool no_pileup_check = false;
  //static bool qc_only = false;
}

static const char* shortopts = "hv:qejb:i:o:w:n:p:s:r:m:z:c:k:u:x:";
static const struct option longopts[] = {
  { "help",                       no_argument, NULL, 'h' },
  { "verbose",                    required_argument, NULL, 'v' },
  { "input-bam",                  required_argument, NULL, 'i' },
  { "output-bam",                 required_argument, NULL, 'o' },
  //{ "min-mapq",                 required_argument, NULL, 'w' },
  //{ "min-clip",                 required_argument, NULL, 'c' },
  //{ "nm-limit",                 required_argument, NULL, 'n' },
  //{ "perc-limit",                 required_argument, NULL, 'p' },
  //{ "isize",                 required_argument, NULL, 's' },
  { "qc-only",                 no_argument, NULL, 'q' },
  { "rules-file",                 required_argument, NULL, 'r' },
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

  // define what is a valid condition
  valid["duplicate"]  = true;
  valid["supplementary"]       = true;
  valid["qcfail"]     = true;
  valid["hardclip"]   = true;
  valid["fwd_strand"] = true;
  valid["rev_strand"] = true;
  valid["mate_fwd_strand"] = true;
  valid["mate_rev_strand"] = true;
  valid["mapped"]          = true;
  valid["mate_mapped"]     = true;
  valid["isize"] = true;
  valid["clip"] = true;
  valid["phred"] = true;
  valid["len"] = true;
  valid["nm"] = true;
  valid["mapq"] = true;

  // start the timer
  clock_gettime(CLOCK_MONOTONIC, &start);

  parseVarOptions(argc, argv);

  if (opt::verbose > 0) {
    cout << "Input BAM:  " << opt::bam << endl;
    cout << "Output BAM: " << opt::out << endl;
    cout << "Input rules script (file or script): " << opt::rules_file << endl;
    cout << "Input proc regions file: " << opt::proc_regions << endl;
  }

  // make the mini rules collection from the rules file
  // this also calls function to parse the BED files
  MiniRulesCollection * mr = new MiniRulesCollection(opt::rules_file);
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
  if (opt::proc_regions.size() > 0)
    GenomicRegion::regionFileToGRV(opt::proc_regions, 0);
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
    grv_proc_regions = GenomicRegion::getWholeGenome();
    if (opt::verbose > 0)
      cout << "Processing whole genome (including centromeres)" << endl;
  }
  // if there is no proc_region file, take the mini collection as a region set
  else if (grv_proc_regions.size() == 0 && !runWholeGenome) {
    grv_proc_regions = mr->sendToGrv();
  }

  // check that there are some regions to run
  if (grv_proc_regions.size() == 0) {
    cerr << "No regions to run. Did you specify a file as a region" << endl;
    exit(EXIT_FAILURE);
  }


  // make sure the global mask is sorted
  if (grv_proc_regions.size() > 0)
    sort(grv_proc_regions.begin(), grv_proc_regions.end());      

  BamAlignment a;
  //  BamQC qc; 

  VariantBamReader sv_reader(opt::bam, opt::out, mr, opt::verbose);

  // loop through the process regions and extract files
  for (auto it = grv_proc_regions.begin(); it != grv_proc_regions.end(); it++) {

    sv_reader.setBamRegion(*it);

    if (opt::verbose > 0)
      cout << "Running region: "  << (*it) << endl;
    //sv_reader.writeVariantBam(qc, opt::qc_only);
    sv_reader.writeVariantBam();
    
    if (opt::verbose > 0) {
      VarUtils::displayRuntime(start);
      cout << endl;
    }

  }

  // index it
  sv_reader.MakeIndex();

  // write out the stats
  //ofstream stats_out;
  //stats_out.open("./qcreport.txt");
  //stats_out << qc;

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
    case 'i': arg >> opt::bam; break;
    case 'o': arg >> opt::out; break;
    case 'w': arg >> opt::mapq; break;
      //case 'n': arg >> opt::nmlim; break;
      //case 'u': arg >> opt::mutect_callstats; break;
      //case 'x': arg >> opt::mutect2_regions; break;
      //case 'p': arg >> opt::perc_limit; break;
      //case 's': arg >> opt::isize; break;
      //case 'q': opt::qc_only = true; break;
    case 'r': 
      if (opt::rules_file.length() > 0) {
	string tmp;
        arg >> tmp;
	if (tmp.length() > 0)
	  opt::rules_file += "%" + tmp;
	break;
      } else {
	arg >> opt::rules_file; 
	break;
      }
    case 'k': arg >> opt::proc_regions; break;
      //case 'm': arg >> opt::min_length; break;
      //case 'z': arg >> opt::min_phred; break;
      //case 'c': arg >> opt::min_clip; break;
      //case 'e': opt::exclude_n = true; break;
    case 'b': opt::skip_cent = true; break;
    }
  }

  // dont stop the run for bad bams for quality checking only
  //opt::perc_limit = opt::qc_only ? 101 : opt::perc_limit;

  // something went wrong, kill
  if (die) {
    cout << "\n" << VARIANT_BAM_USAGE_MESSAGE;
    exit(1);
  }

}


