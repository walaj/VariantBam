#include <string>
#include <getopt.h>
#include <iostream>

#include "SnowTools/gzstream.h"
#include "SnowTools/SnowUtils.h"
#include "SnowTools/GenomicRegionCollection.h"
#include "SnowTools/GenomicRegion.h"
#include "SnowTools/SnowToolsCommon.h"

#include "VariantBamWalker.h"

using SnowTools::GenomicRegion;
using SnowTools::GenomicRegionCollection;
using SnowTools::GRC;

static const char *VARIANT_BAM_USAGE_MESSAGE =
"Usage: variant <input.bam> -g <regions> -r <rules> [OPTIONS] \n\n"
"  Description: Filter a BAM/CRAM file according to hierarchical rules\n"
"\n"
" General options\n"
"      --help                           Display this help and exit\n"
"  -v, --verbose                        Verbose output\n"
"  -c, --counts-file                    File to place read counts per rule / region\n"
"  -x, --counts-file-only               Same as -c, but does counting only (no output BAM)\n"
" Output options\n"
"  -o, --output-bam                     Output BAM file to write instead of SAM-format stdout\n"
"  -C, --cram                           Output file should be in CRAM format\n"
"  -T, --reference                      Path to reference. Required for reading/writing CRAM\n"
"  -h, --include-header                 When outputting to stdout, include the header.\n"
"  -s, --strip-tags                     Remove the specified tags, separated by commas. eg. -s RG,MD\n"
"  -S, --strip-all-tags                 Remove all alignment tags\n"
" Filtering options\n"
"  -q, --qc-only                        Loop through the BAM, but only to make the QC file\n"
"  -g, --region                         Regions (e.g. myvcf.vcf or WG for whole genome) or newline seperated subsequence file.  Applied in same order as -r for multiple\n"
"  -l, --linked-region                  Same as -g, but turns on mate-linking\n"
"  -r, --rules                          Script for the rules. If specified multiple times, will be applied in same order as -g\n"
"  -k, --proc-regions-file              BED file of regions to proess reads from\n"
"\n";

namespace opt {

  static std::string bam;
  static std::string out;
  static bool verbose = false;
  static std::string rules = "";
  static std::string proc_regions = "";
  static bool to_stdout = false;
  static bool cram = false;
  static bool header = false;
  static std::string reference = SnowTools::REFHG19;
  static bool strip_all_tags = false;
  static std::string tag_list = "";
  static std::string counts_file = "";
  static bool counts_only = false;
}

enum {
  OPT_HELP
};

static const char* shortopts = "hvqji:o:r:k:g:Cf:s:ST:l:c:x:";
static const struct option longopts[] = {
  { "help",                       no_argument, NULL, OPT_HELP },
  { "linked-region",              required_argument, NULL, 'l' },
  { "counts-file",                required_argument, NULL, 'c' },
  { "counts-file-only",           required_argument, NULL, 'x' },
  { "cram",                       no_argument, NULL, 'C' },
  { "strip-all-tags",             no_argument, NULL, 'S' },
  { "strip-tags",                 required_argument, NULL, 's' },
  { "reference",                  required_argument, NULL, 'T' },
  { "verbose",                    no_argument, NULL, 'v' },
  { "include-header",             no_argument, NULL, 'h' },
  { "input",                      required_argument, NULL, 'i' },
  { "output-bam",                 required_argument, NULL, 'o' },
  { "qc-only",                    no_argument, NULL, 'q' },
  { "rules",                      required_argument, NULL, 'r' },
  { "region",                     required_argument, NULL, 'g' },
  { "region-with-mates",          required_argument, NULL, 'c' },
  { "proc-regions-file",          required_argument, NULL, 'k' },
  { "no-pileup-check",            no_argument, NULL, 'j' },
  { NULL, 0, NULL, 0 }
};

static struct timespec start;


// forward declare
void parseVarOptions(int argc, char** argv);

int main(int argc, char** argv) {

#ifndef __APPLE__
  // start the timer
  clock_gettime(CLOCK_MONOTONIC, &start);
#endif

  // parse the command line
  parseVarOptions(argc, argv);

  if (opt::verbose) {
    //std::cerr << "Input BAM:  " << opt::bam << std::endl;
    //std::cerr << "Output BAM: " << opt::out << std::endl;
    //std::cerr << "Input rules and regions: " << opt::rules << std::endl;
    //std::cerr << "Input proc regions file: " << opt::proc_regions << std::endl;
    //std::cerr << "TWO-PASS solution?:      " << (opt::twopass ? "ON" : "OFF") << std::endl;
  }

  // setup the walker
  VariantBamWalker walk(opt::bam);

  // set which regions to run
  GRC grv_proc_regions;
  if (opt::proc_regions.length())
    grv_proc_regions.regionFileToGRV(opt::proc_regions, 0, walk.header()); // 0 is pad
  grv_proc_regions.createTreeMap();
  
  // should it print to stdout?
  if (opt::to_stdout) {
    walk.setStdout();
    if (opt::header)
      walk.setPrintHeader();
  }
  // should we print to cram
  else if (opt::cram) {
    walk.setCram(opt::out, opt::reference);
  }

  // should we clear tags?
  if (opt::strip_all_tags)
    walk.setStripAllTags();
  else if (opt::tag_list.length())
    walk.setStripTags(opt::tag_list);

  // make the mini rules collection from the rules file
  // this also calls function to parse the BED files
  if (opt::verbose)
    std::cerr << "...rules: " << opt::rules << std::endl;
  walk.SetMiniRulesCollection(opt::rules);

  // set the regions to run
  if (grv_proc_regions.size())
    walk.setBamWalkerRegions(grv_proc_regions.asGenomicRegionVector());

  /*SnowTools::GRC rules_rg = walk.GetMiniRulesCollection().getAllRegions();
  rules_rg.createTreeMap();

  if (grv_proc_regions.size() && rules_rg.size()) // intersect rules regions with mask regions
    rules_rg = rules_rg.intersection(grv_proc_regions, true); // true -> ignore_strand
  else if (grv_proc_regions.size())
    rules_rg = grv_proc_regions; // rules is whole genome, so just make mask instead

  if (rules_rg.size()) {
    walk.setBamWalkerRegions(rules_rg.asGenomicRegionVector());
  } else if (!rules_rg.size() && grv_proc_regions.size() > 0) {
    std::cerr << "No regions with possibility of reads. This error occurs if no regions in -g are in -k." << std::endl;
    return 1;
    }*/

  // should we count all rules (slower)
  if (opt::counts_only || opt::counts_file.length())
    walk.setCountAllRules();

  // open the output file
  if (!opt::counts_only)
    walk.OpenWriteBam(opt::out);

  // print out some info
  if (opt::verbose) 
    std::cerr << walk << std::endl;

  // set verbosity of walker
  if (opt::verbose)
    walk.setVerbose();

  // do the filtering
  if (opt::verbose)
    std::cerr << "...starting filtering" << std::endl;
  walk.writeVariantBam();

  // display the rule counts
  walk.MiniRulesToFile(opt::counts_file);

  // make a bed file
  //if (opt::verbose > 0)
  //  std::cerr << "...sending merged regions to BED file" << std::endl;
  //mr->sendToBed("merged_rules.bed");

  // index it
#ifdef HAVE_BAMTOOLS
  walk.MakeIndex();
#endif

  if (opt::verbose) {
    SnowTools::displayRuntime(start);
    std::cerr << std::endl;
  }
  
  return 0;
}

void parseVarOptions(int argc, char** argv) {

  bool die = false;

  if (argc < 2) {
    std::cerr << "\n" << VARIANT_BAM_USAGE_MESSAGE;
    exit(EXIT_FAILURE);
  }

  opt::bam = std::string(argv[1]);

  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case OPT_HELP: die = true; break;
    case 'v': opt::verbose = true; break;
    case 's': arg >> opt::tag_list; break;
    case 'S': opt::strip_all_tags = true; break;
    case 'T': arg >> opt::reference; break;
    case 'C': opt::cram = true; break;
    case 'h': opt::header = true;
      //case 't': opt::twopass = true; break;
    case 'i': arg >> opt::bam; break;
    case 'o': arg >> opt::out; break;
    case 'l': 
      {
	std::string tmp;
	arg >> tmp;
	if (tmp.length() == 0)
	  break;
	if (opt::rules.length())
	  opt::rules += "%";
	opt::rules += "mlregion@" + tmp;
      }
      break;
    case 'g': 
      {
	std::string tmp;
	arg >> tmp;
	if (tmp.length() == 0)
	  break;
	if (opt::rules.length())
	  opt::rules += "%";
	opt::rules += "region@" + tmp;
      }
      break;
    case 'c': arg >> opt::counts_file; break;
    case 'x': arg >> opt::counts_file; opt::counts_only = true; break;
    case 'r': 
      {
	std::string tmp;
	arg >> tmp;
	if (tmp.length() == 0)
	  break;

	// check if it's a file
	if (SnowTools::read_access_test(tmp)) 
	  {
	    std::ifstream iss(tmp);
	    std::string val;

	    while(std::getline(iss, val))
	      {
		if (val.find("#") == std::string::npos)
		  {
		    opt::rules += val;
		    opt::rules += "%";
		  }
	      }
	    //trim off the last %
	    opt::rules = SnowTools::cutLastChar(opt::rules); 
	  }
	else if (opt::rules.length()) {
	  opt::rules += "%"; // adding a new line
	  opt::rules += tmp;
	}
	else { // region is empty, so add one first
	  opt::rules = "region@WG%" + tmp; // need to specify a region
	}
	//opt::rules += tmp;
      }
      break;
    case 'k': arg >> opt::proc_regions; break;
      /*    case 'f': 
      {
	string file;
	arg >> file;
	
	string line;
	igzstream iss(file.c_str());
	while(getline(iss, line, '\n')) {
	  if (opt::rules.length() && line.length())
	    opt::rules += "%";
	  if (line.length())
	    opt::rules += line;
	}
      }
      break;*/
    }
  }

  if (opt::bam == "")
    die = true;
  if (opt::out == "")
    opt::to_stdout = true;

  // dont stop the run for bad bams for quality checking only
  //opt::perc_limit = opt::qc_only ? 101 : opt::perc_limit;

  // something went wrong, kill
  if (die) {
    std::cerr << "\n" << VARIANT_BAM_USAGE_MESSAGE;
    exit(1);
  }

}

