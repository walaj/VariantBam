/**
Copyright 2016 Jeremiah A. Wala (jwala@broadinstitute.org)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

#include <string>
#include <getopt.h>
#include <iostream>
#include <fstream>

#include "SeqLib/SeqLibUtils.h"
#include "SeqLib/GenomicRegionCollection.h"
#include "SeqLib/GenomicRegion.h"
#include "SeqLib/SeqLibCommon.h"

#include "VariantBamWalker.h"
#include "CommandLineRegion.h"

using SeqLib::GenomicRegion;
using SeqLib::GenomicRegionCollection;
using SeqLib::GRC;

static const char *VARIANT_BAM_USAGE_MESSAGE =
"Usage: variant <input.bam> [OPTIONS] \n\n"
"  Description: Filter a BAM/SAM/CRAM/STDIN according to hierarchical rules\n"
"\n"
" General options\n"
"  -h, --help                           Display this help and exit\n"
"  -v, --verbose                        Verbose output\n"
  //"  -c, --counts-file                    File to place read counts per rule / region\n"
"  -x, --no-output                      Don't output reads (used for profiling with -q)\n"
"  -r, --rules                          JSON ecript for the rules.\n"
"  -k, --proc-regions-file              Samtools-style region string (e.g. 1:1,000-2,000) or BED/VCF of regions to process. -k UN iterates over unmapped-unmapped reads\n"
"  -Q, --mark-as-qc-fail                Flag reads that don't pass VariantBam with the failed QC flag, rather than deleting the read.\n"
" Output options\n"
"  -o, --output                         Output file to write to (BAM/SAM/CRAM) file instead of stdout\n"
"  -C, --cram                           Output file should be in CRAM format\n"
"  -b, --bam                            Output should be in binary BAM format\n"
"  -T, --reference                      Path to reference. Required for reading/writing CRAM\n"
"  -s, --strip-tags                     Remove the specified tags, separated by commas. eg. -s RG,MD\n"
"  -S, --strip-all-tags                 Remove all alignment tags\n"
"  -Z, --write-trimmed                  Output the base-quality trimmed sequence rather than the original sequence. Also removes quality scores\n"
" Filtering options\n"
"  -q, --qc-file                        Output a qc file that contains information about BAM\n"
"  -m, --max-coverage                   Maximum coverage of output file. BAM must be sorted. Negative values enforce a minimum coverage\n"
"  -p, --min-phred                      Set the minimum base quality score considered to be high-quality\n"
" Region specifiers\n"
"  -g, --region                         Regions (e.g. myvcf.vcf or WG for whole genome) or newline seperated subsequence file.\n"
"  -G, --exclude-region                 Same as -g, but for region where satisfying a rule EXCLUDES this read.\n"
"  -l, --linked-region                  Same as -g, but turns on mate-linking\n"
"  -L, --linked-exclude-region          Same as -l, but for mate-linked region where satisfying this rule EXCLUDES this read.\n"
"  -P, --region-pad                     Apply a padding to each region supplied with the region flags (specify after region flag)\n"
" Command line rules shortcuts (to be used without supplying a -r script)\n"
"      --min-clip                       Minimum number of quality clipped bases\n"
"      --max-nbases                     Maximum number of N bases\n"
"      --min-mapq                       Minimum mapping quality\n"
"      --min-del                        Minimum number of deleted bases\n"
"      --min-ins                        Minimum number of inserted bases\n"
"      --min-length                     Minimum read length (after base-quality trimming)\n"
"      --motif                          Motif file\n"
"  -R, --read-group                     Limit to just a single read group\n"
"  -f, --include-aln-flag               Flags to include (like samtools -f)\n"
"  -F, --exclude-aln-flag               Flags to exclude (like samtools -F)\n"
"\n";

std::vector<CommandLineRegion> command_line_regions;

void __check_command_line(std::vector<CommandLineRegion>& c) {

  if (!c.size())
    c.push_back(CommandLineRegion("WG", -1)); // add whole-genome ALL rule

}

namespace opt {

  static int phred = -1;
  static std::string blacklist;
  static std::string bam;
  static std::string out;
  static int max_cov = 0;
  static bool verbose = false;
  static std::string rules;
  static std::string proc_regions;
  static bool cram = false;
  static std::string reference;
  static bool strip_all_tags = false;
  static std::string tag_list;
  static std::string counts_file;
  static bool noop = false;
  static std::string bam_qcfile;
  static bool bam_output = false;
  static bool write_trimmed = false; // write the quality trimmed read?

  static bool mark_as_qcfail = false; // mark failed reads with QC fail flag, instead of deleting
}

enum {
  OPT_HELP,
  OPT_LENGTH,
  OPT_MAPQ,
  OPT_NBASES,
  OPT_CLIP,
  OPT_MOTIF,
  OPT_INS, 
  OPT_DEL
};

static const char* shortopts = "hvbxi:o:r:k:g:Cf:s:ST:l:c:q:m:L:G:P:F:R:p:QZ";
static const struct option longopts[] = {
  { "help",                       no_argument, NULL, 'h' },
  { "bam",                        no_argument, NULL, 'b' },
  { "linked-region",              required_argument, NULL, 'l' },
  { "write-trimmed",              no_argument, NULL, 'Z'} ,
  { "mark-as-qc-fail",              no_argument, NULL, 'Q'} ,
  { "min-length",              required_argument, NULL, OPT_LENGTH },
  { "motif",              required_argument, NULL, OPT_MOTIF },
  { "min-ins",              required_argument, NULL, OPT_INS},
  { "min-del",              required_argument, NULL, OPT_DEL },
  { "min-phred",              required_argument, NULL, 'p' },
  { "min-mapq",              required_argument, NULL, OPT_MAPQ },
  { "min-clip",              required_argument, NULL, OPT_CLIP },
  { "max-nbases",              required_argument, NULL, OPT_NBASES },
  { "read-group",              required_argument, NULL, 'R' },
  { "include-aln-flag",             required_argument, NULL, 'f' },
  { "exclude-aln-flag",             required_argument, NULL, 'F' },
  { "exclude-region",             required_argument, NULL, 'G' },
  { "linked-exclude-region",      required_argument, NULL, 'L' },
  { "max-coverage",               required_argument, NULL, 'm' },
  //  { "counts-file",                required_argument, NULL, 'c' },
  { "no-output",                  required_argument, NULL, 'x' },
  { "cram",                       no_argument, NULL, 'C' },
  { "strip-all-tags",             no_argument, NULL, 'S' },
  { "strip-tags",                 required_argument, NULL, 's' },
  { "reference",                  required_argument, NULL, 'T' },
  { "verbose",                    no_argument, NULL, 'v' },
  { "input",                      required_argument, NULL, 'i' },
  { "output",                 required_argument, NULL, 'o' },
  { "qc-file",                    no_argument, NULL, 'q' },
  { "rules",                      required_argument, NULL, 'r' },
  { "region",                     required_argument, NULL, 'g' },
  { "region-pad",                 required_argument, NULL, 'P' },
  //  { "region-with-mates",          required_argument, NULL, 'c' },
  { "proc-regions-file",          required_argument, NULL, 'k' },
  { NULL, 0, NULL, 0 }
};

std::string myreplace(std::string &s,
                      std::string toReplace,
                      std::string replaceWith)
{
  if (s.find(toReplace) == std::string::npos)
    return (s);
  return(s.replace(s.find(toReplace), toReplace.length(), replaceWith));
}

// forward declare
void parseVarOptions(int argc, char** argv);

// helper for formatting rules script string with no whitespace
// http://stackoverflow.com/questions/83439/remove-spaces-from-stdstring-in-c
/*  template<typename T, typename P>
    T remove_if(T beg, T end, P pred)
  {
    T dest = beg;
    for (T itr = beg;itr != end; ++itr)
      if (!pred(*itr))
	*(dest++) = *itr;
    return dest;
  }
*/

int main(int argc, char** argv) {

#ifndef __APPLE__
  // start the timer
  // clock_gettime(CLOCK_MONOTONIC, &start);
#endif

  // parse the command line
  parseVarOptions(argc, argv);

  bool has_ml_region = opt::rules.find("mlregion") != std::string::npos;
  
  if (opt::verbose) {
    //std::cerr << "Input BAM:  " << opt::bam << std::endl;
    //std::cerr << "Output BAM: " << opt::out << std::endl;
    //std::cerr << "Input rules and regions: " << opt::rules << std::endl;
    //std::cerr << "Input proc regions file: " << opt::proc_regions << std::endl;
    //std::cerr << "TWO-PASS solution?:      " << (opt::twopass ? "ON" : "OFF") << std::endl;
  }

  if (has_ml_region && opt::verbose) {
    std::cerr << "...mate-linked region supplied. Defaulting to whole BAM run unless trimmed explicitly with -k flag" << std::endl;
  }

  // setup the walker
  VariantBamWalker reader;
  if (!reader.Open(opt::bam))  {
    std::cerr << "ERROR: could not open file " << opt::bam << std::endl;
    exit(EXIT_FAILURE);
  }
    
  // set whether to mark failed as QC fail, or just delete (default)
  reader.m_mark_qc_fail = opt::mark_as_qcfail;
  
  // set the phred trim limit
  reader.phred = opt::phred;

  GRC grv_proc_regions;
  if (opt::proc_regions.length()) {
    if (SeqLib::read_access_test(opt::proc_regions)) {
      grv_proc_regions = GRC(opt::proc_regions, reader.Header());
    } else if (opt::proc_regions.find(":") != std::string::npos) {
      grv_proc_regions.add(SeqLib::GenomicRegion(opt::proc_regions, reader.Header()));
    } else if (opt::proc_regions == "-1" || opt::proc_regions == "UN") {
      grv_proc_regions.add(SeqLib::GenomicRegion(-2, 0, 0));
    } else {
      std::cerr << "...unexpected region format or could not read file" << std::endl;
      exit(EXIT_FAILURE);
    }
    grv_proc_regions.CreateTreeMap();
  }

  // open for writing
  if (!opt::noop) {
    if (opt::out.empty()) {
      reader.m_writer = SeqLib::BamWriter(opt::bam_output ? SeqLib::BAM : SeqLib::SAM);
      reader.m_writer.SetHeader(reader.Header());
      reader.m_writer.Open("-");
    }
    // should we print to cram
    else if (opt::cram) {
      reader.m_writer = SeqLib::BamWriter(SeqLib::CRAM);
      reader.m_writer.SetHeader(reader.Header());

      if (!reader.m_writer.Open(opt::out)) {
	std::cerr << "ERROR: could not open output CRAM " << opt::out << std::endl;
	exit(EXIT_FAILURE);
      }
      if (!reader.m_writer.SetCramReference(opt::reference)) {
	std::cerr << "Failed to set CRAM reference file: " << opt::reference << std::endl;
	exit(EXIT_FAILURE);
      }
    } else {
      reader.m_writer = SeqLib::BamWriter(opt::bam_output ? SeqLib::BAM : SeqLib::SAM);
      reader.m_writer.SetHeader(reader.Header());
      if (!reader.m_writer.Open(opt::out)) {
	std::cerr << "ERROR: could not open output " << (opt::bam_output ? "BAM" : "SAM") << opt::out << std::endl;
	exit(EXIT_FAILURE);
      }
    }
    reader.m_writer.WriteHeader();
  }

  // should we clear tags?
  if (opt::strip_all_tags)
    reader.m_strip_all_tags = true; //.setStripAllTags();
  else if (opt::tag_list.length()) {
    std::istringstream iss(opt::tag_list);
    std::string val;
    while(std::getline(iss, val, ',')) {
      reader.m_tags_to_strip.push_back(val);
    }
  }

  // make the mini rules collection from the rules file
  // this also calls function to parse the BED files
  if (opt::verbose) {
    std::string str = opt::rules;
    str.erase(std::remove_if(str.begin(), str.end(), [](char x){return std::isspace(x);}),str.end());
    std::cerr << "Rules script: " << str << std::endl;
  }

  SeqLib::Filter::ReadFilterCollection rfc;
  
  if (!opt::rules.empty())
    rfc = SeqLib::Filter::ReadFilterCollection(opt::rules, reader.Header());

  // make sure command_line_reigons makes sense
  if (command_line_regions.size() == 2 && command_line_regions[1].all()) {
    std::cerr << "***************************************************" << std::endl
              << "  Region (-l, -L, -g, -G) supplied after rule flags"
              << "  Did you mean to set region flag before rule flags"
              << "***************************************************" << std::endl;
    exit(EXIT_FAILURE);
  }
    

  if (opt::verbose && command_line_regions.size())
    std::cerr << "...building rules from command line" << std::endl;

  // add specific mini rules from command-line
  for (auto& i : command_line_regions) {
    SeqLib::Filter::ReadFilter rf = BuildReadFilterFromCommandLineRegion(i, reader.Header());
    //SeqLib::MiniRules mr(i, walk.header());
    //SeqLib::ReadFilter rf(i, reader.Header());
    //mr.pad = i.pad;
    //mr.mrc = &mrc;
    rfc.AddReadFilter(rf);
    //mrc.m_regions.push_back(mr);
  }
  
  reader.m_mr = rfc;

  if (opt::verbose)
    std::cerr << rfc << std::endl;

  // set max coverage
  reader.max_cov = opt::max_cov;
  if (opt::max_cov > 0 && opt::verbose)
    std::cerr << "--- Setting MAX coverage to: " << opt::max_cov << std::endl;

  // set the regions to run
  if (grv_proc_regions.size()) {
    if (opt::verbose)
       std::cerr << "...from -g flag will run on " << grv_proc_regions.size() << " regions" << std::endl;
    //walk.setBamWalkerRegions(grv_proc_regions.asGenomicRegionVector());
    reader.SetMultipleRegions(grv_proc_regions);
  }

  SeqLib::GRC rules_rg = grv_proc_regions; //reader.GetMiniRulesCollection().getAllRegions();

  rules_rg.CreateTreeMap();

  if (grv_proc_regions.size() && rules_rg.size()) { // intersect rules regions with mask regions. 
    // dont incorporate rules regions if there are any mate-linked regions
    rules_rg = rules_rg.Intersection(grv_proc_regions, true); // true -> ignore_strand
    if (opt::verbose)
      std::cerr << "rules region " << rules_rg.size() << std::endl;
  } else if (grv_proc_regions.size()) {
    rules_rg = grv_proc_regions; // rules is whole genome, so just make mask instead
  }

  if (grv_proc_regions.size() > 0 && (rules_rg.size() || has_ml_region )) // explicitly gave regions
    //walk.setBamWalkerRegions(grv_proc_regions.asGenomicRegionVector());
    reader.SetMultipleRegions(grv_proc_regions);
  /*  else if (rules_rg.size() && !has_ml_region && grv_proc_regions.size() == 0) {
    walk.setBamWalkerRegions(rules_rg.asGenomicRegionVector());
    if (opt::verbose)
      std::cerr << "...from rules, will run on " << rules_rg.size() << " regions" << std::endl;
  } else if (!rules_rg.size() && grv_proc_regions.size() > 0) {
    std::cerr << "No regions with possibility of reads. This error occurs if no regions in -g are in -k." << std::endl;
    return 1;
    }*/

  // should we count all rules (slower)
  //if (opt::counts_file.length())
  //  reader.m_mr.setCountAllRules();

  // print out some info
  if (opt::verbose) 
    std::cerr << reader << std::endl;

  // set verbosity of walker
  reader.m_verbose = opt::verbose;

  // do the filtering
  if (opt::verbose)
    std::cerr << "...starting filtering" << std::endl;

  // set the trim writer opeion
  reader.m_write_trimmed = opt::write_trimmed;

  ////////////
  /// RUN THE WALKER
  ////////////
  reader.writeVariantBam();

  // dump the stats file
  if (!opt::bam_qcfile.empty()) {
    std::ofstream ofs;
    ofs.open(opt::bam_qcfile);
    ofs << reader.m_stats;
    ofs.close();
  }

  // display the rule counts
  /*
  std::ofstream cfile;
  cfile.open(opt::counts_file.c_str());
  cfile << reader.m_mr.EmitCounts(); //MiniRulesToFile(opt::counts_file);
  cfile.close();
  */

  // make a bed file
  //if (opt::verbose > 0)
  //  std::cerr << "...sending merged regions to BED file" << std::endl;
  //mr->sendToBed("merged_rules.bed");
  
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
    
    std::string tmp;
    switch (c) {
    case 'h': die = true; break;
    case 'v': opt::verbose = true; break;
    case 's': arg >> opt::tag_list; break;
    case 'S': opt::strip_all_tags = true; break;
    case 'T': arg >> opt::reference; break;
    case 'Z': opt::write_trimmed = true; break;
    case 'Q': opt::mark_as_qcfail = true; break;
    case 'C': opt::cram = true; break;
    case 'i': arg >> opt::bam; break;
    case 'o': arg >> opt::out; break;
    case 'm': arg >> opt::max_cov; break;
    case 'b': opt::bam_output = true; break;
    case 'l': 
	arg >> tmp;
	command_line_regions.push_back(CommandLineRegion(tmp, MINIRULES_MATE_LINKED));
	break;
    case 'L': 
	arg >> tmp;
	command_line_regions.push_back(CommandLineRegion(tmp, MINIRULES_MATE_LINKED_EXCLUDE));
	break;
    case 'g': 
	arg >> tmp;
	command_line_regions.push_back(CommandLineRegion(tmp, MINIRULES_REGION));
	break;
    case 'G': 
	arg >> tmp;
	command_line_regions.push_back(CommandLineRegion(tmp, MINIRULES_REGION_EXCLUDE));	
	break;
	//    case 'c': arg >> opt::counts_file; break;
    case 'R':
      __check_command_line(command_line_regions);
      arg >> command_line_regions.back().rg;
      break;
    case OPT_MAPQ:
      __check_command_line(command_line_regions);
      arg >> command_line_regions.back().mapq;
      break;
    case OPT_MOTIF:
      __check_command_line(command_line_regions);
      arg >> command_line_regions.back().motif;
      break;
    case OPT_CLIP:
      __check_command_line(command_line_regions);
      arg >> command_line_regions.back().clip;
      break;
    case OPT_LENGTH:
      __check_command_line(command_line_regions);
      arg >> command_line_regions.back().len;
      break;
    case 'p':
      arg >> opt::phred; break;
      break;
    case OPT_NBASES:
      __check_command_line(command_line_regions);
      arg >> command_line_regions.back().nbases;
      break;
    case OPT_INS:
      __check_command_line(command_line_regions);
      arg >> command_line_regions.back().ins;
      break;
    case OPT_DEL:
      __check_command_line(command_line_regions);
      arg >> command_line_regions.back().del;
      break;
    case 'x': opt::noop = true; break;
    case 'q': arg >> opt::bam_qcfile; break;
    case 'P': 
      if (!command_line_regions.size()) {
	std::cerr << "Error: Must input padding *after* specifying a region via -l, -L, -g, -G" << std::endl;
	exit(EXIT_FAILURE);
      }
      arg >> command_line_regions.back().pad;
      break;
    case 'f': 
      __check_command_line(command_line_regions);
      arg >> command_line_regions.back().i_flag;
      break;
    case 'F': 
      __check_command_line(command_line_regions);
      arg >> command_line_regions.back().e_flag;
      break;
    case 'r': 
	arg >> tmp;

	// check if it's a file
	if (SeqLib::read_access_test(tmp)) 
	  {
	    std::ifstream iss(tmp);
	    std::string val;

	    while(std::getline(iss, val))
	       opt::rules += val;
	  }
	else {
	  opt::rules = tmp;
	}
      break;
    case 'k': arg >> opt::proc_regions; break;
    }
  }

  if (opt::bam == "")
    die = true;

  // dont stop the run for bad bams for quality checking only
  //opt::perc_limit = opt::qc_only ? 101 : opt::perc_limit;

  // something went wrong, kill
  if (die) {
    std::cerr << "\n" << VARIANT_BAM_USAGE_MESSAGE;
    exit(1);
  }

}

