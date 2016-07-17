[![Build Status](https://travis-ci.org/jwalabroad/VariantBam.svg?branch=master)](https://travis-ci.org/jwalabroad/VariantBam)

VariantBam: Filtering and profiling of next-generational sequencing data using region-specific rules
======================================================================================================

<div style="text-align:center"><img src="https://raw.githubusercontent.com/jwalabroad/VariantBam/master/labels_vb.png" width="200"></div>

**License:** [GNU GPLv3][license]

[Bioinformatics Paper][biop]
Wala, J., C. Zhang, M. Meyerson, R. Beroukhim. VariantBam: filtering and profiling of nextgenerational sequencing data using region-specific rules. 2016. Bioinformatics, doi: 10.1093/bioinformatics/btw111 

**NOTE:** VariantBam recently was updated to use the more universal JSON syntax, and to remove all dependencies on Boost to make installation easier.


Table of contents
=================

  * [Installation](#gh-md-toc)
  * [Description](#description)
  * [Examples](#examples)
  * [Rules script syntax](#rules-script-syntax)
    * [Region](#region)
    * [Global region](#global-region)
    * [Rules](#rules)
      * [Range rules](#range-rules)
      * [Flag rules](#flag-rules)
      * [Motif rules](#motif-rules)
      * [Tag rules](#tag-rules)
      * [CIGAR rules](#cigar-rules)
      * [Subsample rules](#subsample-rules)            
  * [Command line usage](#command-line-usage)
    * [Full list of options](#full-list-of-options)
  * [Full list of available JSON rules](#full-list-of-available-json-rules)
  * [Attributions](#attributions)

Installation
============
I have succesfully built on Unix with GCC-4.8, 4.9 and 5.1

```
git clone --recursive https://github.com/jwalabroad/VariantBam.git
cd VariantBam
./configure
make
```

Quick Start
===========
```
## using the included test BAM (HCC1143)
VariantBam/src/variant test/small.bam -g 'X:1,000,000-1,100,000' --min-mapq 10 -c counts.tsv -o mini.bam -v

## get help
VariantBam/src/variant --help

## TL;DR examples

## extract all reads and their pair-mates that overlap SNP sites within 100 bp
rfile=<BED file, samtools-style string (e.g. "1:1,000,000-1,000,010"), or VCF>
variant <bam> -l $rfile -P 100 -o mini.bam -v

## mask regions (exclude reads and their pair-mates that overlap)
variant <bam> -L $rfile -o mini.bam -v 

## extract high-quality clipped reads (where clip length account for low quality bases)
variant <bam> --min-phred 4 --min-clip 5 -o mini.bam -v

## extract reads with high mapq that also contain a large insertion or deletion
variant <bam> --min-mapq 20 --min-ins 10 --min-del 10 -v -o mini.bam

## subsample to max-coverage. BAM must be sorted
variant <bam> -m 100 -o mini.bam -v
```

Description
===========

VariantBam is a tool to extract/count specific sets of sequencing reads from next-generational sequencing files. To save money, 
disk space and I/O, one may not want to store an entire BAM on disk. In many cases, it would be more efficient to store only those read-pairs or
reads who intersect some region around the variant locations. Alternatively, if your scientific question is focused on only one aspect of the data (e.g. breakpoints), many 
reads can be removed without losing the information relevant to the problem, and enriching for the signal you are interested in.

##### Tool comparison

VariantBam packages into a single executable a number of filtering features not easily found using ``samtools`` + ``awk``:

> 1. Filter specifically on read clipping, orientation and insert size (all important for structural variation)
> 2. Support for considering only high-quality bases when determining read length or clip count
> 3. [Interval tree][ekg] to efficiently determine if a read overlaps a region
> 4. Ability to link reads to a genomic region if their mate intersects that region.
> 5. Provide different rules for different arbitrarily-sized regions, and to provide these regions as common variant files (VCF, MAF, BED)
> 6. Select reads by matching motifs against a large dictionary using [Aho-Corasick implementation][aho]
> 7. Count reads that satisfy any number of user-defined properties
> 8. Selectively strip alignment tags
> 9. Support for sub-sampling to obtain a BAM file with a coverage limit

VariantBam is implemented in C++ and uses [HTSlib][hlib], a highly optimized C library used as the core of [Samtools][samtools] and [BCFtools][bcf].

To get a full list of options, run ``variant --help``.

Examples
========

##### Example Use 1
Whole-genome analysis has been conducted on a BAM, generating VCF and MAF files. Ideally, these regions could be manually inspected
or reanalyzed without having to keep the entire BAM. Running VariantBam to extract only read-pairs that overlap these events will allow
these regions to be rapidly queried, without having to keep the full BAM record.
```
### Extract all read PAIRS that interset with a variant from a VCF
variant $bam -l myvcf.vcf -o mini.bam
```

##### Example Use 2
In situations where the sequencing or library preparation quality is low, it may be advantageous
to remove poor quality reads before starting the analysis train. VariantBam handles this by optionally taking into
account base-qualities when making a decision whether to keep a sequencing read. For instance, one might 
only be interested in high quality MAPQ 0 or clipped reads. VariantBam can be 
setup to apply unique base-quality filters to different regions or across the entire genome, all with one-pass. 
```
### Extract only high quality reads with >= 50 bases of phred >=4 and MAPQ >= 1 and not duplicated/hardclip/qcfail
### json
{
  "region1" : {
  "region"  : "WG",
  "rules" : [{
      "phred" : [4, 100],
      "length" : [50,1000],
      "mapq" : [1, 60],
      "duplicate" : false,
      "hardclip" : false,
      "qcfail" : false
   }
  ]
   }
}
###
variant $bam -r example2.json -o mini.bam
```
##### Example Use 3
An NGS tool operates only on a subset of the reads (eg. structural variant caller using only clipped/discordant reads). Running VariantBam
to keep only these reads allows the tool to run much faster. This is particurlaly useful for facilitating a more rapid "build/test" cycle.
```
### Extract clipped, discordant, unmapped and indel reads
### json
{
  "global" : { "nbases" : [0,0], "hardclip" : false, "supplementary" : false, "qcfail" : false, "phred" : [4,100] },
  "region_wg" : {"region" : "WG", "rules" : [ 
         { "mapq" : [0, 1000], "clip" : [5,1000] }, {"ic" : true}, {"ff" : true}, {"rf" : true}, {"rr" : true}, { "ins" : [1,1000], "mapq" : [1,100] }, { "del" : [1,1000], "mapq" : [1,1000] } 
  ]}
}
###
variant $bam -r example3.json
```
##### Example Use 4
A user wants to profile a BAM for quality. They would like to count the number of clipped reads in a BAM file, so long
as those reads have sufficient optical quality and mapping quality. VariantBam run with the -x flag for "counting only" 
will accomplish this. Let's try an example of this, just for chromsome 22
```
## example4.json
{
"example4": {
  "region" : "22",
  "rules": [{"clip": [5,1000],
    	     "phred": [4, 1000],
             "length": [20, 1000]}]
}
}									  
##

### 
variant $bam -g 22 --min-clip 5 --min-phred 4 --min-mapq 10 -c counts.tsv
variant $bam -r example4.json -c counts.tsv ## using JSON
```
##### Example Use 5
A team is only interested in variants in known cancer genes, and would like to analyze thousands of exomes and genomes. Running 
VariantBam to extract reads from only these genes, and sending the BAM files to compressed CRAM provides sufficient data reduction
to allow all of the relevant data to be stored on disk.
```
### Grab only reads from predefined regions. Strip unneccessary tags and convert to CRAM for maximum compression
variant $bam -l mygenes.bed -C -o mini.cram -s BI,OQ
```
##### Example Use 6
A research team would like to extract only reads matching a certain motifs, but only if they have high optical quality. 
VariantBam with the ``motif`` rule will accomplish this with rapid O(n) efficiency for an arbitrarily large motif dictionary (where ``n`` is
the length of a read)
```
## example6.json
{
"example6": {
  "rules": [{"motif": "mymotifs.txt",
    	     "phred": 4,
             "length": 20 }]
}
}									  
##

### 
variant $bam -r example6.json ## input as a JSON
variant $bam --min-phred 4 --min-length 20 --motif mymotifs.txt ## using command line shortcuts
```

##### Example Use 7
To reduce the size of the BAM, reads can be removed from centromeric and satellite repeat regions. These reads are rarely helpful for variant calling.
To remove reads that intersect a region, set the region as an inverse-region. In a VariantBam script, use ``"exclude" : true```. For 
quick use on the command line, use ``-L`` or ``-G`` (opposites of ``-l`` and ``-g``).
```
### json 
{
  "" : {
         "region" : "bad.bed",
         "exclude" : true,
	 "matelink" : true
       }
}
###
variant $bam -L bad.bed -o mini.bam -v
```

##### Example Use 8
Massive read-pileups can occur at repetitive regions. These can reduced with VariantBam by subsampling to a max-coverage.
```
### BAM must be sorted
variant $bam -m 100 -o mini.bam -v
```

##### Example Use 9
Obtain basic QC stats from a BAM file, or profile how many reads were accepted by each rule
```
### get QC stats on whole bam AND find how many reads are clipped with high-quality clipped bases
### use the -x flag to produce no output (profiling only)
variant <bam> --min-clip 10 --min-phred 5 -q qcreport.txt -c clipcounts.txt -x
Rscript VariantBam/R/BamQCPlot.R -i qcreport.txt -o qcreport.pdf
```

##### Example Use 10
A user would like to extract only those reads supporting a particular allele at a variant site. This can be done by combining a small 
point-region at the variant site with a motif dictionary. 
Consider two alleles G and A at a site (e.g. 1:143250877), along with their adjacent sequences: 
GCAGAAT and GCAAAAT. To extract variant reads supporting the A allele:

```
## make the motifs file (include reverse complements) 
printf "GCAAAAT\nATTTTGC" > motifs.txt
## just look near the variant
k="1:143,250,677-143,251,077" 
r='{"":{"rules":[{"motif":"motifs.txt"}]}}'
g=1:143250877
variant <bam> -k $k -g $g -r $r -o mini.bam ## with JSON script
variant <bam> -k $k -g $g --motif motifs.txt -o mini.bam ## using command line shortcut
```

Because sequence information is required to match a motif, and reads do not contain the sequence information of their pair-mates, 
extracting all read pairs supporting a particular allele requires a two-pass solution:

```
## two pass solution
variant <bam> -k $k -g $g -r $r | cut -f1 | uniq > q.txt
printf "^@\n" >> q.txt ## keep the sam header too
samtools view <bam> $k -h | grep -f q.txt | samtools view - -b > mini.bam
```
                                             
This can be expanded for an arbitrary number of heterozygous sites, 
for instance to capture reads from a single haplotype:

```
### het.json
{
  "A" : {
      "region" : "1:132,250,677"
      "rules" : {[ "motif" : "motifsA.txt" ]
       }
  },
  "B" : {
      "region" : "1:182,250,325"
      "rules" : {[ "motif" : "motifsB.txt" ]
       }
  },
  
}
###
variant <bam> -r het.json -o mini.bam ## using JSON
variant <bam> -g 1:132,250,677 --motif motifsA.txt -g 1:182,250,325 --motif motifsB.txt ##using command line shortcut

```

Note that for the allele-specific extraction, there could be false negatives (reads not extracted) if a read has a sequencing error within the motif. 


Rules Script Syntax
===================

This section will describe the syntax used by VariantBam to specify the cascades of rules and regions 
applied to a BAM. Below is an example of a valid VariantBam JSON script:

```bash
   { 
     "reg1" : {
          "region" : "WG",
          "rules" : [{RULE_A, RULE_B}, {RULE_C, RULE_D, RULE_E}]
           }
   }
```

This can be read as "Accept a read that passes (RULE_A && RULE_B) OR (RULE_C && RULE_D && RULE_E)".

If no "rules" is supplied, it will default to "accept every".

### Region

The ``region`` keyword marks that what follows is a genomic region, 
which is either the keyword ``WG`` for whole genome, or a VCF, MAF, Callstats, BED file, or samtools-style string. If not specified, 
the default "WG" is applied. Regions are 
treated such that they will include any read who overlaps it, even partially. Optionally,
you can specify that your region of interest is a bit bigger than is actually in the file. You can do this by "padding"
the regions around the sites. For example:

```
### json
{ 
"" : {
       "region" : "myvcf.vcf",
       "pad" : 1000,
     }
}
### command line short-cut
variant $bam -g myvcf.vcf -P 1000
```

Alternatively, if supplying a region directly with the -l, -L, -g or -G flag, you can specify a padding with the -P flag. Note that this padding must be supplied
after a region flag is provided, and will be applied to the last supplied region flag (but with multiple regions, you can provide multiple ``-P`` flags). 

``variant <in.bam> -g myvcf.vcf -P 100``

You can also state that the region applies to reads who don't necessarily overlap the region, but their pair-mate does (called "mate-linking"). 
Note that rules that involve pair-mate information not located within the view read (e.g. mate mapq) are not considered.
Mate-linking is particularly useful for extracting all read PAIRS that cover a variant site.

```
### json
{ 
"" : {
       "region" : "myvcf.vcf",
       "pad" : 1000,
       "matelinked" : true,
     }
}
### command line shortcut
variant $bam -l myvcf.vcf -P 1000
```

```
### json to remove low (<= 10) MAPQ reads in bad region
{ 
"" : {
       "region" : "blacklist.bed",
       "pad" : 1000,
       "rules" : [{"mapq" : [0, 10]}]
       "exclude" : true
     }
}
### command line shortcut (for pure blacklisting)
variant $bam -G blacklist.bed -P 1000
### command line shortcut (for pure blacklisting, with mate linking)
variant $bam -L blacklist.bed -P 1000
```

### Global region

To reduce redundancy, you can name a region-rule set \"global\" anywhere in the stack,
and it will append that rule to everything below. For example, to exclude hardclipped, duplicate, qcfail and 
supplementary reads in every region, you would do:

```bash
{
  "global" : {
  	      "rules" : [{"hardclip" : false, "duplicate" : false, "qcfail" : false, "supplementary" : false}]
             },
  "A" : {
        "rules" : [{"isize" : [600, 0]}, {"clip" : [10,101], "mapq" : [1,60]}]
       },
  "B" : {
  	"region" : myvcf.vcf
   }
}
```

which keeps hardclipped reads, etc out from both of the subsequent regions (A and B)

### Rules

Rules are supplied as a list of criteria that a read must satisfy. VariantBam handles multiple rules in the following way. For each read, VariantBam 
will cycle through the rules within a region until the read satisfies a rule. When it 
does, it includes the read in the output and stops checking. The logic for the entire collection of 
rules is then as follows:

* On a given rule line, the read must satisfy ALL conditions (logical AND)

* Across different rules, the read nead only satisfy ONE rule (logical OR)

Below is an example which uses the ability of VariantBam to interpret VCFs and BED files,
and apply rules separately to them.

```bash
{
   "reg1" : {
              "region" : "myvcf.vcf",
              "pad" : 1000,
	      "matelinked" : true
            },
   "reg2" : {
              "region" : "myexonlist.bed",
              "matelinked" : true,
              "rules" : [{"isize" : [600,0], "mapped" : true, "mate_mapped" : false, "rg" : "H01PE.2"},
                         {"mapped" : false, "mate_mapped" : false},
                         {"hardclip" : true},
			 {"nm" : [1,101], "mapq" : [30, 100]}]
            }
}

The above JSON can be interpreted as a rule-cascade that, in one-pass of the BAM:
--Near VCF sites (to within 1000 bases)
   Keep reads interesecting region OR reads with mate-pairs that intersect region
--In exons:
   keep reads with: (isize outside of 0, 600 && with readgroup H01PE.2) OR (mapepd and mate unmapped) OR hardclipped OR (NM >= 1 && MAPQ >= 30)

```

###### Range rules
Range rules can be input as a two-element JSON array (``"mapq" : [30, 100]``) or a single value specificying the 
minimum accepted value (``"mapq" : 30``). To instead reject reads in this boundary, switch the order. Thus, 
``"mapq" : [100, 30]`` accepts only reads with MAPQ < 30 || MAPQ > 100.

###### Flag rules
Flag rules can be input using keywords like ``"mapped" : true`` or more versatily using the raw alginmetn flag ``"!flag" : 4``. Use
``flag`` to set all of the flags that must be turned on, and ``!flag`` for that must be turned off. Thus, ``"!flag" : 4``
requires that the "unmapped" bit be turned off, and so accepts only mapped reads.

###### Motif rules
A set of motifs can be supplied, so that only reads with (or without) the motif are accepted. A motif file is just a list of sequences in upper case, separated by newlines.
Reverse complements are not automatically considered,
so these must be explicitly provided. To include reads with a motif, use the ``--motif`` flag for simple rules, or in JSON specify ``"motif" : "motiffile.txt"``. To
exclude reads with a motif, use JSON key-value pair: ``"!motif" : "motiffile.txt"``.

Variant BAM can also filter based on the number of ``N`` bases in a read, with the ``nbases`` key, input as a range rule (``"nbaes" : [0,3]``)

###### Tag rules
Filters can be made based on the value of an alignment tag. Supported tags include "rg" (read-group), "nm" (number mismatches), and "xp" (number of supplementary alignments).

###### Cigar rules
Filters can be supplied to enforce a min (or max) number of insertions or deletions, or number of clipped reads. Note that this refers to the max element size. e.g. a CIGAR
of ``10M4D20M2D20M`` will pass the filter ``"del" : [0,4]`` but fail the filter ``"del" : [5, 100]``. The number of clipped bases is consider after trimming low-quality bases (if ``phred`` is supplied).

###### Subsample rule
Region-specific subsampling rates can be applied. For example, in region A you can set ``"subsample" : 0.5``, which will remove half of all reads that otherwise passed the other filters. If a
read falls into two regions, the region listed first in the JSON will apply.

""

Command Line Usage
==================

You can specify simple scripts directly on the command line:
```bash
### keep only paired reads, but not duplicates, in mate-linked region around 100 bp window at VCF sites
variant $bam -l myvcf.vcf -f 1 -F 1024 -P 100 -o out.bam
```

JSON scripts can also be supplied directly, just make sure to encase in single quotes and remove spaces:
```bash
variant $bam -g WG -r '{"":{"rules":[{"motif":"mymotifs.txt"}]}}' -o out.bam
```

### Full list of options

```
Usage: variant <input.bam> [OPTIONS] 

Description: Filter a BAM/CRAM file according to hierarchical rules

 General options
      --help                           Display this help and exit
  -v, --verbose                        Verbose output
  -c, --counts-file                    File to place read counts per rule / region
  -x, --counts-file-only               Same as -c, but does counting only (no output BAM)
  -r, --rules                          JSON script for the rules.
  -k, --proc-regions-file              Samtools-style region string (e.g. 1:1,000-2,000) or BED of regions to process
 Output options
  -o, --output-bam                     Output BAM file to write instead of SAM-format stdout
  -C, --cram                           Output file should be in CRAM format
  -T, --reference                      Path to reference. Required for reading/writing CRAM
  -h, --include-header                 When outputting to stdout, include the header.
  -s, --strip-tags                     Remove the specified tags, separated by commas. eg. -s RG,MD
  -S, --strip-all-tags                 Remove all alignment tags
 Filtering options
  -q, --qc-file                        Output a qc file that contains information about BAM
  -m, --max-coverage                   Maximum coverage of output. BAM must be sorted. Negative vals enforce min coverage.
 Region specifiers
  -g, --region                         Regions (e.g. myvcf.vcf or WG for whole genome) or newline seperated subsequence file.
  -G, --exclude-region                 Same as -g, but for region where satisfying a rule EXCLUDES this read. 
  -l, --linked-region                  Same as -g, but turns on mate-linking
  -L, --linked-exclude-region          Same as -l, but for mate-linked region where satisfying this rule EXCLUDES this read.
  -P, --region-pad                     Apply a padding to each region supplied with the region flags (specify after region flag)
 Command line rules shortcuts (to be used without supplying a -r script)
      --min-phred                      Set the minimum base quality score considered to be high-quality
      --min-clip                       Minimum number of quality clipped bases
      --max-nbases                     Maximum number of N bases
      --min-mapq                       Minimum mapping quality
      --min-del                        Minimum number of inserted bases
      --min-ins                        Minimum number of deleted bases
      --min-readlength                 Minimum read length (after base-quality trimming)
      --motif                          Motif file
  -R, --read-group                     Limit to just a single read group
  -f, --include-aln-flag               Flags to include (like samtools -f)
  -F, --exclude-aln-flag               Flags to exclude (like samtools -F)
```

Full list of available JSON rules
=================================

```
    #RULE           #EXAMPLE                   #DESCRIPTION OF EXAMPLE / FLAG
    flag            "flag" : 4                 Set the flag bits that must be ON
    !flag           "!flag" : 4                Set the flag bits that must be OFF
    motif           "motif" : seqs.txt         File containing substrings that must be present in the sequence.
    !motif          "!motif" : seqs.txt        File containing substrings that must NOT be present in the sequence.
    rg              "rg" : "H01PE.2"           Limit to just a single read-group 
    duplicate       "duplicate" : true         Read must be marked as optical duplicate
    supp            "supp" : false             Read must be primary alignment
    qcfail          "qcfail" : false           Read must note be marked as QC Fail
    fwd_strand      "fwd_strand" : true        Read must be mapped to forward strand
    rev_strand      "rev_strand" : true        Read must be mapped to reverse strand
    mate_fwd_strand "mate_fwd_strand" : true   Mate of read must be mapped to forward strand
    mate_rev_strand "mate_rev_strand" : true   Mate of read must be mapped to reverse strand
    mapped          "mapped" : true            Read must be unmapped
    mate_mapped     "mate_mapped" : true       Mate must be mapped
    subsample       "subsample" : 0.4          Subsample this region to at a certain rate
    ff              "ff" true                  Read pair must have forward-forward orientation
    rr              "rr" : true                Read pair must have reverse-reverse orientation
    fr              "fr" : true                Read pair must have forward-reverse orientation (proper)
    rf              "rf" : true                Read pair must have reverse-forward orientation
    ic              "ic" : true                Read pair must have inter-chromosomal mapping
    ... ALL RANGE RULES FOLLOW THE 3 INPUT OPTIONS ILLUSTRATED BELOW ... 
    ins             "ins"  : [5,101]           Number of inserted bases on the reads (from parsed CIGAR string)
                    "ins" : 5                  ... Take only reads with max insertion size of >= 5
                    "ins" : [101,5]            ... Take only reads with max insertion size NOT in [5,101] (e.g. 0-4)
    del             "del"  : [5,101]           Number of deleted bases relative to reference (from parsed CIGAR string). 
    nm              "nm" : [0,4]               NM tag from BAM (number of mismatches). e.g. must be 0-4 inclusive
    xp              "xp" : [0,4]               Number of supplementary aligments, with XP or XA tag from BAM (hold identity of supplementary alignments)
    isize           "isize" : [100,500]        Insert size, where all insert sizes are converted to positive.
    len             "len" : [80,101]           Length of the read following phred trimming. If phred trimming, don't count hardclips. If not, then HC count to length
    clip            "clip" : [0,5]             Number of clipped bases following phred trimming
    nbases          "nbases" : [0,5]           Removed reads that have within this range of N bases.
    phred           "phred" : [4,100]          Range of phred scores that are considered 'high-quality'
```

Attributions
============

VariantBam is developed and maintained by Jeremiah Wala (jwala@broadinstitute.org) --  Rameen Berkoukhim's lab -- Dana Farber Cancer Institute, Boston, MA. 

This project was developed in collaboration with the Cancer Genome Analysis team at the Broad Institute. Particular thanks to:
* Cheng-Zhong Zhang (Matthew Meyerson Lab)
* Marcin Imielinski
* Gad Getz
* Mara Rosenberg
* Esther Rheinbay 
* Gordon Saksena

[license]: https://github.com/jwalabroad/VariantBam/blob/master/LICENSE

[hlib]: https://github.com/samtools/htslib

[aho]: http://sourceforge.net/projects/multifast/

[ekg]: https://github.com/ekg/intervaltree

[samtools]: https://github.com/samtools/samtools

[bcf]: https://github.com/samtools/bcftools

[biop]: http://bioinformatics.oxfordjournals.org/content/early/2016/03/15/bioinformatics.btw111.full.pdf
