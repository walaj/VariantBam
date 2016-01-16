[![Build Status](https://travis-ci.org/jwalabroad/VariantBam.svg?branch=master)](https://travis-ci.org/jwalabroad/VariantBam)

VariantBam: Filtering and profiling of next-generational sequencing data using region-specific rules
======================================================================================================

<div style="text-align:center"><img src="https://raw.githubusercontent.com/jwalabroad/VariantBam/master/labels_vb.png" width="200"></div>

**License:** [GNU GPLv3][license]

Installation
------------
I have succesfully built on Unix with GCC-4.8+

```
### if on Broad Institute servers, add GCC-4.9
reuse -q GCC-4.9

############## DOWNLOAD AND INSTALL BOOST ###############
############## (only if not already installed) ##########
git clone --recursive https://github.com/boostorg/boost.git
cd boost
./bootstrap.sh --with-libraries=regex,test
./b2

############### DOWNLOAD VARIANT BAM ############### 
git clone --recursive https://github.com/jwalabroad/VariantBam.git
cd VariantBam

############### COMPILE AND INSTALL ###############
./configure --with-boost=<path_to_boost>
make

############### QUICK START ############### 
mkdir -p tmp && cd tmp

## using the included test BAM (HCC1143)
VariantBam/src/variant test/small.bam -g 'X:1,000,000-1,100,000' -r mapq[10,100] -c counts.tsv -o mini.bam -v

## get help
VariantBam/src/variant --help

############## TL;DR ###################

## extract all reads and their mate-pairs that overlap SNP sites within 100 bp
rfile=<BED file, samtools-style string (e.g. "1:1,000,000-1,000,010"), or VCF>
variant <bam> -l $rfile -o mini.bam -P 100 -v

## mask regions (exclude reads and their pair-mates that overlap)
variant <bam> -L $rfile -o mini.bam -v 

## extract high-quality clipped reads
variant <bam> -r 'phred[4,100];clip[5,1000]' -o mini.bam -v

## subsample to max-coverage. BAM must be sorted
variant <bam> -m 100 -o mini.bam -v
```

Description
-----------

VariantBam is a tool to extract/count specific sets of sequencing reads from NGS sequencing files. To save money, 
disk space and I/O, one may not want to store an entire BAM on disk. In many cases, it would be more efficient to store only those read-pairs or
reads who intersect some region around the variant locations. Alternatively, if your scientific question is focused on only one aspect of the data (eg breakpoints), many 
reads can be removed without losing the information relevant to the problem. 

Additionally, VariantBam can be used as part of analysis pipelines that 

##### Example Use 1
Whole-genome analysis has been conducted on a BAM, generating VCF and MAF files. Ideally, these regions could be manually inspected
or reanalyzed without having to keep the entire BAM. Running VariantBam to extract only reads that overlap these events will allow
these regions to be rapidly queried, without having to keep the full BAM record.
```
### Extract all read PAIRS that interset with a variant from a VCF
variant $bam -l myvcf.vcf -o mini.bam
```

##### Example Use 2
In situations where the sequencing or library preparation quality is low, it may be advantageous
to remove poor quality reads before starting the analysis train. VariantBam handles this by optionally taking into
account Phred base-qualities when making a decision whether to keep a sequencing read. For instance, one might 
only be interested in high quality MAPQ 0 or clipped reads. VariantBam can be 
setup to apply unique Phred filters to different regions or across the entire genome, all with one-pass. 
```
### Extract only high quality reads with >= 50 bases of phred >=4 and MAPQ >= 1 and not duplicated/hardclip/qcfail
variant $bam -r 'phred[4,100];length[50,1000];mapq[1,60];!duplicate;!hardclip!;qcfail' -o mini.bam
```
##### Example Use 3
An NGS tool operates only on a subset of the reads (eg. structural variant caller using only clipped/discordant reads). Running VariantBam
to keep only these reads allows the tool to run much faster. This is particurlaly useful for facilitating a more rapid "build/test" cycle.
```
### Extract clipped, discordant, unmapped and indel reads
variant $bam  -r 'global@nbases[0,0];!hardclip;!supplementary;!duplicate;!qcfail;phred[4,100];%region@WG%discordant[0,1000];mapq[1,1000]%mapq[1,1000];clip[5,1000]%ins[1,1000];mapq[1,100]%del[1,1000];mapq[1,1000]' -o mini.bam
```
##### Example Use 4
A user wants to profile a BAM for quality. They would like to count the number of clipped reads in a BAM file, so long
as those reads have sufficient optical quality and mapping quality. VariantBam run with the -x flag for "counting only" 
will accomplish this.
```
### 
variant $bam -r 'clip[5,100];phred[4,100];mapq[10,100]' -x counts.tsv
```
##### Example Use 5
A team is only interested in variants in known cancer genes, and would like to analyze thousands of exomes and genomes. Running 
VariantBam to extract reads from only these genes, and sending the BAM files to compressed CRAM provides sufficient data reduction
to allow all of the relevant data to be stored on disk.
```
### Grab only reads from predefined regions. Strip unneccessary tags and convert to CRAM for maximum compression
variant $bam -g mygenes.bed -r all -C -o mini.cram -s BI,OQ
```
##### Example Use 6
A research team would like to extract only reads matching a certain motifs, but only if they have high optical quality. 
VariantBam with the ``motif`` rule will accomplish this with rapid O(n) efficiency for an arbitrarily large motif dictionary (where ``n`` is
the length of a read)
```
### 
variant $bam -r 'motif[mymotifs.txt];phred[4,100];length[20,1000]' -o mini.bam
```

##### Example Use 7
To reduce the size of the BAM, reads can be removed from centromeric and satellite repeat regions. These reads are rarely helpful for variant calling.
To remove reads that intersect a region, set the region as an inverse-region. In a VariantBam script, use ``!region`` or ``!mlregion``. For 
quick use on the command line, use ``-L`` or ``-G`` (opposites of ``-l`` and ``-g``).
```
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
A user would like to extract only those reads supporting a particular allele at a variant site. This can be done by combining a small 
point-region at the variant site with a motif dictionary. 
Consider two alleles G and A at a site (e.g. 1:143250877), along with their adjacent sequences: 
GCAGAAT and GCAAAAT. To extract variant reads supporting the A allele:

```
## make the motifs file (include reverse complements) 
printf "GCAAAAT\nATTTTGC" > motifs.txt
## just look near the variant
k="1:143,250,677-143,251,077" 
r='motif[motifs.txt]'
g=1:143250877
variant <bam> -k $k -g $g -r $r -o mini.bam
```

Because sequence information is required to match a motif, and reads do not contain the sequence information of their pair-mates, extracting all read pairs supporting a particular allele requires a two-pass solution:

```
## two pass solution                                                                                                                                                                                                                                                            
variant <bam> -k $k -g $g -r $r | cut -f1 | uniq > q.txt                                                                                                                                                                        
printf "^@\n" >> q.txt ## keep the sam header too                                                                                                                                                                                                                           
samtools view <bam> $k -h | grep -f q.txt | samtools view - -b > mini.bam                                                                                                                                                                                                                                  
```
                                             
This can be expanded for an arbitrary number of heterozygous sites, for instance to capture reads from a single haplotype:

```
r='region@1:143,250,677%motif[motifsA.txt]%region@1:183,250,677%motif[motifsB.txt]'
variant <bam> -r $r -o mini.bam
```

Note that for the allele-specific extraction, there will be false negatives (reads not extracted) if a read has a sequencing error within the motif. 

Tool comparison
---------------

In comparing with other avaiable BAM filtering tools, VariantBam provides the following additional features:

> 1. The ability to filter specifically on read clipping, orientation and insert size (all important for structural variation), while taking into account the per-base phred quality.
> 2. Use of interval trees to efficiently determine if a read or read mate overlaps a region.
> 3. The ability to provide different rules for different regions, and the ability to provide these regions as common variant files (VCF, MAF, BED)
> 4. Selecting reads by motif matching against a dictionary of any size
> 5. Ability to count numbers of reads that satisfy any number of user-defined properties
> 6. Read and write CRAM files
> 7. Selectively strip alignment tags
> 8. Support for sub-sampling to obtain a BAM file with a coverage limit

VariantBam is implemented in C++ and uses the HTSlibrary from Heng Li, a highly optimized C library used as the core of Samtools and BCFtools.

To get a full list of options, run ``variant --help``.

Rules Script Syntax
------

This section will describe the syntax used by VariantBam to specify the cascades of rules and regions 
applied to a BAM. Below is an example of a valid VariantBam script:

```bash
    ### this is a comment. The line code below defines filters to be applied to each region/rule
    region@WG
    !hardclip;mapped;mapped_mate;isize[0,600];!mapq[10,100]
    !hardclip;mapped;mapped_mate;clip[10,101]
```

##### Region

Let's look first at the ``region`` tag. The region@ keyword marks that what follows is a genomic region, 
which is either the keyword ``WG`` for whole genome, or a VCF, MAF, Callstats or BED file. Regions are 
treated such that they will include any read who overlaps it, even partially. Optionally,
you can specify that your region of interest is a bit bigger than is actually in the file. You can do this by "padding"
the regions around the sites. For example:

``pad[1000];region@myvcf.vcf``

Alternatively, if supplying a region directly with the -l, -L, -g or -G flag, you can specify a padding with the -P flag. Note that this padding will be applied to every region provided with one of these flags:

``variant <in.bam> -l myvcf.vcf -P 100``

You can also state that the region applies to reads who don't necessarily overlap the region, but their pair-mate does (called "mate-linking"). Note that this only applies to the ``all`` target.
This is particularly useful for extracting all read PAIRS that cover a variant site.

``pad[1000];mlregion@myvcf``

Note that the syntax is such that you must specify the file immediately after the @. 

A region can be inverted, so as to exclude any read that satisfies a rule on that region. In the rules-script, simply prepend the ``region`` with a !.
```bash
## exclude mapq = 0 reads and their mates in blacklist regions
pad[1000];!mlregion@blacklist.bed
mapq[0,0]
```

##### Rules

Rules are supplied as a list of criteria that a read must satisfy.
Note that you can take the complement of a condition 
by prefixing with a ``!``. For example:

```bash
    # do not include hardclipped reads, reads with isize > 600, or reads with mapq between 10 and 100.
    !hardclip;isize[0,600];!mapq[10,100]
    
    # an equivalent specification would be
    !hardclip;mapped;!isize[601,250000000];mapq[0,9]``
```

VariantBam handles multiple rules in the following way. For each read, VariantBam 
will cycle through the rules within a region until the read satisfies a rule. When it 
does, it includes the read in the output and stops checking. The logic for the entire collection of 
rules is then as follows:

On a given rule line, the read must satisfy ALL conditions (logical AND)

Across different rules, the read nead only satisfy ONE rule (logical OR)

To illustrate this, note that there is a small discrepancy in the first rule of the above. In the BAM format, 
unmapped reads and reads with unmapped mates are given an insert size of 0. However, in the same rule 
a condition is described to keep all reads with insert sizes 0-600 inclusive. Recalling the AND logic
within a rule, VariantBam will exclude the read, because it fails the ``mapped`` criteria.

Below is another example which uses the ability of VariantBam to interpret VCFs and BED files,
and apply rules separately to them.

```bash
    ### declare that region is a VCF file with pads of 1000 on either side of the variant.
    ### use the "mlregion" keyword to specify that pairs whose mate falls in the region belong to this rule
    mlregion@/home/unix/jwala/myvcf.vcf;pad[1000]
    #### I want to keep all the reads (this the default). Ill be explicit with the "all" keyword
    all
    #### A BED file which gives a list of exons. In here, I just want to keep "variant" reads pairs
    mlregion@/home/unix/jwala/myexonlist.bed 
    ## keep discordant reads
    !isize[0,600];
    ## keep only unmapped reads and their mates
    !mapped;!mapped_mate
    ## or keep if it is hardclipped
    hardclip
    ## keep reads with a mismatch to reference, but with high mapq
    nm[1,101];mapq[30,100]
```

##### Global

To reduce redundancy, you can also type a ``global@`` rule anywhere in the stack,
and it will append that rule to everything below. For example, to exclude hardclipped, duplicate, qcfail and 
supplementary reads in every region, you would do:

```bash
    global@!hardclip;!duplicate;!qcfail;!supplementary
    region@WG
    !isize[0,600]
    clip[10,101];mapq[1,60]
    region@myvcf.vcf
```

which is equivalent to

```bash
    region@WG
    !isize[0,600];!hardclip;!duplicate;!qcfail;!supplementary
    clip[10,101];mapq[1,60];!hardclip;!duplicate;!qcfail;!supplementary
    region@myvcf.vcf
    !hardclip;!duplicate;!qcfail;!supplementary
```
	
The global tag will apply through all of the regions. If you want to reset it for everything, just add ``global@all`` 
back onto the stack.

To make things run a little faster, you can set the order so that the more inclusive regions / rules are first. This only
applies if there is an overlap among regions. This is because VariantBam will move down the list of regions
that apply to this read and stop as soon as it meets an inclusion criteria. I prefer to start with a whole-genome region / rule
set, and then add more fine-mapped regions later.

##### Command Line Script

The usual method of inputing rules is with a VariantBam script as a text file (passed to
VariantBam with the ``-r`` flag). However, sometimes it is useful to not have to write an intermediate
file and just feed rules directly in. In that case, just pass a string literal to the ``-r, -g, -l`` flags, and VariantBam
will parse directly. ``-r`` will append a new rule, ``-g`` will append a new region and ``-l`` will append a new mate-linke regions. 
You can separate rule lines with either a new ``-r`` flag or with a ``%``. For instance, you might run something like the following:

```bash
variant big.bam -o small.bam -r 'global@!hardclip' -g WG -r '!isize[0,600];%clip[10,101];mapq[1,60]' -l 'myvcf.vcf' 
```

Note the single quotes so that it is interpreted as a string literal in BASH.

Full list of options
--------------------
```
Usage: variant <input.bam> [OPTIONS] 

Description: Filter a BAM/CRAM file according to hierarchical rules

 General options
      --help                           Display this help and exit
  -v, --verbose                        Verbose output
  -c, --counts-file                    File to place read counts per rule / region
  -x, --counts-file-only               Same as -c, but does counting only (no output BAM)
 Output options
  -o, --output-bam                     Output BAM file to write instead of SAM-format stdout
  -C, --cram                           Output file should be in CRAM format
  -T, --reference                      Path to reference. Required for reading/writing CRAM
  -h, --include-header                 When outputting to stdout, include the header.
  -s, --strip-tags                     Remove the specified tags, separated by commas. eg. -s RG,MD
  -S, --strip-all-tags                 Remove all alignment tags
 Filtering options
  -q, --qc-file                        Output a qc file that contains information about BAM
  -m, --max-coverage                   Maximum coverage of output file. BAM must be sorted. Negative values enforce a minimum coverage.
  -g, --region                         Regions (e.g. myvcf.vcf or WG for whole genome) or newline seperated subsequence file.  Applied in same order as -r for multiple
  -G, --exclude-region                 Same as -g, but for region where satisfying a rule EXCLUDES this read. Applied in same order as -r for multiple
  -l, --linked-region                  Same as -g, but turns on mate-linking
  -L, --linked-exclue-region           Same as -l, but for mate-linked region where satisfying this rule EXCLUDES this read.
  -r, --rules                          Script for the rules. If specified multiple times, will be applied in same order as -g
  -k, --proc-regions-file              Samtools-style region string (e.g. 1:1,000,000-2,000,000) or BED file of regions to proess reads from
  -P, --region-pad                     Apply a padding to each region supplied to variantBam with the -l, -L, -g or -G flags
```

Full list of available rules
----------------------------

```
    #RULE           #EXAMPLE             #DESCRIPTION OF EXAMPLE / FLAG 
    motif           motif[seqs.txt]  	 File containing substrings that must be present in the sequence.
    ins             ins[5,101]           Number of inserted bases on the reads (from parsed CIGAR string)
    del             del[10,101]          Number of deleted bases relative to reference (from parsed CIGAR string). 
    nm              nm[0,4]              NM tag from BAM (number of mismatches). e.g. must be 0-4 inclusive
    xp              xp[0,4]              Number of supplementary aligments, with XP or XA tag from BAM (hold identity of supplementary alignments)
    isize           isize[100,500]       Insert size, where all insert sizes are converted to positive.
    len             len[80,101]          Length of the read following phred trimming
    clip            clip[0,5]            Number of clipped bases following phred trimming
    nbases          nbases[0,5]          Removed reads that have within this range of N bases.
    phred           phred[4,100]         Range of phred scores that are "quality" 
    duplicate       duplicate            Read must be marked as optical duplicate 
    supp            !supp                Read must be primary alignment
    qcfail          !qcfail              Read must note be marked as QC Fail
    fwd_strand      fwd_strand           Read must be mapped to forward strand
    rev_strand      rev_strand           Read must be mapped to reverse strand
    mate_fwd_strand mate_fwd_strand      Mate of read must be mapped to forward strand
    mate_rev_strand mate_rev_strand      Mate of read must be mapped to reverse strand  
    mapped          !mapped              Read must be unmapped
    mapped_mate     mapped_mate          Mate must be mapped
    subsample       subsample[0.4]       Subsample this region to at a certain rate
    ff              ff                   Read pair must have forward-forward orientation
    rr              rr                   Read pair must have reverse-reverse orientation
    fr              fr                   Read pair must have forward-reverse orientation (proper)
    rf              rf                   Read pair must have reverse-forward orientation
    ic              ic                   Read pair must have inter-chromosomal mapping
    discordant      discordant[100,600]  Shortcut for !isize[100,600] || rr || ff || rf || ic (!discordant gives "proper" pairs)
```

Longer example rules script
---------------------------
```
##################################
### filters to apply to every read
##################################
global@!duplicate;!hardclip;!qcfail;!supplementary

#######################################################
## specify rules that apply to reads anywhere in genome
#######################################################
region@WG

#### Reads with abs(isize) outside of this range are kept. Also reads
## with non FR orientation. They must be both mapped though, and have
## min mapq of 1
discordant[0,1200];mapped;mate_mapped;mapq[1,100]

#### Reads with 1+ secondary alignments (as stored in XA or XP tag) are kept
xp[1,100]

### Reads without N in sequence, and with clip of 5+ bases (where each
### base must have phred of 4+ are kept). Must also be sufficiently
### long and have min mapq of 1
nbases[0,0];clip[5,101];phred[4,100];!length[0,20];mapq[1,100]

### Any read with non-zero mapq and an insertion (by CIGAR) of 1+ are kept
mapq[1,100];ins[1,100]

### Same, but for dels
mapq[1,100];del[1,100]

### any unmapped mate, and sufficient optical quality (>40 bases with
### phred 4+), is kept
!mapped;mate_mapped;phred[4,100];!length[0,40]

### any read with unmapped mate and mapq >= 10 is kept
!mate_mapped;mapped;mapq[10,100]

### can match sequences against a database of motifs. If read has
### sequence motif in the provided dictionary, it is kept
motif[VariantBam/test/motifs.txt]

################################################################
### can take ALL reads that overlap, or mate overlaps, a region. 
### the mlregion tag is for "mate-linked", which means a read passes
### if either it or its mate overlaps regions
#################################################################
mlregion@example.bed
all

### can specify MuTect regions with KEEP fields as input
#mlregion@LU-A08-43-Tumor.call_stats.txt
#all

### can specify VCF regions to keep. Here we also pad the regions by
### 500bp on each side
pad[500];mlregion@VariantBam/test/test.vcf
all
```

Attributions
------------
VariantBam is developed and maintained by Jeremiah Wala (jwala@broadinstitute.org) --  Rameen Berkoukhim's lab -- Dana Farber Cancer Institute, Boston, MA. 

This project was developed in collaboration with the Cancer Genome Analysis team at the Broad Institute. Particular thanks to:
* Cheng-Zhong Zhang (Matthew Meyerson Lab)
* Marcin Imielinski
* Gad Getz
* Mara Rosenberg
* Esther Rheinbay 
* Gordon Saksena

[license]: https://github.com/broadinstitute/variant-bam/blob/master/LICENSE

[BamTools]: https://raw.githubusercontent.com/wiki/pezmaster31/bamtools/Tutorial_Toolkit_BamTools-1.0.pdf

[API]: http://pezmaster31.github.io/bamtools/annotated.html

[hlib]: https://github.com/samtools/htslib

[snowt]: https://github.com/jwalabroad/SnowTools

[aho]: http://sourceforge.net/projects/multifast/
