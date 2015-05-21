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
#motif[/cga/fh/pcawg_pipeline3/modules/VariantBam/data/abc_38_both_pm_update.uniq]
#motif[/cga/fh/pcawg_pipeline3/modules/VariantBam/data/abc_v14.uniq]

################################################################
### can take ALL reads that overlap, or mate overlaps, a region. 
### the mlregion tag is for "mate-linked", which means a read passes
### if either it or its mate overlaps regions
#################################################################
#mlregion@/cga/fh/pcawg_pipeline3/modules/VariantBam/data/hla.bed
#all

### can specify MuTect regions with KEEP fields as input
#mlregion@LU-A08-43-Tumor.call_stats.txt
#all

### can specify VCF regions to keep. Here we also pad the regions by
### 500bp on each side
#pad[500];mlregion@/home/unix/jwala/myvcf.vcf
#all
