#!/bin/bash

COUNTS_ONLY=1
COUNTS_AND_BAM=0

## example to get counts only
if [[ $COUNTS_ONLY -ne 0 ]]; then
options=(
    --input             /seq/picard_aggregation/G14856/TCGA-05-4432-01A-01D-1931-08/v2/TCGA-05-4432-01A-01D-1931-08.bam
    --output-bam        small.bam
    --rules             example.rules.txt
    --proc-regions-file example.bed
    --counts-file-only  counts.tsv
    --verbose
)
echo "...Counting reads only"
sleep 1
variant ${options[*]}
fi

## example to get counts and get mini-bam
## can uncomment below to strip some aligment tags or all, to make BAM smaller
if [[ $COUNTS_AND_BAM -ne 0 ]]; then
options=(
    --input             /seq/picard_aggregation/G14856/TCGA-05-4432-01A-01D-1931-08/v2/TCGA-05-4432-01A-01D-1931-08.bam
    --output-bam        small.bam
    --rules             example.rules.txt
    --proc-regions-file example.bed
    --counts-file       counts.tsv
    ##--strip-tags        OQ,AM,SM
    ##--strip-all-tags    
    --verbose
)
echo "...Counting and outputting BAM"
sleep 1
variant ${options[*]}
fi
