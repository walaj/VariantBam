Introduction
------------

VariantBam is a tool to extract interesting sequencing reads from a BAM file. VariantBam 
was developed to be a one-pass solution for the various needs of different NGS tools. For instance,
an indel tool might be interested in MAPQ > 0 reads in 30,000 candidate regions of interest, 
a SNP tool might find a different 20,000, and an SV tool might be interested in only discordant or high-quality 
clipped reads across the BAM (where high-quality means they are not clipped to do low Phred quality). Alternatively, 
to save money/space one may not want to store an entire BAM on disk after all the relevant VCF, MAFs, etc have been created. 
Instead it would be more efficient to store only those read-pairs who intersect some region around the variant locations. 
VariantBam is designed to handle all of these situations with a single pass through the BAM.

VariantBam is implemented in C++ and relies on the BamTools API (Derek Barnett, (c) 2009) for BAM I/O. 
It is worth mentioning the capabilities of the BamTools command line ``bamtools filter`` here, 
which may provide a solution more to your needs than VariantBam. ``bamtools filter`` allows you to 
specify a set of filters in JSON format to be applied to a BAM. See the Bamtools documentation_ for more detail. 
Under what situations would you use ``bamtools filter``, and when would you use VariantBam?

1. Extract all MAPQ 0 reads from a BAM - Either tool (prefer ``bamtools filter``)
2. Extract all reads in read group A - ``bamtools filter``
3. Extract all reads with NM tag >= 4 - Either tool (prefer ``bamtools filter``)
4. Extract all reads with NM tag >= 4 in exons - VariantBam.
5. Remove all reads with insert sizes between 100 and 600 bp - VariantBam
6. Extract all reads and mate within 1000bp of a variant or set of genes - VariantBam
7. Extract only high-quality reads with N bases beyong phred score X - VariantBam
8. Reduce a BAM to only high quality reads around your MAFs, VCFs and BED files - VariantBam

A manuscript for VariantBam is currently under preparation.

.. _documentation https://raw.githubusercontent.com/wiki/pezmaster31/bamtools/Tutorial_Toolkit_BamTools-1.0.pdf

