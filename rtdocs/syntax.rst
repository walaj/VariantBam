Syntax
------

This section will describe the syntax used by VariantBam to specify the cascades of rules and regions 
that are applied to the BAM[*]_. Below is an example of a valid VariantBam script:

.. code:: bash

    ### this is a comment. The line code below defines filters to be applied to each region/rule
    region@WG
    rule@!hardclip;!unmapped;!unmapped_mate;isize:[0,600];!mapq:[10,100]
    rule@!hardclip;!unmapped;!unmapped_mate;clip:[10,101]

Region
~~~~~~

Let's look first at the ``region`` tag. The region@ keyword marks that what follows is a region, 
which is either the keyword ``WG`` for whole genome, or a VCF, MAF, Callstats or BED file. Optionally,
you can specify that a region is a bit bigger than is actually in the file. You can do this by "padding"
the regions around the sites. For example:

``region@myvcf.vcf,pad:1000``

You can also state that the region applies to reads who don't necessarily overlap the region, but their pair-mate does.

``region@myvcf,pad:1000,mate``

Note that the syntax is such that you must specify the file immediately after the @, following by other options
in any order.

Rules
~~~~~

The next two lines specify a pair of rules, marked with the ``rule@`` tag. 
The default rule is to include every read, and the conditions are meant to be 
thought of as exclusion criteria. You can take the "opposite" of a condition by prefixing
with a ``!``. For example, the first rule in the above example states:

Keep all reads EXCEPT any read that satisfies the following: Hardclipped, is unmapped, has unmapped mate,
has insert size greater than 600, does NOT have mapping quality between 10 and 100. Thus, we are going to get low mapping 
quality discordant reads from this query. And equivalent specification would be:

``rule@!hardclip;!unmapped;!unmapped_mate;isize:[0,600];mapq:[0,9]``

VariantBam handles multiple rules in the following way. For each read, VariantBam 
will cycle through the rules within a region until the read satisfies a rule. When it 
does, it includes the reads and stops checking. The logic for the entire collectoin is then as follows:

On a given rule line, the read must satisfy ALL conditions (logical AND)

Across different rules, the read nead only satisfy ONE rule (logical OR)

To illustrate this, note that there is a small discrepancy in the first rule of the above. In the BAM format, 
unmapped reads and reads with unmapped mates are given an insert size of 0. However, in the same rule 
a condition is described to keep all reads with insert sizes 0-600 inclusive. Recalling the AND logic
within a rule, VariantBam will exclude the read, because it fails the ``!unmapped`` criteria.

Below is another example which uses the ability of VariantBam to interpret VCFs and BED files,
and apply rules separately to them.

.. code:: bash

    ### declare that my region is a VCF file with pads of 1000 on either side of the variant.
    ### use the "mate" keyword to specify that pairs whose mate falls in the region belong to this rule
    region@/home/unix/jwala/myvcf.vcf,mate,pad:1000
    #### I want to keep all the reads (this the default). Ill be explicit with the "every" keyword
    rule@every
    #### I might also have a BED file which gives a list of exons. In here, I just want to keep "variant" reads
    #### so I can specify something like:
    region@/home/unix/jwala/myexonlist.bed 
    rule@y!isize:[100,600];!unmapped;!unmapped_mate

Global
~~~~~~

To make things more clear and reduce redundancy, you can also type a ``global@`` rule anywhere in the stack,
and it will append that rule to everything below. For example, to exclude hardclipped, duplicate, qcfail and 
supplementary reads in every region, you would do:

.. code:: bash

    global@!hardclip;!duplicate;!qcfail;!supplementary
    region@WG
    rule@!isize:[0,600]
    rule@clip:[10,101];mapq:[1,60]
    region@myvcf.vcf

is equivalent to

.. code:: bash

    region@WG
    rule@!isize:[0,600];!hardclip;!duplicate;!qcfail;!supplementary
    rule@clip:[10,101];mapq:[1,60];!hardclip;!duplicate;!qcfail;!supplementary
    region@myvcf.vcf
    rule@!hardclip;!duplicate;!qcfail;!supplementary
	
The global tag will apply through all of the regions. If you want to reset it for everything, just add ``global@every`` 
back onto the stack.

To make things run a little faster, you can set the order so that the more inclusive regions / rules are first. This only
applies if there is an overlap among regions. This is because VariantBam will move down the list of regions
that apply to this read and stop as soon as it meets an inclusion criteria. I prefer to start with a whole-genome region / rule
set, and then add more fine-mapped regions later.

Command Line Script
~~~~~~~~~~~~~~~~~~~

The usual method of inputing rules is with a VariantBam script as a text file (passed to
VariantBam with the ``-r`` flag). However, sometimes it is useful to not have to write an intermediate
file and just feed rules directly in. In that case, just pass a string literal to the -r flag, and VariantBam
will parse directly. Just separate lines with a ``%``. For instance, you might run something like the following:

``variant -i big.bam -o small.bam -r 'global@!hardclip%region@WG%rule@!isize:[0,600];%rule@clip:[10,101];mapq:[1,60]%region@myvcf.vcf'``

Note the single quotes so that it is interpreted as a string literal in BASH.

.. [*] A standard format like JSON would be better and may be implemented in the future.
