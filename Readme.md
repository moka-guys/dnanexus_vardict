# DNAnexus VarDict v1.2

## What does this app do?

This app performs variant calling using Vardict after sequencing alignment.
Vardict requires a BAM file and a bed file defining which regions are to be assessed. 
Vardict calls SNV, MNV, indels (<120 bp default), and complex variants.

This app can take mulitiple bam files as input, preforming variant calling on multiple samples. 
BED files dictate the genomic regions inside which you want the variant analysis to be performed


## What are typical use cases for this app?

VarDict is used to detect variation (CNV and SNV) in somatic cancer testing. 
This is used to test samples being processed by the SWIFT amplicon panels. 
A vcf file of variants identified (within region specified in a bed file) is generated for each input sample. 
The output vcf will be uploaded into Ingenuity for annotation and filtering.


## What inputs are required for this app to run?

This app the following data:

- Compressed reference genome including `*.fa` and `*.fa.fai` (`*.tar.gz`)
- BAM file(s) (`*.bam`)
- BED file of regions of intrest, for filtering output vcf (`*.bed`)

This app accepts the following inputs:

- -f 	The threshold for allele frequency, default: 0.01 or 1%
- -r	The minimum # of variance reads, default: 2
- -B	The minimum # of reads to determine strand bias, default 2
- -c	The bed file column for chromosome, default 1
- -S	The bed file column for region start, e.g. gene start, default 2
- -E	The bed file column for region end, e.g. gene end, default 3
- -g	The bed file column for gene name, or segment annotation, default 4
- -k 	Indicate whether to perform local realignment.  Default: 1 or true (For Ion or PacBio, 0 or False is recommended).
- -N	(optional) The sample name to be used directly.  Will overwrite naming derived from bam file.

Extra Options Advanced inputs:

- Extra options should be entered as a string and incude the option and value in the following format e.g -I 200 (to specify an indel size of 200bp)
- See Vardict help for full list of additional options.


## What does this app output?

This app outputs a vcf per sample detaling all variants identified. For detailed information about the analysis, consult the Vardict readme at:

https://github.com/AstraZeneca-NGS/VarDict 


## How does this app work?

The app loops through the array of input bam files and does the following steps for each sample: 
- The app uses Samtools to index each bam file 
- The app then uses VarDict, to preform local realignment and call variants from the indexed bam file for the genomic regions specified in the supplied bedfile. 
- A vcf file is output of each sample.

