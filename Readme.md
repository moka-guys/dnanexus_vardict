# DNAnexus VarDict

## What does this app do?

This app is variant caller, Vardict. Exporting a vcf file per sample.

VarDict is an ultra sensitive variant caller for both single and paired sample variant calling from BAM files.
The app here uses single sample variant calling mode. 
VarDict is a variant calling program for SNV, MNV, indels (<120 bp default), and complex variants.  

This app accepts any BAM file. Bam file will be indexing via samtools prior to variant calling. 
By defult this app preforms localrealignment over indels on for more accurate allele frequencies of indels. 

BED files dictate the genomic regions inside which you want the variant analysis to be performed

## What are typical use cases for this app?

This app is used after sequencing alignment, to generate a vcf of variants identified (within region specified in a bed file).

## What inputs are required for this app to run?

This app requires Samtools and the following data:
- Compressed reference genome including `*.fa` and `*.fa.fai` (`*.tar.gz`)
- BAM files (`*.bam`)
- BED file of regions of intrest, for filtering output vcf (`*.bed`)

This app accepts the following inputs:
-f 	The threshold for allele frequency, default: 0.01 or 1%
-r	The minimum # of variance reads, default: 2
-B	The minimum # of reads to determine strand bias, default 2
-c	The column for chromosome
-S	The column for region start, e.g. gene start
-E	The column for region end, e.g. gene end
-g	The column for gene name, or segment annotation
-k 	Indicate whether to perform local realignment.  Default: 1 or true (For Ion or PacBio, 0 or False is recommended).
-N	(optional) The sample name to be used directly.  Will overwrite naming derived from bam file.

Extra Options Advanced inputs:
See Vardict help for full list of additional options.


## What does this app output?

This app outputs a vcf per sample detaling all variants identified. For detailed information about the analysis, consult the Vardict readme at:

https://github.com/AstraZeneca-NGS/VarDict 


## How does this app work?

VarDict, preforms local realignment and makes consensus calls of variants(SNP/Indel/Reference) from an indexed bam file for rgions specified in the supplied bedfile. A vcf file is output of each sample.

