# DNAnexus VarDict v1.3.2

## What does this app do?
This app performs variant calling using the VarDict variant caller, calling SNV, MNV, indels (<120 bp default), and complex variants.

## What are typical use cases for this app?
VarDict is used to detect variation (CNV and SNV) in somatic cancer testing. 
This app can be used to test samples being processed by the SWIFT amplicon panels. 
The output vcf will be uploaded into Ingenuity for annotation and filtering.


## What inputs are required for this app to run?
This app requires the following inputs:

- Compressed reference genome including `*.fa` and `*.fa.fai` (`*.tar.gz`).  The file name should contain the genome build number (37, 19, 38, 20) so that this can be parsed and the appropriate metadata can be added to the VCF header as required when processing with QCI.
- BAM file(s) (`*.bam`). Multiple BAM files can be provided, producing one VCF per sample (multisample variant calling is not performed).
- BED file of regions of interest, for filtering output vcf (`*.bed`)

This following parameters can be passed to the app:
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
- Extra options can be entered as a string and include the option and value in the following format e.g `-I 200 (to specify an indel size of 200bp)`
- See VarDict help for full list of additional options.


## How does this app work?
The app loops through the array of input BAM files and for each sample: 
- The app uses Samtools to index each BAM file 
- The app then uses VarDict to perform local realignment and call variants from the indexed bam file for the genomic regions specified in the supplied bed file. [This VarDict repository was cloned at this point into the app](https://github.com/AstraZeneca-NGS/VarDict/tree/328e00a1166abe4406020a9af12ca816a93517be).A number of scripts are applied in this process:
  - vardict.pl
  - teststrandbias.R 
  - var2vcf_valid.pl - Convert the output into validated VCF file
  - The reference build used is parsed from the reference FASTA file name and added to the VCF header.

- In addition to the parameters stated VarDict applies additional filters including:
  - -I The indel size. Default =  120bp
  - -q The phred score for a base to be considered a good call.  Default = 25 (for Illumina)

## What does this app output?
This app outputs one uncompressed vcf file (.vcf) per sample detailing all called variants within the regions specified in the BED file. 

vcf files are output to `/output`

For detailed information about the analysis, consult the [VarDict readme](https://github.com/AstraZeneca-NGS/VarDict)
