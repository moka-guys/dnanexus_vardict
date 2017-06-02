#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

#Grab inputs
dx-download-all-inputs --except ref_genome --parallel

# Move inputs to home
# mv ~/in/bam_file/* ~/*

echo $allele_freq
echo $min_reads2
echo $sample_name
echo $col_chr
echo $col_start
echo $col_end
echo $col_gene
echo $min_reads
echo $reads_bias
echo $local_realignment
echo $extra_options

opts=" -f $allele_freq -c $col_chr -S $col_start -E $col_end -g $col_gene "

if [ "$min_reads" != "" ]; then
  opts="$opts -r $min_reads"
fi

if [ "$reads_bias" != "" ]; then
	opts="$opts -B $reads_bias" 
fi

if [ "$sample_name" != "" ]; then
	opts="$opts -N $sample_name" 
fi

if [ "$local_realignment" == false ]; then
	opts="$opts -k 0" 
fi

# un-package reference genome
mkdir genome
dx cat "$ref_genome" | tar zxvf - -C genome  
# => genome/<ref>, genome/<ref>.ann, genome/<ref>.bwt, etc.

# rename genome files to grch37 so that the VCF header states the reference to be grch37.fa, which then allows Ingenuity to accept the VCFs (otherwise VCF header would have reference as genome.fa which Ingenuity won't accept)
mv  genome/*.fa  genome/grch37.fa
mv  genome/*.fa.fai  genome/grch37.fa.fai
# mv genome.dict grch37.dict
genome_file=`ls genome/*.fa`


# Run variant annotator for each Bam
mark-section "Run VarDict VariantAnnotator"
for (( i=0; i<${#bam_file_path[@]}; i++ )); 
do echo ${bam_file_prefix[i]}
samtools index ${bam_file_path[i]}
/usr/bin/vardict/vardict -G $genome_file -b ${bam_file_path[i]} $opts $bedfile_path | /usr/bin/vardict/teststrandbias.R | /usr/bin/vardict/var2vcf_valid.pl -E -f $allele_freq > ${bam_file_prefix[i]}.vardict.vcf
done

# Send output back to DNAnexus project
mark-section "Upload output"
mkdir -p ~/out/vardict_vcf/
mv ./*.vardict.vcf ~/out/vardict_vcf/

dx-upload-all-outputs --parallel

mark-success
