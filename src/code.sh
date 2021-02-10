#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

#Grab inputs
dx-download-all-inputs --except ref_genome --parallel

# make output folder
mkdir -p ~/out/vardict_vcf/output

# Move inputs to home
# mv ~/in/bam_file/* ~/*

echo $allele_freq
echo $sample_name
echo $col_chr
echo $col_start
echo $col_end
echo $col_gene
echo $min_reads
echo $reads_bias
echo $local_realignment
echo $extra_options

#Construct opts sting
opts=" -f $allele_freq -c $col_chr -S $col_start -E $col_end -g $col_gene "
# add non-optional arguments to opts string
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

if [ "$extra_options" != "" ]; then
	opts="$opts $extra_options"
fi

# make folder for reference genome
mkdir genome
# download and untar reference genome
dx cat "$ref_genome" | tar zxvf - -C genome

# Parse the reference genome file name for the genome build used and set $genomebuild appropriately.
# This reference build will be added to the VCF header out put by VarDict.

if [[ $ref_genome_name =~ .*37.* ]]
then
	genomebuild="grch37"
elif [[ $ref_genome_name =~ .*19.* ]]
then
	genomebuild="hg19"
elif [[ $ref_genome_name =~ .*38.* ]]
then
	genomebuild="grch38"
elif [[ $ref_genome_name =~ .*20.* ]]
then
	genomebuild="hg20"
else
	echo "$ref_genome_name does not contain a parsable reference genome"
fi

cd genome
for file in *
  do rsync -a "$file" "${file/genome/$genomebuild}"
  done
cd ..

genome_file=$(ls genome/*.fa)


# Run variant annotator for each Bam
mark-section "Run VarDict Variant Caller"
# loop through array of all bam files input, run VarDict for each bam file.
for (( i=0; i<${#bam_file_path[@]}; i++ ));
# show name of current bam file be run
do echo ${bam_file_prefix[i]}
# Index bam file input
samtools index ${bam_file_path[i]}
# run Vardict
/usr/bin/vardict/vardict -G $genome_file -b ${bam_file_path[i]} $opts $bedfile_path | /usr/bin/vardict/teststrandbias.R | /usr/bin/vardict/var2vcf_valid.pl -E -f $allele_freq > ~/out/vardict_vcf/output/${bam_file_prefix[i]}.vardict.vcf
# add the reference genome into the last line of the header
sed -i "s/#CHROM/##REFERENCE=$genomebuild\n#CHROM/" ~/out/vardict_vcf/output/${bam_file_prefix[i]}.vardict.vcf

done

# Send output back to DNAnexus project
mark-section "Upload output"
dx-upload-all-outputs --parallel

mark-success
