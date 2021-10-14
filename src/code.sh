#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

#Grab inputs
dx-download-all-inputs --except ref_genome --parallel


# make output folder
mkdir -p ~/out/vardict_vcf/output
mkdir -p ~/out/vardict_vcf/interim
mkdir -p ~/bcftools_stats/QC/

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

# Get dockerised BCFtools from 001_ToolsReferenceData/Apps/Docker

BCFTOOLS_DOCKER_FILE_ID=project-ByfFPz00jy1fk6PjpZ95F27J:file-G5Z3Yk006yv6gGp6G7zFQG8j
BCFTOOLS_DOCKER_IMAGE_FILE=$(dx describe ${BCFTOOLS_DOCKER_FILE_ID} --name)
BCFTOOLS_DOCKER_IMAGE_NAME=$(tar xfO "${BCFTOOLS_DOCKER_IMAGE_FILE}" manifest.json | sed -E 's/.*"RepoTags":\["?([^"]*)"?.*/\1/')

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

# get the name of the reference genome, and take all before first "."
# should end up with GRCh38, hs37d5 etc - this should ensure reference build is described in the VCf header if not placed there by variant caller.
genomebuild=$(echo $ref_genome_name | cut -d. -f1)


# => genome/<ref>, genome/<ref>.ann, genome/<ref>.bwt, etc.
cd genome
for file in *
  do mv "$file" "${file/genome/$genomebuild}"
  done
cd ..
genome_file=`ls genome/*.fa`


# Run variant annotator for each Bam
mark-section "Run VarDict Variant Caller"
# loop through array of all bam files input, run VarDict for each bam file. 
for (( i=0; i<${#bam_file_path[@]}; i++ )); 
# show name of current bam file be run
do echo ${bam_file_prefix[i]}
# Index bam file input
samtools index ${bam_file_path[i]}
# run Vardict
/usr/bin/VarDict-1.8.2/bin/VarDict -th -G $genome_file -b ${bam_file_path[i]} $opts $bedfile_path | tee ~/out/vardict_vcf/interim/${bam_file_prefix[i]}.vardict.csv | /usr/bin/VarDict-1.8.2/bin/teststrandbias.R | /usr/bin/VarDict-1.8.2/bin/var2vcf_valid.pl -E -f $allele_freq $var2vcf_opts > ~/out/vardict_vcf/output/${bam_file_prefix[i]}.vardict.vcf
# add the reference genome into the last line of the header
sed -i "s/#CHROM/##REFERENCE=$genomebuild\n#CHROM/" ~/out/vardict_vcf/output/${bam_file_prefix[i1]}.vardict.vcf

# Run bcftools stats to produce stats for each file
sudo docker run -v /home/dnanexus:/home --rm ${BCFTOOLS_DOCKER_IMAGE_NAME} bcftools stats ~/out/vardict_vcf/output/${bam_file_prefix[i]}.vardict.vcf > ~/out/bcftools_stats/QC/${bam_file_prefix[i]}.vardict.stats

done

# Send output back to DNAnexus project
mark-section "Upload output"
dx-upload-all-outputs --parallel

mark-success
