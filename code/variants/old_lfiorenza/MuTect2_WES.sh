#!/bin/sh
#PBS -S /bin/sh
#PBS -l select=1:ncpus=4:mem=20g
#PBS -M laura.fiorenza@ieo.it
#PBS -m abe
#PBS -j oe

PATH=/hpcnfs/software/anaconda/anaconda3/bin:$PATH
source activate /hpcnfs/data/PGP/sharedCondaEnv/exome_v2.2/

cd $PBS_O_WORKDIR
mkdir ${sample_name_tum}/MuTect2

gatk Mutect2 \
-R /hpcnfs/data/PGP/reference_genomes/UCSC/hg38/Indexes/bwa_0.7.8/hg38.fa \
-I ./${sample_name_tum}/realignment/${sample_name_tum}_MarkDup_recal.bam \
-I ./${sample_name_norm}/realignment/${sample_name_norm}_MarkDup_recal.bam \
-tumor ${sample_name_tum} \
-normal ${sample_name_norm} \
-L /hpcnfs/data/PGP/exome/referenceBed/hg19/ss_v7/hg38/ss_v7_regions_hg38.bed  \
-O ./${sample_name_tum}/MuTect2/${sample_name_tum}_MuTect2.vcf.gz
