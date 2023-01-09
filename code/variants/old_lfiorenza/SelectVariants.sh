#!/bin/sh
#PBS -S /bin/sh
#PBS -l select=1:ncpus=1:mem=16g
#PBS -M laura.fiorenza@ieo.it
#PBS -m abe
#PBS -j oe

PATH=/hpcnfs/software/anaconda/anaconda3/bin:$PATH
source activate /hpcnfs/data/PGP/sharedCondaEnv/exome_v2.2/

cd $PBS_O_WORKDIR
cd ${sample_name_tum}
cd MuTect2

gatk SelectVariants \
-R /hpcnfs/data/PGP/reference_genomes/UCSC/hg38/Indexes/bwa_0.7.8/hg38.fa \
-L /hpcnfs/data/PGP/exome/referenceBed/hg19/ss_v7/hg38/ss_v7_regions_hg38.bed  \
-V ./${sample_name_tum}_MuTect2.vcf.flagged.gz \
-O ./${sample_name_tum}_MuTect2.vcf.selected.gz \
--exclude-filtered true 
