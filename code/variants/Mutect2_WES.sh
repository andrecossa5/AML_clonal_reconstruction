#!/bin/sh
#PBS -N re_Mutect2
#PBS -l select=1:ncpus=4:mem=20g
#PBS -M andrea.cossa@ieo.it
#PBS -m ae
#PBS -e /hpcnfs/scratch/PGP/acossa/AML_clonal_reconstruction/code/variants/mutect_2_e.txt
#PBS -o /hpcnfs/scratch/PGP/acossa/AML_clonal_reconstruction/code/variants/mutect_2_o.txt
#PBS -q nocg_workq

# Env
PATH=/hpcnfs/software/anaconda/anaconda3/bin:$PATH
source activate /hpcnfs/data/PGP/sharedCondaEnv/exome_v2.2/

# Params
sample_name=sAML1
path_main=/hpcnfs/scratch/PGP/acossa/AML_clonal_reconstruction/
cd $path_main/code/variants/ 
mkdir $path_main/results_and_plots/variants/$sample_name/

gatk Mutect2 \
-R /hpcnfs/data/PGP/reference_genomes/UCSC/hg38/Indexes/bwa_0.7.8/hg38.fa \
-I $path_main/data/bams/$sample_name/D/${sample_name}_D_MarkDup_recal.bam \
-I $path_main/data/bams/$sample_name/N/${sample_name}_N_MarkDup_recal.bam \
-tumor ${sample_name}_D \
-normal ${sample_name}_N \
-L /hpcnfs/data/PGP/exome/referenceBed/hg19/ss_v7/hg38/ss_v7_regions_hg38.bed  \
-O $path_main/results_and_plots/variants/$sample_name/MuTect2.vcf.gz
