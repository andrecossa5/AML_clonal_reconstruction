"""
This takes some .vcf file in input and assess its properites after some filtering procedures. 
Finally selected variants are saved.  
"""

import os
import numpy as np
import pandas as pd
import pysam 

from Cellula.plotting._plotting_base import *


##


# Utils  
def process_VCF(path):
    """
    Process a VCF file.
    """

    vcf = pysam.VariantFile(path)

    mut_name_l = []
    ref_l = []
    alt_l = [] 
    tumor_AD_l = []
    tumor_AF_l = []
    tumor_DP_l = []
    normal_AD_l = []
    normal_AF_l = []
    normal_DP_l = []
    chrom_l = []
    gen_range_l = []
    length_l = []

    i = 0

    for i, r in enumerate(vcf):
        
        if len(r.alleles) == 2:
            ref, alt = r.alleles
            ref_l.append(ref)
            alt_l.append(alt)
            tumor_r = r.samples['sAML1_D']
            normal_r = r.samples['sAML1_N']
            tumor_AD_l.append(tumor_r["AD"])
            tumor_AF_l.append(tumor_r["AF"])
            tumor_DP_l.append(tumor_r["DP"])
            normal_AD_l.append(normal_r["AD"])
            normal_AF_l.append(normal_r["AF"])
            normal_DP_l.append(normal_r["DP"])
            chrom_l.append(r.chrom) 
            gen_range_l.append((r.start, r.stop))
            length_l.append(r.rlen)
            xpos = r.start if r.rlen == 1 else f'{r.start}-{r.stop}'
            mut_name_l.append(f'{ref}>{alt},{r.chrom},{xpos}')
        
        else:
            pass

        i += 1

    D = {
        'id' : mut_name_l,
        'ref' : ref_l,
        'alt' : alt_l,
        'chrom' : chrom_l,
        'gen_range' : gen_range_l,
        'length' : length_l,
        'tumor_AD' : tumor_AD_l,
        'tumor_AF' : tumor_AF_l,
        'tumor_DP' : tumor_DP_l,
        'normal_AD' : normal_AD_l,
        'normal_AF' : normal_AF_l,
        'normal_DP' : normal_DP_l
    }

    print(f'Total number of mutations called: {i}')

    vcf.close()

    return pd.DataFrame(D)


##


# Set paths 
path_main = '/Users/IEO5505/Desktop/AML_clonal_reconstruction/'
path_data = path_main + 'data/'
path_results = path_main + 'results_and_plots/variants/'
path_viz = path_main + 'results_and_plots/visualization/'

# Read data 
old = process_VCF(path_data + 'vcfs/sAML1/sAML1_D_MuTect2.vcf.gz')
new = process_VCF(path_results + 'sAML1/MuTect2.vcf.gz')
processed = process_VCF(path_results + 'sAML1/MuTect2.vcf.selected.gz')

old = old.set_index('id')
new = new.set_index('id')

# Filter and intersection
# vars_old = set(old.query('length == 1 and tumor_AF > 0.1').index)
# vars_new = set(new.query('length == 1 and tumor_AF > 0.1').index)
vars_old = set(old.loc[old['length'] == 1 & old['tumor_AF'].map(lambda x: x[0] > 0.1)].index)
vars_new = set(new.loc[new['length'] == 1 & new['tumor_AF'].map(lambda x: x[0] > 0.1)].index)

len(vars_old) 
len(vars_new)
len(vars_old & vars_new)

# Processed filter
processed = processed.set_index('id')
processed.columns

# Read old output from chiara
df = pd.read_excel(path_data + 'processed/sAML1.xlsx')
df = df.loc[:, 
                [
                    'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 
                    'Tumor_Seq_Allele2', 't_alt_count', 't_ref_count', 'n_alt_count', 
                    'n_ref_count', 'tumor_f'
                ]
    ].rename(columns={'Reference_Allele':'ref', 'Tumor_Seq_Allele2':'alt', 'Chromosome':'chrom', 'tumor_f':'tumor_AF'})
df['gen_range'] = [ (x, y) for x, y in zip(df['Start_Position'], df['End_Position']) ]
df = df.drop(columns=['Start_Position', 'End_Position'])
df['tumor_AD'] = [ (x, y) for x, y in zip(df['t_ref_count'], df['t_alt_count']) ]
df['normal_AD'] = [ (x, y) for x, y in zip(df['n_ref_count'], df['n_alt_count']) ]
df = df.loc[:, ~df.columns.str.contains('count')]

xpos = []
for i in range(df.shape[0]):
    start, stop = df['gen_range'][i]
    if stop == start:
        xpos.append(str(start-1))
    else:
        xpos.append(f'{start}-{stop}')

df['id'] = [ f"{df['ref'][i]}>{df['alt'][i]},{df['chrom'][i]},{xpos[i]}"  for i in range(df.shape[0]) ]
df = df.set_index('id')

# Intersect
old = set(df.index)
new = set(processed.index)

len(old)
len(new)
common = old & new
len(common)

# View VAF 
processed = processed.loc[list(common), :]
df = df.loc[list(common), :]

vafs_processed = processed['tumor_AF'].map(lambda x: x[0])
vafs_df = df['tumor_AF']

np.sum((vafs_processed - vafs_df) < 0.01)


##


# Join CNVs info 
df_cnvs = pd.read_csv(
    path_data + 'processed/copy_number_from_vcf_prova.csv', 
    sep='\t', 
    index_col=0
)
processed.columns

# Process to retrieve PyClone input
df_ = processed.loc[:, ['tumor_AD', 'gen_range']]
df_['ref_counts'] = df_['tumor_AD'].map(lambda x: x[0])
df_['alt_counts'] = df_['tumor_AD'].map(lambda x: x[1])
df_ = df_.drop(columns='tumor_AD')
df_['sample_id'] = 'sAML1'
df_['tumor_content'] = 0.81

df_cnvs = df_cnvs.drop(columns='gene')

major = []
minor = []

for i in range(df_.shape[0]):
    chrom = df_.index[i].split(',')[1]
    start, end = df_['gen_range'][i]
    if df_cnvs.query('chromosome == @chrom').shape[0] != 0:
        d_chrom = df_cnvs.query('chromosome == @chrom')
        if d_chrom.query('start >= @start and end <= @end').shape[0] != 0:
            print(df_.index[i])
            cns = df_cnvs.loc[:, ['cn1, cn2']]
            print('Hello')
        else
            print('Muts is not in an available genomic range...')



