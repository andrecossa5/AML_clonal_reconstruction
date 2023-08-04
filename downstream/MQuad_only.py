#!/usr/bin/python

# MQuad feature selection

########################################################################

# Code
import argparse
import os

##

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='vireoSNP',
    description='Wrapper script within MQuad.'
)

# Add arguments

# Path_main
my_parser.add_argument(
    '-p', 
    '--path_main', 
    type=str,
    default='..',
    help='Path to project dir. Default: .. .'
)

# Sample
my_parser.add_argument(
    '--sample', 
    type=str,
    default='MDA_clones',
    help='Sample to use. Default: MDA_clones.'
)

# GBC
my_parser.add_argument(
    '--GBC', 
    action='store_true',
    help='Read and use GBC information. Default: False.'
)

# ncores
my_parser.add_argument(
    '--ncores', 
    type=int,
    default=None,
    help='n cores to use.'
)

# Parse arguments
args = my_parser.parse_args()


##

# path_main = '/Users/IEO5505/Desktop/mito_bench/'
# sample = 'MDA_clones'
# range_clones = np.arange(2, 15)
# chosen = None
# filtering = 'MQuad'
# min_cov_treshold = 50
# min_cell_number = 10
# p_treshold = 0.8
# GBC = True
# ncores = 4

##

path_main = args.path_main
sample = args.sample
GBC =  args.GBC
ncores = args.ncores


##

####################################################################

# Preparing run: import code, set paths

# Code
from mito_utils.utils import *
from mito_utils.preprocessing import *
from mito_utils.plotting_base import *
from mito_utils.embeddings_plots import *
from mito_utils.dimred import *
from mito_utils.clustering import *

# Paths 
path_data = os.path.join(path_main, 'data')
path_results = os.path.join(path_main, 'results', 'vireo')

####################################################################

# Main
def main():

    # Create folder
    make_folder(path_results, sample, overwrite=True)
    
    if GBC:
        afm = read_one_sample(path_data, sample=sample, with_GBC=True)
        n_all_clones = len(afm.obs['GBC'].unique())
    else:
        afm = read_one_sample(path_data, sample=sample, with_GBC=False)

    # Filter cells and vars
    _, a = filter_cells_and_vars(
        afm,
        sample=sample,
        filtering='MQuad', 
        min_cell_number=0,
        min_cov_treshold=50,
        nproc=ncores, 
        path_=os.path.join(path_results, sample)
    )

    # Save vars
    pd.Series(a.var_names).to_csv(os.path.join(path_results, f'{sample}_MQuad_vars.csv'))



