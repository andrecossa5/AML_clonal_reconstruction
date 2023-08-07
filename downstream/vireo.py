#!/usr/bin/python

# Clonal inference with VireoSNP

########################################################################

# Code
import argparse
import os

##

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='vireoSNP',
    description='Wrapper script within vireoSNP vignette.'
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
    default='sAML1',
    help='Sample to use. Default: sAML1.'
)

# Range n
my_parser.add_argument(
    '--range', 
    type=str,
    default='2:15',
    help='Range of n_clones to search from. Default: 2:15.'
)

# min_cov_treshold
my_parser.add_argument(
    '--min_cov_treshold', 
    type=int,
    default=50,
    help='min_cov_treshold.'
)

# p_treshold
my_parser.add_argument(
    '--p_treshold', 
    type=float,
    default=0.8,
    help='Treshold use to convert clonal assignment to crisp labels.'
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

# path_main = '/Users/IEO5505/Desktop/AML_clonal_reconstruction'
# sample = 'sAML1'
# range_clones = np.arange(2, 15)
# min_cov_treshold = 50
# p_treshold = 0.8
# ncores = 8

##

path_main = args.path_main
sample = args.sample
start, stop = args.range.split(':')
range_clones = range(int(start), int(stop))
min_cov_treshold = args.min_cov_treshold
p_treshold = args.p_treshold
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
from mito_utils._vireo import *

# Paths 
path_data = os.path.join(path_main, 'data')
path_results = os.path.join(path_main, 'results', 'clustering')
path_variants = os.path.join(path_main, 'results', 'variants_selection')

####################################################################

# Main
def main():

    # Load data
    T = Timer()
    T.start()
    
    t = Timer()

    # Create folder
    make_folder(path_results, sample, overwrite=False)
    logger = set_logger(os.path.join(path_results, sample), 'log.txt')
    
    # Set logger
    logger.info(
        f""" 
        Execute clones unsupervised, vireoSNP: \n
        --sample {sample}
        --min_cov_treshold {min_cov_treshold}
        --p_treshold {p_treshold}
        """
    )
    
    ##
    
    # Read data
    t.start()
    afm = read_one_sample(path_data, sample=sample, with_GBC=False)

    # Filter cells and vars, create a mutational embedding
    _path = os.path.join(path_variants, sample, 'MQuad_vars.csv')
    if os.path.exists(_path):
        vois = pd.read_csv(_path, index_col=0)['0']
    else:
        raise ValueError(f'Cannot find vois dataframe at {_path}')
    
    logger.info(f'Found {vois.size} variants for sampe {sample}.')
    
    _, a = filter_cells_and_vars(
        afm,
        sample=sample,
        min_cov_treshold=min_cov_treshold,
        variants=vois
    )

    # UMAP
    a = nans_as_zeros(a) # For sklearn APIs compatibility
    embs, _ = reduce_dimensions(a, method='UMAP', n_comps=2, sqrt=False)
    embs = pd.DataFrame(embs, index=a.obs_names, columns=['UMAP1', 'UMAP2'])
    embs.to_csv(os.path.join(path_results, sample, 'MT_SNVs_umap.csv'))

    logger.info(f'Data preparation and UMAP computation {t.stop()}')

    ##

    # VireoSNP clustering
    t.start()
    df, n_best, fig = vireo_clustering(a, with_knee_plot=True)
    logger.info(f'Finished vireoSNP clustering {t.stop()}')
    logger.info(f'Best solution yiels {n_best} MT_clones.')

    # Save
    df.to_csv(os.path.join(path_results, sample, 'vireo_clones.csv'))
    fig.savefig(
        os.path.join(path_results, sample, 'vireo_ELBO.csv'),
        dpi=300
    )

    # Exit
    logger.info(f'Execution was completed successfully in total {T.stop()}')

####################################################################

# Run program
if __name__ == "__main__":
    main()