"""
Simple assessment of coverage differences across MT-genome, with and without masking of NUMTs regions on nuclear chromosomes.
"""

import os
import pandas as pd
import numpy as np
from mito_utils.plotting_base import *
matplotlib.use('macOSX')


##


# Paths
path_main = '/Users/IEO5505/Desktop/AML_clonal_reconstruction'
path_data = os.path.join(path_main, 'data', 'other')
os.listdir(path_data)

# Read coverages
wi = pd.read_csv(os.path.join(path_data, 'cov_NUMTS.csv'), index_col=0)['0']
wo = pd.read_csv(os.path.join(path_data, 'cov_no_NUMTS.csv'), index_col=0)['0']
wi.describe()
wo.describe()

# Viz
fig, axs = plt.subplots(1,2,figsize=(10,5))
axs[0].plot(wi.values, 'k-')
axs[0].set(title='Still NUMTs')
axs[1].plot(wo.values, 'r-')
axs[1].set(title='Masked NUMTs')
fig.tight_layout()
plt.show()


##