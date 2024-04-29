"""
Collapse pipeline output.
"""

# Code
import os
import re
import pandas as pd


##


## Helper
def get_files_l(path, pattern=None, print_all=False, from_tokey=None):

    path = path_data
    L = []
    d = {}
    for root, _, files in os.walk(path):
        for file in files:
            if print_all:
                print(os.path.join(root, file))
            else:
                if bool(re.search(pattern, file)):
                    if from_tokey is not None:
                        key = root.split('/')[-from_tokey:]
                        if len(key)>1:
                            key = tuple(key)
                        else:
                            key = key[0]
                        d[key] = os.path.join(root, file)
                    else:
                        L.append(os.path.join(root, file))

    if from_tokey is not None:
        return d
    else:          
        return L


# Paths
path_data = '/Users/IEO5505/Desktop/AML_clonal_reconstruction/results/trees/output'


##


# Collapse outputs in summary dfs
get_files_l(path_data, print_all=True)

# Tree outputs
df = pd.concat([ 
    pd.read_csv(x, index_col=0) 
    for x in get_files_l(os.path.join(path_data), pattern='supp')
])
df.to_csv(os.path.join(path_data, 'supports_df.csv'))
df = pd.concat([ 
    pd.read_csv(x, index_col=0) 
    for x in get_files_l(os.path.join(path_data), pattern='resul')
])
df.to_csv(os.path.join(path_data, 'filtering_df.csv'))

# Muts outputs
# df = pd.concat([ 
#     pd.read_csv(v, index_col=0).assign(
#         sample=k[0], filtering=k[1], solver=k[2], metric=k[3]
#     )
#     for k, v in get_files_l(os.path.join(path_data), pattern='mut', from_tokey=4).items()
# ])
# df.to_csv(os.path.join(path_data, 'muts_df.csv'))

# Clones outputs
# df = pd.concat([ 
#     pd.read_csv(v, index_col=0).assign(
#         sample=k[0], filtering=k[1], solver=k[2], metric=k[3]
#     )
#     for k, v in get_files_l(os.path.join(path_data), pattern='clone', from_tokey=4).items()
# ])
# df.to_csv(os.path.join(path_data, 'clones_df.csv'))


##

