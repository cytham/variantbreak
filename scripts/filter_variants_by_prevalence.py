import os
import pandas as pd
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import ScalarFormatter
import scipy.stats as stats
import fastcluster

filename = "data.csv"
df0 = pd.read_csv(filename, delimiter='\t', index_col = 0)
df = df0[~((df0['Filter3'] == 1) & (df0['End'] - df0['Start'] < 100))]  # Only remove SV in homopolyer if SV < 100 bp

sample_names = ['S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 'S9', 'S10', 'S11', 'S12', 'S13', 'S14', 'S15', 'S16',
                'S17', 'S18', 'S19', 'S20', 'S21', 'S22', 'S23', 'S24']

depth = [
    7.55,
    8.14,
    7.66,
    3.68,
    9.38,
    3.42,
    4.9,
    7.18,
    4.38,
    4.44,
    4.23,
    9.36,
    11.18,
    9.85,
    8.09,
    6.61,
    9.88,
    9.1,
    9.04,
    7.0,
    8.25,
    11.64,
    8.06,
    7.81
]

depth_dict = {}
for i,j in zip(sample_names, depth):
    depth_dict[i] = j

# Create list of only samples with depth > 6
sample_list = []
for i in sample_names:
    if depth_dict[i] > 6:
        sample_list.append(i)

# Create list of non-AML sample names
normal = ['S18', 'S19', 'S20', 'S21', 'S22', 'S23', 'S24']

# Create lists of AML sample names
aml = []
for name in sample_list:
    if name in normal:
        pass
    else:
        aml.append(name)

# Get non_AML array from dataframe
non_aml_arr = df.loc[:, normal].values

# Convert all SV hits to 1 or misses to 0
non_aml_arr[non_aml_arr != '-'] = 1
non_aml_arr[non_aml_arr == '-'] = 0

# Count SV hits and misses
non_aml_out = []
for i in non_aml_arr:
    non_aml_out.append([np.count_nonzero(i == 0), np.count_nonzero(i == 1)])

# Get AML array from dataframe
aml_arr = df.loc[:, aml].values

# Convert all SV hits to 1 or misses to 0
aml_arr[aml_arr != '-'] = 1
aml_arr[aml_arr == '-'] = 0

# Count SV hits and misses
aml_out = []
for i in aml_arr:
    aml_out.append([np.count_nonzero(i == 0), np.count_nonzero(i == 1)])

# Combine AML and non-AML and tranpose list, generating list of 2D arrays
out = []
for i in zip(aml_out, non_aml_out):
    out.append(list(map(list, zip(*i))))

total_aml = len(aml)
total_norm = len(normal)



# A) Detect AML-associated variants by prevalence [CHOOSE A OR B]
aml_thres = 0.8  # Change alpha values to adjust stringency
norm_thres = 0.2  # Change beta values to adjust stringency
a = 0
indexes = []
memory = []
for i in out:
    if i[1][0]/total_aml >= aml_thres and i[1][1]/total_norm <= norm_thres:
        indexes.append(a)
        memory.append([a, i[1][0]/total_aml, i[1][1]/total_norm])
    a += 1


# B) Detect normal-associated variants by prevalence [CHOOSE A OR B]
aml_thres = 0.2  # Change beta values to adjust stringency
norm_thres = 0.8 # Change alpha values to adjust stringency
a = 0
indexes_norm = []
for i in out:
    if i[1][0]/total_aml <= aml_thres and i[1][1]/total_norm >= norm_thres:
        indexes_norm.append(a)
        memory.append([a, i[1][0] / total_aml, i[1][1] / total_norm])
    a += 1



select_df = df.iloc[indexes + indexes_norm, :]

# Cluster using linkage hclustering
# Create numpy array for selected samples
arr = select_df.loc[:, aml + normal].values

# Change all SV hits to value of 1 and misses to 0
arr[arr != '-'] = 1
arr[arr == '-'] = 0

# Change astype to int
arr = arr.astype(int)
z = fastcluster.linkage(arr, 'ward')
n = len(z) + 1
cache = dict()
for k in range(len(z)):
    c1, c2 = int(z[k][0]), int(z[k][1])
    c1 = [c1] if c1 < n else cache.pop(c1)
    c2 = [c2] if c2 < n else cache.pop(c2)
    cache[n + k] = c1 + c2

index = cache[2 * len(z)]
sv_order = select_df.index.values.tolist()
indexes1 = []
for i in index:
    indexes1.append(sv_order[i])

# Get default metadata
with pd.HDFStore('data.h5', mode="r") as store:
    df_default = store['dataset']
    metadata_default = store.get_storer('dataset').attrs.metadata

new_df = df_default.reindex(index=indexes1)

# Save dataframe to CSV
new_df.to_csv('new.csv', index=True,
              sep='\t',
              index_label='Index')

# Save dataframe and metadata to HDF5 for viewing on VariantMap
store = pd.HDFStore('new.h5', mode='w')
store.put('dataset', new_df)
store.get_storer('dataset').attrs.metadata = metadata_default
store.close()


