# ------------------------------------------------------------------------------
# Imports
# ------------------------------------------------------------------------------

import argparse
from SDP_interaction_inference.constraints import Constraint
from SDP_interaction_inference.optimization import ModelFreeOptimization
from SDP_interaction_inference.dataset import SparseDataset
import pandas as pd
import numpy as np
import anndata as ad
import json

# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------

def construct_dataset(mirna_sample, mrna_dataset, beta, resamples):

    # size
    cells, gene_pairs = mrna_dataset.shape

    # construct paired count dataframe
    counts_df = pd.DataFrame(
        index = [f"Gene-pair-{i}" for i in range(gene_pairs)],
        columns = [f"Cell-{j}" for j in range(cells)]
    )

    # fill with pairs
    for i in range(gene_pairs):
        gene_i = mirna_sample[:, 0]
        gene_j = mrna_dataset[:, i]
        gene_pair_ij = list(zip(gene_i, gene_j))
        counts_df.iloc[i] = gene_pair_ij

    # construct dataset object
    data = Dataset()
    data.count_dataset = counts_df
    data.cells = cells
    data.gene_pairs = gene_pairs

    # settings
    data.resamples = resamples

    # set capture
    data.beta = beta

    return data

# ------------------------------------------------------------------------------
# Arguments
# ------------------------------------------------------------------------------

# script arguments
parser = argparse.ArgumentParser()

# dataset arguments
parser.add_argument("--array_index", type=int)
#parset.add_argument("--splits", default=1, type=int)

# bootstrap arguments
parser.add_argument("--confidence", default=0.95, type=float)
parser.add_argument("--resamples", default=1000, type=int)

# optimization arguments
parser.add_argument("--d", default=3, type=int)
parser.add_argument("--time_limit", default=30, type=float)
parser.add_argument("--total_time_limit", default=30, type=float)
parser.add_argument("--cut_limit", default=100, type=int)

# result arguments
parser.add_argument("--results", nargs="*", choices=["status", "time", "cuts"], default=["status"])
#parser.add_argument("--outfile", required=True)

# parse
args = parser.parse_args()

# ------------------------------------------------------------------------------
# Data processing
# ------------------------------------------------------------------------------

# load pcRNA
adata_pcRNA = ad.read_h5ad("TotalX_HEK293T_pcRNA.h5ad")

# load miRNA
adata_miRNA = ad.read_h5ad("TotalX_HEK293T_miRNA.h5ad")

# load capture
beta = np.loadtxt("TotalX_HEK293T_capture.txt")

# names
pcRNA_names = adata_pcRNA.var['GeneName'].tolist()[:100] # restrict for testing
miRNA_name = adata_miRNA.var['GeneName'].iloc[args.array_index]

# ------------------------------------------------------------------------------
# Running
# ------------------------------------------------------------------------------

# dataframe to store model free results
result_df = pd.DataFrame(
    index=pcRNA_names
)

# setup dataset
dataset_SDP = SparseDataset()

# construct dataset of miRNA paired with mRNA
dataset_SDP.construct_dataset(
    adata_miRNA[:, args.array_index],
    adata_pcRNA[:, :100], # restrict for testing
    beta
)

# bootstrap
dataset_SDP.bootstrap(
    d=args.d,
    tqdm_disable=False
)

# constraints
constraints = Constraint(
    moment_bounds=True,
    moment_matrices=True,
    factorization=True
)

# construct optimizer
opt_MF_ind = ModelFreeOptimization(
    dataset_SDP,
    d_bd=args.d,
    d_me=0,
    d_sd=args.d,
    constraints=constraints,
    printing=False,
    silent=True,
    time_limit=args.time_limit,
    total_time_limit=args.total_time_limit,
    cut_limit=args.cut_limit
)

# optimize
opt_MF_ind.analyse_dataset()

# for each desired result
for res in args.results:

    # extract results
    res_list = [solution[res] for solution in opt_MF_ind.result_dict.values()]

    # store
    result_df[f'{miRNA_name}_d{args.d}_c{int(args.confidence * 100)}_{res}'] = res_list

# write results
result_df.to_csv(f"ind_MF_{args.array_index}.csv")
