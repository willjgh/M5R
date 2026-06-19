# ------------------------------------------------------------------------------
# Imports
# ------------------------------------------------------------------------------

import argparse
from SDP_interaction_inference.constraints import Constraint
from SDP_interaction_inference.optimization import BirthDeathOptimization
from SDP_interaction_inference.dataset import SparseDataset
import pandas as pd
import numpy as np
import anndata as ad
import json

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
parser.add_argument("--fixed", nargs="*", default=[(3, 1)])
parser.add_argument("--K", default=None)
parser.add_argument("--time_limit", default=30, type=int)
parser.add_argument("--total_time_limit", default=30, type=int)
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

# load independent results
ind_df = pd.read_csv(f"Results/ind_MF_{args.array_index}.csv", index_col=0)

# select infeasible pcRNA
pcRNA_names = ind_df[ind_df.iloc[:, 0] == "INFEASIBLE"].index.tolist()

# convert to adata indices
mask_ad = adata_pcRNA.var['GeneName'].isin(pcRNA_names)
pcRNA_idxs = mask_ad.index[mask_ad]

# names
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
    adata_pcRNA[:, pcRNA_idxs],
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
    moment_equations=True,
    factorization=False
)

# construct optimizer
opt_BD = BirthDeathOptimization(
    dataset_SDP,
    d_bd=args.d,
    d_me=args.d,
    d_sd=args.d,
    constraints=constraints,
    printing=False,
    silent=True,
    time_limit=args.time_limit,
    total_time_limit=args.total_time_limit,
    cut_limit=args.cut_limit,
    fixed=args.fixed,
    K=args.K
)

# optimize
opt_BD.analyse_dataset()

# for each desired result
for res in args.results:

    # extract results
    res_list = [solution[res] for solution in opt_BD.result_dict.values()]

    # store
    result_df[f'{miRNA_name}_d{args.d}_c{int(args.confidence * 100)}_{res}'] = res_list

# write results
result_df.to_csv(f"Results/BD_{args.array_index}.csv")
