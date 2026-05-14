# ------------------------------------------------------------------------------
# Imports
# ------------------------------------------------------------------------------

import argparse
from SDP_interaction_inference.constraints import Constraint
from SDP_interaction_inference.optimization import ModelFreeOptimization
from SDP_interaction_inference.dataset import Dataset
import pandas as pd
import numpy as np
import json


# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------

def construct_dataset(mirna_sample, mrna_dataset, beta, resamples=1000):

    # size
    gene_pairs, cells = mrna_dataset.shape

    # construct paired count dataframe
    counts_df = pd.DataFrame(
        index = [f"Gene-pair-{i}" for i in range(gene_pairs)],
        columns = [f"Cell-{j}" for j in range(cells)]
    )

    # fill with pairs
    for i in range(gene_pairs):
        gene_i = mirna_sample
        gene_j = mrna_dataset.iloc[i]
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
parser.add_argument("--miRNA")
parser.add_argument("--mRNA", nargs="*")
parser.add_argument("--phase", choices=["All", "G1", "S", "G2M"])

# bootstrap arguments
parser.add_argument("--confidence", default=0.95)
parser.add_argument("--resamples", default=1000)

# optimization arguments
parser.add_argument("--d", default=3)
parser.add_argument("--time_limit", default=30)
parser.add_argument("--total_time_limit", default=30)
parser.add_argument("--cut_limit", default=100)

# result arguments
parser.add_argument("--results", nargs="*", choices=["status", "time", "cuts"], default=["status"])
parser.add_argument("--outfile", required=True)

# parse
args = parser.parse_args()

# ------------------------------------------------------------------------------
# Data processing
# ------------------------------------------------------------------------------

# load counts
data = pd.read_csv("../../Moment-equations/Real-Data-2/Data/GSE151334_FIB_counts_thresh.csv", index_col=0)

# load capture
beta = np.loadtxt("../../Moment-equations/Real-Data-2/Capture/beta_FIB.txt")

# load cell cycle phases
phases = pd.read_csv("../R-Notebooks/Cell-Cycle/phase.csv", index_col=0)

# reduce to FIB cells
phases = phases.loc[data.columns]

# load RNA types
biotypes_dict = json.load(open("../../Moment-equations/Real-Data-2/Biotypes/biotypes_FIB.json"))

# select indices of protein coding mRNA and non-coding miRNA
pcRNA_indices = [idx for idx, btype in enumerate(biotypes_dict.values()) if btype == "protein_coding"]
miRNA_indices = [idx for idx, btype in enumerate(biotypes_dict.values()) if btype == "miRNA"]

# select all cells
if args.phase == "All":

    # split into mRNA & miRNA
    data_pcRNA = data.iloc[pcRNA_indices]
    data_miRNA = data.iloc[miRNA_indices]

else:

    # select cells of phase
    phase_mask = (phases['Phase_per_type'] == args.phase)
    data = data.loc[:, phase_mask]
    beta = beta[phase_mask]

    # split into mRNA & miRNA
    data_pcRNA = data.iloc[pcRNA_indices]
    data_miRNA = data.iloc[miRNA_indices]

# ------------------------------------------------------------------------------
# Running
# ------------------------------------------------------------------------------

# dataframe to store model free results
result_df = pd.DataFrame(
    index=args.mRNA
)

# construct dataset of miRNA paired with mRNA
dataset_SDP = construct_dataset(
    data_miRNA.loc[args.miRNA],
    data_pcRNA.loc[args.mRNA],
    beta,
    resamples=args.resamples
)

# bootstrap
dataset_SDP.confidence = args.confidence
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
    result_df[f'{args.miRNA}_{args.phase}_d{args.d}_c{int(args.confidence * 100)}_{res}'] = res_list

# write results
result_df.to_csv(args.outfile)
