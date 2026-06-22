# ------------------------------------------------------------------------------
# Imports
# ------------------------------------------------------------------------------

import argparse
from SDP_interaction_inference.constraints import Constraint
from SDP_interaction_inference.optimization_MOSEK import MOSEKModelFreeInteracting
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
parser.add_argument("--N", default=1000, type=int)
parser.add_argument("--time_limit", default=30, type=float)

# result arguments
parser.add_argument("--method_confidence", nargs="*", default=["95"])
#parser.add_argument("--outfile", required=True)

# parse
args = parser.parse_args()

# ------------------------------------------------------------------------------
# Data processing
# ------------------------------------------------------------------------------

# load pcRNA
adata_pcRNA = ad.read_h5ad("../TotalX_HEK293T_pcRNA.h5ad")

# load miRNA
adata_miRNA = ad.read_h5ad("../TotalX_HEK293T_miRNA.h5ad")

# load capture
beta = np.loadtxt("TotalX_HEK293T_capture_10.txt")

# load independent results
ind_df = pd.read_csv(f"Results-capture/ind_MF_capture_10_{args.array_index}.csv", index_col=0)

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

# construct optimizer
opt_har = MOSEKModelFreeInteracting(
    dataset_SDP,
    d=args.d,
    N=args.N,
    time_limit=args.time_limit
)

# optimize
opt_har.analyse_dataset()

# collect results
har_correlations_dataset = [solution['feasible_correlations'] for solution in opt_har.result_dict_HAR.values()]
har_fano_1_dataset = [solution['feasible_fano_factors_1'] for solution in opt_har.result_dict_HAR.values()]
har_fano_2_dataset = [solution['feasible_fano_factors_2'] for solution in opt_har.result_dict_HAR.values()]

# for each desired confidence
for mc in args.method_confidence:

    # convert to float
    mc = float(mc) / 100

    # set significance
    alpha_method = 1 - mc

    # extract correlations
    min_correlations = []
    max_correlations = []

    for har_correlations in har_correlations_dataset:

        # exception leads to no results
        if har_correlations == None:

            min_correlations.append(np.nan)
            max_correlations.append(np.nan)

        else:

            # tranform None to nan for quantiles
            har_correlations = np.array(har_correlations, dtype=float)
            interval = np.nanquantile(har_correlations, [(alpha_method / 2), 1 - (alpha_method / 2)])

            min_correlations.append(interval[0])
            max_correlations.append(interval[1])

    # extract fano factors
    min_fano_1 = []
    min_fano_2 = []
    max_fano_1 = []
    max_fano_2 = []

    for har_fano_1 in har_fano_1_dataset:

        # exception leads to no results
        if har_fano_1 == None:

            min_fano_1.append(np.nan)
            max_fano_1.append(np.nan)

        else:

            # tranform None to nan for quantiles
            har_fano_1 = np.array(har_fano_1, dtype=float)
            interval_1 = np.nanquantile(har_fano_1, [(alpha_method / 2), 1 - (alpha_method / 2)])

            min_fano_1.append(interval_1[0])
            max_fano_1.append(interval_1[1])

    for har_fano_2 in har_fano_2_dataset:

        # exception leads to no results
        if har_fano_2 == None:

            min_fano_2.append(np.nan)
            max_fano_2.append(np.nan)

        else:

            # tranform None to nan for quantiles
            har_fano_2 = np.array(har_fano_2, dtype=float)
            interval_2 = np.nanquantile(har_fano_2, [(alpha_method / 2), 1 - (alpha_method / 2)])

            min_fano_2.append(interval_2[0])
            max_fano_2.append(interval_2[1])

    # store
    result_df[f'{miRNA_name}_d{args.d}_N{args.N}_c{int(args.confidence * 100)}_mc{int(100 * mc)}_HAR_corr_min'] = min_correlations
    result_df[f'{miRNA_name}_d{args.d}_N{args.N}_c{int(args.confidence * 100)}_mc{int(100 * mc)}_HAR_corr_max'] = max_correlations

    result_df[f'{miRNA_name}_d{args.d}_N{args.N}_c{int(args.confidence * 100)}_mc{int(100 * mc)}_HAR_fano_1_min'] = min_fano_1
    result_df[f'{miRNA_name}_d{args.d}_N{args.N}_c{int(args.confidence * 100)}_mc{int(100 * mc)}_HAR_fano_1_max'] = max_fano_1

    result_df[f'{miRNA_name}_d{args.d}_N{args.N}_c{int(args.confidence * 100)}_mc{int(100 * mc)}_HAR_fano_2_min'] = min_fano_2
    result_df[f'{miRNA_name}_d{args.d}_N{args.N}_c{int(args.confidence * 100)}_mc{int(100 * mc)}_HAR_fano_2_max'] = max_fano_2

# write results
result_df.to_csv(f"Results-capture/int_MF_capture_10_{args.array_index}.csv")
