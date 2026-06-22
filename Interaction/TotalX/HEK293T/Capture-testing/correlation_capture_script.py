# ------------------------------------------------------------------------------
# Imports
# ------------------------------------------------------------------------------

import argparse
from SDP_interaction_inference.dataset import SparseDataset
import pandas as pd
import numpy as np
import anndata as ad
import json
import tqdm
import scipy

# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------

def correlation_bootstrap_sample(rng, sample, beta, confidence=None, resamples=None):

    # get sample size
    n = sample.shape[0]

    # get bootstrap size: default to 1000
    if resamples is None:
        resamples = 1000
    # confidence level: default to 95%
    if confidence is None:
        confidence = 0.95

    # bootstrap to N x n x 2 array
    boot = rng.choice(sample, size=(resamples, n))

    # capture moments
    E_beta = np.mean(beta)
    E_beta2 = np.mean(beta**2)

    # compute correlations
    estimates = np.zeros(resamples)

    b1 = boot[:, :, 0]
    b2 = boot[:, :, 1]

    # OB moments
    E_xy_OB = np.mean(b1 * b2, axis=1)
    E_x_OB = np.mean(b1, axis=1)
    E_y_OB = np.mean(b2, axis=1)
    E_x2_OB = np.mean(b1**2, axis=1)
    E_y2_OB = np.mean(b2**2, axis=1)

    # OG moments
    E_xy_OG = E_xy_OB / E_beta2
    E_x_OG = E_x_OB / E_beta
    E_y_OG = E_y_OB / E_beta
    E_x2_OG = (1 / E_beta2)*E_x2_OB + (1 / E_beta)*E_x_OB - (1 / E_beta2)*E_x_OB
    E_y2_OG = (1 / E_beta2)*E_y2_OB + (1 / E_beta)*E_y_OB - (1 / E_beta2)*E_y_OB

    varx_OG = E_x2_OG - E_x_OG**2
    vary_OG = E_y2_OG - E_y_OG**2

    mask = (varx_OG > 0.0) & (vary_OG > 0.0)
    estimates[~mask] = np.nan
    estimates[mask] = (E_xy_OG[mask] - E_x_OG[mask]*E_y_OG[mask]) / (np.sqrt(varx_OG[mask]) * np.sqrt(vary_OG[mask]))

    # take quantiles
    alpha = 1 - confidence
    interval = np.nanquantile(estimates, [(alpha / 2), 1 - (alpha / 2)])

    # compute point estimate from original sample
    b1 = sample[:, 0]
    b2 = sample[:, 1]

    # OB moments
    E_xy_OB = np.mean(b1 * b2)
    E_x_OB = np.mean(b1)
    E_y_OB = np.mean(b2)
    E_x2_OB = np.mean(b1**2)
    E_y2_OB = np.mean(b2**2)

    # OG moments
    E_xy_OG = E_xy_OB / E_beta2
    E_x_OG = E_x_OB / E_beta
    E_y_OG = E_y_OB / E_beta
    E_x2_OG = (1 / E_beta2)*E_x2_OB + (1 / E_beta)*E_x_OB - (1 / E_beta2)*E_x_OB
    E_y2_OG = (1 / E_beta2)*E_y2_OB + (1 / E_beta)*E_y_OB - (1 / E_beta2)*E_y_OB

    varx_OG = E_x2_OG - E_x_OG**2
    vary_OG = E_y2_OG - E_y_OG**2

    if varx_OG <= 0.0 or vary_OG <= 0.0:
        corr = np.nan
    else:
        corr = (E_xy_OG - E_x_OG*E_y_OG) / (np.sqrt(varx_OG) * np.sqrt(vary_OG))

    # collect results
    result = np.array([
        corr,
        interval[0],
        interval[1]
    ])

    return result

def correlation_bootstrap_dataset(rng, dataset, confidence=None, resamples=None):
    '''Bootstrap dataset correlations.'''

    # convert to sparse column array for faster slicing
    Xmi = dataset.miRNA_dataset.tocsc()
    Xpc = dataset.pcRNA_dataset.tocsc()

    # record results
    results = np.empty((dataset.gene_pairs, 3))
    
    # bootstrap
    for i in tqdm.tqdm(range(dataset.gene_pairs)):

        # select pair
        pair = dataset.pair_indices[i]

        # select samples
        miRNA_sample = Xmi[:, pair[0]]
        pcRNA_sample = Xpc[:, pair[1]]

        # combine
        sample = scipy.sparse.hstack([miRNA_sample, pcRNA_sample])
        sample = sample.toarray()
    
        # compute
        results[i, :] = correlation_bootstrap_sample(rng, sample, np.ones(dataset.cells), confidence, resamples)

    return results

def analytic_correlation_bootstrap_dataset(rng, dataset, confidence=None, resamples=None):
    '''Bootstrap dataset correlations adjusted for capture efficiency.'''

    # convert to sparse column array for faster slicing
    Xmi = dataset.miRNA_dataset.tocsc()
    Xpc = dataset.pcRNA_dataset.tocsc()

    # record results
    results = np.empty((dataset.gene_pairs, 3))
    
    # bootstrap
    for i in tqdm.tqdm(range(dataset.gene_pairs)):

        # select pair
        pair = dataset.pair_indices[i]

        # select samples
        miRNA_sample = Xmi[:, pair[0]]
        pcRNA_sample = Xpc[:, pair[1]]

        # combine
        sample = scipy.sparse.hstack([miRNA_sample, pcRNA_sample])
        sample = sample.toarray()

        # compute
        results[i, :] = correlation_bootstrap_sample(rng, sample, dataset.beta, confidence, resamples)

    return results

# ------------------------------------------------------------------------------
# Arguments
# ------------------------------------------------------------------------------

# script arguments
parser = argparse.ArgumentParser()

# dataset arguments
parser.add_argument("--array_index", type=int)
parser.add_argument("--test_limit", default=-1, type=int)
#parset.add_argument("--splits", default=1, type=int)

# bootstrap arguments
parser.add_argument("--confidence", default=0.95, type=float)
parser.add_argument("--resamples", default=1000, type=int)

# result arguments
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

# test limit on pcRNA
if args.test_limit == -1:
    G = None
else:
    G = args.test_limit

# names
pcRNA_names = adata_pcRNA.var['GeneName'].tolist()[:G]
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
    adata_pcRNA[:, :G],
    beta
)

# initialize random generator
rng = np.random.default_rng()

# bootstrap observed correlation
observed_correlation = correlation_bootstrap_dataset(rng, dataset_SDP, confidence=args.confidence, resamples=args.resamples)

# bootstrap analytic recovered correlation
analytic_correlation = analytic_correlation_bootstrap_dataset(rng, dataset_SDP, confidence=args.confidence, resamples=args.resamples)

# store
result_df[f'{miRNA_name}_c{int(100*args.confidence)}_OB'] = observed_correlation[:, 0]
result_df[f'{miRNA_name}_c{int(100*args.confidence)}_OB_min'] = observed_correlation[:, 1]
result_df[f'{miRNA_name}_c{int(100*args.confidence)}_OB_max'] = observed_correlation[:, 2]
result_df[f'{miRNA_name}_c{int(100*args.confidence)}_AL'] = analytic_correlation[:, 0]
result_df[f'{miRNA_name}_c{int(100*args.confidence)}_AL_min'] = analytic_correlation[:, 1]
result_df[f'{miRNA_name}_c{int(100*args.confidence)}_AL_max'] = analytic_correlation[:, 2]

# write results
result_df.to_csv(f"Results-capture/corr_capture_10_{args.array_index}.csv")
