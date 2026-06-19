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

def fano_bootstrap_sample(rng, sample, beta, confidence=None, resamples=None):

    # get sample size
    n = sample.shape[0]

    # get bootstrap size: default to 1000
    if resamples is None:
        resamples = 1000
    # confidence level: default to 95%
    if confidence is None:
        confidence = 0.95

    # bootstrap to N x n array
    boot = rng.choice(sample, size=(resamples, n))

    # capture moments
    E_beta = np.mean(beta)
    E_beta2 = np.mean(beta**2)

    # compute fano
    estimates = np.zeros(resamples)

    # OB moments
    E_x_OB = np.mean(boot, axis=1)
    E_x2_OB = np.mean(boot**2, axis=1)

    # OG moments
    E_x_OG = E_x_OB / E_beta
    E_x2_OG = (1 / E_beta2)*E_x2_OB + (1 / E_beta)*E_x_OB - (1 / E_beta2)*E_x_OB

    varx_OG = E_x2_OG - E_x_OG**2

    mask = (E_x_OG == 0.0) | (varx_OG < 0.0)
    estimates[mask] = np.nan
    estimates[~mask] = varx_OG[~mask] / E_x_OG[~mask]

    # take quantiles
    alpha = 1 - confidence
    interval = np.nanquantile(estimates, [(alpha / 2), 1 - (alpha / 2)])

    # compute point estimate from original sample

    # OB moments
    E_x_OB = np.mean(sample)
    E_x2_OB = np.mean(sample**2)

    # OG moments
    E_x_OG = E_x_OB / E_beta2
    E_x2_OG = (1 / E_beta2)*E_x2_OB + (1 / E_beta)*E_x_OB - (1 / E_beta2)*E_x_OB

    varx_OG = E_x2_OG - E_x_OG**2

    if E_x_OG == 0.0 or varx_OG < 0.0:
        fano = np.nan
    else:
        fano = varx_OG / E_x_OG

    # collect results
    result = np.array([
        fano,
        interval[0],
        interval[1]
    ])

    return result

def fano_bootstrap_dataset(rng, adata, confidence=None, resamples=None):
    '''Bootstrap data fano factors.'''

    # size
    C = adata.n_obs
    G = adata.n_vars

    # convert to sparse column array for faster slicing
    Xad = adata.X.tocsc()

    # record results
    results = np.empty((G, 3))
    
    # bootstrap
    for i in tqdm.tqdm(range(G)):

        # select sample
        sample = Xad[:, i].toarray().squeeze()
    
        # compute
        results[i, :] = fano_bootstrap_sample(rng, sample, np.ones(C), confidence, resamples)

    return results

def analytic_fano_bootstrap_dataset(rng, adata, beta, confidence=None, resamples=None):
    '''Bootstrap dataset fano factors adjusted for capture efficiency.'''

    # size
    G = adata.n_vars

    # convert to sparse column array for faster slicing
    Xad = adata.X.tocsc()

    # record results
    results = np.empty((G, 3))
    
    # bootstrap
    for i in tqdm.tqdm(range(G)):

        # select sample
        sample = Xad[:, i].toarray().squeeze()
    
        # compute
        results[i, :] = fano_bootstrap_sample(rng, sample, beta, confidence, resamples)

    return results

# ------------------------------------------------------------------------------
# Arguments
# ------------------------------------------------------------------------------

# script arguments
parser = argparse.ArgumentParser()

# dataset arguments
parser.add_argument("--RNA", choices=["pcRNA", "miRNA"])
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

if args.RNA == "pcRNA":

    # load pcRNA
    adata = ad.read_h5ad("TotalX_HEK293T_pcRNA.h5ad")

if args.RNA == "miRNA":

    # load miRNA
    adata = ad.read_h5ad("TotalX_HEK293T_miRNA.h5ad")

# load capture
beta = np.loadtxt("TotalX_HEK293T_capture.txt")

# test limit
if args.test_limit == -1:
    G = None
else:
    G = args.test_limit

adata = adata[:, :G]

# names
RNA_names = adata.var['GeneName'].tolist()[:G]

# ------------------------------------------------------------------------------
# Running
# ------------------------------------------------------------------------------

# dataframe to store model free results
result_df = pd.DataFrame(
    index=RNA_names
)

# initialize random generator
rng = np.random.default_rng()

# bootstrap observed fano factor
observed_fano = fano_bootstrap_dataset(rng, adata, confidence=args.confidence, resamples=args.resamples)

# bootstrap analytic recovered fano factor
analytic_fano = analytic_fano_bootstrap_dataset(rng, adata, beta, confidence=args.confidence, resamples=args.resamples)

# store
result_df[f'c{int(100*args.confidence)}_OB'] = observed_fano[:, 0]
result_df[f'c{int(100*args.confidence)}_OB_min'] = observed_fano[:, 1]
result_df[f'c{int(100*args.confidence)}_OB_max'] = observed_fano[:, 2]
result_df[f'c{int(100*args.confidence)}_AL'] = analytic_fano[:, 0]
result_df[f'c{int(100*args.confidence)}_AL_min'] = analytic_fano[:, 1]
result_df[f'c{int(100*args.confidence)}_AL_max'] = analytic_fano[:, 2]

# write results
result_df.to_csv(f"Results/fano_{args.RNA}.csv")
