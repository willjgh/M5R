# ------------------------------------------------------------------------------
# Imports
# ------------------------------------------------------------------------------

import argparse
from SDP_interaction_inference.optimization_MOSEK import MOSEKModelFreeInteracting
from SDP_interaction_inference.dataset import Dataset
import pandas as pd
import numpy as np
import json
import tqdm

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

def correlation_bootstrap_sample(rng, sample, beta, confidence=None, resamples=None):

    # get sample size
    n = sample.shape[0]

    # get bootstrap size: default to sample size
    if resamples is None:
        resamples = n
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

    # record results
    results = np.empty((dataset.gene_pairs, 3))
    
    # bootstrap
    for i in tqdm.tqdm(range(dataset.gene_pairs)):
        results[i, :] = correlation_bootstrap_sample(rng, np.array([*dataset.count_dataset.iloc[i].values]), np.ones(dataset.cells), confidence, resamples)

    return results

def analytic_correlation_bootstrap_dataset(rng, dataset, confidence=None, resamples=None):
    '''Bootstrap dataset correlations adjusted for capture efficiency.'''

    # record results
    results = np.empty((dataset.gene_pairs, 3))
    
    # bootstrap
    for i in tqdm.tqdm(range(dataset.gene_pairs)):
        results[i, :] = correlation_bootstrap_sample(rng, np.array([*dataset.count_dataset.iloc[i].values]), dataset.beta, confidence, resamples)

    return results

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

# result arguments
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

# initialize random generator
rng = np.random.default_rng()

# bootstrap observed correlation
observed_correlation = correlation_bootstrap_dataset(rng, dataset_SDP, confidence=args.confidence, resamples=args.resamples)

# bootstrap analytic recovered correlation
analytic_correlation = analytic_correlation_bootstrap_dataset(rng, dataset_SDP, confidence=args.confidence, resamples=args.resamples)

# store
result_df[f'{args.miRNA}_{args.phase}_c{int(100*args.confidence)}_OB'] = observed_correlation[:, 0]
result_df[f'{args.miRNA}_{args.phase}_c{int(100*args.confidence)}_OB_min'] = observed_correlation[:, 1]
result_df[f'{args.miRNA}_{args.phase}_c{int(100*args.confidence)}_OB_max'] = observed_correlation[:, 2]
result_df[f'{args.miRNA}_{args.phase}_c{int(100*args.confidence)}_AL'] = analytic_correlation[:, 0]
result_df[f'{args.miRNA}_{args.phase}_c{int(100*args.confidence)}_AL_min'] = analytic_correlation[:, 1]
result_df[f'{args.miRNA}_{args.phase}_c{int(100*args.confidence)}_AL_max'] = analytic_correlation[:, 2]

# write results
result_df.to_csv(args.outfile)
