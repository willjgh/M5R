# M5R

Extension of M4R work investigating the application of optimization inference methods to different areas of genomics

# Interaction

## Overview

Genes do not act independently and can interact with each other to regulate transcription behaviour, even forming large and complex gene regulatory networks. We investigate the problem of detecting if two genes interact using samples of transcript counts from each.

## Model Experiments

Investigate stochastic reaction network models of gene expression and the performance of optimization inference

### birth death interaction

We first consider the reduced problem of detecting interaction in data simulated from a fixed reaction network:

$$ \varnothing \stackrel{k_{tx, 1}}\longrightarrow X_{1} \qquad \varnothing \stackrel{k_{tx, 2}}\longrightarrow X_{2} $$
$$ X_{1} \stackrel{k_{deg, 1}}\longrightarrow \varnothing \qquad X_{2} \stackrel{k_{deg, 2}}\longrightarrow \varnothing $$
$$ X_{1} + X_{2} \stackrel{k_{reg}}\longrightarrow \varnothing $$

where the genes follow a birth-death model of transcription and interact by regulating each others transcription. The problem of detecting interaction then reduces to determining if the parameter $k_{reg}$ is non-zero: using our optimization based inference methods we can produce bounds on the parameter or test infeasibility under this assumption.

Investigate:
- method 1 (min): optimize for min of $k_{reg}$, concluding interaction is present if non-zero
- method 2 (hyp): assuming no interaction / $k_{reg}$ = 0, test feasiblity of data under the reduced model
- compare results to simple correlation tests (corr) for a range of interaction strengths and capture efficiency (measurement error) [Results data]
- new approaches to truncation of distributions

### telegraph interaction

Consider the more general telegraph model structure for the assumed reaction network:

$$ G_{1, off} \underset{k_{1, off}}{\stackrel{k_{1, on}}{\rightleftharpoons}} G_{1, on}  \quad\quad G_{2, off} \underset{k_{2, off}}{\stackrel{k_{2, on}}{\rightleftharpoons}} G_{2, on} $$
$$ G_{1, on} \stackrel{k_{1, tx}}\longrightarrow G_{1, on} + X_{1} \quad\quad G_{2, on} \stackrel{k_{2, tx}}\longrightarrow G_{2, on} + X_{2} $$
$$ X_{1} \stackrel{k_{1, deg}}\longrightarrow \varnothing \quad\quad X_{2} \stackrel{k_{2, deg}}\longrightarrow \varnothing $$
$$ X_{1} + X_{2} \stackrel{k_{reg}}\longrightarrow \varnothing $$

This is a more general model that can capture a wider variety of transcriptional behaviours, at a cost of more computationally expensive and difficult inference.

Investigate:
- performance of min, hyp and corr methods
- GUROBI optimization settings to reduce computation time

TODO:
- revisit using improved hyp methods

### cross model interaction

So far we have considered data simulated from a known model and optimized using this known structure, but with real data we do not know 'which model the data follows' as in reality the models are only approximations used to describe the process of gene expression. So which model should be used? We investigate performance of methods when applied to data simulated from different models e.g. 'crossing' telegraph data with birth death optimization.

Investigate:
- cross model performnce in different parameter settings (telegraph data close to birth-death / very different) [cross model results]
- effect of confidence interval level across interaction strengths (pseudeo p-value for detection conclusion) [strength confidence results]

TODO:
- revisit using improved hyp methods
- ROC / AUC classification metrics using confidence level data

## Marginal Experiments

Investigate the use of marginal probabilities and the performance of the Bootstrap for non-parametric estimation of confidence intervals

### bootstrap investigation

Compare pairwise to individual bootstrap resampling for marginal probability confidence intervals and produce memory efficient code for large sample and bootstrap sizes e.g. 100,000+

### marginal optimization

Compare hypothesis method using only marginal information to full hypothesis method and correlation tests for increasing sample sizes. Suggests marginal results converge to full results with increasing sample sizes but lags behind correlation performance.

### perfect information optimization

Remove bootstrap estimation uncertainity by using exact probabilites / extreme sample sizes and confirm that optimization can produce exact bounds on model parameters in this 'perfect' setting

## Realistic Data

Investigate method performance on data simulated with realistic model parameters (informed by observed data) to get an idea of true performance.

### Simulation

Simulation of realistic model parameters and corresonding samples of transcript counts, producing gene x cell datasets of counts for analysis. 

### Dataset Analysis

Full analysis pipeline for dataset of single cell data, given optimization method and settings, and code to analyse results, producing classification metrics and scatter plots of results.

### Datasets 'Easy' 'Hard' + Results 'Easy' 'Hard'

Simulated data and corresponding results for 'Easy' and 'Hard' datasets i.e. cases where interaction is relatively easy to detect from samples vs cases where interaction is very difficult to detect

### Testing

Testing the performance of optimization methods on these more challenging parameter settings, making changes to heuristics used for data manipulation and optimization.

Investigate:
- classification metrics, confusion matrices, scatter plots and other illustrations of results
- results produced using basic setup
- alternate truncation methods for the range of observed / original counts considered during bootstrap / optimization
- illustration of truncation ranges

## Truncations

Code to compute truncation limits for use in optimization:

[image]

# Perturbation

## Overview

Experiments often 'perturbe' cells, making changes to the environment such as temperature, applying drugs, etc which can affect how genes are expressed. Using data from 'original' and 'perturbed' sets of cells we want to identify if and how the expression behaviour of different genes changes.

Differential gene expression analysis methods such as the package 'DEseq2' can be use to identify perturbations, producing volcano plots that show the statistical significance (p-value) and magnitude (log fold change) of changes in expression.

We invesitgate the use of optimization based inference to identify perturbations: using a telegraph model of expression we can use non-linear optimization to estimate bounds on model parameters using single-cell data. Applying the method to 'original' and 'perturbed' data for the same gene we compare estimated bounds and conclude if there are statistically significant changes:

![screenshot](Perturbation/Plots/example_results_plot.png)

In this example the estimated bounds for all parameters overlap suggesting there is not a stat. sig. change in expression behaviour of the gene, a reasonable conclusion given that the true perturbation effect was a slight decrease of the gene's 'on rate'.

## perturbation simulation

To test the performance of methods we simulate syntheic perturbation data. Includes code to produce 'original' and 'perturbed' datasets of transcript counts and parameter values.

For each gene in a dataset:
- sample original model parameters from log-uniform
- sample original transcript counts from telegraph model stationary distribution for each cell (can also sample and apply capture efficiency)
- simulate perturbation effects to produce perturbed model parameters
- sample perturbed transcript counts

## perturbation identification

Code implementing the optimization method to identify perturbation.

For each gene:
- compute bootstrap confidence intervals on distribution of 'original' and 'perturbed' counts using samples
- optimize for each to estimate bounds on model parameters
- disjoint intervals suggest stat. sig. change in parameter value (sig. level according to confidence of bootstrap intervals)

## DESeq2

R script to produce results of DESeq2 perturbation identification to compare

## Plots

Example plots of results from optimization and differential analysis methods

## Simulated data

Example of datasets simulated using 'perturbation_simulation.ipynb'
