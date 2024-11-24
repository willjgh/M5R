'''
Script to compute dictionary of original truncation bounds given a vector
of capture efficincies and threshold value

Settings:
- name of output file
- capture efficiency vector
- threshold for truncation
- number of values to compute

Output:
- json file of dictionary of truncations
'''

# ------------------------------------------------
# Dependencies
# ------------------------------------------------

import numpy as np
import scipy
import json
import tqdm

# ------------------------------------------------
# Settings
# ------------------------------------------------

output_filename = ".\Truncations\truncations_high"
M = 10
N = 10
beta = np.array([1.0 for j in range(1000)])
thresh_OG = 10**-6

# ------------------------------------------------
# Functions
# ------------------------------------------------

def B(x1_OB, x2_OB, x1_OG, x2_OG, beta):
    '''Compute (1 / n) sum j = 1 to n of P(X1_OB, X2_OB | X1_OG, X2_OG, Beta_j): product of binomials.'''

    if (x1_OB <= x1_OG) and (x2_OB <= x2_OG):

        val = 0.0

        bin_coef_x1 = scipy.special.comb(x1_OG, x1_OB)
        bin_coef_x2 = scipy.special.comb(x2_OG, x2_OB)

        for beta_j in beta:

            val += beta_j**(x1_OB + x2_OB) * (1 - beta_j)**(x1_OG - x1_OB + x2_OG - x2_OB)

        n = beta.shape[0]
        val = float(bin_coef_x1 * bin_coef_x2 * val / n)

        return val
    
    # asssume we cannot observe more counts than are present originally
    else:
        return 0.0

def findTrunc(x1_OB, x2_OB, beta, thresh_OG):
    '''
    Compute box truncation around states (x1_OG, x2_OG) which have
    B(x1_OB, x2_OB, x1_OG, x2_OG, beta) >= thresh_OG

    returns: m_OG, M_OG, n_OG, N_OG
    '''

    trunc_start = False
    trunc_end = False
    m_OG, M_OG, n_OG, N_OG = np.inf, 0, np.inf, 0
    diag = 0
    while (not trunc_start) or (not trunc_end):

        # start at top of grid
        x1_OG = x1_OB
        x2_OG = x2_OB + diag

        # flag if at least one coeff > thresh in diagonal
        trunc_diag = False

        # compute coeffs along diagonal
        while x2_OG >= x2_OB:

            # compute coeff
            coeff = B(x1_OB, x2_OB, x1_OG, x2_OG, beta)

            # above thresh
            if coeff >= thresh_OG:

                # update truncations
                if x1_OG < m_OG:
                    m_OG = x1_OG
                if x2_OG < n_OG:
                    n_OG = x2_OG
                if x1_OG > M_OG:
                    M_OG = x1_OG
                if x2_OG > N_OG:
                    N_OG = x2_OG

                # at least one coeff > thresh (overall)
                trunc_start = True

                # at least one coeff > thresh (in diag)
                trunc_diag = True

            # move down diagonal
            x2_OG -= 1
            x1_OG += 1

        # if NO coeff > thresh (in diag) AND at least one coeff > thresh (overall)
        if (not trunc_diag) and trunc_start:

            # end
            trunc_end = True

        # increment diagonal
        diag += 1

    return m_OG, M_OG, n_OG, N_OG

def preComputeTruncation(M, N, beta, thresh_OG):
    '''
    Compute M x N values of original truncations

    M, N: shape of observed pairs that truncations are computed for
    beta: capture efficiency vector
    thresh_OG: threshold for trunction
    '''
    # store in dictionary (lookup table)
    truncations = {}

    # for each pair of observed counts
    for x1_OB in tqdm.tqdm(range(M)):
        for x2_OB in range(N):

            # compute truncation bounds
            m_OG, M_OG, n_OG, N_OG = findTrunc(x1_OB, x2_OB, beta, thresh_OG)

            # store
            truncations[(x1_OB, x2_OB)] = (m_OG, M_OG, n_OG, N_OG)

    return truncations

# ------------------------------------------------
# Computation
# ------------------------------------------------

truncations = preComputeTruncation(M, N, beta, thresh_OG)

# ------------------------------------------------
# Output
# ------------------------------------------------

json.dump(solution_dict, open(output_filename, 'w'))
