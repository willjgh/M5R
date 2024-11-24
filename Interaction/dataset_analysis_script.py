'''
Script to analyse simulated datasets, using methods to detect interaction between gene-pairs from a stationary sample

Settings:
- file name of dataset of counts to be analysed
- file name of json file containing results
- threshold for bootstrap (confidence interval truncation)
- method used to analyse the dataset: min, hyp, corr
- capture efficiency vector
- pre-computed truncations (for same capture efficiency)
- other settings for chosen method (time limit, MIPGap, BestBdThresh, etc)

Output:
- Stores a dictionary (as json with given name) containing the results: interaction detected, time taken, etc for each gene-pair
'''

# ------------------------------------------------
# Dependencies
# ------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import gurobipy as gp
from gurobipy import GRB
import scipy
from ast import literal_eval
import json

# ------------------------------------------------
# Settings
# ------------------------------------------------

input_filename = ".\Datasets-Easy-Hard\"
output_filename = ".\Results-Easy-Hard\"
truncation_filename = ".\Truncations\"

thresh_OB = 10

beta = np.array([1.0 for j in range(1000)])
truncations = json.load(open(truncation_filename))

method = "hyp"

rng = np.random.default_rng(3844)

# ------------------------------------------------
# Functions
# ------------------------------------------------

def bootstrap(samples, BS=1000, thresh_OB=10, plot=False, printing=False):

    # sample size
    n = len(samples)

    # compute maximum x1 and x2 values
    M, N = np.max(samples, axis=0)
    M, N = int(M), int(N)

    # map (x1, x2) pairs to integers: x2 + (N + 1) * x1
    integer_samples = [x[1] + (N + 1)*x[0] for x in samples]

    # maxiumum of integer sample
    D = (M + 1)*(N + 1) - 1

    # B bootstrap samples as B x n array
    bootstrap = rng.choice(integer_samples, size=(BS, n))

    # offset row i by (D + 1)i
    bootstrap_offset = bootstrap + np.arange(BS)[:, None]*(D + 1)

    # flatten, count occurances of each state and reshape, reversing map to give counts of each (x1, x2) pair
    counts = np.bincount(bootstrap_offset.ravel(), minlength=BS*(D + 1)).reshape(-1, M + 1, N + 1)

    # counts to probability
    counts = counts / n

    # compute 2.5% and 97.5% quantiles for each p(x1, x2)
    bounds = np.quantile(counts, [0.025, 0.975], axis=0)

    # count occurances per (x1, x2) in the in original sample
    sample_counts = np.bincount(integer_samples, minlength=D + 1).reshape(M + 1, N + 1)

    # set truncation bounds
    m_OB, M_OB, n_OB, N_OB = M, 0, N, 0

    # set flag for changes
    thresh_flag = False

    # replace CI's for states below threshold occurances by [0, 1] bounds
    for x1 in range(M + 1):
        for x2 in range(N + 1):
            # below: replace
            if sample_counts[x1, x2] < thresh_OB:
                bounds[:, x1, x2] = [0.0, 1.0]
            # above: update truncation
            else:
                # check if smaller than current min
                if x1 < m_OB:
                    m_OB = x1
                    thresh_flag = True
                if x2 < n_OB:
                    n_OB = x2
                    thresh_flag = True
                # check if larger than current max
                if x1 > M_OB:
                    M_OB = x1
                    thresh_flag = True
                if x2 > N_OB:
                    N_OB = x2
                    thresh_flag = True

    # if no states were above threshold: default to max range, report
    if not thresh_flag:
        m_OB, M_OB, n_OB, N_OB = 0, M, 0, N

    # plotting
    if plot:
        fig, axs = plt.subplots(M + 1, N + 1, figsize=(10, 10))
        fig.tight_layout()
        for x1 in range(M + 1):
            for x2 in range(N + 1):
                # within truncation: green CI lines
                if (x1 >= m_OB) and (x2 >= n_OB) and (x1 <= M_OB) and (x2 <= N_OB):
                    color = "green"
                else:
                    color = "red"
                axs[x1, x2].hist(counts[:, x1, x2]);
                axs[x1, x2].set_title(f"p({x1}, {x2})")
                axs[x1, x2].axvline(bounds[0, x1, x2], color=color)
                axs[x1, x2].axvline(bounds[1, x1, x2], color=color)

        plt.suptitle("X1 X2 Confidence Intervals")
        plt.show()

    if printing:
        print(f"Box truncation: [{m_OB}, {M_OB}] x [{n_OB}, {N_OB}]")

    results =  {
        'samples': samples,
        'sample_counts': sample_counts,
        'joint': bounds,
        'm_OB': m_OB,
        'M_OB': M_OB,
        'n_OB': n_OB,
        'N_OB': N_OB,
        'thresh_flag': thresh_flag
    }

    return results

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

def optimization_hyp(bounds, beta, truncations, K=100, silent=True,
                     print_solution=True, print_truncation=True, thresh_OG=10**-6,
                     time_limit=300):

    # load WLS license credentials
    options = json.load(open("../../../WLS_credentials.json"))
    
    # environment context
    with gp.Env(params=options) as env:

        # model context
        with gp.Model('birth-death-regulation-capture-efficiency-hyp', env=env) as md:

            # set options
            if silent:
                md.Params.LogToConsole = 0

            # set time limit: 5 minute default
            md.Params.TimeLimit = time_limit

            # State space truncations

            # observed truncations: computed during bootstrap
            m_OB = bounds['m_OB']
            M_OB = bounds['M_OB']
            n_OB = bounds['n_OB']
            N_OB = bounds['N_OB']

            # original truncations: find largest original states needed (to define variables)
            min_x1_OG, max_x1_OG, min_x2_OG, max_x2_OG = np.inf, 0, np.inf, 0

            # for each pair of observed states used
            for x1_OB in range(m_OB, M_OB + 1):
                for x2_OB in range(n_OB, N_OB + 1):

                    try:
                        # lookup original truncation
                        m_OG, M_OG, n_OG, N_OG = truncations[(x1_OB, x2_OB)]

                    except KeyError:
                        # compute if not available
                        m_OG, M_OG, n_OG, N_OG = findTrunc(x1_OB, x2_OB, beta, thresh_OG)

                        # store
                        truncations[(x1_OB, x2_OB)] = (m_OG, M_OG, n_OG, N_OG)

                    # if larger than current maximum states: update
                    if M_OG > max_x1_OG:
                        max_x1_OG = M_OG
                    if N_OG > max_x2_OG:
                        max_x2_OG = N_OG

                    # if smaller than current minimum states: update
                    if m_OG < min_x1_OG:
                        min_x1_OG = m_OG
                    if n_OG < min_x2_OG:
                        min_x2_OG = n_OG
            
            if print_truncation:
                print(f"Observed counts: [{m_OB}, {M_OB}] x [{n_OB}, {N_OB}] \nOriginal counts: [{min_x1_OG}, {max_x1_OG}] x [{min_x2_OG}, {max_x2_OG}]")

            # variables

            # marginal stationary distributions: original counts (size = largest original state + 1)
            p1 = md.addMVar(shape=(max_x1_OG + 1), vtype=GRB.CONTINUOUS, name="p1", lb=0, ub=1)
            p2 = md.addMVar(shape=(max_x2_OG + 1), vtype=GRB.CONTINUOUS, name="p2", lb=0, ub=1)

            # dummy joint variable to avoid triple products (as not supported by GUROBI): should be removed by presolve
            p_dummy = md.addMVar(shape=(max_x1_OG + 1, max_x2_OG + 1), vtype=GRB.CONTINUOUS, name="p_dummy", lb=0, ub=1)

            '''aggressive presolve to hopefully ensure this'''
            md.Params.Presolve = 2

            # reaction rate constants
            rate_names = ['k_tx_1', 'k_tx_2', 'k_deg_1', 'k_deg_2']
            rates = md.addVars(rate_names, vtype=GRB.CONTINUOUS, lb=0, ub=K, name=rate_names)

            # constraints

            # fix k_deg_1 = 1, k_deg = 2 for identifiability
            md.addConstr(rates['k_deg_1'] == 1)
            md.addConstr(rates['k_deg_2'] == 1)

            # stationary distribution bounds: for each observed count pair
            for x1_OB in range(m_OB, M_OB + 1):
                for x2_OB in range(n_OB, N_OB + 1):
                    
                    # individual truncation: lookup from pre-computed dict
                    m_OG, M_OG, n_OG, N_OG = truncations[(x1_OB, x2_OB)]
                    
                    # sum over truncation range (INCLUSIVE): drop terms with coefficients < thresh
                    sum_expr = gp.quicksum([B(x1_OB, x2_OB, x1_OG, x2_OG, beta) * p1[x1_OG] * p2[x2_OG] for x1_OG in range(m_OG, M_OG + 1) for x2_OG in range(n_OG, N_OG + 1) if B(x1_OB, x2_OB, x1_OG, x2_OG, beta) >= thresh_OG])
                    
                    md.addConstr(sum_expr >= bounds['joint'][0, x1_OB, x2_OB], name=f"B lb {x1_OB}, {x2_OB}")
                    md.addConstr(sum_expr <= bounds['joint'][1, x1_OB, x2_OB], name=f"B ub {x1_OB}, {x2_OB}")

            # distributional constraints
            md.addConstr(p1.sum() <= 1, name="Distribution x1")
            md.addConstr(p2.sum() <= 1, name="Distribution x2")

            # equate dummy joint variable to product of marginals: all original states
            for x1_OG in range(max_x1_OG + 1):
                for x2_OG in range(max_x2_OG + 1):

                    md.addConstr(p_dummy[x1_OG, x2_OG] == p1[x1_OG] * p2[x2_OG], name=f"Dummy joint definition {x1_OG}, {x2_OG}")

            # CME: use dummy joint variable to avoid triple products: k_[] * p1[] * p2[]
            for x1_OG in range(max_x1_OG):
                for x2_OG in range(max_x2_OG):

                    # remove terms when x's = 0 as not present in equation
                    if x1_OG == 0:
                        x1_zero = 0
                    else:
                        x1_zero = 1
                    if x2_OG == 0:
                        x2_zero = 0
                    else:
                        x2_zero = 1

                    md.addConstr(
                        rates['k_tx_1'] * x1_zero * p_dummy[x1_OG - 1, x2_OG] + \
                        rates['k_tx_2'] * x2_zero * p_dummy[x1_OG, x2_OG - 1] + \
                        rates['k_deg_1'] * (x1_OG + 1) * p_dummy[x1_OG + 1, x2_OG] + \
                        rates['k_deg_2'] * (x2_OG + 1) * p_dummy[x1_OG, x2_OG + 1] - \
                        (rates['k_tx_1'] + rates['k_tx_2'] + \
                        rates['k_deg_1'] * x1_OG + rates['k_deg_2'] * x2_OG) * p_dummy[x1_OG, x2_OG] == 0,
                        name=f"CME {x1_OG}, {x2_OG}"
                        )

            # status of optimization
            status_codes = {1: 'LOADED',
                            2: 'OPTIMAL',
                            3: 'INFEASIBLE',
                            4: 'INF_OR_UNBD',
                            5: 'UNBOUNDED',
                            6: 'CUTOFF',
                            7: 'ITERATION_LIMIT',
                            8: 'NODE_LIMIT',
                            9: 'TIME_LIMIT',
                            10: 'SOLUTION_LIMIT',
                            11: 'INTERRUPTED',
                            12: 'NUMERIC',
                            13: 'SUBOPTIMAL',
                            14: 'INPROGRESS',
                            15: 'USER_OBJ_LIMIT'}

            # solution dict
            solution = {
                'status': None,
                'k_tx_1': "v",
                'k_tx_2': "v",
                'k_deg_1': 1,
                'k_deg_2': 1
            }

            # testing feasibility: simply optimize 0
            md.setObjective(0, GRB.MINIMIZE)

            # set parameter (prevents 'infeasible or unbounded' ambiguity)
            md.Params.DualReductions = 0

            # set solution limit (stop after finding 1 feasible solution)
            md.Params.SolutionLimit = 1

            try:
                md.optimize()
                status_code = md.status
            except:
                status_code = md.status

            # store result
            solution['status'] = status_codes[status_code]

            # print
            if print_solution:
                for key, val in solution.items():
                    if key == "status":
                        print(f"Model is {val}")
                    elif val == "v":
                        print(f"{key} variable")
                    else:
                        print(f"{key} = {val}")

            # return model for IIS, etc
            solution['model'] = md

    return solution

def optimization_min(bounds, beta, truncations, K=100, silent=True,
                     print_solution=True, print_truncation=True, thresh_OG=10**-6,
                     MIPGap=0.05, time_limit=300, BestBdThresh=0.0001):

    # load WLS license credentials
    options = json.load(open("../../WLS_credentials.json"))
    
    # environment context
    with gp.Env(params=options) as env:

        # model context
        with gp.Model('birth-death-regulation-capture-efficiency-min', env=env) as md:

            # set options
            if silent:
                md.Params.LogToConsole = 0

            # set time limit: 5 minute default
            md.Params.TimeLimit = time_limit

            # optimization settings
            md.Params.MIPGap = MIPGap

            '''experimental'''
            # aggressive presolve
            md.Params.Presolve = 2
            # focus on lower bound of objective: allows early termination
            md.Params.MIPFocus = 3

            # set threshold on BestBd for termination
            md.Params.BestBdStop = BestBdThresh

            # State space truncations

            # observed truncations: computed during bootstrap
            m_OB = bounds['m_OB']
            M_OB = bounds['M_OB']
            n_OB = bounds['n_OB']
            N_OB = bounds['N_OB']

            # original truncations: find largest original states needed (to define variables)
            min_x1_OG, max_x1_OG, min_x2_OG, max_x2_OG = np.inf, 0, np.inf, 0

            # for each pair of observed states used
            for x1_OB in range(m_OB, M_OB + 1):
                for x2_OB in range(n_OB, N_OB + 1):

                    try:
                        # lookup original truncation
                        m_OG, M_OG, n_OG, N_OG = truncations[(x1_OB, x2_OB)]

                    except KeyError:
                        # compute if not available
                        m_OG, M_OG, n_OG, N_OG = findTrunc(x1_OB, x2_OB, beta, thresh_OG)

                        # store
                        truncations[(x1_OB, x2_OB)] = (m_OG, M_OG, n_OG, N_OG)

                    # if larger than current maximum states: update
                    if M_OG > max_x1_OG:
                        max_x1_OG = M_OG
                    if N_OG > max_x2_OG:
                        max_x2_OG = N_OG

                    # if smaller than current minimum states: update
                    if m_OG < min_x1_OG:
                        min_x1_OG = m_OG
                    if n_OG < min_x2_OG:
                        min_x2_OG = n_OG
            
            if print_truncation:
                print(f"Observed counts: [{m_OB}, {M_OB}] x [{n_OB}, {N_OB}] \nOriginal counts: [{min_x1_OG}, {max_x1_OG}] x [{min_x2_OG}, {max_x2_OG}]")

            # variables

            # stationary distribution: original counts (size = largest truncation)
            p = md.addMVar(shape=(max_x1_OG + 1, max_x2_OG + 1), vtype=GRB.CONTINUOUS, name="p", lb=0, ub=1)

            # reaction rate constants
            rate_names = ['k_tx_1', 'k_tx_2', 'k_deg_1', 'k_deg_2', 'k_reg']
            rates = md.addVars(rate_names, vtype=GRB.CONTINUOUS, lb=0, ub=K, name=rate_names)

            # constraints

            # fix k_deg_2 = 1 for identifiability
            md.addConstr(rates['k_deg_2'] == 1)

            # stationary distribution bounds: for each observed count
            for x1_OB in range(m_OB, M_OB + 1):
                for x2_OB in range(n_OB, N_OB + 1):

                    # individual truncation: lookup from pre-computed dict
                    m_OG, M_OG, n_OG, N_OG = truncations[(x1_OB, x2_OB)]
                    
                    # sum over truncation range (INCLUSIVE): drop terms with coefficients < thresh
                    sum_expr = gp.quicksum([B(x1_OB, x2_OB, x1_OG, x2_OG, beta) * p[x1_OG, x2_OG] for x1_OG in range(m_OG, M_OG + 1) for x2_OG in range(n_OG, N_OG + 1) if B(x1_OB, x2_OB, x1_OG, x2_OG, beta) >= thresh_OG])
                    
                    md.addConstr(sum_expr >= bounds['joint'][0, x1_OB, x2_OB], name=f"B lb {x1_OB}, {x2_OB}")
                    md.addConstr(sum_expr <= bounds['joint'][1, x1_OB, x2_OB], name=f"B ub {x1_OB}, {x2_OB}")

            # distribution
            md.addConstr(p.sum() <= 1, name="Distribution")
            
            # stationary Qp=0 equations for all relevant variables
            for x1_OG in range(max_x1_OG):
                for x2_OG in range(max_x2_OG):

                    # remove terms when x's = 0 as not present in equation
                    if x1_OG == 0:
                        x1_zero = 0
                    else:
                        x1_zero = 1
                    if x2_OG == 0:
                        x2_zero = 0
                    else:
                        x2_zero = 1

                    md.addConstr(
                        rates['k_tx_1'] * x1_zero * p[x1_OG - 1, x2_OG] + \
                        rates['k_tx_2'] * x2_zero * p[x1_OG, x2_OG - 1] + \
                        rates['k_deg_1'] * (x1_OG + 1) * p[x1_OG + 1, x2_OG] + \
                        rates['k_deg_2'] * (x2_OG + 1) * p[x1_OG, x2_OG + 1] + \
                        rates['k_reg'] * (x1_OG + 1) * (x2_OG + 1) * p[x1_OG + 1, x2_OG + 1] - \
                        (rates['k_tx_1'] + rates['k_tx_2'] + \
                        rates['k_deg_1'] * x1_OG + rates['k_deg_2'] * x2_OG + \
                        rates['k_reg'] * x1_OG * x2_OG) * p[x1_OG, x2_OG] == 0,
                        name=f"Equation {x1_OG}, {x2_OG}"
                        )

            # status of optimization
            status_codes = {1: 'LOADED',
                            2: 'OPTIMAL',
                            3: 'INFEASIBLE',
                            4: 'INF_OR_UNBD',
                            5: 'UNBOUNDED',
                            6: 'CUTOFF',
                            7: 'ITERATION_LIMIT',
                            8: 'NODE_LIMIT',
                            9: 'TIME_LIMIT',
                            10: 'SOLUTION_LIMIT',
                            11: 'INTERRUPTED',
                            12: 'NUMERIC',
                            13: 'SUBOPTIMAL',
                            14: 'INPROGRESS',
                            15: 'USER_OBJ_LIMIT'}

            # solution dict
            solution = {
                'k_tx_1': "v",
                'k_tx_2': "v",
                'k_deg_1': "v",
                'k_deg_2': 1,
                'k_reg': None
            }

            # Optimize

            # set objective: minimize interaction parameter
            md.setObjective(rates['k_reg'], GRB.MINIMIZE)

            # attempt to optimize
            try:
                md.optimize()
                min_val = md.ObjVal
            except:
                min_val = None

            # report status
            status_min = status_codes[md.status]

            # record
            solution['k_reg'] = [min_val, status_min]

            # print
            if print_solution:
                for key, val in solution.items():
                    if key == 'k_reg':
                        if val[1] == 'USER_OBJ_LIMIT':
                            print(f"{key} non-zero lower bound found, early termination status {val[1]}")
                        else:
                            print(f"{key} lower bound {val[0]}, status {val[1]}")
                    elif val == "v":
                        print(f"{key} variable, not optimized")
                    else:
                        print(f"{key} = {val}")

            # return model for IIS, etc
            solution['model'] = md

    return solution

# ------------------------------------------------
# Input
# ------------------------------------------------

counts_df = pd.read_csv(input_filename, index_col=0, converters={f'Cell-{j}': literal_eval for j in range(1000)})

# ------------------------------------------------
# Computation
# ------------------------------------------------

solution_dict = {}

# loop over dataset
for i in range(100):

    # display progress
    print(i)

    # select sample
    samples = list(counts_df.loc[f'Gene-pair-{i}'])

    if method == "hyp":

        # bootstrap
        bounds = bootstrap(samples, BS=1000, thresh=thresh_OB, printing=False)

        # optimize: hyp
        solution = optimization_B_hyp(bounds, beta, truncations, K=100, silent=True,
                            print_solution=True, print_truncation=True, thresh_trunc=10**-6,
                            time_limit=300)

        # store result
        solution_dict[i] = {'status': solution['status'], 'time': solution['model'].Runtime}

    elif method == "min":

        # bootstrap
        bounds = bootstrap(samples, BS=1000, thresh=thresh_OB, printing=False)

        # optimize: min
        solution = optimization_min(bounds, beta, truncations, K=100, silent=True,
                        print_solution=True, print_truncation=True, thresh_trunc=10**-6,
                        MIPGap=0.05, time_limit=300, BestBdThresh=0.0001)

        # store result
        solution_dict[i] = {'bound': solution['k_reg'][0], 'status': solution['k_reg'][1], 'time': solution['model'].Runtime}

    elif method = "pearson":

        # select individual samples
        x1_samples = [x[0] for x in samples]
        x2_samples = [x[1] for x in samples]

        # test
        pearson = scipy.stats.pearsonr(x1_samples, x2_samples)

        # store result
        solution_dict[i] = {'pvalue': float(pearson.statistic), 'statistic': float(pearson.pvalue)}

    elif method == "spearman":

        # select individual samples
        x1_samples = [x[0] for x in samples]
        x2_samples = [x[1] for x in samples]

        # test
        spearman = scipy.stats.spearmanr(x1_samples, x2_samples)

        # store result
        solution_dict[i] = {'pvalue': float(spearman.statistic), 'statistic': float(spearman.pvalue)}


# ------------------------------------------------
# Output
# ------------------------------------------------

json.dump(solution_dict, open(output_filename, 'w'))
