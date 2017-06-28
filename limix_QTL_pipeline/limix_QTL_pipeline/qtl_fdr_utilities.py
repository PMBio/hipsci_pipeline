import numpy as np
import scipy.stats
from joblib import Parallel

def estimate_beta_function_paras(top_pvalues_perm):
    mean = np.mean(top_pvalues_perm)
    variance = np.var(top_pvalues_perm)
    alpha_para = mean * (mean * (1 - mean ) / variance - 1)
    beta_para = alpha_para * (1 / mean - 1)
    return alpha_para,beta_para

def calculate_corrected_pvalues(top_pvalues_perm,nominal_pvalues):
    alpha_para,beta_para = estimate_beta_function_paras(top_pvalues_perm)
    beta_dist = scipy.stats.beta(alpha_para,beta_para)
    #apply correction to nominal p_values - potentially slow
    corrected_pvalues = np.array([beta_dist.cdf(x) for x in nominal_pvalues])
    #corrected_pvalues = Parallel(n_jobs=-2)(np.array([beta_dist.cdf(x) for x in nominal_pvalues]))
    return corrected_pvalues
