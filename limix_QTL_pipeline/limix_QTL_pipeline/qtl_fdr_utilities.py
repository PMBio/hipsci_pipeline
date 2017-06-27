import numpy as np
import scipy.stats

def estimate_beta_function_paras(top_pvalues):
    mean = np.mean(top_pvalues)
    variance = np.var(top_pvalues)
    alpha_para = mean * (mean * (1 - mean ) / variance - 1)
    beta_para = alpha_para * (1 / mean - 1)
    return alpha_para,beta_para

def calculate_corrected_pvalues(nominal_pvalues,top_pvalues):
    alpha_para,beta_para = estimate_beta_function_paras(top_pvalues)
    beta_dist = scipy.stats.beta(alpha_para,beta_para)
    corrected_pvalues = [beta_dist.cdf(x) for x in nominal_pvalues]
    return corrected_pvalues