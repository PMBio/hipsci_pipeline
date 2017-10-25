import numpy as np
import scipy.stats
from scipy.stats import beta
#from joblib import Parallel

#V0.1.1

def estimate_beta_function_paras(top_pvalues_perm):
    mean = np.mean(top_pvalues_perm)
    variance = np.var(top_pvalues_perm)
    alpha_para = mean * (mean * (1 - mean ) / variance - 1)
    beta_para = alpha_para * (1 / mean - 1)
    return alpha_para,beta_para

def define_correction_function(top_pvalues_perm):
    #Always try to use the MLE estimator, new default to 10 permutations.
    #If the MLE estimator fails we go back to the cruder estimation of the beta distribution.
    offset = (np.finfo(np.double).tiny*100)
    ##Replace lowest value with smallest number not 0.
    top_pvalues_perm[top_pvalues_perm <= 0] = offset
    ##Replace highest value with highest number not 1.
    top_pvalues_perm[top_pvalues_perm >= 1] = 1-offset
    try :
        alpha_para,beta_para,loc,fscale =  beta.fit(top_pvalues_perm,floc=0,fscale=1)
    except (scipy.stats._continuous_distns.FitSolverError):
        alpha_para,beta_para = estimate_beta_function_paras(top_pvalues_perm)
    except (scipy.stats._continuous_distns.FitDataError):
        alpha_para,beta_para = estimate_beta_function_paras(top_pvalues_perm)
    beta_dist = scipy.stats.beta(alpha_para,beta_para)
    correction_function = lambda x: beta_dist.cdf(x)
    return correction_function

def calculate_corrected_pvalues(top_pvalues_perm,nominal_pvalues):
    correction_function = define_correction_function(top_pvalues_perm)
    #apply correction to nominal p_values - potentially slow
    corrected_pvalues = np.array([correction_function(x) for x in nominal_pvalues])
    #corrected_pvalues = Parallel(n_jobs=-2)(np.array([correction_function(x) for x in nominal_pvalues]))
    return corrected_pvalues
