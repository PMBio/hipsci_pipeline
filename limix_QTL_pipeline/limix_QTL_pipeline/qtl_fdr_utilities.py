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
    if(len(top_pvalues_perm)<5){
        #If only a small number of features don't use the MLE estimator
        alpha_para,beta_para = estimate_beta_function_paras(top_pvalues_perm)
    } else {
        #Use the MLE estimator
        alpha_para,beta_para,loc,fscale =  beta.fit(top_pvalues_perm,floc=0,fscale=1)
    }
    beta_dist = scipy.stats.beta(alpha_para,beta_para)
    correction_function = lambda x: beta_dist.cdf(x)
    return correction_function

def calculate_corrected_pvalues(top_pvalues_perm,nominal_pvalues):
    if(np.mean(top_pvalues_perm)==1):
        return nominal_pvalues
    correction_function = define_correction_function(top_pvalues_perm)
    #apply correction to nominal p_values - potentially slow
    corrected_pvalues = np.array([correction_function(x) for x in nominal_pvalues])
    #corrected_pvalues = Parallel(n_jobs=-2)(np.array([correction_function(x) for x in nominal_pvalues]))
    return corrected_pvalues
