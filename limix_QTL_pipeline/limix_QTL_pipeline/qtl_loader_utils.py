import limix
import pandas as pd
import numpy as np
import os

##loader functions
def ensure_dir(file_path):
    '''Check if directory exists for output, and create it if it doesn't.'''
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)

def get_kinship_df(kinship_filename):
    if kinship_filename:
        kinship_df = pd.read_csv(kinship_filename,sep='\t',index_col=0)
    else:
        kinship_df = None
    return kinship_df

def get_samplemapping_df(sample_mapping_filename,sample_labels,key_from):
    assert(key_from in ['iid','sample'])
    if sample_mapping_filename:
        mapping_df = pd.read_csv(sample_mapping_filename,sep='\t',header=None,names=['iid','sample'])
        mapping_df.set_index(key_from,inplace=True)
    else:
        #assume the mapping is the identity mapping
        identifiers = sample_labels
        mapping_df = pd.DataFrame(data=identifiers,index=identifiers,columns=['iid','sample'])
    return mapping_df

def get_covariate_df(covariates_filename):
    if covariates_filename:
        covariate_df = pd.read_csv(covariates_filename,sep='\t',index_col=0)
    else:
        covariate_df = None
    return covariate_df

def get_genotype_data(geno_prefix):
    bim,fam,bed = limix.io.read_plink(geno_prefix,verbose=False)
    fam.set_index('iid',inplace=True)
    return bim,fam,bed

def get_annotation_df(anno_filename):
    annotation_col_dtypes = {'feature_id':np.object,
                         'gene_id':np.object,
                         'gene_name':np.object,
                         'chromosome':np.object,
                         'start':np.int64,
                         'end':np.int64,
                         'strand':np.object}
    annotation_df = pd.read_csv(anno_filename,sep='\t',index_col=0,dtype=annotation_col_dtypes)
    return annotation_df

def get_phenotype_df(pheno_filename):
    return pd.read_csv(pheno_filename,sep='\t',index_col=0)
