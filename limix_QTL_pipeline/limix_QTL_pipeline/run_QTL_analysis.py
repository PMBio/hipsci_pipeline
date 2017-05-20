import pandas as pd
import numpy as np
import limix
import qtl_output
import glob
import os
from sklearn.preprocessing import Imputer
import argparse


def get_args():
    parser = argparse.ArgumentParser(description='Run QTL analysis given genotype, phenotype, and annotation.')
    parser.add_argument('-geno_prefix','--geno_prefix', required=True)
    parser.add_argument('-anno_file','--anno_file', required=True)
    parser.add_argument('-pheno_file','--pheno_file', required=True)
    parser.add_argument('-output_dir','--output_dir', required=True)
    parser.add_argument('-cis_window_kb','--cis_window_kb', required=True,
                        help=
                        'The size of the cis window to take SNPs from, in kb.'
                        'The window will extend between:                     '
                        '    (feature_start - (cis_window_kb/2))             '
                        ' and:                                               '
                        '    (feature_end + (cis_window_kb/2))               ')
    parser.add_argument('-chromosome','--chromosome',required=True)
    parser.add_argument('-covariates_file','--covariates_file',required=False)
    parser.add_argument('-kinship_file','--kinship_file',required=False)
    parser.add_argument('-samplemap_file','--samplemap_file',required=False)
    parser.add_argument("--cis",
                        action="store_true",
                        help="Run cis analysis.")
    parser.add_argument("--trans",
                        action="store_true",
                        help="Run trans analysis.")

    args = parser.parse_args()

    return args


def run_QTL_analysis(pheno_filename,anno_filename,geno_prefix,chromosome,window_size,output_dir,
                     covariates_filename=None,kinship_filename=None,sample_mapping_filename=None):
    '''Core function to take input and run QTL tests on a given chromosome.'''
    
    
    #Load input data files    
    phenotype_df = _get_phenotype_df(pheno_filename)
    annotation_df = _get_annotation_df(anno_filename)
    bim,fam,bed = _get_genotype_data(geno_prefix)
    kinship_df = _get_kinship_df(kinship_filename)    
    covariate_df = _get_covariate_df(covariates_filename)
    individual2sample_df = _get_samplemapping_df(sample_mapping_filename,list(phenotype_df.columns))

    #Open output files
    _ensure_dir(output_dir)
    output_writer = qtl_output.hdf5_writer(output_dir+'qtl_results_{}.h5'.format(chromosome))

    #Determine features to be tested
    feature_list = list(set(annotation_df[annotation_df['chromosome']==chromosome].index)&set(phenotype_df.index))
    
    #Array to store indices of snps tested
    tested_snp_idxs = []
    
    #Test features
    for feature_id in feature_list:
        
        chrom = str(annotation_df.loc[feature_id,'chromosome'])
        start = annotation_df.loc[feature_id,'start']
        end = annotation_df.loc[feature_id,'end']
        #make robust to features specified back-to-front
        lowest = min([start,end])
        highest = max([start,end])
        half_window_size = window_size/2
        cis = bim.query("chrom == '%s' & pos > %d & pos < %d" % (chrom, lowest-half_window_size, highest+half_window_size))
        snp_idxs = cis['i'].values
        snp_names = cis['snp'].values
        
        tested_snp_idxs.extend(snp_idxs)
        
        #indices for relevant individuals in genotype matrix
        individual_ids = list(set(fam.index)&set(individual2sample_df.index))
        individual_idxs = fam.loc[individual_ids,'i'].values
    
        #subset genotype matrix and kinship matrix
        snps = bed[snp_idxs,:].compute()
        if len(snps) == 0:
            continue
        snp_matrix = snps[:,individual_idxs].transpose()
        snp_matrix = Imputer(missing_values=np.nan,strategy='mean',axis=0,copy=False).fit_transform(snp_matrix)
    
        if kinship_df is not None:
            kinship_mat = kinship_df.loc[individual_ids,individual_ids].as_matrix()
        else:
            kinship_mat = None
    
        #map individual_ids to samples
        sample_ids = individual2sample_df.loc[individual_ids,'sample'].values
        phenotype = phenotype_df.loc[feature_id,sample_ids].as_matrix()
        
        #generate covariate matrix
        if covariate_df is not None:
            cov_matrix = np.concatenate([np.ones((len(sample_ids),1)),covariate_df.loc[sample_ids,:].as_matrix()],axis=1)
        else:
            cov_matrix = None
        
        #fit modelrun
        LMM = limix.qtl.qtl_test_lmm(snp_matrix, phenotype,K=kinship_mat,covs=cov_matrix)
        
        #add these results to qtl_results
    
        temp_df = pd.DataFrame(index = range(len(snp_names)),columns=['feature_id','snp_id','p_value','beta','n_samples'])
        temp_df['snp_id'] = snp_names
        temp_df['feature_id'] = feature_id
        temp_df['beta'] = LMM.getBetaSNP()[0]
        temp_df['p_value'] = LMM.getPv()[0]
        temp_df['n_samples'] = sum(~np.isnan(phenotype))
        output_writer.add_result_df(temp_df)
    
    output_writer.close()
    
    #gather unique indexes of tested snps
    tested_snp_idxs = list(set(tested_snp_idxs))
    #write annotation and snp data to file
    snp_df = pd.DataFrame()
    snp_df['snp_id'] = bim['snp']
    snp_df['chromosome'] = bim['chrom']
    snp_df['position'] = bim['pos']
    snp_df['assessed_allele'] = bim['a1']
    
    snp_df.ix[tested_snp_idxs,:].to_csv(output_dir+'/snp_metadata_{}.txt'.format(chromosome),sep='\t',index=False)
    annotation_df.loc[feature_list,:].to_csv(output_dir+'/feature_metadata_{}.txt'.format(chromosome),sep='\t')

def merge_QTL_results(results_dir):
    '''Merge QTL results for individual chromosomes into a combined, indexed
    hdf5 file.'''
    qtl_results_files = sorted(glob.glob(results_dir+'qtl_results_*.txt'))
    
    hdf5_outfile = qtl_output.hdf5_writer(results_dir+'qtl_results.h5')
    
    for filename in qtl_results_files:
        df = pd.read_csv(filename,sep='\t')
        hdf5_outfile.add_result_df(df)
    
    hdf5_outfile.close()

def _ensure_dir(file_path):
    '''Check if directory exists for output, and create it if it doesn't.'''
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)

def _get_kinship_df(kinship_filename):
    if kinship_filename:
        kinship_df = pd.read_csv(kinship_filename,sep='\t',index_col=0)
    else:
        kinship_df = None
    return kinship_df

def _get_samplemapping_df(sample_mapping_filename,sample_labels):
    if sample_mapping_filename:
        individual2sample_df = pd.read_csv(sample_mapping_filename,sep='\t',header=None,names=['iid','sample'],index_col=0)
    else:
        #assume the mapping is the identity mapping
        identifiers = sample_labels
        individual2sample_df = pd.DataFrame(data=identifiers,index=identifiers,columns=['sample'])
    return individual2sample_df

def _get_covariate_df(covariates_filename):
    if covariates_filename:
        covariate_df = pd.read_csv(covariates_filename,sep='\t',index_col=0)
    else:
        covariate_df = None
    return covariate_df

def _get_genotype_data(geno_prefix):
    bim,fam,bed = limix.io.read_plink(geno_prefix,verbose=False)
    fam.set_index('iid',inplace=True)
    return bim,fam,bed

def _get_annotation_df(anno_filename):
    annotation_col_dtypes = {'feature_id':np.object,
                         'gene_id':np.object,
                         'gene_name':np.object,
                         'chromosome':np.object,
                         'start':np.int64,
                         'end':np.int64,
                         'strand':np.object}
    annotation_df = pd.read_csv(anno_filename,sep='\t',index_col=0,dtype=annotation_col_dtypes)
    return annotation_df

def _get_phenotype_df(pheno_filename):
    return pd.read_csv(pheno_filename,sep='\t',index_col=0)

if __name__=='__main__':
    args = get_args()
    geno_prefix = args.geno_prefix
    anno_file = args.anno_file
    pheno_file = args.pheno_file
    output_dir = args.output_dir
    cis_window_kb = args.cis_window_kb
    chromosome = args.chromosome
    covariates_file = args.covariates_file
    kinship_file = args.kinship_file
    samplemap_file = args.samplemap_file
    cis = args.cis
    trans = args.trans
    
    if (cis and trans):
        raise ValueError("cis and trans cannot be specified simultaneously")
    if not (cis or trans):
        raise ValueError("One of cis and trans must be specified")

    window_size = int(cis_window_kb)*1000

    run_QTL_analysis(pheno_file,anno_file,geno_prefix,chromosome,window_size,output_dir,
                     covariates_filename=covariates_file,
                     kinship_filename=kinship_file,
                     sample_mapping_filename=samplemap_file)