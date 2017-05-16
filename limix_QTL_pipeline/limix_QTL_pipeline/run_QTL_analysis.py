import pandas as pd
import numpy as np
import limix
import qtl_output
import glob
import os
from sklearn.preprocessing import Imputer

def run_QTL_analysis(pheno_filename,anno_filename,geno_prefix,chromsome,window_size,output_dir,
                     covariates_filename=None,kinship_filename=None,sample_mapping_filename=None):
    '''Core function to take input and run QTL tests on a given chromosome.'''
    
    
    #Load input data files    
    phenotype_df = pd.read_csv(pheno_filename,sep='\t',index_col=0)
    annotation_col_dtypes = {'feature_id':np.object,
                             'gene_id':np.object,
                             'gene_name':np.object,
                             'chromosome':np.object,
                             'start':np.int64,
                             'end':np.int64,
                             'strand':np.object}
    annotation_df = pd.read_csv(anno_filename,sep='\t',index_col=0,dtype=annotation_col_dtypes)
        
    bim,fam,bed = limix.io.read_plink(geno_prefix,verbose=False)
    fam.set_index('iid',inplace=True)
    
    if kinship_filename:
        kinship_df = pd.read_csv(kinship_filename,sep='\t',index_col=0)
    else:
        kinship_df = None
    
    if covariates_filename:
        covariate_df = pd.read_csv(covariates_filename,sep='\t',index_col=0)
    else:
        covariate_df = None
    
    if sample_mapping_filename:
        individual2sample_df = pd.read_csv(sample_mapping_filename,sep='\t',header=None,names=['iid','sample'],index_col=0)
    else:
        #assume the mapping is the identity mapping
        identifiers = list(phenotype_df.columns)
        individual2sample_df = pd.DataFrame(data=identifiers,index=identifiers,columns=['sample'])

    #Open output files
    _ensure_dir(output_dir)
    output_writer = qtl_output.text_writer(output_dir+'qtl_results_{}.txt'.format(chromosome))

    
    #Determine features to be tested
    feature_list = list(set(annotation_df[annotation_df['chromosome']==chromosome].index)&set(phenotype_df.index))
    
    #Test features
    for feature_id in feature_list:
        
        chrom = str(annotation_df.loc[feature_id,'chromosome'])
        start = annotation_df.loc[feature_id,'start']
        end = annotation_df.loc[feature_id,'end']
        center_pos = start + (start-end)/2
        cis = bim.query("chrom == '%s' & pos > %d & pos < %d" % (chrom, center_pos-window_size, center_pos+window_size))
        snp_idxs = cis['i'].values
        snp_names = cis['snp'].values
    
        #indices for relevant individuals in genotype matrix
        individual_ids = list(set(fam.index)&set(individual2sample_df.index))
        individual_idxs = fam.loc[individual_ids,'i'].values
    
        #subset genotype matrix and kinship matrix
        snps = bed[snp_idxs,:].compute()
        if len(snps) == 0:
            continue
        snp_matrix = snps[:,individual_idxs].transpose()
        snp_matrix = Imputer(missing_values=3.,strategy='mean',axis=0,copy=False).fit_transform(snp_matrix)
    
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
    
        temp_df = pd.DataFrame(index = range(len(snp_names)),columns=['feature_id','snp_id','p_value','beta'])
        temp_df['snp_id'] = snp_names
        temp_df['feature_id'] = feature_id
        temp_df['beta'] = LMM.getBetaSNP()[0]
        temp_df['p_value'] = LMM.getPv()[0]
        output_writer.add_result_df(temp_df)
    
    output_writer.close()
    
    #write annotation and snp data to file
    snp_df = pd.DataFrame()
    snp_df['snp_id'] = bim['snp']
    snp_df['chromosome'] = bim['chrom']
    snp_df['position'] = bim['pos']
    
    snp_df.to_csv(output_dir+'/snp_metadata.txt',sep='\t',index=False)
    annotation_df.to_csv(output_dir+'/feature_metadata.txt',sep='\t')

def merge_QTL_results(results_dir):
    '''Merge QTL results for individual chromosomes into a combined, indexed
    hdf5 file.'''
    qtl_results_files = glob.glob(results_dir+'qtl_results_*.txt')
    
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
    
if __name__=='__main__':
    '''Run a test case'''
    data_path = '../data/geuvadis_CEU_YRI_test_data/'
    covariates_filename = data_path+'Geuvadis_CEU_YRI_covariates.txt'
    geno_prefix = data_path+'Geuvadis_chr1'
    pheno_filename = data_path+'Geuvadis_CEU_YRI_Expr.txt'
    anno_filename = data_path+'Geuvadis_CEU_YRI_formatted_annotation_data.txt'
    kinship_filename= data_path+'Geuvadis_chr1_kinship.txt'
    individual2sample_filename = data_path + 'Geuvadis_CEU_gte.txt'
    
    output_dir = data_path+'limix_QTL_results_kinship_covs/'
    
    chromosome = '1'
    
    ws = 250000
    
    run_QTL_analysis(pheno_filename,anno_filename,geno_prefix,chromosome,ws,output_dir,
                     covariates_filename=covariates_filename,
                     kinship_filename=kinship_filename,
                     sample_mapping_filename=individual2sample_filename)

    data_path = '../data/geuvadis_CEU_test_data/'
    geno_prefix = data_path+'Genotypes/Geuvadis'
    pheno_filename = data_path+'Expression/Geuvadis_CEU_Expr.txt'
    anno_filename = data_path+'Expression/Geuvadis_CEU_formatted_annotation_data.txt'
    
    output_dir = data_path+'TestOutput/limix_QTL_results/'
        
    ws = 250000
    
    for chromosome in ['1','2']:
        run_QTL_analysis(pheno_filename,anno_filename,geno_prefix,chromosome,ws,output_dir)
    merge_QTL_results(output_dir)
