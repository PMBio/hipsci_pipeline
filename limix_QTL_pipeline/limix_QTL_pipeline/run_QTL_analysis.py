import pandas as pd
import numpy as np
import limix


#def run_QTL_analysis(pheno_filename,anno_filename,geno_prefix,covariates_filename=None,kinship_filename=None,sample_mapping_filename=None):

data_path = '../data/geuvadis_CEU_YRI_test_data/'
covariates_filename = data_path+'Geuvadis_CEU_YRI_covariates.txt'
geno_prefix = data_path+'Geuvadis_chr1'
pheno_filename = data_path+'Geuvadis_CEU_YRI_Expr.txt'
anno_filename = data_path+'Geuvadis_CEU_YRI_formatted_annotation_data.txt'
kinship_filename= data_path+'Geuvadis_chr1_kinship.txt'

output_dir = data_path+'limix_QTL_results'


phenotype_df = pd.read_csv(pheno_filename,sep='\t',index_col=0)
annotation_df = pd.read_csv(anno_filename,sep='\t',index_col=0)





bim,fam,bed = limix.io.read_plink(geno_prefix,verbose=False)
fam.set_index('iid',inplace=True)

if kinship_filename:
    kinship_df = pd.read_csv(kinship_filename,sep='\t',index_col=0)

if covariates_filename:
    covariate_df = pd.read_csv(covariates_filename,sep='\t',index_col=0)
 
individual2sample_df = pd.read_csv('../data/geuvadis_CEU_YRI_test_data/Geuvadis_CEU_gte.txt',sep='\t',header=None,names=['iid','sample'],index_col=0)

qtl_results_df = pd.DataFrame(columns=['feature_id','snp_id','p_value','beta'])


ws = 250000
feature_list = list(set(annotation_df.index)&set(phenotype_df.index))

for feature_id in feature_list:
    
    chrom = str(annotation_df.loc[feature_id,'chromosome'])
    start = annotation_df.loc[feature_id,'start']
    end = annotation_df.loc[feature_id,'end']
    center_pos = start + (start-end)/2
    cis = bim.query("chrom == '%s' & pos > %d & pos < %d" % (chrom, center_pos-ws, center_pos+ws))
    snp_idxs = cis['i'].values
    snp_names = cis['snp'].values

    #make sure to only look up the non NA valued samples for this feature.
    
    #indices for relevant individuals in genotype matrix
    individual_ids = list(set(fam.index)&set(individual2sample_df.index))
    individual_idxs = fam.loc[individual_ids,'i'].values

    #subset genotype matrix and kinship matrix
    snps = bed[snp_idxs,:].compute().transpose()
    if len(snps) == 0:
        continue
    snps = snps[individual_idxs,:]
    
    kinship_mat = kinship_df.loc[individual_ids,individual_ids].as_matrix()

    #map individual_ids to samples
    sample_ids = individual2sample_df.loc[individual_ids,'sample'].values
    phenotype = phenotype_df.loc[feature_id,sample_ids].as_matrix()
    
    
    #generate covariate matrix
    cov_matrix = np.concatenate([np.ones((len(sample_ids),1)),covariate_df.loc[sample_ids,:].as_matrix()],axis=1)
    
    #fit modelrun
    LMM = limix.qtl.qtl_test_lmm(snps, phenotype,K=kinship_mat)
    
    #add these results to qtl_results

    temp_df = pd.DataFrame(index = range(len(snp_names)),columns=['feature_id','snp_id','p_value','beta'])
    temp_df['snp_id'] = snp_names
    temp_df['feature_id'] = feature_id
    temp_df['beta'] = LMM.getBetaSNP()[0]
    temp_df['p_value'] = LMM.getPv()[0]
    qtl_results_df = qtl_results_df.append(temp_df, ignore_index=True)

snp_df = pd.DataFrame()
snp_df['snp_id'] = bim['snp']
snp_df['chromosome'] = bim['chrom']
snp_df['position'] = bim['pos']

snp_df.to_csv(output_dir+'/snp_metadata.txt',sep='\t',index=False)
annotation_df.to_csv(output_dir+'/feature_metadata.txt',sep='\t')
qtl_results_df.to_csv(output_dir+'/qtl_results.txt',sep='\t',index=False)