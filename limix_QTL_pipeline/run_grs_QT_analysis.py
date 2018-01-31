import pandas as pd
import numpy as np
import limix
import qtl_output
import qtl_loader_utils
import qtl_parse_args
import qtl_utilities as utils
import glob
from sklearn.preprocessing import Imputer
import scipy.stats as scst
import sys

#V0.1.1

def run_PrsQtl_analysis(pheno_filename, anno_filename, prsFile, output_dir, blocksize=1000,
                     skipAutosomeFiltering = False, gaussianize_method=None, minimum_test_samples= 10, seed=np.random.randint(40000), n_perm=0, write_permutations = False, relatedness_score=0.95, feature_variant_covariate_filename = None, snps_filename=None, feature_filename=None, snp_feature_filename=None, genetic_range='all',
                     covariates_filename=None, kinship_filename=None, sample_mapping_filename=None):
    fill_NaN = Imputer(missing_values=np.nan, strategy='mean', axis=0)
    print('Running GRS QT analysis.')
    '''Core function to take input and run QTL tests on a given chromosome.'''

    [phenotype_df, kinship_df, covariate_df, sample2individual_df, annotation_df, snp_filter_df, snp_feature_filter_df, geneticaly_unique_individuals, minimum_test_samples, feature_list, risk_df, chromosome, selectionStart, selectionEnd, feature_variant_covariate_df]=\
    utils.run_PrsQtl_analysis_load_intersect_phenotype_covariates_kinship_sample_mapping(pheno_filename=pheno_filename, anno_filename=anno_filename, prsFile=prsFile, skipAutosomeFiltering = skipAutosomeFiltering,
                      minimum_test_samples= minimum_test_samples,  relatedness_score=relatedness_score, snps_filename=snps_filename, feature_filename=feature_filename, snp_feature_filename=snp_feature_filename, selection=genetic_range,
                     covariates_filename=covariates_filename, kinship_filename=kinship_filename, sample_mapping_filename=sample_mapping_filename, feature_variant_covariate_filename=feature_variant_covariate_filename)
    
    if(feature_list==None or len(feature_list)==0):
        print ('No features to be tested.')
        sys.exit()
    
    #Open output files
    qtl_loader_utils.ensure_dir(output_dir)
    if not selectionStart is None :
        output_writer = qtl_output.hdf5_writer(output_dir+'qtl_results_{}_{}_{}.h5'.format(chromosome,selectionStart,selectionEnd))
    else :
        output_writer = qtl_output.hdf5_writer(output_dir+'qtl_results_{}.h5'.format(chromosome))
    if(write_permutations):
        if not selectionStart is None :
            permutation_writer = qtl_output.hdf5_permutations_writer(output_dir+'perm_results_{}_{}_{}.h5'.format(chromosome,selectionStart,selectionEnd),n_perm)
        else :
            permutation_writer = qtl_output.hdf5_permutations_writer(output_dir+'perm_results_{}.h5'.format(chromosome),n_perm)

    #Arrays to store indices of snps tested and pass and fail QC SNPs for features without missingness.
    tested_snp_names = []
    fail_qc_features = []
    alpha_params = []
    beta_params = []
    n_samples = []
    n_e_samples = []
    currentFeatureNumber = 0
    for feature_id in feature_list:
        currentFeatureNumber+= 1
        if (len(phenotype_df.loc[feature_id,:]))<minimum_test_samples:
            print("Feature: "+feature_id+" not tested not enough samples do QTL test.")
            continue
        data_written = False

        snpQuery = risk_df.index
        snp_cov_df = None
        if(feature_variant_covariate_df is not None):
            if(feature_id in feature_variant_covariate_df['feature'].values):
                covariateSnp = feature_variant_covariate_df['snp_id'].values[feature_variant_covariate_df['feature']==feature_id]
                if(any(i in  risk_df.index.values for i in covariateSnp)):
                    snpQuery_cov = risk_df.index.loc[risk_df.index.map(lambda x: x in list(covariateSnp)),:]
                    snp_cov_df = risk_df.index.loc[risk_df.index.map(lambda x: x in list(covariateSnp)),:].transpose()

        if (len(snpQuery) != 0) and (snp_filter_df is not None):
            snpQuery = snpQuery.loc[snpQuery['snp'].map(lambda x: x in list(map(str, snp_filter_df.index)))]

        if (len(snpQuery) != 0) and (snp_feature_filter_df is not None):
            snpQuery = snpQuery.loc[snpQuery['snp'].map(lambda x: x in list(snp_feature_filter_df['snp_id'].loc[snp_feature_filter_df['feature']==feature_id]))]

        if len(snpQuery) != 0:

            phenotype_ds = phenotype_df.loc[feature_id]
            contains_missing_samples = any(~np.isfinite(phenotype_ds))
            if(contains_missing_samples):
                print ('Feature: ' + feature_id + ' contains missing data.')
                phenotype_ds.dropna(inplace=True)
#   
            '''select indices for relevant individuals in genotype matrix
            These are not unique. NOT to be used to access phenotype/covariates data
            '''
            individual_ids = sample2individual_df.loc[phenotype_ds.index,'iid'].values
            sample2individual_feature= sample2individual_df.loc[phenotype_ds.index]

            if phenotype_ds.empty or len(geneticaly_unique_individuals)<minimum_test_samples :
                print("Feature: "+feature_id+" not tested not enough samples do QTL test.")
                fail_qc_features.append(feature_id)
                continue
            elif np.var(phenotype_ds.values) == 0:
                print("Feature: "+feature_id+" has no variance in selected individuals.")
                fail_qc_features.append(feature_id)
                continue
            
            if contains_missing_samples:
                tmp_unique_individuals = geneticaly_unique_individuals
                geneticaly_unique_individuals = utils.get_unique_genetic_samples(kinship_df.loc[individual_ids,individual_ids], relatedness_score);
            
            n_samples.append(phenotype_ds.size)
            n_e_samples.append(len(geneticaly_unique_individuals))
            
            print ('For feature: ' +str(currentFeatureNumber)+ '/'+str(len(feature_list))+ ' (' + feature_id + '): ' + str(snpQuery.shape[0]) + ' risk scores will be tested.\n Please stand by.')
            if(n_perm!=0):
                bestPermutationPval = np.ones((n_perm), dtype=np.float)
            for snpGroup in utils.chunker(snpQuery, blocksize):
                snp_names = snpGroup.values
                
                tested_snp_names.extend(snp_names)
                snp_matrix_DF = risk_df.loc[snp_names,individual_ids].transpose()
                ##GRS var QC
                snp_matrix_DF = snp_matrix_DF.loc[(np.std(snp_df,axis=0)>0),:]
#               test if the covariates, kinship, snp and phenotype are in the same order
                if ((all(snp_matrix_DF.index==kinship_df.loc[individual_ids,individual_ids].index) if kinship_df is not None else True) &\
                     (all(phenotype_ds.index==covariate_df.loc[sample2individual_feature['sample'],:].index)if covariate_df is not None else True)&\
                     all(snp_matrix_DF.index==sample2individual_feature.loc[phenotype_ds.index]['iid'])):
                    '''
                    if all lines are in order put in arrays the correct genotype and phenotype
                    x=a if cond1 else b <---> equivalent to if cond1: x=a else x=b;                 better readability of the code
                     '''
                    kinship_mat = kinship_df.loc[individual_ids,individual_ids].values if kinship_df is not None else None
                    cov_matrix =  covariate_df.loc[sample2individual_feature['sample'],:].values if covariate_df is not None else None
#                    cov_matrix =  covariate_df[covariate_df.columns.values[np.array([('peer' in c)|(c==feature_id) for c in  covariate_df.columns.values])]].loc[sample2individual_feature['sample'],:].values if covariate_df is not None else None
                    
                    if(snp_cov_df is not None and cov_matrix is not None):
                        snp_cov_df_tmp = snp_cov_df.loc[individual_ids,:]
                        snp_cov_df_tmp.index=sample2individual_feature['sample']
                        cov_matrix = np.concatenate((cov_matrix,snp_cov_df_tmp.values),1)
                    elif snp_cov_df is not None :
                        snp_cov_df_tmp = snp_cov_df.loc[individual_ids,:]
                        snp_cov_df_tmp.index=sample2individual_feature['sample']
                        cov_matrix = snp_cov_df_tmp.values
                        #cov_matrix = np.concatenate((np.ones(snp_cov_df_tmp.shape[0]).reshape(np.ones(snp_cov_df_tmp.shape[0]).shape[0],1),snp_cov_df_tmp.values),1)

                    phenotype = utils.force_normal_distribution(phenotype_ds.values,method=gaussianize_method) if gaussianize_method is not None else phenotype_ds.values
                else:
                    print ('There is an issue in mapping phenotypes and genotypes')
                    sys.exit()
                
                #For limix 1.1 we need to switch to lm our selfs if there is no K.
                #return[snp_matrix_DF,phenotype, kinship_mat,cov_matrix]
                #sys.exit()
                try: 
                    LMM = limix.qtl.scan(snp_matrix_DF.values, phenotype, 'Normal', K=kinship_mat,M=cov_matrix,verbose=False)
                except: 
                    print (feature_id)
                    print ('LMM failed')
                
                if(n_perm!=0):
                    pValueBuffer = []
                    totalSnpsToBeTested = (snp_matrix_DF.shape[1]*n_perm)
                    perm = 0;
                    if(write_permutations):
                        perm_df = pd.DataFrame(index = range(len(snp_matrix_DF.columns)),columns=['snp_id'] + ['permutation_'+str(x) for x in range(n_perm)])
                        perm_df['snp_id'] = snp_matrix_DF.columns
                    for currentNperm in utils.chunker(list(range(1, n_perm+1)), np.floor(totalSnpsToBeTested/blocksize)):
                        if kinship_df is not None:
                            temp = utils.get_shuffeld_genotypes_preserving_kinship(geneticaly_unique_individuals, relatedness_score, snp_matrix_DF,kinship_df.loc[individual_ids,individual_ids], len(currentNperm))
                        else :
                            temp = utils.get_shuffeld_genotypes(snp_matrix_DF, len(currentNperm))
                        LMM_perm = limix.qtl.scan(temp, phenotype, 'Normal',K=kinship_mat,M=cov_matrix,verbose=False)
                        pValueBuffer.extend(np.asarray(LMM_perm.variant_pvalues))
                    if(not(len(pValueBuffer)==totalSnpsToBeTested)):
                        #print(len(pValueBuffer))
                        #print(pValueBuffer)
                        #print(totalSnpsToBeTested)
                        print('Error in blocking logic for permutations.')
                        sys.exit()
                    for relevantOutput in utils.chunker(pValueBuffer,snp_matrix_DF.shape[1]) :
                        if(write_permutations):
                            perm_df['permutation_'+str(perm)] = relevantOutput
                        if(bestPermutationPval[perm] > min(relevantOutput)):
                            bestPermutationPval[perm] = min(relevantOutput)
                        perm+=1

                #add these results to qtl_results
                temp_df = pd.DataFrame(index = range(len(snp_matrix_DF.columns)),columns=['feature_id','snp_id','p_value','beta','beta_se','empirical_feature_p_value'])
                temp_df['snp_id'] = snp_matrix_DF.columns
                temp_df['feature_id'] = feature_id
                temp_df['beta'] = np.asarray(LMM.variant_effsizes)
                temp_df['p_value'] = np.asarray(LMM.variant_pvalues)
                temp_df['beta_se'] = np.asarray(LMM.variant_effsizes_se)
                #insert default dummy value
                temp_df['empirical_feature_p_value'] = -1.0
                
                if not temp_df.empty :
                    data_written = True
                    output_writer.add_result_df(temp_df)
                    if(write_permutations):
                        permutation_writer.add_permutation_results_df(perm_df,feature_id)
            if contains_missing_samples:
                geneticaly_unique_individuals = tmp_unique_individuals
        
        if(n_perm>1 and data_written):
            #updated_permuted_p_in_hdf5(bestPermutationPval, feature_id);
            alpha_para, beta_para = output_writer.apply_pval_correction(feature_id,bestPermutationPval)
            alpha_params.append(alpha_para)
            beta_params.append(beta_para)
        if not data_written :
            fail_qc_features.append(feature_id)
        #print('step 5')
    output_writer.close()
    if(write_permutations):
        permutation_writer.close()

    #gather unique indexes of tested snps
    #write annotation and snp data to file
    snp_df = pd.DataFrame()
    snp_df['snp_id'] = list(set(tested_snp_names))

    feature_list = [x for x in feature_list if x not in fail_qc_features]
    annotation_df = annotation_df.loc[feature_list,:]
    annotation_df['n_samples'] = n_samples
    annotation_df['n_e_samples'] = n_e_samples
    if(n_perm>1 and data_written):
        annotation_df['alpha_param'] = alpha_params
        annotation_df['beta_param'] = beta_params
    if not selectionStart is None :
        snp_df.to_csv(output_dir+'/snp_metadata_{}_{}_{}.txt'.format(chromosome,selectionStart,selectionEnd),sep='\t',index=False)
        annotation_df.to_csv(output_dir+'/feature_metadata_{}_{}_{}.txt'.format(chromosome,selectionStart,selectionEnd),sep='\t')
    else :
        snp_df.to_csv(output_dir+'/snp_metadata_{}.txt'.format(chromosome),sep='\t',index=False)
        annotation_df.to_csv(output_dir+'/feature_metadata_{}.txt'.format(chromosome),sep='\t')

if __name__=='__main__':
    args = qtl_parse_args.get_grsQtl_args()
    grsFile  = args.genetic_risk_scores
    anno_file = args.annotation_file
    pheno_file = args.phenotype_file
    output_dir = args.output_directory
    genetic_range = args.genomic_range
    covariates_file = args.covariates_file
    kinship_file = args.kinship_file
    samplemap_file = args.sample_mapping_file
    block_size = args.block_size
    n_perm = int(args.number_of_permutations)
    snps_filename = args.variant_filter
    snp_feature_filename = args.feature_variant_filter
    feature_variant_covariate_filename = args.feature_variant_covariate
    random_seed = args.seed
    feature_filename = args.feature_filter
    relatedness_score = args.relatedness_score
    minimum_test_samples = args.minimum_test_samples
    gaussianize = args.gaussianize_method
    write_permutations = args.write_permutations
    includeAllChromsomes = args.no_chromosome_filter

    if ((grsFile is None)):
        raise ValueError("No risk scores provided. Either specify a path to Genetic risc score flatfile.")

    if (random_seed is None):
        random_seed = np.random.randint(40000)

    if(n_perm==0 and write_permutations):
        write_permutations=False

    if(n_perm>1 and n_perm<10):
        n_perm=10
        print("Defaults to 10 permutations, if permutations are only used for calibration please give in 1.")

    run_PrsQtl_analysis(pheno_file, anno_file, grsFile, output_dir, blocksize=int(block_size), skipAutosomeFiltering= includeAllChromsomes, gaussianize_method = gaussianize,
                     minimum_test_samples= int(minimum_test_samples), seed=int(random_seed), n_perm=int(n_perm), write_permutations = write_permutations, relatedness_score=float(relatedness_score), 
                     feature_variant_covariate_filename = feature_variant_covariate_filename, snps_filename=snps_filename, feature_filename=feature_filename, snp_feature_filename=snp_feature_filename, 
                     genetic_range=genetic_range, covariates_filename=covariates_file, kinship_filename=kinship_file, sample_mapping_filename=samplemap_file)