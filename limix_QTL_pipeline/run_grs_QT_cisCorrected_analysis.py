from __future__ import division
import pdb
import pandas as pd
import numpy as np
import qtl_output
import qtl_loader_utils
import qtl_parse_args
import qtl_utilities as utils
import glob
from sklearn.preprocessing import Imputer
import sys
from numpy_sugar.linalg import economic_qs, economic_svd
from limix.stats import effsizes_se, lrt_pvalues
from glimix_core.lmm import LMM

#V0.1.2

def run_PrsQtl_analysis(pheno_filename, anno_filename, prsFile, geno_prefix, plinkGenotype, output_dir, blocksize=1000,skipAutosomeFiltering = False, gaussianize_method=None, minimum_test_samples= 10, 
                        seed=np.random.randint(40000), n_perm=0, write_permutations = False, relatedness_score=0.95, feature_variant_covariate_filename = None, snps_filename=None,
                        feature_filename=None, snp_feature_filename=None, genetic_range='all', covariates_filename=None, kinship_filename=None, sample_mapping_filename=None):
    fill_NaN = Imputer(missing_values=np.nan, strategy='mean', axis=0)
    print('Running GRS QT analysis.')
    lik = 'normal'
    '''Core function to take input and run QTL tests on a given chromosome.'''
    [phenotype_df, kinship_df, covariate_df, sample2individual_df, annotation_df, snp_filter_df, snp_feature_filter_df, geneticaly_unique_individuals, minimum_test_samples, feature_list,bim,fam,bed, risk_df, chromosome, selectionStart, selectionEnd, feature_variant_covariate_df]=\
    utils.run_PrsQtl_CisCorrected_analysis_load_intersect_phenotype_covariates_kinship_sample_mapping(pheno_filename=pheno_filename, anno_filename=anno_filename, geno_prefix=geno_prefix, plinkGenotype=plinkGenotype, prsFile=prsFile, skipAutosomeFiltering = skipAutosomeFiltering,
                      minimum_test_samples= minimum_test_samples,  relatedness_score=relatedness_score, snps_filename=snps_filename, feature_filename=feature_filename, snp_feature_filename=snp_feature_filename, selection=genetic_range,
                     covariates_filename=covariates_filename, kinship_filename=kinship_filename, sample_mapping_filename=sample_mapping_filename, feature_variant_covariate_filename=feature_variant_covariate_filename)
    
    mixed = kinship_df is not None
    QS = None
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
                if(any(i in  bim['snp'].values for i in covariateSnp)):
                    snpQuery_cov = bim.loc[bim['snp'].map(lambda x: x in list(covariateSnp)),:]
                    snp_cov_df_t = pd.DataFrame(data=bed[snpQuery_cov['i'].values,:].compute().transpose(),index=fam.index,columns=snpQuery_cov['snp'],)
                    snp_cov_df = pd.DataFrame(fill_NaN.fit_transform(snp_cov_df_t))
                    snp_cov_df.index=snp_cov_df_t.index
                    snp_cov_df.columns=snp_cov_df_t.columns
                    snp_cov_df_t = None
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
            
            '''select indices for relevant individuals in genotype matrix
            These are not unique. NOT to be used to access phenotype/covariates data
            '''
            individual_ids = sample2individual_df.loc[phenotype_ds.index,'iid'].values
            sample2individual_feature= sample2individual_df.loc[phenotype_ds.index]
            
            if contains_missing_samples:
                tmp_unique_individuals = geneticaly_unique_individuals
                geneticaly_unique_individuals = utils.get_unique_genetic_samples(kinship_df.loc[individual_ids,individual_ids], relatedness_score);
            
            if phenotype_ds.empty or len(geneticaly_unique_individuals)<minimum_test_samples :
                print("Feature: "+feature_id+" not tested not enough samples do QTL test.")
                fail_qc_features.append(feature_id)
                if contains_missing_samples:
                    geneticaly_unique_individuals = tmp_unique_individuals
                continue
            elif np.var(phenotype_ds.values) == 0:
                print("Feature: "+feature_id+" has no variance in selected individuals.")
                fail_qc_features.append(feature_id)
                if contains_missing_samples:
                    geneticaly_unique_individuals = tmp_unique_individuals
                continue
            
            print ('For feature: ' +str(currentFeatureNumber)+ '/'+str(len(feature_list))+ ' (' + feature_id + '): ' + str(snpQuery.shape[0]) + ' risk scores will be tested.\n Please stand by.')
            if(n_perm!=0):
                bestPermutationPval = np.ones((n_perm), dtype=np.float)
            
            #Here we need to start preparing the LMM, can use the fam for sample IDS in SNP matrix.
#                test if the covariates, kinship, snp and phenotype are in the same order
            if ((all(kinship_df.loc[individual_ids,individual_ids].index==sample2individual_feature.loc[phenotype_ds.index]['iid']) if kinship_df is not None else True) &\
                 (all(phenotype_ds.index==covariate_df.loc[sample2individual_feature['sample'],:].index)if covariate_df is not None else True)):
                '''
                if all lines are in order put in arrays the correct genotype and phenotype
                x=a if cond1 else b <---> equivalent to if cond1: x=a else x=b;                 better readability of the code
                 '''
                if kinship_df is not None:
                    kinship_mat = kinship_df.loc[individual_ids,individual_ids].values
                    kinship_mat = kinship_mat.astype(float)
                    ##GOWER normalization of Kinship matrix.
                    kinship_mat *= (kinship_mat.shape[0] - 1) / (kinship_mat.trace() - kinship_mat.mean(0).sum())
                    ## This needs to go with the subselection stuff.
                    if(QS is None and not contains_missing_samples):
                        QS = economic_qs(kinship_mat)
                    elif (contains_missing_samples):
                        QS_tmp = QS
                        QS = economic_qs(kinship_mat)
                if kinship_df is None:
                    K = eye(len(phenotype_ds.index))
                    if(QS is None and not contains_missing_samples):
                        QS = economic_qs(kinship_mat)
                    elif (contains_missing_samples):
                        QS_tmp = QS
                        QS = economic_qs(kinship_mat)
                cov_matrix =  covariate_df.loc[sample2individual_feature['sample'],:].values if covariate_df is not None else None
                if covariate_df is None:
                    cov_matrix = ones((len(individual_ids), 1))
                if snp_cov_df is not None:
                    snp_cov_df_tmp = snp_cov_df.loc[individual_ids,:]
                    snp_cov_df_tmp.index=sample2individual_feature['sample']
                    cov_matrix = np.concatenate((cov_matrix,snp_cov_df_tmp.values),1)
                    print("added snp information as covariate.")
                cov_matrix = cov_matrix.astype(float)
            else:
                print ('There is an issue in mapping phenotypes vs covariates and/or kinship')
                sys.exit()
            
            phenotype = utils.force_normal_distribution(phenotype_ds.values,method=gaussianize_method) if gaussianize_method is not None else phenotype_ds.values
            
            #Prepare LMM
            phenotype = phenotype.astype(float)
            
            
            ##Mixed and test.
            ##This is a future change so we don't need to decompose the COVs every time. 
            ##Like QS this needs to happen when genetic unique individuals is the same.
            #svd_cov = economic_svd(cov_matrix)
            #lmm = LMM(phenotype, cov_matrix, QS, SVD=svd_cov)
            #These steps need to happen only once per phenotype.
            #print(QS)
            lmm = LMM(phenotype, cov_matrix, QS)
            if not mixed:
                lmm.delta = 1
                lmm.fix('delta')
            #Prepare null model.
            lmm.fit(verbose=False)
            null_lml = lmm.lml()
            flmm = lmm.get_fast_scanner()
            
            for snpGroup in utils.chunker(snpQuery, blocksize):
                #Fix seed at the start of the first chunker so all permutations are based on the same random first split.
                np.random.seed(seed)
                
                snp_names = snpGroup.values
                
                tested_snp_names.extend(snp_names)
                snp_matrix_DF = risk_df.loc[snp_names,individual_ids].transpose()
                ##GRS var QC
                snp_matrix_DF = snp_matrix_DF.loc[:,(np.std(snp_matrix_DF,axis=0)>0).values]
#               test if the covariates, kinship, snp and phenotype are in the same order
                if (len(snp_matrix_DF.index) != len(sample2individual_feature.loc[phenotype_ds.index]['iid']) or not all(snp_matrix_DF.index==sample2individual_feature.loc[phenotype_ds.index]['iid'])):
                    print ('There is an issue in mapping phenotypes and genotypes')
                    sys.exit()
                
                G = snp_matrix_DF.values
                G = G.astype(float)
                G_index = snp_matrix_DF.columns
                
                alt_lmls, effsizes = flmm.fast_scan(G, verbose=False)
                var_pvalues = lrt_pvalues(null_lml, alt_lmls)
                var_effsizes_se = effsizes_se(effsizes, var_pvalues)
                
                #add these results to qtl_results
                temp_df = pd.DataFrame(index = range(len(G_index)),columns=['feature_id','snp_id','p_value','beta','beta_se','empirical_feature_p_value'])
                temp_df['snp_id'] = G_index
                temp_df['feature_id'] = feature_id
                temp_df['beta'] = np.asarray(effsizes)
                temp_df['p_value'] = np.asarray(var_pvalues)
                temp_df['beta_se'] = np.asarray(var_effsizes_se)
                #insert default dummy value
                temp_df['empirical_feature_p_value'] = -1.0
                
                if(n_perm!=0):
                    pValueBuffer = []
                    totalSnpsToBeTested = (G.shape[1]*n_perm)
                    permutationStepSize = np.floor(totalSnpsToBeTested/blocksize)
                    if(permutationStepSize==0):
                        permutationStepSize=1
                    
                    if(write_permutations):
                        perm_df = pd.DataFrame(index = range(len(G_index)),columns=['snp_id'] + ['permutation_'+str(x) for x in range(n_perm)])
                        perm_df['snp_id'] = G_index
                    for currentNperm in utils.chunker(list(range(1, n_perm+1)), permutationStepSize):
                        if kinship_df is not None:
                            temp = utils.get_shuffeld_genotypes_preserving_kinship(geneticaly_unique_individuals, relatedness_score, snp_matrix_DF,kinship_df.loc[individual_ids,individual_ids], len(currentNperm))
                        else :
                            temp = utils.get_shuffeld_genotypes(snp_matrix_DF, len(currentNperm))
                        temp = temp.astype(float)
                        alt_lmls_p, effsizes_p = flmm.fast_scan(temp, verbose=False)
                        var_pvalues_p = lrt_pvalues(null_lml, alt_lmls_p)
                        pValueBuffer.extend(np.asarray(var_pvalues_p))
                    if(not(len(pValueBuffer)==totalSnpsToBeTested)):
                        #print(len(pValueBuffer))
                        #print(pValueBuffer)
                        #print(totalSnpsToBeTested)
                        print('Error in blocking logic for permutations.')
                        sys.exit()
                    perm = 0
                    for relevantOutput in utils.chunker(pValueBuffer,G.shape[1]) :
                        if(write_permutations):
                            perm_df['permutation_'+str(perm)] = relevantOutput
                        if(bestPermutationPval[perm] > min(relevantOutput)):
                            bestPermutationPval[perm] = min(relevantOutput)
                        perm = perm+1
                        #print(relevantOutput)
                        #print('permutation_'+str(perm))
                
                if not temp_df.empty :
                    data_written = True
                    output_writer.add_result_df(temp_df)
                    if(write_permutations):
                        permutation_writer.add_permutation_results_df(perm_df,feature_id)
            #This we need to change in the written file.
            if(n_perm>1 and data_written):
                #updated_permuted_p_in_hdf5(bestPermutationPval, feature_id);
                alpha_para, beta_para = output_writer.apply_pval_correction(feature_id,bestPermutationPval)
                alpha_params.append(alpha_para)
                beta_params.append(beta_para)
            if not data_written :
                fail_qc_features.append(feature_id)
            else:
                n_samples.append(phenotype_ds.size)
                n_e_samples.append(len(geneticaly_unique_individuals))
            if contains_missing_samples:
                QS = QS_tmp
                geneticaly_unique_individuals = tmp_unique_individuals
                del QS_tmp
                del tmp_unique_individuals
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
    args = qtl_parse_args.get_grsQtl_cis_corrected_args()
    plink  = args.plink
    bgen = args.bgen
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
    
    if ((plink is None) and (bgen is None)):
        raise ValueError("No genotypes provided. Either specify a path to a binary plink genotype file or a bgen file.")
    if ((plink is not None) and (bgen is not None)):
        raise ValueError("Only one genotype file can be provided at once, not both plink and bgen")

    if (bgen is not None) :
        plinkGenotype=False
        geno_prefix = bgen
        raise ValueError("Not supported")
    else:
        plinkGenotype=True
        geno_prefix = plink
    
    if ((grsFile is None)):
        raise ValueError("No risk scores provided. Either specify a path to Genetic risk score flat-file.")
    
    if (random_seed is None):
        random_seed = np.random.randint(40000)
    
    if(n_perm==0 and write_permutations):
        write_permutations=False
    
    if(n_perm==1):
        print("Warning: With only 1 permutation P-value correction is not performed.")
    if(n_perm<50):
        print("Warning: With less than 50 permutations P-values correction is not very accurate.")
    run_PrsQtl_analysis(pheno_file, anno_file, grsFile,geno_prefix, plinkGenotype, output_dir, blocksize=int(block_size), skipAutosomeFiltering= includeAllChromsomes, gaussianize_method = gaussianize,
                     minimum_test_samples= int(minimum_test_samples), seed=int(random_seed), n_perm=int(n_perm), write_permutations = write_permutations, relatedness_score=float(relatedness_score), 
                     feature_variant_covariate_filename = feature_variant_covariate_filename, snps_filename=snps_filename, feature_filename=feature_filename, snp_feature_filename=snp_feature_filename, 
                     genetic_range=genetic_range, covariates_filename=covariates_file, kinship_filename=kinship_file, sample_mapping_filename=samplemap_file)