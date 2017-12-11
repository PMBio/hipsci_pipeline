import pandas as pd
import numpy as np
import limix
import qtl_output
import qtl_loader_utils
import qtl_utilities as utils
from qtl_snp_qc import do_snp_qc
import glob
from sklearn.preprocessing import Imputer
import argparse
import scipy.stats as scst
import sys

#V0.1.1

def run_QTL_analysis(pheno_filename, anno_filename, geno_prefix, plinkGenotype, output_dir, window_size=250000, min_maf=0.05, min_hwe_P=0.001, min_call_rate=0.95, blocksize=1000,
                     cis_mode=True, gaussianize_method=None, minimum_test_samples= 10, seed=np.random.randint(40000), n_perm=0, write_permutations = False, relatedness_score=0.95, feature_variant_covariate_filename = None, snps_filename=None, feature_filename=None, snp_feature_filename=None, genetic_range='all',
                     covariates_filename=None, kinship_filename=None, sample_mapping_filename=None, extended_anno_filename=None):
    fill_NaN = Imputer(missing_values=np.nan, strategy='mean', axis=0)
    '''Core function to take input and run QTL tests on a given chromosome.'''

    [phenotype_df, kinship_df, covariate_df, sample2individual_df,complete_annotation_df, annotation_df, snp_filter_df, snp_feature_filter_df, geneticaly_unique_individuals, minimum_test_samples, feature_list,bim,fam,bed, chromosome, selectionStart, selectionEnd, feature_variant_covariate_df]=\
    utils.run_QTL_analysis_load_intersect_phenotype_covariates_kinship_sample_mapping(pheno_filename=pheno_filename, anno_filename=anno_filename, geno_prefix=geno_prefix, plinkGenotype=plinkGenotype, cis_mode=cis_mode,
                      minimum_test_samples= minimum_test_samples,  relatedness_score=relatedness_score, snps_filename=snps_filename, feature_filename=feature_filename, snp_feature_filename=snp_feature_filename, selection=genetic_range,
                     covariates_filename=covariates_filename, kinship_filename=kinship_filename, sample_mapping_filename=sample_mapping_filename, extended_anno_filename=extended_anno_filename, feature_variant_covariate_filename=feature_variant_covariate_filename)
    #
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
    tested_snp_idxs = []
    pass_qc_snps_all = []
    fail_qc_snps_all = []
    fail_qc_features = []
    for feature_id in feature_list:
        if (len(phenotype_df.loc[feature_id,:]))<minimum_test_samples:
            print("Feature: "+feature_id+" not tested not enough samples do QTL test.")
            continue
        data_written = False
        
        snpQuery = utils.do_snp_selection(feature_id, complete_annotation_df, bim, cis_mode, window_size)
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
            if phenotype_ds.empty | len(phenotype_ds)<minimum_test_samples :
                print("Feature: "+feature_id+" not tested not enough samples do QTL test.")
                fail_qc_features.append(feature_id)
                continue
            elif np.var(phenotype_ds.values) == 0:
                print("Feature: "+feature_id+" has no variance in selected individuals.")
                fail_qc_features.append(feature_id)
                continue

            #If no missing samples we can use the previous SNP Qc information before actually loading data.
            #This allowes for more efficient blocking and retreaving of data
            if not contains_missing_samples:
                snpQuery = snpQuery.loc[snpQuery['snp'].map(lambda x: x not in list(map(str, fail_qc_snps_all)))]
            
            print ('For, feature: ' + feature_id + ' ' + str(snpQuery.shape[0]) + ' SNPs need to be tested.\n Please stand by.')
            if(n_perm!=0):
                bestPermutationPval = np.ones((n_perm), dtype=np.float)
            for snpGroup in utils.chunker(snpQuery, blocksize):
                snp_idxs = snpGroup['i'].values
                snp_names = snpGroup['snp'].values

                tested_snp_idxs.extend(snp_idxs)

                #subset genotype matrix, we cannot subselect at the same time, do in two steps.
                snp_df = pd.DataFrame(data=bed[snp_idxs,:].compute().transpose(),index=fam.index,columns=snp_names)
                snp_df = snp_df.loc[individual_ids,:]

                #SNP QC.
                    #Now we do more proper QC on non-identical samples.
                    #However, we do not use it when checking for missingness.
                    #That could be extended but gives alot of overhead.
                if not contains_missing_samples:
                    #remove snps from snp_df if they fail QC
                    snp_df = snp_df.loc[:,snp_df.columns[~snp_df.columns.isin(fail_qc_snps_all)]]
                    if snp_df.shape[1] == 0:
                        continue
                    snps_to_test_df = snp_df.loc[:,snp_df.columns[~snp_df.columns.isin(pass_qc_snps_all)]]
                    if snps_to_test_df.shape[1] > 0:
                        #Only do QC on relevant SNPs. join pre-QCed list and new QCed list.
                        if kinship_df is not None and len(geneticaly_unique_individuals)<snps_to_test_df.shape[0]:
                            passed_snp_names,failed_snp_names = do_snp_qc(snps_to_test_df.iloc[np.unique(snps_to_test_df.index,return_index=1)[1]].loc[geneticaly_unique_individuals,:], min_call_rate, min_maf, min_hwe_P)
                        else:
                            passed_snp_names,failed_snp_names = do_snp_qc(snps_to_test_df, min_call_rate, min_maf, min_hwe_P)
                        snps_to_test_df = None
                        #append snp_names and failed_snp_names
                        pass_qc_snps_all.extend(passed_snp_names)
                        fail_qc_snps_all.extend(failed_snp_names)
                    snp_df = snp_df.loc[:,snp_df.columns[snp_df.columns.isin(pass_qc_snps_all)]]
                else:
                    tmp_unique_individuals = geneticaly_unique_individuals
                    geneticaly_unique_individuals = utils.get_unique_genetic_samples(kinship_df.loc[individual_ids,individual_ids], relatedness_score);
                    #Do snp QC for relevant section.
                    #Get relevant slice from: phenotype_ds
                    if kinship_df is not None and len(geneticaly_unique_individuals)>1 and len(geneticaly_unique_individuals)<snp_df.shape[0]:
                        passed_snp_names,failed_snp_names = do_snp_qc(snp_df.iloc[np.unique(snp_df.index,return_index=1)[1]].loc[geneticaly_unique_individuals,:], min_call_rate, min_maf, min_hwe_P) 
                    else:
                        passed_snp_names,failed_snp_names = do_snp_qc(snp_df, min_call_rate, min_maf, min_hwe_P)
                    snp_df = snp_df.loc[:,snp_df.columns[snp_df.columns.isin(passed_snp_names)]]
                #print('step 0')
                if len(snp_df.columns) == 0:
                    continue
                #We could make use of relatedness when imputing.
                snp_matrix_DF = pd.DataFrame(fill_NaN.fit_transform(snp_df))
                snp_matrix_DF.columns = snp_df.columns
                snp_matrix_DF.index = snp_df.index
                snp_df = None


#                test if the covariates, kinship, snp and phenotype are in the same order
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

                    ### select discrete covariates with at least 6 lines & add vector of one's if necessary:
                    if cov_matrix is not None :
                        #if np.unique(cov_matrix).shape[0]==2:
                        #    cov_matrix=cov_matrix[:,np.nansum(cov_matrix==1,0)>6]
                        if np.isfinite(cov_matrix.sum(0)/cov_matrix.std(0)).all():
                            cov_matrix = np.concatenate((np.ones(cov_matrix.shape[0]).reshape(np.ones(cov_matrix.shape[0]).shape[0],1),cov_matrix.values),1)
                    phenotype = utils.force_normal_distribution(phenotype_ds.values,method=gaussianize_method) if gaussianize_method is not None else phenotype_ds.values
                else:
                    print ('there is an issue in mapping phenotypes and genotypes')
                    sys.exit()
                LMM = limix.qtl.qtl_test_lmm(snp_matrix_DF.values, phenotype,K=kinship_mat,covs=cov_matrix)

                if(n_perm!=0):
                    if(write_permutations):
                        perm_df = pd.DataFrame(index = range(len(snp_matrix_DF.columns)),columns=['snp_id'] + ['permutation_'+str(x+1) for x in range(n_perm)])
                        perm_df['snp_id'] = snp_matrix_DF.columns
                    if kinship_df is not None and len(geneticaly_unique_individuals)<snp_matrix_DF.shape[0]:
                        temp = utils.get_shuffeld_genotypes_preserving_kinship(geneticaly_unique_individuals, relatedness_score, snp_matrix_DF,kinship_df.loc[individual_ids,individual_ids], n_perm)
                        LMM_perm = limix.qtl.qtl_test_lmm(temp, phenotype,K=kinship_mat,covs=cov_matrix)
                        perm = 0;
                        for relevantOutput in utils.chunker(LMM_perm.getPv()[0],snp_matrix_DF.shape[1]) :
                            if(write_permutations):
                                perm_df['permutation_'+str(perm)] = relevantOutput
                            if(bestPermutationPval[perm] > min(relevantOutput)):
                                bestPermutationPval[perm] = min(relevantOutput)
                            perm+=1
                    else :
                        temp = utils.get_shuffeld_genotypes(snp_matrix_DF,kinship_df, n_perm)
                        LMM_perm = limix.qtl.qtl_test_lmm(temp, phenotype,K=kinship_mat,covs=cov_matrix)
                        perm = 0;
                        for relevantOutput in utils.chunker(LMM_perm.getPv()[0],snp_matrix_DF.shape[1]) :
                            if(write_permutations):
                                perm_df['permutation_'+str(perm)] = relevantOutput
                            if(bestPermutationPval[perm] > min(relevantOutput)):
                                bestPermutationPval[perm] = min(relevantOutput)
                            perm+=1
                #print('step 3')
                #add these results to qtl_results
                temp_df = pd.DataFrame(index = range(len(snp_matrix_DF.columns)),columns=['feature_id','snp_id','p_value','beta','n_samples','empirical_feature_p_value'])
                temp_df['snp_id'] = snp_matrix_DF.columns
                temp_df['feature_id'] = feature_id
                temp_df['beta'] = LMM.getBetaSNP()[0]
                temp_df['p_value'] = LMM.getPv()[0]
                temp_df['n_samples'] = sum(~np.isnan(phenotype))
                temp_df['beta_se'] = LMM.getBetaSNPste()[0]
                #insert default dummy value
                temp_df['empirical_feature_p_value'] = -1.0
                if not temp_df.empty :
                    data_written = True
                    output_writer.add_result_df(temp_df)
                    if(write_permutations):
                        permutation_writer.add_permutation_results_df(perm_df,feature_id)
                if contains_missing_samples:
                    geneticaly_unique_individuals = tmp_unique_individuals
                #print('step 4')
            #This we need to change in the written file.
        if(n_perm>1 and data_written):
            #updated_permuted_p_in_hdf5(bestPermutationPval, feature_id);
            output_writer.apply_pval_correction(feature_id,bestPermutationPval)
        if not data_written :
            fail_qc_features.append(feature_id)
        #print('step 5')
    output_writer.close()
    if(write_permutations):
        permutation_writer.close()

    #gather unique indexes of tested snps
    tested_snp_idxs = list(set(tested_snp_idxs))
    #write annotation and snp data to file
    snp_df = pd.DataFrame()
    snp_df['snp_id'] = bim['snp']
    snp_df['chromosome'] = bim['chrom']
    snp_df['position'] = bim['pos']
    snp_df['assessed_allele'] = bim['a1']

    feature_list = [x for x in feature_list if x not in fail_qc_features]
    if not selectionStart is None :
        snp_df.ix[tested_snp_idxs,:].to_csv(output_dir+'/snp_metadata_{}_{}_{}.txt'.format(chromosome,selectionStart,selectionEnd),sep='\t',index=False)
        annotation_df.loc[feature_list,:].to_csv(output_dir+'/feature_metadata_{}_{}_{}.txt'.format(chromosome,selectionStart,selectionEnd),sep='\t')
    else :
        snp_df.ix[tested_snp_idxs,:].to_csv(output_dir+'/snp_metadata_{}.txt'.format(chromosome),sep='\t',index=False)
        annotation_df.loc[feature_list,:].to_csv(output_dir+'/feature_metadata_{}.txt'.format(chromosome),sep='\t')

if __name__=='__main__':
    args = qtl_parse_args.get_args()
    plink  = args.plink
    bgen = args.bgen
    anno_file = args.annotation_file
    extended_anno_file = args.extended_annotation_file
    pheno_file = args.phenotype_file
    output_dir = args.output_directory
    window_size = args.window
    genetic_range = args.genomic_range
    covariates_file = args.covariates_file
    kinship_file = args.kinship_file
    samplemap_file = args.sample_mapping_file
    min_maf = args.minor_allel_frequency
    min_hwe_P = args.hardy_weinberg_equilibrium
    min_call_rate = args.call_rate
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
    cis = args.cis
    trans = args.trans
    write_permutations = args.write_permutations

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

    if (cis and trans):
        raise ValueError("cis and trans cannot be specified simultaneously")
    if (not cis and not trans):
        cis = True
    if (random_seed is None):
        random_seed = np.random.randint(40000)

    if(n_perm==0 and write_permutations):
        write_permutations=False

    if(n_perm>1 and n_perm<10):
        n_perm=10
        print("Defaults to 10 permutations, if permutations are only used for calibration please give in 1.")

    run_QTL_analysis(pheno_file, anno_file,geno_prefix, plinkGenotype, output_dir, int(window_size),
                     min_maf=float(min_maf), min_hwe_P=float(min_hwe_P), min_call_rate=float(min_call_rate), blocksize=int(block_size),
                     cis_mode=cis, gaussianize_method = gaussianize, minimum_test_samples= int(minimum_test_samples), seed=int(random_seed), 
                     n_perm=int(n_perm), write_permutations = write_permutations, relatedness_score=float(relatedness_score), feature_variant_covariate_filename = feature_variant_covariate_filename,
                     snps_filename=snps_filename, feature_filename=feature_filename, snp_feature_filename=snp_feature_filename, genetic_range=genetic_range, covariates_filename=covariates_file,
                     kinship_filename=kinship_file, sample_mapping_filename=samplemap_file, extended_anno_filename=extended_anno_file)
