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

def get_args():
    parser = argparse.ArgumentParser(description='Run QTL analysis given genotype, phenotype, and annotation.')
    parser.add_argument('-bgen','--bgen',required=False)
    parser.add_argument('-plink','--plink',required=False)
    parser.add_argument('-anno_file','--anno_file', required=True)
    parser.add_argument('-pheno_file','--pheno_file', required=True)
    parser.add_argument('-output_dir','--output_dir', required=True)
    parser.add_argument('-interaction_terms','--interaction_terms',required=True,default=None)
    parser.add_argument('-window','--window', required=False,
                        help=
                        'The size of the cis window to take SNPs from.'
                        'The window will extend between:                     '
                        '    (feature_start - (window))             '
                        ' and:                                               '
                        '    (feature_end + (window))               ',default=250000)
    parser.add_argument('-chromosome','--chromosome',required=False,
                        help=
                        'The chromosome or selection on a chromsome to make during analysis. Selection is based on features.'
                        'Either a chromsome or chromosome:start-end is accepted.',default='all')
    parser.add_argument('-covariates_file','--covariates_file',required=False,default=None)
    parser.add_argument('-kinship_file','--kinship_file',required=False,default=None)
    parser.add_argument('-samplemap_file','--samplemap_file',required=False,default=None)
    parser.add_argument('-maf','--maf',required=False,default=0.05)
    parser.add_argument('-hwe','--hwe',required=False,default=0.0001)
    parser.add_argument('-cr','--cr',required=False,default=0.95)
    parser.add_argument('-block_size','--block_size',required=False,default=500)
    parser.add_argument('-n_perm','--n_perm',required=False,default=0)
    parser.add_argument('-snps','--snps',required=False,default=None)
    parser.add_argument('-features','--features',required=False,default=None)
    parser.add_argument('-seed','--seed',required=False)
    parser.add_argument('-extended_anno_file','--extended_anno_file'
                        help=
                        'Secondary annotation file, to add a multiple locations to one feature.'
                        'This can be used to either link multiple test regions to one feature or exclude multiple regions while testing a feature.', required=False)
    parser.add_argument('-relatedness_score','--relatedness_score',required=False,default=0.95)
    parser.add_argument('-write_permutations','--write_permutations',action="store_true",required=False,default=False)
    parser.add_argument('-minimum_test_samples','--minimum_test_samples',
                    help="Force a minimal number of samples to test a phenotype, automaticaly adds number of covariates to this number.",required=False,default=10)
    parser.add_argument("--gaussianize_method",
                    help="Force normal distribution on phenotypes.", default=None)
    parser.add_argument("--cis",
                        action="store_true",
                        help="Run cis analysis.", default=False)
    parser.add_argument("--trans",
                        action="store_true",
                        help="Run trans analysis.", default=False)

    args = parser.parse_args()

    return args

def run_interaction_QTL_analysis(pheno_filename, anno_filename, geno_prefix, plinkGenotype, output_dir, interaction_terms, window_size=250000, min_maf=0.05, min_hwe_P=0.001, min_call_rate=0.95, blocksize=1000,
                     cis_mode=True, gaussianize_method=None, minimum_test_samples= 10, seed=np.random.randint(40000), n_perm=0, write_permutations = False, relatedness_score=0.95, snps_filename=None, feature_filename=None, chromosome='all',
                     covariates_filename=None, kinship_filename=None, sample_mapping_filename=None, extended_anno_filename=None):
    '''Core function to take input and run QTL tests on a given chromosome.'''
    
    [phenotype_df, kinship_df, covariate_df, sample2individual_df,annotation_df,snp_filter_df, geneticaly_unique_individuals, minimum_test_samples, feature_list,bim,fam,bed, chromosome, selectionStart, selectionEnd]=\
    utils.run_QTL_analysis_load_intersect_phenotype_covariates_kinship_sample_mapping(pheno_filename=pheno_filename, anno_filename=anno_filename, geno_prefix=geno_prefix, plinkGenotype=plinkGenotype, cis_mode=cis_mode,
                      minimum_test_samples= minimum_test_samples,  relatedness_score=relatedness_score, snps_filename=snps_filename, feature_filename=feature_filename, chromosome=chromosome,
                     covariates_filename=covariates_filename, kinship_filename=kinship_filename, sample_mapping_filename=sample_mapping_filename, extended_anno_filename=extended_anno_filename)
    #Open output files
    qtl_loader_utils.ensure_dir(output_dir)
    if not selectionStart is None :
        output_writer = qtl_output.hdf5_writer(output_dir+'qtl_results_{}_{}_{}.h5'.format(chromosome,selectionStart,selectionEnd))
    else :
        output_writer = qtl_output.hdf5_writer(output_dir+'qtl_results_{}.h5'.format(chromosome))
    if(write_permutations):
        if not selectionStart is None :
            permutation_writer = qtl_output.hdf5_writer(output_dir+'perm_results_{}_{}_{}.h5'.format(chromosome,selectionStart,selectionEnd))
        else :
            permutation_writer = qtl_output.hdf5_writer(output_dir+'perm_results_{}.h5'.format(chromosome))
    
    if(not set(interaction_terms).issubset(set(covariate_df.columns))):
        print ('Interaction terms are not found in the covariates')
        sys.exit()
    
    if(write_permutations):
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
        chrom = str(annotation_df.loc[feature_id,'chromosome'])
        start = annotation_df.loc[feature_id,'start']
        end = annotation_df.loc[feature_id,'end']
        #make robust to features specified back-to-front
        lowest = min([start,end])
        highest = max([start,end])
        if (cis_mode) :
            snpQuery = bim.query("chrom == '%s' & pos > %d & pos < %d" % (chrom, lowest-window_size, highest+window_size))
        else :
            snpQuery = bim.query("(chrom == '%s' & (pos < %d | pos > %d))|chrom != '%s'" % (chrom, lowest-window_size, highest+window_size,chrom))
            #Filtering for sites on non allosomes.
            snpQuery = snpQuery.loc[snpQuery['chrom'].map(lambda x: x in list(map(str, range(1, 23))))]

        if (len(snpQuery) != 0) and (snp_filter_df is not None):
            snpQuery = snpQuery.loc[snpQuery['snp'].map(lambda x: x in list(map(str, snp_filter_df.index)))]

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
            else :
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

                    if kinship_df is not None and len(geneticaly_unique_individuals)>1 and len(geneticaly_unique_individuals)<snp_df.shape[0]:
                        passed_snp_names,failed_snp_names = do_snp_qc(snp_df.iloc[np.unique(snp_df.index,return_index=1)[1]].loc[geneticaly_unique_individuals,:], min_call_rate, min_maf, min_hwe_P) 
                    else:
                        passed_snp_names,failed_snp_names = do_snp_qc(snp_df, min_call_rate, min_maf, min_hwe_P)
                    snp_df = snp_df.loc[:,snp_df.columns[snp_df.columns.isin(passed_snp_names)]]
                #print('step 0')
                if len(snp_df.columns) == 0:
                    continue
                #We could make use of relatedness when imputing.
                fill_NaN = Imputer(missing_values=np.nan, strategy='mean', axis=0)
                snp_matrix_DF = pd.DataFrame(fill_NaN.fit_transform(snp_df))
                snp_matrix_DF.columns = snp_df.columns
                snp_matrix_DF.index = snp_df.index
                snp_df = None

                inter = covariate_df.loc[:,interaction_terms]
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
                    ### select discrete covariates with at least 6 lines:
                    if covariate_df is not None :
                        if np.unique(cov_matrix).shape[0]==2:
                            cov_matrix=cov_matrix[:,np.nansum(cov_matrix==1,0)>6]
                    phenotype = utils.force_normal_distribution(phenotype_ds.values,method=gaussianize_method) if gaussianize_method is not None else phenotype_ds.values
                else:
                    print ('there is an issue in mapping phenotypes and genotypes')
                    sys.exit()
                M = np.asarray(cov_matrix, float)
                LMM = limix.qtl.iscan(snp_matrix_DF.values, phenotype, 'normal', inter.values, K=kinship_mat, M=M)
                if(n_perm!=0):
                    if(write_permutations):
                        perm_df = pd.DataFrame(index = range(len(snp_matrix_DF.columns)),columns=['snp_id'] + ['permutation_'+str(x+1) for x in range(n_perm)])
                        perm_df['snp_id'] = snp_matrix_DF.columns
                    if kinship_df is not None and len(geneticaly_unique_individuals)<snp_matrix_DF.shape[0]:
                        temp = utils.get_shuffeld_genotypes_preserving_kinship(geneticaly_unique_individuals, relatedness_score, snp_matrix_DF,kinship_df.loc[individual_ids,individual_ids], n_perm)
                        LMM_perm = limix.qtl.iscan(temp, phenotype, 'normal', inter.values, K=kinship_mat, M=M)
                        perm = 0;
                        for relevantOutput in utils.chunker(LMM_perm.variant_pvalues.values,snp_matrix_DF.shape[1]) :
                            if(write_permutations):
                                perm_df['permutation_'+str(perm)] = relevantOutput
                            if(bestPermutationPval[perm] > min(relevantOutput)):
                                bestPermutationPval[perm] = min(relevantOutput)
                            perm+=1
                    else :
                        temp = utils.get_shuffeld_genotypes(snp_matrix_DF,kinship_df, n_perm)
                        LMM_perm = limix.qtl.iscan(temp, phenotype, 'normal', inter.values, K=kinship_mat, M=M)
                        perm = 0;
                        for relevantOutput in utils.chunker(LMM_perm.variant_pvalues.values,snp_matrix_DF.shape[1]) :
                            if(write_permutations):
                                perm_df['permutation_'+str(perm)] = relevantOutput
                            if(bestPermutationPval[perm] > min(relevantOutput)):
                                bestPermutationPval[perm] = min(relevantOutput)
                            perm+=1
                #print('step 3')
                #add these results to qtl_results
                temp_df = pd.DataFrame(index = range(len(snp_matrix_DF.columns)),columns=['feature_id','snp_id','p_value','beta','n_samples','corr_p_value'])
                temp_df['snp_id'] = snp_matrix_DF.columns
                temp_df['feature_id'] = feature_id
                #temp_df['beta'] = LMM.getBetaSNP()[0]
                temp_df['beta'] = LMM.variant_effsizes.values
                #temp_df['p_value'] = LMM.getPv()[0]
                temp_df['p_value'] = LMM.variant_pvalues.values
                temp_df['n_samples'] = sum(~np.isnan(phenotype))
                #temp_df['beta_se'] = LMM.getBetaSNPste()[0]
                temp_df['beta_se'] = LMM.variant_effsizes_se.values
                #insert default dummy value
                temp_df['corr_p_value'] = -1.0
                if not temp_df.empty :
                    data_written = True
                    output_writer.add_result_df(temp_df)
                    if(write_permutations):
                        permutation_writer.add_permutation_results_df(perm_df,feature_id)
                if contains_missing_samples:
                    geneticaly_unique_individuals = tmp_unique_individuals
                #print('step 4')
            #This we need to change in the written file.
        if(n_perm!=0 and data_written):
            #updated_permuted_p_in_hdf5(bestPermutationPval, feature_id);
            output_writer.apply_pval_correction(feature_id,bestPermutationPval)
        else :
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

    snp_df.ix[tested_snp_idxs,:].to_csv(output_dir+'/snp_metadata_{}.txt'.format(chromosome),sep='\t',index=False)
    feature_list = [x for x in feature_list if x not in fail_qc_features]
    annotation_df.loc[feature_list,:].to_csv(output_dir+'/feature_metadata_{}.txt'.format(chromosome),sep='\t')

if __name__=='__main__':
    args = get_args()
    plink  = args.plink
    bgen = args.bgen
    anno_file = args.anno_file
    extended_anno_file = args.extended_anno_file
    pheno_file = args.pheno_file
    output_dir = args.output_dir
    window_size = args.window
    chromosome = args.chromosome
    covariates_file = args.covariates_file
    kinship_file = args.kinship_file
    samplemap_file = args.samplemap_file
    min_maf = args.maf
    min_hwe_P = args.hwe
    min_call_rate = args.cr
    block_size = args.block_size
    n_perm = args.n_perm
    snps_filename = args.snps
    random_seed = args.seed
    feature_filename = args.features
    relatedness_score = args.relatedness_score
    minimum_test_samples = args.minimum_test_samples
    gaussianize = args.gaussianize_method
    cis = args.cis
    trans = args.trans
    interaction_terms = args.interaction_terms
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

    run_interaction_QTL_analysis(pheno_file, anno_file,geno_prefix, plinkGenotype, output_dir, interaction_terms, int(window_size),
                     min_maf=float(min_maf), min_hwe_P=float(min_hwe_P), min_call_rate=float(min_call_rate), blocksize=int(block_size),
                     cis_mode=cis, gaussianize_method = gaussianize, minimum_test_samples= int(minimum_test_samples), seed=int(random_seed), n_perm=int(n_perm), write_permutations = write_permutations, relatedness_score=float(relatedness_score),
                     snps_filename=snps_filename, feature_filename=feature_filename, chromosome=chromosome, covariates_filename=covariates_file,
                     kinship_filename=kinship_file, sample_mapping_filename=samplemap_file, extended_anno_filename=extended_anno_file)
