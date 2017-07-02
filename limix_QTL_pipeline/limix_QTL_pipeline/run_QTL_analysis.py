import pandas as pd
import numpy as np
import limix
import qtl_output
import qtl_loader_utils
from qtl_snp_qc import do_snp_qc
import glob
import os
from sklearn.preprocessing import Imputer
import argparse
import math

#V0.1

def get_args():
    parser = argparse.ArgumentParser(description='Run QTL analysis given genotype, phenotype, and annotation.')
    parser.add_argument('-bgen','--bgen',required=False)
    parser.add_argument('-plink','--plink',required=False)
    parser.add_argument('-anno_file','--anno_file', required=True)
    parser.add_argument('-pheno_file','--pheno_file', required=True)
    parser.add_argument('-output_dir','--output_dir', required=True)
    parser.add_argument('-window','--window', required=True,
                        help=
                        'The size of the cis window to take SNPs from.'
                        'The window will extend between:                     '
                        '    (feature_start - (window))             '
                        ' and:                                               '
                        '    (feature_end + (window))               ',default=250000)
    parser.add_argument('-chromosome','--chromosome',required=False,default='all')
    parser.add_argument('-covariates_file','--covariates_file',required=False,default=None)
    parser.add_argument('-kinship_file','--kinship_file',required=False,default=None)
    parser.add_argument('-samplemap_file','--samplemap_file',required=False,default=None)
    parser.add_argument('-maf','--maf',required=False,default=0.05)
    parser.add_argument('-hwe','--hwe',required=False,default=0.001)
    parser.add_argument('-cr','--cr',required=False,default=0.95)
    parser.add_argument('-block_size','--block_size',required=False,default=500)
    parser.add_argument('-n_perm','--n_perm',required=False,default=0)
    parser.add_argument('-snps','--snps',required=False,default=None)
    parser.add_argument('-features','--features',required=False,default=None)
    parser.add_argument('-seed','--seed',required=False)
    parser.add_argument('-relatedness_score','--relatedness_score',required=False,default=0.95)
    parser.add_argument('-minimum_test_samples','--minimum_test_samples',
                    help="Force a minimal number of samples to test a phenotype, automaticaly adds number of covariates to this number.",required=False,default=10)
    parser.add_argument("--gaussianize",
                    action="store_true",
                    help="Force normal distribution on phenotypes.", default=False)
    parser.add_argument("--cis",
                        action="store_true",
                        help="Run cis analysis.", default=False)
    parser.add_argument("--trans",
                        action="store_true",
                        help="Run trans analysis.", default=False)

    args = parser.parse_args()

    return args


def run_QTL_analysis(pheno_filename, anno_filename, geno_prefix, plinkGenotype, output_dir, window_size=250000, min_maf=0.05, min_hwe_P=0.001, min_call_rate=0.95, blocksize=1000,
                     cis_mode=True, gaussianize=True, minimum_test_samples= 10, seed=np.random.randint(40000), n_perm=0, relatedness_score=0.95, snps_filename=None, feature_filename=None, chromosome='all',
                     covariates_filename=None, kinship_filename=None, sample_mapping_filename=None):
    '''Core function to take input and run QTL tests on a given chromosome.'''
    #Load input data files & filter for relevant data
    #Load input data filesf
    if(plinkGenotype):
        bim,fam,bed = qtl_loader_utils.get_genotype_data(geno_prefix)
    else :
        geno_prefix+='.bgen'
        print(geno_prefix)
    print("Intersecting data.")
    phenotype_df = qtl_loader_utils.get_phenotype_df(pheno_filename)
    annotation_df = qtl_loader_utils.get_annotation_df(anno_filename)

    individual2sample_df = qtl_loader_utils.get_samplemapping_df(sample_mapping_filename,list(phenotype_df.columns),'iid')
    sample2individual_df = qtl_loader_utils.get_samplemapping_df(sample_mapping_filename,list(phenotype_df.columns),'sample')

    ##Filter first the linking files!
    #Subset linking to relevant genotypes.
    individual2sample_df = individual2sample_df.loc[list(set(individual2sample_df.index) & set(fam.index)),:]
    sample2individual_df = sample2individual_df.loc[sample2individual_df['iid'].map(lambda x: x in list(map(str, fam.index))),:]

    #Subset linking to relevant phenotypes.
    sample2individual_df = sample2individual_df.loc[list(set(sample2individual_df.index) & set(phenotype_df.columns)),:]
    individual2sample_df = individual2sample_df[individual2sample_df['sample'].map(lambda x: x in list(map(str, phenotype_df.columns)))]

    #Subset linking vs kinship.
    kinship_df = qtl_loader_utils.get_kinship_df(kinship_filename)
    if kinship_df is not None:
        #Filter from individual2sample_df & sample2individual_df since we don't want to filter from the genotypes.
        individual2sample_df = individual2sample_df.loc[list(set(individual2sample_df.index) & set(kinship_df.index)),:]
        sample2individual_df = sample2individual_df[sample2individual_df['iid'].map(lambda x: x in list(map(str, kinship_df.index)))]

    #Subset linking vs covariates.
    covariate_df = qtl_loader_utils.get_covariate_df(covariates_filename)
    if covariate_df is not None:
        sample2individual_df = sample2individual_df.loc[list(set(sample2individual_df.index) & set(covariate_df.index)),:]
        individual2sample_df = individual2sample_df[individual2sample_df['sample'].map(lambda x: x in list(map(str, covariate_df.index)))]
        minimum_test_samples += covariate_df.shape[1]
    ###

    ##Filter now the actual data!
    #Filter phenotype data based on the linking files.
    phenotype_df = phenotype_df.loc[list(set(phenotype_df.index)&set(annotation_df.index)),sample2individual_df.index.values]

    #Filter kinship data based on the linking files.
    if kinship_df is not None:
        kinship_df = kinship_df.loc[list(set(kinship_df.index)&set(individual2sample_df.index)),list(set(kinship_df.index)&set(individual2sample_df.index))]
        geneticaly_unique_individuals = get_unique_genetic_samples(kinship_df, relatedness_score);

    #Filter covariate data based on the linking files.
    if covariate_df is not None:
        covariate_df = covariate_df.loc[sample2individual_df.index.values,:]
        minimum_test_samples += covariate_df.shape[1]
    ###
    ##Filtering on features and SNPs to test.
    #Do filtering on features.
    feature_filter_df = qtl_loader_utils.get_snp_df(feature_filename)
    if(feature_filter_df is not None):
        phenotype_df = phenotype_df.loc[feature_filter_df.index,:]
    #Prepare to filter on snps.
    snp_filter_df = qtl_loader_utils.get_snp_df(snps_filename)
    ###

    print("Number of samples with genotype & phenotype data: " + str(phenotype_df.shape[1]))
    #Open output files
    qtl_loader_utils.ensure_dir(output_dir)
    output_writer = qtl_output.hdf5_writer(output_dir+'qtl_results_{}.h5'.format(chromosome))

    if(cis_mode):
        #Remove features from the annotation that are on chromosomes which are not present anyway.
        annotation_df = annotation_df = annotation_df[np.in1d(annotation_df['chromosome'],list(set(bim['chrom'])))]
        #Crude filtering for sites on non allosomes.
        annotation_df = annotation_df[annotation_df['chromosome'].map(lambda x: x in list(map(str, range(1, 23))))]
    
    #Determine features to be tested
    if chromosome=='all':
        feature_list = list(set(annotation_df.index)&set(phenotype_df.index))
    else:
        feature_list = list(set(annotation_df[annotation_df['chromosome']==chromosome].index)&set(phenotype_df.index))
    
    #Arrays to store indices of snps tested and pass and fail QC SNPs for features without missingness.
    tested_snp_idxs = []
    pass_qc_snps_all = []
    fail_qc_snps_all = []
    fail_qc_features = []
    #Test features
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
            #Crude filtering for sites on non allosomes.
            snpQuery = snpQuery.loc[snpQuery['chrom'].map(lambda x: x in list(map(str, range(1, 23))))]

        if (len(snpQuery) != 0) and (snp_filter_df is not None):
            snpQuery = snpQuery.loc[snpQuery['snp'].map(lambda x: x in list(map(str, snp_filter_df.index)))]
   
        if len(snpQuery) != 0:
            sample_ids = individual2sample_df.loc[:,'sample'].values
            phenotype_ds = phenotype_df.loc[feature_id,sample_ids]
            contains_missing_samples = any(phenotype_ds.isnull().values)

            if(contains_missing_samples):
                print ('Feature: ' + feature_id + ' contains missing data.')
            phenotype_ds.dropna(inplace=True)
            if phenotype_ds.empty | len(phenotype_ds)<minimum_test_samples :
                print("Feature: "+feature_id+" not tested not enough samples do QTL test.")
                fail_qc_features.append(feature_id)
                continue
            else :
                print ('For, feature: ' + feature_id + ' ' + str(snpQuery.shape[0]) + ' SNPs need to be tested.\n Please stand by.')
            
            if(n_perm!=0):
                bestPermutationPval = np.ones((n_perm), dtype=np.float)
            for snpGroup in chunker(snpQuery, blocksize):
                snp_idxs = snpGroup['i'].values
                snp_names = snpGroup['snp'].values
                
                tested_snp_idxs.extend(snp_idxs)
                #indices for relevant individuals in genotype matrix
                individual_ids = list(set(fam.index)&set(sample2individual_df.loc[phenotype_ds.index,'iid']))
                individual_idxs = fam.loc[individual_ids,'i'].values
                
                #subset genotype matrix, we cannot subselect at the same time, do in two steps.
                snp_df = pd.DataFrame(data=bed[snp_idxs,:].compute().transpose(),index=fam.index,columns=snp_names)
                snp_df = snp_df.loc[individual_ids,:]
                #print('step -1')
                #SNP QC.
                    #Now we do more proper QC on non-identical samples. 
                    #However, we do not use it when checking for missingness. 
                    #That could be extended but gives alot of overhead.
                if not contains_missing_samples:
                    #remove snps from snp_df if they fail QC
                    snp_df = snp_df.loc[:,snp_df.columns[~snp_df.columns.isin(fail_qc_snps_all)]]
                    snps_to_test_df = snp_df.loc[:,snp_df.columns[~snp_df.columns.isin(pass_qc_snps_all)]]
                    
                    #Only do QC on relevant SNPs. join pre-QCed list and new QCed list.
                    if kinship_df is not None and len(geneticaly_unique_individuals)<snps_to_test_df.shape[0]:
                        passed_snp_names,failed_snp_names = do_snp_qc(snps_to_test_df.loc[geneticaly_unique_individuals,:], min_call_rate, min_maf, min_hwe_P)
                    else:
                        passed_snp_names,failed_snp_names = do_snp_qc(snps_to_test_df, min_call_rate, min_maf, min_hwe_P)
                    snps_to_test_df = None
                    
                    #append snp_names and failed_snp_names
                    pass_qc_snps_all.extend(passed_snp_names)
                    fail_qc_snps_all.extend(failed_snp_names)
                    
                    snp_df = snp_df.loc[:,snp_df.columns[snp_df.columns.isin(pass_qc_snps_all)]]
                else:
                    #Do snp QC for relevant section.
                    
                    if kinship_df is not None and len(geneticaly_unique_individuals)<snp_df.shape[0]:
                        passed_snp_names,failed_snp_names = do_snp_qc(snp_df.loc[geneticaly_unique_individuals,:], min_call_rate, min_maf, min_hwe_P)
                    else:
                        passed_snp_names,failed_snp_names = do_snp_qc(snp_df, min_call_rate, min_maf, min_hwe_P)
                    snp_df = snp_df.loc[:,snp_df.columns[snp_df.columns.isin(pass_qc_snps_all)]]
                #print('step 0')
                if len(snp_df.columns) == 0:
                    print("failed: "+''.join(failed_snp_names))
                    print("passed: "+''.join(passed_snp_names))
                    continue
                #We could make use of relatedness when imputing.
                fill_NaN = Imputer(missing_values=np.nan, strategy='mean', axis=0)
                snp_matrix_DF = pd.DataFrame(fill_NaN.fit_transform(snp_df))
                snp_matrix_DF.columns = snp_df.columns
                snp_matrix_DF.index = snp_df.index
                snp_df = None
                
                if kinship_df is not None:
                    kinship_mat = kinship_df.loc[individual_ids,individual_ids].as_matrix()
                else:
                    kinship_mat = None
                
                #map individual_ids to samples
                sample_ids = individual2sample_df.loc[individual_ids,'sample'].values
                phenotype = phenotype_ds.loc[sample_ids].values
                #print('step 1')
                #generate covariate matrix
                if covariate_df is not None:
                    cov_matrix = np.concatenate([np.ones((len(sample_ids),1)),covariate_df.loc[sample_ids,:].values],axis=1)
                else:
                    cov_matrix = None
                if(gaussianize):
                    phenotype = force_normal_distribution(phenotype)
                #fit modelrun
                LMM = limix.qtl.qtl_test_lmm(snp_matrix_DF.values, phenotype,K=kinship_mat,covs=cov_matrix)
                #print('step 2')
                if(n_perm!=0):
                    if kinship_df is not None and len(geneticaly_unique_individuals)<snp_matrix_DF.shape[0]:
                        
                        temp = get_shuffeld_genotypes_preserving_kinship(geneticaly_unique_individuals, relatedness_score, snp_matrix_DF,kinship_df, n_perm)
                        #print(temp.shape)
                        LMM_perm = limix.qtl.qtl_test_lmm(temp, phenotype,K=kinship_mat,covs=cov_matrix)
                        perm = 0;
                        for relevantOutput in chunker(LMM_perm.getPv()[0],snp_matrix_DF.shape[1]) :
                            if(bestPermutationPval[perm] > min(relevantOutput)): 
                                bestPermutationPval[perm] = min(relevantOutput)
                            perm+=1
                    else : 
                        u_snp_matrix = snp_matrix_DF
                        #Shuffle genotypes.
                        index_samples = np.arange(u_snp_matrix.shape[0])
                        np.random.shuffle(index_samples)
                        temp = snp_matrix_DF.iloc[index_samples,:].values
                        for perm_id in range(1,n_perm) :
                            np.random.shuffle(index_samples)
                            temp = np.concatenate((temp, snp_matrix_DF.iloc[index_samples,:].values),axis=1)
                        #print(temp.shape)
                        LMM_perm = limix.qtl.qtl_test_lmm(temp, phenotype,K=kinship_mat,covs=cov_matrix)
                        perm = 0;
                        for relevantOutput in chunker(LMM_perm.getPv()[0],snp_matrix_DF.shape[1]) :
                            if(bestPermutationPval[perm] > min(relevantOutput)): 
                                bestPermutationPval[perm] = min(relevantOutput)
                            perm+=1
                #print('step 3')
                #add these results to qtl_results
                temp_df = pd.DataFrame(index = range(len(snp_matrix_DF.columns)),columns=['feature_id','snp_id','p_value','beta','n_samples','corr_p_value'])
                temp_df['snp_id'] = snp_matrix_DF.columns
                temp_df['feature_id'] = feature_id
                temp_df['beta'] = LMM.getBetaSNP()[0]
                temp_df['p_value'] = LMM.getPv()[0]
                temp_df['n_samples'] = sum(~np.isnan(phenotype))
                #insert default dummy value
                temp_df['corr_p_value'] = -1.0
                if not temp_df.empty :
                    data_written = True
                    output_writer.add_result_df(temp_df)
                #print('step 4')
            #This we need to change in the written file.
        if(n_perm!=0 and data_written):
            #updated_permuted_p_in_hdf5(bestPermutationPval, feature_id);
            output_writer.apply_pval_correction(feature_id,bestPermutationPval)
        else :
            fail_qc_features.append(feature_id)
        print('step 5')
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
    feature_list = [x for x in feature_list if x not in fail_qc_features]
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

def chunker(seq, size):
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))

def get_unique_genetic_samples(kinship_df, identityScore):
    kinship_df_copy = kinship_df.copy(deep=True)
    kinship_df_copy.values[[np.arange(kinship_df.shape[0])]*2] = 0

    processMatrix = True
    s_index=0
    while processMatrix is True:
        selection = kinship_df_copy.iloc[s_index,].values>=identityScore
        #print(selection)
        if(selection.sum()>0):  
            selection_names = kinship_df_copy.columns[~selection]
            kinship_df_copy = kinship_df_copy.loc[selection_names,selection_names]
        s_index+=1
        if(s_index==kinship_df_copy.shape[0]):
            processMatrix = False
    return(kinship_df_copy.columns)

def force_normal_distribution(phenotype, method='gaussnorm', reference=None):
    _doc='rank transform x into ref/ gaussian;keep the range; keep ties'
    
    indextoupdate=np.isfinite(phenotype);
    y1=phenotype[indextoupdate]
    yuni,yindex=np.unique(y1, return_inverse=True)
    if method=='gaussnorm':
        std=np.nanstd(y1)
        mean=np.nanmean(y1)
        sref = scst.norm.isf(np.linspace(scst.norm.cdf((np.nanmin(y1)-mean)/std), scst.norm.cdf((np.nanmax(y1)-mean)/std),num=yuni.shape[0])[::-1])
        
    elif method=='ranknorm':
        try:
            sref=np.sort(reference[np.isfinite(reference)])[np.linspace(0,reference.shape[0]-0.001, num=yuni.shape[0]).astype(int)]
        except: print ('reference missing. provide reference to force_normal_distribution or choose gaussnorm')
    phenotypenorm=phenotype.copy()
    phenotypenorm[indextoupdate]=sref[yindex]
    
    return phenotypenorm

def get_shuffeld_genotypes_preserving_kinship(geneticaly_unique_individuals, identityScore, snp_matrix_DF,kinship_df, n_perm):
    u_snp_matrix = snp_matrix_DF.loc[geneticaly_unique_individuals,:]
        
    #Shuffle and reinflate
    index_samples = np.arange(u_snp_matrix.shape[0])
    locationBuffer = []
    #Prepare location search for permuted snp_matrix_df.
    first_entry = True
    for current_name in geneticaly_unique_individuals :
        selection = kinship_df.loc[current_name].values>=identityScore
        locationBuffer += np.where(selection)
    
    snp_matrix_DF_copy = np.zeros_like(snp_matrix_DF)
    for perm_id in range(0,n_perm) :
        np.random.shuffle(index_samples)
        u_snp_matrix.index = u_snp_matrix.index[index_samples]
        
        #Re-flate genotype matrix
        for index in np.arange(len(geneticaly_unique_individuals)) :
            snp_matrix_DF_copy[locationBuffer[index],] = u_snp_matrix.loc[geneticaly_unique_individuals[index]].values
        if perm_id != 0:
            temp = np.concatenate((temp, snp_matrix_DF_copy),axis=1)
        else : 
            temp = snp_matrix_DF_copy
    return(temp)

if __name__=='__main__':
    args = get_args()
    plink  = args.plink
    bgen = args.bgen
    anno_file = args.anno_file
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
    gaussianize = args.gaussianize
    cis = args.cis
    trans = args.trans

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

    
    run_QTL_analysis(pheno_file, anno_file,geno_prefix, plinkGenotype, output_dir, int(window_size),
                     min_maf=float(min_maf), min_hwe_P=float(min_hwe_P), min_call_rate=float(min_call_rate), blocksize=int(block_size), 
                     cis_mode=cis, gaussianize = gaussianize, minimum_test_samples= 10, seed=int(random_seed), n_perm=int(n_perm), relatedness_score=float(relatedness_score), 
                     snps_filename=snps_filename, feature_filename=feature_filename, chromosome=chromosome, covariates_filename=covariates_file,
                     kinship_filename=kinship_file, sample_mapping_filename=samplemap_file)
