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
                        '    (feature_end + (window))               ')
    parser.add_argument('-chromosome','--chromosome',required=False,default='all')
    parser.add_argument('-covariates_file','--covariates_file',required=False,default=None)
    parser.add_argument('-kinship_file','--kinship_file',required=False,default=None)
    parser.add_argument('-samplemap_file','--samplemap_file',required=False,default=None)
    parser.add_argument('-maf','--maf',required=False,default=0.05)
    parser.add_argument('-hwe','--hwe',required=False,default=0.001)
    parser.add_argument('-cr','--cr',required=False,default=0.95)
    parser.add_argument('-block_size','--block_size',required=False,default=1000)
    parser.add_argument('-n_perm','--n_perm',required=False,default=0)
    parser.add_argument('-snps','--snps',required=False,default=None)
    parser.add_argument('-features','--features',required=False,default=None)
    parser.add_argument('-seed','--seed',required=False)
    parser.add_argument("--cis",
                        action="store_true",
                        help="Run cis analysis.", default=False)
    parser.add_argument("--trans",
                        action="store_true",
                        help="Run trans analysis.", default=False)

    args = parser.parse_args()

    return args


def run_QTL_analysis(pheno_filename,anno_filename,geno_prefix,plinkGenotype,window_size,output_dir, min_maf, min_hwe_P,min_call_rate,blocksize,cis_mode,seed,n_perm=0,snps_filename=None,feature_filename=None,chromosome='all',
                     covariates_filename=None,kinship_filename=None,sample_mapping_filename=None):
    '''Core function to take input and run QTL tests on a given chromosome.'''
    #Load input data files & filter for relevant data
    #Load input data filesf
    if(plinkGenotype):
        bim,fam,bed = qtl_loader_utils.get_genotype_data(geno_prefix)
    else :
        geno_prefix+='.bgen'
        print(geno_prefix)
    print("Intersecting data.")
    feature_filter_df = qtl_loader_utils.get_snp_df(feature_filename)
    phenotype_df = qtl_loader_utils.get_phenotype_df(pheno_filename)
    annotation_df = qtl_loader_utils.get_annotation_df(anno_filename)

    individual2sample_df = qtl_loader_utils.get_samplemapping_df(sample_mapping_filename,list(phenotype_df.columns),'iid')
    sample2individual_df = qtl_loader_utils.get_samplemapping_df(sample_mapping_filename,list(phenotype_df.columns),'sample')

    kinship_df = qtl_loader_utils.get_kinship_df(kinship_filename) 
    if kinship_df is not None:
        kinship_df = kinship_df.loc[list(set(fam.index)&set(individual2sample_df.index)),list(set(fam.index)&set(individual2sample_df.index))]
        #Filter from individual2sample_df & sample2individual_df since we don't want to filter from the genotypes.
        individual2sample_df = individual2sample_df.loc[kinship_df.index,:]
        sample2individual_df = sample2individual_df[sample2individual_df['iid'].map(lambda x: x in list(map(str, kinship_df.index)))]
        
    covariate_df = qtl_loader_utils.get_covariate_df(covariates_filename)
    phenotype_df = phenotype_df.loc[annotation_df.index.values,individual2sample_df.loc[list(set(fam.index)&set(individual2sample_df.index)),'sample'].values]
    
    if(feature_filter_df is not None):
        phenotype_df = phenotype_df.loc[feature_filter_df.index,:]
    
    if covariate_df is not None:
        phenotype_df = phenotype_df.loc[:,covariate_df.index]
        covariate_df = covariate_df.loc[phenotype_df.columns,:]
    snp_filter_df = qtl_loader_utils.get_snp_df(snps_filename)
    
    print("Number of samples with genotype & phenotype data: " + str(phenotype_df.shape[0]))
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
            print ('For, feature: ' + feature_id + ' ' + str(snpQuery.shape[0]) + ' SNPs need to be tested.\n Please stand by.')
            if(n_perm!=0):
                bestPermutationPval = np.ones((n_perm), dtype=np.float)
                actualPval = []
            for snpGroup in chunker(snpQuery, blocksize):
                snp_idxs = snpGroup['i'].values
                snp_names = snpGroup['snp'].values
                
                tested_snp_idxs.extend(snp_idxs)
                
                sample_ids = individual2sample_df.loc[:,'sample'].values

                phenotype_ds = phenotype_df.loc[feature_id,sample_ids]
                contains_missing_samples = any(phenotype_ds.isnull().values)
                phenotype_ds.dropna(inplace=True)
                if phenotype_ds.empty :
                    fail_qc_features.append(feature_id)
                    continue
                #indices for relevant individuals in genotype matrix
                individual_ids = list(set(fam.index)&set(sample2individual_df.loc[phenotype_ds.index,'iid']))
                individual_idxs = fam.loc[individual_ids,'i'].values
                
                #subset genotype matrix, we cannot subselect at the same time, do in two steps.
                snp_df = pd.DataFrame(data=bed[snp_idxs,:].compute().transpose(),index=fam.index,columns=snp_names)
                snp_df = snp_df.loc[individual_ids,:]
                
                #SNP QC.
                if not contains_missing_samples:
                    #remove snps from snp_df if they fail QC
                    snp_df = snp_df.loc[:,snp_df.columns[~snp_df.columns.isin(fail_qc_snps_all)]]
                    snps_to_test_df = snp_df.loc[:,snp_df.columns[~snp_df.columns.isin(pass_qc_snps_all)]]
                    
                    #Only do QC on relevant SNPs. join pre-QCed list and new QCed list.
                    passed_snp_names,failed_snp_names = do_snp_qc(snps_to_test_df, min_call_rate, min_maf, min_hwe_P)
                    snps_to_test_df = None
                    
                    #append snp_names and failed_snp_names
                    pass_qc_snps_all.extend(passed_snp_names)
                    fail_qc_snps_all.extend(failed_snp_names)
                    
                    snp_df = snp_df.loc[:,snp_df.columns[snp_df.columns.isin(pass_qc_snps_all)]]
                else:
                    #Do snp QC for relevant section.
                    passed_snp_names,failed_snp_names = do_snp_qc(snp_df, min_call_rate, min_maf, min_hwe_P)
                    snp_df = snp_df.loc[:,snp_df.columns[snp_df.columns.isin(pass_qc_snps_all)]]
                
                if len(snp_df.columns) == 0:
                    continue
                snp_matrix = snp_df.values
                snp_names = snp_df.columns
                snp_df = None
                snp_matrix = Imputer(missing_values=np.nan,strategy='mean',axis=0,copy=False).fit_transform(snp_matrix)
                
                if kinship_df is not None:
                    kinship_mat = kinship_df.loc[individual_ids,individual_ids].as_matrix()
                else:
                    kinship_mat = None
                
                #map individual_ids to samples
                sample_ids = individual2sample_df.loc[individual_ids,'sample'].values
                phenotype = phenotype_ds.loc[sample_ids].values

                #generate covariate matrix
                if covariate_df is not None:
                    cov_matrix = np.concatenate([np.ones((len(sample_ids),1)),covariate_df.loc[sample_ids,:].values],axis=1)
                else:
                    cov_matrix = None

                #fit modelrun
                LMM = limix.qtl.qtl_test_lmm(snp_matrix, phenotype,K=kinship_mat,covs=cov_matrix)
                if(n_perm!=0):
                    #countPermutations = np.zeros((snp_matrix.shape[1]), dtype=np.int)
                    #nBetterCorrelation = np.zeros((snp_matrix.shape[1]), dtype=np.int)
                        #If we are going to need another permutation based on breaking the genotypes (Not the phenotypes).
                        #In this case we want to add the snps to one matrix, to speed up the tests and calculate later which come from what permutaiton.
                    np.random.seed(seed)
                    for perm in range(0,n_perm) :
                        LMM_perm = limix.qtl.qtl_test_lmm(snp_matrix, np.random.permutation(phenotype),K=kinship_mat,covs=cov_matrix)
                        #print(np.random.permutation(phenotype))
                        if(bestPermutationPval[perm] > min(LMM_perm.getPv()[0])):
                                bestPermutationPval[perm] = min(LMM_perm.getPv()[0])
                        #print(bestPermutationPval[perm])
                        #for snp in range(0,len(LMM_perm.getPv()[0])):
                            #print(LMM_perm.getPv()[0][snp])
                            #countPermutations[snp]+=1;
                            #if(bestPermutationPval[perm]<=LMM.getPv()[0][snp]):
                                #nBetterCorrelation[snp]+=1

                #add these results to qtl_results
                temp_df = pd.DataFrame(index = range(len(snp_names)),columns=['feature_id','snp_id','p_value','beta','n_samples','corr_p_value'])
                temp_df['snp_id'] = snp_names
                temp_df['feature_id'] = feature_id
                temp_df['beta'] = LMM.getBetaSNP()[0]
                temp_df['p_value'] = LMM.getPv()[0]
                if(n_perm!=0):
                    actualPval.extend(list(LMM.getPv()[0]))
                temp_df['n_samples'] = sum(~np.isnan(phenotype))
                #insert default dummy value
                temp_df['corr_p_value'] = 1.0
                if not temp_df.empty :
                    data_written = True
                    output_writer.add_result_df(temp_df)
            #This we need to change in the written file.
        if(n_perm!=0 and data_written):
            #updated_permuted_p_in_hdf5(bestPermutationPval, feature_id);
            output_writer.apply_pval_correction(feature_id,bestPermutationPval)
        else :
            fail_qc_features.append(feature_id)
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
    if (seed is None):
        seed = np.random.randint(40000)
    
    run_QTL_analysis(pheno_file,anno_file,geno_prefix,plinkGenotype,int(window_size),output_dir,
                     min_maf=float(min_maf), min_hwe_P=float(min_hwe_P),
                     min_call_rate=float(min_call_rate),blocksize=int(block_size), cis_mode=cis,
                     seed=int(random_seed), n_perm=int(n_perm), snps_filename=snps_filename, feature_filename=feature_filename,
                     chromosome=chromosome, covariates_filename=covariates_file,
                     kinship_filename=kinship_file,
                     sample_mapping_filename=samplemap_file)
