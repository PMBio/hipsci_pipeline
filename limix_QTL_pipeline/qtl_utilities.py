import sys
import gc
import numpy as np
import pandas as pd
import scipy.stats as scst
import scipy as sp
import scipy.linalg as la
from bgen_reader import read_bgen
import qtl_loader_utils
import pdb

def run_QTL_analysis_load_intersect_phenotype_covariates_kinship_sample_mapping\
    (pheno_filename, anno_filename, geno_prefix, plinkGenotype,minimum_test_samples= 10, relatedness_score=0.95,cis_mode=True, skipAutosomeFiltering = False, snps_filename=None,
         feature_filename=None, snp_feature_filename=None, selection='all', covariates_filename=None, kinship_filename=None, sample_mapping_filename=None,
         extended_anno_filename=None, feature_variant_covariate_filename=None):
    selectionStart = None
    selectionEnd = None
    if(":" in selection):
        parts = selection.split(":")
        if("-" not in parts[1]):
            print("No correct sub selection.")
            print("Given in: "+selection)
            print("Expected format: (chr number):(start location)-(stop location)")
            sys.exit()
        chromosome = parts[0]
        if("-" in parts[1]):
            parts2 = parts[1].split("-") 
            selectionStart = int(parts2[0])
            selectionEnd = int(parts2[1])
    else :
        chromosome=selection

    ''' function to take input and intersect sample and genotype.'''
    #Load input data files & filter for relevant data
    #Load input data filesf

    phenotype_df = qtl_loader_utils.get_phenotype_df(pheno_filename)
    annotation_df = qtl_loader_utils.get_annotation_df(anno_filename)

    phenotype_df.columns = phenotype_df.columns.astype("str")
    phenotype_df.index = phenotype_df.index.astype("str")
    annotation_df.columns = annotation_df.columns.astype("str")
    annotation_df.index = annotation_df.index.astype("str")

    if(plinkGenotype):
        bim,fam,bed = qtl_loader_utils.get_genotype_data(geno_prefix)
        bgen=None
    else :
        bgen = read_bgen(geno_prefix+'.bgen', verbose=False)
        bim = bgen['variants']
        bim = bim.assign(i = range(bim.shape[0]))
        bim = bim.rename(index=str, columns={"id": "snp"})
        bim['a1'] = bim['allele_ids'].str.split(",", expand=True)[0]
        fam = bgen['samples']
        fam = fam.rename(index=str, columns={"id": "iid"})
        fam.index=fam["iid"]
        bed=None
        print("Warning, the current software only supports biallelic SNPs and plody 2")
        bim['chrom'].replace(['X', 'Y', 'XY', 'MT'], ['23', '24', '25', '26'],inplace=True)
        bim.loc[np.logical_and(bim['nalleles']<3,bim['nalleles']>0),:]

    annotation_df.replace(['X', 'Y', 'XY', 'MT'], ['23', '24', '25', '26'],inplace=True)
    if chromosome=='X' :
        chromosome = '23'
    elif chromosome=='Y':
        chromosome = '24'
    elif chromosome=='XY':
        chromosome='25'
    elif chromosome=='MT':
        chromosome='26'

    print("Intersecting data.")
    if(annotation_df.shape[0] != annotation_df.groupby(annotation_df.index).first().shape[0]): 
        print("Only one location per feature supported. If multiple locations are needed please look at: --extended_anno_file")
        sys.exit()

    ##Make sure that there is only one entry per feature id!.

    sample2individual_df = qtl_loader_utils.get_samplemapping_df(sample_mapping_filename,list(phenotype_df.columns),'sample')
    sample2individual_df.index = sample2individual_df.index.astype('str')
    sample2individual_df = sample2individual_df.astype('str')
    sample2individual_df['sample']=sample2individual_df.index
    sample2individual_df = sample2individual_df.drop_duplicates();
    ##Filter first the linking files!
    #Subset linking to relevant genotypes.
    orgSize = sample2individual_df.shape[0]
    sample2individual_df = sample2individual_df.loc[sample2individual_df['iid'].map(lambda x: x in list(map(str, fam.index))),:]
    diff = orgSize- sample2individual_df.shape[0]
    orgSize = sample2individual_df.shape[0]
    print("Dropped: "+str(diff)+" samples because they are not present in the genotype file.")
    
    #Subset linking to relevant phenotypes.
    sample2individual_df = sample2individual_df.loc[np.intersect1d(sample2individual_df.index,phenotype_df.columns),:]
    diff = orgSize- sample2individual_df.shape[0]
    orgSize = sample2individual_df.shape[0]
    print("Dropped: "+str(diff)+" samples because they are not present in the phenotype file.")
    #Subset linking vs kinship.
    kinship_df = qtl_loader_utils.get_kinship_df(kinship_filename)
    if kinship_df is not None:
        #Filter from individual2sample_df & sample2individual_df since we don't want to filter from the genotypes.
        sample2individual_df = sample2individual_df[sample2individual_df['iid'].map(lambda x: x in list(map(str, kinship_df.index)))]
        diff = orgSize- sample2individual_df.shape[0]
        orgSize = sample2individual_df.shape[0]
        print("Dropped: "+str(diff)+" samples because they are not present in the kinship file.")
    #Subset linking vs covariates.
    covariate_df = qtl_loader_utils.get_covariate_df(covariates_filename)
    if covariate_df is not None:
        if np.nansum(covariate_df==1,0).max()<covariate_df.shape[0]: covariate_df.insert(0, 'ones', np.ones(covariate_df.shape[0]))
        sample2individual_df = sample2individual_df.loc[list(set(sample2individual_df.index) & set(covariate_df.index)),:]
        diff = orgSize- sample2individual_df.shape[0]
        orgSize = sample2individual_df.shape[0]
        print("Dropped: "+str(diff)+" samples because they are not present in the covariate file.")

    ###
    print("Number of samples with genotype & phenotype data: " + str(sample2individual_df.shape[0]))
    if(sample2individual_df.shape[0]<minimum_test_samples):
        print("Not enough samples with both genotype & phenotype data.")
        sys.exit()

    ##Filter now the actual data!
    #Filter phenotype data based on the linking files.
    phenotype_df = phenotype_df.loc[list(set(phenotype_df.index)&set(annotation_df.index)),sample2individual_df.index.values]

    #Filter kinship data based on the linking files.
    genetically_unique_individuals = None
    if kinship_df is not None:
        kinship_df = kinship_df.loc[np.intersect1d(kinship_df.index,sample2individual_df['iid']),np.intersect1d(kinship_df.index,sample2individual_df['iid'])]
        genetically_unique_individuals = get_unique_genetic_samples(kinship_df, relatedness_score);

    #Filter covariate data based on the linking files.
    
    snp_feature_filter_df= qtl_loader_utils.get_snp_feature_df(snp_feature_filename)
    try:
        feature_filter_df = qtl_loader_utils.get_snp_df(feature_filename)
    except:
        if feature_filename  is not None:
            feature_filter_df=pd.DataFrame(index=feature_filename)
    #Do filtering on features.
    if feature_filter_df is not None:
        phenotype_df = phenotype_df.loc[feature_filter_df.index,:]
        ##Filtering on features to test.
    if snp_feature_filter_df is not None:
        phenotype_df = phenotype_df.loc[np.unique(snp_feature_filter_df['feature']),:]
        ##Filtering on features  to test from the combined feature snp filter.

    if ((not cis_mode) and len(set(bim['chrom']))<22) :
        print("Warning, running a trans-analysis on snp data from less than 22 chromosomes.\nTo merge data later the permutation P-values need to be written out.")

    if(cis_mode):
        #Remove features from the annotation that are on chromosomes which are not present anyway.
        annotation_df = annotation_df[np.in1d(annotation_df['chromosome'],list(set(bim['chrom'])))]

    #Prepare to filter on snps.
    snp_filter_df = qtl_loader_utils.get_snp_df(snps_filename)
    if snp_filter_df is not None:
        toSelect = set(snp_filter_df.index).intersection(set(bim['snp']))
        bim = bim.loc[bim['snp'].isin(toSelect)]
        ##Filtering on SNPs to test from the snp filter.

    if snp_feature_filter_df is not None:
        toSelect = set(np.unique(snp_feature_filter_df['snp_id'])).intersection(set(bim['snp']))
        bim = bim.loc[bim['snp'].isin(toSelect)]
        ##Filtering on features  to test from the combined feature snp filter.
    
    #Filtering for sites on non allosomes.
    if not skipAutosomeFiltering :
        annotation_df = annotation_df[annotation_df['chromosome'].map(lambda x: x in list(map(str, range(1, 23))))]
    
    #Determine features to be tested
    if chromosome=='all':
        feature_list = list(set(annotation_df.index)&set(phenotype_df.index))
    else:
        if not selectionStart is None :
            lowest = min([selectionStart,selectionEnd])
            highest = max([selectionStart,selectionEnd])
            annotation_df['mean'] = ((annotation_df["start"] + annotation_df["end"])/2)
            feature_list = list(set(annotation_df.iloc[(annotation_df['chromosome'].values==chromosome) & (annotation_df['mean'].values>=lowest) & (annotation_df["mean"].values<highest)].index.values)&set(phenotype_df.index))
            del annotation_df['mean']
        else :
            feature_list = list(set(annotation_df[annotation_df['chromosome']==chromosome].index)&set(phenotype_df.index))
    #Drop not used feature information.
    phenotype_df = phenotype_df.loc[feature_list,:]
    gc.collect()
    print("Number of features to be tested: " + str(len(feature_list)))
    print("Total number of variants to be considered, before variante QC and feature intersection: " + str(bim.shape[0]))
    
    if(phenotype_df.shape[1]<minimum_test_samples):
        print("Not enough samples with both genotype & phenotype data, for current number of covariates.")
        sys.exit()
    
    if extended_anno_filename is not None:
        complete_annotation_df = pd.read_csv(extended_anno_filename,sep='\t',index_col=0)
        annotation_df['index']=annotation_df.index
        complete_annotation_df['index']=complete_annotation_df.index
        complete_annotation_df = pd.concat([annotation_df,complete_annotation_df]).drop_duplicates()
        del complete_annotation_df['index']
    else:
        complete_annotation_df = annotation_df

    feature_variant_covariate_df = qtl_loader_utils.get_snp_feature_df(feature_variant_covariate_filename) 

    return [phenotype_df, kinship_df, covariate_df, sample2individual_df, complete_annotation_df, annotation_df, snp_filter_df, 
        snp_feature_filter_df, genetically_unique_individuals, minimum_test_samples, feature_list, bim, fam, bed, bgen, chromosome, 
        selectionStart, selectionEnd, feature_variant_covariate_df]

def run_structLMM_QTL_analysis_load_intersect_phenotype_environments_covariates_kinship_sample_mapping\
        (pheno_filename, anno_filename, env_filename, geno_prefix, plinkGenotype,  
            cis_mode = True, association_mode = True, skipAutosomeFiltering = False, minimum_test_samples = 10, 
            relatedness_score = 0.95, snps_filename = None, feature_filename = None, 
            snp_feature_filename = None, selection = 'all', covariates_filename = None, kinship_filename = None, 
            sample_mapping_filename = None, extended_anno_filename = None, feature_variant_covariate_filename = None):  
    selectionStart = None
    selectionEnd = None
    if(":" in selection):
        parts = selection.split(":")
        if("-" not in parts[1]):
            print("No correct sub selection.")
            print("Given in: "+selection)
            print("Expected format: (chr number):(start location)-(stop location)")
            sys.exit()
        chromosome = parts[0]
        if("-" in parts[1]):
            parts2 = parts[1].split("-") 
            selectionStart = int(parts2[0])
            selectionEnd = int(parts2[1])
    else :
        chromosome=selection

    ''' function to take input and intersect sample and genotype.'''
    #Load input data files & filter for relevant data
    #Load input data filesf

    phenotype_df = qtl_loader_utils.get_phenotype_df(pheno_filename)
    annotation_df = qtl_loader_utils.get_annotation_df(anno_filename)

    if(plinkGenotype):
        bim,fam,bed = qtl_loader_utils.get_genotype_data(geno_prefix)
        annotation_df.replace(['X', 'Y', 'XY', 'MT'], ['23', '24', '25', '26'],inplace=True)
        if chromosome=='X' :
            chromosome = '23'
        elif chromosome=='Y':
            chromosome = '24'
        elif chromosome=='XY':
            chromosome='25'
        elif chromosome=='MT':
            chromosome='26'
         #X  -> 23
         #Y  -> 24
         #XY -> 25
         #MT -> 26

    else :
        geno_prefix+='.bgen'
        print(geno_prefix)
    print("Intersecting data.")

    if(annotation_df.shape[0] != annotation_df.groupby(annotation_df.index).first().shape[0]): 
        print("Only one location per feature supported. If multiple locations are needed please look at: --extended_anno_file")
        sys.exit()

    ##Make sure that there is only one entry per feature id!.
    sample2individual_df = qtl_loader_utils.get_samplemapping_df(sample_mapping_filename,list(phenotype_df.columns),'sample')
    sample2individual_df['sample']=sample2individual_df.index
    sample2individual_df = sample2individual_df.drop_duplicates();


    ##Filter first the linking files!
    #Subset linking to relevant genotypes.
    orgSize = sample2individual_df.shape[0]
    sample2individual_df = sample2individual_df.loc[sample2individual_df['iid'].map(lambda x: x in list(map(str, fam.index))),:]
    diff = orgSize- sample2individual_df.shape[0]
    orgSize = sample2individual_df.shape[0]
    print("Dropped: "+str(diff)+" samples because they are not present in the genotype file.")
    
    #Subset linking to relevant phenotypes.
    sample2individual_df = sample2individual_df.loc[np.intersect1d(sample2individual_df.index,phenotype_df.columns),:]
    diff = orgSize- sample2individual_df.shape[0]
    orgSize = sample2individual_df.shape[0]
    print("Dropped: "+str(diff)+" samples because they are not present in the phenotype file.")
    #Subset linking vs kinship.
    kinship_df = qtl_loader_utils.get_kinship_df(kinship_filename)
    if kinship_df is not None:
        #Filter from individual2sample_df & sample2individual_df since we don't want to filter from the genotypes.
        sample2individual_df = sample2individual_df[sample2individual_df['iid'].map(lambda x: x in list(map(str, kinship_df.index)))]
        diff = orgSize- sample2individual_df.shape[0]
        orgSize = sample2individual_df.shape[0]
        print("Dropped: "+str(diff)+" samples because they are not present in the kinship file.")
    #Subset linking vs covariates.
    covariate_df = qtl_loader_utils.get_covariate_df(covariates_filename)
    if covariate_df is not None:
        if np.nansum(covariate_df==1,0).max()<covariate_df.shape[0]: covariate_df.insert(0, 'ones',np.ones(covariate_df.shape[0]))
        sample2individual_df = sample2individual_df.loc[list(set(sample2individual_df.index) & set(covariate_df.index)),:]
        diff = orgSize- sample2individual_df.shape[0]
        orgSize = sample2individual_df.shape[0]
        print("Dropped: "+str(diff)+" samples because they are not present in the covariate file.")
    #Subset linking vs environments.
    environment_df = qtl_loader_utils.get_env_df(env_filename)
    if np.nansum(environment_df==1,0).max()<environment_df.shape[0]: environment_df.insert(0, 'ones',np.ones(environment_df.shape[0]))
    sample2individual_df = sample2individual_df.loc[list(set(sample2individual_df.index) & set(environment_df.index)),:]
    diff = orgSize - sample2individual_df.shape[0]
    orgSize = sample2individual_df.shape[0]
    print("Dropped: "+str(diff)+" samples because they are not present in the environment file.")

    ###
    print("Number of samples with genotype & phenotype data: " + str(sample2individual_df.shape[0]))
    if(sample2individual_df.shape[0]<minimum_test_samples):
        print("Not enough samples with both genotype & phenotype data.")
        sys.exit()

    ##Filter now the actual data!
    #Filter phenotype data based on the linking files.
    phenotype_df = phenotype_df.loc[list(set(phenotype_df.index)&set(annotation_df.index)),sample2individual_df.index.values]

    #Filter kinship data based on the linking files.
    genetically_unique_individuals = None
    if kinship_df is not None:
        kinship_df = kinship_df.loc[np.intersect1d(kinship_df.index,sample2individual_df['iid']),np.intersect1d(kinship_df.index,sample2individual_df['iid'])]
        genetically_unique_individuals = get_unique_genetic_samples(kinship_df, relatedness_score);

    #Filter covariate data based on the linking files.
    if covariate_df is not None:
        covariate_df = covariate_df.loc[np.intersect1d(covariate_df.index,sample2individual_df.index.values),:]
    
    snp_feature_filter_df= qtl_loader_utils.get_snp_feature_df(snp_feature_filename)
    try:
        feature_filter_df = qtl_loader_utils.get_snp_df(feature_filename)
    except:
        if feature_filename  is not None:
            feature_filter_df=pd.DataFrame(index=feature_filename)
    #Do filtering on features.
    if feature_filter_df is not None:
        phenotype_df = phenotype_df.loc[feature_filter_df.index,:]
        ##Filtering on features to test.
    if snp_feature_filter_df is not None:
        phenotype_df = phenotype_df.loc[np.unique(snp_feature_filter_df['feature']),:]
        ##Filtering on features  to test from the combined feature snp filter.

    if ((not cis_mode) and len(set(bim['chrom']))<22) :
        print("Warning, running a trans-analysis on snp data from less than 22 chromosomes.\nTo merge data later the permutation P-values need to be written out.")

    if(cis_mode):
        #Remove features from the annotation that are on chromosomes which are not present anyway.
        annotation_df = annotation_df[np.in1d(annotation_df['chromosome'],list(set(bim['chrom'])))]

    #Prepare to filter on snps.
    snp_filter_df = qtl_loader_utils.get_snp_df(snps_filename)
    if snp_filter_df is not None:
        toSelect = set(snp_filter_df.index).intersection(set(bim['snp']))
        bim = bim.loc[bim['snp'].isin(toSelect)]
        ##Filtering on SNPs to test from the snp filter.

    if snp_feature_filter_df is not None:
        toSelect = set(np.unique(snp_feature_filter_df['snp_id'])).intersection(set(bim['snp']))
        bim = bim.loc[bim['snp'].isin(toSelect)]
        ##Filtering on features  to test from the combined feature snp filter.
    
    #Filtering for sites on non allosomes.
    if not skipAutosomeFiltering :
        annotation_df = annotation_df[annotation_df['chromosome'].map(lambda x: x in list(map(str, range(1, 23))))]
    
    #Determine features to be tested
    if chromosome=='all':
        feature_list = list(set(annotation_df.index)&set(phenotype_df.index))
    else:
        if not selectionStart is None :
            lowest = min([selectionStart,selectionEnd])
            highest = max([selectionStart,selectionEnd])
            annotation_df['mean'] = ((annotation_df["start"] + annotation_df["end"])/2)
            feature_list = list(set(annotation_df.iloc[(annotation_df['chromosome'].values==chromosome) & (annotation_df['mean'].values>=lowest) & (annotation_df["mean"].values<highest)].index.values)&set(phenotype_df.index))
            del annotation_df['mean']
        else :
            feature_list = list(set(annotation_df[annotation_df['chromosome']==chromosome].index)&set(phenotype_df.index))

    print("Number of features to be tested: " + str(len(feature_list)))
    print("Total number of variants to be considered, before variante QC and feature intersection: " + str(bim.shape[0]))
    
    if(phenotype_df.shape[1]<minimum_test_samples):
        print("Not enough samples with both genotype & phenotype data, for current number of covariates.")
        sys.exit()
    
    if extended_anno_filename is not None:
        complete_annotation_df = pd.read_csv(extended_anno_filename,sep='\t',index_col=0)
        annotation_df['index']=annotation_df.index
        complete_annotation_df['index']=complete_annotation_df.index
        complete_annotation_df = pd.concat([annotation_df,complete_annotation_df]).drop_duplicates()
        del complete_annotation_df['index']
    else:
        complete_annotation_df = annotation_df

    feature_variant_covariate_df = qtl_loader_utils.get_snp_feature_df(feature_variant_covariate_filename) 

    return [phenotype_df, kinship_df, covariate_df, environment_df, sample2individual_df, complete_annotation_df, annotation_df, snp_filter_df, snp_feature_filter_df, genetically_unique_individuals, minimum_test_samples, feature_list,bim,fam,bed, chromosome, selectionStart, selectionEnd, feature_variant_covariate_df]

def run_PrsQtl_analysis_load_intersect_phenotype_covariates_kinship_sample_mapping\
        (pheno_filename, anno_filename, prsFile, minimum_test_samples= 10, relatedness_score=0.95, skipAutosomeFiltering = False, snps_filename=None,
         feature_filename=None, snp_feature_filename=None, selection='all', covariates_filename=None, kinship_filename=None, sample_mapping_filename=None, feature_variant_covariate_filename=None):
    
    selectionStart = None
    selectionEnd = None
    if(":" in selection):
        parts = selection.split(":")
        if("-" not in parts[1]):
            print("No correct sub selection.")
            print("Given in: "+selection)
            print("Expected format: (chr number):(start location)-(stop location)")
            sys.exit()
        chromosome = parts[0]
        if("-" in parts[1]):
            parts2 = parts[1].split("-") 
            selectionStart = int(parts2[0])
            selectionEnd = int(parts2[1])
    else :
        chromosome=selection

    ''' function to take input and intersect sample and genotype.'''
    #Load input data files & filter for relevant data
    #Load input data filesf
    #import pdb; pdb.set_trace();
    phenotype_df = qtl_loader_utils.get_phenotype_df(pheno_filename)
    annotation_df = qtl_loader_utils.get_annotation_df(anno_filename)

    phenotype_df.columns = phenotype_df.columns.astype("str")
    phenotype_df.index = phenotype_df.index.astype("str")
    annotation_df.columns = annotation_df.columns.astype("str")
    annotation_df.index = annotation_df.index.astype("str")
    
    if(annotation_df.shape[0] != annotation_df.groupby(annotation_df.index).first().shape[0]): 
        print("Only one location per feature supported. If multiple locations are needed please look at: --extended_anno_file")
        sys.exit()
    
    #Determine features to be tested
    if chromosome!='all':
        if not selectionStart is None :
            lowest = min([selectionStart,selectionEnd])
            highest = max([selectionStart,selectionEnd])
            annotation_df['mean'] = ((annotation_df["start"] + annotation_df["end"])/2)
            feature_list = list(set(annotation_df.iloc[(annotation_df['chromosome'].values==chromosome) & (annotation_df['mean'].values>=lowest) & (annotation_df["mean"].values<highest)].index.values))
            annotation_df = annotation_df.loc[feature_list,]
            del annotation_df['mean']
        else :
            feature_list = list(annotation_df[annotation_df['chromosome']==chromosome].index)
            annotation_df = annotation_df.loc[feature_list,]
    
    #To be able to read variants from a large file we change the loading here.
    #First we subset the genes to the chunk and get the relevant SNPs based on that.
    
    snp_feature_filter_df= qtl_loader_utils.get_snp_feature_df(snp_feature_filename)
    feature_filter_df = qtl_loader_utils.get_snp_df(feature_filename)
    snp_filter_df = qtl_loader_utils.get_snp_df(snps_filename)
    
    #import pdb; pdb.set_trace()
    #Do filtering on variants and features first stage.
    if snp_feature_filter_df is not None:
        if feature_filter_df is not None:
            toSelect = set(feature_filter_df.index.values).intersection(set(annotation_df.index.values))
            annotation_df = annotation_df.loc[toSelect,]
        toSelect = list(set(snp_feature_filter_df['feature'].values).intersection(set(annotation_df.index.values)))
        snp_feature_filter_df = snp_feature_filter_df.loc[snp_feature_filter_df['feature'].isin(toSelect)]
        relSnps = np.unique(snp_feature_filter_df['snp_id'].values)
        
        if snp_filter_df is not None:
            relSnps = set(snp_filter_df.index).intersection(set(relSnps))
        risk_df = qtl_loader_utils.get_grs_subset_df(prsFile, relSnps)
        if risk_df is None:
            print("No variants selected during SNP reading.")
            sys.exit()
        risk_df = risk_df.assign(SnpId=risk_df.index.values)
        risk_df = risk_df.drop_duplicates(keep='first')
        risk_df = risk_df.drop(['SnpId'], axis='columns')
        risk_df = risk_df.loc[risk_df.isnull().sum(axis=1)!=risk_df.shape[1],]
    elif snp_filter_df is not None:
        relSnps = unique(snp_filter_df.index)
        risk_df = qtl_loader_utils.get_grs_subset_df(prsFile, relSnps)
        if risk_df is None:
            print("No variants selected during SNP reading.")
            sys.exit()
        risk_df = risk_df.assign(SnpId=risk_df.index.values)
        risk_df = risk_df.drop_duplicates(keep='first')
        risk_df = risk_df.drop(['SnpId'], axis='columns')
        risk_df = risk_df.loc[risk_df.isnull().sum(axis=1)!=risk_df.shape[1],]
    else :
        risk_df = qtl_loader_utils.get_phenotype_df(prsFile)
    print("Intersecting data.")
    risk_df =  risk_df.astype(float)
    #pdb.set_trace();
    ##Make sure that there is only one entry per feature id!.

    sample2individual_df = qtl_loader_utils.get_samplemapping_df(sample_mapping_filename,list(phenotype_df.columns),'sample')
    sample2individual_df.index = sample2individual_df.index.astype('str')
    sample2individual_df = sample2individual_df.astype('str')
    sample2individual_df['sample']=sample2individual_df.index
    sample2individual_df = sample2individual_df.drop_duplicates();
    ##Filter first the linking files!
    #Subset linking to relevant genotypes.
    orgSize = sample2individual_df.shape[0]
    sample2individual_df = sample2individual_df.loc[sample2individual_df['iid'].map(lambda x: x in list(map(str, risk_df.columns))),:]
    diff = orgSize- sample2individual_df.shape[0]
    orgSize = sample2individual_df.shape[0]
    print("Dropped: "+str(diff)+" samples because they are not present in the genotype file.")
    
    #Subset linking to relevant phenotypes.
    sample2individual_df = sample2individual_df.loc[np.intersect1d(sample2individual_df.index,phenotype_df.columns),:]
    diff = orgSize- sample2individual_df.shape[0]
    orgSize = sample2individual_df.shape[0]
    print("Dropped: "+str(diff)+" samples because they are not present in the phenotype file.")
    #Subset linking vs kinship.
    kinship_df = qtl_loader_utils.get_kinship_df(kinship_filename)
    if kinship_df is not None:
        #Filter from individual2sample_df & sample2individual_df since we don't want to filter from the genotypes.
        sample2individual_df = sample2individual_df[sample2individual_df['iid'].map(lambda x: x in list(map(str, kinship_df.index)))]
        diff = orgSize- sample2individual_df.shape[0]
        orgSize = sample2individual_df.shape[0]
        print("Dropped: "+str(diff)+" samples because they are not present in the kinship file.")
    #Subset linking vs covariates.
    covariate_df = qtl_loader_utils.get_covariate_df(covariates_filename)
    if covariate_df is not None:
        if np.nansum(covariate_df==1,0).max()<covariate_df.shape[0]: covariate_df.insert(0, 'ones',np.ones(covariate_df.shape[0]))
        sample2individual_df = sample2individual_df.loc[list(set(sample2individual_df.index) & set(covariate_df.index)),:]
        diff = orgSize- sample2individual_df.shape[0]
        orgSize = sample2individual_df.shape[0]
        print("Dropped: "+str(diff)+" samples because they are not present in the covariate file.")

    ###
    print("Number of samples with genotype & phenotype data: " + str(sample2individual_df.shape[0]))
    if(sample2individual_df.shape[0]<minimum_test_samples):
        print("Not enough samples with both genotype & phenotype data.")
        sys.exit()
    #import pdb; pdb.set_trace()
    ##Filter now the actual data!
    #Filter phenotype data based on the linking files.
    phenotype_df = phenotype_df.loc[list(set(phenotype_df.index)&set(annotation_df.index)),sample2individual_df.index.values]

    #Filter kinship data based on the linking files.
    genetically_unique_individuals = None
    if kinship_df is not None:
        kinship_df = kinship_df.loc[np.intersect1d(kinship_df.index,sample2individual_df['iid']),np.intersect1d(kinship_df.index,sample2individual_df['iid'])]
        genetically_unique_individuals = get_unique_genetic_samples(kinship_df, relatedness_score);

    #Filter covariate data based on the linking files.
    
    #Do filtering on features.
    if feature_filter_df is not None:
        toSelect = set(feature_filter_df.index.values).intersection(set(phenotype_df.index.values))
        phenotype_df = phenotype_df.loc[toSelect,:]
        ##Filtering on features to test.
    if snp_feature_filter_df is not None:
        toSelect = set(snp_feature_filter_df['feature'].values).intersection(set(phenotype_df.index.values))
        phenotype_df = phenotype_df.loc[toSelect,:]
        if feature_filter_df is not None:
            snp_feature_filter_df = snp_feature_filter_df.loc[snp_feature_filter_df['feature'].isin(toSelect)]
        ##Filtering on features  to test from the combined feature snp filter.

    #Prepare to filter on SNPs.
    if snp_filter_df is not None:
        toSelect = set(snp_filter_df.index).intersection(set(risk_df.index.values))
        risk_df=risk_df.loc[toSelect,:]
        ##Filtering on SNPs to test from the snp filter.
    
    if snp_feature_filter_df is not None:
        toSelect = set(np.unique(snp_feature_filter_df['snp_id'])).intersection(set(risk_df.index.values))
        risk_df=risk_df.loc[toSelect,:]
        ##Filtering on features to test from the combined feature snp filter.
    
    #Filtering for sites on non allosomes.
    if not skipAutosomeFiltering :
        annotation_df = annotation_df[annotation_df['chromosome'].map(lambda x: x in list(map(str, range(1, 23))))]
    
    feature_list = list(set(annotation_df.index)&set(phenotype_df.index))
    print("Number of features to be tested: " + str(len(feature_list)))
    print("Total number of variants to be considered, before variante QC and feature intersection: " + str(risk_df.shape[0]))
    
    if(phenotype_df.shape[1]<minimum_test_samples):
        print("Not enough samples with both genotype & phenotype data, for current number of covariates.")
        sys.exit()
    
    feature_variant_covariate_df = qtl_loader_utils.get_snp_feature_df(feature_variant_covariate_filename)
    
    return [phenotype_df, kinship_df, covariate_df, sample2individual_df, annotation_df, snp_filter_df, snp_feature_filter_df, genetically_unique_individuals, minimum_test_samples, feature_list, risk_df, chromosome, selectionStart, selectionEnd, feature_variant_covariate_df]

def run_PrsQtl_CisCorrected_analysis_load_intersect_phenotype_covariates_kinship_sample_mapping\
        (pheno_filename, anno_filename, geno_prefix, plinkGenotype, prsFile, minimum_test_samples= 10, relatedness_score=0.95, skipAutosomeFiltering = False, snps_filename=None,
         feature_filename=None, snp_feature_filename=None, selection='all', covariates_filename=None, kinship_filename=None, sample_mapping_filename=None, feature_variant_covariate_filename=None):
    
    selectionStart = None
    selectionEnd = None
    if(":" in selection):
        parts = selection.split(":")
        if("-" not in parts[1]):
            print("No correct sub selection.")
            print("Given in: "+selection)
            print("Expected format: (chr number):(start location)-(stop location)")
            sys.exit()
        chromosome = parts[0]
        if("-" in parts[1]):
            parts2 = parts[1].split("-") 
            selectionStart = int(parts2[0])
            selectionEnd = int(parts2[1])
    else :
        chromosome=selection

    ''' function to take input and intersect sample and genotype.'''
    #Load input data files & filter for relevant data
    #Load input data filesf

    phenotype_df = qtl_loader_utils.get_phenotype_df(pheno_filename)
    annotation_df = qtl_loader_utils.get_annotation_df(anno_filename)

    phenotype_df.columns = phenotype_df.columns.astype("str")
    phenotype_df.index = phenotype_df.index.astype("str")
    annotation_df.columns = annotation_df.columns.astype("str")
    annotation_df.index = annotation_df.index.astype("str")

    if(plinkGenotype):
        bim,fam,bed = qtl_loader_utils.get_genotype_data(geno_prefix)
        annotation_df.replace(['X', 'Y', 'XY', 'MT'], ['23', '24', '25', '26'],inplace=True)
        if chromosome=='X' :
            chromosome = '23'
        elif chromosome=='Y':
            chromosome = '24'
        elif chromosome=='XY':
            chromosome=='25'
        elif chromosome=='MT':
            chromosome=='26'
         #X  -> 23
         #Y  -> 24
         #XY -> 25
         #MT -> 26

    else :
        geno_prefix+='.bgen'
        print(geno_prefix)
    
    #Load risk data.
    risk_df = qtl_loader_utils.get_phenotype_df(prsFile)
    print("Intersecting data.")

    if(annotation_df.shape[0] != annotation_df.groupby(annotation_df.index).first().shape[0]): 
        print("Only one location per feature supported. If multiple locations are needed please look at: --extended_anno_file")
        sys.exit()

    ##Make sure that there is only one entry per feature id!.

    sample2individual_df = qtl_loader_utils.get_samplemapping_df(sample_mapping_filename,list(phenotype_df.columns),'sample')
    sample2individual_df.index = sample2individual_df.index.astype('str')
    sample2individual_df = sample2individual_df.astype('str')
    sample2individual_df['sample']=sample2individual_df.index
    sample2individual_df = sample2individual_df.drop_duplicates();
    ##Filter first the linking files!
    #Subset linking to relevant genotypes.
    orgSize = sample2individual_df.shape[0]
    ##subset based on plink
    sample2individual_df = sample2individual_df.loc[sample2individual_df['iid'].map(lambda x: x in list(map(str, fam.index))),:]
    ##subset based on risk scores
    sample2individual_df = sample2individual_df.loc[sample2individual_df['iid'].map(lambda x: x in list(map(str, risk_df.columns))),:]
    diff = orgSize- sample2individual_df.shape[0]
    orgSize = sample2individual_df.shape[0]
    print("Dropped: "+str(diff)+" samples because they are not present in the genotype and or risk file.")
    
    #Subset linking to relevant phenotypes.
    sample2individual_df = sample2individual_df.loc[np.intersect1d(sample2individual_df.index,phenotype_df.columns),:]
    diff = orgSize- sample2individual_df.shape[0]
    orgSize = sample2individual_df.shape[0]
    print("Dropped: "+str(diff)+" samples because they are not present in the phenotype file.")
    #Subset linking vs kinship.
    kinship_df = qtl_loader_utils.get_kinship_df(kinship_filename)
    if kinship_df is not None:
        #Filter from individual2sample_df & sample2individual_df since we don't want to filter from the genotypes.
        sample2individual_df = sample2individual_df[sample2individual_df['iid'].map(lambda x: x in list(map(str, kinship_df.index)))]
        diff = orgSize- sample2individual_df.shape[0]
        orgSize = sample2individual_df.shape[0]
        print("Dropped: "+str(diff)+" samples because they are not present in the kinship file.")
    #Subset linking vs covariates.
    covariate_df = qtl_loader_utils.get_covariate_df(covariates_filename)
    if covariate_df is not None:
        if np.nansum(covariate_df==1,0).max()<covariate_df.shape[0]: covariate_df.insert(0, 'ones',np.ones(covariate_df.shape[0]))
        sample2individual_df = sample2individual_df.loc[list(set(sample2individual_df.index) & set(covariate_df.index)),:]
        diff = orgSize- sample2individual_df.shape[0]
        orgSize = sample2individual_df.shape[0]
        print("Dropped: "+str(diff)+" samples because they are not present in the covariate file.")

    ###
    print("Number of samples with genotype & phenotype data: " + str(sample2individual_df.shape[0]))
    if(sample2individual_df.shape[0]<minimum_test_samples):
        print("Not enough samples with both genotype & phenotype data.")
        sys.exit()

    ##Filter now the actual data!
    #Filter phenotype data based on the linking files.
    phenotype_df = phenotype_df.loc[list(set(phenotype_df.index)&set(annotation_df.index)),sample2individual_df.index.values]

    #Filter kinship data based on the linking files.
    genetically_unique_individuals = None
    if kinship_df is not None:
        kinship_df = kinship_df.loc[np.intersect1d(kinship_df.index,sample2individual_df['iid']),np.intersect1d(kinship_df.index,sample2individual_df['iid'])]
        genetically_unique_individuals = get_unique_genetic_samples(kinship_df, relatedness_score);

    #Filter covariate data based on the linking files.
    
    snp_feature_filter_df= qtl_loader_utils.get_snp_feature_df(snp_feature_filename)
    try:
        feature_filter_df = qtl_loader_utils.get_snp_df(feature_filename)
    except:
        if feature_filename  is not None:
            feature_filter_df=pd.DataFrame(index=feature_filename)
    #Do filtering on features.
    if feature_filter_df is not None:
        phenotype_df = phenotype_df.loc[feature_filter_df.index,:]
        ##Filtering on features to test.
    if snp_feature_filter_df is not None:
        phenotype_df = phenotype_df.loc[np.unique(snp_feature_filter_df['feature']),:]
        ##Filtering on features  to test from the combined feature snp filter.

    #Prepare to filter on snps.
    snp_filter_df = qtl_loader_utils.get_snp_df(snps_filename)
    if snp_filter_df is not None:
        toSelect = set(snp_filter_df.index).intersection(set(risk_df.index.values))
        risk_df=risk_df.loc[toSelect,:]
        ##Filtering on SNPs to test from the snp filter.
    
    if snp_feature_filter_df is not None:
        toSelect = set(np.unique(snp_feature_filter_df['snp_id'])).intersection(set(risk_df.index.values))
        risk_df=risk_df.loc[toSelect,:]
        ##Filtering on features  to test from the combined feature snp filter.
    
    #Filtering for sites on non allosomes.
    if not skipAutosomeFiltering :
        annotation_df = annotation_df[annotation_df['chromosome'].map(lambda x: x in list(map(str, range(1, 23))))]
    
    #Determine features to be tested
    if chromosome=='all':
        feature_list = list(set(annotation_df.index)&set(phenotype_df.index))
    else:
        if not selectionStart is None :
            lowest = min([selectionStart,selectionEnd])
            highest = max([selectionStart,selectionEnd])
            annotation_df['mean'] = ((annotation_df["start"] + annotation_df["end"])/2)
            feature_list = list(set(annotation_df.iloc[(annotation_df['chromosome'].values==chromosome) & (annotation_df['mean'].values>=lowest) & (annotation_df["mean"].values<highest)].index.values)&set(phenotype_df.index))
            del annotation_df['mean']
        else :
            feature_list = list(set(annotation_df[annotation_df['chromosome']==chromosome].index)&set(phenotype_df.index))

    print("Number of features to be tested: " + str(len(feature_list)))
    print("Total number of variants to be considered, before variante QC and feature intersection: " + str(bim.shape[0]))
    
    if(phenotype_df.shape[1]<minimum_test_samples):
        print("Not enough samples with both genotype & phenotype data, for current number of covariates.")
        sys.exit()

    feature_variant_covariate_df = qtl_loader_utils.get_snp_feature_df(feature_variant_covariate_filename) 
    
    return [phenotype_df, kinship_df, covariate_df, sample2individual_df, annotation_df, snp_filter_df, snp_feature_filter_df, genetically_unique_individuals, minimum_test_samples, feature_list,bim,fam,bed, risk_df, chromosome, selectionStart, selectionEnd, feature_variant_covariate_df]

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
    return (seq[pos:pos + np.int(size)] for pos in range(0, len(seq), np.int(size)))

def get_unique_genetic_samples(kinship_df, relatedness_score):
#    tril returns the lower triungular.
#    if two lines are > identity then  kinship_df>=relatedness_score should have an offdiagonal 1.
#    if there is one 1  in the tril then it means that  in the upper triul there was a line withe identical genotype
    return (kinship_df.index[(np.tril(kinship_df>=relatedness_score,k=-1)).sum(1)==0])

def force_normal_distribution(phenotype, method='gaussnorm', reference=None):
    _doc='rank transform x into ref/ gaussian;keep the range; keep ties'
    
    if method=='log':
        return np.log(1+phenotype)
    
    if method=='log_standardize':
        temp=np.log(1+phenotype)
        return (temp-np.nanmean(temp))/np.nanstd(temp)

    if method=='arcsin':
        return np.arcsin(np.sqrt(phenotype))
    
    if method=='arcsin_standardize':
        temp=np.arcsin(np.sqrt(phenotype))
        return (temp-np.nanmean(temp))/np.nanstd(temp)
    
    if method=='standardize':
        return (phenotype-np.nanmean(phenotype))/np.nanstd(phenotype)
                        
    indextoupdate = np.isfinite(phenotype)
    y1 = phenotype[indextoupdate]
    yuni,yindex=np.unique(y1, return_inverse=True)
    phenotypenorm=phenotype.copy()
    
    if method =='gaussnorm':

        sref = scst.norm.isf(np.linspace(0.001, 0.999,num=yuni.shape[0])[::-1])
        phenotypenorm[indextoupdate]=sref[yindex]
        return phenotypenorm
    
    elif method=='ranknorm':
        try:
            xref1=np.unique(reference[np.isfinite(reference)])
            sref=np.sort(xref1)[np.linspace(0,xref1.shape[0]-0.001, num=y1.shape[0]).astype(int)]
        except:
            print ('reference missing. provide reference to force_normal_distribution or choose gaussnorm')
            return 1
        phenotypenorm[indextoupdate]=sref[np.argsort(np.argsort(y1))]
        return phenotypenorm
    
    elif method=='ranknorm_duplicates':
        try:
            xref1=np.unique(reference[np.isfinite(reference)])### unique values from reference
            sref=np.sort(xref1)[np.linspace(0,xref1.shape[0]-0.001, num=yuni.shape[0]).astype(int)]
        except: 
            print ('reference missing. provide reference to force_normal_distribution or choose gaussnorm')
            return 1
        phenotypenorm[indextoupdate]=sref[yindex]
        return phenotypenorm  
    
    else:
        print ('methods are: log, log_standardize, standardize, gaussnorm, ranknorm, ranknorm_duplicates, arcsin, arcsin_standardize')

#get_shuffeld_genotypes_preserving_kinship(genetically_unique_individuals, relatedness_score, snp_matrix_DF,kinship_df.loc[individual_ids,individual_ids], n_perm)
def get_shuffeld_genotypes_preserving_kinship(genetically_unique_individuals, relatedness_score, snp_matrix_DF,kinship_df1,n_perm):

#    snp_matrix_DF.iloc[np.unique(snp_matrix_DF.index,return_index=1)[1]].shape
    '''take only one line for replicates (those with the same name)'''
    temp=snp_matrix_DF.iloc[np.unique(snp_matrix_DF.index,return_index=1)[1]].copy(deep=True)
    u_snp_matrix = temp.loc[genetically_unique_individuals,:]
    kinship_df1=kinship_df1.iloc[np.unique(kinship_df1.index,return_index=1)[1],np.unique(kinship_df1.index,return_index=1)[1]]
    '''has replicates but not same lines form donor (np.setdiff1d(individual_ids,genetically_unique_individuals))'''
    #Shuffle and reinflate
    locationBuffer = np.zeros(snp_matrix_DF.shape[0], dtype=np.int)
    #Prepare location search for permuted snp_matrix_df.
    index = 0
    index_samples = np.arange(u_snp_matrix.shape[0])
    for current_name in genetically_unique_individuals :
        selection = kinship_df1.loc[current_name].values>=relatedness_score
        locationBuffer[np.where(selection)] = index
        index +=1
    snp_matrix_copy = np.zeros((snp_matrix_DF.shape[0],snp_matrix_DF.shape[1]*n_perm))
    counter = 0
    end = (snp_matrix_DF.shape[1])
    for perm_id in range(0,n_perm) :
        np.random.shuffle(index_samples)
        temp_u = u_snp_matrix.values[index_samples,:]
        snp_matrix_copy[:,counter:end] = temp_u[locationBuffer,:]
        counter+= snp_matrix_DF.shape[1]
        end+= snp_matrix_DF.shape[1]
    return(snp_matrix_copy)

def get_shuffeld_genotypes(snp_matrix_DF,n_perm):
    snp_matrix_copy = np.zeros((snp_matrix_DF.shape[0],snp_matrix_DF.shape[1]*n_perm))
    counter = 0
    end = (snp_matrix_DF.shape[1])

    index_samples = np.arange(snp_matrix_DF.shape[0])
    for perm_id in range(0,n_perm) :
        np.random.shuffle(index_samples)
        snp_matrix_copy[:,counter:end] = snp_matrix_DF.values[index_samples,:]
        counter+= snp_matrix_DF.shape[1]
        end+= snp_matrix_DF.shape[1]
    return(snp_matrix_copy)

def qtl_plot(snp_matrix_DF, phenotype,K=None, M=None,LMM=None,snp=None,show_reg_cov=True):
    if LMM is None:
        LMM = limix.qtl.qtl_test_lmm(snp_matrix_DF.values, phenotype,K=K,M=M,verbose=False)
        
    if snp is None:
        snp=snp_matrix_DF.values[:,np.argmin(LMM.variant_effsizes)]           
    
    cov=LMM.null_covariate_effsizes;betacov=np.array([cov[key] for key in cov.keys()])


    temp=pd.DataFrame(data=np.hstack([np.vstack([snp,phenotype,np.zeros(phenotype.shape[0]).astype(bool)]),\
                                          np.vstack([snp,phenotype-np.dot(M,betacov[:,None])[:,0],np.ones(phenotype.shape[0]).astype(bool)])]).T,\
                                            columns=['Variant','Phenotype donor','Covariates regressed'])
    
    temp['Covariates regressed']=temp['Covariates regressed'].astype(bool)
 
    indexunique=np.unique(np.array([l.split('-')[1] for l in snp_matrix_DF.index]),return_index=1)[1]
    ax = sb.boxplot(y=temp['Phenotype donor'],x=temp['Variant'])
    for patch in ax.artists:
        r, g, b, a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, .03))
 
#    sb.swarmplot(y=temp['Phenotype donor'],x=temp['Variant'],color='grey')
    if show_reg_cov:
        index1=np.hstack([indexunique,indexunique+phenotype.shape[0]])
        sb.swarmplot(y=temp['Phenotype donor'].iloc[index1],x=temp['Variant'].iloc[index1],  hue=temp['Covariates regressed'].iloc[index1], split=True)
    else:
        index1= indexunique 
        sb.swarmplot(y=temp['Phenotype donor'].iloc[index1],x=temp['Variant'].iloc[index1])

def do_snp_selection(feature_id, annotation_df, bim, cis_mode, window_size, skipAutosomeFiltering = False):
    annotation_sub_df = annotation_df.loc[[feature_id],:]
    list_of_snp_dfs = []
    snpQuery = None
    for feature_id,annotation_ds in annotation_sub_df.iterrows():
        chrom = str(annotation_ds.loc['chromosome'])
        start = annotation_ds.loc['start']
        end = annotation_ds.loc['end']
        # make robust to features selfpecified back-to-front
        lowest = min([start,end])
        highest = max([start,end])
        if (cis_mode) :
            # for cis, we sequentially add snps that fall within each region
            snpQuery = bim.query("chrom == '%s' & pos > %d & pos < %d" % (chrom, lowest-window_size, highest+window_size))
            list_of_snp_dfs.append(snpQuery)
        else :
            # for trans, we sequentially exclude snps that fall within each region
            if snpQuery is None :
                #Initial search for trans window.
                snpQuery = bim.query("(chrom == '%s' & (pos < %d | pos > %d))|chrom != '%s'" % (chrom, lowest-window_size, highest+window_size,chrom))
            else :
                #Here we start intersecting regions we want to exclude.
                snpQuery = snpQuery.query("(chrom == '%s' & (pos < %d | pos > %d))|chrom != '%s'" % (chrom, lowest-window_size, highest+window_size,chrom))
    
    if (cis_mode):
        selected_snp_df = pd.concat(list_of_snp_dfs).drop_duplicates()
    else:
        selected_snp_df = snpQuery
    if not skipAutosomeFiltering :
        # filtering for sites on non allosomes.
        selected_snp_df = selected_snp_df.loc[selected_snp_df['chrom'].map(lambda x: x in list(map(str, range(1, 23)))).values]
    
    return selected_snp_df

def reduce_snp(snp_df, threshold=0.975):
    ''' input a snp df  samples(rows) x snps( columns)'''
    ''' returns a df with columns: 'lead_snp_id' the pruned snp names and snp_id'''
    allsnps=snp_df.columns
    cor=abs(np.corrcoef(snp_df.T))>=threshold
    cor1=pd.DataFrame(data=np.triu(cor,k=1),index=allsnps,columns=allsnps).sum(1)
    cor2=pd.DataFrame(data=np.triu(cor,k=1),index=allsnps,columns=allsnps)
    uniquesnps=cor1[cor1==0].index.values
    duplicatedsnps=cor1[cor1>0].index.values
    '''the unique ones'''
    rez=pd.DataFrame(data=uniquesnps,index=uniquesnps,columns=['snp_id'])
    ''' the ones that are in LD.. if a,c and a,b but not b,c returns a,b)'''
    rez2=pd.DataFrame(data=duplicatedsnps,index=allsnps[np.argmax(cor2.loc[duplicatedsnps].values*cor2.loc[duplicatedsnps].values.sum(0),1)],columns=['snp_id'])
    rez=pd.concat([rez,rez2])
    rez['lead_snp_id']=rez.index
    return(rez)

def regressOut(Y, X):
    """
    regresses out X from Y
    """
    Xd = la.pinv(X)
    Y_out = Y-X.dot(Xd.dot(Y))
    return Y_out
