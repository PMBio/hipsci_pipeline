#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from run_QTL_analysis import run_QTL_analysis,merge_QTL_results
from hashlib import md5

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

    results_md5_dict = {output_dir+'qtl_results_1.txt':'4e269da9c707b7b5a73cfc6effcdf45f'}
    for f in results_md5_dict.keys():
        assert(md5(f).hexdigest()==results_md5_dict[f])


    data_path = '../data/geuvadis_CEU_test_data/'
    geno_prefix = data_path+'Genotypes/Geuvadis'
    pheno_filename = data_path+'Expression/Geuvadis_CEU_Expr.txt'
    anno_filename = data_path+'Expression/Geuvadis_CEU_formatted_annotation_data.txt'
    
    output_dir = data_path+'TestOutput/limix_QTL_results/'
        
    ws = 250000
    
    for chromosome in ['1','2']:
        run_QTL_analysis(pheno_filename,anno_filename,geno_prefix,chromosome,ws,output_dir)
    merge_QTL_results(output_dir)
    
    results_md5_dict = {output_dir+'qtl_results_1.txt':'a7429ca50ac531c79f3626616f079f91',
                        output_dir+'qtl_results_2.txt':'f2423bbe8605a98d614372d62e3671c3'}
    for f in results_md5_dict.keys():
        assert(md5(f).hexdigest()==results_md5_dict[f])