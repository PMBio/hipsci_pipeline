#import pandas as pd
from load_genotypes import load_genotypes_plink
from generate_kinship import generate_kinship

geno_prefix = '../data/geuvadis_test_data/Geuvadis_chr1'
pheno_filename = '../data/geuvadis_test_data/Geuvadis_CEU_YRI_Expr.txt'
anno_filename = '../data/geuvadis_test_data/Geuvadis_CEU_YRI_Annot.txt'

#bim,genotype_mat = load_genotypes_plink(geno_prefix)
#kinship = generate_kinship(genotype_mat)

phenotypes = pd.read_csv(pheno_filename,sep='\t',index_col=0)
annotation = pd.read_csv(anno_filename,sep='\t',index_col=1)