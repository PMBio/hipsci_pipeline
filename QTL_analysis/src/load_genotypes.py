import limix.io
import numpy as np
import math
import scipy

(bim, fam, bed) = limix.io.read_plink('../data/geuvadis_test_data/Geuvadis_chr1')

#make missing values 0 for now (i.e. homozygous for first allele)
#needs to be replaced with a snp-specific value
missing_value = 0

def bed_to_genotype_value(input_value):
	value_mapping = [0,1,2,missing_value]
	return value_mapping[input_value]

#rows are individuals, SNPs are columns
chrsnps = np.array(np.vectorize(bed_to_genotype_value)(bed.compute()),dtype=float).transpose()
print(chrsnps[0:4,0:8])

kchr = chrsnps
#average by column and subtract
kchr -= chrsnps.mean(axis=0)
print(kchr)
kchr /= kchr.std(axis=0)
print(kchr)


nS = chrsnps.shape[0]
nN = chrsnps.shape[1]

kinship = scipy.dot(kchr, kchr.T)

print(kinship.shape)
print(kinship)
