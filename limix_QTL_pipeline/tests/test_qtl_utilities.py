from context import qtl_utilities
import pandas as pd
import numpy as np
import importlib
importlib.reload(qtl_utilities)

##### set up test

np.random.seed(0)

relatedness_score = 0.95
n_perm = 1

n_snps = 9
n_individuals = 5
snps = ['snp{}'.format(x) for x in range(n_snps)]
individuals = ['donor{}'.format(x) for x in range(n_individuals)]
samples = individuals + individuals[:n_individuals//2]
n_samples = len(samples)

matrix = np.random.binomial(2,0.4,(n_individuals, n_snps))
snp_matrix_DF = pd.DataFrame(columns=snps, index=individuals, data=matrix)

# extend with replicated individuals
snp_matrix_DF_2 = snp_matrix_DF.iloc[:n_individuals//2,:]
snp_matrix_DF_2.index = [x+'_2' for x in snp_matrix_DF_2.index]
snp_matrix_DF = pd.concat([snp_matrix_DF, snp_matrix_DF_2])

# sample mapping to individuals, to allow building the kinship
sm_df = pd.DataFrame(index=snp_matrix_DF.index, columns=['individual'])
sm_df['individual'] = sm_df.index
sm_df['individual'] = sm_df['individual'].apply(lambda x: x.split('_')[0])


#### initialise block diagonal kinship matrix
dummy_df = pd.get_dummies(sm_df['individual'])
kinship_matrix = np.dot(dummy_df.values, dummy_df.values.T)
kinship_df1 = pd.DataFrame(index=dummy_df.index, columns=dummy_df.index, data=kinship_matrix)

# # use an unsual kinship
# kinship_matrix = np.diag([0.001 for _ in snp_matrix_DF.index])
# kinship_df1 = pd.DataFrame(index=snp_matrix_DF.index, columns=snp_matrix_DF.index, data=kinship_matrix)


genetically_unique_individuals = qtl_utilities.get_unique_genetic_samples(kinship_df1, relatedness_score)

snp_matrix_copy = qtl_utilities.get_shuffeld_genotypes_preserving_kinship(genetically_unique_individuals, relatedness_score, snp_matrix_DF,kinship_df1,n_perm)

# put snp_matrix_copy in a dataframe:
perm_df = pd.DataFrame(data=snp_matrix_copy, index=snp_matrix_DF.index, columns=snp_matrix_DF.columns)

print('Original dataframe:')
print(snp_matrix_DF)

print('')
print('Permuted results:')
print(perm_df)
print('')
print('Equality of permuted genotypes for repeated samples from donor0 and donor1:')
print((perm_df.loc['donor0',:]==perm_df.loc['donor0_2']).all())
print((perm_df.loc['donor1',:]==perm_df.loc['donor1_2']).all())
