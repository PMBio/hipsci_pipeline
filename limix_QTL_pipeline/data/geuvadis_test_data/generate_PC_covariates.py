from sklearn.decomposition import PCA
import pandas as pd
import numpy as np

data = pd.read_csv('Geuvadis_CEU_YRI_Expr.txt',sep='\t',index_col=0)

nPCA_components = 5

pca = PCA(n_components=nPCA_components)
pca.fit(data.as_matrix())

covariate_df = pd.DataFrame(data=pca.components_.transpose(),index=data.columns,columns=np.arange(nPCA_components)+1)

covariate_df.to_csv('Geuvadis_CEU_YRI_covariates.txt',sep='\t')