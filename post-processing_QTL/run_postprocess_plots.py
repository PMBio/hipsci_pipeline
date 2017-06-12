from matplotlib.patches import Rectangle
import matplotlib.gridspec as gridspec
#import sys
from matplotlib_venn import venn2, venn2_circles
from matplotlib_venn import venn3, venn3_circles
import matplotlib.pyplot as plt

from matplotlib import colors as mcolors
import matplotlib as mpl
#import h5py
import numpy as np
import pandas

#import limix.stats.fdr as FDR
#import os
#import scipy.stats as scst
#import statsmodels.sandbox.stats.multicomp as scst2
from postprocess_functions_plots import *
from postprocess_functions_genome_results import *

folder_data ='/Users/mirauta/Data/MS/hipsci/TMT/'
folder_destination='/Users/mirauta/Results/hipsci/'
traits=['peptide_lines_filtered_unique_genes_filtered_intensity_scalarcorrected','protein_lines_filtered_unique_genes_filtered_intensity_scalarcorrected']
trait_labels=['Peptide QTL','Protein QTL (method 2)']

for trait in traits:
    summary_gene_feature(output_file='qtl_results_genome', feature_report='ensembl_gene_id',chr_list=np.arange(10,12),folder_data='/Users/mirauta/Data/MS/hipsci/TMT/'+trait)

plot_summary(folder_name='/Users/mirauta/Data/MS/hipsci/TMT/',\
             folder_destination=folder_destination+'Images_pipeline/',plot_name='qtl_summary',\
                 traits=traits, trait_labels=trait_labels,
                 file_name_results_genome='ensembl_gene_id_qtl_results_genome',   qtl_results_file='qtl_results_', \
                 colors=np.array(['orange','darkblue','green']),cis=1.5*10**5,figsize=(10,10),gene_list=None,plot_calibration_flag=True)

rez_pep_pro=plot_replication(rez=None,folder_data ='/Users/mirauta/Data/MS/hipsci/TMT/',    traits=traits,trait_labels= trait_labels,\
    qtl_results_file='qtl_results_',    snp_metadata_file='snp_metadata_',    feature_metadata_file='feature_metadata_',\
    results_genome_file='qtl_results_genome',    feature_report='ensembl_gene_id',\
    folder_destination='/Users/mirauta/Results/hipsci/Images_pipeline/', figsize=7,\
    red_dots_features=None, plot_name='peptide_protein_scaled') 


names=['p_value', 'replicated_p_value', 'replicated_self_p_value', 'feature_id', 'gene_name', 'snp_id', 'chromosome', 'strand', 'position', 'ensembl_gene_id']
df= pandas.DataFrame(data=np.array([rez_pep_pro[key] for key in names ]).T, index=rez_pep_pro['ensembl_gene_id'],columns=names)
df.to_csv(path_or_buf=folder_data+traits[0]+'_'+traits[1]+'_qtl_results.txt',mode='w', sep='\t', columns=None, header=True, index=True)

genes1=df.index[(df['p_value'].values.astype(float)<10**-8)&(df['replicated_p_value'].values.astype(float)<10**-4)]
print (genes1.shape)
for i, g in enumerate(genes1 ):
    plot_manhatan_alone( gene_ensembl_id= genes1[i],folder_name='/Users/mirauta/Data/MS/hipsci/TMT/',\
                        folder_destination='/Users/mirauta/Results/hipsci/Images_pipeline/manhattan/',\
                    plot_name='manhattan_pepQTLonly',traits=traits,trait_labels=trait_labels,file_name_results_genome='ensembl_gene_id_qtl_results_genome',\
                    qtl_results_file='qtl_results_',colors=np.array(['black','green','orange','blue']), figsize=4)
    plt.show()
