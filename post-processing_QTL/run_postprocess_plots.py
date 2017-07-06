import sys
from matplotlib.patches import Rectangle
import matplotlib.gridspec as gridspec
from matplotlib_venn import venn2, venn2_circles
from matplotlib_venn import venn3, venn3_circles
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import matplotlib as mpl
sys.path.append('../../hipsci_pipeline/post-processing_QTL/')
import numpy as np
import pandas

sys.path.append('../')
from scripts.postprocess_functions_genome_results import *
from scripts.postprocess_functions_plots import *


folder_data ='/Users/mirauta/Results/hipsci/QTL/'
folder_destination='/Users/mirauta/Results/hipsci/'
traits=['New2_protein', 'New2_protein_perm' ];trait_labels=traits

trait=traits[0];summary_gene_feature(output_file='qtl_results_genome', feature_report='ensembl_gene_id',\
                         chr_list=np.array([1,20,21]),p_value_field='p_value', folder_data=folder_data+trait,local_adjustment_method='Bonferroni')
trait=traits[1];summary_gene_feature(output_file='qtl_results_genome', feature_report='ensembl_gene_id',\
                         chr_list=np.array([1,20,21]),p_value_field='corr_p_value', folder_data=folder_data+trait,local_adjustment_method=None)
for trait in traits:
    summary_gene_feature(output_file='qtl_results_genome', feature_report='ensembl_gene_id',\
                         chr_list=np.array([1,20,21]),p_value_field='corr_p_value', folder_data=folder_data+trait,local_adjustment_method=None)

traits=['New2_protein', 'New2_protein' ];trait_labels=traits
plot_summary(folder_name=folder_data,\
             folder_destination=folder_destination+'Images_pipeline/',plot_name='qtl_summary2',\
                 traits=traits, trait_labels=trait_labels,
                 file_name_results_genome='ensembl_gene_id_qtl_results_genome',   qtl_results_file='qtl_results_', \
                 colors=np.array(['orange','darkblue','green','red']),cis=2.5*10**5,figsize=(8,8),gene_list=None,\
                 plot_calibration_flag=True)
#plt.show()


#traits=['New2_protein','New2_protein']
#trait_labels=traits
#
rez_pro_pep=plot_replication(rez=None,folder_data =folder_data,\
    traits=traits[::-1],trait_labels= trait_labels, qtl_results_file='qtl_results_',    snp_metadata_file='snp_metadata_',    feature_metadata_file='feature_metadata_',\
    results_genome_file='qtl_results_genome',    feature_report='ensembl_gene_id',\
    folder_destination=folder_destination+'Images_pipeline/', figsize=7,red_dots_features=None, \
                             plot_name='protein_peptide_replication') 
plt.show()

#
#names=['p_value', 'replicated_p_value', 'replicated_self_p_value', 'feature_id', 'gene_name', 'snp_id', 'chromosome', 'strand', 'position', 'ensembl_gene_id']
#df= pandas.DataFrame(data=np.array([rez_pro_pep[key] for key in names ]).T, index=rez_pro_pep['ensembl_gene_id'],columns=names)
#df.to_csv(path_or_buf=folder_data+traits[0]+'_'+traits[1]+'_qtl_results.txt',mode='w', sep='\t', columns=None, header=True, index=True)
##
#genes1=df.index[(df['p_value'].values.astype(float)<10**-4)&(df['replicated_p_value'].values.astype(float)<10**-4)]
#print (genes1.shape)
genes1=['ENSG00000101160']
for i, g in enumerate(genes1 ):
    plot_manhatan_alone( gene_ensembl_id= genes1[i],folder_name=folder_data,\
                        folder_destination='/Users/mirauta/Results/hipsci/Images_pipeline/manhattan/',\
                    plot_name='manhattan_pepQTLonly',traits=traits,trait_labels=trait_labels,file_name_results_genome='ensembl_gene_id_qtl_results_genome',\
                    qtl_results_file='qtl_results_',colors=np.array(['black','green','orange','blue']), figsize=4)
    plt.show()
