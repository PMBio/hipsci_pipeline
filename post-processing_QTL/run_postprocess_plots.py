import sys
import numpy as np
import pandas
import argparse
from matplotlib.patches import Rectangle
import matplotlib.gridspec as gridspec
from matplotlib_venn import venn2, venn2_circles
from matplotlib_venn import venn3, venn3_circles
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import matplotlib as mpl
#sys.path.append('../../hipsci_pipeline/post-processing_QTL/')

from scripts.postprocess_functions_genome_results import *
from scripts.postprocess_functions_plots import *


def get_args():
    parser = argparse.ArgumentParser(description='Run QTL analysis given genotype, phenotype, and annotation.')
    parser.add_argument('-folder_data','--folder_data',required=True)
    parser.add_argument('-folder_destination','--folder_destination',required=False,default=None)
    parser.add_argument('-traits','--traits',required=True)
    parser.add_argument('-trait_labels','--trait_labels',required=False,default=None)
    parser.add_argument('-chromosomes','--chromosomes',required=False,default=np.arange(1,23).astype('U'))
    args = parser.parse_args()
    
    return args


if 1==1:
    if __name__=='__main__':
        print('reading arguments')
        args = get_args()
        print(args)
        chromosomes = args.chromosomes.split(',')
        folder_data =args.folder_data
        print('step1')
        folder_destination = args.folder_destination if args.folder_destination is not None else folder_data
        traits =args.traits.split(',')
        print('step2')
        trait_labels =args.trait_labels.split(',') if args.trait_labels is not None else traits
else:
    folder_data ='/Users/mirauta/Results/hipsci/QTL/'
    folder_destination='/Users/mirauta/Results/hipsci/'
    traits=['New2_protein', 'New2_protein_perm' ];
    trait_labels=traits
    chromosomes=np.array(['21'])

for trait in traits:
    summary_gene_feature(output_file='qtl_results_genome', feature_report='ensembl_gene_id', chr_list=chromosomes,p_value_field='corr_p_value', folder_data=folder_data+trait,local_adjustment_method=None)

plot_summary(folder_name=folder_data,\
             folder_destination=folder_destination+'Images_pipeline/',plot_name='qtl_summary',\
                 traits=traits, trait_labels=trait_labels,
                 file_name_results_genome='ensembl_gene_id_qtl_results_genome',   qtl_results_file='qtl_results_', \
                 colors=np.array(['orange','darkblue','green','red']),cis=2.5*10**5,figsize=(8,8),gene_list=None,\
                 plot_calibration_flag=True)

rez_pro_pep=plot_replication(rez=None,folder_data =folder_data,\
    traits=traits[::-1],trait_labels= trait_labels, qtl_results_file='qtl_results_',    snp_metadata_file='snp_metadata_',    feature_metadata_file='feature_metadata_',\
    results_genome_file='qtl_results_genome',    feature_report='ensembl_gene_id',\
    folder_destination=folder_destination+'Images_pipeline/', figsize=7,red_dots_features=None, \
                             plot_name='protein_peptide_replication') 
#names=['p_value', 'replicated_p_value', 'replicated_self_p_value', 'feature_id', 'gene_name', 'snp_id', 'chromosome', 'strand', 'position', 'ensembl_gene_id']
df= pandas.DataFrame(data=np.array([rez_pro_pep[key] for key in names ]).T, index=rez_pro_pep['ensembl_gene_id'],columns=names)
#df.to_csv(path_or_buf=folder_data+traits[0]+'_'+traits[1]+'_qtl_results.txt',mode='w', sep='\t', columns=None, header=True, index=True)
##
genes1=df.index[(df['p_value'].values.astype(float)<10**-4)&(df['replicated_p_value'].values.astype(float)<10**-4)]
#print (genes1.shape)

for i, g in enumerate(genes1 ):
    plot_manhatan_alone( gene_ensembl_id= genes1[i],folder_name=folder_data,\
                        folder_destination=folder_destination+'/manhattan/',\
                    plot_name='manhattan_pepQTLonly',traits=traits,trait_labels=trait_labels,file_name_results_genome='ensembl_gene_id_qtl_results_genome',\
                    qtl_results_file='qtl_results_',colors=np.array(['black','green','orange','blue']), figsize=4)
    plt.show()
