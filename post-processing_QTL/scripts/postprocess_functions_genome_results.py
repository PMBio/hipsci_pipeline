
import sys

#import statsmodels.api as sm
import h5py
import numpy as np
import limix.stats.fdr as FDR
from  Utilities import *
import pandas

def local_adjustment(pv, N=1,  method=''):
    if method=='': return pv
    if method=='Bonferroni': return pv*N
  
    
#def summary_gene_feature(
qtl_results_file='qtl_results_';snp_metadata_file='snp_metadata_'; feature_metadata_file='feature_metadata_';output_file='qtl_results_genome';\
feature_report='gene_ensembl_id';chr_list=[9];\
folder_data = '/Users/mirauta/Git/hipsci_pipeline/limix_QTL_pipeline/data/geuvadis_CEU_YRI_test_data/\
protein_test_perm2'
folder_data='/Users/mirauta/Data/MS/hipsci/TMT/protein_perm';\
rewrite=False; local_adjustment_method='Bonferroni'
#):

#    _doc=" returns two files: 1) per ferature and 2) per gene\
#    min per gene; min adjusted per gene; \
#    min per feature; min adjusted per ferature "

for ichr,chr in enumerate(chr_list):
    print ('chromosome: '+str(chr))

    frez=h5py.File(folder_data+'/'+qtl_results_file+str(chr)+'.h5','r')
    frezkeys= np.array([k.replace('_i_','') for k in list(frez.keys())])
    ffea= pandas.read_table(folder_data+'/'+feature_metadata_file+ str(chr)+'.txt', sep='\t')
    fsnp= pandas.read_table(folder_data+'/'+snp_metadata_file+ str(chr)+'.txt', sep='\t',index_col=0)
    
    indexnan=np.where((ffea[feature_report])!=(ffea[feature_report]))[0]
    for i in indexnan:
        ffea[feature_report][i]='gene_'+str(ffea['chromosome'][i])+'_'+str(ffea['start'][i])
    ffea_feature=ffea.set_index('feature_id', drop=False).transpose()
    ffea_report_feature=ffea.set_index(feature_report, drop=False).transpose()
    
                         
    if ichr==0:
        fOut=h5py.File(folder_data+'/'+feature_report+'_'+output_file+'.h5','w')
    else:
        fOut=h5py.File(folder_data+'/'+feature_report+'_'+output_file+'.h5','r+')
    #
    # for each report_feature  create h5 groups
    count=0
    for ifea,report_feature in enumerate(np.unique(list(ffea_report_feature))):
        print (report_feature)
#        print (ffea_report_feature[report_feature]['feature_strand'])
#        print (ffea_report_feature[report_feature]['chromosome'])
#        
        #select features for which qtl was computed
        features=np.intersect1d(np.array( ffea_report_feature[report_feature].transpose()['feature_id'] ), frezkeys)
        if len(features) >=1:
            print ('results for '+ report_feature)
    
    #        try: 
            fg=fOut.create_group(report_feature)
            
            pv= np.array([frez[f]['p_value'] for f in  features ])
            print (pv.shape)
            for i, ppv in  enumerate(pv): ppv[ppv!=ppv]=1
            
    #        fgd=fg.create_group('data')
    #        for key in  list(frez[features[0]].dtype.fields):
    #            fgd.create_dataset(key,data=np.array([frez[f][key] for f in  features ]))
    #        
            fgm=fg.create_group('metadata')
            for key in ffea_feature[features[0]].keys():
     
                if isinstance(ffea_feature[features[0]][key],int) :
                    fgm.create_dataset(key,data=np.array([ffea_feature[f][key] for f in  features ]))
                else:
                    fgm.create_dataset(key,data=np.array([ffea_feature[f][key] for f in  features ]).astype('S'))
            
            fgd=fg.create_group('data')
            fgd.create_dataset('p_value',data=  np.hstack(pv) [None]) 
                             
            fgs=fg.create_group('summary_data')
            fgs.create_dataset('min_p_value',data= np.nanmin(np.hstack(pv)) [None])
    
            fgs.create_dataset('min_p_value_local_adjusted',data=np.nanmin(local_adjustment(np.hstack(pv),method=local_adjustment_method))[None])
            minfeatureindex=np.argmin([np.nanmin(ppv) *len(ppv) for ppv in pv])
            
            fgs.create_dataset('min_p_value_feature_id',data= np.array(fgm ['feature_id'][minfeatureindex].astype('S'))[None])
    
            fgs.create_dataset('min_p_value_snp_id',data= np.array(frez[features[minfeatureindex]]['snp_id'][np.nanargmin(pv[minfeatureindex])].astype('S'))[None])
    
            count+=1
    #        except: print (report_feature+' is duplicated ')
            
            
    frez.close()
    ffea 
    fsnp 
    fOut.close()
print (count)  
fOut=h5py.File(folder_data+'/'+feature_report+'_'+output_file+'.h5','r')
list(fOut.keys())
ppv=np.array([fOut[g]['summary_data']['min_p_value_local_adjusted'][:] for g in list(fOut.keys())])
ppv
fOut.close()
#fOut['CRYZ/summary_data']['min_p_value_local_adjusted'][:]
#list(h5py.File(folder_name+'/'+feature_name+'/'+feature_report+'_results_genome.h5','r').keys())
#
#frez=h5py.File(folder_name+'/'+feature_name+'/'+qtl_results_file+str(chr)+'.h5','r')['Q9Y2G5']
#frez[:]

#    
#fIn=h5py.File(folder_name+'/'+feature_name+'/'+file_name_feature+'.h5','r')
#fOut=h5py.File(folder_name+'/'+feature_name+'/'+file_name_gene+'.h5','w')
#genes=np.array([fIn[feature]['metadata/gene_id'][0] for feature in fIn.keys()])
#features=np.array([feature for feature in fIn.keys()])
#
#for gene in np.unique(genes):
#    
#    gene_features=features[genes==gene]
#    fg=fOut.create_group(gene)
#    
#    fgm=fg.create_group('metadata')
#    for key in fIn[gene_features[0]]['metadata'].keys():
#        fgm.create_dataset(key,data=fIn[gene_features[0]]['metadata'][key])
#        
#    fgs=fg.create_group('summary_data')
#    for key in ['min_pv','min_pv_local_adjusted']:
#        fgs.create_dataset(key,data=np.nanmin(np.array([fIn[feature+'/summary_data/'+key][:] for feature in gene_features]))[None])    
#    fgs.create_dataset('min_pv_feature_id',data=np.array(gene_features[np.nanargmin(np.array([fIn[feature+'/summary_data/min_pv_local_adjusted'][:] for feature in gene_features]))]).astype('S')[None])
#
#fIn.close()
#fOut.close()
#        
#    return [folder_name+feature_name+file_name_feature+'.h5',folder_name+feature_name+file_name_gene+'.h5']
# 
#summary_gene_feature(chr_list=[21], folder_name='/Users/mirauta/Data/MS/hipsci/TMT',\
#                     feature_name='protein_test',local_adjustment_method='')
 