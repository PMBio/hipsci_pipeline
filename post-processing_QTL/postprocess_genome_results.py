
import sys

#import statsmodels.api as sm
import h5py
import numpy as np
import limix.stats.fdr as FDR


def local_adjustment(pv, N=1,  method=''):
    if method=='': return pv
    if method=='Bonferroni': return pv*N
  
    
def summary_gene_feature(qtl_results_file='QTL_results_', snp_metadata_file='SNP_metadata_',\
                         feature_metadata_file='Feature_metadata_',\
                         chr_list=np.arange(1,23),\
                         folder_name='/Users/mirauta/Scripts/QTL_pipeline_plots/',\
                         feature_name='Exon_1234',file_name_feature='Feature_results_genome',\
                         file_name_gene='Feature_summarized_results_genome',\
                         rewrite=False, local_adjustment_method='Bonferroni'):
    _doc=" returns two files: 1) per ferature and 2) per gene\
    min per gene; min adjusted per gene; \
    min per feature; min adjusted per ferature "
    
    for ichr,chr in enumerate(chr_list):
        print ('chromosome: '+str(chr))
        
        frez=h5py.File(folder_name+'/'+feature_name+'/'+qtl_results_file+'chr'+str(chr)+'.h5','r')
        ffea=h5py.File(folder_name+'/'+feature_name+'/'+feature_metadata_file+'chr'+str(chr)+'.h5','r')
        fsnp=h5py.File(folder_name+'/'+feature_name+'/'+snp_metadata_file+'chr'+str(chr)+'.h5','r')
    
        if ichr==0:
            fOut=h5py.File(folder_name+'/'+feature_name+'/'+file_name_feature+'.h5','w')
        else:
            fOut=h5py.File(folder_name+'/'+feature_name+'/'+file_name_feature+'.h5','r+')
        
        for feature in frez.keys():

            try: 
                fg=fOut.create_group(feature)
                fgm=fg.create_group('metadata')
                for key in ffea[feature].keys():
                    fgm.create_dataset(key,data=ffea[feature][key])
                
                fgd=fg.create_group('data')
                for key in frez[feature].keys():
                    fgd.create_dataset(key,data=frez[feature][key])
                                               
                fgs=fg.create_group('summary_data')
                feature_pv=frez[feature]['pv'][:]
                fgs.create_dataset('min_pv',data=np.nanmin(feature_pv) [None])
                fgs.create_dataset('min_pv_local_adjusted',data=np.nanmin(local_adjustment(feature_pv,method=local_adjustment_method))[None])
                fgs.create_dataset('min_pv_snp_id',data=np.array(frez[feature]['snp_id'][:][np.nanargmin(feature_pv)].astype('S'))[None])
            except: print (feature+' is duplicated ')
        frez.close()
        ffea.close()
        fsnp.close()
        fOut.close()
        
    fIn=h5py.File(folder_name+'/'+feature_name+'/'+file_name_feature+'.h5','r')
    fOut=h5py.File(folder_name+'/'+feature_name+'/'+file_name_gene+'.h5','w')
    genes=np.array([fIn[feature]['metadata/gene_id'][0] for feature in fIn.keys()])
    features=np.array([feature for feature in fIn.keys()])
    
    for gene in np.unique(genes):
        
        gene_features=features[genes==gene]
        fg=fOut.create_group(gene)
        
        fgm=fg.create_group('metadata')
        for key in fIn[gene_features[0]]['metadata'].keys():
            fgm.create_dataset(key,data=fIn[gene_features[0]]['metadata'][key])
            
        fgs=fg.create_group('summary_data')
        for key in ['min_pv','min_pv_local_adjusted']:
            fgs.create_dataset(key,data=np.nanmin(np.array([fIn[feature+'/summary_data/'+key][:] for feature in gene_features]))[None])    
        fgs.create_dataset('min_pv_feature_id',data=np.array(gene_features[np.nanargmin(np.array([fIn[feature+'/summary_data/min_pv_local_adjusted'][:] for feature in gene_features]))]).astype('S')[None])

    fIn.close()
    fOut.close()
        
    return [folder_name+feature_name+file_name_feature+'.h5',folder_name+feature_name+file_name_gene+'.h5']
 
summary_gene_feature(chr_list=np.arange(1,3), folder_name='/Users/mirauta/Scripts/QTL_pipeline_plots/', feature_name='Exon_1234',local_adjustment_method='')
 