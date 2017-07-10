import sys
import h5py
import numpy as np
import limix.stats.fdr as FDR
import pandas

def local_adjustment(pv, N=1,  method=''):
    if method is None:
        return pv
    if method=='Bonferroni': 
        N=np.hstack(pv).shape[0]
        return pv*N
    else:
        print('Valid multiple testing correction  methods are: None;Bonferroni')


    
def summary_gene_feature(qtl_results_file='qtl_results_',snp_metadata_file='snp_metadata_', feature_metadata_file='feature_metadata_',output_file='qtl_results_genome',\
                         feature_report='ensembl_gene_id',chr_list=[9],folder_data=None,trait=None, \
                            p_value_field='p_value',p_value_raw_field='p_value',local_adjustment_method='Bonferroni'):

    _doc=" aggregates qtl results to feature_report level"

    for ichr,chr in enumerate(chr_list):
        print ('chromosome: '+str(chr))
    
        try:
            frez=h5py.File(folder_data+'/'+trait+'/'+qtl_results_file+str(chr)+'.h5','r')
            frezkeys= np.array([k.replace('_i_','') for k in list(frez.keys())])
            ffea= pandas.read_table(folder_data+'/'+trait+'/'+feature_metadata_file+ str(chr)+'.txt', sep='\t')
            fsnp= pandas.read_table(folder_data+'/'+trait+'/'+snp_metadata_file+ str(chr)+'.txt', sep='\t').set_index('snp_id',drop=False).transpose()
        except:
            print('chromosome'+str(chr)+' missing')
            continue
        indexnan=np.where((ffea[feature_report])!=(ffea[feature_report]))[0]
        for i in indexnan:
            ffea[feature_report][i]='gene_'+str(ffea['chromosome'][i])+'_'+str(ffea['start'][i])
        ffea_feature=ffea.set_index('feature_id', drop=False).transpose()
        ffea_report_feature=ffea.set_index(feature_report, drop=False).transpose()
        
                             
        if ichr==0:
            fOut=h5py.File(folder_data+'/'+trait+'_'+feature_report+'_'+output_file+'.h5','w')
        else:
            fOut=h5py.File(folder_data+'/'+trait+'_'+feature_report+'_'+output_file+'.h5','r+')
    
        # for each report_feature  create h5 groups
        count=0
        for ifea,report_feature in enumerate(np.unique(list(ffea_report_feature))):
#            print (report_feature)
       
            #select features for which qtl was computed
            features=np.intersect1d(np.array( ffea_report_feature[report_feature].transpose()['feature_id']), frezkeys)
            if len(features) >=1:
                
                fg=fOut.create_group(report_feature)
                
                pv= np.array([frez[f][p_value_field] for f in  features ])
                for i, ppv in  enumerate(pv): ppv[ppv!=ppv]=1
                pv2= np.array([frez[f][p_value_raw_field] for f in  features ])
                for i, ppv in  enumerate(pv2): ppv[ppv!=ppv]=1
                beta= np.array([frez[f]['beta'] for f in  features ])
                for i, b in  enumerate(beta): b [b!=b]=1
                snp_id= np.array([frez[f]['snp_id'] for f in  features ])
                position= np.array([fsnp[snp_id[indf].astype('U')].transpose()['position'] for indf, f in  enumerate(features) ])
                for i, p in  enumerate(position): p [p!=p]=1
    
                fgm=fg.create_group('metadata')
                for key in ffea_feature[features[0]].keys():
         
                    if isinstance(ffea_feature[features[0]][key],int) :
                        fgm.create_dataset(key,data=np.array([ffea_feature[f][key] for f in  features ]))
                    else:
                        fgm.create_dataset(key,data=np.array([ffea_feature[f][key] for f in  features ]).astype('S'))
                
                fgd=fg.create_group('data')
                fgd.create_dataset('features',data= features.astype('S'))
                fgdp=fgd.create_group('p_value') 
                for indf,f in enumerate(features): fgdp.create_dataset(f,data= pv[indf])

                fgdp2=fgd.create_group('p_value_raw')
                for indf,f in enumerate(features): fgdp2.create_dataset(f,data= pv2[indf])
                
                fgdb=fgd.create_group('beta') 
                for indf,f in enumerate(features): fgdb.create_dataset(f,data=beta[indf])
                fgdpo=fgd.create_group('position') 
                for indf,f in enumerate(features): fgdpo.create_dataset(f,data=position[indf].astype(int))
                fgds=fgd.create_group('snp_id') 
                for indf,f in enumerate(features): fgds.create_dataset(f,data=snp_id[indf])
                                                                          
                                                 
                fgs=fg.create_group('summary_data')
                fgs.create_dataset('min_p_value',data= np.nanmin(np.hstack(pv)) [None])
                p_bonf=np.nanmin(local_adjustment(np.hstack(pv),method=local_adjustment_method));
                if p_bonf>1:p_bonf=np.array(1)
                fgs.create_dataset('min_p_value_local_adjusted',data=p_bonf[None])
                minfeatureindex=np.argmin([np.nanmin(ppv) *len(ppv) for ppv in pv])
                
                fgs.create_dataset('min_p_value_feature_id',data= np.array(fgm ['feature_id'][minfeatureindex].astype('S'))[None])
        
                min_snp_id=frez[features[minfeatureindex]]['snp_id'][np.nanargmin(pv[minfeatureindex])].astype('U')
                fgs.create_dataset('min_p_value_snp_id',data= np.array(min_snp_id).astype('S')[None])
                fgs.create_dataset('min_p_value_position',data= np.array(fsnp[min_snp_id]['position'])[None])
        
                count+=1

        frez.close()
        fOut.close()

def replication_two_features(folder_data ='/Users/mirauta/Data/MS/hipsci/TMT/',  folder_data2=None,  traits=['peptide_test','protein_test'],
    qtl_results_file='qtl_results_',    snp_metadata_file='snp_metadata_',    feature_metadata_file='feature_metadata_',
    results_genome_file='qtl_results_genome',    feature_report='ensembl_gene_id',p_value_field='p_value'):
    
    _doc=" aggregates qtl results from two traits at feature_report level; return replication of pvalues for trait1  signigicant snps in trait2 "

    if folder_data2 is None:
        folder_data2 =folder_data 
    featureh5=[h5py.File(folder_data+'/'+traits[0]+'_'+feature_report+'_'+results_genome_file+'.h5','r'),h5py.File(folder_data2+'/'+traits[1]+'_'+feature_report+'_'+results_genome_file+'.h5','r')]
    
    feature_ids=np.intersect1d(list(featureh5[0].keys()),list(featureh5[1].keys()))
 
    rez={}
    rez[p_value_field]=np.zeros(len(feature_ids))+np.nan
    rez['replicated_p_value']=np.zeros(len(feature_ids))+np.nan
    rez['replicated_self_p_value']=np.zeros(len(feature_ids))+np.nan
    rez['feature_id']=np.zeros(len(feature_ids),dtype='|S32')
    rez['gene_name']=np.zeros(len(feature_ids),dtype='|S32')
    rez['snp_id']=np.zeros(len(feature_ids),dtype='|S32')
    rez['chromosome']=np.zeros(len(feature_ids),dtype='|S5')
    rez['strand']=np.zeros(len(feature_ids),dtype='|S5')
    rez['position']=np.zeros(len(feature_ids))+np.nan
 
    for indf, feature in enumerate(feature_ids):
    
        temp=featureh5[0][feature]['summary_data/min_p_value'][0]
        if temp<0.001:
            rez[p_value_field][indf]=temp
            rez['feature_id'][indf]=featureh5[0][feature]['summary_data/min_p_value_feature_id'][:][0]
            rez['chromosome'][indf]=featureh5[0][feature]['metadata/chromosome'][:][0]
            rez['strand'][indf]=featureh5[0][feature]['metadata/feature_strand'][:][0]
            rez['gene_name'][indf]=featureh5[0][feature]['metadata/gene_name'][:][0]
            
            rez['snp_id'][indf]=featureh5[0][feature]['summary_data/min_p_value_snp_id'][:][0]
            rez['position'][indf]=featureh5[0][feature]['summary_data/min_p_value_position'][:][0]
            rez['replicated_p_value'][indf]=np.min( [featureh5[1][feature]['data/p_value'][f1][:][featureh5[1][feature]['data/snp_id'][f1][:]==featureh5[0][feature]['summary_data/min_p_value_snp_id'][0]] for f1 in featureh5[1][feature]['data/features']])     
            
            try:
                rez['replicated_self_p_value'][indf]=np.min([featureh5[0][feature]['data/p_value'][f1][:][featureh5[0][feature]['data/snp_id'][f1][:]==featureh5[0][feature]['summary_data/min_p_value_snp_id'][0]]\
                    for f1 in np.setdiff1d(featureh5[0][feature]['data/features'],featureh5[0][feature]['summary_data/min_p_value_feature_id'][0])]) 
            except:
                print (feature)
                1
    for f in featureh5: f.close()
    rez['ensembl_gene_id']=feature_ids.astype('U')
    for key in ['feature_id','chromosome','strand','snp_id','gene_name']:
        rez[key]=rez[key].astype('U')
    rez['traits']=traits
    return rez
 
