from matplotlib.patches import Rectangle
import matplotlib.gridspec as gridspec
import sys

from matplotlib_venn import venn2, venn2_circles
from matplotlib_venn import venn3, venn3_circles
import matplotlib.pyplot as plt
import seaborn as sb
from matplotlib import colors as mcolors
import matplotlib as mpl
import h5py
import numpy as np
#
#import limix.stats.fdr as FDR
import os
import scipy.stats as scst
import statsmodels.sandbox.stats.multicomp as scst2
sys.path.append('../')
from scripts.postprocess_functions_genome_results import *


#
#
#
#folder_name='/Users/mirauta/Results/hipsci/QTL1/'
#file_name_qtl='qtl_results_'
#file_name_perm='perm_results_'
#traits=['param_protein_scaled_peer_gaussnorm_test'];
#chromosome='22'
def plot_summary_onechr_perm(plot_name='qtl_summary_onechr',folder_name=None,folder_destination=None,\
                 traits=[''],chromosome='21',
                 file_name_qtl='qtl_results_',file_name_perm='perm_results_', \
                 colors=np.array(['orange','darkblue','green','m']),cis=2.5*10**5, figsize=(8,5)):
    
    if folder_destination is None:
        folder_destination =folder_name
    if not os.path.exists(folder_destination):
        os.makedirs(folder_destination)
 
    featureh5=[h5py.File(folder_name+'/'+trait+'/'+file_name_qtl+chromosome+'.h5','r') for trait in traits]
    featurepermh5=[h5py.File(folder_name+'/'+trait+'/'+file_name_perm+chromosome+'.h5','r') for trait in traits]
    gene_list_common=np.array([np.array(list(fh5.keys()))  for indf,fh5 in enumerate(featureh5)])
#        return gene_list_common
    temp=np.unique(np.hstack(gene_list_common),return_counts=1)
    gene_list_common=temp[0][temp[1]==gene_list_common.shape[0]]
#    local_adjusted_common=np.array([np.array([fh5[gene]['summary_data/min_p_value_local_adjusted'][:][0] for gene in np.intersect1d(gene_list_common,np.array(list(fh5.keys())))])for indf,fh5 in enumerate(featureh5)])
    fh5=featurepermh5[0]
    perm_fields=np.array(fh5[list(fh5.keys())[0]].dtype.names)[['permutation' in field for field in np.array(fh5[list(fh5.keys())[0]].dtype.names)]]
    
    
    rez={}
    
    rez['p_value']=np.array([np.hstack([fh5[gene]['p_value']for gene in gene_list_common])  for indf,fh5 in enumerate(featureh5)])
    for field in perm_fields:
        rez[field]=np.array([np.hstack([fh5[gene][field] for gene in gene_list_common])  for indf,fh5 in enumerate(featurepermh5)])
    rez['p_value_min_bonferroni']=np.array([np.hstack([np.nanmin(fh5[gene]['p_value'])*fh5[gene]['p_value'].shape[0] for gene in gene_list_common])  for indf,fh5 in enumerate(featureh5)])[0]
    for field in perm_fields:
        rez[field+'_min_bonferroni']=np.array([np.hstack([np.nanmin(fh5[gene][field])*fh5[gene][field].shape[0] for gene in gene_list_common])  for indf,fh5 in enumerate(featurepermh5)])[0]
#    else:
#        local_adjusted=np.array([np.array([fh5[gene]['summary_data/min_p_value_local_adjusted'][:][0] for gene in np.intersect1d(gene_list,np.array(list(fh5.keys())))])for indf,fh5 in enumerate(featureh5)])

 
    fig=plt.figure(figsize=figsize)
    mpl.rcParams.update(mpl.rcParamsDefault)
    fig.patch.set_facecolor('white')
    axes = fig.add_subplot(1, 2, 1, axisbg='white')
    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)
    axes.yaxis.set_ticks_position('left')
    axes.xaxis.set_ticks_position('bottom')
    fdr=np.array([5,4,3,2,1])[::-1]
    trait_labels=['qval']+list(perm_fields)
    local_adjusted=[rez['p_value_min_bonferroni']]+[rez[field+'_min_bonferroni'] for field in perm_fields]
    for iy,yy in enumerate(local_adjusted):
#        print(yy)
 
#        plt.plot(fdr[::-1] ,[(-np.log10(FDR.qvalues1(yy))>thr).sum() for thr in fdr],color=colors[iy],markersize=3,label=trait_labels[iy])
        plt.plot(fdr[::-1] ,[(-np.log10(scst2.fdrcorrection0(yy)[1])>thr).sum() for thr in fdr],color=colors[iy],label=trait_labels[iy],lw=3)
    axes.plot((fdr[-1],fdr[-1]),(0,np.max([(a<10**-2).sum() for a in local_adjusted])) ,'--',color='k',lw=0.5)
    axes.plot((fdr[0],fdr[-1]),(np.max([(a<10**-2).sum() for a in local_adjusted]),np.max([(a<10**-2).sum() for a in local_adjusted])) ,'--',color='k',lw=0.5)
    plt.xlabel('-log10 FDR',fontsize=13);
    plt.ylabel('#pGenes',fontsize=13)
    plt.xticks(fdr,fdr[::-1])
    if trait_labels is not None: plt.legend(loc=2,fontsize=9)
    
#==============================================================================
#     power plot commmon
#==============================================================================
    plt.subplot(1,2,2)
    
    plt.title('Calibration',fontsize=10)
    
    for iy,yy in enumerate([rez['p_value']]+[rez[field] for field in perm_fields]):
        yy=np.sort(yy[(yy==yy)&(yy!=1)])
        plt.plot(-np.log10((0.5+np.arange(len(yy)))/len(yy))[::-1],-np.log10(yy)[::-1],'o',color=colors[iy],markersize=4,label=trait_labels[iy])
                
    plt.plot(-np.log10((0.5+np.arange(len(yy)))/len(yy))[::-1],-np.log10((0.5+np.arange(len(yy)))/len(yy))[::-1],'k',lw=1)
    plt.xlabel('-log10 PV',fontsize=13);
    plt.ylabel('-log10 random',fontsize=13)
    plt.legend(loc=2,fontsize=8)
    for f in featureh5: f.close()
    plt.tight_layout()
    plt.savefig(folder_destination+plot_name+'_'+'_'.join(trait_labels)+'.png',dpi=600)
    return [gene_list_common,rez]


def plot_summary(plot_name='qtl_summary',folder_name=None,folder_destination=None,\
                 traits=[''],trait_labels=None,
                 file_name_results_genome='Feature_results_genome', qtl_results_file='QTL_results_',\
                 colors=np.array(['orange','darkblue','green','m']),cis=2.5*10**5, figsize=(12,12), gene_list=None,\
                 plot_tss_distance_flag=False,plot_calibration_flag=False):
    
    if folder_destination is None:
        folder_destination =folder_name
    if not os.path.exists(folder_destination):
        os.makedirs(folder_destination)
    
    if trait_labels is None: trait_labels=traits
#    
    fig=plt.figure(figsize=figsize)
    mpl.rcParams.update(mpl.rcParamsDefault)
    fig.patch.set_facecolor('white')
 
    featureh5=[h5py.File(folder_name+'/'+trait+'_'+file_name_results_genome+'.h5','r') for trait in traits]
 
    if gene_list is None:
        local_adjusted=np.array([np.array([fh5[gene]['summary_data/min_p_value_local_adjusted'][:][0] for gene in fh5.keys()]) \
                                 for indf,fh5 in enumerate(featureh5)])
        gene_list_common=np.array([np.array(list(fh5.keys()))  for indf,fh5 in enumerate(featureh5)])
        temp=np.unique(np.hstack(gene_list_common),return_counts=1)
        gene_list_common=temp[0][temp[1]==gene_list_common.shape[0]]
        local_adjusted_common=np.array([np.array([fh5[gene]['summary_data/min_p_value_local_adjusted'][:][0] for gene in np.intersect1d(gene_list_common,np.array(list(fh5.keys())))])for indf,fh5 in enumerate(featureh5)])
        
 
    else:
        local_adjusted=np.array([np.array([fh5[gene]['summary_data/min_p_value_local_adjusted'][:][0] for gene in np.intersect1d(gene_list,np.array(list(fh5.keys())))])for indf,fh5 in enumerate(featureh5)])

 
 
    axes = fig.add_subplot(2, 2, 1, axisbg='white')
    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)
    axes.yaxis.set_ticks_position('left')
    axes.xaxis.set_ticks_position('bottom')
    fdr=np.array([5,4,3,2,1])[::-1]
 
    for iy,yy in enumerate(local_adjusted):
#        print(yy)
 
#        plt.plot(fdr[::-1] ,[(-np.log10(FDR.qvalues1(yy))>thr).sum() for thr in fdr],color=colors[iy],markersize=3,label=trait_labels[iy])
        plt.plot(fdr[::-1] ,[(-np.log10(scst2.fdrcorrection0(yy)[1])>thr).sum() for thr in fdr],color=colors[iy],label=trait_labels[iy],lw=3)
    axes.plot((fdr[-1],fdr[-1]),(0,np.max([(a<10**-2).sum() for a in local_adjusted])) ,'--',color='k',lw=0.5)
    axes.plot((fdr[0],fdr[-1]),(np.max([(a<10**-2).sum() for a in local_adjusted]),np.max([(a<10**-2).sum() for a in local_adjusted])) ,'--',color='k',lw=0.5)
    plt.xlabel('-log10 FDR',fontsize=13);
    plt.ylabel('#pGenes',fontsize=13)
    plt.xticks(fdr,fdr[::-1])
    if trait_labels is not None: plt.legend(loc=2,fontsize=9)
    
#==============================================================================
#     power plot commmon
#==============================================================================
    axes = fig.add_subplot(2, 2, 2, axisbg='white')
    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)
    axes.yaxis.set_ticks_position('left')
    axes.xaxis.set_ticks_position('bottom')
    fdr=np.array([5,4,3,2,1])[::-1]
 
    for iy,yy in enumerate(local_adjusted_common):
#        print(yy)
 
#        plt.plot(fdr[::-1] ,[(-np.log10(FDR.qvalues1(yy))>thr).sum() for thr in fdr],color=colors[iy],markersize=3,label=trait_labels[iy])
        plt.plot(fdr[::-1] ,[(-np.log10(scst2.fdrcorrection0(yy)[1])>thr).sum() for thr in fdr],color=colors[iy],label=trait_labels[iy],lw=3)
    axes.plot((fdr[-1],fdr[-1]),(0,np.max([(a<10**-2).sum() for a in local_adjusted_common])) ,'--',color='k',lw=0.5)
    axes.plot((fdr[0],fdr[-1]),(np.max([(a<10**-2).sum() for a in local_adjusted_common]),np.max([(a<10**-2).sum() for a in local_adjusted_common])) ,'--',color='k',lw=0.5)
    plt.xlabel('-log10 FDR',fontsize=13);
    plt.ylabel('#pGenes',fontsize=13)
    plt.xticks(fdr,fdr[::-1])
    plt.title('Common genes')
    if trait_labels is not None: plt.legend(loc=2,fontsize=9)


    if plot_tss_distance_flag:
        axes = fig.add_subplot(2, 2, 3, axisbg='white')
        axes.spines['top'].set_visible(False)
        axes.spines['right'].set_visible(False)
        axes.yaxis.set_ticks_position('left')
        axes.xaxis.set_ticks_position('bottom')
        def get_distance_tss(gene,fh5):
            if ( fh5[gene]['metadata/feature_strand'][:].astype('U')[0]=='+'):
                return fh5[gene]['summary_data/min_p_value_position'][:][0]-fh5[gene]['metadata/start'][:][0]
            else:
                return -fh5[gene]['summary_data/min_p_value_position'][:][0]+fh5[gene]['metadata/end'][:][0]
    
        if gene_list is None:
            distance_tss=  [np.array([get_distance_tss(gene,fh5) for gene in  np.array(list(fh5.keys()))[ local_adjusted[indf]<10**-3]]) for indf,fh5 in enumerate(featureh5)]
        else:
            distance_tss=  [np.array([get_distance_tss(gene,fh5) for gene in np.intersect1d(gene_list,np.array(list(fh5.keys())))[ local_adjusted[indf]<10**-3]]) for indf,fh5 in enumerate(featureh5)]
        axes.hist([np.clip(d[np.isfinite(d)],-cis,cis) for d in distance_tss], bins=7,width=15000,label=trait_labels,color=colors[:len(traits)])
        plt.ylabel('#pGenes \n(PV<0.001)',fontsize=13);plt.xlabel('distance from TSS',fontsize=13);
        plt.legend(loc=2,fontsize=9)
        plt.xticks(np.linspace(-cis,cis,5),np.array([str(int(l))+'k' for l in np.linspace(-cis,cis,5)/1e3]))
            

    if plot_calibration_flag:
        plt.subplot(2,2,4)
    
        featureh5=[h5py.File(folder_name+'/'+trait+'_'+file_name_results_genome+'.h5','r') for trait in traits]
    
        y=[np.hstack([np.hstack([fh5[feature_id]['data/p_value_raw'][f][:] for f in fh5[feature_id]['data/p_value_raw'].keys()]).flatten()  for feature_id in list(fh5.keys())]) for indf,fh5 in enumerate(featureh5)]

        plt.title('Calibration',fontsize=10)
        
        for iy,yy in enumerate(y):
            yy=np.sort(yy[(yy==yy)&(yy!=1)])
            try:
                plt.plot(-np.log10((0.5+np.arange(len(yy)))/len(yy))[::-1],-np.log10(yy)[::-1],'o',color=colors[iy],markersize=4,label=feature_label[iy])
            except:
                plt.plot(-np.log10((0.5+np.arange(len(yy)))/len(yy))[::-1],-np.log10(yy)[::-1],'o',color=colors[iy],markersize=4)
                    
        plt.plot(-np.log10((0.5+np.arange(len(yy)))/len(yy))[::-1],-np.log10((0.5+np.arange(len(yy)))/len(yy))[::-1],'k',lw=1)
        plt.legend(loc=2,fontsize=8)

    for f in featureh5: f.close()
    plt.tight_layout()
    plt.savefig(folder_destination+plot_name+'_'+'_'.join(trait_labels)+'.png',dpi=600)
    return local_adjusted



def plot_manhatan_alone(folder_name='/Users/mirauta/Data/MS/hipsci/TMT/',folder_destination='/Users/mirauta/Data/MS/hipsci/TMT/Images',plot_name='manhattan',\
                        traits=None,trait_labels=None,file_name_results_genome='ensembl_gene_id_qtl_results_genome',   qtl_results_file='qtl_results_',colors=np.array(['b','k','g','m']), figsize=4, gene_ensembl_id= 'ENSG00000182154',\
                        p_value_field='p_value'):
    if folder_destination is None:
        folder_destination =folder_name+'/manhattan/'
    if not os.path.exists(folder_destination):
        os.makedirs(folder_destination)

    featureh5=[h5py.File(folder_name+'/'+trait+'_'+file_name_results_genome+'.h5','r')[gene_ensembl_id] for trait in traits]
    
    rez={}
    temppos=[np.hstack([fh5['data']['position'][f][:] for f in fh5['data']['position'].keys()]) for fh5 in featureh5]
    temp=np.unique(np.hstack(temppos).flatten(),return_counts=1); 
    commonpos=temp[0][temp[1]==np.max(temp[1])]
    
    for dat in [p_value_field,'position']:
        rez[dat]=[np.array([fh5['data'][dat][f][:][np.in1d(fh5['data']['position'][f][:],commonpos)] for f in fh5['data'][dat].keys()]) for fh5 in featureh5]
 
    #==============================================================================
    #  modifty positions for plot
    #==============================================================================
    xpos=(commonpos-commonpos.min())/10.0**6;
    axxpos=np.array([i for i in np.linspace(min(xpos),max(xpos),5)])
    showxpos=np.array([int(i) for i in np.linspace(min(commonpos),max(commonpos),5)])
    startstop=np.array([featureh5[0]['metadata']['start'][:][0],featureh5[0]['metadata']['end'][:][0]]);
    startstop=(startstop-commonpos.min())/10.0**6
    
    fig=plt.figure(figsize=(figsize*2,figsize*len(featureh5)))
    f = gridspec.GridSpec(len(featureh5)*6,8)
    f.update(hspace=10, wspace=10)
    ax = [plt.subplot(f[(i*3):(i*3+3),:6]) for i in np.arange(len(featureh5))]
    fig.set_facecolor('white')
    colors2=colors[1:]
    for indf, a in enumerate(ax):
        a.spines['top'].set_visible(False);    a.spines['right'].set_visible(False);a.yaxis.set_ticks_position('left'); a.xaxis.set_ticks_position('bottom');a.set_axis_bgcolor('white');
        a.add_patch(Rectangle((startstop[0], 0), startstop[1]-startstop[0], 5, facecolor="grey",    alpha=0.35))
        a.set_xticks( axxpos)
        a.set_xticklabels([str(np.around(i,1))+' Mb' for i in showxpos/10.0**6])
        for indt,t in  enumerate(rez[p_value_field][indf][np.argsort([tt.min() for tt in rez[p_value_field][indf]])[:4]]):
            a.plot(xpos,-np.log10(t),'o',color=colors[indt],markersize=1.25)
        a.set_ylabel(trait_labels[indf]+"            \n QTL              \n -log10PV            ",fontsize=10,rotation=0,horizontalalignment= 'center' ,verticalalignment= 'center')
    
        for indt,t in enumerate(rez['p_value'][indf][np.argsort([tt.min() for tt in rez['p_value'][indf]])[:4]]): 
            a.plot(xpos[np.argsort(t)[:5]],-np.log10(t[np.argsort(t)[:5]]),'*',color=colors[indt],markersize=5)
            for indf2, a2 in enumerate(ax):
                if indf!=indf2:
                    for indt2,t2 in enumerate(rez[p_value_field][indf2][np.argsort([tt2.min() for tt2 in rez['p_value'][indf2]])[:4]]):
                        a.plot(xpos[np.argsort(t2)[:3]],-np.log10(t[np.argsort(t2)[:3]],),'ro',markersize=2.5)
     
    print (str(featureh5[0]['metadata']['gene_name'][:].astype('U')[0]))
    ax[0].annotate(str(featureh5[0]['metadata']['gene_name'][:].astype('U')[0]), xy=(startstop[1], 4),fontsize=11)
 
    plt.savefig(folder_destination+plot_name+str(featureh5[0]['metadata']['gene_name'][:].astype('U')[0])+'.pdf',dpi=600)

#==============================================================================
#==============================================================================
# # 
#==============================================================================
#==============================================================================
 

def plot_replication_beta(rez=None,folder_data ='/Users/mirauta/Data/MS/hipsci/TMT/', folder_data2 =None,   traits=['protein_test','peptide_test'],trait_labels=['protein_test','peptide_test'],\
    qtl_results_file='qtl_results_',    snp_metadata_file='snp_metadata_',    feature_metadata_file='feature_metadata_',\
    results_genome_file='qtl_results_genome',    feature_report='ensembl_gene_id',folder_destination='/Users/mirauta/Results/hipsci/Images_pipeline',\
    figsize=5, red_dots_features=None,red_dot='ro',plot_name=''):
    if rez is None:
        rez=replication_two_features(folder_data =folder_data, folder_data2 =folder_data2,    traits=np.array(traits), qtl_results_file=qtl_results_file,    snp_metadata_file=snp_metadata_file,    feature_metadata_file=feature_metadata_file, results_genome_file=results_genome_file,    feature_report= feature_report)
   
    fig=plt.figure(figsize=(figsize,figsize))
    mpl.rcParams.update(mpl.rcParamsDefault)
    axes = fig.add_subplot(1, 1, 1, axisbg='white')
    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)
    axes.spines['top'].set_visible(False)
    axes.spines['left'].set_visible(True)
    axes.yaxis.set_ticks_position('left')
    axes.xaxis.set_ticks_position('bottom')
    fig.patch.set_facecolor('white')
    plt.plot((rez['beta']),(rez['replicated_beta']),'o',color='grey',markersize=6,label=scst.spearmanr(rez['beta'][np.isfinite(rez['replicated_beta']+rez['beta'])],rez['replicated_beta'][np.isfinite(rez['replicated_beta']+rez['beta'])])[0])
    thrs=0.001;
    plt.plot((rez['beta'][rez['replicated_self_beta']<thrs]),(rez['replicated_beta'][rez['replicated_self_beta']<thrs]),'.',markersize=8,color='darkorange')
    
    if red_dots_features is not None:        
        plt.plot((rez['beta'][np.in1d(rez['feature_id'],red_dots_features)]), (rez['replicated_beta'][np.in1d(rez['feature_id'],red_dots_features)]),red_dot,markersize=6)#, mfc='none')
        
    plt.plot((np.nanmin(rez['beta']),(np.nanmax(rez['beta']))),(np.nanmin(rez['beta']),(np.nanmax(rez['beta']))),'k--',lw=0.25)

    plt.xlabel(trait_labels[0]+'\n - log10 PV');plt.ylabel(trait_labels[1]+'\n - log10 PV',rotation=90 )
    plt.legend(loc=1)
    plt.savefig(folder_destination+'replication_'+traits[0]+'_'+traits[1]+plot_name+'.png')
    return rez

def plot_replication_pv(rez=None,folder_data ='/Users/mirauta/Data/MS/hipsci/TMT/', folder_data2 =None,   traits=['protein_test','peptide_test'],trait_labels=['protein_test','peptide_test'],\
    qtl_results_file='qtl_results_',    snp_metadata_file='snp_metadata_',    feature_metadata_file='feature_metadata_',\
    results_genome_file='qtl_results_genome',    feature_report='ensembl_gene_id',folder_destination='/Users/mirauta/Results/hipsci/Images_pipeline',\
    figsize=5, red_dots_features=None,red_dot='ro',plot_name=''):
    if rez is None:
        rez=replication_two_features(folder_data =folder_data, folder_data2 =folder_data2,    traits=np.array(traits), qtl_results_file=qtl_results_file,    snp_metadata_file=snp_metadata_file,    feature_metadata_file=feature_metadata_file, results_genome_file=results_genome_file,    feature_report= feature_report)
   
    fig=plt.figure(figsize=(figsize,figsize))
    mpl.rcParams.update(mpl.rcParamsDefault)
    axes = fig.add_subplot(1, 1, 1, axisbg='white')
    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)
    axes.spines['top'].set_visible(False)
    axes.spines['left'].set_visible(True)
    axes.yaxis.set_ticks_position('left')
    axes.xaxis.set_ticks_position('bottom')
    fig.patch.set_facecolor('white')
    plt.plot(-np.log10(rez['p_value']),-np.log10(rez['replicated_p_value']),'o',color='grey',markersize=6)
    thrs=0.001;
    plt.plot(-np.log10(rez['p_value'][rez['replicated_self_p_value']<thrs]),-np.log10(rez['replicated_p_value'][rez['replicated_self_p_value']<thrs]),'.',markersize=8,color='darkorange')
    
    if red_dots_features is not None:        
        plt.plot(-np.log10(rez['p_value'][np.in1d(rez['feature_id'],red_dots_features)]), -np.log10(rez['replicated_p_value'][np.in1d(rez['feature_id'],red_dots_features)]),red_dot,markersize=6)#, mfc='none')
        
    plt.plot((2,-np.log10(np.nanmin(rez['p_value']))),(2,-np.log10(np.nanmin(rez['p_value']))),'k--',lw=0.25)

    plt.xlabel(trait_labels[0]+'\n - log10 PV');plt.ylabel(trait_labels[1]+'\n - log10 PV',rotation=90 )

    plt.savefig(folder_destination+'replication_'+traits[0]+'_'+traits[1]+plot_name+'.png')
    return rez
