from matplotlib.patches import Rectangle
import matplotlib.gridspec as gridspec
import sys
from matplotlib_venn import venn2, venn2_circles
from matplotlib_venn import venn3, venn3_circles
import matplotlib.pyplot as plt
import seaborn as sb
from matplotlib import colors as mcolors
import h5py
import numpy as np
#
import limix.stats.fdr as FDR
import os
import scipy.stats as scst

def plot_summary(folder_name='/Users/mirauta/Scripts/QTL_pipeline_plots/',folder_destination='',plot_name='qtl_summary',\
                 feature_name_list=['Exon_1234'],feature_label=None,
                 file_name_results_genome='Feature_results_genome', qtl_results_file='QTL_results_',\
                 thresholds=[0.05,0.01],thrrep=0.05,colors=np.array(['b','r','g','m']),cis=2.5*10**5, figsize=(12,12)):
    
    if not os.path.exists(folder_destination):
        os.makedirs(folder_destination)
    
    if feature_label is None: feature_label=feature_name_list
    
    fig=plt.figure(figsize=figsize)
    #    mpl.rcParams.update(mpl.rcParamsDefault)
    axes = fig.add_subplot(1, 1, 1, axisbg='white')
    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)
    axes.yaxis.set_ticks_position('left')
    axes.xaxis.set_ticks_position('bottom')
    fig.patch.set_facecolor('white')
    
#    ###    #
    featureh5=[h5py.File(folder_name+feature+'/'+file_name_results_genome+'.h5','r') for feature in feature_name_list]
    y=[np.hstack([fh5[feature_id]['summary_data/min_p_value_local_adjusted'][:] for feature_id in list(fh5.keys())]) for indf,fh5 in enumerate(featureh5)]  
    for f in featureh5: f.close()
    plt.subplot(2,2,1)    
    fdr=np.array([7,6,5,4,3,2])[::-1]
    plt.title('Number of significant pGenes after local correction',fontsize=10)
    for iy,yy in enumerate(y):
        try:
            plt.plot(fdr[::-1] ,[(-np.log10(yy)>thr).sum() for thr in fdr],color=colors[iy],markersize=3,label=feature_label[iy])
        except:
            plt.plot(fdr[::-1] ,[(-np.log10(yy)>thr).sum() for thr in fdr],color=colors[iy],markersize=3)
    plt.xlabel('-log10 PV Bonferroni',fontsize=10);
    plt.ylabel('genes with a significant pQTL',fontsize=9)
    plt.xticks(fdr,fdr[::-1])
    if feature_label is not None: plt.legend(loc=2,fontsize=8)

    plt.subplot(2,2,2)
    fdr=np.array([5,4,3,2])[::-1]
    plt.title('Number of significant pGenes /Local Bonf correction/ global BH',fontsize=10)
    for iy,yy in enumerate(y):
        print(yy)
        plt.plot(fdr[::-1] ,[(-np.log10(FDR.qvalues1(yy))>thr).sum() for thr in fdr],color=colors[iy],markersize=3,label=feature_label[iy])
    plt.xlabel('-log10 QV',fontsize=10);
    plt.ylabel('genes with a significant pQTL',fontsize=9)
    plt.xticks(fdr,fdr[::-1])
    if feature_label is not None: plt.legend(loc=2,fontsize=8)
    y_test=np.copy(yy)+1-1

#    for re in summary_file: re.close()

    plt.subplot(2,2,3)
     
#    plt.hist([np.squeeze(ppv[da+'_tss'].values())[np.argsort(np.squeeze(ppv[da+'_Bonfeff'].values()))[:len(ppv[da+'_Bonfeff'].values())/20]] for ida,da in enumerate(datas2)], bins=7,width=15000,color=colors2[:len(datas2)],label=datas2)
#     
#    plt.ylabel('p-Genes (signif. PV)');plt.xlabel('distance from TSS');plt.title('distance from TSS',fontsize=7)
#    plt.legend(loc=2,fontsize=6)
#    plt.xticks(np.linspace(-2e5,2e5,5),np.array([str(int(l))+'k' for l in np.linspace(-2e5,2e5,5)/1e3]))
#    
## 
    plt.subplot(2,2,4)
    
    featureh5=[h5py.File(folder_name+feature+'/'+file_name_results_genome+'.h5','r') for feature in feature_name_list]
#    print (featureh5[0][list(featureh5[0].keys())[0]]['summary_data/p_value'])
#    sys.exit()
    print (">>>>>");     print (list(featureh5[0][np.array(list(featureh5[0].keys()))[0]]['metadata'].keys()))
    y=[np.hstack([fh5[feature_id]['data/p_value'][:] for feature_id in list(fh5.keys())]) for indf,fh5 in enumerate(featureh5)]  
    
    for f in featureh5: f.close()
    
    print (len(y))
    print (y[0].shape)
    plt.title('Calibration',fontsize=10)
    
    for iy,yy in enumerate(y):
        yy=np.sort(yy[(yy==yy)&(yy!=1)])
        try:
            plt.plot(-np.log10((0.5+np.arange(len(yy)))/len(yy))[::-1],-np.log10(yy)[::-1],colors[iy]+'o',markersize=4,label=feature_label[iy])
        except:
            plt.plot(-np.log10((0.5+np.arange(len(yy)))/len(yy))[::-1],-np.log10(yy)[::-1],colors[iy]+'o',markersize=4)
                
#    x=np.sort(np.hstack(ppv[da+'_pv0'].values()));x=x[(x==x)&(x!=1)]
#    plt.plot(-np.log10((0.5+np.arange(len(x)))/len(x))[::-1],-np.log10(x)[::-1],colors2[ida]+'--',markersize=0.51,label=da+'0')
    plt.plot(-np.log10((0.5+np.arange(len(yy)))/len(yy))[::-1],-np.log10((0.5+np.arange(len(yy)))/len(yy))[::-1],'k',lw=1)
    plt.legend(loc=2,fontsize=8)


    plt.tight_layout()
    plt.savefig(folder_destination+plot_name+'.png',dpi=600)
    return y_test
#    plt.show()

#featureh5= h5py.File(protein_test/gene_name_results_genome.h5','r')

#yy=plot_summary(folder_name='/Users/mirauta/Git/hipsci_pipeline/limix_QTL_pipeline/data/geuvadis_CEU_YRI_test_data/',folder_destination='/Users/mirauta/Data/MS/hipsci/TMT/Images',plot_name='qtl_summary',\
#                 feature_name_list=['protein_test_test','protein_test_perm','protein_test_perm2'],\
#                                   feature_label=['protein_test', 'protein_perm','protein_perm2'],
#                 file_name_results_genome='feature_id_qtl_results_genome',   qtl_results_file='QTL_results_',\
#                 thresholds=[0.05,0.01],thrrep=0.05,colors=np.array(['b','r','g','m']),cis=2.5*10**5,figsize=(12,12))
#(-np.log10(FDR.qvalues1(yy))>-0.0005).sum()
#(-np.log10(FDR.qvalues1(yy))>1).sum()
##
#
yy=plot_summary(folder_name='/Users/mirauta/Data/MS/hipsci/TMT/',folder_destination='/Users/mirauta/Data/MS/hipsci/TMT/Images',plot_name='qtl_summary',\
                 feature_name_list=['protein_test','protein_perm'],feature_label=['protein_test','protein_perm'],
                 file_name_results_genome='feature_id_qtl_results_genome',   qtl_results_file='QTL_results_',\
                 thresholds=[0.05,0.01],thrrep=0.05,colors=np.array(['b','r','g','m']),cis=2.5*10**5,figsize=(6,6))
(-np.log10(FDR.qvalues1(np.arange(1,100)/1000.0))>2).sum()
((-np.log10(FDR.qvalues(yy)))>11).sum()
plt.plot(yy)
#10**-1.301029995
#
#rez=h5py.File(nameh5,'r')
#testh5=h5py.File('/Users/mirauta/Scripts/QTL_pipeline_plots/test.h5','r')
 

#def plot_manhatan_alone(folder_destination='',plot_name='manhattan_',\
#                        gene,feature_chromosome_file=[[]],gene_metadata_genome =[],snp_metadata_genome =[],\
#                                                     featurelabel=None,colors=np.array(['b','r','g','m']),cis=2.5*10**5):
##def plot_manhatan_alone(folder_data,comment,gene_data_processed=None,gene_data=None,egene_data=None,dat_type=None,dat_type2=None,\
##                        peptides=None,peptide_intensity=None,figsize=5):
#
#    Nfeatures=len(feature_chromosome_file)
#    
#    try: genepos=snp_metadata_genome[][gene][]
#
#    if gene_data_processed is not None:
#        try:
#            genepqtlpv=gene_data_processed[dat_type+'_pv'][gene_data['ensembl_ID'][:][0]]
#
#         
#        except: 1
#    else:
#        genepqtlpv=gene_data['lmm/'+dat_type]['pv'][:][0][np.unique(genepos,return_index=1)[1]]
#        if dat_type2 is not None:
#            genepqtlpv2=gene_data['lmm/'+dat_type2]['pv'][:][0][np.unique(genepos,return_index=1)[1]]
#        genepqtlqv=gene_data['lmm/'+dat_type]['qv'][:][np.unique(genepos,return_index=1)[1]]
#    if peptides is not None: pepqtlpv=np.array([gene_data['lmm/'+pep]['pv'][:][0][np.unique(genepos,return_index=1)[1]] for pep in peptides])
#
#    genepos=genepos[np.unique(genepos,return_index=1)[1]]
#
#    egenepos= egene_data['pos'][:]
#    geneeqtlpv=egene_data['lmm_gene_peer_k30_vsn_covar_gauss_pv'][:][0][np.unique(egenepos,return_index=1)[1]]
#    geneeqtlqv=egene_data['lmm_gene_peer_k30_vsn_covar_gauss_qv'][:][np.unique(egenepos,return_index=1)[1]]
#    egenepos=egenepos[np.unique(egenepos,return_index=1)[1]]
#
#    geneeqtlpv,geneeqtlqv =[g[np.in1d(egenepos,genepos)] for g in [geneeqtlpv,geneeqtlqv]]
#    genepqtlpv,genepqtlqv =[g[np.in1d(genepos,egenepos)] for g in [genepqtlpv,genepqtlqv]]
#    if peptides is not None:
#        pepqtlpv =pepqtlpv[:,np.in1d(genepos,egenepos)]
#        genepqtlpv2 =genepqtlpv2[np.in1d(genepos,egenepos)]
#    genepos=genepos[np.in1d(genepos,egenepos)]
#
#    xpos=(genepos-genepos.min())/10.0**6;
#    axxpos=np.array([i for i in np.linspace(min(xpos),max(xpos),5)])
#    showxpos=np.array([int(i) for i in np.linspace(min(genepos),max(genepos),5)])
#    startstop=np.array([gene_data['gene_start'][:].astype(float)[0],gene_data['gene_end'][:].astype(float)[0]]);startstop=(startstop-genepos.min())/10.0**6
#    if peptides is not None:
#        fig=plt.figure(figsize=(figsize*2,figsize*1.5))
#        f = gridspec.GridSpec(8, 6)
#        f.update(hspace=.8, wspace=1.05)
#        a0 = plt.subplot(f[:2,:4])
#        a1 = plt.subplot(f[2:4,:4])
#        a3 = plt.subplot(f[4:6,:4])
#        a2 = plt.subplot(f[6:,:4])
#
#        a01= plt.subplot(f[:2,4:])
#        a03= plt.subplot(f[4:6,4:])
#    else:
#        fig=plt.figure(figsize=(figsize*2,figsize))
#        f = gridspec.GridSpec(6, 4)
#        f.update(hspace=.8, wspace=1.05)
#        a0 = plt.subplot(f[:3,:4])
#        a2 = plt.subplot(f[3:,:4])
#    fig.set_facecolor('white')
#
#    for a in [a0,a2]:
#        a.spines['top'].set_visible(False);    a.spines['right'].set_visible(False);a.yaxis.set_ticks_position('left'); a.xaxis.set_ticks_position('bottom');a.set_axis_bgcolor('white');
#        a.add_patch(Rectangle((startstop[0], 0), startstop[1]-startstop[0], np.max(-np.log10(genepqtlpv)), facecolor="grey",    alpha=0.35))
#        a.set_xticks( axxpos)
#        a.set_xticklabels([str(np.around(i,1))+' Mb' for i in showxpos/10.0**6])
#
#    if peptides is not None:
#        for a in [a1,a3]:
#            a.spines['top'].set_visible(False);    a.spines['right'].set_visible(False);a.yaxis.set_ticks_position('left'); a.xaxis.set_ticks_position('bottom');  a.set_axis_bgcolor('white');
#            a.set_xticks( axxpos); a.set_xticklabels([str(np.around(i,1))+' Mb' for i in showxpos/10.0**6]);
#            a.add_patch(Rectangle((startstop[0], 0), startstop[1]-startstop[0], np.max(-np.log10(genepqtlpv)), facecolor="grey",alpha=0.35))
#        a3.set_ylabel("peptide\n QTL\n -log10PV ",fontsize=12,rotation=0,horizontalalignment= 'right')
#        a1.set_ylabel("protein model_2\n QTL\n -log10PV ",fontsize=12,rotation=0,horizontalalignment= 'right')
#
#    a0.plot(xpos,-np.log10(genepqtlpv),'ko',markersize=1.5)
#    if (genepqtlqv<0.01).sum()>0:
#        a0.plot((np.min(xpos),np.max(xpos)),(-np.log10(np.max(genepqtlpv[genepqtlqv<0.01])),-np.log10(np.max(genepqtlpv[genepqtlqv<0.01]))),'k--',lw=0.5)
#        a0.annotate("FDR 0.01", xy=(np.min(xpos)+0.05,-np.log10(np.max(genepqtlpv[genepqtlqv<0.01]))+0.05 ),fontsize=10)
#    a0.annotate(gene, xy=(startstop[1]+0.05, -0.5+np.max(-np.log10(genepqtlpv))),fontsize=14)
#    a0.plot(xpos[-np.log10(genepqtlpv)>np.percentile(-np.log10(genepqtlpv),99.5)],-np.log10(genepqtlpv)[-np.log10(genepqtlpv)>np.percentile(-np.log10(genepqtlpv),99.5)],'r*',markersize=8)
#    a0.plot(xpos[-np.log10(geneeqtlpv)>np.percentile(-np.log10(geneeqtlpv),99.5)],-np.log10(genepqtlpv)[-np.log10(geneeqtlpv)>np.percentile(-np.log10(geneeqtlpv),99.5)],'go',markersize=5)
#    a0.set_ylabel("pQTL\n -log10PV ",fontsize=12,rotation=0,horizontalalignment= 'right')
#
#    if peptides is not None:
#        colorspep0=np.array(['b','m','y','c','g','r']*10)
#        colorspep0=np.array(mcolors.cnames.keys())[np.array([('pink' not in p)&('white' not in p)&('light' not in p) for p in np.array(mcolors.cnames.keys())])]
#        box=a01.boxplot(np.log(peptide_intensity+1).T,patch_artist=True);
#        for patch, color in zip(box['boxes'], colorspep0[:peptide_intensity.shape[0]]):    patch.set_facecolor(color)
#        a01.set_ylabel('log intensity')
#  
#
#        a1.plot(xpos,-np.log10(genepqtlpv2),'ko',markersize=1.5)
#
#        for ippp,ppp in enumerate(pepqtlpv):
#
#            a3.plot(xpos,-np.log10(ppp),color=colorspep0[ippp],marker='d',markersize=2,lw=0,label=peptides[ippp].split('_')[2])
#            a3.plot(xpos[-np.log10(ppp)>np.percentile(-np.log10(ppp),99)],-np.log10(ppp)[-np.log10(ppp)>np.percentile(-np.log10(ppp),99)],colorspep0[ippp],markersize=5,lw=0,marker='d')
#        a3.legend(loc=1,fontsize=5,numpoints=1)
#        mg=sb.heatmap(scst.spearmanr(peptide_intensity.T)[0][::-1],ax=a03)
#        a03.set_xticklabels(np.array([pppp.split('_')[2] for pppp in peptides]),rotation=90,fontsize=5)
#        a03.set_ylabel('peptide correlation')
#
#    a2.plot(xpos,-np.log10(geneeqtlpv),'ko',markersize=1.5)
#    if (geneeqtlqv<0.01).sum()>0:
#        a2.plot((np.min(xpos),np.max(xpos)),(-np.log10(np.max(geneeqtlpv[geneeqtlqv<0.01])),-np.log10(np.max(geneeqtlpv[geneeqtlqv<0.01]))),'k--',lw=0.5)
#        a2.annotate("FDR 0.01", xy=(np.min(xpos)+0.05,-np.log10(np.max(geneeqtlpv[geneeqtlqv<0.01]))+0.05 ),fontsize=10)
#    a2.set_ylabel("eQTL\n -log10PV ",fontsize=12,rotation=0,horizontalalignment= 'right')
#    a2.plot(xpos[-np.log10(genepqtlpv)>np.percentile(-np.log10(genepqtlpv),99.5)],-np.log10(geneeqtlpv)[-np.log10(genepqtlpv)>np.percentile(-np.log10(genepqtlpv),99.5)],'ro',markersize=5)
#    a2.plot(xpos[-np.log10(geneeqtlpv)>np.percentile(-np.log10(geneeqtlpv),99.5)],-np.log10(geneeqtlpv)[-np.log10(geneeqtlpv)>np.percentile(-np.log10(geneeqtlpv),99.5)],'g*',markersize=8)
#
#
#
#    plt.savefig(folder_destination+plot_name+gene+'.png',dpi=600)
#    
#==============================================================================
#==============================================================================
# # 
#==============================================================================
#==============================================================================

#
#print (list(file.keys()))
#file.close()

#def fun95(x):return np.nanpercentile(x,95)
#def fun99(x):return np.nanpercentile(x,99)
#def study_relationship_batery(x,y,xlab='',ylab=''):
#    percentiles=np.nanpercentile(x,np.arange(10)*10)
#    perc=np.array([j*(x>=i) for j,i in enumerate(percentiles)]).max(axis=0)
#    rrp=np.squeeze(npi.group_by(perc, y,fun95)).T
##    plt.figure(figsize=(5,5))
#    plt.plot(x,y,'o',markersize=3,label='corS = '+str(np.around(scst.spearmanr(x,y)[0],2)),color='grey')
#    rrp=np.squeeze(npi.group_by(perc, y,fun95)).T
#    plt.plot(percentiles[[int(i) for i in rrp[0]]],rrp[1],'ro-',markersize=5,label='95 perc')
#    rrp=np.squeeze(npi.group_by(perc, y,fun99)).T
#    plt.plot(percentiles[[int(i) for i in rrp[0]]],rrp[1],'bo-',markersize=5,label='99 perc')
#    plt.xlabel(xlab)
#    plt.ylabel(ylab)
##    plt.title(np.corcoef)
#    plt.legend(loc=4,fontsize=8)
#    return [percentiles,rrp]
#def plot_eqtl_pqtl(comment,eqtlpv2,pqtlpv2,eqtlqv2,pqtlqv2,ebeta2,pbeta2,ann2,genome2,genomed2,index_ips_specific2,de=None):
#
#    index=np.arange(len(eqtlpv2))
#
#    plt.figure(figsize=(12, 12))
#
#    plt.subplot(3,2,2);
#    ee=eqtlpv2[index];pp=pqtlpv2[index];
#    plt.plot(-np.log10(ee),-np.log10(pp),'ko',markersize=4)
#    plt.plot(-np.log10(ee[index_ips_specific2]),-np.log10(pp[index_ips_specific2]),'go',markersize=4)
##    plt.plot(-np.log10(ee[index_ips_specific2]),-np.log10(pp[index_ips_specific2]),'ro',markersize=4)
#    plt.title('-log10 pvalue: pQTL vs eQTL; any >2')
#    plt.xlabel('eQTL')
#    plt.ylabel('pQTL')
#
#    plt.subplot(3,2,2); ee=eqtlqv2[index];pp=pqtlqv2[index];
#    plt.plot(-np.log10(ee),-np.log10(pp),'o',markersize=4)
#    plt.title('-log10 qvalue: pQTL vs eQTL ')
#    plt.xlabel('eQTL')
#    plt.ylabel('pQTL')
#
#    plt.subplot(3,2,3); ee=ebeta2[index];pp=pbeta2[index];
#    plt.plot((ee),(pp),'o',markersize=4)
#    plt.title('Beta: pQTL vs eQTL ')
#    plt.xlabel('eQTL')
#    plt.ylabel('pQTL')
#
#    plt.subplot(3,2,5); ee=np.array([float(a) for a in ann2[:,2]])[index];pp=pqtlpv2[index];
#    plt.plot((ee),-np.log10(pp),'o',markersize=4)
#    plt.title('# cisSNP vs pQTL_pv (log10) ')
#    plt.xlabel('# p_cisSNP')
#    plt.ylabel('pQTL')
#
#    plt.subplot(3,2,6); ee=np.array([float(a) for a in ann2[:,3]])[index];pp=eqtlpv2[index];
#    plt.plot((ee),-np.log10(pp),'o',markersize=4)
#    plt.title('# cisSNP vs eQTL_pv (log10) ')
#    plt.xlabel('# e_cisSNP')
#    plt.ylabel('eQTL')
#    plt.tight_layout()#
#


#for gene in ['CRYZ']:#    if de1[ann1[:6]==gene]>0.5:
#    plot_manhatan_alone(comment='_black_egenes_',gene=gene,genome2=genomed1, ann2=ann1,    yp=np.squeeze(ypd[genespepr==gene]),\
#    yr=np.squeeze(yrd[genespepr==gene]),edata=eqtlpv1,pdata=pqtlpv1,edata2=eqtlqv1,pdata2=pqtlqv1)

#summary_gene_feature(chr_list=np.arange(1,3), folder_name='/Users/mirauta/Scripts/QTL_pipeline_plots/', feature_name='Exon_1234',local_adjustment_method='')

#fOut=h5py.File(folder_data+'/'+feature_report+'_'+output_file+'.h5','r')
#
#eqtl=h5py.File('/Users/mirauta/Results/hipsci/QTL/eQTLcis.results.full.ipsc.chr21.hdf5','r')
#
#genes=np.intersect1d(np.array(list(fOut.keys())),np.array(list(eqtl.keys())))
#
#epv=np.array([np.nanmin(eqtl[g]['lmm_peer_k30_covar_residuals_gauss_pv'][:])*\
#              eqtl[g]['lmm_peer_k30_covar_residuals_gauss_pv'].shape[1] for g in genes])
#
#ppv=np.array([fOut[g]['summary_data']['min_p_value_local_adjusted'][:] for g in genes])
#
#plt.plot(-np.log10(epv),-np.log10(ppv),'o')
#scst.spearmanr(-np.log10(epv),-np.log10(ppv))

#    print(list(eqtl[g].keys()))
#    print()
##(eqtl[list(eqtl.keys())[9]]['gdid'])