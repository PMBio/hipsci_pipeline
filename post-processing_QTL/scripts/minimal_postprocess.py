import h5py
import glob
import os.path
import numpy as np
import pandas as pd

def minimal_qtl_processing(QTL_Dir, OutputDir, writeToOneFile=True, compressed = False, overWrite=True, minimalPValue = 1, minimalFeaturePValue = 1):
    qtl_results_file='qtl_results_'
    snp_metadata_file='snp_metadata_'
    feature_metadata_file='feature_metadata_'
    output_file='qtl_results_'


    h5FilesToProcess = (glob.glob(QTL_Dir+"/qtl_*.h5"))
    #print(h5FilesToProcess)
    #print(os.path.dirname(h5FilesToProcess[1]))
    for file in h5FilesToProcess :
        partTmp = os.path.basename(file).replace(qtl_results_file,"").replace(".h5","")
        
        if(writeToOneFile): 
            outputFile = OutputDir+output_file+"all.txt"
        else:
            outputFile = OutputDir+output_file+partTmp+".txt"
        
        #print(outputFile)
        if(((os.path.isfile(outputFile) or os.path.isfile(outputFile+".gz")) and not overWrite) and not writeToOneFile):
            #print("Skipping: "+partTmp)
            continue
        #else :
            #print('Processing: '+partTmp)
        #print(partTmp)
        if not os.path.isfile(QTL_Dir+"/"+snp_metadata_file+partTmp+".txt"):
            print("Skipping: " +partTmp + " not all necessary files are present.")
            continue
        if not os.path.isfile(QTL_Dir+"/"+feature_metadata_file+partTmp+".txt"):
            print("Skipping: " +partTmp + " not all necessary files are present.")
            continue
        try :
            #print(QTL_Dir+"/"+feature_metadata_file+partTmp+".txt")
            #print(QTL_Dir+"/"+snp_metadata_file+partTmp+".txt")
            ffea= pd.read_table(QTL_Dir+"/"+feature_metadata_file+partTmp+".txt", sep='\t')
            fsnp= pd.read_table(QTL_Dir+"/"+snp_metadata_file+partTmp+".txt", sep='\t')
        except:
            print("Issue in features or snp annotation.\n Skipping: "+partTmp)
            continue

        ffea = ffea.rename(index=str, columns={"chromosome": "feature_chromosome", "start": "feature_start", "end": "feature_end"})
        fsnp = fsnp.rename(index=str, columns={"chromosome": "snp_chromosome", "position": "snp_position"})

        frez=h5py.File(file,'r')
        frezkeys= np.array([k.replace('_i_','') for k in list(frez.keys())])

        data={}
        for key in ['feature_id','snp_id','p_value','beta','beta_se','empirical_feature_p_value']:
            data[key]=np.zeros(len(np.unique(list(frezkeys))),dtype='object')+np.nan

        for ifea,report_feature in enumerate(np.unique(list(frezkeys))):
            for key in ['snp_id','p_value','beta','beta_se','empirical_feature_p_value']:
                temp = np.array(frez[report_feature][key])
                data[key][ifea]=np.hstack(temp).astype('U')
            data['feature_id'][ifea]=np.hstack(np.repeat(report_feature,len(frez[report_feature][key])))

        for key in data.keys():
            data[key]=np.hstack(data[key])

        temp=pd.DataFrame(data)
        #print(temp.head())

        temp = pd.merge(temp, ffea, on='feature_id')
        #print(temp.head())
        temp = pd.merge(temp, fsnp, on='snp_id')
        #print(temp.head())
        temp['empirical_feature_p_value'] = temp['empirical_feature_p_value'].astype(float)
        temp['p_value'] = temp['p_value'].astype(float)
        if(minimalPValue<1):
            temp=pd.DataFrame(data).iloc[data['p_value'].astype(float)<minimalPValue]
        
        if(minimalFeaturePValue<1):
            temp=pd.DataFrame(data).iloc[data['empirical_feature_p_value'].astype(float)<minimalFeaturePValue]
        
        temp = temp.sort_values(by =['empirical_feature_p_value',"p_value"], ascending=[True,True])
        
        #print(outputFile)#
        #temp.to_csv(path_or_buf=outputFile, mode='w', sep='\t', columns=None,index=None)
        #print('w'if not os.path.isfile(outputFile) else 'a')
        if(not compressed):
            temp.to_csv(path_or_buf=outputFile, mode='w'if not os.path.isfile(outputFile) else 'a', sep='\t', columns=None,index=None)    
        else:
            temp.to_csv(path_or_buf=outputFile+".gz", mode='w'if not os.path.isfile(outputFile) else 'a', sep='\t', columns=None,index=None,compression='gzip')

def minimal_qtl_processing_top_snp(QTL_Dir, OutputDir, writeToOneFile=True, compressed = False, overWrite=True, minimalPValue = 1, minimalFeaturePValue = 1):
    qtl_results_file='qtl_results_'
    snp_metadata_file='snp_metadata_'
    feature_metadata_file='feature_metadata_'
    output_file='top_qtl_results_'


    h5FilesToProcess = (glob.glob(QTL_Dir+"/qtl_*.h5"))
    #print(h5FilesToProcess)
    #print(os.path.dirname(h5FilesToProcess[1]))
    for file in h5FilesToProcess :
        partTmp = os.path.basename(file).replace(qtl_results_file,"").replace(".h5","")
        #print(partTmp)
        if(writeToOneFile): 
            outputFile = OutputDir+output_file+"all.txt"
        else:
            outputFile = OutputDir+output_file+partTmp+".txt"
        #print(outputFile)
        if(((os.path.isfile(outputFile) or os.path.isfile(outputFile+".gz")) and not overWrite) and not writeToOneFile):
            #print("Skipping: "+partTmp)
            continue
        #else :
            #print('Processing: '+partTmp)
        #print(partTmp)
        if not os.path.isfile(QTL_Dir+"/"+snp_metadata_file+partTmp+".txt"):
            print("Skipping: " +partTmp + " not all necessary files are present.")
            continue
        if not os.path.isfile(QTL_Dir+"/"+feature_metadata_file+partTmp+".txt"):
            print("Skipping: " +partTmp + " not all necessary files are present.")
            continue
        try :
            #print(QTL_Dir+"/"+feature_metadata_file+partTmp+".txt")
            #print(QTL_Dir+"/"+snp_metadata_file+partTmp+".txt")
            ffea= pd.read_table(QTL_Dir+"/"+feature_metadata_file+partTmp+".txt", sep='\t')
            fsnp= pd.read_table(QTL_Dir+"/"+snp_metadata_file+partTmp+".txt", sep='\t')
        except:
            print("Issue in features or snp annotation.\n Skipping: "+partTmp)
            continue

        ffea = ffea.rename(index=str, columns={"chromosome": "feature_chromosome", "start": "feature_start", "end": "feature_end"})
        fsnp = fsnp.rename(index=str, columns={"chromosome": "snp_chromosome", "position": "snp_position"})

        frez=h5py.File(file,'r')
        frezkeys= np.array([k.replace('_i_','') for k in list(frez.keys())])

        data={}
        for key in ['feature_id','snp_id','p_value','beta','beta_se','empirical_feature_p_value']:
            data[key]=np.zeros(len(np.unique(list(frezkeys))),dtype='object')+np.nan

        for ifea,report_feature in enumerate(np.unique(list(frezkeys))):
            for key in ['snp_id','p_value','beta','beta_se','empirical_feature_p_value']:
                temp = np.array(frez[report_feature][key])
                data[key][ifea]=np.hstack(temp).astype('U')
            data['feature_id'][ifea]=np.hstack(np.repeat(report_feature,len(frez[report_feature][key])))

        for key in data.keys():
            data[key]=np.hstack(data[key])

        temp=pd.DataFrame(data)
        #print(temp.head())

        temp = pd.merge(temp, ffea, on='feature_id')
        #print(temp.head())
        temp = pd.merge(temp, fsnp, on='snp_id')
        #print(temp.head())
        temp['empirical_feature_p_value'] = temp['empirical_feature_p_value'].astype(float)
        temp['p_value'] = temp['p_value'].astype(float)
        if(minimalPValue<1):
            temp=pd.DataFrame(data).iloc[data['p_value'].astype(float)<minimalPValue]
        
        if(minimalFeaturePValue<1):
            temp=pd.DataFrame(data).iloc[data['empirical_feature_p_value'].astype(float)<minimalFeaturePValue]
        
        temp = temp.sort_values(by =['empirical_feature_p_value',"p_value"], ascending=[True,True])
        #print(temp)
        temp = temp.groupby(temp['feature_id']).first()
        #print(temp)#
        #temp.to_csv(path_or_buf=outputFile, mode='w', sep='\t', columns=None,index=None)
        #print('w'if not os.path.isfile(outputFile) else 'a')
        if(not compressed):
            temp.to_csv(path_or_buf=outputFile, mode='w'if not os.path.isfile(outputFile) else 'a', sep='\t', columns=None , header = True if not os.path.isfile(outputFile) else False)
        else:
            temp.to_csv(path_or_buf=outputFile+".gz", mode='w'if not os.path.isfile(outputFile) else 'a', sep='\t', columns=None,compression='gzip', header = True if not os.path.isfile(outputFile) else False)