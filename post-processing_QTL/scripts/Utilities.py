
import csv as csv
import numpy as np
csv.field_size_limit(2147483647)


def read_parameters(parameters_file):
    with open(parameters_file, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter='=', quotechar='"')
        parameters = np.array([row for row in reader])
    print("read all parameters from parameters file")
    
    try:    folder_dest=parameters[np.array(['folder_dest' in p for p in parameters[:,0]]),1][0]
    except: print('no folder dest')
    try:    folder_data=parameters[np.array(['folder_data' in p for p in parameters[:,0]]),1][0]
    except: print('no folder data ')
    try:	 batch_name_tmt=np.array([parameters[i,1] for i,p in enumerate(parameters[:,0]) if 'batch_name_tmt' in p])
    except: print('no tmt data')
    try:    file_protein_quantification=parameters[np.array(['file_protein_quantification' in p for p in parameters[:,0]]),1][0]
    except: print('no tmt data')
    
    return [ folder_dest,folder_data,batch_name_tmt,file_protein_quantification]

#
#
#def write_data_csv(file_name,data,delimiter=','):
#    dim=len(data)
#    with open(file_name, 'w') as csvfile:
#        writer = csv.writer(csvfile, delimiter=delimiter, quotechar='"')
#        for n in range(data[0].shape[0]):
#            writer.writerow([data[c][n] for c in range (dim)])
#
# 
#
def load_tab_toh5(file_name,delimiter=',', quotechar='"',fields_to_load=None, skip_before=0,header=0,only_header=0,skip_after=0):
    with open(file_name, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=delimiter, quotechar=quotechar,quoting=csv.QUOTE_NONE)
        for ss in range(skip_before):
            next(reader, None)
        if header:
            headers = next(reader, None)
        else:
            headers=[]
        if only_header:
            if fields_to_load!=None:
                return [np.array(headers)[fields_to_load]]
            return [np.array(headers)]
        
        for ss in range(skip_after):
            next(reader, None)
        if fields_to_load is None:
            dat=[row for row in reader]
        else:		
            dat=[[row[fnr]  for fnr in  fields_to_load] for row in reader]
        
        dat=[np.array([dat[i][j] for i in np.arange(len(dat))]) for j in np.arange(len(dat[0]))]
    return [np.array(headers)[fields_to_load],dat]
#

def write_data(file_name,mat0=None,mat=None,header=None,delim='\t'):
    try: dim=mat.shape[1]
    except: print ('mat needs to have dimension 2')
    
    f = open(file_name,'w')
    if header is not None: 
        for c in header[:-1]:
            f.write(c+delim) # python will convert \n to os.linesep
        f.write(str(header[-1])+'\n')
    for n in range(mat.shape[0]):
        if mat0 is not None: 
            for c in mat0[n]:
                f.write(c+delim) # python will convert \n to os.linesep
        for c in mat[n][:-1]:
            f.write(str(c)+delim) 
        f.write(str(mat[n][-1])+'\n')
    f.close() #
    
def load_data(file_name,delimiter=',', quotechar='"',fields_to_load=None, skip_before=0,header=0,only_header=0,skip_after=0):
    with open(file_name, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=delimiter, quotechar=quotechar,quoting=csv.QUOTE_NONE)
        for ss in range(skip_before):
            next(reader, None)
        if header:
            headers = next(reader, None)
        else:
            headers=[]
        if only_header:
            if fields_to_load!=None:
                return [np.array(headers)[fields_to_load]]
            return [np.array(headers)]
        
        for ss in range(skip_after):
            next(reader, None)
        if fields_to_load is None:
            dat=[row for row in reader]
        else:		
            dat=[[row[fnr]  for fnr in  fields_to_load] for row in reader]
        
        dat=[np.array([dat[i][j] for i in np.arange(len(dat))]) for j in np.arange(len(dat[0]))]
    return [np.array(headers)[fields_to_load],dat]
#

        
def load_data_from_maxquant_output(folder_data='',  name='',data_file='.proteinGroups.txt', field_data='',float_data=1,\
                                   field_metadata=[]):
    
    data_path=folder_data+name+data_file
    
    allfields=np.array(load_data(file_name=data_path,delimiter='\t', header=1,only_header=1)[0])

    metadata={'fields': allfields[np.in1d(allfields,field_metadata)]}
    data={'fields': allfields[np.where([field_data in field for field in allfields])[0]]}

    alldata=load_data(file_name=data_path,fields_to_load=np.where(np.in1d(allfields,data['fields']))[0],delimiter='\t', header=1,       only_header=0)
    allmetadata=load_data(file_name=data_path,fields_to_load=np.where(np.in1d(allfields,metadata['fields']))[0],delimiter='\t',         header=1,only_header=0)

                                                                     #if float_data:
                                                                               #                  data['data']=np.squeeze(np.array(alldata[1])).astype(float)
                                                                                       #                             else:
                                                                                                #                      data['data']=np.squeeze(np.array(alldata[1])).astype('S3')
                                                                                                        # metadata['data']=np.squeeze(np.array(allmetadata[1])).astype('S60')
    if float_data: data[field_data]=np.squeeze(np.array(alldata[1])).astype(float)
    else:    data[field_data]=alldata[1]
    metadata['data']=allmetadata[1]
    rez={}
    rez['allfields']=allfields
    rez['meta']={}
    rez['data']={}
    for ikey,key in enumerate(metadata['fields']):
        rez['meta'][key]=metadata['data'][ikey]
    rez['data']['fields']=data['fields']
    rez['data'][field_data]=data[field_data]
    return rez


#
#
#def dumpDictHdf5(RV,f5):
#        """ Dump a dictionary where each page is a list or an array """
#        for key in RV.keys():
#                f5.create_dataset(name=key,data=np.array(RV[key]),chunks=True,compression='gzip')
#
#
#
#
#def rank_transform(x=None, ref=None):
#	_doc='rank transform x into ref;keep the range'
#        refn0=ref[ref>0]
#        xn0=x[x>0]
#        sref=np.sort(refn0)
#        index= np.linspace(0,len(refn0)-0.1, num=len(xn0)).astype(int);
#
#        orderx=np.argsort(xn0)
#
#        xtr=x*0
#        xtr[x>0]=sref[index][np.argsort(orderx)]
#        return xtr
#
#
#
#
#def normalize_reference(data, ntype='',referencedata=None):
#	if ntype=='no':
#		return data
#
#	elif ntype=='qnorm':
#        	return np.array([rank_transform(x=data[:,l],ref=referencedata) for l in range(data.shape[1])]).T
#
#	elif ntype=='snorm':
#		return np.array([data[:,l]/np.nansum(data[:,l])*np.nansum(referencedata) for l in range(data.shape[1])]).T
#
#
#def matchID(id1,id2):
#    """ match id1 and id2 """
#    idx1 = []
#    idx2 = []
#    for i in range(id1.shape[0]):
#        if id1[i] in id2.tolist():
#            idx1.append(i)
#            idx2.append(np.where(id2==id1[i])[0][0])
#    idx1 = np.array(idx1)
#    idx2 = np.array(idx2)
#    print (id1[idx1]==id2[idx2]).all()  # assert that's right
#    return idx1, idx2
#
#def matchID_order(id1,id2):
#    idx1=np.where(np.in1d(id1,id2))[0]
#    idx1=idx1[np.argsort(id1[idx1])]
#    return 1
#
