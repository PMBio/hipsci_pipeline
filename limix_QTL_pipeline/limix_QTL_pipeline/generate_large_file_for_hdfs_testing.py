import tables
import pandas as pd
import numpy as np
import random

random.seed(0)

class QTL_result(tables.IsDescription):
    feature_id = tables.StringCol(16)   # 16-character String
    snp_id  = tables.StringCol(16)   # 16-character String
    p_value = tables.Float64Col()    # double (double-precision)
    beta = tables.Float64Col()    # double (double-precision)


datafile = '../data/geuvadis_CEU_YRI_test_data/limix_QTL_results.txt'
qtl_results_df = pd.read_csv(datafile,sep='\t')


string_randomiser = lambda s:''.join(random.sample(s,len(s)))



h5file = tables.open_file("../data/qtl_dummy_test_data_20m_rows.h5",'w')
group = h5file.create_group('/','qtl_analysis')
table = h5file.create_table(group,'results',QTL_result,"QTL analysis results")
qtl_result = table.row

column_names = ['feature_id','snp_id','p_value','beta']

#for 1m rows
nMultiplier = 50
#for 1bn rows
nMultiplier = 50*1000
##for 10bn rows
#nMultiplier = 50*10000

for idx in range(nMultiplier):
    feature_ids = list(set(qtl_results_df['feature_id']))
    shuffled_id_dict = dict([(x,string_randomiser(x)) for x in feature_ids])
    dummy_qtl_results_df = qtl_results_df
    dummy_qtl_results_df['feature_id'] = qtl_results_df['feature_id'].map(lambda x: shuffled_id_dict[x])
    for idx,df_row in dummy_qtl_results_df.iterrows():
        for col_name in column_names:
            qtl_result[col_name] = df_row[col_name]
        qtl_result.append()
    table.flush()

indexrows = table.cols.feature_id.create_index()

h5file.close()