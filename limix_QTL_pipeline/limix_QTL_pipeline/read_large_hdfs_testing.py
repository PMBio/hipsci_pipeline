import tables
import pandas as pd
import numpy as np
import time


h5file = tables.open_file("../data/qtl_dummy_test_data_20m_rows.h5",'r')

table = h5file.root.qtl_analysis.results

feature_ids = list(set([row['feature_id'] for row in table]))

record_count = 0

t0 = time.time()

for fid in feature_ids[:50000]:
    snp_ids = [row['snp_id'] for row in table.where(
            '''feature_id==fid''')]
    record_count += len(snp_ids)

t1 = time.time()-t0
    
h5file.close()