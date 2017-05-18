import tables

class hdf5_writer:

    def __init__(self,output_filename):
        self.h5file = tables.open_file(output_filename,'w')
        self.group = self.h5file.create_group('/','qtl_results')

    def close(self):
        self.h5file.close()

    def add_result_df(self,qtl_results_df):
        assert(len(set(qtl_results_df['feature_id']==1)))
        feature_id = qtl_results_df['feature_id'].values[0]
        column_names = ['snp_id','p_value','beta','n_samples']
        group = self.h5file.create_group('/qtl_results/',feature_id)
        table = self.h5file.create_table(group,'results',QTL_result_hdf5,"QTL analysis results")
        qtl_result = table.row
        for idx,df_row in qtl_results_df.iterrows():
            for col_name in column_names:
                qtl_result[col_name] = df_row[col_name]
            qtl_result.append()
        table.flush()
        table.cols.snp_id.create_index()
 
class text_writer:

    def __init__(self,output_filename):
        self.column_names = ['feature_id','snp_id','p_value','beta','n_samples']
        with open(output_filename,'w') as f:
            header = '\t'.join(self.column_names)
            f.write(header+'\n')
        self.outfile = open(output_filename,'a')

    def close(self):
        self.outfile.close()

    def add_result_df(self,qtl_results_df):
        qtl_results_df.loc[:,self.column_names].to_csv(self.outfile,header=None,mode='a',index=False,sep='\t')

       
        
class QTL_result_hdf5(tables.IsDescription):
    snp_id  = tables.StringCol(16)   # 16-character String
    p_value = tables.Float64Col()    # double (double-precision)
    beta = tables.Float64Col()    # double (double-precision)
    n_samples = tables.Int32Col()    # integer
