import tables

class hdf5_writer:

    def __init__(self,output_filename):
        self.h5file = tables.open_file(output_filename,'w')
        self.group = self.h5file.create_group('/','qtl_analysis')
        self.table = self.h5file.create_table(self.group,'results',QTL_result_hdf5,"QTL analysis results")

    def close(self):
        self.h5file.close()

    def add_result_df(self,qtl_results_df):
        column_names = ['feature_id','snp_id','p_value','beta']
        qtl_result = self.table.row
        for idx,df_row in qtl_results_df.iterrows():
            for col_name in column_names:
                qtl_result[col_name] = df_row[col_name]
            qtl_result.append()
        self.table.flush()
        
        
class QTL_result_hdf5(tables.IsDescription):
    feature_id = tables.StringCol(16)   # 16-character String
    snp_id  = tables.StringCol(16)   # 16-character String
    p_value = tables.Float64Col()    # double (double-precision)
    beta = tables.Float64Col()    # double (double-precision)
