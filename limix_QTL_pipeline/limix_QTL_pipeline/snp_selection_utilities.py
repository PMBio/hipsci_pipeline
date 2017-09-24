import pandas as pd

def get_snps(feature_id, annotation_df, bim, cis_mode, window_size):
    annotation_sub_df = annotation_df.loc[[feature_id],:]
    list_of_snp_dfs = []
    for annotation_ds in annotation_sub_df.itterrows():
        chrom = str(annotation_ds.loc['chromosome'])
        start = annotation_ds.loc['start']
        end = annotation_ds.loc['end']
        # make robust to features selfpecified back-to-front
        lowest = min([start,end])
        highest = max([start,end])
        if (cis_mode) :
            snpQuery = bim.query("chrom == '%s' & pos > %d & pos < %d" % (chrom, lowest-window_size, highest+window_size))
        else :
            snpQuery = bim.query("(chrom == '%s' & (pos < %d | pos > %d))|chrom != '%s'" % (chrom, lowest-window_size, highest+window_size,chrom))
            #Filtering for sites on non allosomes.
            snpQuery = selfnpQuery.loc[snpQuery['chrom'].map(lambda x: x in list(map(str, range(1, 23))))]
        list_of_snp_dfs.append(snpQuery)
    selected_snp_df = pd.concat(list_of_snp_dfs).drop_duplicates()
    
    return selected_snp_df
