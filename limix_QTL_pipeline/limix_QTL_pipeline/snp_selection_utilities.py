import pandas as pd

def get_snps(feature_id, annotation_df, bim, cis_mode, window_size):
    annotation_sub_df = annotation_df.loc[[feature_id],:]
    if len(annotation_sub_df['chromosome'].unique())>1:
        print('Warning! Multiple regions specified across multiple chromosomes.\n'
              'If running chromosome-by-chromosome, each feature will be tested \n'
              ' against all relevant SNPs for each chromosome.')
    list_of_snp_dfs = []
    for annotation_ds in annotation_sub_df.itterrows():
        chrom = str(annotation_ds.loc['chromosome'])
        start = annotation_ds.loc['start']
        end = annotation_ds.loc['end']
        # make robust to features selfpecified back-to-front
        lowest = min([start,end])
        highest = max([start,end])
        if (cis_mode) :
            # for cis, we sequentially add snps that fall within each region
            snpQuery = bim.query("chrom == '%s' & pos > %d & pos < %d" % (chrom, lowest-window_size, highest+window_size))
            list_of_snp_dfs.append(snpQuery)
        else :
            # for trans, we sequentially exclude snps that fall within each region
            snpQuery = snpQuery.query("(chrom == '%s' & (pos < %d | pos > %d))|chrom != '%s'" % (chrom, lowest-window_size, highest+window_size,chrom))
    if (cis_mode):
        selected_snp_df = pd.concat(list_of_snp_dfs).drop_duplicates()
    else:
        selected_snp_df = snpQuery
    # filtering for sites on non allosomes.
    selected_snp_df = selected_snp_df.loc[selected_snp_df['chrom'].map(lambda x: x in list(map(str, range(1, 23))))]
    
    return selected_snp_df
