import pandas as pd
column_mapping={'feature_id':'Ensembl Gene ID','gene_id':'Ensembl Gene ID','gene_name':'Associated Gene Name','chromosome':'Chromosome Name','start':'Gene Start (bp)','end':'Gene End (bp)'}

geuvadis_annotation_df = pd.read_csv('Geuvadis_CEU_YRI_Annot.txt',sep='\t')

annotation_df = pd.DataFrame(index=geuvadis_annotation_df.index,columns=['feature_id','gene_id','gene_name','chromosome','start','end'])

for col in column_mapping.keys():
    annotation_df[col] = geuvadis_annotation_df[column_mapping[col]]

annotation_df['strand'] = '.'

annotation_df.to_csv('Geuvadis_CEU_YRI_formatted_annotation_data.txt',sep='\t',index=False)
