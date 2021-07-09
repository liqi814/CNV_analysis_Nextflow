import csv
import pandas as pd

sample_file_name = 'data/IPF_short_telo_noqv.csv'

sample_internal_name_list = []
with open(sample_file_name, 'r') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter='\n')
    for row in csv_reader:
        sample_internal_name_list.extend(row)
sample_internal_name_list = [x.strip() for x in sample_internal_name_list]
sample_internal_name_list_string = "('" + "','".join(sample_internal_name_list) + "')"

sample_already_done = ['IPF0573', 'IPF1916', 'IPF3411', 'IPF-ACE041', 'IPF-ACE088', 'IPF-4310355', 'IPF3353', 'IPF3579',
                       'IPF3578', 'IPF3611', 'IPF0603', ' IPF1773']
sample_internal_name_list = [x for x in sample_internal_name_list if x not in sample_already_done]
# 276 left

with open('data/IPF_short_telo_noqv_finalList.csv', 'w') as csvfile:
    writer = csv.writer(csvfile, delimiter=',')
    writer.writerow(sample_internal_name_list)

# sample_result_16940136_1623356605.csv

sample_result = pd.read_csv("sample_result_16940136_1623356605.csv")
sample_result = sample_result[sample_result.SubProject.str.contains("Garcia")]

sample_result['rawBamFileLoc'] = sample_result.AlignSeqFileLoc + '/' + sample_result.sample_internal_name + '.' + sample_result.experimentID.astype(str) + '/' + sample_result.sample_internal_name + '.' + sample_result.experimentID.astype(str) + '.bam'

sample_result[sample_result.SelfDeclGender == 'F'].rawBamFileLoc.to_csv("WGS_BAM_File/Female_bam_file_path.list", sep="\t", header=False, index=False)
sample_result[sample_result.SelfDeclGender == 'M'].rawBamFileLoc.to_csv("WGS_BAM_File/Male_bam_file_path.list", sep="\t", header=False, index=False)