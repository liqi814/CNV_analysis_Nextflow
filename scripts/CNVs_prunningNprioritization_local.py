import os
import glob
import datetime
import pandas as pd
import numpy as np

os.chdir("D:\\CNV_analysis\\2021-6-28")
output_file_prefix = str(datetime.datetime.now())[:10] + "_sample_called_seg"

calledSeg_List = glob.glob('samples/*/CopyRatioSegments/*')
for calledSegfile in calledSeg_List:
    called_seg = []
    with open(calledSegfile, 'r') as file:
        # called_seg = file.readlines()
        for line in file:
            if line.startswith("@RG") :
                sample_name = line.strip().split(sep="\t")[2].lstrip('SM:')
            elif line.startswith("@") :
                continue
            else:
                called_seg.append(line.strip().split(sep="\t"))

    called_seg = pd.DataFrame(called_seg)
    called_seg.rename(columns=called_seg.iloc[0], inplace=True)
    called_seg = called_seg[1:]
    called_seg = called_seg[called_seg['CALL'] != '0']
    called_seg['SAMPLE_NAME'] = sample_name
    called_seg[['START', 'END', 'NUM_POINTS_COPY_RATIO', 'MEAN_LOG2_COPY_RATIO']] = called_seg[
        ['START', 'END', 'NUM_POINTS_COPY_RATIO', 'MEAN_LOG2_COPY_RATIO']].apply(pd.to_numeric)
    called_seg['SV_type'] = called_seg.apply(lambda row: 'DEL' if row.CALL == '-' else 'DUP', axis=1)

    if not os.path.exists(output_file_prefix + '.bed'):
        called_seg[['CONTIG', 'START', 'END', 'SV_type', 'NUM_POINTS_COPY_RATIO', 'MEAN_LOG2_COPY_RATIO', 'CALL',
                    'SAMPLE_NAME']] \
            .to_csv(output_file_prefix + '.bed', sep='\t', header=True, index=False)
    else:
        called_seg[
            ['CONTIG', 'START', 'END', 'SV_type', 'NUM_POINTS_COPY_RATIO', 'MEAN_LOG2_COPY_RATIO', 'CALL',
             'SAMPLE_NAME']] \
            .to_csv(output_file_prefix + '.bed', sep='\t', mode='a', header=False, index=False)


# reading annotated file and RLCR file
# 7ITrJqgNMY
RLCRs_no_Repeat_Masker_file = "C:\\Users\\ql2387\\Desktop\\CNV_analysis_Nextflow\\data\\\RLCRs\\RLCRs_no_Repeat_Masker.txt"
RLCRs_no_Repeat_Masker = pd.read_csv(RLCRs_no_Repeat_Masker_file, sep="\t", names=['CONTIG', 'REGION_START', 'REGION_END', 'REGION'])
RLCRs_no_Repeat_Masker['CONTIG'] = RLCRs_no_Repeat_Masker['CONTIG'].map(lambda x: x.lstrip('chr'))
# RLCRs_no_Repeat_Masker.CONTIG.str.lstrip("chr")
RLCRs_no_Repeat_Masker[['REGION_START', 'REGION_END']] = RLCRs_no_Repeat_Masker[['REGION_START', 'REGION_END']].apply(pd.to_numeric)
annotated_file = pd.read_csv("AnnotSV_7ITrJqgNMY.tsv", sep="\t", low_memory=False)
annotated_file.rename({'user#1': 'NUM_POINTS_COPY_RATIO', 'user#2': 'MEAN_LOG2_COPY_RATIO', 'user#3': 'CALL', 'user#4': 'SAMPLE_NAME'}, axis=1, inplace=True)

for index, segment in annotated_file.iterrows():
    if segment.Annotation_mode == "full":
        startPOS = segment.SV_start
        endPOS = segment.SV_end
    else:
        startPOS = segment.Intersect_start
        endPOS = segment.Intersect_end

    overlapped_segments = RLCRs_no_Repeat_Masker.loc[(RLCRs_no_Repeat_Masker.CONTIG == segment.SV_chrom) &
                                     (((RLCRs_no_Repeat_Masker.REGION_START <= startPOS) & (startPOS <= RLCRs_no_Repeat_Masker.REGION_END)) |
                                     ((RLCRs_no_Repeat_Masker.REGION_START <= endPOS) & (endPOS <= RLCRs_no_Repeat_Masker.REGION_END)) |
                                     ((RLCRs_no_Repeat_Masker.REGION_START >= startPOS) & (endPOS >= RLCRs_no_Repeat_Masker.REGION_END)))]
    if len(overlapped_segments.REGION)  == 0 :
        annotated_file.loc[index, 'overlapped_RLCR'] = None
    else:
        arr = sorted(list(overlapped_segments.REGION_START) + list(overlapped_segments.REGION_END) + [startPOS, endPOS])
        if (arr[-2] - arr[1] - sum([a -b for a, b in zip(overlapped_segments.REGION_START[1:],overlapped_segments.REGION_END[:-1])])) / (endPOS - startPOS) >= 0.7 :
            annotated_file.loc[index, 'overlapped_RLCR'] = overlapped_segments.REGION.str.cat(sep=';')
        else:
            annotated_file.loc[index, 'overlapped_RLCR'] = None

annotated_file.to_excel(output_file_prefix + '_annotated_RLCR.xlsx', header=True, index=False)