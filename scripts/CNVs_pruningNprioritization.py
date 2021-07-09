#!/usr/bin/python
# usage: test.py -i <inputfile> -o <outputfile>

import pandas as pd
import numpy as np
import os
import sys
import datetime

input_file = sys.argv[1]
sample_name = sys.argv[2]

# input_file = "C:\\Users\\ql2387\\Desktop\\2021-05_TERTp_CNV\\results\\IPF0284.194193\\CopyRatioSegments\\IPF0284.194193.1000.called.seg"
output_file_prefix = "/nfs/external/az-ipf-garcia/CNVanalysis/" + str(datetime.datetime.now())[:10] + "_sample_called_seg"
# sample_name = "IPF0284.194193"
RLCRs_no_Repeat_Masker_file = "/nfs/external/az-ipf-garcia/reference/RLCRs_no_Repeat_Masker.txt"

RLCRs_no_Repeat_Masker = pd.read_csv(RLCRs_no_Repeat_Masker_file, sep="\t", names=['CONTIG', 'REGION_START', 'REGION_END', 'REGION'])
RLCRs_no_Repeat_Masker['CONTIG'] = RLCRs_no_Repeat_Masker['CONTIG'].map(lambda x: x.lstrip('chr'))
# RLCRs_no_Repeat_Masker.CONTIG.str.lstrip("chr")
RLCRs_no_Repeat_Masker[['REGION_START', 'REGION_END']] = RLCRs_no_Repeat_Masker[['REGION_START', 'REGION_END']].apply(pd.to_numeric)

called_seg = []
with open (input_file, 'r') as file:
 # called_seg = file.readlines()
    for line in file:
        if line.startswith("@"):
            continue
        else:
            called_seg.append(line.strip().split(sep="\t"))

called_seg = pd.DataFrame(called_seg)
called_seg.rename(columns=called_seg.iloc[0], inplace = True)
called_seg = called_seg[1:]
called_seg = called_seg[called_seg['CALL'] != '0']
called_seg['SAMPLE_NAME'] = sample_name
called_seg[['START', 'END', 'NUM_POINTS_COPY_RATIO', 'MEAN_LOG2_COPY_RATIO']] = called_seg[[ 'START', 'END', 'NUM_POINTS_COPY_RATIO', 'MEAN_LOG2_COPY_RATIO']].apply(pd.to_numeric)

for index, segment in called_seg.iterrows():
    overlapped_segments = RLCRs_no_Repeat_Masker.loc[(RLCRs_no_Repeat_Masker.CONTIG == segment.CONTIG) &
                                     (((RLCRs_no_Repeat_Masker.REGION_START <= segment.START) & (segment.START <= RLCRs_no_Repeat_Masker.REGION_END)) |
                                     ((RLCRs_no_Repeat_Masker.REGION_START <= segment.END) & (segment.END <= RLCRs_no_Repeat_Masker.REGION_END)) |
                                     ((RLCRs_no_Repeat_Masker.REGION_START >= segment.START) & (segment.END >= RLCRs_no_Repeat_Masker.REGION_END)))]
    if len(overlapped_segments.REGION)  == 0 :
        called_seg.loc[index, 'overlapped_RLCR'] = None
    elif len(overlapped_segments.REGION)  == 1 :
        arr = [int(overlapped_segments.REGION_START), int(overlapped_segments.REGION_END), segment.START, segment.END]
        if (np.percentile(arr, 75) - np.percentile(arr, 25)) / (segment.END -segment.START) >= 0.7 :
            called_seg.loc[index, 'overlapped_RLCR'] = overlapped_segments.REGION.str.cat(sep=';')
        else:
            called_seg.loc[index, 'overlapped_RLCR'] = None
    else:
        called_seg.loc[index, 'overlapped_RLCR'] = overlapped_segments.REGION.str.cat(sep=';')

if not os.path.exists(output_file_prefix + '_wRLCRinfo.csv'):
    called_seg.to_csv(output_file_prefix + '_wRLCRinfo.csv', header=True, index=False)
else:
    called_seg.to_csv(output_file_prefix + '_wRLCRinfo.csv', mode='a', header=False, index=False)

called_seg['SV_type'] = called_seg.apply(lambda row : 'DEL' if row.CALL == '-' else 'DUP', axis=1)

if not os.path.exists(output_file_prefix + '.bed'):
    called_seg[['CONTIG', 'START', 'END', 'SV_type', 'NUM_POINTS_COPY_RATIO', 'MEAN_LOG2_COPY_RATIO', 'CALL', 'SAMPLE_NAME']]\
    .to_csv(output_file_prefix + '.bed', sep='\t', header=True, index=False)
else:
    called_seg[
        ['CONTIG', 'START', 'END', 'SV_type', 'NUM_POINTS_COPY_RATIO', 'MEAN_LOG2_COPY_RATIO', 'CALL', 'SAMPLE_NAME',
         'overlapped_RLCR']] \
        .to_csv(output_file_prefix + '.bed', sep='\t', mode='a', header=False, index=False)
