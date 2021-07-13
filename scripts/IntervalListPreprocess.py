import os
import glob
import datetime
import pandas as pd

os.chdir("D:\\CNV_analysis\\2021-7-13\\IntervalList")


with open("final_centromere_hg19.seg", 'r') as file:
    centromere = []
    for line in file:
        if line.startswith("@") :
            continue
        else:
            centromere.append(line.strip().split(sep="\t"))

centromere_df = pd.DataFrame(centromere[1:], columns=centromere[0])

ENCODE = pd.read_table("hg19-blacklist.v2.bed", sep="\t", header=None)
ENCODE[0] = ENCODE[0].str.lstrip("chr")
ENCODE.columns = ["CONTIG","START","END", "Annotation"]
PAR = pd.read_table("PAR.bed", header=None)
PAR = PAR[0].str.split("   ", expand=True)
PAR.columns = ["CONTIG","START","END"]

centromere_df[['START', 'END']] = centromere_df[['START', 'END']].astype('int32')
ENCODE.loc[:,1:2] = ENCODE.loc[:,1:2].astype('int32')
PAR.loc[:,1:2] =  PAR.loc[:,1:2].astype('int32')

text_file = open('hg19-blacklist_centromere_ENCODE_PAR.bed', "wt")

for chr in [str(x) for x in list(range(1,23))] + ['X', 'Y']:
    new_df = pd.concat([ENCODE[ENCODE.CONTIG == chr][["CONTIG","START","END"]], PAR[PAR.CONTIG == chr][["CONTIG","START","END"]],
                        centromere_df[centromere_df.CONTIG == chr][["CONTIG","START","END"]]]).sort_values(by=['START'], ignore_index=True)
    last_index = 0
    for index in range(1,len(new_df)) :
        if new_df.iloc[index, 1] <= new_df.iloc[last_index, 2]:
            if new_df.iloc[index, 2] > new_df.iloc[last_index, 2]:
                new_df.iloc[last_index, 2] = new_df.iloc[index, 2]
        else:
            text_file.write(str(new_df.iloc[last_index,0]) + ":" + str(new_df.iloc[last_index,1]) + "-" + str(new_df.iloc[last_index,2]) + "\n")
            last_index = index
    text_file.write(str(new_df.iloc[last_index, 0]) + ":" + str(new_df.iloc[last_index, 1]) + "-" + str(
        new_df.iloc[last_index, 2]) + "\n")

text_file.close()