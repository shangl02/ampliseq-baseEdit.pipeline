#Jiangchuan Ye 02/14/2023

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import regex
import re
import os
from Bio import SeqIO
import Bio
from Bio import motifs
from Bio import pairwise2 
from Bio.pairwise2 import format_alignment
import collections
import itertools

# set up working directory
os.chdir(r'X:\users\shangl02\Test\ampliCan.test2\AS15\fastq')
directory = os.getcwd()
filenames = os.listdir(directory)

# generate metadata dataframe from submitted sequence information sheet including ‘Sample #Number’, ‘Guide Sequence’ and ‘Full Amplicon’. Add two columns indicating reverse complement #sequences
meta = pd.read_excel(r'X:\users\shangl02\Test\ampliCan.test2\AS15\AS15_submissionform for Eva.xlsx')
meta['RC_Full Amplicon'] = [Bio.Seq.reverse_complement(seq) for seq in meta['Full Amplicon']]
meta['RC_Guide Sequence'] = [Bio.Seq.reverse_complement(seq) if isinstance(seq, str) else seq for seq in meta['Guide Sequence']]

#Add a column to the metadata dataframe indicating the orientation of guide
guide_in_sense = []
for index, row in meta.iterrows():
    if isinstance(row['Guide Sequence'], float):
        check = 'NA'
    elif row['Guide Sequence'] in row['Full Amplicon']:
        check = True
    elif row['Guide Sequence'] in row['RC_Full Amplicon']:
        check = False
    guide_in_sense.append(check)
meta['guide_in_sense'] = guide_in_sense

#Add two columns indicating the 10bp upstream and downstream of protospacer, respectively
flank1 = {}
flank2 = {}
for index, row in meta.iterrows():
    if isinstance(row['Guide Sequence'], float):
        continue
    if row['Guide Sequence'] in row['Full Amplicon']:
        flank1[index] = row['Full Amplicon'].split(row['Guide Sequence'])[0][-10:]
        flank2[index] = row['Full Amplicon'].split(row['Guide Sequence'])[1][:10]       
    elif row['Guide Sequence'] in row['RC_Full Amplicon']:
        flank1[index] = row['Full Amplicon'].split(row['RC_Guide Sequence'])[0][-10:]
        flank2[index] = row['Full Amplicon'].split(row['RC_Guide Sequence'])[1][:10]
flank_df = pd.DataFrame([flank1, flank2], index = ['flank1', 'flank2']).T
meta1 = meta.merge(flank_df, left_index = True, right_index = True)

#Main base calling function
def basecall(directory, filename, guide, flank1, flank2, guide_in_sense):
    seqs = {}
    seq_parse = list(SeqIO.parse(os.path.join(directory, filename), 'fastq'))
    for record in seq_parse:
        recordqual = [x > 31 for x in record.letter_annotations['phred_quality']]
        if float(sum(recordqual)) / float(len(recordqual)) > 0.5:
            recordseq = ''.join([y if x else 'N' for (x,y) in zip(recordqual, record.seq)])
            split1 = recordseq.split(flank1)
            if len(split1) == 2:
                split2 = split1[1].split(flank2)[0]
                if len(split2) == 20:
                    if guide_in_sense == False:
                        seqs[record.id] = Bio.Seq.reverse_complement(split2)
                    else:
                        seqs[record.id] = split2
    frame = pd.DataFrame(list(seqs.items()), columns = ['ID', 'spacer'])
    frame1 = frame[['N' not in x for x in frame.spacer]]
    Motif = motifs.create(frame1.spacer.values)
    raw = pd.DataFrame(Motif.counts, index = [str(i+1) for i in range(20)]).transpose()
    normalized = raw / raw.sum(axis = 0) * 100
    normalized = normalized.round(2)
    guide_split = [*guide]
    per_base = pd.concat([pd.DataFrame(guide_split, index = [str(i+1) for i in range(20)]).transpose(), normalized])
    counts = frame1.groupby('spacer').count().sort_values('ID', ascending = False)
    counts['percentage'] = counts['ID'] / counts['ID'].sum() * 100
    return [per_base, counts]

#Run every fastq file in working directory
res = {}
for index,row in meta1.iterrows():
    for file in filenames:
        if row['Sample Number'] in file:
            res[file] = basecall(directory, file, row['Guide Sequence'], \
                                 row['flank1'], row['flank2'], row['guide_in_sense'])

#there are two outputs from scripts implemented above
# 1. per base percentage of each nucleotide
 

# 2. Counts and percentage of each sequence pattern
 

#Get a single value of base editing efficiency from all samples
be_eff_all = {}
for key, value in res.items():
    sample = key[:7]
    guide = meta1.loc[meta1['Sample Number'] == sample, 'Guide Sequence'].tolist()[0]
    df = value[1][value[1].index != guide]
    eff = df['percentage'][0:5].sum()
    be_eff_all[key[:7]] = eff
be_eff_df = pd.DataFrame(be_eff_all.items(), columns = ['sample', 'base_editing_efficiency'])
#the output
 
