#!/usr/bin/env python3

#run python script with following command: 
# singularity exec -B /lustre/workspace/home/faste/:/mnt/data /home/faste/singularity/crispr_1.1.sif python3 /home/faste/github/faste-ampliseq/AS14/03_prepoutput.py
'''
singularity exec -B /lustre/workspace/home/faste/:/mnt/data /home/faste/singularity/crispr_1.1.sif python3 /home/faste/github/faste-ampliseq/AS15/03_prepoutput.py --experiment_id AS09.2 --suffix vs15
dfasdfas
vs15 

Lei@02/24/2023
singularity exec -B /hpc/grid/wip_drm_targetsciences/users/shangl02/Test/ampliCan.test/:/mnt/data /lustre/workspace/home/faste/singularity/crispr_1.1.sif python3 03_prepoutput.py --experiment_id AS15 --suffix vs1 

'''

import pandas as pd
import numpy as np
import os
import shutil
import difflib
import warnings
import argparse
import ast
import logging
import math
from logging.handlers import MemoryHandler

pd.set_option('mode.chained_assignment', None)

def init_logger(logger_name, hpc_base_path, exp_id, suffix, loglevel):
    """
    This function will initialise the logger
    :param logger_name: This is the name to be used during logging in the code.
    :type logger_name: string
    :param hpc_base_path: path to the location where output is stored
    :type hpc_base_path: string
    :param exp_id: an ID for the experiment, most common should be something like 'ASXX'
    :type exp_id: string
    :param suffix: a version for the analysis run, e.g. 'vs1'
    :type suffix: string
    :param loglevel: level of logging to be used for console and .log file
    :type loglevel: string

    :return: logger object, memory handler 
    """
    # Initialise Logger
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.DEBUG)

    # this essentially converts 'DEBUG' into '10', the log level into a number
    level=getattr(logging, loglevel)

    # formats the string of logging
    format_string = "%(asctime)s %(name)s [%(levelname)s] %(message)s"
    formatter = logging.Formatter(format_string)

    # Defines the file handler 
    file_handler = logging.FileHandler(f'{hpc_base_path}/{exp_id}/{suffix}_{logger_name}.log', mode='w')
    file_handler.setFormatter(formatter)
    file_handler.setLevel(level)
    logger.addHandler(file_handler)
    
    # Setting up stream handler to display messages on the console
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    stream_handler.setLevel(level)
    logger.addHandler(stream_handler)

    # Memory handler to keep track of error messages that occur
    memory_handler = MemoryHandler(capacity=0)
    memory_handler.setLevel(logging.ERROR)
    logger.addHandler(memory_handler)

    return logger, memory_handler

    df_sum = pd.read_csv(f'{hpc_base_path}/{exp_id}/output/{suffix}/config_summary.csv')
    df_reads = pd.read_csv(f'{hpc_base_path}/{exp_id}/output/{suffix}/barcode_reads_filters.csv')

def generic_output_file(df_sum, df_reads):
    """
    This function generates as simplified output dataframe from the Amplican output. 
    
    :param df_sum: Amplican generated dataframe that contains a summary of the run with overall edits
    :type df_sum: pandas dataframe
    :param df_reads: Amplican generated dataframe that contains summary of read information and filtering
    :type df_reads: pandas dataframe
    
    :return: Pandas dataframe
    """
    #make into one dataframe
    df = pd.merge(df_sum, df_reads, how = 'inner', on='Barcode')

    # calculate percentages
    df['Percent.Reads_Filtered'] = (df['Reads_Filtered']/df['read_count'])*100
    df['Percent.Reads_Edited'] = (df['Reads_Edited']/df['Reads_Filtered'])*100
    df['Percent.Reads_Frameshifted'] = (df['Reads_Frameshifted']/df['Reads_Filtered'])*100
    df['Percent.HDR'] = (df['HDR']/df['Reads_Filtered'])*100
    
    df = df[['ID', 'Barcode', 'Control','read_count','Reads_Filtered','Percent.Reads_Filtered',
         'Reads_Del','Reads_In','Reads_Edited','Reads_Frameshifted','HDR','Percent.Reads_Edited',
         'Percent.Reads_Frameshifted','Percent.HDR']]

    df.columns = ['Sample Name', 'Amplicon', 'Control', 'Total_Reads', 'Reads_Filtered', 'Percent.Reads_Filtered',
                  'Reads_Del', 'Reads_In', 'Reads_Edited', 'Reads_Frameshifted', 'HDR', 'Percent.Reads_Edited',
                  'Percent.Reads_Frameshifted', 'Percent.HDR']
    return df

def replace_nan(row, columnid):
    """
    Function that fills in the NaN with number of unedited reads for each nucleotide and postion
    
    :param columnid: a nucleotide
    :type columnid: string
    :param sample: row or position within the editing window
    :type exp_id: index? not sure
    
    :return: applied as a lambda function 
    """
    # define all nucleotides as a set, calulate set that is not the columnid
    rem_nuc = set(['A', 'C', 'T', 'G']) - set([columnid])
    
    if (math.isnan(row[columnid])):
        return row['total'] - (row[list(rem_nuc)[0]] + row[list(rem_nuc)[1]] + row[list(rem_nuc)[2]])
    else:
        return row[columnid]

def save_df_baseedit(df_norm, df_conf, df_sum, hpc_base_path, exp_id, suffix):
    """
    This function generates + saves a .csv for each base editing sample (saved within a new folder) with the editing purity whichi is the % of each nucleotide at each position accross the entire region that was assayed for edting.
    
    :param df_norm: Amplican generated dataframe that contains all individual editing events and read counts per sample in single table, events are normalized and shifted (first postion at 5' end of editing window is '0')
    :type df_norm: pandas dataframe
    :param df_conf: Config dataframe used for Amplican submission
    :type df_conf: pandas dataframe
    :param df_sum: Amplican generated dataframe that contains a summary of the run with overall edits
    :type df_sum: pandas dataframe
    :param hpc_base_path: path to the location where output is stored
    :type hpc_base_path: string
    :param exp_id: an ID for the experiment, most common should be something like 'ASXX'
    :type exp_id: string
    :param suffix: a version for the analysis run, e.g. 'vs1'
    :type suffix: string
    
    :return: None
    """
    # remove the non base edting samples from the input list
    samples = df_conf.dropna(subset=['Donor'])['ID'].tolist()

    # create new directory for the edit percentages
    newpath_output = os.path.join(hpc_base_path, exp_id, 'output', suffix, 'edit_percent') 

        # create this output path in case it doesn't exist    
    if not os.path.exists(newpath_output):
        os.makedirs(newpath_output)
        logger.info("New output path was created")

    # loop through the samples and create table for each sample
    for sample in samples:
        #select only one sample
        dftmp = df_norm[(df_norm['overlaps'] == True) & \
        (df_norm['seqnames'] == sample) & \
        (df_norm['type'] == 'mismatch')].copy()

        #pivot table into right format
        dftmp = pd.pivot_table(dftmp, values='counts', index=['start', 'originally'],
                    columns=['replacement'], aggfunc=np.sum)
       
        dftmp['total'] = df_sum[df_sum['ID'] == sample]['Reads_Filtered'].tolist()[0]

        cols = ['A', 'C', 'T', 'G']

        for col in cols:
            dftmp[col] = dftmp.apply(lambda row: replace_nan(row, col), axis=1)
            
        dftmp['%A'] = dftmp['A']/dftmp['total']
        dftmp['%C'] = dftmp['C']/dftmp['total']
        dftmp['%G'] = dftmp['G']/dftmp['total']
        dftmp['%T'] = dftmp['T']/dftmp['total']

        dftmp = dftmp[['%A', '%C', '%G', '%T']]

        dftmp.to_csv(f'{hpc_base_path}/{exp_id}/output/{suffix}/edit_percent/{sample}.csv', index=True)

def extract_indels(df_norm, sample, position, targetbase, nontargetbase):
    """
    This function generates two dataframes depending if they contain indels or not for a specific sample
    
    :param df_norm: Amplican generated dataframe that contains all individual editing events and read counts per sample in single table, events are normalized and shifted (first postion at 5' end of editing window is '0')
    :type df_norm: pandas dataframe
    :param sample: Sample name
    :type exp_id: string
    :param position: position of SNP
    :type position: integer
    :param targetbase: nucleotide that is being targeted by gRNA
    :type targetbase: string
    :param nontargetbase: nucleotide that is NOT being targeted by gRNA
    :type nontargetbase: string
    
    :return: Pandas dataframes
    """
    # filter to particular sample and also only events overlapping the window and that have consensus between forward and reverse read
    df = df_norm[(df_norm['seqnames'] == sample) & \
                (df_norm['overlaps'] == True) & \
                (df_norm['consensus'] == True)]
    
    #extract the read_ids that are either deletions or insertions
    indellist = df[(df['type'] == 'deletion')|(df['type'] == 'insertion')]['read_id'].tolist()
    
    # create two new dfs with events containing indels (dfi) and events not containing indels (dfni)
    dfi = df[df['read_id'].isin(indellist)]
    dfni = df[~df['read_id'].isin(indellist)]
    
    return dfi, dfni


def allelespec_output_file(df_norm, df_conf, df_sum):
    """
    This function generates a dataframe that contains allele specific editing events in 3 categoreis, on target indels (1), nontarget indels (2) and non assignable indels (3). Both raw numbers and percentages are provided. 
    
    :param df_norm: Amplican generated dataframe that contains all individual editing events and read counts per sample in single table, events are normalized and shifted (first postion at 5' end of editing window is '0')
    :type df_norm: pandas dataframe
    :param df_conf: Config dataframe used for Amplican submission
    :type df_conf: pandas dataframe
    :param df_sum: Amplican generated dataframe that contains a summary of the run with overall edits
    :type df_sum: pandas dataframe
    
    :return: Pandas dataframe
    """

    replacement_dictionary_caps = {'A':'C', 'G':'T', 'C':'A', 'T':'G'}
    snp_dictionary = {'A' : 'G', 'G' : 'A', 'C' : 'T', 'T': 'C'}

    samplelist = df_conf.dropna(subset=['SNP_position'])['ID'].tolist()

    results_lists=[]
    
    for sample in samplelist:
        # get snp position for that particular sample (subtract one so Biologists don't have to worry about 0 indexing)
        position = int(df_conf[df_conf['ID']==sample]['SNP_position'].tolist()[0]) - 1
        
        # change back to original nucleotide and then generate nontarget using dictionaries
        replacement_nuc = df_conf[df_conf['ID']==sample]['guideRNA'].tolist()[0][position]
        
        targetbase = replacement_dictionary_caps[replacement_nuc]
        nontargetbase = snp_dictionary[targetbase]
        
        dfi, dfni = extract_indels(df_norm, sample, position, targetbase, nontargetbase)
        
        targetSNP_noedit = dfni[(dfni['start'] == position) & (dfni['replacement'] == targetbase)]['counts'].sum()
        nontargetSNP_noedit = dfni[(dfni['start'] == position) & (dfni['replacement'] == nontargetbase)]['counts'].sum()
        notassigned_noedit = df_sum[df_sum['ID']==sample]['Reads_Filtered'].tolist()[0] -\
                            ((dfi.drop_duplicates(subset=['read_id'])['counts'].sum()) +\
                            (dfni.drop_duplicates(subset=['read_id'])['counts'].sum()))
        targetSNP_indel = dfi[(dfi['start'] == position) & (dfi['replacement'] == targetbase)]['counts'].sum()
        nontargetSNP_indel = dfi[(dfi['start'] == position) & (dfi['replacement'] == nontargetbase)]['counts'].sum()
        notassigned_indel = dfi.drop_duplicates(subset=['read_id'])['counts'].sum() -\
                        (dfi[(dfi['start'] == position) & (dfi['replacement'] == targetbase)]['counts'].sum() +\
                            dfi[(dfi['start'] == position) & (dfi['replacement'] == nontargetbase)]['counts'].sum())
        
        prop_targetSNP_noedit = targetSNP_noedit/df_sum[df_sum['ID']==sample]['Reads_Filtered'].tolist()[0]
        prop_nontargetSNP_noedit = nontargetSNP_noedit/df_sum[df_sum['ID']==sample]['Reads_Filtered'].tolist()[0]
        prop_notassigned_noedit = notassigned_noedit/df_sum[df_sum['ID']==sample]['Reads_Filtered'].tolist()[0]
        prop_targetSNP_indel = targetSNP_indel/df_sum[df_sum['ID']==sample]['Reads_Filtered'].tolist()[0]
        prop_nontargetSNP_indel = nontargetSNP_indel/df_sum[df_sum['ID']==sample]['Reads_Filtered'].tolist()[0]
        prop_notassigned_indel = notassigned_indel/df_sum[df_sum['ID']==sample]['Reads_Filtered'].tolist()[0]
        
        results_dic = {'ID': sample, 'targetSNP_noedit':targetSNP_noedit, 'nontargetSNP_noedit':nontargetSNP_noedit, 'notassigned_noedit':notassigned_noedit, 
                                    'targetSNP_indel':targetSNP_indel, 'nontargetSNP_indel':nontargetSNP_indel, 'notassigned_indel':notassigned_indel,
                                    'proportion_targetSNP_noedit':prop_targetSNP_noedit, 'proportion_nontargetSNP_noedit':prop_nontargetSNP_noedit, 'proportion_notassigned_noedit':prop_notassigned_noedit, 
                                    'proportion_targetSNP_indel':prop_targetSNP_indel, 'proportion_nontargetSNP_indel':prop_nontargetSNP_indel, 'proportion_notassigned_indel':prop_notassigned_indel}
        results_lists.append(results_dic)
    
    results_df = pd.DataFrame(results_lists)    
    results_df = results_df.set_index('ID') 

    return results_df

def parse_args():
    """
    Use the argparse library to provide special help text and simplify argument handling

    :return: tuple
    """
    parser = argparse.ArgumentParser(
        description='Check submission file. This can take multiple arguments')
    parser.add_argument('--base_path', dest='hpc_base_path', required=False,
                        default='/mnt/data/ampliseq',
                        help="The path where experiments are stored")
    parser.add_argument('--experiment_id', dest='exp_id', required=True,
                        help="The experiment id should start with 'AS', also corresponds to file names")
    parser.add_argument('--suffix', dest='suffix', required=False,
                        default='vs1',
                        help="Suffix if you are rerunning the analysis")
    parser.add_argument('--log', dest='loglevel', required=False,
                        default='DEBUG', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        help="Set level of logging")

    args = parser.parse_args()

    return args.hpc_base_path, args.exp_id, args.suffix, args.loglevel

if __name__ == "__main__":
    try:
        # Reading the arguments and setting up variables
        hpc_base_path, exp_id, suffix, loglevel = parse_args()

        # debug 
        #hpc_base_path = r'X:\projects\p073_GeneEditingSupport\test\test.03192023'
        #exp_id='AS15'
        #suffix='vs1'
        #loglevel='DEBUG'

        
        # Initiate logger
        logger, memory_handler = init_logger("output-processing-logs", hpc_base_path, exp_id, suffix, loglevel)
        logger.info("Starting output-processing step.")

        logger.info("The path of the experiment is : %s ." % hpc_base_path) 
        logger.info("The experiment ID is : %s ." % exp_id)
        logger.info("The suffix for this run is: %s ." % suffix)
 
        # read in files
        logger.info("Read in files.")

        df_norm = pd.read_csv(os.path.join(hpc_base_path, exp_id, 'output', suffix, 'alignments', 'events_filtered_shifted_normalized.csv'))
        df_conf = pd.read_csv(os.path.join(hpc_base_path, exp_id, f'config_{exp_id}{suffix}.csv'))
        df_sum  = pd.read_csv(os.path.join(hpc_base_path, exp_id, 'output', suffix, 'config_summary.csv'))
        df_reads = pd.read_csv(os.path.join(hpc_base_path, exp_id, 'output', suffix, 'barcode_reads_filters.csv'))

        # get general output
        df_gen_output = generic_output_file(df_sum, df_reads)
        # df_gen_output.to_csv(f'{hpc_base_path}/{exp_id}/output/{suffix}/{exp_id}_output_summary_{suffix}.csv', index=False)
        df_gen_output.to_csv(os.path.join(hpc_base_path, exp_id, 'output', suffix, f'{exp_id}_output_summary_{suffix}.csv'), index=False)


        # if there are any base editing samples save edit percentages
        if df_conf['Donor'].any():
            save_df_baseedit(df_norm, df_conf, df_sum, hpc_base_path, exp_id, suffix)

        # save allele specfic editing
        if df_conf['SNP_position'].any():
            df_allelespec_output = allelespec_output_file(df_norm, df_conf, df_sum)
            df_allelespec_output.to_csv(f'{hpc_base_path}/{exp_id}/output/{suffix}/{exp_id}_output_summary_allelespecificedit_{suffix}.csv', index=True)
    
        # this Exception get's raised if there are any samples that fail the formatting criteria   
        for record in memory_handler.buffer:
            if record.levelno >= logging.ERROR:
                raise Exception(f"One or more error or higher level messages logged: {record.getMessage()}")
    

        logger.info("The output file generation is complete.")  
    
    except:
        logger.error(
            "The output file generation has failed. Kindly check.", exc_info=True)