#!/usr/bin/env python3

# testing
'''
singularity exec -B /lustre/workspace/home/faste/:/mnt/data\
 /home/faste/singularity/crispr_1.1.sif\
 python3 /home/faste/github/faste-ampliseq/AS15/01_prepinput.py --experiment_id AS13 --inputfile_id AS13\ 20230112\ AmpSeq\ CMT1A\ Schwann-like\ cells_FINAL --suffix vstest --log ERROR


singularity exec -B /lustre/workspace/home/faste/:/mnt/data\
 /home/faste/singularity/crispr_1.1.sif\
 python3 /home/faste/github/faste-ampliseq/AS15/01_prepinput.py --experiment_id AS14 --inputfile_id AS14_submissionform\ for\ Eva --suffix vstest --log ERROR

singularity exec -B /lustre/workspace/home/faste/:/mnt/data\
 /home/faste/singularity/crispr_1.1.sif\
 python3 /home/faste/github/faste-ampliseq/AS15/01_prepinput.py --experiment_id AS14 --inputfile_id AS14_submissionform\ for\ Eva_230207 --suffix vstest --log ERROR

singularity exec -B /lustre/workspace/home/faste/:/mnt/data\
 /home/faste/singularity/crispr_1.1.sif\
 python3 /home/faste/github/faste-ampliseq/AS15/01_prepinput.py --experiment_id AS15 --inputfile_id AS15_submissionform\ for\ Eva --suffix vs1 --log ERROR


AS15_submissionform\ for\ Eva

singularity exec -B /lustre/workspace/home/faste/:/mnt/data\
 /home/faste/singularity/crispr_1.1.sif\
 python3 /home/faste/github/faste-ampliseq/AS15/01_prepinput.py --experiment_id AS14 --inputfile_id piplinetesting2 --suffix vstest

piplinetesting.xlsx
AS13\ 20230112\ AmpSeq\ CMT1A\ Schwann-like\ cells_FINAL


Lei@2/22/2023
singularity exec -B /hpc/grid/wip_drm_targetsciences/users/shangl02/Test/ampliCan.test/:/mnt/data /lustre/workspace/home/faste/singularity/crispr_1.1.sif python3 01_prepinput.py --experiment_id AS15 --inputfile_id AS15_submissionform\ for\ Eva --suffix vs1 --log ERROR


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
from logging.handlers import MemoryHandler

pd.set_option('mode.chained_assignment', None)


# function to reverse complement a sequence
def complement(seq):
    """
    This function complements a sequence
    :param seq: a nucleotide sequence
    :type seq: string

    :return: reverse of a nucleotide sequence.
    """
    complement = {'a': 't', 'c': 'g', 'g': 'c', 't': 'a'} 
    bases = list(seq) 
    bases = [complement[base] for base in bases] 
    return ''.join(bases)

def reverse_complement(s):
    """
    This function reverse complements a sequence
    :param s: a nucleotide sequence
    :type s: string

    :return: reverse complement of a nucleotide sequence.
    """
    return complement(s[::-1])

def splittooligos(fullstring, substring):
    """
    This function splits a nucleotide sequence (fullstring) in shorter nucleotide sequences based on length of substring
    :param fullstring: a nucleotide sequence
    :type fullstring: string
    :param substring: a nucleotide sequence, shorter than fullstring
    :type substring: string

    :return: a list of all possible substrings of fullstring with the length of substring
    """
    oligo_list = []
    
    for i in range(len(fullstring)-len(substring)):
        counter = i
        oligo = fullstring[counter:counter+len(substring)]
        oligo_list.append(oligo)
    return oligo_list

def save_output_file(hpc_base_path, exp_id, exp_version):
    """
    This function saves a simplified outpout table (.csv) of the CRISPR editing experiment
    :param hpc_base_path: path to the location where output is stored
    :type hpc_base_path: string
    :param exp_id: an ID for the experiment, most common should be something like 'ASXX'
    :type exp_id: string
    :param exp_version: a version for the analysis run, e.g. 'vs1'
    :type exp_version: string

    :return: df
    """
    # read in output tables from Amplican and Crispresso
    df_r = pd.read_csv(f'{hpc_base_path}/{exp_id}/output/{exp_version}/config_summary.csv')
    df_s = pd.read_csv(f'{hpc_base_path}/{exp_id}/output/{exp_version}/barcode_reads_filters.csv')

    #make into one dataframe
    df = pd.merge(df_r, df_s, how = 'inner', on='Barcode')

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

    df.to_csv(f'{hpc_base_path}/{exp_id}/output/{exp_version}/{exp_id}_output_summary_{exp_version}.csv', index=False)


def reformat_config(df, exp_id, hpc_base_path, suffix):
    """
    This function generates a config file (dataframe that can be saved as .csv) to be used for Amplican. It takes a simple input file submitted by the gene editing team and transforms it into a config.csv. 
    Specifically tranformations include:
    - it harmonizes the column headers
    - generates the mapping between the sample identifier and the .fastq file names
    - generates Amplican compatible input format of the Amplicon reference sequence in lower case and the gRNA in upper case letters
    - if a single control is being used for multiple gRNAs it generates 'pseudoduplicate' samples so that background editing is assayed in the same region than the treated samples
    - for allele specific editing it substitutes the target base with a 'pseudobase' that does not match either reference or alternate allele which allows for allele specific edit quantification

    :param df: a dataframe that contains the submission information
    :type df: pandas dataframe
    :param exp_id: an ID for the experiment, most common should be something like 'ASXX'
    :type exp_id: string
    :param hpc_base_path: path to the location where output is stored
    :type hpc_base_path: string
    :param suffix: a version for the analysis run, e.g. 'vs1'
    :type suffix: string

    :return: Pandas dataframe
    """
    # make the input .csv for amplican
    #def reformat_config(df, exp_id, hpc_base_path, suffix, split_fastq = False):

    filenames = os.listdir(f'{hpc_base_path}/{exp_id}/fastq/')

    #select relevant columns
    ## this is new
    df = df[['Sample Number', 'Sample Name', 'Experiment #', 'Full Amplicon','Expected Sequence (if base editor)','Guide Sequence','Fwd primer bind','Rev primer bind', 'SNP_position']]

    #split df in na and not na (control samples are the ones that do not have a gRNA)
    dff_all = df[df['Guide Sequence'].notna()]
    dfna_all = df[~df['Guide Sequence'].notna()]

    dff_all_copy = dff_all.copy()

    dff_all.loc[:, 'Sample Number new'] = dff_all_copy['Sample Number']
    dff_all.loc[:, 'Control'] = 0

    # create results list
    results_lists = []

    experiments = dfna_all['Experiment #'].unique().tolist()
    # loop through experiments so that only within experiment comparisons are being performed
    for experiment in experiments:

        # filter samples to only experiment of interest
        dfna = dfna_all[dfna_all['Experiment #'] == experiment]
        dff = dff_all[dff_all['Experiment #'] == experiment]

        ctsamples = dfna['Sample Number'].tolist()

        # loop through amplicons matching the control samples to the corresponding edited samples by the amplicon (samples with the same amplicon should be matched)
        for ctsample in ctsamples:
            amplicon = dfna[dfna['Sample Number'] == ctsample]['Full Amplicon'].tolist()[0]

            # scenario for non base editors
            if pd.isna(dff[dff['Full Amplicon']== amplicon]['Expected Sequence (if base editor)'].unique().tolist()[0]) == True:
                grnas = dff[dff['Full Amplicon']== amplicon]['Guide Sequence'].unique().tolist()
                for grna in grnas:
                    sample_number = ctsample
                    sample_number_new = f'{ctsample}_{grna}'
                    sample_name = f'{dfna[dfna["Sample Number"] == ctsample]["Sample Name"].tolist()[0]}_{grna}_{ctsample}'
                    fwd_primer = dff[(dff['Full Amplicon']== amplicon) & (dff['Guide Sequence']== grna)]['Fwd primer bind'].unique().tolist()[0]
                    rev_primer = dff[(dff['Full Amplicon']== amplicon) & (dff['Guide Sequence']== grna)]['Rev primer bind'].unique().tolist()[0]
                    snp_pos = dff[(dff['Full Amplicon']== amplicon) & (dff['Guide Sequence']== grna)]['SNP_position'].unique().tolist()[0]
                    expected_amplicon = ""
                    control = 1

                    # make a dictonary
                    results_dic = {'Sample Number': sample_number, 'Sample Number new': sample_number_new, 'Sample Name': sample_name,'Experiment #': experiment, 'Full Amplicon': amplicon,
                           'Expected Sequence (if base editor)': expected_amplicon, 'Guide Sequence': grna,'Fwd primer bind': fwd_primer,'Rev primer bind': rev_primer, 'SNP_position': snp_pos, 'Control': 1}
                    results_lists.append(results_dic)

            # scenario for base editors    
            else:
                grnas = dff[dff['Full Amplicon']== amplicon]['Guide Sequence'].unique().tolist()
                for grna in grnas:

                    sample_number = ctsample
                    sample_number_new = f'{ctsample}_{grna}'
                    sample_name = f'{dfna[dfna["Sample Number"] == ctsample]["Sample Name"].tolist()[0]}_{grna}_{ctsample}'
                    fwd_primer = dff[(dff['Full Amplicon']== amplicon) & (dff['Guide Sequence']== grna)]['Fwd primer bind'].unique().tolist()[0]
                    rev_primer = dff[(dff['Full Amplicon']== amplicon) & (dff['Guide Sequence']== grna)]['Rev primer bind'].unique().tolist()[0]
                    expected_amplicon = dff[dff['Guide Sequence']== grna]['Expected Sequence (if base editor)'].tolist()[0]
                    snp_pos = dff[(dff['Full Amplicon']== amplicon) & (dff['Guide Sequence']== grna)]['SNP_position'].unique().tolist()[0]
                    control = 1

                    # make a dictonary
                    results_dic = {'Sample Number': sample_number, 'Sample Number new': sample_number_new, 'Sample Name': sample_name,'Experiment #': experiment, 'Full Amplicon': amplicon,
                           'Expected Sequence (if base editor)': expected_amplicon, 'Guide Sequence': grna,'Fwd primer bind': fwd_primer,'Rev primer bind': rev_primer, 'SNP_position': snp_pos, 'Control': 1}
                    results_lists.append(results_dic)
                    
    # create a DataFrame from the list of dictionaries
    results_dfna = pd.DataFrame(results_lists)

    # combine new df of control samples with edited samples
    df = pd.concat([dff_all, results_dfna])

    # replace NaNs with empty values
    df['Expected Sequence (if base editor)'].fillna('', inplace=True)

    amplicon_list = []
    baseedit_list = []
    direction_list = []
    f_read_list = []
    r_read_list = []
    grna_list = []
    control_list = []

    #loop through the rows of the new dataframe and perform several changes like UPPER case gRNA within Amplicon sequence and updating column names
    for index, row in df.iterrows():
        amplicon = 0
        baseedit = 0
        direction = 0
        f_read = 0
        rev_read = 0
        grna = 0
        control = 0

        # look for gRNA in amplicon sequence and revcomplement, make it upper case
        fullstring = row['Full Amplicon'].lower()
        substring = row['Guide Sequence'].lower()
        grna = substring.upper()

        # change SNP if have a position list in submission form, allele specific editing:
        if pd.isna(row['SNP_position']) == False:
            
            # this exchanges the SNP if there is a SNP position + makes it into integer
            position = int(row['SNP_position']) - 1
            replacement_dictionary = {'a':'c', 'g':'t', 'c':'a', 't':'g'}
            new_character = replacement_dictionary[substring[position]]
            new_string = substring[:position] + new_character + substring[position+1:]
            grna = new_string.upper()
            
            if substring in fullstring:
                amplicon = fullstring.replace(substring, new_string.upper())
                direction = 0    
            else:
                new_string_revco = reverse_complement(new_string)
                substring_revco = reverse_complement(substring)
                amplicon = fullstring.replace(substring_revco, new_string_revco.upper())
                direction = 1   
            baseedit = ""
            # all samples even controls are set to '0' so that SNPs don't get substracted from samples
            control = 0

        # make sequence depending if have baseeditor      
        elif row['Expected Sequence (if base editor)'] != '':         
            fullstring_baseedit = row['Expected Sequence (if base editor)'].lower()
                
            if substring in fullstring:
                amplicon = fullstring.replace(substring, substring.upper())
                direction = 0
                oligo_list = splittooligos(fullstring_baseedit, substring)
                gRNA_match = difflib.get_close_matches(substring, oligo_list, n=1, cutoff=0.6)    
                baseedit = fullstring_baseedit.replace(gRNA_match[0], gRNA_match[0].upper())
            else:
                substring_revco = reverse_complement(substring)
                amplicon = fullstring.replace(substring_revco, substring_revco.upper())
                direction = 1
                oligo_list = splittooligos(fullstring_baseedit, substring_revco)
                gRNA_match = difflib.get_close_matches(substring_revco, oligo_list, n=1, cutoff=0.6)
                baseedit = fullstring_baseedit.replace(gRNA_match[0], gRNA_match[0].upper()) 
            control = row['Control']
            
        else:  
            if substring in fullstring:
                amplicon = fullstring.replace(substring, substring.upper())
                direction = 0
            else:
                substring_revco = reverse_complement(substring)
                amplicon = fullstring.replace(substring_revco, substring_revco.upper())
                direction = 1
            baseedit = ""
            control = row['Control']
        # make columns with sequencing file names
        idn = row['Sample Number'].split('_')[1]
        idname = row['Sample Number'].split('_')[0]
        wantedfilename = f'{idname}-{idn}_S'

        result = list(filter(lambda x: x.startswith(wantedfilename), filenames))
        result.sort() # sort so that the R1 is first

        f_read = result[0]
        r_read = result[1]

        f_read_list.append(f_read)
        r_read_list.append(r_read)

        # this sets the controls in allele specific editing to 0 (so that the SNPs are not subtracted)
        control_list.append(control)
        
        grna_list.append(grna)
        amplicon_list.append(amplicon)
        direction_list.append(direction)

        baseedit_list.append(baseedit)

        # some input checks on the rows to make sure there are no mistakes in submission form
        # check if gRNA is part of the Amplicon sequence
        if grna.lower() not in amplicon.lower():
            grna_revco = reverse_complement(grna.lower())
            if grna_revco not in amplicon.lower():
                logger.error("For sample %s the guideRNA can not be found in the Amplicon sequence" % row['Sample Number'])
        
        # check if primers are part of the Amplicon sequence
        if row['Fwd primer bind'].lower() not in amplicon:
            logger.error("For sample %s the forward primer cannot be found in the Amplicon sequence" % row['Sample Number'])

        if reverse_complement(row['Rev primer bind'].lower()) not in amplicon:
            logger.error("For sample %s the reverse primer cannot be found in the Amplicon sequence" % row['Sample Number'])

        # check if primers are actually in the beginning and end of the Amplicon sequence
        #forward primer
        string1 = row['Fwd primer bind'].lower()
        string2 = amplicon[:len(row['Fwd primer bind'])]
        if string1 != string2:
            logger.warning("For sample %s the forward primer is not at the start of the Amplicon sequence" % row['Sample Number'])

        #reverse primer
        string1 = reverse_complement(row['Rev primer bind'].lower())
        string2 = amplicon[-len(row['Rev primer bind']):]
        if string1 != string2:
            logger.warning("For sample %s the reverse primer is not at the end of the Amplicon sequence" % row['Sample Number'])

        # check if Amplicon and HDR sequence are the same lenght  
        if row['Expected Sequence (if base editor)'] != '':
            string1 = len(amplicon)
            string2 = len(row['Expected Sequence (if base editor)'])
            if string1 != string2:
                logger.error("The amplicon sequence and the HDR template have unequal lenghts for sample %s - please check input" % row['Sample Number'])

    df['Amplicon'] = amplicon_list
    df['Direction'] = direction_list

    df['Forward_Reads'] = f_read_list
    df['Reverse_Reads'] = r_read_list

    df['Forward_Primer'] = df['Fwd primer bind']
    df['Reverse_Primer'] = df['Rev primer bind']
    df['ID'] = df['Sample Name']
    df['Barcode'] = df['Sample Number new']
    df['Group'] = df['Experiment #']

    df['Control'] = control_list

    df['guideRNA'] = grna_list
    df['Donor'] = baseedit_list

    df = df[['ID','Barcode', 'Forward_Reads','Reverse_Reads','Group','Control','guideRNA','Forward_Primer',
    'Reverse_Primer','Direction','Amplicon','Donor', 'SNP_position']]

    return df

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
                        help="The experiment id should start with AS, also corresponds to file names")
    parser.add_argument('--inputfile_id', dest='nameinputfile', required=True,
                        help="Name of the input file, should be an Excel file")
    parser.add_argument('--suffix', dest='suffix', required=False,
                        default='vs1',
                        help="Suffix if you are rerunning the analysis")
    parser.add_argument('--log', dest='loglevel', required=False,
                        default='DEBUG', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        help="Set level of logging")

    args = parser.parse_args()

    return args.hpc_base_path, args.exp_id, args.nameinputfile, args.suffix, args.loglevel

if __name__ == "__main__":
    try:
        # Reading the arguments and setting up variables
        hpc_base_path, exp_id, nameinputfile, suffix, loglevel = parse_args()
        
        # Initiate logger
        logger, memory_handler = init_logger("preprocessing-logs", hpc_base_path, exp_id, suffix, loglevel)
        logger.info("Starting pre-processing step.")

        logger.info("The path of the experiment is : %s ." % hpc_base_path)
        logger.info("The experiment ID is : %s ." % exp_id)
        logger.info("The name of the input file is %s.xlsx" % nameinputfile)
        logger.info("The suffix for this run is: %s ." % suffix)
 
       
        # create new output path
        newpath_output = f'{hpc_base_path}/{exp_id}/output' 

        # create this output path in case it doesn't exist    
        if not os.path.exists(newpath_output):
            os.makedirs(newpath_output)
            logger.info("New output path was created")

        newpath_output_suffix = f'{hpc_base_path}/{exp_id}/output/{suffix}' 
            
        if not os.path.exists(newpath_output_suffix):
            os.makedirs(newpath_output_suffix)
            logger.info("New output path with suffix was created")

        df = pd.read_excel(f'{hpc_base_path}/{exp_id}/{nameinputfile}.xlsx', engine='openpyxl')
        logger.info("File has been read in")

        df = reformat_config(df, exp_id, hpc_base_path, suffix) 

        # this Exception get's raised if there are any samples that fail the formatting criteria   
        for record in memory_handler.buffer:
            if record.levelno >= logging.ERROR:
                raise Exception(f"One or more error or higher level messages logged: {record.getMessage()}")
        

        df.to_csv(f'{hpc_base_path}/{exp_id}/config_{exp_id}{suffix}.csv', index=False)

        logger.info("The pre-processing is complete.")  
    
    except:
        logger.error(
            "The pre-processing has failed. Kindly check.", exc_info=True)