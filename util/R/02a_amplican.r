#!/usr/bin/env Rscript
suppressPackageStartupMessages(
    library('amplican')
)

args=commandArgs(trailingOnly=TRUE)
config=args[1]
fastq_folder=args[2]
results_folder=args[3]

# # path to example config file
# config <- '/mnt/data/ampliseq/AS14/config_AS14vs2.csv'
# # path to example fastq files
# fastq_folder <- '/mnt/data/ampliseq/AS14/fastq/'
# # output folder, a full path
# results_folder <- '/mnt/data/ampliseq/AS14/output/vs2'


#  run amplican - normal version
amplicanPipeline(config, fastq_folder, results_folder)

#  run amplican only R1 - with higher threshold for control - mitigate index switichn
#amplicanPipeline(config, fastq_folder, results_folder, fastqfiles = 1)
# amplicanPipelineConservative(config, fastq_folder, results_folder, fastqfiles = 1)

# results of the analysis can be found at
message(results_folder)