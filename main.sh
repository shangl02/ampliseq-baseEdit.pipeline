#!/bin/bash

## args
expID=$1
dataPath=$2
formPath=$3
outDir=$4

suffix="vs1"

## code base directory
dir_codeBase=$( dirname -- "$0" )
#echo "code base dir: ${dir_codeBase}"

## module load
ml ib fastqc/0.11.9

## dir setting
outDir2="${outDir}/${expID}"
dir_fastq="${expID}/fastq"
dir_fastq_abs="${outDir}/${dir_fastq}"
dir_fastqc="${expID}/fastqc"
dir_fastqc_abs="${outDir}/${dir_fastqc}"

echo "[Info] Setting up directories"
mkdir "${outDir2}"
mkdir "${dir_fastq_abs}"
mkdir "${dir_fastqc_abs}"
echo -e '[Success] fastq and QC directories created\n'

## copy raw data
echo "[Info] Copying raw data"
scp -r "${dataPath}/Fastq/." "${dir_fastq_abs}"
echo -e "[Success] fastq data are copied\n"

## fastqc
echo "[Info] fastQC"
for i in "${dir_fastq_abs}/*.fastq.gz"; do
    echo $i;
    fastqc -o ${dir_fastqc_abs} $i -t 4 -q;
done
echo -e "[Success] fastQC is complete\n"

## multiqc
echo "[Info] multiQC"
dir_multiqc="${expID}/multiqc"
dir_multiqc_abs="${outDir}/${dir_multiqc}"
mkdir "${dir_multiqc_abs}"
singularity exec -B ${outDir2}:/mnt/data ${dir_codeBase}/singularity/multiqc_latest.sif multiqc "/mnt/data/fastqc/" -o "/mnt/data/multiqc/"
echo -e "[Success] multiQC is complete\n"

## prepare input
echo "[Info] Preparing input for ampliCan pipeline"
cp "${formPath}" "${outDir2}"
formFile=$( basename -- "${formPath}" )
formName=${formFile%.*}
singularity exec -B \
    ${outDir}:/mnt/data,${dir_codeBase}/util/python:/mnt/code ${dir_codeBase}/singularity/crispr_1.1.sif \
    python3 /mnt/code/01_prepinput.py \
        --base_path /mnt/data \
        --experiment_id "${expID}" \
        --inputfile_id "${formName}" \
        --suffix ${suffix} \
        --log ERROR
echo -e "[Success] Preparation of ampliCan is complete\n"

## run ampliCan
echo "[Info] run ampliCan pipeline"
bsub -q medium -Is singularity exec -B \
    ${outDir}:/mnt/data,${dir_codeBase}/util/R:/mnt/code ${dir_codeBase}/singularity/amplican_1.16.0.sif \
    Rscript /mnt/code/02a_amplican.r \
        /mnt/data/${expID}/config_${expID}${suffix}.csv \
        /mnt/data/${dir_fastq} \
        /mnt/data/${expID}/output/${suffix}
echo -e "[Success] ampliCan pipeline running is complete\n"

## prepare output
echo "[Info] Prepare output"
singularity exec -B \
    ${outDir}:/mnt/data,${dir_codeBase}/util/python:/mnt/code ${dir_codeBase}/singularity/crispr_1.1.sif \
    python3 /mnt/code/03_prepoutput.py \
        --base_path /mnt/data \
        --experiment_id "${expID}" \
        --suffix "${suffix}" 
echo -e "[Success] Output preparation is complete\n"

## base editing advanced module
## will be implemented later

## base editing freq plot
echo "[Info] Plotting"
singularity exec -B \
    ${outDir2}:/mnt/data,${dir_codeBase}/util/python:/mnt/code ${dir_codeBase}/singularity/crispr_1.1.sif \
    python3 /mnt/code/baseEditFreqPlot.py \
        "/mnt/data/output/${suffix}/edit_percent"
echo -e "[Success] Plotting is complete\n" 