#!/bin/bash
#$ -o $HOME/bud_rnaseq/scripts/Shell/filter_data.stdout.txt
#$ -e $HOME/bud_rnaseq/scripts/Shell/filter_data.stderr.txt
#$ -V

# Filtering the combined files of low quality data:
# qsub $HOME/bud_rnaseq/scripts/Shell/filter_data.sh

# Change directory to where we want to be
cd $HOME

# Create an output directory
outdir="$HOME/bud_rnaseq/reads/filtered"
mkdir -p ${outdir}

# Create the file for the program below to use
i=1
while read line
do
    # Skip the first line of the .csv file
    test ${i} -eq 1 && i=$(expr ${i} + 1) && continue

    # Want to retrive the sample name from the .csv file
    sample_name=$(echo ${line} | cut -d "," -f 8)

    # Create the .txt file
    printf "$HOME/bud_rnaseq/reads/raw/${sample_name}\n" >> "${outdir}/raw_read_file_prefix.txt"

done < $HOME/bud_rnaseq/sample_information.csv

# Apply the filtering

# This script runs fastqc and cutadapt on fastq files, it creates two directories, one with the FastQC reports for the raw and filtered reads and the other with the filtered reads.
# The code can be found at: github.com/tavareshugo/BioPipelines/blob/master/qc_fastq.sh

# Specifiying the number of cores, output, input and filtering used
/home/hugot/code/bioPipelines/qc_fastq.sh -c 6 \
-o "${outdir}" \
-p "${outdir}/raw_read_file_prefix.txt" \
-1 "_R1.fq.gz" -2 "_R2.fq.gz" \
-f "--quality-cutoff 20 --minimum-length 50 --max-n 1"

# Apply multiqc to get a visulisation of all the data in one go (using FastQC) and to see what was cut using cutadapt
multiqc ${outdir} -o ${outdir} -n "multiqc"
# multiqc ${outdir}/fastq/filtered -o ${outdir} -n "multiqc_filtered"
# multiqc ${outdir}/fastq/filtered -o ${outdir} -n "multiqc_cutadapt"
