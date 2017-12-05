#!/bin/bash
#$ -o $HOME/bud_rnaseq/scripts/Shell/combining_data.stdout.txt
#$ -e $HOME/bud_rnaseq/scripts/Shell/combining_data.stderr.txt

#Combining the 4 files for each sample into one:
# qsub $HOME/bud_rnaseq/scripts/Shell/combining_data.sh

#Change directory to where we want to be
cd $HOME

#Create an output directory
outdir="$HOME/bud_rnaseq/reads/raw/"
mkdir -p ${outdir}

#Combining the files, ${line} is the variable we are running through and set up the if command
i=1
while read line
do
    echo "Starting program"
    #else
    test ${i} -eq 1 && i=$(expr ${i} + 1) && continue

    #Getting the file name
    file_name=$(echo ${line} | cut -d "," -f 3)

    #Create the basis for the file name
    sample_name=$(echo ${line} | cut -d "," -f 8)

    #Check to see what sample is currently running
    echo ${sample_name}

    #Go to the directory with the files in
    cd $HOME/${file_name}

    #Separate R1 and R2
    read1=$(ls *R1*)
    read2=$(ls *R2*)

    #Check R1 and R2 correspond to the same section
    echo ${read1}
    echo ${read2}

    #Combine the data
    zcat ${read1} | gzip -c > ${outdir}/${sample_name}_R1.fq.gz
    zcat ${read2} | gzip -c > ${outdir}/${sample_name}_R2.fq.gz

done < bud_rnaseq/sample_information.csv
