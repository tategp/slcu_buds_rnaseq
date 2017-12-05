#!/bin/bash

# bash /home/gt346/bud_rnaseq/scripts/Shell/04_quantify_genes_raw.sh

############# Setting up the relevant directories ############

# Working directocry
wd="/home/gt346/bud_rnaseq"
cd ${wd}

# Path to the filtered reads
read_path="/home/gt346/bud_rnaseq/reads/filtered/filtered_reads"

# Where the indexed genome will be stored
index_folder="/home/hugot/reference/arabidopsis/tair10/ensembl_release-37/annotation/salmon_index/fmd/"

# Create an output directory
outdir="/home/gt346/bud_rnaseq/gene_quantification"
mkdir -p ${outdir}


######### Applying the program 04_quantify_genes_raw_pipeline.sh to each of our samples #########

i=1
while read line
do
    # Skip the first line of the .csv file
    test ${i} -eq 1 && i=$(expr ${i} + 1) && continue

    # Want to retrive relavant information from the .csv file
    sample_name=$(echo ${line} | cut -d "," -f 8)
    read1=${read_path}/${sample_name}.filtered._R1.fq.gz
    read2=${read_path}/${sample_name}.filtered._R2.fq.gz

    # Make output directory
    mkdir -p ${outdir}/${sample_name}

    # Applying the program mapping_data.sh with the given inputs
    echo \
    "bash /home/gt346/bud_rnaseq/scripts/Shell/04_quantify_genes_raw_pipeline.sh \
    -o ${outdir}/${sample_name} \
    -1 ${read1} \
    -2 ${read2} \
    -i ${index_folder} \
    -c 3 " | \
    qsub -N ${sample_name}_quantify_genes -o ${outdir}/${sample_name}/${sample_name}.stdout.txt \
    -e ${outdir}/${sample_name}/${sample_name}.stderr.txt \
    -S /bin/bash -pe smp 3 -V \
    -wd "${wd}"

done < sample_information.csv
