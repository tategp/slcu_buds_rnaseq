#!/bin/bash

# bash /home/gt346/temp/scripts/Shell/04_mapping.sh

############# Setting up the relevant directories ############

# Working directocry
wd="/home/gt346/bud_rnaseq"
cd ${wd}

# Path to the filtered reads
map="/home/gt346/bud_rnaseq/mapping"

# Reference genome
ref="/home/hugot/reference/arabidopsis/tair10/ensembl_release-37/genome/Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa"

# Create an output directory
outdir="/home/gt346/bud_rnaseq/gene_quantification"


######### Applying the program 04_quantify_genes_pipeline.sh to each of our samples #########

i=1
while read line
do
    # Skip the first line of the .csv file
    test ${i} -eq 1 && i=$(expr ${i} + 1) && continue

    # Applying the program mapping_data.sh with the given inputs
    echo \
    "bash /home/gt346/bud_rnaseq/scripts/Shell/04_quantify_genes_pipeline.sh \
    -o ${outdir}/${sample_name} \
    -p ${sample_name} \
    -m ${map}/${sample_name} \
    -r ${ref}
    -c 3 " | \
    qsub -N ${sample_name}_map -o ${outdir}/${sample_name}/${sample_name}.stdout.txt \
    -e ${outdir}/${sample_name}/${sample_name}.stderr.txt \
    -S /bin/bash -pe smp 3 -V \
    -wd "${wd}"

done < sample_information.csv
