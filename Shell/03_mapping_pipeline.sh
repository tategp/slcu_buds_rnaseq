#!/bin/bash

######## Calling options from the input ###########
while getopts “o:p:1:2:r:a:i:c:” OPTION
do
  case $OPTION in
    o)  outdir=$OPTARG;;
    p)  sample_name=$OPTARG;;
    1)  read1=$OPTARG;;
    2)  read2=$OPTARG;;
    r)  ref=$OPTARG;;
    a)  annotation=$OPTARG;;
    i)  index_folder=$OPTARG;;
    c)  threads=$OPTARG;;
    ?)  usage; exit;;
  esac
done

############ Mapping using STAR ############

STAR \
--runMode "alignReads" \
--genomeDir "${index_folder}" \
--readFilesIn ${read1} ${read2} \
--readFilesCommand zcat \
--outFileNamePrefix ${outdir}/${sample_name}.star. \
--outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts \
--alignIntronMax 2000 \
--runThreadN ${threads}

mv ${outdir}/${sample_name}.star.Aligned.sortedByCoord.out.bam ${outdir}/${sample_name}.star.bam


############ Mapping using HISAT2 ############

hisat2 \
-x "${ref}" \
-1 ${read1} -2 ${read2} \
--rna-strandness RF \
--max-intronlen 2000 \
--n-ceil L,0,0.10 \
--threads ${threads} \
 | samtools sort -o ${outdir}/${sample_name}.histat2.bam


 #### Run through the created .bam files to mark duplicates and index using Picard tools, then create summary of statistics using samtools flagstat and quantify transcripts using Salmon #####q

for file in ${outdir}/*.bam
do
  file_path=$(echo ${file} | cut -f 1 -d '.')
  mapper=$(echo ${file} | cut -f 2 -d '.')
  java -jar /home/hugot/software/picard-tools/1.135/install/picard.jar MarkDuplicates \
      I=${file} \
      O=${file_path}.${mapper}.markdup.bam \
      M=${file_path}.${mapper}.marked_dup_metrics.txt \
      CREATE_INDEX=true

  samtools flagstat ${file_path}.${mapper}.markdup.bam > ${file_path}.${mapper}.markdup.flagstat.txt

done
