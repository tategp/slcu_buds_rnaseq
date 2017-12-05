#!/bin/bash

######## Calling options from the input ###########
while getopts “o:p:m:r:c:” OPTION
do
  case $OPTION in
    o)  outdir=$OPTARG;;
    p)  sample_name=$OPTARG;;
    m)  map=$OPTARG;;
    r)  ref=$OPTARG;;
    c)  threads=$OPTARG;;
    ?)  usage; exit;;
  esac
done

###### Running over the two types of aligned read

for file in ${map}/*.something.bam ##### Need to change the files running over
do
  ###### Shuffiling order of reads ######
  samtools collate ${file} ${outdir}/then however we want to lable file

  ####### Applying Salmon SMEM quantifcation ########
  salmon quant -t ${ref} \
  -l  ISR \
  -a Output from above \
  --output ${outdir} \
  -- seqBias --gcBias --posBias \
  --threads ${threads}
done
