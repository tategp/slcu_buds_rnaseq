#!/bin/bash

######## Calling options from the input ###########
while getopts “o:p:1:2:i:c:” OPTION
do
  case $OPTION in
    o)  outdir=$OPTARG;;
    1)  read1=$OPTARG;;
    2)  read2=$OPTARG;;
    i)  index_folder=$OPTARG;;
    c)  threads=$OPTARG;;
    ?)  usage; exit;;
  esac
done

####### Applying Salmon SMEM quantifcation ########

salmon quant -i ${index_folder} \
-l  ISR \
-1 <(gunzip -c ${read1}) -2 <(gunzip -c ${read2}) \
--output ${outdir} \
-- seqBias --gcBias --posBias \
--threads ${threads}
