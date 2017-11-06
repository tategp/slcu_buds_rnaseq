#02_combining_data.sh - Combining the 4 files for each sample

#Change directory to where we want to be
cd //slcu.cam.ac.uk/data/TeamOL/personal_folders/gregt

#Getting the read counts, ${line} is the variable we are running through and set up the if command
i=1
while read line
do
    if [[ ${i} -eq 1 ]]
    then
      i=0
    else

    #Getting the file name
    #Getting a variable of just the column needed
    file_path=$(echo ${line} | cut -d "," -f 3)

    #Going to the corresponding file
    file_name="//slcu.cam.ac.uk/data/TeamOL/${file_path}"

    #Check to see how far the program is through running
    base_name="$(basename ${file_name})"
    echo ${base_name}

    #Create sample names
    sample_name_R1="${base_name}_R1"
    sample_name_R2="${base_name}_R2"

    #Go to the directory with the files in
    cd ${file_name}

    #Seperate R1 and R2
    read1=$(ls *R1*)
    read2=$(ls *R2*)

    #Check R1 and R2 correspond to the same section
    echo ${read1}
    echo ${read2}

    #Combine the data
    zcat ${read1} | gzip -c > ${sample_name_R1}.fq.gz
    zcat ${read2} | gzip -c > ${sample_name_R2}.fq.gz
  fi
done < sample_information.csv
