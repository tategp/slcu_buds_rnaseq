#01_read_counts.sh - Count the number of reads for each sample in each library and output them as a csv file with sample name and count number so that we can then plot them.

#Change directory to where we want to be
cd //slcu.cam.ac.uk/data/TeamOL/personal_folders/gregt

#Getting the read counts, ${line} is the variable we are running through and set up the if command
i=1
while read line
do
    if [[ ${i} -eq 1 ]]
    then
      i=0
      echo "sample_name,file_name,count" > read_counts.csv
    else

    #Getting the sample name
    sample=$(echo ${line} | cut -d "," -f 8)
    sample_name="${sample}"

    #Getting the file name
    #Getting a variable of just the column needed
    file_path=$(echo ${line} | cut -d "," -f 3)

    #Going to the corresponding file
    file_name="//slcu.cam.ac.uk/data/TeamOL/${file_path}"

    #Getting the count number
    #Just considering the R1 files
    fastq_files=$(ls ${file_name}/*R1*)

    #Looping over these files, ${f} is the variable we are running through
    for f in $fastq_files
    do
      file=$(basename ${f})
      echo "Processing ${file}"
      count4=$(zcat ${f} | wc -l )
      count=$(expr $count4 / 4)
      #Output
      printf "${sample_name},${file},${count}\n" >> read_counts.csv
    done
  fi
done < sample_information.csv
