Instructions on the order in which to run the programs.

The raw data can be found in /home/gt346/sequence_data/RNAseq/20170713_ssd_buds_nitrate_plasticity

01_combining_data.sh
This program takes the raw data files, which are split up into 4 files per sample and combines them.

02_filter_data.sh
This program applies FastQC, cutadapt and multiqc to the newly created files. Firstly it applies FastQC which gives data on the quality scores, GC content and many other factors. It then uses cutadapt to filter the data for minimum length, maximum number of missing values, minimum quality and removes adapters. After this the program then applies FastQC on the filtered reads and uses multiqc to compile the results into graphs and tables.

03_quantify_genes_raw.sh
This program runs through the .csv file taking out information to then submit 04_quantify_genes_raw_pipeline.sh for each sample to the scheduler.

03_quantify_genes_raw_pipeline.sh
Pipeline to apply Salmon-SMEM to the filtered reads to obtain counts of the transcripts. The program includes options to remove GC bias, sequencing bias, position bias and accounts for the strandedness of the reads.
