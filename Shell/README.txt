Instructions on the order in which to run the programs.

The raw data can be found in /home/gt346/sequence_data/RNAseq/20170713_ssd_buds_nitrate_plasticity

01_combining_data.sh
This program takes the raw data files, which are split up into 4 files per sample, and combines these to make 1 file per sample.

02_filter_data.sh
This program applies FastQC, cutadapt and multiqc to the newly created files. Firstly it applies FastQC which gives data on the quality scores, GC content and many other factors. Then it uses cutadapt to filter the data for minimum length, maximum number, minimum quality of missing values and removes adapters. After this the program then applies FastQC on the filtered reads and then uses multiqc to compile the results into graphs and tables and outputs them as an HTML file.

03_mapping.sh
This program runs through the .csv file taking out information to then submit 03_mapping_pipeline.sh for each sample to the scheduler.

03_mapping_pipeline.sh
Pipeline to run mapping of the data to the given reference. It aligns using STAR and HISAT2 and then marks duplicates and indexes using Picard Tools and creates output statistics using star flagstat. It also quantifies the gene count using STAR.

04_quantify_genes_raw.sh
This program runs through the .csv file taking out information to then submit 04_quantify_genes_raw_pipeline.sh for each sample to the scheduler.

04_quantify_genes_raw_pipeline.sh
Pipeline to apply Salmon-SMEM to the filtered reads to obtain counts of the transcripts.

04_quantify_genes_aligned.sh
This program runs through the .csv file taking out information to then submit 04_quantify_genes_aligned_pipeline.sh for each sample to the scheduler.

04_quantify_genes_aligned_pipeline.sh
Pipeline to apply Salmon-Aln to the previously aligned reads from STAR and HISAT2 to obtain counts of the transcripts.
