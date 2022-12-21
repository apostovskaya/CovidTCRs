#!/bin/sh

# Script to assemble consensuses by UMIs in in-house generated "split" dataset (IMSEQ study).
# We have 2 FASTQ files with R1 and R2.
# We want to assemble consensuses by UMI that is 12 first nucleotides in R1 file.
# There must be only 1 consensus for each UMI.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# https://minnn.readthedocs.io/routines.html
# Consensus assembly consists of 6 stages:
# 1. Extract barcodes (UMIs) from raw sequences.
# 2. Sort sequences by UMI values to group them for further correction.
# 3. Correct mismatches and indels in UMIs.
# 4. Sort sequences by UMI values to group them for further consensus assembly.
# 5. Assembly one consensuses for each UMI.
# 6. Export calculated consensuses to FASTQ format.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# BEFORE RUNNING:
# reminder to make sure there are no such folders in the ../ to not overwrite
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cd /Users/apost/Documents/CloudMail/PhD_2020/IMSEQ/Merged_runs
# at the end, emptied all folders, except 3 and 6 (too heavy to keep)
mkdir ../Minnn_output
mkdir ../Minnn_output/1_Minnn_extracted
mkdir ../Minnn_output/2_Minnn_sorted
mkdir ../Minnn_output/3_Minnn_corrected # compressed at the end
mkdir ../Minnn_output/4_Minnn_filtered
mkdir ../Minnn_output/5_Minnn_sorted-2
mkdir ../Minnn_output/6_Fasta # compressed at the end

for seqfile in $(find ../Merged_runs -type f -name '*.gz' -exec basename "{}" \; \
  | cut -d_ -f1-4 | sort --unique); do
    # R1 and R2 switched in order on purpose, because UMI is supposed to be in the beginning of the 2nd read
    # and Minnn looks for the pattern in the order files are provided for input
    input_file1=../Merged_runs/"${seqfile}"_R2.fastq.gz
    input_file2=../Merged_runs/"${seqfile}"_R1.fastq.gz
    echo "##########################################################"
    echo Starting Minnn for the "${seqfile}": # would be good to add progress bar like file n/x
    echo "##########################################################"
    echo Starting Minnn extraction for the "${seqfile}":
    echo "##########################################################"
    # We have input data where UMI is first 12 nucleotides of R1
    minnn extract --pattern "^(UMI:N{12})\*" -f --input "${input_file1}" "${input_file2}" \
    --output ../Minnn_output/1_Minnn_extracted/"${seqfile}"_extracted.mif \
    --report ../Minnn_output/1_Minnn_extracted/report_"${seqfile}".txt
    minnn sort --groups UMI --input ../Minnn_output/1_Minnn_extracted/"${seqfile}"_extracted.mif \
    --output ../Minnn_output/2_Minnn_sorted/"${seqfile}"_sorted.mif

    echo "##########################################################"
    echo Starting to correct errors in barcodes in "${seqfile}":
    echo "##########################################################"
    minnn correct --groups UMI --input ../Minnn_output/2_Minnn_sorted/"${seqfile}"_sorted.mif \
    --output ../Minnn_output/3_Minnn_corrected/"${seqfile}"_corrected.mif \
    --report ../Minnn_output/3_Minnn_corrected/report_"${seqfile}".txt

    echo "##############################################################################"
    echo Removing reads with UMIs that have a count of less than 10 from "${seqfile}":
    echo "##############################################################################"
    minnn filter-by-count --groups UMI --input ../Minnn_output/3_Minnn_corrected/"${seqfile}"_corrected.mif \
    --output ../Minnn_output/4_Minnn_filtered/"${seqfile}"_filtered.mif --min-count 10 \
    --report ../Minnn_output/4_Minnn_filtered/report_"${seqfile}".txt
    minnn sort --groups UMI --input ../Minnn_output/4_Minnn_filtered/"${seqfile}"_filtered.mif \
    --output ../Minnn_output/5_Minnn_sorted-2/"${seqfile}"_sorted-2.mif

    echo "##########################################################"
    echo Exporting "${seqfile}" files to fastq format
    echo "##########################################################"
    # R1 and R2 switched because that's how they were put in, because UMI is in R2
    minnn mif2fastq --copy-original-headers --input ../Minnn_output/5_Minnn_sorted-2/"${seqfile}"_sorted-2.mif \
    --group R1=../Minnn_output/6_Fasta/"${seqfile}"_10umi_filter_R2.fastq \
    --group R2=../Minnn_output/6_Fasta/"${seqfile}"_10umi_filter_R1.fastq

    echo "##########################################################"
    echo Compressing "${seqfile}"_correcetd.mif and both .fastq files:
    echo "##########################################################"
    gzip ../Minnn_output/3_Minnn_corrected/"${seqfile}"_corrected.mif
    gzip ../Minnn_output/6_Fasta/"${seqfile}"_10umi_filter_R2.fastq
    gzip ../Minnn_output/6_Fasta/"${seqfile}"_10umi_filter_R1.fastq
    # delete all other output files here
done

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# run manually for 15_cd8_3_S19, 19_cd8_3_S22, 5_cd8_2_S48; 9_cd8_3_S16?
#input_file1=../Merged_runs/5_cd8_2_S48_R2.fastq.gz # because UMI is supposed to be in the beginning of the 2nd read
#input_file2=../Merged_runs/5_cd8_2_S48_R1.fastq.gz
#minnn extract --pattern "^(UMI:N{12})\*" -f --input "${input_file1}" "${input_file2}" \
#    --output ../Minnn_output/1_Minnn_extracted/5_cd8_2_S48_extracted.mif
#minnn sort --groups UMI -f --input ../Minnn_output/1_Minnn_extracted/5_cd8_2_S48_extracted.mif \
#    --output ../Minnn_output/2_Minnn_sorted/5_cd8_2_S48_sorted.mif
#minnn correct --groups UMI -f --input ../Minnn_output/2_Minnn_sorted/5_cd8_2_S48_sorted.mif \
#    --output ../Minnn_output/3_Minnn_corrected/5_cd8_2_S48_corrected.mif \
#    --report ../Minnn_output/3_Minnn_corrected/report_5_cd8_2_S48.txt
#minnn filter-by-count --groups UMI -f --input ../Minnn_output/3_Minnn_corrected/5_cd8_2_S48_corrected.mif \
#    --output ../Minnn_output/4_Minnn_filtered/5_cd8_2_S48_filtered.mif --min-count 10 \
#    --report ../Minnn_output/4_Minnn_filtered/report_5_cd8_2_S48.txt
#minnn sort --groups UMI -f --input ../Minnn_output/4_Minnn_filtered/5_cd8_2_S48_filtered.mif \
#    --output ../Minnn_output/5_Minnn_sorted-2/5_cd8_2_S48_sorted-2.mif
#minnn mif2fastq --copy-original-headers -f --input ../Minnn_output/5_Minnn_sorted-2/5_cd8_2_S48_sorted-2.mif \
#    --group R1=../Minnn_output/6_Fasta/5_cd8_2_S48_10umi_filter_R2.fastq \
#    --group R2=../Minnn_output/6_Fasta/5_cd8_2_S48_10umi_filter_R1.fastq

#gzip ../Minnn_output/6_Fasta/5_cd8_2_S48_10umi_filter_R2.fastq
#gzip ../Minnn_output/6_Fasta/5_cd8_2_S48_10umi_filter_R1.fastq

# batch compress all fastq outputs in 6_Fasta folder into gz
#for seqfile in $(find ./Minnn_output/6_Fasta -type f -name '*.fastq' -exec basename "{}" \; \
#  | cut -d_ -f1-4 | sort --unique); do
#    gzip ./Minnn_output/6_Fasta/"${seqfile}"_10umi_filter_R2.fastq
#    gzip ./Minnn_output/6_Fasta/"${seqfile}"_10umi_filter_R1.fastq
#done

# the same for .mif files in 3_Minnn_corrected
#for seqfile in $(find ./Minnn_output/3_Minnn_corrected -type f -name '*.mif' -exec basename "{}" \; \
#  | cut -d_ -f1-4 | sort --unique); do
#    gzip ./Minnn_output/Minnn_corrected/"${seqfile}"_corrected.mif
#done