#!/bin/bash

# 1.
cd ../IMSEQ/Runs

# 2. creating a folder for each sample and copying all relevant files there
# to extract unique file names prefixes, where -d_ specifies the delimeter (_) in the name, -f1-4 says to take fields 1-4
# (fields - sth in the name separated with _, so 1_CD4_1_S1)
for F in $(ls *.gz| cut -d_ -f1-4 | uniq); do
    echo Creating a directory for samples "${F}" 2>&1 | tee -a ./MoveLog.txt
    mkdir "${F}" && mv "${F}"*.gz "${F}"
    echo Files have been copied to "${F}" 2>&1 | tee -a ./MoveLog.txt
done

# 3.
mkdir ../Merged_runs

# 4. merging lanes into one file and putting into separate folder
for folder in */; do
  cat ./"${folder}"?*R1_001.fastq.gz > ../Merged_runs/"${folder%/}"_R1.fastq.gz
  echo "${folder%/}"_R1 files have been merged 2>&1 | tee -a ../Merged_runs/MergeLog.txt
  cat ./"${folder}"?*R2_001.fastq.gz > ../Merged_runs/"${folder%/}"_R2.fastq.gz
  echo "${folder%/}"_R2 files have been merged 2>&1 | tee -a ../Merged_runs/MergeLog.txt
done


# reading list of files from file 'Sample_ID_list.tsv' (for larger number of files)
#cat Sample_ID_list.tsv | while read LINE;do
#  echo "${LINE}";
#done
##### didn't use it

# for i in 1 2 3 4; do
#    for file in NA24694_GCCAAT_L001_R${i}_*fastq.gz; do
#        cat "$file" >> EA00694_GCCAAT_L001_R${i}.fastq.gz

