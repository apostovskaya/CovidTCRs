#!/bin/sh

# Analysis of targeted TCR libraries
# The 'analyze amplicon' pipeline includes alignment of raw sequencing reads using 'align',
# assembly of aligned sequences into clonotypes using assemble
# and exporting the resulting clonotypes into tab-delimited file using 'export'.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# BEFORE RUNNING:
# reminder to chmod u+rwx ./MiXCR_batch.sh
# reminder to check file names before running (modifications before R#.fastq)
# reminder to make sure there are no such folders in the ../ to not overwrite
# https://mixcr.readthedocs.io/en/master/analyze.html#ref-analyze
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# in-house generated 'split' dataset (IMSEQ study)
cd /Users/apost/Documents/CloudMail/PhD_2020/IMSEQ/tmp_wd

mkdir ../MiXCR
mkdir ../MiXCR/MiXCR_reports
mkdir ../MiXCR/MiXCR_outputs
for seqfile in $(find ../Minnn_output/6_Fasta -type f -name '*.fastq.gz' -exec basename "{}" \; \
  | cut -d_ -f1-4 | sort --unique); do
    input_file1=../Minnn_output/6_Fasta/"${seqfile}"_10umi_filter_R1.fastq.gz
    input_file2=../Minnn_output/6_Fasta/"${seqfile}"_10umi_filter_R2.fastq.gz
    mixcr analyze amplicon \
        --species hsa \
        --starting-material rna \
        --5-end no-v-primers --3-end c-primers \
        --adapters no-adapters \
        --report ../MiXCR/MiXCR_reports/"${seqfile}"_10umi_filter_mixcrReport \
        --receptor-type tcr \
        --threads 8 \
        --verbose\
        "${input_file1}" "${input_file2}" ../MiXCR/MiXCR_outputs/"${seqfile}"_10umi_filter_mixcrOut
done


# public 'mixed' dataset (IR-Binder, by Schultheiß et al. 2020)
cd /Users/apost/Documents/CloudMail/PhD_2020/IMSEQ/IR-Binder_PRJEB38339
mkdir ./MiXCR  # the same for MiXCR_new - added nt column to be compatible with VDJtools
mkdir ./MiXCR/MiXCR_reports
mkdir ./MiXCR/MiXCR_outputs
# delete IGH files
for d in ../IR-Binder_Fastq/*/; do
  find "$d" -type f -name 'IGH*' -delete
done
# delete some other files (although TRB, but not typical patients or HD)
for d in ../IR-Binder_Fastq/*/; do
  find "$d" -type f -name 'ChristophS*' -delete
done
# one more
for d in ../IR-Binder_Fastq/*/; do
  find "$d" -type f -name 'CMW*' -delete
done
# delete empty folders (after files were removed from them)
find ../IR-Binder_Fastq/ -empty -type d -delete
# The find command is used to search for files/directories matching a particular search criteria from the specified path,
# in this case the current directory (hence the .). The -empty option holds true for any file and directory that is empty.
#The -type d option holds true for the file type specified; in this case d stands for the file type directory.
#The -delete option is the action to perform, and holds true for all files found in the search.

# Most of the files kept have name format:
# TRB-HD01_S102_L001_R1_001.fastq.gz; TRB-Pt-9-2-250ng-15-04-2020-gDNA_S21_L001_R1_001.fastq.gz
# Some - TCRb-HD38-PB-gDNA_S58_L001_R1_001.fastq.gz; HD22-Jan2017-TCRb_S44_L001_R1_001.fastq.gz

long_arg="-count -fraction -nFeature CDR3 -aaFeature CDR3 -vHit -dHit -jHit -vGene -dGene -jGene"\
" -vAlignment -dAlignment -jAlignment -vHitScore -dHitScore -jHitScore"\
" -vBestIdentityPercent -dBestIdentityPercent -jBestIdentityPercent"
for d in ../IR-Binder_Fastq/*/; do
  for seqfile in $(find "$d" -type f -name '*.fastq.gz' -exec basename "{}" \; \
  | cut -d_ -f1-2 | sort --unique); do
    input_file1="$d""${seqfile}"_L001_R1_001.fastq.gz
    input_file2="$d""${seqfile}"_L001_R2_001.fastq.gz
    # to list all file names:
    # echo "$seqfile"
    mixcr analyze amplicon \
        --species hsa \
        --starting-material dna \
        --5-end v-primers --3-end j-primers \
        --receptor-type trb \
        --adapters adapters-present \
        --export "$long_arg" \
        --report ./MiXCR/MiXCR_reports/"${seqfile}"_mixcrReport \
        --threads 8 \
        --verbose\
        "${input_file1}" "${input_file2}" ./MiXCR/MiXCR_outputs/"${seqfile}"_mixcrOut
  done
done

