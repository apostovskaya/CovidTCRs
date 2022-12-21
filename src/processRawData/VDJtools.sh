#!/bin/sh

# cd /Users/apost/Documents/CloudMail/PhD_2020/IMSEQ

# META_FILE1=$1 # metadata_TRB_10umi_filter_IMSEQ_MiXCR.txt

# if NOT installed with Homebrew
# set variables
GETJAVA="/Library/Java/JavaVirtualMachines/adoptopenjdk-8.jdk/Contents/Home/jre/bin/java"
VDJTOOLS="./vdjtools-1.2.1/vdjtools-1.2.1.jar"

# set path to the old java version and .jar file for execution, zsh version, sth didn't work
#alias -g VDJTOOLS="/Library/Java/JavaVirtualMachines/adoptopenjdk-8.jdk/Contents/Home/jre/bin/java -jar -Xmx20G ./vdjtools-1.2.1/vdjtools-1.2.1.jar"

# tests
# $GETJAVA -jar -Xmx20G $VDJTOOLS Convert -S mixcr ./MiXCR_outputs/10_cd4_3_S49_mixcr.clonotypes.TRB.txt ./VDJtools_format
# $GETJAVA -jar -Xmx20G $VDJTOOLS Convert -S mixcr -m ./metadata_test_IMSEQ_MiXCR.txt ./VDJtools_format/test

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# for TRB and TRA, substitute manually
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mkdir ./VDJtools_format_TRA
# convert mixcr files to vdjtools format, I manually renamed output metadata file after (IMSEQ)
$GETJAVA -jar -Xmx20G $VDJTOOLS Convert -S mixcr -m ./Metadata/metadata_TRA_10umi_filter_IMSEQ_MiXCR.txt ./VDJtools_format_TRA
# for IR-Binder data
$GETJAVA -jar -Xmx20G $VDJTOOLS Convert -S mixcr -m ./MiXCR/metadata_TRB_IR_Binder_MiXCR.txt ./VDJtools_format_TRB

# calculate basic statistics -> file
$GETJAVA -jar -Xmx20G $VDJTOOLS CalcBasicStats -m ./metadata_TRA_10umi_filter_MiXCR_to_VDJtools.txt \
./VDJtools_TRA_stats/IMSEQ_TRA

# re-install R to produce plots, didn't fully finish though, interrupted it as it got stuck, but worked
#$GETJAVA -jar -Xmx20G $VDJTOOLS Rinstall

# segment usage
# actual usage calculation
mkdir ./VDJtools_Segment_Usage_TRA
$GETJAVA -jar -Xmx20G $VDJTOOLS CalcSegmentUsage -p -m ./metadata_TRA_10umi_filter_MiXCR_to_VDJtools.txt \
./VDJtools_Segment_Usage_TRA/IMSEQ_TRA
$GETJAVA -jar -Xmx20G $VDJTOOLS CalcSegmentUsage -p -n -f "patient_id" -l "s_id" \
-m ./metadata_TRB_10umi_filter_MiXCR_to_VDJtools.txt ./VDJtools_Segment_Usage_TRB/IMSEQ_TRB_by_patient
#$GETJAVA -jar -Xmx20G $VDJTOOLS CalcSegmentUsage -p -n -f "T-cell_type" \
#-m ./metadata_TRB_10umi_filter_MiXCR_to_VDJtools.txt ./VDJtools_Segment_Usage_TRB/IMSEQ_TRB_by_cell_type
# Error in seq.int(rx[1L], rx[2L], length.out = nb) :
#  'from' must be a finite number
#Calls: cut -> cut.default
#In addition: Warning messages:
#1: NAs introduced by coercion
#2: In min(x) : no non-missing arguments to min; returning Inf
#3: In max(x) : no non-missing arguments to max; returning -Inf
#Execution halted

# diversity measures
$GETJAVA -jar -Xmx20G $VDJTOOLS CalcDiversityStats -m ./metadata_TRA_10umi_filter_MiXCR_to_VDJtools.txt ./IMSEQ_TRA

# rarefaction plots
$GETJAVA -jar -Xmx20G $VDJTOOLS RarefactionPlot -m ./metadata_TRA_10umi_filter_MiXCR_to_VDJtools.txt \
-f "patient_id" -l "patient_id" --wide-plot --label-exact ./IMSEQ_TRA
# plot manually in RStudio

# Rscript rarefaction_curve.r ./IMSEQ_VDJtools.rarefaction.strict.txt 1 1 F F F F ./IMSEQ_VDJtools.rarefaction.strict.pdf
# open R in terminal, then run
# install.packages("ggplot2", repos="http://cran.r-project.org", lib="/Users/apost/Documents/CloudMail/PhD_2020/IMSEQ/vdjtools-1.2.1/Rpackages/")
# wget http://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_2.0.0.tar.gz

$GETJAVA -jar -Xmx20G $VDJTOOLS CalcSpectratype -m ./metadata_TRA_10umi_filter_MiXCR_to_VDJtools.txt ./IMSEQ_TRA
# $GETJAVA -jar -Xmx20G $VDJTOOLS PlotFancySpectratype -m ./metadata_MiXCR_to_VDJtools.txt ./IMSEQ_VDJtools

$GETJAVA -jar -Xmx20G $VDJTOOLS PlotFancySpectratype

