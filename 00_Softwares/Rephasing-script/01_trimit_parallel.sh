#!/bin/bash


helpFunction()
{
   echo ""
   echo "Usage: $(basename -- "$0") -I <PATH_TO_RAWReads_Dir> -O <PATH_TO_CleanReads_Dir>"
   echo ""
   echo -e "\t-I Path to Directory containing Raw Reads fastq files."
   echo -e "\t-O Path to Output Directory."
   echo ""
   exit 1 # Exit script after printing help
}

# Parse the parameters
while getopts "I:O:T" opt
do
   case "$opt" in
      I ) parameterI="$OPTARG" ;;
      O ) parameterO="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

mkdir -p ${parameterO}


RAW_READS=`realpath ${parameterI}`
CLEAN_READS=`realpath ${parameterO}`

ls ${RAW_READS} | grep "f1.fq.gz" | sed -e 's/_f1.fq.gz//g' >accessions.txt

cat accessions.txt | parallel -j 5 fastp -i ${RAW_READS}/{}_f1.fq.gz -I ${RAW_READS}/{}_r2.fq.gz -o ${CLEAN_READS}/{}_trim_f1.fq.gz -O ${CLEAN_READS}/{}_trim_r2.fq.gz -w 16
rm accessions.txt
