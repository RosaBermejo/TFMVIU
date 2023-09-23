#!/bin/bash
dir=$1

cd $dir
echo "sample-id,absolute-filepath,direction" > qiime_manifest.csv
for muestra in *R1_001.fastq.gz;
do
# temp=${muestra#*-*}
    samp_name=${muestra%%_*}
    submuestra1=${muestra%*_R1_*}
    submuestra2=${muestra#*_R1_*}

    echo -e "${samp_name}\t${submuestra1}_R1_${submuestra2}\t${submuestra1}_R2_${submuestra2}" >> samples_fastq_table.tsv
    echo "${samp_name},${PWD}/fastp/${samp_name}_R1_qc.fastq.gz,forward" >> qiime_manifest.csv
    echo "${samp_name},${PWD}/fastp/${samp_name}_R2_qc.fastq.gz,reverse" >> qiime_manifest.csv
done
