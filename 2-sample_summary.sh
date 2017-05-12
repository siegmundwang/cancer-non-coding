#!/bin/bash

#/////////////////////
# Author: SiegmundWANG
#/////////////////////

#/////////////////////////////////////////////////////////////
# Notification:
# 1. you should always fix in mind the necessity of "uniq"
# 2.
#/////////////////////////////////////////////////////////////

#/////////////////////////////////
# This script:
# 1. generate sample summary file.
# 2. 
#////////////////////////////////

dir="../data"

#//////////////////
# sample-project
#///////////////////
zcat ${dir}raw_data.tsv.gz|awk 'BEGIN{FS="\t";OFS="\t"}$14=="single base substitution"{print $2,$3}'|\
uniq > ${dir}sample-project.tsv

#//////////////////////////////////////
# count mutation freqency of all sample
#//////////////////////////////////////
awk 'BEGIN{FS="\t";OFS="\t"}{num[$4]+=1}END{for(idx in num)print idx,num[idx]}' ${dir}/sorted_sample_range.bed > \
${dir}/sample-mutation.tsv

#//////////////////////////////////////////////
# count freqency of CDS of all filtered1 sample
#//////////////////////////////////////////////
bedtools intersect -a ${dir}/sorted_sample_range.bed -b ${dir}/Hg19_cds.bed |sort|uniq|\
awk 'BEGIN{FS="\t";OFS="\t"}{num[$4]+=1}END{for(idx in num)print idx,num[idx]}' > ${dir}/sample-cds.tsv


#/////////////////////////////////////////////
# summary mutation rate and percentage of CDS 
#/////////////////////////////////////////////
awk 'BEGIN{FS="\t";OFS="\t"}{if(NR==FNR)num[$1]=$2;else{print $1,num[$1],$2,$2/num[$1],num[$1]/3095693981;}}' \
${dir}/sample-mutation.tsv ${dir}/sample-cds.tsv > ${dir}/sample-summary.tsv

#///////////////////////////////////////////////////////////////////////////////////////////
# make final table: sample-id, project, tissue, mutation number, mutation cds, mutaion rate
#///////////////////////////////////////////////////////////////////////////////////////////
awk 'BEGIN{FS="\t";OFS="\t"}{if(NR==FNR)project[$1]=$2;else{print $1,project[$1],$2,$3,$4,$5}}' ${dir}/sample-project.tsv \
${dir}/sample-summary.tsv ${dir}/prj-summ.tsv
awk 'BEGIN{FS="\t";OFS="\t"}{if(NR==FNR)tissue[$1]=$2;else{print $1,tissue[$2],$2,$3,$4,$5,$6}}' ${dir}/tissue_annotation.txt \
${dir}/prj-summ.tsv > ${dir}/sample-final-summary.tsv



