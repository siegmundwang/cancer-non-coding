#!/bin/bash

# Author: SiegmundWANG

# This script:                                           
# 1. get CDS region and covert *.gtf to *.bed            
# 2. summarize the SNV percentage in CDS region          
# 3. remove all the abnormal samples(donors) in raw data.


#/////////////////////////////////////////////////////////////
# Notification:
# 1. you should figure out "R_raw.bed" and "for_R_analysis.bed"
# 2. raw raw data are: "sample_range.bed"
#//////////////////////////////////////////////////////////////

dir="../data"

#/////////////////////////////////////////////////
# extract useful regions in Hg19 annotation file
# 1. CDS region: each extend by 5
# 2. IG|TR region: each extend by 50e3
#/////////////////////////////////////////////////
awk 'BEGIN{FS="\t";OFS="\t"}$3=="CDS"{print $1,$4-1,$5}' ${dir}/Hg19_ant.gtf > ${dir}/cds0.bed
bedtools slop -i ${dir}/cds0.bed -g ${dir}/hg19.genome -b 5 |sort -k1,1 -k2,2n > ${dir}/cds1.bed
bedtools merge -i ${dir}/cds1.bed > ${dir}/Hg19_cds.bed

awk 'BEGIN{FS="\t";OFS="\t"}$2~/IG|TR/{print $1,$4-1,$5}' ${dir}/Hg19_ant.gtf > ${dir}/igtr0.bed
bedtools slop -i ${dir}/igtr0.bed -g ${dir}/hg19.genome -b 50e3 |sort -k1,1 -k2,2n > ${dir}/igtr1.bed
bedtools merge -i ${dir}/igtr1.bed > ${dir}/Hg19_igtr.bed
#/////////////////////////////////////////////
# 1. sort the raw data for further analaysis
# 2. you shoule know the importance of 'uniq'
#/////////////////////////////////////////////

sort -k1,1 -k2,2n ${dir}/sample_range.bed |uniq > ${dir}/sorted_sample_range.bed

#//////////////////////////////////////////////////////////////////////////////////////
# 1. We should re-calculate mutation rate, one sample should not more than 300000 SNV!
# 2. length of the total human genome = 3095693981
# 3. threshold value = 300000/3095693981 = 9.69088e-05
#//////////////////////////////////////////////////////////////////////////////////////

awk 'BEGIN{FS="\t";OFS="\t"}{donr[$4]+=1;}END{for(x in donr)print x,donr[x]/3095693981;}' \
${dir}/sorted_sample_range.bed > ${dir}/sample_prop.tsv
awk 'BEGIN{FS="\t"}$2<=9.69088e-05{print $1}' ${dir}/sample_prop.tsv > ${dir}/P_avl_donor.tsv

#///////////////////////////////////////////
# filter hypermutated sample mutation data
#//////////////////////////////////////////
awk 'BEGIN{FS="\t";OFS="\t";}{if(NR==FNR){sample[$1]=0;}else{if($4 in sample)print $0;}}' ${dir}/P_avl_donor.tsv \
${dir}/sorted_sample_range.bed > ${dir}/filtered1_sample_range.bed

# count Freq of all filtered1 sample
awk 'BEGIN{FS="\t";OFS="\t"}{num[$4]+=1}END{for(idx in num)print idx,num[idx]}' ${dir}/filtered1_sample_range.bed > \
${dir}/Freq_donor.tsv

# count Freq of CDS of all filtered1 sample
bedtools intersect -a ${dir}/filtered1_sample_range.bed -b ${dir}/Hg19_cds.bed |sort|uniq|\
awk 'BEGIN{FS="\t";OFS="\t"}{num[$4]+=1}END{for(idx in num)print idx,num[idx]}' > ${dir}/Freq_CDS_donor.tsv


# get Freq available data
awk 'BEGIN{FS="\t";OFS="\t"}{if(NR==FNR)num[$1]=$2;else{if((num[$1]/$2)<=0.1)print $1,num[$1]/$2}}' \
${dir}/Freq_CDS_donor.tsv ${dir}/Freq_donor.tsv > ${dir}/P_Freq_avl_donor.tsv

#//////////////////////
# filter by CDS data
#//////////////////////

awk 'BEGIN{FS="\t";OFS="\t";}{if(NR==FNR){sample[$1]=0;}else{if($4 in sample)print $0;}}' ${dir}/P_Freq_avl_donor.tsv \
${dir}/filtered1_sample_range.bed > ${dir}/filtered2_sample_range.bed

#///////////////////////////////////////
# filter by intransic property
# 1. remove all mutation in CDS
# 2. remove all mutation in IG|TR region
#///////////////////////////////////////
bedtools intersect -a ${dir}/filtered2_sample_range.bed -b ${dir}/Hg19_cds.bed -v |sort|uniq|\
bedtools intersect -a stdin -b ${dir}/Hg19_igtr.bed -v |sort -k1,1 -k2,2n > ${dir}/filtered3_sample_range.bed


#//////////////////////////////////////////////////////////////////////////////
# merge samples to cluster, one cluster should be more than two mutation.
#//////////////////////////////////////////////////////////////////////////////
bedtools merge -i ${dir}/filtered3_sample_range.bed -d 50 -c 4 -o count,count_distinct,distinct |sort|uniq|\
awk 'BEGIN{FS="\t";OFS="\t";}$5>2{print $0}' > ${dir}/for_R_analysis.bed

#///////////////////////////////////////
# make row data only three cols
#///////////////////////////////////////
cut -f1-3 ${dir}/filtered3_sample_range.bed > ${dir}/R_raw.bed
