#!/bin/bash

fam=$1
vt normalize -q -r $GENOMEREF /mnt/ceph/asalomatov/SSC_Eichler/rerun/ssc${fam}/${fam}-HC-vars-flr.vcf.gz | vcfintersect -b ../b37.exome-pm50.bed | bgzip -c > ${fam}-HC-pm50.vcf.gz
vt normalize -q -r $GENOMEREF /mnt/ceph/asalomatov/SSC_Eichler/rerun/ssc${fam}/${fam}-JHC-vars.vcf.gz | vcfintersect -b ../b37.exome-pm50.bed | bgzip -c > ${fam}-JHC-pm50.vcf.gz
vt normalize -q -r $GENOMEREF /mnt/ceph/asalomatov/SSC_Eichler/rerun/ssc${fam}/${fam}-FB-vars.vcf.gz | vcfintersect -b ../b37.exome-pm50.bed | bgzip -c > ${fam}-FB-pm50.vcf.gz
vt normalize -q -r $GENOMEREF /mnt/ceph/asalomatov/SSC_Eichler/rerun/ssc${fam}/${fam}-PL-vars.vcf.gz | vcfintersect -b ../b37.exome-pm50.bed | bgzip -c > ${fam}-PL-pm50.vcf.gz
tabix -p vcf ${fam}-HC-pm50.vcf.gz
tabix -p vcf ${fam}-JHC-pm50.vcf.gz
tabix -p vcf ${fam}-FB-pm50.vcf.gz
tabix -p vcf ${fam}-PL-pm50.vcf.gz
