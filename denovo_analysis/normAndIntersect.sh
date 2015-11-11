#!/bin/bash

fam=$1
vt normalize -q -r $GENOMEREF ${fam}-ios.vcf.gz | vcfintersect -b ../b37.exome-pm50.bed | bgzip -c > ${fam}-ios-pm50.vcf.gz
tabix -p vcf ${fam}-ios-pm50.vcf.gz
