#!/bin/bash

fam=$1
cp ios_vcf_header.vcf ${fam}-ios.vcf
erep ^$fam SSCexome.csv | ./ios2vcf.py $GENOMEREF >> ${fam}-ios.vcf
bgzip ${fam}-ios.vcf
tabix -p vcf ${fam}-ios.vcf.gz
mkdir -p 200fam
mv ${fam}-ios.vcf* 200fam/
