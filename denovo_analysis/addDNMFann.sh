#!/bin/bash

fam=$1
memb='p1' #$2 #p1 or s1
suffix='-ann' #3 #-ann -ann-dnmfp1

vcfannotate -b dnmfilt/bed/${memb}/${fam}-${memb}-DNMFilter.bed -k DNMFilt${memb} -d NA ${fam}-HC-pm50${suffix}.vcf.gz | bgzip -c > ${fam}-HC-pm50${suffix}-dnmf${memb}.vcf.gz ; tabix -p vcf ${fam}-HC-pm50${suffix}-dnmf${memb}.vcf.gz &
vcfannotate -b dnmfilt/bed/${memb}/${fam}-${memb}-DNMFilter.bed -k DNMFilt${memb} -d NA ${fam}-JHC-pm50${suffix}.vcf.gz | bgzip -c  > ${fam}-JHC-pm50${suffix}-dnmf${memb}.vcf.gz ; tabix -p vcf ${fam}-JHC-pm50${suffix}-dnmf${memb}.vcf.gz &
vcfannotate -b dnmfilt/bed/${memb}/${fam}-${memb}-DNMFilter.bed -k DNMFilt${memb} -d NA ${fam}-FB-pm50${suffix}.vcf.gz | bgzip -c  > ${fam}-FB-pm50${suffix}-dnmf${memb}.vcf.gz; tabix -p vcf ${fam}-FB-pm50${suffix}-dnmf${memb}.vcf.gz &
vcfannotate -b dnmfilt/bed/${memb}/${fam}-${memb}-DNMFilter.bed -k DNMFilt${memb} -d NA ${fam}-PL-pm50${suffix}.vcf.gz | bgzip -c  > ${fam}-PL-pm50${suffix}-dnmf${memb}.vcf.gz; tabix -p vcf ${fam}-PL-pm50${suffix}-dnmf${memb}.vcf.gz &
