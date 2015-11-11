#!/bin/bash

inpdir=/mnt/scratch/asalomatov/SSC_rerun/denovo_analysis/rerun200fam
snpsift='java -jar /mnt/xfs1/bioinfo/software/installs/bcbio_nextgen/150607/Cellar/snpeff/4.1g/libexec/SnpSift.jar'
fam=$1
for m in $(echo HC JHC FB)
do
    zcat ${inpdir}/${fam}-${m}-pm50-ann.vcf.gz | \
        $snpsift filter -c /mnt/scratch/asalomatov/data/b37/snpEff/snpEff.config \
        " ( GEN[0].GT != '0/0' ) & ( GEN[0].GT != '0|0' ) & ( GEN[0].GT != './.' ) & ( GEN[0].GT != '.|.' )" | \
        $snpsift filter -c /mnt/scratch/asalomatov/data/b37/snpEff/snpEff.config \
        " ( GEN[1].GT != '0/0' ) & ( GEN[1].GT != '0|0' ) & ( GEN[1].GT != './.' ) & ( GEN[1].GT != '.|.' )" | \
        $snpsift filter -c /mnt/scratch/asalomatov/data/b37/snpEff/snpEff.config \
        " ( GEN[3].GT = '0/0' ) | ( GEN[3].GT = '0|0' ) " | \
        $snpsift filter -c /mnt/scratch/asalomatov/data/b37/snpEff/snpEff.config \
        " ( GEN[0].DP > 5 ) & ( GEN[1].DP > 5 )  & ( GEN[3].DP > 5 ) " > ${fam}-${m}-pm50-ann-s1.vcf
#        bgzip -c > ${fam}-${m}-pm50-ann-s1.vcf.gz
#        tabix -p ${fam}-${m}-pm50-ann-s1.vcf.gz
done

m='PL'
zcat ${inpdir}/${fam}-${m}-pm50-ann.vcf.gz | \
        $snpsift filter -c /mnt/scratch/asalomatov/data/b37/snpEff/snpEff.config \
        " ( GEN[0].GT != '0/0' ) & ( GEN[0].GT != '0|0' ) & ( GEN[0].GT != './.' ) & ( GEN[0].GT != '.|.' )" | \
        $snpsift filter -c /mnt/scratch/asalomatov/data/b37/snpEff/snpEff.config \
        " ( GEN[1].GT != '0/0' ) & ( GEN[1].GT != '0|0' ) & ( GEN[1].GT != './.' ) & ( GEN[1].GT != '.|.' )" | \
        $snpsift filter -c /mnt/scratch/asalomatov/data/b37/snpEff/snpEff.config \
        " ( GEN[3].GT = '0/0' ) | ( GEN[3].GT = '0|0' ) " | \
        $snpsift filter -c /mnt/scratch/asalomatov/data/b37/snpEff/snpEff.config \
        " ( GEN[0].NR > 5 ) & ( GEN[1].NR > 5 )  & ( GEN[3].NR > 5 ) " > ${fam}-${m}-pm50-ann-s1.vcf
