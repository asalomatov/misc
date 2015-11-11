#!/bin/bash

inpdir=/mnt/ceph/asalomatov/data/SSCexome/rerun200fam
snpsift='java -jar /mnt/xfs1/bioinfo/software/installs/bcbio_nextgen/150607/Cellar/snpeff/4.1g/libexec/SnpSift.jar'
fam=$1
for m in $(echo HC JHC FB)
do
    zcat ${inpdir}/${fam}-${m}-pm50-ann.vcf.gz | \
        $snpsift filter -c /mnt/ceph/asalomatov/data/b37/snpEff/snpEff.config \
        " ( GEN[0].GT = '0/0' ) | ( GEN[0].GT = '0|0' ) " | \
        $snpsift filter -c /mnt/ceph/asalomatov/data/b37/snpEff/snpEff.config \
        " ( GEN[1].GT = '0/0' ) | ( GEN[1].GT = '0|0' ) " | \
        $snpsift filter -c /mnt/ceph/asalomatov/data/b37/snpEff/snpEff.config \
        " ( GEN[2].GT != '0/0' ) & ( GEN[2].GT != '0|0' ) & ( GEN[2].GT != './.' ) " | \
        $snpsift filter -c /mnt/ceph/asalomatov/data/b37/snpEff/snpEff.config \
        " ( GEN[0].DP > 5 ) & ( GEN[1].DP > 5 )  & ( GEN[2].DP > 5 ) " > ${fam}-${m}-pm50-ann-p1.vcf
#        bgzip -c > ${fam}-${m}-pm50-ann-p1.vcf.gz
#        tabix -p ${fam}-${m}-pm50-ann-p1.vcf.gz
done

m='PL'
zcat ${inpdir}/${fam}-${m}-pm50-ann.vcf.gz | \
    $snpsift filter -c /mnt/ceph/asalomatov/data/b37/snpEff/snpEff.config \
    " ( GEN[0].GT = '0/0' ) | ( GEN[0].GT = '0|0' ) " | \
    $snpsift filter -c /mnt/ceph/asalomatov/data/b37/snpEff/snpEff.config \
    " ( GEN[1].GT = '0/0' ) | ( GEN[1].GT = '0|0' ) " | \
    $snpsift filter -c /mnt/ceph/asalomatov/data/b37/snpEff/snpEff.config \
    " ( GEN[2].GT != '0/0' ) & ( GEN[2].GT != '0|0' ) & ( GEN[2].GT != './.' ) " | \
    $snpsift filter -c /mnt/ceph/asalomatov/data/b37/snpEff/snpEff.config \
    " ( GEN[0].NR > 5 ) & ( GEN[1].NR > 5 )  & ( GEN[2].NR > 5 ) " > ${fam}-${m}-pm50-ann-p1.vcf
