#!/bin/bash

infile=$1
outfile=$2

snpsift=/mnt/xfs1/bioinfo/software/installs/bcbio_nextgen/150607/Cellar/snpeff/4.1g/libexec/SnpSift.jar
snpeff=/mnt/xfs1/bioinfo/software/installs/bcbio_nextgen/150607/Cellar/snpeff/4.1g/libexec/snpEff.jar
snpeffconfig=/mnt/ceph/asalomatov/data/b37/snpEff/snpEff.config


java -Xmx5G -jar $snpsift extractFields -s "," -e "." -c $snpeffconfig -noLog $infile \
    CHROM POS ID REF ALT QUAL VARTYPE DNMFilt \
    "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].FEATURE" \
    "GEN[*].GT[*]" "GEN[*].AD[*]"  "GEN[*].DP[*]" "GEN[*].GQ[*]" "GEN[*].PL[*]" \
    > $outfile
