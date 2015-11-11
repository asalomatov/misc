#!/bin/bash

infile=$1
outfile=$2


snpsift=/mnt/xfs1/bioinfo/software/installs/bcbio_nextgen/150607/Cellar/snpeff/4.1g/libexec/SnpSift.jar
snpeff=/mnt/xfs1/bioinfo/software/installs/bcbio_nextgen/150607/Cellar/snpeff/4.1g/libexec/snpEff.jar
snpeffconfig=/mnt/ceph/asalomatov/data/b37/snpEff/snpEff.config


java -Xmx5G -jar $snpsift annotate -id -c $snpeffconfig -noLog -dbsnp $infile | \
    java -Xmx5G -jar $snpsift varType -c $snpeffconfig -noLog - | \
    java -Xmx5G -jar $snpeff ann -noLog -c $snpeffconfig GRCh37.75 -csvStats -s ${outfile}.summary.csv | \
    bgzip -c > $outfil
tabix -p vcf $outfile
