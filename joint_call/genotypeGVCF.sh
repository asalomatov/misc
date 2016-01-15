#!/bin/bash

java -Xms750m -Xmx3500m -XX:+UseSerialGC -Djava.io.tmpdir=/tmp/asalomatov/mktest -jar \
    /bioinfo/software/builds/GATK/fromGitZip150607/gatk-protected-master/target/executable/GenomeAnalysisTK.jar \
    -T GenotypeGVCFs \
    --variant  \
/tmp/asalomatov_working_11006_qrCOunNCWv/11006.mo.5284.merged.sorted.nodups.realigned.all_reads-re-fxgr-flr-dp-23-rlgn-rclb-30-bin.g.vcf
--variant
/tmp/asalomatov_working_11006_qrCOunNCWv/11006.p1.5316.merged.sorted.nodups.realigned.all_reads-re-fxgr-flr-dp-23-rlgn-rclb-30-bin.g.vcf
--variant
/tmp/asalomatov_working_11006_qrCOunNCWv/11006.fa.5252.merged.sorted.nodups.realigned.all_reads-re-fxgr-flr-dp-23-rlgn-rclb-30-bin.g.vcf
-L /tmp/asalomatov_working_11006_qrCOunNCWv/30__bin__11006-uni-mrg.bed -o
/tmp/asalomatov_working_11006_qrCOunNCWv/11006-JHC-30-bin.vcf -R
/bioinfo/data/bcbio_nextgen/150607/genomes/Hsapiens/GRCh37/seq/GRCh37.fa --dbsnp
/bioinfo/data/bcbio_nextgen/150607/genomes/Hsapiens/GRCh37/variation/dbsnp_138.vcf.gz
--standard_min_confidence_threshold_for_calling 10.0   --standard_min_confidence_threshold_for_emitting 10.0  

