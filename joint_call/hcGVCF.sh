#!/bin/bash

java -Xms750m -Xmx5g -XX:+UseSerialGC -Djava.io.tmpdir=/tmp -jar \
    /bioinfo/software/builds/GATK/fromGitZip150607/gatk-protected-master/target/executable/GenomeAnalysisTK.jar \
    -T HaplotypeCaller \
    -I $1 \
    -o ${1}.g.vcf \
    -R /bioinfo/data/bcbio_nextgen/150607/genomes/Hsapiens/GRCh37/seq/GRCh37.fa \
    --dbsnp /bioinfo/data/bcbio_nextgen/150607/genomes/Hsapiens/GRCh37/variation/dbsnp_138.vcf.gz \
    -L /mnt/scratch/asalomatov/data/b37/b37.exome-pm100.bed \
    --emitRefConfidence GVCF \
    --standard_min_confidence_threshold_for_calling 30.0 \
    --standard_min_confidence_threshold_for_emitting 10.0 \
    --annotation BaseQualityRankSumTest \
    --annotation FisherStrand \
    --annotation GCContent \
    --annotation HaplotypeScore \
    --annotation HomopolymerRun \
    --annotation MappingQualityRankSumTest \
    --annotation MappingQualityZero \
    --annotation QualByDepth \
    --annotation ReadPosRankSumTest \
    --annotation RMSMappingQuality \
    --annotation DepthPerAlleleBySample \
    --annotation Coverage --interval_set_rule INTERSECTION --annotation ClippingRankSumTest \
    --annotation DepthPerSampleHC   \
    --pair_hmm_implementation VECTOR_LOGLESS_CACHING \
    -U LENIENT_VCF_PROCESSING  \
    --read_filter BadCigar   \
    --read_filter NotPrimaryAlignment \
    -nct 1

#bgzip ${1}.g.vcf
#tabix -p vcf ${1}.g.vcf.gz
