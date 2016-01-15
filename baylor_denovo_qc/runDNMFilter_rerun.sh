#!/bin/bash

fam=$1
work_dir=/mnt/scratch/asalomatov/baylor_qc/dnm_filter
output_dir=./

###
dnm_dir=/mnt/scratch/asalomatov/software/src/DNMFilter-0.1.0
dnm_conf=/mnt/scratch/asalomatov/software/src/DNMFilter-0.1.0/Features.conf
train_set=/mnt/scratch/asalomatov/data/SSC/SSCdeNovoCalls/dnm_filter_train_sets/combined.csv
ref_gen=/mnt/xfs1/bioinfo/data/bcbio_nextgen/150607/genomes/Hsapiens/GRCh37/seq/GRCh37.fa
###

java -jar ${dnm_dir}/DNMFilter.jar gbm            \
    --reference     $ref_gen                    \
    --pedigree      ${work_dir}/${fam}.ped        \
    --bam           ${work_dir}/baylor-bams.csv        \
    --training      $train_set                      \
    --candidate     ${work_dir}/dnm_files/${fam}-test.csv   \
    --configuration $dnm_conf                               \
    --cutoff        0.0                                     \
    --output        ${ouput_dir}/output/${fam}-dnm.csv

