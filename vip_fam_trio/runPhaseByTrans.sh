#!/bin/bash
echo $1
echo $2
echo $3
echo $4
echo $5

java -jar /bioinfo/software/builds/GATK/fromGit140923/gatk-protected/protected/gatk-package-distribution/target/gatk-package-distribution-3.2.jar -R /bioinfo/data/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa -T PhaseByTransmission -V $1 -ped $2 -o $3 -et NO_ET --pedigreeValidationType SILENT --DeNovoPrior $4 --MendelianViolationsFile $5

