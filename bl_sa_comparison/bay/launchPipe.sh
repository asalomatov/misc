#!/bin/bash

protoc=$1

for i in $(echo 1 2 4)
do
    sbatch -J  ${i}_$protoc -o ${i}_${protoc}.out -e ${i}_${protoc}.err -N 1 --exclusive -D ./ /mnt/xfs1/home/asalomatov/projects/pipeline/ppln/pipe03.sh  /mnt/scratch/asalomatov/baylor_qc/inputs_$protoc ./$protoc $i WG 0 work /mnt/xfs1/home/asalomatov/projects/pipeline/ppln/include.mk 0  ,FixGroups,HaplotypeCaller,Freebayes,Platypus,HaplotypeCallerGVCF, 0 /mnt/xfs1/home/asalomatov/projects/pipeline/ppln/ 20 all
done
