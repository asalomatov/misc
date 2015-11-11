#!/bin/bash

protoc=$1

for i in $(echo 1 2 4)
do
    sbatch -J  ${i}_$protoc -o ${i}_${protoc}.out -e ${i}_${protoc}.err -N 1 --exclusive -D ./ /mnt/xfs1/home/asalomatov/projects/pipeline/ppln/pipe03.sh /mnt/xfs1/scratch/asalomatov/BloodVsSaliva/inputs_$protoc ./$protoc $i WG 0 tmp /mnt/xfs1/home/asalomatov/projects/pipeline/ppln/include.mk 1 ,FixGroups,HaplotypeCaller,Freebayes,Platypus,HaplotypeCallerGVCF, 1 /mnt/xfs1/home/asalomatov/projects/pipeline/ppln/ 20 all
done
