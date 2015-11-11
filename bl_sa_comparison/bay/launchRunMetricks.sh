#!/bin/bash

fam_ids='GEN15-40-D1 GEN15-40F-D1 GEN15-40M-D1 GEN15-41-D1 GEN15-41F-D1 GEN15-41M-D1 GEN15-42-D1 GEN15-42F-D1 GEN15-42M-D1 GEN15-43-D1 GEN15-43F-D1 GEN15-43M-D1 GEN15-44-D1 GEN15-44F-D1 GEN15-44M-D1 GEN15-45-D1 GEN15-45F-D1 GEN15-46-D1 GEN15-46F-D1 GEN15-46M-D1 GEN15-46R-D1 GEN15-54-D1 GEN15-54F-D1 GEN15-54M-D1 GEN15-55-D1 GEN15-56-D1'

datadir='/mnt/xfs1/scratch/asalomatov/BloodVsSaliva'

for f in $fam_ids
do
    echo "sbatch -J $f --ntasks-per-node=4 -D ./ -o ${f}.out -e ${f}.err ${datadir}/scripts/runMetrics.sh $datadir \
        ${datadir}/output $f "
    sbatch -J $f --ntasks-per-node=4 -D ${datadir}/slurm -o ${f}.out -e ${f}.err ${datadir}/scripts/runMetrics.sh $datadir \
        ${datadir}/output $f
done

