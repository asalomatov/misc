#!/bin/bash

for f in $fams; do sbatch -J $f -D ./ -o ${f}.out -e ${f}.err --exclusive
    /mnt/xfs1/home/asalomatov/projects/misc_scripts/baylor_denovo_qc/runDNMFilter.sh $f; done
