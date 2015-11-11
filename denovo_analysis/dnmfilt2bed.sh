#!/bin/bash
fam=$1
cat ${fam}-DNMFilter.csv | awk '{split($0,a,","); print a[2]"\t"a[3]-1"\t"a[3]"\t"a[4]}' | intersectBed -a stdin -b /mnt/ceph/asalomatov/data/SSCexome/b37.exome-pm50.bed > ${fam}-DNMFilter.bed
