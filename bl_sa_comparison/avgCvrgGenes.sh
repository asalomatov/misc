#!/bin/bash

for x in $(ls *-genes.bed); do echo $x;  cat $x | grep -v ^all | awk '{print $4"\t"$6*$9}' | sort -V -k1,1 | bedtools groupby -g 1 -c 2 -o sum > ${x}.txt; done &
