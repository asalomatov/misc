#!/bin/bash

for m in $(echo p1 s1)
do
    for t in $(echo snp ins del)
    do
    sbatch -J SSC${m}$t -e SSC${m}${t}.err -o SSC${m}${t}.out --exclusive -N 1 \
        /mnt/ceph/asalomatov/dnmfilt/runDNMfilter1.sh \
        ../ped/all-${m}.ped ../bam/all-${m}-bams.txt ../train_sets/SSC_comb.csv \
        ../test_sets/all_187_${m}_${t}_test.csv ../output/SSC_comb/all-${m}-${t}-DNMFilter.csv
    done
done

for m in $(echo p1 s1)
do
    for t in $(echo snp ins del)
    do
    sbatch -J comb${m}$t -e comb${m}${t}.err -o comb${m}${t}.out --exclusive -N 1 \
        /mnt/ceph/asalomatov/dnmfilt/runDNMfilter1.sh \
        ../ped/all-${m}.ped ../bam/all-${m}-bams.txt ../train_sets/combined.csv \
        ../test_sets/all_187_${m}_${t}_test.csv ../output/combined/all-${m}-${t}-DNMFilter.csv
    done
done
