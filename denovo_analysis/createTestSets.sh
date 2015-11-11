#!/bin/bash

fam=$(cat /mnt/ceph/asalomatov/data/SSCexome/DNM_traning/fam_187.txt)
fam_quad=$(cat /mnt/ceph/asalomatov/data/SSCexome/DNM_traning/fam_quad_62.txt)

echo snp for s1
for f in $fam_quad; do cat $f-*-pm50-ann-s1.txt | grep -P '\s+SNP\s+' | awk -v fml=$f '{print fml" "$1" "$2}'  >> all_187_s1_snp_test.txt; done
cat all_187_s1_snp_test.txt | sort -V -u -k1,1 -k2,2 -k3,3 | awk '{print $1","$2","$3}' > all_187_s1_snp_test.csv
wc -l all_187_s1_snp_test.txt
wc -l all_187_s1_snp_test.csv
echo " "

echo ins for s1
for f in $fam_quad; do cat $f-*-pm50-ann-s1.txt | grep -P '\s+INS\s+' | awk -v fml=$f '{for (i=1;i<length($5);i++) print fml" "$1" "$2+i}'  >> all_187_s1_ins_test.txt; done
cat all_187_s1_ins_test.txt | sort -V -u -k1,1 -k2,2 -k3,3 | awk '{print $1","$2","$3}' > all_187_s1_ins_test.csv
wc -l all_187_s1_ins_test.txt
wc -l all_187_s1_ins_test.csv
echo " "

echo del for s1
for f in $fam_quad; do cat $f-*-pm50-ann-s1.txt | grep -P '\s+DEL\s+' | awk -v fml=$f '{for (i=1;i<length($4);i++) print fml" "$1" "$2+i}'  >> all_187_s1_del_test.txt; done
cat all_187_s1_del_test.txt | sort -V -u -k1,1 -k2,2 -k3,3 | awk '{print $1","$2","$3}' > all_187_s1_del_test.csv
wc -l all_187_s1_del_test.txt
wc -l all_187_s1_del_test.csv
echo " "

echo snp for p1
for f in $fam; do cat $f-*-pm50-ann-p1.txt | grep -P '\s+SNP\s+' | awk -v fml=$f '{print fml" "$1" "$2}'  >> all_187_p1_snp_test.txt; done
cat all_187_p1_snp_test.txt | sort -V -u -k1,1 -k2,2 -k3,3 | awk '{print $1","$2","$3}' > all_187_p1_snp_test.csv
wc -l all_187_p1_snp_test.txt
wc -l all_187_p1_snp_test.csv
echo " "

echo ins for p1
for f in $fam; do cat $f-*-pm50-ann-p1.txt | grep -P '\s+INS\s+' | awk -v fml=$f '{for (i=1;i<length($5);i++) print fml" "$1" "$2+i}'  >> all_187_p1_ins_test.txt; done
cat all_187_p1_ins_test.txt | sort -V -u -k1,1 -k2,2 -k3,3 | awk '{print $1","$2","$3}' > all_187_p1_ins_test.csv
wc -l all_187_p1_ins_test.txt
wc -l all_187_p1_ins_test.csv
echo " "

echo del for p1
for f in $fam; do cat $f-*-pm50-ann-p1.txt | grep -P '\s+DEL\s+' | awk -v fml=$f '{for (i=1;i<length($4);i++) print fml" "$1" "$2+i}'  >> all_187_p1_del_test.txt; done
cat all_187_p1_del_test.txt | sort -V -u -k1,1 -k2,2 -k3,3 | awk '{print $1","$2","$3}' > all_187_p1_del_test.csv
wc -l all_187_p1_del_test.txt
wc -l all_187_p1_del_test.csv
echo " "
