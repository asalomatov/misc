#!/bin/bash
cp ../ios_vcf_header.vcf 12011-wgs.vcf
cat /mnt/ceph/NYG/Project_REI_01344_WGS_2014_11_18/Quad.Analysis/Family.12011/SSC05559_SSC05555_SSC05551_SSC05560.annotated.consize.txt | grep -v ^CHROM | awk '{split($0,a,"\t"); x="";for (i=1;i<=length(a);i++) x=x"__"a[i]; print $1"\t"$2"\t"".""\t"$3"\t"$4"\t"".""\t"".""\t""FAM=12011"";MISC_INFO="x}' >> 12011-wgs.vcf
