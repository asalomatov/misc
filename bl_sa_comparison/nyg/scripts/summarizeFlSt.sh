#!/bin/bash

fam_ids=$(echo 'GEN15-40-D1 GEN15-40F-D1 GEN15-40M-D1 GEN15-41-D1 GEN15-41F-D1 GEN15-41M-D1 GEN15-42-D1 GEN15-42F-D1 GEN15-42M-D1 GEN15-43-D1 GEN15-43F-D1 GEN15-43M-D1 GEN15-44-D1 GEN15-44F-D1 GEN15-44M-D1 GEN15-45-D1 GEN15-45F-D1 GEN15-46-D1 GEN15-46F-D1 GEN15-46M-D1 GEN15-46R-D1 GEN15-54-D1 GEN15-54F-D1 GEN15-54M-D1 GEN15-55-D1 GEN15-56-D1')

for f in $fam_ids ; do descr='final'; x=$(grep total ${f}*${descr}.bam.FlSt); echo "$f $descr $x" >> \
    BlVsSlv_Fl_total_reads.txt; done
for f in $fam_ids ; do descr='final'; x=$(grep duplicates ${f}*${descr}.bam.FlSt); echo "$f $descr $x" >> \
    BlVsSlv_Fl_duplicate_reads.txt; done
for f in $fam_ids ; do descr='final'; x=$(cat ${f}*${descr}.bam.FlSt | grep mapped | grep -v with); echo "$f $descr $x" >> \
    BlVsSlv_Fl_mapped_reads.txt; done

for f in $fam_ids ; do descr='exome'; x=$(grep total ${f}*${descr}.bam.FlSt); echo "$f $descr $x" >> \
    BlVsSlv_Fl_total_reads.txt; done
for f in $fam_ids ; do descr='exome'; x=$(grep duplicates ${f}*${descr}.bam.FlSt); echo "$f $descr $x" >> \
    BlVsSlv_Fl_duplicate_reads.txt; done
for f in $fam_ids ; do descr='exome'; x=$(cat ${f}*${descr}.bam.FlSt | grep mapped | grep -v with); echo "$f $descr $x" >> \
    BlVsSlv_Fl_mapped_reads.txt; done

for f in $fam_ids ; do descr='exome-extra'; x=$(grep total ${f}*${descr}.bam.FlSt); echo "$f $descr $x" >> \
    BlVsSlv_Fl_total_reads.txt; done
for f in $fam_ids ; do descr='exome-extra'; x=$(grep duplicates ${f}*${descr}.bam.FlSt); echo "$f $descr $x" >> \
    BlVsSlv_Fl_duplicate_reads.txt; done
for f in $fam_ids ; do descr='exome-extra'; x=$(cat ${f}*${descr}.bam.FlSt | grep mapped | grep -v with); echo "$f $descr $x" >> \
    BlVsSlv_Fl_mapped_reads.txt; done

for f in $fam_ids ; do descr='genes'; x=$(grep total ${f}*${descr}.bam.FlSt); echo "$f $descr $x" >> \
    BlVsSlv_Fl_total_reads.txt; done
for f in $fam_ids ; do descr='genes'; x=$(grep duplicates ${f}*${descr}.bam.FlSt); echo "$f $descr $x" >> \
    BlVsSlv_Fl_duplicate_reads.txt; done
for f in $fam_ids ; do descr='genes'; x=$(cat ${f}*${descr}.bam.FlSt | grep mapped | grep -v with); echo "$f $descr $x" >> \
    BlVsSlv_Fl_mapped_reads.txt; done

for f in $fam_ids ; do descr='genes-extra'; x=$(grep total ${f}*${descr}.bam.FlSt); echo "$f $descr $x" >> \
    BlVsSlv_Fl_total_reads.txt; done
for f in $fam_ids ; do descr='genes-extra'; x=$(grep duplicates ${f}*${descr}.bam.FlSt); echo "$f $descr $x" >> \
    BlVsSlv_Fl_duplicate_reads.txt; done
for f in $fam_ids ; do descr='genes-extra'; x=$(cat ${f}*${descr}.bam.FlSt | grep mapped | grep -v with); echo "$f $descr $x" >> \
    BlVsSlv_Fl_mapped_reads.txt; done


