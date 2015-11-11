#!/bin/bash

inpfile=$1
outpfile=$2

cat $inpfile | sort -u -V -k1,1 -k2,2 > ${inpfile}__temp
grep ^MT ${inpfile}__temp > ${inpfile}__temp_MT
grep -v ^MT ${inpfile}__temp > ${outpfile}
cat ${inpfile}__temp_MT >> ${outpfile}
rm ${inpfile}__temp*
