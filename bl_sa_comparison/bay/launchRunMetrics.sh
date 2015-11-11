#!/bin/bash

fam_ids='HM7G3ADXX-1-ID04 HM7G3ADXX-1-ID05 HM7G3ADXX-1-ID06 HM7G3ADXX-2-ID07 HM7G3ADXX-2-ID08 HM7G3ADXX-2-ID09 HM7NHADXX-1-ID10 HM7NHADXX-1-ID11 HM7NHADXX-1-ID12 HM7NHADXX-2-ID01 HM7NHADXX-2-ID02 HM7NHADXX-2-ID12 HMFKKADXX-1-ID09 HMFKKADXX-1-ID10 HMFKKADXX-1-ID11 HMFKKADXX-2-ID01 HMFKKADXX-2-ID02 HTFM5ADXX-1-ID03 HTFM5ADXX-1-ID04 HTFM5ADXX-1-ID05 HTFM5ADXX-2-ID06 HTFM5ADXX-2-ID07 HTFM5ADXX-2-ID08 HTFW7ADXX-2-ID01 HTFW7ADXX-2-ID02 HTFW7ADXX-2-ID03'

datadir='/mnt/xfs1/scratch/asalomatov/baylor_qc'

for f in $fam_ids
do
    echo "sbatch -J $f --ntasks-per-node=4 -D ./ -o ${f}.out -e ${f}.err ${datadir}/scripts/runMetrics.sh $datadir \
        ${datadir}/output $f "
    sbatch -J $f --ntasks-per-node=4 -D ${datadir}/slurm -o ${f}.out -e ${f}.err ${datadir}/scripts/runMetrics.sh $datadir \
        ${datadir}/output $f
done

