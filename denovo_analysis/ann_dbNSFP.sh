#!/bin/bash

#java -jar /mnt/xfs1/bioinfo/software/installs/bcbio_nextgen/150607/Cellar/snpeff/4.1g/libexec/SnpSift.jar dbnsfp -f chr,ref,alt,aaref,aaalt,rs_dbSNP142,genename,Ensembl_geneid,Ensembl_transcriptid,Ensembl_proteinid,SIFT_score SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,MutationTaster_score,MutationTaster_pred,MutationAssessor_score,MutationAssessor_pred,PROVEAN_score,PROVEAN_pred,CADD_raw,CADD_raw_rankscore,CADD_phred,GERP++_NR,GERP++_RS,GERP++_RS_rankscore,1000Gp3_AC,1000Gp3_AF,ExAC_AC,ExAC_AF -c /mnt/scratch/asalomatov/data/snpEff/snpEff.config -v -db /mnt/scratch/asalomatov/data/dbNSFP/hg19/dbNSFP3.0_hg19_sorted.txt.gz $1
filename=$(basename $1)
echo $filename
filenameWOext="${filename%.*}"
echo $filenameWOext
echo "$1 > ${filenameWOext}-nsfp.vcf"
java -jar /mnt/xfs1/bioinfo/software/installs/bcbio_nextgen/150607/Cellar/snpeff/4.1g/libexec/SnpSift.jar dbnsfp \
    -a -c /mnt/scratch/asalomatov/data/snpEff/snpEff.config \
    -db /mnt/scratch/asalomatov/data/dbNSFP/hg19/dbNSFP3.0_hg19_sorted.txt.gz \
    -f 'chr,pos(1-coor),ref,alt,aaref,aaalt,rs_dbSNP142,genename,Ensembl_geneid,Ensembl_transcriptid,Ensembl_proteinid,SIFT_score,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,MutationTaster_score,MutationTaster_pred,MutationAssessor_score,MutationAssessor_pred,PROVEAN_score,PROVEAN_pred,CADD_raw,CADD_raw_rankscore,CADD_phred,GERP++_NR,GERP++_RS,GERP++_RS_rankscore,1000Gp3_AC,1000Gp3_AF,ExAC_AC,ExAC_AF' $1 | sed "s/dbNSFP_GERP++/dbNSFP_GERP/g" | sed "s/pos(1-coor)/pos/g" > ${filenameWOext}-nsfp.vcf
#'pos(1-based)',
