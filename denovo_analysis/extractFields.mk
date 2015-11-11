### 
default: all
SHELL = /bin/bash
USR = $(shell whoami)
#INCLMK = ~/projects/pipeline/ppln/include.mk
include $(INCLMK)
### may override on cl
snpsift = /mnt/xfs1/bioinfo/software/installs/bcbio_nextgen/150607/Cellar/snpeff/4.1g/libexec/SnpSift.jar
snpeff = /mnt/xfs1/bioinfo/software/installs/bcbio_nextgen/150607/Cellar/snpeff/4.1g/libexec/snpEff.jar
snpeffconfig = /mnt/scratch/asalomatov/data/snpEff/snpEff.config
DATASOURCE = SF_JHC_p1
PREFIX = 1
SUFFIX = 
EXT = .vcf
INDIR = .
OUTDIR = .
LOGDIR = $(OUTDIR)
TMPDIR = /tmp/$(USR)
###
inFile = $(wildcard $(INDIR)/$(PREFIX)*$(SUFFIX)$(EXT))
#inFile = $(wildcard $(INDIR)/$(PREFIX)*$(SUFFIX).vcf)
$(info $(inFile))
#outFile = $(addprefix $(OUTDIR)/, $(patsubst %$(SUFFIX).vcf.gz,%$(SUFFIX).txt,$(notdir $(inFile))))
#$(info $(outFile))

define buildTarget_JHC_p1

$(eval targ = $(addprefix $(OUTDIR)/, $(patsubst %$(SUFFIX)$(EXT),%$(SUFFIX).txt,$(notdir $(1)))))
targs += $(targ)
$(eval dep = $(1))

$(targ): $(dep)
	mkdir -p $(OUTDIR)
	mkdir -p $(TMPDIR)
	mkdir -p $(OUTDIR)
	java -Xmx5G -jar $(snpsift) extractFields -s "__" -e "." -c $(snpeffconfig) -noLog $$< CHROM POS ID REF ALT QUAL AC AF AN BaseQRankSum ClippingRankSum DP DS END FS GC HRun HaplotypeScore InbreedingCoeff MLEAC MLEAF MQ MQ0 MQRankSum QD ReadPosRankSum SOR LOF NMD VARTYPE "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].FEATURE" "GEN[0].GT[*]" "GEN[1].GT[*]" "GEN[2].GT[*]" "GEN[0].AD[*]" "GEN[1].AD[*]" "GEN[2].AD[*]" "GEN[0].DP[*]" "GEN[1].DP[*]" "GEN[2].DP[*]" "GEN[0].GQ[*]" "GEN[1].GQ[*]" "GEN[2].GQ[*]" "GEN[0].PL[*]" "GEN[1].PL[*]" "GEN[2].PL[*]" "GEN[0].RGQ[*]" "GEN[1].RGQ[*]" "GEN[2].RGQ[*]" "GEN[0].SB[*]" "GEN[1].SB[*]" "GEN[2].SB[*]" "GEN[*]" dbNSFP_chr dbNSFP_pos dbNSFP_ref dbNSFP_alt dbNSFP_aaref dbNSFP_aaalt dbNSFP_rs_dbSNP142 dbNSFP_genename dbNSFP_Ensembl_geneid dbNSFP_Ensembl_transcriptid dbNSFP_Ensembl_proteinid dbNSFP_SIFT_score dbNSFP_SIFT_pred dbNSFP_Polyphen2_HDIV_score dbNSFP_Polyphen2_HDIV_pred dbNSFP_Polyphen2_HVAR_score dbNSFP_Polyphen2_HVAR_pred dbNSFP_MutationTaster_score dbNSFP_MutationTaster_pred dbNSFP_MutationAssessor_score dbNSFP_MutationAssessor_pred dbNSFP_PROVEAN_score dbNSFP_PROVEAN_pred dbNSFP_CADD_raw dbNSFP_CADD_raw_rankscore dbNSFP_CADD_phred dbNSFP_GERP_NR dbNSFP_GERP_RS dbNSFP_GERP_RS_rankscore dbNSFP_1000Gp3_AC dbNSFP_1000Gp3_AF dbNSFP_ExAC_AC dbNSFP_ExAC_AF > $$@


endef      

define buildTarget_JHC_s1

$(eval targ = $(addprefix $(OUTDIR)/, $(patsubst %$(SUFFIX)$(EXT),%$(SUFFIX).txt,$(notdir $(1)))))
targs += $(targ)
$(eval dep = $(1))

$(targ): $(dep)
	mkdir -p $(OUTDIR)
	mkdir -p $(TMPDIR)
	mkdir -p $(OUTDIR)
	java -Xmx5G -jar $(snpsift) extractFields -s "__" -e "." -c $(snpeffconfig) -noLog $$< CHROM POS ID REF ALT QUAL AC AF AN BaseQRankSum ClippingRankSum DP DS END FS GC HRun HaplotypeScore InbreedingCoeff MLEAC MLEAF MQ MQ0 MQRankSum QD ReadPosRankSum SOR LOF NMD VARTYPE "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].FEATURE" "GEN[0].GT[*]" "GEN[1].GT[*]" "GEN[3].GT[*]" "GEN[0].AD[*]" "GEN[1].AD[*]" "GEN[3].AD[*]" "GEN[0].DP[*]" "GEN[1].DP[*]" "GEN[3].DP[*]" "GEN[0].GQ[*]" "GEN[1].GQ[*]" "GEN[3].GQ[*]" "GEN[0].PL[*]" "GEN[1].PL[*]" "GEN[3].PL[*]" "GEN[0].RGQ[*]" "GEN[1].RGQ[*]" "GEN[3].RGQ[*]" "GEN[0].SB[*]" "GEN[1].SB[*]" "GEN[3].SB[*]" "GEN[*]" dbNSFP_chr dbNSFP_pos dbNSFP_ref dbNSFP_alt dbNSFP_aaref dbNSFP_aaalt dbNSFP_rs_dbSNP142 dbNSFP_genename dbNSFP_Ensembl_geneid dbNSFP_Ensembl_transcriptid dbNSFP_Ensembl_proteinid dbNSFP_SIFT_score dbNSFP_SIFT_pred dbNSFP_Polyphen2_HDIV_score dbNSFP_Polyphen2_HDIV_pred dbNSFP_Polyphen2_HVAR_score dbNSFP_Polyphen2_HVAR_pred dbNSFP_MutationTaster_score dbNSFP_MutationTaster_pred dbNSFP_MutationAssessor_score dbNSFP_MutationAssessor_pred dbNSFP_PROVEAN_score dbNSFP_PROVEAN_pred dbNSFP_CADD_raw dbNSFP_CADD_raw_rankscore dbNSFP_CADD_phred dbNSFP_GERP_NR dbNSFP_GERP_RS dbNSFP_GERP_RS_rankscore dbNSFP_1000Gp3_AC dbNSFP_1000Gp3_AF dbNSFP_ExAC_AC dbNSFP_ExAC_AF > $$@

endef      

define buildTarget_HC_p1

$(eval targ = $(addprefix $(OUTDIR)/, $(patsubst %$(SUFFIX)$(EXT),%$(SUFFIX).txt,$(notdir $(1)))))
targs += $(targ)
$(eval dep = $(1))

$(targ): $(dep)
	mkdir -p $(OUTDIR)
	mkdir -p $(TMPDIR)
	mkdir -p $(OUTDIR) 
	java -Xmx5G -jar $(snpsift) extractFields -s "__" -e "." -c $(snpeffconfig) -noLog $$< CHROM POS ID REF ALT QUAL AC AF AN BaseQRankSum ClippingRankSum DP DS FS GC HRun HaplotypeScore InbreedingCoeff MLEAC MLEAF MQ MQ0 MQRankSum QD ReadPosRankSum SOR LOF NMD VARTYPE "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].FEATURE" "GEN[0].GT[*]" "GEN[1].GT[*]" "GEN[2].GT[*]" "GEN[0].AD[*]" "GEN[1].AD[*]" "GEN[2].AD[*]" "GEN[0].DP[*]" "GEN[1].DP[*]" "GEN[2].DP[*]" "GEN[0].GQ[*]" "GEN[1].GQ[*]" "GEN[2].GQ[*]" "GEN[0].PL[*]" "GEN[1].PL[*]" "GEN[2].PL[*]" "GEN[*]" dbNSFP_chr dbNSFP_pos dbNSFP_ref dbNSFP_alt dbNSFP_aaref dbNSFP_aaalt dbNSFP_rs_dbSNP142 dbNSFP_genename dbNSFP_Ensembl_geneid dbNSFP_Ensembl_transcriptid dbNSFP_Ensembl_proteinid dbNSFP_SIFT_score dbNSFP_SIFT_pred dbNSFP_Polyphen2_HDIV_score dbNSFP_Polyphen2_HDIV_pred dbNSFP_Polyphen2_HVAR_score dbNSFP_Polyphen2_HVAR_pred dbNSFP_MutationTaster_score dbNSFP_MutationTaster_pred dbNSFP_MutationAssessor_score dbNSFP_MutationAssessor_pred dbNSFP_PROVEAN_score dbNSFP_PROVEAN_pred dbNSFP_CADD_raw dbNSFP_CADD_raw_rankscore dbNSFP_CADD_phred dbNSFP_GERP_NR dbNSFP_GERP_RS dbNSFP_GERP_RS_rankscore dbNSFP_1000Gp3_AC dbNSFP_1000Gp3_AF dbNSFP_ExAC_AC dbNSFP_ExAC_AF > $$@


endef      

define buildTarget_HC_s1

$(eval targ = $(addprefix $(OUTDIR)/, $(patsubst %$(SUFFIX)$(EXT),%$(SUFFIX).txt,$(notdir $(1)))))
targs += $(targ)
$(eval dep = $(1))

$(targ): $(dep)
	mkdir -p $(OUTDIR)
	mkdir -p $(TMPDIR)
	mkdir -p $(OUTDIR)
	java -Xmx5G -jar $(snpsift) extractFields -s "__" -e "." -c $(snpeffconfig) -noLog $$< CHROM POS ID REF ALT QUAL AC AF AN BaseQRankSum ClippingRankSum DP DS FS GC HRun HaplotypeScore InbreedingCoeff MLEAC MLEAF MQ MQ0 MQRankSum QD ReadPosRankSum SOR LOF NMD VARTYPE "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].FEATURE" "GEN[0].GT[*]" "GEN[1].GT[*]" "GEN[3].GT[*]" "GEN[0].AD[*]" "GEN[1].AD[*]" "GEN[3].AD[*]" "GEN[0].DP[*]" "GEN[1].DP[*]" "GEN[3].DP[*]" "GEN[0].GQ[*]" "GEN[1].GQ[*]" "GEN[3].GQ[*]" "GEN[0].PL[*]" "GEN[1].PL[*]" "GEN[3].PL[*]" "GEN[*]" dbNSFP_chr dbNSFP_pos dbNSFP_ref dbNSFP_alt dbNSFP_aaref dbNSFP_aaalt dbNSFP_rs_dbSNP142 dbNSFP_genename dbNSFP_Ensembl_geneid dbNSFP_Ensembl_transcriptid dbNSFP_Ensembl_proteinid dbNSFP_SIFT_score dbNSFP_SIFT_pred dbNSFP_Polyphen2_HDIV_score dbNSFP_Polyphen2_HDIV_pred dbNSFP_Polyphen2_HVAR_score dbNSFP_Polyphen2_HVAR_pred dbNSFP_MutationTaster_score dbNSFP_MutationTaster_pred dbNSFP_MutationAssessor_score dbNSFP_MutationAssessor_pred dbNSFP_PROVEAN_score dbNSFP_PROVEAN_pred dbNSFP_CADD_raw dbNSFP_CADD_raw_rankscore dbNSFP_CADD_phred dbNSFP_GERP_NR dbNSFP_GERP_RS dbNSFP_GERP_RS_rankscore dbNSFP_1000Gp3_AC dbNSFP_1000Gp3_AF dbNSFP_ExAC_AC dbNSFP_ExAC_AF > $$@


endef      

define buildTarget_FB_p1

$(eval targ = $(addprefix $(OUTDIR)/, $(patsubst %$(SUFFIX)$(EXT),%$(SUFFIX).txt,$(notdir $(1)))))
targs += $(targ)
$(eval dep = $(1))

$(targ): $(dep)
	mkdir -p $(OUTDIR)
	mkdir -p $(TMPDIR)
	mkdir -p $(OUTDIR) 
	java -Xmx5G -jar $(snpsift) extractFields -s "__" -e "." -c $(snpeffconfig) -noLog $$< CHROM POS ID REF ALT QUAL VARTYPE "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].FEATURE" "GEN[0].GT[*]" "GEN[1].GT[*]" "GEN[2].GT[*]" "GEN[0].DP[*]" "GEN[1].DP[*]" "GEN[2].DP[*]" "GEN[0].GQ[*]" "GEN[1].GQ[*]" "GEN[2].GQ[*]" "GEN[*]" dbNSFP_chr dbNSFP_pos dbNSFP_ref dbNSFP_alt dbNSFP_aaref dbNSFP_aaalt dbNSFP_rs_dbSNP142 dbNSFP_genename dbNSFP_Ensembl_geneid dbNSFP_Ensembl_transcriptid dbNSFP_Ensembl_proteinid dbNSFP_SIFT_score dbNSFP_SIFT_pred dbNSFP_Polyphen2_HDIV_score dbNSFP_Polyphen2_HDIV_pred dbNSFP_Polyphen2_HVAR_score dbNSFP_Polyphen2_HVAR_pred dbNSFP_MutationTaster_score dbNSFP_MutationTaster_pred dbNSFP_MutationAssessor_score dbNSFP_MutationAssessor_pred dbNSFP_PROVEAN_score dbNSFP_PROVEAN_pred dbNSFP_CADD_raw dbNSFP_CADD_raw_rankscore dbNSFP_CADD_phred dbNSFP_GERP_NR dbNSFP_GERP_RS dbNSFP_GERP_RS_rankscore dbNSFP_1000Gp3_AC dbNSFP_1000Gp3_AF dbNSFP_ExAC_AC dbNSFP_ExAC_AF > $$@


endef      

define buildTarget_FB_s1

$(eval targ = $(addprefix $(OUTDIR)/, $(patsubst %$(SUFFIX)$(EXT),%$(SUFFIX).txt,$(notdir $(1)))))
targs += $(targ)
$(eval dep = $(1))

$(targ): $(dep)
	mkdir -p $(OUTDIR)
	mkdir -p $(TMPDIR)
	mkdir -p $(OUTDIR) 
	java -Xmx5G -jar $(snpsift) extractFields -s "__" -e "." -c $(snpeffconfig) -noLog $$< CHROM POS ID REF ALT QUAL VARTYPE "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].FEATURE" "GEN[0].GT[*]" "GEN[1].GT[*]" "GEN[3].GT[*]" "GEN[0].DP[*]" "GEN[1].DP[*]" "GEN[3].DP[*]" "GEN[0].GQ[*]" "GEN[1].GQ[*]" "GEN[3].GQ[*]" "GEN[*]" dbNSFP_chr dbNSFP_pos dbNSFP_ref dbNSFP_alt dbNSFP_aaref dbNSFP_aaalt dbNSFP_rs_dbSNP142 dbNSFP_genename dbNSFP_Ensembl_geneid dbNSFP_Ensembl_transcriptid dbNSFP_Ensembl_proteinid dbNSFP_SIFT_score dbNSFP_SIFT_pred dbNSFP_Polyphen2_HDIV_score dbNSFP_Polyphen2_HDIV_pred dbNSFP_Polyphen2_HVAR_score dbNSFP_Polyphen2_HVAR_pred dbNSFP_MutationTaster_score dbNSFP_MutationTaster_pred dbNSFP_MutationAssessor_score dbNSFP_MutationAssessor_pred dbNSFP_PROVEAN_score dbNSFP_PROVEAN_pred dbNSFP_CADD_raw dbNSFP_CADD_raw_rankscore dbNSFP_CADD_phred dbNSFP_GERP_NR dbNSFP_GERP_RS dbNSFP_GERP_RS_rankscore dbNSFP_1000Gp3_AC dbNSFP_1000Gp3_AF dbNSFP_ExAC_AC dbNSFP_ExAC_AF > $$@


endef      

define buildTarget_PL_p1

$(eval targ = $(addprefix $(OUTDIR)/, $(patsubst %$(SUFFIX)$(EXT),%$(SUFFIX).txt,$(notdir $(1)))))
targs += $(targ)
$(eval dep = $(1))

$(targ): $(dep)
	mkdir -p $(OUTDIR)
	mkdir -p $(TMPDIR)
	mkdir -p $(OUTDIR) 
	java -Xmx5G -jar $(snpsift) extractFields -s "__" -e "." -c $(snpeffconfig) -noLog $$< CHROM POS ID REF ALT QUAL VARTYPE "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].FEATURE" "GEN[0].GT[*]" "GEN[1].GT[*]" "GEN[2].GT[*]" "GEN[0].NR[*]" "GEN[1].NR[*]" "GEN[2].NR[*]" "GEN[0].GQ[*]" "GEN[1].GQ[*]" "GEN[2].GQ[*]" "GEN[*]" dbNSFP_chr dbNSFP_pos dbNSFP_ref dbNSFP_alt dbNSFP_aaref dbNSFP_aaalt dbNSFP_rs_dbSNP142 dbNSFP_genename dbNSFP_Ensembl_geneid dbNSFP_Ensembl_transcriptid dbNSFP_Ensembl_proteinid dbNSFP_SIFT_score dbNSFP_SIFT_pred dbNSFP_Polyphen2_HDIV_score dbNSFP_Polyphen2_HDIV_pred dbNSFP_Polyphen2_HVAR_score dbNSFP_Polyphen2_HVAR_pred dbNSFP_MutationTaster_score dbNSFP_MutationTaster_pred dbNSFP_MutationAssessor_score dbNSFP_MutationAssessor_pred dbNSFP_PROVEAN_score dbNSFP_PROVEAN_pred dbNSFP_CADD_raw dbNSFP_CADD_raw_rankscore dbNSFP_CADD_phred dbNSFP_GERP_NR dbNSFP_GERP_RS dbNSFP_GERP_RS_rankscore dbNSFP_1000Gp3_AC dbNSFP_1000Gp3_AF dbNSFP_ExAC_AC dbNSFP_ExAC_AF > $$@


endef      

define buildTarget_PL_s1

$(eval targ = $(addprefix $(OUTDIR)/, $(patsubst %$(SUFFIX)$(EXT),%$(SUFFIX).txt,$(notdir $(1)))))
targs += $(targ)
$(eval dep = $(1))

$(targ): $(dep)
	mkdir -p $(OUTDIR)
	mkdir -p $(TMPDIR)
	mkdir -p $(OUTDIR) 
	java -Xmx5G -jar $(snpsift) extractFields -s "__" -e "." -c $(snpeffconfig) -noLog $$< CHROM POS ID REF ALT QUAL VARTYPE "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].FEATURE" "GEN[0].GT[*]" "GEN[1].GT[*]" "GEN[3].GT[*]" "GEN[0].NR[*]" "GEN[1].NR[*]" "GEN[3].NR[*]" "GEN[0].GQ[*]" "GEN[1].GQ[*]" "GEN[3].GQ[*]" "GEN[*]" dbNSFP_chr dbNSFP_pos dbNSFP_ref dbNSFP_alt dbNSFP_aaref dbNSFP_aaalt dbNSFP_rs_dbSNP142 dbNSFP_genename dbNSFP_Ensembl_geneid dbNSFP_Ensembl_transcriptid dbNSFP_Ensembl_proteinid dbNSFP_SIFT_score dbNSFP_SIFT_pred dbNSFP_Polyphen2_HDIV_score dbNSFP_Polyphen2_HDIV_pred dbNSFP_Polyphen2_HVAR_score dbNSFP_Polyphen2_HVAR_pred dbNSFP_MutationTaster_score dbNSFP_MutationTaster_pred dbNSFP_MutationAssessor_score dbNSFP_MutationAssessor_pred dbNSFP_PROVEAN_score dbNSFP_PROVEAN_pred dbNSFP_CADD_raw dbNSFP_CADD_raw_rankscore dbNSFP_CADD_phred dbNSFP_GERP_NR dbNSFP_GERP_RS dbNSFP_GERP_RS_rankscore dbNSFP_1000Gp3_AC dbNSFP_1000Gp3_AF dbNSFP_ExAC_AC dbNSFP_ExAC_AF > $$@

endef      

#"GEN[*].AD[*]"  "GEN[*].DP[*]" "GEN[*].GQ[*]" "GEN[*].PL[*]" \ 
define buildTarget_p1

$(eval targ = $(addprefix $(OUTDIR)/, $(patsubst %$(SUFFIX)$(EXT),%$(SUFFIX).txt,$(notdir $(1)))))
targs += $(targ)
$(eval dep = $(1))

$(targ): $(dep)
	mkdir -p $(OUTDIR)
	mkdir -p $(TMPDIR)
	mkdir -p $(OUTDIR)
	java -Xmx5G -jar $(snpsift) extractFields -s "__" -e "." -c $(snpeffconfig) -noLog $$< CHROM POS ID REF ALT QUAL VARTYPE "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].FEATURE" "GEN[0].GT[*]" "GEN[1].GT[*]" "GEN[2].GT[*]" "GEN[0].DP[*]" "GEN[1].DP[*]" "GEN[2].DP[*]" "GEN[*]" > $$@

endef      

define buildTarget_s1

$(eval targ = $(addprefix $(OUTDIR)/, $(patsubst %$(SUFFIX)$(EXT),%$(SUFFIX).txt,$(notdir $(1)))))
targs += $(targ)
$(eval dep = $(1))

$(targ): $(dep)
	mkdir -p $(OUTDIR)
	mkdir -p $(TMPDIR)
	mkdir -p $(OUTDIR)
	java -Xmx5G -jar $(snpsift) extractFields -s "__" -e "." -c $(snpeffconfig) -noLog $$< CHROM POS ID REF ALT QUAL VARTYPE "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].FEATURE" "GEN[0].GT[*]" "GEN[1].GT[*]" "GEN[3].GT[*]" "GEN[0].DP[*]" "GEN[1].DP[*]" "GEN[3].DP[*]" "GEN[*]" > $$@

endef      


define buildTargetIos_denovo

$(eval targ = $(addprefix $(OUTDIR)/, $(patsubst %$(SUFFIX)$(EXT),%$(SUFFIX).txt,$(notdir $(1)))))
targs += $(targ)
$(eval dep = $(1))

$(targ): $(dep)
	mkdir -p $(OUTDIR)
	mkdir -p $(TMPDIR)
	mkdir -p $(OUTDIR)
	java -Xmx5G -jar $(snpsift) extractFields -s "," -e "." -c $(snpeffconfig) -noLog $$< CHROM POS ID REF ALT QUAL VARTYPE FAM CHILD "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].FEATURE" > $$@

endef      

define buildTargetMin

$(eval targ = $(addprefix $(OUTDIR)/, $(patsubst %$(SUFFIX)$(EXT),%$(SUFFIX).txt,$(notdir $(1)))))
targs += $(targ)
$(eval dep = $(1))

$(targ): $(dep)
	mkdir -p $(OUTDIR)
	mkdir -p $(TMPDIR)
	mkdir -p $(OUTDIR)
	java -Xmx5G -jar $(snpsift) extractFields -s "__" -e "." -c $(snpeffconfig) -noLog $$< CHROM POS ID REF ALT QUAL VARTYPE "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].FEATURE" > $$@

endef      

define buildTarget_trdnv

$(eval targ = $(addprefix $(OUTDIR)/, $(patsubst %$(SUFFIX)$(EXT),%$(SUFFIX).txt,$(notdir $(1)))))
targs += $(targ)
$(eval dep = $(1))

$(targ): $(dep)
	mkdir -p $(OUTDIR)
	mkdir -p $(TMPDIR)
	mkdir -p $(OUTDIR)
	java -Xmx5G -jar $(snpsift) extractFields -s "__" -e "." -c $(snpeffconfig) -noLog $$< CHROM POS ID REF ALT QUAL FILTER "GEN[0].DQ" "GEN[0].DGQ" > $$@

endef      

ifeq ($(DATASOURCE),SF_JHC_p1)
$(foreach f,$(inFile),$(eval $(call buildTarget_JHC_p1,$(f)))) 
endif
ifeq ($(DATASOURCE),SF_HC_p1)
$(foreach f,$(inFile),$(eval $(call buildTarget_HC_p1,$(f)))) 
endif
ifeq ($(DATASOURCE),SF_FB_p1)
$(foreach f,$(inFile),$(eval $(call buildTarget_FB_p1,$(f)))) 
endif
ifeq ($(DATASOURCE),SF_PL_p1)
$(foreach f,$(inFile),$(eval $(call buildTarget_PL_p1,$(f)))) 
endif

ifeq ($(DATASOURCE),SF_JHC_s1)
$(foreach f,$(inFile),$(eval $(call buildTarget_JHC_s1,$(f)))) 
endif
ifeq ($(DATASOURCE),SF_HC_s1)
$(foreach f,$(inFile),$(eval $(call buildTarget_HC_s1,$(f)))) 
endif
ifeq ($(DATASOURCE),SF_FB_s1)
$(foreach f,$(inFile),$(eval $(call buildTarget_FB_s1,$(f)))) 
endif
ifeq ($(DATASOURCE),SF_PL_s1)
$(foreach f,$(inFile),$(eval $(call buildTarget_PL_s1,$(f)))) 
endif

ifeq ($(DATASOURCE),trdnv)
$(foreach f,$(inFile),$(eval $(call buildTarget_trdnv,$(f)))) 
endif

ifeq ($(DATASOURCE),SF_s1)
$(foreach f,$(inFile),$(eval $(call buildTarget_s1,$(f)))) 
endif
ifeq ($(DATASOURCE),ios_denovo)
$(foreach f,$(inFile),$(eval $(call buildTargetIos_denovo,$(f)))) 
endif
ifeq ($(DATASOURCE),min)
$(foreach f,$(inFile),$(eval $(call buildTargetMin,$(f)))) 
endif

all: $(targs)
