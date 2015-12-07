### 
SHELL = /bin/bash
USR = $(shell whoami)
INCLMK = ~/projects/pipeline/ppln/include.mk
include $(INCLMK)
### may override on cl
PREFIX = 1
SUFFIX = 
INDIR = .
OUTDIR = .
LOGDIR = $(OUTDIR)
TMPDIR = /tmp/$(USR)
###
inFile = $(wildcard $(INDIR)/$(PREFIX)*$(SUFFIX))
$(info $(inFile))
outFile = $(addprefix $(OUTDIR)/, $(patsubst %$(SUFFIX),%-ann$(SUFFIX),$(notdir $(inFile))))
sumFile = $(addprefix $(OUTDIR)/, $(patsubst %$(SUFFIX),%-ann$(SUFFIX).summary.html,$(notdir $(inFile))))
$(info $(outFile))

all: $(outFile)

$(outFile): $(inFile)
	mkdir -p $(OUTDIR)
	mkdir -p $(TMPDIR)
	mkdir -p $(OUTDIR)
	$(JAVA) -Xmx5G -jar $(SNPSIFTJAR) annotate -id -c $(SNPEFFCONF) -noLog -dbsnp $< | \
        $(JAVA) -Xmx5G -jar $(SNPSIFTJAR) varType -c $(SNPEFFCONF) -noLog - | \
		$(JAVA) -Xmx5G -jar $(SNPEFFJAR) ann -noLog -c $(SNPEFFCONF) $(SNPEFFGENOME) -v -csvStats -s $(sumFile) | \
		$(JAVA) -Xmx5G -jar $(SNPSIFTJAR)  dbnsfp - -a -c $(SNPEFFCONF) -db $(DBNSFP) \
		-f ref,alt,aaref,aaalt,rs_dbSNP142,genename,Ensembl_geneid,Ensembl_transcriptid,Ensembl_proteinid,SIFT_score,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,MutationTaster_score,MutationTaster_pred,MutationAssessor_score,MutationAssessor_pred,PROVEAN_score,PROVEAN_pred,CADD_raw,CADD_raw_rankscore,CADD_phred,GERP++_NR,GERP++_RS,GERP++_RS_rankscore,1000Gp3_AC,1000Gp3_AF,ExAC_AC,ExAC_AF | sed "s/dbNSFP_GERP++/dbNSFP_GERP/g" | $(BGZIP) -c > $@
	$(TABIX) -f -p vcf $@


#        $(JAVA) -Xmx5G -jar $(SNPSIFTJAR) dbnsfp -v -c $(SNPEFFCONF) -noLog -a - | \
