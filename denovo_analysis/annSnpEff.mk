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
		$(BGZIP) -c > $@
	$(TABIX) -f -p vcf $@


#        $(JAVA) -Xmx5G -jar $(SNPSIFTJAR) dbnsfp -v -c $(SNPEFFCONF) -noLog -a - | \
