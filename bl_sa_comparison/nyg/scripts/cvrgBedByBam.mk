### 
#default: all
SHELL = /bin/bash
USR = $(shell whoami)
include $(INCLMK)
### may override on cl
PREFIX =
SUFFIX =
NEWSUFFIX = -targ
TARGBED = 
PROJ = mktest
INDIR = .
OUTDIR = .
LOGDIR = $(OUTDIR)
TMPDIR = /tmp/$(USR)/$(PROJ)
###
inFile = $(wildcard $(INDIR)/$(PREFIX)*$(SUFFIX).bam)
$(info $(inBam))
o0 = $(addprefix $(OUTDIR)/, $(patsubst %$(SUFFIX).bam,%$(SUFFIX)$(NEWSUFFIX).bed,$(notdir $(inFile))))
$(info $(o0))

all: $(o0)

$(OUTDIR)/%$(SUFFIX)$(NEWSUFFIX).bed: $(INDIR)/%$(SUFFIX).bam
	mkdir -p $(OUTDIR)
	mkdir -p $(TMPDIR)
	mkdir -p $(OUTDIR)
	bedtools coverage -sorted -hist -a $(TARGBED) -b $< > $@

