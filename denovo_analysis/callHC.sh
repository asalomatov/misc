#!/bin/bash 
FAMCODE=14011
BEDREGIONS=/mnt/ceph/asalomatov/data/SSCexome/b37.exome-pm50.bed
DBSNP=/mnt/xfs1/bioinfo/data/bcbio_nextgen/150607/genomes/Hsapiens/GRCh37/variation/dbsnp_138.vcf.gz
GENOMEREF=/mnt/xfs1/bioinfoCentos7/data/bcbio_nextgen/150617/genomes/Hsapiens/GRCh37/seq/GRCh37.fa
JAVA=/usr/java/jdk1.7.0_65/bin/java

###
$(eval bams = $(wildcard $(INDIR)/$(FAMCODE)*$(SUFFIX).bam))
$(info $(bams))

$(eval targ = $(OUTDIR)/$(FAMCODE)-HC-$(1)-bin.vcf.gz)
targs += $(targ)
$(eval dep1 = $(INDIR)/$(1)__bin__$(FAMCODE)-uni-mrg.bed)
$(eval dep2 = $(wildcard $(INDIR)/$(FAMCODE)*$(SUFFIX).bam))

$(targ): $(dep1) $(dep2)
	python $(SRCDIR)/gatkHaplotypeCaller.py $(GENOMEREF) $(TMPDIR) $(GATK) $(DBSNP) $(GAPS) $(LOGDIR) $$@ $$^

endef      
#$(info $(targ))

define runCallerSplitBam

$(eval targ = $(OUTDIR)/$(FAMCODE)-HC-$(1)-bin.vcf.gz)
targs += $(targ)
$(eval dep1 = $(INDIR)/$(1)__bin__$(FAMCODE)-uni-mrg.bed)
$(eval dep2 = $(wildcard $(INDIR)/$(FAMCODE)*-$(1)$(SUFFIX).bam))

$(targ): $(dep1) $(dep2)
	python $(SRCDIR)/gatkHaplotypeCaller.py $(GENOMEREF) $(TMPDIR) $(GATK) $(DBSNP) $(GAPS) $(LOGDIR) $$@ $$^

endef      

#$(info $(targ))
ifeq ($(SUFFIX),-bin)
$(foreach bin,$(num_bins),$(eval $(call runCallerSplitBam,$(bin)))) 
else
$(foreach bin,$(num_bins),$(eval $(call runCaller,$(bin)))) 
endif

all: $(targs)

