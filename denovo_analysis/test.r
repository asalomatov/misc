source("http://bioconductor.org/biocLite.R")
biocLite("VariantAnnotation")
library(VariantAnnotation)
library(cgdv17)
iosd <- system.file("vcf", "ios_denovo-pm50-ann.vcf.gza")
hdr <- scanVcfHeader(iosd)
info(hdr)

