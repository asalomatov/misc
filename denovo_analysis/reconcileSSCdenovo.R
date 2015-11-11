load("denovo.RData")
ls()
str(denovo)
dim(denovo)
str(denovo.187)
library(VariantAnnotation)
fl <- system.file("tab", "ios_mut_status_norm.vcf.gz", package="VariantAnnotator")
fl
list.files(pattern="^ios")
mut.status <- readVcf("ios_mut_status_norm.vcf.gz", "hg19")
?readVcf
