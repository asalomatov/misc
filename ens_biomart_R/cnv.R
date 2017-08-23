library('biomaRt')
listMarts()
listMarts(host='grch37.ensembl.org', archive=FALSE)
nsmbl=useMart("ENSEMBL_MART_ENSEMBL", host='grch37.ensembl.org')
listDatasets(nsmbl)
nsmbl = useDataset("hsapiens_gene_ensembl", mart=nsmbl)
head(listFilters(nsmbl))
head(listAttributes(nsmbl), 20)
atr <- (listAttributes(nsmbl))
dim(atr)
grep('name', atr$name, ignore.case=TRUE, value=TRUE)
grep('mim', atr$name, ignore.case=TRUE, value=TRUE)
grep('HGVSp', atr$name, ignore.case=TRUE, value=TRUE)

# read the data
cnv <- read.csv('b1-9.cnv.csv')
head(cnv)
str(cnv[1, c('Chr', 'Start', 'End')])
str(as.list(cnv[1:5, c('Chr', 'Start', 'End')]))
x <- as.list(cnv[1, c('Chr', 'Start', 'End')])
x
x <- list(10, 100008694, 101280279)
# may need to convert some cols from factors to strings later
x <- head(cnv[, c('Chr', 'Start', 'End')])
x

y <- getBM(attributes = c('external_gene_name',
                     'mim_gene_description',
                     'mim_gene_accession'),
      filters = c('chromosome_name', 'start', 'end'),
      values = list(x$Chr, x$Start, x$End),
      mart = nsmbl)
y
head(y)

type(

genesAtInterval <- function(chrom, startpos, endpos){
   y <- getBM(attributes = c('external_gene_name',
                     'mim_gene_description',
                     'mim_gene_accession'),
      filters = c('chromosome_name', 'start', 'end'),
      values = list(chrom, startpos, endpos),
      mart = nsmbl)
   if(nrow(y) == 0) return(data.frame())
   y$Chr <- chrom
   y$Start <- startpos
   y$End <- endpos
   y
}
genesAtInterval(5, 139931633, 139931740)
genesAtInterval("X", 2700157, 3265143)
genesAtRow <- function(i, mydf){
   print(i)
   genesAtInterval(mydf$CHROM[i],
                   mydf$Start[i],
                   mydf$End[i])

}
genesAtRow(98, cnv)
cnv[98,]
zzz <- sapply(1:nrow(cnv), genesAtRow, mydf=cnv, simplify=FALSE)
zzz
anno  <- do.call(rbind, zzz)
head(anno)
dim(anno)
write.csv(anno, row.names=FALSE, quote=FALSE, file = 'cnv_genes_omim.csv')

z <- mapply(genesAtInterval, as.character(x$Chr), x$Start, x$End)
?mapply
str(z)
as.data.frame(do.call(cbind, z))
head(as.data.frame(do.call(cbind, z)))

# try with plyr
install.packages('plyr')
library(plyr)
zzz <- mdply(data.frame(chrom = cnv$Chr,
                 startpos = cnv$Start,
                 endpos = cnv$End),
      genesAtInterval)

sapply

dim(cnv)
names(cnv)
length(unique(cnv$SampleID))
ids.by.batch <- read.csv('/mnt/xfs1/scratch/asalomatov/data/SPARK/info/ids_by_batch.csv')
head(ids.by.batch)
nocnv <- ids.by.batch[!(ids.by.batch$lab_id %in% unique(cnv$SampleID)) & ids.by.batch$lab_id %in% children, c('lab_id', 'batch')]
length(nocnv)

ped <- read.table('/mnt/xfs1/scratch/asalomatov/data/SPARK/info/b1-9_lab_id_ped.tsv', header=TRUE)
head(ped)
children <- ped$personId[(ped$dadId != '0') & (ped$momId != '0')]
children
nocnv <- nocnv[nocnv %in% children]
nocnv
length(nocnv)
write.csv(nocnv, row.names=FALSE, quote=FALSE, file = '/mnt/xfs1/scratch/asalomatov/data/SPARK/info/no_cnv_report.csv')
ig  <- strsplit(ig, ",")
y$external_gene_name %in% ig
write.csv(y, 'chr10-100008694-101280279.csv')


# do the same for VIP vars

x <- read.csv('/mnt/ceph/users/asalomatov/VIP/denovo/SF_lofdmis_anno.csv')
head(x$REF)
length(as.character(x$REF))
length(sapply(as.character(x$REF), nchar))
lenref <- sapply(as.character(x$REF), nchar)
names(lenref)  <- NULL
lenalt <- sapply(as.character(x$ALT), nchar)
names(lenalt)  <- NULL
x$Start <- x$POS
x$End <- x$POS + pmax(lenref, lenalt) - 1
x$chrSE <- paste(x$CHROM, x$Start, x$End, sep='_')
zzz <- sapply(1:nrow(x), genesAtRow, mydf=x, simplify=FALSE)
zzz
anno  <- do.call(rbind, zzz)
head(anno)
dim(anno)
anno$chrSE <- paste(anno$Chr, anno$Start, anno$End, sep='_')
library(plyr)
z <- ddply(anno, .(chrSE), summarize,
      genes = paste(external_gene_name, collapse='|'),
      omim = paste(mim_gene_description, collapse='|')
      )
head(z)
head(x)
dim(z)
length(unique(z$chrSE))
dim(x)
dim(merge(x, z, by='chrSE'))
x <- merge(x, z, by='chrSE', all.x=T, all.y=F)
pheno3 <- read.csv('/mnt/ceph/users/asalomatov/VIP/info/pheno_report_3.csv')
pheno4 <- read.csv('/mnt/ceph/users/asalomatov/VIP/info/pheno_report_4.csv')
vipsmpl <- read.csv('/mnt/ceph/users/asalomatov/VIP/info/svip_pedigrees.csv')

names(pheno3)
names(x)
head(merge(x, pheno3[,c('person_id', 'svip_neuro_exam.measure.eval_age_months')], by.x='lab_id', by.y='person_id'), all.x=T, all.y=F)
x <- merge(x, pheno3[,c('person_id', 'svip_neuro_exam.measure.eval_age_months')], by.x='lab_id', by.y='person_id', all.x=T, all.y=F)
dim(x)
head(x)
head(pheno4)
colSums(is.na(pheno4))
x <- merge(x, pheno4[,c('person_id', 'diagnosis_summary.best_full_scale_iq')], by.x='lab_id', by.y='person_id', all.x=T, all.y=F)
x$lab_id %in% pheno3$person_id
head(vipsmpl)
x <- merge(x, vipsmpl[,c('individual', 'genetic_status_16p')], by.x='lab_id', by.y='individual', all.x=T, all.y=F)
x <- merge(x, vipsmpl[,c('individual', 'father')], by.x='lab_id', by.y='individual', all.x=T, all.y=F)
x <- merge(x, vipsmpl[,c('individual', 'genetic_status_16p')], by.x='father', by.y='individual', all.x=T, all.y=F, suffixes=c('','_fa'))
x <- merge(x, vipsmpl[,c('individual', 'mother')], by.x='lab_id', by.y='individual', all.x=T, all.y=F)
x <- merge(x, vipsmpl[,c('individual', 'genetic_status_16p')], by.x='mother', by.y='individual', all.x=T, all.y=F, suffixes=c('','_mo'))
head(x)
head(x$father)
colSums(is.na(x))
colSums(is.na(vipsmpl))
table(vipsmpl$genetic_status_16p)
names(x)

write.csv(x[, names(x)[!(names(x) %in% c('chrSE','SP_id','SF', 'BCM','Cdfd','Clmb', 'allele_frac.1', 'effect_cat.1','status','batch','pos_a_tr','Start','End'))]], row.names=FALSE, quote=FALSE, file = '/mnt/ceph/users/asalomatov/VIP/denovo/SF_lofdmis_anno_omim.tsv')
write.table(x[, names(x)[!(names(x) %in% c('chrSE','SP_id','SF', 'BCM','Cdfd','Clmb', 'allele_frac.1', 'effect_cat.1','status','batch','pos_a_tr','Start','End'))]], row.names=FALSE, quote=FALSE, file = '/mnt/ceph/users/asalomatov/VIP/denovo/SF_lofdmis_anno_omim.tsv')
