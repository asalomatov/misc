save.image("denovo.RData")
load("denovo.RData")
ls()
effects.lof <- strsplit("exon_loss_variant|frameshift_variant|stop_gained|stop_lost|start_lost|splice_acceptor_variant|splice_donor_variant|splice_region_variant", split="|", fixed=TRUE)[[1]]
effects.miss <- "missense_variant"
effects.syn <- "synonymous_variant"
print(effects.lof)
print(fnInsertGenes)
print(fnInsertEffects)
print(genes)
fnInsertGenes <- function(mydf, genes, splt, other.name="other"){
    x <- sapply(mydf$gene, function(x){ 
        g <- unlist(strsplit(x, split=splt)) 
        g.coord <- match(g, genes)
        g.coord <- g.coord[!is.na(g.coord)][1]
        if (is.na(g.coord))
            return(other.name)
        else
            return(genes[g <- coord])
                  })
}
fnInsertEffects <- function(mydf,field, lof, miss, splt, other.name="other"){
   x <- sapply(mydf[,field], function(x){ 
       e <- unlist(strsplit(x, split=splt)) 
       e_coord <- match(e, lof)
       e_coord <- e_coord[!is.na(e_coord)][1]
       if (!is.na(e_coord))
           return(lof[e_coord])
         #return("lof")
       e_coord <- match(e, miss)
       e_coord <- e_coord[!is.na(e_coord)][1]
       if (!is.na(e_coord))
           return(miss[e_coord])
           #return("miss")
       return(other.name)
      })
}
addScores <- function(mydf, fromfile){
    tempdf <- read.csv(fromfile, header=FALSE, stringsAsFactors=F)
    names(tempdf) <- c("fam","chrom","pos","score")
    tempdf$fam.pos <- do.call(paste, c(tempdf[,c("fam","chrom", "pos")], sep="_"))
    return(rbind(mydf, tempdf))
}
fnShare <- function(x) x/sum(x)

load("denovo.RData")
require(plyr)
install.packages("reshape2")

### use splitstackshape to split annotation into different rows
install.packages("splitstackshape")
install.packages("~/Downloads/splitstackshape_1.4.2.tar.gz", repos = NULL, type="source")
require(splitstackshape)
testdf <- head(sf.cand.p1)
str(testdf[,c("gene", "effect")])
testdf$gene
testdf$effect
testdf$impact
testdf$feature
a <- as.data.frame(cSplit(testdf, c("gene", "effect","impact"), stripWhite=TRUE, direction="long", sep="__"))
i <- sapply(a, is.factor)
i
a[i] <- lapply(a[i], as.character)
str(a)
a

#for full data
names(sf.cand.p1)
sf.cand.p1.long <- as.data.frame(cSplit(sf.cand.p1, splitCols=c("gene", "effect","impact"), makeEqual=TRUE, direction="long", sep="__"))
i <- sapply(sf.cand.p1.long, is.factor)
i
sf.cand.p1.long[i] <- lapply(sf.cand.p1.long[i], as.character)
fnSplitAnno <- function(mydf, mycols, mysep="__"){
    tempdf <- as.data.frame(cSplit(mydf, splitCols=mycols, drop=FALSE, direction="long", sep=mysep))
    i <- sapply(tempdf, is.factor)
    tempdf[i] <- lapply(tempdf[i], as.character)
    tempdf <- tempdf[!is.na(tempdf$gene),]
    return(tempdf)
}
fnSplitAnno(testdf,c("gene", "effect","impact","feature")) 
sf.cand.p1.long <- fnSplitAnno(sf.cand.p1, c("gene", "effect","impact","feature")) 
sf.cand.p1.new.long <- fnSplitAnno(sf.cand.p1.new, c("gene", "effect","impact","feature")) 
sf.cand.p1.new.long$effects <- "other"
sf.cand.p1.new.long$effects[sf.cand.p1.new.long$effect%in%effects.lof] <- "lof"
sf.cand.p1.new.long$effects[sf.cand.p1.new.long$effect%in%effects.miss] <- "mis"
sf.cand.p1.new.long$effects[sf.cand.p1.new.long$effect%in%effects.syn] <- "syn"
sf.cand.p1.new.long$dmg <- "other"
sf.cand.p1.new.long$dmg1 <- "other"
sf.cand.p1.new.long$dmg[(sf.cand.p1.new.long$dbNSFP_Polyphen2_HDIV_pred=="D" | sf.cand.p1.new.long$dbNSFP_Polyphen2_HVAR_pred=="D") & sf.cand.p1.new.long$dbNSFP_SIFT_pred=="D"] <- "damaging"
sf.cand.p1.new.long$dmg1[(sf.cand.p1.new.long$dbNSFP_Polyphen2_HDIV_pred=="D" | sf.cand.p1.new.long$dbNSFP_Polyphen2_HVAR_pred=="D") | sf.cand.p1.new.long$dbNSFP_SIFT_pred=="D"] <- "damaging"
sf.cand.p1.new.long$gene.fam.var <- do.call(paste, c(sf.cand.p1.new.long[,c("gene", "fam.var")], sep="_"))
sf.cand.p1.new.long$gene.fam.var.effect <- do.call(paste, c(sf.cand.p1.new.long[,c("gene.fam.var", "effects")], sep="_"))
#sf.cand.p1.new.long$gene.fam.varsf.cand.p1.new.long <- sf.cand.p1.new.long[!duplicated(sf.cand.p1.new.long$gene.fam.var),]
head(sf.cand.p1.new.long$gene.fam.var)

names(sf.cand.p1.new.long)
dim(sf.cand.p1)
dim(sf.cand.p1.long)
dim(sf.cand.p1.new)
dim(sf.cand.p1.new.long)
str(sf.cand.p1.new.long)
table(sf.cand.p1.new.long$gene)
table(sf.cand.p1.new.long$effect)
table(sf.cand.p1.new.long$effects)
table(sf.cand.p1.new.long$dmg)
length(unique(sf.cand.p1.new.long$gene))
dnm.cut <- .75
c1 <- sf.cand.p1.new.long$score.comb > dnm.cut & !duplicated(sf.cand.p1.new.long$gene.fam.var.effect)
x <- ddply(sf.cand.p1.new.long[c1,], .(gene), summarize,
      n.mis=sum(effects=="mis"),
      n.lof=sum(effects=="lof"),
      n.syn=sum(effects=="syn"),
      n.other=sum(effects=="other")
      )
dim(x)
x <- x[order(x$n.lof, decreasing=T),]
x[x$n.lof > 0,]
x$gene[x$n.lof > 0]
genes.lof <- x$gene[x$n.lof > 0]
tempdf <- sf.cand.p1.new.long[c1,]
lof.mut <- sf.cand.p1.new.long[c1 & sf.cand.p1.new.long$gene%in%genes.lof & sf.cand.p1.new.long$effects=="lof" , 
                    c("fam","chrom", "pos", "ref","alt", "gene","effect")]
genes.lof.function <- list()
genes.lof.function["TBC1D12"] <- "GTPase-activating protein for Rab family protein(s)"
genes.lof.function["MAP3K2"] <- "mitogen-activated protein kinase"
genes.lof.function["COL4A2"] <- "collagen"
genes.lof.function["SSRP1"] <- "structure-specific recognition protein"
genes.lof.function["ESR2"] <- "estrogen receptor"
genes.lof.function["ZNF394"] <- "zinc finger protein"
genes.lof.function["BBS1"] <- "Bardet-Biedl syndrome"
genes.lof.function["APOF"] <- "apolipoprotein F"
genes.lof.function["MGAT3"] <- "mannosyl"
genes.lof.function <- unlist(genes.lof.function)

lof.mut$func <- genes.lof.function[match(lof.mut$gene, names(genes.lof.function))]
lof.mut
write.table(lof.mut, "/mnt/scratch/asalomatov/SSC_rerun/denovo_analysis/lof_mut.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

dnm.cut <- .75
c1 <- sf.cand.p1.new.long$dmg=="damaging" & sf.cand.p1.new.long$score.comb > dnm.cut  & !duplicated(sf.cand.p1.new.long$gene.fam.var.effect)
x <- ddply(sf.cand.p1.new.long[c1,], .(gene), summarize,
      n.mis=sum(effects=="mis"),
      n.lof=sum(effects=="lof"),
      n.syn=sum(effects=="syn"),
      n.other=sum(effects=="other")
      )
dim(x)
x
x <- x[order(x$n.mis, decreasing=T),]
rownames(x) <- NULL
head(x)
x[x$n.mis > 0,]
mis.mut <- x[x$n.mis > 0,]
intersect(lof.mut$gene, mis.mut$gene)
x$gene[x$n.mis > 0]
paste(x$gene[x$n.mis > 0], sep=",", collapse=",")
write.table(x[x$n.mis > 1,], "/mnt/scratch/asalomatov/SSC_rerun/denovo_analysis/num__dmgn_mut_by_gene_dmn75.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
sum(x$n.mis > 1)
write.table(x$gene[x$n.mis > 1], "/mnt/scratch/asalomatov/SSC_rerun/denovo_analysis/num__dmgn_mut_genes_dmn75.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

sf.cand.p1.new.long[with(sf.cand.p1.new.long, gene%in%c("BBS1") & effects=="lof"),  ]
y <- sf.cand.p1.new.long[with(sf.cand.p1.new.long, gene%in%c("GPRC6A","IDH2", "OSTF1","SMTN") & c1), !names(sf.cand.p1.new.long)%in%namesNotToOutput ]
rownames(y) <- NULL
write.table(y, "/mnt/scratch/asalomatov/SSC_rerun/denovo_analysis/mult_dmg_missense.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
t(y)
data.frame(y=y)
ssc.genes <- read.csv("~/Downloads/gene-report.csv", header=TRUE, stringsAsFactors=F)
head(ssc.genes)
intersect(ssc.genes$Gene.symbol, x$gene[x$n.mis > 1])




fam.187 <- strsplit("11006 11009 11013 11023 11029 11043 11055 11056 11064 11069 11083 11093 11096 11109 11120 11124 11141 11148 11172 11184 11190 11193 11195 11198 11205 11218 11224 11229 11246 11257 11262 11303 11346 11375 11388 11390 11396 11398 11425 11452 11459 11469 11471 11472 11479 11480 11491 11498 11504 11506 11510 11518 11523 11526 11545 11556 11569 11587 11599 11610 11629 11653 11659 11660 11691 11696 11707 11711 11715 11722 11734 11753 11773 11788 11827 11834 11843 11863 11872 11928 11942 11947 11948 11959 11964 11969 11989 12011 12015 12036 12073 12086 12106 12114 12118 12130 12152 12153 12157 12161 12185 12198 12212 12225 12233 12237 12238 12249 12285 12296 12300 12304 12335 12341 12346 12358 12373 12381 12390 12430 12437 12444 12521 12532 12555 12578 12603 12621 12630 12641 12667 12674 12703 12741 12744 12752 12810 12905 12933 13008 13031 13048 13116 13158 13169 13188 13207 13274 13314 13333 13335 13346 13415 13447 13494 13517 13530 13532 13533 13557 13593 13606 13610 13629 13668 13678 13701 13733 13741 13742 13757 13793 13812 13815 13822 13835 13844 13857 13863 13890 13914 13926 14006 14011 14020 14201 14292", split=" ")
fam.187  <- unlist(fam.187)
fam.62.quad <- strsplit("13188 14011 11964 13048 11491 13793 11190 13890 13835 12810 12390 13169 12905 11569 11629 11469 12106 11773 13447 12161 13116 11013 11872 11172 11711 11715 12011 14201 12741 11390 11959 13926 13335 11942 13815 12373 12285 13593 12703 11029 11659 11472 11459 11610 11788 13606 11229 13346 11452 11479 11722 13629 12152 12153 12630 12578 11696 12304 13533 12358 12233 11691", split=" ")
fam.62.quad  <- unlist(fam.62.quad)
length(fam.62.quad)
var.187.pm50 <- readLines("/mnt/ceph/asalomatov/data/SSCexome/var_187-pm50.txt")
str(var.187.pm50)
denovo <- read.delim("/mnt/ceph/asalomatov/data/SSCexome/SSC_exome_denovo_Ios_Krumm.csv", stringsAsFactors=F)
str(denovo)
denovo$VARTYPE <- "SNP"
x <- sapply(denovo$Reference, nchar)
table(x)
denovo$VARTYPE[x > 1] <- "DEL"
x <- sapply(denovo$Alternate, nchar)
table(x)
denovo$VARTYPE[x > 1] <- "INS"
table(denovo$VARTYPE)
names(denovo) <- c("chrom","pos","ref","alt","sample","gene","effect","group","Amino.Acid.Change","Exon","Transcript","Validation","Study","VARTYPE")
denovo$fam <- sapply(denovo$sample, function(x) as.integer(strsplit(x, split=".", fixed=TRUE)[[1]][1]))
denovo$memb <- sapply(denovo$sample, function(x) strsplit(x, split=".", fixed=TRUE)[[1]][2])
denovo$var <- do.call(paste, c(denovo[,c("chrom", "pos", "ref", "alt")], sep="_"))
denovo$fam.var <- do.call(paste, c(denovo[,c("fam", "chrom", "pos", "ref", "alt")], sep="_"))
denovo$var.pos <- do.call(paste, c(denovo[,c("chrom", "pos")], sep="_"))
head(denovo)
denovo.187 <- denovo[denovo$fam %in% fam.187 & denovo$memb == "p1",]
denovo.62 <- denovo[denovo$fam %in% fam.62.quad & denovo$memb == "s1",]
denovo.187.pos <- denovo.187[denovo.187$VARTYPE=="SNP" & denovo.187$Validation=="Y", c("fam","chrom","pos", "var")]
denovo.187.pos <- denovo.187.pos[!duplicated(denovo.187.pos[,"var"]),c("fam","chrom","pos")]
dim(denovo.187.pos)
denovo.187.neg <- denovo.187[denovo.187$VARTYPE=="SNP" & denovo.187$Validation=="N", c("fam","chrom","pos", "var")]
denovo.187.neg <- denovo.187.neg[!duplicated(denovo.187.neg[,"var"]),c("fam","chrom","pos")]
dim(denovo.187.neg)
write.table(denovo.187.pos[,c("fam","chrom","pos")], "/mnt/ceph/asalomatov/data/SSCexome/DNM_traning/dnmfilt/denovo-187-pos.csv", quote=FALSE, row.names=FALSE, col.names=FALSE, sep=",")
write.table(denovo.187.neg[,c("fam","chrom","pos")], "/mnt/ceph/asalomatov/data/SSCexome/DNM_traning/dnmfilt/denovo-187-neg.csv", quote=FALSE, row.names=FALSE, col.names=FALSE, sep=",")
denovo.62.pos <- denovo.62[denovo.62$VARTYPE=="SNP" & denovo.62$Validation=="Y", c("fam","chrom","pos", "var")]
denovo.62.pos <- denovo.62.pos[!duplicated(denovo.62.pos[,"var"]),c("fam","chrom","pos")]
dim(denovo.62.pos)
denovo.62.neg <- denovo.62[denovo.62$VARTYPE=="SNP" & denovo.62$Validation=="N", c("fam","chrom","pos", "var")]
denovo.62.neg <- denovo.62.neg[!duplicated(denovo.62.neg[,"var"]),c("fam","chrom","pos")]
dim(denovo.62.neg)
write.table(denovo.62.pos[,c("fam","chrom","pos")], "/mnt/ceph/asalomatov/data/SSCexome/DNM_traning/dnmfilt/denovo-62-pos.csv", quote=FALSE, row.names=FALSE, col.names=FALSE, sep=",")
write.table(denovo.62.neg[,c("fam","chrom","pos")], "/mnt/ceph/asalomatov/data/SSCexome/DNM_traning/dnmfilt/denovo-62-neg.csv", quote=FALSE, row.names=FALSE, col.names=FALSE, sep=",")
?write.csv


denovo.187.pm50 <- denovo.187[denovo.187$var %in% var.187.pm50,]
denovo.62.pm50 <- denovo.62[denovo.62$var %in% var.187.pm50,]
dim(denovo.187.pm50)
dim(denovo.62.pm50)
denovo.187.pm50.uniq <- denovo.187.pm50[!duplicated(denovo.187.pm50[,"fam.var"]),]
denovo.62.pm50.uniq <- denovo.62.pm50[!duplicated(denovo.62.pm50[,"fam.var"]),]
dim(denovo.187.pm50.uniq)
dim(denovo.62.pm50.uniq)
table(denovo.187.pm50$fam)
table(denovo.187.pm50$fam)
table(denovo.187$fam)
table(denovo.187$fam)
dim(denovo.187.pm50.uniq)
table(denovo.187.pm50.uniq$Validation)
table(denovo.187.pm50.uniq$Study)
length(unique(denovo.187.pm50.uniq$fam))
table(denovo.187.pm50.uniq$memb)
ddply(denovo.187.pm50.uniq, .(Study), summarize,
      N = length(var),
      share = round(length(var)/nrow(denovo.187.pm50.uniq), 2)
      )
ddply(denovo.187.pm50.uniq, .(Study, VARTYPE, Validation), summarize,
      N = length(var)
      )
strsplit("10.p1", split=".", fixed=TRUE)[[1]][1]
?strsplit

fam.187  <- unlist(fam.187)
str(fam.187)
sf.cand.p1 <- data.frame(chrom=character(), pos=integer(),id=character(),ref=character(),alt=character(),qual=double(),VARTYPE=character(),effect=character(),
                         impact=character(),gene=character(),feature=character(),gt.fa=character(),gt.mo=character(),gt.p1=character(),dp.fa=integer(),
                         dp.mo=integer(),dp.p1=integer(),GEN=character(),fam=integer(),memb=character(),caller=character(), stringsAsFactors=F)
sf.cand.s1 <- data.frame(chrom=character(), pos=integer(),id=character(),ref=character(),alt=character(),qual=double(),VARTYPE=character(),
                         effect=character(),impact=character(),gene=character(),feature=character(),gt.fa=character(),gt.mo=character(),
                         gt.p1=character(),gt.s1=character(),dp.fa=integer(),dp.mo=integer(),dp.p1=integer(),dp.s1=integer(),GEN=character(),
                         fam=integer(),memb=character(),caller=character(), stringsAsFactors=F)
str(sf.cand.p1)
str(sf.cand.s1)
memb <- "s1"
i <- 0
mpath <- "/mnt/ceph/asalomatov/data/SSCexome/DNM_traning/potential_calls/"
callers <- c("HC","JHC","FB","PL")
readPossibleCallsAllFields <- function(memb, cllr, families, mypath){
    mylist  <- list()
    for (fam in families){
        i <- i + 1
        print("processing family number:")
        print(i)
        mfile <- paste0(mypath,fam,"-",cllr,"-pm50-ann-",memb,"-nsfp",".txt")
        print(mfile)
        tempdf <- read.delim(mfile, stringsAsFactors=F)
        if(cllr == "JHC")
            names(tempdf) <- JHC.fields
        if(cllr == "HC")
            names(tempdf) <- HC.fields
        if(cllr == "FB")
            names(tempdf) <- FB.fields
        if(cllr == "PL")
            names(tempdf) <- PL.fields
        if(nrow(tempdf) == 0) next
        tempdf$fam <- as.integer(fam)
        tempdf$fam.memb <- paste(fam, memb, sep=".")
        tempdf$memb <- memb
        tempdf$caller <- cllr
        tempdf$var <- do.call(paste, c(tempdf[,c("CHROM", "POS", "REF", "ALT")], sep="_"))
        tempdf$var.pos <- do.call(paste, c(tempdf[,c("CHROM", "POS")], sep="_"))
        tempdf$fam.pos <- do.call(paste, c(tempdf[,c("fam", "CHROM", "POS")], sep="_"))
        tempdf$fam.var <- do.call(paste, c(tempdf[,c("fam", "CHROM", "POS", "REF", "ALT")], sep="_"))
        if (is.null(mylist[[as.character(fam)]])) {
            mylist[[as.character(fam)]] <-  tempdf
        }else{
            c1 <- tempdf$var %in% mylist[[as.character(fam)]][["var"]]
            print(paste0("num of new vars : ", sum(!c1)))
            mylist[[as.character(fam)]] <-  rbind.fill(mylist[[as.character(fam)]], tempdf[!c1,])
            mylist[[as.character(fam)]] <- mylist[[as.character(fam)]][!duplicated(mylist[[as.character(fam)]][,"fam.var"]),]
        }
        #c1 <- tempdf$var %in% get(paste0("sf.cand.", memb))
        #assign(paste0("sf.cand.", memb), rbind.fill(get(paste0("sf.cand.", memb)), tempdf[!c1,]))
        print(length(mylist))
        print(dim(mylist[[as.character(fam)]]))
    }
    return(mylist)
}
sf.jhc.p1.list <- readPossibleCallsAllFields("p1", "JHC", fam.187,"/mnt/scratch/asalomatov/SSC_rerun/denovo_analysis/DNM_traning/potential_calls_dp6/dbnsfp/tables/")
sf.jhc.p1 <- rbind.fill(sf.jhc.p1.list)
sf.jhc.s1.list <- readPossibleCallsAllFields("s1", "JHC", fam.62.quad,"/mnt/scratch/asalomatov/SSC_rerun/denovo_analysis/DNM_traning/potential_calls_dp6/dbnsfp/tables/")
sf.jhc.s1 <- rbind.fill(sf.jhc.s1.list)

sf.hc.p1.list <- readPossibleCallsAllFields("p1", "HC", fam.187,"/mnt/scratch/asalomatov/SSC_rerun/denovo_analysis/DNM_traning/potential_calls_dp6/dbnsfp/tables/")
sf.hc.p1 <- rbind.fill(sf.hc.p1.list)
sf.hc.s1.list <- readPossibleCallsAllFields("s1", "HC", fam.62.quad,"/mnt/scratch/asalomatov/SSC_rerun/denovo_analysis/DNM_traning/potential_calls_dp6/dbnsfp/tables/")
sf.hc.s1 <- rbind.fill(sf.hc.s1.list)

sf.fb.p1.list <- readPossibleCallsAllFields("p1", "FB", fam.187,"/mnt/scratch/asalomatov/SSC_rerun/denovo_analysis/DNM_traning/potential_calls_dp6/dbnsfp/tables/")
sf.fb.p1 <- rbind.fill(sf.fb.p1.list)
sf.fb.s1.list <- readPossibleCallsAllFields("s1", "FB", fam.62.quad,"/mnt/scratch/asalomatov/SSC_rerun/denovo_analysis/DNM_traning/potential_calls_dp6/dbnsfp/tables/")
sf.fb.s1 <- rbind.fill(sf.fb.s1.list)

sf.pl.p1.list <- readPossibleCallsAllFields("p1", "PL", fam.187,"/mnt/scratch/asalomatov/SSC_rerun/denovo_analysis/DNM_traning/potential_calls_dp6/dbnsfp/tables/")
sf.pl.p1 <- rbind.fill(sf.pl.p1.list)
sf.pl.s1.list <- readPossibleCallsAllFields("s1", "PL", fam.62.quad,"/mnt/scratch/asalomatov/SSC_rerun/denovo_analysis/DNM_traning/potential_calls_dp6/dbnsfp/tables/")
sf.pl.s1 <- rbind.fill(sf.pl.s1.list)


memb <- "s1"
i <- 0
mylist  <- list()
mpath <- "/mnt/ceph/asalomatov/data/SSCexome/DNM_traning/potential_calls/"
callers <- c("HC","JHC","FB","PL")
readPossibleCalls <- function(memb, families, mypath){
    for (fam in families){
        i <- i + 1
        print("processing family number:")
        print(i)
        for(meth in callers){
            mfile <- paste0(mypath,fam,"-",meth,"-pm50-ann-",memb,".txt")
            print(mfile)
            tempdf <- read.delim(mfile, stringsAsFactors=F)
            names(tempdf) <- c("chrom","pos","id","ref","alt","qual","VARTYPE","effect","impact","gene","feature","gt.fa","gt.mo","gt.p1","dp.fa","dp.mo","dp.p1","GEN")
            if(nrow(tempdf) == 0) next
            tempdf$fam <- as.integer(fam)
            tempdf$memb <- paste(fam, memb, sep=".")
            tempdf$caller <- meth
            tempdf$var <- do.call(paste, c(tempdf[,c("chrom", "pos", "ref", "alt")], sep="_"))
            tempdf$var.pos <- do.call(paste, c(tempdf[,c("chrom", "pos")], sep="_"))
            tempdf$fam.pos <- do.call(paste, c(tempdf[,c("fam", "chrom", "pos")], sep="_"))
            tempdf$fam.var <- do.call(paste, c(tempdf[,c("fam", "chrom", "pos", "ref", "alt")], sep="_"))
            if (is.null(mylist[[as.character(fam)]])) {
                mylist[[as.character(fam)]] <-  tempdf
            }else{
                c1 <- tempdf$var %in% mylist[[as.character(fam)]][["var"]]
                print(paste0("num of new vars : ", sum(!c1)))
                mylist[[as.character(fam)]] <-  rbind.fill(mylist[[as.character(fam)]], tempdf[!c1,])
                mylist[[as.character(fam)]] <- mylist[[as.character(fam)]][!duplicated(mylist[[as.character(fam)]][,"fam.var"]),]
            }
            #c1 <- tempdf$var %in% get(paste0("sf.cand.", memb))
            #assign(paste0("sf.cand.", memb), rbind.fill(get(paste0("sf.cand.", memb)), tempdf[!c1,]))
            print(length(mylist))
            print(dim(mylist[[as.character(fam)]]))
    }
    }
    return(mylist)
}
sf.cand.p1 <- rbind.fill(readPossibleCalls("p1", fam.187, "/mnt/ceph/asalomatov/data/SSCexome/DNM_traning/potential_calls_dp6/"))
#sf.cand.p1$var.pos <- do.call(paste, c(sf.cand.p1[,c("chrom", "pos")], sep="_"))
#sf.cand.p1$fam.pos <- do.call(paste, c(sf.cand.p1[,c("fam","chrom", "pos")], sep="_"))
#sf.cand.p1$fam.var <- do.call(paste, c(sf.cand.p1[,c("fam","chrom", "pos","ref", "alt")], sep="_"))
sf.cand.p1 <- sf.cand.p1[!duplicated(sf.cand.p1[,"fam.var"]),]
sf.cand.p1$dp.fa <- as.integer(sf.cand.p1$dp.fa)
sf.cand.p1$dp.mo <- as.integer(sf.cand.p1$dp.mo)
sf.cand.p1$dp.p1 <- as.integer(sf.cand.p1$dp.p1)
sum(!(is.na(sf.cand.p1$dp.fa) | is.na(sf.cand.p1$dp.mo) | is.na(sf.cand.p1$dp.p1)))
dim(sf.cand.p1)
head(sf.cand.p1)
str(sf.cand.p1)
sf.cand.s1 <- rbind.fill(readPossibleCalls("s1", fam.62.quad, "/mnt/ceph/asalomatov/data/SSCexome/DNM_traning/potential_calls_dp6/"))
#sf.cand.s1$var.pos <- do.call(paste, c(sf.cand.s1[,c("chrom", "pos")], sep="_"))
#sf.cand.s1$fam.pos <- do.call(paste, c(sf.cand.s1[,c("fam","chrom", "pos")], sep="_"))
#sf.cand.s1$fam.var <- do.call(paste, c(sf.cand.s1[,c("fam","chrom", "pos","ref", "alt")], sep="_"))
sf.cand.s1 <- sf.cand.s1[!duplicated(sf.cand.s1[,"fam.var"]),]
dim(sf.cand.s1)
table(sf.cand.s1$dp.p1)
var.sf.cand <- unique(c(sf.cand.p1$var, sf.cand.s1$var))
length(var.sf.cand)
names(sf.cand.p1)
table(sf.cand.p1$gt.p1)
table(sf.cand.p1$gt.fa)
table(sf.cand.p1$gt.mo)
table(sf.cand.p1$caller)
table(sf.cand.p1$caller[sf.cand.p1$dp.fa>30 & sf.cand.p1$dp.mo>30 & sf.cand.p1$dp.p1>30 ])
ddply(sf.cand.p1, .(caller), summarize,
      fa.dp = mean(dp.fa),
      mo.dp = mean(dp.fa),
      p1.dp = mean(dp.fa)
      )


head(sf.cand.p1)
head(sf.cand.p1[sf.cand.p1$caller == 'FB', ])
tail(sf.cand.p1)
dim(sf.cand.p1)
dim(sf.cand.p1)
table(sf.cand.p1$caller)
sf.cand.p1.hc <- sf.cand.p1[sf.cand.p1$caller == "HC",]
dim(sf.cand.p1.hc)
length(unique(sf.cand.p1$var))
#install.packages("dplyr")
#install.packages("~/Downloads/dplyr_0.4.2.tar.gz", repos = NULL, type="source")
sf.cand.p1.hc <- sf.cand.p1.hc[!duplicated(sf.cand.p1.hc[,"var"]),]
dim(sf.cand.p1.hc)
ls()
#

length(denovo.187$var)
length(unique(denovo.187$var))
length(unique(denovo.187$fam))
sum(denovo.187$var %in% sf.cand.p1$var)
sum(denovo.187.pm50.uniq$var %in% sf.cand.p1$var)/length(denovo.187.pm50.uniq$var)
sum(denovo.187.pm50.uniq$var %in% sf.cand.p1$var)
ddply(denovo.187.pm50.uniq, .(memb,Study, Validation), summarize,
      sens = sum(var %in% sf.cand.p1$var)/length(var)
      )
ddply(denovo.187.pm50.uniq, .(memb,Study, Validation), summarize,
      n = sum(var %in% var.sf.cand),
      N = length(var),
      sens = sum(var %in% var.sf.cand)/length(var)
      )
### summary for p1
dnm.cut <- .5
var.type <- "SNP"
min.dp <- 7
known.denovo <- denovo.187.pm50.uniq[denovo.187.pm50.uniq$VARTYPE==var.type,]
#known.denovo <- denovo.187.pm50.uniq
cond1 <-  sf.cand.p1$score.comb > dnm.cut  & sf.cand.p1$id == "." & sf.cand.p1$VARTYPE==var.type
sum(is.na(cond1))
SFcond1<-  sf.cand.p1$score.SFcomb > dnm.cut  & sf.cand.p1$id == "." & sf.cand.p1$VARTYPE==var.type
cond2 <-  sf.cand.p1$dp.fa > min.dp & sf.cand.p1$dp.mo > min.dp & sf.cand.p1$dp.p1 > min.dp
str(sf.cand.p1$dp.fa)
sum(cond2, na.rm=T)
sum(is.na(sf.cand.p1$dp.fa))
ddply(known.denovo, .(VARTYPE, Validation ), summarize,
      N = length(fam.var),
      n = sum(fam.var %in% sf.cand.p1$fam.var[sf.cand.p1$score.comb > dnm.cut & cond2]),
      dif = N - n,
      sens = round(sum(fam.var %in% sf.cand.p1$fam.var[sf.cand.p1$score.comb > dnm.cut & cond2])/length(fam.var),2),
      extra = sum(! sf.cand.p1$fam.var[cond1 & cond2] %in% known.denovo$fam.var),
      n.SF = sum(fam.var %in% sf.cand.p1$fam.var[sf.cand.p1$score.SFcomb > dnm.cut & cond2]),
      N.SF = length(fam.var),
      sens.SF = round(sum(fam.var %in% sf.cand.p1$fam.var[sf.cand.p1$score.SFcomb > dnm.cut & cond2])/length(fam.var),2),
      extra.SF = sum(! unique(sf.cand.p1$fam.var[SFcond1 & cond2]) %in% known.denovo$fam.var),
      not.captured = sum(! fam.var %in% sf.cand.p1$fam.var[cond2])
      )
known.denovo[known.denovo$Validation=="N" & known.denovo$var%in%sf.cand.p1$var[sf.cand.p1$score.comb > .9], ]
known.denovo[known.denovo$Validation=="Y" & known.denovo$var%in%sf.cand.p1$var[sf.cand.p1$score.comb > .9], ]
sf.cand.p1[sf.cand.p1$var %in% known.denovo$var[known.denovo$Validation=="Y"], ]
min(sf.cand.p1$dp.p1[sf.cand.p1$var %in% known.denovo$var[known.denovo$Validation=="Y"]])
table(sf.cand.p1$caller[sf.cand.p1$var %in% known.denovo$var[known.denovo$Validation=="Y"]])
length(sf.cand.p1$fam.var[cond1 & cond2])
length(unique(sf.cand.p1$fam.var[cond1 & cond2]))
sum(!unique(sf.cand.p1$fam.var[cond1 & cond2]) %in% known.denovo$fam.var)

sf.cand.p1.subset <- sf.cand.p1[cond1 & cond2,]
sf.cand.p1.subset <- sf.cand.p1.subset[!is.na(sf.cand.p1.subset$chrom),]
known <- sf.cand.p1.subset$fam.var %in% known.denovo$fam.var
sum(is.na(known))
sum(is.na(sf.cand.p1.subset$chrom))
sf.cand.p1.new <- sf.cand.p1.subset[! known, ]
sf.cand.p1.new <- sf.cand.p1.new[!duplicated(sf.cand.p1.new$fam.var), ]
dim(sf.cand.p1.new)
sum(!duplicated(sf.cand.p1.new$fam.var))
asFnsfpAnn <- function(f.v, data.set, field){
    myind <- match(f.v, data.set$fam.var)
    res <- data.set[myind, field]
    res[res=="."] <- NA
    if(field=="dbNSFP_Polyphen2_HVAR_pred" | field=="dbNSFP_Polyphen2_HDIV_pred" ){
        res1 <- sapply(res, function(x){
               if(grepl("D", x))
                   return("D")
               if(grepl("P", x))
                   return("P")
               if(grepl("B", x))
                   return("B")
               return(NA)
          })
        names(res1) <- NULL
        return(res1)
    }
    if(field=="dbNSFP_SIFT_pred"){
        res1 <- sapply(res, function(x){
               if(grepl("D", x))
                   return("D")
               if(grepl("T", x))
                   return("T")
               return(NA)
          })
        names(res1) <- NULL
        return(res1)
    }
    if(field=="dbNSFP_MutationTaster_pred"){
        res1 <- sapply(res, function(x){
               if(grepl("A", x))
                   return("A")
               if(grepl("D", x))
                   return("D")
               if(grepl("N", x))
                   return("N")
               if(grepl("P", x))
                   return("P")
               return(NA)
          })
        names(res1) <- NULL
        return(res1)
    }
    if(field=="dbNSFP_MutationAssessor_pred"){
        res1 <- sapply(res, function(x){
               if(grepl("H", x))
                   return("H")
               if(grepl("M", x))
                   return("M")
               if(grepl("L", x))
                   return("L")
               if(grepl("N", x))
                   return("N")
               return(NA)
          })
        names(res1) <- NULL
        return(res1)
    }
    if(field=="dbNSFP_PROVEAN_pred"){
        res1 <- sapply(res, function(x){
               if(grepl("D", x))
                   return("D")
               if(grepl("N", x))
                   return("N")
               return(NA)
          })
        names(res1) <- NULL
        return(res1)
    }
    if(field=="dbNSFP_CADD_phred" | grepl("GERP",field) | grepl("1000G",field) | grepl("ExAC",field)){
        res1 <- sapply(res, function(x){
               y <- strsplit(x, split=",")[[1]]
               return(as.double(y[1]))
          })
        names(res1) <- NULL
        return(res1)
    }
    res1 <- sapply(res, function(x){
           y <- strsplit(x, split=",")[[1]]
           return(y[1])
      })
    names(res1) <- NULL
    return(res1)
}
asFnsfpAnn(sf.cand.p1.new$fam.var, sf.jhc.p1,"dbNSFP_Polyphen2_HVAR_pred")     
asFnsfpAnn(sf.cand.p1.new$fam.var, sf.jhc.p1,"dbNSFP_Polyphen2_HDIV_pred")     
asFnsfpAnn(sf.cand.p1.new$fam.var, sf.jhc.p1,"dbNSFP_SIFT_pred")     
asFnsfpAnn(sf.cand.p1.new$fam.var, sf.jhc.p1,"dbNSFP_MutationTaster_pred")     
asFnsfpAnn(sf.cand.p1.new$fam.var, sf.jhc.p1,"dbNSFP_MutationAssessor_pred")     
asFnsfpAnn(sf.cand.p1.new$fam.var, sf.jhc.p1,"dbNSFP_PROVEAN_pred")     
asFnsfpAnn(sf.cand.p1.new$fam.var, sf.jhc.p1,"dbNSFP_CADD_phred")     
asFnsfpAnn(sf.cand.p1.new$fam.var, sf.jhc.p1,"dbNSFP_GERP_RS")     
asFnsfpAnn(sf.cand.p1.new$fam.var, sf.jhc.p1,"dbNSFP_GERP_RS_rankscore")     
asFnsfpAnn(sf.cand.p1.new$fam.var, sf.jhc.p1,"dbNSFP_GERP_NR")     
asFnsfpAnn(sf.cand.p1.new$fam.var, sf.jhc.p1,"dbNSFP_ExAC_AC")     
asFnsfpAnn(sf.cand.p1.new$fam.var, sf.jhc.p1,"dbNSFP_ExAC_AF")     
asFnsfpAnn(sf.cand.p1.new$fam.var, sf.jhc.p1,"dbNSFP_1000Gp3_AC")     
asFnsfpAnn(sf.cand.p1.new$fam.var, sf.jhc.p1,"dbNSFP_1000Gp3_AF")     
asFnsfpAnn(sf.cand.p1.new$fam.var, sf.jhc.p1,"dbNSFP_ref")     
names(sf.jhc.p1)
table(sf.jhc.p1$dbNSFP_GERP_RS)
table(sf.jhc.p1$dbNSFP_1000Gp3_AC)

sf.jhc.p1$dbNSFP_genename[match(sf.cand.p1.new$fam.var, sf.jhc.p1$fam.var)]
match(sf.cand.p1.new$fam.var, sf.hc.p1$fam.var)
match(sf.cand.p1.new$fam.var, sf.hc.p1$fam.var)
asFmergeInd <- function(a,b){
    if(length(a)!=length(b))
        stop("vectors must be of equal length")
    res <- rep(NA, length(a))
    for(i in 1:length(a)){
        if(!is.na(a[i])){
           res[i] <- a[i]
           next
        }
        if(!is.na(b[i])){
           res[i] <- b[i]
           next
        }
       res[i] <- NA
    }
    return(res)
}
asFmergeInd(asFnsfpAnn(sf.cand.p1.new$fam.var, sf.jhc.p1,"dbNSFP_Polyphen2_HVAR_pred"),
                 asFnsfpAnn(sf.cand.p1.new$fam.var, sf.fb.p1,"dbNSFP_Polyphen2_HVAR_pred"))
###add 1000G and ExAc annotations to mail df
for(f in grep("dbNSFP", names(sf.jhc.p1), value=T)){
    print(f)
    a1 <- asFnsfpAnn(sf.jhc.p1$fam.var, sf.jhc.p1,f)     
    a2 <- asFnsfpAnn(sf.jhc$fam.var, sf.hc.p1,f)     
    a3 <- asFnsfpAnn(sf.jhc$fam.var, sf.fb.p1,f)     
    a4 <- asFnsfpAnn(sf.jhc$fam.var, sf.pl.p1,f)     
    i1 <- asFmergeInd(a1,a2)
    i2 <- asFmergeInd(i1,a3)
    sf.cand.p1[[f]] <- asFmergeInd(i2,a4)
}
### Ash's filter
ash.df <- sf.jhc.p1
table(ash.df$VARTYPE)
names(ash.df)
ash.df$dbNSFP_ExAC_AF <- as.double(ash.df$dbNSFP_ExAC_AF)
ash.df <- ash.df[is.na(ash.df$dbNSFP_ExAC_AF) | ash.df$dbNSFP_ExAC_AF <= .001,] 
dim(ash.df)
ash.df$dbNSFP_1000Gp3_AF <- as.double(ash.df$dbNSFP_1000Gp3_AF)
ash.df <- ash.df[is.na(ash.df$dbNSFP_1000Gp3_AF) | ash.df$dbNSFP_1000Gp3_AF <= .001,] 
dim(ash.df)
table(ash.df$VARTYPE)
ash.df <- ash.df[as.double(ash.df$QUAL) >= 30,] 
dim(ash.df)
table(ash.df$VARTYPE)
ash.df$FS <- as.double(ash.df$FS)
ash.df <- ash.df[ash.df$VARTYPE!="SNP" | (ash.df$VARTYPE=="SNP" & ash.df$FS < 45),] 
dim(ash.df)
table(ash.df$VARTYPE)
ash.df$QD <- as.double(ash.df$QD)
ash.df <- ash.df[ash.df$VARTYPE!="SNP" | (ash.df$VARTYPE=="SNP" & ash.df$QD >= 2.5),] 
dim(ash.df)
table(ash.df$VARTYPE)
ash.df <- ash.df[ash.df$VARTYPE=="SNP" | ((ash.df$VARTYPE=="INS" | ash.df$VARTYPE=="DEL") & ash.df$QD >= 1),] 
table(ash.df$VARTYPE)
dim(ash.df)
table(ash.df$VARTYPE)
ash.df <- ash.df[ash.df$VARTYPE!="SNP" | (ash.df$VARTYPE=="SNP" & !(ash.df$FS < 40 & ash.df$FS >= 25 & ash.df$QD < 4 & ash.df$FS >= 2.5)),] 
dim(ash.df)
table(ash.df$VARTYPE)
names(ash.df)
#biased Indels
ash.df <- ash.df[ash.df$VARTYPE=="SNP" | ((ash.df$VARTYPE=="INS" | ash.df$VARTYPE=="DEL") & (ash.df$FS < 25 | ash.df$ReadPosRankSum > -3)),] 
table(ash.df$VARTYPE)
dim(ash.df)
###add dnm scores
ash.df$score.comb <- sf.cand.p1$score.comb[match(ash.df$fam.var, sf.cand.p1$fam.var)]
#summarize
### summary for p1
for (var.type in  c("INS", "DEL")){
var.type <- "SNP"
var.type1 <- "SNP"
min.dp <- 7
dnm.cut <- .25
known.denovo <- denovo.187.pm50.uniq[denovo.187.pm50.uniq$VARTYPE==var.type | denovo.187.pm50.uniq$VARTYPE==var.type1,]
names(known.denovo)
#temp.df <- sf.jhc.p1[sf.jhc.p1$VARTYPE==var.type | sf.jhc.p1$VARTYPE==var.type1,]
temp.df <- ash.df[ash.df$VARTYPE==var.type | ash.df$VARTYPE==var.type1,]
#known.denovo <- denovo.187.pm50.uniq
cond1 <- temp.df$ID == "."# & temp.df$id == "." & temp.df$VARTYPE==var.type
cond2 <-  temp.df$fa.DP > min.dp & temp.df$mo.DP > min.dp & temp.df$ch.DP > min.dp &  temp.df$score.comb > dnm.cut 
Fcond1 <- sf.jhc.p1$VARTYPE==var.type | sf.jhc.p1$VARTYPE==var.type1 #sf.jhc.p1$ID == "." # sf.jhc.p1$score.comb > dnm.cut  & sf.jhc.p1$id == "." & sf.jhc.p1$VARTYPE==var.type
Fcond2 <-  sf.jhc.p1$fa.DP > min.dp & sf.jhc.p1$mo.DP > min.dp & sf.jhc.p1$ch.DP > min.dp
X <- ddply(known.denovo, .(Validation ), summarize,
      N = length(fam.var),
      n = sum(fam.var %in% temp.df$fam.var[cond2]),
      dif = N - n,
      sens = round(sum(fam.var %in% temp.df$fam.var[cond2])/length(fam.var),2),
      extra = sum(! temp.df$fam.var[cond1 & cond2] %in% known.denovo$fam.var),
      not.captured = sum(! fam.var %in% sf.jhc.p1$fam.var[Fcond1 & Fcond2])
      )
X
#sf.jhc.p1.all <- X
#sf.jhc.p1.ash <- X
#sf.jhc.p1.ash1 <- X
#sf.jhc.p1.ash2 <- X

sf.jhc.p1.ash2 <- rbind(sf.jhc.p1.ash2, X)
}
sf.jhc.p1.all
sf.jhc.p1.ash
sf.jhc.p1.ash1
sf.jhc.p1.ash2
write.table(sf.jhc.p1.all, "/mnt/scratch/asalomatov/SSC_rerun/denovo_analysis/sf.jhc.p1.all.csv", quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")
write.table(sf.jhc.p1.ash, "/mnt/scratch/asalomatov/SSC_rerun/denovo_analysis/sf.jhc.p1.ash.csv", quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")
write.table(sf.jhc.p1.ash1, "/mnt/scratch/asalomatov/SSC_rerun/denovo_analysis/sf.jhc.p1.ash1.csv", quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")
write.table(sf.jhc.p1.ash2, "/mnt/scratch/asalomatov/SSC_rerun/denovo_analysis/sf.jhc.p1.ash2.csv", quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")
table(sf.jhc.p1$VARTYPE)

### find validated variant with weak allel ballance

sf.jhc.p1$status <- known.denovo$Validation[match(sf.jhc.p1$fam.var, known.denovo$fam.var)]
str(sf.jhc.p1$status)
table(sf.jhc.p1$ch.AD)
sf.jhc.p1$ch.AD[sf.jhc.p1$status=="Y" & !is.na(sf.jhc.p1$status)]
sf.jhc.p1$fam.var[sf.jhc.p1$ch.AD=="35__7" & sf.jhc.p1$status=="Y" & !is.na(sf.jhc.p1$status)]
sum(sf.jhc.p1$status=="Y", na.rm=T)
table(sf.jhc.p1$ch.AD[sf.jhc.p1$status=="Y"])
x <- sapply(sf.jhc.p1$ch.AD[sf.jhc.p1$status=="Y" & !is.na(sf.jhc.p1$status)], function(x) {
      ad <- strsplit(x, split="__")[[1]]
      print(ad)
      ad <- as.numeric(ad)
      ad[2]/ad[1]
      })
summary(x)
sort(x)


## add dbNSFP annotations to possible new denovos                
for(f in grep("dbNSFP", names(sf.jhc.p1), value=T)){
    print(f)
    a1 <- asFnsfpAnn(sf.cand.p1.new$fam.var, sf.jhc.p1,f)     
    a2 <- asFnsfpAnn(sf.cand.p1.new$fam.var, sf.hc.p1,f)     
    a3 <- asFnsfpAnn(sf.cand.p1.new$fam.var, sf.fb.p1,f)     
    a4 <- asFnsfpAnn(sf.cand.p1.new$fam.var, sf.pl.p1,f)     
    i1 <- asFmergeInd(a1,a2)
    i2 <- asFmergeInd(i1,a3)
    sf.cand.p1.new[[f]] <- asFmergeInd(i2,a4)
}
dim(sf.cand.p1.new)
names(sf.cand.p1.new)

namesNotToOutput <- c("GEN","effect", "impact", "gene", "feature","var","var.pos","fam.pos","fam.var")
head(sf.cand.p1.new[,!names(sf.cand.p1.new)%in%namesNotToOutput])
write.table(sf.cand.p1.new, "/mnt/scratch/asalomatov/SSC_rerun/denovo_analysis/p1_187_extra_dnmcut0.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")


### for s1
dnm.cut <- .75
var.type <- "SNP"
known.denovo <- denovo.62.pm50.uniq[denovo.62.pm50.uniq$VARTYPE==var.type,]
#known.denovo <- denovo.187.pm50.uniq
cond1 <-  sf.cand.s1$score.comb > dnm.cut  & sf.cand.s1$id == "."  & sf.cand.s1$VARTYPE==var.type
ddply(known.denovo, .(VARTYPE, Validation ), summarize,
      n = sum(fam.var %in% sf.cand.s1$fam.var[sf.cand.s1$score.comb > dnm.cut]),
      N = length(fam.var),
      sens = round(sum(fam.var %in% sf.cand.s1$fam.var[sf.cand.s1$score.comb > dnm.cut])/length(fam.var),2),
      Tot = sum(cond1, na.rm=TRUE),
      Tot.u = length(unique(sf.cand.s1$fam.var[cond1])),
      extra = sum(! unique(sf.cand.s1$fam.var[cond1]) %in% known.denovo$fam.var),
      not.captured = sum(! fam.var %in% sf.cand.s1$fam.var)
      )
dim(sf.cand.s1)
sum(is.na(sf.cand.s1$score.comb))


ddply(denovo.187.pm50.uniq, .(VARTYPE, Validation), summarize,
      N = length(var)
      )

c1 <- (denovo.187$var %in% sf.cand.p1$var)
table(sf.cand.p1$VARTYPE)
table(denovo.187.pm50.uniq$Validation)
table(denovo.187.pm50.uniq$memb)
table(denovo.187$Study[!c1])
table(denovo.187$Validation[!c1])
sort(table(denovo.187$fam[!c1]))
length(sort(table(denovo.187$fam[!c1])))
denovo.187[!c1 & denovo.187$memb == "p1",]
denovo.187[c1 & denovo.187$fam == 14011 ,]
denovo.187[c1 & denovo.187$memb == "p1",]
denovo.187$fam %in% fam.187

### load dnm scores

#tr176.p1.dnm <- read.csv("/mnt/ceph/asalomatov/data/SSCexome/DNM_traning/dnmfilt/output/tr176trio/all-p1-DNMFilter.csv", header=FALSE, stringsAsFactors=F)
#names(tr176.p1.dnm) <- c("fam","chrom","pos","score")
#tr176.p1.dnm$fam.pos <- do.call(paste, c(tr176.p1.dnm[,c("fam","chrom", "pos")], sep="_"))
#head(tr176.p1.dnm)
#tr264.p1.dnm <- read.csv("/mnt/ceph/asalomatov/data/SSCexome/DNM_traning/dnmfilt/output/tr264epi/all-p1-DNMFilter.csv", header=FALSE, stringsAsFactors=F)
#names(tr264.p1.dnm) <- c("fam","chrom","pos","score")
#tr264.p1.dnm$fam.pos <- do.call(paste, c(tr264.p1.dnm[,c("fam","chrom", "pos")], sep="_"))
#head(tr264.p1.dnm)

######## p1
### SNPs
comb.p1.snp.dnm <- read.csv("/mnt/ceph/asalomatov/data/SSCexome/DNM_traning/dnmfilt/output/combined/all-p1-snp-DNMFilter.csv", header=FALSE, stringsAsFactors=F)
names(comb.p1.snp.dnm) <- c("fam","chrom","pos","score")
comb.p1.snp.dnm$fam.pos <- do.call(paste, c(comb.p1.snp.dnm[,c("fam","chrom", "pos")], sep="_"))
dim(comb.p1.snp.dnm)
comb.p1.snp.dnm <- addScores(comb.p1.snp.dnm, "/mnt/ceph/asalomatov/data/SSCexome/DNM_traning/dnmfilt_6_9/output/combined/all-p1-snp-DNMFilter.csv")
dim(comb.p1.snp.dnm)

SFcomb.p1.snp.dnm <- read.csv("/mnt/ceph/asalomatov/data/SSCexome/DNM_traning/dnmfilt/output/SSC_comb/all-p1-snp-DNMFilter.csv", header=FALSE, stringsAsFactors=F)
names(SFcomb.p1.snp.dnm) <- c("fam","chrom","pos","score")
SFcomb.p1.snp.dnm$fam.pos <- do.call(paste, c(SFcomb.p1.snp.dnm[,c("fam","chrom", "pos")], sep="_"))
dim(SFcomb.p1.snp.dnm)
SFcomb.p1.snp.dnm <- addScores(SFcomb.p1.snp.dnm, "/mnt/ceph/asalomatov/data/SSCexome/DNM_traning/dnmfilt_6_9/output/SSC_comb/all-p1-snp-DNMFilter.csv")
dim(SFcomb.p1.snp.dnm)

### INS and DEL
comb.p1.ins.dnm <- read.csv("/mnt/ceph/asalomatov/data/SSCexome/DNM_traning/dnmfilt/output/combined/all-p1-ins-DNMFilter.csv", header=FALSE, stringsAsFactors=F)
names(comb.p1.ins.dnm) <- c("fam","chrom","pos","dnm.score")
comb.p1.ins.dnm$pos <- comb.p1.ins.dnm$pos - 1
comb.p1.ins.dnm$st <- comb.p1.ins.dnm$pos
comb.p1.ins.dnm$st[-1] <- diff(comb.p1.ins.dnm$pos)
sum(comb.p1.ins.dnm$st == 1)
comb.p1.ins.dnm$st[comb.p1.ins.dnm$st != 1] <- comb.p1.ins.dnm$pos[comb.p1.ins.dnm$st != 1]
st.1 <- which(comb.p1.ins.dnm$st == 1)
x <- st.1*0
j <- 0
for (i in st.1){
    j <- j + 1
    if (comb.p1.ins.dnm$st[i-1] != 1)
        x[j] <- comb.p1.ins.dnm$st[i-1]
    else
        x[j] <- x[j-1]
}
comb.p1.ins.dnm$st[st.1] <- x
comb.p1.ins.dnm$fam.pos <- do.call(paste, c(comb.p1.ins.dnm[,c("fam","chrom", "st")], sep="_"))
head(st.1,30)
head(comb.p1.ins.dnm,40)
comb.p1.ins.dnm.aggr <- ddply(comb.p1.ins.dnm, .(fam.pos), summarize,
                              fam=fam[1],
                              chrom=chrom[1],
                              pos=pos[1],
                              score=mean(dnm.score),
                              score.min=min(dnm.score),
                              score.max=max(dnm.score)
)
dim(comb.p1.ins.dnm.aggr)
head(comb.p1.ins.dnm.aggr)

comb.p1.del.dnm <- read.csv("/mnt/ceph/asalomatov/data/SSCexome/DNM_traning/dnmfilt/output/combined/all-p1-del-DNMFilter.csv", header=FALSE, stringsAsFactors=F)
names(comb.p1.del.dnm) <- c("fam","chrom","pos","dnm.score")
comb.p1.del.dnm$pos <- comb.p1.del.dnm$pos - 1
comb.p1.del.dnm$st <- comb.p1.del.dnm$pos
comb.p1.del.dnm$st[-1] <- diff(comb.p1.del.dnm$pos)
sum(comb.p1.del.dnm$st == 1)
comb.p1.del.dnm$st[comb.p1.del.dnm$st != 1] <- comb.p1.del.dnm$pos[comb.p1.del.dnm$st != 1]
st.1 <- which(comb.p1.del.dnm$st == 1)
x <- st.1*0
j <- 0
for (i in st.1){
    j <- j + 1
    if (comb.p1.del.dnm$st[i-1] != 1)
        x[j] <- comb.p1.del.dnm$st[i-1]
    else
        x[j] <- x[j-1]
}
comb.p1.del.dnm$st[st.1] <- x
comb.p1.del.dnm$fam.pos <- do.call(paste, c(comb.p1.del.dnm[,c("fam","chrom", "st")], sep="_"))
head(st.1,30)
head(comb.p1.del.dnm,40)
comb.p1.del.dnm.aggr <- ddply(comb.p1.del.dnm, .(fam.pos), summarize,
                              fam=fam[1],
                              chrom=chrom[1],
                              pos=pos[1],
                              score=mean(dnm.score),
                              score.min=min(dnm.score),
                              score.max=max(dnm.score)
)
dim(comb.p1.del.dnm.aggr)
head(comb.p1.del.dnm.aggr)

###### the same for s1
comb.s1.snp.dnm <- read.csv("/mnt/ceph/asalomatov/data/SSCexome/DNM_traning/dnmfilt/output/combined/all-s1-snp-DNMFilter.csv", header=FALSE, stringsAsFactors=F)
names(comb.s1.snp.dnm) <- c("fam","chrom","pos","score")
comb.s1.snp.dnm$fam.pos <- do.call(paste, c(comb.s1.snp.dnm[,c("fam","chrom", "pos")], sep="_"))
head(comb.s1.snp.dnm)
dim(comb.s1.snp.dnm)
comb.s1.snp.dnm <- addScores(comb.s1.snp.dnm, "/mnt/ceph/asalomatov/data/SSCexome/DNM_traning/dnmfilt_6_9/output/combined/all-s1-snp-DNMFilter.csv")

SFcomb.s1.snp.dnm <- read.csv("/mnt/ceph/asalomatov/data/SSCexome/DNM_traning/dnmfilt/output/SSC_comb/all-s1-snp-DNMFilter.csv", header=FALSE, stringsAsFactors=F)
names(SFcomb.s1.snp.dnm) <- c("fam","chrom","pos","score")
SFcomb.s1.snp.dnm$fam.pos <- do.call(paste, c(SFcomb.s1.snp.dnm[,c("fam","chrom", "pos")], sep="_"))
dim(SFcomb.s1.snp.dnm)
SFcomb.s1.snp.dnm <- addScores(SFcomb.s1.snp.dnm, "/mnt/ceph/asalomatov/data/SSCexome/DNM_traning/dnmfilt_6_9/output/SSC_comb/all-s1-snp-DNMFilter.csv")
dim(SFcomb.s1.snp.dnm)


### INS and DEL
comb.s1.ins.dnm <- read.csv("/mnt/ceph/asalomatov/data/SSCexome/DNM_traning/dnmfilt/output/combined/all-s1-ins-DNMFilter.csv", header=FALSE, stringsAsFactors=F)
names(comb.s1.ins.dnm) <- c("fam","chrom","pos","dnm.score")
comb.s1.ins.dnm$pos <- comb.s1.ins.dnm$pos - 1
comb.s1.ins.dnm$st <- comb.s1.ins.dnm$pos
comb.s1.ins.dnm$st[-1] <- diff(comb.s1.ins.dnm$pos)
sum(comb.s1.ins.dnm$st == 1)
comb.s1.ins.dnm$st[comb.s1.ins.dnm$st != 1] <- comb.s1.ins.dnm$pos[comb.s1.ins.dnm$st != 1]
st.1 <- which(comb.s1.ins.dnm$st == 1)
x <- st.1*0
j <- 0
for (i in st.1){
    j <- j + 1
    if (comb.s1.ins.dnm$st[i-1] != 1)
        x[j] <- comb.s1.ins.dnm$st[i-1]
    else
        x[j] <- x[j-1]
}
comb.s1.ins.dnm$st[st.1] <- x
comb.s1.ins.dnm$fam.pos <- do.call(paste, c(comb.s1.ins.dnm[,c("fam","chrom", "st")], sep="_"))
head(st.1,30)
head(comb.s1.ins.dnm,40)
comb.s1.ins.dnm.aggr <- ddply(comb.s1.ins.dnm, .(fam.pos), summarize,
                              fam=fam[1],
                              chrom=chrom[1],
                              pos=pos[1],
                              score=mean(dnm.score),
                              score.min=min(dnm.score),
                              score.max=max(dnm.score)
)
dim(comb.s1.ins.dnm.aggr)
head(comb.s1.ins.dnm.aggr)

comb.s1.del.dnm <- read.csv("/mnt/ceph/asalomatov/data/SSCexome/DNM_traning/dnmfilt/output/combined/all-s1-del-DNMFilter.csv", header=FALSE, stringsAsFactors=F)
names(comb.s1.del.dnm) <- c("fam","chrom","pos","dnm.score")
comb.s1.del.dnm$pos <- comb.s1.del.dnm$pos - 1
comb.s1.del.dnm$st <- comb.s1.del.dnm$pos
comb.s1.del.dnm$st[-1] <- diff(comb.s1.del.dnm$pos)
sum(comb.s1.del.dnm$st == 1)
comb.s1.del.dnm$st[comb.s1.del.dnm$st != 1] <- comb.s1.del.dnm$pos[comb.s1.del.dnm$st != 1]
st.1 <- which(comb.s1.del.dnm$st == 1)
x <- st.1*0
j <- 0
for (i in st.1){
    j <- j + 1
    if (comb.s1.del.dnm$st[i-1] != 1)
        x[j] <- comb.s1.del.dnm$st[i-1]
    else
        x[j] <- x[j-1]
}
comb.s1.del.dnm$st[st.1] <- x
comb.s1.del.dnm$fam.pos <- do.call(paste, c(comb.s1.del.dnm[,c("fam","chrom", "st")], sep="_"))
head(st.1,30)
head(comb.s1.del.dnm,40)
comb.s1.del.dnm.aggr <- ddply(comb.s1.del.dnm, .(fam.pos), summarize,
                              fam=fam[1],
                              chrom=chrom[1],
                              pos=pos[1],
                              score=mean(dnm.score),
                              score.min=min(dnm.score),
                              score.max=max(dnm.score)
)
dim(comb.s1.del.dnm.aggr)
head(comb.s1.del.dnm.aggr)
intersect(comb.s1.del.dnm.aggr$fam.pos, comb.s1.ins.dnm.aggr$fam.pos)
intersect(comb.s1.del.dnm.aggr$fam.pos, comb.s1.snp.dnm$fam.pos)
intersect(comb.p1.del.dnm.aggr$fam.pos, comb.p1.ins.dnm.aggr$fam.pos)
intersect(comb.p1.del.dnm.aggr$fam.pos, comb.p1.snp.dnm$fam.pos)
intersect(comb.p1.del.dnm$fam.pos, comb.s1.snp.dnm$fam.pos)
head(comb.p1.del.dnm)



##### annotate variants with scores

sf.cand.p1$score.comb <- 0
c1 <- sf.cand.p1$VARTYPE == "SNP"
sum(c1)
match.pos <- match(sf.cand.p1$fam.pos[c1], comb.p1.snp.dnm$fam.pos)
sf.cand.p1$score.comb[c1] <- comb.p1.snp.dnm$score[match.pos]
sum(is.na(sf.cand.p1$score.comb[c1]))
head(sf.cand.p1[!c1,])

sf.cand.p1$score.SFcomb <- 0
c1 <- sf.cand.p1$VARTYPE == "SNP"
sum(c1)
match.pos <- match(sf.cand.p1$fam.pos[c1], SFcomb.p1.snp.dnm$fam.pos)
sf.cand.p1$score.SFcomb[c1] <- SFcomb.p1.snp.dnm$score[match.pos]
sum(is.na(sf.cand.p1$score.SFcomb[c1]))
head(sf.cand.p1[!c1,])

c1 <- sf.cand.p1$VARTYPE == "INS"
sum(c1)
match.pos <- match(sf.cand.p1$fam.pos[c1], comb.p1.ins.dnm.aggr$fam.pos)
sf.cand.p1$score.comb[c1] <- comb.p1.ins.dnm.aggr$score[match.pos]
sum(is.na(sf.cand.p1$score.comb[c1]))
head(sf.cand.p1[!c1,])
c1 <- sf.cand.p1$VARTYPE == "DEL"
sum(c1)
match.pos <- match(sf.cand.p1$fam.pos[c1], comb.p1.del.dnm.aggr$fam.pos)
sf.cand.p1$score.comb[c1] <- comb.p1.del.dnm.aggr$score[match.pos]
sum(is.na(sf.cand.p1$score.comb[c1]))
head(sf.cand.p1[!c1,])


sf.cand.s1$score.comb <- 0
c1 <- sf.cand.s1$VARTYPE == "SNP"
sum(c1)
match.pos <- match(sf.cand.s1$fam.pos[c1], comb.s1.snp.dnm$fam.pos)
sf.cand.s1$score.comb[c1] <- comb.s1.snp.dnm$score[match.pos]
sum(is.na(sf.cand.s1$score.comb[c1]))
sum(sf.cand.s1$score.comb[c1]==0)
head(sf.cand.s1[c1,])
head(sf.cand.s1[!c1,])

sf.cand.s1$score.SFcomb <- 0
c1 <- sf.cand.s1$VARTYPE == "SNP"
sum(c1)
match.pos <- match(sf.cand.s1$fam.pos[c1], SFcomb.s1.snp.dnm$fam.pos)
sf.cand.s1$score.SFcomb[c1] <- SFcomb.s1.snp.dnm$score[match.pos]
sum(is.na(sf.cand.s1$score.SFcomb[c1]))
head(sf.cand.s1[!c1,])

c1 <- sf.cand.s1$VARTYPE == "INS"
sum(c1)
match.pos <- match(sf.cand.s1$fam.pos[c1], comb.s1.ins.dnm.aggr$fam.pos)
sf.cand.s1$score.comb[c1] <- comb.s1.ins.dnm.aggr$score[match.pos]
sum(is.na(sf.cand.s1$score.comb[c1]))
head(sf.cand.s1[!c1,])
c1 <- sf.cand.s1$VARTYPE == "DEL"
sum(c1)
match.pos <- match(sf.cand.s1$fam.pos[c1], comb.s1.del.dnm.aggr$fam.pos)
sf.cand.s1$score.comb[c1] <- comb.s1.del.dnm.aggr$score[match.pos]
sum(is.na(sf.cand.s1$score.comb[c1]))
head(sf.cand.s1[!c1,])




