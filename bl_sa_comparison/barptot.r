require(ggplot2)
genes <- read.csv("nov_10_nyg_num_genes_lt15.txt", sep="\t", header=FALSE, stringsAsFactors=F)
genes.bay <- read.csv("nov_10_baylor_num_genes_lt15.txt", sep="\t", header=FALSE, stringsAsFactors=F)
head(genes)
names(genes) <- c("sample", "count")
names(genes.bay) <- c("sample", "count")
head(genes.bay)
genes.bay$sample <- sapply(genes.bay$sample, function(x){
           a <- strsplit(x,"-")[[1]]
           res <- c(a[1], tolower(a[3]), tolower(a[2]))
           print(res)
           return(paste(res, collapse="-"))
})
genes$center <- "NYG"
genes.bay$center <- "Baylor"
genes <- rbind(genes, genes.bay)
genes <- genes[with(genes, order(center, sample)),]
genes <- genes[genes$count != 0,]
genes$id <- paste(genes$center, genes$sample, sep="-")
genes$protocol <- sapply(genes$sample, function(x){
           a <- strsplit(x,"-")[[1]]
           res <-  tolower(a[2])
           return(res)
})
genes$family <- sapply(genes$sample, function(x){
           a <- strsplit(x,"-")[[1]]
           return(a[1])
})
genes$memb <- sapply(genes$sample, function(x){
           a <- strsplit(x,"-")[[1]]
           return(a[3])
})
genes$center_protocol <- paste(genes$center, genes$protocol, sep="-")
genes <- genes[genes$family%in%c("1","2","4"),]
genes <- genes[genes$memb%in%c("proband", "son1", "father","mother"),]
genes
png("~/nov10_genes_barchart.png")
ggplot(data=genes, aes(x=id, y=count, fill=center_protocol)) + 
    geom_bar(color="black", stat="identity") + theme(axis.ticks = element_blank(), axis.text.x = element_blank())
dev.off()


read.csv?
