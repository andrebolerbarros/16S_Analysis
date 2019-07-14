library(phyloseq)
library(dplyr)
library(tidyr)
###################################
dir.create(path = "variables/", showWarnings = FALSE)

###################################
otu_col <- read.table("otu_table.txt",nrow=1,comment.char = "",skip=1,stringsAsFactors = F)
otu_col<-otu_col[-c(1,2)]

otu<-read.table("otu_table.txt",header=F,row.names = 1)
colnames(otu)<-otu_col

tax	<- read.table("taxonomy.tsv",sep = '\t', header = F,row.names = 1,skip=1,stringsAsFactors = F)

tab<-data.frame(Taxon=tax$V2)

tab<-tab %>% separate(Taxon, c("Kingdom",	"Phylum",	"Class",	"Order", "Family", "Genus",	"Species",	"L7",	"L8",	"L9",	"L10",	"L11",	"L12",	"L13",	"L14"),sep = ";")

tax<-data.frame(tab,Confidence=tax$V3,row.names = rownames(tax))
tax[,1]<-gsub("D_[0-9]__","",tax[,1])

for (i in 2:7) {
  j<-i-1
  tax[,i]<-gsub("D_[0-9]__","",tax[,i])
  tax[is.na(tax[,i]),i]<-paste0(tax[is.na(tax[,i]),j],"_unassigned",sep="")
  tax[tax[,i]=="uncultured",i]<-paste0(tax[tax[,i]=="uncultured",j],"_uncultured",sep="")
}

all(rownames(otu)==rownames(tax))

tax<-as.matrix.data.frame(tax)
phy_tree <- read_tree("tree.nwk")

meta_col<-read.table("metadata.tsv",sep="\t",nrow=1,comment.char = "",stringsAsFactors = F)
meta_col<-meta_col[-1]

meta<-read.table("metadata.tsv",header=F,row.names = 1)
colnames(meta)<-meta_col

all(colnames(otu)==rownames(meta))

meta<-meta[order(match(rownames(meta),colnames(otu))),]
all(colnames(otu)==rownames(meta))

OTU <- otu_table(otu, taxa_are_rows = TRUE)
TAX <- tax_table(tax)
META <- sample_data(meta)

ps <- phyloseq(OTU, TAX, META, phy_tree)

saveRDS(ps, "variables/ps.rds") #remember that you need to use readRDS function to read this file
rm(list=ls())

ps<-readRDS("variables/ps.rds")
