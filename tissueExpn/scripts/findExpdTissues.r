# tissue specific expression from modencode. Data describe in this paper:
#Genome Res. (2011) 21(2):325-41.
#A spatial and temporal map of C. elegans gene expression.
#Spencer WC .... Miller DM 3rd
# Gene lists (rather than raw data) downloaded from
#http://www.vanderbilt.edu/wormmap/Core_enriched_genes/ (link is burried in supplementary text)

setwd("~/Documents/MeisterLab/dSMF/tissueExpn/scripts")

library(rtracklayer)
library(GenomicFeatures)

dataFiles<-list.files("../Spencer2011/",pattern="*.txt",recursive=TRUE)
#take tissue data that has been filtered by whole animal reference level
dataFiles<-dataFiles[grep("expressedGenes_WS199",dataFiles)]
#system("ls ")
#list.dirs("../Spencer2011")
geneList=list()
for (t in seq_along(dataFiles)){
   tissueName<-strsplit(strsplit(dataFiles[t],"/",fixed=TRUE)[[1]][2],"_ref")[[1]][1]
   geneList[tissueName]<-read.table(paste0("../Spencer2011/",dataFiles[t]),stringsAsFactors=FALSE)
}

allGenes<-unique(unlist(geneList))

expnTissues<-data.frame(WBgeneID=allGenes,stringsAsFactors=FALSE)
for (t in seq_along(geneList)){
   i<-match(geneList[[t]],expnTissues$WBgeneID)
   expnTissues[i,names(geneList[t])]<-1
}
#15559 genes

numSamples<-dim(expnTissues)[2]

expnTissues$allTissues<-rowSums(expnTissues[2:numSamples],na.rm=TRUE)

L3tissues<-grep("L3",names(expnTissues))
L2tissues<-grep("L2",names(expnTissues))
larval<-c(L2tissues,L3tissues)
Adtissues<-grep("YA",names(expnTissues))
EembTissues<-grep("EE",names(expnTissues))
LembTissues<-grep("LE",names(expnTissues))

expnTissues$larvalTissues<-rowSums(expnTissues[c(2:numSamples)[larval-1]],na.rm=TRUE)

hist(expnTissues$allTissues,breaks=25,col="grey")
hist(expnTissues$larvalTissues,breaks=11,col="grey")

plot(jitter(expnTissues$allTissues),jitter(expnTissues$larvalTissues),pch=16,col="#0000ff33")


###########3### what about house keeping genes?############################
dataFiles<-list.files("../Spencer2011/",pattern="*.txt",recursive=TRUE)
#take tissue data that has been filtered by whole animal reference level
dataFiles<-dataFiles[grep("expressedGenes_unfiltered",dataFiles)]
#system("ls ")
#list.dirs("../Spencer2011")
geneList=list()
for (t in seq_along(dataFiles)){
   tissueName<-strsplit(strsplit(dataFiles[t],"/",fixed=TRUE)[[1]][2],".txt")[[1]][1]
   geneList[tissueName]<-read.table(paste0("../Spencer2011/",dataFiles[t]),stringsAsFactors=FALSE)
}

allGenes<-unique(unlist(geneList))

expnHK<-data.frame(WBgeneID=allGenes,stringsAsFactors=FALSE)
for (t in seq_along(geneList)){
   i<-match(geneList[[t]],expnHK$WBgeneID)
   expnHK[i,names(geneList[t])]<-1
}

numSamples<-dim(expnHK)[2]

#define as housekeeping genes that are expressed in every single tissue in the unfiltered set.
expnHK$allTissues<-rowSums(expnHK[2:numSamples],na.rm=TRUE)
i<-which(expnHK$allTissues==(numSamples-1))
HKgenes<-expnHK$WBgeneID[i]
#remove genes already in the filtered dataset above
i<-which(HKgenes %in% expnTissues$WBgeneID)
HKgenes<-HKgenes[-i]
#4353 genes

#remove reference columns
i<-grep("reference",names(expnHK))
expnHK<-expnHK[,-i]

#subset to only housekeeping genes
i<-which(expnHK$WBgeneID %in% HKgenes)
expnHK<-expnHK[i,]

#sample(expnHK$WBgeneID,10)
#looked at 10 manually. two were very HK (e.g. WBGene00015156), many of others had much weaker expn and 
#expn in particular stage more than other (some dauer, some emb, some L3). but the fact that it caught some
#genuinely widely expressed genes shows that housekeeping genes are no properly included

numSamples=26
#redo allTissue and larvalTissues expn counts
expnHK$allTissues<-rowSums(expnHK[2:numSamples],na.rm=TRUE)
expnHK$larvalTissues<-rowSums(expnHK[c(2:numSamples)[larval-1]],na.rm=TRUE)

expnBreadth<-rbind(expnTissues,expnHK)
sum(duplicated(expnBreadth$WBgeneID))
dim(expnBreadth)
#19912 genes