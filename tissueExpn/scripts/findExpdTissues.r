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

expnTissues<-data.frame(WBgeneID=allGenes)
for (t in seq_along(geneList)){
   i<-match(geneList[[t]],expnTissues$WBgeneID)
   expnTissues[i,names(geneList[t])]<-1
}

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

