setwd("~/Documents/MeisterLab/dSMF/tissueExpn/scripts")

tissues<-readRDS("tissueExpressionBreadth.RDS")
genes<-read.table("newGenes.txt",stringsAsFactors=F,header=F)

getTissueBreadth<-function(genes,tissues){
   i<-match(genes$V1,tissues$WBgeneID)
   genes<-tissues[i,c("WBgeneID","allTissues","larvalTissues","fromUnfilt")]
   return(genes)
}

genesTissueBreadth<-getTissueBreadth(genes,tissues)

write.table(genesTissueBreadth,"newGenes_tissues.txt",row.names=F,quote=F)
