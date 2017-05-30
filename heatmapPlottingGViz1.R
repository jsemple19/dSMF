library(Gviz)
source("~/Documents/MeisterLab/dSMF/usefulFunctions.R")


designT=read.csv('~/Documents/MeisterLab/dSMF/PromoterPrimerDesign/scripts/finalChosenList_tss_dc_tissues.csv',
                 stringsAsFactors=FALSE,header=TRUE)
design=import.bed('~/Documents/MeisterLab/dSMF/PromoterPrimerDesign/scripts/amplicons_stranded.bed')
tss=import.bed('~/Documents/MeisterLab/dSMF/PromoterPrimerDesign/scripts/ampliconTSS.bed')

#most5p=import.bed('~/Documents/MeisterLab/dSMF/PromoterPrimerDesign/scripts/most5p.bed')
#start(most5p)<-end(most5p)

amplicons=design
strand(amplicons)<-strand(tss)
amplicons=resize(amplicons,500,fix='center')

TSS=tss
TSS2=resize(TSS,500,fix='center')
strand(TSS2)='+'

#meth_gr<-readRDS(paste0(path,'/Amplicon_raw_methylation_PolII_covered_amplicons.rds'))

#cO=10
#PolII_in=call_context_methylation(meth_gr,cO,Celegans)

#PolII_inm=unlist(GRangesList(PolII_in))
# saveRDS(PolII_inm,'/work2/gschub/arnaud/nomeseq/WG/DM/amplicon/methylation_calls/PolII_inhib_av10.rds')
#saveRDS(PolII_inm,paste0(path,'PolII_inhib_av10.rds'))
#naix=apply(as.matrix(elementMetadata(PolII_inm[,1])),1,function(x){sum(is.na(x))==1})
#PolII_inm<-PolII_inm[!naix]
#saveRDS(PolII_inm,paste0(path,'PolII_inhib_av10.rds'))

PolII_inm<-readRDS(paste0(path,'PolII_inhib_av10.rds'))

amplicon<-amplicons[84]

getMethMatrix<-function(amplicons,designT,ampNum,proj,sampleName,genome=Celegans) {
   amplicon<-amplicons[ampNum]
   geneName<-mcols(amplicon)$name
   chr<-strsplit(designT[designT$Gene_WB_ID==geneName,]$Amplicon,":")[[1]][1]
   tss<-designT[designT$Gene_WB_ID==geneName,]$maxTSS
   allCs<-getCMethMatrix(proj,amplicon,sampleName)
   onlyCG_GC<-getGCMatrix(allCs,chr=chr,genome=genome,conv.rate=80,destrand=FALSE)
   mCG_GC<-mergeGC_CGmat(onlyCG_GC)
   #reads<-dim(mCG_GC)[1]
   mCG_GC.o<-mCG_GC[hclust(dist(mCG_GC))$order,]
   return(mCG_GC.o)
}


MethMat2gr<-function(mat,amplicon) {
   Cpos<-as.numeric(colnames(mat))
   gr<-GRanges(seqnames=seqnames(amplicon),ranges=IRanges(start=Cpos,width=1),strand=strand(amplicon))
   mcols(gr)<-DataFrame(t(mat))
   return(gr)
}

i=40
a<-getMethMatrix(amplicons, designT,i,NOMEproj,"N2_DE_ampV001")
amplicon.wb<-ucsc2wb(amplicon)
agr<-MethMat2gr(a,amplicon.wb)
#image(t(a),axes=FALSE,col=grey(seq(0,0.8,length=2)),ylab="molecule",
#      xlab="methylation site")


library(Gviz)
# #plotting heatmaps in Gviz
dtrack<-DataTrack(range=agr,type=c("heatmap"), gradient=c("grey","black"), ncolor=3)
# plotTracks(list(gtrack,atrack,dtrack1,dtrack3,dtrack5,dtrack7),chromosome=seqnames(maxol[i]),
#            from=start(maxol[i]),to=end(maxol[i]))
gtrack<-GenomeAxisTrack()
options(ucscChromosomeNames=FALSE)
amplicon<-amplicons[i]
amplicon.wb<-ucsc2wb(amplicon)
st=ucsc2wb(tss[i])
atrack<-AnnotationTrack(range=st)

# txdb<-makeTxDbFromGFF("../../../GenomeVer/annotations/c_elegans.PRJNA13758.WS250.annotations.gff3")
txTr<-GeneRegionTrack(txdb,showFeatureId=TRUE,name="txdb")

pdf("./gvizMethPlots",paper="a4r",height=8,width=11)
plotTracks(list(gtrack,dtrack),
            chromosome=seqnames(amplicon.wb),from=start(amplicon.wb),to=end(amplicon.wb))

dev.off()
# error message:
#Error in .local(GdObject, ...) :
#   Too many stacks to draw. Either increase the device size or limit the drawing to a smaller region.

#need to plot Arnaud's way
