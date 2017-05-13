library(BSgenome.Celegans.UCSC.ce11)
library(Biostrings)
genome <- BSgenome.Celegans.UCSC.ce11

getMedianGap<-function(gr,genome) {
   gr<-sortSeqlevels(gr)
   gr<-gr[strand(gr)=="+"]
   gr<-sort(gr)
   grdist<-gaps(gr)
   medianGap<-median(width(grdist))
   print(medianGap)
   return(medianGap)
}

genome <- BSgenome.Celegans.UCSC.ce11
gpc<-vmatchPattern("GC",genome)
cpg<-vmatchPattern("CG",genome)
ol<-findOverlaps(gpc,cpg)
cggc<-c(gpc[-queryHits(ol)],cpg)

getMedianGap(gpc,genome)
#17
getMedianGap(cpg,genome)
#18
getMedianGap(cggc,genome)
#10


#countCGGC(genome,6)
# "There is 1 CpG every 32 bp in the genome"
# "There is 1 GpC every 30 bp in the genome"
# "A CpG or Gpc occur every 17 bp in the genome"

#countCGGCds(genome,6)
#"There is 1 CpG every 16 bp in the genome"
#"There is 1 GpC every 15 bp in the genome"
#"A CpG or Gpc occur every 8 bp in the genome"

# countCGGC<-function(genome,chrNum) {
#    CpG="CG"
#    GpC="GC"
# 
#    CGcounts=vcountPattern(CpG,genome)
#    CGtotal=sum(CGcounts[CGcounts$strand=="+","count"][1:chrNum])
#    GCcounts=vcountPattern(GpC,genome)
#    GCtotal=sum(GCcounts[GCcounts$strand=="+","count"][1:chrNum])
# 
#    #remove duplicate counts of Cs between two Gs
#    GCGcounts=vcountPattern("GCG",genome)
#    GCGtotal=sum(GCGcounts[GCGcounts$strand=="+","count"][1:chrNum])
# 
#    lengthTotal=sum(seqlengths(genome)[1:chrNum])
# 
#    print(paste("There is 1 CpG every",round(lengthTotal/CGtotal,0),"bp in the genome"))
#    print(paste("There is 1 GpC every",round(lengthTotal/GCtotal,0),"bp in the genome"))
#    print(paste("A CpG or GpC occur every",round(lengthTotal/(CGtotal+GCtotal-GCGtotal),0), "bp in the genome"))
# }
# 
# countCGGCds<-function(genome,chrNum) {
#    CpG="CG"
#    GpC="GC"
#    
#    CGcounts=vcountPattern(CpG,genome)
#    CGtotal=2*sum(CGcounts[CGcounts$strand=="+","count"][1:chrNum])
#    GCcounts=vcountPattern(GpC,genome)
#    GCtotal=2*sum(GCcounts[GCcounts$strand=="+","count"][1:chrNum])
#    
#    #remove duplicate counts of Cs between two Gs
#    GCGcounts=vcountPattern("GCG",genome)
#    GCGtotal=sum(GCGcounts[GCGcounts$strand=="+","count"][1:chrNum])
#    CGCcounts=vcountPattern("CGC",genome)
#    CGCtotal=sum(CGCcounts[CGCcounts$strand=="+","count"][1:chrNum])
#    
#    lengthTotal=sum(seqlengths(genome)[1:chrNum])
#    
#    print(paste("There is 1 CpG every",round(lengthTotal/CGtotal,0),"bp in the genome"))
#    print(paste("There is 1 GpC every",round(lengthTotal/GCtotal,0),"bp in the genome"))
#    print(paste("A CpG or GpC occur every",round(lengthTotal/(CGtotal+GCtotal-GCGtotal-CGCtotal),0), "bp in the genome"))
# }


library(BSgenome.Dmelanogaster.UCSC.dm6)

genome <- BSgenome.Dmelanogaster.UCSC.dm6
gpc<-vmatchPattern("GC",genome)
cpg<-vmatchPattern("CG",genome)
ol<-findOverlaps(gpc,cpg)
cggc<-c(gpc[-queryHits(ol)],cpg)

getMedianGap(gpc,genome)
#10
getMedianGap(cpg,genome)
#13
getMedianGap(cggc,genome)
#7

#countCGGC(genome,7)
#"There is 1 CpG every 24 bp in the genome"
#"There is 1 GpC every 18 bp in the genome"
#"A CpG or Gpc occur every 12 bp in the genome"

#countCGGCds(genome,7)
#"There is 1 CpG every 12 bp in the genome"
#"There is 1 GpC every 9 bp in the genome"
#"A CpG or Gpc occur every 6 bp in the genome"