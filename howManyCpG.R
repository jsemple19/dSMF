library(BSgenome.Celegans.UCSC.ce11)
library(Biostrings)
genome <- BSgenome.Celegans.UCSC.ce11

CpG="CG"
GpC="GC"

# #seqnames(genome)
# CGtotal=0
# GCtotal=0
# GCGtotal=0
# lengthTotal=0
# 
# for (chr in seqnames(genome)[1:6]){
#    CGtotal=CGtotal+countPattern(CpG,genome[[chr]])
#    GCtotal=GCtotal+countPattern(GpC,genome[[chr]])
#    GCGtotal=GCGtotal+countPattern("GCG",genome[[chr]])
#    lengthTotal=lengthTotal+length(genome[[chr]])
# }


CGcounts=vcountPattern(CpG,genome)
CGtotal=sum(CGcounts[CGcounts$strand=="+","count"][1:6])
GCcounts=vcountPattern(GpC,genome)
GCtotal=sum(GCcounts[GCcounts$strand=="+","count"][1:6])

#remove duplicate counts of Cs between to Gs
GCGcounts=vcountPattern("GCG",genome)
GCGtotal=sum(GCGcounts[GCGcounts$strand=="+","count"][1:6])

lengthTotal=sum(seqlengths(genome)[1:6])


print(paste("There is 1 CpG every",round(lengthTotal/CGtotal,0),"bp in the genome"))
print(paste("There is 1 GpC every",round(lengthTotal/GCtotal,0),"bp in the genome"))
print(paste("A CpG or Gpc occur every",round(lengthTotal/(CGtotal+GCtotal-GCGtotal),0), "bp in the genome"))
