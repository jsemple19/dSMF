library(BSgenome.Celegans.UCSC.ce11)
library(Biostrings)
genome <- BSgenome.Celegans.UCSC.ce11

CpG="CG"
GpC="GC"

#seqnames(genome)
CGtotal=0
GCtotal=0
lengthTotal=0

for (chr in seqnames(genome)[1:6]){
   CGtotal=CGtotal+countPattern(CpG,genome[[chr]])
   GCtotal=GCtotal+countPattern(GpC,genome[[chr]])
   lengthTotal=lengthTotal+length(genome[[chr]])
}

print(paste("There is 1 CpG every",round(lengthTotal/CGtotal,0),"bp in the genome"))
print(paste("There is 1 GpC every",round(lengthTotal/GCtotal,0),"bp in the genome"))

print(paste("A CpG or Gpc occur every",round(lengthTotal/(CGtotal+GCtotal),0), "bp in the genome"))
