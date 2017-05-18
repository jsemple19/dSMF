#functions to go with findMotifs2.R file
############################# functions ################################3
findBestMotif<-function(pwm,amplicon,relTSS,expectedPos,rc=TRUE,pattern=FALSE) {
   if (pattern==FALSE) {
      matches<-matchPWM(pwm,amplicon)
   } else {
      matches<-matchPattern(pwm,amplicon,max.mismatch=1,fixed=FALSE)
   }
   if (rc==TRUE) {
      if (pattern==FALSE) {
         matchesRC<-matchPWM(reverseComplement(pwm),amplicon)
      } else {
         matchesRC<-matchPattern(reverseComplement(DNAString(pwm)),amplicon,max.mismatch=1,fixed=FALSE)
      }
      matches<-c(matches,matchesRC)
   }
   bestMatch<-which.min(abs(start(matches)-expectedPos-relTSS))
   return(matches[bestMatch])
}

findBestMotif_genomeCoord<-function(pwm,amplicon,ampMeta,relTSS,expectedPos,rc=TRUE,pattern=FALSE) {
   if (pattern==FALSE) {
      matches<-matchPWM(pwm,amplicon)
   } else {
      matches<-matchPattern(pwm,amplicon,max.mismatch=1,fixed=FALSE)
   }
   if (rc==TRUE) {
      if (pattern==FALSE) {
         matchesRC<-matchPWM(reverseComplement(pwm),amplicon)
      } else {
         matchesRC<-matchPattern(reverseComplement(DNAString(pwm)),amplicon,max.mismatch=1,fixed=FALSE)
      }
      matches<-c(matches,matchesRC)
   }
   bestMatchInd<-which.min(abs(start(matches)-expectedPos-relTSS))
   bestMatch<-matches[bestMatchInd]
   names(bestMatch)<-ampMeta$Gene_WB_ID
   #mcols(bestMatch)<-DataFrame("motifMatch"=as.character(matches),
   #                            "name"=names(as.character(matches)))
   bestMatchGR<-matches2GR(ampMeta,bestMatch)
   mcols(bestMatchGR)$relPos<-start(bestMatch)-relTSS
   return(bestMatchGR)
}


findAllMotifs_genomeCoord<-function(pwm,amplicon,ampMeta,relTSS,rc=TRUE,pattern=FALSE) {
   if (pattern==FALSE) {
      matches<-matchPWM(pwm,amplicon)
   } else {
      matches<-matchPattern(pwm,amplicon,max.mismatch=1,fixed=FALSE)
   }
   if (rc==TRUE) {
      if (pattern==FALSE) {
         matchesRC<-matchPWM(reverseComplement(pwm),amplicon)
      } else {
         matchesRC<-matchPattern(reverseComplement(DNAString(pwm)),amplicon,max.mismatch=1,fixed=FALSE)
      }
      matches<-c(matches,matchesRC)
   }
   relPos<-start(matches)-relTSS
   motifName<-getOriginalVarName(pwm)
   names(matches)<-paste(ampMeta$Gene_WB_ID,motifName,seq(1,length(matches)),sep="_")
   motifMatch<-matches2GR(ampMeta,matches)
   mcols(motifMatch)$relPos<-start(matches)-relTSS
   return(motifMatch)
}

getOriginalVarName<-function(x) {
   my.call<-quote(substitute(x))
   var.name<-eval(my.call)
   for(i in rev(head(sys.frames(), -1L))) { # First frame doesn't matter since we already substituted for first level, reverse since sys.frames is in order of evaluation, and we want to go in reverse order
      my.call[[2]] <- var.name         # this is where we re-use it, modified to replace the variable
      var.name <- eval(my.call, i)
   }
   return(var.name)
}

matches2GR<-function(ampMeta,matches) {
   matchStart<-ifelse(rep(as.character(ampMeta$strand)=="+",length(matches)),
                      (ampMeta$start+start(matches)-1),
                      (ampMeta$end-end(matches)-1))
   matchWidth<-width(matches)
   matchGR<-GRanges(seqnames=rep(ampMeta$chr,length(matches)),
                    IRanges(start=matchStart,width=matchWidth),
                    strand=rep(ampMeta$strand,length(matches)))
   mcols(matchGR)<-DataFrame("motifMatch"=as.character(matches),
                             "name"=names(as.character(matches)))
   return(matchGR)
}

GR2IGV<-function(GR) {
   igv_coord<-paste0(seqnames(GR),":",start(GR),"-",end(GR))
   return(igv_coord)
}

readMemeChipMotif<-function(motifFile) {
   counts<-as.matrix(read.table(motifFile))
   counts<-t(counts)
   rownames(counts)<-c("A","C","G","T")
   colnames(counts)=c(1:dim(counts)[2])
   probs<-counts/colSums(counts)
   return(probs)
}

findPositionalPattern<-function(TSS, expPos,width,pattern) {
   rangePos<-resize(TSS+expPos,width,fix="center",ignore.strand=FALSE)
   #library("BSgenome.Celegans.UCSC.ce11")
   #library(Biostrings)
   genome<-BSgenome.Celegans.UCSC.ce11
   mySeqs<-getSeq(genome,rangePos)
   posPattern<-vmatchPattern(pattern,mySeqs)
   return(posPattern)
}

getRelTSS<-function(strand,tss,ampStart,ampEnd) {
   if (strand=="+"){
      relDistance<-tss-ampStart
   } else {
      relDistance<-ampEnd-tss
   }
   return(relDistance)
}

getPosRel2TSS<-function(strand,tss,pos) {
   if (strand=="+"){
      relDistance<-pos-tss
   } else {
      relDistance<-tss-pos
   }
   return(relDistance)
}

searchAmplicons4pwm<-function(amplicons,ampMeta,pwm,relTSS,expectedPos,rc=TRUE,pattern=FALSE) {
   allBestMatches<-cbind(rep(NA,length(amplicons)),rep(NA,length(amplicons)))
   for (i in 1:length(amplicons)) {
      motifMatches<-findBestMotif(pwm,amplicons[[i]],relTSS[i],expectedPos,rc,pattern)  
      if (length(motifMatches)==0){
         next
      } else {
         ampMeta<-mcols(amplicons)[i,]
         best<-findBestMotif_genomeCoord(pwm, amplicons[[i]], ampMeta,relTSS[i],
                                         expectedPos,rc,pattern)
         allBestMatches[i,1]<-GR2IGV(best)
         allBestMatches[i,2]<-best$relPos
         newGR<-findAllMotifs_genomeCoord(pwm,amplicons[[i]],ampMeta,relTSS[i],rc,pattern)
         if(exists("allMatches")) {
            allMatches<-c(allMatches,newGR)
         } else {
            allMatches<-newGR
         }
      }
   }
   motifType<-paste0(deparse(substitute(pwm)),"_",deparse(substitute(relTSS)))
   export.bed(allMatches, paste("all",motifType,"amplicons.bed",sep="_"))
   allBestMatches<-as.data.frame(allBestMatches,stringsAsFactors=FALSE)
   names(allBestMatches)<-c(motifType,"relPos")
   allBestMatches$relPos<-as.numeric(allBestMatches$relPos)
   return(allBestMatches)
}


#################################################################

