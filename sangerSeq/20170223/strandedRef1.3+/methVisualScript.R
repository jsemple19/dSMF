library(methVisual)
setwd("~/Documents/MeisterLab/dSMF/sangerSeq/20170223/strandedRef1.3+")

#need to trim 5' to facilitate alignment. (and revcomp if necessary)
methData <-MethDataInput("./fileList.txt")

#refseq <- selectRefSeq("./refSeq.txt")
refseq<-selectRefSeq(file.path("/home/jenny/Documents/MeisterLab/dSMF/sangerSeq/20170223/strandedRef1.3+","WBGene00016114.fasta"))

QCdata <- MethylQC(refseq, methData,makeChange=TRUE,identity=80,conversion=90)
methData <- MethAlignNW( refseq , QCdata)
methData

pdf("sangerSeqMethResults.pdf",paper="a4",width=8,height=11)
par(mfrow=c(2,1))
plotAbsMethyl(methData,real=TRUE)

a<-MethLollipops(methData)
dev.off()

TSS<-217
flip=TRUE
#modifiction of MethLollipops function from methVisual package. requires package to be loaded.
methLollipopsReal<-function (methData,TSS,flip=FALSE) {
   lRef <- methData$lengthRef
   cgPos <- methData$positionCGIRef
   if (flip==TRUE) {
      cgPos<-lRef-rev(cgPos)
      TSS<-lRef-TSS
   }
   ncgPos <- length(cgPos)
   cgPos_real <- cgPos
   cgPos <- c(1:ncgPos)
   plot(c(0, lRef*1.05), c(0, 10), yaxt = "n", 
        type = "n", xlab = "position in amplicon (bp)",ylab=NA,bty="n")
   polygon(c(0,lRef,lRef,0),c(0,0,0.5,0.5))
   for (i in 1:ncgPos) {
      lines(rep(cgPos_real[i],2),c(0,0.5),col="black",lwd=1)
   }
   lines(rep(TSS,2),c(0.5,1),col="red",lwd=2)
   text(TSS,1.5,labels="TSS",cex=0.8,col="red")
}


methLollipopsScale<-function (methData,TSS,flip=FALSE) {
   lRef <- methData$lengthRef
   nrow <- lRef
   ncol <- lRef
   methMatrix <- findNonAligned(methData)
   cgPos <- methData$positionCGIRef
   if (flip==TRUE) {
      cgPos<-lRef-rev(cgPos)
      TSS<-lRef-TSS
   }
   TSS<-max(which(TSS>cgPos))+0.5
   ncgPos <- length(cgPos)
   lengthData <- length(methData$seqName)
   methyl <- c()
   summeryMatrix_temp <- matrix(0, nrow = ncgPos, ncol = ncgPos)
   summeryMatrix <- matrix(0, nrow = ncgPos, ncol = ncgPos)
   #cgPos_real <- cgPos
   cgPos <- c(1:ncgPos)
   plot(c(0, ncgPos+1), c(0, lengthData + 1), yaxt = "n", 
        type = "n", xlab = "index of CpG methylation", ylab = "index of clone sequences")
   grid((ncgPos-1), NA, lwd = 1)
   methInAllSeq <- seq(0, 0, length.out = lRef)
   vec1 <- seq(1, 1, length.out = ncgPos)
   methInAllSeq[cgPos] <- vec1
   y_methyl <- seq(lengthData + 1, lengthData + 1, length.out = ncgPos)
   #spacedPos<-lRef*(cgPos)/ncgPos
   lines(cgPos, y_methyl, type = "p", col = "black", lwd = "2")
   abline(v=TSS,col="red")
   methInAllSeq <- seq(3, 3, length.out = lRef)
   for (i in 1:lengthData) {
      seq_temp <- as.array(methMatrix[i, ])
      methInAllSeq_temp <- methInAllSeq
      methInAllSeq_temp[cgPos] <- seq_temp
      non_align <- which(methInAllSeq_temp == 2)
      methyl <- which(methInAllSeq_temp == 1)
      non_methyl <- which(methInAllSeq_temp == 0)
      if (flip==TRUE) {
         non_align <- ncgPos+1-rev(which(methInAllSeq_temp == 2))
         methyl <- ncgPos+1-rev(which(methInAllSeq_temp == 1))
         non_methyl <- ncgPos+1-rev(which(methInAllSeq_temp == 0))
      }
      y_non_align <- seq(i, i, length.out = length(non_align))
      y_methyl <- seq(i, i, length.out = length(methyl))
      y_non_methyl <- seq(i, i, length.out = length(non_methyl))
      axis(side = 2, at = c(1:(lengthData + 1)), labels = c(1:lengthData, 
                                                            "refSeq"), tick = TRUE, par(las = 2))
      lines(non_align, y_non_align, yaxt = "n", type = "p", 
            col = NULL, lwd = "2")
      lines(methyl, y_methyl, type = "p", pch = 19, col = "black", 
            lwd = "5")
      lines(non_methyl, y_non_methyl, type = "p", col = "black", 
            lwd = "2")
   }
   LABEL_Y_AXIS <- c(1:(lengthData), "refSeq")
   Experiment <- c(methData$seqName, "referenceSequence")
   numberExp <- data.frame(LABEL_Y_AXIS = LABEL_Y_AXIS, Experiment = Experiment)
   return(numberExp)
   gc(FALSE)
}



##########################################

#modifiction of MethLollipops function from methVisual package. requires package to be loaded.

methLollipopsAll<-function (methData,TSS,flip=FALSE) {
   lRef <- methData$lengthRef
   nrow <- lRef
   ncol <- lRef
   methMatrix <- findNonAligned(methData)
   cgPos <- methData$positionCGIRef
   if (flip==TRUE) {
      cgPos<-lRef-rev(cgPos)
      TSS_real<-lRef-TSS
   }
   TSS<-max(which(TSS_real>cgPos))+0.5
   ncgPos <- length(cgPos)
   lengthData <- length(methData$seqName)
   methyl <- c()
   summeryMatrix_temp <- matrix(0, nrow = ncgPos, ncol = ncgPos)
   summeryMatrix <- matrix(0, nrow = ncgPos, ncol = ncgPos)
   cgPos_real <- cgPos
   cgPos <- c(1:ncgPos)
   par(mar=c(5,5,5,2))
   plot(c(0, ncgPos+1), c(0, lengthData + 5), bty="n", ylim=c(0,lengthData+5),yaxt="n",
        type = "n", xlab = "index of CpG methylation", ylab = "index of clone sequences",
        las=1,cex.axis=0.8,xaxs="i")
   #grid((ncgPos-1), NA, lwd = 1)
   methInAllSeq <- seq(0, 0, length.out = lRef)
   vec1 <- seq(1, 1, length.out = ncgPos)
   methInAllSeq[cgPos] <- vec1
   y_methyl <- seq(lengthData + 1, lengthData + 1, length.out = ncgPos)
   lines(cgPos, y_methyl, type = "p", col = "black", lwd = "2")
   lines(c(TSS,TSS),c(0.6,lengthData+1.3),col="red")
   text(TSS,0.3,labels="TSS",cex=0.8,col="red")
   methInAllSeq <- seq(3, 3, length.out = lRef)
   for (i in 1:lengthData) {
      seq_temp <- as.array(methMatrix[i, ])
      methInAllSeq_temp <- methInAllSeq
      methInAllSeq_temp[cgPos] <- seq_temp
      non_align <- which(methInAllSeq_temp == 2)
      methyl <- which(methInAllSeq_temp == 1)
      non_methyl <- which(methInAllSeq_temp == 0)
      if (flip==TRUE) {
         non_align <- ncgPos+1-rev(which(methInAllSeq_temp == 2))
         methyl <- ncgPos+1-rev(which(methInAllSeq_temp == 1))
         non_methyl <- ncgPos+1-rev(which(methInAllSeq_temp == 0))
      }
      y_non_align <- seq(i, i, length.out = length(non_align))
      y_methyl <- seq(i, i, length.out = length(methyl))
      y_non_methyl <- seq(i, i, length.out = length(non_methyl))
      axis(side = 2, at = c(1:(lengthData + 1)), 
           labels = c(1:lengthData,"refSeq"), tick = TRUE, par(las = 2),cex.axis=0.8)
      lines(non_align, y_non_align, yaxt = "n", type = "p", 
            col = NULL, lwd = "2")
      lines(methyl, y_methyl, type = "p", pch = 19, col = "black", 
            lwd = "5")
      lines(non_methyl, y_non_methyl, type = "p", col = "black", 
            lwd = "2")
   }
   ##
   par(new=T)
   plot(c(0, lRef*1.05), c(0, lengthData+5), axes=FALSE, type = "n",xlab=NA ,ylab=NA,bty="n",
        ylim=c(0,lengthData+5),xaxs="i")
   #xlab = "position in amplicon (bp)"
   polygon(c(0,lRef,lRef,0),c(lengthData+3,lengthData+3,lengthData+3.5,lengthData+3.5))
   spacedPos<-0.97*lRef*(cgPos)/(ncgPos)
   for (i in 1:ncgPos) {
      lines(rep(cgPos_real[i],2),c(lengthData+3,lengthData+3.5),col="black",lwd=1)
      lines(c(cgPos_real[i],spacedPos[i]),c(lengthData+3,lengthData+1.2),lty=3,lwd=1,
            col="grey")
   }
   lines(rep(TSS_real,2),c(lengthData+3.5,lengthData+3.9),col="red",lwd=2)
   text(TSS_real,lengthData+4.2,labels="TSS",cex=0.8,col="red")
   axis(side=3,cex.axis=0.8,las=1)
   mtext(side=3,line=3,"position in amplicon (bp)",las=1)
   ##
   LABEL_Y_AXIS <- c(1:(lengthData), "refSeq")
   Experiment <- c(methData$seqName, "referenceSequence")
   numberExp <- data.frame(LABEL_Y_AXIS = LABEL_Y_AXIS, Experiment = Experiment)
   return(numberExp)
   gc(FALSE)
}

##### modify to plot from matrix from MiSeq data:

methLollipopsAll<-function (methData,TSS,flip=FALSE) {
   lRef <- methData$lengthRef
   nrow <- lRef
   ncol <- lRef
   methMatrix <- findNonAligned(methData)
   cgPos <- methData$positionCGIRef
   if (flip==TRUE) {
      cgPos<-lRef-rev(cgPos)
      TSS_real<-lRef-TSS
   }
   TSS<-max(which(TSS_real>cgPos))+0.5
   ncgPos <- length(cgPos)
   lengthData <- length(methData$seqName)
   methyl <- c()
   summeryMatrix_temp <- matrix(0, nrow = ncgPos, ncol = ncgPos)
   summeryMatrix <- matrix(0, nrow = ncgPos, ncol = ncgPos)
   cgPos_real <- cgPos
   cgPos <- c(1:ncgPos)
   par(mar=c(5,5,5,2))
   plot(c(0, ncgPos+1), c(0, lengthData + 5), bty="n", ylim=c(0,lengthData+5),yaxt="n",
        type = "n", xlab = "index of CpG methylation", ylab = "index of clone sequences",
        las=1,cex.axis=0.8,xaxs="i")
   #grid((ncgPos-1), NA, lwd = 1)
   methInAllSeq <- seq(0, 0, length.out = lRef)
   vec1 <- seq(1, 1, length.out = ncgPos)
   methInAllSeq[cgPos] <- vec1
   y_methyl <- seq(lengthData + 1, lengthData + 1, length.out = ncgPos)
   lines(cgPos, y_methyl, type = "p", col = "black", lwd = "2")
   lines(c(TSS,TSS),c(0.6,lengthData+1.3),col="red")
   text(TSS,0.3,labels="TSS",cex=0.8,col="red")
   methInAllSeq <- seq(3, 3, length.out = lRef)
   for (i in 1:lengthData) {
      seq_temp <- as.array(methMatrix[i, ])
      methInAllSeq_temp <- methInAllSeq
      methInAllSeq_temp[cgPos] <- seq_temp
      non_align <- which(methInAllSeq_temp == 2)
      methyl <- which(methInAllSeq_temp == 1)
      non_methyl <- which(methInAllSeq_temp == 0)
      if (flip==TRUE) {
         non_align <- ncgPos+1-rev(which(methInAllSeq_temp == 2))
         methyl <- ncgPos+1-rev(which(methInAllSeq_temp == 1))
         non_methyl <- ncgPos+1-rev(which(methInAllSeq_temp == 0))
      }
      y_non_align <- seq(i, i, length.out = length(non_align))
      y_methyl <- seq(i, i, length.out = length(methyl))
      y_non_methyl <- seq(i, i, length.out = length(non_methyl))
      axis(side = 2, at = c(1:(lengthData + 1)), 
           labels = c(1:lengthData,"refSeq"), tick = TRUE, par(las = 2),cex.axis=0.8)
      lines(non_align, y_non_align, yaxt = "n", type = "p", 
            col = NULL, lwd = "2")
      lines(methyl, y_methyl, type = "p", pch = 19, col = "black", 
            lwd = "5")
      lines(non_methyl, y_non_methyl, type = "p", col = "black", 
            lwd = "2")
   }
   ##
   par(new=T)
   plot(c(0, lRef*1.05), c(0, lengthData+5), axes=FALSE, type = "n",xlab=NA ,ylab=NA,bty="n",
        ylim=c(0,lengthData+5),xaxs="i")
   #xlab = "position in amplicon (bp)"
   polygon(c(0,lRef,lRef,0),c(lengthData+3,lengthData+3,lengthData+3.5,lengthData+3.5))
   spacedPos<-0.97*lRef*(cgPos)/(ncgPos)
   for (i in 1:ncgPos) {
      lines(rep(cgPos_real[i],2),c(lengthData+3,lengthData+3.5),col="black",lwd=1)
      lines(c(cgPos_real[i],spacedPos[i]),c(lengthData+3,lengthData+1.2),lty=3,lwd=1,
            col="grey")
   }
   lines(rep(TSS_real,2),c(lengthData+3.5,lengthData+3.9),col="red",lwd=2)
   text(TSS_real,lengthData+4.2,labels="TSS",cex=0.8,col="red")
   axis(side=3,cex.axis=0.8,las=1)
   mtext(side=3,line=3,"position in amplicon (bp)",las=1)
   ##
   LABEL_Y_AXIS <- c(1:(lengthData), "refSeq")
   Experiment <- c(methData$seqName, "referenceSequence")
   numberExp <- data.frame(LABEL_Y_AXIS = LABEL_Y_AXIS, Experiment = Experiment)
   return(numberExp)
   gc(FALSE)
}



pdf("sangerSeqMethResults_flipped.pdf",paper="a4",width=8,height=11)
par(mfrow=c(2,1))
TSS<-217
flip=TRUE

#plotAbsMethyl(methData,real=TRUE)
methLollipopsReal(methData,TSS,flip)
methLollipopsScale(methData,TSS,flip)
methLollipopsAll(methData,TSS,flip=TRUE)

dev.off()
