# 2016-12-14
# Adc_design_primers_template.R
# modification of Arnaud Kreb's script. Used to design primers for autosomal genes from
# two subgroups (dosage compensated and not dosage compensated)

setwd("/SharedDocuments/MeisterLab/dSMF/PromoterPrimerDesign/scripts/")
library(QuasR)

library("BSgenome.Celegans.UCSC.ce11")

source('/SharedDocuments/MeisterLab/dSMF/PromoterPrimerDesign/scripts/useful_functionsV1.r') #load the ranges
source('/SharedDocuments/MeisterLab/dSMF/PromoterPrimerDesign/scripts/Functions_batch_bis_primer_designV201014_dev.r')

genomeVer="WS250"

#added by Jenny
outputDir="./AdcPrimers"
if (!dir.exists(outputDir)) {
  dir.create(outputDir)
}

  ##########
  chosenTSS<-readRDS("chosenTSS_A.Rds")
  TSS<-readRDS("../../TSS/scripts/confirmedTSS_Chen&Kreus_WS235.Rds")
  names(TSS)<-mcols(TSS)$WBGeneID
  
  chosenGR<-lapply(chosenTSS,function(x){ 
    TSS[x]
  })
  
  #########################
	#general priming parametets
	###########################
		size_range='151-500'
		Tm=c(50,52,53)
		
		Trange=c(100,100)

		maxD=2000 #max number of regions considered
	#########################
	#Design primers
	###########################
		
		
		maxS=as.numeric(string.split(size_range,'-',2))	+200
		
		grts=chosenGR#list(target_regions) #,target_regions2) #list of regions
	  #grts=targets
	  #seqlevels(targets)<-sort(seqnames(Celegans))
	    #c("chrI","chrII","chrIII","chrIV","chrM","chrV","chrX")
	  
		design=lapply(sl(grts),function(j){
			#j=1
			#print(paste('j=',j,sep=''))
			igr=grts[[j]]
	
			if(length(igr)>maxD){igr=igr[sample(sl(igr),maxD)]} #randomly select N elements 

	
      target=igr#ectg.tss[ectg.tss.i]
 
			designgr=GRanges(seqnames(target),IRanges(start(target)-Trange[1],start(target)+Trange[2]))
		  names(designgr)<-names(target)
	
			sample.seq=getSeq(Celegans,resize(target,maxS,fix='center'))
			names(sample.seq)=grToNames(target)
			#get terget regions in relative space
			target_regions=paste(round((maxS/2)- Trange[1]),width(designgr),sep=',') #get fw strand
			target_regions.rc=paste(round((maxS/2)- Trange[2]),width(designgr),sep=',') #get reverse strand
	
			#design primers
	
			primers.fw=lapply(seq(length(sample.seq)),function(i){
			  design_bis_primers_GC_dev(sample.seq[i], Tm=Tm, size_range=size_range,target_region=target_regions[i],
			                            primer3_settings_file="/SharedDocuments/MeisterLab/dSMF/PromoterPrimerDesign/scripts/default_settings_10.txt")})
			names(primers.fw)=names(sample.seq)
	
			primers.rc=lapply(seq(length(sample.seq)),function(i){
			  design_bis_primers_GC_dev(reverseComplement(sample.seq[i]),Tm=Tm, size_range=size_range,target_region=target_regions.rc[i],
			                            primer3_settings_file="/SharedDocuments/MeisterLab/dSMF/PromoterPrimerDesign/scripts/default_settings_10.txt")})
			names(primers.rc)=names(sample.seq)
	
	
			#retrive longest amplicon (incl from both strands)
			lprimers=lapply(sl(sample.seq),function(i)
				if(length(primers.fw[[i]])==14 | length(primers.rc[[i]])==14){	
					fw=rep(NA,14)
					rc=rep(NA,14)	
					if(length(primers.fw[[i]])==14 ){
					fw=primers.fw[[i]]
					fw[,1]='fw'
						}
					if(length(primers.rc[[i]])==14 ){
					rc=primers.rc[[i]]
					rc[,1]='rc'}
					
					primers=(rbind(fw,rc))
					primers[order(primers$fragLen,decreasing=T),][1,]
	
				}else(rep(NA,14))
				)
			#data frame
			p.df=do.call(rbind,lprimers)
		
			#reference GR
			

			dfr=p.df#prim_df_list[[i]]
			IDs=dfr$PrimerID 
			snps=resize(designgr,1,fix='center')#=#ref_gr_list[[i]][IDs]
			stP=start(snps)-(maxS/2)
			enP=start(snps)+(maxS/2)
			nai=!is.na(dfr$Fwseq)
			strand=ifelse(p.df$PrimerID=='fw','+','-' )
			refGR=GRanges(as.character(seqnames(snps))[nai],
				      IRanges(ifelse(strand[nai]=='+',stP[nai]+dfr$FwPos[nai]-1,
				                enP[nai]-dfr$FwPos[nai]-dfr$fragLen[nai]-1),
				              ifelse(strand[nai]=='+',stP[nai]+dfr$FwPos[nai]+dfr$fragLen[nai]-1 ,
				                enP[nai]-dfr$FwPos[nai]-1)),
				      primer.class=igr$name[nai],strand=strand[nai])
			names(refGR)<-names(snps[nai]) #added by Jenny to retain WBgeneIDs in output
			list(p.df[nai,],refGR)
		})

	names(design)=names(grts)	
	
	lapply(sl(design), function(i){length(design[[i]][[2]])} )
	GRangesList(lapply(sl(design), function(i){(design[[i]][[2]])} ))[1]
	reg=unlist(GRangesList(lapply(design, '[[',2)))
	#names(reg)=sl(reg) # commented out by Jenny
	cbind(as.data.frame(reg)[,1:3],
	      "name"=row.names(as.data.frame(reg)), # added by Jenny
	      #rep('.',length(reg)), # commented out by Jenny
	      "score"=rep('.',length(reg)),
	      "strand"=as.data.frame(reg)[,5])

	#added by Jenny
	design<-lapply(design,function(x) {
	  x[[1]]<-cbind("Gene_WB_ID"=names(x[[2]]),x[[1]])
	  return(x)
	})	
	
	#export the list containing both  primers and reference GR
	saveRDS(design ,paste0(outputDir,'/tmp_design_test.rds'))
	primer.df=do.call(rbind,lapply(design, '[[',1))
	#remove duplicates
	
	#export the primer data frame
	write.table(primer.df[!(duplicated(primer.df$Fwseq)&duplicated(primer.df$Rvseq)),],
	            paste0(outputDir,'/DM_DE_promoter_test_primers_',genomeVer,'.tab'),sep='\t',quote=F)

	#export the bed file for UCSC visualization
	write.table(cbind(as.data.frame(reg)[,1:3],
	                  row.names(as.data.frame(reg)), #added by Jenny
	                  #rep('',length(reg)), # commented out by Jenny
	                  rep('',length(reg)),
	                  as.data.frame(reg)[,5]),
	            paste0(outputDir,'/test_primers_bed_',genomeVer,'.bed'),sep='\t',quote=F,row.names=F,col.names=F)

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	