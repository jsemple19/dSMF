###########  useful functions

ucsc2wb<-function(BCobject) {
   seqlevels(BCobject)<-gsub("chr","",seqlevels(BCobject))
   seqlevels(BCobject)<-gsub("M","MtDNA",seqlevels(BCobject))
   return(BCobject)
}

wb2ucsc<-function(BCobject) {
   seqlevels(BCobject)<-gsub("MtDNA","M",seqlevels(BCobject))
   seqlevels(BCobject)<-paste0("chr",seqlevels(BCobject))
   return(BCobject)
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


### mode function
Mode <- function(x) {
   ux <- unique(x)
   ux[which.max(tabulate(match(x,ux)))]
}



######### example of tryCatch when some subscripts are NAs
getStartTSS<-function(modes,dataName) {
   result<-tryCatch( 
      {
         start(ol[modes[dataName]])
      }, error=function(cond) {
         return(NA)
      })
   return(result)
}
