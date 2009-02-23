# Miscellaneous support
smoothed.intensity<-function(snpdata,smooth.lambda=4,tau=0.5) {
  chroms<-unique(fData(snpdata)[,"CHR"])
  snpdata<-sortGenomic(snpdata)
  il.smoothed<-matrix(NA,nrow=nrow(snpdata),ncol=ncol(snpdata),dimnames=list(featureNames(snpdata),sampleNames(snpdata)))
  for (i1 in 1:ncol(snpdata) ) {
    for (chrom in chroms) {
      probes<-fData(snpdata)[,"CHR"] == chrom
      if (sum(!is.na(assayData(snpdata)$intensity[probes,i1]))>10){
        il.smoothed[probes,i1]<-quantsmooth(assayData(snpdata)$intensity[probes,i1],tau=tau,smooth.lambda=smooth.lambda,smooth.na=FALSE,segment=100)
      }
    }
  }
  il.smoothed
}

smoothed.probes<-function(snpdata,sProbes=10) {
  
}

smoothed.size<-function(snpdata,sSize=1e6) {

}

setMethod("calculateSmooth", c("SnpSetIllumina","character"), function(object,smoothType,...) {
  smoothType<-match.arg(smoothType,c("quantsmooth","probes","size"))
  smoothed<-switch(smoothType,
    quantsmooth=smoothed.intensity(object,...)
  )
  res<-object
  assayData(res)$smoothed<-smoothed
  res
})

renameOPA<-function(snpdata,newOPA) {
  annotation(snpdata)<-newOPA
  fData(snpdata)$OPA<-rep(newOPA,length(featureNames(snpdata)))
  snpdata
}



