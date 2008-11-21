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

smooth.probes<-function(snpdata,nProbes) {

}

calculateSmooth<-function(snpdata,smoothType=c("quantsmooth","probes","size"),...) {
  smoothType<-match.arg(smoothType)
  smoothed<-switch(smoothType,
    quantsmooth=smoothed.intensity(snpdata,...)
  )
  res<-object
  assayData(res)$smoothed<-smoothed
  res
}

renameOPA<-function(snpdata,newOPA) {
  annotation(snpdata)<-newOPA
  fData(snpdata)$OPA<-rep(newOPA,length(featureNames(snpdata)))
  snpdata
}

sortGenomic<-function(snpdata) {
  snpdata[order(numericCHR(fData(snpdata)[,"CHR"]),fData(snpdata)[,"MapInfo"]),]
}


