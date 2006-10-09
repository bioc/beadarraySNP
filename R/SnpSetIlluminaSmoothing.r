# Miscellaneous support
smoothed.intensity<-function(snpdata,smooth.lambda=4,tau=0.5) {
  chroms<-unique(reporterInfo(snpdata)[,"CHR"])
  ind<-order(numericCHR(reporterInfo(snpdata)[,"CHR"]),reporterInfo(snpdata)[,"MapInfo"])
  snpdata<-snpdata[ind,]
  il.smoothed<-matrix(NA,nrow=nrow(snpdata),ncol=ncol(snpdata),dimnames=list(featureNames(snpdata),sampleNames(snpdata)))
  for (i1 in 1:ncol(snpdata) ) {
    for (chrom in chroms) {
      probes<-reporterInfo(snpdata)[,"CHR"] == chrom
      if (sum(!is.na(assayData(snpdata)$intensity[probes,i1]))>10){
        il.smoothed[probes,i1]<-quantsmooth(assayData(snpdata)$intensity[probes,i1],tau=tau,smooth.lambda=smooth.lambda,smooth.na=FALSE,segment=100)
      }
    }
  }
  il.smoothed
}

renameOPA<-function(snpdata,newOPA) {
  annotation(snpdata)<-newOPA
  pData(featureData(snpdata))$OPA<-rep(newOPA,length(featureNames(snpdata)))
  snpdata
}

