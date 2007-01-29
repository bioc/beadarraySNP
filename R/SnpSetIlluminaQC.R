# QC
removeLowQualityProbes<-function(object, cutoff=0.25) {
  # disable probes that have low overall intensity, ie fraction below overal median
  if ( !("intensity" %in% assayDataElementNames(object))) object<-RG2polar(object)
  probes<-numericCHR(reporterInfo(object)$CHR)>90
  int.med<-median(assayData(object)[["intensity"]][probes,],na.rm=TRUE)
  int.probe.med<-apply(assayData(object)[["intensity"]],1,median,na.rm=TRUE)
  probes<-int.probe.med>(int.med*cutoff) | probes
  cat(sum(!probes),"probes removed with median intensity below",int.med*cutoff,"\n")
  object[probes,]
}

removeLowQualitySamples<-function(object, min.intensity=1500, min.gt=100, subsample="OPA"){
  # remove samples that show median intensity below min.intensity in either red or green channel
  # remove samples that have less than min.gt genotyped probes
  subsample<-getSubsample(object,subsample)

  R<-assayDataElement(object,"R")
  G<-assayDataElement(object,"G")
  my_call<-assayDataElement(object,"call")
  callProbability<-assayDataElement(object,"callProbability")
  for (OPA in levels(subsample)) {
    probes<-subsample == OPA # select all probes of this OPAset
    greenmed<-apply(G[probes,],2,median,na.rm=TRUE)
    redmed<-apply(R[probes,],2,median,na.rm=TRUE)
    untyped<-apply(my_call[probes,],2,function(x) sum(x != "U",na.rm=TRUE))
    low<- greenmed<min.intensity | redmed<min.intensity | untyped<min.gt
    cat(OPA,sum(!is.na(low)),"subsamples, ")
    low[is.na(low)]<-FALSE
    cat(sum(low),"removed\n")
    G[probes,low]<-NA
    R[probes,low]<-NA
    my_call[probes,low]<-NA
    callProbability[probes,low]<-NA
  }
  object<-assayDataElementReplace(object,"R",R)
  object<-assayDataElementReplace(object,"G",G)
  object<-assayDataElementReplace(object,"call",my_call)
  object<-assayDataElementReplace(object,"callProbability",callProbability)
  allmissing<-apply(G,2,function(x) all(is.na(x)))
  if (any(allmissing)) {
    warning(paste("Sample(s)",toString(sampleNames(object)[allmissing]),"removed from result because no data remained"))
    object<-object[,!allmissing]
  }
  object
}

calculateQCarray<-function(object,QCobject=NULL) {
  # object should be SnpSetIllumina
  if (class(object)!="SnpSetIllumina") stop("object not usable for this function")
  # object should not be combined 
  if (length(annotation(object))>1) stop("function does not work with combined datasets")
  if (is.null(QCobject)) {
     QCobject<-new("QCIllumina")
     arrayID(QCobject)<-pData(object)$Sentrix_ID[1]
  }
  R<-assayDataElement(object,"R")
  G<-assayDataElement(object,"G")
  int<-R+G
  for (smp in sampleNames(object)) {
    co<-pData(object)[smp,"Col"]
    ro<-pData(object)[smp,"Row"]
    if (pData(object)[smp,"Sentrix_ID"]!=arrayID(QCobject)) stop("data in QCobject cannot be from different arrays")
		QCobject@samples[ro,co]<-smp
		QCobject@annotation[ro,co]<-annotation(object)
		QCobject@validn[ro,co]<-sum(!is.na(int[,smp]))
		QCobject@intensityMed[ro,co]<-median(int[,smp],na.rm=TRUE)
		QCobject@greenMed[ro,co]<-median(G[,smp],na.rm=TRUE)
		QCobject@redMed[ro,co]<-median(R[,smp],na.rm=TRUE)
	
  }
  QCobject
}

