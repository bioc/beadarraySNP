# segmentation on SNPSetIllimuna objects
# author: J. Oosting



segmentate.old<-function(object, method=c("DNACopy","HMM","BioHMM","GLAD"), normalizedTo=2, doLog=TRUE, doMerge= FALSE, subsample="OPA") {
  # perform segmentation,
  # add states and predicted copy numbers to assayData
  method<-match.arg(method)
  if (!require(snapCGH))
    stop("package snapCGH is not installed")
  sl<-convert2SegList(object,normalizedTo=normalizedTo,doLog=doLog)
  subsample<-getSubsample(object,subsample)
  sl$genes<-cbind(sl$genes,subsample=subsample,
      pos.old=sl$genes[,"Position"],chr.old=sl$genes[,"Chr"])
  # create new mapping with 1 chromosome
  sl$genes[,"Position"]<-sl$genes[,"Position"]+(sl$genes[,"Chr"]*1e9)
  sl$genes[,"Chr"]<-rep(1,nrow(sl))
  subsamples<-unique(as.character(subsample))
  states<-matrix(NA,ncol=ncol(sl),nrow=nrow(sl),dimnames=dimnames(sl))
  predicted<-matrix(NA,ncol=ncol(sl),nrow=nrow(sl),dimnames=dimnames(sl))
  for (smp in subsamples) {
    selection<-subsample==smp
    switch(method, DNACopy = {
      segm<-runDNAcopy(sl[selection,])
    }, HMM = {
      segm<-runHomHMM(sl[selection,], criteria = "AIC")
    }, BioHMM = {
      segm<-runBioHMM(sl[selection,])
    }, GLAD = {
      segm<-runGLAD(sl[selection,])
    })
    if (doMerge)
      segm<-mergeStates(segm)
    states[selection,]<-segm$state
    predicted[selection,]<-segm$M.predicted
  }
  res<-object
  assayData(res)$observed<-sl$M.observed # in case of transformations intensity slot is not enough
  assayData(res)$states<-states
  assayData(res)$predicted<-predicted
  res
}

segmentate<-function(object, method=c("DNACopy","HMM","BioHMM","GLAD"), normalizedTo=2, doLog=TRUE, doMerge= FALSE, useLair=FALSE, subsample="OPA") {
  method<-match.arg(method)
  if (method=="new") {
    object<-sortGenomic(object)
    if (!require(DNAcopy))
      stop("package `DNAcopy` is not installed")
    observed<-assayData(object)$intensity/normalizedTo
    if (doLog) 
      observed<-log2(observed)  
    states<-matrix(NA,ncol=ncol(observed),nrow=nrow(observed),dimnames=dimnames(observed))
    predicted<-states
    cna.intensity<-smooth.CNA(CNA(observed,numericCHR(fData(object)$CHR),fData(object)$MapInfo,sampleid=sampleNames(object)))
    seg.intensity<-DNAcopy::segment(cna.intensity)
    samples<-unique(seg.intensity$output$ID)
    seg.intensity<-split(seg.intensity$output,seg.intensity$output$ID)
    res<-object
    assayData(res)$observed<-observed # in case of transformations intensity slot is not enough
    assayData(res)$states<-states
    assayData(res)$predicted<-predicted
    
    if (useLair) {
      lair.states<-matrix(NA,ncol=ncol(observed),nrow=nrow(observed),dimnames=dimnames(observed))
      lair.predicted<-lair.states
      lair<-assayData(object)$lair
      lair[!assayData(object)$nor.gt]<-NA
      cna.lair<-CNA(lair,numericCHR(fData(object)$CHR),fData(object)$MapInfo,sampleid=sampleNames(object))
      seg.lair<-DNAcopy::segment(cna.lair)
      assayData(res)$lair.states<-lair.states
      assayData(res)$lair.predicted<-lair.predicted
    }
    
    res
  } else segmentate.old(object, method, normalizedTo, doLog, doMerge, subsample)
}
