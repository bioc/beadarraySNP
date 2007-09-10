# segmentation on SNPSetIllimuna objects
# author: J. Oosting
segmentate<-function(object, method=c("DNACopy","HMM","BioHMM","GLAD"), normalizedTo=2, doLog=TRUE, doMerge= FALSE, subsample="OPA") {
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
