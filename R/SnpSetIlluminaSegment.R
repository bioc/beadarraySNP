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
  #sl$genes[,"Position"]<-sl$genes[,"Position"]+(sl$genes[,"Chr"]*1e9)
  #sl$genes[,"Chr"]<-rep(1,nrow(sl))
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
  if (method=="DNACopy") {
    object<-sortGenomic(object)
    if (!require(DNAcopy))
      stop("package `DNAcopy` is not installed")
    observed<-assayData(object)$intensity/normalizedTo
    if (doLog) 
      observed<-log2(observed)  
    states<-matrix(NA,ncol=ncol(observed),nrow=nrow(observed),dimnames=dimnames(observed))
    predicted<-states
    if (useLair) {
      lair.states<-matrix(NA,ncol=ncol(observed),nrow=nrow(observed),dimnames=dimnames(observed))
      lair.predicted<-lair.states
      lair<-assayData(object)$lair
      lair[!assayData(object)$nor.gt]<-NA
    }
    chrom<-numericCHR(fData(object)$CHR)
    maploc<-fData(object)$MapInfo
    cna.intensity<-smooth.CNA(CNA(observed,chrom,maploc,sampleid=sampleNames(object)))
    seg.intensity<-DNAcopy::segment(cna.intensity)
    assays<-unique(seg.intensity$output$ID)
    seg.intensity<-split(seg.intensity$output,seg.intensity$output$ID)
    for (assay in 1:length(assays)) {
      seg.smp<-seg.intensity[[assays[assay]]]
      for(state in 1:nrow(seg.smp)) {
        states[chrom==seg.smp$chrom[state] & maploc>=seg.smp$loc.start[state] & maploc<=seg.smp$loc.end[state],assay]<-state
        predicted[chrom==seg.smp$chrom[state] & maploc>=seg.smp$loc.start[state] & maploc<=seg.smp$loc.end[state],assay]<-seg.smp$seg.mean[state]
      }
      if (useLair) {
        # Exclude segments that have no valid values
        all.na<-aggregate(lair[,assay],by=list(states[,assay]),FUN=function(x) all(is.na(x)))
        selection<-! states %in% all.na[all.na[,2],1]
        # 
        cna.lair<-CNA(lair[selection,assay],states[selection,assay],maploc[selection],sampleid=paste(sampleNames(object)[assay],"lair"))
        seg.lair<-DNAcopy::segment(cna.lair)
        seg.lair<-seg.lair$output
        if (any(all.na[,2])) { # add the excluded segments again
          exc.seg<-seg.smp[all.na[,2],]                                   
          exc.seg$chrom<-all.na[all.na[,2],1]
          exc.seg$seg.mean<-NA
          seg.lair<-rbind(seg.lair,exc.seg)
          idx<-order(seg.lair$chrom,seg.lair$loc.start)
          seg.lair<-seg.lair[idx,]
        }
        # fill in missing probes (lots of them because lair only contains hetrezygote normal)
        seg.lair$loc.start[1]<-seg.smp$loc.start[1]
        prevstate<-seg.lair$chrom[1]
        for (lair.state in 2:nrow(seg.lair)) {
          if(seg.lair$chrom[lair.state]==prevstate) {
            # put non-overlapping split in middle between 2 segments
            seg.lair$loc.end[lair.state-1]<- (seg.lair$loc.end[lair.state-1]+seg.lair$loc.start[lair.state] )%/% 2
            seg.lair$loc.start[lair.state]<- seg.lair$loc.end[lair.state-1]+1
          }else{
            seg.lair$loc.end[lair.state-1]<-seg.smp$loc.end[prevstate]
            seg.lair$loc.start[lair.state]<-seg.smp$loc.start[seg.lair$chrom[lair.state]]
          }
          prevstate<-seg.lair$chrom[lair.state]
        }
        seg.lair$loc.end[nrow(seg.lair)]<-seg.smp$loc.end[prevstate]
        #
        for(state in 1:nrow(seg.lair)) {
          lair.states[states[,assay]==seg.lair$chrom[state] & maploc>=seg.lair$loc.start[state] & maploc<=seg.lair$loc.end[state],assay]<-state
          lair.predicted[states[,assay]==seg.lair$chrom[state] & maploc>=seg.lair$loc.start[state] & maploc<=seg.lair$loc.end[state],assay]<-seg.lair$seg.mean[state]
        }
        
      }
    }
    res<-object
    assayData(res)$observed<-observed # in case of transformations intensity slot is not enough
    if (useLair)
      assayData(res)$states<-lair.states
    else
      assayData(res)$states<-states
    assayData(res)$predicted<-predicted
    if (useLair) {
      assayData(res)$lair.predicted<-lair.predicted
    }
    res
  } else segmentate.old(object, method, normalizedTo, doLog, doMerge, subsample)
}

segmentate(sorteddata.seg,doLog=FALSE,useLair=TRUE)
