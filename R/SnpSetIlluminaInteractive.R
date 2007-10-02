# Interactive plots
interactiveCNselect<-function(object,sample=1,dnaIndex) {

  getIdGoldenGate4OPA<-function(opa,xpos) {
    # retriev the probe that belongs to a certain position within the plot
    subsample<-getSubsample(object,"OPA")
    subsmp<-levels(subsample)[opa]
    selection<-subsample == subsmp
    chrs<-summary(as.factor(featureData(object)$CHR[selection]))
    selection<-featureNames(object)[selection & featureData(object)$CHR %in% names(chrs[chrs>=10])]
    selection[ceiling(xpos*length(selection))]
  }

  plotRaw<-TRUE
  cn.res<-createCNSummary(object,sample=sample,dnaIndex=dnaIndex)
  cn.res<-plotGoldenGate4OPA(object,cn.res,sample=sample,plotRaw=plotRaw,interact=TRUE)
  repeat {
    location<-locator(n=1)
    if (is.null(location)) break
    if (location$x>0.9 & location$y>0 & location$y<1) {
      cn.res<-plotGoldenGate4OPA(object,cn.res,sample=sample,plotRaw=plotRaw,interact=TRUE)
    } else if (location$x>0.9 & location$y>1 & location$y<2) {
      plotRaw<-!plotRaw
      cn.res<-plotGoldenGate4OPA(object,cn.res,sample=sample,plotRaw=plotRaw,interact=TRUE)
    } else {
      opa<-floor((17-location$y)/3)
      if (opa %in% 1:4) {
        snpid<-getIdGoldenGate4OPA(opa,location$x)
        increase<-((16-location$y) - opa*3) < 0
        # rect(0,0,0.9,2,col="white")
        cn.res<-alterCN(cn.res,opa,assayData(object)$predicted[snpid,sample],increase)
        cn.res<-plotGoldenGate4OPA(object,cn.res,sample=sample,plotRaw=plotRaw,interact=TRUE)
        text(0,1.6,paste("Chr",reporterInfo(object)[snpid,"CHR"],"position",reporterInfo(object)[snpid,"MapInfo"]),adj=0)
        text(0,1.2,paste("CN",assayData(object)$predicted[snpid,sample]),adj=0)
        text(0,0.8,paste("increase",increase),adj=0)
      }
    }
  }
  cn.res
}

plotGoldenGate4OPA<-function(object,cn.sum=NULL,sample=1,plotRaw=FALSE,main=NULL,interact=FALSE,...) {
  if (!all(c("lair","nor.gt","loh") %in% assayDataElementNames(object)))
    stop("'calculateLOH' should be performed before making this plot")
  if (!all(c("observed","states","predicted") %in% assayDataElementNames(object)))
    stop("'segmentate' should be performed before making this plot")
  subsample<-getSubsample(object,"OPA")
  if (length(sample)!=1) stop("Only 1 sample can be handled at a time")
  if (is.numeric(sample)) {
    if(sample<1 | sample>ncol(object)) stop(paste("sample",sample,"is not between 1 and ",ncol(object)))
  } else {
    if (!(sample %in% sampleNames(object))) stop(paste(sample,"is not a sample name in the object"))
  }
  if (is.null(cn.sum)) cn.sum<-createCNSummary(object,1)
  if (is.null(main)) {
    if (is.numeric(sample)) main=sampleNames(object)[sample]
    else main=sample
  }
  # setup plot regions
  par(mar=c(1,3,3,1))
  plot(0,type="n",main=main,ylab="",yaxt="none",xlab="",xaxt="none",...)
  if (interact) y.bottom<-0 else y.bottom<-1.2
  par(usr=c(0,1,y.bottom,14))
  y.base<-12
  for (subsmp in levels(subsample)) {
    selection<-subsample == subsmp
    chrs<-summary(as.factor(featureData(object)$CHR[selection]))
    selection<-selection & featureData(object)$CHR %in% names(chrs[chrs>=10])
    selected.n<-sum(selection)
    # Horizontal section in subplot
    abline(h=y.base+c(-0.25,0.25,1,2),col=c("black","black","grey","black"))
    #
    if (plotRaw) points((1:selected.n)/selected.n,y.base+assayData(object)$observed[selection,sample],pch=".")
    # Chromosome annotation
    chr<-featureData(object)$CHR[selection]
    xax<-getMidMaxIdx(chr)
    text(xax$midpos/selected.n,y.base-1,row.names(xax),pos=3)
    axis(2,y.base+c(0.5,1,1.5),c(0.5,1,1.5))
    segments(xax$maxpos/selected.n,y.base-1,xax$maxpos/selected.n,y.base+2)
    # Values from segmentation function
    predicted<-assayData(object)$predicted[selection,sample]
    points((1:selected.n)/selected.n,y.base+predicted,pch="-",col="red")
    # cn from user
    st.sel<-which(cn.sum$states$opa == subsmp)
    st.pred<-cn.sum$states$intensity[st.sel]
    st.tcn<-cn.sum$states$copynumber[st.sel]
    trueCN<-st.tcn[match(assayData(object)$predicted[selection,sample],st.pred)]
    for (i in unique(trueCN)) {
      probes<-trueCN == i
      predicted[probes]<-mean(predicted[probes])
    }
    points((1:selected.n)/selected.n,y.base+predicted,pch="-",col="green")
    # Annotation of segments. Also subdivide at chromosomal boundaries
    idx<-which(c(trueCN,trueCN[length(trueCN)])!=c(-100,trueCN) | c(chr[1],chr)!=c(chr,chr[length(chr)]))
    posi<-(idx+c(idx[-1],length(trueCN)+1))/2 -0.5
    text(posi/selected.n,y.base-0.7,trueCN[idx],pos=3,col="green" )
    segments(idx/selected.n,y.base-1,idx/selected.n,y.base-0.25,col="green")
    # lesser allele intensity ratio
    lair.offset<- -0.15
    lair.range<- 0.40
    het.nrm<-assayData(object)$nor.gt[selection,sample]=="TRUE"
    het.nrm<-names(het.nrm)[het.nrm]
    het.nrm<-het.nrm[!is.na(het.nrm)]
    idx<-which(featureNames(object)[selection] %in% het.nrm)
    points(idx/selected.n,y.base+lair.offset+assayData(object)$lair[het.nrm,sample]*lair.range,pch="-",col="blue")
    # LOH + quality
    loh.offset<- -0.25
    loh.range<- 0.10
    loh.width<- 1.5
    q.col<-ifelse(assayData(object)$GSR[het.nrm,sample]<0.8,"mediumblue","green")
    col<-ifelse(assayData(object)$call[het.nrm,sample]=="H",q.col,"red")
    segments(idx/selected.n,y.base+loh.offset,idx/selected.n,y.base+loh.offset+loh.range,lwd=loh.width,col=col)
    #
    y.base<-y.base-3
  }
  # plot info in lower panel
  abline(h=2)
  if (interact) {
    rect(c(0.9,0.9),c(0,1),c(1,1),c(1,2))
    text(0.9,c(0.5,1.5),c("Redraw","Raw data"),pos=4)
  }
  text(0.02,y.bottom+0.4,paste("Target Index",cn.sum$dnaIndex,"Current index",format(getDNAindex(cn.sum),digits=2,nsmall=2)),adj=0)
  invisible(cn.sum)
}

createCNSummary<-function(object,sample,dnaIndex=1,subsample="OPA"){
  deftarget<-round(dnaIndex*2)
  subsample<-getSubsample(object,subsample)
  res<-NULL
  getExisting<-FALSE
  if (!is.null(assayData(object)$inferred)) if(!any(is.na(assayData(object)$inferred[,sample]))) getExisting<-TRUE
  for (subsmp in levels(subsample)) {
    selection<-subsample == subsmp
    tmp<-summary(as.factor(assayData(object)$predicted[selection,sample]))
    if (getExisting) deftarget<-aggregate(assayData(object)$inferred[selection,sample],by=list(assayData(object)$predicted[selection,sample]),FUN=mean)[,2]
    res<-rbind(res,data.frame(opa=rep(subsmp,length(tmp)),count=tmp,intensity=names(tmp),copynumber=deftarget))
  }
  # determine total copynumber of normal individual in this dataset (each SNP has same weight)
  copynumber.total.nrm<-sum(numericCHR(featureData(object)$CHR)<90)*2
  if (!is.null(pData(object)$"Gender")) gender<-pData(object)[sample,"Gender"]
  else gender<-"F"
  if (gender=="F") copynumber.total.nrm<-copynumber.total.nrm + sum(numericCHR(featureData(object)$CHR)==98)*2 else
     copynumber.total.nrm<-copynumber.total.nrm + sum(numericCHR(featureData(object)$CHR) %in% 98:99)
  list(dnaIndex=dnaIndex,CN.total.nrm=copynumber.total.nrm,states=res)
}

alterCN<-function(cn.sum,opa,value,updown) {
  selection<-which(cn.sum$states$opa == levels(cn.sum$states$opa)[opa])
  idx<-which(cn.sum$states$intensity[selection]==value)
  if (idx>0) {
    if (updown) {
      cn.sum$states$copynumber[selection[idx]]<-cn.sum$states$copynumber[selection[idx]]+1
      for (i in idx:length(selection))
        if (cn.sum$states$copynumber[selection[i]]<cn.sum$states$copynumber[selection[idx]])
          cn.sum$states$copynumber[selection[i]]<-cn.sum$states$copynumber[selection[idx]]
    } else {
      cn.sum$states$copynumber[selection[idx]]<-cn.sum$states$copynumber[selection[idx]]-1
      for (i in 1:idx)
        if (cn.sum$states$copynumber[selection[i]]>cn.sum$states$copynumber[selection[idx]])
          cn.sum$states$copynumber[selection[i]]<-cn.sum$states$copynumber[selection[idx]]
    }
  }
  cn.sum
}

getDNAindex<-function(cn.sum) {
  sum(cn.sum$states$count*cn.sum$states$copynumber)/cn.sum$CN.total.nrm
}

setRealCN<-function(object,sample,cn.sum) {
  if (!("inferred" %in% assayDataElementNames(object))) {
    inferred<-matrix(NA,nrow=nrow(object),ncol=ncol(object),dimnames=dimnames(assayData(object)$G))
  } else {
    inferred<-assayData(object)$inferred
  }
  subsample<-getSubsample(object,"OPA")
  for (subsmp in levels(subsample)) {
    selection<-subsample == subsmp
    # cn from user
    st.sel<-which(cn.sum$states$opa == subsmp)
    st.pred<-cn.sum$states$intensity[st.sel]
    st.tcn<-cn.sum$states$copynumber[st.sel]
    trueCN<-st.tcn[match(assayData(object)$predicted[selection,sample],st.pred)]
    inferred[selection,sample]<-trueCN
  }
  res<-object
  assayData(res)$inferred<-inferred
  res
}