# Interactive plots
interactiveCNselect<-function(object,sample=1,dnaIndex) {
  plotRaw<-TRUE
  cn.res<-createCNSummary(object,sample=sample,dnaIndex=dnaIndex)
  cn.res<-plotGoldenGate4OPA(object,cn.res,sample=sample,plotRaw=plotRaw)
  repeat {
    location<-locator(n=1)
    if (is.null(location)) break
    if (location$x>0.9 & location$y>0 & location$y<1) {
      cn.res<-plotGoldenGate4OPA(object,cn.res,sample=sample,plotRaw=plotRaw)
    } else if (location$x>0.9 & location$y>1 & location$y<2) {
      plotRaw<-!plotRaw
      cn.res<-plotGoldenGate4OPA(object,cn.res,sample=sample,plotRaw=plotRaw)
    } else {
      opa<-floor((17-location$y)/3)
      if (opa %in% 1:4) {
        snpid<-getIdGoldenGate4OPA(object,opa,location$x)
        increase<-((16-location$y) - opa*3) < 0
        # rect(0,0,0.9,2,col="white")
        cn.res<-alterCN(cn.res,opa,assayData(object)$predicted[snpid,sample],increase)
        cn.res<-plotGoldenGate4OPA(object,cn.res,sample=sample,plotRaw=plotRaw)
        text(0,1.6,paste("Chr",reporterInfo(object)[snpid,"CHR"],"position",reporterInfo(object)[snpid,"MapInfo"]),adj=0)
        text(0,1.2,paste("CN",assayData(object)$predicted[snpid,sample]),adj=0)
        text(0,0.8,paste("increase",increase),adj=0)
      }
    }
  }
  cn.res
}

plotGoldenGate4OPA<-function(object,cn.sum=NULL,sample=1,plotRaw=FALSE,main=NULL,...) {
  if (!all(c("lai","nor.gt","loh") %in% assayDataElementNames(object)))
    stop("'calculateLOH' should be performed before making this plot")
  if (!all(c("observed","states","predicted") %in% assayDataElementNames(object)))
    stop("'segmentate' should be performed before making this plot")
  subsample<-beadarraySNP:::getSubsample(object,"OPA")
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
  plot(0,type="n",main=main,ylab="",yaxt="none",xlab="",xaxt="none")
  par(usr=c(0,1,0,14))
  y.base<-12
  for (subsmp in levels(subsample)) {
    selection<-subsample == subsmp
    chrs<-summary(as.factor(featureData(object)$CHR[selection]))
    selection<-selection & featureData(object)$CHR %in% names(chrs[chrs>=10])
    selected.n<-sum(selection)
    if (plotRaw) points((1:selected.n)/selected.n,y.base+assayData(object)$observed[selection,sample],pch=".")
    chr<-featureData(object)$CHR[selection]
    xax<-getMidMaxIdx(chr)
    text(xax$midpos/selected.n,y.base-1,row.names(xax),pos=3)
    axis(2,y.base+c(0.5,1,1.5),c(0.5,1,1.5))
    segments(xax$maxpos/selected.n,y.base-1,xax$maxpos/selected.n,y.base+2)
    abline(h=y.base+c(-0.25,0.25,1,2),col=c("black","black","grey","black"))
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
    idx<-which(c(trueCN,trueCN[length(trueCN)])!=c(-100,trueCN))
    posi<-(idx+c(idx[-1],length(trueCN)+1))/2 -0.5
    text(posi/selected.n,y.base-0.7,trueCN[idx],pos=3,col="green" )
    segments(idx/selected.n,y.base-1,idx/selected.n,y.base-0.25,col="green")
    #
    het.nrm<-assayData(object)$nor.gt[selection,sample]=="TRUE"
    het.nrm<-names(het.nrm)[het.nrm]
    het.nrm<-het.nrm[!is.na(het.nrm)]
    idx<-which(featureNames(object)[selection] %in% het.nrm)
    points(idx/selected.n,y.base+assayData(object)$lai[het.nrm,sample]*0.5,pch="-",col="blue")
    q.col<-ifelse(assayData(object)$GSR[het.nrm,sample]<0.8,"mediumblue","yellow")
    col<-ifelse(assayData(object)$call[het.nrm,sample]=="H",q.col,"red")
    points(idx/selected.n,y.base+(assayData(object)$GSR[het.nrm,sample]*0.25)-0.25,pch="-",col=col)

    y.base<-y.base-3
  }
  # plot info in lower panel
  abline(h=2)
  rect(c(0.9,0.9),c(0,1),c(1,1),c(1,2))
  text(0.9,c(0.5,1.5),c("Redraw","Raw data"),pos=4)
  text(0.02,0.4,paste("Target Index",cn.sum$dnaIndex,"Current index",format(getDNAindex(cn.sum))),adj=0)

  invisible(cn.sum)
}

getIdGoldenGate4OPA<-function(object,opa,xpos) {
  subsample<-beadarraySNP:::getSubsample(object,"OPA")
  subsmp<-levels(subsample)[opa]
  selection<-subsample == subsmp
  chrs<-summary(as.factor(featureData(object)$CHR[selection]))
  selection<-featureNames(object)[selection & featureData(object)$CHR %in% names(chrs[chrs>=10])]
  selection[ceiling(xpos*length(selection))]
}

createCNSummary<-function(object,sample,dnaIndex,subsample="OPA"){
  deftarget<-round(dnaIndex*2)
  subsample<-beadarraySNP:::getSubsample(object,subsample)
  res<-NULL
  for (subsmp in levels(subsample)) {
    selection<-subsample == subsmp
    tmp<-summary(as.factor(assayData(object)$predicted[selection,sample]))
    res<-rbind(res,data.frame(opa=rep(subsmp,length(tmp)),count=tmp,intensity=names(tmp)))
  }
  # determine total copynumber of normal individual in this dataset (each SNP has same weight)
  copynumber.total.nrm<-sum(numericCHR(featureData(object)$CHR)<90)*2
  if (!is.null(pData(object)$"Gender")) gender<-pData(object)[sample,"Gender"]
  else gender<-"F"
  if (gender=="F") copynumber.total.nrm<-copynumber.total.nrm + sum(numericCHR(featureData(object)$CHR)==98)*2 else
     copynumber.total.nrm<-copynumber.total.nrm + sum(numericCHR(featureData(object)$CHR) %in% 98:99)
  res$copynumber<-deftarget
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
  subsample<-beadarraySNP:::getSubsample(object,"OPA")
  for (subsmp in levels(subsample)) {
    selection<-subsample == subsmp
    # cn from user
    st.sel<-which(cn.sum$states$opa == subsmp)
    st.pred<-cn.sum$states$intensity[st.sel]
    st.tcn<-cn.sum$states$copynumber[st.sel]
    trueCN<-st.tcn[match(assayData(object)$predicted[selection,sample],st.pred)]
    inferred[selection,smp]<-trueCN
  }
  res<-object
  assayData(res)$inferred<-inferred
  res
}