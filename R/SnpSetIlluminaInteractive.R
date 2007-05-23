# Interactive plots
plotGoldenGate4OPA<-function(object,cn.sum=NULL,array=1,plotRaw=FALSE,main=sampleNames(object)[array],...) {
  if (!all(c("lai","nor.gt","loh") %in% assayDataElementNames(object)))
    stop("'calculateLOH' should be performed before making this plot")
  if (!all(c("observed","states","predicted") %in% assayDataElementNames(object)))
    stop("'segmentate' should be performed before making this plot")
  subsample<-beadarraySNP:::getSubsample(object,"OPA")
  if (is.null(cn.sum)) cn.sum<-createCNSummary(object,1)
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
    if (plotRaw) points((1:selected.n)/selected.n,y.base+assayData(object)$observed[selection,array],pch=".")
    chr<-featureData(object)$CHR[selection]
    xax<-getMidMaxIdx(chr)
    text(xax$midpos/selected.n,y.base-1,row.names(xax),pos=3)
    axis(2,y.base+c(0.5,1,1.5),c(0.5,1,1.5))
    segments(xax$maxpos/selected.n,y.base-1,xax$maxpos/selected.n,y.base+2)
    abline(h=y.base+c(-0.25,0.25,1,2),col=c("black","black","grey","black"))
    predicted<-assayData(object)$predicted[selection,array]
    points((1:selected.n)/selected.n,y.base+predicted,pch="-",col="red")
    # cn from user
    st.sel<-which(cn.sum$states$opa == subsmp)
    st.pred<-cn.sum$states$intensity[st.sel]
    st.tcn<-cn.sum$states$copynumber[st.sel]
    trueCN<-st.tcn[match(assayData(object)$predicted[selection,array],st.pred)]
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
    het.nrm<-assayData(object)$nor.gt[selection,array]=="TRUE"
    het.nrm<-names(het.nrm)[het.nrm]
    het.nrm<-het.nrm[!is.na(het.nrm)]
    idx<-which(featureNames(object)[selection] %in% het.nrm)
    points(idx/selected.n,y.base+assayData(object)$lai[het.nrm,array]*0.5,pch="-",col="blue")
    q.col<-ifelse(assayData(object)$GSR[het.nrm,array]<0.8,"mediumblue","yellow")
    col<-ifelse(assayData(object)$call[het.nrm,array]=="H",q.col,"red")
    points(idx/selected.n,y.base+(assayData(object)$GSR[het.nrm,array]*0.25)-0.25,pch="-",col=col)

    y.base<-y.base-3
  }
  # plot info in lower panel
  abline(h=2)
  rect(c(0.9,0.9),c(0,1),c(1,1),c(1,2))
  text(0.9,c(0.5,1.5),c("Redraw","Raw data"),pos=4)
  text(0.02,0.4,paste("Target Index",cn.sum$dnaIndex,"Current index",getDNAindex(cn.sum)),adj=0)

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

createCNSummary<-function(object,smp,dnaIndex,subsample="OPA"){
  deftarget<-round(dnaIndex*2)
  subsample<-beadarraySNP:::getSubsample(object,subsample)
  res<-NULL
  for (subsmp in levels(subsample)) {
    selection<-subsample == subsmp
    tmp<-summary(as.factor(assayData(object)$predicted[selection,smp]))
    res<-rbind(res,data.frame(opa=rep(subsmp,length(tmp)),count=tmp,intensity=names(tmp)))
  }
  # determine total copynumber of normal individual in this dataset (each SNP has same weight)
  copynumber.total.nrm<-sum(numericCHR(featureData(object)$CHR)<90)*2
  if (!is.null(pData(object)$"Gender")) gender<-pData(object)[smp,"Gender"]
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
