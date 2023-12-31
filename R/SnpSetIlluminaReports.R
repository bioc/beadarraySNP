## snp reports
## functions to report for all available data in a SnpSetIllumina object
reportSamplesSmoothCopyNumber<-function(snpdata, grouping=NULL, normalizedTo=2, smooth.lambda=2, ridge.kappa=0, 
    plotLOH=c("none","marker","line","NorTum"), sample.colors=NULL, ...){
  # default grouping is by 4 in sequence of samples in snpdata
  plotLOH<-match.arg(plotLOH)
  grouping<-getGrouping(snpdata,grouping,by=4)
  # make sure chromosmes are sorted
  CHR<-numericCHR(pData(featureData(snpdata))[,"CHR"])
  ind<-order(CHR,pData(featureData(snpdata))[,"MapInfo"])
  snpdata<-snpdata[ind,]
  CHR<-CHR[ind]
  chroms<-unique(CHR)

  if (is.null(sample.colors)) sample.colors<-c("red","green","blue","orange","brown","turquoise","yellow","purple","pink","magenta")
  
  intensities<-assayData(snpdata)[["intensity"]]
  for (pageID in levels(factor(grouping))){
    samples<-sampleNames(snpdata)[grouping == pageID]
    if (length(samples)>0) {
      dchrompos<-prepareGenomePlot(pData(featureData(snpdata))[,c("CHR","MapInfo")],...)
      marker.width<-(par("usr")[2]-par("usr")[1])/1000
      tum.n<-sum(pData(snpdata)[samples,"NorTum"]!="N")
      tum.i<-0
      for (i1 in 1:length(samples) ) {
        if(tum.now <- (pData(snpdata)[samples[i1],"NorTum"]!="N")) tum.i<-tum.i+1
        
        for (chrom in chroms) {
          probes<- CHR == chrom
          if (sum(!is.na(intensities[probes,samples[i1]]))>10 ) {
            smoothed<-quantsmooth(intensities[probes,samples[i1]],smooth.lambda=smooth.lambda,ridge.kappa=ridge.kappa,segment=151)
            lines(dchrompos[probes,2],dchrompos[probes,1]+(smoothed-normalizedTo)/normalizedTo,col=sample.colors[i1],lwd=1.5)
            if (plotLOH!="none") {
              probeNames<-featureNames(snpdata)[probes]
              if (plotLOH=="marker") {
                if ("nor.gt" %in% assayDataElementNames(snpdata)) {
                  if (tum.now) {
                    het.nrm<-assayData(snpdata)$nor.gt[probes,samples[i1]]=="TRUE"
                    het.nrm<-names(het.nrm)[het.nrm]
                    het.nrm<-het.nrm[!is.na(het.nrm)]
                    if (length(het.nrm)>0) {
                      idx<-match(het.nrm,featureNames(snpdata))
                      # LOH + quality
                      loh.range<- 0.15/tum.n
                      loh.offset<- -0.62 + (loh.range+0.02)*tum.i
                      loh.width<- 1.5
                      ypos<-mean(dchrompos[idx,1],na.rm=TRUE)+loh.offset+loh.range/2
                      lines(c(min(dchrompos[probes,2],na.rm=TRUE),max(dchrompos[probes,2],na.rm=TRUE)),c(ypos,ypos),lwd=loh.width,col=sample.colors[i1])
                      q.col<-ifelse(assayData(snpdata)$GSR[idx,samples[i1]]<0.8,"mediumblue","green")
                      q.col<-ifelse(assayData(snpdata)$call[idx,samples[i1]]=="H",q.col,"red")
                      rect(dchrompos[idx,2]-marker.width,dchrompos[idx,1]+loh.offset,dchrompos[idx,2]+marker.width,dchrompos[idx,1]+loh.offset+loh.range,border=NA,col=q.col)
                    }
                  }
                } else {
                  chromhet<-heterozygosity(assayData(snpdata)[["call"]][probes,samples[i1]])
                  LOH<-probeNames[chromhet>20]
                  if (length(LOH)>0) points(dchrompos[LOH,2],dchrompos[LOH,1]-0.3-(i1*0.05),pch="-",col=sample.colors[i1])
                }
              }
              if (plotLOH=="line") {
                chromhet<-heterozygosity(assayData(snpdata)[["call"]][probes,samples[i1]])
                lines(dchrompos[probes,2],dchrompos[probes,1]+scaleto(chromhet,c(10,40),c(0.1,-0.4)),col=sample.colors[i1],lty=2)
              }
              if (plotLOH=="NorTum" && tum.now) {
                ## check availability of normal to compare with
                n1<-grep("N",as.character(pData(snpdata)[samples,"NorTum"]))
                if (length(n1)>0) {
                  compGenotype<-compareGenotypes(assayData(snpdata)[["call"]][probes,samples[i1]],assayData(snpdata)[["call"]][probes,samples[n1[1]]])
                  LOH<-probeNames[compGenotype=="l"] # loss 
                  if(length(LOH)>0) points(dchrompos[LOH,2],dchrompos[LOH,1]-0.3-(i1*0.05),pch="|",col=sample.colors[i1])
                  LOH<-probeNames[compGenotype=="i"] # heterozygous normal
                  if(length(LOH)>0) points(dchrompos[LOH,2],dchrompos[LOH,1]-0.3-(i1*0.05),pch="'",col=sample.colors[i1])
                }
              }
            }
          }
        }
      }
      legend("topleft",samples,col=c(sample.colors[1:length(samples)]),pch=18,ncol=3)
    }
  }
}
pdfSamplesSmoothCopyNumber<-function(object,filename,...) {
  pdf(filename,width=7.2,height=11)
	reportSamplesSmoothCopyNumber(object,...)
	dev.off() 
}

reportChromosomesSmoothCopyNumber<-function(snpdata, grouping=NULL, normalizedTo=2, smooth.lambda=2, ridge.kappa=0, 
    plotLOH=c("none","marker","line","NorTum"), sample.colors=NULL, ideo.bleach=0.25, ...){
  ideo.width<-0.15
  ideo.ypos<-normalizedTo+(ideo.width/2)
  plotLOH<-match.arg(plotLOH)
  grouping<-getGrouping(snpdata,grouping,by=4)
  if (is.null(sample.colors)) sample.colors<-c("red","green","blue","orange","brown","turquoise","yellow","purple","pink","magenta")
  # make sure chromosmes are sorted
  CHR<-numericCHR(pData(featureData(snpdata))[,"CHR"])
  ind<-order(CHR,pData(featureData(snpdata))[,"MapInfo"])
  snpdata<-snpdata[ind,]
  CHR<-CHR[ind]
  chroms<-unique(CHR)

  intensities<-assayData(snpdata)[["intensity"]]
  for (group in levels(factor(grouping))){
    samples<-sampleNames(snpdata)[grouping == group]
    if (length(samples)>0) {
      par(mfrow=par()$mfrow) # Start on a new page
      for (chrom in chroms) {
        probes<- CHR == chrom
        if (any(apply(intensities[probes,samples,drop=FALSE],2,function(x) sum(!is.na(x))>10))) {
          plot(c(0,max(lengthChromosome(chrom,"bases"),pData(featureData(snpdata))[probes,"MapInfo"])),c(1,3),main=paste(group,"Chromosome",characterCHR(chrom)),type="n",ylab="intensity",xlab="",xaxt="n")
          paintCytobands(chrom,pos=c(0,ideo.ypos),units="bases",width=ideo.width,legend=FALSE,bleach=ideo.bleach)
          plotSmoothed(intensities[probes,samples,drop=FALSE],pData(featureData(snpdata))[probes,"MapInfo"],smooth.lambda=smooth.lambda,plotnew=FALSE,cols=sample.colors,...)
          legend("topleft",samples,col=sample.colors[1:length(samples)],lty=1,lwd=2,ncol=length(samples))
          if (plotLOH!="none") {
            probeNames<-featureNames(snpdata)[probes]
            markerbase<-par("yaxp")[1]
            markerinterval<-(normalizedTo-markerbase)/13
            if (plotLOH=="marker") {
              if ("nor.gt" %in% assayDataElementNames(snpdata)) {
                tum.i<- 0
                marker.width<-(par("usr")[2]-par("usr")[1])/1000
                for (samp in 1:length(samples)) {
                  if (pData(snpdata)[samples[samp],"NorTum"]!="N") {
                    tum.i<-tum.i+1
                    het.nrm<-assayData(snpdata)$nor.gt[probes,samples[samp]]=="TRUE"
                    het.nrm<-names(het.nrm)[het.nrm]
                    het.nrm<-het.nrm[!is.na(het.nrm)]
                    if (length(het.nrm)>0) {
                      idx<-match(het.nrm,featureNames(snpdata))
                      # LOH + quality
                      abline(h=markerbase+tum.i*markerinterval,lwd=1.5,col=sample.colors[samp])
                      q.col<-ifelse(assayData(snpdata)$GSR[idx,samples[samp]]<0.8,"mediumblue","green")
                      q.col<-ifelse(assayData(snpdata)$call[idx,samples[samp]]=="H",q.col,"red")
                      xpos<-pData(featureData(snpdata))[idx,"MapInfo"]
                      rect(xpos-marker.width,markerbase+(tum.i-0.4)*markerinterval,xpos+marker.width,markerbase+(tum.i+0.4)*markerinterval,border=NA,col=q.col)
                    }
                    
                  }
                }
              } else {
                for (samp in 1:length(samples)) {
                  chromhet<-heterozygosity(assayData(snpdata)[["call"]][probes,samples[samp]])
                  LOH<-probeNames[chromhet>20]
                  if (length(LOH)>0) points(pData(featureData(snpdata))[LOH,"MapInfo"],markerbase+(samp*markerinterval),pch="-",col=sample.colors[samp])
                }
              }
            }
            if (plotLOH=="line") {
              for (samp in 1:length(samples)) {
                chromhet<-heterozygosity(assayData(snpdata)[["call"]][probes,samples[samp]])
                lines(pData(featureData(snpdata))[probes,"MapInfo"],scaleto(chromhet,c(10,40),c(normalizedTo,markerbase)),col=sample.colors[samp],lty=2)
              }
            }
            if (plotLOH=="NorTum") {
              ## check availability of normal to compare with
              n1<-grep("N",as.character(pData(snpdata)[samples,"NorTum"]))
              t1<-samples[-n1]
              if (length(n1)>0 & length(t1)>0) { # at least 1 normal and 1 tumor sample in group
                for (tum in 1:length(t1)) {
                  compGenotype<-compareGenotypes(assayData(snpdata)[["call"]][probes,t1[tum]],assayData(snpdata)[["call"]][probes,samples[n1[1]]])
                  LOH<-probeNames[compGenotype=="l"] # loss
                  if(length(LOH)>0) points(pData(featureData(snpdata))[LOH,"MapInfo"],rep(markerbase+(tum*markerinterval),length(LOH)),pch="|",col=sample.colors[match(t1[tum],samples)])
                  LOH<-probeNames[compGenotype=="i"] # heterozygous normal)
                  if(length(LOH)>0) points(pData(featureData(snpdata))[LOH,"MapInfo"],rep(markerbase+(tum*markerinterval),length(LOH)),pch="'",col=sample.colors[match(t1[tum],samples)])
                }
              }
            }
          }
        }
      }
    }
  }
}

pdfChromosomesSmoothCopyNumber<-function(object,filename,...) {
  pdf(filename,width=7.2,height=11)
  par(mfrow=c(4,1),mar=c(1.5,2,2,0))
	reportChromosomesSmoothCopyNumber(object,...)
	dev.off() 
}

getMidMaxIdx<-function(groups){
  groups<-as.character(groups)
  lvls<-levels(factor(groups))
  midpos<-NULL
  maxpos<-NULL
  for(lvl in 1:length(lvls)) {
    midpos[lvl]<-mean(grep(paste("^",lvls[lvl],"$",sep=""),groups))-0.5
    maxpos[lvl]<-max(grep(paste("^",lvls[lvl],"$",sep=""),groups))
  }
  idx<-order(midpos)
  data.frame(midpos,maxpos,row.names=lvls)[idx,]
}

reportGenomeGainLossLOH<-function(snpdata,grouping=NULL,plotSampleNames=FALSE,sizeSampleNames=4,distance.min,
  upcolor="red",downcolor="blue",lohcolor="grey",hetcolor="lightgrey",lohwidth=1,segment=101,
  orientation=c("V","H"),...) {
  orientation<-match.arg(orientation)
  ind<-order(numericCHR(fData(snpdata)$CHR),fData(snpdata)$MapInfo)
  snpdata<-snpdata[ind,]
  grouping<-getGrouping(snpdata,grouping)
  ind<-order(grouping)
  snpdata<-snpdata[,ind]
  grouping<-grouping[ind]
  if (missing(distance.min)) distance.min=1e+9
  
  # determine gained and lost probes
  if (orientation=="V") ylim<-c(nrow(snpdata),0)  else ylim<-c(0,nrow(snpdata))
  plot(0,xlim=c(0,ncol(snpdata)),ylim=ylim,type="n",xaxt="n",yaxt="n",xlab="",ylab="")
  par(usr=c(0,ncol(snpdata),ylim))
  for (smp in 1:ncol(snpdata)) {
    regions<-getChangedRegions(assayData(snpdata)$intensity[,smp],segment=segment,...)
    if (!is.null(regions)) rect(smp-1,regions[,"start"]-1,smp-0.5,regions[,"end"],col=ifelse(regions[,"up"],upcolor,downcolor),border=NA)
    het<-which(assayData(snpdata)$call[,smp]=="H")
    if (length(het)>0){
      rect(smp-0.25,het-1-lohwidth,smp,het+lohwidth,col=hetcolor,border=NA)
    }

    loh<-which(assayData(snpdata)$loh[,smp])
    if (length(loh)>0){
      position<-pData(featureData(snpdata))$MapInfo[loh]
      distance<-abs(c(position,0)-c(0,position))
      min.distance<-apply(cbind(distance[-1],distance[-length(distance)]),1,min)
      loh<-loh[min.distance<distance.min]
      if (length(loh)>0) rect(smp-0.5,loh-1-lohwidth,smp,loh+lohwidth,col=lohcolor,border=NA)
    }
    #
    #pLOH<-ifelse(assayData(snpdata)$call[,smp]=="H",(1-assayData(snpdata)$GRS[,smp]),1)
    #het.nrm<-which(assayData(snpdata)$nor.gt[,smp]=="H") # & assayData(snpdata)$nor.qs[,smp]>0.9
    #points(smp-0.5*pLOH[het.nrm],het.nrm,pch=".")
    #points(smp-0.5*assayData(snpdata)$nor.qs[het.nrm,smp],het.nrm,pch=".",col="cyan")

  }
  abline(v=1:(ncol(snpdata)-1),col="grey")
  samplenameAxis<-ifelse(orientation=="V",1,1)
  groupingAxis<-ifelse(orientation=="V",3,1)
  groupingLas<-ifelse(orientation=="V",0,2)
  groupingLine<-ifelse(orientation=="V",NA,ifelse(plotSampleNames,sizeSampleNames,NA))
  if (!missing(grouping)) {
    xax<-getMidMaxIdx(grouping)
    axis(groupingAxis,xax$midpos,row.names(xax),line=groupingLine,las=groupingLas)
    axis(groupingAxis,c(0,xax$maxpos),rep("",length(xax$maxpos)+1),line=groupingLine,las=groupingLas,tcl=0.25)
    abline(v=xax$maxpos)
  }

  if (plotSampleNames) {
    axis(samplenameAxis,(1:ncol(snpdata))-0.5,sampleNames(snpdata),las=2,cex.axis=0.6)
  }
  yax<-getMidMaxIdx(fData(snpdata)$CHR)
  axis(2,yax$midpos,row.names(yax))
  abline(h=yax$maxpos)

}

reportChromosomeGainLossLOH<-function(snpdata,grouping=NULL,plotSampleNames=FALSE,distance.min,
  upcolor="red",downcolor="blue",lohcolor="grey",hetcolor="lightgrey",proportion=0.2,plotLOH=TRUE,
  segment=101,...) {
  ind<-order(numericCHR(fData(snpdata)$CHR),fData(snpdata)$MapInfo)
  if (missing(distance.min)) distance.min=1e+9
  snpdata<-snpdata[ind,]
  grouping<-getGrouping(snpdata,grouping)
  ind<-order(grouping)
  snpdata<-snpdata[,ind]
  grouping<-grouping[ind]
  xmin<- -ncol(snpdata)*proportion
  cb.x<-xmin*0.6
  cb.w<- -xmin*0.2
  if (plotLOH) cn.w<-0.5 else cn.w<-1
  for (chrom in 1:22) {
    probes<-fData(snpdata)$CHR == chrom
    lengthchrom<-max(fData(snpdata)$MapInfo[probes],na.rm=TRUE)
    plot(0,xlim=c(xmin,ncol(snpdata)),ylim=c(0,lengthchrom),xlab="",ylab="",xaxt="n",yaxt="n",main=paste("chromosome",chrom),type="n")
    myusr<-par()$usr
    myusr[1]<-xmin
    myusr[2]<-ncol(snpdata)
    par(usr=myusr)
    paintCytobands(chrom,c(cb.x,lengthchrom),units="bases",width=cb.w,orientation="v",legend=TRUE)
    for (smp in 1:ncol(snpdata)) {
      updown<-getChangedRegions(assayData(snpdata)$intensity[probes,smp],fData(snpdata)$MapInfo[probes],
                                segment=segment,...)
      if (!is.null(updown)) {
        rect(smp-1,lengthchrom-updown[,"start"],smp-1+cn.w,lengthchrom-updown[,"end"],col=ifelse(updown[,"up"],upcolor,downcolor),border=NA)
      }
      probe.w<-lengthchrom/1000
      if (plotLOH) {
        gt<-assayData(snpdata)$call[probes,smp]
        het<-names(gt)[gt=="H"]
        if (length(het)>0){
          rect(smp-0.25,lengthchrom-fData(snpdata)[het,"MapInfo"]-probe.w,smp,lengthchrom-fData(snpdata)[het,"MapInfo"]+probe.w,col=hetcolor,border=NA)
        }
        loh<-featureNames(snpdata)[assayData(snpdata)$loh[,smp] & probes]
        if (length(loh)>0) {
          position<-pData(featureData(snpdata))[loh,"MapInfo"]
          distance<-abs(c(position,0)-c(0,position))
          min.distance<-apply(cbind(distance[-1],distance[-length(distance)]),1,min)
          loh<-loh[min.distance<distance.min]
          if (length(loh)>0) rect(smp-0.5,lengthchrom-fData(snpdata)[loh,"MapInfo"]-probe.w,smp,lengthchrom-fData(snpdata)[loh,"MapInfo"]+probe.w,col=lohcolor,border=NA)
        }
      }
    }
    if (plotSampleNames) {
      axis(1,(1:ncol(snpdata))-0.5,sampleNames(snpdata),las=2,cex.axis=0.6)
    }
    abline(v=0:(ncol(snpdata)-1),col="grey")
    if (!missing(grouping)) {
      xax<-getMidMaxIdx(grouping)
      axis(3,xax$midpos,row.names(xax))
      abline(v=xax$maxpos)
    }
  }
}

pdfChromosomeGainLossLOH<-function(object,filename,...) {
  pdf(filename,width=7.2,height=11)
  par(mfrow=c(4,1),mar=c(1.5,2,2,0))
  reportChromosomeGainLossLOH(object,...)
  dev.off()
}

reportGenomeIntensityPlot<-function(snpdata,normalizedTo=NULL,subsample=NULL,smoothing=c("mean","quant"),dot.col="black",smooth.col="red",...) {
  if ( !("intensity" %in% assayDataElementNames(snpdata))) snpdata<-RG2polar(snpdata)
  if (is.null(subsample)) {
    ind<-order(numericCHR(pData(featureData(snpdata))[,"CHR"]),pData(featureData(snpdata))[,"MapInfo"])
    snpdata<-snpdata[ind,]
    # use chromosomes as subsamples
    subsample<-numericCHR(pData(featureData(snpdata))[,"CHR"])
  } else {
    subsample<-getSubsample(snpdata,subsample)
    ind<-order(subsample)
    subsample<-subsample[ind]
    snpdata<-snpdata[ind,]
  }
  # find boundaries between subsamples
  xax<-getMidMaxIdx(subsample)
  if (smoothing=="mean") {
    avg.intensity<-aggregate(assayData(snpdata)$intensity,by=list(subsample),FUN=mean,na.rm=TRUE)
    idx<-match(rownames(xax),avg.intensity[,1])
    avg.intensity<-avg.intensity[idx,]
  }
  for (sample in 1:length(sampleNames(snpdata))) {
    plot(1,ylab="intensity",xlab="",main=sampleNames(snpdata)[sample],xaxt="n",yaxt="n",type="n")
    par(usr=c(1,dim(snpdata)[1],min(assayData(snpdata)[["intensity"]][,sample],na.rm=TRUE),max(assayData(snpdata)[["intensity"]][,sample],na.rm=TRUE)))
    if (!is.null(normalizedTo)) abline(h=normalizedTo,col="grey")
    axis(1,xax$midpos,row.names(xax))
    axis(2) # par(usr) does not recompute/relocate axis
    abline(v=xax$maxpos)
    points(assayData(snpdata)[["intensity"]][,sample],col=dot.col,pch=".",xlab="",...)
    if (smoothing=="mean") {
      segments(c(0,xax$maxpos[-length(xax$maxpos)]),avg.intensity[,sample+1],xax$maxpos,avg.intensity[,sample+1],col=smooth.col)
    } else if (smoothing=="quant") {
      lines(1:dim(snpdata)[1],quantsmooth(assayData(snpdata)[["intensity"]][,sample],smooth.lambda=4,segment=101),col=smooth.col)
    }
  }
}

reportGenotypeSegmentation<-function(object,plotRaw=TRUE,subsample=NULL,panels=0,minProbes=10,maxY=2,...) {
  if (!all(c("lair","nor.gt","loh") %in% assayDataElementNames(object)))
    stop("'calculateLOH' should be performed before making this report")
  if (!all(c("observed","states","predicted") %in% assayDataElementNames(object)))
    stop("'segmentate' should be performed before making this report")
  subsample<-getSubsample(object,subsample)
  if (panels == 0)
    panels<-length(levels(subsample))
  par(mfrow=c(panels,1),mar=c(2.6,4.1,3.6,1.6))
  for (smp in 1:ncol(object)) {
    for (subsmp in levels(subsample)) {
      selection<-subsample == subsmp
      chrs<-summary(as.factor(featureData(object)$CHR[selection]))
      selection<-selection & featureData(object)$CHR %in% names(chrs[chrs>=10])      
      if (is.null(maxY) | maxY<0) maxY.opa<-max(assayData(object)$observed[selection,smp],na.rm=TRUE)
      else maxY.opa<-maxY

      plot(0,type="n",main=sampleNames(object)[smp],ylab="intensity",yaxt="none",xlab="",xaxt="none",...)
      par(usr=c(0,sum(selection),-0.25, maxY.opa))
      if (plotRaw) points(assayData(object)$observed[selection,smp],pch=".")
      chr<-featureData(object)$CHR[selection]
      xax<-getMidMaxIdx(chr)
      axis(1,xax$midpos,row.names(xax))
      axis(2)
      abline(v=xax$maxpos)
      points(assayData(object)$predicted[selection,smp],pch="-",col="red")
      #
      # lesser allele intensity ratio
      lair.offset<- -0.15
      lair.range<- 0.40
      abline(h=lair.offset+lair.range)
      het.nrm<-assayData(object)$nor.gt[selection,smp]=="TRUE"
      het.nrm<-names(het.nrm)[het.nrm]
      het.nrm<-het.nrm[!is.na(het.nrm)]
      idx<-which(featureNames(object)[selection] %in% het.nrm)
      points(idx,lair.offset+assayData(object)$lair[het.nrm,smp]*lair.range,pch="_",col="blue")
      # 
      # LOH + quality
      loh.offset<- -0.24
      loh.range<- 0.09
      loh.width<- 1.5
      q.col<-ifelse(assayData(object)$GSR[het.nrm,smp]<0.8,"mediumblue","green")
      col<-ifelse(assayData(object)$call[het.nrm,smp]=="H",q.col,"red")
      segments(idx,loh.offset,idx,loh.offset+loh.range,lwd=loh.width,col=col)
    }
  }
  invisible()
}

plotGroupZygosity <- function(Green,Red,GenCall,Grouping,NorTum,NormalizedTo=1,...) {
  # Plot alleles. Green to X-axis, Red to Y-axis
  # GenCall: A,AA = -
  #          H,AB = +
  #          B,BB = |
  #          Not-recognized = .
  # Grouping : colors
  # NorTum : Size
  # NormalizedTo : <=0 Fitted to largest spread in Green/Red; >0 used to set limits (min(Red,Green), 2*Normalized)
  # ... : Transferred to plot()
  chars<-c("-","+","|",".")
  sample.colors<-rep(c("red","green","blue","yellow","orange","brown","purple","turquoise","pink","magenta"),length.out=length(Green))
  if (NormalizedTo<=0) plotsize<-c(min(Green,Red,na.rm=TRUE),max(Green,Red,na.rm=TRUE))
  else plotsize<-c(min(Green,Red,na.rm=TRUE),2*NormalizedTo)
  plot(plotsize,plotsize,type="n",cex.axis=0.6,...)
  char<-chars[(GenCall=="A" | GenCall=="AA") * 1 + (GenCall=="H" | GenCall=="AB")*2 + (GenCall=="B" | GenCall=="BB") * 3]
  char[char==0]<-4
  # now draw in colors for sampleids, shapes for snpcall, size for tum/norm (tum is smaller)
  charsize<-1.5 + ((NorTum=="N")*0.5)
  charcolor<-sample.colors[as.numeric(factor(Grouping))]
  for (chip in 1:length(char)) {
    points(Red[chip],Green[chip],col=charcolor[chip],pch=char[chip],cex=charsize[chip])
  }

}

reportGroupZygosity<-function (snpdata,snps,Grouping,NorTum,NormalizedTo=1){
  pdf("GroupZygosity.pdf",width=8,height=11)
  par(mfrow=c(4,3),mar=c(2,1,2,1))
  for (snp in snps) plotGroupZygosity(snpdata$Green[snp,],snpdata$Red[snp,],snpdata$GenCall[snp,],Grouping,NorTum,NormalizedTo,main=snp)  
  dev.off()
}
