## snp reports
## functions to report for all available data in a SnpSetIllumina object
reportSamplesSmoothCopyNumber<-function(snpdata, grouping, normalizedTo=2, smooth.lambda=2, ridge.kappa=0, plotLOH=c("none","marker","line","NorTum"), ...){
  # default grouping is by 4 in sequence of samples in snpdata
  plotLOH<-match.arg(plotLOH)
  if (missing(grouping)) grouping<-floor(seq(along.with=sampleNames(snpdata),by=0.25))
  # make sure chromosmes are sorted
  ind<-order(reporterInfo(snpdata)[,"CHR"],reporterInfo(snpdata)[,"MapInfo"])
  snpdata<-snpdata[ind,]
  sample.colors<-c("red","green","blue","orange","brown","turquoise","yellow","purple","pink","magenta")
  chroms<-unique(reporterInfo(snpdata)[,"CHR"])
  intensities<-assayData(snpdata)[["intensity"]]
  for (pageID in levels(factor(grouping))){
    samples<-sampleNames(snpdata)[grouping == pageID]
    if (length(samples)>0) {
      dchrompos<-prepareGenomePlot(reporterInfo(snpdata)[,c("CHR","MapInfo")],...)
      for (i1 in 1:length(samples) ) {
        for (chrom in chroms) {
          probes<-reporterInfo(snpdata)[,"CHR"] == chrom
          if (sum(!is.na(intensities[probes,samples[i1]]))>10 ) {
            smoothed<-quantsmooth(intensities[probes,samples[i1]],smooth.lambda=smooth.lambda,ridge.kappa=ridge.kappa)
            lines(dchrompos[probes,2],dchrompos[probes,1]+(smoothed-normalizedTo)/normalizedTo,col=sample.colors[i1],lwd=1.5)
            if (plotLOH!="none") {
              probeNames<-featureNames(snpdata)[probes]
              if (plotLOH=="marker") {
                chromhet<-heterozygosity(assayData(snpdata)[["call"]][probes,samples[i1]])
                LOH<-probeNames[chromhet>20]
                if (length(LOH)>0) points(dchrompos[LOH,2],dchrompos[LOH,1]-0.3-(i1*0.05),pch="-",col=sample.colors[i1])
              }
              if (plotLOH=="line") {
                chromhet<-heterozygosity(assayData(snpdata)[["call"]][probes,samples[i1]])
                lines(dchrompos[probes,2],dchrompos[probes,1]+scaleto(chromhet,c(10,40),c(0.1,-0.4)),col=sample.colors[i1],lty=2)
              }
              if (plotLOH=="NorTum" && pData(snpdata)[samples[i1],"NorTum"]!="N") {
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
  pdf(filename,paper="a4",width=7.2,height=11)
	reportSamplesSmoothCopyNumber(object,...)
	dev.off() 
}

reportChromosomesSmoothCopyNumber<-function(snpdata, grouping, normalizedTo=2, smooth.lambda=2, ridge.kappa=0, plotLOH=c("none","marker","line","NorTum"), ...){
  ideo.width<-0.15
  ideo.ypos<-normalizedTo+(ideo.width/2)
  ideo.bleach<-0.25
  plotLOH<-match.arg(plotLOH)
  if (missing(grouping)) grouping<-floor(seq(along.with=sampleNames(snpdata),by=0.25))
  sample.colors<-c("red","green","blue","orange","brown","turquoise","yellow","purple","pink","magenta")
  # make sure chromosmes are sorted
  ind<-order(numericCHR(reporterInfo(snpdata)[,"CHR"]),reporterInfo(snpdata)[,"MapInfo"])
  snpdata<-snpdata[ind,]
  chroms<-unique(numericCHR(reporterInfo(snpdata)[,"CHR"]))
  intensities<-assayData(snpdata)[["intensity"]]
  for (group in levels(factor(grouping))){
    samples<-sampleNames(snpdata)[grouping == group]
    if (length(samples)>0) {
      for (chrom in chroms) {
        probes<-numericCHR(reporterInfo(snpdata)[,"CHR"]) == chrom
        if (any(apply(intensities[probes,samples,drop=FALSE],2,function(x) sum(!is.na(x))>10))) {
          plot(c(0,lengthChromosome(chrom,"bases")),c(1,3),main=paste(group,"Chromosome",characterCHR(chrom)),type="n",ylab="intensity",xlab="",xaxt="n")
          paintCytobands(chrom,pos=c(0,ideo.ypos),units="bases",width=ideo.width,legend=FALSE,bleach=ideo.bleach)
          plotSmoothed(intensities[probes,samples,drop=FALSE],reporterInfo(snpdata)[probes,"MapInfo"],smooth.lambda=smooth.lambda,plotnew=FALSE,...)
          legend("topleft",samples,col=1:length(samples)+1,lty=1,lwd=2,ncol=length(samples))
          if (plotLOH!="none") {
            probeNames<-rownames(snpdata)[probes]
            chromhet<-heterozygosity(assayData(snpdata)[["call"]][probes,samples[i1]])
            if (plotLOH=="marker") {
              LOH<-probeNames[chromhet>20]
              #if (length(LOH)>0) points(dchrompos[LOH,2],dchrompos[LOH,1]-0.3-(i1*0.05),pch="-",col=sample.colors[i1])
            }
            if (plotLOH=="line") {
              #lines(dchrompos[probes,2],dchrompos[probes,1]+scaleto(chromhet,c(10,40),c(0.1,-0.4)),col=sample.colors[i1],lty=2)
            }
            if (plotLOH=="NorTum" && pData(snpdata)[samples[i1],"NorTum"]=="T") {
              ## check availability of normal to compare with
              n1<-grep("N",as.character(pData(snpdata)[samples,"NorTum"]))
              if (length(n1)>0) {
                compGenotype<-compareGenotypes(assayData(snpdata)[["call"]][probes,samples[i1]],assayData(snpdata)[["call"]][probes,samples[n1[1]]])
                LOH<-probeNames[compGenotype=="l"] # loss 
                if(length(LOH)>0) points(reporterInfo(snpdata)[LOH,"MapInfo"],0.4-(i1*0.05),pch="'",col=sample.colors[i1])
                LOH<-probeNames[compGenotype=="i"] # heterozygous normal)
                if(length(LOH)>0) points(reporterInfo(snpdata)[LOH,"MapInfo"],0.4-(i1*0.05),pch="|",col=sample.colors[i1])
              }
            }
          }
        }
      }
    }
  }
}

pdfChromosomesSmoothCopyNumber<-function(object,filename,...) {
  pdf(filename,paper="a4",width=7.2,height=11)
  par(mfrow=c(4,1),mar=c(1.5,2,2,0))
	reportChromosomesSmoothCopyNumber(object,...)
	dev.off() 
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

