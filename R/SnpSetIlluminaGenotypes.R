## Functions to deal with genotypes in SNPList objects
## Initial version 29 may 2006
heterozygousSNPs<-function(object,
  threshold=0.9, 
  useQuality=TRUE, 
  relative=TRUE, 
  percentile=FALSE) {
  ## retrieve a logical from SNPcall, together with quality measure
  ## consider low-quality SNPs to be non-heterozygous
  heterozyg<-assayData(object)[["call"]] == "AB" | assayData(object)[["call"]] == "H"
  if (useQuality) {
    if (!is.null(assayData(object)[["callProbability"]])){
      GenScore<-assayData(object)[["callProbability"]]
      if (relative) {
        if ("GSR" %in% assayDataElementNames(object)) {
          GenScore<-assayData(object)[["GSR"]]
        } else {
	        if ("GTS" %in% varLabels(featureData(object)))
	          GenScore<-sweep(GenScore,1,pData(featureData(object))[,"GTS"],"/")
	        else warning("GenTrainScore is missing; heterozygous SNPs derived from absolute quality score")
        }
      }
      if (percentile) {
        threshold.col<-apply(GenScore,2,quantile,probs=threshold,na.rm=TRUE)
        heterozyg<-heterozyg & sweep(GenScore,2,threshold.col,function(x,y) x>y)

      } else heterozyg<-heterozyg & (GenScore>=threshold)
    } else warning("Quality info missing; heterozygous SNPs derived from genotypes only")
  }
  heterozyg
}

calculateLOH<-function(object,grouping,NorTum="NorTum",...) {
  if (length(NorTum)!=ncol(object)) {
    if (is.null(NorTum) | !(NorTum %in% colnames(pData(object)))) stop("Invalid NorTum argument")
    else NorTum<-pData(object)[,NorTum]
  }
  if (!is.logical(NorTum)) NorTum<-NorTum=="N"
  names(NorTum)<-sampleNames(object)
	#
  hetSNPs<-heterozygousSNPs(object,...)
	loh<-matrix(FALSE,nrow=nrow(object),ncol=ncol(object),dimnames=list(featureNames(object),sampleNames(object)))
	nor.gt<-matrix("",nrow=nrow(object),ncol=ncol(object),dimnames=dimnames(loh))
	nor.qs<-matrix(0,nrow=nrow(object),ncol=ncol(object),dimnames=dimnames(loh))
	for (pageID in levels(factor(grouping))) {
	  samples <- sampleNames(object)[grouping == pageID]
	  n1 <- which(NorTum[samples])
	  t1<-samples[-n1]
	  n1<-samples[n1[1]]
	  if (length(n1)>0 & length(t1)>0) { # at least 1 normal and 1 tumor sample in group
	    for (tum in 1:length(t1)) {
	      loh[,t1[tum]]<-hetSNPs[,n1] & !hetSNPs[,t1[tum]]
	      nor.gt[,t1[tum]]<-hetSNPs[,n1]
	    }
	  }
	}
	assayData(object)[["loh"]]<-loh
	assayData(object)[["nor.gt"]]<-nor.gt
	object
}



  
heterozygosity <-function(genotype,decay=0.8,threshold=0.1) {
  # find stretches of homozygosity
  # depending on decay and treshold low-frequent heterozygous loci are skipped
  if (is.logical(genotype)) heterozygous<-genotype
  else heterozygous<-genotype=="AB" | genotype=="H"
  m<-length(heterozygous)
  heterozygous.dist<-rep(0,m)
  for (i in 1:m) {
    # up
    updist<-0
    found<-0  # decaying number of heterozygous loci in range 
              # a single heterozygous will not interrupt range
              # if at the next heterozygous found is below treshold then
              # that one is skipped also
              # at the default settings (0.8 and 0.1)this can skip about 1 heterozygous 
              # every 10 homozygous
    while(i-updist>1 & (!heterozygous[i-updist] | found<threshold)) {
      if (heterozygous[i-updist]) found<-found+1
      updist<-updist+1
      found<-found*decay
    }
    downdist<-0
    found<-0
    while(i+downdist<m & (!heterozygous[i+downdist] | found<threshold)) {     
      if (heterozygous[i+downdist]) found<-found+1
      downdist<-downdist+1
      found<-found*decay
    }
    heterozygous.dist[i]<-sum(updist,downdist)
  }
  heterozygous.dist
}

compareGenotypes<-function(genotypeT,genotypeN) {
  # Returns character vector 
  # u(ninformative) - both tumor and normal homozygous/not heterozygous
  # i(nformative) - both tumor and normal heterozygous
  # l(oss) - normal heterozygous, tumor homozygous
  # a(rtefact) - normal homozygous, tumor heterozygous
  if (is.logical(genotypeT)) heterozygousT<-genotypeT
  else heterozygousT<-genotypeT=="AB" | genotypeT=="H"
  if (is.logical(genotypeN)) heterozygousN<-genotypeN
  else heterozygousN<-genotypeN=="AB" | genotypeN=="H"
  same<-heterozygousT==heterozygousN
  cmp<-rep("u",length(genotypeT))
  cmp[same & heterozygousN]<-"i"
  cmp[!same & !heterozygousN]<-"a"
  cmp[!same & heterozygousN]<-"l"
  cmp
}

