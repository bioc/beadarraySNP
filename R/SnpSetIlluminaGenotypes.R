## Functions to deal with genotypes in SNPList objects
## Initial version 29 may 2006
heterozygousSNPs<-function(object,
  threshold=0.9, 
  useQuality=TRUE, 
  relative=TRUE, 
  percentile=FALSE) {
  ## retrieve a logical from SNPcall, together with quality measure
  ## consider low-quality SNPs to be non-heterozygous
  heterozyg<-assayData(object)$call == "AB" | assayData(object)$call == "H"
  if (useQuality) {
    if (!is.null(assayData(object)$callProbability)){
      GenScore<-assayData(object)$callProbability
      if (relative) {
        if ("GSR" %in% assayDataElementNames(object)) {
          GenScore<-assayData(object)$GSR
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
  NorTum<-getNorTum(object,NorTum)
	#
	grouping<-getGrouping(object,grouping)
  hetSNPs<-heterozygousSNPs(object,...)
	loh<-matrix(FALSE,nrow=nrow(object),ncol=ncol(object),dimnames=list(featureNames(object),sampleNames(object)))
	nor.gt<-hetSNPs
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
  #
	res<-calculateLair(object,grouping,NorTum)
	assayData(res)$loh<-loh
	assayData(res)$nor.gt<-nor.gt
  res
}

calculateLair<-function(object,grouping=NULL,NorTum="NorTum",min.intensity=NULL,use.homozygous.avg=FALSE) {
  NorTum<-getNorTum(object,NorTum)
	# compute lesser allele intensity ratio value between 0 and 1 (relative to own normal)
	if (is.null(grouping)) {
    snpdata.n<-RG2polar(object) # make sure theta and intensity exist
    theta<-assayData(snpdata.n)$theta
    theta<-1-(theta/(pi/2)) # conform to beadstudio way AA~0, BB~1
    intensity<-assayData(snpdata.n)$intensity[,NorTum]
    gt<-assayData(snpdata.n)$call[,NorTum]
    aa<-ab<-bb<-theta[,NorTum]
    if (is.null(min.intensity)) {
      # this filters extremes 
      min.intensity<-quantile(apply(intensity,1,mean,na.rm=TRUE),probs=0.01)/10
    }
    # Compute average value of AA alleles
    aa[gt!="A" | intensity<min.intensity]<-NA
    aa.avg<-apply(aa,1,mean,na.rm=TRUE)
    # fill in missing values (SNP with no AA normals)
    aa.avgavg<-mean(aa.avg,na.rm=TRUE)
    aa.avg[is.na(aa.avg)]<-aa.avgavg
    #
    bb[gt!="B" | intensity<min.intensity]<-NA
    bb.avg<-apply(bb,1,mean,na.rm=TRUE)
    # fill in missing values (SNP with no BB normals)
    bb.avgavg<-mean(bb.avg,na.rm=TRUE)
    bb.avg[is.na(bb.avg)]<-bb.avgavg
    #
    ab[gt!="H" | intensity<min.intensity]<-NA
    ab.avg<-apply(ab,1,mean,na.rm=TRUE)
    # fill in missing values (SNP with no AB normals) Use middle between values for AA and BB
    ab.avgavg<-(aa.avg+bb.avg)/2
    ab.avg[is.na(ab.avg)]<-ab.avgavg[is.na(ab.avg)]
    # Now compute lair from these values: True heterozygotes should have 1, homozygotes + LOH will have ~0
    if (use.homozygous.avg) {
      # make sure that ab.avg is between aa.avg and bb.avg by zero-ing the homozygous 
      aa.avg[aa.avg-ab.avg > -0.05]<-0
      bb.avg[bb.avg-ab.avg <  0.05]<-1
      # Use average of A and B to define endpoints
      lair<-(theta-aa.avg)/(ab.avg-aa.avg)
      lair2<-(bb.avg-theta)/(bb.avg-ab.avg)
    } else {
      # Use 0 and 1 as endpoints
      lair<-theta/ab.avg
      lair2<-(1-theta)/(1-ab.avg)
    }
    lair<-ifelse(theta-ab.avg<0,lair,lair2)
    # make sure value is between 0 and 1
    lair[lair<0]<-0
    lair[lair>1]<-1
	} else {
	  grouping<-getGrouping(object,grouping)
  	lair<-assayData(object)$G/(assayData(object)$G+assayData(object)$R)
  	lair.n<-lair
    for (smp in unique(grouping)) {
      idx<-which((grouping == smp) & NorTum)[1]
      lair.n[,grouping == smp]<-lair[,idx]
    }
    lair<-ifelse(lair-lair.n<0,(lair/lair.n),((1-lair)/(1-lair.n)))
  }
	res<-object
	assayData(res)$lair<-lair
  res
}

heterozygosity <-function(genotype,decay=0.8,threshold=0.1) {
  # find stretches of homozygosity
  # depending on decay and treshold low-frequent heterozygous loci are skipped
  nas<-is.na(genotype)
  genotype<-genotype[!nas]
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
  res<-rep(NA,length(nas))
  res[!nas]<-heterozygous.dist
  res
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

