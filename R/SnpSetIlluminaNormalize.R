# Normalization procedures for SnpSetIllumina class
# Initial version 12-may-2006
getSubsample<-function(object,subsample) {
  if (length(subsample)!=nrow(object)) {
    if (is.null(subsample) || subsample=="") subsample<-rep(1,nrow(object))
    else subsample<-pData(featureData(object))[,subsample]
  }
  as.factor(subsample)
}

standardNormalization<-function(snpdata) {
  normalizeLoci.SNP(normalizeWithinArrays.SNP(normalizeBetweenAlleles.SNP(snpdata),
       callscore=0.8,relative=TRUE,fixed=FALSE,quantilepersample=TRUE),normalizeTo=2)
}

backgroundEstimate<-function(object,
  method=c("minimum","mode","intmin","anglemode"),
  maxmode=3000,
  bincount=40,
  maxangle=0.3,
  subsample="OPA") {
  ## object   : class SnpSetIllumina; Green and Red not NULL
  ## subsample: character with column name in featureData
  ##            "" or NULL
  ##            vector/factor of length nrow
  ## method   : minimum
  ##            mode   : determine mode by density function
  ##            intmin : determine minimum dependent on other allele
  ##            anglemode : determine mode of angle and set bg accordingly ( needs polar coordinates in snpdata)
  intminBG<-function(allele1,allele2) {
    binsize<-max(allele2,na.rm=TRUE)/bincount
    x<-NULL
    y<-NULL
    ## determine the minima along the intensities of the other allele
    for (bin in 1:bincount) {
      x[bin]<-(bin-0.5) * binsize
      y[bin]<-min(allele1[allele2>(bin -1)*binsize & allele2< bin*binsize])
    }
    y[y == Inf]<-NA
    # fit through minimum as origin
    ym<-y - min(y,na.rm=TRUE)
    fit<-lm(ym ~ x-1)
    
    min(y,na.rm=TRUE) + allele2 * fit$coef

  }

  method<-match.arg(method)
  
  if (method=="anglemode") {
    if ( !all(c("theta","intensity") %in% assayDataElementNames(object))) object<-RG2polar(object)
  }
  subsample<-getSubsample(object,subsample)

  Gb<-matrix(NA,nrow=nrow(object),ncol=ncol(object),dimnames=list(featureNames(object),sampleNames(object)))
  Rb<-matrix(NA,nrow=nrow(object),ncol=ncol(object),dimnames=list(featureNames(object),sampleNames(object)))

  for (OPA in levels(subsample)) {
    probes<-subsample == OPA # select all probes of this OPAset
    for (smp in 1:ncol(object)) {
      if (!all(is.na(assayDataElement(object,"G")[probes,smp]))) {
        switch(method, minimum = {
          min.intensity<-min(assayDataElement(object,"G")[probes,smp],na.rm=TRUE)
          Gb[probes,smp]<-min.intensity-1
          min.intensity<-min(assayDataElement(object,"R")[probes,smp],na.rm=TRUE)
          Rb[probes,smp]<-min.intensity-1
        }, mode = {
          int.density<-density(assayDataElement(object,"G")[probes,smp],na.rm=TRUE)
          dens.probes<-int.density$x<maxmode
          Gb[probes,smp]<-int.density$x[dens.probes][which.max(int.density$y[dens.probes])]
          int.density<-density(assayDataElement(object,"R")[probes,smp],na.rm=TRUE)
          dens.probes<-int.density$x<maxmode
          Rb[probes,smp]<-int.density$x[dens.probes][which.max(int.density$y[dens.probes])]
        }, intmin = {
          Gb[probes,smp]<-intminBG(assayDataElement(object,"G")[probes,smp],assayDataElement(object,"R")[probes,smp])
          Rb[probes,smp]<-intminBG(assayDataElement(object,"R")[probes,smp],assayDataElement(object,"G")[probes,smp])
        }, anglemode = {
          density.p<-density(assayDataElement(object,"theta")[probes,smp],na.rm=TRUE)
          probes.max<-density.p$x<maxangle
          density.p.max<-density.p$x[probes.max]
          mode.p.r<-max(0,density.p.max[which.max(density.p$y[probes.max])],na.rm=TRUE)
          Rb[probes,smp]<-assayDataElement(object,"G")[probes,smp]*sin(mode.p.r)
          probes.max<-density.p$x>(pi/2 -maxangle)
          density.p.max<-density.p$x[probes.max]
          mode.p.g<-min(pi/2,density.p.max[which.max(density.p$y[probes.max])],na.rm=TRUE)
          Gb[probes,smp]<-assayDataElement(object,"R")[probes,smp]*cos(mode.p.g)
        })
      }
    }
  }
  object<-assayDataElementReplace(object,"Gb",Gb)
  object<-assayDataElementReplace(object,"Rb",Rb)
  object
}

backgroundCorrect.SNP<- function(
  object,
  method=c("none", "subtract", "half", "minimum", "edwards", "normexp", "rma"),
  offset = 0){
  # Background correction adapted from limma
  method<-match.arg(method)

  if ( !all(c("Gb","Rb") %in% assayDataElementNames(object))) method<-"none"
  R<-assayDataElement(object,"R")
  G<-assayDataElement(object,"G")
  Rb<-assayDataElement(object,"Rb")
  Gb<-assayDataElement(object,"Gb")

  switch(method, subtract = {
    R<-R-Rb
    G<-G-Gb
  }, half = {
    R<-pmax(R-Rb,0.5)
    G<-pmax(G-Gb,0.5)
  }, minimum = {
    for (slide in 1:ncol(R)) {
      i <- R[, slide] < 1e-18
      if (any(i, na.rm = TRUE)) {
        m <- min(R[!i, slide], na.rm = TRUE)
        R[i, slide] <- m/2
      }
      i <- G[, slide] < 1e-18
      if (any(i, na.rm = TRUE)) {
        m <- min(G[!i, slide], na.rm = TRUE)
        G[i, slide] <- m/2
      }
    }
  }, edwards = {
    one <- matrix(1, nrow(object), 1)
    delta.vec <- function(d, f = 0.1) {
      if (!all(is.na(d))) {quantile(d, mean(d < 1.00000000000000e-16, na.rm = TRUE) *
          (1 + f), na.rm = TRUE)}
      else NA
    }
    sub <- as.matrix(R - Rb)
    delta <- one %*% apply(sub, 2, delta.vec)
    R <- ifelse(sub < delta, delta * exp(1 - (Rb +
        delta)/R), sub)
    sub <- as.matrix(G - Gb)
    delta <- one %*% apply(sub, 2, delta.vec)
    G <- ifelse(sub < delta, delta * exp(1 - (Gb +
        delta)/G), sub)
  }, normexp = {
      for (j in 1:ncol(R)) {
        if(!all(is.na(R[,j]))) {
          x <- G[, j] - Gb[, j]
          out <- normexp.fit(x)
          G[, j] <- normexp.signal(out$par, x)
          x <- R[, j] - Rb[, j]
          out <- normexp.fit(x)
          R[, j] <- normexp.signal(out$par, x)
        }
      }
  }, rma = {
      require("affy")
      R <- apply(R - Rb, 2, bg.adjust)
      G <- apply(G - Gb, 2, bg.adjust)
  })
  if (offset) {
      R <- R + offset
      G <- G + offset
  }
  # remove Rb and Gb from assaydata, and put in new G and R
  #THIS DOESN'T WORK. object<-assayDataElementReplace(object,"Gb",NULL)
  #object<-assayDataElementReplace(object,"Rb",NULL)
  object<-assayDataElementReplace(object,"R",R)
  object<-assayDataElementReplace(object,"G",G)
  object
}
  
normalizeBetweenAlleles.SNP<-function(
  object,
  method=c("quantile"),
  subsample="OPA"){

  method<-match.arg(method)

  subsample<-getSubsample(object,subsample)

  R<-assayDataElement(object,"R")
  G<-assayDataElement(object,"G")
  for (OPA in levels(subsample)) {
    probes<-subsample == OPA # select all probes of this OPAset
    for (smp in 1:ncol(object)) {
      if (!all(is.na(R[probes,smp]))) switch(method,
        quantile = {
          nrm<-normalizeQuantiles(cbind(R[probes,smp],G[probes,smp]))
          R[probes,smp]<-nrm[,1]
          G[probes,smp]<-nrm[,2]
        })
    }
  }
  object<-assayDataElementReplace(object,"R",R)
  object<-assayDataElementReplace(object,"G",G)
  object
}

normalizeWithinArrays.SNP<-function(
  object,
  callscore=0.5,
  normprob=0.5,
  quantilepersample=FALSE,
  relative=FALSE,
  fixed=FALSE,
  useAll=FALSE,
  subsample="OPA") {

  if (useAll) heterozyg<-assayDataElement(object,"call") != ""
  else heterozyg<-assayDataElement(object,"call") == "AB" | assayDataElement(object,"call") == "H"

  # exclude sex chromosomes from normalization
  if ("CHR" %in% varLabels(featureData(object))) {
    numChrom<-numericCHR(pData(featureData(object))[,"CHR"])
    heterozyg[numChrom>90,]<-FALSE
  }
  subsample<-getSubsample(object,subsample)

  GenScore<-assayDataElement(object,"callProbability")
  if (relative) {
    if ("GTS" %in% varLabels(featureData(object))) GenScore<-sweep(GenScore,1,pData(featureData(object))[,"GTS"],"/")
    else stop("featureData does not contain a 'GTS' column, relative=TRUE can not be used")
  }
  if (fixed) gc.min<-callscore
  else if (!quantilepersample) gc.min<-quantile(GenScore[heterozyg],probs=callscore,na.rm=TRUE)

  R<-assayDataElement(object,"R")
  G<-assayDataElement(object,"G")

  for (OPA in levels(subsample)) {
    probes<-subsample == OPA
    green.med<-rep(1,ncol(object))
    red.med<-rep(1,ncol(object))
    for (smp in 1:ncol(object)) {

      if (!all(is.na(G[probes,smp]))) {
        if (quantilepersample & !fixed) gc.min<-quantile(GenScore[heterozyg[,smp] & probes,smp],probs=callscore,na.rm=TRUE)
        green.med[smp]<-quantile(G[heterozyg[,smp] & (GenScore[,smp]>=gc.min) & probes,smp],probs=normprob,na.rm=TRUE)
        red.med[smp]<-quantile(R[heterozyg[,smp] & (GenScore[,smp]>=gc.min) & probes,smp],probs=normprob,na.rm=TRUE)
      }
    }
    G[probes,]<-sweep(G[probes,],2,green.med,FUN="/")
    R[probes,]<-sweep(R[probes,],2,red.med,FUN="/")
  }
  object<-assayDataElementReplace(object,"R",R)
  object<-assayDataElementReplace(object,"G",G)
  object
}

normalizeLoci.SNP <- function(
  object,
  method=c("normals","paired"),
  NorTum="NorTum",
  Gender="Gender",
  Subject="Subject",
  normalizeTo=2) {

  method<-match.arg(method)

  if (length(NorTum)!=ncol(object)) {
    if (is.null(NorTum) | !(NorTum %in% colnames(pData(object)))) NorTum<-rep(TRUE,ncol(object))
    else NorTum<-pData(object)[,NorTum]
  }
  if (!is.logical(NorTum)) NorTum<-NorTum=="N"

  if (length(Gender)!=ncol(object)) {
    if (is.null(Gender) | !(Gender %in% colnames(pData(object)))) Gender<-rep(TRUE,ncol(object))
    else Gender<-pData(object)[,Gender]
  }
  if (!is.logical(Gender)) Gender<-Gender=="F"

  if (length(Subject)!=ncol(object)) {
    if (is.null(Subject) | !(Subject %in% colnames(pData(object)))) Subject<-rep(1,ncol(object))
    else Subject<-pData(object)[,Subject]
  }
  Subject<-as.factor(Subject)

  R<-assayDataElement(object,"R")
  G<-assayDataElement(object,"G")

  switch(method, normals = {
    probe.med<-apply(R[,NorTum,drop=FALSE]+G[,NorTum,drop=FALSE],1,median,na.rm=TRUE)
    if (!all(Gender)) {
      if ("CHR" %in% varLabels(featureData(object))) {
       # Handle Sexchromosomes
        numChrom<-numericCHR(pData(featureData(object))[,"CHR"])
        # Use only normal females for X
        probes<-numChrom==98
        probe.med[probes]<-apply(R[probes,NorTum & Gender,drop=FALSE]+G[probes,NorTum & Gender,drop=FALSE],1,median,na.rm=TRUE)

        # Use only normal males for Y
        probes<-numChrom==99
        probe.med[probes]<-apply(R[probes,NorTum & !Gender,drop=FALSE]+G[probes,NorTum & !Gender,drop=FALSE],1,median,na.rm=TRUE)*2
      }
    }
    R<-sweep(R,1,probe.med,FUN="/")*normalizeTo
    G<-sweep(G,1,probe.med,FUN="/")*normalizeTo

  }, paired = {
    # use Subject to make Tumor/Normal pairs
  })

  object<-assayDataElementReplace(object,"R",R)
  object<-assayDataElementReplace(object,"G",G)
  RG2polar(object)
}


RG2polar <-function(object) {
  R<-assayDataElement(object,"R")
  G<-assayDataElement(object,"G")
  object<-assayDataElementReplace(object,"theta",atan2(G, R))
  object<-assayDataElementReplace(object,"intensity",G+R) # real polar would be r=sqrt(G ^ 2 + R ^ 2)
  object
}

polar2RG <-function(object) {
  theta<-assayDataElement(object,"theta")
  intensity<-assayDataElement(object,"intensity")
  snpdata.cos<-cos(theta)
  Red.perc<-snpdata.cos/(snpdata.cos+sin(theta))

  object<-assayDataElementReplace(object,"R",intensity * Red.perc)
  object<-assayDataElementReplace(object,"G",intensity * (1-Red.perc))
  object
}

