# QC
removeLowQualityProbes<-function(object, cutoff=0.25) {
  # disable probes that have low overall intensity, ie fraction below overal median
  if ( !("intensity" %in% assayDataElementNames(object))) object<-RG2polar(object)
  probes<-numericCHR(fData(object)$CHR)>90
  int.med<-median(assayData(object)[["intensity"]][probes,],na.rm=TRUE)
  int.probe.med<-apply(assayData(object)[["intensity"]],1,median,na.rm=TRUE)
  probes<-int.probe.med>(int.med*cutoff) | probes
  cat(sum(!probes),"probes removed with median intensity below",int.med*cutoff,"\n")
  object[probes,]
}

removeLowQualitySamples<-function(object, min.intensity=1500, min.gt=100, subsample="OPA"){
  # remove samples that show median intensity below min.intensity in either red or green channel
  # remove samples that have less than min.gt genotyped probes
  subsample<-getSubsample(object,subsample)

  R<-assayDataElement(object,"R")
  G<-assayDataElement(object,"G")
  my_call<-assayDataElement(object,"call")
  callProbability<-assayDataElement(object,"callProbability")
  for (OPA in levels(subsample)) {
    probes<-subsample == OPA # select all probes of this OPAset
    greenmed<-apply(G[probes,],2,median,na.rm=TRUE)
    redmed<-apply(R[probes,],2,median,na.rm=TRUE)
    untyped<-apply(my_call[probes,],2,function(x) sum(x != "U",na.rm=TRUE))
    low<- greenmed<min.intensity | redmed<min.intensity | untyped<min.gt
    cat(OPA,sum(!is.na(low)),"subsamples, ")
    low[is.na(low)]<-FALSE
    cat(sum(low),"removed\n")
    G[probes,low]<-NA
    R[probes,low]<-NA
    my_call[probes,low]<-NA
    callProbability[probes,low]<-NA
  }
  object<-assayDataElementReplace(object,"R",R)
  object<-assayDataElementReplace(object,"G",G)
  object<-assayDataElementReplace(object,"call",my_call)
  object<-assayDataElementReplace(object,"callProbability",callProbability)
  allmissing<-apply(G,2,function(x) all(is.na(x)))
  if (any(allmissing)) {
    warning(paste("Sample(s)",toString(sampleNames(object)[allmissing]),"removed from result because no data remained"))
    object<-object[,!allmissing]
  }
  object
}


calculateQCarray<-function(object,QCobject=NULL,arrayType="Sentrix96") {
  # object should be SnpSetIllumina
  if (class(object)!="SnpSetIllumina") stop("object not usable for this function")
  # test for beadstudio combined samplesheet
  barcodesCol<-grep("^SentrixBarcode_",colnames(pData(object)))
  if ((length(barcodesCol) > 0) | is.list(QCobject)) {
    # make a (list of) qcobject for all opapanels in samplesheet
    if (is.null(QCobject)) qcobjects<-list()
    else qcobjects<-QCobject
    opapanels<-sub("^SentrixBarcode_","",colnames(pData(object))[barcodesCol])
    opaAnnot<-annotation(object)[match(opapanels,LETTERS)]
    for (opa in 1:length(barcodesCol)) {
      probes<-pData(featureData(object))$OPA==opaAnnot[opa]
      barcodes<-unique(pData(object)[,barcodesCol[opa]])
      barcodes<-barcodes[barcodes!=""]
      for (barcode in barcodes) {
        samples<-pData(object)[,barcodesCol[opa]]==barcode
        objectsub<-object[probes,samples]
        smpPosition<-pData(objectsub)[,paste("SentrixPosition_",opapanels[opa],sep="")]
        pData(objectsub)<-cbind(pData(objectsub)[,c("Sample_Name","Sample_ID")],Sentrix_ID=pData(objectsub)[,barcodesCol[opa]],Col=as.numeric(substr(smpPosition,7,9)),Row=as.numeric(substr(smpPosition,2,4)))
        annotation(objectsub)<-opaAnnot[opa]
        qcobjects[[barcode]]<-calculateQCarray(objectsub,qcobjects[[barcode]],arrayType)
      }
    }
    if (length(qcobjects)==1) return(qcobjects[[1]])
    else return(qcobjects)
  } else {
    # object should not be combined 
    if (length(annotation(object))>1) stop("function does not work with combined datasets")
    if (is.null(QCobject)) {
       QCobject<-new("QCIllumina")
       arrayID(QCobject)<-as.character(pData(object)$Sentrix_ID[1])
       arrayType(QCobject)<-arrayType
    }
    R<-assayDataElement(object,"R")
    G<-assayDataElement(object,"G")
    conf<-assayDataElement(object,"callProbability")
    calls<-assayDataElement(object,"call")
    int<-R+G
    idx<-order(numericCHR(pData(featureData(object))$CHR),pData(featureData(object))$MapInfo)
    int<-int[idx,,drop=FALSE]
    ptpdiff<-abs(int[-1,,drop=FALSE]-int[-nrow(int),,drop=FALSE])/(int[-1,,drop=FALSE]+int[-nrow(int),,drop=FALSE])
    
    for (smp in sampleNames(object)) {
      if (arrayType(QCobject) %in% c("Sentrix96","Sentrix16")) {
        co<-pData(object)[smp,"Col"]
        ro<-pData(object)[smp,"Row"]
    	} else if (arrayType(QCobject)=="Slide12"){
        co<-match(pData(object)[smp,"Sentrix_Position"],LETTERS)
        ro<-1
    	}
      if (pData(object)[smp,"Sentrix_ID"]!=arrayID(QCobject)) stop("data in QCobject cannot be from different arrays")
  		QCobject@samples[ro,co]<-smp
  		QCobject@annotation[ro,co]<-annotation(object)
  		QCobject@validn[ro,co]<-sum(!is.na(int[,smp]))
  		QCobject@intensityMed[ro,co]<-median(int[,smp],na.rm=TRUE)
  		QCobject@greenMed[ro,co]<-median(G[,smp],na.rm=TRUE)
  		QCobject@redMed[ro,co]<-median(R[,smp],na.rm=TRUE)
  		QCobject@ptpdiff[ro,co]<-median(ptpdiff[,smp],na.rm=TRUE)
  		probes<-conf[,smp]>0.25
  		probes[is.na(probes)]<-FALSE
  	  QCobject@callrate[ro,co]<-sum(probes)/length(probes)
  	  if (sum(probes)>0)  QCobject@hetPerc[ro,co]<-sum(calls[probes,smp]=="H")/sum(probes)
  	  else QCobject@hetPerc[ro,co]<-0
    }
    QCobject
  }
}

BeadstudioQC<-function(object,QClist=list(),arrayType="Sentrix96") {
  # reconstruct chips from beadstudio samplesheet 
  # create a list with QCIllumina objects for each chip
  # split in OPA panels and chips
  OPAs<-annotation(object)
  for (OPA in 1:length(OPAs)) {
    snps<-pData(featureData(object))$OPA == OPAs[OPA]
    SentrixBarcode<-pData(object)[,paste("SentrixBarcode",LETTERS[OPA],sep="_")]
    chips<-unique(SentrixBarcode)
    for (chip in 1:length(chips)) {
      if(chips[chip] %in% names(QClist)) qcObj<-QClist[[chips[chip]]] else qcObj<-NULL
      curObj<-object[snps,chips == chips[chip]]
      SentrixBarcode<-pData(curObj)[,paste("SentrixBarcode",LETTERS[OPA],sep="_")]
      SentrixPosition<-pData(curObj)[,paste("SentrixPosition",LETTERS[OPA],sep="_")]
      annotation(curObj)<-OPAs[OPA]
      if (nchar(SentrixPosition[1])==9) {# Sentrix arrays on 96 well sample plate
        Row<-as.numeric(substr(as.character(SentrixPosition),3,4))
        Col<-as.numeric(substr(as.character(SentrixPosition),8,9))
        phenoData(curObj)<-new("AnnotatedDataFrame",data=cbind(pData(curObj),Row=Row,Col=Col,Sentrix_Position=SentrixPosition,Sentrix_ID=SentrixBarcode))
      } else {
        phenoData(curObj)<-new("AnnotatedDataFrame",data=cbind(pData(curObj),Sentrix_Position=SentrixPosition,Sentrix_ID=SentrixBarcode))
      }
      
      QClist[[chips[chip]]]<-calculateQCarray(curObj,qcObj,arrayType)
    }
  }
  return(QClist)
}

pdfBeadstudioQC<-function(QClist,basename="beadstudio",by=10) {
  for (qc in names(QClist)) 
    pdfQC(QClist[[qc]],filename=paste(basename,"_",qc,".pdf",sep=""),by=by)

}

dist.GT<-function(object) {
	n <- ncol(object)
	m <- nrow(object)
	H <- matrix(rep(0,n^2),n,n)
	for (i in 1:n) {
		for (j in 1:i) {
			if (i==j) H[i,i] <- 1
			else {
			  probes<-!is.na(exprs(object)[,i]) & !is.na(exprs(object)[,j])
				H[i,j] <- 1- sum(exprs(object)[probes,i]==exprs(object)[probes,j])/sum(probes)
				H[j,i] <- H[i,j]
			}
		}
  }
  colnames(H)<-rownames(H)<-sampleNames(object)
  return(as.dist(H))
}