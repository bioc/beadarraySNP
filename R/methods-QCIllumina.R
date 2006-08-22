setMethod("initialize", "QCIllumina",
          function(.Object,
									 arrayType = "Sentrix",
									 arrayID = character(),
									 intensityMed = new("matrix"),
									 greenMed = new("matrix"),
									 redMed = new("matrix"),
									 validn = new("matrix"),
									 annotation = new("matrix"),
									 samples = new("matrix")
									) {
									.Object@arrayID<-arrayID
									.Object@intensityMed<-intensityMed
									.Object@greenMed<-greenMed
									.Object@redMed<-redMed
									.Object@validn<-validn
									.Object@annotation<-annotation
									.Object@samples<-samples
									
									arrayType(.Object) <-arrayType
									
								  .Object
})

setReplaceMethod("arrayType", "QCIllumina", function(object, value) {
	getmatrix<-function(mat) {
		# copy in old values, probably only useful to preserve initialized values
		matn<-matrix(NA,ncol=ncols,nrow=nrows)
		if (!any(dim(mat)==0)){
			matn[1:min(nrow(mat),nrows),1:min(ncol(mat),ncols)]<-mat[1:min(nrow(mat),nrows),1:min(ncol(mat),ncols)]
		}
		matn
	}

	value<-match.arg(value,c("Sentrix"))
	object@arrayType <- value
	if (value=="Sentrix") {
		ncols<-12
		nrows<-8
	}
	
	object@intensityMed<-getmatrix(object@intensityMed)
	object@greenMed<-getmatrix(object@greenMed)
	object@redMed<-getmatrix(object@redMed)
	object@validn<-getmatrix(object@validn)
	object@annotation<-getmatrix(object@annotation)
	object@samples<-getmatrix(object@samples)
  object
})

setMethod("arrayType", "QCIllumina", function(object) object@arrayType)

setReplaceMethod("arrayID", "QCIllumina", function(object, value) {
	object@arrayID<-value
	object
})

setMethod("arrayID", "QCIllumina", function(object) object@arrayID)

setMethod("plotQC", "QCIllumina", function(object,type=c("intensityMed","greenMed","redMed","validn","annotation","samples")) {

  image.plate<-function(z,xdim=dim(z)[2],ydim=dim(z)[1],col = gray (0:99/ 99), zlim=c(0,max(z,na.rm=TRUE)),...) {
    z[is.na(z)]<-0
    image(1:xdim,1:ydim,t(z), zlim=zlim, col = col ,xlab="Col",ylab="Row", ...)
  }

	checkerboard<-function(z,...) {
	  H<-matrix(rep(c(1,0.9),length(z)/2),nrow(z),ncol(z))
    H[,seq(2,ncol(z),by=2)]<-1.9-H[,seq(2,ncol(z),by=2)]
    image.plate(H,zlim=c(0,1),...)
    text(col(z),row(z),labels=z,cex=0.66)
	 }
	 
	 type<-match.arg(type)
	 switch(type, 
	    intensityMed = image.plate(object@intensityMed, main="median Intensity"),
			greenMed = image.plate(object@greenMed, main="median Green", col= rgb(0,0:255,0,max=255)), 
			redMed = image.plate(object@redMed, main="median Red", col= rgb(0:255,0,0,max=255)), 
			validn = image.plate(object@validn,main="valid probes"),
			annotation = checkerboard(object@annotation,main="annotation"),
			samples = checkerboard(object@samples,main="samples"))
	 invisible()
})

setMethod("reportSamplePanelQC", "QCIllumina", function(object, by=10, legend=TRUE, ...) {
  samples<-unique(as.vector(object@samples))
  samples<-samples[!is.na(samples)]
  panels<-unique(as.vector(object@annotation))
  panels<-panels[!is.na(panels)]
  greenmed<-matrix(0,ncol=length(samples),nrow=length(panels),dimnames=list(panels,samples))
  redmed<-matrix(0,ncol=length(samples),nrow=length(panels),dimnames=list(panels,samples))
  for (r in 1:nrow(object@greenMed)) {
    for (c in 1:ncol(object@greenMed)) {
      if (!is.na(object@annotation[r,c])){
        greenmed[object@annotation[r,c],object@samples[r,c]]<-object@greenMed[r,c]
        redmed[object@annotation[r,c],object@samples[r,c]]<-object@redMed[r,c]
      }
    }
  }
  colstart<-3
  colend<-2+length(panels)
  opacol<-c(rgb(r=0,b=0,g=colstart:colend,maxColorValue=colend),rgb(r=colstart:colend,b=0,g=0,maxColorValue=colend))

  greenmed<-rbind(greenmed,redmed)
  for (i in seq(1,length(samples),by=by)) {
    barplot(greenmed[,i:min(length(samples),i+by-1)],beside=TRUE,las=3,col=opacol,...)
  }
  if (legend) {
    plot(c(0,1),type="n",xaxt="n",yaxt="n",xlab="",ylab="",main="",bty="n")
    legend("topleft",legend=c(paste(panels,"Green"),paste(panels,"Red")),fill=opacol,bty="n",cex=1.5)
  }
  invisible()
})

pdfQC<-function(object,filename="arrayQC.pdf",by=10) {
  pdf(filename,paper="a4",width=7.2,height=11)
  par(mfrow=c(3,2),mar=c(2,2,3,1))
  plotQC(object,"greenMed")
  plotQC(object,"redMed")
  plotQC(object,"intensityMed")
  plotQC(object,"validn")
  plotQC(object,"annotation")
  plotQC(object,"samples")
  par(mfrow=c(4,1))
  reportSamplePanelQC(object,by=by)
  dev.off()
}



