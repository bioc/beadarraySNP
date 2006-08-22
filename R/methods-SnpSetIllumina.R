setMethod("initialize", "SnpSetIllumina",
          function(.Object,
                   phenoData = new("AnnotatedDataFrame"),
                   experimentData = new("MIAME"),
                   annotation = character(),
                   call = new("matrix"),
                   callProbability = new("matrix"),
                   G = new("matrix"),
                   R = new("matrix"),
                   reporterInfo = new("data.frameOrNULL"),
                   ... ) {
            .Object<-callNextMethod(.Object,
                           assayData = assayDataNew(
                             call = call,
                             callProbability = callProbability,
                             G = G,
                             R = R,
                             ...),
                           phenoData = phenoData,
                           experimentData = experimentData,
                           annotation = annotation)
            .Object@reporterInfo<-reporterInfo
            if (!is.null(.Object@reporterInfo)) {
              rpNames<-featureNames(.Object)
              if (is.null(rownames(.Object@reporterInfo))) {
                if (dim(.Object)[[1]] != nrow(.Object@reporterInfo)) stop("Conflicting number of rows in 'reporterInfo'")
                rownames(.Object@reporterInfo)<-rpNames
              } else {
                if (!all(rpNames == rownames(.Object@reporterInfo))) stop("Conflicting row names in 'reporterInfo'")
              }

            }
            validObject(.Object)
            .Object

          })

setMethod("[", "SnpSetIllumina", function(x, i, j, ..., drop = FALSE) {
          x<-callNextMethod(x, i, j, ..., drop=drop)
          if(!is.null(reporterInfo(x)) && !missing(i)) reporterInfo(x)<-reporterInfo(x)[i,, ..., drop = drop]
          x
})

setMethod("reporterInfo", "SnpSetIllumina", function(object) object@reporterInfo)

setReplaceMethod("reporterInfo", "SnpSetIllumina", function(object, value) {
  object@reporterInfo <- value
  object
})

setValidity("SnpSetIllumina", function(object) {
  assayDataValidMembers(assayData(object), c("call", "callProbability"))
})

setMethod("exprs", c("SnpSetIllumina"), function(object) assayDataElement(object, "call"))

setReplaceMethod("exprs", c("SnpSetIllumina", "matrix"), function(object, value) {
  assayDataElementReplace(object, "call", value)
})

.mergeAssayData<-function(x, y, newdimnames) {
  # this is derived from assayData combine method
  # differences:
  # - allows different number of reporters/features
  # - will merge data from identical column names into 1 column ie rbind()) 
  # - only works on 2-dimensional assayData elements
  combineElement <- function(x, y) {
    outarr<-array(NA,dim=c(length(newdimnames[[1]]),length(newdimnames[[2]])),newdimnames)
    mode(outarr)<-mode(x)
    outarr[rownames(x),colnames(x)]<-x
    # make sure that values that are NA in y but set in x are not overwritten
    x<-outarr[rownames(y),colnames(y)]
    y[is.na(y)]<-x[is.na(y)]
    #
    outarr[rownames(y),colnames(y)]<-y
    outarr
  }
  storage.mode <- storageMode(x)
  nmfunc <- assayDataElementNames

  if (storageMode(y) != storage.mode)
    stop(paste("assayData must have same storage, but are ",
               storage.mode, ", ", storageMode(y), sep=""))
  if (length(nmfunc(x)) != length(nmfunc(y)))
    stop("assayData have different numbers of elements:\n\t",
         paste(nmfunc(x), collapse=" "), "\n\t",
         paste(nmfunc(y), collapse=" "))
  if (!all(nmfunc(x) == nmfunc(y)))
    stop(paste("assayData have different element names:",
               paste(nmfunc(x), collapse=" "),
               paste(nmfunc(y), collapse=" "), sep="\n\t"))
               
  for (nm in nmfunc(x)) {
    x<-assayDataElementReplace(x,nm,combineElement(assayDataElement(x,nm),assayDataElement(y,nm)))
  }
  x
}

.mergePhenodata<-function(x , y, samples) {
  variables<-union(colnames(pData(x)),colnames(pData(y)))
  outarr<-array(data=NA,dim=c(length(samples),length(variables)),dimnames=list(samples,variables))
  outarr[sampleNames(y),colnames(pData(y))]<-as.matrix(pData(y))
  outarr[sampleNames(x),colnames(pData(x))]<-as.matrix(pData(x))
  pd<-data.frame(outarr)
  vardescs<-union(colnames(varMetadata(x)),colnames(varMetadata(y)))
  outarr<-array(data=NA,dim=c(length(variables),length(vardescs)),dimnames=list(variables,vardescs))
  outarr[colnames(pData(y)),colnames(varMetadata(y))]<-as.matrix(varMetadata(y))
  outarr[colnames(pData(x)),colnames(varMetadata(x))]<-as.matrix(varMetadata(x))
  vd<-data.frame(outarr)
  new("AnnotatedDataFrame", data=pd, varMetadata=vd)
}

.mergeReporterInfo<-function(x, y, reporters) {
  ## merge here will reproduce rows
  if (dim(x)[2] == dim(y)[2] && all(names(x)==names(y)))
     ri <- rbind(x, y)
  else {
     alln <- union(nx <- names(dx <- x), ny <- names(dy <- y))
     if (length(xx <- setdiff(alln,nx))>0)
        for (i in 1:length(xx))
           dx[[ xx[i] ]] <- NA
     if (length(xx <- setdiff(alln,ny))>0)
        for (i in 1:length(xx))
           dy[[ xx[i] ]] <- NA
     ri <- rbind(dx,dy)
  }
  ri[reporters,]
}



setMethod("combine", c("SnpSetIllumina", "SnpSetIllumina"), function(x, y, ...) {

  if (class(x) != class(y))
    stop(paste("objects must be the same class, but are ",
               class(x), ", ", class(y), sep=""))
  ## we need a kind of merge functionality in order to combine OPA panels
  newdimnames<-list(union(featureNames(x),featureNames(y)),union(sampleNames(x),sampleNames(y)))
  x <- .mergeAssayData(x, y, newdimnames)
  # a bit of a hack to only keep the union, and discard double entries
  phenoData(x) <- .mergePhenodata(x, y, newdimnames[[2]])
  experimentData(x) <- combine(experimentData(x),experimentData(y))
  reporterInfo(x)<-.mergeReporterInfo(reporterInfo(x), reporterInfo(y), newdimnames[[1]])
    
  ## annotation -- constant
  if (any(annotation(x) != annotation(y))) {
    warning("objects have different annotations: ",
         annotation(x), ", ", annotation(y))
    annotation(x)<-unique(c(annotation(x),annotation(y)))
  }
  x
})

read.SnpSetIllumina<-function(samplesheet, manifestpath=NULL, reportpath=NULL, rawdatapath=NULL, reportfile=NULL, briefOPAinfo=TRUE, verbose=FALSE) {
  if (verbose) cat("Samplesheet:",ifelse(is.data.frame(samplesheet),"<data.frame>",samplesheet),"\n")
  if (is.data.frame(samplesheet)) {
    samples<-samplesheet
    for (i in 1:length(samples)) samples[[i]]<-as.character(samples[[i]])
    path<-""
  } else {
    path<-dirname(samplesheet)
    firstfield <- scan(samplesheet, what = "", sep = ",", flush = TRUE,
            quiet = TRUE, blank.lines.skip = FALSE, multi.line = FALSE)
    skip <- grep("[Data]", firstfield, fixed=TRUE)
    if (length(skip) == 0) stop("Cannot find \"[Data]\" in samplesheet file")
    samples<-read.table(samplesheet, skip=skip, header = TRUE, sep = ",", as.is = TRUE, check.names = FALSE,colClasses="character")
  }
  samples<-cbind(samples,validn=0,Row=as.numeric(substr(as.character(samples[,"Sentrix_Position"]),3,4)),Col=as.numeric(substr(as.character(samples[,"Sentrix_Position"]),8,9)))
  if (length(unique(samples[,"Sample_Plate"]))>1 | length(unique(samples[,"Pool_ID"]))>1 | length(unique(samples[,"Sentrix_ID"]))>1) stop("Either Sample_Plate or Pool_ID or Sentrix_ID values in samplesheet are not same for all samples")
  OPAname<-as.character(samples[1,"Pool_ID"])
  if (verbose) cat(nrow(samples),"samples in sheet\n")
  #
  if (is.null(manifestpath)) manifestpath<-path
  if (is.null(reportpath)) reportpath<-path
  if (is.null(rawdatapath)) rawdatapath<-path
  # SNPInfo
	SNPinfo<-IlluminaGetOPAinfo(OPAname,manifestpath,brief=briefOPAinfo)
  # PhenoData
  rownames(samples)<-samples[,"Sample_Name"]
  
	if (is.null(reportfile)) { 
    # Import data from GenCall software
    # GenCall, GenScore
    gencalls<-IlluminaGetGencalls(reportpath,OPAname)
    # make SNPinfo and GenCalls contain same ids in case data has less datapoints then annotationfile
    # Order by CHR, MAPinfo
    ind<-order(numericCHR(SNPinfo$CHR),SNPinfo$MapInfo)
    SNPinfo<-SNPinfo[ind,]
    impGenCall<-gencalls$genotypes[ind,]
    impGenScore<-gencalls$callscores[ind,]
    if (verbose) cat(ncol(impGenCall),"samples in report\n")
  	#if (!is.null(gencalls)) locusinfo<-gencalls$"locusinfo"
  	G<-NULL
  	R<-NULL
  	GDev<-NULL
  	RDev<-NULL
  	GenCall<-NULL
  	GenScore<-NULL
    for (sample in 1:nrow(samples)) {
      # read data, sort by rsnumber, new data only has illumnicode
      # drop data that has no rs-codes
      colname<-paste(samples[sample,"Sentrix_ID"],samples[sample,"Sentrix_Position"],sep="_")
      if (verbose) cat(colname)
      GenCall<-cbind(GenCall,impGenCall[,colname])
      GenScore<-cbind(GenScore,impGenScore[,colname])
      sampledata<-read.table(file.path(rawdatapath,paste(paste(samples[sample,"Sentrix_ID"],samples[sample,"Sentrix_Position"],sep="_"),".csv",sep="")),header=TRUE,sep=",",row.names=1,as.is=TRUE)
      sampledata<-sampledata[as.character(SNPinfo[,"IllCode"]),]
      G<-cbind(G,sampledata[,"Mean.GRN"])
      R<-cbind(R,sampledata[,"Mean.RED"])
      GDev<-cbind(GDev,sampledata[,"Dev.GRN"])
      RDev<-cbind(RDev,sampledata[,"Dev.RED"])
      samples[sample,"validn"]<-sum(!is.na(sampledata[,"Mean.GRN"]))
    }
    if (verbose) cat("\n")
    # set all names
    colnames(G)<-samples[,"Sample_Name"]
    rownames(G)<-SNPinfo[,"snpid"]
    dimnames(R)<-dimnames(G)
    dimnames(GDev)<-dimnames(G)
    dimnames(RDev)<-dimnames(G)
    dimnames(GenCall)<-dimnames(G)
    dimnames(GenScore)<-dimnames(G)
    rownames(SNPinfo)<-SNPinfo[,"snpid"]
    SNPinfo<-cbind(SNPinfo,GTS=as.numeric(gencalls$locusinfo[rownames(SNPinfo),"GTS"]))
  } else {
    # Import data from BeadStudio report file
    firstfield <- scan(reportfile, what = "", sep = ",", flush = TRUE,
            quiet = TRUE, blank.lines.skip = FALSE, multi.line = FALSE,nlines=100)
    skip <- grep("[Data]", firstfield, fixed=TRUE)
    if (length(skip) == 0) stop("Cannot find \"[Data]\" in report file")
    alldata<-read.table(reportfile, skip=skip, header = TRUE, sep = "\t", as.is = TRUE, check.names = FALSE,colClasses="character")
    # Integrity checks
    essentialcols<-c("SNP Name","Sample ID","GC Score","Allele1 - AB","Allele2 - AB","GT Score","Cluster Sep","Theta","R","X Raw","Y Raw")
    foundcols<-essentialcols %in% colnames(alldata)
    if (!all(foundcols)) stop ("Columns:",paste(essentialcols[!foundcols])," are missing in the report file")
    m<-length(unique(alldata[,"SNP Name"]))
    datasamples<-unique(alldata[,"Sample ID"])
    n<-length(datasamples)
    if (nrow(alldata) != m*n) stop("Datarows are missing in report file")
    ind<-order(alldata[,"Sample ID"],alldata[,"SNP Name"])
    alldata<-alldata[ind,]
    # Determine col and row names
    tmp1<-matrix(alldata[,"SNP Name"],nrow=m,ncol=n,byrow=FALSE)
    tmp2<-matrix(alldata[,"Sample ID"],nrow=m,ncol=n,byrow=FALSE)
    newdimnames<-list(tmp1[,1],tmp2[1,])
    # Extract data
    G<-matrix(as.numeric(alldata[,"X Raw"]),nrow=m,ncol=n,byrow=FALSE,dimnames=newdimnames)
    R<-matrix(as.numeric(alldata[,"Y Raw"]),nrow=m,ncol=n,byrow=FALSE,dimnames=newdimnames)
    GDev<-NULL
    RDev<-NULL
    GenScore<-matrix(as.numeric(alldata[,"GC Score"]),nrow=m,ncol=n,byrow=FALSE,dimnames=newdimnames)
    GTS<-matrix(as.numeric(alldata[,"GT Score"]),nrow=m,ncol=n,byrow=FALSE,dimnames=newdimnames)
    GTS<-apply(GTS,1,max,na.rm=TRUE)
    GenCall<-paste(alldata[,"Allele1 - AB"],alldata[,"Allele2 - AB"],sep="")
    GenCall<-matrix(sub("AA","A",sub("BB","B",sub("AB","H",GenCall))),nrow=m,ncol=n,byrow=FALSE,dimnames=newdimnames)
    # Select only data from samplesheet
    selected<-paste(samples[,"Sentrix_ID"],samples[,"Sentrix_Position"],sep="_")
    G<-G[,selected]
    R<-R[,selected]
    GenScore<-GenScore[,selected]
    GenCall<-GenCall[,selected]
    colnames(G)<-samples[,"Sample_Name"]
    colnames(R)<-samples[,"Sample_Name"]
    colnames(GenScore)<-samples[,"Sample_Name"]
    colnames(GenCall)<-samples[,"Sample_Name"]
    SNPinfo<-SNPinfo[rownames(G),]
    SNPinfo<-cbind(SNPinfo,GTS=GTS)
    samples[,"validn"]<-apply(G,2,function(x) sum(!is.na(x)))
  }
  new("SnpSetIllumina",phenoData=new("AnnotatedDataFrame",samples,data.frame(labelDescription=colnames(samples),row.names=colnames(samples))), annotation=OPAname, call=GenCall, callProbability=GenScore, G=G, R=R,
             reporterInfo=SNPinfo,storage.mode="list")
}

getExperiments <- function(file="experiments.txt",path=NULL)
{
  if (!is.null(path))
    file <- file.path(path, file)
  tab <- scan(file, what="character")
  if (!is.null(path))
    tab<-file.path(path,tab)
  as.vector(tab)
}

IlluminaGetOPAinfo<-function(OPAname,OPAinfoPath,brief=TRUE) {
	# find .opa file
	# Will be used for SNPinfo field in SNPlist
	# Some columns with special meaning
  # CHR : chromosome
  # MapInfo : position on chrom
  # OPA : linkage panel
  # snpid : international snp id rsxxxxxx (used as row names)
  # IllCode : numeric within linkage panel to connect to snpid
	OPAfile<-list.files(OPAinfoPath,pattern=paste(OPAname,".*\.opa$",sep=""),full.names=TRUE)
	if (length(OPAfile) != 1) stop(paste("OPA info file could not be (uniquely) identified for",OPAset))
	# import it to the database
  firstfield <- scan(OPAfile, what = "", sep = ",", flush = TRUE, quiet = TRUE, blank.lines.skip = FALSE, multi.line = FALSE)
  skip <- grep("Ilmn ID", firstfield, fixed=TRUE)
  if (length(skip) == 0) stop("Cannot find \"Ilmn ID\" in OPA info file")
	enddata<- grep("[Gentrain Request]", firstfield, fixed=TRUE)
  if (length(enddata) == 0) stop("Cannot find \"[Gentrain Request]\" in OPA info file")
	#OPAmetaInfo<-read.table(OPAfile,sep=",",nrows=skip[1]-1,fill=TRUE,as.is=TRUE)
	#OPAtestName<-OPAmetaInfo[OPAmetaInfo[,1]=="Test Name",2]
	#OPAversion<-OPAmetaInfo[OPAmetaInfo[,1]=="Test Version",2]
	#OPAdate<-OPAmetaInfo[OPAmetaInfo[,1]=="Date Manufactured",2]
	OPAinfo<-read.table(OPAfile, skip=skip[1]-1, header = TRUE, sep = ",", as.is = TRUE, check.names = FALSE, nrows=enddata[1]-skip[1]-1)
  colnames(OPAinfo)<-c("Illname","snpid","oligo1","oligo2","oligo3","IllCode","IllOligo","IllStrand","snpbases","CHR","Ploidy","Species","MapInfo","TopGenomicSeq","CustomerStrand")
  rownames(OPAinfo)<-OPAinfo[,"snpid"]
  OPAnames<-rep(OPAname,nrow(OPAinfo))
  if (brief) cbind(OPA=I(OPAnames),OPAinfo[,c("snpid","IllCode","CHR","MapInfo")])
  else cbind(OPA=I(OPAnames),OPAinfo)
}


IlluminaGetGencalls<- function(path,OPAname) {
  gencallfile<-list.files(path,pattern=paste(OPAname,".*LocusByDNA.*csv",sep=""),full.names=TRUE)
  # exclude some common report types *DNA_Report.csv/*Locus_Report.csv/*Final.csv
  exclude<-grep("DNA_Report\.csv|Locus_Report\.csv|Final\.csv",gencallfile)
  if (length(exclude)>0) gencallfile<-gencallfile[-exclude]
  if (length(gencallfile)==1) {
    firstfield <- scan(gencallfile, what = "", sep = ",", flush = TRUE,quiet = TRUE, blank.lines.skip = FALSE, multi.line = FALSE)
    skip <- grep("instituteLabel", firstfield, fixed=TRUE)
    gencalls<-t(read.table(gencallfile, skip=skip[1], header = FALSE, sep = ",", as.is = TRUE, check.names = FALSE))
    skip <- grep("oligoPoolId", firstfield, fixed=TRUE)
    locusinfo<-t(read.table(gencallfile, skip=skip[1], header = FALSE, row.names=2, sep = ",", as.is = TRUE, check.names = FALSE, nrows=7))[-2:-1,c(1,5,6)]
    colnames(locusinfo)<-c("GTS","snpid","IllCode")
    alternate<- rep(c(TRUE,FALSE),dim(gencalls)[2] / 2)
    genotypes<-gencalls[,alternate]
    callscores<-gencalls[,!alternate]
    colnames(genotypes)<-genotypes[5,]    # set to sentrixid_col_row
    colnames(callscores)<-callscores[5,]
    genotypes<-genotypes[-8:-1,]
    callscores<-callscores[-8:-1,]
    callscore.dimnames<-dimnames(callscores)
    callscores<-matrix(as.numeric(callscores),nrow=nrow(genotypes),dimnames=callscore.dimnames)
    rownames(locusinfo)<-locusinfo[,"snpid"]
    rownames(callscores)<-locusinfo[,"snpid"]
    rownames(genotypes)<-locusinfo[,"snpid"]
    list("genotypes"=genotypes,"callscores"=callscores,"locusinfo"=locusinfo)
  } else NULL
}
