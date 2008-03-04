setMethod("initialize", "SnpSetIllumina",
          function(.Object,
                   phenoData = new("AnnotatedDataFrame"),
                   experimentData = new("MIAME"),
                   annotation = character(),
                   call = new("matrix"),
                   callProbability = new("matrix"),
                   G = new("matrix"),
                   R = new("matrix"),
                   featureData = new("AnnotatedDataFrame"),
                   extraData = NULL,
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
                           annotation = annotation,
													 featureData = featureData)
            if (!is.null(extraData)) {
              for (m in names(extraData))
                .Object<-assayDataElementReplace(.Object, m, extraData[[m]])
            }
            validObject(.Object)
            .Object

          })

setValidity("SnpSetIllumina", function(object) {
  assayDataValidMembers(assayData(object), c("call", "callProbability"))
})

setMethod("exprs", c("SnpSetIllumina"), function(object) assayDataElement(object, "call"))

setReplaceMethod("exprs", c("SnpSetIllumina", "matrix"), function(object, value) {
  assayDataElementReplace(object, "call", value)
})

setMethod("fData", c("SnpSetIllumina"), function(object) pData(featureData(object)))

setReplaceMethod("fData", c("SnpSetIllumina", "data.frame"), function(object, value) {
  pData(featureData(object))<- value
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
  storMode <- storageMode(x)
  nmfunc <- function(es) if (storageMode(es) == "list") names(assayData(es)) else ls(assayData(es))

  if (storageMode(y) != storMode)
    stop(paste("assayData must have same storage, but are ",
               storMode, ", ", storageMode(y), sep=""))
  if (length(nmfunc(x)) != length(nmfunc(y)))
    stop("assayData have different numbers of elements:\n\t",
         paste(nmfunc(x), collapse=" "), "\n\t",
         paste(nmfunc(y), collapse=" "))
  if (!all(nmfunc(x) == nmfunc(y)))
    stop(paste("assayData have different element names:",
               paste(nmfunc(x), collapse=" "),
               paste(nmfunc(y), collapse=" "), sep="\n\t"))
  res<-x
  for (nm in nmfunc(res)) {
    res<-assayDataElementReplace(res,nm,combineElement(assayDataElement(res,nm),assayDataElement(y,nm)))
  }
  res
}

.mergeAnnotateddata<-function(x , y, samples) {
  allvariables<-union(colnames(pData(x)),colnames(pData(y)))
  commonvariables<-intersect(colnames(pData(x)),colnames(pData(y)))
  yexclusiverows<-setdiff(rownames(pData(y)),rownames(pData(x)))
  # start with core data.frame, and build up
  outarr<-rbind(pData(x)[,commonvariables],pData(y)[yexclusiverows,commonvariables])
  xonly<-setdiff(colnames(pData(x)),commonvariables)
  if (length(xonly)>0) {
    outarr<-cbind(outarr,rbind(pData(x)[,xonly],array(NA,dim=c(length(yexclusiverows),length(xonly)),
            dimnames=list(yexclusiverows,xonly))))
  }
  yonly<-setdiff(colnames(pData(y)),commonvariables)
  if (length(yonly)>0) {
    outarr<-cbind(outarr,rbind(array(NA,dim=c(dim(x)[1],length(yonly)),dimnames=list(featureNames(x),yonly)
       ),pData(y)[yexclusiverows,yonly]))
  }
  pd<-outarr[samples,]
  vardescs<-union(colnames(varMetadata(x)),colnames(varMetadata(y)))
  outarr<-array(data=NA,dim=c(length(allvariables),length(vardescs)),dimnames=list(allvariables,vardescs))
  outarr[colnames(pData(y)),colnames(varMetadata(y))]<-as.matrix(varMetadata(y))
  outarr[colnames(pData(x)),colnames(varMetadata(x))]<-as.matrix(varMetadata(x))
  vd<-data.frame(outarr[colnames(pd),])
  new("AnnotatedDataFrame", data=pd, varMetadata=vd)
}

.mergefeatureData<-function(x, y, features) {
	## first do data
  ## merge here will reproduce rows
  xd<-pData(x)
  yd<-pData(y)
  if (dim(xd)[2] == dim(yd)[2] && all(names(xd)==names(yd)))
     ri <- rbind(xd, yd)
  else {
     alln <- union(nx <- names(dx <- xd), ny <- names(dy <- yd))
     if (length(xx <- setdiff(alln,nx))>0)
        for (i in 1:length(xx))
           dx[[ xx[i] ]] <- NA
     if (length(xx <- setdiff(alln,ny))>0)
        for (i in 1:length(xx))
           dy[[ xx[i] ]] <- NA
     ri <- rbind(dx,dy)
  }
  ri<-ri[features,]
  # varMetadata
  ld<-c(as.character(varMetadata(x)$labelDescription),as.character(varMetadata(y)$labelDescription))
  rn<-c(rownames(varMetadata(x)),rownames(varMetadata(y)))
}



setMethod("combine", c("SnpSetIllumina", "SnpSetIllumina"), function(x, y, ...) {

  if (class(x) != class(y))
    stop(paste("objects must be the same class, but are ",
               class(x), ", ", class(y), sep=""))
  ## we need a kind of merge functionality in order to combine OPA panels
  newdimnames<-list(union(featureNames(x),featureNames(y)),union(sampleNames(x),sampleNames(y)))
  x <- .mergeAssayData(x, y, newdimnames)
  # a bit of a hack to only keep the union, and discard double entries
  phenoData(x) <- .mergeAnnotateddata(phenoData(x), phenoData(y), newdimnames[[2]])
  experimentData(x) <- combine(experimentData(x),experimentData(y))
  featureData(x)<-.mergeAnnotateddata(featureData(x), featureData(y), newdimnames[[1]])
    
  ## annotation -- constant
  if (any(annotation(x) != annotation(y))) {
    warning("objects have different annotations: ",
         annotation(x), ", ", annotation(y))
    annotation(x)<-unique(c(annotation(x),annotation(y)))
  }
  x
})

setMethod("calculateGSR", "SnpSetIllumina", function(object) {
  assayDataElementReplace(object,"GSR",sweep(assayData(object)[["callProbability"]],1,pData(featureData(object))[,"GTS"],"/"))
})

read.SnpSetIllumina<-function(samplesheet, manifestpath=NULL, reportpath=NULL, rawdatapath=NULL,
  reportfile=NULL, briefOPAinfo=TRUE, verbose=FALSE, readTIF=FALSE, ...) {
  if (verbose) cat("Samplesheet:",ifelse(is.data.frame(samplesheet),"<data.frame>",samplesheet),"\n")
  if (is.data.frame(samplesheet)) {
    samples<-samplesheet
    for (i in 1:length(samples)) samples[[i]]<-as.character(samples[[i]])
    path<-"."
    beadstudio<-FALSE
  } else {
    path<-dirname(samplesheet)
    firstfield <- scan(samplesheet, what = "", sep = ",", flush = TRUE,
            quiet = TRUE, blank.lines.skip = FALSE, multi.line = FALSE)
    # test if beadstudio sample sheet
    manif <- grep("[Manifests]", firstfield, fixed=TRUE)
    beadstudio<-length(manif)>0
    skip <- grep("[Data]", firstfield, fixed=TRUE)
    if (length(skip) == 0) stop("Cannot find \"[Data]\" in samplesheet file")
    if (beadstudio) manifests<-read.table(samplesheet, skip=manif, header = FALSE, sep = ",", as.is = TRUE, check.names = FALSE,colClasses="character",row.names=1,nrows=skip-manif-1)

    samples<-read.table(samplesheet, skip=skip, header = TRUE, sep = ",", as.is = TRUE, check.names = FALSE,colClasses="character")
  }
  if (beadstudio && !is.null(reportfile)) {
    # this only works for beadstudio final report files
    rownames(samples)<-samples[,"Sample_ID"]
    req_cols<-c(paste("SentrixBarcode",rownames(manifests),sep="_"),paste("SentrixPosition",rownames(manifests),sep="_"))
    smpcol<-!(req_cols %in% colnames(samples))
    if (any(smpcol))
      stop("Sample sheet error, column(s) ",paste(req_cols[smpcol],collapse=", ")," are not available")
    OPAname<-as.character(manifests[,1])

  } else {
    req_cols<-c("Sample_Name","Sentrix_Position","Sample_Plate","Pool_ID","Sentrix_ID" )
    smpcol<-!(req_cols %in% colnames(samples))
    if (any(smpcol))
      stop("Sample sheet error, column(s) ",paste(req_cols[smpcol],collapse=", ")," are not available")
    if (nchar(samples[1,"Sentrix_Position"])==9) # Sentrix arrays on 96 well sample plate
      samples<-cbind(samples,validn=0,Row=as.numeric(substr(as.character(samples[,"Sentrix_Position"]),3,4)),Col=as.numeric(substr(as.character(samples[,"Sentrix_Position"]),8,9)))
    if (length(unique(samples[,"Sample_Plate"]))>1 | length(unique(samples[,"Pool_ID"]))>1 | length(unique(samples[,"Sentrix_ID"]))>1) stop("Either Sample_Plate or Pool_ID or Sentrix_ID values in samplesheet are not same for all samples")
    OPAname<-as.character(samples[1,"Pool_ID"])
    rownames(samples)<-samples[,"Sample_Name"]
  }
  if (verbose) cat(nrow(samples),"samples in sheet\n")
  #
  if (is.null(manifestpath)) manifestpath<-path
  if (is.null(reportpath)) reportpath<-path
  if (is.null(rawdatapath)) rawdatapath<-path
  # SNPInfo
	SNPinfo<-IlluminaGetOPAinfo(OPAname,manifestpath,brief=briefOPAinfo)

	if (is.null(reportfile)) {
    if(is.null(SNPinfo)) stop(paste("OPA info file could not be (uniquely) identified for",OPAname))
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
  	BeadData<-list()
    if (readTIF) {
      if (!require(beadarray)) stop("beadarray package is needed to load data from tif files")
      BL<-readIllumina(paste(samples[,"Sentrix_ID"],samples[,"Sentrix_Position"],sep="_"),textType=".txt",rawdatapath,...)
      BS<-createBeadSummaryData(BL,imagesPerArray=1)
      ind<-as.character(SNPinfo$IllCode)
      G<-assayData(BS)$G[ind,]
      R<-assayData(BS)$R[ind,]
      GDev<-assayData(BS)$GBeadStDev[ind,]
      RDev<-assayData(BS)$RBeadStDev[ind,]
      samples[,"validn"]<-apply(G,2,function(x) sum(!is.na(x)))
      rm(BS,BL)
    }
    for (sample in 1:nrow(samples)) {
      # read data, sort by rsnumber, new data only has illumnicode
      # drop data that has no rs-codes
      colname<-paste(samples[sample,"Sentrix_ID"],samples[sample,"Sentrix_Position"],sep="_")
      if (verbose) cat(colname)
      if (colname %in% colnames(impGenCall)) {
        GenCall<-cbind(GenCall,impGenCall[,colname])
        GenScore<-cbind(GenScore,impGenScore[,colname])
        if (!readTIF) {
          beadfile<-list.files(rawdatapath,pattern=paste(paste(samples[sample,"Sentrix_ID"],samples[sample,"Sentrix_Position"],sep="_"),".txt",sep=""),full.names=TRUE)
          if (length(beadfile)==1) {
            BeadData[[colname]]<-read.table(beadfile,header=TRUE,sep="\t",as.is=TRUE)
          }
          summaryfile<-list.files(rawdatapath,pattern=paste(paste(samples[sample,"Sentrix_ID"],samples[sample,"Sentrix_Position"],sep="_"),".csv",sep=""),full.names=TRUE)
          if (length(summaryfile)==1) {
            sampledata<-read.table(summaryfile,header=TRUE,sep=",",row.names=1,as.is=TRUE)
            sampledata<-sampledata[as.character(SNPinfo[,"IllCode"]),]
          } else {
            # calculate from beaddata
            if (length(beadfile)!=1) stop(paste("Missing files to determine intensity for sample",colname))
            sampledata<-aggregate(BeadData[[colname]][,c("Grn","Red")],by=list(BeadData[[colname]]$Code),FUN=median)
            rownames(sampledata)<-sampledata[,1]
            sampledev<-aggregate(BeadData[[colname]][,c("Grn","Red")],by=list(BeadData[[colname]]$Code),FUN=sd)
            sampledata<-cbind(sampledata[,2:3],sampledev[,2:3])
            sampledata<-sampledata[as.character(SNPinfo[,"IllCode"]),]
            colnames(sampledata)<-c("Mean.GRN","Mean.RED","Dev.GRN","Dev.RED")
          }
          G<-cbind(G,sampledata[,"Mean.GRN"])
          R<-cbind(R,sampledata[,"Mean.RED"])
          GDev<-cbind(GDev,sampledata[,"Dev.GRN"])
          RDev<-cbind(RDev,sampledata[,"Dev.RED"])
          samples[sample,"validn"]<-sum(!is.na(sampledata[,"Mean.GRN"]))
        }
      } else {
        warning(paste("Sample",rownames(samples)[sample],"is defined in samplesheet, but is not in the reportfile"))
      }
    }
    if (verbose) cat("\n")
    # set all names
    colnames(G)<-rownames(samples)
    rownames(G)<-SNPinfo[,"snpid"]
    dimnames(R)<-dimnames(G)
    dimnames(GDev)<-dimnames(G)
    dimnames(RDev)<-dimnames(G)
    dimnames(GenCall)<-dimnames(G)
    dimnames(GenScore)<-dimnames(G)
    rownames(SNPinfo)<-SNPinfo[,"snpid"]
    SNPinfo<-cbind(SNPinfo,GTS=as.numeric(gencalls$locusinfo[rownames(SNPinfo),"GTS"]))
    extraData<-NULL
  } else {
    chrompos.fromReport<-FALSE
    if(is.null(SNPinfo)) {
      warning(paste("OPA info file could not be (uniquely) identified for",OPAname,"Using chromosomal position from report file"))
      chrompos.fromReport<-TRUE
    }
    # Import data from BeadStudio report file
    firstfield <- scan(reportfile, what = "", sep = ",", flush = TRUE,
            quiet = TRUE, blank.lines.skip = FALSE, multi.line = FALSE,nlines=100)
    skip <- grep("[Data]", firstfield, fixed=TRUE)
    if (length(skip) == 0) stop("Cannot find \"[Data]\" in report file")
    alldata<-read.table(reportfile, skip=skip, header = TRUE, sep = "\t", as.is = TRUE, check.names = FALSE,colClasses="character")
    # Integrity checks
    essentialcols<-c("SNP Name","Sample ID","GC Score","GT Score","X Raw","Y Raw")
    if (chrompos.fromReport) essentialcols<-c(essentialcols,c("Chr","Position"))
    foundcols<-essentialcols %in% colnames(alldata)
    if (!all(foundcols)) stop ("Columns:",paste(essentialcols[!foundcols],collapse=", ")," are missing in the report file")
    # Test availability of Genotyping information
    ab.report<-all(c("Allele1 - AB","Allele2 - AB") %in% colnames(alldata))
    if (!ab.report){
      alleliccols<-c("Allele1 - Top","Allele2 - Top")
      if (chrompos.fromReport) alleliccols<-c(alleliccols,"ILMN Strand","SNP")
      if (!all(alleliccols %in% colnames(alldata))) stop("Insufficient information on allelic state in the reportfile, see 'read.SnpSetIllumina' help page")
    }
    #TODO: Could compare this to information from reportfile header
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
    # Combine results from single alleles
    if (ab.report) {
      GenCall<-paste(alldata[,"Allele1 - AB"],alldata[,"Allele2 - AB"],sep="")
    } else {
      # Convert SNPs with ACGT notation
      GenCall1<-alldata[,"Allele1 - Top"]
      GenCall2<-alldata[,"Allele2 - Top"]
      if (chrompos.fromReport) {
        TopBot<-alldata[,"ILMN Strand"]=="BOT"
        Snp1<-substr(alldata[,"SNP"],2,2)
      } else {
        SNPinfo<-SNPinfo[rownames(G),]
        TopBot<-rep(SNPinfo[,"IllStrand"]=="BOT",n)
        Snp1<-rep(substr(SNPinfo[,"snpbases"],2,2),n)
      }
      # Generate complement for SNPs on bottom strand
      GenCall1[TopBot]<-c("A","C","G","T","-")[match(GenCall1[TopBot],c("T","G","C","A","-"))]
      GenCall2[TopBot]<-c("A","C","G","T","-")[match(GenCall2[TopBot],c("T","G","C","A","-"))]
      # convert ACGT to AB
      GenCall1<-ifelse(GenCall1==Snp1,"A","B")
      GenCall1[GenCall2=="-"]<-"-"  # preserve no call
      #
      GenCall2<-ifelse(GenCall2==Snp1,"A","B")
      GenCall2[GenCall1=="-"]<-"-"
      #
      GenCall<-paste(GenCall1,GenCall2,sep="")
    }
    GenCall<-matrix(c("A","H","B","-")[match(GenCall,c("AA","AB","BB","--"))],nrow=m,ncol=n,byrow=FALSE,dimnames=newdimnames)
    extraData<-list()
    extraSNPinfo<-NULL
    extraColumns<-colnames(alldata)[!colnames(alldata) %in% c(essentialcols,"Allele1 - AB","Allele2 - AB")]
    for (nam in extraColumns) {
      mat<-matrix(type.convert(alldata[,nam],as.is=TRUE),nrow=m,ncol=n,byrow=FALSE,dimnames=newdimnames)
      # if all columns have same information -> this is an annotation column, put in SNPinfo
      if (all(apply(mat,1,function(x) length(unique(x))==1))) {
        # just use the first column
        extraSNPinfo<-cbind(extraSNPinfo,mat[,1])
        colnames(extraSNPinfo)[ncol(extraSNPinfo)]<-nam
      } else {
        extraData[[nam]]<-mat
      }
    }
    if (is.null(names(extraData))) extraData<-NULL
    if (!is.null(extraSNPinfo)) rownames(extraSNPinfo)<-newdimnames[[1]]

    # Select only data from samplesheet, try a few different ways
    if (all(rownames(samples) %in% colnames(G))) {
      selected<-rownames(samples)
    } else if (all(samples[, "Sample_Name"] %in% colnames(G))){
      selected<-samples[, "Sample_Name"]
    } else if (any(rownames(samples) %in% colnames(G))){
      selected<-rownames(samples)[rownames(samples) %in% colnames(G)]
    } else if (any(samples[, "Sample_Name"] %in% colnames(G))){
      selected<-samples[samples[, "Sample_Name"] %in% colnames(G), "Sample_Name"]
    } else { # default naming of columns in beadstudio
      selected <- paste(samples[, "Sentrix_ID"], samples[,"Sentrix_Position"], sep = "_")
    }
    #
    if (length(selected)==0) stop("None of the samples in the samplesheet match with the samples in the reportfile")
    if (length(selected)<nrow(samples)) {
      samples<-samples[selected,]
      warning("Only a subset of the samples in the samplesheet could be found in the reportfile")
    }
    G<-G[,selected]
    R<-R[,selected]
    GenScore<-GenScore[,selected]
    GenCall<-GenCall[,selected]
    if (!is.null(extraData)) {
      extraData<-lapply(extraData, function(obj) obj[, selected])
    }
    colnames(G)<-rownames(samples)
    colnames(R)<-rownames(samples)
    colnames(GenScore)<-rownames(samples)
    colnames(GenCall)<-rownames(samples)
    if (chrompos.fromReport) {
      CHR<-matrix(alldata[,"Chr"],nrow=m,ncol=n,byrow=FALSE,dimnames=newdimnames)
      MapInfo<-matrix(as.numeric(alldata[,"Position"]),nrow=m,ncol=n,byrow=FALSE,dimnames=newdimnames)
      SNPinfo<-data.frame(CHR[,1,drop=FALSE],MapInfo[,1],GTS,OPA=rep(OPAname[1],length.out=nrow))
      # assigning colnames in previous line (data.frame() ) does not work because the first column retains its column name through drop=FALSE
      # drop=FALSE is used to retain the rownames
      colnames(SNPinfo)<-c("CHR","MapInfo","GTS","OPA")
    } else {
      SNPinfo<-SNPinfo[rownames(G),]
      SNPinfo<-cbind(SNPinfo,GTS=GTS)
    }
    if (!is.null(extraSNPinfo)) SNPinfo<-cbind(SNPinfo,extraSNPinfo)
    samples[,"validn"]<-apply(G,2,function(x) sum(!is.na(x)))
  }
  new("SnpSetIllumina",phenoData=new("AnnotatedDataFrame",samples,data.frame(labelDescription=colnames(samples),row.names=colnames(samples))), annotation=OPAname, call=GenCall, callProbability=GenScore, G=G, R=R,
             featureData=new("AnnotatedDataFrame",SNPinfo,data.frame(labelDescription=colnames(SNPinfo),row.names=colnames(SNPinfo))),extraData=extraData,storage.mode="list")
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

IlluminaGetOPAinfo<-function(OPAnames,OPAinfoPath,brief=TRUE) {
	# find .opa file
	# Will be used for SNPinfo field in SNPlist
	# Some columns with special meaning
  # CHR : chromosome
  # MapInfo : position on chrom
  # OPA : linkage panel
  # snpid : international snp id rsxxxxxx (used as row names)
  # IllCode : numeric within linkage panel to connect to snpid
  # Combine multiple OPA files
  illOPA<-NULL
  for (OPAname in OPAnames) {
  	OPAfile<-list.files(OPAinfoPath,pattern=paste(OPAname,".*\\.opa$",sep=""),full.names=TRUE)
  	if (length(OPAfile) == 1) {
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
      OPA<-rep(OPAname,nrow(OPAinfo))
      #TODO check existing snp_ids in illOPA
      if (brief) illOPA<-rbind(illOPA,cbind(OPA=I(OPA),OPAinfo[,c("snpid","IllCode","CHR","MapInfo","IllStrand","snpbases")]))
      else illOPA<-rbind(illOPA,cbind(OPA=I(OPA),OPAinfo))
    }
  }
  return(illOPA)
}

IlluminaGetGencalls<- function(path,OPAname) {
  gencallfile<-list.files(path,pattern=paste(OPAname,".*LocusByDNA.*csv",sep=""),full.names=TRUE)
  # exclude some common report types *DNA_Report.csv/*Locus_Report.csv/*Final.csv
  exclude<-grep("DNA_Report\\.csv|Locus_Report\\.csv|Final\\.csv",gencallfile)
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

Sample_Map2Samplesheet<-function(samplemapfile,saveas="") {
  # read a SampleMAP file created by beadstudio and convert it to a samplesheet that is usable for read.SnpSetIllumina
  # just assume this is tab-separated
  ss<-read.table(samplemapfile,sep="\t",as.is=TRUE,header = TRUE, check.names = FALSE,colClasses="character")
  req_cols<-c("Index","Name","ID","Plate","SentrixPosition")
  smpcol<-!(req_cols %in% colnames(ss))
  if (any(smpcol)) 
    stop("Sample Map error, column(s) ",paste(req_cols[smpcol],collapse=", ")," are not available")
  colnames(ss)[match(c("Name","Plate","SentrixPosition"),colnames(ss))]<-c("Sample_Name","Sample_Plate","Sentrix_Position")
  # add some required columns
  ss<-cbind(ss,Pool_ID="pool_id",Sentrix_ID="sentrix_id")
  # fill in empty sample_names
  empty_sampleName<-ss$Sample_Name==""
  if (any(empty_sampleName)) ss$Sample_Name[empty_sampleName]<-ss$ID
  if (! saveas=="") {
    cat("[Header]\nInvestigator Name,\nProject Name,\nExperiment Name,\nDate,\n[Data]\n",file=saveas)
    write.table(ss,row.names=FALSE,sep=",",quote=FALSE,file=saveas,append=TRUE)
  }
  return(ss)
}

readClusteringData<-function(object,SNPMap) {

}