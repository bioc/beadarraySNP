# Conversion to other Bioconductor objects
# beadarraySNP package
# Author: J. Oosting
# aCGH: Class aCGH
convert2aCGH<- function(object,normalizedTo=2,doLog=TRUE,organism="hsa") {
  # Intensities to log2.ratios
  # reporterInfo to clones.info
  # pData to phenotype
  if (require(aCGH)) {
    if ( !("intensity" %in% assayDataElementNames(object))) object<-RG2polar(object)
    if (doLog)
      log2.ratios<-log2(assayData(object)[["intensity"]]/normalizedTo)
    else
      log2.ratios<-assayData(object)[["intensity"]]/normalizedTo
    clones.info<-reporterInfo(object)[,c("snpid","IllCode","CHR","MapInfo")]
    names(clones.info)<-c("Clone","Target","Chrom","kb")
    # convert Chrom
    clones.info$Chrom<-numericCHR(clones.info$Chrom)
    # convert sex chromosomes
    if (organism=="hsa") clones.info$Chrom[clones.info$Chrom >= 98] <- clones.info$Chrom[clones.info$Chrom >= 98] - 75
    if (organism=="mmu") clones.info$Chrom[clones.info$Chrom >= 98] <- clones.info$Chrom[clones.info$Chrom >= 98] - 78
    # convert pos to kb
    clones.info$kb<-clones.info$kb/1000
    create.aCGH(log2.ratios, clones.info, pData(object))
  } else stop("Package 'aCGH' is needed for conversion to a aCGH object")
}

convert2SegList <- function(object,normalizedTo=2,doLog=TRUE,organism="hsa") {
  # Intensities to log2.ratios
  # reporterInfo to clones.info
  # pData to phenotype
  if (require(snapCGH)){
    if ( !("intensity" %in% assayDataElementNames(object))) object<-RG2polar(object)
    if (doLog)
      log2.ratios<-log2(assayData(object)[["intensity"]]/normalizedTo)
    else
      log2.ratios<-assayData(object)[["intensity"]]/normalizedTo
    clones.info<-reporterInfo(object)[,c("snpid","IllCode","CHR","MapInfo")]
    names(clones.info)<-c("ID","Target","Chr","Position")
    # convert Chrom
    clones.info$Chr<-numericCHR(clones.info$Chr)
    # convert sex chromosomes
    if (organism=="hsa") clones.info$Chr[clones.info$Chr >= 98] <- clones.info$Chr[clones.info$Chr >= 98] - 75
    if (organism=="mmu") clones.info$Chr[clones.info$Chr >= 98] <- clones.info$Chr[clones.info$Chr >= 98] - 78
    segList<-list()
    segList$M.observed<-log2.ratios
    segList$genes<-clones.info
    new("SegList",segList)
  } else stop("Package 'snapCGH' is needed for conversion to a SegList object")
}

