\name{SnpSetIllumina}
\docType{class}
\alias{class:SnpSetIllumina}
\alias{SnpSetIllumina}
\alias{SnpSetIllumina-class}
\alias{exprs<-,SnpSetIllumina,matrix-method}
\alias{reporterInfo,SnpSetIllumina-method}
\alias{reporterInfo<-,SnpSetIllumina,data.frame-method}
\alias{exprs,SnpSetIllumina-method}
\alias{[,SnpSetIllumina-method}
\alias{initialize,SnpSetIllumina-method}
\alias{combine,SnpSetIllumina,SnpSetIllumina-method}
\alias{combine,SnpSetIllumina,ANY-method}
\alias{featureData<-,SnpSetIllumina-method}
\alias{featureData,SnpSetIllumina-method}
\alias{calculateGSR,SnpSetIllumina-method}
\alias{calculateGSR}

\title{
  Class to Contain Objects Describing High-Throughput SNP Assays.
}

\description{
  Container for high-throughput assays and experimental
  metadata. \code{SnpSetIllumina} class is derived from
  \code{\link[Biobase:eSet]{eSet}}, and requires matrices \code{R}, \code{G}, \code{call},
  \code{callProbability} as assay data members.  
  
  It supports featureData. Several visualization methods use columns \code{CHR} 
  and \code{MapInfo}. The \code{CHR} column is used to handle sex chromosomes in 
  a specific way. The \code{OPA} column is the default way to specify subsamples. } 


\section{Extends}{
   Directly extends class \code{\link[Biobase:eSet]{eSet}}.
}
\section{Creating Objects}{

  \code{new('SnpSetIllumina',
    phenoData = [AnnotatedDataFrame],
    experimentData = [MIAME],
    annotation = [character],
    call = [matrix],
    callProbability = [matrix],
    G = [matrix],
    R = [matrix],
    featureData = [data.frameOrNULL],
    ...)
  }

  \code{SnpSetIllumina} instances are usually created through
  \code{new("SnpSetIllumina", ...)}. Arguments to \code{new}
  include \code{call} (a matrix of gentoypic calls, with features (SNPs)
  corresponding to rows and samples to columns), \code{callProbability},
  \code{G}, \code{R}, \code{phenoData},
  \code{experimentData}, and \code{annotation}. 
  \code{phenoData}, \code{experimentData}, and \code{annotation} can be
  missing, in which case they are assigned default values.
}
\section{Slots}{
  Inherited from \code{Biobase:eSet}:
  \describe{
    \item{\code{assayData}:}{Contains matrices with equal dimensions, and with 
      column number equal to \code{nrow(phenoData)}. \code{assayData} must 
      contain a matrix \code{call} with rows representing features (e.g., SNPs) 
      and columns representing samples, a matrix \code{callProbability} 
      describing the certainty of the call, and matrices \code{R} and \code{G} 
      to describe allele specific intensities. The contents of these matrices 
      are	not enforced by the class. The \code{assayData} matrices \code{Gb, Rb, 
      intensity, theta} are optional, but are either results or input for 
      several methods of the class. Additional matrices of identical size may 
      also be included in \code{assayData}. Class:\code{\link[Biobase:AssayData]{AssayData}}}
    \item{\code{phenoData}:}{See \code{\link[Biobase:eSet]{eSet}}}
    \item{\code{experimentData}:}{See \code{\link[Biobase:eSet]{eSet}}}
    \item{\code{annotation}:}{See \code{\link[Biobase:eSet]{eSet}}}
    \item{\code{featureData}:}{annotation for SNPs, usually will contain a 
      \code{CHR} and a \code{MapInfo} column for genomic localization }
  }
}
\section{Methods}{

  Class-specific methods:
  \describe{
     \item{\code{exprs(SnpSetIllumina)}, \code{exprs(SnpSetIllumina,matrix)<-}}{Access and
       set elements named \code{call} in the \code{AssayData}
       slot.}
     \item{\code{combine(SnpSetIllumina,SnpSetIllumina)}:}{performs \code{union}-like 
		   combination in both dimensions of SnpSetIllumina objects}
     \item{\code{reporterInfo(SnpSetIllumina)}, \code{reporterInfo(SnpSetIllumina,data.frame)<-}}{Access and
       set the \code{pData} in the featureData slot}
     \item{\code{calculateGSR(SnpSetIllumina)}}{calculate ratio of Gentrain score and Gencall score. Creates 
       GSR matrix in \code{assayData}. Should be performed before combining datasets}
  }
  
   Derived from \code{\link[Biobase:eSet]{eSet}}:
   \describe{
     \item{\code{sampleNames(SnpSetIllumina)} and \code{sampleNames(SnpSetIllumina)<-}:}{See \code{\link[Biobase:eSet]{eSet}}}
     \item{\code{featureNames(SnpSetIllumina)}, \code{featureNames(SnpSetIllumina, value)<-}:}{See \code{\link[Biobase:eSet]{eSet}}}
     \item{\code{dims(SnpSetIllumina)}:}{See \code{\link[Biobase:eSet]{eSet}}}
     \item{\code{phenoData(SnpSetIllumina)}, \code{phenoData(SnpSetIllumina,value)<-}:}{See \code{\link[Biobase:eSet]{eSet}}}
     \item{\code{varLabels(SnpSetIllumina)}, \code{varLabels(SnpSetIllumina, value)<-}:}{See \code{\link[Biobase:eSet]{eSet}}}
     \item{\code{varMetadata(SnpSetIllumina)}, \code{varMetadata(SnpSetIllumina,value)<-}:}{See \code{\link[Biobase:eSet]{eSet}}}
     \item{\code{pData(SnpSetIllumina)}, \code{pData(SnpSetIllumina,value)<-}:}{See \code{\link[Biobase:eSet]{eSet}}}
     \item{\code{varMetadata(SnpSetIllumina)}, \code{varMetadata(SnpSetIllumina,value)}}{See \code{\link[Biobase:eSet]{eSet}}}
     \item{\code{experimentData(SnpSetIllumina)},\code{experimentData(SnpSetIllumina,value)<-}:}{See \code{\link[Biobase:eSet]{eSet}}}
     \item{\code{pubMedIds(SnpSetIllumina)}, \code{pubMedIds(SnpSetIllumina,value)}}{See \code{\link[Biobase:eSet]{eSet}}}
     \item{\code{abstract(SnpSetIllumina)}:}{See \code{\link[Biobase:eSet]{eSet}}}
     \item{\code{annotation(SnpSetIllumina)}, \code{annotation(SnpSetIllumina,value)<-}}{See \code{\link[Biobase:eSet]{eSet}}}
     \item{\code{storageMode(eSet)}, \code{storageMode(eSet,character)<-}:}{See \code{\link[Biobase:eSet]{eSet}}}
     \item{\code{reporterNames(SnpSetIllumina)}, \code{reporterNames(SnpSetIllumina,value)<-}:}{DEPRECATED}
     \item{\code{featureData(SnpSetIllumina)}, \code{featureData(SnpSetIllumina,AnnotatedDataFrame)<-}}{See \code{\link[Biobase:eSet]{eSet}}}
     \item{\code{object[(index)}:}{Conducts subsetting of matrices and phenoData and featureData components }
  }

  Standard generic methods:
  \describe{
    \item{\code{initialize(SnpSetIllumina)}:}{Object instantiation, used
      by \code{new}; not to be called directly by the user.}
    \item{\code{validObject(SnpSetIllumina)}:}{Validity-checking method, ensuring
      that \code{call}, \code{callProbability}, \code{G}, and \code{R} are members of
      \code{assayData}. \code{checkValidity(SnpSetIllumina)} imposes this
      validity check, and the validity checks of \code{Biobase:eSet}.}
    \item{\code{show(SnpSetIllumina)}}{See \code{\link[Biobase:eSet]{eSet}}}
    \item{\code{dim(SnpSetIllumina)}, \code{ncol}}{See \code{\link[Biobase:eSet]{eSet}}}
    \item{\code{SnpSetIllumina[(index)}:}{See \code{\link[Biobase:eSet]{eSet}}}
    \item{\code{SnpSetIllumina$}, \code{SnpSetIllumina$<-}}{See \code{\link[Biobase:eSet]{eSet}}}
  }
}
\author{J. Oosting, based on Biobase eSet class }

\seealso{
  \code{\link[Biobase:eSet]{eSet}}}
}

\keyword{classes}
