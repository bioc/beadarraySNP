\name{read.SnpSetIllumina}
\alias{read.SnpSetIllumina}
\concept{SnpSetIllumina}

\title{Read Experimental Data and into an 'SnpSetIllumina' Object}

\description{
   A SnpSetIllumina object is created from the textfiles created by the Illumina 
	 GenCall or BeadStudio software.
}

\usage{
   read.SnpSetIllumina(samplesheet, manifestpath=NULL, reportpath=NULL,
     rawdatapath=NULL, reportfile=NULL, briefOPAinfo=TRUE, readTIF=FALSE, 
     nochecks=FALSE, sepreport="\t", ...)
}
\arguments{
   \item{samplesheet}{a data.frame or filename, contains the sample sheet}
   \item{manifestpath}{a character string for the path containing the 
     manifests / OPA definition files, defaults to path of samplesheet}
   \item{reportpath}{a character string for the path containing the 
     report files, defaults to path of samplesheet}
   \item{rawdatapath}{a character string for the path containing the 
     intensity data files, defaults to path of samplesheet}
   \item{reportfile}{a character string for the name of BeadStudio reportfile}
   \item{briefOPAinfo}{logical, if \code{TRUE} then only the SNP name, Illumi 
     code, chromosome and basepair position are put into the \code{featureData} 
     slot of the result, else all information from the OPA file is put into the 
     \code{featureData} slot}
   \item{readTIF}{logical, uses \code{beadarray} package and raw TIF files to 
     read data}
   \item{nochecks}{logical, limited validity checks on beadstudio report files. See details}
   \item{sepreport}{character, field separator character for beadstudio report files}
   \item{\dots}{arguments are forwarded to \code{\link[beadarray:readIllumina]{readIllumina}} 
     and can be used to perform bead-level normalization} 
}

\details{
  The text files from Illumina software are imported to a SnpSetIllumina object. 
	Both result files from GenCall and BeadStudio can be used. 
  In both cases the sample sheets from the experiments are used to select the 
	proper data from the report or data files. The following columns from the 
	sample sheet file are used for this purpose: \sQuote{\code{Sample_Name}}, 
	\sQuote{\code{Sentrix_Position}}, and \sQuote{\code{Pool_ID}}. The values in 
  columns \sQuote{\code{Sample_Plate}}, \sQuote{\code{Pool_ID}}, and 
  \sQuote{\code{Sentrix_ID}} should be the same for all samples in the file, as 
  this is the case for processed experiments. The contents of the sample sheet 
  are put into the \code{phenoData} slot. 
  
  Zero values in the raw data signals are set to NA
  
	Ideally the OPA manifest file containing SNP annotation should be available, 
	these files are provided by Illumina. Columns \sQuote{\code{IllCode}}, 
  \sQuote{\code{CHR}}, and \sQuote{\code{MapInfo}} are put into the 
  \code{featureData} slot.
  
  \emph{GenCall Data}
  
  In order to process experiments that were genotyped using the GenCall software, 
	the arrays should be scanned with the setting 
  \code{<SaveTextFiles>true</SaveTextFiles>} in the Illumina configuration file 
  \code{Settings.XML}. 3 Types of files need to be present in the same folder: 
  The sample sheet, .csv files containing signal intensity data, and the report 
  file that contains the genotype information. For each sample in the sample 
  sheet there should be a .csv file with the following file mask: 
  \code{[sam_id]_R00[yy]_C00[xx].csv}, where \code{sam_id} is the Illumina ID 
  for the SAM, and \code{xx} and \code{yy} are the column and row number 
  respectively. From the report files the file with mask 
  \code{[Pool_ID]_LocusByDNA[_ExpName].csv} is used. \sQuote{\code{Pool_ID}} is 
  the OPA panel used, and \sQuote{\code{_ExpName}} is optional.
  
  \emph{BeadStudio Data}
  
  To process experiments that were processed with BeadStudio, only two files are 
	needed. The sample sheet and the Final Report file. The sample sheet must 
	contain the same columns as for GenCall, the report file should contain the 
	following columns: \sQuote{\code{SNP Name}}, \sQuote{\code{Sample ID}}, 
  \sQuote{\code{GC Score}}, \sQuote{\code{Allele1 - AB}}, 
  \sQuote{\code{Allele2 - AB}}, \sQuote{\code{GT Score}}, \sQuote{\code{X Raw}}, 
	and \sQuote{\code{Y Raw}}. \sQuote{\code{SNP Name}} and 
  \sQuote{\code{Sample ID}} are used to form rows and columns in the 
  experimental data, \sQuote{\code{GC Score}} is put in the 
  \code{callProbability} matrix, \sQuote{\code{Allele1 - AB}} and 
  \sQuote{\code{Allele2 - AB}} are combined into the \code{call} matrix, 
  \sQuote{\code{GT Score}} is added to the \code{featureData} slot, 
  \sQuote{\code{X Raw}} is put in the \code{R} matrix and \sQuote{\code{Y Raw}} 
  in the \code{G} matrix. Other columns in the report file are added as matrices 
  in the \code{assayData} slot, or columns in the \code{featureData} slot if 
  values are identical for all samples in the reportfile.
  When \code{nochecks} is \code{TRUE} then only the \sQuote{\code{SNP Name}} and
  \sQuote{\code{Sample ID}} columns are required. The resulting object is now of 
  class \code{\link[Biobase:class.MultiSet]{MultiSet}}
  
  \emph{Sample sheets}
  
  To help generate a sample sheet for BeadStudio data a \code{Sample_Map.txt} 
  file can be converted to a sample sheet with the 
  \code{\link{Sample_Map2Samplesheet}} function. For Beadstudio reportfiles it is
  also possible to set samplesheet=\code{NULL}. In this case the \code{phenoData}
  slot will be fabricated from the sample names in the reportfile.
  
  \emph{Manifest/OPA/annotation files}
  
  For BeadStudio reportfiles it is not necessary to have a Manifest file if the
  columns \sQuote{\code{Chr}} and \sQuote{\code{Position}} are available in the 
  report file. Currently this is the only way to import data from Infinium 
  arrays, because Illumina does not supply Manifest files for these arrays.
}
\value{
  This function returns an \code{SnpSetIllumina} object, or a \code{MultiSet} object 
  when \code{nochecks} is \code{TRUE}.
}
\author{Jan Oosting}

\seealso{
   \code{\link{SnpSetIllumina-class}}, \code{\link{Sample_Map2Samplesheet}}, 
   \code{\link[beadarray:readIllumina]{readIllumina}}
}
\examples{
# read a SnpSetIllumina object using example textfiles in data directory
datadir <- system.file("testdata", package="beadarraySNP")
SNPdata <- read.SnpSetIllumina(paste(datadir,"4samples_opa4.csv",sep="/"),datadir)
}
\keyword{file}
\keyword{manip}
