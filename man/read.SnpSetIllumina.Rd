\name{read.SnpSetIllumina}
\alias{read.SnpSetIllumina}
\concept{SnpSetIllumina}

\title{Read Experimental Data and into an 'SnpSetIllumina' Object}

\description{
   A SnpSetIllumina object is created from the textfiles created by the Illumina 
	 GenCall or BeadStudio software.
}

\usage{
   read.SnpSetIllumina(samplesheet, manifestpath = NULL, reportpath = NULL,
     rawdatapath = NULL, reportfile = NULL, briefOPAinfo=TRUE, verbose=FALSE,
     readTIF=FALSE, ...)
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
   \item{briefOPAinfo}{logical, if \code{TRUE} then only the SNP name, Illumi code, 
     chromosome and basepair position are put into the \code{featureData} slot of
     the result, else all information from the OPA file is put into the \code{featureData}
     slot}
   \item{verbose}{logical, if \code{TRUE} then some extra information is given 
     during the import}
   \item{readTIF}{logical, uses \code{beadarray} package and raw TIF files to read data}
   \item{\dots}{arguments are forwarded to \code{readIllumina} and can be used to perform
     bead-level normalization} 
}

\details{
  The text files from Illumina software are imported to a SnpSetIllumina object. 
	Both result files from GenCall and BeadStudio can be used. 
  In both cases the sample sheets from the experiments are used to select the 
	proper data from the report or data files. The following columns from the 
	sample sheet file are used for this purpose: \code{'Sample_Name'}, 
	\code{'Sentrix_Position'}, and \code{'Pool_ID'}. The values in columns 
  \code{'Sample_Plate'}, \code{'Pool_ID'}, and \code{'Sentrix_ID'} should be the 
	same for all samples in the file, as this is the case for processed 
  experiments. The contents of the sample sheet are put into the phenoData slot. 
	Also the OPA definition file containing SNP annotation should be available, 
	these files are provided by Illumina. Columns \code{'IllCode'}, \code{'CHR'}, 
	and \code{'MapInfo'} are put into the \code{featureData} slot.
  
  In order to process experiments that were genotyped using the GenCall software, 
	the arrays should be scanned with the setting \code{<SaveTextFiles>true</SaveTextFiles>} 
	in the Illumina configuration file Settings.XML.
  3 Types of files need to be present in the same folder: The sample sheet, 
	.csv files containing signal intensity data, and the report file that contains 
	the genotype information. For each sample in the sample sheet there should be 
	a .csv file with the following file mask: \code{[sam_id]_R00[yy]_C00[xx].csv}, 
  where \code{sam_id} is the Illumina ID for the SAM, and \code{xx} and \code{yy} 
	are the column and row number respectively. From the report files the file 
	with mask \code{[Pool_ID]_LocusByDNA[_ExpName].csv} is used. \code{'Pool_ID'} 
	is the OPA panel used, and \code{'_ExpName'} is optional.
  
  To process experiments that were processed with BeadStudio, only two files are 
	needed. The sample sheet and the Final Report file. The sample sheet must 
	contain the same columns as for GenCall, the report file should contain the 
	following columns: \code{'SNP Name'}, \code{'Sample ID'}, \code{'GC Score'}, 
	\code{'Allele1 - AB'}, \code{'Allele2 - AB'}, \code{'GT Score'}, \code{'X Raw'}, 
	and \code{'Y Raw'}. \code{'SNP Name'} and \code{'Sample ID'} are used to form 
	rows and columns in the experimental data, \code{'GC Score'} is put in the 
	\code{callProbability} matrix, \code{'Allele1 - AB'} and \code{'Allele2 - AB'} 
	are combined into the \code{call} matrix, \code{'GT Score'} is added to the 
	\code{featureData} slot, \code{'X Raw'} is put in the \code{R} matrix and 
	\code{'Y Raw'} in the \code{G} matrix. Other columns in the report file are 
	added as matrices in the \code{assayData} slot.
  
}
\value{
  This function returns an \code{SnpSetIllumina} object.
}
\author{Jan Oosting}

\seealso{
   \code{\link{SnpSetIllumina-class}}
}
\examples{
# read a SnpSetIllumina object using example textfiles in data directory
datadir <- system.file("testdata", package="beadarraySNP")
SNPdata <- read.SnpSetIllumina(paste(datadir,"4samples_opa4.csv",sep="/"),datadir)
}
\keyword{file}
\keyword{manip}
