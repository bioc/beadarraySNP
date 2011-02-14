# QCIllumina
setGeneric("arrayType",           function(object) standardGeneric("arrayType"))
setGeneric("arrayType<-",         function(object, value) standardGeneric("arrayType<-"))
setGeneric("arrayID",             function(object) standardGeneric("arrayID"))
setGeneric("arrayID<-",           function(object, value) standardGeneric("arrayID<-"))
setGeneric("plotQC",              function(object, type) standardGeneric("plotQC"))
setGeneric("reportSamplePanelQC", function(object, by=10, legend=TRUE, ...) standardGeneric("reportSamplePanelQC"))
# SnpSetIllumina
setGeneric("calculateGSR",        function(object) standardGeneric("calculateGSR"))
setGeneric("calculateSmooth",     function(object, smoothType, ...) standardGeneric("calculateSmooth"))
setGeneric("sortGenomic",         function(object) standardGeneric("sortGenomic"))
#SnpSetSegments
setGeneric("cn.segments",         function(object) standardGeneric("cn.segments"))
setGeneric("cn.segments<-",       function(object, value) standardGeneric("cn.segments<-"))

