setMethod("initialize", "SnpSetSegments",
          function(.Object,
                  assayData = assayDataNew(call = call,
                                           callProbability = callProbability,
                                           G = G, R = R, ...),
                   phenoData = new("AnnotatedDataFrame"),
                   experimentData = new("MIAME"),
                   annotation = character(),
                   protocolData = phenoData[,integer(0)],
                   call = new("matrix"),
                   callProbability = new("matrix"),
                   G = new("matrix"),
                   R = new("matrix"),
                   cn.segments = list(),
                   featureData = new("AnnotatedDataFrame"),
                   extraData = NULL,
                   ...) {
            .Object<-callNextMethod(.Object,
                           assayData = assayData,
                           phenoData = phenoData,
                           experimentData = experimentData,
                           annotation = annotation,
                           featureData = featureData,
                           protocolData = protocolData,
                           extraData = extraData)
            cn.segments(.Object)<-cn.segments
            validObject(.Object)
            .Object

          })

setValidity("SnpSetSegments",
  function(object) {
    if (length(cn.segments(object))>0) {

    }
  }
)


setReplaceMethod("cn.segments",
  signature=signature(
    object="SnpSetSegments",
    value="list"
  ),
  function(object, value) {
  	object@cn.segments<-value
	  object
  }
)

setMethod("cn.segments", "SnpSetSegments", function(object) object@cn.segments)

