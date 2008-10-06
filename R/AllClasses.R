setClass("SnpSetIllumina", 
   contains = "eSet"
)

setClass("QCIllumina",
   representation(
     arrayType = "character", 
     arrayID = "character",
     intensityMed = "matrix",
     greenMed = "matrix",
     redMed = "matrix",
     validn = "matrix",
     annotation = "matrix",
     samples = "matrix" ,
     ptpdiff = "matrix" ,
     callrate = "matrix" ,
     hetPerc = "matrix"
   )
)  
   