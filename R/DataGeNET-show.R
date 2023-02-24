setMethod( "show",
           signature = "DataGeNET.DGN",
           definition = function( object ) {
             cat( "Object of class 'DataGeNET.DGN'\n" )
             cat( " . Search:     ", object@search, "\n" )
             cat( " . Type:       ", object@type, "\n" )
             cat( " . Database:    ", object@database, "\n" )
             cat( " . Score:       ", object@scoreRange, "\n" )
             if( object@search == "single" ) {
               cat( " . Term:       ", object@term, "\n" )
             } else {
               n <- length( object@term )
               cat( " . Term:      ", object@term[ 1 ], "...", object@term[ n ], "\n" )
             }
             if ( object@type == "gene-disease"){
               cat( " . Results: ", length( unique( object@qresult$diseaseid ) ), "\n" )
             } else if (  object@type == "variant-disease"){
               cat( " . Results: ", length( unique( object@qresult$diseaseid ) ), "\n" )
             }
             else if ( object@type == "disease-variant" ){
               cat( " . Results: ", length( unique( object@qresult$variantid)), "\n" )
             }
             else if( object@type == "disease-gene"){
               cat( " . Results: ", length( unique( object@qresult$geneid ) ), "\n" )
             }
             else if( object@type == "disease-disease-gene" | object@type == "disease-disease-variant"){
               cat( " . Results: ",  dim( unique( object@qresult ) )[1], "\n" )
             }
           }
)

setMethod( "show",
           signature = "DataGeNET.RDF",
           definition = function( object ) {
             cat( "Object of class 'DataGeNET.RDF'\n" )
             cat( " . Term:      ", object@input, "\n" )
             cat( " . Type:     ", object@search, "\n" )
             cat( " . Database:    ", object@database, "\n" )
             cat( " . Score:       ", object@scoreRange, "\n" )
             cat( " . Selection:  ", object@selection, "\n" )
             cat( " . Results:    ", object@mapping, "\n" )
           }
)

