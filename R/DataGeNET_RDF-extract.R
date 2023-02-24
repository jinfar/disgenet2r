#' Obtain the raw query from DisGeNET from a \code{DataGeNET.RDF} object.
#'
#' @name extract
#' @rdname extract-methods
#' @aliases DataGeNET.Dis-methods
#' @param object Object of class \code{DataGeNET.RDF}
#' @return A \code{data.frame} containing the raw result from DisGeNET
#' @examples
#' \dontrun{
#' #Being x an DataGeNET.RDF
#' qr <- extract(x) 
#' }
#' @export
setMethod( "extract",
           signature = "DataGeNET.RDF",
           definition = function( object ) {
             return( object@qresult )
           }
)