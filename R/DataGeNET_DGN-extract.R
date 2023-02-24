#' Obtain the raw query from DisGeNET from a \code{DataGeNET.DGN} object.
#'
#' @name extract
#' @rdname extract-methods
#' @aliases DataGeNET.DGN-methods
#' @param object Object of class \code{DataGeNET.DGN}
#' @return A \code{data.frame} containing the raw result from DisGeNET
#' @examples
#' \dontrun{
#' #Being x an DataGeNET.DGN
#' qr <- extract(x) # Get number of unique diseases
#' }
#' @export
setMethod( "extract",
   signature = "DataGeNET.DGN",
   definition = function( object ) {
     return( object@qresult )
   }
)
