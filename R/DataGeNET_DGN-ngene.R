#' Obtain the number of unique genes in a \code{DataGeNET.DGN}.
#'
#' @name ngene
#' @rdname ngene-methods
#' @aliases DataGeNET.DGN-methods
#' @param object Object of class \code{DataGeNET.DGN}
#' @return The number of unique genes
#' @examples
#' \dontrun{
#' #Being x an DataGeNET.DGN
#' nd <- ngene(x) # Get number of unique genes
#' }
#' @export
setMethod( "ngene",
  signature = "DataGeNET.DGN",
  definition = function( object ) {
    return( length( unique( object@qresult$geneid ) ) )
  }
)
