#' Obtain the number of unique diseases in a \code{DataGeNET.DGN}.
#'
#' @name ndisease
#' @rdname ndisease-methods
#' @aliases DataGeNET.DGN-methods
#' @param object Object of class \code{DataGeNET.DGN}
#' @return The number of unique diseases
#' @examples
#' \dontrun{
#' #Being x an DataGeNET.DGN
#' nd <- ndisease(x) # Get number of unique diseases
#' }
#' @export
setMethod( "ndisease",
  signature = "DataGeNET.DGN",
  definition = function( object ) {
    return( length( unique( object@qresult$diseaseid ) ) )
  }
)
