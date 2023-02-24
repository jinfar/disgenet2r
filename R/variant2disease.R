#' Retrieves diseases associated to a variant or list of variants and generates an \code{DataGeNET.DGN}
#'
#' Given a variant or a list of variants, it retrieves the variant-disease associations
#' related to diseases in DisGeNET. It creates an object of type \code{DataGeNET.DGN}.
#'
#' @param variant a variant, or list of variants
#' @param database Name of the database that will be queried. It can take the values:
#' \code{'UNIPROT'} to use Universal Protein Resource;
#' \code{'CLINVAR'} to use ClinVar, a public archive of relationships
#' among sequence variation and human phenotype;
#' \code{'GWASCAT'} to use the NHGRI-EBI GWAS Catalog;
#' \code{'GWASDB'} to use the GWAS Database GWASdb;
#' \code{'CURATED'} to use expert curated, human databases;
#' \code{'ALL'} to use all these databases. Default \code{'CURATED'}.
#' \code{'BEFREE'} to use text mining data, generated using BeFree System;
#' @param score A vector with two elements: 1) initial value of score 2) final value of score
#' @param verbose By default \code{FALSE}. Change it to \code{TRUE} to get a
#' on-time log from the function.
#' @param warnings By default \code{TRUE}. Change it to \code{FALSE} to not see
#' the warnings.
#' @return An object of class \code{DataGeNET.DGN}
#' @examples
#' vd <- variant2disease(variant = "rs121909211", database = "CURATED",  verbose = FALSE, warnings = TRUE)
#' @export variant2disease


variant2disease <- function( variant, database = "CURATED", score=c(0,1),  api_key=NULL,  verbose = FALSE, warnings = TRUE ) {
  if(is.null(api_key)){
    api_key = Sys.getenv('DISGENET_API_KEY')
    if(api_key == ""){
      stop("This is not a valid API KEY! Please, use the function `get_disgenet_api_key` to get your API key from DisGeNET")
    }
  }
  check_disgenet_sources(database)

  if (  database == "CTD_human" | database == "CGI" | database == "CLINGEN"
        | database == "GENOMICS_ENGLAND" | database == "PSYGENET" |  database == "CTD_rat" |
        database == "CTD_mouse" |database == "RGD" | database == "MGD" |
        database == "ANIMAL_MODELS" | database == "HPO" |  database == "INFERRED" |
        database == "LHGDN"  ) {
    stop("This datasource does not annotate variants.")
  }
  if(length(score) != 2) {
    stop("Invalid argument 'score'. It must have two elements: initial value of the score, and final value of the score")
  }

  list_of_variants <- paste(variant,collapse=",")
  result <- NULL

  url <- paste0( get_url_disgenet(), "vda/variant/", list_of_variants ,"?source=",database, "&min_score=",score[1],"&max_score=", score[2],  "&format=tsv")
  r <- get_api_connection(url = url, api_key=api_key)
  if (r$status_code == 200) {
    myTextConnection <- textConnection(rawToChar(r$content ))
    result <- read.csv( myTextConnection, header = TRUE, sep = "\t" )
    close(myTextConnection)
  }else{
    print(httr::http_status(r))
    print(httr::content(r, "text"))
  }

  if(length( variant )==1){
    cls <- "single"
  }else{
    cls <- "list"
  }

  scR <- paste0( score[1],"-", score[2])

  if (!is.null(result)){
    wVariants <- setdiff(variant,result$variantid)
    eVariants <- intersect(variant, result$variantid)
    if( length( wVariants ) != 0 ) {
      variants <- paste( paste( "   -", wVariants ), collapse = "\n" )
      if( warnings ) {
        warning( "\n One or more of the variants in the list is not in DisGeNET ( '", database, "' ):\n", variants )
      }
    }

    result <- new( "DataGeNET.DGN",
                   type = "variant-disease",
                   search = cls,
                   term = as.character( eVariants ),
                   scoreRange = scR,
                   database = database,
                   qresult = result )

    return( result )
  } else {
    print("no results for the query")
  }
}
