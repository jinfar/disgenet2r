#' Retrieves evidences supporting a gene-disease association, for a gene, or list of genes and generates an \code{DataGeNET.DGN}
#'
#' @param variant a variant or a list of variants
#' @param database Name of the database that will be queried. It can take the values:
#' \code{'CTD_human'} to use Comparative Toxicogenomics Database, human data;
#' \code{'UNIPROT'} to use Universal Protein Resource;
#' \code{'CLINGEN'} to use Clinical Genome Resource;
#' \code{'CGI'} to use Cancer Genome Interpreter;
#' \code{'ORPHANET'}, to use Orphanet, the portal for rare diseases and orphan drugs;
#' \code{'PSYGENET'} to use PSYGENET;
#' \code{'GENOMICS_ENGLAND'} to use Genomics England PanelApp;
#' \code{'CURATED'} to use expert curated, human databases;
#' \code{'HPO'} to use HPO;
#' \code{'INFERRED'} to use inferred data from HPO, GWASDB, GWASCAT, and CLINVAR;
#' \code{'CTD_rat'} to use Comparative Toxicogenomics Database, rat data;
#' \code{'CTD_mouse'} to use Comparative Toxicogenomics Database, mouse data;
#' \code{'RGD'}, to use Rat Genome Database;
#' \code{'MGD'}, to use the Mouse Genome Database;
#' \code{'ANIMAL_MODELS'} to use the expert curated, animal models data;
#' \code{'GWASCAT'} to use the NHGRI-EBI GWAS Catalog;
#' \code{'GWASDB'} to use the GWAS Database GWASdb;
#' \code{'CLINVAR'} to use ClinVar, a public archive of relationships
#' among sequence variation and human phenotype;
#' \code{'BEFREE'} to use text mining data, generated using BeFree System;
#' \code{'ALL'} to use all these databases. Default \code{'CURATED'}.
#' @param score A vector with two elements: 1) initial value of score 2) final value of score
#' @param verbose By default \code{FALSE}. Change it to \code{TRUE} to get a
#' on-time log from the function.
#' @param warnings By default \code{TRUE}. Change it to \code{FALSE} to not see
#' the warnings.
#' @return An object of class \code{DataGeNET.Dis}
#' @examples
#' g <- variant2evidence( variant = "rs121913279", "CURATED", score=c(0,1)  )
#' @export variant2evidence


variant2evidence <- function( variant,
                           disease = NULL, year=c(0,2050),
                           database = "CURATED", score=c(0,1), api_key=NULL,  verbose = FALSE, warnings = TRUE ) {
  if(is.null(api_key)){
    api_key = Sys.getenv('DISGENET_API_KEY')
    if(api_key == ""){
      stop("This is not a valid API KEY! Please, use the function `get_disgenet_api_key` to get your API key from DisGeNET")
    }
  }
  check_disgenet_sources( database )

  if( length( variant ) != length( unique( variant ) ) ) {
    variant <- unique( variant )
    warning(
      "Removing duplicates from input genes list."
    )
  }


  if(length(score) != 2) {
    stop("Invalid argument 'score'. It must have two elements: initial value of the score, and final value of the score")
  }
  if(length(year) != 2) {
    stop("Invalid argument 'year'. It must have two elements: initial year and final year")
  }

  list_of_variants <- paste(variant,collapse=",")
  list_of_diseases <- paste(disease,collapse=",")

  result <- NULL
  if(!is.null(disease)){
    url <- paste0( get_url_disgenet(), "vda/evidences/variant/", list_of_variants , "?disease=",list_of_diseases , "&source=",database, "&min_score=",score[1],"&max_score=", score[2], "&min_year=",year[1],"&max_year=", year[2])
  }
  else{
    url <- paste0( get_url_disgenet(), "vda/evidences/variant/", list_of_variants , "?source=",database, "&min_score=",score[1],"&max_score=", score[2], "&min_year=",year[1],"&max_year=", year[2])

  }
  r <- get_api_connection(url = url, api_key=api_key)
  if (r$status_code == 200) {

    res<-jsonlite::fromJSON(httr::content(r, as = "text", encoding = "UTF-8"), flatten = F)
    result <- res$results

    while (is.null(res$`next` )==F) {
      r <- get_api_connection(url = res$`next`, api_key=api_key)
      if (r$status_code == 200) {

        res <- jsonlite::fromJSON(httr::content(r, as = "text", encoding = "UTF-8"), flatten = F)
        result <- rbind(result, res$results)
      }
    }

  } else{
    print(httr::http_status(r))
    print(httr::content(r, "text"))
  }



  if(length( variant )==1){
    cls <- "single"
  }else{
    cls <- "list"
  }

  scR <- paste0( score[1],"-", score[2])

  if (! is.null(result) > 0){
    missing_variants <- setdiff(variant, result$variant_id)
    found_variants <- intersect(variant, result$variant_id)

    if( length( missing_variants) > 0 ) {
      missing_variants <- paste( paste( "   -", missing_variants ), collapse = "\n" )
      if( warnings ) {
        warning( "\n One or more of the variants in the list is not in DisGeNET ( '", database, "' ):\n", missing_variants )
      }
    }

    result <- new( "DataGeNET.DGN",
                   type = "evidence",
                   search = cls,
                   term = as.character( found_variants ),
                   scoreRange = scR,
                   database = database,
                   qresult = result )

    return( result )
  } else {
    print("no results for the query")
  }
}
