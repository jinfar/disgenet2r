#' Retrieves genes associated to a disease, or list of diseases and generates an \code{DataGeNET.DGN}
#'
#' Given the name of one or multiple diseases retrieves their asssociated genes
#' from DisGeNET and creates an object of type \code{DataGeNET.DGN}.
#'
#' @param disease  A disease or a list of disease identifiers (CUIs, MeSH, OMIMs...)
#' @param vocabulary   The vocabulary of the disease identifier(s)
#' Select one of the available:  \code{UMLS} (UMLS), \code{OMIM} (OMIM),
#' \code{MESH} (MeSH), \code{DO} (Disease Ontology),
#' \code{NCI} (NCI thesaurus), \code{ORDO} (Orphanet),
#' \code{ICD9CM} (ICD9-CM) or \code{EFO} (EFO). Default \code{'UMLS'}.
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
#' @param verbose By default \code{FALSE}. Change it to \code{TRUE} to get a
#' on-time log from the function.
#' @param warnings By default \code{TRUE}. Change it to \code{FALSE} to don't see
#' the warnings.
#' @return An object of class \code{DataGeNET.DGN}
#' @examples
#' dis_res <- disease2attribute( "C0028754", database = "CURATED", vocabulary="UMLS")
#' @export disease2attribute



disease2attribute <- function( disease, vocabulary= "UMLS", database = "CURATED",  verbose = FALSE, warnings = TRUE ) {

  vocabularies <- c("UMLS", "OMIM", "NCI", "MESH", "ICD9CM", "EFO", "DO", "HPO", "ORDO")
  vv <- intersect(vocabulary, vocabularies)
  if (length(vv) == 0) {
    stop(paste0(vocabulary,  " is a wrong vocabulary! Allowed vocabularies are: \n", paste(vocabularies, collapse = "\n")))
  }
  if (vv != "UMLS") {
    vv <- paste0(tolower(vv), "_id")
  }
  check_disgenet_sources( database )
  if( length( disease ) != length( unique( disease ) ) ) {
    disease <- unique( disease )
    if( warnings ) {
      warning(
        "Removing duplicates from input diseases list."
      )
    }
  }

  list_of_diseases <- paste(disease,collapse=",")
  result <- NULL


  if (vocabulary== "UMLS") {
    #http://localhost:5555/disgenetv6/api/disease/C0004153,C0002395
    url <- paste0( get_url_disgenet(), "disease/", list_of_diseases , "?format=tsv")
    #url <- paste0( getUrlDGN(), "disease/", list_of_diseases , "?source=",database, "&format=tsv")

  } else {
    url <- paste0( get_url_disgenet(), "disease/",tolower(vocabulary),"/", list_of_diseases ,"?source=",database,"&format=tsv")
  }

  #print(url)
  r <- httr::GET(url)
  if (r$status_code == 200) {
    myTextConnection <- textConnection(rawToChar(r$content ))
    result <- read.csv( myTextConnection, header = TRUE, sep = "\t" )
    close(myTextConnection)
  }else{
    print(httr::http_status(r))
    print(httr::content(r, "text"))
  }
  # x =RCurl::getURL(url)
  # if (grepl("status_code",x) == F) {
  #   dataTsv <- RCurl::getURLContent( url     )
  #   myTextConnection <- textConnection( dataTsv )
  #   result <- read.csv( myTextConnection, header = TRUE, sep = "\t", colClasses=c("character"))
  #   close(myTextConnection)
  # }
  # else{
  #   x <- jsonlite::fromJSON(x)
  #   print(x$detail)
  # }

  if(length( disease )==1){
    cls <- "single"
  }else {
    cls <- "list"
  }

  if (!is.null(result)){
    if(vocabulary == "UMLS"){
      wDiseases <- setdiff(disease, result$diseaseid)
      eDiseases <- intersect(disease, result$diseaseid)
    }else {
      wDiseases <- setdiff(disease,result[, vv])
      eDiseases <- intersect(disease, result[, vv])
    }

    if( length( wDiseases ) != 0 ) {
      diseases <- paste( paste( "   -", wDiseases ), collapse = "\n" )
      if( warnings ) {
        warning( "One or more of the diseases in the list is not in DisGeNET ( '", database, "' ):\n", diseases  )
      }
    }
    result <- new( "DataGeNET.DGN",
                   type       = "disease",
                   search     = cls,
                   term       = as.character( eDiseases ),
                   database   = database,
                   scoreRange = "",
                   qresult    = result
    )
    return( result )
  } else {
    print("no results for the query")
  }
}

