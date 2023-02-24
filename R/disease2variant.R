#' Retrieves variants associated to a disease, or list of diseases and generates an \code{DataGeNET.DGN}
#'
#' Given the name of one or multiple diseases and retrieves the associated variant(s)
#'  and creates an object of type \code{DataGeNET.RDF}.
#' @param disease A disease or a list of disease identifiers (CUIs, MeSH, OMIMs...)
#' @param vocabulary   The vocabulary of the disease identifier(s)
#' Select one of the available: \code{OMIM} (OMIM), \code{MSH} (MeSH), \code{DO} (Disease Ontology),
#' \code{NCI} (NCI thesaurus), \code{ORDO} (Orphanet),
#' \code{ICD9CM} (ICD9-CM) or \code{EFO} (EFO). Default \code{'UMLS'}.
#' @param score A vector with two elements: 1) initial value of score 2) final value of score
#' @param verbose By default \code{FALSE}. Change it to \code{TRUE} to get a
#' on-time log from the function.
#' @param warnings By default \code{TRUE}. Change it to \code{FALSE} to hide the warnings.
#' @param database Name of the database that will be queried. It can take the values:
#' \code{'UNIPROT'} to use Universal Protein Resource;
#' \code{'CLINVAR'} to use ClinVar, a public archive of relationships
#' among sequence variation and human phenotype;
#' \code{'GWASCAT'} to use the NHGRI-EBI GWAS Catalog;
#' \code{'GWASDB'} to use the GWAS Database GWASdb;
#' \code{'CURATED'} to use expert curated, human databases;
#' \code{'BEFREE'} to use text mining data, generated using BeFree System;
#' \code{'ALL'} to use all these databases. Default \code{'CURATED'}.
#' @return An object of class \code{DataGeNET.DGN}
#' @examples
#' dis2var <- disease2variant( disease = "C0751955",vocabulary = "UMLS", database = "ALL", score=c(0,1) )
#' @export disease2variant


disease2variant <- function( disease = "C0751955", vocabulary = "UMLS", database = "CURATED", score = c(0,1), api_key=NULL,  verbose = FALSE, warnings = TRUE ) {
  if(is.null(api_key)){
    api_key = Sys.getenv('DISGENET_API_KEY')
    if(api_key == ""){
      stop("This is not a valid API KEY! Please, use the function `get_disgenet_api_key` to get your API key from DisGeNET")
    }
  }
  check_disgenet_sources(database)
  vocabularies <- c("UMLS", "OMIM", "NCI", "MESH", "ICD9CM", "EFO", "DO", "HPO", "ORDO")
  vv <- intersect(vocabulary, vocabularies)
  if (length(vv) == 0) {
    stop(paste0(vocabulary,  " is a wrong vocabulary! Allowed vocabularies are: \n", paste(vocabularies, collapse = "\n")))
  }
  if (vv != "UMLS") {
    vv <- paste0(tolower(vv), "_id")
  }
  if (  database == "CTD_human" | database == "CTD_rat" |
        database == "CTD_mouse" |database == "RGD" | database == "MGD" |
        database == "ANIMAL_MODELS" | database == "HPO" | database == "PSYGENET" |
        database == "LHGDN"  ) {

    stop("This resource does not have variant-disease associations.")
  }


  if( length( disease ) != length( unique( disease ) ) ) {
    disease <- unique( disease )
    warning(
      "Removing duplicates from input list."
    )
  }

  if(length(score) != 2) {
    stop("Invalid argument 'score'. It must have two elements: initial value of the score, and final value of the score")
  }

  list_of_diseases <- paste(disease,collapse=",")

  result <- NULL
  if (vocabulary== "UMLS") {
    url <- paste0( get_url_disgenet(), "vda/disease/", list_of_diseases , "?source=",database, "&min_score=", score[1],"&max_score=", score[2],  "&format=tsv")
  } else {
    url <- paste0( get_url_disgenet(), "vda/disease/",tolower(vocabulary),"/", list_of_diseases , "?source=", database, "&min_score=",score[1],"&max_score=", score[2],  "&format=tsv")
  }
  #print(url)
  r <- get_api_connection(url = url, api_key=api_key)

  #r <- httr::GET(url)
  if (r$status_code == 200) {
    myTextConnection <- textConnection(rawToChar(r$content ))
    result <- read.csv( myTextConnection, header = TRUE, sep = "\t" )
    close(myTextConnection)
  }else{
    print(httr::http_status(r))
    print(httr::content(r, "text"))
  }
  # x = RCurl::getURL(url)
  # if (grepl("status_code",x) == F) {
  #   dataTsv <- RCurl::getURLContent( url     )
  #   myTextConnection <- textConnection( dataTsv )
  #   result <- read.csv( myTextConnection, header = TRUE, sep = "\t" )
  #
  #   close(myTextConnection)
  # }
  # else{
  #   x <- jsonlite::fromJSON(x)
  #   print(x$detail)
  # }

  if(length( disease )==1){
    cls <- "single"
  }else{
    cls <- "list"
  }

  scR <- paste0( score[1],"-", score[2])
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

    result$score <- as.numeric(as.character(result$score))
    result <- new( "DataGeNET.DGN",
                   type = "disease-variant",
                   search = cls,
                   term = as.character( disease ),
                   scoreRange = scR,
                   database = database,
                   qresult = result )

    return( result )
  } else {
    print("no results for the query")
  }
}
