#' Retrieves genes associated to a disease, or list of diseases and generates an \code{DataGeNET.DGN}
#'
#' Given the disease identifier for one or multiple diseases retrieves their associated genes
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
#' @param score A vector with two elements: 1) initial value of score 2) final value of score
#' @param verbose By default \code{FALSE}. Change it to \code{TRUE} to get a
#' on-time log from the function.
#' @param warnings By default \code{TRUE}. Change it to \code{FALSE} to don't see
#' the warnings.
#' @return An object of class \code{DataGeNET.DGN}
#' @examples
#' dis_res <- disease2gene( "C0028754", database = "CURATED" , score = c(0,1))
#' @export disease2gene



disease2gene <- function( disease, vocabulary= "UMLS", database = "CURATED", score=c(0,1),api_key=NULL,  verbose = FALSE, warnings = TRUE ) {

  if(is.null(api_key)){
    api_key = Sys.getenv('DISGENET_API_KEY')
    if(api_key == ""){
      stop("This is not a valid API KEY! Please, use the function `get_disgenet_api_key` to get your API key from DisGeNET")
    }
  }
  vocabularies <- c("UMLS", "OMIM", "NCI", "MESH", "ICD9CM", "ICD10", "EFO", "DO", "HPO", "ORDO")
  vv <- intersect(vocabulary, vocabularies)
  if (length(vv) == 0) {
    stop(paste0(vocabulary,  " is an invalid vocabulary! ! Allowed vocabularies are: \n", paste(vocabularies, collapse = "\n")))
  }
  if (vv != "UMLS") {
    vv <- paste0(tolower(vv), "id")
  }
  check_disgenet_sources( database )
  if( length( disease ) != length( unique( disease ) ) ) {
    disease <- unique( disease )
    if( warnings ) {get_umls_from_vocabulary
      warning(
        "Removing duplicates from input diseases list."
      )
    }
  }


  if(length(score) != 2) {
    stop("Invalid argument 'score'. It must have two elements: initial value of the score, and final value of the score")
  }

  list_of_diseases <- paste(disease,collapse=",")
  result <- NULL

  if (vocabulary== "UMLS") {
      url <- paste0( get_url_disgenet(), "gda/disease/", list_of_diseases , "?source=",database, "&min_score=",score[1],"&max_score=", score[2],  "&format=tsv")
  } else {
    url <- paste0( get_url_disgenet(), "gda/disease/",tolower(vocabulary),"/", list_of_diseases ,"?source=",database, "&min_score=",score[1],"&max_score=", score[2],  "&format=tsv")
  }

  # r <- httr::GET(url)
  r <- get_api_connection(url = url, api_key=api_key)
  if (r$status_code == 200) {
    myTextConnection <- textConnection(rawToChar(r$content ))
    result <- read.csv( myTextConnection, header = TRUE, sep = "\t" )
    close(myTextConnection)
  }else{
    print(httr::http_status(r))
    print(httr::content(r, "text"))
  }


  if(length( disease )==1){
    cls <- "single"
  }else {
    cls <- "list"
  }

  scR <- paste0( score[1],"-", score[2])

  if (!is.null(result)){
    if(vocabulary == "UMLS"){
      wDiseases <- setdiff(disease, result$diseaseid)
      eDiseases <- intersect(disease, result$diseaseid)
      if( length( wDiseases ) != 0 ) {
        diseases <- paste( paste( "   -", wDiseases ), collapse = "\n" )
        if( warnings ) {
          warning( "One or more of the diseases in the list is not in DisGeNET ( '", database, "' ):\n", diseases  )
        }
      }
     }


    result$score <- as.numeric(as.character(result$score))
    result <- new( "DataGeNET.DGN",
                   type       = "disease-gene",
                   search     = cls,
                   term       = as.character( list_of_diseases ),
                   database   = database,
                   scoreRange = scR,
                   qresult    = result
    )
    return( result )
  } else {
    print("no results for the query")
  }
}
