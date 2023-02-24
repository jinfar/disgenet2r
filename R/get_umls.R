#' Retrieves mappings from different vocabularies associated to a disease, or list of diseases and generates an \code{DataGeNET.DGN}
#'
#' Given the disease identifier for one or multiple diseases retrieves the mappings to other vocabularies
#' from DisGeNET and creates an object of type \code{DataGeNET.DGN}.
#'
#' @param disease  A disease or a list of disease identifiers (CUIs, MeSH, OMIMs...)
#' @param vocabulary   The vocabulary of the disease identifier(s)
#' Select one of the available:  \code{UMLS} (UMLS), \code{OMIM} (OMIM),
#' \code{MESH} (MeSH), \code{DO} (Disease Ontology),
#' \code{NCI} (NCI thesaurus), \code{ORDO} (Orphanet),
#' \code{ICD9CM} (ICD9-CM) or \code{EFO} (EFO). Default \code{'UMLS'}.
#' @param verbose By default \code{FALSE}. Change it to \code{TRUE} to get a
#' on-time log from the function.
#' @param warnings By default \code{TRUE}. Change it to \code{FALSE} to don't see
#' the warnings.
#' @return An object of class \code{DataGeNET.DGN}
#' @examples
#' dis_res <- get_umls_from_vocabulary( "	D000130", vocabulary = "MESH"  )
#' @export get_umls_from_vocabulary



get_umls_from_vocabulary <- function(disease, vocabulary,  api_key=NULL,
                                     allowed_vocabularies = c("UMLS",
                                                              "OMIM",
                                                              "NCI",
                                                              "MESH",
                                                              "ICD9CM",
                                                              "ICD10",
                                                              "EFO",
                                                              "DO",
                                                              "HPO",
                                                              "ORDO", "NAME")) {

  if(is.null(api_key)){
    api_key = Sys.getenv('DISGENET_API_KEY')
    if(api_key == ""){
      stop("This is not a valid API KEY! Please, use the function `get_disgenet_api_key` to get your API key from DisGeNET")
    }
  }
  if (!vocabulary %in% allowed_vocabularies) {
      stop(
        paste0(  vocabulary, " is an invalid vocabulary! Allowed vocabularies are: \n",
          paste(allowed_vocabularies, collapse = "\n")  )
      )
  }
  else{

    list_of_diseases <- paste(unique(disease),collapse=",")

    if (vocabulary!= "NAME"){
      url <-   paste0(get_url_disgenet(), "disease/mappings/", tolower(vocabulary),"/", list_of_diseases    )
    }
    else{
      disease <- gsub(" ","%20", disease)
      url <-   paste0(get_url_disgenet(), "disease/mappings/name/", disease    )
      #print(url)
    }

    r <- get_api_connection(url = url, api_key=api_key)
    if (r$status_code == 200) {
      result <-  jsonlite::fromJSON(httr::content(r, as = "text", encoding = "UTF-8"), flatten = F)
      return(result)
    } else{
      print(httr::http_status(r))
      print(httr::content(r, "text"))
    }
  }

}
