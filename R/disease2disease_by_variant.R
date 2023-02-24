#' Retrieves disease-disease associations via shared genes and generates an \code{DataGeNET.DGN}
#'
#' Given the CUI of one or multiple diseases, retrieves the associated diseases
#' from DisGeNET and creates an object of type \code{DataGeNET.Dis}.
#'
#' @param disease a CUI or a vector of CUIs
#' The diseases non existing in DisGeNET will
#' be removed from the output.
#' @param vocabulary   The vocabulary of the disease identifier(s)
#' Select one of the available:  \code{UMLS} (UMLS), \code{OMIM} (OMIM), \code{MSH} (MeSH), \code{DO} (Disease Ontology),
#' \code{NCI} (NCI thesaurus), \code{ORDO} (Orphanet),
#' \code{ICD9CM} (ICD9-CM) or \code{EFO} (EFO). Default \code{'UMLS'}.
#' @param database Name of the database that will be queried. It can take the values:
#' \code{'UNIPROT'} to use Universal Protein Resource;
#' \code{'CURATED'} to use expert curated, human databases;
#' \code{'GWASCAT'} to use the NHGRI-EBI GWAS Catalog;
#' \code{'GWASDB'} to use the GWAS Database GWASdb;
#' \code{'CLINVAR'} to use ClinVar, a public archive of relationships
#' among sequence variation and human phenotype;
#' \code{'BEFREE'} to use text mining data, generated using BeFree System;
#' \code{'ALL'} to use all these databases. Default \code{'CURATED'}.
#' @param ndiseases The number of associated diseases to retrieve. By default 10. Set to XX to retrieve all diseases
#' @param pvalue p value of the bootstrap. See more info at www.disgenet.org. By default 1
#' @param verbose By default \code{FALSE}. Set to \code{TRUE} to get a
#' on-time log from the function.
#' @param warnings By default \code{TRUE}. Set to \code{FALSE} to hide the warnings.
#' @return An object of class \code{DataGeNET.Dis}
#' @examples
#' dd1 <- disease2disease_by_variant( "C0028754", database = "CURATED" , pvalue = 0.05)
#' @export disease2disease_by_variant


disease2disease_by_variant <- function( disease, vocabulary= "UMLS",  ndiseases =10, database = "CURATED", pvalue = 1,  api_key=NULL,  verbose = FALSE, warnings = TRUE ) {
  if(is.null(api_key)){
    api_key = Sys.getenv('DISGENET_API_KEY')
    if(api_key == ""){
      stop("This is not a valid API KEY! Please, use the function `get_disgenet_api_key` to get your API key from DisGeNET")
    }
  }
  check_disgenet_sources( database )
  if (  database == "CTD_human" | database == "CTD_rat" |
        database == "CTD_mouse" |database == "RGD" | database == "MGD" |
        database == "ANIMAL_MODELS" | database == "HPO" | database == "PSYGENET" |
        database == "LHGDN"  ) {

    stop("This resource does not have variant-disease associations.")
  }


  vocabularies <- c("UMLS", "OMIM", "NCI", "MESH", "ICD9CM", "EFO", "DO", "HPO", "ORDO")
  vv <- intersect(vocabulary, vocabularies)
  if (length(vv) == 0) {
    stop(paste0(vocabulary,  " is a wrong vocabulary! Allowed vocabularies are: \n", paste(vocabularies, collapse = "\n")))
  }
  if (vv != "UMLS") {
    vv <- paste0(tolower(vv), "_id")
  }

  if( length( disease ) != length( unique( disease ) ) ) {
    disease <- unique( disease )
    if( warnings ) {
      warning(
        "Removing duplicates from input diseases list."
      )
    }
  }

  result    <- data.frame()
  for (ii in 1:length(disease) ){

    if( verbose ) {
      message(
        "Retrieving data for ", disease[ii], " from DisGeNET"
      )
    }
    if (vocabulary== "UMLS") {
      url <- paste0( get_url_disgenet(), "dda/variants/disease/", disease[ii] ,"?limit=", ndiseases, "&source=", database, "&format=tsv")
    } else {
      url <- paste0( get_url_disgenet(), "dda/variants/disease/",tolower(vocabulary),"/", disease[ii] ,"?limit=", ndiseases, "&source=", database, "&format=tsv")
    }
    r <- get_api_connection(url = url, api_key=api_key)

    if (r$status_code == 200) {
      myTextConnection <- textConnection(rawToChar(r$content ))
      dataNew <- read.csv( myTextConnection, header = TRUE, sep = "\t" )
      result <- rbind(result, dataNew)
      close(myTextConnection)
    } else if (r$status_code == 404){
      message(
        "The disease ", disease[ii], " is not present in DisGeNET, data source ", database
      )
    } else{
      #print(httr::http_status(r))
      print(httr::content(r, "text"))
    }

  }

  if(length( disease )==1){
    cls <- "single"
  }else {
    cls <- "list"
  }
  # print(result)
  result <- result[result$nvariants > 0,]
  # if (!is.null(result)){
  if ( dim(result)[2]> 0){
    if(vocabulary == "UMLS"){
      wDiseases <- setdiff(disease, result$diseaseid1)
      eDiseases <- intersect(disease, result$diseaseid1)
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
    result <- pvalue_estimation(result, database = database, nboot = 1000, api_key = api_key)
     result$pvalue <- as.numeric(as.character(result$pvalue))
    result$pvalue <- round(result$pvalue, 3)
    result$jaccard_variants <- as.numeric(as.character(result$jaccard_variants))
    result$jaccard_variants <- round(result$jaccard_variants, 3)
    result <- result[ result$pvalue < pvalue  , ]

    result <- result[c( "diseaseid1", "disease1_name", "diseaseid2", "disease2_name",
                        "disease1_disease_class", "disease1_disease_class_name",
                        "disease2_disease_class", "disease2_disease_class_name",
                        "jaccard_variants", "pvalue", "nvariants",
                        "nvariants1" ,"nvariants2" ,  "source"  )]

    result <- new( "DataGeNET.DGN",
                   type     = "disease-disease-variant",
                   search   = cls,
                   term       = as.character( eDiseases ),
                   database = database,
                   scoreRange =  "",
                   qresult  = result
    )

    return( result )
  } else{
    message(
      "There are no results for those parameters "
    )
  }

}
