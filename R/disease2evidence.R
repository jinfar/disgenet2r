#' Retrieves evidences supporting the associations, for a disease, or a list of diseases,
#' and generates an \code{DataGeNET.DGN}
#'
#' @param disease  A disease or a list of diseases identified by the Concept Unique Identifiers
#' @param variant a variant or a list of variants
#' @param gene a gene, or a list of genes identified by NCBI identifiers, or HGNC symbols.
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
#' @param type one of two possible values: "GDA" or "VDA" depending if the user wants to retrieve gene-disease associations or
#' variant-disease associations.
#' @param year A vector with two elements: 1) initial value of year 2) final value of year
#' @param verbose By default \code{FALSE}. Change it to \code{TRUE} to get a
#' on-time log from the function.
#' @param warnings By default \code{TRUE}. Change it to \code{FALSE} to not see
#' the warnings.
#' @return An object of class \code{DataGeNET.Dis}
#' @examples
#' g <- disease2evidence( disease="C0002395", "CURATED", score=c(0,1)  )
#' @export disease2evidence


disease2evidence <- function( disease=disease, type = 'GDA',
                           gene = NULL, variant=NULL, year=c(0, 2050), limit = 10,
                           database = "CURATED", score=c(0,1), api_key=NULL,  verbose = FALSE, warnings = TRUE ) {
  if(is.null(api_key)){
    api_key = Sys.getenv('DISGENET_API_KEY')
    if(api_key == ""){
      stop("This is not a valid API KEY! Please, use the function `get_disgenet_api_key` to get your API key from DisGeNET")
    }
  }
  check_disgenet_sources( database )

  if( length( disease ) != length( unique( disease ) ) ) {
    disease <- unique( disease )
    warning(
      "Removing duplicates from input disease list."
    )
  }

  # type <- "symbol"
  # if( vocabulary == "HGNC" ){
  #   #do nothing
  # } else if( vocabulary == "ENTREZ" ){
  #   type <- "entrez"
  # } else if( class( gene[1] ) == "factor"){
  #   message("Your input genes are in factor format.
  #   Please, revise your input genes and save them as numeric or as character.")
  #   stop()
  # } else {
  #   message(paste0("Your input genes are a wrong vocabulary ", vocabulary,
  #                  "Please, revise your input vocabulary. Remember that genes should be identified
  #                  using HGNC gene symbols or NCBI Entrez identifiers"))
  #   stop()
  # }
  if( type =="GDA" & !is.null(variant)) {
    stop("You must select 'type=VDA' if you are going to query disease2evidence with a variant or a list of variants")
  }
  if( type =="VDA" & !is.null(gene)) {
    message("You are using the argument 'type=VDA' and querying with genes. This will retrieve the variants in the gene(s) associated to diseases.
            To obtain gene-disease associations, change the argument 'type' to 'GDA'")
  }


  if(length(score) != 2) {
    stop("Invalid argument 'score'. It must have two elements: initial value of the score, and final value of the score")
  }

  if(length(year) != 2) {
    stop("Invalid argument 'year'. It must have two elements: initial year and final year")
  }
  list_of_genes <- paste(gene,collapse=",")
  list_of_diseases <- paste(disease,collapse=",")
  list_of_variants <- paste(variant,collapse=",")

  result <- NULL

  if(type == "GDA"){
    if(!is.null(gene)){
      url <- paste0( get_url_disgenet(), "gda/evidences/disease/", list_of_diseases , "?gene=", list_of_genes, "&source=",database, "&min_score=",score[1],"&max_score=", score[2], "&min_year=",year[1],"&max_year=", year[2])
    }
    else{
      url <- paste0( get_url_disgenet(), "gda/evidences/disease/", list_of_diseases , "?source=",database, "&min_score=",score[1],"&max_score=", score[2], "&min_year=",year[1],"&max_year=", year[2])

    }
  } else if(type == "VDA"){
    if (!is.null(variant) & !is.null(gene)){
      url <- paste0( get_url_disgenet(), "vda/evidences/disease/", list_of_diseases , "?variant=",list_of_variants , "&gene=",list_of_genes, "&source=",database, "&min_score=",score[1],"&max_score=", score[2], "&min_year=",year[1],"&max_year=", year[2])
    }
    else  if(!is.null(variant)){
      url <- paste0( get_url_disgenet(), "vda/evidences/disease/", list_of_diseases , "?variant=",list_of_variants , "&source=",database, "&min_score=",score[1],"&max_score=", score[2], "&min_year=",year[1],"&max_year=", year[2])
    } else   if(!is.null(gene)){
      url <- paste0( get_url_disgenet(), "vda/evidences/disease/", list_of_diseases , "?gene=",list_of_genes , "&source=",database, "&min_score=",score[1],"&max_score=", score[2], "&min_year=",year[1],"&max_year=", year[2])
    }
    else{
      url <- paste0( get_url_disgenet(), "vda/evidences/disease/", list_of_diseases , "?source=",database, "&min_score=",score[1],"&max_score=", score[2], "&min_year=",year[1],"&max_year=", year[2])
    }
  }
  else {
    message(paste0("The argument 'type' is invalid: ", type,
                   "Please, revise your input type Remember that it should be VDA or GDA"))
    stop()
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



  if(length( disease )==1){
    cls <- "single"
  }else{
    cls <- "list"
  }

  scR <- paste0( score[1],"-", score[2])

  if (! is.null(result) > 0){
   missing_diseases <- setdiff(disease, result$disease_id)
   found_diseases <- intersect(disease, result$disease_id)

    if( length( missing_diseases ) > 0 ) {
      diseases <- paste( paste( "   -", missing_diseases ), collapse = "\n" )
      if( warnings ) {
        warning( "\n One or more of the genes in the list is not in DisGeNET ( '", database, "' ):\n", diseases )
      }
    }

    result <- new( "DataGeNET.DGN",
                   type = "evidence",
                   search = cls,
                   term = as.character( found_diseases ),
                   scoreRange = scR,
                   database = database,
                   qresult = result )

    return( result )
  } else {
    print("no results for the query")
  }
}
