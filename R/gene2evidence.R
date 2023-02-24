#' Retrieves evidences supporting a gene-disease association, for a gene, or list of genes and generates an \code{DataGeNET.DGN}
#'
#' Given the NCBI identifier, or HGNC symbol of one or multiple genes, retrieves their associated diseases
#' from DisGeNET and creates an object of type \code{DataGeNET.Dis}.
#'
#' @param gene a gene, or a list of genes identified by NCBI identifiers, or HGNC symbols.
#' The genes non contained in DisGeNET will
#' be removed from the output.
#' @param vocabulary   The vocabulary of the gene identifier(s)
#' Select one of the available:  \code{HGNC} (Gene Symbols), \code{ENTREZ} (NCBI Gene Identifiers)
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
#' @param year A vector with two elements: 1) initial value of year 2) final value of year
#' @param verbose By default \code{FALSE}. Change it to \code{TRUE} to get a
#' on-time log from the function.
#' @param warnings By default \code{TRUE}. Change it to \code{FALSE} to not see
#' the warnings.
#' @return An object of class \code{DataGeNET.Dis}
#' @examples
#' g <- gene2evidence( 3953, "CURATED", score=c(0,1), vocabulary = "ENTREZ" )
#' @export gene2evidence


gene2evidence <- function( gene, vocabulary= "HGNC",
                           disease = NULL, year=c(0, 2050), limit = 10,
                           database = "CURATED", score=c(0,1), api_key=NULL,  verbose = FALSE, warnings = TRUE ) {
  if(is.null(api_key)){
    api_key = Sys.getenv('DISGENET_API_KEY')
    if(api_key == ""){
      stop("This is not a valid API KEY! Please, use the function `get_disgenet_api_key` to get your API key from DisGeNET")
    }
  }
  check_disgenet_sources( database )

  if( length( gene ) != length( unique( gene ) ) ) {
    gene <- unique( gene )
    warning(
      "Removing duplicates from input genes list."
    )
  }

  type <- "symbol"
  if( vocabulary == "HGNC" ){
    #do nothing
  } else if( vocabulary == "ENTREZ" ){
    type <- "entrez"
  } else if( class( gene[1] ) == "factor"){
    message("Your input genes are in factor format.
    Please, revise your input genes and save them as numeric or as character.")
    stop()
  } else {
    message(paste0("Your input genes are a wrong vocabulary ", vocabulary,
                   "Please, revise your input vocabulary. Remember that genes should be identified
                   using HGNC gene symbols or NCBI Entrez identifiers"))
    stop()
  }

  if(length(score) != 2) {
    stop("Invalid argument 'score'. It must have two elements: initial value of the score, and final value of the score")
  }
  if(length(year) != 2) {
    stop("Invalid argument 'year'. It must have two elements: initial year and final year")
  }
  list_of_genes <- paste(gene,collapse=",")
  list_of_diseases <- paste(disease,collapse=",")

  result <- NULL
  if(!is.null(disease)){
    url <- paste0( get_url_disgenet(), "gda/evidences/gene/", list_of_genes , "?disease=",list_of_diseases , "&source=",database, "&min_score=",score[1],"&max_score=", score[2], "&min_year=",year[1],"&max_year=", year[2])
  }
  else{
    url <- paste0( get_url_disgenet(), "gda/evidences/gene/", list_of_genes , "?source=",database, "&min_score=",score[1],"&max_score=", score[2], "&min_year=",year[1],"&max_year=", year[2])

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



  if(length( gene )==1){
    cls <- "single"
  }else{
    cls <- "list"
  }

  scR <- paste0( score[1],"-", score[2])

  if (! is.null(result) > 0){
    if(type == "symbol"){
      wGenes <- setdiff(gene, result$gene_symbol)
      eGenes <- intersect(gene, result$gene_symbol)
    }else if (type == "entrez"){
      wGenes <- setdiff(gene,result$geneid)
      eGenes <- intersect(gene, result$geneid)
    }

    if( length( wGenes ) != 0 ) {
      genes <- paste( paste( "   -", wGenes ), collapse = "\n" )
      if( warnings ) {
        warning( "\n One or more of the genes in the list is not in DisGeNET ( '", database, "' ):\n", genes )
      }
    }

    result <- new( "DataGeNET.DGN",
                   type = "evidence",
                   search = cls,
                   term = as.character( eGenes ),
                   scoreRange = scR,
                   database = database,
                   qresult = result )

    return( result )
  } else {
    print("no results for the query")
  }
}
