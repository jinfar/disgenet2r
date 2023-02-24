#' Query for a disease or a list of diseases and generates an \code{DataGeNET.RDF}
#'
#' Given a disease or a list of diseases, retrieves the GO BP of the gene products
#' associated to the disease(s) and creates an object of type \code{DataGeNET.RDF}.
#'
#' @param disease a CUI or a vector of CUIs
#' @param database Name of the database that will be queried.
#' It can take the values \code{'CTD_human'} to use Comparative
#' Toxicogenomics Database, human data; \code{'UNIPROT'} to use Universal
#' Protein Resource;\code{'CLINVAR'} to use ClinVar, a public archive of relationships
#' among sequence variation and human phenotype; \code{'GWASCAT'} to use
#' the NHGRI-EBI GWAS Catalog; \code{'ORPHANET'}, to use
#' Orphanet, the portal for rare diseases and orphan drugs;
#' \code{'CURATED'} to use expert curated, human databases;
#' \code{'RGD'}, to use Rat Genome Database; \code{'MGD'}, to use the Mouse Genome Database;
#' \code{'CTD_rat'} to use Comparative Toxicogenomics Database, rat data;
#' \code{'CTD_mouse'} to use Comparative Toxicogenomics Database, mouse data;
#' \code{'PREDICTED'} to use the expert curated, animal models data;
#' \code{'ALL'} to use all these databases. Default \code{'CURATED'}.
#' @param score A vector with two elements: 1) character with greather
#' \code{'>'} or with lower \code{'<'} meaing greather or equal and lower or
#' equal; 2) the score to be compared. By default: \code{c('>', 0)}.
#' @param verbose By default \code{FALSE}. Change it to \code{TRUE} to get a
#' on-time log from the function.
#' @param warnings By default \code{TRUE}. Change it to \code{FALSE} to don't see
#' the warnings.
#' @return An object of class \code{DataGeNET.RDF}
#' @examples
#' disease2bp <- disease2molecularfunction( disease = "C1859062" )
#' @export disease2molecularfunction


disease2molecularfunction <- function( disease, vocabulary = "umls", database = "CURATED", score = c(0,1), verbose = FALSE, warnings = TRUE ) {
  check_disgenet_sources(database);

  if(length(score) != 2) {
    stop("Invalid argument 'score'. It must have two elements: initial value of the score, and final value of the score")
  }
  gdas <- disease2gene(disease, vocabulary= "UMLS", database =  database , score = score)
  genes <- unique(gdas@qresult$geneid)
  message(paste0("There are ",length(genes) ," genes to be queried on the
                 Wikidata SPARQL endopoint which may take some time, please wait..."));



  results <- data.frame()
  for (gene in genes ){
    if( verbose ) {
      message(
        "Retrieving data for gene ", gene, " from Wikidata"
      )
    }

    g2mf <- gene2molecularfunction(gene = gene )
    if(!is.null(g2mf)){
      results <- rbind(results, g2mf@qresult)
    } else{
      warning(paste0("There were no annotations for ", gene, " in Wikidata"))
    }
  }


  if (dim(results)[1] > 0){
    result <- new( "DataGeNET.RDF",
                   input     = as.character(disease),
                   search    = "GOMF",
                   selection = "disease",
                   mapping   = nrow(results),
                   qresult   = results
    )
    return( result )
  }else {
    if( verbose ) {
      message("The genes associated to disease ", disease, " have no annotations in Wikidata.")
    }
  }

}
