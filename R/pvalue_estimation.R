#' Estimates the p value of the Jaccard Coefficient for pairs of diseases
#'
#' This function estimates the statistic significance of the Jaccard Coefficient with a bootstrap
#' @param object receives dataframe produced by function disgenetDisDis
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
#' @param nboot Number of iterations sued to compute the pvalue associted
#' to the calculated Jaccard Index. Default: 1000.
#' @param ncores Number of cores used to calculate the pvalue associated to
#' the computed Jaccard Index. Default: 1.
#' @param verbose By default \code{FALSE}. Change it to \code{TRUE} to get a
#' on-time log from the function.
#' @examples
#' ji <- pvalue_estimation( object = dis, "CURATED" )
#' @export pvalue_estimation



pvalue_estimation <- function(input, database="CURATED", api_key=NULL,  nboot = 100, ncores = 1, verbose = FALSE) {
  if(is.null(api_key)){
    api_key = Sys.getenv('DISGENET_API_KEY')
    if(api_key == ""){
      stop("This is not a valid API KEY! Please, use the function `get_disgenet_api_key` to get your API key from DisGeNET")
    }
  }
  check_disgenet_sources(database)

  colnames(input)<- tolower(colnames(input))
  if ( length(intersect("jaccard_genes", colnames(input))) > 0 ){
    type <- "gene"
    url <- paste0( get_url_disgenet(), "gene/source/",database, "?format=tsv")
    r <- get_api_connection(url = url, api_key=api_key)

    if (r$status_code == 200) {
      myTextConnection <- textConnection(rawToChar(r$content ))
      result <- read.csv( myTextConnection, header = TRUE, sep = "\t" )
      close(myTextConnection)
    }else{
      #print(httr::http_status(r))
      print(httr::content(r, "text"))
    }
    input$jaccard_index <- input$jaccard_genes
  }   else if ( length(intersect("jaccard_variants", colnames(input))) > 0 ){
    type <- "variant"
    load(system.file("extdata", "universe.rda", package="disgenet2r"))

    result <-unique(universe_vdas[,c("variantid", "source")])
    result <- result[ result$source==database,  ]

    input$jaccard_index <- input$jaccard_variants
  }


  if (type =="gene"){
    universe <- as.character(unique(result$geneid))
  } else {
    universe <- as.character(unique(result$variantid))
  }

  message(paste0("A total of ", length(universe), " ", type,
            "s obtained from DisGeNET database ", database, " are
                 being used for the bootstrap process\n")
  )

    input$pvalue <- NA
    message("Pvalue estimation")

    for(i in 1:nrow( input )){
      if( verbose ) {
        message(paste0("\n\t->Disease pair ", i, " of ", nrow(input), " total diseases' pairs."))
      }
      if (type == "gene"){
        bb <- ji.internal(input$ngenes1[i], input$ngenes2[i], universe, nboot, ncores)
      } else {
        bb <- ji.internal(input$nvariants1[i], input$nvariants2[i], universe, nboot, ncores)
      }
      pvalue <- (sum(bb > input$jaccard_index[i]) * 1.0) / (nboot)
      input$pvalue[i] <- round(pvalue, 5)
    }
    return(input)

  }

  ji.internal <- function(len1, len2, universe, nboot, ncores) {
    if (!requireNamespace("parallel", quietly = TRUE)) {
      pfun <- lapply
    } else {
      pfun <- parallel::mclapply
    }

    unlist(pfun(1:nboot, function(ii) {
      g1 <- sample( universe, len1 )
      g2 <- sample( universe, len2 )
      ja.coefr <- length(intersect(g1, g2)) / length(union(g1, g2))
    }, mc.cores = ncores))
  }
