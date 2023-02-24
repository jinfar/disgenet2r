#' Performs an enrichment analysis on a list of genes/variants and generates an object \code{DataGeNET.DGN}
#'
#' @param entities Name or vector of genes (Entrez identifiers or HGNC symbols) or variants (DBSNP identifiers) to
#' for the enrichment
#' @param vocabulary The vocabulary for the genes (ENTREZ or HGNC) and for the variants (DBSNP)
#' @param database Name of the database that will be queried to retrieve the gene/variant disease annotations.
#' The possible values are:
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
#' @param universe The universe to be used for the Fisher test.
#' The possible values are:
#' \code{DISGENET} - All genes in DisGeNET
#' \code{HUMAN} - All genes according to the NCBI
#' \code{HUMAN_CODING} - All protein coding genes to the NCBI
#' \code{CUSTOM} - A list of genes or variants submitted by the user
#' \code{custom_universe} a list of genes/variants as Entrez or symbols/dbSNP identifiers to be used
#' as universe when setting the universe="CUSTOM"
#' @param verbose By default \code{FALSE}. Change it to \code{TRUE} to get a
#' on-time log from the function.
#' @param warnings By default \code{TRUE}. Change it to \code{FALSE} to not see
#' the warnings.
#' @return An object of class \code{DataGeNET.Dis}
#' @examples
#' res <- disease_enrichment( entities =c("APP", "PSEN","APOE"), vocabulary = "HGNC", database = "CURATED", universe = "HUMAN_CODING")
#' @export disease_enrichment


disease_enrichment <- function(entities=entities, universe = "DISGENET", custom_universe ="", vocabulary = "HGNC",
                               verbose = TRUE, database = "CURATED", warnings = TRUE) {
  check_disgenet_sources( database )

  if( length(entities) != length( unique( entities ) ) ) {
    entities <- unique( entities )
    warning(
      "Removing duplicates from input list."
    )
  }

  type <- "symbol"
  if( vocabulary == "HGNC" ){
    #do nothing
  } else if( vocabulary == "ENTREZ" ){
    type <- "entrezid"
  } else if( vocabulary == "DBSNP" ){
    type <- "dbsnp"
  }else if( class( entities[1] ) == "factor"){
    message("Your input entities are in factor format.
    Please, revise your input entities and save them as numeric or as character.")
    stop()
  } else {
    message(paste0("Your input entities are a wrong vocabulary ", vocabulary,
                   ". Please, revise your input vocabulary. Remember that genes should be identified
                   using HGNC gene symbols or NCBI Entrez identifiersm and variants as DBSNP identifiers"))
    stop()
  }

  # if (!  universe %in% c("DISGENET", "HUMAN", "HUMAN_CODING", "CUSTOM")){
  #   message(paste0("Your input universe  ", universe,
  #                  "  should be one of the following:
  #                 - DISGENET - All genes/variants in DisGeNET \\n
  #                 - HUMAN - All genes according to the NCBI or \\n
  #                 - HUMAN - All variants according to CLINVAR, UNIPROT, all variants from GWASCAT and GWASDB without filtering by p-value. \\n
  #                 - HUMAN_CODING - All protein coding genes to the NCBI \\n
  #                 - CUSTOM - a list of genes supplied by the user \\n
  #                  "))
  #   stop()
  #
  # }
  # if (universe == "CUSTOM" & length(custom_universe ) < 2  ){
  #   message("To use the CUSTOM universe you must supply your own list of genes or variants!")
  #   stop()
  # }

  list_of_entities = paste(entities, collapse = ",")

  if (vocabulary == "DBSNP"){
      url <- paste0( get_url_disgenet(), "enrichment/variants")
      body <-  paste0("variants=",list_of_entities, "&source=",database)
      # if (universe== "CUSTOM"){
      #   custom_universe = paste(custom_universe, collapse = ",")
      #   body <-  paste0("variants=",list_of_entities, "&universe=CUSTOM&source=",database, "&custom_universe=", custom_universe)
      # } else {
      #   body <-  paste0("variants=",list_of_entities, "&universe=",universe, "&source=",database)
      # }

  } else if (vocabulary %in% c("HGNC", "ENTREZ" ) ){
    url <- paste0( get_url_disgenet(), "enrichment/genes")
    body <-  paste0("genes=",list_of_entities, "&typeid=",type ,  "&source=",database)
    # if (universe== "CUSTOM"){
    #   custom_universe = paste(custom_universe, collapse = ",")
    #   body <-  paste0("genes=",list_of_entities, "&typeid=",type ,"&universe=CUSTOM&source=",database, "&custom_universe=", custom_universe)
    # }
    # else{
    #       body <-  paste0("genes=",list_of_entities, "&typeid=",type ,"&universe=",universe, "&source=",database)
    #
    # }
  } else{
    message(paste0("Your input vocabulary  ", vocabulary,
                   "  should be one of the following:
                  - HGNC \\n
                  - ENTREZ  \\n
                  - DBSNP \\n"))
    stop()
  }

  r  <- httr::POST(url = url,
                     httr::add_headers('Content-Type' = "application/x-www-form-urlencoded"),
                     body = body)

  if (r$status_code == 200) {
    res <-jsonlite::fromJSON(httr::content(r, as = "text", encoding = "UTF-8"), flatten = F)
  }else{
    print(httr::http_status(r))
    print(httr::content(r, "text"))
    stop()

  }

  if (type == "dbsnp") {
    data <- res$results[,c("diseaseid", "disease_name", "source", "variant_ratio", "bg_ratio", "pvalue", "adjusted_pvalue")]
    # shared_entities <-   cbind(res$results$diseaseid, purrr::map_dfr(res$results$intersection, ~as_tibble(t(.), .name_repair = 'unique')))
    shared_entities <-  sapply(res$results$intersection, '[', seq(max(sapply(res$results$intersection, length))))
    shared_entities <- as.data.frame(t(shared_entities),stringsAsFactors =F)
    shared_entities <- cbind(res$results$diseaseid, shared_entities)
    shared_entities <- reshape2::melt(shared_entities, id.vars="res$results$diseaseid")
    colnames(shared_entities) <- c("diseaseid", "v1", "shared_variant")
    shared_entities <- aggregate(shared_variant~diseaseid, data= shared_entities , function(x)paste(x,collapse=";"))

  } else{

    l <- res$background$not_found_genes
    if(length(l) > 0){
      message(paste0(l, " gene(s) from the input list not found in DisGeNET ", database))
    }
    data <- res$results[,c("diseaseid", "disease_name", "source", "gene_ratio", "bg_ratio", "pvalue", "adjusted_pvalue")]
    shared_entities <-   cbind(res$results$diseaseid, res$results$intersection)
    shared_entities <- reshape2::melt(shared_entities, id.vars="res$results$diseaseid")
    colnames(shared_entities) <- c("diseaseid", "shared_geneid", "shared_symbol")
    shared_entities <- subset(shared_entities, !is.na(shared_symbol))
    ag1 <- aggregate(shared_geneid~diseaseid, data= shared_entities , function(x)paste(x,collapse=";"))
    ag2 <- aggregate(shared_symbol~diseaseid, data= shared_entities , function(x)paste(x,collapse=";"))
    shared_entities <- merge(ag1, ag2, by = "diseaseid")

  }

  if( dim( res$results$disease_class)[2] > 0){
    disclass <-   cbind(res$results$diseaseid, res$results$disease_class)
    colnames(disclass)[1] <-"diseaseid"
    disclass <- reshape2::melt(disclass, id = "diseaseid")
    disclass <- unique(disclass[! is.na(disclass$value), ])
    disclass$variable <- as.character(disclass$variable)
    disclass$value <- as.character(disclass$value)
    disclass <- disclass[ order(disclass$variable), ]
    disclass1 <- aggregate( variable ~ diseaseid, data= disclass ,  function(x)paste(x,collapse=";"))
    disclass2 <- aggregate( value ~ diseaseid, data= disclass ,  function(x)paste(x,collapse=";"))
    disclass <- merge(disclass1, disclass2, by = "diseaseid")
    colnames(disclass) <- c("diseaseid", "disease_class","disease_class_name")
    data <- merge(data, disclass, by = "diseaseid",all.x = T)
  } else {
    data$disease_class <- NA
    data$disease_class_name <- NA
  }

  disclass <-   cbind(res$results$diseaseid, res$results$semantic_type)
  colnames(disclass)[1] <-"diseaseid"
  disclass <- reshape2::melt(disclass, id = "diseaseid")
  disclass <- unique(disclass[! is.na(disclass$value), ])
  disclass$variable <- as.character(disclass$variable)
  disclass$value <- as.character(disclass$value)
  disclass <- disclass[ order(disclass$variable), ]
  disclass <- aggregate( value ~ diseaseid, data= disclass ,  function(x)paste(x,collapse=";"))
  colnames(disclass) <- c("diseaseid", "disease_semantic_type")

  data <- merge(data, disclass, by = "diseaseid",all.x = T)
  data <- merge(data, shared_entities, by = "diseaseid",all.x = T)

  colnames(data)[1:7] <- c("ID", "Description", "source",  "Ratio", "BgRatio", "pvalue", "FDR")
  data$Count <- sapply(strsplit(data$Ratio, "/"), `[`, 1)
  data$Count <- as.numeric(as.character(data$Count))

  data$gg <- data$Count/length(entities)
  data <- data[order(as.numeric(as.character(data$FDR))), ]
  result <- new( "DataGeNET.DGN",
                 type       = "disease-enrichment",
                 search     = "list",
                 term       = as.character( entities ),
                 database   = database,

                 qresult    = data
  )
  return( result )

}
