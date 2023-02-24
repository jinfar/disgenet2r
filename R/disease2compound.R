#' Query for given disease(s) and generates an \code{DataGeNET.RDF}
#'
#' Given the UMLS CUI for a disease, it retrieves the information
#' related to the pharmacology in ChEMBL and creates an object of type \code{DataGeNET.RDF}.
#'
#' @param disease UMLS ID (CUI) for a disease
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
#' @return An object of class \code{DataGeNET.RDF}
#' @examples
#' dis2com <- disease2compound( disease = "C3542021", database = "ALL" )
#' @export disease2compound


disease2compound <- function( disease = "C3542021", database = "CURATED", score = c(0,1), verbose = FALSE, warnings = TRUE ) {

  cat(paste0("Started at ", Sys.time(),"\n"));
  message("Please, be aware of ChEMBL restrictions for their SPARQL endpoint \"https://www.ebi.ac.uk/rdf/what-are-limitations-sparql-endpoints\". \nIf 1000 compounds are retrieved for a single protein it probably means that there are more than 1000.");

  # this function, gives you the pharmacological profile from ChEMBL for an specific disease

  check_disgenet_sources(database);

  if(length(score) != 2) {
    stop("Invalid argument 'score'. It must have two elements: initial value of the score, and final value of the score")
  }
  gdas <- disease2gene(disease, vocabulary= "UMLS", database =  database , score = score)
  unis <- unique(gdas@qresult$uniprotid)
  message(paste0("There are ",length(unis) ," proteins to be queried on the
                 ChEMBL SPARQL endopoint which may take some time, please wait..."));

    # endpoint
    endpoint <- "https://www.ebi.ac.uk/rdf/services/sparql";
    results <- data.frame(uniprot=character(),molecule = character(),moleculeLabel=character())


    for( uniprot in  unis){

            query <-    "PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
            PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
            PREFIX owl: <http://www.w3.org/2002/07/owl#>
            PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
            PREFIX dc: <http://purl.org/dc/elements/1.1/>
            PREFIX dcterms: <http://purl.org/dc/terms/>
            PREFIX foaf: <http://xmlns.com/foaf/0.1/>
            PREFIX skos: <http://www.w3.org/2004/02/skos/core#>
            PREFIX cco: <http://rdf.ebi.ac.uk/terms/chembl#>

            SELECT DISTINCT ?uniprot
            ?molecule
            ?moleculeLabel
            WHERE {
            ?targetcmpt cco:targetCmptXref ?uniprot .
            ?target cco:hasTargetComponent ?targetcmpt .
            ?molecule rdfs:label ?moleculeLabel .
            ?assay cco:hasTarget ?target .
            ?activity a cco:Activity ;
            cco:hasMolecule ?molecule ;
            cco:hasAssay ?assay .
            FILTER ( ?uniprot = <http://purl.uniprot.org/uniprot/%s> ) }";


          query <- sprintf(query,  uniprot)

            # prefix <- c("rdfs","http://www.w3.org/2000/01/rdf-schema#","cco","http://rdf.ebi.ac.uk/terms/chembl#");
            cat(paste0("Querying ChEMBL for protein ", uniprot, " ... "));

            res <- tryCatch({
              SPARQL(url=endpoint,query=query )
            },
            error = function(err){
              print(err)
            })
            if (length(res$results) == 0){
              cat(paste0("no compounds found for uniprot ", uniprot, " \n"));

            } else {
              cat(paste0(nrow(res$results)," compounds retrieved for uniprot ", uniprot, " \n"));
              results <- rbind(results, res$results)
            }
          }
      cat(paste0("Ended at ", Sys.time(),"\n"));

      results$uniprot <- gsub("<http://purl.uniprot.org/uniprot/", "", results$uniprot)
      results$uniprot <- gsub(">", "", results$uniprot)
      results <- merge(results, gdas@qresult, by.x = "uniprot", by.y = "uniprotid", all = T)
      scR <- paste0( score[1],"-", score[2])
      result <- new( "DataGeNET.RDF",
                     input     = disease,
                     search    = "compound",
                     database = database,
                     scoreRange = scR,
                     selection = "ChEMBL",
                     mapping   = nrow(results),
                     qresult   = results
      );
    return( result )

}





federated_disease2compound <- function( disease = "C0751955", database = "ALL", score = c(">", 0), verbose = FALSE, warnings = TRUE ) {

  # this function, gives you the pharmacological profile from ChEMBL for an specific disease
  if( verbose ) {
    message( "Starting querying DisGeNET for the disease code ", disease, " and cross its disease genes with ChEMBL data." )
  }
  checkDisGeNET_dbs(database)
  diseaseValidation <-  diseaseMapping(disease, inVocabulary = "umls")
  diseaseValidation <- diseaseValidation@qresult
  disease <- diseaseValidation [diseaseValidation$isFormatOK == T & diseaseValidation$existsInDGN == T, "inputID"]

  if (length(disease)> 0){
    DB_ID <- "\"\""
    if (database == "CURATED") {
      DB_ID <- " \"uniprot|ctd_human|clinvar|orphanet|gwascat\" "
    }
    else if (database == "PREDICTED") {
      DB_ID <- "\"ctd_mouse|mgd|ctd_rat|rgd\""
    }
    else if (database == "ALL") {
      DB_ID <- "\"\""
    }
    else{
      DB_ID <- paste0("\"", tolower(database), "\"")
    }

    if (length(score) != 2) {
      stop("Invalid argument 'score'. It must have two elements.")
    } else if (!score[1] %in% c('>', '<')) {
      stop("Invalid argument 'score'. First element must be '>' or '<'.")
    }

    SCORE <- 0
    if (score[2] > 0) {
      SCORE <- score[2]
    }
    SIGN <- ">"
    if (score[1] ==  "<") {
      SIGN <- "<"
    }
    # endpoint
    endpoint <- "http://rdf.disgenet.org/sparql/"

    # oql is the query that let us look for compounds in ChEMBL
    oql <- "PREFIX cco: <http://rdf.ebi.ac.uk/terms/chembl#>
    SELECT DISTINCT ?disease str(?diseaseName) as ?diseaseName ?gene str(?geneName) as ?geneName  ?uniprot ?source ?score  ?molecule ?moleculeLabel WHERE {
    ?gda sio:SIO_000628 ?gene,?disease .
    ?gda sio:SIO_000253 ?source . FILTER regex(?source, DB_ID )
    ?gda sio:SIO_000216 ?scoreIRI .
    ?scoreIRI sio:SIO_000300 ?score . FILTER (?score SIGN \"SCORE\"^^xsd:decimal)

    ?gene rdf:type ncit:C16612 .
    ?gene dcterms:title ?geneName .
    ?gene sio:SIO_010078 ?protein .
    ?disease rdf:type ncit:C7057 ;
    dcterms:title ?diseaseName ;
    dcterms:identifier ?inputID .
    FILTER regex(?inputID,concat('umls:',?cui))
    values ?cui {'UMLS_CUI'} .
    ?protein skos:exactMatch ?uniprot .
    FILTER regex(?uniprot, 'http://purl.uniprot.org/uniprot/')
    # Query ChEMBL for active molecules
{SELECT ?uniprot ?molecule ?moleculeLabel WHERE {
    SERVICE <http://www.ebi.ac.uk/rdf/services/chembl/sparql> {
    ?molecule rdfs:label ?moleculeLabel .
    ?activity a cco:Activity ;
    cco:hasMolecule ?molecule ;
    cco:hasAssay ?assay .
    ?assay cco:hasTarget ?target .
    ?target cco:hasTargetComponent ?targetcmpt .
    ?targetcmpt cco:targetCmptXref ?uniprot .
    ?uniprot a cco:UniprotRef
    } # end of service
} # end of select
} # end of subquery
    }"

    oql <- stringr::str_replace(string = oql,
                                pattern = "DB_ID",
                                replacement = DB_ID)
    oql <- stringr::str_replace(string = oql,
                                pattern = "SCORE",
                                replacement = SCORE)
    oql <- stringr::str_replace(string = oql,
                                pattern = "SIGN",
                                replacement = SIGN)
    # disease
    # input: query=single, namespace=UMLS. Parse input to substract 'namespace:'.
    if (substr(disease,1,4) == 'umls' & nchar(disease) == 13 & substr(disease,6,6) == 'C' & !is.na(as.numeric(substring(disease,7)))) {
      disease = substr(disease,6,13)
      oql <- stringr::str_replace(
        string = oql,
        pattern = "UMLS_CUI",
        replacement = disease
      )
    } else if (substr(disease,1,1) == 'C' & nchar(substring(disease,2)) == 7 & !is.na(as.numeric(substring(disease,2)))) {
      oql <- stringr::str_replace(
        string = oql,
        pattern = "UMLS_CUI",
        replacement = disease
      )
    }else{
      stop( "\nThe disease introduced: '", disease, "' is not a valid input.
            Remember that only one type of ontology is valid for this query: UMLS CUI (inputted as \"umls:C0028754\" or \"C0028754\").
            Please, check that the disease introduced is a valid UMLS CUI and it is inputted in the correct format, and try again.\n\n" )
    }

    # SPARQL query
    # perform the SPARQL query to DisGeNET
    prefix <- c("sio","http://semanticscience.org/resource/","skos","http://www.w3.org/2004/02/skos/core#")
    message("Performing the query. This is a federated query with ChEMBL and may take some time, please wait...")
    res <- tryCatch({
      SPARQL(url=endpoint,query=oql,ns=prefix)
    },
    error = function(err){
      print(paste0("The query did not work. The error: "),err)
    })
    if (length(res$results) == 0){
      warning("There are no genes for disease ", disease, " and the input parameters that have annotations in ChEMBL.")
    } else {

      result <- new( "DataGeNET.RDF",
                     input     = disease,
                     search    = "compound",
                     selection = "umls",
                     mapping   = nrow(res$result),
                     qresult   = res$result
      )
      return( result )
    }

  }else{
    warning("The  disease ", disease, " is not in DisGeNET (", database, ", ", score, ")")
  }
}
