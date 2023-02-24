#' Query for given disease(s) and generates an \code{DataGeNET.RDF}
#'
#' Given the UMLS ID (CUI) for a disease, retrieves the associated genes,
#' and the pathway(s) in WikiPathways to which the genes are annotated.
#' It creates an object of type \code{DataGeNET.RDF}.
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
#' dis2path <- disease2pathway( disease = "C0751955" )
#' @export disease2pathway


disease2pathway <- function( disease , database = "CURATED", score = c(0,1), verbose = FALSE, warnings = TRUE ) {

  # this function, gives you a list of pathway(s) from WikiPathways for an specific disease

  check_disgenet_sources( database )

  DB_ID <- "\"\""
  if (database == "CURATED") {
    DB_ID <- toupper(" \"uniprot|ctd_human|clinvar|orphanet|CLINGEN|GENOMICS_ENGLAND|CGI|psygenet\" ")
  }
  else if (database == "ANIMAL_MODELS") {
    DB_ID <- toupper("\"ctd_mouse|mgd|ctd_rat|rgd\"")
  }
  else if (database == "ALL") {
    DB_ID <- "\"\""
  }
  else{
    DB_ID <- paste0("\"", toupper(database), "\"")
  }


  if(length(score) != 2) {
    stop("Invalid argument 'score'. It must have two elements: initial value of the score, and final value of the score")
  }
  # endpoint
  endpoint <- "http://rdf.disgenet.org/sparql/"



  message("Performing the query. This is a federated query with WikiPathways and may take some time, please wait...")
  query <- "PREFIX wp:      <http://vocabularies.wikipathways.org/wp#>

      SELECT DISTINCT ?diseaseid str(?diseaseName) as ?diseaseName
      ?geneid str(?geneSymbol) as ?geneSymbol str(?geneName) as ?geneName
       ?score ?source ?pathway str(?pathwayName) as ?pathwayName WHERE {
      # Query DisGeNET for disease-genes
      ?disease skos:exactMatch <http://bio2rdf.org/umls:%s> .
      ?gda sio:SIO_000628 ?gene,?disease .
      ?gda sio:SIO_000253 ?source .
     ?gene rdf:type ncit:C16612 ;
      dcterms:identifier ?geneid;
      sio:SIO_000205 ?gSymbol;
      dcterms:title ?geneName .
      ?gSymbol dcterms:title ?geneSymbol .
      ?disease rdf:type ncit:C7057;
      dcterms:identifier ?diseaseid;
      dcterms:title ?diseaseName .
      ?gda sio:SIO_000253 ?source .

      ?gda sio:SIO_000216 ?scoreIRI .
      ?scoreIRI sio:SIO_000300 ?score .
      FILTER (?score >= %f && ?score<= %f)
      FILTER regex( ?source, %s )

      # Query WikiPathways for gene-pathways
      SERVICE <http://sparql.wikipathways.org/sparql> {
      ?geneProduct a wp:GeneProduct .
      ?geneProduct wp:bdbEntrezGene ?gene .
      ?geneProduct rdfs:label ?GeneLabel .
      ?geneProduct dcterms:isPartOf ?pathwayid .
      ?pathwayid dc:identifier ?pathway .
      ?pathwayid dc:title ?pathwayName .
      }
      }
      ORDER BY DESC(?pathwayName)"

      query <- sprintf(query, disease, score[1],score[2],DB_ID)
      #print(    query)
      # SPARQL query
      # perform the SPARQL query to DisGeNET
      #print(oql)
      prefix <- c("sio","http://semanticscience.org/resource/","skos","http://www.w3.org/2004/02/skos/core#")
      res <- tryCatch({
        SPARQL(url=endpoint,query=query,ns=prefix)
      },
      error = function(err){
        print(err)
      })



   if (length(res$results) > 0) {
    scR <- paste0( score[1],"-", score[2])

    results <- res$result
    results<- getFormattedDataFrame(results)
    result <- new( "DataGeNET.RDF",
                   input     = disease,
                   search    = "pathway",
                   database = database,
                   scoreRange = scR,
                   selection = "wikipathways",
                   mapping   = nrow(results),
                   qresult   = results
    )
    return( result )
  }
  else{
    print(paste0("There were no results for this query"))
  }

}
