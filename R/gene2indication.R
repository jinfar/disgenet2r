#' Query for given gene(s) and generates an \code{DataGeNET.RDF}
#
#' Given the NCBI gene identifier, retrieves the information
#' related to drugs which interacts with the gene product and
#' their indications and creates an object of type \code{DataGeNET.RDF}.
#
#' @param gene NCBI gene ID
#' @param verbose By default \code{FALSE}. Change it to \code{TRUE} to get a
#' on-time log from the function.
#' @param warnings By default \code{TRUE}. Change it to \code{FALSE} to don't see
#' the warnings.
#' @return An object of class \code{DataGeNET.RDF}
# @examples
# gene2indi <- gene2indication( gene = 1588 )
#' @export gene2indication


gene2indication <- function( gene,  verbose = FALSE, warnings = TRUE ) {

    # endpoint
  endpoint <- "http://rdf.disgenet.org/sparql/"

  # oql is the query that let us look for drug indications in Wikidata
  query <- "PREFIX wd: <http://www.wikidata.org/entity/>
  PREFIX wdt: <http://www.wikidata.org/prop/direct/>

  SELECT str(?ncbigene_id) as ?ncbi_geneid ?gene_iri str(?gene_product_label) as ?gene_product_label str(?chembl_id) as ?chembl_id str(?drug_label) as ?drug_label ?drug_iri str(?indication_label) as ?indication_label ?indication_iri WHERE {

    # Wikidata query
    SERVICE <https://query.wikidata.org/sparql> {
      ?gene wdt:P351 ?ncbigene_id . FILTER (?ncbigene_id = \"%s\")
      ?gene wdt:P688 ?gene_product .  # gene_product (usually a protein) is a product of a gene (a region of DNA)
      ?gene_product rdfs:label ?gene_product_label .
      filter (lang(?gene_product_label) = 'en')
      ?drug_iri wdt:P129 ?gene_product ;   # drug interacts with a gene_product
      wdt:P2175 ?indication_iri ;
      wdt:P592 ?chembl_id .
      ?drug_iri rdfs:label ?drug_label . filter (lang(?drug_label) = 'en')
      ?indication_iri rdfs:label ?indication_label . filter (lang(?indication_label) = 'en')
  }

  # DisGeNET query
  BIND(IRI(CONCAT(\"http://identifiers.org/ncbigene/\",?ncbigene_id)) as ?gene_iri)
  ?gene_iri rdf:type ncit:C16612 .
  }"

  query <- sprintf(query,  as.character(gene))

  # SPARQL query
  # perform the SPARQL query to DisGeNET
  prefix <- c("sio","http://semanticscience.org/resource/","skos","http://www.w3.org/2004/02/skos/core#")
  message("Performing the query. This is a federated query with Wikidata and may take some time, please wait...")
  res <- tryCatch({
    SPARQL(url=endpoint,query=query,ns=prefix)
  },
  error = function(err){
    print(paste0("The query did not work. The error: "),err)
  })

  if (length(res$results) == 0){
    warning("There is no indications data in Wikidata for gene ", gene)
    #warning("There are no genes ", gene, " and the input parameters that have annotations in Wikidata.")
  } else {

    result <- new( "DataGeNET.RDF",
                   input     = as.character(gene),
                   search    = "indication",
                   selection = "gene",
                   mapping   = nrow(res$result),
                   qresult   = res$result
    )
    return( result )
  }

}
