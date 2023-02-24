#' Query for given gene(s) and generates an \code{DataGeNET.RDF}
#'
#' Given a NCBI gene identifier, retrieves the GO (Gene Ontology) information
#' related to biological processes where the gene product is involved and creates
#' an object of type \code{DataGeNET.RDF}.
#'
#' @param gene NCBI gene ID
#' @param verbose By default \code{FALSE}. Change it to \code{TRUE} to get a
#' on-time log from the function.
#' @param warnings By default \code{TRUE}. Change it to \code{FALSE} to don't see
#' the warnings.
#' @return An object of class \code{DataGeNET.RDF}
#' @examples
#' gene2bp <- gene2biologicalprocess( gene = 1588 )
#' @export gene2biologicalprocess


gene2biologicalprocess <- function( gene, verbose = FALSE, warnings = TRUE ) {

  # endpoint
  endpoint <- "http://rdf.disgenet.org/sparql/"
  message("Performing the query. This is a federated query with Wikidata and may take some time, please wait...")

    # oql is the query that let us look for drug indications in Wikidata
    query <- "PREFIX wd: <http://www.wikidata.org/entity/>
      PREFIX wdt: <http://www.wikidata.org/prop/direct/>

      SELECT DISTINCT str(?ncbigene_id) as ?ncbi_geneid ?gene_iri
str(?gene_product_label) as ?gene_product_label
str(?go_id) as ?go_id str(?go_label) as ?go_label ?biological_process_iri WHERE {

      # Wikidata query
      SERVICE <https://query.wikidata.org/sparql> {
        ?gene wdt:P351 ?ncbigene_id . FILTER (?ncbigene_id =\"%s\")
        ?gene wdt:P688 ?gene_product .  # gene_product (usually a protein) is a product of a gene (a region of DNA)
        ?gene_product rdfs:label ?gene_product_label .
        filter (lang(?gene_product_label) = 'en')
        ?drug wdt:P129 ?gene_product .   # drug interacts with a gene_product
        ?gene_product wdt:P682 ?biological_process_iri . #add information about the GO cell component
        ?biological_process_iri rdfs:label ?go_label ;
                        wdt:P686 ?go_id .
        filter (lang(?go_label) = 'en')
      }

      # DisGeNET query
      BIND(IRI(CONCAT(\"http://identifiers.org/ncbigene/\",?ncbigene_id)) as ?gene_iri)
      ?gene_iri rdf:type ncit:C16612 .
      }"


    query <- sprintf(query,  as.character(gene))
    #print(query)

    # SPARQL query
    # perform the SPARQL query to DisGeNET
    prefix <- c("sio","http://semanticscience.org/resource/","skos","http://www.w3.org/2004/02/skos/core#")
    #message("Performing the query. This is a federated query with Wikidata and may take some time, please wait...")
    res <- tryCatch({
      SPARQL(url=endpoint,query= query,ns=prefix)
    },
    error = function(err){
      print(paste0("The query did not work. The error: "),err)
    })

    if (length(res$results) > 0){
      results <- res$results
      result <- new( "DataGeNET.RDF",
                     input     = as.character(gene),
                     search    = "GOBP",
                     selection = "gene",
                     mapping   = nrow(results),
                     qresult   = results
      )
      return( result )

    } else {
      if( verbose ) {
        message("The gene ",gene, " has no annotations in Wikidata.")
      }
    }



}
