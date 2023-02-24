# Query for given compound(s) and generates an \code{DataGeNET.RDF}
#
# Given a compound identifier in ChEMBL, it retrieves the related genes and diseases
# in DisGeNET. It creates an object of type \code{DataGeNET.RDF}.
#
# @param CHEMBL ID or vector of CHEMBL IDs
# @param verbose By default \code{FALSE}. Change it to \code{TRUE} to get a
# on-time log from the function.
# @param warnings By default \code{TRUE}. Change it to \code{FALSE} to don't see
# the warnings.
# @return An object of class \code{DataGeNET.RDF}
# @examples
# com2dis <- compound2disease( compound = "CHEMBL281124" )
### @export compound2disease


compound2disease <- function( compound , verbose = FALSE, warnings = TRUE ) {

  # this function, gives you the disease profile from DisGeNET for an specific compound in ChEMBL
  message("Please, be aware of ChEMBL restrictions for their SPARQL endpoint \"https://www.ebi.ac.uk/rdf/what-are-limitations-sparql-endpoints\"");
  if( verbose ) {
    message( "Starting querying ChEMBL for the compound ID ", compound, " and cross its compounds with DisGeNET data." )
  }

  # endpoint
  endpoint <- "http://rdf.disgenet.org/sparql/"

  # oql is the query that let us look for diseases in DisGeNET
  oql <- "PREFIX cco: <http://rdf.ebi.ac.uk/terms/chembl#>
          SELECT DISTINCT
            ?molecule
            ?moleculeLabel
            ?uniprot
            ?gene
            str(?geneName) as ?name
            str(?geneSymbol) as ?geneSymbol
            ?disease
            str(?diseaseName) as ?diseaseName
            WHERE {
            # Query ChEMBL for activity data for the input molecule
              {
                SELECT
                  ?molecule
                  ?moleculeLabel
                  ?uniprot
                WHERE {
                  SERVICE <http://www.ebi.ac.uk/rdf/services/chembl/sparql> {
                    SELECT
                      ?activity
                      ?molecule
                      ?moleculeLabel
                      ?assay
                      ?target
                      ?targetcmpt
                      ?uniprot
                    WHERE{
                      ?activity a cco:Activity ;
                      cco:hasMolecule ?molecule ;
                      cco:hasAssay ?assay .
                      ?molecule rdfs:label ?moleculeLabel.
                      VALUES ?moleculeLabel {'CHEMBL_ID'} .
                      ?assay cco:hasTarget ?target .
                      ?target cco:hasTargetComponent ?targetcmpt .
                      ?targetcmpt cco:targetCmptXref ?uniprot .
                      ?uniprot a cco:UniprotRef
                    } # end chembl query
                  } # end of service
                } # end of select
              } # end of subquery
            ?gda sio:SIO_000628 ?gene,?disease .
            ?gene rdf:type ncit:C16612 ;
            dcterms:title ?geneName ;
            sio:SIO_010078 ?protein ;
            sio:SIO_000205 ?geneSymbolURI .
            ?geneSymbolURI dcterms:title ?geneSymbol .
            ?disease rdf:type ncit:C7057 ;
            dcterms:title ?diseaseName .
            ?protein skos:exactMatch ?uniprot .
            FILTER regex(?uniprot, 'http://purl.uniprot.org/uniprot/')
            }"

  # compound
  # http://rdf.ebi.ac.uk/resource/chembl/molecule/CHEMBL281124
  # input: query=single, namespace=CHEMBL. Parse input to check ID format.
  if (substr(compound,1,6) == 'CHEMBL' & !is.na(as.numeric(substring(compound,7)))) {
    oql <- stringr::str_replace(
      string = oql,
      pattern = "CHEMBL_ID",
      replacement = compound
    )
  } else {
    stop( "\nThe compound ID introduced: '", compound, "' is not a valid input.
          Remember that the valid input is: CHEMBLID (inputted as \"CHEMBL281124\").
          Please, check that the compound introduced is a valid CHEMBL ID and it is inputted in the correct format, and try again.\n\n" )
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
    warning("The compound input with ID ", compound, " does not have drug target(s) encoded by disease genes in the DisGeNET database.")
  } else {
    result <- res$result
    result<- getFormattedDataFrame(result)
    result <- new( "DataGeNET.RDF",
                   input     = compound,
                   search    = "disease",
                   selection = "chembl",
                   mapping   = nrow(result),
                   qresult   = result
    )
    return( result )

  }
}
