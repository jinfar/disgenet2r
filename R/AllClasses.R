#' Class DataGeNET.DGN
#'
#' Class \code{DataGeNET.DGN} is one of the two objects used in \code{disgenet2r}
#' package. It is the main data container to using the different functions to
#' query DisGeNET database and generate their output.
#'
#' @name DataGeNET.DGN-class
#' @rdname DataGeNET.DGN-class
#' @exportClass DataGeNET.DGN
#' @slot type type of query: \code{'gene-disease'} or  \code{'variant-disease'}
#'  or \code{'disease-gene'},  or \code{'disease-variant'} or \code{'disease-disease'} or \code{'disease-enrichment'}.
#' @slot search Character containing \code{'single'} of \code{'list'}.It is
#' used to perform the correct query to DisGeNET
#' @slot database Character containing the name of the database that will be
#' queried. It can take the values
#' \code{'CTD_human'} to use Comparative Toxicogenomics Database, human data;
#' \code{'UNIPROT'} to use Universal Protein Resource;
#' \code{'CLINGEN'} to use Clinical Genome Resource;
#' \code{'CGI'} to use Cancer Genome Interpreter;
#' \code{'ORPHANET'}, to use Orphanet, the portal for rare diseases and orphan drugs;
#' \code{'PSYGENET'} to use PSYGENET;
#' \code{'GENOMICS_ENGLAND'} to use GENOMICS England;
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
#' @slot scoreRange the disGeNET score range to retrieved the GDAs or VDAs
#' @slot term Character with the term(s) to search into the database(s).
#' @slot qresult \code{data.frame} with the obtained result
#' @seealso disease2gene, gene2disease, DataGeNET.DGN-methods


setClass( "DataGeNET.DGN",
          representation =
              representation( type       = "character",  # type of query
                              search     = "character",  # single or list
                              database   = "character",  # where to search
                              scoreRange = "character",  # score
                              term       = "character",  # what to search
                              qresult    = "data.frame"  # result
              ),
          prototype =
              prototype( type       = "",
                         search     = "",
                         database   = "",
                         scoreRange = "",
                         term       = "",
                         qresult    = data.frame()
              )
)

#' Class DataGeNET.RDF
#'
#' Class \code{DataGeNET.RDF} is the basic object generated when RDF is needed.
#' It can be passed to other disgenet2r package functions.
#'
#' @name DataGeNET.RDF-class
#' @rdname DataGeNET.RDF-class
#' @exportClass DataGeNET.RDF
#' @slot input Character containing \code{'ontology'}, \code{'phenotype'} or
#' \code{'disease'}. used as input to the query.
#' @slot search Character containing \code{'ontology'},  \code{'phenotype'} or
#' \code{'disease'} it is looking for.It is used to perform the correct query
#' @slot database Character containing the name of the database that will be
#' queried. It can take the values
#' \code{'CTD_human'} to use Comparative Toxicogenomics Database, human data;
#' \code{'UNIPROT'} to use Universal Protein Resource;
#' \code{'CLINGEN'} to use Clinical Genome Resource;
#' \code{'CGI'} to use Cancer Genome Interpreter;
#' \code{'ORPHANET'}, to use Orphanet, the portal for rare diseases and orphan drugs;
#' \code{'PSYGENET'} to use PSYGENET;
#' \code{'GENOMICS_ENGLAND'} to use GENOMICS England;
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
#' @slot scoreRange the disGeNET score range to retrieved the GDAs/VDAs
#' @slot selection Character containing the name of the ontoloyg selected
#' @slot mapping Number with the mappings obtained.
#' @slot qresult \code{data.frame} with the obtained result
#' @seealso ontology2disease, disease2ontology, phenotype2disease, disease2phenotype

setClass( "DataGeNET.RDF",
          representation =
            representation( input     = "character",  # disease, ontology, phenotype, gene
                            search    = "character",  # what is looked for
                            database = "character",  # where to search
                            scoreRange = "character", # score
                            selection = "character",  # resource/namespace where the search is performed
                            mapping   = "numeric",    # number of mappings
                            qresult   = "data.frame"  # result
            ),
          prototype =
            prototype( input     = "",
                       search    = "",
                       database = "",
                       scoreRange = "",
                       selection = "",
                       mapping   = 0,
                       qresult   = data.frame()
            )
)
