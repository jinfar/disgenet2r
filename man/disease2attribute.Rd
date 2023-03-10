% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/disease2attribute.R
\name{disease2attribute}
\alias{disease2attribute}
\title{Retrieves genes associated to a disease, or list of diseases and generates an \code{DataGeNET.DGN}}
\usage{
disease2attribute(
  disease,
  vocabulary = "UMLS",
  database = "CURATED",
  verbose = FALSE,
  warnings = TRUE
)
}
\arguments{
\item{disease}{A disease or a list of disease identifiers (CUIs, MeSH, OMIMs...)}

\item{vocabulary}{The vocabulary of the disease identifier(s)
Select one of the available:  \code{UMLS} (UMLS), \code{OMIM} (OMIM),
\code{MESH} (MeSH), \code{DO} (Disease Ontology),
\code{NCI} (NCI thesaurus), \code{ORDO} (Orphanet),
\code{ICD9CM} (ICD9-CM) or \code{EFO} (EFO). Default \code{'UMLS'}.}

\item{database}{Name of the database that will be queried. It can take the values:
\code{'CTD_human'} to use Comparative Toxicogenomics Database, human data;
\code{'UNIPROT'} to use Universal Protein Resource;
\code{'CLINGEN'} to use Clinical Genome Resource;
\code{'CGI'} to use Cancer Genome Interpreter;
\code{'ORPHANET'}, to use Orphanet, the portal for rare diseases and orphan drugs;
\code{'PSYGENET'} to use PSYGENET;
\code{'GENOMICS_ENGLAND'} to use Genomics England PanelApp;
\code{'CURATED'} to use expert curated, human databases;
\code{'HPO'} to use HPO;
\code{'INFERRED'} to use inferred data from HPO, GWASDB, GWASCAT, and CLINVAR;
\code{'CTD_rat'} to use Comparative Toxicogenomics Database, rat data;
\code{'CTD_mouse'} to use Comparative Toxicogenomics Database, mouse data;
\code{'RGD'}, to use Rat Genome Database;
\code{'MGD'}, to use the Mouse Genome Database;
\code{'ANIMAL_MODELS'} to use the expert curated, animal models data;
\code{'GWASCAT'} to use the NHGRI-EBI GWAS Catalog;
\code{'GWASDB'} to use the GWAS Database GWASdb;
\code{'CLINVAR'} to use ClinVar, a public archive of relationships
among sequence variation and human phenotype;
\code{'BEFREE'} to use text mining data, generated using BeFree System;
\code{'ALL'} to use all these databases. Default \code{'CURATED'}.}

\item{verbose}{By default \code{FALSE}. Change it to \code{TRUE} to get a
on-time log from the function.}

\item{warnings}{By default \code{TRUE}. Change it to \code{FALSE} to don't see
the warnings.}
}
\value{
An object of class \code{DataGeNET.DGN}
}
\description{
Given the name of one or multiple diseases retrieves their asssociated genes
from DisGeNET and creates an object of type \code{DataGeNET.DGN}.
}
\examples{
dis_res <- disease2attribute( "C0028754", database = "CURATED", vocabulary="UMLS")
}
