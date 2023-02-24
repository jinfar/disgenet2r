check_disgenet_sources <- function( database,
                            allowed_db = c( "CTD_HUMAN", "UNIPROT", "HPO",
                                            "ORPHANET", "PSYGENET", "CGI",
                                            "GENOMICS_ENGLAND", "CLINGEN",
                                            "CLINVAR", "GWASCAT", "GWASDB",
                                            "CURATED",
                                            "CTD_MOUSE", "CTD_RAT","MGD", "RGD",
                                            "ANIMAL_MODELS", "INFERRED",
                                             "LHGDN","BEFREE", "ALL" ) ) {
  database <- toupper( database )
  if( sum( database %in% allowed_db ) == 1 ) {
    return( database )
  } else {
    stop("Invalid given 'database' ('", database, '). Try ?disease2gene or ',
         '?gene2disease to check the available databases.')
  }
}


