# Private function. Given an RDF data frame, returns the formatted data frame
# input: inID = input ID, InVocabulary = input vocabulary
# output: uri = the resulting URI. If the format is incorrect it returns NA
getFormattedDataFrame <- function( df, verbose = FALSE, warnings = FALSE  ) {

  #df <- as.data.frame(apply(df, 1:2, function(x) strsplit(gsub(">", "", x), "\\/")[[1]][length(strsplit(x, "\\/")[[1]])] ))

  # df$disease <- gsub("<http://linkedlifedata.com/resource/umls/id/", "", df$disease)
  # df$disease <- gsub(">", "", df$disease)
  # df["gene"] <- apply(df["gene"], 1:2, function(x) strsplit(gsub(">", "", x), "\\/")[[1]][length(strsplit(x, "\\/")[[1]])] )
  colnames(df)[which(colnames(df) == "disease")] <- "diseaseid"
  colnames(df)[which(colnames(df) == "gene")] <- "geneid"
  colnames(df)[which(colnames(df) == "geneSymbol")] <- "geneSymbol"
  colnames(df)[which(colnames(df) == "geneName")] <- "description"
  df$source <- gsub("<http://rdf.disgenet.org/v6.0.0/void/","", df$source)
  df$source <- gsub(">","", df$source)
  df$diseaseid <- gsub("umls:","", df$diseaseid)
  df$geneid <- gsub("ncbigene:","", df$geneid)
  if (length(df[, colnames(df)=="molecule"]) == 0){
    #df <- as.data.frame(apply(df, 1:2, function(x) strsplit(gsub(">", "", x), "\\/")[[1]][length(strsplit(x, "\\/")[[1]])] ))
    df["source"] <- apply(df["source"], 1:2, function(x) strsplit(gsub(">", "", x), "\\/")[[1]][length(strsplit(x, "\\/")[[1]])] )

    #colnames(df)[which(colnames(df) == "source")] <- "sourceId"
    #df$sourceId <- gsub("-2016", "", df$sourceId)
    df$source <- toupper(df$source )
    df$source <- gsub("CTD_HUMAN", "CTD_human", df$source)
    df$source <- gsub("CTD_RAT", "CTD_rat", df$source)
    df$source <- gsub("CTD_MOUSE", "CTD_mouse", df$source)
  }
  else{
    df["uniprot"] <- apply(df["uniprot"], 1:2, function(x) strsplit(gsub(">", "", x), "\\/")[[1]][length(strsplit(x, "\\/")[[1]])] )

  }
  if (length(df[, colnames(df)=="snp"]) >  0){
    df["snp"] <- apply(df["snp"], 1:2, function(x) strsplit(gsub(">", "", x), "\\/")[[1]][length(strsplit(x, "\\/")[[1]])] )
  }
  return(df);
}
