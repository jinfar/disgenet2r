#' Retrieves diseases similar to an input disease and generates an \code{DataGeNET.DGN}
#'
#' Given the disease identifier for one or multiple diseases retrieves their asssociated genes
#' from DisGeNET and creates an object of type \code{DataGeNET.DGN}.
#'
#' @param disease  A disease or a list of disease identifiers (CUIs, MeSH, OMIMs...)
#' specific genes from DisGeNET. The genes non contained in DisGeNET will
#' be removed from the output.
#' @param verbose By default \code{FALSE}. Change it to \code{TRUE} to get a
#' on-time log from the function.
#' @param warnings By default \code{TRUE}. Change it to \code{FALSE} to not see
#' the warnings.
#' @return An object of class \code{DataGeNET.Dis}
#' @examples
#' dis_res <- get_similar_diseases( "C0028754" )
#' @export get_similar_diseases

get_similar_diseases <- function(disease, api_key=NULL, limit = 5   ) {

  if(is.null(api_key)){
    api_key = Sys.getenv('DISGENET_API_KEY')
    if(api_key == ""){
      stop("This is not a valid API KEY! Please, use the function `get_disgenet_api_key` to get your API key from DisGeNET")
    }
  }
  list_of_diseases <- paste(unique(disease),collapse=",")
  url <-   paste0(get_url_disgenet(), "disease/similarity/", list_of_diseases , "?source=ALL&format=json&limit=", limit   )
  #print(url)
  # r <- httr::GET(url) http://miranda/disgenetv6/api/disease/similarity/C0002395?format=json
  r <- get_api_connection(url = url, api_key=api_key)
  if (r$status_code == 200) {
    res<-jsonlite::fromJSON(httr::content(r, as = "text", encoding = "UTF-8"), flatten = F)

    result <-res$similar_diseases[[1]]
    result$ref_disease <- res$ref_disease[1]
    result$name <- res$disease_name[1]
    if (length(res$similar_diseases)> 1) {
      for (i in 2:length(res$similar_diseases)) {
        rr <-res$similar_diseases[[i]]
        rr$ref_disease <- res$ref_disease[i]
        rr$name <- res$disease_name[i]
        result <- rbind(result, rr)
      }

    }

    colnames(result)[1]<- "similar_disease"
    colnames(result)[2]<- "similar_disease_name"
    return(result)
  } else{
    print(httr::http_status(r))
    print(httr::content(r, "text"))
  }

}
