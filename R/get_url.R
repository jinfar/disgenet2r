get_url_disgenet <- function() {
  url <- "https://www.disgenet.org/api/"
  # url <- "http://miranda/disgenetv6/api/"
  # url <- "http://localhost:5555/disgenetv6/api/"
  return( url )
}


get_api_connection <- function(url = url, api_key = NULL) {

  if(is.null(api_key)){
      api_key = Sys.getenv('DISGENET_API_KEY')
      if(api_key == ""){
        stop("This is not a valid API KEY! Please, use the function `get_disgenet_api_key` to get your API key from DisGeNET")
      }
  }

  api_key <- paste0('Bearer ', api_key )
  res <- httr::GET(url = url,
                   httr::add_headers(.headers=c(  'accept' = "application/json",
                                                  'Authorization' = api_key)
                   ))

  return( res )
}
