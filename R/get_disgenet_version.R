#' Retrieves DisGeNET data version
#' @examples
#' dgn_version <- get_disgenet_version()
#' @export get_disgenet_version


get_disgenet_version <- function() {
  url <- paste0( get_url_disgenet(), "version/")
  #x = RCurl::getURL(url)
  #x <- jsonlite::fromJSON(x)
  r <- httr::GET(url)
  r <- httr::content(r, "text", encoding = "UTF-8")
  r <- gsub('\\"', " ", r)
  print(r)
}
