#' Retrieves the DisGeNET API key for a user
#'
#' Given the email and password, retrieves the users API Key
#' @param email the email of the user
#' @param password the password of the user to connect to DisGeNET
#' @param verbose By default \code{FALSE}. Change it to \code{TRUE} to get a
#' on-time log from the function.
#' @param warnings By default \code{TRUE}. Change it to \code{FALSE} to don't see
#' the warnings.
#' @return A string
#' @examples
#' disgenet_api_key <- get_disgenet_api_key( email = "user@gmail.com", password = "myspwd" )
#' @export get_disgenet_api_key



get_disgenet_api_key <- function( email, password,  verbose = FALSE, warnings = TRUE ) {
  url <- paste0(get_url_disgenet(), "auth/")
  r  <- httr::POST(url = url, httr::add_headers(.headers=c(
    `accept` =  '*/*' )), body = list(    `email` = email,    `password` = password  )  )

  if (r$status_code == 200) {
    api_key <- jsonlite::fromJSON(httr::content( r, as = "text", encoding = "UTF-8"), flatten = F)$token
    #Sys.setenv('DISGENET_API_KEY' = api_key)
  }else{
    print(httr::http_status(r))
    print(httr::content(r, "text"))
    api_key <-  NULL
  }
   return(api_key)
}
