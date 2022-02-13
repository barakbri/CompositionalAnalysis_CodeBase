#' Shiny GUI for ancom.R
#'
#' @description Opens a graphical user interface to run \pkg{ancom.R}, built from the \pkg{shiny} package.
#'
#' @rdname shiny_ancom
#'
#' @param input  input from GUI.
#' @param output output to GUI.
#' 
#' @details 
#' This GUI is the primary purpose for the package \pkg{ancom.R}.
#' 
#' 
#' @examples
#' \dontrun{ shiny_ancom() }
#' 
#' @import shiny
#' @export
#' 


shiny_ancom <- function(){

  runApp(
    list(
      ui     = shiny_ancomUI,
      server = shiny_ancomServer
    )
  )
}



